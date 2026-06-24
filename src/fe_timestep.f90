module fe_timestep
   !! Adaptive time-stepping strategies for the viscoelastic + sea-level model.
   !! Isolated from fe_coupling so the model can carry more than one stepping strategy
   !! (fixed sub-steps, adaptive, …) behind a common driver.
   !!
   !! `adaptive_stepper` advances the VE+SLE model across a coupling interval [t0,t1]
   !! with the ice load LINEARLY INTERPOLATED between its endpoints, choosing Δt by a
   !! step-doubling local-error estimate on the Maxwell memory (the §3c adaptive-dt
   !! controller). It rests on the 2nd-order trapezoidal memory rule (SCHEME_TRAP) +
   !! the SLE↔memory co-convergence (§3c 3b): each candidate step is taken once at Δt
   !! and once as two Δt/2 sub-steps; the coarse/fine memory difference is the local
   !! error a controller needs (∝ Δt^{p+1}). A rejected step rolls the model state back
   !! (resp%save_state/restore_state) and retries smaller; an accepted step keeps the
   !! more accurate fine state and grows Δt. Δt enters the response only through
   !! Mk = (μ/η)Δt, so changing it (resp%set_dt) is a cheap rescale — no operator
   !! re-factorization (the band LU is Δt-independent).
   use fe_precision,    only: wp
   use fe_sht,          only: sht_grid
   use fe_response,     only: ve_response
   use fe_sle,          only: sle_solver, sle_result
   use fe_viscoelastic, only: scheme_order, scheme_is_implicit
   implicit none
   private

   public :: adaptive_stepper

   type :: adaptive_stepper
      !! Embedded step-doubling controller. Tolerances are on the memory ∞-norm (the
      !! surface observable is algebraic in the memory τ, so this also bounds the
      !! observable's local error). Defaults are placeholders — set per problem.
      real(wp) :: rtol      = 1.0e-4_wp   !! relative local-error tolerance (memory norm)
      real(wp) :: atol      = 1.0e-3_wp   !! absolute local-error floor (memory units)
      real(wp) :: safety    = 0.9_wp      !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp      !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp     !! min Δt shrink per step (reject or grow)
      real(wp) :: dt_min    = 0.0_wp      !! Δt floor (0 = none); a step at the floor is
                                          !! accepted even if over tolerance (no infinite
                                          !! subdivision) and tallied in n_floor
      real(wp) :: dt_max    = huge(1.0_wp)!! Δt ceiling
      real(wp) :: dt_try    = 0.0_wp      !! next-step Δt suggestion (carried across
                                          !! intervals; 0 ⇒ first guess = whole interval)
      integer  :: n_accept  = 0           !! accepted steps (cumulative)
      integer  :: n_reject  = 0           !! rejected attempts (cumulative)
      integer  :: n_floor   = 0           !! steps accepted at dt_min over tolerance
      integer  :: n_solve   = 0           !! SLE solves issued (the dominant cost unit)
      real(wp) :: worst_mass_resid = 0.0_wp  !! worst SLE mass residual over the LAST advance()
   contains
      procedure :: advance => stepper_advance
   end type adaptive_stepper

contains

   subroutine stepper_advance(self, sht, resp, sle, topo0, ice0, ice1, ice_ref, &
                              t0, t1, rsl, C)
      !! Advance the VE+SLE model from t0 to t1 with the ice load interpolated linearly
      !! ice(t) = ice0 + (t−t0)/(t1−t0)·(ice1−ice0), absolute; the SLE load is
      !! d_ice(t) = ice(t) − ice_ref (total change from the reference state). Δt is
      !! chosen adaptively; on return resp holds the end-of-interval memory at t1 and
      !! rsl/C hold the converged sea-level fields there.
      class(adaptive_stepper),  intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      class(ve_response),       intent(inout) :: resp
      type(sle_solver),         intent(inout) :: sle
      real(wp),                 intent(in)    :: topo0(:,:)
      real(wp),                 intent(in)    :: ice0(:,:), ice1(:,:)  !! abs. ice at t0,t1
      real(wp),                 intent(in)    :: ice_ref(:,:)          !! reference ice
      real(wp),                 intent(in)    :: t0, t1
      real(wp),                 intent(inout) :: rsl(:,:)
      real(wp),                 intent(out)   :: C(:,:)

      type(sle_result)      :: res
      real(wp), allocatable :: rsl_n(:,:), ice_now(:,:), dice_now(:,:)
      complex(wp), allocatable :: sig0(:)
      real(wp) :: t, dt, err_inf, tau_inf, errsc, fac, ricfac, expo, span
      integer  :: p, np, nl
      logical  :: at_floor, accept

      span = t1 - t0
      if (span <= 0.0_wp) return
      self%worst_mass_resid = 0.0_wp             ! per-interval diagnostic (reset each advance)
      if (.not. scheme_is_implicit(resp%scheme)) then
         ! 1st-order / explicit response carries no step-doubling order signal here;
         ! the adaptive controller is meaningful only for the trapezoidal scheme.
         ! Fall back to a single fixed step over the interval.
         np = sht%nphi;  nl = sht%nlat
         allocate(ice_now(np,nl), dice_now(np,nl))
         call resp%set_dt(span)
         ice_now = ice1;  dice_now = ice1 - ice_ref
         call sle%solve(sht, resp, dice_now, ice_now, topo0, rsl, C, res)
         self%worst_mass_resid = max(self%worst_mass_resid, res%mass_resid)
         self%n_accept = self%n_accept + 1
         return
      end if

      np = sht%nphi;  nl = sht%nlat
      allocate(rsl_n(np,nl), ice_now(np,nl), dice_now(np,nl))

      ! Seed the trapezoidal start-of-step load σ_0 once, if not already tracked: at the
      ! interval start the response carries its entering memory, so a report-only SLE
      ! solve (no memory/time advance) gives the load consistent with it — at t=0 that
      ! is the elastic-consistent load. Without this the first step falls back to the
      ! σ_{n+1} proxy for σ_n (O(Δt) for the relaxing SLE load → 2nd order on that step).
      if (.not. resp%sigma_primed) then
         allocate(sig0(sht%nlm))
         ice_now = ice0;  dice_now = ice0 - ice_ref
         call sle%solve(sht, resp, dice_now, ice_now, topo0, rsl, C, res, &
                        report_only=.true., sigma_lm=sig0)
         call resp%prime_sigma(sig0)
      end if

      p      = scheme_order(resp%scheme)          ! global order (2 for trapezoidal)
      ricfac = real(2**p - 1, wp)                 ! Richardson factor for the estimate
      expo   = 1.0_wp/real(p + 1, wp)             ! step-size exponent (local err ∝ Δt^{p+1})

      t = t0
      if (self%dt_try <= 0.0_wp) self%dt_try = span
      do
         if (t >= t1 - 1.0e-9_wp*span) exit
         dt = min(self%dt_try, t1 - t)
         at_floor = (self%dt_min > 0.0_wp .and. dt <= self%dt_min*(1.0_wp + 1.0e-9_wp))

         ! --- field step-doubling around [t, t+dt] -------------------------------
         rsl_n = rsl
         call resp%save_state()                   ! buffer A = τ_n
         call resp%set_dt(dt)
         call solve_at(t + dt)                     ! coarse: one Δt
         call resp%stash_coarse()                  ! buffer B = τ_coarse
         call resp%restore_state();  rsl = rsl_n   ! back to τ_n (and its rsl seed)
         call resp%set_dt(0.5_wp*dt)
         call solve_at(t + 0.5_wp*dt)              ! fine sub-step 1
         call solve_at(t + dt)                     ! fine sub-step 2 → τ_fine
         call resp%coarse_fine_error(err_inf, tau_inf)
         call resp%set_dt(dt)                       ! leave resp%dt at the step size
         errsc = (err_inf/ricfac) / (self%atol + self%rtol*tau_inf)

         ! --- accept / reject ----------------------------------------------------
         accept = (errsc <= 1.0_wp) .or. at_floor
         if (accept) then
            t = t + dt;  self%n_accept = self%n_accept + 1
            if (at_floor .and. errsc > 1.0_wp) self%n_floor = self%n_floor + 1
         else
            call resp%restore_state();  rsl = rsl_n
            self%n_reject = self%n_reject + 1
         end if

         ! --- step-size update ---------------------------------------------------
         if (errsc <= 0.0_wp) then
            fac = self%grow_max
         else
            fac = self%safety * errsc**(-expo)
         end if
         fac = min(self%grow_max, max(self%shrink_min, fac))
         self%dt_try = min(self%dt_max, max(self%dt_min, dt*fac))
      end do

   contains

      subroutine solve_at(t_eval)
         !! One SLE solve (co-converged memory advance) with the ice load interpolated
         !! at t_eval. Advances resp by the current resp%dt (set by the caller).
         real(wp), intent(in) :: t_eval
         real(wp) :: frac
         frac = (t_eval - t0)/span
         ice_now  = ice0 + frac*(ice1 - ice0)
         dice_now = ice_now - ice_ref
         self%n_solve = self%n_solve + 1
         call sle%solve(sht, resp, dice_now, ice_now, topo0, rsl, C, res)
         self%worst_mass_resid = max(self%worst_mass_resid, res%mass_resid)
      end subroutine solve_at

   end subroutine stepper_advance

end module fe_timestep
