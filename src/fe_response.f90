module fe_response
   !! Surface-load response operator: the abstraction the sea-level equation
   !! (fe_sle) is built on. Given a spectral surface mass-density load σ_lm
   !! [kg m^-2], it returns the two fields the SLE needs,
   !!
   !!     u_lm  — radial displacement of the solid surface  [m]
   !!     n_lm  — geoid / sea-surface-equipotential height   [m]
   !!
   !! per spherical-harmonic coefficient. The SLE depends only on this interface,
   !! so the elastic and (later) viscoelastic earth responses are swappable.
   !!
   !! Geoid mapping. The per-degree solve returns U(a) and F(a) = φ₁(a), the
   !! surface coefficients of radial displacement and the perturbed gravitational
   !! potential. Martinec's φ₁ carries the load's own direct potential with the
   !! sign OPPOSITE to φ^L (φ₁ → −φ^L for a rigid sphere, k = −F/φ^L − 1; see
   !! fe_radial_fe%loading_love). The geopotential perturbation is therefore −F,
   !! and Bruns' formula gives the geoid height
   !!
   !!     N(a) = −F(a)/g .
   !!
   !! This uses only U and F — NOT the horizontal Love number l, whose sign /
   !! normalization is still being calibrated — so the SLE is not blocked by
   !! that open item. Both U and the −F/g geoid are pinned by the validated
   !! rigid (U→0, 1+k→1) and fluid (U→−(2j+1)/3·φ^L/g, 1+k→0) limits.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_mesh, radial_operator, &
                                 idx_u, idx_v, idx_f, ndof_of
   use fe_viscoelastic,    only: NLAM, ve_strain_constants, dissipative_rhs, &
                                 advance_memory, scheme_is_implicit, &
                                 SCHEME_FE, SCHEME_TRAP
   use fe_sht,             only: sht_grid
   implicit none
   private

   public :: response_operator, elastic_response, null_response, ve_response

   type, abstract :: response_operator
      !! Maps a spectral surface load to surface displacement + geoid.
      !!
      !! Stateful (viscoelastic) responses are AFFINE in the current load at a
      !! fixed time: apply() returns gain(l)·σ_lm + drift_lm, where drift_lm is
      !! frozen from the past-relaxation memory. The SLE fixed point may call
      !! apply() many times per time step with different trial loads (safe,
      !! pure); the surrounding step is bracketed by begin_step (freeze the
      !! drift) and commit_step (advance the memory with the converged load).
      !! For elastic / null responses these brackets are no-ops.
   contains
      procedure(apply_if), deferred :: apply
      procedure :: begin_step  => response_begin_default
      procedure :: commit_step => response_commit_default
      !! Co-convergence hooks for the SLE driver (§3c part 3b): when the load
      !! evolves fast within a step, the load σ and the relaxation memory τ_{n+1}
      !! are a mutual fixed point, so the driver iterates them together. It calls
      !! prepare_endpoint once (snapshot τ_n), then repeatedly: converge σ against
      !! the current τ_{n+1} estimate, advance_endpoint (one trapezoid pass that
      !! also refreshes the report drift to τ_{n+1}), until endpoint_converged;
      !! finalize_step then advances time. For stateless / 1st-order responses
      !! these reduce to a single memory advance (no co-iteration), so they are
      !! no-ops here except advance_endpoint, which stateful drivers override.
      procedure :: prepare_endpoint   => response_prepare_default
      procedure :: advance_endpoint   => response_advance_default
      procedure :: endpoint_converged => response_converged_default
      procedure :: finalize_step      => response_finalize_default
      !! horizontal() returns the spheroidal scalar v_lm = V(a)·σ_lm (+ frozen
      !! drift, for stateful responses) whose surface gradient ∇₁(Σ v_lm Y_lm) is
      !! the horizontal displacement (u_θ, u_φ). It is a pure DIAGNOSTIC — the SLE
      !! fixed point needs only u and N (apply), so it never calls this; a caller
      !! evaluates it ONCE on the converged load. Default (rigid / null) ⇒ zero.
      procedure :: horizontal  => response_horizontal_default
   end type response_operator

   abstract interface
      subroutine apply_if(self, sht, sigma_lm, u_lm, n_lm)
         import :: response_operator, sht_grid, wp
         class(response_operator), intent(inout) :: self
         type(sht_grid),           intent(in)    :: sht
         complex(wp),              intent(in)    :: sigma_lm(:)  !! load [kg m^-2]
         complex(wp),              intent(out)   :: u_lm(:)      !! uplift  [m]
         complex(wp),              intent(out)   :: n_lm(:)      !! geoid   [m]
      end subroutine apply_if
   end interface

   type, extends(response_operator) :: elastic_response
      !! Time-independent (elastic) response: per-degree surface response to a
      !! unit load, precomputed once. Linear and order-independent, so a single
      !! gain per degree multiplies every coefficient of that degree.
      integer  :: lmax = 0
      real(wp) :: g    = 0.0_wp          !! surface gravity [m s^-2]
      real(wp) :: a    = 0.0_wp          !! surface radius  [m]
      real(wp), allocatable :: ugain(:)  !! (0:lmax) U(a) per unit σ_l  [m / (kg m^-2)]
      real(wp), allocatable :: ngain(:)  !! (0:lmax) N(a)=−F(a)/g per unit σ_l
      real(wp), allocatable :: vgain(:)  !! (0:lmax) V(a) per unit σ_l (horizontal)
   contains
      procedure :: init       => elastic_response_init
      procedure :: apply      => elastic_response_apply
      procedure :: horizontal => elastic_response_horizontal
      procedure :: destroy    => elastic_response_destroy
   end type elastic_response

   type, extends(response_operator) :: null_response
      !! Rigid, non-self-gravitating Earth: u ≡ 0 and N ≡ 0. The SLE then
      !! reduces to a uniform (eustatic/barystatic) ocean response, which is the
      !! textbook limit used to check mass conservation and the uniform term.
   contains
      procedure :: apply => null_response_apply
   end type null_response

   type, extends(response_operator) :: ve_response
      !! Viscoelastic field driver. Holds one per-degree saddle-point operator
      !! (assembled + factored once) shared across all orders m, plus an
      !! independent Maxwell memory-stress history per spectral coefficient
      !! (l,m) — each (l,m) load has its own time history, so memory cannot be
      !! collapsed across m. Because the operator and the M = μΔt/η factors are
      !! real, each complex (l,m) history is two real histories (re/im).
      !!
      !! The per-step response is affine: solving the unit load once per degree
      !! gives the elastic gains gu(l), gn(l) AND the nodal field used to update
      !! memory; begin_step solves the frozen memory forcing per (l,m) for the
      !! drift; apply combines them; commit_step advances the memory.
      integer  :: lmax = 0, nr = 0, ne = 0, ndof = 0, nlm = 0
      real(wp) :: g = 0.0_wp, a = 0.0_wp, dt = 0.0_wp, time = 0.0_wp
      !! begin_step skips the drift solve for a coefficient whose Maxwell memory is
      !! below skip_tol × (the largest memory over all coefficients): its drift is
      !! negligible, so it is set to zero rather than solved. Self-consistent (all
      !! memory is zero at t=0 ⇒ all skipped ⇒ exact elastic first step) and cheap
      !! to gate. skip_tol = 0 disables skipping (solve every coefficient).
      real(wp) :: skip_tol = 1.0e-4_wp
      real(wp), allocatable :: mnorm(:)                 !! (nk) max|memory| per slot
      type(radial_operator), allocatable :: ops(:)      !! (1:lmax) per-degree operator
      ! degree-independent element fields
      real(wp), allocatable :: r(:)                     !! node radii (nr)
      real(wp), allocatable :: mu(:), Mk(:)             !! (ne) shear, M=μΔt/η
      ! per-degree constants and unit-load response
      real(wp), allocatable :: Jr(:)                    !! (1:lmax) l(l+1)
      real(wp), allocatable :: nrmc(:,:)                !! (NLAM,1:lmax) Z:Z norms
      real(wp), allocatable :: sa(:,:,:), sb(:,:,:), sc(:,:,:)  !! (4,NLAM,1:lmax)
      real(wp), allocatable :: gu(:), gn(:), gv(:)      !! (0:lmax) elastic gains (gv=V(a))
      real(wp), allocatable :: xUn(:,:), xVn(:,:)       !! (nr,1:lmax) unit-load nodal U,V
      ! Per-(l,m) state is stored in DEGREE-GROUPED order k = 1..nk (all orders m of
      ! a degree l contiguous, l ascending; degree 0 carries no memory and is
      ! excluded). This makes begin_step/commit_step iterate k contiguously AND
      ! reuse the per-degree operator ops(l) across its orders (cache-hot), instead
      ! of the SHTns m-major lm order which switches operator on nearly every solve.
      integer :: nk = 0                                 !! # deforming coeffs (l>=1)
      integer, allocatable :: k2lm(:)                   !! (nk) slot k -> SHTns lm index
      integer, allocatable :: kdeg(:)                   !! (nk) degree l of slot k
      integer, allocatable :: kbeg(:)                   !! (1:lmax+1) first k of each degree
      ! per-(l,m) memory stress, split into real/imag (NLAM,ne,nk)
      real(wp), allocatable :: Are(:,:,:), Aim(:,:,:)
      real(wp), allocatable :: Bre(:,:,:), Bim(:,:,:)
      real(wp), allocatable :: Cre(:,:,:), Cim(:,:,:)
      ! frozen per-step drift (from the memory τ_n), set by begin_step
      complex(wp), allocatable :: dUa(:), dFa(:), dVa(:) !! (nk) surface drift (U,F,V)
      real(wp), allocatable :: dUn_re(:,:), dUn_im(:,:) !! (nr,nk) nodal drift U from τ_n (ε_n)
      real(wp), allocatable :: dVn_re(:,:), dVn_im(:,:) !! (nr,nk) nodal drift V from τ_n
      ! Memory time-integration scheme (see fe_viscoelastic). FE (default) advances
      ! the memory once in commit_step from the report strain — the historical path,
      ! bit-identical. TRAP (2nd-order) is implicit in the end-of-step strain, so
      ! commit_step Picard-iterates the endpoint (frozen load: σ held at the converged
      ! value) — re-solving the drift against the trial τ_{n+1} each pass. §3c part 3a.
      integer  :: scheme = SCHEME_FE
      integer  :: max_couple_iter = 1        !! coupling-iteration cap for implicit schemes
      real(wp) :: couple_tol = 1.0e-6_wp     !! relative surface-drift change to stop iterating
      integer  :: couple_iters_last = 0      !! iterations taken last commit (diagnostic)
      ! Co-convergence state for the SLE driver's σ<->τ fixed point (§3c 3b). Each
      ! advance_endpoint does ONE trapezoid pass and refreshes the report drift to
      ! the new τ_{n+1}; the driver re-converges σ against it and repeats until done.
      logical  :: couple_done = .false.      !! co-convergence reached (drift settled)
      integer  :: couple_pass = 0            !! co-convergence passes taken this step
      ! Scratch for the implicit commit (allocated lazily on first TRAP commit):
      real(wp), allocatable :: Are0(:,:,:), Aim0(:,:,:)  !! (NLAM,ne,nk) τ_n snapshot
      real(wp), allocatable :: Bre0(:,:,:), Bim0(:,:,:)
      real(wp), allocatable :: Cre0(:,:,:), Cim0(:,:,:)
      real(wp), allocatable :: edUn_re(:,:), edUn_im(:,:) !! (nr,nk) endpoint nodal drift U (ε_{n+1})
      real(wp), allocatable :: edVn_re(:,:), edVn_im(:,:) !! (nr,nk) endpoint nodal drift V
      complex(wp), allocatable :: dUa_prev(:)             !! (nk) prev-iterate surface drift (convergence)
      ! Start-of-step load σ_n for the trapezoidal ε_n term. The rule is ½(ε_n+ε_{n+1})
      ! with ε_n = σ_n·xUn + drift(τ_n) and ε_{n+1} = σ_{n+1}·xUn + drift(τ_{n+1}); ε_n
      ! must use the load at t_n, not the current σ_{n+1}. For a held load σ_n=σ_{n+1}
      ! so this is invisible (the historical/3a path), but for a fast-evolving load
      ! (the SLE 3b driver) it sets the order. sigma_n carries the previous step's
      ! converged load; sigma_next stages σ_{n+1} until finalize commits it. Until the
      ! first step finalizes (sigma_primed=.false.) ε_n falls back to σ_{n+1}, which
      ! reproduces the historical first step exactly for a load present at t=0.
      complex(wp), allocatable :: sigma_n(:), sigma_next(:)  !! (nlm)
      logical :: sigma_primed = .false.
      ! Δt-invariant memory rate Mk/Δt = μ/η (set once in init): set_dt rescales the
      ! memory factor Mk = MkPerDt·Δt exactly, so the adaptive controller can change
      ! Δt with no operator re-factor and no drift from repeated rescaling.
      real(wp), allocatable :: MkPerDt(:)                 !! (ne)
      ! Controller state snapshots (§3c controller, lazily allocated): buffer A holds
      ! the entering state τ_n that a rejected/fine step restores to; buffer B holds the
      ! coarse τ_{n+1} for the step-doubling error estimate. Distinct from Are0 (which
      ! sle_solve overwrites internally each step). See save_state/stash_coarse.
      real(wp), allocatable :: Are_s(:,:,:), Aim_s(:,:,:), Bre_s(:,:,:), Bim_s(:,:,:)
      real(wp), allocatable :: Cre_s(:,:,:), Cim_s(:,:,:)            !! buffer A (τ_n)
      real(wp), allocatable :: Are_c(:,:,:), Aim_c(:,:,:), Bre_c(:,:,:), Bim_c(:,:,:)
      real(wp), allocatable :: Cre_c(:,:,:), Cim_c(:,:,:)            !! buffer B (τ_coarse)
      real(wp)    :: time_s = 0.0_wp                       !! saved time (buffer A)
      complex(wp), allocatable :: sigma_n_s(:)             !! saved σ_n (buffer A)
      logical     :: sigma_primed_s = .false.
   contains
      procedure :: init        => ve_response_init
      procedure :: begin_step  => ve_response_begin
      procedure :: apply       => ve_response_apply
      procedure :: horizontal  => ve_response_horizontal
      procedure :: commit_step => ve_response_commit
      procedure :: prepare_endpoint   => ve_response_prepare_endpoint
      procedure :: advance_endpoint   => ve_response_advance_endpoint
      procedure :: endpoint_converged => ve_response_endpoint_converged
      procedure :: finalize_step      => ve_response_finalize_step
      ! Adaptive-Δt controller primitives (§3c controller; see fe_timestep)
      procedure :: prime_sigma       => ve_response_prime_sigma
      procedure :: set_dt            => ve_response_set_dt
      procedure :: save_state        => ve_response_save_state
      procedure :: restore_state     => ve_response_restore_state
      procedure :: stash_coarse      => ve_response_stash_coarse
      procedure :: coarse_fine_error => ve_response_coarse_fine_error
      procedure :: destroy     => ve_response_destroy
   end type ve_response

contains

   subroutine response_begin_default(self, sht)
      !! No-op step bracket for stateless responses.
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_begin_default

   subroutine response_commit_default(self, sht, sigma_lm)
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_commit_default

   subroutine response_prepare_default(self, sht)
      !! No-op endpoint bracket for stateless / 1st-order responses.
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_prepare_default

   subroutine response_advance_default(self, sht, sigma_lm)
      !! No memory to advance (stateless response): nothing to do.
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_advance_default

   logical function response_converged_default(self) result(done)
      !! Stateless / 1st-order responses converge in a single pass.
      class(response_operator), intent(in) :: self
      done = .true.
   end function response_converged_default

   subroutine response_finalize_default(self, sht)
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_finalize_default

   subroutine response_horizontal_default(self, sht, sigma_lm, v_lm)
      !! No horizontal displacement for a rigid / non-deforming response.
      class(response_operator), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)  !! load [kg m^-2]
      complex(wp),              intent(out)   :: v_lm(:)      !! spheroidal V(a) [m]
      v_lm = (0.0_wp, 0.0_wp)
   end subroutine response_horizontal_default

   subroutine null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      class(null_response), intent(inout) :: self
      type(sht_grid),       intent(in)    :: sht
      complex(wp),          intent(in)    :: sigma_lm(:)
      complex(wp),          intent(out)   :: u_lm(:)
      complex(wp),          intent(out)   :: n_lm(:)
      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
   end subroutine null_response_apply

   subroutine elastic_response_init(self, earth, lmax)
      !! Precompute the per-degree elastic surface gains for degrees 0..lmax.
      !!
      !!   l = 0 : incompressibility (Div u = 0) forbids degree-0 radial
      !!           deformation, so U(0)=0; the geoid feels only the load's own
      !!           monopole potential, N(0) = φ^L_0/g = 4πGa/g per unit σ.
      !!   l ≥ 1 : assemble the per-degree saddle-point operator, solve a unit
      !!           surface load, store U(a) and N(a) = −F(a)/g.
      class(elastic_response), intent(inout) :: self
      type(earth_model),       intent(in)    :: earth
      integer,                 intent(in)    :: lmax
      type(radial_mesh)     :: mesh
      type(radial_operator) :: op
      integer  :: l
      real(wp) :: ua, va, fa

      call self%destroy()
      self%lmax = lmax
      self%a    = earth%r_earth
      self%g    = earth%gravity_at(earth%r_earth)
      allocate(self%ugain(0:lmax), self%ngain(0:lmax), self%vgain(0:lmax))

      ! degree 0: no deformation, pure monopole geoid offset
      self%ugain(0) = 0.0_wp
      self%ngain(0) = 4.0_wp*pi*grav_G*self%a / self%g
      self%vgain(0) = 0.0_wp                    ! no horizontal at degree 0

      call mesh%build(earth)
      do l = 1, lmax
         call op%assemble(earth, mesh, l)
         call op%solve(1.0_wp, ua, va, fa)     ! unit surface load coefficient
         self%ugain(l) = ua
         self%ngain(l) = -fa / self%g
         self%vgain(l) = va
         call op%destroy()
      end do

      ! degree-1 geoid frame: the per-degree solve fixes the displacement gauge
      ! (wᵀd=0, geocenter/CE-like, h₁≈0), but the geoid (sea surface) is referenced
      ! to the CM frame, in which the degree-1 external potential vanishes ⇒ N₁≡0.
      ! (The benchmark M3-L70-V01 table has k₁=−1 exactly, i.e. N₁=(1+k₁)φ^L/g=0;
      ! validated against the Spada-2011 disc n_disc, which matches once N₁ is
      ! dropped.) Displacement degree-1 (ugain(1)) is left as solved.
      if (lmax >= 1) self%ngain(1) = 0.0_wp
   end subroutine elastic_response_init

   subroutine elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Spectral multiply: u_lm = ugain(l)·σ_lm, n_lm = ngain(l)·σ_lm. Degrees
      !! above the precomputed lmax are zeroed.
      class(elastic_response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: u_lm(:)
      complex(wp),             intent(out)   :: n_lm(:)
      integer :: l, m, lm, lcap

      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht%lmidx(l, m)
            u_lm(lm) = self%ugain(l) * sigma_lm(lm)
            n_lm(lm) = self%ngain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_apply

   subroutine elastic_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal multiply: v_lm = vgain(l)·σ_lm (degree-1 left as solved, like
      !! ugain — the horizontal displacement is in the CE-like gauge, not the geoid
      !! CM frame). Synthesize ∇₁(Σ v_lm Y_lm) for (u_θ, u_φ).
      class(elastic_response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: v_lm(:)
      integer :: l, m, lm, lcap
      v_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht%lmidx(l, m)
            v_lm(lm) = self%vgain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_horizontal

   subroutine elastic_response_destroy(self)
      class(elastic_response), intent(inout) :: self
      if (allocated(self%ugain)) deallocate(self%ugain)
      if (allocated(self%ngain)) deallocate(self%ngain)
      if (allocated(self%vgain)) deallocate(self%vgain)
      self%lmax = 0
   end subroutine elastic_response_destroy

   ! --- viscoelastic field driver ---------------------------------------------

   subroutine ve_response_init(self, earth, sht, dt)
      !! Assemble the per-degree operators, precompute the unit-load response and
      !! Maxwell constants, and zero the per-(l,m) memory. Tied to the grid sht
      !! (sets lmax = sht%lmax and the coefficient layout).
      class(ve_response), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: dt
      type(radial_mesh) :: mesh
      real(wp), allocatable :: x(:)
      real(wp) :: eta_e
      integer  :: l, m, e, lay, node, k

      call self%destroy()
      call mesh%build(earth)
      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%dt = dt;  self%time = 0.0_wp
      self%a  = earth%r_earth;  self%g = earth%gravity_at(earth%r_earth)

      ! degree-independent element fields: node radii, shear, Maxwell factor
      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne), self%MkPerDt(self%ne))
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         if (eta_e > 0.0_wp) then
            self%Mk(e)      = self%mu(e)*dt/eta_e     ! historical form (byte-identical)
            self%MkPerDt(e) = self%mu(e)/eta_e        ! Δt-invariant; set_dt rescales Mk
         else
            self%Mk(e)      = 0.0_wp                  ! elastic layer: no memory
            self%MkPerDt(e) = 0.0_wp
         end if
      end do

      ! per-degree constants, operators, and unit-load response
      allocate(self%Jr(self%lmax), self%nrmc(NLAM,self%lmax))
      allocate(self%sa(4,NLAM,self%lmax), self%sb(4,NLAM,self%lmax), &
               self%sc(4,NLAM,self%lmax))
      allocate(self%gu(0:self%lmax), self%gn(0:self%lmax))
      allocate(self%xUn(self%nr,self%lmax), self%xVn(self%nr,self%lmax))
      allocate(self%gv(0:self%lmax))
      allocate(self%ops(self%lmax), x(self%ndof))

      ! degree 0: monopole geoid, no deformation, no memory (no operator).
      ! degree 1: geocenter motion — carried in the CM frame. The sparse KKT
      ! border in radial_operator removes the rigid-translation null space
      ! (wᵀd = 0, Blewitt 2003), so j=1 assembles and steps like any other
      ! degree; it joins the l-loop below.
      self%gu(0) = 0.0_wp
      self%gn(0) = 4.0_wp*pi*grav_G*self%a/self%g
      self%gv(0) = 0.0_wp                       ! no horizontal at degree 0

      do l = 1, self%lmax
         self%Jr(l) = real(l, wp)*real(l+1, wp)
         call ve_strain_constants(self%Jr(l), self%nrmc(:,l), &
                                  self%sa(:,:,l), self%sb(:,:,l), self%sc(:,:,l))
         call self%ops(l)%assemble(earth, mesh, l)
         call self%ops(l)%solve_vec(self%ops(l)%load_rhs(1.0_wp), x)
         self%gu(l) = x(idx_u(self%nr))
         self%gn(l) = -x(idx_f(self%nr))/self%g
         self%gv(l) = x(idx_v(self%nr))         ! surface horizontal V(a)
         do node = 1, self%nr
            self%xUn(node,l) = x(idx_u(node))
            self%xVn(node,l) = x(idx_v(node))
         end do
      end do

      ! degree-1 geoid frame (see elastic_response_init): the geoid is referenced
      ! to the CM frame ⇒ N₁≡0. Zero the degree-1 geoid gain here; the degree-1
      ! relaxation drift is likewise zeroed in begin_step. Displacement (gu(1),
      ! xUn/xVn) is left as solved (CE-like geocenter, h₁≈0).
      if (self%lmax >= 1) self%gn(1) = 0.0_wp

      ! Degree-grouped coefficient map: slot k = 1..nk over (l>=1, m=0..min(l,mmax)),
      ! l ascending then m ascending. k2lm bridges back to the SHTns lm index for
      ! the load/uplift/geoid spectra; kdeg gives the degree (operator) per slot.
      self%nk = 0
      do l = 1, self%lmax
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            self%nk = self%nk + 1
         end do
      end do
      allocate(self%k2lm(self%nk), self%kdeg(self%nk), self%kbeg(self%lmax+1))
      k = 0
      do l = 1, self%lmax
         self%kbeg(l) = k + 1                  ! first slot of degree l (contiguous)
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            k = k + 1
            self%k2lm(k) = sht%lmidx(l, m)
            self%kdeg(k) = l
         end do
      end do
      self%kbeg(self%lmax+1) = k + 1           ! sentinel (one past the last slot)

      ! per-(l,m) memory + drift, all in degree-grouped k order (zeroed)
      allocate(self%Are(NLAM,self%ne,self%nk), self%Aim(NLAM,self%ne,self%nk))
      allocate(self%Bre(NLAM,self%ne,self%nk), self%Bim(NLAM,self%ne,self%nk))
      allocate(self%Cre(NLAM,self%ne,self%nk), self%Cim(NLAM,self%ne,self%nk))
      self%Are = 0.0_wp; self%Aim = 0.0_wp; self%Bre = 0.0_wp
      self%Bim = 0.0_wp; self%Cre = 0.0_wp; self%Cim = 0.0_wp
      allocate(self%dUa(self%nk), self%dFa(self%nk), self%dVa(self%nk))
      allocate(self%dUn_re(self%nr,self%nk), self%dUn_im(self%nr,self%nk))
      allocate(self%dVn_re(self%nr,self%nk), self%dVn_im(self%nr,self%nk))
      allocate(self%mnorm(self%nk))
      self%dUa = (0.0_wp,0.0_wp); self%dFa = (0.0_wp,0.0_wp); self%dVa = (0.0_wp,0.0_wp)
      self%dUn_re = 0.0_wp; self%dUn_im = 0.0_wp
      self%dVn_re = 0.0_wp; self%dVn_im = 0.0_wp
      self%mnorm = 0.0_wp                       ! zero memory ⇒ all slots skipped initially
   end subroutine ve_response_init

   subroutine ve_response_begin(self, sht)
      !! Freeze the per-(l,m) drift from the entering memory τ_n: solve the memory
      !! forcing −∫τ^V:δε with the load held at zero, storing surface + nodal drift.
      !! The nodal drift (the ε_n term) lands in self%dUn_*/dVn_*; the implicit
      !! commit re-uses the same solver against the trial τ_{n+1} (see solve_drift).
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      call solve_drift(self, sht, self%dUn_re, self%dUn_im, self%dVn_re, self%dVn_im)
   end subroutine ve_response_begin

   subroutine solve_drift(self, sht, Un_re, Un_im, Vn_re, Vn_im)
      !! Solve the per-(l,m) drift (load=0 memory forcing) from self's CURRENT memory
      !! arrays (self%Are…): surface drift → self%dUa/dFa/dVa, nodal drift → the four
      !! target arrays. begin_step passes self%dUn_* (drift from τ_n); the implicit
      !! commit passes self%edUn_* (drift from the trial τ_{n+1}). Targets are distinct
      !! components from everything read via self, so there is no argument aliasing.
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(out)   :: Un_re(:,:), Un_im(:,:), Vn_re(:,:), Vn_im(:,:)
      real(wp), allocatable :: fre(:), fim(:), xre(:), xim(:)
      integer :: k, l, node, e, mm
      real(wp) :: thr, mk

      ! Refresh the per-slot memory magnitude from the current memory (a cheap pass
      ! vs the solves) so it is always consistent — including after a restart, which
      ! reloads the memory arrays but not this derived cache. Explicit loop (no
      ! abs(slice) temporaries, which would each heap-allocate).
      !$omp parallel do default(shared) private(k, e, mm, mk) schedule(static)
      do k = 1, self%nk
         mk = 0.0_wp
         do e = 1, self%ne
            do mm = 1, NLAM
               mk = max(mk, abs(self%Are(mm,e,k)), abs(self%Bre(mm,e,k)), &
                           abs(self%Cre(mm,e,k)), abs(self%Aim(mm,e,k)), &
                           abs(self%Bim(mm,e,k)), abs(self%Cim(mm,e,k)))
            end do
         end do
         self%mnorm(k) = mk
      end do
      !$omp end parallel do
      ! Skip the drift solve for coefficients with negligible memory (their drift is
      ! negligible too): zero their drift instead of solving. thr is relative to the
      ! largest memory present.
      thr = self%skip_tol * maxval(self%mnorm)

      ! Solve for the drift, PARALLEL OVER DEGREE l so each per-degree operator
      ! ops(l) is touched by a single thread. Safe because j>=2 uses the re-entrant
      ! banded LU (fe_band); degree 1 (the lone LIS solver, not re-entrant) is a
      ! single iteration, hence run by a single thread — no concurrent LIS call. The
      ! scratch vectors are per-thread; dynamic schedule balances the rising work
      ! per degree (l+1 orders). Inactive (ordinary serial loop) unless openmp=1.
      !$omp parallel default(shared) private(l, k, node, fre, fim, xre, xim)
      allocate(fre(self%ndof), fim(self%ndof), xre(self%ndof), xim(self%ndof))
      !$omp do schedule(dynamic)
      do l = 1, self%lmax
         do k = self%kbeg(l), self%kbeg(l+1) - 1
            if (self%mnorm(k) <= thr) then       ! negligible memory ⇒ negligible drift
               self%dUa(k) = (0.0_wp,0.0_wp);  self%dFa(k) = (0.0_wp,0.0_wp)
               self%dVa(k) = (0.0_wp,0.0_wp)
               Un_re(:,k) = 0.0_wp;  Un_im(:,k) = 0.0_wp
               Vn_re(:,k) = 0.0_wp;  Vn_im(:,k) = 0.0_wp
               cycle
            end if
            fre = 0.0_wp;  fim = 0.0_wp
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Are(:,:,k), self%Bre(:,:,k), &
                 self%Cre(:,:,k), fre)
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Aim(:,:,k), self%Bim(:,:,k), &
                 self%Cim(:,:,k), fim)
            call self%ops(l)%solve_vec(fre, xre)
            call self%ops(l)%solve_vec(fim, xim)
            self%dUa(k) = cmplx(xre(idx_u(self%nr)), xim(idx_u(self%nr)), wp)
            self%dFa(k) = cmplx(xre(idx_f(self%nr)), xim(idx_f(self%nr)), wp)
            self%dVa(k) = cmplx(xre(idx_v(self%nr)), xim(idx_v(self%nr)), wp)
            if (l == 1) self%dFa(k) = (0.0_wp, 0.0_wp)   ! N₁≡0 (CM frame; see init)
            do node = 1, self%nr
               Un_re(node,k) = xre(idx_u(node))
               Un_im(node,k) = xim(idx_u(node))
               Vn_re(node,k) = xre(idx_v(node))
               Vn_im(node,k) = xim(idx_v(node))
            end do
         end do
      end do
      !$omp end do
      deallocate(fre, fim, xre, xim)
      !$omp end parallel
   end subroutine solve_drift

   subroutine ve_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Affine response at the frozen time: u = gu(l)·σ + drift_U,
      !! N = gn(l)·σ − drift_F/g.
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: u_lm(:)
      complex(wp),        intent(out)   :: n_lm(:)
      integer :: k, l, lm, lm0

      u_lm = (0.0_wp,0.0_wp);  n_lm = (0.0_wp,0.0_wp)
      ! degree 0: monopole geoid, no deformation, no memory
      lm0 = sht%lmidx(0, 0)
      u_lm(lm0) = self%gu(0)*sigma_lm(lm0)
      n_lm(lm0) = self%gn(0)*sigma_lm(lm0)
      ! degrees l>=1, in degree-grouped k order (gn(1)=0 and dFa(k)=0 give N₁≡0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         u_lm(lm) = self%gu(l)*sigma_lm(lm) + self%dUa(k)
         n_lm(lm) = self%gn(l)*sigma_lm(lm) - self%dFa(k)/self%g
      end do
   end subroutine ve_response_apply

   subroutine ve_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal V at the frozen time: v_lm = gv(l)·σ + drift_V. Uses the same
      !! frozen drift (dVa) as the last begin_step, so calling it after a converged
      !! step gives the horizontal consistent with apply()'s u/N. Degree 1 left as
      !! solved (CE-like gauge, like u — not the geoid CM frame).
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: v_lm(:)
      integer :: k, l, lm
      v_lm = (0.0_wp,0.0_wp)               ! degree 0 has no horizontal (gv(0)=0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         v_lm(lm) = self%gv(l)*sigma_lm(lm) + self%dVa(k)
      end do
   end subroutine ve_response_horizontal

   subroutine ve_response_commit(self, sht, sigma_lm)
      !! Advance the memory with the converged load, frozen σ (held/slow-load step;
      !! §3c part 3a). Explicit (FE): total nodal strain = σ·(unit-load nodal) +
      !! drift(τ_n), one Maxwell update per (l,m). Implicit (TRAP): the endpoint is
      !! solved by Picard iteration — re-solve the drift against the trial τ_{n+1},
      !! form the endpoint strain, trapezoid-advance from τ_n, repeat to couple_tol.
      !! Advances time by Δt. For fast-evolving loads the SLE driver instead iterates
      !! prepare_endpoint/advance_endpoint/finalize_step so σ co-converges (3b).
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm
      integer  :: iter

      if (.not. scheme_is_implicit(self%scheme)) then
         call fe_advance(self, sigma_lm)            ! explicit: byte-for-byte historical
         self%couple_iters_last = 1
         self%time = self%time + self%dt
         return
      end if

      ! --- implicit (TRAP): iterate the endpoint to a consistent τ_{n+1} ------------
      call ensure_commit_scratch(self)
      call snapshot_taun(self)                      ! τ_n base for every trapezoid pass
      self%dUa_prev = self%dUa
      do iter = 1, self%max_couple_iter
         ! Endpoint drift from the current τ_{n+1} estimate (self%Are…); also refreshes
         ! self%dUa (the surface drift, used as the convergence signal). Then reset to
         ! τ_n and trapezoid-advance with (ε_n, ε_{n+1}).
         call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im)
         call trapezoid_advance_all(self, sigma_lm)
         ! iter 1 re-solves drift from τ_n and so reproduces begin_step's drift; the
         ! fixed point only moves at iter 2 (never exit on the first pass).
         cnorm = maxval(abs(self%dUa - self%dUa_prev))
         snorm = maxval(abs(self%dUa))
         self%dUa_prev = self%dUa
         if (iter >= 2 .and. cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp))) exit
      end do
      self%couple_iters_last = min(iter, self%max_couple_iter)
      self%sigma_n = sigma_lm;  self%sigma_primed = .true.   ! σ_n for the next step's ε_n
      self%time = self%time + self%dt
   end subroutine ve_response_commit

   subroutine fe_advance(self, sigma_lm)
      !! Explicit forward-Euler memory advance: one Maxwell update per (l,m) from the
      !! report strain σ·(unit-load nodal) + drift(τ_n). The historical path, shared by
      !! commit_step and advance_endpoint so both stay byte-identical for FE.
      class(ve_response), intent(inout) :: self
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp) :: sre, sim
      integer  :: k, l, lm, node
      !$omp parallel default(shared) private(k, l, lm, node, sre, sim, Ure, Uim, Vre, Vim)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk                          ! degree-grouped order
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))
         do node = 1, self%nr
            Ure(node) = sre*self%xUn(node,l) + self%dUn_re(node,k)
            Uim(node) = sim*self%xUn(node,l) + self%dUn_im(node,k)
            Vre(node) = sre*self%xVn(node,l) + self%dVn_re(node,k)
            Vim(node) = sim*self%xVn(node,l) + self%dVn_im(node,k)
         end do
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k))
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k))
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim)
      !$omp end parallel
   end subroutine fe_advance

   subroutine snapshot_taun(self)
      !! Snapshot the entering memory τ_n into the *0 arrays so every trapezoid pass
      !! advances from τ_n (not compounding). ε_n nodal drift is in self%dUn_* (set by
      !! begin_step) and stays fixed; ε_{n+1} drift is re-solved into self%edUn_*.
      class(ve_response), intent(inout) :: self
      self%Are0 = self%Are;  self%Aim0 = self%Aim
      self%Bre0 = self%Bre;  self%Bim0 = self%Bim
      self%Cre0 = self%Cre;  self%Cim0 = self%Cim
   end subroutine snapshot_taun

   subroutine trapezoid_advance_all(self, sigma_lm)
      !! One trapezoid endpoint advance for all (l,m): reset memory to τ_n (the *0
      !! snapshot), then advance with the entering strain ε_n (σ·xUn + dUn) and the
      !! endpoint strain ε_{n+1} (σ·xUn + edUn). Reads self%edUn_*/edVn_* (the current
      !! τ_{n+1} drift estimate); writes self%Are…. Does not touch self%dUa.
      class(ve_response), intent(inout) :: self
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp), allocatable :: Ure_n(:), Uim_n(:), Vre_n(:), Vim_n(:)
      real(wp) :: sre, sim, srn, sin
      integer  :: k, l, lm, node
      !$omp parallel default(shared) &
      !$omp   private(k, l, lm, node, sre, sim, srn, sin, Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr), &
               Ure_n(self%nr), Uim_n(self%nr), Vre_n(self%nr), Vim_n(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))          ! σ_{n+1}
         if (self%sigma_primed) then
            srn = real(self%sigma_n(lm), wp);  sin = aimag(self%sigma_n(lm))   ! σ_n
         else
            srn = sre;  sin = sim     ! first step: ε_n uses σ_{n+1} (load present at t=0)
         end if
         do node = 1, self%nr
            Ure_n(node) = srn*self%xUn(node,l) + self%dUn_re(node,k)   ! ε_n (σ_n)
            Vre_n(node) = srn*self%xVn(node,l) + self%dVn_re(node,k)
            Uim_n(node) = sin*self%xUn(node,l) + self%dUn_im(node,k)
            Vim_n(node) = sin*self%xVn(node,l) + self%dVn_im(node,k)
            Ure(node)   = sre*self%xUn(node,l) + self%edUn_re(node,k)  ! ε_{n+1}
            Vre(node)   = sre*self%xVn(node,l) + self%edVn_re(node,k)
            Uim(node)   = sim*self%xUn(node,l) + self%edUn_im(node,k)
            Vim(node)   = sim*self%xVn(node,l) + self%edVn_im(node,k)
         end do
         ! reset to τ_n, then trapezoid-advance with (ε_n, ε_{n+1})
         self%Are(:,:,k) = self%Are0(:,:,k);  self%Bre(:,:,k) = self%Bre0(:,:,k)
         self%Cre(:,:,k) = self%Cre0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Ure_n, Vn_prev=Vre_n)
         self%Aim(:,:,k) = self%Aim0(:,:,k);  self%Bim(:,:,k) = self%Bim0(:,:,k)
         self%Cim(:,:,k) = self%Cim0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Uim_n, Vn_prev=Vim_n)
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      !$omp end parallel
   end subroutine trapezoid_advance_all

   subroutine ve_response_prepare_endpoint(self, sht)
      !! Open a co-converging step (§3c 3b): snapshot τ_n and seed the endpoint drift
      !! ε_{n+1} with the entering τ_n drift (begin_step's dUn). The SLE driver then
      !! converges σ against the current report drift (dUa, = drift(τ_n) on entry) and
      !! calls advance_endpoint, which refreshes the report drift to τ_{n+1}.
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_pass = 0
      self%couple_done = .false.
      if (.not. scheme_is_implicit(self%scheme)) return    ! FE: nothing to snapshot
      call ensure_commit_scratch(self)
      call snapshot_taun(self)
      self%edUn_re = self%dUn_re;  self%edUn_im = self%dUn_im
      self%edVn_re = self%dVn_re;  self%edVn_im = self%dVn_im
      self%dUa_prev = self%dUa
   end subroutine ve_response_prepare_endpoint

   subroutine ve_response_advance_endpoint(self, sht, sigma_lm)
      !! One co-convergence pass with the SLE-converged load σ (§3c 3b). FE: a single
      !! Maxwell update (1st-order; no co-iteration). TRAP: trapezoid-advance τ_n→τ_{n+1}
      !! using the current endpoint-drift estimate, THEN refresh the report drift dUa
      !! (and ε_{n+1}=edUn) from the new τ_{n+1} — so the driver's next σ-convergence,
      !! and the next pass's endpoint strain, see the advanced memory. Sets couple_done
      !! when the surface drift settles to couple_tol. Does NOT advance time.
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm

      if (.not. scheme_is_implicit(self%scheme)) then
         call fe_advance(self, sigma_lm)
         self%couple_pass = 1
         self%couple_done = .true.
         return
      end if

      call trapezoid_advance_all(self, sigma_lm)   ! ε_{n+1} from the previous estimate
      self%sigma_next = sigma_lm                   ! stage σ_{n+1}; finalize commits it to σ_n
      ! refresh dUa + ε_{n+1} drift from the new τ_{n+1} (Are…), ready for the next pass
      call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im)
      self%couple_pass = self%couple_pass + 1
      cnorm = maxval(abs(self%dUa - self%dUa_prev))
      snorm = maxval(abs(self%dUa))
      self%dUa_prev = self%dUa
      self%couple_done = (self%couple_pass >= 2 .and. &
                          cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp)))
   end subroutine ve_response_advance_endpoint

   logical function ve_response_endpoint_converged(self) result(done)
      class(ve_response), intent(in) :: self
      done = self%couple_done
   end function ve_response_endpoint_converged

   subroutine ve_response_finalize_step(self, sht)
      !! Close a co-converging step: the memory already holds the converged τ_{n+1}
      !! (advance_endpoint left it there), so only advance time.
      class(ve_response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_iters_last = self%couple_pass
      if (allocated(self%sigma_next)) then          ! TRAP: commit σ_{n+1} as next σ_n
         self%sigma_n = self%sigma_next;  self%sigma_primed = .true.
      end if
      self%time = self%time + self%dt
   end subroutine ve_response_finalize_step

   subroutine ve_response_set_dt(self, dt)
      !! Change the step size. Δt enters only through Mk = (μ/η)·Δt, so rescale Mk from
      !! the Δt-invariant rate MkPerDt — exact (no drift from repeated halving/restoring)
      !! and no operator re-factor (the band LU is Δt-independent). The adaptive
      !! controller uses this to try a step, halve for the fine sub-steps, and restore.
      class(ve_response), intent(inout) :: self
      real(wp),           intent(in)    :: dt
      self%dt = dt
      self%Mk = self%MkPerDt * dt
   end subroutine ve_response_set_dt

   subroutine ensure_state_scratch(self)
      !! Lazily allocate the controller's state buffers A (τ_n) and B (τ_coarse).
      class(ve_response), intent(inout) :: self
      if (allocated(self%Are_s)) return
      allocate(self%Are_s(NLAM,self%ne,self%nk), self%Aim_s(NLAM,self%ne,self%nk))
      allocate(self%Bre_s(NLAM,self%ne,self%nk), self%Bim_s(NLAM,self%ne,self%nk))
      allocate(self%Cre_s(NLAM,self%ne,self%nk), self%Cim_s(NLAM,self%ne,self%nk))
      allocate(self%Are_c(NLAM,self%ne,self%nk), self%Aim_c(NLAM,self%ne,self%nk))
      allocate(self%Bre_c(NLAM,self%ne,self%nk), self%Bim_c(NLAM,self%ne,self%nk))
      allocate(self%Cre_c(NLAM,self%ne,self%nk), self%Cim_c(NLAM,self%ne,self%nk))
      allocate(self%sigma_n_s(self%nlm))
   end subroutine ensure_state_scratch

   subroutine ve_response_save_state(self)
      !! Snapshot the entering prognostic state (memory τ_n + time + σ_n) into buffer A.
      !! A rejected step or the fine sub-step path restores to this with restore_state.
      class(ve_response), intent(inout) :: self
      call ensure_state_scratch(self)
      self%Are_s = self%Are;  self%Aim_s = self%Aim
      self%Bre_s = self%Bre;  self%Bim_s = self%Bim
      self%Cre_s = self%Cre;  self%Cim_s = self%Cim
      self%time_s = self%time
      if (allocated(self%sigma_n)) self%sigma_n_s = self%sigma_n
      self%sigma_primed_s = self%sigma_primed
   end subroutine ve_response_save_state

   subroutine ve_response_restore_state(self)
      !! Restore the prognostic state saved by save_state (buffer A).
      class(ve_response), intent(inout) :: self
      self%Are = self%Are_s;  self%Aim = self%Aim_s
      self%Bre = self%Bre_s;  self%Bim = self%Bim_s
      self%Cre = self%Cre_s;  self%Cim = self%Cim_s
      self%time = self%time_s
      if (allocated(self%sigma_n)) self%sigma_n = self%sigma_n_s
      self%sigma_primed = self%sigma_primed_s
   end subroutine ve_response_restore_state

   subroutine ve_response_stash_coarse(self)
      !! Snapshot the current memory (the coarse one-Δt τ_{n+1}) into buffer B for the
      !! step-doubling error estimate, to be compared against the fine result.
      class(ve_response), intent(inout) :: self
      call ensure_state_scratch(self)
      self%Are_c = self%Are;  self%Aim_c = self%Aim
      self%Bre_c = self%Bre;  self%Bim_c = self%Bim
      self%Cre_c = self%Cre;  self%Cim_c = self%Cim
   end subroutine ve_response_stash_coarse

   subroutine ve_response_coarse_fine_error(self, err_inf, tau_inf)
      !! After the fine path, return the coarse↔fine memory difference ‖τ_fine−τ_coarse‖∞
      !! (buffer B is τ_coarse, self%Are… is τ_fine) and the memory magnitude ‖τ_fine‖∞,
      !! for the controller's scaled local-error estimate.
      class(ve_response), intent(in)  :: self
      real(wp),           intent(out) :: err_inf, tau_inf
      err_inf = max(maxval(abs(self%Are - self%Are_c)), maxval(abs(self%Aim - self%Aim_c)), &
                    maxval(abs(self%Bre - self%Bre_c)), maxval(abs(self%Bim - self%Bim_c)), &
                    maxval(abs(self%Cre - self%Cre_c)), maxval(abs(self%Cim - self%Cim_c)))
      tau_inf = max(maxval(abs(self%Are)), maxval(abs(self%Aim)), &
                    maxval(abs(self%Bre)), maxval(abs(self%Bim)), &
                    maxval(abs(self%Cre)), maxval(abs(self%Cim)))
   end subroutine ve_response_coarse_fine_error

   subroutine ensure_commit_scratch(self)
      !! Lazily allocate the implicit-commit scratch (the τ_n memory snapshot doubles
      !! the memory footprint, so it is only paid when a TRAP commit is first used).
      class(ve_response), intent(inout) :: self
      if (allocated(self%Are0)) return
      allocate(self%Are0(NLAM,self%ne,self%nk), self%Aim0(NLAM,self%ne,self%nk))
      allocate(self%Bre0(NLAM,self%ne,self%nk), self%Bim0(NLAM,self%ne,self%nk))
      allocate(self%Cre0(NLAM,self%ne,self%nk), self%Cim0(NLAM,self%ne,self%nk))
      allocate(self%edUn_re(self%nr,self%nk), self%edUn_im(self%nr,self%nk))
      allocate(self%edVn_re(self%nr,self%nk), self%edVn_im(self%nr,self%nk))
      allocate(self%dUa_prev(self%nk))
      call ensure_sigma(self)
   end subroutine ensure_commit_scratch

   subroutine ensure_sigma(self)
      !! Lazily allocate the start-of-step load buffers (σ_n / σ_next). Separate from
      !! ensure_commit_scratch so prime_sigma can seed σ_0 before the first commit.
      class(ve_response), intent(inout) :: self
      if (allocated(self%sigma_n)) return
      allocate(self%sigma_n(self%nlm), self%sigma_next(self%nlm))
      self%sigma_n = (0.0_wp, 0.0_wp)     ! σ at t=0; primed by prime_sigma or step 1
      self%sigma_primed = .false.
   end subroutine ensure_sigma

   subroutine ve_response_prime_sigma(self, sigma_lm)
      !! Seed the start-of-step load σ_n with a known load (the elastic-consistent SLE
      !! load at t=0) and mark it tracked, so the trapezoidal ε_n on the FIRST step uses
      !! the true σ_0 rather than the σ_{n+1} proxy — making that step 2nd→3rd order.
      class(ve_response), intent(inout) :: self
      complex(wp),        intent(in)    :: sigma_lm(:)
      call ensure_sigma(self)
      self%sigma_n = sigma_lm
      self%sigma_primed = .true.
   end subroutine ve_response_prime_sigma

   subroutine ve_response_destroy(self)
      class(ve_response), intent(inout) :: self
      integer :: l
      if (allocated(self%ops)) then
         do l = 1, size(self%ops);  call self%ops(l)%destroy();  end do
         deallocate(self%ops)
      end if
      if (allocated(self%r))      deallocate(self%r)
      if (allocated(self%mu))     deallocate(self%mu)
      if (allocated(self%Mk))     deallocate(self%Mk)
      if (allocated(self%Jr))     deallocate(self%Jr)
      if (allocated(self%nrmc))   deallocate(self%nrmc)
      if (allocated(self%sa))     deallocate(self%sa)
      if (allocated(self%sb))     deallocate(self%sb)
      if (allocated(self%sc))     deallocate(self%sc)
      if (allocated(self%gu))     deallocate(self%gu)
      if (allocated(self%gn))     deallocate(self%gn)
      if (allocated(self%gv))     deallocate(self%gv)
      if (allocated(self%xUn))    deallocate(self%xUn)
      if (allocated(self%xVn))    deallocate(self%xVn)
      if (allocated(self%Are))    deallocate(self%Are)
      if (allocated(self%Aim))    deallocate(self%Aim)
      if (allocated(self%Bre))    deallocate(self%Bre)
      if (allocated(self%Bim))    deallocate(self%Bim)
      if (allocated(self%Cre))    deallocate(self%Cre)
      if (allocated(self%Cim))    deallocate(self%Cim)
      if (allocated(self%dUa))    deallocate(self%dUa)
      if (allocated(self%dFa))    deallocate(self%dFa)
      if (allocated(self%dVa))    deallocate(self%dVa)
      if (allocated(self%dUn_re)) deallocate(self%dUn_re)
      if (allocated(self%dUn_im)) deallocate(self%dUn_im)
      if (allocated(self%dVn_re)) deallocate(self%dVn_re)
      if (allocated(self%dVn_im)) deallocate(self%dVn_im)
      if (allocated(self%k2lm))   deallocate(self%k2lm)
      if (allocated(self%kdeg))   deallocate(self%kdeg)
      if (allocated(self%kbeg))   deallocate(self%kbeg)
      if (allocated(self%mnorm))  deallocate(self%mnorm)
      if (allocated(self%Are0))   deallocate(self%Are0)
      if (allocated(self%Aim0))   deallocate(self%Aim0)
      if (allocated(self%Bre0))   deallocate(self%Bre0)
      if (allocated(self%Bim0))   deallocate(self%Bim0)
      if (allocated(self%Cre0))   deallocate(self%Cre0)
      if (allocated(self%Cim0))   deallocate(self%Cim0)
      if (allocated(self%edUn_re)) deallocate(self%edUn_re)
      if (allocated(self%edUn_im)) deallocate(self%edUn_im)
      if (allocated(self%edVn_re)) deallocate(self%edVn_re)
      if (allocated(self%edVn_im)) deallocate(self%edVn_im)
      if (allocated(self%dUa_prev)) deallocate(self%dUa_prev)
      if (allocated(self%sigma_n))    deallocate(self%sigma_n)
      if (allocated(self%sigma_next)) deallocate(self%sigma_next)
      if (allocated(self%MkPerDt))    deallocate(self%MkPerDt)
      if (allocated(self%Are_s))      deallocate(self%Are_s, self%Aim_s, self%Bre_s, &
                                                 self%Bim_s, self%Cre_s, self%Cim_s)
      if (allocated(self%Are_c))      deallocate(self%Are_c, self%Aim_c, self%Bre_c, &
                                                 self%Bim_c, self%Cre_c, self%Cim_c)
      if (allocated(self%sigma_n_s))  deallocate(self%sigma_n_s)
      self%lmax = 0;  self%nlm = 0;  self%nk = 0
      self%sigma_primed = .false.;  self%sigma_primed_s = .false.
   end subroutine ve_response_destroy

end module fe_response
