program test_timestep
   !! Adaptive-Δt controller (§3c) on the VE+SLE field model. Two checks:
   !!
   !!  (A) the field step-doubling LOCAL-ERROR ESTIMATE is order-faithful: for the
   !!      trapezoidal scheme (p=2) it scales as Δt^{p+1}=Δt^3. Measured exactly like
   !!      the 1-D step_doubling_check, but through the full SLE driver, using the
   !!      ve_response controller primitives (save_state/set_dt/stash_coarse/
   !!      coarse_fine_error). Held load → isolates the integrator order.
   !!
   !!  (B) the adaptive_stepper CONTROLS the global error: driving a fast ice ramp to
   !!      a fixed end time, the adaptive result converges to a fine fixed-Δt reference,
   !!      tightening rtol lowers the error monotonically, and it reaches the reference
   !!      accuracy in FAR fewer steps than the fixed reference (the payoff that
   !!      amortizes the ~6× trapezoidal per-step cost).
   !!
   !! Same homogeneous Maxwell sphere + fixed ocean as test_sle_couple_order.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi, rho_ice, rho_water, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: ve_response
   use fe_viscoelastic,    only: SCHEME_TRAP
   use fe_sht,             only: sht_grid
   use fe_sle,             only: sle_solver, sle_result
   use fe_timestep,        only: adaptive_stepper
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp
   integer,  parameter :: LMAX = 16, MAXIT = 30
   real(wp), parameter :: ICE_MAX = 2000.0_wp
   real(wp), parameter :: T_END = 2.0_wp*kyr

   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp), allocatable :: topo0(:,:), iceF(:,:), ice0(:,:), zero(:,:)
   logical :: ok

   ok = .true.
   call sht%init(LMAX, nlat=2*LMAX, nphi=4*LMAX)
   call mk_earth(e)
   allocate(topo0(sht%nphi,sht%nlat), iceF(sht%nphi,sht%nlat), &
            ice0(sht%nphi,sht%nlat), zero(sht%nphi,sht%nlat))
   call make_topo(topo0)
   call ice_field(ICE_MAX, iceF)         ! full grounded cap
   ice0 = 0.0_wp;  zero = 0.0_wp

   write(*,'(a)') ' Adaptive-Δt controller (§3c) -- VE+SLE field model'

   call part_A_estimate_order()
   call part_B_controller()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: field step-doubling estimate is ~3rd order and the adaptive'
      write(*,'(a)') '       controller converges to the reference with far fewer steps.'
   else
      write(*,'(a)') ' FAIL: adaptive-Δt controller checks did not all pass'
      call sht%destroy();  call radial_fe_finalize();  error stop 1
   end if
   call sht%destroy();  call radial_fe_finalize()

contains

   subroutine part_A_estimate_order()
      !! Step-doubling local-error estimate vs Δt (held load), measured AFTER one warmup
      !! step so the start-of-step load σ_n is tracked (the controller's normal regime;
      !! the very first step from rest uses the σ_{n+1} proxy for σ_n, which the SLE
      !! ocean load makes O(Δt) and so caps that one step at 2nd order). Must scale as
      !! Δt^{p+1}=Δt^3 (p=2 trapezoidal).
      integer,  parameter :: NS = 4
      real(wp), parameter :: dts(NS) = [160.0_wp, 80.0_wp, 40.0_wp, 20.0_wp]  ! yr
      real(wp) :: est(NS), p
      integer  :: i
      write(*,'(a)') ''
      write(*,'(a)') '  (A) step-doubling estimate vs Δt (held load, after a warmup step)'
      write(*,'(a)') '      Δt[yr]      est ‖τ_f−τ_c‖∞/3'
      do i = 1, NS
         est(i) = one_estimate(dts(i)*yr)
         write(*,'(a,f8.1,es18.3)') '   ', dts(i), est(i)
      end do
      p = log(est(NS-1)/est(NS))/log(2.0_wp)
      write(*,'(a,f6.2)') '      estimate scaling order (p+1, expect ~3) = ', p
      if (p < 2.6_wp .or. p > 3.4_wp) then
         write(*,'(a,f5.2)') '      FAIL: step-doubling estimate not ~Δt^3, p+1=', p; ok = .false.
      end if
   end subroutine part_A_estimate_order

   function one_estimate(dt) result(est)
      !! One field step-doubling estimate of the first step from rest under a held load,
      !! exercising the controller primitives directly.
      real(wp), intent(in) :: dt
      real(wp), parameter  :: dt_warm = 40.0_wp*yr   ! fixed → identical warmed state per Δt
      real(wp) :: est
      type(ve_response) :: resp
      type(sle_solver)  :: sle
      type(sle_result)  :: res
      real(wp), allocatable :: rsl(:,:), C(:,:), rsl_n(:,:)
      real(wp) :: err_inf, tau_inf
      allocate(rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), rsl_n(sht%nphi,sht%nlat))
      call setup(resp, sle, dt)
      rsl = 0.0_wp
      ! one warmup step so σ_n is tracked (τ≠0, sigma_primed) before measuring
      call resp%set_dt(dt_warm)
      call sle%solve(sht, resp, iceF, iceF, topo0, rsl, C, res)
      call resp%save_state();  rsl_n = rsl
      call resp%set_dt(dt)
      call sle%solve(sht, resp, iceF, iceF, topo0, rsl, C, res)     ! coarse (held load)
      call resp%stash_coarse()
      call resp%restore_state();  rsl = rsl_n
      call resp%set_dt(0.5_wp*dt)
      call sle%solve(sht, resp, iceF, iceF, topo0, rsl, C, res)     ! fine 1
      call sle%solve(sht, resp, iceF, iceF, topo0, rsl, C, res)     ! fine 2
      call resp%coarse_fine_error(err_inf, tau_inf)
      est = err_inf/3.0_wp                                          ! 2^p−1 = 3
      call resp%destroy()
      deallocate(rsl, C, rsl_n)
   end function one_estimate

   subroutine part_B_controller()
      !! Drive the ice ramp 0→full over [0,T_END] and compare adaptive runs to a fine
      !! fixed-Δt reference.
      real(wp), allocatable :: rref(:,:), rad(:,:), C(:,:)
      real(wp), parameter :: tols(3) = [1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]
      real(wp) :: errs(3), sig
      integer  :: nacc(3), nrej(3), i, nref
      allocate(rref(sht%nphi,sht%nlat), rad(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))

      ! Fine fixed-Δt reference: forced Δt = 12.5 yr (the truth).
      call run_fixed(12.5_wp*yr, rref, nref)
      sig = maxval(abs(rref))
      write(*,'(a)') ''
      write(*,'(a,f9.3,a,i4,a)') '  (B) reference (fixed Δt=12.5 yr): max|rsl|=', sig, &
                                 ' m  in ', nref, ' steps'
      write(*,'(a)') '      rtol        adaptive err [m]   rel    n_accept  n_reject'
      do i = 1, 3
         call run_adaptive(tols(i), rad, nacc(i), nrej(i))
         errs(i) = maxval(abs(rad - rref))
         write(*,'(a,es9.1,es16.3,f9.4,i9,i9)') '   ', tols(i), errs(i), errs(i)/sig, &
                                                nacc(i), nrej(i)
      end do

      ! (1) tightening rtol lowers the global error monotonically
      if (.not. (errs(1) > errs(2) .and. errs(2) > errs(3))) then
         write(*,'(a)') '      FAIL: error not monotone-decreasing in rtol'; ok = .false.
      end if
      ! (2) the tight-tol adaptive result matches the reference closely
      if (errs(3) > 0.01_wp*sig) then
         write(*,'(a)') '      FAIL: tight-rtol adaptive run not close to the reference'; ok = .false.
      end if
      ! (3) the payoff: adaptive reaches reference accuracy in far fewer steps
      if (nacc(2) >= nref/2) then
         write(*,'(a)') '      FAIL: adaptive not cheaper than the fixed reference'; ok = .false.
      end if
      deallocate(rref, rad, C)
   end subroutine part_B_controller

   subroutine run_fixed(dt_fixed, rsl, nsteps)
      real(wp), intent(in)  :: dt_fixed
      real(wp), intent(out) :: rsl(:,:)
      integer,  intent(out) :: nsteps
      type(ve_response)     :: resp
      type(sle_solver)      :: sle
      type(adaptive_stepper):: st
      real(wp), allocatable :: C(:,:)
      allocate(C(sht%nphi,sht%nlat))
      call setup(resp, sle, dt_fixed)
      st%dt_try = dt_fixed;  st%dt_min = dt_fixed;  st%dt_max = dt_fixed
      rsl = 0.0_wp
      call st%advance(sht, resp, sle, topo0, ice0, iceF, zero, 0.0_wp, T_END, rsl, C)
      nsteps = st%n_accept
      call resp%destroy();  deallocate(C)
   end subroutine run_fixed

   subroutine run_adaptive(rtol, rsl, nacc, nrej)
      real(wp), intent(in)  :: rtol
      real(wp), intent(out) :: rsl(:,:)
      integer,  intent(out) :: nacc, nrej
      type(ve_response)     :: resp
      type(sle_solver)      :: sle
      type(adaptive_stepper):: st
      real(wp), allocatable :: C(:,:)
      allocate(C(sht%nphi,sht%nlat))
      call setup(resp, sle, T_END)
      st%rtol = rtol;  st%atol = 1.0e-3_wp
      st%dt_try = 0.01_wp*T_END;  st%dt_min = 0.0_wp;  st%dt_max = T_END
      rsl = 0.0_wp
      call st%advance(sht, resp, sle, topo0, ice0, iceF, zero, 0.0_wp, T_END, rsl, C)
      nacc = st%n_accept;  nrej = st%n_reject
      call resp%destroy();  deallocate(C)
   end subroutine run_adaptive

   subroutine setup(resp, sle, dt)
      type(ve_response), intent(out) :: resp
      type(sle_solver),  intent(out) :: sle
      real(wp),          intent(in)  :: dt
      call resp%init(e, sht, dt)
      resp%scheme = SCHEME_TRAP;  resp%couple_tol = 1.0e-9_wp
      sle%fixed_ocean = .true.;  sle%subgrid = .false.;  sle%max_mem_iter = MAXIT
      sle%warm_start  = .true.
   end subroutine setup

   subroutine make_topo(topo0)
      real(wp), intent(out) :: topo0(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            topo0(i,j) = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine make_topo

   subroutine ice_field(h, ice)
      real(wp), intent(in)  :: h
      real(wp), intent(out) :: ice(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            ice(i,j) = merge(h, 0.0_wp, th < 40.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine ice_field

   subroutine mk_earth(e)
      type(earth_model), intent(out) :: e
      e%name = "maxwell";  e%r_earth = a;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_timestep
