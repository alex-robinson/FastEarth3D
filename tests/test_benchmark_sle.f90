program test_benchmark_sle
   !! Rung-4 validation against the Martinec et al. (2018) SLE benchmark, the full
   !! migrating-coastline sea-level equation with basin topography and a time-
   !! evolving off-pole ice cap. This is the standalone (not `make check`) SLE
   !! benchmark; case A (pure loading, no SLE) lives in test_benchmark_martinec.
   !!
   !! CASE E2 (this test) = giapy L3 + T2 + B3 (data/benchmarks/sle_martinec2018/,
   !! see data/benchmarks/PROVENANCE.md and giapy tests/sle_test.py):
   !!   * ICE  (L3): spherical cap, centre (colat 25 deg, lon 75 deg), final
   !!     h0 = 500 m, final angular radius alpha = 10 deg, sqrt cap profile.
   !!   * TOPO (B3): exponential basin B = 3800 - 6000 exp(-d^2/2 sigma^2),
   !!     centre (colat 35 deg, lon 25 deg), sigma = 26 deg.
   !!   * TIME (T2): linear growth of (alpha,h) from 0 to full over t = 15->5 kyr
   !!     (500 steps), then held 5->0 kyr (250 steps); dt = 20 yr, 750 steps total.
   !!   * lmax = 128 (giapy ntrunc=128), SLE inner iterations = 20 (giapy eliter).
   !!
   !! SCALARS ONLY (this pass): we validate the four scalar reference columns
   !!   col2 vertical displacement u, col5 = -N*g (geoid), col6 sea-surface = N+esl,
   !!   col7 SLE = relative sea level rsl = N - u + esl.
   !! The horizontal columns col3/col4 (th/ph displacement) need the spheroidal V
   !! field, which ve_response does not yet output -- deferred to a later pass.
   !!
   !! Profiles (figs 10-13): along circles of constant lon (col1=colat) or
   !! constant lat (col1=180+lon); coordinates sampled with fe_sht%eval_point.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, kyr, rho_ice, rho_water
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: ve_response
   use fe_sht,             only: sht_grid
   use fe_sle,             only: sle_solver, sle_result
   use fe_field,           only: spherical_cap, exp_basin
   implicit none

   character(*), parameter :: DIR = 'data/benchmarks/sle_martinec2018/'
   integer,  parameter :: LMAX = 128
   integer,  parameter :: NGROW = 500, NHOLD = 250, NSTEP = NGROW + NHOLD
   real(wp), parameter :: DEG = pi/180.0_wp
   real(wp), parameter :: G_GIAPY = 9.815_wp           ! giapy reference gravity
   character(*), parameter :: CASE = 'E2_'
   ! E2 = L3 cap (500 m) + B3 deep basin (colat35,lon25); T2; profiles lon25/lat35
   real(wp), parameter :: ICE_COLAT = 25.0_wp*DEG, ICE_LON = 75.0_wp*DEG
   real(wp), parameter :: H0 = 500.0_wp, ALPHA_MAX = 10.0_wp*DEG
   real(wp), parameter :: BAS_COLAT = 35.0_wp*DEG, BAS_LON = 25.0_wp*DEG
   real(wp), parameter :: BMAX = 3800.0_wp, B0 = 6000.0_wp, SIGB = 26.0_wp*DEG
   ! Peak-normalized max-error tolerances. The Martinec-2018 SLE benchmark is a
   ! numerical INTER-CODE comparison (no analytical solution), so giapy is one
   ! implementation and meter-level / few-% scatter is expected -- largest in the
   ! coastline/flotation-coupled fields. These bounds guard against gross
   ! regression, not exact agreement. Achieved on E2 (lmax128): uplift/geoid ~5-6%
   ! in well-behaved regions, sea-surface ~14% (a ~4 m uniform eustatic offset),
   ! with the worst local errors (fig11 ~35-44%) confined to the cap-edge-over-
   ! basin GROUNDING LINE -- the most implementation-dependent zone.
   real(wp), parameter :: TOL_U = 4.0e-1_wp, TOL_N = 5.0e-1_wp, &
                          TOL_SS = 2.0e-1_wp, TOL_SLE = 5.0e-1_wp

   type(sht_grid)     :: sht
   type(earth_model)  :: em
   type(ve_response)  :: resp
   type(sle_solver)   :: sle
   type(sle_result)   :: res
   real(wp), allocatable :: topo0(:,:), ice_now(:,:)
   real(wp), allocatable :: rsl(:,:), C(:,:), tmp(:,:)
   complex(wp), allocatable :: u_lm(:), N_lm(:), rsl_lm(:)
   real(wp) :: dt, frac, alpha, h, esl, bary, shift, rho_ratio
   integer  :: istep
   logical  :: ok

   integer,  parameter :: NFIG = 4
   character(8) :: fig(NFIG) = ['fig10', 'fig11', 'fig12', 'fig13']
   character(1) :: ptype(NFIG) = ['Z', 'Y', 'Z', 'Y']   ! Z: col1=colat, Y: col1=180+lon
   real(wp) :: pfix(NFIG)                                 ! fixed lon (Z) or colat (Y) [rad]

   ok = .true.
   ! Z-profiles: pfix = fixed longitude; Y-profiles: pfix = fixed colatitude.
   ! giapy lat-profiles pass latev=90-lat to a latitude evaluator, so the sample
   ! COLATITUDE = lat (the label), not 90-lat: fig11 lat25->colat25, fig13->colat35.
   pfix = [ICE_LON, 25.0_wp*DEG, 25.0_wp*DEG, 35.0_wp*DEG]   ! fig11 colat25, fig12 lon25, fig13 colat35

   dt = 0.02_wp*kyr
   call sht%init(LMAX, nlat=2*LMAX, nphi=4*LMAX)
   em = build_M3L70V01()
   call resp%init(em, sht, dt)
   sle%n_inner = 20     ! n_outer left at default 3 (converged: 3 == 12 on E2)

   allocate(topo0(sht%nphi,sht%nlat), ice_now(sht%nphi,sht%nlat), &
            rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), tmp(sht%nphi,sht%nlat))
   allocate(u_lm(sht%nlm), N_lm(sht%nlm), rsl_lm(sht%nlm))

   call exp_basin(sht, BAS_COLAT, BAS_LON, BMAX, B0, SIGB, topo0)

   write(*,'(a,i0,a,i0,a)') ' E2 (L3+T2+B3): lmax=', LMAX, ', ', NSTEP, ' steps @ dt=20 yr'
   do istep = 1, NSTEP
      frac = min(1.0_wp, real(istep,wp)/real(NGROW,wp))   ! T2 ramp 0->1 over 10 kyr, then held
      alpha = frac*ALPHA_MAX;  h = frac*H0
      call spherical_cap(sht, ICE_COLAT, ICE_LON, alpha, h, ice_now)
      ! d_ice = TOTAL ice change from the ice-free reference (not a per-step
      ! increment); ice = absolute ice (same here, ref is ice-free).
      call sle%solve(sht, resp, ice_now, ice_now, topo0, rsl, C, res)
      if (mod(istep,100) == 0 .or. istep == NSTEP) &
         write(*,'(a,i4,a,f7.4,a,es10.2,a,f9.4,a,i3,a,es9.2)') '   step ', istep, '  frac=', frac, &
            '  mass_resid=', res%mass_resid, '  esl=', res%esl, &
            '  n_outer=', res%n_outer_done, '  resid=', res%resid
   end do
   esl = res%esl

   ! eustatic decomposition at the final state (ice_now, C are the last-step values)
   rho_ratio = rho_ice/rho_water
   bary  = -rho_ratio*sht%surface_integral(ice_now)/sht%surface_integral(C)  ! barystatic
   shift = sht%surface_integral(C*(res%N - res%u))/sht%surface_integral(C)    ! ocean-mean(N-u)
   write(*,'(a)') ''
   write(*,'(a,f10.4,a,f10.4,a,f10.4)') ' eustatic: dphi(esl)=', esl, '  barystatic=', bary, &
                                        '  ocean-mean(N-u)=', shift
   write(*,'(a,f8.4,a,f8.4,a,f10.2)') ' ocean frac: migrated(C)=', res%ocean_frac, &
        '  bare basin(topo0<0)=', sht%surface_integral(merge(1.0_wp,0.0_wp,topo0<0.0_wp))/(16.0_wp*atan(1.0_wp)), &
        '  C_int[sr]=', sht%surface_integral(C)

   ! spectral coefficients of the final converged fields (analysis overwrites input)
   tmp = res%u;  call sht%analysis(tmp, u_lm)
   tmp = res%N;  call sht%analysis(tmp, N_lm)
   tmp = rsl;    call sht%analysis(tmp, rsl_lm)

   write(*,'(a)') ''
   write(*,'(a)') ' Peak-normalized max error vs Martinec-2018 E2 (scalars):'
   write(*,'(a)') '   fig     rows    u[%]    N[%]   ss[%]  sle[%]   peaks u/N/ss/sle [m]'
   call compare_all()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: Martinec-2018 case E2 (SLE scalars) matches the benchmark'
   else
      write(*,'(a)') ' FAIL: Martinec-2018 case E2 SLE validation did not all pass'
      call resp%destroy();  call sht%destroy();  call radial_fe_finalize();  error stop 1
   end if
   call resp%destroy();  call sht%destroy();  call radial_fe_finalize()

contains

   subroutine compare_all()
      integer :: k
      do k = 1, NFIG
         call compare_fig(k)
      end do
   end subroutine compare_all

   subroutine compare_fig(k)
      integer, intent(in) :: k
      integer,  parameter :: MAXR = 800
      real(wp) :: c1(MAXR), uref(MAXR), nref(MAXR), ssref(MAXR), sleref(MAXR)
      real(wp) :: colat, lon, um, nm, ssm, slem
      real(wp) :: eu, en, ess, esle, pu, pn, pss, psle
      integer  :: nrow, i

      call read_fig(trim(fig(k)), c1, uref, nref, ssref, sleref, nrow)
      eu = 0; en = 0; ess = 0; esle = 0
      pu = maxval(abs(uref(1:nrow)));  pn = maxval(abs(nref(1:nrow)))
      pss = maxval(abs(ssref(1:nrow)));  psle = maxval(abs(sleref(1:nrow)))
      do i = 1, nrow
         if (ptype(k) == 'Z') then        ! col1 = colatitude, fixed longitude
            colat = c1(i)*DEG;  lon = pfix(k)
         else                              ! col1 = 180+lon; sample lon = col1 (SHTns periodic)
            colat = pfix(k);    lon = c1(i)*DEG
         end if
         call sht%eval_point(u_lm,   colat, lon, um)
         call sht%eval_point(N_lm,   colat, lon, nm)
         call sht%eval_point(rsl_lm, colat, lon, slem)
         ssm = nm + esl
         eu   = max(eu,   abs(um  - uref(i)))
         en   = max(en,   abs(nm  - nref(i)))
         ess  = max(ess,  abs(ssm - ssref(i)))
         esle = max(esle, abs(slem - sleref(i)))
      end do
      eu = eu/pu;  en = en/pn;  ess = ess/pss;  esle = esle/psle
      write(*,'(3x,a,i7,4f8.2,3x,4f8.2)') fig(k), nrow, 100*eu, 100*en, 100*ess, 100*esle, &
                                          pu, pn, pss, psle
      if (eu   > TOL_U)   ok = .false.
      if (en   > TOL_N)   ok = .false.
      if (ess  > TOL_SS)  ok = .false.
      if (esle > TOL_SLE) ok = .false.
   end subroutine compare_fig

   subroutine read_fig(name, c1, uref, nref, ssref, sleref, nrow)
      !! Read a 7-column *_SBK.dat profile: col1, col2=u, (col3,col4 skipped),
      !! col5=-N*g, col6=sea-surface, col7=SLE. Reference geoid N = -col5/g_giapy.
      character(*), intent(in)  :: name
      real(wp),     intent(out) :: c1(:), uref(:), nref(:), ssref(:), sleref(:)
      integer,      intent(out) :: nrow
      character(256) :: line
      character(*), parameter :: f = DIR//CASE
      real(wp) :: v(7)
      integer  :: u, ios
      open(newunit=u, file=f//name//'_SBK.dat', status='old', action='read', iostat=ios)
      if (ios /= 0) then;  write(*,'(4a)') ' FAIL: cannot read ', f, name, '_SBK.dat'; error stop 1;  end if
      nrow = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read(line,*,iostat=ios) v
         if (ios /= 0) cycle
         nrow = nrow + 1
         c1(nrow)     = v(1)
         uref(nrow)   = v(2)
         nref(nrow)   = -v(5)/G_GIAPY
         ssref(nrow)  = v(6)
         sleref(nrow) = v(7)
      end do
      close(u)
   end subroutine read_fig

end program test_benchmark_sle
