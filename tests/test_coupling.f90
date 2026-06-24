program test_coupling
   !! CLIMBER-X coupling contract (design.md item 1): init → update → finalize.
   !! Drives the solid_earth derived type with an ice load and checks that
   !!
   !!   (0) at the reference ice (h_ice = h_ice_ref) nothing moves: rsl ≈ 0 and
   !!       z_bed ≈ z_bed_eq (the supplied reference IS the relaxed state);
   !!   (1) ocean mass is conserved every internal Maxwell sub-step;
   !!   (2) a steady grounded-ice load relaxes sensibly: the bed subsides under
   !!       the ice immediately (elastic) and keeps subsiding with decreasing
   !!       increments toward an isostatic limit (viscoelastic relaxation), and
   !!       the ocean draws down (adding land ice removes ocean water).
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_params,          only: fe_param_class
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_sht,             only: sht_grid
   use fe_coupling,        only: solid_earth
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 15
   type(sht_grid), target :: sht
   type(fe_param_class)   :: p
   type(solid_earth)      :: se
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_ref(:,:), h_ice(:,:)
   real(wp) :: dt_couple, bed_eq_ice, bed_prev, bed_now, d_first, d_last
   real(wp) :: ocean_rsl, worst_mass, net_subs
   integer  :: i, jice, jocean, step
   logical  :: ok, monotone

   ok = .true.
   call sht%init(LMAX, nlat=2*LMAX, nphi=4*LMAX)

   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_ref(sht%nphi,sht%nlat), &
            h_ice(sht%nphi,sht%nlat))
   call make_reference(z_bed_eq, h_ice_ref)

   ! representative latitude rows: under the ice cap (colat ~ 15°) and open ocean
   jice   = nearest_row(15.0_wp)
   jocean = nearest_row(70.0_wp)
   bed_eq_ice = z_bed_eq(1,jice)

   dt_couple    = 2.0_wp*kyr                        ! interval per se%update call
   p%dt_couple  = dt_couple                         ! default cadence (M3-L70-V01, trapezoidal)
   call se%init(p, sht, z_bed_eq, h_ice_ref)

   write(*,'(a,i0,a)') ' coupling: lmax=', LMAX, &
                       '  (dt_couple=2 kyr, adaptive trapezoidal)'

   ! --- (0) reference consistency: no ice change -> no motion -----------------
   call se%update(h_ice_ref, dt_couple)
   write(*,'(a,es10.2,a,es10.2)') '   ref: max|rsl|=', maxval(abs(se%rsl)), &
                                  '   max|z_bed-z_bed_eq|=', maxval(abs(se%z_bed - z_bed_eq))
   if (maxval(abs(se%rsl)) > 1.0e-9_wp .or. maxval(abs(se%z_bed - z_bed_eq)) > 1.0e-9_wp) then
      write(*,'(a)') '   FAIL: reference state is not stationary'; ok = .false.
   end if

   ! --- (1)+(2) hold a 2 km grounded ice cap and watch it relax ---------------
   call make_load(h_ice)                           ! cap on colat<30 (on land)
   write(*,'(a)') '       step   z_bed under ice [m]   d(subs) [m]   worst mass resid'
   bed_prev = bed_eq_ice;  worst_mass = 0.0_wp;  monotone = .true.
   do step = 1, NSTEP
      call se%update(h_ice, dt_couple)
      bed_now = se%z_bed(1,jice)
      worst_mass = max(worst_mass, se%worst_mass_resid)
      if (step == 1)     d_first = bed_prev - bed_now
      if (step == NSTEP) d_last  = bed_prev - bed_now
      if (bed_now > bed_prev + 1.0e-6_wp) monotone = .false.   ! must not rebound
      if (mod(step,3) == 1 .or. step == NSTEP) &
         write(*,'(i9,f18.4,f15.4,es18.2)') step, bed_now, bed_prev-bed_now, &
                                            se%worst_mass_resid
      bed_prev = bed_now
   end do

   ocean_rsl = se%rsl(1,jocean)
   net_subs  = bed_eq_ice - se%z_bed(1,jice)

   write(*,'(a)') ''
   write(*,'(a,f12.4,a)')  '      net bed subsidence under ice =', net_subs, ' m'
   write(*,'(a,f12.4,a,f12.4,a)') '      subsidence increment: first =', d_first, &
                                  '   last =', d_last, ' m'
   write(*,'(a,f12.4,a)')  '      ocean rsl (drawdown, should be < 0) =', ocean_rsl, ' m'
   write(*,'(a,es11.2)')   '      worst ocean-mass residual over the run =', worst_mass

   ! --- assertions ------------------------------------------------------------
   if (worst_mass > 1.0e-10_wp) then
      write(*,'(a)') '      FAIL: ocean mass not conserved during stepping'; ok = .false.
   end if
   if (.not. monotone) then
      write(*,'(a)') '      FAIL: bed rebounds under a steady load'; ok = .false.
   end if
   if (net_subs < 10.0_wp) then
      write(*,'(a)') '      FAIL: no appreciable subsidence under the ice'; ok = .false.
   end if
   if (d_last >= d_first) then
      write(*,'(a)') '      FAIL: subsidence not slowing (no relaxation toward a limit)'
      ok = .false.
   end if
   if (ocean_rsl >= 0.0_wp) then
      write(*,'(a)') '      FAIL: building land ice should draw the ocean down'; ok = .false.
   end if

   call se%finalize()
   call sht%destroy()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: coupling init/update/finalize conserves ocean mass,'
      write(*,'(a)') '       holds the reference state, and relaxes a steady load'
   else
      write(*,'(a)') ' FAIL: coupling validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine make_reference(z_bed_eq, h_ice_ref)
      !! Reference (equilibrium) topography: a polar land cap (colat<50°, +500 m)
      !! over deep ocean (−4000 m); no reference ice.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_ref(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
         end do
      end do
      h_ice_ref = 0.0_wp
   end subroutine make_reference

   subroutine make_load(h_ice)
      !! A 2 km grounded ice cap on the land interior (colat<30°), held fixed.
      real(wp), intent(out) :: h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            h_ice(i,j) = merge(2000.0_wp, 0.0_wp, thd < 30.0_wp)
         end do
      end do
   end subroutine make_load

   integer function nearest_row(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*180.0_wp/pi - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest_row

end program test_coupling
