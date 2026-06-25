program test_earth
   !! Validate the M3-L70-V01 Earth-structure dataset for internal consistency.
   !!
   !! The strongest available check: the unperturbed gravity computed from the
   !! density profile, g(r) = G*M(<r)/r^2, must reproduce the gravity values
   !! Spada et al. (2011) Table 3 states at each interface. This ties the
   !! hard-coded densities and radii to independently published numbers. Mass and
   !! the moment-of-inertia factor are checked for physical plausibility.
   use fe_precision,       only: wp
   use fe_constants,       only: m_earth
   use fe_earth_structure, only: earth_moi, earth_total_mass, earth_gravity_at, earth_n_layers, earth_model, build_M3L70V01
   implicit none

   type(earth_model) :: em
   integer  :: i
   logical  :: ok
   real(wp) :: mass, moi_factor, g, gtol
   ! Benchmark interface radii [m] and stated gravities [m s^-2] (Spada 2011 Tab 3)
   real(wp), parameter :: r_if(5) = [6371.0_wp, 6301.0_wp, 5951.0_wp, &
                                     5701.0_wp, 3480.0_wp]*1.0e3_wp
   real(wp), parameter :: g_if(5) = [9.815_wp, 9.854_wp, 9.978_wp, &
                                     10.024_wp, 10.457_wp]

   ok = .true.
   em = build_M3L70V01()
   print '(a,a,a,i0,a)', ' model: ', em%name, '  (', earth_n_layers(em), ' layers)'

   ! --- Gravity profile vs benchmark interface values -------------------------
   ! 0.02 m/s^2 tolerance absorbs the benchmark's 4-sig-fig rounding and any
   ! difference in the gravitational constant convention.
   gtol = 0.02_wp
   print '(a)', '   interface       r [km]      g model     g bench      err'
   do i = 1, 5
      g = earth_gravity_at(em, r_if(i))
      print '(i10,4f12.4)', i, r_if(i)/1.0e3_wp, g, g_if(i), abs(g - g_if(i))
      if (abs(g - g_if(i)) >= gtol) ok = .false.
   end do

   ! --- Total mass: plausible vs real Earth (simplified model, so ~few %) -----
   mass = earth_total_mass(em)
   print '(a,es12.5,a,f7.4)', ' total mass = ', mass, ' kg   (M_earth ratio = ', &
        mass/m_earth, ')'
   if (mass/m_earth < 0.95_wp .or. mass/m_earth > 1.05_wp) ok = .false.

   ! --- Moment-of-inertia factor I/(M R^2): physically ~0.30-0.40 -------------
   moi_factor = earth_moi(em)/(mass*em%r_earth**2)
   print '(a,f8.5,a)', ' MOI factor I/(M R^2) = ', moi_factor, '  (Earth ~0.3307)'
   if (moi_factor < 0.30_wp .or. moi_factor > 0.40_wp) ok = .false.

   if (ok) then
      print '(a)', ' PASS: M3-L70-V01 gravity profile, mass, and MOI consistent'
   else
      print '(a)', ' FAIL: Earth-structure consistency check exceeded tolerance'
      error stop 1
   end if
end program test_earth
