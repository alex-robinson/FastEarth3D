program test_love
   !! Rung-2 validation: solve the per-degree saddle-point system with LIS and
   !! check the surface loading Love numbers. Two analytic limits of a
   !! homogeneous incompressible self-gravitating sphere pin the solver with no
   !! external benchmark table:
   !!
   !!   fluid (μ→0):  h_j → −(2j+1)/3  and  k_j → −1  (full compensation)
   !!   rigid (μ→∞):  h_j, l_j, k_j → 0 (the sphere does not deform)
   !!
   !! and the elastic benchmark model M3-L70-V01 must converge to physical
   !! values (printed for comparison against Spada et al. 2011 Test 2/1).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, earth_layer, build_M3L70V01, &
                                 RHEOL_ELASTIC, RHEOL_FLUID
   use fe_radial_fe,       only: radial_mesh, radial_operator, loading_love, &
                                 radial_fe_finalize
   implicit none

   real(wp), parameter :: km = 1.0e3_wp
   integer  :: j, it
   real(wp) :: h, l, k, ua, va, fa, rsd, hf
   logical  :: ok
   type(radial_operator) :: op

   ok = .true.

   ! --- 1. fluid limit: h → −(2j+1)/3, k → −1 ---------------------------------
   write(*,'(a)') ' (1) homogeneous near-fluid sphere -> fluid limits'
   write(*,'(a)') '      j      h        h=-(2j+1)/3      k         resid'
   do j = 2, 6
      call solve_homogeneous(1.0_wp, j, ua, va, fa, it, rsd)   ! mu ~ 0
      call love_for(j, ua, va, fa, h, l, k)
      hf = -real(2*j+1, wp)/3.0_wp
      write(*,'(i7,2f12.5,f12.5,es11.2)') j, h, hf, k, rsd
      if (abs(h - hf) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: h off the fluid limit'; ok = .false.
      end if
      if (abs(k + 1.0_wp) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: k off the fluid limit (-1)'; ok = .false.
      end if
      if (rsd > 1.0e-8_wp) then
         write(*,'(a)') '      FAIL: LIS did not converge'; ok = .false.
      end if
   end do

   ! --- 2. rigid limit: h, l, k → 0 -------------------------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) near-rigid sphere (mu=1e18) -> h,l,k -> 0'
   write(*,'(a)') '      j       h           l           k'
   do j = 2, 6
      call solve_homogeneous(1.0e18_wp, j, ua, va, fa, it, rsd)
      call love_for(j, ua, va, fa, h, l, k)
      write(*,'(i7,3es13.4)') j, h, l, k
      if (abs(h) > 1.0e-4_wp .or. abs(l) > 1.0e-4_wp .or. abs(k) > 1.0e-4_wp) then
         write(*,'(a)') '      FAIL: not approaching the rigid limit'; ok = .false.
      end if
   end do

   ! --- 3. elastic M3-L70-V01 (vs Spada 2011 Test 2/1) ------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (3) elastic M3-L70-V01 loading Love numbers'
   write(*,'(a)') '      j       h           l           k        iters  resid'
   call elastic_M3()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: LIS saddle-point solve reproduces the analytic'
      write(*,'(a)') '       fluid and rigid Love-number limits'
   else
      write(*,'(a)') ' FAIL: Love-number validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine solve_homogeneous(mu, j, ua, va, fa, it, rsd)
      !! Single-layer incompressible sphere of rigidity mu, solved at degree j.
      real(wp), intent(in)  :: mu
      integer,  intent(in)  :: j
      real(wp), intent(out) :: ua, va, fa, rsd
      integer,  intent(out) :: it
      type(earth_model) :: e
      type(radial_mesh) :: m
      e%name = "homog";  e%r_earth = 6371.0_wp*km;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, 6371.0_wp*km, 5511.0_wp, mu, &
                                huge(1.0_wp), RHEOL_ELASTIC)
      call m%build(e)
      call op%assemble(e, m, j)
      call op%solve(1.0_wp, ua, va, fa, iters=it, resid=rsd)
   end subroutine solve_homogeneous

   subroutine love_for(j, ua, va, fa, h, l, k)
      integer,  intent(in)  :: j
      real(wp), intent(in)  :: ua, va, fa
      real(wp), intent(out) :: h, l, k
      type(earth_model) :: e
      e%r_earth = 6371.0_wp*km;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, 6371.0_wp*km, 5511.0_wp, 1.0_wp, &
                                huge(1.0_wp), RHEOL_ELASTIC)
      call loading_love(e, j, 1.0_wp, ua, va, fa, h, l, k)
   end subroutine love_for

   subroutine elastic_M3()
      type(earth_model) :: e
      type(radial_mesh) :: m
      integer  :: jj, iter
      real(wp) :: uu, vv, ff, hh, ll, kk, rr
      e = build_M3L70V01()
      call m%build(e)
      do jj = 2, 8
         call op%assemble(e, m, jj)
         call op%solve(1.0_wp, uu, vv, ff, iters=iter, resid=rr)
         call loading_love(e, jj, 1.0_wp, uu, vv, ff, hh, ll, kk)
         write(*,'(i7,3es13.5,i6,es10.2)') jj, hh, ll, kk, iter, rr
         if (rr > 1.0e-6_wp) then
            write(*,'(a)') '      FAIL: M3 solve did not converge'; ok = .false.
         end if
         if (jj == 2) then
            if (kk > -0.25_wp .or. kk < -0.55_wp .or. hh > -0.3_wp) then
               write(*,'(a)') '      FAIL: M3 degree-2 Love numbers unphysical'
               ok = .false.
            end if
         end if
      end do
   end subroutine elastic_M3

end program test_love
