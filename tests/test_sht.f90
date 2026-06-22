program test_sht
   !! Round-trip test for the SHTns transform wrapper (fe_sht): synthesize a
   !! known band-limited spectrum to the spatial grid, analyze it back, and check
   !! the coefficients are recovered to machine precision. This proves the build
   !! system, the SHTns linkage, and the spectral<->spatial kernel before any
   !! physics is added.
   use, intrinsic :: iso_c_binding, only: c_double
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid
   implicit none

   type(sht_grid) :: g
   integer  :: lmax, l, m, lm
   real(wp),    allocatable :: sh(:,:)
   complex(wp), allocatable :: slm(:), slm0(:)
   real(wp) :: err, tol

   lmax = 24
   call g%init(lmax, nlat=48, nphi=96)
   print '(a,i0,a,i0,a,i0,a,i0)', ' SHTns grid: lmax=', g%lmax, &
        '  nlat=', g%nlat, '  nphi=', g%nphi, '  nlm=', g%nlm

   allocate(sh(g%nphi, g%nlat), slm(g%nlm), slm0(g%nlm))

   ! Deterministic band-limited spectrum. m=0 coefficients must be real for a
   ! real-valued spatial field; m>0 may be complex.
   slm0 = (0.0_wp, 0.0_wp)
   do l = 0, lmax
      do m = 0, l
         lm = g%lmidx(l, m)
         if (m == 0) then
            slm0(lm) = cmplx(1.0_wp/real(l+1, wp), 0.0_wp, wp)
         else
            slm0(lm) = cmplx(1.0_wp/real(l+1, wp), &
                             0.1_wp*real(m, wp)/real(l+1, wp), wp)
         end if
      end do
   end do

   slm = slm0
   call g%synthesis(slm, sh)    ! spectral -> spatial
   call g%analysis(sh, slm)     ! spatial  -> spectral (round trip)

   err = maxval(abs(slm - slm0))
   tol = 1.0e-11_wp
   print '(a,es12.4,a,es12.4)', ' round-trip max error = ', err, '   tol = ', tol

   call g%destroy()

   if (err < tol) then
      print '(a)', ' PASS: SHTns synthesis/analysis round trip'
   else
      print '(a)', ' FAIL: round-trip error exceeds tolerance'
      error stop 1
   end if
end program test_sht
