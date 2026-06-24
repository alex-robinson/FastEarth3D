program test_tensor_sh
   !! Rung 6 — axisymmetric tensor-SH dyadic transforms (fe_tensor_sh). Two checks,
   !! no external data:
   !!
   !!   (1) round trip — analysis∘synth = identity on random coefficients (the
   !!       per-channel orthogonal bases must invert exactly under Gauss quadrature);
   !!   (2) physical correctness — the tensor double-dot ∫τ:ε dΩ computed from the
   !!       reconstructed dyadic grid fields (weights [1,½,½,½] on [rr,rθ,tr,df])
   !!       must equal the spectral form Σ_λ norm_λ τ^λ ε^λ with the Martinec (2000)
   !!       B13 norms [1, j(j+1)/2, 2j²(j+1)², 2(j−1)j(j+1)(j+2)]. Round-trip alone
   !!       is convention-blind; this ties the bases to the validated spectral norms.
   use fe_precision, only: wp
   use fe_constants, only: pi
   use fe_sht,       only: sht_grid
   use fe_tensor_sh, only: tensor_sh, TLAM, DY_RR, DY_RT, DY_TR, DY_DF
   implicit none

   integer, parameter :: LMAX = 12
   type(sht_grid)  :: sht
   type(tensor_sh) :: tsh
   real(wp), allocatable :: ctau(:,:), ceps(:,:), c2(:,:), dtau(:,:), deps(:,:)
   real(wp) :: rt_err, dd_phys, dd_spec, jj, nrm(TLAM)
   integer  :: lam, j, seed
   logical  :: ok

   ok = .true.
   call sht%init(LMAX, nlat=3*LMAX, nphi=2, mmax=0)
   call tsh%init(sht)
   allocate(ctau(TLAM,LMAX), ceps(TLAM,LMAX), c2(TLAM,LMAX))
   allocate(dtau(sht%nlat,4), deps(sht%nlat,4))

   ! Deterministic pseudo-random coefficients; zero the degree-1 Z⁶ (no such harmonic).
   seed = 1
   do j = 1, LMAX
      do lam = 1, TLAM
         ctau(lam,j) = frand(seed)
         ceps(lam,j) = frand(seed)
      end do
   end do
   ctau(4,1) = 0.0_wp;  ceps(4,1) = 0.0_wp

   ! (1) round trip
   call tsh%synth(ctau, dtau)
   call tsh%analysis(dtau, c2)
   rt_err = maxval(abs(c2 - ctau))
   write(*,'(a,es11.2)') ' (1) round-trip max|analysis(synth(c)) - c| =', rt_err
   if (rt_err > 1.0e-10_wp) then
      write(*,'(a)') '     FAIL: dyadic transform does not round-trip'
      ok = .false.
   end if

   ! (2) physical double-dot == spectral double-dot
   call tsh%synth(ceps, deps)
   dd_phys = 2.0_wp*pi*sum( tsh%w * ( dtau(:,DY_RR)*deps(:,DY_RR) &
                + 0.5_wp*dtau(:,DY_RT)*deps(:,DY_RT) &
                + 0.5_wp*dtau(:,DY_TR)*deps(:,DY_TR) &
                + 0.5_wp*dtau(:,DY_DF)*deps(:,DY_DF) ) )
   dd_spec = 0.0_wp
   do j = 1, LMAX
      jj = real(j,wp)*real(j+1,wp)
      nrm = [ 1.0_wp, 0.5_wp*jj, 2.0_wp*jj*jj, &
              2.0_wp*real(j-1,wp)*real(j,wp)*real(j+1,wp)*real(j+2,wp) ]
      do lam = 1, TLAM
         dd_spec = dd_spec + nrm(lam)*ctau(lam,j)*ceps(lam,j)
      end do
   end do
   write(*,'(a,es14.6)') ' (2) physical  ∫τ:ε dΩ =', dd_phys
   write(*,'(a,es14.6)') '     spectral  Σ norm_λ τ^λ ε^λ =', dd_spec
   write(*,'(a,es11.2)') '     relative difference =', abs(dd_phys-dd_spec)/abs(dd_spec)
   if (abs(dd_phys-dd_spec)/abs(dd_spec) > 1.0e-10_wp) then
      write(*,'(a)') '     FAIL: dyadic reconstruction inconsistent with the B13 norms'
      ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: axisymmetric tensor-SH dyadic transforms validated'
   else
      write(*,'(a)') ' FAIL: tensor-SH dyadic transforms did not all pass'
      call tsh%destroy();  call sht%destroy();  error stop 1
   end if
   call tsh%destroy();  call sht%destroy()

contains

   real(wp) function frand(s) result(r)
      !! Tiny LCG in [-1,1), deterministic, no Math.random/state dependence.
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      r = 2.0_wp*real(s,wp)/2147483647.0_wp - 1.0_wp
   end function frand

end program test_tensor_sh
