program test_ve_response
   !! Rung-4 viscoelastic field driver (ve_response). Two checks pin it to
   !! already-validated code with no new reference data:
   !!
   !!   (1) elastic limit — with zero memory (the first step) the affine drift is
   !!       zero, so the per-degree gains gu(l), gn(l) must equal the standalone
   !!       elastic_response gains to machine precision;
   !!   (2) 1-D stepper agreement — driving a single held degree-2 coefficient
   !!       through begin/apply/commit must reproduce the validated 1-D
   !!       ve_degree relaxation history U(a,t), F(a,t), confirming the
   !!       per-(l,m) memory path is the same Maxwell kernel.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_mesh, radial_fe_finalize
   use fe_viscoelastic,    only: ve_degree
   use fe_response,        only: elastic_response, ve_response
   use fe_sht,             only: sht_grid
   implicit none

   integer, parameter :: LMAX = 8
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   real(wp) :: dt
   logical  :: ok

   ok = .true.
   dt = 0.02_wp*kyr                 ! 20 yr explicit step (VEGA)
   call sht%init(LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()

   write(*,'(a)') ' (1) elastic limit: ve_response gains == elastic_response'
   call elastic_limit(ok)

   write(*,'(a)') ''
   write(*,'(a)') ' (2) held degree-2 load reproduces the 1-D ve_degree stepper'
   call stepper_agreement(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: ve_response matches the elastic gains and the 1-D'
      write(*,'(a)') '       Maxwell stepper through the per-(l,m) memory path'
   else
      write(*,'(a)') ' FAIL: ve_response validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call sht%destroy()
   call radial_fe_finalize()

contains

   subroutine elastic_limit(ok)
      logical, intent(inout) :: ok
      type(ve_response)      :: ve
      type(elastic_response) :: el
      integer  :: l
      real(wp) :: du, dn, dmax
      call ve%init(e, sht, dt)
      call el%init(e, lmax=LMAX)
      ! The field driver carries degrees l>=2 (degree 0/1 handled separately:
      ! degree-1 deformation is deferred to the CM-frame treatment), so compare
      ! the gains there.
      dmax = 0.0_wp
      do l = 2, LMAX
         du = abs(ve%gu(l) - el%ugain(l))
         dn = abs(ve%gn(l) - el%ngain(l))
         dmax = max(dmax, du, dn)
      end do
      write(*,'(a,es11.2)') '      max |gain difference| (l>=2) =', dmax
      if (dmax > 1.0e-12_wp) then
         write(*,'(a)') '      FAIL: ve gains differ from the elastic response'
         ok = .false.
      end if
      if (ve%gu(1) /= 0.0_wp .or. ve%gn(1) /= 0.0_wp) then
         write(*,'(a)') '      FAIL: degree-1 deformation should be zeroed'
         ok = .false.
      end if
      call ve%destroy();  call el%destroy()
   end subroutine elastic_limit

   subroutine stepper_agreement(ok)
      logical, intent(inout) :: ok
      type(ve_response)  :: ve
      type(ve_degree)    :: vd
      type(radial_mesh)  :: mesh
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
      integer, parameter :: NSTEP = 25
      integer  :: lm, i
      real(wp) :: sigma, t1, ua1, va1, fa1, ua2, fa2, g
      real(wp) :: emax_u, emax_f, ref_u

      sigma = 1.0_wp
      ! reference: validated 1-D stepper at degree 2
      call mesh%build(e)
      call vd%init(e, mesh, j=2, dt=dt)

      ! field driver, single held (l=2,m=0) coefficient
      call ve%init(e, sht, dt)
      g  = ve%g
      lm = sht%lmidx(2, 0)
      allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
      slm = (0.0_wp, 0.0_wp);  slm(lm) = cmplx(sigma, 0.0_wp, wp)

      emax_u = 0.0_wp;  emax_f = 0.0_wp;  ref_u = 0.0_wp
      write(*,'(a)') '       step   U_a(1D)      U_a(field)     F_a(1D)      F_a(field)'
      do i = 1, NSTEP
         call vd%step(sigma, t1, ua1, va1, fa1)          ! 1-D: state then advance
         call ve%begin_step(sht)
         call ve%apply(sht, slm, ulm, nlm)               ! field: state at this t
         ua2 = real(ulm(lm), wp)
         fa2 = -real(nlm(lm), wp)*g                       ! N = -F/g  ->  F = -N g
         call ve%commit_step(sht, slm)                    ! advance memory
         emax_u = max(emax_u, abs(ua1 - ua2))
         emax_f = max(emax_f, abs(fa1 - fa2))
         ref_u  = max(ref_u, abs(ua1))
         if (mod(i,5) == 1 .or. i == NSTEP) &
            write(*,'(i9,4es14.5)') i, ua1, ua2, fa1, fa2
         ! imaginary channel must stay exactly zero
         if (abs(aimag(ulm(lm))) > 1.0e-14_wp) then
            write(*,'(a)') '      FAIL: imaginary channel leaked'; ok = .false.
         end if
      end do
      write(*,'(a,es11.2,a,es11.2)') '      max|ΔU_a| =', emax_u, &
                                     '   max|ΔF_a| =', emax_f
      write(*,'(a,es11.2)')          '      relative U_a error =', emax_u/ref_u
      if (emax_u/ref_u > 1.0e-6_wp .or. emax_f/abs(fa1) > 1.0e-6_wp) then
         write(*,'(a)') '      FAIL: field driver disagrees with the 1-D stepper'
         ok = .false.
      end if
      call vd%destroy();  call ve%destroy()
   end subroutine stepper_agreement

end program test_ve_response
