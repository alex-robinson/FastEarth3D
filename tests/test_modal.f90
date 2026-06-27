program test_modal
   !! Foundation test for fe_modal (RESP_MODAL). Verifies that the modal
   !! relaxation propagator `modal_apply_p` reproduces fe_viscoelastic's
   !! homogeneous (σ=0) `ve_step` relaxation exactly — it is the same kernel,
   !! wrapped as a matrix-free linear map for the subspace eigensolve. Checks:
   !!   (1) FE propagator == reference ve_step(σ=0), step for step, to ~ULP;
   !!   (2) the relaxation actually decays the memory (and stays finite);
   !!   (3) the backward-Euler propagator (used by the eigensolve) runs, stays
   !!       finite, and relaxes (A-stable, large dt).
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_mesh, radial_mesh_build, radial_fe_finalize
   use fe_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, &
                                 SCHEME_FE, SCHEME_BE
   use fe_modal,           only: modal_degree, modal_degree_init, &
                                 modal_degree_destroy, modal_apply_p
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp
   integer,  parameter :: j = 2

   type(earth_model) :: e
   type(radial_mesh) :: m
   type(ve_degree)   :: ref
   type(modal_degree):: md
   real(wp), allocatable :: t0A(:,:), t0B(:,:), t0C(:,:)        ! seed memory τ0
   real(wp), allocatable :: cA(:,:), cB(:,:), cC(:,:)           ! modal current
   real(wp), allocatable :: oA(:,:), oB(:,:), oC(:,:)           ! modal next
   real(wp) :: dt, t, ua, va, fa, n0, maxdiff, dlast, mdone, mbe
   integer  :: ne, istep, nstep
   logical  :: ok, finite

   ok = .true.
   dt = 10.0_wp*yr
   call mk_earth(e)
   call radial_mesh_build(m, e)
   ne = m%ne
   allocate(t0A(NLAM,ne), t0B(NLAM,ne), t0C(NLAM,ne))
   allocate(cA(NLAM,ne), cB(NLAM,ne), cC(NLAM,ne))
   allocate(oA(NLAM,ne), oB(NLAM,ne), oC(NLAM,ne))

   ! --- seed a non-trivial memory state τ0 by holding a unit load for a while ---
   call ve_init(ref, e, m, j, dt)
   do istep = 1, 50
      call ve_step(ref, 1.0_wp, t, ua, va, fa)
   end do
   t0A = ref%Am;  t0B = ref%Bm;  t0C = ref%Cm
   call ve_destroy(ref)
   n0 = norm3(t0A, t0B, t0C)

   ! --- (1) FE propagator vs reference homogeneous ve_step(σ=0) -----------------
   call ve_init(ref, e, m, j, dt)              ! FE scheme (default)
   ref%Am = t0A;  ref%Bm = t0B;  ref%Cm = t0C
   ref%Un_prev = 0.0_wp;  ref%Vn_prev = 0.0_wp;  ref%time = 0.0_wp
   call modal_degree_init(md, e, m, j, dt, SCHEME_FE)
   cA = t0A;  cB = t0B;  cC = t0C
   nstep = 100;  maxdiff = 0.0_wp;  finite = .true.
   do istep = 1, nstep
      call ve_step(ref, 0.0_wp, t, ua, va, fa)         ! reference relaxation step
      call modal_apply_p(md, cA, cB, cC, oA, oB, oC)   ! P · current
      cA = oA;  cB = oB;  cC = oC
      maxdiff = max(maxdiff, norm3(cA-ref%Am, cB-ref%Bm, cC-ref%Cm))
      if (any(cA /= cA) .or. any(cB /= cB) .or. any(cC /= cC)) finite = .false.
   end do
   mdone = norm3(cA, cB, cC)
   call ve_destroy(ref)
   call modal_degree_destroy(md)

   write(*,'(a)')         ' Modal propagator foundation (homogeneous sphere, degree 2)'
   write(*,'(a,es10.3)')  '   seed memory norm  |tau0|        = ', n0
   write(*,'(a,es10.3)')  '   FE prop vs ve_step max|diff|     = ', maxdiff
   write(*,'(a,es10.3)')  '   memory after 100 steps |tau|     = ', mdone

   if (.not. finite) then
      write(*,'(a)') '   FAIL: modal FE propagator not finite';  ok = .false.
   end if
   if (maxdiff > 1.0e-12_wp*n0) then
      write(*,'(a)') '   FAIL: FE propagator /= ve_step(sigma=0)';  ok = .false.
   end if
   if (mdone >= n0) then
      write(*,'(a)') '   FAIL: homogeneous relaxation did not decay memory';  ok = .false.
   end if

   ! --- (3) backward-Euler propagator (the eigensolve uses this) ----------------
   call modal_degree_init(md, e, m, j, 1.0e3_wp*yr, SCHEME_BE)   ! large dt, A-stable
   call modal_apply_p(md, t0A, t0B, t0C, oA, oB, oC)
   mbe = norm3(oA, oB, oC)
   finite = .not. (any(oA /= oA) .or. any(oB /= oB) .or. any(oC /= oC))
   call modal_degree_destroy(md)
   write(*,'(a,es10.3)')  '   BE prop (dt=1kyr) |P.tau0|       = ', mbe
   if (.not. finite) then
      write(*,'(a)') '   FAIL: BE propagator not finite';  ok = .false.
   end if
   if (mbe >= n0) then
      write(*,'(a)') '   FAIL: BE propagator did not relax memory';  ok = .false.
   end if

   call radial_fe_finalize()
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: modal propagator matches ve_step and relaxes (FE + BE)'
   else
      write(*,'(a)') ' FAIL: modal propagator foundation checks failed'
      error stop 1
   end if

contains

   pure real(wp) function norm3(A, B, C) result(nrm)
      real(wp), intent(in) :: A(:,:), B(:,:), C(:,:)
      nrm = max(maxval(abs(A)), maxval(abs(B)), maxval(abs(C)))
   end function norm3

   subroutine mk_earth(em)
      type(earth_model), intent(out) :: em
      em%name = "maxwell";  em%r_earth = a;  em%r_core = 0.0_wp
      allocate(em%layers(1))
      em%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_modal
