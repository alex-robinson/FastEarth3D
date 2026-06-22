module fe_viscoelastic
   !! Time-domain viscoelastic relaxation: the heart of the method.
   !!
   !! Incompressible Maxwell rheology integrated explicitly in time (Martinec
   !! 2000 §3, the ω=1 explicit scheme, eqs 23-25). Each step the deviatoric
   !! stress splits into an instantaneous elastic part — the SAME per-degree
   !! saddle-point operator as the elastic problem — plus a viscous *memory*
   !! stress τ^V carried from the previous step:
   !!
   !!     τ^{i+1} = τ^{E,i+1} + τ^{V,i},   τ^{V,i} = (1−M)τ^{V,i-1} − 2μ M ε^i,
   !!
   !! with M = μΔt/η (eq 17). The memory stress contributes the dissipative
   !! forcing −∫ τ^{V,i}:δε dV (eq 35) to the right-hand side; the left-hand
   !! operator never changes, so it is assembled/equilibrated once. The explicit
   !! scheme is conditionally stable: Δt ≲ 2η_min/μ ⇒ a viscosity floor.
   !!
   !! This module provides the per-degree (1-D) stepper `ve_degree`. For a
   !! radially symmetric viscosity the memory stress evolves directly on the
   !! tensor spherical-harmonic coefficients (§9, eqs 105-110) — no spatial grid;
   !! the spheroidal strain keeps the four tensor components λ ∈ {1,2,5,6}. The
   !! laterally varying (3-D) case re-uses this same path with the memory update
   !! done pointwise on the Gauss grid (the project's 3-D goal, built on top).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_operator, radial_mesh, &
                                 idx_u, idx_v, idx_f, ndof_of
   implicit none
   private

   public :: ve_degree
   ! Per-element Maxwell kernel, shared with the field driver (fe_response):
   public :: NLAM, strain_coeffs, ve_strain_constants, dissipative_rhs, &
             advance_memory

   ! Spheroidal strain keeps four tensor-harmonic components; LAM maps the local
   ! index 1..4 to Martinec's λ ∈ {1,2,5,6} (λ=3,4 are toroidal, dropped).
   integer, parameter :: NLAM = 4

   type :: ve_degree
      !! Explicit Maxwell time stepper for a single spherical-harmonic degree j.
      integer  :: j  = -1
      integer  :: nr = 0, ne = 0, ndof = 0
      real(wp) :: dt = 0.0_wp, time = 0.0_wp, Jr = 0.0_wp
      type(radial_operator) :: op            !! elastic operator (assembled once)
      real(wp), allocatable :: r(:)          !! node radii (nr)
      real(wp), allocatable :: mu(:)         !! element shear modulus (ne)
      real(wp), allocatable :: Mk(:)         !! element M = μΔt/η (ne)
      real(wp), allocatable :: norm(:)       !! Z^λ:Z^λ norms (NLAM)
      real(wp) :: sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)  !! test-dof strain coeffs
      ! memory-stress coefficients per element, per λ (eq 109: A/h+Bψ_k/r+Cψ_{k+1}/r)
      real(wp), allocatable :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      real(wp), allocatable :: Un(:), Vn(:)  !! current nodal U, V (nr)
   contains
      procedure :: init  => ve_init
      procedure :: step  => ve_step
      procedure :: destroy => ve_destroy
   end type ve_degree

contains

   subroutine ve_init(self, earth, mesh, j, dt)
      !! Assemble the elastic operator and set up the per-element Maxwell factors
      !! M = μΔt/η. Elastic layers (η→∞) get M→0 (frozen, never relax); fluid
      !! layers (μ=0) carry no memory stress.
      class(ve_degree),  intent(inout) :: self
      type(earth_model), intent(in)    :: earth
      type(radial_mesh), intent(in)    :: mesh
      integer,           intent(in)    :: j
      real(wp),          intent(in)    :: dt
      integer  :: e, lay
      real(wp) :: eta_e

      call self%destroy()
      self%j  = j;  self%dt = dt;  self%time = 0.0_wp
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%Jr = real(j, wp)*real(j+1, wp)
      call self%op%assemble(earth, mesh, j)

      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne))
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         if (eta_e > 0.0_wp) then
            self%Mk(e) = self%mu(e)*dt/eta_e        ! ~0 for elastic (η huge)
         else
            self%Mk(e) = 0.0_wp                     ! fluid (μ=0): no memory
         end if
      end do

      ! Z^λ:Z^λ orthogonality norms (eqs B13/110) and the strain coefficients of
      ! the four unit test dofs (δU_k, δU_{k+1}, δV_k, δV_{k+1}); both r-
      ! independent, so precompute once. Shared with the field driver.
      allocate(self%norm(NLAM))
      call ve_strain_constants(self%Jr, self%norm, self%sa, self%sb, self%sc)

      allocate(self%Am(NLAM,self%ne), self%Bm(NLAM,self%ne), self%Cm(NLAM,self%ne))
      self%Am = 0.0_wp;  self%Bm = 0.0_wp;  self%Cm = 0.0_wp
      allocate(self%Un(self%nr), self%Vn(self%nr))
      self%Un = 0.0_wp;  self%Vn = 0.0_wp
   end subroutine ve_init

   subroutine ve_step(self, sigma, t_now, U_a, V_a, F_a)
      !! Advance the held degree-j load of coefficient `sigma` by one Δt and
      !! return the surface response at the time BEFORE advancing (so the first
      !! call returns the elastic t=0 state). The memory stress is updated from
      !! the new strain, ready for the next step.
      class(ve_degree), intent(inout) :: self
      real(wp),         intent(in)    :: sigma
      real(wp),         intent(out)   :: t_now, U_a, V_a, F_a
      real(wp), allocatable :: f(:), x(:)
      integer :: node

      t_now = self%time
      allocate(f(self%ndof), x(self%ndof))
      f = self%op%load_rhs(sigma)            ! elastic load forcing (eq 84)
      call add_dissipative(self, f)          ! + memory forcing (eqs 94,110)
      call self%op%solve_vec(f, x)

      do node = 1, self%nr
         self%Un(node) = x(idx_u(node))
         self%Vn(node) = x(idx_v(node))
      end do
      U_a = x(idx_u(self%nr));  V_a = x(idx_v(self%nr));  F_a = x(idx_f(self%nr))

      call update_memory(self)               ! τ^{V,i} <- (1-M)τ^{V,i-1} - 2μM ε^i
      self%time = self%time + self%dt
   end subroutine ve_step

   ! --- internals -------------------------------------------------------------

   subroutine add_dissipative(self, f)
      !! Add the dissipative memory forcing −∫ τ^{V}:δε dV to this degree's RHS.
      class(ve_degree), intent(in)    :: self
      real(wp),         intent(inout) :: f(:)
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, &
                           self%norm, self%Am, self%Bm, self%Cm, f)
   end subroutine add_dissipative

   subroutine update_memory(self)
      !! Explicit Maxwell update of this degree's memory stress from the new
      !! nodal strain.
      class(ve_degree), intent(inout) :: self
      call advance_memory(self%ne, self%mu, self%Mk, self%Un, self%Vn, &
                          self%Jr, self%Am, self%Bm, self%Cm)
   end subroutine update_memory

   ! --- shared per-element Maxwell kernel (reused by the field driver) --------

   pure subroutine strain_coeffs(u1, u2, v1, v2, Jr, a, b, c)
      !! Spheroidal strain-tensor coefficients (a,b,c) for the four kept tensor
      !! components λ = 1,2,5,6 in terms of an element's nodal U,V (Martinec eq
      !! 87, W dropped). The strain is ε = a/h + b ψ_k/r + c ψ_{k+1}/r (eq 88).
      real(wp), intent(in)  :: u1, u2, v1, v2, Jr
      real(wp), intent(out) :: a(NLAM), b(NLAM), c(NLAM)
      a(1) = -u1 + u2;        b(1) = 0.0_wp;            c(1) = 0.0_wp           ! λ=1
      a(2) = -v1 + v2;        b(2) = u1 - v1;           c(2) = u2 - v2          ! λ=2
      a(3) = 0.0_wp;          b(3) = -u1/Jr + 0.5_wp*v1; c(3) = -u2/Jr + 0.5_wp*v2 ! λ=5
      a(4) = 0.0_wp;          b(4) = 0.5_wp*v1;         c(4) = 0.5_wp*v2        ! λ=6
   end subroutine strain_coeffs

   pure subroutine ve_strain_constants(Jr, norm, sa, sb, sc)
      !! Per-degree, r-independent constants: the Z^λ:Z^λ norms (eqs B13/110) and
      !! the strain coefficients of the four unit test dofs.
      real(wp), intent(in)  :: Jr
      real(wp), intent(out) :: norm(NLAM), sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)
      norm = [ 1.0_wp, 0.5_wp*Jr, 2.0_wp*Jr**2, 2.0_wp*Jr*(Jr - 2.0_wp) ]
      call strain_coeffs(1.0_wp,0.0_wp,0.0_wp,0.0_wp, Jr, sa(1,:), sb(1,:), sc(1,:))
      call strain_coeffs(0.0_wp,1.0_wp,0.0_wp,0.0_wp, Jr, sa(2,:), sb(2,:), sc(2,:))
      call strain_coeffs(0.0_wp,0.0_wp,1.0_wp,0.0_wp, Jr, sa(3,:), sb(3,:), sc(3,:))
      call strain_coeffs(0.0_wp,0.0_wp,0.0_wp,1.0_wp, Jr, sa(4,:), sb(4,:), sc(4,:))
   end subroutine ve_strain_constants

   pure subroutine dissipative_rhs(ne, r, sa, sb, sc, norm, Am, Bm, Cm, f)
      !! Accumulate the dissipative memory forcing −∫ τ^{V}:δε dV into the RHS f
      !! (length ndof). Per element: 2-point radial Gauss quadrature (eqs 94-95)
      !! of the spectral double-dot Σ_λ norm_λ τ^{V,λ}(r) δε^λ(r) (eq 110), with
      !! τ^{V,λ} from the stored memory coefficients (eq 109). Operates on plain
      !! arrays so both the 1-D stepper and the per-(l,m) field driver share it.
      integer,  intent(in)    :: ne
      real(wp), intent(in)    :: r(:), sa(:,:), sb(:,:), sc(:,:), norm(:)
      real(wp), intent(in)    :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      real(wp), intent(inout) :: f(:)
      real(wp), parameter :: xg = 0.5773502691896257_wp   ! 1/√3
      real(wp) :: gp(2)
      real(wp) :: rk, rk1, h, ra, psik, psik1, tauV(NLAM), deps, D, floc(4)
      integer  :: e, ig, t, m
      gp = [ -xg, xg ]
      do e = 1, ne
         rk = r(e);  rk1 = r(e+1);  h = rk1 - rk
         floc = 0.0_wp
         do ig = 1, 2
            ra    = 0.5_wp*(h*gp(ig) + rk + rk1)        ! Gauss node (eq 95)
            psik  = (rk1 - ra)/h
            psik1 = (ra - rk)/h
            do m = 1, NLAM                               ! τ^{V,λ}(ra) (eq 109)
               tauV(m) = Am(m,e)/h + Bm(m,e)*psik/ra + Cm(m,e)*psik1/ra
            end do
            do t = 1, 4                                  ! the 4 local test dofs
               D = 0.0_wp
               do m = 1, NLAM
                  deps = sa(t,m)/h + sb(t,m)*psik/ra + sc(t,m)*psik1/ra
                  D = D + norm(m)*tauV(m)*deps
               end do
               floc(t) = floc(t) - D*ra*ra*h*0.5_wp      ! −D r² h/2 (eq 94, w=1)
            end do
         end do
         f(idx_u(e))   = f(idx_u(e))   + floc(1)
         f(idx_u(e+1)) = f(idx_u(e+1)) + floc(2)
         f(idx_v(e))   = f(idx_v(e))   + floc(3)
         f(idx_v(e+1)) = f(idx_v(e+1)) + floc(4)
      end do
   end subroutine dissipative_rhs

   pure subroutine advance_memory(ne, mu, Mk, Un, Vn, Jr, Am, Bm, Cm)
      !! Explicit Maxwell update of the memory stress from the new nodal strain
      !! (eq 102/107): [A,B,C] ← (1−M)[A,B,C] − 2μM [a,b,c](uⁱ). At the first
      !! step the memory starts at 0, so this sets τ^{V,0} = −2μM ε⁰ (eq 25).
      integer,  intent(in)    :: ne
      real(wp), intent(in)    :: mu(:), Mk(:), Un(:), Vn(:), Jr
      real(wp), intent(inout) :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      real(wp) :: a(NLAM), b(NLAM), c(NLAM), om, two_muM
      integer  :: e, m
      do e = 1, ne
         call strain_coeffs(Un(e), Un(e+1), Vn(e), Vn(e+1), Jr, a, b, c)
         om      = 1.0_wp - Mk(e)
         two_muM = 2.0_wp*mu(e)*Mk(e)
         do m = 1, NLAM
            Am(m,e) = om*Am(m,e) - two_muM*a(m)
            Bm(m,e) = om*Bm(m,e) - two_muM*b(m)
            Cm(m,e) = om*Cm(m,e) - two_muM*c(m)
         end do
      end do
   end subroutine advance_memory

   subroutine ve_destroy(self)
      class(ve_degree), intent(inout) :: self
      call self%op%destroy()
      if (allocated(self%r))    deallocate(self%r)
      if (allocated(self%mu))   deallocate(self%mu)
      if (allocated(self%Mk))   deallocate(self%Mk)
      if (allocated(self%norm)) deallocate(self%norm)
      if (allocated(self%Am))   deallocate(self%Am)
      if (allocated(self%Bm))   deallocate(self%Bm)
      if (allocated(self%Cm))   deallocate(self%Cm)
      if (allocated(self%Un))   deallocate(self%Un)
      if (allocated(self%Vn))   deallocate(self%Vn)
   end subroutine ve_destroy

end module fe_viscoelastic
