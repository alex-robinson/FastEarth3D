module fe_rotation
   !! Rotational feedback / true polar wander (Spada et al. 2011 §2.1.1; the
   !! time-domain rotational theory of Martinec & Hagedoorn 2014, as in VILMA).
   !!
   !! A surface load and the deformation it drives perturb the off-diagonal
   !! inertia (I₁₃, I₂₃); the linearized Liouville equation maps that to equatorial
   !! polar motion m = m₁ + i m₂; the shifted pole perturbs the centrifugal
   !! potential — a degree-2, order-1 potential — which deforms the Earth back
   !! (tidal response) and feeds the geoid / sea-level equation. We integrate the
   !! GIA (quasi-static) Liouville equation with the Chandler wobble neglected
   !! (Spada eq. 7), the regime the GIA timescale lives in:
   !!
   !!     [1 − k^T(t)/k_s] ∗ m(t) = Ψ_L(t),                              (Spada 7)
   !!     Ψ_L(t) = I(t)/(C−A),  I(t) = [δ(t)+k^L(t)] ∗ I_rigid(t),       (Spada 20)
   !!
   !! with k^T, k^L the degree-2 tidal / loading Love numbers and k_s ≡ k^T_f the
   !! secular (fluid) tidal Love number (Spada eq. 11). NO explicit Ω appears —
   !! the centrifugal scaling is absorbed into m (= ω/Ω) and k_s by construction.
   !!
   !! Method (composes with the time-domain Maxwell machinery, no normal modes,
   !! no convolution quadrature). Two compact degree-2 viscoelastic channels carry
   !! the convolutions as Maxwell memory (reusing the per-element kernel of
   !! fe_viscoelastic):
   !!   - a LOADING channel: forced by the rigid inertia I_rigid, returns
   !!     I(t) = [1+k^L]∗I_rigid  ⇒  Ψ_L = I/(C−A);
   !!   - a TIDAL channel: forced by the centrifugal potential ∝ m, returns the
   !!     induced potential k^T∗m  ⇒  the rotational feedback.
   !! Each step the feedback makes the Liouville equation ALGEBRAIC in m (the affine
   !! begin/apply/commit structure of the field driver):
   !!
   !!     m_n = [ Ψ_L,n − dF_tidal/k_s ] / [ 1 − k^T_e/k_s ],
   !!
   !! where k^T_e is the elastic tidal Love number and dF_tidal is the tidal
   !! channel's frozen memory (the m-history). The memory is then advanced with the
   !! converged m_n. Rotation is purely degree 2, so this is a single complex
   !! coefficient per channel — self-contained, decoupled from the SLE field driver
   !! (the centrifugal potential is fed back into the SLE in a later step).
   !!
   !! 3-D ready: I_rigid is a direct Gauss-grid quadrature of the actual load
   !! (any field), so no axisymmetric assumption enters; only the (1+k^L)/k^T
   !! channels use the 1-D radial relaxation (laterally-varying η is rung 6).
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G
   use fe_earth_structure, only: earth_model, RHEOL_FLUID
   use fe_radial_fe,       only: radial_mesh, radial_operator, tidal_love, &
                                 idx_u, idx_v, idx_f, ndof_of
   use fe_viscoelastic,    only: NLAM, ve_strain_constants, dissipative_rhs, &
                                 advance_memory, SCHEME_FE
   use fe_sht,             only: sht_grid
   implicit none
   private

   public :: rotation_state

   integer, parameter :: JROT = 2          !! rotation is purely degree 2

   type :: deg2_channel
      !! A single degree-2 viscoelastic response to a COMPLEX forcing coefficient,
      !! carrying Maxwell memory. `tidal` selects the forcing kind (external
      !! potential via tidal_rhs, vs. surface load via load_rhs); the surface
      !! perturbed-potential coefficient F(a) is the readout. Affine in the current
      !! forcing: F = Fe·coeff + dF, with dF frozen from the entering memory each
      !! step (begin_step) and the memory advanced with the converged coeff (commit).
      logical  :: tidal = .false.
      integer  :: nr = 0, ne = 0, ndof = 0
      real(wp) :: Jr = 0.0_wp, dt = 0.0_wp
      type(radial_operator) :: op                 !! degree-2 operator (assembled once)
      real(wp), allocatable :: r(:), mu(:), Mk(:), MkPerDt(:)
      real(wp) :: norm(NLAM), sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)
      real(wp) :: Fe = 0.0_wp                      !! elastic surface-F per unit forcing
      real(wp), allocatable :: xUn(:), xVn(:)      !! unit-forcing nodal U,V
      ! per-element memory stress, real/imag (NLAM, ne)
      real(wp), allocatable :: Are(:,:), Aim(:,:), Bre(:,:), Bim(:,:), Cre(:,:), Cim(:,:)
      ! frozen drift from the entering memory τ_n (begin_step)
      complex(wp) :: dF = (0.0_wp, 0.0_wp)         !! surface F drift
      real(wp), allocatable :: dUn_re(:), dUn_im(:), dVn_re(:), dVn_im(:)  !! nodal ε_n drift
   contains
      procedure :: init       => channel_init
      procedure :: set_dt     => channel_set_dt
      procedure :: begin_step => channel_begin
      procedure :: commit     => channel_commit
      procedure :: destroy    => channel_destroy
   end type deg2_channel

   type :: rotation_state
      logical     :: enabled = .false.       !! off until validated against benchmark
      complex(wp) :: m = (0.0_wp, 0.0_wp)     !! polar motion m₁ + i m₂ [rad]
      real(wp)    :: time = 0.0_wp            !! model time [s]
      ! physics constants (defaults: Spada 2011 Table 2; overridable for deep time)
      real(wp)    :: a       = 6.371e6_wp     !! Earth radius [m]
      real(wp)    :: CminusA = 2.63e35_wp     !! C − A [kg m²]
      real(wp)    :: k_s     = 0.0_wp         !! secular (fluid) tidal Love number k^T_f
      logical     :: k_s_set = .false.        !! .true. if k_s was supplied (deep-time override)
      real(wp)    :: kTe     = 0.0_wp         !! elastic tidal Love number k^T_e (degree 2)
      type(deg2_channel) :: load_ch          !! (1+k^L)∗ channel
      type(deg2_channel) :: tidal_ch         !! k^T∗ channel
   contains
      procedure :: init    => rotation_init
      procedure :: update  => rotation_update
      procedure :: destroy => rotation_destroy
   end type rotation_state

contains

   ! === rotation_state ========================================================

   subroutine rotation_init(self, earth, sht, dt, k_s)
      !! Build the two degree-2 channels, the elastic tidal Love number k^T_e, and
      !! the secular k_s = k^T_f (the relaxed tidal limit = elastic tidal solve of
      !! the model with every Maxwell layer fluidized). Pass k_s to override it with
      !! an observed-flattening value (Mitrovica et al. 2005) for deep-time runs.
      class(rotation_state), intent(inout) :: self
      type(earth_model),     intent(in)    :: earth
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: dt
      real(wp), optional,    intent(in)    :: k_s
      type(radial_mesh) :: mesh

      call self%destroy()
      call mesh%build(earth)
      self%a = earth%r_earth
      call self%load_ch%init(earth, mesh, dt, tidal=.false.)
      call self%tidal_ch%init(earth, mesh, dt, tidal=.true.)
      ! elastic tidal Love number from the tidal channel's unit response:
      ! k^T = −F(a)/φ_t − 1 (tidal_love convention), φ_t = 1 ⇒ k^T_e = −Fe − 1.
      self%kTe = -self%tidal_ch%Fe - 1.0_wp
      if (present(k_s)) then
         self%k_s = k_s;  self%k_s_set = .true.
      else
         self%k_s = fluid_tidal_k(earth, mesh)
         self%k_s_set = .false.
      end if
      self%m = (0.0_wp, 0.0_wp);  self%time = 0.0_wp
   end subroutine rotation_init

   subroutine rotation_update(self, sht, load, dt)
      !! Advance the polar motion one step under the surface mass load `load`
      !! [kg m⁻²] on the Gauss grid, returning the time BEFORE advancing in
      !! self%time-aligned fashion (first call ⇒ elastic m₀). Solves the algebraic
      !! Liouville equation with the tidal feedback frozen, then advances both
      !! channels' memory.
      class(rotation_state), intent(inout) :: self
      type(sht_grid),        intent(in)    :: sht
      real(wp),              intent(in)    :: load(:,:)   !! surface mass load [kg m⁻²]
      real(wp),              intent(in)    :: dt
      complex(wp) :: Irig, cload, Itot, psiL, m_new
      real(wp)    :: scl

      if (.not. self%enabled) return
      if (dt /= self%load_ch%dt) then
         call self%load_ch%set_dt(dt);  call self%tidal_ch%set_dt(dt)
      end if

      ! (2,1) rigid inertia of the current load (direct grid quadrature).
      Irig = inertia21(sht, load, self%a)

      ! Freeze both channels' relaxation drift from the entering memory τ_n.
      call self%load_ch%begin_step()
      call self%tidal_ch%begin_step()

      ! LOADING: feed the channel a load whose own degree-2 potential equals I_rigid
      ! (φ^L = 4πGa σ/(2j+1) ⇒ σ = I_rigid·(2j+1)/(4πGa)); then −F = [1+k^L]∗I_rigid
      ! = I(t), and Ψ_L = I/(C−A).
      scl   = real(2*JROT+1, wp)/(4.0_wp*pi*grav_G*self%a)
      cload = Irig*scl
      Itot  = -(self%load_ch%Fe*cload + self%load_ch%dF)
      psiL  = Itot/self%CminusA

      ! Liouville (algebraic): the tidal feedback's induced potential is
      ! k^T∗m = k^T_e·m − dF_tidal (units cancel: feed φ_t = m, read −F − φ_t).
      ! m = Ψ_L + (1/k_s)(k^T_e m − dF_tidal) ⇒ solve for m.
      m_new = (psiL - self%tidal_ch%dF/self%k_s) / (1.0_wp - self%kTe/self%k_s)

      ! Advance memory: loading with cload (load-driven), tidal with the converged m.
      call self%load_ch%commit(cload)
      call self%tidal_ch%commit(m_new)

      self%m = m_new
      self%time = self%time + dt
   end subroutine rotation_update

   subroutine rotation_destroy(self)
      class(rotation_state), intent(inout) :: self
      call self%load_ch%destroy()
      call self%tidal_ch%destroy()
      self%m = (0.0_wp, 0.0_wp);  self%time = 0.0_wp
      self%k_s = 0.0_wp;  self%kTe = 0.0_wp;  self%k_s_set = .false.
   end subroutine rotation_destroy

   ! === degree-2 inertia from the load (3-D-ready grid quadrature) =============

   complex(wp) function inertia21(sht, load, a) result(I21)
      !! Off-diagonal inertia perturbation I₁₃ + i I₂₃ of a surface mass load on the
      !! Gauss grid: I₁₃ = −a⁴∫σ sinθcosθ cosφ dΩ, I₂₃ = −a⁴∫σ sinθcosθ sinφ dΩ
      !! (= −a⁴∫σ sinθcosθ e^{iφ} dΩ, the degree-2 order-1 mass moment). Direct
      !! quadrature — no spherical-harmonic normalization enters, so it is exact for
      !! any (3-D) load field.
      type(sht_grid), intent(in) :: sht
      real(wp),       intent(in) :: load(:,:)   !! (nphi, nlat) [kg m⁻²]
      real(wp),       intent(in) :: a
      real(wp), allocatable :: w13(:,:), w23(:,:)
      real(wp) :: st, ct, sc2
      integer  :: il, ip
      allocate(w13(sht%nphi, sht%nlat), w23(sht%nphi, sht%nlat))
      do il = 1, sht%nlat
         st  = sin(sht%colat(il));  ct = cos(sht%colat(il))
         sc2 = st*ct                                   ! sinθ cosθ
         do ip = 1, sht%nphi
            w13(ip,il) = load(ip,il)*sc2*cos(sht%lon(ip))
            w23(ip,il) = load(ip,il)*sc2*sin(sht%lon(ip))
         end do
      end do
      I21 = cmplx(-a**4*sht%surface_integral(w13), &
                  -a**4*sht%surface_integral(w23), wp)
   end function inertia21

   ! === secular (fluid) tidal Love number =====================================

   real(wp) function fluid_tidal_k(earth, mesh) result(k_s)
      !! k_s = k^T_f: the relaxed (t→∞) degree-2 tidal Love number = the ELASTIC
      !! tidal solve of the model with every MAXWELL (viscous) layer fluidized (μ=0).
      !! The elastic lithosphere (RHEOL_ELASTIC, η→∞) is kept elastic and the inviscid
      !! core (RHEOL_FLUID) is unchanged — same construction as the loading fluid limit
      !! in test_benchmark_love. This is Spada eq. 11's secular Love number, and the
      !! rotational secular slope is pathologically sensitive to it (the lithosphere-
      !! thickness paradox, Mitrovica et al. 2005): fluidizing the lithosphere by
      !! mistake inflates k^T_f and badly under-drives the late-time polar motion.
      use fe_earth_structure, only: RHEOL_MAXWELL
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      type(earth_model)     :: ef
      type(radial_operator) :: op
      real(wp), allocatable :: x(:)
      real(wp) :: h, l, ua, va, fa
      integer  :: lay
      ef = earth
      do lay = 1, ef%n_layers()
         if (ef%layers(lay)%rheology == RHEOL_MAXWELL) then
            ef%layers(lay)%mu = 0.0_wp;  ef%layers(lay)%rheology = RHEOL_FLUID
         end if
      end do
      call op%assemble(ef, mesh, JROT)
      allocate(x(op%ndof))
      call op%solve_vec(op%tidal_rhs(1.0_wp), x)
      ua = x(idx_u(mesh%nr));  va = x(idx_v(mesh%nr));  fa = x(idx_f(mesh%nr))
      call tidal_love(ef, JROT, 1.0_wp, ua, va, fa, h, l, k_s)
      call op%destroy()
   end function fluid_tidal_k

   ! === deg2_channel ==========================================================

   subroutine channel_init(self, earth, mesh, dt, tidal)
      !! Assemble the degree-2 operator, the unit-forcing response (Fe + nodal U,V),
      !! the per-element Maxwell factors, and zero the memory. `tidal` picks the
      !! forcing kind (tidal_rhs vs load_rhs).
      class(deg2_channel), intent(inout) :: self
      type(earth_model),   intent(in)    :: earth
      type(radial_mesh),   intent(in)    :: mesh
      real(wp),            intent(in)    :: dt
      logical,             intent(in)    :: tidal
      real(wp), allocatable :: x(:)
      real(wp) :: eta_e
      integer  :: e, lay, node

      call self%destroy()
      self%tidal = tidal
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%Jr = real(JROT, wp)*real(JROT+1, wp);  self%dt = dt

      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne), self%MkPerDt(self%ne))
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         if (eta_e > 0.0_wp) then
            self%MkPerDt(e) = self%mu(e)/eta_e
         else
            self%MkPerDt(e) = 0.0_wp
         end if
      end do
      self%Mk = self%MkPerDt*dt

      call ve_strain_constants(self%Jr, self%norm, self%sa, self%sb, self%sc)
      call self%op%assemble(earth, mesh, JROT)

      ! unit-forcing response: Fe (surface F) + nodal U,V for the memory forcing
      allocate(x(self%ndof), self%xUn(self%nr), self%xVn(self%nr))
      if (tidal) then
         call self%op%solve_vec(self%op%tidal_rhs(1.0_wp), x)
      else
         call self%op%solve_vec(self%op%load_rhs(1.0_wp), x)
      end if
      self%Fe = x(idx_f(self%nr))
      do node = 1, self%nr
         self%xUn(node) = x(idx_u(node));  self%xVn(node) = x(idx_v(node))
      end do

      allocate(self%Are(NLAM,self%ne), self%Aim(NLAM,self%ne))
      allocate(self%Bre(NLAM,self%ne), self%Bim(NLAM,self%ne))
      allocate(self%Cre(NLAM,self%ne), self%Cim(NLAM,self%ne))
      self%Are = 0.0_wp; self%Aim = 0.0_wp; self%Bre = 0.0_wp
      self%Bim = 0.0_wp; self%Cre = 0.0_wp; self%Cim = 0.0_wp
      allocate(self%dUn_re(self%nr), self%dUn_im(self%nr), &
               self%dVn_re(self%nr), self%dVn_im(self%nr))
      self%dUn_re = 0.0_wp; self%dUn_im = 0.0_wp
      self%dVn_re = 0.0_wp; self%dVn_im = 0.0_wp
      self%dF = (0.0_wp, 0.0_wp)
   end subroutine channel_init

   subroutine channel_set_dt(self, dt)
      !! Rescale the Maxwell factor for a new Δt (Mk = (μ/η)·Δt); no re-factor.
      class(deg2_channel), intent(inout) :: self
      real(wp),            intent(in)    :: dt
      self%dt = dt;  self%Mk = self%MkPerDt*dt
   end subroutine channel_set_dt

   subroutine channel_begin(self)
      !! Freeze the drift from the entering memory τ_n: solve the load-free memory
      !! forcing −∫τ^V:δε for the real and imaginary parts, storing the surface F
      !! drift (self%dF) and the nodal strain drift (ε_n term for the commit).
      class(deg2_channel), intent(inout) :: self
      real(wp), allocatable :: fre(:), fim(:), xre(:), xim(:)
      integer :: node
      allocate(fre(self%ndof), fim(self%ndof), xre(self%ndof), xim(self%ndof))
      fre = 0.0_wp;  fim = 0.0_wp
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, self%norm, &
                           self%Are, self%Bre, self%Cre, fre)
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, self%norm, &
                           self%Aim, self%Bim, self%Cim, fim)
      call self%op%solve_vec(fre, xre)
      call self%op%solve_vec(fim, xim)
      self%dF = cmplx(xre(idx_f(self%nr)), xim(idx_f(self%nr)), wp)
      do node = 1, self%nr
         self%dUn_re(node) = xre(idx_u(node));  self%dUn_im(node) = xim(idx_u(node))
         self%dVn_re(node) = xre(idx_v(node));  self%dVn_im(node) = xim(idx_v(node))
      end do
   end subroutine channel_begin

   subroutine channel_commit(self, coeff)
      !! Advance the Maxwell memory (forward-Euler) with the converged forcing
      !! coefficient: total nodal strain = coeff·(unit response) + drift(τ_n).
      class(deg2_channel), intent(inout) :: self
      complex(wp),         intent(in)    :: coeff
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp) :: cr, ci
      integer  :: node
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr))
      cr = real(coeff, wp);  ci = aimag(coeff)
      do node = 1, self%nr
         Ure(node) = cr*self%xUn(node) + self%dUn_re(node)
         Uim(node) = ci*self%xUn(node) + self%dUn_im(node)
         Vre(node) = cr*self%xVn(node) + self%dVn_re(node)
         Vim(node) = ci*self%xVn(node) + self%dVn_im(node)
      end do
      call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr, &
                          self%Are, self%Bre, self%Cre)
      call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr, &
                          self%Aim, self%Bim, self%Cim)
   end subroutine channel_commit

   subroutine channel_destroy(self)
      class(deg2_channel), intent(inout) :: self
      call self%op%destroy()
      if (allocated(self%r))      deallocate(self%r)
      if (allocated(self%mu))     deallocate(self%mu)
      if (allocated(self%Mk))     deallocate(self%Mk)
      if (allocated(self%MkPerDt))deallocate(self%MkPerDt)
      if (allocated(self%xUn))    deallocate(self%xUn)
      if (allocated(self%xVn))    deallocate(self%xVn)
      if (allocated(self%Are))    deallocate(self%Are, self%Aim, self%Bre, &
                                             self%Bim, self%Cre, self%Cim)
      if (allocated(self%dUn_re)) deallocate(self%dUn_re, self%dUn_im, &
                                             self%dVn_re, self%dVn_im)
      self%Fe = 0.0_wp;  self%dF = (0.0_wp, 0.0_wp)
      self%nr = 0;  self%ne = 0;  self%ndof = 0
   end subroutine channel_destroy

end module fe_rotation
