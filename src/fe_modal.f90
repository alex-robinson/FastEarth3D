module fe_modal
   !! Modal reduction of the per-degree viscoelastic relaxation (RESP_MODAL).
   !!
   !! For one spherical-harmonic degree the FE memory-stress state relaxes as a
   !! sum of exponential normal modes. This module extracts the few dominant
   !! modes per degree, so the response can be carried as K scalar modal
   !! amplitudes per (l,m) instead of the full radial memory tensor. See
   !! doc/design-modal.md.
   !!
   !! Propagator. The per-degree stepper run with the load OFF maps a memory state
   !! τ → Pτ (fe_viscoelastic's ve_step, σ=0). The BACKWARD-EULER propagator
   !! (SCHEME_BE) has eigenvalues λ = 1/(1+Δt/τ), placing the slowest (dominant)
   !! modes at the largest |λ| → 1, fast modes → 0, so subspace iteration
   !! converges to the physical modes. τ = Δt·λ/(1−λ) recovers the relaxation time
   !! exactly, so Δt affects only conditioning, not accuracy.
   !!
   !! Physical subspace. The FE memory is over-parametrised (3·NLAM coeffs per
   !! element) versus the ~4 strain DOFs an element can carry, so the full memory
   !! operator has spurious (marginally-unstable, zero-residue) eigenmodes that
   !! would swamp the largest-|λ| end. The eigensolve is therefore confined to the
   !! physical (strain-reachable) subspace by iterating in the NODAL (U,V)
   !! displacement parametrisation: a nodal field w maps to memory via
   !! strain_coeffs (StrainGen), and back by least squares (StrainGen⁺ =
   !! Gm⁻¹·StrainGenᵀ, Gm = StrainGenᵀStrainGen). The reduced propagator
   !! P_nodal = StrainGen⁺·P·StrainGen has only the physical modes. This is exact
   !! for uniform viscosity and a benign Galerkin reduction across layer
   !! boundaries.
   !!
   !! Modal form. A step load σ relaxes each amplitude φ_k toward σ with time τ_k,
   !! and u(t) = u_el·σ + Σ_k C^u_k·φ_k, C^u_k = r^u_k·b_k·τ_k (likewise N, V).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_mesh, radial_mesh_build, &
                                 radial_operator_solve_vec, radial_operator_load_rhs, &
                                 idx_u, idx_v, idx_f, ndof_of
   use fe_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, &
                                 strain_coeffs, dissipative_rhs, SCHEME_FE, SCHEME_BE
   implicit none
   private

   public :: modal_degree, modal_spectrum
   public :: modal_degree_init, modal_degree_destroy, modal_apply_p
   public :: modal_solve, modal_spectrum_destroy
   public :: RANK_ISOSTATIC, RANK_RATE, RANK_RESIDUE

   integer, parameter :: RANK_ISOSTATIC = 1   !! |C^u_k| = |r^u·b·τ| (final relaxed uplift)
   integer, parameter :: RANK_RATE      = 2   !! |r^u·b| (initial relaxation rate)
   integer, parameter :: RANK_RESIDUE   = 3   !! |r^u_k| (pure surface coupling)

   type :: modal_degree
      !! Per-degree propagator engine + nodal (strain-subspace) infrastructure.
      !! Wraps a `ve_degree` stepper used as the matrix-free map P (load held at
      !! zero) and the StrainGen / StrainGen⁺ machinery that confines the
      !! eigensolve to the physical subspace. Relaxation lives only on Maxwell
      !! elements; elastic/fluid elements are forced memory-free (Mk = 0).
      integer  :: j  = -1                 !! spherical-harmonic degree
      integer  :: nr = 0, ne = 0          !! radial nodes / elements
      integer  :: nmax = 0                !! number of Maxwell elements
      integer  :: nact = 0, ndn = 0       !! active nodes, nodal DOFs (ndn = 2*nact)
      real(wp) :: dt = 0.0_wp             !! propagator step [s]
      real(wp) :: Jr = 0.0_wp             !! l(l+1)
      logical,  allocatable :: maxwell(:) !! (ne) element carries memory
      integer,  allocatable :: nodepos(:) !! (nr) active-node index, or 0 if inactive
      real(wp), allocatable :: Gm(:,:)    !! (ndn,ndn) LU of StrainGenᵀStrainGen
      integer,  allocatable :: piv(:)     !! (ndn) LU pivots
      type(ve_degree) :: eng              !! the reused per-degree VE stepper
   end type modal_degree

   type :: modal_spectrum
      !! Result of modal_solve for one degree: the selected relaxation modes plus
      !! the instantaneous elastic surface gains. The response to a load history
      !! is u = gu·σ + Σ_k Cu_k·φ_k (and likewise N, V), φ_k the first-order lag of
      !! time constant tau(k) driven by σ.
      integer  :: j = -1, nmode = 0
      real(wp) :: gu = 0.0_wp, gn = 0.0_wp, gv = 0.0_wp
      real(wp), allocatable :: tau(:)
      real(wp), allocatable :: Cu(:), Cn(:), Cv(:)
   end type modal_spectrum

contains

   ! ===================================================================
   ! Propagator + nodal infrastructure
   ! ===================================================================

   subroutine modal_degree_init(self, earth, mesh, j, dt, scheme)
      !! Set up the per-degree propagator and the nodal strain-subspace machinery
      !! (active nodes, Gm = StrainGenᵀStrainGen factored). `scheme` selects the
      !! propagator: SCHEME_BE (default) for the eigensolve, SCHEME_FE for the
      !! load-forcing probe. Non-Maxwell elements are forced memory-free.
      type(modal_degree), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(radial_mesh),  intent(in)    :: mesh
      integer,            intent(in)    :: j
      real(wp),           intent(in)    :: dt
      integer, optional,  intent(in)    :: scheme
      integer  :: sch, e, n, ia
      logical, allocatable :: act(:)

      sch = SCHEME_BE;  if (present(scheme)) sch = scheme
      call modal_degree_destroy(self)
      self%j  = j;  self%dt = dt
      self%nr = mesh%nr;  self%ne = mesh%ne
      call ve_init(self%eng, earth, mesh, j, dt)
      self%eng%scheme = sch
      self%eng%max_couple_iter = 100
      self%eng%couple_tol      = 1.0e-13_wp
      self%Jr = self%eng%Jr

      ! Maxwell mask (by rheology); zero Mk on non-Maxwell (frozen, no memory).
      allocate(self%maxwell(self%ne))
      self%nmax = 0
      do e = 1, self%ne
         self%maxwell(e) = (earth%layers(mesh%elem_layer(e))%rheology == RHEOL_MAXWELL)
         if (self%maxwell(e)) then
            self%nmax = self%nmax + 1
         else
            self%eng%Mk(e) = 0.0_wp
         end if
      end do

      ! Active nodes: endpoints of any Maxwell element. Number them for the nodal w.
      allocate(act(self%nr));  act = .false.
      do e = 1, self%ne
         if (self%maxwell(e)) then;  act(e) = .true.;  act(e+1) = .true.;  end if
      end do
      allocate(self%nodepos(self%nr));  self%nodepos = 0
      ia = 0
      do n = 1, self%nr
         if (act(n)) then;  ia = ia + 1;  self%nodepos(n) = ia;  end if
      end do
      self%nact = ia;  self%ndn = 2*ia

      call build_gm(self)
   end subroutine modal_degree_init

   subroutine modal_apply_p(self, inA, inB, inC, outA, outB, outC, sigma)
      !! Apply the memory-space relaxation propagator: out = P·in (load `sigma`,
      !! default 0). One step of the VE stepper — reuses ve_step verbatim.
      type(modal_degree), intent(inout) :: self
      real(wp),           intent(in)    :: inA(:,:), inB(:,:), inC(:,:)
      real(wp),           intent(out)   :: outA(:,:), outB(:,:), outC(:,:)
      real(wp), optional, intent(in)    :: sigma
      real(wp) :: t, ua, va, fa, sig
      sig = 0.0_wp;  if (present(sigma)) sig = sigma
      self%eng%Am = inA;  self%eng%Bm = inB;  self%eng%Cm = inC
      self%eng%time = 0.0_wp
      self%eng%Un_prev = 0.0_wp;  self%eng%Vn_prev = 0.0_wp
      call ve_step(self%eng, sig, t, ua, va, fa)
      outA = self%eng%Am;  outB = self%eng%Bm;  outC = self%eng%Cm
   end subroutine modal_apply_p

   subroutine modal_degree_destroy(self)
      type(modal_degree), intent(inout) :: self
      call ve_destroy(self%eng)
      if (allocated(self%maxwell)) deallocate(self%maxwell)
      if (allocated(self%nodepos)) deallocate(self%nodepos)
      if (allocated(self%Gm))      deallocate(self%Gm)
      if (allocated(self%piv))     deallocate(self%piv)
      self%j = -1;  self%nr = 0;  self%ne = 0;  self%nmax = 0
      self%nact = 0;  self%ndn = 0;  self%dt = 0.0_wp
   end subroutine modal_degree_destroy

   ! --- StrainGen: nodal field <-> memory --------------------------------------

   subroutine strain_fwd(self, w, A, B, C)
      !! StrainGen: nodal (U,V) field w (ndn) → memory coefficients (NLAM,ne) via
      !! strain_coeffs on each Maxwell element. Non-Maxwell elements zeroed.
      type(modal_degree), intent(in)  :: self
      real(wp),           intent(in)  :: w(:)
      real(wp),           intent(out) :: A(:,:), B(:,:), C(:,:)
      real(wp) :: aa(NLAM), bb(NLAM), cc(NLAM), u1, u2, v1, v2
      integer  :: e, pe, pe1
      A = 0.0_wp;  B = 0.0_wp;  C = 0.0_wp
      do e = 1, self%ne
         if (.not. self%maxwell(e)) cycle
         pe = self%nodepos(e);  pe1 = self%nodepos(e+1)
         u1 = w(2*pe-1);  u2 = w(2*pe1-1);  v1 = w(2*pe);  v2 = w(2*pe1)
         call strain_coeffs(u1, u2, v1, v2, self%Jr, aa, bb, cc)
         A(:,e) = aa;  B(:,e) = bb;  C(:,e) = cc
      end do
   end subroutine strain_fwd

   subroutine strain_adj(self, A, B, C, g)
      !! StrainGenᵀ: memory (NLAM,ne) → nodal field g (ndn). Exact transpose of
      !! strain_fwd (S_eᵀ scattered over Maxwell elements).
      type(modal_degree), intent(in)  :: self
      real(wp),           intent(in)  :: A(:,:), B(:,:), C(:,:)
      real(wp),           intent(out) :: g(:)
      real(wp) :: g1, g2, g3, g4, Jr
      integer  :: e, pe, pe1
      Jr = self%Jr;  g = 0.0_wp
      do e = 1, self%ne
         if (.not. self%maxwell(e)) cycle
         pe = self%nodepos(e);  pe1 = self%nodepos(e+1)
         g1 = -A(1,e) + B(2,e) - B(3,e)/Jr
         g2 =  A(1,e) + C(2,e) - C(3,e)/Jr
         g3 = -A(2,e) - B(2,e) + 0.5_wp*(B(3,e) + B(4,e))
         g4 =  A(2,e) - C(2,e) + 0.5_wp*(C(3,e) + C(4,e))
         g(2*pe-1)  = g(2*pe-1)  + g1
         g(2*pe1-1) = g(2*pe1-1) + g2
         g(2*pe)    = g(2*pe)    + g3
         g(2*pe1)   = g(2*pe1)   + g4
      end do
   end subroutine strain_adj

   subroutine build_gm(self)
      !! Assemble Gm = StrainGenᵀStrainGen (ndn×ndn) and LU-factor it. Built by
      !! applying strain_adj∘strain_fwd to each nodal basis vector, so it is the
      !! exact Gram matrix of the StrainGen used by the iteration.
      type(modal_degree), intent(inout) :: self
      real(wp), allocatable :: w(:), A(:,:), B(:,:), C(:,:), g(:)
      integer :: i
      allocate(self%Gm(self%ndn, self%ndn), self%piv(self%ndn))
      allocate(w(self%ndn), A(NLAM,self%ne), B(NLAM,self%ne), C(NLAM,self%ne), g(self%ndn))
      do i = 1, self%ndn
         w = 0.0_wp;  w(i) = 1.0_wp
         call strain_fwd(self, w, A, B, C)
         call strain_adj(self, A, B, C, g)
         self%Gm(:,i) = g
      end do
      call lu_factor(self%Gm, self%piv, self%ndn)
   end subroutine build_gm

   subroutine strain_pinv(self, A, B, C, w)
      !! StrainGen⁺ = Gm⁻¹·StrainGenᵀ : memory → nodal field (least-squares).
      type(modal_degree), intent(in)  :: self
      real(wp),           intent(in)  :: A(:,:), B(:,:), C(:,:)
      real(wp),           intent(out) :: w(:)
      call strain_adj(self, A, B, C, w)
      call lu_solve(self%Gm, self%piv, w, self%ndn)
   end subroutine strain_pinv

   subroutine apply_p_nodal(self, scrA, scrB, scrC, oA, oB, oC, win, wout)
      !! Reduced propagator P_nodal·w = StrainGen⁺·P·StrainGen·w (dim ndn).
      type(modal_degree), intent(inout) :: self
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:), oA(:,:), oB(:,:), oC(:,:)
      real(wp), intent(in)    :: win(:)
      real(wp), intent(out)   :: wout(:)
      call strain_fwd(self, win, scrA, scrB, scrC)
      call modal_apply_p(self, scrA, scrB, scrC, oA, oB, oC)   ! BE, σ=0
      call strain_pinv(self, oA, oB, oC, wout)
   end subroutine apply_p_nodal

   ! ===================================================================
   ! Main solve
   ! ===================================================================

   subroutine modal_solve(spec, earth, mesh, j, n_modes, mode_rank, dt_be, &
                          p_block, tol, maxit)
      !! Extract the dominant relaxation modes for degree j (≥1) into `spec` by
      !! block subspace iteration on the nodal BE propagator + Rayleigh–Ritz, then
      !! rank and truncate. n_modes ≤ 0 keeps all significant modes.
      type(modal_spectrum), intent(out) :: spec
      type(earth_model),    intent(in)  :: earth
      type(radial_mesh),    intent(in)  :: mesh
      integer,              intent(in)  :: j
      integer,  optional,   intent(in)  :: n_modes, mode_rank, p_block, maxit
      real(wp), optional,   intent(in)  :: dt_be, tol

      type(modal_degree) :: md
      real(wp), allocatable :: Q(:,:), Hm(:,:), Bnod(:), Bhat(:), bcoef(:)
      real(wp), allocatable :: evr(:), Wh(:,:), Whinv(:,:), wk(:)
      real(wp), allocatable :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), allocatable :: BmA(:,:), BmB(:,:), BmC(:,:)
      real(wp), allocatable :: tau_a(:), ru_a(:), cu_a(:), cn_a(:), cv_a(:), strength(:)
      integer,  allocatable :: ord(:)
      logical,  allocatable :: keep(:)
      real(wp) :: dtbe, lam, ru, rn, rv, bk, g, dtfe
      integer  :: nmreq, rank, pb, ndn, p, nm, k, nkeep, kk, nout

      nmreq = -1;             if (present(n_modes))   nmreq = n_modes
      rank  = RANK_ISOSTATIC; if (present(mode_rank)) rank  = mode_rank
      pb    = 16;             if (present(p_block))   pb    = p_block
      dtbe  = 1.0e3_wp*3.15576e7_wp
      if (present(dt_be)) dtbe = dt_be
      if (present(maxit) .or. present(tol)) continue   ! reserved (Arnoldi is direct)

      call modal_degree_init(md, earth, mesh, j, dtbe, SCHEME_BE)
      ndn = md%ndn;  g = md%eng%op%g_surf
      p = min(pb, ndn)
      allocate(scrA(NLAM,md%ne), scrB(NLAM,md%ne), scrC(NLAM,md%ne))

      ! load forcing B (nodal) — the Krylov seed and the modal load projection.
      ! B_mem = τ₁/Δt from one FE step at τ=0, σ=1 (exact); project to the nodal
      ! strain subspace.
      dtfe = dtbe
      allocate(BmA(NLAM,md%ne), BmB(NLAM,md%ne), BmC(NLAM,md%ne), Bnod(ndn))
      call load_forcing(earth, mesh, j, dtfe, BmA, BmB, BmC)
      call strain_pinv(md, BmA, BmB, BmC, Bnod)

      ! Krylov (Arnoldi) subspace from the load forcing under P_nodal: targets the
      ! load-controllable modes — the surface response only sees these, so the
      ! unobservable near-marginal modes that swamp plain subspace iteration drop
      ! out. Hm is the Arnoldi Hessenberg = QᵀP_nodalQ; nm ≤ p is the Krylov dim.
      allocate(Q(ndn,p), Hm(p,p))
      call arnoldi(md, scrA, scrB, scrC, Bnod, Q, Hm, p, nm)

      ! Rayleigh–Ritz: real eigenpairs of the small Hm
      allocate(evr(nm), Wh(nm,nm), Whinv(nm,nm), Bhat(nm), bcoef(nm))
      call small_eigen(Hm(1:nm,1:nm), nm, evr, Wh)
      call invert_matrix(Wh, Whinv, nm)
      Bhat  = matmul(transpose(Q(:,1:nm)), Bnod)
      bcoef = matmul(Whinv, Bhat)

      ! per-mode τ, residues, strengths
      allocate(wk(ndn), tau_a(nm), ru_a(nm), cu_a(nm), cn_a(nm), cv_a(nm), keep(nm))
      do k = 1, nm
         lam = evr(k)
         keep(k) = (lam > 1.0e-12_wp .and. lam < 1.0_wp - 1.0e-12_wp)
         if (.not. keep(k)) cycle
         tau_a(k) = dtbe * lam / (1.0_wp - lam)
         wk = matmul(Q(:,1:nm), Wh(:,k))
         call mode_residue(md, scrA, scrB, scrC, wk, g, j, ru, rn, rv)
         bk = bcoef(k)
         ru_a(k) = ru
         cu_a(k) = ru * bk * tau_a(k)
         cn_a(k) = rn * bk * tau_a(k)
         cv_a(k) = rv * bk * tau_a(k)
      end do

      ! rank and select
      allocate(strength(nm), ord(nm))
      do k = 1, nm
         if (.not. keep(k)) then;  strength(k) = -1.0_wp;  cycle;  end if
         select case (rank)
         case (RANK_RATE);    strength(k) = abs(ru_a(k)*bcoef(k))
         case (RANK_RESIDUE); strength(k) = abs(ru_a(k))
         case default;        strength(k) = abs(cu_a(k))
         end select
      end do
      call sort_desc(strength, ord, nm)

      nkeep = count(keep)
      if (nmreq > 0) nkeep = min(nkeep, nmreq)
      nout = 0
      do k = 1, nkeep
         kk = ord(k)
         if (strength(kk) <= 0.0_wp) exit
         if (k > 1 .and. strength(kk) < 1.0e-6_wp*strength(ord(1))) exit
         nout = nout + 1
      end do

      call elastic_gains(md, j, g, spec%gu, spec%gn, spec%gv)
      spec%j = j;  spec%nmode = nout
      allocate(spec%tau(nout), spec%Cu(nout), spec%Cn(nout), spec%Cv(nout))
      do k = 1, nout
         kk = ord(k)
         spec%tau(k) = tau_a(kk);  spec%Cu(k) = cu_a(kk)
         spec%Cn(k)  = cn_a(kk);   spec%Cv(k) = cv_a(kk)
      end do

      call modal_degree_destroy(md)
   end subroutine modal_solve

   subroutine modal_spectrum_destroy(spec)
      type(modal_spectrum), intent(inout) :: spec
      if (allocated(spec%tau)) deallocate(spec%tau)
      if (allocated(spec%Cu))  deallocate(spec%Cu)
      if (allocated(spec%Cn))  deallocate(spec%Cn)
      if (allocated(spec%Cv))  deallocate(spec%Cv)
      spec%nmode = 0;  spec%j = -1
   end subroutine modal_spectrum_destroy

   ! ===================================================================
   ! Solve internals
   ! ===================================================================

   subroutine arnoldi(md, scrA, scrB, scrC, b0, Q, Hm, p, nm)
      !! Arnoldi process: build the Krylov basis Q of span{b0, P_nodal·b0, …} up to
      !! p vectors, with Hm the upper-Hessenberg projection QᵀP_nodalQ. Returns the
      !! actual dimension nm (≤ p; smaller if an invariant subspace is reached).
      !! Modified Gram–Schmidt with one re-orthogonalisation pass for stability.
      type(modal_degree), intent(inout) :: md
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), intent(in)    :: b0(:)
      real(wp), intent(out)   :: Q(:,:), Hm(:,:)
      integer,  intent(in)    :: p
      integer,  intent(out)   :: nm
      real(wp), allocatable :: w(:), oA(:,:), oB(:,:), oC(:,:)
      real(wp) :: beta, h, r
      integer  :: ndn, k, i, pass
      ndn = size(Q,1)
      allocate(w(ndn), oA(NLAM,md%ne), oB(NLAM,md%ne), oC(NLAM,md%ne))
      Hm = 0.0_wp
      beta = sqrt(dot_product(b0, b0))
      if (beta <= 1.0e-300_wp) then;  nm = 1;  Q(:,1) = 0.0_wp;  return;  end if
      Q(:,1) = b0/beta
      nm = p
      do k = 1, p
         call apply_p_nodal(md, scrA, scrB, scrC, oA, oB, oC, Q(:,k), w)
         do pass = 1, 2
            do i = 1, k
               r = dot_product(Q(:,i), w)
               Hm(i,k) = Hm(i,k) + r
               w = w - r*Q(:,i)
            end do
         end do
         if (k == p) exit
         h = sqrt(dot_product(w, w))
         if (h <= 1.0e-12_wp*beta) then;  nm = k;  exit;  end if
         Hm(k+1,k) = h
         Q(:,k+1) = w/h
      end do
   end subroutine arnoldi

   subroutine load_forcing(earth, mesh, j, dtfe, BmA, BmB, BmC)
      !! Memory-space load-forcing B = τ₁/Δt from one FE step at τ=0, σ=1 (exact:
      !! from τ=0, FE gives τ₁ = Δt·B·σ).
      type(earth_model), intent(in)  :: earth
      type(radial_mesh), intent(in)  :: mesh
      integer,           intent(in)  :: j
      real(wp),          intent(in)  :: dtfe
      real(wp),          intent(out) :: BmA(:,:), BmB(:,:), BmC(:,:)
      type(modal_degree) :: fe
      real(wp), allocatable :: zA(:,:), zB(:,:), zC(:,:)
      call modal_degree_init(fe, earth, mesh, j, dtfe, SCHEME_FE)
      allocate(zA(NLAM,fe%ne), zB(NLAM,fe%ne), zC(NLAM,fe%ne))
      zA = 0.0_wp;  zB = 0.0_wp;  zC = 0.0_wp
      call modal_apply_p(fe, zA, zB, zC, BmA, BmB, BmC, sigma=1.0_wp)
      BmA = BmA/dtfe;  BmB = BmB/dtfe;  BmC = BmC/dtfe
      call modal_degree_destroy(fe)
   end subroutine load_forcing

   subroutine mode_residue(md, scrA, scrB, scrC, wk, g, j, ru, rn, rv)
      !! Surface drift produced by the nodal mode shape wk: build its memory
      !! (StrainGen), solve K⁻¹·(dissipative forcing), read surface coefficients.
      type(modal_degree), intent(inout) :: md
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), intent(in)    :: wk(:), g
      integer,  intent(in)    :: j
      real(wp), intent(out)   :: ru, rn, rv
      real(wp), allocatable :: f(:), x(:)
      integer :: nd
      nd = ndof_of(md%nr)
      allocate(f(nd), x(nd))
      call strain_fwd(md, wk, scrA, scrB, scrC)
      f = 0.0_wp
      call dissipative_rhs(md%eng%ne, md%eng%r, md%eng%sa, md%eng%sb, md%eng%sc, &
                           md%eng%norm, scrA, scrB, scrC, f)
      call radial_operator_solve_vec(md%eng%op, f, x)
      ru = x(idx_u(md%nr))
      rn = -x(idx_f(md%nr))/g
      rv = x(idx_v(md%nr))
      if (j == 1) rn = 0.0_wp
   end subroutine mode_residue

   subroutine elastic_gains(md, j, g, gu, gn, gv)
      !! Instantaneous elastic surface gains: unit-load solve of the operator.
      type(modal_degree), intent(inout) :: md
      integer,  intent(in)  :: j
      real(wp), intent(in)  :: g
      real(wp), intent(out) :: gu, gn, gv
      real(wp), allocatable :: x(:)
      allocate(x(ndof_of(md%nr)))
      call radial_operator_solve_vec(md%eng%op, radial_operator_load_rhs(md%eng%op, 1.0_wp), x)
      gu = x(idx_u(md%nr))
      gn = -x(idx_f(md%nr))/g
      gv = x(idx_v(md%nr))
      if (j == 1) gn = 0.0_wp
   end subroutine elastic_gains

   ! ===================================================================
   ! Small dense linear algebra (in-house; runs once per degree)
   ! ===================================================================

   subroutine small_eigen(Hin, n, evr, evec)
      !! Real eigenvalues (evr, descending) and right eigenvectors (evec) of the
      !! small real matrix Hin: shifted-QR with deflation, then inverse iteration.
      real(wp), intent(in)  :: Hin(:,:)
      integer,  intent(in)  :: n
      real(wp), intent(out) :: evr(:)
      real(wp), intent(out) :: evec(:,:)
      real(wp), allocatable :: H(:,:), Qf(:,:), Rf(:,:)
      real(wp) :: mu, aa, bb, cc, dd, tr, det, disc, r1, r2, sub
      integer  :: m, sweep, i, k
      integer, parameter :: MAXSWEEP = 3000
      allocate(H(n,n), Qf(n,n), Rf(n,n))
      H = Hin(1:n,1:n)
      m = n
      do while (m > 1)
         do sweep = 1, MAXSWEEP
            sub = abs(H(m,m-1))
            if (sub <= 1.0e-14_wp*(abs(H(m-1,m-1)) + abs(H(m,m)) + tiny(1.0_wp))) exit
            aa = H(m-1,m-1);  bb = H(m-1,m);  cc = H(m,m-1);  dd = H(m,m)
            tr = aa + dd;  det = aa*dd - bb*cc;  disc = 0.25_wp*tr*tr - det
            if (disc >= 0.0_wp) then
               r1 = 0.5_wp*tr + sqrt(disc);  r2 = 0.5_wp*tr - sqrt(disc)
               if (abs(r1-dd) <= abs(r2-dd)) then;  mu = r1;  else;  mu = r2;  end if
            else
               mu = dd
            end if
            do i = 1, m;  H(i,i) = H(i,i) - mu;  end do
            call qr_square(H(1:m,1:m), Qf(1:m,1:m), Rf(1:m,1:m), m)
            H(1:m,1:m) = matmul(Rf(1:m,1:m), Qf(1:m,1:m))
            do i = 1, m;  H(i,i) = H(i,i) + mu;  end do
         end do
         evr(m) = H(m,m)
         m = m - 1
      end do
      evr(1) = H(1,1)
      call sort_real_desc(evr, n)
      do k = 1, n
         call inverse_iteration(Hin, n, evr(k), evec(:,k))
      end do
   end subroutine small_eigen

   subroutine qr_square(A, Q, R, n)
      real(wp), intent(in)  :: A(:,:)
      real(wp), intent(out) :: Q(:,:), R(:,:)
      integer,  intent(in)  :: n
      integer :: i, jj
      real(wp) :: nrm
      Q(1:n,1:n) = A(1:n,1:n);  R(1:n,1:n) = 0.0_wp
      do jj = 1, n
         do i = 1, jj-1
            R(i,jj) = dot_product(Q(1:n,i), Q(1:n,jj))
            Q(1:n,jj) = Q(1:n,jj) - R(i,jj)*Q(1:n,i)
         end do
         nrm = sqrt(dot_product(Q(1:n,jj), Q(1:n,jj)))
         if (nrm <= 1.0e-300_wp) nrm = 1.0_wp
         R(jj,jj) = nrm;  Q(1:n,jj) = Q(1:n,jj)/nrm
      end do
   end subroutine qr_square

   subroutine inverse_iteration(H, n, lambda, v)
      real(wp), intent(in)  :: H(:,:)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: lambda
      real(wp), intent(out) :: v(:)
      real(wp), allocatable :: M(:,:), b(:), w(:)
      integer,  allocatable :: piv(:)
      real(wp) :: eps, nrm
      integer  :: i, it
      allocate(M(n,n), b(n), w(n), piv(n))
      eps = 1.0e-10_wp*(abs(lambda) + 1.0_wp)
      M = H(1:n,1:n)
      do i = 1, n;  M(i,i) = M(i,i) - (lambda + eps);  end do
      call lu_factor(M, piv, n)
      b = 1.0_wp/sqrt(real(n,wp))
      do it = 1, 5
         w = b
         call lu_solve(M, piv, w, n)
         nrm = sqrt(dot_product(w, w))
         if (nrm <= 1.0e-300_wp) nrm = 1.0_wp
         b = w/nrm
      end do
      v(1:n) = b
   end subroutine inverse_iteration

   subroutine invert_matrix(A, Ainv, n)
      real(wp), intent(in)  :: A(:,:)
      real(wp), intent(out) :: Ainv(:,:)
      integer,  intent(in)  :: n
      real(wp), allocatable :: M(:,:), e(:)
      integer,  allocatable :: piv(:)
      integer :: k
      allocate(M(n,n), e(n), piv(n))
      M = A(1:n,1:n)
      call lu_factor(M, piv, n)
      do k = 1, n
         e = 0.0_wp;  e(k) = 1.0_wp
         call lu_solve(M, piv, e, n)
         Ainv(1:n,k) = e
      end do
   end subroutine invert_matrix

   subroutine lu_factor(A, piv, n)
      !! In-place LU with partial pivoting (Doolittle).
      real(wp), intent(inout) :: A(:,:)
      integer,  intent(out)   :: piv(:)
      integer,  intent(in)    :: n
      integer  :: i, k, ip
      real(wp) :: big, f
      real(wp), allocatable :: tmp(:)
      allocate(tmp(n))
      do k = 1, n
         ip = k;  big = abs(A(k,k))
         do i = k+1, n
            if (abs(A(i,k)) > big) then;  big = abs(A(i,k));  ip = i;  end if
         end do
         piv(k) = ip
         if (ip /= k) then;  tmp = A(k,1:n);  A(k,1:n) = A(ip,1:n);  A(ip,1:n) = tmp;  end if
         if (abs(A(k,k)) <= 1.0e-300_wp) A(k,k) = 1.0e-300_wp
         do i = k+1, n
            f = A(i,k)/A(k,k);  A(i,k) = f
            A(i,k+1:n) = A(i,k+1:n) - f*A(k,k+1:n)
         end do
      end do
   end subroutine lu_factor

   subroutine lu_solve(A, piv, b, n)
      !! Solve LU x = b in place (b → x) given lu_factor output.
      real(wp), intent(in)    :: A(:,:)
      integer,  intent(in)    :: piv(:)
      real(wp), intent(inout) :: b(:)
      integer,  intent(in)    :: n
      integer  :: i, k
      real(wp) :: s
      do k = 1, n
         if (piv(k) /= k) then;  s = b(k);  b(k) = b(piv(k));  b(piv(k)) = s;  end if
      end do
      do i = 2, n
         b(i) = b(i) - dot_product(A(i,1:i-1), b(1:i-1))
      end do
      do i = n, 1, -1
         s = b(i)
         if (i < n) s = s - dot_product(A(i,i+1:n), b(i+1:n))
         b(i) = s/A(i,i)
      end do
   end subroutine lu_solve

   subroutine sort_real_desc(a, n)
      real(wp), intent(inout) :: a(:)
      integer,  intent(in)    :: n
      integer :: i, k
      real(wp) :: v
      do i = 2, n
         v = a(i);  k = i-1
         do while (k >= 1)
            if (a(k) >= v) exit
            a(k+1) = a(k);  k = k-1
         end do
         a(k+1) = v
      end do
   end subroutine sort_real_desc

   subroutine sort_desc(val, ord, n)
      real(wp), intent(in)  :: val(:)
      integer,  intent(out) :: ord(:)
      integer,  intent(in)  :: n
      integer :: i, k, t
      do i = 1, n;  ord(i) = i;  end do
      do i = 2, n
         t = ord(i);  k = i-1
         do while (k >= 1)
            if (val(ord(k)) >= val(t)) exit
            ord(k+1) = ord(k);  k = k-1
         end do
         ord(k+1) = t
      end do
   end subroutine sort_desc

end module fe_modal
