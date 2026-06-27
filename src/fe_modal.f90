module fe_modal
   !! Modal reduction of the per-degree viscoelastic relaxation (RESP_MODAL).
   !!
   !! For one spherical-harmonic degree the FE memory-stress state relaxes as a
   !! sum of exponential normal modes. This module extracts the few dominant
   !! modes {τ_k, residues} per degree by block subspace iteration on the FE
   !! relaxation propagator, so the response can be carried as K scalar modal
   !! amplitudes per (l,m) instead of the full radial memory tensor. See
   !! doc/design-modal.md.
   !!
   !! The propagator is the existing, validated per-degree stepper run with the
   !! load OFF: one homogeneous step of fe_viscoelastic's `ve_step` maps a memory
   !! state τ → Pτ. Subspace iteration uses the BACKWARD-EULER propagator
   !! (SCHEME_BE), whose eigenvalues λ_k = 1/(1+Δt/τ_k) place the slowest
   !! (physically dominant) modes at the largest |λ| → 1, well separated from the
   !! fast modes (λ → 0). The eigenvectors are those of the generator A, and the
   !! relaxation time is recovered exactly as τ_k = Δt·λ_k/(1−λ_k).
   !!
   !! This file currently provides the propagator foundation (modal_degree +
   !! modal_apply_p); the eigensolve and residue extraction build on it.
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_mesh, radial_mesh_build
   use fe_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, &
                                 SCHEME_FE, SCHEME_BE
   implicit none
   private

   public :: modal_degree
   public :: modal_degree_init, modal_degree_destroy, modal_apply_p

   type :: modal_degree
      !! Per-degree relaxation propagator engine. Wraps a `ve_degree` stepper
      !! (assembled operator + Maxwell element data) used as the matrix-free
      !! linear map P: memory → memory with the surface load held at zero. The
      !! memory state per degree is the three radial shape-coefficients A,B,C of
      !! the Maxwell stress, each (NLAM, ne) — relaxation lives only on Maxwell
      !! elements (MkPerDt /= 0); elastic/fluid elements stay identically zero
      !! under P, so the iteration is confined to the Maxwell subspace.
      integer  :: j  = -1                 !! spherical-harmonic degree
      integer  :: nr = 0, ne = 0          !! radial nodes / elements
      real(wp) :: dt = 0.0_wp             !! propagator step [s] (sets the λ↔τ map)
      type(ve_degree) :: eng              !! the reused per-degree VE stepper
   end type modal_degree

contains

   subroutine modal_degree_init(self, earth, mesh, j, dt, scheme)
      !! Set up the per-degree propagator engine for degree j at step dt. The
      !! scheme selects the propagator: SCHEME_BE (default) for the shift-invert
      !! subspace iteration; SCHEME_FE for the plain forward map (testing). The
      !! engine is configured for a held (σ=0) implicit advance — Picard-iterated
      !! to a consistent endpoint for BE.
      type(modal_degree), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(radial_mesh),  intent(in)    :: mesh
      integer,            intent(in)    :: j
      real(wp),           intent(in)    :: dt
      integer, optional,  intent(in)    :: scheme
      integer :: sch

      sch = SCHEME_BE;  if (present(scheme)) sch = scheme
      call modal_degree_destroy(self)
      self%j  = j;  self%dt = dt
      self%nr = mesh%nr;  self%ne = mesh%ne
      call ve_init(self%eng, earth, mesh, j, dt)
      self%eng%scheme = sch
      ! Implicit (BE) propagator needs the within-step coupling iteration to reach
      ! the consistent endpoint memory; explicit (FE) ignores it. A generous cap
      ! with a tight tolerance — this runs once at init, so converge it well.
      self%eng%max_couple_iter = 100
      self%eng%couple_tol      = 1.0e-12_wp
   end subroutine modal_degree_init

   subroutine modal_apply_p(self, inA, inB, inC, outA, outB, outC)
      !! Apply the relaxation propagator: out = P·in, where `in`/`out` are the
      !! Maxwell memory coefficients (NLAM, ne). One homogeneous (σ=0) step of the
      !! VE stepper: the surface report balances against the input memory and the
      !! memory advances by the configured scheme. Reuses `ve_step` verbatim, so
      !! P is exactly the discrete relaxation the full RESP_VE model would apply.
      type(modal_degree), intent(inout) :: self
      real(wp),           intent(in)    :: inA(:,:), inB(:,:), inC(:,:)
      real(wp),           intent(out)   :: outA(:,:), outB(:,:), outC(:,:)
      real(wp) :: t, ua, va, fa

      self%eng%Am = inA;  self%eng%Bm = inB;  self%eng%Cm = inC
      self%eng%time     = 0.0_wp
      self%eng%Un_prev  = 0.0_wp;  self%eng%Vn_prev = 0.0_wp
      call ve_step(self%eng, 0.0_wp, t, ua, va, fa)   ! σ = 0: pure relaxation
      outA = self%eng%Am;  outB = self%eng%Bm;  outC = self%eng%Cm
   end subroutine modal_apply_p

   subroutine modal_degree_destroy(self)
      type(modal_degree), intent(inout) :: self
      call ve_destroy(self%eng)
      self%j = -1;  self%nr = 0;  self%ne = 0;  self%dt = 0.0_wp
   end subroutine modal_degree_destroy

end module fe_modal
