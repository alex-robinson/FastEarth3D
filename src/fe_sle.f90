module fe_sle
   !! Sea-level equation: gravitationally self-consistent, mass-conserving
   !! redistribution of ocean water over a deforming solid Earth and geoid, with
   !! migrating coastlines (Kendall, Mitrovica & Milne 2005; Martinec et al.
   !! 2018 benchmark).
   !!
   !! Solved pseudo-spectrally: the load → (uplift, geoid) convolution is done in
   !! spectral space (via a response_operator, fe_response), while the ocean
   !! function multiply C·S is a pointwise product on the spatial Gauss grid —
   !! this is what avoids Gibbs ringing at coastlines.
   !!
   !! The equation is a Fredholm equation of the second kind, so it is iterated.
   !! For a fixed coastline the change in relative sea level over the ocean is
   !!
   !!     S = C·( N − u + Δφ ) ,
   !!
   !! where N (geoid) and u (uplift) are the response to the TOTAL surface load
   !! L = ρ_i ΔI + ρ_w (C·S), and the spatial constant Δφ (a uniform shift of the
   !! equipotential) is fixed each iteration by ocean-mass conservation,
   !!
   !!     ρ_w ∫ C·S dA = −ρ_i ∫ ΔI dA   (melt water volume) ,
   !!     Δφ = [ −(ρ_i/ρ_w) ∫ΔI dΩ − ∫ C(N−u) dΩ ] / ∫ C dΩ .
   !!
   !! Because Δφ is built to satisfy that balance, mass is conserved to machine
   !! precision at every iteration. The inner loop iterates S (the water load
   !! feeds back through the response); the outer loop rebuilds the ocean
   !! function C from the migrated topography topo0 − S (moving shorelines).
   use fe_precision, only: wp
   use fe_constants, only: rho_ice, rho_water
   use fe_sht,       only: sht_grid
   use fe_response,  only: response_operator
   implicit none
   private

   public :: sle_solver, sle_result

   type :: sle_result
      !! Diagnostics returned by a solve.
      integer  :: n_outer_done = 0      !! coastline iterations performed
      integer  :: n_inner_last = 0      !! inner iterations in the last outer pass
      real(wp) :: resid        = 0.0_wp !! last inner max|ΔS| [m]
      real(wp) :: mass_resid   = 0.0_wp !! relative ocean-mass-conservation error
      real(wp) :: ocean_frac   = 0.0_wp !! ∫C dΩ / 4π
      real(wp) :: esl          = 0.0_wp !! eustatic offset Δφ [m] (uniform sea-surface shift)
      real(wp), allocatable :: u(:,:)   !! converged solid uplift [m] (nphi,nlat)
      real(wp), allocatable :: N(:,:)   !! converged geoid rise [m] (nphi,nlat)
   end type sle_result

   type :: sle_solver
      integer  :: n_outer = 3          !! paleotopography / coastline iterations
      integer  :: n_inner = 20         !! water-load fixed-point iterations
      real(wp) :: tol     = 1.0e-7_wp  !! inner convergence on max|ΔS| / max|S| [-]
      !! Ocean geometry. .false. (default) = time-varying: the coastline migrates
      !! each outer pass from the deformed surface topo0 − rsl (Martinec 2018 §2.2,
      !! the SLE2 suite). .true. = fixed: the ocean function is held at the initial
      !! O⁽⁰⁾ = (topo0 < 0) throughout (§2.1, eq 1, the SLE1 suite) — no coastline
      !! migration, so a single outer pass converges the inner water load.
      logical  :: fixed_ocean = .false.
      !! Sloping-coast ("subgrid") ocean water load. .true. (default) uses the
      !! actual water-column change (Martinec 2018 §2.2, eqs 15-19): s = C·rsl −
      !! ζ⁽⁰⁾·(C − C⁽⁰⁾), which equals rsl over permanent ocean but rsl − ζ⁽⁰⁾ over
      !! newly-flooded cells, tapering to zero at the coast where the bed meets the
      !! sea surface — the mass-correct load on a migrating, sloping coastline.
      !! .false. is the simpler binary load ρ_w·C·rsl (the full sea-level change
      !! wherever C = 1, a sharp coastline), kept for comparison; it agrees with
      !! subgrid when the coastline does not move (deep basins, fixed_ocean).
      logical  :: subgrid = .true.
   contains
      procedure :: solve => sle_solve
   end type sle_solver

contains

   subroutine sle_solve(self, sht, resp, d_ice, ice, topo0, rsl, C, res)
      !! Solve for the relative-sea-level change rsl [m] driven by a grounded-ice
      !! thickness change d_ice [m], on a reference topography topo0 [m] (solid
      !! surface relative to the reference sea surface; ocean where < 0).
      !!
      !! rsl is the FULL relative-sea-level change field, N − u + Δφ (geoid rise
      !! minus solid uplift, plus the mass-conservation offset), defined over the
      !! WHOLE grid — not just the ocean. On land it is the bedrock-vs-sea-surface
      !! change that drives bedrock motion under grounded ice; over the ocean it
      !! is the sea-level change. The ocean-masked sea level is simply C·rsl, so
      !! it is not returned separately. The bedrock relative to the sea surface is
      !! topo0 − rsl everywhere.
      !!
      !! The absolute grounded-ice thickness ice [m] is passed alongside the
      !! change d_ice: d_ice drives the surface load, while ice enters the
      !! coastline test in ocean_function so that ice thick enough to ground on
      !! the bed is excluded from the ocean (it bears on the solid surface, it
      !! does not float — even where the bed has subsided below the sea surface).
      !! Both are needed because the load is incremental but flotation is an
      !! absolute condition.
      class(sle_solver),        intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      class(response_operator), intent(inout) :: resp
      real(wp),                 intent(in)    :: d_ice(:,:)  !! ice CHANGE [m] (load)
      real(wp),                 intent(in)    :: ice(:,:)    !! abs. ice [m] (flotation)
      real(wp),                 intent(in)    :: topo0(:,:)  !! (nphi,nlat) [m]
      real(wp),                 intent(out)   :: rsl(:,:)    !! full RSL change [m]
      real(wp),                 intent(out)   :: C(:,:)      !! (nphi,nlat) ocean fn
      type(sle_result),         intent(out)   :: res

      real(wp), allocatable :: load(:,:), u(:,:), N(:,:), Sraw(:,:), rsl_new(:,:)
      real(wp), allocatable :: C0(:,:), wcorr(:,:)
      complex(wp), allocatable :: load_lm(:), u_lm(:), N_lm(:)
      real(wp) :: rho_ratio, ice_int, dphi, C_int, Cs_int, zeta_int, smax, dmax
      integer  :: io, ii, np, nl

      np = sht%nphi;  nl = sht%nlat
      allocate(load(np,nl), u(np,nl), N(np,nl), Sraw(np,nl), rsl_new(np,nl))
      allocate(C0(np,nl), wcorr(np,nl))
      allocate(load_lm(sht%nlm), u_lm(sht%nlm), N_lm(sht%nlm))

      rho_ratio = rho_ice/rho_water
      rsl = 0.0_wp;  u = 0.0_wp;  N = 0.0_wp;  dphi = 0.0_wp;  ice_int = 0.0_wp
      zeta_int = 0.0_wp;  wcorr = 0.0_wp
      res%n_inner_last = 0;  res%resid = 0.0_wp;  res%n_outer_done = 0

      ! Initial ocean function O⁽⁰⁾ — the reference (t0) coastline against which the
      ! subgrid term measures newly flooded / emerged cells, and the held coastline
      ! in fixed-ocean mode. It is the flotation-aware ocean function of the
      ! reference state: bathymetry ζ⁽⁰⁾ = topo0 (rsl = 0 at t0) and the reference
      ! ice ice − d_ice (= 0 for an ice-free reference ⇒ O⁽⁰⁾ = (topo0 < 0); = ice
      ! when d_ice = 0 ⇒ grounded reference ice is excluded, so a no-change solve
      ! leaves C ≡ C⁽⁰⁾ and the subgrid term vanishes).
      call ocean_function(topo0, ice - d_ice, C0)

      ! Freeze the response's relaxation drift for this time step; for elastic /
      ! null responses this is a no-op.
      call resp%begin_step(sht)

      do io = 1, self%n_outer
         if (self%fixed_ocean) then
            ! Fixed ocean geometry (Martinec 2018 §2.1): hold the coastline at the
            ! reference O⁽⁰⁾ for all time (no migration). Computed once.
            if (io == 1) C = C0
         else
            ! migrate the coastline using the current (full-field) sea level: ocean
            ! where the deformed solid surface topo0 − rsl is below the sea surface
            ! AND the ice there floats rather than grounds (grounded ice keeps a
            ! subsided cell as land).
            call ocean_function(topo0 - rsl, ice, C)
         end if
         C_int = sht%surface_integral(C)
         if (C_int <= 0.0_wp) exit          ! no ocean: nothing to redistribute

         ! water-equivalent melt source ∝ −(ρ_i/ρ_w)∫ΔI dΩ over GROUNDED ice only:
         ! floating ice (C=1) is already in the ocean, so it does not change the
         ! ocean-water budget. Recomputed per coastline pass (grounded set shifts).
         ice_int = -rho_ratio * sht%surface_integral(d_ice*(1.0_wp - C))

         ! Subgrid sloping-coast correction (Martinec 2018 eq 17): the ocean water
         ! column change is C·rsl − ζ⁽⁰⁾·(C − C⁽⁰⁾), not C·rsl. The −ζ⁽⁰⁾(C−C⁽⁰⁾)
         ! piece accounts for the bed elevation of cells that crossed the coastline
         ! since t0 (a newly flooded cell fills from its bed, ζ⁽⁰⁾, not from the
         ! reference sea surface), so the load tapers to zero at the moving coast.
         ! Constant over the inner loop (depends only on the coastline C, C⁽⁰⁾). It
         ! also enters the mass balance via ζ̄⁽⁰⁾ = ∫ζ⁽⁰⁾(C−C⁽⁰⁾) (eqs 19-20).
         if (self%subgrid) then
            wcorr    = -rho_water * topo0 * (C - C0)
            zeta_int = sht%surface_integral(topo0*(C - C0))
         else
            wcorr = 0.0_wp;  zeta_int = 0.0_wp
         end if

         do ii = 1, self%n_inner
            ! total surface mass load = GROUNDED ice + ocean water. Ice over ocean
            ! cells (C=1: open ocean or floating ice) does not press its full weight
            ! on the bed -- it is borne by buoyancy and carried by the ocean term
            ! ρ_w·C·rsl. The (1−C) mask keeps the ice load only where it grounds
            ! (C=0). Without it, ice overhanging a deep basin over-loads the bed.
            ! wcorr is the subgrid sloping-coast term (zero unless self%subgrid).
            load = rho_ice*d_ice*(1.0_wp - C) + rho_water*(C*rsl) + wcorr
            call sht%analysis(load, load_lm)            ! analysis overwrites load
            call resp%apply(sht, load_lm, u_lm, N_lm)
            call sht%synthesis(u_lm, u)
            call sht%synthesis(N_lm, N)

            Sraw = N - u
            Cs_int = sht%surface_integral(C*Sraw)
            dphi   = (ice_int - Cs_int + zeta_int)/C_int ! mass-conservation offset
            rsl_new = Sraw + dphi                        ! full field, everywhere

            dmax = maxval(abs(C*(rsl_new - rsl)))        ! converge on the ocean part
            rsl  = rsl_new
            res%n_inner_last = ii;  res%resid = dmax
            smax = maxval(abs(C*rsl))
            if (dmax <= self%tol*max(smax, tiny(1.0_wp))) exit
         end do

         res%n_outer_done = io
         if (self%fixed_ocean) exit         ! C is fixed: one coastline pass converges
      end do

      ! Commit the relaxation memory using the converged total load (no-op for
      ! elastic / null), advancing the response by one time step. Same grounded-
      ! ice masking + subgrid sloping-coast term (wcorr) as the inner load.
      load = rho_ice*d_ice*(1.0_wp - C) + rho_water*(C*rsl) + wcorr
      call sht%analysis(load, load_lm)
      call resp%commit_step(sht, load_lm)

      ! diagnostics. The conserved ocean-water volume is ∫s dΩ = ∫C·rsl − ζ̄⁽⁰⁾
      ! (the subgrid sloping-coast term; ζ̄⁽⁰⁾ = 0 in the binary case), which must
      ! balance the melt source ice_int.
      Cs_int = sht%surface_integral(C*rsl) - zeta_int
      res%ocean_frac = C_int/(16.0_wp*atan(1.0_wp))      ! ∫C dΩ / 4π
      if (abs(ice_int) > 0.0_wp) then
         res%mass_resid = abs(Cs_int - ice_int)/abs(ice_int)
      else
         res%mass_resid = abs(Cs_int)
      end if
      res%u = u;  res%N = N;  res%esl = dphi             ! converged fields + offset
   end subroutine sle_solve

   subroutine ocean_function(topo, ice, C)
      !! Migrating-coastline ocean function with grounded-ice flotation. A cell
      !! is ocean (C = 1) only where BOTH
      !!   (a) the solid surface is below the sea surface, topo < 0, and
      !!   (b) the ice column is thin enough to float rather than ground:
      !!       ρ_i·I < −ρ_w·topo  (−topo > 0 is the water depth; the inequality
      !!       compares the ice draft to the column it would displace).
      !! Where ice grounds (ρ_i·I ≥ −ρ_w·topo) the cell is land (C = 0): the ice
      !! rests on and bears on the bed, so that column is not free ocean. With
      !! ice = 0 this reduces to the bare bathymetry test topo < 0.
      real(wp), intent(in)  :: topo(:,:)   !! solid surface vs. sea surface [m]
      real(wp), intent(in)  :: ice(:,:)    !! absolute grounded-ice thickness [m]
      real(wp), intent(out) :: C(:,:)
      where (topo < 0.0_wp .and. rho_ice*ice < -rho_water*topo)
         C = 1.0_wp
      elsewhere
         C = 0.0_wp
      end where
   end subroutine ocean_function

end module fe_sle
