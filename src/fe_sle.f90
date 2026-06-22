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
   end type sle_result

   type :: sle_solver
      integer  :: n_outer = 3          !! paleotopography / coastline iterations
      integer  :: n_inner = 20         !! water-load fixed-point iterations
      real(wp) :: tol     = 1.0e-7_wp  !! inner convergence on max|ΔS| / max|S| [-]
   contains
      procedure :: solve => sle_solve
   end type sle_solver

contains

   subroutine sle_solve(self, sht, resp, d_ice, topo0, S, C, res)
      !! Solve for the relative-sea-level change S [m] driven by a grounded-ice
      !! thickness change d_ice [m], on a reference topography topo0 [m] (solid
      !! surface relative to the reference sea surface; ocean where < 0).
      class(sle_solver),        intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      class(response_operator), intent(inout) :: resp
      real(wp),                 intent(in)    :: d_ice(:,:)  !! (nphi,nlat) [m]
      real(wp),                 intent(in)    :: topo0(:,:)  !! (nphi,nlat) [m]
      real(wp),                 intent(out)   :: S(:,:)      !! (nphi,nlat) [m]
      real(wp),                 intent(out)   :: C(:,:)      !! (nphi,nlat) ocean fn
      type(sle_result),         intent(out)   :: res

      real(wp), allocatable :: load(:,:), u(:,:), N(:,:), Sraw(:,:), Snew(:,:)
      complex(wp), allocatable :: load_lm(:), u_lm(:), N_lm(:)
      real(wp) :: rho_ratio, ice_int, dphi, C_int, Cs_int, smax, dmax
      integer  :: io, ii, np, nl

      np = sht%nphi;  nl = sht%nlat
      allocate(load(np,nl), u(np,nl), N(np,nl), Sraw(np,nl), Snew(np,nl))
      allocate(load_lm(sht%nlm), u_lm(sht%nlm), N_lm(sht%nlm))

      rho_ratio = rho_ice/rho_water
      ! water-equivalent melt source ∝ −(ρ_i/ρ_w)∫ΔI dΩ (a² cancels in Δφ)
      ice_int = -rho_ratio * sht%surface_integral(d_ice)

      S = 0.0_wp
      res%n_inner_last = 0;  res%resid = 0.0_wp;  res%n_outer_done = 0

      ! Freeze the response's relaxation drift for this time step; for elastic /
      ! null responses this is a no-op.
      call resp%begin_step(sht)

      do io = 1, self%n_outer
         ! migrate the coastline: ocean where the current solid surface is below
         ! the (reference) sea surface, topo0 − S < 0.
         call ocean_function(topo0 - S, C)
         C_int = sht%surface_integral(C)
         if (C_int <= 0.0_wp) exit          ! no ocean: nothing to redistribute

         do ii = 1, self%n_inner
            ! total surface mass load = grounded ice + ocean water
            load = rho_ice*d_ice + rho_water*(C*S)
            call sht%analysis(load, load_lm)            ! analysis overwrites load
            call resp%apply(sht, load_lm, u_lm, N_lm)
            call sht%synthesis(u_lm, u)
            call sht%synthesis(N_lm, N)

            Sraw = N - u
            Cs_int = sht%surface_integral(C*Sraw)
            dphi   = (ice_int - Cs_int)/C_int           ! mass-conservation offset
            Snew   = C*(Sraw + dphi)

            dmax = maxval(abs(Snew - S))
            S    = Snew
            res%n_inner_last = ii;  res%resid = dmax
            smax = maxval(abs(S))
            if (dmax <= self%tol*max(smax, tiny(1.0_wp))) exit
         end do

         res%n_outer_done = io
      end do

      ! Commit the relaxation memory using the converged total load (no-op for
      ! elastic / null), advancing the response by one time step.
      load = rho_ice*d_ice + rho_water*(C*S)
      call sht%analysis(load, load_lm)
      call resp%commit_step(sht, load_lm)

      ! diagnostics
      Cs_int = sht%surface_integral(C*S)
      res%ocean_frac = C_int/(16.0_wp*atan(1.0_wp))      ! ∫C dΩ / 4π
      if (abs(ice_int) > 0.0_wp) then
         res%mass_resid = abs(Cs_int - ice_int)/abs(ice_int)
      else
         res%mass_resid = abs(Cs_int)
      end if
   end subroutine sle_solve

   subroutine ocean_function(topo, C)
      !! Migrating-coastline ocean function: C = 1 where the solid surface is
      !! below the sea surface (topo < 0), else 0. (Ice grounding is a later
      !! refinement; grounded ice belongs to the bedrock load, not the ocean.)
      real(wp), intent(in)  :: topo(:,:)
      real(wp), intent(out) :: C(:,:)
      where (topo < 0.0_wp)
         C = 1.0_wp
      elsewhere
         C = 0.0_wp
      end where
   end subroutine ocean_function

end module fe_sle
