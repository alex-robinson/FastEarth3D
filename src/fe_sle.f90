module fe_sle
   !! Sea-level equation: gravitationally self-consistent, mass-conserving
   !! redistribution of ocean water over a deforming solid Earth and geoid, with
   !! migrating coastlines (Kendall, Mitrovica & Milne 2005; Martinec et al.
   !! 2018 benchmark).
   !!
   !! Solved pseudo-spectrally: convolutions in spectral space, the ocean-
   !! function product `C*beta` as a pointwise multiply on the spatial grid (this
   !! is what avoids Gibbs ringing at coastlines). It is a Fredholm equation of
   !! the second kind, so it is iterated; outer iterations regenerate the
   !! paleotopography and ocean function (see doc/design.md, section D).
   !!
   !! STATUS: scaffold — ocean-function container + solver interface only.
   use fe_precision, only: wp
   use fe_constants, only: rho_ice, rho_water
   use fe_sht,       only: sht_grid
   implicit none
   private

   public :: sle_solver

   type :: sle_solver
      integer :: n_outer = 3   !! paleotopography / coastline iterations
      integer :: n_inner = 3   !! water-load fixed-point iterations
      ! TODO: ocean function C (water present), grounded mask beta, present-day
      ! topography, and the spatial-grid work arrays.
   contains
      procedure :: solve => sle_solve
   end type sle_solver

contains

   subroutine sle_solve(self, sht, ice_thickness, topo, rsl)
      !! Given an ice load and (paleo)topography, return relative sea level.
      !! STATUS: stub — interface matches the intended coupling fields.
      class(sle_solver), intent(inout) :: self
      type(sht_grid),    intent(in)     :: sht
      real(wp),          intent(in)     :: ice_thickness(:,:)  !! (nphi, nlat) [m]
      real(wp),          intent(in)     :: topo(:,:)           !! (nphi, nlat) [m]
      real(wp),          intent(out)    :: rsl(:,:)            !! (nphi, nlat) [m]
      rsl = 0.0_wp
      ! TODO:
      !   - partition ice into grounded (loads bed) vs floating (loads ocean)
      !     via the flotation criterion (rho_ice/rho_water)
      !   - outer loop: rebuild ocean function C and grounded mask beta
      !   - inner loop: solve the water-load increment by Love/Green convolution
      !   - add the uniform (eustatic) term pinned by ocean mass conservation
   end subroutine sle_solve

end module fe_sle
