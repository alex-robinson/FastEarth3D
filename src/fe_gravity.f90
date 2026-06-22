module fe_gravity
   !! Self-gravitation: the perturbed gravitational potential consistent with the
   !! deforming, redistributed mass (Poisson's equation), coupled into the radial
   !! momentum balance. Solved together with the FE radial operator degree by
   !! degree; kept as its own module so the rotational centrifugal-potential
   !! source term (fe_rotation) and the sea-level equation (fe_sle) can add their
   !! contributions through one interface.
   !!
   !! STATUS: scaffold — interface only.
   use fe_precision, only: wp
   use fe_earth_structure, only: earth_model
   implicit none
   private

   public :: potential_perturbation

contains

   subroutine potential_perturbation(earth, l, density_lm, phi_lm)
      !! Degree-l perturbed potential from a spectral density load. STATUS: stub.
      type(earth_model), intent(in)  :: earth
      integer,           intent(in)  :: l
      complex(wp),       intent(in)  :: density_lm(:)
      complex(wp),       intent(out) :: phi_lm(:)
      phi_lm = (0.0_wp, 0.0_wp)
      ! TODO: radial Green's-function / FE coupling for Poisson's equation.
   end subroutine potential_perturbation

end module fe_gravity
