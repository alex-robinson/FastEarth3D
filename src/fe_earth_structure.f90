module fe_earth_structure
   !! Reference Earth structure: radial layering plus an optional 3D (laterally
   !! varying) viscosity field.
   !!
   !! The solver is built 3D-ready from the start (project goal): the radial
   !! profile is the spherically symmetric *background*, and `visc_3d`, when
   !! allocated, carries log10-viscosity perturbations on the Gauss-Legendre
   !! spatial grid at each radial node — exactly how VILMA injects lateral
   !! heterogeneity (Albrecht et al. 2024). 1D runs simply leave `visc_3d`
   !! unallocated and the same code path reduces to the spherically symmetric
   !! case. See doc/design.md.
   !!
   !! STATUS: scaffold — types defined, loaders are stubs.
   use fe_precision, only: wp
   use fe_constants, only: r_earth
   implicit none
   private

   public :: earth_model

   type :: earth_model
      !! Radially discretized, incompressible Maxwell Earth.
      integer  :: nr = 0                      !! number of radial nodes
      real(wp), allocatable :: r(:)           !! radii of nodes [m] (nr)
      real(wp), allocatable :: rho(:)         !! density [kg m^-3] (nr)
      real(wp), allocatable :: mu(:)          !! elastic shear modulus [Pa] (nr)
      real(wp), allocatable :: eta(:)         !! background viscosity [Pa s] (nr)
      real(wp) :: r_lith = 0.0_wp             !! lithosphere base radius [m]
      ! 3D lateral viscosity: log10 perturbation on the (nlat*nphi, nr) grid.
      ! Unallocated => 1D run.
      real(wp), allocatable :: visc_3d(:,:)
   contains
      procedure :: init_1d => earth_model_init_1d
      procedure :: is_3d   => earth_model_is_3d
   end type earth_model

contains

   subroutine earth_model_init_1d(self, r, rho, mu, eta, r_lith)
      !! Populate a 1D (spherically symmetric) Earth from radial profiles.
      !! STATUS: stores the profiles; FE node insertion / interpolation TODO.
      class(earth_model), intent(inout) :: self
      real(wp), intent(in) :: r(:), rho(:), mu(:), eta(:)
      real(wp), intent(in) :: r_lith
      self%nr = size(r)
      self%r   = r
      self%rho = rho
      self%mu  = mu
      self%eta = eta
      self%r_lith = r_lith
   end subroutine earth_model_init_1d

   logical function earth_model_is_3d(self) result(yes)
      class(earth_model), intent(in) :: self
      yes = allocated(self%visc_3d)
   end function earth_model_is_3d

end module fe_earth_structure
