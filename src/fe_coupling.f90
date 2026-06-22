module fe_coupling
   !! Top-level coupling API — the contract a host climate/ice model (CLIMBER-X)
   !! drives the solid-Earth model through. Mirrors the existing VILMA wrapper in
   !! CLIMBER-X (src/geo/vilma.F90): ice thickness goes in, relative sea level
   !! and bedrock elevation come out, on the model's own grid.
   !!
   !!   call se%init(...)                       ! once, with the reference state
   !!   call se%update(h_ice, rsl, z_bed)       ! every coupling step (e.g. 10 yr)
   !!   call se%finalize()
   !!
   !! The host is responsible for mapping between its grid and the model's
   !! Gauss-Legendre grid (CLIMBER-X uses conservative/bilinear SCRIP weights);
   !! all spherical-harmonic work stays inside this model.
   !!
   !! STATUS: scaffold — assembles the components and defines the API; the step
   !! is a no-op until the solver is implemented.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_sht,             only: sht_grid
   use fe_earth_structure, only: earth_model
   use fe_viscoelastic,    only: viscoelastic_state
   use fe_sle,             only: sle_solver
   use fe_rotation,        only: rotation_state
   implicit none
   private

   public :: solid_earth

   type :: solid_earth
      type(sht_grid)           :: sht       !! horizontal transform engine
      type(earth_model)        :: earth     !! radial (+ optional 3D) structure
      type(viscoelastic_state) :: visco     !! time-domain relaxation state
      type(sle_solver)         :: sle       !! sea-level equation
      type(rotation_state)     :: rotation  !! TPW feedback (off by default)
      real(wp), allocatable    :: z_bed_eq(:,:)  !! relaxed bedrock [m] (nphi,nlat)
      real(wp)                 :: dt_couple = 10.0_wp*kyr/1000.0_wp  !! [s], default 10 yr
   contains
      procedure :: init     => solid_earth_init
      procedure :: update   => solid_earth_update
      procedure :: finalize => solid_earth_finalize
   end type solid_earth

contains

   subroutine solid_earth_init(self, lmax, nlat, nphi)
      !! Configure the transform grid and the sub-solvers. STATUS: partial — the
      !! transform engine is real; structure/solver initialization is TODO.
      class(solid_earth), intent(inout) :: self
      integer,            intent(in)     :: lmax
      integer,            intent(in), optional :: nlat, nphi
      call self%sht%init(lmax, nlat=nlat, nphi=nphi)
      allocate(self%z_bed_eq(self%sht%nphi, self%sht%nlat), source=0.0_wp)
      ! TODO: load earth structure, assemble per-degree radial operators,
      ! initialize the viscoelastic state and the sea-level solver.
   end subroutine solid_earth_init

   subroutine solid_earth_update(self, h_ice, rsl, z_bed)
      !! Advance one coupling interval and return outputs. STATUS: stub.
      class(solid_earth), intent(inout) :: self
      real(wp), intent(in)  :: h_ice(:,:)   !! ice thickness [m] (nphi, nlat)
      real(wp), intent(out) :: rsl(:,:)     !! relative sea level [m]
      real(wp), intent(out) :: z_bed(:,:)   !! bedrock elevation [m]
      rsl   = 0.0_wp
      z_bed = self%z_bed_eq - rsl
      ! TODO: drive viscoelastic time steps over the interval under the
      ! ice+water load, solve the SLE for rsl, derive z_bed = z_bed_eq - rsl.
   end subroutine solid_earth_update

   subroutine solid_earth_finalize(self)
      class(solid_earth), intent(inout) :: self
      call self%sht%destroy()
      if (allocated(self%z_bed_eq)) deallocate(self%z_bed_eq)
   end subroutine solid_earth_finalize

end module fe_coupling
