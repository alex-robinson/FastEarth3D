module fe_coupling
   !! Top-level coupling API — the contract a host climate/ice model (CLIMBER-X)
   !! drives the solid-Earth model through. Mirrors the VILMA wrapper in
   !! CLIMBER-X (src/geo/vilma.F90): ice thickness goes in, relative sea level
   !! and bedrock elevation come out, on the model's own Gauss-Legendre grid.
   !!
   !!   call se%init(earth, sht, z_bed_eq, h_ice_ref, dt_couple, dt_step)
   !!   ...
   !!   call se%update(h_ice)       ! every coupling step; mutates se%rsl, se%z_bed
   !!   ...
   !!   call se%finalize()
   !!
   !! The host maps between its grid and the model's Gauss grid (CLIMBER-X uses
   !! conservative/bilinear SCRIP weights); all spherical-harmonic work stays
   !! inside this model. update() advances the viscoelastic state and the
   !! sea-level equation, then stores the results IN the derived type — the host
   !! reads se%rsl (relative sea level) and se%z_bed (bedrock) afterwards.
   !!
   !! Reference (equilibrium) state. The supplied (z_bed_eq, h_ice_ref) IS the
   !! relaxed state: the Maxwell memory starts at zero, so with h_ice = h_ice_ref
   !! we get d_ice = 0, rsl = 0, z_bed = z_bed_eq. Departures from the reference
   !! ice load drive the deformation and sea-level change incrementally. The
   !! reference topography z_bed_eq doubles as the SLE's reference topo0.
   use fe_precision,       only: wp
   use fe_sht,             only: sht_grid
   use fe_earth_structure, only: earth_model
   use fe_response,        only: ve_response
   use fe_sle,             only: sle_solver, sle_result
   use fe_rotation,        only: rotation_state
   implicit none
   private

   public :: solid_earth

   type :: solid_earth
      type(sht_grid), pointer  :: sht => null()  !! transform grid (borrowed, host-owned)
      type(earth_model)        :: earth     !! radial (+ optional 3D) structure
      type(ve_response)        :: resp      !! viscoelastic field driver (load → u, N)
      type(sle_solver)         :: sle       !! sea-level equation
      type(rotation_state)     :: rotation  !! TPW feedback (off by default)
      ! reference (equilibrium) state — set once at init
      real(wp), allocatable    :: z_bed_eq(:,:)   !! relaxed bedrock [m] (nphi,nlat)
      real(wp), allocatable    :: h_ice_ref(:,:)  !! reference grounded ice [m]
      ! current state — updated each coupling step
      real(wp), allocatable    :: h_ice(:,:)  !! current grounded-ice thickness [m]
      real(wp), allocatable    :: rsl(:,:)    !! relative sea level change [m] (full field)
      real(wp), allocatable    :: z_bed(:,:)  !! bedrock = z_bed_eq − rsl [m]
      real(wp), allocatable    :: C(:,:)      !! ocean function (1 ocean / 0 land)
      ! configuration
      real(wp) :: dt_couple = 0.0_wp   !! coupling interval [s]
      real(wp) :: dt_step   = 0.0_wp   !! internal Maxwell time step [s]
      integer  :: n_sub     = 0        !! steps per coupling interval = dt_couple/dt_step
      real(wp) :: time      = 0.0_wp   !! model time [s] (mirrors resp%time)
      ! diagnostics from the last update
      type(sle_result) :: last_res
      real(wp) :: worst_mass_resid = 0.0_wp  !! worst SLE mass residual over the interval
   contains
      procedure :: init     => solid_earth_init
      procedure :: update   => solid_earth_update
      procedure :: finalize => solid_earth_finalize
   end type solid_earth

contains

   subroutine solid_earth_init(self, earth, sht, z_bed_eq, h_ice_ref, &
                               dt_couple, dt_step)
      !! Wire the sub-solvers to a (host-owned) transform grid and set the
      !! reference state. The grid is borrowed by pointer — the host keeps it
      !! alive for the model's lifetime and is responsible for destroying it.
      class(solid_earth), intent(inout)        :: self
      type(earth_model),  intent(in)           :: earth
      type(sht_grid),     intent(in),   target :: sht
      real(wp),           intent(in)           :: z_bed_eq(:,:)   !! (nphi,nlat) [m]
      real(wp),           intent(in)           :: h_ice_ref(:,:)  !! (nphi,nlat) [m]
      real(wp),           intent(in)           :: dt_couple       !! coupling interval [s]
      real(wp),           intent(in)           :: dt_step         !! Maxwell step [s]
      integer :: np, nl

      call self%finalize()                       ! clean slate (safe on a fresh object)

      np = sht%nphi;  nl = sht%nlat
      if (size(z_bed_eq,1) /= np .or. size(z_bed_eq,2) /= nl .or. &
          size(h_ice_ref,1) /= np .or. size(h_ice_ref,2) /= nl) &
         error stop 'solid_earth_init: reference fields do not match the grid'

      self%n_sub = nint(dt_couple/dt_step)
      if (self%n_sub < 1 .or. &
          abs(self%n_sub*dt_step - dt_couple) > 1.0e-6_wp*dt_couple) &
         error stop 'solid_earth_init: dt_couple must be a positive multiple of dt_step'

      self%sht       => sht
      self%earth     = earth
      self%dt_couple = dt_couple
      self%dt_step   = dt_step
      self%time      = 0.0_wp

      ! reference (equilibrium) state; z_bed_eq doubles as the SLE topo0
      allocate(self%z_bed_eq,  source=z_bed_eq)
      allocate(self%h_ice_ref, source=h_ice_ref)

      ! current state — at the reference, nothing has moved yet
      allocate(self%h_ice(np,nl), source=h_ice_ref)
      allocate(self%rsl(np,nl),   source=0.0_wp)
      allocate(self%z_bed,        source=z_bed_eq)
      allocate(self%C(np,nl),     source=0.0_wp)

      ! viscoelastic driver: operators assembled + factored once, memory zeroed
      call self%resp%init(earth, sht, dt_step)
   end subroutine solid_earth_init

   subroutine solid_earth_update(self, h_ice)
      !! Advance one coupling interval under the ice thickness h_ice and store the
      !! results in the derived type (se%rsl, se%z_bed, se%C, se%h_ice). The ice
      !! load change relative to the reference, d_ice = h_ice − h_ice_ref, is held
      !! constant across the n_sub internal Maxwell steps; each step solves the
      !! SLE against the current relaxation state and advances the memory by Δt.
      class(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice(:,:)   !! grounded-ice thickness [m]
      real(wp), allocatable :: d_ice(:,:)
      integer :: k

      allocate(d_ice, source = h_ice - self%h_ice_ref)
      self%h_ice = h_ice
      self%worst_mass_resid = 0.0_wp

      do k = 1, self%n_sub
         call self%sle%solve(self%sht, self%resp, d_ice, h_ice, self%z_bed_eq, &
                             self%rsl, self%C, self%last_res)
         self%worst_mass_resid = max(self%worst_mass_resid, self%last_res%mass_resid)
      end do

      ! bedrock relative to the sea surface; rsl is the full RSL change field, so
      ! the bed subsides under grounded ice (land) as well as in the ocean.
      self%z_bed = self%z_bed_eq - self%rsl
      self%time  = self%resp%time
   end subroutine solid_earth_update

   subroutine solid_earth_finalize(self)
      class(solid_earth), intent(inout) :: self
      call self%resp%destroy()
      self%sht => null()                 ! borrowed grid — the host destroys it
      if (allocated(self%z_bed_eq))  deallocate(self%z_bed_eq)
      if (allocated(self%h_ice_ref)) deallocate(self%h_ice_ref)
      if (allocated(self%h_ice))     deallocate(self%h_ice)
      if (allocated(self%rsl))       deallocate(self%rsl)
      if (allocated(self%z_bed))     deallocate(self%z_bed)
      if (allocated(self%C))         deallocate(self%C)
      self%n_sub = 0;  self%time = 0.0_wp
   end subroutine solid_earth_finalize

end module fe_coupling
