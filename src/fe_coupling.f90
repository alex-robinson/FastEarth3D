module fe_coupling
   !! Top-level coupling API — the contract a host climate/ice model (CLIMBER-X)
   !! drives the solid-Earth model through. Mirrors the VILMA wrapper in
   !! CLIMBER-X (src/geo/vilma.F90): ice thickness goes in, relative sea level
   !! and bedrock elevation come out, on the model's own Gauss-Legendre grid.
   !!
   !!   call se%init(p, sht, z_bed_eq, h_ice_ref)   ! p = fe_param_class (one &fe3d group)
   !!   ...
   !!   call se%update(h_ice, dt)   ! advance time -> time+dt; mutates se%rsl, se%z_bed
   !!   ...
   !!   call se%finalize()
   !!
   !! The host maps between its grid and the model's Gauss grid (CLIMBER-X uses
   !! conservative/bilinear SCRIP weights); all spherical-harmonic work stays
   !! inside this model. update() advances the viscoelastic state and the
   !! sea-level equation across the interval [time, time+dt] — the ice load is
   !! ramped linearly between the previous and current h_ice and the internal Δt
   !! is chosen adaptively (fe_timestep) — then stores the results IN the derived
   !! type — the host reads se%rsl (relative sea level) and se%z_bed (bedrock).
   !!
   !! Reference (equilibrium) state. The supplied (z_bed_eq, h_ice_ref) IS the
   !! relaxed state: the Maxwell memory starts at zero, so with h_ice = h_ice_ref
   !! we get d_ice = 0, rsl = 0, z_bed = z_bed_eq. Departures from the reference
   !! ice load drive the deformation and sea-level change incrementally. The
   !! reference topography z_bed_eq doubles as the SLE's reference topo0.
   use fe_precision,       only: wp
   use fe_params,          only: fe_param_class
   use fe_sht,             only: sht_grid
   use fe_earth_structure, only: earth_model, build_earth
   use fe_viscoelastic,    only: scheme_from_name
   use fe_response,        only: ve_response
   use fe_sle,             only: sle_solver, sle_result
   use fe_timestep,        only: adaptive_stepper
   use fe_rotation,        only: rotation_state
   implicit none
   private

   public :: solid_earth

   type :: solid_earth
      type(sht_grid), pointer  :: sht => null()  !! transform grid (borrowed, host-owned)
      type(earth_model)        :: earth     !! radial (+ optional 3D) structure
      type(ve_response)        :: resp      !! viscoelastic field driver (load → u, N)
      type(sle_solver)         :: sle       !! sea-level equation
      type(adaptive_stepper)   :: stepper   !! adaptive-Δt controller (fe_timestep)
      type(rotation_state)     :: rotation  !! TPW feedback (off by default)
      ! reference (equilibrium) state — set once at init
      real(wp), allocatable    :: z_bed_eq(:,:)   !! relaxed bedrock [m] (nphi,nlat)
      real(wp), allocatable    :: h_ice_ref(:,:)  !! reference grounded ice [m]
      ! current state — updated each coupling step
      real(wp), allocatable    :: h_ice(:,:)  !! current grounded-ice thickness [m]
      real(wp), allocatable    :: rsl(:,:)    !! relative sea level change [m] (full field)
      real(wp), allocatable    :: z_bed(:,:)  !! bedrock = z_bed_eq − rsl [m]
      real(wp), allocatable    :: C(:,:)      !! ocean function (1 ocean / 0 land)
      ! configuration / clock
      real(wp) :: time      = 0.0_wp   !! model time [s] (mirrors resp%time)
      ! diagnostics from the last update
      real(wp) :: worst_mass_resid = 0.0_wp  !! worst SLE mass residual over the last interval
   contains
      procedure :: init     => solid_earth_init
      procedure :: update   => solid_earth_update
      procedure :: finalize => solid_earth_finalize
   end type solid_earth

contains

   subroutine solid_earth_init(self, p, sht, z_bed_eq, h_ice_ref)
      !! Build the model from the parameter record and wire it to a (host-owned)
      !! transform grid, setting the reference state. The earth structure is built
      !! from p (named built-in or custom layers); the SLE, adaptive-Δt and memory-
      !! scheme knobs are distributed from p to the sub-solvers. The grid is
      !! borrowed by pointer — the host keeps it alive for the model's lifetime and
      !! is responsible for destroying it.
      class(solid_earth),   intent(inout)        :: self
      type(fe_param_class), intent(in)           :: p
      type(sht_grid),       intent(in),   target :: sht
      real(wp),             intent(in)           :: z_bed_eq(:,:)   !! (nphi,nlat) [m]
      real(wp),             intent(in)           :: h_ice_ref(:,:)  !! (nphi,nlat) [m]
      integer  :: np, nl
      real(wp) :: dt0

      call self%finalize()                       ! clean slate (safe on a fresh object)

      np = sht%nphi;  nl = sht%nlat
      if (size(z_bed_eq,1) /= np .or. size(z_bed_eq,2) /= nl .or. &
          size(h_ice_ref,1) /= np .or. size(h_ice_ref,2) /= nl) &
         error stop 'solid_earth_init: reference fields do not match the grid'

      self%sht   => sht
      self%earth = build_earth(p)
      self%time  = 0.0_wp

      ! viscoelastic driver: operators assembled + factored once, memory zeroed.
      ! Δt enters only through Mk = (μ/η)Δt, which the adaptive stepper rescales
      ! per sub-step (resp%set_dt) — so the init Δt is just a nominal seed.
      dt0 = p%dt_init;  if (dt0 <= 0.0_wp) dt0 = p%dt_couple
      call self%resp%init(self%earth, sht, dt0)
      self%resp%scheme          = scheme_from_name(p%scheme)
      self%resp%max_couple_iter = p%max_couple_iter

      ! sea-level equation knobs. Warm-start the fixed point across steps: self%rsl
      ! persists (seeded to 0 below), and between adjacent steps the coastline
      ! barely moves, so the previous solution is a near-converged seed — sharply
      ! cutting the inner iteration count (the dominant per-step cost is the SLE's
      ! spherical-harmonic transforms, ~linear in that count).
      self%sle%n_outer      = p%sle_n_outer
      self%sle%n_inner      = p%sle_n_inner
      self%sle%tol          = p%sle_tol
      self%sle%max_mem_iter = p%sle_max_mem_iter
      self%sle%fixed_ocean  = p%sle_fixed_ocean
      self%sle%subgrid      = p%sle_subgrid
      self%sle%warm_start   = .true.

      ! adaptive-Δt controller (fe_timestep)
      self%stepper%rtol       = p%rtol
      self%stepper%atol       = p%atol
      self%stepper%safety     = p%safety
      self%stepper%grow_max   = p%grow_max
      self%stepper%shrink_min = p%shrink_min
      self%stepper%dt_min     = p%dt_min
      self%stepper%dt_max     = p%dt_max
      self%stepper%dt_try     = p%dt_init       ! 0 => first guess = whole interval

      ! rotational feedback
      self%rotation%enabled = p%rotation

      ! reference (equilibrium) state; z_bed_eq doubles as the SLE topo0
      allocate(self%z_bed_eq,  source=z_bed_eq)
      allocate(self%h_ice_ref, source=h_ice_ref)

      ! current state — at the reference, nothing has moved yet
      allocate(self%h_ice(np,nl), source=h_ice_ref)
      allocate(self%rsl(np,nl),   source=0.0_wp)
      allocate(self%z_bed,        source=z_bed_eq)
      allocate(self%C(np,nl),     source=0.0_wp)
   end subroutine solid_earth_init

   subroutine solid_earth_update(self, h_ice, dt)
      !! Advance the model from time to time+dt under the ice thickness h_ice and
      !! store the results in the derived type (se%rsl, se%z_bed, se%C, se%h_ice).
      !! The ice load is ramped linearly from the previous h_ice to the new one
      !! across the interval; the adaptive stepper (fe_timestep) chooses the
      !! internal Δt, solving the SLE against the current relaxation state and
      !! advancing the Maxwell memory at each sub-step.
      class(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice(:,:)   !! grounded-ice thickness [m]
      real(wp),           intent(in)    :: dt           !! interval to advance [s]
      real(wp) :: t0, t1

      t0 = self%time;  t1 = t0 + dt
      call self%stepper%advance(self%sht, self%resp, self%sle, self%z_bed_eq, &
                                self%h_ice, h_ice, self%h_ice_ref, t0, t1, &
                                self%rsl, self%C)
      self%h_ice            = h_ice
      self%worst_mass_resid = self%stepper%worst_mass_resid

      ! bedrock relative to the sea surface; rsl is the full RSL change field, so
      ! the bed subsides under grounded ice (land) as well as in the ocean.
      self%z_bed = self%z_bed_eq - self%rsl
      self%time  = t1
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
      self%time = 0.0_wp
   end subroutine solid_earth_finalize

end module fe_coupling
