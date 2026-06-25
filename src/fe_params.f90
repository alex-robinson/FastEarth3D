module fe_params
   !! Single configuration record for the whole solid-Earth model, loaded from one
   !! namelist group `&fe3d` (yelmo convention: a flat parameter type filled by
   !! nml_read, see fesm-utils/utils/src/nml.f90). Every runtime knob lives here;
   !! the high-level system init (solid_earth%init) consumes the whole record and
   !! distributes the values to the sub-solvers, while the specific component inits
   !! keep their direct-argument signatures.
   !!
   !! The standalone driver (program fastearth) reads `&fe3d` from a file and runs
   !! a forced simulation; a host model (CLIMBER-X) can instead fill the record in
   !! memory. The IO fields (file_*, name_*, time_*) are used only by the driver.
   use fe_precision, only: wp
   use fe_constants, only: kyr, sec_per_year
   use nml
   implicit none
   private

   public :: fe_param_class, fe_par_load, fe_par_print
   public :: MAX_LAYER

   integer, parameter :: MAX_LAYER = 16   !! cap on custom earth-structure layers

   type :: fe_param_class
      ! --- spectral grid (fe_sht) ------------------------------------------------
      integer  :: lmax  = 0       !! maximum spherical-harmonic degree (required)
      integer  :: nlat  = 0       !! Gauss latitudes  (0 => SHTns default = lmax+2)
      integer  :: nphi  = 0       !! longitudes       (0 => SHTns default)
      integer  :: mmax  = -1      !! maximum order    (<0 => = lmax)
      integer  :: mres  = 1       !! order stride
      real(wp) :: eps_polar = -1.0_wp   !! polar-optimization threshold (<0 => library default)

      ! --- earth structure (fe_earth_structure) ---------------------------------
      character(len=64) :: earth = "M3-L70-V01"
         !! named built-in model, or "custom" to build from the layer arrays below
      integer  :: n_layer = 0                  !! # custom layers (surface-first)
      real(wp) :: r_earth = 6371.0e3_wp        !! custom: surface radius [m]
      real(wp) :: r_core  = 3480.0e3_wp        !! custom: core-mantle boundary radius [m]
      real(wp) :: r_bot(MAX_LAYER) = 0.0_wp    !! custom: layer inner radii [m]
      real(wp) :: r_top(MAX_LAYER) = 0.0_wp    !! custom: layer outer radii [m]
      real(wp) :: rho(MAX_LAYER)   = 0.0_wp    !! custom: layer densities [kg m^-3]
      real(wp) :: mu(MAX_LAYER)    = 0.0_wp    !! custom: layer shear moduli [Pa]
      real(wp) :: eta(MAX_LAYER)   = 0.0_wp    !! custom: layer viscosities [Pa s]
      integer  :: rheology(MAX_LAYER) = 1      !! custom: 0=elastic 1=Maxwell 2=fluid

      ! --- viscoelastic memory scheme (fe_response / fe_viscoelastic) ------------
      character(len=8) :: scheme = "trap"      !! fe | etd1 | trap | be
      integer  :: max_couple_iter = 20         !! SLE<->memory co-convergence cap (implicit schemes)

      ! --- sea-level equation (fe_sle) ------------------------------------------
      integer  :: sle_n_outer      = 3
      integer  :: sle_n_inner      = 20
      real(wp) :: sle_tol          = 1.0e-7_wp
      integer  :: sle_max_mem_iter = 20
      logical  :: sle_fixed_ocean  = .false.
      logical  :: sle_subgrid      = .true.

      ! --- adaptive time stepping (fe_timestep) ----------------------------------
      ! Time fields are SI [s] in the record; the nml supplies them in YEARS and
      ! fe_par_load converts on read (so the in-memory record is uniformly SI).
      real(wp) :: dt_couple = kyr            !! default coupling interval [s] (driver cadence)
      real(wp) :: dt_init   = 0.0_wp         !! first trial Δt [s] (0 => try the whole interval)
      real(wp) :: dt_min    = 0.0_wp         !! Δt floor [s] (0 => none)
      real(wp) :: dt_max    = huge(1.0_wp)   !! Δt ceiling [s]
      real(wp) :: rtol      = 1.0e-4_wp      !! relative local-error tolerance (memory ∞-norm)
      real(wp) :: atol      = 1.0e-3_wp      !! absolute local-error floor
      real(wp) :: safety    = 0.9_wp         !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp         !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp        !! min Δt shrink per step

      ! --- rotational feedback (fe_rotation) ------------------------------------
      logical  :: rotation = .false.         !! TPW feedback (off until validated)

      ! --- driver I/O (program fastearth only) ----------------------------------
      character(len=512) :: file_forcing = ""            !! ice-thickness forcing (lon,lat,time)
      character(len=64)  :: name_ice     = "h_ice"       !! ice variable in file_forcing
      character(len=64)  :: name_time    = "time"        !! time axis [years] in file_forcing
      character(len=512) :: file_ref     = ""            !! reference state (z_bed_eq, h_ice_ref)
      character(len=64)  :: name_zbed_eq = "z_bed_eq"    !! bed var (legacy 2D ref; or 3D bed in forcing)
      character(len=64)  :: name_hice_ref = "h_ice_ref"  !! ice var (legacy 2D ref)
      character(len=512) :: file_out     = "fastearth_out.nc"  !! step output
      real(wp) :: time_init = -huge(1.0_wp)  !! start time [years] (default: first forcing slice)
      real(wp) :: time_end  =  huge(1.0_wp)  !! end time   [years] (default: last  forcing slice)
      ! --- online lon-lat -> Gauss remap of the forcing (program fastearth) ------
      logical  :: remap_input = .true.   !! .true.: forcing is lon-lat, remap on the fly;
                                         !! .false.: forcing already on the Gauss grid (legacy)
      character(len=64) :: name_lon = "lon"  !! source longitude axis variable [deg]
      character(len=64) :: name_lat = "lat"  !! source latitude  axis variable [deg]
      ! --- reference / equilibration (program fastearth) ------------------------
      ! In remap mode the reference is taken from the forcing file at the first
      ! in-window slice: bed = name_zbed_eq, ice = name_ice (both 3D lon-lat-time).
      integer  :: i_eq     = 1           !! 0: declare the start slice as equilibrium (memory 0);
                                         !! 1: ice-free reference, spin up under the start load
                                         !! (paleotopo fixed point) for dt_equil, then transient.
      real(wp) :: dt_equil = 50.0_wp*kyr !! spin-up hold per equilibration pass [s] (i_eq=1)
   end type fe_param_class

contains

   subroutine fe_par_load(p, filename, defaults_file, group)
      !! Fill the whole parameter record from the `&fe3d` group of `filename`,
      !! overlaid on a complete `defaults_file` (yelmo convention): every parameter
      !! must exist in the defaults file, but the user `filename` may set only the
      !! subset it wants to override. If `defaults_file` is omitted, `filename` IS
      !! its own defaults — so it must then be complete. Override `group` to read a
      !! differently-named namelist.
      type(fe_param_class), intent(inout) :: p
      character(len=*),     intent(in)    :: filename
      character(len=*),     intent(in), optional :: defaults_file
      character(len=*),     intent(in), optional :: group
      character(len=64)  :: g
      character(len=512) :: df
      real(wp) :: dt_couple_yr, dt_init_yr, dt_min_yr, dt_max_yr, time_init_yr, time_end_yr
      real(wp) :: dt_equil_yr

      g  = "fe3d";      if (present(group))         g  = group
      df = filename;    if (present(defaults_file)) df = defaults_file
      call nml_set_verbose(.false.)             ! fe_par_print echoes a concise summary instead

      ! grid
      call nml_read(filename, g, "lmax",      p%lmax,      defaults_file=df)
      call nml_read(filename, g, "nlat",      p%nlat,      defaults_file=df)
      call nml_read(filename, g, "nphi",      p%nphi,      defaults_file=df)
      call nml_read(filename, g, "mmax",      p%mmax,      defaults_file=df)
      call nml_read(filename, g, "mres",      p%mres,      defaults_file=df)
      call nml_read(filename, g, "eps_polar", p%eps_polar, defaults_file=df)

      ! earth structure
      call nml_read(filename, g, "earth",     p%earth,     defaults_file=df)
      call nml_read(filename, g, "n_layer",   p%n_layer,   defaults_file=df)
      call nml_read(filename, g, "r_earth",   p%r_earth,   defaults_file=df)
      call nml_read(filename, g, "r_core",    p%r_core,    defaults_file=df)
      call nml_read(filename, g, "r_bot",     p%r_bot,     defaults_file=df)
      call nml_read(filename, g, "r_top",     p%r_top,     defaults_file=df)
      call nml_read(filename, g, "rho",       p%rho,       defaults_file=df)
      call nml_read(filename, g, "mu",        p%mu,        defaults_file=df)
      call nml_read(filename, g, "eta",       p%eta,       defaults_file=df)
      call nml_read(filename, g, "rheology",  p%rheology,  defaults_file=df)

      ! memory scheme
      call nml_read(filename, g, "scheme",          p%scheme,          defaults_file=df)
      call nml_read(filename, g, "max_couple_iter", p%max_couple_iter, defaults_file=df)

      ! sea-level equation
      call nml_read(filename, g, "sle_n_outer",      p%sle_n_outer,      defaults_file=df)
      call nml_read(filename, g, "sle_n_inner",      p%sle_n_inner,      defaults_file=df)
      call nml_read(filename, g, "sle_tol",          p%sle_tol,          defaults_file=df)
      call nml_read(filename, g, "sle_max_mem_iter", p%sle_max_mem_iter, defaults_file=df)
      call nml_read(filename, g, "sle_fixed_ocean",  p%sle_fixed_ocean,  defaults_file=df)
      call nml_read(filename, g, "sle_subgrid",      p%sle_subgrid,      defaults_file=df)

      ! adaptive time stepping. The Δt / time fields are given in YEARS in the nml
      ! and converted to SI seconds here (the record is uniformly SI internally).
      dt_couple_yr = p%dt_couple/sec_per_year
      dt_init_yr   = p%dt_init  /sec_per_year
      dt_min_yr    = p%dt_min   /sec_per_year
      dt_max_yr    = p%dt_max   /sec_per_year
      call nml_read(filename, g, "dt_couple",  dt_couple_yr, defaults_file=df)
      call nml_read(filename, g, "dt_init",    dt_init_yr,   defaults_file=df)
      call nml_read(filename, g, "dt_min",     dt_min_yr,    defaults_file=df)
      call nml_read(filename, g, "dt_max",     dt_max_yr,    defaults_file=df)
      p%dt_couple = dt_couple_yr*sec_per_year
      p%dt_init   = dt_init_yr  *sec_per_year
      p%dt_min    = dt_min_yr   *sec_per_year
      p%dt_max    = dt_max_yr   *sec_per_year
      call nml_read(filename, g, "rtol",       p%rtol,       defaults_file=df)
      call nml_read(filename, g, "atol",       p%atol,       defaults_file=df)
      call nml_read(filename, g, "safety",     p%safety,     defaults_file=df)
      call nml_read(filename, g, "grow_max",   p%grow_max,   defaults_file=df)
      call nml_read(filename, g, "shrink_min", p%shrink_min, defaults_file=df)

      ! rotation
      call nml_read(filename, g, "rotation",   p%rotation,   defaults_file=df)

      ! driver I/O
      call nml_read(filename, g, "file_forcing",  p%file_forcing,  defaults_file=df)
      call nml_read(filename, g, "name_ice",      p%name_ice,      defaults_file=df)
      call nml_read(filename, g, "name_time",     p%name_time,     defaults_file=df)
      call nml_read(filename, g, "file_ref",      p%file_ref,      defaults_file=df)
      call nml_read(filename, g, "name_zbed_eq",  p%name_zbed_eq,  defaults_file=df)
      call nml_read(filename, g, "name_hice_ref", p%name_hice_ref, defaults_file=df)
      call nml_read(filename, g, "file_out",      p%file_out,      defaults_file=df)
      time_init_yr = p%time_init/sec_per_year
      time_end_yr  = p%time_end /sec_per_year
      call nml_read(filename, g, "time_init",     time_init_yr,    defaults_file=df)
      call nml_read(filename, g, "time_end",      time_end_yr,     defaults_file=df)
      p%time_init = time_init_yr*sec_per_year
      p%time_end  = time_end_yr *sec_per_year

      ! remap + equilibration
      call nml_read(filename, g, "remap_input",   p%remap_input,   defaults_file=df)
      call nml_read(filename, g, "name_lon",      p%name_lon,      defaults_file=df)
      call nml_read(filename, g, "name_lat",      p%name_lat,      defaults_file=df)
      call nml_read(filename, g, "i_eq",          p%i_eq,          defaults_file=df)
      dt_equil_yr = p%dt_equil/sec_per_year
      call nml_read(filename, g, "dt_equil",      dt_equil_yr,     defaults_file=df)
      p%dt_equil = dt_equil_yr*sec_per_year
   end subroutine fe_par_load

   subroutine fe_par_print(p, unit)
      !! Echo the active configuration (to stdout, or `unit` if given).
      type(fe_param_class), intent(in) :: p
      integer, intent(in), optional :: unit
      integer :: u, k
      u = 6;  if (present(unit)) u = unit

      write(u,'(a)')          ' [fe3d] configuration'
      write(u,'(a,i0,a,i0,a,i0)') '   grid:   lmax=', p%lmax, '  nlat=', p%nlat, '  nphi=', p%nphi
      write(u,'(a,a)')        '   earth:  ', trim(p%earth)
      if (trim(p%earth) == "custom") then
         do k = 1, p%n_layer
            write(u,'(a,i0,a,es9.2,a,es9.2,a,f8.1,a,es9.2,a,es9.2,a,i0)') &
                 '     layer ', k, ': r=[', p%r_bot(k), ',', p%r_top(k), &
                 ']  rho=', p%rho(k), '  mu=', p%mu(k), '  eta=', p%eta(k), &
                 '  rheol=', p%rheology(k)
         end do
      end if
      write(u,'(a,a,a,i0)')   '   scheme: ', trim(p%scheme), '   max_couple_iter=', p%max_couple_iter
      write(u,'(a,i0,a,i0,a,es8.1,a,l1,a,l1)') &
           '   sle:    n_outer=', p%sle_n_outer, '  n_inner=', p%sle_n_inner, &
           '  tol=', p%sle_tol, '  fixed_ocean=', p%sle_fixed_ocean, '  subgrid=', p%sle_subgrid
      write(u,'(a,es9.2,a,es9.2,a,es8.1,a,es8.1)') &
           '   dt:     couple=', p%dt_couple, '  init=', p%dt_init, &
           '  rtol=', p%rtol, '  atol=', p%atol
      write(u,'(a,l1)')       '   rotation: ', p%rotation
      write(u,'(a,l1,a,i0,a,es9.2)') '   forcing: remap_input=', p%remap_input, &
           '  i_eq=', p%i_eq, '  dt_equil=', p%dt_equil
   end subroutine fe_par_print

end module fe_params
