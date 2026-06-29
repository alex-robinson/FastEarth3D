module fe_drive
   !! Standalone forced-run driver. Given a &fe3d configuration (fe_params), build
   !! the transform grid and the solid-Earth model, read an ice-thickness forcing
   !! time series (and a bedrock reference) from netCDF, optionally remap it onto the
   !! Gauss grid on the fly, set up the reference / equilibration state, and march the
   !! model across the forcing, writing the diagnostic surface fields each step.
   !!
   !! General forced-transient driver: fields are read per (file, var, time-index) so
   !! experiments can mix and match. Two input modes:
   !!   remap_input=.true.  (default): the forcing is on a regular lon-lat grid and is
   !!                       conservatively remapped onto the model Gauss grid per slice
   !!                       (fe_remap); raw datasets stay on disk, no preprocessing.
   !!   remap_input=.false.: the forcing is already on the Gauss grid (host remapped,
   !!                       or produced by the offline fastearth_remap tool).
   !!
   !! Reference state (i_eq), mirroring CLIMBER-X i_equilibrium. RSL is measured
   !! against the relaxed reference z_bed_eq, so i_eq=1 yields rsl ~0 at present day.
   !!   i_eq=0  start slice is the equilibrium (z_bed_eq=bed[k0], h_ice_eq=ice[k0]).
   !!   i_eq=1  (default) present-day reference: z_bed_eq / h_ice_eq from
   !!           z_bed_ref_file / h_ice_ref_file (e.g. RTopo), remapped online.
   !!   i_eq=2  equilibrium read from z_bed_eq_file / h_ice_eq_file.
   !!   i_eq=3  z_bed_eq = present-day reference bed + rsl from rsl_restart_file.
   !! Optional (non-default) ice-free paleotopo spin-up when dt_equil>0: find the
   !! ice-free relaxed bed by a fixed point so the EQUILIBRIUM under the start-slice
   !! ice reproduces z_bed_eq, holding the load dt_equil per pass; this supersedes
   !! the i_eq-selected bed and leaves the model spun up (bed AND viscous memory).
   use fe_precision, only: wp
   use fe_constants, only: sec_per_year
   use fe_params,    only: fe_param_class, fe_par_load, fe_par_print
   use fe_sht,       only: sht_grid, sht_grid_destroy, sht_grid_surface_integral, sht_grid_init
   use fe_coupling,  only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth, &
                           solid_earth_enable_visc_3d
   use fe_response,  only: RESP_MODAL
   use fe_remap,     only: ll2gauss_map, ll2gauss_init, ll2gauss_apply
   use fe_io,        only: fe_write_step, fe_restart_write, fe_restart_read
   use ncio,         only: nc_read, nc_size, nc_exists_var
   use iso_fortran_env, only: error_unit
   implicit none
   private

   public :: fastearth_run

   ! diagnostic surface fields written each output step (the prognostic memory is
   ! written separately via fe_restart_write when a restart is wanted)
   character(len=8), parameter :: OUT_VARS(5) = &
        [character(len=8) :: "h_ice", "rsl", "z_bed", "C_ocean", "bsl"]

   integer,  parameter :: MAX_EQ  = 12       !! cap on LGM-memory spin-up passes
   real(wp), parameter :: EQ_TOL  = 1.0_wp   !! spin-up convergence: mean|Δz_bed/pass| [m]

contains

   subroutine fastearth_run(cfg_file, defaults_file)
      !! Run a forced simulation defined by the &fe3d group of cfg_file (overlaid on
      !! defaults_file if given; otherwise cfg_file must be complete).
      character(len=*), intent(in)           :: cfg_file
      character(len=*), intent(in), optional :: defaults_file

      type(fe_param_class)   :: p
      type(sht_grid), target :: sht
      type(solid_earth)      :: se
      type(ll2gauss_map)     :: rmap
      real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:), ice_lgm(:,:)
      real(wp), allocatable  :: src(:,:), tyr(:)
      real(wp) :: t0, dt
      integer  :: nt, k, k0, k1, np, nl, nlon, nls
      logical  :: remap
      character(len=:), allocatable :: rundir
      integer(kind=8) :: pc0, pc1, prate          ! PROFILE: per-step phase timers
      real(wp) :: t_read = 0.0_wp, t_upd = 0.0_wp, t_wrt = 0.0_wp
      real(wp) :: t_dr, t_mm                      ! PROFILE: solid_earth_update sub-phases
      integer  :: nstep = 0

      ! --- configuration --------------------------------------------------------
      if (present(defaults_file)) then
         call fe_par_load(p, cfg_file, defaults_file=defaults_file)
      else
         call fe_par_load(p, cfg_file)
      end if
      call fe_par_print(p)

      ! --- validate inputs up front (fail loudly on a config typo) --------------
      ! ncio's nc_read stop's with exit status 0 on a missing variable, so a typo
      ! in a name_* / file_* config would otherwise abort silently mid-run; check
      ! every (file, var) the run will read before doing any work.
      call validate_inputs(p)

      ! --- transform grid + work arrays -----------------------------------------
      call build_grid(p, sht)
      np = sht%nphi;  nl = sht%nlat
      allocate(z_bed_eq(np,nl), h_ice_eq(np,nl), h_ice(np,nl), ice_lgm(np,nl))

      ! --- build the conservative lon-lat -> Gauss map (online remap mode) -------
      remap = p%remap_input
      call system_clock(pc0, prate)
      if (remap) then
         call build_remap(p, sht, rmap, nlon, nls)
         allocate(src(nlon, nls))
         write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' remap: lon-lat (', nlon, 'x', nls, &
              ') -> Gauss (', np, 'x', nl, ')'
      end if
      call system_clock(pc1)
      write(*,'(a,f8.2,a)') ' [PROFILE setup] build_remap =', real(pc1-pc0,wp)/prate, ' s'

      ! --- forcing time axis (years) + window -----------------------------------
      nt = nc_size(p%file_forcing, trim(p%name_time))
      allocate(tyr(nt));  call nc_read(p%file_forcing, trim(p%name_time), tyr)
      call select_window(tyr*sec_per_year, p%time_init, p%time_end, k0, k1)
      if (k1 <= k0) error stop 'fastearth_run: forcing window contains < 2 time slices'
      write(*,'(a,i0,a,f0.1,a,f0.1,a)') ' fastearth: ', k1-k0+1, ' slices, t = ', &
           tyr(k0), ' -> ', tyr(k1), ' yr'

      ! --- start-slice ice + reference state (z_bed_eq, h_ice_eq) per i_eq ------
      call read_ice(p, rmap, sht, remap, k0, src, ice_lgm)      ! start-slice (LGM) ice
      call setup_reference(p, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_eq)

      ! --- initialise; restart resume and/or optional LGM-memory spin-up --------
      ! Restarts are written alongside the run output: <rundir>/spinup after the spin-up
      ! and <rundir>/final at the end, where rundir is the directory of file_out.
      rundir = dir_of(p%file_out)
      call system_clock(pc0, prate)
      if (len_trim(p%restart_in_file) > 0) then
         ! resume from a saved full state (memory + clock). init at the reference, then
         ! restore; a lower-resolution restart is interpolated up to the model grid.
         se%par = p; call solid_earth_init(se, sht, z_bed_eq, h_ice_eq)
         call fe_restart_read(se, trim(p%restart_in_file))
         if (p%dt_equil > 0.0_wp) then            ! optional further equilibration
            call equilibrate(p, sht, se, z_bed_eq, h_ice_eq, ice_lgm, from_restart=.true.)
            call fe_restart_write(se, se%time, folder=trim(rundir)//"/spinup")
         else
            ! seed the entering (LGM) ice and re-solve the diagnostics from the restored
            ! memory (also fills rsl/z_bed/C after a cross-resolution restart).
            call solid_earth_update(se, ice_lgm, 0.0_wp)
         end if
      else if (p%dt_equil > 0.0_wp) then
         ! spin up se to isostatic equilibrium under the start-slice (LGM) ice while
         ! HOLDING the present-day reference (z_bed_eq, h_ice_eq) as the datum, so the
         ! transient measures rsl/bsl against today (-> 0 as ice -> h_ice_eq).
         call equilibrate(p, sht, se, z_bed_eq, h_ice_eq, ice_lgm, spinup_1d=p%spinup_1d)
         call fe_restart_write(se, se%time, folder=trim(rundir)//"/spinup")
      else
         se%par = p; call solid_earth_init(se, sht, z_bed_eq, h_ice_eq)
         call solid_earth_update(se, ice_lgm, 0.0_wp)                        ! seed entering ice, no integration
      end if
      call system_clock(pc1)
      write(*,'(a,f8.2,a)') ' [PROFILE setup] solid_earth_init+visc3d+seed =', real(pc1-pc0,wp)/prate, ' s'

      ! Lateral-viscosity diagnostic. RESP_VE: how many radial elements are genuinely
      ! 3-D (pay the tensor-SH advance) vs collapse to the cheap 1-D path. RESP_MODAL:
      ! how many within-degree mode ranks carry a lateral anomaly (the K scalar-SHT path).
      if (p%l_visc_3d) then
         if (se%resp%kind == RESP_MODAL) then
            write(*,'(a,i0,a,i0,a)') ' modal lateral viscosity: ', se%resp%nrank3d, &
                 ' of ', se%resp%maxmode, ' mode ranks carry a lateral anomaly'
         else
            write(*,'(a,i0,a,i0,a)') ' visc3d split: ', se%resp%ne3d, &
                 ' of ', se%resp%ne, ' radial elements laterally 3-D (rest advance as 1-D)'
         end if
      end if

      ! --- march the transient --------------------------------------------------
      t0 = tyr(k0)*sec_per_year
      se%time = t0;  se%resp%time = t0
      call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.true.)

      se%resp%t_drift = 0.0_wp;  se%resp%t_mem = 0.0_wp   ! PROFILE: time the transient only
      se%resp%n_drift = 0;       se%resp%n_mem = 0
      do k = k0, k1-1
         dt = (tyr(k+1) - tyr(k))*sec_per_year
         call system_clock(pc0, prate)
         call read_ice(p, rmap, sht, remap, k+1, src, h_ice)
         call system_clock(pc1);  t_read = t_read + real(pc1-pc0,wp)/prate
         call system_clock(pc0)
         call solid_earth_update(se, h_ice, dt)
         call system_clock(pc1);  t_upd = t_upd + real(pc1-pc0,wp)/prate
         call system_clock(pc0)
         call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.false.)
         call system_clock(pc1);  t_wrt = t_wrt + real(pc1-pc0,wp)/prate
         nstep = nstep + 1
         write(*,'(a,f12.2,a,es10.2,a,es9.2)') '   t=', se%time/sec_per_year, &
              ' yr   max|rsl|=', maxval(abs(se%rsl)), '   mass_resid=', se%worst_mass_resid
      end do

      if (nstep > 0) write(*,'(a,i0,a,/,3(a,f8.1,a,f5.1,a,/))') &
         ' [PROFILE] per coupling step (mean over ', nstep, ' steps):', &
         '   read_ice (remap+IO) =', 1.0e3_wp*t_read/nstep, ' ms (', &
            100.0_wp*t_read/(t_read+t_upd+t_wrt), ' %)', &
         '   solid_earth_update  =', 1.0e3_wp*t_upd /nstep, ' ms (', &
            100.0_wp*t_upd /(t_read+t_upd+t_wrt), ' %)', &
         '   fe_write_step (out) =', 1.0e3_wp*t_wrt /nstep, ' ms (', &
            100.0_wp*t_wrt /(t_read+t_upd+t_wrt), ' %)'
      ! solid_earth_update internal split (wall-clock, accumulated in the response).
      ! drift = per-degree band LU; memory advance = the Maxwell update (3-D dyadic SHT
      ! round-trip when laterally 3-D); rest = SLE iteration + load/geoid SHTs.
      if (nstep > 0) then
         t_dr = se%resp%t_drift;  t_mm = se%resp%t_mem
         write(*,'(a,/,3(a,f8.1,a,f5.1,a,/))') &
            ' [PROFILE] solid_earth_update breakdown (per step, wall-clock):', &
            '   drift solve (band LU) =', 1.0e3_wp*t_dr/nstep, ' ms (', &
               100.0_wp*t_dr/max(t_upd,tiny(1.0_wp)), ' % of update)', &
            '   memory advance        =', 1.0e3_wp*t_mm/nstep, ' ms (', &
               100.0_wp*t_mm/max(t_upd,tiny(1.0_wp)), ' % of update)', &
            '   SLE + coupling (rest) =', 1.0e3_wp*(t_upd-t_dr-t_mm)/nstep, ' ms (', &
               100.0_wp*(t_upd-t_dr-t_mm)/max(t_upd,tiny(1.0_wp)), ' % of update)'
      end if
      if (nstep > 0) write(*,'(a,f7.1,a,f7.1,a)') &
         '   sub-steps/interval: n_accept=', real(se%stepper%n_accept,wp)/nstep, &
         '  n_solve=', real(se%stepper%n_solve,wp)/nstep, '  (per coupling step)'
      write(*,'(a,a)') ' fastearth: wrote ', trim(p%file_out)
      call fe_restart_write(se, se%time, folder=trim(rundir)//"/final")
      write(*,'(a,a)') ' fastearth: wrote restart ', trim(rundir)//'/final/fe_restart.nc'
      call solid_earth_finalize(se);  call sht_grid_destroy(sht)
   end subroutine fastearth_run

   ! --- equilibration ----------------------------------------------------------

   subroutine equilibrate(p, sht, se, z_bed_eq, h_ice_eq, ice_lgm, spinup_1d, from_restart)
      !! LGM-memory spin-up (i_eq=1): bring se to isostatic equilibrium under the
      !! start-slice ice ice_lgm while HOLDING the present-day reference (z_bed_eq =
      !! observed bed, h_ice_eq = observed ice) as the fixed datum. The model is
      !! initialised at the reference (memory 0, rsl 0), the entering ice is set to
      !! ice_lgm, then ice_lgm is held for dt_equil per pass until the bed stops moving
      !! between passes (the slow Maxwell modes have relaxed). On return se carries the
      !! relaxed LGM deflection rsl = equilibrium response to the load anomaly
      !! (ice_lgm - h_ice_eq) and its consistent memory state; the reference datum is
      !! UNCHANGED, so the transient that follows measures rsl/bsl against today
      !! (rsl, bsl -> 0 as ice -> h_ice_eq).
      !!
      !! spinup_1d:    run the relaxation with 1-D (radial) viscosity, then enable the
      !!               3-D field for the transient (cheap seed; ignored if not l_visc_3d).
      !! from_restart: se is already initialised + state-restored, so skip the init and
      !!               continue relaxing from the restored memory (further-equilibration).
      type(fe_param_class), intent(in)            :: p
      type(sht_grid),       intent(in),   target  :: sht
      type(solid_earth),    intent(inout)         :: se
      real(wp),             intent(in)            :: z_bed_eq(:,:)   !! PD reference bed (datum, held)
      real(wp),             intent(in)            :: h_ice_eq(:,:)  !! PD reference ice (datum, held)
      real(wp),             intent(in)            :: ice_lgm(:,:)
      logical, optional,    intent(in)            :: spinup_1d, from_restart
      real(wp), allocatable :: z_prev(:,:), resid(:,:)
      real(wp) :: rmean, rmax, fourpi
      integer  :: it, np, nl
      logical  :: l1d, resume

      l1d    = .false.;  if (present(spinup_1d))    l1d    = spinup_1d
      resume = .false.;  if (present(from_restart)) resume = from_restart
      np = sht%nphi;  nl = sht%nlat
      fourpi = 16.0_wp*atan(1.0_wp)
      allocate(z_prev(np,nl), resid(np,nl))

      write(*,'(a,a,a,es9.2,a)') ' equilibration (i_eq=1): spin up LGM memory vs PD reference', &
           merge(' [1-D]', '      ', l1d .and. p%l_visc_3d), ', dt_equil=', p%dt_equil/sec_per_year, ' yr/pass'
      if (.not. resume) then
         se%par = p
         call solid_earth_init(se, sht, z_bed_eq, h_ice_eq, defer_visc_3d=l1d)  ! reference, memory 0
      end if
      call solid_earth_update(se, ice_lgm, 0.0_wp)                   ! set entering ice = start (LGM) ice
      z_prev = se%z_bed
      do it = 1, MAX_EQ
         call solid_earth_update(se, ice_lgm, p%dt_equil)            ! hold LGM ice -> relax further
         resid = se%z_bed - z_prev                                   ! bed motion since the previous pass
         rmean = sht_grid_surface_integral(sht, abs(resid))/fourpi
         rmax  = maxval(abs(resid))
         write(*,'(a,i2,a,f10.4,a,f10.2,a)') '   pass ', it, ':  d|z_bed|(vs prev pass)=', &
              rmean, ' m   max=', rmax, ' m'
         if (rmean < EQ_TOL) exit
         z_prev = se%z_bed
      end do

      ! spinup_1d: the relaxation ran 1-D; enable the 3-D field for the transient,
      ! preserving the spun-up memory as the seed.
      if (l1d .and. p%l_visc_3d) call solid_earth_enable_visc_3d(se, sht)
   end subroutine equilibrate

   pure function dir_of(path) result(d)
      !! Directory part of a path (everything before the last "/"); "." if none.
      character(len=*), intent(in)  :: path
      character(len=:), allocatable :: d
      integer :: i
      i = index(trim(path), "/", back=.true.)
      if (i > 0) then;  d = path(1:i-1);  else;  d = "."; end if
   end function dir_of

   ! --- input helpers ----------------------------------------------------------

   subroutine validate_inputs(p)
      !! Check every (file, variable) the run will actually read exists, before any
      !! work begins. Mirrors the conditional read paths (remap_input, i_eq) so a
      !! mistyped name_* / file_* config fails loudly with nonzero status instead of
      !! aborting silently mid-run (ncio's nc_read stop's with exit status 0 on a
      !! missing variable). 'what' names the offending config key in the message.
      type(fe_param_class), intent(in) :: p
      logical :: remap
      remap = p%remap_input

      ! forcing file: time axis + ice always; lon/lat only when remapping online.
      call require_var(p%file_forcing, p%name_time, 'file_forcing/name_time')
      call require_var(p%file_forcing, p%name_ice,  'file_forcing/name_ice')
      if (remap) then
         call require_var(p%file_forcing, p%name_lon, 'file_forcing/name_lon')
         call require_var(p%file_forcing, p%name_lat, 'file_forcing/name_lat')
      end if

      ! relaxed reference state (z_bed_eq, h_ice_eq), per i_eq -- see read_bed /
      ! setup_reference for which (file, var) each case reads.
      select case (p%i_eq)
      case (0)
         if (remap) then
            call require_var(p%file_forcing, p%name_zbed_eq, 'file_forcing/name_zbed_eq')
         else
            call require_var(p%file_ref, p%name_zbed_eq, 'file_ref/name_zbed_eq')
         end if
      case (1)
         call require_var(p%z_bed_ref_file, p%name_z_bed_ref, 'z_bed_ref_file/name_z_bed_ref')
         call require_var(p%h_ice_ref_file, p%name_h_ice_ref, 'h_ice_ref_file/name_h_ice_ref')
      case (2)
         call require_var(p%z_bed_eq_file, p%name_z_bed_ref, 'z_bed_eq_file/name_z_bed_ref')
         call require_var(p%h_ice_eq_file, p%name_h_ice_ref, 'h_ice_eq_file/name_h_ice_ref')
      case (3)
         call require_var(p%z_bed_ref_file,   p%name_z_bed_ref, 'z_bed_ref_file/name_z_bed_ref')
         call require_var(p%rsl_restart_file, p%name_rsl,       'rsl_restart_file/name_rsl')
         call require_var(p%h_ice_ref_file,   p%name_h_ice_ref, 'h_ice_ref_file/name_h_ice_ref')
      case default
         error stop 'fastearth_run: i_eq must be 0, 1, 2, or 3'
      end select
   end subroutine validate_inputs

   subroutine require_var(file, var, what)
      !! error stop (nonzero, flushed) if 'file' is unset/missing or lacks netCDF
      !! variable 'var'. 'what' is the config key(s) involved, named in the message.
      character(len=*), intent(in) :: file, var, what
      logical :: ok
      if (len_trim(file) == 0) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": file not set for the chosen i_eq'
         flush(error_unit)
         error stop 'fastearth_run: required input file not set (see message above)'
      end if
      inquire(file=trim(file), exist=ok)
      if (.not. ok) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": file not found'
         write(error_unit,'(a)') '   file: '//trim(file)
         flush(error_unit)
         error stop 'fastearth_run: input file not found (see message above)'
      end if
      if (.not. nc_exists_var(trim(file), trim(var))) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": variable "'//trim(var)//'" not found'
         write(error_unit,'(a)') '   file: '//trim(file)
         flush(error_unit)
         error stop 'fastearth_run: required netCDF variable not found (see message above)'
      end if
   end subroutine require_var

   subroutine build_remap(p, sht, rmap, nlon, nls)
      !! Read the source lon/lat axes from the forcing file and build the conservative
      !! lon-lat -> Gauss map. Returns the source dimensions for the work buffer.
      type(fe_param_class), intent(in)    :: p
      type(sht_grid),       intent(in)    :: sht
      type(ll2gauss_map),   intent(out)   :: rmap
      integer,              intent(out)   :: nlon, nls
      real(wp), allocatable :: lon_s(:), lat_s(:)
      nlon = nc_size(p%file_forcing, trim(p%name_lon))
      nls  = nc_size(p%file_forcing, trim(p%name_lat))
      allocate(lon_s(nlon), lat_s(nls))
      call nc_read(p%file_forcing, trim(p%name_lon), lon_s)
      call nc_read(p%file_forcing, trim(p%name_lat), lat_s)
      call ll2gauss_init(rmap, sht, lon_s, lat_s)
   end subroutine build_remap

   subroutine read_ice(p, rmap, sht, remap, k, src, h_ice)
      !! Read ice slice k onto the Gauss grid. remap: read lon-lat then conservatively
      !! remap (mass-conserving); else read the Gauss-grid slice directly.
      type(fe_param_class), intent(in)    :: p
      type(ll2gauss_map),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(out)   :: h_ice(:,:)
      if (remap) then
         call nc_read(p%file_forcing, trim(p%name_ice), src, &
                      start=[1,1,k], count=[size(src,1), size(src,2), 1])
         call ll2gauss_apply(rmap, sht, src, h_ice, conserve_mass=.true.)
      else
         call nc_read(p%file_forcing, trim(p%name_ice), h_ice, &
                      start=[1,1,k], count=[size(h_ice,1), size(h_ice,2), 1])
      end if
   end subroutine read_ice

   subroutine read_bed(p, rmap, sht, remap, k, src, z_bed)
      !! Read the reference bedrock. remap: slice k of name_zbed_eq from the forcing
      !! file (lon-lat-time), remapped (no mass rescale -- bed is geometry, not mass).
      !! else: the 2D name_zbed_eq from file_ref (legacy Gauss-grid reference file).
      type(fe_param_class), intent(in)    :: p
      type(ll2gauss_map),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(out)   :: z_bed(:,:)
      if (remap) then
         call nc_read(p%file_forcing, trim(p%name_zbed_eq), src, &
                      start=[1,1,k], count=[size(src,1), size(src,2), 1])
         call ll2gauss_apply(rmap, sht, src, z_bed, conserve_mass=.false.)
      else
         if (len_trim(p%file_ref) == 0) error stop 'fastearth_run: file_ref not set (remap_input=.false.)'
         call nc_read(p%file_ref, trim(p%name_zbed_eq), z_bed)
      end if
   end subroutine read_bed

   subroutine setup_reference(p, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_eq)
      !! Fill the relaxed reference (z_bed_eq = SLE topo0, h_ice_eq) per p%i_eq,
      !! mirroring CLIMBER-X i_equilibrium. Reference files are lon-lat and remapped
      !! online with a per-file conservative map (their grid differs from the forcing).
      type(fe_param_class), intent(in)    :: p
      type(ll2gauss_map),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k0
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(in)    :: ice_lgm(:,:)
      real(wp),             intent(out)   :: z_bed_eq(:,:), h_ice_eq(:,:)
      real(wp), allocatable :: rsl_r(:,:)
      select case (p%i_eq)
      case (0)
         call read_bed(p, rmap, sht, remap, k0, src, z_bed_eq)   ! data bed at the start slice
         h_ice_eq = ice_lgm
      case (1)
         call read_ref2d(p, sht, remap, p%z_bed_ref_file, p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%h_ice_ref_file, p%name_h_ice_ref, .true.,  h_ice_eq)
      case (2)
         call read_ref2d(p, sht, remap, p%z_bed_eq_file, p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%h_ice_eq_file, p%name_h_ice_ref, .true.,  h_ice_eq)
      case (3)
         allocate(rsl_r(size(z_bed_eq,1), size(z_bed_eq,2)))
         call read_ref2d(p, sht, remap, p%z_bed_ref_file,   p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%rsl_restart_file, p%name_rsl,       .false., rsl_r)
         call read_ref2d(p, sht, remap, p%h_ice_ref_file,   p%name_h_ice_ref, .true.,  h_ice_eq)
         z_bed_eq = z_bed_eq + rsl_r
      case default
         error stop 'fastearth_run: i_eq must be 0, 1, 2, or 3'
      end select
   end subroutine setup_reference

   subroutine read_ref2d(p, sht, remap, file, varname, conserve, field)
      !! Read a 2D lon-lat reference field onto the Gauss grid. If the file is
      !! already on the Gauss grid (its lon/lat dims match nphi/nlat — e.g. a
      !! prebaked reference from fastearth_mkref) it is read directly, skipping the
      !! expensive conservative-map build. Otherwise a per-file conservative map is
      !! built from the file's own axes and applied (conserve=.true. for ice). In
      !! legacy (remap=.false.) mode the field is read directly.
      type(fe_param_class), intent(in)  :: p
      type(sht_grid),       intent(in)  :: sht
      logical,              intent(in)  :: remap, conserve
      character(len=*),     intent(in)  :: file, varname
      real(wp),             intent(out) :: field(:,:)
      type(ll2gauss_map)    :: m
      real(wp), allocatable :: lon_s(:), lat_s(:), buf(:,:)
      integer :: nlon, nls
      if (len_trim(file) == 0) &
         error stop 'fastearth_run: reference file not set for the chosen i_eq'
      if (remap) then
         nlon = nc_size(file, trim(p%name_lon));  nls = nc_size(file, trim(p%name_lat))
         if (nlon == sht%nphi .and. nls == sht%nlat) then
            call nc_read(file, trim(varname), field)   ! already on the Gauss grid
            return
         end if
         allocate(lon_s(nlon), lat_s(nls), buf(nlon, nls))
         call nc_read(file, trim(p%name_lon), lon_s)
         call nc_read(file, trim(p%name_lat), lat_s)
         call nc_read(file, trim(varname),    buf)
         call ll2gauss_init(m, sht, lon_s, lat_s)
         call ll2gauss_apply(m, sht, buf, field, conserve_mass=conserve)
      else
         call nc_read(file, trim(varname), field)
      end if
   end subroutine read_ref2d

   ! --- grid + window ----------------------------------------------------------

   subroutine build_grid(p, sht)
      !! Build the Gauss-Legendre transform grid from p. When nlat/nphi are unset
      !! (<=0) default to a de-aliased grid (nlat=2 lmax+2, nphi=4 lmax) sized for
      !! the SLE's quadratic ocean-function product.
      type(fe_param_class), intent(in)            :: p
      type(sht_grid),       intent(inout), target :: sht
      integer :: nlat, nphi
      nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
      nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
      if (p%mmax >= 0 .and. p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
      else if (p%mmax >= 0) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
      else if (p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
      else
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
      end if
   end subroutine build_grid

   subroutine select_window(tsec, t_init, t_end, k0, k1)
      !! First/last forcing indices inside [t_init, t_end] (all in seconds).
      real(wp), intent(in)  :: tsec(:), t_init, t_end
      integer,  intent(out) :: k0, k1
      integer :: k
      k0 = 0;  k1 = 0
      do k = 1, size(tsec)
         if (tsec(k) >= t_init .and. tsec(k) <= t_end) then
            if (k0 == 0) k0 = k
            k1 = k
         end if
      end do
      if (k0 == 0) error stop 'fastearth_run: no forcing times within [time_init, time_end]'
   end subroutine select_window

end module fe_drive
