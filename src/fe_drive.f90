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
   !!   i_eq=0  start slice is the equilibrium (z_bed_eq=bed[k0], h_ice_ref=ice[k0]).
   !!   i_eq=1  (default) present-day reference: z_bed_eq / h_ice_ref from
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
   use fe_sht,       only: sht_grid
   use fe_coupling,  only: solid_earth
   use fe_remap,     only: ll2gauss_map, ll2gauss_init, ll2gauss_apply
   use fe_io,        only: fe_write_step
   use ncio,         only: nc_read, nc_size
   implicit none
   private

   public :: fastearth_run

   ! diagnostic surface fields written each output step (the prognostic memory is
   ! written separately via fe_restart_write when a restart is wanted)
   character(len=8), parameter :: OUT_VARS(5) = &
        [character(len=8) :: "h_ice", "rsl", "z_bed", "C_ocean", "bsl"]

   integer,  parameter :: MAX_EQ  = 12       !! cap on paleotopo fixed-point passes
   real(wp), parameter :: EQ_TOL  = 1.0_wp   !! equilibration residual mean|.| [m]

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
      real(wp), allocatable  :: z_bed_eq(:,:), h_ice_ref(:,:), h_ice(:,:), ice_lgm(:,:)
      real(wp), allocatable  :: src(:,:), tyr(:)
      real(wp) :: t0, dt
      integer  :: nt, k, k0, k1, np, nl, nlon, nls
      logical  :: remap
      integer(kind=8) :: pc0, pc1, prate          ! PROFILE: per-step phase timers
      real(wp) :: t_read = 0.0_wp, t_upd = 0.0_wp, t_wrt = 0.0_wp
      integer  :: nstep = 0

      ! --- configuration --------------------------------------------------------
      if (present(defaults_file)) then
         call fe_par_load(p, cfg_file, defaults_file=defaults_file)
      else
         call fe_par_load(p, cfg_file)
      end if
      call fe_par_print(p)

      ! --- transform grid + work arrays -----------------------------------------
      call build_grid(p, sht)
      np = sht%nphi;  nl = sht%nlat
      allocate(z_bed_eq(np,nl), h_ice_ref(np,nl), h_ice(np,nl), ice_lgm(np,nl))

      if (len_trim(p%file_forcing) == 0) error stop 'fastearth_run: file_forcing not set'

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

      ! --- start-slice ice + reference state (z_bed_eq, h_ice_ref) per i_eq ------
      call read_ice(p, rmap, sht, remap, k0, src, ice_lgm)      ! start-slice (LGM) ice
      call setup_reference(p, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_ref)

      ! --- initialise; optional ice-free paleotopo spin-up (dt_equil>0) ---------
      call system_clock(pc0, prate)
      if (p%dt_equil > 0.0_wp) then
         ! non-default: replace z_bed_eq with the ice-free relaxed bed (paleotopo
         ! fixed point) and leave se spun up to the start-ice equilibrium.
         call equilibrate(p, sht, se, z_bed_eq, ice_lgm)
      else
         call se%init(p, sht, z_bed_eq, h_ice_ref)
         call se%update(ice_lgm, 0.0_wp)                        ! seed entering ice, no integration
      end if
      call system_clock(pc1)
      write(*,'(a,f8.2,a)') ' [PROFILE setup] se%init+visc3d+seed =', real(pc1-pc0,wp)/prate, ' s'

      ! 1-D/3-D layer split diagnostic: how many elements are genuinely laterally 3-D
      ! (pay the pseudo-spectral tensor-SH advance) vs collapse to the cheap 1-D path.
      if (p%l_visc_3d) write(*,'(a,i0,a,i0,a)') ' visc3d split: ', se%resp%ne3d, &
           ' of ', se%resp%ne, ' radial elements laterally 3-D (rest advance as 1-D)'

      ! --- march the transient --------------------------------------------------
      t0 = tyr(k0)*sec_per_year
      se%time = t0;  se%resp%time = t0
      call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.true.)

      do k = k0, k1-1
         dt = (tyr(k+1) - tyr(k))*sec_per_year
         call system_clock(pc0, prate)
         call read_ice(p, rmap, sht, remap, k+1, src, h_ice)
         call system_clock(pc1);  t_read = t_read + real(pc1-pc0,wp)/prate
         call system_clock(pc0)
         call se%update(h_ice, dt)
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
         '   se%update (SLE+adv) =', 1.0e3_wp*t_upd /nstep, ' ms (', &
            100.0_wp*t_upd /(t_read+t_upd+t_wrt), ' %)', &
         '   fe_write_step (out) =', 1.0e3_wp*t_wrt /nstep, ' ms (', &
            100.0_wp*t_wrt /(t_read+t_upd+t_wrt), ' %)'
      if (nstep > 0) write(*,'(a,f7.1,a,f7.1,a)') &
         '   sub-steps/interval: n_accept=', real(se%stepper%n_accept,wp)/nstep, &
         '  n_solve=', real(se%stepper%n_solve,wp)/nstep, '  (per coupling step)'
      write(*,'(a,a)') ' fastearth: wrote ', trim(p%file_out)
      call se%finalize();  call sht%destroy()
   end subroutine fastearth_run

   ! --- equilibration ----------------------------------------------------------

   subroutine equilibrate(p, sht, se, z_bed_eq, ice_lgm)
      !! Paleotopo fixed point (i_eq=1): find the ice-free relaxed bed z_bed_eq whose
      !! viscoelastic EQUILIBRIUM under the start-slice ice ice_lgm reproduces the data
      !! bed (the entering z_bed_eq). Each pass re-inits the model ice-free (memory 0),
      !! holds ice_lgm for dt_equil to relax, and Newton-updates z_bed_eq by the bed
      !! residual (the deflection is ~independent of the datum, so this converges in a
      !! couple of passes). On return z_bed_eq is the ice-free relaxed bed and se holds
      !! the spun-up equilibrium state (bed + memory) ready for the transient.
      type(fe_param_class), intent(in)            :: p
      type(sht_grid),       intent(in),   target  :: sht
      type(solid_earth),    intent(inout)         :: se
      real(wp),             intent(inout)         :: z_bed_eq(:,:)  !! in: data bed; out: ice-free bed
      real(wp),             intent(in)            :: ice_lgm(:,:)
      real(wp), allocatable :: bed_target(:,:), h0(:,:), resid(:,:)
      real(wp) :: rmean, rmax, fourpi
      integer  :: it, np, nl

      np = sht%nphi;  nl = sht%nlat
      fourpi = 16.0_wp*atan(1.0_wp)
      allocate(bed_target(np,nl), source=z_bed_eq)
      allocate(h0(np,nl), source=0.0_wp)
      allocate(resid(np,nl))

      write(*,'(a,es9.2,a)') ' equilibration (i_eq=1): spin up under start ice, dt_equil=', &
           p%dt_equil/sec_per_year, ' yr/pass'
      do it = 1, MAX_EQ
         call se%init(p, sht, z_bed_eq, h0)              ! ice-free reference, memory 0
         call se%update(ice_lgm, 0.0_wp)                 ! set entering ice = start ice
         call se%update(ice_lgm, p%dt_equil)             ! hold -> relax to equilibrium
         resid = se%z_bed - bed_target
         rmean = sht%surface_integral(abs(resid))/fourpi
         rmax  = maxval(abs(resid))
         write(*,'(a,i2,a,f10.4,a,f10.2,a)') '   pass ', it, ':  mean|z_bed-bed_data|=', &
              rmean, ' m   max=', rmax, ' m'
         if (rmean < EQ_TOL) exit
         z_bed_eq = z_bed_eq - resid                     ! Newton step (deflection ~ datum-independent)
      end do
   end subroutine equilibrate

   ! --- input helpers ----------------------------------------------------------

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

   subroutine setup_reference(p, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_ref)
      !! Fill the relaxed reference (z_bed_eq = SLE topo0, h_ice_ref) per p%i_eq,
      !! mirroring CLIMBER-X i_equilibrium. Reference files are lon-lat and remapped
      !! online with a per-file conservative map (their grid differs from the forcing).
      type(fe_param_class), intent(in)    :: p
      type(ll2gauss_map),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k0
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(in)    :: ice_lgm(:,:)
      real(wp),             intent(out)   :: z_bed_eq(:,:), h_ice_ref(:,:)
      real(wp), allocatable :: rsl_r(:,:)
      select case (p%i_eq)
      case (0)
         call read_bed(p, rmap, sht, remap, k0, src, z_bed_eq)   ! data bed at the start slice
         h_ice_ref = ice_lgm
      case (1)
         call read_ref2d(p, sht, remap, p%z_bed_ref_file, p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%h_ice_ref_file, p%name_h_ice_ref, .true.,  h_ice_ref)
      case (2)
         call read_ref2d(p, sht, remap, p%z_bed_eq_file, p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%h_ice_eq_file, p%name_h_ice_ref, .true.,  h_ice_ref)
      case (3)
         allocate(rsl_r(size(z_bed_eq,1), size(z_bed_eq,2)))
         call read_ref2d(p, sht, remap, p%z_bed_ref_file,   p%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(p, sht, remap, p%rsl_restart_file, p%name_rsl,       .false., rsl_r)
         call read_ref2d(p, sht, remap, p%h_ice_ref_file,   p%name_h_ice_ref, .true.,  h_ice_ref)
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
         call sht%init(p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
      else if (p%mmax >= 0) then
         call sht%init(p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
      else if (p%eps_polar > 0.0_wp) then
         call sht%init(p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
      else
         call sht%init(p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
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
