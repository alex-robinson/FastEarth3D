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
   !! Reference / equilibration (i_eq):
   !!   i_eq=0  declare the first in-window slice as the relaxed reference (z_bed_eq =
   !!           bed[k0], h_ice_ref = ice[k0], memory 0); the transient load is the ice
   !!           change relative to that slice.
   !!   i_eq=1  (default) ice-free reference; find the ice-free relaxed bed by a
   !!           paleotopo fixed point so the EQUILIBRIUM under the start-slice ice
   !!           reproduces the data bed[k0], holding the load dt_equil per pass. Leaves
   !!           the model at that equilibrium (bed AND viscous memory), then deglaciates
   !!           with the absolute ice as load. The reported residual ‖z_bed-bed[k0]‖ is
   !!           the consistency check vs i_eq=0.
   use fe_precision, only: wp
   use fe_constants, only: sec_per_year
   use fe_params,    only: fe_param_class, fe_par_load, fe_par_print
   use fe_sht,       only: sht_grid
   use fe_coupling,  only: solid_earth
   use fe_remap,     only: ll2gauss_map
   use fe_io,        only: fe_write_step
   use ncio,         only: nc_read, nc_size
   implicit none
   private

   public :: fastearth_run

   ! diagnostic surface fields written each output step (the prognostic memory is
   ! written separately via fe_restart_write when a restart is wanted)
   character(len=8), parameter :: OUT_VARS(4) = &
        [character(len=8) :: "h_ice", "rsl", "z_bed", "C_ocean"]

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
      if (remap) then
         call build_remap(p, sht, rmap, nlon, nls)
         allocate(src(nlon, nls))
         write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' remap: lon-lat (', nlon, 'x', nls, &
              ') -> Gauss (', np, 'x', nl, ')'
      end if

      ! --- forcing time axis (years) + window -----------------------------------
      nt = nc_size(p%file_forcing, trim(p%name_time))
      allocate(tyr(nt));  call nc_read(p%file_forcing, trim(p%name_time), tyr)
      call select_window(tyr*sec_per_year, p%time_init, p%time_end, k0, k1)
      if (k1 <= k0) error stop 'fastearth_run: forcing window contains < 2 time slices'
      write(*,'(a,i0,a,f0.1,a,f0.1,a)') ' fastearth: ', k1-k0+1, ' slices, t = ', &
           tyr(k0), ' -> ', tyr(k1), ' yr'

      ! --- reference bed + start-slice ice --------------------------------------
      call read_bed(p, rmap, sht, remap, k0, src, z_bed_eq)     ! data bed at the start slice
      call read_ice(p, rmap, sht, remap, k0, src, ice_lgm)      ! start-slice (LGM) ice

      ! --- reference / equilibration --------------------------------------------
      if (p%i_eq == 0) then
         ! declare the start slice as the relaxed reference (memory 0)
         h_ice_ref = ice_lgm
         call se%init(p, sht, z_bed_eq, h_ice_ref)
         call se%update(ice_lgm, 0.0_wp)                        ! seed entering ice, no integration
      else
         ! ice-free reference + paleotopo spin-up; leaves se at the LGM equilibrium
         h_ice_ref = 0.0_wp
         call equilibrate(p, sht, se, z_bed_eq, ice_lgm)        ! z_bed_eq -> ice-free relaxed bed
      end if

      ! --- march the transient --------------------------------------------------
      t0 = tyr(k0)*sec_per_year
      se%time = t0;  se%resp%time = t0
      call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.true.)

      do k = k0, k1-1
         dt = (tyr(k+1) - tyr(k))*sec_per_year
         call read_ice(p, rmap, sht, remap, k+1, src, h_ice)
         call se%update(h_ice, dt)
         call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.false.)
         write(*,'(a,f12.2,a,es10.2,a,es9.2)') '   t=', se%time/sec_per_year, &
              ' yr   max|rsl|=', maxval(abs(se%rsl)), '   mass_resid=', se%worst_mass_resid
      end do

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
      call rmap%init(sht, lon_s, lat_s)
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
         call rmap%apply(sht, src, h_ice, conserve_mass=.true.)
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
         call rmap%apply(sht, src, z_bed, conserve_mass=.false.)
      else
         if (len_trim(p%file_ref) == 0) error stop 'fastearth_run: file_ref not set (remap_input=.false.)'
         call nc_read(p%file_ref, trim(p%name_zbed_eq), z_bed)
      end if
   end subroutine read_bed

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
