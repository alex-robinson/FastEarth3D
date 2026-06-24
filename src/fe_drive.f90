module fe_drive
   !! Standalone forced-run driver. Given a &fe3d configuration (fe_params), build
   !! the transform grid and the solid-Earth model, read a reference equilibrium
   !! state and an ice-thickness forcing time series from netCDF, and march the
   !! model across the forcing, writing the diagnostic surface fields each step.
   !!
   !! This is the same model a host (CLIMBER-X) drives through fe_coupling — only
   !! the I/O and the time loop live here, so `program fastearth` is a thin wrapper
   !! around fastearth_run. Inputs are assumed already on the model's Gauss-Legendre
   !! grid (the host does any remapping); the forcing time axis (in years) defines
   !! the coupling cadence — one se%update per slice, advancing by the slice spacing.
   use fe_precision, only: wp
   use fe_constants, only: sec_per_year
   use fe_params,    only: fe_param_class, fe_par_load, fe_par_print
   use fe_sht,       only: sht_grid
   use fe_coupling,  only: solid_earth
   use fe_io,        only: fe_write_step
   use ncio,         only: nc_read, nc_size
   implicit none
   private

   public :: fastearth_run

   ! diagnostic surface fields written each output step (the prognostic memory is
   ! written separately via fe_restart_write when a restart is wanted)
   character(len=8), parameter :: OUT_VARS(4) = &
        [character(len=8) :: "h_ice", "rsl", "z_bed", "C_ocean"]

contains

   subroutine fastearth_run(cfg_file, defaults_file)
      !! Run a forced simulation defined by the &fe3d group of cfg_file (overlaid on
      !! defaults_file if given; otherwise cfg_file must be complete).
      character(len=*), intent(in)           :: cfg_file
      character(len=*), intent(in), optional :: defaults_file

      type(fe_param_class)   :: p
      type(sht_grid), target :: sht
      type(solid_earth)      :: se
      real(wp), allocatable  :: z_bed_eq(:,:), h_ice_ref(:,:), h_ice(:,:), tyr(:)
      real(wp) :: t0, dt
      integer  :: nt, k, k0, k1, np, nl

      ! --- configuration --------------------------------------------------------
      if (present(defaults_file)) then
         call fe_par_load(p, cfg_file, defaults_file=defaults_file)
      else
         call fe_par_load(p, cfg_file)
      end if
      call fe_par_print(p)

      ! --- transform grid + model ----------------------------------------------
      call build_grid(p, sht)
      np = sht%nphi;  nl = sht%nlat
      allocate(z_bed_eq(np,nl), h_ice_ref(np,nl), h_ice(np,nl))

      if (len_trim(p%file_ref) == 0)     error stop 'fastearth_run: file_ref not set'
      if (len_trim(p%file_forcing) == 0) error stop 'fastearth_run: file_forcing not set'
      call nc_read(p%file_ref, trim(p%name_zbed_eq),  z_bed_eq)
      call nc_read(p%file_ref, trim(p%name_hice_ref), h_ice_ref)

      call se%init(p, sht, z_bed_eq, h_ice_ref)

      ! --- forcing time axis (years) + window -----------------------------------
      nt = nc_size(p%file_forcing, trim(p%name_time))
      allocate(tyr(nt));  call nc_read(p%file_forcing, trim(p%name_time), tyr)
      call select_window(tyr*sec_per_year, p%time_init, p%time_end, k0, k1)
      if (k1 <= k0) error stop 'fastearth_run: forcing window contains < 2 time slices'

      write(*,'(a,i0,a,f0.1,a,f0.1,a)') ' fastearth: ', k1-k0+1, ' slices, t = ', &
           tyr(k0), ' -> ', tyr(k1), ' yr'

      ! --- march the model ------------------------------------------------------
      ! Align the model clock to the first forcing time (absolute years) and seed
      ! the entering load with slice k0 (a zero-length update sets se%h_ice without
      ! integrating), then write the initial state.
      t0 = tyr(k0)*sec_per_year
      se%time = t0;  se%resp%time = t0
      call read_slice(p, k0, h_ice)
      call se%update(h_ice, 0.0_wp)
      call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.true.)

      do k = k0, k1-1
         dt = (tyr(k+1) - tyr(k))*sec_per_year
         call read_slice(p, k+1, h_ice)
         call se%update(h_ice, dt)
         call fe_write_step(se, p%file_out, se%time, nms=OUT_VARS, init=.false.)
         write(*,'(a,f12.2,a,es10.2,a,es9.2)') '   t=', se%time/sec_per_year, &
              ' yr   max|rsl|=', maxval(abs(se%rsl)), '   mass_resid=', se%worst_mass_resid
      end do

      write(*,'(a,a)') ' fastearth: wrote ', trim(p%file_out)
      call se%finalize();  call sht%destroy()
   end subroutine fastearth_run

   ! --- helpers ----------------------------------------------------------------

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

   subroutine read_slice(p, k, h_ice)
      !! Read forcing ice-thickness slice k, h_ice(lon,lat) at time index k.
      type(fe_param_class), intent(in)  :: p
      integer,              intent(in)  :: k
      real(wp),             intent(out) :: h_ice(:,:)
      call nc_read(p%file_forcing, trim(p%name_ice), h_ice, &
                   start=[1,1,k], count=[size(h_ice,1), size(h_ice,2), 1])
   end subroutine read_slice

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
