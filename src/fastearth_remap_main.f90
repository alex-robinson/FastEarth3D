program fastearth_remap
   !! Offline lon-lat -> Gauss-Legendre remapper. Reads a &fe3d config (the grid
   !! knobs lmax/nlat/nphi, the source file_forcing + name_ice/name_lon/name_lat/
   !! name_time, and file_out), conservatively remaps the ice variable over every
   !! time slice onto the model Gauss grid, and writes a Gauss-grid forcing file the
   !! standalone driver can consume with remap_input=.false.
   !!
   !!   ./bin/fastearth_remap.x config.nml [defaults.nml]
   !!
   !! The online path (fe_drive, remap_input=.true.) is the default and needs no
   !! preprocessing; this tool is for pre-baking a remapped forcing when that is
   !! preferred (e.g. reusing one remap across many runs). Same conservative engine
   !! (fe_remap), so results are identical.
   use fe_precision, only: wp
   use fe_params,    only: fe_param_class, fe_par_load
   use fe_sht,       only: sht_grid
   use fe_remap,     only: ll2gauss_map
   use ncio,         only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
   implicit none

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   character(len=512)  :: cfg, defs
   type(fe_param_class) :: p
   type(sht_grid)       :: sht
   type(ll2gauss_map)   :: rmap
   real(wp), allocatable :: lon_s(:), lat_s(:), src(:,:), gauss(:,:)
   real(wp), allocatable :: lon_g(:), lat_g(:), tyr(:)
   integer :: nlon, nls, np, nl, nt, k, nlat, nphi

   ! --- config ---------------------------------------------------------------
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "fastearth.nml"
   end if
   if (command_argument_count() >= 2) then
      call get_command_argument(2, defs)
      call fe_par_load(p, cfg, defaults_file=trim(defs))
   else
      call fe_par_load(p, cfg)
   end if
   if (len_trim(p%file_forcing) == 0) error stop 'fastearth_remap: file_forcing not set'

   ! --- Gauss grid (same defaulting as the driver) ---------------------------
   nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
   nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
   call sht%init(p%lmax, nlat=nlat, nphi=nphi)
   np = sht%nphi;  nl = sht%nlat

   ! --- conservative map from the source axes --------------------------------
   nlon = nc_size(p%file_forcing, trim(p%name_lon))
   nls  = nc_size(p%file_forcing, trim(p%name_lat))
   allocate(lon_s(nlon), lat_s(nls))
   call nc_read(p%file_forcing, trim(p%name_lon), lon_s)
   call nc_read(p%file_forcing, trim(p%name_lat), lat_s)
   call rmap%init(sht, lon_s, lat_s)

   nt = nc_size(p%file_forcing, trim(p%name_time))
   allocate(tyr(nt));  call nc_read(p%file_forcing, trim(p%name_time), tyr)
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') ' fastearth_remap: ', nlon, 'x', nls, &
        ' lon-lat -> ', np, 'x', nl, ' Gauss, ', nt, ' slices'

   ! --- Gauss output axes ----------------------------------------------------
   allocate(lon_g(np), lat_g(nl))
   lon_g = sht%lon*RAD2DEG
   do k = 1, nl;  lat_g(k) = 90.0_wp - sht%colat(k)*RAD2DEG;  end do

   call nc_create(p%file_out)
   call nc_write_dim(p%file_out, "lon",  x=lon_g, units="degrees_east")
   call nc_write_dim(p%file_out, "lat",  x=lat_g, units="degrees_north")
   call nc_write_dim(p%file_out, "time", x=tyr,   units="years", unlimited=.true.)

   allocate(src(nlon,nls), gauss(np,nl))
   do k = 1, nt
      call nc_read(p%file_forcing, trim(p%name_ice), src, &
                   start=[1,1,k], count=[nlon, nls, 1])
      call rmap%apply(sht, src, gauss, conserve_mass=.true.)
      call nc_write(p%file_out, trim(p%name_ice), gauss, dim1="lon", dim2="lat", &
                    dim3="time", start=[1,1,k], count=[np, nl, 1])
   end do

   call sht%destroy()
   write(*,'(a,a)') ' fastearth_remap: wrote ', trim(p%file_out)
end program fastearth_remap
