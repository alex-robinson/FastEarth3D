program fastearth_mkref
   !! Offline prebake of the 2D reference state (bedrock + ice) onto the model
   !! Gauss grid. Reads a &fe3d config (grid knobs lmax/nlat/nphi, the reference
   !! source z_bed_ref_file / h_ice_ref_file + name_z_bed_ref / name_h_ice_ref /
   !! name_lon / name_lat, and file_out), conservatively remaps bed (as-is) and ice
   !! (mass-conserving) once, and writes a small Gauss-grid reference file.
   !!
   !!   ./bin/fastearth_mkref.x prebake.nml [defaults.nml]
   !!
   !! The driver (fe_drive, i_eq=1/2/3) reads such a file directly — its lon/lat
   !! dims match the Gauss grid, so read_ref2d skips the conservative-map build. The
   !! same fe_remap engine is used here, so the prebaked field is identical to the
   !! online remap. Both reference vars are assumed on a common source grid (RTopo).
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
   real(wp), allocatable :: lon_s(:), lat_s(:), buf(:,:), bed(:,:), ice(:,:)
   real(wp), allocatable :: lon_g(:), lat_g(:)
   integer :: nlon, nls, np, nl, k, nlat, nphi

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
   if (len_trim(p%z_bed_ref_file) == 0) error stop 'fastearth_mkref: z_bed_ref_file not set'
   if (len_trim(p%h_ice_ref_file) == 0) error stop 'fastearth_mkref: h_ice_ref_file not set'

   ! --- Gauss grid (same defaulting as the driver) ---------------------------
   nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
   nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
   call sht%init(p%lmax, nlat=nlat, nphi=nphi)
   np = sht%nphi;  nl = sht%nlat

   ! --- conservative map from the reference source axes ----------------------
   nlon = nc_size(p%z_bed_ref_file, trim(p%name_lon))
   nls  = nc_size(p%z_bed_ref_file, trim(p%name_lat))
   allocate(lon_s(nlon), lat_s(nls), buf(nlon,nls), bed(np,nl), ice(np,nl))
   call nc_read(p%z_bed_ref_file, trim(p%name_lon), lon_s)
   call nc_read(p%z_bed_ref_file, trim(p%name_lat), lat_s)
   write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' fastearth_mkref: ', nlon, 'x', nls, &
        ' lon-lat -> ', np, 'x', nl, ' Gauss (building conservative map)'
   call rmap%init(sht, lon_s, lat_s)

   call nc_read(p%z_bed_ref_file, trim(p%name_z_bed_ref), buf)
   call rmap%apply(sht, buf, bed, conserve_mass=.false.)
   call nc_read(p%h_ice_ref_file, trim(p%name_h_ice_ref), buf)
   call rmap%apply(sht, buf, ice, conserve_mass=.true.)

   ! --- Gauss output ---------------------------------------------------------
   allocate(lon_g(np), lat_g(nl))
   lon_g = sht%lon*RAD2DEG
   do k = 1, nl;  lat_g(k) = 90.0_wp - sht%colat(k)*RAD2DEG;  end do

   call nc_create(p%file_out)
   call nc_write_dim(p%file_out, "lon", x=lon_g, units="degrees_east")
   call nc_write_dim(p%file_out, "lat", x=lat_g, units="degrees_north")
   call nc_write(p%file_out, trim(p%name_z_bed_ref), bed, dim1="lon", dim2="lat", &
                 units="m", long_name="reference bedrock elevation (Gauss grid)")
   call nc_write(p%file_out, trim(p%name_h_ice_ref), ice, dim1="lon", dim2="lat", &
                 units="m", long_name="reference ice thickness (Gauss grid)")

   call sht%destroy()
   write(*,'(a,a)') ' fastearth_mkref: wrote ', trim(p%file_out)
end program fastearth_mkref
