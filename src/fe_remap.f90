module fe_remap
   !! Conservative lon-lat -> Gauss-Legendre remapping for the standalone driver
   !! (fe_drive). Real forcing (e.g. CLIMBER-X geo_ice_tarasov_deglac.nc) lives on a
   !! regular lon-lat grid; the model runs on the SHTns Gauss grid. This wraps the
   !! fesm-utils `coords` library (great-circle polygon clipping; the same SCRIP-style
   !! apparatus CLIMBER-X uses host-side) to build a conservative source->Gauss map
   !! once, then apply it per field / time slice.
   !!
   !! The coupling contract (fe_coupling) assumes the HOST remaps onto the Gauss grid,
   !! so this is a driver-only concern; nothing in the coupled path depends on it.
   !!
   !! Mass note: `coords` derives target cell boundaries from axis midpoints, so the
   !! Gauss-cell areas it uses are not exactly the SHTns Gauss quadrature weights. The
   !! map is conservative w.r.t. its own areas, but when SHTns re-integrates the
   !! remapped field its quadrature total differs by O(1/nlat^2). For a mass-bearing
   !! field (ice) pass conserve_mass=.true.: a single global factor rescales the result
   !! so its SHTns surface integral equals the source area-integral exactly (the whole
   !! sphere is 4*pi sr, used to convert the conserved coords-area total to steradians
   !! without needing the planet radius). Geometry fields (bed) are remapped as-is.
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid, sht_grid_surface_integral
   use coords,       only: grid_class, grid_init, map_class, &
                           map_init_conservative, map_field
   implicit none
   private

   public :: ll2gauss_map
   public :: ll2gauss_init, ll2gauss_apply

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   real(wp), parameter :: FOURPI  = 12.566370614359172_wp

   type :: ll2gauss_map
      !! A precomputed conservative map from a lon-lat source grid onto the model's
      !! SHTns Gauss grid. Build once with %init, then %apply to each field.
      type(grid_class) :: src                 !! lon-lat source
      type(grid_class) :: dst                 !! Gauss target (matches sht), lat ascending
      type(map_class)  :: map                 !! conservative weights src -> dst
      integer :: nlon = 0, nlat_src = 0       !! source dimensions
      integer :: nphi = 0, nlat = 0           !! target (SHTns) dimensions
   end type ll2gauss_map

contains

   subroutine ll2gauss_init(self, sht, lon_src, lat_src)
      !! Build the conservative map for source cell-centre axes lon_src(:), lat_src(:)
      !! [degrees, ascending] onto the SHTns Gauss grid. The target grid is built with
      !! ASCENDING latitude (south first), which coords expects; %apply flips back into
      !! the SHTns (north-first) row order.
      type(ll2gauss_map), intent(out) :: self
      type(sht_grid),      intent(in)  :: sht
      real(wp),            intent(in)  :: lon_src(:)   !! source longitudes [deg]
      real(wp),            intent(in)  :: lat_src(:)   !! source latitudes  [deg]
      real(wp), allocatable :: lon_t(:), lat_t(:)
      integer :: j

      self%nlon = size(lon_src);  self%nlat_src = size(lat_src)
      self%nphi = sht%nphi;       self%nlat = sht%nlat

      allocate(lon_t(self%nphi), lat_t(self%nlat))
      do j = 1, self%nphi
         lon_t(j) = sht%lon(j)*RAD2DEG               ! SHTns longitudes [0,360)
      end do
      do j = 1, self%nlat
         ! SHTns colat is ascending (north -> south), so lat descends; reverse it so
         ! the target axis is ascending (south -> north).
         lat_t(j) = 90.0_wp - sht%colat(self%nlat - j + 1)*RAD2DEG
      end do

      call grid_init(self%src, name="ll-src", mtype="latlon", units="degrees", &
                     lon180=.true.,  x=lon_src, y=lat_src)
      call grid_init(self%dst, name="gauss",  mtype="latlon", units="degrees", &
                     lon180=.false., x=lon_t,   y=lat_t)
      call map_init_conservative(self%map, self%src, self%dst)
   end subroutine ll2gauss_init

   subroutine ll2gauss_apply(self, sht, fsrc, fdst, conserve_mass)
      !! Remap fsrc(nlon, nlat_src) onto fdst(nphi, nlat) in the SHTns spatial layout.
      !! conserve_mass (default .false.): rescale fdst so its SHTns surface integral
      !! equals the source area-integral (use for ice thickness; leave off for bed).
      type(ll2gauss_map), intent(in)  :: self
      type(sht_grid),      intent(in)  :: sht
      real(wp),            intent(in)  :: fsrc(:,:)    !! (nlon, nlat_src)
      real(wp),            intent(out) :: fdst(:,:)    !! (nphi, nlat)
      logical, optional,   intent(in)  :: conserve_mass
      real(wp), allocatable :: vt(:,:)
      logical,  allocatable :: m2(:,:)
      real(wp) :: src_int_sr, gauss_int
      integer  :: j
      logical  :: do_mass

      if (size(fsrc,1) /= self%nlon .or. size(fsrc,2) /= self%nlat_src) &
         error stop 'fe_remap: source field shape /= map source grid'
      if (size(fdst,1) /= self%nphi .or. size(fdst,2) /= self%nlat) &
         error stop 'fe_remap: target field shape /= SHTns grid'

      allocate(vt(self%nphi, self%nlat), m2(self%nphi, self%nlat))
      call map_field(self%map, "f", fsrc, vt, method="mean", mask2=m2)
      where (.not. m2) vt = 0.0_wp                    ! uncovered cells (global src: none)

      do j = 1, self%nlat                              ! ascending-lat -> SHTns north-first
         fdst(:, j) = vt(:, self%nlat - j + 1)
      end do

      do_mass = .false.;  if (present(conserve_mass)) do_mass = conserve_mass
      if (do_mass) then
         ! conserved source total as a steradian integral: the map preserves
         ! sum(fsrc*src_area); convert to sr via 4*pi/total_area (full-sphere grid).
         src_int_sr = sum(fsrc * self%src%area) / sum(self%src%area) * FOURPI
         gauss_int  = sht_grid_surface_integral(sht, fdst)
         if (abs(gauss_int) > tiny(1.0_wp)) fdst = fdst * (src_int_sr/gauss_int)
      end if
   end subroutine ll2gauss_apply

end module fe_remap
