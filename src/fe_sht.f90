module fe_sht
   !! Thin Fortran wrapper around the SHTns spherical-harmonic transform library.
   !!
   !! Horizontal transforms are the performance kernel of the spectral-finite-
   !! element solver (Martinec 2000): every time step maps fields between the
   !! Gauss-Legendre spatial grid (where lateral viscosity and the sea-level
   !! equation live) and the spectral coefficients (where the radial solves are
   !! per-degree and decoupled). This module isolates the SHTns C API behind a
   !! small derived type so the rest of FastEarth3D never touches iso_c_binding.
   !!
   !! Convention: fully-normalized real spherical harmonics, no Condon-Shortley
   !! phase (SHT_ORTHONORMAL + SHT_NO_CS_PHASE), on a Gauss grid with a
   !! phi-contiguous spatial layout. Spectral arrays hold m >= 0 only (the field
   !! is real). See doc/design.md for the rationale and the degree-1 / frame
   !! caveats that the load Love numbers will later depend on.
   use, intrinsic :: iso_c_binding
   use fe_precision, only: wp
   implicit none
   private

   ! SHTns Fortran 2003 interface (installed by `make install` into the SHTns
   ! prefix include dir; found via -I$(SHTNSROOT)/include).
   include 'shtns.f03'

   real(wp), parameter :: pi = acos(-1.0_wp)

   public :: sht_grid

   type :: sht_grid
      !! One configured SHTns transform (a spectral truncation + matching grid).
      type(c_ptr) :: cfg = c_null_ptr   !! opaque SHTns configuration handle
      integer :: lmax = 0               !! max spherical-harmonic degree
      integer :: mmax = 0               !! max order (in units of mres)
      integer :: mres = 1              !! order step
      integer :: nlat = 0              !! # latitudes (Gauss nodes)
      integer :: nphi = 0              !! # longitudes
      integer :: nlm  = 0              !! # spectral coefficients (m >= 0)
      ! Physical grid geometry, cached at init. Spatial fields are (nphi, nlat).
      real(wp), allocatable :: colat(:)    !! colatitude of each latitude row [rad] (nlat)
      real(wp), allocatable :: lon(:)      !! longitude of each column [rad] (nphi)
      real(wp), allocatable :: gauss_w(:)  !! Gauss quadrature weights [-] (nlat), sum = 2
   contains
      procedure :: init      => sht_grid_init
      procedure :: destroy   => sht_grid_destroy
      procedure :: synthesis => sht_grid_synthesis   !! spectral -> spatial
      procedure :: analysis  => sht_grid_analysis    !! spatial  -> spectral
      procedure :: sph_synthesis => sht_grid_sph_synthesis !! spheroidal (gradient) vector synth
      procedure :: eval_point => sht_grid_eval_point !! scalar field at (colat,lon)
      procedure :: eval_point_horiz => sht_grid_eval_point_horiz !! ∇₁ field at (colat,lon)
      procedure :: lmidx     => sht_grid_lmidx       !! (l,m) -> coefficient index
      procedure :: surface_integral => sht_grid_surface_integral  !! ∫ f dΩ
   end type sht_grid

contains

   subroutine sht_grid_init(self, lmax, nlat, nphi, mmax, mres, eps_polar)
      !! Create the SHTns config and its Gauss grid.
      !!
      !! Grid defaults to the smallest that resolves degree `lmax` exactly with
      !! Gauss quadrature (nlat >= lmax+1, nphi >= 2*mmax+1). For quadratic
      !! products (the sea-level equation's ocean-function multiply) the caller
      !! should pass a de-aliased grid: nlat ~ 3*lmax/2, nphi ~ 3*lmax.
      class(sht_grid), intent(inout) :: self
      integer,  intent(in)           :: lmax
      integer,  intent(in), optional :: nlat, nphi, mmax, mres
      real(wp), intent(in), optional :: eps_polar

      type(shtns_info), pointer :: info
      real(c_double) :: eps
      integer :: mmax_, mres_, nlat_, nphi_, norm, layout

      mres_ = 1;          if (present(mres)) mres_ = mres
      mmax_ = lmax/mres_; if (present(mmax)) mmax_ = mmax
      eps   = 1.0e-10_c_double
      if (present(eps_polar)) eps = real(eps_polar, c_double)

      nlat_ = lmax + 2;       if (present(nlat)) nlat_ = nlat
      nphi_ = 2*(mmax_ + 1);  if (present(nphi)) nphi_ = nphi

      norm   = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
      layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS

      self%lmax = lmax
      self%mmax = mmax_
      self%mres = mres_

      self%cfg = shtns_create(lmax, mmax_, mres_, norm)
      call shtns_set_grid(self%cfg, layout, eps, nlat_, nphi_)

      ! Read back the realized grid sizes from the SHTns info struct.
      call c_f_pointer(self%cfg, info)
      self%nlat = info%nlat
      self%nphi = info%nphi
      self%nlm  = info%nlm

      call cache_geometry(self, info)
   end subroutine sht_grid_init

   subroutine cache_geometry(self, info)
      !! Cache the physical grid coordinates and quadrature weights from SHTns.
      class(sht_grid),  intent(inout) :: self
      type(shtns_info), intent(in)    :: info
      real(c_double), pointer     :: cos_theta(:)
      real(c_double), allocatable :: wts_half(:)
      integer :: j, nh

      ! Colatitudes from the cos(theta) array SHTns exposes on its struct.
      call c_f_pointer(info%ct, cos_theta, [self%nlat])
      self%colat = acos(real(cos_theta, wp))

      ! Uniform longitudes.
      allocate(self%lon(self%nphi))
      do j = 1, self%nphi
         self%lon(j) = 2.0_wp*pi*real(j-1, wp)/real(self%nphi, wp)
      end do

      ! Gauss weights: SHTns returns one hemisphere (nlat_2); mirror to full grid.
      nh = info%nlat_2
      allocate(wts_half(nh))
      call shtns_gauss_wts(self%cfg, wts_half)
      allocate(self%gauss_w(self%nlat))
      self%gauss_w(1:nh) = real(wts_half, wp)
      self%gauss_w(self%nlat:self%nlat-nh+1:-1) = real(wts_half, wp)
      deallocate(wts_half)
   end subroutine cache_geometry

   subroutine sht_grid_destroy(self)
      !! Release the SHTns configuration and cached geometry.
      class(sht_grid), intent(inout) :: self
      if (c_associated(self%cfg)) call shtns_destroy(self%cfg)
      self%cfg = c_null_ptr
      self%lmax = 0; self%mmax = 0; self%nlat = 0; self%nphi = 0; self%nlm = 0
      if (allocated(self%colat))   deallocate(self%colat)
      if (allocated(self%lon))     deallocate(self%lon)
      if (allocated(self%gauss_w)) deallocate(self%gauss_w)
   end subroutine sht_grid_destroy

   real(wp) function sht_grid_surface_integral(self, f) result(total)
      !! Surface integral ∫ f dΩ over the unit sphere by Gauss-Legendre
      !! quadrature in latitude and the uniform rule in longitude. Used for the
      !! mass/area integrals the sea-level equation depends on; for f≡1 it
      !! returns 4π.
      class(sht_grid), intent(in) :: self
      real(wp),        intent(in) :: f(:,:)   !! (nphi, nlat)
      real(wp) :: dlon, row
      integer  :: i, j
      dlon  = 2.0_wp*pi/real(self%nphi, wp)
      total = 0.0_wp
      do j = 1, self%nlat
         row = 0.0_wp
         do i = 1, self%nphi
            row = row + f(i, j)
         end do
         total = total + self%gauss_w(j)*row*dlon
      end do
   end function sht_grid_surface_integral

   subroutine sht_grid_synthesis(self, slm, sh)
      !! Spectral coefficients -> spatial field, shape (nphi, nlat).
      class(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: slm(:)   !! length nlm
      real(wp),        intent(out) :: sh(:,:)  !! (nphi, nlat)
      call SH_to_spat(self%cfg, slm, sh)
   end subroutine sht_grid_synthesis

   subroutine sht_grid_analysis(self, sh, slm)
      !! Spatial field -> spectral coefficients.
      !! NB: SHTns overwrites the input spatial array, hence intent(inout).
      class(sht_grid), intent(in)    :: self
      real(wp),        intent(inout) :: sh(:,:)  !! (nphi, nlat)
      complex(wp),     intent(out)   :: slm(:)   !! length nlm
      call spat_to_SH(self%cfg, sh, slm)
   end subroutine sht_grid_analysis

   subroutine sht_grid_sph_synthesis(self, slm, vth, vph)
      !! Spheroidal (gradient) vector synthesis: from the spectral coefficients slm
      !! of a scalar potential, return the surface-gradient field on the grid,
      !!   vth = ∂_θ(Σ slm Y_lm),   vph = (1/sinθ) ∂_φ(Σ slm Y_lm),
      !! i.e. the E_lm and F_lm tensor-harmonic building blocks (Martinec 2000 B11).
      !! Used to bootstrap the tensor-SH dyadic basis (fe_tensor_sh).
      class(sht_grid), intent(in)  :: self
      complex(wp),     intent(in)  :: slm(:)      !! length nlm
      real(wp),        intent(out) :: vth(:,:), vph(:,:)  !! (nphi, nlat)
      call SHsph_to_spat(self%cfg, slm, vth, vph)
   end subroutine sht_grid_sph_synthesis

   subroutine sht_grid_eval_point(self, f_lm, colat, lon, val)
      !! Evaluate a scalar field (spectral coefficients f_lm) at an ARBITRARY point
      !! (colat, lon) [rad] -- not restricted to a grid node. Uses SHTns'
      !! SHqst_to_point with zero spheroidal/toroidal parts, so the normalization
      !! matches analysis/synthesis exactly. Used to sample model fields along the
      !! Martinec-2018 benchmark profiles (circles of constant lon or lat).
      class(sht_grid), intent(in) :: self
      complex(wp),     intent(in) :: f_lm(:)    !! length nlm
      real(wp),        intent(in) :: colat, lon !! [rad]
      real(wp),        intent(out):: val
      complex(wp) :: q(self%nlm), s(self%nlm), t(self%nlm)
      real(wp)    :: vr(1), vt(1), vp(1)
      q = f_lm;  s = (0.0_wp, 0.0_wp);  t = (0.0_wp, 0.0_wp)
      call SHqst_to_point(self%cfg, q, s, t, cos(colat), lon, vr, vt, vp)
      val = vr(1)
   end subroutine sht_grid_eval_point

   subroutine sht_grid_eval_point_horiz(self, s_lm, colat, lon, vth, vph)
      !! Evaluate the HORIZONTAL field of a spheroidal scalar (coefficients s_lm)
      !! at an arbitrary point (colat, lon) [rad]: ∇₁(Σ s_lm Y_lm) =
      !! (∂_θ, (1/sinθ)∂_φ)(Σ s_lm Y_lm), via SHqst_to_point with q = t = 0 and
      !! spheroidal = s_lm. Returns vth = θ-component, vph = φ-component. Used for
      !! the horizontal-displacement (u_θ, u_φ) columns of the Martinec benchmark;
      !! s_lm = response%horizontal output (the per-degree V(a) coefficients).
      class(sht_grid), intent(in) :: self
      complex(wp),     intent(in) :: s_lm(:)    !! length nlm (spheroidal V(a))
      real(wp),        intent(in) :: colat, lon !! [rad]
      real(wp),        intent(out):: vth, vph   !! θ-, φ-components
      complex(wp) :: q(self%nlm), s(self%nlm), t(self%nlm)
      real(wp)    :: vr(1), vt(1), vp(1)
      q = (0.0_wp, 0.0_wp);  s = s_lm;  t = (0.0_wp, 0.0_wp)
      call SHqst_to_point(self%cfg, q, s, t, cos(colat), lon, vr, vt, vp)
      vth = vt(1);  vph = vp(1)
   end subroutine sht_grid_eval_point_horiz

   integer function sht_grid_lmidx(self, l, m) result(lm)
      !! 1-based index into the spectral array for harmonic (l, m).
      class(sht_grid), intent(in) :: self
      integer,         intent(in) :: l, m
      lm = shtns_lmidx(self%cfg, l, m)
   end function sht_grid_lmidx

end module fe_sht
