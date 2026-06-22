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
   contains
      procedure :: init      => sht_grid_init
      procedure :: destroy   => sht_grid_destroy
      procedure :: synthesis => sht_grid_synthesis   !! spectral -> spatial
      procedure :: analysis  => sht_grid_analysis    !! spatial  -> spectral
      procedure :: lmidx     => sht_grid_lmidx       !! (l,m) -> coefficient index
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
   end subroutine sht_grid_init

   subroutine sht_grid_destroy(self)
      !! Release the SHTns configuration.
      class(sht_grid), intent(inout) :: self
      if (c_associated(self%cfg)) call shtns_destroy(self%cfg)
      self%cfg = c_null_ptr
      self%lmax = 0; self%mmax = 0; self%nlat = 0; self%nphi = 0; self%nlm = 0
   end subroutine sht_grid_destroy

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

   integer function sht_grid_lmidx(self, l, m) result(lm)
      !! 1-based index into the spectral array for harmonic (l, m).
      class(sht_grid), intent(in) :: self
      integer,         intent(in) :: l, m
      lm = shtns_lmidx(self%cfg, l, m)
   end function sht_grid_lmidx

end module fe_sht
