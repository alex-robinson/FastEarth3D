module fe_response
   !! Surface-load response operator: the abstraction the sea-level equation
   !! (fe_sle) is built on. Given a spectral surface mass-density load σ_lm
   !! [kg m^-2], it returns the two fields the SLE needs,
   !!
   !!     u_lm  — radial displacement of the solid surface  [m]
   !!     n_lm  — geoid / sea-surface-equipotential height   [m]
   !!
   !! per spherical-harmonic coefficient. The SLE depends only on this interface,
   !! so the elastic and (later) viscoelastic earth responses are swappable.
   !!
   !! Geoid mapping. The per-degree solve returns U(a) and F(a) = φ₁(a), the
   !! surface coefficients of radial displacement and the perturbed gravitational
   !! potential. Martinec's φ₁ carries the load's own direct potential with the
   !! sign OPPOSITE to φ^L (φ₁ → −φ^L for a rigid sphere, k = −F/φ^L − 1; see
   !! fe_radial_fe%loading_love). The geopotential perturbation is therefore −F,
   !! and Bruns' formula gives the geoid height
   !!
   !!     N(a) = −F(a)/g .
   !!
   !! This uses only U and F — NOT the horizontal Love number l, whose sign /
   !! normalization is still being calibrated — so the SLE is not blocked by
   !! that open item. Both U and the −F/g geoid are pinned by the validated
   !! rigid (U→0, 1+k→1) and fluid (U→−(2j+1)/3·φ^L/g, 1+k→0) limits.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_mesh, radial_operator
   use fe_sht,             only: sht_grid
   implicit none
   private

   public :: response_operator, elastic_response, null_response

   type, abstract :: response_operator
      !! Maps a spectral surface load to surface displacement + geoid.
   contains
      procedure(apply_if), deferred :: apply
   end type response_operator

   abstract interface
      subroutine apply_if(self, sht, sigma_lm, u_lm, n_lm)
         import :: response_operator, sht_grid, wp
         class(response_operator), intent(inout) :: self
         type(sht_grid),           intent(in)    :: sht
         complex(wp),              intent(in)    :: sigma_lm(:)  !! load [kg m^-2]
         complex(wp),              intent(out)   :: u_lm(:)      !! uplift  [m]
         complex(wp),              intent(out)   :: n_lm(:)      !! geoid   [m]
      end subroutine apply_if
   end interface

   type, extends(response_operator) :: elastic_response
      !! Time-independent (elastic) response: per-degree surface response to a
      !! unit load, precomputed once. Linear and order-independent, so a single
      !! gain per degree multiplies every coefficient of that degree.
      integer  :: lmax = 0
      real(wp) :: g    = 0.0_wp          !! surface gravity [m s^-2]
      real(wp) :: a    = 0.0_wp          !! surface radius  [m]
      real(wp), allocatable :: ugain(:)  !! (0:lmax) U(a) per unit σ_l  [m / (kg m^-2)]
      real(wp), allocatable :: ngain(:)  !! (0:lmax) N(a)=−F(a)/g per unit σ_l
   contains
      procedure :: init    => elastic_response_init
      procedure :: apply   => elastic_response_apply
      procedure :: destroy => elastic_response_destroy
   end type elastic_response

   type, extends(response_operator) :: null_response
      !! Rigid, non-self-gravitating Earth: u ≡ 0 and N ≡ 0. The SLE then
      !! reduces to a uniform (eustatic/barystatic) ocean response, which is the
      !! textbook limit used to check mass conservation and the uniform term.
   contains
      procedure :: apply => null_response_apply
   end type null_response

contains

   subroutine null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      class(null_response), intent(inout) :: self
      type(sht_grid),       intent(in)    :: sht
      complex(wp),          intent(in)    :: sigma_lm(:)
      complex(wp),          intent(out)   :: u_lm(:)
      complex(wp),          intent(out)   :: n_lm(:)
      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
   end subroutine null_response_apply

   subroutine elastic_response_init(self, earth, lmax)
      !! Precompute the per-degree elastic surface gains for degrees 0..lmax.
      !!
      !!   l = 0 : incompressibility (Div u = 0) forbids degree-0 radial
      !!           deformation, so U(0)=0; the geoid feels only the load's own
      !!           monopole potential, N(0) = φ^L_0/g = 4πGa/g per unit σ.
      !!   l ≥ 1 : assemble the per-degree saddle-point operator, solve a unit
      !!           surface load, store U(a) and N(a) = −F(a)/g.
      class(elastic_response), intent(inout) :: self
      type(earth_model),       intent(in)    :: earth
      integer,                 intent(in)    :: lmax
      type(radial_mesh)     :: mesh
      type(radial_operator) :: op
      integer  :: l
      real(wp) :: ua, va, fa

      call self%destroy()
      self%lmax = lmax
      self%a    = earth%r_earth
      self%g    = earth%gravity_at(earth%r_earth)
      allocate(self%ugain(0:lmax), self%ngain(0:lmax))

      ! degree 0: no deformation, pure monopole geoid offset
      self%ugain(0) = 0.0_wp
      self%ngain(0) = 4.0_wp*pi*grav_G*self%a / self%g

      call mesh%build(earth)
      do l = 1, lmax
         call op%assemble(earth, mesh, l)
         call op%solve(1.0_wp, ua, va, fa)     ! unit surface load coefficient
         self%ugain(l) = ua
         self%ngain(l) = -fa / self%g
         call op%destroy()
      end do
   end subroutine elastic_response_init

   subroutine elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Spectral multiply: u_lm = ugain(l)·σ_lm, n_lm = ngain(l)·σ_lm. Degrees
      !! above the precomputed lmax are zeroed.
      class(elastic_response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: u_lm(:)
      complex(wp),             intent(out)   :: n_lm(:)
      integer :: l, m, lm, lcap

      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht%lmidx(l, m)
            u_lm(lm) = self%ugain(l) * sigma_lm(lm)
            n_lm(lm) = self%ngain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_apply

   subroutine elastic_response_destroy(self)
      class(elastic_response), intent(inout) :: self
      if (allocated(self%ugain)) deallocate(self%ugain)
      if (allocated(self%ngain)) deallocate(self%ngain)
      self%lmax = 0
   end subroutine elastic_response_destroy

end module fe_response
