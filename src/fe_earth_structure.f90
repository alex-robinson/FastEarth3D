module fe_earth_structure
   !! Reference Earth structure: a radially layered, incompressible Maxwell Earth,
   !! plus an optional 3D (laterally varying) viscosity field.
   !!
   !! The layered description is the *physical* model (piecewise-constant density,
   !! shear modulus, and viscosity between interface radii) — e.g. the benchmark
   !! model M3-L70-V01. The radial finite-element mesh that the solver discretizes
   !! onto is built from this by fe_radial_fe.
   !!
   !! 3D-ready (project goal): `visc_3d`, when allocated, carries absolute
   !! log10-viscosity on the Gauss-Legendre spatial grid per radial node — how
   !! VILMA injects lateral heterogeneity (Albrecht et al. 2024). 1D runs leave it
   !! unallocated; the same solver path reduces to the spherically symmetric case.
   use fe_precision, only: wp
   use fe_constants, only: pi, grav_G
   use fe_params,    only: fe_param_class, MAX_LAYER
   implicit none
   private

   ! Rheology of a layer.
   integer, parameter, public :: RHEOL_ELASTIC = 0   !! elastic (viscosity -> inf)
   integer, parameter, public :: RHEOL_MAXWELL = 1   !! Maxwell viscoelastic
   integer, parameter, public :: RHEOL_FLUID   = 2   !! inviscid fluid (mu = 0)

   public :: earth_layer, earth_model, build_M3L70V01, build_earth

   type :: earth_layer
      !! A homogeneous spherical shell, r_bot <= r <= r_top.
      real(wp) :: r_bot = 0.0_wp   !! inner radius [m]
      real(wp) :: r_top = 0.0_wp   !! outer radius [m]
      real(wp) :: rho   = 0.0_wp   !! density [kg m^-3]
      real(wp) :: mu    = 0.0_wp   !! shear modulus / rigidity [Pa]
      real(wp) :: eta   = 0.0_wp   !! viscosity [Pa s] (huge for elastic, 0 fluid)
      integer  :: rheology = RHEOL_MAXWELL
   end type earth_layer

   type :: earth_model
      character(len=:), allocatable :: name
      real(wp) :: r_earth = 0.0_wp        !! surface radius [m]
      real(wp) :: r_core  = 0.0_wp        !! core-mantle boundary radius [m]
      type(earth_layer), allocatable :: layers(:)   !! surface-first (index 1 = top)
      ! Optional lateral viscosity: ABSOLUTE log10(η [Pa·s]) on the Gauss grid ×
      ! FE radial nodes, (nphi*nlat, nr). Populated by fe_read_visc_3d from a real
      ! lon-lat-r field (rung 6c). The node→element bridge (ve_response) takes the
      ! log10-mean of the two bracketing nodes and forms the perturbation against
      ! the element's radial reference η — storing the absolute field here avoids a
      ! per-node reference ambiguity at layer interfaces.
      real(wp), allocatable :: visc_3d(:,:)
   contains
      procedure :: n_layers     => earth_n_layers
      procedure :: rho_at       => earth_rho_at
      procedure :: mu_at        => earth_mu_at
      procedure :: eta_at       => earth_eta_at
      procedure :: mass_below   => earth_mass_below
      procedure :: total_mass   => earth_total_mass
      procedure :: gravity_at   => earth_gravity_at
      procedure :: moi          => earth_moi
      procedure :: is_3d        => earth_is_3d
   end type earth_model

contains

   ! --- Construction ----------------------------------------------------------

   function build_earth(p) result(em)
      !! Build the earth model selected by p%earth: a named built-in, or "custom"
      !! assembled from the per-layer arrays (surface-first) in the parameter record.
      type(fe_param_class), intent(in) :: p
      type(earth_model) :: em
      integer :: k, n

      select case (trim(p%earth))
      case ("M3-L70-V01")
         em = build_M3L70V01()
      case ("custom")
         n = p%n_layer
         if (n < 1 .or. n > MAX_LAYER) &
            error stop 'build_earth: n_layer out of range for a custom earth model'
         em%name    = "custom"
         em%r_earth = p%r_earth
         em%r_core  = p%r_core
         allocate(em%layers(n))
         do k = 1, n
            em%layers(k) = earth_layer(p%r_bot(k), p%r_top(k), p%rho(k), &
                                       p%mu(k), p%eta(k), p%rheology(k))
         end do
      case default
         error stop 'build_earth: unknown earth model "'//trim(p%earth)//'"'
      end select
   end function build_earth

   function build_M3L70V01() result(em)
      !! Benchmark Earth model M3-L70-V01 (Spada et al. 2011, GJI 185:106,
      !! Table 3; values cross-checked against ALMA3 MODELS/M3L70V01.dat).
      !! Incompressible, self-gravitating, Maxwell. 70 km elastic lithosphere,
      !! 3 Maxwell mantle layers, inviscid fluid core.
      type(earth_model) :: em
      real(wp), parameter :: km = 1.0e3_wp

      em%name    = "M3-L70-V01"
      em%r_earth = 6371.0_wp*km
      em%r_core  = 3480.0_wp*km
      allocate(em%layers(5))
      !                       r_bot[m]      r_top[m]      rho      mu[Pa]      eta[Pa s]      rheology
      em%layers(1) = earth_layer(6301.0_wp*km, 6371.0_wp*km, 3037.0_wp,  5.0605e10_wp, huge(1.0_wp), RHEOL_ELASTIC) ! lithosphere
      em%layers(2) = earth_layer(5951.0_wp*km, 6301.0_wp*km, 3438.0_wp,  7.0363e10_wp, 1.0e21_wp,    RHEOL_MAXWELL) ! upper mantle
      em%layers(3) = earth_layer(5701.0_wp*km, 5951.0_wp*km, 3871.0_wp,  1.0549e11_wp, 1.0e21_wp,    RHEOL_MAXWELL) ! transition
      em%layers(4) = earth_layer(3480.0_wp*km, 5701.0_wp*km, 4978.0_wp,  2.2834e11_wp, 2.0e21_wp,    RHEOL_MAXWELL) ! lower mantle
      em%layers(5) = earth_layer(   0.0_wp,    3480.0_wp*km, 10750.0_wp, 0.0_wp,       0.0_wp,       RHEOL_FLUID)   ! core
   end function build_M3L70V01

   ! --- Queries ---------------------------------------------------------------

   integer function earth_n_layers(self) result(n)
      class(earth_model), intent(in) :: self
      n = 0
      if (allocated(self%layers)) n = size(self%layers)
   end function earth_n_layers

   integer function layer_at(self, r) result(k)
      !! Index of the layer containing radius r (-1 if outside the model).
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: i
      k = -1
      do i = 1, self%n_layers()
         if (r >= self%layers(i)%r_bot .and. r <= self%layers(i)%r_top) then
            k = i
            return
         end if
      end do
   end function layer_at

   real(wp) function earth_rho_at(self, r) result(val)
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%rho
   end function earth_rho_at

   real(wp) function earth_mu_at(self, r) result(val)
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%mu
   end function earth_mu_at

   real(wp) function earth_eta_at(self, r) result(val)
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%eta
   end function earth_eta_at

   real(wp) function earth_mass_below(self, r) result(m)
      !! Mass enclosed within radius r [kg], integrating the piecewise-constant
      !! density profile (partial shells handled).
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      real(wp) :: lo, hi
      integer  :: i
      m = 0.0_wp
      do i = 1, self%n_layers()
         lo = self%layers(i)%r_bot
         hi = min(self%layers(i)%r_top, r)
         if (hi > lo) m = m + (4.0_wp/3.0_wp)*pi*self%layers(i)%rho*(hi**3 - lo**3)
      end do
   end function earth_mass_below

   real(wp) function earth_total_mass(self) result(m)
      class(earth_model), intent(in) :: self
      m = self%mass_below(self%r_earth)
   end function earth_total_mass

   real(wp) function earth_gravity_at(self, r) result(g)
      !! Unperturbed gravity at radius r [m s^-2]: G * M(<r) / r^2.
      class(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      g = 0.0_wp
      if (r > 0.0_wp) g = grav_G*self%mass_below(r)/(r*r)
   end function earth_gravity_at

   real(wp) function earth_moi(self) result(inertia)
      !! Moment of inertia about a diameter [kg m^2] for the spherically
      !! symmetric profile: I = (8 pi / 15) * sum rho * (r_top^5 - r_bot^5).
      class(earth_model), intent(in) :: self
      integer :: i
      inertia = 0.0_wp
      do i = 1, self%n_layers()
         inertia = inertia + (8.0_wp*pi/15.0_wp)*self%layers(i)%rho* &
                   (self%layers(i)%r_top**5 - self%layers(i)%r_bot**5)
      end do
   end function earth_moi

   logical function earth_is_3d(self) result(yes)
      class(earth_model), intent(in) :: self
      yes = allocated(self%visc_3d)
   end function earth_is_3d

end module fe_earth_structure
