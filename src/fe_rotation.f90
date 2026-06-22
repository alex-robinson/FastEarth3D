module fe_rotation
   !! Rotational feedback / true polar wander (Martinec & Hagedoorn 2014; in
   !! VILMA, so included here). Surface loading and deformation perturb the
   !! inertia tensor; a linearized Liouville equation maps that to polar motion,
   !! which feeds back as a centrifugal-potential perturbation forcing the
   !! sea-level equation.
   !!
   !! Key conventions (doc/design.md, section D4): the rotational forcing is
   !! purely degree 2 (NOT degree 1), it uses *tidal* Love numbers, and the fluid
   !! Love number k_f must be fixed from the observed flattening (Mitrovica et
   !! al. 2005) to avoid the spurious lithosphere-thickness sensitivity.
   !!
   !! STATUS: scaffold — interface only.
   use fe_precision, only: wp
   use fe_constants, only: omega_earth
   use fe_sht,       only: sht_grid
   implicit none
   private

   public :: rotation_state

   type :: rotation_state
      logical  :: enabled = .false.    !! off until validated against benchmark
      real(wp) :: m(3) = 0.0_wp        !! polar-motion / spin-rate perturbations
   contains
      procedure :: update => rotation_update
   end type rotation_state

contains

   subroutine rotation_update(self, inertia_perturbation, centrifugal_lm)
      !! Liouville update: inertia-tensor change -> polar motion -> centrifugal
      !! potential perturbation (degree-2 spectral forcing). STATUS: stub.
      class(rotation_state), intent(inout) :: self
      real(wp),              intent(in)    :: inertia_perturbation(3,3)
      complex(wp),           intent(out)   :: centrifugal_lm(:)
      centrifugal_lm = (0.0_wp, 0.0_wp)
      if (.not. self%enabled) return
      ! TODO: solve linearized Liouville equation for self%m, then build the
      ! degree-2 centrifugal potential from m using tidal Love numbers.
   end subroutine rotation_update

end module fe_rotation
