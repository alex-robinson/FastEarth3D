module fe_viscoelastic
   !! Time-domain viscoelastic relaxation: the heart of the method.
   !!
   !! Incompressible Maxwell rheology integrated explicitly in the time domain
   !! (Martinec 2000). Each step the stress splits into an instantaneous elastic
   !! part plus a viscous *memory* stress carried from the previous step and
   !! relaxed by a factor set by Δt / Maxwell-time. The memory stress is the only
   !! place lateral viscosity enters, evaluated on the Gauss-Legendre grid and
   !! transformed to spectral each step (pseudo-spectral). The explicit scheme is
   !! conditionally stable, hence a viscosity floor (~1e19 Pa s).
   !!
   !! STATUS: scaffold — state container + time-loop interface only.
   use fe_precision, only: wp
   use fe_radial_fe, only: radial_operator
   use fe_sht,       only: sht_grid
   implicit none
   private

   public :: viscoelastic_state

   type :: viscoelastic_state
      real(wp) :: time = 0.0_wp     !! current model time [s]
      real(wp) :: dt   = 0.0_wp     !! time step [s]
      ! TODO: spectral displacement, gravitational potential, and the radial
      ! memory-stress field carried between steps.
   contains
      procedure :: step => viscoelastic_step
   end type viscoelastic_state

contains

   subroutine viscoelastic_step(self, ops, sht, load_lm)
      !! Advance one explicit time step under the surface load `load_lm`.
      !! STATUS: stub — wires the intended dependencies, no physics yet.
      class(viscoelastic_state), intent(inout) :: self
      type(radial_operator),     intent(in)    :: ops(:)   !! one per degree
      type(sht_grid),            intent(in)     :: sht
      complex(wp),               intent(in)     :: load_lm(:)
      self%time = self%time + self%dt
      ! TODO:
      !   1. evaluate memory stress on the spatial grid (lateral viscosity here)
      !   2. forward transform to spectral (sht%analysis)
      !   3. per-degree banded solve (ops(l)%solve) for the elastic-like update
      !   4. update displacement / potential / memory stress
   end subroutine viscoelastic_step

end module fe_viscoelastic
