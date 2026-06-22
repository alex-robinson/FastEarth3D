module fe_radial_fe
   !! Radial finite-element discretization of the viscoelastic field equations.
   !!
   !! For the spherically symmetric background each spherical-harmonic degree `l`
   !! decouples into a small banded linear system in radius (piecewise-linear
   !! FE, Galerkin weak form; Martinec 2000). These per-degree "stiffness"
   !! operators are assembled once and reused every time step — the explicit
   !! memory-stress scheme means the angular orders never couple in the solve,
   !! which is what makes the method fast and trivially parallel over (l,m).
   !!
   !! STATUS: scaffold — interface only.
   use fe_precision, only: wp
   use fe_earth_structure, only: earth_model
   implicit none
   private

   public :: radial_operator

   type :: radial_operator
      !! Banded radial system for one degree l (factored, ready to solve).
      integer :: l  = -1
      integer :: nr = 0
      ! TODO: banded matrix storage + factorization for the (u, v, potential)
      ! degrees of freedom per node.
   contains
      procedure :: assemble => radial_operator_assemble
      procedure :: solve    => radial_operator_solve
   end type radial_operator

contains

   subroutine radial_operator_assemble(self, earth, l)
      !! Assemble + factor the banded radial system for degree l. STATUS: stub.
      class(radial_operator), intent(inout) :: self
      type(earth_model),      intent(in)    :: earth
      integer,                intent(in)    :: l
      self%l  = l
      self%nr = earth%nr
      ! TODO: build element matrices from earth%{rho,mu}, apply self-gravity
      ! coupling and surface/CMB boundary conditions, factor.
   end subroutine radial_operator_assemble

   subroutine radial_operator_solve(self, rhs, sol)
      !! Solve the assembled system for one (l,m) right-hand side. STATUS: stub.
      class(radial_operator), intent(in)  :: self
      complex(wp),            intent(in)  :: rhs(:)
      complex(wp),            intent(out) :: sol(:)
      sol = (0.0_wp, 0.0_wp)
      ! TODO: banded back-substitution using the stored factorization.
   end subroutine radial_operator_solve

end module fe_radial_fe
