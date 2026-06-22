module fe_radial_fe
   !! Radial finite-element discretization of the viscoelastic field equations.
   !!
   !! For the spherically symmetric background each spherical-harmonic degree `l`
   !! decouples into a small banded linear system in radius (piecewise-linear
   !! "tent" FE, Galerkin weak form; Martinec 2000). These per-degree stiffness
   !! operators are assembled once and reused every time step — the explicit
   !! memory-stress scheme means the angular orders never couple in the solve.
   !!
   !! This module provides:
   !!   - radial_mesh:      the P1 mesh over the mantle shell [r_core, r_earth]
   !!                       (the inviscid core is a CMB boundary condition, not
   !!                       meshed — Wu & Peltier 1982 / VEGA, Martinec et al.
   !!                       2018). VERIFIED and implemented.
   !!   - radial_operator:  the per-degree banded system. STATUS: interface only;
   !!                       the incompressible weak form awaits confirmation of
   !!                       Martinec (2000)'s mixed element pair (see doc/design.md).
   use fe_precision, only: wp
   use fe_earth_structure, only: earth_model
   implicit none
   private

   public :: radial_mesh, radial_operator

   ! Default radial element-size targets [m], after VEGA (Martinec et al. 2018):
   ! 5 km in the lithosphere, 10 km in the upper mantle, 40 km in the lower
   ! mantle. Selected by depth thresholds below.
   real(wp), parameter :: DR_LITHO = 5.0e3_wp
   real(wp), parameter :: DR_UPPER = 10.0e3_wp
   real(wp), parameter :: DR_LOWER = 40.0e3_wp
   real(wp), parameter :: DEPTH_LITHO =  70.0e3_wp   !! lithosphere base depth
   real(wp), parameter :: DEPTH_UPPER = 670.0e3_wp   !! upper/lower mantle divide

   type :: radial_mesh
      !! Piecewise-linear radial mesh over the meshed (solid) shell.
      integer :: nr = 0                       !! number of nodes
      integer :: ne = 0                       !! number of elements (= nr - 1)
      real(wp), allocatable :: r(:)           !! node radii [m], strictly ascending (nr)
      integer,  allocatable :: elem_layer(:)  !! source earth layer per element (ne)
   contains
      procedure :: build => radial_mesh_build
   end type radial_mesh

   type :: radial_operator
      !! Banded radial system for one degree l (factored, ready to solve).
      integer :: l  = -1
      integer :: nr = 0
      ! TODO: banded matrix storage + factorization for the spheroidal degrees of
      ! freedom (U, V, potential F, pressure Π) per node, once B4 is confirmed.
   contains
      procedure :: assemble => radial_operator_assemble
      procedure :: solve    => radial_operator_solve
   end type radial_operator

contains

   ! --- Mesh ------------------------------------------------------------------

   pure function dr_target(depth) result(dr)
      !! Target radial element size at a given depth [m] (VEGA spacing).
      real(wp), intent(in) :: depth
      real(wp) :: dr
      if (depth <= DEPTH_LITHO) then
         dr = DR_LITHO
      else if (depth <= DEPTH_UPPER) then
         dr = DR_UPPER
      else
         dr = DR_LOWER
      end if
   end function dr_target

   subroutine radial_mesh_build(self, earth)
      !! Build the radial mesh over the SOLID shell [r_core, r_earth]. The core
      !! (a fluid layer reaching r=0) is excluded — it enters as a CMB boundary
      !! condition. Each solid layer is subdivided into uniform elements no larger
      !! than the depth-dependent target, with nodes pinned to every material
      !! interface so no element straddles a density/rigidity jump.
      class(radial_mesh), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      real(wp), allocatable :: r(:)
      integer,  allocatable :: lay(:)
      integer :: i, k, ne_layer, ntot, off
      real(wp) :: r0, r1, depth_mid, dr, h

      ! Count elements per solid layer first (a layer is "solid" if its top is
      ! above the CMB; this excludes the fluid core layer).
      ntot = 0
      do i = 1, earth%n_layers()
         r0 = earth%layers(i)%r_bot
         r1 = earth%layers(i)%r_top
         if (r1 <= earth%r_core) cycle            ! core (or below CMB): skip
         r0 = max(r0, earth%r_core)               ! clip to the meshed shell
         depth_mid = earth%r_earth - 0.5_wp*(r0 + r1)
         dr = dr_target(depth_mid)
         ne_layer = max(1, ceiling((r1 - r0)/dr))
         ntot = ntot + ne_layer
      end do

      self%ne = ntot
      self%nr = ntot + 1
      allocate(r(self%nr), lay(self%ne))

      ! Lay down nodes layer by layer, sharing interface nodes between layers.
      r(1) = earth%r_core
      off  = 1            ! index of the last node placed
      do i = earth%n_layers(), 1, -1   ! innermost solid layer first (ascending r)
         r0 = earth%layers(i)%r_bot
         r1 = earth%layers(i)%r_top
         if (r1 <= earth%r_core) cycle
         r0 = max(r0, earth%r_core)
         depth_mid = earth%r_earth - 0.5_wp*(r0 + r1)
         dr = dr_target(depth_mid)
         ne_layer = max(1, ceiling((r1 - r0)/dr))
         h = (r1 - r0)/real(ne_layer, wp)
         do k = 1, ne_layer
            r(off + k)   = r0 + real(k, wp)*h
            lay(off + k - 1) = i
         end do
         off = off + ne_layer
      end do
      r(self%nr) = earth%r_earth   ! pin the surface exactly

      call move_alloc(r,   self%r)
      call move_alloc(lay, self%elem_layer)
   end subroutine radial_mesh_build

   ! --- Operator (stub) -------------------------------------------------------

   subroutine radial_operator_assemble(self, earth, mesh, l)
      !! Assemble + factor the banded radial system for degree l on `mesh`.
      !! STATUS: stub — awaits the confirmed Martinec (2000) weak form (B4).
      class(radial_operator), intent(inout) :: self
      type(earth_model),      intent(in)    :: earth
      type(radial_mesh),      intent(in)    :: mesh
      integer,                intent(in)    :: l
      self%l  = l
      self%nr = mesh%nr
      ! TODO: build element matrices by sampling earth%{rho_at,mu_at,gravity_at}
      ! on the mesh, couple self-gravity (Poisson), apply the surface-load and
      ! CMB fluid boundary conditions, enforce incompressibility, factor.
      if (earth%n_layers() == 0) return
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
