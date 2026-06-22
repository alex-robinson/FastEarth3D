module fe_lis
   !! Thin Fortran wrapper around LIS (Library of Iterative Solvers), isolating
   !! every LIS call (and the lisf.h preprocessor interface) from the physics
   !! modules. fe_radial_fe hands this module a sparse system in coordinate
   !! (COO) form and gets back the solution — no LIS types leak outward.
   !!
   !! The per-degree saddle-point operator is band-diagonal, indefinite and
   !! non-symmetric, so the solver is a nonsymmetric Krylov method (GMRES) with
   !! an ILU preconditioner. The operator is FIXED across all loads, orders and
   !! time steps of a degree, so `fe_lis_system` builds the matrix and factors
   !! the ILU preconditioner ONCE and reuses them for every right-hand side
   !! (`lis_solve_kernel`) — essential for the viscoelastic time stepper, which
   !! solves the same operator thousands of times.
#include "lisf.h"
   use fe_precision, only: wp
   implicit none
   private

   public :: fe_lis_initialize, fe_lis_finalize
   public :: fe_lis_system

   logical, save :: initialized = .false.

   type :: fe_lis_system
      !! A built + preconditioned LIS system, ready for repeated RHS solves.
      private
      LIS_MATRIX  :: A
      LIS_VECTOR  :: bb, xx
      LIS_SOLVER  :: solver
      LIS_PRECON  :: precon
      integer     :: n = 0
      logical     :: built = .false.
   contains
      procedure :: build   => lis_system_build
      procedure :: solve   => lis_system_solve
      procedure :: destroy => lis_system_destroy
   end type fe_lis_system

contains

   subroutine fe_lis_initialize()
      !! Initialize the LIS runtime once (idempotent). Must precede any solve.
      LIS_INTEGER :: ierr
      if (initialized) return
      call lis_initialize(ierr)
      initialized = .true.
   end subroutine fe_lis_initialize

   subroutine fe_lis_finalize()
      !! Shut the LIS runtime down (idempotent).
      LIS_INTEGER :: ierr
      if (.not. initialized) return
      call lis_finalize(ierr)
      initialized = .false.
   end subroutine fe_lis_finalize

   subroutine lis_system_build(self, n, rows, cols, vals, options)
      !! Build the CSR matrix from COO (1-based, duplicates summed) and factor
      !! the preconditioner once. `options` is the LIS solver/preconditioner
      !! string (e.g. '-i gmres -p ilu ...'); the matrix and precon are then
      !! reused for every solve().
      class(fe_lis_system), intent(inout) :: self
      integer,              intent(in)    :: n
      integer,              intent(in)    :: rows(:), cols(:)
      real(wp),             intent(in)    :: vals(:)
      character(len=*),     intent(in)    :: options
      LIS_INTEGER :: ierr, i, nnz

      call fe_lis_initialize()
      call self%destroy()
      self%n = n
      nnz = size(vals)

      call lis_matrix_create(LIS_COMM_WORLD, self%A, ierr)
      call lis_matrix_set_size(self%A, 0, n, ierr)
      do i = 1, nnz
         call lis_matrix_set_value(LIS_ADD_VALUE, rows(i), cols(i), vals(i), self%A, ierr)
      end do
      call lis_matrix_set_type(self%A, LIS_MATRIX_CSR, ierr)
      call lis_matrix_assemble(self%A, ierr)

      call lis_vector_duplicate(self%A, self%bb, ierr)
      call lis_vector_duplicate(self%A, self%xx, ierr)

      ! Build the solver + preconditioner ONCE (test8f.F90 reuse pattern):
      ! set options and the matrix, then create/update the preconditioner.
      call lis_solver_create(self%solver, ierr)
      call lis_solver_set_option(trim(options), self%solver, ierr)
      call lis_solver_set_matrix(self%A, self%solver, ierr)
      call lis_precon_psd_create(self%solver, self%precon, ierr)
      call lis_precon_psd_update(self%solver, self%precon, ierr)
      self%built = .true.
   end subroutine lis_system_build

   subroutine lis_system_solve(self, b, x, iters, resid, info)
      !! Solve A x = b for a new RHS, reusing the stored matrix + preconditioner.
      class(fe_lis_system), intent(in)  :: self
      real(wp),             intent(in)  :: b(:)
      real(wp),             intent(out) :: x(:)
      integer,  optional,   intent(out) :: iters, info
      real(wp), optional,   intent(out) :: resid
      LIS_INTEGER :: ierr, i, lis_iter
      LIS_REAL    :: lis_resid
      real(wp)    :: val
      do i = 1, self%n
         call lis_vector_set_value(LIS_INS_VALUE, i, b(i), self%bb, ierr)
      end do
      call lis_solve_kernel(self%A, self%bb, self%xx, self%solver, self%precon, ierr)
      if (present(info)) info = int(ierr)
      call lis_solver_get_iter(self%solver, lis_iter, ierr)
      call lis_solver_get_residualnorm(self%solver, lis_resid, ierr)
      if (present(iters)) iters = int(lis_iter)
      if (present(resid)) resid = real(lis_resid, wp)
      do i = 1, self%n
         call lis_vector_get_value(self%xx, i, val, ierr)
         x(i) = val
      end do
   end subroutine lis_system_solve

   subroutine lis_system_destroy(self)
      class(fe_lis_system), intent(inout) :: self
      LIS_INTEGER :: ierr
      if (.not. self%built) return
      call lis_precon_destroy(self%precon, ierr)
      call lis_solver_destroy(self%solver, ierr)
      call lis_vector_destroy(self%bb, ierr)
      call lis_vector_destroy(self%xx, ierr)
      call lis_matrix_destroy(self%A, ierr)
      self%built = .false.
      self%n = 0
   end subroutine lis_system_destroy

end module fe_lis
