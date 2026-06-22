module fe_lis
   !! Thin Fortran wrapper around LIS (Library of Iterative Solvers), isolating
   !! every LIS call (and the lisf.h preprocessor interface) from the physics
   !! modules. fe_radial_fe hands this module a sparse system in coordinate
   !! (COO) form and gets back the solution — no LIS types leak outward.
   !!
   !! The per-degree saddle-point operator is band-diagonal, indefinite and
   !! non-symmetric, so the default solver here is a nonsymmetric Krylov method
   !! (BiCGSTAB/GMRES) with an ILU preconditioner. Callers pass the LIS option
   !! string so the solver/preconditioner can be tuned per problem.
#include "lisf.h"
   use fe_precision, only: wp
   implicit none
   private

   public :: fe_lis_initialize, fe_lis_finalize, fe_lis_solve_coo

   logical, save :: initialized = .false.

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

   subroutine fe_lis_solve_coo(n, rows, cols, vals, b, x, options, &
                               iters, resid, info)
      !! Solve A x = b for a sparse A given in COO form (1-based rows/cols).
      !! Duplicate (i,j) entries are summed by LIS during assembly, so the COO
      !! list may carry the raw element-by-element contributions.
      integer,          intent(in)  :: n             !! system size
      integer,          intent(in)  :: rows(:)        !! COO row indices (1-based)
      integer,          intent(in)  :: cols(:)        !! COO column indices (1-based)
      real(wp),         intent(in)  :: vals(:)        !! COO values
      real(wp),         intent(in)  :: b(:)           !! right-hand side (n)
      real(wp),         intent(out) :: x(:)           !! solution (n)
      character(len=*), intent(in)  :: options        !! LIS option string
      integer,          intent(out) :: iters          !! iterations taken
      real(wp),         intent(out) :: resid          !! final relative residual
      integer,          intent(out) :: info           !! LIS error code (0 = OK)

      LIS_MATRIX  :: A
      LIS_VECTOR  :: bb, xx
      LIS_SOLVER  :: solver
      LIS_INTEGER :: ierr, i, nnz, lis_iter
      LIS_REAL    :: lis_resid
      real(wp)    :: val

      call fe_lis_initialize()
      nnz = size(vals)

      ! --- matrix: COO -> CSR --------------------------------------------------
      call lis_matrix_create(LIS_COMM_WORLD, A, ierr)
      call lis_matrix_set_size(A, 0, n, ierr)
      do i = 1, nnz
         call lis_matrix_set_value(LIS_ADD_VALUE, rows(i), cols(i), vals(i), A, ierr)
      end do
      call lis_matrix_set_type(A, LIS_MATRIX_CSR, ierr)
      call lis_matrix_assemble(A, ierr)

      ! --- vectors -------------------------------------------------------------
      call lis_vector_duplicate(A, bb, ierr)
      call lis_vector_duplicate(A, xx, ierr)
      do i = 1, n
         call lis_vector_set_value(LIS_INS_VALUE, i, b(i), bb, ierr)
      end do

      ! --- solve ---------------------------------------------------------------
      call lis_solver_create(solver, ierr)
      call lis_solver_set_option(trim(options), solver, ierr)
      call lis_solve(A, bb, xx, solver, ierr)
      info = ierr

      call lis_solver_get_iter(solver, lis_iter, ierr)
      call lis_solver_get_residualnorm(solver, lis_resid, ierr)
      iters = int(lis_iter)
      resid = real(lis_resid, wp)

      do i = 1, n
         call lis_vector_get_value(xx, i, val, ierr)
         x(i) = val
      end do

      call lis_solver_destroy(solver, ierr)
      call lis_vector_destroy(bb, ierr)
      call lis_vector_destroy(xx, ierr)
      call lis_matrix_destroy(A, ierr)
   end subroutine fe_lis_solve_coo

end module fe_lis
