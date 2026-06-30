! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Solvers of linear systems.
module regrid_solvers

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

public :: solve_linear_system, linear_solver, solve_tridiagonal_system, solve_diag_dominant_tridiag


  interface
module subroutine solve_linear_system( A, R, X, N, answer_date )
  integer,              intent(in)    :: N  !< The size of the system
  real, dimension(N,N), intent(inout) :: A  !< The matrix being inverted in arbitrary units [A] on
                                            !! input, but internally modified to become nondimensional
                                            !! during the solver.
  real, dimension(N),   intent(inout) :: R  !< system right-hand side in arbitrary units [A B] on
                                            !! input, but internally modified to have units of [B]
                                            !! during the solver
  real, dimension(N),   intent(inout) :: X  !< solution vector in arbitrary units [B]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use
  ! Local variables

end subroutine solve_linear_system
module subroutine linear_solver( N, A, R, X )
  integer,              intent(in)    :: N  !< The size of the system
  real, dimension(N,N), intent(inout) :: A  !< The matrix being inverted in arbitrary units [A] on
                                            !! input, but internally modified to become nondimensional
                                            !! during the solver.
  real, dimension(N),   intent(inout) :: R  !< system right-hand side in [A B] on input, but internally
                                            !! modified to have units of [B] during the solver
  real, dimension(N),   intent(inout) :: X  !< solution vector [B]

  ! Local variables

  ! Loop on rows to transform the problem into multiplication by an upper-right matrix.
end subroutine linear_solver
module subroutine solve_tridiagonal_system( Al, Ad, Au, R, X, N, answer_date )
  integer,            intent(in)  :: N   !< The size of the system
  real, dimension(N), intent(in)  :: Ad  !< Matrix center diagonal in arbitrary units [A]
  real, dimension(N), intent(in)  :: Al  !< Matrix lower diagonal [A]
  real, dimension(N), intent(in)  :: Au  !< Matrix upper diagonal [A]
  real, dimension(N), intent(in)  :: R   !< system right-hand side in arbitrary units [A B]
  real, dimension(N), intent(out) :: X   !< solution vector in arbitrary units [B]
  integer,  optional, intent(in)  :: answer_date  !< The vintage of the expressions to use
  ! Local variables

end subroutine solve_tridiagonal_system
module subroutine solve_diag_dominant_tridiag( Al, Ac, Au, R, X, N )
  integer,            intent(in)  :: N   !< The size of the system
  real, dimension(N), intent(in)  :: Ac  !< Matrix center diagonal offset from Al + Au in arbitrary units [A]
  real, dimension(N), intent(in)  :: Al  !< Matrix lower diagonal [A]
  real, dimension(N), intent(in)  :: Au  !< Matrix upper diagonal [A]
  real, dimension(N), intent(in)  :: R   !< system right-hand side in arbitrary units [A B]
  real, dimension(N), intent(out) :: X   !< solution vector in arbitrary units [B]
  ! Local variables

  ! Factorization and forward sweep, in a form that will never give a division by a
  ! zero pivot for positive definite Ac, Al, and Au.
end subroutine solve_diag_dominant_tridiag
  end interface

end module regrid_solvers
