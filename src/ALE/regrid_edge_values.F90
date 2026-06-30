! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Edge value estimation for high-order reconstruction
module regrid_edge_values

use MOM_error_handler, only : MOM_error, FATAL
use regrid_solvers, only : solve_linear_system, linear_solver
use regrid_solvers, only : solve_tridiagonal_system, solve_diag_dominant_tridiag
use polynomial_functions, only : evaluation_polynomial

implicit none ; private

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public bound_edge_values, average_discontinuous_edge_values, check_discontinuous_edge_values
public edge_values_explicit_h2, edge_values_explicit_h4, edge_values_explicit_h4cw
public edge_values_implicit_h4, edge_values_implicit_h6
public edge_slopes_implicit_h3, edge_slopes_implicit_h5

! The following parameter is used to avoid singular matrices for boundary
! extrapolation. It is needed only in the case where thicknesses vanish
! to a small enough values such that the eigenvalues of the matrix can not
! be separated.
real, parameter :: hMinFrac      = 1.e-5  !< A minimum fraction for min(h)/sum(h) [nondim]


  interface
module subroutine bound_edge_values( N, h, u, edge_val, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Potentially modified edge values [A]; the
                                           !! second index is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

end subroutine bound_edge_values
module subroutine average_discontinuous_edge_values( N, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values that may be modified [A]; the
                                           !! second index is for the two edges of each cell.
  ! Local variables

  ! Loop on interior edges
end subroutine average_discontinuous_edge_values
module subroutine check_discontinuous_edge_values( N, u, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: u !< cell averages in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Cell edge values [A]; the
                                           !! second index is for the two edges of each cell.
  ! Local variables

end subroutine check_discontinuous_edge_values
module subroutine edge_values_explicit_h2( N, h, u, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the
                                           !! second index is for the two edges of each cell.

  ! Local variables

  ! Boundary edge values are simply equal to the boundary cell averages
end subroutine edge_values_explicit_h2
module subroutine edge_values_explicit_h4( N, h, u, edge_val, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
                                 ! in units that vary with the second (j) index as [H^j]
                                 ! with the index (j) as [A H^(j-1)]

end subroutine edge_values_explicit_h4
module subroutine edge_values_explicit_h4cw( N, h, u, edge_val, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val  !< Returned edge values [A]; the second index
                                                   !! is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]

  ! Local variables

  ! Set the thicknesses for very thin layers to some minimum value.
end subroutine edge_values_explicit_h4cw
module subroutine edge_values_implicit_h4( N, h, u, edge_val, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
                                        ! in units that vary with the second (j) index as [H^j]
                                        ! with the index (j) as [A H^(j-1)]

end subroutine edge_values_implicit_h4
module subroutine end_value_h4(dz, u, Csys)
  real, dimension(4), intent(in)  :: dz    !< The thicknesses of 4 layers, starting at the edge [H].
                                           !! The values of dz must be positive.
  real, dimension(4), intent(in)  :: u     !< The average properties of 4 layers, starting at the edge [A]
  real, dimension(4), intent(out) :: Csys  !< The four coefficients of a 4th order polynomial fit
                                           !! of u as a function of z [A H-(n-1)]

  ! Local variables
                          ! The units of Wt vary with the second index as [H-(n-1)].
  ! real :: I_h1          ! The inverse of the a thickness [H-1]

  ! These are only used for code verification
  ! real, dimension(4) :: Atest  ! The  coefficients of an expression that is being tested.
  ! real :: zavg, u_mag, c_mag
  ! character(len=128) :: mesg
  ! real, parameter :: C1_12 = 1.0 / 12.0

 ! if ((dz(1) == dz(2)) .and. (dz(1) == dz(3)) .and. (dz(1) == dz(4))) then
 !   ! There are simple closed-form expressions in this case
 !   I_h1 = 0.0 ; if (dz(1) > 0.0) I_h1 = 1.0 / dz(1)
 !   Csys(1) = u(1) + (-13.0 * (u(2)-u(1)) + 10.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25*C1_3)
 !   Csys(2) = (35.0 * (u(2)-u(1)) - 34.0 * (u(3)-u(2)) + 11.0 * (u(4)-u(3))) * (0.25*C1_3 * I_h1)
 !   Csys(3) = (-5.0 * (u(2)-u(1)) + 8.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25 * I_h1**2)
 !   Csys(4) = ((u(2)-u(1)) - 2.0 * (u(3)-u(2)) + (u(4)-u(3))) * (0.5*C1_3)
 ! else

  ! Express the coefficients as sums of the differences between properties of successive layers.

end subroutine end_value_h4
module subroutine edge_slopes_implicit_h3( N, h, u, edge_slopes, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_slopes !< Returned edge slopes [A H-1]; the
                                           !! second index is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
                                        ! in units that vary with the second (j) index as [H^j]
                                        ! index (j) as [A H^(j-1)]
                                        ! in units that vary with the index (j) as [A H^(j-2)]

end subroutine edge_slopes_implicit_h3
module subroutine edge_slopes_implicit_h5( N, h, u, edge_slopes, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_slopes !< Returned edge slopes [A H-1]; the
                                           !! second index is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

! -----------------------------------------------------------------------------
! Fifth-order implicit estimates of edge slopes are based on a four-cell,
! three-edge stencil. A tridiagonal system is set up and is based on
! expressing the edge slopes in terms of neighboring cell averages.
!
! The generic relationship is
!
! \alpha u'_{i-1/2} + u'_{i+1/2} + \beta u'_{i+3/2} =
! a \bar{u}_{i-1} + b \bar{u}_i + c \bar{u}_{i+1} + d \bar{u}_{i+2}
!
! and the stencil looks like this
!
!         i-1     i     i+1    i+2
!   ..--o------o------o------o------o--..
!            i-1/2  i+1/2  i+3/2
!
! In this routine, the coefficients \alpha, \beta, a, b, c and d are
! computed, the tridiagonal system is built, boundary conditions are
! prescribed and the system is solved to yield edge-value estimates.
!
! Note that the centered stencil only applies to edges 3 to N-1 (edges are
! numbered 1 to n+1), which yields N-3 equations for N+1 unknowns. Two other
! equations are written by using a right-biased stencil for edge 2 and a
! left-biased stencil for edge N. The prescription of boundary conditions
! (using sixth-order polynomials) closes the system.
!
! CAUTION: For each edge, in order to determine the coefficients of the
!          implicit expression, a 6x6 linear system is solved. This may
!          become computationally expensive if regridding is carried out
!          often. Figuring out closed-form expressions for these coefficients
!          on nonuniform meshes turned out to be intractable.
! -----------------------------------------------------------------------------

  ! Local variables
                                        ! in units that might vary with the second (j) index as [H^j]
                                        ! units that sometimes vary with the intex (j) as [H^(j-1)] or [H^j]
                                        ! or might be [A]

  ! Loop on cells (except the first and last ones)
end subroutine edge_slopes_implicit_h5
module subroutine edge_values_implicit_h6( N, h, u, edge_val, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val  !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
                                        ! in units that might vary with the second (j) index as [H^j]
                                        ! units that sometimes vary with the intex (j) as [H^(j-1)] or [H^j]
                                        ! or might be [A]
                                        ! coefficients of a fit polynomial in units that vary with the
                                        ! index (j) as [A H^(j-1)]

  ! Loop on interior cells
end subroutine edge_values_implicit_h6
module subroutine test_line(msg, N, A, C, R, mag, tol)
  character(len=*),   intent(in) :: msg  !< An identifying message for this test
  integer,            intent(in) :: N    !< The number of points in the system
  real, dimension(4), intent(in) :: A    !< One of the two vectors being multiplied in arbitrary units [A]
  real, dimension(4), intent(in) :: C    !< One of the two vectors being multiplied in arbitrary units [B]
  real,               intent(in) :: R    !< The expected solution of the equation [A B]
  real,               intent(in) :: mag  !< The magnitude of leading order terms in this line [A B]
  real, optional,     intent(in) :: tol  !< The fractional tolerance for the sums [nondim]


end subroutine test_line
  end interface

end module regrid_edge_values
