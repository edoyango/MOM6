! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Vertical interpolation for regridding
module regrid_interp

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase

use regrid_edge_values, only : edge_values_explicit_h2, edge_values_explicit_h4
use regrid_edge_values, only : edge_values_explicit_h4cw
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_values, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5

use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PPM_functions, only : PPM_monotonicity
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

use P1M_functions, only : P1M_interpolation, P1M_boundary_extrapolation
use P3M_functions, only : P3M_interpolation, P3M_boundary_extrapolation

implicit none ; private

!> Control structure for regrid_interp module
type, public :: interp_CS_type ; private

  !> The following parameter is only relevant when used with the target
  !! interface densities regridding scheme. It indicates which interpolation
  !! to use to determine the grid.
  integer :: interpolation_scheme = -1

  !> Indicate whether high-order boundary extrapolation should be used within
  !! boundary cells
  logical :: boundary_extrapolation

  !> The vintage of the expressions to use for regridding
  integer :: answer_date = 99991231
end type interp_CS_type

public regridding_set_ppolys, build_and_interpolate_grid
public set_interp_scheme, set_interp_extrap, set_interp_answer_date

! List of interpolation schemes
integer, parameter :: INTERPOLATION_P1M_H2     = 0 !< O(h^2)
integer, parameter :: INTERPOLATION_P1M_H4     = 1 !< O(h^2)
integer, parameter :: INTERPOLATION_P1M_IH4    = 2 !< O(h^2)
integer, parameter :: INTERPOLATION_PLM        = 3 !< O(h^2)
integer, parameter :: INTERPOLATION_PPM_CW     =10 !< O(h^3)
integer, parameter :: INTERPOLATION_PPM_H4     = 4 !< O(h^3)
integer, parameter :: INTERPOLATION_PPM_IH4    = 5 !< O(h^3)
integer, parameter :: INTERPOLATION_P3M_IH4IH3 = 6 !< O(h^4)
integer, parameter :: INTERPOLATION_P3M_IH6IH5 = 7 !< O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH4IH3 = 8 !< O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH6IH5 = 9 !< O(h^5)

!>@{ Interpolant degrees
integer, parameter :: DEGREE_1 = 1, DEGREE_2 = 2, DEGREE_3 = 3, DEGREE_4 = 4
integer, public, parameter :: DEGREE_MAX = 5
!>@}

!> When the N-R algorithm produces an estimate that lies outside [0,1], the
!! estimate is set to be equal to the boundary location, 0 or 1, plus or minus
!! an offset, respectively, when the derivative is zero at the boundary [nondim].
real, public, parameter    :: NR_OFFSET = 1e-6
!> Maximum number of Newton-Raphson iterations. Newton-Raphson iterations are
!! used to build the new grid by finding the coordinates associated with
!! target densities and interpolations of degree larger than 1.
integer, public, parameter :: NR_ITERATIONS = 8
!> Tolerance for Newton-Raphson iterations (stop when increment falls below this) [nondim]
real, public, parameter    :: NR_TOLERANCE = 1e-12


  interface
module subroutine regridding_set_ppolys(CS, densities, n0, h0, ppoly0_E, ppoly0_S, &
               ppoly0_coefs, degree, h_neglect, h_neglect_edge)
  type(interp_CS_type),  intent(in)    :: CS !< Interpolation control structure
  integer,               intent(in)    :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)    :: densities !< Actual cell densities [A]
  real, dimension(n0),   intent(in)    :: h0 !< cell widths on source grid [H]
  real, dimension(n0,2), intent(inout) :: ppoly0_E  !< Edge value of polynomial [A]
  real, dimension(n0,2), intent(inout) :: ppoly0_S  !< Edge slope of polynomial [A H-1]
  real, dimension(n0,DEGREE_MAX+1), intent(inout) :: ppoly0_coefs !< Coefficients of polynomial [A]
  integer,               intent(inout) :: degree    !< The degree of the polynomials
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions [H]
                                             !! in the same units as h0.
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width
                                             !! for the purpose of edge value calculations [H]
                                             !! in the same units as h0.
  ! Local variables
                      ! calculations in the same units as h0 [H]

end subroutine regridding_set_ppolys
module subroutine interpolate_grid( n0, h0, x0, ppoly0_E, ppoly0_coefs, &
                             target_values, degree, n1, h1, x1, answer_date )
  integer,               intent(in)     :: n0            !< Number of points on source grid
  integer,               intent(in)     :: n1            !< Number of points on target grid
  real, dimension(n0),   intent(in)     :: h0            !< Thicknesses of source grid cells [H]
  real, dimension(n0+1), intent(in)     :: x0            !< Source interface positions [H]
  real, dimension(n0,2), intent(in)     :: ppoly0_E      !< Edge values of interpolating polynomials [A]
  real, dimension(n0,DEGREE_MAX+1), &
                          intent(in)    :: ppoly0_coefs  !< Coefficients of interpolating polynomials [A]
  real, dimension(n1+1),  intent(in)    :: target_values !< Target values of interfaces [A]
  integer,                intent(in)    :: degree        !< Degree of interpolating polynomials
  real, dimension(n1),    intent(inout) :: h1            !< Thicknesses of target grid cells [H]
  real, dimension(n1+1),  intent(inout) :: x1            !< Target interface positions [H]
  integer,      optional, intent(in)    :: answer_date   !< The vintage of the expressions to use

  ! Local variables

  ! Make sure boundary coordinates of new grid coincide with boundary
  ! coordinates of previous grid
end subroutine interpolate_grid
module subroutine build_and_interpolate_grid(CS, densities, n0, h0, x0, target_values, &
                                      n1, h1, x1, h_neglect, h_neglect_edge)
  type(interp_CS_type),  intent(in)    :: CS  !< A control structure for regrid_interp
  integer,               intent(in)    :: n0  !< The number of points on the input grid
  integer,               intent(in)    :: n1  !< The number of points on the output grid
  real, dimension(n0),   intent(in)    :: densities !< Input cell densities [R ~> kg m-3]
  real, dimension(n1+1), intent(in)    :: target_values !< Target values of interfaces [R ~> kg m-3]
  real, dimension(n0),   intent(in)    :: h0  !< Initial cell widths usually in [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(n0+1), intent(in)    :: x0  !< Source interface positions [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(n1),   intent(inout) :: h1  !< Output cell widths [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(n1+1), intent(inout) :: x1  !< Target interface positions [H ~> m or kg m-2] or [Z ~> m]
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the
                                           !! purpose of cell reconstructions in the same
                                           !! units as h0 [H ~> m or kg m-2] or [Z ~> m].
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the
                                           !! purpose of edge value calculations in the same
                                           !! units as h0 [H ~> m or kg m-2] or [Z ~> m]


end subroutine build_and_interpolate_grid
module function get_polynomial_coordinate( N, h, x_g, edge_values, ppoly_coefs, &
                                    target_value, degree, answer_date ) result ( x_tgt )
  ! Arguments
  integer,              intent(in) :: N            !< Number of grid cells
  real, dimension(N),   intent(in) :: h            !< Grid cell thicknesses [H]
  real, dimension(N+1), intent(in) :: x_g          !< Grid interface locations [H]
  real, dimension(N,2), intent(in) :: edge_values  !< Edge values of interpolating polynomials [A]
  real, dimension(N,DEGREE_MAX+1), intent(in) :: ppoly_coefs  !< Coefficients of interpolating polynomials [A]
  real,                 intent(in) :: target_value !< Target value to find position for [A]
  integer,              intent(in) :: degree       !< Degree of the interpolating polynomials
  integer,              intent(in) :: answer_date  !< The vintage of the expressions to use
  real                             :: x_tgt        !< The position of x_g at which target_value is found [H]

  ! Local variables
!   real                        :: x           ! global target coordinate [nondim]

end function get_polynomial_coordinate
integer module function interpolation_scheme(interp_scheme)
  character(len=*), intent(in) :: interp_scheme !< Name of the interpolation scheme
        !! Valid values include "P1M_H2", "P1M_H4", "P1M_IH2", "PLM", "PPM_CW", "PPM_H4",
        !!   "PPM_IH4", "P3M_IH4IH3", "P3M_IH6IH5", "PQM_IH4IH3", and "PQM_IH6IH5"

end function interpolation_scheme
module subroutine set_interp_scheme(CS, interp_scheme)
  type(interp_CS_type), intent(inout) :: CS  !< A control structure for regrid_interp
  character(len=*),     intent(in) :: interp_scheme !< Name of the interpolation scheme
        !! Valid values include "P1M_H2", "P1M_H4", "P1M_IH2", "PLM", "PPM_CW", "PPM_H4",
        !!   "PPM_IH4", "P3M_IH4IH3", "P3M_IH6IH5", "PQM_IH4IH3", and "PQM_IH6IH5"

end subroutine set_interp_scheme
module subroutine set_interp_extrap(CS, extrap)
  type(interp_CS_type), intent(inout) :: CS  !< A control structure for regrid_interp
  logical,              intent(in)    :: extrap !< Indicate whether high-order boundary
                                             !! extrapolation should be used in boundary cells

end subroutine set_interp_extrap
module subroutine set_interp_answer_date(CS, answer_date)
  type(interp_CS_type), intent(inout) :: CS  !< A control structure for regrid_interp
  integer,              intent(in)    :: answer_date !< An integer encoding the vintage of
                                             !! the expressions to use for regridding

end subroutine set_interp_answer_date
  end interface

end module regrid_interp
