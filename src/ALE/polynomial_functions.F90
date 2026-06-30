! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Polynomial functions
module polynomial_functions

implicit none ; private

public :: evaluation_polynomial, integration_polynomial, first_derivative_polynomial


  interface
real module function evaluation_polynomial( coeff, ncoef, x )
  real, dimension(:), intent(in) :: coeff !< The coefficients of the polynomial, in units that
                                          !! vary with the index k as [A H^(k-1)]
  integer,            intent(in) :: ncoef !< The number of polynomial coefficients
  real,               intent(in) :: x     !< The position at which to evaluate the polynomial
                                          !! in arbitrary thickness units [H]
  ! Local variables

end function evaluation_polynomial
real module function first_derivative_polynomial( coeff, ncoef, x )
  real, dimension(:), intent(in) :: coeff !< The coefficients of the polynomial, in units that
                                          !! vary with the index k as [A H^(k-1)]
  integer,            intent(in) :: ncoef !< The number of polynomial coefficients
  real, intent(in)               :: x     !< The position at which to evaluate the derivative
                                          !! in arbitrary thickness units [H]
  ! Local variables

end function first_derivative_polynomial
real module function integration_polynomial( xi0, xi1, Coeff, npoly )
  real,               intent(in) :: xi0   !< The lower bound of the integral in arbitrary
                                          !! thickness units [H]
  real,               intent(in) :: xi1   !< The upper bound of the integral in arbitrary
                                          !! thickness units [H]
  real, dimension(:), intent(in) :: Coeff !< The coefficients of the polynomial, in units that
                                          !! vary with the index k as [A H^(k-1)]
  integer,            intent(in) :: npoly !< The degree of the polynomial
  ! Local variables

end function integration_polynomial
  end interface

end module polynomial_functions
