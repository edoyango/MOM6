! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A module with intrinsic functions that are used by MOM but are not supported
!! by some compilers.
module MOM_intrinsic_functions

use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
use iso_fortran_env, only : int64, real64

implicit none ; private

public :: invcosh, cuberoot
public :: intrinsic_functions_unit_tests

! Floating point model, if bit layout from high to low is (sign, exp, frac)

integer, parameter :: bias = maxexponent(1.) - 1
  !< The double precision exponent offset
integer, parameter :: signbit = storage_size(1.) - 1
  !< Position of sign bit
integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
  !< Bit size of exponent
integer, parameter :: expbit = signbit - explen
  !< Position of lowest exponent bit
integer, parameter :: fraclen = expbit
  !< Length of fractional part


  interface
module function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh [nondim].  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh  ! The inverse of cosh of x [nondim]

end function invcosh
elemental module function cuberoot(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
              ! in arbitrary units that can grow or shrink with each iteration [B C]
              ! in arbitrary units that can grow or shrink with each iteration [C]
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]


end function cuberoot
pure module subroutine rescale_cbrt(a, x, e_r, s_a)
  real, intent(in) :: a
    !< The number to be rescaled for cube-root computation [A3]
  real, intent(out) :: x
    !< The rescaled value of `a` in the range [0.125, 1) [B3]
  integer(kind=int64), intent(out) :: e_r
    !< The integral component of the cube-root exponent of `a`.
  integer(kind=int64), intent(out) :: s_a
    !< Sign bit of `a`.  A nonzero value indicates negative sign.

    ! Floating point integer representation of `a`
    ! Exponent of `a`
    ! Exponent of `x`

  ! Pack bits of a into xb and extract its exponent and sign.
end subroutine rescale_cbrt
pure module function descale(x, e_a, s_a) result(a)
  real, intent(in) :: x
    !< The rescaled value which is to be restored in ambiguous units [B]
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the unscaled value
  real :: a
    !< Restored value with the corrected exponent and sign in arbitrary units [A]

    ! Bit-packed real number into integer form
    ! Biased exponent of x

  ! Apply the corrected exponent and sign to x.
end function descale
module function intrinsic_functions_unit_tests(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  ! Local variables

end function intrinsic_functions_unit_tests
logical module function Test_cuberoot(verbose, val)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
  ! Local variables

end function Test_cuberoot
  end interface

end module MOM_intrinsic_functions
