! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction with h4 interpolation for edges
!!
!! This implementation of PPM follows White and Adcroft 2008 \cite white2008, with cells
!! resorting to PCM for extrema including first and last cells in column.
!! This scheme differs from Colella and Woodward, 1984 \cite colella1984, in the method
!! of first estimating the fourth-order accurate edge values.
!! This uses numerical expressions refactored at the beginning of 2019.
!! The first and last cells are always limited to PCM.
module Recon1d_PPM_H4_2019

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PPM_H4_2019, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PPM_H4_2019

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]

contains
  !> Implementation of the PPM_H4_2019 initialization
  procedure :: init => init
  !> Implementation of the PPM_H4_2019 reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PPM_H4_2019 average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PPM_H4_2019 reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PPM_H4_2019 reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PPM_H4_2019
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PPM_H4_2019 reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_H4_2019 reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PPM_H4_2019


  interface
module subroutine init(this, n, h_neglect, check)
  class(PPM_H4_2019),     intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(PPM_H4_2019), intent(inout) :: this !< This reconstruction
  real,          intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                                    ! differences across the cell [A]
                     ! in units that vary with the second (j) index as [H^j]
                     ! with the index (j) as [A H^(j-1)]

end subroutine reconstruct
module subroutine end_value_h4(dz, u, Csys)
  real, intent(in)  :: dz(4)    !< The thicknesses of 4 layers, starting at the edge [H].
                                !! The values of dz must be positive.
  real, intent(in)  :: u(4)     !< The average properties of 4 layers, starting at the edge [A]
  real, intent(out) :: Csys(4)  !< The four coefficients of a 4th order polynomial fit
                                !! of u as a function of z [A H-(n-1)]

  ! Local variables
                          ! The units of Wt vary with the second index as [H-(n-1)].
  ! real :: I_h1          ! The inverse of the a thickness [H-1]

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
real module function f(this, k, x)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,            intent(in) :: k    !< Cell number
  real,               intent(in) :: x    !< Non-dimensional position within element [nondim]

end function f
real module function dfdx(this, k, x)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,            intent(in) :: k    !< Cell number
  real,               intent(in) :: x    !< Non-dimensional position within element [nondim]

end function dfdx
real module function average(this, k, xa, xb)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,          intent(in) :: xb   !< End of averaging interval on element (0 to 1)

end function average
module subroutine destroy(this)
  class(PPM_H4_2019), intent(inout) :: this !< This reconstruction

end subroutine destroy
logical module function check_reconstruction(this, h, u)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  real,          intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PPM_H4_2019), intent(inout) :: this    !< This reconstruction
  logical,       intent(in)    :: verbose !< True, if verbose
  integer,       intent(in)    :: stdout  !< I/O channel for stdout
  integer,       intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PPM_H4_2019
