! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows Colella and Woodward, 1984 \cite colella1984, with cells
!! resorting to PCM for extrema including the first and last cells in column.
!! The cell-wise reconstructions are limited so that the edge values (which are also the extrema
!! in a cell) are bounded by the neighboring cell means.
!! This does not yield monotonic profiles for the general remapping problem.
module Recon1d_PLM_CW

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PLM_CW, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - x()                       *locally defined
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PLM_CW

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]

contains
  !> Implementation of the PLM_CW initialization
  procedure :: init => init
  !> Implementation of the PLM_CW reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PLM_CW average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PLM_CW reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PLM_CW reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of solver for x: f(x)=t
  procedure :: x => x
  !> Implementation of deallocation for PLM_CW
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PLM_CW reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PLM_CW reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PLM_CW


  interface
module subroutine init(this, n, h_neglect, check)
  class(PLM_CW),     intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(PLM_CW), intent(inout) :: this !< This reconstruction
  real,          intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                                    ! differences across the cell [A]

end subroutine reconstruct
real module function f(this, k, x)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: x    !< Non-dimensional position within element [nondim]

end function f
real module function dfdx(this, k, x)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: x    !< Non-dimensional position within element [nondim]

end function dfdx
real module function x(this, k, t)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: t    !< Value to solve for [A]

end function x
real module function average(this, k, xa, xb)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,          intent(in) :: xb   !< End of averaging interval on element (0 to 1)

  ! This form is not guaranteed to be bounded by {ul,ur}
! u_a = this%ul(k) * ( 1. - xa ) + this%ur(k) * xa
! u_b = this%ul(k) * ( 1. - xb ) + this%ur(k) * xb
! average = 0.5 * ( u_a + u_b )

  ! Mid-point between xa and xb
end function average
module subroutine destroy(this)
  class(PLM_CW), intent(inout) :: this !< This reconstruction

end subroutine destroy
logical module function check_reconstruction(this, h, u)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  real,          intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PLM_CW), intent(inout) :: this    !< This reconstruction
  logical,       intent(in)    :: verbose !< True, if verbose
  integer,       intent(in)    :: stdout  !< I/O channel for stdout
  integer,       intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PLM_CW
