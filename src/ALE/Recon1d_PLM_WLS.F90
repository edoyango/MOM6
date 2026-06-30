!> Piecewise Linear Method using Weighted Conservative Least Squares 1D reconstruction
module Recon1d_PLM_WLS

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PLM_WLS, testing

!> PLM reconstruction using Weighted Least Squares constrained to conserve for central cell
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - x()                    -> recon1d_type.x()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PLM_WLS

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]
  real, allocatable, private :: slp(:) !< Difference across cell, ur - ul [A].
                              !! This is redundant with ul and ur and not used
                              !! in any evaluations, but is needed for testing.

contains
  !> Implementation of the PLM_WLS initialization
  procedure :: init => init
  !> Implementation of the PLM_WLS reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PLM_WLS average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PLM_WLS reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PLM_WLS reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PLM_WLS
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PLM_WLS reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PLM_WLS reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PLM_WLS


  interface
module subroutine init(this, n, h_neglect, check)
  class(PLM_WLS),    intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(PLM_WLS), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
real module function f(this, k, x)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]

end function f
real module function dfdx(this, k, x)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]

end function dfdx
real module function average(this, k, xa, xb)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)

  ! Mid-point between xa and xb
end function average
module subroutine destroy(this)
  class(PLM_WLS), intent(inout) :: this !< This reconstruction

end subroutine destroy
logical module function check_reconstruction(this, h, u)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
real module function LS_error(this, k, h, u)
  type(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function LS_error
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PLM_WLS), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PLM_WLS
