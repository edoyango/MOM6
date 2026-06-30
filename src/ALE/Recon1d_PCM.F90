! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> 1D reconstructions using the Piecewise Constant Method (PCM)
module Recon1d_PCM

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PCM

!> PCM (piecewise constant) reconstruction
!!
!! The source for the methods ultimately used by this class are:
!!   init()                    *locally defined
!!   reconstruct()             *locally defined
!!   average()                 *locally defined
!!   f()                       *locally defined
!!   dfdx()                    *locally defined
!! - x()                       *locally defined
!!   check_reconstruction()    *locally defined
!!   unit_tests()              *locally defined
!!   destroy()                 *locally defined
!!   remap_to_sub_grid()    -> Recon1d%remap_to_sub_grid()
!!   init_parent()          -> init()
!!   reconstruct_parent()   -> parent()
type, extends (Recon1d) :: PCM

contains
  !> Implementation of the PCM initialization
  procedure :: init => init
  !> Implementation of the PCM reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PCM average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PCM reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PCM reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of solver for x: f(x)=t
  procedure :: x => x
  !> Implementation of deallocation for PCM
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PCM reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PCM reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PCM


  interface
module subroutine init(this, n, h_neglect, check)
  class(PCM),        intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H].
                                              !! Not used by PCM.
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(PCM), intent(inout) :: this !< This reconstruction
  real,       intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,       intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
real module function f(this, k, x)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: x    !< Non-dimensional position within element [nondim]

end function f
real module function dfdx(this, k, x)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: x    !< Non-dimensional position within element [nondim]

end function dfdx
real module function x(this, k, t)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: t    !< Value to solve for [A]

end function x
real module function average(this, k, xa, xb)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,       intent(in) :: xb   !< End of averaging interval on element (0 to 1)

end function average
module subroutine destroy(this)
  class(PCM), intent(inout) :: this !< This reconstruction

end subroutine destroy
logical module function check_reconstruction(this, h, u)
  class(PCM), intent(in) :: this !< This reconstruction
  real,       intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,       intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PCM), intent(inout) :: this    !< This reconstruction
  logical,    intent(in)    :: verbose !< True, if verbose
  integer,    intent(in)    :: stdout  !< I/O channel for stdout
  integer,    intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PCM
