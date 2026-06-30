! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction in model index space
!!
!! This implementation of PPM follows Colella and Woodward, 1984, using uniform thickness
!! and with cells resorting to PCM for local extrema including the first and last cells.
!!
!! "Fourth order" estimates of edge values use PLM also calculated in index space
!! (i.e. with no grid dependence). First and last PLM slopes are extrapolated.
!! Limiting follows Colella and Woodward thereafter. The high accuracy of this scheme is
!! realized only when the grid-spacing is exactly uniform. This scheme deviates from CW84
!! when the grid spacing is variable.
module Recon1d_PPM_CWK

use Recon1d_type, only : Recon1d, testing
use Recon1d_PLM_CWK, only : PLM_CWK

implicit none ; private

public PPM_CWK, testing

!> PPM reconstruction in index space (no grid dependence).
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
type, extends (Recon1d) :: PPM_CWK

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]
  type(PLM_CWK) :: PLM !< The PLM reconstruction used to estimate edge values

contains
  !> Implementation of the PPM_CWK initialization
  procedure :: init => init
  !> Implementation of the PPM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PPM_CWK average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PPM_CWK reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PPM_CWK reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of solver for x: f(x)=t
  procedure :: x => x
  !> Implementation of deallocation for PPM_CWK
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PPM_CWK reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_CWK reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PPM_CWK


  interface
module subroutine init(this, n, h_neglect, check)
  class(PPM_CWK),     intent(out) :: this      !< This reconstruction
  integer,            intent(in)  :: n         !< Number of cells in this column
  real, optional,     intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional,  intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(PPM_CWK), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
real module function f(this, k, x)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]

end function f
real module function dfdx(this, k, x)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]

end function dfdx
real module function x(this, k, t)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: t    !< Value to solve for [A]

  ! The PPM profile is the quadratic profile: f(x) = ul + (slp+a6)*x - a6*x^2.
  ! Setting f(x)=t gives: -a6*x^2 + (slp+a6)*x + (ul-t) = 0.
  ! In the common parlance of solving a*x^2 + b*x + c = 0, this means
  !  a = -a6;  b = slp+a6; c = ul-t
  ! The quadratic formula x = ( -b +/- sD ) / ( 2a ) with sD = sqrt(b^2-4*a*c)
  ! can suffer from catastrophic cancellation in some scenarios.
  ! A mathematically equivalent form of x = 2c / ( -b -/+ sD ) also can fail.
  ! Usually, to avoid catastrophic cancellation, we use the rule
  !   If b>0 then the two roots are
  !     ra = -(b+sD)/(2a)
  !     rc = -2c/(b+sD)
  !   otherwise if b<0 then the two roots are
  !     ra = (-b+sD)/(2a)
  !     rc = 2c/(-b+sD)
  ! In all expressions, sD and b do not have cancelling contributions due to the signs.
  ! Note that here, if b>0 then c<0, and vice versa, because we are looking
  ! for f(x)=t which shifts "c" by t so that the root we are interested in
  ! falls in the range 0 <= x <= 1 (assuming t falls in ul...ur).
  ! When b>0 and a>0 then -b/(2a)<0 and ra<0<rc, so we need rc
  ! When b>0 and a<0 then -b/(2a)>0 and ra>rc, so we need rc
  ! When b<0 and a>0 then -b/(2a)>0 and ra>rc, so we need rc
  ! When b<0 and a<0 then -b/(2a)<0 and ra<0<rc, so we need rc
  ! So a form that always gives us the root that we want is
  !   x = -2c/(b+sgn(b)*sD)
end function x
real module function average(this, k, xa, xb)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)

end function average
module subroutine destroy(this)
  class(PPM_CWK), intent(inout) :: this !< This reconstruction

end subroutine destroy
logical module function check_reconstruction(this, h, u)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PPM_CWK), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PPM_CWK
