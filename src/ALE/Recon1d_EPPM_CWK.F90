! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction in model index space with linear
!! extrapolation for first and last cells
!!
!! This implementation of PPM follows Colella and Woodward, 1984, using uniform thickness
!! and with cells resorting to PCM for local extrema. First and last cells use a PLM
!! representation with slope set by matching the edge of the first interior cell.
module Recon1d_EPPM_CWK

use Recon1d_type, only : Recon1d, testing
use Recon1d_PPM_CWK, only : PPM_CWK

implicit none ; private

public EPPM_CWK, testing

!> PPM reconstruction in index space (no grid dependence) with linear extrapolation
!! for first and last cells.
!!
!! Implemented by extending recon1d_ppm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_cwk.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_ppm_cwk.average()
!! - f()                    -> recon1d_ppm_cwk.f()
!! - dfdx()                 -> recon1d_ppm_cwk.dfdx()
!! - check_reconstruction() -> recon1d_ppm_cwk.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_ppm_cwk.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_ppm_cwk.init()
!! - reconstruct_parent()   -> recon1d_ppm_cwk.reconstruct()
type, extends (PPM_CWK) :: EPPM_CWK

contains
  !> Implementation of the EPPM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EPPM_CWK reconstruction
  procedure :: unit_tests => unit_tests

end type EPPM_CWK


  interface
module subroutine reconstruct(this, h, u)
  class(EPPM_CWK), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
logical module function unit_tests(this, verbose, stdout, stderr)
  class(EPPM_CWK), intent(inout) :: this    !< This reconstruction
  logical,         intent(in)    :: verbose !< True, if verbose
  integer,         intent(in)    :: stdout  !< I/O channel for stdout
  integer,         intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_EPPM_CWK
