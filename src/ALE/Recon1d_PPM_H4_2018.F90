! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction with h4 interpolation for edges (2018 version)
!!
!! This implementation of PPM follows White and Adcroft 2008 \cite white2008, with cells
!! resorting to PCM for extrema including first and last cells in column.
!! This scheme differs from Colella and Woodward, 1984 \cite colella1984, in the method
!! of first estimating the fourth-order accurate edge values.
!! This uses numerical expressions that predate a 2019 refactoring.
!! The first and last cells are always limited to PCM.
module Recon1d_PPM_H4_2018

use Recon1d_PPM_H4_2019, only : PPM_H4_2019, testing
use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values
use regrid_solvers, only :  solve_linear_system

implicit none ; private

public PPM_H4_2018, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! Implemented by extending recon1d_ppm_h4_2019.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_h4_2019.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_ppm_h4_2019.average()
!! - f()                    -> recon1d_ppm_h4_2019.f()
!! - dfdx()                 -> recon1d_ppm_h4_2019.dfdx()
!! - check_reconstruction() -> recon1d_ppm_h4_2019.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_ppm_h4_2019.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_ppm_h4_2019.init()
!! - reconstruct_parent()   -> recon1d_ppm_h4_2019.reconstruct()
type, extends (PPM_H4_2019) :: PPM_H4_2018

contains
  !> Implementation of the PPM_H4_2018 reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the PPM_H4_2018 reconstruction
  procedure :: unit_tests => unit_tests

end type PPM_H4_2018


  interface
module subroutine reconstruct(this, h, u)
  class(PPM_H4_2018), intent(inout) :: this !< This reconstruction
  real,               intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,               intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                            ! in units that vary with the second (j) index as [H^j]
                            ! with the index (j) as [A H^(j-1)]

end subroutine reconstruct
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PPM_H4_2018), intent(inout) :: this    !< This reconstruction
  logical,            intent(in)    :: verbose !< True, if verbose
  integer,            intent(in)    :: stdout  !< I/O channel for stdout
  integer,            intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_PPM_H4_2018
