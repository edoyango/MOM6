! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction following Colella and Woodward, 1984
!!
!! This implementation of PPM follows Colella and Woodward, 1984 \cite colella1984, with
!! cells resorting to PCM for extrema including first and last cells in column. The algorithm was
!! first ported from Hycom as hybgen_ppm_coefs() in the mom_hybgen_remap module. This module is
!! a refactor to facilitate more complete testing and evaluation.
!!
!! The mom_hybgen_remap.hybgen_ppm_coefs() function (reached with "PPM_HYGEN"),
!! regrid_edge_values.edge_values_explicit_h4cw() function followed by ppm_functions.ppm_reconstruction()
!! (reached with "PPM_CW"), are equivalent. Similarly recon1d_ppm_hybgen (this implementation) is equivalent also.
module Recon1d_PPM_hybgen

use Recon1d_type, only : testing
use Recon1d_PPM_CW, only : PPM_CW

implicit none ; private

public PPM_hybgen, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! Implemented by extending recon1d_ppm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_cw.init()
!! - reconstruct()             *locally defined
!! - average()                 *locally defined but calls recon1d_ppm_cw.average()
!! - f()                    -> recon1d_ppm_cw.f()
!! - dfdx()                 -> recon1d_ppm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()           -> recon1d_ppm_cw.unit_tests()
!! - destroy()              -> recon1d_ppm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PPM_CW) :: PPM_hybgen

contains
  !> Implementation of the PPM_hybgen reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PPM_hybgen average over an interval [A]
  procedure :: average => average
  !> Implementation of check reconstruction for the PPM_hybgen reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_hybgen reconstruction
  procedure :: unit_tests => unit_tests

end type PPM_hybgen


  interface
module subroutine reconstruct(this, h, u)
  class(PPM_hybgen), intent(inout) :: this !< This reconstruction
  real,              intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,              intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
real module function average(this, k, xa, xb)
  class(PPM_hybgen), intent(in) :: this !< This reconstruction
  integer,           intent(in) :: k    !< Cell number
  real,              intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,              intent(in) :: xb   !< End of averaging interval on element (0 to 1)

end function average
logical module function check_reconstruction(this, h, u)
  class(PPM_hybgen), intent(in) :: this !< This reconstruction
  real,              intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,              intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(PPM_hybgen), intent(inout) :: this    !< This reconstruction
  logical,           intent(in)    :: verbose !< True, if verbose
  integer,           intent(in)    :: stdout  !< I/O channel for stdout
  integer,           intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
module subroutine hybgen_ppm_coefs(s, h_src, edges, nk, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  real,    intent(in)  :: s(nk)     !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2) !< The PPM interpolation edge values [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping


end subroutine hybgen_ppm_coefs
module subroutine bound_edge_values(N, h, u, edge_val, h_neglect, answer_date)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< Cell widths [H]
  real, dimension(N),   intent(in)    :: u !< Cell averages [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date !< The vintage of the expressions to use


end subroutine bound_edge_values
module subroutine check_discontinuous_edge_values(N, u, edge_val)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: u !< Cell averages [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values [A]


end subroutine check_discontinuous_edge_values
  end interface

end module Recon1d_PPM_hybgen
