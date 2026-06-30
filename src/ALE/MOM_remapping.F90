! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides column-wise vertical remapping functions
module MOM_remapping

! Original module written by Laurent White, 2008.06.09

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase
use numerical_testing_type, only : testing
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_explicit_h4cw
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_values, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PPM_functions, only : PPM_monotonicity
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1
use MOM_hybgen_remap, only : hybgen_plm_coefs, hybgen_ppm_coefs, hybgen_weno_coefs

use Recon1d_type, only : Recon1d
use Recon1d_PCM, only : PCM
use Recon1d_PLM_CW, only : PLM_CW
use Recon1d_PLM_hybgen, only : PLM_hybgen
use Recon1d_PLM_CWK, only : PLM_CWK
use Recon1d_MPLM_CWK, only : MPLM_CWK
use Recon1d_EMPLM_CWK, only : EMPLM_CWK
use Recon1d_MPLM_WA, only : MPLM_WA
use Recon1d_EMPLM_WA, only : EMPLM_WA
use Recon1d_MPLM_WA_poly, only : MPLM_WA_poly
use Recon1d_EMPLM_WA_poly, only : EMPLM_WA_poly
use Recon1d_PPM_CW, only : PPM_CW
use Recon1d_PPM_hybgen, only : PPM_hybgen
use Recon1d_PPM_CWK, only : PPM_CWK
use Recon1d_EPPM_CWK, only : EPPM_CWK
use Recon1d_PPM_H4_2019, only : PPM_H4_2019
use Recon1d_PPM_H4_2018, only : PPM_H4_2018
use Recon1d_PLM_WLS, only : PLM_WLS

implicit none ; private

!> Container for remapping parameters
type, public :: remapping_CS ; private
  !> Determines which reconstruction to use
  integer :: remapping_scheme = -911
  !> Degree of polynomial reconstruction
  integer :: degree = 0
  !> If true, extrapolate boundaries
  logical :: boundary_extrapolation = .true.
  !> If true, reconstructions are checked for consistency.
  logical :: check_reconstruction = .false.
  !> If true, the result of remapping are checked for conservation and bounds.
  logical :: check_remapping = .false.
  !> If true, the intermediate values used in remapping are forced to be bounded.
  logical :: force_bounds_in_subcell = .false.
  !> If true, impose bounds on the remapping from sub-cells to target grid
  logical :: force_bounds_in_target = .true.
  !> If true, impose bounds on the remapping from non-vanished sub-cells to target grid
  logical :: better_force_bounds_in_target = .false.
  !> If true, calculate and use an offset when summing sub-cells to the target grid
  logical :: offset_tgt_summation = .false.
  !> The vintage of the expressions to use for remapping. Values below 20190101 result
  !! in the use of older, less accurate expressions.
  integer :: answer_date = 99991231
  !> If true, use the OM4 version of the remapping algorithm that makes poor assumptions
  !! about the reconstructions in top and bottom layers of the source grid
  logical :: om4_remap_via_sub_cells = .false.

  !> A negligibly small width for the purpose of cell reconstructions in the same units
  !! as the h0 argument to remapping_core_h [H]
  real :: h_neglect
  !> A negligibly small width for the purpose of edge value calculations in the same units
  !! as the h0 argument to remapping_core_h [H]
  real :: h_neglect_edge

  !> If true, do some debugging as operations proceed
  logical :: debug = .false.

  !> The instance of the actual equation of state
  class(Recon1d), pointer :: reconstruction => Null()
end type

! The following routines are visible to the outside world
public remapping_core_h, remapping_core_w
public initialize_remapping, end_remapping, remapping_set_param, extract_member_remapping_CS
public remapping_unit_tests, build_reconstructions_1d, average_value_ppoly
public interpolate_column, reintegrate_column, dzFromH1H2

! The following are private parameter constants
integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
integer, parameter  :: REMAPPING_PLM        = 2 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PLM_HYBGEN = 3 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PPM_CW     =10 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_H4     = 4 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_IH4    = 5 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_HYBGEN = 6 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_WENO_HYBGEN= 7 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 8 !< O(h^4) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 9 !< O(h^5) remapping scheme
integer, parameter  :: REMAPPING_VIA_CLASS  =99 !< Scheme is controlled by Recon1d class

integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

!> Documentation for external callers
character(len=360), public :: remappingSchemesDoc = &
                 "PCM         (1st-order accurate)\n"//&
                 "PLM         (2nd-order accurate)\n"//&
                 "PLM_HYBGEN  (2nd-order accurate)\n"//&
                 "PPM_H4      (3rd-order accurate)\n"//&
                 "PPM_IH4     (3rd-order accurate)\n"//&
                 "PPM_HYBGEN  (3rd-order accurate)\n"//&
                 "WENO_HYBGEN (3rd-order accurate)\n"//&
                 "PQM_IH4IH3  (4th-order accurate)\n"//&
                 "PQM_IH6IH5  (5th-order accurate)\n"
character(len=3), public :: remappingDefaultScheme = "PLM" !< Default remapping method


  interface
module subroutine remapping_set_param(CS, remapping_scheme, boundary_extrapolation,  &
               check_reconstruction, check_remapping, force_bounds_in_subcell, &
               force_bounds_in_target, better_force_bounds_in_target, offset_tgt_summation, &
               om4_remap_via_sub_cells, answers_2018, answer_date, nk, &
               h_neglect, h_neglect_edge)
  type(remapping_CS),         intent(inout) :: CS !< Remapping control structure
  character(len=*), optional, intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,          intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,          intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,          intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,          intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,          intent(in)    :: force_bounds_in_target !< Force target values to be bounded
  logical, optional,          intent(in)    :: better_force_bounds_in_target !< Force target values to be bounded
  logical, optional,          intent(in)    :: offset_tgt_summation !< Use an offset when summing sub-cells
  logical, optional,          intent(in)    :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm
  logical, optional,          intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,          intent(in)    :: answer_date  !< The vintage of the expressions to use
  real,    optional,          intent(in)    :: h_neglect !< A negligibly small width for the purpose of cell
                                                         !! reconstructions in the same units as the h0 argument
                                                         !! to remapping_core_h [H]
  real,    optional,          intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of edge
                                                         !! value calculations in the same units as as the h0
                                                         !! argument to remapping_core_h [H]
  integer, optional,          intent(in)    :: nk !< Number of levels to initialize reconstruction class with

end subroutine remapping_set_param
module subroutine extract_member_remapping_CS(CS, remapping_scheme, degree, boundary_extrapolation, check_reconstruction, &
                                       check_remapping, force_bounds_in_subcell, force_bounds_in_target, &
                                       better_force_bounds_in_target, offset_tgt_summation)
  type(remapping_CS), intent(in) :: CS !< Control structure for remapping module
  integer, optional, intent(out) :: remapping_scheme        !< Determines which reconstruction scheme to use
  integer, optional, intent(out) :: degree                  !< Degree of polynomial reconstruction
  logical, optional, intent(out) :: boundary_extrapolation  !< If true, extrapolate boundaries
  logical, optional, intent(out) :: check_reconstruction    !< If true, reconstructions are checked for consistency.
  logical, optional, intent(out) :: check_remapping         !< If true, the result of remapping are checked
                                                            !!  for conservation and bounds.
  logical, optional, intent(out) :: force_bounds_in_subcell !< If true, the intermediate values used in
                                                            !! remapping are forced to be bounded.
  logical, optional, intent(out) :: force_bounds_in_target  !< Force target values to be bounded
  logical, optional, intent(out) :: better_force_bounds_in_target  !< Force target values to be bounded
  logical, optional, intent(out) :: offset_tgt_summation    !< Use an offset when summing sub-cells

end subroutine extract_member_remapping_CS
module subroutine remapping_core_h(CS, n0, h0, u0, n1, h1, u1, net_err, PCM_cell)
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1), intent(in)  :: h1 !< Cell widths on target grid [H]
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid [A]
  real, optional,      intent(out) :: net_err !< Error in total column [A H]
  logical, dimension(n0), optional, intent(in) :: PCM_cell !< If present, use PCM remapping for
                                         !! cells in the source grid where this is true.
  ! Local variables
  ! For error checking/debugging

  ! Calculate sub-layer thicknesses and indices connecting sub-layers to source and target grids
  ! Sets: h_sub, h0_eff, isrc_start, isrc_end, isrc_max, isub_src, itgt_start, itgt_end
end subroutine remapping_core_h
module subroutine remapping_core_w( CS, n0, h0, u0, n1, dx, u1)
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,               intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1+1), intent(in)  :: dx !< Cell widths on target grid [H]
  real, dimension(n1),   intent(out) :: u1 !< Cell averages on target grid [A]

  ! Local variables
  ! For error checking/debugging

end subroutine remapping_core_w
module subroutine build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, &
                                     ppoly_r_E, ppoly_r_S, iMethod, h_neglect, &
                                     h_neglect_edge, PCM_cell, debug )
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,CS%degree+1), &
                         intent(out) :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_E !< Edge value of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_S !< Edge slope of polynomial [A H-1]
  integer,               intent(out) :: iMethod !< Integration method
  real,                  intent(in)  :: h_neglect !< A negligibly small width for the
                                         !! purpose of cell reconstructions
                                         !! in the same units as h0 [H]
  real, optional,        intent(in)  :: h_neglect_edge !< A negligibly small width for the purpose
                                         !! of edge value calculations in the same units as h0 [H].
                                         !! The default is h_neglect.
  logical, optional,     intent(in)  :: PCM_cell(n0) !< If present, use PCM remapping for
                                         !! cells from the source grid where this is true.
  logical, optional,     intent(in) :: debug !< If true, enable debugging

  ! Local variables
                      ! calculations in the same units as h0 [H]

end subroutine build_reconstructions_1d
module subroutine check_reconstructions_1d(n0, h0, u0, deg, boundary_extrapolation, &
                                    ppoly_r_coefs, ppoly_r_E)
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(in)  :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2),    intent(in)  :: ppoly_r_E !< Edge value of polynomial [A]
  ! Local variables

end subroutine check_reconstructions_1d
module subroutine intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                    isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h1(n1)  !< Target grid widths (size n1) [H]
  real,    intent(out) :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(out) :: h0_eff(n0) !< Effective thickness of source cells [H]
  integer, intent(out) :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(out) :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(out) :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(out) :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(out) :: itgt_end(n1) !< Index of last sub-cell within each target cell
  integer, intent(out) :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  ! Local variables
  ! For error checking/debugging

  ! Initialize algorithm
end subroutine intersect_src_tgt_grids
module subroutine remap_src_to_sub_grid_om4(n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, h_sub, &
                                 h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                                 method, force_bounds_in_subcell, u_sub, uh_sub, u02_err)
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,    intent(in)  :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: h0_eff(n0) !< Effective thickness of source cells [H]
  integer, intent(in)  :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  integer, intent(in)  :: method  !< Remapping scheme to use
  logical, intent(in)  :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,    intent(out) :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables
  ! For error checking/debugging

end subroutine remap_src_to_sub_grid_om4
module subroutine remap_src_to_sub_grid(n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, h_sub, &
                                 isrc_start, isrc_end, isrc_max, isub_src, &
                                 method, force_bounds_in_subcell, u_sub, uh_sub, u02_err)
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,    intent(in)  :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  integer, intent(in)  :: method  !< Remapping scheme to use
  logical, intent(in)  :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,    intent(out) :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables
  ! For error checking/debugging

end subroutine remap_src_to_sub_grid
module subroutine remap_sub_to_tgt_grid_om4(n0, n1, h1, h_sub, u_sub, uh_sub, &
                                 itgt_start, itgt_end, force_bounds_in_target, u1, uh_err)
  integer, intent(in)  :: n0     !< Number of cells in source grid
  integer, intent(in)  :: n1     !< Number of cells in target grid
  real,    intent(in)  :: h1(n1) !< Target grid widths (size n1) [H]
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(in)  :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  integer, intent(in)  :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(in)  :: itgt_end(n1) !< Index of last sub-cell within each target cell
  logical, intent(in)  :: force_bounds_in_target !< Force sub-cell values to be bounded
  real,    intent(out) :: u1(n1) !< Target cell averages (size n1) [A]
  real,    intent(out) :: uh_err !< Estimate of bound on error in sum of u*h [A H]
  ! Local variables

end subroutine remap_sub_to_tgt_grid_om4
module subroutine remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, &
                                 itgt_start, itgt_end, force_bounds_in_target, &
                                 better_force_bounds_in_target, offset_summation, u1, uh_err)
  integer, intent(in)  :: n0     !< Number of cells in source grid
  integer, intent(in)  :: n1     !< Number of cells in target grid
  real,    intent(in)  :: h1(n1) !< Target grid widths (size n1) [H]
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(in)  :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  integer, intent(in)  :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(in)  :: itgt_end(n1) !< Index of last sub-cell within each target cell
  logical, intent(in)  :: force_bounds_in_target !< Force sub-cell values to be bounded
  logical, intent(in)  :: better_force_bounds_in_target !< Force sub-cell values to be bounded
  logical, intent(in)  :: offset_summation !< Offset values in summation for accuracy
  real,    intent(out) :: u1(n1) !< Target cell averages (size n1) [A]
  real,    intent(out) :: uh_err !< Estimate of bound on error in sum of u*h [A H]
  ! Local variables

end subroutine remap_sub_to_tgt_grid
module subroutine interpolate_column(nsrc, h_src, u_src, ndest, h_dest, u_dest, mask_edges)
  integer,                  intent(in)    :: nsrc   !< Number of source cells
  real, dimension(nsrc),    intent(in)    :: h_src  !< Thickness of source cells [H]
  real, dimension(nsrc+1),  intent(in)    :: u_src  !< Values at source cell interfaces [A]
  integer,                  intent(in)    :: ndest  !< Number of destination cells
  real, dimension(ndest),   intent(in)    :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest+1), intent(inout) :: u_dest !< Interpolated value at destination cell interfaces [A]
  logical,                  intent(in)    :: mask_edges !< If true, mask the values outside of massless
                                                    !! layers at the top and bottom of the column.

  ! Local variables
                            ! within the source layer [nondim], 0 <= frac_pos <= 1.

  ! The following forces the "do while" loop to do one cycle that will set u1, u2, dh.
end subroutine interpolate_column
module subroutine reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, uh_dest)
  integer,                intent(in)    :: nsrc    !< Number of source cells
  real, dimension(nsrc),  intent(in)    :: h_src   !< Thickness of source cells [H]
  real, dimension(nsrc),  intent(in)    :: uh_src  !< Values at source cell interfaces [A H]
  integer,                intent(in)    :: ndest   !< Number of destination cells
  real, dimension(ndest), intent(in)    :: h_dest  !< Thickness of destination cells [H]
  real, dimension(ndest), intent(inout) :: uh_dest !< Interpolated value at destination cell interfaces [A H]

  ! Local variables

end subroutine reintegrate_column
real module function average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefs, method, i0, xa, xb)
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: u0(n0) !< Cell means [A]
  real,          intent(in)    :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,          intent(in)    :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer,       intent(in)    :: method !< Remapping scheme to use
  integer,       intent(in)    :: i0     !< Source cell index
  real,          intent(in)    :: xa     !< Non-dimensional start position within source cell [nondim]
  real,          intent(in)    :: xb     !< Non-dimensional end position within source cell [nondim]
  ! Local variables

end function average_value_ppoly
module subroutine check_remapped_values(n0, h0, u0, ppoly_r_E, deg, ppoly_r_coefs, &
                                 n1, h1, u1, iMethod, uh_err, caller)
  integer,               intent(in) :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in) :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in) :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,2), intent(in) :: ppoly_r_E  !< Edge values of polynomial fits [A]
  integer,               intent(in) :: deg !< Degree of the piecewise polynomial reconstrution
  real, dimension(n0,deg+1), intent(in) :: ppoly_r_coefs !< Coefficients of the piecewise
                                          !! polynomial reconstructions [A]
  integer,               intent(in) :: n1 !< Number of cells on target grid
  real, dimension(n1),   intent(in) :: h1 !< Cell widths on target grid [H]
  real, dimension(n1),   intent(in) :: u1 !< Cell averages on target grid [A]
  integer,               intent(in) :: iMethod !< An integer indicating the integration method used
  real,                  intent(in) :: uh_err  !< A bound on the error in the sum of u*h as
                                               !! estimated by the remapping code [H A]
  character(len=*),      intent(in) :: caller  !< The name of the calling routine.

  ! Local variables

  ! Check errors and bounds
end subroutine check_remapped_values
module subroutine measure_input_bounds( n0, h0, u0, edge_values, h0tot, h0err, u0tot, u0err, u0min, u0max )
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,2), intent(in)  :: edge_values !< Cell edge values on source grid [A]
  real,                  intent(out) :: h0tot !< Sum of cell widths [H]
  real,                  intent(out) :: h0err !< Magnitude of round-off error in h0tot [H]
  real,                  intent(out) :: u0tot !< Sum of cell widths times values [H A]
  real,                  intent(out) :: u0err !< Magnitude of round-off error in u0tot [H A]
  real,                  intent(out) :: u0min !< Minimum value in reconstructions of u0 [A]
  real,                  intent(out) :: u0max !< Maximum value in reconstructions of u0 [A]
  ! Local variables

end subroutine measure_input_bounds
module subroutine measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
  integer,               intent(in)  :: n1 !< Number of cells on destination grid
  real, dimension(n1),   intent(in)  :: h1 !< Cell widths on destination grid [H]
  real, dimension(n1),   intent(in)  :: u1 !< Cell averages on destination grid [A]
  real,                  intent(out) :: h1tot !< Sum of cell widths [H]
  real,                  intent(out) :: h1err !< Magnitude of round-off error in h1tot [H]
  real,                  intent(out) :: u1tot !< Sum of cell widths times values [H A]
  real,                  intent(out) :: u1err !< Magnitude of round-off error in u1tot [H A]
  real,                  intent(out) :: u1min !< Minimum value in reconstructions of u1 [A]
  real,                  intent(out) :: u1max !< Maximum value in reconstructions of u1 [A]
  ! Local variables

end subroutine measure_output_bounds
module subroutine dzFromH1H2( n1, h1, n2, h2, dx )
  integer,            intent(in)  :: n1 !< Number of cells on source grid
  real, dimension(:), intent(in)  :: h1 !< Cell widths of source grid (size n1) [H]
  integer,            intent(in)  :: n2 !< Number of cells on target grid
  real, dimension(:), intent(in)  :: h2 !< Cell widths of target grid (size n2) [H]
  real, dimension(:), intent(out) :: dx !< Change in interface position (size n2+1) [H]
  ! Local variables

end subroutine dzFromH1H2
module subroutine initialize_remapping( CS, remapping_scheme, boundary_extrapolation, &
                check_reconstruction, check_remapping, force_bounds_in_subcell, &
                force_bounds_in_target, better_force_bounds_in_target, offset_tgt_summation, &
                om4_remap_via_sub_cells, answers_2018, answer_date, nk, &
                h_neglect, h_neglect_edge)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  character(len=*),   intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,  intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,  intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,  intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,  intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,  intent(in)    :: force_bounds_in_target !< Force target values to be bounded
  logical, optional,  intent(in)    :: better_force_bounds_in_target !< Force target values to be bounded
  logical, optional,  intent(in)    :: offset_tgt_summation !< Use an offset when summing sub-cells
  logical, optional,  intent(in)    :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm
  logical, optional,  intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,  intent(in)    :: answer_date  !< The vintage of the expressions to use
  real,    optional,  intent(in)    :: h_neglect !< A negligibly small width for the purpose of cell
                                                 !! reconstructions in the same units as h0 [H]
  real,    optional,  intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of edge
                                                      !! value calculations in the same units as h0 [H].
  integer, optional,  intent(in)    :: nk !< Number of levels to initialize reconstruction class with

  ! Note that remapping_scheme is mandatory for initialize_remapping()
end subroutine initialize_remapping
module subroutine setReconstructionType(string,CS)
  character(len=*),   intent(in)    :: string !< String to parse for method
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  ! Local variables
end subroutine setReconstructionType
module subroutine end_remapping(CS)
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure

end subroutine end_remapping
module subroutine test_interp(test, msg, nsrc, h_src, u_src, ndest, h_dest, u_true)
  type(testing),         intent(inout) :: test   !< Unit testing convenience functions
  character(len=*),         intent(in) :: msg    !< Message to label test
  integer,                  intent(in) :: nsrc   !< Number of source cells
  real, dimension(nsrc),    intent(in) :: h_src  !< Thickness of source cells [H]
  real, dimension(nsrc+1),  intent(in) :: u_src  !< Values at source cell interfaces [A]
  integer,                  intent(in) :: ndest  !< Number of destination cells
  real, dimension(ndest),   intent(in) :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest+1), intent(in) :: u_true !< Correct value at destination cell interfaces [A]
  ! Local variables

  ! Interpolate from src to dest
end subroutine test_interp
module subroutine test_reintegrate(test, msg, nsrc, h_src, uh_src, ndest, h_dest, uh_true)
  type(testing),       intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),       intent(in) :: msg     !< Message to label test
  integer,                intent(in) :: nsrc    !< Number of source cells
  real, dimension(nsrc),  intent(in) :: h_src   !< Thickness of source cells [H]
  real, dimension(nsrc),  intent(in) :: uh_src  !< Values of source cell stuff [A H]
  integer,                intent(in) :: ndest   !< Number of destination cells
  real, dimension(ndest), intent(in) :: h_dest  !< Thickness of destination cells [H]
  real, dimension(ndest), intent(in) :: uh_true !< Correct value of destination cell stuff [A H]
  ! Local variables

  ! Interpolate from src to dest
end subroutine test_reintegrate
module subroutine test_recon_consistency(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local

end subroutine test_recon_consistency
module subroutine test_preserve_uniform(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local

end subroutine test_preserve_uniform
module subroutine test_unchanged_grid(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local

end subroutine test_unchanged_grid
module subroutine compare_two_schemes(test, CS1, CS2, n0, n1, niter, msg)
  type(testing),      intent(inout) :: test  !< Unit testing convenience functions
  type(remapping_CS), intent(inout) :: CS1   !< Remapping control structure configured for
                                             !! original implementation
  type(remapping_CS), intent(inout) :: CS2   !< Remapping control structure configured for
                                             !! class-based implementation
  integer,            intent(in)    :: n0    !< Number of source cells
  integer,            intent(in)    :: n1    !< Number of destination cells
  integer,            intent(in)    :: niter !< Number of randomized columns to try
  character(len=*),   intent(in)    :: msg   !< Message to label test
  ! Local

end subroutine compare_two_schemes
logical module function remapping_unit_tests(verbose, num_comp_samp)
  logical, intent(in) :: verbose !< If true, write results to stdout
  integer, optional, intent(in) :: num_comp_samp !< If present, number of samples to
                                 !! try comparing class-based cade against OM4 code
  ! Local variables

end function remapping_unit_tests
  end interface

end module MOM_remapping
