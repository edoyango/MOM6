! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> provides runtime remapping of diagnostics to z star, sigma and
!! rho vertical coordinates.
!!
!! The diag_remap_ctrl type represents a remapping of diagnostics to a particular
!! vertical coordinate. The module is used by the diag mediator module in the
!! following way:
!! 1. diag_remap_init() is called to initialize a diag_remap_ctrl instance.
!! 2. diag_remap_configure_axes() is called to read the configuration file and set up the
!!    vertical coordinate / axes definitions.
!! 3. diag_remap_get_axes_info() returns information needed for the diag mediator to
!!    define new axes for the remapped diagnostics.
!! 4. diag_remap_update() is called periodically (whenever h, T or S change) to either
!!    create or update the target remapping grids.
!! 5. diag_remap_do_remap() is called from within a diag post() to do the remapping before
!!    the diagnostic is written out.


! NOTE: In the following functions, the fields are initially passed using 1-based
! indexing, which are then passed to separate private internal routines that shift
! the indexing to use the same indexing conventions used elsewhere in the MOM6 code.
!
!   * diag_remap_do_remap, which calls do_remap
!   * vertically_reintegrate_diag_field, which calls vertically_reintegrate_field
!   * vertically_interpolate_diag_field, which calls vertically_interpolate_field
!   * horizontally_average_diag_field, which calls horizontally_average_field


module MOM_diag_remap

use MOM_coms,             only : reproducing_sum_EFP, EFP_to_real
use MOM_coms,             only : EFP_type, assignment(=), EFP_sum_across_PEs
use MOM_error_handler,    only : MOM_error, FATAL, assert, WARNING
use MOM_debugging,        only : check_column_integrals
use MOM_diag_manager_infra,only : MOM_diag_axis_init
use MOM_file_parser,      only : get_param, log_param, param_file_type
use MOM_string_functions, only : lowercase, extractWord
use MOM_grid,             only : ocean_grid_type
use MOM_unit_scaling,     only : unit_scale_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : EOS_type
use MOM_remapping,        only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_remapping,        only : interpolate_column, reintegrate_column
use MOM_regridding,       only : regridding_CS, initialize_regridding, end_regridding
use MOM_regridding,       only : set_regrid_params, get_regrid_size
use MOM_regridding,       only : getCoordinateInterfaces, set_h_neglect, set_dz_neglect
use MOM_regridding,       only : get_zlike_CS, get_sigma_CS, get_rho_CS
use regrid_consts,        only : coordinateMode
use coord_zlike,          only : build_zstar_column
use coord_sigma,          only : build_sigma_column
use coord_rho,            only : build_rho_column


implicit none ; private

#include "MOM_memory.h"

public diag_remap_ctrl
public diag_remap_init, diag_remap_end, diag_remap_update, diag_remap_do_remap
public diag_remap_configure_axes, diag_remap_axes_configured
public diag_remap_calc_hmask
public diag_remap_get_axes_info, diag_remap_set_active
public diag_remap_diag_registration_closed
public vertically_reintegrate_diag_field
public vertically_interpolate_diag_field
public horizontally_average_diag_field

!> Represents remapping of diagnostics to a particular vertical coordinate.
!!
!! There is one of these types for each vertical coordinate. The vertical axes
!! of a diagnostic will reference an instance of this type indicating how (or
!! if) the diagnostic should be vertically remapped when being posted.
type :: diag_remap_ctrl
  logical :: configured = .false. !< Whether vertical coordinate has been configured
  logical :: initialized = .false.  !< Whether remappping initialized
  logical :: used = .false.  !< Whether this coordinate actually gets used.
  integer :: vertical_coord = 0 !< The vertical coordinate that we remap to
  character(len=10) :: vertical_coord_name ='' !< The coordinate name as understood by ALE
  logical :: Z_based_coord = .false.  !< If true, this coordinate is based on remapping of
                                      !! geometric distances across layers (in [Z ~> m]) rather
                                      !! than layer thicknesses (in [H ~> m or kg m-2]).  This
                                      !! distinction only matters in non-Boussinesq mode.
  character(len=16) :: diag_coord_name = '' !< A name for the purpose of run-time parameters
  character(len=8) :: diag_module_suffix = '' !< The suffix for the module to appear in diag_table
  type(remapping_CS) :: remap_cs !< Remapping control structure use for this axes
  type(regridding_CS) :: regrid_cs !< Regridding control structure that defines the coordinates for this axes
  integer :: nz = 0 !< Number of vertical levels used for remapping
  real, dimension(:,:,:), allocatable :: h !< Remap grid thicknesses in [H ~> m or kg m-2] or
                                      !! vertical extents in [Z ~> m], depending on the setting of Z_based_coord.
  real, dimension(:,:,:), allocatable :: h_extensive !< Remap grid thicknesses in [H ~> m or kg m-2] or
                                      !! vertical extents in [Z ~> m] for remapping extensive variables
  integer :: interface_axes_id = 0 !< Vertical axes id for remapping at interfaces
  integer :: layer_axes_id = 0 !< Vertical axes id for remapping on layers
  logical :: om4_remap_via_sub_cells !< Use the OM4-era ramap_via_sub_cells
  integer :: answer_date      !< The vintage of the order of arithmetic and expressions
                              !! to use for remapping.  Values below 20190101 recover
                              !! the answers from 2018, while higher values use more
                              !! robust forms of the same remapping expressions.

end type diag_remap_ctrl


  interface
module subroutine diag_remap_init(remap_cs, coord_tuple, om4_remap_via_sub_cells, answer_date, GV)
  type(diag_remap_ctrl), intent(inout) :: remap_cs    !< Diag remapping control structure
  character(len=*),      intent(in)    :: coord_tuple !< A string in form of
                                                      !! MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME
  logical,               intent(in)    :: om4_remap_via_sub_cells !< Use the OM4-era ramap_via_sub_cells
  integer,               intent(in)    :: answer_date !< The vintage of the order of arithmetic and expressions
                                                      !! to use for remapping.  Values below 20190101 recover
                                                      !! the answers from 2018, while higher values use more
                                                      !! robust forms of the same remapping expressions.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean vertical grid structure, used here to evaluate
                                                      !! whether the model is in non-Boussinesq mode.

end subroutine diag_remap_init
module subroutine diag_remap_end(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

end subroutine diag_remap_end
module subroutine diag_remap_diag_registration_closed(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

end subroutine diag_remap_diag_registration_closed
module subroutine diag_remap_set_active(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

end subroutine diag_remap_set_active
module subroutine diag_remap_configure_axes(remap_cs, G, GV, US, param_file)
  type(diag_remap_ctrl),   intent(inout) :: remap_cs !< Diag remap control structure
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid type
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file structure

  ! Local variables
end subroutine diag_remap_configure_axes
module subroutine diag_remap_get_axes_info(remap_cs, nz, id_layer, id_interface)
  type(diag_remap_ctrl), intent(in) :: remap_cs !< Diagnostic coordinate control structure
  integer, intent(out) :: nz !< Number of vertical levels for the coordinate
  integer, intent(out) :: id_layer !< 1D-axes id for layer points
  integer, intent(out) :: id_interface !< 1D-axes id for interface points

end subroutine diag_remap_get_axes_info
module function diag_remap_axes_configured(remap_cs)
  type(diag_remap_ctrl), intent(in) :: remap_cs !< Diagnostic coordinate control structure
  logical :: diag_remap_axes_configured

end function diag_remap_axes_configured
module subroutine diag_remap_update(remap_cs, G, GV, US, h, T, S, eqn_of_state, h_target)
  type(diag_remap_ctrl),   intent(inout) :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),   pointer    :: G  !< The ocean's grid type
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h  !< New thickness in [H ~> m or kg m-2] or [Z ~> m], depending
                                            !! on the value of remap_cs%Z_based_coord
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: T  !< New temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: S  !< New salinities [S ~> ppt]
  type(EOS_type),          intent(in) :: eqn_of_state !< A pointer to the equation of state
  real, dimension(SZI_(G),SZJ_(G),remap_cs%nz), &
                        intent(inout) :: h_target  !< The new diagnostic thicknesses in [H ~> m or kg m-2]
                                            !! or [Z ~> m], depending on the value of remap_cs%Z_based_coord

  ! Local variables
                         ! in units of [H Z-1 ~> 1 or kg m-3] or [nondim], depending on remap_cs%Z_based_coord.

  ! Note that coordinateMode('LAYER') is never 'configured' so will always return here.
end subroutine diag_remap_update
module subroutine diag_remap_do_remap(remap_cs, G, GV, US, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                               mask, field, remapped_field)
  type(diag_remap_ctrl),   intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),   intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(:,:,:),  intent(in)  :: h  !< The current thicknesses [H ~> m or kg m-2] or [Z ~> m],
                                             !! depending on the value of remap_CS%Z_based_coord
  integer, dimension(:,:), intent(in)  :: OBC_u !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at u-points,
                                             !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                             !! or -1 for a Western OBC
  integer, dimension(:,:), intent(in)  :: OBC_v !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at v-points,
                                             !! with a value of 0 for no OBC, 1 for a Northern OBC
                                             !! or -1 for a Southern OBC
  logical,                 intent(in)  :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                 intent(in)  :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:),  pointer     :: mask !< A mask for the field [nondim].
  real, dimension(:,:,:),  intent(in)  :: field(:,:,:) !< The diagnostic field to be remapped [A]
  real, dimension(:,:,:),  intent(out) :: remapped_field !< Field remapped to new coordinate [A]

  ! Local variables

end subroutine diag_remap_do_remap
module subroutine do_remap(remap_cs, G, GV, US, isdf, jsdf, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                    field, remapped_field, mask)
  type(diag_remap_ctrl),   intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),   intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  integer,                 intent(in)  :: isdf !< The starting i-index in memory for field
  integer,                 intent(in)  :: jsdf !< The starting j-index in memory for field
  real, dimension(G%isd:,G%jsd:,:), &
                           intent(in)  :: h  !< The current thicknesses [H ~> m or kg m-2] or [Z ~> m],
                                             !! depending on the value of remap_CS%Z_based_coord
  integer, dimension(G%IsdB:,G%jsd:), &
                           intent(in) :: OBC_u !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at u-points,
                                             !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                             !! or -1 for a Western OBC
  integer, dimension(G%isd:,G%JsdB:), &
                           intent(in) :: OBC_v !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at v-points,
                                             !! with a value of 0 for no OBC, 1 for a Northern OBC
                                             !! or -1 for a Southern OBC
  logical,                 intent(in)  :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                 intent(in)  :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(isdf:,jsdf:,:), &
                           intent(in)  :: field !< The diagnostic field to be remapped [A]
  real, dimension(isdf:,jsdf:,:), &
                           intent(out) :: remapped_field !< Field remapped to new coordinate [A]
  real, dimension(isdf:,jsdf:), &
                 optional, intent(in)  :: mask !< A mask for the field [nondim]

  ! Local variables

end subroutine do_remap
module subroutine diag_remap_calc_hmask(remap_cs, G, mask)
  type(diag_remap_ctrl),  intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),  intent(in)  :: G    !< Ocean grid structure
  real, dimension(G%isd:,G%jsd:,:), &
                          intent(out) :: mask !< h-point mask for target grid [nondim]

  ! Local variables

end subroutine diag_remap_calc_hmask
module subroutine vertically_reintegrate_diag_field(remap_cs, G, h, h_target, OBC_u, OBC_v, &
                                             staggered_in_x, staggered_in_y, mask, field, reintegrated_field)
  type(diag_remap_ctrl),  intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),  intent(in)  :: G        !< Ocean grid structure
  real, dimension(:,:,:), intent(in)  :: h        !< The thicknesses of the source grid [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(:,:,:), intent(in)  :: h_target !< The thicknesses of the target grid [H ~> m or kg m-2] or [Z ~> m]
  integer, dimension(:,:), intent(in) :: OBC_u    !< An array that indicates the presence and direction
                                                  !! of any open boundary conditions at u-points,
                                                  !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                                  !! or -1 for a Western OBC
  integer, dimension(:,:), intent(in) :: OBC_v    !< An array that indicates the presence and direction
                                                  !! of any open boundary conditions at v-points,
                                                  !! with a value of 0 for no OBC, 1 for a Northern OBC
                                                  !! or -1 for a Southern OBC
  logical,                intent(in)  :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in)  :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer     :: mask     !< A mask for the field [nondim].  Note that because this
                                                  !! is a pointer it retains its declared indexing conventions.
  real, dimension(:,:,:), intent(in)  :: field    !<  The diagnostic field to be remapped [A]
  real, dimension(:,:,:), intent(out) :: reintegrated_field !< Field argument remapped to alternative coordinate [A]

  ! Local variables

end subroutine vertically_reintegrate_diag_field
module subroutine vertically_reintegrate_field(remap_cs, G, isdf, jsdf, h, h_target, OBC_u, OBC_v, &
                                        staggered_in_x, staggered_in_y, field, reintegrated_field, mask)
  type(diag_remap_ctrl),  intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),  intent(in)  :: G        !< Ocean grid structure
  integer,                intent(in)  :: isdf     !< The starting i-index in memory for field
  integer,                intent(in)  :: jsdf     !< The starting j-index in memory for field
  real, dimension(G%isd:,G%jsd:,:), &
                          intent(in)  :: h        !< The thicknesses of the source grid [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(G%isd:,G%jsd:,:), &
                          intent(in)  :: h_target !< The thicknesses of the target grid [H ~> m or kg m-2] or [Z ~> m]
  integer, dimension(G%IsdB:,G%jsd:),  &
                           intent(in) :: OBC_u    !< An array that indicates the presence and direction
                                                  !! of any open boundary conditions at u-points,
                                                  !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                                  !! or -1 for a Western OBC
  integer, dimension(G%isd:,G%JsdB:), &
                           intent(in) :: OBC_v    !< An array that indicates the presence and direction
                                                  !! of any open boundary conditions at v-points,
                                                  !! with a value of 0 for no OBC, 1 for a Northern OBC
                                                  !! or -1 for a Southern OBC
  logical,                intent(in)  :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in)  :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(isdf:,jsdf:,:), &
                          intent(in)  :: field   !< The diagnostic field to be remapped [A]
  real, dimension(isdf:,jsdf:,:), &
                          intent(out) :: reintegrated_field !< Field argument remapped to alternative coordinate [A]
  real, dimension(isdf:,jsdf:), &
                optional, intent(in)  :: mask !< A mask for the field [nondim]

  ! Local variables

end subroutine vertically_reintegrate_field
module subroutine vertically_interpolate_diag_field(remap_cs, G, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                                             mask, field, interpolated_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),  intent(in) :: G   !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !< The current thicknesses [H ~> m or kg m-2] or [Z ~> m],
                                            !! depending on the value of remap_cs%Z_based_coord
  integer, dimension(:,:), intent(in) :: OBC_u !< An array that indicates the presence and direction
                                            !! of any open boundary conditions at u-points,
                                            !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                            !! or -1 for a Western OBC
  integer, dimension(:,:), intent(in) :: OBC_v !< An array that indicates the presence and direction
                                            !! of any open boundary conditions at v-points,
                                            !! with a value of 0 for no OBC, 1 for a Northern OBC
                                            !! or -1 for a Southern OBC
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask !< A mask for the field [nondim].  Note that because this
                                             !! is a pointer it retains its declared indexing conventions.
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped [A]
  real, dimension(:,:,:), intent(inout) :: interpolated_field !< Field argument remapped to alternative coordinate [A]

  ! Local variables

end subroutine vertically_interpolate_diag_field
module subroutine vertically_interpolate_field(remap_cs, G, isdf, jsdf, h, OBC_u, OBC_v, &
                                        staggered_in_x, staggered_in_y, field, interpolated_field, mask)
  type(diag_remap_ctrl),  intent(in)  :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),  intent(in)  :: G    !< Ocean grid structure
  integer,                intent(in)  :: isdf !< The starting i-index in memory for field
  integer,                intent(in)  :: jsdf !< The starting j-index in memory for field
  real, dimension(G%isd:,G%jsd:,:), &
                          intent(in)  :: h    !< The current thicknesses [H ~> m or kg m-2] or [Z ~> m],
                                              !! depending on the value of remap_cs%Z_based_coord
  integer, dimension(G%IsdB:,G%jsd:),  &
                           intent(in) :: OBC_u !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at u-points,
                                             !! with a value of 0 for no OBC, 1 for an Eastern OBC
                                             !! or -1 for a Western OBC
  integer, dimension(G%isd:,G%JsdB:), &
                           intent(in) :: OBC_v !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at v-points,
                                             !! with a value of 0 for no OBC, 1 for a Northern OBC
                                             !! or -1 for a Southern OBC
  logical,                intent(in)  :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in)  :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(isdf:,jsdf:,:), &
                          intent(in)  :: field !< The diagnostic field to be remapped [A]
  real, dimension(isdf:,jsdf:,:), &
                          intent(out) :: interpolated_field !< Field argument remapped to alternative coordinate [A]
  real, dimension(isdf:,jsdf:), &
                optional, intent(in)  :: mask !< A mask for the field [nondim]

  ! Local variables

end subroutine vertically_interpolate_field
module subroutine horizontally_average_diag_field(G, GV, h, staggered_in_x, staggered_in_y, &
                                           is_layer, is_extensive, &
                                           field, averaged_field, &
                                           averaged_mask)
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean vertical grid structure
  real, dimension(:,:,:), intent(in) :: h !< The current thicknesses [H ~> m or kg m-2]
  logical,                intent(in) :: staggered_in_x !< True if the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True if the y-axis location is at v or q points
  logical,                intent(in) :: is_layer !< True if the z-axis location is at h points
  logical,                intent(in) :: is_extensive !< True if the z-direction is spatially integrated (over layers)
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped [A]
  real, dimension(:),    intent(out) :: averaged_field !< Field argument horizontally averaged [A]
  logical, dimension(:), intent(out) :: averaged_mask  !< Mask for horizontally averaged field [nondim]

  ! Local variables

end subroutine horizontally_average_diag_field
module subroutine horizontally_average_field(G, GV, isdf, jsdf, h, staggered_in_x, staggered_in_y, &
                                      is_layer, is_extensive, field, averaged_field, averaged_mask)
  type(ocean_grid_type),   intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean vertical grid structure
  integer,                 intent(in)  :: isdf !< The starting i-index in memory for field
  integer,                 intent(in)  :: jsdf !< The starting j-index in memory for field
  real, dimension(G%isd:,G%jsd:,:), &
                           intent(in)  :: h  !< The current thicknesses [H ~> m or kg m-2]
  logical,                 intent(in)  :: staggered_in_x !< True if the x-axis location is at u or q points
  logical,                 intent(in)  :: staggered_in_y !< True if the y-axis location is at v or q points
  logical,                 intent(in)  :: is_layer !< True if the z-axis location is at h points
  logical,                 intent(in)  :: is_extensive !< True if the z-direction is spatially integrated (over layers)
  real, dimension(isdf:,jsdf:,:), &
                           intent(in)  :: field !<  The diagnostic field to be remapped [A]
  real, dimension(:),      intent(out) :: averaged_field !< Field argument horizontally averaged [A]
  logical, dimension(:),   intent(out) :: averaged_mask  !< Mask for horizontally averaged field [nondim]

  ! Local variables
                                             ! or mass [L2 kg m-2 ~> kg] of each cell.
                                             ! field being averaged in each cell, in [L2 a ~> m2 A],
                                             ! [L2 m a ~> m3 A] or [L2 kg m-2 A ~> kg A],
                                             ! depending on the weighting for the averages and whether the
                                             ! model makes the Boussinesq approximation.
                                               ! in the cells that used in the weighted averages.
                                               ! [A m2], [A m3] or [A kg]

end subroutine horizontally_average_field
  end interface

end module MOM_diag_remap
