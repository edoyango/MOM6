! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initializes fixed aspects of the model, such as horizontal grid metrics,
!! topography and Coriolis.
module MOM_fixed_initialization

use MOM_debugging, only : hchksum, qchksum, uvchksum
use MOM_domains, only : pass_var
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_io, only : slasher
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_open_boundary, only : ocean_OBC_type
use MOM_open_boundary, only : open_boundary_config, open_boundary_query
use MOM_open_boundary, only : open_boundary_impose_normal_slope
use MOM_open_boundary, only : open_boundary_impose_land_mask
use MOM_shared_initialization, only : MOM_initialize_rotation, MOM_calculate_grad_Coriolis
use MOM_shared_initialization, only : initialize_topography_from_file, apply_topography_edits_from_file
use MOM_shared_initialization, only : initialize_topography_named, limit_topography, diagnoseMaximumDepth
use MOM_shared_initialization, only : set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
use MOM_shared_initialization, only : reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
use MOM_shared_initialization, only : read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
use MOM_shared_initialization, only : set_subgrid_topo_at_vel_from_file
use MOM_shared_initialization, only : compute_global_grid_integrals
use MOM_shared_initialization, only : set_meanSL_from_file
use MOM_unit_scaling, only : unit_scale_type

use user_initialization, only : user_initialize_topography
use DOME_initialization, only : DOME_initialize_topography
use ISOMIP_initialization, only : ISOMIP_initialize_topography
use basin_builder, only : basin_builder_topography
use benchmark_initialization, only : benchmark_initialize_topography
use Neverworld_initialization, only : Neverworld_initialize_topography
use DOME2d_initialization, only : DOME2d_initialize_topography
use Kelvin_initialization, only : Kelvin_initialize_topography
use sloshing_initialization, only : sloshing_initialize_topography
use seamount_initialization, only : seamount_initialize_topography
use dumbbell_initialization, only : dumbbell_initialize_topography
use shelfwave_initialization, only : shelfwave_initialize_topography
use Phillips_initialization, only : Phillips_initialize_topography
use dense_water_initialization, only : dense_water_initialize_topography

implicit none ; private

public MOM_initialize_fixed, MOM_initialize_rotation, MOM_initialize_topography


  interface
module subroutine MOM_initialize_fixed(G, US, OBC, PF)
  type(dyn_horgrid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundary structure.
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine MOM_initialize_fixed
module subroutine MOM_initialize_topography(D, max_depth, G, PF, US, meanSL)
  type(dyn_horgrid_type),           intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D  !< Ocean bottom depth [Z ~> m]
  type(param_file_type),            intent(in)  :: PF !< Parameter file structure
  real,                             intent(out) :: max_depth !< Maximum depth or geometric thickness,
                                                             !! with meanSL present, of model [Z ~> m]
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                          optional, intent(in)  :: meanSL !< Mean sea level [Z ~> m]

  ! This subroutine makes the appropriate call to set up the bottom depth.
  ! This is a separate subroutine so that it can be made public and shared with
  ! the ice-sheet code or other components.

  ! Local variables
                                        ! to meanSL. A temporary field used to diagnose maximum
                                        ! static column thickness. D_meanSL = D + meanSL [Z ~> m].

end subroutine MOM_initialize_topography
  end interface

end module MOM_fixed_initialization
