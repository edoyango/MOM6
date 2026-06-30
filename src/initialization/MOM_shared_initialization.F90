! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Code that initializes fixed aspects of the model grid, such as horizontal
!! grid metrics, topography and Coriolis, and can be shared between components.
module MOM_shared_initialization

use MOM_coms, only : max_across_PEs, reproducing_sum
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, param_file_type, log_version
use MOM_io, only : create_MOM_file, file_exists, field_size, get_filename_appendix
use MOM_io, only : MOM_infra_file, MOM_field
use MOM_io, only : MOM_read_data, MOM_read_vector, read_variable, stdout
use MOM_io, only : open_file_to_read, close_file_to_read, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, MOM_write_field, var_desc
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

public MOM_shared_init_init
public MOM_initialize_rotation, MOM_calculate_grad_Coriolis
public initialize_topography_from_file, apply_topography_edits_from_file
public initialize_topography_named, limit_topography, diagnoseMaximumDepth
public set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
public reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
public read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
public set_subgrid_topo_at_vel_from_file
public compute_global_grid_integrals, write_ocean_geometry_file
public set_meanSL_from_file

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine MOM_shared_init_init(PF)
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.


! This include declares and sets the variable "version".
end subroutine MOM_shared_init_init
module subroutine MOM_initialize_rotation(f, G, PF, US)
  type(dyn_horgrid_type),                       intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f  !< The Coriolis parameter [T-1 ~> s-1]
  type(param_file_type),                        intent(in)  :: PF !< Parameter file structure
  type(unit_scale_type),                        intent(in)  :: US !< A dimensional unit scaling type

!   This subroutine makes the appropriate call to set up the Coriolis parameter.
! This is a separate subroutine so that it can be made public and shared with
! the ice-sheet code or other components.
! Set up the Coriolis parameter, f, either analytically or from file.

end subroutine MOM_initialize_rotation
module subroutine MOM_calculate_grad_Coriolis(dF_dx, dF_dy, G, US)
  type(dyn_horgrid_type),             intent(inout) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dx !< x-component of grad f [T-1 L-1 ~> s-1 m-1]
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dy !< y-component of grad f [T-1 L-1 ~> s-1 m-1]
  type(unit_scale_type),    optional, intent(in)    :: US !< A dimensional unit scaling type
  ! Local variables

end subroutine MOM_calculate_grad_Coriolis
module function diagnoseMaximumDepth(D, G)
  type(dyn_horgrid_type),  intent(in) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(in) :: D !< Ocean bottom depth in [m] or [Z ~> m]
  real :: diagnoseMaximumDepth             !< The global maximum ocean bottom depth in [m] or [Z ~> m]
  ! Local variables
end function diagnoseMaximumDepth
module subroutine set_meanSL_from_file(meanSL, G, param_file, US)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: meanSL !< Mean sea level referenced to a zero
                                                          !! reference height at tracer points [Z ~> m].
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  ! Local variables

end subroutine set_meanSL_from_file
module subroutine initialize_topography_from_file(D, G, param_file, US)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  ! Local variables

end subroutine initialize_topography_from_file
module subroutine apply_topography_edits_from_file(D, G, param_file, US)
  type(dyn_horgrid_type),           intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(inout) :: D !< Ocean bottom depth [m] or [Z ~> m] if
                                                       !! US is present
  type(param_file_type),            intent(in)    :: param_file !< Parameter file structure
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine apply_topography_edits_from_file
module subroutine initialize_topography_named(D, G, param_file, topog_config, max_depth, US)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  character(len=*),                 intent(in)  :: topog_config !< The name of an idealized
                                                              !! topographic configuration
  real,                             intent(in)  :: max_depth  !< Maximum depth [Z ~> m]
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type

  ! This subroutine places the bottom depth in m into D(:,:), shaped according to the named config.

  ! Local variables
end subroutine initialize_topography_named
module subroutine limit_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type), intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                          intent(inout) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),  intent(in)    :: param_file !< Parameter file structure
  real,                   intent(in)    :: max_depth  !< Maximum depth of model [Z ~> m]
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type

  ! Local variables

end subroutine limit_topography
module subroutine set_rotation_planetary(f, G, param_file, US)
  type(dyn_horgrid_type), intent(in)  :: G  !< The dynamic horizontal grid
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(out) :: f  !< Coriolis parameter (vertical component) [T-1 ~> s-1]
  type(param_file_type),  intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type),  intent(in)  :: US !< A dimensional unit scaling type

! This subroutine sets up the Coriolis parameter for a sphere

end subroutine set_rotation_planetary
module subroutine set_rotation_beta_plane(f, G, param_file, US)
  type(dyn_horgrid_type), intent(in)  :: G  !< The dynamic horizontal grid
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(out) :: f  !< Coriolis parameter (vertical component) [T-1 ~> s-1]
  type(param_file_type),  intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type),  intent(in)  :: US !< A dimensional unit scaling type

! This subroutine sets up the Coriolis parameter for a beta-plane

end subroutine set_rotation_beta_plane
module subroutine initialize_grid_rotation_angle(G, PF)
  type(dyn_horgrid_type), intent(inout) :: G   !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: PF  !< A structure indicating the open file
                                               !! to parse for model parameter values.

                        ! to equivalent distances in latitudes [nondim]

end subroutine initialize_grid_rotation_angle
module function modulo_around_point(x, xc, Lx) result(x_mod)
  real, intent(in) :: x  !< Value to which to apply modulo arithmetic [A]
  real, intent(in) :: xc !< Center of modulo range [A]
  real, intent(in) :: Lx !< Modulo range width [A]
  real :: x_mod          !< x shifted by an integer multiple of Lx to be close to xc [A].

end function modulo_around_point
module subroutine reset_face_lengths_named(G, param_file, name, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),       intent(in)    :: name !< The name for the set of face lengths. Only "global_1deg"
                                                !! is currently implemented.
  type(unit_scale_type),  intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables
end subroutine reset_face_lengths_named
module subroutine reset_face_lengths_file(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type),  intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine reset_face_lengths_file
module subroutine reset_face_lengths_list(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type),  intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables
                          ! +/- 360 degrees from the specified range of values.

end subroutine reset_face_lengths_list
module subroutine read_face_length_list(iounit, filename, num_lines, lines)
  integer,                          intent(in)  :: iounit    !< An open I/O unit number for the file
  character(len=*),                 intent(in)  :: filename  !< The name of the face-length file to read
  integer,                          intent(out) :: num_lines !< The number of non-blank lines in the file
  character(len=120), dimension(:), pointer     :: lines  !< The non-blank lines, after removing comments

  !   This subroutine reads and counts the non-blank lines in the face length
  ! list file, after removing comments.

end subroutine read_face_length_list
module subroutine set_subgrid_topo_at_vel_from_file(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G          !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: param_file !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US         !< A dimensional unit scaling type

  ! Local variables

end subroutine set_subgrid_topo_at_vel_from_file
module subroutine set_velocity_depth_max(G)
  type(dyn_horgrid_type), intent(inout) :: G   !< The dynamic horizontal grid
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! maximum of the adjacent depths.

end subroutine set_velocity_depth_max
module subroutine set_velocity_depth_min(G)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! minimum of the adjacent depths.

end subroutine set_velocity_depth_min
module subroutine compute_global_grid_integrals(G, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(unit_scale_type),  intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine compute_global_grid_integrals
module subroutine write_ocean_geometry_file(G, param_file, directory, US, geom_file)
  type(dyn_horgrid_type),       intent(inout) :: G         !< The dynamic horizontal grid
  type(param_file_type),        intent(in)    :: param_file !< Parameter file structure
  character(len=*),             intent(in)    :: directory !< The directory into which to place the geometry file.
  type(unit_scale_type),        intent(in)    :: US        !< A dimensional unit scaling type
  character(len=*),   optional, intent(in)    :: geom_file !< If present, the name of the geometry file
                                                           !! (otherwise the file is "ocean_geometry")

  ! Local variables.

end subroutine write_ocean_geometry_file
  end interface

end module MOM_shared_initialization
