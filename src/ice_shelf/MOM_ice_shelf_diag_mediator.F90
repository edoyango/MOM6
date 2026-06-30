! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The subroutines here provide convenient wrappers to the FMS diag_manager
!! interfaces with additional diagnostic capabilities.
module MOM_IS_diag_mediator

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,          only : chksum0, hchksum, uchksum, vchksum, Bchksum
use MOM_coms,               only : PE_here
use MOM_cpu_clock,          only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,          only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_manager_infra, only : MOM_diag_manager_init
use MOM_diag_manager_infra, only : MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : send_data_infra, MOM_diag_field_add_attribute, EAST, NORTH
use MOM_diag_manager_infra, only : register_diag_field_infra, register_static_field_infra
use MOM_diag_manager_infra, only : get_MOM_diag_field_id, DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : diag_send_complete_infra
use MOM_error_handler,      only : MOM_error, FATAL, is_root_pe, assert, callTree_showQuery
use MOM_error_handler,      only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,        only : get_param, log_version, param_file_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : get_filename_appendix
use MOM_safe_alloc,         only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions,   only : lowercase, uppercase, slasher, ints_to_string, trim_trailing_commas
use MOM_time_manager,       only : time_type, get_time
use MOM_unit_scaling,       only : unit_scale_type

implicit none ; private

public MOM_IS_diag_mediator_infrastructure_init
public MOM_IS_diag_mediator_init, MOM_IS_diag_mediator_end, set_IS_diag_mediator_grid
public set_IS_axes_info, MOM_diag_axis_init
public register_MOM_IS_diag_field, register_MOM_IS_static_field, register_MOM_IS_scalar_field
public post_IS_data, post_IS_data_0d, MOM_IS_diag_send_complete
public safe_alloc_ptr, safe_alloc_alloc, time_type
public enable_averaging, enable_averages, disable_averaging, query_averaging_enabled
public MOM_IS_diag_mediator_close_registration, get_diag_time_end
public define_axes_group, diag_masks_set
public diag_register_area_ids, found_in_diagtable

!> Make a diagnostic available for averaging or output.
interface post_IS_data
  module procedure post_IS_data_2d, post_IS_data_0d
end interface post_IS_data

!> Registers a non-array scalar diagnostic, returning an integer handle
interface register_MOM_IS_scalar_field
  module procedure register_scalar_field_CS, register_scalar_field_axes
end interface register_MOM_IS_scalar_field

!> A group of 1D axes that comprise a 1D/2D/3D mesh
type, public :: axes_grp
  character(len=15) :: id   !< The id string for this particular combination of handles.
  integer           :: rank !< Number of dimensions in the list of axes.
  integer, dimension(:), allocatable :: handles !< Handles to 1D axes.
  type(diag_ctrl), pointer :: diag_cs => null() !< Circular link back to the main diagnostics control structure
                                                !! (Used to avoid passing said structure into every possible call).
  ! ID's for cell_methods
  character(len=9) :: x_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes x-direction.
  character(len=9) :: y_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes y-direction.
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.

  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes [nondim]
  real, pointer, dimension(:,:)   :: mask2d_comp => null() !< Mask for 2-d axes on the computational
                                                      !! domain for this diagnostic [nondim]
end type axes_grp

!> This type is used to represent a diagnostic at the diag_mediator level.
!!
!! There can be both 'primary' and 'secondary' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use              !< True if this entry is being used.
  integer :: fms_diag_id         !< Underlying FMS diag_manager id.
  character(len=64) :: debug_str = '' !< The diagnostic name and module for FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< If non-zero, a factor to multiply data by before posting to FMS,
                                 !! often including factors to undo internal scaling in units of [a A-1 ~> 1]
end type diag_type

!>   The diag_ctrl data type contains times to regulate diagnostics along with masks and
!! axes to use with diagnostics, and a list of structures with data about each diagnostic.
type, public :: diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_iounit = -1           !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum  !< If true, log chksums in a text file instead of posting diagnostics
  logical :: show_call_tree  !< Display the call tree while running. Set by VERBOSITY level.
  logical :: index_space_axes !< If true, diagnostic horizontal coordinates axes are in index space.

  ! The following fields are used for the output of the data.
  ! These give the computational-domain sizes, and are relative to a start value
  ! of 1 in memory for the tracer-point arrays.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain
  ! These give the memory-domain sizes, and can be start at any value on each PE.
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval for any fields
                                !! that are offered for averaging [s].
  type(time_type) :: time_end   !< The end time of the valid interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names indicate
  !! the horizontal locations (B, T, Cu, or Cv) and vertical locations (here just 1).
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  !>@}
  type(axes_grp) :: axesNull !< An axis group for scalars

  ! Mask arrays for 2D diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points [nondim]
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points [nondim]
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dT_comp => null() !< 2D cell-center mask on the computational domain [nondim]

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 15
  type(diag_type), dimension(:), allocatable :: diags !< The list of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID

  !> default missing value to be sent to ALL diagnostics registrations [various]
  real :: missing_value = -1.0e34

  type(ocean_grid_type), pointer :: G => null()  !< The ocean grid type
  type(unit_scale_type), pointer :: US => null() !< A dimensional unit scaling type

  !> Number of checksum-only diagnostics
  integer :: num_chksum_diags

end type diag_ctrl

!>@{ CPU clocks
integer :: id_clock_diag_mediator
!>@}


  interface
module subroutine set_IS_axes_info(G, diag_cs, axes_set_name)
  type(ocean_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(diag_ctrl),       intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in) :: axes_set_name !<  A name to use for this set of axes.
                                                !! The default is "ice".
!   This subroutine sets up the grid and axis information for use by the ice shelf model.

  ! Local variables

end subroutine set_IS_axes_info
module subroutine diag_register_area_ids(diag_cs, id_area_t, id_area_q)
  type(diag_ctrl), intent(inout) :: diag_cs   !< Diagnostics control structure
  integer,   optional, intent(in)    :: id_area_t !< Diag_mediator id for area of h-cells
  integer,   optional, intent(in)    :: id_area_q !< Diag_mediator id for area of q-cells
  ! Local variables
end subroutine diag_register_area_ids
module subroutine define_axes_group(diag_cs, handles, axes, &
                             x_cell_method, y_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Structure used to regulate diagnostic output
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles that define the axis group
  type(axes_grp),             intent(out) :: axes    !< The group of axes that is set up here
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields

  ! Local variables

end subroutine define_axes_group
module subroutine set_IS_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(diag_ctrl),     intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

end subroutine set_IS_diag_mediator_grid
module subroutine post_IS_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_MOM_IS_diag_field.
  real,              intent(in) :: field         !< real value being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables


end subroutine post_IS_data_0d
module subroutine post_IS_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_MOM_IS_diag_field.
  real,      target, intent(in) :: field(:,:)    !< 2-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables

end subroutine post_IS_data_2d
module subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:) !< 2-d array being offered for output or averaging
                                              !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS   !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real, optional, target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables

end subroutine post_data_2d_low
module subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,            intent(in)    :: time_int_in !< The time interval [s] over which any
                                                !! values that are offered are valid.
  type(time_type), intent(in)    :: time_end_in !< The end time of the valid interval
  type(diag_ctrl), intent(inout) :: diag_cs     !< Structure used to regulate diagnostic output
  ! This subroutine enables the accumulation of time averages over the specified time interval.

!  if (num_file==0) return
end subroutine enable_averaging
module subroutine enable_averages(time_int, time_end, diag_CS, T_to_s)
  real,            intent(in)    :: time_int !< The time interval over which any values
                                             !! that are offered are valid [T ~> s].
  type(time_type), intent(in)    :: time_end !< The end time of the valid interval.
  type(diag_ctrl), intent(inout) :: diag_CS  !< A structure that is used to regulate diagnostic output
  real,  optional, intent(in)    :: T_to_s   !< A conversion factor for time_int to seconds [s T-1 ~> 1].
  ! This subroutine enables the accumulation of time averages over the specified time interval.

end subroutine enable_averages
module subroutine disable_averaging(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

end subroutine disable_averaging
logical module function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_cs  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag_cs%time_int [s]
  type(time_type), optional, intent(out) :: time_end !< Current setting of diag_cs%time_end

end function query_averaging_enabled
module subroutine MOM_IS_diag_mediator_infrastructure_init(err_msg)
  character(len=*), optional, intent(out)   :: err_msg !< An error message

end subroutine MOM_IS_diag_mediator_infrastructure_init
module function get_diag_time_end(diag_cs)
  type(diag_ctrl), intent(in)  :: diag_cs !< Structure used to regulate diagnostic output
  type(time_type) :: get_diag_time_end
  !   This function returns the valid end time for diagnostics that are handled
  ! outside of the MOM6 infrastructure, such as via the generic tracer code.

end function get_diag_time_end
module function register_MOM_IS_diag_field(module_name, field_name, axes_in, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion) result (register_diag_field)
  integer :: register_diag_field  !< The returned diagnostic handle
  character(len=*),           intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*),           intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),     target, intent(in) :: axes_in   !< Container with up to 3 integer handles that
                                                      !! indicates axes for this field
  type(time_type),            intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_IS_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                    !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                          !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count    !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to
                                                         !! have no attribute.  If present, this overrides the
                                                         !! default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables

end function register_MOM_IS_diag_field
logical module function register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion)
  integer,          intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model" or "ice_model_fast"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),   intent(in) :: axes !< Container with up to 3 integer handles that indicates axes
                                             !! for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]
  ! Local variables

end function register_diag_field_expand_cmor
integer module function register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,  &
            verbose, do_not_log, err_msg, interp_method, tile_count)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that indicates
                                             !! axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something
                                                       !! (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  ! Local variables

  ! This gets the cell area associated with the grid location of this variable
end function register_diag_field_expand_axes
module subroutine add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
  type(diag_ctrl),        pointer       :: diag_cs !< Diagnostics mediator control structure
  integer,                intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  integer,                intent(in)    :: fms_id !< The FMS diag_manager ID for this diagnostic
  type(diag_type),        pointer       :: this_diag !< This diagnostic
  type(axes_grp), target, intent(in)    :: axes !< Container with up to 3 integer handles that
                                                !! indicates axes for this field
  character(len=*),       intent(in)    :: module_name !< Name of this module, usually
                                                       !! "ocean_model" or "ice_shelf_model"
  character(len=*),       intent(in)    :: field_name !< Name of diagnostic

  ! If the diagnostic is needed obtain a diag_mediator ID (if needed)
end subroutine add_diag_to_list
module subroutine attach_cell_methods(id, axes, ostring, cell_methods, x_cell_method, y_cell_method)
  integer,                    intent(in)  :: id !< Handle to diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container with up to 3 integer handles that indicates
                                                  !! axes for this field
  character(len=*),           intent(out) :: ostring !< The cell_methods strings that would appear in the file
  character(len=*), optional, intent(in)  :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  ! Local variables

end subroutine attach_cell_methods
module function register_scalar_field_axes(module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion) result (register_scalar_field)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that
                                             !! indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

end function register_scalar_field_axes
module function register_scalar_field_CS(module_name, field_name, init_time, diag_cs, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion) result (register_scalar_field)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  type(diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables

end function register_scalar_field_CS
module function register_MOM_IS_static_field(module_name, field_name, axes, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            do_not_log, interp_method, tile_count, &
            cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
            x_cell_method, y_cell_method, area_cell_method, conversion) result(register_static_field)
  integer :: register_static_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that
                                             !! indicates axes for this field
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_IS_data calls (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  integer,          optional, intent(in) :: area !< fms_id for area_t
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
  character(len=*), optional, intent(in) :: area_cell_method !< Specifies the cell method for area
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables

end function register_MOM_IS_static_field
module subroutine describe_option(opt_name, value, diag_CS)
  character(len=*),    intent(in) :: opt_name !< The name of the option
  character(len=*),    intent(in) :: value    !< The value of the option
  type(diag_ctrl), intent(in) :: diag_CS  !< Structure used to regulate diagnostic output

  ! Local variables

end subroutine describe_option
module subroutine MOM_IS_diag_mediator_init(G, US, param_file, diag_cs, component, err_msg, &
                                  doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The horizontal grid type
  type(unit_scale_type), target, intent(in) :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),            intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in)    :: component !< An optional component name
  character(len=*), optional, intent(out)   :: err_msg !< A string for a returned error message
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine MOM_IS_diag_mediator_init
module subroutine diag_masks_set(G, missing_value, diag_cs)
  type(ocean_grid_type), target, intent(in)    :: G   !< The horizontal grid type
  real,                          intent(in)    :: missing_value !< A fill value for missing points
  type(diag_ctrl),               intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

  ! Local variables

  ! 2d masks point to the model masks since they are identical
end subroutine diag_masks_set
module subroutine MOM_IS_diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

end subroutine MOM_IS_diag_mediator_close_registration
module subroutine MOM_IS_diag_mediator_end(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  ! Local variables

end subroutine MOM_IS_diag_mediator_end
integer module function get_new_diag_id(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< Diagnostics control structure
  ! Local variables

end function get_new_diag_id
module subroutine initialize_diag_type(diag)
  type(diag_type), intent(inout) :: diag !< diag_type to be initialized

end subroutine initialize_diag_type
module subroutine alloc_diag_with_id(diag_id, diag_cs, diag)
  integer,                 intent(in   ) :: diag_id !< id for the diagnostic
  type(diag_ctrl), target, intent(inout) :: diag_cs !< structure used to regulate diagnostic output
  type(diag_type),         pointer       :: diag    !< structure representing a diagnostic (inout)


end subroutine alloc_diag_with_id
module subroutine log_available_diag(used, module_name, field_name, cell_methods_string, comment, &
                              diag_CS, long_name, units, standard_name, variants, dimensions)
  logical,          intent(in) :: used !< Whether this diagnostic was in the diag_table or not
  character(len=*), intent(in) :: module_name !< Name of the diagnostic module
  character(len=*), intent(in) :: field_name !< Name of this diagnostic field
  character(len=*), intent(in) :: cell_methods_string !< The spatial component of the CF cell_methods attribute
  character(len=*), intent(in) :: comment !< A comment to append after [Used|Unused]
  type(diag_ctrl), intent(in) :: diag_CS  !< The diagnotics control structure
  character(len=*), optional, intent(in) :: dimensions !< Descriptor of the horizontal and vertical dimensions
  character(len=*), optional, intent(in) :: long_name !< CF long name of diagnostic
  character(len=*), optional, intent(in) :: units !< Units for diagnostic
  character(len=*), optional, intent(in) :: standard_name !< CF standardized name of diagnostic
  character(len=*), optional, intent(in) :: variants !< Alternate modules and variable names for
                                                     !! this diagnostic and derived diagnostics
  ! Local variables

end subroutine log_available_diag
module subroutine log_chksum_diag(docunit, description, chksum)
  integer,          intent(in) :: docunit     !< Handle of the log file
  character(len=*), intent(in) :: description !< Name of the diagnostic module
  integer,          intent(in) :: chksum      !< chksum of the diagnostic

end subroutine log_chksum_diag
logical module function found_in_diagtable(diag, varName)
  type(diag_ctrl), intent(in) :: diag     !< A structure used to control diagnostics.
  character(len=*),    intent(in) :: varName  !< The obsolete diagnostic name
  ! Local

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
end function found_in_diagtable
module subroutine MOM_IS_diag_send_complete()
end subroutine MOM_IS_diag_send_complete
  end interface

end module MOM_IS_diag_mediator
