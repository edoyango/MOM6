! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The subroutines here provide convenient wrappers to the FMS diag_manager
!! interfaces with additional diagnostic capabilities.
module MOM_diag_mediator

use MOM_array_transform,  only : symmetric_sum
use MOM_checksums,        only : chksum0, zchksum, hchksum, uchksum, vchksum, Bchksum
use MOM_coms,             only : PE_here
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_buffers,     only : diag_buffer_2d, diag_buffer_3d
use MOM_diag_manager_infra, only : MOM_diag_manager_init, MOM_diag_manager_end
use MOM_diag_manager_infra, only : diag_axis_init=>MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : send_data_infra, MOM_diag_field_add_attribute, EAST, NORTH
use MOM_diag_manager_infra, only : register_diag_field_infra, register_static_field_infra
use MOM_diag_manager_infra, only : get_MOM_diag_field_id, DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : diag_send_complete_infra
use MOM_diag_remap,       only : diag_remap_ctrl, diag_remap_update, diag_remap_calc_hmask
use MOM_diag_remap,       only : diag_remap_init, diag_remap_end, diag_remap_do_remap
use MOM_diag_remap,       only : vertically_reintegrate_diag_field, vertically_interpolate_diag_field
use MOM_diag_remap,       only : horizontally_average_diag_field, diag_remap_get_axes_info
use MOM_diag_remap,       only : diag_remap_configure_axes, diag_remap_axes_configured
use MOM_diag_remap,       only : diag_remap_diag_registration_closed, diag_remap_set_active
use MOM_EOS,              only : EOS_type
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, is_root_pe, assert, callTree_showQuery
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_io,               only : vardesc, query_vardesc
use MOM_io,               only : get_filename_appendix
use MOM_safe_alloc,       only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase, slasher, ints_to_string, trim_trailing_commas
use MOM_time_manager,     only : time_type, get_time
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_domains,          only : get_domain_extent, clone_MOM_domain

implicit none ; private

#undef __DO_SAFETY_CHECKS__
#define IMPLIES(A, B) ((.not. (A)) .or. (B))

public set_axes_info, post_data, register_diag_field, time_type
public post_data_3d_by_column, post_data_3d_final
public post_product_u, post_product_sum_u, post_product_v, post_product_sum_v
public set_masks_for_axes, MOM_diag_send_complete
! post_data_1d_k is a deprecated interface that can be replaced by a call to post_data, but
! it is being retained for backward compatibility to older versions of the ocean_BGC code.
public post_data_1d_k
public safe_alloc_ptr, safe_alloc_alloc
public enable_averaging, enable_averages, disable_averaging, query_averaging_enabled
public diag_mediator_init, diag_mediator_end, set_diag_mediator_grid
public diag_mediator_infrastructure_init, diag_mediator_set_OBC_info
public diag_mediator_close_registration, get_diag_time_end
public diag_axis_init, ocean_register_diag, register_static_field
public register_scalar_field
public define_axes_group, diag_masks_set
public set_piecemeal_extents
public diag_register_area_ids
public register_cell_measure, diag_associate_volume_cell_measure
public diag_get_volume_cell_measure_dm_id
public diag_set_state_ptrs, diag_update_remap_grids
public diag_grid_storage_init, diag_grid_storage_end
public diag_copy_diag_to_storage, diag_copy_storage_to_diag
public diag_save_grids, diag_restore_grids
public found_in_diagtable

!> Make a diagnostic available for averaging or output.
interface post_data
  module procedure post_data_3d, post_data_2d, post_data_1d_k, post_data_0d
end interface post_data

!> Registers a non-array scalar diagnostic, returning an integer handle
interface register_scalar_field
  module procedure register_scalar_field_CS, register_scalar_field_axes
end interface register_scalar_field

!> Down sample a field
interface downsample_field
  module procedure downsample_field_2d, downsample_field_3d
end interface downsample_field

!> Down sample the mask of a field
interface downsample_mask
  module procedure downsample_mask_2d, downsample_mask_3d
end interface downsample_mask

!> Down sample a diagnostic field
interface downsample_diag_field
  module procedure downsample_diag_field_2d, downsample_diag_field_3d
end interface downsample_diag_field

!> Contained for down sampled masks
type, private :: diag_dsamp
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes [nondim]
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes [nondim]
end type diag_dsamp

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
  character(len=9) :: v_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes vertical direction.
  ! For remapping
  integer :: nz = 0 !< Vertical dimension of diagnostic
  integer :: vertical_coordinate_number = 0 !< Index of the corresponding diag_remap_ctrl for this axis group
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.
  logical :: is_layer = .false. !< If true, indicates that this axes group is for a layer vertically-located field.
  logical :: is_interface = .false. !< If true, indicates that this axes group is for an interface
                                    !! vertically-located field.
  logical :: is_native = .true. !< If true, indicates that this axes group is for a native model grid.
                                !! False for any other grid. Used for rank>2.
  logical :: needs_remapping = .false. !< If true, indicates that this axes group is for a intensive layer-located
                                       !! field that must be remapped to these axes. Used for rank>2.
  logical :: needs_interpolating = .false. !< If true, indicates that this axes group is for a sampled
                                         !! interface-located field that must be interpolated to
                                         !! these axes. Used for rank>2.
  integer :: downsample_level_factor = 1 !< If greater than 1, the factor by which this diagnostic will be downsampled
  integer :: downsample_level_index = 0 !< If greater than 0, the index for the downsample level for this diagnostic
                                        !! in the diag_cs%dsamp array.
  ! For horizontally averaged diagnostics (applies to 2d and 3d fields only)
  type(axes_grp), pointer :: xyave_axes => null() !< The associated 1d axes for horizontally area-averaged diagnostics
  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  integer :: id_volume = -1 !< The diag_manager id for volume to be used for cell_measure of variables
                            !! with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes [nondim]
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes [nondim]
  type(diag_dsamp), dimension(:), allocatable :: dsamp !< Downsample container

  ! For diagnostics posted piecemeal
  type(diag_buffer_2d) :: piecemeal_2d !< A dynamically reallocated buffer for 2d piecemeal diagnostics
  type(diag_buffer_3d) :: piecemeal_3d !< A dynamically reallocated buffer for 3d piecemeal diagnostics
end type axes_grp

!> Contains an array to store a diagnostic target grid
type, private :: diag_grids_type
  real, dimension(:,:,:), allocatable :: h  !< Target grid for remapped coordinate [H ~> m or kg m-2] or [Z ~> m]
end type diag_grids_type

!> Stores all the remapping grids and the model's native space thicknesses
type, public :: diag_grid_storage
  integer                                          :: num_diag_coords !< Number of target coordinates
  real, dimension(:,:,:), allocatable              :: h_state         !< Layer thicknesses in native
                                                                      !! space [H ~> m or kg m-2]
  type(diag_grids_type), dimension(:), allocatable :: diag_grids      !< Primarily empty, except h field
end type diag_grid_storage

! Integers to encode the total cell methods
! Note that vorticity points (the PPP and PPM methods) are not fully dealt with for downsampling.
integer :: PPP=111  !< x:point,y:point,z:point
!integer :: PPS=112 ! x:point,y:point,z:sum  , this kind of diagnostic is not currently present in diag_table.MOM6
integer :: PPM=113  !< x:point,y:point,z:mean
integer :: PSP=121  !< x:point,y:sum,z:point
integer :: PSS=122  !< x:point,y:sum,z:point
integer :: PSM=123  !< x:point,y:sum,z:mean
integer :: PMP=131  !< x:point,y:mean,z:point
integer :: PMM=133  !< x:point,y:mean,z:mean
integer :: SPP=211  !< x:sum,y:point,z:point
integer :: SPS=212  !< x:sum,y:point,z:sum
integer :: SSP=221  !< x:sum,y:sum,z:point
integer :: MPP=311  !< x:mean,y:point,z:point
integer :: MPM=313  !< x:mean,y:point,z:mean
integer :: MMP=331  !< x:mean,y:mean,z:point
integer :: MMS=332  !< x:mean,y:mean,z:sum
integer :: SSS=222  !< x:sum,y:sum,z:sum
integer :: MMM=333  !< x:mean,y:mean,z:mean

!> This type is used to represent a diagnostic at the diag_mediator level.
!!
!! There can be both 'primary' and 'secondary' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use !< True if this entry is being used.
  integer :: fms_diag_id !< Underlying FMS diag_manager id.
  integer :: fms_xyave_diag_id = -1 !< For a horizontally area-averaged diagnostic.
  integer :: downsample_diag_id = -1 !< For a horizontally area-downsampled diagnostic.
  character(len=64) :: debug_str = '' !< The diagnostic name and module for FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< If non-zero, a factor to multiply data by before posting to FMS,
                                 !! often including factors to undo internal scaling in units of [a A-1 ~> 1]
  logical :: v_extensive = .false. !< True for vertically extensive fields (vertically integrated).
                                   !! False for intensive (concentrations).
  integer :: xyz_method = 0 !< A 3 digit integer encoding the diagnostics cell method
                            !! It can be used to determine the downsample algorithm
end type diag_type

!> Container for down sampling information
type diagcs_dsamp
  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain
  integer :: isgB !< The start i-index of cell corners within the global domain
  integer :: iegB !< The end i-index of cell corners within the global domain
  integer :: jsgB !< The start j-index of cell corners within the global domain
  integer :: jegB !< The end j-index of cell corners within the global domain

  !>@{ Axes for each location on a diagnostic grid
  type(axes_grp)  :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp)  :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp)  :: axesB1, axesT1, axesCu1, axesCv1
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi
  !>@}

  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points [nondim]
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points [nondim]
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points [nondim]
  !>@{ 3D mask arrays for diagnostics at layers (mask...L) and interfaces (mask...i), all [nondim]
  real, dimension(:,:,:), pointer :: mask3dTL  => null()
  real, dimension(:,:,:), pointer :: mask3dBL  => null()
  real, dimension(:,:,:), pointer :: mask3dCuL => null()
  real, dimension(:,:,:), pointer :: mask3dCvL => null()
  real, dimension(:,:,:), pointer :: mask3dTi  => null()
  real, dimension(:,:,:), pointer :: mask3dBi  => null()
  real, dimension(:,:,:), pointer :: mask3dCui => null()
  real, dimension(:,:,:), pointer :: mask3dCvi => null()
  !>@}
end type diagcs_dsamp

!> The following data type a list of diagnostic fields an their variants,
!! as well as variables that control the handling of model output.
type, public :: diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_iounit = -1           !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum  !< If true, log chksums in a text file instead of posting diagnostics
  logical :: show_call_tree  !< Display the call tree while running. Set by VERBOSITY level.
  logical :: index_space_axes !< If true, diagnostic horizontal coordinates axes are in index space.
  logical :: symmetric_downsample_sums !< If true, use rotationally symmetric sums when downsampling diagnostics.

  ! The following fields are used for the output of the data.
  ! These give the computational-domain sizes, and are relative to a start value
  ! of 1 in memory for the tracer-point arrays.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain
  ! These give the memory-domain sizes, and can start at any value on each PE.
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval for any fields
                                !! that are offered for averaging [s].
  type(time_type) :: time_end   !< The end time of the valid interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names
  !! indicate the horizontal (B, T, Cu, or Cv) and vertical (L, i, or 1) locations.
  type(axes_grp) :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp) :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  !>@}
  type(axes_grp) :: axesZi !< A 1-D z-space axis at interfaces
  type(axes_grp) :: axesZL !< A 1-D z-space axis at layer centers
  type(axes_grp) :: axesNull !< An axis group for scalars

  ! Mask arrays for 2D diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points [nondim]
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points [nondim]
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points [nondim]
  !>@{ 3D mask arrays for diagnostics at layers (mask...L) and interfaces (mask...i) all [nondim]
  real, dimension(:,:,:), pointer :: mask3dTL  => null()
  real, dimension(:,:,:), pointer :: mask3dBL  => null()
  real, dimension(:,:,:), pointer :: mask3dCuL => null()
  real, dimension(:,:,:), pointer :: mask3dCvL => null()
  real, dimension(:,:,:), pointer :: mask3dTi  => null()
  real, dimension(:,:,:), pointer :: mask3dBi  => null()
  real, dimension(:,:,:), pointer :: mask3dCui => null()
  real, dimension(:,:,:), pointer :: mask3dCvi => null()

  integer :: num_diag_dsamp_levels !< The number of downsampled levels requested in the parameters files (default 0)
  integer, dimension(:), allocatable :: diag_dsamp_levels !< The downsample levels requested by diag registrations
  type(diagcs_dsamp), dimension(:), allocatable :: dsamp !< An array of downsampling control containers
                                                         !! for each level of downsampling that is being used,
                                                         !! with a size determined at runtime via NUM_DIAG_DOWNSAMP_LEV
  !>@}

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 100
  type(diag_type), dimension(:), allocatable :: diags !< The list of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID

  !> default missing value to be sent to ALL diagnostics registrations [various]
  real :: missing_value = -1.0e34

  !> Number of diagnostic vertical coordinates (remapped)
  integer :: num_diag_coords
  !> Control structure for each possible coordinate
  type(diag_remap_ctrl), dimension(:), allocatable :: diag_remap_cs
  type(diag_grid_storage) :: diag_grid_temp !< Stores the remapped diagnostic grid
  logical :: diag_grid_overridden = .false. !< True if the diagnostic grids have been overriden

  type(axes_grp), dimension(:), allocatable :: &
    remap_axesZL, &  !< The 1-D z-space cell-centered axis for remapping
    remap_axesZi     !< The 1-D z-space interface axis for remapping
  !>@{ Axes used for remapping
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi
  !>@}

  ! Pointer to H, G and T&S needed for remapping
  real, dimension(:,:,:), pointer :: h => null() !< The thicknesses needed for remapping [H ~> m or kg m-2]
  real, dimension(:,:,:), pointer :: T => null() !< The temperatures needed for remapping [C ~> degC]
  real, dimension(:,:,:), pointer :: S => null() !< The salinities needed for remapping [S ~> ppt]
  type(EOS_type),         pointer :: eqn_of_state => null() !< The equation of state type
  type(thermo_var_ptrs),  pointer :: tv => null()   !< A structure with thermodynamic variables that are
                                                    !! used to convert thicknesses to vertical extents
  type(ocean_grid_type), pointer :: G => null()  !< The ocean grid type
  type(verticalGrid_type), pointer :: GV => null()  !< The model's vertical ocean grid
  type(unit_scale_type), pointer :: US => null() !< A dimensional unit scaling type

  !> The volume cell measure (special diagnostic) manager id
  integer :: volume_cell_measure_dm_id = -1

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of h so that we know whether it has changed [H ~> m or kg m-2]. If it has then
  ! need the target grid for vertical remapping needs to have been updated.
  real, dimension(:,:,:), allocatable :: h_old
#endif

  !> Number of checksum-only diagnostics
  integer :: num_chksum_diags

  integer, dimension(:,:), allocatable :: OBC_u !< An array that indicates the presence and direction
                                                !! of any open boundary conditions at u-points,
                                                !! with a value of 0 for no OBC, 1 for an
                                                !! Eastern OBC or -1 for a Western OBC
  integer, dimension(:,:), allocatable :: OBC_v !< An array that indicates the presence and direction
                                                !! of any open boundary conditions at v-points,
                                                !! with a value of 0 for no OBC, 1 for a Northern OBC
                                                !! or -1 for a Southern OBC
  real, dimension(:,:,:), allocatable :: h_begin !< Layer thicknesses at the beginning of the timestep used
                                                 !! for remapping of extensive variables [H ~> m or kg m-2]
  real, dimension(:,:,:), allocatable :: dz_begin !< Layer vertical extents at the beginning of the timestep used
                                                 !! for remapping of extensive variables [Z ~> m]

end type diag_ctrl

!>@{ CPU clocks
integer :: id_clock_diag_mediator, id_clock_diag_remap, id_clock_diag_grid_updates
!>@}


  interface
module subroutine set_axes_info(G, GV, US, param_file, diag_cs, set_vertical)
  type(ocean_grid_type),   intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),         intent(inout) :: diag_cs !< Diagnostics control structure
  logical,       optional, intent(in)    :: set_vertical !< If true or missing, set up
                                                       !! vertical axes
  ! Local variables
                          ! that might be [m], [kg m-3] or [nondim], depending on the coordinate.
                          ! that might be [m], [kg m-3] or [nondim], depending on the coordinate.

end subroutine set_axes_info
module subroutine set_axes_info_dsamp(G, GV, param_file, diag_cs, id_zl_native, id_zi_native)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),       intent(inout) :: diag_cs !< Diagnostics control structure
  integer,               intent(in)    :: id_zl_native !< ID of native layers
  integer,               intent(in)    :: id_zi_native !< ID of native interfaces

  ! Local variables
                                                         ! the output axes, often in units of [degrees_N] or
                                                         ! [km] or [m] or [gridpoints].
                                                         ! the output axes, often in units of [degrees_N] or
                                                         ! [km] or [m] or [gridpoints].
                                                         ! the output axes, often in units of [degrees_N] or
                                                         ! [km] or [m] or [gridpoints].
                                                         ! the output axes, often in units of [degrees_N] or
                                                         ! [km] or [m] or [gridpoints].


  ! Axes group for native downsampled diagnostics
  !Loop over the downsampling levels requested in parameters.
end subroutine set_axes_info_dsamp
module subroutine set_masks_for_axes(G, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G !< The ocean grid type.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables

end subroutine set_masks_for_axes
module subroutine set_masks_for_axes_dsamp(G, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G !< The ocean grid type.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables

  ! Each downsampled axis needs both downsampled and non-downsampled masks.
  ! The downsampled mask is needed for sending out the diagnostics output via diag_manager.
  ! The non-downsampled mask is needed for downsampling the diagnostics field.
end subroutine set_masks_for_axes_dsamp
module subroutine diag_register_area_ids(diag_cs, id_area_t, id_area_q)
  type(diag_ctrl),   intent(inout) :: diag_cs   !< Diagnostics control structure
  integer, optional, intent(in)    :: id_area_t !< Diag_mediator id for area of h-cells
  integer, optional, intent(in)    :: id_area_q !< Diag_mediator id for area of q-cells
  ! Local variables
end subroutine diag_register_area_ids
module subroutine register_cell_measure(G, diag, Time)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< Regulates diagnostic output
  type(time_type),         intent(in)    :: Time !< Model time
  ! Local variables
end subroutine register_cell_measure
module subroutine diag_associate_volume_cell_measure(diag_cs, id_h_volume)
  type(diag_ctrl),   intent(inout) :: diag_cs     !< Diagnostics control structure
  integer,           intent(in)    :: id_h_volume !< Diag_manager id for volume of h-cells
  ! Local variables

end subroutine diag_associate_volume_cell_measure
integer module function diag_get_volume_cell_measure_dm_id(diag_cs)
  type(diag_ctrl),   intent(in) :: diag_cs   !< Diagnostics control structure

end function diag_get_volume_cell_measure_dm_id
module subroutine define_axes_group(diag_cs, handles, axes, nz, vertical_coordinate_number, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point, &
                             is_layer, is_interface, &
                             is_native, needs_remapping, needs_interpolating, &
                             xyave_axes)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Diagnostics control structure
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles
  type(axes_grp),             intent(out) :: axes    !< The group of 1D axes
  integer,          optional, intent(in)  :: nz      !< Number of layers in this diagnostic grid
  integer,          optional, intent(in)  :: vertical_coordinate_number !< Index number for vertical coordinate
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: v_cell_method !< A vertical direction cell method used to construct
                                                        !! the "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields
  logical,          optional, intent(in)  :: is_layer   !< If true, indicates that this axes group is
                                                        !! for a layer vertically-located field.
  logical,          optional, intent(in)  :: is_interface !< If true, indicates that this axes group
                                                        !! is for an interface vertically-located field.
  logical,          optional, intent(in)  :: is_native  !< If true, indicates that this axes group is
                                                        !! for a native model grid. False for any other grid.
  logical,          optional, intent(in)  :: needs_remapping !< If true, indicates that this axes group is
                                                        !! for a intensive layer-located field that must
                                                        !! be remapped to these axes. Used for rank>2.
  logical,          optional, intent(in)  :: needs_interpolating !< If true, indicates that this axes group
                                                        !! is for a sampled interface-located field that must
                                                        !! be interpolated to these axes. Used for rank>2.
  type(axes_grp),   optional, target      :: xyave_axes !< The corresponding axes group for horizontally
                                                        !! area-average diagnostics
  ! Local variables

end subroutine define_axes_group
module subroutine define_axes_group_dsamp(diag_cs, handles, axes, dl, nz, vertical_coordinate_number, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point, &
                             is_layer, is_interface, &
                             is_native, needs_remapping, needs_interpolating, &
                             xyave_axes)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Diagnostics control structure
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles
  type(axes_grp),             intent(out) :: axes    !< The group of 1D axes
  integer,                    intent(in)  :: dl      !< Downsample level index
  integer,          optional, intent(in)  :: nz      !< Number of layers in this diagnostic grid
  integer,          optional, intent(in)  :: vertical_coordinate_number !< Index number for vertical coordinate
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: v_cell_method !< A vertical direction cell method used to construct
                                                        !! the "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields
  logical,          optional, intent(in)  :: is_layer   !< If true, indicates that this axes group is
                                                        !! for a layer vertically-located field.
  logical,          optional, intent(in)  :: is_interface !< If true, indicates that this axes group
                                                        !! is for an interface vertically-located field.
  logical,          optional, intent(in)  :: is_native  !< If true, indicates that this axes group is
                                                        !! for a native model grid. False for any other grid.
  logical,          optional, intent(in)  :: needs_remapping !< If true, indicates that this axes group is
                                                        !! for a intensive layer-located field that must
                                                        !! be remapped to these axes. Used for rank>2.
  logical,          optional, intent(in)  :: needs_interpolating !< If true, indicates that this axes group
                                                        !! is for a sampled interface-located field that must
                                                        !! be interpolated to these axes. Used for rank>2.
  type(axes_grp),   optional, target      :: xyave_axes !< The corresponding axes group for horizontally
                                                        !! area-average diagnostics
  ! Local variables

end subroutine define_axes_group_dsamp
module subroutine set_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

end subroutine set_diag_mediator_grid
module subroutine post_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field         !< real value being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables


end subroutine post_data_0d
module subroutine post_data_1d_k(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real, target,      intent(in) :: field(:)      !< 1-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables


end subroutine post_data_1d_k
module subroutine post_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field(:,:)    !< 2-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables

end subroutine post_data_2d
module subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:) !< 2-d array being offered for output or averaging
                                              !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real, optional, target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables

end subroutine post_data_2d_low
module subroutine post_data_3d(diag_field_id, field, diag_cs, is_static, mask, alt_h)

  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field(:,:,:)  !< 3-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask [nondim]
  real, dimension(:,:,:), &
         target, optional, intent(in) :: alt_h  !< An alternate thickness to use for vertically
                                                !! remapping this diagnostic [H ~> m or kg m-2].

  ! Local variables
                                                !! remapping this diagnostic [H ~> m or kg m-2].


end subroutine post_data_3d
module subroutine post_data_3d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:,:) !< 3-d array being offered for output or averaging
                                                !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional,target, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables


end subroutine post_data_3d_low
module subroutine post_data_3d_by_column(diag_field_id, field, diag_cs, i, j)
  integer,            intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                  !! previous call to register_diag_field.
  real, dimension(:), intent(in) :: field         !< 3-d array being offered for output or averaging
                                                  !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS  !< Structure used to regulate diagnostic output
  integer,            intent(in) :: i             !< The i-index to post the data in the buffer
  integer,            intent(in) :: j             !< The j-index to post the data in the buffer


end subroutine post_data_3d_by_column
module subroutine post_data_3d_by_point(diag_field_id, field, diag_cs, i, j, k)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field         !< 3-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  integer,           intent(in) :: i             !< The i-index to post the data in the buffer
  integer,           intent(in) :: j             !< The j-index to post the data in the buffer
  integer,           intent(in) :: k             !< The k-index to post the data in the buffer


end subroutine post_data_3d_by_point
module subroutine post_data_3d_final(diag_field_id, diag_cs)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output


end subroutine post_data_3d_final
module subroutine post_product_u(id, u_a, u_b, G, nz, diag, mask, alt_h)
  integer,                  intent(in) :: id   !< The ID for this diagnostic
  type(ocean_grid_type),    intent(in) :: G    !< ocean grid structure
  integer,                  intent(in) :: nz   !< The size of the arrays in the vertical
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, nz), &
                            intent(in) :: u_a  !< The first u-point array in arbitrary units [A]
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, nz), &
                            intent(in) :: u_b  !< The second u-point array in arbitrary units [B]
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output
  real,           optional, intent(in) :: mask(:,:,:)  !< If present, use this real array as the data mask [nondim]
  real,   target, optional, intent(in) :: alt_h(:,:,:) !< An alternate thickness to use for vertically
                                               !! remapping this diagnostic [H ~> m or kg m-2]

  ! Local variables

end subroutine post_product_u
module subroutine post_product_sum_u(id, u_a, u_b, G, nz, diag)
  integer,                  intent(in) :: id   !< The ID for this diagnostic
  type(ocean_grid_type),    intent(in) :: G    !< ocean grid structure
  integer,                  intent(in) :: nz   !< The size of the arrays in the vertical
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, nz), &
                            intent(in) :: u_a  !< The first u-point array in arbitrary units [A]
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, nz), &
                            intent(in) :: u_b  !< The second u-point array in arbitrary units [B]
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output


end subroutine post_product_sum_u
module subroutine post_product_v(id, v_a, v_b, G, nz, diag, mask, alt_h)
  integer,                  intent(in) :: id   !< The ID for this diagnostic
  type(ocean_grid_type),    intent(in) :: G    !< ocean grid structure
  integer,                  intent(in) :: nz   !< The size of the arrays in the vertical
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB, nz), &
                            intent(in) :: v_a  !< The first v-point array in arbitrary units [A]
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB, nz), &
                            intent(in) :: v_b  !< The second v-point array in arbitrary units [B]
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output
  real,           optional, intent(in) :: mask(:,:,:)  !< If present, use this real array as the data mask [nondim]
  real,   target, optional, intent(in) :: alt_h(:,:,:) !< An alternate thickness to use for vertically
                                               !! remapping this diagnostic [H ~> m or kg m-2]

  ! Local variables

end subroutine post_product_v
module subroutine post_product_sum_v(id, v_a, v_b, G, nz, diag)
  integer,                  intent(in) :: id   !< The ID for this diagnostic
  type(ocean_grid_type),    intent(in) :: G    !< ocean grid structure
  integer,                  intent(in) :: nz   !< The size of the arrays in the vertical
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB, nz), &
                            intent(in) :: v_a  !< The first v-point array in arbitrary units [A]
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB, nz), &
                            intent(in) :: v_b  !< The second v-point array in arbitrary units [B]
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output


end subroutine post_product_sum_v
module subroutine post_xy_average(diag_cs, diag, field)
  type(diag_type),   intent(in) :: diag !< This diagnostic
  real,    target,   intent(in) :: field(:,:,:) !< Diagnostic field in arbitrary units [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_cs !< Diagnostics mediator control structure
  ! Local variable


end subroutine post_xy_average
module subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,            intent(in)    :: time_int_in !< The time interval [s] over which any
                                                !! values that are offered are valid.
  type(time_type), intent(in)    :: time_end_in !< The end time of the valid interval
  type(diag_ctrl), intent(inout) :: diag_CS     !< Structure used to regulate diagnostic output

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
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

end subroutine disable_averaging
module function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_CS  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag%time_int [s]
  type(time_type), optional, intent(out) :: time_end !< Current setting of diag%time_end
  logical :: query_averaging_enabled

end function query_averaging_enabled
module function get_diag_time_end(diag_cs)
  type(diag_ctrl), intent(in)  :: diag_CS !< Structure used to regulate diagnostic output
  type(time_type) :: get_diag_time_end
  !   This function returns the valid end time for diagnostics that are handled
  ! outside of the MOM6 infrastructure, such as via the generic tracer code.

end function get_diag_time_end
integer module function register_diag_field(module_name, field_name, axes_in, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  character(len=*),           intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                                        !! or "ice_shelf_model"
  character(len=*),           intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),     target, intent(in) :: axes_in   !< Container w/ up to 3 integer handles that
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
                                                         !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
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
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]
  logical,          optional, intent(in) :: v_extensive  !< True for vertically extensive fields (vertically
                                                         !! integrated). Default/absent for intensive.
  ! Local variables

end function register_diag_field
logical module function register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  integer,          intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
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
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]
  logical,          optional, intent(in) :: v_extensive !< True for vertically extensive fields (vertically
                                                         !! integrated). Default/absent for intensive.
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
integer module function xyz_method(axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
  type(axes_grp),             intent(in)  :: axes !< Container w/ up to 3 integer handles that indicates
                                                  !! axes for this field
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  logical,          optional, intent(in)  :: v_extensive !< True for vertically extensive fields
                                                         !! (vertically integrated). Default/absent for intensive.


  ! This is a simple way to encode the cell method information made from 3 strings
  ! (x_cell_method,y_cell_method,v_cell_method) in a 3 digit integer xyz
  ! x_cell_method,y_cell_method,v_cell_method can each be 'point' or 'sum' or 'mean'
  ! We can encode these with setting  1 for 'point', 2 for 'sum, 3 for 'mean' in
  ! the 100s position for x, 10s position for y, 1s position for z
  ! E.g., x:sum,y:point,z:mean is 213

end function xyz_method
module subroutine attach_cell_methods(id, axes, ostring, cell_methods, &
                               x_cell_method, y_cell_method, v_cell_method, v_extensive)
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
  character(len=*), optional, intent(in)  :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  logical,          optional, intent(in)  :: v_extensive !< True for vertically extensive fields
                                                         !! (vertically integrated). Default/absent for intensive.
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
module function register_static_field(module_name, field_name, axes, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            do_not_log, interp_method, tile_count, &
            cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
            x_cell_method, y_cell_method, area_cell_method, conversion)
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
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
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

end function register_static_field
module subroutine describe_option(opt_name, value, diag_CS)
  character(len=*), intent(in) :: opt_name !< The name of the option
  character(len=*), intent(in) :: value   !< A character string with the setting of the option.
  type(diag_ctrl),  intent(in) :: diag_CS !< Structure used to regulate diagnostic output


end subroutine describe_option
module function ocean_register_diag(var_desc, G, diag_CS, day)
  integer :: ocean_register_diag  !< An integer handle to this diagnostic.
  type(vardesc),         intent(in) :: var_desc !< The vardesc type describing the diagnostic
  type(ocean_grid_type), intent(in) :: G        !< The ocean's grid type
  type(diag_ctrl), intent(in), target :: diag_CS  !< The diagnostic control structure
  type(time_type),       intent(in) :: day      !< The current model time

                     ! as might be needed to convert from intensive to extensive
                     ! or for dimensional consistency testing [various] or [a A-1 ~> 1]

end function ocean_register_diag
module subroutine diag_mediator_infrastructure_init(err_msg)
  ! This subroutine initializes the FMS diag_manager.
  character(len=*), optional, intent(out)   :: err_msg !< An error message

end subroutine diag_mediator_infrastructure_init
module subroutine diag_mediator_init(G, GV, US, nz, param_file, diag_cs, doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The ocean grid type.
  type(verticalGrid_type), target, intent(in)  :: GV !< The ocean vertical grid structure
  type(unit_scale_type),   target, intent(in)  :: US !< A dimensional unit scaling type
  integer,                    intent(in)    :: nz    !< The number of layers in the model's native grid.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),            intent(inout) :: diag_cs !< A pointer to a type with many variables
                                                     !! used for diagnostics
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the
                                                     !! file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
  ! This include declares and sets the variable "version".

end subroutine diag_mediator_init
module subroutine diag_mediator_set_OBC_info(G, OBC_seg_u, OBC_seg_v, diag_cs)
  type(ocean_grid_type), intent(inout) :: G  !< The ocean grid type.
  integer, dimension(G%IsdB:G%IedB,G%jsd:G%jed),  &
                         intent(in) :: OBC_seg_u !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at u-points,
                                             !! with a value of 0 for no OBC, a positive value for an
                                             !! Eastern OBC or a negative value for a Western OBC
  integer, dimension(G%isd:G%ied,G%JsdB:G%JedB), &
                         intent(in) :: OBC_seg_v !< An array that indicates the presence and direction
                                             !! of any open boundary conditions at v-points,
                                             !! with a value of 0 for no OBC, a positive value for a
                                             !! Northern OBC or a negative value for a Southern OBC
  type(diag_ctrl),       intent(inout) :: diag_cs !< A defined type used to regulate diagnostics


end subroutine diag_mediator_set_OBC_info
module subroutine diag_set_state_ptrs(h, tv, diag_cs)
  real, dimension(:,:,:), target, intent(in   ) :: h !< the model thickness array [H ~> m or kg m-2]
  type(thermo_var_ptrs),  target, intent(in   ) :: tv !< A structure with thermodynamic variables that are
                                                      !! used to convert thicknesses to vertical extents
  type(diag_ctrl),                intent(inout) :: diag_cs !< diag mediator control structure

  ! Keep pointers to h, T, S needed for the diagnostic remapping
end subroutine diag_set_state_ptrs
module subroutine diag_update_remap_grids(diag_cs, alt_h, alt_T, alt_S, update_intensive, update_extensive )
  type(diag_ctrl),        intent(inout) :: diag_cs      !< Diagnostics control structure
  real, target, optional, intent(in   ) :: alt_h(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current thicknesses [H ~> m or kg m-2]
  real, target, optional, intent(in   ) :: alt_T(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current temperatures [C ~> degC]
  real, target, optional, intent(in   ) :: alt_S(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current salinity [S ~> ppt]
  logical, optional,      intent(in   ) :: update_intensive !< If true (default), update the grids used for
                                                            !! intensive diagnostics
  logical, optional,      intent(in   ) :: update_extensive !< If true (not default), update the grids used for
                                                            !! intensive diagnostics
  ! Local variables

end subroutine diag_update_remap_grids
module subroutine diag_masks_set(G, nz, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G  !< The ocean grid type.
  integer,                       intent(in) :: nz !< The number of layers in the model's native grid.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables

  ! 2d masks point to the model masks since they are identical
end subroutine diag_masks_set
module subroutine set_piecemeal_extents(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics

  ! Piecemeal buffers for 2d axes
end subroutine set_piecemeal_extents
module subroutine diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output


end subroutine diag_mediator_close_registration
module subroutine axes_grp_end(axes)
  type(axes_grp), intent(inout) :: axes   !< Axes group to be destroyed

end subroutine axes_grp_end
module subroutine diag_mediator_end(time, diag_CS, end_diag_manager)
  type(time_type),   intent(in)  :: time !< The current model time
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in)  :: end_diag_manager !< If true, call diag_manager_end()

  ! Local variables

end subroutine diag_mediator_end
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
  type(diag_ctrl),  intent(in) :: diag_CS  !< The diagnostics control structure
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
module subroutine diag_grid_storage_init(grid_storage, G, GV, diag)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  type(ocean_grid_type),   intent(in)    :: G           !< Horizontal grid
  type(verticalGrid_type), intent(in)    :: GV          !< ocean vertical grid structure
  type(diag_ctrl),         intent(in)    :: diag        !< Diagnostic control structure used as the constructor
                                                        !! template for this routine

end subroutine diag_grid_storage_init
module subroutine diag_copy_diag_to_storage(grid_storage, h_state, diag)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  real, dimension(:,:,:),  intent(in)    :: h_state     !< Current model thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),         intent(in)    :: diag     !< Diagnostic control structure used as the constructor


  ! Don't do anything else if there are no remapped coordinates
end subroutine diag_copy_diag_to_storage
module subroutine diag_copy_storage_to_diag(diag, grid_storage)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the constructor
  type(diag_grid_storage), intent(in)    :: grid_storage !< Structure containing a snapshot of the target grids


  ! Don't do anything else if there are no remapped coordinates
end subroutine diag_copy_storage_to_diag
module subroutine diag_save_grids(diag)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the constructor


  ! Don't do anything else if there are no remapped coordinates
end subroutine diag_save_grids
module subroutine diag_restore_grids(diag)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the constructor


  ! Don't do anything else if there are no remapped coordinates
end subroutine diag_restore_grids
module subroutine diag_grid_storage_end(grid_storage)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  ! Local variables

  ! Don't do anything else if there are no remapped coordinates
end subroutine diag_grid_storage_end
module subroutine downsample_diag_masks_set(G, nz, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G  !< The ocean grid type.
  integer,                       intent(in) :: nz !< The number of layers in the model's native grid.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables

!print*,'original c extents ',G%isc,G%iec,G%jsc,G%jec
!print*,'original c extents ',G%iscb,G%iecb,G%jscb,G%jecb
!print*,'coarse   c extents ',G%HId2%isc,G%HId2%iec,G%HId2%jsc,G%HId2%jec
!print*,'original d extents ',G%isd,G%ied,G%jsd,G%jed
!print*,'coarse   d extents ',G%HId2%isd,G%HId2%ied,G%HId2%jsd,G%HId2%jed
! original c  extents           5          52           5          52
! original cB-nonsym extents    5          52           5          52
! original cB-sym    extents    4          52           4          52
! coarse   c extents            3          26           3          26
! original d extents            1          56           1          56
! original dB-nonsym extents    1          56           1          56
! original dB-sym extents       0          56           0          56
! coarse   d extents            1          28           1          28

end subroutine downsample_diag_masks_set
module subroutine downsample_diag_indices_get(fo1, fo2, dl, diag_cs, isv, iev, jsv, jev)
  integer,           intent(in)  :: fo1     !< The size of the original diag field in x on data domain including halos
  integer,           intent(in)  :: fo2     !< The size of the original diag field in y on data domain including halos
  integer,           intent(in)  :: dl      !< Index of downsample level
  type(diag_ctrl),   intent(in)  :: diag_CS !< Structure used to regulate diagnostic output
  integer,           intent(out) :: isv     !< i-start index for diagnostics
  integer,           intent(out) :: iev     !< i-end index for diagnostics
  integer,           intent(out) :: jsv     !< j-start index for diagnostics
  integer,           intent(out) :: jev     !< j-end index for diagnostics
  ! Local variables

  !   The current implementation of the downsampled diagnostics assumes that the tracer-point
  ! computational domain on each processor can be evenly divided by dL in each direction, which
  ! avoids the need for halo updates or checks that the halo regions are up-to-date.  The following
  ! check that this assumption is true is only relevant if there are in fact downsampled diagnostics,
  ! which is why it occurs during the first call to this routine instead of during initialization.
end subroutine downsample_diag_indices_get
module subroutine downsample_diag_field_3d(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
  real, dimension(:,:,:), pointer :: locfield  !< Input array pointer in arbitrary units [A ~> a]
  real, dimension(:,:,:), allocatable, intent(inout) :: locfield_dsamp !< Output (downsampled) array [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: dl                !< Index of Level of down sampling
  integer, intent(inout) :: isv            !< i-start index for diagnostics
  integer, intent(inout) :: iev            !< i-end index for diagnostics
  integer, intent(inout) :: jsv            !< j-start index for diagnostics
  integer, intent(inout) :: jev            !< j-end index for diagnostics
  real,    optional,target, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask [nondim]
  ! Local variables

end subroutine downsample_diag_field_3d
module subroutine downsample_diag_field_2d(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
  real, dimension(:,:), pointer :: locfield !< Input array pointer in arbitrary units [A ~> a]
  real, dimension(:,:), allocatable, intent(inout) :: locfield_dsamp !< Output (downsampled) array [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: dl                !< Index of Level of down sampling
  integer, intent(inout) :: isv            !< i-start index for diagnostics
  integer, intent(inout) :: iev            !< i-end index for diagnostics
  integer, intent(inout) :: jsv            !< j-start index for diagnostics
  integer, intent(inout) :: jev            !< j-end index for diagnostics
  real,    optional,target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim].
  ! Local variables

end subroutine downsample_diag_field_2d
module subroutine downsample_field_3d(field_in, field_out, dL, method, mask, diag_cs, diag, &
                               isv_o, jsv_o, isv_d, iev_d, jsv_d, jev_d)
  real, dimension(:,:,:), pointer :: field_in      !< Original field to be downsampled in arbitrary units [A ~> a]
  real, dimension(:,:,:), allocatable :: field_out !< Downsampled field in the same arbitrary units [A ~> a]
  integer, intent(in) :: dL                !< Level of down sampling
  integer,  intent(in) :: method           !< Sampling method
  real,  dimension(:,:,:), pointer :: mask !< Mask for input field [nondim]
  type(diag_ctrl), intent(in) :: diag_CS   !< Structure used to regulate diagnostic output
  type(diag_type), intent(in) :: diag      !< A structure describing the diagnostic to post
  integer, intent(in) :: isv_o             !< Original i-start index
  integer, intent(in) :: jsv_o             !< Original j-start index
  integer, intent(in) :: isv_d             !< i-start index of down sampled data
  integer, intent(in) :: iev_d             !< i-end index of down sampled data
  integer, intent(in) :: jsv_d             !< j-start index of down sampled data
  integer, intent(in) :: jev_d             !< j-end index of down sampled data

  ! Local variables
                             ! indices when i or j is 0.
                    ! value [nondim], [L2 ~> m2], [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]
                    ! [A H L ~> a m2 or a kg m-1] or [A H L2 ~> a m3 or a kg]
                    ! value [nondim], [L2 ~> m2], [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]
                    ! [A H L ~> a m2 or a kg m-1] or [A H L2 ~> a m3 or a kg]
                    ! [A H L ~> a m2 or a kg m-1] or [A H L2 ~> a m3 or a kg]
                    ! value [nondim], [L2 ~> m2], [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]
                    ! [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]

end subroutine downsample_field_3d
module subroutine downsample_field_2d(field_in, field_out, dl, method, mask, diag_cs, diag, &
                               isv_o, jsv_o, isv_d, iev_d, jsv_d, jev_d)
  real, dimension(:,:), pointer :: field_in      !< Original field to be downsampled in arbitrary units [A ~> a]
  real, dimension(:,:), allocatable :: field_out !< Downsampled field in the same arbitrary units [A ~> a]
  integer, intent(in) :: dl                !< Level of down sampling
  integer,  intent(in) :: method           !< Sampling method
  real, dimension(:,:), pointer :: mask    !< Mask for input field [nondim]
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: isv_o             !< Original i-start index
  integer, intent(in) :: jsv_o             !< Original j-start index
  integer, intent(in) :: isv_d             !< i-start index of down sampled data
  integer, intent(in) :: iev_d             !< i-end index of down sampled data
  integer, intent(in) :: jsv_d             !< j-start index of down sampled data
  integer, intent(in) :: jev_d             !< j-end index of down sampled data

  ! Local variables
                             ! indices when i or j is 0.
                    ! value [nondim], [L2 ~> m2], [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]
                    ! [A H L ~> a m2 or a kg m-1] or [A H L2 ~> a m3 or a kg]
                    ! value [nondim], [L2 ~> m2], [H L ~> m2 or kg m-1] or [H L2 ~> m3 or kg]
                    ! [A H L ~> a m2 or a kg m-1] or [A H L2 ~> a m3 or a kg]

end subroutine downsample_field_2d
module function sum_1d(field, sz) result(sum)
  integer, intent(in) :: sz        !<  The size of the array to sum
  real,    intent(in) :: field(sz) !< The field to sum in arbitrary units [A ~> a]
  real :: sum !< The rotationally symmetric sum of the entries in field [A ~> a]

  ! Local variables

end function sum_1d
module function square_sum(field, sz, naive_sum) result(sum)
  integer,           intent(in) :: sz         !< The size of the array along each axis
  real,              intent(in) :: field(sz, sz) !< The field to sum in arbitrary units [A ~> a]
  logical, optional, intent(in) :: naive_sum !< If true, sum the elements in the order they appear in memory.
  real :: sum !< The sum of the entries in field [A ~> a]

  ! Local variables

end function square_sum
module subroutine downsample_mask_2d(mask_in, mask_dsamp, dL, method, isc_o, jsc_o, isd_o, jsd_o, &
                              isc_d, iec_d, jsc_d, jec_d, isd_d, ied_d, jsd_d, jed_d)
  integer, intent(in) :: isd_o !< Original data domain i-start index
  integer, intent(in) :: jsd_o !< Original data domain j-start index
  real, dimension(isd_o:,jsd_o:), intent(in) :: mask_in !< Original mask to be down sampled [nondim]
  real, dimension(:,:), pointer :: mask_dsamp   !< Down-sampled mask [nondim]
  integer, intent(in) :: method !< Sampling method
  integer, intent(in) :: dL    !< Level of down sampling
  integer, intent(in) :: isc_o !< Original i-start index
  integer, intent(in) :: jsc_o !< Original j-start index
  integer, intent(in) :: isc_d !< Computational i-start index of down sampled data
  integer, intent(in) :: iec_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsc_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jec_d !< Computational j-end index of down sampled data
  integer, intent(in) :: isd_d !< Data domain i-start index of down sampled data
  integer, intent(in) :: ied_d !< Data domain i-end index of down sampled data
  integer, intent(in) :: jsd_d !< Data domain j-start index of down sampled data
  integer, intent(in) :: jed_d !< Data domain j-end index of down sampled data

  ! Local variables
                             ! indices when i or j is 0.

  ! down sampled mask = 0 unless the mask value of one of the down sampling cells is 1
end subroutine downsample_mask_2d
module subroutine downsample_mask_3d(mask_in, mask_dsamp, dL, method, isc_o, jsc_o, isd_o, jsd_o, &
                              isc_d, iec_d, jsc_d, jec_d, isd_d, ied_d, jsd_d, jed_d)
  integer, intent(in) :: isd_o !< Original data domain i-start index
  integer, intent(in) :: jsd_o !< Original data domain j-start index
  real, dimension(isd_o:,jsd_o:,:), intent(in) :: mask_in !< Original mask to be down sampled [nondim]
  real, dimension(:,:,:), pointer :: mask_dsamp   !< Down-sampled mask [nondim]
  integer, intent(in) :: dL    !< Level of down sampling
  integer, intent(in) :: method !< Sampling method
  integer, intent(in) :: isc_o !< Original i-start index
  integer, intent(in) :: jsc_o !< Original j-start index
  integer, intent(in) :: isc_d !< Computational i-start index of down sampled data
  integer, intent(in) :: iec_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsc_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jec_d !< Computational j-end index of down sampled data
  integer, intent(in) :: isd_d !< Computational i-start index of down sampled data
  integer, intent(in) :: ied_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsd_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jed_d !< Computational j-end index of down sampled data

  ! Local variables
                             ! indices when i or j is 0.

  ! down sampled mask = 0 unless the mask value of one of the down sampling cells is 1
end subroutine downsample_mask_3d
logical module function found_in_diagtable(diag, varName)
  type(diag_ctrl),            intent(in) :: diag       !< A structure used to control diagnostics.
  character(len=*),           intent(in) :: varName    !< The obsolete diagnostic name
  ! Local

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
end function found_in_diagtable
module subroutine MOM_diag_send_complete()
end subroutine MOM_diag_send_complete
  end interface

end module MOM_diag_mediator
