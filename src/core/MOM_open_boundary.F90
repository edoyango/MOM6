! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Controls where open boundary conditions are applied
module MOM_open_boundary

use MOM_array_transform,      only : rotate_array, rotate_array_pair
use MOM_coms,                 only : sum_across_PEs, any_across_PEs
use MOM_coms,                 only : Set_PElist, Get_PElist, PE_here, num_PEs
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_debugging,            only : hchksum, uvchksum, chksum
use MOM_diag_mediator,        only : diag_ctrl, time_type
use MOM_domains,              only : pass_var, pass_vector
use MOM_domains,              only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,              only : To_All, EAST_FACE, NORTH_FACE, SCALAR_PAIR, CGRID_NE, CORNER
use MOM_dyn_horgrid,          only : dyn_horgrid_type
use MOM_error_handler,        only : MOM_mesg, MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser,          only : get_param, log_version, param_file_type, read_param
use MOM_grid,                 only : ocean_grid_type, hor_index_type
use MOM_interface_heights,    only : thickness_to_dz
use MOM_interpolate,          only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate,          only : external_field
use MOM_io,                   only : slasher, field_size, file_exists, stderr, SINGLE_FILE
use MOM_io,                   only : vardesc, query_vardesc, var_desc
use MOM_regridding,           only : regridding_CS
use MOM_remapping,            only : remappingSchemesDoc, remappingDefaultScheme, remapping_CS
use MOM_remapping,            only : initialize_remapping, remapping_core_h, end_remapping
use MOM_restart,              only : register_restart_field, register_restart_pair
use MOM_restart,              only : query_initialized, set_initialized, MOM_restart_CS
use MOM_string_functions,     only : extract_word, remove_spaces, uppercase, lowercase
use MOM_tidal_forcing,        only : astro_longitudes, astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency
use MOM_time_manager,         only : set_date, time_type, time_minus_signed
use MOM_tracer_registry,      only : tracer_type, tracer_registry_type, tracer_name_lookup
use MOM_unit_scaling,         only : unit_scale_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public open_boundary_apply_normal_flow
public open_boundary_config
public open_boundary_setup_vert
public open_boundary_halo_update
public open_boundary_query
public open_boundary_end
public open_boundary_impose_normal_slope
public open_boundary_impose_land_mask
public radiation_open_bdry_conds
public read_OBC_segment_data
public update_OBC_segment_data
public initialize_OBC_segment_reservoirs
public open_boundary_test_extern_uv
public open_boundary_test_extern_h
public open_boundary_zero_normal_flow
public parse_segment_str
public register_OBC, OBC_registry_init
public register_file_OBC, file_OBC_end
public segment_tracer_registry_init
public segment_tracer_registry_end
public segment_thickness_reservoir_init
public register_segment_tracer
public register_temp_salt_segments
public register_obgc_segments
public fill_temp_salt_segments
public fill_obgc_segments
public fill_thickness_segments
public set_obgc_segments_props
public setup_OBC_tracer_reservoirs
public setup_OBC_thickness_reservoirs
public open_boundary_register_restarts
public copy_thickness_reservoirs
public update_segment_tracer_reservoirs
public update_segment_thickness_reservoirs
public set_initialized_OBC_tracer_reservoirs
public update_OBC_ramp
public remap_OBC_fields
public rotate_OBC_config
public rotate_OBC_segment_direction
public write_OBC_info, chksum_OBC_segments
public initialize_segment_data
public flood_fill
public flood_fill2

integer, parameter, public :: OBC_NONE = 0      !< Indicates the use of no open boundary
integer, parameter, public :: OBC_DIRECTION_N = 100 !< Indicates the boundary is an effective northern boundary
integer, parameter, public :: OBC_DIRECTION_S = 200 !< Indicates the boundary is an effective southern boundary
integer, parameter, public :: OBC_DIRECTION_E = 300 !< Indicates the boundary is an effective eastern boundary
integer, parameter, public :: OBC_DIRECTION_W = 400 !< Indicates the boundary is an effective western boundary
!>@{ Enumeration values for OBC relative vorticity configurations
integer, parameter, public :: OBC_VORTICITY_NONE = 0
integer, parameter, public :: OBC_VORTICITY_ZERO = 1
integer, parameter, public :: OBC_VORTICITY_FREESLIP = 2
integer, parameter, public :: OBC_VORTICITY_COMPUTED = 3
integer, parameter, public :: OBC_VORTICITY_SPECIFIED = 4
!>@}
!>@{ Enumeration values for OBC strain configurations
integer, parameter, public :: OBC_STRAIN_NONE = 0
integer, parameter, public :: OBC_STRAIN_ZERO = 1
integer, parameter, public :: OBC_STRAIN_FREESLIP = 2
integer, parameter, public :: OBC_STRAIN_COMPUTED = 3
integer, parameter, public :: OBC_STRAIN_SPECIFIED = 4
!>@}
integer, parameter :: NUM_PHYS_FIELDS = 13  !< Number of physical fields
!>@{ Indices of physical field positions in segment%field array
integer, parameter :: &
    F_U = 1, F_V = 2, F_VX = 3, F_UY = 4, F_Z = 5, F_UAMP = 6, F_UPHASE = 7, &
    F_VAMP = 8, F_VPHASE = 9, F_ZAMP = 10, F_ZPHASE = 11, F_T = 12, F_S = 13
!>@}
character(len=8), parameter :: PHYS_FIELD_NAMES(NUM_PHYS_FIELDS) = &
    [character(len=8) :: 'U', 'V', 'DVDX', 'DUDY', 'SSH', 'Uamp', &
     'Uphase', 'Vamp', 'Vphase', 'SSHamp', 'SSHphase', 'TEMP', 'SALT']  !< Physical field name
                                                            !! strings used by input parameter

!> Open boundary segment data from files (mostly).
type, public :: OBC_segment_data_type
  type(external_field) :: handle            !< handle from FMS associated with segment data on disk
  type(external_field) :: dz_handle         !< handle from FMS associated with segment thicknesses on disk
  logical           :: required = .false.   !< True if this field is required
  logical           :: use_IO = .false.     !< True if segment data is based on file input
  character(len=32) :: name                 !< A name identifier for the segment data.  When there is grid
                                            !! rotation, this is the name on the rotated internal grid.
  integer           :: tr_index = -1        !< If this field is a tracer, its index in registry is stored here.
  logical           :: bgc_tracer           !< True if this field is a BGC tracer
  logical           :: on_face              !< If true, this field is discretized on the OBC segment
                                            !! (velocity-point) faces, or if false it as the vorticiy points
  real              :: scale                !< A scaling factor for converting input data to
                                            !! the internal units of this field.  For salinity this would
                                            !! be in units of [S ppt-1 ~> 1]
  real, allocatable :: buffer_src(:,:,:)    !< buffer for segment data located at cell faces and on
                                            !! the original vertical grid in the internally scaled
                                            !! units for the field in question, such as [L T-1 ~> m s-1]
                                            !! for a velocity or [S ~> ppt] for salinity.
  integer           :: nk_src               !< Number of vertical levels in the source data
  real, allocatable :: dz_src(:,:,:)        !< vertical grid cell spacing of the incoming segment
                                            !! data in [Z ~> m].
  real, allocatable :: buffer_dst(:,:,:)    !< buffer src data remapped to the target vertical grid
                                            !! in the internally scaled units for the field in
                                            !! question, such as [L T-1 ~> m s-1] for a velocity or
                                            !! [S ~> ppt] for salinity.
  real              :: value                !< A constant value for the inflow concentration if not read
                                            !! from file, in the internal units of a field, such as [S ~> ppt]
                                            !! for salinity.
  real              :: resrv_lfac_in = 1.   !< The reservoir inverse length scale factor for the inward
                                            !! direction per field [nondim].  The general 1/Lscale_in is
                                            !! multiplied by this factor for a specific tracer or thickness.
  real              :: resrv_lfac_out= 1.   !< The reservoir inverse length scale factor for the outward
                                            !! direction per field [nondim].  The general 1/Lscale_out is
                                            !! multiplied by this factor for a specific tracer or thickness.
end type OBC_segment_data_type

!> Tracer on OBC segment data structure, for putting into a segment tracer registry.
type, public :: OBC_segment_tracer_type
  real, allocatable          :: t(:,:,:)              !< tracer concentration array in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  real                       :: OBC_inflow_conc = 0.0 !< tracer concentration for generic inflows in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  character(len=32)          :: name                  !< tracer name used for error messages
  type(tracer_type), pointer :: Tr => NULL()          !< metadata describing the tracer
  real, allocatable          :: tres(:,:,:)           !< tracer reservoir array in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  real                       :: scale                 !< A scaling factor for converting the units of input
                                                      !! data, like [S ppt-1 ~> 1] for salinity.
  logical                    :: is_initialized        !< reservoir values have been set when True
  integer                    :: ntr_index = -1        !< index of segment tracer in the global tracer registry
  real              :: resrv_lfac_in = 1.   !< The reservoir inverse length scale factor for the inward
                                            !! direction per field [nondim].  The general 1/Lscale_in is
                                            !! multiplied by this factor for a specific tracer or thickness.
  real              :: resrv_lfac_out= 1.   !< The reservoir inverse length scale factor for the outward
                                            !! direction per field [nondim].  The general 1/Lscale_out is
                                            !! multiplied by this factor for a specific tracer or thickness.
end type OBC_segment_tracer_type

!> Thickness on OBC segment data structure, with a reservoir
type, public :: OBC_segment_thickness_type
  real, allocatable          :: h(:,:,:)              !< layer thickness array in rescaled units, [Z ~> m].
  real                       :: OBC_inflow_conc = 0.0 !< layer thickness for generic inflows in rescaled units,
                                                      !! [Z ~> m].
  character(len=32)          :: name                  !< thickness name used for error messages
  real, allocatable          :: h_res(:,:,:)          !< thickness reservoir array in rescaled units,
                                                      !! [Z ~> m].
  real                       :: scale                 !< A scaling factor for converting the units of input
                                                      !! data, [Z m-1 ~> 1].
  logical                    :: is_initialized        !< reservoir values have been set when True
  integer                    :: fd_index = -1         !< index of segment thickness in the input fields
end type OBC_segment_thickness_type

!> Registry type for tracers on segments
type, public :: segment_tracer_registry_type
  integer                       :: ntseg = 0         !< number of registered tracer segments
  type(OBC_segment_tracer_type) :: Tr(MAX_FIELDS_)   !< array of registered tracers
  logical                       :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                     !! When locked=.true.,no more tracers can be registered.
                                                     !! Not sure who should lock it or when...
end type segment_tracer_registry_type

!> Open boundary segment data structure.  Unless otherwise noted, 2-d and 3-d arrays are discretized
!! at the same position as normal velocity points in the middle of the OBC segments.
type, public :: OBC_segment_type
  logical :: Flather        !< If true, applies Flather + Chapman radiation of barotropic gravity waves.
  logical :: radiation      !< If true, 1D Orlanksi radiation boundary conditions are applied.
                            !! If False, a gradient condition is applied.
  logical :: radiation_tan  !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! tangential flows.
  logical :: radiation_grad !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! dudv and dvdx.
  logical :: oblique        !< Oblique waves supported at radiation boundary.
  logical :: oblique_tan    !< If true, 2D radiation boundary conditions are applied to
                            !! tangential flows.
  logical :: oblique_grad   !< If true, 2D radiation boundary conditions are applied to
                            !! dudv and dvdx.
  logical :: nudged         !< Optional supplement to radiation boundary.
  logical :: nudged_tan     !< Optional supplement to nudge tangential velocity.
  logical :: nudged_grad    !< Optional supplement to nudge normal gradient of tangential velocity.
  logical :: specified      !< Boundary normal velocity fixed to external value.
  logical :: specified_tan  !< Boundary tangential velocity fixed to external value.
  logical :: specified_grad !< Boundary gradient of tangential velocity fixed to external value.
  logical :: open           !< Boundary is open for continuity solver, and there are no other
                            !! parameterized mass fluxes at the open boundary.
  logical :: gradient       !< Zero gradient at boundary.
  integer :: direction      !< Boundary faces one of the four directions.
  logical :: is_N_or_S      !< True if the OB is facing North or South and exists on this PE.
  logical :: is_E_or_W      !< True if the OB is facing East or West and exists on this PE.
  logical :: is_E_or_W_2    !< True if the OB is facing East or West anywhere.
  type(OBC_segment_data_type), pointer :: field(:) => NULL()  !< OBC data
  integer :: num_fields     !< number of OBC data fields (e.g. u_normal,u_parallel and eta for Flather)
  integer :: Is_obc         !< Starting local i-index of boundary segment, this may be outside of the local PE.
  integer :: Ie_obc         !< Ending local i-index of boundary segment, this may be outside of the local PE.
  integer :: Js_obc         !< Starting local j-index of boundary segment, this may be outside of the local PE.
  integer :: Je_obc         !< Ending local j-index of boundary segment, this may be outside of the local PE.
  real :: Velocity_nudging_timescale_in  !< Nudging timescale on inflow [T ~> s].
  real :: Velocity_nudging_timescale_out !< Nudging timescale on outflow [T ~> s].
  logical :: on_pe          !< true if any portion of the segment is located in this PE's data domain
  logical :: temp_segment_data_exists !< true if temperature data arrays are present
  logical :: salt_segment_data_exists !< true if salinity data arrays are present
  real, allocatable :: Htot(:,:)  !< The total column thickness [H ~> m or kg m-2] at OBC-points.
  real, allocatable :: dZtot(:,:) !< The total column vertical extent [Z ~> m] at OBC segment faces.
  real, allocatable :: normal_vel(:,:,:)      !< The layer velocity normal to the OB
                                              !! segment [L T-1 ~> m s-1].
  real, allocatable :: tangential_vel(:,:,:)  !< The layer velocity tangential to the OB segment
                                              !! [L T-1 ~> m s-1], discretized at the corner points.
  real, allocatable :: tangential_grad(:,:,:) !< The gradient of the velocity tangential to the OB
                                              !! segment [T-1 ~> s-1], discretized at the corner points.
  real, allocatable :: normal_trans(:,:,:)    !< The layer transport normal to the OB
                                              !! segment [H L2 T-1 ~> m3 s-1].
  real, allocatable :: normal_vel_bt(:,:)     !< The barotropic velocity normal to
                                              !! the OB segment [L T-1 ~> m s-1].
  real, allocatable :: normal_trans_bt(:,:)   !< The barotropic transport normal
                                              !! the OB segment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, allocatable :: tidal_vn(:,:)          !< The barotropic tidal velocity normal to
                                              !! the OB segment [L T-1 ~> m s-1].
  real, allocatable :: tidal_vt(:,:)          !< The barotropic tidal velocity tangential to
                                              !! the OB segment [L T-1 ~> m s-1].
  real, allocatable :: SSH(:,:)               !< The sea-surface elevation along the
                                              !! segment [Z ~> m].
  real, allocatable :: tidal_elev(:,:)        !< Tidal elevation at the OBC points [Z ~> m]
  real, allocatable :: grad_normal(:,:,:)     !< The gradient of the normal flow along the
                                              !! segment times the grid spacing [L T-1 ~> m s-1],
                                              !! with the first index being the corner-point index
                                              !! along the segment, and the second index being 1 (for
                                              !! values one point into the domain) or 2 (for values
                                              !! along the OBC itself)
  real, allocatable :: grad_tan(:,:,:)        !< The gradient of the tangential flow along the
                                              !! segment times the grid spacing [L T-1 ~> m s-1], with the
                                              !! first index being the velocity/tracer point index along the
                                              !! segment, and the second being 1 for the value 1.5 points
                                              !! inside the domain and 2 for the value half a point
                                              !! inside the domain.
  real, allocatable :: grad_gradient(:,:,:)   !< The gradient normal to the segment of the gradient
                                              !! tangetial to the segment of tangential flow along the segment
                                              !! times the grid spacing [T-1 ~> s-1], with the first
                                              !! index being the velocity/tracer point index along the segment,
                                              !! and the second being 1 for the value 2 points into the domain
                                              !! and 2 for the value 1 point into the domain.
  real, allocatable :: rx_norm_rad(:,:,:)     !< The previous normal phase speed use for EW radiation
                                              !! OBC, in grid points per timestep [nondim]
  real, allocatable :: ry_norm_rad(:,:,:)     !< The previous normal phase speed use for NS radiation
                                              !! OBC, in grid points per timestep [nondim]
  real, allocatable :: rx_norm_obl(:,:,:)     !< The previous x-direction normalized radiation coefficient
                                              !! for either EW or NS oblique OBCs [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_norm_obl(:,:,:)     !< The previous y-direction normalized radiation coefficient
                                              !! for either EW or NS oblique OBCs [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal(:,:,:)      !< The denominator for oblique radiation of the normal
                                              !! velocity [L2 T-2 ~> m2 s-2]
  real, allocatable :: nudged_normal_vel(:,:,:) !< The layer velocity normal to the OB segment
                                              !! that values should be nudged towards [L T-1 ~> m s-1].
  real, allocatable :: nudged_tangential_vel(:,:,:) !< The layer velocity tangential to the OB segment
                                              !! that values should be nudged towards [L T-1 ~> m s-1],
                                              !! discretized at the corner (PV) points.
  real, allocatable :: nudged_tangential_grad(:,:,:)  !< The layer dvdx or dudy towards which nudging
                                              !! can occur [T-1 ~> s-1].
  type(OBC_segment_thickness_type), pointer  :: h_Reg=> NULL()!< A pointer to the thickness for the segment.
  type(segment_tracer_registry_type), pointer  :: tr_Reg=> NULL()!< A pointer to the tracer registry for the segment.
  type(hor_index_type) :: HI !< Horizontal index ranges
  real :: Tr_InvLscale_out                                  !< An effective inverse length scale for restoring
                                                            !! the tracer concentration in a fictitious
                                                            !! reservoir towards interior values when flow
                                                            !! is exiting the domain [L-1 ~> m-1]
  real :: Tr_InvLscale_in                                   !< An effective inverse length scale for restoring
                                                            !! the tracer concentration towards an externally
                                                            !! imposed value when flow is entering [L-1 ~> m-1]
  real :: Th_InvLscale_out                                  !< An effective inverse length scale for restoring
                                                            !! the layer thickness in a fictitious
                                                            !! reservoir towards interior values when flow
                                                            !! is exiting the domain [L-1 ~> m-1]
  real :: Th_InvLscale_in                                   !< An effective inverse length scale for restoring
                                                            !! the layer thickness towards an externally
                                                            !! imposed value when flow is entering [L-1 ~> m-1]
end type OBC_segment_type

!> Open-boundary data
type, public :: ocean_OBC_type
  integer :: number_of_segments = 0                   !< The number of open-boundary segments.
  logical :: reverse_segment_order = .false.          !< If true, store the segments internally in the reversed order.
  integer :: ke = 0                                   !< The number of model layers
  logical :: open_u_BCs_exist_globally = .false.      !< True if any zonal velocity points
                                                      !! in the global domain use open BCs.
  logical :: open_v_BCs_exist_globally = .false.      !< True if any meridional velocity points
                                                      !! in the global domain use open BCs.
  logical :: Flather_u_BCs_exist_globally = .false.   !< True if any zonal velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: Flather_v_BCs_exist_globally = .false.   !< True if any meridional velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: oblique_BCs_exist_globally = .false.     !< True if any velocity points
                                                      !! in the global domain use oblique BCs.
  logical :: nudged_u_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: nudged_v_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: specified_u_BCs_exist_globally = .false. !< True if any zonal velocity points
                                                      !! in the global domain use specified BCs.
  logical :: specified_v_BCs_exist_globally = .false. !< True if any meridional velocity points
                                                      !! in the global domain use specified BCs.
  logical :: radiation_BCs_exist_globally = .false.   !< True if radiations BCs are in use anywhere.
  logical :: user_BCs_set_globally = .false.          !< True if any OBC_USER_CONFIG is set
                                                      !! for input from user directory.
  logical :: update_OBC = .false.                     !< Is OBC data time-dependent
  logical :: update_OBC_seg_data = .false.            !< Is it the time for OBC segment data update for fields that
                                                      !! require less frequent update
  logical :: any_needs_IO_for_data = .false.          !< Is any i/o needed for OBCs globally
  integer :: vorticity_config                         !< An integer indicating OBC relative vorticity configuration
  integer :: strain_config                            !< An integer indicating OBC strain configuration
  logical :: zero_biharmonic = .false.                !< If True, zeros the Laplacian of flow on open boundaries for
                                                      !! use in the biharmonic viscosity term.
  logical :: brushcutter_mode = .false.               !< If True, read data on supergrid.
  logical, allocatable :: tracer_x_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with x reservoirs (needed for restarts).
  logical, allocatable :: tracer_y_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with y reservoirs (needed for restarts).
  logical :: thickness_x_reservoirs_used = .false.    !< True for thichness reservoirs in x (needed for restarts).
  logical :: thickness_y_reservoirs_used = .false.    !< True for thichness reservoirs in y (needed for restarts).
  integer                       :: ntr = 0            !< number of tracers
  integer :: n_tide_constituents = 0                  !< Number of tidal constituents to add to the boundary.
  logical :: add_tide_constituents = .false.          !< If true, add tidal constituents to the boundary elevation
                                                      !! and velocity. Will be set to true if n_tide_constituents > 0.
  character(len=2), allocatable, dimension(:) :: tide_names  !< Names of tidal constituents to add to the boundary data.
  real, allocatable, dimension(:) :: tide_frequencies !< Angular frequencies of chosen tidal
                                                      !! constituents [rad T-1 ~> rad s-1].
  real, allocatable, dimension(:) :: tide_eq_phases   !< Equilibrium phases of chosen tidal constituents [rad].
  real, allocatable, dimension(:) :: tide_fn          !< Amplitude modulation of boundary tides by nodal cycle [nondim].
  real, allocatable, dimension(:) :: tide_un          !< Phase modulation of boundary tides by nodal cycle [rad].
  logical :: add_eq_phase = .false.                   !< If true, add the equilibrium phase argument
                                                      !! to the specified boundary tidal phase.
  logical :: add_nodal_terms = .false.                !< If true, insert terms for the 18.6 year modulation when
                                                      !! calculating tidal boundary conditions.
  type(time_type) :: time_ref                         !< Reference date (t = 0) for tidal forcing.
  type(astro_longitudes) :: tidal_longitudes          !< Lunar and solar longitudes used to calculate tidal forcing.
  ! Properties of the segments used.
  type(OBC_segment_type), allocatable :: segment(:)   !< List of segment objects.
  ! Which segment object describes the current point.
  integer, allocatable :: segnum_u(:,:) !< The absolute value gives the segment number of any OBCs at u-points,
                                        !! while the sign indicates whether they are Eastern (> 0) or Western (< 0)
                                        !! OBCs, with 0 for velocities that are not on an OBC.
  integer, allocatable :: segnum_v(:,:) !< The absolute value gives the segment number of any OBCs at v-points,
                                        !! while the sign indicates whether they are Northern (> 0) or Southern (< 0)
                                        !! OBCs, with 0 for velocities that are not on an OBC.
  ! Keep the OBC segment properties for external BGC tracers
  type(external_tracers_segments_props), pointer :: obgc_segments_props => NULL() !< obgc segment properties
  integer :: num_obgc_tracers = 0       !< The total number of obgc tracers

  ! The following parameters are used in the baroclinic radiation code:
  real :: gamma_uv !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for velocities [nondim].
                   !! Valid values range from 0 to 1, with a default of 0.3.
  real :: rx_max   !< The maximum magnitude of the baroclinic radiation velocity (or speed of
                   !! characteristics) in units of grid points per timestep [nondim].
  logical :: OBC_pe !< Is there an open boundary on this tile?
  logical :: u_OBCs_on_PE   !< True if there are any u-point OBCs on this PE, including in its halos.
  logical :: v_OBCs_on_PE   !< True if there are any v-point OBCs on this PE, including in its halos.
  logical :: v_N_OBCs_on_PE !< True if there are any northern v-point OBCs on this PE, including in its halos.
  logical :: v_S_OBCs_on_PE !< True if there are any southern v-point OBCs on this PE, including in its halos.
  logical :: u_E_OBCs_on_PE !< True if there are any eastern u-point OBCs on this PE, including in its halos.
  logical :: u_W_OBCs_on_PE !< True if there are any western u-point OBCs on this PE, including in its halos.
  !>@{ Index ranges on the local PE for the open boundary conditions in various directions
  integer :: Is_u_W_obc, Ie_u_W_obc, js_u_W_obc, je_u_W_obc
  integer :: Is_u_E_obc, Ie_u_E_obc, js_u_E_obc, je_u_E_obc
  integer :: is_v_S_obc, ie_v_S_obc, Js_v_S_obc, Je_v_S_obc
  integer :: is_v_N_obc, ie_v_N_obc, Js_v_N_obc, Je_v_N_obc
  !>@}
  type(remapping_CS), pointer :: remap_z_CS => NULL() !< ALE remapping control structure for
                                                      !! z-space data on segments
  type(remapping_CS), pointer :: remap_h_CS => NULL() !< ALE remapping control structure for
                                                      !! thickness-based fields on segments
  type(OBC_registry_type), pointer :: OBC_Reg => NULL()  !< Registry type for boundaries
  real, allocatable :: rx_normal(:,:,:)     !< Array storage for normal phase speed for EW radiation OBCs
                                            !! in units of grid points per timestep [nondim]
  real, allocatable :: ry_normal(:,:,:)     !< Array storage for normal phase speed for NS radiation OBCs
                                            !! in units of grid points per timestep [nondim]
  real, allocatable :: rx_oblique_u(:,:,:)  !< X-direction oblique boundary condition radiation speeds
                                            !! squared at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_oblique_u(:,:,:)  !< Y-direction oblique boundary condition radiation speeds
                                            !! squared at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: rx_oblique_v(:,:,:)  !< X-direction oblique boundary condition radiation speeds
                                            !! squared at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_oblique_v(:,:,:)  !< Y-direction oblique boundary condition radiation speeds
                                            !! squared at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal_u(:,:,:)  !< Denominator for normalizing EW oblique boundary condition
                                            !! radiation rates at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal_v(:,:,:)  !< Denominator for normalizing NS oblique boundary condition
                                            !! radiation rates at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: tres_x(:,:,:,:)      !< Array storage of tracer reservoirs for restarts,
                                            !! in unscaled units [conc]
  real, allocatable :: tres_y(:,:,:,:)      !< Array storage of tracer reservoirs for restarts,
                                            !! in unscaled units [conc]
  real, allocatable :: h_res_x(:,:,:)       !< Array storage of thickness reservoirs for restarts,
                                            !! [Z ~> m]
  real, allocatable :: h_res_y(:,:,:)       !< Array storage of thickness reservoirs for restarts,
                                            !! [Z ~> m]
  logical :: use_h_res = .false.            !< If true, use thickness reservoirs
  logical :: debug                          !< If true, write verbose checksums for debugging purposes.
  integer :: nk_OBC_debug = 0               !< The number of layers of OBC segment data to write out
                                            !! in full when DEBUG_OBCS is true.
  real :: silly_h  !< A silly value of thickness outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [Z ~> m].
  real :: silly_u  !< A silly value of velocity outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [L T-1 ~> m s-1].
  logical :: ramp = .false.                 !< If True, ramp from zero to the external values for SSH.
  logical :: ramping_is_activated = .false. !< True if the ramping has been initialized
  real :: ramp_timescale                    !< If ramp is True, use this timescale for ramping [T ~> s].
  real :: trunc_ramp_time                   !< If ramp is True, time after which ramp is done [T ~> s].
  real :: ramp_value                        !< If ramp is True, where we are on the ramp from
                                            !! zero to one [nondim].
  type(time_type) :: ramp_start_time        !< Time when model was started.
  integer :: remap_answer_date  !< The vintage of the order of arithmetic and expressions to use
                                !! for remapping.  Values below 20190101 recover the remapping
                                !! answers from 2018, while higher values use more robust
                                !! forms of the same remapping expressions.
  logical :: check_reconstruction !< Flag for remapping to run checks on reconstruction
  logical :: check_remapping      !< Flag for remapping to run internal checks
  logical :: force_bounds_in_subcell !< Flag for remapping to hide overshoot using bounds
  logical :: om4_remap_via_sub_cells !< If true, use the OM4 remapping algorithm
  character(40) :: remappingScheme !< String selecting the vertical remapping scheme
  type(group_pass_type) :: pass_oblique  !< Structure for group halo pass
  logical :: exterior_OBC_bug   !< If true, use incorrect form of tracers exterior to OBCs.
  logical :: hor_index_bug      !< If true, recover set of a horizontal indexing bugs in the OBC code.
  logical :: reservoir_init_bug !< If true, set the OBC tracer reservoirs at the startup of a new
                                !! run from the interior tracer concentrations regardless of
                                !! properties that may be explicitly specified for the reservoir
                                !! concentrations.
  logical :: ts_needed_bug      !< If true, recover a bug that temperature and salinity can be ignored
                                !! even if they are registered tracers in the rest of the model.
end type ocean_OBC_type

!> Control structure for open boundaries that read from files.
!! Probably lots to update here.
type, public :: file_OBC_CS ; private
  logical :: OBC_file_used = .false.     !< Placeholder for now to avoid an empty type.
end type file_OBC_CS

!> Type to carry something (what??) for the OBC registry.
type, public :: OBC_struct_type
  character(len=32)               :: name             !< OBC name used for error messages
end type OBC_struct_type

!> Type to carry basic OBC information needed for updating values.
type, public :: OBC_registry_type
  integer               :: nobc = 0          !< number of registered open boundary types.
  type(OBC_struct_type) :: OB(MAX_FIELDS_)   !< array of registered boundary types.
  logical               :: locked = .false.  !< New OBC types may be registered if locked=.false.
                                             !! When locked=.true.,no more boundaries can be registered.
end type OBC_registry_type

!> Type to carry OBC information needed for setting segments for OBGC tracers
type, private :: external_tracers_segments_props
   type(external_tracers_segments_props), pointer :: next => NULL() !< pointer to the next node
   character(len=128) :: tracer_name      !< tracer name
   character(len=128) :: tracer_src_file  !< tracer source file for BC
   character(len=128) :: tracer_src_field !< name of the field in source file to extract BC
   real               :: lfac_in  !< multiplicative factor for inbound  tracer reservoir length scale [nondim]
   real               :: lfac_out !< multiplicative factor for outbound tracer reservoir length scale [nondim]
end type external_tracers_segments_props
integer :: id_clock_pass !< A CPU time clock

character(len=40)  :: mdl = "MOM_open_boundary" !< This module's name.


  interface
module subroutine open_boundary_config(G, US, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure

  ! Local variables
                             ! recreate the bugs, or if false bugs are only used if actively selected.
                             ! of the open boundary condition code.
  ! This include declares and sets the variable "version".

end subroutine open_boundary_config
module subroutine open_boundary_setup_vert(GV, US, OBC)
  type(verticalGrid_type), intent(in)    :: GV  !< Container for vertical grid information
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure

  ! Local variables

end subroutine open_boundary_setup_vert
module subroutine segment_determine_required_fields(segment, tides, temp_salt)
  type(OBC_segment_type), intent(inout) :: segment !< OBC segment
  logical, optional, intent(in) :: tides         !< Switch for tidal variables
  logical, optional, intent(in) :: temp_salt     !< Switch for thermodynamic variables

  ! Local variables

end subroutine segment_determine_required_fields
integer module function find_phys_field_index(name)
  character(len=*), intent(in) :: name !< Field name

  ! Local variables

end function find_phys_field_index
module subroutine OBC_any_IO(OBC)
  type(ocean_OBC_type), intent(inout) :: OBC !< Open boundary control structure

  ! Local variables

end subroutine OBC_any_IO
module subroutine allocate_segment_field_data(field, OBC, segment, US, inputdir, filename, varname, &
                                       suffix, value, turns, nz)
  type(OBC_segment_data_type), &
                          intent(inout) :: field    !< A field of the segment
  type(ocean_OBC_type),   intent(in)    :: OBC      !< Open boundary control structure
  type(OBC_segment_type), intent(inout) :: segment  !< Segment to work on
  type(unit_scale_type),  intent(in)    :: US       !< A dimensional unit scaling type
  character(len=*),       intent(in)    :: inputdir !< The directory of input files
  character(len=*),       intent(in)    :: filename !< Input file name
  character(len=*),       intent(in)    :: varname  !< Variable name in the input file
  character(len=*),       intent(in)    :: suffix   !< Variable name suffix, "_segment_xxx"
  real,                   intent(in)    :: value    !< Unscaled specified value of the field [a]
  integer,                intent(in)    :: turns    !< Number of quarter turns of the grid
  integer,                intent(in)    :: nz       !< Default k-axis size in buffer_dst

  ! Local variables

end subroutine allocate_segment_field_data
module subroutine initialize_segment_data(GV, US, OBC, PF, turns, use_temperature)
  type(verticalGrid_type),      intent(in)    :: GV  !< Container for vertical grid information
  type(unit_scale_type),        intent(in)    :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type), target, intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),        intent(in)    :: PF  !< Parameter file handle
  integer,                      intent(in)    :: turns !< Number of quarter turns of the grid
  logical,                      intent(in)    :: use_temperature !< If true, temperature and
                                                 !! salinity used as state variables.

  ! Local variables

end subroutine initialize_segment_data
logical module function field_is_on_face(name, is_E_or_W)
  character(len=*), intent(in) :: name       !< The OBC segment data name to interpret
  logical,          intent(in) :: is_E_or_W  !< This is true for an eastern or western open boundary condition

end function field_is_on_face
logical module function field_is_tidal(name)
  character(len=*), intent(in) :: name       !< The OBC segment data name to interpret

end function field_is_tidal
module subroutine set_segnum_signs(OBC, G)
  type(ocean_OBC_type),   intent(inout) :: OBC !< Open boundary control structure, perhaps on a rotated grid.
  type(dyn_horgrid_type), intent(in)    :: G   !< Ocean grid structure used by OBC


end subroutine set_segnum_signs
real module function scale_factor_from_name(name, US, Tr_Reg)
  character(len=*),        intent(in) :: name  !< The OBC segment data name to interpret
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(segment_tracer_registry_type), pointer :: Tr_Reg  !< pointer to tracer registry for this segment


end function scale_factor_from_name
module subroutine initialize_obc_tides(OBC, US, param_file)
  type(ocean_OBC_type), intent(inout) :: OBC  !< Open boundary control structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< Parameter file handle

end subroutine initialize_obc_tides
module subroutine setup_segment_indices(G, seg, Is_obc, Ie_obc, Js_obc, Je_obc)
  type(dyn_horgrid_type), intent(in) :: G !< grid type
  type(OBC_segment_type), intent(inout) :: seg  !< Open boundary segment
  integer, intent(in) :: Is_obc !< Q-point global i-index of start of segment
  integer, intent(in) :: Ie_obc !< Q-point global i-index of end of segment
  integer, intent(in) :: Js_obc !< Q-point global j-index of start of segment
  integer, intent(in) :: Je_obc !< Q-point global j-index of end of segment
  ! Local variables

  ! Isg, Ieg will be I*_obc in global space
end subroutine setup_segment_indices
module subroutine setup_u_point_obc(OBC, G, US, segment_str, l_seg, l_seg_io, PF, reentrant_y)
  type(ocean_OBC_type),    intent(inout) :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "I=%,J=%:%,string"
  integer,                 intent(in) :: l_seg !< The internal segment number
  integer,                 intent(in) :: l_seg_io !< The segment number used for reading parameters
  type(param_file_type), intent(in)   :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_y !< is the domain reentrant in y?
  ! Local variables
  ! This returns the global indices for the segment
end subroutine setup_u_point_obc
module subroutine setup_v_point_obc(OBC, G, US, segment_str, l_seg, l_seg_io, PF, reentrant_x)
  type(ocean_OBC_type),    intent(inout) :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "J=%,I=%:%,string"
  integer,                 intent(in) :: l_seg !< The internal segment number
  integer,                 intent(in) :: l_seg_io !< The segment number used for reading parameters
  type(param_file_type),   intent(in) :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_x !< is the domain reentrant in x?
  ! Local variables

  ! This returns the global indices for the segment
end subroutine setup_v_point_obc
module subroutine parse_segment_str(ni_global, nj_global, segment_str, l, m, n, action_str, reentrant)
  integer,          intent(in)  :: ni_global !< Number of h-points in zonal direction
  integer,          intent(in)  :: nj_global !< Number of h-points in meridional direction
  character(len=*), intent(in)  :: segment_str !< A string in form of "I=l,J=m:n,string" or "J=l,I=m,n,string"
  integer,          intent(out) :: l !< The value of I=l, if segment_str begins with I=l, or the value of J=l
  integer,          intent(out) :: m !< The value of J=m, if segment_str begins with I=, or the value of I=m
  integer,          intent(out) :: n !< The value of J=n, if segment_str begins with I=, or the value of I=n
  character(len=*), intent(out) :: action_str(:) !< The "string" part of segment_str
  logical,          intent(in)  :: reentrant !< is domain reentrant in relevant direction?
  ! Local variables
                                                    !! "I=%,J=%:%,string"

  ! Process first word which will started with either 'I=' or 'J='
end subroutine parse_segment_str
module subroutine parse_segment_manifest_str(segment_str, num_fields, fields)
  character(len=*), intent(in) :: segment_str   !< A string in form of
                                        !< "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
  integer, intent(out) :: num_fields    !< The number of fields in the segment data
  character(len=*), dimension(NUM_PHYS_FIELDS), intent(out) :: fields
                                        !< List of fieldnames for each segment

  ! Local variables

end subroutine parse_segment_manifest_str
module subroutine parse_segment_data_str(segment_str, idx, var, value, filename, fieldname)
  character(len=*), intent(in) :: segment_str   !< A string in form of
      !! "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
  integer, intent(in) :: idx                    !< Index of segment_str record
  character(len=*), intent(in) :: var           !< The name of the variable for which parameters are needed
  character(len=*), intent(out) :: filename     !< The name of the input file if using "file" method
  character(len=*), intent(out) :: fieldname    !< The name of the variable in the input file if using
                                                !! "file" method
  real, optional, intent(out)  :: value         !< A constant value if using the "value" method in various
                                                !! units but without the internal rescaling [various units]

  ! Local variables

  ! Process first word which will start with the fieldname
end subroutine parse_segment_data_str
module subroutine parse_for_tracer_reservoirs(OBC, PF, use_temperature)
  type(ocean_OBC_type), target, intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),  intent(in) :: PF  !< Parameter file handle
  logical,                intent(in) :: use_temperature !< If true, T and S are used

  ! Local variables

end subroutine parse_for_tracer_reservoirs
module subroutine open_boundary_halo_update(G, OBC)
  type(ocean_grid_type),   intent(in) :: G   !< Ocean grid structure
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure

  ! Local variables

end subroutine open_boundary_halo_update
logical module function open_boundary_query(OBC, apply_open_OBC, apply_specified_OBC, apply_Flather_OBC, &
                                     apply_nudged_OBC, needs_ext_seg_data)
  type(ocean_OBC_type), pointer    :: OBC !< Open boundary control structure
  logical, optional,    intent(in) :: apply_open_OBC      !< Returns True if open_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_specified_OBC !< Returns True if specified_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_Flather_OBC   !< Returns True if Flather_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_nudged_OBC    !< Returns True if nudged_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: needs_ext_seg_data  !< Returns True if external segment data needed
end function open_boundary_query
module subroutine open_boundary_dealloc(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure

end subroutine open_boundary_dealloc
module subroutine open_boundary_end(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
end subroutine open_boundary_end
module subroutine open_boundary_impose_normal_slope(OBC, G, depth)
  type(ocean_OBC_type),             pointer       :: OBC   !< Open boundary control structure
  type(dyn_horgrid_type),           intent(in)    :: G     !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: depth !< Bathymetry at h-points, in [Z ~> m] or other units
  ! Local variables

end subroutine open_boundary_impose_normal_slope
module subroutine open_boundary_impose_land_mask(OBC, G, areaCu, areaCv, US)
  type(ocean_OBC_type),              pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),             intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: areaCu !< Area of a u-cell [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: areaCv !< Area of a u-cell [L2 ~> m2]
  ! Local variables

end subroutine open_boundary_impose_land_mask
module subroutine setup_OBC_tracer_reservoirs(G, GV, OBC, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),        intent(in)    :: GV  !< The ocean's vertical grid structure
  type(ocean_OBC_type), target,   intent(inout) :: OBC !< Open boundary control structure
  type(MOM_restart_CS), optional, intent(in)    :: restart_CS !< MOM restart control structure

  ! Local variables
                          ! For salinity the units would be [ppt S-1 ~> 1]

end subroutine setup_OBC_tracer_reservoirs
module subroutine setup_OBC_thickness_reservoirs(G, GV, OBC, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),        intent(in)    :: GV  !< The ocean's vertical grid structure
  type(ocean_OBC_type), target,   intent(inout) :: OBC !< Open boundary control structure
  type(MOM_restart_CS), optional, intent(in)    :: restart_CS !< MOM restart control structure

  ! Local variables
                          ! [m Z-1 ~> 1]

end subroutine setup_OBC_thickness_reservoirs
module subroutine set_initialized_OBC_tracer_reservoirs(G, OBC, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G   !< Ocean grid structure
  type(ocean_OBC_type),           intent(in)    :: OBC !< Open boundary control structure
  type(MOM_restart_CS),           intent(inout) :: restart_CS !< MOM restart control structure

end subroutine set_initialized_OBC_tracer_reservoirs
module subroutine copy_thickness_reservoirs(OBC, G, GV)
  type(ocean_grid_type),          intent(inout) :: G     !< Ocean grid structure
  type(verticalGrid_type),        intent(in)    :: GV    !< The ocean's vertical grid structure
  type(ocean_OBC_type),           pointer       :: OBC   !< Open boundary control structure
  ! Local variables

end subroutine copy_thickness_reservoirs
module subroutine radiation_open_bdry_conds(OBC, u_new, u_old, v_new, v_old, G, GV, US, dt)
  type(ocean_grid_type),                      intent(inout) :: G     !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV    !< The ocean's vertical grid structure
  type(ocean_OBC_type),                       pointer       :: OBC   !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u_new !< On exit, new u values on open boundaries
                                                                     !! On entry, the old time-level u but including
                                                                     !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u_old !< Original unadjusted u [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v_new !< On exit, new v values on open boundaries.
                                                                     !! On entry, the old time-level v but including
                                                                     !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v_old !< Original unadjusted v  [L T-1 ~> m s-1]
  type(unit_scale_type),                      intent(in)    :: US    !< A dimensional unit scaling type
  real,                                       intent(in)    :: dt    !< Appropriate timestep [T ~> s]
  ! Local variables
                   ! discretized at the corner (PV) points.
end subroutine radiation_open_bdry_conds
module subroutine open_boundary_apply_normal_flow(OBC, G, GV, u, v)
  ! Arguments
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                     intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< u field to update on open
                                                                  !! boundaries [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< v field to update on open
                                                                  !! boundaries [L T-1 ~> m s-1]
  ! Local variables

end subroutine open_boundary_apply_normal_flow
module subroutine open_boundary_zero_normal_flow(OBC, G, GV, u, v)
  ! Arguments
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                      intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV  !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< u field to update on open boundaries [arbitrary]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< v field to update on open boundaries [arbitrary]
  ! Local variables

end subroutine open_boundary_zero_normal_flow
module subroutine gradient_at_q_points(G, GV, segment, uvel, vvel)
  type(ocean_grid_type),   intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  type(OBC_segment_type), intent(inout) :: segment !< OBC segment structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uvel !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vvel !< meridional velocity [L T-1 ~> m s-1]

end subroutine gradient_at_q_points
module function lookup_seg_field(OBC_seg, field)
  type(OBC_segment_type), intent(in) :: OBC_seg !< OBC segment
  character(len=32), intent(in) :: field !< The field name
  integer :: lookup_seg_field
  ! Local variables

end function lookup_seg_field
module function get_tracer_index(OBC_seg,tr_name)
  type(OBC_segment_type), pointer :: OBC_seg !< OBC segment
  character(len=*), intent(in) :: tr_name   !< The field name
  integer :: get_tracer_index
end function get_tracer_index
module subroutine allocate_OBC_segment_data(OBC, segment)
  type(ocean_OBC_type),   intent(in)    :: OBC     !< Open boundary structure
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment
  ! Local variables

end subroutine allocate_OBC_segment_data
module subroutine deallocate_OBC_segment_data(segment)
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment

end subroutine deallocate_OBC_segment_data
module subroutine open_boundary_test_extern_uv(G, GV, OBC, u, v)
  type(ocean_grid_type),                     intent(in)    :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< The ocean's vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),intent(inout) :: u !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),intent(inout) :: v !< Meridional velocity [L T-1 ~> m s-1]
  ! Local variables

end subroutine open_boundary_test_extern_uv
module subroutine open_boundary_test_extern_h(G, GV, OBC, h)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !<  Ocean vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)),intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2]
  ! Local variables

end subroutine open_boundary_test_extern_h
module subroutine read_OBC_segment_data(G, GV, US, OBC, tv, h, Time)
  type(ocean_grid_type),                     intent(in) :: G    !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in) :: US   !< A dimensional unit scaling type
  type(ocean_OBC_type),                      pointer    :: OBC  !< Open boundary structure
  type(thermo_var_ptrs),                     intent(in) :: tv   !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Thickness [H ~> m or kg m-2]
  type(time_type),                           intent(in) :: Time !< Model time

  ! Local variables

end subroutine read_OBC_segment_data
module subroutine update_OBC_segment_data(G, GV, US, OBC, h, Time)
  type(ocean_grid_type),                     intent(in) :: G    !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in) :: US   !< A dimensional unit scaling type
  type(ocean_OBC_type),                      pointer    :: OBC  !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Thickness [H ~> m or kg m-2]
  type(time_type),                           intent(in) :: Time !< Model time

  ! Local variables

end subroutine update_OBC_segment_data
module subroutine initialize_OBC_segment_reservoirs(GV, OBC)
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary structure

  ! Local variables

end subroutine initialize_OBC_segment_reservoirs
module subroutine update_OBC_ramp(Time, OBC, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(ocean_OBC_type),    intent(inout) :: OBC      !< Open boundary structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables

end subroutine update_OBC_ramp
module subroutine register_OBC(name, param_file, Reg)
  character(len=32),     intent(in)  :: name        !< OBC name used for error messages
  type(param_file_type), intent(in)  :: param_file  !< file to parse for  model parameter values
  type(OBC_registry_type), pointer   :: Reg         !< pointer to the tracer registry

end subroutine register_OBC
module subroutine OBC_registry_init(param_file, Reg)
  type(param_file_type),   intent(in) :: param_file !< open file to parse for model parameters
  type(OBC_registry_type), pointer    :: Reg        !< pointer to OBC registry



end subroutine OBC_registry_init
module function register_file_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(file_OBC_CS),        pointer    :: CS         !< file control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_file_OBC

end function register_file_OBC
module subroutine file_OBC_end(CS)
  type(file_OBC_CS), pointer    :: CS   !< OBC file control structure.

end subroutine file_OBC_end
module subroutine segment_tracer_registry_init(param_file, segment)
  type(param_file_type),      intent(in)      :: param_file !< open file to parse for model parameters
  type(OBC_segment_type), intent(inout)       :: segment    !<  the segment


! This include declares and sets the variable "version".
  !character(len=256) :: mesg    ! Message for error messages.

end subroutine segment_tracer_registry_init
module subroutine segment_thickness_reservoir_init(GV, US, OBC, param_file)
  type(param_file_type),  intent(in)    :: param_file !< open file to parse for model parameters
  type(verticalGrid_type), intent(in)   :: GV         !< ocean vertical grid structure
  type(unit_scale_type),  intent(in)    :: US         !< Unit scaling type
  type(ocean_OBC_type),   pointer       :: OBC        !< Open boundary structure
! real,         optional, intent(in)    :: OBC_scalar !< If present, use scalar value for segment tracer
!                                                     !! inflow concentration, including any rescaling to
!                                                     !! put the tracer concentration into its internal units,
!                                                     !! like [S ~> ppt] for salinity.
! logical,      optional, intent(in)    :: OBC_array  !< If true, use array values for segment tracer
!                                                     !! inflow concentration.
! Local variables
                  ! salinity, or other various units depending on what rescaling has occurred previously.

! This include declares and sets the variable "version".

end subroutine segment_thickness_reservoir_init
module subroutine register_segment_tracer(tr_ptr, ntr_index, param_file, GV, segment, &
                                   OBC_scalar, OBC_array, scale, resrv_lfac_in, resrv_lfac_out)
  type(verticalGrid_type), intent(in)   :: GV         !< ocean vertical grid structure
  type(tracer_type), target             :: tr_ptr     !< A target that can be used to set a pointer to the
                                                      !! stored value of tr. This target must be
                                                      !! an enduring part of the control structure,
                                                      !! because the tracer registry will use this memory,
                                                      !! but it also means that any updates to this
                                                      !! structure in the calling module will be
                                                      !! available subsequently to the tracer registry.
  integer, intent(in)                   :: ntr_index  !< index of segment tracer in the global tracer registry
  type(param_file_type),  intent(in)    :: param_file !< file to parse for model parameter values
  type(OBC_segment_type), intent(inout) :: segment    !< current segment data structure
  real,         optional, intent(in)    :: OBC_scalar !< If present, use scalar value for segment tracer
                                                      !! inflow concentration, including any rescaling to
                                                      !! put the tracer concentration into its internal units,
                                                      !! like [S ~> ppt] for salinity.
  logical,      optional, intent(in)    :: OBC_array  !< If true, use array values for segment tracer
                                                      !! inflow concentration.
  real,         optional, intent(in)    :: scale      !< A scaling factor that should be used with any
                                                      !! data that is read in to convert it to the internal
                                                      !! units of this tracer, in units like [S ppt-1 ~> 1]
                                                      !! for salinity.
  real,         optional, intent(in)    :: resrv_lfac_in   !< The reservoir inverse length scale factor
  real,         optional, intent(in)    :: resrv_lfac_out  !< The reservoir inverse length scale factor

! Local variables
                  ! salinity, or other various units depending on what rescaling has occurred previously.

end subroutine register_segment_tracer
module subroutine segment_tracer_registry_end(Reg)
  type(segment_tracer_registry_type), pointer :: Reg        !< pointer to tracer registry

! Local variables

end subroutine segment_tracer_registry_end
module subroutine segment_thickness_registry_end(Reg)
  type(OBC_segment_thickness_type), pointer :: Reg        !< pointer to thickness reservoir

! Local variables

end subroutine segment_thickness_registry_end
module subroutine register_temp_salt_segments(GV, US, OBC, tr_Reg, param_file)
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US         !< Unit scaling type
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer       :: tr_Reg     !< Tracer registry
  type(param_file_type),      intent(in)    :: param_file !< file to parse for  model parameter values

  ! Local variables

end subroutine register_temp_salt_segments
module subroutine set_obgc_segments_props(OBC,tr_name,obc_src_file_name,obc_src_field_name,lfac_in,lfac_out)
  type(ocean_OBC_type),pointer  :: OBC                !< Open boundary structure
  character(len=*),  intent(in) :: tr_name            !< Tracer name
  character(len=*),  intent(in) :: obc_src_file_name  !< OBC source file name
  character(len=*),  intent(in) :: obc_src_field_name !< name of the field in the source file
  real,              intent(in) :: lfac_in            !< factors for tracer reservoir inbound length scales [nondim]
  real,              intent(in) :: lfac_out           !< factors for tracer reservoir outbound length scales [nondim]

                                                                    ! the tracer segment properties
end subroutine set_obgc_segments_props
module subroutine get_obgc_segments_props(node, tr_name,obc_src_file_name,obc_src_field_name,lfac_in,lfac_out)
  type(external_tracers_segments_props),pointer :: node !< pointer to tracer segment properties
  character(len=*), intent(out) :: tr_name            !< Tracer name
  character(len=*), intent(out) :: obc_src_file_name  !< OBC source file name
  character(len=*), intent(out) :: obc_src_field_name !< name of the field in the source file
  real,             intent(out) :: lfac_in   !< multiplicative factor for inbound  reservoir length scale [nondim]
  real,             intent(out) :: lfac_out  !< multiplicative factor for outbound reservoir length scale [nondim]
end subroutine get_obgc_segments_props
module subroutine register_obgc_segments(GV, OBC, tr_Reg, param_file, tr_name)
  type(verticalGrid_type),    intent(in) :: GV         !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer    :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer    :: tr_Reg     !< Tracer registry
  type(param_file_type),      intent(in) :: param_file !< file to parse for model parameter values
  character(len=*),           intent(in) :: tr_name    !< Tracer name

  ! Local variables
                                                    ! relaxation length scale [nondim].
                                                    ! relaxation length scale [nondim].

end subroutine register_obgc_segments
module subroutine fill_obgc_segments(G, GV, OBC, tr_ptr, tr_name)
  type(ocean_grid_type),      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV     !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer       :: OBC    !< Open boundary structure
  real, dimension(:,:,:),     pointer       :: tr_ptr !< Pointer to tracer field in scaled concentration
                                                      !! units, like [S ~> ppt] for salinity.
  character(len=*),           intent(in)    :: tr_name !< Tracer name
! Local variables

end subroutine fill_obgc_segments
module subroutine fill_temp_salt_segments(G, GV, US, OBC, tv)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< Unit scaling
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure


end subroutine fill_temp_salt_segments
module subroutine fill_thickness_segments(G, GV, US, OBC, h)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< Unit scaling
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]


end subroutine fill_thickness_segments
module subroutine mask_outside_OBCs(G, US, param_file, OBC)
  type(dyn_horgrid_type),       intent(inout) :: G          !< Ocean grid structure
  type(param_file_type),        intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),         pointer       :: OBC        !< Open boundary structure
  type(unit_scale_type),        intent(in)    :: US         !< A dimensional unit scaling type

  ! Local variables
                                                      ! two different ways [nondim]

end subroutine mask_outside_OBCs
module subroutine flood_fill(G, color, cin, cout, cland)
  type(dyn_horgrid_type), intent(inout) :: G      !< Ocean grid structure
  real, dimension(:,:),   intent(inout) :: color  !< For sorting inside from outside [nondim]
  integer, intent(in) :: cin    !< color for inside the domain
  integer, intent(in) :: cout   !< color for outside the domain
  integer, intent(in) :: cland  !< color for inside the land mask

! Local variables

end subroutine flood_fill
module subroutine flood_fill2(G, color, cin, cout, cland)
  type(dyn_horgrid_type), intent(inout) :: G       !< Ocean grid structure
  real, dimension(:,:),   intent(inout) :: color   !< For sorting inside from outside [nondim]
  integer, intent(in) :: cin    !< color for inside the domain
  integer, intent(in) :: cout   !< color for outside the domain
  integer, intent(in) :: cland  !< color for inside the land mask

! Local variables

end subroutine flood_fill2
module subroutine open_boundary_register_restarts(HI, GV, US, OBC, Reg, param_file, restart_CS, &
                                           use_temperature)
  type(hor_index_type),    intent(in) :: HI !< Horizontal indices
  type(verticalGrid_type), pointer    :: GV !< Container for vertical grid information
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer    :: OBC !< OBC data structure, data intent(inout)
  type(tracer_registry_type), pointer :: Reg !< pointer to tracer registry
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure
  logical,                 intent(in) :: use_temperature !< If true, T and S are used
  ! Local variables

end subroutine open_boundary_register_restarts
module subroutine update_segment_tracer_reservoirs(G, GV, uhr, vhr, h, OBC, Reg)
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in) :: GV  !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhr !< accumulated volume/mass flux through
                                                                !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vhr !< accumulated volume/mass flux through
                                                                !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< layer thickness after advection
                                                                !! [H ~> m or kg m-2]
  type(ocean_OBC_type),                       pointer    :: OBC !< Open boundary structure
  type(tracer_registry_type),                 pointer    :: Reg !< pointer to tracer registry

  ! Local variables
                            ! flow is from the interior to the reservoir [H L2 ~> m3 or kg].
                            ! direction per field [nondim]
                            ! direction per field [nondim]
                            ! length scale [nondim]
                            ! 1 if the length scale of reservoir is zero [nondim]
                            ! e.g. a_in is -1 only if b_in ==1 and uhr or vhr is inward
                            ! e.g. a_out is 1 only if b_out==1 and uhr or vhr is outward
                            ! It's clear that a_in and a_out cannot be both non-zero [nondim]
                            ! For salinity the units would be [ppt S-1 ~> 1]

end subroutine update_segment_tracer_reservoirs
module subroutine update_segment_thickness_reservoirs(G, GV, uhr, vhr, h, OBC)
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in) :: GV  !<  Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhr !< accumulated volume/mass flux through
                                                                !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vhr !< accumulated volume/mass flux through
                                                                !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< layer thickness after advection
                                                                !! [H ~> m or kg m-2]
  type(ocean_OBC_type),                       pointer    :: OBC !< Open boundary structure

  ! Local variable
                          ! length scale [nondim]
                          ! length scale [nondim]
                          ! For salinity the units would be [ppt S-1 ~> 1]
                          ! direction per field [nondim]
                          ! direction per field [nondim]
                          ! 1 if the length scale of reservoir is zero [nondim]
                          ! e.g. a_in is -1 only if b_in ==1 and uhr or vhr is inward
                          ! e.g. a_out is 1 only if b_out==1 and uhr or vhr is outward
                          ! It's clear that a_in and a_out cannot be both non-zero [nondim]
end subroutine update_segment_thickness_reservoirs
module subroutine remap_OBC_fields(G, GV, h_old, h_new, OBC, PCM_cell)
  type(ocean_grid_type),                     intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in) :: GV    !<  Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old !< Thickness of source grid [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_new !< Thickness of destination grid [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer    :: OBC   !< Open boundary structure
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in) :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables

                        ! For salinity the units would be [S ~> ppt].
                        ! For salinity the units would be [ppt S-1 ~> 1].

end subroutine remap_OBC_fields
module subroutine adjustSegmentEtaToFitBathymetry(G, GV, US, segment, fld, at_node)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(OBC_segment_type),  intent(inout) :: segment !< OBC segment
  integer,                 intent(in)    :: fld  !< field index to adjust thickness
  logical,                 intent(in)    :: at_node !< True this point is at the OBC nodes rather than the faces

  ! real :: dilate      ! A factor by which to dilate the water column [nondim]
  !character(len=100) :: mesg

end subroutine adjustSegmentEtaToFitBathymetry
module subroutine rotate_OBC_config(OBC_in, G_in, OBC, G, turns)
  type(ocean_OBC_type), pointer, intent(in)    :: OBC_in !< Input OBC
  type(dyn_horgrid_type),        intent(in)    :: G_in  !< Input grid
  type(ocean_OBC_type), pointer, intent(inout) :: OBC   !< Rotated OBC
  type(dyn_horgrid_type),        intent(in)    :: G     !< Rotated grid
  integer,                       intent(in)    :: turns !< Number of quarter turns


end subroutine rotate_OBC_config
module subroutine rotate_OBC_segment_config(segment_in, G_in, segment, G, turns)
  type(OBC_segment_type), intent(in) :: segment_in  !< Input OBC segment
  type(dyn_horgrid_type),  intent(in) :: G_in       !< Input grid metric
  type(OBC_segment_type), intent(inout) :: segment  !< Rotated OBC segment
  type(dyn_horgrid_type),  intent(in) :: G          !< Rotated grid metric
  integer, intent(in) :: turns                      !< Number of quarter turns

  ! Global segment indices

  ! NOTE: A "rotation" of the OBC segment string would allow us to use
  !   setup_[uv]_point_obc to set up most of this.  For now, we just copy/swap
  !   flags and manually rotate the indices.

  ! This is set if the segment is in the local grid
end subroutine rotate_OBC_segment_config
module function rotate_OBC_segment_direction(direction, turns) result(rotated_dir)
  integer, intent(in) :: direction  !< The orientation of an OBC segment on the original grid
  integer, intent(in) :: turns      !< Number of quarter turns
  integer :: rotated_dir  !< An integer encoding the new rotated segment direction


end function rotate_OBC_segment_direction
module function rotated_field_name(input_name, turns)
  character(len=*),    intent(in) :: input_name !< The unrotated field name
  integer,             intent(in) :: turns !< Number of quarter turns of the grid
  character(len=len(input_name)) :: rotated_field_name !< The rotated field name

end function rotated_field_name
module subroutine allocate_rotated_seg_data(src_array, HI_in, tgt_array, segment)
  real, dimension(:,:,:), intent(in) :: src_array !< The segment data on the unrotated source grid
  type(hor_index_type),   intent(in) :: HI_in !< Horizontal indices on the source grid
  real, dimension(:,:,:), allocatable, intent(inout) :: tgt_array !< The segment data that is being allocated
  type(OBC_segment_type), intent(inout) :: segment !< OBC segment on the target grid

  ! Local variables

end subroutine allocate_rotated_seg_data
module subroutine write_OBC_info(OBC, G, GV, US)
  type(ocean_OBC_type),    pointer    :: OBC   !< An open boundary condition control structure
  type(ocean_grid_type),   intent(in) :: G     !< Rotated grid metric
  type(verticalGrid_type), intent(in) :: GV    !< Vertical grid
  type(unit_scale_type),   intent(in) :: US    !< Unit scaling

  ! Local variables
                        ! without grid rotation

end subroutine write_OBC_info
module subroutine chksum_OBC_segments(OBC, G, GV, US, nk)
  type(ocean_OBC_type),    intent(in) :: OBC   !< An open boundary condition control structure
  type(ocean_grid_type),   intent(in) :: G     !< Rotated grid metric
  type(verticalGrid_type), intent(in) :: GV    !< Vertical grid
  type(unit_scale_type),   intent(in) :: US    !< Unit scaling
  integer,                 intent(in) :: nk    !< The number of layers to print

  ! Local variables

end subroutine chksum_OBC_segments
module subroutine chksum_OBC_segment_data(segment, GV, US, nk, nseg_out)
  type(OBC_segment_type),  intent(in) :: segment !< Segment type to checksum
  type(verticalGrid_type), intent(in) :: GV    !< Vertical grid
  type(unit_scale_type),   intent(in) :: US    !< Unit scaling
  integer,                 intent(in) :: nk    !< The number of layers to print
  integer,                 intent(in) :: nseg_out !< The segment number reported in output

  ! Local variables

end subroutine chksum_OBC_segment_data
  end interface

end module MOM_open_boundary
