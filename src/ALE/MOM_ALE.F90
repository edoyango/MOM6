! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the main regridding routines.
!!
!! Regridding comprises two steps:
!! 1. Interpolation and creation of a new grid based on target interface
!!    densities (or any other criterion).
!! 2. Remapping of quantities between old grid and new grid.
!!
!! Original module written by Laurent White, 2008.06.09
module MOM_ALE

use MOM_debugging,        only : check_column_integrals
use MOM_diag_mediator,    only : register_diag_field, post_data, diag_ctrl
use MOM_diag_mediator,    only : time_type, diag_update_remap_grids, query_averaging_enabled
use MOM_domains,          only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,    only : MOM_error, FATAL, WARNING
use MOM_error_handler,    only : callTree_showQuery
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_hybgen_unmix,     only : hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix, hybgen_unmix_CS
use MOM_hybgen_regrid,    only : hybgen_regrid_CS
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_interface_heights,only : find_eta, calc_derived_thermo
use MOM_open_boundary,    only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,    only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_regridding,       only : initialize_regridding, regridding_main, end_regridding
use MOM_regridding,       only : uniformResolution
use MOM_regridding,       only : inflate_vanished_layers_old
use MOM_regridding,       only : regridding_preadjust_reqs, convective_adjustment
use MOM_regridding,       only : set_target_densities_from_GV, set_target_densities
use MOM_regridding,       only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding,       only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding,       only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding,       only : regriddingDefaultMinThickness
use MOM_regridding,       only : regridding_CS, set_regrid_params, write_regrid_file
use MOM_regridding,       only : getCoordinateInterfaces
use MOM_regridding,       only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding,       only : getStaticThickness
use MOM_remapping,        only : initialize_remapping, end_remapping
use MOM_remapping,        only : remapping_core_h, remapping_core_w
use MOM_remapping,        only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,        only : interpolate_column, reintegrate_column
use MOM_remapping,        only : remapping_CS, dzFromH1H2, remapping_set_param
use MOM_string_functions, only : uppercase, extractWord, extract_integer
use MOM_tracer_registry,  only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,     only : get_thickness_units, verticalGrid_type

!use regrid_consts,       only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts,        only : coordinateUnits, coordinateMode, state_dependent
use regrid_edge_values,   only : edge_values_implicit_h4
use PLM_functions,        only : PLM_reconstruction, PLM_boundary_extrapolation
use PLM_functions,        only : PLM_extrapolate_slope, PLM_monotonized_slope, PLM_slope_wa
use PPM_functions,        only : PPM_reconstruction, PPM_boundary_extrapolation
use Recon1d_PLM_WLS,      only : PLM_WLS

implicit none ; private
#include <MOM_memory.h>


!> ALE control structure
type, public :: ALE_CS ; private
  logical :: remap_uv_using_old_alg !< If true, uses the old "remapping via a delta z"
                                    !! method. If False, uses the new method that
                                    !! remaps between grids described by h.
  logical :: partial_cell_vel_remap !< If true, use partial cell thicknesses at velocity points
                                    !! that are masked out where they extend below the shallower
                                    !! of the neighboring bathymetry for remapping velocity.

  real :: regrid_time_scale !< The time-scale used in blending between the current (old) grid
                            !! and the target (new) grid [T ~> s]

  type(regridding_CS) :: regridCS !< Regridding parameters and work arrays
  type(remapping_CS)  :: remapCS  !< Remapping parameters and work arrays
  type(remapping_CS)  :: vel_remapCS  !< Remapping parameters for velocities and work arrays

  type(hybgen_unmix_CS), pointer :: hybgen_unmixCS => NULL() !< Parameters for hybgen remapping

  logical :: use_hybgen_unmix   !< If true, use the hybgen unmixing code before regridding
  logical :: do_conv_adj        !< If true, do convective adjustment before regridding

  integer :: nk             !< Used only for queries, not directly by this module
  real :: BBL_h_vel_mask    !< The thickness of a bottom boundary layer within which velocities in
                            !! thin layers are zeroed out after remapping, following practice with
                            !! Hybgen remapping, or a negative value to avoid such filtering
                            !! altogether, in [H ~> m or kg m-2].
  real :: h_vel_mask        !< A thickness at velocity points below which near-bottom layers are
                            !! zeroed out after remapping, following the practice with Hybgen
                            !! remapping, or a negative value to avoid such filtering altogether,
                            !! in [H ~> m or kg m-2].

  logical :: remap_after_initialization !< Indicates whether to regrid/remap after initializing the state.

  integer :: answer_date    !< The vintage of the expressions and order of arithmetic to use for
                            !! remapping. Values below 20190101 result in the use of older, less
                            !! accurate expressions that were in use at the end of 2018.  Higher
                            !! values result in the use of more robust and accurate forms of
                            !! mathematically equivalent expressions.

  logical :: conserve_ke    !< Apply a correction to the baroclinic velocity after remapping to
                            !! conserve KE.

  logical :: debug   !< If true, write verbose checksums for debugging purposes.
  logical :: show_call_tree !< For debugging

  ! for diagnostics
  type(diag_ctrl), pointer           :: diag                          !< structure to regulate output
  integer, dimension(:), allocatable :: id_tracer_remap_tendency      !< diagnostic id
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency     !< diagnostic id
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency_2d  !< diagnostic id
  logical, dimension(:), allocatable :: do_tendency_diag              !< flag for doing diagnostics
  integer                            :: id_dzRegrid = -1              !< diagnostic id

  ! diagnostic for fields prior to applying ALE remapping
  integer :: id_u_preale = -1 !< diagnostic id for zonal velocity before ALE.
  integer :: id_v_preale = -1 !< diagnostic id for meridional velocity before ALE.
  integer :: id_h_preale = -1 !< diagnostic id for layer thicknesses before ALE.
  integer :: id_T_preale = -1 !< diagnostic id for temperatures before ALE.
  integer :: id_S_preale = -1 !< diagnostic id for salinities before ALE.
  integer :: id_e_preale = -1 !< diagnostic id for interface heights before ALE.
  integer :: id_vert_remap_h = -1      !< diagnostic id for layer thicknesses used for remapping
  integer :: id_vert_remap_h_tendency = -1 !< diagnostic id for layer thickness tendency due to ALE
  integer :: id_remap_delta_integ_u2 = -1  !< Change in depth-integrated rho0*u**2/2
  integer :: id_remap_delta_integ_v2 = -1  !< Change in depth-integrated rho0*v**2/2

end type

! Publicly available functions
public ALE_init
public ALE_end
public ALE_regrid
public ALE_offline_inputs
public ALE_regrid_accelerated
public ALE_remap_scalar
public ALE_remap_tracers
public ALE_remap_velocities
public ALE_remap_set_h_vel, ALE_remap_set_h_vel_via_dz
public ALE_remap_interface_vals
public ALE_remap_vertex_vals
public ALE_PLM_edge_values
public TS_PLM_edge_values
public TS_PPM_edge_values
public TS_PLM_WLS_edge_values
public adjustGridForIntegrity
public ALE_initRegridding
public ALE_getCoordinate
public ALE_getCoordinateUnits
public ALE_writeCoordinateFile
public ALE_updateVerticalGridType
public ALE_initThicknessToCoord
public ALE_update_regrid_weights
public pre_ALE_diagnostics
public pre_ALE_adjustments
public ALE_remap_init_conds
public ALE_register_diags
public ALE_set_extrap_boundaries

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine ALE_init( param_file, G, GV, US, max_depth, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(ocean_grid_type),   intent(in) :: G          !< Grid structure
  type(verticalGrid_type), intent(in) :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  real,                    intent(in) :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(ALE_CS),            pointer    :: CS         !< Module control structure

  ! Local variables
                                                         ! for sharing parameters.

end subroutine ALE_init
module subroutine ALE_set_extrap_boundaries( param_file, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(ALE_CS),            pointer    :: CS         !< Module control structure

end subroutine ALE_set_extrap_boundaries
module subroutine ALE_set_OM4_remap_algorithm( CS, om4_remap_via_sub_cells )
  type(ALE_CS), pointer :: CS !< Module control structure
  logical, intent(in)   :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm

end subroutine ALE_set_OM4_remap_algorithm
module subroutine ALE_register_diags(Time, G, GV, US, diag, CS)
  type(time_type),target,     intent(in)  :: Time  !< Time structure
  type(ocean_grid_type),      intent(in)  :: G     !< Grid structure
  type(unit_scale_type),      intent(in)  :: US    !< A dimensional unit scaling type
  type(verticalGrid_type),    intent(in)  :: GV    !< Ocean vertical grid structure
  type(diag_ctrl), target,    intent(in)  :: diag  !< Diagnostics control structure
  type(ALE_CS), pointer                   :: CS    !< Module control structure

  ! Local variables

end subroutine ALE_register_diags
module subroutine adjustGridForIntegrity( CS, G, GV, h )
  type(ALE_CS),                              intent(in)    :: CS  !< Regridding parameters and options
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                   intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h   !< Current 3D grid thickness that
                                                                  !! are to be adjusted [H ~> m or kg m-2]
end subroutine adjustGridForIntegrity
module subroutine ALE_end(CS)
  type(ALE_CS), pointer :: CS  !< module control structure

  ! Deallocate memory used for the regridding
end subroutine ALE_end
module subroutine pre_ALE_diagnostics(G, GV, US, h, u, v, tv, CS)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options

  ! Local variables

end subroutine pre_ALE_diagnostics
module subroutine pre_ALE_adjustments(G, GV, US, h, tv, Reg, CS, u, v)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]


  ! Do column-wise convective adjustment.
  ! Tracers and velocities should probably also undergo consistent adjustments.
end subroutine pre_ALE_adjustments
module subroutine ALE_regrid( G, GV, US, h, h_new, dzRegrid, tv, CS, frac_shelf_h, PCM_cell)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h   !< Layer thicknesses in 3D grid before
                                                                   !! regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out)   :: h_new !< Layer thicknesses in 3D grid after
                                                                   !! regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out)  :: dzRegrid !< The change in grid interface positions
                                                                   !! due to regridding, in the same units as
                                                                   !! thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: frac_shelf_h !< Fractional ice shelf coverage [nondim]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out)   :: PCM_cell !< If true, use PCM remapping in a cell.

  ! Local variables

end subroutine ALE_regrid
module subroutine ALE_offline_inputs(CS, G, GV, US, h, tv, Reg, uhtr, vhtr, Kd, debug, OBC)
  type(ALE_CS),                                 pointer       :: CS    !< Regridding parameters and options
  type(ocean_grid_type),                        intent(in   ) :: G     !< Ocean grid informations
  type(verticalGrid_type),                      intent(in   ) :: GV    !< Ocean vertical grid structure
  type(unit_scale_type),                        intent(in   ) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(inout) :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(inout) :: tv    !< Thermodynamic variable structure
  type(tracer_registry_type),                   pointer       :: Reg   !< Tracer registry structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: uhtr  !< Zonal mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(inout) :: vhtr  !< Meridional mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(inout) :: Kd    !< Input diffusivities
                                                                       !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  logical,                                      intent(in   ) :: debug !< If true, then turn checksums
  type(ocean_OBC_type),                         pointer       :: OBC   !< Open boundary structure
  ! Local variables

end subroutine ALE_offline_inputs
module subroutine ALE_regrid_accelerated(CS, G, GV, US, h, tv, n_itt, u, v, OBC, Reg, dt, &
                                  dzRegrid, initial)
  type(ALE_CS),            pointer       :: CS     !< ALE control structure
  type(ocean_grid_type),   intent(inout) :: G      !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV     !< Vertical grid
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h      !< Original thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv     !< Thermo vars (T/S/EOS)
  integer,                 intent(in)    :: n_itt  !< Number of times to regrid
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary structure
  type(tracer_registry_type), &
                 optional, pointer       :: Reg    !< Tracer registry to remap onto new grid
  real,          optional, intent(in)    :: dt     !< Model timestep to provide a timescale for regridding [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                 optional, intent(inout) :: dzRegrid !< Final change in interface positions [H ~> m or kg m-2]
  logical,       optional, intent(in)    :: initial !< Whether we're being called from an initialization
                                                    !! routine (and expect diagnostics to work)

  ! Local variables
                                                               ! velocity points [H ~> m or kg m-2]
                                                               ! velocity points [H ~> m or kg m-2]
                                                               ! velocity points [H ~> m or kg m-2]
                                                               ! velocity points [H ~> m or kg m-2]

  ! we have to keep track of the total dzInterface if for some reason
  ! we're using the old remapping algorithm for u/v
                                                             ! an iteration [H ~> m or kg m-2]

end subroutine ALE_regrid_accelerated
module subroutine ALE_remap_tracers(CS, G, GV, h_old, h_new, Reg, debug, dt, PCM_cell)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old        !< Thickness of source grid
                                                                           !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness of destination grid
                                                                           !! [H ~> m or kg m-2]
  type(tracer_registry_type),                pointer       :: Reg          !< Tracer registry structure
  logical,                         optional, intent(in)    :: debug  !< If true, show the call tree
  real,                            optional, intent(in)    :: dt     !< time step for diagnostics [T ~> s]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in)    :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables
                                                       ! content [Conc H T-1 ~> Conc m s-1 or Conc kg m-2 s-1] or
                                                       ! cell thickness [H T-1 ~> m s-1 or kg m-2 s-1]
                                                       ! content [Conc H T-1 ~> Conc m s-1 or Conc kg m-2 s-1]

end subroutine ALE_remap_tracers
module subroutine ALE_remap_set_h_vel(CS, G, GV, h_new, h_u, h_v, OBC, debug)
  type(ALE_CS),                              intent(in)    :: CS      !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G       !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV      !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new   !< Thickness at tracer points of the
                                                                      !! grid being interpolated to velocity
                                                                      !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u     !< Grid thickness at zonal velocity
                                                                      !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v     !< Grid thickness at meridional velocity
                                                                      !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC     !< Open boundary structure
  logical,                         optional, intent(in)    :: debug   !< If true, show the call tree

  ! Local variables

end subroutine ALE_remap_set_h_vel
module subroutine ALE_remap_set_h_vel_via_dz(CS, G, GV, h_new, h_u, h_v, OBC, h_old, dzInterface, debug)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness at tracer points of the
                                                                           !! grid being interpolated to velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC          !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_old        !< Thickness of source grid when generating
                                                                           !! the destination grid via the old
                                                                           !! algorithm [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(in)    :: dzInterface  !< Change in interface position
                                                                           !! [H ~> m or kg m-2]
  logical,                         optional, intent(in)    :: debug        !< If true, show the call tree

  ! Local variables

end subroutine ALE_remap_set_h_vel_via_dz
module subroutine ALE_remap_set_h_vel_partial(CS, G, GV, h_mask, h_u, h_v)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_mask       !< Thickness at tracer points
                                                                           !! used to apply the partial
                                                                           !! cell masking [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  ! Local variables

end subroutine ALE_remap_set_h_vel_partial
module subroutine ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness at tracer points of the
                                                                           !! grid being interpolated to velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC          !< Open boundary structure

  ! Local variables

end subroutine ALE_remap_set_h_vel_OBC
module subroutine ALE_remap_velocities(CS, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, u, v, debug, &
                                dt, allow_preserve_variance)
  type(ALE_CS),                              intent(in)    :: CS        !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV        !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_old_u   !< Source grid thickness at zonal
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(in)    :: h_old_v   !< Source grid thickness at meridional
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_new_u   !< Destination grid thickness at zonal
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(in)    :: h_new_v   !< Destination grid thickness at meridional
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: u         !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: v         !< Meridional velocity [L T-1 ~> m s-1]
  logical,                         optional, intent(in)    :: debug     !< If true, show the call tree
  real,                            optional, intent(in)    :: dt        !< time step for diagnostics [T ~> s]
  logical,                         optional, intent(in)    :: allow_preserve_variance !< If true, enables ke-conserving
                                                                                      !! correction

  ! Local variables
                                                 ! 0.5 * rho0 *  u**2 [R Z L2 T-3 ~> W m-2]
                                                 ! 0.5 * rho0 *  v**2 [R Z L2 T-3 ~> W m-2]

end subroutine ALE_remap_velocities
module subroutine ALE_remap_interface_vals(CS, G, GV, h_old, h_new, int_val)
  type(ALE_CS),                              intent(in)    :: CS       !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G        !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV       !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old    !< Thickness of source grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new    !< Thickness of destination grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(inout) :: int_val  !< The interface values to interpolate [A]


end subroutine ALE_remap_interface_vals
module subroutine ALE_remap_vertex_vals(CS, G, GV, h_old, h_new, vert_val)
  type(ALE_CS),                              intent(in)    :: CS       !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G        !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV       !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old    !< Thickness of source grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new    !< Thickness of destination grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                                             intent(inout) :: vert_val  !< The interface values to interpolate [A]


end subroutine ALE_remap_vertex_vals
module subroutine apply_partial_cell_mask(h1, h_mask)
  real, dimension(:), intent(inout) :: h1 !< A column of thicknesses to be masked out after their
                                          !! running vertical sum exceeds h_mask [H ~> m or kg m-2]
  real,               intent(in)    :: h_mask !< The depth after which the thicknesses in h1 are
                                          !! masked out [H ~> m or kg m-2]
  ! Local variables

end subroutine apply_partial_cell_mask
module subroutine mask_near_bottom_vel(vel, h, h_BBL, h_thin, nk)
  integer, intent(in)    :: nk      !< The number of layers in this column
  real,    intent(inout) :: vel(nk) !< The velocity component being zeroed out [L T-1 ~> m s-1]
  real,    intent(in)    :: h(nk)   !< The layer thicknesses at velocity points  [H ~> m or kg m-2]
  real,    intent(in)    :: h_BBL   !< The thickness of the near-bottom region over which to apply
                                    !! the filtering [H ~> m or kg m-2]
  real,    intent(in)    :: h_thin  !< A layer thickness below which the filtering is applied [H ~> m or kg m-2]

  ! Local variables

end subroutine mask_near_bottom_vel
module subroutine ALE_remap_scalar(CS, G, GV, nk_src, h_src, s_src, h_dst, s_dst, all_cells, old_remap)
  type(remapping_CS),                      intent(in)    :: CS        !< Remapping control structure
  type(ocean_grid_type),                   intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                 intent(in)    :: GV        !< Ocean vertical grid structure
  integer,                                 intent(in)    :: nk_src    !< Number of levels on source grid
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: h_src     !< Level thickness of source grid
                                                                      !! [H ~> m or kg m-2] or other units
                                                                      !! if H_neglect is provided
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: s_src     !< Scalar on source grid, in arbitrary units [A]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(in)   :: h_dst     !< Level thickness of destination grid in the
                                                                      !! same units as h_src, often [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(inout) :: s_dst    !< Scalar on destination grid, in the same
                                                                      !! arbitrary units as s_src [A]
  logical, optional,                       intent(in)    :: all_cells !< If false, only reconstruct for
                                                                      !! non-vanished cells. Use all vanished
                                                                      !! layers otherwise (default).
  logical, optional,                       intent(in)    :: old_remap !< If true, use the old "remapping_core_w"
                                                                      !! method, otherwise use "remapping_core_h".
   ! Local variables

end subroutine ALE_remap_scalar
module subroutine TS_PLM_edge_values( CS, S_t, S_b, T_t, T_b, G, GV, tv, h, bdry_extrap )
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),            intent(inout) :: CS   !< module control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer [C ~> degC]
  type(thermo_var_ptrs),   intent(in)    :: tv   !< thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thickness [H ~> m or kg m-2]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells

end subroutine TS_PLM_edge_values
module subroutine ALE_PLM_edge_values( CS, G, GV, h, Q, bdry_extrap, Q_t, Q_b )
  type(ALE_CS),            intent(in)    :: CS   !< module control structure
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Q    !< 3d scalar array, in arbitrary units [A]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_t  !< Scalar at the top edge of each layer [A]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_b  !< Scalar at the bottom edge of each layer [A]
  ! Local variables

end subroutine ALE_PLM_edge_values
module subroutine TS_PPM_edge_values( CS, S_t, S_b, T_t, T_b, G, GV, tv, h, bdry_extrap )
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),            intent(inout) :: CS   !< module control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer [C ~> degC]
  type(thermo_var_ptrs),   intent(in)    :: tv   !< thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thicknesses [H ~> m or kg m-2]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells

  ! Local variables

end subroutine TS_PPM_edge_values
module subroutine TS_PLM_WLS_edge_values(CS, S_t, S_b, T_t, T_b, G, GV, tv, h)
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),            intent(inout) :: CS   !< module control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer [C ~> degC]
  type(thermo_var_ptrs),   intent(in)    :: tv   !< thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thickness [H ~> m or kg m-2]
  ! Local variables

end subroutine TS_PLM_WLS_edge_values
module subroutine ALE_initRegridding(G, GV, US, max_depth, param_file, mdl, regridCS)
  type(ocean_grid_type),   intent(in)  :: G          !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  real,                    intent(in)  :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(param_file_type),   intent(in)  :: param_file !< parameter file
  character(len=*),        intent(in)  :: mdl        !< Name of calling module
  type(regridding_CS),     intent(out) :: regridCS   !< Regridding parameters and work arrays
  ! Local variables

end subroutine ALE_initRegridding
module function ALE_getCoordinate( CS )
  type(ALE_CS), pointer    :: CS                  !< module control structure

  real, dimension(CS%nk+1) :: ALE_getCoordinate !< The coordinate positions, in the appropriate units
                                                !! of the target coordinate, e.g. [Z ~> m] for z*,
                                                !! non-dimensional for sigma, etc.
end function ALE_getCoordinate
module function ALE_getCoordinateUnits( CS )
  type(ALE_CS), pointer :: CS   !< module control structure

  character(len=20)     :: ALE_getCoordinateUnits

end function ALE_getCoordinateUnits
logical module function ALE_remap_init_conds( CS )
  type(ALE_CS), pointer :: CS   !< module control structure

end function ALE_remap_init_conds
module subroutine ALE_update_regrid_weights( dt, CS )
  real,         intent(in) :: dt !< Time-step used between ALE calls [T ~> s]
  type(ALE_CS), pointer    :: CS !< ALE control structure
  ! Local variables

end subroutine ALE_update_regrid_weights
module subroutine ALE_updateVerticalGridType(CS, GV)
  type(ALE_CS),            pointer :: CS  !< ALE control structure
  type(verticalGrid_type), pointer :: GV  !< vertical grid information


end subroutine ALE_updateVerticalGridType
module subroutine ALE_writeCoordinateFile( CS, GV, directory )
  type(ALE_CS),            pointer     :: CS         !< module control structure
  type(verticalGrid_type), intent(in)  :: GV         !< ocean vertical grid structure
  character(len=*),        intent(in)  :: directory  !< directory for writing grid info


end subroutine ALE_writeCoordinateFile
module subroutine ALE_initThicknessToCoord( CS, G, GV, h, height_units )
  type(ALE_CS), intent(inout)                            :: CS  !< module control structure
  type(ocean_grid_type), intent(in)                      :: G   !< module grid structure
  type(verticalGrid_type), intent(in)                    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h   !< layer thickness in thickness units
                                                                !! [H ~> m or kg m-2] or height units [Z ~> m]
  logical,                          optional, intent(in) :: height_units !< If present and true, the
                                                                !! thicknesses are in height units

  ! Local variables

end subroutine ALE_initThicknessToCoord
  end interface

end module MOM_ALE
