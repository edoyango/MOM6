! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization functions for state variables, u, v, h, T and S.
module MOM_state_initialization

use MOM_debugging, only : hchksum, qchksum, uvchksum
use MOM_density_integrals, only : int_specific_vol_dp
use MOM_density_integrals, only : find_depth_of_pressure_in_cell
use MOM_coms, only : max_across_PEs, min_across_PEs, reproducing_sum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only :  CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : find_eta, dz_to_thickness, dz_to_thickness_simple
use MOM_interface_heights, only : calc_derived_thermo
use MOM_io, only : file_exists, field_size, MOM_read_data, MOM_read_vector, slasher
use MOM_open_boundary, only : ocean_OBC_type, open_boundary_test_extern_h
use MOM_open_boundary, only : fill_temp_salt_segments, setup_OBC_tracer_reservoirs
use MOM_open_boundary, only : fill_thickness_segments
use MOM_open_boundary, only : set_initialized_OBC_tracer_reservoirs
use MOM_restart, only : restore_state, is_new_run, copy_restart_var, copy_restart_vector
use MOM_restart, only : restart_registry_lock, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, set_up_sponge_ML_density
use MOM_sponge, only : initialize_sponge, sponge_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, set_up_ALE_sponge_vel_field
use MOM_ALE_sponge, only : ALE_sponge_CS, initialize_ALE_sponge
use MOM_string_functions, only : uppercase, lowercase
use MOM_time_manager, only : time_type, operator(/=)
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : setVerticalGridAxes, verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type, EOS_domain
use MOM_EOS, only : convert_temp_salt_for_TEOS10
use user_initialization, only : user_initialize_thickness, user_initialize_velocity
use user_initialization, only : user_init_temperature_salinity, user_set_OBC_data
use user_initialization, only : user_initialize_sponges
use DOME_initialization, only : DOME_initialize_thickness
use DOME_initialization, only : DOME_set_OBC_data
use DOME_initialization, only : DOME_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_thickness
use ISOMIP_initialization, only : ISOMIP_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_temperature_salinity
use RGC_initialization, only : RGC_initialize_sponges
use baroclinic_zone_initialization, only : baroclinic_zone_init_temperature_salinity
use benchmark_initialization, only : benchmark_initialize_thickness
use benchmark_initialization, only : benchmark_init_temperature_salinity
use Neverworld_initialization, only : Neverworld_initialize_thickness
use circle_obcs_initialization, only : circle_obcs_initialize_thickness
use lock_exchange_initialization, only : lock_exchange_initialize_thickness
use external_gwave_initialization, only : external_gwave_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_temperature_salinity
use DOME2d_initialization, only : DOME2d_initialize_sponges
use adjustment_initialization, only : adjustment_initialize_thickness
use adjustment_initialization, only : adjustment_initialize_temperature_salinity
use sloshing_initialization, only : sloshing_initialize_thickness
use sloshing_initialization, only : sloshing_initialize_temperature_salinity
use seamount_initialization, only : seamount_initialize_thickness
use seamount_initialization, only : seamount_initialize_temperature_salinity
use dumbbell_initialization, only : dumbbell_initialize_thickness
use dumbbell_initialization, only : dumbbell_initialize_temperature_salinity
use Phillips_initialization, only : Phillips_initialize_thickness
use Phillips_initialization, only : Phillips_initialize_velocity
use Phillips_initialization, only : Phillips_initialize_sponges
use Rossby_front_2d_initialization, only : Rossby_front_initialize_thickness
use Rossby_front_2d_initialization, only : Rossby_front_initialize_temperature_salinity
use Rossby_front_2d_initialization, only : Rossby_front_initialize_velocity
use SCM_CVMix_tests, only: SCM_CVMix_tests_TS_init
use dyed_channel_initialization, only : dyed_channel_set_OBC_tracer_data
use dyed_obcs_initialization, only : dyed_obcs_set_OBC_data
use supercritical_initialization, only : supercritical_set_OBC_data
use soliton_initialization, only : soliton_initialize_velocity
use soliton_initialization, only : soliton_initialize_thickness
use BFB_initialization, only : BFB_initialize_sponges_southonly
use dense_water_initialization, only : dense_water_initialize_TS
use dense_water_initialization, only : dense_water_initialize_sponges
use dumbbell_initialization, only : dumbbell_initialize_sponges
use MOM_tracer_Z_init, only : tracer_Z_init_array, determine_temperature
use MOM_ALE, only : ALE_initRegridding, ALE_CS, ALE_initThicknessToCoord
use MOM_ALE, only : ALE_remap_scalar, ALE_regrid_accelerated, TS_PLM_edge_values
use MOM_regridding, only : regridding_CS, set_regrid_params, getCoordinateResolution
use MOM_regridding, only : regridding_main, regridding_preadjust_reqs, convective_adjustment
use MOM_regridding, only : set_dz_neglect, set_h_neglect
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer, homogenize_field
use MOM_oda_incupd, only: oda_incupd_CS, initialize_oda_incupd_fixed, initialize_oda_incupd
use MOM_oda_incupd, only: set_up_oda_incupd_field, set_up_oda_incupd_vel_field
use MOM_oda_incupd, only: calc_oda_increments, output_oda_incupd_inc

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize_state, MOM_initialize_OBCs

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40)  :: mdl = "MOM_state_initialization" !< This module's name.


  interface
module subroutine MOM_initialize_state(u, v, h, tv, Time, G, GV, US, PF, dirs, &
                                restart_CS, ALE_CSp, tracer_Reg, sponge_CSp, &
                                ALE_sponge_CSp, oda_incupd_CSp, OBC_for_remap, &
                                Time_in, frac_shelf_h, mass_shelf, OBC_for_bug)
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(out)   :: u    !< The zonal velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(out)   :: v    !< The meridional velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(out)   :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),      intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                    !! variables
  type(time_type),            intent(inout) :: Time !< Time at the start of the run segment.
  type(param_file_type),      intent(in)    :: PF   !< A structure indicating the open file to parse
                                                    !! for model parameter values.
  type(directories),          intent(in)    :: dirs !< A structure containing several relevant
                                                    !! directory paths.
  type(MOM_restart_CS),       intent(inout) :: restart_CS !< MOM restart control structure
  type(ALE_CS),               pointer       :: ALE_CSp !< The ALE control structure for remapping
  type(tracer_registry_type), pointer       :: tracer_Reg !< A pointer to the tracer registry
  type(sponge_CS),            pointer       :: sponge_CSp !< The layerwise sponge control structure.
  type(ALE_sponge_CS),        pointer       :: ALE_sponge_CSp !< The ALE sponge control structure.
  type(ocean_OBC_type),       pointer       :: OBC_for_remap !< The open boundary condition control
                                                    !! structure that may be used for remapping velocities.
                                                    !! This must be on the unrotated grid, but only the
                                                    !! position and directions of the OBC faces are used.
  type(oda_incupd_CS),        pointer       :: oda_incupd_CSp !< The oda_incupd control structure.
  type(time_type), optional,  intent(in)    :: Time_in !< Time at the start of the run segment.
  real, dimension(SZI_(G),SZJ_(G)), &
                     optional, intent(in)   :: frac_shelf_h    !< The fraction of the grid cell covered
                                                               !! by a floating ice shelf [nondim].
  real, dimension(SZI_(G),SZJ_(G)), &
                     optional, intent(in)   :: mass_shelf      !< The mass per unit area of the overlying
                                                               !! ice shelf [R Z ~> kg m-2]
  type(ocean_OBC_type), optional, pointer   :: OBC_for_bug  !< An open boundary condition control structure
                                                    !! that might be used to store OBC temperatures and
                                                    !! salinities if OBC_RESERVOIR_INIT_BUG is true.
  ! Local variables

                         ! run from the interior tracer concentrations regardless of properties that
                         ! may be explicitly specified for the reservoir concentrations.
                         ! by a large surface pressure by squeezing the column.
                         ! by a large surface pressure, such as with an ice sheet.
                        ! is a run from a restart file; this option
                        ! allows the use of Fatal unused parameters.
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  ! This include declares and sets the variable "version".

end subroutine MOM_initialize_state
module subroutine MOM_initialize_OBCs(h, tv, OBC, Time, G, GV, US, PF, restart_CS, tracer_Reg)
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),      intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                    !! variables
  type(ocean_OBC_type),       pointer       :: OBC   !< The open boundary condition control structure.
  type(time_type),            intent(in)    :: Time !< Time at the start of the run segment.
  type(param_file_type),      intent(in)    :: PF   !< A structure indicating the open file to parse
                                                    !! for model parameter values.
  type(MOM_restart_CS),       intent(inout) :: restart_CS !< MOM restart control structure
  type(tracer_registry_type), pointer       :: tracer_Reg !< A pointer to the tracer registry

  ! Local variables
                          ! recreate the bugs, or if false bugs are only used if actively selected.
                        ! of the open boundary condition code.
                        ! run from the interior tracer concentrations regardless of properties that
                        ! may be explicitly specified for the reservoir concentrations.

end subroutine MOM_initialize_OBCs
module subroutine initialize_thickness_from_file(h, depth_tot, G, GV, US, param_file, file_has_thickness, &
                                          just_read, mass_file)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h    !< The thickness that is being initialized, in height
                                               !! or thickness units, depending on the value of
                                               !! mass_file [Z ~> m] or [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,                 intent(in)  :: file_has_thickness !< If true, this file contains layer
                                               !! thicknesses; otherwise it contains
                                               !! interface heights.
  logical,                 intent(in)  :: just_read !< If true, this call will only read
                                               !! parameters without changing h.
  logical,                 intent(in)  :: mass_file !< If true, this file contains layer thicknesses in
                                               !! units of mass per unit area.

  ! Local variables
                      ! file to convert it to units of m [various]
                      ! them to units of m or correct sign conventions to positive upward [various]
                      ! thickness to fit the bathymetry [Z ~> m].
                      ! correct_thickness is false [Z ~> m]

end subroutine initialize_thickness_from_file
module subroutine adjustEtaToFitBathymetry(G, GV, US, eta, h, ht, dZ_ref_eta)
  type(ocean_grid_type),                       intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: eta !< Interface heights [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: h   !< Layer thicknesses [Z ~> m]
  real,                                        intent(in)    :: ht  !< Tolerance to exceed adjustment
                                                                    !! criteria [Z ~> m]
  real,                              optional, intent(in)    :: dZ_ref_eta !< The difference between the
                                                                    !! reference heights for bathyT and
                                                                    !! eta [Z ~> m], 0 by default.
  ! Local variables

end subroutine adjustEtaToFitBathymetry
module subroutine initialize_thickness_uniform(h, depth_tot, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                          ! negative because it is positive upward.
                          ! positive upward [Z ~> m].

end subroutine initialize_thickness_uniform
module subroutine initialize_thickness_list(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                          ! usually negative because it is positive upward.
                          ! positive upward, in depth units [Z ~> m].

end subroutine initialize_thickness_list
module subroutine initialize_thickness_param(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                          ! negative because it is positive upward.
                          ! positive upward [Z ~> m].

end subroutine initialize_thickness_param
module subroutine initialize_thickness_search
end subroutine initialize_thickness_search
module subroutine depress_surface(h, G, GV, US, param_file, tv, just_read, z_top_shelf)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                      !! parameters without changing h.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: z_top_shelf    !< Top interface position under ice shelf [Z ~> m]
  ! Local variables
                       ! which can be used to change units, for example, often [Z m-1 ~> 1].

end subroutine depress_surface
module subroutine trim_for_ice(PF, G, GV, US, ALE_CSp, tv, h, just_read)
  type(param_file_type),   intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(ALE_CS),            pointer       :: ALE_CSp !< ALE control structure
  type(thermo_var_ptrs),   intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                                                        ! of salinity within each layer [S ~> ppt]
                                                        ! of temperature within each layer [C ~> degC]
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.

end subroutine trim_for_ice
module subroutine calc_sfc_displacement(PF, G, GV, US, mass_shelf, tv, h)
  type(param_file_type),   intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: mass_shelf  !< Ice shelf mass [R Z ~> kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]

  ! temporary arrays
                    ! mass [R Z ~> kg m-2]

end subroutine calc_sfc_displacement
module subroutine cut_off_column_top(nk, tv, GV, US, G_earth, depth, min_thickness, T, T_t, T_b, &
                              S, S_t, S_b, p_surf, h, remap_CS, z_tol, frac_dp_bugfix)
  integer,               intent(in)    :: nk  !< Number of layers
  type(thermo_var_ptrs), intent(in)    :: tv  !< Thermodynamics structure
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  real,                  intent(in)    :: G_earth !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                  intent(in)    :: depth !< Depth of ocean column [Z ~> m].
  real,                  intent(in)    :: min_thickness !< Smallest thickness allowed [H ~> m or kg m-2].
  real, dimension(nk),   intent(inout) :: T   !< Layer mean temperature [C ~> degC]
  real, dimension(nk),   intent(in)    :: T_t !< Temperature at top of layer [C ~> degC]
  real, dimension(nk),   intent(in)    :: T_b !< Temperature at bottom of layer [C ~> degC]
  real, dimension(nk),   intent(inout) :: S   !< Layer mean salinity [S ~> ppt]
  real, dimension(nk),   intent(in)    :: S_t !< Salinity at top of layer [S ~> ppt]
  real, dimension(nk),   intent(in)    :: S_b !< Salinity at bottom of layer [S ~> ppt]
  real,                  intent(in)    :: p_surf !< Imposed pressure on ocean at surface [R L2 T-2 ~> Pa]
  real, dimension(nk),   intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(remapping_CS),    pointer       :: remap_CS !< Remapping structure for remapping T and S,
                                                   !! if associated
  real,                  intent(in)    :: z_tol !< The tolerance with which to find the depth
                                                !! matching the specified pressure [Z ~> m].
  logical,               intent(in)    :: frac_dp_bugfix !< If true, use bugfix in frac_dp_at_pos

  ! Local variables

  ! Keep a copy of the initial thicknesses in reverse order to use in remapping
end subroutine cut_off_column_top
module subroutine initialize_velocity_from_file(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables

end subroutine initialize_velocity_from_file
module subroutine initialize_velocity_zero(u, v, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables

end subroutine initialize_velocity_zero
module subroutine initialize_velocity_uniform(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables

end subroutine initialize_velocity_uniform
module subroutine initialize_velocity_circular(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables
end subroutine initialize_velocity_circular
module subroutine initialize_temp_salt_from_file(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                                   intent(in)  :: just_read !< If true, this call will only
                                                           !! read parameters without changing T or S.
  ! Local variables

end subroutine initialize_temp_salt_from_file
module subroutine initialize_temp_salt_from_profile(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                                   intent(in)  :: just_read !< If true, this call will only read
                                                               !! parameters without changing T or S.
  ! Local variables

end subroutine initialize_temp_salt_from_profile
module subroutine initialize_temp_salt_fit(T, S, G, GV, US, param_file, eqn_of_state, P_Ref, just_read)
  type(ocean_grid_type),   intent(in)  :: G            !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV           !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< The potential temperature that is
                                                       !! being initialized [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< The salinity that is being
                                                       !! initialized [S ~> ppt].
  type(unit_scale_type),   intent(in)  :: US           !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file   !< A structure to parse for run-time
                                                       !! parameters.
  type(EOS_type),          intent(in)  :: eqn_of_state !< Equation of state structure
  real,                    intent(in)  :: P_Ref        !< The coordinate-density reference pressure
                                                       !! [R L2 T-2 ~> Pa].
  logical,                 intent(in)  :: just_read    !< If true, this call will only read
                                                       !! parameters without changing T or S.
  ! Local variables
end subroutine initialize_temp_salt_fit
module subroutine initialize_temp_salt_linear(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for
                                                               !! run-time parameters
  logical,                                   intent(in)  :: just_read !< If present and true,
                                                               !! this call will only read parameters
                                                               !! without changing T or S.

  ! Local variables

end subroutine initialize_temp_salt_linear
module subroutine initialize_sponges_file(G, GV, US, use_temperature, tv, u, v, depth_tot, param_file, &
                                   Layer_CSp, ALE_CSp, Time)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  logical,                 intent(in) :: use_temperature !< If true, T & S are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic
                                              !! variables.
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity that is being
                                              !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity that is being
                                              !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(sponge_CS),         pointer    :: Layer_CSp  !< A pointer that is set to point to the control
                                              !! structure for this module (in layered mode).
  type(ALE_sponge_CS),     pointer    :: ALE_CSp  !< A pointer that is set to point to the control
                                                  !! structure for this module (in ALE mode).
  type(time_type),         intent(in) :: Time !< Time at the start of the run segment. Time_in
                                              !! overrides any value set for Time.
  ! Local variables

                                    ! on the vertical grid of the input file  [C ~> degC]
                                    ! on the vertical grid of the input file [S ~> ppt]
                                    ! velocities on the vertical grid of the input file [L T-1 ~> m s-1]
                                    ! velocities on the vertical grid of the input file [L T-1 ~> m s-1]



                              ! the horizontal dimension and in time prior to vertical remapping.

end subroutine initialize_sponges_file
module subroutine initialize_oda_incupd_file(G, GV, US, use_temperature, tv, h, u, v, param_file, &
                                      oda_incupd_CSp, restart_CS, Time)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  logical,                 intent(in)    :: use_temperature !< If true, T & S are state variables.
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic
                                                 !! variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: u     !< The zonal velocity that is being
                                                 !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             intent(in) :: v     !< The meridional velocity that is being
                                                 !! initialized [L T-1 ~> m s-1]
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(oda_incupd_CS),     pointer    :: oda_incupd_CSp  !< A pointer that is set to point to the control
                                                 !! structure for this module.
  type(MOM_restart_CS),    intent(in) :: restart_CS !< MOM restart control structure
  type(time_type),         intent(in) :: Time    !< Time at the start of the run segment. Time_in
                                                 !! overrides any value set for Time.
  ! Local variables
                                    ! on the vertical grid of the input file, used for both
                                    ! temperatures [C ~> degC] and salinities [S ~> ppt]
                                    ! increments on the vertical grid of the input file [L T-1 ~> m s-1]
                                    ! increments on the vertical grid of the input file [L T-1 ~> m s-1]




!  logical :: use_ALE ! True if ALE is being used, False if in layered mode

end subroutine initialize_oda_incupd_file
module subroutine set_velocity_depth_max(G)
  type(ocean_grid_type), intent(inout) :: G !< The ocean's grid structure
  ! Local variables

end subroutine set_velocity_depth_max
module subroutine set_velocity_depth_min(G)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  ! Local variables

end subroutine set_velocity_depth_min
module subroutine MOM_temp_salt_initialize_from_Z(h, tv, depth_tot, G, GV, US, PF, just_read, frac_shelf_h)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h    !< Layer thicknesses being initialized [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                 !! variables including temperature and salinity
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                 !! parameters without changing T or S.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: frac_shelf_h  !< The fraction of the grid cell covered
                                                 !! by a floating ice shelf [nondim].

  ! Local variables
                                   !! and salinity in z-space; by default it is also used for ice shelf area.

  ! This include declares and sets the variable "version".



                      ! interpolation from an input dataset [C ~> degC]
                      ! interpolation from an input dataset [S ~> ppt]
                         ! thickness to fit the bathymetry [Z ~> m].
                         ! correct_thickness is false [Z ~> m]


  ! data arrays
                                                ! relative to the surface [Z ~> m].

  ! Local variables for ALE remapping
                                    ! regridding [H ~> m or kg m-2]
                                    ! remapping cell reconstructions [Z ~> m]
                                    ! remapping edge value calculations [Z ~> m]

                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
                                  ! for horizontal regridding.  Values below 20190101 recover the
                                  ! answers from 2018, while higher values use expressions that have
                                  ! been rearranged for rotational invariance.
                                   ! extrapolating the densities at the bottom of unstable profiles
                                   ! from data when finding the initial interface locations in
                                   ! layered mode from a dataset of T and S.

end subroutine MOM_temp_salt_initialize_from_Z
module subroutine find_interfaces(rho, zin, nk_data, Rb, Z_bot, zi, G, GV, US, nlevs, nkml, hml, &
                           eps_z, eps_rho, density_extrap_bug)
  type(ocean_grid_type),      intent(in)  :: G     !< The ocean's grid structure
  type(verticalGrid_type),    intent(in)  :: GV    !< The ocean's vertical grid structure
  integer,                    intent(in)  :: nk_data !< The number of levels in the input data
  real, dimension(SZI_(G),SZJ_(G),nk_data), &
                              intent(in)  :: rho   !< Potential density in z-space [R ~> kg m-3]
  real, dimension(nk_data),   intent(in)  :: zin   !< Input data levels [Z ~> m].
  real, dimension(SZK_(GV)+1), intent(in) :: Rb    !< target interface densities [R ~> kg m-3]
  real, dimension(SZI_(G),SZJ_(G)), &
                              intent(in)  :: Z_bot !< The (usually negative) height of the seafloor
                                                   !! relative to the surface [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(out) :: zi    !< The returned interface heights [Z ~> m]
  type(unit_scale_type),      intent(in)  :: US    !< A dimensional unit scaling type
  integer, dimension(SZI_(G),SZJ_(G)), &
                              intent(in)  :: nlevs !< number of valid points in each column
  integer,                    intent(in)  :: nkml  !< number of mixed layer pieces to distribute over
                                                   !! a depth of hml.
  real,                       intent(in)  :: hml   !< mixed layer depth [Z ~> m].
  real,                       intent(in)  :: eps_z !< A negligibly small layer thickness [Z ~> m].
  real,                       intent(in)  :: eps_rho !< A negligibly small density difference [R ~> kg m-3].
  logical,                    intent(in)  :: density_extrap_bug !< If true use an expression with an
                                                   !! indexing bug for projecting the densities at
                                                   !! the bottom of unstable profiles from data when
                                                   !! finding the initial interface locations in
                                                   !! layered mode from a dataset of T and S.

  ! Local variables

end subroutine find_interfaces
module subroutine MOM_state_init_tests(G, GV, US, tv)
  type(ocean_grid_type),     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),     intent(in)    :: tv   !< Thermodynamics structure.

  ! Local variables

end subroutine MOM_state_init_tests
  end interface

end module MOM_state_initialization
