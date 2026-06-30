! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The central module of the MOM6 ocean model
module MOM

! Infrastructure modules
use MOM_array_transform,      only : rotate_array, rotate_vector
use MOM_debugging,            only : MOM_debugging_init, hchksum, uvchksum, totalTandS
use MOM_debugging,            only : check_redundant, query_debugging_checks
use MOM_checksum_packages,    only : MOM_thermo_chksum, MOM_state_chksum
use MOM_checksum_packages,    only : MOM_accel_chksum, MOM_surface_chksum
use MOM_coms,                 only : num_PEs
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,            only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,            only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,        only : diag_mediator_init, enable_averaging, enable_averages
use MOM_diag_mediator,        only : diag_mediator_infrastructure_init, diag_mediator_set_OBC_info
use MOM_diag_mediator,        only : diag_set_state_ptrs, diag_update_remap_grids
use MOM_diag_mediator,        only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator,        only : register_diag_field, register_cell_measure
use MOM_diag_mediator,        only : set_axes_info, diag_ctrl, diag_masks_set
use MOM_diag_mediator,        only : set_masks_for_axes
use MOM_diag_mediator,        only : diag_grid_storage, diag_grid_storage_init
use MOM_diag_mediator,        only : diag_save_grids, diag_restore_grids
use MOM_diag_mediator,        only : diag_copy_storage_to_diag, diag_copy_diag_to_storage
use MOM_domains,              only : MOM_domains_init, MOM_domain_type
use MOM_domains,              only : sum_across_PEs, pass_var, pass_vector
use MOM_domains,              only : clone_MOM_domain, deallocate_MOM_domain
use MOM_domains,              only : To_North, To_East, To_South, To_West
use MOM_domains,              only : To_All, Omit_corners, CGRID_NE, SCALAR_PAIR
use MOM_domains,              only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,              only : start_group_pass, complete_group_pass, Omit_Corners
use MOM_error_handler,        only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,        only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,        only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,          only : read_param, get_param, log_version, param_file_type
use MOM_forcing_type,         only : forcing, mech_forcing, find_ustar
use MOM_forcing_type,         only : MOM_forcing_chksum, MOM_mech_forcing_chksum
use MOM_get_input,            only : Get_MOM_Input, directories
use MOM_io,                   only : MOM_io_init, vardesc, var_desc
use MOM_io,                   only : slasher, file_exists, MOM_read_data
use MOM_obsolete_params,      only : find_obsolete_params
use MOM_restart,              only : register_restart_field, register_restart_pair, save_restart
use MOM_restart,              only : query_initialized, set_initialized, restart_registry_lock
use MOM_restart,              only : restart_init, is_new_run, determine_is_new_run, MOM_restart_CS
use MOM_spatial_means,        only : global_mass_integral
use MOM_time_manager,         only : time_type, real_to_time, operator(+)
use MOM_time_manager,         only : operator(-), operator(>), operator(*), operator(/)
use MOM_time_manager,         only : operator(>=), operator(==), increment_date
use MOM_unit_tests,           only : unit_tests

! MOM core modules
use MOM_ALE,                   only : ALE_init, ALE_end, ALE_regrid, ALE_CS, adjustGridForIntegrity
use MOM_ALE,                   only : ALE_getCoordinate, ALE_getCoordinateUnits, ALE_writeCoordinateFile
use MOM_ALE,                   only : ALE_updateVerticalGridType, ALE_remap_init_conds, pre_ALE_adjustments
use MOM_ALE,                   only : ALE_remap_tracers, ALE_remap_velocities
use MOM_ALE,                   only : ALE_remap_set_h_vel, ALE_remap_set_h_vel_via_dz
use MOM_ALE,                   only : ALE_update_regrid_weights, pre_ALE_diagnostics, ALE_register_diags
use MOM_ALE,                   only : ALE_set_extrap_boundaries
use MOM_ALE_sponge,            only : rotate_ALE_sponge, update_ALE_sponge_field
use MOM_barotropic,            only : Barotropic_CS
use MOM_boundary_update,       only : call_OBC_register, OBC_register_end, update_OBC_CS
use MOM_check_scaling,         only : check_MOM6_scaling_factors
use MOM_coord_initialization,  only : MOM_initialize_coord, write_vertgrid_file
use MOM_diabatic_driver,       only : diabatic, diabatic_driver_init, diabatic_CS, extract_diabatic_member
use MOM_diabatic_driver,       only : adiabatic, adiabatic_driver_init, diabatic_driver_end
use MOM_diabatic_driver,       only : register_diabatic_restarts
use MOM_stochastics,           only : stochastics_init, update_stochastics, stochastic_CS, apply_skeb
use MOM_diagnostics,           only : calculate_diagnostic_fields, MOM_diagnostics_init
use MOM_diagnostics,           only : register_transport_diags, post_transport_diagnostics
use MOM_diagnostics,           only : register_surface_diags, write_static_fields
use MOM_diagnostics,           only : post_surface_dyn_diags, post_surface_thermo_diags
use MOM_diagnostics,           only : diagnostics_CS, surface_diag_IDs, transport_diag_IDs
use MOM_diagnostics,           only : MOM_diagnostics_end
use MOM_dynamics_unsplit,      only : step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
use MOM_dynamics_unsplit,      only : initialize_dyn_unsplit, end_dyn_unsplit
use MOM_dynamics_unsplit,      only : MOM_dyn_unsplit_CS
use MOM_dynamics_split_RK2,    only : step_MOM_dyn_split_RK2, register_restarts_dyn_split_RK2
use MOM_dynamics_split_RK2,    only : initialize_dyn_split_RK2, end_dyn_split_RK2
use MOM_dynamics_split_RK2,    only : MOM_dyn_split_RK2_CS, remap_dyn_split_rk2_aux_vars
use MOM_dynamics_split_RK2,    only : init_dyn_split_RK2_diabatic
use MOM_dynamics_split_RK2b,   only : step_MOM_dyn_split_RK2b, register_restarts_dyn_split_RK2b
use MOM_dynamics_split_RK2b,   only : initialize_dyn_split_RK2b, end_dyn_split_RK2b
use MOM_dynamics_split_RK2b,   only : MOM_dyn_split_RK2b_CS, remap_dyn_split_RK2b_aux_vars
use MOM_dynamics_unsplit_RK2,  only : step_MOM_dyn_unsplit_RK2, register_restarts_dyn_unsplit_RK2
use MOM_dynamics_unsplit_RK2,  only : initialize_dyn_unsplit_RK2, end_dyn_unsplit_RK2
use MOM_dynamics_unsplit_RK2,  only : MOM_dyn_unsplit_RK2_CS
use MOM_dyn_horgrid,           only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_dyn_horgrid,           only : rotate_dyn_horgrid
use MOM_EOS,                   only : EOS_init, calculate_density, calculate_TFreeze, EOS_domain
use MOM_fixed_initialization,  only : MOM_initialize_fixed
use MOM_forcing_type,          only : allocate_forcing_type, allocate_mech_forcing
use MOM_forcing_type,          only : deallocate_mech_forcing, deallocate_forcing_type
use MOM_forcing_type,          only : rotate_forcing, rotate_mech_forcing
use MOM_forcing_type,          only : copy_common_forcing_fields, set_derived_forcing_fields
use MOM_forcing_type,          only : homogenize_forcing, homogenize_mech_forcing
use MOM_grid,                  only : ocean_grid_type, MOM_grid_init, MOM_grid_end
use MOM_grid,                  only : set_first_direction
use MOM_harmonic_analysis,     only : HA_accum, harmonic_analysis_CS
use MOM_hor_index,             only : hor_index_type, hor_index_init
use MOM_hor_index,             only : rotate_hor_index
use MOM_interface_heights,     only : find_eta, calc_derived_thermo, thickness_to_dz
use MOM_interface_filter,      only : interface_filter, interface_filter_init, interface_filter_end
use MOM_interface_filter,      only : interface_filter_CS
use MOM_internal_tides,        only : int_tide_CS
use MOM_kappa_shear,           only : kappa_shear_at_vertex
use MOM_lateral_mixing_coeffs, only : calc_slope_functions, VarMix_init, VarMix_end
use MOM_lateral_mixing_coeffs, only : calc_resoln_function, calc_depth_function, VarMix_CS
use MOM_MEKE,                  only : MEKE_alloc_register_restart, step_forward_MEKE
use MOM_MEKE,                  only : MEKE_CS, MEKE_init, MEKE_end
use MOM_MEKE_types,            only : MEKE_type
use MOM_mixed_layer_restrat,   only : mixedlayer_restrat, mixedlayer_restrat_init, mixedlayer_restrat_CS
use MOM_mixed_layer_restrat,   only : mixedlayer_restrat_register_restarts
use MOM_obsolete_diagnostics,  only : register_obsolete_diagnostics
use MOM_open_boundary,         only : ocean_OBC_type, open_boundary_end
use MOM_open_boundary,         only : register_temp_salt_segments, update_segment_tracer_reservoirs
use MOM_open_boundary,         only : read_OBC_segment_data, initialize_OBC_segment_reservoirs
use MOM_open_boundary,         only : setup_OBC_tracer_reservoirs
use MOM_open_boundary,         only : setup_OBC_thickness_reservoirs
use MOM_open_boundary,         only : open_boundary_register_restarts, remap_OBC_fields
use MOM_open_boundary,         only : open_boundary_setup_vert, initialize_segment_data
use MOM_open_boundary,         only : update_OBC_segment_data, rotate_OBC_config
use MOM_open_boundary,         only : open_boundary_halo_update, write_OBC_info, chksum_OBC_segments
use MOM_open_boundary,         only : segment_thickness_reservoir_init
use MOM_porous_barriers,       only : porous_widths_layer, porous_widths_interface, porous_barriers_init
use MOM_porous_barriers,       only : porous_barrier_CS
use MOM_set_visc,              only : set_viscous_BBL, set_viscous_ML, set_visc_CS
use MOM_set_visc,              only : set_visc_register_restarts, remap_vertvisc_aux_vars
use MOM_set_visc,              only : set_visc_init, set_visc_end
use MOM_shared_initialization, only : write_ocean_geometry_file
use MOM_sponge,                only : init_sponge_diags, sponge_CS
use MOM_state_initialization,  only : MOM_initialize_state, MOM_initialize_OBCs
use MOM_stoch_eos,             only : MOM_stoch_eos_init, MOM_stoch_eos_run, MOM_stoch_eos_CS
use MOM_stoch_eos,             only : stoch_EOS_register_restarts, post_stoch_EOS_diags, mom_calc_varT
use MOM_sum_output,            only : write_energy, accumulate_net_input
use MOM_sum_output,            only : MOM_sum_output_init, MOM_sum_output_end
use MOM_sum_output,            only : sum_output_CS
use MOM_ALE_sponge,            only : init_ALE_sponge_diags, ALE_sponge_CS
use MOM_thickness_diffuse,     only : thickness_diffuse, thickness_diffuse_init
use MOM_thickness_diffuse,     only : thickness_diffuse_end, thickness_diffuse_CS
use MOM_tracer_advect,         only : advect_tracer, tracer_advect_init
use MOM_tracer_advect,         only : tracer_advect_end, tracer_advect_CS
use MOM_tracer_hor_diff,       only : tracer_hordiff, tracer_hor_diff_init
use MOM_tracer_hor_diff,       only : tracer_hor_diff_end, tracer_hor_diff_CS
use MOM_tracer_registry,       only : tracer_registry_type, register_tracer, tracer_registry_init
use MOM_tracer_registry,       only : register_tracer_diagnostics, post_tracer_diagnostics_at_sync
use MOM_tracer_registry,       only : post_tracer_transport_diagnostics, MOM_tracer_chksum
use MOM_tracer_registry,       only : preALE_tracer_diagnostics, postALE_tracer_diagnostics
use MOM_tracer_registry,       only : lock_tracer_registry, tracer_registry_end
use MOM_tracer_flow_control,   only : call_tracer_register, tracer_flow_control_CS
use MOM_tracer_flow_control,   only : tracer_flow_control_init, call_tracer_surface_state
use MOM_tracer_flow_control,   only : tracer_flow_control_end, call_tracer_register_obc_segments
use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_unit_scaling,          only : unit_scale_type, unit_scaling_init, unit_scaling_end
use MOM_variables,             only : surface, allocate_surface_state, deallocate_surface_state
use MOM_variables,             only : thermo_var_ptrs, vertvisc_type, porous_barrier_type
use MOM_variables,             only : accel_diag_ptrs, cont_diag_ptrs, ocean_internal_state
use MOM_variables,             only : rotate_surface_state
use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd
use MOM_verticalGrid,          only : get_thickness_units, get_flux_units, get_tr_flux_units
use MOM_wave_interface,        only : wave_parameters_CS, waves_end, waves_register_restarts
use MOM_wave_interface,        only : Update_Stokes_Drift

! Database client used for machine-learning interface
use MOM_database_comms,       only : dbcomms_CS_type, database_comms_init, dbclient_type

! ODA modules
use MOM_oda_driver_mod,        only : ODA_CS, oda, init_oda, oda_end
use MOM_oda_driver_mod,        only : set_prior_tracer, set_analysis_time, apply_oda_tracer_increments
use MOM_oda_incupd,            only : oda_incupd_CS, init_oda_incupd_diags

! Offline modules
use MOM_offline_main,          only : offline_transport_CS, offline_transport_init, update_offline_fields
use MOM_offline_main,          only : insert_offline_main, extract_offline_main, post_offline_convergence_diags
use MOM_offline_main,          only : register_diags_offline_transport, offline_advection_ale
use MOM_offline_main,          only : offline_redistribute_residual, offline_diabatic_ale
use MOM_offline_main,          only : offline_fw_fluxes_into_ocean, offline_fw_fluxes_out_ocean
use MOM_offline_main,          only : offline_advection_layer, offline_transport_end
use MOM_ice_shelf,             only : ice_shelf_CS, ice_shelf_query, initialize_ice_shelf
use MOM_particles_mod,         only : particles, particles_init, particles_run, particles_save_restart, particles_end
use MOM_particles_mod,         only : particles_to_k_space, particles_to_z_space
implicit none ; private

#include <MOM_memory.h>

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure with diagnostic IDs of the state variables
type MOM_diag_IDs
  !>@{ 3-d state field diagnostic IDs
  integer :: id_u  = -1, id_v  = -1, id_h  = -1
  !>@}
  !> 2-d state field diagnostic ID
  integer :: id_ssh_inst = -1
end type MOM_diag_IDs

!> Control structure for the MOM module, including the variables that describe
!! the state of the ocean.
type, public :: MOM_control_struct ; private
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
    h, &            !< layer thickness [H ~> m or kg m-2]
    T, &            !< potential temperature [C ~> degC]
    S               !< salinity [S ~> ppt]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    u,  &           !< zonal velocity component [L T-1 ~> m s-1]
    uh, &           !< uh = u * h * dy at u grid points [H L2 T-1 ~> m3 s-1 or kg s-1]
    uhtr            !< accumulated zonal thickness fluxes to advect tracers [H L2 ~> m3 or kg]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    v,  &           !< meridional velocity [L T-1 ~> m s-1]
    vh, &           !< vh = v * h * dx at v grid points [H L2 T-1 ~> m3 s-1 or kg s-1]
    vhtr            !< accumulated meridional thickness fluxes to advect tracers [H L2 ~> m3 or kg]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: ssh_rint
                    !< A running time integral of the sea surface height [T Z ~> s m].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: ave_ssh_ibc
                    !< time-averaged (over a forcing time step) sea surface height
                    !! with a correction for the inverse barometer [Z ~> m]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: eta_av_bc
                    !< free surface height or column mass time averaged over the last
                    !! baroclinic dynamics time step [H ~> m or kg m-2]
  real, dimension(:,:), pointer :: Hml => NULL()
                    !< active mixed layer depth, or 0 if there is no boundary layer scheme [Z ~> m]
  real :: time_in_cycle !< The running time of the current time-stepping cycle
                    !! in calls that step the dynamics, and also the length of
                    !! the time integral of ssh_rint [T ~> s].
  real :: time_in_thermo_cycle !< The running time of the current time-stepping
                    !! cycle in calls that step the thermodynamics [T ~> s].

  type(ocean_grid_type) :: G_in                   !< Input grid metric
  type(ocean_grid_type), pointer :: G => NULL()   !< Model grid metric
  logical :: rotate_index = .false.   !< True if index map is rotated
  logical :: homogenize_forcings = .false. !< True if all inputs are homogenized
  logical :: update_ustar = .false.   !< True to update ustar from homogenized tau
  logical :: vertex_shear = .false. !< True if vertex shear is on

  type(verticalGrid_type), pointer :: &
    GV => NULL()    !< structure containing vertical grid info
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors
  type(thermo_var_ptrs) :: tv !< structure containing pointers to available thermodynamic fields
  real :: t_dyn_rel_adv !< The time of the dynamics relative to tracer advection and lateral mixing
                    !! [T ~> s], or equivalently the elapsed time since advectively updating the
                    !! tracers.  t_dyn_rel_adv is invariably positive and may span multiple coupling timesteps.
  integer :: n_dyn_steps_in_adv !< The number of dynamics time steps that contributed to uhtr
                    !! and vhtr since the last time tracer advection occured.
  real :: t_dyn_rel_thermo  !< The time of the dynamics relative to diabatic  processes and remapping
                    !! [T ~> s].  t_dyn_rel_thermo can be negative or positive depending on whether
                    !! the diabatic processes are applied before or after the dynamics and may span
                    !! multiple coupling timesteps.
  real :: t_dyn_rel_diag !< The time of the diagnostics relative to diabatic processes and remapping
                    !!  [T ~> s].  t_dyn_rel_diag is always positive, since the diagnostics must lag.
  logical :: preadv_h_stored = .false. !< If true, the thicknesses from before the advective cycle
                    !! have been stored for use in diagnostics.

  type(diag_ctrl)     :: diag !< structure to regulate diagnostic output timing
  type(vertvisc_type) :: visc !< structure containing vertical viscosities,
                    !! bottom drag viscosities, and related fields
  type(MEKE_type) :: MEKE   !< Fields related to the Mesoscale Eddy Kinetic Energy
  logical :: adiabatic !< If true, there are no diapycnal mass fluxes, and no calls
                    !! to routines to calculate or apply diapycnal fluxes.
  logical :: diabatic_first !< If true, apply diabatic and thermodynamic processes before time
                    !! stepping the dynamics.
  logical :: use_ALE_algorithm  !< If true, use the ALE algorithm rather than layered
                    !! isopycnal/stacked shallow water mode. This logical is set by calling the
                    !! function useRegridding() from the MOM_regridding module.
  logical :: remap_aux_vars     !< If true, apply ALE remapping to all of the auxiliary 3-D
                    !! variables that are needed to reproduce across restarts,
                    !! similarly to what is done with the primary state variables.
  logical :: remap_uv_using_old_alg !< If true, use the old "remapping via a delta z" method for
                    !! velocities.  If false, remap between two grids described by thicknesses.

  type(MOM_stoch_eos_CS) :: stoch_eos_CS !< structure containing random pattern for stoch EOS
  logical :: alternate_first_direction !< If true, alternate whether the x- or y-direction
                    !! updates occur first in directionally split parts of the calculation.
  real    :: first_dir_restart = -1.0 !< A real copy of G%first_direction for use in restart files [nondim]
  logical :: offline_tracer_mode = .false.
                    !< If true, step_offline() is called instead of step_MOM().
                    !! This is intended for running MOM6 in offline tracer mode
  logical :: MEKE_in_dynamics !< If .true. (default), MEKE is called in the dynamics routine otherwise
                              !! it is called during the tracer dynamics

  type(time_type), pointer :: Time   !< pointer to the ocean clock
  real    :: dt                      !< (baroclinic) dynamics time step [T ~> s]
  real    :: dt_therm                !< diabatic time step [T ~> s]
  real    :: dt_tr_adv               !< tracer advection time step [T ~> s]
  logical :: thermo_spans_coupling   !< If true, thermodynamic and tracer time
                                     !! steps can span multiple coupled time steps.
  logical :: tradv_spans_coupling    !< If true, thermodynamic and tracer time
  integer :: nstep_tot = 0           !< The total number of dynamic timesteps taken
                                     !! so far in this run segment
  logical :: count_calls = .false.   !< If true, count the calls to step_MOM, rather than the
                                     !! number of dynamics steps in nstep_tot
  logical :: debug                   !< If true, write verbose checksums for debugging purposes.
  logical :: debug_OBCs              !< If true, write verbose OBC values for debugging purposes.
  integer :: ntrunc                  !< number u,v truncations since last call to write_energy

  integer :: cont_stencil            !< The stencil for thickness from the continuity solver.
  integer :: dyn_h_stencil           !< The stencil for thickness for the dynamics based on
                                     !! the continuity solver and Coriolis schemes.
  ! These elements are used to control the dynamics updates.
  logical :: do_dynamics             !< If false, does not call step_MOM_dyn_*. This is an
                                     !! undocumented run-time flag that is fragile.
  logical :: split                   !< If true, use the split time stepping scheme.
  logical :: use_alt_split           !< If true, use a version of the split explicit time stepping
                                     !! scheme that exchanges velocities with step_MOM that have the
                                     !! average barotropic phase over a baroclinic timestep rather
                                     !! than the instantaneous barotropic phase.
  logical :: use_RK2                 !< If true, use RK2 instead of RK3 in unsplit mode
                                     !! (i.e., no split between barotropic and baroclinic).
  logical :: interface_filter        !< If true, apply an interface height filter immediately
                                     !! after any calls to thickness_diffuse.
  logical :: thickness_diffuse       !< If true, diffuse interface height w/ a diffusivity KHTH.
  logical :: thickness_diffuse_first !< If true, diffuse thickness before dynamics.
  logical :: interface_filter_dt_bug !< If true, uses the wrong time interval in
                                     !! calls to interface_filter and thickness_diffuse.
  logical :: mixedlayer_restrat      !< If true, use submesoscale mixed layer restratifying scheme.
  logical :: useMEKE                 !< If true, call the MEKE parameterization.
  logical :: use_stochastic_EOS      !< If true, use the stochastic EOS parameterizations.
  logical :: useWaves                !< If true, update Stokes drift
  real :: dtbt_reset_period          !< The time interval between dynamic recalculation of the
                                     !! barotropic time step [T ~> s]. If this is negative dtbt is never
                                     !! calculated, and if it is 0, dtbt is calculated every step.
  type(time_type) :: dtbt_reset_interval !< A time_time representation of dtbt_reset_period.
  type(time_type) :: dtbt_reset_time     !< The next time DTBT should be calculated.
  real            :: dt_obc_seg_period   !< The time interval between OBC segment updates for OBGC
                                         !! tracers [T ~> s], or a negative value if the segment
                                         !! data are time-invarant, or zero to update the OBGC
                                         !! segment data with every call to update_OBC_segment_data.
  type(time_type) :: dt_obc_seg_interval !< A time_time representation of dt_obc_seg_period.
  type(time_type) :: dt_obc_seg_time     !< The next time OBC segment update is applied to OBGC tracers.

  real, dimension(:,:), pointer :: frac_shelf_h => NULL() !< fraction of total area occupied
  !! by ice shelf [nondim]
  real, dimension(:,:), pointer :: mass_shelf => NULL() !< Mass of ice shelf [R Z ~> kg m-2]
  type(accel_diag_ptrs) :: ADp  !< structure containing pointers to accelerations,
                                !! for derived diagnostics (e.g., energy budgets)
  type(cont_diag_ptrs)  :: CDp  !< structure containing pointers to continuity equation
                                !! terms, for derived diagnostics (e.g., energy budgets)
  real, dimension(:,:,:), pointer :: &
    u_prev => NULL(), &         !< previous value of u stored for diagnostics [L T-1 ~> m s-1]
    v_prev => NULL()            !< previous value of v stored for diagnostics [L T-1 ~> m s-1]

  logical :: interp_p_surf      !< If true, linearly interpolate surface pressure
                                !! over the coupling time step, using specified value
                                !! at the end of the coupling step. False by default.
  logical :: p_surf_prev_set    !< If true, p_surf_prev has been properly set from
                                !! a previous time-step or the ocean restart file.
                                !! This is only valid when interp_p_surf is true.
  real, dimension(:,:), pointer :: &
    p_surf_prev  => NULL(), &   !< surface pressure [R L2 T-2 ~> Pa] at end  previous call to step_MOM
    p_surf_begin => NULL(), &   !< surface pressure [R L2 T-2 ~> Pa] at start of step_MOM_dyn_...
    p_surf_end   => NULL()      !< surface pressure [R L2 T-2 ~> Pa] at end   of step_MOM_dyn_...

  ! Variables needed to reach between start and finish phases of initialization
  logical :: write_IC           !< If true, then the initial conditions will be written to file
  character(len=120) :: IC_file !< A file into which the initial conditions are
                                !! written in a new run if SAVE_INITIAL_CONDS is true.

  logical :: calc_rho_for_sea_lev !< If true, calculate rho to convert pressure to sea level

  ! These elements are used to control the calculation and error checking of the surface state
  real :: Hmix                  !< Diagnostic mixed layer thickness over which to
                                !! average surface tracer properties when a bulk
                                !! mixed layer is not used [H ~> m or kg m-2], or a negative value
                                !! if a bulk mixed layer is being used.
  real :: HFrz                  !< If HFrz > 0, the nominal depth over which melt potential is computed
                                !! [H ~> m or kg m-2].  The actual depth over which melt potential is
                                !! computed is min(HFrz, OBLD), where OBLD is the boundary layer depth.
                                !! If HFrz <= 0 (default), melt potential will not be computed.
  real :: Hmix_UV               !< Depth scale over which to average surface flow to
                                !! feedback to the coupler/driver [H ~> m or kg m-2] when
                                !! bulk mixed layer is not used, or a negative value
                                !! if a bulk mixed layer is being used.
  logical :: check_bad_sfc_vals !< If true, scan surface state for ridiculous values.
  real    :: bad_val_ssh_max    !< Maximum SSH before triggering bad value message [Z ~> m]
  real    :: bad_val_sst_max    !< Maximum SST before triggering bad value message [C ~> degC]
  real    :: bad_val_sst_min    !< Minimum SST before triggering bad value message [C ~> degC]
  real    :: bad_val_sss_max    !< Maximum SSS before triggering bad value message [S ~> ppt]
  real    :: bad_val_col_thick  !< Minimum column thickness before triggering bad value message [Z ~> m]
  integer :: answer_date        !< The vintage of the expressions for the surface properties.  Values
                                !! below 20190101 recover the answers from the end of 2018, while
                                !! higher values use more appropriate expressions that differ at
                                !! roundoff for non-Boussinesq cases.
  logical :: use_particles      !< Turns on the particles package
  logical :: use_uh_particles   !< particles are advected by uh/h
  logical :: uh_particles_bug   !< If true, uses an inconsistent timestep for particle advection
  logical :: use_dbclient       !< Turns on the database client used for ML inference/analysis
  character(len=10) :: particle_type !< Particle types include: surface(default), profiling and sail drone.

  type(MOM_diag_IDs)       :: IDs      !<  Handles used for diagnostics.
  type(transport_diag_IDs) :: transport_IDs  !< Handles used for transport diagnostics.
  type(surface_diag_IDs)   :: sfc_IDs  !< Handles used for surface diagnostics.
  type(diag_grid_storage)  :: diag_pre_sync !< The grid (thicknesses) before remapping
  type(diag_grid_storage)  :: diag_pre_dyn  !< The grid (thicknesses) before dynamics

  ! The remainder of this type provides pointers to child module control structures.

  type(MOM_dyn_unsplit_CS),      pointer :: dyn_unsplit_CSp => NULL()
    !< Pointer to the control structure used for the unsplit dynamics
  type(MOM_dyn_unsplit_RK2_CS),  pointer :: dyn_unsplit_RK2_CSp => NULL()
    !< Pointer to the control structure used for the unsplit RK2 dynamics
  type(MOM_dyn_split_RK2_CS),    pointer :: dyn_split_RK2_CSp => NULL()
    !< Pointer to the control structure used for the mode-split RK2 dynamics
  type(MOM_dyn_split_RK2b_CS),    pointer :: dyn_split_RK2b_CSp => NULL()
    !< Pointer to the control structure used for an alternate version of the mode-split RK2 dynamics
  type(harmonic_analysis_CS),    pointer :: HA_CSp => NULL()
    !< Pointer to the control structure for harmonic analysis
  type(thickness_diffuse_CS) :: thickness_diffuse_CSp
    !< Pointer to the control structure used for the isopycnal height diffusive transport.
    !! This is also common referred to as Gent-McWilliams diffusion
  type(interface_filter_CS) :: interface_filter_CSp
    !< Control structure used for the interface height smoothing operator.
  type(mixedlayer_restrat_CS) :: mixedlayer_restrat_CSp
    !< Pointer to the control structure used for the mixed layer restratification
  type(set_visc_CS)           :: set_visc_CSp
    !< Pointer to the control structure used to set viscosities
  type(diabatic_CS),             pointer :: diabatic_CSp => NULL()
    !< Pointer to the control structure for the diabatic driver
  type(MEKE_CS) :: MEKE_CSp
    !< Pointer to the control structure for the MEKE updates
  type(VarMix_CS) :: VarMix
    !< Control structure for the variable mixing module
  type(tracer_registry_type),    pointer :: tracer_Reg => NULL()
    !< Pointer to the MOM tracer registry
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp => NULL()
    !< Pointer to the MOM tracer advection control structure
  type(tracer_hor_diff_CS),      pointer :: tracer_diff_CSp => NULL()
    !< Pointer to the MOM along-isopycnal tracer diffusion control structure
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp => NULL()
    !< Pointer to the control structure that orchestrates the calling of tracer packages
    ! Although update_OBC_CS is not used directly outside of initialization, other modules
    ! set pointers to this type, so it should be kept for the duration of the run.
  type(update_OBC_CS),           pointer :: update_OBC_CSp => NULL()
    !< Pointer to the control structure for updating open boundary condition properties
  type(ocean_OBC_type),          pointer :: OBC => NULL()
    !< Pointer to the MOM open boundary condition type
  type(sponge_CS),               pointer :: sponge_CSp => NULL()
    !< Pointer to the layered-mode sponge control structure
  type(ALE_sponge_CS),           pointer :: ALE_sponge_CSp => NULL()
    !< Pointer to the ALE-mode sponge control structure
  type(oda_incupd_CS),           pointer :: oda_incupd_CSp => NULL()
    !< Pointer to the oda incremental update control structure
  type(int_tide_CS),             pointer :: int_tide_CSp => NULL()
    !< Pointer to the internal tides control structure
  type(ALE_CS),                  pointer :: ALE_CSp => NULL()
    !< Pointer to the Arbitrary Lagrangian Eulerian (ALE) vertical coordinate control structure

  ! Pointers to control structures used for diagnostics
  type(sum_output_CS),           pointer :: sum_output_CSp => NULL()
    !< Pointer to the globally summed output control structure
  type(diagnostics_CS) :: diagnostics_CSp
    !< Pointer to the MOM diagnostics control structure
  type(offline_transport_CS),    pointer :: offline_CSp => NULL()
    !< Pointer to the offline tracer transport control structure
  type(porous_barrier_CS)                :: por_bar_CS
    !< Control structure for porous barrier

  logical               :: ensemble_ocean !< if true, this run is part of a
                                !! larger ensemble for the purpose of data assimilation
                                !! or statistical analysis.
  type(ODA_CS), pointer :: odaCS => NULL() !< a pointer to the control structure for handling
                                !! ensemble model state vectors and data assimilation
                                !! increments and priors
  type(dbcomms_CS_type)   :: dbcomms_CS !< Control structure for database client used for online ML/AI
  logical :: use_porbar !< If true, use porous barrier to constrain the widths and face areas
                        !! at the edges of the grid cells.
  type(porous_barrier_type) :: pbv !< porous barrier fractional cell metrics
  type(particles), pointer :: particles => NULL() !<Lagrangian particles
  type(stochastic_CS), pointer :: stoch_CS => NULL() !< a pointer to the stochastics control structure
  type(MOM_restart_CS), pointer :: restart_CS => NULL()
    !< Pointer to MOM's restart control structure
end type MOM_control_struct

public initialize_MOM, finish_MOM_initialization, MOM_end
public step_MOM, step_offline
public extract_surface_state, get_ocean_stocks
public get_MOM_state_elements, MOM_state_is_synchronized
public allocate_surface_state, deallocate_surface_state
public save_MOM_restart

!>@{ CPU time clock IDs
integer :: id_clock_ocean
integer :: id_clock_dynamics
integer :: id_clock_thermo
integer :: id_clock_MOM_end
integer :: id_clock_remap
integer :: id_clock_tracer
integer :: id_clock_diabatic
integer :: id_clock_adiabatic
integer :: id_clock_continuity  ! also in dynamics s/r
integer :: id_clock_thick_diff
integer :: id_clock_int_filter
integer :: id_clock_BBL_visc
integer :: id_clock_ml_restrat
integer :: id_clock_diagnostics
integer :: id_clock_Z_diag
integer :: id_clock_init
integer :: id_clock_MOM_init
integer :: id_clock_pass       ! also in dynamics d/r
integer :: id_clock_pass_init  ! also in dynamics d/r
integer :: id_clock_ALE
integer :: id_clock_other
integer :: id_clock_offline_tracer
integer :: id_clock_save_restart
integer :: id_clock_unit_tests
integer :: id_clock_stoch
integer :: id_clock_varT
!>@}


  interface
module subroutine step_MOM(forces_in, fluxes_in, sfc_state, Time_start, time_int_in, CS, &
                    Waves, do_dynamics, do_thermodynamics, start_cycle, &
                    end_cycle, cycle_length, reset_therm)
  type(mech_forcing), target, intent(inout) :: forces_in !< A structure with the driving mechanical forces
  type(forcing), target, intent(inout) :: fluxes_in  !< A structure with pointers to themodynamic,
                                                     !! tracer and mass exchange forcing fields
  type(surface), target, intent(inout) :: sfc_state  !< surface ocean state
  type(time_type),    intent(in)    :: Time_start    !< starting time of a segment, as a time type
  real,               intent(in)    :: time_int_in   !< time interval covered by this run segment [T ~> s].
  type(MOM_control_struct), intent(inout), target :: CS   !< control structure from initialize_MOM
  type(Wave_parameters_CS), &
            optional, pointer       :: Waves         !< An optional pointer to a wave property CS
  logical,  optional, intent(in)    :: do_dynamics   !< Present and false, do not do updates due
                                                     !! to the dynamics.
  logical,  optional, intent(in)    :: do_thermodynamics  !< Present and false, do not do updates due
                                                     !! to the thermodynamics or remapping.
  logical,  optional, intent(in)    :: start_cycle   !< This indicates whether this call is to be
                                                     !! treated as the first call to step_MOM in a
                                                     !! time-stepping cycle; missing is like true.
  logical,  optional, intent(in)    :: end_cycle     !< This indicates whether this call is to be
                                                     !! treated as the last call to step_MOM in a
                                                     !! time-stepping cycle; missing is like true.
  real,     optional, intent(in)    :: cycle_length  !< The amount of time in a coupled time
                                                     !! stepping cycle [T ~> s].
  logical,  optional, intent(in)    :: reset_therm   !< This indicates whether the running sums of
                                                     !! thermodynamic quantities should be reset.
                                                     !! If missing, this is like start_cycle.

  ! local variables
                                                   ! metrics and related information
                                                   ! various unit conversion factors



                          ! and beginning of the current time step [nondim]
                          ! properties will apply, for use in diagnostics, or 0
                          ! if it is not to be calculated anew [T ~> s].

                                       ! multiple coupling timesteps.
                                       ! multiple coupling timesteps.
                        ! can use nonblocking halo updates
                        ! a stepping cycle (whatever that may mean).
                        ! the end of a stepping cycle (whatever that may mean).
                ! the time-evolving surface density in non-Boussinesq mode [Z T-1 ~> m s-1]




  ! External forcing fields on the model index map

end subroutine step_MOM
module subroutine step_MOM_dynamics(forces, p_surf_begin, p_surf_end, dt, dt_tr_adv, &
                             bbl_time_int, CS, Time_local, Waves)
  type(mech_forcing), intent(in)    :: forces     !< A structure with the driving mechanical forces
  real, dimension(:,:), pointer     :: p_surf_begin !< A pointer (perhaps NULL) to the surface
                                                  !! pressure at the beginning of this dynamic
                                                  !! step, intent in [R L2 T-2 ~> Pa].
  real, dimension(:,:), pointer     :: p_surf_end !< A pointer (perhaps NULL) to the surface
                                                  !! pressure at the end of this dynamic step,
                                                  !! intent in [R L2 T-2 ~> Pa].
  real,               intent(in)    :: dt         !< time interval covered by this call [T ~> s].
  real,               intent(in)    :: dt_tr_adv  !< time interval covered by any updates that may
                                                  !! span multiple dynamics steps [T ~> s].
  real,               intent(in)    :: bbl_time_int !< time interval over which updates to the
                                                  !! bottom boundary layer properties will apply [T ~> s],
                                                  !! or zero not to update the properties.
  type(MOM_control_struct), intent(inout), target :: CS   !< control structure from initialize_MOM
  type(time_type),    intent(in)    :: Time_local !< End time of a segment, as a time type
  type(wave_parameters_CS), &
            optional, pointer       :: Waves      !< Container for wave related parameters; the
                                                  !! fields in Waves are intent in here.

  ! local variables
                                                   ! metrics and related information
                                                   ! various unit conversion factors

                        ! barotropic time step needs to be updated.


end subroutine step_MOM_dynamics
module subroutine step_MOM_tracer_dyn(CS, G, GV, US, h, Time_local)
  type(MOM_control_struct), intent(inout) :: CS     !< control structure
  type(ocean_grid_type),    intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h      !< layer thicknesses after the transports [H ~> m or kg m-2]
  type(time_type),          intent(in)    :: Time_local !< The model time at the end
                                                    !! of the time step.
end subroutine step_MOM_tracer_dyn
module subroutine step_MOM_thermo(CS, G, GV, US, u, v, h, tv, fluxes, dtdia, &
                           Time_end_thermo, update_BBL, Waves)
  type(MOM_control_struct), intent(inout) :: CS     !< Master MOM control structure
  type(ocean_grid_type),    intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),  intent(inout) :: GV     !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: u      !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(inout) :: v      !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: h      !< layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables
  type(forcing),            intent(inout) :: fluxes !< pointers to forcing fields
  real,                     intent(in)    :: dtdia  !< The time interval over which to advance [T ~> s]
  type(time_type),          intent(in)    :: Time_end_thermo !< End of averaging interval for thermo diags
  logical,                  intent(in)    :: update_BBL !< If true, calculate the bottom boundary layer properties.
  type(wave_parameters_CS), &
                  optional, pointer       :: Waves  !< Container for wave related parameters
                                                    !! the fields in Waves are intent in here.

                               ! in the dynamic core.

end subroutine step_MOM_thermo
module subroutine ALE_regridding_and_remapping(CS, G, GV, US, u, v, h, tv, dtdia, Time_end_thermo)
  type(MOM_control_struct), intent(inout) :: CS     !< Master MOM control structure
  type(ocean_grid_type),    intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),  intent(inout) :: GV     !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: u      !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(inout) :: v      !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: h      !< layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables
  real,                     intent(in)    :: dtdia  !< The time interval over which to advance [T ~> s]
  type(time_type),          intent(in)    :: Time_end_thermo !< End of averaging interval for thermo diags

                                               ! in the same units as thicknesses [H ~> m or kg m-2]
                                               ! velocity points [H ~> m or kg m-2]
                                               ! velocity points [H ~> m or kg m-2]
                                               ! velocity points [H ~> m or kg m-2]
                                               ! velocity points [H ~> m or kg m-2]

end subroutine ALE_regridding_and_remapping
module subroutine post_diabatic_halo_updates(CS, G, GV, US, u, v, h, tv)
  type(MOM_control_struct), intent(inout) :: CS     !< Master MOM control structure
  type(ocean_grid_type),    intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),  intent(inout) :: GV     !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: u      !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(inout) :: v      !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: h      !< layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables

                               ! in the dynamic core.

end subroutine post_diabatic_halo_updates
module subroutine step_offline(forces, fluxes, sfc_state, Time_start, time_interval, CS)
  type(mech_forcing), intent(in)    :: forces        !< A structure with the driving mechanical forces
  type(forcing),      intent(inout) :: fluxes        !< pointers to forcing fields
  type(surface),      intent(inout) :: sfc_state     !< surface ocean state
  type(time_type),    intent(in)    :: Time_start    !< starting time of a segment, as a time type
  real,               intent(in)    :: time_interval !< time interval [T ~> s]
  type(MOM_control_struct), intent(inout) :: CS      !< control structure from initialize_MOM

  ! Local pointers
                                                      ! metrics and related information
                                                      ! about the vertical grid
                                                      ! various unit conversion factors


                                                  ! in the same units as thicknesses [H ~> m or kg m-2]

                               ! in the dynamic core.

  ! 3D pointers


  ! Grid-related pointer assignments
end subroutine step_offline
module subroutine initialize_MOM(Time, Time_init, param_file, dirs, CS, &
                          Time_in, offline_tracer_mode, input_restart_file, diag_ptr, &
                          count_calls, tracer_flow_CSp,  ice_shelf_CSp, waves_CSp, ensemble_num, &
                          calve_ice_shelf_bergs)
  type(time_type), target,   intent(inout) :: Time        !< model time, set in this routine
  type(time_type),           intent(in)    :: Time_init   !< The start time for the coupled model's calendar
  type(param_file_type),     intent(out)   :: param_file  !< structure indicating parameter file to parse
  type(directories),         intent(out)   :: dirs        !< structure with directory paths
  type(MOM_control_struct),  intent(inout), target :: CS  !< pointer set in this routine to MOM control structure
  type(time_type), optional, intent(in)    :: Time_in     !< time passed to MOM_initialize_state when
                                                          !! model is not being started from a restart file
  logical,         optional, intent(out)   :: offline_tracer_mode !< True is returned if tracers are being run offline
  character(len=*),optional, intent(in)    :: input_restart_file !< If present, name of restart file to read
  type(diag_ctrl), optional, pointer       :: diag_ptr    !< A pointer set in this routine to the diagnostic
                                                          !! regulatory structure
  type(tracer_flow_control_CS), &
                   optional, pointer       :: tracer_flow_CSp !< A pointer set in this routine to
                                                          !! the tracer flow control structure.
  logical,         optional, intent(in)    :: count_calls !< If true, nstep_tot counts the number of
                                                          !! calls to step_MOM instead of the number of
                                                          !! dynamics timesteps.
  type(ice_shelf_CS), optional,     pointer :: ice_shelf_CSp !< A pointer to an ice shelf control structure
  type(Wave_parameters_CS), &
                   optional, pointer       :: Waves_CSp   !< An optional pointer to a wave property CS
  integer, optional :: ensemble_num                       !< Ensemble index provided by the cap (instead of FMS
                                                          !! ensemble manager)
  logical, optional :: calve_ice_shelf_bergs !< If true, will add point iceberg calving variables to the ice
                                             !! shelf restart
  ! local variables

  ! Initial state on the input index map
                                                  ! by an ice shelf [nondim]
                                                  ! [R Z ~> kg m-2]

  ! This include declares and sets the variable "version".

                               ! of the maximum stable value [nondim].

                                                  ! in the same units as thicknesses [H ~> m or kg m-2]
                                                  ! points [H ~> m or kg m-2]
                                                  ! velocity points [H ~> m or kg m-2]
                                                  ! velocity points [H ~> m or kg m-2]


                               ! with nkml sublayers and nkbl buffer layer.
                               ! in equation of state calculations.
                               ! with accumulated heat deficit returned to surface ocean.
                               ! a minimum value, and the deficit is reported.
                               ! and absolute salinity. Care should be taken to convert them
                               ! to potential temperature and practical salinity before
                               ! exchanging them with the coupler and/or reporting T&S diagnostics.
                               ! and salnity is performed
                               ! of having the data domain on each processor start at 1.
                               ! the velocity points.
                               ! time step needs to be updated before it is used.
                               ! updated first in directionally split parts of the
                               ! calculation.
                               ! set to recreate the bugs so that the code can be moved forward
                               ! without changing answers for existing configurations.  When this is
                               ! false, bugs are only used if they are actively selected.
                               ! in the dynamic core.
                               ! fluxes [J m-2 H-1 C-1 ~> J m-3 degC-1 or J kg-1 degC-1]

                                                                ! (To be used for writing out ocean geometry)

end subroutine initialize_MOM
module subroutine finish_MOM_initialization(Time, dirs, CS)
  type(time_type),          intent(in)    :: Time        !< model time, used in this routine
  type(directories),        intent(in)    :: dirs        !< structure with directory paths
  type(MOM_control_struct), intent(inout) :: CS          !< MOM control structure

                                                   ! metrics and related information
                                                   ! various unit conversion factors

end subroutine finish_MOM_initialization
module subroutine register_diags(Time, G, GV, US, IDs, diag)
  type(time_type),         intent(in)    :: Time  !< current model time
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(inout) :: US    !< A dimensional unit scaling type
  type(MOM_diag_IDs),      intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),         intent(inout) :: diag  !< regulates diagnostic output


end subroutine register_diags
module subroutine MOM_timing_init(CS)
  type(MOM_control_struct), intent(in) :: CS  !< control structure set up by initialize_MOM.

end subroutine MOM_timing_init
module subroutine set_restart_fields(GV, US, param_file, CS, restart_CSp)
  type(verticalGrid_type),  intent(inout) :: GV         !< ocean vertical grid structure
  type(unit_scale_type),    intent(inout) :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in) :: param_file    !< opened file for parsing to get parameters
  type(MOM_control_struct), intent(in) :: CS            !< control structure set up by initialize_MOM
  type(MOM_restart_CS),     pointer    :: restart_CSp   !< pointer to the restart control
                                                        !! structure that will be used for MOM.
  ! Local variables

end subroutine set_restart_fields
module subroutine adjust_ssh_for_p_atm(tv, G, GV, US, ssh, p_atm, use_EOS)
  type(thermo_var_ptrs),             intent(in)    :: tv  !< A structure pointing to various thermodynamic variables
  type(ocean_grid_type),             intent(in)    :: G   !< ocean grid structure
  type(verticalGrid_type),           intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),             intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: ssh !< time mean surface height [Z ~> m]
  real, dimension(:,:),              pointer       :: p_atm !< Ocean surface pressure [R L2 T-2 ~> Pa]
  logical,                           intent(in)    :: use_EOS !< If true, calculate the density for
                                                       !! the SSH correction using the equation of state.

                      ! a corrected effective SSH [R ~> kg m-3].

end subroutine adjust_ssh_for_p_atm
module subroutine extract_surface_state(CS, sfc_state_in)
  type(MOM_control_struct), intent(inout), target :: CS   !< Master MOM control structure
  type(surface), target, intent(inout) :: sfc_state_in !< transparent ocean surface state
                                             !! structure shared with the calling routine
                                             !! data in this structure is intent out.

  ! Local variables
                                                  !! metrics and related information
                             !! layer properties [Z ~> m] or [H ~> m or kg m-2]
                             !! calculation of properties of the uppermost ocean [nondim] or [Z H-1 ~> 1 or m3 kg-1]
                             !  After the ANSWERS_2018 flag has been obsoleted, H_rescale will be 1.

end subroutine extract_surface_state
module subroutine rotate_initial_state(u_in, v_in, h_in, T_in, S_in, &
    use_temperature, turns, u, v, h, T, S)
  real, dimension(:,:,:), intent(in)  :: u_in  !< Zonal velocity on the initial grid [L T-1 ~> m s-1]
  real, dimension(:,:,:), intent(in)  :: v_in  !< Meridional velocity on the initial grid [L T-1 ~> m s-1]
  real, dimension(:,:,:), intent(in)  :: h_in  !< Layer thickness on the initial grid [H ~> m or kg m-2]
  real, dimension(:,:,:), intent(in)  :: T_in  !< Temperature on the initial grid [C ~> degC]
  real, dimension(:,:,:), intent(in)  :: S_in  !< Salinity on the initial grid [S ~> ppt]
  logical,                intent(in)  :: use_temperature !< If true, temperature and salinity are active
  integer,                intent(in)  :: turns !< The number quarter-turns to apply
  real, dimension(:,:,:), intent(out) :: u     !< Zonal velocity on the rotated grid [L T-1 ~> m s-1]
  real, dimension(:,:,:), intent(out) :: v     !< Meridional velocity on the rotated grid [L T-1 ~> m s-1]
  real, dimension(:,:,:), intent(out) :: h     !< Layer thickness on the rotated grid [H ~> m or kg m-2]
  real, dimension(:,:,:), intent(out) :: T     !< Temperature on the rotated grid [C ~> degC]
  real, dimension(:,:,:), intent(out) :: S     !< Salinity on the rotated grid [S ~> ppt]

end subroutine rotate_initial_state
module function MOM_state_is_synchronized(CS, adv_dyn) result(in_synch)
  type(MOM_control_struct), intent(inout) :: CS !< MOM control structure
  logical,        optional, intent(in) :: adv_dyn  !< If present and true, only check
                                          !! whether the advection is up-to-date with
                                          !! the dynamics.
  logical :: in_synch !< True if all phases of the update are synchronized.


end function MOM_state_is_synchronized
module subroutine get_MOM_state_elements(CS, G, GV, US, C_p, C_p_scaled, use_temp)
  type(MOM_control_struct), intent(inout), target :: CS  !< MOM control structure
  type(ocean_grid_type),   optional, pointer     :: G    !< structure containing metrics and grid info
  type(verticalGrid_type), optional, pointer     :: GV   !< structure containing vertical grid info
  type(unit_scale_type),   optional, pointer     :: US   !< A dimensional unit scaling type
  real,                    optional, intent(out) :: C_p  !< The heat capacity [J kg degC-1]
  real,                    optional, intent(out) :: C_p_scaled !< The heat capacity in scaled
                                                         !! units [Q C-1 ~> J kg-1 degC-1]
  logical,                 optional, intent(out) :: use_temp !< True if temperature is a state variable

end subroutine get_MOM_state_elements
module subroutine get_ocean_stocks(CS, mass, heat, salt, on_PE_only)
  type(MOM_control_struct), intent(inout) :: CS !< MOM control structure
  real,    optional, intent(out) :: heat  !< The globally integrated integrated ocean heat [J].
  real,    optional, intent(out) :: salt  !< The globally integrated integrated ocean salt [kg].
  real,    optional, intent(out) :: mass  !< The globally integrated integrated ocean mass [kg].
  logical, optional, intent(in)  :: on_PE_only !< If present and true, only sum on the local PE.

end subroutine get_ocean_stocks
module subroutine save_MOM_restart(CS, directory, time, G, time_stamped, filename, &
    GV, num_rest_files, write_IC)
  type(MOM_control_struct), intent(inout) :: CS
    !< MOM control structure
  character(len=*), intent(in) :: directory
    !< The directory where the restart files are to be written
  type(time_type), intent(in) :: time
    !< The current model time
  type(ocean_grid_type), intent(inout) :: G
    !< The ocean's grid structure
  logical, optional, intent(in) :: time_stamped
    !< If present and true, add time-stamp to the restart file names
  character(len=*), optional, intent(in) :: filename
    !< A filename that overrides the name in CS%restartfile
  type(verticalGrid_type), optional, intent(in) :: GV
    !< The ocean's vertical grid structure
  integer, optional, intent(out) :: num_rest_files
    !< number of restart files written
  logical, optional, intent(in) :: write_IC
    !< If present and true, initial conditions are being written

end subroutine save_MOM_restart
module subroutine MOM_end(CS)
  type(MOM_control_struct), intent(inout) :: CS   !< MOM control structure

end subroutine MOM_end
  end interface

end module MOM
