! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Time step the adiabatic dynamic core of MOM using RK2 method with greater use of the
!! time-filtered velocities and less inheritance of tedencies from the previous step in the
!! predictor step than in the original MOM_dyanmics_split_RK2.
module MOM_dynamics_split_RK2b

use MOM_variables,    only : vertvisc_type, thermo_var_ptrs, porous_barrier_type
use MOM_variables,    only : BT_cont_type, alloc_bt_cont_type, dealloc_bt_cont_type
use MOM_variables,    only : accel_diag_ptrs, ocean_internal_state, cont_diag_ptrs
use MOM_forcing_type, only : mech_forcing

use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,         only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,         only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,     only : diag_mediator_init, enable_averages
use MOM_diag_mediator,     only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator,     only : post_product_u, post_product_sum_u
use MOM_diag_mediator,     only : post_product_v, post_product_sum_v
use MOM_diag_mediator,     only : register_diag_field, register_static_field
use MOM_diag_mediator,     only : set_diag_mediator_grid, diag_ctrl, diag_update_remap_grids
use MOM_domains,           only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_domains,           only : To_North, To_East, Omit_Corners
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : start_group_pass, complete_group_pass, pass_var, pass_vector
use MOM_debugging,         only : hchksum, uvchksum, query_debugging_checks
use MOM_error_handler,     only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,     only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,     only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_get_input,         only : directories
use MOM_io,                only : vardesc, var_desc, EAST_FACE, NORTH_FACE
use MOM_restart,           only : register_restart_field, register_restart_pair
use MOM_restart,           only : query_initialized, set_initialized, save_restart
use MOM_restart,           only : only_read_from_restarts
use MOM_restart,           only : restart_init, is_new_run, MOM_restart_CS
use MOM_time_manager,      only : time_type, operator(+)
use MOM_time_manager,      only : operator(-), operator(>), operator(*), operator(/)

use MOM_ALE,                   only : ALE_CS, ALE_remap_velocities
use MOM_barotropic,            only : barotropic_init, btstep, btcalc, bt_mass_source
use MOM_barotropic,            only : register_barotropic_restarts, set_dtbt, barotropic_CS
use MOM_barotropic,            only : barotropic_end
use MOM_boundary_update,       only : update_OBC_data, update_OBC_CS
use MOM_continuity,            only : continuity, continuity_CS
use MOM_continuity,            only : continuity_init, continuity_stencil
use MOM_CoriolisAdv,           only : CorAdCalc, CoriolisAdv_CS
use MOM_CoriolisAdv,           only : CoriolisAdv_init, CoriolisAdv_end, CoriolisAdv_stencil
use MOM_debugging,             only : check_redundant
use MOM_grid,                  only : ocean_grid_type
use MOM_harmonic_analysis,     only : HA_init, harmonic_analysis_CS
use MOM_hor_index,             only : hor_index_type
use MOM_hor_visc,              only : horizontal_viscosity, hor_visc_CS
use MOM_hor_visc,              only : hor_visc_init, hor_visc_end
use MOM_interface_heights,     only : thickness_to_dz, find_col_avg_SpV
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, radiation_open_bdry_conds
use MOM_open_boundary,         only : open_boundary_zero_normal_flow, open_boundary_query
use MOM_open_boundary,         only : open_boundary_test_extern_h, update_OBC_ramp
use MOM_open_boundary,         only : copy_thickness_reservoirs
use MOM_open_boundary,         only : update_segment_thickness_reservoirs
use MOM_PressureForce,         only : PressureForce, PressureForce_CS
use MOM_PressureForce,         only : PressureForce_init
use MOM_set_visc,              only : set_viscous_ML, set_visc_CS
use MOM_thickness_diffuse,     only : thickness_diffuse_CS
use MOM_self_attr_load,        only : SAL_CS
use MOM_self_attr_load,        only : SAL_init, SAL_end
use MOM_tidal_forcing,         only : tidal_forcing_CS
use MOM_tidal_forcing,         only : tidal_forcing_init, tidal_forcing_end
use MOM_unit_scaling,          only : unit_scale_type
use MOM_vert_friction,         only : vertvisc, vertvisc_coef, vertvisc_remnant
use MOM_vert_friction,         only : vertvisc_init, vertvisc_end, vertvisc_CS
use MOM_vert_friction,         only : updateCFLtruncationValue, vertFPmix
use MOM_verticalGrid,          only : verticalGrid_type, get_thickness_units
use MOM_verticalGrid,          only : get_flux_units, get_tr_flux_units
use MOM_wave_interface,        only : wave_parameters_CS, Stokes_PGF

implicit none ; private

#include <MOM_memory.h>

!> MOM_dynamics_split_RK2b module control structure
type, public :: MOM_dyn_split_RK2b_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    !< CAu = f*v - u.grad(u) [L T-2 ~> m s-2]
    CAu_pred, & !< The predictor step value of CAu = f*v - u.grad(u) [L T-2 ~> m s-2]
    PFu, &    !< PFu = -dM/dx [L T-2 ~> m s-2]
    PFu_Stokes, & !< PFu_Stokes = -d/dx int_r (u_L*duS/dr) [L T-2 ~> m s-2]
    diffu     !< Zonal acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    !< CAv = -f*u - u.grad(v) [L T-2 ~> m s-2]
    CAv_pred, & !< The predictor step value of CAv = -f*u - u.grad(v) [L T-2 ~> m s-2]
    PFv, &    !< PFv = -dM/dy [L T-2 ~> m s-2]
    PFv_Stokes, & !< PFv_Stokes = -d/dy int_r (v_L*dvS/dr) [L T-2 ~> m s-2]
    diffv     !< Meridional acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: visc_rem_u
              !< Both the fraction of the zonal momentum originally in a
              !! layer that remains after a time-step of viscosity, and the
              !! fraction of a time-step worth of a barotropic acceleration
              !! that a layer experiences after viscosity is applied [nondim].
              !! Nondimensional between 0 (at the bottom) and 1 (far above).
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: u_accel_bt
              !< The zonal layer accelerations due to the difference between
              !! the barotropic accelerations and the baroclinic accelerations
              !! that were fed into the barotopic calculation [L T-2 ~> m s-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: visc_rem_v
              !< Both the fraction of the meridional momentum originally in
              !! a layer that remains after a time-step of viscosity, and the
              !! fraction of a time-step worth of a barotropic acceleration
              !! that a layer experiences after viscosity is applied [nondim].
              !! Nondimensional between 0 (at the bottom) and 1 (far above).
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: v_accel_bt
              !< The meridional layer accelerations due to the difference between
              !! the barotropic accelerations and the baroclinic accelerations
              !! that were fed into the barotopic calculation [L T-2 ~> m s-2]

  ! The following variables are only used with the split time stepping scheme.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_)             :: eta    !< Instantaneous free surface height (in Boussinesq
                                                                  !! mode) or column mass anomaly (in non-Boussinesq
                                                                  !! mode) [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: u_av   !< layer x-velocity with vertical mean replaced by
                                                                  !! time-mean barotropic velocity over a baroclinic
                                                                  !! timestep [L T-1 ~> m s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: v_av   !< layer y-velocity with vertical mean replaced by
                                                                  !! time-mean barotropic velocity over a baroclinic
                                                                  !! timestep [L T-1 ~> m s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)      :: h_av   !< arithmetic mean of two successive layer
                                                                  !! thicknesses [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_)             :: eta_PF !< instantaneous SSH used in calculating PFu and
                                                                  !! PFv [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_)        :: uhbt   !< average x-volume or mass flux determined by the
                                                                  !! barotropic solver [H L2 T-1 ~> m3 s-1 or kg s-1].
                                                                  !! uhbt is roughly equal to the vertical sum of uh.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_)        :: vhbt   !< average y-volume or mass flux determined by the
                                                                  !! barotropic solver [H L2 T-1 ~> m3 s-1 or kg s-1].
                                                                  !! vhbt is roughly equal to vertical sum of vh.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)      :: pbce   !< pbce times eta gives the baroclinic pressure
                                                                  !! anomaly in each layer due to free surface height
                                                                  !! anomalies [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_)        :: du_av_inst !< The barotropic zonal velocity increment
                                                                  !! between filtered and instantaneous velocities
                                                                  !! [L T-1 ~> m s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_)        :: dv_av_inst !< The barotropic meridional velocity increment
                                                                  !! between filtered and instantaneous velocities
                                                                  !! [L T-1 ~> m s-1]

  real, pointer, dimension(:,:) :: taux_bot => NULL() !< frictional x-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]
  real, pointer, dimension(:,:) :: tauy_bot => NULL() !< frictional y-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]
  type(BT_cont_type), pointer   :: BT_cont  => NULL() !<  A structure with elements that describe the
                                                      !! effective summed open face areas as a function
                                                      !! of barotropic flow.

  logical :: BT_adj_corr_mass_src !< If true, recalculates the barotropic mass source after
                                  !! predictor step. This should make little difference in the
                                  !! deep ocean but appears to help for vanished layers.
  logical :: split_bottom_stress  !< If true, provide the bottom stress
                                  !! calculated by the vertical viscosity to the
                                  !! barotropic solver.
  logical :: dtbt_use_bt_cont     !< If true, use BT_cont to calculate DTBT.
  logical :: calculate_SAL        !< If true, calculate self-attraction and loading.
  logical :: use_tides            !< If true, tidal forcing is enabled.
  logical :: use_HA               !< If true, perform inline harmonic analysis.
  logical :: remap_aux            !< If true, apply ALE remapping to all of the auxiliary 3-D
                                  !! variables that are needed to reproduce across restarts,
                                  !! similarly to what is done with the primary state variables.

  real    :: be      !< A nondimensional number from 0.5 to 1 that controls
                     !! the backward weighting of the time stepping scheme [nondim]
  real    :: begw    !< A nondimensional number from 0 to 1 that controls
                     !! the extent to which the treatment of gravity waves
                     !! is forward-backward (0) or simulated backward
                     !! Euler (1) [nondim].  0 is often used.
  logical :: debug   !< If true, write verbose checksums for debugging purposes.
  logical :: debug_OBC !< If true, do additional calls resetting values to help verify the correctness
                       !! of the open boundary condition code.
  logical :: fpmix = .false.                 !< If true, applies profiles of momentum flux magnitude and direction.
  logical :: module_is_initialized = .false. !< Record whether this module has been initialized.
  logical :: visc_rem_dt_bug = .true. !< If true, recover a bug that uses dt_pred rather than dt for vertvisc_rem
                                      !! at the end of predictor.

  !>@{ Diagnostic IDs
  !  integer :: id_uold   = -1, id_vold   = -1
  integer :: id_uh     = -1, id_vh     = -1
  integer :: id_umo    = -1, id_vmo    = -1
  integer :: id_umo_2d = -1, id_vmo_2d = -1
  integer :: id_PFu    = -1, id_PFv    = -1
  integer :: id_CAu    = -1, id_CAv    = -1
  integer :: id_ueffA = -1, id_veffA = -1
  ! integer :: id_hf_PFu    = -1, id_hf_PFv    = -1
  integer :: id_h_PFu    = -1, id_h_PFv    = -1
  integer :: id_hf_PFu_2d = -1, id_hf_PFv_2d = -1
  integer :: id_intz_PFu_2d = -1, id_intz_PFv_2d = -1
  integer :: id_PFu_visc_rem = -1, id_PFv_visc_rem = -1
  ! integer :: id_hf_CAu    = -1, id_hf_CAv    = -1
  integer :: id_h_CAu    = -1, id_h_CAv    = -1
  integer :: id_hf_CAu_2d = -1, id_hf_CAv_2d = -1
  integer :: id_intz_CAu_2d = -1, id_intz_CAv_2d = -1
  integer :: id_CAu_visc_rem = -1, id_CAv_visc_rem = -1
  integer :: id_deta_dt = -1

  ! Split scheme only.
  integer :: id_uav        = -1, id_vav        = -1
  integer :: id_u_BT_accel = -1, id_v_BT_accel = -1
  ! integer :: id_hf_u_BT_accel    = -1, id_hf_v_BT_accel    = -1
  integer :: id_h_u_BT_accel    = -1, id_h_v_BT_accel    = -1
  integer :: id_hf_u_BT_accel_2d = -1, id_hf_v_BT_accel_2d = -1
  integer :: id_intz_u_BT_accel_2d = -1, id_intz_v_BT_accel_2d = -1
  integer :: id_u_BT_accel_visc_rem    = -1, id_v_BT_accel_visc_rem    = -1
  !>@}

  type(diag_ctrl), pointer       :: diag => NULL() !< A structure that is used to regulate the
                                         !! timing of diagnostic output.
  type(accel_diag_ptrs), pointer :: ADp => NULL()  !< A structure pointing to the various
                                         !! accelerations in the momentum equations,
                                         !! which can later be used to calculate
                                         !! derived diagnostics like energy budgets.
  type(accel_diag_ptrs), pointer :: AD_pred => NULL() !< A structure pointing to the various
                                         !! predictor step accelerations in the momentum equations,
                                         !! which can be used to debug truncations.
  type(cont_diag_ptrs), pointer  :: CDp => NULL()  !< A structure with pointers to various
                                         !! terms in the continuity equations,
                                         !! which can later be used to calculate
                                         !! derived diagnostics like energy budgets.

  ! The remainder of the structure points to child subroutines' control structures.
  !> A pointer to the horizontal viscosity control structure
  type(hor_visc_CS) :: hor_visc
  !> A pointer to the continuity control structure
  type(continuity_CS) :: continuity_CSp
  !> The CoriolisAdv control structure
  type(CoriolisAdv_CS) :: CoriolisAdv
  !> A pointer to the PressureForce control structure
  type(PressureForce_CS) :: PressureForce_CSp
  !> A pointer to a structure containing interface height diffusivities
  type(vertvisc_CS),      pointer :: vertvisc_CSp      => NULL()
  !> A pointer to the set_visc control structure
  type(set_visc_CS),      pointer :: set_visc_CSp      => NULL()
  !> A pointer to the barotropic stepping control structure
  type(barotropic_CS) :: barotropic_CSp
  !> A pointer to the SAL control structure
  type(SAL_CS) :: SAL_CSp
  !> A pointer to the tidal forcing control structure
  type(tidal_forcing_CS) :: tides_CSp
  !> A pointer to the harmonic analysis control structure
  type(harmonic_analysis_CS) :: HA_CSp
  !> A pointer to the ALE control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

  type(ocean_OBC_type),   pointer :: OBC => NULL() !< A pointer to an open boundary
     !! condition type that specifies whether, where, and  what open boundary
     !! conditions are used.  If no open BCs are used, this pointer stays
     !! nullified.  Flather OBCs use open boundary_CS as well.
  !> A pointer to the update_OBC control structure
  type(update_OBC_CS),    pointer :: update_OBC_CSp => NULL()

  type(group_pass_type) :: pass_eta      !< Structure for group halo pass
  type(group_pass_type) :: pass_visc_rem !< Structure for group halo pass
  type(group_pass_type) :: pass_uvp      !< Structure for group halo pass
  type(group_pass_type) :: pass_uv_inst  !< Structure for group halo pass
  type(group_pass_type) :: pass_hp_uv    !< Structure for group halo pass
  type(group_pass_type) :: pass_hp_uhvh  !< Structure for group halo pass
  type(group_pass_type) :: pass_h_uv     !< Structure for group halo pass

end type MOM_dyn_split_RK2b_CS


public step_MOM_dyn_split_RK2b
public register_restarts_dyn_split_RK2b
public initialize_dyn_split_RK2b
public remap_dyn_split_RK2b_aux_vars
public end_dyn_split_RK2b

!>@{ CPU time clock IDs
integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_mom_update
integer :: id_clock_continuity, id_clock_thick_diff
integer :: id_clock_btstep, id_clock_btcalc, id_clock_btforce
integer :: id_clock_pass
!>@}


  interface
module subroutine step_MOM_dyn_split_RK2b(u_av, v_av, h, tv, visc, Time_local, dt, forces, &
                                   p_surf_begin, p_surf_end, uh, vh, uhtr, vhtr, eta_av, &
                                   G, GV, US, CS, calc_dtbt, VarMix, MEKE, thickness_diffuse_CSp, pbv, Waves)
  type(ocean_grid_type),             intent(inout) :: G            !< Ocean grid structure
  type(verticalGrid_type),           intent(in)    :: GV           !< Ocean vertical grid structure
  type(unit_scale_type),             intent(in)    :: US           !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             target, intent(inout) :: u_av         !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             target, intent(inout) :: v_av         !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                     intent(inout) :: h            !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),             intent(in)    :: tv           !< Thermodynamic type
  type(vertvisc_type),               intent(inout) :: visc         !< Vertical visc, bottom drag, and related
  type(time_type),                   intent(in)    :: Time_local   !< Model time at end of time step
  real,                              intent(in)    :: dt           !< Baroclinic dynamics time step [T ~> s]
  type(mech_forcing),                intent(in)    :: forces       !< A structure with the driving mechanical forces
  real, dimension(:,:),              pointer       :: p_surf_begin !< Surface pressure at the start of this dynamic
                                                                   !! time step [R L2 T-2 ~> Pa]
  real, dimension(:,:),              pointer       :: p_surf_end   !< Surface pressure at the end of this dynamic
                                                                   !! time step [R L2 T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             target, intent(inout) :: uh           !< Zonal volume or mass transport
                                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             target, intent(inout) :: vh           !< Meridional volume or mass transport
                                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                     intent(inout) :: uhtr         !< Accumulated zonal volume or mass transport
                                                                   !! since last tracer advection [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                     intent(inout) :: vhtr         !< Accumulated meridional volume or mass transport
                                                                   !! since last tracer advection [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: eta_av       !< Free surface height or column mass
                                                                   !! averaged over time step [H ~> m or kg m-2]
  type(MOM_dyn_split_RK2b_CS),       pointer       :: CS           !< Module control structure
  logical,                           intent(in)    :: calc_dtbt    !< If true, recalculate the barotropic time step
  type(VarMix_CS),                   intent(inout) :: VarMix       !< Variable mixing control structure
  type(MEKE_type),                   intent(inout) :: MEKE         !< MEKE fields
  type(thickness_diffuse_CS),        intent(inout) :: thickness_diffuse_CSp !< Pointer to a structure containing
                                                                   !! interface height diffusivities
  type(porous_barrier_type),         intent(in)    :: pbv          !< porous barrier fractional cell metrics
  type(wave_parameters_CS), optional, pointer      :: Waves        !< A pointer to a structure containing
                                                                   !! fields related to the surface wave conditions

  ! local variables

                                                        ! of each layer calculated by the non-barotropic
                                                        ! part of the model [L T-2 ~> m s-2]
                                                        ! of each layer calculated by the non-barotropic
                                                        ! part of the model [L T-2 ~> m s-2]

                                                     ! and end of a time step [H ~> m or kg m-2]

                                ! obtained using the initial velocities [H L2 T-1 ~> m3 s-1 or kg s-1]
                                ! obtained using the initial velocities [H L2 T-1 ~> m3 s-1 or kg s-1]

                                               ! or column mass [H ~> m or kg m-2]
                                               ! height or column mass [H T-1 ~> m s-1 or kg m-2 s-1]

                                ! saved for use in the Flather open boundary condition code [L T-1 ~> m s-1]
                                ! saved for use in the Flather open boundary condition code [L T-1 ~> m s-1]

  ! GMM, TODO: make these allocatable?
  ! real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uold ! u-velocity before vert_visc is applied, for fpmix
  !                                                    !                                      [L T-1 ~> m s-1]
  ! real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vold ! v-velocity before vert_visc is applied, for fpmix
  !                                                    !                                      [L T-1 ~> m s-1]
                      ! [H T2 R-1 L-2 ~> m Pa-1 or kg m-2 Pa-1]
end subroutine step_MOM_dyn_split_RK2b
module subroutine register_restarts_dyn_split_RK2b(HI, GV, US, param_file, CS, restart_CS, uh, vh)
  type(hor_index_type),          intent(in)    :: HI         !< Horizontal index structure
  type(verticalGrid_type),       intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< parameter file
  type(MOM_dyn_split_RK2b_CS),   pointer       :: CS         !< module control structure
  type(MOM_restart_CS),          intent(inout) :: restart_CS !< MOM restart control structure
  real, dimension(SZIB_(HI),SZJ_(HI),SZK_(GV)), &
                         target, intent(inout) :: uh !< zonal volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(HI),SZJB_(HI),SZK_(GV)), &
                         target, intent(inout) :: vh !< merid volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]



end subroutine register_restarts_dyn_split_RK2b
module subroutine remap_dyn_split_RK2b_aux_vars(G, GV, CS, h_old_u, h_old_v, h_new_u, h_new_v, ALE_CSp)
  type(ocean_grid_type),            intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV       !< ocean vertical grid structure
  type(MOM_dyn_split_RK2b_CS),      pointer       :: CS       !< module control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_old_u  !< Source grid thickness at zonal
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: h_old_v  !< Source grid thickness at meridional
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_new_u  !< Destination grid thickness at zonal
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: h_new_v  !< Destination grid thickness at meridional
                                                              !! velocity points [H ~> m or kg m-2]
  type(ALE_CS),                     pointer       :: ALE_CSp  !< ALE control structure to use when remapping

end subroutine remap_dyn_split_RK2b_aux_vars
module subroutine initialize_dyn_split_RK2b(u, v, h, tv, uh, vh, eta, Time, G, GV, US, param_file, &
                      diag, CS, HA_CSp, restart_CS, dt, Accel_diag, Cont_diag, MIS, &
                      VarMix, MEKE, thickness_diffuse_CSp,                  &
                      OBC, update_OBC_CSp, ALE_CSp, set_visc, &
                      visc, dirs, ntrunc, pbv, calc_dtbt, cont_stencil, dyn_h_stencil)
  type(ocean_grid_type),            intent(inout) :: G          !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US         !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(inout) :: u          !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(inout) :: v          !< merid velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                    intent(inout) :: h          !< layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv         !< Thermodynamic type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            target, intent(inout) :: uh    !< zonal volume/mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            target, intent(inout) :: vh    !< merid volume/mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: eta        !< free surface height or column mass [H ~> m or kg m-2]
  type(time_type),          target, intent(in)    :: Time       !< current model time
  type(param_file_type),            intent(in)    :: param_file !< parameter file for parsing
  type(diag_ctrl),          target, intent(inout) :: diag       !< to control diagnostics
  type(MOM_dyn_split_RK2b_CS),      pointer       :: CS         !< module control structure
  type(harmonic_analysis_CS),       pointer       :: HA_CSp     !< A pointer to the control structure of the
                                                                !! harmonic analysis module
  type(MOM_restart_CS),             intent(inout) :: restart_CS !< MOM restart control structure
  real,                             intent(in)    :: dt         !< time step [T ~> s]
  type(accel_diag_ptrs),    target, intent(inout) :: Accel_diag !< points to momentum equation terms for
                                                                !! budget analysis
  type(cont_diag_ptrs),     target, intent(inout) :: Cont_diag  !< points to terms in continuity equation
  type(ocean_internal_state),       intent(inout) :: MIS        !< "MOM6 internal state" used to pass
                                                                !! diagnostic pointers
  type(VarMix_CS),                  intent(inout) :: VarMix     !< points to spatially variable viscosities
  type(MEKE_type),                  intent(inout) :: MEKE       !< MEKE fields
  type(thickness_diffuse_CS),       intent(inout) :: thickness_diffuse_CSp !< Pointer to the control structure
                                                                !! used for the isopycnal height diffusive transport.
  type(ocean_OBC_type),             pointer       :: OBC        !< points to OBC related fields
  type(update_OBC_CS),              pointer       :: update_OBC_CSp !< points to OBC update related fields
  type(ALE_CS),                     pointer       :: ALE_CSp    !< points to ALE control structure
  type(set_visc_CS),        target, intent(in)    :: set_visc   !< set_visc control structure
  type(vertvisc_type),              intent(inout) :: visc       !< vertical viscosities, bottom drag, and related
  type(directories),                intent(in)    :: dirs       !< contains directory paths
  integer, target,                  intent(inout) :: ntrunc     !< A target for the variable that records
                                                                !! the number of times the velocity is
                                                                !! truncated (this should be 0).
  logical,                          intent(out)   :: calc_dtbt  !< If true, recalculate the barotropic time step
  type(porous_barrier_type),        intent(in)    :: pbv        !< porous barrier fractional cell metrics
  integer,                          intent(out)   :: cont_stencil !< The stencil for thickness
                                                                !! from the continuity solver.
  integer,                          intent(out)   :: dyn_h_stencil !< The stencil for thickness for the
                                                                !! dynamics based on the continuity
                                                                !! solver and Coriolis scheme.

  ! local variables
  ! This include declares and sets the variable "version".
                          ! recreate the bugs, or if false bugs are only used if actively selected.

end subroutine initialize_dyn_split_RK2b
module subroutine end_dyn_split_RK2b(CS)
  type(MOM_dyn_split_RK2b_CS), pointer :: CS  !< module control structure

end subroutine end_dyn_split_RK2b
  end interface

end module MOM_dynamics_split_RK2b
