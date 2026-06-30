! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Barotropic solver
module MOM_barotropic

use MOM_checksums, only : chksum0
use MOM_coms,      only : any_across_PEs
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_debugging, only : hchksum, uvchksum, Bchksum
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl, enable_averaging, enable_averages
use MOM_domains, only : min_across_PEs, clone_MOM_domain, deallocate_MOM_domain
use MOM_domains, only : To_All, Scalar_Pair, AGRID, CORNER, MOM_domain_type
use MOM_domains, only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains, only : start_group_pass, complete_group_pass, pass_var, pass_vector
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_harmonic_analysis, only : HA_accum, harmonic_analysis_CS
use MOM_hor_index, only : hor_index_type
use MOM_io, only : vardesc, var_desc, MOM_read_data, slasher, NORTH_FACE, EAST_FACE
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, open_boundary_query
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_segment_type
use MOM_restart, only : register_restart_field, register_restart_pair
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_self_attr_load, only : scalar_SAL_sensitivity
use MOM_self_attr_load, only : SAL_CS
use MOM_streaming_filter, only : Filt_register, Filt_init, Filt_accum, Filter_CS
use MOM_time_manager, only : time_type, real_to_time, operator(+), operator(-)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, alloc_bt_cont_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_variables, only : accel_diag_ptrs
use MOM_wave_drag, only : wave_drag_init, wave_drag_calc, wave_drag_CS

implicit none ; private

#include <MOM_memory.h>
#ifdef STATIC_MEMORY_
#  ifndef BTHALO_
#    define BTHALO_ 0
#  endif
#  define WHALOI_ MAX(BTHALO_-NIHALO_,0)
#  define WHALOJ_ MAX(BTHALO_-NJHALO_,0)
#  define NIMEMW_   1-WHALOI_:NIMEM_+WHALOI_
#  define NJMEMW_   1-WHALOJ_:NJMEM_+WHALOJ_
#  define NIMEMBW_  -WHALOI_:NIMEM_+WHALOI_
#  define NJMEMBW_  -WHALOJ_:NJMEM_+WHALOJ_
#  define SZIW_(G)  NIMEMW_
#  define SZJW_(G)  NJMEMW_
#  define SZIBW_(G) NIMEMBW_
#  define SZJBW_(G) NJMEMBW_
#else
#  define NIMEMW_   :
#  define NJMEMW_   :
#  define NIMEMBW_  :
#  define NJMEMBW_  :
#  define SZIW_(G)  G%isdw:G%iedw
#  define SZJW_(G)  G%jsdw:G%jedw
#  define SZIBW_(G) G%isdw-1:G%iedw
#  define SZJBW_(G) G%jsdw-1:G%jedw
#endif

public btcalc, bt_mass_source, btstep, barotropic_init, barotropic_end
public register_barotropic_restarts, set_dtbt, barotropic_get_tav

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The barotropic stepping open boundary condition type
type, private :: BT_OBC_type
  real, allocatable :: Cg_u(:,:)  !< The external wave speed at u-points [L T-1 ~> m s-1].
  real, allocatable :: Cg_v(:,:)  !< The external wave speed at u-points [L T-1 ~> m s-1].
  real, allocatable :: dZ_u(:,:)  !< The total vertical column extent at the u-points [Z ~> m].
  real, allocatable :: dZ_v(:,:)  !< The total vertical column extent at the v-points [Z ~> m].
  real, allocatable :: uhbt(:,:)  !< The zonal barotropic thickness fluxes specified
                                  !! for open boundary conditions (if any) [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, allocatable :: vhbt(:,:)  !< The meridional barotropic thickness fluxes specified
                                  !! for open boundary conditions (if any) [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, allocatable :: ubt_outer(:,:) !< The zonal velocities just outside the domain,
                                  !! as set by the open boundary conditions [L T-1 ~> m s-1].
  real, allocatable :: vbt_outer(:,:) !< The meridional velocities just outside the domain,
                                  !! as set by the open boundary conditions [L T-1 ~> m s-1].
  real, allocatable :: SSH_outer_u(:,:) !< The surface height outside of the domain
                                  !! at a u-point with an open boundary condition [Z ~> m].
  real, allocatable :: SSH_outer_v(:,:) !< The surface height outside of the domain
                                  !! at a v-point with an open boundary condition [Z ~> m].
  integer, allocatable :: u_OBC_type(:,:) !< An integer encoding the type and direction of u-point OBCs
  integer, allocatable :: v_OBC_type(:,:) !< An integer encoding the type and direction of v-point OBCs
  logical :: u_OBCs_on_PE !< True if this PE has an open boundary at any u-points.
  logical :: v_OBCs_on_PE !< True if this PE has an open boundary at any v-points.
  !>@{ Index ranges on the local PE for the open boundary conditions in various directions
  integer :: Is_u_W_obc, Ie_u_W_obc, js_u_W_obc, je_u_W_obc
  integer :: Is_u_E_obc, Ie_u_E_obc, js_u_E_obc, je_u_E_obc
  integer :: is_v_S_obc, ie_v_S_obc, Js_v_S_obc, Je_v_S_obc
  integer :: is_v_N_obc, ie_v_N_obc, Js_v_N_obc, Je_v_N_obc
  !>@}

  type(group_pass_type) :: pass_uv   !< Structure for group halo pass of vectors
  type(group_pass_type) :: scalar_pass  !< Structure for group halo pass of scalars
end type BT_OBC_type

integer, parameter :: SPECIFIED_OBC = 1 !< An integer used to encode a specified OBC point
integer, parameter :: FLATHER_OBC = 2   !< An integer used to encode a Flather OBC point
integer, parameter :: GRADIENT_OBC = 4  !< An integer used to encode a gradient OBC point

!> The barotropic stepping control structure
type, public :: barotropic_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: frhatu
          !< The fraction of the total column thickness interpolated to u grid points in each layer [nondim].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: frhatv
          !< The fraction of the total column thickness interpolated to v grid points in each layer [nondim].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: IDatu
          !< Inverse of the total thickness at u grid points [H-1 ~> m-1 or m2 kg-1].
  real, allocatable, dimension(:,:) :: lin_drag_u
          !< A spatially varying linear drag coefficient acting on the zonal barotropic flow
          !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real, allocatable, dimension(:,:) :: ubt_IC
          !< The barotropic solvers estimate of the zonal velocity that will be the initial
          !! condition for the next call to btstep [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: ubtav
          !< The barotropic zonal velocity averaged over the baroclinic time step [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: IDatv
          !< Inverse of the basin depth at v grid points [Z-1 ~> m-1].
  real, allocatable, dimension(:,:) :: lin_drag_v
          !< A spatially varying linear drag coefficient acting on the zonal barotropic flow
          !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real, allocatable, dimension(:,:) :: vbt_IC
          !< The barotropic solvers estimate of the zonal velocity that will be the initial
          !! condition for the next call to btstep [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: vbtav
          !< The barotropic meridional velocity averaged over the  baroclinic time step [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: eta_cor
          !< The difference between the free surface height from the barotropic calculation and the sum
          !! of the layer thicknesses. This difference is imposed as a forcing term in the barotropic
          !! calculation over a baroclinic timestep [H ~> m or kg m-2].
  real, allocatable, dimension(:,:) :: eta_cor_bound
          !< A limit on the rate at which eta_cor can be applied while avoiding instability
          !! [H T-1 ~> m s-1 or kg m-2 s-1]. This is only used if CS%bound_BT_corr is true.
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: &
    ua_polarity, &  !< Test vector components for checking grid polarity [nondim]
    va_polarity, &  !< Test vector components for checking grid polarity [nondim]
    bathyT          !< A copy of bathyT (ocean bottom depth) with wide halos [Z ~> m]
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: IareaT
                    !<   This is a copy of G%IareaT with wide halos, but will
                    !! still utilize the macro IareaT when referenced, [L-2 ~> m-2].
  real ALLOCABLE_, dimension(NIMEMBW_,NJMEMW_) :: &
    dy_Cu, &        !<   A copy of G%dy_Cu with wide halos [L ~> m].
    IdxCu, &        !<   A copy of G%IdxCu with wide halos [L-1 ~> m-1].
    OBCmask_u       !< An array to multiplicatively mask out changes at OBC points, 0 or 1 [nondim]
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMBW_) :: &
    dx_Cv, &        !<   A copy of G%dx_Cv with wide halos [L ~> m].
    IdyCv, &        !<   A copy of G%IdyCv with wide halos [L-1 ~> m-1].
    OBCmask_v       !< An array to multiplicatively mask out changes at OBC points, 0 or 1 [nondim]
  real, allocatable, dimension(:,:) :: &
    D_u_Cor, &      !<   A simply averaged depth at u points recast as a thickness [H ~> m or kg m-2]
    D_v_Cor, &      !<   A simply averaged depth at v points recast as a thickness [H ~> m or kg m-2]
    q_D             !< f / D at PV points [T-1 H-1 ~> s-1 m-1 or m2 s-1 kg-1]
  real, allocatable, dimension(:,:,:) :: &
    q_wt            !< The area weights for the thicknesses around a corner point to be used when
                    !! calculating PV for use in the Coriolis term, taking OBCs into account [L2 ~> m2].
                    !! The order of the 4 values at a point is the order in which the neighboring
                    !! tracer points occur in memory, i.e. SW, SE, NW then NE.
  real, allocatable :: frhatu1(:,:,:)  !< Predictor step values of frhatu stored for diagnostics [nondim]
  real, allocatable :: frhatv1(:,:,:)  !< Predictor step values of frhatv stored for diagnostics [nondim]
  real, allocatable :: IareaT_OBCmask(:,:)  !< If non-zero, work on given points [L-2 ~> m-2].

  type(BT_OBC_type) :: BT_OBC !< A structure with all of this modules fields
                              !! for applying open boundary conditions.

  real    :: dtbt            !< The barotropic time step [T ~> s].
  real    :: dtbt_fraction   !<   The fraction of the maximum time-step that
                             !! should used [nondim].  The default is 0.98.
  real    :: dtbt_max        !<   The maximum stable barotropic time step [T ~> s].
  real    :: dt_bt_filter    !<   The time-scale over which the barotropic mode solutions are
                             !! filtered [T ~> s] if positive, or as a fraction of DT if
                             !! negative [nondim].  This can never be taken to be longer than 2*dt.
                             !! Set this to 0 to apply no filtering.
  integer :: nstep_last = 0  !< The number of barotropic timesteps per baroclinic
                             !! time step the last time btstep was called.
  real    :: bebt            !< A nondimensional number, from 0 to 1, that
                             !! determines the gravity wave time stepping scheme [nondim].
                             !! 0.0 gives a forward-backward scheme, while 1.0
                             !! give backward Euler. In practice, bebt should be
                             !! of order 0.2 or greater.
  real    :: Rho_BT_lin      !< A density that is used to convert total water column thicknesses
                             !! into mass in non-Boussinesq mode with linearized options in the
                             !! barotropic solver or when estimating the stable barotropic timestep
                             !! without access to the full baroclinic model state [R ~> kg m-3]
  logical :: split           !< If true, use the split time stepping scheme.
  logical :: bound_BT_corr   !< If true, the magnitude of the fake mass source
                             !! in the barotropic equation that drives the two
                             !! estimates of the free surface height toward each
                             !! other is bounded to avoid driving corrective
                             !! velocities that exceed MAXCFL_BT_CONT.
  logical :: gradual_BT_ICs  !< If true, adjust the initial conditions for the
                             !! barotropic solver to the values from the layered
                             !! solution over a whole timestep instead of
                             !! instantly.  This is a decent approximation to the
                             !! inclusion of sum(u dh_dt) while also correcting
                             !! for truncation errors.
  logical :: Sadourny        !< If true, the Coriolis terms are discretized
                             !! with Sadourny's energy conserving scheme,
                             !! otherwise the Arakawa & Hsu scheme is used.  If
                             !! the deformation radius is not resolved Sadourny's
                             !! scheme should probably be used.
  logical :: integral_bt_cont !< If true, use the time-integrated velocity over the barotropic steps
                             !! to determine the integrated transports used to update the continuity
                             !! equation.  Otherwise the transports are the sum of the transports
                             !! based on a series of instantaneous velocities and the BT_CONT_TYPE
                             !! for transports.  This is only valid if a BT_CONT_TYPE is used.
  logical :: bt_adjust_src_for_filter !< If true, increases the rate at which BT mass sources are
                             !! applied so that they are all used up before the steps within the
                             !! filtering period start. This avoids the mass sink driving the SSH
                             !! below the bottom during the period of filtering.
  logical :: bt_limit_integral_transport !< If true, limit the time-integrated transports by the
                             !! initial volume accounting for sinks of mass.
  logical :: integral_OBCs   !< This is true if integral_bt_cont is true and there are open boundary
                             !! conditions being applied somewhere in the global domain.
  logical :: Nonlinear_continuity !< If true, the barotropic continuity equation
                             !! uses the full ocean thickness for transport.
  integer :: Nonlin_cont_update_period !< The number of barotropic time steps
                             !! between updates to the face area, or 0 only to
                             !! update at the start of a call to btstep.  The
                             !! default is 1.
  logical :: BT_project_velocity !< If true, step the barotropic velocity first
                             !! and project out the velocity tendency by 1+BEBT
                             !! when calculating the transport.  The default
                             !! (false) is to use a predictor continuity step to
                             !! find the pressure field, and then do a corrector
                             !! continuity step using a weighted average of the
                             !! old and new velocities, with weights of (1-BEBT) and BEBT.
  logical :: nonlin_stress   !< If true, use the full depth of the ocean at the start of the
                             !! barotropic step when calculating the surface stress contribution to
                             !! the barotropic accelerations.  Otherwise use the depth based on bathyT.
  real    :: BT_Coriolis_scale !< A factor by which the barotropic Coriolis acceleration anomaly
                             !! terms are scaled [nondim].
  integer :: answer_date     !< The vintage of the expressions in the barotropic solver.
                             !! Values below 20190101 recover the answers from the end of 2018,
                             !! while higher values use more efficient or general expressions.

  logical :: dynamic_psurf   !< If true, add a dynamic pressure due to a viscous
                             !! ice shelf, for instance.
  real    :: Dmin_dyn_psurf  !< The minimum total thickness to use in limiting the size
                             !! of the dynamic surface pressure for stability [H ~> m or kg m-2].
  real    :: ice_strength_length  !< The length scale at which the damping rate
                             !! due to the ice strength should be the same as if
                             !! a Laplacian were applied [L ~> m].
  real    :: const_dyn_psurf !< The constant that scales the dynamic surface
                             !! pressure [nondim].  Stable values are < ~1.0.
                             !! The default is 0.9.
  logical :: calculate_SAL   !< If true, calculate self-attraction and loading.
  logical :: tidal_sal_bug   !< If true, the tidal self-attraction and loading anomaly in the
                             !! barotropic solver has the wrong sign, replicating a long-standing
                             !! bug.
  real    :: G_extra         !< A nondimensional factor by which gtot is enhanced [nondim].
  integer :: hvel_scheme     !< An integer indicating how the thicknesses at
                             !! velocity points are calculated. Valid values are
                             !! given by the parameters defined below:
                             !!   HARMONIC, ARITHMETIC, HYBRID, and FROM_BT_CONT
  logical :: strong_drag     !< If true, use a stronger estimate of the retarding
                             !! effects of strong bottom drag.
  logical :: rescale_strong_drag !< If true, reduce the barotropic contribution to the layer
                             !! accelerations to account for the difference between the forces that
                             !! can be counteracted  by the stronger drag with BT_STRONG_DRAG and the
                             !! average of the layer viscous remnants after a baroclinic timestep.
  logical :: linear_wave_drag  !< If true, apply a linear drag to the barotropic
                             !! velocities, using rates set by lin_drag_u & _v
                             !! divided by the depth of the ocean.
  logical :: linearized_BT_PV  !< If true, the PV and interface thicknesses used
                             !! in the barotropic Coriolis calculation is time
                             !! invariant and linearized.
  logical :: use_filter      !< If true, use streaming band-pass filter to detect the
                             !! instantaneous tidal signals in the simulation.
  logical :: linear_freq_drag  !< If true, apply a linear frequency-dependent drag to the tidal
                             !! velocities. The streaming band-pass filter must be turned on.
  logical :: use_wide_halos  !< If true, use wide halos and march in during the
                             !! barotropic time stepping for efficiency.
  integer :: min_stencil     !< The minimum stencil width to use with the wide halo iterations.
                             !! A nonzero value may reflect the distribution of OBC faces or it
                             !! may be useful for debugging purposes.
  logical :: clip_velocity   !< If true, limit any velocity components that are
                             !! are large enough for a CFL number to exceed
                             !! CFL_trunc.  This should only be used as a
                             !! desperate debugging measure.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  logical :: debug_bt        !< If true, write verbose checksums from within the barotropic
                             !! time-stepping loop for debugging purposes.
  logical :: debug_wide_halos !< If true, write the checksums on the full wide halos.   Otherwise
                             !! only the output for the final computational domain is written.
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: maxvel          !< Velocity components greater than maxvel are
                             !! truncated to maxvel [L T-1 ~> m s-1].
  real    :: CFL_trunc       !< If clip_velocity is true, velocity components will
                             !! be truncated when they are large enough that the
                             !! corresponding CFL number exceeds this value [nondim].
  real    :: maxCFL_BT_cont  !< The maximum permitted CFL number associated with the
                             !! barotropic accelerations from the summed velocities
                             !! times the time-derivatives of thicknesses [nondim].  The
                             !! default is 0.1, and there will probably be real
                             !! problems if this were set close to 1.
  logical :: BT_cont_bounds  !< If true, use the BT_cont_type variables to set limits
                             !! on the magnitude of the corrective mass fluxes.
  logical :: visc_rem_u_uh0  !< If true, use the viscous remnants when estimating
                             !! the barotropic velocities that were used to
                             !! calculate uh0 and vh0.  False is probably the
                             !! better choice.
  logical :: adjust_BT_cont  !< If true, adjust the curve fit to the BT_cont type
                             !! that is used by the barotropic solver to match the
                             !! transport about which the flow is being linearized.
  logical :: use_old_coriolis_bracket_bug !< If True, use an order of operations
                             !! that is not bitwise rotationally symmetric in the
                             !! meridional Coriolis term of the barotropic solver.
  logical :: tidal_sal_flather !< Apply adjustment to external gravity wave speed
                             !! consistent with tidal self-attraction and loading
                             !! used within the barotropic solver
  logical :: wt_uv_bug = .true. !< If true, recover a bug that wt_[uv] that is not normalized.
  logical :: exterior_OBC_bug = .true. !< If true, recover a bug with boundary conditions
                             !! inside the domain.
  logical :: interior_OBC_PV !< If true, use only interior ocean points at OBCs to specify the PV
                             !!  used in the barotropic Coriolis anomalies.  Otherwise the
                             !! calculation relies on bathymetry and eta being projected outward
                             !! across OBCs.  Unfortunately, this option does change answers near
                             !! convex (peninsula-type) pairs of OBC segments.
  type(time_type), pointer :: Time  => NULL() !< A pointer to the ocean models clock.
  type(diag_ctrl), pointer :: diag => NULL()  !< A structure that is used to regulate
                             !! the timing of diagnostic output.
  type(MOM_domain_type), pointer :: BT_Domain => NULL()  !< Barotropic MOM domain
  type(hor_index_type), pointer :: debug_BT_HI => NULL() !< debugging copy of horizontal index_type
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< Control structure for SAL
  type(harmonic_analysis_CS), pointer :: HA_CSp => NULL() !< Control structure for harmonic analysis
  type(Filter_CS) :: Filt_CS_u, & !< Control structures for the streaming band-pass filter of ubt
                     Filt_CS_v    !< Control structures for the streaming band-pass filter of vbt
  type(wave_drag_CS) :: Drag_CS !< Control structures for the frequency-dependent drag
  logical :: module_is_initialized = .false.  !< If true, module has been initialized

  integer :: isdw !< The lower i-memory limit for the wide halo arrays.
  integer :: iedw !< The upper i-memory limit for the wide halo arrays.
  integer :: jsdw !< The lower j-memory limit for the wide halo arrays.
  integer :: jedw !< The upper j-memory limit for the wide halo arrays.

  type(group_pass_type) :: pass_q_DCor !< Handle for a group halo pass
  type(group_pass_type) :: pass_gtot !< Handle for a group halo pass
  type(group_pass_type) :: pass_tmp_uv !< Handle for a group halo pass
  type(group_pass_type) :: pass_eta_bt_rem !< Handle for a group halo pass
  type(group_pass_type) :: pass_force_hbt0_Cor_ref !< Handle for a group halo pass
  type(group_pass_type) :: pass_Dat_uv !< Handle for a group halo pass
  type(group_pass_type) :: pass_eta_ubt !< Handle for a group halo pass
  type(group_pass_type) :: pass_etaav !< Handle for a group halo pass
  type(group_pass_type) :: pass_ubt_Cor !< Handle for a group halo pass
  type(group_pass_type) :: pass_ubta_uhbta !< Handle for a group halo pass
  type(group_pass_type) :: pass_e_anom !< Handle for a group halo pass
  type(group_pass_type) :: pass_SpV_avg !< Handle for a group halo pass

  !>@{ Diagnostic IDs
  integer :: id_PFu_bt = -1, id_PFv_bt = -1, id_Coru_bt = -1, id_Corv_bt = -1
  integer :: id_LDu_bt = -1, id_LDv_bt = -1, id_eta_cor = -1
  integer :: id_ubtforce = -1, id_vbtforce = -1, id_uaccel = -1, id_vaccel = -1
  integer :: id_visc_rem_u = -1, id_visc_rem_v = -1, id_bt_rem_u = -1, id_bt_rem_v = -1
  integer :: id_ubt = -1, id_vbt = -1, id_eta_bt = -1, id_ubtav = -1, id_vbtav = -1
  integer :: id_ubt_st = -1, id_vbt_st = -1, id_eta_st = -1
  integer :: id_ubtdt = -1, id_vbtdt = -1
  integer :: id_ubt_hifreq = -1, id_vbt_hifreq = -1, id_eta_hifreq = -1
  integer :: id_uhbt_hifreq = -1, id_vhbt_hifreq = -1, id_eta_pred_hifreq = -1
  integer :: id_etaPF_hifreq = -1, id_etaPF_anom = -1
  integer :: id_gtotn = -1, id_gtots = -1, id_gtote = -1, id_gtotw = -1
  integer :: id_uhbt = -1, id_frhatu = -1, id_vhbt = -1, id_frhatv = -1
  integer :: id_frhatu1 = -1, id_frhatv1 = -1

  integer :: id_BTC_FA_u_EE = -1, id_BTC_FA_u_E0 = -1, id_BTC_FA_u_W0 = -1, id_BTC_FA_u_WW = -1
  integer :: id_BTC_ubt_EE = -1, id_BTC_ubt_WW = -1
  integer :: id_BTC_FA_v_NN = -1, id_BTC_FA_v_N0 = -1, id_BTC_FA_v_S0 = -1, id_BTC_FA_v_SS = -1
  integer :: id_BTC_vbt_NN = -1, id_BTC_vbt_SS = -1
  integer :: id_BTC_FA_u_rat0 = -1, id_BTC_FA_v_rat0 = -1, id_BTC_FA_h_rat0 = -1
  integer :: id_uhbt0 = -1, id_vhbt0 = -1
  integer :: id_SSH_u_OBC = -1, id_SSH_v_OBC = -1, id_ubt_OBC = -1, id_vbt_OBC = -1
  !>@}

end type barotropic_CS

!> A description of the functional dependence of transport at a u-point
type, private :: local_BT_cont_u_type
  real :: FA_u_EE !< The effective open face area for zonal barotropic transport
                  !! drawing from locations far to the east [H L ~> m2 or kg m-1].
  real :: FA_u_E0 !< The effective open face area for zonal barotropic transport
                  !! drawing from nearby to the east [H L ~> m2 or kg m-1].
  real :: FA_u_W0 !< The effective open face area for zonal barotropic transport
                  !! drawing from nearby to the west [H L ~> m2 or kg m-1].
  real :: FA_u_WW !< The effective open face area for zonal barotropic transport
                  !! drawing from locations far to the west [H L ~> m2 or kg m-1].
  real :: uBT_WW  !< uBT_WW is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_u_WW.  uBT_WW must be non-negative.
  real :: uBT_EE  !< uBT_EE is a barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_u_EE. uBT_EE must be non-positive.
  real :: uh_crvW !< The curvature of face area with velocity for flow from the west [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> nondim or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: uh_crvE !< The curvature of face area with velocity for flow from the east [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> nondim or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: uh_WW   !< The zonal transport when ubt=ubt_WW [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
  real :: uh_EE   !< The zonal transport when ubt=ubt_EE [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
end type local_BT_cont_u_type

!> A description of the functional dependence of transport at a v-point
type, private :: local_BT_cont_v_type
  real :: FA_v_NN !< The effective open face area for meridional barotropic transport
                  !! drawing from locations far to the north [H L ~> m2 or kg m-1].
  real :: FA_v_N0 !< The effective open face area for meridional barotropic transport
                  !! drawing from nearby to the north [H L ~> m2 or kg m-1].
  real :: FA_v_S0 !< The effective open face area for meridional barotropic transport
                  !! drawing from nearby to the south [H L ~> m2 or kg m-1].
  real :: FA_v_SS !< The effective open face area for meridional barotropic transport
                  !! drawing from locations far to the south [H L ~> m2 or kg m-1].
  real :: vBT_SS  !< vBT_SS is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_v_SS. vBT_SS must be non-negative.
  real :: vBT_NN  !< vBT_NN is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_v_NN.  vBT_NN must be non-positive.
  real :: vh_crvS !< The curvature of face area with velocity for flow from the south [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> nondim or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: vh_crvN !< The curvature of face area with velocity for flow from the north [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> nondim or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: vh_SS   !< The meridional transport when vbt=vbt_SS [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
  real :: vh_NN   !< The meridional transport when vbt=vbt_NN [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
end type local_BT_cont_v_type

!> A container for passing around active tracer point memory limits
type, private :: memory_size_type
  !>@{ Currently active memory limits
  integer :: isdw, iedw, jsdw, jedw ! The memory limits of the wide halo arrays.
  !>@}
end type memory_size_type

!>@{ CPU time clock IDs
integer :: id_clock_sync=-1, id_clock_calc=-1
integer :: id_clock_calc_pre=-1, id_clock_calc_post=-1
integer :: id_clock_pass_step=-1, id_clock_pass_pre=-1, id_clock_pass_post=-1
!>@}

!>@{ Enumeration values for various schemes
integer, parameter :: HARMONIC        = 1
integer, parameter :: ARITHMETIC      = 2
integer, parameter :: HYBRID          = 3
integer, parameter :: FROM_BT_CONT    = 4
integer, parameter :: HYBRID_BT_CONT  = 5
character*(20), parameter :: HYBRID_STRING = "HYBRID"
character*(20), parameter :: HARMONIC_STRING = "HARMONIC"
character*(20), parameter :: ARITHMETIC_STRING = "ARITHMETIC"
character*(20), parameter :: BT_CONT_STRING = "FROM_BT_CONT"
!>@}

!> A negligible parameter which avoids division by zero, but is too small to
!! modify physical values [nondim].
real, parameter :: subroundoff = 1e-30


  interface
module subroutine btstep(U_in, V_in, eta_in, dt, bc_accel_u, bc_accel_v, forces, pbce, &
                  eta_PF_in, U_Cor, V_Cor, accel_layer_u, accel_layer_v, &
                  eta_out, uhbtav, vhbtav, G, GV, US, CS, &
                  visc_rem_u, visc_rem_v, SpV_avg, ADp, OBC, BT_cont, eta_PF_start, &
                  taux_bot, tauy_bot, uh0, vh0, u_uh0, v_vh0, etaav)
  type(ocean_grid_type),                   intent(inout) :: G       !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US      !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: U_in    !< The initial (3-D) zonal
                                                                    !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: V_in    !< The initial (3-D) meridional
                                                                    !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: eta_in  !< The initial barotropic free surface height
                                                         !! anomaly or column mass anomaly [H ~> m or kg m-2].
  real,                                      intent(in)  :: dt      !< The time increment to integrate over [T ~> s].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: bc_accel_u !< The zonal baroclinic accelerations,
                                                                       !! [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: bc_accel_v !< The meridional baroclinic accelerations,
                                                                       !! [L T-2 ~> m s-2].
  type(mech_forcing),                        intent(in)  :: forces     !< A structure with the driving mechanical forces
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: pbce       !< The baroclinic pressure anomaly in each layer
                                                         !! due to free surface height anomalies
                                                         !! [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: eta_PF_in  !< The 2-D eta field (either SSH anomaly or
                                                         !! column mass anomaly) that was used to calculate the input
                                                         !! pressure gradient accelerations (or its final value if
                                                         !! eta_PF_start is provided [H ~> m or kg m-2].
                                                         !! Note: eta_in, pbce, and eta_PF_in must have up-to-date
                                                         !! values in the first point of their halos.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: U_Cor      !< The (3-D) zonal velocities used to
                                                         !! calculate the Coriolis terms in bc_accel_u [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: V_Cor      !< The (3-D) meridional velocities used to
                                                         !! calculate the Coriolis terms in bc_accel_u [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: accel_layer_u !< The zonal acceleration of each layer due
                                                         !! to the barotropic calculation [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: accel_layer_v !< The meridional acceleration of each layer
                                                         !! due to the barotropic calculation [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: eta_out       !< The final barotropic free surface
                                                         !! height anomaly or column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(out) :: uhbtav        !< the barotropic zonal volume or mass
                                                         !! fluxes averaged through the barotropic steps
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)),         intent(out) :: vhbtav        !< the barotropic meridional volume or mass
                                                         !! fluxes averaged through the barotropic steps
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(barotropic_CS),                       intent(inout) :: CS           !< Barotropic control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: visc_rem_u    !< Both the fraction of the momentum
                                                         !! originally in a layer that remains after a time-step of
                                                         !! viscosity, and the fraction of a time-step's worth of a
                                                         !! barotropic acceleration that a layer experiences after
                                                         !! viscosity is applied, in the zonal direction [nondim].
                                                         !! Visc_rem_u is between 0 (at the bottom) and 1 (far above).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: visc_rem_v    !< Ditto for meridional direction [nondim].
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)  :: SpV_avg     !< The column average specific volume, used
                                                         !! in non-Boussinesq OBC calculations [R-1 ~> m3 kg-1]
  type(accel_diag_ptrs),                      pointer    :: ADp          !< Acceleration diagnostic pointers
  type(ocean_OBC_type),                       pointer    :: OBC          !< The open boundary condition structure.
  type(BT_cont_type),                         pointer    :: BT_cont      !< A structure with elements that describe
                                                         !! the effective open face areas as a function of barotropic
                                                         !! flow.
  real, dimension(:,:),                       pointer    :: eta_PF_start !< The eta field consistent with the pressure
                                                         !! gradient at the start of the barotropic stepping
                                                         !! [H ~> m or kg m-2].
  real, dimension(:,:),                       pointer    :: taux_bot     !< The zonal bottom frictional stress from
                                                         !! ocean to the seafloor [R L Z T-2 ~> Pa].
  real, dimension(:,:),                       pointer    :: tauy_bot     !< The meridional bottom frictional stress
                                                         !! from ocean to the seafloor [R L Z T-2 ~> Pa].
  real, dimension(:,:,:),                     pointer    :: uh0     !< The zonal layer transports at reference
                                                                    !! velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:,:),                     pointer    :: u_uh0   !< The velocities used to calculate
                                                                    !! uh0 [L T-1 ~> m s-1]
  real, dimension(:,:,:),                     pointer    :: vh0     !< The zonal layer transports at reference
                                                                    !! velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:,:),                     pointer    :: v_vh0   !< The velocities used to calculate
                                                                    !! vh0 [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: etaav        !< The free surface height or column mass
                                                         !! averaged over the barotropic integration [H ~> m or kg m-2].

  ! Local variables
                                    ! terms [L T-1 ~> m s-1].
                ! be used in calculating barotropic velocities, possibly with
                ! sums less than one due to viscous losses [nondim]
                ! used to normalize wt_u and wt_v [nondim]
                  ! averaged between the beginning and end of the time step,
                  ! relative to eta_PF, with SAL effects included [H ~> m or kg m-2].

  ! These are always allocated with symmetric memory and wide halos.
                  ! bt_rem_u is between 0 and 1.
end subroutine btstep
module subroutine btstep_timeloop(eta, ubt, vbt, uhbt0, Datu, BTCL_u, vhbt0, Datv, BTCL_v, eta_IC, &
                eta_PF_1, d_eta_PF, eta_src, dyn_coef_eta, uhbtav, vhbtav, u_accel_bt, v_accel_bt, &
                f_4_u, f_4_v, bt_rem_u, bt_rem_v, &
                BT_force_u, BT_force_v, Cor_ref_u, Cor_ref_v, Rayleigh_u, Rayleigh_v, &
                eta_PF, gtot_E, gtot_W, gtot_N, gtot_S, SpV_col_avg, dgeo_de, &
                eta_sum, eta_wtd, ubt_wtd, vbt_wtd, Coru_avg, PFu_avg, LDu_avg, Corv_avg, PFv_avg, &
                LDv_avg, use_BT_cont, interp_eta_PF, find_etaav, dt, dtbt, nstep, nfilter, &
                wt_vel, wt_eta, wt_accel, wt_trans, wt_accel2, ADp, BT_OBC, CS, G, MS, GV, US)

  type(barotropic_CS),    intent(inout) :: CS    !< Barotropic control structure
  type(ocean_grid_type),  intent(inout) :: G     !< The ocean's grid structure (inout to allow for halo updates)
  type(memory_size_type), intent(in)    :: MS    !< A type that describes the memory sizes of
                                                 !! the argument arrays.
  real, dimension(SZIW_(CS),SZJW_(CS)), target, intent(inout) :: &
    eta           !< The barotropic free surface height anomaly or column mass anomaly [H ~> m or kg m-2]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    ubt           !< The zonal barotropic velocity [L T-1 ~> m s-1]
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    vbt           !< The meridional barotropic velocity [L T-1 ~> m s-1]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    uhbt0         !< The difference between the sum of the layer zonal thickness flux and the
                  !! barotropic thickness flux using the same velocity [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    Datu          !< Basin depth at u-velocity grid points times the y-grid spacing [H L ~> m2 or kg m-1]
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: &
    BTCL_u        !< Structure of information used for a dynamic estimate of the face areas at u-points.
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    vhbt0         !< The difference between the sum of the layer meridional thickness flux and the
                  !! barotropic thickness flux using the same velocity [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    Datv          !< Basin depth at v-velocity grid points times the x-grid spacing [H L ~> m2 or kg m-1]
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: &
    BTCL_v        !< Structure of information used for a dynamic estimate of the face areas at v-points
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_IC        !< A local copy of the initial 2-D eta field (eta_in) [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_PF_1      !< The initial value of eta_PF, when interp_eta_PF is true [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    d_eta_PF      !< The change in eta_PF over the barotropic time stepping when
                  !! interp_eta_PF is true [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_src       !< The source of eta per barotropic timestep [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    dyn_coef_eta  !< The coefficient relating the changes in eta to the dynamic surface pressure
                  !! under rigid ice [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: &
    uhbtav        !< the barotropic zonal volume or mass fluxes averaged through the barotropic
                  !! steps [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: &
    vhbtav        !< the barotropic meridional volume or mass fluxes averaged through the barotropic
                  !! steps [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    u_accel_bt    !! The difference between the zonal acceleration from the
                  !< barotropic calculation and BT_force_v [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    v_accel_bt    !< The difference between the meridional acceleration from the
                  !! barotropic calculation and BT_force_v [L T-2 ~> m s-2].
  real, dimension(4,SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    f_4_u         !< The terms giving the contribution to the Coriolis acceleration at a zonal
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at v points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! v-velocities to the southwest, southeast, northwest and northeast.
  real, dimension(4,SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    f_4_v         !< The terms giving the contribution to the Coriolis acceleration at a meridional
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at u points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! u-velocities to the southwest, southeast, northwest and northeast.
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    bt_rem_u      !< The fraction of the barotropic zonal velocity that remains after a time step,
                  !! the rest being lost to bottom drag [nondim].  bt_rem_v is between 0 and 1.
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    bt_rem_v      !< The fraction of the barotropic meridional velocity that remains after a time step,
                  !! the rest being lost to bottom drag [nondim].  bt_rem_v is between 0 and 1.
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    BT_force_u    !< The vertical average of all of the v-accelerations that are
                  !! not explicitly included in the barotropic equation [L T-2 ~> m s-2]
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    BT_force_v    !< The vertical average of all of the v-accelerations that are
                  !! not explicitly included in the barotropic equation [L T-2 ~> m s-2]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    Cor_ref_u     !< The meridional barotropic Coriolis acceleration due
                  !! to the reference velocities [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    Cor_ref_v     !< The meridional barotropic Coriolis acceleration due
                  !! to the reference velocities [L T-2 ~> m s-2].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    Rayleigh_u    !< A Rayleigh drag timescale operating at u-points for drag parameterizations
                  !! that introduced directly into the barotropic solver rather than coming
                  !! in via the visc_rem_u arrays from the layered equations [T-1 ~> s-1]
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    Rayleigh_v    !< A Rayleigh drag timescale operating at v-points for drag parameterizations
                  !! that introduced directly into the barotropic solver rather than coming
                  !! in via the visc_rem_v arrays from the layered equations [T-1 ~> s-1]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(inout) :: &
    eta_PF        !< The 2-D eta field (either SSH anomaly or column mass anomaly) that was used to
                  !! calculate the input pressure gradient accelerations [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_E        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the east of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_W        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the west of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
                  !! (See Hallberg, J Comp Phys 1997 for a discussion of gtot_E and gtot_W.)
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_N        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the north of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_S        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the south of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
                  !! (See Hallberg, J Comp Phys 1997 for a discussion of gtot_E and gtot_W.)
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in) :: &
    SpV_col_avg   !< The column average specific volume [R-1 ~> m3 kg-1]
  real,    intent(in) :: dgeo_de !< The constant of proportionality between geopotential and
                  !! sea surface height [nondim].  It is of order 1, but for stability this
                  !! may be made larger than the physical problem would suggest.
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(out) :: &
    eta_sum       !< eta summed across the timesteps [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(out) :: &
    eta_wtd       !< A weighted estimate used to calculate eta_out [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: &
    ubt_wtd       !< A weighted sum used to find the filtered final ubt [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: &
    vbt_wtd       !< A weighted sum used to find the filtered final vbt [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: &
    Coru_avg      !< The average zonal barotropic Coriolis acceleration [L T-2 ~> m s-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: &
    PFu_avg       !< The average zonal barotropic pressure gradient force [L T-2 ~> m s-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: &
    LDu_avg       !< The average zonal barotropic linear wave drag acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: &
    Corv_avg      !< The average meridional barotropic Coriolis acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: &
    PFv_avg       !< The average meridional barotropic pressure gradient force [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: &
    LDv_avg       !< The average meridional barotropic linear wave drag acceleration [L T-2 ~> m s-2]
  logical, intent(in) :: use_BT_cont  !< If true, use the information in the bt_cont_types to
                  !! calculate the mass transports
  logical, intent(in) :: interp_eta_PF !< If true, interpolate the reference value of eta used
                  !! to calculate the pressure force with time.
  logical, intent(in) :: find_etaav !< If true, diagnose the time mean value of eta
  real,    intent(in) :: dt       !< The time increment to integrate over [T ~> s]
  real,    intent(in) :: dtbt     !< The barotropic time step [T ~> s]
  integer, intent(in) :: nstep    !< The number of barotropic time steps to take to cover the specified time interval
  integer, intent(in) :: nfilter  !< The number of extra barotropic steps to take to allow for time filtering
  real, dimension(nstep+nfilter), intent(in) :: &
    wt_vel        !< The raw or relative weights of each of the barotropic timesteps
                  !! in determining the average velocities [nondim]
  real, dimension(nstep+nfilter), intent(in) :: &
    wt_eta        !< The raw or relative weights of each of the barotropic timesteps
                  !! in determining the average eta [nondim]
  real, dimension(nstep+nfilter+1), intent(in) :: &
    wt_accel      !< The raw or relative weights of each of the barotropic timesteps
                  !! in determining the average accelerations [nondim]
  real, dimension(nstep+nfilter+1), intent(in) :: &
    wt_trans      !< The raw or relative weights of each of the barotropic timesteps
                  !! in determining the average transports [nondim]
  real, dimension(nstep+nfilter+1), intent(in) :: &
    wt_accel2     !< Potentially un-normalized relative weights of each of the
                  !! barotropic timesteps in determining the average accelerations [nondim]
  type(accel_diag_ptrs),    pointer       :: ADp     !< Acceleration diagnostic pointers
  type(BT_OBC_type),        intent(in)    :: BT_OBC  !< A structure with the private barotropic arrays
                                                     !! related to the open boundary conditions,
                                                     !! with time evolving data stored via set_up_BT_OBC
  type(verticalGrid_type),  intent(in)    :: GV      !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)    :: US      !< A dimensional unit scaling type

  ! Local variables
                      ! dtbt_diag = dt/(nstep+nfilter)
                      ! when project_velocity is true [nondim]. For now be_proj is set
                      ! to equal bebt, as they have similar roles and meanings.
                      ! source is all used up by the beginning of the filtering [nondim]
                      ! from the initial condition using the time-integrated barotropic velocity.

end subroutine btstep_timeloop
module subroutine btstep_find_Cor(q, DCor_u, DCor_v, f_4_u, f_4_v, isvf, ievf, jsvf, jevf, CS)
  type(barotropic_CS), intent(inout) :: CS      !< Barotropic control structure
  real, intent(in) :: q(SZIBW_(CS),SZJBW_(CS))  !< A pseudo potential vorticity [T-1 Z-1 ~> s-1 m-1]
                  !! or [T-1 H-1 ~> s-1 m-1 or m2 s-1 kg-1]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    DCor_u        !< An averaged depth or total thickness at u points [Z ~> m] or [H ~> m or kg m-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    DCor_v        !< An averaged depth or total thickness at v points [Z ~> m] or [H ~> m or kg m-2].
  real, dimension(4,SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    f_4_u         !< The terms giving the contribution to the Coriolis acceleration at a zonal
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at v points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! v-velocities to the southwest, southeast, northwest and northeast.
  real, dimension(4,SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    f_4_v         !< The terms giving the contribution to the Coriolis acceleration at a meridional
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at u points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! u-velocities to the southwest, southeast, northwest and northeast.
  integer, intent(in) :: isvf  !< The starting i-index of the largest valid range for tracer points
  integer, intent(in) :: ievf  !< The ending i-index of the largest valid range for tracer points
  integer, intent(in) :: jsvf  !< The starting j-index of the largest valid range for tracer points
  integer, intent(in) :: jevf  !< The ending j-index of the largest valid range for tracer points

  ! real :: C1_3 ! One third [nondim]

end subroutine btstep_find_Cor
module subroutine truncate_velocities(ubt, vbt, dt, G, CS, isv, iev, jsv, jev)
  type(ocean_grid_type), intent(inout) :: G  !< The ocean's grid structure.
  type(barotropic_CS),   intent(inout) :: CS !< Barotropic control structure
  real,    intent(inout) :: ubt(SZIBW_(CS),SZJW_(CS)) !< The zonal barotropic velocity [L T-1 ~> m s-1]
  real,    intent(inout) :: vbt(SZIW_(CS),SZJBW_(CS)) !< The meridional barotropic velocity [L T-1 ~> m s-1]
  real,    intent(in)    :: dt  !< The time increment to integrate over [T ~> s].
  integer, intent(in)    :: isv !< The starting valid tracer array i-index that is being worked on
  integer, intent(in)    :: iev !< The ending valid tracer array i-index that is being worked on
  integer, intent(in)    :: jsv !< The starting valid tracer array j-index that is being worked on
  integer, intent(in)    :: jev !< The ending valid tracer array j-index being that is worked on


end subroutine truncate_velocities
module subroutine btloop_eta_predictor(n, dtbt, ubt, vbt, eta, ubt_int, vbt_int, uhbt, vhbt, uhbt0, vhbt0, &
                        uhbt_int, vhbt_int, BTCL_u, BTCL_v, Datu, Datv, &
                        eta_IC, eta_src, eta_pred, isv, iev, jsv, jev, &
                        integral_BT_cont, use_BT_cont, G, US, CS)
  type(ocean_grid_type), intent(in)  :: G     !< The ocean's grid structure
  type(barotropic_CS),   intent(in)  :: CS    !< Barotropic control structure
  integer,               intent(in)  :: n     !< The current step in loop of timesteps
  real,                  intent(in)  :: dtbt  !< The barotropic time step [T ~> s]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    ubt           !< The zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    vbt           !< The zonal barotropic velocity [L T-1 ~> m s-1].
  real, target, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta           !< The barotropic free surface height anomaly or column mass
                  !! anomaly [H ~> m or kg m-2]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    ubt_int       !< The running time integral of ubt over the time steps [L ~> m].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    vbt_int       !< The running time integral of vbt over the time steps [L ~> m].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    uhbt0         !< The difference between the sum of the layer zonal thickness
                  !! fluxes and the barotropic thickness flux using the same
                  !! velocity [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    vhbt0         !< The difference between the sum of the layer meridional
                  !! thickness fluxes and the barotropic thickness flux using
                  !! the same velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    BTCL_u        !< A repackaged version of the u-point information in BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    BTCL_v        !< A repackaged version of the v-point information in BT_cont.
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    Datu          !< Basin depth at u-velocity grid points times the y-grid
                  !! spacing [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    Datv          !< Basin depth at v-velocity grid points times the x-grid
                  !! spacing [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_IC        !< A local copy of the initial 2-D eta field (eta_in) [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_src       !< The source of eta per barotropic timestep [H ~> m or kg m-2].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    uhbt          !< The zonal barotropic thickness fluxes [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    vhbt          !< The meridional barotropic thickness fluxes [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    uhbt_int      !< The running time integral of uhbt over the time steps [H L2 ~> m3 or kg].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    vhbt_int      !< The running time integral of vhbt over the time steps [H L2 ~> m3 or kg].
  real, target, dimension(SZIW_(CS),SZJW_(CS)), intent(inout) :: &
    eta_pred      !< A predictor value of eta [H ~> m or kg m-2] like eta.
  integer, intent(in)  :: isv         !< The starting i-index of eta_pred to calculate
  integer, intent(in)  :: iev         !< The ending i-index of eta_pred to calculate
  integer, intent(in)  :: jsv         !< The starting j-index of eta_pred to calculate
  integer, intent(in)  :: jev         !< The ending j-index of eta_pred to calculate
  logical, intent(in)  :: integral_BT_cont !< If true, update the barotropic continuity equation directly
                                      !! from the initial condition using the time-integrated barotropic velocity.
  logical, intent(in)  :: use_BT_cont !< If true, use the information in the BT_cont_type to determine
                                      !! barotropic transports as a function of the barotropic velocities.
  type(unit_scale_type), intent(in)  :: US  !< A dimensional unit scaling type


  !$OMP parallel default(shared)
end subroutine btloop_eta_predictor
module subroutine btloop_find_PF(PFu, PFv, isv, iev, jsv, jev, eta_PF_BT, eta_PF, &
                          gtot_N, gtot_S, gtot_E, gtot_W, dgeo_de, find_etaav, &
                          wt_accel2_n, eta_sum, v_first, G, US, CS)
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure.
  type(barotropic_CS),     intent(inout) :: CS    !< Barotropic control structure
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    PFu           !< The anomalous zonal pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    PFv           !< The meridional pressure force acceleration [L T-2 ~> m s-2].
  integer, intent(in)  :: isv         !< The starting i-index of eta being set in ths loop
  integer, intent(in)  :: iev         !< The ending i-index of eta_pred being set in ths loop
  integer, intent(in)  :: jsv         !< The starting j-index of eta_pred being set in ths loop
  integer, intent(in)  :: jev         !< The ending j-index of eta_pred being set in ths loop
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_PF_BT     !< The eta array (either the SSH anomaly or column mass anomaly) that
                  !! determines the barotropic pressure force [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_PF        !< The input 2-D eta field (either SSH anomaly or column mass anomaly)
                  !! that was used to calculate the input pressure gradient
                  !! accelerations [H ~> m or kg m-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_N        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the north of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_S        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the south of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  !! (See Hallberg, J Comp Phys 1997 for a discussion of gtot_E and gtot_W.)
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_E        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the east of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_W        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the west of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  !! (See Hallberg, J Comp Phys 1997 for a discussion of gtot_E and gtot_W.)
  real,    intent(in) :: dgeo_de !< The constant of proportionality between geopotential and
                  !! sea surface height [nondim].  It is of order 1, but for stability this
                  !! may be made larger than the physical  problem would suggest.
  logical, intent(in) :: find_etaav !< If true, diagnose the time mean value of eta
  real,    intent(in) :: wt_accel2_n !< The weighting value of wt_accel2 at step n.
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(inout) :: &
    eta_sum       !< A weighted running sum of eta summed across the timesteps [H ~> m or kg m-2]
  logical, intent(in) :: v_first !< If true, update the v-velocity first with the present loop iteration
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type

  ! Local variables

  ! Ensure that the extra points used for the temporally staggered Coriolis terms are updated.
end subroutine btloop_find_PF
module subroutine btloop_add_dyn_PF(PFu, PFv, eta_pred, eta, dyn_coef_eta, p_surf_dyn, &
                             isv, iev, jsv, jev, v_first, G, US, CS)
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure.
  type(barotropic_CS),     intent(inout) :: CS    !< Barotropic control structure
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    PFu           !< The anomalous zonal pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    PFv           !< The meridional pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta_pred      !< The updated eta field (either SSH anomaly or column mass anomaly) that is
                  !! used to estimate the divergence that is to be damped [H ~> m or kg m-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    eta           !< The previous eta field (either SSH anomaly or column mass anomaly) that is
                  !! used to estimate the divergence that is to be damped [H ~> m or kg m-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    dyn_coef_eta  !< The coefficient relating the changes in eta to the dynamic surface pressure
                  !! under rigid ice [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(inout) :: &
    p_surf_dyn    !< A dynamic surface pressure under rigid ice [L2 T-2 ~> m2 s-2].
  integer, intent(in)  :: isv         !< The starting i-index of eta being set in ths loop
  integer, intent(in)  :: iev         !< The ending i-index of eta_pred being set in ths loop
  integer, intent(in)  :: jsv         !< The starting j-index of eta_pred being set in ths loop
  integer, intent(in)  :: jev         !< The ending j-index of eta_pred being set in ths loop
  logical, intent(in) :: v_first !< If true, update the v-velocity first with the present loop iteration
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type

  ! Local variables

  ! Ensure that the extra points used for the temporally staggered Coriolis terms are updated.
end subroutine btloop_add_dyn_PF
module subroutine btloop_update_v(dtbt, ubt, vbt, v_accel_bt, &
                           Cor_v, PFv, is_v, ie_v, Js_v, Je_v, f_4_v, &
                           bt_rem_v, BT_force_v, Cor_ref_v, Rayleigh_v, &
                           wt_accel_n, G, US, CS, Cor_bracket_bug)
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure.
  type(barotropic_CS),     intent(inout) :: CS    !< Barotropic control structure
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    ubt           !< The zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    vbt           !< The meridional barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    v_accel_bt    !< The difference between the meridional acceleration from the
                  !! barotropic calculation and BT_force_v [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(inout) :: &
    Cor_v         !< The meridional Coriolis acceleration [L T-2 ~> m s-2]
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    PFv           !< The meridional pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(4,SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    f_4_v         !< The terms giving the contribution to the Coriolis acceleration at a meridional
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at u points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! u-velocities to the southwest, southeast, northwest and northeast.
  integer, intent(in)  :: is_v !< The starting i-index of the range of v-point values to calculate
  integer, intent(in)  :: ie_v !< The ending i-index of the range of v-point values to calculate
  integer, intent(in)  :: Js_v !< The starting j-index of the range of v-point values to calculate
  integer, intent(in)  :: Je_v !< The ending j-index of the range of v-point values to calculate
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    bt_rem_v      !< The fraction of the barotropic meridional velocity that
                  !! remains after a time step, the rest being lost to bottom
                  !! drag [nondim].  bt_rem_v is between 0 and 1.
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    BT_force_v    !< The vertical average of all of the v-accelerations that are
                  !! not explicitly included in the barotropic equation [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    Cor_ref_v     !< The meridional barotropic Coriolis acceleration due
                  !! to the reference velocities [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    Rayleigh_v    !< A Rayleigh drag timescale operating at v-points for drag parameterizations
                  !! that introduced directly into the barotropic solver rather than coming
                  !! in via the visc_rem_v arrays from the layered equations [T-1 ~> s-1]
  real,    intent(in) :: wt_accel_n  !< The raw or relative weights of each of the barotropic timesteps
                  !! in determining the average accelerations [nondim]
  real,    intent(in) :: dtbt !< The barotropic time step [T ~> s].
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  logical, optional, intent(in) :: Cor_bracket_bug !< If present and true, use an order of operations that is
                  !! not bitwise rotationally symmetric in the meridional Coriolis term

  ! Local variables

end subroutine btloop_update_v
module subroutine btloop_update_u(dtbt, ubt, vbt, u_accel_bt, &
                           Cor_u, PFu, Is_u, Ie_u, js_u, je_u, f_4_u, &
                           bt_rem_u, BT_force_u, Cor_ref_u, Rayleigh_u, &
                           wt_accel_n, G, US, CS)
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure.
  type(barotropic_CS),     intent(inout) :: CS    !< Barotropic control structure
  real,    intent(in) :: dtbt     !< The barotropic time step [T ~> s].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    ubt           !< The zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    vbt           !< The meridional barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    u_accel_bt    !! The difference between the zonal acceleration from the
                  !< barotropic calculation and BT_force_v [L T-2 ~> m s-2].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(inout) :: &
    Cor_u         !< The anomalous zonal Coriolis acceleration [L T-2 ~> m s-2]
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    PFu           !< The anomalous zonal pressure force acceleration [L T-2 ~> m s-2].
  integer, intent(in)  :: Is_u !< The starting i-index of the range of u-point values to calculate
  integer, intent(in)  :: Ie_u !< The ending i-index of the range of u-point values to calculate
  integer, intent(in)  :: js_u !< The starting j-index of the range of u-point values to calculate
  integer, intent(in)  :: je_u !< The ending j-index of the range of u-point values to calculate
  real, dimension(4,SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    f_4_u         !< The terms giving the contribution to the Coriolis acceleration at a zonal
                  !! velocity point from the neighboring meridional velocity anomalies [T-1 ~> s-1].
                  !! These are the products of thicknesses at v points and appropriately staggered
                  !! averaged pseudo potential vorticities, but with sufficiently smooth topography
                  !! they are approximately f / 4.  The 4 values on the innermost loop are for
                  !! v-velocities to the southwest, southeast, northwest and northeast.
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    bt_rem_u      !< The fraction of the barotropic meridional velocity that
                  !! remains after a time step, the rest being lost to bottom
                  !! drag [nondim].  bt_rem_v is between 0 and 1.
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    BT_force_u    !< The vertical average of all of the v-accelerations that are
                  !! not explicitly included in the barotropic equation [L T-2 ~> m s-2].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    Cor_ref_u     !< The meridional barotropic Coriolis acceleration due
                  !! to the reference velocities [L T-2 ~> m s-2].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    Rayleigh_u    !< A Rayleigh drag timescale operating at u-points for drag parameterizations
                  !! that introduced directly into the barotropic solver rather than coming
                  !! in via the visc_rem_u arrays from the layered equations [T-1 ~> s-1].
  real,    intent(in) :: wt_accel_n  !< The raw or relative weights of each of the barotropic timesteps
                                  !! in determining the average accelerations [nondim]
  type(unit_scale_type),   intent(in)  :: US      !< A dimensional unit scaling type

  ! Local variables

  !$OMP do schedule(static)
end subroutine btloop_update_u
module subroutine btstep_ubt_from_layer(U_in, V_in, wt_u, wt_v, ubt, vbt,  G, GV, CS)
  type(verticalGrid_type), intent(in)  :: GV      !< The ocean's vertical grid structure.
  type(barotropic_CS),     intent(inout) :: CS    !< Barotropic control structure
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure.
  real, intent(in)  :: U_in(SZIB_(G),SZJ_(G),SZK_(GV)) !< The initial (3-D) zonal velocity [L T-1 ~> m s-1]
  real, intent(in)  :: V_in(SZI_(G),SZJB_(G),SZK_(GV)) !< The initial (3-D) meridional velocity [L T-1 ~> m s-1]
  real, intent(in)  :: wt_u(SZIB_(G),SZJ_(G),SZK_(GV)) !< The normalized weights to be used in calculating
                                                  !! zonal barotropic velocities, possibly with sums
                                                  !! less than one due to viscous losses [nondim]
  real, intent(in)  :: wt_v(SZI_(G),SZJB_(G),SZK_(GV)) !< The normalized weights to be used in calculating
                                                  !! meridional barotropic velocities, possibly with
                                                  !! sums less than one due to viscous losses [nondim]
  real, intent(out) :: ubt(SZIBW_(CS),SZJW_(CS))  !< The zonal barotropic velocity [L T-1 ~> m s-1]
  real, intent(out) :: vbt(SZIW_(CS),SZJBW_(CS))  !< The meridional barotropic velocity [L T-1 ~> m s-1]

  ! Local variables

end subroutine btstep_ubt_from_layer
module subroutine btstep_layer_accel(dt, u_accel_bt, v_accel_bt, pbce, gtot_E, gtot_W, gtot_N, gtot_S, &
                              e_anom, G, GV, CS, accel_layer_u, accel_layer_v)
  type(barotropic_CS),      intent(inout) :: CS !< Barotropic control structure
  type(ocean_grid_type),    intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV !< The ocean's vertical grid structure.
  real, intent(in)  :: dt      !< The time increment to integrate over [T ~> s].
  real, dimension(SZIBW_(CS),SZJW_(CS)), intent(in) :: &
    u_accel_bt  !< The difference between the zonal acceleration from the
                !! barotropic calculation and BT_force_u [L T-2 ~> m s-2].
  real, dimension(SZIW_(CS),SZJBW_(CS)), intent(in) :: &
    v_accel_bt  !< The difference between the meridional acceleration from the
                !! barotropic calculation and BT_force_v [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: pbce !< The baroclinic pressure anomaly in each layer
                                                         !! due to free surface height anomalies
                                                         !! [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_E        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the east of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_W        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the west of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_N        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the north of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZIW_(CS),SZJW_(CS)), intent(in) :: &
    gtot_S        !< The effective total reduced gravity used to relate free surface height
                  !! deviations to pressure forces (including GFS and baroclinic contributions)
                  !! in the barotropic momentum equations half a grid-point to the south of a
                  !! thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  !! (See Hallberg, J Comp Phys 1997 for a discussion of gtot_E, etc.)
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_anom        !< The anomaly in the sea surface height or column mass
                  !! averaged between the beginning and end of the time step,
                  !! relative to eta_PF, with SAL effects included [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: accel_layer_u !< The zonal acceleration of each layer due
                                                         !! to the barotropic calculation [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: accel_layer_v !< The meridional acceleration of each layer
                                                         !! due to the barotropic calculation [L T-2 ~> m s-2].

  ! Local variables

end subroutine btstep_layer_accel
module subroutine set_dtbt(G, GV, US, CS, pbce, gtot_est, BT_cont, eta, SSH_add)
  type(ocean_grid_type),        intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),      intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(barotropic_CS),          intent(inout) :: CS   !< Barotropic control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                      optional, intent(in)    :: pbce !< The baroclinic pressure anomaly in each layer due to free
                                                      !! surface height anomalies [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real,               optional, intent(in)    :: gtot_est !< An estimate of the total gravitational acceleration
                                                      !! [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  type(BT_cont_type), optional, pointer       :: BT_cont  !< A structure with elements that describe the effective open
                                                      !! face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZJ_(G)), &
                      optional, intent(in)    :: eta  !< The barotropic free surface height anomaly or  column mass
                                                      !! anomaly [H ~> m or kg m-2].
  real,               optional, intent(in)    :: SSH_add !< An additional contribution to SSH to provide a margin of
                                                      !! error when calculating the external wave speed [Z ~> m].

  ! Local variables
                  ! from the thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  ! (See Hallberg, J Comp Phys 1997 for a discussion.)
                  ! spacing [H L ~> m2 or kg m-1].
                  ! spacing [H L ~> m2 or kg m-1].
                  ! of the reference geopotential with the sea surface height [nondim].
                  ! This is typically ~0.09 or less.
                  ! sea surface height [nondim].  It is a nondimensional number of
                  ! order 1.  For stability, this may be made larger
                  ! than physical problem would suggest.
                  ! when calculating the external wave speed [Z ~> m].
                      ! timesteps [T2 ~> s2]
                      ! barotropic time step [T-2 ~> s-2].


end subroutine set_dtbt
module subroutine apply_u_velocity_OBCs(ubt, uhbt, ubt_trans, eta, SpV_avg, ubt_old, BT_OBC, G, MS, &
                               GV, US, CS, halo, dtbt, bebt, use_BT_cont, integral_BT_cont, dt_elapsed, &
                              Datu, BTCL_u, uhbt0, ubt_int, ubt_int_prev, uhbt_int, uhbt_int_prev)
  type(ocean_grid_type),                 intent(in)    :: G       !< The ocean's grid structure.
  type(memory_size_type),                intent(in)    :: MS      !< A type that describes the memory sizes of
                                                                  !! the argument arrays.
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt     !< the zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: uhbt    !< the zonal barotropic transport
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt_trans !< The zonal barotropic velocity used in
                                                                  !! transport [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta     !< The barotropic free surface height anomaly or
                                                                  !! column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: SpV_avg !< The column average specific volume [R-1 ~> m3 kg-1]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: ubt_old !< The starting value of ubt in a barotropic
                                                                  !! step [L T-1 ~> m s-1].
  type(BT_OBC_type),                     intent(in)    :: BT_OBC  !< A structure with the private barotropic arrays
                                                                  !! related to the open boundary conditions,
                                                                  !! set by set_up_BT_OBC.
  type(verticalGrid_type),               intent(in)    :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),                 intent(in)    :: US      !< A dimensional unit scaling type
  type(barotropic_CS),                   intent(in)    :: CS      !< Barotropic control structure
  integer,                               intent(in)    :: halo    !< The extra halo size to use here.
  real,                                  intent(in)    :: dtbt    !< The time step [T ~> s].
  real,                                  intent(in)    :: bebt    !< The fractional weighting of the future velocity
                                                                  !! in determining the transport [nondim]
  logical,                               intent(in)    :: use_BT_cont !< If true, use the BT_cont_types to calculate
                                                                  !! transports.
  logical,                               intent(in)    :: integral_BT_cont !< If true, update the barotropic continuity
                                                                  !! equation directly from the initial condition
                                                                  !! using the time-integrated barotropic velocity.
  real,                                  intent(in)    :: dt_elapsed !< The amount of time in the barotropic stepping
                                                                  !! that will have elapsed [T ~> s].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: Datu    !< A fixed estimate of the face areas at u points
                                                                  !! [H L ~> m2 or kg m-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u !< Structure of information used
                                                                  !! for a dynamic estimate of the face areas at
                                                                  !! u-points.
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: uhbt0   !< A correction to the zonal transport so that
                                                                  !! the barotropic functions agree with the sum
                                                                  !! of the layer transports
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt_int !< The time-integrated zonal barotropic
                                                                  !! velocity after this update [L T-1 ~> m s-1]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: ubt_int_prev  !< The time-integrated zonal barotropic
                                                                  !! velocity before this update [L T-1 ~> m s-1]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: uhbt_int !< The time-integrated zonal barotropic transport
                                                                  !! after this update [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: uhbt_int_prev !< The time-integrated zonal barotropic
                                                                  !! transport before this update
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]

  ! Local variables

end subroutine apply_u_velocity_OBCs
module subroutine apply_v_velocity_OBCs(vbt, vhbt, vbt_trans, eta, SpV_avg, vbt_old, BT_OBC, &
                               G, MS, GV, US, CS, halo, dtbt, bebt, use_BT_cont, integral_BT_cont, dt_elapsed, &
                               Datv, BTCL_v, vhbt0, vbt_int, vbt_int_prev, vhbt_int, vhbt_int_prev)
  type(ocean_grid_type),                 intent(in)    :: G       !< The ocean's grid structure.
  type(memory_size_type),                intent(in)    :: MS      !< A type that describes the memory sizes of
                                                                  !! the argument arrays.
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt     !< The meridional barotropic velocity
                                                                  !! [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vhbt    !< the meridional barotropic transport
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt_trans !< the meridional BT velocity used in
                                                                  !! transports [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta     !< The barotropic free surface height anomaly or
                                                                  !! column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: SpV_avg !< The column average specific volume [R-1 ~> m3 kg-1]
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vbt_old !< The starting value of vbt in a barotropic
                                                                  !! step [L T-1 ~> m s-1].
  type(BT_OBC_type),                     intent(in)    :: BT_OBC  !< A structure with the private barotropic arrays
                                                                  !! related to the open boundary conditions,
                                                                  !! set by set_up_BT_OBC.
  type(verticalGrid_type),               intent(in)    :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),                 intent(in)    :: US      !< A dimensional unit scaling type
  type(barotropic_CS),                   intent(in)    :: CS      !< Barotropic control structure
  integer,                               intent(in)    :: halo    !< The extra halo size to use here.
  real,                                  intent(in)    :: dtbt    !< The time step [T ~> s].
  real,                                  intent(in)    :: bebt    !< The fractional weighting of the future velocity
                                                                  !! in determining the transport [nondim]
  logical,                               intent(in)    :: use_BT_cont !< If true, use the BT_cont_types to calculate
                                                                  !! transports.
  logical,                               intent(in)    :: integral_BT_cont !< If true, update the barotropic continuity
                                                                  !! equation directly from the initial condition
                                                                  !! using the time-integrated barotropic velocity.
  real,                                  intent(in)    :: dt_elapsed !< The amount of time in the barotropic stepping
                                                                  !! that will have elapsed [T ~> s].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: Datv    !< A fixed estimate of the face areas at v points
                                                                  !! [H L ~> m2 or kg m-1].
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v !< Structure of information used
                                                                  !! for a dynamic estimate of the face areas at
                                                                  !! v-points.
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vhbt0   !< A correction to the meridional transport so that
                                                                  !! the barotropic functions agree with the sum
                                                                  !! of the layer transports
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt_int !< The time-integrated meridional barotropic
                                                                  !! velocity after this update [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vbt_int_prev !< The time-integrated meridional barotropic
                                                                  !! velocity before this update [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vhbt_int !< The time-integrated meridional barotropic
                                                                  !! transport after this update
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vhbt_int_prev !< The time-integrated meridional barotropic
                                                                  !! transport before this update
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]

  ! Local variables

end subroutine apply_v_velocity_OBCs
module subroutine initialize_BT_OBC(OBC, BT_OBC, G, CS)
  type(ocean_OBC_type), target,          intent(inout) :: OBC    !< An associated pointer to an OBC type.
  type(BT_OBC_type),                     intent(inout) :: BT_OBC !< A structure with the private barotropic arrays
                                                                 !! related to the open boundary conditions,
                                                                 !! set by set_up_BT_OBC.
  type(ocean_grid_type),                 intent(inout) :: G      !< The ocean's grid structure.
  type(barotropic_CS),                   intent(inout) :: CS     !< Barotropic control structure

  ! Local variables
                 ! converted to real numbers to work with the MOM6 halo update code [nondim]
                 ! converted to real numbers to work with the MOM6 halo update code [nondim]
                            ! with a southern OBC in a northern halo.

end subroutine initialize_BT_OBC
module subroutine set_up_BT_OBC(OBC, eta, SpV_avg, BT_OBC, BT_Domain, G, GV, US, CS, MS, halo, use_BT_cont, &
                         integral_BT_cont, dt_baroclinic, Datu, Datv, BTCL_u, BTCL_v, dgeo_de)
  type(ocean_OBC_type), target,          intent(inout) :: OBC    !< An associated pointer to an OBC type.
  type(memory_size_type),                intent(in)    :: MS     !< A type that describes the memory sizes of the
                                                                 !! argument arrays.
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta    !< The barotropic free surface height anomaly or
                                                                 !! column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: SpV_avg !< The column average specific volume [R-1 ~> m3 kg-1]
  type(BT_OBC_type),                     intent(inout) :: BT_OBC !< A structure with the private barotropic arrays
                                                                 !! related to the open boundary conditions,
                                                                 !! set by set_up_BT_OBC.
  type(MOM_domain_type),                 intent(inout) :: BT_Domain !< MOM_domain_type associated with wide arrays
  type(ocean_grid_type),                 intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),                 intent(in)    :: US     !< A dimensional unit scaling type
  type(barotropic_CS),                   intent(inout) :: CS     !< Barotropic control structure
  integer,                               intent(in)    :: halo   !< The extra halo size to use here.
  logical,                               intent(in)    :: use_BT_cont !< If true, use the BT_cont_types to calculate
                                                                 !! transports.
  logical,                               intent(in)    :: integral_BT_cont !< If true, update the barotropic continuity
                                                                 !! equation directly from the initial condition
                                                                 !! using the time-integrated barotropic velocity.
  real,                                  intent(in)    :: dt_baroclinic !< The baroclinic timestep for this cycle of
                                                                 !! updates to the barotropic solver [T ~> s]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: Datu   !< A fixed estimate of the face areas at u points
                                                                 !! [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: Datv   !< A fixed estimate of the face areas at v points
                                                                 !! [H L ~> m2 or kg m-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u !< Structure of information used
                                                                 !! for a dynamic estimate of the face areas at
                                                                 !! u-points.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v !< Structure of information used
                                                                 !! for a dynamic estimate of the face areas at
                                                                 !! v-points.
  real,                                  intent(in)    :: dgeo_de  !< The constant of proportionality between
                                                                 !! geopotential and sea surface height [nondim].
  ! Local variables

end subroutine set_up_BT_OBC
module subroutine destroy_BT_OBC(BT_OBC)
  type(BT_OBC_type), intent(inout) :: BT_OBC !< A structure with the private barotropic arrays
                                             !! related to the open boundary conditions,
                                             !! set by set_up_BT_OBC.

end subroutine destroy_BT_OBC
module subroutine btcalc(h, G, GV, CS, h_u, h_v, may_use_default, OBC)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(barotropic_CS),     intent(inout) :: CS   !< Barotropic control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: h_u  !< The specified effective thicknesses at u-points,
                                                 !! perhaps scaled down to account for viscosity and
                                                 !! fractional open areas [H ~> m or kg m-2].  These
                                                 !! are used here as non-normalized weights for each
                                                 !! layer that are converted the normalized weights
                                                 !! for determining the barotropic accelerations.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: h_v  !< The specified effective thicknesses at v-points,
                                                 !! perhaps scaled down to account for viscosity and
                                                 !! fractional open areas [H ~> m or kg m-2].  These
                                                 !! are used here as non-normalized weights for each
                                                 !! layer that are converted the normalized weights
                                                 !! for determining the barotropic accelerations.
  logical,       optional, intent(in)    :: may_use_default !< An optional logical argument
                                                 !! to indicate that the default velocity point
                                                 !! thicknesses may be used for this particular
                                                 !! calculation, even though the setting of
                                                 !! CS%hvel_scheme would usually require that h_u
                                                 !! and h_v be passed in.
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundary control structure.

  ! Local variables
                               ! in roundoff and can be neglected [H ~> m or kg m-2].
                               ! The harmonic mean uses a weight of (1 - wt_arith).
                               ! around a u-point (positive upward) [H ~> m or kg m-2]
                               ! around a v-point (positive upward) [H ~> m or kg m-2]


end subroutine btcalc
module function find_uhbt(u, BTC) result(uhbt)
  real, intent(in) :: u    !< The local zonal velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.

  real :: uhbt !< The zonal barotropic transport [L2 H T-1 ~> m3 s-1] or time integrated transport [L2 H ~> m3]

end function find_uhbt
module function find_duhbt_du(u, BTC) result(duhbt_du)
  real, intent(in) :: u    !< The local zonal velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: duhbt_du !< The zonal barotropic face area [L H ~> m2 or kg m-1]

end function find_duhbt_du
module function uhbt_to_ubt(uhbt, BTC) result(ubt)
  real, intent(in) :: uhbt                      !< The barotropic zonal transport that should be inverted for,
                                                !! [H L2 T-1 ~> m3 s-1 or kg s-1] or the time-integrated
                                                !! transport [H L2 ~> m3 or kg].
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that allow the
                                                !! barotropic transports to be calculated consistently with the
                                                !! layers' continuity equations.  The dimensions of some
                                                !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: ubt                                   !< The result - The velocity that gives uhbt transport [L T-1 ~> m s-1]
                                                !! or the time-integrated velocity [L ~> m].

  ! Local variables
                                 ! or [H L2 ~> m3 or kg].
                                 ! maximum increase of vs2, both [nondim].

  ! Find the value of ubt that gives uhbt.
end function uhbt_to_ubt
module function find_vhbt(v, BTC) result(vhbt)
  real, intent(in) :: v    !< The local meridional velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: vhbt !< The meridional barotropic transport [L2 H T-1 ~> m3 s-1] or time integrated transport [L2 H ~> m3]

end function find_vhbt
module function find_dvhbt_dv(v, BTC) result(dvhbt_dv)
  real, intent(in) :: v    !< The local meridional velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: dvhbt_dv !< The meridional barotropic face area [L H ~> m2 or kg m-1]

end function find_dvhbt_dv
module function vhbt_to_vbt(vhbt, BTC) result(vbt)
  real, intent(in) :: vhbt                      !< The barotropic meridional transport that should be
                                                !! inverted for [H L2 T-1 ~> m3 s-1 or kg s-1] or the
                                                !! time-integrated transport [H L2 ~> m3 or kg].
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that allow the
                                                !! barotropic transports to be calculated consistently
                                                !! with the layers' continuity equations.  The dimensions of some
                                                !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: vbt                                   !< The result - The velocity that gives vhbt transport [L T-1 ~> m s-1]
                                                !! or the time-integrated velocity [L ~> m].

  ! Local variables
                                 ! or [H L2 ~> m3 or kg].
                                 ! maximum increase of vs2, both [nondim].

  ! Find the value of vbt that gives vhbt.
end function vhbt_to_vbt
module subroutine set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, US, MS, BT_Domain, halo, dt_baroclinic)
  type(BT_cont_type),     intent(inout) :: BT_cont    !< The BT_cont_type input to the barotropic solver
  type(memory_size_type), intent(in)    :: MS         !< A type that describes the memory sizes of
                                                      !! the argument arrays
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(out) :: BTCL_u       !< A structure with the u information from BT_cont
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(out) :: BTCL_v       !< A structure with the v information from BT_cont
  type(ocean_grid_type),  intent(in)    :: G          !< The ocean's grid structure
  type(unit_scale_type),  intent(in)    :: US         !< A dimensional unit scaling type
  type(MOM_domain_type),  intent(inout) :: BT_Domain  !< The domain to use for updating the halos
                                                      !! of wide arrays
  integer,                intent(in)    :: halo       !< The extra halo size to use here
  real,         optional, intent(in)    :: dt_baroclinic !< The baroclinic time step [T ~> s], which
                                                      !! is provided if INTEGRAL_BT_CONTINUITY is true.

  ! Local variables

end subroutine set_local_BT_cont_types
module subroutine adjust_local_BT_cont_types(ubt, uhbt, vbt, vhbt, BTCL_u, BTCL_v, &
                                      G, US, MS, halo, dt_baroclinic)
  type(memory_size_type), intent(in)  :: MS   !< A type that describes the memory sizes of the argument arrays.
  real, dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(in)  :: ubt  !< The linearization zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(in)  :: uhbt !< The linearization zonal barotropic transport
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(in)  :: vbt  !< The linearization meridional barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(in)  :: vhbt !< The linearization meridional barotropic transport
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(out) :: BTCL_u !< A structure with the u information from BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(out) :: BTCL_v !< A structure with the v information from BT_cont.
  type(ocean_grid_type),  intent(in)  :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)  :: US   !< A dimensional unit scaling type
  integer,                intent(in)  :: halo !< The extra halo size to use here.
  real,         optional, intent(in)  :: dt_baroclinic !< The baroclinic time step [T ~> s], which is
                                                       !! provided if INTEGRAL_BT_CONTINUITY is true.

  ! Local variables

end subroutine adjust_local_BT_cont_types
module subroutine BT_cont_to_face_areas(BT_cont, Datu, Datv, G, US, MS, halo)
  type(BT_cont_type),     intent(inout) :: BT_cont    !< The BT_cont_type input to the
                                                      !! barotropic solver.
  type(memory_size_type), intent(in)    :: MS         !< A type that describes the memory
                                                      !! sizes of the argument arrays.
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), &
                          intent(out)   :: Datu       !< The effective zonal face area [H L ~> m2 or kg m-1].
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), &
                          intent(out)   :: Datv       !< The effective meridional face area [H L ~> m2 or kg m-1].
  type(ocean_grid_type),  intent(in)    :: G          !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US         !< A dimensional unit scaling type
  integer,      optional, intent(in)    :: halo       !< The extra halo size to use here.

  ! Local variables
end subroutine BT_cont_to_face_areas
module subroutine swap(a,b)
  real, intent(inout) :: a !< The first variable to be swapped [arbitrary units]
  real, intent(inout) :: b !< The second variable to be swapped [arbitrary units]
end subroutine swap
module subroutine find_face_areas(Datu, Datv, G, GV, US, CS, MS, halo, eta, add_max)
  type(memory_size_type),  intent(in) :: MS    !< A type that describes the memory sizes of the argument arrays.
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), &
                           intent(out) :: Datu !< The open zonal face area [H L ~> m2 or kg m-1].
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), &
                           intent(out) :: Datv !< The open meridional face area [H L ~> m2 or kg m-1].
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(barotropic_CS),     intent(in)  :: CS   !< Barotropic control structure
  integer,                 intent(in)  :: halo !< The halo size to use, default = 1.
  real, dimension(MS%isdw:MS%iedw,MS%jsdw:MS%jedw), &
                 optional, intent(in)  :: eta  !< The barotropic free surface height anomaly
                                               !! or column mass anomaly [H ~> m or kg m-2].
  real,          optional, intent(in)  :: add_max !< A value to add to the maximum depth (used
                                               !! to overestimate the external wave speed) [Z ~> m].

  ! Local variables
end subroutine find_face_areas
module subroutine bt_mass_source(h, eta, set_cor, G, GV, CS)
  type(ocean_grid_type),              intent(in) :: G        !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in) :: GV       !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),   intent(in) :: eta      !< The free surface height that is to be
                                                             !! corrected [H ~> m or kg m-2].
  logical,                            intent(in) :: set_cor  !< A flag to indicate whether to set the corrective
                                                             !! fluxes (and update the slowly varying part of eta_cor)
                                                             !! (.true.) or whether to incrementally update the
                                                             !! corrective fluxes.
  type(barotropic_CS),                intent(inout) :: CS    !< Barotropic control structure

  ! Local variables
                              ! the sum of the layer thicknesses [H ~> m or kg m-2].
                              ! thicknesses [H ~> m or kg m-2].

end subroutine bt_mass_source
module subroutine barotropic_init(u, v, h, Time, G, GV, US, param_file, diag, CS, &
                           restart_CS, calc_dtbt, BT_cont, OBC, SAL_CSp, HA_CSp)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(barotropic_CS),     intent(inout) :: CS   !< Barotropic control structure
  type(MOM_restart_CS),    intent(in)    :: restart_CS !< MOM restart control structure
  logical,                 intent(out)   :: calc_dtbt  !< If true, the barotropic time step must
                                                 !! be recalculated before stepping.
  type(BT_cont_type),      pointer       :: BT_cont    !< A structure with elements that describe the
                                                 !! effective open face areas as a function of
                                                 !! barotropic flow.
  type(ocean_OBC_type),    pointer       :: OBC  !< The open boundary condition structure.
  type(SAL_CS), target,    optional      :: SAL_CSp  !< A pointer to the control structure of the
                                                 !! SAL module.
  type(harmonic_analysis_CS), target, optional :: HA_CSp !< A pointer to the control structure of the
                                                 !! harmonic analysis module

  ! This include declares and sets the variable "version".
  ! Local variables
                        ! upper-bound estimate for pbce.
                        ! in calculating the safe external wave speed [Z ~> m].
                          ! piston velocities [nondim].
                                       ! drag piston velocity.
                                       ! name in wave_drag_file.
                                       ! name in wave_drag_file.
                                       ! name in wave_drag_file.
                      ! geopotential with the sea surface height when scalar SAL are enabled [nondim].
                      ! This is typically ~0.09 or less.
                      ! in roundoff and can be neglected [H L2 ~> m3 or kg]
                                        ! that acts on the barotropic flow [H T-1 ~> m s-1 or kg m-2 s-1].

                          ! so that diagnosed barotropic pressure gradient forces are zero at
                          ! land, coastal or OBC points.
                          ! recreate the bugs, or if false bugs are only used if actively selected.
end subroutine barotropic_init
module subroutine barotropic_get_tav(CS, ubtav, vbtav, G, US)
  type(barotropic_CS),               intent(in)    :: CS    !< Barotropic control structure
  type(ocean_grid_type),             intent(in)    :: G     !< Grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: ubtav !< Zonal barotropic velocity averaged
                                                            !! over a baroclinic timestep [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vbtav !< Meridional barotropic velocity averaged
                                                            !! over a baroclinic timestep [L T-1 ~> m s-1]
  type(unit_scale_type),             intent(in)    :: US    !< A dimensional unit scaling type
  ! Local variables

end subroutine barotropic_get_tav
module subroutine barotropic_end(CS)
  type(barotropic_CS), intent(inout) :: CS  !< Control structure to clear out.

end subroutine barotropic_end
module subroutine register_barotropic_restarts(HI, GV, US, param_file, CS, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(barotropic_CS),     intent(inout) :: CS      !< Barotropic control structure
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables

end subroutine register_barotropic_restarts
  end interface

end module MOM_barotropic
