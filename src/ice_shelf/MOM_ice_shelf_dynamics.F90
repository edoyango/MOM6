! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements a crude placeholder for a later implementation of full
!! ice shelf dynamics.
module MOM_ice_shelf_dynamics

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_IS_diag_mediator, only : post_data=>post_IS_data
use MOM_IS_diag_mediator, only : register_diag_field=>register_MOM_IS_diag_field, safe_alloc_ptr
!use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_init, set_IS_diag_mediator_grid
use MOM_IS_diag_mediator, only : diag_ctrl, time_type, enable_averages, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE, AGRID, CORNER, CENTER
use MOM_domains, only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_io, only : file_exists, slasher, MOM_read_data
use MOM_io, only : open_ASCII_file, get_filename_appendix
use MOM_io, only : APPEND_FILE, WRITEONLY_FILE
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_time_manager, only : time_type, get_time, set_time, time_type_to_real, operator(>)
use MOM_time_manager,  only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,  only : operator(/=), operator(<=), operator(>=), operator(<)
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state
use MOM_coms, only : reproducing_sum, max_across_PEs, min_across_PEs
use MOM_checksums, only : hchksum, qchksum
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_channel,initialize_ice_flow_from_file
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_from_file,initialize_ice_C_basal_friction
use MOM_ice_shelf_initialize, only : initialize_ice_AGlen
implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf, IS_dynamics_post_data
public ice_time_step_CFL, ice_shelf_dyn_end, change_in_draft, write_ice_shelf_energy
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask, volume_above_floatation
public masked_var_grounded

! SSA inner solver flags
integer, parameter :: INNER_CG = 1       !< Conjugate gradient (default)
integer, parameter :: INNER_MINRES = 2   !< MINRES
integer, parameter :: INNER_CR = 3       !< Conjugate residual

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the ice shelf dynamics.
type, public :: ice_shelf_dyn_CS ; private
  real, pointer, dimension(:,:) :: u_shelf => NULL() !< the zonal velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: v_shelf => NULL() !< the meridional velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: taudx_shelf => NULL() !< the zonal driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: taudy_shelf => NULL() !< the meridional driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: sx_shelf => NULL() !< the zonal surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: sy_shelf => NULL() !< the meridional surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: u_face_mask => NULL() !< mask for velocity boundary conditions on the C-grid
                                       !! u-face - this is because the FEM cares about FACES THAT GET INTEGRATED OVER,
                                       !! not vertices. Will represent boundary conditions on computational boundary
                                       !! (or permanent boundary between fast-moving and near-stagnant ice
                                       !! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition,
                                       !! 3=inhomogeneous Dirichlet boundary for u and v, 4=flux boundary: at these
                                       !! faces a flux will be specified which will override velocities; a homogeneous
                                       !! velocity condition will be specified (this seems to give the solver less
                                       !! difficulty)  5=inhomogenous Dirichlet boundary for u only. 6=inhomogenous
                                       !! Dirichlet boundary for v only
  real, pointer, dimension(:,:) :: v_face_mask => NULL()  !< A mask for velocity boundary conditions on the C-grid
                                       !! v-face, with valued defined similarly to u_face_mask, but 5 is Dirichlet for v
                                       !! and 6 is Dirichlet for u
  real, pointer, dimension(:,:) :: u_face_mask_bdry => NULL() !< A duplicate copy of u_face_mask?
  real, pointer, dimension(:,:) :: v_face_mask_bdry => NULL() !< A duplicate copy of v_face_mask?
  real, pointer, dimension(:,:) :: u_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary u-faces (where u_face_mask=4) [Z L T-1 ~> m2 s-1]
  real, pointer, dimension(:,:) :: v_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary v-faces (where v_face_mask=4) [Z L T-1 ~> m2 s-1]??
   ! needed where u_face_mask is equal to 4, similarly for v_face_mask
  real, pointer, dimension(:,:) :: umask => NULL()      !< u-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: vmask => NULL()      !< v-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: calve_mask => NULL() !< a mask to prevent the ice shelf front from
                                          !! advancing past its initial position (but it may retreat)
  real, pointer, dimension(:,:) :: t_shelf => NULL() !< Vertically integrated temperature in the ice shelf/stream,
                                                     !! on corner-points (B grid) [C ~> degC]
  real, pointer, dimension(:,:) :: tmask => NULL()   !< A mask on tracer points that is 1 where there is ice.
  real, pointer, dimension(:,:,:) :: ice_visc => NULL() !< Area and depth-integrated Glen's law ice viscosity
                                                        !!  (Pa m3 s) in [R L4 Z T-1 ~> kg m2 s-1].
                                                        !!  at either 1 (cell-centered) or 4 quadrature points per cell
  real, pointer, dimension(:,:,:) :: newton_visc_factor => NULL() !< Newton tangent stiffness coefficient:
                                                      !!  (1/n_glen - 1)/2 * ice_visc / eps_e2 at each
                                                      !!  viscosity quadrature point [R L4 Z T ~> kg m2 s]
  real, pointer, dimension(:,:,:) :: newton_str_ux => NULL() !< Longitudinal x-strain-rate ux at each viscosity
                                                      !!  quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:,:) :: newton_str_vy => NULL() !< Longitudinal y-strain-rate vy at each viscosity
                                                      !!  quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:,:) :: newton_str_sh => NULL() !< Engineering shear strain-rate uy+vx at each
                                                      !!  viscosity quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:) :: AGlen_visc => NULL() !< Ice-stiffness parameter in Glen's law ice viscosity,
                                                      !! often in [Pa-3 s-1] if n_Glen is 3.
  real, pointer, dimension(:,:) :: u_bdry_val => NULL() !< The zonal ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: v_bdry_val => NULL() !< The meridional ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: h_bdry_val => NULL() !< The ice thickness at inflowing boundaries [Z ~> m].
  real, pointer, dimension(:,:) :: t_bdry_val => NULL() !< The ice temperature at inflowing boundaries [C ~> degC].

  real, pointer, dimension(:,:) :: bed_elev => NULL()  !< The bed elevation used for ice dynamics [Z ~> m],
                                                       !! relative to mean sea-level.  This is
                                                       !! the same as G%bathyT+Z_ref, when below sea-level.
                                                       !! Sign convention: positive below sea-level, negative above.

  real, pointer, dimension(:,:) :: C_basal_friction => NULL()!< Coefficient in sliding law tau_b = C u^(n_basal_fric),
                               !! units of [R L Z T-2 (s m-1)^(n_basal_fric) ~> Pa (s m-1)^(n_basal_fric)]
  real, pointer, dimension(:,:) :: coef_prefactor => NULL() !< Pre-computed area*C_basal_friction*L_T_to_m_s for
                               !! basal friction quadrature evaluation [R L2 Z T-1 ~> kg s-1].
  real, pointer, dimension(:,:) :: fB_elem => NULL()        !< Pre-computed element-level Coulomb fB parameter
                               !! [(T L-1)^CF_PostPeak]; 0 for Weertman.
                               !! Updated each outer iteration by calc_shelf_basal_prefactors.
  real :: alpha_coulomb = 1.0  !< Coulomb prefactor (CF_PostPeak-1)^(CF_PostPeak-1)/CF_PostPeak^CF_PostPeak [nondim]
  real :: coulomb_pp_n         !< CF_PostPeak/n_basal_fric [nondim]
  real, pointer, dimension(:,:) :: OD_rt => NULL()         !< A running total for calculating OD_av [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac_rt => NULL() !< A running total for calculating ground_frac.
  real, pointer, dimension(:,:) :: OD_av => NULL()         !< The time average open ocean depth [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac => NULL()   !< Fraction of the time a cell is "exposed", i.e. the column
                               !! thickness is below a threshold and interacting with the rock [nondim].  When this
                               !! is 1, the ice-shelf is grounded
  real, pointer, dimension(:,:) :: float_cond => NULL()   !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, pointer, dimension(:,:,:,:) :: Phi => NULL() !< The gradients of bilinear basis elements at Gaussian
                                                !! 4 quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, pointer, dimension(:,:,:) :: PhiC => NULL()  !< The gradients of bilinear basis elements at 1 cell-centered
                                                !! quadrature point per cell [L-1 ~> m-1].
  real, pointer, dimension(:,:,:) :: Jac => NULL()   !< Jacobian determinant |J_q| = a_q*d_q of the element
                                                !! mapping at each of the 4 Gaussian quadrature points [L2 ~> m2].
                                                !! Equal to G%areaT only for rectangular elements; differs when
                                                !! opposite cell edges have unequal lengths (non-rectangular quads).
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub => NULL() !< Quadrature structure weights at subgridscale
                                                !!  locations for finite element calculations [nondim]
  integer :: OD_rt_counter = 0 !< A counter of the number of contributions to OD_rt.

  real :: velocity_update_time_step !< The time interval over which to update the ice shelf velocity
                    !! using the nonlinear elliptic equation, or 0 to update every timestep [T ~> s].
                    ! DNGoldberg thinks this should be done no more often than about once a day
                    ! (maybe longer) because it will depend on ocean values  that are averaged over
                    ! this time interval, and solving for the equilibrated flow will begin to lose
                    ! meaning if it is done too frequently.
  real :: elapsed_velocity_time  !< The elapsed time since the ice velocities were last updated [T ~> s].

  real :: g_Earth      !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: density_ice  !< A typical density of ice [R ~> kg m-3].
  real :: Cp_ice       !< The heat capacity of fresh ice [Q C-1 ~> J kg-1 degC-1].

  logical :: advect_shelf !< If true (default), advect ice shelf and evolve thickness
  logical :: reentrant_x !< If true, the domain is zonally reentrant
  logical :: reentrant_y !< If true, the domain is meridionally reentrant
  logical :: alternate_first_direction_IS !< If true, alternate whether the x- or y-direction
                                          !! updates occur first in directionally split parts of the calculation.
  integer :: first_direction_IS !< An integer that indicates which direction is
                                !! to be updated first in directionally split
                                !! parts of the ice sheet calculation (e.g. advection).
  real    :: first_dir_restart_IS = -1.0 !< A real copy of CS%first_direction_IS for use in restart files
  integer :: visc_qps !< The number of quadrature points per cell (1 or 4) on which to calculate ice viscosity.
  character(len=40) :: ice_viscosity_compute !< Specifies whether the ice viscosity is computed internally
                                   !! according to Glen's flow law; is constant (for debugging purposes)
                                   !! or using observed strain rates and read from a file
  logical :: shelf_top_slope_bugs !< If true, use directionally inconsistent estimates of the grid
                            !! spacing when calculating the ice shelf surface slope, and underestimate
                            !! slopes near the edge of the ice shelf by a factor of 2.
  logical :: GL_regularize  !< Specifies whether to regularize the floatation condition
                            !! at the grounding line as in Goldberg Holland Schoof 2009
  integer :: n_sub_regularize
                            !< partition of cell over which to integrate for
                            !! interpolated grounding line the (rectangular) is
                            !! divided into nxn equally-sized rectangles, over which
                            !!  basal contribution is integrated (iterative quadrature)
  logical :: GL_couple      !< whether to let the floatation condition be
                            !! determined by ocean column thickness means update_OD_ffrac
                            !! will be called (note: GL_regularize and GL_couple
                            !! should be exclusive)

  real    :: CFL_factor     !< A factor used to limit subcycled advective timestep in uncoupled runs
                            !! i.e. dt <= CFL_factor * min(dx / u) [nondim]

  real :: min_h_shelf !< The minimum ice thickness used during ice dynamics [Z ~> m].
  real :: min_basal_traction !< The minimum basal traction for grounded ice (Pa m-1 s) [R Z T-1 ~> kg m-2 s-1]
  real :: max_surface_slope !< The maximum allowed ice-sheet surface slope (to ignore, set to zero) [nondim]
  real :: min_ice_visc !< The minimum allowed Glen's law ice viscosity (Pa s), in [R L2 T-1 ~> kg m-1 s-1].

  real :: n_glen            !< Nonlinearity exponent in Glen's Law [nondim]
  real :: eps_glen_min      !< Min. strain rate to avoid infinite Glen's law viscosity, [T-1 ~> s-1].
  real :: n_basal_fric      !< Exponent in sliding law tau_b = C u^(m_slide) [nondim]
  logical :: CoulombFriction !< Use Coulomb friction law (Schoof 2005, Gagliardini et al 2007)
  real :: CF_MinN           !< Minimum Coulomb friction effective pressure [R Z L T-2 ~> Pa]
  real :: CF_PostPeak       !< Coulomb friction post peak exponent [nondim]
  real :: CF_Max            !< Coulomb friction maximum coefficient [nondim]
  real :: density_ocean_avg !< A typical ocean density [R ~> kg m-3].  This does not affect ocean
                            !! circulation or thermodynamics.  It is used to estimate the
                            !! gravitational driving force at the shelf front (until we think of
                            !! a better way to do it, but any difference will be negligible).
  real :: rhoi_rhow         !< The density of ice divided by a typical water density [nondim]
  real :: rhow_rhoi         !< A typical water density divided by the density of ice [nondim]
  real :: thresh_float_col_depth !< The water column depth over which the shelf if considered to be floating
  logical :: moving_shelf_front  !< Specify whether to advance shelf front (and calve).
  logical :: calve_to_mask       !< If true, calve off the ice shelf when it passes the edge of a mask.
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].
  real :: T_shelf_missing   !< An ice shelf temperature to use where there is no ice shelf [C ~> degC]
  real :: cg_tolerance !< For Picard iterations, the tolerance in the CG solver, relative to initial residual, that
                       !! determines when to stop the conjugate gradient iterations [nondim].
  real :: cg_newton_tolerance !< For inexact Newton iterations, the initial tolerance in the CG solver, relative to
                              !!  initial residual, that determines when to stop the CG iterations [nondim].
  real :: cg_tol_current !< Working CG tolerance for the current inner solve [nondim].
  real :: nonlinear_tolerance !< The fractional nonlinear tolerance, relative to the initial error,
                              !! that sets when to stop the iterative velocity solver [nondim]
  real :: newton_after_tolerance !< The fractional nonlinear tolerance, relative to the initial error, at
                                 !! which to switch from Picard to Newton iterations in the velocity solver
                                 !! If set to <= 0, no Picard [nondim]
  type(group_pass_type) :: pass_visc_and_newton !< Handle for Newton-and-viscosity-related group passes
  type(group_pass_type) :: pass_newton !< Handle for Newton-related group passes
  logical :: newton_adapt_cg_tol !< Use an adaptive CG tolerance during Newton iterations
  real :: ew_gamma !< Gamma in Eisenstat-Walker adaptive Newton tolerance [nondim].
  real :: ew_alpha !< Alpha in Eisenstat-Walker adaptive Newton tolerance [nondim].
  integer :: ew_safety !< Safeguard Eisenstat-Walker using:
                    !!(0) no safeguard, (1) EW choice 2 threshold or (2) PETSc option 3 (Chacon 2008)
  real :: ew_1_thres !< Threshold for Eisenstat-Walker version 1 [nondim]
  real :: ew_eta_max !< Maximum allowed Eisenstat-Walker eta [nondim]
  integer :: cg_max_iterations !< The maximum number of iterations that can be used in the CG solver
  integer :: nonlin_solve_err_mode  !< 1: exit based on nonlin residual | F | / | F_0 | where | | is infty-norm
                    !! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol) where | | is infty-norm
                    !! 3: exit based on change of solution norm 2*abs(|u|-|u_last|)/(|u|+|u_last|) where | | is L2-norm
                    !! 4: exit based on nonlin residual  | F | / | F_0 | where | | is L2-norm
                    !! 5: exit based on relative residual | F | / | tau | where | | is L2-norm
  logical :: ssa_add_rel_resid !< Nonlinear error in velocity solve will also depend on the
                                   !! L2 residual norm relative to RHS norm
  real :: rr_nonlinear_tolerance !< If ssa_add_rel_resid, the additional nonlin tolerance in the iterative
                    !! velocity solve used for the relative residual [nondim]
  ! for write_ice_shelf_energy
  type(time_type) :: energysavedays            !< The interval between writing the energies
                                               !! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  !< The starting interval for computing a geometric
                                               !! progression of time deltas between calls to
                                               !! write_energy. This interval will increase by a factor of 2.
                                               !! after each call to write_energy.
  logical         :: energysave_geometric      !< Logical to control whether calls to write_energy should
                                               !! follow a geometric progression
  type(time_type) :: write_energy_time         !< The next time to write to the energy file.
  type(time_type) :: geometric_end_time        !< Time at which to stop the geometric progression
                                               !! of calls to write_energy and revert to the standard
                                               !! energysavedays interval
  real    :: timeunit           !< The length of the units for the time axis and certain input parameters
                                !! including ENERGYSAVEDAYS [s].
  type(time_type) :: Start_time !< The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  integer :: prev_IS_energy_calls = 0 !< The number of times write_ice_shelf_energy has been called.
  integer :: IS_fileenergy_ascii   !< The unit number of the ascii version of the energy file.
  character(len=200) :: IS_energyfile  !< The name of the ice sheet energy file with path.

  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
  logical :: doing_newton = .false. !< If true, the outer iteration is using Newton (tangent) linearization
                                    !! instead of Picard (secant) linearization for the ice viscosity
  integer :: inner_solver !< The inner linear solver: INNER_CG (1),INNER_MINRES (2), or INNER_CR (3)
  logical :: cg_halo_shrink = .true. !< If true, CG uses halo-shrinking to defer pass_vector calls;
                                     !! if false, uses fixed CG_action range with 1 pass_vector per iteration
  logical :: module_is_initialized = .false. !< True if this module has been initialized.

  !>@{ Diagnostic handles
  integer :: id_u_shelf = -1, id_v_shelf = -1, id_shelf_speed, id_t_shelf = -1, &
             id_taudx_shelf = -1, id_taudy_shelf = -1, id_taud_shelf = -1, id_bed_elev = -1, &
             id_ground_frac = -1, id_col_thick = -1, id_OD_av = -1, id_float_cond = -1, &
             id_u_mask = -1, id_v_mask = -1, id_ufb_mask =-1, id_vfb_mask = -1, id_t_mask = -1, &
             id_sx_shelf = -1, id_sy_shelf = -1, id_surf_slope_mag_shelf, &
             id_duHdx = -1, id_dvHdy = -1, id_fluxdiv = -1, &
             id_strainrate_xx = -1, id_strainrate_yy = -1, id_strainrate_xy = -1, &
             id_pstrainrate_1 = -1, id_pstrainrate_2, &
             id_devstress_xx = -1, id_devstress_yy = -1, id_devstress_xy = -1, &
             id_pdevstress_1 = -1, id_pdevstress_2 = -1

  !>@}
  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !>@{ Diagnostic handles for debugging
  integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1, &
             id_visc_shelf = -1, id_taub = -1
  !>@}
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to control diagnostic output.

end type ice_shelf_dyn_CS

!> A container for loop bounds
type :: loop_bounds_type ; private
  integer :: ish !< Starting i-index of the computational domain [nondim]
  integer :: ieh !< Ending i-index of the computational domain [nondim]
  integer :: jsh !< Starting j-index of the computational domain [nondim]
  integer :: jeh !< Ending j-index of the computational domain [nondim]
end type loop_bounds_type


  interface
module function slope_limiter(num, denom)
  real, intent(in)    :: num   !< The numerator of the ratio used in the Van Leer slope limiter
  real, intent(in)    :: denom !< The denominator of the ratio used in the Van Leer slope limiter
  real :: slope_limiter ! The slope limiter value, between 0 and 2 [nondim].

end function slope_limiter
module function quad_area (X, Y)
  real, dimension(4), intent(in) :: X !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4), intent(in) :: Y !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real :: quad_area ! Computed area [L2 ~> m2]

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

end function quad_area
module subroutine register_ice_shelf_dyn_restarts(G, US, param_file, CS, restart_CS)
  type(ocean_grid_type),  intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_dyn_CS), pointer       :: CS !< A pointer to the ice shelf dynamics control structure
  type(MOM_restart_CS),   intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables

end subroutine register_ice_shelf_dyn_restarts
module subroutine initialize_ice_shelf_dyn(param_file, Time, ISS, CS, G, US, diag, new_sim, Cp_ice, &
                                    Input_start_time, directory, solo_ice_sheet_in)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(time_type),         intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_state),   intent(in)    :: ISS  !< A structure with elements that describe
                                                 !! the ice-shelf state
  type(ice_shelf_dyn_CS),  pointer       :: CS   !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type),   intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),   intent(in)    :: US   !< A structure containing unit conversion factors
  type(diag_ctrl), target, intent(in)    :: diag !< A structure that is used to regulate the diagnostic output.
  logical,                 intent(in)    :: new_sim !< If true this is a new simulation, otherwise
                                                 !! has been started from a restart file.
  real,                    intent(in)    :: Cp_ice !< Heat capacity of ice [Q C-1 ~> J kg-1 degC-1]
  type(time_type),         intent(in)    :: Input_start_time !< The start time of the simulation.
  character(len=*),        intent(in)    :: directory  !< The directory where the ice sheet energy file goes.
  logical,       optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                 !! a solo ice-sheet driver.

  ! Local variables
                          ! in through open boundaries [C ~> degC]
  !This include declares and sets the variable "version".
                          ! recreate the bugs, or if false bugs are only used if actively selected.

end subroutine initialize_ice_shelf_dyn
module subroutine initialize_diagnostic_fields(CS, ISS, G, US, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(time_type),        intent(in)    :: Time !< The current model time

!
end subroutine initialize_diagnostic_fields
module function ice_time_step_CFL(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  real :: ice_time_step_CFL !< The maximum permitted timestep based on the ice velocities [T ~> s].


end function ice_time_step_CFL
module subroutine update_ice_shelf(CS, ISS, G, US, time_step, Time, calve_ice_shelf_bergs, &
                            ocean_mass, coupled_grounding, must_update_vel)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< To convert ice flux through front
                                                                 !! to bergs
  real, dimension(SZDI_(G),SZDJ_(G)), &
                optional, intent(in)    :: ocean_mass !< If present this is the mass per unit area
                                              !! of the ocean [R Z ~> kg m-2].
  logical,      optional, intent(in)    :: coupled_grounding !< If true, the grounding line is
                                              !! determined by coupled ice-ocean dynamics
  logical,      optional, intent(in)    :: must_update_vel !< Always update the ice velocities if true.

end subroutine update_ice_shelf
module subroutine volume_above_floatation(CS, G, ISS, vaf, hemisphere)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  type(ice_shelf_state),  intent(in) :: ISS !< A structure with elements that describe
                                            !! the ice-shelf state
  real, intent(out) :: vaf !< area integrated volume above floatation [Z L2 ~> m3]
  integer, optional, intent(in) :: hemisphere !< 0 for Antarctica only, 1 for Greenland only. Otherwise, all ice sheets

end subroutine volume_above_floatation
module subroutine masked_var_grounded(G,CS,var,varout)
  type(ocean_grid_type), intent(in) :: G !< The grid structure used by the ice shelf.
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var !< variable in
  real, dimension(SZI_(G),SZJ_(G)), intent(out)  :: varout !<variable out
end subroutine masked_var_grounded
module subroutine IS_dynamics_post_data(time_step, Time, CS, ISS, G)
  real :: time_step !< Length of time for post data averaging [T ~> s].
  type(time_type),        intent(in)    :: Time !< The current model time
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
                                                  !! [R L2 Z T-1 ~> Pa s m]
                                                  !! [R L T-1 ~> Pa s m-1]


end subroutine IS_dynamics_post_data
module subroutine ice_visc_diag(CS,G,ice_visc)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(out)  :: ice_visc !< area-averaged vertically integrated ice viscosity
                                                               !! [R L2 Z T-1 ~> Pa s m]

end subroutine ice_visc_diag
module subroutine write_ice_shelf_energy(CS, G, US, mass, area, day, time_step)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: mass !< The mass per unit area of the ice shelf
                                                !! or sheet [R Z ~> kg m-2]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(in)    :: area !< The ice shelf or ice sheet area [L2 ~> m2]
  type(time_type),         intent(in)    :: day !< The current model time.
  type(time_type),  optional, intent(in) :: time_step !< The current time step
  ! Local variables

  ! write_energy_time is the next integral multiple of energysavedays.
end subroutine write_ice_shelf_energy
module subroutine ice_shelf_advect(CS, ISS, G, time_step, Time, calve_ice_shelf_bergs)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< If true, track ice shelf flux through a
                                               !! static ice shelf, so that it can be converted into icebergs

! 3/8/11 DNG
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!    ADDITIONALLY, it will update the volume of ice in partially-filled cells, and update
!        hmask accordingly
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells


end subroutine ice_shelf_advect
module subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, iters, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time

  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudx !< Driving x-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudy !< Driving y-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: u_bdry_cont ! Boundary u-stress contribution [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: v_bdry_cont ! Boundary v-stress contribution [R L3 Z T-2 ~> kg m s-2]
                                                ! the grounding line (float_cond=1) or not (float_cond=0)

end subroutine ice_shelf_solve_outer
module subroutine ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, float_cond, &
                                 hmask, conv_flag, iters, time, Phi, Phisub)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudx !< The x-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudy  !< The y-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal (corner) points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< If GL_regularize=true, indicates cells containing
                                                      !! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                                 !! partly or fully covered by an ice-shelf
  integer,                intent(out)   :: conv_flag !< A flag indicating whether (1) or not (0) the
                                                     !! iterations have converged to the specified tolerance
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                               !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                                  !! locations for finite element calculations [nondim]

                           ! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]

end subroutine ice_shelf_solve_inner
module subroutine ice_shelf_solve_inner_CG(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                     IDIAGu, IDIAGv, H_node, float_cond, hmask, &
                                     rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                     Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< Grounding line indicator [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays

                                                   !! at quadrature points [L T-1 ~> m s-1]
                                                   !! at quadrature points [L T-1 ~> m s-1]
                                                  ! [kg m2 s-3]
                                                       ! sum_vec_3d(:,:,1) [kg m2 s-3]
                                                       ! sum_vec_3d(:,:,2) [kg2 m2 s-4]
                         ! sv3dsums(1) [kg m2 s-3]
                         ! sv3dsums(2) [kg2 m2 s-4]
                         ! iteration, scaled by resid_scale [kg m2 s-3]
                          ! [T4 kg2 m2 R-2 Z-2 L-6 s-4 ~> 1]

end subroutine ice_shelf_solve_inner_CG
module subroutine ice_shelf_solve_inner_MINRES(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                         IDIAGu, IDIAGv, H_node, float_cond, hmask, &
                                         rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                         Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< Grounding line indicator [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays

                                                   !! at quadrature points [L T-1 ~> m s-1]
                                                   !! at quadrature points [L T-1 ~> m s-1]
                                                     ! [kg m2 s-3] before normalization;
                                                     ! [nondim] inside loop (after Lanczos normalization)
                         ! initial value [kg^1/2 m s^-3/2], then [nondim] after iter 1
                         ! [kg m2 s-3] before normalization, [nondim] inside loop

end subroutine ice_shelf_solve_inner_MINRES
module subroutine ice_shelf_solve_inner_CR(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                     IDIAGu, IDIAGv, H_node, float_cond, hmask, &
                                     rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                     Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< Grounding line indicator [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays

                                                   !! at quadrature points [L T-1 ~> m s-1]
                                                   !! at quadrature points [L T-1 ~> m s-1]
                                ! sum_vec_3d(:,:,1): r^2 [kg2 m2 s-4] or z·q [kg m2 s-3] (context-dependent)
                                ! sum_vec_3d(:,:,2): z·w or q·(M^-1 q) [kg m2 s-3]
                          ! sv3dsums(1): r^2 or z·q [kg2 m2 s-4 or kg m2 s-3] (context-dependent)
                          ! sv3dsums(2): z·w or q·M^-1 q [kg m2 s-3]

end subroutine ice_shelf_solve_inner_CR
module subroutine ice_shelf_advect_thickness_x(CS, G, LB, time_step, hmask, h0, h_after_uflux, uh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses after
                                              !! the zonal mass fluxes [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization



!  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

end subroutine ice_shelf_advect_thickness_x
module subroutine ice_shelf_advect_thickness_y(CS, G, LB, time_step, hmask, h0, h_after_vflux, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses after
                                              !! the meridional mass fluxes [Z ~> m].
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization



end subroutine ice_shelf_advect_thickness_y
module subroutine shelf_advance_front(CS, ISS, G, hmask, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                           !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! in this subroutine we go through the computational cells only and, if they are empty or partial cells,
  ! we find the reference thickness and update the shelf mass and partial area fraction and the hmask if necessary

  ! if any cells go from partial to complete, we then must set the thickness, update hmask accordingly,
  ! and divide the overflow across the adjacent EMPTY (not partly-covered) cells.
  ! (it is highly unlikely there will not be any; in which case this will need to be rethought.)

  ! most likely there will only be one "overflow". If not, though, a pass_var of all relevant variables
  ! is done; there will therefore be a loop which, in practice, will hopefully not have to go through
  ! many iterations

  ! when 3d advected scalars are introduced, they will be impacted by what is done here

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from eastern neighbor:  flux_enter(:,:,1)
  !   from western neighbor:  flux_enter(:,:,2)
  !   from southern neighbor: flux_enter(:,:,3)
  !   from northern neighbor: flux_enter(:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !


                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].
                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].

end subroutine shelf_advance_front
module subroutine ice_shelf_min_thickness_calve(G, h_shelf, area_shelf_h, hmask, thickness_calve, halo)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                  intent(in)    :: thickness_calve !< The thickness at which to trigger calving [Z ~> m].
  integer,     optional, intent(in)    :: halo  !< The number of halo points to use.  If not present,
                                                !! work on the entire data domain.

end subroutine ice_shelf_min_thickness_calve
module subroutine calve_to_mask(G, h_shelf, area_shelf_h, hmask, calve_mask)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)    :: calve_mask !< A mask that indicates where the ice
                                                             !! shelf can exist, and where it will calve.


end subroutine calve_to_mask
module subroutine calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(in)   :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: OD  !< ocean floor depth at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]


! driving stress!

! ! taudx and taudy will hold driving stress in the x- and y- directions when done.
!    they will sit on the BGrid, and so their size depends on whether the grid is symmetric
!
! Since this is a finite element solve, they will actually have the form \int \Phi_i rho g h \nabla s
!
! OD -this is important and we do not yet know where (in MOM) it will come from. It represents
!     "average" ocean depth -- and is needed to find surface elevation
!    (it is assumed that base_ice = bed + OD)


end subroutine calc_shelf_driving_stress
module subroutine CG_action(CS, uret, vret, u_shlf, v_shlf, Phi, Phisub, umask, vmask, hmask, H_node, &
                     ice_visc, float_cond, bathyT, u_curr, v_curr, G, US, is, ie, js, je, dens_ratio, use_newton_in)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: uret !< The retarding stresses working at u-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: vret !< The retarding stresses working at v-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                         intent(in)   :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                         intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                         intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                               !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: float_cond !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points
                                                 !! relative to sea-level [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]

  real,                  intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater, nondimensional
  type(unit_scale_type), intent(in)    :: US  !< A structure containing unit conversion factors
  integer,               intent(in)    :: is  !< The starting i-index to work on
  integer,               intent(in)    :: ie  !< The ending i-index to work on
  integer,               intent(in)    :: js  !< The starting j-index to work on
  integer,               intent(in)    :: je  !< The ending j-index to work on
  logical, optional,     intent(in)    :: use_newton_in !< If present, overrides CS%doing_newton for Newton correction

! the linear action of the matrix on (u,v) with bilinear finite elements
! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
! but this may change pursuant to conversations with others
!
! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!     in order to make less frequent halo updates

! the linear action of the matrix on (u,v) with bilinear finite elements
! Phi has the form
! Phi(k,q,i,j) - applies to cell i,j

    !  3 - 4
    !  |   |
    !  1 - 2

! Phi(2*k-1,q,i,j) gives d(Phi_k)/dx at quadrature point q
! Phi(2*k,q,i,j) gives d(Phi_k)/dy at quadrature point q
! Phi_k is equal to 1 at vertex k, and 0 at vertex l /= k, and bilinear


end subroutine CG_action
module subroutine CG_action_subgrid_basal(CS, G, US, Phisub, H, U_curr, V_curr, U_delta, V_delta, &
                                   bathyT, dens_ratio, i_elem, j_elem, fB_e, use_newton, Ucontr, Vcontr, &
                                   dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2),   intent(in) :: H       !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: U_delta !< Search direction δu at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_delta !< Search direction δv at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  logical,                intent(in) :: use_newton !< If true, include Newton basal drag correction
  real, dimension(2,2),   intent(out) :: Ucontr !< Nodal u-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2),   intent(out) :: Vcontr !< Nodal v-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real,                   intent(in) :: dxCv_S !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in) :: dxCv_N !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_W !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_E !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in) :: IareaT !< The inverse of the cell area at the tracer point [L-2 ~> m-2]

                                                                                !! at each sub-cell
                                                ! accumulated then pair-summed for rotation invariance

end subroutine CG_action_subgrid_basal
module subroutine compute_basal_coef(unorm2_qp, coef_prefactor, min_trac_area, fB_e, &
    n_basal_fric, CoulombFriction, CF_PostPeak, L_T_to_m_s, use_newton, &
    basal_coef, drag_newt)
  real,    intent(in)  :: unorm2_qp      !< Regularized |u^k|^2 > 0 at quadrature point [L2 T-2 ~> m2 s-2]
  real,    intent(in)  :: coef_prefactor !< Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real,    intent(in)  :: min_trac_area  !< Pre-computed min_basal_traction * areaT floor [R L2 Z T-1 ~> kg s-1]
  real,    intent(in)  :: fB_e           !< Element-level Coulomb fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  real,    intent(in)  :: n_basal_fric   !< Friction sliding exponent m [nondim]
  logical, intent(in)  :: CoulombFriction !< True if using Coulomb friction
  real,    intent(in)  :: CF_PostPeak    !< Coulomb post-peak exponent q [nondim]
  real,    intent(in)  :: L_T_to_m_s    !< Unit conversion factor from internal [L T-1] to [m s-1]
  logical, intent(in)  :: use_newton     !< If true, evaluate drag_newt; otherwise set to 0
  real,    intent(out) :: basal_coef     !< Picard friction coefficient at quadrature point [R L2 Z T-1 ~> kg s-1]
  real,    intent(out) :: drag_newt      !< Newton drag coefficient [R Z T ~> kg m-2 s]; 0 without Newton


end subroutine compute_basal_coef
module subroutine sum_square_matrix(sum_out, mat_in, n)
  integer, intent(in) :: n !< The length and width of each matrix in mat_in
  real, dimension(n,n), intent(in) :: mat_in !< The n x n matrix whose elements will be summed
  real, intent(out) :: sum_out !< The sum of the elements of matrix mat_in

end subroutine sum_square_matrix
module subroutine matrix_diagonal(CS, G, US, float_cond, H_node, ice_visc, u_curr, v_curr, &
                           hmask, dens_ratio, Phi, Phisub, u_diagonal, v_diagonal)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal
                                                 !! (corner) points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                          intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                                !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                   intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater [nondim]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_diagonal !< The diagonal elements of the u-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_diagonal  !< The diagonal elements of the v-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]


! returns the diagonal entries of the matrix for a Jacobi preconditioning


end subroutine matrix_diagonal
module subroutine CG_diagonal_subgrid_basal(CS, G, US, Phisub, H_node, U_curr, V_curr, &
                                     bathyT, dens_ratio, i_elem, j_elem, fB_e, u_diag, v_diag, &
                                     dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub  !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2),   intent(in) :: H_node  !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  real, dimension(2,2),   intent(out) :: u_diag !< Nodal u-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real, dimension(2,2),   intent(out) :: v_diag !< Nodal v-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real,                   intent(in)  :: dxCv_S !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in)  :: dxCv_N !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in)  :: dyCu_W !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in)  :: dyCu_E !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in)  :: IareaT !< The inverse of the cell area at the tracer point [L-2 ~> m-2]

                                                         ! pair-summed for rotation invariance

end subroutine CG_diagonal_subgrid_basal
module subroutine IS_dynamics_post_data_2(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.

  !Allocate the gradient basis functions for 1 cell-centered quadrature point per cell
end subroutine IS_dynamics_post_data_2
module subroutine calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: u_shlf !< The zonal ice shelf velocity [L T-1 ~> m s-1].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: v_shlf !< The meridional ice shelf velocity [L T-1 ~> m s-1].

! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for bilinear FEM solve


! this may be subject to change later... to make it "hybrid"
!  real, dimension(SZDIB_(G),SZDJB_(G)) ::  eII, ux, uy, vx, vy

end subroutine calc_shelf_visc
module subroutine calc_shelf_basal_prefactors(CS, ISS, G, US)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< Ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< Ice shelf state (hmask, h_shelf)
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure
  type(unit_scale_type),  intent(in)    :: US  !< Unit conversion factors


end subroutine calc_shelf_basal_prefactors
module subroutine calc_shelf_taub(CS, ISS, G, basal_tr)
  type(ice_shelf_dyn_CS), intent(in)  :: CS  !< Ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)  :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type),  intent(in)  :: G   !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(out) :: basal_tr !< Area-averaged basal traction [R L T-1 ~> Pa s m-1]


end subroutine calc_shelf_taub
module subroutine update_OD_ffrac(CS, G, US, ocean_mass, find_avg)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)     :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: ocean_mass !< The mass per unit area of the ocean [R Z ~> kg m-2].
  logical,                intent(in)    :: find_avg !< If true, find the average of OD and ffrac, and
                                              !! reset the underlying running sums to 0.


end subroutine update_OD_ffrac
module subroutine update_OD_ffrac_uncoupled(CS, G, h_shelf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf !< the thickness of the ice shelf [Z ~> m].


end subroutine update_OD_ffrac_uncoupled
module subroutine change_in_draft(CS, G, h_shelf0, h_shelf1, ddraft)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf0 !< the previous thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf1 !< the current thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout)    :: ddraft !< the change in shelf draft thickness

end subroutine change_in_draft
module subroutine bilinear_shape_functions (X, Y, Phi, area)
  real, dimension(4),   intent(in)    :: X   !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4),   intent(in)    :: Y   !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(8,4), intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real,                 intent(out)   :: area !< The quadrilateral cell area [L2 ~> m2].

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

! this subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell. values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?
! ... will all cells have the same shape and dimension?


end subroutine bilinear_shape_functions
module subroutine bilinear_shape_fn_grid(G, i, j, Phi, Jac)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8,4),  intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(4), optional, intent(out) :: Jac !< Jacobian determinant |J_q| = a_q*d_q at each
                                              !! Gaussian quadrature point [L2 ~> m2].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?

  ! Mirror lookups: xquad_m(qp) == 1 - xquad(qp), yquad_m(qp) == 1 - yquad(qp) mathematically,
  ! but each mirror entry is the stored value at the x- or y-mirrored quadrature point. This
  ! ensures rotation-paired QPs read bit-identical operand values.

end subroutine bilinear_shape_fn_grid
module subroutine bilinear_shape_fn_grid_1qp(G, i, j, Phi)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8),    intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! a cell-cented point of gaussian quadrature. (in 1D: .5 for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1) gives d(Phi_i)/dx at the quadrature point
! Phi(2*i) gives d(Phi_i)/dy at the quadrature point
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear


    ! d(x)/d(x*)
end subroutine bilinear_shape_fn_grid_1qp
module subroutine bilinear_shape_functions_subgrid(Phisub, nsub)
  integer, intent(in)    :: nsub   !< The number of subgridscale quadrature locations in each direction
  real, dimension(2,2,nsub,nsub,2,2), &
           intent(inout) :: Phisub !< Quadrature structure weights at subgridscale
                                   !! locations for finite element calculations [nondim]

  ! this subroutine is a helper for interpolation of floatation condition
  ! for the purposes of evaluating the terms \int (u,v) \phi_i dx dy in a cell that is
  !     in partial floatation
  ! the array Phisub contains the values of \phi_i (where i is a node of the cell)
  !     at quad point j
  ! i think this general approach may not work for nonrectangular elements...
  !

  ! Phisub(q1,q2,i,j,k,l)
  !  q1: quad point x-index
  !  q2: quad point y-index
  !  i: subgrid index in x-direction
  !  j: subgrid index in y-direction
  !  k: basis function x-index
  !  l: basis function y-index

  ! e.g. k=1,l=1 => node 1
  !      q1=2,q2=1 => quad point 2

    !  3 - 4
    !  |   |
    !  1 - 2

  ! Mirror-symmetric per-direction node weights: a_left == 1-x_global, a_right == x_global
  ! mathematically, but constructed so that a_right(qx,i) is computed by exactly the same
  ! operand sequence as a_left(3-qx, nsub+1-i). This guarantees bit-exact rotation symmetry.

end subroutine bilinear_shape_functions_subgrid
module subroutine update_velocity_masks(CS, G, hmask, umask, vmask, u_face_mask, v_face_mask)
  type(ice_shelf_dyn_CS),intent(in)    :: CS !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: u_face_mask !< A coded mask for velocities at the C-grid u-face
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: v_face_mask !< A coded mask for velocities at the C-grid v-face
  ! sets masks for velocity solve
  ! ignores the fact that their might be ice-free cells - this only considers the computational boundary

  ! !!!IMPORTANT!!! relies on thickness mask - assumed that this is called after hmask has been updated & halo-updated


end subroutine update_velocity_masks
module subroutine interpolate_H_to_B(G, h_shelf, hmask, H_node, min_h_shelf)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: h_shelf !< The ice shelf thickness at tracer points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, intent(in) :: min_h_shelf !< The minimum ice thickness used during ice dynamics [Z ~> m].


end subroutine interpolate_H_to_B
module subroutine ice_shelf_dyn_end(CS)
  type(ice_shelf_dyn_CS), pointer   :: CS !< A pointer to the ice shelf dynamics control structure

end subroutine ice_shelf_dyn_end
module subroutine ice_shelf_temp(CS, ISS, G, US, time_step, melt_rate, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: melt_rate !< basal melt rate [R Z T-1 ~> kg m-2 s-1]
  type(time_type),        intent(in)    :: Time !< The current model time

!    This subroutine takes the velocity (on the Bgrid) and timesteps
!      (HT)_t = - div (uHT) + (adot Tsurf -bdot Tbot) once and then calculates T=HT/H
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells



  ! For now adot and Tsurf are defined here adot=surf acc 0.1m/yr, Tsurf=-20oC, vary them later
end subroutine ice_shelf_temp
module subroutine ice_shelf_advect_temp_x(CS, G, time_step, hmask, h0, h_after_uflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses times temperature [C Z ~> degC m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


end subroutine ice_shelf_advect_temp_x
module subroutine ice_shelf_advect_temp_y(CS, G, time_step, hmask, h_after_uflux, h_after_vflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses times temperature after
                                              !! the meridional mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


end subroutine ice_shelf_advect_temp_y
  end interface

end module MOM_ice_shelf_dynamics
