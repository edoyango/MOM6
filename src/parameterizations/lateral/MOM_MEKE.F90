! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the Mesoscale Eddy Kinetic Energy framework
!! with topographic beta effect included in computing beta in Rhines scale

module MOM_MEKE

use iso_fortran_env,       only : real32

use MOM_coms,              only : PE_here
use MOM_database_comms,    only : dbclient_type, dbcomms_CS_type
use MOM_debugging,         only : hchksum, uvchksum
use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,     only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : pass_vector, pass_var
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg, is_root_pe
use MOM_file_parser,       only : read_param, get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_hor_index,         only : hor_index_type
use MOM_interface_heights, only : find_eta
use MOM_interpolate,       only : init_external_field, time_interp_external
use MOM_interpolate,       only : time_interp_external_init
use MOM_interpolate,       only : external_field
use MOM_io,                only : vardesc, var_desc, slasher
use MOM_isopycnal_slopes,  only : calc_isoneutral_slopes
use MOM_restart,           only : MOM_restart_CS, register_restart_field, query_initialized
use MOM_string_functions,  only : lowercase
use MOM_time_manager,      only : time_type_to_real
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : vertvisc_type, thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_MEKE_types,        only : MEKE_type


implicit none ; private

#include <MOM_memory.h>

public step_forward_MEKE, MEKE_init, MEKE_alloc_register_restart, MEKE_end

! Constants for this module
integer, parameter :: NUM_FEATURES = 4 !< How many features used to predict EKE
integer, parameter :: MKE_IDX = 1     !< Index of mean kinetic energy in the feature array
integer, parameter :: SLOPE_Z_IDX = 2 !< Index of vertically averaged isopycnal slope in the feature array
integer, parameter :: RV_IDX = 3      !< Index of surface relative vorticity in the feature array
integer, parameter :: RD_DX_Z_IDX = 4 !< Index of the radius of deformation over the grid size in the feature array

integer, parameter :: EKE_PROG = 1     !< Use prognostic equation to calculate EKE
integer, parameter :: EKE_FILE = 2     !< Read in EKE from a file
integer, parameter :: EKE_DBCLIENT = 3 !< Infer EKE using a neural network

!> Control structure that contains MEKE parameters and diagnostics handles
type, public :: MEKE_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  ! Parameters
  real :: MEKE_FrCoeff  !< Efficiency of conversion of ME into MEKE [nondim]
  real :: MEKE_bhFrCoeff!< Efficiency of conversion of ME into MEKE by the biharmonic dissipation [nondim]
  real :: MEKE_GMcoeff  !< Efficiency of conversion of PE into MEKE [nondim]
  real :: MEKE_GMECoeff !< Efficiency of conversion of MEKE into ME by GME [nondim]
  real :: MEKE_damping  !< Local depth-independent MEKE dissipation rate [T-1 ~> s-1].
  real :: MEKE_Cd_scale !< The ratio of the bottom eddy velocity to the column mean
                        !! eddy velocity, i.e. sqrt(2*MEKE), [nondim]. This should be less than 1
                        !! to account for the surface intensification of MEKE.
  real :: MEKE_Cb       !< Coefficient in the \f$\gamma_{bot}\f$ expression [nondim]
  real :: MEKE_min_gamma!< Minimum value of gamma_b^2 allowed [nondim]
  real :: MEKE_Ct       !< Coefficient in the \f$\gamma_{bt}\f$ expression [nondim]
  logical :: visc_drag  !< If true use the vertvisc_type to calculate bottom drag.
  logical :: MEKE_GEOMETRIC !< If true, uses the GM coefficient formulation from the GEOMETRIC
                        !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha !< The nondimensional coefficient governing the efficiency of the
                        !! GEOMETRIC thickness diffusion [nondim].
  logical :: MEKE_equilibrium_alt !< If true, use an alternative calculation for the
                        !! equilibrium value of MEKE.
  logical :: MEKE_equilibrium_restoring !< If true, restore MEKE back to its equilibrium value,
                        !!  which is calculated at each time step.
  logical :: GM_src_alt !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                        !! than the streamfunction for the MEKE GM source term.
  real    :: MEKE_min_depth_tot  !< The minimum total thickness over which to distribute MEKE energy
                        !! sources from GM energy conversion [H ~> m or kg m-2].  When the total
                        !! thickness is less than this, the sources are scaled away.
  logical :: Rd_as_max_scale !< If true the length scale can not exceed the
                        !! first baroclinic deformation radius.
  logical :: use_old_lscale !< Use the old formula for mixing length scale.
  logical :: use_min_lscale !< Use simple minimum for mixing length scale.
  logical :: MEKE_positive  !< If true, it guarantees that MEKE will always be >= 0.
  real :: lscale_maxval !< The ceiling on the MEKE mixing length scale when use_min_lscale is true [L ~> m].
  real :: cdrag         !< The bottom drag coefficient for MEKE, times rescaling factors [H L-1 ~> nondim or kg m-3]
  real :: MEKE_BGsrc    !< Background energy source for MEKE [L2 T-3 ~> W kg-1] (= m2 s-3).
  real :: MEKE_dtScale  !< Scale factor to accelerate time-stepping [nondim]
  real :: MEKE_KhCoeff  !< Scaling factor to convert MEKE into Kh [nondim]
  real :: MEKE_Uscale   !< MEKE velocity scale for bottom drag [L T-1 ~> m s-1]
  real :: MEKE_KH       !< Background lateral diffusion of MEKE [L2 T-1 ~> m2 s-1]
  real :: MEKE_K4       !< Background bi-harmonic diffusivity (of MEKE) [L4 T-1 ~> m4 s-1]
  real :: KhMEKE_Fac    !< A factor relating MEKE%Kh to the diffusivity used for
                        !! MEKE itself [nondim].
  real :: viscosity_coeff_Ku !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral harmonic momentum mixing
                        !! by unresolved eddies represented by MEKE [nondim].
  real :: viscosity_coeff_Au !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral biharmonic momentum mixing
                        !! by unresolved eddies represented by MEKE [nondim].
  real :: Lfixed        !< Fixed mixing length scale [L ~> m].
  real :: aDeform       !< Weighting towards deformation scale of mixing length [nondim]
  real :: aRhines       !< Weighting towards Rhines scale of mixing length [nondim]
  real :: aFrict        !< Weighting towards frictional arrest scale of mixing length [nondim]
  real :: aEady         !< Weighting towards Eady scale of mixing length [nondim]
  real :: aGrid         !< Weighting towards grid scale of mixing length [nondim]
  real :: MEKE_advection_factor !< A scaling in front of the advection of MEKE [nondim]
  real :: MEKE_topographic_beta !< Weight for how much topographic beta is considered
                                !! when computing beta in Rhines scale [nondim]
  real :: MEKE_restoring_rate !< Inverse of the timescale used to nudge MEKE toward its
                        !! equilibrium value [T-1 ~> s-1].
  logical :: MEKE_advection_bug !< If true, recover a bug in the calculation of the barotropic
                        !! transport for the advection of MEKE, wherein only the transports in the
                        !! deepest layer are used.
  logical :: fixed_total_depth  !< If true, use the nominal bathymetric depth as the estimate of
                        !! the time-varying ocean depth.  Otherwise base the depth on the total
                        !! ocean mass per unit area.
  real :: rho_fixed_total_depth !< A density used to translate the nominal bathymetric depth into an
                        !! estimate of the total ocean mass per unit area when MEKE_FIXED_TOTAL_DEPTH
                        !! is true [R ~> kg m-3]
  logical :: kh_flux_enabled !< If true, lateral diffusive MEKE flux is enabled.
  logical :: initialize !< If True, invokes a steady state solver to calculate MEKE.
  logical :: debug      !< If true, write out checksums of data for debugging
  integer :: eke_src !< Enum specifying whether EKE is stepped forward prognostically (default),
                     !! read in from a file, or inferred via a neural network
  logical :: sqg_use_MEKE !< If True, use MEKE%Le for the SQG vertical structure.
  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_MEKE = -1, id_Ue = -1, id_Kh = -1, id_src = -1
  integer :: id_src_adv = -1, id_src_mom_K4 = -1, id_src_btm_drag = -1
  integer :: id_src_GM = -1, id_src_mom_lp = -1, id_src_mom_bh = -1
  integer :: id_Ub = -1, id_Ut = -1
  integer :: id_GM_src = -1, id_mom_src = -1, id_mom_src_bh = -1, id_GME_snk = -1, id_decay = -1
  integer :: id_KhMEKE_u = -1, id_KhMEKE_v = -1, id_Ku = -1, id_Au = -1
  integer :: id_Le = -1, id_gamma_b = -1, id_gamma_t = -1
  integer :: id_Lrhines = -1, id_Leady = -1
  integer :: id_MEKE_equilibrium = -1
  !>@}
  type(external_field) :: eke_handle   !< Handle for reading in EKE from a file
  ! Infrastructure
  integer :: id_clock_pass !< Clock for group pass calls
  type(group_pass_type) :: pass_MEKE !< Group halo pass handle for MEKE%MEKE and maybe MEKE%Kh_diff
  type(group_pass_type) :: pass_Kh   !< Group halo pass handle for MEKE%Kh, MEKE%Ku, and/or MEKE%Au

  ! MEKE via Machine Learning
  type(dbclient_type), pointer :: client => NULL() !< Pointer to the database client

  logical :: online_analysis !< If true, post the EKE used in MOM6 at every timestep
  character(len=5) :: model_key  = 'mleke'  !< Key where the ML-model is stored
  character(len=7) :: key_suffix !< Suffix appended to every key sent to Redis
  real :: eke_max !< The maximum value of EKE considered physically reasonable [L2 T-2 ~> m2 s-2]

  ! Clock ids
  integer :: id_client_init   !< Clock id to time initialization of the client
  integer :: id_put_tensor    !< Clock id to time put_tensor routine
  integer :: id_run_model     !< Clock id to time running of the ML model
  integer :: id_unpack_tensor !< Clock id to time retrieval of EKE prediction

  ! Diagnostic ids
  integer :: id_mke     = -1 !< Diagnostic id for surface mean kinetic energy
  integer :: id_slope_z = -1 !< Diagnostic id for vertically averaged horizontal slope magnitude
  integer :: id_slope_x = -1 !< Diagnostic id for isopycnal slope in the x-direction
  integer :: id_slope_y = -1 !< Diagnostic id for isopycnal slope in the y-direction
  integer :: id_rv      = -1 !< Diagnostic id for surface relative vorticity

end type MEKE_CS


  interface
module subroutine step_forward_MEKE(MEKE, h, SN_u, SN_v, visc, dt, G, GV, US, CS, hu, hv, u, v, tv, Time)
  type(MEKE_type),                          intent(inout) :: MEKE !< MEKE data.
  type(ocean_grid_type),                    intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),                  intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),                    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: h    !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  type(vertvisc_type),                      intent(in)    :: visc !< The vertical viscosity type.
  real,                                     intent(in)    :: dt   !< Model(baroclinic) time-step [T ~> s].
  type(MEKE_CS),                            intent(inout) :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: hu   !< Accumulated zonal mass flux [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: hv   !< Accumulated meridional mass flux [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                    intent(in)    :: tv   !< Type containing thermodynamic variables
  type(time_type),                          intent(in)    :: Time !< The time used for interpolating EKE

  ! Local variables
end subroutine step_forward_MEKE
module subroutine MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass, depth_tot)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  type(MEKE_type),                   intent(inout) :: MEKE !< MEKE fields
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: drag_rate_visc !< Mean flow velocity contribution
                                                           !! to the MEKE drag rate [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: I_mass  !< Inverse of column mass [R-1 Z-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].

  ! Local variables

end subroutine MEKE_equilibrium
module subroutine MEKE_equilibrium_restoring(CS, G, GV, US, SN_u, SN_v, depth_tot, &
                                      equilibrium_value)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type.
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: equilibrium_value
      !< Equilbrium value of MEKE to be calculated at each time step [L2 T-2 ~> m2 s-2]

  ! Local variables

end subroutine MEKE_equilibrium_restoring
module subroutine MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, EKE, depth_tot, &
                             bottomFac2, barotrFac2, LmixScale)
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  type(MEKE_type),                   intent(in)    :: MEKE !< MEKE field
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: EKE  !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: bottomFac2 !< gamma_b^2 [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: barotrFac2 !< gamma_t^2 [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: LmixScale !< Eddy mixing length [L ~> m].
  ! Local variables

end subroutine MEKE_lengthScales
module subroutine MEKE_lengthScales_0d(CS, US, area, beta, depth_tot, Rd_dx, SN, EKE, &
                                bottomFac2, barotrFac2, LmixScale, Lrhines, Leady)
  type(MEKE_CS), intent(in)    :: CS         !< MEKE control structure.
  type(unit_scale_type), intent(in) :: US    !< A dimensional unit scaling type
  real,          intent(in)    :: area       !< Grid cell area [L2 ~> m2]
  real,          intent(in)    :: beta       !< Planetary beta = \f$ \nabla f\f$  [T-1 L-1 ~> s-1 m-1]
  real,          intent(in)    :: depth_tot  !< The total thickness of the water column [H ~> m or kg m-2]
  real,          intent(in)    :: Rd_dx      !< Resolution Ld/dx [nondim].
  real,          intent(in)    :: SN         !< Eady growth rate [T-1 ~> s-1].
  real,          intent(in)    :: EKE        !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real,          intent(out)   :: bottomFac2 !< gamma_b^2 [nondim]
  real,          intent(out)   :: barotrFac2 !< gamma_t^2 [nondim]
  real,          intent(out)   :: LmixScale  !< Eddy mixing length [L ~> m].
  real,          intent(out)   :: Lrhines    !< Rhines length scale [L ~> m].
  real,          intent(out)   :: Leady      !< Eady length scale [L ~> m].
  ! Local variables

  ! Length scale for MEKE derived diffusivity
end subroutine MEKE_lengthScales_0d
logical module function MEKE_init(Time, G, GV, US, param_file, diag, dbcomms_CS, CS, MEKE, restart_CS, meke_in_dynamics)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Ocean vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(dbcomms_CS_type),   intent(in)    :: dbcomms_CS !< Database communications control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(MEKE_CS),           intent(inout) :: CS         !< MEKE control structure.
  type(MEKE_type),         intent(inout) :: MEKE       !< MEKE fields
  type(MOM_restart_CS),    intent(in)    :: restart_CS !< MOM restart control structure
  logical,                 intent(  out) :: meke_in_dynamics !< If true, MEKE is stepped forward in dynamics
                                                             !! otherwise in tracer dynamics

  ! Local variables
  ! This include declares and sets the variable "version".

end function MEKE_init
module subroutine ML_MEKE_init(diag, G, US, Time, param_file, dbcomms_CS, CS)
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ocean_grid_type),         intent(inout) :: G           !< The ocean's grid structure.
  type(unit_scale_type),         intent(in)    :: US          !< A dimensional unit scaling type
  type(time_type),               intent(in)    :: Time        !< The current model time.
  type(param_file_type),         intent(in)    :: param_file  !< Parameter file parser structure.
  type(dbcomms_CS_type),         intent(in)    :: dbcomms_CS  !< Control structure for database communication
  type(MEKE_CS),                 intent(inout) :: CS          !< Control structure for this module


  ! Store pointers in control structure
end subroutine ML_MEKE_init
module subroutine ML_MEKE_calculate_features(G, GV, US, CS, Rd_dx_h, u, v, tv, h, dt, features_array)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(MEKE_CS),                             intent(in)    :: CS !< Control structure for MEKE
  real, dimension(SZI_(G),SZJ_(G)),          intent(in   ) :: Rd_dx_h !< Rossby radius of deformation over
                                                                 !! the grid length scale [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Type containing thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt !< Model(baroclinic) time-step [T ~> s].
  real(kind=real32), dimension(SIZE(h),num_features), intent(  out) :: features_array
                                                                 !< The array of features needed for machine
                                                                 !! learning inference, with different units
                                                                 !! for the various subarrays [various]




end subroutine ML_MEKE_calculate_features
module subroutine predict_MEKE(G, US, CS, npts, Time, features_array, MEKE)
  type(ocean_grid_type),                                 intent(inout) :: G  !< Ocean grid
  type(unit_scale_type),                                 intent(in)    :: US   !< A dimensional unit scaling type
  type(MEKE_CS),                                         intent(in   ) :: CS !< Control structure for MEKE
  integer,                                               intent(in   ) :: npts !< Number of T-grid cells on the local
                                                                               !! domain
  type(time_type),                                       intent(in   ) :: Time !< The current model time
  real(kind=real32), dimension(npts,num_features),       intent(in   ) :: features_array
                                                                          !< The array of features needed for machine
                                                                          !! learning inference, with different units
                                                                          !! for the various subarrays [various]
  real, dimension(SZI_(G),SZJ_(G)),                      intent(  out) :: MEKE !< Eddy kinetic energy [L2 T-2 ~> m2 s-2]

  ! Local variables
                                                       ! energy in mks units [m2 s-2]
                                                       ! in mks units [m2 s-2]

end subroutine predict_MEKE
real module function vertical_average_interface(h, w, h_min)

  real, dimension(:), intent(in) :: h  !< Layer Thicknesses [H ~> m or kg m-2]
  real, dimension(:), intent(in) :: w  !< Quantity to average [arbitrary]
  real, intent(in) :: h_min !< The vanishingly small layer thickness [H ~> m or kg m-2]


end function vertical_average_interface
module subroutine MEKE_alloc_register_restart(HI, US, param_file, MEKE, restart_CS)
! Arguments
  type(hor_index_type),  intent(in)    :: HI         !< Horizontal index structure
  type(unit_scale_type), intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in)    :: param_file !< Parameter file parser structure.
  type(MEKE_type),       intent(inout) :: MEKE       !< MEKE fields
  type(MOM_restart_CS),  intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables

! Determine whether this module will be used
end subroutine MEKE_alloc_register_restart
module subroutine MEKE_end(MEKE)
  type(MEKE_type), intent(inout) :: MEKE !< A structure with MEKE-related fields.

  ! NOTE: MEKE will always be allocated by MEKE_init, even if MEKE is disabled.
  !  So these must all be conditional, even though MEKE%MEKE and MEKE%Rd_dx_h
  !  are always allocated (when MEKE is enabled)

end subroutine MEKE_end
  end interface

end module MOM_MEKE
