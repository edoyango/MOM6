! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Energetically consistent planetary boundary layer parameterization
module MOM_energetic_PBL

use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_coms,           only : EFP_type, real_to_EFP, EFP_to_real, operator(+), assignment(=), EFP_sum_across_PEs
use MOM_debugging,      only : hchksum
use MOM_diag_mediator,  only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator,  only : post_data_3d_by_column, post_data_3d_final
use MOM_diag_mediator,  only : time_type, diag_ctrl
use MOM_domains,        only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,  only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser,    only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,   only : forcing
use MOM_grid,           only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_intrinsic_functions, only : cuberoot
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS, Get_Langmuir_Number
use MOM_stochastics,    only : stochastic_CS

implicit none ; private

#include <MOM_memory.h>

public energetic_PBL, energetic_PBL_init, energetic_PBL_end
public energetic_PBL_get_MLD

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters for the MOM_energetic_PBL module
type, public :: energetic_PBL_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.

  !/ Constants
  real    :: VonKar          !< The von Karman coefficient as used in the ePBL module [nondim]
  real    :: omega           !< The Earth's rotation rate [T-1 ~> s-1].
  real    :: omega_frac      !< When setting the decay scale for turbulence, use this fraction of
                             !! the absolute rotation rate blended with the local value of f, as
                             !! sqrt((1-omega_frac)*f^2 + omega_frac*4*omega^2) [nondim].

  !/ Convection related terms
  real    :: nstar           !< The fraction of the TKE input to the mixed layer available to drive
                             !! entrainment [nondim]. This quantity is the vertically integrated
                             !! buoyancy production minus the vertically integrated dissipation of
                             !! TKE produced by buoyancy.

  !/ Mixing Length terms
  logical :: Use_MLD_iteration !< If true, use the proximity to the bottom of the actively turbulent
                             !! surface boundary layer to constrain the mixing lengths.
  logical :: MLD_iteration_guess !< False to default to guessing half the
                             !! ocean depth for the first iteration.
  logical :: MLD_bisection   !< If true, use bisection with the iterative determination of the
                             !! self-consistent mixed layer depth.  Otherwise use the false position
                             !! after a maximum and minimum bound have been evaluated and the
                             !! returned value from the previous guess or bisection before this.
  logical :: MLD_iter_bug    !< If true use buggy logic that gives the wrong bounds for the next
                             !! iteration when successive guesses increase by exactly EPBL_MLD_TOLERANCE.
  integer :: max_MLD_its     !< The maximum number of iterations that can be used to find a
                             !! self-consistent mixed layer depth with Use_MLD_iteration.
  real    :: MixLenExponent  !< Exponent in the mixing length shape-function [nondim].
                             !! 1 is law-of-the-wall at top and bottom,
                             !! 2 is more KPP like.
  real    :: MKE_to_TKE_effic !< The efficiency with which mean kinetic energy released by
                             !!  mechanically forced entrainment of the mixed layer is converted to
                             !!  TKE, times conversion factors between the natural units of mean
                             !!  kinetic energy and those used for TKE [Z2 L-2 ~> nondim].
  logical :: direct_calc     !< If true and there is no conversion from mean kinetic energy to ePBL
                             !! turbulent kinetic energy, use a direct calculation of the
                             !! diffusivity that is supported by a given energy input instead of the
                             !! more general but slower iterative solver.
  real    :: ustar_min       !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
                             !! If the value is small enough, this should not affect the solution.
  real    :: Ekman_scale_coef !< A nondimensional scaling factor controlling the inhibition of the
                             !! diffusive length scale by rotation [nondim].  Making this larger decreases
                             !! the diffusivity in the planetary boundary layer.
  real    :: transLay_scale  !< A scale for the mixing length in the transition layer
                             !! at the edge of the boundary layer as a fraction of the
                             !! boundary layer thickness [nondim].  The default is 0, but a
                             !! value of 0.1 might be better justified by observations.
  real    :: MLD_tol         !< A tolerance for determining the boundary layer thickness when
                             !! Use_MLD_iteration is true [Z ~> m].
  real    :: min_mix_len     !< The minimum mixing length scale that will be used by ePBL [Z ~> m].
                             !! The default (0) does not set a minimum.

  !/ Velocity scale terms
  integer :: wT_scheme       !< An enumerated value indicating the method for finding the turbulent
                             !! velocity scale.  There are currently two options:
                             !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
                             !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real    :: wstar_ustar_coef !< A ratio relating the efficiency with which convectively released
                             !! energy is converted to a turbulent velocity, relative to
                             !! mechanically forced turbulent kinetic energy [nondim].
                             !! Making this larger increases the diffusivity.
  real    :: vstar_surf_fac  !< If (wT_scheme == wT_from_RH18) this is the proportionality coefficient between
                             !! ustar and the surface mechanical contribution to vstar [nondim]
  real    :: vstar_scale_fac !< An overall nondimensional scaling factor for vstar [nondim].  Making
                             !! this larger increases the diffusivity.

  !mstar related options
  integer :: mstar_scheme    !< An encoded integer to determine which formula is used to set mstar
  integer :: BBL_mstar_scheme !< An encoded integer to determine which formula is used to set mstar
  real    :: mstar_cap       !< Since mstar is restoring undissipated energy to mixing,
                             !! there must be a cap on how large it can be [nondim].  This
                             !! is definitely a function of latitude (Ekman limit),
                             !! but will be taken as constant for now.

  !/ vertical decay related options
  real    :: TKE_decay       !< The ratio of the natural Ekman depth to the TKE decay scale [nondim].

  !/ mstar_scheme == 0
  real    :: fixed_mstar     !< mstar is the ratio of the friction velocity cubed to the TKE available to
                             !! drive entrainment [nondim]. This quantity is the vertically
                             !! integrated shear production minus the vertically integrated
                             !! dissipation of TKE produced by shear.  This value is used if the option
                             !! for using a fixed mstar is used.
  real    :: BBL_fixed_mstar !< Similar to fixed_mstar, but for the bottom boundary layer

  !/ mstar_scheme == 2
  real :: C_Ek = 0.17        !< mstar Coefficient in rotation limit for EPBL_MSTAR_SCHEME=OM4 [nondim]
  real :: mstar_coef = 0.3   !< mstar coefficient in rotation/stabilizing balance for EPBL_MSTAR_SCHEME=OM4 [nondim]

  !/ mstar_scheme == 3
  real    :: RH18_mstar_cN1  !< mstar_N coefficient 1 (outer-most coefficient for fit) [nondim].
                             !! Value of 0.275 in RH18.  Increasing this
                             !! coefficient increases mechanical mixing for all values of Hf/ust,
                             !! but is most effective at low values (weakly developed OSBLs).
  real    :: RH18_mstar_cN2  !< mstar_N coefficient 2 (coefficient outside of exponential decay) [nondim].
                             !! Value of 8.0 in RH18.  Increasing this coefficient increases mstar
                             !! for all values of HF/ust, with a consistent affect across
                             !! a wide range of Hf/ust.
  real    :: RH18_mstar_cN3  !< mstar_N coefficient 3 (exponential decay coefficient) [nondim]. Value of
                             !! -5.0 in RH18.  Increasing this increases how quickly the value
                             !! of mstar decreases as Hf/ust increases.
  real    :: RH18_mstar_cS1  !< mstar_S coefficient for RH18 in stabilizing limit [nondim].
                             !! Value of 0.2 in RH18.
  real    :: RH18_mstar_cS2  !< mstar_S exponent for RH18 in stabilizing limit [nondim].
                             !! Value of 0.4 in RH18.

  !/ Coefficient for shear/convective turbulence interaction
  real :: mstar_convect_coef !< Factor to reduce mstar when statically unstable [nondim].

  !/ Langmuir turbulence related parameters
  logical :: Use_LT = .false. !< Flag for using LT in Energy calculation
  integer :: LT_enhance_form !< Integer for Enhancement functional form (various options)
  real    :: LT_enhance_coef !< Coefficient in fit for Langmuir Enhancement [nondim]
  real    :: LT_enhance_exp  !< Exponent in fit for Langmuir Enhancement [nondim]
  real :: LaC_MLD_Ek         !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Ekman depth [nondim].
  real :: LaC_MLD_Ob_stab    !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_Ek_Ob_stab     !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_MLD_Ob_un      !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: LaC_Ek_Ob_un       !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: Max_Enhance_M = 5. !< The maximum allowed LT enhancement to the mixing [nondim].

  !/ Machine learned equation discovery model paramters
  logical :: eqdisc       !< Uses machine learned shape function
  logical :: eqdisc_v0    !< Uses machine learned velocity scale
  logical :: eqdisc_v0h   !< Uses machine learned velocity scale that uses boundary layer depth as input
  real :: v0_lower_cap    !< Lower cap to prevent v0 from attaining anomlously low values [Z T-1 ~> m s-1]
  real :: v0_upper_cap    !< Upper cap to prevent v0 from attaining anomlously high values [Z T-1 ~> m s-1]
  real :: f_lower !< Lower cap of |f| i.e. absolute of Coriolis parameter [T-1 ~> s-1]
                  !! Used only in get_eqdisc_v0 subroutine. Default is 0.1deg Lat
  real :: bflux_lower_cap !< Lower cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: bflux_upper_cap !< Upper cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: sigma_max_lower_cap    !< Lower cap to prevent sigma_max from attaining low values [nondim]
  real :: sigma_max_upper_cap    !< Upper cap to prevent sigma_max from attaining high values [nondim]
  real :: Eh_upper_cap !< Upper cap to prevent Eh = hf/(u__*) from attaining high values [nondim]
  real :: Lh_cap       !< Cap to prevent Lh = h/Monin_Obukhov_depth from attaining beyond extreme values [nondim]
  real, allocatable, dimension(:) :: shape_function !< shape function used in machine learned diffusivity [nondim]
  !/ Coefficients used for Machine learned diffusivity
  real :: ML_c(18) !< Array of non-dimensional constants used in machine learned (ML) diffusivity [nondim]
  real :: shape_function_epsilon !< An small value of shape_function below the boundary layer depth [nondim]

  !/ Bottom boundary layer mixing related options
  real :: ePBL_BBL_effic     !< The efficiency of bottom boundary layer mixing via ePBL driven by
                             !! the bottom drag dissipation of mean kinetic energy, times
                             !! conversion factors between the natural units of mean kinetic energy
                             !! and those used for TKE [Z2 L-2 ~> nondim].
  real :: ePBL_tidal_effic   !< The efficiency of bottom boundary layer mixing via ePBL driven by
                             !! the bottom drag dissipation of tides, times conversion factors
                             !! between the natural units of mean kinetic energy and those used for
                             !! TKE [Z2 L-2 ~> nondim].
  logical :: Use_BBLD_iteration !< If true, use the proximity to the top of the actively turbulent
                             !! bottom boundary layer to constrain the mixing lengths.
  real    :: TKE_decay_BBL   !< The ratio of the natural Ekman depth to the TKE decay scale for
                             !! bottom boundary layer mixing [nondim]
  real    :: min_BBL_mix_len !< The minimum mixing length scale that will be used by ePBL in the bottom
                             !! boundary layer mixing [Z ~> m].  The default (0) does not set a minimum.
  real    :: MixLenExponent_BBL !< Exponent in the bottom boundary layer mixing length shape-function [nondim].
                             !! 1 is law-of-the-wall at top and bottom,
                             !! 2 is more KPP like.
  real    :: BBLD_tol        !< The tolerance for the iteratively determined bottom boundary layer depth [Z ~> m].
                             !! This is only used with USE_MLD_ITERATION.
  integer :: max_BBLD_its    !< The maximum number of iterations that can be used to find a self-consistent
                             !! bottom boundary layer depth.
  integer :: wT_scheme_BBL   !< An enumerated value indicating the method for finding the bottom boundary
                             !! layer turbulent velocity scale.  There are currently two options:
                             !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
                             !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real :: vstar_scale_fac_BBL !< An overall nondimensional scaling factor for wT in the bottom boundary layer [nondim].
                             !! Making this larger increases the bottom boundary layer diffusivity.", &
  real :: vstar_surf_fac_BBL !< If (wT_scheme_BBL == wT_from_RH18) this is the proportionality coefficient between
                             !! ustar and the bottom boundayer layer mechanical contribution to vstar [nondim]
  real :: Ekman_scale_coef_BBL !< A nondimensional scaling factor controlling the inhibition of the
                             !! diffusive length scale by rotation in the bottom boundary layer [nondim].
                             !! Making this larger decreases the bottom boundary layer diffusivity.
  logical :: decay_adjusted_BBL_TKE !< If true, include an adjustment factor in the bottom boundary layer
                             !! energetics that accounts for an exponential decay of TKE from a
                             !! near-bottom source and an assumed piecewise linear linear profile
                             !! of the buoyancy flux response to a change in a diffusivity.
  logical :: BBL_effic_bug   !< If true, overestimate the efficiency of the non-tidal ePBL bottom boundary
                             !! layer diffusivity by a factor of 1/sqrt(CDRAG), which is often a factor of
                             !! about 18.3.
  logical :: ePBL_BBL_use_mstar !< If true, use an mstar*ustar^3 paramaterization to get the TKE available
                             !! to drive mixing in the bottom boundary layer version of ePBL.  Otherwise,
                             !! use the meanflow energy loss to bottom drag scaled by a constant efficiency.

  !/ Options for documenting differences from parameter choices
  integer :: options_diff    !< If positive, this is a coded integer indicating a pair of
                             !! settings whose differences are diagnosed in a passive diagnostic mode
                             !! via extra calls to ePBL_column.  If this is 0 or negative no extra
                             !! calls occur.

  !/ Others
  type(time_type), pointer :: Time=>NULL() !< A pointer to the ocean model's clock.

  logical :: TKE_diagnostics = .false. !< If true, diagnostics of the TKE budget are being calculated.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in the ePBL
                             !! calculations.  Values below 20190101 recover the answers from the
                             !! end of 2018, while higher values use updated and more robust forms
                             !! of the same expressions.  Values below 20240101 use A**(1./3.) to
                             !! estimate the cube root of A in several expressions, while higher
                             !! values use the integer root function cuberoot(A) and therefore
                             !! can work with scaled variables.
  logical :: orig_PE_calc    !< If true, the ePBL code uses the original form of the
                             !! potential energy change code.  Otherwise, it uses a newer version
                             !! that can work with successive increments to the diffusivity in
                             !! upward or downward passes.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.

  real, allocatable, dimension(:,:) :: &
    ML_depth                 !< The mixed layer depth determined by active mixing in ePBL, which may
                             !! be used for the first guess in the next time step [H ~> m or kg m-2]
  real, allocatable, dimension(:,:) :: &
    BBL_depth                !< The bottom boundary layer depth determined by active mixing in ePBL [H ~> m or kg m-2]

  type(EFP_type), dimension(2) :: sum_its !< The total number of iterations and columns worked on
  type(EFP_type), dimension(2) :: sum_its_BBL !< The total number of iterations and columns worked on

  !>@{ Diagnostic IDs
  integer :: id_Kd_ePBL_col_by_col = -1
  integer :: id_ML_depth = -1, id_hML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_ustar_ePBL = -1, id_bflx_ePBL = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_Kd_BBL = -1, id_BBL_Mix_Length = -1, id_BBL_Vel_Scale = -1
  integer :: id_TKE_BBL = -1, id_TKE_BBL_mixing = -1, id_TKE_BBL_decay = -1
  integer :: id_ustar_BBL = -1, id_bflx_BBL = -1, id_BBL_decay_scale = -1, id_BBL_depth = -1
  integer :: id_mstar_sfc = -1, id_mstar_BBL = -1, id_LA_mod = -1, id_LA = -1, id_mstar_LT = -1
  ! The next options are used when passively diagnosing sensitivities from parameter choices
  integer :: id_opt_diff_Kd_ePBL = -1, id_opt_maxdiff_Kd_ePBL = -1, id_opt_diff_hML_depth = -1
  !>@}
end type energetic_PBL_CS

!>@{ Enumeration values for mstar_scheme
integer, parameter :: Use_Fixed_mstar = 0  !< The value of mstar_scheme to use a constant mstar
integer, parameter :: mstar_from_Ekman = 2 !< The value of mstar_scheme to base mstar on the ratio
                                           !! of the Ekman layer depth to the Obukhov depth
integer, parameter :: mstar_from_RH18 = 3  !< The value of mstar_scheme to base mstar of of RH18
integer, parameter :: No_Langmuir = 0      !< The value of LT_enhance_form not use Langmuir turbulence.
integer, parameter :: Langmuir_rescale = 2 !< The value of LT_enhance_form to use a multiplicative
                                           !! rescaling of mstar to account for Langmuir turbulence.
integer, parameter :: Langmuir_add = 3     !< The value of LT_enhance_form to add a contribution to
                                           !! mstar from Langmuir turbulence to other contributions.
integer, parameter :: wT_from_cRoot_TKE = 0 !< Use a constant times the cube root of remaining TKE
                                           !! to calculate the turbulent velocity.
integer, parameter :: wT_from_RH18 = 1     !< Use a scheme based on a combination of w* and v* as
                                           !! documented in Reichl & Hallberg (2018) to calculate
                                           !! the turbulent velocity.
character*(20), parameter :: CONSTANT_STRING = "CONSTANT"
character*(20), parameter :: OM4_STRING = "OM4"
character*(20), parameter :: RH18_STRING = "REICHL_H18"
character*(20), parameter :: ROOT_TKE_STRING = "CUBE_ROOT_TKE"
character*(20), parameter :: NONE_STRING = "NONE"
character*(20), parameter :: RESCALED_STRING = "RESCALE"
character*(20), parameter :: ADDITIVE_STRING = "ADDITIVE"
!>@}

logical :: report_avg_its = .false.  !< Report the average number of ePBL iterations for debugging.

!> A type for conveniently passing around ePBL diagnostics for a column.
type, public :: ePBL_column_diags ; private
  !>@{ Local column copies of energy change diagnostics, all in [R Z3 T-3 ~> W m-2].
  real :: dTKE_conv, dTKE_forcing, dTKE_wind, dTKE_mixing ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  real :: dTKE_MKE, dTKE_mech_decay, dTKE_conv_decay      ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  real :: dTKE_BBL, dTKE_BBL_decay, dTKE_BBL_mixing       ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  !>@}
  real :: LA        !< The value of the Langmuir number [nondim]
  real :: LAmod     !< The modified Langmuir number by convection [nondim]
  real :: mstar     !< The value of mstar used in ePBL [nondim]
  real :: mstar_BBL !< The value of mstar used in ePBL BBL [nondim]
  real :: mstar_LT  !< The portion of mstar due to Langmuir turbulence [nondim]
  integer :: OBL_its !< The number of iterations used to find a self-consistent surface boundary layer depth
  integer :: BBL_its !< The number of iterations used to find a self-consistent bottom boundary layer depth
end type ePBL_column_diags


  interface
module subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, visc, dt, Kd_int, G, GV, US, CS, &
                         stoch_CS, dSV_dT, dSV_dS, TKE_forced, buoy_flux, BBL_buoy_flux, Waves )
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_3d   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u_3d   !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: v_3d   !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: TKE_forced !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! [R Z3 T-2 ~> J m-2].
  type(thermo_var_ptrs),   intent(inout) :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  type(forcing),           intent(inout) :: fluxes !< A structure containing pointers to any
                                                   !! possible forcing fields. Unused fields have
                                                   !! NULL ptrs.
  type(vertvisc_type),     intent(in)    :: visc   !< Structure with vertical viscosities,
                                                   !! BBL properties and related fields
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: Kd_int !< The diagnosed diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(energetic_PBL_CS),  intent(inout) :: CS     !< Energetic PBL control structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: buoy_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: BBL_buoy_flux !< The bottom buoyancy flux [Z2 T-3 ~> m2 s-3].
  type(wave_parameters_CS), pointer      :: Waves  !< Waves control structure for Langmuir turbulence
  type(stochastic_CS),     pointer       :: stoch_CS  !< The control structure returned by a previous

!    This subroutine determines the diffusivities from the integrated energetics
!  mixed layer model.  It assumes that heating, cooling and freshwater fluxes
!  have already been applied.  All calculations are done implicitly, and there
!  is no stability limit on the time step.
!
!    For each interior interface, first discard the TKE to account for mixing
! of shortwave radiation through the next denser cell.  Next drive mixing based
! on the local? values of ustar + wstar, subject to available energy.  This
! step sets the value of Kd(K).  Any remaining energy is then subject to decay
! before being handed off to the next interface.  mech_TKE and conv_PErel are treated
! separately for the purposes of decay, but are used proportionately to drive
! mixing.
!
!   The key parameters for the mixed layer are found in the control structure.
!   To use the classic constant mstar mixed layers choose EPBL_MSTAR_SCHEME=CONSTANT.
! The key parameters then include mstar, nstar, TKE_decay, and conv_decay.
! For the Oberhuber (1993) mixed layer,the values of these are:
!      mstar = 1.25,  nstar = 1, TKE_decay = 2.5, conv_decay = 0.5
! TKE_decay is 1/kappa in eq. 28 of Oberhuber (1993), while conv_decay is 1/mu.
! For a traditional Kraus-Turner mixed layer, the values are:
!      mstar = 1.25, nstar = 0.4, TKE_decay = 0.0, conv_decay = 0.0

  ! Local variables
end subroutine energetic_PBL
module subroutine ePBL_column(h, dz, u, v, T0, S0, dSV_dT, dSV_dS, SpV_dt, TKE_forcing, B_flux, absf, &
                       u_star, u_star_mean, mech_TKE_in, dt, MLD_io, Kd, mixvel, mixlen, GV, US, CS, eCD, &
                       Waves, G, i, j, TKE_gen_stoch, TKE_diss_stoch)
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m].
  real, dimension(SZK_(GV)), intent(in)  :: u      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: v      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [C ~> degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [S ~> ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1), intent(in) :: SpV_dt !< Specific volume interpolated to interfaces
                                                   !! divided by dt or 1.0 / (dt * Rho0), times conversion
                                                   !! factors for answer dates before 20240101 in
                                                   !! [m3 Z-3 R-1 T2 s-3 ~> m3 kg-1 s-1] or without
                                                   !! the conversion factors for answer dates of
                                                   !! 20240101 and later in [R-1 T-1 ~> m3 kg-1 s-1],
                                                   !! used to convert local TKE into a turbulence
                                                   !! velocity cubed.
  real, dimension(SZK_(GV)), intent(in)  :: TKE_forcing !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! [R Z3 T-2 ~> J m-2].
  real,                    intent(in)    :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                    intent(in)    :: absf   !< The absolute value of the Coriolis parameter [T-1 ~> s-1].
  real,                    intent(in)    :: u_star !< The surface friction velocity [Z T-1 ~> m s-1].
  real,                    intent(in)    :: u_star_mean !< The surface friction velocity without any
                                                   !! contribution from unresolved gustiness  [Z T-1 ~> m s-1].
  real,                    intent(in)    :: mech_TKE_in !< The mechanically generated turbulent
                                                   !! kinetic energy available for mixing over a time
                                                   !! step before the application of the efficiency
                                                   !! in mstar. [R Z3 T-2 ~> J m-2].
  real,                    intent(inout) :: MLD_io !< A first guess at the mixed layer depth on input, and
                                                   !! the calculated mixed layer depth on output [Z ~> m]
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd     !< The diagnosed diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel !< The mixing velocity scale used in Kd
                                                   !! [Z T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen !< The mixing length scale used in Kd [Z ~> m].
  type(energetic_PBL_CS),  intent(in)    :: CS     !< Energetic PBL control structure
  type(ePBL_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.
  type(wave_parameters_CS), pointer      :: Waves  !< Waves control structure for Langmuir turbulence
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  integer,                 intent(in)    :: i      !< The i-index to work on (used for Waves)
  integer,                 intent(in)    :: j      !< The j-index to work on (used for Waves)
  real,          optional, intent(in)    :: TKE_gen_stoch  !< random factor used to perturb TKE generation [nondim]
  real,          optional, intent(in)    :: TKE_diss_stoch !< random factor used to perturb TKE dissipation [nondim]

!    This subroutine determines the diffusivities in a single column from the integrated energetics
!  planetary boundary layer (ePBL) model.  It assumes that heating, cooling and freshwater fluxes
!  have already been applied.  All calculations are done implicitly, and there
!  is no stability limit on the time step.
!
!    For each interior interface, first discard the TKE to account for mixing
! of shortwave radiation through the next denser cell.  Next drive mixing based
! on the local? values of ustar + wstar, subject to available energy.  This
! step sets the value of Kd(K).  Any remaining energy is then subject to decay
! before being handed off to the next interface.  mech_TKE and conv_PErel are treated
! separately for the purposes of decay, but are used proportionately to drive
! mixing.

  ! Local variables
end subroutine ePBL_column
module subroutine ePBL_BBL_column(h, dz, u, v, T0, S0, dSV_dT, dSV_dS, SpV_dt, absf, &
                           dt, Kd, BBL_TKE_in, u_star_BBL, u_star_BBL_z_t, b_flux_BBL, Kd_BBL, BBLD_io, mixvel_BBL, &
                           mixlen_BBL, GV, US, CS, eCD)
  type(verticalGrid_type),   intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m].
  real, dimension(SZK_(GV)), intent(in)  :: u      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: v      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [C ~> degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [S ~> ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1), intent(in) :: SpV_dt !< Specific volume interpolated to interfaces
                                                   !! divided by dt (if non-Boussinesq) or
                                                   !! 1.0 / (dt * Rho0), in [R-1 T-1 ~> m3 kg-1 s-1],
                                                   !! used to convert local TKE into a turbulence
                                                   !! velocity cubed.
  real,                    intent(in)    :: absf   !< The absolute value of the Coriolis parameter [T-1 ~> s-1].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(in)    :: Kd     !< The diffusivities at interfaces due to previously
                                                   !! applied mixing processes [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real,                    intent(in)    :: BBL_TKE_in !< The mechanically generated turbulent
                                                   !! kinetic energy available for bottom boundary
                                                   !! layer mixing within a time step [R Z3 T-2 ~> J m-2].
  real,                    intent(in)    :: u_star_BBL !< The bottom boundary layer friction velocity
                                                       !! in thickness flux units [H T-1 ~> m s-1 or kg m-2 s-1]
  real,                    intent(in)    :: u_star_BBL_z_t !< The bottom boundary layer friction velocity
                                                       !! converted to length flux units [Z T-1 ~> m s-1]
  real,                    intent(in)    :: b_flux_BBL !< The bottom boundary layer buoyancy flux
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd_BBL !< The bottom boundary layer contribution to diffusivities
                                                   !! at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real,                    intent(inout) :: BBLD_io !< A first guess at the bottom boundary layer depth on input, and
                                                   !! the calculated bottom boundary layer depth on output [Z ~> m]
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel_BBL !< The profile of boundary layer turbulent mixing
                                                   !! velocities [Z T-1 ~> m s-1]
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen_BBL !< The profile of bottom boundary layer turbulent
                                                   !! mixing lengths [Z ~> m]
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(energetic_PBL_CS),  intent(in)    :: CS     !< Energetic PBL control structure
  type(ePBL_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.

!    This subroutine determines the contributions from diffusivities in a single column from a
!  bottom-boundary layer adaptation of the integrated energetics planetary boundary layer (ePBL)
!  model.  It accounts for the possibility that the surface boundary diffusivities have already
!  been determined.  All calculations are done implicitly, and there is no stability limit on the
!  time step.  Only mechanical mixing in the bottom boundary layer is considered.  (Geothermal heat
!  fluxes are addressed elsewhere in the MOM6 code, and convection throughout the water column is
!  handled by the surface version of ePBL.)  There is no conversion of released mean kinetic energy
!  into bottom boundary layer turbulent kinetic energy (at least for now), apart from the explicit
!  energy that is supplied as an argument to this routine.

  ! Local variables
end subroutine ePBL_BBL_column
module subroutine kappa_eqdisc(shape_func, CS, GV, dz, absf, B_flux, u_star, MLD_guess)

  type(verticalGrid_type), intent(in) :: GV     !< The ocean's vertical grid structure.
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, dimension(SZK_(GV)+1), intent(inout) :: shape_func  !< shape function, [nondim]
  real, intent(in) :: absf      !< The absolute value of f [T-1 ~> s-1]
  real, intent(in) :: u_star    !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: B_Flux    !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m]
  real, intent(in) :: MLD_guess !< Mixing Layer depth guessed/found for iteration [Z ~> m].

  ! local variables for this subroutine

  ! variables used for optimizing computations:

end subroutine kappa_eqdisc
module subroutine get_eqdisc_v0(CS, absf, B_flux, u_star, v0_dummy)
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, intent(in) :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, intent(in) :: u_star !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: absf  !< The absolute value of f [T-1 ~> s-1].
  real, intent(inout) :: v0_dummy   !< velocity scale v0, local variable [Z T-1 ~> m s-1]

  ! local variables for this subroutine

end subroutine get_eqdisc_v0
module subroutine get_eqdisc_v0h(CS, B_flux, u_star, MLD_guess, v0_dummy)
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, intent(in) :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, intent(in) :: u_star !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: MLD_guess !< boundary layer depth guessed/found for iteration [Z ~> m]

  real, intent(inout) :: v0_dummy   !< velocity scale v0, local variable [Z T-1 ~> m s-1]

  ! local variables for this subroutine

end subroutine get_eqdisc_v0h
module function exp_decay_TKE_adjust(hb, ha, Idecay) result(TKE_to_PE_scale)
  real, intent(in) :: hb   !< The thickness over which the buoyancy flux varies on the
                           !! near-boundary side of an interface (e.g., a well-mixed bottom
                           !! boundary layer thickness) [H ~> m or kg m-2]
  real, intent(in) :: ha   !< The thickness of the layer on the opposite side of an interface from
                           !! the boundary [H ~> m or kg m-2]
  real, intent(in) :: Idecay !< The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1]
  real             :: TKE_to_PE_scale !< The effective fractional change in energy available to
                           !! drive mixing at this interface once the exponential decay of TKE
                           !! is accounted for [nondim].  TKE_to_PE_scale is always positive.


end function exp_decay_TKE_adjust
module subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres_Z, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, PE_ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface below [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].

  real, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                !! dKddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with dKddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realized by applying a huge value of dKddt_h at the
                                          !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with dKddt_h in the
                                            !! limit where dKddt_h = 0 [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: PE_ColHt_cor !< The correction to PE_chg that is made due to a net
                                            !! change in the column height [R Z3 T-2 ~> J m-2].

  ! Local variables
                   ! for the potential energy changes [R Z2 T-2 ~> J m-3].
                     ! for the column height changes [H Z ~> m2 or kg m-1].

  !   The expression for the change in potential energy used here is derived
  ! from the expression for the final estimates of the changes in temperature
  ! and salinities, and then extensively manipulated to get it into its most
  ! succinct form. The derivation is not necessarily obvious, but it demonstrably
  ! works by comparison with separate calculations of the energy changes after
  ! the tridiagonal solver for the final changes in temperature and salinity are
  ! applied.

end subroutine find_PE_chg
module subroutine find_Kd_from_PE_chg(Kd_prev, dKd_max, dt_h, max_PE_chg, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, pres_Z, &
                       dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       Kd_add, PE_chg, dPE_max, frac_dKd_max_PE)
  real, intent(in)  :: Kd_prev  !< The previously used diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(in)  :: dKd_max  !< The maximum change in the diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(in)  :: dt_h     !< The time step and divided by the average of the
                                !! thicknesses around the interface [T Z-1 ~> s m-1].
  real, intent(in)  :: max_PE_chg !< The maximum change in the column potential energy due to
                                !! additional mixing at an interface [R Z3 T-2 ~> J m-2].

  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface below [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].
  real, intent(out) :: Kd_add   !< The additional diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(out) :: PE_chg   !< The realized change in the column potential energy due to
                                !! additional mixing at an interface [R Z3 T-2 ~> J m-2].
  real, optional, &
        intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                !! be realized by applying a huge value of dKddt_h at the
                                !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, &
        intent(out) :: frac_dKd_max_PE !< The fraction of the energy required to support dKd_max
                                !! that is supplied by max_PE_chg [nondim]

  ! Local variables
                    ! and divided by the average of the thicknesses around the
                    ! interface [H ~> m or kg m-2].
                    ! the time step and divided by the average of the thicknesses around
                    ! the interface [H ~> m or kg m-2].
                    ! for the potential energy changes [R Z2 T-2 ~> J m-3].
                     ! for the column height changes [H Z ~> m2 or kg m-1].

  ! The expression for the change in potential energy used here is derived from the expression
  ! for the final estimates of the changes in temperature and salinities, which is then
  ! extensively manipulated to get it into its most succinct form.  It is the same as the
  ! expression that appears in find_PE_chg.

end subroutine find_Kd_from_PE_chg
module subroutine find_PE_chg_orig(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                       dT_to_dPEa, dS_to_dPEa, pres_Z, dT_to_dColHt_k, &
                       dS_to_dColHt_k, dT_to_dColHta, dS_to_dColHta, PE_chg, &
                       dPEc_dKd, dPE_max, dPEc_dKd_0)
  real, intent(in)  :: Kddt_h   !< The diffusivity at an interface times the time step and
                                !! divided by the average of the thicknesses around the
                                !! interface [H ~> m or kg m-2].
  real, intent(in)  :: h_k      !< The thickness of the layer below the interface [H ~> m or kg m-2].
  real, intent(in)  :: b_den_1  !< The first term in the denominator of the pivot
                                !! for the tridiagonal solver, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: dTe_term !< A diffusivity-independent term related to the temperature change
                                !! in the layer below the interface [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: dSe_term !< A diffusivity-independent term related to the salinity change
                                !! in the layer below the interface [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_km1_t2 !< A diffusivity-independent term related to the
                                 !! temperature change in the layer above the interface [C ~> degC].
  real, intent(in)  :: dS_km1_t2 !< A diffusivity-independent term related to the
                                 !! salinity change in the layer above the interface [S ~> ppt].
  real, intent(in)  :: pres_Z    !< The rescaled hydrostatic interface pressure, which relates
                                 !! the changes in column thickness to the energy that is radiated
                                 !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! in the salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].

  real, intent(out) :: PE_chg    !< The change in column potential energy from applying
                                 !! Kddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realized by applying a huge value of Kddt_h at the
                                          !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                          !! limit where Kddt_h = 0 [R Z3 T-2 H-1 ~> J m-3 or J kg-1].

!   This subroutine determines the total potential energy change due to mixing
! at an interface, including all of the implicit effects of the prescribed
! mixing at interfaces above.  Everything here is derived by careful manipulation
! of the robust tridiagonal solvers used for tracers by MOM6.  The results are
! positive for mixing in a stably stratified environment.
!   The comments describing these arguments are for a downward mixing pass, but
! this routine can also be used for an upward pass with the sense of direction
! reversed.

  ! Local variables
                        ! per unit change in Kddt_h [C H-1 ~> degC m-1 or degC m2 kg-1]
                        ! per unit change in Kddt_h [S H-1 ~> ppt m-1 or ppt m2 kg-1]

end subroutine find_PE_chg_orig
module subroutine find_mstar(CS, US, Buoyancy_Flux, UStar, &
                      BLD, Abs_Coriolis, Is_BBL, mstar, &
                      Langmuir_Number, mstar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), intent(in) :: CS    !< Energetic PBL control structure
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: UStar !< ustar including gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: Abs_Coriolis !< absolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  logical,               intent(in)  :: Is_BBL !< Logcal flag to indicate if bottom boundary layer mode
  real,                  intent(out) :: mstar !< Output mstar (Mixing/ustar**3) [nondim]
  real,        optional, intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,        optional, intent(out) :: mstar_LT !< mstar increase due to Langmuir turbulence [nondim]
  real,        optional, intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/ Variables used in computing mstar

  !/  Integer options for how to find mstar

  !/

end subroutine find_mstar
module subroutine mstar_Langmuir(CS, US, Abs_Coriolis, Buoyancy_Flux, UStar, BLD, Langmuir_Number, &
                          mstar, mstar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), intent(in) :: CS    !< Energetic PBL control structure
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: Abs_Coriolis !< Absolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: UStar !< Surface friction velocity with? gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  real,                  intent(inout) :: mstar !< Input/output mstar (Mixing/ustar**3) [nondim]
  real,                  intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,                  intent(out) :: mstar_LT !< mstar increase due to Langmuir turbulence [nondim]
  real,                  intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/
                             ! conditions or 0 under unstable conditions [nondim].
                             ! conditions or 0 under unstable conditions [nondim].
                             ! conditions or 0 under stable conditions [nondim].
                             ! conditions or 0 under stable conditions [nondim].

  ! Set default values for no Langmuir effects.
end subroutine mstar_Langmuir
module subroutine energetic_PBL_get_MLD(CS, MLD, G, US, m_to_MLD_units)
  type(energetic_PBL_CS),           intent(in)  :: CS  !< Energetic PBL control structure
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: MLD !< Depth of ePBL active mixing layer [Z ~> m]
                                                       !! or other units
  real,                   optional, intent(in)  :: m_to_MLD_units !< A conversion factor from meters
                                                       !! to the desired units for MLD, sometimes [Z m-1 ~> 1]
  ! Local variables

end subroutine energetic_PBL_get_MLD
module subroutine energetic_PBL_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(energetic_PBL_CS),  intent(inout) :: CS   !< Energetic PBL control structure

  ! Local variables
  ! This include declares and sets the variable "version".
                          ! recreate the bugs, or if false bugs are only used if actively selected.
                     ! bottom boundary layer mixing is not enabled.
end subroutine energetic_PBL_init
module subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), intent(inout) :: CS !< Energetic_PBL control structure


end subroutine energetic_PBL_end
  end interface

end module MOM_energetic_PBL
