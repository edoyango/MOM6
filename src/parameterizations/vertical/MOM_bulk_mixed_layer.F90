! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Build mixed layer parameterization
module MOM_bulk_mixed_layer

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type, diag_ctrl, diag_update_remap_grids
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_EOS,           only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,           only : average_specific_vol, calculate_density_derivs
use MOM_EOS,           only : calculate_spec_vol, calculate_specific_vol_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : extractFluxes1d, forcing, find_ustar
use MOM_grid,          only : ocean_grid_type
use MOM_opacity,       only : absorbRemainingSW, optics_type, extract_optics_slice
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public bulkmixedlayer, bulkmixedlayer_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters for the MOM_bulk_mixed_layer module
type, public :: bulkmixedlayer_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  integer :: nkml            !< The number of layers in the mixed layer.
  integer :: nkbl            !< The number of buffer layers.
  integer :: nsw             !< The number of bands of penetrating shortwave radiation.
  real    :: mstar           !< The ratio of the friction velocity cubed to the
                             !! TKE input to the mixed layer [nondim].
  real    :: nstar           !< The fraction of the TKE input to the mixed layer
                             !! available to drive entrainment [nondim].
  real    :: nstar2          !< The fraction of potential energy released by
                             !! convective adjustment that drives entrainment [nondim].
  logical :: absorb_all_SW   !< If true, all shortwave radiation is absorbed by the
                             !! ocean, instead of passing through to the bottom mud.
  real    :: TKE_decay       !< The ratio of the natural Ekman depth to the TKE
                             !! decay scale [nondim].
  real    :: bulk_Ri_ML      !< The efficiency with which mean kinetic energy released by
                             !! mechanically forced entrainment of the mixed layer is
                             !! converted to TKE, times conversion factors between the
                             !! natural units of mean kinetic energy and TKE [Z2 L-2 ~> nondim]
  real    :: bulk_Ri_convective !< The efficiency with which convectively released mean kinetic
                             !! energy becomes TKE, times conversion factors between the natural
                             !! units of mean kinetic energy and TKE [Z2 L-2 ~> nondim]
  real    :: vonKar          !< The von Karman constant as used for mixed layer viscosity [nondim]
  real    :: Hmix_min        !< The minimum mixed layer thickness [H ~> m or kg m-2].
  real    :: mech_TKE_floor  !< A tiny floor on the amount of turbulent kinetic energy that is
                             !! used when the mixed layer does not yet contain HMIX_MIN fluid
                             !! [H Z2 T-2 ~> m3 s-2 or J m-2].  The default is so small that its actual
                             !! value is irrelevant, but it is detectably greater than 0.
  real    :: H_limit_fluxes  !< When the total ocean depth is less than this
                             !! value [H ~> m or kg m-2], scale away all surface forcing to
                             !! avoid boiling the ocean.
  real    :: ustar_min       !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
                             !! If the value is small enough, this should not affect the solution.
  real    :: omega           !<   The Earth's rotation rate [T-1 ~> s-1].
  real    :: dT_dS_wt        !<   When forced to extrapolate T & S to match the
                             !! layer densities, this factor [C S-1 ~> degC ppt-1] is
                             !! combined with the derivatives of density with T & S
                             !! to determines what direction is orthogonal to
                             !! density contours.  It should be a typical value of
                             !! (dR/dS) / (dR/dT) in oceanic profiles.
                             !! 6 degC ppt-1 might be reasonable.
  real    :: Hbuffer_min     !< The minimum buffer layer thickness when the mixed layer
                             !! is very large [H ~> m or kg m-2].
  real    :: Hbuffer_rel_min !< The minimum buffer layer thickness relative to the combined
                             !! mixed and buffer layer thicknesses when they are thin [nondim]
  real    :: BL_detrain_time !< A timescale that characterizes buffer layer detrainment
                             !! events [T ~> s].
  real    :: BL_extrap_lim   !< A limit on the density range over which
                             !! extrapolation can occur when detraining from the
                             !! buffer layers, relative to the density range
                             !! within the mixed and buffer layers, when the
                             !! detrainment is going into the lightest interior
                             !! layer  [nondim].
  real :: BL_split_rho_tol   !< The fractional tolerance for matching layer target densities
                             !! when splitting layers to deal with massive interior layers
                             !! that are lighter than one of the mixed or buffer layers [nondim].
  logical :: ML_resort       !<   If true, resort the layers by density, rather than
                             !! doing convective adjustment.
  integer :: ML_presort_nz_conv_adj !< If ML_resort is true, do convective
                             !! adjustment on this many layers (starting from the
                             !! top) before sorting the remaining layers.
  real    :: omega_frac      !<   When setting the decay scale for turbulence, use this fraction
                             !! of the absolute rotation rate blended with the local value of f,
                             !! as sqrt((1-of)*f^2 + of*4*omega^2) [nondim].
  logical :: correct_absorption !< If true, the depth at which penetrating
                             !! shortwave radiation is absorbed is corrected by
                             !! moving some of the heating upward in the water
                             !! column.  The default is false.
  logical :: nonBous_energetics  !< If true, use non-Boussinesq expressions for the energetic
                             !! calculations used in the bulk mixed layer calculations.
  logical :: Resolve_Ekman   !<   If true, the nkml layers in the mixed layer are
                             !! chosen to optimally represent the impact of the
                             !! Ekman transport on the mixed layer TKE budget.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  logical :: TKE_diagnostics = .false. !< If true, calculate extensive diagnostics of the TKE budget
  logical :: do_rivermix = .false. !< Provide additional TKE to mix river runoff
                             !! at the river mouths to rivermix_depth
  real    :: rivermix_depth = 0.0  !< The depth of mixing if do_rivermix is true [H ~> m or kg m-2].
  logical :: limit_det       !< If true, limit the extent of buffer layer
                             !! detrainment to be consistent with neighbors.
  real    :: lim_det_dH_sfc  !< The fractional limit in the change between grid
                             !! points of the surface region (mixed & buffer
                             !! layer) thickness [nondim].  0.5 by default.
  real    :: lim_det_dH_bathy !< The fraction of the total depth by which the
                             !! thickness of the surface region (mixed & buffer layers) is allowed
                             !! to change between grid points [nondim].  0.2 by default.
  logical :: use_river_heat_content !< If true, use the fluxes%runoff_Hflx field
                             !! to set the heat carried by runoff, instead of
                             !! using SST for temperature of liq_runoff
  logical :: use_calving_heat_content !< Use SST for temperature of froz_runoff
  logical :: convect_mom_bug !< If true, use code with a bug that causes a loss of momentum
                             !! conservation during mixedlayer convection.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  real    :: Allowed_T_chg   !< The amount by which temperature is allowed
                             !! to exceed previous values during detrainment [C ~> degC]
  real    :: Allowed_S_chg   !< The amount by which salinity is allowed
                             !! to exceed previous values during detrainment [S ~> ppt]

  ! These are terms in the mixed layer TKE budget, all in [H Z2 T-3 ~> m3 s-3 or W m-2] except as noted.
  real, allocatable, dimension(:,:) :: &
    ML_depth, &        !< The mixed layer depth [H ~> m or kg m-2].
    diag_TKE_wind, &   !< The wind source of TKE [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_RiBulk, & !< The resolved KE source of TKE [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_conv, &   !< The convective source of TKE [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_pen_SW, & !< The TKE sink required to mix penetrating shortwave heating [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_mech_decay, & !< The decay of mechanical TKE [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_conv_decay, & !< The decay of convective TKE [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_mixing, & !< The work done by TKE to deepen the mixed layer [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_TKE_conv_s2, & !< The convective source of TKE due to to mixing in sigma2 [H Z2 T-3 ~> m3 s-3 or W m-2].
    diag_PE_detrain, & !< The spurious source of potential energy due to mixed layer
                       !! detrainment [R Z3 T-3 ~> W m-2].
    diag_PE_detrain2   !< The spurious source of potential energy due to mixed layer only
                       !! detrainment [R Z3 T-3 ~> W m-2].
  type(group_pass_type) :: pass_h_sum_hmbl_prev !< For group halo pass

  !>@{ Diagnostic IDs
  integer :: id_ML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_TKE_RiBulk = -1, id_TKE_conv = -1, id_TKE_pen_SW = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1, id_TKE_conv_s2 = -1
  integer :: id_PE_detrain = -1, id_PE_detrain2 = -1, id_h_mismatch = -1
  integer :: id_Hsfc_used = -1, id_Hsfc_max = -1, id_Hsfc_min = -1
  !>@}
end type bulkmixedlayer_CS

!>@{ CPU clock IDs
integer :: id_clock_pass=0
!>@}


  interface
module subroutine bulkmixedlayer(h_3d, u_3d, v_3d, tv, fluxes, dt, ea, eb, G, GV, US, CS, &
                          optics, BLD, H_ml, aggregate_FW_forcing, dt_diag, last_call)
  type(ocean_grid_type),      intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_3d   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: u_3d   !< Zonal velocities interpolated to h points
                                                      !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: v_3d   !< Zonal velocities interpolated to h points
                                                      !! [L T-1 ~> m s-1].
  type(thermo_var_ptrs),      intent(inout) :: tv     !< A structure containing pointers to any
                                                      !! available thermodynamic fields. Absent
                                                      !! fields have NULL pointers.
  type(forcing),              intent(inout) :: fluxes !< A structure containing pointers to any
                                                      !! possible forcing fields.  Unused fields
                                                      !! have NULL pointers.
  real,                       intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea     !< The amount of fluid moved downward into a
                                                      !! layer; this should be increased due to
                                                      !! mixed layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb     !< The amount of fluid moved upward into a
                                                      !! layer; this should be increased due to
                                                      !! mixed layer entrainment [H ~> m or kg m-2].
  type(bulkmixedlayer_CS),    intent(inout) :: CS     !< Bulk mixed layer control structure
  type(optics_type),          pointer       :: optics !< The structure that can be queried for the
                                                      !! inverse of the vertical absorption decay
                                                      !! scale for penetrating shortwave radiation.
  real, dimension(SZI_(G),SZJ_(G)), &
                              intent(inout) :: BLD    !< Active mixed layer depth [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                              intent(inout) :: H_ml   !< Active mixed layer thickness [H ~> m or kg m-2].
  logical,                    intent(in)    :: aggregate_FW_forcing !< If true, the net incoming and
                                                     !! outgoing surface freshwater fluxes are
                                                     !! combined before being applied, instead of
                                                     !! being applied separately.
  real,             optional, intent(in)    :: dt_diag  !< The diagnostic time step,
                                                      !! which may be less than dt if there are
                                                      !! two calls to mixedlayer [T ~> s].
  logical,          optional, intent(in)    :: last_call !< if true, this is the last call
                                                      !! to mixedlayer in the current time step, so
                                                      !! diagnostics will be written. The default is
                                                      !! .true.

  ! Local variables
end subroutine bulkmixedlayer
module subroutine convective_adjustment(h, u, v, R0, SpV0, Rcv, T, S, eps, d_eb, &
                                 dKE_CA, cTKE, j, G, GV, US, CS, nz_conv)
  type(ocean_grid_type),              intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV  !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2].
                                                           !! The units of h are referred to as H below.
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: u   !< Zonal velocities interpolated to h
                                                           !! points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: v   !< Zonal velocities interpolated to h
                                                           !! points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: R0  !< Potential density referenced to
                                                           !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: SpV0 !< Specific volume referenced to
                                                           !! surface pressure [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: Rcv !< The coordinate defining potential
                                                           !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: T   !< Layer temperatures [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: S   !< Layer salinities [S ~> ppt].
  real, dimension(SZI_(G),SZK_(GV)),  intent(in)    :: eps !< The negligibly small amount of water
                                                           !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_eb !< The downward increase across a layer
                                                           !! in the entrainment from below [H ~> m or kg m-2].
                                                           !! Positive values go with mass gain by
                                                           !! a layer.
  real, dimension(SZI_(G),SZK_(GV)),  intent(out)   :: dKE_CA !< The vertically integrated change in
                                                           !! kinetic energy due to convective
                                                           !! adjustment [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZK_(GV)),  intent(out)   :: cTKE !< The buoyant turbulent kinetic energy
                                                           !! source due to convective adjustment
                                                           !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  integer,                            intent(in)    :: j   !< The j-index to work on.
  type(unit_scale_type),              intent(in)    :: US  !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),            intent(in)    :: CS  !< Bulk mixed layer control structure
  integer,                  optional, intent(in)    :: nz_conv !< If present, the number of layers
                                                           !! over which to do convective adjustment
                                                           !! (perhaps CS%nkml).

  ! Local variables
end subroutine convective_adjustment
module subroutine mixedlayer_convection(h, d_eb, htot, Ttot, Stot, uhtot, vhtot,      &
                                 R0_tot, SpV0_tot, Rcv_tot, u, v, T, S, R0, SpV0, Rcv, eps,    &
                                 dR0_dT, dSpV0_dT, dRcv_dT, dR0_dS, dSpV0_dS, dRcv_dS,             &
                                 netMassInOut, netMassOut, Net_heat, Net_salt, &
                                 nsw, Pen_SW_bnd, opacity_band, Conv_En,       &
                                 dKE_FC, j, ksort, G, GV, US, CS, tv, fluxes, dt,      &
                                 aggregate_FW_forcing)
  type(ocean_grid_type),    intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV    !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(inout) :: h     !< Layer thickness [H ~> m or kg m-2].
                                                   !! The units of h are referred to as H below.
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: d_eb  !< The downward increase across a layer in the
                                                   !! layer in the entrainment from below [H ~> m or kg m-2].
                                                   !! Positive values go with mass gain by a layer.
  real, dimension(SZI_(G)), intent(out)   :: htot  !< The accumulated mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: Ttot  !< The depth integrated mixed layer temperature
                                                   !! [C H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: Stot  !< The depth integrated mixed layer salinity
                                                   !! [S H ~> ppt m or ppt kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: uhtot !< The depth integrated mixed layer zonal
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(out)   :: vhtot !< The integrated mixed layer meridional
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(out)   :: R0_tot !< The integrated mixed layer potential density referenced
                                                   !! to 0 pressure [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G)), intent(out)   :: SpV0_tot !< The integrated mixed layer specific volume referenced
                                                   !! to 0 pressure [H R-1 ~> m4 kg-1 or m].
  real, dimension(SZI_(G)), intent(out)   :: Rcv_tot !< The integrated mixed layer coordinate
                                                   !! variable potential density [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: u     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: v     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: T     !< Layer temperatures [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: S     !< Layer salinities [S ~> ppt].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: R0    !< Potential density referenced to
                                                   !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: SpV0  !< Specific volume referenced to
                                                   !! surface pressure [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: Rcv   !< The coordinate defining potential
                                                   !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: eps   !< The negligibly small amount of water
                                                   !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dT  !< The partial derivative of R0 with respect to
                                                   !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dSpV0_dT  !< The partial derivative of SpV0 with respect to
                                                   !! temperature [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dT !< The partial derivative of Rcv with respect to
                                                   !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dS  !< The partial derivative of R0 with respect to
                                                   !! salinity [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)), intent(in)    :: dSpV0_dS  !< The partial derivative of SpV0 with respect to
                                                   !! salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dS !< The partial derivative of Rcv with respect to
                                                   !! salinity [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)), intent(in)    :: netMassInOut !< The net mass flux (if non-Boussinesq)
                                                   !! or volume flux (if Boussinesq) into the ocean
                                                   !! within a time step [H ~> m or kg m-2]. (I.e. P+R-E.)
  real, dimension(SZI_(G)), intent(in)    :: netMassOut !< The mass or volume flux out of the ocean
                                                   !! within a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: Net_heat !< The net heating at the surface over a time
                                                   !! step [C H ~> degC m or degC kg m-2].  Any penetrating
                                                   !! shortwave radiation is not included in Net_heat.
  real, dimension(SZI_(G)), intent(in)    :: Net_salt !< The net surface salt flux into the ocean
                                                   !! over a time step [S H ~> ppt m or ppt kg m-2].
  integer,                  intent(in)    :: nsw   !< The number of bands of penetrating
                                                   !! shortwave radiation.
  real, dimension(max(nsw,1),SZI_(G)), intent(inout) :: Pen_SW_bnd !< The penetrating shortwave
                                                   !! heating at the sea surface in each penetrating
                                                   !! band [C H ~> degC m or degC kg m-2].
  real, dimension(max(nsw,1),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< The opacity in each band of
                                                   !! penetrating shortwave radiation [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)), intent(out)   :: Conv_En !< The buoyant turbulent kinetic energy source
                                                   !! due to free convection [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)), intent(out)   :: dKE_FC !< The vertically integrated change in kinetic
                                                   !! energy due to free convection [H Z2 T-2 ~> m3 s-2 or J m-2].
  integer,                  intent(in)    :: j     !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: ksort !< The density-sorted k-indices.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),  intent(in)    :: CS    !< Bulk mixed layer control structure
  type(thermo_var_ptrs),    intent(inout) :: tv    !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent
                                                   !! fields have NULL pointers.
  type(forcing),            intent(inout) :: fluxes  !< A structure containing pointers to any
                                                   !! possible forcing fields.  Unused fields
                                                   !! have NULL pointers.
  real,                     intent(in)    :: dt    !< Time increment [T ~> s].
  logical,                  intent(in)    :: aggregate_FW_forcing !< If true, the net incoming and
                                                   !! outgoing surface freshwater fluxes are
                                                   !! combined before being applied, instead of
                                                   !! being applied separately.

!   This subroutine causes the mixed layer to entrain to the depth of free
! convection.  The depth of free convection is the shallowest depth at which the
! fluid is denser than the average of the fluid above.

  ! Local variables
                       ! that is not absorbed in a layer [nondim].
                       ! that is absorbed in a layer [C H ~> degC m or degC kg m-2].
                       ! entrainment [H ~> m or kg m-2].
                       ! h_ent between iterations [H ~> m or kg m-2].
                       ! the conversion from H to Z divided by the mean density,
                       ! [Z2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
                       ! shortwave radiation, integrated over a layer
                       ! [H R ~> kg m-2 or kg2 m-5].

end subroutine mixedlayer_convection
module subroutine find_starting_TKE(htot, h_CA, fluxes, U_star_2d, Conv_En, cTKE, dKE_FC, dKE_CA, &
                             TKE, TKE_river, Idecay_len_TKE, cMKE, tv, dt, Idt_diag, &
                             j, ksort, G, GV, US, CS)
  type(ocean_grid_type),      intent(in)    :: G       !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  real, dimension(SZI_(G)),   intent(in)    :: htot    !< The accumulated mixed layer thickness
                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G)),   intent(in)    :: h_CA    !< The mixed layer depth after convective
                                                       !! adjustment [H ~> m or kg m-2].
  type(forcing),              intent(in)    :: fluxes  !< A structure containing pointers to any
                                                       !! possible forcing fields.  Unused fields
                                                       !! have NULL pointers.
  real, dimension(SZI_(G),SZJ_(G)), intent(in) ::  U_star_2d !< The wind friction velocity, calculated
                                                       !! using the Boussinesq reference density or
                                                       !! the time-evolving surface density in
                                                       !! non-Boussinesq mode [Z T-1 ~> m s-1]
  real, dimension(SZI_(G)),   intent(inout) :: Conv_En !< The buoyant turbulent kinetic energy source
                                                       !! due to free convection [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),   intent(in)    :: dKE_FC  !< The vertically integrated change in
                                                       !! kinetic energy due to free convection
                                                       !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: cTKE    !< The buoyant turbulent kinetic energy
                                                       !! source due to convective adjustment
                                                       !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: dKE_CA  !< The vertically integrated change in
                                                       !! kinetic energy due to convective
                                                       !! adjustment [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),   intent(out)   :: TKE     !< The turbulent kinetic energy available for
                                                       !! mixing over a time step [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G)),   intent(out)   :: Idecay_len_TKE !< The inverse of the vertical decay
                                                       !! scale for TKE [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)),   intent(in)    :: TKE_river !< The source of turbulent kinetic energy
                                                       !! available for driving mixing at river mouths
                                                       !! [H Z2 T-3 ~> m3 s-3 or W m-2].
  real, dimension(2,SZI_(G)), intent(out)   :: cMKE    !< Coefficients of HpE and HpE^2 in
                                                       !! calculating the denominator of MKE_rate,
                                                       !! [H-1 ~> m-1 or m2 kg-1] and [H-2 ~> m-2 or m4 kg-2].
  type(thermo_var_ptrs),      intent(inout) :: tv      !< A structure containing pointers to any
                                                       !! available thermodynamic fields.
  real,                       intent(in)    :: dt      !< The time step [T ~> s].
  real,                       intent(in)    :: Idt_diag !< The inverse of the accumulated diagnostic
                                                       !! time interval [T-1 ~> s-1].
  integer,                    intent(in)    :: j       !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: ksort   !< The density-sorted k-indices.
  type(bulkmixedlayer_CS),    intent(inout) :: CS      !< Bulk mixed layer control structure

!   This subroutine determines the TKE available at the depth of free
! convection to drive mechanical entrainment.

  ! Local variables
                    ! free convection is converted to TKE, often ~0.2 [nondim].
                    ! convective adjustment is converted to TKE, often ~0.2 [nondim].
                    ! that release is positive [H Z2 T-2 ~> m3 s-2 or J m-2].
                    ! timestep (which may include 2 calls) [nondim].
                    ! based on the layer-averaged specific volume [Z H-1 ~> nondim or m3 kg-1]

end subroutine find_starting_TKE
module subroutine mechanical_entrainment(h, d_eb, htot, Ttot, Stot, uhtot, vhtot, &
                                  R0_tot, SpV0_tot, Rcv_tot, u, v, T, S, R0, SpV0, Rcv, eps, &
                                  dR0_dT, dSpV0_dT, dRcv_dT, cMKE, Idt_diag, nsw, &
                                  Pen_SW_bnd, opacity_band, TKE, &
                                  Idecay_len_TKE, j, ksort, G, GV, US, CS)
  type(ocean_grid_type),    intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV    !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(inout) :: h     !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: d_eb  !< The downward increase across a layer in the
                                                   !! layer in the entrainment from below [H ~> m or kg m-2].
                                                   !! Positive values go with mass gain by a layer.
  real, dimension(SZI_(G)), intent(inout) :: htot  !< The accumulated mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: Ttot  !< The depth integrated mixed layer temperature
                                                   !! [C H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: Stot  !< The depth integrated mixed layer salinity
                                                   !! [S H ~> ppt m or ppt kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: uhtot !< The depth integrated mixed layer zonal
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(inout) :: vhtot !< The integrated mixed layer meridional
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(inout) :: R0_tot !< The integrated mixed layer potential density
                                                   !! referenced to 0 pressure [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G)), intent(inout) :: SpV0_tot !< The integrated mixed layer specific volume referenced
                                                   !! to 0 pressure [H R-1 ~> m4 kg-1 or m].
  real, dimension(SZI_(G)), intent(inout) :: Rcv_tot !< The integrated mixed layer coordinate variable
                                                   !! potential density [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: u     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: v     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: T     !< Layer temperatures [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: S     !< Layer salinities [S ~> ppt].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: R0    !< Potential density referenced to
                                                   !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: SpV0  !< Specific volume referenced to
                                                   !! surface pressure [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G),SZK0_(GV)), &
                            intent(in)    :: Rcv   !< The coordinate defining potential
                                                   !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: eps   !< The negligibly small amount of water
                                                   !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dT  !< The partial derivative of R0 with respect to
                                                   !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dSpV0_dT  !< The partial derivative of SpV0 with respect to
                                                   !! temperature [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dT !< The partial derivative of Rcv with respect to
                                                   !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(2,SZI_(G)), intent(in)  :: cMKE  !< Coefficients of HpE and HpE^2 used in calculating the
                                                   !! denominator of MKE_rate; the two elements have differing
                                                   !! units of [H-1 ~> m-1 or m2 kg-1] and [H-2 ~> m-2 or m4 kg-2].
  real,                     intent(in)    :: Idt_diag !< The inverse of the accumulated diagnostic
                                                   !! time interval [T-1 ~> s-1].
  integer,                  intent(in)    :: nsw   !< The number of bands of penetrating
                                                   !! shortwave radiation.
  real, dimension(max(nsw,1),SZI_(G)), intent(inout) :: Pen_SW_bnd !< The penetrating shortwave
                                                   !! heating at the sea surface in each penetrating
                                                   !! band [C H ~> degC m or degC kg m-2].
  real, dimension(max(nsw,1),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< The opacity in each band of
                                                   !! penetrating shortwave radiation [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)), intent(inout) :: TKE   !< The turbulent kinetic energy
                                                   !! available for mixing over a time
                                                   !! step [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)), intent(inout) :: Idecay_len_TKE !< The vertical TKE decay rate [H-1 ~> m-1 or m2 kg-1].
  integer,                  intent(in)    :: j     !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: ksort !< The density-sorted k-indices.
  type(bulkmixedlayer_CS),  intent(inout) :: CS    !< Bulk mixed layer control structure

! This subroutine calculates mechanically driven entrainment.

  ! Local variables
                    ! absorbed in a layer [nondim].
                        ! that is absorbed in a layer [C H ~> degC m or degC kg m-2].
                       ! h_ent between iterations [H ~> m or kg m-2].
                    ! within the mixed layer that will be eliminated
                    ! within a timestep [nondim], 0 to 1.
                      ! conversion from H to m divided by the mean density,
                      ! in [Z2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
                        ! [H Z2 T-2 ~> m3 s-2 or J m-2].
                    ! across the mixed layer [Z2 T-2 ~> m2 s-2].
                          ! TKE, divided by layer thickness in m [Z2 T-2 ~> m2 s-2].
                    ! kinetic energy [H2 Z2 T-2 ~> m4 s-2 or kg2 m-2 s-2]
                    ! release of mean kinetic energy [H Z2 T-2 ~> m3 s-2 or J m-2]
                    ! dTKE_dh [Z2 T-2 ~> m2 s-2].
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
                    ! fractional decay of TKE across a layer [nondim].
                    ! of TKE and SW radiation across a layer [nondim]
                    ! thicknesses [H-1 ~> m-1 or m2 kg-1].

end subroutine mechanical_entrainment
module subroutine sort_ML(h, R0, SpV0, eps, G, GV, CS, ksort)
  type(ocean_grid_type),                intent(in)  :: G     !< The ocean's grid structure.
  type(verticalGrid_type),              intent(in)  :: GV    !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)),   intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(in)  :: R0    !< The potential density used to sort
                                                             !! the layers [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(in)  :: SpV0  !< Specific volume referenced to
                                                             !! surface pressure [R-1 ~> m3 kg-1]
  real, dimension(SZI_(G),SZK_(GV)),    intent(in)  :: eps   !< The (small) thickness that must
                                                             !! remain in each layer [H ~> m or kg m-2].
  type(bulkmixedlayer_CS),              intent(in)  :: CS    !< Bulk mixed layer control structure
  integer, dimension(SZI_(G),SZK_(GV)), intent(out) :: ksort !< The k-index to use in the sort.

  ! Local variables

end subroutine sort_ML
module subroutine resort_ML(h, T, S, R0, SpV0, Rcv, RcvTgt, eps, d_ea, d_eb, ksort, G, GV, CS, &
                     dR0_dT, dR0_dS, dSpV0_dT, dSpV0_dS, dRcv_dT, dRcv_dS)
  type(ocean_grid_type),                intent(in)    :: G       !< The ocean's grid structure.
  type(verticalGrid_type),              intent(in)    :: GV      !< The ocean's vertical grid
                                                                 !! structure.
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: h       !< Layer thickness [H ~> m or kg m-2].
                                                                 !! Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: T       !< Layer temperatures [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: S       !< Layer salinities [S ~> ppt].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: R0      !< Potential density referenced to
                                                                 !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: SpV0    !< Specific volume referenced to
                                                                 !! surface pressure [R-1 ~> m3 kg-1]
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: Rcv     !< The coordinate defining
                                                                 !! potential density [R ~> kg m-3].
  real, dimension(SZK_(GV)),            intent(in)    :: RcvTgt  !< The target value of Rcv for each
                                                                 !! layer [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: eps     !< The (small) thickness that must
                                                                 !! remain in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: d_ea    !< The upward increase across a
                                                                 !! layer in the entrainment from
                                                                 !! above [H ~> m or kg m-2].
                                                                 !!  Positive d_ea goes with layer
                                                                 !! thickness increases.
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: d_eb    !< The downward increase across a
                                                                 !! layer in the entrainment from
                                                                 !! below [H ~> m or kg m-2]. Positive values go
                                                                 !! with mass gain by a layer.
  integer, dimension(SZI_(G),SZK_(GV)), intent(in)    :: ksort   !< The density-sorted k-indices.
  type(bulkmixedlayer_CS),              intent(in)    :: CS      !< Bulk mixed layer control structure
  real, dimension(SZI_(G)),             intent(in)    :: dR0_dT  !< The partial derivative of
                                                                 !! potential density referenced
                                                                 !! to the surface with potential
                                                                 !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),             intent(in)    :: dR0_dS  !< The partial derivative of
                                                                 !! potential density referenced
                                                                 !! to the surface with salinity,
                                                                 !! [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),             intent(in)    :: dSpV0_dT !< The partial derivative of SpV0 with respect
                                                                 !! to temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G)),             intent(in)    :: dSpV0_dS !< The partial derivative of SpV0 with respect
                                                                 !! to salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
  real, dimension(SZI_(G)),             intent(in)    :: dRcv_dT !< The partial derivative of
                                                                 !! coordinate defining potential
                                                                 !! density with potential
                                                                 !! temperature [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),             intent(in)    :: dRcv_dS !< The partial derivative of
                                                                 !! coordinate defining potential
                                                                 !! density with salinity,
                                                                 !! [R S-1 ~> kg m-3 ppt-1].

!   If there are no massive light layers above the deepest of the mixed- and
! buffer layers, do nothing (except perhaps to reshuffle these layers).
!   If there are nkbl or fewer layers above the deepest mixed- or buffer-
! layers, move them (in sorted order) into the buffer layers, even if they
! were previously interior layers.
!   If there are interior layers that are intermediate in density (both in-situ
! and the coordinate density (sigma-2)) between the newly forming mixed layer
! and a residual buffer- or mixed layer, and the number of massive layers above
! the deepest massive buffer or mixed layer is greater than nkbl, then split
! those buffer layers into pieces that match the target density of the two
! nearest interior layers.
!   Otherwise, if there are more than nkbl+1 remaining massive layers

  ! Local variables
                        ! when extrapolating to match a target density [C2 S-2 ~> degC2 ppt-2]
                        ! extrapolating [C R-1 ~> degC m3 kg-1]
                        ! extrapolating [S R-1 ~> ppt m3 kg-1]
                        ! densities of two layers [R ~> kg m-3]
                        ! densities of two layers [R-1 ~> m3 kg-1]
                        ! pair of layers [R H2 ~> kg m-1 or kg3 m-7] or [R-1 H2 ~> m5 kg-1 or kg m-1]

end subroutine resort_ML
module subroutine mixedlayer_detrain_2(h, T, S, R0, Spv0, Rcv, RcvTgt, dt, dt_diag, d_ea, j, G, GV, US, CS, &
                                dR0_dT, dR0_dS, dSpV0_dT, dSpV0_dS, dRcv_dT, dRcv_dS, max_BL_det)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: h    !< Layer thickness [H ~> m or kg m-2].
                                                            !!  Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: T    !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: S    !< Salinity [S ~> ppt].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: R0   !< Potential density referenced to
                                                            !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: SpV0 !< Specific volume referenced to
                                                            !! surface pressure [R-1 ~> m3 kg-1]
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: Rcv  !< The coordinate defining potential
                                                            !! density [R ~> kg m-3].
  real, dimension(SZK_(GV)),          intent(in)    :: RcvTgt  !< The target value of Rcv for each
                                                            !! layer [R ~> kg m-3].
  real,                               intent(in)    :: dt   !< Time increment [T ~> s].
  real,                               intent(in)    :: dt_diag !< The diagnostic time step [T ~> s].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_ea !< The upward increase across a layer in
                                                            !! the entrainment from above
                                                            !! [H ~> m or kg m-2]. Positive d_ea
                                                            !! goes with layer thickness increases.
  integer,                            intent(in)    :: j    !< The meridional row to work on.
  type(unit_scale_type),              intent(in)    :: US   !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),            intent(inout) :: CS   !< Bulk mixed layer control structure
  real, dimension(SZI_(G)),           intent(in)    :: dR0_dT  !< The partial derivative of
                                                            !! potential density referenced to the
                                                            !! surface with potential temperature,
                                                            !! [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dR0_dS  !< The partial derivative of
                                                            !! potential density referenced to the
                                                            !! surface with salinity
                                                            !! [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: dSpV0_dT !< The partial derivative of specific
                                                            !! volume with respect to temeprature
                                                            !! [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G)),           intent(in)    :: dSpV0_dS  !< The partial derivative of specific
                                                            !! volume with respect to salinity
                                                            !! [R-1 S-1 ~> m3 kg-1 ppt-1]
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dT !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with potential temperature,
                                                            !! [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dS !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with salinity [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: max_BL_det !< If non-negative, the maximum
                                                            !! detrainment permitted from the buffer
                                                            !! layers [H ~> m or kg m-2].

! This subroutine moves any water left in the former mixed layers into the
! two buffer layers and may also move buffer layer water into the interior
! isopycnal layers.

  ! Local variables
                                  ! layers [H ~> m or kg m-2].
                                  ! buffer layer [H R ~> kg m-2 or kg2 m-5]
                                  ! buffer layer [H R-1 ~> m4 kg-1 or m]
                                  ! buffer layer [H R ~> kg m-2 or kg2 m-5]
                                  ! buffer layer [C H ~> degC m or degC kg m-2]
                                  ! buffer layer [S H ~> ppt m or ppt kg m-2]

                                  ! h(i,CS%nkml+1) and h(i,CS%nkml+2) [H ~> m or kg m-2].
                                  ! available to move into the lower buffer
                                  ! layer [H ~> m or kg m-2].
                                  ! layer that remains there [H ~> m or kg m-2].
                                  ! stays [H ~> m or kg m-2].

                                  ! between the water in kb2 and the water being detrained.
                                  ! buffer layers and create water that matches
                                  ! the target density of an interior layer.
                                  ! stays_merge is the thickness of the upper
                                  ! layer that remains [H ~> m or kg m-2].

!  real :: dT_2dz                 ! Half the vertical gradient of T [C H-1 ~> degC m-1 or degC m2 kg-1]
!  real :: dS_2dz                 ! Half the vertical gradient of S [S H-1 ~> ppt m-1 or ppt m2 kg-1]
                                  ! the slope within the upper buffer layer when
                                  ! water MUST be detrained to the lower layer [nondim].

                                  ! advection or mixing layers, divided by
                                  ! rho_0*g [H2 ~> m2 or kg2 m-4].
                                  ! mixing layers [R Z3 T-2 ~> J m-2].
                                  ! into the buffer layer or the merge the two
                                  ! buffer layers [R H2 Z T-2 ~> J m-2 or J kg2 m-8].
                                  ! into the buffer layer or the merge the two
                                  ! buffer layers [R Z3 T-2 ~> J m-2].

                                  ! drawn from the mixed layer [H ~> m or kg m-2].
                                  ! water that will go directly into the lower
                                  ! buffer layer [H ~> m or kg m-2].

                                  ! the lower buffer layer [H ~> m or kg m-2].
                                  ! the upper buffer layer [H ~> m or kg m-2].
                                  ! and to an interior layer that is just denser than the lower
                                  ! buffer layer [H ~> m or kg m-2].
                                  ! is just denser than the lower buffer layer [H ~> m or kg m-2].

                                  ! buffer layer and the water that moves into
                                  ! an interior layer or that stays in that
                                  ! layer [R ~> kg m-3].
                                  ! the lower buffer layer and the water that
                                  ! moves into an interior layer [R ~> kg m-3].
                                  ! advection [R H-1 ~> kg m-4 or m-1].
                                  ! buffer layer and the water that stays in that layer [R-1 ~> m3 kg-1]
                                  ! between the lower buffer layer and the water that
                                  ! moves into an interior layer [R-1 ~> m3 kg-1]
                                  ! permitted - here (detrainment_per_day/dt)*30
                                  ! days? [nondim]
                                  ! to prefer merging the buffer layers [nondim].
                                  ! salinity changes in defining spiciness, in
                                  ! [C S-1 ~> degC ppt-1] and [S C-1 ~> ppt degC-1].

                                  ! divided by the time step [Z2 H-2 T-1 ~> s-1 or m6 kg-2 s-1].
                                  ! respect to the coordinate potential density.
                    ! in roundoff and can be neglected [H ~> m or kg m-2].



end subroutine mixedlayer_detrain_2
module subroutine mixedlayer_detrain_1(h, T, S, R0, SpV0, Rcv, RcvTgt, dt, dt_diag, d_ea, d_eb, &
                                j, G, GV, US, CS, dRcv_dT, dRcv_dS, max_BL_det)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: h    !< Layer thickness [H ~> m or kg m-2].
                                                            !! Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: T    !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: S    !< Salinity [S ~> ppt].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: R0   !< Potential density referenced to
                                                            !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: SpV0 !< Specific volume referenced to
                                                            !! surface pressure [R-1 ~> m3 kg-1]
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: Rcv  !< The coordinate defining potential
                                                            !! density [R ~> kg m-3].
  real, dimension(SZK_(GV)),          intent(in)    :: RcvTgt !< The target value of Rcv for each
                                                            !! layer [R ~> kg m-3].
  real,                               intent(in)    :: dt   !< Time increment [T ~> s].
  real,                               intent(in)    :: dt_diag !< The accumulated time interval for
                                                            !! diagnostics [T ~> s].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_ea !< The upward increase across a layer in
                                                            !! the entrainment from above
                                                            !! [H ~> m or kg m-2]. Positive d_ea
                                                            !! goes with layer thickness increases.
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_eb !< The downward increase across a layer
                                                            !! in the entrainment from below [H ~> m or kg m-2].
                                                            !! Positive values go with mass gain by
                                                            !! a layer.
  integer,                            intent(in)    :: j    !< The meridional row to work on.
  type(unit_scale_type),              intent(in)    :: US   !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),            intent(inout) :: CS   !< Bulk mixed layer control structure
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dT !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with potential temperature
                                                            !! [R C-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dS    !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with salinity [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: max_BL_det !< If non-negative, the maximum
                                                            !! detrainment permitted from the buffer
                                                            !! layers [H ~> m or kg m-2].

  ! Local variables
                              ! entrained [H ~> m or kg m-2].
                              ! from the mixed layer [H ~> m or kg m-2].
                     ! when extraploating to match a target density [C2 S-2 ~> degC2 ppt-2]
                     ! extrapolating [C R-1 ~> degC m3 kg-1]
                     ! extrapolating [S R-1 ~> ppt m3 kg-1]
                              ! conversion from H to m divided by the mean density times the time
                              ! step [Z2 T-3 H-1 R-1 ~> m4 s-3 kg-1 or m7 s-3 kg-2].
                              ! conversion from H to Z divided by the diagnostic time step
                              ! [Z3 H-2 T-3 ~> m s-3 or m7 kg-2 s-3].
                              ! H to RZ divided by the diagnostic time step
                              ! [R Z2 H-1 T-3 ~> kg m-2 s-3 or m s-3].
                              ! H to RZ squared divided by the diagnostic time step
                              ! [R2 Z3 H-2 T-3 ~> kg2 m-5 s-3 or m s-3]

end subroutine mixedlayer_detrain_1
module subroutine bulkmixedlayer_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The model's clock with the current time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(bulkmixedlayer_CS), intent(inout) :: CS   !< Bulk mixed layer control structure

  ! This include declares and sets the variable "version".
end subroutine bulkmixedlayer_init
module function EF4(Ht, En, I_L, dR_de)
  real,           intent(in)    :: Ht  !< Total thickness [H ~> m or kg m-2].
  real,           intent(in)    :: En  !< Entrainment [H ~> m or kg m-2].
  real,           intent(in)    :: I_L !< The e-folding scale [H-1 ~> m-1 or m2 kg-1]
  real, optional, intent(inout) :: dR_de !< The partial derivative of the result R with E [H-2 ~> m-2 or m4 kg-2].
  real :: EF4 !< The integral [H-1 ~> m-1 or m2 kg-1].

  ! Local variables

end function EF4
  end interface

end module MOM_bulk_mixed_layer
