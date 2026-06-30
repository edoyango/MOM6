! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This routine drives the diabatic/dianeutral physics for MOM
module MOM_diabatic_driver

use MOM_bulk_mixed_layer,    only : bulkmixedlayer, bulkmixedlayer_init, bulkmixedlayer_CS
use MOM_debugging,           only : hchksum
use MOM_checksum_packages,   only : MOM_state_chksum, MOM_state_stats
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_CVMix_shear,         only : CVMix_shear_is_used
use MOM_CVMix_ddiff,         only : CVMix_ddiff_is_used
use MOM_diabatic_aux,        only : diabatic_aux_init, diabatic_aux_end, diabatic_aux_CS
use MOM_diabatic_aux,        only : make_frazil, adjust_salt, differential_diffuse_T_S, triDiagTS
use MOM_diabatic_aux,        only : triDiagTS_Eulerian, find_uv_at_h
use MOM_diabatic_aux,        only : applyBoundaryFluxesInOut, set_pen_shortwave
use MOM_diag_mediator,       only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,       only : post_product_sum_u, post_product_sum_v
use MOM_diag_mediator,       only : diag_ctrl, time_type, diag_update_remap_grids
use MOM_diag_mediator,       only : diag_ctrl, query_averaging_enabled, enable_averages, disable_averaging
use MOM_diag_mediator,       only : diag_grid_storage, diag_grid_storage_init, diag_grid_storage_end
use MOM_diag_mediator,       only : diag_copy_diag_to_storage, diag_copy_storage_to_diag
use MOM_diag_mediator,       only : diag_save_grids, diag_restore_grids
use MOM_diagnose_mld,        only : diagnoseMLDbyDensityDifference, diagnoseMLDbyEnergy
use MOM_diagnose_kdwork,     only : vbf_CS, KdWork_init, KdWork_end, KdWork_diagnostics
use MOM_diagnose_kdwork,     only : Allocate_VBF_CS, Deallocate_VBF_CS
use MOM_diapyc_energy_req,   only : diapyc_energy_req_init, diapyc_energy_req_end
use MOM_diapyc_energy_req,   only : diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_CS
use MOM_CVMix_conv,          only : CVMix_conv_init, CVMix_conv_cs
use MOM_CVMix_conv,          only : calculate_CVMix_conv
use MOM_domains,             only : pass_var, To_West, To_South, To_All, Omit_Corners
use MOM_domains,             only : create_group_pass, do_group_pass, group_pass_type
use MOM_energetic_PBL,       only : energetic_PBL, energetic_PBL_init
use MOM_energetic_PBL,       only : energetic_PBL_end, energetic_PBL_CS
use MOM_energetic_PBL,       only : energetic_PBL_get_MLD
use MOM_entrain_diffusive,   only : entrainment_diffusive, entrain_diffusive_init
use MOM_entrain_diffusive,   only : entrain_diffusive_CS
use MOM_EOS,                 only : calculate_density, calculate_density_derivs, calculate_TFreeze
use MOM_EOS,                 only : calculate_specific_vol_derivs, EOS_domain
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, callTree_showQuery,MOM_mesg
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_version, param_file_type, read_param
use MOM_forcing_type,        only : forcing, MOM_forcing_chksum, find_ustar
use MOM_forcing_type,        only : calculateBuoyancyFlux2d, forcing_SinglePointPrint
use MOM_geothermal,          only : geothermal_entraining, geothermal_in_place
use MOM_geothermal,          only : geothermal_init, geothermal_end, geothermal_CS
use MOM_grid,                only : ocean_grid_type
use MOM_int_tide_input,      only : set_int_tide_input, int_tide_input_init
use MOM_int_tide_input,      only : int_tide_input_end, int_tide_input_CS, int_tide_input_type
use MOM_interface_heights,   only : find_eta, calc_derived_thermo, thickness_to_dz
use MOM_interface_heights,   only : convert_MLD_to_ML_thickness
use MOM_internal_tides,      only : propagate_int_tide, register_int_tide_restarts
use MOM_internal_tides,      only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear,         only : kappa_shear_is_used
use MOM_CVMix_KPP,           only : KPP_CS, KPP_init, KPP_compute_BLD, KPP_calculate
use MOM_CVMix_KPP,           only : KPP_end, KPP_get_BLD, register_KPP_restarts
use MOM_CVMix_KPP,           only : KPP_NonLocalTransport_temp, KPP_NonLocalTransport_saln
use MOM_oda_incupd,          only : apply_oda_incupd, oda_incupd_CS
use MOM_opacity,             only : opacity_init, opacity_end, opacity_CS
use MOM_opacity,             only : absorbRemainingSW, optics_type, optics_nbands
use MOM_open_boundary,       only : ocean_OBC_type
use MOM_regularize_layers,   only : regularize_layers, regularize_layers_init, regularize_layers_CS
use MOM_restart,             only : MOM_restart_CS
use MOM_set_diffusivity,     only : set_diffusivity, set_BBL_TKE
use MOM_set_diffusivity,     only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity,     only : set_diffusivity_CS
use MOM_sponge,              only : apply_sponge, sponge_CS
use MOM_ALE_sponge,          only : apply_ALE_sponge, ALE_sponge_CS
use MOM_time_manager,        only : time_type, real_to_time, operator(-), operator(<=)
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_tracer_diabatic,     only : tracer_vertdiff, tracer_vertdiff_Eulerian
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables,           only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_verticalGrid,        only : verticalGrid_type, get_thickness_units
use MOM_wave_interface,      only : wave_parameters_CS
use MOM_stochastics,         only : stochastic_CS

implicit none ; private

#include <MOM_memory.h>

public diabatic
public diabatic_driver_init
public diabatic_driver_end
public extract_diabatic_member
public adiabatic
public adiabatic_driver_init
public register_diabatic_restarts

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for this module
type, public :: diabatic_CS ; private
  logical :: initialized = .false.   !< True if this control structure has been initialized.

  logical :: use_legacy_diabatic     !< If true (default), use a legacy version of the diabatic
                                     !! algorithm. This is temporary and is needed to avoid change
                                     !! in answers.
  logical :: use_bulkmixedlayer      !< If true, a refined bulk mixed layer is used with
                                     !! nkml sublayers (and additional buffer layers).
  logical :: use_energetic_PBL       !< If true, use the implicit energetics planetary
                                     !! boundary layer scheme to determine the diffusivity
                                     !! in the surface boundary layer.
  logical :: use_KPP                 !< If true, use CVMix/KPP boundary layer scheme to determine the
                                     !! OBLD and the diffusivities within this layer.
  logical :: use_kappa_shear         !< If true, use the kappa_shear module to find the
                                     !! shear-driven diapycnal diffusivity.
  logical :: use_CVMix_shear         !< If true, use the CVMix module to find the
                                     !! shear-driven diapycnal diffusivity.
  logical :: use_CVMix_ddiff         !< If true, use the CVMix double diffusion module.
  logical :: use_CVMix_conv          !< If true, use the CVMix module to get enhanced
                                     !! mixing due to convection.
  logical :: double_diffuse          !< If true, some form of double-diffusive mixing is used.
  logical :: use_sponge              !< If true, sponges may be applied anywhere in the
                                     !! domain.  The exact location and properties of
                                     !! those sponges are set by calls to
                                     !! initialize_sponge and set_up_sponge_field.
  logical :: use_oda_incupd          !< If True, DA incremental update is
                                     !! applied everywhere
  logical :: use_geothermal          !< If true, apply geothermal heating.
  logical :: use_int_tides           !< If true, use the code that advances a separate set
                                     !! of equations for the internal tide energy density.
  logical :: ePBL_is_additive        !< If true, the diffusivity from ePBL is added to all
                                     !! other diffusivities. Otherwise, the larger of kappa-
                                     !! shear and ePBL diffusivities are used.
  real    :: ePBL_Prandtl            !< The Prandtl number used by ePBL to convert vertical
                                     !! diffusivities into viscosities [nondim].
  logical :: useALEalgorithm         !< If true, use the ALE algorithm rather than layered
                                     !! isopycnal/stacked shallow water mode. This logical
                                     !! passed by argument to diabatic_driver_init.
  logical :: aggregate_FW_forcing    !< Determines whether net incoming/outgoing surface
                                     !! FW fluxes are applied separately or combined before
                                     !! being applied.
  real    :: ML_mix_first            !< The nondimensional fraction of the mixed layer
                                     !! algorithm that is applied before diffusive mixing [nondim].
                                     !! The default is 0, while 0.5 gives Strang splitting
                                     !! and 1 is a sensible value too.  Note that if there
                                     !! are convective instabilities in the initial state,
                                     !! the first call may do much more than the second.
  integer :: NKBL                    !< The number of buffer layers (if bulk_mixed_layer)
  logical :: massless_match_targets  !< If true (the default), keep the T & S
                                     !! consistent with the target values.
  logical :: mix_boundary_tracers    !< If true, mix the passive tracers in massless layers at the
                                     !! bottom into the interior as though a diffusivity of
                                     !! Kd_min_tr (see below) were operating.
  logical :: mix_boundary_tracer_ALE !< If true, in ALE mode mix the passive tracers in massless
                                     !! layers at the bottom into the interior as though a
                                     !! diffusivity of Kd_min_tr (see below) were operating.
  real    :: Kd_BBL_tr               !< A bottom boundary layer tracer diffusivity that
                                     !! will allow for explicitly specified bottom fluxes
                                     !! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].  The entrainment at the
                                     !! bottom is at least sqrt(Kd_BBL_tr*dt) over the same distance.
  real    :: Kd_min_tr               !< A minimal diffusivity that should always be
                                     !! applied to tracers, especially in massless layers
                                     !! near the bottom [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: minimum_forcing_depth   !< The smallest depth over which heat and freshwater
                                     !! fluxes are applied [H ~> m or kg m-2].
  real    :: evap_CFL_limit = 0.8    !< The largest fraction of a layer that can be
                                     !! evaporated in one time-step [nondim].
  integer :: halo_TS_diff = 0        !< The temperature, salinity and thickness halo size that
                                     !! must be valid for the diffusivity calculations.
  integer :: halo_diabatic = 0       !< The temperature, salinity, specific volume and thickness
                                     !! halo size that must be valid for the diabatic calculations,
                                     !! including vertical mixing and internal tide propagation.
  logical :: useKPP = .false.        !< use CVMix/KPP diffusivities and non-local transport
  logical :: KPPisPassive            !< If true, KPP is in passive mode, not changing answers.
  logical :: debug                   !< If true, write verbose checksums for debugging purposes.
  logical :: debugConservation       !< If true, monitor conservation and extrema.
  logical :: tracer_tridiag          !< If true, use tracer_vertdiff instead of tridiagTS for
                                     !< vertical diffusion of T and S
  logical :: debug_energy_req        !< If true, test the mixing energy requirement code.
  type(diag_ctrl), pointer :: diag   !< structure used to regulate timing of diagnostic output
  real    :: MLDdensityDifference    !< Density difference used to determine MLD_user [R ~> kg m-3]
  real    :: dz_subML_N2             !< The distance over which to calculate a diagnostic of the
                                     !! average stratification at the base of the mixed layer [Z ~> m].
  real    :: MLD_En_vals(3)          !< Energy values for energy mixed layer diagnostics [R Z3 T-2 ~> J m-2]
  real    :: BMLD_En_vals(3)         !< Energy values for energy bottom mixed layer diagnostics [R Z3 T-2 ~> J m-2]
  logical :: use_OM4_MLD_En_iter     !< If true, uses an older iteration in the energetics MLD calculation to bitwise
                                     !! reproduce OM4 era models
  real    :: ref_h_mld = 0.0         !< The depth of the "surface"  density used in a difference mixed based
                                     !! MLD calculation [Z ~> m].
  logical :: Use_KdWork_diag = .false.  !< Logical flag to indicate if any Kd_work diagnostics are on.
  logical :: Use_N2_diag = .false.   !< Logical flag to indicate if any N2 diagnostics are on.
  logical :: MLD_param_003 = .false. !< Logical flag if MLD in brine plume should use the 0.03 mixed layer depth
  logical :: MLD_param_EN1 = .false. !< Logical flag if MLD in brine plume should use the EN1 mixed layer depth
  logical :: MLD_param_ePBL = .false.!< Logical flag if MLD in brine plume should use the ePBL boundary layer depth

  !>@{ Diagnostic IDs
  integer :: id_ea       = -1, id_eb       = -1 ! used by layer diabatic
  integer :: id_ea_t     = -1, id_eb_t     = -1, id_ea_s   = -1, id_eb_s     = -1
  integer :: id_Kd_heat  = -1, id_Kd_salt  = -1, id_Kd_int = -1, id_Kd_ePBL  = -1
  integer :: id_Tdif     = -1, id_Sdif     = -1, id_Tadv   = -1, id_Sadv     = -1
  integer :: id_N2_dd    = -1, id_N2_salt_dd = -1, id_N2_temp_dd
  ! These are handles to diagnostics related to the mixed layer properties.
  integer :: id_MLD_003 = -1, id_MLD_0125 = -1, id_MLD_user = -1, id_mlotstsq = -1
  integer :: id_MLD_003_zr = -1, id_MLD_003_rr = -1
  integer :: id_MLD_EN1 = -1, id_MLD_EN2  = -1, id_MLD_EN3  = -1, id_subMLN2  = -1
  integer :: id_BMLD_EN1 = -1, id_BMLD_EN2  = -1, id_BMLD_EN3  = -1

  ! These are handles to diagnostics that are only available in non-ALE layered mode.
  integer :: id_wd       = -1
  integer :: id_dudt_dia = -1, id_dvdt_dia = -1
  integer :: id_hf_dudt_dia_2d = -1, id_hf_dvdt_dia_2d = -1

  ! diagnostic for fields prior to applying diapycnal physics
  integer :: id_u_predia = -1, id_v_predia = -1, id_h_predia = -1
  integer :: id_T_predia = -1, id_S_predia = -1, id_e_predia = -1

  integer :: id_diabatic_diff_temp_tend     = -1
  integer :: id_diabatic_diff_saln_tend     = -1
  integer :: id_diabatic_diff_heat_tend     = -1
  integer :: id_diabatic_diff_salt_tend     = -1
  integer :: id_diabatic_diff_heat_tend_2d  = -1
  integer :: id_diabatic_diff_salt_tend_2d  = -1
  integer :: id_diabatic_diff_h = -1

  integer :: id_boundary_forcing_h       = -1
  integer :: id_boundary_forcing_h_tendency   = -1
  integer :: id_boundary_forcing_temp_tend    = -1
  integer :: id_boundary_forcing_saln_tend    = -1
  integer :: id_boundary_forcing_heat_tend    = -1
  integer :: id_boundary_forcing_salt_tend    = -1
  integer :: id_boundary_forcing_heat_tend_2d = -1
  integer :: id_boundary_forcing_salt_tend_2d = -1

  integer :: id_frazil_h    = -1
  integer :: id_frazil_temp_tend    = -1
  integer :: id_frazil_heat_tend    = -1
  integer :: id_frazil_heat_tend_2d = -1
  !>@}

  logical :: diabatic_diff_tendency_diag = .false. !< If true calculate diffusive tendency diagnostics
  logical :: boundary_forcing_tendency_diag = .false. !< If true calculate frazil diagnostics
  logical :: frazil_tendency_diag = .false. !< If true calculate frazil tendency diagnostics

  type(diabatic_aux_CS),        pointer :: diabatic_aux_CSp      => NULL() !< Control structure for a child module
  type(int_tide_input_CS),      pointer :: int_tide_input_CSp    => NULL() !< Control structure for a child module
  type(int_tide_input_type),    pointer :: int_tide_input        => NULL() !< Control structure for a child module
  type(set_diffusivity_CS),     pointer :: set_diff_CSp          => NULL() !< Control structure for a child module
  type(sponge_CS),              pointer :: sponge_CSp            => NULL() !< Control structure for a child module
  type(ALE_sponge_CS),          pointer :: ALE_sponge_CSp        => NULL() !< Control structure for a child module
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp       => NULL() !< Control structure for a child module
  type(optics_type),            pointer :: optics                => NULL() !< Control structure for a child module
  type(KPP_CS),                 pointer :: KPP_CSp               => NULL() !< Control structure for a child module
  type(diapyc_energy_req_CS),   pointer :: diapyc_en_rec_CSp     => NULL() !< Control structure for a child module
  type(oda_incupd_CS),          pointer :: oda_incupd_CSp        => NULL() !< Control structure for a child module
  type(int_tide_CS),            pointer :: int_tide_CSp          => NULL() !< Control structure for a child module
  type(vbf_CS),                 pointer :: VBF                   => NULL() !< Control structure for a child module


  type(bulkmixedlayer_CS) :: bulkmixedlayer         !< Bulk mixed layer control structure
  type(CVMix_conv_CS) :: CVMix_conv                 !< CVMix convection control structure
  type(energetic_PBL_CS) :: ePBL                    !< Energetic PBL control structure
  type(entrain_diffusive_CS) :: entrain_diffusive   !< Diffusive entrainment control structure
  type(geothermal_CS) :: geothermal                 !< Geothermal control structure
  type(opacity_CS) :: opacity                       !< Opacity control structure
  type(regularize_layers_CS) :: regularize_layers   !< Regularize layer control structure

  type(group_pass_type) :: pass_hold_eb_ea !< For group halo pass
  type(diag_grid_storage) :: diag_grids_prev !< Stores diagnostic grids at some previous point in the algorithm

  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
end type diabatic_CS

!>@{ clock ids
integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_differential_diff, id_clock_remap
integer :: id_clock_kpp, id_clock_oda_incupd
!>@}


  interface
module subroutine diabatic(u, v, h, tv, BLD, fluxes, visc, ADp, CDp, dt, Time_end, &
                    G, GV, US, CS, stoch_CS, OBC, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: BLD      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< Structure with vertical viscosities,
                                                                        !! BBL properties and related fields
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< Points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(stochastic_CS),                        pointer       :: stoch_CS !< stochastic control structure
  type(ocean_OBC_type),                       pointer       :: OBC      !< Open boundaries control structure.
  type(Wave_parameters_CS),                   pointer       :: Waves    !< Surface gravity waves

  ! local variables


end subroutine diabatic
module subroutine diabatic_ALE_legacy(u, v, h, tv, BLD, fluxes, visc, ADp, CDp, dt, Time_end, &
                           G, GV, US, CS, stoch_CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: BLD      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< Structure with vertical viscosities,
                                                                        !! BBL properties and related fields
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< Points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(stochastic_CS),                        pointer       :: stoch_CS !< stochastic control structure
  type(Wave_parameters_CS),                   pointer       :: Waves    !< Surface gravity waves

  ! local variables

end subroutine diabatic_ALE_legacy
module subroutine diabatic_ALE(u, v, h, tv, BLD, fluxes, visc, ADp, CDp, dt, Time_end, &
                        G, GV, US, CS, stoch_CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: BLD      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< Structure with vertical viscosities,
                                                                        !! BBL properties and related fields
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< Points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(stochastic_CS),                        pointer       :: stoch_CS !< stochastic control structure
  type(Wave_parameters_CS),                   pointer       :: Waves    !< Surface gravity waves

  ! local variables

end subroutine diabatic_ALE
module subroutine layered_diabatic(u, v, h, tv, BLD, fluxes, visc, ADp, CDp, dt, Time_end, &
                            G, GV, US, CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: BLD      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< Structure with vertical viscosities,
                                                                        !! BBL properties and related fields
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< Points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(Wave_parameters_CS),                   pointer       :: Waves    !< Surface gravity waves

end subroutine layered_diabatic
module subroutine extract_diabatic_member(CS, opacity_CSp, optics_CSp, evap_CFL_limit, minimum_forcing_depth, &
                                   KPP_CSp, energetic_PBL_CSp, diabatic_aux_CSp, diabatic_halo, use_KPP)
  type(diabatic_CS), target, intent(in)      :: CS !< module control structure
  ! All output arguments are optional
  type(opacity_CS),  optional, pointer       :: opacity_CSp !< A pointer to be set to the opacity control structure
  type(optics_type), optional, pointer       :: optics_CSp  !< A pointer to be set to the optics control structure
  type(KPP_CS),      optional, pointer       :: KPP_CSp     !< A pointer to be set to the KPP CS
  type(energetic_PBL_CS), optional, pointer  :: energetic_PBL_CSp !< A pointer to be set to the ePBL CS
  real,              optional, intent(  out) :: evap_CFL_limit !<The largest fraction of a layer that can be
                                                            !! evaporated in one time-step [nondim].
  real,              optional, intent(  out) :: minimum_forcing_depth !< The smallest depth over which heat
                                                            !! and freshwater fluxes are applied [H ~> m or kg m-2].
  type(diabatic_aux_CS), optional, pointer   :: diabatic_aux_CSp !< A pointer to be set to the diabatic_aux
                                                            !! control structure
  integer,           optional, intent(  out) :: diabatic_halo !< The halo size where the diabatic algorithms
                                                            !! assume thermodynamics properties are valid.
  logical,           optional, intent(  out) :: use_KPP       !< If true, diabatic is using KPP vertical mixing

  ! Pointers to control structures
end subroutine extract_diabatic_member
module subroutine adiabatic(h, tv, fluxes, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h      !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv     !< points to thermodynamic fields
  type(forcing),           intent(inout) :: fluxes !< boundary fluxes
  real,                    intent(in)    :: dt     !< time step [T ~> s]
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(diabatic_CS),       pointer       :: CS     !< module control structure


end subroutine adiabatic
module subroutine diagnose_diabatic_diff_tendency(tv, h, temp_old, saln_old, dt, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                      intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h        !< thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: temp_old !< temperature prior to diabatic
                                                                     !! physics [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: saln_old !< salinity prior to diabatic physics [S ~> ppt]
  real,                                       intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),                      intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                          pointer    :: CS       !< module control structure

  ! Local variables

end subroutine diagnose_diabatic_diff_tendency
module subroutine diagnose_boundary_forcing_tendency(tv, h, temp_old, saln_old, h_old, &
                                              dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),   intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h        !< thickness after boundary flux application [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: temp_old !< temperature prior to boundary flux application [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: saln_old !< salinity prior to boundary flux application [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old    !< thickness prior to boundary flux application [H ~> m or kg m-2]
  real,                    intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),   intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),       pointer    :: CS       !< module control structure

  ! Local variables

end subroutine diagnose_boundary_forcing_tendency
module subroutine diagnose_frazil_tendency(tv, h, temp_old, dt, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                     intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h        !< thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: temp_old !< temperature prior to frazil formation [C ~> degC]
  real,                                      intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),                     intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                         pointer    :: CS       !< module control structure


end subroutine diagnose_frazil_tendency
module subroutine adiabatic_driver_init(Time, G, param_file, diag, CS, &
                                tracer_flow_CSp)
  type(time_type),         intent(in)    :: Time             !< current model time
  type(ocean_grid_type),   intent(in)    :: G                !< model grid structure
  type(param_file_type),   intent(in)    :: param_file       !< the file to parse for parameter values
  type(diag_ctrl), target, intent(inout) :: diag             !< regulates diagnostic output
  type(diabatic_CS),       pointer       :: CS               !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp  !< pointer to control structure of the
                                                             !! tracer flow control module

  ! This "include" declares and sets the variable "version".

end subroutine adiabatic_driver_init
module subroutine diabatic_driver_init(Time, G, GV, US, param_file, useALEalgorithm, diag, &
                                ADp, CDp, CS, tracer_flow_CSp, sponge_CSp, &
                                ALE_sponge_CSp, oda_incupd_CSp, int_tide_CSp)
  type(time_type), target                :: Time             !< model time
  type(ocean_grid_type),   intent(inout) :: G                !< model grid structure
  type(verticalGrid_type), intent(in)    :: GV               !< model vertical grid structure
  type(unit_scale_type),   intent(in)    :: US               !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file       !< file to parse for parameter values
  logical,                 intent(in)    :: useALEalgorithm  !< logical for whether to use ALE remapping
  type(diag_ctrl), target, intent(inout) :: diag             !< structure to regulate diagnostic output
  type(accel_diag_ptrs),   intent(inout) :: ADp              !< pointers to accelerations in momentum equations,
                                                             !! to enable diagnostics, like energy budgets
  type(cont_diag_ptrs),    intent(inout) :: CDp              !< pointers to terms in continuity equations
  type(diabatic_CS),       pointer       :: CS               !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp  !< pointer to control structure of the
                                                             !! tracer flow control module
  type(sponge_CS),         pointer       :: sponge_CSp       !< pointer to the sponge module control structure
  type(ALE_sponge_CS),     pointer       :: ALE_sponge_CSp   !< pointer to the ALE sponge module control structure
  type(oda_incupd_CS),     pointer       :: oda_incupd_CSp   !< pointer to the ocean data assimilation incremental
                                                             !! update module control structure
  type(int_tide_CS),       pointer       :: int_tide_CSp     !< pointer to the internal tide structure

  ! Local variables

  ! This "include" declares and sets the variable "version".
end subroutine diabatic_driver_init
module subroutine register_diabatic_restarts(G, GV, US, param_file, int_tide_CSp, restart_CSp, CS)
  type(ocean_grid_type), intent(in)    :: G           !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure
  type(unit_scale_type), intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type), intent(in)    :: param_file  !< A structure to parse for run-time parameters
  type(int_tide_CS),     pointer       :: int_tide_CSp !< Internal tide control structure
  type(MOM_restart_CS),  pointer       :: restart_CSp  !< MOM restart control structure
  type(diabatic_CS),     pointer       :: CS           !< module control structure


end subroutine register_diabatic_restarts
module subroutine diabatic_driver_end(CS)
  type(diabatic_CS), intent(inout) :: CS  !< module control structure

end subroutine diabatic_driver_end
  end interface

end module MOM_diabatic_driver
