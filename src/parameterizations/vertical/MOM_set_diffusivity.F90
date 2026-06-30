! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculate vertical diffusivity from all mixing processes
module MOM_set_diffusivity

use MOM_bkgnd_mixing,        only : calculate_bkgnd_mixing, bkgnd_mixing_init, bkgnd_mixing_cs
use MOM_bkgnd_mixing,        only : bkgnd_mixing_end
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_CVMix_ddiff,         only : CVMix_ddiff_init, CVMix_ddiff_end, CVMix_ddiff_cs
use MOM_CVMix_ddiff,         only : compute_ddiff_coeffs
use MOM_CVMix_shear,         only : calculate_CVMix_shear, CVMix_shear_init, CVMix_shear_cs
use MOM_CVMix_shear,         only : CVMix_shear_end
use MOM_diag_mediator,       only : diag_ctrl, time_type
use MOM_diag_mediator,       only : post_data, register_diag_field
use MOM_diagnose_kdwork,     only : vbf_CS
use MOM_debugging,           only : hchksum, uvchksum, Bchksum, hchksum_pair
use MOM_EOS,                 only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_error_handler,       only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_error_handler,       only : callTree_showQuery
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,        only : forcing, optics_type
use MOM_full_convection,     only : full_convection
use MOM_grid,                only : ocean_grid_type
use MOM_interface_heights,   only : thickness_to_dz, find_rho_bottom
use MOM_internal_tides,      only : int_tide_CS, get_lowmode_loss, get_lowmode_diffusivity
use MOM_intrinsic_functions, only : invcosh
use MOM_io,                  only : slasher, MOM_read_data
use MOM_isopycnal_slopes,    only : vert_fill_TS
use MOM_kappa_shear,         only : calculate_kappa_shear, kappa_shear_init, Kappa_shear_CS
use MOM_kappa_shear,         only : calc_kappa_shear_vertex, kappa_shear_at_vertex
use MOM_open_boundary,       only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary,       only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_string_functions,    only : uppercase
use MOM_tidal_mixing,        only : tidal_mixing_CS, calculate_tidal_mixing, tidal_mixing_h_amp
use MOM_tidal_mixing,        only : setup_tidal_diagnostics, post_tidal_diagnostics
use MOM_tidal_mixing,        only : tidal_mixing_init, tidal_mixing_end
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid,        only : verticalGrid_type
use user_change_diffusivity, only : user_change_diff, user_change_diff_init
use user_change_diffusivity, only : user_change_diff_end, user_change_diff_CS

implicit none ; private

#include <MOM_memory.h>

public set_diffusivity
public set_BBL_TKE
public set_diffusivity_init
public set_diffusivity_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure contains parameters for MOM_set_diffusivity.
type, public :: set_diffusivity_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: debug           !< If true, write verbose checksums for debugging.

  logical :: bulkmixedlayer  !< If true, a refined bulk mixed layer is used with
                             !! GV%nk_rho_varies variable density mixed & buffer layers.
  real    :: FluxRi_max      !< The flux Richardson number where the stratification is
                             !! large enough that N2 > omega2 [nondim].  The full expression
                             !! for the Flux Richardson number is usually
                             !! FLUX_RI_MAX*N2/(N2+OMEGA2). The default is 0.2.
  logical :: bottomdraglaw   !< If true, the  bottom stress is calculated with a
                             !! drag law c_drag*|u|*u.
  logical :: BBL_mixing_as_max !<  If true, take the maximum of the diffusivity
                             !! from the BBL mixing and the other diffusivities.
                             !! Otherwise, diffusivities from the BBL_mixing is added.
  logical :: use_LOTW_BBL_diffusivity !< If true, use simpler/less precise, BBL diffusivity.
  logical :: LOTW_BBL_use_omega !< If true, use simpler/less precise, BBL diffusivity.
  real    :: Von_Karm        !< The von Karman constant as used in the BBL diffusivity calculation
                             !! [nondim].  See (http://en.wikipedia.org/wiki/Von_Karman_constant)
  real    :: BBL_effic       !< Efficiency with which the energy extracted
                             !! by bottom drag drives BBL diffusion in the original BBL scheme, times
                             !! conversion factors between the natural units of mean kinetic energy
                             !! and those those used for TKE [Z2 L-2 ~> nondim].
  real    :: ePBL_BBL_effic  !< efficiency with which the energy extracted
                             !! by bottom drag drives BBL diffusion in the ePBL BBL scheme [nondim]
  logical :: ePBL_BBL_mstar  !< logical if the bottom boundary layer uses an mstar x ustar^3 formulation
                             !! needed here to know whether or not to populate the bottom ustar
  real    :: cdrag           !< quadratic drag coefficient [nondim]
  real    :: dz_BBL_avg_min  !< A minimal distance over which to average to determine the average
                             !! bottom boundary layer density [Z ~> m]
  real    :: IMax_decay      !< Inverse of a maximum decay scale for
                             !! bottom-drag driven turbulence [H-1 ~> m-1 or m2 kg-1].
  real    :: Kv              !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: Kd              !< interior diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kd_min          !< minimum diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kd_max          !< maximum increment for diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                             !! Set to a negative value to have no limit.
  real    :: Kd_add          !< uniform diffusivity added everywhere without
                             !! filtering or scaling [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kd_smooth       !< Vertical diffusivity used to interpolate more
                             !! sensible values of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(diag_ctrl), pointer :: diag => NULL() !< structure to regulate diagnostic output timing

  logical :: limit_dissipation !< If enabled, dissipation is limited to be larger
                               !! than the following:
  real :: dissip_min    !< Minimum dissipation [R Z2 T-3 ~> W m-3]
  real :: dissip_N0     !< Coefficient a in minimum dissipation = a+b*N [R Z2 T-3 ~> W m-3]
  real :: dissip_N1     !< Coefficient b in minimum dissipation = a+b*N [R Z2 T-2 ~> J m-3]
  real :: dissip_N2     !< Coefficient c in minimum dissipation = c*N2 [R Z2 T-1 ~> J s m-3]
  real :: dissip_Kd_min !< Minimum Kd [H Z T-1 ~> m2 s-1 or kg m-1 s-1], with dissipation Rho0*Kd_min*N^2

  real :: omega         !< Earth's rotation frequency [T-1 ~> s-1]
  logical :: ML_radiation !< allow a fraction of TKE available from wind work
                          !! to penetrate below mixed layer base with a vertical
                          !! decay scale determined by the minimum of
                          !! (1) The depth of the mixed layer, or
                          !! (2) An Ekman length scale.
                          !! Energy available to drive mixing below the mixed layer is
                          !! given by E = ML_RAD_COEFF*MSTAR*USTAR**3.  Optionally, if
                          !! ML_rad_TKE_decay is true, this is further reduced by a factor
                          !! of exp(-h_ML*Idecay_len_TkE), where Idecay_len_TKE is
                          !! calculated the same way as in the mixed layer code.
                          !! The diapycnal diffusivity is KD(k) = E/(N2(k)+OMEGA2),
                          !! where N2 is the squared buoyancy frequency [T-2 ~> s-2] and OMEGA2
                          !! is the rotation rate of the earth squared.
  real :: ML_rad_kd_max   !< Maximum diapycnal diffusivity due to turbulence radiated from
                          !! the base of the mixed layer [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: ML_rad_efold_coeff  !< Coefficient to scale penetration depth [nondim]
  real :: ML_rad_coeff        !< Coefficient which scales MSTAR*USTAR^3 to obtain energy
                              !! available for mixing below mixed layer base [nondim]
  logical :: ML_rad_bug       !< If true use code with a bug that reduces the energy available
                              !! in the transition layer by a factor of the inverse of the energy
                              !! deposition lenthscale (in m).
  logical :: ML_rad_TKE_decay !< If true, apply same exponential decay
                              !! to ML_rad as applied to the other surface
                              !! sources of TKE in the mixed layer code.
  real    :: ustar_min        !< A minimum value of ustar to avoid numerical
                              !! problems [Z T-1 ~> m s-1].  If the value is small enough,
                              !! this parameter should not affect the solution.
  real    :: TKE_decay        !< ratio of natural Ekman depth to TKE decay scale [nondim]
  real    :: mstar            !< ratio of friction velocity cubed to
                              !! TKE input to the mixed layer [nondim]
  logical :: ML_use_omega     !< If true, use absolute rotation rate instead
                              !! of the vertical component of rotation when
                              !! setting the decay scale for mixed layer turbulence.
  real    :: ML_omega_frac    !<   When setting the decay scale for turbulence, use
                              !! this fraction [nondim] of the absolute rotation rate blended
                              !! with the local value of f, as f^2 ~= (1-of)*f^2 + of*4*omega^2.
  logical :: user_change_diff !< If true, call user-defined code to change diffusivity.
  logical :: useKappaShear    !< If true, use the kappa_shear module to find the
                              !! shear-driven diapycnal diffusivity.
  logical :: Vertex_Shear     !< If true, do the calculations of the shear-driven mixing
                              !! at the cell vertices (i.e., the vorticity points).
  logical :: use_CVMix_shear  !< If true, use one of the CVMix modules to find
                              !! shear-driven diapycnal diffusivity.
  logical :: double_diffusion !< If true, enable double-diffusive mixing using an old method.
  logical :: use_CVMix_ddiff  !< If true, enable double-diffusive mixing via CVMix.
  logical :: use_tidal_mixing !< If true, activate tidal mixing diffusivity.
  logical :: use_int_tides    !< If true, use internal tides ray tracing
  logical :: simple_TKE_to_Kd !< If true, uses a simple estimate of Kd/TKE that
                              !! does not rely on a layer-formulation.
  real    :: Max_Rrho_salt_fingers      !< max density ratio for salt fingering [nondim]
  real    :: Max_salt_diff_salt_fingers !< max salt diffusivity for salt fingers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kv_molecular     !< Molecular viscosity for double diffusive convection [H Z T-1 ~> m2 s-1 or Pa s]

  integer :: answer_date      !< The vintage of the order of arithmetic and expressions in this module's
                              !! calculations.  Values below 20190101 recover the answers from the
                              !! end of 2018, while higher values use updated and more robust forms
                              !! of the same expressions.  Values above 20240630 use more accurate
                              !! expressions for cases where USE_LOTW_BBL_DIFFUSIVITY is true.  Values
                              !! above 20250301 use less confusing expressions to set the bottom-drag
                              !! generated diffusivity when USE_LOTW_BBL_DIFFUSIVITY is false.
  integer :: LOTW_BBL_answer_date !< The vintage of the order of arithmetic and expressions
                              !! in the LOTW_BBL calculations.  Values below 20240630 recover the
                              !! original answers, while higher values use more accurate expressions.
                              !! This only applies when USE_LOTW_BBL_DIFFUSIVITY is true.
  integer :: drag_diff_answer_date !< The vintage of the order of arithmetic in the drag diffusivity
                              !! calculations.  Values above 20250301 use less confusing expressions
                              !! to set the bottom-drag generated diffusivity when
                              !! USE_LOTW_BBL_DIFFUSIVITY is false.

  character(len=200) :: inputdir !< The directory in which input files are found
  type(user_change_diff_CS), pointer :: user_change_diff_CSp => NULL() !< Control structure for a child module
  type(Kappa_shear_CS),      pointer :: kappaShear_CSp       => NULL() !< Control structure for a child module
  type(CVMix_shear_cs),      pointer :: CVMix_shear_csp      => NULL() !< Control structure for a child module
  type(CVMix_ddiff_cs),      pointer :: CVMix_ddiff_csp      => NULL() !< Control structure for a child module
  type(bkgnd_mixing_cs),     pointer :: bkgnd_mixing_csp     => NULL() !< Control structure for a child module
  type(int_tide_CS),         pointer :: int_tide_CSp         => NULL() !< Control structure for a child module
  type(tidal_mixing_cs) :: tidal_mixing   !< Control structure for a child module

  !>@{ Diagnostic IDs
  integer :: id_maxTKE     = -1, id_TKE_to_Kd   = -1, id_Kd_user    = -1
  integer :: id_Kd_layer   = -1, id_Kd_BBL      = -1, id_N2         = -1
  integer :: id_Kd_Work    = -1, id_KT_extra    = -1, id_KS_extra   = -1, id_R_rho    = -1
  integer :: id_Kd_bkgnd   = -1, id_Kv_bkgnd    = -1, id_Kd_leak    = -1
  integer :: id_Kd_quad    = -1, id_Kd_itidal   = -1, id_Kd_Froude  = -1, id_Kd_slope = -1
  integer :: id_prof_leak  = -1, id_prof_quad   = -1, id_prof_itidal= -1
  integer :: id_prof_Froude= -1, id_prof_slope  = -1, id_bbl_thick = -1, id_kbbl = -1
  integer :: id_Kd_Work_added = -1
  !>@}

end type set_diffusivity_CS

!> This structure has memory for used in calculating diagnostics of diffusivity
type diffusivity_diags
  real, pointer, dimension(:,:,:) :: &
    N2_3d     => NULL(), & !< squared buoyancy frequency at interfaces [T-2 ~> s-2]
    Kd_user   => NULL(), & !< user-added diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_BBL    => NULL(), & !< BBL diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Work   => NULL(), & !< layer integrated work by diapycnal mixing [R Z3 T-3 ~> W m-2]
    Kd_Work_added   => NULL(), & !< layer integrated work by added mixing [R Z3 T-3 ~> W m-2]
    maxTKE    => NULL(), & !< energy required to entrain to h_max [H Z2 T-3 ~> m3 s-3 or W m-2]
    Kd_bkgnd  => NULL(), & !< Background diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kv_bkgnd  => NULL(), & !< Viscosity from background diffusivity at interfaces [H Z T-1 ~> m2 s-1 or Pa s]
    KT_extra  => NULL(), & !< Double diffusion diffusivity for temperature [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    KS_extra  => NULL(), & !< Double diffusion diffusivity for salinity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    drho_rat  => NULL(), & !< The density difference ratio used in double diffusion [nondim].
    Kd_leak   => NULL(), & !< internal tides leakage diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_quad   => NULL(), & !< internal tides bottom drag diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_itidal => NULL(), & !< internal tides wave drag diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Froude => NULL(), & !< internal tides high Froude diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_slope  => NULL(), & !< internal tides critical slopes diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    prof_leak   => NULL(), & !< vertical profile for leakage [H-1 ~> m-1 or m2 kg-1]
    prof_quad   => NULL(), & !< vertical profile for bottom drag [H-1 ~> m-1 or m2 kg-1]
    prof_itidal => NULL(), & !< vertical profile for wave drag [H-1 ~> m-1 or m2 kg-1]
    prof_Froude => NULL(), & !< vertical profile for Froude drag [H-1 ~> m-1 or m2 kg-1]
    prof_slope  => NULL()    !< vertical profile for critical slopes [H-1 ~> m-1 or m2 kg-1]
  real, pointer, dimension(:,:) ::  bbl_thick => NULL(), & !< bottom boundary layer thickness [H ~> m or kg m-2]
                                    kbbl => NULL() !< top of bottom boundary layer

  real, pointer, dimension(:,:,:) :: TKE_to_Kd => NULL()
                          !< conversion rate (~1.0 / (G_Earth + dRho_lay)) between TKE
                          !! dissipated within a layer and Kd in that layer [T2 Z-1 ~> s2 m-1]

end type diffusivity_diags

!>@{ CPU time clocks
integer :: id_clock_kappaShear, id_clock_CVMix_ddiff
!>@}


  interface
module subroutine set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, optics, visc, dt, Kd_int, &
                           G, GV, US, CS, VBF, Kd_lay, Kd_extra_T, Kd_extra_S)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: u_h  !< Zonal velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: v_h  !< Meridional velocity interpolated to h points [L T-1 ~> m s-1].
  type(thermo_var_ptrs),     intent(inout) :: tv   !< Structure with pointers to thermodynamic
                                                   !! fields. Out is for tv%TempxPmE.
  type(forcing),             intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(optics_type),         pointer       :: optics !< A structure describing the optical
                                                   !!  properties of the ocean.
  type(vertvisc_type),       intent(inout) :: visc !< Structure containing vertical viscosities, bottom
                                                   !! boundary layer properties and related fields.
  real,                      intent(in)    :: dt   !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                             intent(out)   :: Kd_int !< Diapycnal diffusivity at each interface
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(set_diffusivity_CS),  pointer       :: CS   !< Module control structure.
  type(vbf_CS),           pointer          :: VBF  !< A diagnostic control structure for vertical buoyancy fluxes
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   optional, intent(out)   :: Kd_lay !< Diapycnal diffusivity of each layer
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                   optional, intent(out)   :: Kd_extra_T !< The extra diffusivity at interfaces of
                                                   !! temperature due to double diffusion relative
                                                   !! to the diffusivity of density
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                   optional, intent(out)   :: Kd_extra_S !< The extra diffusivity at interfaces of
                                                   !! salinity due to double diffusion relative
                                                   !! to the diffusivity of density
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! local variables



                  ! filled vertically by diffusion or the properties after full convective adjustment.

                  !< TKE dissipated within a layer and Kd in that layer [T2 Z-1 ~> s2 m-1]



                            ! buffer layer, or -1 without a bulk mixed layer.



end subroutine set_diffusivity
module subroutine find_TKE_to_Kd(h, tv, dRho_int, N2_lay, j, dt, G, GV, US, CS, &
                          TKE_to_Kd, maxTKE, kb)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv   !< Structure containing pointers to any available
                                                          !! thermodynamic fields.
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: dRho_int !< Change in locally referenced potential density
                                                          !! across each interface [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: N2_lay !< The squared buoyancy frequency of the
                                                          !! layers [T-2 ~> s-2].
  integer,                          intent(in)    :: j    !< j-index of row to work on
  real,                             intent(in)    :: dt   !< Time increment [T ~> s].
  type(set_diffusivity_CS),         pointer       :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: TKE_to_Kd !< The conversion rate between the
                                                          !! TKE dissipated within a layer and the
                                                          !! diapycnal diffusivity within that layer,
                                                          !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                          !! [T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: maxTKE !< The energy required to for a layer to entrain to its
                                                          !! maximum realizable thickness [H Z2 T-3 ~> m3 s-3 or W m-2]
  integer, dimension(SZI_(G)),      intent(out)   :: kb   !< Index of lightest layer denser than the buffer
                                                          !! layer, or -1 without a bulk mixed layer.
  ! Local variables
                  ! below it [nondim]
end subroutine find_TKE_to_Kd
module subroutine find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, US, CS, dRho_int, &
                   N2_lay, N2_int, N2_bot, Rho_bot, h_bot, k_bot)
  type(ocean_grid_type),    intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(in)  :: tv   !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: T_f  !< layer temperature with the values in massless layers
                                                !! filled vertically by diffusion [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: S_f  !< Layer salinities with values in massless
                                                !! layers filled vertically by diffusion [S ~> ppt].
  type(forcing),            intent(in)  :: fluxes !< A structure of thermodynamic surface fluxes
  integer,                  intent(in)  :: j    !< j-index of row to work on
  type(set_diffusivity_CS), pointer     :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(out) :: dRho_int !< Change in locally referenced potential density
                                                !! across each interface [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(out) :: N2_int !< The squared buoyancy frequency at the interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(out) :: N2_lay !< The squared buoyancy frequency of the layers [T-2 ~> s-2].
  real, dimension(SZI_(G)), intent(out) :: N2_bot !< The near-bottom squared buoyancy frequency [T-2 ~> s-2].
  real, dimension(SZI_(G)), intent(out) :: Rho_bot !< Near-bottom density [R ~> kg m-3].
  real, dimension(SZI_(G)), optional, intent(out) :: h_bot !< Bottom boundary layer thickness [H ~> m or kg m-2].
  integer, dimension(SZI_(G)), optional, intent(out) :: k_bot !< Bottom boundary layer top layer index.

  ! Local variables

                    ! times some unit conversion factors [H T-2 R-1 ~> m4 s-2 kg-1 or m s-2].


end subroutine find_N2
module subroutine double_diffusion(tv, h, T_f, S_f, j, G, GV, US, CS, Kd_T_dd, Kd_S_dd)
  type(ocean_grid_type),    intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)  :: US  !< A dimensional unit scaling type
  type(thermo_var_ptrs),    intent(in)  :: tv  !< Structure containing pointers to any available
                                               !! thermodynamic fields; absent fields have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: T_f !< layer temperatures with the values in massless layers
                                               !! filled vertically by diffusion [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: S_f !< Layer salinities with values in massless
                                               !! layers filled vertically by diffusion [S ~> ppt].
  integer,                  intent(in)  :: j   !< Meridional index upon which to work.
  type(set_diffusivity_CS), pointer     :: CS  !< Module control structure.
  real, dimension(SZI_(G),SZK_(GV)+1),       &
                            intent(out) :: Kd_T_dd !< Interface double diffusion diapycnal
                                               !! diffusivity for temp [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1),       &
                            intent(out) :: Kd_S_dd !< Interface double diffusion diapycnal
                                               !! diffusivity for saln [H Z T-1 ~> m2 s-1 or kg m-1 s-1]




end subroutine double_diffusion
module subroutine add_drag_diffusivity(h, u, v, tv, fluxes, visc, j, TKE_to_Kd, maxTKE, &
                                kb, rho_bot, G, GV, US, CS, Kd_lay, Kd_int, Kd_BBL)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv   !< Structure containing pointers to any available
                                                          !! thermodynamic fields.
  type(forcing),                    intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(vertvisc_type),              intent(in)    :: visc !< Structure containing vertical viscosities, bottom
                                                          !! boundary layer properties and related fields
  integer,                          intent(in)    :: j    !< j-index of row to work on
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                          !! TKE dissipated within a layer and the
                                                          !! diapycnal diffusivity within that layer,
                                                          !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                          !! [T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: maxTKE !< The energy required to for a layer to entrain to its
                                                          !! maximum-realizable thickness [H Z2 T-3 ~> m3 s-3 or W m-2]
  integer, dimension(SZI_(G)),      intent(in)    :: kb   !< Index of lightest layer denser than the buffer
                                                          !! layer, or -1 without a bulk mixed layer
  real, dimension(SZI_(G)),         intent(in)    :: rho_bot !< In situ density averaged over a near-bottom
                                                          !! region [R ~> kg m-3]
  type(set_diffusivity_CS),         pointer       :: CS   !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers,
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1), intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces,
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(:,:,:),           pointer       :: Kd_BBL !< Interface BBL diffusivity
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

! This routine adds diffusion sustained by flow energy extracted by bottom drag.

end subroutine add_drag_diffusivity
module subroutine add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, Rho_bot, Kd_int, &
                                    G, GV, US, CS, Kd_BBL, Kd_lay)
  type(ocean_grid_type),    intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u  !< u component of flow [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v  !< v component of flow [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(in)    :: tv !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  type(forcing),            intent(in)    :: fluxes !< Surface fluxes structure
  type(vertvisc_type),      intent(in)    :: visc !< Structure containing vertical viscosities, bottom
                                                  !! boundary layer properties and related fields.
  integer,                  intent(in)    :: j  !< j-index of row to work on
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(in)    :: N2_int !< Square of Brunt-Vaisala at interfaces [T-2 ~> s-2]
  real, dimension(SZI_(G)), intent(in)    :: rho_bot !< In situ density averaged over a near-bottom
                                                     !! region [R ~> kg m-3]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(inout) :: Kd_int !< Interface net diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(set_diffusivity_CS), pointer       :: CS !< Diffusivity control structure
  real, dimension(:,:,:),   pointer       :: Kd_BBL !< Interface BBL diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)), &
                  optional, intent(inout) :: Kd_lay !< Layer net diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! Local variables
                           ! can act as a source of TKE [H L2 T-3 ~> m3 s-3 or W m-2]
                           ! height [H-1 ~> m-1 or m2 kg-1].
                           ! the assumption that this extracted energy also drives diapycnal mixing.

end subroutine add_LOTW_BBL_diffusivity
module subroutine add_MLrad_diffusivity(dz, fluxes, tv, j, Kd_int, G, GV, US, CS, TKE_to_Kd, Kd_lay)
  type(ocean_grid_type),            intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: dz     !< Height change across layers [Z ~> m]
  type(forcing),                    intent(in)    :: fluxes !< Surface fluxes structure
  type(thermo_var_ptrs),            intent(in)    :: tv     !< Structure containing pointers to any available
                                                            !! thermodynamic fields.
  integer,                          intent(in)    :: j      !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)+1), intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(set_diffusivity_CS),         pointer       :: CS     !< Diffusivity control structure
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                            !! TKE dissipated within  a layer and the
                                                            !! diapycnal diffusivity witin that layer,
                                                            !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                            !! [T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].

! This routine adds effects of mixed layer radiation to the layer diffusivities.

                                        ! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

                            ! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                            ! layer code [Z-2 ~> m-2]

end subroutine add_MLrad_diffusivity
module subroutine set_BBL_TKE(u, v, h, tv, fluxes, visc, G, GV, US, CS, OBC)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),    intent(in)    :: tv   !< Structure with pointers to thermodynamic fields
  type(forcing),            intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(vertvisc_type),      intent(inout) :: visc !< Structure containing vertical viscosities, bottom
                                                  !! boundary layer properties and related fields.
  type(set_diffusivity_CS), pointer       :: CS   !< Diffusivity control structure
  type(ocean_OBC_type),     pointer       :: OBC  !< Open boundaries control structure.

  ! This subroutine calculates several properties related to bottom
  ! boundary layer turbulence.







end subroutine set_BBL_TKE
module subroutine set_density_ratios(h, tv, kb, G, GV, US, CS, j, ds_dsp1, rho_0)
  type(ocean_grid_type),            intent(in)   :: G  !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)   :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)   :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),            intent(in)   :: tv !< Structure containing pointers to any
                                                       !! available thermodynamic fields; absent
                                                       !! fields have NULL ptrs.
  integer, dimension(SZI_(G)),      intent(in)   :: kb !< Index of lightest layer denser than the buffer
                                                       !! layer, or -1 without a bulk mixed layer.
  type(unit_scale_type),            intent(in)   :: US !< A dimensional unit scaling type
  type(set_diffusivity_CS),         pointer      :: CS !< Control structure returned by previous
                                                       !! call to diabatic_entrain_init.
  integer,                          intent(in)   :: j  !< Meridional index upon which to work.
  real, dimension(SZI_(G),SZK_(GV)), intent(out) :: ds_dsp1 !< Coordinate variable (sigma-2)
                                                       !! difference across an interface divided by
                                                       !! the difference across the interface below
                                                       !! it [nondim]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(in)   :: rho_0 !< Layer potential densities relative to
                                                       !! surface press [R ~> kg m-3].

  ! Local variables
                                   ! layers [R-1 ~> m3 kg-1]

end subroutine set_density_ratios
module subroutine set_diffusivity_init(Time, G, GV, US, param_file, diag, CS, int_tide_CSp, halo_TS, &
                                double_diffuse, physical_OBL_scheme)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time
                                                  !! parameters.
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure used to regulate diagnostic output.
  type(set_diffusivity_CS), pointer       :: CS   !< pointer set to point to the module control
                                                  !! structure.
  type(int_tide_CS),        pointer       :: int_tide_CSp !< Internal tide control structure
  integer,                  intent(out)   :: halo_TS !< The halo size of tracer points that must be
                                                  !! valid for the calculations in set_diffusivity.
  logical,                  intent(out)   :: double_diffuse !< This indicates whether some version
                                                  !! of double diffusion is being used.
  logical,                  intent(in)    :: physical_OBL_scheme !< If true, a physically based
                                                  !! parameterization (like KPP or ePBL or a bulk mixed
                                                  !! layer) is used outside of set_diffusivity to
                                                  !! specify the mixing that occurs in the ocean's
                                                  !! surface boundary layer.

  ! Local variables
  ! This include declares and sets the variable "version".
                             ! in setting the default for other diffusivities.
                             ! that is used in place of the absolute value of the local Coriolis
                             ! parameter in the denominator of some expressions [nondim]
                                     ! the Bryan-Lewis (1979) style tanh profile.
                             ! isopycnal or stacked shallow water mode.

end subroutine set_diffusivity_init
module subroutine set_diffusivity_end(CS)
  type(set_diffusivity_CS), intent(inout) :: CS !< Control structure for this module

end subroutine set_diffusivity_end
  end interface

end module MOM_set_diffusivity
