! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface to vertical tidal mixing schemes including CVMix tidal mixing.
module MOM_tidal_mixing

use MOM_diag_mediator,      only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,      only : safe_alloc_ptr, post_data
use MOM_diagnose_Kdwork,    only : vbf_CS
use MOM_debugging,          only : hchksum
use MOM_error_handler,      only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,        only : openParameterBlock, closeParameterBlock
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : slasher, MOM_read_data, field_size
use MOM_io,                 only : read_netCDF_data
use MOM_internal_tides,     only : int_tide_CS, get_lowmode_loss
use MOM_remapping,          only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_string_functions,   only : uppercase, lowercase
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : thermo_var_ptrs, p3d
use MOM_verticalGrid,       only : verticalGrid_type
use CVMix_tidal,            only : CVMix_init_tidal, CVMix_compute_Simmons_invariant
use CVMix_tidal,            only : CVMix_coeffs_tidal, CVMix_tidal_params_type
use CVMix_tidal,            only : CVMix_compute_Schmittner_invariant, CVMix_compute_SchmittnerCoeff
use CVMix_tidal,            only : CVMix_coeffs_tidal_schmittner
use CVMix_kinds_and_types,  only : CVMix_global_params_type
use CVMix_put_get,          only : CVMix_put

implicit none ; private

#include <MOM_memory.h>

public tidal_mixing_init
public setup_tidal_diagnostics
public calculate_tidal_mixing
public post_tidal_diagnostics
public tidal_mixing_h_amp
public tidal_mixing_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Containers for tidal mixing diagnostics
type, public :: tidal_mixing_diags ; private
  real, allocatable :: Kd_itidal(:,:,:)       !< internal tide diffusivity at interfaces
                                              !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, allocatable :: Fl_itidal(:,:,:)       !< vertical flux of tidal turbulent dissipation
                                              !! [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable :: Kd_Niku(:,:,:)         !< lee-wave diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, allocatable :: Kd_Niku_work(:,:,:)    !< layer integrated work by lee-wave driven mixing [R Z3 T-3 ~> W m-2]
  real, allocatable :: Kd_Itidal_Work(:,:,:)  !< layer integrated work by int tide driven mixing [R Z3 T-3 ~> W m-2]
  real, allocatable :: Kd_Lowmode_Work(:,:,:) !< layer integrated work by low mode driven mixing [R Z3 T-3 ~> W m-2]
  real, allocatable :: N2_int(:,:,:)          !< Buoyancy frequency squared at interfaces [T-2 ~> s-2]
  real, allocatable :: vert_dep_3d(:,:,:)     !< The 3-d mixing energy deposition vertical fraction [nondim]?
  real, allocatable :: Schmittner_coeff_3d(:,:,:) !< The coefficient in the Schmittner et al mixing scheme [nondim]
  real, allocatable :: tidal_qe_md(:,:,:)     !< Input tidal energy dissipated locally,
                                              !! interpolated to model vertical coordinate [R Z3 T-3 ~> W m-2]
  real, allocatable :: Kd_lowmode(:,:,:)      !< internal tide diffusivity at interfaces
                                              !! due to propagating low modes [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, allocatable :: Fl_lowmode(:,:,:)      !< vertical flux of tidal turbulent
                                              !! dissipation due to propagating low modes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable :: TKE_itidal_used(:,:)   !< internal tide TKE input at ocean bottom [R Z3 T-3 ~> W m-2]
  real, allocatable :: N2_bot(:,:)            !< bottom squared buoyancy frequency [T-2 ~> s-2]
  real, allocatable :: N2_meanz(:,:)          !< vertically averaged buoyancy frequency [T-2 ~> s-2]
  real, allocatable :: Polzin_decay_scale_scaled(:,:) !< Vertical scale of decay for tidal dissipation [Z ~> m]
  real, allocatable :: Polzin_decay_scale(:,:)  !< Vertical decay scale for tidal dissipation with Polzin [Z ~> m]
  real, allocatable :: Simmons_coeff_2d(:,:)  !< The Simmons et al mixing coefficient [nondim]
end type

!> Control structure with parameters for the tidal mixing module.
type, public :: tidal_mixing_cs ; private
  logical :: debug = .true.   !< If true, do more extensive debugging checks.  This is hard-coded.

  ! Parameters
  logical :: int_tide_dissipation = .false. !< Internal tide conversion (from barotropic)
                              !! with the schemes of St Laurent et al (2002) & Simmons et al (2004)

  integer :: Int_tide_profile !< A coded integer indicating the vertical profile
                              !! for dissipation of the internal waves.  Schemes that are
                              !! currently encoded are St Laurent et al (2002) and Polzin (2009).
  logical :: Lee_wave_dissipation = .false. !< Enable lee-wave driven mixing, following
                              !! Nikurashin (2010), with a vertical energy
                              !! deposition profile specified by Lee_wave_profile to be
                              !! St Laurent et al (2002) or Simmons et al (2004) scheme

  integer :: Lee_wave_profile !< A coded integer indicating the vertical profile
                              !! for dissipation of the lee waves.  Schemes that are
                              !! currently encoded are St Laurent et al (2002) and
                              !! Polzin (2009).
  real :: Int_tide_decay_scale !< decay scale for internal wave TKE [Z ~> m]

  real :: Mu_itides           !< efficiency for conversion of dissipation
                              !! to potential energy [nondim]

  real :: Gamma_itides        !< fraction of local dissipation [nondim]

  real :: Gamma_lee           !< fraction of local dissipation for lee waves
                              !! (Nikurashin's energy input) [nondim]
  real :: Decay_scale_factor_lee !< Scaling factor for the decay scale of lee
                              !! wave energy dissipation [nondim]

  real :: min_zbot_itides     !< minimum depth for internal tide conversion [Z ~> m].
  logical :: Lowmode_itidal_dissipation = .false.  !< If true, consider mixing due to breaking low
                              !! modes that have been remotely generated using an internal tidal
                              !! dissipation scheme to specify the vertical profile of the energy
                              !! input to drive diapycnal mixing, along the lines of St. Laurent
                              !! et al. (2002) and Simmons et al. (2004).

  real :: Nu_Polzin           !< The non-dimensional constant used in Polzin form of
                              !! the vertical scale of decay of tidal dissipation [nondim]

  real :: Nbotref_Polzin      !< Reference value for the buoyancy frequency at the
                              !! ocean bottom used in Polzin formulation of the
                              !! vertical scale of decay of tidal dissipation [T-1 ~> s-1]
  real :: Polzin_decay_scale_factor !< Scaling factor for the decay length scale
                              !! of the tidal dissipation profile in Polzin [nondim]
  real :: Polzin_decay_scale_max_factor !< The decay length scale of tidal dissipation
                              !! profile in Polzin formulation should not exceed
                              !! Polzin_decay_scale_max_factor * depth of the ocean [nondim].
  real :: Polzin_min_decay_scale !< minimum decay scale of the tidal dissipation
                              !! profile in Polzin formulation [Z ~> m]

  real :: TKE_itide_max       !< maximum internal tide conversion [R Z3 T-3 ~> W m-2]
                              !! available to mix above the BBL

  real :: utide               !< constant tidal amplitude [Z T-1 ~> m s-1] if READ_TIDEAMP is false.
  real :: kappa_itides        !< topographic wavenumber and non-dimensional scaling [Z-1 ~> m-1].
  real :: kappa_h2_factor     !< factor for the product of wavenumber * rms sgs height [nondim]
  character(len=200) :: inputdir !< The directory in which to find input files

  logical :: use_CVMix_tidal = .false. !< true if CVMix is to be used for determining
                              !! diffusivity due to tidal mixing

  real :: min_thickness       !< Minimum thickness allowed [Z ~> m]

  ! CVMix-specific parameters
  integer                         :: CVMix_tidal_scheme = -1  !< 1 for Simmons, 2 for Schmittner
  type(CVMix_tidal_params_type)   :: CVMix_tidal_params !< A CVMix-specific type with parameters for tidal mixing
  type(CVMix_global_params_type)  :: CVMix_glb_params   !< CVMix-specific for Prandtl number only
  real                            :: tidal_max_coef     !< CVMix-specific maximum allowable tidal
                                                        !! diffusivity. [Z2 T-1 ~> m2 s-1]
  real                            :: tidal_diss_lim_tc  !< CVMix-specific dissipation limit depth for
                                                        !! tidal-energy-constituent data [Z ~> m].
  type(remapping_CS)              :: remap_CS           !< The control structure for remapping
  integer :: remap_answer_date  !< The vintage of the order of arithmetic and expressions to use
                                !! for remapping.  Values below 20190101 recover the remapping
                                !! answers from 2018, while higher values use more robust
                                !! forms of the same remapping expressions.
  integer :: tidal_answer_date  !< The vintage of the order of arithmetic and expressions in the tidal
                                !! mixing calculations.  Values below 20190101 recover the answers
                                !! from the end of 2018, while higher values use updated and more robust
                                !! forms of the same expressions.

  type(int_tide_CS), pointer    :: int_tide_CSp=> NULL() !< Control structure for a child module

  ! Data containers
  real, allocatable :: TKE_Niku(:,:)    !< Lee wave driven Turbulent Kinetic Energy input
                                        !! [R Z3 T-3 ~> W m-2]
  real, allocatable :: TKE_itidal(:,:)  !< The internal Turbulent Kinetic Energy input divided by
                                        !! the bottom stratification and in non-Boussinesq mode by
                                        !! the near-bottom density [R Z4 H-1 T-2 ~> J m-2 or J m kg-1]
  real, allocatable :: Nb(:,:)          !< The near bottom buoyancy frequency [T-1 ~> s-1].
  real, allocatable :: mask_itidal(:,:) !< A mask of where internal tide energy is input [nondim]
  real, allocatable :: h2(:,:)          !< Squared bottom depth variance [Z2 ~> m2].
  real, allocatable :: tideamp(:,:)     !< RMS tidal amplitude [Z T-1 ~> m s-1]
  real, allocatable :: h_src(:)         !< tidal constituent input layer thickness [m]
  real, allocatable :: tidal_qe_2d(:,:) !< Tidal energy input times the local dissipation
                                        !! fraction, q*E(x,y), with the CVMix implementation
                                        !! of Jayne et al tidal mixing [R Z3 T-3 ~> W m-2].
                                        !! TODO: make this E(x,y) only
  real, allocatable :: tidal_qe_3d_in(:,:,:) !< q*E(x,y,z) with the Schmittner parameterization [R Z3 T-3 ~> W m-2]


  ! Diagnostics
  type(diag_ctrl),          pointer :: diag => NULL() !< structure to regulate diagnostic output timing
  type(tidal_mixing_diags) :: dd        !< Tidal mixing diagnostic arrays

  !>@{ Diagnostic identifiers
  integer :: id_TKE_itidal                = -1
  integer :: id_TKE_leewave               = -1
  integer :: id_Kd_itidal                 = -1
  integer :: id_Kd_Niku                   = -1
  integer :: id_Kd_lowmode                = -1
  integer :: id_Kd_Itidal_Work            = -1
  integer :: id_Kd_Niku_Work              = -1
  integer :: id_Kd_Lowmode_Work           = -1
  integer :: id_Nb                        = -1
  integer :: id_N2_bot                    = -1
  integer :: id_N2_meanz                  = -1
  integer :: id_Fl_itidal                 = -1
  integer :: id_Fl_lowmode                = -1
  integer :: id_Polzin_decay_scale        = -1
  integer :: id_Polzin_decay_scale_scaled = -1
  integer :: id_N2_int                    = -1
  integer :: id_Simmons_coeff             = -1
  integer :: id_Schmittner_coeff          = -1
  integer :: id_tidal_qe_md               = -1
  integer :: id_vert_dep                  = -1
  !>@}

end type tidal_mixing_cs

!>@{ Coded parmameters for specifying mixing schemes
character*(20), parameter :: STLAURENT_PROFILE_STRING   = "STLAURENT_02"
character*(20), parameter :: POLZIN_PROFILE_STRING      = "POLZIN_09"
integer,        parameter :: STLAURENT_02 = 1
integer,        parameter :: POLZIN_09    = 2
character*(20), parameter :: SIMMONS_SCHEME_STRING      = "SIMMONS"
character*(20), parameter :: SCHMITTNER_SCHEME_STRING   = "SCHMITTNER"
integer,        parameter :: SIMMONS   = 1
integer,        parameter :: SCHMITTNER   = 2
!>@}


  interface
logical module function tidal_mixing_init(Time, G, GV, US, param_file, int_tide_CSp, diag, CS)
  type(time_type),          intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),    intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type),  intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< Run-time parameter file handle
  type(int_tide_CS),        pointer       :: int_tide_CSp !< A pointer to the internal tides control structure
  type(diag_ctrl), target,  intent(inout) :: diag       !< Diagnostics control structure.
  type(tidal_mixing_cs),    intent(inout) :: CS         !< This module's control structure.

  ! Local variables
                        ! diffusivities into viscosities [nondim]

  ! This include declares and sets the variable "version".

end function tidal_mixing_init
module subroutine calculate_tidal_mixing(dz, j, N2_bot, Rho_bot, N2_lay, N2_int, TKE_to_Kd, max_TKE, &
                                  G, GV, US, CS, Kd_max, Kv, Kd_lay, Kd_int, VBF)
  type(ocean_grid_type),            intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: dz     !< The vertical distance across layers [Z ~> m]
  integer,                          intent(in)    :: j      !< The j-index to work on
  real, dimension(SZI_(G)),         intent(in)    :: N2_bot !< The near-bottom squared buoyancy
                                                            !! frequency [T-2 ~> s-2].
  real, dimension(SZI_(G)),         intent(in)    :: Rho_bot !< The near-bottom in situ density [R ~> kg m-3]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: N2_lay !< The squared buoyancy frequency of the
                                                            !! layers [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: N2_int !< The squared buoyancy frequency at the
                                                            !! interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                            !! dissipated within a layer and the
                                                            !! diapycnal diffusivity within that layer,
                                                            !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                            !! [T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: max_TKE !< The energy required for a layer to
                                                            !! entrain to its maximum realizable
                                                            !! thickness [H Z2 T-3 ~> m3 s-3 or W m-2]
  type(tidal_mixing_cs),            intent(inout) :: CS     !< The control structure for this module
  real,                             intent(in)    :: Kd_max !< The maximum increment for diapycnal
                                                            !! diffusivity due to TKE-based processes,
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                                                            !! Set this to a negative value to have no limit.
  real, dimension(:,:,:),           pointer       :: Kv     !< The "slow" vertical viscosity at each interface
                                                            !! (not layer!) [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                          optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces
                                                            !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(vbf_CS), pointer                           :: VBF    !< A diagnostic structure for vertical buoyancy fluxes

end subroutine calculate_tidal_mixing
module subroutine calculate_CVMix_tidal(dz, j, N2_int, G, GV, US, CS, Kv, Kd_lay, Kd_int)
  type(ocean_grid_type),   intent(in)    :: G     !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(tidal_mixing_cs),   intent(inout) :: CS    !< This module's control structure.
  real, dimension(SZI_(G),SZK_(GV)),   intent(in) :: dz     !< The vertical distance across layers [Z ~> m]
  integer,                 intent(in)    :: j     !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: N2_int !< The squared buoyancy
                                                  !! frequency at the interfaces [T-2 ~> s-2].
  real, dimension(:,:,:),  pointer       :: Kv    !< The "slow" vertical viscosity at each interface
                                                  !! (not layer!) [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZI_(G),SZK_(GV)), &
                 optional, intent(inout) :: Kd_lay!< The diapycnal diffusivity in the layers
                                                  !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                 optional, intent(inout) :: Kd_int!< The diapycnal diffusivity at interfaces
                                                  !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  ! Local variables
                                             ! use in the Southern Ocean [nondim].  If this is smaller
                                             ! than Schmittner_coeff, that standard value is used.
                                             ! to model coordinates [R Z3 T-3 ~> W m-2]
                                             ! parameterization [nondim]
                         ! related to the distribution of tidal mixing energy, with unusual array
                         ! extents that are not explained, that is set and used by the CVMix
                         ! tidal mixing schemes, perhaps in [m3 kg-1]?

                                     ! TODO: when coupled, get this from CESM (SHR_CONST_RHOFW)

end subroutine calculate_CVMix_tidal
module subroutine add_int_tide_diffusivity(dz, j, N2_bot, Rho_bot, N2_lay, TKE_to_Kd, max_TKE, &
                                    G, GV, US, CS, Kd_max, Kd_lay, Kd_int, VBF)
  type(ocean_grid_type),             intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),           intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),             intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: dz     !< The vertical distance across layers [Z ~> m]
  integer,                           intent(in)    :: j      !< The j-index to work on
  real, dimension(SZI_(G)),          intent(in)    :: N2_bot !< The near-bottom squared buoyancy frequency
                                                             !! frequency [T-2 ~> s-2].
  real, dimension(SZI_(G)),          intent(in)    :: Rho_bot !< The near-bottom in situ density [R ~> kg m-3]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: N2_lay !< The squared buoyancy frequency of the
                                                             !! layers [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: TKE_to_Kd !< The conversion rate between the TKE
                                                             !! dissipated within a layer and the
                                                             !! diapycnal diffusivity within that layer,
                                                             !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                             !! [T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: max_TKE !< The energy required for a layer
                                                             !! to entrain to its maximum realizable
                                                             !! thickness [H Z2 T-3 ~> m3 s-3 or W m-2]
  type(tidal_mixing_cs),             intent(inout) :: CS     !< The control structure for this module
  real,                              intent(in)    :: Kd_max !< The maximum increment for diapycnal
                                                             !! diffusivity due to TKE-based processes
                                                             !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
                                                             !! Set this to a negative value to have no limit.
  real, dimension(SZI_(G),SZK_(GV)), &
                           optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers
                                                             !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                           optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces
                                                             !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(vbf_CS), pointer                            :: VBF    !< A diagnostics structure for vertical buoyancy fluxes

  ! local

                        ! z*=int(N2/N2_bot) * N2_bot/N2_meanz = int(N2/N2_meanz)
                        ! z0_Polzin_scaled = z0_Polzin * N2_bot/N2_meanz
end subroutine add_int_tide_diffusivity
module subroutine setup_tidal_diagnostics(G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure
  type(tidal_mixing_cs),   intent(inout) :: CS !< The control structure for this module

  ! local

end subroutine setup_tidal_diagnostics
module subroutine post_tidal_diagnostics(G, GV, h ,CS)
  type(ocean_grid_type),    intent(in)   :: G   !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)   :: GV  !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                            intent(in)   :: h   !< Layer thicknesses [H ~> m or kg m-2].
  type(tidal_mixing_cs),    intent(inout) :: CS !< The control structure for this module

end subroutine post_tidal_diagnostics
module subroutine tidal_mixing_h_amp(h_amp, G, j, CS)
  type(ocean_grid_type),    intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G)), intent(out) :: h_amp !< The topographic roughness amplitude [Z ~> m]
  integer,                  intent(in)  :: j     !< j-index of the row to work on
  type(tidal_mixing_cs),    intent(in)  :: CS    !< The control structure for this module


end subroutine tidal_mixing_h_amp
module subroutine read_tidal_energy(G, GV, US, tidal_energy_type, param_file, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  character(len=20),       intent(in) :: tidal_energy_type !< The type of tidal energy inputs to read
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(tidal_mixing_cs),   intent(inout) :: CS   !< The control structure for this module

  ! local variables

end subroutine read_tidal_energy
module subroutine read_tidal_constituents(G, GV, US, tidal_energy_file, param_file, CS)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(unit_scale_type), intent(in) :: US   !< A dimensional unit scaling type
  character(len=200),    intent(in) :: tidal_energy_file !< The file from which to read tidal energy inputs
  type(param_file_type), intent(in)    :: param_file !< Run-time parameter file handle
  type(tidal_mixing_cs), intent(inout) :: CS   !< The control structure for this module

  ! local variables

end subroutine read_tidal_constituents
module subroutine tidal_mixing_end(CS)
  type(tidal_mixing_cs), intent(inout) :: CS !< This module's control structure, which
                                             !! will be deallocated in this routine.

  ! TODO: deallocate all the dynamically allocated members here ...
end subroutine tidal_mixing_end
  end interface

end module MOM_tidal_mixing
