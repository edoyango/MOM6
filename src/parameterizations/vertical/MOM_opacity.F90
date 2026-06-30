! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines used to calculate the opacity of the ocean.
module MOM_opacity

use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public set_opacity, opacity_init, opacity_end
public extract_optics_slice, extract_optics_fields, optics_nbands
public absorbRemainingSW, sumSWoverBands

!> This type is used to store information about ocean optical properties
type, public :: optics_type
  integer :: nbands     !< The number of penetrating bands of SW radiation

  real, allocatable :: opacity_band(:,:,:,:) !< SW optical depth per unit thickness [Z-1 ~> m-1]
                        !! The number of radiation bands is most rapidly varying (first) index.

  real, allocatable :: sw_pen_band(:,:,:) !< shortwave radiation [Q R Z T-1 ~> W m-2]
                        !! at the surface in each of the nbands bands that penetrates beyond the surface.
                        !! The most rapidly varying dimension is the band.

  real, allocatable :: min_wavelength_band(:)
      !< The minimum wavelength in each band of penetrating shortwave radiation [nm]
  real, allocatable :: max_wavelength_band(:)
      !< The maximum wavelength in each band of penetrating shortwave radiation [nm]

  real :: PenSW_flux_absorb !< A heat flux that is small enough to be completely absorbed in the next
                        !! sufficiently thick layer [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
  real :: PenSW_absorb_Invlen !< The inverse of the thickness that is used to absorb the remaining
  !! shortwave heat flux when it drops below PEN_SW_FLUX_ABSORB [H ~> m or kg m-2].

  !! Lookup tables for Ohlmann solar penetration scheme
  !! These would naturally exist as private module variables but that is prohibited in MOM6
  real :: dlog10chl           !< Chl increment within lookup table  [log10 of Chl in mg m-3]
  real :: chl_min             !< Lower bound of Chl in lookup table [mg m-3]
  real :: log10chl_min        !< Lower bound of Chl in lookup table [log10 of Chl in mg m-3]
  real :: log10chl_max        !< Upper bound of Chl in lookup table [log10 of Chl in mg m-3]
  real, allocatable, dimension(:) :: a1_lut,&       !< Coefficient for band 1 [nondim]
       &                             a2_lut,&       !< Coefficient for band 2 [nondim]
       &                             b1_lut,&       !< Exponential decay scale for band 1 [Z-1 ~> m-1]
       &                             b2_lut         !< Exponential decay scale for band 2 [Z-1 ~> m-1]

  integer :: answer_date  !< The vintage of the order of arithmetic and expressions in the optics
                          !! calculations.  Values below 20190101 recover the answers from the
                          !! end of 2018, while higher values use updated and more robust
                          !! forms of the same expressions.

end type optics_type

!> The control structure with parameters for the MOM_opacity module
type, public :: opacity_CS ; private
  logical :: var_pen_sw      !<   If true, use one of the CHL_A schemes (specified by OPACITY_SCHEME) to
                             !! determine the e-folding depth of incoming shortwave radiation.
  integer :: opacity_scheme  !<   An integer indicating which scheme should be used to translate
                             !! water properties into the opacity (i.e., the e-folding depth) and
                             !! (perhaps) the number of bands of penetrating shortwave radiation to use.
  real :: pen_sw_scale       !<   The vertical absorption e-folding depth of the
                             !! penetrating shortwave radiation [Z ~> m].
  real :: pen_sw_scale_2nd   !<   The vertical absorption e-folding depth of the
                             !! (2nd) penetrating shortwave radiation [Z ~> m].
  real :: SW_1ST_EXP_RATIO   !< Ratio for 1st exp decay in Two Exp decay opacity [nondim]
  real :: pen_sw_frac        !<   The fraction of shortwave radiation that is
                             !! penetrating with a constant e-folding approach [nondim]
  real :: blue_frac          !<   The fraction of the penetrating shortwave
                             !! radiation that is in the blue band [nondim].
  real :: opacity_land_value !< The value to use for opacity over land [Z-1 ~> m-1].
                             !! The default is 10 m-1 - a value for muddy water.
  real, allocatable, dimension(:,:) &
       :: opacity_coef       !< Groups of coefficients, in [Z-1 ~> m-1] or [Z ~> m] depending on the
                             !! scheme, in expressions for opacity, with the second index being the
                             !! wavelength band.  For example, when OPACITY_SCHEME = MANIZZA_05,
                             !! these are coef_1 and coef_2 in the
                             !! expression opacity = coef_1 + coef_2 * chl**pow.
  real, allocatable, dimension(:) &
       :: sw_pen_frac_coef   !< Coefficients in the expression for the penetrating shortwave
                             !! fracetion [nondim]
  real, allocatable, dimension(:) &
       :: chl_power          !< Powers of chlorophyll [nondim] for each band for expressions for
                             !! opacity of the form opacity = coef_1 + coef_2 * chl**pow.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  integer :: chl_dep_bands   !< The number of bands that depend on the Chlorophyll concentrations.
  logical :: warning_issued  !< A flag that is used to avoid repetitive warnings.

  !>@{ Diagnostic IDs
  integer :: id_sw_pen = -1, id_sw_vis_pen = -1
  integer, allocatable :: id_opacity(:)
  !>@}
end type opacity_CS

!>@{ Coded integers to specify the opacity scheme
integer, parameter :: NO_SCHEME = 0, MANIZZA_05 = 1, MOREL_88 = 2, SINGLE_EXP = 3, DOUBLE_EXP = 4,&
     &                OHLMANN_03 = 5
!>@}

character*(10), parameter :: MANIZZA_05_STRING = "MANIZZA_05" !< String to specify the opacity scheme
character*(10), parameter :: MOREL_88_STRING   = "MOREL_88"   !< String to specify the opacity scheme
character*(10), parameter :: OHLMANN_03_STRING = "OHLMANN_03" !< String to specify the opacity scheme
character*(10), parameter :: SINGLE_EXP_STRING = "SINGLE_EXP" !< String to specify the opacity scheme
character*(10), parameter :: DOUBLE_EXP_STRING = "DOUBLE_EXP" !< String to specify the opacity scheme


  interface
module subroutine set_opacity(optics, sw_total, sw_vis_dir, sw_vis_dif, sw_nir_dir, sw_nir_dif, &
                       G, GV, US, CS, chl_2d, chl_3d)
  type(optics_type),       intent(inout) :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  real, dimension(:,:),    pointer       :: sw_total !< Total shortwave flux into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dir !< Visible, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dif !< Visible, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dir !< Near-IR, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dif !< Near-IR, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(opacity_CS)                       :: CS     !< The control structure earlier set up by opacity_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: chl_2d !< Vertically uniform chlorophyll-A concentrations [mg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: chl_3d !< The chlorophyll-A concentrations of each layer [mg m-3]

  ! Local variables
                            ! shortwave radiation [nondim]
                                        ! summed across all bands [Q R Z T-1 ~> W m-2].
                            ! from op to 1/op_diag_len * tanh(op * op_diag_len)
end subroutine set_opacity
module subroutine opacity_from_chl(optics, sw_total, sw_vis_dir, sw_vis_dif, sw_nir_dir, sw_nir_dif, &
                            G, GV, US, CS, chl_2d, chl_3d)
  type(optics_type),       intent(inout) :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  real, dimension(:,:),    pointer       :: sw_total !< Total shortwave flux into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dir !< Visible, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_vis_dif !< Visible, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dir !< Near-IR, direct shortwave into the ocean [Q R Z T-1 ~> W m-2]
  real, dimension(:,:),    pointer       :: sw_nir_dif !< Near-IR, diffuse shortwave into the ocean [Q R Z T-1 ~> W m-2]
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(opacity_CS)                       :: CS     !< The control structure.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: chl_2d !< Vertically uniform chlorophyll-A concentrations [mg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: chl_3d !< A 3-d field of chlorophyll-A concentrations [mg m-3]

                            ! shortwave radiation [nondim]
                            ! near-infrared radiation [nondim]
                            ! shortwave radiation [Q R Z T-1 ~> W m-2].
                            ! radiation [Q R Z T-1 ~> W m-2].
                            ! radiation [Q R Z T-1 ~> W m-2].

end subroutine opacity_from_chl
module function opacity_morel(chl_data, CS)
  real, intent(in)  :: chl_data !< The chlorophyll-A concentration in [mg m-3]
  type(opacity_CS)  :: CS       !< Opacity control structure
  real :: opacity_morel !< The returned opacity [Z-1 ~> m-1]


end function opacity_morel
module function SW_pen_frac_morel(chl_data, CS)
  real, intent(in)  :: chl_data !< The chlorophyll-A concentration [mg m-3]
  type(opacity_CS)  :: CS       !< Opacity control structure
  real :: SW_pen_frac_morel     !< The returned penetrating shortwave fraction [nondim]

  !   The following are coefficients for the optical model taken from Morel and
  ! Antoine (1994). These coefficients represent a non uniform distribution of
  ! chlorophyll-a through the water column.  Other approaches may be more
  ! appropriate when using an interactive ecosystem model that predicts
  ! three-dimensional chl-a values.

end function SW_pen_frac_morel
module subroutine extract_optics_slice(optics, j, G, GV, opacity, opacity_scale, penSW_top, penSW_scale, SpV_avg)
  type(optics_type),       intent(in)  :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer,                 intent(in)  :: j      !< j-index to extract
  type(ocean_grid_type),   intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(max(optics%nbands,1),SZI_(G),SZK_(GV)), &
                 optional, intent(out) :: opacity   !< The opacity in each band, i-point, and layer [Z-1 ~> m-1],
                                                    !! but with units that can be altered by opacity_scale
                                                    !! and the presence of SpV_avg to change this to other
                                                    !! units like [H-1 ~> m-1 or m2 kg-1]
  real,          optional, intent(in)  :: opacity_scale !< A factor by which to rescale the opacity [nondim] or
                                                    !! [Z H-1 ~> 1 or m3 kg-1]
  real, dimension(max(optics%nbands,1),SZI_(G)), &
                 optional, intent(out) :: penSW_top !< The shortwave radiation [Q R Z T-1 ~> W m-2]
                                                    !! at the surface in each of the nbands bands
                                                    !! that penetrates beyond the surface skin layer.
  real,          optional, intent(in)  :: penSW_scale !< A factor by which to rescale the shortwave flux [nondim]
                                                    !! or other units.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: SpV_avg   !< The layer-averaged specific volume [R-1 ~> m3 kg-1]
                                                    !! that is used along with opacity_scale in non-Boussinesq
                                                    !! cases to change the opacity from distance based units to
                                                    !! mass-based units

  ! Local variables
                        ! same units as penSW_scale
end subroutine extract_optics_slice
module subroutine extract_optics_fields(optics, nbands)
  type(optics_type),       intent(in)  :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer, optional,       intent(out) :: nbands !< The number of penetrating bands of SW radiation

end subroutine extract_optics_fields
module function optics_nbands(optics)
  type(optics_type),           pointer :: optics !< An optics structure that has values of opacities
                                                 !! and shortwave fluxes.
  integer :: optics_nbands !< The number of penetrating bands of SW radiation

end function optics_nbands
module subroutine absorbRemainingSW(G, GV, US, h, opacity_band, nsw, optics, j, dt, H_limit_fluxes, &
                             adjustAbsorptionProfile, absorbAllSW, T, Pen_SW_bnd, &
                             eps, ksort, htot, Ttot, TKE, dSV_dT)

  type(ocean_grid_type),             intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),           intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  integer,                           intent(in)    :: nsw  !< Number of bands of penetrating
                                                           !! shortwave radiation.
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(max(1,nsw),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< Opacity in each band of penetrating
                                                           !! shortwave radiation [H-1 ~> m-1 or m2 kg-1].
                                                           !! The indices are band, i, k.
  type(optics_type),                 intent(in)    :: optics !< An optics structure that has values of
                                                           !! opacities and shortwave fluxes.
  integer,                           intent(in)    :: j    !< j-index to work on.
  real,                              intent(in)    :: dt   !< Time step [T ~> s].
  real,                              intent(in)    :: H_limit_fluxes !< If the total ocean depth is
                                                           !! less than this, they are scaled away
                                                           !! to avoid numerical instabilities
                                                           !! [H ~> m or kg m-2]. This would
                                                           !! not be necessary if a finite heat
                                                           !! capacity mud-layer were added.
  logical,                          intent(in)    :: adjustAbsorptionProfile !< If true, apply
                                                           !! heating above the layers in which it
                                                           !! should have occurred to get the
                                                           !! correct mean depth (and potential
                                                           !! energy change) of the shortwave that
                                                           !! should be absorbed by each layer.
  logical,                          intent(in)    :: absorbAllSW !< If true, apply heating above the
                                                           !! layers in which it should have occurred
                                                           !! to get the correct mean depth (and
                                                           !! potential energy change) of the
                                                           !! shortwave that should be absorbed by
                                                           !! each layer.
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: T    !< Layer potential/conservative
                                                           !! temperatures [C ~> degC]
  real, dimension(max(1,nsw),SZI_(G)), intent(inout) :: Pen_SW_bnd !< Penetrating shortwave heating in
                                                           !! each band that hits the bottom and will
                                                           !! will be redistributed through the water
                                                           !! column [C H ~> degC m or degC kg m-2],
                                                           !! size nsw x SZI_(G).
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: eps !< Small thickness that must remain in
                                                           !! each layer, and which will not be
                                                           !! subject to heating [H ~> m or kg m-2]
  integer, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: ksort !< Density-sorted k-indices.
  real, dimension(SZI_(G)), optional, intent(in)    :: htot !< Total mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), optional, intent(inout) :: Ttot !< Depth integrated mixed layer
                                                           !! temperature [C H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(in) :: dSV_dT !< The partial derivative of specific volume
                                                           !! with temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G),SZK_(GV)), optional, intent(inout) :: TKE !< The TKE sink from mixing the heating
                                                           !! throughout a layer [R Z3 T-2 ~> J m-2].

  ! Local variables
                   ! layers above a given layer [C ~> degC].  This is only nonzero if
                   ! adjustAbsorptionProfile is true, in which case the net
                   ! change in the temperature of a layer is the sum of the
                   ! direct heating of that layer plus T_chg_above from all of
                   ! the layers below, plus any contribution from absorbing
                   ! radiation that hits the bottom.
end subroutine absorbRemainingSW
module subroutine sumSWoverBands(G, GV, US, h, dz, nsw, optics, j, dt, &
                          H_limit_fluxes, absorbAllSW, iPen_SW_bnd, netPen)
  type(ocean_grid_type),    intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: h   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: dz  !< Layer vertical extent [Z ~> m].
  integer,                  intent(in)    :: nsw !< The number of bands of penetrating shortwave
                                                 !! radiation, perhaps from optics_nbands(optics),
  type(optics_type),        intent(in)    :: optics !< An optics structure that has values
                                                   !! set based on the opacities.
  integer,                  intent(in)    :: j   !< j-index to work on.
  real,                     intent(in)    :: dt  !< Time step [T ~> s].
  real,                     intent(in)    :: H_limit_fluxes !< the total depth at which the
                                                 !! surface fluxes start to be limited to avoid
                                                 !! excessive heating of a thin ocean [H ~> m or kg m-2]
  logical,                  intent(in)    :: absorbAllSW !< If true, ensure that all shortwave
                                                 !! radiation is absorbed in the ocean water column.
  real, dimension(max(nsw,1),SZI_(G)), intent(in) :: iPen_SW_bnd !< The incident penetrating shortwave
                                                 !! in each band at the sea surface; size nsw x SZI_(G)
                                                 !! [C H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                             intent(inout) :: netPen !< Net penetrating shortwave heat flux at each
                                                 !! interface, summed across all bands
                                                 !! [C H ~> degC m or degC kg m-2].
  ! Local variables
                              ! remaining shortwave radiation [H ~> m or kg m-2].
                              ! penetrating shortwave heating that hits the bottom
                              ! and will be redistributed through the water column
                              ! [C H ~> degC m or degC kg m-2]

                          ! in each band, initially iPen_SW_bnd [C H ~> degC m or degC kg m-2]
                          ! absorbed in a layer [nondim]
                          ! not absorbed because the layers are too thin [nondim].
                          ! surface fluxes start to be limited [H-1 ~> m-1 or m2 kg-1]
                          ! absorbed in the next layer for computational efficiency, instead of
                          ! continuing to penetrate [C H ~> degC m or degC kg m-2].
                          ! was not entirely absorbed.

end subroutine sumSWoverBands
module subroutine opacity_init(Time, G, GV, US, param_file, diag, CS, optics)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< model vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(opacity_CS) :: CS                         !< Opacity control structure
  type(optics_type) :: optics                    !< An optics structure that has parameters
                                                 !! set and arrays allocated here.
  ! Local variables
  ! This include declares and sets the variable "version".
                                ! near-infrared radiation with parameterizations following the
                                ! functional form from Manizza et al., GRL 2005, namely in the form
                                ! opacity = coef_1 + coef_2 * chl**pow for each band.
                                ! radiation bands, in expressions for opacity of the form
                                ! opacity = coef_1 + coef_2 * chl**pow.
                                ! radiation in the form proposed by Morel and Antoine (1994), namely
                                ! opacity = 1 / (sum(n=1:6, Coef(n) * log10(Chl)**(n-1)))
                                ! fifth order polynomial fit as a funciton of log10(Chlorophyll).
                                ! flux when that flux drops below PEN_SW_FLUX_ABSORB [H ~> m or kg m-2]
                                ! radiation bands [nm]
end subroutine opacity_init
module subroutine init_ohlmann_table(optics)

  implicit none

  type(optics_type), intent(inout) :: optics

  ! Local variables

  !! These are the data from Ohlmann (2003) Table 1a with additional
  !! values provided by C. Ohlmann and implemented in CESM-POP by B. Briegleb





  !! Make the table big enough so step size is smaller
  !! in log-space that any increment in Table 1a

end subroutine init_ohlmann_table
module function lookup_ohlmann_swpen(chl,optics) result(A)

  implicit none

  real, intent(in) :: chl
  type(optics_type), intent(in) :: optics
  real, dimension(2) :: A

  ! Local variables


  ! Make sure we are in the table
end function lookup_ohlmann_swpen
module function lookup_ohlmann_opacity(chl,optics) result(B)

  implicit none
  real, intent(in) :: chl
  type(optics_type), intent(in) :: optics
  real, dimension(2) :: B

  ! Local variables

  ! Make sure we are in the table
end function lookup_ohlmann_opacity
module subroutine opacity_end(CS, optics)
  type(opacity_CS)  :: CS     !< Opacity control structure
  type(optics_type) :: optics !< An optics type structure that should be deallocated.

end subroutine opacity_end
  end interface

end module MOM_opacity
