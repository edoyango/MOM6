! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Functions and routines to take area, volume, mass-weighted, layerwise, zonal or meridional means
module MOM_spatial_means

use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=)
use MOM_coms, only : EFP_to_real, real_to_EFP, EFP_sum_across_PEs
use MOM_coms, only : reproducing_sum, reproducing_sum_EFP, EFP_to_real
use MOM_coms, only : query_EFP_overflow_error, reset_EFP_overflow_error
use MOM_coms, only : max_across_PEs, min_across_PEs
use MOM_error_handler, only : MOM_error, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public :: global_i_mean, global_j_mean
public :: global_area_mean, global_area_mean_u, global_area_mean_v, global_layer_mean
public :: global_area_integral
public :: global_volume_mean, global_mass_integral, global_mass_int_EFP
public :: adjust_area_mean_to_zero
public :: array_global_min_max

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].


  interface
module function global_area_mean(var, G, scale, tmp_scale, unscale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                                         !! or arbitrary rescaled units [A ~> a] if unscale
                                                         !! or tmp_scale is present
                                                         !! arbitrary, possibly rescaled units [A ~> a]
  real,                    optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                                         !! that is reversed in the return value [a A-1 ~> 1],
                                                         !! or [b B-1 ~> 1] if unscale is also present.
  real,                    optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums, or
                                                         !! a factor converting between rescaled units if
                                                         !! tmp_scale is also present [B A-1 ~> b a-1].
                                                         !! Here scale and unscale are synonymous, but unscale
                                                         !! is preferred and takes precedence.
  real :: global_area_mean  ! The mean of the variable in arbitrary unscaled units [a] or scaled units [A ~> a]
                            ! or [B ~> b], depending on which optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums.
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
                                          ! scaled cell integral in [A L2 ~> a m2] or [B L2 ~> b m2]
end function global_area_mean
module function global_area_mean_v(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: var  !< The variable to average in arbitrary
                                                         !! units [a], or arbitrary rescaled units
                                                         !! [A ~> a] if tmp_scale is present
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                         !! variable that converts it back to unscaled
                                                         !! (e.g., mks) units to enable the use of the
                                                         !! reproducing sums [a A-1 ~> 1], but is reversed
                                                         !! before output so that the return value has
                                                         !! the same units as var

  real :: global_area_mean_v  ! The mean of the variable in the same arbitrary units as var [A ~> a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same unless tmp_scale and unscale are both present.

end function global_area_mean_v
module function global_area_mean_u(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)  :: var  !< The variable to average in arbitrary
                                                         !! units [a], or arbitrary rescaled units
                                                         !! [A ~> a] if tmp_scale is present
  real,                    optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                         !! variable that converts it back to unscaled
                                                         !! (e.g., mks) units to enable the use of the
                                                         !! reproducing sums [a A-1 ~> 1], but is reversed
                                                         !! before output so that the return value has
                                                         !! the same units as var
  real :: global_area_mean_u  ! The mean of the variable in the same arbitrary units as var [A ~> a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums

end function global_area_mean_u
module function global_area_integral(var, G, scale, area, tmp_scale, unscale)
  type(ocean_grid_type),            intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var   !< The variable to integrate in arbitrary units [a],
                                                         !! or arbitrary rescaled units [A ~> a] if unscale
                                                         !! or tmp_scale is present
  real,                   optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: area !< The alternate area to use, including
                                                         !! any required masking [L2 ~> m2].
  real,                   optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                                         !! that is reversed in the return value [a A-1 ~> 1],
                                                         !! or [b B-1 ~> 1] if unscale is also present.
  real,                   optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                         !! that converts it back to unscaled (e.g., mks)
                                                         !! units to enable the use of the reproducing sums, or
                                                         !! a factor converting between rescaled units if
                                                         !! tmp_scale is also present [B A-1 ~> b a-1].
                                                         !! Here scale and unscale are synonymous, but unscale
                                                         !! is preferred and takes precedence if both are present.
  real :: global_area_integral !< The returned area integral, usually in the units of var times an area,
                               !! [a m2] or [A L2 ~> a m2] or [B L2 ~> b m2], depending on which optional
                               !! arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums.
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
                    ! a scaled cell integral in [B L2 ~> b m2] or other units as indicated below
                    ! or [1] or [B m2 A-1 L-2 ~> b a-1] or [B A-1 ~> b a-1] depending on which
                    ! optional arguments are present.
  !_______________________________________________________________________________________________
  ! Table of units of scalefac and tmpForSumming, depending on the presence of optional arguments |
  !_______________________________________________________________________________________________|
  ! present(tmp_scale) | present(unscale) | scalefac units          | tmpForSumming units         |
  !____________________|__________________|_________________________|_____________________________!
  !      True          |      True        | [B A-1 ~> b a-1]        | [B L2 ~> b m2]              |
  !      True          |      False       | [1]                     | [A L2 ~> a m2]              |
  !      False         |      True        | [a m2 A-1 L-2 ~> b a-1] | [a m2]                      |
  !      False         |      False       | [m2 L-2 ~> 1]           | [a m2]                      |
  !____________________|__________________|_________________________|_____________________________!
end function global_area_integral
module function global_layer_mean(var, h, G, GV, scale, tmp_scale, unscale)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                                                 !! or arbitrary rescaled units [A ~> a] if unscale
                                                                 !! or tmp_scale is present
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real,                            optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                                 !! that converts it back to unscaled (e.g., mks)
                                                                 !! units to enable the use of the reproducing sums
  real,                            optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                                                 !! variable for use in the reproducing sums
                                                                 !! that is reversed in the return value [a A-1 ~> 1],
                                                                 !! or [b B-1 ~> 1] if unscale is also present.
  real,                            optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                                 !! that converts it back to unscaled (e.g., mks)
                                                                 !! units to enable the use of the reproducing sums, or
                                                                 !! a factor converting between rescaled units if
                                                                 !! tmp_scale is also present [B A-1 ~> b a-1].
                                                                 !! Here scale and unscale are synonymous, but unscale
                                                                 !! is preferred and takes precedence.
  real, dimension(SZK_(GV)) :: global_layer_mean  !< The mean of the variable in the arbitrary scaled [A ~> a]
                                                  !! or [B ~> b] or unscaled [a] units of var, depending on which
                                                  !! optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
                                                    ! [L2 a kg m-2 ~> a kg] or a scaled cell integral in
                                                    ! [L2 B m ~> b m3] or [L2 B m ~> b m3] or other units
                                                    ! as indicated the table below.
                                                    ! the model is Boussinesq, used as a weight [L2 m ~> m3]
                                                    ! or [L2 kg m-2 ~> kg]
  !__________________________________________________________________________________________________
  ! Units of weight, scalefac and tmpForSumming, depending on the presence of optional arguments    |
  !_________________________________________________________________________________________________|
  ! Boussinesq | tmp_scale | unscale |  weight units     | scalefac units   | tmpForSumming units   |
  !            | present   | present |                   |                  |                       |
  !____________|___________|_________|___________________|__________________|_______________________!
  !   True     |  True     |  True   | [L2 m ~> m3]      | [B A-1 ~> b a-1] | [B L2 m ~> b m3]      |
  !   True     |  True     |  False  | [L2 m ~> m3]      | [1]              | [A L2 m ~> a m3]      |
  !   True     |  False    |  True   | [L2 m ~> m3]      | [a A-1 ~> 1]     | [L2 a m ~> a m3]      |
  !   True     |  False    |  False  | [L2 m ~> m3]      | [1]              | [L2 a m ~> a m3]      |
  !   False    |  True     |  True   | [L2 kg m-2 ~> kg] | [B A-1 ~> b a-1] | [B L2 kg m-2 ~> b kg] |
  !   False    |  True     |  False  | [L2 kg m-2 ~> kg] | [1]              | [A L2 kg m-2 ~> a kg] |
  !   False    |  False    |  True   | [L2 kg m-2 ~> kg] | [a A-1 ~> 1]     | [L2 a kg m-2 ~> a kg] |
  !   False    |  False    |  False  | [L2 kg m-2 ~> kg] | [1]              | [L2 a kg m-2 ~> a kg] |
  !____________|___________|_________|___________________|__________________|_______________________!
                                        ! half being the tracer integrals in [b m3] or [b kg] and the
                                        ! second half being the summed weights in [m3] or [kg]
                               ! layers [L2 a m ~> a m3] or [L2 a kg m-2 ~> a kg]
                               ! layers [L2 m ~> m3] or [L2 kg m-2 ~> kg]
end function global_layer_mean
module function global_volume_mean(var, h, G, GV, scale, tmp_scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: var  !< The variable to average in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! or tmp_scale is present
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the
                                               !! variable that is reversed in the return value [a A-1 ~> 1],
                                               !! or [b B-1 ~> 1] if unscale is also present.
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums, or
                                               !! a factor converting between rescaled units if
                                               !! tmp_scale is also present [B A-1 ~> b a-1].
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  real :: global_volume_mean  !< The thickness-weighted average of var in the arbitrary scaled [A ~> a] or [B ~> b] or
                              !! unscaled [a] units of var, depending on which optional arguments are provided

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
                                                  ! [B L2 m ~> b m3] or [B L2 kg m-2 ~> b kg] or
                                                  ! [L2 a m ~> a m3] or [L2 a kg m-2 ~> a kg]
                                                  ! [L2 m ~> m3] or [L2 kg m-2 ~> kg]
end function global_volume_mean
module function global_mass_integral(h, G, GV, var, on_PE_only, scale, tmp_scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: var  !< The variable to integrate in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! or tmp_scale is present
  logical,       optional, intent(in)  :: on_PE_only  !< If present and true, the sum is only done
                                               !! on the local PE, and it is _not_ order invariant.
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: tmp_scale !< A temporary rescaling factor for the variable
                                               !! that is reversed in the return value [a A-1 ~> 1],
                                               !! or [b B-1 ~> 1] if unscale is also present.
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums, or
                                               !! a factor converting between rescaled units if
                                               !! tmp_scale is also present [B A-1 ~> b a-1].
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  real :: global_mass_integral  !< The mass-weighted integral of var (or 1) in kg times the arbitrary
                                !! units of var [kg a] or in [R Z L2 A ~> kg a] if tmp_scale is present
                                !! or [R Z L2 B ~> kg b] if both unscale and tmp_scale are present

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  ! [A ~> a] and [B ~> b] are the same units unless tmp_scale and unscale are both present.
                     ! [kg a] or [kg] or if tmp_scale is present in [B R Z L2 ~> kg b] or
                     ! [A R Z L2 !> kg m] or [R Z L2 ~> kg]
                     ! or [kg R-1 Z-1 L-2 ~> 1] or [1] or [B A-1 ~> b a-1] if tmp_scale is present.
                     ! [kg a R-1 Z-1 L-2 A-1 ~> 1] or [kg b R-1 Z-1 L-2 B-1 ~> 1] or [kg R-1 Z-1 L-2 ~> 1]
  !_______________________________________________________________________________________
  ! Units of scalefac and tmpForSumming, depending on the presence of optional arguments |
  !______________________________________________________________________________________|
  ! var     | tmp_scale | unscale |   scalefac units            | tmpForSumming units    |
  ! present | present   | present |                             |                        |
  !_________|___________|_________|_____________________________|________________________!
  !  True   |  True     |  True   | [B A-1 ~> b a-1]            | [B R Z L2 ~> b kg]     |
  !  True   |  True     |  False  | [1]                         | [A R Z L2 ~> a kg]     |
  !  True   |  False    |  True   | [a kg A-1 R-1 Z-1 L-2 ~> 1] | [a kg]                 |
  !  True   |  False    |  False  | [kg R-1 Z-1 L-2 ~> 1]       | [a kg]                 |
  !  False  |  True     |  either | [1]                         | [R Z L2 ~> kg]         |
  !  False  |  False    |  either | [kg R-1 Z-1 L-2 ~> 1]       | [kg]                   |
  !_________|___________|_________|_____________________________|________________________!
end function global_mass_integral
module function global_mass_int_EFP(h, G, GV, var, on_PE_only, scale, unscale)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)  :: var  !< The variable to integrate in arbitrary units [a],
                                               !! or arbitrary rescaled units [A ~> a] if unscale
                                               !! is present
  logical,       optional, intent(in)  :: on_PE_only  !< If present and true, the sum is only done
                                               !! on the local PE, but it is still order invariant.
  real,          optional, intent(in)  :: scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums
  real,          optional, intent(in)  :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                               !! that converts it back to unscaled (e.g., mks)
                                               !! units to enable the use of the reproducing sums.
                                               !! Here scale and unscale are synonymous, but unscale
                                               !! is preferred and takes precedence if both are present.
  type(EFP_type) :: global_mass_int_EFP  !< The mass-weighted integral of var (or 1) in
                                         !! kg times the arbitrary units of var [kg a]

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums

end function global_mass_int_EFP
module subroutine global_i_mean(array, i_mean, G, mask, scale, tmp_scale, unscale)
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array !< The variable to integrate in arbitrary units [a],
                                                           !! or arbitrary rescaled units [A ~> a] if unscale
                                                           !! is present
  real, dimension(SZJ_(G)),         intent(out)   :: i_mean !< Global mean of array along its i-axis [a] or [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: mask  !< An array used for weighting the i-mean [nondim]
  real,                   optional, intent(in)    :: scale !< A rescaling factor for the output variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: tmp_scale !< A rescaling factor for the internal
                                                           !! calculations that is removed from the output [a A-1 ~> 1]
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums.
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! is preferred and takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums

end subroutine global_i_mean
module subroutine global_j_mean(array, j_mean, G, mask, scale, tmp_scale, unscale)
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array !< The variable to integrate in arbitrary units [a],
                                                           !! or arbitrary rescaled units [A ~> a] if unscale
                                                           !! is present
  real, dimension(SZI_(G)),         intent(out)   :: j_mean !<  Global mean of array along its j-axis [a] or [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: mask  !< An array used for weighting the j-mean [nondim]
  real,                   optional, intent(in)    :: scale !< A rescaling factor for the output variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: tmp_scale !< A rescaling factor for the internal
                                                           !! calculations that is removed from the output [a A-1 ~> 1]
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                           !! that converts it back to unscaled (e.g., mks)
                                                           !! units to enable the use of the reproducing sums.
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! is preferred and takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums

end subroutine global_j_mean
module subroutine adjust_area_mean_to_zero(array, G, scaling, unit_scale, unscale)
  type(ocean_grid_type),            intent(in)    :: G       !< Grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: array   !< 2D array to be adjusted in  arbitrary units [a],
                                                             !! or arbitrary rescaled units [A ~> a] if unscale
                                                             !! is present
  real, optional,                   intent(out)   :: scaling !< The scaling factor used [nondim]
  real,                   optional, intent(in)    :: unit_scale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                             !! that converts it back to unscaled (e.g., mks)
                                                             !! units to enable the use of the reproducing sums
  real,                   optional, intent(in)    :: unscale !< A rescaling factor for the variable [a A-1 ~> 1]
                                                             !! that converts it back to unscaled (e.g., mks)
                                                             !! units to enable the use of the reproducing sums.
                                                             !! Here unit_scale and unscale are synonymous, but unscale
                                                             !! is preferred and takes precedence if both are present.
  ! Local variables
  ! In the following comments, [A ~> a] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums

end subroutine adjust_area_mean_to_zero
module subroutine array_global_min_max(tr_array, G, nk, g_min, g_max, &
                                xgmin, ygmin, zgmin, xgmax, ygmax, zgmax, unscale)
  integer,                      intent(in)  :: nk    !< The number of vertical levels
  type(ocean_grid_type),        intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),nk), intent(in)  :: tr_array !< The tracer array to search for
                                                     !! extrema in arbitrary concentration units [CU ~> conc]
  real,                         intent(out) :: g_min !< The global minimum of tr_array, either in
                                                     !! the same units as tr_array [CU ~> conc] or in
                                                     !! unscaled units if unscale is present [conc]
  real,                         intent(out) :: g_max !< The global maximum of tr_array, either in
                                                     !! the same units as tr_array [CU ~> conc] or in
                                                     !! unscaled units if unscale is present [conc]
  real,               optional, intent(out) :: xgmin !< The x-position of the global minimum in the
                                                     !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real,               optional, intent(out) :: ygmin !< The y-position of the global minimum in the
                                                     !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real,               optional, intent(out) :: zgmin !< The z-position of the global minimum [layer]
  real,               optional, intent(out) :: xgmax !< The x-position of the global maximum in the
                                                     !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real,               optional, intent(out) :: ygmax !< The y-position of the global maximum in the
                                                     !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real,               optional, intent(out) :: zgmax !< The z-position of the global maximum [layer]
  real,               optional, intent(in)  :: unscale !< A factor to use to undo any scaling of
                                                     !! the input tracer array [conc CU-1 ~> 1]

  ! Local variables
                             ! maximum values in units that vary between the array elements [various]

end subroutine array_global_min_max
module function ijk_loc(i, j, k, nk, HI)
  integer,              intent(in) :: i   !< Local i-index
  integer,              intent(in) :: j   !< Local j-index
  integer,              intent(in) :: k   !< Local k-index
  integer,              intent(in) :: nk  !< Range of k-index, used to pick out a low-k position.
  type(hor_index_type), intent(in) :: HI  !< Horizontal index ranges
  integer :: ijk_loc  ! An integer encoding the cell position in the global grid.

  ! Local variables

  ! These global i-grid positions run from 1 to HI%niglobal, and analogously for jg.
end function ijk_loc
  end interface

end module MOM_spatial_means
