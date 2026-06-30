! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides a transparent unit rescaling type to facilitate dimensional consistency testing
module MOM_unit_scaling

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type

implicit none ; private

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, T, R and Q, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the rescaled
! combination is a nondimensional variable, the notation would be "a slope [Z L-1 ~> nondim]",
! but if (as the case for the variables here), the rescaled combination is exactly 1, the right
! notation would be something like "a dimensional scaling factor [Z m-1 ~> 1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

public unit_scaling_init, unit_no_scaling_init, unit_scaling_end, fix_restart_unit_scaling

!> Describes various unit conversion factors
type, public :: unit_scale_type
  real :: m_to_Z     !< A constant that translates distances in meters to the units of depth              [Z m-1 ~> 1]
  real :: Z_to_m     !< A constant that translates distances in the units of depth to meters              [m Z-1 ~> 1]
  real :: m_to_L     !< A constant that translates lengths in meters to the units of horizontal lengths   [L m-1 ~> 1]
  real :: L_to_m     !< A constant that translates lengths in the units of horizontal lengths to meters   [m L-1 ~> 1]
  real :: s_to_T     !< A constant that translates time intervals in seconds to the units of time         [T s-1 ~> 1]
  real :: T_to_s     !< A constant that translates the units of time to seconds                           [s T-1 ~> 1]
  real :: R_to_kg_m3 !< A constant that translates the units of density to kilograms per meter cubed [kg m-3 R-1 ~> 1]
  real :: kg_m3_to_R !< A constant that translates kilograms per meter cubed to the units of density  [R m3 kg-1 ~> 1]
  real :: Q_to_J_kg  !< A constant that translates the units of enthalpy to Joules per kilogram      [J kg-1 Q-1 ~> 1]
  real :: J_kg_to_Q  !< A constant that translates Joules per kilogram to the units of enthalpy        [Q kg J-1 ~> 1]
  real :: C_to_degC  !< A constant that translates the units of temperature to degrees Celsius         [degC C-1 ~> 1]
  real :: degC_to_C  !< A constant that translates degrees Celsius to the units of temperature         [C degC-1 ~> 1]
  real :: S_to_ppt   !< A constant that translates the units of salinity to parts per thousand          [ppt S-1 ~> 1]
  real :: ppt_to_S   !< A constant that translates parts per thousand to the units of salinity          [S ppt-1 ~> 1]

  ! These are useful combinations of the fundamental scale conversion factors above.
  real :: Z_to_L          !< Convert vertical distances to lateral lengths                                [L Z-1 ~> 1]
  real :: L_to_Z          !< Convert lateral lengths to vertical distances                                [Z L-1 ~> 1]
  real :: L_T_to_m_s      !< Convert lateral velocities from L T-1 to m s-1                         [T m L-1 s-1 ~> 1]
  real :: m_s_to_L_T      !< Convert lateral velocities from m s-1 to L T-1                         [L s T-1 m-1 ~> 1]
  real :: L_T2_to_m_s2    !< Convert lateral accelerations from L T-2 to m s-2                     [L s2 T-2 m-1 ~> 1]
  real :: Z2_T_to_m2_s    !< Convert vertical diffusivities from Z2 T-1 to m2 s-1                  [T m2 Z-2 s-1 ~> 1]
  real :: m2_s_to_Z2_T    !< Convert vertical diffusivities from m2 s-1 to Z2 T-1                  [Z2 s T-1 m-2 ~> 1]
  real :: W_m2_to_QRZ_T   !< Convert heat fluxes from W m-2 to Q R Z T-1                       [Q R Z m2 T-1 W-1 ~> 1]
  real :: QRZ_T_to_W_m2   !< Convert heat fluxes from Q R Z T-1 to W m-2                    [W T Q-1 R-1 Z-1 m-2 ~> 1]
  ! Not used enough:  real :: kg_m2_to_RZ   !< Convert mass loads from kg m-2 to R Z                [R Z m2 kg-1 ~> 1]
  real :: RZ_to_kg_m2     !< Convert mass loads from R Z to kg m-2                               [kg R-1 Z-1 m-2 ~> 1]
  real :: RZL2_to_kg      !< Convert masses from R Z L2 to kg                                    [kg R-1 Z-1 L-2 ~> 1]
  real :: kg_m2s_to_RZ_T  !< Convert mass fluxes from kg m-2 s-1 to R Z T-1                   [R Z m2 s T-1 kg-1 ~> 1]
  real :: RZ_T_to_kg_m2s  !< Convert mass fluxes from R Z T-1 to kg m-2 s-1                [T kg R-1 Z-1 m-2 s-1 ~> 1]
  real :: RZ3_T3_to_W_m2  !< Convert turbulent kinetic energy fluxes from R Z3 T-3 to W m-2    [W T3 R-1 Z-3 m-2 ~> 1]
  real :: W_m2_to_RZ3_T3  !< Convert turbulent kinetic energy fluxes from W m-2 to R Z3 T-3     [R Z3 m2 T-3 W-1 ~> 1]
  real :: RL2_T2_to_Pa    !< Convert pressures from R L2 T-2 to Pa                                [Pa T2 R-1 L-2 ~> 1]
  real :: RLZ_T2_to_Pa    !< Convert wind stresses from R L Z T-2 to Pa                       [Pa T2 R-1 L-1 Z-1 ~> 1]
  real :: Pa_to_RL2_T2    !< Convert pressures from Pa to R L2 T-2                                [R L2 T-2 Pa-1 ~> 1]
  real :: Pa_to_RLZ_T2    !< Convert wind stresses from Pa to R L Z T-2                          [R L Z T-2 Pa-1 ~> 1]

  ! These are no longer used for changing scaling across restarts.
  real :: m_to_Z_restart = 1.0 !< A copy of the m_to_Z that is used in restart files.
  real :: m_to_L_restart = 1.0 !< A copy of the m_to_L that is used in restart files.
  real :: s_to_T_restart = 1.0 !< A copy of the s_to_T that is used in restart files.
  real :: kg_m3_to_R_restart = 1.0 !< A copy of the kg_m3_to_R that is used in restart files.
  real :: J_kg_to_Q_restart = 1.0 !< A copy of the J_kg_to_Q that is used in restart files.
end type unit_scale_type


  interface
module subroutine unit_scaling_init( param_file, US )
  type(param_file_type), intent(in) :: param_file !< Parameter file handle/type
  type(unit_scale_type), pointer    :: US         !< A dimensional unit scaling type

  ! This routine initializes a unit_scale_type structure (US).

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine unit_scaling_init
module subroutine unit_no_scaling_init(US)
  type(unit_scale_type), pointer    :: US         !< A dimensional unit scaling type

end subroutine unit_no_scaling_init
module subroutine set_unit_scaling_combos(US)
  type(unit_scale_type), intent(inout) :: US !< A dimensional unit scaling type

  ! Convert vertical to horizontal length scales and the reverse:
end subroutine set_unit_scaling_combos
module subroutine fix_restart_unit_scaling(US, unscaled)
  type(unit_scale_type), intent(inout) :: US !< A dimensional unit scaling type
  logical,     optional, intent(in)    :: unscaled !< If true, set the restart factors as though the
                                             !! model would be unscaled, which is appropriate if the
                                             !! scaling is undone when writing a restart file.

end subroutine fix_restart_unit_scaling
module subroutine unit_scaling_end( US )
  type(unit_scale_type), pointer :: US !< A dimensional unit scaling type

end subroutine unit_scaling_end
  end interface

end module MOM_unit_scaling
