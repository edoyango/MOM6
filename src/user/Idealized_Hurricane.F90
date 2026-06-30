! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Forcing for the idealized hurricane and SCM_idealized_hurricane examples.
module Idealized_hurricane

! History
!--------
! November 2014: Origination.
! October 2018: Renamed module from SCM_idealized_hurricane to idealized_hurricane
!               This module is no longer exclusively for use in SCM mode.
!               Legacy code that can be deleted is at the bottom (currently maintained
!               only to preserve exact answers in SCM mode).
!               The T/S initializations have been removed since they are redundant
!               w/ T/S initializations in CVMix_tests (which should be moved
!               into the main state_initialization to their utility
!               for multiple example cases).
! December 2024: Removed the legacy subroutine SCM_idealized_hurricane_wind_forcing

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), time_to_real
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public idealized_hurricane_wind_init !Public interface to initialize the idealized
                                     ! hurricane wind profile.
public idealized_hurricane_wind_forcing !Public interface to update the idealized
                                        ! hurricane wind profile.

!> Container for parameters describing idealized wind structure
type, public :: idealized_hurricane_CS ; private

  ! Parameters used to compute Holland radial wind profile
  real    :: rho_a                !< Mean air density [R ~> kg m-3]
  real    :: pressure_ambient     !< Pressure at surface of ambient air [R L2 T-2 ~> Pa]
  real    :: pressure_central     !< Pressure at surface at hurricane center [R L2 T-2 ~> Pa]
  real    :: rad_max_wind         !< Radius of maximum winds [L ~> m]
  real    :: rad_edge             !< Radius of the edge of the hurricane, normalized by
                                  !! the radius of maximum winds [nondim]
  real    :: rad_ambient          !< Radius at which the winds are at their ambient background values,
                                  !! normalized by the radius of maximum winds [nondim]
  real    :: max_windspeed        !< Maximum wind speeds [L T-1 ~> m s-1]
  real    :: hurr_translation_spd !< Hurricane translation speed [L T-1 ~> m s-1]
  real    :: hurr_translation_dir !< Hurricane translation direction [radians]
  real    :: gustiness            !< Gustiness (used in u*) [R Z2 T-2 ~> Pa]
  real    :: Rho0                 !< A reference ocean density [R ~> kg m-3]
  real    :: Hurr_cen_Y0          !< The initial y position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Hurr_cen_X0          !< The initial x position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Holland_B            !< Parameter 'B' from the Holland formula [nondim]
  logical :: relative_tau         !< A logical to take difference between wind
                                  !! and surface currents to compute the stress
  integer :: answer_date          !< The vintage of the expressions in the idealized hurricane
                                  !! test case.  Values below 20190101 recover the answers
                                  !! from the end of 2018, while higher values use expressions
                                  !! that are rescalable and respect rotational symmetry.
  ! Parameters used in a simple wind-speed dependent expression for C_drag
  real :: Cd_calm       !< The drag coefficient with weak relative winds [nondim]
  real :: calm_speed    !< The relative wind speed below which the drag coefficient takes its
                        !! calm value [L T-1 ~> m s-1]
  real :: Cd_windy      !< The drag coefficient with strong relative winds [nondim]
  real :: windy_speed   !< The relative wind speed below which the drag coefficient takes its
                        !! windy value [L T-1 ~> m s-1]
  real :: dCd_dU10      !< The partial derivative of the drag coefficient times 1000 with the 10 m
                        !! wind speed for intermediate wind speeds [T L-1 ~> s m-1]
  real :: Cd_intercept  !< The zero-wind intercept times 1000 of the linear fit for the drag
                        !! coefficient for the intermediate speeds where there is a linear
                        !! dependence on the 10 m wind speed [nondim]

  ! Parameters used to set the inflow angle as a function of radius and maximum wind speed
  real :: A0_0          !< The zero-radius, zero-speed intercept of the axisymmetric inflow angle [degrees]
  real :: A0_Rnorm      !< The normalized radius dependence of the axisymmetric inflow angle [degrees]
  real :: A0_speed      !< The maximum wind speed dependence of the axisymmetric inflow angle
                        !! [degrees T L-1 ~> degrees s m-1]
  real :: A1_0          !< The zero-radius, zero-speed intercept of the normalized inflow angle
                        !! asymmetry [degrees]
  real :: A1_Rnorm      !< The normalized radius dependence of the normalized inflow angle asymmetry [degrees]
  real :: A1_speed      !< The translation speed dependence of the normalized inflow angle asymmetry
                        !! [degrees T L-1 ~> degrees s m-1]
  real :: P1_0          !< The zero-radius, zero-speed intercept of the angle difference between the
                        !! translation direction and the inflow direction [degrees]
  real :: P1_Rnorm      !< The normalized radius dependence of the angle difference between the
                        !! translation direction and the inflow direction [degrees]
  real :: P1_speed      !< The translation speed dependence of the angle difference between the
                        !! translation direction and the inflow direction [degrees T L-1 ~> degrees s m-1]

  ! Parameters used if in SCM (single column model) mode
  logical :: SCM_mode   !< If true this being used in Single Column Model mode
  logical :: edge_taper_bug !< If true and SCM_mode is true, use a bug that does all of the tapering
                        !! and inflow angle calculations for radii between RAD_EDGE and RAD_AMBIENT
                        !! as though they were at RAD_EDGE.
  real :: f_column      !< Coriolis parameter used in the single column mode idealized
                        !! hurricane wind profile [T-1 ~> s-1]
  logical :: BR_Bench   !< A "benchmark" configuration (which is meant to
                        !! provide identical wind to reproduce a previous
                        !! experiment, where that wind formula contained an error)
  real    :: dy_from_center  !< (Fixed) distance in y from storm center path [L ~> m]

  real :: pi      !< The circumference of a circle divided by its diameter [nondim]
  real :: Deg2Rad !< The conversion factor from degrees to radians [radian degree-1]

end type

character(len=40)  :: mdl = "idealized_hurricane" !< This module's name.


  interface
module subroutine idealized_hurricane_wind_init(Time, G, US, param_file, CS)
  type(time_type),               intent(in) :: Time   !< Model time
  type(ocean_grid_type),         intent(in) :: G      !< Grid structure
  type(unit_scale_type),         intent(in) :: US     !< A dimensional unit scaling type
  type(param_file_type),         intent(in) :: param_file !< Input parameter structure
  type(idealized_hurricane_CS),  pointer    :: CS     !< Parameter container for this module

  ! Local variables
                 ! function of wind speed with the idealized hurricane.  When this is false, the
                 ! linear shape for the mid-range wind speeds is specified separately.

  ! This include declares and sets the variable "version".

end subroutine idealized_hurricane_wind_init
module subroutine idealized_hurricane_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                intent(in)    :: sfc_state  !< Surface state structure
  type(mech_forcing),           intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),              intent(in)    :: day    !< Time in days
  type(ocean_grid_type),        intent(inout) :: G      !< Grid structure
  type(unit_scale_type),        intent(in)    :: US     !< A dimensional unit scaling type
  type(idealized_hurricane_CS), pointer       :: CS     !< Container for idealized hurricane parameters

  ! Local variables

                      !!  benchmark 'f' value [nondim]
                      !!  current relative stress calculation [nondim]

  ! Bounds for loops and memory allocation
end subroutine idealized_hurricane_wind_forcing
module subroutine idealized_hurricane_wind_profile(CS, US, absf, YY, XX, UOCN, VOCN, Tx, Ty)
  ! Author: Brandon Reichl
  ! Date: Nov-20-2014
  !       Aug-14-2018 Generalized for non-SCM configuration

  ! Input parameters
  type(idealized_hurricane_CS), pointer       :: CS   !< Container for idealized hurricane parameters
  type(unit_scale_type),        intent(in)    :: US     !< A dimensional unit scaling type
  real, intent(in)  :: absf !< Input Coriolis magnitude [T-1 ~> s-1]
  real, intent(in)  :: YY   !< Location in m relative to center y [L ~> m]
  real, intent(in)  :: XX   !< Location in m relative to center x [L ~> m]
  real, intent(in)  :: UOCN !< X surface current [L T-1 ~> m s-1]
  real, intent(in)  :: VOCN !< Y surface current [L T-1 ~> m s-1]
  real, intent(out) :: Tx   !< X stress [R L Z T-2 ~> Pa]
  real, intent(out) :: Ty   !< Y stress [R L Z T-2 ~> Pa]

  ! Local variables

  ! Wind profile terms
  ! These variables with weird units are only used with pre-20240501 expressions
                         ! for the Holland profile calculation [m^B R L2 T-2 ~> m^B Pa]
  ! These variables are used with expressions from 20240501 or later
                ! to the power of Holland_B [nondim]

  !Wind angle variables

  ! Implementing Holland (1980) parametric wind profile

end subroutine idealized_hurricane_wind_profile
module function simple_wind_scaled_Cd(u10, du10, CS) result(Cd)
  real,                      intent(in) :: U10  !< The 10 m wind speed [L T-1 ~> m s-1]
  real,                      intent(in) :: du10 !< The magnitude of the difference between the 10 m wind
                                                !! and the ocean flow [L T-1 ~> m s-1]
  type(idealized_hurricane_CS), pointer :: CS   !< Container for SCM parameters
  real :: Cd  ! Air-sea drag coefficient [nondim]

  ! Note that these expressions are discontinuous at dU10 = 11 and 20 m s-1.
end function simple_wind_scaled_Cd
  end interface

end module idealized_hurricane
