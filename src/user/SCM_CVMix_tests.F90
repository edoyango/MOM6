! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initial conditions and forcing for the single column model (SCM) CVMix test set.
module SCM_CVMix_tests

use MOM_domains,       only : pass_var, pass_vector, TO_ALL
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing, mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_verticalgrid,  only : verticalGrid_type
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_unit_scaling,  only : unit_scale_type
use MOM_time_manager,  only : time_type, operator(+), operator(/), time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface

implicit none ; private

#include <MOM_memory.h>

public SCM_CVMix_tests_TS_init
public SCM_CVMix_tests_surface_forcing_init
public SCM_CVMix_tests_wind_forcing
public SCM_CVMix_tests_buoyancy_forcing
public SCM_CVMix_tests_CS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Container for surface forcing parameters
type SCM_CVMix_tests_CS ; private
  logical :: UseWindStress  !< True to use wind stress
  logical :: UseHeatFlux    !< True to use heat flux
  logical :: UseEvaporation !< True to use evaporation
  logical :: UseDiurnalSW   !< True to use diurnal sw radiation
  real :: tau_x !< (Constant) Wind stress, X [R L Z T-2 ~> Pa]
  real :: tau_y !< (Constant) Wind stress, Y [R L Z T-2 ~> Pa]
  real :: surf_HF !< (Constant) Heat flux [C Z T-1 ~> m degC s-1]
  real :: surf_evap !< (Constant) Evaporation rate [Z T-1 ~> m s-1]
  real :: Max_sw !< maximum of diurnal sw radiation [C Z T-1 ~> degC m s-1]
  real :: Rho0 !< reference density [R ~> kg m-3]
  real :: rho_restore !< The density that is used to convert piston velocities
                      !! into salt or heat fluxes [R ~> kg m-3]
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "SCM_CVMix_tests" !< This module's name.


  interface
module subroutine SCM_CVMix_tests_TS_init(T, S, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [Z ~> m]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Input parameter structure
  logical,                                   intent(in)  :: just_read !< If present and true, this call
                                                               !! will only read parameters without changing T & S.
  ! Local variables

end subroutine SCM_CVMix_tests_TS_init
module subroutine SCM_CVMix_tests_surface_forcing_init(Time, G, param_file, CS)
  type(time_type),          intent(in) :: Time       !< Model time
  type(ocean_grid_type),    intent(in) :: G          !< Grid structure
  type(param_file_type),    intent(in) :: param_file !< Input parameter structure
  type(SCM_CVMix_tests_CS), pointer    :: CS         !< Parameter container


  ! This include declares and sets the variable "version".

end subroutine SCM_CVMix_tests_surface_forcing_init
module subroutine SCM_CVMix_tests_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(in)    :: sfc_state  !< Surface state structure
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time in days
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(SCM_CVMix_tests_CS), pointer       :: CS     !< Container for SCM parameters
  ! Local variables
  ! Bounds for loops and memory allocation
end subroutine SCM_CVMix_tests_wind_forcing
module subroutine SCM_CVMix_tests_buoyancy_forcing(sfc_state, fluxes, day, G, US, CS)
  type(surface),            intent(in)    :: sfc_state  !< Surface state structure
  type(forcing),            intent(inout) :: fluxes !< Surface fluxes structure
  type(time_type),          intent(in)    :: day    !< Current model time
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(SCM_CVMix_tests_CS), pointer       :: CS     !< Container for SCM parameters

  ! Local variables

end subroutine SCM_CVMix_tests_buoyancy_forcing
  end interface

end module SCM_CVMix_tests
