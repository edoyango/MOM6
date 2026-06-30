! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Template for user to code up surface forcing.
module user_surface_forcing

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public USER_wind_forcing, USER_buoyancy_forcing, USER_surface_forcing_init

!>   This control structure should be used to store any run-time variables
!! associated with the user-specified forcing.
!!
!! It can be readily modified for a specific case, and because it is private there
!! will be no changes needed in other code (although they will have to be recompiled).
type, public :: user_surface_forcing_CS ; private
  !   The variables in the canonical example are used for some common
  ! cases, but do not need to be used.

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: rho_restore        !< The density that is used to convert piston velocities into salt
                             !! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]
  real :: gust_const         !< A constant unresolved background gustiness
                             !! that contributes to ustar [R Z2 T-2 ~> Pa].

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type user_surface_forcing_CS


  interface
module subroutine USER_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),            intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),         intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(user_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned
                                                       !! by a previous call to user_surface_forcing_init

  ! Local variables

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
end subroutine USER_wind_forcing
module subroutine USER_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [T ~> s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(user_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned
                                                       !! by a previous call to user_surface_forcing_init

!    This subroutine specifies the current surface fluxes of buoyancy or
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.

!    When temperature is used, there are long list of fluxes that need to be
!  set - essentially the same as for a full coupled model, but most of these
!  can be simply set to zero.  The net fresh water flux should probably be
!  set in fluxes%evap and fluxes%lprec, with any salinity restoring
!  appearing in fluxes%vprec, and the other water flux components
!  (fprec, lrunoff and frunoff) left as arrays full of zeros.
!  Evap is usually negative and precip is usually positive.  All heat fluxes
!  are in W m-2 and positive for heat going into the ocean.  All fresh water
!  fluxes are in [R Z T-1 ~> kg m-2 s-1] and positive for water moving into the ocean.

  ! Local variables
                         ! toward [R ~> kg m-3].
                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].


end subroutine USER_buoyancy_forcing
module subroutine USER_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),               intent(in) :: Time !< The current model time
  type(ocean_grid_type),         intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),         intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,       intent(in) :: diag !< A structure that is used to regulate diagnostic output.
  type(user_surface_forcing_CS), pointer    :: CS   !< A pointer that is set to point to
                                                    !! the control structure for this module

  ! This include declares and sets the variable "version".

end subroutine USER_surface_forcing_init
  end interface

end module user_surface_forcing
