! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Surface forcing for the boundary-forced-basin (BFB) configuration
module BFB_surface_forcing

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public BFB_buoyancy_forcing, BFB_surface_forcing_init

!> Control structure for BFB_surface_forcing
type, public :: BFB_surface_forcing_CS ; private

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: rho_restore        !< The density that is used to convert piston velocities into salt
                             !! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]
  real :: SST_s              !< SST at the southern edge of the linear forcing ramp [C ~> degC]
  real :: SST_n              !< SST at the northern edge of the linear forcing ramp [C ~> degC]
  real :: S_ref              !< Reference salinity used throughout the domain [S ~> ppt]
  real :: lfrslat            !< Southern latitude where the linear forcing ramp begins [degrees_N] or [km]
  real :: lfrnlat            !< Northern latitude where the linear forcing ramp ends [degrees_N] or [km]
  real :: Rho_T0_S0          !< The density at T=0, S=0 [R ~> kg m-3]
  real :: dRho_dT            !< The partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real :: dRho_dS            !< The partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1]
                             !!   Note that temperature and salinity are being used as dummy variables here.
                             !! All temperatures are converted into density.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
end type BFB_surface_forcing_CS


  interface
module subroutine BFB_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                intent(inout) :: sfc_state  !< A structure containing fields that
                                                      !! describe the surface state of the ocean.
  type(forcing),                intent(inout) :: fluxes !< A structure containing pointers to any
                                                      !! possible forcing fields. Unused fields
                                                      !! have NULL ptrs.
  type(time_type),              intent(in)    :: day  !< Time of the fluxes.
  real,                         intent(in)    :: dt   !< The amount of time over which
                                                      !! the fluxes apply [T ~> s]
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(BFB_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure
                                                      !! returned by a previous call to
                                                      !! BFB_surface_forcing_init.
  ! Local variables
                         ! toward [R ~> kg m-3].
                           ! factors [Q R C-1 ~> J m-3 degC-1]
                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].

end subroutine BFB_buoyancy_forcing
module subroutine BFB_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),              intent(in) :: Time !< The current model time.
  type(ocean_grid_type),        intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,      intent(in) :: diag !< A structure that is used to
                                                   !! regulate diagnostic output.
  type(BFB_surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure for this module

  ! This include declares and sets the variable "version".

end subroutine BFB_surface_forcing_init
  end interface

end module BFB_surface_forcing
