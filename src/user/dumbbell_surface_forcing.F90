! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Surface forcing for the dumbbell test case
module dumbbell_surface_forcing

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public dumbbell_dynamic_forcing, dumbbell_buoyancy_forcing, dumbbell_surface_forcing_init

!> Control structure for the dumbbell test case forcing
type, public :: dumbbell_surface_forcing_CS ; private
  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const         !< The restoring rate at the surface [R Z T-1 ~> kg m-2 s-1].
! real :: gust_const         !< A constant unresolved background gustiness
!                            !! that contributes to ustar [R L Z T-2 ~> Pa].
  real :: slp_amplitude      !< The amplitude of pressure loading [R L2 T-2 ~> Pa] applied
                             !! to the reservoirs
  real :: slp_period         !< Period of sinusoidal pressure wave [days]
  real, dimension(:,:), allocatable :: &
    forcing_mask             !< A mask regulating where forcing occurs [nondim]
  real, dimension(:,:), allocatable :: &
    S_restore                !< The surface salinity field toward which to restore [S ~> ppt].
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type dumbbell_surface_forcing_CS


  interface
module subroutine dumbbell_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state  !< A structure containing fields that
                                                         !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing pointers to any
                                                         !! possible forcing fields. Unused fields
                                                         !! have NULL ptrs.
  type(time_type),               intent(in)    :: day    !< Time of the fluxes.
  real,                          intent(in)    :: dt     !< The amount of time over which
                                                         !! the fluxes apply [T ~> s]
  type(ocean_grid_type),         intent(in)    :: G      !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US     !< A dimensional unit scaling type
  type(dumbbell_surface_forcing_CS),  pointer  :: CS     !< A control structure returned by a previous
                                                         !! call to dumbbell_surface_forcing_init
  ! Local variables

end subroutine dumbbell_buoyancy_forcing
module subroutine dumbbell_dynamic_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state  !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing pointers to any
                                                       !! possible forcing fields. Unused fields
                                                       !! have NULL ptrs.
  type(time_type),               intent(in)    :: day  !< Time of the fluxes.
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [T ~> s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(dumbbell_surface_forcing_CS),  pointer  :: CS   !< A control structure returned by a previous
                                                       !! call to dumbbell_surface_forcing_init
  ! Local variables


end subroutine dumbbell_dynamic_forcing
module subroutine dumbbell_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),              intent(in) :: Time !< The current model time.
  type(ocean_grid_type),        intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),      target, intent(in) :: diag !< A structure that is used to
                                                   !! regulate diagnostic output.
  type(dumbbell_surface_forcing_CS), &
                                pointer    :: CS   !< A pointer to the control structure for this module
  ! Local variables
                        ! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]

end subroutine dumbbell_surface_forcing_init
  end interface

end module dumbbell_surface_forcing
