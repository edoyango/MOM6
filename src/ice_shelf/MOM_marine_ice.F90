! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines incorporating the effects of marine ice (sea-ice and icebergs) into
!! the ocean model dynamics and thermodynamics.
module MOM_marine_ice

use MOM_constants,     only : hlf
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : TO_ALL, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : allocate_forcing_type
use MOM_forcing_type,  only : forcing, mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface

implicit none ; private

#include <MOM_memory.h>

public iceberg_forces, iceberg_fluxes, marine_ice_init

!> Control structure for MOM_marine_ice
type, public :: marine_ice_CS ; private
  real :: kv_iceberg          !< The viscosity of the icebergs [L4 Z-2 T-1 ~> m2 s-1] (for ice rigidity)
  real :: berg_area_threshold !< Fraction of grid cell which iceberg must occupy
                              !! so that fluxes below are set to zero [nondim]. (0.5 is a
                              !! good value to use.) Not applied for negative values.
  real :: latent_heat_fusion  !< Latent heat of fusion [Q ~> J kg-1]
  real :: density_iceberg     !< A typical density of icebergs [R ~> kg m-3] (for ice rigidity)

  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
end type marine_ice_CS


  interface
module subroutine iceberg_forces(G, forces, use_ice_shelf, sfc_state, time_step, CS)
  type(ocean_grid_type), intent(inout) :: G       !< The ocean's grid structure
  type(mech_forcing),    intent(inout) :: forces  !< A structure with the driving mechanical forces
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  logical,               intent(in)    :: use_ice_shelf  !< If true, this configuration uses ice shelves.
  real,                  intent(in)    :: time_step  !< The coupling time step [T ~> s].
  type(marine_ice_CS),   pointer       :: CS      !< Pointer to the control structure for MOM_marine_ice

end subroutine iceberg_forces
module subroutine iceberg_fluxes(G, US, fluxes, use_ice_shelf, sfc_state, time_step, CS)
  type(ocean_grid_type), intent(inout) :: G       !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US      !< A dimensional unit scaling type
  type(forcing),         intent(inout) :: fluxes  !< A structure with pointers to themodynamic,
                                                  !! tracer and mass exchange forcing fields
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  logical,               intent(in)    :: use_ice_shelf  !< If true, this configuration uses ice shelves.
  real,                  intent(in)    :: time_step   !< The coupling time step [T ~> s].
  type(marine_ice_CS),   pointer       :: CS      !< Pointer to the control structure for MOM_marine_ice

end subroutine iceberg_fluxes
module subroutine marine_ice_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(marine_ice_CS),     pointer       :: CS   !< Pointer to the control structure for MOM_marine_ice

  ! This include declares and sets the variable "version".

end subroutine marine_ice_init
  end interface

end module MOM_marine_ice
