! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Sets forcing for the MESO configuration
module MESO_surface_forcing

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public MESO_buoyancy_forcing, MESO_surface_forcing_init

!> This control structure is used to store parameters associated with the MESO forcing.
type, public :: MESO_surface_forcing_CS ; private

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: rho_restore        !< The density that is used to convert piston velocities into salt
                             !! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]
  real :: gust_const         !< A constant unresolved background gustiness
                             !! that contributes to ustar [R L Z T-2 ~> Pa]
  real, dimension(:,:), pointer :: &
    T_Restore(:,:) => NULL(), & !< The temperature to restore the SST toward [C ~> degC].
    S_Restore(:,:) => NULL(), & !< The salinity to restore the sea surface salnity toward [S ~> ppt]
    PmE(:,:) => NULL(), &       !< The prescribed precip minus evap [Z T-1 ~> m s-1].
    Solar(:,:) => NULL()        !< The shortwave forcing into the ocean [Q R Z T-1 ~> W m-2].
  real, dimension(:,:), pointer :: Heat(:,:) => NULL() !< The prescribed longwave, latent and sensible
                                !! heat flux into the ocean [Q R Z T-1 ~> W m-2].
  character(len=200) :: inputdir !< The directory where NetCDF input files are.
  character(len=200) :: salinityrestore_file !< The file with the target sea surface salinity
  character(len=200) :: SSTrestore_file !< The file with the target sea surface temperature
  character(len=200) :: Solar_file !< The file with the shortwave forcing
  character(len=200) :: heating_file !< The file with the longwave, latent, and sensible heating
  character(len=200) :: PmE_file !< The file with precipitation minus evaporation
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type MESO_surface_forcing_CS

logical :: first_call = .true. !< True until after the first call to the MESO forcing routines


  interface
module subroutine MESO_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [T ~> s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(MESO_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned by
                                                       !! a previous call to MESO_surface_forcing_init

!    When temperature is used, there are long list of fluxes that need to be
!  set - essentially the same as for a full coupled model, but most of these
!  can be simply set to zero.  The net fresh water flux should probably be
!  set in fluxes%evap and fluxes%lprec, with any salinity restoring
!  appearing in fluxes%vprec, and the other water flux components
!  (fprec, lrunoff and frunoff) left as arrays full of zeros.
!  Evap is usually negative and precip is usually positive.  All heat fluxes
!  are in W m-2 and positive for heat going into the ocean.  All fresh water
!  fluxes are in kg m-2 s-1 and positive for water moving into the ocean.

                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].


end subroutine MESO_buoyancy_forcing
module subroutine MESO_surface_forcing_init(Time, G, US, param_file, diag, CS)

  type(time_type),               intent(in)    :: Time !< The current model time
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,       intent(inout) :: diag !< structure used to regulate diagnostic output
  type(MESO_surface_forcing_CS), pointer       :: CS   !< A pointer that is set to point to the
                                                       !! control structure for this module

  ! This include declares and sets the variable "version".

end subroutine MESO_surface_forcing_init
  end interface

end module MESO_surface_forcing
