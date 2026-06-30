! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module provides a common datatype to provide forcing for MARBL tracers
!! regardless of driver
module MARBL_forcing_mod

!! This module exists to house code used by multiple drivers in config_src/
!! for passing forcing fields to MARBL
!! (This comment can go in the wiki on the NCAR fork?)

use MOM_diag_mediator,        only : safe_alloc_ptr, diag_ctrl, register_diag_field, post_data
use MOM_time_manager,         only : time_type
use MOM_error_handler,        only : MOM_error, WARNING, FATAL
use MOM_file_parser,          only : get_param, log_param, param_file_type
use MOM_grid,                 only : ocean_grid_type
use MOM_unit_scaling,         only : unit_scale_type
use MOM_interpolate,          only : external_field, init_external_field, time_interp_external
use MOM_io,                   only : slasher
use marbl_constants_mod,      only : molw_Fe
use MOM_forcing_type,         only : forcing

implicit none ; private

#include <MOM_memory.h>

public :: MARBL_forcing_init
public :: convert_driver_fields_to_forcings

!> Data type used to store diagnostic index returned from register_diag_field()
!! For the forcing fields that can be written via post_data()
type, private :: marbl_forcing_diag_ids
  integer :: atm_fine_dust   !< Atmospheric fine dust component of dust_flux
  integer :: atm_coarse_dust !< Atmospheric coarse dust component of dust_flux
  integer :: atm_bc          !< Atmospheric black carbon component of iron_flux
  integer :: ice_dust        !< Sea-ice dust component of dust_flux
  integer :: ice_bc          !< Sea-ice black carbon component of iron_flux
end type marbl_forcing_diag_ids

!> Control structure for this module
type, public :: marbl_forcing_CS ; private
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.

  real    :: dust_ratio_thres               !< coarse/fine dust ratio threshold [1]
  real    :: dust_ratio_to_fe_bioavail_frac !< ratio of dust to iron bioavailability fraction [1]
  real    :: fe_bioavail_frac_offset        !< offset for iron bioavailability fraction [1]
  real    :: atm_fe_to_bc_ratio             !< atmospheric iron to black carbon ratio [1]
  real    :: atm_bc_fe_bioavail_frac        !< atmospheric black carbon to iron bioavailablity fraction ratio [1]
  real    :: seaice_fe_to_bc_ratio          !< sea-ice iron to black carbon ratio [1]
  real    :: seaice_bc_fe_bioavail_frac     !< sea-ice black carbon to iron bioavailablity fraction ratio [1]
  real    :: iron_frac_in_atm_fine_dust     !< Fraction of fine dust from the atmosphere that is iron [1]
  real    :: iron_frac_in_atm_coarse_dust   !< Fraction of coarse dust from the atmosphere that is iron [1]
  real    :: iron_frac_in_seaice_dust       !< Fraction of dust from the sea ice that is iron [1]
  real    :: atm_co2_const                  !< atmospheric CO2 (if specifying a constant value) [ppm]
  real    :: atm_alt_co2_const              !< alternate atmospheric CO2 for _ALT_CO2 tracers
                                            !! (if specifying a constant value) [ppm]

  type(marbl_forcing_diag_ids) :: diag_ids  !< used for registering and posting some MARBL forcing fields as diagnostics

  logical :: use_marbl_tracers    !< most functions can return immediately
                                  !! MARBL tracers are turned off
  integer :: atm_co2_iopt         !< Integer version of atm_co2_opt, which determines source of atm_co2
  integer :: atm_alt_co2_iopt     !< Integer version of atm_alt_co2_opt, which determines source of atm_alt_co2

end type marbl_forcing_CS

! Module parameters
integer, parameter :: atm_co2_constant_iopt = 0     !< module parameter denoting atm_co2_opt = 'constant'
integer, parameter :: atm_co2_prognostic_iopt = 1   !< module parameter denoting atm_co2_opt = 'diagnostic'
integer, parameter :: atm_co2_diagnostic_iopt = 2   !< module parameter denoting atm_co2_opt = 'prognostic'


  interface
  module subroutine MARBL_forcing_init(G, US, param_file, diag, day, inputdir, use_marbl, CS)
    type(ocean_grid_type),           intent(in)    :: G           !< The ocean's grid structure
    type(unit_scale_type),           intent(in)    :: US          !< A dimensional unit scaling type
    type(param_file_type),           intent(in)    :: param_file  !< A structure to parse for run-time parameters
    type(diag_ctrl), target,         intent(in)    :: diag        !< Structure used to regulate diagnostic output.
    type(time_type), target,         intent(in)    :: day         !< Time of the start of the run.
    character(len=*),                intent(in)    :: inputdir    !< Directory containing input files
    logical,                         intent(in)    :: use_marbl   !< Is MARBL tracer package active?
    type(marbl_forcing_CS), pointer, intent(inout) :: CS          !< A pointer that is set to point to control
                                                                  !! structure for MARBL forcing


  end subroutine MARBL_forcing_init
  module subroutine convert_driver_fields_to_forcings(atm_fine_dust_flux, atm_coarse_dust_flux, &
                                               seaice_dust_flux, atm_bc_flux, seaice_bc_flux, &
                                               nhx_dep, noy_dep, atm_co2_prog, atm_co2_diag, &
                                               afracr, swnet_afracr, ifrac_n, &
                                               swpen_ifrac_n, Time, G, US, i0, j0, fluxes, CS)

    real, dimension(:,:),   pointer, intent(in)    :: atm_fine_dust_flux   !< atmosphere fine dust flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_coarse_dust_flux !< atmosphere coarse dust flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: seaice_dust_flux     !< sea ice dust flux from IOB [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_bc_flux          !< atmosphere black carbon flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: seaice_bc_flux       !< sea ice black carbon flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: afracr               !< open ocean fraction [1]
    real, dimension(:,:),   pointer, intent(in)    :: nhx_dep              !< NHx flux from atmosphere [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: noy_dep              !< NOy flux from atmosphere [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_co2_prog         !< Prognostic atmospheric CO2 concentration
                                                                           !! [ppm]
    real, dimension(:,:),   pointer, intent(in)    :: atm_co2_diag         !< Diagnostic atmospheric CO2 concentration
                                                                           !! [ppm]
    real, dimension(:,:),   pointer, intent(in)    :: swnet_afracr         !< shortwave flux * open ocean fraction
                                                                           !! [W m-2]
    real, dimension(:,:,:), pointer, intent(in)    :: ifrac_n              !< per-category ice fraction [1]
    real, dimension(:,:,:), pointer, intent(in)    :: swpen_ifrac_n        !< per-category shortwave flux * ice fraction
                                                                           !! [W m-2]
    type(time_type),                 intent(in)    :: Time                 !< The time of the fluxes, used for
                                                                           !! interpolating the salinity to the
                                                                           !! right time, when it is being
                                                                           !! restored.
    type(ocean_grid_type),           intent(in)    :: G                    !< The ocean's grid structure
    type(unit_scale_type),           intent(in)    :: US                   !< A dimensional unit scaling type
    integer,                         intent(in)    :: i0                   !< i index offset
    integer,                         intent(in)    :: j0                   !< j index offset
    type(forcing),                   intent(inout) :: fluxes               !< MARBL-specific forcing fields
    type(marbl_forcing_CS), pointer, intent(inout) :: CS                   !< A pointer that is set to point to
                                                                           !! control structure for MARBL forcing

    ! Note: following two conversion factors are used to both convert from km m-2 s-1 -> mmol m-2 s-1
    !!      AND cast in MOM6's unique dimensional consistency scaling system [conc Z T-1]
                                     !! [s m2 kg-1 conc Z T-1 ~> mmol kg-1]
                                     !! [s m2 kg-1 conc Z T-1 ~> mmol kg-1]

  end subroutine convert_driver_fields_to_forcings
  end interface

end module MARBL_forcing_mod
