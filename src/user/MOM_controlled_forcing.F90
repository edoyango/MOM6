! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Use control-theory to adjust the surface heat flux and precipitation.
!!
!! Adjustments are based on the time-mean or periodically (seasonally) varying
!! anomalies from the observed state.
!!
!! The techniques behind this are described in Hallberg and Adcroft (2018, in prep.).
module MOM_controlled_forcing

use MOM_diag_mediator, only : post_data, query_averaging_enabled, enable_averages, disable_averaging
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains,       only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,   only : read_param, get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing
use MOM_grid,          only : ocean_grid_type
use MOM_restart,       only : register_restart_field, MOM_restart_CS
use MOM_time_manager,  only : time_type, operator(+), operator(/), operator(-)
use MOM_time_manager,  only : get_date, set_date
use MOM_time_manager,  only : time_type_to_real, real_to_time
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface

implicit none ; private

#include <MOM_memory.h>

public apply_ctrl_forcing, register_ctrl_forcing_restarts
public controlled_forcing_init, controlled_forcing_end

!> Control structure for MOM_controlled_forcing
type, public :: ctrl_forcing_CS ; private
  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: do_integrated  !< If true, use time-integrated anomalies to control the surface state.
  integer :: num_cycle      !< The number of elements in the forcing cycle.
  real    :: heat_int_rate  !< The rate at which heating anomalies accumulate [T-1 ~> s-1]
  real    :: prec_int_rate  !< The rate at which precipitation anomalies accumulate [T-1 ~> s-1]
  real    :: heat_cyc_rate  !< The rate at which cyclical heating anomalies accumulate [T-1 ~> s-1]
  real    :: prec_cyc_rate  !< The rate at which cyclical precipitation anomalies
                            !! accumulate [T-1 ~> s-1]
  real    :: Len2           !< The square of the length scale over which the anomalies
                            !! are smoothed via a Laplacian filter [L2 ~> m2]
  real    :: lam_heat       !< A constant of proportionality between SST anomalies
                            !! and heat fluxes [Q R Z T-1 C-1 ~> W m-2 degC-1]
  real    :: lam_prec       !< A constant of proportionality between SSS anomalies
                            !! (normalised by mean SSS) and precipitation [R Z T-1 ~> kg m-2 s-1]
  real    :: lam_cyc_heat   !< A constant of proportionality between cyclical SST
                            !! anomalies and corrective heat fluxes [Q R Z T-1 C-1 ~> W m-2 degC-1]
  real    :: lam_cyc_prec   !< A constant of proportionality between cyclical SSS
                            !! anomalies (normalised by mean SSS) and corrective
                            !! precipitation [R Z T-1 ~> kg m-2 s-1]

  real, pointer, dimension(:,:) :: &
    heat_0 => NULL(), &     !< The non-periodic integrative corrective heat flux that has been
                            !! evolved to control mean SST anomalies [Q R Z T-1 ~> W m-2]
    precip_0 => NULL()      !< The non-periodic integrative corrective precipitation that has been
                            !! evolved to control mean SSS anomalies [R Z T-1 ~> kg m-2 s-1]

  ! The final dimension of each of the six variables that follow is for the periodic bins.
  real, pointer, dimension(:,:,:) :: &
    heat_cyc => NULL(), &   !< The periodic integrative corrective heat flux that has been evolved
                            !! to control periodic (seasonal) SST anomalies [Q R Z T-1 ~> W m-2].
                            !! The third dimension is the periodic bins.
    precip_cyc => NULL()    !< The non-periodic integrative corrective precipitation that has been
                            !! evolved to control periodic (seasonal) SSS anomalies [R Z T-1 ~> kg m-2 s-1].
                            !! The third dimension is the periodic bins.
  real, pointer, dimension(:) :: &
    avg_time => NULL()      !< The accumulated averaging time in each part of the cycle [T ~> s] or
                            !! a negative value to indicate that the variables like avg_SST_anom are
                            !! the actual averages, and not time integrals.
                            !! The dimension is the periodic bins.
  real, pointer, dimension(:,:,:) :: &
    avg_SST_anom => NULL(), & !< The time-averaged periodic sea surface temperature anomalies [C ~> degC],
                              !! or (at some points in the code), the time-integrated periodic
                              !! temperature anomalies [T C ~> s degC].
                              !! The third dimension is the periodic bins.
    avg_SSS_anom => NULL(), & !< The time-averaged periodic sea surface salinity anomalies [S ~> ppt],
                              !! or (at some points in the code), the time-integrated periodic
                              !! salinity anomalies [T S ~> s ppt].
                              !! The third dimension is the periodic bins.
    avg_SSS => NULL()         !< The time-averaged periodic sea surface salinities [S ~> ppt], or (at
                              !! some points in the code), the time-integrated periodic
                              !! salinities [T S ~> s ppt].
                              !! The third dimension is the periodic bins.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                            !! regulate the timing of diagnostic output.
  integer :: id_heat_0 = -1 !< Diagnostic handle for the steady heat flux
  integer :: id_prec_0 = -1 !< Diagnostic handle for the steady precipitation
end type ctrl_forcing_CS


  interface
module subroutine apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, virt_heat, virt_precip, &
                              day_start, dt, G, US, CS)
  type(ocean_grid_type), intent(inout) :: G         !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SST_anom  !< The sea surface temperature anomalies [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SSS_anom  !< The sea surface salinity anomlies [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SSS_mean  !< The mean sea surface salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: virt_heat !< Virtual (corrective) heat
                                                    !! fluxes that are augmented in this
                                                    !! subroutine [Q R Z T-1 ~> W m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: virt_precip !< Virtual (corrective)
                                                    !! precipitation fluxes that are augmented
                                                    !! in this subroutine [R Z T-1 ~> kg m-2 s-1]
  type(time_type),       intent(in)    :: day_start !< Start time of the fluxes.
  real,                  intent(in)    :: dt        !< Length of time over which these fluxes
                                                    !! will be applied [T ~> s]
  type(unit_scale_type), intent(in)    :: US        !< A dimensional unit scaling type
  type(ctrl_forcing_CS), pointer       :: CS        !< A pointer to the control structure returned
                                                    !! by a previous call to ctrl_forcing_init.

  ! Local variables

end subroutine apply_ctrl_forcing
module function periodic_int(rval, num_period) result (m)
  real,    intent(in) :: rval       !< Input for mapping [nondim]
  integer, intent(in) :: num_period !< Maximum output.
  integer             :: m          !< Return value.

end function periodic_int
module function periodic_real(rval, num_period) result(val_out)
  real,    intent(in) :: rval       !< Input to be shifted into valid range [nondim]
  integer, intent(in) :: num_period !< Maximum valid value.
  real                :: val_out    !< Return value [nondim]

end function periodic_real
module subroutine register_ctrl_forcing_restarts(G, US, param_file, CS, restart_CS)
  type(ocean_grid_type), intent(in) :: G          !< The ocean's grid structure.
  type(unit_scale_type), intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.
  type(ctrl_forcing_CS), pointer :: CS            !< A pointer that is set to point to the
                                                  !! control structure for this module.
  type(MOM_restart_CS), intent(inout) :: restart_CS !< MOM restart control struct

end subroutine register_ctrl_forcing_restarts
module subroutine controlled_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),           intent(in) :: Time       !< The current model time.
  type(ocean_grid_type),     intent(in) :: G          !< The ocean's grid structure.
  type(unit_scale_type),     intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),     intent(in) :: param_file !< A structure indicating the
                                                      !! open file to parse for model
                                                      !! parameter values.
  type(diag_ctrl), target,   intent(in) :: diag       !< A structure that is used to regulate
                                                      !! diagnostic output.
  type(ctrl_forcing_CS),     pointer    :: CS         !< A pointer that is set to point to the
                                                      !! control structure for this module.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine controlled_forcing_init
module subroutine controlled_forcing_end(CS)
  type(ctrl_forcing_CS),    pointer :: CS !< A pointer to the control structure
                                          !! returned by a previous call to
                                          !! controlled_forcing_init, it will be
                                          !! deallocated here.

end subroutine controlled_forcing_end
  end interface

end module MOM_controlled_forcing
