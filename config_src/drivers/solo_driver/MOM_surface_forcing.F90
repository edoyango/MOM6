! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Functions that calculate the surface wind stresses and fluxes of buoyancy
!! or temperature/salinity and fresh water, in ocean-only (solo) mode.
!!
!! These functions are called every time step, even if the wind stresses
!! or buoyancy fluxes are constant in time - in that case these routines
!! return quickly without doing anything.  In addition, any I/O of forcing
!! fields is controlled by surface_forcing_init, located in this file.
module MOM_surface_forcing

use MOM_constants,           only : hlv, hlf
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE
use MOM_data_override,       only : data_override_init, data_override
use MOM_diag_mediator,       only : post_data, query_averaging_enabled
use MOM_diag_mediator,       only : diag_ctrl, safe_alloc_ptr
use MOM_domains,             only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_domains,             only : fill_symmetric_edges, CGRID_NE
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_string_functions,    only : uppercase
use MOM_forcing_type,        only : forcing, mech_forcing
use MOM_forcing_type,        only : set_net_mass_forcing, copy_common_forcing_fields
use MOM_forcing_type,        only : set_derived_forcing_fields
use MOM_forcing_type,        only : forcing_diags, mech_forcing_diags, register_forcing_type_diags
use MOM_forcing_type,        only : allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,        only : allocate_mech_forcing, deallocate_mech_forcing
use MOM_grid,                only : ocean_grid_type
use MOM_get_input,           only : Get_MOM_Input, directories
use MOM_io,                  only : file_exists, MOM_read_data, MOM_read_vector, slasher
use MOM_io,                  only : read_netCDF_data, EAST_FACE, NORTH_FACE, num_timelevels
use MOM_restart,             only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,             only : restart_init_end, save_restart, restore_state
use MOM_time_manager,        only : time_type, operator(+), operator(/), operator(*)
use MOM_time_manager,        only : set_time, get_time, get_date, time_to_real
use MOM_tracer_flow_control, only : call_tracer_set_forcing, tracer_flow_control_CS
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : surface
use MESO_surface_forcing,    only : MESO_buoyancy_forcing
use MESO_surface_forcing,    only : MESO_surface_forcing_init, MESO_surface_forcing_CS
use user_surface_forcing,    only : USER_wind_forcing, USER_buoyancy_forcing
use user_surface_forcing,    only : USER_surface_forcing_init, user_surface_forcing_CS
use user_revise_forcing,     only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,     only : user_revise_forcing_CS
use idealized_hurricane,     only : idealized_hurricane_wind_forcing
use idealized_hurricane,     only : idealized_hurricane_wind_init, idealized_hurricane_CS
use SCM_CVmix_tests,         only : SCM_CVmix_tests_surface_forcing_init
use SCM_CVmix_tests,         only : SCM_CVmix_tests_wind_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_buoyancy_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_CS
use BFB_surface_forcing,     only : BFB_buoyancy_forcing
use BFB_surface_forcing,     only : BFB_surface_forcing_init, BFB_surface_forcing_CS
use dumbbell_surface_forcing, only : dumbbell_surface_forcing_init, dumbbell_surface_forcing_CS
use dumbbell_surface_forcing, only : dumbbell_buoyancy_forcing
use MARBL_forcing_mod,       only : marbl_forcing_CS, MARBL_forcing_init
use MARBL_forcing_mod,       only : convert_driver_fields_to_forcings

implicit none ; private

#include <MOM_memory.h>

public set_forcing
public surface_forcing_init
public forcing_save_restart

!> Structure containing pointers to the forcing fields that may be used to drive MOM.
!!  All fluxes are positive into the ocean.
type, public :: surface_forcing_CS ; private

  logical :: use_temperature    !< if true, temp & salinity used as state variables
  logical :: restorebuoy        !< if true, use restoring surface buoyancy forcing
  logical :: adiabatic          !< if true, no diapycnal mass fluxes or surface buoyancy forcing
  logical :: nonBous            !< If true, this run is fully non-Boussinesq
  logical :: variable_winds     !< if true, wind stresses vary with time
  logical :: variable_buoyforce !< if true, buoyancy forcing varies with time.
  real    :: south_lat          !< southern latitude of the domain [degrees_N] or [km] or [m]
  real    :: len_lat            !< domain length in latitude [degrees_N] or [km] or [m]

  real :: Rho0                  !< Boussinesq reference density [R ~> kg m-3]
  real :: G_Earth               !< gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const = 0.0      !< piston velocity for surface restoring [Z T-1 ~> m s-1]
  real :: Flux_const_T = 0.0    !< piston velocity for surface temperature restoring [Z T-1 ~> m s-1]
  real :: Flux_const_S = 0.0    !< piston velocity for surface salinity restoring [Z T-1 ~> m s-1]
  real :: rho_restore           !< The density that is used to convert piston velocities into salt
                                !! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]
  real :: latent_heat_fusion    !< latent heat of fusion times [Q ~> J kg-1]
  real :: latent_heat_vapor     !< latent heat of vaporization [Q ~> J kg-1]
  real :: tau_x0                !< Constant zonal wind stress used in the WIND_CONFIG="const"
                                !! forcing [R L Z T-2 ~> Pa]
  real :: tau_y0                !< Constant meridional wind stress used in the WIND_CONFIG="const"
                                !! forcing [R L Z T-2 ~> Pa]
  real :: taux_mag              !< Peak magnitude of the zonal wind stress for several analytic
                                !! profiles [R L Z T-2 ~> Pa]

  real    :: gust_const                 !< constant unresolved background gustiness for ustar [R Z2 T-2 ~> Pa]
  logical :: read_gust_2d               !< if true, use 2-dimensional gustiness supplied from a file
  real, pointer :: gust(:,:) => NULL()  !< spatially varying unresolved background gustiness [R L Z T-2 ~> Pa]
                                        !! gust is used when read_gust_2d is true.

  real, pointer :: T_Restore(:,:)    => NULL()  !< temperature to damp (restore) the SST to [C ~> degC]
  real, pointer :: S_Restore(:,:)    => NULL()  !< salinity to damp (restore) the SSS [S ~> ppt]
  real, pointer :: Dens_Restore(:,:) => NULL()  !< density to damp (restore) surface density [R ~> kg m-3]

  ! if WIND_CONFIG=='gyres' then use the following as  = A, B, C and n respectively for
  ! taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L)
  real :: gyres_taux_const   !< A constant wind stress [R L Z T-2 ~> Pa].
  real :: gyres_taux_sin_amp !< The amplitude of cosine wind stress gyres [R L Z T-2 ~> Pa], if WIND_CONFIG=='gyres'
  real :: gyres_taux_cos_amp !< The amplitude of cosine wind stress gyres [R L Z T-2 ~> Pa], if WIND_CONFIG=='gyres'
  real :: gyres_taux_n_pis   !< The number of sine lobes in the basin if WIND_CONFIG=='gyres' [nondim]
  integer :: answer_date     !< This 8-digit integer gives the approximate date with which the order
                             !! of arithmetic and expressions were added to the code.
                             !! Dates before 20190101 use original answers.
                             !! Dates after 20190101 use a form of the gyre wind stresses that are
                             !! rotationally invariant and more likely to be the same between compilers.
  logical :: ustar_gustless_bug   !< If true, include a bug in the time-averaging of the
                                  !! gustless wind friction velocity.
  logical :: use_marbl_tracers              !< If true, allocate memory for forcing needed by MARBL
  ! if WIND_CONFIG=='scurves' then use the following to define a piecewise scurve profile
  real :: scurves_ydata(20) = 90. !< Latitudes of scurve nodes [degreesN]
  real :: scurves_taux(20) = 0.   !< Zonal wind stress values at scurve nodes [R L Z T-2 ~> Pa]

  real :: T_north   !< Target temperatures at north used in buoyancy_forcing_linear [C ~> degC]
  real :: T_south   !< Target temperatures at south used in buoyancy_forcing_linear [C ~> degC]
  real :: S_north   !< Target salinity at north used in buoyancy_forcing_linear [S ~> ppt]
  real :: S_south   !< Target salinity at south used in buoyancy_forcing_linear [S ~> ppt]

  logical :: first_call_set_forcing = .true. !< True until after the first call to set_forcing
  logical :: archaic_OMIP_file = .true. !< If true use the variable names and data fields from
                                        !! a very old version of the OMIP forcing
  logical :: dataOverrideIsInitialized = .false. !< If true, data override has been initialized

  real :: wind_scale          !< value by which wind-stresses are scaled [nondim]
  real :: constantHeatForcing !< value used for sensible heat flux when buoy_config="const" [Q R Z T-1 ~> W m-2]

  character(len=8)   :: wind_stagger !< A character indicating how the wind stress components
                              !! are staggered in WIND_FILE.  Valid values are A or C for now.
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL() !< A pointer to the structure
                              !! that is used to orchestrate the calling of tracer packages
!#CTRL#  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(diag_ctrl), pointer :: diag !< structure used to regulate timing of diagnostic output

  character(len=200) :: inputdir    !< directory where NetCDF input files are.
  character(len=200) :: wind_config !< indicator for wind forcing type (2gyre, USER, FILE..)
  character(len=200) :: wind_file   !< if wind_config is "file", file to use
  character(len=200) :: buoy_config !< indicator for buoyancy forcing type

  character(len=200) :: longwave_file     = '' !< The file from which the longwave heat flux is read
  character(len=200) :: shortwave_file    = '' !< The file from which the shortwave heat flux is read
  character(len=200) :: evaporation_file  = '' !< The file from which the evaporation is read
  character(len=200) :: sensibleheat_file = '' !< The file from which the sensible heat flux is read
  character(len=200) :: latentheat_file   = '' !< The file from which the latent heat flux is read

  character(len=200) :: rain_file   = '' !< The file from which the rainfall is read
  character(len=200) :: snow_file   = '' !< The file from which the snowfall is read
  character(len=200) :: runoff_file = '' !< The file from which the runoff is read

  character(len=200) :: longwaveup_file  = '' !< The file from which the upward longwave heat flux is read
  character(len=200) :: shortwaveup_file = '' !< The file from which the upward shortwave heat flux is read

  character(len=200) :: SSTrestore_file      = '' !< The file from which to read the sea surface
                                                  !! temperature to restore toward
  character(len=200) :: salinityrestore_file = '' !< The file from which to read the sea surface
                                                  !! salinity to restore toward

  character(len=80)  :: stress_x_var  = '' !< X-wind stress variable name in the input file
  character(len=80)  :: stress_y_var  = '' !< Y-wind stress variable name in the input file
  character(len=80)  :: ustar_var     = '' !< ustar variable name in the input file
  character(len=80)  :: LW_var        = '' !< longwave heat flux variable name in the input file
  character(len=80)  :: SW_var        = '' !< shortwave heat flux variable name in the input file
  character(len=80)  :: latent_var    = '' !< latent heat flux variable name in the input file
  character(len=80)  :: sens_var      = '' !< sensible heat flux variable name in the input file
  character(len=80)  :: evap_var      = '' !< evaporation variable name in the input file
  character(len=80)  :: rain_var      = '' !< rainfall variable name in the input file
  character(len=80)  :: snow_var      = '' !< snowfall variable name in the input file
  character(len=80)  :: lrunoff_var   = '' !< liquid runoff variable name in the input file
  character(len=80)  :: frunoff_var   = '' !< frozen runoff variable name in the input file
  character(len=80)  :: SST_restore_var = '' !< target sea surface temperature variable name in the input file
  character(len=80)  :: SSS_restore_var = '' !< target sea surface salinity variable name in the input file

  ! These variables relate model times to time levels in the various forcing files.
  integer :: wind_days_per_rec = 0    !< If positive the number of days of wind stress per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: SW_days_per_rec = 0      !< If positive the number of days shortwave heat flux per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: LW_days_per_rec = 0      !< If positive the number of days longwave heat flux per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: latent_days_per_rec = 0  !< If positive the number of days latent heat flux per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: sens_days_per_rec = 0    !< If positive the number of days sensible heat flux per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: evap_days_per_rec = 0    !< If positive the number of days evaporation per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: precip_days_per_rec = 0  !< If positive the number of days precipitation per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: runoff_days_per_rec = 0  !< If positive the number of days runoff per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: SST_days_per_rec = 0     !< If positive the number of days target SST per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.
  integer :: SSS_days_per_rec = 0     !< If positive the number of days target SSS per forcing file time level,
                                      !! or if negative the number of time levels per day.  If 31 change forcing
                                      !! monthly, or if 0 the model will guess the right value based on the file size.

  ! These variables give the number of time levels in the various forcing files.
  integer :: wind_nlev = -1   !< The number of time levels in the file of wind stress
  integer :: SW_nlev   = -1   !< The number of time levels in the file of shortwave heat flux
  integer :: LW_nlev = -1     !< The number of time levels in the file of longwave heat flux
  integer :: latent_nlev = -1 !< The number of time levels in the file of latent heat flux
  integer :: sens_nlev = -1   !< The number of time levels in the file of sensible heat flux
  integer :: evap_nlev = -1   !< The number of time levels in the file of evaporation
  integer :: precip_nlev = -1 !< The number of time levels in the file of precipitation
  integer :: runoff_nlev = -1 !< The number of time levels in the file of runoff
  integer :: SST_nlev  = -1   !< The number of time levels in the file of target SST
  integer :: SSS_nlev = -1    !< The number of time levels in the file of target SSS

  ! These variables give the last time level read for the various forcing files.
  integer :: wind_last_lev = -1   !< The last time level read of wind stress
  integer :: SW_last_lev   = -1   !< The last time level read of shortwave heat flux
  integer :: LW_last_lev = -1     !< The last time level read of longwave heat flux
  integer :: latent_last_lev = -1 !< The last time level read of latent heat flux
  integer :: sens_last_lev = -1   !< The last time level read of sensible heat flux
  integer :: evap_last_lev = -1   !< The last time level read of evaporation
  integer :: precip_last_lev = -1 !< The last time level read of precipitation
  integer :: runoff_last_lev = -1 !< The last time level read of runoff
  integer :: SST_last_lev  = -1   !< The last time level read of target SST
  integer :: SSS_last_lev = -1    !< The last time level read of target SSS

  type(forcing_diags), public :: handles !< A structure with diagnostics handles

  !>@{ Control structures for named forcing packages
  type(user_revise_forcing_CS),  pointer :: urf_CS => NULL()
  type(user_surface_forcing_CS), pointer :: user_forcing_CSp => NULL()
  type(BFB_surface_forcing_CS), pointer :: BFB_forcing_CSp => NULL()
  type(dumbbell_surface_forcing_CS), pointer :: dumbbell_forcing_CSp => NULL()
  type(MESO_surface_forcing_CS), pointer :: MESO_forcing_CSp => NULL()
  type(idealized_hurricane_CS), pointer :: idealized_hurricane_CSp => NULL()
  type(SCM_CVmix_tests_CS),      pointer :: SCM_CVmix_tests_CSp => NULL()
  type(marbl_forcing_CS), pointer :: marbl_forcing_CSp => NULL()
  !>@}

end type surface_forcing_CS

integer :: id_clock_forcing !< A CPU time clock


  interface
module subroutine set_forcing(sfc_state, forces, fluxes, day_start, day_interval, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),    intent(inout) :: forces !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day_start !< The start time of the fluxes
  type(time_type),       intent(in)    :: day_interval !< Length of time over which these fluxes applied
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
end subroutine set_forcing
module subroutine wind_forcing_const(sfc_state, forces, tau_x0, tau_y0, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  real,                     intent(in)    :: tau_x0 !< The zonal wind stress [R Z L T-2 ~> Pa]
  real,                     intent(in)    :: tau_y0 !< The meridional wind stress [R Z L T-2 ~> Pa]
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables

end subroutine wind_forcing_const
module subroutine wind_forcing_2gyre(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables

end subroutine wind_forcing_2gyre
module subroutine wind_forcing_1gyre(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables

end subroutine wind_forcing_1gyre
module subroutine wind_forcing_gyres(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables

end subroutine wind_forcing_gyres
module subroutine Neverworld_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure.
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control structure returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables

end subroutine Neverworld_wind_forcing
module subroutine scurve_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure.
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control structure returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables
! real :: ydata(7) = (/ -70., -45., -15., 0., 15., 45., 70. /)
! real :: taudt(7) = (/ 0., 0.2, -0.1, -0.02, -0.1, 0.1, 0. /)

  ! Allocate the forcing arrays, if necessary.
end subroutine scurve_wind_forcing
real module function scurve(x,L)
  real , intent(in) :: x       !< non-dimensional position [nondim]
  real , intent(in) :: L       !< non-dimensional width [nondim]

end function scurve
module subroutine wind_forcing_from_file(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
                        ! sub-gridscale variability or gustiness [R Z2 T-2 ~> Pa]

end subroutine wind_forcing_from_file
module subroutine wind_forcing_by_data_override(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
                        ! sub-gridscale variability or gustiness [R Z2 T-2 ~> Pa]

end subroutine wind_forcing_by_data_override
module subroutine stresses_to_ustar(forces, G, US, CS)
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(ocean_grid_type),    intent(in)    :: G      !< Grid structure.
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control structure returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables
                        ! sub-gridscale variability or gustiness [R Z2 T-2 ~> Pa]

end subroutine stresses_to_ustar
module subroutine buoyancy_forcing_from_files(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [T ~> s]
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
                  ! [R Z T-1 ~> kg m-2 s-1]
!#CTRL#  real, dimension(SZI_(G),SZJ_(G)) :: &
!#CTRL#    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
!#CTRL#                  ! target (observed) value [C ~> degC].
!#CTRL#    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
!#CTRL#                  ! (observed) value [S ~> ppt].
!#CTRL#    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
!#CTRL#                  ! anomalies when calculating restorative precipitation anomalies [S ~> ppt].




end subroutine buoyancy_forcing_from_files
module subroutine buoyancy_forcing_from_data_override(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                  !! describe the surface state of the ocean.
  type(forcing),            intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  real,                     intent(in)    :: dt   !< The amount of time over which
                                                  !! the fluxes apply [T ~> s]
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control structure returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
!#CTRL#  real, dimension(SZI_(G),SZJ_(G)) :: &
!#CTRL#    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
!#CTRL#                  ! target (observed) value [C ~> degC].
!#CTRL#    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
!#CTRL#                  ! (observed) value [S ~> ppt].
!#CTRL#    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
!#CTRL#                  ! anomalies when calculating restorative precipitation anomalies [S ~> ppt].

end subroutine buoyancy_forcing_from_data_override
module subroutine buoyancy_forcing_zero(sfc_state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [T ~> s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  ! Local variables

end subroutine buoyancy_forcing_zero
module subroutine buoyancy_forcing_const(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [T ~> s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
end subroutine buoyancy_forcing_const
module subroutine buoyancy_forcing_linear(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [T ~> s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  ! Local variables

end subroutine buoyancy_forcing_linear
module function get_file_time_level(Time, nlev_file, days_per_rec) result (time_lev)
  type(time_type), intent(in) :: Time         !< The time of the fluxes
  integer ,        intent(in) :: nlev_file    !< The number of time records in a forcing file
  integer ,        intent(in) :: days_per_rec !< If positive, the number of days spanned by each
                                              !! time record in a file, if negative the number
                                              !! time records per day, or if 0 determine this
                                              !! by guessing based on the number of records in
                                              !! the file.  If this is 31, the time levels will
                                              !! be based on the months of the calendar.
                                              !! Setting this larger than 1000000 will always
                                              !! cause the time level to be set to 1.
  integer :: time_lev                !< The time level in a file that will be read.

  ! Local variables
                                     ! taking the periodicity of the file into account.

end function get_file_time_level
module subroutine MARBL_forcing_from_data_override(fluxes, day, G, US, CS)
  type(forcing),            intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),          intent(in)    :: day    !< The time of the fluxes
  type(ocean_grid_type),    intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control structure returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]
                                                                 !! [R Z T-1 ~> kg m-2 s-1]

  ! Necessary null pointers for arguments to convert_driver_fields_to_forcings()
  ! Since they are null, MARBL will not use multiple ice categories

end subroutine MARBL_forcing_from_data_override
module subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< pointer to control structure returned by
                                                    !! a previous surface_forcing_init call
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !< model time at this call; needed for mpp_write calls
  character(len=*),           intent(in)    :: directory !< directory into which to write these restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file names
                                                    !! include a unique time stamp; the  default is false.
  character(len=*), optional, intent(in)    :: filename_suffix !< optional suffix (e.g., a time-stamp)
                                                    !! to append to the restart file name

end subroutine forcing_save_restart
module subroutine surface_forcing_init(Time, G, US, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),              intent(in)    :: Time !< The current model time
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,      intent(inout) :: diag !< structure used to regulate diagnostic output
  type(surface_forcing_CS),     pointer       :: CS   !< pointer to control structure returned by
                                                      !! a previous surface_forcing_init call
  type(tracer_flow_control_CS), pointer       :: tracer_flow_CSp !< Forcing for tracers?

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine surface_forcing_init
module subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control structure returned by
                                               !! a previous surface_forcing_init call
  type(forcing), optional,  intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields

end subroutine surface_forcing_end
  end interface

end module MOM_surface_forcing
