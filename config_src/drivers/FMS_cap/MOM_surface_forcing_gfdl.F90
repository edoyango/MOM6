! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module MOM_surface_forcing_gfdl

!#CTRL# use MOM_controlled_forcing, only : apply_ctrl_forcing, register_ctrl_forcing_restarts
!#CTRL# use MOM_controlled_forcing, only : controlled_forcing_init, controlled_forcing_end
!#CTRL# use MOM_controlled_forcing, only : ctrl_forcing_CS
use MOM_coms,             only : reproducing_sum, field_chksum
use MOM_constants,        only : hlv, hlf
use MOM_coupler_types,    only : coupler_2d_bc_type, coupler_type_write_chksums
use MOM_coupler_types,    only : coupler_type_initialized, coupler_type_spawn
use MOM_coupler_types,    only : coupler_type_copy_data
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_SUBCOMPONENT
use MOM_data_override,    only : data_override_init, data_override
use MOM_diag_mediator,    only : diag_ctrl, safe_alloc_ptr, time_type
use MOM_domains,          only : pass_vector, pass_var, fill_symmetric_edges
use MOM_domains,          only : AGRID, BGRID_NE, CGRID_NE, To_All
use MOM_domains,          only : To_North, To_East, Omit_Corners
use MOM_EOS,              only : gsw_sr_from_sp
use MOM_error_handler,    only : MOM_error, WARNING, FATAL, is_root_pe, MOM_mesg
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,     only : forcing, mech_forcing
use MOM_forcing_type,     only : forcing_diags, mech_forcing_diags, register_forcing_type_diags
use MOM_forcing_type,     only : allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,     only : allocate_mech_forcing, deallocate_mech_forcing
use MOM_get_input,        only : Get_MOM_Input, directories
use MOM_grid,             only : ocean_grid_type
use MOM_interpolate,      only : init_external_field, time_interp_external
use MOM_interpolate,      only : time_interp_external_init
use MOM_interpolate,      only : external_field
use MOM_io,               only : slasher, write_version_number, MOM_read_data
use MOM_io,               only : read_netCDF_data
use MOM_io,               only : stdout_if_root
use MOM_restart,          only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,          only : restart_init_end, save_restart, restore_state
use MOM_string_functions, only : uppercase
use MOM_spatial_means,    only : adjust_area_mean_to_zero
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : surface
use user_revise_forcing,  only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,  only : user_revise_forcing_CS
use iso_fortran_env, only : int64

implicit none ; private

#include <MOM_memory.h>

public convert_IOB_to_fluxes, convert_IOB_to_forces
public surface_forcing_init
public ice_ocn_bnd_type_chksum
public forcing_save_restart

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> surface_forcing_CS is a structure containing pointers to the forcing fields
!! which may be used to drive MOM.  All fluxes are positive downward.
type, public :: surface_forcing_CS ; private
  integer :: wind_stagger       !< AGRID, BGRID_NE, or CGRID_NE (integer values
                                !! from MOM_domains) to indicate the staggering of
                                !! the winds that are being provided in calls to
                                !! update_ocean_model.
  logical :: use_temperature    !< If true, temp and saln used as state variables.
  logical :: nonBous            !< If true, this run is fully non-Boussinesq
  real :: wind_stress_multiplier !< A multiplier applied to incoming wind stress [nondim].

  real :: Rho0                  !< Boussinesq reference density [R ~> kg m-3]
  real :: area_surf = -1.0      !< Total ocean surface area [L2 ~> m2]
  real :: latent_heat_fusion    !< Latent heat of fusion [Q ~> J kg-1]
  real :: latent_heat_vapor     !< Latent heat of vaporization [Q ~> J kg-1]

  real :: max_p_surf            !< The maximum surface pressure that can be exerted by
                                !! the atmosphere and floating sea-ice [R L2 T-2 ~> Pa].
                                !! This is needed because the FMS coupling structure
                                !! does not limit the water that can be frozen out
                                !! of the ocean and the ice-ocean heat fluxes are
                                !! treated explicitly.
  logical :: use_limited_P_SSH  !< If true, return the sea surface height with
                                !! the correction for the atmospheric (and sea-ice)
                                !! pressure limited by max_p_surf instead of the
                                !! full atmospheric pressure.  The default is true.
  logical :: approx_net_mass_src !< If true, use the net mass sources from the ice-ocean boundary
                                !! type without any further adjustments to drive the ocean dynamics.
                                !! The actual net mass source may differ due to corrections.

  real :: gust_const            !< Constant unresolved background gustiness for ustar [R Z2 T-2 ~> Pa]
  logical :: read_gust_2d       !< If true, use a 2-dimensional gustiness supplied from an input file.
  real, pointer, dimension(:,:) :: &
    BBL_tidal_dis => NULL()     !< Tidal energy dissipation in the bottom boundary layer that can act as a
                                !! source of energy for bottom boundary layer mixing [R Z L2 T-3 ~> W m-2]
  real, pointer, dimension(:,:) :: &
    gust => NULL()              !< A spatially varying unresolved background gustiness that
                                !! contributes to ustar [R Z2 T-2 ~> Pa].  gust is used when read_gust_2d is true.
  real, pointer, dimension(:,:) :: &
    ustar_tidal => NULL()       !< Tidal contribution to the bottom friction velocity [Z T-1 ~> m s-1]
  real :: cd_tides              !< Drag coefficient that applies to the tides [nondim]
  real :: utide                 !< Constant tidal velocity to use if read_tideamp is false [Z T-1 ~> m s-1].
  logical :: read_tideamp       !< If true, spatially varying tidal amplitude read from a file.

  logical :: rigid_sea_ice      !< If true, sea-ice exerts a rigidity that acts to damp surface
                                !! deflections (especially surface gravity waves).  The default is false.
  real    :: g_Earth            !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real    :: Kv_sea_ice         !< Viscosity in sea-ice that resists sheared vertical motions [L4 Z-2 T-1 ~> m2 s-1]
  real    :: density_sea_ice    !< Typical density of sea-ice [R ~> kg m-3]. The value is only used to convert
                                !! the ice pressure into appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass !< A mass per unit area of sea-ice beyond which sea-ice viscosity
                                !! becomes effective [R Z ~> kg m-2], typically of order 1000 kg m-2.
  logical :: allow_flux_adjustments !< If true, use data_override to obtain flux adjustments
  logical :: allow_carbon_flux_exchange !< If true, allows fluxes and diagnostics of carbon in runoff.

  logical :: restore_salt       !< If true, the coupled MOM driver adds a term to restore surface
                                !! salinity to a specified value.
  logical :: restore_temp       !< If true, the coupled MOM driver adds a term to restore sea
                                !! surface temperature to a specified value.
  real    :: Flux_const_salt    !< Piston velocity for surface salinity restoring [Z T-1 ~> m s-1]
  real    :: Flux_const_temp    !< Piston velocity for surface temperature restoring [Z T-1 ~> m s-1]
  real    :: rho_restore        !< The density that is used to convert piston velocities into salt
                                !! or heat fluxes with salinity or temperature restoring [R ~> kg m-3]
  logical :: trestore_SPEAR_ECDA            !< If true, modify restoring data wrt local SSS
  real    :: SPEAR_dTf_dS                   !< The derivative of the freezing temperature with
                                            !! salinity [C S-1 ~> degC ppt-1].
  logical :: salt_restore_as_sflux          !< If true, SSS restore as salt flux instead of water flux
  logical :: adjust_net_srestore_to_zero    !< Adjust srestore to zero (for both salt_flux or vprec)
  logical :: adjust_net_srestore_by_scaling !< Adjust srestore w/o moving zero contour
  logical :: adjust_net_fresh_water_to_zero !< Adjust net surface fresh-water (with restoring) to zero
  logical :: use_net_FW_adjustment_sign_bug !< Use the wrong sign when adjusting net FW
  logical :: adjust_net_fresh_water_by_scaling !< Adjust net surface fresh-water w/o moving zero contour
  logical :: mask_srestore_under_ice        !< If true, use an ice mask defined by frazil criteria
                                            !! for salinity restoring.
  real    :: ice_salt_concentration         !< Salt concentration for sea ice [kg/kg]
  logical :: mask_srestore_marginal_seas    !< If true, then mask SSS restoring in marginal seas
  real    :: max_delta_srestore             !< Maximum delta salinity used for restoring [S ~> ppt]
  real    :: max_delta_trestore             !< Maximum delta sst used for restoring [C ~> degC]
  real, pointer, dimension(:,:) :: basin_mask => NULL() !< Mask for surface salinity restoring by basin [nondim]
  integer :: answer_date        !< The vintage of the order of arithmetic and expressions in the
                                !! gustiness calculations.  Values below 20190101 recover the answers
                                !! from the end of 2018, while higher values use a simpler expression
                                !! to calculate gustiness.
  logical :: ustar_gustless_bug             !< If true, include a bug in the time-averaging of the
                                            !! gustless wind friction velocity.
  logical :: check_no_land_fluxes           !< Return warning if IOB flux over land is non-zero

  type(diag_ctrl), pointer :: diag => NULL()  !< Structure to regulate diagnostic output timing
  character(len=200) :: inputdir              !< Directory where NetCDF input files are
  character(len=200) :: salt_restore_file     !< Filename for salt restoring data
  character(len=30)  :: salt_restore_var_name !< Name of surface salinity in salt_restore_file
  logical            :: salt_restore_is_practical !< Specifies that the target salinity is practical and not absolute.
  logical            :: mask_srestore         !< If true, apply a 2-dimensional mask to the surface
                                              !! salinity restoring fluxes. The masking file should be
                                              !! in inputdir/salt_restore_mask.nc and the field should
                                              !! be named 'mask'
  real, pointer, dimension(:,:) :: srestore_mask => NULL() !< mask for SSS restoring [nondim]
  character(len=200) :: temp_restore_file     !< Filename for sst restoring data
  character(len=30)  :: temp_restore_var_name !< Name of surface temperature in temp_restore_file
  logical            :: mask_trestore         !< If true, apply a 2-dimensional mask to the surface
                                              !! temperature restoring fluxes. The masking file should be
                                              !! in inputdir/temp_restore_mask.nc and the field should
                                              !! be named 'mask'
  real, pointer, dimension(:,:) :: trestore_mask => NULL() !< Mask for SST restoring [nondim]
  type(external_field) :: srestore_handle
    !< Handle for time-interpolated salt restoration field
  type(external_field) :: trestore_handle
    !< Handle for time-interpolated temperature restoration field

  type(forcing_diags), public :: handles !< Diagnostics handles

!#CTRL#  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL() !< A control structure for user forcing revisions
end type surface_forcing_CS


!> ice_ocean_boundary_type is a structure corresponding to forcing, but with the elements, units,
!! and conventions that exactly conform to the use for MOM6-based coupled models.
type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: u_flux          =>NULL() !< i-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: v_flux          =>NULL() !< j-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: t_flux          =>NULL() !< sensible heat flux [W m-2]
  real, pointer, dimension(:,:) :: q_flux          =>NULL() !< specific humidity flux [kg m-2 s-1]
  real, pointer, dimension(:,:) :: salt_flux       =>NULL() !< salt flux [kg m-2 s-1]
  real, pointer, dimension(:,:) :: excess_salt     =>NULL() !< salt left behind by brine rejection [kg m-2 s-1]
  real, pointer, dimension(:,:) :: lw_flux         =>NULL() !< long wave radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dir =>NULL() !< direct visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dif =>NULL() !< diffuse visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dir =>NULL() !< direct Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dif =>NULL() !< diffuse Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: lprec           =>NULL() !< mass flux of liquid precip [kg m-2 s-1]
  real, pointer, dimension(:,:) :: fprec           =>NULL() !< mass flux of frozen precip [kg m-2 s-1]
  real, pointer, dimension(:,:) :: runoff          =>NULL() !< mass flux of liquid runoff [kg m-2 s-1]
  real, pointer, dimension(:,:) :: runoff_carbon   =>NULL() !< mass flux of carbon in liquid runoff [kg m-2 s-1]
  real, pointer, dimension(:,:) :: calving         =>NULL() !< mass flux of frozen runoff [kg m-2 s-1]
  real, pointer, dimension(:,:) :: stress_mag      =>NULL() !< The time-mean magnitude of the stress on the ocean [Pa]
  real, pointer, dimension(:,:) :: ustar_berg      =>NULL() !< frictional velocity beneath icebergs [m s-1]
  real, pointer, dimension(:,:) :: area_berg       =>NULL() !< fractional area covered by icebergs [m2 m-2]
  real, pointer, dimension(:,:) :: mass_berg       =>NULL() !< mass of icebergs per unit ocean area [kg m-2]
  real, pointer, dimension(:,:) :: runoff_hflx     =>NULL() !< heat content of liquid runoff [W m-2]
  real, pointer, dimension(:,:) :: calving_hflx    =>NULL() !< heat content of frozen runoff [W m-2]
  real, pointer, dimension(:,:) :: p               =>NULL() !< pressure of overlying ice and atmosphere
                                                            !< on ocean surface [Pa]
  real, pointer, dimension(:,:) :: mi              =>NULL() !< mass of ice per unit ocean area [kg m-2]
  real, pointer, dimension(:,:) :: ice_rigidity    =>NULL() !< rigidity of the sea ice, sea-ice and
                                                            !! ice-shelves, expressed as a coefficient
                                                            !! for divergence damping, as determined
                                                            !! outside of the ocean model [m3 s-1]
  real, pointer, dimension(:,:) :: shelf_sfc_mass_flux =>NULL() !< mass flux to surface of ice sheet [kg m-2 s-1]
  integer :: xtype                    !< The type of the exchange - REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type) :: fluxes  !< A structure that may contain an array of named fields
                                      !! used for passive tracer fluxes.
  integer :: wind_stagger = -999      !< A flag indicating the spatial discretization of wind stresses.
                                      !! This flag may be set by the flux-exchange code, based on what
                                      !! the sea-ice model is providing.  Otherwise, the value from
                                      !! the surface_forcing_CS is used.
end type ice_ocean_boundary_type

integer :: id_clock_forcing !< A CPU time clock


  interface
module subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, valid_time, G, US, CS, sfc_state)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
                                                   !! the ocean in a coupled model
  type(forcing),           intent(inout) :: fluxes !< A structure containing pointers to all
                                                   !! possible mass, heat or salt flux forcing fields.
                                                   !!  Unused fields have NULL ptrs.
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time   !< The time of the fluxes, used for interpolating the
                                                   !! salinity to the right time, when it is being restored.
  real,                    intent(in)    :: valid_time !< The amount of time over which these fluxes
                                                   !! should be applied [T ~> s].
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS     !< A pointer to the control structure returned by a
                                                   !! previous call to surface_forcing_init.
  type(surface),           intent(in)    :: sfc_state !< A structure containing fields that describe the
                                                   !! surface state of the ocean.

end subroutine convert_IOB_to_fluxes
module subroutine convert_IOB_to_forces(IOB, forces, index_bounds, Time, G, US, CS, dt_forcing, reset_avg)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
                                                   !! the ocean in a coupled model
  type(mech_forcing),      intent(inout) :: forces !< A structure with the driving mechanical forces
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time   !< The time of the fluxes, used for interpolating the
                                                   !! salinity to the right time, when it is being restored.
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS     !< A pointer to the control structure returned by a
                                                   !! previous call to surface_forcing_init.
  real,          optional, intent(in)    :: dt_forcing !< A time interval over which to apply the
                                                   !! current value of ustar as a weighted running
                                                   !! average [T ~> s], or if 0 do not average ustar.
                                                   !! Missing is equivalent to 0.
  logical,       optional, intent(in)    :: reset_avg !< If true, reset the time average.

  ! Local variables

                              ! mass fluxes [R Z s m2 kg-1 T-1 ~> 1]


end subroutine convert_IOB_to_forces
module subroutine extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux, tauy, ustar, &
                                gustless_ustar, mag_tau, gustless_mag_tau, tau_halo)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB  !< An ice-ocean boundary type with fluxes to drive
                                                 !! the ocean in a coupled model
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time !< The time of the fluxes, used for interpolating the
                                                 !! salinity to the right time, when it is being restored.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS   !< A pointer to the control structure returned by a
                                                 !! previous call to surface_forcing_init.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(inout) :: taux !< The zonal wind stresses on a C-grid [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(inout) :: tauy !< The meridional wind stresses on a C-grid [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(inout) :: ustar !< The surface friction velocity [Z T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out)   :: gustless_ustar !< The surface friction velocity without
                                                 !! any contributions from gustiness [Z T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(inout) :: mag_tau !< The magintude of the wind stress at tracer points
                                                 !! including subgridscale variability and gustiness [R Z2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out) :: gustless_mag_tau !< The magintude of the wind stress at tracer points
                                                 !! without any contributions from gustiness [R Z2 T-2 ~> Pa]
  integer,       optional, intent(in)    :: tau_halo !< The halo size of wind stresses to set, 0 by default.

  ! Local variables



end subroutine extract_IOB_stresses
module subroutine apply_flux_adjustments(G, US, CS, Time, fluxes)
  type(ocean_grid_type),    intent(inout) :: G  !< Ocean grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(forcing),            intent(inout) :: fluxes !< Surface fluxes structure

  ! Local variables
                                                 ! [Q R Z T-1 ~> W m-2] or [R Z T-1 ~> kg m-2 s-1]


end subroutine apply_flux_adjustments
module subroutine apply_force_adjustments(G, US, CS, Time, forces)
  type(ocean_grid_type),    intent(inout) :: G  !< Ocean grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces

  ! Local variables


end subroutine apply_force_adjustments
module subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< A pointer to the control structure returned
                                                    !! by a previous call to surface_forcing_init
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !< The current model time
  character(len=*),           intent(in)    :: directory !< The directory into which to write the
                                                    !! restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file names include
                                                    !! a unique time stamp.  The default is false.
  character(len=*), optional, intent(in)    :: filename_suffix !< An optional suffix (e.g., a time-
                                                    !! stamp) to append to the restart file names.

end subroutine forcing_save_restart
module subroutine surface_forcing_init(Time, G, US, param_file, diag, CS, wind_stagger)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure that is used to regulate
                                                  !! diagnostic output
  type(surface_forcing_CS), pointer       :: CS   !< A pointer that is set to point to the control
                                                  !! structure for this module
  integer,        optional, intent(in)    :: wind_stagger !< If present, the staggering of the winds
                                                  !! that are being provided in calls to update_ocean_model

  ! Local variables
                            ! the tidal bottom TKE input used with INT_TIDE_DISSIPATION, times the
                            ! factor rescaling from the units of TKE to those of mean kinetic
                            ! energy [R L2 Z-2 ~> kg m-3]
                                  ! or other equivalent files.
  ! This include declares and sets the variable "version".

end subroutine surface_forcing_init
module subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS !< A pointer to the control structure returned by
                                                !! a previous call to surface_forcing_init, it will
                                                !! be deallocated here.
  type(forcing), optional,  intent(inout) :: fluxes !< A structure containing pointers to all
                                                !! possible mass, heat or salt flux forcing fields.
                                                !! If present, it will be deallocated here.

end subroutine surface_forcing_end
module subroutine ice_ocn_bnd_type_chksum(id, timestep, iobt)

  character(len=*), intent(in) :: id     !< An identifying string for this call
  integer,          intent(in) :: timestep !< The number of elapsed timesteps
  type(ice_ocean_boundary_type), &
                    intent(in) :: iobt   !< An ice-ocean boundary type with fluxes to drive the
                                         !! ocean in a coupled model whose checksums are reported
  ! Local variables

end subroutine ice_ocn_bnd_type_chksum
module subroutine check_mask_val_consistency(val, mask, i, j, varname, G)

  real, intent(in) :: val  !< value of flux/variable passed by IOB [various]
  real, intent(in) :: mask !< value of ocean mask [nondim]
  integer, intent(in) :: i !< model grid cell indices
  integer, intent(in) :: j !< model grid cell indices
  character(len=*), intent(in) :: varname !< variable name
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure
  ! Local variables

end subroutine check_mask_val_consistency
  end interface

end module MOM_surface_forcing_gfdl
