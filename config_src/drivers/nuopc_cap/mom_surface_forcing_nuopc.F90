! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Converts the input ESMF data (import data) to a MOM-specific data type (surface_forcing_CS).
module MOM_surface_forcing_nuopc

use MOM_coms,             only : reproducing_sum, field_chksum
use MOM_constants,        only : hlv, hlf
use MOM_coupler_types,    only : coupler_2d_bc_type, coupler_type_write_chksums
use MOM_coupler_types,    only : coupler_type_initialized, coupler_type_spawn
use MOM_coupler_types,    only : coupler_type_copy_data
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_SUBCOMPONENT
use MOM_data_override,    only : data_override_init, data_override
use MOM_diag_mediator,    only : diag_ctrl
use MOM_diag_mediator,    only : safe_alloc_ptr, time_type
use MOM_domains,          only : pass_vector, pass_var, fill_symmetric_edges
use MOM_domains,          only : AGRID, BGRID_NE, CGRID_NE, To_All
use MOM_domains,          only : To_North, To_East, Omit_Corners
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
use MOM_io,               only : stdout
use MOM_restart,          only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,          only : restart_init_end, save_restart, restore_state
use MOM_string_functions, only : uppercase
use MOM_spatial_means,    only : adjust_area_mean_to_zero
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : surface
use user_revise_forcing,  only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,  only : user_revise_forcing_CS
use iso_fortran_env,      only : int64
use MARBL_forcing_mod,    only : marbl_forcing_CS, MARBL_forcing_init
use MARBL_forcing_mod,    only : convert_driver_fields_to_forcings

implicit none ; private

#include <MOM_memory.h>

public convert_IOB_to_fluxes
public convert_IOB_to_forces
public surface_forcing_init
public forcing_save_restart
public ice_ocn_bnd_type_chksum

private apply_flux_adjustments
private apply_force_adjustments
private surface_forcing_end

!> Contains pointers to the forcing fields which may be used to drive MOM.
!! All fluxes are positive downward.
type, public :: surface_forcing_CS ; private
  integer :: wind_stagger       !< AGRID, BGRID_NE, or CGRID_NE (integer values
                                !! from MOM_domains) to indicate the staggering of
                                !! the winds that are being provided in calls to
                                !! update_ocean_model.
  logical :: use_temperature    !! If true, temp and saln used as state variables
  real :: wind_stress_multiplier !< A multiplier applied to incoming wind stress (nondim).

  real :: Rho0                  !< Boussinesq reference density [R ~> kg m-3]
  real :: area_surf = -1.0      !< total ocean surface area [L2 ~> m2]
  real :: latent_heat_fusion    !< latent heat of fusion [J/kg]
  real :: latent_heat_vapor     !< latent heat of vaporization [J/kg]

  real :: max_p_surf            !< maximum surface pressure that can be exerted by the
                                !! atmosphere and floating sea-ice [R L2 T-2 ~> Pa].
                                !! This is needed because the FMS coupling
                                !! structure does not limit the water that can be
                                !! frozen out of the ocean and the ice-ocean heat
                                !! fluxes are treated explicitly.
  logical :: use_limited_P_SSH  !< If true, return the sea surface height with
                                !! the correction for the atmospheric (and sea-ice)
                                !! pressure limited by max_p_surf instead of the
                                !! full atmospheric pressure.  The default is true.
  logical :: use_CFC            !< enables the MOM_CFC_cap tracer package.
  logical :: use_marbl_tracers  !< enables the MARBL tracer package.
  logical :: enthalpy_cpl       !< Controls if enthalpy terms are provided by the coupler or computed
                                !! internally.
  real :: gust_const            !< constant unresolved background gustiness for ustar [R L Z T-2 ~> Pa]
  logical :: read_gust_2d       !< If true, use a 2-dimensional gustiness supplied
                                !! from an input file.
  real, pointer, dimension(:,:) :: &
    TKE_tidal => NULL(), &      !< turbulent kinetic energy introduced to the
                                !! bottom boundary layer by drag on the tidal flows [R Z3 T-3 ~> W m-2]
    gust => NULL(), &           !< spatially varying unresolved background
                                !! gustiness that contributes to ustar [R L Z T-2 ~> Pa].
                                !! gust is used when read_gust_2d is true.
    ustar_tidal => NULL()       !< tidal contribution to the bottom friction velocity [Z T-1 ~> m s-1]
  real :: cd_tides              !< drag coefficient that applies to the tides (nondimensional)
  real :: utide                 !< constant tidal velocity to use if read_tideamp
                                !! is false [Z T-1 ~> m s-1]
  logical :: read_tideamp       !< If true, spatially varying tidal amplitude read from a file.

  logical :: rigid_sea_ice      !< If true, sea-ice exerts a rigidity that acts
                                !! to damp surface deflections (especially surface
                                !! gravity waves).  The default is false.
  real    :: g_Earth            !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real    :: Kv_sea_ice         !< Viscosity in sea-ice that resists sheared vertical motions [L4 Z-2 T-1 ~> m2 s-1]
  real    :: density_sea_ice    !< Typical density of sea-ice [R ~> kg m-3]. The value is
                                !! only used to convert the ice pressure into
                                !! appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass !< A mass per unit area of sea-ice beyond which
                                !! sea-ice viscosity becomes effective [R Z ~> kg m-2],
                                !! typically of order 1000 kg m-2.
  logical :: allow_flux_adjustments !< If true, use data_override to obtain flux adjustments
  logical :: liquid_runoff_from_data !< If true, use data_override to obtain liquid runoff

  real    :: Flux_const                     !< piston velocity for surface restoring [Z T-1 ~> m s-1]
  logical :: salt_restore_as_sflux          !< If true, SSS restore as salt flux instead of water flux
  logical :: adjust_net_srestore_to_zero    !< adjust srestore to zero (for both salt_flux or vprec)
  logical :: adjust_net_srestore_by_scaling !< adjust srestore w/o moving zero contour
  logical :: adjust_net_fresh_water_to_zero !< adjust net surface fresh-water (w/ restoring) to zero
  logical :: use_net_FW_adjustment_sign_bug !< use the wrong sign when adjusting net FW
  logical :: adjust_net_fresh_water_by_scaling !< adjust net surface fresh-water w/o moving zero contour
  logical :: mask_srestore_under_ice        !< If true, use an ice mask defined by frazil
                                            !! criteria for salinity restoring.
  real    :: ice_salt_concentration         !< salt concentration for sea ice [kg/kg]
  logical :: mask_srestore_marginal_seas    !< if true, then mask SSS restoring in marginal seas
  real    :: max_delta_srestore             !< maximum delta salinity used for restoring [S ~> ppt]
  real    :: max_delta_trestore             !< maximum delta sst used for restoring [C ~> degC]
  real, pointer, dimension(:,:) :: basin_mask => NULL() !< mask for SSS restoring by basin
  logical :: ustar_gustless_bug             !< If true, include a bug in the time-averaging of the
                                            !! gustless wind friction velocity.

  type(diag_ctrl), pointer :: diag                  !< structure to regulate diagnostic output timing
  character(len=200)       :: inputdir              !< directory where NetCDF input files are
  character(len=200)       :: salt_restore_file     !< filename for salt restoring data
  character(len=30)        :: salt_restore_var_name !< name of surface salinity in salt_restore_file
  logical                  :: mask_srestore         !< if true, apply a 2-dimensional mask to the surface
                                                    !< salinity restoring fluxes. The masking file should be
                                                    !< in inputdir/salt_restore_mask.nc and the field should
                                                    !! be named 'mask'
  real, pointer, dimension(:,:) :: srestore_mask => NULL() !< mask for SSS restoring
  character(len=200)       :: temp_restore_file     !< filename for sst restoring data
  character(len=30)        :: temp_restore_var_name !< name of surface temperature in temp_restore_file
  logical                  :: mask_trestore         !< if true, apply a 2-dimensional mask to the surface
                                                    !! temperature restoring fluxes. The masking file should be
                                                    !! in inputdir/temp_restore_mask.nc and the field should
                                                    !! be named 'mask'
  real, pointer, dimension(:,:) :: trestore_mask => NULL() !< mask for SST restoring
  type(external_field) :: srestore_handle
    !< Handle for time-interpolated salt restoration field
  type(external_field) :: trestore_handle
    !< Handle for time-interpolated temperature restoration field
  ! Diagnostics handles
  type(forcing_diags), public :: handles

  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL()

  type(marbl_forcing_CS), pointer :: marbl_forcing_CSp => NULL() !< parameters for getting MARBL forcing
end type surface_forcing_CS

!> Structure corresponding to forcing, but with the elements, units, and conventions
!! that exactly conform to the use for MOM-based coupled models.
type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: lrunoff           =>NULL() !< liquid runoff [km m-2 s-1]
  real, pointer, dimension(:,:) :: frunoff           =>NULL() !< ice runoff [km m-2 s-1]
  real, pointer, dimension(:,:) :: lrunoff_glc       =>NULL() !< liquid glc runoff via rof [km m-2 s-1]
  real, pointer, dimension(:,:) :: frunoff_glc       =>NULL() !< frozen glc runoff via rof [km m-2 s-1]
  real, pointer, dimension(:,:) :: u_flux            =>NULL() !< i-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: v_flux            =>NULL() !< j-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: t_flux            =>NULL() !< sensible heat flux [W m-2]
  real, pointer, dimension(:,:) :: q_flux            =>NULL() !< specific humidity flux [km m-2 s-1]
  real, pointer, dimension(:,:) :: salt_flux         =>NULL() !< salt flux [km m-2 s-1]
  real, pointer, dimension(:,:) :: seaice_melt_heat  =>NULL() !< sea ice and snow melt heat flux [W m-2]
  real, pointer, dimension(:,:) :: seaice_melt       =>NULL() !< water flux due to sea ice and snow melting [km m-2 s-1]
  real, pointer, dimension(:,:) :: lw_flux           =>NULL() !< long wave radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dir   =>NULL() !< direct visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dif   =>NULL() !< diffuse visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dir   =>NULL() !< direct Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dif   =>NULL() !< diffuse Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: lprec             =>NULL() !< mass flux of liquid precip [km m-2 s-1]
  real, pointer, dimension(:,:) :: fprec             =>NULL() !< mass flux of frozen precip [km m-2 s-1]
  real, pointer, dimension(:,:) :: ustar_berg        =>NULL() !< frictional velocity beneath icebergs [m s-1]
  real, pointer, dimension(:,:) :: area_berg         =>NULL() !< area covered by icebergs[m2 m-2]
  real, pointer, dimension(:,:) :: mass_berg         =>NULL() !< mass of icebergs(kg m-2)
  real, pointer, dimension(:,:) :: hrofl             =>NULL() !< heat content from liquid runoff [W m-2]
  real, pointer, dimension(:,:) :: hrofi             =>NULL() !< heat content from frozen runoff [W m-2]
  real, pointer, dimension(:,:) :: hrofl_glc         =>NULL() !< heat content from liquid glc runoff [W m-2]
  real, pointer, dimension(:,:) :: hrofi_glc         =>NULL() !< heat content from frozen glc runoff [W m-2]
  real, pointer, dimension(:,:) :: hrain             =>NULL() !< heat content from liquid precipitation [W m-2]
  real, pointer, dimension(:,:) :: hsnow             =>NULL() !< heat content from frozen precipitation [W m-2]
  real, pointer, dimension(:,:) :: hevap             =>NULL() !< heat content from evaporation [W m-2]
  real, pointer, dimension(:,:) :: hcond             =>NULL() !< heat content from condensation [W m-2]
  real, pointer, dimension(:,:) :: p                 =>NULL() !< pressure of overlying ice and atmosphere
                                                              !< on ocean surface [Pa]
  real, pointer, dimension(:,:) :: ice_fraction      =>NULL() !< fractional ice area [1]
  real, pointer, dimension(:,:) :: u10_sqr           =>NULL() !< wind speed squared at 10m [m2 s-2]
  real, pointer, dimension(:,:) :: nhx_dep           =>NULL() !< Nitrogen deposition [kg m-2 s-1]
  real, pointer, dimension(:,:) :: noy_dep           =>NULL() !< Nitrogen deposition [kg m-2 s-1]
  real, pointer, dimension(:,:) :: atm_co2_prog      =>NULL() !< Prognostic atmospheric co2 concentration [ppm]
  real, pointer, dimension(:,:) :: atm_co2_diag      =>NULL() !< Diagnostic atmospheric co2 concentration [ppm]
  real, pointer, dimension(:,:) :: atm_fine_dust_flux   =>NULL() !< Fine dust flux from atmosphere [kg m-2 s-1]
  real, pointer, dimension(:,:) :: atm_coarse_dust_flux =>NULL() !< Coarse dust flux from atmosphere [kg m-2 s-1]
  real, pointer, dimension(:,:) :: seaice_dust_flux     =>NULL() !< Dust flux from seaice [kg m-2 s-1]
  real, pointer, dimension(:,:) :: atm_bc_flux          =>NULL() !< Black carbon flux from atmosphere [kg m-2 s-1]
  real, pointer, dimension(:,:) :: seaice_bc_flux       =>NULL() !< Black carbon flux from seaice [kg m-2 s-1]
  real, pointer, dimension(:,:) :: afracr               =>NULL() !< Fractional atmosphere coverage wrt ocean [1]
  real, pointer, dimension(:,:) :: swnet_afracr         =>NULL() !< Net shortwave radiation times atmosphere fraction
                                                                 !! positive => into the ocean [W m-2]
  real, pointer, dimension(:,:,:) :: swpen_ifrac_n      =>NULL() !< Net shortwave radiation penetrating into ice and
                                                                 !! ocean times ice fraction for thickness
                                                                 !! positive => into the ocean [W m-2]
  real, pointer, dimension(:,:,:) :: ifrac_n            =>NULL() !< Ice fraction per category [1]
  real, pointer, dimension(:,:) :: mi                =>NULL() !< mass of ice [km m-2]
  real, pointer, dimension(:,:) :: ice_rigidity      =>NULL() !< rigidity of the sea ice, sea-ice and
                                                              !! ice-shelves, expressed as a coefficient
                                                              !! for divergence damping, as determined
                                                              !! outside of the ocean model in [m3 s-1]
  real, pointer, dimension(:,:)   :: lamult          => NULL() !< Langmuir enhancement factor [1]
  real, pointer, dimension(:)     :: stk_wavenumbers => NULL() !< The central wave number of Stokes bands [rad/m]
  real, pointer, dimension(:,:,:) :: ustkb           => NULL() !< Stokes Drift spectrum, zonal [m s-1]
                                                               !! Horizontal  - u points
                                                               !! 3rd dimension - wavenumber
  real, pointer, dimension(:,:,:) :: vstkb           => NULL() !< Stokes Drift spectrum, meridional [m s-1]
                                                               !! Horizontal  - v points
                                                               !! 3rd dimension - wavenumber
  integer :: num_stk_bands            !< Number of Stokes drift bands passed through the coupler
  integer :: xtype                                            !< The type of the exchange - REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type)      :: fluxes                     !< A structure that may contain an array of
                                                              !! named fields used for passive tracer fluxes.
  integer :: wind_stagger = -999                              !< A flag indicating the spatial discretization of
                                                              !! wind stresses.  This flag may be set by the
                                                              !! flux-exchange code, based on what the sea-ice
                                                              !! model is providing.  Otherwise, the value from
                                                              !! the surface_forcing_CS is used.

  ! Forcing when receiving multiple ice categories from CMEPS
  integer                                      :: ice_ncat            !< Number of ice categories coming from coupler
                                                                      !! (1 => not using separate categories)
end type ice_ocean_boundary_type

integer :: id_clock_forcing


  interface
module subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, valid_time, G, US, CS, &
                                 sfc_state, restore_salt, restore_temp)
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
  logical,       optional, intent(in)    :: restore_salt !< If true, salinity is restored to a target value.
  logical,       optional, intent(in)    :: restore_temp !< If true, temperature is restored to a target value.

  ! local variables
end subroutine convert_IOB_to_fluxes
module subroutine convert_IOB_to_forces(IOB, forces, index_bounds, Time, G, US, CS)
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

  ! local variables




end subroutine convert_IOB_to_forces
module subroutine apply_flux_adjustments(G, US, CS, Time, fluxes)
  type(ocean_grid_type),    intent(inout) :: G  !< Ocean grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(forcing),            intent(inout) :: fluxes !< Surface fluxes structure

  ! Local variables


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
module subroutine surface_forcing_init(Time, G, US, param_file, diag, CS, restore_salt, restore_temp, use_waves)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure that is used to regulate
                                                  !! diagnostic output
  type(surface_forcing_CS), pointer       :: CS   !< A pointer that is set to point to the control
                                                  !! structure for this module
  logical, optional,        intent(in)    :: restore_salt !< If present and true surface salinity
                                                  !! restoring will be applied in this model.
  logical, optional,        intent(in)    :: restore_temp !< If present and true surface temperature
                                                  !! restoring will be applied in this model.
  logical, optional,        intent(in)    :: use_waves !< If present and true, use waves and activate
                                                  !! the corresponding wave forcing diagnostics

  ! Local variables
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
  end interface

end module MOM_surface_forcing_nuopc
