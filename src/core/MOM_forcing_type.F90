! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module implements boundary forcing for MOM6.
module MOM_forcing_type

use MOM_array_transform, only : rotate_array, rotate_vector, rotate_array_pair
use MOM_coupler_types, only : coupler_2d_bc_type, coupler_type_destructor
use MOM_coupler_types, only : coupler_type_increment_data, coupler_type_initialized
use MOM_coupler_types, only : coupler_type_copy_data, coupler_type_spawn
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_debugging,     only : hchksum, uvchksum
use MOM_diag_mediator, only : post_data, register_diag_field, register_scalar_field
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_alloc, query_averaging_enabled
use MOM_diag_mediator, only : enable_averages, disable_averaging
use MOM_EOS,           only : calculate_density_derivs, calculate_specific_vol_derivs, EOS_domain
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_opacity,       only : sumSWoverBands, optics_type, extract_optics_slice, optics_nbands
use MOM_spatial_means, only : global_area_integral, global_area_mean
use MOM_spatial_means, only : global_area_mean_u, global_area_mean_v
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, extractFluxes2d, optics_type
public MOM_forcing_chksum, MOM_mech_forcing_chksum
public calculateBuoyancyFlux1d, calculateBuoyancyFlux2d, find_ustar
public forcing_accumulate, fluxes_accumulate
public forcing_SinglePointPrint, mech_forcing_diags, forcing_diagnostics
public register_forcing_type_diags, allocate_forcing_type, deallocate_forcing_type
public copy_common_forcing_fields, allocate_mech_forcing, deallocate_mech_forcing
public set_derived_forcing_fields, copy_back_forcing_fields
public set_net_mass_forcing, get_net_mass_forcing
public rotate_forcing, rotate_mech_forcing
public homogenize_forcing, homogenize_mech_forcing

!> Allocate the fields of a (flux) forcing type, based on either a set of input
!! flags for each group of fields, or a pre-allocated reference forcing.
interface allocate_forcing_type
  module procedure allocate_forcing_by_group
  module procedure allocate_forcing_by_ref
end interface allocate_forcing_type

!> Allocate the fields of a mechanical forcing type, based on either a set of
!! input flags for each group of fields, or a pre-allocated reference forcing.
interface allocate_mech_forcing
  module procedure allocate_mech_forcing_by_group
  module procedure allocate_mech_forcing_from_ref
end interface allocate_mech_forcing

!> Allocate arrays if optional flag is present and true (works for 2D and 3D)
interface myAlloc
  module procedure myAlloc_2d
  module procedure myAlloc_3d
end interface myAlloc

!> Determine the friction velocity from a forcing type or a mechanical forcing type.
interface find_ustar
  module procedure find_ustar_fluxes
  module procedure find_ustar_mech_forcing
end interface find_ustar

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Structure that contains pointers to the boundary forcing used to drive the
!! liquid ocean simulated by MOM.
!!
!! Data in this type is allocated in the module MOM_surface_forcing.F90, of which there
!! are three: solo, coupled, and ice-shelf. Alternatively, they are allocated in
!! MESO_surface_forcing.F90, which is a special case of solo_driver/MOM_surface_forcing.F90.
type, public :: forcing

  ! surface stress components and turbulent velocity scale
  real, pointer, dimension(:,:) :: &
    omega_w2x     => NULL(), & !< the counter-clockwise angle of the wind stress with respect
    ustar         => NULL(), & !< surface friction velocity scale [Z T-1 ~> m s-1].
    tau_mag       => NULL(), & !< Magnitude of the wind stress averaged over tracer cells,
                               !! including any contributions from sub-gridscale variability
                               !! or gustiness, rescaled to units that are more convenient for
                               !! calculating turbulent fluxes and friction velocities [R Z2 T-2 ~> Pa]
    ustar_gustless => NULL(), & !< surface friction velocity scale without any
                               !! any augmentation for gustiness [Z T-1 ~> m s-1].
    tau_mag_gustless => NULL() !< Magnitude of the wind stress averaged over tracer cells,
                               !! without any augmentation for sub-gridscale variability
                               !! or gustiness, rescaled to units that are more convenient for
                               !! calculating turbulent fluxes and friction velocities [R Z2 T-2 ~> Pa]

  ! surface buoyancy force, used when temperature is not a state variable
  real, pointer, dimension(:,:) :: &
    buoy          => NULL()  !< buoyancy flux [L2 T-3 ~> m2 s-3]

  ! radiative heat fluxes into the ocean [Q R Z T-1 ~> W m-2]
  real, pointer, dimension(:,:) :: &
    sw         => NULL(), & !< shortwave [Q R Z T-1 ~> W m-2]
    sw_vis_dir => NULL(), & !< visible, direct shortwave [Q R Z T-1 ~> W m-2]
    sw_vis_dif => NULL(), & !< visible, diffuse shortwave [Q R Z T-1 ~> W m-2]
    sw_nir_dir => NULL(), & !< near-IR, direct shortwave [Q R Z T-1 ~> W m-2]
    sw_nir_dif => NULL(), & !< near-IR, diffuse shortwave [Q R Z T-1 ~> W m-2]
    lw         => NULL()    !< longwave [Q R Z T-1 ~> W m-2] (typically negative)

  ! turbulent heat fluxes into the ocean [Q R Z T-1 ~> W m-2]
  real, pointer, dimension(:,:) :: &
    latent           => NULL(), & !< latent [Q R Z T-1 ~> W m-2] (typically < 0)
    sens             => NULL(), & !< sensible [Q R Z T-1 ~> W m-2] (typically negative)
    seaice_melt_heat => NULL(), & !< sea ice and snow melt or formation [Q R Z T-1 ~> W m-2] (typically negative)
    heat_added       => NULL()    !< additional heat flux from SST restoring or flux adjustments [Q R Z T-1 ~> W m-2]

  ! components of latent heat fluxes used for diagnostic purposes
  real, pointer, dimension(:,:) :: &
    latent_evap_diag        => NULL(), & !< latent [Q R Z T-1 ~> W m-2] from evaporating liquid water (typically < 0)
    latent_fprec_diag       => NULL(), & !< latent [Q R Z T-1 ~> W m-2] from melting fprec  (typically < 0)
    latent_frunoff_diag     => NULL(), & !< latent [Q R Z T-1 ~> W m-2] from melting frunoff (calving) (typically < 0)
    latent_frunoff_glc_diag => NULL()    !< latent [Q R Z T-1 ~> W m-2] from melting glacier frunoff (typically < 0)

  ! water mass fluxes into the ocean [R Z T-1 ~> kg m-2 s-1]; these fluxes impact the ocean mass
  real, pointer, dimension(:,:) :: &
    evap          => NULL(), & !< (-1)*fresh water flux evaporated out of the ocean [R Z T-1 ~> kg m-2 s-1]
    lprec         => NULL(), & !< precipitating liquid water into the ocean [R Z T-1 ~> kg m-2 s-1]
    fprec         => NULL(), & !< precipitating frozen water into the ocean [R Z T-1 ~> kg m-2 s-1]
    vprec         => NULL(), & !< virtual liquid precip associated w/ SSS restoring [R Z T-1 ~> kg m-2 s-1]
    lrunoff       => NULL(), & !< liquid river runoff entering ocean [R Z T-1 ~> kg m-2 s-1]
    frunoff       => NULL(), & !< frozen river runoff (calving) entering ocean [R Z T-1 ~> kg m-2 s-1]
    lrunoff_glc   => NULL(), & !< liquid river glacier runoff entering ocean [R Z T-1 ~> kg m-2 s-1]
    frunoff_glc   => NULL(), & !< frozen river glacier runoff entering ocean [R Z T-1 ~> kg m-2 s-1]
    seaice_melt   => NULL()    !< snow/seaice melt (positive) or formation (negative) [R Z T-1 ~> kg m-2 s-1]

  ! carbon content associated with water crossing ocean surface
  real, pointer, dimension(:,:) :: &
    carbon_content_lrunoff     => NULL() !< carbon content associated with liquid runoff [R Z T-1 ~> kg m-2 s-1]

  ! Integrated water mass fluxes into the ocean, used for passive tracer sources [H ~> m or kg m-2]
  real, pointer, dimension(:,:) :: &
    netMassIn     => NULL(), & !< Sum of water mass fluxes into the ocean integrated over a
                               !! forcing timestep [H ~> m or kg m-2]
    netMassOut    => NULL()    !< Net water mass flux out of the ocean integrated over a forcing timestep,
                               !! with negative values for water leaving the ocean [H ~> m or kg m-2]

  ! heat associated with water crossing ocean surface
  real, pointer, dimension(:,:) :: &
    heat_content_cond        => NULL(), & !< heat content associated with condensating water [Q R Z T-1 ~> W m-2]
    heat_content_evap        => NULL(), & !< heat content associated with evaporating water  [Q R Z T-1 ~> W m-2]
    heat_content_lprec       => NULL(), & !< heat content associated with liquid >0 precip   [Q R Z T-1 ~> W m-2]
    heat_content_fprec       => NULL(), & !< heat content associated with frozen precip      [Q R Z T-1 ~> W m-2]
    heat_content_vprec       => NULL(), & !< heat content associated with virtual >0 precip  [Q R Z T-1 ~> W m-2]
    heat_content_lrunoff     => NULL(), & !< heat content associated with liquid runoff      [Q R Z T-1 ~> W m-2]
    heat_content_frunoff     => NULL(), & !< heat content associated with frozen runoff      [Q R Z T-1 ~> W m-2]
    heat_content_lrunoff_glc => NULL(), & !< heat content associated with liquid runoff      [Q R Z T-1 ~> W m-2]
    heat_content_frunoff_glc => NULL(), & !< heat content associated with frozen runoff      [Q R Z T-1 ~> W m-2]
    heat_content_massout     => NULL(), & !< heat content associated with mass leaving ocean [Q R Z T-1 ~> W m-2]
    heat_content_massin      => NULL()    !< heat content associated with mass entering ocean [Q R Z T-1 ~> W m-2]

  ! salt mass flux (contributes to ocean mass only if non-Bouss )
  real, pointer, dimension(:,:) :: &
    salt_flux       => NULL(), & !< net salt flux into the ocean [R Z T-1 ~> kgSalt m-2 s-1]
    salt_flux_in    => NULL(), & !< salt flux provided to the ocean from coupler [R Z T-1 ~> kgSalt m-2 s-1]
    salt_flux_added => NULL(), & !< additional salt flux from restoring or flux adjustment before adjustment
                                 !! to net zero [R Z T-1 ~> kgSalt m-2 s-1]
    salt_left_behind => NULL()   !< salt left in ocean at the surface from brine rejection
                                 !! [R Z T-1 ~> kgSalt m-2 s-1]

  ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
  real, pointer, dimension(:,:) :: p_surf_full => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa].
                !! if there is sea-ice, then p_surf_flux is at ice-ocean interface
  real, pointer, dimension(:,:) :: p_surf => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] as used to drive the ocean model.
                !! If p_surf is limited, p_surf may be smaller than p_surf_full, otherwise they are the same.
  real, pointer, dimension(:,:) :: p_surf_SSH => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] that is used in corrections to the sea surface
                !! height field that is passed back to the calling routines.
                !! p_surf_SSH may point to p_surf or to p_surf_full.
  logical :: accumulate_p_surf = .false. !< If true, the surface pressure due to the atmosphere
                                 !! and various types of ice needs to be accumulated, and the
                                 !! surface pressure explicitly reset to zero at the driver level
                                 !! when appropriate.

  ! tide related inputs
  real, pointer, dimension(:,:) :: &
    BBL_tidal_dis => NULL(), & !< Tidal energy dissipation in the bottom boundary layer that can act
                               !! as a source of energy for bottom boundary layer mixing [R Z L2 T-3 ~> W m-2]
    ustar_tidal   => NULL()    !< tidal contribution to bottom ustar [Z T-1 ~> m s-1]

  ! iceberg related inputs
  real, pointer, dimension(:,:) :: &
    ustar_berg => NULL(), &   !< iceberg contribution to top ustar [Z T-1 ~> m s-1].
    area_berg  => NULL(), &   !< fractional area of ocean surface covered by icebergs [nondim]
    mass_berg  => NULL()      !< mass of icebergs [R Z ~> kg m-2]

  ! land ice-shelf related inputs
  real, pointer, dimension(:,:) :: ustar_shelf => NULL()  !< Friction velocity under ice-shelves [Z T-1 ~> m s-1].
                                 !! as computed by the ocean at the previous time step.
  real, pointer, dimension(:,:) :: frac_shelf_h => NULL() !< Fractional ice shelf coverage of
                                 !! h-cells, from 0 to 1 [nondim]. This is only
                                 !! associated if ice shelves are enabled, and are
                                 !! exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: iceshelf_melt => NULL() !< Ice shelf melt rate (positive)
                                 !! or freezing (negative) [R Z T-1 ~> kg m-2 s-1]
  real, pointer, dimension(:,:) :: shelf_sfc_mass_flux => NULL() !< Ice shelf surface mass flux
                                 !! deposition from the atmosphere. [R Z T-1 ~> kg m-2 s-1]

  ! Scalars set by surface forcing modules
  real :: vPrecGlobalAdj = 0.     !< adjustment to restoring vprec to zero out global net [R Z T-1 ~> kg m-2 s-1]
  real :: saltFluxGlobalAdj = 0.  !< adjustment to restoring salt flux to zero out global
                                  !! net [R Z T-1 ~> kgSalt m-2 s-1]
  real :: netFWGlobalAdj = 0.     !< adjustment to net fresh water to zero out global net [R Z T-1 ~> kg m-2 s-1]
  real :: vPrecGlobalScl = 0.     !< scaling of restoring vprec to zero out global net ( -1..1 ) [nondim]
  real :: saltFluxGlobalScl = 0.  !< scaling of restoring salt flux to zero out global net ( -1..1 ) [nondim]
  real :: netFWGlobalScl = 0.     !< scaling of net fresh water to zero out global net ( -1..1 ) [nondim]

  logical :: fluxes_used = .true. !< If true, all of the heat, salt, and mass
                                  !! fluxes have been applied to the ocean.
  real :: dt_buoy_accum = -1.0    !< The amount of time over which the buoyancy fluxes
                                  !! should be applied [T ~> s].  If negative, this forcing
                                  !! type variable has not yet been initialized.
  logical :: gustless_accum_bug = .true. !< If true, use an incorrect expression in the time
                                  !! average of the gustless wind stress.
  real :: C_p                   !< heat capacity of seawater [Q C-1 ~> J kg-1 degC-1].
                                !! C_p is is the same value as in thermovar_ptrs_type.

  ! arrays needed in the some tracer modules, e.g., MOM_CFC_cap
  real, pointer, dimension(:,:) :: &
    ice_fraction  => NULL(), &  !< fraction of sea ice coverage at h-cells, from 0 to 1 [nondim].
    u10_sqr       => NULL()     !< wind magnitude at 10 m squared [L2 T-2 ~> m2 s-2]

  ! Forcing fields required for MARBL
  real, pointer, dimension(:,:) :: &
    noy_dep => NULL(),               & !< NOy Deposition [conc Z T-1 ~> conc m s-1]
    nhx_dep => NULL(),               & !< NHx Deposition [conc Z T-1 ~> conc m s-1]
    atm_co2 => NULL(),               & !< Atmospheric CO2 Concentration [ppm]
    atm_alt_co2 => NULL(),           & !< Alternate atmospheric CO2 Concentration [ppm]
    dust_flux => NULL(),             & !< Flux of dust into the ocean [R Z T-1 ~> kgN m-2 s-1]
    iron_flux => NULL()                !< Flux of dust into the ocean [conc Z T-1 ~> conc m s-1]

  real, pointer, dimension(:,:,:) :: &
    fracr_cat   => NULL(),           & !< per-category ice fraction [nondim]
    qsw_cat     => NULL()              !< per-category shortwave [Q R Z T-1 ~> W m-2]

  real, pointer, dimension(:,:) :: &
    lamult => NULL()            !< Langmuir enhancement factor [nondim]

  ! passive tracer surface fluxes
  type(coupler_2d_bc_type) :: tr_fluxes !< This structure contains arrays of
     !! of named fields used for passive tracer fluxes.
     !! All arrays in tr_fluxes use the coupler indexing, which has no halos.
     !! This is not a convenient convention, but imposed on MOM6 by the coupler.

  ! For internal error tracking
  integer :: num_msg = 0 !< Number of messages issued about excessive SW penetration
  integer :: max_msg = 2 !< Maximum number of messages to issue about excessive SW penetration

end type forcing

!> Structure that contains pointers to the mechanical forcing at the surface
!! used to drive the liquid ocean simulated by MOM.
!! Data in this type is allocated in the module MOM_surface_forcing.F90,
!! of which there are three versions:  solo, coupled, and ice-shelf.
type, public :: mech_forcing
  ! surface stress components and turbulent velocity scale
  real, pointer, dimension(:,:) :: &
    taux  => NULL(), & !< zonal wind stress [R L Z T-2 ~> Pa]
    tauy  => NULL(), & !< meridional wind stress [R L Z T-2 ~> Pa]
    tau_mag => NULL(), & !< Magnitude of the wind stress averaged over tracer cells, including any
                       !! contributions from sub-gridscale variability or gustiness [R L Z T-2 ~> Pa]
    ustar => NULL(), & !< surface friction velocity scale [Z T-1 ~> m s-1].
    net_mass_src => NULL(), & !< The net mass source to the ocean [R Z T-1 ~> kg m-2 s-1]
    omega_w2x    => NULL()    !< the counter-clockwise angle of the wind stress with respect
                              !! to the horizontal abscissa (x-coordinate) at tracer points [rad].

  ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
  real, pointer, dimension(:,:) :: p_surf_full => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa].
                !! if there is sea-ice, then p_surf_flux is at ice-ocean interface
  real, pointer, dimension(:,:) :: p_surf => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] as used to drive the ocean model.
                !! If p_surf is limited, p_surf may be smaller than p_surf_full, otherwise they are the same.
  real, pointer, dimension(:,:) :: p_surf_SSH => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] that is used in corrections
                !! to the sea surface height field that is passed back to the calling routines.
                !! p_surf_SSH may point to p_surf or to p_surf_full.

  ! iceberg related inputs
  real, pointer, dimension(:,:) :: &
    area_berg  => NULL(), &    !< fractional area of ocean surface covered by icebergs [nondim]
    mass_berg  => NULL()       !< mass of icebergs per unit ocean area [R Z ~> kg m-2]

  ! land ice-shelf related inputs
  real, pointer, dimension(:,:) :: frac_shelf_u  => NULL() !< Fractional ice shelf coverage of u-cells,
                !! nondimensional from 0 to 1 [nondim]. This is only associated if ice shelves are enabled,
                !! and is exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: frac_shelf_v  => NULL() !< Fractional ice shelf coverage of v-cells,
                !! nondimensional from 0 to 1 [nondim]. This is only associated if ice shelves are enabled,
                !! and is exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: &
    rigidity_ice_u => NULL(), & !< Depth-integrated lateral viscosity of ice shelves or sea ice at
                                !! u-points [L4 Z-1 T-1 ~> m3 s-1]
    rigidity_ice_v => NULL()    !< Depth-integrated lateral viscosity of ice shelves or sea ice at
                                !! v-points [L4 Z-1 T-1 ~> m3 s-1]
  real :: dt_force_accum = -1.0 !< The amount of time over which the mechanical forcing fluxes
                                !! have been averaged [T ~> s].
  logical :: net_mass_src_set = .false. !< If true, an estimate of net_mass_src has been provided.
  logical :: accumulate_p_surf = .false. !< If true, the surface pressure due to the atmosphere
                                !! and various types of ice needs to be accumulated, and the
                                !! surface pressure explicitly reset to zero at the driver level
                                !! when appropriate.
  logical :: accumulate_rigidity = .false. !< If true, the rigidity due to various types of
                                !! ice needs to be accumulated, and the rigidity explicitly
                                !! reset to zero at the driver level when appropriate.
  real, pointer, dimension(:) :: &
    stk_wavenumbers => NULL()   !< The central wave number of Stokes bands [rad Z-1 ~> rad m-1]
  real, pointer, dimension(:,:,:) :: &
    ustkb => NULL(), &          !< Stokes Drift spectrum, zonal [L T-1 ~> m s-1]
                                !! Horizontal - u points
                                !! 3rd dimension - wavenumber
    vstkb => NULL()             !< Stokes Drift spectrum, meridional [L T-1 ~> m s-1]
                                !! Horizontal - v points
                                !! 3rd dimension - wavenumber

  logical :: initialized = .false. !< This indicates whether the appropriate arrays have been initialized.
end type mech_forcing

!> Structure that defines the id handles for the forcing type
type, public :: forcing_diags ; private

  !>@{ Forcing diagnostic handles
  ! mass flux diagnostic handles
  integer :: id_prcme        = -1, id_evap        = -1
  integer :: id_precip       = -1, id_vprec       = -1
  integer :: id_lprec        = -1, id_fprec       = -1
  integer :: id_lrunoff      = -1, id_frunoff     = -1
  integer :: id_lrunoff_glc  = -1, id_frunoff_glc = -1
  integer :: id_net_massout  = -1, id_net_massin  = -1
  integer :: id_massout_flux = -1, id_massin_flux = -1
  integer :: id_seaice_melt  = -1

  ! global area integrated mass flux diagnostic handles
  integer :: id_total_prcme        = -1, id_total_evap        = -1
  integer :: id_total_precip       = -1, id_total_vprec       = -1
  integer :: id_total_lprec        = -1, id_total_fprec       = -1
  integer :: id_total_lrunoff      = -1, id_total_frunoff     = -1
  integer :: id_total_lrunoff_glc  = -1, id_total_frunoff_glc = -1
  integer :: id_total_net_massout  = -1, id_total_net_massin  = -1
  integer :: id_total_seaice_melt  = -1

  ! global area averaged mass flux diagnostic handles
  integer :: id_prcme_ga  = -1, id_evap_ga = -1
  integer :: id_lprec_ga  = -1, id_fprec_ga= -1
  integer :: id_precip_ga = -1, id_vprec_ga= -1

  ! heat flux diagnostic handles
  integer :: id_net_heat_coupler        = -1, id_net_heat_surface        = -1
  integer :: id_sens                    = -1, id_LwLatSens               = -1
  integer :: id_sw                      = -1, id_lw                      = -1
  integer :: id_sw_vis                  = -1, id_sw_nir                  = -1
  integer :: id_lat_evap                = -1, id_lat_frunoff             = -1
  integer :: id_lat_frunoff_glc         = -1
  integer :: id_lat                     = -1, id_lat_fprec               = -1
  integer :: id_heat_content_lrunoff    = -1, id_heat_content_frunoff    = -1
  integer :: id_heat_content_lrunoff_glc= -1, id_heat_content_frunoff_glc= -1
  integer :: id_heat_content_lprec      = -1, id_heat_content_fprec      = -1
  integer :: id_heat_content_cond       = -1, id_heat_content_surfwater  = -1
  integer :: id_heat_content_evap       = -1
  integer :: id_heat_content_vprec      = -1, id_heat_content_massout    = -1
  integer :: id_heat_added              = -1, id_heat_content_massin     = -1
  integer :: id_hfrainds                = -1, id_hfrunoffds              = -1
  integer :: id_seaice_melt_heat        = -1
  integer :: id_carbon_content_lrunoff  = -1

  ! global area integrated heat flux diagnostic handles
  integer :: id_total_net_heat_coupler        = -1, id_total_net_heat_surface        = -1
  integer :: id_total_sens                    = -1, id_total_LwLatSens               = -1
  integer :: id_total_sw                      = -1, id_total_lw                      = -1
  integer :: id_total_lat_evap                = -1, id_total_lat_frunoff             = -1
  integer :: id_total_lat_frunoff_glc         = -1
  integer :: id_total_lat                     = -1, id_total_lat_fprec               = -1
  integer :: id_total_heat_content_lrunoff    = -1, id_total_heat_content_frunoff    = -1
  integer :: id_total_heat_content_lrunoff_glc= -1, id_total_heat_content_frunoff_glc=-1
  integer :: id_total_heat_content_lprec      = -1, id_total_heat_content_fprec      = -1
  integer :: id_total_heat_content_cond       = -1, id_total_heat_content_surfwater  = -1
  integer :: id_total_heat_content_evap       = -1
  integer :: id_total_heat_content_vprec      = -1, id_total_heat_content_massout    = -1
  integer :: id_total_heat_added              = -1, id_total_heat_content_massin     = -1
  integer :: id_total_seaice_melt_heat        = -1

  ! global area averaged heat flux diagnostic handles
  integer :: id_net_heat_coupler_ga = -1, id_net_heat_surface_ga = -1
  integer :: id_sens_ga             = -1, id_LwLatSens_ga        = -1
  integer :: id_sw_ga               = -1, id_lw_ga               = -1
  integer :: id_lat_ga              = -1

  ! salt flux diagnostic handles
  integer :: id_saltflux          = -1
  integer :: id_saltFluxIn        = -1
  integer :: id_saltFluxAdded     = -1
  integer :: id_saltFluxBehind    = -1

  integer :: id_total_saltflux        = -1
  integer :: id_total_saltFluxIn      = -1
  integer :: id_total_saltFluxAdded   = -1

  integer :: id_vPrecGlobalAdj    = -1
  integer :: id_vPrecGlobalScl    = -1
  integer :: id_saltFluxGlobalAdj = -1
  integer :: id_saltFluxGlobalScl = -1
  integer :: id_netFWGlobalAdj    = -1
  integer :: id_netFWGlobalScl    = -1

  ! momentum flux and forcing diagnostic handles
  integer :: id_taux  = -1
  integer :: id_tauy  = -1
  integer :: id_ustar = -1
  integer :: id_omega_w2x = -1
  integer :: id_tau_mag = -1
  integer :: id_psurf     = -1
  integer :: id_TKE_tidal = -1
  integer :: id_buoy      = -1

  ! tracer surface flux related diagnostics handles
  integer :: id_ice_fraction = -1
  integer :: id_u10_sqr      = -1

  ! iceberg diagnostic handles
  integer :: id_ustar_berg = -1
  integer :: id_area_berg = -1
  integer :: id_mass_berg = -1

  ! Iceberg + Ice shelf diagnostic handles
  integer :: id_ustar_ice_cover = -1
  integer :: id_frac_ice_cover = -1
  ! wave forcing diagnostics handles.
  integer :: id_lamult = -1
  !>@}

  integer :: id_clock_forcing = -1 !< CPU clock id

end type forcing_diags


  interface
module subroutine extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt, &
                  FluxRescaleDepth, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, net_salt, pen_SW_bnd, tv, &
                  aggregate_FW, nonpenSW, netmassInOut_rate, net_Heat_Rate, &
                  net_salt_rate, pen_sw_bnd_Rate)

  type(ocean_grid_type),    intent(in)    :: G              !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),            intent(inout) :: fluxes         !< structure containing pointers to possible
                                                            !! forcing fields. NULL unused fields.
  type(optics_type),        pointer       :: optics         !< pointer to optics
  integer,                  intent(in)    :: nsw            !< number of bands of penetrating SW
  integer,                  intent(in)    :: j              !< j-index to work on
  real,                     intent(in)    :: dt             !< The time step for these fluxes [T ~> s]
  real,                     intent(in)    :: FluxRescaleDepth !< min ocean depth before fluxes
                                                            !! are scaled away [H ~> m or kg m-2]
  logical,                  intent(in)    :: useRiverHeatContent   !< logical for river heat content
  logical,                  intent(in)    :: useCalvingHeatContent !< logical for calving heat content
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: T              !< layer temperatures [C ~> degC]
  real, dimension(SZI_(G)), intent(out)   :: netMassInOut   !< net mass flux (non-Bouss) or volume flux
                                                            !! (if Bouss) of water in/out of ocean over
                                                            !! a time step [H ~> m or kg m-2]
  real, dimension(SZI_(G)), intent(out)   :: netMassOut     !< net mass flux (non-Bouss) or volume flux
                                                            !! (if Bouss) of water leaving ocean surface
                                                            !! over a time step [H ~> m or kg m-2].
                                                            !! netMassOut < 0 means mass leaves ocean.
  real, dimension(SZI_(G)), intent(out)   :: net_heat       !< net heat at the surface accumulated over a
                                                            !! time step for coupler + restoring.
                                                            !! Exclude two terms from net_heat:
                                                            !! (1) downwelling (penetrative) SW,
                                                            !! (2) evaporation heat content,
                                                            !! (since do not yet know evap temperature).
                                                            !! [C H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: net_salt       !< surface salt flux into the ocean
                                                            !! accumulated over a time step
                                                            !! [S H ~> ppt m or ppt kg m-2].
  real, dimension(max(1,nsw),G%isd:G%ied), intent(out) :: pen_SW_bnd !< penetrating SW flux, split into bands.
                                                            !! [C H ~> degC m or degC kg m-2]
                                                            !! and array size nsw x SZI_(G), where
                                                            !! nsw=number of SW bands in pen_SW_bnd.
                                                            !! This heat flux is not part of net_heat.
  type(thermo_var_ptrs),    intent(inout) :: tv             !< structure containing pointers to available
                                                            !! thermodynamic fields. Used to keep
                                                            !! track of the heat flux associated with net
                                                            !! mass fluxes into the ocean.
  logical,                  intent(in)    :: aggregate_FW   !< For determining how to aggregate forcing.
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: nonpenSW       !< Non-penetrating SW used in net_heat
                                                            !! [C H ~> degC m or degC kg m-2].
                                                            !! Summed over SW bands when diagnosing nonpenSW.
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: net_Heat_rate  !< Rate of net surface heating
                                                            !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: net_salt_rate  !< Surface salt flux into the ocean
                                                            !! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1].
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: netmassInOut_rate !< Rate of net mass flux into the ocean
                                                            !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real, dimension(max(1,nsw),G%isd:G%ied), &
                  optional, intent(out)   :: pen_sw_bnd_rate !< Rate of penetrative shortwave heating
                                                             !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].

  ! local
                              ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
                              ! or 0 for no limiting [H-1 ~> m-1 or m2 kg-1]
                              ! [C H R-1 Z-1 Q-1 ~> degC m3 J-1 or kg degC J-1]


  !BGR-Jul 5,2017{
  ! Initializes/sets logicals if 'rates' are requested
  ! These factors are required for legacy reasons
  !  and therefore computed only when optional outputs are requested
end subroutine extractFluxes1d
module subroutine extractFluxes2d(G, GV, US, fluxes, optics, nsw, dt, FluxRescaleDepth, &
                           useRiverHeatContent, useCalvingHeatContent, h, T, &
                           netMassInOut, netMassOut, net_heat, Net_salt, Pen_SW_bnd, tv, &
                           aggregate_FW)

  type(ocean_grid_type),            intent(in)    :: G              !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),                    intent(inout) :: fluxes         !< structure containing pointers to forcing.
  type(optics_type),                pointer       :: optics         !< pointer to optics
  integer,                          intent(in)    :: nsw            !< number of bands of penetrating SW
  real,                             intent(in)    :: dt             !< The time step for these fluxes [T ~> s]
  real,                             intent(in)    :: FluxRescaleDepth !< min ocean depth before fluxes
                                                                    !! are scaled away [H ~> m or kg m-2]
  logical,                          intent(in)    :: useRiverHeatContent   !< logical for river heat content
  logical,                          intent(in)    :: useCalvingHeatContent !< logical for calving heat content
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: T              !< layer temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: netMassInOut   !< net mass flux (non-Bouss) or volume flux
                                                                    !! (if Bouss) of water in/out of ocean over
                                                                    !! a time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: netMassOut     !< net mass flux (non-Bouss) or volume flux
                                                                    !! (if Bouss) of water leaving ocean surface
                                                                    !! over a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: net_heat       !< net heat at the surface accumulated over a
                                                                    !! time step associated with coupler + restore.
                                                                    !! Exclude two terms from net_heat:
                                                                    !! (1) downwelling (penetrative) SW,
                                                                    !! (2) evaporation heat content,
                                                                    !! (since do not yet know temperature of evap).
                                                                    !! [C H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: net_salt       !< surface salt flux into the ocean accumulated
                                                                    !! over a time step [S H ~> ppt m or ppt kg m-2]
  real, dimension(max(1,nsw),G%isd:G%ied,G%jsd:G%jed), intent(out) :: pen_SW_bnd !< penetrating SW flux, by frequency
                                                                    !! band [C H ~> degC m or degC kg m-2] with array
                                                                    !! size nsw x SZI_(G), where nsw=number of SW bands
                                                                    !! in pen_SW_bnd. This heat flux is not in net_heat.
  type(thermo_var_ptrs),            intent(inout) :: tv             !< structure containing pointers to available
                                                                    !! thermodynamic fields. Here it is used to keep
                                                                    !! track of the heat flux associated with net
                                                                    !! mass fluxes into the ocean.
  logical,                          intent(in)    :: aggregate_FW   !< For determining how to aggregate the forcing.

  !$OMP parallel do default(shared)
end subroutine extractFluxes2d
module subroutine calculateBuoyancyFlux1d(G, GV, US, fluxes, optics, nsw, h, Temp, Salt, tv, j, &
                                   buoyancyFlux, netHeatMinusSW, netSalt)
  type(ocean_grid_type),                    intent(in)    :: G              !< ocean grid
  type(verticalGrid_type),                  intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),                    intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),                            intent(inout) :: fluxes         !< surface fluxes
  type(optics_type),                        pointer       :: optics         !< penetrating SW optics
  integer,                                  intent(in)    :: nsw            !< The number of frequency bands of
                                                                            !! penetrating shortwave radiation
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: Temp           !< prognostic temp [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: Salt           !< salinity [S ~> ppt]
  type(thermo_var_ptrs),                    intent(inout) :: tv             !< thermodynamics type
  integer,                                  intent(in)    :: j              !< j-row to work on
  real, dimension(SZI_(G),SZK_(GV)+1),      intent(out)   :: buoyancyFlux   !< buoyancy fluxes [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G)),                 intent(out)   :: netHeatMinusSW !< Surface heat flux excluding shortwave
                                                                          !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(SZI_(G)),                 intent(out)   :: netSalt        !< surface salt flux
                                                                            !! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]

  ! local variables
                                                      ! [H T-1 ~> m s-1 or kg m-2 s-1]
                                                      ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
                                                      ! to temperature [R-1 C-1 ~> m3 kg-1 degC-1]
                                                      ! to salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
                                                      ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]

                    ! unit conversion factor [L2 H-1 R-1 T-2 ~> m4 kg-1 s-2 or m7 kg-2 s-2]
                    ! thickness units to mass per units area [R L2 H-1 T-2 ~> kg m-2 s-2 or m s-2]
                            ! it is necessary to eliminate fluxes [H ~> m or kg m-2]

  !  smg: what do we do when have heat fluxes from calving and river?
end subroutine calculateBuoyancyFlux1d
module subroutine calculateBuoyancyFlux2d(G, GV, US, fluxes, optics, h, Temp, Salt, tv, &
                                   buoyancyFlux, netHeatMinusSW, netSalt)
  type(ocean_grid_type),                      intent(in)    :: G      !< ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  type(forcing),                              intent(inout) :: fluxes !< surface fluxes
  type(optics_type),                          pointer       :: optics !< SW ocean optics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Temp   !< temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Salt   !< salinity [S ~> ppt]
  type(thermo_var_ptrs),                      intent(inout) :: tv     !< thermodynamics type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: buoyancyFlux  !< buoyancy fluxes [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: netHeatMinusSW !< surface heat flux excluding shortwave
                                                                      !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(SZI_(G),SZJ_(G)),           intent(inout) :: netSalt !< Net surface salt flux
                                                                      !! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]

  ! local variables

  !$OMP parallel do default(shared)
end subroutine calculateBuoyancyFlux2d
module subroutine find_ustar_fluxes(fluxes, tv, U_star, G, GV, US, halo, H_T_units)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(forcing),           intent(in)  :: fluxes !< Surface fluxes container
  type(thermo_var_ptrs),   intent(in)  :: tv   !< Structure containing pointers to any
                                               !! available thermodynamic fields.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(out) :: U_star !< The surface friction velocity [Z T-1 ~> m s-1] or
                                               !! [H T-1 ~> m s-1 or kg m-2 s-1], depending on H_T_units.
  integer,       optional, intent(in)  :: halo !< The extra halo size to fill in, 0 by default
  logical,       optional, intent(in)  :: H_T_units !< If present and true, return U_star in units
                                               !! of [H T-1 ~> m s-1 or kg m-2 s-1]

  ! Local variables
                       ! or in some semi-Boussinesq cases the reference
                       ! density [H2 Z-2 R-1 ~> m3 kg-1 or kg m-3]
                       ! returned in [H T-1 ~> m s-1 or kg m-2 s-1]

end subroutine find_ustar_fluxes
module subroutine find_ustar_mech_forcing(forces, tv, U_star, G, GV, US, halo, H_T_units)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(mech_forcing),      intent(in)  :: forces !< Surface forces container
  type(thermo_var_ptrs),   intent(in)  :: tv   !< Structure containing pointers to any
                                               !! available thermodynamic fields.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(out) :: U_star !< The surface friction velocity [Z T-1 ~> m s-1]
  integer,       optional, intent(in)  :: halo !< The extra halo size to fill in, 0 by default
  logical,       optional, intent(in)  :: H_T_units !< If present and true, return U_star in units
                                               !! of [H T-1 ~> m s-1 or kg m-2 s-1]

  ! Local variables
                       ! the rescaled reference density [H2 Z-2 R-1 ~> m3 kg-1 or kg m-3]
                       ! returned in [H T-1 ~> m s-1 or kg m-2 s-1]

end subroutine find_ustar_mech_forcing
module subroutine MOM_forcing_chksum(mesg, fluxes, G, US, haloshift)
  character(len=*),        intent(in) :: mesg      !< message
  type(forcing),           intent(in) :: fluxes    !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in) :: G         !< grid type
  type(unit_scale_type),   intent(in) :: US        !< A dimensional unit scaling type
  integer, optional,       intent(in) :: haloshift !< shift in halo


end subroutine MOM_forcing_chksum
module subroutine MOM_mech_forcing_chksum(mesg, forces, G, US, haloshift)
  character(len=*),        intent(in) :: mesg      !< message
  type(mech_forcing),      intent(in) :: forces    !< A structure with the driving mechanical forces
  type(ocean_grid_type),   intent(in) :: G         !< grid type
  type(unit_scale_type),   intent(in) :: US        !< A dimensional unit scaling type
  integer, optional,       intent(in) :: haloshift !< shift in halo


end subroutine MOM_mech_forcing_chksum
module subroutine mech_forcing_SinglePointPrint(forces, G, i, j, mesg)
  type(mech_forcing),    intent(in) :: forces !< A structure with the driving mechanical forces
  type(ocean_grid_type), intent(in) :: G      !< Grid type
  character(len=*),      intent(in) :: mesg   !< Message
  integer,               intent(in) :: i      !< i-index
  integer,               intent(in) :: j      !< j-index

end subroutine mech_forcing_SinglePointPrint
module subroutine forcing_SinglePointPrint(fluxes, G, i, j, mesg)
  type(forcing),         intent(in) :: fluxes !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type), intent(in) :: G      !< Grid type
  character(len=*),      intent(in) :: mesg   !< Message
  integer,               intent(in) :: i      !< i-index
  integer,               intent(in) :: j      !< j-index

end subroutine forcing_SinglePointPrint
module subroutine register_forcing_type_diags(Time, diag, US, use_temperature, handles, use_berg_fluxes, use_waves, &
                                       use_cfcs, use_glc_runoff, use_carbon_runoff)
  type(time_type),     intent(in)    :: Time            !< time type
  type(diag_ctrl),     intent(inout) :: diag            !< diagnostic control type
  type(unit_scale_type), intent(in)  :: US              !< A dimensional unit scaling type
  logical,             intent(in)    :: use_temperature !< True if T/S are in use
  type(forcing_diags), intent(inout) :: handles         !< handles for diagnostics
  logical, optional,   intent(in)    :: use_berg_fluxes !< If true, allow iceberg flux diagnostics
  logical, optional,   intent(in)    :: use_waves       !< If true, allow wave forcing diagnostics
  logical, optional,   intent(in)    :: use_cfcs        !< If true, allow cfc related diagnostics
  logical, optional,   intent(in)    :: use_glc_runoff  !< If true, allow separate glacial runoff diagnostics
  logical, optional,   intent(in)    :: use_carbon_runoff  !< If true, allow separate carbon runoff diagnostics

  ! Clock for forcing diagnostics
end subroutine register_forcing_type_diags
module subroutine forcing_accumulate(flux_tmp, forces, fluxes, G, wt2)
  type(forcing),         intent(in)    :: flux_tmp !< A temporary structure with current
                                                 !!thermodynamic forcing fields
  type(mech_forcing),    intent(in)    :: forces !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes !< A structure containing time-averaged
                                                 !! thermodynamic forcing fields
  type(ocean_grid_type), intent(inout) :: G      !< The ocean's grid structure
  real,                  intent(out)   :: wt2    !< The relative weight of the new fluxes [nondim]

  ! This subroutine copies mechancal forcing from flux_tmp to fluxes and
  ! stores the time-weighted averages of the various buoyancy fluxes in fluxes,
  ! and increments the amount of time over which the buoyancy forcing should be
  ! applied, all via a call to fluxes accumulate.

end subroutine forcing_accumulate
module subroutine fluxes_accumulate(flux_tmp, fluxes, G, wt2, forces)
  type(forcing),             intent(in)    :: flux_tmp !< A temporary structure with current
                                                     !! thermodynamic forcing fields
  type(forcing),             intent(inout) :: fluxes !< A structure containing time-averaged
                                                     !! thermodynamic forcing fields
  type(ocean_grid_type),     intent(inout) :: G      !< The ocean's grid structure
  real,                      intent(out)   :: wt2    !< The relative weight of the new fluxes [nondim]
  type(mech_forcing), optional, intent(in) :: forces !< A structure with the driving mechanical forces

  ! This subroutine copies mechanical forcing from flux_tmp to fluxes and
  ! stores the time-weighted averages of the various buoyancy fluxes in fluxes,
  ! and increments the amount of time over which the buoyancy forcing in fluxes should be
  ! applied based on the time interval stored in flux_tmp.

end subroutine fluxes_accumulate
module subroutine copy_common_forcing_fields(forces, fluxes, G, skip_pres)
  type(mech_forcing),      intent(in)    :: forces   !< A structure with the driving mechanical forces
  type(forcing),           intent(inout) :: fluxes   !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in)    :: G        !< grid type
  logical,       optional, intent(in)    :: skip_pres !< If present and true, do not copy pressure fields.

end subroutine copy_common_forcing_fields
module subroutine set_derived_forcing_fields(forces, fluxes, G, US, Rho0)
  type(mech_forcing),      intent(in)    :: forces   !< A structure with the driving mechanical forces
  type(forcing),           intent(inout) :: fluxes   !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in)    :: G        !< grid type
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  real,                    intent(in)    :: Rho0     !< A reference density of seawater [R ~> kg m-3],
                                                     !! as used to calculate ustar.

end subroutine set_derived_forcing_fields
module subroutine set_net_mass_forcing(fluxes, forces, G, US)
  type(forcing),           intent(in)    :: fluxes   !< A structure containing thermodynamic forcing fields
  type(mech_forcing),      intent(inout) :: forces   !< A structure with the driving mechanical forces
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  type(ocean_grid_type),   intent(in)    :: G        !< The ocean grid type

end subroutine set_net_mass_forcing
module subroutine get_net_mass_forcing(fluxes, G, US, net_mass_src)
  type(forcing),                    intent(in)  :: fluxes !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),            intent(in)  :: G      !< The ocean grid type
  type(unit_scale_type),            intent(in)  :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: net_mass_src !< The net mass flux of water into the ocean
                                                          !! [R Z T-1 ~> kg m-2 s-1].

end subroutine get_net_mass_forcing
module subroutine copy_back_forcing_fields(fluxes, forces, G)
  type(forcing),           intent(in)    :: fluxes   !< A structure containing thermodynamic forcing fields
  type(mech_forcing),      intent(inout) :: forces   !< A structure with the driving mechanical forces
  type(ocean_grid_type),   intent(in)    :: G        !< grid type

end subroutine copy_back_forcing_fields
module subroutine mech_forcing_diags(forces_in, dt, G, time_end, diag, handles)
  type(mech_forcing), target, intent(in) :: forces_in !< mechanical forcing input fields
  real,                  intent(in)    :: dt       !< time step for the forcing [T ~> s]
  type(ocean_grid_type), intent(in)    :: G        !< grid type
  type(time_type),       intent(in)    :: time_end !< The end time of the diagnostic interval.
  type(diag_ctrl),       intent(inout) :: diag     !< diagnostic type
  type(forcing_diags),   intent(inout) :: handles  !< diagnostic id for diag_manager



end subroutine mech_forcing_diags
module subroutine forcing_diagnostics(fluxes_in, sfc_state, G_in, US, time_end, diag, handles, enthalpy)
  type(forcing), target, intent(in)    :: fluxes_in !< A structure containing thermodynamic forcing fields
  type(surface),         intent(in)    :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(ocean_grid_type), target, intent(in) :: G_in !< Input grid type
  type(unit_scale_type), intent(in)    :: US        !< A dimensional unit scaling type
  type(time_type),       intent(in)    :: time_end  !< The end time of the diagnostic interval.
  type(diag_ctrl),       intent(inout) :: diag      !< diagnostic regulator
  type(forcing_diags),   intent(inout) :: handles   !< diagnostic ids
  logical,     optional, intent(in   ) :: enthalpy  !< If present and true, the heat content associated
                                                    !! with mass entering/leaving the ocean is provided
                                                    !! by the coupler. Diagnostics net_heat_surface and
                                                    !! heat_content_surfwater are computed using
                                                    !! heat_content_evap instead of heat_content_massout.

  ! local variables
                          ! of mass fluxes [R Z T-1 ~> kg m-2 s-1] or heat fluxes [Q R Z T-1 ~> W m-2]

end subroutine forcing_diagnostics
module subroutine allocate_forcing_by_group(G, fluxes, water, heat, ustar, press, &
                                  shelf, iceberg, salt, fix_accum_bug, cfc, marbl, &
                                  waves, shelf_sfc_accumulation, lamult, hevap, &
                                  ice_ncat, tau_mag, carbon)
  type(ocean_grid_type), intent(in) :: G       !< Ocean grid structure
  type(forcing),      intent(inout) :: fluxes  !< A structure containing thermodynamic forcing fields
  logical, optional,     intent(in) :: water   !< If present and true, allocate water fluxes
  logical, optional,     intent(in) :: heat    !< If present and true, allocate heat fluxes
  logical, optional,     intent(in) :: ustar   !< If present and true, allocate ustar and related fields
  logical, optional,     intent(in) :: press   !< If present and true, allocate p_surf and related fields
  logical, optional,     intent(in) :: shelf   !< If present and true, allocate fluxes for ice-shelf
  logical, optional,     intent(in) :: iceberg !< If present and true, allocate fluxes for icebergs
  logical, optional,     intent(in) :: salt    !< If present and true, allocate salt fluxes
  logical, optional,     intent(in) :: fix_accum_bug !< If present and true, avoid using a bug in
                                               !! accumulation of ustar_gustless
  logical, optional,     intent(in) :: cfc     !< If present and true, allocate fields needed
                                               !! for cfc surface fluxes
  logical, optional,     intent(in) :: marbl   !< If present and true, allocate fields needed
                                               !! for MARBL surface fluxes
  logical, optional,     intent(in) :: waves   !< If present and true, allocate wave fields
  logical, optional,     intent(in) :: shelf_sfc_accumulation !< If present and true, and shelf is true,
                                               !! then allocate surface flux deposition from the atmosphere
                                               !! over ice shelves and ice sheets.
  logical, optional,     intent(in) :: lamult  !< If present and true, allocate langmuir enhancement factor
  logical, optional,     intent(in) :: hevap   !< If present and true, allocate heat content evap.
                                               !! This field must be allocated when enthalpy is provided
                                               !! via coupler.
  integer, optional,     intent(in) :: ice_ncat !< number of ice categories
  logical, optional,     intent(in) :: tau_mag !< If present and true, allocate tau_mag and related fields
  logical, optional,     intent(in) :: carbon  !< If present and true, allocate carbon fluxes

  ! Local variables

  ! if true, allocate fluxes needed to calculate enthalpy terms in MOM6
end subroutine allocate_forcing_by_group
module subroutine allocate_forcing_by_ref(fluxes_ref, G, fluxes, turns)
  type(forcing),         intent(in)  :: fluxes_ref !< Reference fluxes
  type(ocean_grid_type), intent(in)  :: G          !< Grid metric of target fluxes
  type(forcing),         intent(out) :: fluxes     !< Target fluxes
  integer,     optional, intent(in)  :: turns      !< If present, the number of counterclockwise
                                                   !! quarter turns to use on the new grid.


end subroutine allocate_forcing_by_ref
module subroutine allocate_mech_forcing_by_group(G, forces, stress, ustar, shelf, &
                                          press, iceberg, waves, num_stk_bands, tau_mag)
  type(ocean_grid_type), intent(in) :: G       !< Ocean grid structure
  type(mech_forcing), intent(inout) :: forces  !< Forcing fields structure

  logical, optional,     intent(in) :: stress  !< If present and true, allocate taux, tauy
  logical, optional,     intent(in) :: ustar   !< If present and true, allocate ustar and related fields
  logical, optional,     intent(in) :: shelf   !< If present and true, allocate forces for ice-shelf
  logical, optional,     intent(in) :: press   !< If present and true, allocate p_surf and related fields
  logical, optional,     intent(in) :: iceberg !< If present and true, allocate forces for icebergs
  logical, optional,     intent(in) :: waves   !< If present and true, allocate wave fields
  integer, optional,     intent(in) :: num_stk_bands !< Number of Stokes bands to allocate
  logical, optional,     intent(in) :: tau_mag !< If present and true, allocate tau_mag

  ! Local variables

end subroutine allocate_mech_forcing_by_group
module subroutine allocate_mech_forcing_from_ref(forces_ref, G, forces)
  type(mech_forcing), intent(in) :: forces_ref  !< Reference forcing fields
  type(ocean_grid_type), intent(in) :: G      !< Grid metric of target forcing
  type(mech_forcing), intent(out) :: forces   !< Mechanical forcing fields


  ! Identify the active fields in the reference forcing
end subroutine allocate_mech_forcing_from_ref
module subroutine get_forcing_groups(fluxes, water, heat, ustar, tau_mag, press, shelf, &
                             iceberg, salt, heat_added, buoy, carbon)
  type(forcing), intent(in) :: fluxes  !< Reference flux fields
  logical, intent(out) :: water   !< True if fluxes contains water-based fluxes
  logical, intent(out) :: heat    !< True if fluxes contains heat-based fluxes
  logical, intent(out) :: ustar   !< True if fluxes contains ustar
  logical, intent(out) :: tau_mag !< True if fluxes contains tau_mag
  logical, intent(out) :: press   !< True if fluxes contains surface pressure
  logical, intent(out) :: shelf   !< True if fluxes contains ice shelf fields
  logical, intent(out) :: iceberg !< True if fluxes contains iceberg fluxes
  logical, intent(out) :: salt    !< True if fluxes contains salt flux
  logical, intent(out) :: heat_added !< True if fluxes contains explicit heat
  logical, intent(out) :: buoy    !< True if fluxes contains buoyancy fluxes
  logical, optional, intent(out) :: carbon  !< True if fluxes contains carbon fluxes

  ! NOTE: heat, salt, heat_added, and buoy would typically depend on each other
  !   to some degree.  But since this would be enforced at the driver level,
  !   we handle them here as independent flags.

end subroutine get_forcing_groups
module subroutine get_mech_forcing_groups(forces, stress, ustar, tau_mag, shelf, press, iceberg)
  type(mech_forcing), intent(in) :: forces  !< Reference forcing fields
  logical, intent(out) :: stress  !< True if forces contains wind stress fields
  logical, intent(out) :: ustar   !< True if forces contains ustar field
  logical, intent(out) :: tau_mag !< True if forces contains tau_mag field
  logical, intent(out) :: shelf   !< True if forces contains ice shelf fields
  logical, intent(out) :: press   !< True if forces contains pressure fields
  logical, intent(out) :: iceberg !< True if forces contains iceberg fields

end subroutine get_mech_forcing_groups
module subroutine myAlloc_2d(array, is, ie, js, je, flag)
  real, dimension(:,:), pointer :: array !< Array to be allocated
  integer,           intent(in) :: is !< Start i-index
  integer,           intent(in) :: ie !< End i-index
  integer,           intent(in) :: js !< Start j-index
  integer,           intent(in) :: je !< End j-index
  logical, optional, intent(in) :: flag !< Flag to indicate to allocate

end subroutine myAlloc_2d
module subroutine myAlloc_3d(array, is, ie, js, je, ks, ke, flag)
  real, dimension(:,:,:), pointer :: array !< Array to be allocated
  integer,             intent(in) :: is !< Start i-index
  integer,             intent(in) :: ie !< End i-index
  integer,             intent(in) :: js !< Start j-index
  integer,             intent(in) :: je !< End j-index
  integer,             intent(in) :: ks !< Start k-index
  integer,             intent(in) :: ke !< End k-index
  logical, optional,   intent(in) :: flag !< Flag to indicate to allocate

end subroutine myAlloc_3d
module subroutine deallocate_forcing_type(fluxes)
  type(forcing), intent(inout) :: fluxes !< Forcing fields structure

end subroutine deallocate_forcing_type
module subroutine deallocate_mech_forcing(forces)
  type(mech_forcing), intent(inout) :: forces  !< Forcing fields structure

end subroutine deallocate_mech_forcing
module subroutine rotate_forcing(fluxes_in, fluxes, turns)
  type(forcing), intent(in)  :: fluxes_in     !< Input forcing structure
  type(forcing), intent(inout) :: fluxes      !< Rotated forcing structure
  integer, intent(in) :: turns                !< Number of quarter turns


end subroutine rotate_forcing
module subroutine rotate_mech_forcing(forces_in, turns, forces)
  type(mech_forcing), intent(in)  :: forces_in  !< Forcing on the input domain
  integer, intent(in) :: turns                  !< Number of quarter-turns
  type(mech_forcing), intent(inout) :: forces   !< Forcing on the rotated domain


end subroutine rotate_mech_forcing
module subroutine homogenize_mech_forcing(forces, G, US, Rho0, UpdateUstar)
  type(mech_forcing),    intent(inout) :: forces !< Forcing on the input domain
  type(ocean_grid_type),    intent(in) :: G      !< Grid metric of target forcing
  type(unit_scale_type),    intent(in) :: US     !< A dimensional unit scaling type
  real,                     intent(in) :: Rho0   !< A reference density of seawater [R ~> kg m-3],
                                                 !! as used to calculate ustar.
  logical, optional,        intent(in) :: UpdateUstar !< A logical to determine if Ustar should be directly averaged
                                                 !! or updated from mean tau.

end subroutine homogenize_mech_forcing
module subroutine homogenize_forcing(fluxes, G, GV, US)
  type(forcing),           intent(inout) :: fluxes !< Input forcing struct
  type(ocean_grid_type),   intent(in)    :: G      !< Grid metric of target forcing
  type(verticalGrid_type), intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type


end subroutine homogenize_forcing
module subroutine homogenize_field_t(var, G, tmp_scale)
  type(ocean_grid_type),            intent(in)    :: G   !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: var !< The variable to homogenize [A ~> a]
  real,                    optional, intent(in)    :: tmp_scale !< A temporary rescaling factor for the
                                                         !! variable that is reversed in the
                                                         !! return value [a A-1 ~> 1]

end subroutine homogenize_field_t
module subroutine homogenize_field_v(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: var  !< The variable to homogenize [A ~> a]
  real,                    optional, intent(in)    :: tmp_scale !< A temporary rescaling factor for the
                                                           !! variable that is reversed in the
                                                           !! return value [a A-1 ~> 1]

end subroutine homogenize_field_v
module subroutine homogenize_field_u(var, G, tmp_scale)
  type(ocean_grid_type),             intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: var  !< The variable to homogenize [A ~> a]
  real,                    optional, intent(in)    :: tmp_scale !< A temporary rescaling factor for the
                                                           !! variable that is reversed in the
                                                           !! return value [a A-1 ~> 1]

end subroutine homogenize_field_u
  end interface

end module MOM_forcing_type
