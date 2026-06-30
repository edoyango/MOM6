! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides functions for some diabatic processes such as frazil, brine rejection,
!! tendency due to surface flux divergence.
module MOM_diabatic_aux

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS,           only : calculate_density, calculate_TFreeze, EOS_domain
use MOM_EOS,           only : calculate_specific_vol_derivs, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, callTree_showQuery
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing, extractFluxes1d, forcing_SinglePointPrint
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_interpolate,   only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate,   only : external_field
use MOM_io,            only : slasher
use MOM_opacity,       only : set_opacity, opacity_CS, extract_optics_slice, extract_optics_fields
use MOM_opacity,       only : optics_type, optics_nbands, absorbRemainingSW, sumSWoverBands
use MOM_tracer_flow_control, only : get_chl_from_model, tracer_flow_control_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public diabatic_aux_init, diabatic_aux_end
public make_frazil, adjust_salt, differential_diffuse_T_S, triDiagTS, triDiagTS_Eulerian
public find_uv_at_h, applyBoundaryFluxesInOut, set_pen_shortwave

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for diabatic_aux
type, public :: diabatic_aux_CS ; private
  logical :: do_rivermix = .false. !< Provide additional TKE to mix river runoff at the
                                   !! river mouths to a depth of "rivermix_depth"
  real    :: rivermix_depth = 0.0  !< The depth to which rivers are mixed if do_rivermix = T [Z ~> m].
  real    :: dSalt_frac_max  !< An upper limit on the fraction of the salt in a layer that can be
                             !! lost to the net surface salt fluxes within a timestep [nondim]
  logical :: reclaim_frazil  !<   If true, try to use any frazil heat deficit to
                             !! to cool the topmost layer down to the freezing
                             !! point.  The default is true.
  logical :: pressure_dependent_frazil  !< If true, use a pressure dependent
                             !! freezing temperature when making frazil.  The
                             !! default is false, which will be faster but is
                             !! inappropriate with ice-shelf cavities.
  logical :: ignore_fluxes_over_land    !< If true, the model does not check
                             !! if fluxes are applied over land points. This
                             !! flag must be used when the ocean is coupled with
                             !! sea ice and ice shelves and use_ePBL = true.
  logical :: use_river_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a river heat content. Otherwise, runoff
                             !! is added with a temperature of the local SST.
  logical :: use_calving_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a calving heat content. Otherwise, calving
                             !! is added with a temperature of the local SST.
  logical :: var_pen_sw      !<   If true, use one of the CHL_A schemes to determine the
                             !! e-folding depth of incoming shortwave radiation.
  type(external_field) :: sbc_chl   !< A handle used in time interpolation of
                             !! chlorophyll read from a file.
  logical :: chl_from_file   !< If true, chl_a is read from a file.
  logical :: do_brine_plume  !< If true, insert salt flux below the surface according to
                             !! a parameterization by \cite Nguyen2009.
  logical :: check_salt_bp   !< A logical to check for salt conservation in the brine plume scheme
  !TODO: Delete DEBUG lines after brine plume is proven to be conservative to numerical precision.
  !DEBUG logical :: check_salt_verbose !< A logical to be verbose when checking salt conservation
  integer :: brine_plume_n   !< The exponent in the brine plume parameterization.
  real :: plume_strength     !< Fraction of the available brine to take to the bottom of the mixed
                             !! layer [nondim].
  real :: plume_mld_fac      !< Proportionality factor between the mixed/mixing layer depth and the
                             !! vertical scale used for the brine plume parameterization [nondim].
  real :: check_salt_threshold!< The maximum relative salt change acceptable in a time step [nondim]

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< Structure used to regulate timing of diagnostic output

  ! Diagnostic handles
  integer :: id_createdH       = -1 !< Diagnostic ID of mass added to avoid grounding
  integer :: id_brine_input    = -1 !< Diagnostic ID of which layer receives the brine salt flux
  integer :: id_penSW_diag     = -1 !< Diagnostic ID of Penetrative shortwave heating (flux convergence)
  integer :: id_penSWflux_diag = -1 !< Diagnostic ID of Penetrative shortwave flux
  integer :: id_nonpenSW_diag  = -1 !< Diagnostic ID of Non-penetrative shortwave heating
  integer :: id_Chl            = -1 !< Diagnostic ID of chlorophyll-A handles for opacity

  ! Optional diagnostic arrays
  real, allocatable, dimension(:,:)   :: createdH       !< The amount of volume added in order to
                                                        !! avoid grounding [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: brine_input    !< Brine input diagnostic indicating
                                                        !! the resulting salt tendency [S T-1 ~> ppt s-1]
  real, allocatable, dimension(:,:,:) :: penSW_diag     !< Heating in a layer from convergence of
                                                        !! penetrative SW [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:,:) :: penSWflux_diag !< Penetrative SW flux at base of grid
                                                        !! layer [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:)   :: nonpenSW_diag  !< Non-downwelling SW radiation at ocean
                                                        !! surface [Q R Z T-1 ~> W m-2]

end type diabatic_aux_CS

!>@{ CPU time clock IDs
integer :: id_clock_uv_at_h, id_clock_frazil
!>@}


  interface
module subroutine make_frazil(h, tv, G, GV, US, CS, p_surf, halo)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Structure containing pointers to any available
                                               !! thermodynamic fields.
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(diabatic_aux_CS),   intent(in)    :: CS !< The control structure returned by a previous
                                               !! call to diabatic_aux_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa].
  integer,       optional, intent(in)    :: halo !< Halo width over which to calculate frazil
  ! Local variables
                       ! row of points.

end subroutine make_frazil
module subroutine differential_diffuse_T_S(h, T, S, Kd_T, Kd_S, tv, dt, G, GV)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T    !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S    !< Salinity [PSU] or [gSalt/kg], generically [S ~> ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)    :: Kd_T !< The extra diffusivity of temperature due to
                                                 !! double diffusion relative to the diffusivity of
                                                 !! density [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)    :: Kd_S !< The extra diffusivity of salinity due to
                                                 !! double diffusion relative to the diffusivity of
                                                 !! density [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< Structure containing pointers to any
                                                 !! available thermodynamic fields.
  real,                    intent(in)    :: dt   !<  Time increment [T ~> s].

  ! local variables
                    ! added to ensure positive definiteness [H ~> m or kg m-2].
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
                    ! in roundoff and can be neglected [Z ~> m].

end subroutine differential_diffuse_T_S
module subroutine adjust_salt(h, tv, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< Structure containing pointers to any
                                                 !! available thermodynamic fields.
  type(diabatic_aux_CS),   intent(in)    :: CS   !< The control structure returned by a previous
                                                 !! call to diabatic_aux_init.

  ! local variables

end subroutine adjust_salt
module subroutine triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, T, S)
  type(ocean_grid_type),                     intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: is !< The start i-index to work on.
  integer,                                   intent(in)    :: ie !< The end i-index to work on.
  integer,                                   intent(in)    :: js !< The start j-index to work on.
  integer,                                   intent(in)    :: je !< The end j-index to work on.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: hold !< The layer thicknesses before entrainment,
                                                                 !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: ea !< The amount of fluid entrained from the layer
                                                                 !! above within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: eb !< The amount of fluid entrained from the layer
                                                                 !! below within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T  !< Layer potential temperatures [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S  !< Layer salinities [S ~> ppt].

  ! Local variables

  !$OMP parallel do default(shared) private(h_tr,b1,d1,c1,b_denom_1)
end subroutine triDiagTS
module subroutine triDiagTS_Eulerian(G, GV, is, ie, js, je, hold, ent, T, S)
  type(ocean_grid_type),                     intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: is   !< The start i-index to work on.
  integer,                                   intent(in)    :: ie   !< The end i-index to work on.
  integer,                                   intent(in)    :: js   !< The start j-index to work on.
  integer,                                   intent(in)    :: je   !< The end j-index to work on.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: hold !< The layer thicknesses before entrainment,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: ent  !< The amount of fluid mixed across an interface
                                                                   !! within this time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T    !< Layer potential temperatures [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S    !< Layer salinities [S ~> ppt].

  ! Local variables

  !$OMP parallel do default(shared) private(h_tr,b1,d1,c1,b_denom_1)
end subroutine triDiagTS_Eulerian
module subroutine find_uv_at_h(u, v, h, u_h, v_h, G, GV, US, ea, eb, zero_mix)
  type(ocean_grid_type),     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)  :: u    !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             intent(in)  :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(out)   :: u_h !< Zonal velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(out)   :: v_h !< Meridional velocity interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)  :: ea !< The amount of fluid entrained from the layer
                                                 !! above within this time step [H ~> m or kg m-2].
                                                 !! Omitting ea is the same as setting it to 0.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)  :: eb !< The amount of fluid entrained from the layer
                                                 !! below within this time step [H ~> m or kg m-2].
                                                 !! Omitting eb is the same as setting it to 0.
  logical,           optional, intent(in)  :: zero_mix !< If true, do the calculation of u_h and
                                                 !! v_h as though ea and eb were being supplied with
                                                 !! uniformly zero values.

  ! Local variables
                       ! in roundoff and can be neglected [H ~> m or kg m-2].
  ! Fractional weights of the neighboring velocity points, ~1/2 in the open ocean.
end subroutine find_uv_at_h
module subroutine set_pen_shortwave(optics, fluxes, G, GV, US, CS, opacity, tracer_flow_CSp)
  type(optics_type),       pointer       :: optics !< An optics structure that has will contain
                                                   !! information about shortwave fluxes and absorption.
  type(forcing),           intent(inout) :: fluxes !< points to forcing fields
                                                   !! unused fields have NULL pointers
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(diabatic_aux_CS),   pointer       :: CS     !< Control structure for diabatic_aux
  type(opacity_CS)                       :: opacity !< The control structure for the opacity module.
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp !< A pointer to the control structure
                                                   !! organizing the tracer modules.

  ! Local variables
end subroutine set_pen_shortwave
module subroutine applyBoundaryFluxesInOut(CS, G, GV, US, dt, fluxes, optics, nsw, h, tv, &
                                    aggregate_FW_forcing, evap_CFL_limit, &
                                    minimum_forcing_depth, cTKE, dSV_dT, dSV_dS, &
                                    SkinBuoyFlux, MLD_h)
  type(diabatic_aux_CS),   pointer       :: CS !< Control structure for diabatic_aux
  type(ocean_grid_type),   intent(in)    :: G  !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real,                    intent(in)    :: dt !< Time-step over which forcing is applied [T ~> s]
  type(forcing),           intent(inout) :: fluxes !< Surface fluxes container
  type(optics_type),       pointer       :: optics !< Optical properties container
  integer,                 intent(in)    :: nsw !< The number of frequency bands of penetrating
                                                !! shortwave radiation
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Structure containing pointers to any
                                               !! available thermodynamic fields.
  logical,                 intent(in)    :: aggregate_FW_forcing !< If False, treat in/out fluxes separately.
  real,                    intent(in)    :: evap_CFL_limit !< The largest fraction of a layer that
                                               !! can be evaporated in one time-step [nondim].
  real,                    intent(in)    :: minimum_forcing_depth !< The smallest depth over which
                                               !! heat and freshwater fluxes is applied [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: cTKE !< Turbulent kinetic energy requirement to mix
                                               !! forcing through each layer [R Z3 T-2 ~> J m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: dSV_dT !< Partial derivative of specific volume with
                                               !! potential temperature [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: dSV_dS !< Partial derivative of specific volume with
                                               !! salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out)   :: SkinBuoyFlux !< Buoyancy flux at surface [Z2 T-3 ~> m2 s-3].
  real, dimension(:,:), &
                 optional, pointer       :: MLD_h !< Mixed layer thickness for brine plumes [H ~> m or kg m-2]

  ! Local variables
                     ! drops below this value [H ~> m or kg m-2]
                     ! shifted to the next deeper layer [H ~> m or kg m-2]
                         ! By default EnthalpyConst = 1.0. If fluxes%heat_content_evap
                         ! is associated enthalpy is provided via coupler and EnthalpyConst = 0.0.
end subroutine applyBoundaryFluxesInOut
module subroutine diabatic_aux_init(Time, G, GV, US, param_file, diag, CS, useALEalgorithm, use_ePBL)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure used to regulate diagnostic output
  type(diabatic_aux_CS),   pointer       :: CS   !< A pointer to the control structure for the
                                                 !! diabatic_aux module, which is initialized here.
  logical,                 intent(in)    :: useALEalgorithm !< If true, use the ALE algorithm rather
                                                 !! than layered mode.
  logical,                 intent(in)    :: use_ePBL !< If true, use the implicit energetics planetary
                                                 !! boundary layer scheme to determine the diffusivity
                                                 !! in the surface boundary layer.

  ! This "include" declares and sets the variable "version".
                                 ! when var_pen_sw is defined and reading from file.
end subroutine diabatic_aux_init
module subroutine diabatic_aux_end(CS)
  type(diabatic_aux_CS), pointer :: CS !< The control structure returned by a previous
                                       !! call to diabatic_aux_init; it is deallocated here.

end subroutine diabatic_aux_end
  end interface

end module MOM_diabatic_aux
