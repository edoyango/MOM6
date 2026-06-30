! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Orchestrates the registration and calling of tracer packages
module MOM_tracer_flow_control

use MOM_coms,          only : EFP_type, assignment(=), EFP_to_real, real_to_EFP, EFP_sum_across_PEs
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type, close_param_file
use MOM_forcing_type,  only : forcing, optics_type
use MOM_get_input,     only : Get_MOM_input
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_interface_heights, only : convert_MLD_to_ML_thickness
use MOM_CVMix_KPP,     only : KPP_CS
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart,       only : MOM_restart_CS
use MOM_sponge,        only : sponge_CS
use MOM_ALE_sponge,    only : ALE_sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
#include <MOM_memory.h>

! Add references to other user-provide tracer modules here.
use USER_tracer_example, only : tracer_column_physics, USER_initialize_tracer, USER_tracer_stock
use USER_tracer_example, only : USER_register_tracer_example, USER_tracer_surface_state
use USER_tracer_example, only : USER_tracer_example_end, USER_tracer_example_CS
use DOME_tracer, only : register_DOME_tracer, initialize_DOME_tracer
use DOME_tracer, only : DOME_tracer_column_physics, DOME_tracer_surface_state
use DOME_tracer, only : DOME_tracer_end, DOME_tracer_CS
use ISOMIP_tracer, only : register_ISOMIP_tracer, initialize_ISOMIP_tracer
use ISOMIP_tracer, only : ISOMIP_tracer_column_physics, ISOMIP_tracer_surface_state
use ISOMIP_tracer, only : ISOMIP_tracer_end, ISOMIP_tracer_CS
use RGC_tracer, only : register_RGC_tracer, initialize_RGC_tracer
use RGC_tracer, only : RGC_tracer_column_physics
use RGC_tracer, only : RGC_tracer_end, RGC_tracer_CS
use ideal_age_example, only : register_ideal_age_tracer, initialize_ideal_age_tracer
use ideal_age_example, only : ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
use ideal_age_example, only : ideal_age_stock, ideal_age_example_end, ideal_age_tracer_CS
use MARBL_tracers, only : register_MARBL_tracers, initialize_MARBL_tracers
use MARBL_tracers, only : MARBL_tracers_column_physics, MARBL_tracers_set_forcing
use MARBL_tracers, only : MARBL_tracers_surface_state, MARBL_tracers_get
use MARBL_tracers, only : MARBL_tracers_stock, MARBL_tracers_end, MARBL_tracers_CS
use regional_dyes, only : register_dye_tracer, initialize_dye_tracer
use regional_dyes, only : dye_tracer_column_physics, dye_tracer_surface_state
use regional_dyes, only : dye_stock, regional_dyes_end, dye_tracer_CS
use MOM_OCMIP2_CFC, only : register_OCMIP2_CFC, initialize_OCMIP2_CFC, flux_init_OCMIP2_CFC
use MOM_OCMIP2_CFC, only : OCMIP2_CFC_column_physics, OCMIP2_CFC_surface_state
use MOM_OCMIP2_CFC, only : OCMIP2_CFC_stock, OCMIP2_CFC_end, OCMIP2_CFC_CS
use MOM_CFC_cap, only : register_CFC_cap, initialize_CFC_cap
use MOM_CFC_cap, only : CFC_cap_column_physics, CFC_cap_set_forcing
use MOM_CFC_cap, only : CFC_cap_stock, CFC_cap_end, CFC_cap_CS
use oil_tracer, only : register_oil_tracer, initialize_oil_tracer
use oil_tracer, only : oil_tracer_column_physics, oil_tracer_surface_state
use oil_tracer, only : oil_stock, oil_tracer_end, oil_tracer_CS
use advection_test_tracer, only : register_advection_test_tracer, initialize_advection_test_tracer
use advection_test_tracer, only : advection_test_tracer_column_physics, advection_test_tracer_surface_state
use advection_test_tracer, only : advection_test_stock, advection_test_tracer_end, advection_test_tracer_CS
use dyed_obc_tracer, only : register_dyed_obc_tracer, initialize_dyed_obc_tracer
use dyed_obc_tracer, only : dyed_obc_tracer_column_physics
use dyed_obc_tracer, only : dyed_obc_tracer_end, dyed_obc_tracer_CS
use MOM_generic_tracer, only : register_MOM_generic_tracer, initialize_MOM_generic_tracer
use MOM_generic_tracer, only : MOM_generic_tracer_column_physics, MOM_generic_tracer_surface_state
use MOM_generic_tracer, only : end_MOM_generic_tracer, MOM_generic_tracer_get, MOM_generic_flux_init
use MOM_generic_tracer, only : MOM_generic_tracer_stock, MOM_generic_tracer_min_max, MOM_generic_tracer_CS
use MOM_generic_tracer, only : register_MOM_generic_tracer_segments
use pseudo_salt_tracer, only : register_pseudo_salt_tracer, initialize_pseudo_salt_tracer
use pseudo_salt_tracer, only : pseudo_salt_tracer_column_physics, pseudo_salt_tracer_surface_state
use pseudo_salt_tracer, only : pseudo_salt_stock, pseudo_salt_tracer_end, pseudo_salt_tracer_CS
use boundary_impulse_tracer, only : register_boundary_impulse_tracer, initialize_boundary_impulse_tracer
use boundary_impulse_tracer, only : boundary_impulse_tracer_column_physics, boundary_impulse_tracer_surface_state
use boundary_impulse_tracer, only : boundary_impulse_stock, boundary_impulse_tracer_end
use boundary_impulse_tracer, only : boundary_impulse_tracer_CS
use nw2_tracers, only : nw2_tracers_CS, register_nw2_tracers, nw2_tracer_column_physics
use nw2_tracers, only : initialize_nw2_tracers, nw2_tracers_end

implicit none ; private

public call_tracer_register, tracer_flow_control_init, call_tracer_set_forcing
public call_tracer_column_fns, call_tracer_surface_state, call_tracer_stocks
public call_tracer_flux_init, get_chl_from_model, tracer_flow_control_end
public call_tracer_register_obc_segments

!> The control structure for orchestrating the calling of tracer packages
type, public :: tracer_flow_control_CS ; private
  logical :: use_USER_tracer_example = .false.     !< If true, use the USER_tracer_example package
  logical :: use_DOME_tracer = .false.             !< If true, use the DOME_tracer package
  logical :: use_ISOMIP_tracer = .false.           !< If true, use the ISOMPE_tracer package
  logical :: use_RGC_tracer =.false.               !< If true, use the RGC_tracer package
  logical :: use_ideal_age = .false.               !< If true, use the ideal age tracer package
  logical :: use_MARBL_tracers = .false.           !< If true, use the MARBL tracer package
  logical :: use_regional_dyes = .false.           !< If true, use the regional dyes tracer package
  logical :: use_oil = .false.                     !< If true, use the oil tracer package
  logical :: use_advection_test_tracer = .false.   !< If true, use the advection_test_tracer package
  logical :: use_OCMIP2_CFC = .false.              !< If true, use the OCMIP2_CFC tracer package
  logical :: use_CFC_cap = .false.                 !< If true, use the CFC_cap tracer package
  logical :: use_MOM_generic_tracer = .false.      !< If true, use the MOM_generic_tracer packages
  logical :: use_pseudo_salt_tracer = .false.      !< If true, use the psuedo_salt tracer  package
  logical :: use_boundary_impulse_tracer = .false. !< If true, use the boundary impulse tracer package
  logical :: use_dyed_obc_tracer = .false.         !< If true, use the dyed OBC tracer package
  logical :: use_nw2_tracers = .false.             !< If true, use the NW2 tracer package
  logical :: get_chl_from_MARBL = .false.          !< If true, use the MARBL-provided Chl for shortwave penetration
  !>@{ Pointers to the control strucures for the tracer packages
  type(USER_tracer_example_CS), pointer :: USER_tracer_example_CSp => NULL()
  type(DOME_tracer_CS), pointer :: DOME_tracer_CSp => NULL()
  type(ISOMIP_tracer_CS), pointer :: ISOMIP_tracer_CSp => NULL()
  type(RGC_tracer_CS), pointer :: RGC_tracer_CSp => NULL()
  type(ideal_age_tracer_CS), pointer :: ideal_age_tracer_CSp => NULL()
  type(MARBL_tracers_CS), pointer :: MARBL_tracers_CSp => NULL()
  type(dye_tracer_CS), pointer :: dye_tracer_CSp => NULL()
  type(oil_tracer_CS), pointer :: oil_tracer_CSp => NULL()
  type(advection_test_tracer_CS), pointer :: advection_test_tracer_CSp => NULL()
  type(OCMIP2_CFC_CS), pointer :: OCMIP2_CFC_CSp => NULL()
  type(CFC_cap_CS),    pointer :: CFC_cap_CSp => NULL()
  type(MOM_generic_tracer_CS), pointer :: MOM_generic_tracer_CSp => NULL()
  type(pseudo_salt_tracer_CS), pointer :: pseudo_salt_tracer_CSp => NULL()
  type(boundary_impulse_tracer_CS), pointer :: boundary_impulse_tracer_CSp => NULL()
  type(dyed_obc_tracer_CS), pointer :: dyed_obc_tracer_CSp => NULL()
  type(nw2_tracers_CS), pointer :: nw2_tracers_CSp => NULL()
  !>@}
end type tracer_flow_control_CS


  interface
module subroutine call_tracer_flux_init(verbosity)
  integer, optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.


  ! Determine which tracer routines with tracer fluxes are to be called.  Note
  ! that not every tracer package is required to have a flux_init call.
end subroutine call_tracer_flux_init
module subroutine call_tracer_register(G, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),        intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),      intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),        intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time
                                                         !! parameters.
  type(tracer_flow_control_CS), pointer    :: CS         !< A pointer that is set to point to the
                                                         !! control structure for this module.
  type(tracer_registry_type),   pointer    :: tr_Reg     !< A pointer that is set to point to the
                                                         !! control structure for the tracer
                                                         !! advection and diffusion module.
  type(MOM_restart_CS), intent(inout) :: restart_CS      !< A pointer to the restart control
                                                         !! structure.

  ! This include declares and sets the variable "version".

end subroutine call_tracer_register
module subroutine tracer_flow_control_init(restart, day, G, GV, US, h, param_file, diag, OBC, &
                                    CS, sponge_CSp, ALE_sponge_CSp, tv)
  logical,                               intent(in)    :: restart !< 1 if the fields have already
                                                                  !! been read from a restart file.
  type(time_type), target,               intent(in)    :: day     !< Time of the start of the run.
  type(ocean_grid_type),                 intent(inout) :: G       !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in)    :: GV      !< The ocean's vertical grid
                                                                  !! structure.
  type(unit_scale_type),                 intent(in)    :: US      !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                         intent(in)    :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),                 intent(in)    :: param_file !< A structure to parse for
                                                                  !! run-time parameters
  type(diag_ctrl), target,               intent(in)    :: diag    !< A structure that is used to
                                                                  !! regulate diagnostic output.
  type(ocean_OBC_type),                  pointer       :: OBC     !< This open boundary condition
                                                                  !! type specifies whether, where,
                                                                  !! and what open boundary
                                                                  !! conditions are used.
  type(tracer_flow_control_CS),          pointer       :: CS      !< The control structure returned
                                                                  !! by a previous call to
                                                                  !! call_tracer_register.
  type(sponge_CS),                       pointer       :: sponge_CSp     !< A pointer to the control
                                               !! structure for the sponges, if they are in use.
                                               !! Otherwise this may be unassociated.
  type(ALE_sponge_CS),                   pointer       :: ALE_sponge_CSp !< A pointer to the control
                                               !! structure for the ALE sponges, if they are in use.
                                               !! Otherwise this may be unassociated.
  type(thermo_var_ptrs),                 intent(in)    :: tv      !< A structure pointing to various
                                                                  !! thermodynamic variables

end subroutine tracer_flow_control_init
module subroutine call_tracer_register_obc_segments(GV, param_file, CS, tr_Reg, OBC)
  type(verticalGrid_type),      intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time
                                                         !! parameters.
  type(tracer_flow_control_CS), pointer    :: CS         !< A pointer that is set to point to the
                                                         !! control structure for this module.
  type(tracer_registry_type),   pointer    :: tr_Reg     !< A pointer that is set to point to the
                                                         !! control structure for the tracer
                                                         !! advection and diffusion module.
  type(ocean_OBC_type),      pointer       :: OBC        !< This open boundary condition
                                                         !! type specifies whether, where,
                                                         !! and what open boundary
                                                         !! conditions are used.

end subroutine call_tracer_register_obc_segments
module subroutine get_chl_from_model(Chl_array, G, GV, CS)
  type(ocean_grid_type),        intent(in)  :: G         !< The ocean's grid structure.
  type(verticalGrid_type),      intent(in)  :: GV        !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                intent(out) :: Chl_array !< The array in which to store the model's
                                                         !! Chlorophyll-A concentrations [mg m-3].
  type(tracer_flow_control_CS), pointer     :: CS        !< The control structure returned by a
                                                         !! previous call to call_tracer_register.

end subroutine get_chl_from_model
module subroutine call_tracer_set_forcing(sfc_state, fluxes, day_start, day_interval, G, US, Rho0, CS)

  type(surface),                intent(inout) :: sfc_state !< A structure containing fields that
                                                           !! describe the surface state of the
                                                           !! ocean.
  type(forcing),                intent(inout) :: fluxes    !< A structure containing pointers to any
                                                           !! possible forcing fields. Unused fields
                                                           !! have NULL ptrs.
  type(time_type),              intent(in)    :: day_start !< Start time of the fluxes.
  type(time_type),              intent(in)    :: day_interval !< Length of time over which these
                                                           !! fluxes will be applied.
  type(ocean_grid_type),        intent(in)    :: G         !< The ocean's grid structure.
  type(unit_scale_type),        intent(in)    :: US        !< A dimensional unit scaling type
  real,                         intent(in)    :: Rho0      !< The mean ocean density [R ~> kg m-3]
  type(tracer_flow_control_CS), pointer       :: CS        !< The control structure returned by a
                                                           !! previous call to call_tracer_register.

end subroutine call_tracer_set_forcing
module subroutine call_tracer_column_fns(h_old, h_new, ea, eb, fluxes, mld, dt, G, GV, US, tv, optics, CS, &
                                  debug, KPP_CSp, nonLocalTrans, evap_CFL_limit, minimum_forcing_depth, h_BL)
  type(ocean_grid_type),                 intent(in) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in) :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old !< Layer thickness before entrainment
                                                              !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_new !< Layer thickness after entrainment
                                                              !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: ea !< an array to which the amount of
                                          !! fluid entrained from the layer above during this call
                                          !! will be added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: eb !< an array to which the amount of
                                          !! fluid entrained from the layer below during this call
                                          !! will be added [H ~> m or kg m-2].
  type(forcing),                         intent(in) :: fluxes !< A structure containing pointers to
                                                              !! any possible forcing fields.
                                                              !! Unused fields have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G)),      intent(in) :: mld    !< Mixed layer depth [Z ~> m]
  real,                                  intent(in) :: dt     !< The amount of time covered by this
                                                              !! call [T ~> s]
  type(unit_scale_type),                 intent(in) :: US     !< A dimensional unit scaling type
  type(thermo_var_ptrs),                 intent(in) :: tv     !< A structure pointing to various
                                                              !! thermodynamic variables.
  type(optics_type),                     pointer    :: optics !< The structure containing optical
                                                              !! properties.
  type(tracer_flow_control_CS),          pointer    :: CS     !< The control structure returned by
                                                              !! a previous call to
                                                              !! call_tracer_register.
  logical,                               intent(in) :: debug  !< If true calculate checksums
  type(KPP_CS),                optional, pointer    :: KPP_CSp  !< KPP control structure
  real,                        optional, intent(in) :: nonLocalTrans(:,:,:) !< Non-local transport [nondim]
  real,                        optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of
                                                              !! the water that can be fluxed out
                                                              !! of the top layer in a timestep [nondim]
  real,                        optional, intent(in) :: minimum_forcing_depth !< The smallest depth over
                                                              !! which fluxes can be applied [H ~> m or kg m-2]
  real, dimension(:,:),        optional, pointer    :: h_BL   !< Thickness of active mixing layer [H ~> m or kg m-2]

  ! Local variables

end subroutine call_tracer_column_fns
module subroutine call_tracer_stocks(h, stock_values, G, GV, US, CS, stock_names, stock_units, &
                              num_stocks, stock_index, got_min_max, global_min, global_max, &
                              xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)
  type(ocean_grid_type),          intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    &
                                  intent(in)  :: h           !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),             intent(out) :: stock_values !< The globally mass-integrated
                                                             !! amount of a tracer [kg conc].
  type(unit_scale_type),          intent(in)  :: US          !< A dimensional unit scaling type
  type(tracer_flow_control_CS),   pointer     :: CS          !< The control structure returned by a
                                                             !! previous call to
                                                             !! call_tracer_register.
  character(len=*), dimension(:), &
                        optional, intent(out) :: stock_names !< Diagnostic names to use for each stock.
  character(len=*), dimension(:), &
                        optional, intent(out) :: stock_units !< Units to use in the metadata for each stock.
  integer,              optional, intent(out) :: num_stocks  !< The number of tracer stocks being returned.
  integer,              optional, intent(in)  :: stock_index !< The integer stock index from
                             !! stocks_constants_mod of the stock to be returned.  If this is
                             !! present and greater than 0, only a single stock can be returned.
  logical, dimension(:), &
                      optional, intent(inout) :: got_min_max !< Indicates whether the global min and
                                                             !! max are found for each tracer
  real, dimension(:), optional, intent(out)   :: global_min  !< The global minimum of each tracer [conc]
  real, dimension(:), optional, intent(out)   :: global_max  !< The global maximum of each tracer [conc]
  real, dimension(:), optional, intent(out)   :: xgmin       !< The x-position of the global minimum in the
                                                             !! units of G%geoLonT, often [degrees_E] or [km]
  real, dimension(:), optional, intent(out)   :: ygmin       !< The y-position of the global minimum in the
                                                             !! units of G%geoLatT, often [degrees_N] or [km]
  real, dimension(:), optional, intent(out)   :: zgmin       !< The z-position of the global minimum [layer]
  real, dimension(:), optional, intent(out)   :: xgmax       !< The x-position of the global maximum in the
                                                             !! units of G%geoLonT, often [degrees_E] or [km]
  real, dimension(:), optional, intent(out)   :: ygmax       !< The y-position of the global maximum in the
                                                             !! units of G%geoLatT, often [degrees_N] or [km]
  real, dimension(:), optional, intent(out)   :: zgmax       !< The z-position of the global maximum [layer]

  ! Local variables
  ! real, dimension(MAX_FIELDS_) :: values ! Globally integrated tracer amounts in a
                                           ! new list for each tracer package [kg conc]
                                                           ! new list for each tracer package [kg conc]
                                                           ! single master list for all tracers [kg conc]

end subroutine call_tracer_stocks
module subroutine store_stocks(pkg_name, ns, names, units, values, index, stock_values, &
                        set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  character(len=*),   intent(in)    :: pkg_name !< The tracer package name
  integer,            intent(in)    :: ns      !< The number of stocks associated with this tracer package
  character(len=*), dimension(:), &
                      intent(in)    :: names   !< Diagnostic names to use for each stock.
  character(len=*), dimension(:), &
                      intent(in)    :: units   !< Units to use in the metadata for each stock.
  type(EFP_type), dimension(:), &
                      intent(in)    :: values  !< The values of the tracer stocks [conc kg]
  integer,            intent(in)    :: index   !< The integer stock index from
                             !! stocks_constants_mod of the stock to be returned.  If this is
                             !! present and greater than 0, only a single stock can be returned.
  type(EFP_type), dimension(:), &
                      intent(inout) :: stock_values !< The master list of stock values [conc kg]
  character(len=*),   intent(inout) :: set_pkg_name !< The name of the last tracer package whose
                                               !! stocks were stored for a specific index.  This is
                                               !! used to trigger an error if there are redundant stocks.
  integer,            intent(in)    :: max_ns  !< The maximum size of the master stock list
  integer,            intent(inout) :: ns_tot  !< The total number of stocks in the master list
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_names !< Diagnostic names to use for each stock in the master list
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_units !< Units to use in the metadata for each stock in the master list

! This routine stores the stocks and does error handling for call_tracer_stocks.

end subroutine store_stocks
module subroutine call_tracer_surface_state(sfc_state, h, G, GV, US, CS)
  type(surface),                intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(ocean_grid_type),        intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),      intent(in)    :: GV    !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(unit_scale_type),        intent(in)    :: US    !< A dimensional unit scaling type
  type(tracer_flow_control_CS), pointer       :: CS    !< The control structure returned by a
                                                       !! previous call to call_tracer_register.

end subroutine call_tracer_surface_state
module subroutine tracer_flow_control_end(CS)
  type(tracer_flow_control_CS), pointer :: CS    !< The control structure returned by a
                                                 !! previous call to call_tracer_register.

end subroutine tracer_flow_control_end
  end interface

end module MOM_tracer_flow_control
