! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A tracer package for using dyes to diagnose regional flows.
module regional_dyes

use MOM_coms,               only : EFP_type
use MOM_coupler_types,      only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,      only : diag_ctrl, post_data, register_diag_field
use MOM_error_handler,      only : MOM_error, FATAL, WARNING
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,       only : forcing
use MOM_grid,               only : ocean_grid_type
use MOM_hor_index,          only : hor_index_type
use MOM_interface_heights,  only : thickness_to_dz
use MOM_io,                 only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,      only : ocean_OBC_type
use MOM_restart,            only : query_initialized, MOM_restart_CS
use MOM_spatial_means,      only : global_mass_int_EFP
use MOM_sponge,             only : set_up_sponge_field, sponge_CS
use MOM_time_manager,       only : time_type
use MOM_tracer_registry,    only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic,    only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,      only : tracer_Z_init
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : surface, thermo_var_ptrs
use MOM_verticalGrid,       only : verticalGrid_type
use MOM_tracer_advect_schemes, only : set_tracer_advect_scheme, TracerAdvectionSchemeDoc

implicit none ; private

#include <MOM_memory.h>

public register_dye_tracer, initialize_dye_tracer
public dye_tracer_column_physics, dye_tracer_surface_state
public dye_stock, regional_dyes_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the regional dyes tracer package
type, public :: dye_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  real, allocatable, dimension(:) :: dye_source_minlon !< Minimum longitude of region dye will be
                                                       !! injected, in [m] or [km] or [degrees_E]
  real, allocatable, dimension(:) :: dye_source_maxlon !< Maximum longitude of region dye will be
                                                       !! injected, in [m] or [km] or [degrees_E]
  real, allocatable, dimension(:) :: dye_source_minlat !< Minimum latitude of region dye will be
                                                       !! injected, in [m] or [km] or [degrees_N]
  real, allocatable, dimension(:) :: dye_source_maxlat !< Maximum latitude of region dye will be
                                                       !! injected, in [m] or [km] or [degrees_N]
  real, allocatable, dimension(:) :: dye_source_mindepth !< Minimum depth of region dye will be injected [Z ~> m].
  real, allocatable, dimension(:) :: dye_source_maxdepth !< Maximum depth of region dye will be injected [Z ~> m].
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine [CU ~> conc]

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  integer, allocatable, dimension(:) :: id_tr_dia_diff !< Diagnostic IDs for vertical tracer fluxes (positive up)

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .true. !< If true the tracers may be initialized if not found in a restart file
end type dye_tracer_CS


  interface
module function register_dye_tracer(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dye_tracer_CS),        pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables
  ! This include declares and sets the variable "version".
  logical :: register_dye_tracer

end function register_dye_tracer
module subroutine initialize_dye_tracer(restart, day, G, GV, US, h, diag, OBC, CS, sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(dye_tracer_CS),                pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_dye_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< A pointer to the control structure
                                                         !! for the sponges, if they are in use.
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  ! Local variables

end subroutine initialize_dye_tracer
module subroutine dye_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, tv, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(dye_tracer_CS),     pointer    :: CS   !< The control structure returned by a previous
                                              !! call to register_dye_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

  ! Local variables
                                              !! [conc H T-1 ~> conc m s-1]

end subroutine dye_tracer_column_physics
module function dye_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                            !! tracer, in kg times concentration units [kg conc]
  type(dye_tracer_CS),                pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to register_dye_tracer.
  character(len=*), dimension(:),     intent(out)   :: names !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: dye_stock   !< Return value: the number of stocks
                                                                   !! calculated here.

  ! Local variables

end function dye_stock
module subroutine dye_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(dye_tracer_CS),     pointer       :: CS !< The control structure returned by a previous
                                               !! call to register_dye_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine dye_tracer_surface_state
module subroutine regional_dyes_end(CS)
  type(dye_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                     !! call to register_dye_tracer.

end subroutine regional_dyes_end
  end interface

end module regional_dyes
