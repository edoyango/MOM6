! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A tracer package that mimics salinity
module pseudo_salt_tracer

use MOM_coms,            only : EFP_type
use MOM_debugging,       only : hchksum
use MOM_diag_mediator,   only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_CVMix_KPP,       only : KPP_NonLocalTransport, KPP_CS
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type, tracer_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_pseudo_salt_tracer, initialize_pseudo_salt_tracer
public pseudo_salt_tracer_column_physics, pseudo_salt_tracer_surface_state
public pseudo_salt_stock, pseudo_salt_tracer_end

!> The control structure for the pseudo-salt tracer
type, public :: pseudo_salt_tracer_CS ; private
  type(tracer_type), pointer :: tr_ptr !< pointer to tracer inside Tr_reg
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: ps(:,:,:) => NULL()   !< The array of pseudo-salt tracer used in this
                                         !! subroutine [ppt]
  real, allocatable :: diff(:,:,:)       !< The difference between the pseudo-salt
                                         !! tracer and the real salt [ppt].
  logical :: pseudo_salt_may_reinit = .true. !< Hard coding since this should not matter

  integer :: id_psd = -1                 !< A diagnostic ID

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                                         !! the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc) :: tr_desc !< A description and metadata for the pseudo-salt tracer
end type pseudo_salt_tracer_CS


  interface
module function register_pseudo_salt_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(pseudo_salt_tracer_CS),  pointer  :: CS   !< The control structure returned by a previous
                                                 !! call to register_pseudo_salt_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and
                                                 !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables
  ! This include declares and sets the variable "version".
  logical :: register_pseudo_salt_tracer
end function register_pseudo_salt_tracer
module subroutine initialize_pseudo_salt_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(pseudo_salt_tracer_CS),        pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_pseudo_salt_tracer
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure containing various thermodynamic variables

  !   This subroutine initializes the tracer fields in CS%ps(:,:,:).

  ! Local variables

end subroutine initialize_pseudo_salt_tracer
module subroutine pseudo_salt_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, debug, &
              KPP_CSp, nonLocalTrans, evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea    !< The amount of fluid entrained from the layer above
                                               !! during this call [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb    !< The amount of fluid entrained from the layer below
                                               !! during this call [H ~> m or kg m-2]
  type(forcing),           intent(in) :: fluxes !< A structure containing thermodynamic and
                                               !! tracer forcing fields
  real,                    intent(in) :: dt    !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  type(pseudo_salt_tracer_CS), pointer :: CS   !< The control structure returned by a previous
                                               !! call to register_pseudo_salt_tracer
  type(thermo_var_ptrs),   intent(in) :: tv    !< A structure pointing to various thermodynamic variables
  logical,                 intent(in) :: debug !< If true calculate checksums
  type(KPP_CS),  optional, pointer    :: KPP_CSp  !< KPP control structure
  real,          optional, intent(in)   :: nonLocalTrans(:,:,:) !< Non-local transport [nondim]
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                               !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                               !! fluxes can be applied [H ~> m or kg m-2]

  !   This subroutine applies diapycnal diffusion and any other column
  ! tracer physics or chemistry to the tracers from this file.

  ! The arguments to this subroutine are redundant in that
  !     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
                              ! [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
                              ! a timestep [ppt H ~> ppt m or ppt kg m-2]
                              ! away [H ~> m or kg m-2]

end subroutine pseudo_salt_tracer_column_physics
module function pseudo_salt_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV     !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc]
  type(pseudo_salt_tracer_CS),        pointer       :: CS     !< The control structure returned by a previous
                                                              !! call to register_pseudo_salt_tracer
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                              !! being sought
  integer                                           :: pseudo_salt_stock !< Return value: the number of
                                                              !! stocks calculated here


end function pseudo_salt_stock
module subroutine pseudo_salt_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(pseudo_salt_tracer_CS),  pointer  :: CS !< The control structure returned by a previous
                                               !! call to register_pseudo_salt_tracer

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine pseudo_salt_tracer_surface_state
module subroutine pseudo_salt_tracer_end(CS)
  type(pseudo_salt_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                             !! call to register_pseudo_salt_tracer

end subroutine pseudo_salt_tracer_end
  end interface

end module pseudo_salt_tracer
