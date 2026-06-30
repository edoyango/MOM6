! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A sample tracer package that has striped initial conditions
module USER_tracer_example

use MOM_coms,            only : EFP_type
use MOM_coupler_types,   only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : file_exists, MOM_read_data, slasher
use MOM_io,              only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public USER_register_tracer_example, USER_initialize_tracer, USER_tracer_stock
public tracer_column_physics, USER_tracer_surface_state, USER_tracer_example_end

integer, parameter :: NTR = 1 !< The number of tracers in this module.

!> The control structure for the USER_tracer_example module
type, public :: USER_tracer_example_CS ; private
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " "
                                       !! to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()  !< The array of tracers used in this subroutine, perhaps in [g kg-1]?
  real :: land_val(NTR) = -1.0 !< The value of tr that is used where land is masked out, perhaps in [g kg-1]?

  real :: stripe_width  !< The Gaussian width of the stripe in the initial condition
                        !! for the tracer_example tracers [L ~> m]
  real :: stripe_lat    !< The central latitude of the stripe in the initial condition
                        !! for the tracer_example tracers, in [degrees_N] or [km] or [m].
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                    !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions of each of the tracers.
end type USER_tracer_example_CS


  interface
module function USER_register_tracer_example(G, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),   intent(in)   :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)   :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)   :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)   :: param_file !< A structure to parse for run-time parameters
  type(USER_tracer_example_CS), pointer :: CS   !< A pointer that is set to point to the control
                                                !! structure for this module
  type(tracer_registry_type), pointer   :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS),   intent(inout) :: restart_CS !< MOM restart control struct

! Local variables
  ! This include declares and sets the variable "version".
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  logical :: USER_register_tracer_example
end function USER_register_tracer_example
module subroutine USER_initialize_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(USER_tracer_example_CS),       pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to USER_register_tracer_example.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables

end subroutine USER_initialize_tracer
module subroutine tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS)
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
  type(USER_tracer_example_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to USER_register_tracer_example.

! Local variables
                               ! with surface mass fluxes added back [H ~> m or kg m-2].
                               ! in roundoff and can be neglected [H ~> m or kg m-2].
                               ! advection of the tracer [nondim]

  ! These are the settings for most "physical" tracers, which
  ! are advected diapycnally in the usual manner.
end subroutine tracer_column_physics
module function USER_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc]
  type(USER_tracer_example_CS),       pointer       :: CS     !< The control structure returned by a
                                                              !! previous call to register_USER_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                              !! being sought.
  integer                                           :: USER_tracer_stock !< Return value: the number of
                                                              !! stocks calculated here.

  ! Local variables

end function USER_tracer_stock
module subroutine USER_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),        intent(in)    :: G     !< The ocean's grid structure
  type(verticalGrid_type),      intent(in)    :: GV    !< The ocean's vertical grid structure
  type(surface),                intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(USER_tracer_example_CS), pointer       :: CS !< The control structure returned by a previous
                                                    !! call to register_USER_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine USER_tracer_surface_state
module subroutine USER_tracer_example_end(CS)
  type(USER_tracer_example_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_USER_tracer.

end subroutine USER_tracer_example_end
  end interface

end module USER_tracer_example
