! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A tracer package of ideal age tracers
module ideal_age_example

use MOM_coms,          only : EFP_type
use MOM_coupler_types, only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means, only : global_mass_int_EFP
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, time_to_real
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_ideal_age_tracer, initialize_ideal_age_tracer
public ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
public ideal_age_stock, ideal_age_example_end
public count_BL_layers

integer, parameter :: NTR_MAX = 4 !< the maximum number of tracers in this module.

!> The control structure for the ideal_age_tracer package
type, public :: ideal_age_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  integer :: nkbl   !< The number of layers in the boundary layer.  The ideal
                    !1 age tracers are reset in the top nkbl layers.
  character(len=200) :: IC_file !< The file in which the age-tracer initial values
                    !! can be found, or an empty string for internal initialization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package [years] or other units
  real, dimension(NTR_MAX) :: IC_val = 0.0    !< The (uniform) initial condition value [years] or other units
  real, dimension(NTR_MAX) :: young_val = 0.0 !< The value assigned to tr at the surface [years] or other units
  real, dimension(NTR_MAX) :: land_val = -1.0 !< The value of tr used where land is masked out [years] or other units
  real, dimension(NTR_MAX) :: growth_rate !< The exponential growth rate for the young value [year-1]
  real, dimension(NTR_MAX) :: tracer_start_year !< The year in which tracers start aging, or at which the
                                              !! surface value equals young_val [years].
  logical :: use_real_BL_depth   !< If true, uses the BL scheme to determine the number of
                                 !! layers above the BL depth instead of the fixed nkbl value.
  integer :: BL_residence_num    !< The tracer number assigned to the BL residence tracer in this module
  logical :: tracers_may_reinit  !< If true, these tracers be set up via the initialization code if
                                 !! they are not found in the restart files.
  logical :: tracer_ages(NTR_MAX) !< Indicates whether each tracer ages.

  integer, dimension(NTR_MAX) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                        !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart controls structure

  type(vardesc) :: tr_desc(NTR_MAX) !< Descriptions and metadata for the tracers

end type ideal_age_tracer_CS


  interface
module function register_ideal_age_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ideal_age_tracer_CS),  pointer    :: CS !< The control structure returned by a previous
                                               !! call to register_ideal_age_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! This include declares and sets the variable "version".
  logical :: register_ideal_age_tracer
end function register_ideal_age_tracer
module subroutine initialize_ideal_age_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
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
  type(ideal_age_tracer_CS),          pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_ideal_age_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.

!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

  ! Local variables

end subroutine initialize_ideal_age_tracer
module subroutine ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth, Hbl)
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
  type(ideal_age_tracer_CS), pointer  :: CS   !< The control structure returned by a previous
                                              !! call to register_ideal_age_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: Hbl !< Boundary layer thickness [H ~> m or kg m-2]

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
  ! Local variables
end subroutine ideal_age_tracer_column_physics
module function ideal_age_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                            !! tracer, in kg times concentration units [kg conc].
  type(ideal_age_tracer_CS),          pointer       :: CS   !< The control structure returned by a previous
                                                            !! call to register_ideal_age_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: ideal_age_stock !< The number of stocks calculated here.

  ! Local variables

end function ideal_age_stock
module subroutine ideal_age_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(ideal_age_tracer_CS), pointer     :: CS !< The control structure returned by a previous
                                               !! call to register_ideal_age_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine ideal_age_tracer_surface_state
module subroutine ideal_age_example_end(CS)
  type(ideal_age_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                           !! call to register_ideal_age_tracer.

end subroutine ideal_age_example_end
module subroutine count_BL_layers(G, GV, h, Hbl, BL_layers)
  type(ocean_grid_type),            intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: Hbl  !< Boundary layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: BL_layers !< Number of model layers in the boundary layer [nondim]

end subroutine count_BL_layers
  end interface

end module ideal_age_example
