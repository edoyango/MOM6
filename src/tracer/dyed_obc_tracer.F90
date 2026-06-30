! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This tracer package dyes flow through open boundaries
module dyed_obc_tracer

use MOM_coupler_types,      only : atmos_ocn_coupler_flux
use MOM_diag_mediator,      only : diag_ctrl
use MOM_error_handler,      only : MOM_error, FATAL, WARNING
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,       only : forcing
use MOM_hor_index,          only : hor_index_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,      only : ocean_OBC_type
use MOM_restart,            only : MOM_restart_CS
use MOM_restart,            only : query_initialized, set_initialized
use MOM_time_manager,       only : time_type
use MOM_tracer_registry,    only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic,    only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : surface
use MOM_verticalGrid,       only : verticalGrid_type
use MOM_tracer_registry,    only : tracer_type
use MOM_tracer_registry,    only : tracer_name_lookup
use MOM_tracer_advect_schemes, only : set_tracer_advect_scheme, TracerAdvectionSchemeDoc

implicit none ; private

#include <MOM_memory.h>

public register_dyed_obc_tracer, initialize_dyed_obc_tracer
public dyed_obc_tracer_column_physics, dyed_obc_tracer_end

!> The control structure for the dyed_obc tracer package
type, public :: dyed_obc_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this subroutine in [conc]

  logical :: tracers_may_reinit  !< If true, these tracers be set up via the initialization code if
                                 !! they are not found in the restart files.

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
end type dyed_obc_tracer_CS


  interface
module function register_dyed_obc_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dyed_obc_tracer_CS),   pointer    :: CS   !< A pointer that is set to point to the
                                                 !! control structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  ! This include declares and sets the variable "version".
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  logical :: register_dyed_obc_tracer

end function register_dyed_obc_tracer
module subroutine initialize_dyed_obc_tracer(restart, day, G, GV, h, diag, OBC, CS)
  type(ocean_grid_type),                 intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in) :: GV   !< The ocean's vertical grid structure
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                               !! been read from a restart file.
  type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,               intent(in) :: diag    !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC     !< Structure specifying open boundary options.
  type(dyed_obc_tracer_CS),              pointer    :: CS      !< The control structure returned by a previous
                                                               !! call to dyed_obc_register_tracer.

! Local variables
                            ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine initialize_dyed_obc_tracer
module subroutine dyed_obc_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
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
  type(dyed_obc_tracer_CS), pointer   :: CS   !< The control structure returned by a previous
                                              !! call to dyed_obc_register_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! Local variables
end subroutine dyed_obc_tracer_column_physics
module subroutine dyed_obc_tracer_end(CS)
  type(dyed_obc_tracer_CS), pointer :: CS   !< The control structure returned by a previous
                                            !! call to dyed_obc_register_tracer.

end subroutine dyed_obc_tracer_end
  end interface

end module dyed_obc_tracer
