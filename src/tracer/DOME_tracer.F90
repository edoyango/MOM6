! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A tracer package that is used as a diagnostic in the DOME experiments
module DOME_tracer

use MOM_coupler_types,   only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_hor_index,       only : hor_index_type
use MOM_grid,            only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_io,              only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type, OBC_segment_tracer_type
use MOM_open_boundary,   only : OBC_segment_type
use MOM_restart,         only : MOM_restart_CS
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_DOME_tracer, initialize_DOME_tracer
public DOME_tracer_column_physics, DOME_tracer_surface_state, DOME_tracer_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

integer, parameter :: ntr = 11 !< The number of tracers in this module.

!> The DOME_tracer control structure
type, public :: DOME_tracer_CS ; private
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, perhaps in [g kg-1]
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out, perhaps in [g kg-1]
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  real :: stripe_width  !< The meridional width of the vertical stripes in the initial condition
                        !! for some of the DOME tracers, in [km] or [degrees_N] or [m].
  real :: stripe_s_lat  !< The southern latitude of the first vertical stripe in the initial condition
                        !! for some of the DOME tracers, in [km] or [degrees_N] or [m].
  real :: sheet_spacing !< The vertical spacing between successive horizontal sheets of tracer in the initial
                        !! conditions for some of the DOME tracers [Z ~> m], and twice the thickness of
                        !! these horizontal tracer sheets

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                    !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers
end type DOME_tracer_CS


  interface
module function register_DOME_tracer(G, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),    intent(in)   :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)   :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)   :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)   :: param_file !< A structure to parse for run-time parameters
  type(DOME_tracer_CS),     pointer      :: CS   !< A pointer that is set to point to the
                                                 !! control structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  ! This include declares and sets the variable "version".
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  logical :: register_DOME_tracer
end function register_DOME_tracer
module subroutine initialize_DOME_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp, tv)
  type(ocean_grid_type),                 intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                 intent(in) :: US   !< A dimensional unit scaling type
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                               !! been read from a restart file.
  type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,               intent(in) :: diag    !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC     !< Structure specifying open boundary options.
  type(DOME_tracer_CS),                  pointer    :: CS      !< The control structure returned by a previous
                                                               !! call to DOME_register_tracer.
  type(sponge_CS),                       pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                     !! for the sponges, if they are in use.
  type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  ! Local variables
                            ! in roundoff and can be neglected [Z ~> m]

end subroutine initialize_DOME_tracer
module subroutine DOME_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
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
  type(DOME_tracer_CS),    pointer    :: CS   !< The control structure returned by a previous
                                              !! call to DOME_register_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! Local variables
end subroutine DOME_tracer_column_physics
module subroutine DOME_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(DOME_tracer_CS),    pointer       :: CS !< The control structure returned by a previous
                                               !! call to DOME_register_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine DOME_tracer_surface_state
module subroutine DOME_tracer_end(CS)
  type(DOME_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                      !! call to DOME_register_tracer.
end subroutine DOME_tracer_end
  end interface

end module DOME_tracer
