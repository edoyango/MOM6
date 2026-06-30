! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines used to set up and use a set of (one for now)
!! dynamically passive tracers in the ISOMIP configuration.
!!
!! For now, just one passive tracer is injected in
!! the sponge layer.
module ISOMIP_tracer

! Original sample tracer package by Robert Hallberg, 2002
! Adapted to the ISOMIP test case by Gustavo Marques, May 2016

use MOM_coms, only : max_across_PEs
use MOM_coupler_types, only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_ISOMIP_tracer, initialize_ISOMIP_tracer
public ISOMIP_tracer_column_physics, ISOMIP_tracer_surface_state, ISOMIP_tracer_end

integer, parameter :: ntr = 1 !< ntr is the number of tracers in this module.

!> ISOMIP tracer package control structure
type, public :: ISOMIP_tracer_CS ; private
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len = 200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in [conc] (g m-3)?
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out [conc].
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux
             !< if it is used and the surface tracer concentrations are to be
             !< provided to the coupler.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers in this package
end type ISOMIP_tracer_CS


  interface
module function register_ISOMIP_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),      intent(in) :: HI    !<A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                       !! to parse for model parameter values.
  type(ISOMIP_tracer_CS),     pointer    :: CS !<A pointer that is set to point to the control
                                                       !! structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !<A pointer to the tracer registry.
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control struct

  ! This include declares and sets the variable "version".
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  logical :: register_ISOMIP_tracer
end function register_ISOMIP_tracer
module subroutine initialize_ISOMIP_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    ALE_sponge_CSp)

  type(ocean_grid_type),                 intent(in) :: G !< Grid structure.
  type(verticalGrid_type),               intent(in) :: GV !< The ocean's vertical grid structure.
  logical,                               intent(in) :: restart !< .true. if the fields have already
                                                       !! been read from a restart file.
  type(time_type), target,               intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2].
  type(diag_ctrl),               target, intent(in) :: diag !< A structure that is used to regulate
                                                       !! diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC !< This open boundary condition type specifies
                                                       !! whether, where, and what open boundary conditions
                                                       !! are used. This is not being used for now.
  type(ISOMIP_tracer_CS),                pointer    :: CS !< The control structure returned by a previous call
                                                       !! to ISOMIP_register_tracer.
  type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp !< A pointer to the control structure for
                                                       !! the sponges, if they are in use.  Otherwise this
                                                       !! may be unassociated.

                            ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine initialize_ISOMIP_tracer
module subroutine ISOMIP_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
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
  type(ISOMIP_tracer_CS),  pointer    :: CS !< The control structure returned by a previous
                                              !! call to ISOMIP_register_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
end subroutine ISOMIP_tracer_column_physics
module subroutine ISOMIP_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(ISOMIP_tracer_CS),  pointer       :: CS !< The control structure returned by a previous
                                               !! call to ISOMIP_register_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine ISOMIP_tracer_surface_state
module subroutine ISOMIP_tracer_end(CS)
  type(ISOMIP_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                        !! call to ISOMIP_register_tracer.

end subroutine ISOMIP_tracer_end
  end interface

end module ISOMIP_tracer
