! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the routines used to set up a
!! dynamically passive tracer.
!! Set up and use passive tracers requires the following:
!! (1) register_RGC_tracer
!! (2) apply diffusion, physics/chemistry and advect the tracer

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Elizabeth Yankovsky, June 2019                                  *
!*********+*********+*********+*********+*********+*********+***********

module RGC_tracer

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_restart, only :  MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS, get_ALE_sponge_nz_data
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_open_boundary, only : ocean_OBC_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_RGC_tracer, initialize_RGC_tracer
public RGC_tracer_column_physics, RGC_tracer_end

integer, parameter :: NTR = 1 !< The number of tracers in this module.

!> tracer control structure
type, public :: RGC_tracer_CS ; private
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len = 200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry.
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package [kg kg-1]
  real, pointer :: tr_aux(:,:,:,:) => NULL() !< The masked tracer concentration  [kg kg-1]
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out [kg kg-1]
  real :: CSL              !< The length of the continental shelf (x direction) [km]
  real :: lensponge        !< the length of the sponge layer [km]
  logical :: mask_tracers  !< If true, tracers are masked out in massless layers.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers.
end type RGC_tracer_CS


  interface
module function register_RGC_tracer(G, GV, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !<A structure indicating the open file to parse
                                                 !! for model parameter values.
  type(RGC_tracer_CS),        pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure

  ! This include declares and sets the variable "version".
  logical :: register_RGC_tracer
end function register_RGC_tracer
module subroutine initialize_RGC_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    layer_CSp, sponge_CSp)

  type(ocean_grid_type),   intent(in) :: G   !< Grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure.
  logical,                 intent(in) :: restart !< .true. if the fields have already
                                             !! been read from a restart file.
  type(time_type), target, intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(diag_ctrl), target, intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),    pointer    :: OBC !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used. This is not being used for now.
  type(RGC_tracer_CS),     pointer    :: CS  !< The control structure returned by a previous
                                             !!   call to RGC_register_tracer.
  type(sponge_CS),         pointer    :: layer_CSp  !< A pointer to the control structure
  type(ALE_sponge_CS),     pointer    :: sponge_CSp !< A pointer to the control structure for the
                                             !! sponges, if they are in use.  Otherwise this may be unassociated.

                            ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine initialize_RGC_tracer
module subroutine RGC_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
                              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
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
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to any possible
                                              !! forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s].
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(RGC_tracer_CS),     pointer    :: CS   !< The control structure returned by a previous call.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can be
                                              !! fluxed out of the top layer in a timestep [nondim].
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which fluxes
                                              !! can be applied [H ~> m or kg m-2].

! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]


end subroutine RGC_tracer_column_physics
module subroutine RGC_tracer_end(CS)
  type(RGC_tracer_CS), pointer :: CS !< The control structure returned by a previous call to RGC_register_tracer.

end subroutine RGC_tracer_end
  end interface

end module RGC_tracer
