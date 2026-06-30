! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Ideal tracers designed to help diagnose a tracer diffusivity tensor in NeverWorld2
module nw2_tracers

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc
use MOM_restart, only : query_initialized, set_initialized, MOM_restart_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_nw2_tracers
public initialize_nw2_tracers
public nw2_tracer_column_physics
public nw2_tracers_end

!> The control structure for the nw2_tracers package
type, public :: nw2_tracers_CS ; private
  integer :: ntr = 0  !< The number of tracers that are actually used.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in [conc] (g m-3)?
  real, allocatable , dimension(:) :: restore_rate !< The rate at which the tracer is damped toward
                                             !! its target profile [T-1 ~> s-1]
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart controls structure
end type nw2_tracers_CS


  interface
logical module function register_nw2_tracers(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(nw2_tracers_CS),       pointer    :: CS !< The control structure returned by a previous
                                               !! call to register_nw2_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! This include declares and sets the variable "version".
end function register_nw2_tracers
module subroutine initialize_nw2_tracers(restart, day, G, GV, US, h, tv, diag, CS)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various
                                                         !! thermodynamic variables
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(nw2_tracers_CS),               pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_nw2_tracer.
  ! Local variables

end subroutine initialize_nw2_tracers
module subroutine nw2_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, tv, CS, &
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
  type(nw2_tracers_CS),    pointer    :: CS   !< The control structure returned by a previous
                                              !! call to register_nw2_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
  ! Local variables

! if (.not.associated(CS)) return

end subroutine nw2_tracer_column_physics
real module function nw2_tracer_dist(m, G, GV, eta, i, j, k)
  integer, intent(in) :: m !< Indicates the NW2 tracer
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in) :: eta !< Interface position [Z ~> m]
  integer, intent(in) :: i !< Cell index i
  integer, intent(in) :: j !< Cell index j
  integer, intent(in) :: k !< Layer index k
  ! Local variables
end function nw2_tracer_dist
module subroutine nw2_tracers_end(CS)
  type(nw2_tracers_CS), pointer :: CS !< The control structure returned by a previous
                                      !! call to register_nw2_tracers.

end subroutine nw2_tracers_end
  end interface

end module nw2_tracers
