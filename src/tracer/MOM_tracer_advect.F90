! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!>  This module contains the subroutines that advect tracers along coordinate surfaces.
module MOM_tracer_advect

use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,       only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,   only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,         only : sum_across_PEs, max_across_PEs
use MOM_domains,         only : create_group_pass, do_group_pass, group_pass_type, pass_var
use MOM_error_handler,   only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,   only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_tracer_advect_schemes, only : ADVECT_PLM, ADVECT_PPMH3, ADVECT_PPM
use MOM_tracer_advect_schemes, only : set_tracer_advect_scheme, TracerAdvectionSchemeDoc
implicit none ; private

#include <MOM_memory.h>

public advect_tracer
public tracer_advect_init
public tracer_advect_end

!> Control structure for this module
type, public :: tracer_advect_CS ; private
  real    :: dt                    !< The baroclinic dynamics time step [T ~> s].
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !< timing of diagnostic output.
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: useHuynhStencilBug = .false. !< If true, use the incorrect stencil width.
                                   !! This is provided for compatibility with legacy simuations.
  type(group_pass_type) :: pass_uhr_vhr_t_hprev !< A structure used for group passes
  integer :: default_advect_scheme = -1 !< Determines which reconstruction to use
end type tracer_advect_CS

!>@{ CPU time clocks
integer :: id_clock_advect
integer :: id_clock_pass
integer :: id_clock_sync
!>@}


  interface
module subroutine advect_tracer(h_end, uhtr, vhtr, OBC, dt, G, GV, US, CS, Reg, x_first_in, &
                         vol_prev, max_iter_in, update_vol_prev, uhr_out, vhr_out)
  type(ocean_grid_type),   intent(inout) :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_end !< Layer thickness after advection [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uhtr  !< Accumulated volume or mass flux through the
                                                  !! zonal faces [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vhtr  !< Accumulated volume or mass flux through the
                                                  !! meridional faces [H L2 ~> m3 or kg]
  type(ocean_OBC_type),    pointer       :: OBC   !< specifies whether, where, and what OBCs are used
  real,                    intent(in)    :: dt    !< time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(tracer_advect_CS),  pointer       :: CS    !< control structure for module
  type(tracer_registry_type), pointer    :: Reg   !< pointer to tracer registry
  logical,       optional, intent(in)    :: x_first_in !< If present, indicate whether to update
                                                  !! first in the x- or y-direction.
  ! The remaining optional arguments are only used in offline tracer mode.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: vol_prev !< Cell volume before advection [H L2 ~> m3 or kg].
                                                  !! If update_vol_prev is true, the returned value is
                                                  !! the cell volume after the transport that was done
                                                  !! by this call, and if all the transport could be
                                                  !! accommodated it should be close to h_end*G%areaT.
  integer,       optional, intent(in)    :: max_iter_in !< The maximum number of iterations
  logical,       optional, intent(in)    :: update_vol_prev !< If present and true, update vol_prev to
                                                  !! return its value after the tracer have been updated.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: uhr_out !< Remaining accumulated volume or mass fluxes
                                                  !! through the zonal faces [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: vhr_out !< Remaining accumulated volume or mass fluxes
                                                  !! through the meridional faces [H L2 ~> m3 or kg]

                                       ! can be simply discarded [H L2 ~> m3 or kg].


end subroutine advect_tracer
module subroutine advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, advect_schemes)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: ntr  !< The number of tracers
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhr !< accumulated volume/mass flux through
                                                                  !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: uh_neglect !< A tiny zonal mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJ_(G),SZK_(GV)),      intent(inout) :: domore_u !< If true, there is more advection to be
                                                                  !! done in this u-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  integer, dimension(ntr),                   intent(in)    :: advect_schemes !< list of advection schemes to use


                        ! part of that volume that might be lost
                        ! due to advection out the other side of
                        ! the grid box, both in [H L2 ~> m3 or kg].
                        ! current iteration [H L2 ~> m3 or kg].
                        ! any of the passes [H ~> m or kg m-2].
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
                        ! the value in the cell whose reconstruction is being found [conc]
                        ! is being found and the minimum of the surrounding values [conc]

  ! keep a local copy of the initial values of domore_u, which is to be used when computing ad2d_x
  ! diagnostic at the end of this subroutine.
end subroutine advect_x
module subroutine advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, advect_schemes)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: ntr !< The number of tracers
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhr !< accumulated volume/mass flux through
                                                                  !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vh_neglect !< A tiny meridional mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJB_(G),SZK_(GV)),     intent(inout) :: domore_v !< If true, there is more advection to be
                                                                  !! done in this v-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  integer, dimension(ntr),                   intent(in)    :: advect_schemes !< list of advection schemes to use

                                ! current iteration [H L2 ~> m3 or kg].
                                ! part of that volume that might be lost
                                ! due to advection out the other side of
                                ! the grid box, both in  [H L2 ~> m3 or kg].
                        ! any of the passes [H ~> m or kg m-2].
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
                        ! the value in the cell whose reconstruction is being found [conc]
                        ! is being found and the minimum of the surrounding values [conc]

end subroutine advect_y
module subroutine tracer_advect_init(Time, G, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time        !< current model time
  type(ocean_grid_type),   intent(in)    :: G           !< ocean grid structure
  type(unit_scale_type),   intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file  !< open file to parse for model parameters
  type(diag_ctrl), target, intent(inout) :: diag        !< regulates diagnostic output
  type(tracer_advect_CS),  pointer       :: CS          !< module control structure

  ! This include declares and sets the variable "version".

end subroutine tracer_advect_init
module subroutine tracer_advect_end(CS)
  type(tracer_advect_CS), pointer :: CS  !< module control structure

end subroutine tracer_advect_end
  end interface

end module MOM_tracer_advect
