! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains routines that implement physical fluxes of tracers (e.g. due
!! to surface fluxes or mixing). These are intended to be called from call_tracer_column_fns
!! in the MOM_tracer_flow_control module.
module MOM_tracer_diabatic

use MOM_grid,             only : ocean_grid_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_forcing_type,     only : forcing
use MOM_error_handler,    only : MOM_error, FATAL, WARNING

implicit none ; private

#include <MOM_memory.h>
public tracer_vertdiff, tracer_vertdiff_Eulerian
public applyTracerBoundaryFluxesInOut


  interface
module subroutine tracer_vertdiff(h_old, ea, eb, dt, tr, G, GV, &
                           sfc_flux, btm_flux, btm_reservoir, sink_rate, convert_flux_in)
  type(ocean_grid_type),                     intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old  !< layer thickness before entrainment
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: ea     !< amount of fluid entrained from the layer
                                                                     !! above [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: eb     !< amount of fluid entrained from the layer
                                                                     !! below [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: tr     !< tracer concentration in concentration units [CU]
  real,                                      intent(in)    :: dt     !< amount of time covered by this call [T ~> s]
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: sfc_flux !< surface flux of the tracer in units of
                                                                     !! [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: btm_flux !< The (negative upward) bottom flux of the
                                                                     !! tracer in [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(inout) :: btm_reservoir !< amount of tracer in a bottom reservoir
                                                                     !! [CU R Z ~> CU kg m-2]
  real,                             optional,intent(in)    :: sink_rate !< rate at which the tracer sinks
                                                                     !! [Z T-1 ~> m s-1]
  logical,                          optional,intent(in)    :: convert_flux_in !< True if the specified sfc_flux needs
                                                                     !! to be integrated in time

  ! local variables
                    !! difference in sinking rates across the layer [H ~> m or kg m-2].
                    !! By construction, 0 <= h_minus_dsink < h_work.
                    !! interfaces, limited to prevent characteristics from
                    !! crossing within a single timestep [H ~> m or kg m-2].
                    !! ensure positive definiteness [H ~> m or kg m-2].
                    !! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine tracer_vertdiff
module subroutine tracer_vertdiff_Eulerian(h_old, ent, dt, tr, G, GV, &
                                    sfc_flux, btm_flux, btm_reservoir, sink_rate, convert_flux_in)
  type(ocean_grid_type),                     intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old  !< layer thickness before entrainment
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: ent    !< Amount of fluid mixed across interfaces
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: tr     !< tracer concentration in concentration units [CU]
  real,                                      intent(in)    :: dt     !< amount of time covered by this call [T ~> s]
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: sfc_flux !< surface flux of the tracer in units of
                                                                     !! [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: btm_flux !< The (negative upward) bottom flux of the
                                                                     !! tracer in [CU kg m-2 T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(inout) :: btm_reservoir !< amount of tracer in a bottom reservoir
                                                                     !! [CU R Z ~> CU kg m-2]
  real,                             optional,intent(in)    :: sink_rate !< rate at which the tracer sinks
                                                                     !! [Z T-1 ~> m s-1]
  logical,                          optional,intent(in)    :: convert_flux_in !< True if the specified sfc_flux needs
                                                                     !! to be integrated in time

  ! local variables
                    !! difference in sinking rates across the layer [H ~> m or kg m-2].
                    !! By construction, 0 <= h_minus_dsink < h_work.
                    !! interfaces, limited to prevent characteristics from
                    !! crossing within a single timestep [H ~> m or kg m-2].
                    !! ensure positive definiteness [H ~> m or kg m-2].
                    !! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine tracer_vertdiff_Eulerian
module subroutine applyTracerBoundaryFluxesInOut(G, GV, Tr, dt, fluxes, h, evap_CFL_limit, minimum_forcing_depth, &
               in_flux_optional, out_flux_optional, update_h_opt)

  type(ocean_grid_type),                      intent(in   ) :: G  !< Grid structure
  type(verticalGrid_type),                    intent(in   ) :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: Tr !< Tracer concentration on T-cell [conc]
  real,                                       intent(in   ) :: dt !< Time-step over which forcing is applied [T ~> s]
  type(forcing),                              intent(in   ) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real,                                       intent(in   ) :: evap_CFL_limit !< Limit on the fraction of the
                                                                  !! water that can be fluxed out of the top
                                                                  !! layer in a timestep [nondim]
  real,                                       intent(in   ) :: minimum_forcing_depth !< The smallest depth over
                                                                  !! which fluxes can be applied [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in   ) :: in_flux_optional !< The total time-integrated
                                                                  !! amount of tracer that enters with freshwater
                                                                  !! [conc H ~> conc m or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: out_flux_optional !< The total time-integrated
                                                                  !! amount of tracer that leaves with freshwater
                                                                  !! [conc H ~> conc m or conc kg m-2]
  logical,                          optional, intent(in) :: update_h_opt  !< Optional flag to determine whether
                                                                  !! h should be updated


                                    ! enters with freshwater [conc H ~> conc m or conc kg m-2]
                                    ! leaves with freshwater [conc H ~> conc m or conc kg m-2]
                               ! the freshwater [conc H ~> conc m or conc kg m-2]
                               ! the freshwater [conc H ~> conc m or conc kg m-2]
                               ! supplied from a column that grounded out [H ~> m or kg m-2]

end subroutine applyTracerBoundaryFluxesInOut
  end interface

end module MOM_tracer_diabatic
