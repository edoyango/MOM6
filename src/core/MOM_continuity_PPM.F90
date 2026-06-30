! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Solve the layer continuity equation using the PPM method for layer fluxes.
module MOM_continuity_PPM

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity_PPM, continuity_PPM_init, continuity_PPM_stencil
public continuity_fluxes, continuity_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_merdional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds

!>@{ CPU time clock IDs
integer :: id_clock_reconstruct, id_clock_update, id_clock_correct
!>@}

!> Control structure for mom_continuity_ppm
type, public :: continuity_PPM_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag !< Diagnostics control structure.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic
                             !! limiter; otherwise use a simple positive
                             !! definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic
                             !! mean) interpolation of the edge values instead
                             !! of the higher order interpolation.
  real :: tol_eta            !< The tolerance for free-surface height
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [H ~> m or kg m-2].
  real :: tol_vel            !< The tolerance for barotropic velocity
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [L T-1 ~> m s-1].
  real :: CFL_limit_adjust   !< The maximum CFL of the adjusted velocities [nondim]
  real :: h_marg_min         !< Negligible floor on h_marg, the marginal thickness
                             !! used to calculate the partial derivative of transports
                             !! with velocities [H ~> m or kg m-2]
  logical :: aggress_adjust  !< If true, allow the adjusted velocities to have a
                             !! relative CFL change up to 0.5.  False by default.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths
                             !! to the tracer cell areas when estimating CFL
                             !! numbers.  Without aggress_adjust, the default is
                             !! false; it is always true with.
  logical :: better_iter     !< If true, stop corrective iterations using a
                             !! velocity-based criterion and only stop if the
                             !! iteration is better than all predecessors.
  logical :: use_visc_rem_max !< If true, use more appropriate limiting bounds
                             !! for corrections in strongly viscous columns.
  logical :: marginal_faces  !< If true, use the marginal face areas from the
                             !! continuity solver for use as the weights in the
                             !! barotropic solver.  Otherwise use the transport
                             !! averaged areas.
end type continuity_PPM_CS

!> A container for loop bounds
type, public :: cont_loop_bounds_type ; private
  !>@{ Loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type cont_loop_bounds_type

!> Finds the thickness fluxes from the continuity solver or their vertical sum without
!! actually updating the layer thicknesses.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes


  interface
module subroutine continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
                          visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont, du_cor, dv_cor)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: hin !< Initial layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Final layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Zonal volume flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Meridional volume flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< pointers to porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: vhbt !< The summed volume flux through meridional faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                             !< The zonal velocities that give uhbt as the depth-integrated transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: v_cor
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(out)   :: dv_cor !< The meridional velocity increments from v that give vhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  ! Local variables

end subroutine continuity_PPM
module subroutine continuity_3d_fluxes(u, v, h, uh, vh, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Thickness fluxes through zonal faces,
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Thickness fluxes through meridional faces,
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables

end subroutine continuity_3d_fluxes
module subroutine continuity_2d_fluxes(u, v, h, uhbt, vhbt, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(out)   :: uhbt !< Vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(out)   :: vhbt !< Vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables

end subroutine continuity_2d_fluxes
module subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, visc_rem_u, visc_rem_v)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u   !< Zonal velocity, which will be adjusted to
                                                !! give uhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v   !< Meridional velocity, which will be adjusted
                                                !! to give vhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uhbt !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: vhbt !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of the meridional momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.

  ! Local variables
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.
end subroutine continuity_adjust_vel
module subroutine continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin, hmin)
  type(ocean_grid_type),       intent(in)    :: G    !< Ocean's grid structure
  type(verticalGrid_type),     intent(in)    :: GV   !< Ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(inout) :: h    !< Final layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                               intent(in)    :: uh   !< Zonal thickness flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(cont_loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)    :: hin  !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]


end subroutine continuity_zonal_convergence
module subroutine continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin, hmin)
  type(ocean_grid_type),       intent(in)    :: G    !< Ocean's grid structure
  type(verticalGrid_type),     intent(in)    :: GV   !< Ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(inout) :: h    !< Final layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                               intent(in)    :: vh   !< Meridional thickness flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(cont_loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)    :: hin  !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]


end subroutine continuity_merdional_convergence
module subroutine zonal_edge_thickness(h_in, h_W, h_E, G, GV, US, CS, OBC, LB_in)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_W  !< Western edge layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_E  !< Eastern edge layer thickness [H ~> m or kg m-2].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.

  ! Local variables

end subroutine zonal_edge_thickness
module subroutine meridional_edge_thickness(h_in, h_S, h_N, G, GV, US, CS, OBC, LB_in)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_S  !< Southern edge layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_N  !< Northern edge layer thickness [H ~> m or kg m-2].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.

  ! Local variables

end subroutine meridional_edge_thickness
module subroutine zonal_mass_flux(u, h_in, h_W, h_E, uh, dt, G, GV, US, CS, OBC, por_face_areaU, &
                           LB_in, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_W !< Western edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_E !< Eastern edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                           intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a barotropic
                     !! acceleration that a layer experiences after viscosity is applied [nondim].
                     !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                     !< The zonal velocities (u with a barotropic correction)
                     !! that give uhbt as the depth-integrated transport [L T-1 ~> m s-1]
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  ! Local variables
                  ! the time step [T-1 ~> s-1].

end subroutine zonal_mass_flux
module subroutine zonal_BT_mass_flux(u, h_in, h_W, h_E, uhbt, dt, G, GV, US, CS, OBC, por_face_areaU, LB_in)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: u    !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_W  !< Western edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_E  !< Eastern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(out) :: uhbt !< The summed volume flux through zonal
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),  intent(in)  :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.

  ! Local variables

end subroutine zonal_BT_mass_flux
module subroutine zonal_flux_layer(u, h, h_W, h_E, uh, duhdu, visc_rem, dt, G, US, j, &
                            ish, ieh, do_I, vol_CFL, por_face_areaU, h_marg_min, OBC)
  type(ocean_grid_type),        intent(in)    :: G        !< Ocean's grid structure.
  real, dimension(SZIB_(G)),    intent(in)    :: u        !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZIB_(G)),    intent(in)    :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G)),     intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)),     intent(in)    :: h_W      !< West edge thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)),     intent(in)    :: h_E      !< East edge thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G)),    intent(inout) :: uh       !< Zonal mass or volume
                                                          !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)),    intent(inout) :: duhdu    !< Partial derivative of uh
                                                          !! with u [H L ~> m2 or kg m-1].
  real,                         intent(in)    :: dt       !< Time increment [T ~> s]
  type(unit_scale_type),        intent(in)    :: US       !< A dimensional unit scaling type
  integer,                      intent(in)    :: j        !< Spatial index.
  integer,                      intent(in)    :: ish      !< Start of index range.
  integer,                      intent(in)    :: ieh      !< End of index range.
  logical, dimension(SZIB_(G)), intent(in)    :: do_I     !< Which i values to work on.
  logical,                      intent(in)    :: vol_CFL  !< If true, rescale the
  real, dimension(SZIB_(G)),    intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
          !! ratio of face areas to the cell areas when estimating the CFL number.
  real,                         intent(in)    :: h_marg_min !< Negligible floor on h_marg [H ~> m or kg m-2]
  type(ocean_OBC_type), optional, pointer     :: OBC !< Open boundaries control structure.
  ! Local variables

end subroutine zonal_flux_layer
module subroutine zonal_flux_thickness(u, h, h_W, h_E, h_u, dt, G, GV, US, LB, vol_CFL, &
                                marginal, OBC, por_face_areaU, visc_rem_u)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h_u !< Effective thickness at zonal faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(cont_loop_bounds_type),               intent(in)    :: LB   !< Loop bounds structure.
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                                   intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in)    :: visc_rem_u
                          !< Both the fraction of the momentum originally in a layer that remains after
                          !! a time-step of viscosity, and the fraction of a time-step's worth of a
                          !! barotropic acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
end subroutine zonal_flux_thickness
module subroutine zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             j, ish, ieh, do_I_in, por_face_areaU, uh_3d, OBC)

  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZK_(GV)),        intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G)),                 intent(in)    :: uhbt !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].

  real, dimension(SZIB_(G)),                 intent(in)    :: du_max_CFL  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: du_min_CFL  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: uh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G)),                 intent(out)   :: du !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  integer,                                   intent(in)    :: j    !< Spatial index.
  integer,                                   intent(in)    :: ish  !< Start of index range.
  integer,                                   intent(in)    :: ieh  !< End of index range.
  logical, dimension(SZIB_(G)),              intent(in)    :: do_I_in     !<
                       !! A logical flag indicating which I values to work on.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                                      intent(in) :: por_face_areaU !< fractional open area of U-faces [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), optional, intent(inout) :: uh_3d !<
                       !! Volume flux through zonal faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type),            optional, pointer       :: OBC !< Open boundaries control structure.
  ! Local variables

end subroutine zonal_flux_adjust
module subroutine set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, j, ish, ieh, do_I, por_face_areaU)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(BT_cont_type),                        intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G)),                 intent(in)    :: uh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: du_max_CFL  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G)),                 intent(in)    :: du_min_CFL  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZIB_(G),SZK_(GV)),        intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G)),                 intent(in)    :: visc_rem_max !< Maximum allowable visc_rem [nondim].
  integer,                                   intent(in)    :: j        !< Spatial index.
  integer,                                   intent(in)    :: ish      !< Start of index range.
  integer,                                   intent(in)    :: ieh      !< End of index range.
  logical, dimension(SZIB_(G)),              intent(in)    :: do_I     !< A logical flag indicating
                       !! which I values to work on.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                                    intent(in) :: por_face_areaU !< fractional open area of U-faces [nondim]
  ! Local variables
end subroutine set_zonal_BT_cont
module subroutine meridional_mass_flux(v, h_in, h_S, h_N, vh, dt, G, GV, US, CS, OBC, por_face_areaV, &
                                LB_in, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v    !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_S  !< South edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_N  !< North edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: vh   !< Volume flux through meridional
                                                                  !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: por_face_areaV !< fractional open area of V-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(in) :: vhbt !< The summed volume flux through meridional
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(in)  :: visc_rem_v !< Both the fraction of the momentum
                                   !! originally in a layer that remains after a time-step of viscosity,
                                   !! and the fraction of a time-step's worth of a barotropic acceleration
                                   !! that a layer experiences after viscosity is applied [nondim].
                                   !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(out) :: v_cor
                                   !< The meridional velocities (v with a barotropic correction)
                                   !! that give vhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(BT_cont_type),               optional, pointer     :: BT_cont !< A structure with elements that describe
                                   !! the effective open face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZJB_(G)), &
                                    optional, intent(out)   :: dv_cor !< The meridional velocity increments from v
                                                                  !! that give vhbt as the depth-integrated
                                                                  !! transports [L T-1 ~> m s-1].

  ! Local variables
                  ! the time step [T-1 ~> s-1].

end subroutine meridional_mass_flux
module subroutine meridional_BT_mass_flux(v, h_in, h_S, h_N, vhbt, dt, G, GV, US, CS, OBC, por_face_areaV, LB_in)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v    !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_S  !< Southern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_N  !< Northern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G)),          intent(out) :: vhbt !< The summed volume flux through meridional
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: por_face_areaV !< fractional open area of V-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.

  ! Local variables

end subroutine meridional_BT_mass_flux
module subroutine merid_flux_layer(v, h, h_S, h_N, vh, dvhdv, visc_rem, dt, G, US, J, &
                            ish, ieh, do_I, vol_CFL, por_face_areaV, h_marg_min, OBC)
  type(ocean_grid_type),        intent(in)    :: G        !< Ocean's grid structure.
  real, dimension(SZI_(G)),     intent(in)    :: v        !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G)),     intent(in)    :: visc_rem !< Both the fraction of the
         !! momentum originally in a layer that remains after a time-step
         !! of viscosity, and the fraction of a time-step's worth of a barotropic
         !! acceleration that a layer experiences after viscosity is applied [nondim].
         !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: h      !< Layer thickness used to calculate fluxes,
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: h_S    !< South edge thickness in the reconstruction
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: h_N    !< North edge thickness in the reconstruction
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G)),     intent(inout) :: vh       !< Meridional mass or volume transport
                                                          !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G)),     intent(inout) :: dvhdv    !< Partial derivative of vh with v
                                                          !! [H L ~> m2 or kg m-1].
  real,                         intent(in)    :: dt       !< Time increment [T ~> s].
  type(unit_scale_type),        intent(in)    :: US       !< A dimensional unit scaling type
  integer,                      intent(in)    :: j        !< Spatial index.
  integer,                      intent(in)    :: ish      !< Start of index range.
  integer,                      intent(in)    :: ieh      !< End of index range.
  logical, dimension(SZI_(G)),  intent(in)    :: do_I     !< Which i values to work on.
  logical,                      intent(in)    :: vol_CFL  !< If true, rescale the
         !! ratio of face areas to the cell areas when estimating the CFL number.
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in) :: por_face_areaV !< fractional open area of V-faces [nondim]
  real,                         intent(in)    :: h_marg_min !< Negligible floor on h_marg [H ~> m or kg m-2]
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables
                 ! with the same units as h, i.e. [H ~> m or kg m-2].

end subroutine merid_flux_layer
module subroutine meridional_flux_thickness(v, h, h_S, h_N, h_v, dt, G, GV, US, LB, vol_CFL, &
                                     marginal, OBC, por_face_areaV, visc_rem_v)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< Layer thickness used to calculate fluxes,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_S  !< South edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_N  !< North edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: h_v !< Effective thickness at meridional faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(cont_loop_bounds_type),               intent(in)    :: LB   !< Loop bounds structure.
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                     intent(in) :: por_face_areaV  !< fractional open area of V-faces [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), optional, intent(in) :: visc_rem_v !< Both the fraction
                          !! of the momentum originally in a layer that remains after a time-step of
                          !! viscosity, and the fraction of a time-step's worth of a barotropic
                          !! acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
                 ! with the same units as h [H ~> m or kg m-2] .
end subroutine meridional_flux_thickness
module subroutine meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, &
                             dv, dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             j, ish, ieh, do_I_in, por_face_areaV, vh_3d, OBC)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),&
                           intent(in)    :: h_S  !< South edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_N  !< North edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in) :: visc_rem
                             !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G)), intent(in)    :: vhbt !< The summed volume flux through meridional faces
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G)), intent(in)    :: dv_max_CFL !< Maximum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G)), intent(in)    :: dv_min_CFL !< Minimum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G)), intent(in)    :: vh_tot_0   !< The summed transport with 0 adjustment
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G)), intent(in)    :: dvhdv_tot_0 !< The partial derivative of dv_err with
                                                  !! dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G)), intent(out)   :: dv   !< The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                     intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),  intent(in)    :: CS   !< This module's control structure.
  integer,                  intent(in)    :: j    !< Spatial index.
  integer,                  intent(in)    :: ish  !< Start of index range.
  integer,                  intent(in)    :: ieh  !< End of index range.
  logical, dimension(SZI_(G)), &
                            intent(in)    :: do_I_in  !< A flag indicating which I values to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                     intent(in) :: por_face_areaV !< fractional open area of V-faces [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                  optional, intent(inout) :: vh_3d !< Volume flux through meridional
                             !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables

end subroutine meridional_flux_adjust
module subroutine set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, j, ish, ieh, do_I, por_face_areaV)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to calculate fluxes,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_S  !< South edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_N  !< North edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  type(BT_cont_type),                        intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZI_(G)),                  intent(in)    :: vh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G)),                  intent(in)    :: dvhdv_tot_0 !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G)),                  intent(in)    :: dv_max_CFL !< Maximum acceptable value
                                                                   !!  of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G)),                  intent(in)    :: dv_min_CFL !< Minimum acceptable value
                                                                   !!  of dv [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZI_(G),SZK_(GV)),         intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G)),                  intent(in)    :: visc_rem_max !< Maximum allowable visc_rem [nondim]
  integer,                                   intent(in)    :: j        !< Spatial index.
  integer,                                   intent(in)    :: ish      !< Start of index range.
  integer,                                   intent(in)    :: ieh      !< End of index range.
  logical, dimension(SZI_(G)),               intent(in)    :: do_I     !< A logical flag indicating
                       !! which I values to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                intent(in) :: por_face_areaV !< fractional open area of V-faces [nondim]
  ! Local variables
end subroutine set_merid_BT_cont
module subroutine PPM_reconstruction_x(h_in, h_W, h_E, G, LB, h_min, monotonic, simple_2nd, OBC, k)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out) :: h_W  !< West edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out) :: h_E  !< East edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  type(cont_loop_bounds_type),       intent(in)  :: LB   !< Active loop bounds structure.
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.
  integer :: k      !< vertical grid index

  ! Local variables with useful mnemonic names.
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]

end subroutine PPM_reconstruction_x
module subroutine PPM_reconstruction_y(h_in, h_S, h_N, G, LB, h_min, monotonic, simple_2nd, OBC, k)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out) :: h_S  !< South edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out) :: h_N  !< North edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  type(cont_loop_bounds_type),       intent(in)  :: LB   !< Active loop bounds structure.
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.
  integer :: k      !< vertical grid index

  ! Local variables with useful mnemonic names.
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]

end subroutine PPM_reconstruction_y
module subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_L !< Left thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_R !< Right thickness in the reconstruction [H ~> m or kg m-2].
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  integer,                           intent(in)  :: iis      !< Start of i index range.
  integer,                           intent(in)  :: iie      !< End of i index range.
  integer,                           intent(in)  :: jis      !< Start of j index range.
  integer,                           intent(in)  :: jie      !< End of j index range.

! Local variables

end subroutine PPM_limit_pos
module subroutine PPM_limit_CW84(h_in, h_L, h_R, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G     !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in  !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_L !< Left thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_R !< Right thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  integer,                           intent(in)  :: iis   !< Start of i index range.
  integer,                           intent(in)  :: iie   !< End of i index range.
  integer,                           intent(in)  :: jis   !< Start of j index range.
  integer,                           intent(in)  :: jie   !< End of j index range.

  ! Local variables

end subroutine PPM_limit_CW84
module function ratio_max(a, b, maxrat) result(ratio)
  real, intent(in) :: a       !< Numerator, in arbitrary units [A]
  real, intent(in) :: b       !< Denominator, in arbitrary units [B]
  real, intent(in) :: maxrat  !< Maximum value of ratio [A B-1]
  real :: ratio               !< Return value [A B-1]

end function ratio_max
module subroutine continuity_PPM_init(Time, G, GV, US, param_file, diag, CS, OBC)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating
                  !! the open file to parse for model parameter values.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to
                  !! regulate diagnostic output.
  type(continuity_PPM_CS), intent(inout) :: CS   !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.

  !> This include declares and sets the variable "version".

end subroutine continuity_PPM_init
module function continuity_PPM_stencil(CS) result(stencil)
  type(continuity_PPM_CS), intent(in) :: CS   !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

end function continuity_PPM_stencil
module function set_continuity_loop_bounds(G, CS, i_stencil, j_stencil) result(LB)
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(continuity_PPM_CS), intent(in) :: CS  !< Module's control structure.
  logical,       optional, intent(in) :: i_stencil !< If present and true, extend the i-loop bounds
                                             !! by the stencil width of the continuity scheme.
  logical,       optional, intent(in) :: j_stencil !< If present and true, extend the j-loop bounds
                                             !! by the stencil width of the continuity scheme.
  type(cont_loop_bounds_type)         :: LB  !< A type storing the array sizes to work on in the continuity routines.

  ! Local variables

end function set_continuity_loop_bounds
  end interface

end module MOM_continuity_PPM
