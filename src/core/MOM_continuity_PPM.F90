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

!> Interfaces to the zonal mass flux routines that are implemented in MOM_continuity_PPM_s.
interface

!> Calculates the mass or volume fluxes through the zonal faces, and other related quantities.
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
end subroutine zonal_mass_flux

!> Sets the effective interface thickness associated with the fluxes at each zonal velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
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
end subroutine zonal_flux_thickness

!> Calculates the mass or volume fluxes through the meridional faces, and other related quantities.
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
end subroutine meridional_mass_flux

!> Sets the effective interface thickness associated with the fluxes at each meridional velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
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
end subroutine meridional_flux_thickness

!> Return the maximum ratio of a/b or maxrat.
module function ratio_max(a, b, maxrat) result(ratio)
  real, intent(in) :: a       !< Numerator, in arbitrary units [A]
  real, intent(in) :: b       !< Denominator, in arbitrary units [B]
  real, intent(in) :: maxrat  !< Maximum value of ratio [A B-1]
  real :: ratio               !< Return value [A B-1]
end function ratio_max

end interface

contains

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
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
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  type(cont_loop_bounds_type) :: LB ! A type indicating the loop range for a phase of the updates
  logical :: x_first

  h_min = GV%Angstrom_H

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_continuity_PPM: Module must be initialized before it is used.")

  x_first = (MOD(G%first_direction,2) == 0)

  if (present(visc_rem_u) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "MOM_continuity_PPM: Either both visc_rem_u and visc_rem_v or neither "// &
      "one must be present in call to continuity_PPM.")

  if (x_first) then
    !  First advect zonally, with loop bounds that accomodate the subsequent meridional advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.true.)
    call zonal_edge_thickness(hin, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, hin, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin)

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hmin=h_min)

  else  ! .not. x_first
    !  First advect meridionally, with loop bounds that accomodate the subsequent zonal advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.true., j_stencil=.false.)
    call meridional_edge_thickness(hin, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, hin, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hmin=h_min)
  endif

end subroutine continuity_PPM

!> Finds the thickness fluxes from the continuity solver without actually updating the
!! layer thicknesses.  Because the fluxes in the two directions are calculated based on the
!! input thicknesses, which are not updated between the direcitons, the fluxes returned here
!! are not the same as those that would be returned by a call to continuity.
subroutine continuity_3d_fluxes(u, v, h, uh, vh, dt, G, GV, US, CS, OBC, pbv)
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
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_mass_flux(u, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_mass_flux(v, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV)

end subroutine continuity_3d_fluxes

!> Find the vertical sum of the thickness fluxes from the continuity solver without actually
!! updating the layer thicknesses.  Because the fluxes in the two directions are calculated
!! based on the input thicknesses, which are not updated between the directions, the fluxes
!! returned here are not the same as those that would be returned by a call to continuity.
subroutine continuity_2d_fluxes(u, v, h, uhbt, vhbt, dt, G, GV, US, CS, OBC, pbv)
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
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_BT_mass_flux(u, h, h_W, h_E, uhbt, dt, G, GV, US, CS, OBC, pbv%por_face_areaU)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_BT_mass_flux(v, h, h_S, h_N, vhbt, dt, G, GV, US, CS, OBC, pbv%por_face_areaV)

end subroutine continuity_2d_fluxes

!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, visc_rem_u, visc_rem_v)
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
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_in  !< Input zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_in  !< Input meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uh  !< Volume flux through zonal faces =
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vh  !< Volume flux through meridional faces =
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.
  u_in(:,:,:) = u(:,:,:)
  v_in(:,:,:) = v(:,:,:)

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_mass_flux(u_in, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                       uhbt=uhbt, visc_rem_u=visc_rem_u, u_cor=u)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_mass_flux(v_in, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                            vhbt=vhbt, visc_rem_v=visc_rem_v, v_cor=v)

end subroutine continuity_adjust_vel


!> Updates the thicknesses due to zonal thickness fluxes.
subroutine continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin, hmin)
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

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  if (present(hin)) then
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo ; enddo ; enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_zonal_convergence

!> Updates the thicknesses due to meridional thickness fluxes.
subroutine continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin, hmin)
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

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  if (present(hin)) then
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo ; enddo ; enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_merdional_convergence


!> Set the reconstructed thicknesses at the eastern and western edges of tracer cells.
subroutine zonal_edge_thickness(h_in, h_W, h_E, G, GV, US, CS, OBC, LB_in)
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
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  call cpu_clock_begin(id_clock_reconstruct)

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  if (CS%upwind_1st) then
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh+1
      h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,nz
      call PPM_reconstruction_x(h_in(:,:,k), h_W(:,:,k), h_E(:,:,k), G, LB, &
                                2.0*GV%Angstrom_H, CS%monotonic, CS%simple_2nd, OBC, k)
    enddo
  endif

  call cpu_clock_end(id_clock_reconstruct)

end subroutine zonal_edge_thickness


!> Set the reconstructed thicknesses at the eastern and western edges of tracer cells.
subroutine meridional_edge_thickness(h_in, h_S, h_N, G, GV, US, CS, OBC, LB_in)
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
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  call cpu_clock_begin(id_clock_reconstruct)

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  if (CS%upwind_1st) then
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh-1,jeh+1 ; do i=ish,ieh
      h_S(i,j,k) = h_in(i,j,k) ; h_N(i,j,k) = h_in(i,j,k)
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,nz
      call PPM_reconstruction_y(h_in(:,:,k), h_S(:,:,k), h_N(:,:,k), G, LB, &
                                2.0*GV%Angstrom_H, CS%monotonic, CS%simple_2nd, OBC, k)
    enddo
  endif

  call cpu_clock_end(id_clock_reconstruct)

end subroutine meridional_edge_thickness



!> Calculates the vertically integrated mass or volume fluxes through the zonal faces.
subroutine zonal_BT_mass_flux(u, h_in, h_W, h_E, uhbt, dt, G, GV, US, CS, OBC, por_face_areaU, LB_in)
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
  real :: uh(SZIB_(G))      ! Volume flux through zonal faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: duhdu(SZIB_(G))   ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  logical, dimension(SZIB_(G)) :: do_I
  real :: ones(SZIB_(G))    ! An array of 1's [nondim]
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC, OBC_in_row

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    ish = LB_in%ish ; ieh = LB_in%ieh ; jsh = LB_in%jsh ; jeh = LB_in%jeh ; nz = GV%ke
  else
    ish = G%isc ; ieh = G%iec ; jsh = G%jsc ; jeh = G%jec ; nz = GV%ke
  endif

  ones(:) = 1.0 ; do_I(:) = .true.

  uhbt(:,:) = 0.0
  !$OMP parallel do default(shared) private(uh,duhdu,OBC_in_row)
  do j=jsh,jeh
    ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
    OBC_in_row = .false.
    if (local_specified_BC) then ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
      if (OBC%segment(abs(OBC%segnum_u(I,j)))%specified) OBC_in_row = .true.
    endif ; enddo ; endif
    do k=1,nz
      ! This sets uh and duhdu.
      call zonal_flux_layer(u(:,j,k), h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), uh, duhdu, ones, &
                            dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, por_face_areaU(:,j,k), &
                            CS%h_marg_min, OBC)
      if (OBC_in_row) then ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
        l_seg = abs(OBC%segnum_u(I,j))
        if (OBC%segment(l_seg)%specified) uh(I) = OBC%segment(l_seg)%normal_trans(I,j,k)
      endif ; enddo ; endif

      ! Accumulate the barotropic transport.
      do I=ish-1,ieh
        uhbt(I,j) = uhbt(I,j) + uh(I)
      enddo
    enddo ! k-loop
  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine zonal_BT_mass_flux


!> Evaluates the zonal mass or volume fluxes in a layer.
subroutine zonal_flux_layer(u, h, h_W, h_E, uh, duhdu, visc_rem, dt, G, US, j, &
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
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i
  logical :: local_open_BC

  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  do I=ish-1,ieh ; if (do_I(I)) then
    ! Set new values of uh and duhdu.
    if (u(I) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i) + h_E(i)) - 2.0*h(i)
      uh(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * u(I) * &
          (h_E(i) + CFL * (0.5*(h_W(i) - h_E(i)) + curv_3*(CFL - 1.5)))
      h_marg = h_E(i) + CFL * ((h_W(i) - h_E(i)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1) + h_E(i+1)) - 2.0*h(i+1)
      uh(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * u(I) * &
          (h_W(i+1) + CFL * (0.5*(h_E(i+1)-h_W(i+1)) + curv_3*(CFL - 1.5)))
      h_marg = h_W(i+1) + CFL * ((h_E(i+1)-h_W(i+1)) + 3.0*curv_3*(CFL - 1.0))
    else
      uh(I) = 0.0
      h_marg = 0.5 * (h_W(i+1) + h_E(i))
    endif
    h_marg = max(h_marg, h_marg_min)
    duhdu(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * h_marg * visc_rem(I)
  endif ; enddo

  if (local_open_BC) then
    do I=ish-1,ieh ; if (do_I(I)) then ; if (OBC%segnum_u(I,j) /= 0) then
      if (OBC%segment(abs(OBC%segnum_u(I,j)))%open) then
        if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
          uh(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * u(I) * h(i)
          duhdu(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * max(h(i), h_marg_min) * visc_rem(I)
        else !  OBC_DIRECTION_W
          uh(I) = (G%dy_Cu(I,j) * por_face_areaU(I)) * u(I) * h(i+1)
          duhdu(I) = (G%dy_Cu(I,j)* por_face_areaU(I)) * max(h(i+1), h_marg_min) * visc_rem(I)
        endif
      endif
    endif ; endif ; enddo
  endif
end subroutine zonal_flux_layer

!> Calculates the vertically integrated mass or volume fluxes through the meridional faces.
subroutine meridional_BT_mass_flux(v, h_in, h_S, h_N, vhbt, dt, G, GV, US, CS, OBC, por_face_areaV, LB_in)
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
  real :: vh(SZI_(G))      ! Volume flux through meridional faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real ::  dvhdv(SZI_(G))  ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  logical, dimension(SZI_(G)) :: do_I
  real :: ones(SZI_(G))    ! An array of 1's [nondim]
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC, OBC_in_row

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    ish = LB_in%ish ; ieh = LB_in%ieh ; jsh = LB_in%jsh ; jeh = LB_in%jeh ; nz = GV%ke
  else
    ish = G%isc ; ieh = G%iec ; jsh = G%jsc ; jeh = G%jec ; nz = GV%ke
  endif

  ones(:) = 1.0 ; do_I(:) = .true.

  vhbt(:,:) = 0.0
  !$OMP parallel do default(shared) private(vh,dvhdv,OBC_in_row)
  do J=jsh-1,jeh
    ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
    OBC_in_row = .false.
    if (local_specified_BC) then ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
      if (OBC%segment(abs(OBC%segnum_v(i,J)))%specified) OBC_in_row = .true.
    endif ; enddo ; endif
    do k=1,nz
      ! This sets vh and dvhdv.
      call merid_flux_layer(v(:,J,k), h_in(:,:,k), h_S(:,:,k), h_N(:,:,k), vh, dvhdv, ones, &
                            dt, G, US, J, ish, ieh, do_I, CS%vol_CFL, por_face_areaV(:,:,k), &
                            CS%h_marg_min, OBC)
      if (OBC_in_row) then ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
        l_seg = abs(OBC%segnum_v(i,J))
        if (OBC%segment(l_seg)%specified) vh(i) = OBC%segment(l_seg)%normal_trans(i,J,k)
      endif ; enddo ; endif

      ! Accumulate the barotropic transport.
      do i=ish,ieh
        vhbt(i,J) = vhbt(i,J) + vh(i)
      enddo
    enddo ! k-loop
  enddo ! j-loop

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_BT_mass_flux


!> Evaluates the meridional mass or volume fluxes in a layer.
subroutine merid_flux_layer(v, h, h_S, h_N, vh, dvhdv, visc_rem, dt, G, US, J, &
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
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h, i.e. [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i
  logical :: local_open_BC

  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  do i=ish,ieh ; if (do_I(i)) then
    if (v(i) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i) * dt * G%IdyT(i,j) ; endif
      curv_3 = (h_S(i,j) + h_N(i,j)) - 2.0*h(i,j)
      vh(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * v(i) * ( h_N(i,j) + CFL * &
          (0.5*(h_S(i,j) - h_N(i,j)) + curv_3*(CFL - 1.5)) )
      h_marg = h_N(i,j) + CFL * ((h_S(i,j) - h_N(i,j)) + &
                                  3.0*curv_3*(CFL - 1.0))
    elseif (v(i) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1) + h_N(i,j+1)) - 2.0*h(i,j+1)
      vh(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * v(i) * ( h_S(i,j+1) + CFL * &
          (0.5*(h_N(i,j+1)-h_S(i,j+1)) + curv_3*(CFL - 1.5)) )
      h_marg = h_S(i,j+1) + CFL * ((h_N(i,j+1)-h_S(i,j+1)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      vh(i) = 0.0
      h_marg = 0.5 * (h_S(i,j+1) + h_N(i,j))
    endif
    h_marg = max(h_marg, h_marg_min)
    dvhdv(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * h_marg * visc_rem(i)
  endif ; enddo

  if (local_open_BC) then
    do i=ish,ieh ; if (do_I(i)) then
      if (OBC%segnum_v(i,J) /= 0) then
        if (OBC%segment(abs(OBC%segnum_v(i,J)))%open) then
          if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
            vh(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * v(i) * h(i,j)
            dvhdv(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * max(h(i,j), h_marg_min) * visc_rem(i)
          else
            vh(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * v(i) * h(i,j+1)
            dvhdv(i) = (G%dx_Cv(i,J)*por_face_areaV(i,J)) * max(h(i,j+1), h_marg_min) * visc_rem(i)
          endif
        endif
      endif
    endif ; enddo
  endif
end subroutine merid_flux_layer

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_x(h_in, h_W, h_E, G, LB, h_min, monotonic, simple_2nd, OBC, k)
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
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, parameter :: oneSixth = 1./6.  ! [nondim]
  real :: h_ip1, h_im1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, n, stencil
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif

  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl-stencil < G%isd) .or. (iel+stencil > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",I0,".")') &
               stencil + max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",I0,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (simple_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_W(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_E(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl,jel ; do i=isl-1,iel+1
      if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
      endif
    enddo ; enddo

    if (local_open_BC) then
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_E_or_W) then
          I=segment%HI%IsdB
          do j=segment%HI%jsd,segment%HI%jed
            slp(i+1,j) = 0.0
            slp(i,j) = 0.0
          enddo
        endif
      enddo
    endif

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_W(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
      h_E(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
    enddo ; enddo
  endif

  if (local_open_BC) then
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_E) then
        I=segment%HI%IsdB
        if (associated(segment%h_Reg)) then
          do j=segment%HI%jsd,segment%HI%jed
            h_W(i+1,j) = segment%h_Reg%h_res(i,j,k)
            h_E(i+1,j) = segment%h_Reg%h_res(i,j,k)
            h_W(i,j) = segment%h_Reg%h_res(i,j,k)
            h_E(i,j) = segment%h_Reg%h_res(i,j,k)
          enddo
        else
          do j=segment%HI%jsd,segment%HI%jed
            h_W(i+1,j) = h_in(i,j)
            h_E(i+1,j) = h_in(i,j)
            h_W(i,j) = h_in(i,j)
            h_E(i,j) = h_in(i,j)
          enddo
        endif
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        if (associated(segment%h_Reg)) then
          do j=segment%HI%jsd,segment%HI%jed
            h_W(i,j) = segment%h_Reg%h_res(i,j,k)
            h_E(i,j) = segment%h_Reg%h_res(i,j,k)
            h_W(i+1,j) = segment%h_Reg%h_res(i,j,k)
            h_E(i+1,j) = segment%h_Reg%h_res(i,j,k)
          enddo
        else
          do j=segment%HI%jsd,segment%HI%jed
            h_W(i,j) = h_in(i+1,j)
            h_E(i,j) = h_in(i+1,j)
            h_W(i+1,j) = h_in(i+1,j)
            h_E(i+1,j) = h_in(i+1,j)
          enddo
        endif
      endif
    enddo
  endif

  if (monotonic) then
    call PPM_limit_CW84(h_in, h_W, h_E, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_W, h_E, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_x

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_y(h_in, h_S, h_N, G, LB, h_min, monotonic, simple_2nd, OBC, k)
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
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, parameter :: oneSixth = 1./6.      ! [nondim]
  real :: h_jp1, h_jm1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, n, stencil
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif

  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",I0,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stencil < G%jsd) .or. (jel+stencil > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",I0,".")') &
                 stencil + max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (simple_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_S(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_N(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo

    if (local_open_BC) then
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_N_or_S) then
          J=segment%HI%JsdB
          do i=segment%HI%isd,segment%HI%ied
            slp(i,j+1) = 0.0
            slp(i,j) = 0.0
          enddo
        endif
      enddo
    endif

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_S(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_N(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  if (local_open_BC) then
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        if (associated(segment%h_Reg)) then
          do i=segment%HI%isd,segment%HI%ied
            h_S(i,j+1) = segment%h_Reg%h_res(i,j,k)
            h_N(i,j+1) = segment%h_Reg%h_res(i,j,k)
            h_S(i,j) = segment%h_Reg%h_res(i,j,k)
            h_N(i,j) = segment%h_Reg%h_res(i,j,k)
          enddo
        else
          do i=segment%HI%isd,segment%HI%ied
            h_S(i,j+1) = h_in(i,j)
            h_N(i,j+1) = h_in(i,j)
            h_S(i,j) = h_in(i,j)
            h_N(i,j) = h_in(i,j)
          enddo
        endif
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        if (associated(segment%h_Reg)) then
          do i=segment%HI%isd,segment%HI%ied
            h_S(i,j) = segment%h_Reg%h_res(i,j,k)
            h_N(i,j) = segment%h_Reg%h_res(i,j,k)
            h_S(i,j+1) = segment%h_Reg%h_res(i,j,k)
            h_N(i,j+1) = segment%h_Reg%h_res(i,j,k)
          enddo
        else
          do i=segment%HI%isd,segment%HI%ied
            h_S(i,j) = h_in(i,j+1)
            h_N(i,j) = h_in(i,j+1)
            h_S(i,j+1) = h_in(i,j+1)
            h_N(i,j+1) = h_in(i,j+1)
          enddo
        endif
      endif
    enddo
  endif

  if (monotonic) then
    call PPM_limit_CW84(h_in, h_S, h_N, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_S, h_N, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_y

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
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
  real    :: curv  ! The grid-normalized curvature of the three thicknesses  [H ~> m or kg m-2]
  real    :: dh    ! The difference between the edge thicknesses             [H ~> m or kg m-2]
  real    :: scale ! A scaling factor to reduce the curvature of the fit               [nondim]
  integer :: i, j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*((h_L(i,j) + h_R(i,j)) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo

end subroutine PPM_limit_pos

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! according to the monotonic prescription of Colella and Woodward, 1984.
subroutine PPM_limit_CW84(h_in, h_L, h_R, G, iis, iie, jis, jie)
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
  real    :: h_i      ! A copy of the cell-average layer thickness                [H ~> m or kg m-2]
  real    :: RLdiff   ! The difference between the input edge values              [H ~> m or kg m-2]
  real    :: RLdiff2  ! The squared difference between the input edge values   [H2 ~> m2 or kg2 m-4]
  real    :: RLmean   ! The average of the input edge thicknesses                 [H ~> m or kg m-2]
  real    :: FunFac   ! A curious product of the thickness slope and curvature [H2 ~> m2 or kg2 m-4]
  integer :: i, j

  do j=jis,jie ; do i=iis,iie
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j)
    if ( ( h_R(i,j) - h_i ) * ( h_i - h_L(i,j) ) <= 0. ) then
      h_L(i,j) = h_i ; h_R(i,j) = h_i
    else
      RLdiff = h_R(i,j) - h_L(i,j)            ! Difference of edge values
      RLmean = 0.5 * ( h_R(i,j) + h_L(i,j) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_L(i,j) = 3. * h_i - 2. * h_R(i,j)
      if ( FunFac < -RLdiff2 ) h_R(i,j) = 3. * h_i - 2. * h_L(i,j)
    endif
  enddo ; enddo

  return
end subroutine PPM_limit_CW84

!> Initializes continuity_ppm_cs
subroutine continuity_PPM_init(Time, G, GV, US, param_file, diag, CS, OBC)
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
  logical :: local_open_BC, use_h_marg_min
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: n

  !> This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_continuity_PPM" ! This module's name.
  character(len=256) :: mesg

  CS%initialized = .true.

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif

! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MONOTONIC_CONTINUITY", CS%monotonic, &
                 "If true, CONTINUITY_PPM uses the Colella and Woodward "//&
                 "monotonic limiter.  The default (false) is to use a "//&
                 "simple positive definite limiter.", default=.false.)
  call get_param(param_file, mdl, "SIMPLE_2ND_PPM_CONTINUITY", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "UPWIND_1ST_CONTINUITY", CS%upwind_1st, &
                 "If true, CONTINUITY_PPM becomes a 1st-order upwind "//&
                 "continuity solver.  This scheme is highly diffusive "//&
                 "but may be useful for debugging or in single-column "//&
                 "mode where its minimal stencil is useful.", default=.false.)
  call get_param(param_file, mdl, "ETA_TOLERANCE", CS%tol_eta, &
                 "The tolerance for the differences between the "//&
                 "barotropic and baroclinic estimates of the sea surface "//&
                 "height due to the fluxes through each face.  The total "//&
                 "tolerance for SSH is 4 times this value.  The default "//&
                 "is 0.5*NK*ANGSTROM, and this should not be set less "//&
                 "than about 10^-15*MAXIMUM_DEPTH.", units="m", scale=GV%m_to_H, &
                 default=0.5*GV%ke*GV%Angstrom_m)

  call get_param(param_file, mdl, "VELOCITY_TOLERANCE", CS%tol_vel, &
                 "The tolerance for barotropic velocity discrepancies "//&
                 "between the barotropic solution and  the sum of the "//&
                 "layer thicknesses.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
                 ! The speed of light is the default.

  call get_param(param_file, mdl, "CONT_PPM_AGGRESS_ADJUST", CS%aggress_adjust,&
                 "If true, allow the adjusted velocities to have a "//&
                 "relative CFL change up to 0.5.", default=.false.)
  CS%vol_CFL = CS%aggress_adjust
  call get_param(param_file, mdl, "CONT_PPM_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers.  The "//&
                 "default is set by CONT_PPM_AGGRESS_ADJUST.", &
                 default=CS%aggress_adjust, do_not_read=CS%aggress_adjust)
  call get_param(param_file, mdl, "CONTINUITY_CFL_LIMIT", CS%CFL_limit_adjust, &
                 "The maximum CFL of the adjusted velocities.", units="nondim", &
                 default=0.5)
  call get_param(param_file, mdl, "CONT_PPM_BETTER_ITER", CS%better_iter, &
                 "If true, stop corrective iterations using a velocity "//&
                 "based criterion and only stop if the iteration is "//&
                 "better than all predecessors.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_USE_VISC_REM_MAX", CS%use_visc_rem_max, &
                 "If true, use more appropriate limiting bounds for "//&
                 "corrections in strongly viscous columns.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_MARGINAL_FACE_AREAS", CS%marginal_faces, &
                 "If true, use the marginal face areas from the continuity "//&
                 "solver for use as the weights in the barotropic solver. "//&
                 "Otherwise use the transport averaged areas.", default=.true.)
  call get_param(param_file, mdl, "CONT_USE_H_MARG_MIN", use_h_marg_min, &
                 "If true, the marginal thickness used and returned from continuity "//&
                 "is bounded from below by a sub-roundoff value.  Otherwise the "//&
                 "minimum is a large negative value to recreate previous answers.", default=.false.)
  CS%diag => diag

  id_clock_reconstruct = cpu_clock_id('(Ocean continuity reconstruction)', grain=CLOCK_ROUTINE)
  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)

  if (use_h_marg_min) then
    CS%h_marg_min = GV%H_subroundoff
  else
    ! The light-year here means there is effectively no floor on negative marginal thicknesses.
    CS%h_marg_min = -1.0e16*GV%m_to_H  ! Negative 1.057 light years
  endif

  if (local_open_BC) then
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (associated(segment%h_Reg)) then
        if (.not. allocated(segment%h_Reg%h_res)) then
          write(mesg,'("In MOM_continuity_PPM, continuity_PPM_init called with ", &
                & "badly configured h_res.")')
          call MOM_error(FATAL, mesg)
        endif
      endif
    enddo
  endif

end subroutine continuity_PPM_init

!> continuity_PPM_stencil returns the continuity solver stencil size
function continuity_PPM_stencil(CS) result(stencil)
  type(continuity_PPM_CS), intent(in) :: CS   !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

  stencil = 3 ; if (CS%simple_2nd) stencil = 2 ; if (CS%upwind_1st) stencil = 1

end function continuity_PPM_stencil

!> Set up a structure that stores the sizes of the i- and j-loops to to work on in the continuity solver.
function set_continuity_loop_bounds(G, CS, i_stencil, j_stencil) result(LB)
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(continuity_PPM_CS), intent(in) :: CS  !< Module's control structure.
  logical,       optional, intent(in) :: i_stencil !< If present and true, extend the i-loop bounds
                                             !! by the stencil width of the continuity scheme.
  logical,       optional, intent(in) :: j_stencil !< If present and true, extend the j-loop bounds
                                             !! by the stencil width of the continuity scheme.
  type(cont_loop_bounds_type)         :: LB  !< A type storing the array sizes to work on in the continuity routines.

  ! Local variables
  logical :: add_i_stencil, add_j_stencil ! Local variables set based on i_stencil and j_stensil
  integer :: stencil    ! The continuity solver stencil size with the current continuity scheme.

  add_i_stencil = .false. ; if (present(i_stencil)) add_i_stencil = i_stencil
  add_j_stencil = .false. ; if (present(j_stencil)) add_j_stencil = j_stencil

  stencil = continuity_PPM_stencil(CS)

  if (add_i_stencil) then
    LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil
  else
    LB%ish = G%isc ; LB%ieh = G%iec
  endif

  if (add_j_stencil) then
    LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
  else
    LB%jsh = G%jsc ; LB%jeh = G%jec
  endif

end function set_continuity_loop_bounds

!> \namespace mom_continuity_ppm
!!
!! This module contains the subroutines that advect layer
!! thickness.  The scheme here uses a Piecewise-Parabolic method with
!! a positive definite limiter.

end module MOM_continuity_PPM
