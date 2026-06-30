! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates various values related to the bottom boundary layer, such as the viscosity and
!! thickness of the BBL (set_viscous_BBL).
module MOM_set_visc

use MOM_ALE,           only : ALE_CS, ALE_remap_velocities, ALE_remap_interface_vals, ALE_remap_vertex_vals
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_cvmix_conv,    only : cvmix_conv_is_used
use MOM_CVMix_ddiff,   only : CVMix_ddiff_is_used
use MOM_cvmix_shear,   only : cvmix_shear_is_used
use MOM_debugging,     only : uvchksum, hchksum
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains,       only : pass_var, CORNER
use MOM_EOS,           only : calculate_density, calculate_density_derivs, calculate_specific_vol_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_forcing_type,  only : forcing, mech_forcing, find_ustar
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_intrinsic_functions, only : cuberoot
use MOM_io,            only : slasher, MOM_read_data, vardesc, var_desc
use MOM_kappa_shear,   only : kappa_shear_is_used, kappa_shear_at_vertex
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_restart,       only : register_restart_field_as_obsolete, register_restart_pair
use MOM_safe_alloc,    only : safe_alloc_ptr, safe_alloc_alloc
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type, porous_barrier_type
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units

implicit none ; private

#include <MOM_memory.h>

public set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_end
public set_visc_register_restarts, set_u_at_v, set_v_at_u
public remap_vertvisc_aux_vars

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for MOM_set_visc
type, public :: set_visc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Hbbl           !< The static bottom boundary layer thickness [H ~> m or kg m-2].
                            !! Runtime parameter `HBBL`.
  real    :: dz_bbl         !< The static bottom boundary layer thickness in height units [Z ~> m].
                            !! Runtime parameter `HBBL`.
  real    :: cdrag          !< The quadratic drag coefficient [nondim].
                            !! Runtime parameter `CDRAG`.
  real    :: c_Smag         !< The Laplacian Smagorinsky coefficient for
                            !! calculating the drag in channels [nondim].
  real    :: drag_bg_vel    !< An assumed unresolved background velocity for
                            !! calculating the bottom drag [L T-1 ~> m s-1].
                            !! Runtime parameter `DRAG_BG_VEL`.
                            !! Should not be used if BBL_USE_TIDAL_BG is True.
  real    :: BBL_thick_min  !< The minimum bottom boundary layer thickness [Z ~> m].
                            !! This might be Kv / (cdrag * drag_bg_vel) to give
                            !! Kv as the minimum near-bottom viscosity.
  real    :: Htbl_shelf     !< A nominal thickness of the surface boundary layer for use
                            !! in calculating the near-surface velocity [H ~> m or kg m-2].
  real    :: Htbl_shelf_min !< The minimum surface boundary layer thickness [Z ~> m].
  real    :: KV_BBL_min     !< The minimum viscosity in the bottom boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: KV_TBL_min     !< The minimum viscosity in the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost `HBBL`, depending
                            !! on whether linear_drag is true.
                            !! Runtime parameter `BOTTOMDRAGLAW`.
  logical :: bottomdragmap  !< If true, apply the spatially varying drag coefficient (cdrag_2d)
                            !! instead of the spatially uniform drag coefficient (cdrag).
  logical :: body_force_drag !< If true, the bottom stress is imposed as an explicit body force
                            !! applied over a fixed distance from the bottom, rather than as an
                            !! implicit calculation based on an enhanced near-bottom viscosity.
  logical :: BBL_use_EOS    !< If true, use the equation of state in determining
                            !! the properties of the bottom boundary layer.
  logical :: linear_drag    !< If true, the drag law is cdrag*`DRAG_BG_VEL`*u.
                            !! Runtime parameter `LINEAR_DRAG`.
  logical :: Channel_drag   !< If true, the drag is exerted directly on each layer
                            !! according to what fraction of the bottom they overlie.
  real    :: Chan_drag_max_vol !< The maximum bottom boundary layer volume within which the
                            !! channel drag is applied, normalized by the full cell area,
                            !! or a negative value to apply no maximum [Z ~> m].
  real    :: channel_break_depth !< When CHANNEL_DRAG is true, the bathymetric depth interpolated
                            !! to the vorticity point is a combination of the harmonic mean of the
                            !! adjacent velocity point depths below this depth [Z ~> m] and the
                            !! arithmetic mean of the adjacent depths above it, to roughly mimic a
                            !! continental shelf break profile.  The internal version of this depth
                            !! uses the same offset (G%Z_ref) as the bathymetry.
  logical :: correct_BBL_bounds !< If true, uses the correct bounds on the BBL thickness and
                            !! viscosity so that the bottom layer feels the intended drag.
  logical :: RiNo_mix       !< If true, use Richardson number dependent mixing.
  logical :: dynamic_viscous_ML !< If true, use a bulk Richardson number criterion to
                            !! determine the mixed layer thickness for viscosity.
  real    :: bulk_Ri_ML     !< The bulk mixed layer used to determine the
                            !! thickness of the viscous mixed layer [nondim]
  real    :: omega          !<   The Earth's rotation rate [T-1 ~> s-1].
  real    :: ustar_min      !< A minimum value of ustar to avoid numerical
                            !! problems [H T-1 ~> m s-1 or kg m-2 s-1].  If the value is
                            !! small enough, this should not affect the solution.
  real    :: TKE_decay      !< The ratio of the natural Ekman depth to the TKE
                            !! decay scale [nondim]
  real    :: omega_frac     !<   When setting the decay scale for turbulence, use this
                            !! fraction of the absolute rotation rate blended with the local
                            !! value of f, as sqrt((1-of)*f^2 + of*4*omega^2) [nondim]
  real    :: tideampfac2    !< A factor to multiply by tideamp to convert to a mean ustar,
                            !! accounts for conversion of amplitude to mean magnitude over
                            !! a time average much longer than the tidal periods and for
                            !! non-commuting conversion of mean tideamp to mean ustar**3 [nondim]
  logical :: concave_trigonometric_L  !< If true, use trigonometric expressions to determine the
                            !! fractional open interface lengths for concave topography.
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the set
                            !! viscosity calculations.  Values below 20190101 recover the answers
                            !! from the end of 2018, while higher values use updated and more robust
                            !! forms of the same expressions.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: BBL_use_tidal_bg !< If true, use a tidal background amplitude for the bottom velocity
                            !! when computing the bottom stress.
  character(len=200) :: inputdir !< The directory for input files.
  type(ocean_OBC_type), pointer :: OBC => NULL() !< Open boundaries control structure
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                            !! regulate the timing of diagnostic output.
  ! Allocatable data arrays
  real, allocatable, dimension(:,:) :: cdrag_u !< The spatially varying quadratic drag coefficient [nondim]
  real, allocatable, dimension(:,:) :: cdrag_v !< The spatially varying quadratic drag coefficient [nondim]
  real, allocatable, dimension(:,:) :: tideamp !< RMS tidal amplitude at h points [L T-1 ~> m s-1]
  ! Diagnostic arrays
  real, allocatable, dimension(:,:) :: bbl_u !< BBL mean U current [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: bbl_v !< BBL mean V current [L T-1 ~> m s-1]
  !>@{ Diagnostics handles
  integer :: id_bbl_thick_u = -1, id_kv_bbl_u = -1, id_bbl_u = -1
  integer :: id_bbl_thick_v = -1, id_kv_bbl_v = -1, id_bbl_v = -1
  integer :: id_Ray_u = -1, id_Ray_v = -1
  integer :: id_nkml_visc_u = -1, id_nkml_visc_v = -1
  !>@}
end type set_visc_CS


  interface
module subroutine set_viscous_BBL(u, v, h, tv, visc, G, GV, US, CS, pbv)
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),    intent(in)    :: tv   !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL pointers.
  type(vertvisc_type),      intent(inout) :: visc !< A structure containing vertical viscosities and
                                                  !! related fields.
  type(set_visc_CS),        intent(inout) :: CS   !< The control structure returned by a previous
                                                  !! call to set_visc_init.
  type(porous_barrier_type),intent(in)    :: pbv  !< porous barrier fractional cell metrics

  ! Local variables
end subroutine set_viscous_BBL
module subroutine find_L_open_uniform_slope(vol_below, Dp, Dm, L, GV)
  type(verticalGrid_type),     intent(in)  :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), intent(in)  :: vol_below !< The volume below each interface, normalized by
                                                   !! the full horizontal area of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dp   !< The larger of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dm   !< The smaller of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real, dimension(SZK_(GV)+1), intent(out) :: L    !< The fraction of the full cell width that is open at
                                                   !! the depth of each interface [nondim]

  ! Local variables

end subroutine find_L_open_uniform_slope
module subroutine find_L_open_concave_trigonometric(vol_below, D_vel, Dp, Dm, L, GV)
  type(verticalGrid_type),     intent(in)  :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), intent(in)  :: vol_below !< The volume below each interface, normalized by
                                                   !! the full horizontal area of a velocity cell [Z ~> m]
  real,                        intent(in)  :: D_vel !< The average bottom depth at a velocity point [Z ~> m]
  real,                        intent(in)  :: Dp   !< The larger of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dm   !< The smaller of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real, dimension(SZK_(GV)+1), intent(out) :: L    !< The fraction of the full cell width that is open at
                                                   !! the depth of each interface [nondim]

  ! Local variables
                           ! cell, times the cell width squared [Z ~> m].
                           ! a cell times the cell width [Z ~> m].
  ! The following "volumes" have units of vertical heights because they are normalized
  ! by the full horizontal area of a velocity cell.
                           ! open areas that must be integrated [Z ~> m].

end subroutine find_L_open_concave_trigonometric
module subroutine find_L_open_concave_iterative(vol_below, D_vel, Dp, Dm, L, GV)
  type(verticalGrid_type),     intent(in)  :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), intent(in)  :: vol_below !< The volume below each interface, normalized by
                                                   !! the full horizontal area of a velocity cell [Z ~> m]
  real,                        intent(in)  :: D_vel !< The average bottom depth at a velocity point [Z ~> m]
  real,                        intent(in)  :: Dp   !< The larger of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dm   !< The smaller of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real, dimension(SZK_(GV)+1), intent(out) :: L    !< The fraction of the full cell width that is open at
                                                   !! the depth of each interface [nondim]

  ! Local variables
                           ! cell, times the cell width squared [Z ~> m].
                           ! a cell times the cell width [Z ~> m].

  ! The following "volumes" have units of vertical heights because they are normalized
  ! by the full horizontal area of a velocity cell.
                           ! open areas that must be integrated [Z ~> m].
                           ! relating L to vol_err when there is a single open region [Z ~> m]
                           ! relating L to vol_err when there are two open regions [Z ~> m]

                           ! relating L to vol_err when there is a single open region [nondim]
                           ! relating L to vol_err when there is are two open regions [nondim]
                           ! of L and the target value [Z ~> m]

  ! The following combinations of slope and crv are reused across layers, and hence are pre-calculated
  ! for efficiency.  All are non-negative.
  ! These are only used if the slope exceeds or matches the curvature.
  ! These are only used if the curvature exceeds the slope.
                           ! divided by the curvature [Z ~> m]


end subroutine find_L_open_concave_iterative
module subroutine test_L_open_concave(vol_below, D_vel, Dp, Dm, L, vol_err, GV)
  type(verticalGrid_type),     intent(in)  :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), intent(in)  :: vol_below !< The volume below each interface, normalized by
                                                   !! the full horizontal area of a velocity cell [Z ~> m]
  real,                        intent(in)  :: D_vel !< The average bottom depth at a velocity point [Z ~> m]
  real,                        intent(in)  :: Dp   !< The larger of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dm   !< The smaller of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real, dimension(SZK_(GV)+1), intent(in)  :: L    !< The fraction of the full cell width that is open at
                                                   !! the depth of each interface [nondim]
  real, dimension(SZK_(GV)+1), intent(out) :: vol_err !< The difference between vol_below and the
                                                   !! value obtained from using L in the cubic equation [Z ~> m]

  ! Local variables
                           ! cell, times the cell width squared [Z ~> m].
                           ! a cell times the cell width [Z ~> m].

  ! The following "volumes" have units of vertical heights because they are normalized
  ! by the full horizontal area of a velocity cell.
                           ! open areas that must be integrated [Z ~> m].

  ! The following combinations of slope and crv are reused across layers, and hence are pre-calculated
  ! for efficiency.  All are non-negative.
  ! These are only used if the curvature exceeds the slope.


end subroutine test_L_open_concave
module subroutine find_L_open_convex(vol_below, D_vel, Dp, Dm, L, GV, US, CS)
  type(verticalGrid_type),     intent(in)  :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), intent(in)  :: vol_below  !< The volume below each interface, normalized by
                                                   !! the full horizontal area of a velocity cell [Z ~> m]
  real,                        intent(in)  :: D_vel !< The average bottom depth at a velocity point [Z ~> m]
  real,                        intent(in)  :: Dp   !< The larger of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real,                        intent(in)  :: Dm   !< The smaller of the two depths at the edge
                                                   !! of a velocity cell [Z ~> m]
  real, dimension(SZK_(GV)+1), intent(out) :: L    !< The fraction of the full cell width that is open at
                                                   !! the depth of each interface [nondim]
  type(unit_scale_type),       intent(in)  :: US   !< A dimensional unit scaling type
  type(set_visc_CS),           intent(in)  :: CS   !< The control structure returned by a previous
                                                   !! call to set_visc_init.

  ! Local variables
                           ! cell, times the cell width squared [Z ~> m].
                           ! a cell times the cell width [Z ~> m].
  ! All of the following "volumes" have units of vertical heights because they are normalized
  ! by the full horizontal area of a velocity cell.
                           ! L, or the error for the interface below [Z ~> m].
                           ! solution of a cubic equation for L.
                           ! evaluated at L=L0 [Z ~> m].
                           ! accuracy of a single L(:) Newton iteration [Z5 ~> m5]

end subroutine find_L_open_convex
module function set_v_at_u(v, h, G, GV, i, j, k, mask2dCv, OBC)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer,                 intent(in) :: i    !< The i-index of the u-location to work on.
  integer,                 intent(in) :: j    !< The j-index of the u-location to work on.
  integer,                 intent(in) :: k    !< The k-index of the u-location to work on.
  real, dimension(SZI_(G),SZJB_(G)),&
                           intent(in) :: mask2dCv !< A multiplicative mask of the v-points [nondim]
  type(ocean_OBC_type),    pointer    :: OBC  !< A pointer to an open boundary condition structure
  real                                :: set_v_at_u !< The return value of v at u points points in the
                                              !! same units as u, i.e. [L T-1 ~> m s-1] or other units.

  ! This subroutine finds a thickness-weighted value of v at the u-points.

end function set_v_at_u
module function set_u_at_v(u, h, G, GV, i, j, k, mask2dCu, OBC)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity [L T-1 ~> m s-1] or other units.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer,                 intent(in) :: i    !< The i-index of the u-location to work on.
  integer,                 intent(in) :: j    !< The j-index of the u-location to work on.
  integer,                 intent(in) :: k    !< The k-index of the u-location to work on.
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in) :: mask2dCu !< A multiplicative mask of the u-points [nondim]
  type(ocean_OBC_type),    pointer    :: OBC  !< A pointer to an open boundary condition structure
  real                                :: set_u_at_v !< The return value of u at v points in the
                                              !! same units as u, i.e. [L T-1 ~> m s-1] or other units.

  ! This subroutine finds a thickness-weighted value of u at the v-points.

end function set_u_at_v
module subroutine set_viscous_ML(u, v, h, tv, forces, visc, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure containing pointers to any available
                                                 !! thermodynamic fields. Absent fields have
                                                 !! NULL pointers.
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(inout) :: visc !< A structure containing vertical viscosities and
                                                 !! related fields.
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(set_visc_CS),       intent(inout) :: CS   !< The control structure returned by a previous
                                                 !! call to set_visc_init.

  ! Local variables
end subroutine set_viscous_ML
module subroutine set_visc_register_restarts(HI, G, GV, US, param_file, visc, restart_CS, use_ice_shelf)
  type(hor_index_type),    intent(in)    :: HI         !< A horizontal index type structure.
  type(ocean_grid_type),   intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                       !! parameters.
  type(vertvisc_type),     intent(inout) :: visc       !< A structure containing vertical
                                                       !! viscosities and related fields.
                                                       !! Allocated here.
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure
  logical,                 intent(in) :: use_ice_shelf !< if true, register tau_shelf restarts
  ! Local variables
end subroutine set_visc_register_restarts
module subroutine remap_vertvisc_aux_vars(G, GV, visc, h_old, h_new, ALE_CSp, OBC)
  type(ocean_grid_type),            intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV       !< ocean vertical grid structure
  type(vertvisc_type),              intent(inout) :: visc     !< A structure containing vertical
                                                              !! viscosities and related fields.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_old    !< Thickness of source grid  [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_new    !< Thickness of destination grid [H ~> m or kg m-2]
  type(ALE_CS),                     pointer       :: ALE_CSp  !< ALE control structure to use when remapping
  type(ocean_OBC_type),             pointer       :: OBC      !< Open boundary structure

end subroutine remap_vertvisc_aux_vars
module subroutine set_visc_init(Time, G, GV, US, param_file, diag, visc, CS, restart_CS, OBC)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(vertvisc_type),     intent(inout) :: visc !< A structure containing vertical viscosities and
                                                 !! related fields.
  type(set_visc_CS),       intent(inout) :: CS   !< Vertical viscosity control structure
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure
  type(ocean_OBC_type),    pointer       :: OBC  !< A pointer to an open boundary condition structure

  ! Local variables
                             ! rate of TKE [nondim]
                             ! is used in place of the absolute value of the local Coriolis
                             ! parameter in the denominator of some expressions [nondim]
                             ! to the vorticity point is a combination of the harmonic mean of the
                             ! adjacent velocity point depths below this depth [Z ~> m] and the
                             ! arithmetic mean of the adjacent depths above it, to roughly mimic a
                             ! continental shelf break profile.

                             ! isopycnal or stacked shallow water mode.
  ! This include declares and sets the variable "version".

end subroutine set_visc_init
module subroutine set_visc_end(visc, CS)
  type(vertvisc_type), intent(inout) :: visc !< A structure containing vertical viscosities and
                                             !! related fields.  Elements are deallocated here.
  type(set_visc_CS),   intent(inout) :: CS   !< The control structure returned by a previous
                                             !! call to set_visc_init.

end subroutine set_visc_end
  end interface

end module MOM_set_visc
