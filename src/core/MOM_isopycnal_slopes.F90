! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculations of isoneutral slopes and stratification.
module MOM_isopycnal_slopes

use MOM_debugging,     only : hchksum, uvchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density_derivs, calculate_density_second_derivs, EOS_domain
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S

implicit none ; private

#include <MOM_memory.h>

public calc_isoneutral_slopes, vert_fill_TS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine calc_isoneutral_slopes(G, GV, US, h, e, tv, dt_kappa_smooth, use_stanley, slope_x, slope_y, &
                                  N2_u, N2_v, dzu, dzv, dzSxN, dzSyN, halo, OBC, OBC_N2, &
                                  drdx_u, drdy_v, drdz_u, drdz_v)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: e    !< Interface heights [Z ~> m]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  real,                                        intent(in)    :: dt_kappa_smooth !< A smoothing vertical
                                                                     !! diffusivity times a smoothing
                                                                     !! timescale [H Z ~> m2 or kg m-1]
  logical,                                     intent(in)    :: use_stanley !< turn on stanley param in slope
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_u !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_v !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between v-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzu  !< Z-thickness at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzv  !< Z-thickness at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSxN !< Z-thickness times zonal slope contribution to
                                                                     !! Eady growth rate at u-points. [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSyN !< Z-thickness times meridional slope contrib. to
                                                                     !! Eady growth rate at v-points. [Z T-1 ~> m s-1]
  integer,                           optional, intent(in)    :: halo !< Halo width over which to compute
  type(ocean_OBC_type),              optional, pointer       :: OBC  !< Open boundaries control structure.
  logical,                           optional, intent(in)    :: OBC_N2 !< If present and true, use interior data
                                                                     !! to calculate stratification at open boundary
                                                                     !! condition faces.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: drdx_u !< Zonal density gradient at u
                                                                       !! along surfaces of constant z
                                                                       !! (not along isopycnals or
                                                                       !! model interfaces) [R L-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: drdy_v !< Meridional density gradient at v
                                                                       !! along surfaces of constant z
                                                                       !! (not along isopycnals or
                                                                       !! model interfaces) [R L-1 ~> kg m-4]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: drdz_u !< Vertical density gradient
                                                                       !! at u [R Z-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: drdz_v !< Vertical density gradient
                                                                       !! at v [R Z-1 ~> kg m-4]

  ! Local variables
end subroutine calc_isoneutral_slopes
module subroutine vert_fill_TS(h, T_in, S_in, kappa_dt, T_f, S_f, G, GV, US, halo_here, larger_h_denom)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: T_in !< Input temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: S_in !< Input salinity [S ~> ppt]
  real,                                      intent(in)  :: kappa_dt !< A vertical diffusivity to use for smoothing
                                                                 !! times a smoothing timescale [H Z ~> m2 or kg m-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T_f  !< Filled temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S_f  !< Filled salinity [S ~> ppt]
  integer,                         optional, intent(in)  :: halo_here !< Number of halo points to work on,
                                                                 !! 0 by default
  logical,                         optional, intent(in)  :: larger_h_denom !< Present and true, add a large
                                                                 !! enough minimal thickness in the denominator of
                                                                 !! the flux calculations so that the fluxes are
                                                                 !! never so large as eliminate the transmission
                                                                 !! of information across groups of massless layers.
  ! Local variables
                                   ! between layers in a timestep [H ~> m or kg m-2].
                                   ! completely decoupling groups of layers [H ~> m or kg m-2].
                                   ! Often 0 < h_neglect << h0.
                                   ! added to ensure positive definiteness [H ~> m or kg m-2].

end subroutine vert_fill_TS
  end interface

end module MOM_isopycnal_slopes
