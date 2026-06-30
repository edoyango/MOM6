! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides routines that do checksums of groups of MOM variables
module MOM_checksum_packages

!   This module provides several routines that do check-sums of groups
! of variables in the various dynamic solver routines.

use MOM_coms, only : min_across_PEs, max_across_PEs, reproducing_sum
use MOM_debugging, only : hchksum, uvchksum
use MOM_error_handler, only : MOM_mesg, is_root_pe
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

public MOM_state_chksum, MOM_thermo_chksum, MOM_accel_chksum
public MOM_state_stats, MOM_surface_chksum

!> Write out checksums of the MOM6 state variables
interface MOM_state_chksum
  module procedure MOM_state_chksum_5arg
  module procedure MOM_state_chksum_3arg
end interface

#include <MOM_memory.h>

!> A type for storing statistica about a variable
type :: stats ; private
  real :: minimum = 1.E34  !< The minimum value [degC] or [ppt] or other units
  real :: maximum = -1.E34 !< The maximum value [degC] or [ppt] or other units
  real :: average = 0.     !< The average value [degC] or [ppt] or other units
end type stats


  interface
module subroutine MOM_state_chksum_5arg(mesg, u, v, h, uh, vh, G, GV, US, haloshift, symmetric, omit_corners, vel_scale)
  character(len=*),                          &
                           intent(in) :: mesg !< A message that appears on the chksum lines.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity [L T-1 ~> m s-1] or other units.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity [L T-1 ~> m s-1] or other units.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: uh   !< Volume flux through zonal faces = u*h*dy
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: vh   !< Volume flux through meridional faces = v*h*dx
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  integer,       optional, intent(in) :: haloshift !< The width of halos to check (default 0).
  logical,       optional, intent(in) :: symmetric !< If true, do checksums on the fully symmetric
                                                   !! computational domain.
  logical,       optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,          optional, intent(in) :: vel_scale !< The scaling factor to convert velocities to [T m L-1 s-1 ~> 1]


  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
end subroutine MOM_state_chksum_5arg
module subroutine MOM_state_chksum_3arg(mesg, u, v, h, G, GV, US, haloshift, symmetric, omit_corners)
  character(len=*),                intent(in) :: mesg !< A message that appears on the chksum lines.
  type(ocean_grid_type),           intent(in) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),         intent(in) :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                   intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1] or [m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                   intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1] or [m s-1]..
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                   intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(unit_scale_type),           intent(in) :: US !< A dimensional unit scaling type, which is
                                                    !! used to rescale u and v if present.
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0).
  logical,               optional, intent(in) :: symmetric !< If true, do checksums on the fully
                                                    !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts


  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
end subroutine MOM_state_chksum_3arg
module subroutine MOM_thermo_chksum(mesg, tv, G, US, haloshift, omit_corners)
  character(len=*),         intent(in) :: mesg !< A message that appears on the chksum lines.
  type(thermo_var_ptrs),    intent(in) :: tv   !< A structure pointing to various
                                               !! thermodynamic variables.
  type(ocean_grid_type),    intent(in) :: G    !< The ocean's grid structure.
  type(unit_scale_type),    intent(in) :: US   !< A dimensional unit scaling type
  integer,        optional, intent(in) :: haloshift !< The width of halos to check (default 0).
  logical,        optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts

end subroutine MOM_thermo_chksum
module subroutine MOM_surface_chksum(mesg, sfc_state, G, US, haloshift, symmetric)
  character(len=*),      intent(in)    :: mesg !< A message that appears on the chksum lines.
  type(surface),         intent(inout) :: sfc_state !< transparent ocean surface state structure
                                               !! shared with the calling routine data in this
                                               !! structure is intent out.
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US    !< A dimensional unit scaling type
  integer,     optional, intent(in)    :: haloshift !< The width of halos to check (default 0).
  logical,     optional, intent(in)    :: symmetric !< If true, do checksums on the fully symmetric
                                               !! computational domain.


end subroutine MOM_surface_chksum
module subroutine MOM_accel_chksum(mesg, CAu, CAv, PFu, PFv, diffu, diffv, G, GV, US, pbce, &
                            u_accel_bt, v_accel_bt, symmetric)
  character(len=*),         intent(in) :: mesg !< A message that appears on the chksum lines.
  type(ocean_grid_type),    intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in) :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in) :: CAu  !< Zonal acceleration due to Coriolis
                                               !! and momentum advection terms [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in) :: CAv  !< Meridional acceleration due to Coriolis
                                               !! and momentum advection terms [L T-2 ~> m s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in) :: PFu  !< Zonal acceleration due to pressure gradients
                                               !! (equal to -dM/dx) [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in) :: PFv  !< Meridional acceleration due to pressure gradients
                                               !! (equal to -dM/dy) [L T-2 ~> m s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in) :: diffu !< Zonal acceleration due to convergence of the
                                                !! along-isopycnal stress tensor [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(in) :: diffv !< Meridional acceleration due to convergence of
                                                !! the along-isopycnal stress tensor [L T-2 ~> m s-2].
  type(unit_scale_type),    intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                  optional, intent(in) :: pbce !< The baroclinic pressure anomaly in each layer
                                               !! due to free surface height anomalies
                                               !! [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                  optional, intent(in) :: u_accel_bt !< The zonal acceleration from terms in the
                                                     !! barotropic solver [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                  optional, intent(in) :: v_accel_bt !< The meridional acceleration from terms in
                                                     !! the barotropic solver [L T-2 ~> m s-2].
  logical,        optional, intent(in) :: symmetric !< If true, do checksums on the fully symmetric
                                                    !! computational domain.


end subroutine MOM_accel_chksum
module subroutine MOM_state_stats(mesg, u, v, h, Temp, Salt, G, GV, US, allowChange, permitDiminishing)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  character(len=*),        intent(in) :: mesg !< A message that appears on the chksum lines.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, pointer, dimension(:,:,:),           &
                           intent(in) :: Temp !< Temperature [C ~> degC].
  real, pointer, dimension(:,:,:),           &
                           intent(in) :: Salt !< Salinity [S ~> ppt].
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  logical,       optional, intent(in) :: allowChange !< do not flag an error
                                                     !! if the statistics change.
  logical,       optional, intent(in) :: permitDiminishing !< do not flag error if the
                                                           !! extrema are diminishing.

  ! Local variables
end subroutine MOM_state_stats
  end interface

end module MOM_checksum_packages
