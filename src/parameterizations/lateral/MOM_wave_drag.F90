! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Frequency-dependent linear wave drag

module MOM_wave_drag

use MOM_domains,       only : pass_vector, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, NOTE
use MOM_file_parser,   only : get_param, log_param, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : MOM_read_data, slasher, EAST_FACE, NORTH_FACE
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public wave_drag_init, wave_drag_calc

#include <MOM_memory.h>

!> Control structure for the MOM_wave_drag module
type, public :: wave_drag_CS ; private
  integer :: nf                                 !< Number of filters to be used in the simulation
  real, allocatable, dimension(:,:,:) :: coef_u !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: coef_v !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: coef_uv !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: coef_vu !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
  logical :: tensor_drag                        !< If true, include the off-diagonal components of the
                                                !! wave drag tensor for computing the wave drag
end type wave_drag_CS


  interface
module subroutine wave_drag_init(param_file, wave_drag_file, G, GV, US, CS)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),        intent(in)    :: wave_drag_file !< The file from which to read drag coefficients
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(wave_drag_CS),      intent(out)   :: CS         !< Control structure of MOM_wave_drag

  ! Local variables

  ! The number and names of drag coefficients should match those of the streaming filters.
end subroutine wave_drag_init
module subroutine wave_drag_calc(u, v, drag_u, drag_v, G, CS)
  type(ocean_grid_type),           intent(in) :: G     !< The ocean's grid structure
  type(wave_drag_CS),              intent(in) :: CS    !< Control structure of MOM_wave_drag
  real, dimension(:,:,:), pointer, intent(in) :: u     !< Zonal velocity from the output of
                                                       !! streaming band-pass filters [L T-1 ~> m s-1]
  real, dimension(:,:,:), pointer, intent(in) :: v     !< Meridional velocity from the output of
                                                       !! streaming band-pass filters [L T-1 ~> m s-1]
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed), intent(out) :: drag_u !< Sum of products of filtered velocities
                                                       !! and scaled frequency-dependent drag [L2 T-2 ~> m2 s-2]
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB), intent(out) :: drag_v !< Sum of products of filtered velocities
                                                       !! and scaled frequency-dependent drag [L2 T-2 ~> m2 s-2]

  ! Local variables

end subroutine wave_drag_calc
  end interface

end module MOM_wave_drag
