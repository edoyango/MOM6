! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Module for calculating curve fit for porous topography.
!written by sjd
module MOM_porous_barriers

use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_error_handler,     only : MOM_error, FATAL
use MOM_grid,              only : ocean_grid_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs, porous_barrier_type
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_interface_heights, only : find_eta
use MOM_time_manager,      only : time_type
use MOM_diag_mediator,     only : register_diag_field, diag_ctrl, post_data
use MOM_file_parser,       only : param_file_type, get_param, log_version
use MOM_unit_scaling,      only : unit_scale_type
use MOM_debugging,         only : hchksum, uvchksum

implicit none ; private

public porous_widths_layer, porous_widths_interface, porous_barriers_init

#include <MOM_memory.h>

!> The control structure for the MOM_porous_barriers module
type, public :: porous_barrier_CS ; private
  logical :: initialized = .false.  !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: &
      diag => Null()                !< A structure to regulate diagnostic output timing
  logical :: debug                  !< If true, write verbose checksums for debugging purposes.
  real    :: mask_depth             !< The depth shallower than which porous barrier is not applied [Z ~> m]
  integer :: eta_interp             !< An integer indicating how the interface heights at the velocity
                                    !! points are calculated. Valid values are given by the parameters
                                    !! defined below: MAX, MIN, ARITHMETIC and HARMONIC.
  integer :: answer_date            !< The vintage of the porous barrier weight function calculations.
                                    !! Values below 20220806 recover the old answers in which the layer
                                    !! averaged weights are not strictly limited by an upper-bound of 1.0 .
  !>@{ Diagnostic IDs
  integer :: id_por_layer_widthU = -1, id_por_layer_widthV = -1, &
             id_por_face_areaU = -1, id_por_face_areaV = -1
  !>@}
end type porous_barrier_CS

integer :: id_clock_porous_barrier !< CPU clock for porous barrier

!>@{ Enumeration values for eta interpolation schemes
integer, parameter :: ETA_INTERP_MAX   = 1
integer, parameter :: ETA_INTERP_MIN   = 2
integer, parameter :: ETA_INTERP_ARITH = 3
integer, parameter :: ETA_INTERP_HARM  = 4
character(len=20), parameter :: ETA_INTERP_MAX_STRING = "MAX"
character(len=20), parameter :: ETA_INTERP_MIN_STRING = "MIN"
character(len=20), parameter :: ETA_INTERP_ARITH_STRING = "ARITHMETIC"
character(len=20), parameter :: ETA_INTERP_HARM_STRING = "HARMONIC"
!>@}


  interface
module subroutine porous_widths_layer(h, tv, G, GV, US, pbv, CS, eta_bt)
  ! Note: eta_bt is not currently used
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in) :: tv  !< A structure pointing to various
                                                                !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  type(porous_barrier_type),                  intent(inout) :: pbv !< porous barrier fractional cell metrics
  type(porous_barrier_CS),                    intent(in) :: CS     !< Control structure for porous barrier

  !local variables
                                                     ! to the previous layer at u or v points [Z ~> m]
                                                ! updated while moving up layers

end subroutine porous_widths_layer
module subroutine porous_widths_interface(h, tv, G, GV, US, pbv, CS, eta_bt)
  ! Note: eta_bt is not currently used
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in) :: tv  !< A structure pointing to various
                                                                !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  type(porous_barrier_type),                  intent(inout) :: pbv  !< porous barrier fractional cell metrics
  type(porous_barrier_CS),                    intent(in) :: CS !< Control structure for porous barrier

  !local variables
                                                ! updated while moving up layers

end subroutine porous_widths_interface
module subroutine calc_eta_at_uv(eta_u, eta_v, interp, dmask, h, tv, G, GV, US, eta_bt)
  !variables needed to call find_eta
  type(ocean_grid_type),                        intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                      intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                        intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(in) :: tv  !< A structure pointing to various
                                                                  !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), optional,   intent(in) :: eta_bt !< optional barotropic variable
                                                                   !! used to dilate the layer thicknesses
                                                                   !! [H ~> m or kg m-2].
  real,                                         intent(in) :: dmask !< The depth shallower than which
                                                                    !! porous barrier is not applied [Z ~> m]
  integer,                                      intent(in) :: interp !< eta interpolation method
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: eta_u !< Layer interface heights at u points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: eta_v !< Layer interface heights at v points [Z ~> m]

  ! local variables

end subroutine calc_eta_at_uv
module subroutine calc_por_layer(D_min, D_max, D_avg, eta_layer, A_layer, do_next)
  real,    intent(in)  :: D_min     !< minimum topographic height (deepest) [Z ~> m]
  real,    intent(in)  :: D_max     !< maximum topographic height (shallowest) [Z ~> m]
  real,    intent(in)  :: D_avg     !< mean topographic height [Z ~> m]
  real,    intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real,    intent(out) :: A_layer   !< frac. open face area of below eta_layer [Z ~> m]
  logical, intent(out) :: do_next   !< False if eta_layer>D_max

  ! local variables

end subroutine calc_por_layer
module subroutine calc_por_interface(D_min, D_max, D_avg, eta_layer, w_layer, do_next)
  real,    intent(in)  :: D_min     !< minimum topographic height (deepest) [Z ~> m]
  real,    intent(in)  :: D_max     !< maximum topographic height (shallowest) [Z ~> m]
  real,    intent(in)  :: D_avg     !< mean topographic height [Z ~> m]
  real,    intent(in)  :: eta_layer !< height of interface [Z ~> m]
  real,    intent(out) :: w_layer   !< frac. open interface width at eta_layer [nondim]
  logical, intent(out) :: do_next   !< False if eta_layer>D_max

  ! local variables

end subroutine calc_por_interface
module subroutine porous_barriers_init(Time, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time       !< Current model time
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure
  type(porous_barrier_CS), intent(inout) :: CS         !< Module control structure

  ! local variables
  !> This include declares and sets the variable "version".

end subroutine porous_barriers_init
  end interface

end module MOM_porous_barriers
