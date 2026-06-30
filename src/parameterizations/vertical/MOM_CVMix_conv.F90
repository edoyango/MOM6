! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface to CVMix convection scheme.
module MOM_CVMix_conv

use MOM_debugging,      only : hchksum
use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use CVMix_convection,   only : CVMix_init_conv, CVMix_coeffs_conv
use CVMix_kpp,          only : CVMix_kpp_compute_kOBL_depth

implicit none ; private

#include <MOM_memory.h>

public CVMix_conv_init, calculate_CVMix_conv, CVMix_conv_is_used

!> Control structure including parameters for CVMix convection.
type, public :: CVMix_conv_cs ; private

  ! Parameters
  real    :: kd_conv_const !< diffusivity constant used in convective regime [Z2 T-1 ~> m2 s-1]
  real    :: kv_conv_const !< viscosity constant used in convective regime [Z2 T-1 ~> m2 s-1]
  real    :: bv_sqr_conv   !< Threshold for squared buoyancy frequency
                           !! needed to trigger Brunt-Vaisala parameterization [T-2 ~> s-2]
  real    :: min_thickness !< Minimum thickness allowed [Z ~> m]
  logical :: debug         !< If true, turn on debugging

  ! Diagnostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to diagnostics control structure
  !>@{ Diagnostics handles
  integer :: id_N2 = -1, id_kd_conv = -1, id_kv_conv = -1
  !>@}

end type CVMix_conv_cs

character(len=40)  :: mdl = "MOM_CVMix_conv"     !< This module's name.


  interface
logical module function CVMix_conv_init(Time, G, GV, US, param_file, diag, CS)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(CVMix_conv_cs),     intent(inout) :: CS         !< CVMix convection control structure


  ! This include declares and sets the variable "version".

  ! Read parameters
end function CVMix_conv_init
module subroutine calculate_CVMix_conv(h, tv, G, GV, US, CS, hbl, Kd, Kv, Kd_aux)

  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                     intent(in)  :: tv !< Thermodynamics structure.
  type(CVMix_conv_cs),                       intent(in)  :: CS !< CVMix convection control structure
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: hbl !< Depth of ocean boundary layer [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(inout) :: Kd !< Diapycnal diffusivity at each interface
                                                                 !! that will be incremented here
                                                                 !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(inout) :: Kv !< Viscosity at each interface that will be
                                                                 !! incremented here [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                   optional, intent(inout) :: Kd_aux !< A second diapycnal diffusivity at each
                                                                 !! interface that will also be incremented
                                                                 !! here [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! local variables
                                       !! variable since here convection is always
                                       !! computed based on Brunt Vaisala.
                                       !! a dummy variable, same reason as above.
                    ! [H s-2 R-1 ~> m4 s-2 kg-1 or m s-2]

end subroutine calculate_CVMix_conv
logical module function CVMix_conv_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
end function CVMix_conv_is_used
  end interface

end module MOM_CVMix_conv
