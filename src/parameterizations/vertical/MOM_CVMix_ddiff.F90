! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface to CVMix double diffusion scheme.
module MOM_CVMix_ddiff

use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density_derivs
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_debugging,      only : hchksum
use MOM_grid,           only : ocean_grid_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use cvmix_ddiff,        only : cvmix_init_ddiff, CVMix_coeffs_ddiff
use cvmix_kpp,          only : CVmix_kpp_compute_kOBL_depth
implicit none ; private

#include <MOM_memory.h>

public CVMix_ddiff_init, CVMix_ddiff_end, CVMix_ddiff_is_used, compute_ddiff_coeffs

!> Control structure including parameters for CVMix double diffusion.
type, public :: CVMix_ddiff_cs ; private

  ! Parameters
  real    :: strat_param_max !< maximum value for the stratification parameter [nondim]
  real    :: kappa_ddiff_s   !< leading coefficient in formula for salt-fingering regime
                             !! for salinity diffusion [Z2 T-1 ~> m2 s-1]
  real    :: ddiff_exp1      !< interior exponent in salt-fingering regime formula [nondim]
  real    :: ddiff_exp2      !< exterior exponent in salt-fingering regime formula [nondim]
  real    :: mol_diff        !< molecular diffusivity [Z2 T-1 ~> m2 s-1]
  real    :: kappa_ddiff_param1 !< exterior coefficient in diffusive convection regime [nondim]
  real    :: kappa_ddiff_param2 !< middle coefficient in diffusive convection regime [nondim]
  real    :: kappa_ddiff_param3 !< interior coefficient in diffusive convection regime [nondim]
  real    :: min_thickness      !< Minimum thickness allowed [H ~> m or kg m-2]
  character(len=4) :: diff_conv_type !< type of diffusive convection to use. Options are Marmorino &
                                !! Caldwell 1976 ("MC76"; default) and Kelley 1988, 1990 ("K90")
  logical :: debug              !< If true, turn on debugging

end type CVMix_ddiff_cs

character(len=40)  :: mdl = "MOM_CVMix_ddiff"     !< This module's name.


  interface
logical module function CVMix_ddiff_init(Time, G, GV, US, param_file, diag, CS)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(CVMix_ddiff_cs),    pointer       :: CS         !< This module's control structure.

  ! This include declares and sets the variable "version".

end function CVMix_ddiff_init
module subroutine compute_ddiff_coeffs(h, tv, G, GV, US, j, Kd_T, Kd_S, CS, R_rho)

  type(ocean_grid_type),                      intent(in)    :: G    !< Grid structure.
  type(verticalGrid_type),                    intent(in)    :: GV   !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                      intent(in)    :: tv   !< Thermodynamics structure.
  type(unit_scale_type),                      intent(in)    :: US   !< A dimensional unit scaling type
  integer,                                    intent(in)    :: j    !< Meridional grid index to work on.
  ! Kd_T and Kd_S are intent inout because only one j-row is set here, but they are essentially outputs.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kd_T !< Interface double diffusion diapycnal
                                                                    !! diffusivity for temperature
                                                                    !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kd_S !< Interface double diffusion diapycnal
                                                                    !! diffusivity for salinity
                                                                    !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(CVMix_ddiff_cs),                       pointer       :: CS   !< The control structure returned
                                                                    !! by a previous call to CVMix_ddiff_init.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                    optional, intent(inout) :: R_rho !< The density ratios at interfaces [nondim].

  ! Local variables


  ! initialize dummy variables
end subroutine compute_ddiff_coeffs
logical module function CVMix_ddiff_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
end function CVMix_ddiff_is_used
module subroutine CVMix_ddiff_end(CS)
  type(CVMix_ddiff_cs), pointer :: CS !< Control structure for this module that
                                      !! will be deallocated in this subroutine
end subroutine CVMix_ddiff_end
  end interface

end module MOM_CVMix_ddiff
