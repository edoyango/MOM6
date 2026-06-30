! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface to CVMix interior shear schemes
module MOM_CVMix_shear

!> \author Brandon Reichl

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density
use CVMix_shear, only : CVMix_init_shear, CVMix_coeffs_shear
use MOM_kappa_shear, only : kappa_shear_is_used
implicit none ; private

#include <MOM_memory.h>

public calculate_CVMix_shear, CVMix_shear_init, CVMix_shear_is_used, CVMix_shear_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure including parameters for CVMix interior shear schemes.
type, public :: CVMix_shear_cs ; private
  logical :: use_LMD94                      !< Flags to use the LMD94 scheme
  logical :: use_PP81                       !< Flags to use Pacanowski and Philander (JPO 1981)
  integer :: n_smooth_ri                    !< Number of times to smooth Ri using a 1-2-1 filter
  real    :: Ri_zero                        !< LMD94 critical Richardson number [nondim]
  real    :: Nu_zero                        !< LMD94 maximum interior diffusivity [Z2 T-1 ~> m2 s-1]
  real    :: KPP_exp                        !< Exponent of unitless factor of diffusivities
                                            !! for KPP internal shear mixing scheme [nondim]
  real, allocatable, dimension(:,:,:) :: N2 !< Squared Brunt-Vaisala frequency [T-2 ~> s-2]
  real, allocatable, dimension(:,:,:) :: S2 !< Squared shear frequency [T-2 ~> s-2]
  real, allocatable, dimension(:,:,:) :: ri_grad !< Gradient Richardson number [nondim]
  real, allocatable, dimension(:,:,:) :: ri_grad_orig !< Gradient Richardson number
                                                      !! after smoothing [nondim]
  character(10) :: Mix_Scheme               !< Mixing scheme name (string)

  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to the diagnostics control structure
  !>@{ Diagnostic handles
  integer :: id_N2 = -1, id_S2 = -1, id_ri_grad = -1, id_kv = -1, id_kd = -1
  integer :: id_ri_grad_orig = -1
  !>@}

end type CVMix_shear_cs

character(len=40)  :: mdl = "MOM_CVMix_shear"  !< This module's name.


  interface
module subroutine calculate_CVMix_shear(u_H, v_H, h, tv, kd, kv, G, GV, US, CS )
  type(ocean_grid_type),                      intent(in)  :: G   !< Grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: u_H !< Initial zonal velocity on T points [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: v_H !< Initial meridional velocity on T
                                                                 !! points [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: kd !< The vertical diffusivity at each interface
                                                                 !! (not layer!) [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: kv !< The vertical viscosity at each interface
                                                                 !! (not layer!) [H Z T-1 ~> m2 s-1 or Pa s]
  type(CVMix_shear_cs),                       pointer     :: CS  !< The control structure returned by a previous
                                                                 !! call to CVMix_shear_init.
  ! Local variables

  ! some constants
end subroutine calculate_CVMix_shear
logical module function CVMix_shear_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< The current time.
  type(ocean_grid_type),   intent(in)    :: G  !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure.
  type(CVMix_shear_cs),    pointer       :: CS !< This module's control structure.
  ! Local variables

! This include declares and sets the variable "version".

end function CVMix_shear_init
logical module function CVMix_shear_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< Run-time parameter files handle.
  ! Local variables
end function CVMix_shear_is_used
module subroutine CVMix_shear_end(CS)
  type(CVMix_shear_cs), intent(inout) :: CS !< Control structure for this module that
                                            !! will be deallocated in this subroutine
end subroutine CVMix_shear_end
  end interface

end module MOM_CVMix_shear
