! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Increments the diapycnal diffusivity in a specified band of latitudes and densities.
module user_change_diffusivity

use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public user_change_diff, user_change_diff_init, user_change_diff_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for user_change_diffusivity
type, public :: user_change_diff_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real :: Kd_add        !< The scale of a diffusivity that is added everywhere without
                        !! any filtering or scaling [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: lat_range(4)  !< 4 values that define the latitude range over which
                        !! a diffusivity scaled by Kd_add is added [degrees_N].
  real :: rho_range(4)  !< 4 values that define the coordinate potential
                        !! density range over which a diffusivity scaled by
                        !! Kd_add is added [R ~> kg m-3].
  logical :: use_abs_lat  !< If true, use the absolute value of latitude when
                          !! setting lat_range.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                          !! regulate the timing of diagnostic output.
end type user_change_diff_CS


  interface
module subroutine user_change_diff(h, tv, G, GV, US, CS, Kd_lay, Kd_int, T_f, S_f, Kd_int_add)
  type(ocean_grid_type),                    intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                  intent(in)    :: GV  !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: h   !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                    intent(in)    :: tv  !< A structure containing pointers
                                                                 !! to any available thermodynamic
                                                                 !! fields. Absent fields have NULL ptrs.
  type(unit_scale_type),                    intent(in)    :: US  !< A dimensional unit scaling type
  type(user_change_diff_CS),                pointer       :: CS  !< This module's control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity of each
                                                                  !! layer [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at each
                                                                  !! interface [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   optional, intent(in)    :: T_f !< Temperature with massless
                                                                  !! layers filled in vertically [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   optional, intent(in)    :: S_f !< Salinity with massless
                                                                  !! layers filled in vertically [S ~> ppt].
  real, dimension(:,:,:),                      optional, pointer       :: Kd_int_add !< The diapycnal
                                                                  !! diffusivity that is being added at
                                                                  !! each interface [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  ! Local variables
                      ! equation of state.


end subroutine user_change_diff
module function range_OK(range) result(OK)
  real, dimension(4), intent(in) :: range  !< Four values to check [arbitrary]
  logical                        :: OK     !< Return value.

end function range_OK
module function val_weights(val, range) result(ans)
  real,               intent(in) :: val    !< Value for which we need an answer [arbitrary units].
  real, dimension(4), intent(in) :: range  !< Range over which the answer is non-zero [arbitrary units].
  real                           :: ans    !< Return value [nondim].
  ! Local variables

end function val_weights
module subroutine user_change_diff_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),           intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),     intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),     intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure indicating the
                                                         !! open file to parse for
                                                         !! model parameter values.
  type(diag_ctrl), target,   intent(inout) :: diag       !< A structure that is used to
                                                         !! regulate diagnostic output.
  type(user_change_diff_CS), pointer       :: CS         !< A pointer that is set to
                                                         !! point to the control
                                                         !! structure for this module.

  ! This include declares and sets the variable "version".

end subroutine user_change_diff_init
module subroutine user_change_diff_end(CS)
  type(user_change_diff_CS), pointer :: CS !< A pointer that is set to point to the control
                                           !! structure for this module.

end subroutine user_change_diff_end
  end interface

end module user_change_diffusivity
