! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides regularization of layers in isopycnal mode
module MOM_regularize_layers

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public regularize_layers, regularize_layers_init

!> This control structure holds parameters used by the MOM_regularize_layers module
type, public :: regularize_layers_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: regularize_surface_layers !< If true, vertically restructure the
                             !! near-surface layers when they have too much
                             !! lateral variations to allow for sensible lateral
                             !! barotropic transports.
  logical :: reg_sfc_detrain !< If true, allow the buffer layers to detrain into the
                             !! interior as a part of the restructuring when
                             !! regularize_surface_layers is true
  real    :: density_match_tol !< A relative tolerance for how well the densities must match
                             !! with the target densities during detrainment when regularizing
                             !! the near-surface layers [nondim]
  real    :: sufficient_adjustment !< The fraction of the target entrainment of mass to the mixed
                             !! and buffer layers that is enough for one timestep when regularizing
                             !! the near-surface layers [nondim].  No more mass will be sought from
                             !! deeper layers in the interior after this fraction is exceeded.
  real    :: h_def_tol1      !< The value of the relative thickness deficit at
                             !! which to start modifying the structure, 0.5 by
                             !! default (or a thickness ratio of 5.83) [nondim].
  real    :: h_def_tol2      !< The value of the relative thickness deficit at
                             !! which to the structure modification is in full
                             !! force, now 20% of the way from h_def_tol1 to 1 [nondim].
  real    :: h_def_tol3      !< The value of the relative thickness deficit at which to start
                             !! detrainment from the buffer layers to the interior, now 30% of
                             !! the way from h_def_tol1 to 1 [nondim].
  real    :: h_def_tol4      !< The value of the relative thickness deficit at which to do
                             !! detrainment from the buffer layers to the interior at full
                             !! force, now 50% of the way from h_def_tol1 to 1 [nondim].
  real    :: Hmix_min        !< The minimum mixed layer thickness [H ~> m or kg m-2].
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in this module's
                             !! calculations.  Values below 20190101 recover the answers from the
                             !! end of 2018, while higher values use updated and more robust forms
                             !! of the same expressions.
  logical :: debug           !< If true, do more thorough checks for debugging purposes.

  integer :: id_def_rat = -1 !< A diagnostic ID
end type regularize_layers_CS

!>@{ Clock IDs
!! \todo Should these be global?
integer :: id_clock_pass
!>@}


  interface
module subroutine regularize_layers(h, tv, dt, ea, eb, G, GV, US, CS)
  type(ocean_grid_type),      intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),      intent(inout) :: tv !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL pointers.
  real,                       intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea !< The amount of fluid moved downward into a
                                                  !! layer; this should be increased due to mixed
                                                  !! layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb !< The amount of fluid moved upward into a layer
                                                  !! this should be increased due to mixed layer
                                                  !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(regularize_layers_CS), intent(in)    :: CS !< Regularize layer control structure

end subroutine regularize_layers
module subroutine regularize_surface(h, tv, dt, ea, eb, G, GV, US, CS)
  type(ocean_grid_type),      intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),      intent(inout) :: tv !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL pointers.
  real,                       intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea !< The amount of fluid moved downward into a
                                                  !! layer; this should be increased due to mixed
                                                  !! layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb !< The amount of fluid moved upward into a layer
                                                  !! this should be increased due to mixed layer
                                                  !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(regularize_layers_CS), intent(in)    :: CS !< Regularize layer control structure

  ! Local variables

                ! d_eb correspond to a gain in mass by a layer by upward motion.
end subroutine regularize_surface
module subroutine find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, h)
  type(ocean_grid_type),      intent(in)  :: G         !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)  :: GV        !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)  :: e         !< Interface depths [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),          &
                              intent(out) :: def_rat_u !< The thickness deficit ratio at u points,
                                                       !! [nondim].
  real, dimension(SZI_(G),SZJB_(G)),          &
                              intent(out) :: def_rat_v !< The thickness deficit ratio at v points,
                                                       !! [nondim].
  type(regularize_layers_CS), intent(in)  :: CS        !< Regularize layer control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                              intent(in)  :: h         !< Layer thicknesses [H ~> m or kg m-2].

  ! Local variables
                ! h_def_u is normalized [H ~> m or kg m-2].
                ! h_def_v is normalized [H ~> m or kg m-2].
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine find_deficit_ratios
module subroutine regularize_layers_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for
                                                 !! run-time parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate
                                                 !! diagnostic output.
  type(regularize_layers_CS), intent(inout) :: CS !< Regularize layer control structure

end subroutine regularize_layers_init
  end interface

end module MOM_regularize_layers
