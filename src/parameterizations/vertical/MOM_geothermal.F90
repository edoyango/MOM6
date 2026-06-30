! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implemented geothermal heating at the ocean bottom.
module MOM_geothermal

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,           only : calculate_specific_vol_derivs

implicit none ; private

#include <MOM_memory.h>

public geothermal_entraining, geothermal_in_place, geothermal_init, geothermal_end

!> Control structure for geothermal heating
type, public :: geothermal_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: dRcv_dT_inplace  !< The value of dRcv_dT above which (dRcv_dT is negative) the
                              !! water is heated in place instead of moving upward between
                              !! layers in non-ALE layered mode [R C-1 ~> kg m-3 degC-1]
  real, allocatable, dimension(:,:) :: geo_heat !< The geothermal heat flux [Q R Z T-1 ~> W m-2]
  real    :: geothermal_thick !< The thickness over which geothermal heating is
                              !! applied [H ~> m or kg m-2]
  logical :: apply_geothermal !< If true, geothermal heating will be applied.  This is false if
                              !! GEOTHERMAL_SCALE is 0 and there is no heat to apply.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing
                                             !! timing of diagnostic output
  integer :: id_internal_heat_heat_tendency = -1  !< ID for diagnostic of heat tendency
  integer :: id_internal_heat_temp_tendency = -1  !< ID for diagnostic of temperature tendency
  integer :: id_internal_heat_h_tendency = -1     !< ID for diagnostic of thickness tendency
  integer :: id_geothermal_buoyancy_flux = -1     !< ID for diagnostic of bottom buoyancy flux
end type geothermal_CS


  interface
module subroutine geothermal_entraining(h, tv, dt, ea, eb, G, GV, US, CS, halo)
  type(ocean_grid_type),                     intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: ea !< The amount of fluid moved
                                                                 !! downward into a layer; this
                                                                 !! should be increased due to mixed
                                                                 !! layer detrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: eb !< The amount of fluid moved upward
                                                                 !! into a layer; this should be
                                                                 !! increased due to mixed layer
                                                                 !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(geothermal_CS),                       intent(in)    :: CS !< The control structure returned by
                                                                 !! a previous call to
                                                                 !! geothermal_init.
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work
  ! Local variables


                        ! in the present layer [R C-1 ~> kg m-3 degC-1]; usually negative
                        ! [C H ~> degC m or degC kg m-2]
                        ! layer [C H ~> degC m or degC kg m-2]
                        ! [C H ~> degC m or degC kg m-2]
                        ! 0 <= heating <= heat_trans
                        ! [C H Q-1 R-1 Z-1 ~> degC m3 J-1 or degC kg J-1]



end subroutine geothermal_entraining
module subroutine geothermal_in_place(h, tv, dt, G, GV, US, CS, BFlx_geothermal, halo)
  type(ocean_grid_type),                     intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(geothermal_CS),                       intent(in)    :: CS !< Geothermal heating control struct
  real, dimension(SZI_(G), SZJ_(G)),         intent(out)   :: BFlx_geothermal !< Geothermal buoyancy flux
                                                                 !! in [Z2 T-3 ~> m2 s-3]
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work


  ! Local variables

                        ! [C H Q-1 R-1 Z-1 ~> degC m3 J-1 or degC kg J-1]

                        ! converted into a layer-integrated heat tendency [Q R Z T-1 ~> W m-2]

end subroutine geothermal_in_place
module subroutine geothermal_init(Time, G, GV, US, param_file, diag, CS, useALEalgorithm)
  type(time_type), target, intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure used to regulate diagnostic output.
  type(geothermal_CS),     intent(inout) :: CS   !< Geothermal heating control struct
  logical,       optional, intent(in)    :: useALEalgorithm  !< logical for whether to use ALE remapping

! This include declares and sets the variable "version".
  ! Local variables
                     ! [Q R Z T-1 ~> W m-2] or [Q R Z m2 s J-1 T-1 ~> nondim]
end subroutine geothermal_init
module subroutine geothermal_end(CS)
  type(geothermal_CS), intent(inout) :: CS !< Geothermal heating control struct
end subroutine geothermal_end
  end interface

end module MOM_geothermal
