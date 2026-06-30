! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization for the "bench mark" configuration
module benchmark_initialization

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public benchmark_initialize_topography
public benchmark_initialize_thickness
public benchmark_init_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine benchmark_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  ! This include declares and sets the variable "version".
end subroutine benchmark_initialize_topography
module subroutine benchmark_initialize_thickness(h, depth_tot, G, GV, US, param_file, eqn_of_state, &
                                          P_Ref, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  type(EOS_type),          intent(in)  :: eqn_of_state !< Equation of state structure
  real,                    intent(in)  :: P_Ref       !< The coordinate-density
                                                      !! reference pressure [R L2 T-2 ~> Pa].
  logical,                 intent(in)  :: just_read   !< If true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
                             ! usually negative because it is positive upward.
                             ! in depth units [Z ~> m].
                             ! positive upward, in depth units [Z ~> m].
                    ! between SST and the bottom temperature [nondim].
                    ! interface temperature for a given z [nondim]
                    ! temperature and the interface temperature with z [Z-1 ~> m-1]
  ! This include declares and sets the variable "version".

end subroutine benchmark_initialize_thickness
module subroutine benchmark_init_temperature_salinity(T, S, G, GV, US, param_file, &
               eqn_of_state, P_Ref, just_read)
  type(ocean_grid_type),               intent(in)  :: G            !< The ocean's grid structure
  type(verticalGrid_type),             intent(in)  :: GV           !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T      !< The potential temperature
                                                                   !! that is being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S      !< The salinity that is being
                                                                   !! initialized [S ~> ppt]
  type(unit_scale_type),               intent(in)  :: US           !< A dimensional unit scaling type
  type(param_file_type),               intent(in)  :: param_file   !< A structure indicating the
                                                                   !! open file to parse for
                                                                   !! model parameter values.
  type(EOS_type),                      intent(in)  :: eqn_of_state !< Equation of state structure
  real,                                intent(in)  :: P_Ref        !< The coordinate-density
                                                                   !! reference pressure [R L2 T-2 ~> Pa]
  logical,                             intent(in)  :: just_read    !< If true, this call will only read
                                                                   !! parameters without changing T & S.
  ! Local variables

end subroutine benchmark_init_temperature_salinity
  end interface

end module benchmark_initialization
