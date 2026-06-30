! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization for the "Phillips" channel configuration
module Phillips_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public Phillips_initialize_thickness
public Phillips_initialize_velocity
public Phillips_initialize_sponges
public Phillips_initialize_topography

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

! This include declares and sets the variable "version".
#include "version_variable.h"


  interface
module subroutine Phillips_initialize_thickness(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h          !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file
                                                     !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read  !< If true, this call will only read
                                                     !! parameters without changing h.

                          ! geolat, often [km]
                          ! but this could be 1000 [m km-1]

end subroutine Phillips_initialize_thickness
module subroutine Phillips_initialize_velocity(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file to
                                                     !! parse for modelparameter values.
  logical,                 intent(in)  :: just_read  !< If true, this call will only read
                                                     !! parameters without changing u & v.

                          ! domain width [nondim]
                          ! as geolat, often [km]
                          ! domain width [nondim]
                          ! but this could be 1000 [m km-1]
                          ! Values below 20250101 recover the answers from the end of 2018, while
                          ! higher values use mathematically equivalent expressions that are fully
                          ! rescalable.
end subroutine Phillips_initialize_velocity
module subroutine Phillips_initialize_sponges(G, GV, US, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(unit_scale_type), intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure containing pointers
                                            !! to any available thermodynamic
                                            !! fields, potential temperature and
                                            !! salinity or mixed layer density.
                                            !! Absent fields have NULL ptrs.
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                            !! open file to parse for model
                                            !! parameter values.
  type(sponge_CS),   pointer    :: CSp      !< A pointer that is set to point to
                                            !! the control structure for the
                                            !! sponge module.
  real, intent(in), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< Thickness field [H ~> m or kg m-2].

  ! Local variables
                       ! geolat, often [km]
                          ! but this could be 1000 [m km-1]


end subroutine Phillips_initialize_sponges
module function sech(x)
  real, intent(in) :: x    !< Input value [nondim].
  real             :: sech !< Result [nondim].

  ! This is here to prevent overflows or underflows.
end function sech
module subroutine Phillips_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine Phillips_initialize_topography
  end interface

end module Phillips_initialization
