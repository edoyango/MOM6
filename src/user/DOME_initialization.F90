! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the model for the "DOME" experiment.
!! DOME = Dynamics of Overflows and Mixing Experiment
module DOME_initialization

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_tracer_registry, only : tracer_name_lookup
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public DOME_initialize_topography
public DOME_initialize_thickness
public DOME_initialize_sponges
public DOME_set_OBC_data, register_DOME_OBC

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine DOME_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
                          ! but this could be 1000 [m km-1]
  ! This include declares and sets the variable "version".
end subroutine DOME_initialize_topography
module subroutine DOME_initialize_thickness(h, depth_tot, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

                            ! negative because it is positive upward [Z ~> m].
                            ! positive upward [Z ~> m].

end subroutine DOME_initialize_thickness
module subroutine DOME_initialize_sponges(G, GV, US, tv, depth_tot, PF, CSp)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure containing any available
                                              !! thermodynamic fields, including potential
                                              !! temperature and salinity or mixed layer density.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: PF   !< A structure indicating the open file to
                                              !! parse for model parameter values.
  type(sponge_CS),         pointer    :: CSp  !< A pointer that is set to point to the control
                                              !! structure for this module.



end subroutine DOME_initialize_sponges
module subroutine register_DOME_OBC(param_file, US, OBC, tr_Reg)
  type(param_file_type),      intent(in) :: param_file !< parameter file.
  type(unit_scale_type),      intent(in) :: US       !< A dimensional unit scaling type
  type(ocean_OBC_type),       pointer    :: OBC      !< OBC registry.
  type(tracer_registry_type), pointer    :: tr_Reg   !< Tracer registry.

end subroutine register_DOME_OBC
module subroutine DOME_set_OBC_data(OBC, tv, G, GV, US, PF, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv  !< A structure containing pointers to any
                              !! available thermodynamic fields, including potential
                              !! temperature and salinity or mixed layer density. Absent
                              !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: PF  !< A structure indicating the open file
                              !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

  ! Local variables
  ! The following variables are used to set up the transport in the DOME example.
                            ! top and bottom of a layer [nondim]
                            ! with a range from 0 for the densest water to -1 for the lightest
                            ! v-velocity face, in the same units as G%geoLon [km]
                            ! inner edge of the inflow
                            ! properties [T-1 ~> s-1]
                            ! region of the specified shear profile [nondim]
                            ! but this could be 1000 [m km-1]

end subroutine DOME_set_OBC_data
  end interface

end module DOME_initialization
