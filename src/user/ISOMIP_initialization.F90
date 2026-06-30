! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the ISOMIP test case.
module ISOMIP_initialization

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_SIGMA_SHELF_ZSTAR
implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "ISOMIP_initialization" !< This module's name.

! The following routines are visible to the outside world
public ISOMIP_initialize_topography
public ISOMIP_initialize_thickness
public ISOMIP_initialize_temperature_salinity
public ISOMIP_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine ISOMIP_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  ! The following variables are used to set up the bathymetry in the ISOMIP example.
  ! This include declares and sets the variable "version".
end subroutine ISOMIP_initialize_topography
module subroutine ISOMIP_initialize_thickness ( h, depth_tot, G, GV, US, param_file, tv, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure to parse for model parameter values
  type(thermo_var_ptrs),   intent(in)  :: tv          !< A structure containing pointers to any
                                                      !! available thermodynamic fields, including
                                                      !! the eqn. of state.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                          !  usually negative because it is positive upward.
                            ! positive upward, in depth units [Z ~> m].
  !character(len=256) :: mesg  ! The text of an error message

end subroutine ISOMIP_initialize_thickness
module subroutine ISOMIP_initialize_temperature_salinity ( T, S, h, depth_tot, G, GV, US, param_file, &
                                                    eqn_of_state, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: depth_tot  !< The nominal total bottom-to-top
                                                               !! depth of the ocean [Z ~> m]
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            intent(in)  :: eqn_of_state !< Equation of state structure
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables
  !real :: rho_tmp    ! A temporary density used for debugging [R ~> kg m-3]
  !character(len=256) :: mesg ! The text of an error message

end subroutine ISOMIP_initialize_temperature_salinity
module subroutine ISOMIP_initialize_sponges(G, GV, US, tv, depth_tot, PF, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure containing pointers to any available
                                              !! thermodynamic fields, potential temperature and
                                              !! salinity or mixed layer density.
                                              !! Absent fields have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: PF   !< A structure to parse for model parameter values
  logical,                 intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),         pointer    :: CSp  !< Layer-mode sponge structure
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE-mode sponge structure
  ! Local variables
  ! real :: RHO(SZI_(G),SZJ_(G),SZK_(GV)) ! A temporary array for RHO [R ~> kg m-3]

                                    ! negative because it is positive upward.
  !real :: rho_tmp                   ! A temporary density used for debugging [R ~> kg m-3]


end subroutine ISOMIP_initialize_sponges
  end interface

end module ISOMIP_initialization
