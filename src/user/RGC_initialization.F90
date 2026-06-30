! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the models sponges for the Rotating Gravity Current (RGC) experiment.
module RGC_initialization

!***********************************************************************
!* By Elizabeth Yankovsky, May 2018                                    *
!***********************************************************************

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_ALE_sponge, only : set_up_ALE_sponge_vel_field
use MOM_domains, only : pass_var
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_sponge, only : set_up_sponge_ML_density
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, EOS_domain
implicit none ; private

#include <MOM_memory.h>

public RGC_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine RGC_initialize_sponges(G, GV, US, tv, u, v, depth_tot, PF, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< A structure containing pointers
                                            !! to any available thermodynamic
                                            !! fields, potential temperature and
                                            !! salinity or mixed layer density.
                                            !! Absent fields have NULL pointers.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 target, intent(in) :: u    !< Array with the u velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 target, intent(in) :: v    !< Array with the v velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type), intent(in) :: PF   !< A structure to parse for model parameter values.
  logical,               intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),       pointer    :: CSp  !< Layer-mode sponge structure
  type(ALE_sponge_CS),   pointer    :: ACSp !< ALE-mode sponge structure

  ! Local variables


end subroutine RGC_initialize_sponges
  end interface

end module RGC_initialization
