! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initializes hydrography from z-coordinate climatology files
module MOM_tracer_initialization_from_Z

use MOM_debugging,     only : hchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, param_file_type, log_version
use MOM_grid,          only : ocean_grid_type
use MOM_horizontal_regridding, only : myStats, horiz_interp_and_extrap_tracer
use MOM_interface_heights, only : dz_to_thickness_simple
use MOM_regridding,    only : set_dz_neglect, set_h_neglect
use MOM_remapping,     only : remapping_CS, initialize_remapping
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_ALE,           only : ALE_remap_scalar

implicit none ; private

#include <MOM_memory.h>

public :: MOM_initialize_tracer_from_Z

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40)  :: mdl = "MOM_tracer_initialization_from_Z" !< This module's name.


  interface
module subroutine MOM_initialize_tracer_from_Z(h, tr, G, GV, US, PF, src_file, src_var_nam, &
                          src_var_unit_conversion, src_var_record, homogenize, &
                          useALEremapping, remappingScheme, src_var_gridspec, h_in_Z_units, &
                          ongrid)
  type(ocean_grid_type),      intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type),    intent(in)    :: GV  !< Ocean vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h   !< Layer thicknesses, in [H ~> m or kg m-2] or
                                                   !! [Z ~> m] depending on the value of h_in_Z_units.
  real, dimension(:,:,:),     pointer       :: tr  !< Pointer to array to be initialized [CU ~> conc]
  type(param_file_type),      intent(in)    :: PF  !< parameter file
  character(len=*),           intent(in)    :: src_file !< source filename
  character(len=*),           intent(in)    :: src_var_nam !< variable name in file
  real,             optional, intent(in)    :: src_var_unit_conversion !< optional multiplicative unit conversion,
                                                   !! often used for rescaling into model units [CU conc-1 ~> 1]
  integer,          optional, intent(in)    :: src_var_record  !< record to read for multiple time-level files
  logical,          optional, intent(in)    :: homogenize !< optionally homogenize to mean value
  logical,          optional, intent(in)    :: useALEremapping !< to remap or not (optional)
  character(len=*), optional, intent(in)    :: remappingScheme !< remapping scheme to use.
  character(len=*), optional, intent(in)    :: src_var_gridspec !< Source variable name in a gridspec file.
                                                                !! This is not implemented yet.
  logical,          optional, intent(in)    :: h_in_Z_units !< If present and true, the input grid
                                                            !! thicknesses are in the units of height
                                                            !! ([Z ~> m]) instead of the usual units of
                                                            !! thicknesses ([H ~> m or kg m-2])
  logical,          optional, intent(in)    :: ongrid       !< If true, then data are assumed to have been
                                                            !! interpolated to the model horizontal grid. In this case,
                                                            !! only extrapolation is performed by
                                                            !! horiz_interp_and_extrap_tracer()
  ! Local variables

  ! This include declares and sets the variable "version".

                                                        ! and input-file vertical levels [CU ~> conc]
                                                        ! and input-file vertical levels [nondim]

  ! Local variables for ALE remapping

                                  ! remapping cell reconstructions [Z ~> m] or [H ~> m or kg m-2]
                                  ! remapping edge value calculations [Z ~> m] or [H ~> m or kg m-2]
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
                                  ! for horizontal regridding.  Values below 20190101 recover the
                                  ! answers from 2018, while higher values use expressions that have
                                  ! been rearranged for rotational invariance.

end subroutine MOM_initialize_tracer_from_Z
  end interface

end module MOM_tracer_initialization_from_Z
