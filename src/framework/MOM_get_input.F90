! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> \brief Reads the only Fortran name list needed to boot-strap the model.
!!
!! The name list parameters indicate which directories to use for
!! certain types of input and output, and which files to look in for
!! the full parsable input parameter file(s).
module MOM_get_input

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : open_param_file, param_file_type
use MOM_io, only : file_exists, close_file, slasher, ensembler
use MOM_io, only : open_namelist_file, check_nml_error
use posix, only : mkdir, stat, stat_buf

implicit none ; private

public get_MOM_input

!> Container for paths and parameter file names.
type, public :: directories
  character(len=240) :: &
    restart_input_dir = ' ',& !< The directory to read restart and input files.
    restart_output_dir = ' ',&!< The directory into which to write restart files.
    output_directory = ' '    !< The directory to use to write the model output.
  character(len=2048) :: &
    input_filename  = ' '     !< A string that indicates the input files or how
                              !! the run segment should be started.
end type directories


  interface
module subroutine get_MOM_input(param_file, dirs, check_params, default_input_filename, ensemble_num)
  type(param_file_type), optional, intent(out) :: param_file   !< A structure to parse for run-time parameters.
  type(directories),     optional, intent(out) :: dirs         !< Container for paths and parameter file names.
  logical,               optional, intent(in)  :: check_params !< If present and False will stop error checking for
                                                               !! run-time parameters.
  character(len=*),      optional, intent(in)  :: default_input_filename !< If present, is the value assumed for
                                                               !! input_filename if input_filename is not listed
                                                               !! in the namelist MOM_input_nml.
  integer, optional, intent(in) :: ensemble_num !< The ensemble id of the current member
  ! Local variables

                               ! the run segment should be started.



  ! Default values in case parameter is not set in file input.nml
end subroutine get_MOM_input
  end interface

end module MOM_get_input
