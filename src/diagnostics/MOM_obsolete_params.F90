! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Methods for testing for, and list of, obsolete run-time parameters.
module MOM_obsolete_params

! This module was first conceived and written by Robert Hallberg, July 2010.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, log_version, param_file_type

implicit none ; private

#include <MOM_memory.h>

public find_obsolete_params
public obsolete_logical, obsolete_int, obsolete_real, obsolete_char


  interface
module subroutine find_obsolete_params(param_file)
  type(param_file_type), intent(in) :: param_file !< Structure containing parameter file data.
  ! Local variables
! This include declares and sets the variable "version".

end subroutine find_obsolete_params
module subroutine obsolete_logical(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete LOGICAL parameter.
  logical,     optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  ! Local variables

end subroutine obsolete_logical
module subroutine obsolete_char(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname    !< Name of obsolete STRING parameter.
  character(len=*), optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint  !< A hint to the user about what to do.
  ! Local variables

end subroutine obsolete_char
module subroutine obsolete_real(param_file, varname, warning_val, hint, only_warn)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete REAL parameter.
  real,        optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  logical,     optional, intent(in) :: only_warn   !< If present and true, issue warnings instead of fatal errors.

  ! Local variables

end subroutine obsolete_real
module subroutine obsolete_int(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete INTEGER parameter.
  integer,     optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  ! Local variables

end subroutine obsolete_int
  end interface

end module MOM_obsolete_params
