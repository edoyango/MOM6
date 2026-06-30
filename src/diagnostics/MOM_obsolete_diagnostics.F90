! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides a mechanism for recording diagnostic variables that are no longer
!! valid, along with their replacement name if appropriate.
module MOM_obsolete_diagnostics

use MOM_diag_mediator, only : diag_ctrl, found_in_diagtable
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : param_file_type, log_version, get_param

implicit none ; private

#include <MOM_memory.h>

public register_obsolete_diagnostics


  interface
module subroutine register_obsolete_diagnostics(param_file, diag)
  type(param_file_type),       intent(in)    :: param_file !< The parameter file handle.
  type(diag_ctrl),             intent(in)    :: diag       !< A structure used to control diagnostics.
! This include declares and sets the variable "version".
  ! Local variables

end subroutine register_obsolete_diagnostics
logical module function diag_found(diag, varName, newVarName)
  type(diag_ctrl),            intent(in) :: diag       !< A structure used to control diagnostics.
  character(len=*),           intent(in) :: varName    !< The obsolete diagnostic name
  character(len=*), optional, intent(in) :: newVarName !< The valid name of this diagnostic

end function diag_found
  end interface

end module MOM_obsolete_diagnostics
