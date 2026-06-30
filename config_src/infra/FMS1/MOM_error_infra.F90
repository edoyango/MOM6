! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines for error handling and I/O management
module MOM_error_infra

use mpp_mod, only : mpp_error, mpp_pe, mpp_root_pe, mpp_stdlog=>stdlog, mpp_stdout=>stdout
use mpp_mod, only : NOTE, WARNING, FATAL

implicit none ; private

public :: MOM_err, is_root_pe, stdlog, stdout
!> Integer parameters encoding the severity of an error message
public :: NOTE, WARNING, FATAL


  interface
module subroutine MOM_err(severity, message)
  integer,           intent(in) :: severity !< The severity level of this error
  character(len=*),  intent(in) :: message  !< A message to write out

end subroutine MOM_err
integer module function stdout()
end function stdout
integer module function stdlog()
end function stdlog
logical module function is_root_pe()
end function is_root_pe
  end interface

end module MOM_error_infra
