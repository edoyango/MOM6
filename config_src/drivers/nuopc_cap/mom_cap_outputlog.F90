!> This module contains a set of subroutines that check if MOM restart and history files
!! have been written and closed. This file is specific to UWM operational requirements
!! and configurations (eg specific output frequencies in hours) and may break if used outside
!! the scope of intended use.
!! This module is a stub when CESMCOUPLED is defined
module MOM_cap_outputlog

#ifdef CESMCOUPLED
use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_SUCCESS
implicit none ; private

public :: outputlog_init, outputlog_run, outputlog_restart

  interface
module subroutine outputlog_init(gcomp, mclock, rc)
  type(ESMF_GridComp)  :: gcomp  !< an ESMF_GridComp object
  type(ESMF_Clock)     :: mclock !< the ESMF_clock for the model
  integer, intent(out) :: rc     !< return code
end subroutine outputlog_init
module subroutine outputlog_run(mclock, atStopTime, rc)
  type(ESMF_Clock)              :: mclock     !< the ESMF_clock for the model
  logical, intent(in), optional :: atStopTime !< if true, checks for final output file
  integer, intent(out)          :: rc         !< return code
end subroutine outputlog_run
module subroutine outputlog_restart(mclock, num_rest_files, rc)
  type(ESMF_Clock)     :: mclock         !< the ESMF_clock for the model
  integer, intent(in)  :: num_rest_files !< the number of restart files
  integer, intent(out) :: rc             !< return code
end subroutine outputlog_restart
module subroutine outputlog_init(gcomp, mclock, rc)
  type(ESMF_GridComp)  :: gcomp
  type(ESMF_Clock)     :: mclock
  integer, intent(out) :: rc

  ! local variables
  !----------------------------------------------------------------------------

end subroutine outputlog_init
module subroutine outputlog_run(mclock, atStopTime, rc)
  type(ESMF_Clock)              :: mclock
  logical, intent(in), optional :: atStopTime
  integer, intent(out)          :: rc

  ! local variables
  !----------------------------------------------------------------------------

end subroutine outputlog_run
module subroutine outputlog_restart(mclock, num_rest_files, rc)
  type(ESMF_Clock)     :: mclock
  integer, intent(in)  :: num_rest_files
  integer, intent(out) :: rc

  ! local variables
  !----------------------------------------------------------------------------

end subroutine outputlog_restart
logical module function file_is_complete(fname, chk4size, createsize, rc) result(filecomplete)
  character(len=*), intent(in)  :: fname
  logical,          intent(in)  :: chk4size
  integer,          intent(in)  :: createsize
  integer,          intent(out) :: rc

  !----------------------------------------------------------------------------

end function file_is_complete
integer module function get_unlimited_len(fname) result(unlen)
  character(len=*), intent(in) :: fname

  !----------------------------------------------------------------------------

end function get_unlimited_len
character(len=16) module function get_timestr(MyTime, rc) result(timestr)
  type(ESMF_Time), intent(in)  :: MyTime
  integer,         intent(out) :: rc

  !----------------------------------------------------------------------------

end function get_timestr
character(len=40) module function get_importexport(currTime, nextTime, rc) result(importexport)

  type(ESMF_Time), intent(in)  :: currTime, nextTime
  integer,         intent(out) :: rc

  !----------------------------------------------------------------------------

end function get_importexport
module subroutine debug_info(tag,fname,chkflag,filesize,timestring)
  character(len=*), intent(in) :: tag
  character(len=*), intent(in) :: fname
  integer,          intent(in) :: filesize
  logical,          intent(in) :: chkflag
  character(len=*), intent(in) :: timestring

  !----------------------------------------------------------------------------

end subroutine debug_info
module subroutine nf90_err(ierr, string)
  integer,          intent(in) :: ierr
  character(len=*), intent(in) :: string
  !----------------------------------------------------------------------------

end subroutine nf90_err
  end interface

end module MOM_cap_outputlog
