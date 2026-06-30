! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A module to monitor the overall CPU time used by MOM6 and project when to stop the model
module MOM_write_cputime

use MOM_coms,          only : sum_across_PEs, num_pes
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_io,            only : open_ASCII_file, close_file, APPEND_FILE, WRITEONLY_FILE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_time_manager,  only : time_type, get_time, operator(>)

implicit none ; private

public write_cputime, MOM_write_cputime_init, MOM_write_cputime_end, write_cputime_start_clock

!-----------------------------------------------------------------------

integer :: CLOCKS_PER_SEC = 1000 !< The number of clock cycles per second, used by the system clock
integer :: MAX_TICKS      = 1000 !< The number of ticks per second, used by the system clock

!> A control structure that regulates the writing of CPU time
type, public :: write_cputime_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real :: maxcpu                !<   The maximum amount of CPU time per processor
                                !! for which MOM should run before saving a restart
                                !! file and quitting with a return value that
                                !! indicates that further execution is required to
                                !! complete the simulation [wall-clock seconds].
  type(time_type) :: Start_time !< The start time of the simulation.
                                !! Start_time is set in MOM_initialization.F90
  real :: startup_cputime       !< The CPU time used in the startup phase of the model [clock_cycles].
  real :: prev_cputime = 0.0    !< The last measured CPU time [clock_cycles].
  real :: dn_dcpu_min = -1.0    !< The minimum derivative of timestep with CPU time [steps clock_cycles-1].
  real :: cputime2 = 0.0        !< The accumulated CPU time [clock_cycles].
  integer :: previous_calls = 0 !< The number of times write_CPUtime has been called.
  integer :: prev_n = 0         !< The value of n from the last call.
  integer :: fileCPU_ascii= -1  !< The unit number of the CPU time file.
  character(len=200) :: CPUfile !< The name of the CPU time file.
end type write_cputime_CS


  interface
module subroutine write_cputime_start_clock(CS)
  type(write_cputime_CS), pointer :: CS !< The control structure set up by a previous
                                        !! call to MOM_write_cputime_init.
end subroutine write_cputime_start_clock
module subroutine MOM_write_cputime_init(param_file, directory, Input_start_time, CS)
  type(param_file_type),  intent(in) :: param_file !< A structure to parse for run-time parameters
  character(len=*),       intent(in) :: directory  !< The directory where the CPU time file goes.
  type(time_type),        intent(in) :: Input_start_time !< The start model time of the simulation.
  type(write_cputime_CS), pointer    :: CS         !< A pointer that may be set to point to the
                                                   !! control structure for this module.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine MOM_write_cputime_init
module subroutine MOM_write_cputime_end(CS)
  type(write_cputime_CS), pointer    :: CS    !< The control structure set up by a previous
                                              !! call to MOM_write_cputime_init.

end subroutine MOM_write_cputime_end
module subroutine write_cputime(day, n, CS, nmax, call_end)
  type(time_type),        intent(inout) :: day  !< The current model time.
  integer,                intent(in)    :: n    !< The time step number of the current execution.
  type(write_cputime_CS), pointer       :: CS   !< The control structure set up by a previous
                                                !! call to MOM_write_cputime_init.
  integer,      optional, intent(inout) :: nmax !< The number of iterations after which to stop so
                                                !! that the simulation will not run out of CPU time.
  logical,      optional, intent(in)    :: call_end !< If true, also call MOM_write_cputime_end.

  ! Local variables
                           ! this subroutine [clock_cycles]

end subroutine write_cputime
  end interface

end module MOM_write_cputime
