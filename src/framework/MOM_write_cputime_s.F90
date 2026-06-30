submodule (MOM_write_cputime) MOM_write_cputime_s
  implicit none
contains
module procedure write_cputime_start_clock
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  if (.not.associated(CS)) allocate(CS)

  call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
  CS%prev_cputime = new_cputime
end procedure write_cputime_start_clock
module procedure MOM_write_cputime_init
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
# include "version_variable.h"
  character(len=40)  :: mdl = 'MOM_write_cputime'  ! This module's name.
  logical :: all_default   ! If true, all parameters are using their default values.
  if (.not.associated(CS)) then
    allocate(CS)
    call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
    CS%prev_cputime = new_cputime
  endif

  CS%initialized = .true.

  ! Read all relevant parameters and write them to the model log.

  ! Determine whether all parameters are set to their default values.
  call get_param(param_file, mdl, "MAXCPU", CS%maxcpu, units="wall-clock seconds", default=-1.0, do_not_log=.true.)
  call get_param(param_file, mdl, "CPU_TIME_FILE", CS%CPUfile, default="CPU_stats", do_not_log=.true.)
  all_default = (CS%maxcpu == -1.0) .and. (trim(CS%CPUfile) == trim("CPU_stats"))

  call log_version(param_file, mdl, version, "", all_default=all_default)
  call get_param(param_file, mdl, "MAXCPU", CS%maxcpu, &
                 "The maximum amount of cpu time per processor for which "//&
                 "MOM should run before saving a restart file and "//&
                 "quitting with a return value that indicates that a "//&
                 "further run is required to complete the simulation. "//&
                 "If automatic restarts are not desired, use a negative "//&
                 "value for MAXCPU.  MAXCPU has units of wall-clock "//&
                 "seconds, so the actual CPU time used is larger by a "//&
                 "factor of the number of processors used.", &
                 units="wall-clock seconds", default=-1.0)
  call get_param(param_file, mdl, "CPU_TIME_FILE", CS%CPUfile, &
                 "The file into which CPU time is written.",default="CPU_stats")
  CS%CPUfile = trim(directory)//trim(CS%CPUfile)
  call log_param(param_file, mdl, "directory/CPU_TIME_FILE", CS%CPUfile)
#ifdef STATSLABEL
  CS%CPUfile = trim(CS%CPUfile)//"."//trim(adjustl(STATSLABEL))
#endif

  CS%Start_time = Input_start_time

end procedure MOM_write_cputime_init
module procedure MOM_write_cputime_end
  if (.not.associated(CS)) return

  ! Flush and close the output files.
  if (is_root_pe() .and. CS%fileCPU_ascii > 0) then
    flush(CS%fileCPU_ascii)
    call close_file(CS%fileCPU_ascii)
  endif

  deallocate(CS)

end procedure MOM_write_cputime_end
module procedure write_cputime
  real    :: d_cputime     ! The change in CPU time since the last call
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK [clock_cycles]
  real    :: reday         ! The time in days, including fractional days [days]
  integer :: start_of_day  ! The number of seconds since the start of the day
  integer :: num_days      ! The number of days in the time
  if (.not.associated(CS)) call MOM_error(FATAL, &
         "write_energy: Module must be initialized before it is used.")

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "write_cputime: Module must be initialized before it is used.")

  call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
!   The following lines extract useful information even if the clock has rolled
! over, assuming a 32-bit SYSTEM_CLOCK.  With more bits, rollover is essentially
! impossible. Negative fluctuations of less than 10 seconds are not interpreted
! as the clock rolling over.  This should be unnecessary but is sometimes needed
! on the GFDL SGI/O3k.
  if (new_cputime < CS%prev_cputime-(10.0*CLOCKS_PER_SEC)) then
    d_cputime = new_cputime - CS%prev_cputime + MAX_TICKS
  else
    d_cputime = new_cputime - CS%prev_cputime
  endif

  call sum_across_PEs(d_cputime)
  if (CS%previous_calls == 0) CS%startup_cputime = d_cputime

  CS%cputime2 = CS%cputime2 + d_cputime

  if ((CS%previous_calls >= 1) .and. (CS%maxcpu > 0.0)) then
    ! Determine the slowest rate at which time steps are executed.
    if ((n > CS%prev_n) .and. (d_cputime > 0.0) .and. &
        ((CS%dn_dcpu_min*d_cputime < (n - CS%prev_n)) .or. &
         (CS%dn_dcpu_min < 0.0))) &
      CS%dn_dcpu_min = (n - CS%prev_n) / d_cputime
    if (present(nmax) .and. (CS%dn_dcpu_min >= 0.0)) then
      ! Have the model stop itself after 95% of the CPU time has been used.
      nmax = n + INT( CS%dn_dcpu_min * &
          (0.95*CS%maxcpu * REAL(num_pes())*CLOCKS_PER_SEC - &
           (CS%startup_cputime + CS%cputime2)) )
!     write(mesg,*) "Resetting nmax to ",nmax," at day",reday
!     call MOM_mesg(mesg)
    endif
  endif
  CS%prev_cputime = new_cputime ; CS%prev_n = n

  call get_time(day, start_of_day, num_days)
  reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)

  !  Reopen or create a text output file.
  if ((CS%previous_calls == 0) .and. (is_root_pe())) then
    if (day > CS%Start_time) then
      call open_ASCII_file(CS%fileCPU_ascii, trim(CS%CPUfile), action=APPEND_FILE)
    else
      call open_ASCII_file(CS%fileCPU_ascii, trim(CS%CPUfile), action=WRITEONLY_FILE)
    endif
  endif

  if (is_root_pe()) then
    if (CS%previous_calls == 0) then
      write(CS%fileCPU_ascii, &
        '("Startup CPU time: ", F12.3, " sec summed across", I5, " PEs.")') &
                            (CS%startup_cputime / CLOCKS_PER_SEC), num_pes()
      write(CS%fileCPU_ascii,*)"        Day, Step number,     CPU time, CPU time change"
    endif
    write(CS%fileCPU_ascii,'(F12.3,", ",I11,", ",F12.3,", ",F12.3)') &
           reday, n, (CS%cputime2 / real(CLOCKS_PER_SEC)), &
           d_cputime / real(CLOCKS_PER_SEC)

    flush(CS%fileCPU_ascii)
  endif
  CS%previous_calls = CS%previous_calls + 1

  if (present(call_end)) then
    if (call_end) call MOM_write_cputime_end(CS)
  endif

end procedure write_cputime
end submodule MOM_write_cputime_s
