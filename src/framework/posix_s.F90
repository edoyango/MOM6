submodule (posix) posix_s
  implicit none
contains
module procedure chmod
  integer(kind=c_int) :: mode_c
  integer(kind=c_int) :: rc_c
  mode_c = int(mode, kind=c_int)
  rc_c = chmod_posix(path//c_null_char, mode_c)
  rc = int(rc_c)
end procedure chmod
module procedure mkdir
  integer(kind=c_int) :: mode_c
  integer(kind=c_int) :: rc_c
  mode_c = int(mode, kind=c_int)
  rc_c = mkdir_posix(path//c_null_char, mode_c)
  rc = int(rc_c)
end procedure mkdir
module procedure stat
  integer(kind=c_int) :: rc_c
  rc_c = stat_posix(path//c_null_char, buf)

  rc = int(rc_c)
end procedure stat
module procedure signal
  integer(kind=c_int) :: sig_c
  type(c_funptr) :: handle_c
  sig_c = int(sig, kind=c_int)
  handle_c = signal_posix(sig_c, c_funloc(func))
  call c_f_procpointer(handle_c, handle)
end procedure signal
module procedure kill
  integer(kind=c_int) :: pid_c, sig_c, rc_c
  pid_c = int(pid, kind=c_int)
  sig_c = int(sig, kind=c_int)
  rc_c = kill_posix(pid_c, sig_c)
  rc = int(rc_c)
end procedure kill
module procedure getpid
  integer(kind=c_long) :: pid_c
  pid_c = getpid_posix()
  pid = int(pid_c)
end procedure getpid
module procedure getppid
  integer(kind=c_long) :: pid_c
  pid_c = getppid_posix()
  pid = int(pid_c)
end procedure getppid
module procedure sleep
  integer(kind=c_int) :: seconds_c
  integer(kind=c_int) :: rc_c
  seconds_c = int(seconds, kind=c_int)
  rc_c = sleep_posix(seconds_c)
  rc = int(rc_c)
end procedure sleep
module procedure longjmp
  integer(kind=c_int) :: val_c
  val_c = int(val, kind=c_int)
  call longjmp_posix(env, val_c)
end procedure longjmp
module procedure siglongjmp
  integer(kind=c_int) :: val_c
  val_c = int(val, kind=c_int)
  call siglongjmp_posix(env, val_c)
end procedure siglongjmp
module procedure setjmp_missing
  print '(a)', 'ERROR: setjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSETJMP_NAME=\"<symbol name>\".'
  error stop

  ! NOTE: compilers may expect a return value, even if it is unreachable
  read env%state
  rc = -1
end procedure setjmp_missing
module procedure longjmp_missing
  print '(a)', 'ERROR: longjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DLONGJMP_NAME=\"<symbol name>\".'
  error stop

  read env%state
  read char(val)
end procedure longjmp_missing
module procedure sigsetjmp_missing
  print '(a)', 'ERROR: sigsetjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSIGSETJMP_NAME=\"<symbol name>\".'
  error stop

  ! NOTE: compilers may expect a return value, even if it is unreachable
  read env%state
  read char(savesigs)
  rc = -1
end procedure sigsetjmp_missing
module procedure siglongjmp_missing
  print '(a)', 'ERROR: siglongjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSIGLONGJMP_NAME=\"<symbol name>\".'
  read env%state
  read char(val)
  error stop
end procedure siglongjmp_missing
end submodule posix_s
