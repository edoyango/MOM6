!> Contains wrapper routines that call the ufs tracing routines
module mom_cap_profiling

#ifdef UFS_TRACING
  use ufs_trace_mod, only: ufs_trace_init, ufs_trace, ufs_trace_finalize
#endif

  implicit none

  private

  public cap_profiling_init
  public cap_profiling
  public cap_profiling_finalize


  interface
  module subroutine cap_profiling_init()
  end subroutine cap_profiling_init
  module subroutine cap_profiling(component, routine, ph)
    character(len=*), intent(in) :: component !< Name of the component, 'mom'
    character(len=*), intent(in) :: routine   !< Name of the profiled subroutine
    character(len=*), intent(in) :: ph        !< Duration event phase type. 'B' or 'E' for begin/end
  end subroutine cap_profiling
  module subroutine cap_profiling_finalize()
  end subroutine cap_profiling_finalize
  end interface

end module mom_cap_profiling
