submodule (mom_cap_profiling) mom_cap_profiling_s
  implicit none
contains
  module procedure cap_profiling_init
#ifdef UFS_TRACING
    call ufs_trace_init()
#endif
    return
  end procedure cap_profiling_init
  module procedure cap_profiling
#ifdef UFS_TRACING
    call ufs_trace(component, routine, ph)
#endif
    return
  end procedure cap_profiling
  module procedure cap_profiling_finalize
#ifdef UFS_TRACING
    call ufs_trace_finalize()
#endif
    return
  end procedure cap_profiling_finalize
end submodule mom_cap_profiling_s
