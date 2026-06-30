submodule (MOM_error_infra) MOM_error_infra_s
  implicit none
contains
module procedure MOM_err
  call mpp_error(severity, message)
end procedure MOM_err
module procedure stdout
  stdout = mpp_stdout()
end procedure stdout
module procedure stdlog
  stdlog = mpp_stdlog()
end procedure stdlog
module procedure is_root_pe
  is_root_pe = .false.
  if (mpp_pe() == mpp_root_pe()) is_root_pe = .true.
end procedure is_root_pe
end submodule MOM_error_infra_s
