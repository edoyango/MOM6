submodule (MOM_cpu_clock_infra) MOM_cpu_clock_infra_s
  implicit none
contains
module procedure cpu_clock_begin
  call mpp_clock_begin(id)

end procedure cpu_clock_begin
module procedure cpu_clock_end
  call mpp_clock_end(id)

end procedure cpu_clock_end
module procedure cpu_clock_id
  integer :: clock_flags
  clock_flags = clock_flag_default
  if (present(sync)) then
    if (sync) then
      clock_flags = ibset(clock_flags, 0)
    else
      clock_flags = ibclr(clock_flags, 0)
    endif
  endif

  cpu_clock_id = mpp_clock_id(name, flags=clock_flags, grain=grain)
end procedure cpu_clock_id
end submodule MOM_cpu_clock_infra_s
