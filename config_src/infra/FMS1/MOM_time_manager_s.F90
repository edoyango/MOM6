submodule (MOM_time_manager) MOM_time_manager_s
  implicit none
contains
module procedure real_to_time
  real             :: x  ! The time in real seconds [s]
  real             :: real_subsecond_remainder ! The fractional seconds from time_in [s]
  integer          :: seconds, days, ticks
  x = time_in ; if (present(unscale)) x = unscale*time_in

  days = floor(x/86400.)
  seconds = floor(x - 86400.*days)
  real_subsecond_remainder = x - (days*86400. + seconds)
  ticks = nint(real_subsecond_remainder * get_ticks_per_second())

  real_to_time = set_time(seconds=seconds, days=days, ticks=ticks, err_msg=err_msg)

end procedure real_to_time
module procedure time_to_real
  time_to_real = time_type_to_real(time)
  if (present(scale)) time_to_real = scale * time_to_real

end procedure time_to_real
module procedure time_minus_signed
  real :: abs_diff ! The absolute value of the difference in times [s] or [T ~> s]
  abs_diff = time_to_real(time_a - time_b, scale)

  ! Add the sign back by comparing time_a and time_b
  time_minus_signed = merge(abs_diff, -abs_diff, time_a >= time_b)

end procedure time_minus_signed
end submodule MOM_time_manager_s
