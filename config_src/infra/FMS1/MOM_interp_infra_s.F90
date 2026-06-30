submodule (MOM_interp_infra) MOM_interp_infra_s
  implicit none
contains
module procedure horizontal_interp_init
  call horiz_interp_init()
end procedure horizontal_interp_init
module procedure time_interp_extern_init
  call time_interp_external_init()
end procedure time_interp_extern_init
module procedure horiz_interp_from_weights_field2d
  call horiz_interp(Interp, data_in, data_out, verbose, &
                    mask_in, mask_out, missing_value, missing_permit, &
                    err_msg, new_missing_handle=.true. )

end procedure horiz_interp_from_weights_field2d
module procedure horiz_interp_from_weights_field3d
  call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, &
                    missing_value, missing_permit, err_msg)

end procedure horiz_interp_from_weights_field3d
module procedure build_horiz_interp_weights_2d_to_2d
  call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, &
                        verbose, interp_method, num_nbrs, max_dist, &
                        src_modulo, mask_in, mask_out, &
                        is_latlon_in, is_latlon_out)

end procedure build_horiz_interp_weights_2d_to_2d
module procedure get_axis_data
  call mpp_get_axis_data( axis, dat )
end procedure get_axis_data
module procedure get_extern_field_size
  get_extern_field_size = get_external_field_size(index)

end procedure get_extern_field_size
module procedure get_extern_field_axes
  get_extern_field_axes = get_external_field_axes(index)
end procedure get_extern_field_axes
module procedure get_extern_field_missing
  get_extern_field_missing = get_external_field_missing(index)

end procedure get_extern_field_missing
module procedure get_external_field_info
  if (present(size)) then
    size(:) = get_extern_field_size(field%id)
  endif

  if (present(axes)) then
    axes(:) = get_extern_field_axes(field%id)
  endif

  if (present(missing)) then
    missing = get_extern_field_missing(field%id)
  endif
end procedure get_external_field_info
module procedure time_interp_extern_0d
  call time_interp_external(field%id, time, data_in, verbose=verbose)
end procedure time_interp_extern_0d
module procedure time_interp_extern_2d
  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end procedure time_interp_extern_2d
module procedure time_interp_extern_3d
  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end procedure time_interp_extern_3d
module procedure init_extern_field
  if (present(MOM_Domain)) then
    field%id = init_external_field(file, fieldname, domain=MOM_domain%mpp_domain, &
             verbose=verbose, threading=threading, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency)
  else
    field%id = init_external_field(file, fieldname, domain=domain, &
             verbose=verbose, threading=threading, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency)
  endif
end procedure init_extern_field
end submodule MOM_interp_infra_s
