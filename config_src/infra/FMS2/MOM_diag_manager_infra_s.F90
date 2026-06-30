submodule (MOM_diag_manager_infra) MOM_diag_manager_infra_s
  implicit none
contains
module procedure MOM_diag_axis_init
  integer :: coarsening ! The degree of grid coarsening, this is the index of an array of coarsening levels
  if (present(null_axis)) then ; if (null_axis) then
    ! Return the special null axis id for scalars
    MOM_diag_axis_init = null_axis_id
    return
  endif ; endif

  if (present(MOM_domain)) then
    coarsening = 0 ; if (present(coarsen)) coarsening = coarsen
    if (coarsening == 0) then
      MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain, domain_position=position)
    else
      MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain_d(coarsening), domain_position=position)
    endif
  else
    if (present(coarsen)) then ; if (coarsen /= 1) then
      call MOM_error(FATAL, "diag_axis_init does not support grid coarsening without a MOM_domain.")
    endif ; endif
    MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
            direction=direction, set_name=set_name, edges=edges)
  endif

end procedure MOM_diag_axis_init
module procedure get_MOM_diag_axis_name
  call fms_get_diag_axis_name(id, name)

end procedure get_MOM_diag_axis_name
module procedure get_MOM_diag_field_id
  get_MOM_diag_field_id = -1
  get_MOM_diag_field_id = get_diag_field_id_fms(module_name, field_name)

end procedure get_MOM_diag_field_id
module procedure MOM_diag_manager_init
  call FMS_diag_manager_init(diag_model_subset, time_init, err_msg)

end procedure MOM_diag_manager_init
module procedure MOM_diag_manager_end
  call FMS_diag_manager_end(time)

end procedure MOM_diag_manager_end
module procedure register_diag_field_infra_scalar
  register_diag_field_infra_scalar = register_diag_field_fms(module_name, field_name, init_time, &
        long_name, units, missing_value, range, standard_name, do_not_log, err_msg, area, volume)

end procedure register_diag_field_infra_scalar
module procedure register_diag_field_infra_array
  register_diag_field_infra_array = register_diag_field_fms(module_name, field_name, axes, init_time, &
        long_name, units, missing_value, range, mask_variant, standard_name, verbose, do_not_log, &
        err_msg, interp_method, tile_count, area, volume)

end procedure register_diag_field_infra_array
module procedure register_static_field_infra
  if(present(missing_value) .or. present(range)) then
    register_static_field_infra = register_static_field_fms(module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant=mask_variant, standard_name=standard_name, dynamic=.false.,&
       do_not_log=do_not_log, interp_method=interp_method,tile_count=tile_count, area=area, volume=volume)
  else
    register_static_field_infra = register_static_field_fms(module_name, field_name, axes, long_name, units,&
       &  mask_variant=mask_variant, standard_name=standard_name, dynamic=.false.,do_not_log=do_not_log, &
       interp_method=interp_method,tile_count=tile_count, area=area, volume=volume)
  endif
end procedure register_static_field_infra
module procedure send_data_infra_0d
  send_data_infra_0d = send_data_fms(diag_field_id, field, time, err_msg)
end procedure send_data_infra_0d
module procedure send_data_infra_1d
  if(present(rmask) .or. present(weight)) then
   if(present(rmask) .and. present(weight)) then
  send_data_infra_1d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, mask=mask, rmask=rmask, ie_in=ie_in,&
                                     weight=weight, err_msg=err_msg)
   elseif(present(rmask)) then
  send_data_infra_1d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, mask=mask, rmask=rmask, ie_in=ie_in,&
                                     err_msg=err_msg)
   elseif(present(weight)) then
  send_data_infra_1d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, ie_in=ie_in, weight=weight,&
                                     err_msg=err_msg)
   endif
  else
  send_data_infra_1d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, ie_in=ie_in, err_msg=err_msg)
  endif

end procedure send_data_infra_1d
module procedure send_data_infra_2d
  if(present(rmask) .or. present(weight)) then
   if(present(rmask) .and. present(weight)) then
    send_data_infra_2d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, js_in=js_in, mask=mask, &
                                rmask=rmask, ie_in=ie_in, je_in=je_in, weight=weight, err_msg=err_msg)
   elseif(present(rmask)) then
    send_data_infra_2d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, js_in=js_in, mask=mask, &
                                rmask=rmask, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
   elseif(present(weight)) then
    send_data_infra_2d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, js_in=js_in, mask=mask, &
                                ie_in=ie_in, je_in=je_in, weight=weight, err_msg=err_msg)
   endif
  else
    send_data_infra_2d = send_data_fms(diag_field_id, field, time=time, is_in=is_in, js_in=js_in, mask=mask, &
                                ie_in=ie_in, je_in=je_in, err_msg=err_msg)
  endif
end procedure send_data_infra_2d
module procedure send_data_infra_3d
  send_data_infra_3d = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, &
                               rmask, ie_in, je_in, ke_in, weight, err_msg)

end procedure send_data_infra_3d
module procedure send_data_infra_2d_r8
  send_data_infra_2d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, mask, &
                                   rmask, ie_in, je_in, weight, err_msg)

end procedure send_data_infra_2d_r8
module procedure send_data_infra_3d_r8
  send_data_infra_3d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, rmask, &
                                ie_in, je_in, ke_in, weight, err_msg)

end procedure send_data_infra_3d_r8
module procedure MOM_diag_field_add_attribute_scalar_r
  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end procedure MOM_diag_field_add_attribute_scalar_r
module procedure MOM_diag_field_add_attribute_scalar_i
  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end procedure MOM_diag_field_add_attribute_scalar_i
module procedure MOM_diag_field_add_attribute_scalar_c
  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end procedure MOM_diag_field_add_attribute_scalar_c
module procedure MOM_diag_field_add_attribute_r1d
  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end procedure MOM_diag_field_add_attribute_r1d
module procedure MOM_diag_field_add_attribute_i1d
  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end procedure MOM_diag_field_add_attribute_i1d
module procedure diag_send_complete_infra
  call diag_send_complete (set_time(0))
end procedure diag_send_complete_infra
module procedure diag_manager_set_time_end_infra
  call diag_manager_set_time_end(time)
end procedure diag_manager_set_time_end_infra
end submodule MOM_diag_manager_infra_s
