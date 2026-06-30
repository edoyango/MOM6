submodule (MOM_couplertype_infra) MOM_couplertype_infra_s
  implicit none
contains
module procedure atmos_ocn_coupler_flux
  coupler_index = aof_set_coupler_flux(name, flux_type, implementation,      &
                              param=param, mol_wt=mol_wt, ice_restart_file=ice_restart_file, &
                              ocean_restart_file=ocean_restart_file, &
                              units=units, caller=caller, verbosity=verbosity)

end procedure atmos_ocn_coupler_flux
module procedure CT_spawn_1d_2d
  call coupler_type_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_1d_2d
module procedure CT_spawn_1d_3d
  call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_1d_3d
module procedure CT_spawn_2d_2d
  call coupler_type_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_2d_2d
module procedure CT_spawn_2d_3d
  call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_2d_3d
module procedure CT_spawn_3d_2d
  call coupler_type_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_3d_2d
module procedure CT_spawn_3d_3d
  call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_3d_3d
module procedure CT_copy_data_2d
  call coupler_type_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end procedure CT_copy_data_2d
module procedure CT_copy_data_3d
  call coupler_type_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end procedure CT_copy_data_3d
module procedure CT_copy_data_2d_3d
  call coupler_type_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice, ind3_start, ind3_end)
end procedure CT_copy_data_2d_3d
module procedure CT_increment_data_2d
  call coupler_type_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                         scale_prev=scale_prev)

end procedure CT_increment_data_2d
module procedure CT_increment_data_3d
  call coupler_type_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                                   scale_prev=scale_prev, exclude_flux_type=exclude_flux_type, &
                                   only_flux_type=only_flux_type)

end procedure CT_increment_data_3d
module procedure CT_rescale_data_2d
  call coupler_type_rescale_data(var, scale)

end procedure CT_rescale_data_2d
module procedure CT_rescale_data_3d
  call coupler_type_rescale_data(var, scale)

end procedure CT_rescale_data_3d
module procedure CT_increment_data_2d_3d
  call coupler_type_increment_data(var_in, weights, var, halo_size=halo_size)

end procedure CT_increment_data_2d_3d
module procedure CT_redistribute_data_2d
  call coupler_type_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end procedure CT_redistribute_data_2d
module procedure CT_redistribute_data_3d
  call coupler_type_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end procedure CT_redistribute_data_3d
module procedure CT_extract_data
  call coupler_type_extract_data(var_in, bc_index, field_index, array_out, scale_factor, halo_size, idim, jdim)

end procedure CT_extract_data
module procedure CT_set_data
  call coupler_type_set_data(array_in, bc_index, field_index, var, scale_factor, halo_size, idim, jdim)

end procedure CT_set_data
module procedure CT_data_override
  call coupler_type_data_override(gridname, var, time)
end procedure CT_data_override
module procedure CT_set_diags
  call coupler_type_set_diags(var, diag_name, axes, time)

end procedure CT_set_diags
module procedure CT_send_data
  call coupler_type_send_data(var, Time)
end procedure CT_send_data
module procedure CT_write_chksums_2d
  call coupler_type_write_chksums(var, outunit, name_lead)

end procedure CT_write_chksums_2d
module procedure CT_write_chksums_3d
  call coupler_type_write_chksums(var, outunit, name_lead)

end procedure CT_write_chksums_3d
module procedure CT_initialized_1d
  CT_initialized_1d = coupler_type_initialized(var)
end procedure CT_initialized_1d
module procedure CT_initialized_2d
  CT_initialized_2d = coupler_type_initialized(var)
end procedure CT_initialized_2d
module procedure CT_initialized_3d
  CT_initialized_3d = coupler_type_initialized(var)
end procedure CT_initialized_3d
module procedure CT_destructor_1d
  call coupler_type_destructor(var)

end procedure CT_destructor_1d
module procedure CT_destructor_2d
  call coupler_type_destructor(var)

end procedure CT_destructor_2d
end submodule MOM_couplertype_infra_s
