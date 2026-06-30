submodule (MOM_coupler_types) MOM_coupler_types_s
  implicit none
contains
module procedure CT_spawn_1d_2d
  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_1d_2d
module procedure CT_spawn_1d_3d
  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_1d_3d
module procedure CT_spawn_2d_2d
  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_2d_2d
module procedure CT_spawn_2d_3d
  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_2d_3d
module procedure CT_spawn_3d_2d
  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_3d_2d
module procedure CT_spawn_3d_3d
  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end procedure CT_spawn_3d_3d
module procedure CT_copy_data_2d
  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end procedure CT_copy_data_2d
module procedure CT_copy_data_3d
  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end procedure CT_copy_data_3d
module procedure CT_copy_data_2d_3d
  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice, ind3_start, ind3_end)
end procedure CT_copy_data_2d_3d
module procedure CT_increment_data_2d
  call CT_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                         scale_prev=scale_prev)

end procedure CT_increment_data_2d
module procedure CT_increment_data_3d
  call CT_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                         scale_prev=scale_prev, exclude_flux_type=exclude_flux_type, &
                         only_flux_type=only_flux_type)

end procedure CT_increment_data_3d
module procedure CT_increment_data_2d_3d
  call CT_increment_data(var_in, weights, var, halo_size=halo_size)

end procedure CT_increment_data_2d_3d
module procedure CT_rescale_data_2d
  call CT_rescale_data(var, scale)

end procedure CT_rescale_data_2d
module procedure CT_rescale_data_3d
  call CT_rescale_data(var, scale)

end procedure CT_rescale_data_3d
module procedure CT_redistribute_data_2d
  call CT_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end procedure CT_redistribute_data_2d
module procedure CT_redistribute_data_3d
  call CT_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end procedure CT_redistribute_data_3d
module procedure coupler_type_data_override
  call CT_data_override(gridname, var, time)
end procedure coupler_type_data_override
module procedure extract_coupler_type_data
  real, allocatable :: array_unrot(:,:)  ! Array on the unrotated grid in arbitrary units [A]
  integer :: q_turns ! The number of quarter turns through which array_out is to be rotated
  integer :: index
  index = ind_flux ; if (present(field_index)) index = field_index
  q_turns = 0 ; if (present(turns)) q_turns = modulo(turns, 4)

  ! The case with non-trivial grid rotation is complicated by the fact that the data fields
  ! in the coupler_2d_bc_type are never rotated, so they need to be handled separately.
  if (q_turns == 0) then
    call CT_extract_data(var_in, bc_index, index, array_out, &
        scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)
  elseif (present(idim) .and. present(jdim)) then
    call allocate_rotated_array(array_out, [1,1], -q_turns, array_unrot)

    if (modulo(q_turns, 2) /= 0) then
      call CT_extract_data(var_in, bc_index, index, array_unrot, &
          idim=jdim, jdim=idim, scale_factor=scale_factor, halo_size=halo_size)
    else
      call CT_extract_data(var_in, bc_index, index, array_unrot, &
          idim=idim, jdim=jdim, scale_factor=scale_factor, halo_size=halo_size)
    endif

    call rotate_array(array_unrot, q_turns, array_out)
    deallocate(array_unrot)
  else
    call allocate_rotated_array(array_out, [1,1], -q_turns, array_unrot)
    call CT_extract_data(var_in, bc_index, index, array_unrot, &
        scale_factor=scale_factor, halo_size=halo_size)
    call rotate_array(array_unrot, q_turns, array_out)
    deallocate(array_unrot)
  endif

end procedure extract_coupler_type_data
module procedure set_coupler_type_data
  real, allocatable :: array_unrot(:,:)  ! Array on the unrotated grid in the same arbitrary units
  integer :: subfield ! An integer indicating which field to set.
  integer :: q_turns ! The number of quarter turns through which array_in is rotated
  q_turns = 0 ; if (present(turns)) q_turns = modulo(turns, 4)

  subfield = ind_csurf
  if (present(solubility)) then ; if (solubility) subfield = ind_alpha ; endif
  if (present(field_index)) subfield = field_index

  ! The case with non-trivial grid rotation is complicated by the fact that the data fields
  ! in the coupler_2d_bc_type are never rotated, so they need to be handled separately.
  if (q_turns == 0) then
    call CT_set_data(array_in, bc_index, subfield, var, &
                     scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)
  elseif (present(idim) .and. present(jdim)) then
    call allocate_rotated_array(array_in, [1,1], -q_turns, array_unrot)
    call rotate_array(array_in, -q_turns, array_unrot)

    if (modulo(q_turns, 2) /= 0) then
      call CT_set_data(array_unrot, bc_index, subfield, var, &
          idim=jdim, jdim=idim, &
          scale_factor=scale_factor, halo_size=halo_size)
    else
      call CT_set_data(array_unrot, bc_index, subfield, var, &
          idim=idim, jdim=jdim, &
          scale_factor=scale_factor, halo_size=halo_size)
    endif

    deallocate(array_unrot)
  else
    call allocate_rotated_array(array_in, [1,1], -q_turns, array_unrot)
    call rotate_array(array_in, -q_turns, array_unrot)
    call CT_set_data(array_in, bc_index, subfield, var, &
                     scale_factor=scale_factor, halo_size=halo_size)
    deallocate(array_unrot)
  endif

end procedure set_coupler_type_data
module procedure coupler_type_set_diags
  call CT_set_diags(var, diag_name, axes, time)

end procedure coupler_type_set_diags
module procedure coupler_type_send_data
  call CT_send_data(var, Time)
end procedure coupler_type_send_data
module procedure CT_write_chksums_2d
  call CT_write_chksums(var, outunit, name_lead)

end procedure CT_write_chksums_2d
module procedure CT_write_chksums_3d
  call CT_write_chksums(var, outunit, name_lead)

end procedure CT_write_chksums_3d
module procedure CT_initialized_1d
  CT_initialized_1d = CT_initialized(var)
end procedure CT_initialized_1d
module procedure CT_initialized_2d
  CT_initialized_2d = CT_initialized(var)
end procedure CT_initialized_2d
module procedure CT_initialized_3d
  CT_initialized_3d = CT_initialized(var)
end procedure CT_initialized_3d
module procedure CT_destructor_1d
  call CT_destructor(var)

end procedure CT_destructor_1d
module procedure CT_destructor_2d
  call CT_destructor(var)

end procedure CT_destructor_2d
end submodule MOM_coupler_types_s
