submodule (MOM_data_override_infra) MOM_data_override_infra_s
  implicit none
contains
module procedure impose_data_init
  if (present(MOM_domain_in)) then
    call data_override_init(Ocean_domain_in=MOM_domain_in%mpp_domain, Ice_domain_in=Ice_domain_in)
  else
    call data_override_init(Ocean_domain_in=Ocean_domain_in, Ice_domain_in=Ice_domain_in)
  endif
end procedure impose_data_init
module procedure data_override_MD
  logical :: overridden, is_ocean
  integer :: i, j, is, ie, js, je
  overridden = .false.
  is_ocean = .true. ; if (present(is_ice)) is_ocean = .not.is_ice
  if (is_ocean) then
    call data_override('OCN', fieldname, data_2D, time, override=overridden)
  else
    call data_override('ICE', fieldname, data_2D, time, override=overridden)
  endif

  if (overridden .and. present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data in the computational domain if the data override has occurred.
    call get_simple_array_i_ind(domain, size(data_2D,1), is, ie)
    call get_simple_array_j_ind(domain, size(data_2D,2), js, je)
    do j=js,je ; do i=is,ie
      data_2D(i,j) = scale*data_2D(i,j)
    enddo ; enddo
  endif ; endif

  if (present(override)) override = overridden

end procedure data_override_MD
module procedure data_override_2d
  call data_override(gridname, fieldname, data_2D, time, override)

end procedure data_override_2d
module procedure impose_data_unset_domains
  call data_override_unset_domains(unset_Ocean=unset_Ocean, unset_Ice=unset_Ice, &
                                   must_be_set=must_be_set)
end procedure impose_data_unset_domains
end submodule MOM_data_override_infra_s
