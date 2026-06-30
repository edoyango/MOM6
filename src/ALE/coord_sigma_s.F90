submodule (coord_sigma) coord_sigma_s
  implicit none
contains
module procedure init_coord_sigma
  if (associated(CS)) call MOM_error(FATAL, "init_coord_sigma: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
end procedure init_coord_sigma
module procedure end_coord_sigma
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end procedure end_coord_sigma
module procedure set_sigma_params
  if (.not. associated(CS)) call MOM_error(FATAL, "set_sigma_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end procedure set_sigma_params
module procedure build_sigma_column
  integer :: k
  zInterface(CS%nk+1) = -depth
  do k = CS%nk,1,-1
    zInterface(k) = zInterface(k+1) + (totalThickness * CS%coordinateResolution(k))
    ! Adjust interface position to accommodate inflating layers
    ! without disturbing the interface above
    if (zInterface(k) < (zInterface(k+1) + CS%min_thickness)) then
      zInterface(k) = zInterface(k+1) + CS%min_thickness
    endif
  enddo
end procedure build_sigma_column
end submodule coord_sigma_s
