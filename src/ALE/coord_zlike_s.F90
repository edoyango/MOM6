submodule (coord_zlike) coord_zlike_s
  implicit none
contains
module procedure init_coord_zlike
  if (associated(CS)) call MOM_error(FATAL, "init_coord_zlike: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
end procedure init_coord_zlike
module procedure end_coord_zlike
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end procedure end_coord_zlike
module procedure set_zlike_params
  if (.not. associated(CS)) call MOM_error(FATAL, "set_zlike_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end procedure set_zlike_params
module procedure build_zstar_column
  real :: eta   ! Free surface height [Z ~> m] or [H ~> m or kg m-2]
  real :: stretching ! A stretching factor for the coordinate [nondim]
  real :: dh, min_thickness, z0_top, z_star, z_scale ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]
  integer :: k
  logical :: new_zstar_def
  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  new_zstar_def = .false.
  min_thickness = min( CS%min_thickness, total_thickness/real(CS%nk) )
  z0_top = 0.
  if (present(z_rigid_top)) then
    z0_top = z_rigid_top
    new_zstar_def = .true.
  endif

  ! Position of free-surface (or the rigid top, for which eta ~ z0_top)
  eta = total_thickness - depth
  if (present(eta_orig)) eta = eta_orig

  ! Conventional z* coordinate:
  !   z* = (z-eta) / stretching   where stretching = (H+eta)/H
  !   z = eta + stretching * z*
  ! The above gives z*(z=eta) = 0, z*(z=-H) = -H.
  ! With a rigid top boundary at eta = z0_top then
  !   z* = z0 + (z-eta) / stretching   where stretching = (H+eta)/(H+z0)
  !   z = eta + stretching * (z*-z0) * stretching
  stretching = total_thickness / ( depth + z0_top )

  if (new_zstar_def) then
    ! z_star is the notional z* coordinate in absence of upper/lower topography
    z_star = 0. ! z*=0 at the free-surface
    zInterface(1) = eta ! The actual position of the top of the column
    do k = 2,CS%nk
      z_star = z_star - CS%coordinateResolution(k-1)*z_scale
      ! This ensures that z is below a rigid upper surface (ice shelf bottom)
      zInterface(k) = min( eta + stretching * ( z_star - z0_top ), z0_top )
      ! This ensures that the layer in inflated
      zInterface(k) = min( zInterface(k), zInterface(k-1) - min_thickness )
      ! This ensures that z is above or at the topography
      zInterface(k) = max( zInterface(k), -depth + real(CS%nk+1-k) * min_thickness )
    enddo
    zInterface(CS%nk+1) = -depth

  else
    ! Integrate down from the top for a notional new grid, ignoring topography
    ! The starting position is offset by z0_top which, if z0_top<0, will place
    ! interfaces above the rigid boundary.
    zInterface(1) = eta
    do k = 1,CS%nk
      dh = stretching * CS%coordinateResolution(k)*z_scale ! Notional grid spacing
      zInterface(k+1) = zInterface(k) - dh
    enddo

    ! Integrating up from the bottom adjusting interface position to accommodate
    ! inflating layers without disturbing the interface above
    zInterface(CS%nk+1) = -depth
    do k = CS%nk,1,-1
      if ( zInterface(k) < (zInterface(k+1) + min_thickness) ) then
        zInterface(k) = zInterface(k+1) + min_thickness
      endif
    enddo
  endif

end procedure build_zstar_column
end submodule coord_zlike_s
