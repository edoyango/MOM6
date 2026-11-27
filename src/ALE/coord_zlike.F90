! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
#include "do_concurrent_compat.h"

!> Regrid columns for a z-like coordinate (z-star, z-level)
module coord_zlike

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

!> Control structure containing required parameters for a z-like coordinate
type, public :: zlike_CS ; private

  !> Number of levels to be generated
  integer :: nk

  !> Minimum thickness allowed for layers, in the same thickness units (perhaps [H ~> m or kg m-2])
  !! that will be used in all subsequent calls to build_zstar_column with this structure.
  real :: min_thickness

  !> Target coordinate resolution, usually in [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution
end type zlike_CS

public init_coord_zlike, set_zlike_params, build_zstar_column, end_coord_zlike
public build_zstar_block

contains

!> Initialise a zlike_CS with pointers to parameters
subroutine init_coord_zlike(CS, nk, coordinateResolution)
  type(zlike_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of levels in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Target coordinate resolution [Z ~> m]

  if (associated(CS)) call MOM_error(FATAL, "init_coord_zlike: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
end subroutine init_coord_zlike

!> Deallocates the zlike control structure
subroutine end_coord_zlike(CS)
  type(zlike_CS), pointer :: CS !< Coordinate control structure

  ! Nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end subroutine end_coord_zlike

!> Set parameters in the zlike structure
subroutine set_zlike_params(CS, min_thickness)
  type(zlike_CS), pointer    :: CS !< Coordinate control structure
  real, optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_zlike_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_zlike_params

!> Builds a z* coordinate with a minimum thickness
subroutine build_zstar_column(CS, depth, total_thickness, zInterface, &
                              z_rigid_top, eta_orig, zScale)
  type(zlike_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive downward in the
                                                   !! output units), units may be [Z ~> m] or [H ~> m or kg m-2]
  real,                     intent(in)    :: total_thickness !< Column thickness (positive definite in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: z_rigid_top !< The height of a rigid top (positive upward in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: eta_orig !< The actual original height of the top (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: zScale !< Scaling factor from the target coordinate resolution
                                                    !! in Z to desired units for zInterface, perhaps Z_to_H,
                                                    !! often [nondim] or [H Z-1 ~> 1 or kg m-3]
  ! Local variables
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

end subroutine build_zstar_column

!> Builds a z* coordinate with a minimum thickness for each column of a block, as in build_zstar_column
!! without a rigid top or an original surface height, with the interfaces in the last array index.
subroutine build_zstar_block(CS, nib, njb, iie, jje, do_ij, depth, total_thickness, zInterface, zScale)
  type(zlike_CS), intent(in) :: CS  !< Coordinate control structure
  integer,        intent(in) :: nib !< The first dimension of the block arrays
  integer,        intent(in) :: njb !< The second dimension of the block arrays
  integer,        intent(in) :: iie !< The number of columns of this block in the first dimension
  integer,        intent(in) :: jje !< The number of columns of this block in the second dimension
  logical, dimension(nib,njb),         intent(in)    :: do_ij  !< True for the columns to work on
  real,    dimension(nib,njb),         intent(in)    :: depth  !< Depth of ocean bottom (positive downward
                                                   !! in the output units), units may be [Z ~> m] or [H ~> m or kg m-2]
  real,    dimension(nib,njb),         intent(in)    :: total_thickness !< Column thickness (positive definite in
                                                   !! the same units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real,    dimension(nib,njb,CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,                      intent(in)    :: zScale !< Scaling factor from the target coordinate
                                                    !! resolution in Z to desired units for zInterface, perhaps
                                                    !! Z_to_H, often [nondim] or [H Z-1 ~> 1 or kg m-3]
  ! Local variables
  real :: eta   ! Free surface height [Z ~> m] or [H ~> m or kg m-2]
  real :: stretching(nib,njb) ! A stretching factor for the coordinate [nondim]
  real :: min_thickness(nib,njb) ! The minimum layer thickness [Z ~> m] or [H ~> m or kg m-2]
  real :: dh, z0_top, z_scale ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]
  real :: CS_min_thickness ! A copy of CS%min_thickness [Z ~> m] or [H ~> m or kg m-2]
  real :: resolution(CS%nk) ! A copy of CS%coordinateResolution [Z ~> m]
  integer :: k, nk, ii, jj

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale
  nk = CS%nk ; CS_min_thickness = CS%min_thickness
  do k=1,nk ; resolution(k) = CS%coordinateResolution(k) ; enddo
  z0_top = 0.

  !$omp target enter data map(to: resolution) map(alloc: stretching, min_thickness)

  do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(eta))
    min_thickness(ii,jj) = min( CS_min_thickness, total_thickness(ii,jj)/real(nk) )

    ! Position of free-surface
    eta = total_thickness(ii,jj) - depth(ii,jj)

    ! Conventional z* coordinate:
    !   z* = (z-eta) / stretching   where stretching = (H+eta)/H
    !   z = eta + stretching * z*
    ! The above gives z*(z=eta) = 0, z*(z=-H) = -H.
    stretching(ii,jj) = total_thickness(ii,jj) / ( depth(ii,jj) + z0_top )

    ! Integrate down from the top for a notional new grid, ignoring topography
    zInterface(ii,jj,1) = eta
  enddo
  do concurrent (jj=1:jje)
    do k = 1,nk
      do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(dh))
        dh = stretching(ii,jj) * resolution(k)*z_scale ! Notional grid spacing
        zInterface(ii,jj,k+1) = zInterface(ii,jj,k) - dh
      enddo
    enddo
  enddo

  ! Integrating up from the bottom adjusting interface position to accommodate
  ! inflating layers without disturbing the interface above
  do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj))
    zInterface(ii,jj,nk+1) = -depth(ii,jj)
  enddo
  do concurrent (jj=1:jje)
    do k = nk,1,-1
      do concurrent (ii=1:iie, do_ij(ii,jj))
        if ( zInterface(ii,jj,k) < (zInterface(ii,jj,k+1) + min_thickness(ii,jj)) ) then
          zInterface(ii,jj,k) = zInterface(ii,jj,k+1) + min_thickness(ii,jj)
        endif
      enddo
    enddo
  enddo

  !$omp target exit data map(release: resolution, stretching, min_thickness)

end subroutine build_zstar_block

end module coord_zlike
