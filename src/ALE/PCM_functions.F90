! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise constant reconstruction functions
module PCM_functions

implicit none ; private

public PCM_reconstruction


  interface
module subroutine PCM_reconstruction( N, u, edge_values, ppoly_coef )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: u !< cell averages in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial,
                                           !! with the same units as u [A].
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial,
                                           !! with the same units as u [A].

  ! Local variables

  ! The coefficients of the piecewise constant polynomial are simply
  ! the cell averages.
end subroutine PCM_reconstruction
  end interface

end module PCM_functions
