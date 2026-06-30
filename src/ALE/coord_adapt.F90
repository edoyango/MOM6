! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for the adaptive coordinate
module coord_adapt

use MOM_EOS,           only : calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure for adaptive coordinates (coord_adapt).
type, public :: adapt_CS ; private

  !> Number of layers/levels
  integer :: nk

  !> Nominal near-surface resolution [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: coordinateResolution

  !> Ratio of optimisation and diffusion timescales [nondim]
  real :: adaptTimeRatio

  !> Nondimensional coefficient determining how much optimisation to apply [nondim]
  real :: adaptAlpha

  !> Near-surface zooming depth [H ~> m or kg m-2]
  real :: adaptZoom

  !> Near-surface zooming coefficient [nondim]
  real :: adaptZoomCoeff

  !> Stratification-dependent diffusion coefficient [nondim]
  real :: adaptBuoyCoeff

  !> Reference density difference for stratification-dependent diffusion [R ~> kg m-3]
  real :: adaptDrho0

  !> If true, form a HYCOM1-like mixed layet by preventing interfaces
  !! from becoming shallower than the depths set by coordinateResolution
  logical :: adaptDoMin  = .false.
end type adapt_CS

public init_coord_adapt, set_adapt_params, build_adapt_column, end_coord_adapt


  interface
module subroutine init_coord_adapt(CS, nk, coordinateResolution, m_to_H, kg_m3_to_R)
  type(adapt_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of layers in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Nominal near-surface resolution [m] or
                                       !! other units specified with m_to_H
  real,               intent(in) :: m_to_H !< A conversion factor from m to the units of thicknesses,
                                       !! perhaps in units of [H m-1 ~> 1 or kg m-3]
  real,               intent(in) :: kg_m3_to_R !< A conversion factor from kg m-3 to the units of density,
                                       !! perhaps in units of [R m3 kg-1 ~> 1]

end subroutine init_coord_adapt
module subroutine end_coord_adapt(CS)
  type(adapt_CS), pointer :: CS  !< The control structure for this module

  ! nothing to do
end subroutine end_coord_adapt
module subroutine set_adapt_params(CS, adaptTimeRatio, adaptAlpha, adaptZoom, adaptZoomCoeff, &
                            adaptBuoyCoeff, adaptDrho0, adaptDoMin)
  type(adapt_CS),    pointer    :: CS  !< The control structure for this module
  real,    optional, intent(in) :: adaptTimeRatio !< Ratio of optimisation and diffusion timescales [nondim]
  real,    optional, intent(in) :: adaptAlpha     !< Nondimensional coefficient determining
                                                  !! how much optimisation to apply [nondim]
  real,    optional, intent(in) :: adaptZoom      !< Near-surface zooming depth [H ~> m or kg m-2]
  real,    optional, intent(in) :: adaptZoomCoeff !< Near-surface zooming coefficient [nondim]
  real,    optional, intent(in) :: adaptBuoyCoeff !< Stratification-dependent diffusion coefficient [nondim]
  real,    optional, intent(in) :: adaptDrho0  !< Reference density difference for
                                               !! stratification-dependent diffusion [R ~> kg m-3]
  logical, optional, intent(in) :: adaptDoMin  !< If true, form a HYCOM1-like mixed layer by
                                               !! preventing interfaces from becoming shallower than
                                               !! the depths set by coordinateResolution

end subroutine set_adapt_params
module subroutine build_adapt_column(CS, G, GV, US, tv, i, j, zInt, tInt, sInt, h, nom_depth_H, zNext)
  type(adapt_CS),                              intent(in)    :: CS   !< The control structure for this module
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  integer,                                     intent(in)    :: i    !< The i-index of the column to work on
  integer,                                     intent(in)    :: j    !< The j-index of the column to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: zInt !< Interface heights [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: tInt !< Interface temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: sInt !< Interface salinities [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                     !! relative to mean sea level or another locally
                                                                     !! valid reference height, converted to thickness
                                                                     !! units [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1),                 intent(inout) :: zNext !< updated interface positions [H ~> m or kg m-2]

  ! Local variables
                      ! adjustive fluxes [H ~> m or kg m-2]

end subroutine build_adapt_column
  end interface

end module coord_adapt
