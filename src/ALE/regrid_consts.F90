! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Contains constants for interpreting input parameters that control regridding.
module regrid_consts

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase

implicit none ; public

! List of regridding types. These should be consecutive and starting at 1.
! This allows them to be used as array indices.
integer, parameter :: REGRIDDING_LAYER     = 1      !< Layer mode identifier
integer, parameter :: REGRIDDING_ZSTAR     = 2      !< z* coordinates identifier
integer, parameter :: REGRIDDING_RHO       = 3      !< Density coordinates identifier
integer, parameter :: REGRIDDING_SIGMA     = 4      !< Sigma coordinates identifier
integer, parameter :: REGRIDDING_ARBITRARY = 5      !< Arbitrary coordinates identifier
integer, parameter :: REGRIDDING_HYCOM1    = 6      !< Simple HyCOM coordinates without BBL
integer, parameter :: REGRIDDING_SIGMA_SHELF_ZSTAR = 8 !< Identifiered for z* coordinates at the bottom,
                                                    !!  sigma-near the top
integer, parameter :: REGRIDDING_ADAPTIVE = 9       !< Adaptive coordinate mode identifier
integer, parameter :: REGRIDDING_HYBGEN    = 10     !< Hybgen coordinates identifier

character(len=*), parameter :: REGRIDDING_LAYER_STRING = "LAYER"   !< Layer string
character(len=*), parameter :: REGRIDDING_ZSTAR_STRING_OLD = "Z*"  !< z* string (legacy name)
character(len=*), parameter :: REGRIDDING_ZSTAR_STRING = "ZSTAR"   !< z* string
character(len=*), parameter :: REGRIDDING_RHO_STRING   = "RHO"     !< Rho string
character(len=*), parameter :: REGRIDDING_SIGMA_STRING = "SIGMA"   !< Sigma string
character(len=*), parameter :: REGRIDDING_ARBITRARY_STRING = "ARB" !< Arbitrary coordinates
character(len=*), parameter :: REGRIDDING_HYCOM1_STRING = "HYCOM1" !< Hycom string
character(len=*), parameter :: REGRIDDING_HYBGEN_STRING = "HYBGEN" !< Hybgen string
character(len=*), parameter :: REGRIDDING_SIGMA_SHELF_ZSTAR_STRING = "SIGMA_SHELF_ZSTAR" !< Hybrid z*/sigma
character(len=*), parameter :: REGRIDDING_ADAPTIVE_STRING = "ADAPTIVE" !< Adaptive coordinate string
character(len=*), parameter :: DEFAULT_COORDINATE_MODE = REGRIDDING_LAYER_STRING !< Default coordinate mode

!> Returns a string with the coordinate units associated with the coordinate mode.
interface coordinateUnits
  module procedure coordinateUnitsI
  module procedure coordinateUnitsS
end interface

!> Returns true if the coordinate is dependent on the state density, returns false otherwise.
interface state_dependent
  module procedure state_dependent_char
  module procedure state_dependent_int
end interface


  interface
module function coordinateMode(string)
  integer :: coordinateMode !< Enumerated integer indicating coordinate mode
  character(len=*), intent(in) :: string !< String to indicate coordinate mode
end function coordinateMode
module function coordinateUnitsI(coordMode)
  character(len=16) :: coordinateUnitsI !< Units of coordinate
  integer, intent(in) :: coordMode !< Coordinate mode
end function coordinateUnitsI
module function coordinateUnitsS(string)
  character(len=16) :: coordinateUnitsS !< Units of coordinate
  character(len=*), intent(in) :: string !< Coordinate mode
end function coordinateUnitsS
logical module function state_dependent_char(string)
  character(len=*), intent(in) :: string !< String to indicate coordinate mode

end function state_dependent_char
logical module function state_dependent_int(mode)
  integer, intent(in) :: mode !< Coordinate mode
end function state_dependent_int
  end interface

end module regrid_consts
