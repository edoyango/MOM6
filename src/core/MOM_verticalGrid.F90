! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides a transparent vertical ocean grid type and supporting routines
module MOM_verticalGrid

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public verticalGridInit, verticalGridEnd
public setVerticalGridAxes
public get_flux_units, get_thickness_units, get_tr_flux_units

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Describes the vertical ocean grid, including unit conversion factors
type, public :: verticalGrid_type

  ! Commonly used parameters
  integer :: ke     !< The number of layers/levels in the vertical
  real :: max_depth !< The maximum depth of the ocean [Z ~> m].
!  real :: mks_g_Earth !< The gravitational acceleration in unscaled MKS units [m s-2].  This might not be used.
  real :: g_Earth   !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: g_Earth_Z_T2 !< The gravitational acceleration in alternatively rescaled units [Z T-2 ~> m s-2]
  real :: Rho0      !< The density used in the Boussinesq approximation or nominal
                    !! density used to convert depths into mass units [R ~> kg m-3].

  ! Vertical coordinate descriptions for diagnostics and I/O
  character(len=40) :: zAxisUnits !< The units that vertical coordinates are written in
  character(len=40) :: zAxisLongName !< Coordinate name to appear in files,
                                  !! e.g. "Target Potential Density" or "Height"
  real, allocatable, dimension(:) :: sLayer !< Coordinate values of layer centers, in unscaled
                        !! units that depend on the vertical coordinate, such as [kg m-3] for an
                        !! isopycnal or some hybrid coordinates, [m] for a Z* coordinate,
                        !! or [nondim] for a sigma coordinate.
  real, allocatable, dimension(:) :: sInterface !< Coordinate values on interfaces, in the same
                        !! unscale units as sLayer [various].
  integer :: direction = 1 !< Direction defaults to 1, positive up.

  ! The following variables give information about the vertical grid.
  logical :: Boussinesq !< If true, make the Boussinesq approximation.
  logical :: semi_Boussinesq !< If true, do non-Boussinesq pressure force calculations and
                        !! use mass-based "thicknesses, but use Rho0 to convert layer thicknesses
                        !! into certain height changes.  This only applies if BOUSSINESQ is false.
  real :: Angstrom_H    !< A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2].
  real :: Angstrom_Z    !< A one-Angstrom thickness in the model depth units [Z ~> m].
  real :: Angstrom_m    !< A one-Angstrom thickness [m].
  real :: H_subroundoff !< A thickness that is so small that it can be added to a thickness of
                        !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
                        !! If Angstrom is 0 or exceedingly small, this is negligible compared to 1e-17 m.
  real :: dZ_subroundoff !< A thickness in height units that is so small that it can be added to a
                        !! vertical distance of Angstrom_Z or 1e-17 m without changing it at the bit
                        !! level [Z ~> m].  This is the height equivalent of H_subroundoff.
  real, allocatable, dimension(:) :: &
    g_prime, &          !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
    Rlay                !< The target coordinate value (potential density) in each layer [R ~> kg m-3].
  integer :: nkml = 0   !< The number of layers at the top that should be treated
                        !! as parts of a homogeneous region.
  integer :: nk_rho_varies = 0 !< The number of layers at the top where the
                        !! density does not track any target density.
  real :: H_to_kg_m2    !< A constant that translates thicknesses from the units of thickness
                        !! to kg m-2 [kg m-2 H-1 ~> kg m-3 or 1].
  real :: kg_m2_to_H    !< A constant that translates thicknesses from kg m-2 to the units
                        !! of thickness [H m2 kg-1 ~> m3 kg-1 or 1].
  real :: m_to_H        !< A constant that translates distances in m to the units of
                        !! thickness [H m-1 ~> 1 or kg m-3].
  real :: H_to_m        !< A constant that translates distances in the units of thickness
                        !! to m [m H-1 ~> 1 or m3 kg-1].
  real :: H_to_Pa       !< A constant that translates the units of thickness to pressure
                        !! [Pa H-1 ~> kg m-2 s-2 or m s-2].
  real :: H_to_Z        !< A constant that translates thickness units to the units of
                        !! depth [Z H-1 ~> 1 or m3 kg-1].
  real :: Z_to_H        !< A constant that translates depth units to thickness units
                        !! depth [H Z-1 ~> 1 or kg m-3].
  real :: H_to_RZ       !< A constant that translates thickness units to the units of
                        !! mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real :: RZ_to_H       !< A constant that translates mass per unit area units to
                        !! thickness units [H R-1 Z-1 ~> m3 kg-2 or 1].
  real :: H_to_MKS      !< A constant that translates thickness units to its MKS unit
                        !! (m or kg m-2) based on GV%Boussinesq [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  real :: m2_s_to_HZ_T  !< The combination of conversion factors that converts kinematic viscosities
                        !! in m2 s-1 to the internal units of the kinematic (in Boussinesq mode)
                        !! or dynamic viscosity [H Z s T-1 m-2 ~> 1 or kg m-3]
  real :: HZ_T_to_m2_s  !< The combination of conversion factors that converts the viscosities from
                        !! their internal representation into a kinematic viscosity in m2 s-1
                        !! [T m2 H-1 Z-1 s-1 ~> 1 or m3 kg-1]
  real :: HZ_T_to_MKS   !< The combination of conversion factors that converts the viscosities from
                        !! their internal representation into their unnscaled MKS units
                        !! (m2 s-1 or Pa s), depending on whether the model is Boussinesq
                        !! [T m2 H-1 Z-1 s-1 ~> 1] or [T Pa s H-1 Z-1 ~> 1]

end type verticalGrid_type


  interface
module subroutine verticalGridInit( param_file, GV, US )
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle/type
  type(verticalGrid_type), pointer    :: GV         !< The container for vertical grid data
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  ! This routine initializes the verticalGrid_type structure (GV).
  ! All memory is allocated but not necessarily set to meaningful values until later.

  ! Local variables
                     ! when in non-Boussinesq mode [R ~> kg m-3]
  ! This include declares and sets the variable "version".

end subroutine verticalGridInit
module function get_thickness_units(GV)
  character(len=48)                   :: get_thickness_units !< The vertical thickness units
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  !   This subroutine returns the appropriate units for thicknesses,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.

end function get_thickness_units
module function get_flux_units(GV)
  character(len=48)                   :: get_flux_units !< The thickness flux units
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  !   This subroutine returns the appropriate units for thickness fluxes,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.

end function get_flux_units
module function get_tr_flux_units(GV, tr_units, tr_vol_conc_units, tr_mass_conc_units)
  character(len=48)                      :: get_tr_flux_units !< The model's flux units
                                                              !! for a tracer.
  type(verticalGrid_type),    intent(in) :: GV                !< The ocean's vertical
                                                              !! grid structure.
  character(len=*), optional, intent(in) :: tr_units          !< Units for a tracer, for example
                                                              !! Celsius or PSU.
  character(len=*), optional, intent(in) :: tr_vol_conc_units !< The concentration units per unit
                                                              !! volume, for example if the units are
                                                              !! umol m-3, tr_vol_conc_units would
                                                              !! be umol.
  character(len=*), optional, intent(in) :: tr_mass_conc_units !< The concentration units per unit
                                                              !! mass of sea water, for example if
                                                              !! the units are mol kg-1,
                                                              !! tr_vol_conc_units would be mol.

  !   This subroutine returns the appropriate units for thicknesses and fluxes,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.

end function get_tr_flux_units
module subroutine setVerticalGridAxes( Rlay, GV, scale )
  type(verticalGrid_type), intent(inout) :: GV    !< The container for vertical grid data
  real, dimension(GV%ke),  intent(in)    :: Rlay  !< The layer target density [R ~> kg m-3]
  real,                    intent(in)    :: scale !< A unit scaling factor for Rlay to convert
                                                  !! it into the units of sInterface, usually
                                                  !! [kg m-3 R-1 ~> 1] when used in layer mode.
  ! Local variables

end subroutine setVerticalGridAxes
module subroutine verticalGridEnd( GV )
  type(verticalGrid_type), pointer :: GV !< The ocean's vertical grid structure

end subroutine verticalGridEnd
  end interface

end module MOM_verticalGrid
