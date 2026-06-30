! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for the HyCOM coordinate
module coord_hycom

use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, NOTE
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_remapping,     only : remapping_CS, remapping_core_h
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, regridding_set_ppolys
use regrid_interp,     only : DEGREE_MAX

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom_CS ; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Maximum depths of interfaces [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_layer_thickness

  !> If true, an interface only moves if it improves the density fit
  logical :: only_improves = .false.

  !> If true, use 3-D control fields
  logical :: use_3d = .false.

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: target_density_3d

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:,:,:) :: coordinateResolution_3d

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom_CS

public init_coord_hycom, init_3d_coord_hycom, set_hycom_params, build_hycom1_column, end_coord_hycom


  interface
module subroutine init_coord_hycom(CS, nk, coordinateResolution, target_density, interp_CS)
  type(hycom_CS),       pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  real, dimension(nk),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(nk+1),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

end subroutine init_coord_hycom
module subroutine init_3d_coord_hycom(CS, G, nk, coordinateResolution, target_density, interp_CS)
  type(hycom_CS),       pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(ocean_grid_type),intent(in) :: G   !< Ocean grid structure
  integer,              intent(in) :: nk  !< Number of layers in generated grid
  real, dimension(SZI_(G),SZJ_(G),nk), intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),nk+1), intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation
  ! Local variables

end subroutine init_3d_coord_hycom
module subroutine end_coord_hycom(CS)
  type(hycom_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
end subroutine end_coord_hycom
module subroutine set_hycom_params(CS, max_interface_depths, max_layer_thickness, only_improves, interp_CS)
  type(hycom_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces [H ~> m or kg m-2]
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  logical, optional, intent(in) :: only_improves !< If true, an interface only moves if it improves the density fit
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

end subroutine set_hycom_params
module subroutine build_hycom1_column(CS, remapCS, eqn_of_state, nz, ix, jy, depth, h, T, S, p_col, &
                               z_col, z_col_new, zScale, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)    :: CS    !< Coordinate control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)    :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz    !< Number of levels
  integer,               intent(in)    :: ix    !< x direction array index
  integer,               intent(in)    :: jy    !< y direction array index
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T     !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)    :: S     !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: z_col_new !< Absolute positions of interfaces [H ~> m or kg m-2]
  real, optional,        intent(in)    :: zScale !< Scaling factor from the input coordinate thicknesses in [Z ~> m]
                                                !! to desired units for zInterface, perhaps GV%Z_to_H in which
                                                !! case this has units of [H Z-1 ~> nondim or kg m-3]
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                !! cell reconstruction [H ~> m or kg m-2]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]

  ! Local variables
                                        ! interface target densities [R ~> kg m-3]
                                        ! interface target densities [R ~> kg m-3]
                     ! perhaps 1 or a factor in [H Z-1 ~> 1 or kg m-3]

end subroutine build_hycom1_column
module subroutine build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, nz, ix, jy, depth, h, T, S, &
                                       p_col, R, RiAnom, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)  :: CS     !< Coordinate control structure
  type(remapping_CS),    intent(in)  :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)  :: eqn_of_state !< Equation of state structure
  integer,               intent(in)  :: nz     !< Number of levels
  integer,               intent(in)  :: ix     !< x direction array index
  integer,               intent(in)  :: jy     !< y direction array index
  real,                  intent(in)  :: depth  !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)  :: T      !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)  :: S      !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)  :: p_col  !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz),   intent(out) :: R      !< Layer density [R ~> kg m-3]
  real, dimension(nz+1), intent(out) :: RiAnom !< The interface density anomaly
                                               !! w.r.t. the interface target
                                               !! densities [R ~> kg m-3]
  real,                  intent(in)  :: h_neglect !< A negligibly small width for the purpose of
                                               !! cell reconstruction [H ~> m or kg m-2]
  real,        optional, intent(in)  :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]
  ! Local variables

  ! Work bottom recording potential density
end subroutine build_hycom1_target_anomaly
  end interface

end module coord_hycom
