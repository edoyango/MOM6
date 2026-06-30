! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for the continuous isopycnal (rho) coordinate
module coord_rho

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the rho coordinate
type, public :: rho_CS ; private

  !> Number of layers
  integer :: nk

  !> Minimum thickness allowed for layers, often in [H ~> m or kg m-2]
  real :: min_thickness = 0.

  !> Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .false.

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type rho_CS

public init_coord_rho, set_rho_params, build_rho_column, old_inflate_layers_1d, end_coord_rho


  interface
module subroutine init_coord_rho(CS, nk, ref_pressure, target_density, interp_CS)
  type(rho_CS),         pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in the grid
  real,                 intent(in) :: ref_pressure !< Coordinate reference pressure [R L2 T-2 ~> Pa]
  real, dimension(:),   intent(in) :: target_density !< Nominal density of interfaces [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

end subroutine init_coord_rho
module subroutine end_coord_rho(CS)
  type(rho_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
end subroutine end_coord_rho
module subroutine set_rho_params(CS, min_thickness, integrate_downward_for_e, interp_CS, ref_pressure)
  type(rho_CS),      pointer    :: CS !< Coordinate control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface
                                      !! positions from the top downward.  If false, integrate
                                      !! from the bottom upward, as does the rest of the model.
  real,    optional, intent(in) :: ref_pressure     !< The reference pressure for density-dependent
                                                    !! coordinates [R L2 T-2 ~> Pa]

  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

end subroutine set_rho_params
module subroutine build_rho_column(CS, nz, depth, h, T, S, eqn_of_state, z_interface, &
                            z_rigid_top, eta_orig, h_neglect, h_neglect_edge)
  type(rho_CS),        intent(in)    :: CS !< coord_rho control structure
  integer,             intent(in)    :: nz !< Number of levels on source grid (i.e. length of  h, T, S)
  real,                intent(in)    :: depth !< Depth of ocean bottom (positive downward) [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: T  !< Temperature for source column [C ~> degC]
  real, dimension(nz), intent(in)    :: S  !< Salinity for source column [S ~> ppt]
  type(EOS_type),      intent(in)    :: eqn_of_state !< Equation of state structure
  real, dimension(CS%nk+1), &
                       intent(inout) :: z_interface !< Absolute positions of interfaces [H ~> m or kg m-2]
  real, optional,      intent(in)    :: z_rigid_top !< The height of a rigid top (positive upward in the same
                                             !! units as depth) [H ~> m or kg m-2]
  real, optional,      intent(in)    :: eta_orig !< The actual original height of the top in the same
                                                   !! units as depth) [H ~> m or kg m-2]
  real,                intent(in)    :: h_neglect !< A negligibly small width for the purpose
                                             !! of cell reconstructions [H ~> m or kg m-2]
  real,      optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                             !! of edge value calculations [H ~> m or kg m-2]

  ! Local variables

  ! Construct source column with vanished layers removed (stored in h_nv)
end subroutine build_rho_column
module subroutine build_rho_column_iteratively(CS, remapCS, nz, depth, h, T, S, eqn_of_state, &
                                        zInterface, h_neglect, h_neglect_edge, dev_tol)
  type(rho_CS),          intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom [Z ~> m]
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses in Z coordinates [Z ~> m]
  real, dimension(nz),   intent(in)    :: T  !< T for column [C ~> degC]
  real, dimension(nz),   intent(in)    :: S  !< S for column [S ~> ppt]
  type(EOS_type),        intent(in)    :: eqn_of_state !< Equation of state structure
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces [Z ~> m]
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions
                                             !! in the same units as h [Z ~> m]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width
                                             !! for the purpose of edge value calculations
                                             !! in the same units as h [Z ~> m]
  real,        optional, intent(in)    :: dev_tol !< The tolerance for the deviation between
                                             !! successive grids for determining when the
                                             !! iterative solver has converged [Z ~> m]

  ! Local variables
                               ! deviation between two successive grids [Z ~> m].
                               ! regridding iterations [Z ~> m]

  !  Maximum number of regridding iterations

end subroutine build_rho_column_iteratively
module subroutine copy_finite_thicknesses(nk, h_in, thresh, nout, h_out, mapping)
  integer,                intent(in)  :: nk      !< Number of layer for h_in, T_in, S_in
  real, dimension(nk),    intent(in)  :: h_in    !< Thickness of input column [H ~> m or kg m-2] or [Z ~> m]
  real,                   intent(in)  :: thresh  !< Thickness threshold defining vanished
                                                 !! layers [H ~> m or kg m-2] or [Z ~> m]
  integer,                intent(out) :: nout    !< Number of non-vanished layers
  real, dimension(nk),    intent(out) :: h_out   !< Thickness of output column [H ~> m or kg m-2] or [Z ~> m]
  integer, dimension(nk), intent(out) :: mapping !< Index of k-out corresponding to k-in
  ! Local variables

  ! Build up new grid
end subroutine copy_finite_thicknesses
module subroutine old_inflate_layers_1d( min_thickness, nk, h )

  ! Argument
  real,               intent(in)    :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2] or other units
  integer,            intent(in)    :: nk  !< Number of layers in the grid
  real, dimension(:), intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2] or other units

  ! Local variable
                             ! same units as h, often [H ~> m or kg m-2]
                             ! to give mass conservation in the same units as h, often [H ~> m or kg m-2]

  ! Count number of nonzero layers
end subroutine old_inflate_layers_1d
  end interface

end module coord_rho
