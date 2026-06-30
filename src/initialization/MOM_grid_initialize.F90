! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initializes horizontal grid
module MOM_grid_initialize

use MOM_checksums,     only : hchksum, Bchksum, uvchksum, hchksum_pair, Bchksum_pair
use MOM_domains,       only : pass_var, pass_vector, pe_here, root_PE, broadcast
use MOM_domains,       only : AGRID, BGRID_NE, CGRID_NE, To_All, Scalar_Pair
use MOM_domains,       only : To_North, To_South, To_East, To_West
use MOM_domains,       only : MOM_domain_type, clone_MOM_domain, deallocate_MOM_domain
use MOM_dyn_horgrid,   only : dyn_horgrid_type, set_derived_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher, file_exists, stdout
use MOM_io,            only : CORNER, NORTH_FACE, EAST_FACE
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public set_grid_metrics, initialize_masks, Adcroft_reciprocal

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Global positioning system (aka container for information to describe the grid)
type, private :: GPS ; private
  real :: len_lon  !< The longitudinal or x-direction length of the domain [degrees_E] or [km] or [m].
  real :: len_lat  !< The latitudinal or y-direction length of the domain [degrees_N] or [km] or [m].
  real :: west_lon !< The western longitude of the domain or the equivalent
                   !! starting value for the x-axis [degrees_E] or [km] or [m].
  real :: south_lat  !< The southern latitude of the domain or the equivalent
                   !! starting value for the y-axis [degrees_N] or [km] or [m].
  real :: Rad_Earth_L !< The radius of the Earth in rescaled units [L ~> m]
  real :: Lat_enhance_factor  !< The amount by which the meridional resolution
                   !! is enhanced within LAT_EQ_ENHANCE of the equator [nondim]
  real :: Lat_eq_enhance !< The latitude range to the north and south of the equator
                   !! over which the resolution is enhanced [degrees_N]
  logical :: isotropic !< If true, an isotropic grid on a sphere (also known as a Mercator grid)
                   !! is used. With an isotropic grid, the meridional extent of the domain
                   !! (LENLAT), the zonal extent (LENLON), and the number of grid points in each
                   !! direction are _not_ independent. In MOM the meridional extent is determined
                   !! to fit the zonal extent and the number of grid points, while grid is
                   !! perfectly isotropic.
  logical :: equator_reference !< If true, the grid is defined to have the equator at the
                   !!  nearest q- or h- grid point to (-LOWLAT*NJGLOBAL/LENLAT).
  integer :: niglobal !< The number of i-points in the global grid computational domain
  integer :: njglobal !< The number of j-points in the global grid computational domain
end type GPS


  interface
module subroutine set_grid_metrics(G, param_file, US)
  type(dyn_horgrid_type),          intent(inout) :: G  !< The dynamic horizontal grid type
  type(param_file_type),           intent(in)    :: param_file !< Parameter file structure
  type(unit_scale_type),           intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine set_grid_metrics
module subroutine grid_metrics_chksum(parent, G, US)
  character(len=*),       intent(in) :: parent !< A string identifying the caller
  type(dyn_horgrid_type), intent(in) :: G      !< The dynamic horizontal grid type
  type(unit_scale_type),  intent(in) :: US !< A dimensional unit scaling type


end subroutine grid_metrics_chksum
module subroutine set_grid_metrics_from_mosaic(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G           !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: param_file  !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US          !< A dimensional unit scaling type

  ! Local variables
  ! These are symmetric arrays, corresponding to the data in the mosaic file
                                                                   ! longitudes [degrees_E]

end subroutine set_grid_metrics_from_mosaic
module subroutine set_grid_metrics_cartesian(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G           !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: param_file  !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US    !< A dimensional unit scaling type
  ! Local variables

end subroutine set_grid_metrics_cartesian
module subroutine set_grid_metrics_spherical(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G           !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: param_file  !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US    !< A dimensional unit scaling type
  ! Local variables

end subroutine set_grid_metrics_spherical
module subroutine set_grid_metrics_mercator(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G           !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: param_file  !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US    !< A dimensional unit scaling type
  ! Local variables
                          ! Int_dj_dy at a latitude or longitude that is
                          ! being set to be at grid index jRef or iRef [gridpoints]

  !   All of the metric terms should be defined over the domain from
  ! isd to ied.  Outside of the physical domain, both the metrics
  ! and their inverses may be set to zero.
end subroutine set_grid_metrics_mercator
module function ds_di(x, y, GP)
  real, intent(in) :: x  !< The longitude in question [radians]
  real, intent(in) :: y  !< The latitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters

  real :: ds_di  ! The returned grid spacing [L ~> m]

end function ds_di
module function ds_dj(x, y, GP)
  real, intent(in) :: x  !< The longitude in question [radians]
  real, intent(in) :: y  !< The latitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters

  real :: ds_dj  ! The returned grid spacing [L ~> m]

end function ds_dj
module function  dL(x1, x2, y1, y2)
  real, intent(in) :: x1 !< Segment starting longitude [radians]
  real, intent(in) :: x2 !< Segment ending longitude [radians]
  real, intent(in) :: y1 !< Segment starting latitude [radians]
  real, intent(in) :: y2 !< Segment ending latitude [radians]
  ! Local variables
  real :: dL ! A contribution to the spanned area the surface of the sphere [radian2]

end function dL
module function find_root( fn, dy_df, GP, fnval, y1, ymin, ymax, ittmax)
  real :: find_root !< The value of y where fn(y) = fnval that will be returned [radians]
  real,      external    :: fn    !< The external function whose root is being sought [gridpoints]
  real,      external    :: dy_df !< The inverse of the derivative of that function [radian gridpoint-1]
  type(GPS), intent(in)  :: GP    !< A structure of grid parameters
  real,      intent(in)  :: fnval !< The value of fn being sought [gridpoints]
  real,      intent(in)  :: y1    !< A first guess for y [radians]
  real,      intent(in)  :: ymin  !< The minimum permitted value of y [radians]
  real,      intent(in)  :: ymax  !< The maximum permitted value of y [radians]
  integer,   intent(out) :: ittmax !< The number of iterations used to polish the root
  ! Local variables

!  Bracket the root.  Do not use the bounding values because the value at the
! function at the bounds could be infinite, as is the case for the Mercator
! grid recursion relation. (I.e., this is a search on an open interval.)
end function find_root
module function dx_di(x, GP)
  real, intent(in) :: x !< The longitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters
  real :: dx_di         ! The derivative of zonal position with the grid index [radian gridpoint-1]

end function dx_di
module function Int_di_dx(x, GP)
  real, intent(in) :: x  !< The longitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters
  real :: Int_di_dx   ! A position in the global i-index space [gridpoints]

end function Int_di_dx
module function dy_dj(y, GP)
  real, intent(in) :: y !< The latitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters
  real :: dy_dj         ! The derivative of meridional position with the grid index [radian gridpoint-1]
  ! Local variables
                        ! gridpoints to the nominal spacing in Radians [radian gridpoint-1]
                        ! is enhanced [radians]
end function dy_dj
module function Int_dj_dy(y, GP)
  real, intent(in) :: y  !< The latitude in question [radians]
  type(GPS), intent(in) :: GP  !< A structure of grid parameters
  real :: Int_dj_dy        ! The grid position of latitude y [gridpoints]
  ! Local variables
                           ! nominal spacing in gridpoints to the nominal
                           ! spacing in Radians [gridpoint radian-1]
                           ! grid spacing is enhanced by a factor of GP%lat_enhance_factor [radians]

end function Int_dj_dy
module subroutine extrapolate_metric(var, jh, missing)
  real, dimension(:,:), intent(inout) :: var     !< The array in which to fill in halos in arbitrary units [A]
  integer,              intent(in)    :: jh      !< The size of the halos to be filled
  real,       optional, intent(in)    :: missing !< The missing data fill value, 0 by default [A]
  ! Local variables

end subroutine extrapolate_metric
module function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted in arbitrary units [A]
  real :: I_val            !< The Adcroft reciprocal of val [A-1]

end function Adcroft_reciprocal
module subroutine initialize_masks(G, PF, US, OBC_dir_u, OBC_dir_v, open_corner_OBCs, maskT)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid type
  type(param_file_type),  intent(in)    :: PF !< Parameter file structure
  type(unit_scale_type),  intent(in)    :: US !< A dimensional unit scaling type
  integer, dimension(G%IsdB:G%IedB,G%jsd:G%jed), &
                optional, intent(in)    :: OBC_dir_u  !< Trinary values that indicate whether there
                                              !! is an open boundary condition at zonal velocity
                                              !! faces and their orientation, with 0 for no OBC,
                                              !! a positive value for an Eastern OBC and
                                              !! a negative value for a Western OBC.
  integer, dimension(G%isd:G%ied,G%JsdB:G%JedB), &
                optional, intent(in)    :: OBC_dir_v  !< Trinary values that indicate whether there
                                              !! is an open boundary condition at zonal velocity
                                              !! faces and their orientation, with 0 for no OBC,
                                              !! a positive value for a Northern OBC and
                                              !! a negative value for a Southern OBC.
  logical,      optional, intent(in)   :: open_corner_OBCs  !< If present and true, the bay-like corner
                                              !! between two orthogonal open boundary segments is open,
                                              !! otherwise it is closed.
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                optional, intent(in)   :: maskT !< If present, this array is used to set the
                                              !! the mask at tracer points instead of using the
                                              !! bathymetry to determine the masks [nondim]

  ! Local variables

end subroutine initialize_masks
  end interface

end module MOM_grid_initialize
