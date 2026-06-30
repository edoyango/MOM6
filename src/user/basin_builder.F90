! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> An idealized topography building system
module basin_builder

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_string_functions, only : lowercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public basin_builder_topography

! This include declares and sets the variable "version".
# include "version_variable.h"
character(len=40) :: mdl = "basin_builder" !< This module's name.


  interface
module subroutine basin_builder_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max [A]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units [A]
  ! Local variables

end subroutine basin_builder_topography
real module function cone(x, x0, L, clip)
  real,           intent(in) :: x    !< Coordinate in arbitrary units [A]
  real,           intent(in) :: x0   !< position of peak in arbitrary units [A]
  real,           intent(in) :: L    !< half-width of base of cone in arbitrary units [A]
  real, optional, intent(in) :: clip !< clipping height of cone [nondim]

end function cone
real module function scurve(x, x0, L)
  real, intent(in) :: x       !< Coordinate in arbitrary units [A]
  real, intent(in) :: x0      !< position of peak in arbitrary units [A]
  real, intent(in) :: L       !< half-width of base of cone in arbitrary units [A]

end function scurve
real module function cstprof(x, x0, L, lf, bf, sf, sh)
  real, intent(in) :: x       !< Coordinate in arbitrary units [A]
  real, intent(in) :: x0      !< position of peak in arbitrary units [A]
  real, intent(in) :: L       !< width of profile in arbitrary units [A]
  real, intent(in) :: lf      !< fraction of width that is "land" [nondim]
  real, intent(in) :: bf      !< fraction of width that is "beach" [nondim]
  real, intent(in) :: sf      !< fraction of width that is "continental slope" [nondim]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function cstprof
real module function dist_line_fixed_x(x, y, x0, y0, y1)
  real, intent(in) :: x       !< X-coordinate in arbitrary units [A]
  real, intent(in) :: y       !< Y-coordinate in arbitrary units [A]
  real, intent(in) :: x0      !< x-position of line segment in arbitrary units [A]
  real, intent(in) :: y0      !< y-position of line segment end in arbitrary units [A]
  real, intent(in) :: y1      !< y-position of line segment end in arbitrary units [A]

end function dist_line_fixed_x
real module function dist_line_fixed_y(x, y, x0, x1, y0)
  real, intent(in) :: x       !< X-coordinate in arbitrary units [A]
  real, intent(in) :: y       !< Y-coordinate in arbitrary units [A]
  real, intent(in) :: x0      !< x-position of line segment end in arbitrary units [A]
  real, intent(in) :: x1      !< x-position of line segment end in arbitrary units [A]
  real, intent(in) :: y0      !< y-position of line segment in arbitrary units [A]

end function dist_line_fixed_y
real module function angled_coast(lon, lat, lon_eq, lat_mer, dr, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon_eq  !< Longitude intersection with Equator [degrees_E]
  real, intent(in) :: lat_mer !< Latitude intersection with Prime Meridian [degrees_N]
  real, intent(in) :: dr      !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function angled_coast
real module function NS_coast(lon, lat, lonC, lat0, lat1, dlon, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of coast [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of coast end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function NS_coast
real module function EW_coast(lon, lat, latC, lon0, lon1, dlat, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: latC    !< Latitude of coast [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lon1    !< Longitude of coast end [degrees_E]
  real, intent(in) :: dlat    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function EW_coast
real module function NS_conic_ridge(lon, lat, lonC, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]

end function NS_conic_ridge
real module function NS_scurve_ridge(lon, lat, lonC, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]

end function NS_scurve_ridge
real module function circ_conic_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
  real, intent(in) :: lon            !< Longitude [degrees_E]
  real, intent(in) :: lat            !< Latitude [degrees_N]
  real, intent(in) :: lon0           !< Longitude of center of ring [degrees_E]
  real, intent(in) :: lat0           !< Latitude of center of ring [degrees_N]
  real, intent(in) :: ring_radius    !< Radius of ring [degrees]
  real, intent(in) :: ring_thickness !< Radial thickness of ring [degrees]
  real, intent(in) :: ridge_height   !< Ridge height as fraction of full depth [nondim]

end function circ_conic_ridge
real module function circ_scurve_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
  real, intent(in) :: lon            !< Longitude [degrees_E]
  real, intent(in) :: lat            !< Latitude [degrees_N]
  real, intent(in) :: lon0           !< Longitude of center of ring [degrees_E]
  real, intent(in) :: lat0           !< Latitude of center of ring [degrees_N]
  real, intent(in) :: ring_radius    !< Radius of ring [degrees]
  real, intent(in) :: ring_thickness !< Radial thickness of ring [degrees]
  real, intent(in) :: ridge_height   !< Ridge height as fraction of full depth [nondim]

end function circ_scurve_ridge
  end interface

end module basin_builder
