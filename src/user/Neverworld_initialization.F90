! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization for the "Neverworld" configuration
module Neverworld_initialization

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream

implicit none ; private

#include <MOM_memory.h>

public Neverworld_initialize_topography
public Neverworld_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine Neverworld_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max [A]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units [A]

  ! Local variables
  ! This include declares and sets the variable "version".
end subroutine Neverworld_initialize_topography
real module function cosbell(x, L)
  real , intent(in) :: x       !< Position in arbitrary units [A]
  real , intent(in) :: L       !< Width in arbitrary units [A]

end function cosbell
real module function spike(x, L)

  real , intent(in) :: x       !< Position in arbitrary units [A]
  real , intent(in) :: L       !< Width in arbitrary units [A]

end function spike
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
real module function NS_coast(lon, lat, lon0, lat0, lat1, dlon, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of coast end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function NS_coast
real module function EW_coast(lon, lat, lon0, lon1, lat0, dlat, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lon1    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast [degrees_N]
  real, intent(in) :: dlat    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]

end function EW_coast
real module function NS_ridge(lon, lat, lon0, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]

end function NS_ridge
real module function circ_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
  real, intent(in) :: lon            !< Longitude [degrees_E]
  real, intent(in) :: lat            !< Latitude [degrees_N]
  real, intent(in) :: lon0           !< Longitude of center of ring [degrees_E]
  real, intent(in) :: lat0           !< Latitude of center of ring [degrees_N]
  real, intent(in) :: ring_radius    !< Radius of ring [degrees]
  real, intent(in) :: ring_thickness !< Radial thickness of ring [degrees]
  real, intent(in) :: ridge_height   !< Ridge height as fraction of full depth [nondim]

end function circ_ridge
module subroutine Neverworld_initialize_thickness(h, depth_tot, G, GV, US, param_file, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US                   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h !< The thickness that is being
                                                              !! initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure [R L2 T-2 ~> Pa].
  ! Local variables
                            ! usually negative because it is positive upward.
                  ! by the domain sizes [nondim]
                  ! by the domain sizes [nondim]

end subroutine Neverworld_initialize_thickness
  end interface

end module Neverworld_initialization
