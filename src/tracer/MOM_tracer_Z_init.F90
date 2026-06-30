! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Used to initialize tracers from a depth- (or z*-) space file.
module MOM_tracer_Z_init

use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : MOM_read_data, get_var_sizes, read_attribute, read_variable
use MOM_io, only : open_file_to_read, close_file_to_read
use MOM_EOS, only : EOS_type, calculate_density, calculate_density_derivs, EOS_domain
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public tracer_Z_init, read_Z_edges, tracer_Z_init_array, determine_temperature

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module function tracer_Z_init(tr, h, filename, tr_name, G, GV, US, missing_val, land_val, scale)
  logical :: tracer_Z_init !< A return code indicating if the initialization has been successful
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                         intent(out)   :: tr   !< The tracer to initialize [CU ~> conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                         intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2] or other
                                               !! arbitrary units such as [Z ~> m]
  character(len=*),      intent(in)    :: filename  !< The name of the file to read from
  character(len=*),      intent(in)    :: tr_name   !< The name of the tracer in the file
  real,        optional, intent(in)    :: missing_val !< The missing value for the tracer [CU ~> conc]
  real,        optional, intent(in)    :: land_val  !< A value to use to fill in land points [CU ~> conc]
  real,        optional, intent(in)    :: scale     !< A factor by which to scale the output tracers from the
                                                    !! their units in the file [CU conc-1 ~> 1]

  ! Local variables
end function tracer_Z_init
module subroutine tracer_z_init_array(tr_in, z_edges, nk_data, e, land_fill, G, nlay, nlevs, &
                               eps_z, tr, scale)
  type(ocean_grid_type),      intent(in)  :: G     !< The ocean's grid structure
  integer,                    intent(in)  :: nk_data !< The number of levels in the input data
  real, dimension(SZI_(G),SZJ_(G),nk_data), &
                              intent(in)  :: tr_in !< The z-space array of tracer concentrations
                                                   !! that is read in [A]
  real, dimension(nk_data+1), intent(in)  :: z_edges !< The depths of the cell edges in the input z* data
                                                   !! [Z ~> m] or [m]
  integer,                    intent(in)  :: nlay  !< The number of vertical layers in the target grid
  real, dimension(SZI_(G),SZJ_(G),nlay+1), &
                              intent(in)  :: e     !< The depths of the target layer interfaces [Z ~> m] or [m]
  real,                       intent(in)  :: land_fill !< fill in data over land [B]
  integer, dimension(SZI_(G),SZJ_(G)), &
                              intent(in)  :: nlevs !< The number of input levels with valid data
  real,                       intent(in)  :: eps_z !< A negligibly thin layer thickness [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),nlay), &
                              intent(out) :: tr    !< tracers in model space [B]
  real,             optional, intent(in)  :: scale !< A factor by which to scale the output tracers from the
                                                   !! input tracers [B A-1 ~> 1]

  ! Local variables
                         ! a layer, relative to the cell center and normalized by the cell thickness [nondim].
                         ! a layer, relative to the cell center and normalized by the cell thickness [nondim].
                         ! Note that -1/2 <= z1 <= z2 <= 1/2.

end subroutine tracer_z_init_array
module subroutine read_Z_edges(filename, tr_name, z_edges, nz_out, has_edges, &
                        use_missing, missing, scale, missing_scale)
  character(len=*), intent(in)    :: filename !< The name of the file to read from.
  character(len=*), intent(in)    :: tr_name !< The name of the tracer in the file.
  real, dimension(:), allocatable, &
                    intent(out)   :: z_edges !< The depths of the vertical edges of the tracer array [Z ~> m]
  integer,          intent(out)   :: nz_out  !< The number of vertical layers in the tracer array
  logical,          intent(out)   :: has_edges !< If true the values in z_edges are the edges of the
                                             !! tracer cells, otherwise they are the cell centers
  logical,          intent(inout) :: use_missing !< If false on input, see whether the tracer has a
                                             !! missing value, and if so return true
  real,             intent(inout) :: missing !< The missing value, if one has been found [CU ~> conc]
  real,             intent(in)    :: scale   !< A scaling factor for z_edges into new units [Z m-1 ~> 1]
  real,             intent(in)    :: missing_scale  !< A scaling factor to use to convert the
                                             !! tracers and their missing value from the units in
                                             !! the file into their internal units [CU conc-1 ~> 1]

  !   This subroutine reads the vertical coordinate data for a field from a
  ! NetCDF file.  It also might read the missing value attribute for that same field.

end subroutine read_Z_edges
module subroutine find_overlap(e, Z_top, Z_bot, k_max, k_start, k_top, k_bot, wt, z1, z2)
  real, dimension(:), intent(in)  :: e      !< Column interface heights, [Z ~> m] or other units.
  real,               intent(in)  :: Z_top  !< Top of range being mapped to, in the units of e [Z ~> m].
  real,               intent(in)  :: Z_bot  !< Bottom of range being mapped to, in the units of e [Z ~> m].
  integer,            intent(in)  :: k_max  !< Number of valid layers.
  integer,            intent(in)  :: k_start !< Layer at which to start searching.
  integer,            intent(out) :: k_top  !< Indices of top layers that overlap with the depth range.
  integer,            intent(out) :: k_bot  !< Indices of bottom layers that overlap with the depth range.
  real, dimension(:), intent(out) :: wt     !< Relative weights of each layer from k_top to k_bot [nondim].
  real, dimension(:), intent(out) :: z1     !< Depth of the top limits of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness [nondim].  Note that -1/2 <= z1 < z2 <= 1/2.
  real, dimension(:), intent(out)   :: z2     !< Depths of the bottom limit of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness [nondim].  Note that -1/2 <= z1 < z2 <= 1/2.

  ! Local variables

end subroutine find_overlap
module function find_limited_slope(val, e, k) result(slope)
  real, dimension(:), intent(in) :: val !< A column of the values that are being interpolated, in arbitrary units [A]
  real, dimension(:), intent(in) :: e   !< A column's interface heights [Z ~> m] or other units.
  integer,            intent(in) :: k   !< The layer whose slope is being determined.
  real :: slope !< The normalized slope in the intracell distribution of val [A]
  ! Local variables

end function find_limited_slope
module subroutine determine_temperature(temp, salt, R_tgt, EOS, p_ref, niter, k_start, G, GV, US, PF, &
                                 just_read)
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),       intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: temp !< potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: salt !< salinity [S ~> ppt]
  real, dimension(SZK_(GV)),     intent(in)    :: R_tgt !< desired potential density [R ~> kg m-3].
  type(EOS_type),                intent(in)    :: EOS   !< seawater equation of state control structure
  real,                          intent(in)    :: p_ref !< reference pressure [R L2 T-2 ~> Pa].
  integer,                       intent(in)    :: niter !< maximum number of iterations
  integer,                       intent(in)    :: k_start !< starting index (i.e. below the buffer layer)
  type(unit_scale_type),         intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: PF  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                       intent(in)    :: just_read !< If true, this call will only read
                                                      !! parameters without changing T or S.

  ! Local variables (All of which need documentation!)
                        ! minimizing property changes while correcting density [C S-1 ~> degC ppt-1].
                        ! T-S space when stretched with dT_dS_gauge [S2 R-2 ~> ppt2 m6 kg-2]
                    ! when old_fit is true [C ~> degC]
                    ! when old_fit is true [S ~> ppt]
  ! This include declares and sets the variable "version".

end subroutine determine_temperature
  end interface

end module MOM_tracer_Z_init
