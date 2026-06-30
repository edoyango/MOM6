! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Horizontal interpolation
module MOM_horizontal_regridding

use MOM_debugging,     only : hchksum
use MOM_coms,          only : max_across_PEs, min_across_PEs, sum_across_PEs, broadcast
use MOM_coms,          only : reproducing_sum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_LOOP
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_error_handler, only : MOM_get_verbosity
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_interpolate,   only : time_interp_external
use MOM_interp_infra,  only : run_horiz_interp, build_horiz_interp_weights
use MOM_interp_infra,  only : horiz_interp_type, horizontal_interp_init
use MOM_interpolate,   only : get_external_field_info
use MOM_interp_infra,  only : external_field
use MOM_time_manager,  only : time_type
use MOM_io,            only : axis_info, get_axis_info, get_var_axes_info, MOM_read_data
use MOM_io,            only : read_attribute, read_variable

implicit none ; private

#include <MOM_memory.h>

public :: horiz_interp_and_extrap_tracer, myStats, homogenize_field

!> Extrapolate and interpolate data
interface horiz_interp_and_extrap_tracer
  module procedure horiz_interp_and_extrap_tracer_record
  module procedure horiz_interp_and_extrap_tracer_fms_id
end interface

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].


  interface
module subroutine myStats(array, missing, G, k, mesg, unscale, full_halo)
  type(ocean_grid_type), intent(in) :: G     !< Ocean grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in) :: array !< input array in arbitrary units [A ~> a]
  real,                  intent(in) :: missing !< missing value in arbitrary units [A ~> a]
  integer,               intent(in) :: k     !< Level to calculate statistics for
  character(len=*),      intent(in) :: mesg  !< Label to use in message
  real,        optional, intent(in) :: unscale !< A scaling factor for output that countacts
                                             !! any internal dimesional scaling [a A-1 ~> 1]
  logical,     optional, intent(in) :: full_halo !< If present and true, test values on the whole
                                             !! array rather than just the computational domain.
  ! Local variables

end subroutine myStats
module subroutine fill_miss_2d(aout, good, fill, prev, G, acrit, num_pass, relc, debug, answer_date)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(inout) :: aout !< The array with missing values to fill [arbitrary]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: good !< Valid data mask for incoming array
                                               !! (1==good data; 0==missing data) [nondim].
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: fill !< Same shape array of points which need
                                               !! filling (1==fill;0==dont fill) [nondim]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: prev !< First guess where isolated holes exist [arbitrary]
  real,                  intent(in)    :: acrit !< A minimal value for deltas between iterations that
                                               !! determines when the smoothing has converged [arbitrary].
  integer,     optional, intent(in)    :: num_pass !< The maximum number of iterations
  real,        optional, intent(in)    :: relc !< A relaxation coefficient for Laplacian [nondim]
  logical,     optional, intent(in)    :: debug !< If true, write verbose debugging messages.
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                !! Dates before 20190101 give the same  answers
                                                !! as the code did in late 2018, while later versions
                                                !! add parentheses for rotational symmetry.




end subroutine fill_miss_2d
module subroutine horiz_interp_and_extrap_tracer_record(filename, varnam, recnum, G, tr_z, mask_z, &
                                                 z_in, z_edges_in, missing_value, scale, &
                                                 homogenize, m_to_Z, answers_2018, ongrid, tr_iter_tol, answer_date)

  character(len=*),      intent(in)    :: filename   !< Path to file containing tracer to be
                                                     !! interpolated.
  character(len=*),      intent(in)    :: varnam     !< Name of tracer in file.
  integer,               intent(in)    :: recnum     !< Record number of tracer to be read.
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:), intent(out) :: tr_z
                                                     !< Allocatable tracer array on the horizontal
                                                     !! model grid and input-file vertical levels
                                                     !! in arbitrary units [A ~> a]
  real, allocatable, dimension(:,:,:), intent(out) :: mask_z
                                                     !< Allocatable tracer mask array on the horizontal
                                                     !! model grid and input-file vertical levels [nondim]
  real, allocatable, dimension(:), intent(out) :: z_in
                                                     !< Cell grid values for input data [Z ~> m]
  real, allocatable, dimension(:), intent(out) :: z_edges_in
                                                     !< Cell grid edge values for input data [Z ~> m]
  real,                  intent(out)   :: missing_value !< The missing value in the returned array, scaled
                                                     !! to avoid accidentally having valid values match
                                                     !! missing values in the same units as tr_z [A ~> a]
  real,                  intent(in)    :: scale      !< Scaling factor for tracer into the internal
                                                     !! units of the model for the units in the file [A a-1 ~> 1]
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth [Z m-1 ~> 1].  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.
  logical,     optional, intent(in)    :: ongrid     !< If true, then data are assumed to have been interpolated
                                                     !! to the model horizontal grid. In this case, only
                                                     !! extrapolation is performed by this routine
  real,        optional, intent(in)    :: tr_iter_tol !< The tolerance for changes in tracer concentrations
                                                     !! between smoothing iterations that determines when to
                                                     !! stop iterating in the same units as tr_z [A ~> a]
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                     !! Dates before 20190101 give the same  answers
                                                     !! as the code did in late 2018, while later versions
                                                     !! add parentheses for rotational symmetry.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
                                                     !! native horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
                                                     !! model horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
                                                     !! with units that change as the input data is
                                                     !! interpreted [a] then [A ~> a]

                                ! iterations that determines when to stop iterating [A ~> a]

end subroutine horiz_interp_and_extrap_tracer_record
module subroutine horiz_interp_and_extrap_tracer_fms_id(field, Time, G, tr_z, mask_z, &
                                                 z_in, z_edges_in, missing_value, scale, &
                                                 homogenize, spongeOngrid, m_to_Z, &
                                                 answers_2018, tr_iter_tol, answer_date, &
                                                 axes)

  type(external_field), intent(in)     :: field      !< Handle for the time interpolated field
  type(time_type),       intent(in)    :: Time       !< A FMS time type
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:), intent(out) :: tr_z
                                                     !< Allocatable tracer array on the horizontal
                                                     !! model grid and input-file vertical levels
                                                     !! in arbitrary units [A ~> a]
  real, allocatable, dimension(:,:,:), intent(out) :: mask_z
                                                     !< Allocatable tracer mask array on the horizontal
                                                     !! model grid and input-file vertical levels [nondim]
  real, allocatable, dimension(:), intent(out) :: z_in
                                                     !< Cell grid values for input data [Z ~> m]
  real, allocatable, dimension(:), intent(out) :: z_edges_in
                                                     !< Cell grid edge values for input data [Z ~> m]
  real,                  intent(out)   :: missing_value !< The missing value in the returned array, scaled
                                                     !! to avoid accidentally having valid values match
                                                     !! missing values, in the same arbitrary units as tr_z [A ~> a]
  real,                  intent(in)    :: scale      !< Scaling factor for tracer into the internal
                                                     !! units of the model [A a-1 ~> 1]
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  logical,     optional, intent(in)    :: spongeOngrid !< If present and true, the sponge data are on the model grid
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth [Z m-1 ~> 1].  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.
  real,        optional, intent(in)    :: tr_iter_tol !< The tolerance for changes in tracer concentrations
                                                     !! between smoothing iterations that determines when to
                                                     !! stop iterating, in the same arbitrary units as tr_z [A ~> a]
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                     !! Dates before 20190101 give the same  answers
                                                     !! as the code did in late 2018, while later versions
                                                     !! add parentheses for rotational symmetry.
  type(axis_info), allocatable, dimension(:), optional, intent(inout) :: axes !< Axis types for the input data

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
                                                     !! native horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
                                                     !! with units that change as the input data is
                                                     !! interpreted [a] then [A ~> a]
                                                     !! on the original grid [a]

                                ! iterations that determines when to stop iterating [A ~> a]

end subroutine horiz_interp_and_extrap_tracer_fms_id
module subroutine homogenize_field(field, G, tmp_scale, weights, answer_date, wt_unscale)
  type(ocean_grid_type),            intent(inout) :: G      !< Ocean grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: field  !< The tracer on the model grid in arbitrary units [A ~> a]
  real,                   optional, intent(in)    :: tmp_scale !< A temporary rescaling factor for the
                                                            !! variable that is reversed in the
                                                            !! return value [a A-1 ~> 1]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: weights !< The weights for the tracer in arbitrary units that
                                                            !! typically differ from those used by field [B ~> b]
  integer,                optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                            !! Dates before 20230101 use non-reproducing sums
                                                            !! in their averages, while later versions use
                                                            !! reproducing sums for rotational symmetry and
                                                            !! consistency across PE layouts.
  real,                   optional, intent(in)    :: wt_unscale !< A factor that undoes any dimensional scaling
                                                            !! of the weights so that they can be used with
                                                            !! reproducing sums [b B-1 ~> 1]

  ! Local variables
  ! In the following comments, [A] and [B] are used to indicate the arbitrary, possibly rescaled
  ! units of the input field and the weighting array, while [a] and [b] indicate the corresponding
  ! unscaled (e.g., mks) units that can be used with the reproducing sums
                      ! tracer-point grid mask if it weights is absent [B ~> b]

end subroutine homogenize_field
module subroutine meshgrid(x, y, x_T, y_T)
  real, dimension(:),                   intent(in)    :: x  !< input 1-dimensional vector [arbitrary]
  real, dimension(:),                   intent(in)    :: y  !< input 1-dimensional vector [arbitrary]
  real, dimension(size(x,1),size(y,1)), intent(inout) :: x_T !< output 2-dimensional array [arbitrary]
  real, dimension(size(x,1),size(y,1)), intent(inout) :: y_T !< output 2-dimensional array [arbitrary]


end subroutine meshgrid
  end interface

end module MOM_horizontal_regridding
