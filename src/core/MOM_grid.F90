! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides the ocean grid type
module MOM_grid

use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_domains, only : MOM_domain_type, get_domain_extent, compute_block_extent
use MOM_domains, only : get_global_shape, deallocate_MOM_domain
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public MOM_grid_init, MOM_grid_end, set_derived_metrics, set_first_direction
public isPointInCell, hor_index_type, get_global_grid_size

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Ocean grid type. See mom_grid for details.
type, public :: ocean_grid_type
  type(MOM_domain_type), pointer :: Domain => NULL() !< Ocean model domain
  type(MOM_domain_type), pointer :: Domain_aux => NULL() !< A non-symmetric auxiliary domain type.
  type(hor_index_type) :: HI !< Horizontal index ranges
  type(hor_index_type), allocatable :: HId(:) !< Horizontal index ranges for downsampling

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: isd_global !< The value of isd in the global index space (decomposition invariant).
  integer :: jsd_global !< The value of isd in the global index space (decomposition invariant).
  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  integer :: ke         !< The number of layers in the vertical.
  logical :: symmetric  !< True if symmetric memory is used.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                  !! allowed.  The default is .false. (for now).
  integer :: first_direction !< An integer that indicates which direction is
                             !! to be updated first in directionally split
                             !! parts of the calculation.  This can be altered
                             !! during the course of the run via calls to
                             !! set_first_direction.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    mask2dT, &   !< 0 for land points and 1 for ocean points on the h-grid [nondim].
    geoLatT, &   !< The geographic latitude at tracer (h) points [degrees_N] or [km] or [m]
    geoLonT, &   !< The geographic longitude at tracer (h) points [degrees_E] or [km] or [m]
    dxT, &       !< dxT is delta x at h points [L ~> m].
    IdxT, &      !< 1/dxT [L-1 ~> m-1].
    dyT, &       !< dyT is delta y at h points [L ~> m].
    IdyT, &      !< IdyT is 1/dyT [L-1 ~> m-1].
    areaT, &     !< The area of an h-cell [L2 ~> m2].
    IareaT, &    !< 1/areaT [L-2 ~> m-2].
    sin_rot, &   !< The sine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].
    cos_rot      !< The cosine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    mask2dCu, &  !< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    OBCmaskCu, & !< 0 for boundary or OBC points and 1 for ocean points on the u grid [nondim].
    geoLatCu, &  !< The geographic latitude at u points [degrees_N] or [km] or [m]
    geoLonCu, &  !< The geographic longitude at u points [degrees_E] or [km] or [m].
    dxCu, &      !< dxCu is delta x at u points [L ~> m].
    IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
    IdxCu_OBCmask, & !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
    dyCu, &      !< dyCu is delta y at u points [L ~> m].
    IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
    dy_Cu, &     !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    IareaCu, &   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    areaCu       !< The areas of the u-grid cells [L2 ~> m2].

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    mask2dCv, &  !< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    OBCmaskCv, & !< 0 for boundary or OBC points and 1 for ocean points on the v grid [nondim].
    geoLatCv, &  !< The geographic latitude at v points [degrees_N] or [km] or [m]
    geoLonCv, &  !< The geographic longitude at v points [degrees_E] or [km] or [m].
    dxCv, &      !< dxCv is delta x at v points [L ~> m].
    IdxCv, &     !< 1/dxCv [L-1 ~> m-1].
    dyCv, &      !< dyCv is delta y at v points [L ~> m].
    IdyCv, &     !< 1/dyCv [L-1 ~> m-1].
    IdyCv_OBCmask, & !< 1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
    dx_Cv, &     !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    IareaCv, &   !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    areaCv       !< The areas of the v-grid cells [L2 ~> m2].

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    porous_DminU, & !< minimum topographic height (deepest) of U-face [Z ~> m]
    porous_DmaxU, & !< maximum topographic height (shallowest) of U-face [Z ~> m]
    porous_DavgU    !< average topographic height of U-face [Z ~> m]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    porous_DminV, & !< minimum topographic height (deepest) of V-face [Z ~> m]
    porous_DmaxV, & !< maximum topographic height (shallowest) of V-face [Z ~> m]
    porous_DavgV    !< average topographic height of V-face [Z ~> m]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    mask2dBu, &  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    geoLatBu, &  !< The geographic latitude at q points [degrees_N] or [km] or [m]
    geoLonBu, &  !< The geographic longitude at q points [degrees_E] or [km] or [m].
    dxBu, &      !< dxBu is delta x at q points [L ~> m].
    IdxBu, &     !< 1/dxBu [L-1 ~> m-1].
    dyBu, &      !< dyBu is delta y at q points [L ~> m].
    IdyBu, &     !< 1/dyBu [L-1 ~> m-1].
    areaBu, &    !< areaBu is the area of a q-cell [L2 ~> m2]
    IareaBu      !< IareaBu = 1/areaBu [L-2 ~> m-2].

  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & !< The latitude of T points for the purpose of labeling the output axes,
                          !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLatT.
    gridLatB => NULL()    !< The latitude of B points for the purpose of labeling the output axes,
                          !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLatBu.
  real, pointer, dimension(:) :: &
    gridLonT => NULL(), & !< The longitude of T points for the purpose of labeling the output axes,
                          !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLonT.
    gridLonB => NULL()    !< The longitude of B points for the purpose of labeling the output axes,
                          !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLonBu.
  character(len=40) :: &
    ! Except on a Cartesian grid, these are usually some variant of "degrees".
    x_axis_units, &     !< The units that are used in labeling the x coordinate axes.
    y_axis_units, &     !< The units that are used in labeling the y coordinate axes.
    ! These are internally generated names, including "m", "km", "deg_E" and "deg_N".
    x_ax_unit_short, &  !< A short description of the x-axis units for documenting parameter units
    y_ax_unit_short     !< A short description of the y-axis units for documenting parameter units

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    bathyT           !< Ocean bottom depth, referenced to Z_ref at tracer points. bathyT is in
                     !! depth units and positive *below* Z_ref [Z ~> m].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    meanSL           !< Spatially varying time mean sea level, referenced to Z_ref at tracer points.
                     !! meanSL is in height units and positive *above* Z_ref. It is used
                     !! a) as the height where p = p_atm or zero;
                     !! b) to calculate time mean thickness of the water column, where
                     !!    mean thickness = max(meanSL + bathyT, 0.0).
                     !! meanSL is 2D for the consideration of a domain with spatically varying mean
                     !! height, e.g. the Great Lakes system [Z ~> m].
  real    :: Z_ref   !< A reference value for all geometric height fields, such as bathyT [Z ~> m].

  logical :: bathymetry_at_vel  !< If true, there are separate values for the
                  !! basin depths at velocity points.  Otherwise the effects of
                  !! of topography are entirely determined from thickness points.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Dblock_u, &   !< Topographic depths at u-points at which the flow is blocked [Z ~> m].
    Dopen_u       !< Topographic depths at u-points at which the flow is open at width dy_Cu [Z ~> m].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Dblock_v, &   !< Topographic depths at v-points at which the flow is blocked [Z ~> m].
    Dopen_v       !< Topographic depths at v-points at which the flow is open at width dx_Cv [Z ~> m].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    CoriolisBu, & !< The Coriolis parameter at corner points [T-1 ~> s-1].
    Coriolis2Bu   !< The square of the Coriolis parameter at corner points [T-2 ~> s-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    df_dx, &      !< Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].

  ! These variables are global sums that are useful for 1-d diagnostics.
  real :: areaT_global  !< Global sum of h-cell area [L2 ~> m2]
  real :: IareaT_global !< Global sum of inverse h-cell area (1/areaT_global) [L-2 ~> m-2].

  type(unit_scale_type), pointer :: US => NULL() !< A dimensional unit scaling type


  ! These variables are for block structures.
  integer :: nblocks  !< The number of sub-PE blocks on this PE
  type(hor_index_type), pointer :: Block(:) => NULL() !< Index ranges for each block

  ! These parameters are run-time parameters that are used during some
  ! initialization routines (but not all)
  real :: grid_unit_to_L !< A factor that converts a the geoLat and geoLon variables and related
                        !! variables like len_lat and len_lon into rescaled horizontal distance
                        !! units on a Cartesian grid, in [L km ~> 1000] or [L m-1 ~> 1] or
                        !! is 0 for a non-Cartesian grid.
  real :: south_lat     !< The latitude (or y-coordinate) of the first v-line [degrees_N] or [km] or [m]
  real :: west_lon      !< The longitude (or x-coordinate) of the first u-line [degrees_E] or [km] or [m]
  real :: len_lat       !< The latitudinal (or y-coord) extent of physical domain [degrees_N] or [km] or [m]
  real :: len_lon       !< The longitudinal (or x-coord) extent of physical domain [degrees_E] or [km] or [m]
  real :: Rad_Earth_L   !< The radius of the planet in rescaled units [L ~> m]
  real :: max_depth     !< The maximum depth of the ocean in depth units [Z ~> m]
end type ocean_grid_type


  interface
module subroutine MOM_grid_init(G, param_file, US, HI, global_indexing, bathymetry_at_vel)
  type(ocean_grid_type), intent(inout) :: G          !< The horizontal grid type
  type(param_file_type), intent(in)    :: param_file !< Parameter file handle
  type(unit_scale_type), optional, pointer :: US !< A dimensional unit scaling type
  type(hor_index_type), &
                  optional, intent(in) :: HI !< A hor_index_type for array extents
  logical,        optional, intent(in) :: global_indexing !< If true use global index
                             !! values instead of having the data domain on each
                             !! processor start at 1.
  logical,        optional, intent(in) :: bathymetry_at_vel !< If true, there are
                             !! separate values for the ocean bottom depths at
                             !! velocity points.  Otherwise the effects of topography
                             !! are entirely determined from thickness points.

  ! Local variables
                             ! the data domain on each processor start at 1.
  ! This include declares and sets the variable "version".


end subroutine MOM_grid_init
module subroutine set_derived_metrics(G, US)
  type(ocean_grid_type), intent(inout) :: G  !< The horizontal grid structure
  type(unit_scale_type), intent(in)    :: US !< A dimensional unit scaling type
!    Various inverse grid spacings and derived areas are calculated within this
!  subroutine.

end subroutine set_derived_metrics
module function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted [A].
  real :: I_val            !< The Adcroft reciprocal of val [A-1].

end function Adcroft_reciprocal
logical module function isPointInCell(G, i, j, x, y)
  type(ocean_grid_type), intent(in) :: G !< Grid type
  integer,               intent(in) :: i !< i index of cell to test
  integer,               intent(in) :: j !< j index of cell to test
  real,                  intent(in) :: x !< x coordinate of point [degrees_E]
  real,                  intent(in) :: y !< y coordinate of point [degrees_N]
  ! Local variables
end function isPointInCell
module subroutine set_first_direction(G, y_first)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  integer,               intent(in) :: y_first !< The first direction to store

end subroutine set_first_direction
module subroutine get_global_grid_size(G, niglobal, njglobal)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type
  integer,               intent(out)   :: niglobal !< i-index global size of grid
  integer,               intent(out)   :: njglobal !< j-index global size of grid

end subroutine get_global_grid_size
module subroutine allocate_metrics(G)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type

  ! This subroutine allocates the lateral elements of the ocean_grid_type that
  ! are always used and zeros them out.

end subroutine allocate_metrics
module subroutine MOM_grid_end(G)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type

end subroutine MOM_grid_end
  end interface

end module MOM_grid
