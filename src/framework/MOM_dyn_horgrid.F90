! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Contains a shareable dynamic type for describing horizontal grids and metric data
!! and utilty routines that work on this type.
module MOM_dyn_horgrid

use MOM_array_transform, only : rotate_array, rotate_array_pair
use MOM_domains,         only : MOM_domain_type, deallocate_MOM_domain
use MOM_error_handler,   only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_hor_index,       only : hor_index_type
use MOM_unit_scaling,    only : unit_scale_type

implicit none ; private

public create_dyn_horgrid, destroy_dyn_horgrid, set_derived_dyn_horgrid
public rescale_dyn_horgrid_bathymetry, rotate_dyn_horgrid

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Describes the horizontal ocean grid with only dynamic memory arrays
type, public :: dyn_horgrid_type
  type(MOM_domain_type), pointer :: Domain => NULL() !< Ocean model domain
  type(MOM_domain_type), pointer :: Domain_aux => NULL() !< A non-symmetric auxiliary domain type.
  type(hor_index_type) :: HI !< Horizontal index ranges

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

  integer :: isd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: jsd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  logical :: symmetric  !< True if symmetric memory is used.

  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                  !! allowed.  The default is .false. (for now).
  integer :: first_direction !< An integer that indicates which direction is to be updated first in
                             !! directionally split parts of the calculation.  This can be altered
                             !! during the course of the run via calls to set_first_direction.

  real, allocatable, dimension(:,:) :: &
    mask2dT, &   !< 0 for land points and 1 for ocean points on the h-grid [nondim].
    geoLatT, &   !< The geographic latitude at q points [degrees of latitude] or [m].
    geoLonT, &   !< The geographic longitude at q points [degrees of longitude] or [m].
    dxT, &       !< dxT is delta x at h points [L ~> m].
    IdxT, &      !< 1/dxT [L-1 ~> m-1].
    dyT, &       !< dyT is delta y at h points [L ~> m].
    IdyT, &      !< IdyT is 1/dyT [L-1 ~> m-1].
    areaT, &     !< The area of an h-cell [L2 ~> m2].
    IareaT       !< 1/areaT [L-2 ~> m-2].
  real, allocatable, dimension(:,:) :: sin_rot
                 !< The sine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].
  real, allocatable, dimension(:,:) :: cos_rot
                 !< The cosine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].

  real, allocatable, dimension(:,:) :: &
    mask2dCu, &  !< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    OBCmaskCu, & !< 0 for boundary or OBC points and 1 for ocean points on the u grid [nondim].
    geoLatCu, &  !< The geographic latitude at u points [degrees of latitude] or [m].
    geoLonCu, &  !< The geographic longitude at u points [degrees of longitude] or [m].
    dxCu, &      !< dxCu is delta x at u points [L ~> m].
    IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
    IdxCu_OBCmask, & !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
    dyCu, &      !< dyCu is delta y at u points [L ~> m].
    IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
    dy_Cu, &     !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    IareaCu, &   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    areaCu       !< The areas of the u-grid cells [L2 ~> m2].

  real, allocatable, dimension(:,:) :: &
    mask2dCv, &  !< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    OBCmaskCv, & !< 0 for boundary or OBC points and 1 for ocean points on the v grid [nondim].
    geoLatCv, &  !< The geographic latitude at v points [degrees of latitude] or [m].
    geoLonCv, &  !< The geographic longitude at v points [degrees of longitude] or [m].
    dxCv, &      !< dxCv is delta x at v points [L ~> m].
    IdxCv, &     !< 1/dxCv [L-1 ~> m-1].
    dyCv, &      !< dyCv is delta y at v points [L ~> m].
    IdyCv, &     !< 1/dyCv [L-1 ~> m-1].
    IdyCv_OBCmask, & !< 1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
    dx_Cv, &     !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    IareaCv, &   !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    areaCv       !< The areas of the v-grid cells [L2 ~> m2].

  real, allocatable, dimension(:,:) :: &
    porous_DminU, & !< minimum topographic height (deepest) of U-face [Z ~> m]
    porous_DmaxU, & !< maximum topographic height (shallowest) of U-face [Z ~> m]
    porous_DavgU    !< average topographic height of U-face [Z ~> m]

  real, allocatable, dimension(:,:) :: &
    porous_DminV, & !< minimum topographic height (deepest) of V-face [Z ~> m]
    porous_DmaxV, & !< maximum topographic height (shallowest) of V-face [Z ~> m]
    porous_DavgV    !< average topographic height of V-face [Z ~> m]

  real, allocatable, dimension(:,:) :: &
    mask2dBu, &  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    geoLatBu, &  !< The geographic latitude at q points [degrees of latitude] or [m].
    geoLonBu, &  !< The geographic longitude at q points [degrees of longitude] or [m].
    dxBu, &      !< dxBu is delta x at q points [L ~> m].
    IdxBu, &     !< 1/dxBu [L-1 ~> m-1].
    dyBu, &      !< dyBu is delta y at q points [L ~> m].
    IdyBu, &     !< 1/dyBu [L-1 ~> m-1].
    areaBu, &    !< areaBu is the area of a q-cell [L ~> m]
    IareaBu      !< IareaBu = 1/areaBu [L-2 ~> m-2].

  real, pointer, dimension(:) :: gridLatT => NULL()
        !< The latitude of T points for the purpose of labeling the output axes,
        !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLatT.
  real, pointer, dimension(:) :: gridLatB => NULL()
        !< The latitude of B points for the purpose of labeling the output axes,
        !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLatBu.
  real, pointer, dimension(:) :: gridLonT => NULL()
        !< The longitude of T points for the purpose of labeling the output axes,
        !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLonT.
  real, pointer, dimension(:) :: gridLonB => NULL()
        !< The longitude of B points for the purpose of labeling the output axes,
        !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLonBu.
  character(len=40) :: &
    ! Except on a Cartesian grid, these are usually some variant of "degrees".
    x_axis_units, &     !< The units that are used in labeling the x coordinate axes.
    y_axis_units, &     !< The units that are used in labeling the y coordinate axes.
    ! These are internally generated names, including "m", "km", "deg_E" and "deg_N".
    x_ax_unit_short, &  !< A short description of the x-axis units for documenting parameter units
    y_ax_unit_short     !< A short description of the y-axis units for documenting parameter units

  real, allocatable, dimension(:,:) :: &
    bathyT        !< Ocean bottom depth, referenced to a zero reference height at tracer points.
                  !! bathyT is in depth units and positive *below* the reference height [Z ~> m].
  real, allocatable, dimension(:,:) :: &
    meanSL        !< Spatially varying time mean sea level, referenced to a zero reference height
                  !! at tracer points. meanSL is in height units and positive *above* zero. It is used
                  !! a) as the height where p = p_atm or zero;
                  !! b) to calculate time mean thickness of the water column, where
                  !!    mean thickness = max(meanSL + bathyT, 0.0).
                  !! meanSL is 2D for the consideration of a domain with spatically varying mean
                  !! height, e.g. the Great Lakes system [Z ~> m].

  logical :: bathymetry_at_vel  !< If true, there are separate values for the
                  !! basin depths at velocity points.  Otherwise the effects of
                  !! of topography are entirely determined from thickness points.
  real, allocatable, dimension(:,:) :: &
    Dblock_u, &   !< Topographic depths at u-points at which the flow is blocked [Z ~> m].
    Dopen_u       !< Topographic depths at u-points at which the flow is open at width dy_Cu [Z ~> m].
  real, allocatable, dimension(:,:) :: &
    Dblock_v, &   !< Topographic depths at v-points at which the flow is blocked [Z ~> m].
    Dopen_v       !< Topographic depths at v-points at which the flow is open at width dx_Cv [Z ~> m].
  real, allocatable, dimension(:,:) :: &
    CoriolisBu, & !< The Coriolis parameter at corner points [T-1 ~> s-1].
    Coriolis2Bu   !< The square of the Coriolis parameter at corner points [T-2 ~> s-2].
  real, allocatable, dimension(:,:) :: &
    df_dx, &      !< Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].

  ! These variables are global sums that are useful for 1-d diagnostics.
  real :: areaT_global  !< Global sum of h-cell area [L2 ~> m2]
  real :: IareaT_global !< Global sum of inverse h-cell area (1/areaT_global) [L-2 ~> m-2]

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
  real :: max_depth     !< The maximum depth of the ocean [Z ~> m]
end type dyn_horgrid_type


  interface
module subroutine create_dyn_horgrid(G, HI, bathymetry_at_vel)
  type(dyn_horgrid_type), pointer, intent(inout) :: G  !< A pointer to the dynamic horizontal grid type
  type(hor_index_type),   intent(in) :: HI !< A hor_index_type for array extents
  logical,        optional, intent(in) :: bathymetry_at_vel !< If true, there are
                             !! separate values for the basin depths at velocity
                             !! points.  Otherwise the effects of topography are
                             !! entirely determined from thickness points.

  ! This subroutine allocates the lateral elements of the dyn_horgrid_type that
  ! are always used and zeros them out.

end subroutine create_dyn_horgrid
module subroutine rotate_dyn_horgrid(G_in, G, US, turns)
  type(dyn_horgrid_type), intent(in)    :: G_in   !< The input horizontal grid type
  type(dyn_horgrid_type), intent(inout) :: G      !< An output rotated horizontal grid type
                                                  !! that has already been allocated, but whose
                                                  !! contents are largely replaced here.
  type(unit_scale_type),  intent(in)    :: US     !< A dimensional unit scaling type
  integer, intent(in) :: turns                    !< Number of quarter turns

  ! Center point
end subroutine rotate_dyn_horgrid
module subroutine rescale_dyn_horgrid_bathymetry(G, m_in_new_units)
  type(dyn_horgrid_type), intent(inout) :: G !< The dynamic horizontal grid type
  real,                   intent(in)    :: m_in_new_units !< The new internal representation of 1 m depth [m Z-1 ~> 1]

  ! Local variables

end subroutine rescale_dyn_horgrid_bathymetry
module subroutine set_derived_dyn_horgrid(G, US)
  type(dyn_horgrid_type), intent(inout) :: G !< The dynamic horizontal grid type
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type
!    Various inverse grid spacings and derived areas are calculated within this
!  subroutine.

end subroutine set_derived_dyn_horgrid
module function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted in arbitrary units [A ~> a]
  real :: I_val            !< The Adcroft reciprocal of val [A-1 ~> a-1].

end function Adcroft_reciprocal
module subroutine destroy_dyn_horgrid(G)
  type(dyn_horgrid_type), pointer :: G !< The dynamic horizontal grid type

end subroutine destroy_dyn_horgrid
  end interface

end module MOM_dyn_horgrid
