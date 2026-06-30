! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the routines used to apply sponge layers when using
!! the ALE mode.
!!
!! Applying sponges requires the following:
!! 1. initialize_ALE_sponge
!! 2. set_up_ALE_sponge_field (tracers) and set_up_ALE_sponge_vel_field (vel)
!! 3. apply_ALE_sponge
!! 4. init_ALE_sponge_diags (not being used for now)
!! 5. ALE_sponge_end (not being used for now)

module MOM_ALE_sponge


use MOM_array_transform, only: rotate_array
use MOM_coms,          only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains,       only : pass_var, To_ALL, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer
use MOM_interpolate,   only : init_external_field, time_interp_external_init
use MOM_interpolate,   only : get_external_field_info
use MOM_interpolate,   only : external_field
use MOM_io,            only : axis_info
use MOM_remapping,     only : remapping_cs, remapping_core_h, initialize_remapping
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Store the reference profile at h points for a variable
interface set_up_ALE_sponge_field
  module procedure set_up_ALE_sponge_field_fixed
  module procedure set_up_ALE_sponge_field_varying
end interface

!> This subroutine stores the reference profile at u and v points for a vector
interface set_up_ALE_sponge_vel_field
  module procedure set_up_ALE_sponge_vel_field_fixed
  module procedure set_up_ALE_sponge_vel_field_varying
end interface

!> Determine the number of points which are within sponges in this computational domain.
!!
!! Only points that have positive values of Iresttime and which mask2dT indicates are ocean
!! points are included in the sponges.  It also stores the target interface heights.
interface initialize_ALE_sponge
  module procedure initialize_ALE_sponge_fixed
  module procedure initialize_ALE_sponge_varying
end interface

!  Publicly available functions
public set_up_ALE_sponge_field, set_up_ALE_sponge_vel_field
public get_ALE_sponge_thicknesses, get_ALE_sponge_nz_data
public initialize_ALE_sponge, apply_ALE_sponge, ALE_sponge_end, init_ALE_sponge_diags
public rotate_ALE_sponge, update_ALE_sponge_field

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays with extra gridding information
type :: p3d ; private
  !integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field.
  integer :: num_tlevs !< The number of time records contained in the file
  real, dimension(:,:,:), pointer :: p => NULL() !< pointer to the data [various]
  real, dimension(:,:,:), pointer :: dz => NULL() !< pointer to the data grid spacing [Z ~> m]
end type p3d

!> A structure for creating arrays of pointers to 2D arrays with extra gridding information
type :: p2d ; private
  type(external_field) :: field !< Time interpolator field handle
  integer :: nz_data !< The number of vertical levels in the input field
  integer :: num_tlevs !< The number of time records contained in the file
  real :: scale = 1.0  !< A multiplicative factor by which to rescale input data [various]
  real, dimension(:,:), pointer :: p => NULL() !< pointer to the data [various]
  real, dimension(:,:), pointer :: dz => NULL() !< pointer to the data grid spacing [Z ~> m]
  character(len=:), allocatable  :: name  !< The name of the input field
  character(len=:), allocatable  :: long_name !< The long name of the input field
  character(len=:), allocatable  :: unit !< The unit of the input field
  type(axis_info),  allocatable  :: axes_data(:) !< Axis types for the input field
end type p2d

!> ALE sponge control structure
type, public :: ALE_sponge_CS ; private
  integer :: nz        !< The total number of layers.
  integer :: nz_data   !< The total number of arbitrary layers (used by older code).
  integer :: num_col   !< The number of sponge tracer points within the computational domain.
  integer :: num_col_u !< The number of sponge u-points within the computational domain.
  integer :: num_col_v !< The number of sponge v-points within the computational domain.
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_sponge_field
  logical :: sponge_uv !< Control whether u and v are included in sponge
  integer, allocatable :: col_i(:)    !< Array of the i-indices of each tracer column being damped
  integer, allocatable :: col_j(:)    !< Array of the j-indices of each tracer column being damped
  integer, allocatable :: col_i_u(:)  !< Array of the i-indices of each u-column being damped
  integer, allocatable :: col_j_u(:)  !< Array of the j-indices of each u-column being damped
  integer, allocatable :: col_i_v(:)  !< Array of the i-indices of each v-column being damped
  integer, allocatable :: col_j_v(:)  !< Array of the j-indices of each v-column being damped

  real, allocatable :: Iresttime_col(:)   !< The inverse restoring time of each tracer column [T-1 ~> s-1]
  real, allocatable :: Iresttime_col_u(:) !< The inverse restoring time of each u-column [T-1 ~> s-1]
  real, allocatable :: Iresttime_col_v(:) !< The inverse restoring time of each v-column [T-1 ~> s-1]

  type(p3d) :: var(MAX_FIELDS_)      !< Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_) !< The values to which the fields are damped.
  type(p2d) :: Ref_val_u  !< The values to which the u-velocities are damped.
  type(p2d) :: Ref_val_v  !< The values to which the v-velocities are damped.
  type(p3d) :: var_u  !< Pointer to the u velocities that are being damped.
  type(p3d) :: var_v  !< Pointer to the v velocities that are being damped.
  type(p2d) :: Ref_dz !< Grid on which reference data is provided (older code).
  type(p2d) :: Ref_dzu !< u-point grid on which reference data is provided (older code).
  type(p2d) :: Ref_dzv !< v-point grid on which reference data is provided (older code).

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays
  integer :: remap_answer_date     !< The vintage of the order of arithmetic and expressions to use
                                   !! for remapping.  Values below 20190101 recover the remapping
                                   !! answers from 2018, while higher values use more robust
                                   !! forms of the same remapping expressions.
  integer :: hor_regrid_answer_date !< The vintage of the order of arithmetic and expressions to use
                                   !! for horizontal regridding.  Values below 20190101 recover the
                                   !! answers from 2018, while higher values use expressions that have
                                   !! been rearranged for rotational invariance.

  logical :: time_varying_sponges  !< True if using newer sponge code
  logical :: spongeDataOngrid      !< True if the sponge data are on the model horizontal grid
  real :: varying_input_dz_mask    !< An input file thickness below which the target values with time-varying
                                   !! sponges are replaced by the value above [Z ~> m].
                                   !! It is not clear why this needs to be greater than 0.

  !>@{ Diagnostic IDs
  integer, dimension(MAX_FIELDS_) :: id_sp_tendency = reshape([-1], [MAX_FIELDS_], [-1]) !< Diagnostic ids for tracer
                                                                                         !! tendencies due to sponges.
                                                                                         !! Init all to -1.
  integer :: id_sp_u_tendency                  !< Diagnostic id for zonal momentum tendency due to
                                               !! Rayleigh damping
  integer :: id_sp_v_tendency                  !< Diagnostic id for meridional momentum tendency due to
                                               !! Rayleigh damping
end type ALE_sponge_CS


  interface
module subroutine initialize_ALE_sponge_fixed(Iresttime, G, GV, param_file, CS, data_h, nz_data, &
                                       Iresttime_u_in, Iresttime_v_in, data_h_is_Z)

  type(ocean_grid_type),            intent(in) :: G !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in) :: GV !< ocean vertical grid structure
  integer,                          intent(in) :: nz_data !< The total number of sponge input layers.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: Iresttime !< The inverse of the restoring time [T-1 ~> s-1].
  type(param_file_type),            intent(in) :: param_file !< A structure indicating the open file
                                                             !! to parse for model parameter values.
  type(ALE_sponge_CS),              pointer    :: CS !< A pointer that is set to point to the control
                                                     !! structure for this module (in/out).
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(inout) :: data_h !< The thicknesses of the sponge
                                                     !! input layers, in [H ~> m or kg m-2] or [Z ~> m]
                                                     !! depending on data_h_is_Z.
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(in) :: Iresttime_u_in  !< The inverse of the restoring
                                                                             !! time at U-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(in) :: Iresttime_v_in  !< The inverse of the restoring
                                                                             !! time at v-points [T-1 ~> s-1].
  logical,                optional, intent(in) :: data_h_is_Z  !< If present and true data_h is already in
                                                     !! depth units.  Omitting this is the same as setting
                                                     !! it to false.

  ! Local variables
                                                     !! input layers [Z ~> m].
  ! This include declares and sets the variable "version".

end subroutine initialize_ALE_sponge_fixed
module function get_ALE_sponge_nz_data(CS)
  type(ALE_sponge_CS), intent(in) :: CS !< ALE sponge control struct
  integer :: get_ALE_sponge_nz_data  !< The number of layers in the fixed sponge data.

end function get_ALE_sponge_nz_data
module subroutine get_ALE_sponge_thicknesses(G, GV, data_h, sponge_mask, CS, data_h_in_Z)
  type(ocean_grid_type), intent(in)    :: G !< The ocean's grid structure (in).
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  real, allocatable, dimension(:,:,:), &
                         intent(inout) :: data_h !< The thicknesses of the sponge input layers expressed
                                             !! as vertical extents [Z ~> m] or in thickness units
                                             !! [H ~> m or kg m-2], depending on the value of data_h_in_Z.
  logical, dimension(SZI_(G),SZJ_(G)), &
                         intent(out)   :: sponge_mask !< A logical mask that is true where
                                                 !! sponges are being applied.
  type(ALE_sponge_CS),   pointer       :: CS !< A pointer that is set to point to the control
                                             !! structure for the ALE_sponge module.
  logical, optional,     intent(in) :: data_h_in_Z  !< If present and true data_h is returned in
                                                    !! depth units.  Omitting this is the same as setting
                                                    !! it to false.

end subroutine get_ALE_sponge_thicknesses
module subroutine initialize_ALE_sponge_varying(Iresttime, G, GV, US, param_file, CS, Iresttime_u_in, Iresttime_v_in)

  type(ocean_grid_type),            intent(in) :: G !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in) :: GV !< ocean vertical grid structure
  type(unit_scale_type),            intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: Iresttime !< The inverse of the restoring time [T-1 ~> s-1].
  type(param_file_type),            intent(in) :: param_file !< A structure indicating the open file to parse
                                                             !! for model parameter values.
  type(ALE_sponge_CS),              pointer    :: CS !< A pointer that is set to point to the control
                                                     !! structure for this module (in/out).
  real, dimension(SZIB_(G),SZJ_(G)), intent(in), optional :: Iresttime_u_in !< The inverse of the restoring time
                                                                            !! for u [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in), optional :: Iresttime_v_in !< The inverse of the restoring time
                                                                            !! for v [T-1 ~> s-1].

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine initialize_ALE_sponge_varying
module subroutine init_ALE_sponge_diags(Time, G, diag, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(ALE_sponge_CS),     intent(inout) :: CS   !< ALE sponge control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  ! Local Variables

end subroutine init_ALE_sponge_diags
module subroutine set_up_ALE_sponge_field_fixed(sp_val, G, GV, f_ptr, CS,  &
                                           sp_name, sp_long_name, sp_unit, scale)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(ALE_sponge_CS),     pointer    :: CS !< ALE sponge control structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), &
                           intent(in) :: sp_val !< Field to be used in the sponge, it can have an
                                            !! arbitrary number of layers [various]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr !< Pointer to the field to be damped [various]
  character(len=*),        intent(in) :: sp_name  !< The name of the tracer field
  character(len=*),        optional, &
                           intent(in) :: sp_long_name !< The long name of the tracer field
                                                      !! if not given, use the sp_name
  character(len=*),        optional, &
                           intent(in) :: sp_unit !< The unit of the tracer field
                                                 !! if not given, use the none
  real,          optional, intent(in) :: scale !< A factor by which to rescale the input data, including any
                                               !! contributions due to dimensional rescaling [various ~> 1].
                                               !! The default is 1.


end subroutine set_up_ALE_sponge_field_fixed
module subroutine set_up_ALE_sponge_field_varying(filename, fieldname, Time, G, GV, US, f_ptr, CS,  &
                                             sp_name, sp_long_name, sp_unit, scale)
  character(len=*),        intent(in) :: filename !< The name of the file with the
                                                  !! time varying field data
  character(len=*),        intent(in) :: fieldname !< The name of the field in the file
                                                  !! with the time varying field data
  type(time_type),         intent(in) :: Time  !< The current model time
  type(ocean_grid_type),   intent(in) :: G     !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr !< Pointer to the field to be damped (in) [various].
  type(ALE_sponge_CS),     pointer    :: CS    !< Sponge control structure (in/out).
  character(len=*),        intent(in) :: sp_name  !< The name of the tracer field
  character(len=*),        optional,  &
                           intent(in) :: sp_long_name !< The long name of the tracer field
                                                      !! if not given, use the sp_name
  character(len=*),        optional,  &
                           intent(in) :: sp_unit !< The unit of the tracer field
                                                 !! if not given, use 'none'
  real,          optional, intent(in) :: scale !< A factor by which to rescale the input data, including any
                                               !! contributions due to dimensional rescaling [various ~> 1].

  ! Local variables
end subroutine set_up_ALE_sponge_field_varying
module subroutine set_up_ALE_sponge_vel_field_fixed(u_val, v_val, G, GV, u_ptr, v_ptr, CS, scale)
  type(ocean_grid_type),   intent(in) :: G     !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(ALE_sponge_CS),     pointer    :: CS    !< Sponge structure (in/out).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u_val !< u field to be used in the sponge [L T-1 ~> m s-1],
                                               !! it is provided on its own vertical grid that may
                                               !! have fewer layers than the model itself, but not more.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v_val !< v field to be used in the sponge [L T-1 ~> m s-1],
                                               !! it is provided on its own vertical grid that may
                                               !! have fewer layers than the model itself, but not more.
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_ptr !< u-field to be damped [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v_ptr !< v-field to be damped [L T-1 ~> m s-1]
  real,          optional, intent(in) :: scale !< A factor by which to rescale the input data, including any
                                               !! contributions due to dimensional rescaling [various ~> 1].
                                               !! The default is 1.


end subroutine set_up_ALE_sponge_vel_field_fixed
module subroutine set_up_ALE_sponge_vel_field_varying(filename_u, fieldname_u, filename_v, fieldname_v, &
                                               Time, G, GV, US, CS, u_ptr, v_ptr, scale)
  character(len=*), intent(in)    :: filename_u  !< File name for u field
  character(len=*), intent(in)    :: fieldname_u !< Name of u variable in file
  character(len=*), intent(in)    :: filename_v  !< File name for v field
  character(len=*), intent(in)    :: fieldname_v !< Name of v variable in file
  type(time_type),  intent(in)    :: Time        !< Model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid (in)
  type(verticalGrid_type), intent(in) :: GV      !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(ALE_sponge_CS),     pointer    :: CS      !< Sponge structure (in/out).
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_ptr !< u-field to be damped [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v_ptr !< v-field to be damped [L T-1 ~> m s-1]
  real,          optional, intent(in) :: scale   !< A factor by which to rescale the input data, including any
                                                 !! contributions due to dimensional rescaling, often in
                                                 !! [L s T-1 m-1 ~> 1].  For varying velocities the
                                                 !! default is the same as using US%m_s_to_L_T.

  ! Local variables
end subroutine set_up_ALE_sponge_vel_field_varying
module subroutine apply_ALE_sponge(h, tv, dt, G, GV, US, CS, Time)
  type(ocean_grid_type),     intent(inout) :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(thermo_var_ptrs),     intent(in)    :: tv !< A structure pointing to various
                                                 !! thermodynamic variables
  real,                      intent(in)    :: dt !< The amount of time covered by this call [T ~> s].
  type(ALE_sponge_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_ALE_sponge (in).
  type(time_type),           intent(in)    :: Time !< The current model date

  ! Local variables
                                                !! diagnostics [various] then in [various T-1 ~> various s-1]
                                                !! first in [L T-1 ~> m s-1] then in [L T-2 ~> m s-2]
                                                !! first in [L T-1 ~> m s-1] then in [L T-2 ~> m s-2]

  ! Local variables for ALE remapping
                                                        ! edges in the input file [Z ~> m]

end subroutine apply_ALE_sponge
module subroutine rotate_ALE_sponge(sponge_in, G_in, sponge, G, GV, US, turns, param_file)
  type(ALE_sponge_CS),     intent(in) :: sponge_in !< The control structure for this module with the
                                                   !! original grid rotation
  type(ocean_grid_type),   intent(in) :: G_in      !< The ocean's grid structure with the original rotation.
  type(ALE_sponge_CS),     pointer    :: sponge    !< A pointer to the control that will be set up with
                                                   !! the new grid rotation
  type(ocean_grid_type),   intent(in) :: G         !< The ocean's grid structure with the new rotation.
  type(verticalGrid_type), intent(in) :: GV        !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US        !< A dimensional unit scaling type
  integer,                 intent(in) :: turns     !< The number of 90-degree turns between grids
  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file
                                                   !! to parse for model parameter values.

  ! First part: Index construction
  !   1. Reconstruct Iresttime(:,:) from sponge_in
  !   2. rotate Iresttime(:,:)
  !   3. Call initialize_ALE_sponge using new grid and rotated Iresttime(:,:)
  ! All the index adjustment should follow from the Iresttime rotation


end subroutine rotate_ALE_sponge
module subroutine update_ALE_sponge_field(sponge, p_old, G, GV, p_new)
  type(ALE_sponge_CS),     intent(inout) :: sponge !< ALE sponge control struct
  real, dimension(:,:,:), &
                   target, intent(in) :: p_old !< The previous array of target values [various]
  type(ocean_grid_type),   intent(in) :: G     !< The updated ocean grid structure
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: p_new !< The new array of target values [various]


end subroutine update_ALE_sponge_field
module subroutine ALE_sponge_end(CS)
  type(ALE_sponge_CS), pointer :: CS !< A pointer to the control structure that is
                                     !! set by a previous call to initialize_ALE_sponge.


end subroutine ALE_sponge_end
  end interface

end module MOM_ALE_sponge
