! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the routines used to apply incremental updates
!! from data assimilation.
!
!! Applying incremental updates requires the following:
!! 1. initialize_oda_incupd_fixed and initialize_oda_incupd
!! 2. set_up_oda_incupd_field (tracers) and set_up_oda_incupd_vel_field (vel)
!! 3. calc_oda_increments (if using full fields input)
!! 4. apply_oda_incupd
!! 5. output_oda_incupd_inc (output increment if using full fields input)
!! 6. init_oda_incupd_diags (to output increments in diagnostics)
!! 7. oda_incupd_end (not being used for now)

module MOM_oda_incupd


use MOM_array_transform, only : rotate_array
use MOM_coms,            only : sum_across_PEs
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator,   only : diag_ctrl
use MOM_domains,         only : pass_var,pass_vector
use MOM_error_handler,   only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_get_input,       only : directories, Get_MOM_input
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : vardesc, var_desc
use MOM_remapping,       only : remapping_cs, remapping_core_h, initialize_remapping
use MOM_remapping,       only : remappingSchemesDoc
use MOM_restart,         only : register_restart_field, register_restart_pair, MOM_restart_CS
use MOM_restart,         only : restart_init, save_restart, query_initialized
use MOM_spatial_means,   only : global_i_mean
use MOM_time_manager,    only : time_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type, get_thickness_units

implicit none ; private

#include <MOM_memory.h>


!  Publicly available functions
public set_up_oda_incupd_field, set_up_oda_incupd_vel_field
public initialize_oda_incupd_fixed, initialize_oda_incupd, apply_oda_incupd, oda_incupd_end
public init_oda_incupd_diags,calc_oda_increments,output_oda_incupd_inc

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays with extra gridding information
type :: p3d
  integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field.
  real, dimension(:,:,:), pointer :: mask_in => NULL() !< pointer to the data mask (perhaps unused) [nondim]
  real, dimension(:,:,:), pointer :: p => NULL() !< pointer to the data, in units that depend
                                                 !! on the field it refers to [various].
  real, dimension(:,:,:), pointer :: h => NULL() !< pointer to the data grid (perhaps unused)
                                                 !! in [H ~> m or kg m-2]
end type p3d

!> oda incupd control structure
type, public :: oda_incupd_CS ; private
  integer :: nz        !< The total number of layers.
  integer :: nz_data   !< The total number of arbritary layers (used by older code).
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_oda_incupd_field

  type(p3d) :: Inc(MAX_FIELDS_)      !< The increments to be applied to the field
  type(p3d) :: Inc_u  !< The increments to be applied to the u-velocities, with data in [L T-1 ~> m s-1]
  type(p3d) :: Inc_v  !< The increments to be applied to the v-velocities, with data in [L T-1 ~> m s-1]
  type(p3d) :: Ref_h  !< Vertical grid on which the increments are provided, with data in [H ~> m or kg m-2]


  integer :: nstep_incupd          !< number of time step for full update
  real    :: ncount = 0.0          !< increment time step counter [nondim].  This could be an integer
                                   !! but a real variable works better with the existing restarts.
  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays
  logical :: incupdDataOngrid  !< True if the incupd data are on the model horizontal grid
  logical :: uv_inc      !< use u and v increments

  ! for diagnostics
  type(diag_ctrl), pointer           :: diag !<structure to regulate output
  ! diagnostic for inc. fields
  integer :: id_u_oda_inc = -1 !< diagnostic id for zonal velocity inc.
  integer :: id_v_oda_inc = -1 !< diagnostic id for meridional velocity inc.
  integer :: id_h_oda_inc = -1 !< diagnostic id for layer thicknesses inc.
  integer :: id_T_oda_inc = -1 !< diagnostic id for temperature inc.
  integer :: id_S_oda_inc = -1 !< diagnostic id for salinity inc.

end type oda_incupd_CS


  interface
module subroutine initialize_oda_incupd_fixed( G, GV, US, CS, restart_CS)
  type(ocean_grid_type),   intent(in)    :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US          !< A dimensional unit scaling type
  type(oda_incupd_CS),     pointer       :: CS          !< A pointer that is set to point to the control
                                                        !! structure for this module (in/out).
  type(MOM_restart_CS),    intent(inout) :: restart_CS  !< MOM restart control struct

end subroutine initialize_oda_incupd_fixed
module subroutine initialize_oda_incupd( G, GV, US, param_file, CS, data_h, nz_data, restart_CS)
  type(ocean_grid_type),      intent(in) :: G           !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),      intent(in) :: US          !< A dimensional unit scaling type
  integer,                    intent(in) :: nz_data     !< The total number of incr. input layers.
  type(param_file_type),      intent(in) :: param_file  !< A structure indicating the open file
                                                        !! to parse for model parameter values.
  type(oda_incupd_CS),        pointer    :: CS          !< A pointer that is set to point to the control
                                                        !! structure for this module (in/out).
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(in) :: data_h !< The ODA h
                                                                 !! [H ~> m or kg m-2].
  type(MOM_restart_CS),       intent(in) :: restart_CS  !< MOM restart control struct

  ! This include declares and sets the variable "version".

end subroutine initialize_oda_incupd
module subroutine set_up_oda_incupd_field(sp_val, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in) :: GV     !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS     !< oda_incupd control structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), &
                           intent(in) :: sp_val !< increment field, it can have an arbitrary number
                                                !! of layers, in various units depending on the
                                                !! field it refers to [various].


end subroutine set_up_oda_incupd_field
module subroutine set_up_oda_incupd_vel_field(u_val, v_val, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS !< oda incupd structure (in/out).

  real, dimension(SZIB_(G),SZJ_(G),CS%nz_data), &
                          intent(in) :: u_val !< u increment, it has arbritary number of layers but
                                              !! not to exceed the total number of model layers [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),CS%nz_data), &
                          intent(in) :: v_val !< v increment, it has arbritary number of layers but
                                              !! not to exceed the number of model layers [L T-1 ~> m s-1]

end subroutine set_up_oda_incupd_vel_field
module subroutine calc_oda_increments(h, tv, u, v, G, GV, US, CS)

  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(thermo_var_ptrs),     intent(in)    :: tv !< A structure pointing to various thermodynamic variables

  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)      :: u    !< The zonal velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)      :: v    !< The meridional velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  type(oda_incupd_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_oda_incupd (in).


                                               ! like [S ~> ppt] for salinity.
                                               ! like [S ~> ppt] for salinity.



end subroutine calc_oda_increments
module subroutine apply_oda_incupd(h, tv, u, v, dt, G, GV, US, CS)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout)    :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(thermo_var_ptrs),     intent(inout) :: tv !< A structure pointing to various thermodynamic variables

  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout)   :: u    !< The zonal velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout)   :: v    !< The meridional velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]

  real,                      intent(in)    :: dt !< The amount of time covered by this call [T ~> s].
  type(oda_incupd_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_oda_incupd (in).

  ! Local variables
                                               ! like [S ~> ppt] for salinity.
                                               ! like [S ~> ppt] for salinity.



!  integer :: ncount      ! time step counter

end subroutine apply_oda_incupd
module subroutine output_oda_incupd_inc(Time, G, GV, param_file, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating the open file
                                                             !! to parse for
                                                             !model parameter
                                                             !values.
  type(oda_incupd_CS),     pointer       :: CS   !< ODA incupd control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling




end subroutine output_oda_incupd_inc
module subroutine init_oda_incupd_diags(Time, G, GV, diag, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(oda_incupd_CS),     pointer       :: CS   !< ALE sponge control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling

end subroutine init_oda_incupd_diags
module subroutine oda_incupd_end(CS)
  type(oda_incupd_CS), pointer :: CS !< A pointer to the control structure that is
                                     !! set by a previous call to initialize_oda_incupd.


end subroutine oda_incupd_end
  end interface

end module MOM_oda_incupd
