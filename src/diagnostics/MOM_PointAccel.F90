! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Debug accelerations at a given point
!!
!!    The two subroutines in this file write out all of the terms
!! in the u- or v-momentum balance at a given point.  Usually
!! these subroutines are called after the velocities exceed some
!! threshold, in order to determine which term is culpable.
!! often this is done for debugging purposes.
module MOM_PointAccel

use MOM_diag_mediator, only : diag_ctrl
use MOM_domains, only : pe_here
use MOM_error_handler, only : MOM_error, NOTE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : open_ASCII_file, APPEND_FILE, MULTIPLE, SINGLE_FILE
use MOM_time_manager, only : time_type, get_time, get_date, set_date, operator(-)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : ocean_internal_state, accel_diag_ptrs, cont_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public write_u_accel, write_v_accel, PointAccel_init

!> The control structure for the MOM_PointAccel module
type, public :: PointAccel_CS ; private
  character(len=200) :: u_trunc_file !< The complete path to the file in which a column's worth of
                                     !! u-accelerations are written if u-velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to the file in which a column's worth of
                                     !! v-accelerations are written if v-velocity truncations occur.
  integer :: u_file         !< The unit number for an opened u-truncation files, or -1 if it has not yet been opened.
  integer :: v_file         !< The unit number for an opened v-truncation files, or -1 if it has not yet been opened.
  integer :: cols_written   !< The number of columns whose output has been
                            !! written by this PE during the current run.
  integer :: max_writes     !< The maximum number of times any PE can write out
                            !! a column's worth of accelerations during a run.
  logical :: full_column    !< If true, write out the accelerations in all massive layers,
                            !! otherwise just document the ones with large velocities.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
! The following are pointers to many of the state variables and accelerations
! that are used to step the physical model forward.  They all use the same
! names as the variables they point to in MOM.F90
  real, pointer, dimension(:,:,:) :: &
    u_av => NULL(), &       !< Time average u-velocity [L T-1 ~> m s-1]
    v_av => NULL(), &       !< Time average velocity [L T-1 ~> m s-1]
    u_prev => NULL(), &     !< Previous u-velocity [L T-1 ~> m s-1]
    v_prev => NULL(), &     !< Previous v-velocity [L T-1 ~> m s-1]
    T => NULL(), &          !< Temperature [C ~> degC]
    S => NULL(), &          !< Salinity [S ~> ppt]
    u_accel_bt => NULL(), & !< Barotropic u-accelerations [L T-2 ~> m s-2]
    v_accel_bt => NULL()    !< Barotropic v-accelerations [L T-2 ~> m s-2]
end type PointAccel_CS


  interface
module subroutine write_u_accel(I, j, um, hin, ADp, CDp, dt, G, GV, US, CS, vel_rpt, str, a, hv)
  integer,                     intent(in) :: I   !< The zonal index of the column to be documented.
  integer,                     intent(in) :: j   !< The meridional index of the column to be documented.
  type(ocean_grid_type),       intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),     intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),       intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: um  !< The new zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: hin !< The layer thickness [H ~> m or kg m-2].
  type(accel_diag_ptrs),       intent(in) :: ADp !< A structure pointing to the various
                                                 !! accelerations in the momentum equations.
  type(cont_diag_ptrs),        intent(in) :: CDp !<  A structure with pointers to various terms
                                                 !! in the continuity equations.
  real,                        intent(in) :: dt  !< The ocean dynamics time step [T ~> s].
  type(PointAccel_CS),         pointer    :: CS  !< The control structure returned by a previous
                                                 !! call to PointAccel_init.
  real,                        intent(in) :: vel_rpt !< The velocity magnitude that triggers a report [L T-1 ~> m s-1]
  real, optional,              intent(in) :: str !< The surface wind stress [R L Z T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                     optional, intent(in) :: a   !< The layer coupling coefficients from vertvisc
                                                 !! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in) :: hv  !< The layer thicknesses at velocity grid points,
                                                 !! from vertvisc [H ~> m or kg m-2].

  ! Local variables
                              ! or [kg T m-1 s-1 L-1 H-1 ~> 1]

end subroutine write_u_accel
module subroutine write_v_accel(i, J, vm, hin, ADp, CDp, dt, G, GV, US, CS, vel_rpt, str, a, hv)
  integer,                     intent(in) :: i   !< The zonal index of the column to be documented.
  integer,                     intent(in) :: J   !< The meridional index of the column to be documented.
  type(ocean_grid_type),       intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),     intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),       intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                               intent(in) :: vm  !< The new meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: hin !< The layer thickness [H ~> m or kg m-2].
  type(accel_diag_ptrs),       intent(in) :: ADp !< A structure pointing to the various
                                                 !! accelerations in the momentum equations.
  type(cont_diag_ptrs),        intent(in) :: CDp !< A structure with pointers to various terms in
                                                 !! the continuity equations.
  real,                        intent(in) :: dt  !< The ocean dynamics time step [T ~> s].
  type(PointAccel_CS),         pointer    :: CS  !< The control structure returned by a previous
                                                 !! call to PointAccel_init.
  real,                        intent(in) :: vel_rpt !< The velocity magnitude that triggers a report [L T-1 ~> m s-1]
  real, optional,              intent(in) :: str !< The surface wind stress [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                     optional, intent(in) :: a   !< The layer coupling coefficients from vertvisc
                                                 !! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                     optional, intent(in) :: hv  !< The layer thicknesses at velocity grid points,
                                                 !! from vertvisc [H ~> m or kg m-2].

  ! Local variables
                              ! or [kg T m-1 s-1 L-1 H-1 ~> 1]

end subroutine write_v_accel
module subroutine PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS)
  type(ocean_internal_state), &
                        target, intent(in)    :: MIS  !< For "MOM Internal State" a set of pointers
                                                      !! to the fields and accelerations that make
                                                      !! up the ocean's physical state.
  type(time_type),      target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure.
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time
                                                      !! parameters.
  type(diag_ctrl),      target, intent(inout) :: diag !< A structure that is used to regulate
                                                      !! diagnostic output.
  type(directories),            intent(in)    :: dirs !< A structure containing several relevant
                                                      !! directory paths.
  type(PointAccel_CS),          pointer       :: CS   !< A pointer that is set to point to the
                                                      !! control structure for this module.
  ! This include declares and sets the variable "version".

end subroutine PointAccel_init
  end interface

end module MOM_PointAccel
