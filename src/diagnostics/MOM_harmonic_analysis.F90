! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Inline harmonic analysis (conventional)
module MOM_harmonic_analysis

use MOM_time_manager,  only : time_type, real_to_time, time_to_real, time_minus_signed
use MOM_time_manager,  only : set_date, get_date, increment_date
use MOM_time_manager,  only : operator(+), operator(-), operator(<), operator(>), operator(>=)
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_file_parser,   only : param_file_type, get_param
use MOM_io,            only : file_exists, open_ASCII_file, READONLY_FILE, close_file
use MOM_io,            only : MOM_infra_file, vardesc, MOM_field
use MOM_io,            only : var_desc, create_MOM_file, SINGLE_FILE, MOM_write_field
use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE
use MOM_tidal_forcing, only : astro_longitudes, astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency

implicit none ; private

public HA_init, HA_accum

#include <MOM_memory.h>

!> The private control structure for storing the HA info of a particular field
type, private :: HA_type
  character(len=16) :: key = "none"          !< Name of the field of which harmonic analysis is to be performed
  character(len=1)  :: grid                  !< The grid on which the field is defined ('h', 'q', 'u', or 'v')
  real :: old_time = -1.0                    !< The time of the previous accumulating step [T ~> s]
  real, allocatable :: ref(:,:)              !< The initial field in arbitrary units [A]
  real, allocatable :: FtF(:,:)              !< Accumulator of (F' * F) [nondim]
  real, allocatable :: FtSSH(:,:,:)          !< Accumulator of (F' * SSH_in) in arbitrary units [A]
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
end type HA_type

!> A linked list of control structures that store the HA info of different fields
type, private :: HA_node
  type(HA_type)          :: this             !< Control structure of the current field in the list
  type(HA_node), pointer :: next             !< The list of other fields
end type HA_node

!> The public control structure of the MOM_harmonic_analysis module
type, public :: harmonic_analysis_CS ; private
  logical :: HAready = .false.               !< If true, perform harmonic analysis
  type(time_type) :: &
    time_start, &                            !< Start time of harmonic analysis
    time_end, &                              !< End time of harmonic analysis
    time_ref                                 !< Reference time (t = 0) used to calculate tidal forcing
  real, allocatable, dimension(:) :: &
    freq, &                                  !< The frequency of a tidal constituent [T-1 ~> s-1]
    phase0, &                                !< The phase of a tidal constituent at time 0 [rad]
    tide_fn, &                               !< Amplitude modulation of tides by nodal cycle [nondim].
    tide_un                                  !< Phase modulation of tides by nodal cycle [rad].
  integer :: nc                              !< The number of tidal constituents in use
  integer :: length                          !< Number of fields of which harmonic analysis is to be performed
  character(len=4), allocatable, dimension(:) :: const_name !< The name of each constituent
  character(len=255) :: path                 !< Path to directory where output will be written
  type(unit_scale_type)  :: US               !< A dimensional unit scaling type
  type(HA_node), pointer :: list => NULL()   !< A linked list for storing the HA info of different fields
end type harmonic_analysis_CS


  interface
module subroutine HA_init(Time, US, param_file, nc, CS)
  type(time_type),       intent(in)  :: Time        !< The current model time
  type(unit_scale_type), intent(in)  :: US          !< A dimensional unit scaling type
  type(param_file_type), intent(in)  :: param_file  !< A structure to parse for run-time parameters
  integer,               intent(in)  :: nc          !< The number of tidal constituents in use
  type(harmonic_analysis_CS), intent(out) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
                                                    !! equilibrium tide. Set to false if providing tidal phases
                                                    !! that have already been shifted by the
                                                    !! astronomical/equilibrium argument
                                                    !! calculating tidal forcing.
                                                    !! tidal phases at t = 0.


end subroutine HA_init
module subroutine HA_register(key, grid, CS)
  character(len=*),           intent(in)    :: key     !< Name of the current field
  character(len=1),           intent(in)    :: grid    !< The grid on which the key field is defined
  type(harmonic_analysis_CS), intent(inout) :: CS      !< Control structure of the MOM_harmonic_analysis module

  ! Local variables

end subroutine HA_register
module subroutine HA_accum(key, data, Time, G, CS)
  character(len=*),           intent(in) :: key  !< Name of the current field
  real, dimension(:,:),       intent(in) :: data !< Input data of which harmonic analysis is to be performed [A]
  type(time_type),            intent(in) :: Time !< The current model time
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(harmonic_analysis_CS), intent(inout) :: CS   !< Control structure of the MOM_harmonic_analysis module

  ! Local variables

  ! Exit the accumulator in the following cases
end subroutine HA_accum
module subroutine HA_write(ha1, Time, G, CS)
  type(HA_type), pointer,     intent(in) :: ha1    !< Control structure for the current field
  type(time_type),            intent(in) :: Time   !< The current model time
  type(ocean_grid_type),      intent(in) :: G      !< The ocean's grid structure
  type(harmonic_analysis_CS), intent(in) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables


end subroutine HA_write
module subroutine HA_solver(ha1, nc, FtF, x)
  type(HA_type), pointer,              intent(in)  :: ha1    !< Control structure for the current field
  integer,                             intent(in)  :: nc     !< Number of harmonic constituents
  real, dimension(:,:),                intent(in)  :: FtF    !< Accumulator of (F' * F) for all fields [nondim]
  real, dimension(ha1%is:ha1%ie,ha1%js:ha1%je,2*nc+1), &
                                       intent(out) :: x      !< Solution vector of harmonic constants [A]

  ! Local variables

  ! Cholesky decomposition
end subroutine HA_solver
  end interface

end module MOM_harmonic_analysis
