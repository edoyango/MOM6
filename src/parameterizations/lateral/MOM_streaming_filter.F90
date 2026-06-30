! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Streaming band-pass filter for detecting the instantaneous tidal signals in the simulation

module MOM_streaming_filter

use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE, FATAL
use MOM_file_parser,   only : get_param, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : axis_info, set_axis_info
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_tidal_forcing, only : tidal_frequency
use MOM_time_manager,  only : time_type, time_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public Filt_register, Filt_init, Filt_accum

#include <MOM_memory.h>

!> Control structure for the MOM_streaming_filter module
type, public :: Filter_CS ; private
  integer :: nf                        !< Number of filters to be used in the simulation
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
  character(len=8) :: key              !< Identifier of the variable to be filtered
  character(len=2), allocatable, dimension(:) :: filter_names !< Names of filters
  real, allocatable, dimension(:)      :: filter_omega !< Target frequencies of filters [rad T-1 ~> rad s-1]
  real, allocatable, dimension(:)      :: filter_alpha !< Bandwidth parameters of filters [nondim]
  real, allocatable, dimension(:,:,:)  :: s1, &        !< A dummy variable for solving the system of ODEs [A]
                                          u1           !< Filtered data, representing the narrow-band signal
                                                       !< oscillating around the target frequency [A]
  real :: old_time = -1.0              !< The time of the previous accumulating step [T ~> s]
end type Filter_CS


  interface
module subroutine Filt_register(nf, key, grid, HI, CS, restart_CS)
  integer,               intent(in)    :: nf           !< Number of filters to be used in the simulation
  character(len=*),      intent(in)    :: key          !< Identifier of the variable to be filtered
  character(len=*),      intent(in)    :: grid         !< Horizontal grid location: "h", "u", or "v"
  type(hor_index_type),  intent(in)    :: HI           !< Horizontal index type structure
  type(Filter_CS),       intent(out)   :: CS           !< Control structure of MOM_streaming_filter
  type(MOM_restart_CS),  intent(inout) :: restart_CS   !< MOM restart control structure

  ! Local variables

end subroutine Filt_register
module subroutine Filt_init(param_file, US, CS, restart_CS)
  type(param_file_type), intent(in)    :: param_file   !< A structure to parse for run-time parameters
  type(unit_scale_type), intent(in)    :: US           !< A dimensional unit scaling type
  type(Filter_CS),       intent(inout) :: CS           !< Control structure of MOM_streaming_filter
  type(MOM_restart_CS),  intent(in)    :: restart_CS   !< MOM restart control structure

  ! Local variables

end subroutine Filt_init
module subroutine Filt_accum(u, u1, Time, US, CS)
  real, dimension(:,:,:), pointer, intent(out)   :: u1   !< Output of the filter [A]
  type(time_type),                 intent(in)    :: Time !< The current model time
  type(unit_scale_type),           intent(in)    :: US   !< A dimensional unit scaling type
  type(Filter_CS),        target,  intent(inout) :: CS   !< Control structure of MOM_streaming_filter
  real, dimension(CS%is:CS%ie,CS%js:CS%je), intent(in) :: u !< Input into the filter [A]

  ! Local variables

end subroutine Filt_accum
  end interface

end module MOM_streaming_filter
