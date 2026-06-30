! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Tidal contributions to geopotential
module MOM_tidal_forcing

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : field_exists, file_exists, MOM_read_data
use MOM_time_manager,  only : set_date, time_type, time_minus_signed
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public calc_tidal_forcing, tidal_forcing_init, tidal_forcing_end
public calc_tidal_forcing_legacy
! MOM_open_boundary uses the following to set tides on the boundary.
public astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency

#include <MOM_memory.h>

integer, parameter :: MAX_CONSTITUENTS = 10 !< The maximum number of tidal
                                            !! constituents that could be used.
!> Simple type to store astronomical longitudes used to calculate tidal phases.
type, public :: astro_longitudes
  real :: s  !< Mean longitude of moon [rad]
  real :: h  !< Mean longitude of sun [rad]
  real :: p  !< Mean longitude of lunar perigee [rad]
  real :: N  !< Longitude of ascending node [rad]
end type astro_longitudes

!> The control structure for the MOM_tidal_forcing module
type, public :: tidal_forcing_CS ; private
  logical :: use_tidal_sal_file !< If true, Read the tidal self-attraction
                      !! and loading from input files, specified
                      !! by TIDAL_INPUT_FILE.
  logical :: use_tidal_sal_prev !< If true, use the SAL from the previous
                      !! iteration of the tides to facilitate convergence.
  logical :: use_eq_phase !< If true, tidal forcing is phase-shifted to match
                      !! equilibrium tide. Set to false if providing tidal phases
                      !! that have already been shifted by the
                      !! astronomical/equilibrium argument.
  real    :: sal_scalar = 0.0 !< The constant of proportionality between self-attraction and
                      !! loading (SAL) geopotential anomaly and total geopotential geopotential
                      !! anomalies. This is only used if USE_PREVIOUS_TIDES is true. [nondim].
  integer :: nc       !< The number of tidal constituents in use.
  real, dimension(MAX_CONSTITUENTS) :: &
    freq, &           !< The frequency of a tidal constituent [rad T-1 ~> rad s-1].
    phase0, &         !< The phase of a tidal constituent at time 0 [rad].
    amp, &            !< The amplitude of a tidal constituent at time 0 [Z ~> m].
    love_no           !< The Love number of a tidal constituent at time 0 [nondim].
  integer :: struct(MAX_CONSTITUENTS) !< An encoded spatial structure for each constituent
  character (len=16) :: const_name(MAX_CONSTITUENTS) !< The name of each constituent

  type(time_type) :: time_ref !< Reference time (t = 0) used to calculate tidal forcing.
  type(astro_longitudes) :: tidal_longitudes !< Astronomical longitudes used to calculate
                                   !! tidal phases at t = 0.
  real, allocatable :: &
    sin_struct(:,:,:), &    !< The sine based structures that can be associated with
                            !! the astronomical forcing [nondim].
    cos_struct(:,:,:), &    !< The cosine based structures that can be associated with
                            !! the astronomical forcing [nondim].
    cosphasesal(:,:,:), &   !< The cosine of the phase of the self-attraction and loading amphidromes [nondim].
    sinphasesal(:,:,:), &   !< The sine of the phase of the self-attraction and loading amphidromes [nondim].
    ampsal(:,:,:), &        !< The amplitude of the SAL [Z ~> m].
    cosphase_prev(:,:,:), & !< The cosine of the phase of the amphidromes in the previous tidal solutions [nondim].
    sinphase_prev(:,:,:), & !< The sine of the phase of the amphidromes in the previous tidal solutions [nondim].
    amp_prev(:,:,:), &      !< The amplitude of the previous tidal solution [Z ~> m].
    tide_fn(:), &           !< Amplitude modulation of tides by nodal cycle [nondim].
    tide_un(:)              !< Phase modulation of tides by nodal cycle [rad].
end type tidal_forcing_CS

integer :: id_clock_tides !< CPU clock for tides


  interface
module subroutine astro_longitudes_init(time_ref, longitudes)
  type(time_type), intent(in) :: time_ref            !> Time to calculate longitudes for.
  type(astro_longitudes), intent(out) :: longitudes  !> Lunar and solar longitudes at time_ref.

  ! Local variables

  ! Find date at time_ref in days since midnight at the start of 1900-01-01
end subroutine astro_longitudes_init
module function eq_phase(constit, longitudes)
  character (len=2), intent(in) :: constit !> Name of constituent (e.g., M2).
  type(astro_longitudes), intent(in) :: longitudes   !> Mean longitudes calculated using astro_longitudes_init
  real :: eq_phase                         !> The equilibrium phase argument for the constituent [rad].

end function eq_phase
module function tidal_frequency(constit)
  character (len=2), intent(in) :: constit !> Constituent to look up
  real :: tidal_frequency                  !> Angular frequency [rad s-1]

end function tidal_frequency
module subroutine nodal_fu(constit, nodelon, fn, un)
  character (len=2), intent(in)  :: constit !> Tidal constituent to find modulation for.
  real,              intent(in)  :: nodelon !> Longitude of ascending node [rad], which
                                            !! can be calculated using astro_longitudes_init.
  real,              intent(out) :: fn      !> Amplitude modulation [nondim]
  real,              intent(out) :: un      !> Phase modulation [rad]


end subroutine nodal_fu
module subroutine tidal_forcing_init(Time, G, US, param_file, CS)
  type(time_type),        intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(tidal_forcing_CS), intent(inout) :: CS   !< Tidal forcing control structure

  ! Local variables
                                              !! calculating tidal forcing.
  ! This include declares and sets the variable "version".

end subroutine tidal_forcing_init
module subroutine find_in_files(filenames, varname, array, G, scale)
  character(len=*), dimension(:),   intent(in)  :: filenames !< The names of the files to search for the named variable
  character(len=*),                 intent(in)  :: varname   !< The name of the variable to read
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: array     !< The array to fill with the data [arbitrary]
  real,                   optional, intent(in)  :: scale     !< A factor by which to rescale the array to translate it
                                                             !! into its desired units [arbitrary]
  ! Local variables

end subroutine find_in_files
module subroutine calc_tidal_forcing(Time, e_tide_eq, e_tide_sal, G, US, CS)
  type(ocean_grid_type),            intent(in)  :: G          !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time       !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_eq  !< The geopotential height anomalies
                                                              !! due to the equilibrium tides [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_sal !< The geopotential height anomalies
                                                              !! due to the tidal SAL [Z ~> m].
  type(unit_scale_type),            intent(in)  :: US         !< A dimensional unit scaling type
  type(tidal_forcing_CS),           intent(in)  :: CS         !< The control structure returned by a
                                                              !! previous call to tidal_forcing_init.

  ! Local variables
end subroutine calc_tidal_forcing
module subroutine calc_tidal_forcing_legacy(Time, e_sal, e_sal_tide, e_tide_eq, e_tide_sal, G, US, CS)
  type(ocean_grid_type),            intent(in)  :: G          !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time       !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: e_sal      !< The self-attraction and loading fields
                                                              !! calculated previously used to
                                                              !! initialized e_sal_tide [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_sal_tide !< The total geopotential height anomalies
                                                              !! due to both SAL and tidal forcings [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_eq  !< The geopotential height anomalies
                                                              !! due to the equilibrium tides [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_sal !< The geopotential height anomalies
                                                              !! due to the tidal SAL [Z ~> m].
  type(unit_scale_type),            intent(in)  :: US         !< A dimensional unit scaling type
  type(tidal_forcing_CS),           intent(in)  :: CS         !< The control structure returned by a
                                                              !! previous call to tidal_forcing_init.

  ! Local variables
end subroutine calc_tidal_forcing_legacy
module subroutine tidal_forcing_end(CS)
  type(tidal_forcing_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                              !! to tidal_forcing_init; it is deallocated here.

end subroutine tidal_forcing_end
  end interface

end module MOM_tidal_forcing
