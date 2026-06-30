! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Reports integrated quantities for monitoring the model state
module MOM_sum_output

use iso_fortran_env, only : int64
use MOM_checksums,     only : is_NaN, field_checksum
use MOM_coms,          only : sum_across_PEs, PE_here, root_PE, num_PEs, max_across_PEs
use MOM_coms,          only : reproducing_sum, reproducing_sum_EFP, EFP_to_real, real_to_EFP
use MOM_coms,          only : EFP_type, operator(+), operator(-), assignment(=), EFP_sum_across_PEs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_io,            only : create_MOM_file, reopen_MOM_file
use MOM_io,            only : MOM_infra_file, MOM_netcdf_file, MOM_field
use MOM_io,            only : file_exists, slasher, vardesc, var_desc, MOM_write_field
use MOM_io,            only : field_size, read_variable, read_attribute, open_ASCII_file, stdout
use MOM_io,            only : axis_info, set_axis_info, delete_axis_info, get_filename_appendix
use MOM_io,            only : attribute_info, set_attribute_info, delete_attribute_info
use MOM_io,            only : APPEND_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_spatial_means, only : array_global_min_max
use MOM_time_manager,  only : time_type, get_time, get_date, set_time
use MOM_time_manager,  only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,  only : operator(/=), operator(<=), operator(>=), operator(<), operator(>)
use MOM_time_manager,  only : get_calendar_type, time_type_to_real, NO_CALENDAR
use MOM_tracer_flow_control, only : tracer_flow_control_CS, call_tracer_stocks
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public write_energy, accumulate_net_input
public MOM_sum_output_init, MOM_sum_output_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

integer, parameter :: NUM_FIELDS = 17 !< Number of diagnostic fields
character (*), parameter :: depth_chksum_attr = "bathyT_checksum"
                                      !< Checksum attribute name of G%bathyT
                                      !! over the compute domain
character (*), parameter :: area_chksum_attr = "mask2dT_areaT_checksum"
                                      !< Checksum attribute of name of
                                      !! G%mask2dT * G%areaT over the compute
                                      !! domain

!> A list of depths and corresponding globally integrated ocean area at each
!! depth and the ocean volume below each depth.
type :: Depth_List
  integer                         :: listsize  !< length of the list <= niglobal*njglobal + 1
  real, allocatable, dimension(:) :: depth     !< A list of depths [Z ~> m]
  real, allocatable, dimension(:) :: area      !< The cross-sectional area of the ocean at that depth [L2 ~> m2]
  real, allocatable, dimension(:) :: vol_below !< The ocean volume below that depth [Z L2 ~> m3]
end type Depth_List

!> The control structure for the MOM_sum_output module
type, public :: sum_output_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.

  type(Depth_List)              :: DL !< The sorted depth list.

  integer, allocatable, dimension(:) :: lH
                                !< This saves the entry in DL with a volume just
                                !! less than the volume of fluid below the interface.
  logical :: do_APE_calc        !<   If true, calculate the available potential energy of the
                                !! interfaces.  Disabling this reduces the memory footprint of
                                !! high-PE-count models dramatically.
  logical :: read_depth_list    !<   Read the depth list from a file if it exists
                                !! and write it if it doesn't.
  character(len=200) :: depth_list_file  !< The name of the depth list file.
  real    :: D_list_min_inc     !<  The minimum increment [Z ~> m], between the depths of the
                                !! entries in the depth-list file, 0 by default.
  logical :: require_depth_list_chksum
                                !< Require matching checksums in Depth_list.nc when reading
                                !! the file.
  logical :: update_depth_list_chksum
                                !< Automatically update the Depth_list.nc file if the
                                !! checksums are missing or do not match current values.
  logical :: use_temperature    !<   If true, temperature and salinity are state variables.
  type(EFP_type) :: fresh_water_in_EFP !< The total mass of fresh water added by surface fluxes on
                                  !! this PE since the last time that write_energy was called [kg].
  type(EFP_type) :: net_salt_in_EFP !< The total salt added by surface fluxes on this PE since
                                  !! the last time that write_energy was called [ppt kg].
  type(EFP_type) :: net_heat_in_EFP !<  The total heat added by surface fluxes on this PE since
                                  !! the last time that write_energy was called [J].
  type(EFP_type) :: heat_prev_EFP !<  The total amount of heat in the ocean the last
                                  !! time that write_energy was called [J].
  type(EFP_type) :: salt_prev_EFP !< The total amount of salt in the ocean the last
                                  !! time that write_energy was called [ppt kg].
  type(EFP_type) :: mass_prev_EFP !< The total ocean mass the last time that
                                  !! write_energy was called [kg].
  real    :: dt_in_T            !< The baroclinic dynamics time step [T ~> s].

  type(time_type) :: energysavedays            !< The interval between writing the energies
                                               !! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  !< The starting interval for computing a geometric
                                               !! progression of time deltas between calls to
                                               !! write_energy. This interval will increase by a factor of 2.
                                               !! after each call to write_energy.
  logical         :: energysave_geometric      !< Logical to control whether calls to write_energy should
                                               !! follow a geometric progression
  type(time_type) :: write_energy_time         !< The next time to write to the energy file.
  type(time_type) :: geometric_end_time        !< Time at which to stop the geometric progression
                                               !! of calls to write_energy and revert to the standard
                                               !! energysavedays interval

  real    :: timeunit           !< The length of the units for the time axis and certain input parameters
                                !! including ENERGYSAVEDAYS [s].

  logical :: date_stamped_output !< If true, use dates (not times) in messages to stdout.
  logical :: ISO_date_stamped_output !< If true, use ISO formatted dates in messages to stdout.
  type(time_type) :: Start_time !< The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  integer, pointer :: ntrunc => NULL() !< The number of times the velocity has been
                                !! truncated since the last call to write_energy.
  real    :: max_Energy         !< The maximum permitted energy per unit mass.  If there is
                                !! more energy than this, the model should stop [L2 T-2 ~> m2 s-2].
  integer :: maxtrunc           !< The number of truncations per energy save
                                !! interval at which the run is stopped.
  logical :: write_stocks       !< If true, write the integrated tracer amounts
                                !! to stdout when the energy files are written.
  logical :: write_min_max      !< If true, write the maximum and minimum values of temperature,
                                !! salinity and some tracer concentrations to stdout when the energy
                                !! files are written.
  logical :: write_min_max_loc  !< If true, write the locations of the maximum and minimum values
                                !! of temperature, salinity and some tracer concentrations to stdout
                                !! when the energy files are written.
  integer :: previous_calls = 0 !< The number of times write_energy has been called.
  integer :: prev_n = 0         !< The value of n from the last call.
  type(MOM_netcdf_file) :: fileenergy_nc !< The file handle for the netCDF version of the energy file.
  integer :: fileenergy_ascii   !< The unit number of the ascii version of the energy file.
  type(MOM_field), dimension(NUM_FIELDS+MAX_FIELDS_) :: &
             fields             !< fieldtype variables for the output fields.
  character(len=200) :: energyfile  !< The name of the energy file with path.
end type sum_output_CS


  interface
module subroutine MOM_sum_output_init(G, GV, US, param_file, directory, ntrnc, &
                               Input_start_time, CS)
  type(ocean_grid_type),   intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                       !! parameters.
  character(len=*),        intent(in)    :: directory  !< The directory where the energy file goes.
  integer, target,         intent(inout) :: ntrnc      !< The integer that stores the number of times
                                                       !! the velocity has been truncated since the
                                                       !! last call to write_energy.
  type(time_type),         intent(in)    :: Input_start_time !< The start time of the simulation.
  type(Sum_output_CS),     pointer       :: CS         !< A pointer that is set to point to the
                                                       !! control structure for this module.
  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine MOM_sum_output_init
module subroutine MOM_sum_output_end(CS)
  type(Sum_output_CS), pointer :: CS  !< The control structure returned by a
                                      !! previous call to MOM_sum_output_init.
end subroutine MOM_sum_output_end
module subroutine write_energy(u, v, h, tv, day, n, G, GV, US, CS, tracer_CSp, dt_forcing)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h   !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv  !< A structure pointing to various
                                                !! thermodynamic variables.
  type(time_type),         intent(in)    :: day !< The current model time.
  integer,                 intent(in)    :: n   !< The time step number of the
                                                !! current execution.
  type(Sum_output_CS),     pointer       :: CS  !< The control structure returned by a
                                                !! previous call to MOM_sum_output_init.
  type(tracer_flow_control_CS), pointer  :: tracer_CSp !< Control structure with the tree of
                                                !! all registered tracer packages
  type(time_type),  optional, intent(in) :: dt_forcing !< The forcing time step

  ! Local variables
                       ! volume as is below an interface [Z ~> m].
                       ! the total mass of the ocean [L2 T-2 ~> m2 s-2].
                       ! the last call to this subroutine [R Z L2 ~> kg]
                       ! by the surface fluxes [R Z L2 ~> kg]
                       ! to this subroutine [1e-3 R Z L2 ~> g Salt]
                       ! the surface fluxes [1e-3 R Z L2 ~> g Salt]
                       ! the surface fluxes divided by total mass [ppt].
                       ! by the surface fluxes, divided by the total heat
                       ! capacity of the ocean [C ~> degC]
                       ! height of the basin depth over H otherwise [Z ~> m].
                       ! This makes PE only include real fluid.
end subroutine write_energy
module subroutine accumulate_net_input(fluxes, sfc_state, tv, dt, G, US, CS)
  type(forcing),         intent(in) :: fluxes !< A structure containing pointers to any possible
                                              !! forcing fields.  Unused fields are unallocated.
  type(surface),         intent(in) :: sfc_state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  type(thermo_var_ptrs), intent(in) :: tv     !< A structure pointing to various
                                              !! thermodynamic variables.
  real,                  intent(in) :: dt     !< The amount of time over which to average [T ~> s].
  type(ocean_grid_type), intent(in) :: G      !< The ocean's grid structure.
  type(unit_scale_type), intent(in) :: US     !< A dimensional unit scaling type
  type(Sum_output_CS),   pointer    :: CS     !< The control structure returned by a previous call
                                              !! to MOM_sum_output_init.
  ! Local variables
end subroutine accumulate_net_input
module subroutine depth_list_setup(G, GV, US, DL, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(Depth_List),        intent(inout) :: DL  !< The list of depths, areas and volumes to set up
  type(Sum_output_CS),     pointer       :: CS  !< The control structure returned by a
                                                !! previous call to MOM_sum_output_init.
  ! Local variables

end subroutine depth_list_setup
module subroutine create_depth_list(G, DL, min_depth_inc)
  type(ocean_grid_type), intent(in)    :: G  !< The ocean's grid structure.
  type(Depth_List),      intent(inout) :: DL !< The list of depths, areas and volumes to create
  real,                  intent(in)    :: min_depth_inc !< The minimum increment between depths in the list [Z ~> m]

  ! Local variables


end subroutine create_depth_list
module subroutine write_depth_list(G, US, DL, filename)
  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure.
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  type(Depth_List),      intent(in) :: DL  !< The list of depths, areas and volumes to write
  character(len=*),      intent(in) :: filename !< The path to the depth list file to write.

  ! Local variables

  ! All ranks are required to compute the global checksum
end subroutine write_depth_list
module subroutine read_depth_list(G, US, DL, filename, require_chksum, file_matches)
  type(ocean_grid_type), intent(in)    :: G   !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(Depth_List),      intent(inout) :: DL  !< The list of depths, areas and volumes
  character(len=*),      intent(in)    :: filename !< The path to the depth list file to read.
  logical,               intent(in)    :: require_chksum !< If true, missing or mismatched depth
                                              !! and area checksums result in a fatal error.
  logical, optional,     intent(out)   :: file_matches !< If present, this indicates whether the file
                                              !! has been read with matching depth and area checksums

  ! Local variables

  ! Check bathymetric consistency between this configuration and the depth list file.
end subroutine read_depth_list
module subroutine get_depth_list_checksums(G, US, depth_chksum, area_chksum)
  type(ocean_grid_type), intent(in) :: G          !< Ocean grid structure
  type(unit_scale_type), intent(in) :: US         !< A dimensional unit scaling type
  character(len=16), intent(out) :: depth_chksum  !< Depth checksum hexstring
  character(len=16), intent(out) :: area_chksum   !< Area checksum hexstring

  ! Local variables

end subroutine get_depth_list_checksums
  end interface

end module MOM_sum_output
