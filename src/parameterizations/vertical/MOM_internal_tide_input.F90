! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates energy input to the internal tides
module MOM_int_tide_input

use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,    only : diag_ctrl, query_averaging_enabled
use MOM_diag_mediator,    only : disable_averaging, enable_averages
use MOM_diag_mediator,    only : safe_alloc_ptr, post_data, register_diag_field
use MOM_debugging,        only : hchksum
use MOM_error_handler,    only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,      only : read_param
use MOM_forcing_type,     only : forcing
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, vardesc, MOM_read_data
use MOM_interface_heights, only : thickness_to_dz, find_rho_bottom
use MOM_isopycnal_slopes, only : vert_fill_TS
use MOM_string_functions, only : extractWord
use MOM_time_manager,     only : time_type, set_time, operator(+), operator(<=)
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : calculate_density_derivs, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public set_int_tide_input, int_tide_input_init, int_tide_input_end
public get_input_TKE, get_barotropic_tidal_vel

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters that regulate internal tide energy inputs.
type, public :: int_tide_input_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: debug      !< If true, write verbose checksums for debugging.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                        !! regulate the timing of diagnostic output.
  real :: TKE_itide_maxi !< Maximum Internal tide conversion
                        !! available to mix above the BBL [H Z2 T-3 ~> m3 s-3 or W m-2]
  real :: kappa_fill    !< Vertical diffusivity used to interpolate sensible values
                        !! of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  real, allocatable, dimension(:,:,:) :: TKE_itidal_coef
            !< The time-invariant field that enters the TKE_itidal input calculation noting that the
            !! stratification and perhaps density are time-varying [R Z4 H-1 T-2 ~> J m-2 or J m kg-1].
  real, allocatable, dimension(:,:,:) :: &
    TKE_itidal_input, & !< The internal tide TKE input at the bottom of the ocean [H Z2 T-3 ~> m3 s-3 or W m-2].
    tideamp             !< The amplitude of the tidal velocities [L T-1 ~> m s-1].

  character(len=200) :: inputdir !< The directory for input files.

  logical :: int_tide_source_test    !< If true, apply an arbitrary generation site
                                     !! for internal tide testing
  type(time_type) :: time_max_source !< A time for use in testing internal tides
  real    :: int_tide_source_x       !< X Location of generation site
                                     !! for internal tide for testing [degrees_E] or [km]
  real    :: int_tide_source_y       !< Y Location of generation site
                                     !! for internal tide for testing [degrees_N] or [km]
  integer :: int_tide_source_i       !< I Location of generation site
  integer :: int_tide_source_j       !< J Location of generation site
  logical :: int_tide_use_glob_ij    !< Use global indices for generation site
  integer :: nFreq = 0               !< The number of internal tide frequency bands


  !>@{ Diagnostic IDs
  integer, allocatable, dimension(:) :: id_TKE_itidal_itide
  integer :: id_Nb = -1, id_N2_bot = -1
  !>@}
end type int_tide_input_CS

!> This type is used to exchange fields related to the internal tides.
type, public :: int_tide_input_type
  real, allocatable, dimension(:,:) :: &
    h2, &               !< The squared topographic roughness height [Z2 ~> m2].
    Nb, &               !< The bottom stratification [T-1 ~> s-1].
    Rho_bot             !< The bottom density or the Boussinesq reference density [R ~> kg m-3].
end type int_tide_input_type


  interface
module subroutine set_int_tide_input(u, v, h, tv, fluxes, itide, dt, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv !< A structure containing pointers to the
                                                                  !! thermodynamic fields
  type(forcing),                              intent(in)    :: fluxes !< A structure of thermodynamic surface fluxes
  type(int_tide_input_type),                  intent(inout) :: itide !< A structure containing fields related
                                                                  !! to the internal tide sources.
  real,                                       intent(in)    :: dt !< The time increment [T ~> s].
  type(int_tide_input_CS),                    pointer       :: CS !< This module's control structure.

  ! Local variables

                  ! the massless layers filled vertically by diffusion.
                        ! equation of state.
                          ! to mks [T3 kg H-1 Z-2 s-3 ~> kg m-3 or 1]
                          ! units [H Z2 s3 T-3 kg-1 ~> m3 kg-1 or 1]


end subroutine set_int_tide_input
module subroutine find_N2_bottom(G, GV, US, tv, fluxes, h, T_f, S_f, h2, N2_bot, Rho_bot, h_bot, k_bot)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),                     intent(in)  :: tv   !< A structure containing pointers to the
                                                                 !! thermodynamic fields
  type(forcing),                             intent(in)  :: fluxes !< A structure of thermodynamic surface fluxes
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: T_f  !< Temperature after vertical filtering to
                                                                 !! smooth out the values in thin layers [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: S_f  !< Salinity after vertical filtering to
                                                                 !! smooth out the values in thin layers [S ~> ppt].
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: h2   !< Bottom topographic roughness [Z2 ~> m2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: N2_bot !< The squared buoyancy frequency at the
                                                                 !! ocean bottom [T-2 ~> s-2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: Rho_bot !< The average density near the ocean
                                                                 !! bottom [R ~> kg m-3]
  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: h_bot !< Bottom boundary layer thickness [H ~> m or kg m-2]
  integer, dimension(SZI_(G),SZJ_(G)),       intent(out) :: k_bot !< Bottom boundary layer top layer index

  ! Local variables

                  ! density [H T-2 R-1 ~> m4 s-2 kg-1 or m s-2].

end subroutine find_N2_bottom
module subroutine get_input_TKE(G, TKE_itidal_input, nFreq, CS)
  type(ocean_grid_type), intent(in)    :: G !< The ocean's grid structure (in).
  integer, intent(in) :: nFreq !< number of frequencies
  real, dimension(SZI_(G),SZJ_(G),nFreq), &
                         intent(out) :: TKE_itidal_input !< The energy input to the internal waves
                                                         !! [H Z2 T-3 ~> m3 s-3 or W m-2].
  type(int_tide_input_CS),   target       :: CS !< A pointer that is set to point to the control
                                                 !! structure for the internal tide input module.

end subroutine get_input_TKE
module subroutine get_barotropic_tidal_vel(G, vel_btTide, nFreq, CS)
  type(ocean_grid_type), intent(in)    :: G !< The ocean's grid structure (in).
  integer, intent(in) :: nFreq !< number of frequencies
  real, dimension(SZI_(G),SZJ_(G),nFreq), &
                         intent(out) :: vel_btTide !< Barotropic velocity read from file [L T-1 ~> m s-1].
  type(int_tide_input_CS),   target       :: CS !< A pointer that is set to point to the control
                                                 !! structure for the internal tide input module.

end subroutine get_barotropic_tidal_vel
module subroutine int_tide_input_init(Time, G, GV, US, param_file, diag, CS, itide)
  type(time_type),           intent(in)    :: Time !< The current model time
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),   target, intent(inout) :: diag !< structure used to regulate diagnostic output.
  type(int_tide_input_CS),   pointer       :: CS   !< This module's control structure, which is initialized here.
  type(int_tide_input_type), pointer       :: itide !< A structure containing fields related
                                                   !! to the internal tide sources.
  ! Local variables
  ! This include declares and sets the variable "version".

                             ! to the mean depth [nondim]
                             ! tidal amplitude file is not present.
                             ! to mks [T3 kg H-1 Z-2 s-3 ~> kg m-3 or 1]
                             ! units [H Z2 s3 T-3 kg-1 ~> m3 kg-1 or 1]
                             !! for testing internal tides (BDM)

end subroutine int_tide_input_init
module subroutine int_tide_input_end(CS)
  type(int_tide_input_CS), pointer :: CS !< This module's control structure, which is deallocated here.

end subroutine int_tide_input_end
  end interface

end module MOM_int_tide_input
