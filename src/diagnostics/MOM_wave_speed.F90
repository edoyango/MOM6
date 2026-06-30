! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines for calculating baroclinic wave speeds
module MOM_wave_speed

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : log_version
use MOM_grid, only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, interpolate_column
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density_derivs, calculate_specific_vol_derivs

implicit none ; private

#include <MOM_memory.h>

public wave_speed, wave_speeds, wave_speed_init, wave_speed_set_param

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for MOM_wave_speed
type, public :: wave_speed_CS ; private
  logical :: initialized = .false.     !< True if this control structure has been initialized.
  logical :: use_ebt_mode = .false.    !< If true, calculate the equivalent barotropic wave speed instead
                                       !! of the first baroclinic wave speed.
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  logical :: better_cg1_est = .false.  !< If true, use an improved estimate of the first mode
                                       !! internal wave speed.
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent barotropic
                                       !! wave speed [nondim]. This parameter controls the default behavior of
                                       !! wave_speed() which can be overridden by optional arguments.
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [H ~> m or kg m-2].
                                       !! If this parameter is negative, this limiting does not occur.
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  real :: min_speed2 = 0.              !< The minimum mode 1 internal wave speed squared [L2 T-2 ~> m2 s-2]
  real :: wave_speed_tol = 0.001       !< The fractional tolerance with which to solve for the wave
                                       !! speeds [nondim]
  real :: c1_thresh = -1.0             !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but
                                       !! are simply reported as 0 [L T-1 ~> m s-1].  A non-negative
                                       !! value must be specified via a call to wave_speed_init for
                                       !! the subroutine wave_speeds to be used (but not wave_speed).
  type(remapping_CS) :: remap_2018_CS  !< Used for vertical remapping when calculating equivalent barotropic
                                       !! mode structure for answer dates below 20190101.
  type(remapping_CS) :: remap_CS       !< Used for vertical remapping when calculating equivalent barotropic
                                       !! mode structure.
  integer :: remap_answer_date = 99991231 !< The vintage of the order of arithmetic and expressions to use
                                       !! for remapping.  Values below 20190101 recover the remapping
                                       !! answers from 2018, while higher values use more robust
                                       !! forms of the same remapping expressions.
  type(diag_ctrl), pointer :: diag     !< Diagnostics control structure
end type wave_speed_CS


  interface
module subroutine wave_speed(h, tv, G, GV, US, cg1, CS, halo_size, use_ebt_mode, mono_N2_column_fraction, &
                      mono_N2_depth, modal_structure)
  type(ocean_grid_type),            intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),          intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)  :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: cg1 !< First mode internal wave speed [L T-1 ~> m s-1]
  type(wave_speed_CS),              intent(in)  :: CS !< Wave speed control struct
  integer,                optional, intent(in)  :: halo_size !< Width of halo within which to
                                                       !! calculate wave speeds
  logical,                optional, intent(in)  :: use_ebt_mode !< If true, use the equivalent
                                          !! barotropic mode instead of the first baroclinic mode.
  real,                   optional, intent(in)  :: mono_N2_column_fraction !< The lower fraction
                                          !! of water column over which N2 is limited as monotonic
                                          !! for the purposes of calculating vertical modal structure [nondim].
  real,                   optional, intent(in)  :: mono_N2_depth !< A depth below which N2 is limited as
                                          !! monotonic for the purposes of calculating vertical
                                          !! modal structure [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          optional, intent(out) :: modal_structure !< Normalized model structure [nondim]

  ! Local variables
                  ! the thickness of the layer below (Igl) or above (Igu) it, in [T2 L-2 ~> s2 m-2].
                  ! units of the eigenvalue change with the number of layers and because of the
                  ! dynamic rescaling that is used to keep det in a numerically representable range,
                  ! the units of of det are hard to interpret, but det/ddet is always in units
                  ! of [T2 L-2 ~> s2 m-2]
                   ! thicknesses [H R-1 ~> m4 kg-1 or m], negative for stable stratification.
                   ! its derivative with lam between rows of the Thomas algorithm solver [L2 s2 T-2 m-2 ~> nondim].
                   ! The exact value should not matter for the final result if it is an even power of 2.
                    ! the total water column can be merged for efficiency [nondim].
                    ! when deciding to merge layers in the calculation [nondim]
                    ! with each iteration.  Because of all of the dynamic rescaling of the determinant
                    ! between rows, its units are not easily interpretable, but the ratio of det/ddet
                    ! always has units of [T2 L-2 ~> s2 m-2]
                         ! in units of [L2 T-2 ~> m2 s-2] after it is modified inside of tdma6.

end subroutine wave_speed
module subroutine tdma6(n, a, c, lam, y)
  integer,            intent(in)    :: n !< Number of rows of matrix
  real, dimension(:), intent(in)    :: a !< Lower diagonal   [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(in)    :: c !< Upper diagonal   [T2 L-2 ~> s2 m-2]
  real,               intent(in)    :: lam !< Scalar subtracted from leading diagonal [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(inout) :: y !< RHS on entry [A ~> a], result on exit [A L2 T-2 ~> a m2 s-2]

  ! Local variables

end subroutine tdma6
module subroutine wave_speeds(h, tv, G, GV, US, nmodes, cn, CS, w_struct, u_struct, u_struct_max, u_struct_bot, Nb, int_w2, &
                       int_U2, int_N2w2, halo_size)
  type(ocean_grid_type),                           intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                         intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                           intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),       intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                           intent(in)  :: tv !< Thermodynamic variables
  integer,                                         intent(in)  :: nmodes !< Number of modes
  type(wave_speed_CS),                             intent(in)  :: CS !< Wave speed control struct
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1,nmodes),intent(out) :: w_struct !< Wave vertical velocity profile [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),nmodes),intent(out) :: u_struct !< Wave horizontal velocity profile
                                                                     !! [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: cn !< Waves speeds [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: u_struct_max !< Maximum of wave horizontal velocity
                                                                     !! profile [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: u_struct_bot !< Bottom value of wave horizontal
                                                                     !! velocity profile [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G)),                intent(out) :: Nb !< Bottom value of buoyancy freqency
                                                                     !! [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_w2 !< depth-integrated vertical velocity
                                                                     !! profile squared [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_U2 !< depth-integrated horizontal velocity
                                                                     !! profile squared [H Z-2 ~> m-1 or kg m-4]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_N2w2 !< depth-integrated buoyancy frequency
                                                                     !! times vertical velocity profile
                                                                     !! squared [H T-2 ~> m s-2 or kg m-2 s-2]
  integer,                               optional, intent(in)  :: halo_size !< Width of halo within which to
                                                                     !! calculate wave speeds

  ! Local variables
end subroutine wave_speeds
module subroutine tridiag_det(a, c, ks, ke, lam, det, ddet, row_scale)
  real, dimension(:), intent(in) :: a     !< Lower diagonal of matrix (first entry unused) [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(in) :: c     !< Upper diagonal of matrix (last entry unused) [T2 L-2 ~> s2 m-2]
  integer,            intent(in) :: ks    !< Starting index to use in determinant
  integer,            intent(in) :: ke    !< Ending index to use in determinant
  real,               intent(in) :: lam   !< Value subtracted from b [T2 L-2 ~> s2 m-2]
  real,               intent(out):: det   !< Determinant of the matrix in dynamically rescaled units that
                                          !! depend on the number of rows and the cumulative magnitude of
                                          !! det and are therefore difficult to interpret, but the units
                                          !! of det/ddet are always in [T2 L-2 ~> s2 m-2]
  real,               intent(out):: ddet  !< Derivative of determinant with lam in units that are dynamically
                                          !! rescaled along with those of det, such that the units of
                                          !! det/ddet are always in [T2 L-2 ~> s2 m-2]
  real,               intent(in) :: row_scale !< A scaling factor of the rows of the matrix to
                                          !! limit the growth of the determinant [L2 s2 T-2 m-2 ~> 1]
  ! Local variables
                           ! that vary with the number of layers that have been worked on [various]
                           ! layers [various], but the units of detKm1/ddetKm1 are [T2 L-2 ~> s2 m-2]

end subroutine tridiag_det
module subroutine wave_speed_init(CS, GV, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, remap_answers_2018, &
                           remap_answer_date, better_speed_est, om4_remap_via_sub_cells, &
                           min_speed, wave_speed_tol, c1_thresh)
  type(wave_speed_CS), intent(inout) :: CS  !< Wave speed control struct
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                     !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                     !! which N2 is limited as monotonic for the purposes of
                                     !! calculating the vertical modal structure [nondim].
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                     !! as monotonic for the purposes of calculating the
                                     !! vertical modal structure [H ~> m or kg m-2].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                     !! that recover the remapping answers from 2018.  Otherwise
                                     !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date  !< The vintage of the order of arithmetic and expressions
                                      !! to use for remapping.  Values below 20190101 recover the remapping
                                      !! answers from 2018, while higher values use more robust
                                      !! forms of the same remapping expressions.
  logical, optional, intent(in) :: better_speed_est !< If true, use a more robust estimate of the first
                                     !! mode speed as the starting point for iterations.
  logical, optional, intent(in) :: om4_remap_via_sub_cells !< Use the OM4-era ramap_via_sub_cells
                                     !! for calculating the EBT structure
  real,    optional, intent(in) :: min_speed !< If present, set a floor in the first mode speed
                                     !! below which 0 is returned [L T-1 ~> m s-1].
  real,    optional, intent(in) :: wave_speed_tol !< The fractional tolerance for finding the
                                     !! wave speeds [nondim]
  real,    optional, intent(in) :: c1_thresh !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but are
                                       !! simply reported as 0 [L T-1 ~> m s-1].  A non-negative value
                                       !! must be specified for wave_speeds to be used (but not wave_speed).

  ! This include declares and sets the variable "version".

end subroutine wave_speed_init
module subroutine wave_speed_set_param(CS, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, remap_answers_2018, &
                                remap_answer_date, better_speed_est, min_speed, wave_speed_tol, c1_thresh)
  type(wave_speed_CS), intent(inout)  :: CS
                                      !< Control structure for MOM_wave_speed
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                      !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                      !! which N2 is limited as monotonic for the purposes of
                                      !! calculating the vertical modal structure [nondim].
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                      !! as monotonic for the purposes of calculating the
                                      !! vertical modal structure [H ~> m or kg m-2].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                      !! that recover the remapping answers from 2018.  Otherwise
                                      !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date  !< The vintage of the order of arithmetic and expressions
                                      !! to use for remapping.  Values below 20190101 recover the remapping
                                      !! answers from 2018, while higher values use more robust
                                      !! forms of the same remapping expressions.
  logical, optional, intent(in) :: better_speed_est !< If true, use a more robust estimate of the first
                                     !! mode speed as the starting point for iterations.
  real,    optional, intent(in) :: min_speed !< If present, set a floor in the first mode speed
                                     !! below which 0 is returned [L T-1 ~> m s-1].
  real,    optional, intent(in) :: wave_speed_tol !< The fractional tolerance for finding the
                                     !! wave speeds [nondim]
  real,    optional, intent(in) :: c1_thresh !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but are
                                       !! simply reported as 0 [L T-1 ~> m s-1].  A non-negative value
                                       !! must be specified for wave_speeds to be used (but not wave_speed).

end subroutine wave_speed_set_param
  end interface

end module MOM_wave_speed
