! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides the ocean stochastic equation of state
module MOM_stoch_eos

use MOM_diag_mediator,    only : register_diag_field, post_data, diag_ctrl
use MOM_error_handler,    only : MOM_error, FATAL
use MOM_file_parser,      only : get_param, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_hor_index,        only : hor_index_type
use MOM_isopycnal_slopes, only : vert_fill_TS
use MOM_random,           only : PRNG, random_2d_constructor, random_2d_norm
use MOM_restart,          only : MOM_restart_CS, register_restart_field, is_new_run, query_initialized
use MOM_time_manager,     only : time_type
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type
!use random_numbers_mod,  only : getRandomNumbers, initializeRandomNumberStream, randomNumberStream

implicit none ; private
#include <MOM_memory.h>

public MOM_stoch_eos_init
public MOM_stoch_eos_run
public stoch_EOS_register_restarts
public post_stoch_EOS_diags
public MOM_calc_varT

!> Describes parameters of the stochastic component of the EOS
!! correction, described in Stanley et al. JAMES 2020.
type, public :: MOM_stoch_eos_CS ; private
  real, allocatable :: l2_inv(:,:)  !< One over sum of the T cell side side lengths squared [L-2 ~> m-2]
  real, allocatable :: rgauss(:,:)  !< nondimensional random Gaussian [nondim]
  real      :: tfac = 0.27          !< Nondimensional decorrelation time factor, ~1/3.7 [nondim]
  real      :: amplitude = 0.624499 !< Nondimensional standard deviation of Gaussian [nondim]
  integer   :: seed                 !< PRNG seed
  type(PRNG)  ::  rn_CS             !< PRNG control structure
  real, allocatable :: pattern(:,:) !< Random pattern for stochastic EOS [nondim]
  real, allocatable :: phi(:,:)     !< temporal correlation stochastic EOS [nondim]
  logical :: use_stoch_eos!< If true, use the stochastic equation of state (Stanley et al. 2020)
  real :: stanley_coeff   !< Coefficient correlating the temperature gradient
                          !! and SGS T variance [nondim]; if <0, turn off scheme in all codes
  real :: stanley_a       !< a in exp(aX) in stochastic coefficient [nondim]
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  !>@{ Diagnostic IDs
  integer :: id_stoch_eos  = -1, id_stoch_phi  = -1, id_tvar_sgs = -1
  !>@}

end type MOM_stoch_eos_CS


  interface
logical module function MOM_stoch_eos_init(Time, G, GV, US, param_file, diag, CS, restart_CS)
  type(time_type),         intent(in)    :: Time       !< Time for stochastic process
  type(ocean_grid_type),   intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(diag_ctrl), target, intent(inout) :: diag       !< Structure used to control diagnostics
  type(MOM_stoch_eos_CS),  intent(inout) :: CS         !< Stochastic control structure
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.

  ! local variables

end function MOM_stoch_eos_init
module subroutine stoch_EOS_register_restarts(HI, param_file, CS, restart_CS)
  type(hor_index_type),    intent(in)    :: HI         !< Horizontal index structure
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(MOM_stoch_eos_CS),  intent(inout) :: CS         !< Stochastic control structure
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.

end subroutine stoch_EOS_register_restarts
module subroutine MOM_stoch_eos_run(G, u, v, delt, Time, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: delt !< Time step size for AR1 process [T ~> s].
  type(time_type),         intent(in)    :: Time !< Time for stochastic process
  type(MOM_stoch_eos_CS),  intent(inout) :: CS   !< Stochastic control structure

  ! local variables

  ! Return without doing anything if this capability is not enabled.
end subroutine MOM_stoch_eos_run
module subroutine post_stoch_EOS_diags(CS, tv, diag)
  type(MOM_stoch_eos_CS), intent(in) :: CS  !< Stochastic control structure
  type(thermo_var_ptrs),  intent(in) :: tv  !< Thermodynamics structure
  type(diag_ctrl),        intent(inout) :: diag !< Structure to control diagnostics

end subroutine post_stoch_EOS_diags
module subroutine MOM_calc_varT(G, GV, US, h, tv, CS, dt)
  type(ocean_grid_type),   intent(in)   :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)   :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)   :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                          intent(in)    :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),  intent(inout) :: tv  !< Thermodynamics structure
  type(MOM_stoch_eos_CS), intent(inout) :: CS  !< Stochastic control structure
  real,                   intent(in)    :: dt  !< Time increment [T ~> s]

  ! local variables
end subroutine MOM_calc_varT
  end interface

end module MOM_stoch_eos
