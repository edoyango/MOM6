! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Accelerations due to the Coriolis force and momentum advection
module MOM_CoriolisAdv

!> \author Robert Hallberg, April 1994 - June 2002

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary, only : OBC_VORTICITY_ZERO, OBC_VORTICITY_FREESLIP
use MOM_open_boundary, only : OBC_VORTICITY_COMPUTED, OBC_VORTICITY_SPECIFIED
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : accel_diag_ptrs, porous_barrier_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS

implicit none ; private

public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end, CoriolisAdv_stencil

#include <MOM_memory.h>

!> Control structure for mom_coriolisadv
type, public :: CoriolisAdv_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  integer :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! - WENOVI3RD_PV_ENSTRO    - 3rd-order WENO scheme for PV reconstruction
                             !! - WENOVI5TH_PV_ENSTRO    - 5th-order WENO scheme for PV reconstruction
                             !! - WENOVI7TH_PV_ENSTRO    - 7th-order WENO scheme for PV reconstruction
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
  integer :: KE_Scheme       !< KE_SCHEME selects the discretization for
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  logical :: KE_use_limiter  !< If true, use the Koren limiter for KE_UP3 scheme
  integer :: PV_Adv_Scheme   !< PV_ADV_SCHEME selects the discretization for PV advection
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
  real    :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                             !! This must be greater than 2.0, and is 4.0 by default.
  real    :: wt_lin_blend    !< A weighting value beyond which the blending between
                             !! Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                             !! This must be between 1 and 1e-15, often 1/8.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Coriolis  !< If true, the Coriolis terms at u points are
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
  logical :: Coriolis_En_Dis !< If CORIOLIS_EN_DIS is defined, two estimates of
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
  logical :: weno_velocity_smooth !< If true, use velocity to compute the smoothness indicator for WENO
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_rv = -1, id_PV = -1, id_gKEu = -1, id_gKEv = -1
  integer :: id_rvxu = -1, id_rvxv = -1
  ! integer :: id_hf_gKEu    = -1, id_hf_gKEv    = -1
  integer :: id_hf_gKEu_2d = -1, id_hf_gKEv_2d = -1
  integer :: id_intz_gKEu_2d = -1, id_intz_gKEv_2d = -1
  ! integer :: id_hf_rvxu    = -1, id_hf_rvxv    = -1
  integer :: id_hf_rvxu_2d = -1, id_hf_rvxv_2d = -1
  integer :: id_h_gKEu = -1, id_h_gKEv = -1
  integer :: id_h_rvxu = -1, id_h_rvxv = -1
  integer :: id_intz_rvxu_2d = -1, id_intz_rvxv_2d = -1
  integer :: id_CAuS = -1, id_CAvS = -1
  !>@}
end type CoriolisAdv_CS

!>@{ Enumeration values for Coriolis_Scheme
integer, parameter :: SADOURNY75_ENERGY = 1
integer, parameter :: ARAKAWA_HSU90     = 2
integer, parameter :: ROBUST_ENSTRO     = 3
integer, parameter :: SADOURNY75_ENSTRO = 4
integer, parameter :: ARAKAWA_LAMB81    = 5
integer, parameter :: AL_BLEND          = 6
integer, parameter :: wenovi7th_PV_ENSTRO = 7
integer, parameter :: wenovi5th_PV_ENSTRO = 8
integer, parameter :: wenovi3rd_PV_ENSTRO = 9
character*(20), parameter :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
character*(20), parameter :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
character*(20), parameter :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
character*(20), parameter :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
character*(20), parameter :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
character*(20), parameter :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
character*(20), parameter :: WENOVI7TH_PV_ENSTRO_STRING = "WENOVI7TH_PV_ENSTRO"
character*(20), parameter :: WENOVI5TH_PV_ENSTRO_STRING = "WENOVI5TH_PV_ENSTRO"
character*(20), parameter :: WENOVI3RD_PV_ENSTRO_STRING = "WENOVI3RD_PV_ENSTRO"
!>@}
!>@{ Enumeration values for KE_Scheme
integer, parameter :: KE_ARAKAWA        = 10
integer, parameter :: KE_SIMPLE_GUDONOV = 11
integer, parameter :: KE_GUDONOV        = 12
integer, parameter :: KE_UP3            = 13
character*(20), parameter :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
character*(20), parameter :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
character*(20), parameter :: KE_GUDONOV_STRING = "KE_GUDONOV"
character*(20), parameter :: KE_UP3_STRING = "KE_UP3"
!>@}
!>@{ Enumeration values for PV_Adv_Scheme
integer, parameter :: PV_ADV_CENTERED   = 21
integer, parameter :: PV_ADV_UPWIND1    = 22
character*(20), parameter :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
character*(20), parameter :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"
!>@}


  interface
module subroutine CorAdCalc(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS

  ! Local variables

                ! surrounding an h grid point.  At small scales, a = q/4,
                ! b = q/4, etc.  All are in [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1],
                ! and use the indexing of the corresponding u point.

end subroutine CorAdCalc
module subroutine gradKE(u, v, h, KE, KEx, KEy, G, GV, US, CS)
  type(ocean_grid_type),             intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),           intent(in)  :: GV  !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)  :: u   !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: v   !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out) :: KE  !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(out) :: KEx !< Zonal acceleration due to kinetic
                                                        !! energy gradient [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: KEy !< Meridional acceleration due to kinetic
                                                        !! energy gradient [L T-2 ~> m s-2]
  type(unit_scale_type),             intent(in)  :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),              intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  ! Local variables

end subroutine gradKE
module subroutine UP3_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [l t-1 ~> m s-1] or [l2 t-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q at point i-1/2 [A ~> a]

end subroutine UP3_reconstruction
module subroutine UP3_Koren_limiter_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q on point i-1/2 [A ~> a]

end subroutine UP3_Koren_limiter_reconstruction
module function fac_fn(tau, b) result(fac)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: fac               !< The factor for the weight [nondim]

end function fac_fn
module subroutine weno_three_h_weight_reconstruction(q4, h4, u4, &
                                              h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q4(4)   !< Tracer value times thickness on points i-2, i-1, i, i+1 [A ~> a]
    real, intent(in)    :: h4(4)   !< Thickness values on points i-2, i-1, i, i+1 [L ~> m]
    real, optional, intent(in)    :: u4(4) !< Velocity values on points i-2, i-1, i, i+1
                                                    !![L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u              !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    real, intent(inout) :: qr             !< Reconstruction of tracer q on point i-1/2 [A ~> a]
    logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute smoothness indicator

end subroutine weno_three_h_weight_reconstruction
module subroutine weno_three_weight(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Smoothness indicator for this stencil [A2 ~> a2]

end subroutine weno_three_weight
module subroutine weno_three_reconstruction_0(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A2 ~> a2]

end subroutine weno_three_reconstruction_0
module subroutine weno_three_reconstruction_1(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A ~> a]

end subroutine weno_three_reconstruction_1
module subroutine weno_five_h_weight_reconstruction(q6, h6, u6, &
                                             h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q6(6)
    !< Tracer values on points i-3, i-2, i-1, i, i+1, i+2 [A ~> a]
    real, intent(in)    :: h6(6)
    !< Thickness values on points i-3, i-2, i-1, i, i+1, i+2 [L ~> m]
    real, optional, intent(in)    :: u6(6)
    !< Velocity values on points i-3, i-2, i-1, i, i+1, i+2 [L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u                      !< Velocity or thickness flux on point i-1/2
                                                  !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    logical, intent(in) :: velocity_smoothing     !< If ture, use velocity to compute the smoothness indicator
    real, intent(inout) :: qr                     !< Reconstruction of tracer q on point i-1/2 [A ~> a]

end subroutine weno_five_h_weight_reconstruction
module subroutine weno_five_weight_0(q3, w0)
  real, intent(in) :: q3(3)       !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w0       !< Smoothness indicator for this stencil [A2 ~> a2]

end subroutine weno_five_weight_0
module subroutine weno_five_weight_1(q3, w1)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w1        !< Smoothness indicator for this stencil [A2 ~> a2]

end subroutine weno_five_weight_1
module subroutine weno_five_weight_2(q3, w2)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w2        !< Smoothness indicator for this stencil [A2 ~> a2]

end subroutine weno_five_weight_2
module subroutine weno_five_reconstruction_0(q3, p0)
  real, intent(in) :: q3(3)        !< Tracer values on three points [A ~> a]
  real, intent(inout) :: p0        !< Reconstruction of the quantity [A ~> a]

end subroutine weno_five_reconstruction_0
module subroutine weno_five_reconstruction_1(q3, p1)
  real, intent(in) :: q3(3)         !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p1         !< Reconstruction of the quantity [A ~> a]

end subroutine weno_five_reconstruction_1
module subroutine weno_five_reconstruction_2(q3, p2)
  real, intent(in) :: q3(3)          !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p2          !< Reconstruction of the quantity [A ~> a]

end subroutine weno_five_reconstruction_2
module subroutine weno_seven_h_weight_reconstruction(q8, h8, u8, &
                                            h_tiny, u, qr, velocity_smoothing)
  real, intent(in)    :: q8(8)
  !< Tracer values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3
  real, intent(in)    :: h8(8)
  !< Thickness on the same tracer points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L ~> m]
  real, optional, intent(in)    :: u8(8)
  !< Velocity values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L T-1 ~> m s-1]
  real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
  real, intent(in)    :: u    !< Velocity or thickness flux on point i-1/2
                              !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute the smoothness indicator
  real, intent(inout) :: qr   !< Reconstruction of tracer q on point i-1/2 [A ~> a]

end subroutine weno_seven_h_weight_reconstruction
module subroutine weno_seven_weight_0(q4, w0)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w0          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
end subroutine weno_seven_weight_0
module subroutine weno_seven_weight_1(q4, w1)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w1          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
end subroutine weno_seven_weight_1
module subroutine weno_seven_weight_2(q4, w2)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w2           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
end subroutine weno_seven_weight_2
module subroutine weno_seven_weight_3(q4, w3)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w3           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
end subroutine weno_seven_weight_3
module subroutine weno_seven_reconstruction_0(q4, p0)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p0            !< Reconstruction of the quantity [A ~> a]

end subroutine weno_seven_reconstruction_0
module subroutine weno_seven_reconstruction_1(q4, p1)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p1            !< Reconstruction of the quantity [A ~> a]

end subroutine weno_seven_reconstruction_1
module subroutine weno_seven_reconstruction_2(q4, p2)
  real, intent(in) :: q4(4)             !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p2             !< Reconstruction of the quantity [A ~> a]

end subroutine weno_seven_reconstruction_2
module subroutine weno_seven_reconstruction_3(q4, p3)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p3            !< Reconstruction of the quantity [A ~> a]

end subroutine weno_seven_reconstruction_3
module function CoriolisAdv_stencil(CS) result(stencil)
  type(CoriolisAdv_CS), intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  integer :: stencil  !< The halo stencil size for the Coriolis advection scheme

end function CoriolisAdv_stencil
module subroutine CoriolisAdv_init(Time, G, GV, US, param_file, diag, AD, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(accel_diag_ptrs),   target, intent(inout) :: AD !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(inout) :: CS   !< Control structure for MOM_CoriolisAdv
  ! Local variables
! This include declares and sets the variable "version".

end subroutine CoriolisAdv_init
module subroutine CoriolisAdv_end(CS)
  type(CoriolisAdv_CS), intent(inout) :: CS !< Control structure for MOM_CoriolisAdv
end subroutine CoriolisAdv_end
  end interface

end module MOM_CoriolisAdv
