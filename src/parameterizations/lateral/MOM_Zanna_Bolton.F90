! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates Zanna and Bolton 2020 parameterization
!! Implemented by Perezhogin P.A. Contact: pperezhogin@gmail.com
module MOM_Zanna_Bolton

use MOM_grid,          only : ocean_grid_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type, &
                              start_group_pass, complete_group_pass
use MOM_domains,       only : To_North, To_East
use MOM_domains,       only : pass_var, CORNER
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_ANN,           only : ANN_init, ANN_apply_array_sio, ANN_end, ANN_CS

implicit none ; private

#include <MOM_memory.h>

public ZB2020_lateral_stress, ZB2020_init, ZB2020_end, ZB2020_copy_gradient_and_thickness

!> Control structure for Zanna-Bolton-2020 parameterization.
type, public :: ZB2020_CS ; private
  ! Parameters
  real      :: amplitude      !< The nondimensional scaling factor in ZB model,
                              !! typically 0.1 - 10 [nondim].
  integer   :: ZB_type        !< Select how to compute the trace part of ZB model:
                              !! 0 - both deviatoric and trace components are computed
                              !! 1 - only deviatoric component is computed
                              !! 2 - only trace component is computed
  integer   :: ZB_cons        !< Select a discretization scheme for ZB model
                              !! 0 - non-conservative scheme
                              !! 1 - conservative scheme for deviatoric component
  integer   :: HPF_iter       !< Number of sharpening passes for the Velocity Gradient (VG) components
                              !! in ZB model.
  integer   :: Stress_iter    !< Number of smoothing passes for the Stress tensor components
                              !! in ZB model.
  real      :: Klower_R_diss  !< Attenuation of
                              !! the ZB parameterization in the regions of
                              !! geostrophically-unbalanced flows (Klower 2018, Juricke2020,2019)
                              !! Subgrid stress is multiplied by 1/(1+(shear/(f*R_diss)))
                              !! R_diss=-1: attenuation is not used; typical value R_diss=1.0 [nondim]
  integer   :: Klower_shear   !< Type of expression for shear in Klower formula
                              !! 0: sqrt(sh_xx**2 + sh_xy**2)
                              !! 1: sqrt(sh_xx**2 + sh_xy**2 + vort_xy**2)
  integer   :: Marching_halo  !< The number of filter iterations per a single MPI
                              !! exchange

  real, dimension(:,:,:), allocatable :: &
          sh_xx,   & !< Horizontal tension (du/dx - dv/dy) in h (CENTER)
                     !! points including metric terms [T-1 ~> s-1]
          sh_xy,   & !< Horizontal shearing strain (du/dy + dv/dx) in q (CORNER)
                     !! points including metric terms [T-1 ~> s-1]
          vort_xy, & !< Vertical vorticity (dv/dx - du/dy) in q (CORNER)
                     !! points including metric terms [T-1 ~> s-1]
          hq         !< Thickness in CORNER points [H ~> m or kg m-2]

  real, dimension(:,:,:), allocatable :: &
          Txx,     & !< Subgrid stress xx component in h [L2 T-2 ~> m2 s-2]
          Tyy,     & !< Subgrid stress yy component in h [L2 T-2 ~> m2 s-2]
          Txy        !< Subgrid stress xy component in q [L2 T-2 ~> m2 s-2]

  real, dimension(:,:), allocatable :: &
          kappa_h, & !< Scaling coefficient in h points [L2 ~> m2]
          kappa_q    !< Scaling coefficient in q points [L2 ~> m2]

  real, allocatable ::    &
        ICoriolis_h(:,:), &  !< Inverse Coriolis parameter at h points [T ~> s]
        c_diss(:,:,:)        !< Attenuation parameter at h points
                             !! (Klower 2018, Juricke2019,2020) [nondim]

  real, dimension(:,:), allocatable ::    &
        maskw_h,  & !< Mask of land point at h points multiplied by filter weight [nondim]
        maskw_q     !< Same mask but for q points [nondim]

  logical :: use_ann  !< If True, momentum fluxes are inferred with ANN
  integer :: stencil_size  !< Default is 3x3
  type(ANN_CS) :: ann_Tall !< ANN instance for off-diagonal and diagonal stress
  character(len=200) :: ann_file_Tall !< Path to netcdf file with ANN
  real :: subroundoff_shear !< Small dimensional constant for save division by zero [T-1 ~> s-1]

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_ZB2020u = -1, id_ZB2020v = -1, id_KE_ZB2020 = -1
  integer :: id_Txx = -1
  integer :: id_Tyy = -1
  integer :: id_Txy = -1
  integer :: id_cdiss = -1
  !>@}

  !>@{ CPU time clock IDs
  integer :: id_clock_module
  integer :: id_clock_copy
  integer :: id_clock_cdiss
  integer :: id_clock_stress
  integer :: id_clock_stress_ANN
  integer :: id_clock_divergence
  integer :: id_clock_mpi
  integer :: id_clock_filter
  integer :: id_clock_post
  integer :: id_clock_source
  !>@}

  !>@{ MPI group passes
  type(group_pass_type) :: &
      pass_Tq, pass_Th, &        !< handles for halo passes of Txy and Txx, Tyy
      pass_xx, pass_xy           !< handles for halo passes of sh_xx and sh_xy, vort_xy
  integer :: Stress_halo = -1, & !< The halo size in filter of the stress tensor
             HPF_halo = -1       !< The halo size in filter of the velocity gradient
  !>@}

end type ZB2020_CS


  interface
module subroutine ZB2020_init(Time, G, GV, US, param_file, diag, CS, use_ZB2020)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ZB2020_CS),         intent(inout) :: CS         !< ZB2020 control structure.
  logical,                 intent(out)   :: use_ZB2020 !< If true, turns on ZB scheme.

                              ! but small compared to Coriolis parameter [T-1 ~> s-1]


  ! This include declares and sets the variable "version".

end subroutine ZB2020_init
module subroutine ZB2020_end(CS)
  type(ZB2020_CS), intent(inout) :: CS  !< ZB2020 control structure.

end subroutine ZB2020_end
module subroutine ZB2020_copy_gradient_and_thickness(sh_xx, sh_xy, vort_xy, hq, &
                                       G, GV, CS, k)
  type(ocean_grid_type),         intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(inout) :: CS     !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: sh_xy       !< horizontal shearing strain (du/dy + dv/dx)
                              !! including metric terms [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: vort_xy     !< Vertical vorticity (dv/dx - du/dy)
                              !! including metric terms [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: hq          !< harmonic mean of the harmonic means
                              !! of the u- & v point thicknesses [H ~> m or kg m-2]

  real, dimension(SZI_(G),SZJ_(G)), &
    intent(in) :: sh_xx       !< horizontal tension (du/dx - dv/dy)
                              !! including metric terms [T-1 ~> s-1]

  integer, intent(in) :: k    !< The vertical index of the layer to be passed.


end subroutine ZB2020_copy_gradient_and_thickness
module subroutine ZB2020_lateral_stress(u, v, h, diffu, diffv, G, GV, CS, &
                             dx2h, dy2h, dx2q, dy2q)
  type(ocean_grid_type),         intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(inout) :: CS !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                 intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2].

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                        intent(inout) :: diffu   !< Zonal acceleration due to eddy viscosity.
                                                 !! It is updated with ZB closure [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                        intent(inout) :: diffv   !< Meridional acceleration due to eddy viscosity.
                                                 !! It is updated with ZB closure [L T-2 ~> m s-2]

  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: dx2h    !< dx^2 at h points [L2 ~> m2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: dy2h    !< dy^2 at h points [L2 ~> m2]

  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: dx2q    !< dx^2 at q points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: dy2q    !< dy^2 at q points [L2 ~> m2]

end subroutine ZB2020_lateral_stress
module subroutine compute_c_diss(G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(inout) :: CS   !< ZB2020 control structure.



end subroutine compute_c_diss
module subroutine compute_stress(G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(inout) :: CS   !< ZB2020 control structure.



  ! Local variables



end subroutine compute_stress
module subroutine compute_stress_ANN_collocated(G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(inout) :: CS   !< ZB2020 control structure.


                                     ! number of horizontal grid points x
                                     ! (sh_xy, sh_xx, vort_xy) on a stencil    [nondim]
                                     ! output features number of horizontal grid points x
                                     ! (Txy,Txx,Tyy) [nondim]
                                     ! output features (Txy,Txx,Tyy) [L2 T-2 ~> m2 s-2]



end subroutine compute_stress_ANN_collocated
module subroutine compute_stress_divergence(u, v, h, diffu, diffv, dx2h, dy2h, dx2q, dy2q, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(in) :: CS   !< ZB2020 control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
        intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
        intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
        intent(in) :: h             !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
        intent(out) :: diffu           !< Zonal acceleration due to convergence of
                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
        intent(out) :: diffv           !< Meridional acceleration due to convergence
                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJ_(G)),           &
        intent(in) :: dx2h          !< dx^2 at h points [L2 ~> m2]
  real, dimension(SZI_(G),SZJ_(G)),           &
        intent(in) :: dy2h          !< dy^2 at h points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)),         &
        intent(in) :: dx2q          !< dx^2 at q points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)),         &
        intent(in) :: dy2q          !< dy^2 at q points [L2 ~> m2]

  ! Local variables

                          !! along-coordinate stress tensor for ZB model
                          !! [L T-2 ~> m s-2]
                          !! of along-coordinate stress tensor for ZB model
                          !! [L T-2 ~> m s-2]


                       ! roundoff and so neglected [H ~> m or kg m-2]


end subroutine compute_stress_divergence
module subroutine filter_velocity_gradients(G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G       !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV      !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(inout) :: CS   !< ZB2020 control structure.




end subroutine filter_velocity_gradients
module subroutine filter_stress(G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G       !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV      !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(inout) :: CS   !< ZB2020 control structure.


end subroutine filter_stress
module subroutine filter_hq(G, GV, CS, current_halo, remaining_iterations, q, h)
  type(ocean_grid_type),   intent(in) :: G       !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV      !< The ocean's vertical grid structure
  type(ZB2020_CS),         intent(in) :: CS      !< ZB2020 control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional,   &
           intent(inout) :: h !< Input/output array in h points [arbitrary]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), optional, &
           intent(inout) :: q !< Input/output array in q points [arbitrary]
  integer, intent(inout) :: current_halo         !< Currently available halo points
  integer, intent(inout) :: remaining_iterations !< The number of iterations to perform


end subroutine filter_hq
module subroutine filter_3D(x, maskw, isd, ied, jsd, jed, is, ie, js, je, nz, &
                     current_halo, remaining_iterations,               &
                     direction)
  integer, intent(in) :: isd !< Indices of array size
  integer, intent(in) :: ied !< Indices of array size
  integer, intent(in) :: jsd !< Indices of array size
  integer, intent(in) :: jed !< Indices of array size
  integer, intent(in) :: is  !< Indices of owned points
  integer, intent(in) :: ie  !< Indices of owned points
  integer, intent(in) :: js  !< Indices of owned points
  integer, intent(in) :: je  !< Indices of owned points
  integer, intent(in) :: nz  !< Vertical array size
  real, dimension(isd:ied,jsd:jed,nz), &
           intent(inout) :: x !< Input/output array [arbitrary]
  real, dimension(isd:ied,jsd:jed), &
           intent(in) :: maskw !< Mask array of land points divided by 16 [nondim]
  integer, intent(inout) :: current_halo         !< Currently available halo points
  integer, intent(inout) :: remaining_iterations !< The number of iterations to perform
  logical, intent(in)    :: direction            !< The direction of the first 1D filter



  ! Do as many iterations as needed and possible
end subroutine filter_3D
module subroutine compute_energy_source(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(in)  :: CS   !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                 intent(in) :: h     !< Layer thicknesses [H ~> m or kg m-2].

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: fx    !< Zonal acceleration due to convergence of
                                                     !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in) :: fy    !< Meridional acceleration due to convergence
                                                     !! of along-coordinate stress tensor [L T-2 ~> m s-2]

                                            ! [H L2 T-3 ~> m3 s-3 or W m-2]
                                            ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
                                            ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

                                            ! [H L2 T-1 ~> m3 s-1 or kg s-1].
                                            ! [H L2 T-1 ~> m3 s-1 or kg s-1].



end subroutine compute_energy_source
  end interface

end module MOM_Zanna_Bolton
