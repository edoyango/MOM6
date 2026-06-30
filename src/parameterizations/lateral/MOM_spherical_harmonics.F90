! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Laplace's spherical harmonic transforms (SHT)
module MOM_spherical_harmonics
use MOM_coms_infra,    only : sum_across_PEs
use MOM_coms,          only : reproducing_sum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type

implicit none ; private

public spherical_harmonics_init, spherical_harmonics_end, order2index, calc_lmax
public spherical_harmonics_forward, spherical_harmonics_inverse

#include <MOM_memory.h>

!> Control structure for spherical harmonic transforms
type, public :: sht_CS ; private
  logical :: initialized = .False. !< True if this control structure has been initialized.
  integer :: ndegree !< Maximum degree of the spherical harmonics [nondim].
  integer :: lmax !< Number of associated Legendre polynomials of nonnegative m
                  !! [lmax=(ndegree+1)*(ndegree+2)/2] [nondim].
  real, allocatable :: cos_clatT(:,:) !< Precomputed cosine of colatitude at the t-cells [nondim].
  real, allocatable :: Pmm(:,:,:) !< Precomputed associated Legendre polynomials (m=n) at the t-cells [nondim].
  real, allocatable :: cos_lonT(:,:,:), & !< Precomputed cosine factors at the t-cells [nondim].
                       sin_lonT(:,:,:)    !< Precomputed sine factors at the t-cells [nondim].
  real, allocatable :: cos_lonT_wtd(:,:,:), & !< Precomputed area-weighted cosine factors at the t-cells [nondim]
                       sin_lonT_wtd(:,:,:)    !< Precomputed area-weighted sine factors at the t-cells [nondim]
  real, allocatable :: a_recur(:,:), & !< Precomputed recurrence coefficients a [nondim].
                       b_recur(:,:)    !< Precomputed recurrence coefficients b [nondim].
  logical :: reprod_sum !< True if use reproducible global sums
end type sht_CS

integer :: id_clock_sht=-1 !< CPU clock for SHT [MODULE]
integer :: id_clock_sht_forward=-1 !< CPU clock for forward transforms [ROUTINE]
integer :: id_clock_sht_inverse=-1  !< CPU clock for inverse transforms [ROUTINE]
integer :: id_clock_sht_global_sum=-1  !< CPU clock for global summation in forward transforms [LOOP]


  interface
module subroutine spherical_harmonics_forward(G, CS, var, Snm_Re, Snm_Im, Nd, tmp_scale)
  type(ocean_grid_type), intent(in)    :: G            !< The ocean's grid structure.
  type(sht_CS),          intent(inout) :: CS           !< Control structure for SHT
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: var          !< Input 2-D variable in arbitrary mks units [a]
                                                       !! or in arbitrary rescaled units [A ~> a] if
                                                       !! tmp_scale is present
  real,                  intent(out)   :: Snm_Re(:)    !< SHT coefficients for the real modes (cosine) in
                                                       !! the same arbitrary units as var [a] or [A ~> a]
  real,                  intent(out)   :: Snm_Im(:)    !< SHT coefficients for the imaginary modes (sine) in
                                                       !! the same arbitrary units as var [a] or [A ~> a]
  integer,     optional, intent(in)    :: Nd           !< Maximum degree of the spherical harmonics
                                                       !! overriding ndegree in the CS [nondim]
  real,        optional, intent(in)    :: tmp_scale    !< A temporary rescaling factor to convert
                                                       !! var to MKS units during the reproducing
                                                       !! sums [a A-1 ~> 1]
  ! local variables
end subroutine spherical_harmonics_forward
module subroutine spherical_harmonics_inverse(G, CS, Snm_Re, Snm_Im, var, Nd)
  type(ocean_grid_type), intent(in)  :: G            !< The ocean's grid structure.
  type(sht_CS),          intent(in)  :: CS           !< Control structure for SHT
  real,                  intent(in)  :: Snm_Re(:)    !< SHT coefficients for the real modes (cosine)
                                                     !! in arbitrary units [a] or [A ~> a]
  real,                  intent(in)  :: Snm_Im(:)    !< SHT coefficients for the imaginary modes (sine) in
                                                     !! the same arbitrary units as Snm_Re [a] or [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(out) :: var          !< Output 2-D variable in the same arbitrary units
                                                     !! as Snm_Re and Snm_Im [a] or [A ~> a]
  integer,     optional, intent(in)  :: Nd           !< Maximum degree of the spherical harmonics
                                                     !! overriding ndegree in the CS [nondim]
  ! local variables

end subroutine spherical_harmonics_inverse
module subroutine spherical_harmonics_init(G, param_file, CS)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure indicating
  type(sht_CS), intent(inout)       :: CS !< Control structure for spherical harmonic transforms

  ! local variables
  ! This include declares and sets the variable "version".

end subroutine spherical_harmonics_init
module subroutine spherical_harmonics_end(CS)
  type(sht_CS), intent(inout) :: CS !< Control structure for spherical harmonic transforms

end subroutine spherical_harmonics_end
module function calc_lmax(Nd) result(lmax)
  integer :: lmax           !< Number of real spherical harmonic modes [nondim]
  integer, intent(in) :: Nd !< Maximum degree [nondim]

end function calc_lmax
module function order2index(m, Nd) result(l)
  integer :: l              !< One-dimensional index number [nondim]
  integer, intent(in) :: m  !< Current order number [nondim]
  integer, intent(in) :: Nd !< Maximum degree [nondim]

end function order2index
  end interface

end module MOM_spherical_harmonics
