! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements an ANN-based mesoscale streamfunction parameterization for use
!! with isopycnal height diffusion in MOM_thickness_diffuse.
!!
!! The network reads a nondimensionalized stencil of density gradients,
!! strain rate components, and relative vorticity, and returns two density
!! flux components at the cell center. The dimensionalization in
!! meso_sfn_ANN_compute (multiplication by rho_grad_mag * vel_grad_mag *
!! areaT * ann_coeff) must match the nondimensionalization used when the
!! network was trained -- changing one without the other will produce
!! garbage fluxes. The training procedure is the implicit contract.
!!
!! Density fluxes are converted to a velocity-scale streamfunction
!! Upsilon (Ferrari et al. 2010) by dividing by the local 3-D density
!! gradient magnitude; a configurable clamp acts on Upsilon so the cap is
!! grid-independent. The volume-transport streamfunction passed back to
!! thickness_diffuse is Upsilon * dy_Cu (or dx_Cv), matching MOM6's
!! Sfn_unlim convention.
module MOM_meso_sfn_ANN

use MOM_ANN,              only : ANN_init, ANN_apply_array_sio, ANN_end, ANN_CS
use MOM_diag_mediator,    only : post_data, register_diag_field, diag_ctrl, time_type
use MOM_error_handler,    only : MOM_error, FATAL
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_isopycnal_slopes, only : calc_isoneutral_slopes
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_domains,          only : pass_vector

implicit none ; private

#include <MOM_memory.h>

public :: meso_sfn_ANN_init, meso_sfn_ANN_compute, meso_sfn_ANN_end

!> Control structure for meso-scale streamfunction ANN parameterization
type, public :: MESO_SFN_ANN_CS; private
  logical :: initialized = .false. !< If true, the module has been initialized.
  logical :: debug !< if true, write verbose checksums for debugging purposes.

  real :: ann_coeff  !< Coefficient to multiply the ANN output by.
  real    :: kappa_smooth        !< Vertical diffusivity used to interpolate more sensible values
                                 !! of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  integer :: ann_window !< Size of the window used in the ANN model.

  type(ANN_CS) :: ann_rho_flux !< ANN instance for off-diagonal and diagonal stress
  character(len=200) :: ann_file_rho_flux !< Path to netcdf file with ANN
  real :: min_dist_from_boundary  !< Minimum distance from bottom for valid interface [Z ~> m]
  real :: mag_grad_floor  !< Floor for density gradient magnitude [R Z-1 ~> kg m-4]
  real :: flux_clamp  !< Maximum magnitude of ANN output density flux [R L T-1 ~> kg m-2 s-1]
  real :: Upsilon_clamp !< Maximum magnitude of the velocity-scale streamfunction
                        !! Upsilon (Ferrari et al. 2010) [L Z T-1 ~> m2 s-1]
  type(diag_ctrl), pointer :: diag => NULL() !< structure used to regulate timing of diagnostics
  ! Diagnostic identifiers
  integer :: id_drdx_u !< Diagnostic id for zonal density gradient at u-points.
  integer :: id_drdy_v !< Diagnostic id for meridional density gradient at v-points.
  integer :: id_drdz_u !< Diagnostic id for vertical density gradient at u-points.
  integer :: id_drdz_v !< Diagnostic id for vertical density gradient at v-points.
  integer :: id_drdx_c !< Diagnostic id for zonal density gradient at center points.
  integer :: id_drdy_c !< Diagnostic id for meridional density gradient at center points.
  integer :: id_Fx_c   !< Diagnostic id for zonal density flux at center points.
  integer :: id_Fy_c   !< Diagnostic id for meridional density flux at center points.
  integer :: id_Fx_u   !< Diagnostic id for zonal density flux at u-points.
  integer :: id_Fy_v   !< Diagnostic id for meridional density flux at v-points.
  integer :: id_sfn_u  !< Diagnostic id for volume streamfunction at u-points.
  integer :: id_sfn_v  !< Diagnostic id for volume streamfunction at v-points.
end type MESO_SFN_ANN_CS


  interface
module subroutine meso_sfn_ANN_compute(h, e, sfn_u, sfn_v, G, GV, US, tv, CS, dt, u, v)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Layer thickness [Z ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: e      !< Layer thickness [Z ~> m or kg m-2]
  type(MESO_SFN_ANN_CS),                intent(inout) :: CS !< Control structure for thickness_flux_ann
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: sfn_u  !< Mesoscale volume streamfunction
                                                                     !! on u-points [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: sfn_v  !< Mesoscale volume streamfunction
                                                                     !! on v-points [Z L2 T-1 ~> m3 s-1]
  real,                                      intent(in)    :: dt     !< Model time step [T ~> s]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1].

  ! Local variables



                                                          !! at center points [R L-1 ~> kg m-4]

                                                       !! center points [R L T-1 ~> kg m-2 s-1]



                            ! roundoff; used to prevent division by zero [R L-1 ~> kg m-4]
                            ! roundoff; used to prevent division by zero [T-1 ~> s-1]

end subroutine meso_sfn_ANN_compute
module subroutine center_grad_rho(drdx_u, drdy_v, drdx_c, drdy_c, G, GV, CS)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(MESO_SFN_ANN_CS),                intent(inout) :: CS !< Control structure for thickness_flux_ann
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: drdx_u !< Zonal density gradient
                                                                    !! at u-points [R L-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: drdy_v !< Meridional density gradient
                                                                    !! at v-points [R L-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: drdx_c !< Zonal density gradient
                                                                       !! at center [R L-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: drdy_c !< Meridional density gradient
                                                                       !! at center [R L-1 ~> kg m-4]


end subroutine center_grad_rho
module subroutine center2uv(var1_c, var2_c, var1_u, var2_v, G, GV)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: var1_c !< Variable at center points [arbitrary]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: var2_c !< Variable at center points [arbitrary]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout)   :: var1_u !< Variable at u points [arbitrary]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout)   :: var2_v !< Variable at v points [arbitrary]


end subroutine center2uv
module subroutine vel_gradients(u, v, G, GV, dudx, dudy, dvdx, dvdy, CS)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),intent(in)    :: u   !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),intent(in)    :: v   !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: dudx !< du/dx [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: dvdy !< dv/dy [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: dudy !< du/dy [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: dvdx !< dv/dx [T-1 ~> s-1]
  type(MESO_SFN_ANN_CS), intent(in) :: CS !< Control structure for thickness_flux_ann

  ! Corner points

end subroutine vel_gradients
module subroutine calc_layered_density_gradients(G, GV, US, h, e, &
                                          drdx_u, drdy_v, drdz_u, drdz_v, halo, min_dist_from_boundary)
  type(ocean_grid_type),                       intent(in)  :: G
  type(verticalGrid_type),                     intent(in)  :: GV
  type(unit_scale_type),                       intent(in)  :: US
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)  :: h   ! Layer thickness [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: e   ! Interface heights [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: drdx_u ! [R L-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: drdy_v ! [R L-1 ~> kg m-4]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: drdz_u ! [R Z-1 ~> kg m-4]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: drdz_v ! [R Z-1 ~> kg m-4]
  integer,                                      intent(in)  :: halo
  real,                                         intent(in)  :: min_dist_from_boundary  ! Threshold for boundaries [Z]

  ! Local variables


end subroutine calc_layered_density_gradients
module subroutine meso_sfn_ANN_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(MESO_SFN_ANN_CS), intent(inout) :: CS !< Control structure for meso sfn ann

  ! Local variables

end subroutine meso_sfn_ANN_init
module subroutine meso_sfn_ANN_end(CS)
  type(MESO_SFN_ANN_CS), intent(inout) :: CS !< Control structure

  ! Deallocate anything that needs to be.
end subroutine meso_sfn_ANN_end
  end interface

end module MOM_meso_sfn_ANN
