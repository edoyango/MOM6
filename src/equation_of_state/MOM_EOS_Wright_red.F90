! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state using the Wright 1997 expressions with reduced range of data.
module MOM_EOS_Wright_red

use MOM_EOS_base_type, only : EOS_base
use MOM_hor_index, only : hor_index_type

implicit none ; private

public Wright_red_EOS
public int_density_dz_wright_red, int_spec_vol_dp_wright_red
public avg_spec_vol_Wright_red

!>@{ Parameters in the Wright equation of state using the reduced range formula, which is a fit to the UNESCO
!    equation of state for the restricted range: -2 < theta < 30 [degC], 28 < S < 38 [PSU], 0  < p < 5e7 [Pa].

  ! Note that a0/a1 ~= 2028 [degC] ; a0/a2 ~= -6343 [PSU]
  !           b0/b1 ~= 165 [degC]  ; b0/b4 ~= 974 [PSU]
  !           c0/c1 ~= 216 [degC]  ; c0/c4 ~= -740 [PSU]
real, parameter :: a0 = 7.057924e-4  ! A parameter in the Wright alpha_0 fit [m3 kg-1]
real, parameter :: a1 = 3.480336e-7  ! A parameter in the Wright alpha_0 fit [m3 kg-1 degC-1]
real, parameter :: a2 = -1.112733e-7 ! A parameter in the Wright alpha_0 fit [m3 kg-1 PSU-1]
real, parameter :: b0 = 5.790749e8   ! A parameter in the Wright p_0 fit [Pa]
real, parameter :: b1 = 3.516535e6   ! A parameter in the Wright p_0 fit [Pa degC-1]
real, parameter :: b2 = -4.002714e4  ! A parameter in the Wright p_0 fit [Pa degC-2]
real, parameter :: b3 = 2.084372e2   ! A parameter in the Wright p_0 fit [Pa degC-3]
real, parameter :: b4 = 5.944068e5   ! A parameter in the Wright p_0 fit [Pa PSU-1]
real, parameter :: b5 = -9.643486e3  ! A parameter in the Wright p_0 fit [Pa degC-1 PSU-1]
real, parameter :: c0 = 1.704853e5   ! A parameter in the Wright lambda fit [m2 s-2]
real, parameter :: c1 = 7.904722e2   ! A parameter in the Wright lambda fit [m2 s-2 degC-1]
real, parameter :: c2 = -7.984422    ! A parameter in the Wright lambda fit [m2 s-2 degC-2]
real, parameter :: c3 = 5.140652e-2  ! A parameter in the Wright lambda fit [m2 s-2 degC-3]
real, parameter :: c4 = -2.302158e2  ! A parameter in the Wright lambda fit [m2 s-2 PSU-1]
real, parameter :: c5 = -3.079464    ! A parameter in the Wright lambda fit [m2 s-2 degC-1 PSU-1]
!>@}

!> The EOS_base implementation of the reduced range Wright 1997 equation of state
type, extends (EOS_base) :: Wright_red_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_Wright_red
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_Wright_red
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_Wright_red
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_Wright_red
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_Wright_red
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_Wright_red
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_Wright_red
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_Wright_red
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_Wright_red

  !> Local implementation of generic calculate_density_array for efficiency
  procedure :: calculate_density_array => calculate_density_array_Wright_red
  !> Local implementation of generic calculate_spec_vol_array for efficiency
  procedure :: calculate_spec_vol_array => calculate_spec_vol_array_Wright_red

end type Wright_red_EOS


  interface
real elemental module function density_elem_Wright_red(this, T, S, pressure)
  class(Wright_red_EOS), intent(in) :: this !< This EOS
  real,           intent(in) :: T        !< potential temperature relative to the surface [degC].
  real,           intent(in) :: S        !< salinity [PSU].
  real,           intent(in) :: pressure !< pressure [Pa].

  ! Local variables

end function density_elem_Wright_red
real elemental module function density_anomaly_elem_Wright_red(this, T, S, pressure, rho_ref)
  class(Wright_red_EOS), intent(in) :: this !< This EOS
  real, intent(in) :: T        !< potential temperature relative to the surface [degC].
  real, intent(in) :: S        !< salinity [PSU].
  real, intent(in) :: pressure !< pressure [Pa].
  real, intent(in) :: rho_ref  !< A reference density [kg m-3].

  ! Local variables

end function density_anomaly_elem_Wright_red
real elemental module function spec_vol_elem_Wright_red(this, T, S, pressure)
  class(Wright_red_EOS), intent(in) :: this !< This EOS
  real,           intent(in) :: T        !< potential temperature relative to the surface [degC]
  real,           intent(in) :: S        !< salinity [PSU]
  real,           intent(in) :: pressure !< pressure [Pa]

  ! Local variables
                  ! an offset to account for spv_ref

end function spec_vol_elem_Wright_red
real elemental module function spec_vol_anomaly_elem_Wright_red(this, T, S, pressure, spv_ref)
  class(Wright_red_EOS), intent(in) :: this !< This EOS
  real,           intent(in) :: T        !< potential temperature relative to the surface [degC]
  real,           intent(in) :: S        !< salinity [PSU]
  real,           intent(in) :: pressure !< pressure [Pa]
  real,           intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables
                  ! an offset to account for spv_ref

end function spec_vol_anomaly_elem_Wright_red
elemental module subroutine calculate_density_derivs_elem_Wright_red(this, T, S, pressure, drho_dT, drho_dS)
  class(Wright_red_EOS), intent(in) :: this   !< This EOS
  real,               intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,               intent(in)  :: S        !< Salinity [PSU]
  real,               intent(in)  :: pressure !< Pressure [Pa]
  real,               intent(out) :: drho_dT  !< The partial derivative of density with potential
                                              !! temperature [kg m-3 degC-1]
  real,               intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                              !! in [kg m-3 PSU-1]

  ! Local variables

end subroutine calculate_density_derivs_elem_Wright_red
elemental module subroutine calculate_density_second_derivs_elem_Wright_red(this, T, S, pressure, &
                              drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
  class(Wright_red_EOS), intent(in) :: this       !< This EOS
  real,               intent(in)    :: T          !< Potential temperature referenced to 0 dbar [degC]
  real,               intent(in)    :: S          !< Salinity [PSU]
  real,               intent(in)    :: pressure   !< Pressure [Pa]
  real,               intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 PSU-2]
  real,               intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                  !! to T [kg m-3 PSU-1 degC-1]
  real,               intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
  real,               intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real,               intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables
                  ! offset (p0 elsewhere) in the Wright EOS [Pa]

end subroutine calculate_density_second_derivs_elem_Wright_red
elemental module subroutine calculate_specvol_derivs_elem_Wright_red(this, T, S, pressure, dSV_dT, dSV_dS)
  class(Wright_red_EOS), intent(in) :: this     !< This EOS
  real,               intent(in)    :: T        !< Potential temperature [degC]
  real,               intent(in)    :: S        !< Salinity [PSU]
  real,               intent(in)    :: pressure !< Pressure [Pa]
  real,               intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                !! potential temperature [m3 kg-1 degC-1]
  real,               intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                !! salinity [m3 kg-1 PSU-1]

  ! Local variables

  !al0 = a0 + (a1*T + a2*S)
end subroutine calculate_specvol_derivs_elem_Wright_red
elemental module subroutine calculate_compress_elem_Wright_red(this, T, S, pressure, rho, drho_dp)
  class(Wright_red_EOS), intent(in) :: this   !< This EOS
  real,               intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,               intent(in)  :: S        !< Salinity [PSU]
  real,               intent(in)  :: pressure !< Pressure [Pa]
  real,               intent(out) :: rho      !< In situ density [kg m-3]
  real,               intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                              !! (also the inverse of the square of sound speed)
                                              !! [s2 m-2].

  ! Local variables

end subroutine calculate_compress_elem_Wright_red
module subroutine avg_spec_vol_Wright_red(T, S, p_t, dp, SpV_avg, start, npts)
  real, dimension(:), intent(in)    :: T         !< Potential temperature relative to the surface
                                                 !! [degC].
  real, dimension(:), intent(in)    :: S         !< Salinity [PSU].
  real, dimension(:), intent(in)    :: p_t       !< Pressure at the top of the layer [Pa]
  real, dimension(:), intent(in)    :: dp        !< Pressure change in the layer [Pa]
  real, dimension(:), intent(inout) :: SpV_avg   !< The vertical average specific volume
                                                 !! in the layer [m3 kg-1]
  integer,            intent(in)    :: start     !< the starting point in the arrays.
  integer,            intent(in)    :: npts      !< the number of values to calculate.

  ! Local variables

  !  alpha(j) = al0 + lambda / (pressure(j) + p0)
end subroutine avg_spec_vol_Wright_red
module subroutine EoS_fit_range_Wright_red(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(Wright_red_EOS), intent(in) :: this !< This EOS
  real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: S_max !< The maximum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_Wright_red
module subroutine int_density_dz_wright_red(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, dz_neglect, &
                                 MassWghtInterp, rho_scale, pres_scale, temp_scale, saln_scale, Z_0p)
  type(hor_index_type), intent(in)  :: HI       !< The horizontal index type for the arrays.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T        !< Potential temperature relative to the surface
                                                !! [C ~> degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S        !< Salinity [S ~> PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t      !< Height at the top of the layer in depth units [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b      !< Height at the top of the layer [Z ~> m].
  real,                 intent(in)  :: rho_ref  !< A mean density [R ~> kg m-3], that is subtracted
                                                !! out to reduce the magnitude of each of the integrals.
                                                !! (The pressure is calculated as p~=-z*rho_0*G_e.)
  real,                 intent(in)  :: rho_0    !< Density [R ~> kg m-3], that is used
                                                !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                                !! used in the equation of state.
  real,                 intent(in)  :: G_e      !< The Earth's gravitational acceleration
                                                !! [L2 Z-1 T-2 ~> m s-2].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dpa    !< The change in the pressure anomaly across the
                                                !! layer [R L2 T-2 ~> Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer
                                                !! of the pressure anomaly relative to the anomaly
                                                !! at the top of the layer [R Z L2 T-2 ~> Pa m].
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the x grid spacing [R L2 T-2 ~> Pa].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the y grid spacing [R L2 T-2 ~> Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT   !< The depth of the bathymetry [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: SSH      !< The sea surface height [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m].
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                                !! mass weighting to interpolate T/S in integrals
  real,       optional, intent(in)  :: rho_scale !< A multiplicative factor by which to scale density
                                                 !! from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real,       optional, intent(in)  :: pres_scale !< A multiplicative factor to convert pressure
                                                 !! into Pa [Pa T2 R-1 L-2 ~> 1].
  real,       optional, intent(in)  :: temp_scale  !< A multiplicative factor by which to scale
                            !! temperature into degC [degC C-1 ~> 1]
  real,       optional, intent(in)  :: saln_scale !< A multiplicative factor to convert pressure
                            !! into PSU [PSU S-1 ~> 1].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p      !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa].
                       ! pres_scale [R L2 T-2 Pa-1 ~> 1].

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
end subroutine int_density_dz_wright_red
module subroutine int_spec_vol_dp_wright_red(T, S, p_t, p_b, spv_ref, HI, dza, &
                                  intp_dza, intx_dza, inty_dza, halo_size, bathyP, P_surf, dP_neglect, &
                                  MassWghtInterp, SV_scale, pres_scale, temp_scale, saln_scale)
  type(hor_index_type), intent(in)  :: HI        !< The ocean's horizontal index type.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [C ~> degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S         !< Salinity [S ~> PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t       !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b       !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: spv_ref   !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1].
                            !! The calculation is mathematically identical with different values of
                            !! spv_ref, but this reduces the effects of roundoff.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza     !< The change in the geopotential anomaly across
                                                 !! the layer [L2 T-2 ~> m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                                                 !! the geopotential anomaly relative to the anomaly
                                                 !! at the bottom of the layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the x grid spacing
                                                 !! [L2 T-2 ~> m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the y grid spacing
                                                 !! [L2 T-2 ~> m2 s-2].
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate
                                                 !! dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP    !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: P_surf    !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_neglect !< A miniscule pressure change with
                                                 !! the same units as p_t [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                                 !! mass weighting to interpolate T/S in integrals
  real,       optional, intent(in)  :: SV_scale  !< A multiplicative factor by which to scale specific
                            !! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  real,       optional, intent(in)  :: pres_scale !< A multiplicative factor to convert pressure
                            !! into Pa [Pa T2 R-1 L-2 ~> 1].
  real,       optional, intent(in)  :: temp_scale  !< A multiplicative factor by which to scale
                            !! temperature into degC [degC C-1 ~> 1]
  real,       optional, intent(in)  :: saln_scale !< A multiplicative factor to convert pressure
                            !! into PSU [PSU S-1 ~> 1].

  ! Local variables
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2].

end subroutine int_spec_vol_dp_wright_red
module subroutine calculate_density_array_Wright_red(this, T, S, pressure, rho, start, npts, rho_ref)
  class(Wright_red_EOS), intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end subroutine calculate_density_array_Wright_red
module subroutine calculate_spec_vol_array_Wright_red(this, T, S, pressure, specvol, start, npts, spv_ref)
  class(Wright_red_EOS),  intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables

end subroutine calculate_spec_vol_array_Wright_red
  end interface

end module MOM_EOS_Wright_red
