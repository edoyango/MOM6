! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module ocn_cpl_indices

  use mct_mod,              only: mct_avect_init, mct_avect_indexra, mct_aVect_clean, mct_aVect
  use seq_flds_mod,         only: ice_ncat, seq_flds_i2o_per_cat
  use seq_flds_mod,         only: seq_flds_x2o_fields, seq_flds_o2x_fields

  implicit none ; public

  !> Structure with indices needed for MCT attribute vectors
  type cpl_indices_type
    ! ocean to coupler
    integer :: o2x_So_t          !< Surface potential temperature (deg C)
    integer :: o2x_So_u          !< Surface zonal velocity [m s-1]
    integer :: o2x_So_v          !< Surface meridional velocity [m s-1]
    integer :: o2x_So_s          !< Surface salinity (PSU)
    integer :: o2x_So_dhdx       !< Zonal slope in the sea surface height
    integer :: o2x_So_dhdy       !< Meridional lope in the sea surface height
    integer :: o2x_So_bldepth    !< Boundary layer depth (m)
    integer :: o2x_Fioo_q        !< Ocean melt and freeze potential (W/m2)
    integer :: o2x_Faoo_fco2_ocn !< CO2 flux
    integer :: o2x_Faoo_fdms_ocn !< DMS flux

    ! coupler to ocean
    integer :: x2o_Si_ifrac      !< Fractional ice wrt ocean
    integer :: x2o_So_duu10n     !< 10m wind speed squared (m^2/s^2)
    integer :: x2o_Sa_pslv       !< Sea-level pressure (Pa)
    integer :: x2o_Sa_co2prog    !< Bottom atm level prognostic CO2
    integer :: x2o_Sa_co2diag    !< Bottom atm level diagnostic CO2
    integer :: x2o_Sw_lamult     !< Wave model langmuir multiplier
    integer :: x2o_Sw_ustokes    !< Surface Stokes drift, x-component
    integer :: x2o_Sw_vstokes    !< Surface Stokes drift, y-component
    integer :: x2o_Foxx_taux     !< Zonal wind stress (W/m2)
    integer :: x2o_Foxx_tauy     !< Meridonal wind stress (W/m2)
    integer :: x2o_Foxx_swnet    !< Net short-wave heat flux (W/m2)
    integer :: x2o_Foxx_sen      !< Sensible heat flux (W/m2)
    integer :: x2o_Foxx_lat      !< Latent heat flux  (W/m2)
    integer :: x2o_Foxx_lwup     !< Longwave radiation, up (W/m2)
    integer :: x2o_Faxa_lwdn     !< Longwave radiation, down (W/m2)
    integer :: x2o_Faxa_swvdr    !< Visible, direct shortwave (W/m2)
    integer :: x2o_Faxa_swvdf    !< Visible, diffuse shortwave (W/m2)
    integer :: x2o_Faxa_swndr    !< near-IR, direct shortwave (W/m2)
    integer :: x2o_Faxa_swndf    !< near-IR, direct shortwave (W/m2)
    integer :: x2o_Fioi_melth    !< Heat flux from snow & ice melt (W/m2)
    integer :: x2o_Fioi_meltw    !< Water flux from sea ice and snow melt (kg/m2/s)
    integer :: x2o_Fioi_bcpho    !< Black Carbon hydrophobic release from sea ice component
    integer :: x2o_Fioi_bcphi    !< Black Carbon hydrophilic release from sea ice component
    integer :: x2o_Fioi_flxdst   !< Dust release from sea ice component
    integer :: x2o_Fioi_salt     !< Salt flux    (kg(salt)/m2/s)
    integer :: x2o_Foxx_evap     !< Evaporation flux  (kg/m2/s)
    integer :: x2o_Faxa_prec     !< Total precipitation flux (kg/m2/s)
    integer :: x2o_Faxa_snow     !< Water flux due to snow (kg/m2/s)
    integer :: x2o_Faxa_rain     !< Water flux due to rain (kg/m2/s)
    integer :: x2o_Faxa_bcphidry !< Black   Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_bcphodry !< Black   Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_bcphiwet !< Black   Carbon hydrophilic wet deposition
    integer :: x2o_Faxa_ocphidry !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_ocphodry !< Organic Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_ocphiwet !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_dstwet1  !< Size 1 dust -- wet deposition
    integer :: x2o_Faxa_dstwet2  !< Size 2 dust -- wet deposition
    integer :: x2o_Faxa_dstwet3  !< Size 3 dust -- wet deposition
    integer :: x2o_Faxa_dstwet4  !< Size 4 dust -- wet deposition
    integer :: x2o_Faxa_dstdry1  !< Size 1 dust -- dry deposition
    integer :: x2o_Faxa_dstdry2  !< Size 2 dust -- dry deposition
    integer :: x2o_Faxa_dstdry3  !< Size 3 dust -- dry deposition
    integer :: x2o_Faxa_dstdry4  !< Size 4 dust -- dry deposition
    integer :: x2o_Foxx_rofl     !< River runoff flux (kg/m2/s)
    integer :: x2o_Foxx_rofi     !< Ice runoff flux (kg/m2/s)

    ! optional per thickness category fields
    integer, dimension(:), allocatable :: x2o_frac_col      !< Fraction of ocean cell, per column
    integer, dimension(:), allocatable :: x2o_fracr_col     !< Fraction of ocean cell used  in radiation computations,
                                                            !! per column
    integer, dimension(:), allocatable :: x2o_qsw_fracr_col !< qsw * fracr, per column
  end type cpl_indices_type

  public :: cpl_indices_init

!=======================================================================

  interface
  module subroutine cpl_indices_init(ind)
    type(cpl_indices_type), intent(inout) :: ind !< Structure with coupler indices and vectors

    ! Local Variables

    ! create temporary attribute vectors
  end subroutine cpl_indices_init
  end interface

end module ocn_cpl_indices
