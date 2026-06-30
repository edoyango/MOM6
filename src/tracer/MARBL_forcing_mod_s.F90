submodule (MARBL_forcing_mod) MARBL_forcing_mod_s
#include <MOM_memory.h>
  implicit none
contains
  module procedure MARBL_forcing_init
    character(len=40)  :: mdl = "MARBL_forcing_mod"  ! This module's name.
    character(len=15)  :: atm_co2_opt
    character(len=200) :: err_message
    if (associated(CS)) then
      call MOM_error(WARNING, "marbl_forcing_init called with an associated control structure.")
      return
    endif

    allocate(CS)
    CS%diag => diag

    CS%use_marbl_tracers = .true.
    if (.not. use_marbl) then
      CS%use_marbl_tracers = .false.
      return
    endif

    call get_param(param_file, mdl, "DUST_RATIO_THRES", CS%dust_ratio_thres, &
        "coarse/fine dust ratio threshold", units="1", default=69.00594)
    call get_param(param_file, mdl, "DUST_RATIO_TO_FE_BIOAVAIL_FRAC", CS%dust_ratio_to_fe_bioavail_frac, &
        "ratio of dust to iron bioavailability fraction", units="1", default=1./366.314)
    call get_param(param_file, mdl, "FE_BIOAVAIL_FRAC_OFFSET", CS%fe_bioavail_frac_offset, &
        "offset for iron bioavailability fraction", units="1", default=0.0146756)
    call get_param(param_file, mdl, "ATM_FE_TO_BC_RATIO", CS%atm_fe_to_bc_ratio, &
        "atmospheric iron to black carbon ratio", units="1", default=1.)
    call get_param(param_file, mdl, "ATM_BC_FE_BIOAVAIL_FRAC", CS%atm_bc_fe_bioavail_frac, &
        "atmospheric black carbon to iron bioavailablity fraction ratio", units="1", default=0.06)
    call get_param(param_file, mdl, "SEAICE_FE_TO_BC_RATIO", CS%seaice_fe_to_bc_ratio, &
        "sea-ice iron to black carbon ratio", units="1", default=1.)
    call get_param(param_file, mdl, "SEAICE_BC_FE_BIOAVAIL_FRAC", CS%seaice_bc_fe_bioavail_frac, &
        "sea-ice black carbon to iron bioavailablity fraction ratio", units="1", default=0.06)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_FINE_DUST", CS%iron_frac_in_atm_fine_dust, &
        "Fraction of fine dust from the atmosphere that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_COARSE_DUST", CS%iron_frac_in_atm_coarse_dust, &
        "Fraction of coarse dust from the atmosphere that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_SEAICE_DUST", CS%iron_frac_in_seaice_dust, &
        "Fraction of dust from sea ice that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "ATM_CO2_OPT", atm_co2_opt, &
        "Source of atmospheric CO2 [constant, diagnostic, or prognostic]", &
        default="constant")
    select case (trim(atm_co2_opt))
      case("prognostic")
        CS%atm_co2_iopt = atm_co2_prognostic_iopt
      case("diagnostic")
        CS%atm_co2_iopt = atm_co2_diagnostic_iopt
      case("constant")
        CS%atm_co2_iopt = atm_co2_constant_iopt
      case DEFAULT
        write(err_message, "(3A)") "'", trim(atm_co2_opt), "' is not a valid ATM_CO2_OPT value"
        call MOM_error(FATAL, err_message)
    end select
    if (CS%atm_co2_iopt == atm_co2_constant_iopt) then
      call get_param(param_file, mdl, "ATM_CO2_CONST", CS%atm_co2_const, &
          "Value to send to MARBL as xco2", &
          default=284.317, units="ppm")
    endif
    call get_param(param_file, mdl, "ATM_ALT_CO2_OPT", atm_co2_opt, &
        "Source of alternate atmospheric CO2 [constant, diagnostic, or prognostic]", &
        default="constant")
    select case (trim(atm_co2_opt))
      case("prognostic")
        CS%atm_alt_co2_iopt = atm_co2_prognostic_iopt
      case("diagnostic")
        CS%atm_alt_co2_iopt = atm_co2_diagnostic_iopt
      case("constant")
        CS%atm_alt_co2_iopt = atm_co2_constant_iopt
      case DEFAULT
        write(err_message, "(3A)") "'", trim(atm_co2_opt), "' is not a valid ATM_ALT_CO2_OPT value"
        call MOM_error(FATAL, err_message)
    end select
    if (CS%atm_alt_co2_iopt == atm_co2_constant_iopt) then
      call get_param(param_file, mdl, "ATM_ALT_CO2_CONST", CS%atm_alt_co2_const, &
          "Value to send to MARBL as xco2_alt_co2", &
          default=284.317, units="ppm")
    endif

    ! Register diagnostic fields for outputing forcing values
    ! These fields are posted from convert_driver_fields_to_forcings(), and they are received
    ! in physical units so no conversion is necessary here.
    CS%diag_ids%atm_fine_dust = register_diag_field("ocean_model", "ATM_FINE_DUST_FLUX_CPL", &
        CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
        day, "ATM_FINE_DUST_FLUX from cpl", "kg/m^2/s")
    CS%diag_ids%atm_coarse_dust = register_diag_field("ocean_model", "ATM_COARSE_DUST_FLUX_CPL", &
        CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
        day, "ATM_COARSE_DUST_FLUX from cpl", "kg/m^2/s")
    CS%diag_ids%atm_bc = register_diag_field("ocean_model", "ATM_BLACK_CARBON_FLUX_CPL", &
        CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
        day, "ATM_BLACK_CARBON_FLUX from cpl",  "kg/m^2/s")

    CS%diag_ids%ice_dust = register_diag_field("ocean_model", "SEAICE_DUST_FLUX_CPL", &
        CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
        day, "SEAICE_DUST_FLUX from cpl", "kg/m^2/s")
    CS%diag_ids%ice_bc = register_diag_field("ocean_model", "SEAICE_BLACK_CARBON_FLUX_CPL", &
        CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
        day, "SEAICE_BLACK_CARBON_FLUX from cpl", "kg/m^2/s")

  end procedure MARBL_forcing_init
  module procedure convert_driver_fields_to_forcings
    integer :: i, j, is, ie, js, je, m
    real :: atm_fe_bioavail_frac     !< Fraction of iron from the atmosphere available for biological uptake [1]
    real :: seaice_fe_bioavail_frac  !< Fraction of iron from sea ice available for biological uptake [1]
    real :: iron_flux_conversion     !< Factor to convert iron flux from kg m-2 s-1 -> mmol m-3 (m s-1)
    real :: ndep_conversion          !< Factor to convert nitrogen deposition from kg m-2 s-1 -> mmol m-3 (m s-1)
    if (.not. CS%use_marbl_tracers) return

    is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
    ndep_conversion = (1.e6/14.) * (US%m_to_Z * US%T_to_s)
    iron_flux_conversion = (1.e6 / molw_Fe) * (US%m_to_Z * US%T_to_s)

    ! Post fields from coupler to diagnostics
    ! TODO: units from diag register are incorrect; we should be converting these in the cap, I think
    if (CS%diag_ids%atm_fine_dust > 0) &
      call post_data(CS%diag_ids%atm_fine_dust, atm_fine_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
          CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_coarse_dust > 0) &
      call post_data(CS%diag_ids%atm_coarse_dust, atm_coarse_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
          CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_bc > 0) &
      call post_data(CS%diag_ids%atm_bc, atm_bc_flux(is-i0:ie-i0,js-j0:je-j0), CS%diag, &
          mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_dust > 0) &
      call post_data(CS%diag_ids%ice_dust, seaice_dust_flux(is-i0:ie-i0,js-j0:je-j0), CS%diag, &
          mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_bc > 0) &
      call post_data(CS%diag_ids%ice_bc, seaice_bc_flux(is-i0:ie-i0,js-j0:je-j0), CS%diag, &
          mask=G%mask2dT(is:ie,js:je))

    do j=js,je ; do i=is,ie
      ! Nitrogen Deposition
      fluxes%nhx_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * nhx_dep(i-i0,j-j0)
      fluxes%noy_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * noy_dep(i-i0,j-j0)
    enddo ; enddo

    ! Atmospheric CO2
    select case (CS%atm_co2_iopt)
      case (atm_co2_prognostic_iopt)
        if (associated(atm_co2_prog)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_co2(i,j) = G%mask2dT(i,j) * atm_co2_prog(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_CO2_OPT = 'prognostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_diagnostic_iopt)
        if (associated(atm_co2_diag)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_co2(i,j) = G%mask2dT(i,j) * atm_co2_diag(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_CO2_OPT = 'diagnostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_constant_iopt)
        do j=js,je ; do i=is,ie
          fluxes%atm_co2(i,j) = G%mask2dT(i,j) * CS%atm_co2_const
        enddo ; enddo
    end select

    ! Alternate Atmospheric CO2
    select case (CS%atm_alt_co2_iopt)
      case (atm_co2_prognostic_iopt)
        if (associated(atm_co2_prog)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * atm_co2_prog(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_ALT_CO2_OPT = 'prognostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_diagnostic_iopt)
        if (associated(atm_co2_diag)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * atm_co2_diag(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_ALT_CO2_OPT = 'diagnostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_constant_iopt)
        do j=js,je ; do i=is,ie
          fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * CS%atm_co2_const
        enddo ; enddo
    end select

    ! Dust flux
    if (associated(atm_fine_dust_flux)) then
      do j=js,je ; do i=is,ie
        fluxes%dust_flux(i,j) = (US%kg_m2s_to_RZ_T * G%mask2dT(i,j)) * &
            ((atm_fine_dust_flux(i-i0,j-j0) + atm_coarse_dust_flux(i-i0,j-j0)) + &
            seaice_dust_flux(i-i0,j-j0))
      enddo ; enddo
    endif

    if (associated(atm_bc_flux)) then
      do j=js,je ; do i=is,ie
        ! TODO: abort if atm_fine_dust_flux and atm_coarse_dust_flux are not associated?
        ! Contribution of atmospheric dust to iron flux
        if (atm_coarse_dust_flux(i-i0,j-j0) < &
            CS%dust_ratio_thres * atm_fine_dust_flux(i-i0,j-j0)) then
          atm_fe_bioavail_frac = CS%fe_bioavail_frac_offset + CS%dust_ratio_to_fe_bioavail_frac * &
            (CS%dust_ratio_thres - atm_coarse_dust_flux(i-i0,j-j0) / atm_fine_dust_flux(i-i0,j-j0))
        else
          atm_fe_bioavail_frac = CS%fe_bioavail_frac_offset
        endif

        ! Contribution of atmospheric dust to iron flux
        fluxes%iron_flux(i,j) = (atm_fe_bioavail_frac * &
            (CS%iron_frac_in_atm_fine_dust * atm_fine_dust_flux(i-i0,j-j0) + &
            CS%iron_frac_in_atm_coarse_dust * atm_coarse_dust_flux(i-i0,j-j0)))

        ! Contribution of atmospheric black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (CS%atm_bc_fe_bioavail_frac * &
            (CS%atm_fe_to_bc_ratio * atm_bc_flux(i-i0,j-j0)))

        seaice_fe_bioavail_frac = atm_fe_bioavail_frac
        ! Contribution of seaice dust to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (seaice_fe_bioavail_frac * &
            (CS%iron_frac_in_seaice_dust * seaice_dust_flux(i-i0,j-j0)))

        ! Contribution of seaice black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (CS%seaice_bc_fe_bioavail_frac * &
            (CS%seaice_fe_to_bc_ratio * seaice_bc_flux(i-i0,j-j0)))

        ! Unit conversion (kg m-2 s-1 -> conc Z T-1)
        fluxes%iron_flux(i,j) = (G%mask2dT(i,j) * iron_flux_conversion) * fluxes%iron_flux(i,j)

      enddo ; enddo
    endif

      ! Per ice-category forcings
      ! If the cap receives per-category fields, memory should be allocated in fluxes
    if (associated(ifrac_n)) then
      do j=js,je ; do i=is,ie
        fluxes%fracr_cat(i,j,1) = min(1., afracr(i-i0,j-j0))
        fluxes%qsw_cat(i,j,1) = swnet_afracr(i-i0,j-j0)
        do m=1,size(ifrac_n, 3)
          fluxes%fracr_cat(i,j,m+1) = min(1., ifrac_n(i-i0,j-j0,m))
          fluxes%qsw_cat(i,j,m+1)   = swpen_ifrac_n(i-i0,j-j0,m)
        enddo
        where (fluxes%fracr_cat(i,j,:) > 0.)
          fluxes%qsw_cat(i,j,:) = fluxes%qsw_cat(i,j,:) / fluxes%fracr_cat(i,j,:)
        elsewhere
          fluxes%fracr_cat(i,j,:) = 0.
          fluxes%qsw_cat(i,j,:) = 0.
        endwhere
        fluxes%fracr_cat(i,j,:) = G%mask2dT(i,j) * fluxes%fracr_cat(i,j,:)
        fluxes%qsw_cat(i,j,:)   = (US%W_m2_to_QRZ_T * G%mask2dT(i,j)) * fluxes%qsw_cat(i,j,:)
      enddo ; enddo
    endif

  end procedure convert_driver_fields_to_forcings
end submodule MARBL_forcing_mod_s
