submodule (MOM_surface_forcing_nuopc) MOM_surface_forcing_nuopc_s
#include <MOM_memory.h>
  implicit none
contains
module procedure convert_IOB_to_fluxes
  real, dimension(SZI_(G),SZJ_(G)) :: &
    data_restore,  & !< The surface value toward which to restore [S ~> ppt] or [C ~> degC]
    PmE_adj,       & !< The adjustment to PminusE that will cause the salinity
                     !! to be restored toward its target value [kg/(m^2 * s)]
    net_FW,        & !< The area integrated net freshwater flux into the ocean [R Z L2 T-1 ~> kg s-1]
    net_FW2,       & !< The area averaged net freshwater flux into the ocean [R Z T-1 ~> kg m-2 s-1]
    work_sum,      & !< A 2-d array that is used as the work space for a global
                     !! sum, used with units of [L2 ~> m2] or [R Z L2 T-1 ~> kg s-1]
    open_ocn_mask    !< a binary field indicating where ice is present based on frazil criteria

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  logical :: restore_salinity !< local copy of the argument restore_salt, if it
                              !! is present, or false (no restoring) otherwise.
  logical :: restore_sst      !< local copy of the argument restore_temp, if it
                              !! is present, or false (no restoring) otherwise.
  real :: delta_sss           !< temporary storage for sss diff from restoring value [S ~> ppt]
  real :: delta_sst           !< temporary storage for sst diff from restoring value [C ~> degC]
  real :: kg_m2_s_conversion  !< A combination of unit conversion factors for rescaling
                              !! mass fluxes [R Z s m2 kg-1 T-1 ~> 1].

  real :: C_p                 !< heat capacity of seawater [J kg-1 degC-1]
  real :: sign_for_net_FW_bug !< Should be +1. but an old bug can be recovered by using -1.

  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  kg_m2_s_conversion = US%kg_m2s_to_RZ_T
  C_p                    = US%Q_to_J_kg*US%degC_to_C*fluxes%C_p
  open_ocn_mask(:,:)     = 1.0
  pme_adj(:,:)           = 0.0
  fluxes%vPrecGlobalAdj  = 0.0
  fluxes%vPrecGlobalScl  = 0.0
  fluxes%saltFluxGlobalAdj = 0.0
  fluxes%saltFluxGlobalScl = 0.0
  fluxes%netFWGlobalAdj = 0.0
  fluxes%netFWGlobalScl = 0.0

  restore_salinity = .false.
  if (present(restore_salt)) restore_salinity = restore_salt
  restore_sst = .false.
  if (present(restore_temp)) restore_sst = restore_temp

  ! allocation and initialization if this is the first time that this
  ! flux type has been used.
  if (fluxes%dt_buoy_accum < 0) then
    call allocate_forcing_type(G, fluxes, water=.true., heat=.true., ustar=.true., &
                               press=.true., fix_accum_bug=.not.CS%ustar_gustless_bug, &
                               cfc=CS%use_CFC, marbl=CS%use_marbl_tracers, hevap=CS%enthalpy_cpl, &
                               tau_mag=.true., ice_ncat=IOB%ice_ncat)
    call safe_alloc_ptr(fluxes%omega_w2x,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%p_surf     ,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed)
    if (CS%use_limited_P_SSH) then
      fluxes%p_surf_SSH => fluxes%p_surf
    else
      fluxes%p_surf_SSH => fluxes%p_surf_full
    endif

    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_in,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%BBL_tidal_dis,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed)

    if (CS%allow_flux_adjustments) then
      call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
      call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)
    endif

    do j=js-2,je+2 ; do i=is-2,ie+2
      fluxes%BBL_tidal_dis(i,j)   = US%Z_to_L**2*CS%TKE_tidal(i,j)
      fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
    enddo ; enddo

    if (restore_temp) call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)

  endif   ! endif for allocation and initialization

  if (((associated(IOB%ustar_berg) .and. (.not.associated(fluxes%ustar_berg))) &
    .or. (associated(IOB%area_berg) .and. (.not.associated(fluxes%area_berg)))) &
    .or. (associated(IOB%mass_berg) .and. (.not.associated(fluxes%mass_berg)))) &
    call allocate_forcing_type(G, fluxes, iceberg=.true.)

  if ((.not.coupler_type_initialized(fluxes%tr_fluxes)) .and. &
      coupler_type_initialized(IOB%fluxes)) &
    call coupler_type_spawn(IOB%fluxes, fluxes%tr_fluxes, &
                            (/is,is,ie,ie/), (/js,js,je,je/))
  !   It might prove valuable to use the same array extents as the rest of the
  ! ocean model, rather than using haloless arrays, in which case the last line
  ! would be: (             (/isd,is,ie,ied/), (/jsd,js,je,jed/))

  ! allocation and initialization on first call to this routine
  if (CS%area_surf < 0.0) then
    do j=js,je ; do i=is,ie
      work_sum(i,j) = G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    CS%area_surf = reproducing_sum(work_sum, isr, ier, jsr, jer, unscale=US%L_to_m**2)
  endif    ! endif for allocation and initialization


  ! Indicate that there are new unused fluxes.
  fluxes%fluxes_used = .false.
  fluxes%dt_buoy_accum = valid_time

  if (CS%allow_flux_adjustments) then
    fluxes%heat_added(:,:)=0.0
    fluxes%salt_flux_added(:,:)=0.0
  endif

  do j=js,je ; do i=is,ie
    fluxes%salt_flux(i,j) = 0.0
    fluxes%vprec(i,j) = 0.0
  enddo ; enddo

  ! Salinity restoring logic
  if (restore_salinity) then
    call time_interp_external(CS%srestore_handle, Time, data_restore, scale=US%ppt_to_S)
    ! open_ocn_mask indicates where to restore salinity (1 means restore, 0 does not)
    open_ocn_mask(:,:) = 1.0
    if (CS%mask_srestore_under_ice) then ! Do not restore under sea-ice
      do j=js,je ; do i=is,ie
        if (sfc_state%SST(i,j) <= -0.0539*US%degC_to_C*US%S_to_ppt*sfc_state%SSS(i,j)) open_ocn_mask(i,j)=0.0
      enddo ; enddo
    endif
    if (CS%salt_restore_as_sflux) then
      do j=js,je ; do i=is,ie
        delta_sss = data_restore(i,j) - sfc_state%SSS(i,j)
        delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
        fluxes%salt_flux(i,j) = 1.e-3*US%S_to_ppt*G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)* &
                  (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j)) *delta_sss  ! kg Salt m-2 s-1
      enddo ; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%salt_flux, G, fluxes%saltFluxGlobalScl, &
                         unit_scale=US%RZ_T_to_kg_m2s)
          fluxes%saltFluxGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) =  G%areaT(is:ie,js:je) * fluxes%salt_flux(is:ie,js:je)
          fluxes%saltFluxGlobalAdj = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer, unscale=US%RZL2_to_kg*US%s_to_T) &
                                     / CS%area_surf
          fluxes%salt_flux(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) -  fluxes%saltFluxGlobalAdj
        endif
      endif
      fluxes%salt_flux_added(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) ! Diagnostic
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          delta_sss = sfc_state%SSS(i,j) - data_restore(i,j)
          delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
          fluxes%vprec(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j))* &
                      (CS%Rho0*CS%Flux_const) * &
                      delta_sss / (0.5*(sfc_state%SSS(i,j) + data_restore(i,j)))
        endif
      enddo ; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%vprec, G, fluxes%vPrecGlobalScl, &
                                        unit_scale=US%RZ_T_to_kg_m2s)
          fluxes%vPrecGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je) * fluxes%vprec(is:ie,js:je)
          fluxes%vPrecGlobalAdj = reproducing_sum(work_sum(:,:), isr, ier, jsr, jer, unscale=US%RZL2_to_kg*US%s_to_T) &
                                  / CS%area_surf
          do j=js,je ; do i=is,ie
            fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%vPrecGlobalAdj ) * G%mask2dT(i,j)
          enddo ; enddo
        endif
      endif
    endif
  endif

  ! SST restoring logic
  if (restore_sst) then
    call time_interp_external(CS%trestore_handle, Time, data_restore, scale=US%degC_to_C)
    do j=js,je ; do i=is,ie
      delta_sst = data_restore(i,j) - sfc_state%SST(i,j)
      delta_sst = sign(1.0,delta_sst)*min(abs(delta_sst),CS%max_delta_trestore)
      fluxes%heat_added(i,j) = G%mask2dT(i,j) * CS%trestore_mask(i,j) * &
                      (CS%Rho0*fluxes%C_p) * delta_sst * CS%Flux_const   ! Q R Z T-1 ~> W m-2
    enddo ; enddo
  endif


  ! Check that liquid runoff has a place to go
  if (CS%liquid_runoff_from_data .and. .not. associated(IOB%lrunoff)) then
    call MOM_error(FATAL, "liquid runoff is being added via data_override but "// &
                          "there is no associated runoff in the IOB")
    return
  endif
  if (associated(IOB%lrunoff)) then
    if (CS%liquid_runoff_from_data) call data_override('OCN', 'runoff', IOB%lrunoff, Time)
  endif

  ! obtain fluxes from IOB; note the staggering of indices
  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie

    if (associated(IOB%lprec)) &
      fluxes%lprec(i,j) =  kg_m2_s_conversion * IOB%lprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%fprec)) &
      fluxes%fprec(i,j) = kg_m2_s_conversion * IOB%fprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%q_flux)) &
      fluxes%evap(i,j) = kg_m2_s_conversion * IOB%q_flux(i-i0,j-j0) * G%mask2dT(i,j)

    ! liquid runoff flux
    if (associated(IOB%lrunoff)) then
      fluxes%lrunoff(i,j) = kg_m2_s_conversion * IOB%lrunoff(i-i0,j-j0) * G%mask2dT(i,j)
    endif

    ! ice runoff flux
    if (associated(IOB%frunoff)) then
      fluxes%frunoff(i,j) = kg_m2_s_conversion * IOB%frunoff(i-i0,j-j0) * G%mask2dT(i,j)
    endif

    ! add liquid glc runoff flux via rof
    if (associated(IOB%lrunoff_glc)) then
      fluxes%lrunoff_glc(i,j) = kg_m2_s_conversion * IOB%lrunoff_glc(i-i0,j-j0) * G%mask2dT(i,j)
    endif

    ! ice glc runoff flux via rof
    if (associated(IOB%frunoff_glc)) then
      fluxes%frunoff_glc(i,j) = kg_m2_s_conversion * IOB%frunoff_glc(i-i0,j-j0) * G%mask2dT(i,j)
    endif

    if (associated(IOB%ustar_berg)) &
      fluxes%ustar_berg(i,j) = US%m_to_Z*US%T_to_s * IOB%ustar_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%area_berg)) &
      fluxes%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%mass_berg)) &
      fluxes%mass_berg(i,j) = US%m_to_Z*US%kg_m3_to_R * IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%lw_flux)) &
      fluxes%lw(i,j) = US%W_m2_to_QRZ_T * IOB%lw_flux(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%t_flux)) &
      fluxes%sens(i,j) = US%W_m2_to_QRZ_T * IOB%t_flux(i-i0,j-j0) * G%mask2dT(i,j)

    ! sea ice and snow melt heat flux [Q R Z T-1 ~> W m-2]
    if (associated(IOB%seaice_melt_heat)) &
      fluxes%seaice_melt_heat(i,j) = US%W_m2_to_QRZ_T * G%mask2dT(i,j) * IOB%seaice_melt_heat(i-i0,j-j0)

    ! water flux due to sea ice and snow melt [km m-2 s-1]
    if (associated(IOB%seaice_melt)) &
      fluxes%seaice_melt(i,j) = kg_m2_s_conversion * G%mask2dT(i,j) * IOB%seaice_melt(i-i0,j-j0)

    fluxes%latent(i,j) = 0.0
    ! notice minus sign since fprec is positive into the ocean
    if (associated(IOB%fprec)) then
      fluxes%latent(i,j)            = fluxes%latent(i,j) - &
          IOB%fprec(i-i0,j-j0)*US%W_m2_to_QRZ_T*CS%latent_heat_fusion
      fluxes%latent_fprec_diag(i,j) = - G%mask2dT(i,j) * IOB%fprec(i-i0,j-j0)*US%W_m2_to_QRZ_T*CS%latent_heat_fusion
    endif
    ! notice minus sign since frunoff is positive into the ocean
    if (associated(IOB%frunoff)) then
      fluxes%latent(i,j)              = fluxes%latent(i,j) - &
          IOB%frunoff(i-i0,j-j0) * US%W_m2_to_QRZ_T * CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j) = - G%mask2dT(i,j) * &
          IOB%frunoff(i-i0,j-j0) * US%W_m2_to_QRZ_T * CS%latent_heat_fusion
    endif
    ! notice minus sign since frunoff_glc is positive into the ocean
    if (associated(IOB%frunoff_glc)) then
      fluxes%latent(i,j)              = fluxes%latent(i,j) - &
          IOB%frunoff_glc(i-i0,j-j0) * US%W_m2_to_QRZ_T * CS%latent_heat_fusion
      fluxes%latent_frunoff_glc_diag(i,j) = fluxes%latent_frunoff_glc_diag(i,j) - G%mask2dT(i,j) * &
          IOB%frunoff_glc(i-i0,j-j0) * US%W_m2_to_QRZ_T * CS%latent_heat_fusion
    endif
    if (associated(IOB%q_flux)) then
      fluxes%latent(i,j)           = fluxes%latent(i,j) + &
          IOB%q_flux(i-i0,j-j0)*US%W_m2_to_QRZ_T*CS%latent_heat_vapor
      fluxes%latent_evap_diag(i,j) = G%mask2dT(i,j) * IOB%q_flux(i-i0,j-j0)*US%W_m2_to_QRZ_T*CS%latent_heat_vapor
    endif
    fluxes%latent(i,j) = G%mask2dT(i,j) * fluxes%latent(i,j)

    if (associated(IOB%sw_flux_vis_dir)) &
      fluxes%sw_vis_dir(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_vis_dir(i-i0,j-j0)

    if (associated(IOB%sw_flux_vis_dif)) &
      fluxes%sw_vis_dif(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_vis_dif(i-i0,j-j0)

    if (associated(IOB%sw_flux_nir_dir)) &
      fluxes%sw_nir_dir(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_nir_dir(i-i0,j-j0)

    if (associated(IOB%sw_flux_nir_dif)) &
      fluxes%sw_nir_dif(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_nir_dif(i-i0,j-j0)

    fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
                     fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)

    ! enthalpy terms
    if (CS%enthalpy_cpl) then
      if (associated(IOB%hrofl)) &
        fluxes%heat_content_lrunoff(i,j) = US%W_m2_to_QRZ_T * IOB%hrofl(i-i0,j-j0) * G%mask2dT(i,j)

      if (associated(IOB%hrofi)) &
        fluxes%heat_content_frunoff(i,j) = US%W_m2_to_QRZ_T * IOB%hrofi(i-i0,j-j0) * G%mask2dT(i,j)

      if (associated(IOB%hrain)) &
        fluxes%heat_content_lprec(i,j)   = US%W_m2_to_QRZ_T * IOB%hrain(i-i0,j-j0) * G%mask2dT(i,j)

      if (associated(IOB%hsnow)) &
        fluxes%heat_content_fprec(i,j)   = US%W_m2_to_QRZ_T * IOB%hsnow(i-i0,j-j0) * G%mask2dT(i,j)

       if (associated(IOB%hevap)) &
        fluxes%heat_content_evap(i,j)    = US%W_m2_to_QRZ_T * IOB%hevap(i-i0,j-j0) * G%mask2dT(i,j)

       if (associated(IOB%hcond)) &
        fluxes%heat_content_cond(i,j)    = US%W_m2_to_QRZ_T * IOB%hcond(i-i0,j-j0) * G%mask2dT(i,j)

       if (associated(IOB%hrofl_glc)) &
        fluxes%heat_content_lrunoff_glc(i,j) = US%W_m2_to_QRZ_T * IOB%hrofl_glc(i-i0,j-j0) * G%mask2dT(i,j)

       if (associated(IOB%hrofi_glc)) &
        fluxes%heat_content_frunoff_glc(i,j) = US%W_m2_to_QRZ_T * IOB%hrofi_glc(i-i0,j-j0) * G%mask2dT(i,j)
    endif

    ! sea ice fraction [1]
    if (associated(IOB%ice_fraction) .and. associated(fluxes%ice_fraction)) &
      fluxes%ice_fraction(i,j) = G%mask2dT(i,j) * IOB%ice_fraction(i-i0,j-j0)
    ! 10-m wind speed squared [m2 s-2]
    if (associated(IOB%u10_sqr) .and. associated(fluxes%u10_sqr)) &
      fluxes%u10_sqr(i,j) = US%m_to_L**2 * US%T_to_s**2 * G%mask2dT(i,j) * IOB%u10_sqr(i-i0,j-j0)

  enddo ; enddo

  ! Copy MARBL-specific IOB fields into fluxes; also set some MARBL-specific forcings to other values
  ! (constants, values from netCDF, etc)
  if (CS%use_marbl_tracers) &
    call convert_driver_fields_to_forcings(IOB%atm_fine_dust_flux, IOB%atm_coarse_dust_flux, &
                                           IOB%seaice_dust_flux, IOB%atm_bc_flux, IOB%seaice_bc_flux, &
                                           IOB%nhx_dep, IOB%noy_dep, IOB%atm_co2_prog, IOB%atm_co2_diag, &
                                           IOB%afracr, IOB%swnet_afracr, IOB%ifrac_n, IOB%swpen_ifrac_n, &
                                           Time, G, US, i0, j0, fluxes, CS%marbl_forcing_CSp)

  ! wave to ocean coupling
  if ( associated(IOB%lamult)) then
    do j=js,je ; do i=is,ie
      if (IOB%ice_fraction(i-i0,j-j0) <= 0.05 ) then
        fluxes%lamult(i,j) = IOB%lamult(i-i0,j-j0)
      else
        fluxes%lamult(i,j) = 1.0
      endif
    enddo ; enddo
    call pass_var(fluxes%lamult, G%domain, halo=1 )
  endif

  ! applied surface pressure from atmosphere and cryosphere
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * US%kg_m3_to_R*US%m_s_to_L_T**2*IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * US%kg_m3_to_R*US%m_s_to_L_T**2*IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
      enddo ; enddo
    endif
    fluxes%accumulate_p_surf = .true. ! Multiple components may contribute to surface pressure.
  endif

  if (associated(IOB%salt_flux)) then
    do j=js,je ; do i=is,ie
      fluxes%salt_flux(i,j)    = G%mask2dT(i,j)*(fluxes%salt_flux(i,j) + kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0))
      fluxes%salt_flux_in(i,j) = G%mask2dT(i,j)*( kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0) )
    enddo ; enddo
  endif

  ! adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    sign_for_net_FW_bug = 1.
    if (CS%use_net_FW_adjustment_sign_bug) sign_for_net_FW_bug = -1.
    do j=js,je ; do i=is,ie
      net_FW(i,j) = ((((fluxes%lprec(i,j) + fluxes%fprec(i,j)) + fluxes%seaice_melt(i,j)) + &
                      ((fluxes%lrunoff(i,j) + fluxes%frunoff(i,j)) + (fluxes%lrunoff_glc(i,j) + &
                        fluxes%frunoff_glc(i,j)))) + (fluxes%evap(i,j) + fluxes%vprec(i,j)) ) * &
                      G%areaT(i,j)
      net_FW2(i,j) = net_FW(i,j) / G%areaT(i,j)
    enddo ; enddo

    if (CS%adjust_net_fresh_water_by_scaling) then
      call adjust_area_mean_to_zero(net_FW2, G, fluxes%netFWGlobalScl, unscale=US%RZ_T_to_kg_m2s)
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = fluxes%vprec(i,j) + (net_FW2(i,j) - net_FW(i,j) / G%areaT(i,j)) * G%mask2dT(i,j)
      enddo ; enddo
    else
      fluxes%netFWGlobalAdj = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer, unscale=US%RZL2_to_kg*US%s_to_T) / &
                              CS%area_surf
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%netFWGlobalAdj ) * G%mask2dT(i,j)
      enddo ; enddo
    endif

  endif

  if (coupler_type_initialized(fluxes%tr_fluxes) .and. &
      coupler_type_initialized(IOB%fluxes)) &
    call coupler_type_copy_data(IOB%fluxes, fluxes%tr_fluxes)

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to fluxes
    call apply_flux_adjustments(G, US, CS, Time, fluxes)
  endif

  ! Allow for user-written code to alter fluxes after all the above
  call user_alter_forcing(sfc_state, fluxes, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)

end procedure convert_IOB_to_fluxes
module procedure convert_IOB_to_forces
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    taux_at_q, & !< Zonal wind stresses at q points [R Z L T-2 ~> Pa]
    tauy_at_q    !< Meridional wind stresses at q points [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rigidity_at_h, & !< Ice rigidity at tracer points [L4 Z-1 T-1 ~> m3 s-1]
    taux_at_h, & !< Zonal wind stresses at h points [R Z L T-2 ~> Pa]
    tauy_at_h    !< Meridional wind stresses at h points [R Z L T-2 ~> Pa]
  real :: gustiness     !< unresolved gustiness that contributes to ustar [R Z L T-2 ~> Pa]
  real :: Irho0         !< inverse of the mean density in [Z L-1 R-1 ~> m3 kg-1]
  real :: taux2, tauy2  !< squared wind stresses [R2 Z2 L2 T-4 ~> Pa2]
  real :: tau_mag       !< magnitude of the wind stress [R Z L T-2 ~> Pa]
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal wind stress units [R Z L T-2 Pa-1 ~> 1]
  real :: stress_conversion ! A unit conversion factor from Pa times any stress multiplier [R Z L T-2 Pa-1 ~> 1]
  real :: I_GEarth      !< The inverse of the gravitational acceleration [T2 Z L-2 ~> s2 m-1]
  real :: Kv_rho_ice    !< (CS%Kv_sea_ice / CS%density_sea_ice) [L4 Z-2 T-1 R-1 ~> m5 s-1 kg-1]
  real :: mass_ice      !< mass of sea ice at a face [R Z ~> kg m-2]
  real :: mass_eff      !< effective mass of sea ice for rigidity [R Z ~> kg m-2]
  integer :: wind_stagger  !< AGRID, BGRID_NE, or CGRID_NE (integers from MOM_domains)
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0, istk
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1
  i0 = is - isc_bnd ; j0 = js - jsc_bnd

  Irho0 = US%L_to_Z / CS%Rho0
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z
  stress_conversion = Pa_conversion * CS%wind_stress_multiplier

  ! allocation and initialization if this is the first time that this
  ! mechanical forcing type has been used.
  if (.not.forces%initialized) then

    call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., press=.true., tau_mag=.true.)

    call safe_alloc_ptr(forces%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(forces%p_surf_full,isd,ied,jsd,jed)
    call safe_alloc_ptr(forces%omega_w2x,isd,ied,jsd,jed)

    if (CS%rigid_sea_ice) then
      call safe_alloc_ptr(forces%rigidity_ice_u,IsdB,IedB,jsd,jed)
      call safe_alloc_ptr(forces%rigidity_ice_v,isd,ied,JsdB,JedB)
    endif

    forces%initialized = .true.
  endif

  if ( (associated(IOB%area_berg) .and. (.not. associated(forces%area_berg))) .or. &
       (associated(IOB%mass_berg) .and. (.not. associated(forces%mass_berg))) ) &
    call allocate_mech_forcing(G, forces, iceberg=.true.)

  if (associated(IOB%ice_rigidity)) then
    rigidity_at_h(:,:) = 0.0
    call safe_alloc_ptr(forces%rigidity_ice_u,IsdB,IedB,jsd,jed)
    call safe_alloc_ptr(forces%rigidity_ice_v,isd,ied,JsdB,JedB)
  endif

  forces%accumulate_rigidity = .true. ! Multiple components may contribute to rigidity.
  if (associated(forces%rigidity_ice_u)) forces%rigidity_ice_u(:,:) = 0.0
  if (associated(forces%rigidity_ice_v)) forces%rigidity_ice_v(:,:) = 0.0

  if ( associated(IOB%lamult) ) then
    call allocate_mech_forcing(G, forces, waves=.true., num_stk_bands=0)
  elseif ( associated(IOB%ustkb) ) then
    call allocate_mech_forcing(G, forces, waves=.true., num_stk_bands=IOB%num_stk_bands)
  endif

  ! applied surface pressure from atmosphere and cryosphere
  if (CS%use_limited_P_SSH) then
    forces%p_surf_SSH => forces%p_surf
  else
    forces%p_surf_SSH => forces%p_surf_full
  endif
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * US%kg_m3_to_R*US%m_s_to_L_T**2*IOB%p(i-i0,j-j0)
        forces%p_surf(i,j) = MIN(forces%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * US%kg_m3_to_R*US%m_s_to_L_T**2*IOB%p(i-i0,j-j0)
        forces%p_surf(i,j) = forces%p_surf_full(i,j)
      enddo ; enddo
    endif
  else
    do j=js,je ; do i=is,ie
      forces%p_surf_full(i,j) = 0.0
      forces%p_surf(i,j) = 0.0
    enddo ; enddo
  endif
  forces%accumulate_p_surf = .true. ! Multiple components may contribute to surface pressure.

#ifdef CESMCOUPLED
  wind_stagger = AGRID
#else
  wind_stagger = CS%wind_stagger
#endif

  if (wind_stagger == BGRID_NE) then
    ! This is necessary to fill in the halo points.
    taux_at_q(:,:) = 0.0 ; tauy_at_q(:,:) = 0.0
  endif
  if (wind_stagger == AGRID) then
    ! This is necessary to fill in the halo points.
    taux_at_h(:,:) = 0.0 ; tauy_at_h(:,:) = 0.0
  endif

  ! obtain fluxes from IOB; note the staggering of indices
  do j=js,je ; do i=is,ie
    if (associated(IOB%area_berg)) &
      forces%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%mass_berg)) &
      forces%mass_berg(i,j) = US%m_to_Z*US%kg_m3_to_R * IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (associated(IOB%ice_rigidity)) &
      rigidity_at_h(i,j) = US%m_to_L**3*US%Z_to_L*US%T_to_s * IOB%ice_rigidity(i-i0,j-j0) * G%mask2dT(i,j)

    if (wind_stagger == BGRID_NE) then
      if (associated(IOB%u_flux)) taux_at_q(I,J) = IOB%u_flux(i-i0,j-j0) * stress_conversion
      if (associated(IOB%v_flux)) tauy_at_q(I,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
    elseif (wind_stagger == AGRID) then
      if (associated(IOB%u_flux)) taux_at_h(i,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
      if (associated(IOB%v_flux)) tauy_at_h(i,j) = IOB%v_flux(i-i0,j-j0) * stress_conversion
    else ! C-grid wind stresses.
      if (associated(IOB%u_flux)) forces%taux(I,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
      if (associated(IOB%v_flux)) forces%tauy(i,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
    endif

  enddo ; enddo

  ! surface momentum stress related fields as function of staggering
  if (wind_stagger == BGRID_NE) then
    if (G%symmetric) &
      call fill_symmetric_edges(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE)
    call pass_vector(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE, halo=1)

    do j=js,je ; do I=Isq,Ieq
      forces%taux(I,j) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0.0) &
        forces%taux(I,j) = (G%mask2dBu(I,J)*taux_at_q(I,J) + &
                            G%mask2dBu(I,J-1)*taux_at_q(I,J-1)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      forces%tauy(i,J) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0.0) &
        forces%tauy(i,J) = (G%mask2dBu(I,J)*tauy_at_q(I,J) + &
                            G%mask2dBu(I-1,J)*tauy_at_q(I-1,J)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
    enddo ; enddo

    ! ustar is required for the bulk mixed layer formulation. The background value
    ! of 0.02 Pa is a relatively small value intended to give reasonable behavior
    ! in regions of very weak winds.

    do j=js,je ; do i=is,ie
      tau_mag = 0.0 ; gustiness = CS%gust_const
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
           (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0.0) then
        tau_mag = sqrt(((G%mask2dBu(I,J)*((taux_at_q(I,J)**2) + (tauy_at_q(I,J)**2)) + &
            G%mask2dBu(I-1,J-1)*((taux_at_q(I-1,J-1)**2) + (tauy_at_q(I-1,J-1)**2))) + &
           (G%mask2dBu(I,J-1)*((taux_at_q(I,J-1)**2) + (tauy_at_q(I,J-1)**2)) + &
            G%mask2dBu(I-1,J)*((taux_at_q(I-1,J)**2) + (tauy_at_q(I-1,J)**2))) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)
      endif
      forces%tau_mag(i,j) = US%L_to_Z*(gustiness + tau_mag)
      forces%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0*tau_mag)
    enddo ; enddo
    call pass_vector(forces%taux, forces%tauy, G%Domain, halo=1)
  elseif (wind_stagger == AGRID) then
    call pass_vector(taux_at_h, tauy_at_h, G%Domain, To_All+Omit_Corners, stagger=AGRID, halo=1)

    do j=js,je ; do I=Isq,Ieq
      forces%taux(I,j) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i+1,j)) > 0.0) &
        forces%taux(I,j) = (G%mask2dT(i,j)*taux_at_h(i,j) + &
                            G%mask2dT(i+1,j)*taux_at_h(i+1,j)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i+1,j))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      forces%tauy(i,J) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i,j+1)) > 0.0) &
        forces%tauy(i,J) = (G%mask2dT(i,j)*tauy_at_h(i,j) + &
                            G%mask2dT(i,J+1)*tauy_at_h(i,j+1)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i,j+1))
    enddo ; enddo

    do j=js,je ; do i=is,ie
      gustiness = CS%gust_const
      if (CS%read_gust_2d .and. (G%mask2dT(i,j) > 0.0)) gustiness = CS%gust(i,j)
      forces%tau_mag(i,j) = US%L_to_Z*(gustiness + G%mask2dT(i,j) * sqrt((taux_at_h(i,j)**2) + (tauy_at_h(i,j)**2)))
      forces%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0 * G%mask2dT(i,j) * &
                               sqrt((taux_at_h(i,j)**2) + (tauy_at_h(i,j)**2)))
      forces%omega_w2x(i,j) = atan2(tauy_at_h(i,j), taux_at_h(i,j))
    enddo ; enddo
    call pass_vector(forces%taux, forces%tauy, G%Domain, halo=1)
  else ! C-grid wind stresses.
    if (G%symmetric) &
      call fill_symmetric_edges(forces%taux, forces%tauy, G%Domain)
    call pass_vector(forces%taux, forces%tauy, G%Domain, halo=1)

    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0.0) &
        taux2 = (G%mask2dCu(I-1,j)*(forces%taux(I-1,j)**2) + &
                 G%mask2dCu(I,j)*(forces%taux(I,j)**2)) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))

      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0.0) &
        tauy2 = (G%mask2dCv(i,J-1)*(forces%tauy(i,J-1)**2) + &
                 G%mask2dCv(i,J)*(forces%tauy(i,J)**2)) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))

      if (CS%read_gust_2d) then
        forces%tau_mag(i,j) = US%L_to_Z*(CS%gust(i,j) + sqrt(taux2 + tauy2))
        forces%ustar(i,j) = sqrt(CS%gust(i,j)*Irho0 + Irho0*sqrt(taux2 + tauy2))
      else
        forces%tau_mag(i,j) = US%L_to_Z*(CS%gust_const + sqrt(taux2 + tauy2))
        forces%ustar(i,j) = sqrt(CS%gust_const*Irho0 + Irho0*sqrt(taux2 + tauy2))
      endif
    enddo ; enddo

  endif   ! endif for wind related fields

  ! wave to ocean coupling
  if ( associated(IOB%ustkb) ) then

    forces%stk_wavenumbers(:) = IOB%stk_wavenumbers * US%Z_to_m
    do istk = 1,IOB%num_stk_bands
      do j=js,je ; do i=is,ie
        forces%ustkb(i,j,istk) = IOB%ustkb(i-I0,j-J0,istk) * US%m_s_to_L_T
        forces%vstkb(i,j,istk) = IOB%vstkb(i-I0,j-J0,istk) * US%m_s_to_L_T
      enddo ; enddo
      call pass_var(forces%ustkb(:,:,istk), G%domain )
      call pass_var(forces%vstkb(:,:,istk), G%domain )
    enddo
  endif

  ! sea ice related dynamic fields
  if (associated(IOB%ice_rigidity)) then
    call pass_var(rigidity_at_h, G%Domain, halo=1)
    do I=is-1,ie ; do j=js,je
      forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + &
              min(rigidity_at_h(i,j), rigidity_at_h(i+1,j))
    enddo ; enddo
    do i=is,ie ; do J=js-1,je
      forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + &
              min(rigidity_at_h(i,j), rigidity_at_h(i,j+1))
    enddo ; enddo
  endif

  if (CS%rigid_sea_ice) then
    call pass_var(forces%p_surf_full, G%Domain, halo=1)
    I_GEarth = 1.0 / CS%g_Earth
    Kv_rho_ice = (CS%Kv_sea_ice / CS%density_sea_ice)
    do I=is-1,ie ; do j=js,je
      mass_ice = min(forces%p_surf_full(i,j), forces%p_surf_full(i+1,j)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass)**2 / (mass_ice + CS%rigid_sea_ice_mass)
      endif
      forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + Kv_rho_ice * mass_eff
    enddo ; enddo
    do i=is,ie ; do J=js-1,je
      mass_ice = min(forces%p_surf_full(i,j), forces%p_surf_full(i,j+1)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass)**2 / (mass_ice + CS%rigid_sea_ice_mass)
      endif
      forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + Kv_rho_ice * mass_eff
    enddo ; enddo
  endif

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to forces
    call apply_force_adjustments(G, US, CS, Time, forces)
  endif

!###  ! Allow for user-written code to alter fluxes after all the above
!###  call user_alter_mech_forcing(forces, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end procedure convert_IOB_to_forces
module procedure apply_flux_adjustments
  real, dimension(SZI_(G),SZJ_(G)) :: temp_at_h !< Fluxes at h points [W m-2 or kg m-2 s-1]
  integer :: isc, iec, jsc, jec, i, j
  logical :: overrode_h
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  overrode_h = .false.
  call data_override('OCN', 'hflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%heat_added(i,j) = fluxes%heat_added(i,j) + US%W_m2_to_QRZ_T*temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%heat_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'sflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%salt_flux_added(i,j) = fluxes%salt_flux_added(i,j) + &
        US%kg_m2s_to_RZ_T * temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%salt_flux_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'prcme_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%vprec(i,j) = fluxes%vprec(i,j) + US%kg_m2s_to_RZ_T * temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%vprec, G%Domain)
end procedure apply_flux_adjustments
module procedure apply_force_adjustments
  real, dimension(SZI_(G),SZJ_(G)) :: tempx_at_h !< Delta to zonal wind stress at h points [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: tempy_at_h !< Delta to meridional wind stress at h points [R Z L T-2 ~> Pa]
  integer :: isc, iec, jsc, jec, i, j
  real :: dLonDx, dLonDy, rDlon, cosA, sinA, zonal_tau, merid_tau
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal units [R Z L T-2 Pa-1 ~> 1]
  logical :: overrode_x, overrode_y
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

  tempx_at_h(:,:) = 0.0 ; tempy_at_h(:,:) = 0.0
  ! Either reads data or leaves contents unchanged
  overrode_x = .false. ; overrode_y = .false.
  call data_override('OCN', 'taux_adj', tempx_at_h(isc:iec,jsc:jec), Time, override=overrode_x)
  call data_override('OCN', 'tauy_adj', tempy_at_h(isc:iec,jsc:jec), Time, override=overrode_y)

  if (overrode_x .or. overrode_y) then
    if (.not. (overrode_x .and. overrode_y)) call MOM_error(FATAL,"apply_flux_adjustments: "//&
            "Both taux_adj and tauy_adj must be specified, or neither, in data_table")

    ! Rotate winds
    call pass_vector(tempx_at_h, tempy_at_h, G%Domain, To_All, AGRID, halo=1)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      dLonDx = G%geoLonCu(I,j) - G%geoLonCu(I-1,j)
      dLonDy = G%geoLonCv(i,J) - G%geoLonCv(i,J-1)
      rDlon = sqrt( dLonDx * dLonDx + dLonDy * dLonDy )
      if (rDlon > 0.) rDlon = 1. / rDlon
      cosA = dLonDx * rDlon
      sinA = dLonDy * rDlon
      zonal_tau = Pa_conversion * tempx_at_h(i,j)
      merid_tau = Pa_conversion * tempy_at_h(i,j)
      tempx_at_h(i,j) = cosA * zonal_tau - sinA * merid_tau
      tempy_at_h(i,j) = sinA * zonal_tau + cosA * merid_tau
    enddo ; enddo

    ! Average to C-grid locations
    do j=jsc,jec ; do I=isc-1,iec
      forces%taux(I,j) = forces%taux(I,j) + 0.5 * ( tempx_at_h(i,j) + tempx_at_h(i+1,j) )
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      forces%tauy(i,J) = forces%tauy(i,J) + 0.5 * ( tempy_at_h(i,j) + tempy_at_h(i,j+1) )
    enddo ; enddo
  endif ! overrode_x .or. overrode_y

end procedure apply_force_adjustments
module procedure forcing_save_restart
  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return
  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end procedure forcing_save_restart
module procedure surface_forcing_init
  real :: utide  ! The RMS tidal velocity [Z T-1 ~> m s-1].
  type(directories)  :: dirs
  logical            :: new_sim, iceberg_flux_diags, glc_runoff_diags, fix_ustar_gustless_bug
  logical :: test_value  ! This is used to determine whether a logical parameter is being set explicitly.
  logical :: explicit_bug, explicit_fix ! These indicate which parameters are set explicitly.
  type(time_type)    :: Time_frc
  character(len=200) :: TideAmp_file, gust_file, salt_file, temp_file ! Input file names.
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_surface_forcing_nuopc"  ! This module's name.
  character(len=48)  :: stagger
  character(len=48)  :: flnam
  character(len=240) :: basin_file
  integer :: i, j, isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  id_clock_forcing=cpu_clock_id('Ocean surface forcing', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_forcing)

  CS%diag => diag

  call write_version_number(version)
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv)
  call get_param(param_file, mdl, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the "//&
                 "atmosphere and floating sea-ice or ice shelves. This is "//&
                 "needed because the FMS coupling structure does not "//&
                 "limit the water that can be frozen out of the ocean and "//&
                 "the ice-ocean heat fluxes are treated explicitly.  No "//&
                 "limit is applied if a negative value is used.", &
                 units="Pa", default=-1.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_TO_ZERO", &
                 CS%adjust_net_srestore_to_zero, &
                 "If true, adjusts the salinity restoring seen to zero "//&
                 "whether restoring is via a salt flux or virtual precip.",&
                 default=restore_salt)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_BY_SCALING", &
                 CS%adjust_net_srestore_by_scaling, &
                 "If true, adjustments to salt restoring to achieve zero net are "//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_TO_ZERO", &
                 CS%adjust_net_fresh_water_to_zero, &
                 "If true, adjusts the net fresh-water forcing seen "//&
                 "by the ocean (including restoring) to zero.", default=.false.)
  if (CS%adjust_net_fresh_water_to_zero) &
    call get_param(param_file, mdl, "USE_NET_FW_ADJUSTMENT_SIGN_BUG", &
                 CS%use_net_FW_adjustment_sign_bug, &
                   "If true, use the wrong sign for the adjustment to "//&
                   "the net fresh-water.", default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_BY_SCALING", &
                 CS%adjust_net_fresh_water_by_scaling, &
                 "If true, adjustments to net fresh water to achieve zero net are "//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ICE_SALT_CONCENTRATION", &
                 CS%ice_salt_concentration, &
                 "The assumed sea-ice salinity needed to reverse engineer the "//&
                 "melt flux (or ice-ocean fresh-water flux).", &
                 units="kg/kg", default=0.005)
  call get_param(param_file, mdl, "USE_LIMITED_PATM_SSH", CS%use_limited_P_SSH, &
                 "If true, return the sea surface height with the "//&
                 "correction for the atmospheric (and sea-ice) pressure "//&
                 "limited by max_p_surf instead of the full atmospheric "//&
                 "pressure.", default=.true.)

  call get_param(param_file, mdl, "WIND_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the "//&
                 "staggering of the input wind stress field.  Valid "//&
                 "values are 'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; CS%wind_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; CS%wind_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; CS%wind_stagger = CGRID_NE
  else ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif
  call get_param(param_file, mdl, "WIND_STRESS_MULTIPLIER", CS%wind_stress_multiplier, &
                 "A factor multiplying the wind-stress given to the ocean by the "//&
                 "coupler. This is used for testing and should be =1.0 for any "//&
                 "production runs.", units="nondim", default=1.0)

  call get_param(param_file, mdl, "USE_CFC_CAP", CS%use_CFC, &
                 default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "USE_MARBL_TRACERS", CS%use_marbl_tracers, &
                 default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "ENTHALPY_FROM_COUPLER", CS%enthalpy_cpl, &
                 "If True, the heat (enthalpy) associated with mass entering/leaving the "//&
                 "ocean is provided via coupler.", default=.false.)

  if (restore_salt) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0)
    call get_param(param_file, mdl, "SALT_RESTORE_FILE", CS%salt_restore_file, &
                 "A file in which to find the surface salinity to use for restoring.", &
                 default="salt_restore.nc")
    call get_param(param_file, mdl, "SALT_RESTORE_VARIABLE", CS%salt_restore_var_name, &
                 "The name of the surface salinity variable to read from "//&
                 "SALT_RESTORE_FILE for restoring salinity.", &
                 default="salt")
    call get_param(param_file, mdl, "SRESTORE_AS_SFLUX", CS%salt_restore_as_sflux, &
                 "If true, the restoring of salinity is applied as a salt "//&
                 "flux instead of as a freshwater flux.", default=.false.)
    call get_param(param_file, mdl, "MAX_DELTA_SRESTORE", CS%max_delta_srestore, &
                 "The maximum salinity difference used in restoring terms.", &
                 units="PSU or g kg-1", default=999.0, scale=US%ppt_to_S)
    call get_param(param_file, mdl, "MASK_SRESTORE_UNDER_ICE", &
                 CS%mask_srestore_under_ice, &
                 "If true, disables SSS restoring under sea-ice based on a frazil "//&
                 "criteria (SST<=Tf). Only used when RESTORE_SALINITY is True.",      &
                 default=.false.)
    call get_param(param_file, mdl, "MASK_SRESTORE_MARGINAL_SEAS", &
                 CS%mask_srestore_marginal_seas, &
                 "If true, disable SSS restoring in marginal seas. Only used when "//&
                 "RESTORE_SALINITY is True.", default=.false.)
    call get_param(param_file, mdl, "BASIN_FILE", basin_file, &
                 "A file in which to find the basin masks, in variable 'basin'.", &
                 default="basin.nc")
    basin_file = trim(CS%inputdir) // trim(basin_file)
    call safe_alloc_ptr(CS%basin_mask,isd,ied,jsd,jed) ; CS%basin_mask(:,:) = 1.0
    if (CS%mask_srestore_marginal_seas) then
      call MOM_read_data(basin_file,'basin',CS%basin_mask,G%domain, timelevel=1)
      do j=jsd,jed ; do i=isd,ied
        if (CS%basin_mask(i,j) >= 6.0) then ; CS%basin_mask(i,j) = 0.0
        else ; CS%basin_mask(i,j) = 1.0 ; endif
      enddo ; enddo
    endif
    call get_param(param_file, mdl, "MASK_SRESTORE", CS%mask_srestore, &
                 "If true, read a file (salt_restore_mask) containing "//&
                 "a mask for SSS restoring.", default=.false.)
  endif

  if (restore_temp) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0)
    call get_param(param_file, mdl, "SST_RESTORE_FILE", CS%temp_restore_file, &
                 "A file in which to find the surface temperature to use for restoring.", &
                 default="temp_restore.nc")
    call get_param(param_file, mdl, "SST_RESTORE_VARIABLE", CS%temp_restore_var_name, &
                 "The name of the surface temperature variable to read from "//&
                 "SST_RESTORE_FILE for restoring sst.", &
                 default="temp")

    call get_param(param_file, mdl, "MAX_DELTA_TRESTORE", CS%max_delta_trestore, &
                 "The maximum sst difference used in restoring terms.", &
                 units="degC ", default=999.0, scale=US%degC_to_C)
    call get_param(param_file, mdl, "MASK_TRESTORE", CS%mask_trestore, &
                 "If true, read a file (temp_restore_mask) containing "//&
                 "a mask for SST restoring.", default=.false.)

  endif

! Optionally read tidal amplitude from input file (m s-1) on model grid.
! Otherwise use default tidal amplitude for bottom frictionally-generated
! dissipation. Default cd_tides is chosen to yield approx 1 TWatt of
! work done against tides globally using OSU tidal amplitude.
  call get_param(param_file, mdl, "CD_TIDES", CS%cd_tides, &
                 "The drag coefficient that applies to the tides.", &
                 units="nondim", default=1.0e-4)
  call get_param(param_file, mdl, "READ_TIDEAMP", CS%read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing "//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (CS%read_TIDEAMP) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying "//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", &
                 default="tideamp.nc")
    CS%utide=0.0
  else
    call get_param(param_file, mdl, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
  endif

  call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed)

  if (CS%read_TIDEAMP) then
    TideAmp_file = trim(CS%inputdir) // trim(TideAmp_file)
    call MOM_read_data(TideAmp_file,'tideamp',CS%TKE_tidal,G%domain,timelevel=1, scale=US%m_to_Z*US%T_to_s)
    do j=jsd,jed ; do i=isd,ied
      utide = CS%TKE_tidal(i,j)
      CS%TKE_tidal(i,j) = G%mask2dT(i,j)*CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  else
    do j=jsd,jed ; do i=isd,ied
      utide = CS%utide
      CS%TKE_tidal(i,j) = CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  endif

  ! initialize time interpolator module
  call time_interp_external_init()

! Optionally read a x-y gustiness field in place of a global
! constant.

  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
               "The background gustiness in the winds.", &
               units="Pa", default=0.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.")

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed)
    gust_file = trim(CS%inputdir) // trim(gust_file)
    call MOM_read_data(gust_file, 'gustiness', CS%gust, G%domain, timelevel=1, &
               scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z) ! units in file should be Pa
  endif

  call get_param(param_file, mdl, "USTAR_GUSTLESS_BUG", CS%ustar_gustless_bug, &
                 "If true include a bug in the time-averaging of the gustless wind friction velocity", &
                 default=.false., do_not_log=.true.)
  ! This is used to test whether USTAR_GUSTLESS_BUG is being actively set.
  call get_param(param_file, mdl, "USTAR_GUSTLESS_BUG", test_value, default=.true., do_not_log=.true.)
  explicit_bug = CS%ustar_gustless_bug .eqv. test_value
  call get_param(param_file, mdl, "FIX_USTAR_GUSTLESS_BUG", fix_ustar_gustless_bug, &
                 "If true correct a bug in the time-averaging of the gustless wind friction velocity", &
                 default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "FIX_USTAR_GUSTLESS_BUG", test_value, default=.false., do_not_log=.true.)
  explicit_fix = fix_ustar_gustless_bug .eqv. test_value

  if (explicit_bug .and. explicit_fix .and. (fix_ustar_gustless_bug .eqv. CS%ustar_gustless_bug)) then
    ! USTAR_GUSTLESS_BUG is being explicitly set, and should not be changed.
    call MOM_error(FATAL, "USTAR_GUSTLESS_BUG and FIX_USTAR_GUSTLESS_BUG are both being set "//&
                   "with inconsistent values.  FIX_USTAR_GUSTLESS_BUG is an obsolete "//&
                   "parameter and should be removed.")
  elseif (explicit_fix) then
    call MOM_error(WARNING, "FIX_USTAR_GUSTLESS_BUG is an obsolete parameter.  "//&
                   "Use USTAR_GUSTLESS_BUG instead (noting that it has the opposite sense).")
    CS%ustar_gustless_bug = .not.fix_ustar_gustless_bug
  endif
  call log_param(param_file, mdl, "USTAR_GUSTLESS_BUG", CS%ustar_gustless_bug, &
                 "If true include a bug in the time-averaging of the gustless wind friction velocity", &
                 default=.false.)

! See whether sufficiently thick sea ice should be treated as rigid.
  call get_param(param_file, mdl, "USE_RIGID_SEA_ICE", CS%rigid_sea_ice, &
                 "If true, sea-ice is rigid enough to exert a "//&
                 "nonhydrostatic pressure that resist vertical motion.", &
                 default=.false.)
  if (CS%rigid_sea_ice) then
    call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%Z_to_m*US%m_s_to_L_T**2)
    call get_param(param_file, mdl, "SEA_ICE_MEAN_DENSITY", CS%density_sea_ice, &
                 "A typical density of sea ice, used with the kinematic "//&
                 "viscosity, when USE_RIGID_SEA_ICE is true.", &
                 units="kg m-3", default=900.0, scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "SEA_ICE_VISCOSITY", CS%Kv_sea_ice, &
                 "The kinematic viscosity of sufficiently thick sea ice "//&
                 "for use in calculating the rigidity of sea ice.", &
                 units="m2 s-1", default=1.0e9, scale=US%Z_to_L**2*US%m_to_L**2*US%T_to_s)
    call get_param(param_file, mdl, "SEA_ICE_RIGID_MASS", CS%rigid_sea_ice_mass, &
                 "The mass of sea-ice per unit area at which the sea-ice "//&
                 "starts to exhibit rigidity", &
                 units="kg m-2", default=1000.0, scale=US%kg_m3_to_R*US%m_to_Z)
  endif

  call get_param(param_file, mdl, "ALLOW_ICEBERG_FLUX_DIAGNOSTICS", iceberg_flux_diags, &
                 "If true, makes available diagnostics of fluxes from icebergs "//&
                 "as seen by MOM6.", default=.false.)

  call get_param(param_file, mdl, "ALLOW_GLC_RUNOFF_DIAGNOSTICS", glc_runoff_diags, &
                 "If true, makes available diagnostics of separate glacier runoff fluxes"//&
                 "as seen by MOM6.", default=.false.)

  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles, &
                                   use_berg_fluxes=iceberg_flux_diags, use_waves=use_waves, &
                                   use_cfcs=CS%use_CFC, use_glc_runoff=glc_runoff_diags)

  call get_param(param_file, mdl, "ALLOW_FLUX_ADJUSTMENTS", CS%allow_flux_adjustments, &
                 "If true, allows flux adjustments to specified via the "//&
                 "data_table using the component name 'OCN'.", default=.false.)

  call get_param(param_file, mdl, "LIQUID_RUNOFF_FROM_DATA", CS%liquid_runoff_from_data, &
                 "If true, allows liquid river runoff to be specified via the "//&
                 "data_table using the component name 'OCN'.", default=.false.)

  if (CS%allow_flux_adjustments .or. CS%liquid_runoff_from_data) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
  endif

  ! Set up MARBL forcing control structure
  call MARBL_forcing_init(G, US, param_file, diag, Time, CS%inputdir, CS%use_marbl_tracers, &
      CS%marbl_forcing_CSp)

  if (present(restore_salt)) then ; if (restore_salt) then
    salt_file = trim(CS%inputdir) // trim(CS%salt_restore_file)
    CS%srestore_handle = init_external_field(salt_file, CS%salt_restore_var_name, domain=G%Domain%mpp_domain)
    call safe_alloc_ptr(CS%srestore_mask,isd,ied,jsd,jed) ; CS%srestore_mask(:,:) = 1.0
    if (CS%mask_srestore) then ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'salt_restore_mask.nc'
      call MOM_read_data(flnam,'mask', CS%srestore_mask, G%domain, timelevel=1)
    endif
  endif ; endif

  if (present(restore_temp)) then ; if (restore_temp) then
    temp_file = trim(CS%inputdir) // trim(CS%temp_restore_file)
    CS%trestore_handle = init_external_field(temp_file, CS%temp_restore_var_name, domain=G%Domain%mpp_domain)
    call safe_alloc_ptr(CS%trestore_mask,isd,ied,jsd,jed) ; CS%trestore_mask(:,:) = 1.0
    if (CS%mask_trestore) then  ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'temp_restore_mask.nc'
      call MOM_read_data(flnam, 'mask', CS%trestore_mask, G%domain, timelevel=1)
    endif
  endif ; endif

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
  call restart_init_end(CS%restart_CSp)

  if (associated(CS%restart_CSp)) then
    call Get_MOM_Input(dirs=dirs)

    new_sim = .false.
    if ((dirs%input_filename(1:1) == 'n') .and. &
        (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.
    if (.not.new_sim) then
      call restore_state(dirs%input_filename, dirs%restart_input_dir, Time_frc, &
                         G, CS%restart_CSp)
    endif
  endif

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end procedure surface_forcing_init
module procedure surface_forcing_end
  if (present(fluxes)) call deallocate_forcing_type(fluxes)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end procedure surface_forcing_end
module procedure ice_ocn_bnd_type_chksum
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  integer :: outunit ! The output unit to write to
  outunit = stdout
  root = is_root_pe()

  if (root) write(outunit,*) "BEGIN CHECKSUM(ice_ocean_boundary_type):: ", id, timestep
  chks = field_chksum( iobt%u_flux         ) ; if (root) write(outunit,100) 'iobt%u_flux          ', chks
  chks = field_chksum( iobt%v_flux         ) ; if (root) write(outunit,100) 'iobt%v_flux          ', chks
  chks = field_chksum( iobt%t_flux         ) ; if (root) write(outunit,100) 'iobt%t_flux          ', chks
  chks = field_chksum( iobt%q_flux         ) ; if (root) write(outunit,100) 'iobt%q_flux          ', chks
  chks = field_chksum( iobt%seaice_melt_heat) ; if (root) write(outunit,100) 'iobt%seaice_melt_heat', chks
  chks = field_chksum( iobt%seaice_melt)     ; if (root) write(outunit,100) 'iobt%seaice_melt    ', chks
  chks = field_chksum( iobt%salt_flux      ) ; if (root) write(outunit,100) 'iobt%salt_flux      ', chks
  chks = field_chksum( iobt%lw_flux        ) ; if (root) write(outunit,100) 'iobt%lw_flux        ', chks
  chks = field_chksum( iobt%sw_flux_vis_dir) ; if (root) write(outunit,100) 'iobt%sw_flux_vis_dir', chks
  chks = field_chksum( iobt%sw_flux_vis_dif) ; if (root) write(outunit,100) 'iobt%sw_flux_vis_dif', chks
  chks = field_chksum( iobt%sw_flux_nir_dir) ; if (root) write(outunit,100) 'iobt%sw_flux_nir_dir', chks
  chks = field_chksum( iobt%sw_flux_nir_dif) ; if (root) write(outunit,100) 'iobt%sw_flux_nir_dif', chks
  chks = field_chksum( iobt%lprec          ) ; if (root) write(outunit,100) 'iobt%lprec          ', chks
  chks = field_chksum( iobt%fprec          ) ; if (root) write(outunit,100) 'iobt%fprec          ', chks
  chks = field_chksum( iobt%lrunoff        ) ; if (root) write(outunit,100) 'iobt%lrunoff        ', chks
  chks = field_chksum( iobt%frunoff        ) ; if (root) write(outunit,100) 'iobt%frunoff        ', chks
  chks = field_chksum( iobt%lrunoff_glc    ) ; if (root) write(outunit,100) 'iobt%lrunoff_glc    ', chks
  chks = field_chksum( iobt%frunoff_glc    ) ; if (root) write(outunit,100) 'iobt%frunoff_glc    ', chks
  chks = field_chksum( iobt%p              ) ; if (root) write(outunit,100) 'iobt%p              ', chks
  if (associated(iobt%ice_fraction)) then
    chks = field_chksum( iobt%ice_fraction ) ; if (root) write(outunit,100) 'iobt%ice_fraction   ', chks
  endif
  if (associated(iobt%u10_sqr)) then
    chks = field_chksum( iobt%u10_sqr ) ; if (root) write(outunit,100) 'iobt%u10_sqr   ', chks
  endif
  if (associated(iobt%ustar_berg)) then
    chks = field_chksum( iobt%ustar_berg ) ; if (root) write(outunit,100) 'iobt%ustar_berg     ', chks
  endif
  if (associated(iobt%area_berg)) then
    chks = field_chksum( iobt%area_berg  ) ; if (root) write(outunit,100) 'iobt%area_berg      ', chks
  endif
  if (associated(iobt%mass_berg)) then
    chks = field_chksum( iobt%mass_berg  ) ; if (root) write(outunit,100) 'iobt%mass_berg      ', chks
  endif

  ! MARBL forcing
  if (associated(iobt%atm_fine_dust_flux)) then
    chks = field_chksum(iobt%atm_fine_dust_flux)
    if (root) write(outunit,110) 'iobt%atm_fine_dust_flux   ', chks
  endif
  if (associated(iobt%atm_coarse_dust_flux)) then
    chks = field_chksum(iobt%atm_coarse_dust_flux)
    if (root) write(outunit,110) 'iobt%atm_coarse_dust_flux ', chks
  endif
  if (associated(iobt%seaice_dust_flux)) then
    chks = field_chksum(iobt%seaice_dust_flux)
    if (root) write(outunit,110) 'iobt%seaice_dust_flux     ', chks
  endif
  if (associated(iobt%atm_bc_flux)) then
    chks = field_chksum(iobt%atm_bc_flux)
    if (root) write(outunit,110) 'iobt%atm_bc_flux          ', chks
  endif
  if (associated(iobt%seaice_bc_flux)) then
    chks = field_chksum(iobt%seaice_bc_flux)
    if (root) write(outunit,110) 'iobt%seaice_bc_flux       ', chks
  endif
  if (associated(iobt%nhx_dep)) then
    chks = field_chksum(iobt%nhx_dep)
    if (root) write(outunit,100) 'iobt%nhx_dep    ', chks
  endif
  if (associated(iobt%noy_dep)) then
    chks = field_chksum(iobt%noy_dep)
    if (root) write(outunit,100) 'iobt%noy_dep    ', chks
  endif
  if (associated(iobt%atm_co2_prog)) then
    chks = field_chksum(iobt%atm_co2_prog)
    if (root) write(outunit,110) 'iobt%atm_co2_prog         ', chks
  endif
  if (associated(iobt%atm_co2_diag)) then
    chks = field_chksum(iobt%atm_co2_diag)
    if (root) write(outunit,110) 'iobt%atm_co2_diag         ', chks
  endif
  if (associated(iobt%afracr)) then
    chks = field_chksum(iobt%afracr)
    if (root) write(outunit,100) 'iobt%afracr     ', chks
  endif
  if (associated(iobt%swnet_afracr)) then
    chks = field_chksum(iobt%swnet_afracr)
    if (root) write(outunit,110) 'iobt%swnet_afracr         ', chks
  endif
  if (associated(iobt%ifrac_n)) then
    chks = field_chksum(iobt%ifrac_n)
    if (root) write(outunit,100) 'iobt%ifrac_n    ', chks
  endif
  if (associated(iobt%swpen_ifrac_n)) then
    chks = field_chksum(iobt%swpen_ifrac_n)
    if (root) write(outunit,110) 'iobt%swpen_ifrac_n        ', chks
  endif

  ! enthalpy
  if (associated(iobt%hrofl)) then
    chks = field_chksum( iobt%hrofl  ) ; if (root) write(outunit,100) 'iobt%hrofl      ', chks
  endif
  if (associated(iobt%hrofi)) then
    chks = field_chksum( iobt%hrofi  ) ; if (root) write(outunit,100) 'iobt%hrofi      ', chks
  endif
  if (associated(iobt%hrain)) then
    chks = field_chksum( iobt%hrain  ) ; if (root) write(outunit,100) 'iobt%hrain      ', chks
  endif
  if (associated(iobt%hsnow)) then
    chks = field_chksum( iobt%hsnow  ) ; if (root) write(outunit,100) 'iobt%hsnow      ', chks
  endif
  if (associated(iobt%hevap)) then
    chks = field_chksum( iobt%hevap  ) ; if (root) write(outunit,100) 'iobt%hevap      ', chks
  endif
  if (associated(iobt%hcond)) then
    chks = field_chksum( iobt%hcond  ) ; if (root) write(outunit,100) 'iobt%hcond      ', chks
  endif
  if (associated(iobt%hrofl_glc)) then
    chks = field_chksum( iobt%hrofl_glc  ) ; if (root) write(outunit,100) 'iobt%hrofl_glc      ', chks
  endif
  if (associated(iobt%hrofl_glc)) then
    chks = field_chksum( iobt%hrofl_glc  ) ; if (root) write(outunit,100) 'iobt%hrofl_glc      ', chks
  endif

100 FORMAT("   CHECKSUM::",A20," = ",Z20)
110 FORMAT("   CHECKSUM::",A30," = ",Z20)

  call coupler_type_write_chksums(iobt%fluxes, outunit, 'iobt%')

end procedure ice_ocn_bnd_type_chksum
end submodule MOM_surface_forcing_nuopc_s
