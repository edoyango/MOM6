submodule (MOM_surface_forcing_gfdl) MOM_surface_forcing_gfdl_s
#include <MOM_memory.h>
  implicit none
contains
module procedure convert_IOB_to_fluxes
  real, dimension(SZI_(G),SZJ_(G)) :: &
    data_restore,  & ! The surface value toward which to restore [S ~> ppt] or [C ~> degC]
    SST_anom,      & ! Instantaneous sea surface temperature anomalies from a target value [C ~> degC]
    SSS_anom,      & ! Instantaneous sea surface salinity anomalies from a target value [S ~> ppt]
    SSS_mean,      & ! A (mean?) salinity about which to normalize local salinity
                     ! anomalies when calculating restorative precipitation anomalies [S ~> ppt]
    net_FW,        & ! The area integrated net freshwater flux into the ocean [R Z L2 T-1 ~> kg s-1]
    net_FW2,       & ! The area averaged net freshwater flux into the ocean [R Z T-1 ~> kg m-2 s-1]
    work_sum,      & ! A 2-d array that is used as the work space for global sums [L2 ~> m2] or [R Z L2 T-1 ~> kg s-1]
    open_ocn_mask    ! a binary field indicating where ice is present based on frazil criteria [nondim]

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  real :: delta_sss           ! temporary storage for sss diff from restoring value [S ~> ppt]
  real :: delta_sst           ! temporary storage for sst diff from restoring value [C ~> degC]

  real :: kg_m2_s_conversion  ! A combination of unit conversion factors for rescaling
                              ! mass fluxes [R Z s m2 kg-1 T-1 ~> 1]
  real :: rhoXcp              ! Reference density times heat capacity times unit scaling
                              ! factors [Q R C-1 ~> J m-3 degC-1]
  real :: sign_for_net_FW_bug ! Should be +1. but an old bug can be recovered by using -1 [nondim]

  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  kg_m2_s_conversion = US%kg_m2s_to_RZ_T
  if (CS%restore_temp) rhoXcp = CS%rho_restore * fluxes%C_p
  open_ocn_mask(:,:)     = 1.0
  fluxes%vPrecGlobalAdj  = 0.0
  fluxes%vPrecGlobalScl  = 0.0
  fluxes%saltFluxGlobalAdj = 0.0
  fluxes%saltFluxGlobalScl = 0.0
  fluxes%netFWGlobalAdj = 0.0
  fluxes%netFWGlobalScl = 0.0

  ! allocation and initialization if this is the first time that this
  ! flux type has been used.
  if (fluxes%dt_buoy_accum < 0) then
    call allocate_forcing_type(G, fluxes, water=.true., heat=.true., ustar=.not.CS%nonBous, press=.true., &
                               fix_accum_bug=.not.CS%ustar_gustless_bug, tau_mag=CS%nonBous,&
                               carbon=CS%allow_carbon_flux_exchange)

    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed)
    if (CS%use_limited_P_SSH) then
      fluxes%p_surf_SSH => fluxes%p_surf
    else
      fluxes%p_surf_SSH => fluxes%p_surf_full
    endif

    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_in,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%BBL_tidal_dis,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)

    if (associated(IOB%excess_salt)) call safe_alloc_ptr(fluxes%salt_left_behind,isd,ied,jsd,jed)

    do j=js-2,je+2 ; do i=is-2,ie+2
      fluxes%BBL_tidal_dis(i,j)   = CS%BBL_tidal_dis(i,j)
      fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
    enddo ; enddo


  endif   ! endif for allocation and initialization


  if (((associated(IOB%ustar_berg) .and. (.not.associated(fluxes%ustar_berg))) &
    .or. (associated(IOB%area_berg) .and. (.not.associated(fluxes%area_berg)))) &
    .or. (associated(IOB%mass_berg) .and. (.not.associated(fluxes%mass_berg)))) &
    call allocate_forcing_type(G, fluxes, iceberg=.true.)

  if ((.not.coupler_type_initialized(fluxes%tr_fluxes)) .and. &
      coupler_type_initialized(IOB%fluxes)) &
    call coupler_type_spawn(IOB%fluxes, fluxes%tr_fluxes, (/is,is,ie,ie/), (/js,js,je,je/))
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

  fluxes%heat_added(:,:) = 0.0
  fluxes%salt_flux_added(:,:) = 0.0

  do j=js,je ; do i=is,ie
    fluxes%salt_flux(i,j) = 0.0
    fluxes%vprec(i,j) = 0.0
  enddo ; enddo

  ! Salinity restoring logic
  if (CS%restore_salt) then
    call time_interp_external(CS%srestore_handle, Time, data_restore, scale=US%ppt_to_S)
    if (sfc_state%S_is_absS .and. CS%salt_restore_is_practical) then
      !Adjust the salt restoring data to absolute
      do j=js,je
        do i=is,ie
          data_restore(i,j) = gsw_sr_from_sp(data_restore(i,j))
        enddo
      enddo
    endif
    ! open_ocn_mask indicates where to restore salinity (1 means restore, 0 does not)
    open_ocn_mask(:,:) = 1.0
    if (CS%mask_srestore_under_ice) then ! Do not restore under sea-ice
      do j=js,je ; do i=is,ie
        if (sfc_state%SST(i,j) <= CS%SPEAR_dTf_dS*sfc_state%SSS(i,j)) open_ocn_mask(i,j)=0.0
      enddo ; enddo
    endif
    if (CS%salt_restore_as_sflux) then
      do j=js,je ; do i=is,ie
        delta_sss = data_restore(i,j) - sfc_state%SSS(i,j)
        delta_sss = sign(1.0,delta_sss) * min(abs(delta_sss), CS%max_delta_srestore)
        fluxes%salt_flux(i,j) = 1.e-3*US%S_to_ppt*G%mask2dT(i,j) * (CS%rho_restore*CS%Flux_const_salt)* &
             (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j)) * delta_sss  ! R Z T-1 ~> kg Salt m-2 s-1
      enddo ; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%salt_flux, G, fluxes%saltFluxGlobalScl, &
                          unit_scale=US%RZ_T_to_kg_m2s)
          fluxes%saltFluxGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je)*fluxes%salt_flux(is:ie,js:je) * G%mask2dT(is:ie,js:je)
          fluxes%saltFluxGlobalAdj = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer, unscale=US%RZL2_to_kg*US%s_to_T) &
                                     / CS%area_surf
          fluxes%salt_flux(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) - &
                                          fluxes%saltFluxGlobalAdj * G%mask2dT(is:ie,js:je)
        endif
      endif
      fluxes%salt_flux_added(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) ! Diagnostic
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          delta_sss = sfc_state%SSS(i,j) - data_restore(i,j)
          delta_sss = sign(1.0,delta_sss) * min(abs(delta_sss), CS%max_delta_srestore)
          fluxes%vprec(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j))* &
                      (CS%rho_restore*CS%Flux_const_salt) * &
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
  if (CS%restore_temp) then
    call time_interp_external(CS%trestore_handle, Time, data_restore, scale=US%degC_to_C)
    if ( CS%trestore_SPEAR_ECDA ) then
      do j=js,je ; do i=is,ie
        if (abs(data_restore(i,j)+1.8*US%degC_to_C) < 0.0001*US%degC_to_C) then
          data_restore(i,j) = CS%SPEAR_dTf_dS*sfc_state%SSS(i,j)
        endif
      enddo ; enddo
    endif

    do j=js,je ; do i=is,ie
      delta_sst = data_restore(i,j) - sfc_state%SST(i,j)
      delta_sst = sign(1.0,delta_sst) * min(abs(delta_sst), CS%max_delta_trestore)
      fluxes%heat_added(i,j) = G%mask2dT(i,j) * CS%trestore_mask(i,j) * &
                               rhoXcp * delta_sst * CS%Flux_const_temp  ! [Q R Z T-1 ~> W m-2]
    enddo ; enddo
  endif


  ! obtain fluxes from IOB; note the staggering of indices
  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie

    if (associated(IOB%lprec)) then
      fluxes%lprec(i,j) = kg_m2_s_conversion * IOB%lprec(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%lprec(i-i0,j-j0), G%mask2dT(i,j), i, j, 'lprec', G)
    endif

    if (associated(IOB%fprec)) then
      fluxes%fprec(i,j) = kg_m2_s_conversion * IOB%fprec(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%fprec(i-i0,j-j0), G%mask2dT(i,j), i, j, 'fprec', G)
    endif

    if (associated(IOB%q_flux)) then
      fluxes%evap(i,j) = - kg_m2_s_conversion * IOB%q_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%q_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'q_flux', G)
    endif

    if (associated(IOB%runoff)) then
      fluxes%lrunoff(i,j) = kg_m2_s_conversion * IOB%runoff(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%runoff(i-i0,j-j0), G%mask2dT(i,j), i, j, 'runoff', G)
    endif

    if (associated(IOB%calving)) then
      fluxes%frunoff(i,j) = kg_m2_s_conversion * IOB%calving(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%calving(i-i0,j-j0), G%mask2dT(i,j), i, j, 'calving', G)
    endif

    if (associated(IOB%shelf_sfc_mass_flux)) then
       fluxes%shelf_sfc_mass_flux(i,j) = kg_m2_s_conversion * IOB%shelf_sfc_mass_flux(i-i0,j-j0)
    endif

    if (associated(IOB%ustar_berg)) then
      fluxes%ustar_berg(i,j) = US%m_to_Z*US%T_to_s * IOB%ustar_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%ustar_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'ustar_berg', G)
    endif

    if (associated(IOB%area_berg)) then
      fluxes%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%area_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'area_berg', G)
    endif

    if (associated(IOB%mass_berg)) then
      fluxes%mass_berg(i,j) = US%m_to_Z*US%kg_m3_to_R * IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%mass_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'mass_berg', G)
    endif

    if (associated(IOB%runoff_hflx)) then
      fluxes%heat_content_lrunoff(i,j) = US%W_m2_to_QRZ_T * IOB%runoff_hflx(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%runoff_hflx(i-i0,j-j0), G%mask2dT(i,j), i, j, 'runoff_hflx', G)
    endif

    if (associated(IOB%runoff_carbon) .and. CS%allow_carbon_flux_exchange) then
      fluxes%carbon_content_lrunoff(i,j) = US%kg_m2s_to_RZ_T * IOB%runoff_carbon(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%runoff_carbon(i-i0,j-j0), G%mask2dT(i,j), i, j, 'runoff_carbon', G)
    endif

    if (associated(IOB%calving_hflx)) then
      fluxes%heat_content_frunoff(i,j) = US%W_m2_to_QRZ_T * IOB%calving_hflx(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%calving_hflx(i-i0,j-j0), G%mask2dT(i,j), i, j, 'calving_hflx', G)
    endif

    if (associated(IOB%lw_flux)) then
      fluxes%LW(i,j) = US%W_m2_to_QRZ_T * IOB%lw_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%lw_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'lw_flux', G)
    endif

    if (associated(IOB%t_flux)) then
      fluxes%sens(i,j) = -US%W_m2_to_QRZ_T* IOB%t_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%t_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 't_flux', G)
    endif

    fluxes%latent(i,j) = 0.0
    if (associated(IOB%fprec)) then
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%fprec(i-i0,j-j0)*kg_m2_s_conversion * CS%latent_heat_fusion
      fluxes%latent_fprec_diag(i,j) = -G%mask2dT(i,j) * IOB%fprec(i-i0,j-j0)*kg_m2_s_conversion * CS%latent_heat_fusion
    endif
    if (associated(IOB%calving)) then
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%calving(i-i0,j-j0)*kg_m2_s_conversion * CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j) = -G%mask2dT(i,j) * IOB%calving(i-i0,j-j0)*kg_m2_s_conversion * &
                                        CS%latent_heat_fusion
    endif
    if (associated(IOB%q_flux)) then
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%q_flux(i-i0,j-j0)*kg_m2_s_conversion * CS%latent_heat_vapor
      fluxes%latent_evap_diag(i,j) = -G%mask2dT(i,j) * IOB%q_flux(i-i0,j-j0)*kg_m2_s_conversion * CS%latent_heat_vapor
    endif

    fluxes%latent(i,j) = G%mask2dT(i,j) * fluxes%latent(i,j)

    if (associated(IOB%sw_flux_vis_dir)) then
      fluxes%sw_vis_dir(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_vis_dir(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_vis_dir(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_vis_dir', G)
    endif
    if (associated(IOB%sw_flux_vis_dif)) then
      fluxes%sw_vis_dif(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_vis_dif(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_vis_dif(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_vis_dif', G)
    endif
    if (associated(IOB%sw_flux_nir_dir)) then
      fluxes%sw_nir_dir(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_nir_dir(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_nir_dir(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_nir_dir', G)
    endif
    if (associated(IOB%sw_flux_nir_dif)) then
      fluxes%sw_nir_dif(i,j) = G%mask2dT(i,j) * US%W_m2_to_QRZ_T * IOB%sw_flux_nir_dif(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_nir_dif(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_nir_dif', G)
    endif
    if (CS%answer_date < 20190101) then
      fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
                       fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)
    else
      fluxes%sw(i,j) = (fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j)) + &
                       (fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j))
    endif

  enddo ; enddo

  ! applied surface pressure from atmosphere and cryosphere
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * US%Pa_to_RL2_T2*IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
        if (CS%check_no_land_fluxes) &
          call check_mask_val_consistency(IOB%p(i-i0,j-j0), G%mask2dT(i,j), i, j, 'p', G)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * US%Pa_to_RL2_T2*IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
        if (CS%check_no_land_fluxes) &
          call check_mask_val_consistency(IOB%p(i-i0,j-j0), G%mask2dT(i,j), i, j, 'p', G)
      enddo ; enddo
    endif
    fluxes%accumulate_p_surf = .true. ! Multiple components may contribute to surface pressure.
  endif

  ! more salt restoring logic
  if (associated(IOB%salt_flux)) then
    do j=js,je ; do i=is,ie
      fluxes%salt_flux(i,j)    = G%mask2dT(i,j)*(fluxes%salt_flux(i,j) - kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0))
      fluxes%salt_flux_in(i,j) = G%mask2dT(i,j)*( -kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0) )
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%salt_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'salt_flux', G)
    enddo ; enddo
  endif
  if (associated(IOB%excess_salt)) then
    do j=js,je ; do i=is,ie
      fluxes%salt_left_behind(i,j) = G%mask2dT(i,j)*(kg_m2_s_conversion*IOB%excess_salt(i-i0,j-j0))
    enddo ; enddo
  endif

!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!#CTRL#                           fluxes%vprec, day, valid_time, G, US, CS%ctrl_forcing_CSp)
!#CTRL# endif

  ! adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    sign_for_net_FW_bug = 1.
    if (CS%use_net_FW_adjustment_sign_bug) sign_for_net_FW_bug = -1.
    do j=js,je ; do i=is,ie
      net_FW(i,j) =  (((fluxes%lprec(i,j)   + fluxes%fprec(i,j)) + &
                       (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j))) + &
                       (fluxes%evap(i,j)    + fluxes%vprec(i,j)) ) * G%areaT(i,j)
      !   The following contribution appears to be calculating the volume flux of sea-ice
      ! melt. This calculation is clearly WRONG if either sea-ice has variable
      ! salinity or the sea-ice is completely fresh.
      !   Bob thinks this is trying ensure the net fresh-water of the ocean + sea-ice system
      ! is constant.
      !   To do this correctly we will need a sea-ice melt field added to IOB. -AJA
      if (associated(IOB%salt_flux) .and. (CS%ice_salt_concentration>0.0)) &
        net_FW(i,j) = net_FW(i,j) + sign_for_net_FW_bug * G%areaT(i,j) * &
                     (kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0) / CS%ice_salt_concentration)
      net_FW2(i,j) = net_FW(i,j) / G%areaT(i,j)
    enddo ; enddo

    if (CS%adjust_net_fresh_water_by_scaling) then
      call adjust_area_mean_to_zero(net_FW2, G, fluxes%netFWGlobalScl, unscale=US%RZ_T_to_kg_m2s)
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = fluxes%vprec(i,j) + &
            (net_FW2(i,j) - net_FW(i,j)/G%areaT(i,j)) * G%mask2dT(i,j)
      enddo ; enddo
    else
      fluxes%netFWGlobalAdj = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer, unscale=US%RZL2_to_kg*US%s_to_T) / &
                               CS%area_surf
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%netFWGlobalAdj ) * G%mask2dT(i,j)
      enddo ; enddo
    endif

  endif

  ! Set the wind stresses and ustar.
  if (associated(fluxes%ustar) .and. associated(fluxes%ustar_gustless) .and. associated(fluxes%tau_mag) &
      .and. associated(fluxes%tau_mag_gustless) ) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=fluxes%ustar, &
                              mag_tau=fluxes%tau_mag, gustless_ustar=fluxes%ustar_gustless, &
                              gustless_mag_tau=fluxes%tau_mag_gustless)
  else
    if (associated(fluxes%ustar)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=fluxes%ustar)
    if (associated(fluxes%ustar_gustless)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, gustless_ustar=fluxes%ustar_gustless)
    if (associated(fluxes%tau_mag)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, mag_tau=fluxes%tau_mag)
    if (associated(fluxes%tau_mag_gustless)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, gustless_mag_tau=fluxes%tau_mag_gustless)
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
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rigidity_at_h, &  ! Ice rigidity at tracer points [L4 Z-1 T-1 ~> m3 s-1]
    net_mass_src, &   ! A temporary of net mass sources [R Z T-1 ~> kg m-2 s-1].
    ustar_tmp, &      ! A temporary array of ustar values [Z T-1 ~> m s-1].
    tau_mag_tmp       ! A temporary array of surface stress magnitudes [R Z2 T-2 ~> Pa]
  real :: I_GEarth      ! The inverse of the gravitational acceleration [T2 Z L-2 ~> s2 m-1]
  real :: Kv_rho_ice    ! (CS%Kv_sea_ice / CS%density_sea_ice) [L4 Z-2 T-1 R-1 ~> m5 s-1 kg-1]
  real :: mass_ice      ! mass of sea ice at a face [R Z ~> kg m-2]
  real :: mass_eff      ! effective mass of sea ice for rigidity [R Z ~> kg m-2]
  real :: wt1, wt2      ! Relative weights of previous and current values of ustar [nondim].
  real :: kg_m2_s_conversion  ! A combination of unit conversion factors for rescaling
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
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

  kg_m2_s_conversion = US%kg_m2s_to_RZ_T

  ! allocation and initialization if this is the first time that this
  ! mechanical forcing type has been used.
  if (.not.forces%initialized) then
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.not.CS%nonBous, &
                               press=.true., tau_mag=CS%nonBous)

    call safe_alloc_ptr(forces%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(forces%p_surf_full,isd,ied,jsd,jed)
    if (CS%use_limited_P_SSH) then
      forces%p_surf_SSH => forces%p_surf
    else
      forces%p_surf_SSH => forces%p_surf_full
    endif

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

  ! Set the weights for forcing fields that use running time averages.
  if (present(reset_avg)) then ; if (reset_avg) forces%dt_force_accum = 0.0 ; endif
  wt1 = 0.0 ; wt2 = 1.0
  if (present(dt_forcing)) then
    if ((forces%dt_force_accum > 0.0) .and. (dt_forcing > 0.0)) then
      wt1 = forces%dt_force_accum / (forces%dt_force_accum + dt_forcing)
      wt2 = 1.0 - wt1
    endif
    if (dt_forcing > 0.0) then
      forces%dt_force_accum = max(forces%dt_force_accum, 0.0) + dt_forcing
    else
      forces%dt_force_accum = 0.0 ! Reset the averaging time interval.
    endif
  else
    forces%dt_force_accum = 0.0 ! Reset the averaging time interval.
  endif

  ! applied surface pressure from atmosphere and cryosphere
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * US%Pa_to_RL2_T2*IOB%p(i-i0,j-j0)
        forces%p_surf(i,j) = MIN(forces%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * US%Pa_to_RL2_T2*IOB%p(i-i0,j-j0)
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

  ! Set the wind stresses and ustar.
  if (wt1 <= 0.0) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux=forces%taux, tauy=forces%tauy, &
                              tau_halo=1)
    if (associated(forces%ustar)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=forces%ustar)
    if (associated(forces%tau_mag)) &
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, mag_tau=forces%tau_mag)
  else
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux=forces%taux, tauy=forces%tauy, &
                              tau_halo=1)
    if (associated(forces%ustar)) then
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=ustar_tmp)
      do j=js,je ; do i=is,ie
        forces%ustar(i,j) = wt1*forces%ustar(i,j) + wt2*ustar_tmp(i,j)
      enddo ; enddo
    endif
    if (associated(forces%tau_mag)) then
      call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, mag_tau=tau_mag_tmp)
      do j=js,je ; do i=is,ie
        forces%tau_mag(i,j) = wt1*forces%tau_mag(i,j) + wt2*tau_mag_tmp(i,j)
      enddo ; enddo
    endif
  endif

  ! Find the net mass source in the input forcing without other adjustments.
  if (CS%approx_net_mass_src .and. associated(forces%net_mass_src)) then
    net_mass_src(:,:) = 0.0
    i0 = is - isc_bnd ; j0 = js - jsc_bnd
    do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
      if (associated(IOB%lprec)) &
        net_mass_src(i,j) = net_mass_src(i,j) + kg_m2_s_conversion * IOB%lprec(i-i0,j-j0)
      if (associated(IOB%fprec)) &
        net_mass_src(i,j) = net_mass_src(i,j) + kg_m2_s_conversion * IOB%fprec(i-i0,j-j0)
      if (associated(IOB%runoff)) &
        net_mass_src(i,j) = net_mass_src(i,j) + kg_m2_s_conversion * IOB%runoff(i-i0,j-j0)
      if (associated(IOB%calving)) &
        net_mass_src(i,j) = net_mass_src(i,j) + kg_m2_s_conversion * IOB%calving(i-i0,j-j0)
      if (associated(IOB%q_flux)) &
        net_mass_src(i,j) = net_mass_src(i,j) - kg_m2_s_conversion * IOB%q_flux(i-i0,j-j0)
    endif ; enddo ; enddo
    if (wt1 <= 0.0) then
      do j=js,je ; do i=is,ie
        forces%net_mass_src(i,j) = wt2*net_mass_src(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        forces%net_mass_src(i,j) = wt1*forces%net_mass_src(i,j) + wt2*net_mass_src(i,j)
      enddo ; enddo
    endif
    forces%net_mass_src_set = .true.
  else
    forces%net_mass_src_set = .false.
  endif

  ! Obtain optional ice-berg related fluxes from the IOB type:
  if (associated(IOB%area_berg)) then ; do j=js,je ; do i=is,ie
    forces%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)
  enddo ; enddo ; endif

  if (associated(IOB%mass_berg)) then ; do j=js,je ; do i=is,ie
    forces%mass_berg(i,j) = US%m_to_Z*US%kg_m3_to_R * IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)
  enddo ; enddo ; endif

  ! Obtain sea ice related dynamic fields
  if (associated(IOB%ice_rigidity)) then
    do j=js,je ; do i=is,ie
      rigidity_at_h(i,j) = US%m_to_L**3*US%Z_to_L*US%T_to_s * IOB%ice_rigidity(i-i0,j-j0) * G%mask2dT(i,j)
    enddo ; enddo
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
module procedure extract_IOB_stresses
  real, dimension(SZI_(G),SZJ_(G)) :: taux_in_A   ! Zonal wind stresses [R Z L T-2 ~> Pa] at h points
  real, dimension(SZI_(G),SZJ_(G)) :: tauy_in_A   ! Meridional wind stresses [R Z L T-2 ~> Pa] at h points
  real, dimension(SZIB_(G),SZJ_(G)) :: taux_in_C  ! Zonal wind stresses [R Z L T-2 ~> Pa] at u points
  real, dimension(SZI_(G),SZJB_(G)) :: tauy_in_C  ! Meridional wind stresses [R Z L T-2 ~> Pa] at v points
  real, dimension(SZIB_(G),SZJB_(G)) :: taux_in_B ! Zonal wind stresses [R Z L T-2 ~> Pa] at q points
  real, dimension(SZIB_(G),SZJB_(G)) :: tauy_in_B ! Meridional wind stresses [R Z L T-2 ~> Pa] at q points
  real :: gustiness     ! unresolved gustiness that contributes to ustar [R Z2 T-2 ~> Pa]
  real :: Irho0         ! Inverse of the Boussinesq mean density [R-1 ~> m3 kg-1]
  real :: taux2, tauy2  ! squared wind stresses [R2 Z2 L2 T-4 ~> Pa2]
  real :: tau_mag       ! magnitude of the wind stress [R Z2 T-2 ~> Pa]
  real :: stress_conversion ! A unit conversion factor from Pa times any stress multiplier [R Z L T-2 Pa-1 ~> 1]
  real :: Pa_to_RZ2_T2  ! The combination of unit conversion factors used for mag_tau [R Z2 T-2 Pa-1 ~> 1]
  logical :: do_ustar, do_gustless, do_tau_mag, do_gustless_tau_mag
  integer :: wind_stagger  ! AGRID, BGRID_NE, or CGRID_NE (integers from MOM_domains)
  integer :: i, j, is, ie, js, je, ish, ieh, jsh, jeh, Isqh, Ieqh, Jsqh, Jeqh, i0, j0, halo
  halo = 0 ; if (present(tau_halo)) halo = tau_halo
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  ish  = G%isc-halo  ; ieh   = G%iec+halo  ; jsh  = G%jsc-halo  ; jeh  = G%jec+halo
  Isqh = G%IscB-halo ; Ieqh  = G%IecB+halo ; Jsqh = G%JscB-halo ; Jeqh = G%JecB+halo
  i0 = is - index_bounds(1) ; j0 = js - index_bounds(3)

  IRho0 = 1.0 / CS%Rho0
  stress_conversion = US%Pa_to_RLZ_T2 * CS%wind_stress_multiplier

  do_ustar = present(ustar) ; do_gustless = present(gustless_ustar)
  do_tau_mag = present(mag_tau) ; do_gustless_tau_mag = present(gustless_mag_tau)

  wind_stagger = CS%wind_stagger
  if ((IOB%wind_stagger == AGRID) .or. (IOB%wind_stagger == BGRID_NE) .or. &
      (IOB%wind_stagger == CGRID_NE)) wind_stagger = IOB%wind_stagger

  if (associated(IOB%u_flux).neqv.associated(IOB%v_flux)) call MOM_error(FATAL,"extract_IOB_stresses: "//&
            "associated(IOB%u_flux) /= associated(IOB%v_flux !!!")
  if (present(taux).neqv.present(tauy)) call MOM_error(FATAL,"extract_IOB_stresses: "//&
            "present(taux) /= present(tauy) !!!")

  ! Set surface momentum stress related fields as a function of staggering.
  if (present(taux) .or. present(tauy) .or. &
      ((do_ustar .or. do_tau_mag .or. do_gustless .or. do_gustless_tau_mag) &
       .and. .not.associated(IOB%stress_mag)) ) then

    if (wind_stagger == BGRID_NE) then
      taux_in_B(:,:) = 0.0 ; tauy_in_B(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do J=js,je ; do I=is,ie
          taux_in_B(I,J) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_B(I,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (G%symmetric) call fill_symmetric_edges(taux_in_B, tauy_in_B, G%Domain, stagger=BGRID_NE)
      call pass_vector(taux_in_B, tauy_in_B, G%Domain, stagger=BGRID_NE, halo=max(1,halo))

      if (present(taux).and.present(tauy)) then
        do j=jsh,jeh ; do I=Isqh,Ieqh
          taux(I,j) = 0.0
          if ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0.0) &
            taux(I,j) = (G%mask2dBu(I,J)*taux_in_B(I,J) + G%mask2dBu(I,J-1)*taux_in_B(I,J-1)) / &
                        (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
        enddo ; enddo
        do J=Jsqh,Jeqh ; do i=ish,ieh
          tauy(i,J) = 0.0
          if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0.0) &
            tauy(i,J) = (G%mask2dBu(I,J)*tauy_in_B(I,J) + G%mask2dBu(I-1,J)*tauy_in_B(I-1,J)) / &
                        (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
        enddo ; enddo
      endif
    elseif (wind_stagger == AGRID) then
      taux_in_A(:,:) = 0.0 ; tauy_in_A(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do j=js,je ; do i=is,ie
          taux_in_A(i,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_A(i,j) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (halo == 0) then
        call pass_vector(taux_in_A, tauy_in_A, G%Domain, To_All+Omit_Corners, stagger=AGRID, halo=1)
      else
        call pass_vector(taux_in_A, tauy_in_A, G%Domain, stagger=AGRID, halo=max(1,halo))
      endif

      if (present(taux)) then ; do j=jsh,jeh ; do I=Isqh,Ieqh
        taux(I,j) = 0.0
        if ((G%mask2dT(i,j) + G%mask2dT(i+1,j)) > 0.0) &
          taux(I,j) = (G%mask2dT(i,j)*taux_in_A(i,j) + G%mask2dT(i+1,j)*taux_in_A(i+1,j)) / &
                      (G%mask2dT(i,j) + G%mask2dT(i+1,j))
      enddo ; enddo ; endif

      if (present(tauy)) then ; do J=Jsqh,Jeqh ; do i=ish,ieh
        tauy(i,J) = 0.0
        if ((G%mask2dT(i,j) + G%mask2dT(i,j+1)) > 0.0) &
          tauy(i,J) = (G%mask2dT(i,j)*tauy_in_A(i,j) + G%mask2dT(i,J+1)*tauy_in_A(i,j+1)) / &
                      (G%mask2dT(i,j) + G%mask2dT(i,j+1))
      enddo ; enddo ; endif

    else ! C-grid wind stresses.
      taux_in_C(:,:) = 0.0 ; tauy_in_C(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do j=js,je ; do i=is,ie
          taux_in_C(I,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_C(i,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (G%symmetric) call fill_symmetric_edges(taux_in_C, tauy_in_C, G%Domain)
      call pass_vector(taux_in_C, tauy_in_C, G%Domain, halo=max(1,halo))

      if (present(taux).and.present(tauy)) then
        do j=jsh,jeh ; do I=Isqh,Ieqh
          taux(I,j) = G%mask2dCu(I,j)*taux_in_C(I,j)
        enddo ; enddo
        do J=Jsqh,Jeqh ; do i=ish,ieh
          tauy(i,J) = G%mask2dCv(i,J)*tauy_in_C(i,J)
        enddo ; enddo
      endif
    endif   ! endif for extracting wind stress fields with various staggerings
  endif

  if (do_ustar .or. do_tau_mag .or. do_gustless .or. do_gustless_tau_mag) then
    ! Set surface friction velocity directly or as a function of staggering.
    ! ustar is required for the bulk mixed layer formulation and other turbulent mixing
    ! parametizations. The background gustiness (for example with a relatively small value
    ! of 0.02 Pa) is intended to give reasonable behavior in regions of very weak winds.
    if (associated(IOB%stress_mag)) then
      Pa_to_RZ2_T2 = US%Pa_to_RLZ_T2 * US%L_to_Z

      if (do_ustar .or. do_tau_mag) then ; do j=js,je ; do i=is,ie
        gustiness = CS%gust_const
        if (CS%read_gust_2d) then
          if ((wind_stagger == CGRID_NE) .or. &
              ((wind_stagger == AGRID) .and. (G%mask2dT(i,j) > 0.0)) .or. &
              ((wind_stagger == BGRID_NE) .and. &
               (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
                (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0.0)) ) &
            gustiness = CS%gust(i,j)
        endif
        if (do_tau_mag) &
          mag_tau(i,j) = gustiness + Pa_to_RZ2_T2*IOB%stress_mag(i-i0,j-j0)
        if (do_gustless_tau_mag) &
          gustless_mag_tau(i,j) = Pa_to_RZ2_T2*IOB%stress_mag(i-i0,j-j0)
        if (do_ustar) &
          ustar(i,j) = sqrt(gustiness*IRho0 + IRho0*Pa_to_RZ2_T2*IOB%stress_mag(i-i0,j-j0))
      enddo ; enddo ; endif
      if (CS%answer_date < 20190101) then
        if (do_gustless) then ; do j=js,je ; do i=is,ie
          gustless_ustar(i,j) = sqrt(Pa_to_RZ2_T2*IOB%stress_mag(i-i0,j-j0) / CS%Rho0)
        enddo ; enddo ; endif
      else
        if (do_gustless) then ; do j=js,je ; do i=is,ie
          gustless_ustar(i,j) = sqrt(IRho0 * Pa_to_RZ2_T2*IOB%stress_mag(i-i0,j-j0))
        enddo ; enddo ; endif
      endif
    elseif (wind_stagger == BGRID_NE) then
      do j=js,je ; do i=is,ie
        tau_mag = 0.0 ; gustiness = CS%gust_const
        if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
             (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0.0) then
          tau_mag = US%L_to_Z * sqrt(((G%mask2dBu(I,J)*((taux_in_B(I,J)**2) + (tauy_in_B(I,J)**2)) + &
              G%mask2dBu(I-1,J-1)*((taux_in_B(I-1,J-1)**2) + (tauy_in_B(I-1,J-1)**2))) + &
             (G%mask2dBu(I,J-1)*((taux_in_B(I,J-1)**2) + (tauy_in_B(I,J-1)**2)) + &
              G%mask2dBu(I-1,J)*((taux_in_B(I-1,J)**2) + (tauy_in_B(I-1,J)**2))) ) / &
            ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
          if (CS%read_gust_2d) gustiness = CS%gust(i,j)
        endif
        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (do_tau_mag) mag_tau(i,j) = gustiness + tau_mag
        if (do_gustless_tau_mag) gustless_mag_tau(i,j) = tau_mag
        if (CS%answer_date < 20190101) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    elseif (wind_stagger == AGRID) then
      do j=js,je ; do i=is,ie
        tau_mag = G%mask2dT(i,j) * US%L_to_Z * sqrt((taux_in_A(i,j)**2) + (tauy_in_A(i,j)**2))
        gustiness = CS%gust_const
        if (CS%read_gust_2d .and. (G%mask2dT(i,j) > 0.0)) gustiness = CS%gust(i,j)
        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (do_tau_mag) mag_tau(i,j) = gustiness + tau_mag
        if (do_gustless_tau_mag) gustless_mag_tau(i,j) = tau_mag
        if (CS%answer_date < 20190101) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    else  ! C-grid wind stresses.
      do j=js,je ; do i=is,ie
        taux2 = 0.0 ; tauy2 = 0.0
        if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0.0) &
          taux2 = (G%mask2dCu(I-1,j)*(taux_in_C(I-1,j)**2) + G%mask2dCu(I,j)*(taux_in_C(I,j)**2)) / &
                  (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
        if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0.0) &
          tauy2 = (G%mask2dCv(i,J-1)*(tauy_in_C(i,J-1)**2) + G%mask2dCv(i,J)*(tauy_in_C(i,J)**2)) / &
                  (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
        tau_mag = US%L_to_Z * sqrt(taux2 + tauy2)

        gustiness = CS%gust_const
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)

        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (do_tau_mag) mag_tau(i,j) = gustiness + tau_mag
        if (do_gustless_tau_mag) gustless_mag_tau(i,j) = tau_mag
        if (CS%answer_date < 20190101) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    endif ! endif for wind friction velocity fields
  endif

end procedure extract_IOB_stresses
module procedure apply_flux_adjustments
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: temp_at_h ! Various fluxes at h points
  integer :: isc, iec, jsc, jec, i, j
  logical :: overrode_h
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call data_override(G%Domain, 'hflx_adj', temp_at_h, Time, override=overrode_h, &
                     scale=US%W_m2_to_QRZ_T)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%heat_added(i,j) = fluxes%heat_added(i,j) + temp_at_h(i,j) * G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%heat_added, G%Domain)

  call data_override(G%Domain, 'sflx_adj', temp_at_h, Time, override=overrode_h, &
                     scale=US%kg_m2s_to_RZ_T)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%salt_flux_added(i,j) = fluxes%salt_flux_added(i,j) + temp_at_h(i,j) * G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%salt_flux_added, G%Domain)

  call data_override(G%Domain, 'prcme_adj', temp_at_h, Time, override=overrode_h, &
                     scale=US%kg_m2s_to_RZ_T)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%vprec(i,j) = fluxes%vprec(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%vprec, G%Domain)
end procedure apply_flux_adjustments
module procedure apply_force_adjustments
  real, dimension(SZI_(G),SZJ_(G)) :: tempx_at_h ! Delta to zonal wind stress at h points [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: tempy_at_h ! Delta to meridional wind stress at h points [R Z L T-2 ~> Pa]
  integer :: isc, iec, jsc, jec, i, j
  real :: dLonDx, dLonDy ! The change in longitude across the cell in the x- and y-directions [degrees_E]
  real :: rDlon ! The magnitude of the change in longitude [degrees_E] and then its inverse [degrees_E-1]
  real :: cosA, sinA  ! The cosine and sine of the angle between the grid and true north [nondim]
  real :: zonal_tau, merid_tau ! True zonal and meridional wind stresses [R Z L T-2 ~> Pa]
  logical :: overrode_x, overrode_y
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  tempx_at_h(:,:) = 0.0 ; tempy_at_h(:,:) = 0.0
  ! Either reads data or leaves contents unchanged
  overrode_x = .false. ; overrode_y = .false.
  call data_override(G%Domain, 'taux_adj', tempx_at_h(isc:iec,jsc:jec), Time, &
                     override=overrode_x, scale=US%Pa_to_RLZ_T2)
  call data_override(G%Domain, 'tauy_adj', tempy_at_h(isc:iec,jsc:jec), Time, &
                     override=overrode_y, scale=US%Pa_to_RLZ_T2)

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
      zonal_tau = tempx_at_h(i,j)
      merid_tau = tempy_at_h(i,j)
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
  real :: utide             ! The RMS tidal velocity [Z T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    utide_2d                ! A 2d array of RMS tidal velocities [Z T-1 ~> m s-1].
  real :: Flux_const_dflt   ! A default piston velocity for restoring surface properties [m day-1]
  logical :: Boussinesq       ! If true, this run is fully Boussinesq
  logical :: semi_Boussinesq  ! If true, this run is partially non-Boussinesq
  real :: rho_TKE_tidal     ! The constant bottom density used to translate tidal amplitudes into
  logical :: new_sim              ! False if this simulation was started from a restart file
  logical :: iceberg_flux_diags   ! If true, diagnostics of fluxes from icebergs are available.
  logical :: fix_ustar_gustless_bug  ! If false, include a bug using an older run-time parameter.
  logical :: test_value  ! This is used to determine whether a logical parameter is being set explicitly.
  logical :: explicit_bug, explicit_fix ! These indicate which parameters are set explicitly.
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  type(time_type)    :: Time_frc
  type(directories)  :: dirs      ! A structure containing relevant directory paths and input filenames.
  character(len=200) :: TideAmp_file, gust_file, salt_file, temp_file ! Input file names.
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_surface_forcing"  ! This module's name.
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
  call log_version(param_file, mdl, version, "", log_to_all=.true., debugging=.true.)

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "BOUSSINESQ", Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "SEMI_BOUSSINESQ", semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=.true.)
  CS%nonBous = .not.(Boussinesq .or. semi_Boussinesq)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R) ! (, do_not_log=CS%nonBous)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the "//&
                 "atmosphere and floating sea-ice or ice shelves. This is "//&
                 "needed because the FMS coupling structure does not "//&
                 "limit the water that can be frozen out of the ocean and "//&
                 "the ice-ocean heat fluxes are treated explicitly.  No "//&
                 "limit is applied if a negative value is used.", &
                 units="Pa", default=-1.0, scale=US%Pa_to_RL2_T2)
  call get_param(param_file, mdl, "RESTORE_SALINITY", CS%restore_salt, &
                 "If true, the coupled driver will add a globally-balanced "//&
                 "fresh-water flux that drives sea-surface salinity "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RESTORE_TEMPERATURE", CS%restore_temp, &
                 "If true, the coupled driver will add a  "//&
                 "heat flux that drives sea-surface temperature "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_TO_ZERO", &
                 CS%adjust_net_srestore_to_zero, &
                 "If true, adjusts the salinity restoring seen to zero "//&
                 "whether restoring is via a salt flux or virtual precip.",&
                 default=CS%restore_salt)
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
  call get_param(param_file, mdl, "APPROX_NET_MASS_SRC", CS%approx_net_mass_src, &
                 "If true, use the net mass sources from the ice-ocean "//&
                 "boundary type without any further adjustments to drive "//&
                 "the ocean dynamics.  The actual net mass source may differ "//&
                 "due to internal corrections.", default=.false.)

  if (present(wind_stagger)) then
    if     (wind_stagger == AGRID)    then ; stagger = 'AGRID'
    elseif (wind_stagger == BGRID_NE) then ; stagger = 'BGRID_NE'
    elseif (wind_stagger == CGRID_NE) then ; stagger = 'CGRID_NE'
    else ; stagger = 'UNKNOWN' ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                      trim(stagger)// "is invalid.") ; endif
    call log_param(param_file, mdl, "WIND_STAGGER", stagger, &
                   "The staggering of the input wind stress field "//&
                   "from the coupler that is actually used.")
    CS%wind_stagger = wind_stagger
  else
    call get_param(param_file, mdl, "WIND_STAGGER", stagger, &
                   "A case-insensitive character string to indicate the "//&
                   "staggering of the input wind stress field.  Valid "//&
                   "values are 'A', 'B', or 'C'.", default="C")
    if     (uppercase(stagger(1:1)) == 'A') then ; CS%wind_stagger = AGRID
    elseif (uppercase(stagger(1:1)) == 'B') then ; CS%wind_stagger = BGRID_NE
    elseif (uppercase(stagger(1:1)) == 'C') then ; CS%wind_stagger = CGRID_NE
    else ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                          trim(stagger)//" is invalid.") ; endif
  endif

  call get_param(param_file, mdl, "WIND_STRESS_MULTIPLIER", CS%wind_stress_multiplier, &
                 "A factor multiplying the wind-stress given to the ocean by the "//&
                 "coupler. This is used for testing and should be =1.0 for any "//&
                 "production runs.", units="nondim", default=1.0)

  if (CS%restore_salt) then
    call get_param(param_file, mdl, "FLUXCONST", Flux_const_dflt, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", default=0.0)
    call get_param(param_file, mdl, "FLUXCONST_SALT", CS%Flux_const_salt, &
                 "The constant that relates the restoring surface salt fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", default=Flux_const_dflt, scale=US%m_to_Z*US%T_to_s)
    ! Finish converting CS%Flux_const_salt from m day-1 to [Z T-1 ~> m s-1].  Ideally this would be
    ! included in the scale factors above, but doing so would change answers because a/b /= a*(1/b).
    CS%Flux_const_salt = CS%Flux_const_salt / 86400.0
    call get_param(param_file, mdl, "SALT_RESTORE_FILE", CS%salt_restore_file, &
                 "A file in which to find the surface salinity to use for restoring.", &
                 default="salt_restore.nc")
    call get_param(param_file, mdl, "SALT_RESTORE_VARIABLE", CS%salt_restore_var_name, &
                 "The name of the surface salinity variable to read from "//&
                 "SALT_RESTORE_FILE for restoring salinity.", &
                 default="salt")
    call get_param(param_file, mdl, "SALT_RESTORE_PRACTICAL_SALINITY", CS%salt_restore_is_practical, &
                 "Specifies if the restoring surface salinity variable is practical salinity.  If this "//&
                 "flag is set to false it is assumed that the salinity is absolute salinity.", default=.false.)
    call get_param(param_file, mdl, "SRESTORE_AS_SFLUX", CS%salt_restore_as_sflux, &
                 "If true, the restoring of salinity is applied as a salt "//&
                 "flux instead of as a freshwater flux.", default=.false.)
    call get_param(param_file, mdl, "MAX_DELTA_SRESTORE", CS%max_delta_srestore, &
                 "The maximum salinity difference used in restoring terms.", &
                 units="PSU or g kg-1", default=999.0, scale=US%ppt_to_S)
    call get_param(param_file, mdl, "MASK_SRESTORE_UNDER_ICE", CS%mask_srestore_under_ice, &
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

  if (CS%restore_temp) then
    call get_param(param_file, mdl, "FLUXCONST", Flux_const_dflt, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", default=0.0)
    call get_param(param_file, mdl, "FLUXCONST_TEMP", CS%Flux_const_temp, &
                 "The constant that relates the restoring surface temperature fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", default=Flux_const_dflt, scale=US%m_to_Z*US%T_to_s)
    ! Finish converting CS%Flux_const_temp from [m day-1] to [Z T-1 ~> m s-1].  Ideally this would be
    ! included in the scale factors above, but doing so would change answers because a/b /= a*(1/b).
    CS%Flux_const_temp = CS%Flux_const_temp / 86400.0
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

    call get_param(param_file, mdl, "SPEAR_ECDA_SST_RESTORE_TFREEZE", CS%trestore_SPEAR_ECDA, &
                 "If true, modify SST restoring field using SSS state. This only modifies the "//&
                 "restoring data that is within 0.0001degC of -1.8degC.", default=.false.)
  else
    CS%trestore_SPEAR_ECDA = .false. ! Needed to toggle logging of SPEAR_DTFREEZE_DS
  endif
  call get_param(param_file, mdl, "SPEAR_DTFREEZE_DS", CS%SPEAR_dTf_dS, &
                 "The derivative of the freezing temperature with salinity.", &
                 units="degC ppt-1", default=-0.054, scale=US%degC_to_C*US%S_to_ppt, &
                 do_not_log=.not.CS%trestore_SPEAR_ECDA)
  call get_param(param_file, mdl, "RESTORE_FLUX_RHO", CS%rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat "//&
                 "fluxes with RESTORE_SALINITY or RESTORE_TEMPERATURE.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=.not.(CS%restore_temp.or.CS%restore_salt))

  ! Optionally read tidal amplitude from input file [Z T-1 ~> m s-1] on model grid.
  ! Otherwise use default tidal amplitude for bottom frictionally-generated
  ! dissipation. Default cd_tides is chosen to yield approx 1 TWatt of
  ! work done against tides globally using OSU tidal amplitude.
  ! Note that the slightly unusual length scaling is deliberate, because the tidal
  ! amplitudes are used to set the friction velocity.
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
  call get_param(param_file, mdl, "TKE_TIDAL_RHO", rho_TKE_tidal, &
                 "The constant bottom density used to translate tidal amplitudes into the tidal "//&
                 "bottom TKE input used with INT_TIDE_DISSIPATION.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R*US%Z_to_L**2, &
                 do_not_log=.not.(CS%read_TIDEAMP.or.(CS%utide>0.0)))

  call safe_alloc_ptr(CS%BBL_tidal_dis,isd,ied,jsd,jed)
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed)

  if (CS%read_TIDEAMP) then
    TideAmp_file = trim(CS%inputdir) // trim(TideAmp_file)
    ! NOTE: There are certain cases where FMS is unable to read this file, so
    ! we use read_netCDF_data in place of MOM_read_data.
    utide_2d(:,:) = 0.0
    call read_netCDF_data(TideAmp_file, 'tideamp', utide_2d, G%Domain, &
        rescale=US%m_to_Z*US%T_to_s)
    do j=jsd,jed ; do i=isd,ied
      utide = utide_2d(i,j)
      CS%BBL_tidal_dis(i,j) = G%mask2dT(i,j)*rho_TKE_tidal*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  else
    do j=jsd,jed ; do i=isd,ied
      utide = CS%utide
      CS%BBL_tidal_dis(i,j) = rho_TKE_tidal*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  endif

  call time_interp_external_init()

  ! Optionally read a x-y gustiness field in place of a global constant.
  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2*US%L_to_Z)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.", fail_if_missing=.true.)

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed)
    gust_file = trim(CS%inputdir) // trim(gust_file)
    ! NOTE: There are certain cases where FMS is unable to read this file, so
    ! we use read_netCDF_data in place of MOM_read_data.
    call read_netCDF_data(gust_file, 'gustiness', CS%gust, G%Domain, &
                          rescale=US%Pa_to_RLZ_T2*US%L_to_Z) ! units in file should be [Pa]
  endif
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "SURFACE_FORCING_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the gustiness "//&
                 "calculations.  Values below 20190101 recover the answers from the end "//&
                 "of 2018, while higher values use a simpler expression to calculate gustiness.", &
                 default=default_answer_date)

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
                 units="m s-2", default=9.80, scale=US%Z_to_m*US%m_s_to_L_T**2)
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
  call get_param(param_file, mdl, "ALLOW_CARBON_FLUX_EXCHANGE", CS%allow_carbon_flux_exchange, &
                 "If true, makes available fluxes and diagnostics of carbon in runoff "//&
                 "within MOM6.", default=.false.)
  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles, &
                                   use_berg_fluxes=iceberg_flux_diags, &
                                   use_carbon_runoff=CS%allow_carbon_flux_exchange)

  call get_param(param_file, mdl, "ALLOW_FLUX_ADJUSTMENTS", CS%allow_flux_adjustments, &
                 "If true, allows flux adjustments to specified via the "//&
                 "data_table using the component name 'OCN'.", default=.false.)

  call get_param(param_file, mdl, "CHECK_NO_LAND_FLUXES", CS%check_no_land_fluxes, &
                 "If true, checks that values from IOB fluxes are zero "//&
                 "above land points (i.e. G%mask2dT = 0).", default=.false., &
                 debuggingParam=.true.)

  call data_override_init(G%Domain)

  if (CS%restore_salt) then
    salt_file = trim(CS%inputdir) // trim(CS%salt_restore_file)
    CS%srestore_handle = init_external_field(salt_file, CS%salt_restore_var_name, MOM_domain=G%Domain)
    call safe_alloc_ptr(CS%srestore_mask,isd,ied,jsd,jed) ; CS%srestore_mask(:,:) = 1.0
    if (CS%mask_srestore) then ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'salt_restore_mask.nc'
      call MOM_read_data(flnam,'mask', CS%srestore_mask, G%domain, timelevel=1)
    endif
  endif

  if (CS%restore_temp) then
    temp_file = trim(CS%inputdir) // trim(CS%temp_restore_file)
    CS%trestore_handle = init_external_field(temp_file, CS%temp_restore_var_name, MOM_domain=G%Domain)
    call safe_alloc_ptr(CS%trestore_mask,isd,ied,jsd,jed) ; CS%trestore_mask(:,:) = 1.0
    if (CS%mask_trestore) then  ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'temp_restore_mask.nc'
      call MOM_read_data(flnam, 'mask', CS%trestore_mask, G%domain, timelevel=1)
    endif
  endif

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
!#CTRL#  call register_ctrl_forcing_restarts(G, param_file, CS%ctrl_forcing_CSp, &
!#CTRL#                                      CS%restart_CSp)
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

!#CTRL#  call controlled_forcing_init(Time, G, US, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end procedure surface_forcing_init
module procedure surface_forcing_end
  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!#CTRL#  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end procedure surface_forcing_end
module procedure ice_ocn_bnd_type_chksum
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  integer :: outunit ! The output unit to write to
  root = is_root_pe()
  outunit = stdout_if_root()

  if (root) write(outunit,*) "BEGIN CHECKSUM(ice_ocean_boundary_type):: ", id, timestep
  chks = field_chksum( iobt%u_flux         ) ; if (root) write(outunit,100) 'iobt%u_flux         ', chks
  chks = field_chksum( iobt%v_flux         ) ; if (root) write(outunit,100) 'iobt%v_flux         ', chks
  chks = field_chksum( iobt%t_flux         ) ; if (root) write(outunit,100) 'iobt%t_flux         ', chks
  chks = field_chksum( iobt%q_flux         ) ; if (root) write(outunit,100) 'iobt%q_flux         ', chks
  chks = field_chksum( iobt%salt_flux      ) ; if (root) write(outunit,100) 'iobt%salt_flux      ', chks
  chks = field_chksum( iobt%lw_flux        ) ; if (root) write(outunit,100) 'iobt%lw_flux        ', chks
  chks = field_chksum( iobt%sw_flux_vis_dir) ; if (root) write(outunit,100) 'iobt%sw_flux_vis_dir', chks
  chks = field_chksum( iobt%sw_flux_vis_dif) ; if (root) write(outunit,100) 'iobt%sw_flux_vis_dif', chks
  chks = field_chksum( iobt%sw_flux_nir_dir) ; if (root) write(outunit,100) 'iobt%sw_flux_nir_dir', chks
  chks = field_chksum( iobt%sw_flux_nir_dif) ; if (root) write(outunit,100) 'iobt%sw_flux_nir_dif', chks
  chks = field_chksum( iobt%lprec          ) ; if (root) write(outunit,100) 'iobt%lprec          ', chks
  chks = field_chksum( iobt%fprec          ) ; if (root) write(outunit,100) 'iobt%fprec          ', chks
  chks = field_chksum( iobt%runoff         ) ; if (root) write(outunit,100) 'iobt%runoff         ', chks
  chks = field_chksum( iobt%calving        ) ; if (root) write(outunit,100) 'iobt%calving        ', chks
  chks = field_chksum( iobt%p              ) ; if (root) write(outunit,100) 'iobt%p              ', chks
  if (associated(iobt%shelf_sfc_mass_flux)) then
     chks = field_chksum( iobt%shelf_sfc_mass_flux ) ; if (root) write(outunit,100) 'iobt%shelf_sfc_mass_flux     ',&
        chks
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
  if (associated(iobt%excess_salt)) then
    chks = field_chksum( iobt%excess_salt    ) ; if (root) write(outunit,100) 'iobt%excess_salt    ', chks
  endif
100 FORMAT("   CHECKSUM::",A20," = ",Z20)

  call coupler_type_write_chksums(iobt%fluxes, outunit, 'iobt%')

end procedure ice_ocn_bnd_type_chksum
module procedure check_mask_val_consistency
  character(len=48) :: ci, cj !< model local grid cell indices as strings
  character(len=48) :: ciglo, cjglo !< model global grid cell indices as strings
  character(len=48) :: cval !< value to be displayed
  character(len=256) :: error_message !< error message to be displayed
  if ((mask == 0.) .and. (val /= 0.)) then
    write(ci, '(I8)') i
    write(cj, '(I8)') j
    write(ciglo, '(I8)') i + G%HI%idg_offset
    write(cjglo, '(I8)') j + G%HI%jdg_offset
    write(cval, '(E22.16)') val
    error_message = "MOM_surface_forcing: found non-zero value (="//trim(cval)//") over land "//&
                    "for variable "//trim(varname)//" at local point (i, j) = ("//trim(ci)//", "//trim(cj)//&
                    ", global point (iglo, jglo) = ("//trim(ciglo)//", "//trim(cjglo)//")"
    call MOM_error(WARNING, error_message)
  endif

end procedure check_mask_val_consistency
end submodule MOM_surface_forcing_gfdl_s
