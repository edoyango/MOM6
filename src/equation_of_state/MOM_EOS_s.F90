submodule (MOM_EOS) MOM_EOS_s
  implicit none
contains
module procedure calculate_density_elem
  real :: Ta      ! An array of temperatures [degC]
  real :: Sa      ! An array of salinities [ppt]
  real :: pres    ! An mks version of the pressure to use [Pa]
  real :: rho_mks ! An mks version of the density to be returned [kg m-3]
  real :: rho_scale  ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  pres = EOS%RL2_T2_to_Pa * pressure
  Ta = EOS%C_to_degC * T
  Sa = EOS%S_to_ppt * S

  if (present(rho_ref)) then
    rho_mks = EOS%type%density_anomaly_elem(Ta, Sa, pres, EOS%R_to_kg_m3*rho_ref)
  else
    rho_mks = EOS%type%density_elem(Ta, Sa, pres)
  endif

  ! Rescale the output density to the desired units.
  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  calculate_density_elem = rho_scale * rho_mks

end procedure calculate_density_elem
module procedure calculate_density_scalar
  real :: Ta      ! An array of temperatures [degC]
  real :: Sa      ! An array of salinities [ppt]
  real :: pres    ! An mks version of the pressure to use [Pa]
  real :: rho_mks ! An mks version of the density to be returned [kg m-3]
  real :: rho_scale  ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  pres = EOS%RL2_T2_to_Pa * pressure
  Ta = EOS%C_to_degC * T
  Sa = EOS%S_to_ppt * S

  if (present(rho_ref)) then
    rho_mks = EOS%type%density_anomaly_elem(Ta, Sa, pres, EOS%R_to_kg_m3*rho_ref)
  else
    rho_mks = EOS%type%density_elem(Ta, Sa, pres)
  endif

  ! Rescale the output density to the desired units.
  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  rho = rho_scale * rho_mks

end procedure calculate_density_scalar
module procedure calculate_stanley_density_scalar
  real :: d2RdTT   ! Second derivative of density with temperature [R C-2 ~> kg m-3 degC-2]
  real :: d2RdST   ! Second derivative of density with temperature and salinity [R S-1 C-1 ~> kg m-3 degC-1 ppt-1]
  real :: d2RdSS   ! Second derivative of density with salinity [R S-2 ~> kg m-3 ppt-2]
  real :: d2RdSp   ! Second derivative of density with salinity and pressure [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
  real :: d2RdTp   ! Second derivative of density with temperature and pressure [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  call calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref)
  call calculate_density_second_derivs_scalar(T, S, pressure, d2RdSS, d2RdST, d2RdTT, d2RdSp, d2RdTP, EOS)

  ! Equation 25 of Stanley et al., 2020.
  rho = rho + ( 0.5 * d2RdTT * Tvar + ( d2RdST * TScov + 0.5 * d2RdSS * Svar ) )

  if (present(scale)) rho = rho * scale

end procedure calculate_stanley_density_scalar
module procedure calculate_density_1d
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real, dimension(size(rho)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(rho)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(rho)) :: Sa    ! Salinity converted to [ppt]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%R_to_kg_m3 == 1.0) .and. &
      (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call EOS%type%calculate_density_array(T, S, pressure, rho, is, npts, rho_ref=rho_ref)
  else ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    if (present(rho_ref)) then
      call EOS%type%calculate_density_array(Ta, Sa, pres, rho, is, npts, rho_ref=EOS%R_to_kg_m3*rho_ref)
    else
      call EOS%type%calculate_density_array(Ta, Sa, pres, rho, is, npts)
    endif
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    rho(i) = rho_scale * rho(i)
  enddo ; endif

end procedure calculate_density_1d
module procedure calculate_stanley_density_1d
  real, dimension(size(T)) :: &
    d2RdTT, &   ! Second derivative of density with temperature [R C-2 ~> kg m-3 degC-2]
    d2RdST, &   ! Second derivative of density with temperature and salinity [R S-1 C-1 ~> kg m-3 degC-1 ppt-1]
    d2RdSS, &   ! Second derivative of density with salinity [R S-2 ~> kg m-3 ppt-2]
    d2RdSp, &   ! Second derivative of density with salinity and pressure [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
    d2RdTp      ! Second derivative of density with temperature and pressure [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  call calculate_density_1d(T, S, pressure, rho, EOS, dom, rho_ref)
  call calculate_density_second_derivs_1d(T, S, pressure, d2RdSS, d2RdST, d2RdTT, d2RdSp, d2RdTP, EOS, dom)

  ! Equation 25 of Stanley et al., 2020.
  do i=is,ie
    rho(i) = rho(i) + ( 0.5 * d2RdTT(i) * Tvar(i) + ( d2RdST(i) * TScov(i) + 0.5 * d2RdSS(i) * Svar(i) ) )
  enddo

  if (present(scale)) then ; if (scale /= 1.0) then ; do i=is,ie
    rho(i) = scale * rho(i)
  enddo ; endif ; endif

end procedure calculate_stanley_density_1d
module procedure calculate_spec_vol_array
  integer :: j
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
              "calculate_spec_vol_array: EOS%form_of_EOS is not valid.")

  call EOS%type%calculate_spec_vol_array(T, S, pressure, specvol, start, npts, spv_ref)

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    specvol(j) = scale * specvol(j)
  enddo ; endif ; endif

end procedure calculate_spec_vol_array
module procedure calc_spec_vol_scalar
  real, dimension(1) :: Ta   ! Rescaled single element array version of temperature [degC]
  real, dimension(1) :: Sa   ! Rescaled single element array version of salinity [ppt]
  real, dimension(1) :: pres ! Rescaled single element array version of pressure [Pa]
  real, dimension(1) :: spv  ! Rescaled single element array version of specific volume [m3 kg-1]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]
  pres(1) = EOS%RL2_T2_to_Pa * pressure
  Ta(1) = EOS%C_to_degC * T ; Sa(1) = EOS%S_to_ppt * S

  if (present(spv_ref)) then
    call calculate_spec_vol_array(Ta, Sa, pres, spv, 1, 1, EOS, EOS%kg_m3_to_R*spv_ref)
  else
    call calculate_spec_vol_array(Ta, Sa, pres, spv, 1, 1, EOS)
  endif
  specvol = spv(1)

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then
    specvol = spv_scale * specvol
  endif

end procedure calc_spec_vol_scalar
module procedure calc_spec_vol_1d
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(specvol) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%kg_m3_to_R == 1.0) .and. &
      (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_spec_vol_array(T, S, pressure, specvol, is, npts, EOS, spv_ref)
  else ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    if (present(spv_ref)) then
      call calculate_spec_vol_array(Ta, Sa, pres, specvol, is, npts, EOS, EOS%kg_m3_to_R*spv_ref)
    else
      ! There is rescaling of variables, but spv_ref is not present. Passing a 0 value of spv_ref
      ! changes answers at roundoff for some equations of state, like Wright and UNESCO.
      call calculate_spec_vol_array(Ta, Sa, pres, specvol, is, npts, EOS)
    endif
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then ; do i=is,ie
    specvol(i) = spv_scale * specvol(i)
  enddo ; endif

end procedure calc_spec_vol_1d
module procedure calculate_TFreeze_scalar
  real :: p_scale ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: S_scale ! A factor to convert salinity to units of ppt [ppt S-1 ~> 1]
  real :: iS_scale! A factor to convert salinity to units of S [S ppt-1 ~> 1]
  real :: absS    ! A salinity converted to absolute salinity, only used in specific scenarios [ppt]
  real :: TFreeze_S ! The salinity for the freezing equation in model units [S ~> PSU or ppt]
  p_scale = 1.0 ; S_scale = 1.0 ; iS_scale = 1.0
  if (present(pres_scale)) p_scale = pres_scale
  if (present(scale_from_EOS)) then ; if (scale_from_EOS) then
    p_scale = EOS%RL2_T2_to_Pa
    S_scale = EOS%S_to_ppt
    iS_scale = EOS%ppt_to_S
  endif ; endif

  if (EOS%use_conT_absS) then
    ! Otherwise absS is unneeded and therefore unset
    absS = S*S_scale
    if (EOS%TFreeze_S_is_pracS) then
      TFreeze_S = gsw_sp_from_sr(absS)*iS_scale
    else
      TFreeze_S = S
    endif
  else
    TFreeze_S = S
  endif

  select case (EOS%form_of_TFreeze)
    case (TFREEZE_LINEAR)
      call calculate_TFreeze_linear(S_scale*TFreeze_S, p_scale*pressure, T_fr, EOS%TFr_S0_P0, &
                                    EOS%dTFr_dS, EOS%dTFr_dp)
    case (TFREEZE_MILLERO)
      call calculate_TFreeze_Millero(S_scale*TFreeze_S, p_scale*pressure, T_fr)
    case (TFREEZE_TEOSPOLY)
      call calculate_TFreeze_TEOS_poly(S_scale*TFreeze_S, p_scale*pressure, T_fr)
    case (TFREEZE_TEOS10)
      call calculate_TFreeze_teos10(S_scale*TFreeze_S, p_scale*pressure, T_fr)
    case default
      call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
  end select

  if (EOS%use_conT_absS .and. EOS%TFreeze_T_is_potT) then
    ! absS is set only if EOS%use_conT_absS is True
    ! absS and T_fr have physical units here and don't need converted
    T_fr = gsw_ct_from_pt(absS,T_fr)
  endif

  if (present(scale_from_EOS)) then ; if (scale_from_EOS) then
    T_fr = EOS%degC_to_C * T_fr
  endif ; endif

end procedure calculate_TFreeze_scalar
module procedure calculate_TFreeze_array
  real, dimension(size(pressure)) :: pres  ! Pressure converted to [Pa]
  real :: p_scale  ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real, dimension(size(S)) :: absS ! A salinity converted to absolute salinity, only used in specific scenarios [ppt]
  real, dimension(size(S)) :: TFreeze_S ! The salinity for the freezing equation in model units [S ~> PSU or ppt]
  integer :: j
  p_scale = 1.0 ; if (present(pres_scale)) p_scale = pres_scale

  if (EOS%use_conT_absS) then
    ! Otherwise absS is unneeded and therefore unset
    absS(:) = S(:)
    if (EOS%TFreeze_S_is_pracS) then
      TFreeze_S(:) = gsw_sp_from_sr(absS(:))
    else
      TFreeze_S(:) = S(:)
    endif
  else
    TFreeze_S(:) = S(:)
  endif

  if (p_scale == 1.0) then
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(TFreeze_S, pressure, T_fr, start, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(TFreeze_S, pressure, T_fr, start, npts)
      case (TFREEZE_TEOSPOLY)
        call calculate_TFreeze_TEOS_poly(TFreeze_S, pressure, T_fr, start, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(TFreeze_S, pressure, T_fr, start, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  else
    do j=start,start+npts-1 ; pres(j) = p_scale * pressure(j) ; enddo
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(TFreeze_S, pres, T_fr, start, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(TFreeze_S, pres, T_fr, start, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(TFreeze_S, pres, T_fr, start, npts)
      case (TFREEZE_TEOSPOLY)
        call calculate_TFreeze_TEOS_poly(TFreeze_S, pres, T_fr, start, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  endif

  if (EOS%use_conT_absS .and. EOS%TFreeze_T_is_potT) then
    ! absS is set only if EOS%use_conT_absS is True!
    T_fr(:) = gsw_ct_from_pt(absS(:),T_fr(:))
  endif


end procedure calculate_TFreeze_array
module procedure calculate_TFreeze_1d
  real, dimension(size(T_fr)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T_fr)) :: Sa    ! Salinity converted to [ppt]
  real, dimension(size(T_fr)) :: absS  ! Salinity converted to absoluate salinity [ppt]
  real, dimension(size(T_fr)) :: TFreeze_S ! The salinity for the freezing equation in model units [S ~> PSU or ppt]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(T_Fr) ; npts = 1 + ie - is
  endif

  if (EOS%use_conT_absS) then
    ! Otherwise absS is unneeded and therefore unset
    absS(:) = S(:)*EOS%S_to_ppt
    if (EOS%TFreeze_S_is_pracS) then
      TFreeze_S(:) = gsw_sp_from_sr(absS(:))*EOS%ppt_to_S
    else
      TFreeze_S(:) = S(:)
    endif
  else
    TFreeze_S(:) = S(:)
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(TFreeze_S, pressure, T_fr, is, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(TFreeze_S, pressure, T_fr, is, npts)
      case (TFREEZE_TEOSPOLY)
        call calculate_TFreeze_TEOS_poly(TFreeze_S, pressure, T_fr, is, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(TFreeze_S, pressure, T_fr, is, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Sa(i) = EOS%S_to_ppt * TFreeze_S(i)
    enddo
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(Sa, pres, T_fr, is, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(Sa, pres, T_fr, is, npts)
      case (TFREEZE_TEOSPOLY)
        call calculate_TFreeze_TEOS_poly(Sa, pres, T_fr, is, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(Sa, pres, T_fr, is, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  endif

  if (EOS%use_conT_absS .and. EOS%TFreeze_T_is_potT) then
    ! absS is set only if EOS%use_conT_absS is True!
    ! absS is in ppt and T_fr is in degC at this point.
    T_fr(:) = gsw_ct_from_pt(absS(:),T_fr(:))
  endif


  if (EOS%degC_to_C /= 1.0) then
    do i=is,ie ; T_fr(i) = EOS%degC_to_C * T_fr(i) ; enddo
  endif

end procedure calculate_TFreeze_1d
module procedure calculate_density_derivs_array
  integer :: j
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
              "calculate_density_derivs_array: EOS%form_of_EOS is not valid.")

  call EOS%type%calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, start, npts)

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    drho_dT(j) = scale * drho_dT(j)
    drho_dS(j) = scale * drho_dS(j)
  enddo ; endif ; endif

end procedure calculate_density_derivs_array
module procedure calculate_density_derivs_1d
  real, dimension(size(drho_dT)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(drho_dT)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(drho_dT)) :: Sa    ! Salinity converted to [ppt]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(drho_dT) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, is, npts, EOS)
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    call calculate_density_derivs_array(Ta, Sa, pres, drho_dT, drho_dS, is, npts, EOS)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  dRdT_scale = rho_scale * EOS%C_to_degC
  dRdS_scale = rho_scale * EOS%S_to_ppt
  if ((dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0)) then ; do i=is,ie
    drho_dT(i) = dRdT_scale * drho_dT(i)
    drho_dS(i) = dRdS_scale * drho_dS(i)
  enddo ; endif

end procedure calculate_density_derivs_1d
module procedure calculate_density_derivs_scalar
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  real :: pres(1)  ! Pressure converted to [Pa]
  real :: Ta(1)    ! Temperature converted to [degC]
  real :: Sa(1)    ! Salinity converted to [ppt]
  pres(1) = EOS%RL2_T2_to_Pa*pressure
  Ta(1) = EOS%C_to_degC * T
  Sa(1) = EOS%S_to_ppt * S

  call EOS%type%calculate_density_derivs_scalar(Ta(1), Sa(1), pres(1), drho_dT, drho_dS)

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  dRdT_scale = rho_scale * EOS%C_to_degC
  dRdS_scale = rho_scale * EOS%S_to_ppt
  if ((dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0)) then
    drho_dT = dRdT_scale * drho_dT
    drho_dS = dRdS_scale * drho_dS
  endif

end procedure calculate_density_derivs_scalar
module procedure calculate_density_second_derivs_1d
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  integer :: i, is, ie, npts
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
              "calculate_density_second_derivs: EOS%form_of_EOS is not valid.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(T) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call EOS%type%calculate_density_second_derivs_array(T, S, pressure, &
                        drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    call EOS%type%calculate_density_second_derivs_array(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    drho_dS_dS(i) = rho_scale * drho_dS_dS(i)
    drho_dS_dT(i) = rho_scale * drho_dS_dT(i)
    drho_dT_dT(i) = rho_scale * drho_dT_dT(i)
    drho_dS_dP(i) = rho_scale * drho_dS_dP(i)
    drho_dT_dP(i) = rho_scale * drho_dT_dP(i)
  enddo ; endif

  if (EOS%RL2_T2_to_Pa /= 1.0) then ; do i=is,ie
    drho_dS_dP(i) = EOS%RL2_T2_to_Pa * drho_dS_dP(i)
    drho_dT_dP(i) = EOS%RL2_T2_to_Pa * drho_dT_dP(i)
  enddo ; endif

  if (EOS%C_to_degC /= 1.0) then ; do i=is,ie
    drho_dS_dT(i) = EOS%C_to_degC * drho_dS_dT(i)
    drho_dT_dT(i) = EOS%C_to_degC**2 * drho_dT_dT(i)
    drho_dT_dP(i) = EOS%C_to_degC * drho_dT_dP(i)
  enddo ; endif

  if (EOS%S_to_ppt /= 1.0) then ; do i=is,ie
    drho_dS_dS(i) = EOS%S_to_ppt**2 * drho_dS_dS(i)
    drho_dS_dT(i) = EOS%S_to_ppt * drho_dS_dT(i)
    drho_dS_dP(i) = EOS%S_to_ppt * drho_dS_dP(i)
  enddo ; endif

end procedure calculate_density_second_derivs_1d
module procedure calculate_density_second_derivs_scalar
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: pres  ! Pressure converted to [Pa]
  real :: Ta    ! Temperature converted to [degC]
  real :: Sa    ! Salinity converted to [ppt]
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
             "calculate_density_second_derivs: EOS%form_of_EOS is not valid.")

  pres = EOS%RL2_T2_to_Pa*pressure
  Ta = EOS%C_to_degC * T
  Sa = EOS%S_to_ppt * S

  call EOS%type%calculate_density_second_derivs_scalar(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                              drho_dT_dT, drho_dS_dP, drho_dT_dP)

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then
    drho_dS_dS = rho_scale * drho_dS_dS
    drho_dS_dT = rho_scale * drho_dS_dT
    drho_dT_dT = rho_scale * drho_dT_dT
    drho_dS_dP = rho_scale * drho_dS_dP
    drho_dT_dP = rho_scale * drho_dT_dP
  endif

  if (EOS%RL2_T2_to_Pa /= 1.0) then
    drho_dS_dP = EOS%RL2_T2_to_Pa * drho_dS_dP
    drho_dT_dP = EOS%RL2_T2_to_Pa * drho_dT_dP
  endif

  if (EOS%C_to_degC /= 1.0) then
    drho_dS_dT = EOS%C_to_degC * drho_dS_dT
    drho_dT_dT = EOS%C_to_degC**2 * drho_dT_dT
    drho_dT_dP = EOS%C_to_degC * drho_dT_dP
  endif

  if (EOS%S_to_ppt /= 1.0) then
    drho_dS_dS = EOS%S_to_ppt**2 * drho_dS_dS
    drho_dS_dT = EOS%S_to_ppt * drho_dS_dT
    drho_dS_dP = EOS%S_to_ppt * drho_dS_dP
  endif

end procedure calculate_density_second_derivs_scalar
module procedure calculate_spec_vol_derivs_array
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
              "calculate_spec_vol_derivs_array: EOS%form_of_EOS is not valid.")

  call EOS%type%calculate_specvol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, start, npts)

end procedure calculate_spec_vol_derivs_array
module procedure calc_spec_vol_derivs_1d
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]
  real :: dSVdT_scale ! A factor to convert dSV_dT to the desired units [kg degC R-1 C-1 m-3 ~> 1]
  real :: dSVdS_scale ! A factor to convert dSV_dS to the desired units [kg ppt R-1 S-1 m-3 ~> 1]
  integer :: i, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(dSV_dT) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_spec_vol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, is, npts, EOS)
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    call calculate_spec_vol_derivs_array(Ta, Sa, pres, dSV_dT, dSV_dS, is, npts, EOS)
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  dSVdT_scale = spv_scale * EOS%C_to_degC
  dSVdS_scale = spv_scale * EOS%S_to_ppt
  if ((dSVdT_scale /= 1.0) .or. (dSVdS_scale /= 1.0)) then ; do i=is,ie
    dSV_dT(i) = dSVdT_scale * dSV_dT(i)
    dSV_dS(i) = dSVdS_scale * dSV_dS(i)
  enddo ; endif

end procedure calc_spec_vol_derivs_1d
module procedure calculate_compress_1d
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  integer :: i, is, ie, npts
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
              "calculate_compress_1d: EOS%form_of_EOS is not valid.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  do i=is,ie
    pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
    Ta(i) = EOS%C_to_degC * T(i)
    Sa(i) = EOS%S_to_ppt * S(i)
  enddo

  call EOS%type%calculate_compress_array(Ta, Sa, pres, rho, drho_dp, is, npts)

  if (EOS%kg_m3_to_R /= 1.0) then ; do i=is,ie
    rho(i) = EOS%kg_m3_to_R * rho(i)
  enddo ; endif
  if (EOS%L_T_to_m_s /= 1.0) then ; do i=is,ie
    drho_dp(i) = EOS%L_T_to_m_s**2 * drho_dp(i)
  enddo ; endif

end procedure calculate_compress_1d
module procedure calculate_compress_scalar
  real, dimension(1) :: pa    ! Pressure in a size-1 1d array [R L2 T-2 ~> Pa]
  real, dimension(1) :: Ta    ! Temperature in a size-1 1d array [C ~> degC]
  real, dimension(1) :: Sa    ! Salinity in a size-1 1d array [S ~> ppt]
  real, dimension(1) :: rhoa  ! In situ density in a size-1 1d array [R ~> kg m-3]
  real, dimension(1) :: drho_dpa ! The partial derivative of density with pressure (also the
  Ta(1) = T ; Sa(1) = S ; pa(1) = pressure

  call calculate_compress_1d(Ta, Sa, pa, rhoa, drho_dpa, EOS)
  rho = rhoa(1) ; drho_dp = drho_dpa(1)

end procedure calculate_compress_scalar
module procedure average_specific_vol
  real, dimension(size(T)) :: pres  ! Layer-top pressure converted to [Pa]
  real, dimension(size(T)) :: dpres ! Pressure change converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: T5(5)      ! Temperatures at five quadrature points [C ~> degC]
  real :: S5(5)      ! Salinities at five quadrature points [S ~> ppt]
  real :: p5(5)      ! Pressures at five quadrature points [R L2 T-2 ~> Pa]
  real :: a5(5)      ! Specific volumes at five quadrature points [R-1 ~> m3 kg-1]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  integer :: i, n, is, ie, npts
  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(T) ; npts = 1 + ie - is
  endif

  if (EOS%EOS_quadrature) then
    do i=is,ie
      do n=1,5
        T5(n) = T(i) ; S5(n) = S(i)
        p5(n) = p_t(i) + 0.25*real(5-n)*dp(i)
      enddo
      call calculate_spec_vol(T5, S5, p5, a5, EOS)

      ! Use Boole's rule to estimate the average specific volume.
      SpV_avg(i) = C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + 12.0*a5(3))
    enddo
  elseif ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call avg_spec_vol_linear(T, S, p_t, dp, SpV_avg, is, npts, EOS%Rho_T0_S0, &
                                 EOS%dRho_dT, EOS%dRho_dS, EOS%dRho_dp)
      case (EOS_WRIGHT)
        call avg_spec_vol_buggy_wright(T, S, p_t, dp, SpV_avg, is, npts)
      case (EOS_WRIGHT_FULL)
        call avg_spec_vol_wright_full(T, S, p_t, dp, SpV_avg, is, npts)
      case (EOS_WRIGHT_REDUCED)
        call avg_spec_vol_wright_red(T, S, p_t, dp, SpV_avg, is, npts)
      case default
        call MOM_error(FATAL, "No analytic average specific volume option is available with this EOS!")
    end select
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * p_t(i)
      dpres(i) = EOS%RL2_T2_to_Pa * dp(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call avg_spec_vol_linear(Ta, Sa, pres, dpres, SpV_avg, is, npts, EOS%Rho_T0_S0, &
                                 EOS%dRho_dT, EOS%dRho_dS, EOS%dRho_dp)
      case (EOS_WRIGHT)
        call avg_spec_vol_buggy_wright(Ta, Sa, pres, dpres, SpV_avg, is, npts)
      case (EOS_WRIGHT_FULL)
        call avg_spec_vol_wright_full(Ta, Sa, pres, dpres, SpV_avg, is, npts)
      case (EOS_WRIGHT_REDUCED)
        call avg_spec_vol_wright_red(Ta, Sa, pres, dpres, SpV_avg, is, npts)
      case default
        call MOM_error(FATAL, "No analytic average specific volume option is available with this EOS!")
    end select
  endif

  spv_scale = EOS%R_to_kg_m3
  if (EOS%EOS_quadrature) spv_scale = 1.0
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then ; do i=is,ie
    SpV_avg(i) = spv_scale * SpV_avg(i)
  enddo ; endif

end procedure average_specific_vol
module procedure EoS_fit_range
  if (.not. allocated(EOS%type)) call MOM_error(FATAL, &
                  "calculate_compress: EOS%form_of_EOS is not valid.")

  call EOS%type%EoS_fit_range(T_min, T_max, S_min, S_max, p_min, p_max)

end procedure EoS_fit_range
module procedure EOS_domain
  integer :: halo_sz
  halo_sz = 0 ; if (present(halo)) halo_sz = halo

  EOSdom(1) = HI%isc - (HI%isd-1) - halo_sz
  EOSdom(2) = HI%iec - (HI%isd-1) + halo_sz

end procedure EOS_domain
module procedure analytic_int_specific_vol_dp
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  real :: dRdp_scale ! A factor to convert drho_dp to the desired units [T-2 L2 s2 m-2 ~> 1]
  if (EOS%EOS_quadrature) call MOM_error(FATAL, "EOS_quadrature is set!")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      dRdT_scale = EOS%kg_m3_to_R * EOS%C_to_degC
      dRdS_scale = EOS%kg_m3_to_R * EOS%S_to_ppt
      dRdp_scale = EOS%kg_m3_to_R * EOS%RL2_T2_to_Pa
      call int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, HI, EOS%kg_m3_to_R*EOS%Rho_T0_S0, &
                          dRdT_scale*EOS%dRho_dT, dRdS_scale*EOS%dRho_dS, dRdp_scale*EOS%dRho_dp, &
                          dza, intp_dza, intx_dza, inty_dza, halo_size, &
                          bathyP, P_surf, dP_tiny, MassWghtInterp)
    case (EOS_WRIGHT)
      call int_spec_vol_dp_wright(T, S, p_t, p_b, alpha_ref, HI, dza, intp_dza, intx_dza, &
                                  inty_dza, halo_size, bathyP, P_surf, dP_tiny, MassWghtInterp, &
                                  SV_scale=EOS%R_to_kg_m3, pres_scale=EOS%RL2_T2_to_Pa, &
                                  temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt)
    case (EOS_WRIGHT_FULL)
      call int_spec_vol_dp_wright_full(T, S, p_t, p_b, alpha_ref, HI, dza, intp_dza, intx_dza, &
                                  inty_dza, halo_size, bathyP, P_surf, dP_tiny, MassWghtInterp, &
                                  SV_scale=EOS%R_to_kg_m3, pres_scale=EOS%RL2_T2_to_Pa, &
                                  temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt)
    case (EOS_WRIGHT_REDUCED)
      call int_spec_vol_dp_wright_red(T, S, p_t, p_b, alpha_ref, HI, dza, intp_dza, intx_dza, &
                                  inty_dza, halo_size, bathyP, P_surf, dP_tiny, MassWghtInterp, &
                                  SV_scale=EOS%R_to_kg_m3, pres_scale=EOS%RL2_T2_to_Pa, &
                                  temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt)
    case default
      call MOM_error(FATAL, "No analytic integration option is available with this EOS!")
  end select

end procedure analytic_int_specific_vol_dp
module procedure analytic_int_density_dz
  real :: rho_scale  ! A multiplicative factor by which to scale density from kg m-3 to the
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  real :: dRdp_scale ! A factor to convert drho_dp to the desired units [T-2 L2 s2 m-2 ~> 1]
  real :: pres_scale ! A multiplicative factor to convert pressure into Pa [Pa T2 R-1 L-2 ~> 1]
  if (EOS%EOS_quadrature) call MOM_error(FATAL, "EOS_quadrature is set!")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      rho_scale = EOS%kg_m3_to_R
      dRdT_scale = EOS%kg_m3_to_R * EOS%C_to_degC
      dRdS_scale = EOS%kg_m3_to_R * EOS%S_to_ppt
      dRdp_scale = EOS%kg_m3_to_R * EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0) .or. (dRdp_scale /= 1.0)) then
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, rho_scale*EOS%Rho_T0_S0, &
                         dRdT_scale*EOS%dRho_dT, dRdS_scale*EOS%dRho_dS, dRdp_scale*EOS%dRho_dp, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, dz_neglect, MassWghtInterp, Z_0p=Z_0p)
      else
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                         EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, EOS%dRho_dp, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, dz_neglect, MassWghtInterp, Z_0p=Z_0p)
      endif
    case (EOS_WRIGHT)
      rho_scale = EOS%kg_m3_to_R
      pres_scale = EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (pres_scale /= 1.0) .or. (EOS%C_to_degC /= 1.0) .or. (EOS%S_to_ppt /= 1.0)) then
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, rho_scale, pres_scale, &
                                   temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt, Z_0p=Z_0p)
      else
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, Z_0p=Z_0p)
      endif
    case (EOS_WRIGHT_FULL)
      rho_scale = EOS%kg_m3_to_R
      pres_scale = EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (pres_scale /= 1.0) .or. (EOS%C_to_degC /= 1.0) .or. (EOS%S_to_ppt /= 1.0)) then
        call int_density_dz_wright_full(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, rho_scale, pres_scale, &
                                   temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt, Z_0p=Z_0p)
      else
        call int_density_dz_wright_full(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, Z_0p=Z_0p)
      endif
    case (EOS_WRIGHT_REDUCED)
      rho_scale = EOS%kg_m3_to_R
      pres_scale = EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (pres_scale /= 1.0) .or. (EOS%C_to_degC /= 1.0) .or. (EOS%S_to_ppt /= 1.0)) then
        call int_density_dz_wright_red(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, rho_scale, pres_scale, &
                                   temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt, Z_0p=Z_0p)
      else
        call int_density_dz_wright_red(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                   dz_neglect, MassWghtInterp, Z_0p=Z_0p)
      endif
    case default
      call MOM_error(FATAL, "No analytic integration option is available with this EOS!")
  end select

end procedure analytic_int_density_dz
module procedure query_compressible
  query_compressible = EOS%compressible
end procedure query_compressible
module procedure get_EOS_name
  select case (id)
    case (EOS_LINEAR)
      eos_name = EOS_LINEAR_STRING
    case (EOS_UNESCO)
      eos_name = EOS_UNESCO_STRING
    case (EOS_WRIGHT)
      eos_name = EOS_WRIGHT_STRING
    case (EOS_WRIGHT_REDUCED)
      eos_name = EOS_WRIGHT_RED_STRING
    case (EOS_WRIGHT_FULL)
      eos_name = EOS_WRIGHT_FULL_STRING
    case (EOS_TEOS10)
      eos_name = EOS_TEOS10_STRING
    case (EOS_ROQUET_RHO)
      eos_name = EOS_ROQUET_RHO_STRING
    case (EOS_ROQUET_SPV)
      eos_name = EOS_ROQUET_SPV_STRING
    case (EOS_JACKETT06)
      eos_name = EOS_JACKETT06_STRING
    case default
      call MOM_error(FATAL, "get_EOS_name: something went wrong internally - enumeration is not valid.")
  end select

end procedure get_EOS_name
module procedure EOS_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_EOS" ! This module's name.
  character(len=12)  :: TFREEZE_DEFAULT ! The default freezing point expression
  character(len=40)  :: tmpstr
  logical :: EOS_quad_default, EOS_TS_default
  real :: Rho_Tref_Sref ! Density at Tref degC and Sref ppt [kg m-3]
  real :: Tref          ! Reference temperature [degC]
  real :: Sref          ! Reference salinity [psu]
  real :: pref          ! Reference pressure [Pa]
  real :: rho0          ! Density at T=0, S=0 and p=0 [kg m-3]
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "EQN_OF_STATE", tmpstr, &
                 "EQN_OF_STATE determines which ocean equation of state should be used.  "//&
                 'Currently, the valid choices are "LINEAR", "UNESCO", "JACKETT_MCD", '//&
                 '"WRIGHT", "WRIGHT_REDUCED", "WRIGHT_FULL", "NEMO", "ROQUET_RHO", "ROQUET_SPV" '//&
                 'and "TEOS10".  This is only used if USE_EOS is true.', default=EOS_DEFAULT)
  select case (uppercase(tmpstr))
    case (EOS_LINEAR_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_LINEAR)
    case (EOS_UNESCO_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_UNESCO)
    case (EOS_JACKETT_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_UNESCO)
    case (EOS_WRIGHT_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_WRIGHT)
    case (EOS_WRIGHT_RED_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_WRIGHT_REDUCED)
    case (EOS_WRIGHT_FULL_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_WRIGHT_FULL)
    case (EOS_TEOS10_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_TEOS10)
    case (EOS_NEMO_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_ROQUET_RHO)
    case (EOS_ROQUET_RHO_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_ROQUET_RHO)
    case (EOS_ROQUET_SPV_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_ROQUET_SPV)
    case (EOS_JACKETT06_STRING)
      call EOS_manual_init(EOS, form_of_EOS=EOS_JACKETT06)
    case default
      call MOM_error(FATAL, "interpret_eos_selection: EQN_OF_STATE "//&
                              trim(tmpstr) // " in input file is invalid.")
  end select
  call MOM_mesg('interpret_eos_selection: equation of state set to "' // &
                trim(tmpstr)//'"', 5)

  if (EOS%form_of_EOS == EOS_LINEAR) then
    ! RHO(T,S) = RHO_REF + DRHO_DT*(T-T_REF) + DRHO_DS*(S-S_REF) + DRHO_DP*(P-P_REF)
    !          = RHO_REF - (DRHO_DT*T_REF + DRHO_DS*SREF + DRHO_DP*PREF) + (DRHO_DT*T + DRHO_DS*S + DRHO_DP*P)
    !          = RHO_T0_S0 + (DRHO_DT*T + DRHO_DS*S + DRHO_DP*P)
    call get_param(param_file, mdl, "RHO_REF_LINEAR_EOS", Rho_Tref_Sref, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is the density "//&
                   "at T=T_REF_LINEAR_EOS, S=S_REF_LINEAR_EOS and p=P_REF_LINEAR_EOS", &
                   units="kg m-3", default=1000.0, old_name="RHO_TREF_SREF")
    call get_param(param_file, mdl, "T_REF_LINEAR_EOS", Tref, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is the reference "//&
                   "temperature.", units="degC", default=0.0, old_name="TREF")
    call get_param(param_file, mdl, "S_REF_LINEAR_EOS", Sref, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is the reference "//&
                   "salinity.", units="psu", default=0.0, old_name="SREF")
    call get_param(param_file, mdl, "P_REF_LINEAR_EOS", pref, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is the reference "//&
                   "pressure.", units="Pa", default=0.0)
    call get_param(param_file, mdl, "DRHO_DT", EOS%dRho_dT, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is "//&
                   "the partial derivative of density with temperature.", &
                   units="kg m-3 K-1", default=-0.2)
    call get_param(param_file, mdl, "DRHO_DS", EOS%dRho_dS, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is "//&
                   "the partial derivative of density with salinity.", &
                   units="kg m-3 ppt-1", default=0.8)
    call get_param(param_file, mdl, "DRHO_DP", EOS%dRho_dp, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is "//&
                   "the partial derivative of density with pressure (the inverse of "//&
                   "sound speed squared).", units="s2 m-2", default=0.0)
    rho0 = Rho_Tref_Sref - ((EOS%dRho_dT * Tref + EOS%dRho_dS * Sref) + EOS%dRho_dp * pref)
    call get_param(param_file, mdl, "RHO_T0_S0", EOS%Rho_T0_S0, &
                   "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", this is the density "//&
                   "at T=0, S=0 and p=0.  If RHO_TO_SO is specified, RHO_REF_LINEAR_EOS, "//&
                   "T_REF_LINEAR_EOS, S_REF_LINEAR_EOS and P_REF_LINEAR_EOS are not used.", &
                   units="kg m-3", default=rho0)
    EOS%Compressible = (EOS%dRho_dp/=0.0)
    call EOS_manual_init(EOS, form_of_EOS=EOS_LINEAR, Rho_T0_S0=EOS%Rho_T0_S0, &
                         dRho_dT=EOS%dRho_dT, dRho_dS=EOS%dRho_dS, dRho_dp=EOS%dRho_dp)
  endif
  if (EOS%form_of_EOS == EOS_WRIGHT) then
    call get_param(param_file, mdl, "USE_WRIGHT_2ND_DERIV_BUG", EOS%use_Wright_2nd_deriv_bug, &
                 "If true, use a bug in the calculation of the second derivatives of density "//&
                 "with temperature and with temperature and pressure that causes some terms "//&
                 "to be only 2/3 of what they should be.", default=.false.)
    call EOS_manual_init(EOS, form_of_EOS=EOS_WRIGHT, use_Wright_2nd_deriv_bug=EOS%use_Wright_2nd_deriv_bug)
  endif

  if (present(use_conT_absS)) then
    EOS%use_conT_absS = use_conT_absS
  else
    EOS%use_conT_absS = .false. ! Assuming it is not needed, it is set to false
  endif

  EOS_quad_default = .not.((EOS%form_of_EOS == EOS_LINEAR) .or. &
                           (EOS%form_of_EOS == EOS_WRIGHT) .or. &
                           (EOS%form_of_EOS == EOS_WRIGHT_REDUCED) .or. &
                           (EOS%form_of_EOS == EOS_WRIGHT_FULL))
  call get_param(param_file, mdl, "EOS_QUADRATURE", EOS%EOS_quadrature, &
                 "If true, always use the generic (quadrature) code "//&
                 "code for the integrals of density.", default=EOS_quad_default)

  TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING
  if ((EOS%form_of_EOS == EOS_TEOS10 .or. EOS%form_of_EOS == EOS_ROQUET_RHO .or. &
       EOS%form_of_EOS == EOS_ROQUET_SPV)) &
    TFREEZE_DEFAULT = TFREEZE_TEOS10_STRING
  call get_param(param_file, mdl, "TFREEZE_FORM", tmpstr, &
                 "TFREEZE_FORM determines which expression should be "//&
                 "used for the freezing point.  Currently, the valid "//&
                 'choices are "LINEAR", "MILLERO_78", "TEOS_POLY", "TEOS10"', &
                 default=TFREEZE_DEFAULT)
  select case (uppercase(tmpstr))
    case (TFREEZE_LINEAR_STRING)
      EOS%form_of_TFreeze = TFREEZE_LINEAR
    case (TFREEZE_MILLERO_STRING)
      EOS%form_of_TFreeze = TFREEZE_MILLERO
    case (TFREEZE_TEOSPOLY_STRING)
      EOS%form_of_TFreeze = TFREEZE_TEOSPOLY
    case (TFREEZE_TEOS10_STRING)
      EOS%form_of_TFreeze = TFREEZE_TEOS10
    case default
      call MOM_error(FATAL, "interpret_eos_selection:  TFREEZE_FORM "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select

  if (EOS%form_of_TFreeze == TFREEZE_LINEAR) then
    call get_param(param_file, mdl, "TFREEZE_S0_P0",EOS%TFr_S0_P0, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the freezing potential temperature at "//&
                 "S=0, P=0.", units="degC", default=0.0)
    call get_param(param_file, mdl, "DTFREEZE_DS",EOS%dTFr_dS, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with salinity.", &
                 units="degC ppt-1", default=-0.054)
    call get_param(param_file, mdl, "DTFREEZE_DP",EOS%dTFr_dP, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with pressure.", &
                 units="degC Pa-1", default=0.0)
  endif

  if ((EOS%form_of_TFreeze==TFREEZE_TEOSPOLY) .or. (EOS%form_of_TFreeze==TFREEZE_TEOS10)) then
     ! Which default is appropriate for Millero?
     EOS_TS_default = .false.
  else
     EOS_TS_default = .true.
  endif
  call get_param(param_file, mdl, "TFREEZE_S_IS_PRACS", EOS%TFreeze_S_is_pracS, &
               "When True, the model will check if the model internal salinity is "//&
               "practical salinity.  If the model uses absolute salinity, a "//&
               "conversion will be applied.", default=EOS_TS_default)
  call get_param(param_file, mdl, "TFREEZE_T_IS_POTT", EOS%TFreeze_T_is_potT, &
               "When True, the model will check if the model internal temperature is "//&
               "potential temperature.  If the model uses conservative temperature, a "//&
               "conversion will be applied.", default=EOS_TS_default)

  if ((EOS%form_of_EOS == EOS_TEOS10 .or. EOS%form_of_EOS == EOS_ROQUET_RHO .or. &
       EOS%form_of_EOS == EOS_ROQUET_SPV) .and. &
      .not.((EOS%form_of_TFreeze == TFREEZE_TEOS10) .or. (EOS%form_of_TFreeze == TFREEZE_TEOSPOLY)) ) then
    call MOM_error(WARNING, "interpret_eos_selection:  EOS_TEOS10 or EOS_ROQUET_RHO or EOS_ROQUET_SPV "//&
                   "should only be used along with TFREEZE_FORM = TFREEZE_TEOS10 or TFREEZE_TEOSPOLY.")
  endif

  ! Unit conversions
  EOS%m_to_Z = 1. ; if (present(US)) EOS%m_to_Z = US%m_to_Z
  EOS%kg_m3_to_R = 1. ; if (present(US)) EOS%kg_m3_to_R = US%kg_m3_to_R
  EOS%R_to_kg_m3 = 1. ; if (present(US)) EOS%R_to_kg_m3 = US%R_to_kg_m3
  EOS%RL2_T2_to_Pa = 1. ; if (present(US)) EOS%RL2_T2_to_Pa = US%RL2_T2_to_Pa
  EOS%L_T_to_m_s = 1. ; if (present(US)) EOS%L_T_to_m_s = US%L_T_to_m_s
  EOS%degC_to_C = 1. ; if (present(US)) EOS%degC_to_C = US%degC_to_C
  EOS%C_to_degC = 1. ; if (present(US)) EOS%C_to_degC = US%C_to_degC
  EOS%ppt_to_S = 1. ; if (present(US)) EOS%ppt_to_S = US%ppt_to_S
  EOS%S_to_ppt = 1. ; if (present(US)) EOS%S_to_ppt = US%S_to_ppt

end procedure EOS_init
module procedure EOS_manual_init
  if (present(form_of_EOS)) then
    EOS%form_of_EOS     = form_of_EOS
    if (allocated(EOS%type)) deallocate(EOS%type) ! Needed during testing which re-initializes
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        allocate(linear_EOS :: EOS%type)
      case (EOS_UNESCO)
        allocate(UNESCO_EOS :: EOS%type)
      case (EOS_WRIGHT)
        allocate(buggy_Wright_EOS :: EOS%type)
      case (EOS_WRIGHT_FULL)
        allocate(Wright_full_EOS :: EOS%type)
      case (EOS_WRIGHT_REDUCED)
        allocate(Wright_red_EOS :: EOS%type)
      case (EOS_JACKETT06)
        allocate(Jackett06_EOS :: EOS%type)
      case (EOS_TEOS10)
        allocate(TEOS10_EOS :: EOS%type)
      case (EOS_ROQUET_RHO)
        allocate(Roquet_rho_EOS :: EOS%type)
      case (EOS_ROQUET_SPV)
        allocate(Roquet_SpV_EOS :: EOS%type)
    end select
    select type (t => EOS%type)
      type is (linear_EOS)
        call t%set_params_linear(Rho_T0_S0, dRho_dT, dRho_dS, dRho_dp)
      type is (buggy_Wright_EOS)
        call t%set_params_buggy_Wright(use_Wright_2nd_deriv_bug)
    end select
  endif
  if (present(form_of_TFreeze))  EOS%form_of_TFreeze = form_of_TFreeze
  if (present(EOS_quadrature ))  EOS%EOS_quadrature  = EOS_quadrature
  if (present(Compressible   ))  EOS%Compressible    = Compressible
  if (present(Rho_T0_S0      ))  EOS%Rho_T0_S0       = Rho_T0_S0
  if (present(drho_dT        ))  EOS%drho_dT         = drho_dT
  if (present(dRho_dS        ))  EOS%dRho_dS         = dRho_dS
  if (present(dRho_dp        ))  EOS%dRho_dp         = dRho_dp
  if (present(TFr_S0_P0      ))  EOS%TFr_S0_P0       = TFr_S0_P0
  if (present(dTFr_dS        ))  EOS%dTFr_dS         = dTFr_dS
  if (present(dTFr_dp        ))  EOS%dTFr_dp         = dTFr_dp
  if (present(use_Wright_2nd_deriv_bug)) EOS%use_Wright_2nd_deriv_bug = use_Wright_2nd_deriv_bug

end procedure EOS_manual_init
module procedure convert_temp_salt_for_TEOS10
  real, parameter :: Sref_Sprac = (35.16504/35.0) ! The TEOS 10 conversion factor to go from
  integer :: i, j, k
  if ((EOS%form_of_EOS /= EOS_TEOS10) .and. (EOS%form_of_EOS /= EOS_ROQUET_RHO) .and. &
      (EOS%form_of_EOS /= EOS_ROQUET_SPV)) return

  do k=1,kd ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    if (mask_z(i,j,k) >= 1.0) then
      S(i,j,k) = Sref_Sprac * S(i,j,k)
      T(i,j,k) = EOS%degC_to_C*poTemp_to_consTemp(EOS%C_to_degC*T(i,j,k), EOS%S_to_ppt*S(i,j,k))
    endif
  enddo ; enddo ; enddo
end procedure convert_temp_salt_for_TEOS10
module procedure cons_temp_to_pot_temp
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(S)) :: Sa    ! Salinity converted to [ppt]
  real :: T_scale ! A factor to convert potential temperature from degC to the desired units [C degC-1 ~> 1]
  integer :: i, is, ie
  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(T)
  endif

  if ((EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    poTemp(is:ie) = consTemp_to_poTemp(T(is:ie), S(is:ie))
  else
    do i=is,ie
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    poTemp(is:ie) = consTemp_to_poTemp(Ta(is:ie), Sa(is:ie))
  endif

  T_scale = EOS%degC_to_C
  if (present(scale)) T_scale = scale
  if (T_scale /= 1.0) then ; do i=is,ie
    poTemp(i) = T_scale * poTemp(i)
  enddo ; endif

end procedure cons_temp_to_pot_temp
module procedure pot_temp_to_cons_temp
  real, dimension(size(T)) :: Tp    ! Potential temperature converted to [degC]
  real, dimension(size(S)) :: Sa    ! Absolute salinity converted to [ppt]
  real :: T_scale ! A factor to convert potential temperature from degC to the desired units [C degC-1 ~> 1]
  integer :: i, is, ie
  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(T)
  endif


  if ((EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    consTemp(is:ie) = poTemp_to_consTemp(T(is:ie), S(is:ie))
  else
    do i=is,ie
      Tp(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    consTemp(is:ie) = poTemp_to_consTemp(Tp(is:ie), Sa(is:ie))
  endif

  T_scale = EOS%degC_to_C
  if (present(scale)) T_scale = scale
  if (T_scale /= 1.0) then ; do i=is,ie
    consTemp(i) = T_scale * consTemp(i)
  enddo ; endif

end procedure pot_temp_to_cons_temp
module procedure abs_saln_to_prac_saln
  real :: S_scale ! A factor to convert practical salinity from ppt to the desired units [S PSU-1 ~> 1]
  real, parameter :: Sprac_Sref = (35.0/35.16504) ! The TEOS 10 conversion factor to go from
  integer :: i, is, ie
  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(S)
  endif

  if (present(scale)) then
    S_scale = Sprac_Sref * scale
    do i=is,ie
      prSaln(i) = S_scale * S(i)
    enddo
  else
    do i=is,ie
      prSaln(i) = Sprac_Sref * S(i)
    enddo
  endif

end procedure abs_saln_to_prac_saln
module procedure prac_saln_to_abs_saln
  real :: S_scale ! A factor to convert absolute salinity from ppt to the desired units [S ppt-1 ~> 1]
  real, parameter :: Sref_Sprac = (35.16504/35.0) ! The TEOS 10 conversion factor to go from
  integer :: i, is, ie
  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(S)
  endif

  if (present(scale)) then
    S_scale = Sref_Sprac * scale
    do i=is,ie
      absSaln(i) = S_scale * S(i)
    enddo
  else
    do i=is,ie
      absSaln(i) = Sref_Sprac * S(i)
    enddo
  endif

end procedure prac_saln_to_abs_saln
module procedure EOS_quadrature
  EOS_quadrature  = EOS%EOS_quadrature

end procedure EOS_quadrature
module procedure EOS_unit_tests
  type(EOS_type) :: EOS_tmp
  logical :: fail
  if (verbose) write(stdout,*) '==== MOM_EOS: EOS_unit_tests ===='
  EOS_unit_tests = .false. ! Normally return false

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_TEOS10)
  fail = test_TS_conversion_consistency(T_cons=9.989811727177308, S_abs=35.16504, &
                                        T_pot=10.0, S_prac=35.0, EOS=EOS_tmp, verbose=verbose)
  if (verbose .and. fail) call MOM_error(WARNING, "Some EOS variable conversions tests have failed.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_UNESCO)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "UNESCO", &
                              rho_check=1027.54345796120*EOS_tmp%kg_m3_to_R)
  if (verbose .and. fail) call MOM_error(WARNING, "UNESCO EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_WRIGHT_FULL)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "WRIGHT_FULL", &
                              rho_check=1027.55177447616*EOS_tmp%kg_m3_to_R, avg_Sv_check=.true.)
  if (verbose .and. fail) call MOM_error(WARNING, "WRIGHT_FULL EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_WRIGHT_REDUCED)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "WRIGHT_REDUCED", &
                              rho_check=1027.54303596346*EOS_tmp%kg_m3_to_R, avg_Sv_check=.true.)
  if (verbose .and. fail) call MOM_error(WARNING, "WRIGHT_REDUCED EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  ! This test is deliberately outside of the fit range for WRIGHT_REDUCED, and it results in the expected warnings.
  ! call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_WRIGHT_REDUCED)
  ! fail = test_EOS_consistency(25.0, 15.0, 1.0e7, EOS_tmp, verbose, "WRIGHT_REDUCED", &
  !                             rho_check=1012.625699301455*EOS_tmp%kg_m3_to_R)
  ! if (verbose .and. fail) call MOM_error(WARNING, "WRIGHT_REDUCED EOS has failed some self-consistency tests.")
  ! EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_WRIGHT, use_Wright_2nd_deriv_bug=.true.)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "WRIGHT", &
                              rho_check=1027.54303596346*EOS_tmp%kg_m3_to_R, avg_Sv_check=.true.)
  ! These last test is a known failure and since MPI is not necessarily initializaed when running these tests
  ! we need to avoid flagging the fails.
  !if (verbose .and. fail) call MOM_error(WARNING, "WRIGHT EOS has failed some self-consistency tests.")
  !EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_ROQUET_RHO)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "ROQUET_RHO", &
                              rho_check=1027.42385663668*EOS_tmp%kg_m3_to_R)
  if (verbose .and. fail) call MOM_error(WARNING, "ROQUET_RHO EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_ROQUET_SPV)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "ROQUET_SPV", &
                              rho_check=1027.42387475199*EOS_tmp%kg_m3_to_R)
  if (verbose .and. fail) call MOM_error(WARNING, "ROQUET_SPV EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_JACKETT06)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "JACKETT06", &
                              rho_check=1027.539690758425*EOS_tmp%kg_m3_to_R)
  if (verbose .and. fail) call MOM_error(WARNING, "JACKETT06 EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  ! The TEOS10 equation of state is not passing the self consistency tests for dho_dS_dp due
  ! to a bug (a missing division by the square root of offset-salinity) on line 111 of
  ! pkg/GSW-Fortan/toolbox/gsw_specvol_second_derivatives.f90.  This bug has been highlighted in an
  ! issue posted to the TEOS-10/GSW-Fortran page at github.com/TEOS-10/GSW-Fortran/issues/26, and
  ! it will be corrected by github.com/mom-ocean/GSW-Fortran/pull/1 .
  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_TEOS10)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "TEOS10", skip_2nd=.true., &
                              rho_check=1027.42355961492*EOS_tmp%kg_m3_to_R)
  if (verbose .and. fail) call MOM_error(WARNING, "TEOS10 EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_ROQUET_RHO)
  fail = test_EOS_consistency(10.0, 30.0, 1.0e7, EOS_tmp, verbose, "ROQUET_RHO", &
                              rho_check=1027.45140117152*EOS_tmp%kg_m3_to_R)
  ! The corresponding check value published by Roquet et al. (2015) is 1027.45140 [kg m-3].
  if (verbose .and. fail) call MOM_error(WARNING, "Roquet_rho EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_ROQUET_SPV)
  fail = test_EOS_consistency(10.0, 30.0, 1.0e7, EOS_tmp, verbose, "ROQUET_SPV", &
                              spv_check=9.73282046614623e-04*EOS_tmp%R_to_kg_m3)
  ! The corresponding check value here published by Roquet et al. (2015) is 9.732819628e-04 [m3 kg-1],
  ! but the order of arithmetic there was not completely specified with parentheses.
  if (verbose .and. fail) call MOM_error(WARNING, "ROQUET_SPV EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_EOS=EOS_LINEAR, Rho_T0_S0=1000.0, drho_dT=-0.2, dRho_dS=0.8, dRho_dp=5.0e-7)
  fail = test_EOS_consistency(25.0, 35.0, 1.0e7, EOS_tmp, verbose, "LINEAR", &
                              rho_check=1028.0*EOS_tmp%kg_m3_to_R, avg_Sv_check=.true.)
  if (verbose .and. fail) call MOM_error(WARNING, "LINEAR EOS has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  ! Test the freezing point calculations

  call EOS_manual_init(EOS_tmp, form_of_TFreeze=TFREEZE_LINEAR, TFr_S0_P0=0.0, dTFr_dS=-0.054, &
                       dTFr_dP=-7.6e-8)
  fail = test_TFr_consistency(35.0, 1.0e7, EOS_tmp, verbose, "LINEAR", TFr_check=-2.65*EOS_tmp%degC_to_C)
  if (verbose .and. fail) call MOM_error(WARNING, "LINEAR TFr has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_TFreeze=TFREEZE_MILLERO)
  fail = test_TFr_consistency(35.0, 1.0e7, EOS_tmp, verbose, "MILLERO_78", &
                              TFr_check=-2.69730134114106*EOS_tmp%degC_to_C)
  if (verbose .and. fail) call MOM_error(WARNING, "MILLERO_78 TFr has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_TFreeze=TFREEZE_TEOS10)
  fail = test_TFr_consistency(35.0, 1.0e7, EOS_tmp, verbose, "TEOS10", &
                              TFr_check=-2.69099996992861*EOS_tmp%degC_to_C)
  if (verbose .and. fail) call MOM_error(WARNING, "TEOS10 TFr has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  call EOS_manual_init(EOS_tmp, form_of_TFreeze=TFREEZE_TEOSPOLY)
  fail = test_TFr_consistency(35.0, 1.0e7, EOS_tmp, verbose, "TEOS_POLY", &
                              TFr_check=-2.691165259327735*EOS_tmp%degC_to_C)
  if (verbose .and. fail) call MOM_error(WARNING, "TEOS_POLY TFr has failed some self-consistency tests.")
  EOS_unit_tests = EOS_unit_tests .or. fail

  if (EOS_unit_tests) then
    call MOM_error(WARNING, "EOS_unit_tests: One or more EOS tests have failed!")
  else
    if (verbose) call MOM_mesg("EOS_unit_tests: All EOS consistency tests have passed.")
  endif

end procedure EOS_unit_tests
module procedure test_TS_conversion_consistency
  real :: Sabs(1)  ! Absolute or reference salinity [g kg-1]
  real :: Sprac(1) ! Practical salinity [PSU]
  real :: Stest(1) ! A converted salinity [ppt]
  real :: Tcons(1) ! Conservative temperature [degC]
  real :: Tpot(1)  ! Potential temperature [degC]
  real :: Ttest(1) ! A converted temperature [degC]
  real :: Stol     ! Roundoff error on a typical value of salinities [ppt]
  real :: Ttol     ! Roundoff error on a typical value of temperatures [degC]
  logical :: test_OK ! True if a particular test is consistent.
  logical :: OK      ! True if all checks so far are consistent.
  OK = .true.

  ! Copy scalar input values into the corresponding arrays
  Sabs(1) = S_abs ; Sprac(1) = S_prac ; Tcons(1) = T_cons ; Tpot(1) = T_pot

  ! Set tolerances for the conversions.
  Ttol = 2.0 * 400.0*epsilon(Ttol)
  Stol = 35.0 * 400.0*epsilon(Stol)

  ! Check that the converted salinities agree
  call abs_saln_to_prac_saln(Sabs, Stest, EOS)
  test_OK = (abs(Stest(1) - Sprac(1)) <= Stol)
  if (verbose) call write_check_msg("MOM6 Sprac", Stest(1), Sprac(1), Stol, test_OK)
  OK = OK .and. test_OK

  call prac_saln_to_abs_saln(Sprac, Stest, EOS)
  test_OK = (abs(Stest(1) - Sabs(1)) <= Stol)
  if (verbose) call write_check_msg("MOM6 Sabs", Stest(1), Sabs(1), Stol, test_OK)
  OK = OK .and. test_OK

  call cons_temp_to_pot_temp(Tcons, Sabs, Ttest, EOS)
  test_OK = (abs(Ttest(1) - Tpot(1)) <= Ttol)
  if (verbose) call write_check_msg("MOM6 Tpot", Ttest(1), Tpot(1), Ttol, test_OK)
  OK = OK .and. test_OK

  call pot_temp_to_cons_temp(Tpot, Sabs, Ttest, EOS)
  test_OK = (abs(Ttest(1) - Tcons(1)) <= Ttol)
  if (verbose) call write_check_msg("MOM6 Tcons", Ttest(1), Tcons(1), Ttol, test_OK)
  OK = OK .and. test_OK

  inconsistent = .not.OK
end procedure test_TS_conversion_consistency
module procedure test_TFr_consistency
  real, dimension(-3:3,-3:3) :: S ! Salinities at the test value and perturbed points [S ~> ppt]
  real, dimension(-3:3,-3:3) :: P ! Pressures at the test value and perturbed points [R L2 T-2 ~> Pa]
  real, dimension(-3:3,-3:3,2) :: TFr ! Freezing point at the test value and perturbed points [C ~> degC]
  real :: dS         ! Magnitude of salinity perturbations [S ~> ppt]
  real :: dp         ! Magnitude of pressure perturbations [R L2 T-2 ~> Pa]
  real :: TFr_tol    ! Roundoff error on a typical value of TFreeze [C ~> degC]
  logical :: test_OK ! True if a particular test is consistent.
  logical :: OK      ! True if all checks so far are consistent.
  integer :: i, j, n
  OK = .true.

  dS = 0.5*EOS%ppt_to_S      ! Salinity perturbations [S ~> ppt]
  dp = 10.0e4 / EOS%RL2_T2_to_Pa ! Pressure perturbations [R L2 T-2 ~> Pa]

  ! TEOS 10 requires a tolerance that is ~20 times larger than other freezing point
  ! expressions because it lacks parentheses.
  TFr_tol = 2.0*EOS%degC_to_C * 400.0*epsilon(TFr_tol)

  do n=1,2
    ! Calculate density values with a wide enough stencil to estimate first and second derivatives
    ! with up to 6th order accuracy.  Doing this twice with different sizes of perturbations allows
    ! the evaluation of whether the finite differences are converging to the calculated values at a
    ! rate that is consistent with the order of accuracy of the finite difference forms, and hence
    ! the consistency of the calculated values.
    do j=-3,3 ; do i=-3,3
      S(i,j) = max(S_test + n*dS*i, 0.0)
      p(i,j) = max(p_test + n*dp*j, 0.0)
    enddo ; enddo
    do j=-3,3
      call calculate_TFreeze(S(:,j), p(:,j), TFr(:,j,n), EOS)
    enddo
  enddo

  ! Check that the freezing point agrees with the provided check value
  if (present(TFr_check)) then
    test_OK = (abs(TFr_check - TFr(0,0,1)) <= TFr_tol)
    OK = OK .and. test_OK
    if (verbose) call write_check_msg(trim(EOS_name)//" TFr", TFr(0,0,1), TFr_check, Tfr_tol, test_OK)
  endif

  inconsistent = .not.OK
end procedure test_TFr_consistency
module procedure write_check_msg
  character(len=200) :: mesg
  write(mesg, '(ES24.16," vs. ",ES24.16,", diff=",ES12.4,", tol=",ES12.4)') &
        val, val_chk, val-val_chk, val_tol
  if (test_OK) then
    write(stdout,*) trim(var_name)//" agrees with its check value :"//trim(mesg)
  else
    write(stderr,*) trim(var_name)//" disagrees with its check value :"//trim(mesg)
  endif
end procedure write_check_msg
module procedure test_EOS_consistency
  real, dimension(-3:3,-3:3,-3:3) :: T ! Temperatures at the test value and perturbed points [C ~> degC]
  real, dimension(-3:3,-3:3,-3:3) :: S ! Salinities at the test value and perturbed points [S ~> ppt]
  real, dimension(-3:3,-3:3,-3:3) :: P ! Pressures at the test value and perturbed points [R L2 T-2 ~> Pa]
  real, dimension(-3:3,-3:3,-3:3,2) :: rho ! Densities relative to rho_ref at the test value and
  real, dimension(-3:3,-3:3,-3:3,2) :: spv ! Specific volumes relative to spv_ref at the test value and
  real :: dT        ! Magnitude of temperature perturbations [C ~> degC]
  real :: dS        ! Magnitude of salinity perturbations [S ~> ppt]
  real :: dp        ! Magnitude of pressure perturbations [R L2 T-2 ~> Pa]
  real :: rho_ref   ! A reference density that is extracted for greater accuracy [R ~> kg m-3]
  real :: spv_ref   ! A reference specific volume that is extracted for greater accuracy [R-1 ~> m3 kg-1]
  real :: rho_nooff ! Density with no reference offset [R ~> kg m-3]
  real :: spv_nooff ! Specific volume with no reference offset [R-1 ~> m3 kg-1]
  real :: drho_dT   ! The partial derivative of density with potential
  real :: drho_dS   ! The partial derivative of density with salinity
  real :: drho_dp   ! The partial derivative of density with pressure (also the
  real :: dSV_dT(1) ! The partial derivative of specific volume with potential
  real :: dSV_dS(1) ! The partial derivative of specific volume with salinity
  real :: SpV_avg_a(1) ! The pressure-averaged specific volume determined analytically [R-1 ~> m3 kg-1]
  real :: SpV_avg_q(1) ! The pressure-averaged specific volume determined via quadrature [R-1 ~> m3 kg-1]
  real :: drho_dS_dS ! Second derivative of density with respect to S  [R S-2 ~> kg m-3 ppt-2]
  real :: drho_dS_dT ! Second derivative of density with respect to T and S [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
  real :: drho_dT_dT ! Second derivative of density with respect to T [R C-2 ~> kg m-3 degC-2]
  real :: drho_dS_dP ! Second derivative of density with respect to salinity and pressure
  real :: drho_dT_dP ! Second derivative of density with respect to temperature and pressure
  real :: drho_dT_fd(2) ! Two 6th order finite difference estimates of the partial derivative of density
  real :: drho_dS_fd(2) ! Two 6th order finite difference estimates of the partial derivative of density
  real :: drho_dp_fd(2) ! Two 6th order finite difference estimates of the partial derivative of density
  real :: dSV_dT_fd(2)  ! Two 6th order finite difference estimates of the partial derivative of
  real :: dSV_dS_fd(2)  ! Two 6th order finite difference estimates of the partial derivative of
  real :: drho_dS_dS_fd(2)  ! Two 6th order finite difference estimates of the second derivative of
  real :: drho_dS_dT_fd(2)  ! Two 6th order finite difference estimates of the second derivative of density
  real :: drho_dT_dT_fd(2)  ! Two 6th order finite difference estimates of the second derivative of
  real :: drho_dS_dP_fd(2)  ! Two 6th order finite difference estimates of the second derivative of density
  real :: drho_dT_dP_fd(2)  ! Two 6th order finite difference estimates of the second derivative of density
  real :: rho_tmp    ! A temporary copy of the situ density [R ~> kg m-3]
  real :: tol        ! The nondimensional tolerance from roundoff [nondim]
  real :: r_tol      ! Roundoff error on a typical value of density anomaly [R ~> kg m-3]
  real :: sv_tol     ! Roundoff error on a typical value of specific volume anomaly [R-1 ~> m3 kg-1]
  real :: tol_here   ! The tolerance for each check, in various units [various]
  real :: T_min, T_max ! The minimum and maximum temperature over which this EoS is fitted [degC]
  real :: S_min, S_max ! The minimum and maximum temperature over which this EoS is fitted [ppt]
  real :: p_min, p_max ! The minimum and maximum temperature over which this EoS is fitted [Pa]
  real :: count_fac  ! A factor in the roundoff estimates based on the factors in the numerator and
  real :: count_fac2 ! A factor in the roundoff estimates based on the factors in the numerator and
  character(len=200) :: mesg
  type(EOS_type) :: EOS_tmp
  logical :: test_OK ! True if a particular test is consistent.
  logical :: OK      ! True if all checks so far are consistent.
  logical :: test_2nd  ! If true, do tests on the 2nd derivative calculations
  logical :: test_avg_Sv ! If true, compare numerical and analytical estimates of the vertically
  integer :: order   ! The order of accuracy of the centered finite difference estimates (2, 4 or 6).
  integer :: i, j, k, n
  test_2nd = .true. ; if (present(skip_2nd)) test_2nd = .not.skip_2nd
  test_avg_Sv = .false. ; if (present(avg_Sv_check)) test_avg_Sv = avg_Sv_check

  dT = 0.1*EOS%degC_to_C     ! Temperature perturbations [C ~> degC]
  dS = 0.5*EOS%ppt_to_S      ! Salinity perturbations [S ~> ppt]
  dp = 10.0e4 / EOS%RL2_T2_to_Pa ! Pressure perturbations [R L2 T-2 ~> Pa]

  r_tol = 50.0*EOS%kg_m3_to_R * 10.*epsilon(r_tol)
  sv_tol = 5.0e-5*EOS%R_to_kg_m3 * 10.*epsilon(sv_tol)
  rho_ref = 1000.0*EOS%kg_m3_to_R
  spv_ref = 1.0 / rho_ref

  order = 4  ! This should be 2, 4 or 6.

  ! Check whether the consistency test is being applied outside of the value range of this EoS.
  call EoS_fit_range(EOS, T_min, T_max, S_min, S_max, p_min, p_max)
  if ((T_test < T_min) .or. (T_test > T_max)) then
    write(mesg, '(ES12.4," [degC] which is outside of the fit range of ",ES12.4," to ",ES12.4)') T_test, T_min, T_max
    call MOM_error(WARNING, trim(EOS_name)//" is being evaluated at a temperature of "//trim(mesg))
  endif
  if ((S_test < S_min) .or. (S_test > S_max)) then
    write(mesg, '(ES12.4," [ppt] which is outside of the fit range of ",ES12.4," to ",ES12.4)') S_test, S_min, S_max
    call MOM_error(WARNING, trim(EOS_name)//" is being evaluated at a salinity of "//trim(mesg))
  endif
  if ((p_test < p_min) .or. (p_test > p_max)) then
    write(mesg, '(ES12.4," [Pa] which is outside of the fit range of ",ES12.4," to ",ES12.4)') p_test, p_min, p_max
    call MOM_error(WARNING, trim(EOS_name)//" is being evaluated at a pressure of "//trim(mesg))
  endif

  do n=1,2
    ! Calculate density values with a wide enough stencil to estimate first and second derivatives
    ! with up to 6th order accuracy.  Doing this twice with different sizes of perturbations allows
    ! the evaluation of whether the finite differences are converging to the calculated values at a
    ! rate that is consistent with the order of accuracy of the finite difference forms, and hence
    ! the consistency of the calculated values.
    do k=-3,3 ; do j=-3,3 ; do i=-3,3
      T(i,j,k) = T_test + n*dT*i
      S(i,j,k) = S_test + n*dS*j
      p(i,j,k) = p_test + n*dp*k
    enddo ; enddo ; enddo
    do k=-3,3 ; do j=-3,3
      call calculate_density(T(:,j,k), S(:,j,k), p(:,j,k), rho(:,j,k,n), EOS, rho_ref=rho_ref)
      call calculate_spec_vol(T(:,j,k), S(:,j,k), p(:,j,k), spv(:,j,k,n), EOS, spv_ref=spv_ref)
    enddo ; enddo

    drho_dT_fd(n) = first_deriv(rho(:,0,0,n), n*dT, order)
    drho_dS_fd(n) = first_deriv(rho(0,:,0,n), n*dS, order)
    drho_dp_fd(n) = first_deriv(rho(0,0,:,n), n*dp, order)
    dSV_dT_fd(n) = first_deriv(spv(:,0,0,n), n*dT, order)
    dSV_dS_fd(n) = first_deriv(spv(0,:,0,n), n*dS, order)
    if (test_2nd) then
      drho_dT_dT_fd(n) = second_deriv(rho(:,0,0,n), n*dT, order)
      drho_dS_dS_fd(n) = second_deriv(rho(0,:,0,n), n*dS, order)
      drho_dS_dT_fd(n) = derivs_2d(rho(:,:,0,n), n**2*dT*dS, order)
      drho_dT_dP_fd(n) = derivs_2d(rho(:,0,:,n), n**2*dT*dP, order)
      drho_dS_dP_fd(n) = derivs_2d(rho(0,:,:,n), n**2*dS*dP, order)
    endif
  enddo

  call calculate_density_derivs(T(0,0,0), S(0,0,0), p(0,0,0), drho_dT, drho_dS, EOS)
  ! The first indices here are "0:0" because there is no scalar form of calculate_specific_vol_derivs.
  call calculate_specific_vol_derivs(T(0:0,0,0), S(0:0,0,0), p(0:0,0,0), dSV_dT, dSV_dS, EOS)
  if (test_2nd) &
    call calculate_density_second_derivs(T(0,0,0), S(0,0,0), p(0,0,0), &
                                       drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, EOS)
  call calculate_compress(T(0,0,0), S(0,0,0), p(0,0,0), rho_tmp, drho_dp, EOS)

  if (test_avg_Sv) then
    EOS_tmp = EOS
    call EOS_manual_init(EOS_tmp, EOS_quadrature=.false.)
    call average_specific_vol(T(0:0,0,0), S(0:0,0,0), p(0:0,0,0), p(0:0,0,0), SpV_avg_a, EOS_tmp)
    call EOS_manual_init(EOS_tmp, EOS_quadrature=.true.)
    call average_specific_vol(T(0:0,0,0), S(0:0,0,0), p(0:0,0,0), p(0:0,0,0), SpV_avg_q, EOS_tmp)
  endif

  OK = .true.

  tol = 1000.0*epsilon(tol)

  ! Check that the density agrees with the provided check value
  if (present(rho_check)) then
    test_OK = (abs(rho_check - (rho_ref + rho(0,0,0,1))) < tol*(rho_ref + rho(0,0,0,1)))
    OK = OK .and. test_OK
    if (verbose) &
      call write_check_msg(trim(EOS_name)//" rho", rho_ref+rho(0,0,0,1), rho_check, tol*rho(0,0,0,1), test_OK)
  endif

  ! Check that the specific volume agrees with the provided check value or the inverse of density
  if (present(spv_check)) then
    test_OK = (abs(spv_check - (spv_ref + spv(0,0,0,1))) < tol*abs(spv_ref + spv(0,0,0,1)))
    if (verbose) &
      call write_check_msg(trim(EOS_name)//" spv", spv_ref+spv(0,0,0,1), spv_check, tol*spv(0,0,0,1), test_OK)
    OK = OK .and. test_OK
  else
    test_OK = (abs((rho_ref+rho(0,0,0,1)) * (spv_ref + spv(0,0,0,1)) - 1.0) < tol)
    OK = OK .and. test_OK
    if (verbose) then
      write(mesg, '(ES16.8," and ",ES16.8,", ratio - 1 = ",ES16.8)') &
          rho_ref+rho(0,0,0,1), 1.0/(spv_ref + spv(0,0,0,1)), &
          (rho_ref+rho(0,0,0,1)) * (spv_ref + spv(0,0,0,1)) - 1.0
      if (test_OK) then
        write(stdout,*) "The values of "//trim(EOS_name)//" rho and 1/spv agree.  "//trim(mesg)
      else
        write(stderr,*) "The values of "//trim(EOS_name)//" rho and 1/spv disagree.  "//trim(mesg)
      endif
    endif
  endif

  ! Check that the densities are consistent when the reference value is extracted
  call calculate_density(T(0,0,0), S(0,0,0), p(0,0,0), rho_nooff, EOS)
  test_OK = (abs(rho_nooff - (rho_ref + rho(0,0,0,1))) < tol*rho_nooff)
  OK = OK .and. test_OK
  if (verbose .and. .not.test_OK) then
    write(mesg, '(ES24.16," vs. ",ES24.16," with tolerance ",ES12.4)') &
          rho_ref+rho(0,0,0,1), rho_nooff, tol*rho_nooff
    write(stderr,*) "For "//trim(EOS_name)//&
                   " rho with and without a reference value disagree: "//trim(mesg)
  endif

  ! Check that the specific volumes are consistent when the reference value is extracted
  call calculate_spec_vol(T(0,0,0), S(0,0,0), p(0,0,0), spv_nooff, EOS)
  test_OK = (abs(spv_nooff - (spv_ref + spv(0,0,0,1))) < tol*rho_nooff)
  OK = OK .and. test_OK
  if (verbose .and. .not.test_OK) then
    write(mesg, '(ES24.16," vs. ",ES24.16," with tolerance ",ES12.4)') &
          spv_ref + spv(0,0,0,1), spv_nooff, tol*spv_nooff
    write(stderr,*) "For "//trim(EOS_name)//&
                   " spv with and without a reference value disagree: "//trim(mesg)
  endif

  ! Account for the factors of terms in the numerator and denominator when estimating roundoff
  if (order == 6) then
    count_fac = 110.0/60.0 ; count_fac2 = 1088.0/180.0
  elseif (order == 4) then ! Use values appropriate for 4th order schemes.
    count_fac = 18.0/12.0 ; count_fac2 = 64.0/12.0
  else ! Use values appropriate for 2nd order schemes.
    count_fac = 2.0/2.0 ; count_fac2 = 4.0
  endif

  ! Check for the rate of convergence expected with a 4th or 6th order accurate discretization
  ! with a 20% margin of error and a tolerance for contributions from roundoff.
  tol_here = tol*abs(drho_dT) + count_fac*r_tol/dT
  OK = OK .and. check_FD(drho_dT, drho_dT_fd, tol_here, verbose, trim(EOS_name)//" drho_dT", order)
  tol_here = tol*abs(drho_dS) + count_fac*r_tol/dS
  OK = OK .and. check_FD(drho_dS, drho_dS_fd, tol_here, verbose, trim(EOS_name)//" drho_dS", order)
  tol_here = tol*abs(drho_dp) + count_fac*r_tol/dp
  OK = OK .and. check_FD(drho_dp, drho_dp_fd, tol_here, verbose, trim(EOS_name)//" drho_dp", order)
  tol_here = tol*abs(dSV_dT(1)) + count_fac*sv_tol/dT
  OK = OK .and. check_FD(dSV_dT(1), dSV_dT_fd, tol_here, verbose, trim(EOS_name)//" dSV_dT", order)
  tol_here = tol*abs(dSV_dS(1)) + count_fac*sv_tol/dS
  OK = OK .and. check_FD(dSV_dS(1), dSV_dS_fd, tol_here, verbose, trim(EOS_name)//" dSV_dS", order)
  if (test_2nd) then
    tol_here = tol*abs(drho_dT_dT) + count_fac2*r_tol/dT**2
    OK = OK .and. check_FD(drho_dT_dT, drho_dT_dT_fd, tol_here, verbose, trim(EOS_name)//" drho_dT_dT", order)
    ! The curvature in salinity is relatively weak, so looser tolerances are needed for some forms of EOS?
    tol_here = 10.0*(tol*abs(drho_dS_dS) + count_fac2*r_tol/dS**2)
    OK = OK .and. check_FD(drho_dS_dS, drho_dS_dS_fd, tol_here, verbose, trim(EOS_name)//" drho_dS_dS", order)
    tol_here = tol*abs(drho_dS_dT) + count_fac**2*r_tol/(dS*dT)
    OK = OK .and. check_FD(drho_dS_dT, drho_dS_dT_fd, tol_here, verbose, trim(EOS_name)//" drho_dS_dT", order)
    tol_here = tol*abs(drho_dT_dP) + count_fac**2*r_tol/(dT*dp)
    OK = OK .and. check_FD(drho_dT_dP, drho_dT_dP_fd, tol_here, verbose, trim(EOS_name)//" drho_dT_dP", order)
    tol_here = tol*abs(drho_dS_dP) + count_fac**2*r_tol/(dS*dp)
    OK = OK .and. check_FD(drho_dS_dP, drho_dS_dP_fd, tol_here, verbose, trim(EOS_name)//" drho_dS_dP", order)
  endif

  if (test_avg_Sv) then
    tol_here = 0.5*tol*(abs(SpV_avg_a(1)) + abs(SpV_avg_q(1)))
    test_OK = (abs(SpV_avg_a(1) - SpV_avg_q(1)) < tol_here)
    if (verbose) then
      write(mesg, '(ES24.16," and ",ES24.16," differ by ",ES16.8," (",ES10.2,"), tol=",ES16.8)') &
        SpV_avg_a(1), SpV_avg_q(1), SpV_avg_a(1) - SpV_avg_q(1), &
        2.0*(SpV_avg_a(1) - SpV_avg_q(1)) / (abs(SpV_avg_a(1)) + abs(SpV_avg_q(1)) + tiny(SpV_avg_a(1))), &
        tol_here
      if (verbose .and. .not.test_OK) then
        write(stderr,*) "The values of "//trim(EOS_name)//" SpV_avg disagree. "//trim(mesg)
      elseif (verbose) then
        write(stdout,*) "The values of "//trim(EOS_name)//" SpV_avg agree: "//trim(mesg)
      endif
    endif
    OK = OK .and. test_OK
  endif

  inconsistent = .not.OK

  contains

  !> Return a finite difference estimate of the first derivative of a field in arbitrary units [A B-1]
  real function first_deriv(R, dx, order)
    real,    intent(in) :: R(-3:3) !< The field whose derivative is being taken, in arbitrary units [A]
    real,    intent(in) :: dx      !< The spacing in parameter space, in different arbitrary units [B]
    integer, intent(in) :: order   !< The order of accuracy of the centered finite difference estimates (2, 4 or 6)

    if (order == 6) then  ! Find a 6th order accurate first derivative on a regular grid.
      first_deriv = (45.0*(R(1)-R(-1)) + (-9.0*(R(2)-R(-2)) + (R(3)-R(-3))) ) / (60.0 * dx)
    elseif (order == 4) then  ! Find a 4th order accurate first derivative on a regular grid.
      first_deriv = (8.0*(R(1)-R(-1)) - (R(2)-R(-2)) ) / (12.0 * dx)
    else  ! Find a 2nd order accurate first derivative on a regular grid.
      first_deriv = (R(1)-R(-1)) / (2.0 * dx)
    endif
  end function first_deriv

  !> Return a finite difference estimate of the second derivative of a field in arbitrary units [A B-2]
  real function second_deriv(R, dx, order)
    real,    intent(in) :: R(-3:3) !< The field whose derivative is being taken, in arbitrary units [A]
    real,    intent(in) :: dx      !< The spacing in parameter space, in different arbitrary units [B]
    integer, intent(in) :: order   !< The order of accuracy of the centered finite difference estimates (2, 4 or 6)

    if (order == 6) then  ! Find a 6th order accurate second derivative on a regular grid.
      second_deriv = ( -490.0*R(0) + (270.0*(R(1)+R(-1)) + (-27.0*(R(2)+R(-2)) + 2.0*(R(3)+R(-3))) )) / (180.0 * dx**2)
    elseif (order == 4) then   ! Find a 4th order accurate second derivative on a regular grid.
      second_deriv = ( -30.0*R(0) + (16.0*(R(1)+R(-1)) - (R(2)+R(-2))) ) / (12.0 * dx**2)
    else  ! Find a 2nd order accurate second derivative on a regular grid.
      second_deriv = ( -2.0*R(0) + (R(1)+R(-1)) ) / dx**2
    endif
  end function second_deriv

  !> Return a finite difference estimate of the second derivative with respect to two different
  !! parameters of a field in arbitrary units [A B-1 C-1]
  real function derivs_2d(R, dxdy, order)
    real,    intent(in) :: R(-3:3,-3:3) !< The field whose derivative is being taken in arbitrary units [A]
    real,    intent(in) :: dxdy   !< The spacing in two directions in parameter space in different arbitrary units [B C]
    integer, intent(in) :: order  !< The order of accuracy of the centered finite difference estimates (2, 4 or 6)

    real :: dRdx(-3:3)  ! The first derivative in one direction times the grid spacing in that direction [A]
    integer :: i

    do i=-3,3
      dRdx(i) = first_deriv(R(:,i), 1.0, order)
    enddo
    derivs_2d = first_deriv(dRdx, dxdy, order)

  end function derivs_2d

  !> Check for the rate of convergence expected with a finite difference discretization
  !! with a 20% margin of error and a tolerance for contributions from roundoff.
  logical function check_FD(val, val_fd, tol, verbose, field_name, order)
    real, intent(in) :: val       !< The derivative being checked, in arbitrary units [arbitrary]
    real, intent(in) :: val_fd(2) !< Two finite difference estimates of val taken with a spacing
                                  !! in parameter space and twice this spacing, in the same
                                  !! arbitrary units as val [arbitrary]
    real, intent(in) :: tol       !< An estimated fractional tolerance due to roundoff [arbitrary]
    logical, intent(in) :: verbose !< If true, write results to stdout
    character(len=*), intent(in) :: field_name !< A name used to describe the field in error messages
    integer, intent(in) :: order  !< The order of accuracy of the centered finite difference estimates (2, 4 or 6)

    character(len=200) :: mesg

    check_FD = ( abs(val_fd(1) - val) < (1.2*abs(val_fd(2) - val)/2**order + abs(tol)) )

    write(mesg, '(ES24.16," and ",ES24.16," differ by ",ES16.8," (",ES10.2,"), tol=",ES16.8)') &
          val, val_fd(1), val - val_fd(1), &
          2.0*(val - val_fd(1)) / (abs(val) + abs(val_fd(1)) + tiny(val)), &
          (1.2*abs(val_fd(2) - val)/2**order + abs(tol))
    ! This message is useful for debugging the two estimates:
    ! write(mesg, '(ES16.8," and ",ES16.8," or ",ES16.8," differ by ",2ES16.8," (",2ES10.2"), tol=",ES16.8)') &
    !       val, val_fd(1), val_fd(2), val - val_fd(1), val - val_fd(2), &
    !       2.0*(val - val_fd(1)) / (abs(val) + abs(val_fd(1)) + tiny(val)), &
    !       2.0*(val - val_fd(2)) / (abs(val) + abs(val_fd(2)) + tiny(val)), &
    !       (1.2*abs(val_fd(2) - val)/2**order + abs(tol))
    if (verbose .and. .not.check_FD) then
      write(stderr,*) "The values of "//trim(field_name)//" disagree. "//trim(mesg)
    elseif (verbose) then
      write(stdout,*) "The values of "//trim(field_name)//" agree: "//trim(mesg)
    endif
  end function check_FD

end procedure test_EOS_consistency
end submodule MOM_EOS_s
