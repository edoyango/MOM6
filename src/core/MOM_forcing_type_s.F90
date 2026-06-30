submodule (MOM_forcing_type) MOM_forcing_type_s
#include <MOM_memory.h>
  implicit none
contains
module procedure extractFluxes1d
  real :: htot(SZI_(G))       ! total ocean depth [H ~> m or kg m-2]
  real :: Pen_sw_tot(SZI_(G)) ! sum across all bands of Pen_SW [C H ~> degC m or degC kg m-2].
  real :: pen_sw_tot_rate(SZI_(G)) ! Summed rate of shortwave heating across bands
  real :: Ih_limit            ! inverse depth at which surface fluxes start to be limited
  real :: scale               ! scale scales away fluxes if depth < FluxRescaleDepth [nondim]
  real :: I_Cp                ! 1.0 / C_p [C Q-1 ~> kg degC J-1]
  real :: I_Cp_Hconvert       ! Unit conversion factors divided by the heat capacity
  logical :: calculate_diags  ! Indicate to calculate/update diagnostic arrays
  logical :: do_enthalpy      ! If true (default) enthalpy terms are computed in MOM6
  character(len=200) :: mesg
  integer            :: is, ie, nz, i, k, n
  logical :: do_NHR, do_NSR, do_NMIOR, do_PSWBR
  do_NHR = .false.
  do_NSR = .false.
  do_NMIOR = .false.
  do_PSWBR = .false.
  if (present(net_heat_rate)) do_NHR = .true.
  if (present(net_salt_rate)) do_NSR = .true.
  if (present(netmassinout_rate)) do_NMIOR = .true.
  if (present(pen_sw_bnd_rate)) do_PSWBR = .true.
  !}BGR

  ! GMM: by default heat content from mass entering and leaving the ocean (enthalpy)
  ! is diagnosed in this subroutine. When heat_content_evap is associated,
  ! the enthalpy terms are provided via coupler and, therefore, they do not need
  ! to be computed again.
  do_enthalpy = .true.
  if (associated(fluxes%heat_content_evap)) do_enthalpy = .false.

  Ih_limit  = 0.0 ; if (FluxRescaleDepth > 0.0) Ih_limit  = 1.0 / FluxRescaleDepth
  I_Cp      = 1.0 / tv%C_p
  I_Cp_Hconvert = 1.0 / (GV%H_to_RZ * tv%C_p)

  is = G%isc ; ie = G%iec ; nz = GV%ke

  calculate_diags = .true.

  ! error checking

  if (nsw > 0) then ; if (nsw /= optics_nbands(optics)) call MOM_error(WARNING, &
    "mismatch in the number of bands of shortwave radiation in MOM_forcing_type extract_fluxes.")
  endif

  if (.not.associated(fluxes%sw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sw is not associated.")

  if (.not.associated(fluxes%lw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%lw is not associated.")

  if (.not.associated(fluxes%latent)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%latent is not associated.")

  if (.not.associated(fluxes%sens)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sens is not associated.")

  if (.not.associated(fluxes%evap)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No evaporation defined.")

  if (.not.associated(fluxes%vprec)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%vprec not defined.")

  if ((.not.associated(fluxes%lprec)) .or. &
      (.not.associated(fluxes%fprec))) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No precipitation defined.")

  do i=is,ie ; htot(i) = h(i,1) ; enddo
  do k=2,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,k) ; enddo ; enddo

  if (nsw >= 1) then
    call extract_optics_slice(optics, j, G, GV, penSW_top=Pen_SW_bnd)
    if (do_PSWBR) call extract_optics_slice(optics, j, G, GV, penSW_top=Pen_SW_bnd_rate)
  endif

  do i=is,ie

    scale = 1.0 ; if ((Ih_limit > 0.0) .and. (htot(i)*Ih_limit < 1.0)) scale = htot(i)*Ih_limit

    ! Convert the penetrating shortwave forcing to (C * H) and reduce fluxes for shallow depths.
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)
    Pen_sw_tot(i) = 0.0
    if (nsw >= 1) then
      do n=1,nsw
        Pen_SW_bnd(n,i) = I_Cp_Hconvert*scale*dt * max(0.0, Pen_SW_bnd(n,i))
        Pen_sw_tot(i)   = Pen_sw_tot(i) + Pen_SW_bnd(n,i)
      enddo
    else
      Pen_SW_bnd(1,i) = 0.0
    endif

    if (do_PSWBR) then  ! Repeat the above code w/ dt=1s for legacy reasons
      pen_sw_tot_rate(i) = 0.0
      if (nsw >= 1) then
        do n=1,nsw
          Pen_SW_bnd_rate(n,i) = I_Cp_Hconvert*scale * max(0.0, Pen_SW_bnd_rate(n,i))
          pen_sw_tot_rate(i) = pen_sw_tot_rate(i) + pen_sw_bnd_rate(n,i)
        enddo
      else
        pen_sw_bnd_rate(1,i) = 0.0
      endif
    endif

    ! net volume/mass of liquid and solid passing through surface boundary fluxes
    netMassInOut(i) = dt * (scale * &
                                 (((((((( fluxes%lprec(i,j)        &
                                        + fluxes%fprec(i,j)      )  &
                                        + fluxes%evap(i,j)       )  &
                                        + fluxes%lrunoff(i,j)    )  &
                                        + fluxes%lrunoff_glc(i,j))  &
                                        + fluxes%vprec(i,j)      )  &
                                        + fluxes%seaice_melt(i,j))  &
                                        + fluxes%frunoff(i,j)    )  &
                                        + fluxes%frunoff_glc(i,j)))

    if (do_NMIOr) then  ! Repeat the above code without multiplying by a timestep for legacy reasons
      netMassInOut_rate(i) = (scale * &
                                 (((((((( fluxes%lprec(i,j)      &
                                        + fluxes%fprec(i,j)      )  &
                                        + fluxes%evap(i,j)       )  &
                                        + fluxes%lrunoff(i,j)    )  &
                                        + fluxes%lrunoff_glc(i,j))  &
                                        + fluxes%vprec(i,j)      )  &
                                        + fluxes%seaice_melt(i,j))  &
                                        + fluxes%frunoff(i,j)    )  &
                                        + fluxes%frunoff_glc(i,j)))
    endif

    ! smg:
    ! for non-Bouss, we add/remove salt mass to total ocean mass. to conserve
    ! total salt mass ocean+ice, the sea ice model must lose mass when salt mass
    ! is added to the ocean, which may still need to be coded.  Not that the units
    ! of netMassInOut are still [Z R ~> kg m-2], so no conversion to H should occur yet.
    if (.not.GV%Boussinesq .and. associated(fluxes%salt_flux)) then
      netMassInOut(i) = netMassInOut(i) + dt * (scale * fluxes%salt_flux(i,j))
      if (do_NMIOr) netMassInOut_rate(i) = netMassInOut_rate(i) + &
                                               (scale * fluxes%salt_flux(i,j))
    endif

    ! net volume/mass of water leaving the ocean.
    ! check that fluxes are < 0, which means mass is indeed leaving.
    netMassOut(i) = 0.0

    ! evap > 0 means condensating water is added into ocean.
    ! evap < 0 means evaporation of water from the ocean, in
    ! which case heat_content_massout is computed in MOM_diabatic_driver.F90
    if (fluxes%evap(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%evap(i,j)
  !   if (associated(fluxes%heat_content_cond)) fluxes%heat_content_cond(i,j) = 0.0 !??? --AJA

    ! lprec < 0 means sea ice formation taking water from the ocean.
    ! smg: we should split the ice melt/formation from the lprec
    if (fluxes%lprec(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%lprec(i,j)

    ! seaice_melt < 0 means sea ice formation taking water from the ocean.
    if (fluxes%seaice_melt(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%seaice_melt(i,j)

    ! vprec < 0 means virtual evaporation arising from surface salinity restoring,
    ! in which case heat_content_vprec is computed in MOM_diabatic_driver.F90.
    if (fluxes%vprec(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%vprec(i,j)

    netMassOut(i) = dt * scale * netMassOut(i)

    ! convert to H units (Bouss=meter or non-Bouss=kg/m^2)
    netMassInOut(i) = GV%RZ_to_H * netMassInOut(i)
    if (do_NMIOr) netMassInOut_rate(i) = GV%RZ_to_H * netMassInOut_rate(i)
    netMassOut(i)   = GV%RZ_to_H * netMassOut(i)

    ! surface heat fluxes from radiation and turbulent fluxes (K * H)
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)

    ! CIME provides heat flux from snow&ice melt (seaice_melt_heat), so this is added below
    ! Note: this term accounts for the enthalpy associated with water flux due to sea ice melting/freezing
    if (associated(fluxes%seaice_melt_heat)) then
      net_heat(i) = scale * dt * I_Cp_Hconvert * &
                    ( fluxes%sw(i,j) + (((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) + &
                      fluxes%seaice_melt_heat(i,j)) )
      !Repeats above code w/ dt=1. for legacy reason
      if (do_NHR)  net_heat_rate(i) = scale * I_Cp_Hconvert * &
           ( fluxes%sw(i,j) + (((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) + &
             fluxes%seaice_melt_heat(i,j)))
    else
      net_heat(i) = scale * dt * I_Cp_Hconvert * &
                    ( fluxes%sw(i,j) + ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )
      !Repeats above code w/ dt=1. for legacy reason
      if (do_NHR)  net_heat_rate(i) = scale * I_Cp_Hconvert * &
           ( fluxes%sw(i,j) + ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )
    endif

    ! Add heat flux from surface damping (restoring) (K * H) or flux adjustments.
    if (associated(fluxes%heat_added)) then
      net_heat(i) = net_heat(i) + (scale * (dt * I_Cp_Hconvert)) * fluxes%heat_added(i,j)
      if (do_NHR) net_heat_rate(i) = net_heat_rate(i) + (scale * I_Cp_Hconvert) * fluxes%heat_added(i,j)
    endif

    ! Add explicit heat flux for runoff (which is part of the ice-ocean boundary
    ! flux type). Runoff is otherwise added with a temperature of SST.
    if (useRiverHeatContent) then
      ! remove lrunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = (net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_lrunoff(i,j)) - &
                     (GV%RZ_to_H * (scale * dt)) * fluxes%lrunoff(i,j) * T(i,1)
      net_heat(i) = (net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_lrunoff_glc(i,j)) - &
                     (GV%RZ_to_H * (scale * dt)) * fluxes%lrunoff_glc(i,j) * T(i,1)
      !BGR-Jul 5, 2017{
      !Intentionally neglect the following contribution to rate for legacy reasons.
      !if (do_NHR) net_heat_rate(i) = (net_heat_rate(i) + (scale*I_Cp_Hconvert) * fluxes%heat_content_lrunoff(i,j)) - &
      !               (GV%RZ_to_H * (scale)) * fluxes%lrunoff(i,j) * T(i,1)
      !}BGR
      if (calculate_diags .and. associated(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_lrunoff(i,j) - fluxes%lrunoff(i,j)*T(i,1))
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_lrunoff_glc(i,j) - fluxes%lrunoff_glc(i,j)*T(i,1))
      endif
    endif

    ! Add explicit heat flux for calving (which is part of the ice-ocean boundary
    ! flux type). Calving is otherwise added with a temperature of SST.
    if (useCalvingHeatContent) then
      ! remove frunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_frunoff(i,j) - &
                    (GV%RZ_to_H * (scale * dt)) * fluxes%frunoff(i,j) * T(i,1)
      net_heat(i) = net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_frunoff_glc(i,j) - &
                    (GV%RZ_to_H * (scale * dt)) * fluxes%frunoff_glc(i,j) * T(i,1)
      !BGR-Jul 5, 2017{
      !Intentionally neglect the following contribution to rate for legacy reasons.
!      if (do_NHR) net_heat_rate(i) = net_heat_rate(i) + (scale*I_Cp_Hconvert) * fluxes%heat_content_frunoff(i,j) - &
!                    (GV%RZ_to_H * scale) * fluxes%frunoff(i,j) * T(i,1)
      !}BGR
      if (calculate_diags .and. associated(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_frunoff(i,j) - fluxes%frunoff(i,j)*T(i,1))
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_frunoff_glc(i,j) - fluxes%frunoff_glc(i,j)*T(i,1))
      endif
    endif

! smg: new code
    ! add heat from all terms that may add mass to the ocean (K * H).
    ! if evap, lprec, or vprec < 0, then compute their heat content
    ! inside MOM_diabatic_driver.F90 and fill in fluxes%heat_content_massout.
    ! we do so since we do not here know the temperature
    ! of water leaving the ocean, as it could be leaving from more than
    ! one layer of the upper ocean in the case of very thin layers.
    ! When evap, lprec, or vprec > 0, then we know their heat content here
    ! via settings from inside of the appropriate config_src driver files.
!    if (associated(fluxes%heat_content_lprec)) then
!      net_heat(i) = net_heat(i) + scale * dt * I_Cp_Hconvert * &
!     (fluxes%heat_content_lprec(i,j)    + (fluxes%heat_content_fprec(i,j)   + &
!     (fluxes%heat_content_lrunoff(i,j)  + (fluxes%heat_content_frunoff(i,j) + &
!     (fluxes%heat_content_cond(i,j)     +  fluxes%heat_content_vprec(i,j))))))
!    endif

    ! When enthalpy terms are provided via coupler, they must be included in net_heat
    if (.not. do_enthalpy) then
      net_heat(i) = net_heat(i) + (scale * dt * I_Cp_Hconvert * &
                    ((((fluxes%heat_content_lrunoff(i,j) + fluxes%heat_content_frunoff(i,j)) + &
                       (fluxes%heat_content_lrunoff_glc(i,j) + fluxes%heat_content_frunoff_glc(i,j))) + &
                       (fluxes%heat_content_lprec(i,j)   + fluxes%heat_content_fprec(i,j)))   + &
                       (fluxes%heat_content_evap(i,j)    + fluxes%heat_content_cond(i,j))))
    endif

    if (fluxes%num_msg < fluxes%max_msg) then
      if (Pen_SW_tot(i) > 1.000001 * I_Cp_Hconvert*scale*dt*fluxes%sw(i,j)) then
        fluxes%num_msg = fluxes%num_msg + 1
        write(mesg,'("Penetrating shortwave of ",1pe17.10, &
                    &" exceeds total shortwave of ",1pe17.10,&
                    &" at ",1pg11.4,",E,",1pg11.4,"N.")') &
               US%C_to_degC*Pen_SW_tot(i), US%C_to_degC*I_Cp_Hconvert*scale*dt * fluxes%sw(i,j), &
               G%geoLonT(i,j), G%geoLatT(i,j)
        call MOM_error(WARNING,mesg)
      endif
    endif

    ! remove penetrative portion of the SW that is NOT absorbed within a
    ! tiny layer at the top of the ocean.
    net_heat(i) = net_heat(i) - Pen_SW_tot(i)
    !Repeat above code for 'rate' term
    if (do_NHR) net_heat_rate(i) = net_heat_rate(i) - Pen_SW_tot_rate(i)

    ! diagnose non-downwelling SW
    if (present(nonPenSW)) then
      nonPenSW(i) = scale * dt * I_Cp_Hconvert * fluxes%sw(i,j) - Pen_SW_tot(i)
    endif

    ! Salt fluxes
    net_salt(i) = 0.0
    if (do_NSR) net_salt_rate(i) = 0.0
    ! Convert salt_flux from kg (salt)/(m^2 * s) to
    ! Boussinesq: (ppt * m)
    ! non-Bouss:  (g/m^2)
    if (associated(fluxes%salt_flux)) then
      net_salt(i) = (scale * dt * (1000.0*US%ppt_to_S * fluxes%salt_flux(i,j))) * GV%RZ_to_H
      !Repeat above code for 'rate' term
      if (do_NSR) net_salt_rate(i) = (scale * 1. * (1000.0*US%ppt_to_S * fluxes%salt_flux(i,j))) * GV%RZ_to_H
    endif

    ! Diagnostics follow...
    if (calculate_diags .and. do_enthalpy) then

      ! Initialize heat_content_massin that is diagnosed in mixedlayer_convection or
      ! applyBoundaryFluxes such that the meaning is as the sum of all incoming components.
      if (associated(fluxes%heat_content_massin))  then
        if (aggregate_FW) then
          if (netMassInOut(i) > 0.0) then ! net is "in"
            fluxes%heat_content_massin(i,j) = -tv%C_p * netMassOut(i) * T(i,1) * GV%H_to_RZ / dt
          else ! net is "out"
            fluxes%heat_content_massin(i,j) = tv%C_p * ( netMassInout(i) - netMassOut(i) ) * &
                                               T(i,1) * GV%H_to_RZ / dt
          endif
        else
          fluxes%heat_content_massin(i,j) = 0.
        endif
      endif

      ! Initialize heat_content_massout that is diagnosed in mixedlayer_convection or
      ! applyBoundaryFluxes such that the meaning is as the sum of all outgoing components.
      if (associated(fluxes%heat_content_massout)) then
        if (aggregate_FW) then
          if (netMassInOut(i) > 0.0) then ! net is "in"
            fluxes%heat_content_massout(i,j) = tv%C_p * netMassOut(i) * T(i,1) * GV%H_to_RZ / dt
          else ! net is "out"
            fluxes%heat_content_massout(i,j) = -tv%C_p * ( netMassInout(i) - netMassOut(i) ) * &
                                                T(i,1) * GV%H_to_RZ / dt
          endif
        else
          fluxes%heat_content_massout(i,j) = 0.0
        endif
      endif

      ! smg: we should remove sea ice melt from lprec!!!
      ! fluxes%lprec > 0 means ocean gains mass via liquid precipitation and/or sea ice melt.
      ! When atmosphere does not provide heat of this precipitation, the ocean assumes
      ! it enters the ocean at the SST.
      ! fluxes%lprec < 0 means ocean loses mass via sea ice formation. As we do not yet know
      ! the layer at which this mass is removed, we cannot compute it heat content. We must
      ! wait until MOM_diabatic_driver.F90.
      if (associated(fluxes%heat_content_lprec)) then
        if (fluxes%lprec(i,j) > 0.0) then
          fluxes%heat_content_lprec(i,j) = tv%C_p*fluxes%lprec(i,j)*T(i,1)
        else
          fluxes%heat_content_lprec(i,j) = 0.0
        endif
      endif

      ! fprec SHOULD enter ocean at 0degC if atmos model does not provide fprec heat content.
      ! However, we need to adjust netHeat above to reflect the difference between 0decC and SST
      ! and until we do so fprec is treated like lprec and enters at SST. -AJA
      if (associated(fluxes%heat_content_fprec)) then
        if (fluxes%fprec(i,j) > 0.0) then
          fluxes%heat_content_fprec(i,j) = tv%C_p*fluxes%fprec(i,j)*T(i,1)
        else
          fluxes%heat_content_fprec(i,j) = 0.0
        endif
      endif

      ! virtual precip associated with salinity restoring
      ! vprec > 0 means add water to ocean, assumed to be at SST
      ! vprec < 0 means remove water from ocean; set heat_content_vprec in MOM_diabatic_driver.F90
      if (associated(fluxes%heat_content_vprec)) then
        if (fluxes%vprec(i,j) > 0.0) then
          fluxes%heat_content_vprec(i,j) = tv%C_p*fluxes%vprec(i,j)*T(i,1)
        else
          fluxes%heat_content_vprec(i,j) = 0.0
        endif
      endif

      ! fluxes%evap < 0 means ocean loses mass due to evaporation.
      ! Evaporation leaves ocean surface at a temperature that has yet to be determined,
      ! since we do not know the precise layer that the water evaporates.  We therefore
      ! compute fluxes%heat_content_massout at the relevant point inside MOM_diabatic_driver.F90.
      ! fluxes%evap > 0 means ocean gains moisture via condensation.
      ! Condensation is assumed to drop into the ocean at the SST, just like lprec.
      if (associated(fluxes%heat_content_cond)) then
        if (fluxes%evap(i,j) > 0.0) then
          fluxes%heat_content_cond(i,j) = tv%C_p*fluxes%evap(i,j)*T(i,1)
        else
          fluxes%heat_content_cond(i,j) = 0.0
        endif
      endif

      ! Liquid runoff enters ocean at SST if land model does not provide runoff heat content.
      if (.not. useRiverHeatContent) then
        if (associated(fluxes%lrunoff) .and. associated(fluxes%heat_content_lrunoff)) then
          fluxes%heat_content_lrunoff(i,j) = tv%C_p*fluxes%lrunoff(i,j)*T(i,1)
        endif
        if (associated(fluxes%lrunoff_glc) .and. associated(fluxes%heat_content_lrunoff_glc)) then
          fluxes%heat_content_lrunoff_glc(i,j) = tv%C_p*fluxes%lrunoff_glc(i,j)*T(i,1)
        endif
      endif

      ! Icebergs enter ocean at SST if land model does not provide calving heat content.
      if (.not. useCalvingHeatContent) then
        if (associated(fluxes%frunoff) .and. associated(fluxes%heat_content_frunoff)) then
          fluxes%heat_content_frunoff(i,j) = tv%C_p*fluxes%frunoff(i,j)*T(i,1)
        endif
        if (associated(fluxes%frunoff_glc) .and. associated(fluxes%heat_content_frunoff_glc)) then
          fluxes%heat_content_frunoff_glc(i,j) = tv%C_p*fluxes%frunoff_glc(i,j)*T(i,1)
        endif
      endif

    elseif (.not. do_enthalpy) then

      ! virtual precip associated with salinity restoring. Heat content associated with
      ! that is *not* provided by the coupler and must be calculated by MOM6.
      ! vprec > 0 means add water to ocean, assumed to be at SST
      ! vprec < 0 means remove water from ocean; set heat_content_vprec in MOM_diabatic_driver.F90
      if (associated(fluxes%heat_content_vprec)) then
        if (fluxes%vprec(i,j) > 0.0) then
          fluxes%heat_content_vprec(i,j) = fluxes%C_p*fluxes%vprec(i,j)*T(i,1)
        else
          fluxes%heat_content_vprec(i,j) = 0.0
        endif
      endif

      if (associated(tv%TempxPmE)) then
        tv%TempxPmE(i,j) =  (I_Cp*dt*scale) * &
         ((((fluxes%heat_content_lprec(i,j) + fluxes%heat_content_fprec(i,j)) + &
            (fluxes%heat_content_lrunoff(i,j) + fluxes%heat_content_frunoff(i,j))) + &
            (fluxes%heat_content_lrunoff_glc(i,j) + fluxes%heat_content_frunoff_glc(i,j))) + &
            (fluxes%heat_content_evap(i,j) + fluxes%heat_content_cond(i,j)))
      endif

    endif ! calculate_diags and do_enthalpy

  enddo ! i-loop

end procedure extractFluxes1d
module procedure extractFluxes2d
  integer :: j
  do j=G%jsc, G%jec
    call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt, &
            FluxRescaleDepth, useRiverHeatContent, useCalvingHeatContent,&
            h(:,j,:), T(:,j,:), netMassInOut(:,j), netMassOut(:,j),              &
            net_heat(:,j), net_salt(:,j), pen_SW_bnd(:,:,j), tv, aggregate_FW)
  enddo

end procedure extractFluxes2d
module procedure calculateBuoyancyFlux1d
  real, dimension(SZI_(G))              :: netH       ! net FW flux [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G))              :: netEvap    ! net FW flux leaving ocean via evaporation
  real, dimension(SZI_(G))              :: netHeat    ! net temp flux [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(SZI_(G), SZK_(GV))    :: dz         ! Layer thicknesses in depth units [Z ~> m]
  real, dimension(max(nsw,1), SZI_(G))  :: penSWbnd   ! penetrating SW radiation by band
  real, dimension(SZI_(G))              :: pressure   ! pressure at the surface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G))              :: dRhodT     ! density partial derivative wrt temp [R C-1 ~> kg m-3 degC-1]
  real, dimension(SZI_(G))              :: dRhodS     ! density partial derivative wrt saln [R S-1 ~> kg m-3 ppt-1]
  real, dimension(SZI_(G))              :: dSpV_dT    ! Partial derivative of specific volume with respect
  real, dimension(SZI_(G))              :: dSpV_dS    ! Partial derivative of specific volume with respect
  real, dimension(SZI_(G),SZK_(GV)+1)   :: netPen     ! The net penetrating shortwave radiation at each level
  logical :: useRiverHeatContent
  logical :: useCalvingHeatContent
  real    :: GoRho  ! The gravitational acceleration divided by mean density times a
  real    :: g_conv ! The gravitational acceleration times the conversion factors from non-Boussinesq
  real    :: H_limit_fluxes ! A depth scale that specifies when the ocean is shallow that
  integer :: i, k
  useRiverHeatContent   = .False.
  useCalvingHeatContent = .False.

  H_limit_fluxes = max( GV%Angstrom_H, 1.e-30*GV%m_to_H )

  ! The surface forcing is contained in the fluxes type.
  ! We aggregate the thermodynamic forcing for a time step into the following:
  ! netH       = water added/removed via surface fluxes [H T-1 ~> m s-1 or kg m-2 s-1]
  ! netHeat    = heat via surface fluxes [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  ! netSalt    = salt via surface fluxes [S H T-1 ~> ppt m s-1 or gSalt m-2 s-1]
  ! Note that unlike other calls to extractFLuxes1d() that return the time-integrated flux
  ! this call returns the rate because dt=1 (in arbitrary time units)
  call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, 1.0,                        &
                H_limit_fluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netEvap, netHeatMinusSW,                 &
                netSalt, penSWbnd, tv, .false.)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call thickness_to_dz(h, tv, dz, j, G, GV)
  call sumSWoverBands(G, GV, US, h(:,j,:), dz, optics_nbands(optics), optics, j, 1.0, &
                      H_limit_fluxes, .true., penSWbnd, netPen)

  ! Adjust netSalt to reflect dilution effect of FW flux
  ! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
  netSalt(G%isc:G%iec) = netSalt(G%isc:G%iec) - Salt(G%isc:G%iec,j,1) * netH(G%isc:G%iec)

  ! Add in the SW heating for purposes of calculating the net
  ! surface buoyancy flux affecting the top layer.
  ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd, dim=1 )
  netHeat(G%isc:G%iec) = netHeatMinusSW(G%isc:G%iec) + netPen(G%isc:G%iec,1)

  ! Determine the buoyancy flux
  pressure(:) = 0.
  if (associated(tv%p_surf)) then ; do i=G%isc,G%iec ; pressure(i) = tv%p_surf(i,j) ; enddo ; endif

  if ((.not.GV%Boussinesq) .and. (.not.GV%semi_Boussinesq)) then
    g_conv = GV%g_Earth * GV%H_to_RZ

    ! Specific volume derivatives
    call calculate_specific_vol_derivs(Temp(:,j,1), Salt(:,j,1), pressure, dSpV_dT, dSpV_dS, &
                                  tv%eqn_of_state, EOS_domain(G%HI))

    ! Convert to a buoyancy flux [L2 T-3 ~> m2 s-3], first excluding penetrating SW heating
    do i=G%isc,G%iec
      buoyancyFlux(i,1) = g_conv * (dSpV_dS(i) * netSalt(i) + dSpV_dT(i) * netHeat(i))
    enddo
    ! We also have a penetrative buoyancy flux associated with penetrative SW
    do k=2,GV%ke+1 ; do i=G%isc,G%iec
      buoyancyFlux(i,k) = g_conv * ( dSpV_dT(i) * netPen(i,k) ) ! [L2 T-3 ~> m2 s-3]
    enddo ; enddo
  else
    GoRho = (GV%g_Earth * GV%H_to_Z) / GV%Rho0

    ! Density derivatives
    call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, dRhodT, dRhodS, &
                                  tv%eqn_of_state, EOS_domain(G%HI))

    ! Convert to a buoyancy flux [L2 T-3 ~> m2 s-3], excluding penetrating SW heating
    do i=G%isc,G%iec
      buoyancyFlux(i,1) = - GoRho * ( dRhodS(i) * netSalt(i) + dRhodT(i) * netHeat(i) )
    enddo
    ! We also have a penetrative buoyancy flux associated with penetrative SW
    do k=2,GV%ke+1 ; do i=G%isc,G%iec
      buoyancyFlux(i,k) = - GoRho * ( dRhodT(i) * netPen(i,k) ) ! [L2 T-3 ~> m2 s-3]
    enddo ; enddo
  endif

end procedure calculateBuoyancyFlux1d
module procedure calculateBuoyancyFlux2d
  integer :: j
  do j=G%jsc,G%jec
    call calculateBuoyancyFlux1d(G, GV, US, fluxes, optics, optics_nbands(optics), h, Temp, Salt, &
                                 tv, j, buoyancyFlux(:,j,:), netHeatMinusSW(:,j),  netSalt(:,j))
  enddo

end procedure calculateBuoyancyFlux2d
module procedure find_ustar_fluxes
  real :: I_rho        ! The inverse of the reference density [R-1 ~> m3 kg-1]
  logical :: Z_T_units ! If true, U_star is returned in units of [Z T-1 ~> m s-1], otherwise it is
  integer :: i, j, is, ie, js, je, hs
  hs = 0 ; if (present(halo)) hs = max(halo, 0)
  is = G%isc - hs ; ie = G%iec + hs ; js = G%jsc - hs ; je = G%jec + hs

  Z_T_units = .true. ; if (present(H_T_units)) Z_T_units = .not.H_T_units

  if (.not.(associated(fluxes%ustar) .or. associated(fluxes%tau_mag))) &
    call MOM_error(FATAL, "find_ustar_fluxes requires that either ustar or tau_mag be associated.")

  if (associated(fluxes%ustar) .and. (GV%Boussinesq .or. .not.associated(fluxes%tau_mag))) then
    if (Z_T_units) then
      do j=js,je ; do i=is,ie
        U_star(i,j) = fluxes%ustar(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        U_star(i,j) = GV%Z_to_H * fluxes%ustar(i,j)
      enddo ; enddo
    endif
  elseif (allocated(tv%SpV_avg)) then
    if (tv%valid_SpV_halo < 0) call MOM_error(FATAL, &
        "find_ustar_fluxes called in non-Boussinesq mode with invalid values of SpV_avg.")
    if (tv%valid_SpV_halo < hs) call MOM_error(FATAL, &
        "find_ustar_fluxes called in non-Boussinesq mode with insufficient valid values of SpV_avg.")
    if (Z_T_units) then
      do j=js,je ; do i=is,ie
        U_star(i,j) = sqrt(fluxes%tau_mag(i,j) * tv%SpV_avg(i,j,1))
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        U_star(i,j) = GV%RZ_to_H * sqrt(fluxes%tau_mag(i,j) / tv%SpV_avg(i,j,1))
      enddo ; enddo
    endif
  else
    I_rho = GV%Z_to_H * GV%RZ_to_H
    if (Z_T_units) I_rho = GV%H_to_Z * GV%RZ_to_H ! == 1.0 / GV%Rho0
    do j=js,je ; do i=is,ie
      U_star(i,j) = sqrt(fluxes%tau_mag(i,j) * I_rho)
    enddo ; enddo
  endif

end procedure find_ustar_fluxes
module procedure find_ustar_mech_forcing
  real :: I_rho        ! The inverse of the reference density [R-1 ~> m3 kg-1] or in some semi-Boussinesq cases
  logical :: Z_T_units ! If true, U_star is returned in units of [Z T-1 ~> m s-1], otherwise it is
  integer :: i, j, is, ie, js, je, hs
  hs = 0 ; if (present(halo)) hs = max(halo, 0)
  is = G%isc - hs ; ie = G%iec + hs ; js = G%jsc - hs ; je = G%jec + hs

  Z_T_units = .true. ; if (present(H_T_units)) Z_T_units = .not.H_T_units

  if (.not.(associated(forces%ustar) .or. associated(forces%tau_mag))) &
    call MOM_error(FATAL, "find_ustar_mech requires that either ustar or tau_mag be associated.")

  if (associated(forces%ustar) .and. (GV%Boussinesq .or. .not.associated(forces%tau_mag))) then
    if (Z_T_units) then
      do j=js,je ; do i=is,ie
        U_star(i,j) = forces%ustar(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        U_star(i,j) = GV%Z_to_H * forces%ustar(i,j)
      enddo ; enddo
    endif
  elseif (allocated(tv%SpV_avg)) then
    if (tv%valid_SpV_halo < 0) call MOM_error(FATAL, &
        "find_ustar_mech called in non-Boussinesq mode with invalid values of SpV_avg.")
    if (tv%valid_SpV_halo < hs) call MOM_error(FATAL, &
        "find_ustar_mech called in non-Boussinesq mode with insufficient valid values of SpV_avg.")
    if (Z_T_units) then
      do j=js,je ; do i=is,ie
        U_star(i,j) = sqrt(forces%tau_mag(i,j) * tv%SpV_avg(i,j,1))
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        U_star(i,j) = GV%RZ_to_H * sqrt(forces%tau_mag(i,j) / tv%SpV_avg(i,j,1))
      enddo ; enddo
    endif
  else
    I_rho = GV%Z_to_H * GV%RZ_to_H
    if (Z_T_units) I_rho = GV%H_to_Z * GV%RZ_to_H ! == 1.0 / GV%Rho0
    do j=js,je ; do i=is,ie
      U_star(i,j) = sqrt(forces%tau_mag(i,j) * I_rho)
    enddo ; enddo
  endif

end procedure find_ustar_mech_forcing
module procedure MOM_forcing_chksum
  integer :: hshift
  hshift = 1 ; if (present(haloshift)) hshift = haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(fluxes%ustar)) &
    call hchksum(fluxes%ustar, mesg//" fluxes%ustar", G%HI, haloshift=hshift, unscale=US%Z_to_m*US%s_to_T)
  if (associated(fluxes%tau_mag)) &
    call hchksum(fluxes%tau_mag, mesg//" fluxes%tau_mag", G%HI, haloshift=hshift, unscale=US%RLZ_T2_to_Pa*US%Z_to_L)
  if (associated(fluxes%buoy)) &
    call hchksum(fluxes%buoy, mesg//" fluxes%buoy ", G%HI, haloshift=hshift, unscale=US%L_to_m**2*US%s_to_T**3)
  if (associated(fluxes%sw)) &
    call hchksum(fluxes%sw, mesg//" fluxes%sw", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_vis_dir)) &
    call hchksum(fluxes%sw_vis_dir, mesg//" fluxes%sw_vis_dir", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_vis_dif)) &
    call hchksum(fluxes%sw_vis_dif, mesg//" fluxes%sw_vis_dif", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_nir_dir)) &
    call hchksum(fluxes%sw_nir_dir, mesg//" fluxes%sw_nir_dir", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_nir_dif)) &
    call hchksum(fluxes%sw_nir_dif, mesg//" fluxes%sw_nir_dif", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%lw)) &
    call hchksum(fluxes%lw, mesg//" fluxes%lw", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent)) &
    call hchksum(fluxes%latent, mesg//" fluxes%latent", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_evap_diag)) &
    call hchksum(fluxes%latent_evap_diag, mesg//" fluxes%latent_evap_diag", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_fprec_diag)) &
    call hchksum(fluxes%latent_fprec_diag, mesg//" fluxes%latent_fprec_diag", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_frunoff_diag)) &
    call hchksum(fluxes%latent_frunoff_diag, mesg//" fluxes%latent_frunoff_diag", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_frunoff_glc_diag)) &
    call hchksum(fluxes%latent_frunoff_glc_diag, mesg//" fluxes%latent_frunoff_glc_diag", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens", G%HI, haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%lprec)) &
    call hchksum(fluxes%lprec, mesg//" fluxes%lprec", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%fprec)) &
    call hchksum(fluxes%fprec, mesg//" fluxes%fprec", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%vprec)) &
    call hchksum(fluxes%vprec, mesg//" fluxes%vprec", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%seaice_melt)) &
    call hchksum(fluxes%seaice_melt, mesg//" fluxes%seaice_melt", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%seaice_melt_heat)) &
    call hchksum(fluxes%seaice_melt_heat, mesg//" fluxes%seaice_melt_heat", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf", G%HI, haloshift=hshift, unscale=US%RL2_T2_to_Pa)
  if (associated(fluxes%u10_sqr)) &
    call hchksum(fluxes%u10_sqr, mesg//" fluxes%u10_sqr", G%HI, haloshift=hshift, unscale=US%L_to_m**2*US%s_to_T**2)
  if (associated(fluxes%ice_fraction)) &
    call hchksum(fluxes%ice_fraction, mesg//" fluxes%ice_fraction", G%HI, haloshift=hshift)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%BBL_tidal_dis)) &
    call hchksum(fluxes%BBL_tidal_dis, mesg//" fluxes%BBL_tidal_dis", G%HI, haloshift=hshift, &
                 unscale=US%L_to_Z**2*US%RZ3_T3_to_W_m2)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal", G%HI, haloshift=hshift, unscale=US%Z_to_m*US%s_to_T)
  if (associated(fluxes%lrunoff)) &
    call hchksum(fluxes%lrunoff, mesg//" fluxes%lrunoff", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%lrunoff_glc)) &
    call hchksum(fluxes%lrunoff_glc, mesg//" fluxes%lrunoff_glc", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%frunoff)) &
    call hchksum(fluxes%frunoff, mesg//" fluxes%frunoff", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%frunoff_glc)) &
    call hchksum(fluxes%frunoff_glc, mesg//" fluxes%frunoff_glc", G%HI, haloshift=hshift, unscale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%heat_content_lrunoff)) &
    call hchksum(fluxes%heat_content_lrunoff, mesg//" fluxes%heat_content_lrunoff", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_lrunoff_glc)) &
    call hchksum(fluxes%heat_content_lrunoff_glc, mesg//" fluxes%heat_content_lrunoff_glc", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_frunoff)) &
    call hchksum(fluxes%heat_content_frunoff, mesg//" fluxes%heat_content_frunoff", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_frunoff_glc)) &
    call hchksum(fluxes%heat_content_frunoff_glc, mesg//" fluxes%heat_content_frunoff_glc", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_lprec)) &
    call hchksum(fluxes%heat_content_lprec, mesg//" fluxes%heat_content_lprec", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_fprec)) &
    call hchksum(fluxes%heat_content_fprec, mesg//" fluxes%heat_content_fprec", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_cond)) &
    call hchksum(fluxes%heat_content_cond, mesg//" fluxes%heat_content_cond", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_evap)) &
    call hchksum(fluxes%heat_content_evap, mesg//" fluxes%heat_content_evap", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_massout)) &
    call hchksum(fluxes%heat_content_massout, mesg//" fluxes%heat_content_massout", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_massin)) &
    call hchksum(fluxes%heat_content_massin, mesg//" fluxes%heat_content_massin", G%HI, &
                 haloshift=hshift, unscale=US%QRZ_T_to_W_m2)
end procedure MOM_forcing_chksum
module procedure MOM_mech_forcing_chksum
  integer :: hshift
  hshift = 1 ; if (present(haloshift)) hshift = haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(forces%taux) .and. associated(forces%tauy)) &
    call uvchksum(mesg//" forces%tau[xy]", forces%taux, forces%tauy, G%HI, &
                  haloshift=hshift, symmetric=.true., unscale=US%RLZ_T2_to_Pa)
  if (associated(forces%p_surf)) &
    call hchksum(forces%p_surf, mesg//" forces%p_surf", G%HI, haloshift=hshift, unscale=US%RL2_T2_to_Pa)
  if (associated(forces%ustar)) &
    call hchksum(forces%ustar, mesg//" forces%ustar", G%HI, haloshift=hshift, unscale=US%Z_to_m*US%s_to_T)
  if (associated(forces%tau_mag)) &
    call hchksum(forces%tau_mag, mesg//" forces%tau_mag", G%HI, haloshift=hshift, unscale=US%RLZ_T2_to_Pa*US%Z_to_L)
  if (associated(forces%rigidity_ice_u) .and. associated(forces%rigidity_ice_v)) &
    call uvchksum(mesg//" forces%rigidity_ice_[uv]", forces%rigidity_ice_u, &
        forces%rigidity_ice_v, G%HI, haloshift=hshift, symmetric=.true., &
        unscale=US%L_to_m**3*US%L_to_Z*US%s_to_T, scalar_pair=.true.)

end procedure MOM_mech_forcing_chksum
module procedure mech_forcing_SinglePointPrint
  write(0,'(2a)') 'MOM_forcing_type, forcing_SinglePointPrint: Called from ',mesg
  write(0,'(a,2es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: lon,lat = ',G%geoLonT(i,j),G%geoLatT(i,j)
  call locMsg(forces%taux,'taux')
  call locMsg(forces%tauy,'tauy')

  contains
  !> Format and write a message depending on associated state of array
  subroutine locMsg(array,aname)
    real, dimension(:,:), pointer :: array !< Array to write element from
    character(len=*)              :: aname !< Name of array

    if (associated(array)) then
      write(0,'(3a,es15.3)') 'MOM_forcing_type, mech_forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
    else
      write(0,'(4a)') 'MOM_forcing_type, mech_forcing_SinglePointPrint: ',trim(aname),' is not associated.'
    endif
  end subroutine locMsg

end procedure mech_forcing_SinglePointPrint
module procedure forcing_SinglePointPrint
  write(0,'(2a)') 'MOM_forcing_type, forcing_SinglePointPrint: Called from ',mesg
  write(0,'(a,2es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: lon,lat = ',G%geoLonT(i,j),G%geoLatT(i,j)
  call locMsg(fluxes%ustar,'ustar')
  call locMsg(fluxes%tau_mag,'tau_mag')
  call locMsg(fluxes%buoy,'buoy')
  call locMsg(fluxes%sw,'sw')
  call locMsg(fluxes%sw_vis_dir,'sw_vis_dir')
  call locMsg(fluxes%sw_vis_dif,'sw_vis_dif')
  call locMsg(fluxes%sw_nir_dir,'sw_nir_dir')
  call locMsg(fluxes%sw_nir_dif,'sw_nir_dif')
  call locMsg(fluxes%lw,'lw')
  call locMsg(fluxes%latent,'latent')
  call locMsg(fluxes%latent_evap_diag,'latent_evap_diag')
  call locMsg(fluxes%latent_fprec_diag,'latent_fprec_diag')
  call locMsg(fluxes%latent_frunoff_diag,'latent_frunoff_diag')
  call locMsg(fluxes%latent_frunoff_glc_diag,'latent_frunoff_glc_diag')
  call locMsg(fluxes%sens,'sens')
  call locMsg(fluxes%evap,'evap')
  call locMsg(fluxes%lprec,'lprec')
  call locMsg(fluxes%fprec,'fprec')
  call locMsg(fluxes%vprec,'vprec')
  call locMsg(fluxes%seaice_melt,'seaice_melt')
  call locMsg(fluxes%seaice_melt_heat,'seaice_melt_heat')
  call locMsg(fluxes%p_surf,'p_surf')
  call locMsg(fluxes%salt_flux,'salt_flux')
  call locMsg(fluxes%BBL_tidal_dis,'BBL_tidal_dis')
  call locMsg(fluxes%ustar_tidal,'ustar_tidal')
  call locMsg(fluxes%lrunoff,'lrunoff')
  call locMsg(fluxes%lrunoff_glc,'lrunoff_glc')
  call locMsg(fluxes%frunoff,'frunoff')
  call locMsg(fluxes%frunoff_glc,'frunoff_glc')
  call locMsg(fluxes%heat_content_lrunoff,'heat_content_lrunoff')
  call locMsg(fluxes%heat_content_lrunoff_glc,'heat_content_lrunoff_glc')
  call locMsg(fluxes%heat_content_frunoff,'heat_content_frunoff')
  call locMsg(fluxes%heat_content_frunoff_glc,'heat_content_frunoff_glc')
  call locMsg(fluxes%heat_content_lprec,'heat_content_lprec')
  call locMsg(fluxes%heat_content_fprec,'heat_content_fprec')
  call locMsg(fluxes%heat_content_vprec,'heat_content_vprec')
  call locMsg(fluxes%heat_content_cond,'heat_content_cond')
  call locMsg(fluxes%heat_content_cond,'heat_content_massout')
  call locMsg(fluxes%heat_content_evap,'heat_content_evap')
  call locMsg(fluxes%heat_content_massout,'heat_content_massout')
  call locMsg(fluxes%heat_content_massin,'heat_content_massin')

  contains
  !> Format and write a message depending on associated state of array
  subroutine locMsg(array,aname)
    real, dimension(:,:), pointer :: array !< Array to write element from
    character(len=*)              :: aname !< Name of array

    if (associated(array)) then
      write(0,'(3a,es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
    else
      write(0,'(4a)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' is not associated.'
    endif
  end subroutine locMsg

end procedure forcing_SinglePointPrint
module procedure register_forcing_type_diags
  handles%id_clock_forcing=cpu_clock_id('(Ocean forcing diagnostics)', grain=CLOCK_ROUTINE)


  handles%id_taux = register_diag_field('ocean_model', 'taux', diag%axesCu1, Time,  &
        'Zonal surface stress from ocean interactions with atmos and ice',          &
        'Pa', conversion=US%RLZ_T2_to_Pa,                                           &
        standard_name='surface_downward_x_stress', cmor_field_name='tauuo',         &
        cmor_units='N m-2', cmor_long_name='Surface Downward X Stress',             &
        cmor_standard_name='surface_downward_x_stress')

  handles%id_tauy = register_diag_field('ocean_model', 'tauy', diag%axesCv1, Time,  &
        'Meridional surface stress ocean interactions with atmos and ice',         &
        'Pa', conversion=US%RLZ_T2_to_Pa,                                          &
        standard_name='surface_downward_y_stress', cmor_field_name='tauvo',        &
        cmor_units='N m-2', cmor_long_name='Surface Downward Y Stress',            &
        cmor_standard_name='surface_downward_y_stress')

  handles%id_tau_mag = register_diag_field('ocean_model', 'tau_mag', diag%axesT1, Time, &
        'Average magnitude of the wind stress including contributions from gustiness', &
        'Pa', conversion=US%RLZ_T2_to_Pa*US%Z_to_L)

  handles%id_ustar = register_diag_field('ocean_model', 'ustar', diag%axesT1, Time, &
      'Surface friction velocity = [(gustiness + tau_magnitude)/rho0]^(1/2)', &
      'm s-1', conversion=US%Z_to_m*US%s_to_T)

  handles%id_omega_w2x = register_diag_field('ocean_model', 'omega_w2x', diag%axesT1, Time, &
      'Counter-clockwise angle of the wind stress from the horizontal axis.', 'rad', conversion=1.0)

  if (present(use_berg_fluxes)) then
    if (use_berg_fluxes) then
      handles%id_ustar_berg = register_diag_field('ocean_model', 'ustar_berg', diag%axesT1, Time, &
          'Friction velocity below iceberg ', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

      handles%id_area_berg = register_diag_field('ocean_model', 'area_berg', diag%axesT1, Time, &
          'Area of grid cell covered by iceberg ', 'm2 m-2', conversion=1.0)

      handles%id_mass_berg = register_diag_field('ocean_model', 'mass_berg', diag%axesT1, Time, &
          'Mass of icebergs ', 'kg m-2', conversion=US%RZ_to_kg_m2)

      handles%id_ustar_ice_cover = register_diag_field('ocean_model', 'ustar_ice_cover', diag%axesT1, Time, &
          'Friction velocity below iceberg and ice shelf together', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

      handles%id_frac_ice_cover = register_diag_field('ocean_model', 'frac_ice_cover', diag%axesT1, Time, &
          'Area of grid cell below iceberg and ice shelf together ', 'm2 m-2', conversion=1.0)
    endif
  endif

  ! See:
  if (present(use_cfcs)) then
    if (use_cfcs) then
      handles%id_ice_fraction = register_diag_field('ocean_model', 'ice_fraction', diag%axesT1, Time, &
          'Fraction of cell area covered by sea ice', 'm2 m-2', conversion=1.0)

      handles%id_u10_sqr = register_diag_field('ocean_model', 'u10_sqr', diag%axesT1, Time, &
          'Wind magnitude at 10m, squared', 'm2 s-2', conversion=US%L_to_m**2*US%s_to_T**2)
    endif
  endif

  handles%id_psurf = register_diag_field('ocean_model', 'p_surf', diag%axesT1, Time, &
        'Pressure at ice-ocean or atmosphere-ocean interface', &
        'Pa', conversion=US%RL2_T2_to_Pa, cmor_field_name='pso', &
        cmor_long_name='Sea Water Pressure at Sea Water Surface', &
        cmor_standard_name='sea_water_pressure_at_sea_water_surface')

  handles%id_TKE_tidal = register_diag_field('ocean_model', 'TKE_tidal', diag%axesT1, Time, &
        'Tidal source of BBL mixing', 'W m-2', conversion=US%L_to_Z**2*US%RZ3_T3_to_W_m2)

  if (.not. use_temperature) then
    handles%id_buoy = register_diag_field('ocean_model', 'buoy', diag%axesT1, Time, &
          'Buoyancy forcing', 'm2 s-3', conversion=US%L_to_m**2*US%s_to_T**3)
    return
  endif


  !===============================================================
  ! surface mass flux maps

  handles%id_prcme = register_diag_field('ocean_model', 'PRCmE', diag%axesT1, Time, &
        'Net surface water flux (precip+melt+lrunoff+ice calving-evap)',  &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                       &
        standard_name='water_flux_into_sea_water', cmor_field_name='wfo', &
        cmor_standard_name='water_flux_into_sea_water',cmor_long_name='Water Flux Into Sea Water')

  handles%id_evap = register_diag_field('ocean_model', 'evap', diag%axesT1, Time, &
        'Evaporation/condensation at ocean surface (evaporation is negative)', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        standard_name='water_evaporation_flux', cmor_field_name='evs', &
        cmor_standard_name='water_evaporation_flux', &
        cmor_long_name='Water Evaporation Flux Where Ice Free Ocean over Sea')

  ! smg: seaice_melt field requires updates to the sea ice model
  handles%id_seaice_melt = register_diag_field('ocean_model', 'seaice_melt',       &
        diag%axesT1, Time, 'water flux to ocean from snow/sea ice melting(> 0) or formation(< 0)', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                                  &
        standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',     &
        cmor_field_name='fsitherm',                                                  &
        cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',&
        cmor_long_name='water flux to ocean from sea ice melt(> 0) or form(< 0)')

  handles%id_precip = register_diag_field('ocean_model', 'precip', diag%axesT1, Time, &
        'Liquid + frozen precipitation into ocean', 'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_fprec = register_diag_field('ocean_model', 'fprec', diag%axesT1, Time,     &
        'Frozen precipitation into ocean',                                              &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                               &
        standard_name='snowfall_flux', cmor_field_name='prsn',                          &
        cmor_standard_name='snowfall_flux', cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea')

  handles%id_lprec = register_diag_field('ocean_model', 'lprec', diag%axesT1, Time,       &
        'Liquid precipitation into ocean',                                                &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                                 &
        standard_name='rainfall_flux',                                                    &
        cmor_field_name='prlq', cmor_standard_name='rainfall_flux',                       &
        cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea')

  handles%id_vprec = register_diag_field('ocean_model', 'vprec', diag%axesT1, Time, &
        'Virtual liquid precip into ocean due to SSS restoring',                    &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_frunoff = register_diag_field('ocean_model', 'frunoff', diag%axesT1, Time,    &
        'Frozen runoff (calving) and iceberg melt into ocean',                             &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                                  &
        standard_name='water_flux_into_sea_water_from_icebergs',                           &
        cmor_field_name='ficeberg',                                                        &
        cmor_standard_name='water_flux_into_sea_water_from_icebergs',                      &
        cmor_long_name='Water Flux into Seawater from Icebergs')

  handles%id_lrunoff = register_diag_field('ocean_model', 'lrunoff', diag%axesT1, Time, &
        'Liquid runoff (rivers) into ocean',                                                  &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                                     &
        standard_name='water_flux_into_sea_water_from_rivers', cmor_field_name='friver',      &
        cmor_standard_name='water_flux_into_sea_water_from_rivers',                           &
        cmor_long_name='Water Flux into Sea Water From Rivers')

  if (present(use_glc_runoff)) then
    handles%id_frunoff_glc = register_diag_field('ocean_model', 'frunoff_glc', diag%axesT1, Time,    &
          'Frozen glacier runoff (calving) and iceberg melt into ocean', &
          units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
          standard_name='glc_water_flux_into_sea_water_from_icebergs') ! todo: update cmor names

    handles%id_lrunoff_glc = register_diag_field('ocean_model', 'lrunoff_glc', diag%axesT1, Time, &
          'Liquid runoff (glaciers) into ocean', &
          units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
          standard_name='water_flux_into_sea_water_from_glaciers') ! todo: update cmor names
  endif

  handles%id_net_massout = register_diag_field('ocean_model', 'net_massout', diag%axesT1, Time, &
        'Net mass leaving the ocean due to evaporation, seaice formation', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_net_massin  = register_diag_field('ocean_model', 'net_massin', diag%axesT1, Time, &
        'Net mass entering ocean due to precip, runoff, ice melt', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_massout_flux = register_diag_field('ocean_model', 'massout_flux', diag%axesT1, Time, &
        'Net mass flux of freshwater out of the ocean (used in the boundary flux calculation)', &
         'kg m-2', conversion=diag%GV%H_to_kg_m2)

  handles%id_massin_flux  = register_diag_field('ocean_model', 'massin_flux', diag%axesT1, Time, &
        'Net mass flux of freshwater into the ocean (used in boundary flux calculation)', &
        'kg m-2', conversion=diag%GV%H_to_kg_m2)

  !=========================================================================
  ! area integrated surface mass transport, all are rescaled to MKS units before area integration.

  handles%id_total_prcme = register_scalar_field('ocean_model', 'total_PRCmE', Time, diag,         &
      long_name='Area integrated net surface water flux (precip+melt+liq runoff+ice calving-evap)', &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                          &
      standard_name='water_flux_into_sea_water_area_integrated',                                   &
      cmor_field_name='total_wfo',                                                                 &
      cmor_standard_name='water_flux_into_sea_water_area_integrated',                              &
      cmor_long_name='Water Transport Into Sea Water Area Integrated')

  handles%id_total_evap = register_scalar_field('ocean_model', 'total_evap', Time, diag,&
      long_name='Area integrated evap/condense at ocean surface',                       &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                               &
      standard_name='water_evaporation_flux_area_integrated',                           &
      cmor_field_name='total_evs',                                                      &
      cmor_standard_name='water_evaporation_flux_area_integrated',                      &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Integrated')

  ! seaice_melt field requires updates to the sea ice model
  handles%id_total_seaice_melt = register_scalar_field('ocean_model', 'total_icemelt', Time, diag, &
      long_name='Area integrated sea ice melt (>0) or form (<0)',                                      &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                              &
      standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',         &
      cmor_field_name='total_fsitherm',                                                                &
      cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',    &
      cmor_long_name='Water Melt/Form from Sea Ice Area Integrated')

  handles%id_total_precip = register_scalar_field('ocean_model', 'total_precip', Time, diag, &
      long_name='Area integrated liquid+frozen precip into ocean', &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)

  handles%id_total_fprec = register_scalar_field('ocean_model', 'total_fprec', Time, diag,&
      long_name='Area integrated frozen precip into ocean',                               &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                 &
      standard_name='snowfall_flux_area_integrated',                                      &
      cmor_field_name='total_prsn',                                                       &
      cmor_standard_name='snowfall_flux_area_integrated',                                 &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_lprec = register_scalar_field('ocean_model', 'total_lprec', Time, diag,&
      long_name='Area integrated liquid precip into ocean',                               &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                 &
      standard_name='rainfall_flux_area_integrated',                                      &
      cmor_field_name='total_pr',                                                         &
      cmor_standard_name='rainfall_flux_area_integrated',                                 &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_vprec = register_scalar_field('ocean_model', 'total_vprec', Time, diag, &
      long_name='Area integrated virtual liquid precip due to SSS restoring', &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)

  handles%id_total_frunoff = register_scalar_field('ocean_model', 'total_frunoff', Time, diag,    &
      long_name='Area integrated frozen runoff (calving) & iceberg melt into ocean',              &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                         &
      cmor_field_name='total_ficeberg',                                                           &
      cmor_standard_name='water_flux_into_sea_water_from_icebergs_area_integrated',               &
      cmor_long_name='Water Flux into Seawater from Icebergs Area Integrated')

  handles%id_total_lrunoff = register_scalar_field('ocean_model', 'total_lrunoff', Time, diag,&
      long_name='Area integrated liquid runoff into ocean',                                   &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T,                                     &
      cmor_field_name='total_friver',                                                         &
      cmor_standard_name='water_flux_into_sea_water_from_rivers_area_integrated',             &
      cmor_long_name='Water Flux into Sea Water From Rivers Area Integrated')

  if (present(use_glc_runoff)) then
    handles%id_total_frunoff_glc = register_scalar_field('ocean_model', 'total_frunoff_glc', Time, diag, &
        long_name='Area integrated frozen glacier runoff (calving) & iceberg melt into ocean', &
        units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)

    handles%id_total_lrunoff_glc = register_scalar_field('ocean_model', 'total_lrunoff_glc', Time, diag, &
        long_name='Area integrated liquid glacier runoff into ocean', &
        units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)
  endif

  handles%id_total_net_massout = register_scalar_field('ocean_model', 'total_net_massout', Time, diag, &
      long_name='Area integrated mass leaving ocean due to evap and seaice form', &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)

  handles%id_total_net_massin = register_scalar_field('ocean_model', 'total_net_massin', Time, diag, &
      long_name='Area integrated mass entering ocean due to predip, runoff, ice melt', &
      units='kg s-1', conversion=US%RZL2_to_kg*US%s_to_T)

  !=========================================================================
  ! area averaged surface mass transport

  handles%id_prcme_ga = register_scalar_field('ocean_model', 'PRCmE_ga', Time, diag, &
      long_name='Area averaged net surface water flux (precip+melt+liq runoff+ice calving-evap)', &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
      standard_name='water_flux_into_sea_water_area_averaged', &
      cmor_field_name='ave_wfo', cmor_standard_name='rainfall_flux_area_averaged', &
      cmor_long_name='Water Transport Into Sea Water Area Averaged')

  handles%id_evap_ga = register_scalar_field('ocean_model', 'evap_ga', Time, diag, &
      long_name='Area averaged evap/condense at ocean surface', &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
      standard_name='water_evaporation_flux_area_averaged', &
      cmor_field_name='ave_evs', cmor_standard_name='water_evaporation_flux_area_averaged', &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Averaged')

  handles%id_lprec_ga = register_scalar_field('ocean_model', 'lprec_ga', Time, diag,&
      long_name='Area integrated liquid precip into ocean', &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
      standard_name='rainfall_flux_area_averaged', &
      cmor_field_name='ave_pr', cmor_standard_name='rainfall_flux_area_averaged', &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Averaged')

  handles%id_fprec_ga = register_scalar_field('ocean_model', 'fprec_ga', Time, diag, &
      long_name='Area integrated frozen precip into ocean',                        &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                            &
      standard_name='snowfall_flux_area_averaged',                                 &
      cmor_field_name='ave_prsn',cmor_standard_name='snowfall_flux_area_averaged', &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Averaged')

  handles%id_precip_ga = register_scalar_field('ocean_model', 'precip_ga', Time, diag, &
      long_name='Area averaged liquid+frozen precip into ocean', &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_vprec_ga = register_scalar_field('ocean_model', 'vrec_ga', Time, diag, &
      long_name='Area averaged virtual liquid precip due to SSS restoring', &
      units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  !===============================================================
  ! surface heat flux maps

  handles%id_heat_content_frunoff = register_diag_field('ocean_model', 'heat_content_frunoff', &
        diag%axesT1, Time, 'Heat content (relative to 0C) of solid runoff into ocean',         &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lrunoff = register_diag_field('ocean_model', 'heat_content_lrunoff', &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid runoff into ocean',        &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  if (present(use_carbon_runoff)) then
    if (use_carbon_runoff) then
       handles%id_carbon_content_lrunoff = register_diag_field('ocean_model', 'carbon_content_lrunoff', &
             diag%axesT1, Time, 'Carbon content of liquid runoff into ocean',        &
             'kg m-2 s-1', standard_name='carbon_flux_due_to_runoff')
    endif
  endif

  if (present(use_glc_runoff)) then
    handles%id_heat_content_frunoff_glc = register_diag_field('ocean_model', 'heat_content_frunoff_glc', &
          diag%axesT1, Time, 'Heat content (relative to 0C) of solid glacier runoff into ocean',         &
          'W m-2', conversion=US%QRZ_T_to_W_m2)

    handles%id_heat_content_lrunoff_glc = register_diag_field('ocean_model', 'heat_content_lrunoff_glc', &
          diag%axesT1, Time, 'Heat content (relative to 0C) of liquid glacier runoff into ocean',        &
          'W m-2', conversion=US%QRZ_T_to_W_m2)
  endif

  handles%id_hfrunoffds = register_diag_field('ocean_model', 'hfrunoffds',                            &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid+solid runoff into ocean', &
        'W m-2', conversion=US%QRZ_T_to_W_m2,                                                 &
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lprec = register_diag_field('ocean_model', 'heat_content_lprec',             &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid precip entering ocean',           &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_fprec = register_diag_field('ocean_model', 'heat_content_fprec',&
        diag%axesT1,Time,'Heat content (relative to 0degC) of frozen prec entering ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_vprec = register_diag_field('ocean_model', 'heat_content_vprec',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of virtual precip entering ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_cond = register_diag_field('ocean_model', 'heat_content_cond',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of water condensing into ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_evap = register_diag_field('ocean_model', 'heat_content_evap',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of water evaporating from ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_hfrainds = register_diag_field('ocean_model', 'hfrainds',                                 &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid+frozen precip entering ocean',    &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water',&
        cmor_long_name='Heat Content (relative to 0degC) of Liquid + Frozen Precipitation')

  handles%id_heat_content_surfwater = register_diag_field('ocean_model', 'heat_content_surfwater',&
         diag%axesT1, Time,                                                                       &
        'Heat content (relative to 0degC) of net water crossing ocean surface (frozen+liquid)',   &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_massout = register_diag_field('ocean_model', 'heat_content_massout',                      &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass leaving ocean ocean via evap and ice form',&
        'W m-2', conversion=US%QRZ_T_to_W_m2,                                                      &
        cmor_field_name='hfevapds',                                                                                 &
        cmor_standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water',           &
        cmor_long_name='Heat Content (relative to 0degC) of Water Leaving Ocean via Evaporation and Ice Formation')

  handles%id_heat_content_massin = register_diag_field('ocean_model', 'heat_content_massin',   &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass entering ocean ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_net_heat_coupler = register_diag_field('ocean_model', 'net_heat_coupler',          &
        diag%axesT1,Time,'Surface ocean heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_net_heat_surface = register_diag_field('ocean_model', 'net_heat_surface',diag%axesT1, Time,  &
        'Surface ocean heat flux from SW+LW+lat+sens+mass transfer+frazil+restore+seaice_melt_heat or '// &
        'flux adjustments', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='surface_downward_heat_flux_in_sea_water', cmor_field_name='hfds',            &
        cmor_standard_name='surface_downward_heat_flux_in_sea_water',           &
        cmor_long_name='Surface ocean heat flux from SW+LW+latent+sensible+masstransfer+frazil+seaice_melt_heat')

  handles%id_sw = register_diag_field('ocean_model', 'SW', diag%axesT1, Time,  &
        'Shortwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='net_downward_shortwave_flux_at_sea_water_surface',      &
        cmor_field_name='rsntds',                                              &
        cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface', &
        cmor_long_name='Net Downward Shortwave Radiation at Sea Water Surface')
  handles%id_sw_vis = register_diag_field('ocean_model', 'sw_vis', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the visible band', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)
  handles%id_sw_nir = register_diag_field('ocean_model', 'sw_nir', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the near-infrared band', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', diag%axesT1, Time, &
        'Combined longwave, latent, and sensible heating at ocean surface', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_lw = register_diag_field('ocean_model', 'LW', diag%axesT1, Time, &
        'Longwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='surface_net_downward_longwave_flux',                   &
        cmor_field_name='rlntds',                                             &
        cmor_standard_name='surface_net_downward_longwave_flux',              &
        cmor_long_name='Surface Net Downward Longwave Radiation')

  handles%id_lat = register_diag_field('ocean_model', 'latent', diag%axesT1, Time,                    &
        'Latent heat flux into ocean due to fusion and evaporation (negative means ocean heat loss)', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, cmor_field_name='hflso',                                &
        cmor_standard_name='surface_downward_latent_heat_flux',                                       &
        cmor_long_name='Surface Downward Latent Heat Flux due to Evap + Melt Snow/Ice')

  handles%id_lat_evap = register_diag_field('ocean_model', 'latent_evap', diag%axesT1, Time, &
        'Latent heat flux into ocean due to evaporation/condensation', 'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_lat_fprec = register_diag_field('ocean_model', 'latent_fprec_diag', diag%axesT1, Time,&
        'Latent heat flux into ocean due to melting of frozen precipitation',                      &
        'W m-2', conversion=US%QRZ_T_to_W_m2, cmor_field_name='hfsnthermds',                       &
        cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics',                  &
        cmor_long_name='Latent Heat to Melt Frozen Precipitation')

  handles%id_lat_frunoff = register_diag_field('ocean_model', 'latent_frunoff', diag%axesT1, Time, &
        'Latent heat flux into ocean due to melting of icebergs', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        cmor_field_name='hfibthermds',                                                             &
        cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics',               &
        cmor_long_name='Latent Heat to Melt Frozen Runoff/Iceberg')

  if (present(use_glc_runoff)) then
    handles%id_lat_frunoff_glc = register_diag_field('ocean_model', 'latent_frunoff_glc', diag%axesT1, Time, &
          'Latent heat flux into ocean due to melting of frozen glacier runoff', 'W m-2', conversion=US%QRZ_T_to_W_m2)
  endif

  handles%id_sens = register_diag_field('ocean_model', 'sensible', diag%axesT1, Time, &
        'Sensible heat flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2,        &
        standard_name='surface_downward_sensible_heat_flux',                         &
        cmor_field_name='hfsso',                                                     &
        cmor_standard_name='surface_downward_sensible_heat_flux',                    &
        cmor_long_name='Surface Downward Sensible Heat Flux')

  handles%id_seaice_melt_heat = register_diag_field('ocean_model', 'seaice_melt_heat', diag%axesT1, Time,&
        'Heat flux into ocean due to snow and sea ice melt/freeze', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='snow_ice_melt_heat_flux',                         &
  !GMM TODO cmor_field_name='hfsso',                                                     &
        cmor_standard_name='snow_ice_melt_heat_flux',                    &
        cmor_long_name='Heat flux into ocean from snow and sea ice melt')

  handles%id_heat_added = register_diag_field('ocean_model', 'heat_added', diag%axesT1, Time, &
        'Flux Adjustment or restoring surface heat flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2)


  !===============================================================
  ! area integrated surface heat fluxes

  handles%id_total_heat_content_frunoff = register_scalar_field('ocean_model',                     &
      'total_heat_content_frunoff', Time, diag,                                                    &
      long_name='Area integrated heat content (relative to 0C) of solid runoff',                   &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2, cmor_field_name='total_hfsolidrunoffds', &
      cmor_standard_name=                                                                          &
      'temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                              &
      'Temperature Flux due to Solid Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_lrunoff = register_scalar_field('ocean_model',               &
      'total_heat_content_lrunoff', Time, diag,                                              &
      long_name='Area integrated heat content (relative to 0C) of liquid runoff',            &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2, cmor_field_name='total_hfrunoffds', &
      cmor_standard_name=                                                                    &
      'temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                        &
      'Temperature Flux due to Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  if (present(use_glc_runoff)) then
    handles%id_total_heat_content_frunoff_glc = register_scalar_field('ocean_model',                 &
        'total_heat_content_frunoff_glc', Time, diag,                                                &
        long_name='Area integrated heat content (relative to 0C) of solid glacier runoff',           &
        units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2) ! todo: update cmor names

    handles%id_total_heat_content_lrunoff_glc = register_scalar_field('ocean_model',               &
        'total_heat_content_lrunoff_glc', Time, diag,                                              &
        long_name='Area integrated heat content (relative to 0C) of liquid glacier runoff',        &
        units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2) ! todo: update cmor names
  endif

  handles%id_total_heat_content_lprec = register_scalar_field('ocean_model',                   &
      'total_heat_content_lprec', Time, diag,                                                  &
      long_name='Area integrated heat content (relative to 0C) of liquid precip',              &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2, cmor_field_name='total_hfrainds',   &
      cmor_standard_name=                                                                      &
      'temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                          &
      'Temperature Flux due to Rainfall Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_fprec = register_scalar_field('ocean_model',     &
      'total_heat_content_fprec', Time, diag,                                    &
      long_name='Area integrated heat content (relative to 0C) of frozen precip',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_heat_content_vprec = register_scalar_field('ocean_model',      &
      'total_heat_content_vprec', Time, diag,                                     &
      long_name='Area integrated heat content (relative to 0C) of virtual precip',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_heat_content_cond = register_scalar_field('ocean_model',   &
      'total_heat_content_cond', Time, diag,                                  &
      long_name='Area integrated heat content (relative to 0C) of condensate',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_heat_content_evap = register_scalar_field('ocean_model',    &
      'total_heat_content_evap', Time, diag,                                   &
      long_name='Area integrated heat content (relative to 0C) of evaporation',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_heat_content_surfwater = register_scalar_field('ocean_model',          &
      'total_heat_content_surfwater', Time, diag,                                         &
      long_name='Area integrated heat content (relative to 0C) of water crossing surface',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_heat_content_massout = register_scalar_field('ocean_model',                      &
      'total_heat_content_massout', Time, diag,                                                     &
      long_name='Area integrated heat content (relative to 0C) of water leaving ocean',             &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                                          &
      cmor_field_name='total_hfevapds',                                                             &
      cmor_standard_name=                                                                           &
      'temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water_area_integrated',&
      cmor_long_name='Heat Flux Out of Sea Water due to Evaporating Water Area Integrated')

  handles%id_total_heat_content_massin = register_scalar_field('ocean_model',           &
      'total_heat_content_massin', Time, diag,                                          &
      long_name='Area integrated heat content (relative to 0C) of water entering ocean',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_net_heat_coupler = register_scalar_field('ocean_model',                       &
      'total_net_heat_coupler', Time, diag,                                                      &
      long_name='Area integrated surface heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_net_heat_surface = register_scalar_field('ocean_model',                      &
      'total_net_heat_surface', Time, diag,                                                     &
      long_name='Area integrated surface heat flux from SW+LW+lat+sens+mass+frazil+restore or flux adjustments', &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,   &
      cmor_field_name='total_hfds',                                                             &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_integrated',             &
      cmor_long_name=                                                                           &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Integrated')

  handles%id_total_sw = register_scalar_field('ocean_model',                                &
      'total_sw', Time, diag,                                                               &
      long_name='Area integrated net downward shortwave at sea water surface',              &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                                  &
      cmor_field_name='total_rsntds',                                                       &
      cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface_area_integrated',&
      cmor_long_name=                                                                       &
      'Net Downward Shortwave Radiation at Sea Water Surface Area Integrated')

  handles%id_total_LwLatSens = register_scalar_field('ocean_model',&
      'total_LwLatSens', Time, diag,                               &
      long_name='Area integrated longwave+latent+sensible heating',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_lw = register_scalar_field('ocean_model',                  &
      'total_lw', Time, diag,                                                 &
      long_name='Area integrated net downward longwave at sea water surface', &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,      &
      cmor_field_name='total_rlntds',                                         &
      cmor_standard_name='surface_net_downward_longwave_flux_area_integrated',&
      cmor_long_name=                                                         &
      'Surface Net Downward Longwave Radiation Area Integrated')

  handles%id_total_lat = register_scalar_field('ocean_model',                &
      'total_lat', Time, diag,                                               &
      long_name='Area integrated surface downward latent heat flux',         &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                   &
      cmor_field_name='total_hflso',                                         &
      cmor_standard_name='surface_downward_latent_heat_flux_area_integrated',&
      cmor_long_name=                                                        &
      'Surface Downward Latent Heat Flux Area Integrated')

  handles%id_total_lat_evap = register_scalar_field('ocean_model',      &
      'total_lat_evap', Time, diag,                                     &
      long_name='Area integrated latent heat flux due to evap/condense',&
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_lat_fprec = register_scalar_field('ocean_model',                            &
      'total_lat_fprec', Time, diag,                                                           &
      long_name='Area integrated latent heat flux due to melting frozen precip',               &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                                     &
      cmor_field_name='total_hfsnthermds',                                                     &
      cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics_area_integrated',&
      cmor_long_name=                                                                          &
      'Latent Heat to Melt Frozen Precipitation Area Integrated')

  handles%id_total_lat_frunoff = register_scalar_field('ocean_model',                             &
      'total_lat_frunoff', Time, diag,                                                            &
      long_name='Area integrated latent heat flux due to melting icebergs',                       &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                                        &
      cmor_field_name='total_hfibthermds',                                                        &
      cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics_area_integrated',&
      cmor_long_name=                                                                             &
      'Heat Flux into Sea Water due to Iceberg Thermodynamics Area Integrated')

  if (present(use_glc_runoff)) then
    handles%id_total_lat_frunoff_glc = register_scalar_field('ocean_model',                             &
        'total_lat_frunoff_glc', Time, diag,                                                            &
        long_name='Area integrated latent heat flux due to melting frozen glacier runoff',              &
        units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2) ! todo: update cmor names
  endif

  handles%id_total_sens = register_scalar_field('ocean_model',                 &
      'total_sens', Time, diag,                                                &
      long_name='Area integrated downward sensible heat flux',                 &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2,                     &
      cmor_field_name='total_hfsso',                                           &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_integrated',&
      cmor_long_name=                                                          &
      'Surface Downward Sensible Heat Flux Area Integrated')

  handles%id_total_heat_added = register_scalar_field('ocean_model',&
      'total_heat_adjustment', Time, diag,                               &
      long_name='Area integrated surface heat flux from restoring and/or flux adjustment',   &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  handles%id_total_seaice_melt_heat = register_scalar_field('ocean_model',&
      'total_seaice_melt_heat', Time, diag,                               &
      long_name='Area integrated surface heat flux from snow and sea ice melt',   &
      units='W', conversion=US%QRZ_T_to_W_m2*US%L_to_m**2)

  !===============================================================
  ! area averaged surface heat fluxes

  handles%id_net_heat_coupler_ga = register_scalar_field('ocean_model',                       &
      'net_heat_coupler_ga', Time, diag,                                                      &
      long_name='Area averaged surface heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
      units='W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_net_heat_surface_ga = register_scalar_field('ocean_model',                       &
      'net_heat_surface_ga', Time, diag, long_name=                                           &
      'Area averaged surface heat flux from SW+LW+lat+sens+mass+frazil+restore+seaice_melt_heat or flux adjustments', &
      units='W m-2', conversion=US%QRZ_T_to_W_m2,                                             &
      cmor_field_name='ave_hfds',                                                             &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_averaged',             &
      cmor_long_name=                                                                         &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Averaged')

  handles%id_sw_ga = register_scalar_field('ocean_model',                                 &
      'sw_ga', Time, diag,                                                                &
      long_name='Area averaged net downward shortwave at sea water surface',              &
      units='W m-2', conversion=US%QRZ_T_to_W_m2,                                         &
      cmor_field_name='ave_rsntds',                                                       &
      cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface_area_averaged',&
      cmor_long_name=                                                                     &
      'Net Downward Shortwave Radiation at Sea Water Surface Area Averaged')

  handles%id_LwLatSens_ga = register_scalar_field('ocean_model',&
      'LwLatSens_ga', Time, diag,                               &
      long_name='Area averaged longwave+latent+sensible heating',&
      units='W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_lw_ga = register_scalar_field('ocean_model',                   &
      'lw_ga', Time, diag,                                                  &
      long_name='Area averaged net downward longwave at sea water surface', &
      units='W m-2', conversion=US%QRZ_T_to_W_m2,                           &
      cmor_field_name='ave_rlntds',                                         &
      cmor_standard_name='surface_net_downward_longwave_flux_area_averaged',&
      cmor_long_name=                                                       &
      'Surface Net Downward Longwave Radiation Area Averaged')

  handles%id_lat_ga = register_scalar_field('ocean_model',                 &
      'lat_ga', Time, diag,                                                &
      long_name='Area averaged surface downward latent heat flux',         &
      units='W m-2', conversion=US%QRZ_T_to_W_m2,                          &
      cmor_field_name='ave_hflso',                                         &
      cmor_standard_name='surface_downward_latent_heat_flux_area_averaged',&
      cmor_long_name=                                                      &
      'Surface Downward Latent Heat Flux Area Averaged')

  handles%id_sens_ga = register_scalar_field('ocean_model',                  &
      'sens_ga', Time, diag,                                                 &
      long_name='Area averaged downward sensible heat flux',                 &
      units='W m-2', conversion=US%QRZ_T_to_W_m2,                            &
      cmor_field_name='ave_hfsso',                                           &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_averaged',&
      cmor_long_name=                                                        &
      'Surface Downward Sensible Heat Flux Area Averaged')


  !===============================================================
  ! maps of surface salt fluxes, virtual precip fluxes, and adjustments

  handles%id_saltflux = register_diag_field('ocean_model', 'salt_flux', diag%axesT1, Time,&
        'Net salt flux into ocean at surface (restoring + sea-ice)',                      &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                               &
        cmor_field_name='sfdsi', cmor_standard_name='downward_sea_ice_basal_salt_flux', &
        cmor_long_name='Downward Sea Ice Basal Salt Flux')

  handles%id_saltFluxIn = register_diag_field('ocean_model', 'salt_flux_in', diag%axesT1, Time, &
        'Salt flux into ocean at surface from coupler', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxAdded = register_diag_field('ocean_model', 'salt_flux_added', &
        diag%axesT1,Time,'Salt flux into ocean at surface due to restoring or flux adjustment', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxBehind = register_diag_field('ocean_model', 'salt_left_behind', &
        diag%axesT1,Time,'Salt left in ocean at surface due to ice formation', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxGlobalAdj = register_scalar_field('ocean_model',              &
        'salt_flux_global_restoring_adjustment', Time, diag,                       &
        'Adjustment needed to balance net global salt flux into ocean at surface', &
         units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_vPrecGlobalAdj = register_scalar_field('ocean_model',  &
        'vprec_global_adjustment', Time, diag,                      &
        'Adjustment needed to adjust net vprec into ocean to zero', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_netFWGlobalAdj = register_scalar_field('ocean_model',       &
        'net_fresh_water_global_adjustment', Time, diag,                 &
        'Adjustment needed to adjust net fresh water into ocean to zero',&
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxGlobalScl = register_scalar_field('ocean_model',            &
        'salt_flux_global_restoring_scaling', Time, diag,                        &
        'Scaling applied to balance net global salt flux into ocean at surface', &
        'nondim', conversion=1.0)

  handles%id_vPrecGlobalScl = register_scalar_field('ocean_model',&
        'vprec_global_scaling', Time, diag,                       &
        'Scaling applied to adjust net vprec into ocean to zero', &
        'nondim', conversion=1.0)

  handles%id_netFWGlobalScl = register_scalar_field('ocean_model',      &
        'net_fresh_water_global_scaling', Time, diag,                   &
        'Scaling applied to adjust net fresh water into ocean to zero', &
        'nondim', conversion=1.0)

  !===============================================================
  ! area integrals of surface salt fluxes

  handles%id_total_saltflux = register_scalar_field('ocean_model', 'total_salt_flux', &
      Time, diag,  long_name='Area integrated surface salt flux',           &
      units='kg s-1', conversion=1e-3*US%RZL2_to_kg*US%s_to_T,              &
      cmor_field_name='total_sfdsi',                                        &
      cmor_standard_name='downward_sea_ice_basal_salt_flux_area_integrated',&
      cmor_long_name='Downward Sea Ice Basal Salt Flux Area Integrated')

  handles%id_total_saltFluxIn = register_scalar_field('ocean_model', 'total_salt_Flux_In', &
      Time, diag, long_name='Area integrated surface salt flux at surface from coupler', &
      units='kg s-1', conversion=1e-3*US%RZL2_to_kg*US%s_to_T)

  handles%id_total_saltFluxAdded = register_scalar_field('ocean_model', 'total_salt_Flux_Added', &
      Time, diag, long_name='Area integrated surface salt flux due to restoring or flux adjustment', &
      units='kg s-1', conversion=1e-3*US%RZL2_to_kg*US%s_to_T)

  !===============================================================
  ! wave forcing diagnostics
  if (present(use_waves)) then
    if (use_waves) then
      handles%id_lamult = register_diag_field('ocean_model', 'lamult', &
        diag%axesT1, Time, long_name='Langmuir enhancement factor received from WW3', units="nondim", conversion=1.0)
    endif
  endif

end procedure register_forcing_type_diags
module procedure forcing_accumulate
  call fluxes_accumulate(flux_tmp, fluxes, G, wt2, forces)

end procedure forcing_accumulate
module procedure fluxes_accumulate
  real :: wt1  ! The relative weight of the previous fluxes [nondim]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB


  if (fluxes%dt_buoy_accum < 0) call MOM_error(FATAL, "fluxes_accumulate: "//&
     "fluxes must be initialzed before it can be augmented.")

  ! wt1 is the relative weight of the previous fluxes.
  wt1 = fluxes%dt_buoy_accum / (fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum)
  wt2 = 1.0 - wt1 ! = flux_tmp%dt_buoy_accum / (fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum)
  fluxes%dt_buoy_accum = fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum

  ! Copy over the pressure fields and accumulate averages of ustar or tau_mag, either from the forcing
  ! type or from the temporary fluxes type.
  if (present(forces)) then
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = forces%p_surf(i,j)
      fluxes%p_surf_full(i,j) = forces%p_surf_full(i,j)
    enddo ; enddo

    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = wt1*fluxes%ustar(i,j) + wt2*forces%ustar(i,j)
    enddo ; enddo ; endif
    if (associated(fluxes%tau_mag)) then ; do j=js,je ; do i=is,ie
      fluxes%tau_mag(i,j) = wt1*fluxes%tau_mag(i,j) + wt2*forces%tau_mag(i,j)
    enddo ; enddo ; endif
  else
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = flux_tmp%p_surf(i,j)
      fluxes%p_surf_full(i,j) = flux_tmp%p_surf_full(i,j)
    enddo ; enddo

    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = wt1*fluxes%ustar(i,j) + wt2*flux_tmp%ustar(i,j)
    enddo ; enddo ; endif
    if (associated(fluxes%tau_mag)) then ; do j=js,je ; do i=is,ie
      fluxes%tau_mag(i,j) = wt1*fluxes%tau_mag(i,j) + wt2*flux_tmp%tau_mag(i,j)
    enddo ; enddo ; endif
  endif

  ! Average ustar_gustless.
  if (associated(fluxes%ustar_gustless)) then
    if (fluxes%gustless_accum_bug) then
      do j=js,je ; do i=is,ie
        fluxes%ustar_gustless(i,j) = flux_tmp%ustar_gustless(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%ustar_gustless(i,j) = wt1*fluxes%ustar_gustless(i,j) + wt2*flux_tmp%ustar_gustless(i,j)
      enddo ; enddo
    endif
  endif

  if (associated(fluxes%tau_mag_gustless)) then
    do j=js,je ; do i=is,ie
      fluxes%tau_mag_gustless(i,j) = wt1*fluxes%tau_mag_gustless(i,j) + wt2*flux_tmp%tau_mag_gustless(i,j)
    enddo ; enddo
  endif

  ! Average the water, heat, and salt fluxes.
  do j=js,je ; do i=is,ie
    fluxes%evap(i,j) = wt1*fluxes%evap(i,j) + wt2*flux_tmp%evap(i,j)
    fluxes%lprec(i,j) = wt1*fluxes%lprec(i,j) + wt2*flux_tmp%lprec(i,j)
    fluxes%fprec(i,j) = wt1*fluxes%fprec(i,j) + wt2*flux_tmp%fprec(i,j)
    fluxes%vprec(i,j) = wt1*fluxes%vprec(i,j) + wt2*flux_tmp%vprec(i,j)
    fluxes%lrunoff(i,j) = wt1*fluxes%lrunoff(i,j) + wt2*flux_tmp%lrunoff(i,j)
    fluxes%frunoff(i,j) = wt1*fluxes%frunoff(i,j) + wt2*flux_tmp%frunoff(i,j)
    fluxes%lrunoff_glc(i,j) = wt1*fluxes%lrunoff_glc(i,j) + wt2*flux_tmp%lrunoff_glc(i,j)
    fluxes%frunoff_glc(i,j) = wt1*fluxes%frunoff_glc(i,j) + wt2*flux_tmp%frunoff_glc(i,j)
    fluxes%seaice_melt(i,j) = wt1*fluxes%seaice_melt(i,j) + wt2*flux_tmp%seaice_melt(i,j)
    fluxes%sw(i,j) = wt1*fluxes%sw(i,j) + wt2*flux_tmp%sw(i,j)
    fluxes%sw_vis_dir(i,j) = wt1*fluxes%sw_vis_dir(i,j) + wt2*flux_tmp%sw_vis_dir(i,j)
    fluxes%sw_vis_dif(i,j) = wt1*fluxes%sw_vis_dif(i,j) + wt2*flux_tmp%sw_vis_dif(i,j)
    fluxes%sw_nir_dir(i,j) = wt1*fluxes%sw_nir_dir(i,j) + wt2*flux_tmp%sw_nir_dir(i,j)
    fluxes%sw_nir_dif(i,j) = wt1*fluxes%sw_nir_dif(i,j) + wt2*flux_tmp%sw_nir_dif(i,j)
    fluxes%lw(i,j) = wt1*fluxes%lw(i,j) + wt2*flux_tmp%lw(i,j)
    fluxes%latent(i,j) = wt1*fluxes%latent(i,j) + wt2*flux_tmp%latent(i,j)
    fluxes%sens(i,j) = wt1*fluxes%sens(i,j) + wt2*flux_tmp%sens(i,j)

    fluxes%salt_flux(i,j) = wt1*fluxes%salt_flux(i,j) + wt2*flux_tmp%salt_flux(i,j)
  enddo ; enddo
  if (associated(fluxes%heat_added) .and. associated(flux_tmp%heat_added)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_added(i,j) = wt1*fluxes%heat_added(i,j) + wt2*flux_tmp%heat_added(i,j)
    enddo ; enddo
  endif
  ! These might always be associated, in which case they can be combined?
  if (associated(fluxes%heat_content_cond) .and. associated(flux_tmp%heat_content_cond)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_cond(i,j) = wt1*fluxes%heat_content_cond(i,j) + wt2*flux_tmp%heat_content_cond(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_evap) .and. associated(flux_tmp%heat_content_evap)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_evap(i,j) = wt1*fluxes%heat_content_evap(i,j) + wt2*flux_tmp%heat_content_evap(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_lprec) .and. associated(flux_tmp%heat_content_lprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_lprec(i,j) = wt1*fluxes%heat_content_lprec(i,j) + wt2*flux_tmp%heat_content_lprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_fprec) .and. associated(flux_tmp%heat_content_fprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_fprec(i,j) = wt1*fluxes%heat_content_fprec(i,j) + wt2*flux_tmp%heat_content_fprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_vprec) .and. associated(flux_tmp%heat_content_vprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_vprec(i,j) = wt1*fluxes%heat_content_vprec(i,j) + wt2*flux_tmp%heat_content_vprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_lrunoff) .and. associated(flux_tmp%heat_content_lrunoff)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_lrunoff(i,j) = wt1*fluxes%heat_content_lrunoff(i,j) + wt2*flux_tmp%heat_content_lrunoff(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_frunoff) .and. associated(flux_tmp%heat_content_frunoff)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_frunoff(i,j) = wt1*fluxes%heat_content_frunoff(i,j) + wt2*flux_tmp%heat_content_frunoff(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_lrunoff_glc) .and. associated(flux_tmp%heat_content_lrunoff_glc)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_lrunoff_glc(i,j) = wt1*fluxes%heat_content_lrunoff_glc(i,j) + &
                                             wt2*flux_tmp%heat_content_lrunoff_glc(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_frunoff_glc) .and. associated(flux_tmp%heat_content_frunoff_glc)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_frunoff_glc(i,j) = wt1*fluxes%heat_content_frunoff_glc(i,j) + &
                                             wt2*flux_tmp%heat_content_frunoff_glc(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%carbon_content_lrunoff) .and. associated(flux_tmp%carbon_content_lrunoff)) then
    do j=js,je ; do i=is,ie
      fluxes%carbon_content_lrunoff(i,j) = wt1*fluxes%carbon_content_lrunoff(i,j) + &
                                           wt2*flux_tmp%carbon_content_lrunoff(i,j)
    enddo ; enddo
  endif

  if (associated(fluxes%ustar_shelf) .and. associated(flux_tmp%ustar_shelf)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%ustar_shelf(i,j)  = flux_tmp%ustar_shelf(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%iceshelf_melt) .and. associated(flux_tmp%iceshelf_melt)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%iceshelf_melt(i,j)  = flux_tmp%iceshelf_melt(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%shelf_sfc_mass_flux) &
                 .and. associated(flux_tmp%shelf_sfc_mass_flux)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%shelf_sfc_mass_flux(i,j)  = flux_tmp%shelf_sfc_mass_flux(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%frac_shelf_h) .and. associated(flux_tmp%frac_shelf_h)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%frac_shelf_h(i,j)  = flux_tmp%frac_shelf_h(i,j)
    enddo ; enddo
  endif

  if (coupler_type_initialized(fluxes%tr_fluxes) .and. &
      coupler_type_initialized(flux_tmp%tr_fluxes)) &
    call coupler_type_increment_data(flux_tmp%tr_fluxes, fluxes%tr_fluxes, &
                              scale_factor=wt2, scale_prev=wt1)

end procedure fluxes_accumulate
module procedure copy_common_forcing_fields
  logical :: do_pres
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do_pres = .true. ; if (present(skip_pres)) do_pres = .not.skip_pres

  if (associated(forces%ustar) .and. associated(fluxes%ustar)) then
    do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = forces%ustar(i,j)
    enddo ; enddo
  endif
  if (associated(forces%omega_w2x) .and. associated(fluxes%omega_w2x)) then
    do j=js,je ; do i=is,ie
      fluxes%omega_w2x(i,j) = forces%omega_w2x(i,j)
    enddo ; enddo
  endif
  if (associated(forces%tau_mag) .and. associated(fluxes%tau_mag)) then
    do j=js,je ; do i=is,ie
      fluxes%tau_mag(i,j) = forces%tau_mag(i,j)
    enddo ; enddo
  endif

  if (do_pres) then
    if (associated(forces%p_surf) .and. associated(fluxes%p_surf)) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf(i,j) = forces%p_surf(i,j)
      enddo ; enddo
    endif

    if (associated(forces%p_surf_full) .and. associated(fluxes%p_surf_full)) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = forces%p_surf_full(i,j)
      enddo ; enddo
    endif

    if (associated(forces%p_surf_SSH, forces%p_surf_full)) then
      fluxes%p_surf_SSH => fluxes%p_surf_full
    elseif (associated(forces%p_surf_SSH, forces%p_surf)) then
      fluxes%p_surf_SSH => fluxes%p_surf
    endif
  endif

end procedure copy_common_forcing_fields
module procedure set_derived_forcing_fields
  real :: taux2, tauy2 ! Squared wind stress components [R2 L2 Z2 T-4 ~> Pa2].
  real :: Irho0        ! Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  Irho0 = US%L_to_Z / Rho0

  if ( associated(forces%taux) .and. associated(forces%tauy) .and. &
       (associated(fluxes%ustar_gustless) .or. associated(fluxes%tau_mag_gustless)) ) then
    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0.0) &
        taux2 = (G%mask2dCu(I-1,j) * (forces%taux(I-1,j)**2) + &
                 G%mask2dCu(I,j) * (forces%taux(I,j)**2)) / &
                (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0.0) &
        tauy2 = (G%mask2dCv(i,J-1) * (forces%tauy(i,J-1)**2) + &
                 G%mask2dCv(i,J) * (forces%tauy(i,J)**2)) / &
                (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))

      if (associated(fluxes%ustar_gustless)) then
        if (fluxes%gustless_accum_bug) then
          ! This change is just for computational efficiency, but it is wrapped with another change.
          fluxes%ustar_gustless(i,j) = sqrt(US%L_to_Z * sqrt(taux2 + tauy2) / Rho0)
        else
          fluxes%ustar_gustless(i,j) = sqrt(sqrt(taux2 + tauy2) * Irho0)
        endif
      endif
      if (associated(fluxes%tau_mag_gustless)) then
        fluxes%tau_mag_gustless(i,j) = US%L_to_Z*sqrt(taux2 + tauy2)
      endif
    enddo ; enddo
  endif

end procedure set_derived_forcing_fields
module procedure set_net_mass_forcing
  if (associated(forces%net_mass_src)) &
    call get_net_mass_forcing(fluxes, G, US, forces%net_mass_src)

end procedure set_net_mass_forcing
module procedure get_net_mass_forcing
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  net_mass_src(:,:) = 0.0
  if (associated(fluxes%lprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%lprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%fprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%fprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%vprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%vprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%lrunoff)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%lrunoff(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%frunoff)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%frunoff(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%lrunoff_glc)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%lrunoff_glc(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%frunoff_glc)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%frunoff_glc(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%evap)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%evap(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%seaice_melt)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + fluxes%seaice_melt(i,j)
  enddo ; enddo ; endif

end procedure get_net_mass_forcing
module procedure copy_back_forcing_fields
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (associated(forces%ustar) .and. associated(fluxes%ustar)) then
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = fluxes%ustar(i,j)
    enddo ; enddo
  endif
  if (associated(forces%omega_w2x) .and. associated(fluxes%omega_w2x)) then
    do j=js,je ; do i=is,ie
      forces%omega_w2x(i,j) = fluxes%omega_w2x(i,j)
    enddo ; enddo
  endif
  if (associated(forces%tau_mag) .and. associated(fluxes%tau_mag)) then
    do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = fluxes%tau_mag(i,j)
    enddo ; enddo
  endif

end procedure copy_back_forcing_fields
module procedure mech_forcing_diags
  integer :: is, ie, js, je
  type(mech_forcing), pointer :: forces
  integer :: turns
  call cpu_clock_begin(handles%id_clock_forcing)

  ! NOTE: post_data expects data to be on the rotated index map, so any
  !   rotations must be applied before saving the output.
  turns = diag%G%HI%turns
  if (turns /= 0) then
    allocate(forces)
    call allocate_mech_forcing(forces_in, diag%G, forces)
    call rotate_mech_forcing(forces_in, turns, forces)
  else
    forces => forces_in
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call enable_averages(dt, time_end, diag)
  ! if (query_averaging_enabled(diag)) then

    if ((handles%id_taux > 0) .and. associated(forces%taux)) &
      call post_data(handles%id_taux, forces%taux, diag)

    if ((handles%id_tauy > 0) .and. associated(forces%tauy)) &
      call post_data(handles%id_tauy, forces%tauy, diag)

    if ((handles%id_mass_berg > 0) .and. associated(forces%mass_berg)) &
      call post_data(handles%id_mass_berg, forces%mass_berg, diag)

    if ((handles%id_area_berg > 0) .and. associated(forces%area_berg)) &
      call post_data(handles%id_area_berg, forces%area_berg, diag)

  ! endif

  call disable_averaging(diag)

  if (turns /= 0) then
    call deallocate_mech_forcing(forces)
    deallocate(forces)
  endif

  call cpu_clock_end(handles%id_clock_forcing)
end procedure mech_forcing_diags
module procedure forcing_diagnostics
  type(ocean_grid_type), pointer :: G   ! Grid metric on model index map
  type(forcing), pointer :: fluxes      ! Fluxes on the model index map
  real, dimension(SZI_(diag%G),SZJ_(diag%G)) :: res ! A temporary array for combinations
  real :: total_mass_flux ! Diagnostic of an integrated boundary mass flux in [R Z L2 T-1 ~> kg s-1]
  real :: total_heat_flux ! Diagnostic of an integrated boundary heat flux in [Q R Z L2 T-1 ~> W]
  real :: total_salt_flux ! Diagnostic of an integrated boundary salt flux in [R Z L2 T-1 ~> kg s-1]
  real :: ave_mass_flux   ! Diagnostic of the average of a surface mass flux in [R Z T-1 ~> kg m-2 s-1]
  real :: ave_heat_flux   ! Diagnostic of the average of a surface heat flux in [Q R Z T-1 ~> W m-2]
  real :: I_dt            ! inverse time step [T-1 ~> s-1]
  integer :: turns        ! Number of index quarter turns
  logical :: mom_enthalpy ! If true (default) enthalpy terms are computed in MOM6
  integer :: i, j, is, ie, js, je
  call cpu_clock_begin(handles%id_clock_forcing)

  mom_enthalpy = .true.
  if (present(enthalpy)) mom_enthalpy = .not. enthalpy

  ! NOTE: post_data expects data to be on the rotated index map, so any
  !   rotations must be applied before saving the output.
  turns = diag%G%HI%turns
  if (turns /= 0) then
    G => diag%G
    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, G, fluxes, turns=turns)
    call rotate_forcing(fluxes_in, fluxes, turns)
  else
    G => G_in
    fluxes => fluxes_in
  endif

  I_dt    = 1.0 / fluxes%dt_buoy_accum
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call enable_averages(fluxes%dt_buoy_accum, time_end, diag)
  ! if (query_averaging_enabled(diag)) then

    ! post the diagnostics for surface mass fluxes ==================================

    if (handles%id_prcme > 0 .or. handles%id_total_prcme > 0 .or. handles%id_prcme_ga > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%lprec))       res(i,j) = res(i,j) + fluxes%lprec(i,j)
        if (associated(fluxes%fprec))       res(i,j) = res(i,j) + fluxes%fprec(i,j)
        ! fluxes%cond is not needed because it is derived from %evap > 0
        if (associated(fluxes%evap))        res(i,j) = res(i,j) + fluxes%evap(i,j)
        if (associated(fluxes%lrunoff))     res(i,j) = res(i,j) + fluxes%lrunoff(i,j)
        if (associated(fluxes%frunoff))     res(i,j) = res(i,j) + fluxes%frunoff(i,j)
        if (associated(fluxes%lrunoff_glc)) res(i,j) = res(i,j) + fluxes%lrunoff_glc(i,j)
        if (associated(fluxes%frunoff_glc)) res(i,j) = res(i,j) + fluxes%frunoff_glc(i,j)
        if (associated(fluxes%vprec))       res(i,j) = res(i,j) + fluxes%vprec(i,j)
        if (associated(fluxes%seaice_melt)) res(i,j) = res(i,j) + fluxes%seaice_melt(i,j)
      enddo ; enddo
      if (handles%id_prcme > 0) call post_data(handles%id_prcme, res, diag)
      if (handles%id_total_prcme > 0) then
        total_mass_flux = global_area_integral(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_prcme, total_mass_flux, diag)
      endif
      if (handles%id_prcme_ga > 0) then
        ave_mass_flux = global_area_mean(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_prcme_ga, ave_mass_flux, diag)
      endif
    endif

    if (handles%id_net_massout > 0 .or. handles%id_total_net_massout > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%lprec)) then
          if (fluxes%lprec(i,j) < 0.0) res(i,j) = res(i,j) + fluxes%lprec(i,j)
        endif
        if (associated(fluxes%vprec)) then
          if (fluxes%vprec(i,j) < 0.0) res(i,j) = res(i,j) + fluxes%vprec(i,j)
        endif
        if (associated(fluxes%evap)) then
          if (fluxes%evap(i,j) < 0.0) res(i,j) = res(i,j) + fluxes%evap(i,j)
        endif
        if (associated(fluxes%seaice_melt)) then
          if (fluxes%seaice_melt(i,j) < 0.0) res(i,j) = res(i,j) + fluxes%seaice_melt(i,j)
        endif
      enddo ; enddo
      if (handles%id_net_massout > 0) call post_data(handles%id_net_massout, res, diag)
      if (handles%id_total_net_massout > 0) then
        total_mass_flux = global_area_integral(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_net_massout, total_mass_flux, diag)
      endif
    endif

    if (handles%id_massout_flux > 0 .and. associated(fluxes%netMassOut)) &
      call post_data(handles%id_massout_flux, fluxes%netMassOut, diag)

    if (handles%id_net_massin > 0 .or. handles%id_total_net_massin > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%fprec)) res(i,j) = res(i,j) + fluxes%fprec(i,j)
        if (associated(fluxes%lrunoff)) res(i,j) = res(i,j) + fluxes%lrunoff(i,j)
        if (associated(fluxes%frunoff)) res(i,j) = res(i,j) + fluxes%frunoff(i,j)
        if (associated(fluxes%lrunoff_glc)) res(i,j) = res(i,j) + fluxes%lrunoff_glc(i,j)
        if (associated(fluxes%frunoff_glc)) res(i,j) = res(i,j) + fluxes%frunoff_glc(i,j)

        if (associated(fluxes%lprec)) then
          if (fluxes%lprec(i,j) > 0.0) res(i,j) = res(i,j) + fluxes%lprec(i,j)
        endif
        if (associated(fluxes%vprec)) then
          if (fluxes%vprec(i,j) > 0.0) res(i,j) = res(i,j) + fluxes%vprec(i,j)
        endif
        ! fluxes%cond is not needed because it is derived from %evap > 0
        if (associated(fluxes%evap)) then
          if (fluxes%evap(i,j) > 0.0) res(i,j) = res(i,j) + fluxes%evap(i,j)
        endif
        if (associated(fluxes%seaice_melt)) then
          if (fluxes%seaice_melt(i,j) > 0.0) res(i,j) = res(i,j) + fluxes%seaice_melt(i,j)
        endif
      enddo ; enddo
      if (handles%id_net_massin > 0) call post_data(handles%id_net_massin, res, diag)
      if (handles%id_total_net_massin > 0) then
        total_mass_flux = global_area_integral(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_net_massin, total_mass_flux, diag)
      endif
    endif

    if (handles%id_massin_flux > 0 .and. associated(fluxes%netMassIn)) &
      call post_data(handles%id_massin_flux, fluxes%netMassIn, diag)

    if ((handles%id_evap > 0) .and. associated(fluxes%evap)) &
      call post_data(handles%id_evap, fluxes%evap, diag)
    if ((handles%id_total_evap > 0) .and. associated(fluxes%evap)) then
      total_mass_flux = global_area_integral(fluxes%evap, G, tmp_scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_evap, total_mass_flux, diag)
    endif
    if ((handles%id_evap_ga > 0) .and. associated(fluxes%evap)) then
      ave_mass_flux = global_area_mean(fluxes%evap, G, tmp_scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_evap_ga, ave_mass_flux, diag)
    endif

    if (associated(fluxes%lprec) .and. associated(fluxes%fprec)) then
      do j=js,je ; do i=is,ie
        res(i,j) = fluxes%lprec(i,j) + fluxes%fprec(i,j)
      enddo ; enddo
      if (handles%id_precip > 0) call post_data(handles%id_precip, res, diag)
      if (handles%id_total_precip > 0) then
        total_mass_flux = global_area_integral(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_precip, total_mass_flux, diag)
      endif
      if (handles%id_precip_ga > 0) then
        ave_mass_flux = global_area_mean(res, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_precip_ga, ave_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%lprec)) then
      if (handles%id_lprec > 0) call post_data(handles%id_lprec, fluxes%lprec, diag)
      if (handles%id_total_lprec > 0) then
        total_mass_flux = global_area_integral(fluxes%lprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_lprec, total_mass_flux, diag)
      endif
      if (handles%id_lprec_ga > 0) then
        ave_mass_flux = global_area_mean(fluxes%lprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_lprec_ga, ave_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%fprec)) then
      if (handles%id_fprec > 0) call post_data(handles%id_fprec, fluxes%fprec, diag)
      if (handles%id_total_fprec > 0) then
        total_mass_flux = global_area_integral(fluxes%fprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_fprec, total_mass_flux, diag)
      endif
      if (handles%id_fprec_ga > 0) then
        ave_mass_flux = global_area_mean(fluxes%fprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_fprec_ga, ave_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%vprec)) then
      if (handles%id_vprec > 0) call post_data(handles%id_vprec, fluxes%vprec, diag)
      if (handles%id_total_vprec > 0) then
        total_mass_flux = global_area_integral(fluxes%vprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_vprec, total_mass_flux, diag)
      endif
      if (handles%id_vprec_ga > 0) then
        ave_mass_flux = global_area_mean(fluxes%vprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_vprec_ga, ave_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%lrunoff)) then
    if (handles%id_lrunoff > 0) call post_data(handles%id_lrunoff, fluxes%lrunoff, diag)
      if (handles%id_total_lrunoff > 0) then
        total_mass_flux = global_area_integral(fluxes%lrunoff, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_lrunoff, total_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%lrunoff_glc)) then
    if (handles%id_lrunoff_glc > 0) call post_data(handles%id_lrunoff_glc, fluxes%lrunoff_glc, diag)
      if (handles%id_total_lrunoff_glc > 0) then
        total_mass_flux = global_area_integral(fluxes%lrunoff_glc, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_lrunoff_glc, total_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%frunoff)) then
      if (handles%id_frunoff > 0) call post_data(handles%id_frunoff, fluxes%frunoff, diag)
      if (handles%id_total_frunoff > 0) then
        total_mass_flux = global_area_integral(fluxes%frunoff, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_frunoff, total_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%frunoff_glc)) then
      if (handles%id_frunoff_glc > 0) call post_data(handles%id_frunoff_glc, fluxes%frunoff_glc, diag)
      if (handles%id_total_frunoff_glc > 0) then
        total_mass_flux = global_area_integral(fluxes%frunoff_glc, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_frunoff_glc, total_mass_flux, diag)
      endif
    endif

    if (associated(fluxes%seaice_melt)) then
      if (handles%id_seaice_melt > 0) call post_data(handles%id_seaice_melt, fluxes%seaice_melt, diag)
      if (handles%id_total_seaice_melt > 0) then
        total_mass_flux = global_area_integral(fluxes%seaice_melt, G, tmp_scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_seaice_melt, total_mass_flux, diag)
      endif
    endif

    if ((handles%id_carbon_content_lrunoff > 0) .and. associated(fluxes%carbon_content_lrunoff))  &
      call post_data(handles%id_carbon_content_lrunoff, fluxes%carbon_content_lrunoff, diag)

    ! post diagnostics for boundary heat fluxes ====================================

    if ((handles%id_heat_content_lrunoff > 0) .and. associated(fluxes%heat_content_lrunoff))  &
      call post_data(handles%id_heat_content_lrunoff, fluxes%heat_content_lrunoff, diag)
    if ((handles%id_total_heat_content_lrunoff > 0) .and. associated(fluxes%heat_content_lrunoff)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_lrunoff, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_lrunoff, total_heat_flux, diag)
    endif


    if ((handles%id_heat_content_lrunoff_glc > 0) .and. associated(fluxes%heat_content_lrunoff_glc))  &
      call post_data(handles%id_heat_content_lrunoff_glc, fluxes%heat_content_lrunoff_glc, diag)
    if ((handles%id_total_heat_content_lrunoff_glc > 0) .and. associated(fluxes%heat_content_lrunoff_glc)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_lrunoff_glc, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_lrunoff_glc, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_frunoff > 0) .and. associated(fluxes%heat_content_frunoff))  &
      call post_data(handles%id_heat_content_frunoff, fluxes%heat_content_frunoff, diag)
    if ((handles%id_total_heat_content_frunoff > 0) .and. associated(fluxes%heat_content_frunoff)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_frunoff, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_frunoff, total_heat_flux, diag)
    endif
    if ((handles%id_heat_content_frunoff_glc > 0) .and. associated(fluxes%heat_content_frunoff_glc))  &
      call post_data(handles%id_heat_content_frunoff_glc, fluxes%heat_content_frunoff_glc, diag)
    if ((handles%id_total_heat_content_frunoff_glc > 0) .and. associated(fluxes%heat_content_frunoff_glc)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_frunoff_glc, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_frunoff_glc, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_lprec > 0) .and. associated(fluxes%heat_content_lprec))      &
      call post_data(handles%id_heat_content_lprec, fluxes%heat_content_lprec, diag)
    if ((handles%id_total_heat_content_lprec > 0) .and. associated(fluxes%heat_content_lprec)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_lprec, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_lprec, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_fprec > 0) .and. associated(fluxes%heat_content_fprec))      &
      call post_data(handles%id_heat_content_fprec, fluxes%heat_content_fprec, diag)
    if ((handles%id_total_heat_content_fprec > 0) .and. associated(fluxes%heat_content_fprec)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_fprec, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_fprec, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_vprec > 0) .and. associated(fluxes%heat_content_vprec))      &
      call post_data(handles%id_heat_content_vprec, fluxes%heat_content_vprec, diag)
    if ((handles%id_total_heat_content_vprec > 0) .and. associated(fluxes%heat_content_vprec)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_vprec, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_vprec, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_cond > 0) .and. associated(fluxes%heat_content_cond))        &
      call post_data(handles%id_heat_content_cond, fluxes%heat_content_cond, diag)
    if ((handles%id_total_heat_content_cond > 0) .and. associated(fluxes%heat_content_cond)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_cond, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_cond, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_evap > 0) .and. associated(fluxes%heat_content_evap))        &
      call post_data(handles%id_heat_content_evap, fluxes%heat_content_evap, diag)
    if ((handles%id_total_heat_content_evap > 0) .and. associated(fluxes%heat_content_evap)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_evap, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_evap, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_massout > 0) .and. associated(fluxes%heat_content_massout))  &
      call post_data(handles%id_heat_content_massout, fluxes%heat_content_massout, diag)
    if ((handles%id_total_heat_content_massout > 0) .and. associated(fluxes%heat_content_massout)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_massout, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_massout, total_heat_flux, diag)
    endif

    if ((handles%id_heat_content_massin > 0) .and. associated(fluxes%heat_content_massin))  &
      call post_data(handles%id_heat_content_massin, fluxes%heat_content_massin, diag)
    if ((handles%id_total_heat_content_massin > 0) .and. associated(fluxes%heat_content_massin)) then
      total_heat_flux = global_area_integral(fluxes%heat_content_massin, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_massin, total_heat_flux, diag)
    endif

    if (handles%id_net_heat_coupler > 0 .or. handles%id_total_net_heat_coupler > 0 .or. &
        handles%id_net_heat_coupler_ga > 0. ) then
      do j=js,je ; do i=is,ie
      res(i,j) = 0.0
      if (associated(fluxes%LW))               res(i,j) = res(i,j) + fluxes%lw(i,j)
      if (associated(fluxes%latent))           res(i,j) = res(i,j) + fluxes%latent(i,j)
      if (associated(fluxes%sens))             res(i,j) = res(i,j) + fluxes%sens(i,j)
      if (associated(fluxes%SW))               res(i,j) = res(i,j) + fluxes%sw(i,j)
      if (associated(fluxes%seaice_melt_heat)) res(i,j) = res(i,j) + fluxes%seaice_melt_heat(i,j)
      enddo ; enddo
      if (handles%id_net_heat_coupler > 0) call post_data(handles%id_net_heat_coupler, res, diag)
      if (handles%id_total_net_heat_coupler > 0) then
        total_heat_flux = global_area_integral(res, G, tmp_scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_net_heat_coupler, total_heat_flux, diag)
      endif
      if (handles%id_net_heat_coupler_ga > 0) then
        ave_heat_flux = global_area_mean(res, G, tmp_scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_net_heat_coupler_ga, ave_heat_flux, diag)
      endif
    endif

    if (handles%id_net_heat_surface > 0 .or. handles%id_total_net_heat_surface > 0 .or. &
        handles%id_net_heat_surface_ga > 0. ) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%LW))               res(i,j) = res(i,j) + fluxes%lw(i,j)
        if (associated(fluxes%latent))           res(i,j) = res(i,j) + fluxes%latent(i,j)
        if (associated(fluxes%sens))             res(i,j) = res(i,j) + fluxes%sens(i,j)
        if (associated(fluxes%SW))               res(i,j) = res(i,j) + fluxes%sw(i,j)
        if (associated(fluxes%seaice_melt_heat)) res(i,j) = res(i,j) + fluxes%seaice_melt_heat(i,j)
        if (allocated(sfc_state%frazil))         res(i,j) = res(i,j) + sfc_state%frazil(i,j) * I_dt
        if (associated(fluxes%heat_content_lrunoff)) &
          res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
        if (associated(fluxes%heat_content_frunoff)) &
          res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
        if (associated(fluxes%heat_content_lrunoff_glc)) &
          res(i,j) = res(i,j) + fluxes%heat_content_lrunoff_glc(i,j)
        if (associated(fluxes%heat_content_frunoff_glc)) &
          res(i,j) = res(i,j) + fluxes%heat_content_frunoff_glc(i,j)
        if (associated(fluxes%heat_content_lprec)) &
          res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
        if (associated(fluxes%heat_content_fprec)) &
          res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
        if (associated(fluxes%heat_content_vprec)) &
          res(i,j) = res(i,j) + fluxes%heat_content_vprec(i,j)
        if (associated(fluxes%heat_content_cond)) &
          res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
        if (mom_enthalpy) then
          if (associated(fluxes%heat_content_massout)) &
            res(i,j) = res(i,j) + fluxes%heat_content_massout(i,j)
        else
          if (associated(fluxes%heat_content_evap)) &
            res(i,j) = res(i,j) + fluxes%heat_content_evap(i,j)
        endif
        if (associated(fluxes%heat_added)) res(i,j) = res(i,j) + fluxes%heat_added(i,j)
      enddo ; enddo
      if (handles%id_net_heat_surface > 0) call post_data(handles%id_net_heat_surface, res, diag)

      if (handles%id_total_net_heat_surface > 0) then
        total_heat_flux = global_area_integral(res, G, tmp_scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_net_heat_surface, total_heat_flux, diag)
      endif
      if (handles%id_net_heat_surface_ga > 0) then
        ave_heat_flux = global_area_mean(res, G, tmp_scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_net_heat_surface_ga, ave_heat_flux, diag)
      endif
    endif

    if (handles%id_heat_content_surfwater > 0 .or. handles%id_total_heat_content_surfwater > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%heat_content_lrunoff))     res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
        if (associated(fluxes%heat_content_frunoff))     res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
        if (associated(fluxes%heat_content_lrunoff_glc)) res(i,j) = res(i,j) + fluxes%heat_content_lrunoff_glc(i,j)
        if (associated(fluxes%heat_content_frunoff_glc)) res(i,j) = res(i,j) + fluxes%heat_content_frunoff_glc(i,j)
        if (associated(fluxes%heat_content_lprec))       res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
        if (associated(fluxes%heat_content_fprec))       res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
        if (associated(fluxes%heat_content_vprec))       res(i,j) = res(i,j) + fluxes%heat_content_vprec(i,j)
        if (associated(fluxes%heat_content_cond))        res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
        if (mom_enthalpy) then
          if (associated(fluxes%heat_content_massout)) res(i,j) = res(i,j) + fluxes%heat_content_massout(i,j)
        else
          if (associated(fluxes%heat_content_evap))    res(i,j) = res(i,j) + fluxes%heat_content_evap(i,j)
        endif
      enddo ; enddo
      if (handles%id_heat_content_surfwater > 0) call post_data(handles%id_heat_content_surfwater, res, diag)
      if (handles%id_total_heat_content_surfwater > 0) then
        total_heat_flux = global_area_integral(res, G, tmp_scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_heat_content_surfwater, total_heat_flux, diag)
      endif
    endif

    ! for OMIP, hfrunoffds = heat content of liquid plus frozen runoff
    if (handles%id_hfrunoffds > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%heat_content_lrunoff)) res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
        if (associated(fluxes%heat_content_frunoff)) res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
        if (associated(fluxes%heat_content_lrunoff_glc)) res(i,j) = res(i,j) + fluxes%heat_content_lrunoff_glc(i,j)
        if (associated(fluxes%heat_content_frunoff_glc)) res(i,j) = res(i,j) + fluxes%heat_content_frunoff_glc(i,j)
      enddo ; enddo
      call post_data(handles%id_hfrunoffds, res, diag)
    endif

    ! for OMIP, hfrainds = heat content of lprec + fprec + cond
    if (handles%id_hfrainds > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%heat_content_lprec)) res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
        if (associated(fluxes%heat_content_fprec)) res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
        if (associated(fluxes%heat_content_cond)) res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
      enddo ; enddo
      call post_data(handles%id_hfrainds, res, diag)
    endif

    if ((handles%id_LwLatSens > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = (fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)
      enddo ; enddo
      call post_data(handles%id_LwLatSens, res, diag)
    endif

    if ((handles%id_total_LwLatSens > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = (fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)
      enddo ; enddo
      total_heat_flux = global_area_integral(res, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_LwLatSens, total_heat_flux, diag)
    endif

    if ((handles%id_LwLatSens_ga > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j))
      enddo ; enddo
      ave_heat_flux = global_area_mean(res, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_LwLatSens_ga, ave_heat_flux, diag)
    endif

    if ((handles%id_sw > 0) .and. associated(fluxes%sw)) then
      call post_data(handles%id_sw, fluxes%sw, diag)
    endif
    if ((handles%id_sw_vis > 0) .and. associated(fluxes%sw_vis_dir) .and. &
        associated(fluxes%sw_vis_dif)) then
      call post_data(handles%id_sw_vis, fluxes%sw_vis_dir+fluxes%sw_vis_dif, diag)
    endif
    if ((handles%id_sw_nir > 0) .and. associated(fluxes%sw_nir_dir) .and. &
        associated(fluxes%sw_nir_dif)) then
      call post_data(handles%id_sw_nir, fluxes%sw_nir_dir+fluxes%sw_nir_dif, diag)
    endif
    if ((handles%id_total_sw > 0) .and. associated(fluxes%sw)) then
      total_heat_flux = global_area_integral(fluxes%sw, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_sw, total_heat_flux, diag)
    endif
    if ((handles%id_sw_ga > 0) .and. associated(fluxes%sw)) then
      ave_heat_flux = global_area_mean(fluxes%sw, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_sw_ga, ave_heat_flux, diag)
    endif

    if ((handles%id_lw > 0) .and. associated(fluxes%lw)) then
      call post_data(handles%id_lw, fluxes%lw, diag)
    endif
    if ((handles%id_total_lw > 0) .and. associated(fluxes%lw)) then
      total_heat_flux = global_area_integral(fluxes%lw, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lw, total_heat_flux, diag)
    endif
    if ((handles%id_lw_ga > 0) .and. associated(fluxes%lw)) then
      ave_heat_flux = global_area_mean(fluxes%lw, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_lw_ga, ave_heat_flux, diag)
    endif

    if ((handles%id_lat > 0) .and. associated(fluxes%latent)) then
      call post_data(handles%id_lat, fluxes%latent, diag)
    endif
    if ((handles%id_total_lat > 0) .and. associated(fluxes%latent)) then
      total_heat_flux = global_area_integral(fluxes%latent, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat, total_heat_flux, diag)
    endif
    if ((handles%id_lat_ga > 0) .and. associated(fluxes%latent)) then
      ave_heat_flux = global_area_mean(fluxes%latent, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_lat_ga, ave_heat_flux, diag)
    endif

    if ((handles%id_lat_evap > 0) .and. associated(fluxes%latent_evap_diag)) then
      call post_data(handles%id_lat_evap, fluxes%latent_evap_diag, diag)
    endif
    if ((handles%id_total_lat_evap > 0) .and. associated(fluxes%latent_evap_diag)) then
      total_heat_flux = global_area_integral(fluxes%latent_evap_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_evap, total_heat_flux, diag)
    endif

    if ((handles%id_lat_fprec > 0) .and. associated(fluxes%latent_fprec_diag)) then
      call post_data(handles%id_lat_fprec, fluxes%latent_fprec_diag, diag)
    endif
    if ((handles%id_total_lat_fprec > 0) .and. associated(fluxes%latent_fprec_diag)) then
      total_heat_flux = global_area_integral(fluxes%latent_fprec_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_fprec, total_heat_flux, diag)
    endif

    if ((handles%id_lat_frunoff > 0) .and. associated(fluxes%latent_frunoff_diag)) then
      call post_data(handles%id_lat_frunoff, fluxes%latent_frunoff_diag, diag)
    endif
    if (handles%id_total_lat_frunoff > 0 .and. associated(fluxes%latent_frunoff_diag)) then
      total_heat_flux = global_area_integral(fluxes%latent_frunoff_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_frunoff, total_heat_flux, diag)
    endif

    if ((handles%id_lat_frunoff_glc > 0) .and. associated(fluxes%latent_frunoff_glc_diag)) then
      call post_data(handles%id_lat_frunoff_glc, fluxes%latent_frunoff_glc_diag, diag)
    endif
    if (handles%id_total_lat_frunoff_glc > 0 .and. associated(fluxes%latent_frunoff_glc_diag)) then
      total_heat_flux = global_area_integral(fluxes%latent_frunoff_glc_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_frunoff_glc, total_heat_flux, diag)
    endif

    if ((handles%id_sens > 0) .and. associated(fluxes%sens)) then
      call post_data(handles%id_sens, fluxes%sens, diag)
    endif

    if ((handles%id_seaice_melt_heat > 0) .and. associated(fluxes%seaice_melt_heat)) then
      call post_data(handles%id_seaice_melt_heat, fluxes%seaice_melt_heat, diag)
    endif

    if ((handles%id_total_seaice_melt_heat > 0) .and. associated(fluxes%seaice_melt_heat)) then
      total_heat_flux = global_area_integral(fluxes%seaice_melt_heat, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_seaice_melt_heat, total_heat_flux, diag)
    endif

    if ((handles%id_total_sens > 0) .and. associated(fluxes%sens)) then
      total_heat_flux = global_area_integral(fluxes%sens, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_sens, total_heat_flux, diag)
    endif
    if ((handles%id_sens_ga > 0) .and. associated(fluxes%sens)) then
      ave_heat_flux = global_area_mean(fluxes%sens, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_sens_ga, ave_heat_flux, diag)
    endif

    if ((handles%id_heat_added > 0) .and. associated(fluxes%heat_added)) then
      call post_data(handles%id_heat_added, fluxes%heat_added, diag)
    endif

    if ((handles%id_total_heat_added > 0) .and. associated(fluxes%heat_added)) then
      total_heat_flux = global_area_integral(fluxes%heat_added, G, tmp_scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_added, total_heat_flux, diag)
    endif


    ! post the diagnostics for boundary salt fluxes ==========================

    if ((handles%id_saltflux > 0) .and. associated(fluxes%salt_flux)) &
      call post_data(handles%id_saltflux, fluxes%salt_flux, diag)
    if ((handles%id_total_saltflux > 0) .and. associated(fluxes%salt_flux)) then
      total_salt_flux = global_area_integral(fluxes%salt_flux, G, tmp_scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltflux, total_salt_flux, diag)
    endif

    if ((handles%id_saltFluxAdded > 0) .and. associated(fluxes%salt_flux_added)) &
      call post_data(handles%id_saltFluxAdded, fluxes%salt_flux_added, diag)
    if ((handles%id_total_saltFluxAdded > 0) .and. associated(fluxes%salt_flux_added)) then
      total_salt_flux = global_area_integral(fluxes%salt_flux_added, G, tmp_scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltFluxAdded, total_salt_flux, diag)
    endif

    if (handles%id_saltFluxIn > 0 .and. associated(fluxes%salt_flux_in)) &
      call post_data(handles%id_saltFluxIn, fluxes%salt_flux_in, diag)
    if ((handles%id_total_saltFluxIn > 0) .and. associated(fluxes%salt_flux_in)) then
      total_salt_flux = global_area_integral(fluxes%salt_flux_in, G, tmp_scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltFluxIn, total_salt_flux, diag)
    endif

    if (handles%id_saltFluxBehind > 0 .and. associated(fluxes%salt_left_behind)) &
      call post_data(handles%id_saltFluxBehind, fluxes%salt_left_behind, diag)

    if (handles%id_saltFluxGlobalAdj > 0)                                            &
      call post_data(handles%id_saltFluxGlobalAdj, fluxes%saltFluxGlobalAdj, diag)
    if (handles%id_vPrecGlobalAdj > 0)                                               &
      call post_data(handles%id_vPrecGlobalAdj, fluxes%vPrecGlobalAdj, diag)
    if (handles%id_netFWGlobalAdj > 0)                                               &
      call post_data(handles%id_netFWGlobalAdj, fluxes%netFWGlobalAdj, diag)
    if (handles%id_saltFluxGlobalScl > 0)                                            &
      call post_data(handles%id_saltFluxGlobalScl, fluxes%saltFluxGlobalScl, diag)
    if (handles%id_vPrecGlobalScl > 0)                                               &
      call post_data(handles%id_vPrecGlobalScl, fluxes%vPrecGlobalScl, diag)
    if (handles%id_netFWGlobalScl > 0)                                               &
      call post_data(handles%id_netFWGlobalScl, fluxes%netFWGlobalScl, diag)

    ! post diagnostics related to tracer surface fluxes  ========================

    if ((handles%id_ice_fraction > 0) .and. associated(fluxes%ice_fraction)) &
      call post_data(handles%id_ice_fraction, fluxes%ice_fraction, diag)

    if ((handles%id_u10_sqr > 0) .and. associated(fluxes%u10_sqr)) &
      call post_data(handles%id_u10_sqr, fluxes%u10_sqr, diag)

    ! remaining boundary terms ==================================================

    if ((handles%id_psurf > 0) .and. associated(fluxes%p_surf))                      &
      call post_data(handles%id_psurf, fluxes%p_surf, diag)

    if ((handles%id_TKE_tidal > 0) .and. associated(fluxes%BBL_tidal_dis))    &
      call post_data(handles%id_TKE_tidal, fluxes%BBL_tidal_dis, diag)

    if ((handles%id_buoy > 0) .and. associated(fluxes%buoy))                         &
      call post_data(handles%id_buoy, fluxes%buoy, diag)

    if ((handles%id_tau_mag > 0) .and. associated(fluxes%tau_mag)) &
      call post_data(handles%id_tau_mag, fluxes%tau_mag, diag)

    if ((handles%id_ustar > 0) .and. associated(fluxes%ustar)) &
      call post_data(handles%id_ustar, fluxes%ustar, diag)

    if ((handles%id_omega_w2x > 0) .and. associated(fluxes%omega_w2x)) &
      call post_data(handles%id_omega_w2x, fluxes%omega_w2x, diag)

    if ((handles%id_ustar_berg > 0) .and. associated(fluxes%ustar_berg)) &
      call post_data(handles%id_ustar_berg, fluxes%ustar_berg, diag)

    if ((handles%id_frac_ice_cover > 0) .and. associated(fluxes%frac_shelf_h)) &
      call post_data(handles%id_frac_ice_cover, fluxes%frac_shelf_h, diag)

    if ((handles%id_ustar_ice_cover > 0) .and. associated(fluxes%ustar_shelf)) &
      call post_data(handles%id_ustar_ice_cover, fluxes%ustar_shelf, diag)

    ! wave forcing ===============================================================
    if (handles%id_lamult > 0)                                            &
      call post_data(handles%id_lamult, fluxes%lamult, diag)

  ! endif  ! query_averaging_enabled
  call disable_averaging(diag)

  if (turns /= 0) then
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif

  call cpu_clock_end(handles%id_clock_forcing)
end procedure forcing_diagnostics
module procedure allocate_forcing_by_group
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: shelf_sfc_acc, enthalpy_mom
  enthalpy_mom = .true.
  if (present (hevap)) enthalpy_mom = .not. hevap

  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB

  shelf_sfc_acc = .false.
  if (present(shelf_sfc_accumulation)) shelf_sfc_acc = shelf_sfc_accumulation

  call myAlloc(fluxes%ustar,isd,ied,jsd,jed, ustar)
  call myAlloc(fluxes%ustar_gustless,isd,ied,jsd,jed, ustar)
  call myAlloc(fluxes%tau_mag,isd,ied,jsd,jed, ustar)

  ! Note that myAlloc can be called safely multiple times for the same pointer.
  call myAlloc(fluxes%tau_mag,isd,ied,jsd,jed, tau_mag)
  call myAlloc(fluxes%tau_mag_gustless,isd,ied,jsd,jed, tau_mag)

  call myAlloc(fluxes%evap,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%lprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%fprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%vprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%lrunoff,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%frunoff,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%lrunoff_glc,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%frunoff_glc,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%seaice_melt,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%netMassOut,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%netMassIn,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%seaice_melt_heat,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%sw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%lw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%sens,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_evap_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_fprec_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_frunoff_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_frunoff_glc_diag,isd,ied,jsd,jed, heat)

  call myAlloc(fluxes%salt_flux,isd,ied,jsd,jed, salt)
  call myAlloc(fluxes%carbon_content_lrunoff,isd,ied,jsd,jed, carbon)

  if (present(heat) .and. present(water)) then ; if (heat .and. water) then
    call myAlloc(fluxes%heat_content_cond,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_evap,isd,ied,jsd,jed, .not. enthalpy_mom)
    call myAlloc(fluxes%heat_content_lprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_fprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_vprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lrunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_frunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lrunoff_glc,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_frunoff_glc,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_massout,isd,ied,jsd,jed, enthalpy_mom)
    call myAlloc(fluxes%heat_content_massin,isd,ied,jsd,jed,  enthalpy_mom)
  endif ; endif

  call myAlloc(fluxes%p_surf,isd,ied,jsd,jed, press)

  ! These fields should only be allocated if ice shelf is enabled.
  if (present(shelf)) then ; if (shelf) then
    call myAlloc(fluxes%frac_shelf_h,isd,ied,jsd,jed, shelf)
    call myAlloc(fluxes%ustar_shelf,isd,ied,jsd,jed, shelf)
    call myAlloc(fluxes%iceshelf_melt,isd,ied,jsd,jed, shelf)
    if (shelf_sfc_acc) call myAlloc(fluxes%shelf_sfc_mass_flux,isd,ied,jsd,jed, shelf_sfc_acc)
  endif ; endif

  !These fields should only be allocated when iceberg area is being passed through the coupler.
  call myAlloc(fluxes%ustar_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%area_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%mass_berg,isd,ied,jsd,jed, iceberg)

  !These fields should only be allocated when USE_CFC_CAP is activated.
  call myAlloc(fluxes%ice_fraction,isd,ied,jsd,jed, cfc)
  call myAlloc(fluxes%u10_sqr,isd,ied,jsd,jed, cfc)

  !These fields should only be allocated when wave coupling is activated.
  call myAlloc(fluxes%ice_fraction,isd,ied,jsd,jed, waves)
  call myAlloc(fluxes%lamult,isd,ied,jsd,jed, lamult)

  if (present(fix_accum_bug)) fluxes%gustless_accum_bug = .not.fix_accum_bug

  !These fields should only be allocated when USE_MARBL is activated.
  call myAlloc(fluxes%ice_fraction,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%u10_sqr,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%noy_dep,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%nhx_dep,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%atm_co2,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%atm_alt_co2,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%dust_flux,isd,ied,jsd,jed, marbl)
  call myAlloc(fluxes%iron_flux,isd,ied,jsd,jed, marbl)

  ! These fields should only be allocated when receiving multiple ice categories
  if (present(ice_ncat)) then
    call myAlloc(fluxes%fracr_cat,isd,ied,jsd,jed,1,ice_ncat+1, ice_ncat > 0)
    call myAlloc(fluxes%qsw_cat,isd,ied,jsd,jed,1,ice_ncat+1, ice_ncat > 0)
  endif

end procedure allocate_forcing_by_group
module procedure allocate_forcing_by_ref
  logical :: do_ustar, do_taumag, do_water, do_heat, do_salt, do_press, do_shelf
  logical :: do_iceberg, do_heat_added, do_buoy, do_carbon
  logical :: even_turns  ! True if turns is absent or even
  call get_forcing_groups(fluxes_ref, do_water, do_heat, do_ustar, do_taumag, do_press, &
      do_shelf, do_iceberg, do_salt, do_heat_added, do_buoy, do_carbon)

  call allocate_forcing_type(G, fluxes, do_water, do_heat, do_ustar, &
      do_press, do_shelf, do_iceberg, do_salt, tau_mag=do_taumag, carbon=do_carbon)

  ! The following fluxes would typically be allocated by the driver
  call myAlloc(fluxes%sw_vis_dir, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_vis_dir))
  call myAlloc(fluxes%sw_vis_dif, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_vis_dif))
  call myAlloc(fluxes%sw_nir_dir, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_nir_dir))
  call myAlloc(fluxes%sw_nir_dif, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_nir_dif))

  call myAlloc(fluxes%salt_flux_in, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%salt_flux_in))
  call myAlloc(fluxes%salt_flux_added, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%salt_flux_added))

  call myAlloc(fluxes%p_surf_full, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%p_surf_full))

  call myAlloc(fluxes%heat_added, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%heat_added))
  call myAlloc(fluxes%buoy, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%buoy))

  call myAlloc(fluxes%BBL_tidal_dis, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%BBL_tidal_dis))
  call myAlloc(fluxes%ustar_tidal, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%ustar_tidal))

  ! This flag would normally be set by a control flag in allocate_forcing_type.
  ! Here we copy the flag from the reference forcing.
  fluxes%gustless_accum_bug = fluxes_ref%gustless_accum_bug

  if (coupler_type_initialized(fluxes_ref%tr_fluxes)) then
    ! The data fields in the coupler_2d_bc_type are never rotated.
    even_turns = .true. ; if (present(turns)) even_turns = (modulo(turns, 2) == 0)
    if (even_turns) then
      call coupler_type_spawn(fluxes_ref%tr_fluxes, fluxes%tr_fluxes, &
                (/G%isc,G%isc,G%iec,G%iec/), (/G%jsc,G%jsc,G%jec,G%jec/))
    else
      call coupler_type_spawn(fluxes_ref%tr_fluxes, fluxes%tr_fluxes, &
                (/G%jsc,G%jsc,G%jec,G%jec/), (/G%isc,G%isc,G%iec,G%iec/))
    endif
  endif

end procedure allocate_forcing_by_ref
module procedure allocate_mech_forcing_by_group
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB

  call myAlloc(forces%taux,IsdB,IedB,jsd,jed, stress)
  call myAlloc(forces%tauy,isd,ied,JsdB,JedB, stress)

  call myAlloc(forces%ustar,isd,ied,jsd,jed, ustar)
  call myAlloc(forces%tau_mag,isd,ied,jsd,jed, ustar)
  ! Note that myAlloc can be called safely multiple times for the same pointer.
  call myAlloc(forces%tau_mag,isd,ied,jsd,jed, tau_mag)

  call myAlloc(forces%p_surf,isd,ied,jsd,jed, press)
  call myAlloc(forces%p_surf_full,isd,ied,jsd,jed, press)
  call myAlloc(forces%net_mass_src,isd,ied,jsd,jed, press)

  call myAlloc(forces%rigidity_ice_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(forces%rigidity_ice_v,isd,ied,JsdB,JedB, shelf)
  call myAlloc(forces%frac_shelf_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(forces%frac_shelf_v,isd,ied,JsdB,JedB, shelf)

  !These fields should only on allocated when iceberg area is being passed through the coupler.
  call myAlloc(forces%area_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(forces%mass_berg,isd,ied,jsd,jed, iceberg)

  !These fields should only be allocated when waves
  if (present(waves)) then ; if (waves) then
    if (.not. present(num_stk_bands)) then
      call MOM_error(FATAL,"Requested to &
      &initialize with waves, but no waves are present.")
    endif
    if (num_stk_bands > 0) then
      if (.not.associated(forces%ustkb)) then
        allocate(forces%stk_wavenumbers(num_stk_bands), source=0.0)
        allocate(forces%ustkb(isd:ied,jsd:jed,num_stk_bands), source=0.0)
        allocate(forces%vstkb(isd:ied,jsd:jed,num_stk_bands), source=0.0)
      endif
    endif
  endif ; endif

end procedure allocate_mech_forcing_by_group
module procedure allocate_mech_forcing_from_ref
  logical :: do_stress, do_ustar, do_tau_mag, do_shelf, do_press, do_iceberg
  call get_mech_forcing_groups(forces_ref, do_stress, do_ustar, do_tau_mag, do_shelf, &
                               do_press, do_iceberg)

  call allocate_mech_forcing(G, forces, do_stress, do_ustar, do_shelf, &
                             do_press, do_iceberg, tau_mag=do_tau_mag)
end procedure allocate_mech_forcing_from_ref
module procedure get_forcing_groups
  ustar = associated(fluxes%ustar) .and. associated(fluxes%ustar_gustless)
  tau_mag = associated(fluxes%tau_mag) .and. associated(fluxes%tau_mag_gustless)
  ! TODO: Check for all associated fields, but for now just check one as a marker
  water = associated(fluxes%evap)
  heat = associated(fluxes%seaice_melt_heat)
  salt = associated(fluxes%salt_flux)
  press = associated(fluxes%p_surf)
  shelf = associated(fluxes%frac_shelf_h)
  iceberg = associated(fluxes%ustar_berg)
  heat_added = associated(fluxes%heat_added)
  buoy = associated(fluxes%buoy)
  if (present(carbon)) carbon = associated(fluxes%carbon_content_lrunoff)
end procedure get_forcing_groups
module procedure get_mech_forcing_groups
  stress = associated(forces%taux) &
      .and. associated(forces%tauy)
  ustar = associated(forces%ustar)
  tau_mag = associated(forces%tau_mag)
  shelf = associated(forces%rigidity_ice_u) &
      .and. associated(forces%rigidity_ice_v) &
      .and. associated(forces%frac_shelf_u) &
      .and. associated(forces%frac_shelf_v)
  press = associated(forces%p_surf) &
      .and. associated(forces%p_surf_full) &
      .and. associated(forces%net_mass_src)
  iceberg = associated(forces%area_berg) &
      .and. associated(forces%mass_berg)
end procedure get_mech_forcing_groups
module procedure myAlloc_2d
  if (present(flag)) then ; if (flag) then ; if (.not.associated(array)) then
    allocate(array(is:ie,js:je), source=0.0)
  endif ; endif ; endif
end procedure myAlloc_2d
module procedure myAlloc_3d
  if (present(flag)) then ; if (flag) then ; if (.not.associated(array)) then
    allocate(array(is:ie,js:je,ks:ke), source=0.0)
  endif ; endif ; endif
end procedure myAlloc_3d
module procedure deallocate_forcing_type
  if (associated(fluxes%omega_w2x))            deallocate(fluxes%omega_w2x)
  if (associated(fluxes%ustar))                deallocate(fluxes%ustar)
  if (associated(fluxes%ustar_gustless))       deallocate(fluxes%ustar_gustless)
  if (associated(fluxes%tau_mag))              deallocate(fluxes%tau_mag)
  if (associated(fluxes%buoy))                 deallocate(fluxes%buoy)
  if (associated(fluxes%sw))                   deallocate(fluxes%sw)
  if (associated(fluxes%seaice_melt_heat))     deallocate(fluxes%seaice_melt_heat)
  if (associated(fluxes%sw_vis_dir))           deallocate(fluxes%sw_vis_dir)
  if (associated(fluxes%sw_vis_dif))           deallocate(fluxes%sw_vis_dif)
  if (associated(fluxes%sw_nir_dir))           deallocate(fluxes%sw_nir_dir)
  if (associated(fluxes%sw_nir_dif))           deallocate(fluxes%sw_nir_dif)
  if (associated(fluxes%lw))                   deallocate(fluxes%lw)
  if (associated(fluxes%latent))               deallocate(fluxes%latent)
  if (associated(fluxes%latent_evap_diag))     deallocate(fluxes%latent_evap_diag)
  if (associated(fluxes%latent_fprec_diag))    deallocate(fluxes%latent_fprec_diag)
  if (associated(fluxes%latent_frunoff_diag))  deallocate(fluxes%latent_frunoff_diag)
  if (associated(fluxes%latent_frunoff_glc_diag))  deallocate(fluxes%latent_frunoff_glc_diag)
  if (associated(fluxes%sens))                 deallocate(fluxes%sens)
  if (associated(fluxes%carbon_content_lrunoff)) deallocate(fluxes%carbon_content_lrunoff)
  if (associated(fluxes%heat_added))           deallocate(fluxes%heat_added)
  if (associated(fluxes%heat_content_lrunoff)) deallocate(fluxes%heat_content_lrunoff)
  if (associated(fluxes%heat_content_frunoff)) deallocate(fluxes%heat_content_frunoff)
  if (associated(fluxes%heat_content_lrunoff_glc)) deallocate(fluxes%heat_content_lrunoff_glc)
  if (associated(fluxes%heat_content_frunoff_glc)) deallocate(fluxes%heat_content_frunoff_glc)
  if (associated(fluxes%heat_content_lprec))   deallocate(fluxes%heat_content_lprec)
  if (associated(fluxes%heat_content_fprec))   deallocate(fluxes%heat_content_fprec)
  if (associated(fluxes%heat_content_cond))    deallocate(fluxes%heat_content_cond)
  if (associated(fluxes%heat_content_evap))    deallocate(fluxes%heat_content_evap)
  if (associated(fluxes%heat_content_massout)) deallocate(fluxes%heat_content_massout)
  if (associated(fluxes%heat_content_massin))  deallocate(fluxes%heat_content_massin)
  if (associated(fluxes%evap))                 deallocate(fluxes%evap)
  if (associated(fluxes%lprec))                deallocate(fluxes%lprec)
  if (associated(fluxes%fprec))                deallocate(fluxes%fprec)
  if (associated(fluxes%vprec))                deallocate(fluxes%vprec)
  if (associated(fluxes%lrunoff))              deallocate(fluxes%lrunoff)
  if (associated(fluxes%frunoff))              deallocate(fluxes%frunoff)
  if (associated(fluxes%lrunoff_glc))          deallocate(fluxes%lrunoff_glc)
  if (associated(fluxes%frunoff_glc))          deallocate(fluxes%frunoff_glc)
  if (associated(fluxes%seaice_melt))          deallocate(fluxes%seaice_melt)
  if (associated(fluxes%netMassOut))           deallocate(fluxes%netMassOut)
  if (associated(fluxes%netMassIn))            deallocate(fluxes%netMassIn)
  if (associated(fluxes%salt_flux))            deallocate(fluxes%salt_flux)
  if (associated(fluxes%p_surf_full))          deallocate(fluxes%p_surf_full)
  if (associated(fluxes%p_surf))               deallocate(fluxes%p_surf)
  if (associated(fluxes%BBL_tidal_dis))        deallocate(fluxes%BBL_tidal_dis)
  if (associated(fluxes%ustar_tidal))          deallocate(fluxes%ustar_tidal)
  if (associated(fluxes%ustar_shelf))          deallocate(fluxes%ustar_shelf)
  if (associated(fluxes%iceshelf_melt))        deallocate(fluxes%iceshelf_melt)
  if (associated(fluxes%shelf_sfc_mass_flux)) &
                                               deallocate(fluxes%shelf_sfc_mass_flux)
  if (associated(fluxes%frac_shelf_h))         deallocate(fluxes%frac_shelf_h)
  if (associated(fluxes%ustar_berg))           deallocate(fluxes%ustar_berg)
  if (associated(fluxes%area_berg))            deallocate(fluxes%area_berg)
  if (associated(fluxes%mass_berg))            deallocate(fluxes%mass_berg)
  if (associated(fluxes%ice_fraction))         deallocate(fluxes%ice_fraction)
  if (associated(fluxes%u10_sqr))              deallocate(fluxes%u10_sqr)
  if (associated(fluxes%noy_dep))              deallocate(fluxes%noy_dep)
  if (associated(fluxes%nhx_dep))              deallocate(fluxes%nhx_dep)
  if (associated(fluxes%atm_co2))              deallocate(fluxes%atm_co2)
  if (associated(fluxes%atm_alt_co2))          deallocate(fluxes%atm_alt_co2)
  if (associated(fluxes%dust_flux))            deallocate(fluxes%dust_flux)
  if (associated(fluxes%iron_flux))            deallocate(fluxes%iron_flux)
  if (associated(fluxes%fracr_cat))            deallocate(fluxes%fracr_cat)
  if (associated(fluxes%qsw_cat))              deallocate(fluxes%qsw_cat)

  call coupler_type_destructor(fluxes%tr_fluxes)

end procedure deallocate_forcing_type
module procedure deallocate_mech_forcing
  if (associated(forces%omega_w2x))      deallocate(forces%omega_w2x)
  if (associated(forces%taux))           deallocate(forces%taux)
  if (associated(forces%tauy))           deallocate(forces%tauy)
  if (associated(forces%ustar))          deallocate(forces%ustar)
  if (associated(forces%tau_mag))        deallocate(forces%tau_mag)
  if (associated(forces%p_surf))         deallocate(forces%p_surf)
  if (associated(forces%p_surf_full))    deallocate(forces%p_surf_full)
  if (associated(forces%net_mass_src))   deallocate(forces%net_mass_src)
  if (associated(forces%rigidity_ice_u)) deallocate(forces%rigidity_ice_u)
  if (associated(forces%rigidity_ice_v)) deallocate(forces%rigidity_ice_v)
  if (associated(forces%frac_shelf_u))   deallocate(forces%frac_shelf_u)
  if (associated(forces%frac_shelf_v))   deallocate(forces%frac_shelf_v)
  if (associated(forces%area_berg))      deallocate(forces%area_berg)
  if (associated(forces%mass_berg))      deallocate(forces%mass_berg)

end procedure deallocate_mech_forcing
module procedure rotate_forcing
  logical :: do_ustar, do_taumag, do_water, do_heat, do_salt, do_press, do_shelf, &
      do_iceberg, do_heat_added, do_buoy
  call get_forcing_groups(fluxes_in, do_water, do_heat, do_ustar, do_taumag, do_press, &
      do_shelf, do_iceberg, do_salt, do_heat_added, do_buoy)

  if (associated(fluxes_in%ustar)) &
    call rotate_array(fluxes_in%ustar, turns, fluxes%ustar)
  if (associated(fluxes_in%ustar_gustless)) &
    call rotate_array(fluxes_in%ustar_gustless, turns, fluxes%ustar_gustless)

  if (associated(fluxes_in%tau_mag)) &
    call rotate_array(fluxes_in%tau_mag, turns, fluxes%tau_mag)
  if (associated(fluxes_in%tau_mag_gustless)) &
    call rotate_array(fluxes_in%tau_mag_gustless, turns, fluxes%tau_mag_gustless)

  if (do_water) then
    call rotate_array(fluxes_in%evap, turns, fluxes%evap)
    call rotate_array(fluxes_in%lprec, turns, fluxes%lprec)
    call rotate_array(fluxes_in%fprec, turns, fluxes%fprec)
    call rotate_array(fluxes_in%vprec, turns, fluxes%vprec)
    call rotate_array(fluxes_in%lrunoff, turns, fluxes%lrunoff)
    call rotate_array(fluxes_in%frunoff, turns, fluxes%frunoff)
    call rotate_array(fluxes_in%lrunoff_glc, turns, fluxes%lrunoff_glc)
    call rotate_array(fluxes_in%frunoff_glc, turns, fluxes%frunoff_glc)
    call rotate_array(fluxes_in%seaice_melt, turns, fluxes%seaice_melt)
    call rotate_array(fluxes_in%netMassOut, turns, fluxes%netMassOut)
    call rotate_array(fluxes_in%netMassIn, turns, fluxes%netMassIn)
  endif

  if (do_heat) then
    call rotate_array(fluxes_in%seaice_melt_heat, turns, fluxes%seaice_melt_heat)
    call rotate_array(fluxes_in%sw, turns, fluxes%sw)
    call rotate_array(fluxes_in%lw, turns, fluxes%lw)
    call rotate_array(fluxes_in%latent, turns, fluxes%latent)
    call rotate_array(fluxes_in%sens, turns, fluxes%sens)
    call rotate_array(fluxes_in%latent_evap_diag, turns, fluxes%latent_evap_diag)
    call rotate_array(fluxes_in%latent_fprec_diag, turns, fluxes%latent_fprec_diag)
    call rotate_array(fluxes_in%latent_frunoff_diag, turns, fluxes%latent_frunoff_diag)
    call rotate_array(fluxes_in%latent_frunoff_glc_diag, turns, fluxes%latent_frunoff_glc_diag)
  endif

  if (do_salt) then
    call rotate_array(fluxes_in%salt_flux, turns, fluxes%salt_flux)
  endif

  if (do_heat .and. do_water) then
    call rotate_array(fluxes_in%heat_content_cond, turns, fluxes%heat_content_cond)
    call rotate_array(fluxes_in%heat_content_lprec, turns, fluxes%heat_content_lprec)
    call rotate_array(fluxes_in%heat_content_fprec, turns, fluxes%heat_content_fprec)
    call rotate_array(fluxes_in%heat_content_vprec, turns, fluxes%heat_content_vprec)
    call rotate_array(fluxes_in%heat_content_lrunoff, turns, fluxes%heat_content_lrunoff)
    call rotate_array(fluxes_in%heat_content_lrunoff_glc, turns, fluxes%heat_content_lrunoff_glc)
    call rotate_array(fluxes_in%heat_content_frunoff, turns, fluxes%heat_content_frunoff)
    call rotate_array(fluxes_in%heat_content_frunoff_glc, turns, fluxes%heat_content_frunoff_glc)
    if (associated (fluxes_in%heat_content_evap))  then
      call rotate_array(fluxes_in%heat_content_evap, turns, fluxes%heat_content_evap)
    else
      call rotate_array(fluxes_in%heat_content_massout, turns, fluxes%heat_content_massout)
      call rotate_array(fluxes_in%heat_content_massin, turns, fluxes%heat_content_massin)
    endif
  endif

  if (do_press) then
    call rotate_array(fluxes_in%p_surf, turns, fluxes%p_surf)
  endif

  if (do_shelf) then
    call rotate_array(fluxes_in%frac_shelf_h, turns, fluxes%frac_shelf_h)
    call rotate_array(fluxes_in%ustar_shelf, turns, fluxes%ustar_shelf)
    call rotate_array(fluxes_in%iceshelf_melt, turns, fluxes%iceshelf_melt)
    call rotate_array(fluxes_in%shelf_sfc_mass_flux, turns, fluxes%shelf_sfc_mass_flux)
  endif

  if (do_iceberg) then
    call rotate_array(fluxes_in%ustar_berg, turns, fluxes%ustar_berg)
    call rotate_array(fluxes_in%area_berg, turns, fluxes%area_berg)
    !BGR: pretty sure the following line isn't supposed to be here.
    call rotate_array(fluxes_in%iceshelf_melt, turns, fluxes%iceshelf_melt)
  endif

  if (do_heat_added) then
    call rotate_array(fluxes_in%heat_added, turns, fluxes%heat_added)
  endif

  ! The following fields are handled by drivers rather than control flags.
  if (associated(fluxes_in%sw_vis_dir)) &
    call rotate_array(fluxes_in%sw_vis_dir, turns, fluxes%sw_vis_dir)
  if (associated(fluxes_in%sw_vis_dif)) &
    call rotate_array(fluxes_in%sw_vis_dif, turns, fluxes%sw_vis_dif)
  if (associated(fluxes_in%sw_nir_dir)) &
    call rotate_array(fluxes_in%sw_nir_dir, turns, fluxes%sw_nir_dir)
  if (associated(fluxes_in%sw_nir_dif)) &
    call rotate_array(fluxes_in%sw_nir_dif, turns, fluxes%sw_nir_dif)

  if (associated(fluxes_in%salt_flux_in)) &
    call rotate_array(fluxes_in%salt_flux_in, turns, fluxes%salt_flux_in)
  if (associated(fluxes_in%salt_flux_added)) &
    call rotate_array(fluxes_in%salt_flux_added, turns, fluxes%salt_flux_added)

  if (associated(fluxes_in%p_surf_full)) &
    call rotate_array(fluxes_in%p_surf_full, turns, fluxes%p_surf_full)

  if (associated(fluxes_in%buoy)) &
    call rotate_array(fluxes_in%buoy, turns, fluxes%buoy)

  if (associated(fluxes_in%BBL_tidal_dis)) &
    call rotate_array(fluxes_in%BBL_tidal_dis, turns, fluxes%BBL_tidal_dis)
  if (associated(fluxes_in%ustar_tidal)) &
    call rotate_array(fluxes_in%ustar_tidal, turns, fluxes%ustar_tidal)

  ! NOTE: Tracer fields are handled by FMS, so are left unrotated.  Any
  ! reads/writes to tr_fields must be appropriately rotated.
  if (coupler_type_initialized(fluxes%tr_fluxes)) then
    call coupler_type_copy_data(fluxes_in%tr_fluxes, fluxes%tr_fluxes)
  endif

  ! Scalars and flags
  fluxes%accumulate_p_surf = fluxes_in%accumulate_p_surf

  fluxes%vPrecGlobalAdj = fluxes_in%vPrecGlobalAdj
  fluxes%saltFluxGlobalAdj = fluxes_in%saltFluxGlobalAdj
  fluxes%netFWGlobalAdj = fluxes_in%netFWGlobalAdj
  fluxes%vPrecGlobalScl = fluxes_in%vPrecGlobalScl
  fluxes%saltFluxGlobalScl = fluxes_in%saltFluxGlobalScl
  fluxes%netFWGlobalScl = fluxes_in%netFWGlobalScl

  fluxes%fluxes_used = fluxes_in%fluxes_used
  fluxes%dt_buoy_accum = fluxes_in%dt_buoy_accum
  fluxes%C_p = fluxes_in%C_p
  ! NOTE: gustless_accum_bug is set during allocation

  fluxes%num_msg = fluxes_in%num_msg
  fluxes%max_msg = fluxes_in%max_msg
end procedure rotate_forcing
module procedure rotate_mech_forcing
  logical :: do_stress, do_ustar, do_tau_mag, do_shelf, do_press, do_iceberg
  call get_mech_forcing_groups(forces_in, do_stress, do_ustar, do_tau_mag, do_shelf, &
                              do_press, do_iceberg)

  if (do_stress) &
    call rotate_vector(forces_in%taux, forces_in%tauy, turns, &
        forces%taux, forces%tauy)

  if (associated(forces_in%ustar)) &
    call rotate_array(forces_in%ustar, turns, forces%ustar)
  if (associated(forces_in%tau_mag)) &
    call rotate_array(forces_in%tau_mag, turns, forces%tau_mag)

  if (do_shelf) then
    call rotate_array_pair( &
      forces_in%rigidity_ice_u, forces_in%rigidity_ice_v, turns, &
      forces%rigidity_ice_u, forces%rigidity_ice_v &
    )
    call rotate_array_pair( &
      forces_in%frac_shelf_u, forces_in%frac_shelf_v, turns, &
      forces%frac_shelf_u, forces%frac_shelf_v &
    )
  endif

  if (do_press) then
    call rotate_array(forces_in%p_surf, turns, forces%p_surf)
    call rotate_array(forces_in%p_surf_full, turns, forces%p_surf_full)
    call rotate_array(forces_in%net_mass_src, turns, forces%net_mass_src)

    ! p_surf_SSH points to either p_surf or p_surf_full
    if (associated(forces_in%p_surf_SSH, forces_in%p_surf)) then
      forces%p_surf_SSH => forces%p_surf
    else if (associated(forces_in%p_surf_SSH, forces_in%p_surf_full)) then
      forces%p_surf_SSH => forces%p_surf_full
    else
      forces%p_surf_SSH => null()
    endif
  endif

  if (do_iceberg) then
    call rotate_array(forces_in%area_berg, turns, forces%area_berg)
    call rotate_array(forces_in%mass_berg, turns, forces%mass_berg)
  endif

  ! Copy fields
  forces%dt_force_accum = forces_in%dt_force_accum
  forces%net_mass_src_set = forces_in%net_mass_src_set
  forces%accumulate_p_surf = forces_in%accumulate_p_surf
  forces%accumulate_rigidity = forces_in%accumulate_rigidity
  forces%initialized = forces_in%initialized
end procedure rotate_mech_forcing
module procedure homogenize_mech_forcing
  real :: tx_mean, ty_mean ! Mean wind stresses [R L Z T-2 ~> Pa]
  real :: tau_mag      ! The magnitude of the wind stresses [R Z2 T-2 ~> Pa]
  real :: Irho0        ! Inverse of the mean density [R-1 ~> m3 kg-1]
  logical :: do_stress, do_ustar, do_taumag, do_shelf, do_press, do_iceberg, tau2ustar
  integer :: i, j, is, ie, js, je, isB, ieB, jsB, jeB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB

  Irho0 = 1.0 / Rho0

  tau2ustar = .false.
  if (present(UpdateUstar)) tau2ustar = UpdateUstar

  call get_mech_forcing_groups(forces, do_stress, do_ustar, do_taumag, do_shelf, &
                              do_press, do_iceberg)

  if (do_stress) then
    tx_mean = global_area_mean_u(forces%taux, G, tmp_scale=US%RLZ_T2_to_Pa)
    do j=js,je ; do i=isB,ieB
      if (G%mask2dCu(I,j) > 0.0) forces%taux(I,j) = tx_mean
    enddo ; enddo
    ty_mean = global_area_mean_v(forces%tauy, G, tmp_scale=US%RLZ_T2_to_Pa)
    do j=jsB,jeB ; do i=is,ie
      if (G%mask2dCv(i,J) > 0.0) forces%tauy(i,J) = ty_mean
    enddo ; enddo
    if (tau2ustar) then
      tau_mag = US%L_to_Z*sqrt((tx_mean**2) + (ty_mean**2))
      if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        forces%tau_mag(i,j) = tau_mag
      endif ; enddo ; enddo ; endif
      if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        forces%ustar(i,j) = sqrt(tau_mag * Irho0)
      endif ; enddo ; enddo ; endif
    else
      if (associated(forces%ustar)) &
        call homogenize_field_t(forces%ustar, G, tmp_scale=US%Z_to_m*US%s_to_T)
      if (associated(forces%tau_mag)) &
        call homogenize_field_t(forces%tau_mag, G, tmp_scale=US%RLZ_T2_to_Pa*US%Z_to_L)
    endif
  else
    if (associated(forces%ustar)) &
      call homogenize_field_t(forces%ustar, G, tmp_scale=US%Z_to_m*US%s_to_T)
    if (associated(forces%tau_mag)) &
      call homogenize_field_t(forces%tau_mag, G, tmp_scale=US%RLZ_T2_to_Pa*US%Z_to_L)
  endif

  if (do_shelf) then
    call homogenize_field_u(forces%rigidity_ice_u, G, tmp_scale=US%L_T_to_m_s*US%L_to_m**2*US%L_to_Z)
    call homogenize_field_v(forces%rigidity_ice_v, G, tmp_scale=US%L_T_to_m_s*US%L_to_m**2*US%L_to_Z)
    call homogenize_field_u(forces%frac_shelf_u, G)
    call homogenize_field_v(forces%frac_shelf_v, G)
  endif

  if (do_press) then
    ! NOTE: p_surf_SSH either points to p_surf or p_surf_full
    call homogenize_field_t(forces%p_surf, G, tmp_scale=US%RL2_T2_to_Pa)
    call homogenize_field_t(forces%p_surf_full, G, tmp_scale=US%RL2_T2_to_Pa)
    call homogenize_field_t(forces%net_mass_src, G, tmp_scale=US%RZ_T_to_kg_m2s)
  endif

  if (do_iceberg) then
    call homogenize_field_t(forces%area_berg, G)
    call homogenize_field_t(forces%mass_berg, G, tmp_scale=US%RZ_to_kg_m2)
  endif

end procedure homogenize_mech_forcing
module procedure homogenize_forcing
  logical :: do_ustar, do_taumag, do_water, do_heat, do_salt, do_press, do_shelf
  logical :: do_iceberg, do_heat_added, do_buoy
  call get_forcing_groups(fluxes, do_water, do_heat, do_ustar, do_taumag, do_press, &
      do_shelf, do_iceberg, do_salt, do_heat_added, do_buoy)

  if (associated(fluxes%ustar)) &
    call homogenize_field_t(fluxes%ustar, G, tmp_scale=US%Z_to_m*US%s_to_T)
  if (associated(fluxes%ustar_gustless)) &
    call homogenize_field_t(fluxes%ustar_gustless, G, tmp_scale=US%Z_to_m*US%s_to_T)

  if (associated(fluxes%tau_mag)) &
    call homogenize_field_t(fluxes%tau_mag, G, tmp_scale=US%RLZ_T2_to_Pa*US%Z_to_L)
  if (associated(fluxes%tau_mag_gustless)) &
    call homogenize_field_t(fluxes%tau_mag_gustless, G, tmp_scale=US%RLZ_T2_to_Pa*US%Z_to_L)

  if (do_water) then
    call homogenize_field_t(fluxes%evap, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%lprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%fprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%vprec, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%lrunoff, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%frunoff, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%lrunoff_glc, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%frunoff_glc, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%seaice_melt, G, tmp_scale=US%RZ_T_to_kg_m2s)
    !  These two calls might not be needed.
    call homogenize_field_t(fluxes%netMassOut, G, tmp_scale=GV%H_to_mks)
    call homogenize_field_t(fluxes%netMassIn, G, tmp_scale=GV%H_to_mks)
    !This was removed and I don't think replaced. Not needed?
    !call homogenize_field_t(fluxes%netSalt, G)
  endif

  if (do_heat) then
    call homogenize_field_t(fluxes%seaice_melt_heat, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%sw, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%lw, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%latent, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%sens, G, tmp_scale=US%QRZ_T_to_W_m2)
    !### These are for diagnostics only and may not be needed.
    call homogenize_field_t(fluxes%latent_evap_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%latent_fprec_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%latent_frunoff_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%latent_frunoff_glc_diag, G, tmp_scale=US%QRZ_T_to_W_m2)
  endif

  if (do_salt) call homogenize_field_t(fluxes%salt_flux, G, tmp_scale=US%RZ_T_to_kg_m2s)

  if (do_heat .and. do_water) then
    call homogenize_field_t(fluxes%heat_content_cond, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_lprec, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_fprec, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_vprec, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_lrunoff, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_frunoff, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_lrunoff_glc, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_frunoff_glc, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_massout, G, tmp_scale=US%QRZ_T_to_W_m2)
    call homogenize_field_t(fluxes%heat_content_massin, G, tmp_scale=US%QRZ_T_to_W_m2)
  endif

  if (do_press) call homogenize_field_t(fluxes%p_surf, G, tmp_scale=US%RL2_T2_to_Pa)

  if (do_shelf) then
    call homogenize_field_t(fluxes%frac_shelf_h, G)
    call homogenize_field_t(fluxes%ustar_shelf, G, tmp_scale=US%Z_to_m*US%s_to_T)
    call homogenize_field_t(fluxes%iceshelf_melt, G, tmp_scale=US%RZ_T_to_kg_m2s)
    call homogenize_field_t(fluxes%shelf_sfc_mass_flux, G, tmp_scale=US%RZ_T_to_kg_m2s)
  endif

  if (do_iceberg) then
    call homogenize_field_t(fluxes%ustar_berg, G, tmp_scale=US%Z_to_m*US%s_to_T)
    call homogenize_field_t(fluxes%area_berg, G)
  endif

  if (do_heat_added) then
    call homogenize_field_t(fluxes%heat_added, G, tmp_scale=US%QRZ_T_to_W_m2)
  endif

  ! The following fields are handled by drivers rather than control flags.
  if (associated(fluxes%sw_vis_dir)) &
    call homogenize_field_t(fluxes%sw_vis_dir, G, tmp_scale=US%QRZ_T_to_W_m2)

  if (associated(fluxes%sw_vis_dif)) &
    call homogenize_field_t(fluxes%sw_vis_dif, G, tmp_scale=US%QRZ_T_to_W_m2)

  if (associated(fluxes%sw_nir_dir)) &
    call homogenize_field_t(fluxes%sw_nir_dir, G, tmp_scale=US%QRZ_T_to_W_m2)

  if (associated(fluxes%sw_nir_dif)) &
    call homogenize_field_t(fluxes%sw_nir_dif, G, tmp_scale=US%QRZ_T_to_W_m2)

  if (associated(fluxes%salt_flux_in)) &
    call homogenize_field_t(fluxes%salt_flux_in, G, tmp_scale=US%RZ_T_to_kg_m2s)

  if (associated(fluxes%salt_flux_added)) &
    call homogenize_field_t(fluxes%salt_flux_added, G, tmp_scale=US%RZ_T_to_kg_m2s)

  if (associated(fluxes%p_surf_full)) &
    call homogenize_field_t(fluxes%p_surf_full, G, tmp_scale=US%RL2_T2_to_Pa)

  if (associated(fluxes%buoy)) &
    call homogenize_field_t(fluxes%buoy, G, tmp_scale=US%L_to_m**2*US%s_to_T**3)

  if (associated(fluxes%BBL_tidal_dis)) &
    call homogenize_field_t(fluxes%BBL_tidal_dis, G, tmp_scale=US%L_to_Z**2*US%RZ3_T3_to_W_m2)

  if (associated(fluxes%ustar_tidal)) &
    call homogenize_field_t(fluxes%ustar_tidal, G, tmp_scale=US%Z_to_m*US%s_to_T)

  ! TODO: tracer flux homogenization
  ! Having a warning causes a lot of errors (each time step).
  !if (coupler_type_initialized(fluxes%tr_fluxes)) &
  !  call MOM_error(WARNING, "Homogenization of tracer BC fluxes not yet implemented.")

end procedure homogenize_forcing
module procedure homogenize_field_t
  real    :: avg   ! Global average of var, in the same units as var [A ~> a]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  avg = global_area_mean(var, G, tmp_scale=tmp_scale)
  do j=js,je ; do i=is,ie
    if (G%mask2dT(i,j) > 0.0) var(i,j) = avg
  enddo ; enddo

end procedure homogenize_field_t
module procedure homogenize_field_v
  real    :: avg   ! Global average of var, in the same units as var [A ~> a]
  integer :: i, j, is, ie, jsB, jeB
  is = G%isc ; ie = G%iec ; jsB = G%jscB ; jeB = G%jecB

  avg = global_area_mean_v(var, G, tmp_scale=tmp_scale)
  do J=jsB,jeB ; do i=is,ie
    if (G%mask2dCv(i,J) > 0.0) var(i,J) = avg
  enddo ; enddo

end procedure homogenize_field_v
module procedure homogenize_field_u
  real    :: avg   ! Global average of var, in the same units as var [A ~> a]
  integer :: i, j, isB, ieB, js, je
  isB = G%iscB ; ieB = G%iecB ; js = G%jsc ; je = G%jec

  avg = global_area_mean_u(var, G, tmp_scale=tmp_scale)
  do j=js,je ; do I=isB,ieB
    if (G%mask2dCu(I,j) > 0.0) var(I,j) = avg
  enddo ; enddo

end procedure homogenize_field_u
end submodule MOM_forcing_type_s
