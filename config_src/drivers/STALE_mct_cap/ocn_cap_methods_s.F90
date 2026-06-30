submodule (ocn_cap_methods) ocn_cap_methods_s
  implicit none
contains
module procedure ocn_import
  integer         :: i, j, isc, iec, jsc, jec  ! Grid indices
  integer         :: k
  integer         :: day, secs, rc
  type(ESMF_time) :: currTime
  character(*), parameter :: F01  = "('(ocn_import) ',a,4(i6,2x),d21.14)"
  isc = GRID%isc ; iec = GRID%iec ; jsc = GRID%jsc ; jec = GRID%jec

  k = 0
  do j = jsc, jec
    do i = isc, iec
      k = k + 1 ! Increment position within gindex

      ! rotate taux and tauy from true zonal/meridional to local coordinates
      ! taux
      ice_ocean_boundary%u_flux(i,j) = GRID%cos_rot(i,j) * x2o(ind%x2o_Foxx_taux,k) &
                                      - GRID%sin_rot(i,j) * x2o(ind%x2o_Foxx_tauy,k)

      ! tauy
      ice_ocean_boundary%v_flux(i,j) = GRID%cos_rot(i,j) * x2o(ind%x2o_Foxx_tauy,k) &
                                      + GRID%sin_rot(i,j) * x2o(ind%x2o_Foxx_taux,k)

      ! liquid precipitation (rain)
      ice_ocean_boundary%lprec(i,j) = x2o(ind%x2o_Faxa_rain,k)

      ! frozen precipitation (snow)
      ice_ocean_boundary%fprec(i,j) = x2o(ind%x2o_Faxa_snow,k)

      ! longwave radiation, sum up and down (W/m2)
      ice_ocean_boundary%lw_flux(i,j) = (x2o(ind%x2o_Faxa_lwdn,k) + x2o(ind%x2o_Foxx_lwup,k))

      ! specific humitidy flux
      ice_ocean_boundary%q_flux(i,j) = x2o(ind%x2o_Foxx_evap,k)

      ! sensible heat flux (W/m2)
      ice_ocean_boundary%t_flux(i,j) = x2o(ind%x2o_Foxx_sen,k)

      ! snow&ice melt heat flux  (W/m^2)
      ice_ocean_boundary%seaice_melt_heat(i,j) = x2o(ind%x2o_Fioi_melth,k)

      ! water flux from snow&ice melt (kg/m2/s)
      ice_ocean_boundary%seaice_melt(i,j) = x2o(ind%x2o_Fioi_meltw,k)

      ! liquid runoff
      ice_ocean_boundary%rofl_flux(i,j) = x2o(ind%x2o_Foxx_rofl,k) * GRID%mask2dT(i,j)

      ! ice runoff
      ice_ocean_boundary%rofi_flux(i,j) = x2o(ind%x2o_Foxx_rofi,k) * GRID%mask2dT(i,j)

      ! surface pressure
      ice_ocean_boundary%p(i,j) = x2o(ind%x2o_Sa_pslv,k) * GRID%mask2dT(i,j)

      ! salt flux
      ice_ocean_boundary%salt_flux(i,j) = x2o(ind%x2o_Fioi_salt,k) * GRID%mask2dT(i,j)

      ! 1) visible, direct shortwave  (W/m2)
      ! 2) visible, diffuse shortwave (W/m2)
      ! 3) near-IR, direct shortwave  (W/m2)
      ! 4) near-IR, diffuse shortwave (W/m2)
      if (present(c1) .and. present(c2) .and. present(c3) .and. present(c4)) then
        ! Use runtime coefficients to decompose net short-wave heat flux into 4 components
        ice_ocean_boundary%sw_flux_vis_dir(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c1 * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_vis_dif(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c2 * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_nir_dir(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c3 * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_nir_dif(i,j) = x2o(ind%x2o_Foxx_swnet,k) * c4 * GRID%mask2dT(i,j)
      else
        ice_ocean_boundary%sw_flux_vis_dir(i,j) = x2o(ind%x2o_Faxa_swvdr,k) * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_vis_dif(i,j) = x2o(ind%x2o_Faxa_swvdf,k) * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_nir_dir(i,j) = x2o(ind%x2o_Faxa_swndr,k) * GRID%mask2dT(i,j)
        ice_ocean_boundary%sw_flux_nir_dif(i,j) = x2o(ind%x2o_Faxa_swndf,k) * GRID%mask2dT(i,j)
      endif
    enddo
  enddo

  if (debug .and. is_root_pe()) then
    call ESMF_ClockGet(EClock, CurrTime=CurrTime, rc=rc)
    call ESMF_TimeGet(CurrTime, d=day, s=secs, rc=rc)

    do j = GRID%jsc, GRID%jec
      do i = GRID%isc, GRID%iec
        write(logunit,F01)'import: day, secs, j, i, u_flux           = ',day,secs,j,i,ice_ocean_boundary%u_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, v_flux           = ',day,secs,j,i,ice_ocean_boundary%v_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, lprec            = ',day,secs,j,i,ice_ocean_boundary%lprec(i,j)
        write(logunit,F01)'import: day, secs, j, i, lwrad            = ',day,secs,j,i,ice_ocean_boundary%lw_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, q_flux           = ',day,secs,j,i,ice_ocean_boundary%q_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, t_flux           = ',day,secs,j,i,ice_ocean_boundary%t_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, seaice_melt_heat = ',&
                          day,secs,j,i,ice_ocean_boundary%seaice_melt_heat(i,j)
        write(logunit,F01)'import: day, secs, j, i, seaice_melt      = ',&
                          day,secs,j,i,ice_ocean_boundary%seaice_melt(i,j)
        write(logunit,F01)'import: day, secs, j, i, runoff          = ',&
                          day,secs,j,i,ice_ocean_boundary%rofl_flux(i,j) + ice_ocean_boundary%rofi_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, psurf           = ',&
                          day,secs,j,i,ice_ocean_boundary%p(i,j)
        write(logunit,F01)'import: day, secs, j, i, salt_flux       = ',&
                          day,secs,j,i,ice_ocean_boundary%salt_flux(i,j)
        write(logunit,F01)'import: day, secs, j, i, sw_flux_vis_dir = ',&
                          day,secs,j,i,ice_ocean_boundary%sw_flux_vis_dir(i,j)
        write(logunit,F01)'import: day, secs, j, i, sw_flux_vis_dif = ',&
                          day,secs,j,i,ice_ocean_boundary%sw_flux_vis_dif(i,j)
        write(logunit,F01)'import: day, secs, j, i, sw_flux_nir_dir = ',&
                          day,secs,j,i,ice_ocean_boundary%sw_flux_nir_dir(i,j)
        write(logunit,F01)'import: day, secs, j, i, sw_flux_nir_dif = ',&
                          day,secs,j,i,ice_ocean_boundary%sw_flux_nir_dir(i,j)
      enddo
    enddo
  endif

end procedure ocn_import
module procedure ocn_export
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: sshx!< Zonal SSH gradient, local coordinate.
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: sshy!< Meridional SSH gradient, local coordinate.
  integer :: i, j, n, ig, jg  !< Grid indices
  real    :: slp_L, slp_R, slp_C, slope, u_min, u_max
  real :: I_time_int  !< The inverse of coupling time interval [s-1].
  I_time_int = 0.0 ; if (dt_int > 0.0) I_time_int = 1.0 / dt_int

  ! Copy from ocn_public to o2x. ocn_public uses global indexing with no halos.
  ! The mask comes from "grid" that uses the usual MOM domain that has halos
  ! and does not use global indexing.

  n = 0
  do j=grid%jsc, grid%jec
    jg = j + grid%jdg_offset
    do i=grid%isc,grid%iec
      n = n+1
      ig = i + grid%idg_offset
      ! surface temperature in Kelvin
      o2x(ind%o2x_So_t, n) = ocn_public%t_surf(ig,jg) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_s, n) = ocn_public%s_surf(ig,jg) * grid%mask2dT(i,j)
      ! rotate ocn current from local tripolar grid to true zonal/meridional (inverse transformation)
      o2x(ind%o2x_So_u, n) = (grid%cos_rot(i,j) * ocn_public%u_surf(ig,jg) + &
                              grid%sin_rot(i,j) * ocn_public%v_surf(ig,jg)) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_v, n) = (grid%cos_rot(i,j) * ocn_public%v_surf(ig,jg) - &
                              grid%sin_rot(i,j) * ocn_public%u_surf(ig,jg)) * grid%mask2dT(i,j)

      ! boundary layer depth (m)
      o2x(ind%o2x_So_bldepth, n) = ocn_public%OBLD(ig,jg) * grid%mask2dT(i,j)
      ! ocean melt and freeze potential (o2x_Fioo_q), W m-2
      if (ocn_public%frazil(ig,jg) > 0.0) then
        ! Frazil: change from J/m^2 to W/m^2
        o2x(ind%o2x_Fioo_q, n) = ocn_public%frazil(ig,jg) * grid%mask2dT(i,j) * I_time_int
      else
        ! Melt_potential: change from J/m^2 to W/m^2
        o2x(ind%o2x_Fioo_q, n) = -ocn_public%melt_potential(ig,jg) * grid%mask2dT(i,j) * I_time_int !* ncouple_per_day
        ! make sure Melt_potential is always <= 0
        if (o2x(ind%o2x_Fioo_q, n) > 0.0) o2x(ind%o2x_Fioo_q, n) = 0.0
      endif
      ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
      ! in order to update halos. i.e. does not use global indexing.
      ssh(i,j) = ocn_public%sea_lev(ig,jg)
    enddo
  enddo

  ! Update halo of ssh so we can calculate gradients
  call pass_var(ssh, grid%domain)

  ! d/dx ssh
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    ! This is a simple second-order difference
    ! o2x(ind%o2x_So_dhdx, n) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * grid%US%m_to_L*grid%IdxT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = (ssh(I,j) - ssh(I-1,j)) * grid%mask2dCu(I-1,j)
    if (grid%mask2dCu(I-1,j)==0.) slp_L = 0.
    slp_R = (ssh(I+1,j) - ssh(I,j)) * grid%mask2dCu(I,j)
    if (grid%mask2dCu(I+1,j)==0.) slp_R = 0.
    slp_C = 0.5 * (slp_L + slp_R)
    if ( (slp_L * slp_R) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      u_max = max( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    endif
    sshx(i,j) = slope * grid%US%m_to_L*grid%IdxT(i,j) * grid%mask2dT(i,j)
    if (grid%mask2dT(i,j)==0.) sshx(i,j) = 0.0
  enddo ; enddo

  ! d/dy ssh
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    ! This is a simple second-order difference
    ! o2x(ind%o2x_So_dhdy, n) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * grid%US%m_to_L*grid%IdyT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = ssh(i,J) - ssh(i,J-1) * grid%mask2dCv(i,J-1)
    if (grid%mask2dCv(i,J-1)==0.) slp_L = 0.

    slp_R = ssh(i,J+1) - ssh(i,J) * grid%mask2dCv(i,J)
    if (grid%mask2dCv(i,J+1)==0.) slp_R = 0.

    slp_C = 0.5 * (slp_L + slp_R)
    if ((slp_L * slp_R) > 0.0) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      u_max = max( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    endif
    sshy(i,j) = slope * grid%US%m_to_L*grid%IdyT(i,j) * grid%mask2dT(i,j)
    if (grid%mask2dT(i,j)==0.) sshy(i,j) = 0.0
  enddo ; enddo

  ! rotate ssh gradients from local coordinates to true zonal/meridional (inverse transformation)
  n = 0
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    n = n+1
    o2x(ind%o2x_So_dhdx, n) = grid%cos_rot(i,j) * sshx(i,j) + grid%sin_rot(i,j) * sshy(i,j)
    o2x(ind%o2x_So_dhdy, n) = grid%cos_rot(i,j) * sshy(i,j) - grid%sin_rot(i,j) * sshx(i,j)
  enddo ; enddo

end procedure ocn_export
end submodule ocn_cap_methods_s
