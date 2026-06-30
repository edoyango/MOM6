submodule (MOM_int_tide_input) MOM_int_tide_input_s
#include <MOM_memory.h>
  implicit none
contains
module procedure set_int_tide_input
  real, dimension(SZI_(G),SZJ_(G)) :: &
    N2_bot        ! The bottom squared buoyancy frequency [T-2 ~> s-2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Rho_bot, &   ! The average near-bottom density or the Boussinesq reference density [R ~> kg m-3].
    h_bot        ! Bottom boundary layer thickness [H ~> m or kg m-2].
  integer, dimension(SZI_(G),SZJ_(G)) ::  k_bot ! Bottom boundary layer top layer index.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    T_f, S_f      ! The temperature and salinity in [C ~> degC] and [S ~> ppt] with the values in
  logical :: use_EOS    ! If true, density is calculated from T & S using an
  logical :: avg_enabled  ! for testing internal tides (BDM)
  type(time_type) :: time_end        !< For use in testing internal tides (BDM)
  real :: HZ2_T3_to_W_m2  ! unit conversion factor for TKE from internal units
  real :: W_m2_to_HZ2_T3  ! unit conversion factor for TKE from mks to internal
  integer :: i, j, is, ie, js, je, nz, isd, ied, jsd, jed
  integer :: i_global, j_global
  integer :: fr
  HZ2_T3_to_W_m2 = GV%H_to_kg_m2*(US%Z_to_m**2)*(US%s_to_T**3)
  W_m2_to_HZ2_T3 = GV%kg_m2_to_H*(US%m_to_Z**2)*(US%T_to_s**3)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%initialized) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  use_EOS = associated(tv%eqn_of_state)

  ! Smooth the properties through massless layers.
  if (use_EOS) then
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_fill*dt, T_f, S_f, G, GV, US, larger_h_denom=.true.)
  endif

  call find_N2_bottom(G, GV, US, tv, fluxes, h, T_f, S_f, itide%h2, N2_bot, Rho_bot, h_bot, k_bot)

  avg_enabled = query_averaging_enabled(CS%diag, time_end=time_end)

  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    !$OMP parallel do default(shared)
    do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
      itide%Nb(i,j) = G%mask2dT(i,j) * sqrt(N2_bot(i,j))
      CS%TKE_itidal_input(i,j,fr) = min(GV%RZ_to_H*GV%Z_to_H*CS%TKE_itidal_coef(i,j,fr)*itide%Nb(i,j), &
                                        CS%TKE_itide_maxi)
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
      itide%Nb(i,j) = G%mask2dT(i,j) * sqrt(N2_bot(i,j))
      itide%Rho_bot(i,j) = G%mask2dT(i,j) * Rho_bot(i,j)
      CS%TKE_itidal_input(i,j,fr) = min((GV%RZ_to_H*GV%RZ_to_H*Rho_bot(i,j))*CS%TKE_itidal_coef(i,j,fr)*itide%Nb(i,j), &
                                        CS%TKE_itide_maxi)
    enddo ; enddo ; enddo
  endif

  if (CS%int_tide_source_test) then
    CS%TKE_itidal_input(:,:,:) = 0.0
    if (time_end <= CS%time_max_source) then
      if (CS%int_tide_use_glob_ij) then
        do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
          i_global = i + G%idg_offset
          j_global = j + G%jdg_offset
          if ((i_global == CS%int_tide_source_i) .and. (j_global == CS%int_tide_source_j)) then
            CS%TKE_itidal_input(i,j,fr) = 1.0*W_m2_to_HZ2_T3
          endif
        enddo ; enddo ; enddo
      else
        do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
          ! Input  an arbitrary energy point source.id_
          if (((G%geoLonCu(I-1,j)-CS%int_tide_source_x) * (G%geoLonBu(I,j)-CS%int_tide_source_x) <= 0.0) .and. &
              ((G%geoLatCv(i,J-1)-CS%int_tide_source_y) * (G%geoLatCv(i,j)-CS%int_tide_source_y) <= 0.0)) then
            CS%TKE_itidal_input(i,j,fr) = 1.0*W_m2_to_HZ2_T3
          endif
        enddo ; enddo ; enddo
      endif
    endif
  endif

  if (CS%debug) then
    call hchksum(N2_bot, "N2_bot", G%HI, haloshift=0, unscale=US%s_to_T**2)
    call hchksum(CS%TKE_itidal_input,"TKE_itidal_input", G%HI, haloshift=0, &
                 unscale=HZ2_T3_to_W_m2)
  endif

  call enable_averages(dt, time_end, CS%diag)

  do fr=1,CS%nFreq
    if (CS%id_TKE_itidal_itide(fr) > 0) call post_data(CS%id_TKE_itidal_itide(fr), &
                                                       CS%TKE_itidal_input(isd:ied,jsd:jed,fr), CS%diag)
  enddo
  if (CS%id_Nb > 0) call post_data(CS%id_Nb, itide%Nb, CS%diag)
  if (CS%id_N2_bot > 0 ) call post_data(CS%id_N2_bot, N2_bot, CS%diag)

  call disable_averaging(CS%diag)

end procedure set_int_tide_input
module procedure find_N2_bottom
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    pres, &       ! The pressure at each interface [R L2 T-2 ~> Pa].
    dRho_int      ! The unfiltered density differences across interfaces [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)) :: dz ! Layer thicknesses in depth units [Z ~> m]
  real, dimension(SZI_(G)) :: &
    Temp_int, &   ! The temperature at each interface [C ~> degC]
    Salin_int, &  ! The salinity at each interface [S ~> ppt]
    drho_bot, &   ! The density difference at the bottom of a layer [R ~> kg m-3]
    h_amp, &      ! The amplitude of topographic roughness [Z ~> m].
    hb, &         ! The thickness of the water column below the midpoint of a layer [H ~> m or kg m-2]
    z_from_bot, & ! The distance of a layer center from the bottom [Z ~> m]
    dRho_dT, &    ! The partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS       ! The partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: dz_int  ! The vertical extent of water associated with an interface [Z ~> m]
  real :: G_Rho0  ! The gravitational acceleration, sometimes divided by the Boussinesq
  logical :: do_i(SZI_(G)), do_any
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  G_Rho0 = GV%g_Earth_Z_T2 / GV%H_to_RZ
  EOSdom(:) = EOS_domain(G%HI)

  ! Find the (limited) density jump across each interface.
  do i=is,ie
    dRho_int(i,1) = 0.0 ; dRho_int(i,nz+1) = 0.0
  enddo

  !$OMP parallel do default(none) shared(is,ie,js,je,nz,tv,fluxes,G,GV,US,h,T_f,S_f, &
  !$OMP                                  h2,N2_bot,Rho_bot,h_bot,k_bot,G_Rho0,EOSdom) &
  !$OMP                          private(pres,Temp_Int,Salin_Int,dRho_dT,dRho_dS, &
  !$OMP                                  dz,hb,dRho_bot,z_from_bot,do_i,h_amp,do_any,dz_int) &
  !$OMP                     firstprivate(dRho_int)
  do j=js,je

    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dz, j, G, GV)

    if (associated(tv%eqn_of_state)) then
      if (associated(fluxes%p_surf)) then
        do i=is,ie ; pres(i,1) = fluxes%p_surf(i,j) ; enddo
      else
        do i=is,ie ; pres(i,1) = 0.0 ; enddo
      endif
      do K=2,nz
        do i=is,ie
          pres(i,K) = pres(i,K-1) + (GV%g_Earth*GV%H_to_RZ)*h(i,j,k-1)
          Temp_Int(i) = 0.5 * (T_f(i,j,k) + T_f(i,j,k-1))
          Salin_Int(i) = 0.5 * (S_f(i,j,k) + S_f(i,j,k-1))
        enddo
        call calculate_density_derivs(Temp_int, Salin_int, pres(:,K), dRho_dT(:), dRho_dS(:), &
                                      tv%eqn_of_state, EOSdom)
        do i=is,ie
          dRho_int(i,K) = max(dRho_dT(i)*(T_f(i,j,k) - T_f(i,j,k-1)) + &
                              dRho_dS(i)*(S_f(i,j,k) - S_f(i,j,k-1)), 0.0)
        enddo
      enddo
    else
      do K=2,nz ; do i=is,ie
        dRho_int(i,K) = (GV%Rlay(k) - GV%Rlay(k-1))
      enddo ; enddo
    endif

    ! Find the bottom boundary layer stratification.
    do i=is,ie
      hb(i) = 0.0 ; dRho_bot(i) = 0.0
      z_from_bot(i) = 0.5*dz(i,nz)
      do_i(i) = (G%mask2dT(i,j) > 0.0)
      h_amp(i) = sqrt(h2(i,j))
    enddo

    do k=nz,2,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        dz_int = 0.5*(dz(i,k) + dz(i,k-1))
        z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

        hb(i) = hb(i) + 0.5*(h(i,j,k) + h(i,j,k-1))
        dRho_bot(i) = dRho_bot(i) + dRho_int(i,K)

        if (z_from_bot(i) > h_amp(i)) then
          if (k>2) then
            ! Always include at least one full layer.
            hb(i) = hb(i) + 0.5*(h(i,j,k-1) + h(i,j,k-2))
            dRho_bot(i) = dRho_bot(i) + dRho_int(i,K-1)
          endif
          do_i(i) = .false.
        else
          do_any = .true.
        endif
      endif ; enddo
      if (.not.do_any) exit
    enddo

    do i=is,ie
      if (hb(i) > 0.0) then
        N2_bot(i,j) = (G_Rho0 * dRho_bot(i)) / hb(i)
      else ;  N2_bot(i,j) = 0.0 ; endif
    enddo

    if (GV%Boussinesq .or. GV%semi_Boussinesq) then
      do i=is,ie
        rho_bot(i,j) = GV%Rho0
      enddo
    else
      ! Average the density over the envelope of the topography.
      call find_rho_bottom(G, GV, US, tv, h, dz, pres, h_amp, j, Rho_bot(:,j), h_bot(:,j), k_bot(:,j))
    endif
  enddo

end procedure find_N2_bottom
module procedure get_input_TKE
  integer :: i,j,fr
  do fr=1,nFreq ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    TKE_itidal_input(i,j,fr) = CS%TKE_itidal_input(i,j,fr)
  enddo ; enddo ; enddo

end procedure get_input_TKE
module procedure get_barotropic_tidal_vel
  integer :: i,j,fr
  do fr=1,nFreq ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    vel_btTide(i,j,fr) = CS%tideamp(i,j,fr)
  enddo ; enddo ; enddo

end procedure get_barotropic_tidal_vel
module procedure int_tide_input_init
  logical :: read_tideamp
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_int_tide_input"  ! This module's name.
  character(len=200) :: filename, tideamp_file, h2_file ! Input file names or paths
  character(len=80)  :: tideamp_var, rough_var ! Input file variable names
  character(len=80)  :: var_name
  character(len=200) :: var_descript
  character(len=200) :: tidefile_varnames
  real :: mask_itidal        ! A multiplicative land mask, 0 or 1 [nondim]
  real :: max_frac_rough     ! The fraction relating the maximum topographic roughness
  real :: utide              ! constant tidal amplitude [L T-1 ~> m s-1] to be used if
  real :: kappa_h2_factor    ! factor for the product of wavenumber * rms sgs height [nondim].
  real :: kappa_itides       ! topographic wavenumber and non-dimensional scaling [L-1 ~> m-1]
  real :: min_zbot_itides    ! Minimum ocean depth for internal tide conversion [Z ~> m].
  real :: HZ2_T3_to_W_m2     ! unit conversion factor for TKE from internal units
  real :: W_m2_to_HZ2_T3     ! unit conversion factor for TKE from mks to internal
  integer :: tlen_days       !< Time interval from start for adding wave source
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: num_freq, fr
  if (associated(CS)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "control structure.")
    return
  endif
  if (associated(itide)) then
    call MOM_error(WARNING, "int_tide_input_init called with an associated "// &
                            "internal tide input type.")
    return
  endif
  allocate(CS)
  allocate(itide)

  CS%initialized = .true.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag

  HZ2_T3_to_W_m2 = GV%H_to_kg_m2*(US%Z_to_m**2)*(US%s_to_T**3)
  W_m2_to_HZ2_T3 = GV%kg_m2_to_H*(US%m_to_Z**2)*(US%T_to_s**3)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "MIN_ZBOT_ITIDES", min_zbot_itides, &
               "Turn off internal tidal dissipation when the total "//&
               "ocean depth is less than this value.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_fill, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T)

  call get_param(param_file, mdl, "UTIDE", utide, &
               "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
               units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  call read_param(param_file, "INTERNAL_TIDE_FREQS", num_freq)
  CS%nFreq= num_freq

  allocate(itide%Nb(isd:ied,jsd:jed), source=0.0)
  allocate(itide%Rho_bot(isd:ied,jsd:jed), source=0.0)
  allocate(itide%h2(isd:ied,jsd:jed), source=0.0)
  allocate(CS%TKE_itidal_input(isd:ied,jsd:jed,num_freq), source=0.0)
  allocate(CS%tideamp(isd:ied,jsd:jed,num_freq), source=utide)
  allocate(CS%TKE_itidal_coef(isd:ied,jsd:jed, num_freq), source=0.0)

  call get_param(param_file, mdl, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. "//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0), scale=US%L_to_m)

  call get_param(param_file, mdl, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with "//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
  call get_param(param_file, mdl, "TKE_ITIDE_MAX", CS%TKE_itide_maxi, &
               "The maximum internal tide energy source available to mix "//&
               "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
               units="W m-2", default=1.0e3, scale=W_m2_to_HZ2_T3)

  call get_param(param_file, mdl, "READ_TIDEAMP", read_tideamp, &
               "If true, read a file (given by TIDEAMP_FILE) containing "//&
               "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (read_tideamp) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", tideamp_file, &
               "The path to the file containing the spatially varying "//&
               "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
    filename = trim(CS%inputdir) // trim(tideamp_file)
    call log_param(param_file, mdl, "INPUTDIR/TIDEAMP_FILE", filename)

    call read_param(param_file, "INTTIDE_AMP_VARNAMES", tidefile_varnames)
    do fr=1,num_freq
      tideamp_var = extractWord(tidefile_varnames,fr)
      call MOM_read_data(filename, tideamp_var, CS%tideamp(:,:,fr), G%domain, scale=US%m_s_to_L_T)
    enddo

  endif

  call get_param(param_file, mdl, "H2_FILE", h2_file, &
               "The path to the file containing the sub-grid-scale "//&
               "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
               fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
  call get_param(param_file, mdl, "ROUGHNESS_VARNAME", rough_var, &
                 "The name in the input file of the squared sub-grid-scale "//&
                 "topographic roughness amplitude variable.", default="h2")
  call MOM_read_data(filename, rough_var, itide%h2, G%domain, scale=US%m_to_Z**2)

  call get_param(param_file, mdl, "FRACTIONAL_ROUGHNESS_MAX", max_frac_rough, &
                 "The maximum topographic roughness amplitude as a fraction of the mean depth, "//&
                 "or a negative value for no limitations on roughness.", &
                 units="nondim", default=0.1)

  ! The following parameters are used in testing the internal tide code.
  call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_TEST", CS%int_tide_source_test, &
                 "If true, apply an arbitrary generation site for internal tide testing", &
                 default=.false.)
  if (CS%int_tide_source_test)then
    call get_param(param_file, mdl, "INTERNAL_TIDE_USE_GLOB_IJ", CS%int_tide_use_glob_ij, &
                 "Use global IJ for internal tide generation source test", default=.false.)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_X", CS%int_tide_source_x, &
                 "X Location of generation site for internal tide", &
                 units=G%x_ax_unit_short, default=1.0, do_not_log=CS%int_tide_use_glob_ij)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_Y", CS%int_tide_source_y, &
                 "Y Location of generation site for internal tide", &
                 units=G%y_ax_unit_short, default=1.0, do_not_log=CS%int_tide_use_glob_ij)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_I", CS%int_tide_source_i, &
                 "I Location of generation site for internal tide", default=0, &
                 do_not_log=.not.CS%int_tide_use_glob_ij)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_J", CS%int_tide_source_j, &
                 "J Location of generation site for internal tide", default=0, &
                 do_not_log=.not.CS%int_tide_use_glob_ij)
    call get_param(param_file, mdl, "INTERNAL_TIDE_SOURCE_TLEN_DAYS", tlen_days, &
                 "Time interval from start of experiment for adding wave source", &
                 units="days", default=0)
    CS%time_max_source = Time + set_time(0, days=tlen_days)

    if ((CS%int_tide_use_glob_ij) .and. ((CS%int_tide_source_x /= 1.) .or. (CS%int_tide_source_y /= 1.))) then
      call MOM_error(FATAL, "MOM_internal_tide_input: "//&
                     "Internal tide source set to use (i,j) indices hence (x,y) geographical coords are meaningless.")
    endif
    if ((.not.CS%int_tide_use_glob_ij) .and. ((CS%int_tide_source_i /= 0) .or. (CS%int_tide_source_j /= 0))) then
      call MOM_error(FATAL, "MOM_internal_tide_input: "//&
                     "Internal tide source set to use (x,y) geographical coords hence (i,j) indices are meaningless.")
    endif
  endif

  do fr=1,num_freq ; do j=js,je ; do i=is,ie
    mask_itidal = 1.0
    if (G%meanSL(i,j) + G%bathyT(i,j) < min_zbot_itides) mask_itidal = 0.0

    CS%tideamp(i,j,fr) = CS%tideamp(i,j,fr) * mask_itidal * G%mask2dT(i,j)

    ! Restrict rms topo to a fraction (often 10 percent) of the column depth.
    if (max_frac_rough >= 0.0) &
      itide%h2(i,j) = min((max_frac_rough * max(G%meanSL(i,j) + G%bathyT(i,j), 0.0))**2, itide%h2(i,j))

    ! Compute the fixed part of internal tidal forcing; units are [R Z4 H-1 T-2 ~> J m-2 or J m kg-1] here.
    CS%TKE_itidal_coef(i,j,fr) = 0.5*US%L_to_Z*kappa_h2_factor * GV%H_to_RZ * &
         kappa_itides * itide%h2(i,j) * CS%tideamp(i,j,fr)**2
  enddo ; enddo ; enddo


  allocate( CS%id_TKE_itidal_itide(num_freq), source=-1)

  do fr=1,num_freq
    write(var_name, '("TKE_itidal_itide_freq",i1)') fr
    write(var_descript, '("Internal Tide Driven Turbulent Kinetic Energy in frequency ",i1)') fr

    CS%id_TKE_itidal_itide(fr) = register_diag_field('ocean_model',var_name,diag%axesT1,Time, &
                                                     var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
  enddo

  CS%id_Nb = register_diag_field('ocean_model','Nb_itide',diag%axesT1,Time, &
       'Bottom Buoyancy Frequency', 's-1', conversion=US%s_to_T)

  CS%id_N2_bot = register_diag_field('ocean_model','N2_b_itide',diag%axesT1,Time, &
       'Bottom Buoyancy frequency squared', 's-2', conversion=US%s_to_T**2)

end procedure int_tide_input_init
module procedure int_tide_input_end
  if (associated(CS)) deallocate(CS)

end procedure int_tide_input_end
end submodule MOM_int_tide_input_s
