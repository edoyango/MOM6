submodule (MOM_surface_forcing) MOM_surface_forcing_s
#include <MOM_memory.h>
  implicit none
contains
module procedure set_forcing
  real :: dt                     ! length of time over which fluxes applied [T ~> s]
  type(time_type) :: day_center  ! central time of the fluxes.
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call cpu_clock_begin(id_clock_forcing)
  call callTree_enter("set_forcing, MOM_surface_forcing.F90")

  day_center = day_start + day_interval/2
  dt = time_to_real(day_interval, scale=US%s_to_T)

  if (CS%first_call_set_forcing) then
    ! Allocate memory for the mechanical and thermodynamic forcing fields.
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.not.CS%nonBous, press=.true., tau_mag=CS%nonBous)

    call allocate_forcing_type(G, fluxes, ustar=.not.CS%nonBous, marbl=CS%use_marbl_tracers, tau_mag=CS%nonBous, &
                               fix_accum_bug=.not.CS%ustar_gustless_bug)
    if (trim(CS%buoy_config) /= "NONE") then
      if ( CS%use_temperature ) then
        call allocate_forcing_type(G, fluxes, water=.true., heat=.true., press=.true.)
        if (CS%restorebuoy) then
          call safe_alloc_ptr(CS%T_Restore,isd, ied, jsd, jed)
          call safe_alloc_ptr(fluxes%heat_added, isd, ied, jsd, jed)
          call safe_alloc_ptr(CS%S_Restore, isd, ied, jsd, jed)
        endif
      else ! CS%use_temperature false.
        call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)

        if (CS%restorebuoy) call safe_alloc_ptr(CS%Dens_Restore, isd, ied, jsd, jed)
      endif  ! endif for  CS%use_temperature
    endif
  endif

  ! calls to various wind options
  if (CS%variable_winds .or. CS%first_call_set_forcing) then

    if (trim(CS%wind_config) == "file") then
      call wind_forcing_from_file(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "data_override") then
      call wind_forcing_by_data_override(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "2gyre") then
      call wind_forcing_2gyre(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "1gyre") then
      call wind_forcing_1gyre(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "gyres") then
      call wind_forcing_gyres(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "zero") then
      call wind_forcing_const(sfc_state, forces, 0., 0., day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "const") then
      call wind_forcing_const(sfc_state, forces, CS%tau_x0, CS%tau_y0, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "Neverworld" .or. trim(CS%wind_config) == "Neverland") then
      call Neverworld_wind_forcing(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "scurves") then
      call scurve_wind_forcing(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "ideal_hurr") then
      call idealized_hurricane_wind_forcing(sfc_state, forces, day_center, G, US, CS%idealized_hurricane_CSp)
    elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
      call MOM_error(FATAL, "MOM_surface_forcing (set_forcing): "//&
                            'WIND_CONFIG = "SCM_ideal_hurr" is a depricated option.')
    elseif (trim(CS%wind_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_wind_forcing(sfc_state, forces, day_center, G, US, CS%SCM_CVmix_tests_CSp)
    elseif (trim(CS%wind_config) == "USER") then
      call USER_wind_forcing(sfc_state, forces, day_center, G, US, CS%user_forcing_CSp)
    elseif (CS%variable_winds .and. .not.CS%first_call_set_forcing) then
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Variable winds defined with no wind config")
    else
      call MOM_error(FATAL, &
       "MOM_surface_forcing:Unrecognized wind config "//trim(CS%wind_config))
    endif
  endif

  ! calls to various buoyancy forcing options
  if (CS%restorebuoy .and. .not.CS%variable_buoyforce) then
    call MOM_error(FATAL, "With RESTOREBUOY = True, VARIABLE_BUOYFORCE = True should be used. "//&
                          "Otherwise, this can lead to diverging solutions when a simulation "//&
                          "is continued using a restart file.")
  endif

  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    if (trim(CS%buoy_config) == "file") then
      call buoyancy_forcing_from_files(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "data_override") then
      call buoyancy_forcing_from_data_override(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "zero") then
      call buoyancy_forcing_zero(sfc_state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "const") then
      call buoyancy_forcing_const(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "linear") then
      call buoyancy_forcing_linear(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "MESO") then
      call MESO_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%MESO_forcing_CSp)
    elseif (trim(CS%buoy_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_buoyancy_forcing(sfc_state, fluxes, day_center, G, US, CS%SCM_CVmix_tests_CSp)
    elseif (trim(CS%buoy_config) == "USER") then
      call USER_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%user_forcing_CSp)
    elseif (trim(CS%buoy_config) == "BFB") then
      call BFB_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%BFB_forcing_CSp)
    elseif (trim(CS%buoy_config) == "dumbbell") then
      call dumbbell_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%dumbbell_forcing_CSp)
    elseif (trim(CS%buoy_config) == "NONE") then
      call MOM_mesg("MOM_surface_forcing: buoyancy forcing has been set to omitted.")
    elseif (CS%variable_buoyforce .and. .not.CS%first_call_set_forcing) then
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Variable buoy defined with no buoy config.")
    else
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Unrecognized buoy config "//trim(CS%buoy_config))
    endif
  endif

  if (CS%use_marbl_tracers) then
    call MARBL_forcing_from_data_override(fluxes, day_center, G, US, CS)
  endif

  if (associated(CS%tracer_flow_CSp)) then
    call call_tracer_set_forcing(sfc_state, fluxes, day_start, day_interval, G, US, CS%Rho0, &
                                 CS%tracer_flow_CSp)
  endif

  ! Allow for user-written code to alter the fluxes after all the above
  call user_alter_forcing(sfc_state, fluxes, day_center, G, CS%urf_CS)

  ! Fields that exist in both the forcing and mech_forcing types must be copied.
  if (CS%variable_winds .or. CS%first_call_set_forcing) then
    call copy_common_forcing_fields(forces, fluxes, G)
    call set_derived_forcing_fields(forces, fluxes, G, US, CS%Rho0)
  endif

  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    call set_net_mass_forcing(fluxes, forces, G, US)
  endif

  CS%first_call_set_forcing = .false.

  call cpu_clock_end(id_clock_forcing)
  call callTree_leave("set_forcing")

end procedure set_forcing
module procedure wind_forcing_const
  real :: mag_tau  ! Magnitude of the wind stress [R Z2 T-2 ~> Pa]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  call callTree_enter("wind_forcing_const, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  mag_tau = US%L_to_Z * sqrt( tau_x0**2 + tau_y0**2)

  ! Set the steady surface wind stresses, in units of [R L Z T-2 ~> Pa].
  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = tau_x0
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = tau_y0
  enddo ; enddo

  if (CS%read_gust_2d) then
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( ( mag_tau + CS%gust(i,j) ) / CS%Rho0 )
    enddo ; enddo ; endif
    if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = mag_tau + CS%gust(i,j)
    enddo ; enddo ; endif
  else
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( ( mag_tau + CS%gust_const ) / CS%Rho0 )
    enddo ; enddo ; endif
    if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = mag_tau + CS%gust_const
    enddo ; enddo ; endif
  endif

  call callTree_leave("wind_forcing_const")
end procedure wind_forcing_const
module procedure wind_forcing_2gyre
  real :: PI            ! A common irrational number, 3.1415926535... [nondim]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  call callTree_enter("wind_forcing_2gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  PI = 4.0*atan(1.0)

  ! Set the steady surface wind stresses, in units of [R L Z T-2 ~> Pa].
  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = CS%taux_mag * (1.0 - cos(2.0*PI*(G%geoLatCu(I,j)-CS%South_lat) / CS%len_lat))
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  if (associated(forces%ustar)) call stresses_to_ustar(forces, G, US, CS)

  call callTree_leave("wind_forcing_2gyre")
end procedure wind_forcing_2gyre
module procedure wind_forcing_1gyre
  real :: PI            ! A common irrational number, 3.1415926535... [nondim]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  call callTree_enter("wind_forcing_1gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  PI = 4.0*atan(1.0)

  ! Set the steady surface wind stresses, in units of [R Z L T-2 ~> Pa].
  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = CS%taux_mag * cos(PI*(G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat)
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  if (associated(forces%ustar)) call stresses_to_ustar(forces, G, US, CS)

  call callTree_leave("wind_forcing_1gyre")
end procedure wind_forcing_1gyre
module procedure wind_forcing_gyres
  real :: PI            ! A common irrational number, 3.1415926535... [nondim]
  real :: y             ! The latitude relative to the south normalized by the domain extent [nondim]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  call callTree_enter("wind_forcing_gyres, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  PI = 4.0*atan(1.0)

  ! steady surface wind stresses [R L Z T-2 ~> Pa]
  do j=js-1,je+1 ; do I=is-1,Ieq
    y = (G%geoLatCu(I,j)-CS%South_lat) / CS%len_lat
    forces%taux(I,j) = CS%gyres_taux_const + &
                       (   CS%gyres_taux_sin_amp*sin(CS%gyres_taux_n_pis*PI*y) &
                         + CS%gyres_taux_cos_amp*cos(CS%gyres_taux_n_pis*PI*y) )
  enddo ; enddo

  do J=js-1,Jeq ; do i=is-1,ie+1
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  ! set the friction velocity
  if (CS%answer_date < 20190101) then
    if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = CS%gust_const + US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                                                  ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
    enddo ; enddo ; endif
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( (CS%gust_const/CS%Rho0) + &
              US%L_to_Z * sqrt(0.5*((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2) + &
                                    (forces%taux(I-1,j)**2) + (forces%taux(I,j)**2)))/CS%Rho0 )
    enddo ; enddo ; endif
  else
    call stresses_to_ustar(forces, G, US, CS)
  endif

  call callTree_leave("wind_forcing_gyres")
end procedure wind_forcing_gyres
module procedure Neverworld_wind_forcing
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: PI            ! A common irrational number, 3.1415926535... [nondim]
  real :: y             ! The latitude relative to the south normalized by the domain extent [nondim]
  real :: tau_max       ! The magnitude of the wind stress [R Z L T-2 ~> Pa]
  real :: off           ! An offset in the relative latitude [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true.)

  !  Set the surface wind stresses, in units of [R Z L T-2 ~> Pa].  A positive taux
  !  accelerates the ocean to the (pseudo-)east.

  !  The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  PI = 4.0*atan(1.0)

  forces%taux(:,:) = 0.0
  tau_max = CS%taux_mag
  off = 0.02
  do j=js,je ; do I=is-1,Ieq
    y = (G%geoLatT(i,j)-G%south_lat)/G%len_lat

    if (y <= 0.29) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max * ( (1/0.29)*y - ( 1/(2*PI) )*sin( (2*PI*y) / 0.29 ) )
    endif
    if ((y > 0.29) .and. (y <= (0.8-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *(0.35+0.65*cos(PI*(y-0.29)/(0.51-off))  )
    endif
    if ((y > (0.8-off)) .and. (y <= (1-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *( 1.5*( (y-1+off) - (0.1/PI)*sin(10.0*PI*(y-0.8+off)) ) )
    endif
    forces%taux(I,j) = G%mask2dCu(I,j) * forces%taux(I,j)
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0
  enddo ; enddo

  ! Set the surface friction velocity, in units of [Z T-1 ~> m s-1].  ustar is always positive.
  if (associated(forces%ustar)) call stresses_to_ustar(forces, G, US, CS)

end procedure Neverworld_wind_forcing
module procedure scurve_wind_forcing
  integer :: i, j, kseg
  real :: y_curve       ! The latitude relative to the southern end of a curve segment [degreesN]
  real :: L_curve       ! The latitudinal extent of a curve segment [degreesN]
  call allocate_mech_forcing(G, forces, stress=.true.)

  kseg = 1
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    ! Find segment k s.t. ydata(k)<= G%geoLatCu(I,j) < ydata(k+1)
    do while (G%geoLatCu(I,j) >= CS%scurves_ydata(kseg+1) .and. kseg<6) ! Should this be kseg<19?
      kseg = kseg+1
    enddo
    do while (G%geoLatCu(I,j) < CS%scurves_ydata(kseg) .and. kseg>1)
      kseg = kseg-1
    enddo

    y_curve = G%geoLatCu(I,j) - CS%scurves_ydata(kseg)
    L_curve = CS%scurves_ydata(kseg+1) - CS%scurves_ydata(kseg)
    forces%taux(I,j) = CS%scurves_taux(kseg) +  &
                       (CS%scurves_taux(kseg+1) - CS%scurves_taux(kseg)) * scurve(y_curve, L_curve)
    forces%taux(I,j) = G%mask2dCu(I,j) * forces%taux(I,j)
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0
  enddo ; enddo

  ! Set the surface friction velocity, in units of [Z T-1 ~> m s-1].  ustar is always positive.
  if (associated(forces%ustar)) call stresses_to_ustar(forces, G, US, CS)

end procedure scurve_wind_forcing
module procedure scurve
  real :: s  ! The evaluated function value [nondim]
  s = x/L
  scurve = (3. - 2.*s) * (s*s)
end procedure scurve
module procedure wind_forcing_from_file
  character(len=200) :: filename  ! The name of the input file.
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal wind stresses at h-points [R L Z T-2 ~> Pa]
  real    :: temp_y(SZI_(G),SZJ_(G)) ! Pseudo-meridional wind stresses at h-points [R L Z T-2 ~> Pa]
  real    :: ustar_loc(SZI_(G),SZJ_(G)) ! The local value of ustar [Z T-1 ~> m s-1]
  real    :: tau_mag    ! The magnitude of the wind stress including any contributions from
  integer :: time_lev                ! The time level that is used for a field.
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  logical :: read_Ustar
  call callTree_enter("wind_forcing_from_file, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  time_lev = get_file_time_level(day, CS%wind_nlev, CS%wind_days_per_rec)

  if (time_lev /= CS%wind_last_lev) then
    filename = trim(CS%wind_file)
    read_Ustar = (len_trim(CS%ustar_var) > 0)
!    if (is_root_pe()) &
!      write(*,'("Wind_forcing Reading time level ",I," last was ",I,".")')&
!           time_lev-1,CS%wind_last_lev-1
    select case ( uppercase(CS%wind_stagger(1:1)) )
    case ("A")
      temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
      call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                           temp_x(:,:), temp_y(:,:), G%Domain, stagger=AGRID, &
                           timelevel=time_lev, scale=US%Pa_to_RLZ_T2)

      call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
      do j=js,je ; do I=is-1,Ieq
        forces%taux(I,j) = 0.5 * CS%wind_scale * (temp_x(i,j) + temp_x(i+1,j))
      enddo ; enddo
      do J=js-1,Jeq ; do i=is,ie
        forces%tauy(i,J) = 0.5 * CS%wind_scale * (temp_y(i,j) + temp_y(i,j+1))
      enddo ; enddo

      if (.not.read_Ustar) then
        if (CS%read_gust_2d) then
          if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
            forces%tau_mag(i,j) = CS%gust(i,j) + US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2))
          enddo ; enddo ; endif
          if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
            tau_mag = CS%gust(i,j) + US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2))
            forces%ustar(i,j) = sqrt(tau_mag / CS%Rho0)
          enddo ; enddo ; endif
        else
          if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
            forces%tau_mag(i,j) = CS%gust_const + US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2))
          enddo ; enddo ; endif
          if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
            forces%ustar(i,j) = sqrt( CS%gust_const/CS%Rho0 + &
                    US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2)) / CS%Rho0 )
          enddo ; enddo ; endif
        endif
      endif
    case ("C")
      if (G%symmetric) then
        if (.not.associated(G%Domain_aux)) call MOM_error(FATAL, &
          " wind_forcing_from_file with C-grid input and symmetric memory "//&
          " called with a non-associated auxiliary domain in the grid type.")
        !   Read the data as though symmetric memory were not being used, and
        ! then translate it appropriately.
        temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
        call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                             temp_x(:,:), temp_y(:,:), &
                             G%Domain_aux, stagger=CGRID_NE, timelevel=time_lev, &
                             scale=US%Pa_to_RLZ_T2)
        do j=js,je ; do i=is,ie
          forces%taux(I,j) = CS%wind_scale * temp_x(I,j)
          forces%tauy(i,J) = CS%wind_scale * temp_y(i,J)
        enddo ; enddo
        call fill_symmetric_edges(forces%taux, forces%tauy, G%Domain, stagger=CGRID_NE)
      else
        call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                             forces%taux(:,:), forces%tauy(:,:), &
                             G%Domain, stagger=CGRID_NE, timelevel=time_lev, &
                             scale=US%Pa_to_RLZ_T2)

        if (CS%wind_scale /= 1.0) then
          do j=js,je ; do I=Isq,Ieq
            forces%taux(I,j) = CS%wind_scale * forces%taux(I,j)
          enddo ; enddo
          do J=Jsq,Jeq ; do i=is,ie
            forces%tauy(i,J) = CS%wind_scale * forces%tauy(i,J)
          enddo ; enddo
        endif
      endif

      call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)
      if (.not.read_Ustar) then
        if (CS%read_gust_2d) then
          if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
            forces%tau_mag(i,j) = CS%gust(i,j) + &
                    US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                          ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
          enddo ; enddo ; endif
          if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
            tau_mag = CS%gust(i,j) + &
                   US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                         ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
            forces%ustar(i,j) = sqrt( tau_mag / CS%Rho0 )
          enddo ; enddo ; endif
        else
          if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
            forces%tau_mag(i,j) = CS%gust_const + &
                  US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                        ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
          enddo ; enddo ; endif
          if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
             forces%ustar(i,j) = sqrt( CS%gust_const/CS%Rho0 + &
                  US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                        ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))/CS%Rho0 )
          enddo ; enddo ; endif
        endif
      endif
    case default
      call MOM_error(FATAL, "wind_forcing_from_file: Unrecognized stagger "//&
                      trim(CS%wind_stagger)//" is not 'A' or 'C'.")
    end select

    if (read_Ustar) then
      call MOM_read_data(filename, CS%Ustar_var, ustar_loc(:,:), &
                         G%Domain, timelevel=time_lev, scale=US%m_to_Z*US%T_to_s)
      if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
        forces%tau_mag(i,j) = CS%Rho0 * ustar_loc(i,j)**2
      enddo ; enddo ; endif
      if (associated(forces%ustar)) then ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        forces%ustar(i,j) = ustar_loc(i,j)
      enddo ; enddo ; endif
    endif

    CS%wind_last_lev = time_lev

  endif ! time_lev /= CS%wind_last_lev

  call callTree_leave("wind_forcing_from_file")
end procedure wind_forcing_from_file
module procedure wind_forcing_by_data_override
  real :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal wind stresses at h-points [R Z L T-2 ~> Pa].
  real :: temp_y(SZI_(G),SZJ_(G)) ! Pseudo-meridional wind stresses at h-points [R Z L T-2 ~> Pa].
  real :: ustar_prev(SZI_(G),SZJ_(G)) ! The pre-override value of ustar [Z T-1 ~> m s-1]
  real :: ustar_loc(SZI_(G),SZJ_(G)) ! The value of ustar, perhaps altered by data override [Z T-1 ~> m s-1]
  real :: tau_mag       ! The magnitude of the wind stress including any contributions from
  integer :: i, j
  call callTree_enter("wind_forcing_by_data_override, MOM_surface_forcing.F90")

  if (.not.CS%dataOverrideIsInitialized) then
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.not.CS%nonBous, press=.true., tau_mag=CS%nonBous)
    call data_override_init(G%Domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
  ! CS%wind_scale is ignored here because it is not set in this mode.
  call data_override(G%Domain, 'taux', temp_x, day, scale=US%Pa_to_RLZ_T2)
  call data_override(G%Domain, 'tauy', temp_y, day, scale=US%Pa_to_RLZ_T2)
  call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
  do j=G%jsc,G%jec ; do I=G%isc-1,G%IecB
    forces%taux(I,j) = 0.5 * (temp_x(i,j) + temp_x(i+1,j))
  enddo ; enddo
  do J=G%jsc-1,G%JecB ; do i=G%isc,G%iec
    forces%tauy(i,J) = 0.5 * (temp_y(i,j) + temp_y(i,j+1))
  enddo ; enddo

  if (CS%read_gust_2d) then
    call data_override(G%Domain, 'gust', CS%gust, day, scale=US%Pa_to_RLZ_T2*US%L_to_Z)
    if (associated(forces%tau_mag)) then ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      forces%tau_mag(i,j) = US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2)) + CS%gust(i,j)
    enddo ; enddo ; endif
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      tau_mag = US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2)) + CS%gust(i,j)
      ustar_loc(i,j) = sqrt( tau_mag / CS%Rho0 )
    enddo ; enddo
  else
    if (associated(forces%tau_mag)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        forces%tau_mag(i,j) = US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2)) + CS%gust_const
      ! ustar_loc(i,j) = sqrt( forces%tau_mag(i,j) / CS%Rho0 )
      enddo ; enddo
    endif
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ustar_loc(i,j) = sqrt(US%L_to_Z * sqrt((temp_x(i,j)**2) + (temp_y(i,j)**2))/CS%Rho0 + &
                            CS%gust_const/CS%Rho0)
    enddo ; enddo
  endif

  ! Give the data override the option to modify the newly calculated forces%ustar.
  ustar_prev(:,:) = ustar_loc(:,:)
  call data_override(G%Domain, 'ustar', ustar_loc, day, scale=US%m_to_Z*US%T_to_s)

  ! Only reset values where data override of ustar has occurred
  if (associated(forces%tau_mag)) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (ustar_prev(i,j) /= ustar_loc(i,j)) then
      forces%tau_mag(i,j) = CS%Rho0 * ustar_loc(i,j)**2
    endif ; enddo ; enddo
  endif

  if (associated(forces%ustar)) then ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    forces%ustar(i,j) = ustar_loc(i,j)
  enddo ; enddo ; endif

  call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)

  call callTree_leave("wind_forcing_by_data_override")
end procedure wind_forcing_by_data_override
module procedure stresses_to_ustar
  real :: I_rho         ! The inverse of the Boussinesq reference density [R-1 ~> m3 kg-1]
  real :: tau_mag       ! The magnitude of the wind stress including any contributions from
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  I_rho = 1.0 / CS%Rho0

  if (CS%read_gust_2d) then
    if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = CS%gust(i,j) + &
              US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                    ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
    enddo ; enddo ; endif
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      tau_mag = CS%gust(i,j) + &
              US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                    ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
      forces%ustar(i,j) = sqrt( tau_mag * I_rho )
    enddo ; enddo ; endif
  else
    if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
      forces%tau_mag(i,j) = CS%gust_const + &
              US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                    ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
    enddo ; enddo ; endif
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      tau_mag = CS%gust_const + &
              US%L_to_Z * sqrt(0.5*(((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)) + &
                                    ((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2))))
      forces%ustar(i,j) = sqrt( tau_mag * I_rho )
    enddo ; enddo ; endif
  endif

end procedure stresses_to_ustar
module procedure buoyancy_forcing_from_files
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp          ! A 2-d temporary work array in various units of [Q R Z T-1 ~> W m-2] or
  real :: rhoXcp ! reference density times heat capacity [Q R C-1 ~> J m-3 degC-1]
  logical :: fluxes_changed     ! True if any of the fluxes might have been altered
  integer :: time_lev           ! time level that for a field
  integer :: i, j, is, ie, js, je
  call callTree_enter("buoyancy_forcing_from_files, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je = G%jec

  if (CS%use_temperature) rhoXcp = CS%rho_restore * fluxes%C_p

  ! Read the buoyancy forcing file
  fluxes_changed = .false.

  ! longwave
  time_lev = get_file_time_level(day, CS%LW_nlev, CS%LW_days_per_rec)
  if (time_lev /= CS%LW_last_lev) then
    call MOM_read_data(CS%longwave_file, CS%LW_var, fluxes%lw(:,:), &
                       G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%longwaveup_file, "lwup_sfc", temp(:,:), G%Domain, &
                         timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
      do j=js,je ; do i=is,ie ; fluxes%LW(i,j) = fluxes%LW(i,j) - temp(i,j) ; enddo ; enddo
    endif
    CS%LW_last_lev = time_lev ; fluxes_changed = .true.
  endif

  ! evaporation
  if ( (CS%evap_nlev /= CS%LW_nlev) .or. (CS%evap_days_per_rec /= CS%LW_days_per_rec) ) &
    time_lev = get_file_time_level(day, CS%evap_nlev, CS%evap_days_per_rec)
  if (time_lev /= CS%evap_last_lev) then
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%evaporation_file, CS%evap_var, fluxes%evap(:,:), &
                     G%Domain, timelevel=time_lev, scale=-US%kg_m2s_to_RZ_T)
      do j=js,je ; do i=is,ie
        fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
        fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
      enddo ; enddo
    else
      call MOM_read_data(CS%evaporation_file, CS%evap_var, fluxes%evap(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    endif
    CS%evap_last_lev = time_lev ; fluxes_changed = .true.
  endif

  if ( (CS%latent_nlev /= CS%evap_nlev) .or. (CS%latent_days_per_rec /= CS%evap_days_per_rec) ) &
    time_lev = get_file_time_level(day, CS%latent_nlev, CS%latent_days_per_rec)
  if (time_lev /= CS%latent_last_lev) then
    if (.not.CS%archaic_OMIP_file) then
      call MOM_read_data(CS%latentheat_file, CS%latent_var, fluxes%latent(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
      do j=js,je ; do i=is,ie
        fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
      enddo ; enddo
    endif
    CS%latent_last_lev = time_lev ; fluxes_changed = .true.
  endif

  if ( (CS%sens_nlev /= CS%latent_nlev) .or. (CS%sens_days_per_rec /= CS%latent_days_per_rec) )  &
    time_lev = get_file_time_level(day, CS%sens_nlev, CS%sens_days_per_rec)
  if (time_lev /= CS%sens_last_lev) then
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%sensibleheat_file, CS%sens_var, fluxes%sens(:,:), &
                     G%Domain, timelevel=time_lev, scale=-US%W_m2_to_QRZ_T)
    else
      call MOM_read_data(CS%sensibleheat_file, CS%sens_var, fluxes%sens(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    endif
    CS%sens_last_lev = time_lev ; fluxes_changed = .true.
  endif

  if ( (CS%SW_nlev /= CS%sens_nlev) .or. (CS%SW_days_per_rec /= CS%sens_days_per_rec) )  &
    time_lev = get_file_time_level(day, CS%SW_nlev, CS%SW_days_per_rec)
  if (time_lev /= CS%SW_last_lev) then
    call MOM_read_data(CS%shortwave_file, CS%SW_var, fluxes%sw(:,:), G%Domain, &
                       timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%shortwaveup_file, "swup_sfc", temp(:,:), G%Domain, &
                         timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
      do j=js,je ; do i=is,ie
        fluxes%sw(i,j) = fluxes%sw(i,j) - temp(i,j)
      enddo ; enddo
    endif
    CS%SW_last_lev = time_lev ; fluxes_changed = .true.
  endif

  if ( (CS%precip_nlev /= CS%SW_nlev) .or. (CS%precip_days_per_rec /= CS%SW_days_per_rec) )  &
    time_lev = get_file_time_level(day, CS%precip_nlev, CS%precip_days_per_rec)
  if (time_lev /= CS%precip_last_lev) then
    call MOM_read_data(CS%snow_file, CS%snow_var, &
             fluxes%fprec(:,:), G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    call MOM_read_data(CS%rain_file, CS%rain_var, &
             fluxes%lprec(:,:), G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    if (CS%archaic_OMIP_file) then
      do j=js,je ; do i=is,ie
        fluxes%lprec(i,j) = fluxes%lprec(i,j) - fluxes%fprec(i,j)
      enddo ; enddo
    endif
    CS%precip_last_lev = time_lev ; fluxes_changed = .true.
  endif

  if ( (CS%runoff_nlev /= CS%precip_nlev) .or. (CS%runoff_days_per_rec /= CS%precip_days_per_rec) )  &
    time_lev = get_file_time_level(day, CS%runoff_nlev, CS%runoff_days_per_rec)
  if (time_lev /= CS%runoff_last_lev) then
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%runoff_file, CS%lrunoff_var, temp(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T*US%m_to_L**2)
      do j=js,je ; do i=is,ie
        fluxes%lrunoff(i,j) = temp(i,j)*G%IareaT(i,j)
      enddo ; enddo
      call MOM_read_data(CS%runoff_file, CS%frunoff_var, temp(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T*US%m_to_L**2)
      do j=js,je ; do i=is,ie
        fluxes%frunoff(i,j) = temp(i,j)*G%IareaT(i,j)
      enddo ; enddo
    else
      call MOM_read_data(CS%runoff_file, CS%lrunoff_var, fluxes%lrunoff(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
      call MOM_read_data(CS%runoff_file, CS%frunoff_var, fluxes%frunoff(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    endif
    CS%runoff_last_lev = time_lev ; fluxes_changed = .true.
  endif

  !  Read the SST and SSS fields for damping.
  if (CS%restorebuoy) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
    time_lev = get_file_time_level(day, CS%SST_nlev, CS%SST_days_per_rec)
    if (time_lev /= CS%SST_last_lev) then
      call MOM_read_data(CS%SSTrestore_file, CS%SST_restore_var, &
               CS%T_Restore(:,:), G%Domain, timelevel=time_lev, scale=US%degC_to_C)
      CS%SST_last_lev = time_lev
    endif

    if ( (CS%SSS_nlev /= CS%SST_nlev) .or. (CS%SSS_days_per_rec /= CS%SST_days_per_rec) )  &
      time_lev = get_file_time_level(day, CS%SSS_nlev, CS%SSS_days_per_rec)
    if (time_lev /= CS%SSS_last_lev) then
      call MOM_read_data(CS%salinityrestore_file, CS%SSS_restore_var, &
               CS%S_Restore(:,:), G%Domain, timelevel=time_lev, scale=US%ppt_to_S)
      CS%SSS_last_lev = time_lev
    endif
  endif

  if (fluxes_changed) then
    ! mask out land points and compute heat content of water fluxes
    ! assume liquid precipitation enters ocean at SST
    ! assume frozen precipitation enters ocean at 0degC
    ! assume liquid runoff enters ocean at SST
    ! assume solid runoff (calving) enters ocean at 0degC
    ! mass leaving the ocean has heat_content determined in MOM_diabatic_driver.F90
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)    = fluxes%evap(i,j)    * G%mask2dT(i,j)
      fluxes%lprec(i,j)   = fluxes%lprec(i,j)   * G%mask2dT(i,j)
      fluxes%fprec(i,j)   = fluxes%fprec(i,j)   * G%mask2dT(i,j)
      fluxes%lrunoff(i,j) = fluxes%lrunoff(i,j) * G%mask2dT(i,j)
      fluxes%frunoff(i,j) = fluxes%frunoff(i,j) * G%mask2dT(i,j)
      fluxes%lw(i,j)      = fluxes%lw(i,j)      * G%mask2dT(i,j)
      fluxes%sens(i,j)    = fluxes%sens(i,j)    * G%mask2dT(i,j)
      fluxes%sw(i,j)      = fluxes%sw(i,j)      * G%mask2dT(i,j)
      fluxes%latent(i,j)  = fluxes%latent(i,j)  * G%mask2dT(i,j)

      fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
      fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
    enddo ; enddo
  endif ! fluxes have changed and need to be masked


  ! restoring surface boundary fluxes
  if (CS%restorebuoy) then

    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const_T)
          fluxes%vprec(i,j) = - (CS%rho_restore*CS%Flux_const_S) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
                             (CS%G_Earth * CS%Flux_const / CS%rho_restore)
        else
          fluxes%buoy(i,j) = 0.0
        endif
      enddo ; enddo
    endif

  else                                              ! not RESTOREBUOY
    if (.not.CS%use_temperature) then
      call MOM_error(FATAL, "buoyancy_forcing in MOM_surface_forcing: "// &
                     "The fluxes need to be defined without RESTOREBUOY.")
    endif

  endif                                             ! end RESTOREBUOY

!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!#CTRL#                           fluxes%vprec, day, dt, G, US, CS%ctrl_forcing_CSp)
!#CTRL# endif

  call callTree_leave("buoyancy_forcing_from_files")
end procedure buoyancy_forcing_from_files
module procedure buoyancy_forcing_from_data_override
  real :: rhoXcp ! The mean density times the heat capacity [Q R C-1 ~> J m-3 degC-1].
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  call callTree_enter("buoyancy_forcing_from_data_override, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (CS%use_temperature) rhoXcp = CS%rho_restore * fluxes%C_p

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(G%Domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  call data_override(G%Domain, 'lw', fluxes%lw, day, scale=US%W_m2_to_QRZ_T)
  call data_override(G%Domain, 'sw', fluxes%sw, day, scale=US%W_m2_to_QRZ_T)

  ! The normal MOM6 sign conventions are that fluxes%evap and fluxes%sens are positive into the
  ! ocean but evap and sens are normally positive quantities in the files.
  call data_override(G%Domain, 'evap', fluxes%evap, day, scale=-US%kg_m2s_to_RZ_T)
  call data_override(G%Domain, 'sens', fluxes%sens, day, scale=-US%W_m2_to_QRZ_T)

  do j=js,je ; do i=is,ie
    fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
    fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
  enddo ; enddo

  call data_override(G%Domain, 'snow', fluxes%fprec, day, scale=US%kg_m2s_to_RZ_T)
  call data_override(G%Domain, 'rain', fluxes%lprec, day, scale=US%kg_m2s_to_RZ_T)
  call data_override(G%Domain, 'runoff', fluxes%lrunoff, day, scale=US%kg_m2s_to_RZ_T)
  call data_override(G%Domain, 'calving', fluxes%frunoff, day, scale=US%kg_m2s_to_RZ_T)

!     Read the SST and SSS fields for damping.
  if (CS%restorebuoy) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
    call data_override(G%Domain, 'SST_restore', CS%T_restore, day, scale=US%degC_to_C)
    call data_override(G%Domain, 'SSS_restore', CS%S_restore, day, scale=US%ppt_to_S)
  endif

  ! restoring boundary fluxes
  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const_T)
          fluxes%vprec(i,j) = - (CS%rho_restore*CS%Flux_const_S) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
                             (CS%G_Earth * CS%Flux_const / CS%rho_restore)
        else
          fluxes%buoy(i,j) = 0.0
        endif
      enddo ; enddo
    endif
  else                                              ! not RESTOREBUOY
    if (.not.CS%use_temperature) then
      call MOM_error(FATAL, "buoyancy_forcing in MOM_surface_forcing: "// &
                     "The fluxes need to be defined without RESTOREBUOY.")
    endif
  endif                                             ! end RESTOREBUOY


  ! mask out land points and compute heat content of water fluxes
  ! assume liquid precip enters ocean at SST
  ! assume frozen precip enters ocean at 0degC
  ! assume liquid runoff enters ocean at SST
  ! assume solid runoff (calving) enters ocean at 0degC
  ! mass leaving ocean has heat_content determined in MOM_diabatic_driver.F90
  do j=js,je ; do i=is,ie
    fluxes%evap(i,j)    = fluxes%evap(i,j)    * G%mask2dT(i,j)
    fluxes%lprec(i,j)   = fluxes%lprec(i,j)   * G%mask2dT(i,j)
    fluxes%fprec(i,j)   = fluxes%fprec(i,j)   * G%mask2dT(i,j)
    fluxes%lrunoff(i,j) = fluxes%lrunoff(i,j) * G%mask2dT(i,j)
    fluxes%frunoff(i,j) = fluxes%frunoff(i,j) * G%mask2dT(i,j)
    fluxes%lw(i,j)      = fluxes%lw(i,j)      * G%mask2dT(i,j)
    fluxes%latent(i,j)  = fluxes%latent(i,j)  * G%mask2dT(i,j)
    fluxes%sens(i,j)    = fluxes%sens(i,j)    * G%mask2dT(i,j)
    fluxes%sw(i,j)      = fluxes%sw(i,j)      * G%mask2dT(i,j)

    fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
    fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
    fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
  enddo ; enddo

!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!#CTRL#                           fluxes%vprec, day, dt, G, US, CS%ctrl_forcing_CSp)
!#CTRL# endif

  call callTree_leave("buoyancy_forcing_from_data_override")
end procedure buoyancy_forcing_from_data_override
module procedure buoyancy_forcing_zero
  integer :: i, j, is, ie, js, je
  call callTree_enter("buoyancy_forcing_zero, MOM_surface_forcing.F90")
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  call callTree_leave("buoyancy_forcing_zero")
end procedure buoyancy_forcing_zero
module procedure buoyancy_forcing_const
  integer :: i, j, is, ie, js, je
  call callTree_enter("buoyancy_forcing_const, MOM_surface_forcing.F90")
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec

  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = CS%constantHeatForcing * G%mask2dT(i,j)
      fluxes%sw(i,j)                   = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  call callTree_leave("buoyancy_forcing_const")
end procedure buoyancy_forcing_const
module procedure buoyancy_forcing_linear
  real :: y             ! The latitude relative to the south normalized by the domain extent [nondim]
  real :: T_restore     ! The temperature towards which to restore [C ~> degC]
  real :: S_restore     ! The salinity towards which to restore [S ~> ppt]
  integer :: i, j, is, ie, js, je
  call callTree_enter("buoyancy_forcing_linear, MOM_surface_forcing.F90")
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! This case has no surface buoyancy forcing.
  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        y = (G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat
        T_restore = CS%T_south + (CS%T_north-CS%T_south)*y
        S_restore = CS%S_south + (CS%S_north-CS%S_south)*y
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((T_Restore - sfc_state%SST(i,j)) * ((CS%rho_restore * fluxes%C_p) * CS%Flux_const))
          fluxes%vprec(i,j) = - (CS%rho_restore*CS%Flux_const) * &
              (S_Restore - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + S_Restore))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "RESTOREBUOY to linear not written yet.")
     !do j=js,je ; do i=is,ie
     !  if (G%mask2dT(i,j) > 0.0) then
     !    fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
     !                       (CS%G_Earth * CS%Flux_const / CS%rho_restore)
     !  else
     !    fluxes%buoy(i,j) = 0.0
     !  endif
     !enddo ; enddo
    endif
  else                                              ! not RESTOREBUOY
    if (.not.CS%use_temperature) then
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "The fluxes need to be defined without RESTOREBUOY.")
    endif
  endif                                             ! end RESTOREBUOY

  call callTree_leave("buoyancy_forcing_linear")
end procedure buoyancy_forcing_linear
module procedure get_file_time_level
  integer :: days, seconds           ! The number of days and seconds since the start of the calendar
  integer :: year, month, day, hour, minute, second ! The components of the model time
  integer :: recs_per_day            ! The number of file time records per day
  integer :: recs                    ! The number of time levels into the file to read without
  if ( (days_per_rec >= 1000000) .or. &
       ( (days_per_rec == 0) .and. .not.((nlev_file == 12) .or. (nlev_file == 365)) ) ) then
    ! The second condition above is to recreate the existing behavior, but it should perhaps be
    ! phased out.
    time_lev = 1
  elseif ( (days_per_rec == 31) .or. ((days_per_rec == 0) .and. (nlev_file == 12)) ) then
    call get_date(Time, year, month, day, hour, minute, second)
    time_lev = month
  else
    call get_time(Time, seconds, days)
    if ( (days_per_rec == 0) .or. (abs(days_per_rec) == 1) ) then
      recs = days
    elseif (days_per_rec < 0) then
      recs_per_day = -days_per_rec
      recs = days * recs_per_day + ( (recs_per_day*set_time(seconds, 0)) / set_time(0, 1) )
      ! When integer rounding in the time-type arithmetic is considered, the line above is equivalent to:
      !   seconds_per_day = set_time(0, 1) / set_time(1, 0)
      !   recs = days * recs_per_day + floor(real(recs_per_day*seconds) / real(seconds_per_day))
    else
      recs = days / days_per_rec
    endif
    time_lev = recs - nlev_file*floor(real(recs) / real(nlev_file)) + 1
  endif

end procedure get_file_time_level
module procedure MARBL_forcing_from_data_override
  real, pointer, dimension(:,:) :: atm_co2_prog         =>NULL() !< Prognostic atmospheric CO2 concentration [ppm]
  real, pointer, dimension(:,:) :: atm_co2_diag         =>NULL() !< Diagnostic atmospheric CO2 concentration [ppm]
  real, pointer, dimension(:,:) :: atm_fine_dust_flux   =>NULL() !< Fine dust flux from atmosphere
  real, pointer, dimension(:,:) :: atm_coarse_dust_flux =>NULL() !< Coarse dust flux from atmosphere
  real, pointer, dimension(:,:) :: seaice_dust_flux     =>NULL() !< Dust flux from seaice
  real, pointer, dimension(:,:) :: atm_bc_flux          =>NULL() !< Black carbon flux from atmosphere
  real, pointer, dimension(:,:) :: seaice_bc_flux       =>NULL() !< Black carbon flux from seaice
  real, pointer, dimension(:,:) :: nhx_dep              =>NULL() !< Nitrogen deposition
  real, pointer, dimension(:,:) :: noy_dep              =>NULL() !< Nitrogen deposition
  integer :: isc, iec, jsc, jec
  real, pointer, dimension(:,:)   :: afracr        =>NULL()
  real, pointer, dimension(:,:)   :: swnet_afracr  =>NULL()
  real, pointer, dimension(:,:,:) :: swpen_ifrac_n =>NULL()
  real, pointer, dimension(:,:,:) :: ifrac_n       =>NULL()
  call callTree_enter("MARBL_forcing_from_data_override, MOM_surface_forcing.F90")

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(G%Domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  ! Allocate memory for pointers
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  allocate ( atm_co2_prog   (isc:iec,jsc:jec),       &
             atm_co2_diag   (isc:iec,jsc:jec),       &
             atm_fine_dust_flux (isc:iec,jsc:jec),   &
             atm_coarse_dust_flux (isc:iec,jsc:jec), &
             seaice_dust_flux (isc:iec,jsc:jec),     &
             atm_bc_flux (isc:iec,jsc:jec),          &
             seaice_bc_flux (isc:iec,jsc:jec),       &
             nhx_dep (isc:iec,jsc:jec),              &
             noy_dep (isc:iec,jsc:jec),              &
             source=0.0)


  ! fluxes used directly as MARBL inputs
  ! (should be scaled)
  call data_override(G%Domain, 'ice_fraction', fluxes%ice_fraction, day)
  call data_override(G%Domain, 'u10_sqr', fluxes%u10_sqr, day, scale=US%m_s_to_L_T**2)

  ! fluxes used to compute MARBL inputs
  ! These are kept in physical units, and will be scaled appropriately in
  ! convert_driver_fields_to_forcings()
  call data_override(G%Domain, 'atm_co2_prog', atm_co2_prog, day)
  call data_override(G%Domain, 'atm_co2_diag', atm_co2_diag, day)
  call data_override(G%Domain, 'atm_fine_dust_flux', atm_fine_dust_flux, day)
  call data_override(G%Domain, 'atm_coarse_dust_flux', atm_coarse_dust_flux, day)
  call data_override(G%Domain, 'atm_bc_flux', atm_bc_flux, day)
  call data_override(G%Domain, 'seaice_dust_flux', seaice_dust_flux, day)
  call data_override(G%Domain, 'seaice_bc_flux', seaice_bc_flux, day)
  call data_override(G%Domain, 'nhx_dep', nhx_dep, day)
  call data_override(G%Domain, 'noy_dep', noy_dep, day)

  call convert_driver_fields_to_forcings(atm_fine_dust_flux, atm_coarse_dust_flux, &
                                         seaice_dust_flux, atm_bc_flux, seaice_bc_flux, &
                                         nhx_dep, noy_dep, atm_co2_prog, atm_co2_diag, &
                                         afracr, swnet_afracr, ifrac_n, swpen_ifrac_n, &
                                         day, G, US, 0, 0, fluxes, CS%marbl_forcing_CSp)

  deallocate ( atm_co2_prog,         &
               atm_co2_diag,         &
               atm_fine_dust_flux,   &
               atm_coarse_dust_flux, &
               seaice_dust_flux,     &
               atm_bc_flux,          &
               seaice_bc_flux,       &
               nhx_dep,              &
               noy_dep)

  call callTree_leave("MARBL_forcing_from_data_override")

end procedure MARBL_forcing_from_data_override
module procedure forcing_save_restart
  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end procedure forcing_save_restart
module procedure surface_forcing_init
  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
# include "version_variable.h"
  real :: flux_const_default ! The unscaled value of FLUXCONST [m day-1]
  logical :: Boussinesq       ! If true, this run is fully Boussinesq
  logical :: semi_Boussinesq  ! If true, this run is partially non-Boussinesq
  logical :: fix_ustar_gustless_bug  ! If false, include a bug using an older run-time parameter.
  logical :: test_value  ! This is used to determine whether a logical parameter is being set explicitly.
  logical :: explicit_bug, explicit_fix ! These indicate which parameters are set explicitly.
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  character(len=40)  :: mdl = "MOM_surface_forcing" ! This module's name.
  character(len=200) :: filename, gust_file ! The name of the gustiness input file.
  if (associated(CS)) then
    call MOM_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  id_clock_forcing=cpu_clock_id('(Ocean surface forcing)', grain=CLOCK_MODULE)
  call cpu_clock_begin(id_clock_forcing)

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, '')
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, "MOM", "BOUSSINESQ", Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true., do_not_log=.true.)
  call get_param(param_file, "MOM", "SEMI_BOUSSINESQ", semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=.true.)
  CS%nonBous = .not.(Boussinesq .or. semi_Boussinesq)
  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mdl, "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is "//&
                 "true. This assumes that KD = KDML = 0.0 and that "//&
                 "there is no buoyancy forcing, but makes the model "//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, mdl, "VARIABLE_WINDS", CS%variable_winds, &
                 "If true, the winds vary in time after the initialization.", &
                 default=.true.)
  call get_param(param_file, mdl, "VARIABLE_BUOYFORCE", CS%variable_buoyforce, &
                 "If true, the buoyancy forcing varies in time after the "//&
                 "initialization of the model.", default=.true.)

  ! Determine parameters related to the buoyancy forcing.
  call get_param(param_file, mdl, "BUOY_CONFIG", CS%buoy_config, &
                 "The character string that indicates how buoyancy forcing is specified.  Valid "//&
                 "options include (file), (data_override), (zero), (const), (linear), (MESO), "//&
                 "(SCM_CVmix_tests), (BFB), (dumbbell), (USER) and (NONE).", default="zero")
  if (trim(CS%buoy_config) == "file") then
    call get_param(param_file, mdl, "ARCHAIC_OMIP_FORCING_FILE", CS%archaic_OMIP_file, &
                 "If true, use the forcing variable decomposition from "//&
                 "the old German OMIP prescription that predated CORE. If "//&
                 "false, use the variable groupings available from MOM "//&
                 "output diagnostics of forcing variables.", default=.true.)
    if (CS%archaic_OMIP_file) then
      call get_param(param_file, mdl, "LONGWAVEDOWN_FILE", CS%longwave_file, &
                 "The file with the downward longwave heat flux, in "//&
                 "variable lwdn_sfc.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LONGWAVEUP_FILE", CS%longwaveup_file, &
                 "The file with the upward longwave heat flux, in "//&
                 "variable lwup_sfc.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in "//&
                 "variable evap.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in "//&
                 "variable shflx.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SHORTWAVEUP_FILE", CS%shortwaveup_file, &
                 "The file with the upward shortwave heat flux.", &
                 fail_if_missing=.true.)
      call get_param(param_file, mdl, "SHORTWAVEDOWN_FILE", CS%shortwave_file, &
                 "The file with the downward shortwave heat flux.", &
                 fail_if_missing=.true.)
      call get_param(param_file, mdl, "SNOW_FILE", CS%snow_file, &
                 "The file with the downward frozen precip flux, in "//&
                 "variable snow.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "PRECIP_FILE", CS%rain_file, &
                 "The file with the downward total precip flux, in "//&
                 "variable precip.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "FRESHDISCHARGE_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving fluxes, "//&
                 "invariables disch_w and disch_s.", fail_if_missing=.true.)

      ! These variable names are hard-coded, per the archaic OMIP conventions.
      CS%latentheat_file = CS%evaporation_file ; CS%latent_var = "evap"
      CS%LW_var = "lwdn_sfc" ; CS%SW_var = "swdn_sfc" ; CS%sens_var = "shflx"
      CS%evap_var = "evap" ; CS%rain_var = "precip" ; CS%snow_var = "snow"
      CS%lrunoff_var = "disch_w" ; CS%frunoff_var = "disch_s"

    else
      call get_param(param_file, mdl, "LONGWAVE_FILE", CS%longwave_file, &
                 "The file with the longwave heat flux, in the variable "//&
                 "given by LONGWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LONGWAVE_FORCING_VAR", CS%LW_var, &
                 "The variable with the longwave forcing field.", default="LW")
      call get_param(param_file, mdl, "LONGWAVE_FILE_DAYS_PER_RECORD", CS%LW_days_per_rec, &
                 "If positive the number of days of longwave fluxes per time level in LONGWAVE_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=0)

      call get_param(param_file, mdl, "SHORTWAVE_FILE", CS%shortwave_file, &
                 "The file with the shortwave heat flux, in the variable "//&
                 "given by SHORTWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SHORTWAVE_FORCING_VAR", CS%SW_var, &
                 "The variable with the shortwave forcing field.", default="SW")
      call get_param(param_file, mdl, "SHORTWAVE_FILE_DAYS_PER_RECORD", CS%SW_days_per_rec, &
                 "If positive the number of days of shortwave fluxes per time level in SHORTWAVE_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)

      call get_param(param_file, mdl, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in the "//&
                 "variable given by EVAP_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "EVAP_FORCING_VAR", CS%evap_var, &
                 "The variable with the evaporative moisture flux.", &
                 default="evap")
      call get_param(param_file, mdl, "EVAPORATION_FILE_DAYS_PER_RECORD", CS%evap_days_per_rec, &
                 "If positive the number of days of evaporation per time level in EVAPORATION_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)

      call get_param(param_file, mdl, "LATENTHEAT_FILE", CS%latentheat_file, &
                 "The file with the latent heat flux, in the variable "//&
                 "given by LATENT_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LATENT_FORCING_VAR", CS%latent_var, &
                 "The variable with the latent heat flux.", default="latent")
      call get_param(param_file, mdl, "LATENTHEAT_FILE_DAYS_PER_RECORD", CS%latent_days_per_rec, &
                 "If positive the number of days of latent heat fluxes per time level in LATENTHEAT_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)

      call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in the variable "//&
                 "given by SENSIBLE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SENSIBLE_FORCING_VAR", CS%sens_var, &
                 "The variable with the sensible heat flux.", default="sensible")
      call get_param(param_file, mdl, "SENSIBLEHEAT_FILE_DAYS_PER_RECORD", CS%sens_days_per_rec, &
                 "If positive the number of days of sensible heat fluxes per time level in SENSIBLEHEAT_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)

      call get_param(param_file, mdl, "RAIN_FILE", CS%rain_file, &
                 "The file with the liquid precipitation flux, in the "//&
                 "variable given by RAIN_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "RAIN_FORCING_VAR", CS%rain_var, &
                 "The variable with the liquid precipitation flux.", &
                 default="liq_precip")
      call get_param(param_file, mdl, "RAIN_FILE_DAYS_PER_RECORD", CS%precip_days_per_rec, &
                 "If positive the number of days of rain fluxes per time level in RAIN_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)
      call get_param(param_file, mdl, "SNOW_FILE", CS%snow_file, &
                 "The file with the frozen precipitation flux, in the "//&
                 "variable given by SNOW_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SNOW_FORCING_VAR", CS%snow_var, &
                 "The variable with the frozen precipitation flux.", &
                 default="froz_precip")
      call get_param(param_file, mdl, "SHORTWAVE_FILE_DAYS_PER_RECORD", CS%SW_days_per_rec, &
                 "If positive the number of days of shortwave fluxes per time level in SHORTWAVE_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%LW_days_per_rec)

      call get_param(param_file, mdl, "RUNOFF_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving "//&
                 "fluxes, in variables given by LIQ_RUNOFF_FORCING_VAR "//&
                 "and FROZ_RUNOFF_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LIQ_RUNOFF_FORCING_VAR", CS%lrunoff_var, &
                 "The variable with the liquid runoff flux.", &
                 default="liq_runoff")
      call get_param(param_file, mdl, "FROZ_RUNOFF_FORCING_VAR", CS%frunoff_var, &
                 "The variable with the frozen runoff flux.", &
                 default="froz_runoff")
      call get_param(param_file, mdl, "RUNOFF_FILE_DAYS_PER_RECORD", CS%SW_days_per_rec, &
                 "If positive the number of days of runoff per time level in RUNOFF_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=0)
    endif

    call get_param(param_file, mdl, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in the "//&
                 "variable given by SST_RESTORE_VAR.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to "//&
                 "restore in the variable given by SSS_RESTORE_VAR.", &
                 fail_if_missing=.true.)
    if (CS%archaic_OMIP_file) then
      CS%SST_restore_var = "TEMP" ; CS%SSS_restore_var = "SALT"
    else
      call get_param(param_file, mdl, "SST_RESTORE_VAR", CS%SST_restore_var, &
                 "The variable with the SST toward which to restore.", &
                 default="SST")
      call get_param(param_file, mdl, "SSTRESTORE_FILE_DAYS_PER_RECORD", CS%SST_days_per_rec, &
                 "If positive the number of days of SST per time level in SSTRESTORE_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=0)
      call get_param(param_file, mdl, "SSS_RESTORE_VAR", CS%SSS_restore_var, &
                 "The variable with the SSS toward which to restore.", &
                 default="SSS")
      call get_param(param_file, mdl, "SALINITYRESTORE_FILE_DAYS_PER_RECORD", CS%SSS_days_per_rec, &
                 "If positive the number of days of salinity per time level in SALINITYRESTORE_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=CS%SST_days_per_rec)
    endif

    ! Add inputdir to the file names.
    CS%shortwave_file = trim(CS%inputdir)//trim(CS%shortwave_file)
    CS%longwave_file = trim(CS%inputdir)//trim(CS%longwave_file)
    CS%sensibleheat_file = trim(CS%inputdir)//trim(CS%sensibleheat_file)
    CS%latentheat_file = trim(CS%inputdir)//trim(CS%latentheat_file)
    CS%evaporation_file = trim(CS%inputdir)//trim(CS%evaporation_file)
    CS%snow_file = trim(CS%inputdir)//trim(CS%snow_file)
    CS%rain_file = trim(CS%inputdir)//trim(CS%rain_file)
    CS%runoff_file = trim(CS%inputdir)//trim(CS%runoff_file)

    CS%shortwaveup_file = trim(CS%inputdir)//trim(CS%shortwaveup_file)
    CS%longwaveup_file = trim(CS%inputdir)//trim(CS%longwaveup_file)

    CS%SSTrestore_file = trim(CS%inputdir)//trim(CS%SSTrestore_file)
    CS%salinityrestore_file = trim(CS%inputdir)//trim(CS%salinityrestore_file)
  elseif (trim(CS%buoy_config) == "const") then
    call get_param(param_file, mdl, "SENSIBLE_HEAT_FLUX", CS%constantHeatForcing, &
                 "A constant heat forcing (positive into ocean) applied "//&
                 "through the sensible heat flux field. ", &
                 units='W/m2', scale=US%W_m2_to_QRZ_T, fail_if_missing=.true.)
  endif

  ! Determine parameters related to the wind forcing.
  call get_param(param_file, mdl, "WIND_CONFIG", CS%wind_config, &
                 "The character string that indicates how wind forcing is specified.  Valid "//&
                 "options include (file), (data_override), (2gyre), (1gyre), (gyres), (zero), "//&
                 "(const), (Neverworld), (scurves), (ideal_hurr), (SCM_CVmix_tests) and (USER).", &
                 default="zero")
  if (trim(CS%wind_config) == "file") then
    call get_param(param_file, mdl, "WIND_FILE", CS%wind_file, &
                 "The file in which the wind stresses are found in "//&
                 "variables STRESS_X and STRESS_Y.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "WINDSTRESS_X_VAR",CS%stress_x_var, &
                 "The name of the x-wind stress variable in WIND_FILE.", &
                 default="STRESS_X")
    call get_param(param_file, mdl, "WINDSTRESS_Y_VAR", CS%stress_y_var, &
                 "The name of the y-wind stress variable in WIND_FILE.", &
                 default="STRESS_Y")
    call get_param(param_file, mdl, "WIND_STAGGER",CS%wind_stagger, &
                 "A character indicating how the wind stress components "//&
                 "are staggered in WIND_FILE.  This may be A or C for now.", &
                 default="C")
    call get_param(param_file, mdl, "WINDSTRESS_SCALE", CS%wind_scale, &
                 "A value by which the wind stresses in WIND_FILE are rescaled.", &
                 default=1.0, units="nondim")
    call get_param(param_file, mdl, "USTAR_FORCING_VAR", CS%ustar_var, &
                 "The name of the friction velocity variable in WIND_FILE "//&
                 "or blank to get ustar from the wind stresses plus the "//&
                 "gustiness.", default=" ")
    CS%wind_file = trim(CS%inputdir) // trim(CS%wind_file)
    call get_param(param_file, mdl, "WIND_FILE_DAYS_PER_RECORD", CS%wind_days_per_rec, &
                 "If positive the number of days of wind stress per time level in WIND_FILE, "//&
                 "or if negative the number of time levels per day.  If 31 change forcing monthly, "//&
                 "or if 0 the model will guess the right value based on the file size.", &
                 default=0)
  endif
  if (trim(CS%wind_config) == "gyres") then
    call get_param(param_file, mdl, "TAUX_CONST", CS%gyres_taux_const, &
                 "With the gyres wind_config, the constant offset in the "//&
                 "zonal wind stress profile: "//&
                 "  A in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2)
    call get_param(param_file, mdl, "TAUX_SIN_AMP", CS%gyres_taux_sin_amp, &
                 "With the gyres wind_config, the sine amplitude in the "//&
                 "zonal wind stress profile: "//&
                 "  B in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2)
    call get_param(param_file, mdl, "TAUX_COS_AMP", CS%gyres_taux_cos_amp, &
                 "With the gyres wind_config, the cosine amplitude in "//&
                 "the zonal wind stress profile: "//&
                 "  C in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2)
    call get_param(param_file, mdl, "TAUX_N_PIS",CS%gyres_taux_n_pis, &
                 "With the gyres wind_config, the number of gyres in "//&
                 "the zonal wind stress profile: "//&
                 "  n in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="nondim", default=0.0)
    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
    call get_param(param_file, mdl, "WIND_GYRES_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the expressions used to set gyre wind stresses. "//&
                 "Values below 20190101 recover the answers from the end of 2018, "//&
                 "while higher values use a form of the gyre wind stresses that are "//&
                 "rotationally invariant and more likely to be the same between compilers.", &
                 default=default_answer_date)
  else
    CS%answer_date = 20190101
  endif
  if (trim(CS%wind_config) == "scurves") then
    call get_param(param_file, mdl, "WIND_SCURVES_LATS", CS%scurves_ydata, &
                 "A list of latitudes defining a piecewise scurve profile "//&
                 "for zonal wind stress.", &
                 units="degrees N", fail_if_missing=.true.)
    call get_param(param_file, mdl, "WIND_SCURVES_TAUX", CS%scurves_taux, &
                 "A list of zonal wind stress values at latitudes "//&
                 "WIND_SCURVES_LATS defining a piecewise scurve profile.", &
                 units="Pa", scale=US%Pa_to_RLZ_T2, fail_if_missing=.true.)
  endif
  if (trim(CS%wind_config) == "2gyre") then
    call get_param(param_file, mdl, "TAUX_MAGNITUDE", CS%taux_mag, &
                 "The peak zonal wind stress when WIND_CONFIG = 2gyre.", &
                 units="Pa", default=0.1, scale=US%Pa_to_RLZ_T2)
  endif
  if (trim(CS%wind_config) == "1gyre") then
    call get_param(param_file, mdl, "TAUX_MAGNITUDE", CS%taux_mag, &
                 "The peak zonal wind stress when WIND_CONFIG = 1gyre.", &
                 units="Pa", default=-0.2, scale=US%Pa_to_RLZ_T2)
  endif
  if (trim(CS%wind_config) == "Neverworld" .or. trim(CS%wind_config) == "Neverland") then
    call get_param(param_file, mdl, "TAUX_MAGNITUDE", CS%taux_mag, &
                 "The peak zonal wind stress when WIND_CONFIG = Neverworld.", &
                 units="Pa", default=0.2, scale=US%Pa_to_RLZ_T2)
  endif

  if ((trim(CS%wind_config) == "2gyre") .or. &
      (trim(CS%wind_config) == "1gyre") .or. &
      (trim(CS%wind_config) == "gyres") .or. &
      (trim(CS%buoy_config) == "linear")) then
    CS%south_lat = G%south_lat
    CS%len_lat = G%len_lat
  endif

  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R) ! (, do_not_log=CS%nonBous)
  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back toward some "//&
                 "specified surface state with a rate given by FLUXCONST.", default=.false.)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", default=hlf, &
                 units="J/kg", scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", default=hlv, units="J/kg", scale=US%J_kg_to_Q)
  if (CS%restorebuoy) then
    ! These three variables use non-standard time units, but are rescaled as they are read.
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0)

    if (CS%use_temperature) then
      call get_param(param_file, mdl, "FLUXCONST", flux_const_default, &
                 default=0.0, units="m day-1", do_not_log=.true.)
      call get_param(param_file, mdl, "FLUXCONST_T", CS%Flux_const_T, &
                 "The constant that relates the restoring surface temperature flux to the "//&
                 "relative surface anomaly (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, default=flux_const_default)
      call get_param(param_file, mdl, "FLUXCONST_S", CS%Flux_const_S, &
                 "The constant that relates the restoring surface salinity flux to the "//&
                 "relative surface anomaly (akin to a piston velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, default=flux_const_default)
    endif

    if (trim(CS%buoy_config) == "linear") then
      call get_param(param_file, mdl, "SST_NORTH", CS%T_north, &
                 "With buoy_config linear, the sea surface temperature "//&
                 "at the northern end of the domain toward which to "//&
                 "to restore.", units="degC", default=0.0, scale=US%degC_to_C)
      call get_param(param_file, mdl, "SST_SOUTH", CS%T_south, &
                 "With buoy_config linear, the sea surface temperature "//&
                 "at the southern end of the domain toward which to "//&
                 "to restore.", units="degC", default=0.0, scale=US%degC_to_C)
      call get_param(param_file, mdl, "SSS_NORTH", CS%S_north, &
                 "With buoy_config linear, the sea surface salinity "//&
                 "at the northern end of the domain toward which to "//&
                 "to restore.", units="ppt", default=35.0, scale=US%ppt_to_S)
      call get_param(param_file, mdl, "SSS_SOUTH", CS%S_south, &
                 "With buoy_config linear, the sea surface salinity "//&
                 "at the southern end of the domain toward which to "//&
                 "to restore.", units="ppt", default=35.0, scale=US%ppt_to_S)
    endif
    call get_param(param_file, mdl, "RESTORE_FLUX_RHO", CS%rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat "//&
                 "fluxes with RESTORE_SALINITY or RESTORE_TEMPERATURE.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(((CS%Flux_const==0.0).and.(CS%Flux_const_T==0.0).and.(CS%Flux_const_S==0.0))&
                            .or.(.not.CS%restorebuoy)))
  endif
  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)

  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2*US%L_to_Z)

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

  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.", fail_if_missing=.true.)
    call safe_alloc_ptr(CS%gust,G%isd,G%ied,G%jsd,G%jed)
    filename = trim(CS%inputdir) // trim(gust_file)
    ! NOTE: There are certain cases where FMS is unable to read this file, so
    ! we use read_netCDF_data in place of MOM_read_data.
    call read_netCDF_data(filename, 'gustiness', CS%gust, G%Domain, &
                          rescale=US%Pa_to_RLZ_T2*US%L_to_Z) ! units in file should be [Pa]
  endif
  call get_param(param_file, mdl, "USE_MARBL_TRACERS", CS%use_marbl_tracers, &
                  default=.false., do_not_log=.true.)

!  All parameter settings are now known.

  if (trim(CS%wind_config) == "USER" .or. trim(CS%buoy_config) == "USER" ) then
    call USER_surface_forcing_init(Time, G, US, param_file, diag, CS%user_forcing_CSp)
  elseif (trim(CS%buoy_config) == "BFB" ) then
    call BFB_surface_forcing_init(Time, G, US, param_file, diag, CS%BFB_forcing_CSp)
  elseif (trim(CS%buoy_config) == "dumbbell" ) then
    call dumbbell_surface_forcing_init(Time, G, US, param_file, diag, CS%dumbbell_forcing_CSp)
  elseif (trim(CS%wind_config) == "MESO" .or. trim(CS%buoy_config) == "MESO" ) then
    call MESO_surface_forcing_init(Time, G, US, param_file, diag, CS%MESO_forcing_CSp)
  elseif (trim(CS%wind_config) == "ideal_hurr") then
    call idealized_hurricane_wind_init(Time, G, US, param_file, CS%idealized_hurricane_CSp)
  elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
    call MOM_error(FATAL, "MOM_surface_forcing (surface_forcing_init): "//&
          'WIND_CONFIG = "SCM_ideal_hurr" is a depricated option.  '//&
          'To obtain mathematically equivalent results set '//&
          'WIND_CONFIG = "ideal_hurr", IDL_HURR_SCM = True and IDL_HURR_X0 = 6.48e+05.')
  elseif (trim(CS%wind_config) == "const") then
    call get_param(param_file, mdl, "CONST_WIND_TAUX", CS%tau_x0, &
                 "With wind_config const, this is the constant zonal wind-stress", &
                 units="Pa", scale=US%Pa_to_RLZ_T2, fail_if_missing=.true.)
    call get_param(param_file, mdl, "CONST_WIND_TAUY", CS%tau_y0, &
                 "With wind_config const, this is the constant meridional wind-stress", &
                 units="Pa", scale=US%Pa_to_RLZ_T2, fail_if_missing=.true.)
  elseif (trim(CS%wind_config) == "SCM_CVmix_tests" .or. &
          trim(CS%buoy_config) == "SCM_CVmix_tests") then
    call SCM_CVmix_tests_surface_forcing_init(Time, G, param_file, CS%SCM_CVmix_tests_CSp)
  endif

  ! Set up MARBL forcing control structure
  call MARBL_forcing_init(G, US, param_file, diag, Time, CS%inputdir, CS%use_marbl_tracers, &
      CS%marbl_forcing_CSp)

  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles)

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

  ! Determine how many time levels are in each forcing variable.
  if (trim(CS%buoy_config) == "file") then
    CS%SW_nlev = num_timelevels(CS%shortwave_file, CS%SW_var, min_dims=3)
    CS%LW_nlev = num_timelevels(CS%longwave_file, CS%LW_var, min_dims=3)
    CS%latent_nlev = num_timelevels(CS%latentheat_file, CS%latent_var, 3)
    CS%sens_nlev = num_timelevels(CS%sensibleheat_file, CS%sens_var, min_dims=3)

    CS%evap_nlev = num_timelevels(CS%evaporation_file, CS%evap_var, min_dims=3)
    CS%precip_nlev = num_timelevels(CS%rain_file, CS%rain_var, min_dims=3)
    CS%runoff_nlev = num_timelevels(CS%runoff_file, CS%lrunoff_var, 3)

    CS%SST_nlev = num_timelevels(CS%SSTrestore_file, CS%SST_restore_var, 3)
    CS%SSS_nlev = num_timelevels(CS%salinityrestore_file, CS%SSS_restore_var, 3)
  endif

  if (trim(CS%wind_config) == "file") &
    CS%wind_nlev = num_timelevels(CS%wind_file, CS%stress_x_var, min_dims=3)

!#CTRL#  call controlled_forcing_init(Time, G, US, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end procedure surface_forcing_init
module procedure surface_forcing_end
  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!#CTRL#  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

  call callTree_leave("MARBL_forcing_from_data_override, MOM_surface_forcing.F90")
end procedure surface_forcing_end
end submodule MOM_surface_forcing_s
