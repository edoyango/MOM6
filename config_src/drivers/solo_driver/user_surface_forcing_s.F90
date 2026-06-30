submodule (user_surface_forcing) user_surface_forcing_s
  implicit none
contains
module procedure USER_wind_forcing
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  call MOM_error(FATAL, "User_wind_surface_forcing: " // &
     "User forcing routine called without modification." )

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., tau_mag=.true.)

  !  Set the surface wind stresses, in units of [R L Z T-2 ~> Pa].  A positive taux
  !  accelerates the ocean to the (pseudo-)east.

  !  The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    ! Change this to the desired expression.
    forces%taux(I,j) = G%mask2dCu(I,j) * 0.0*US%Pa_to_RLZ_T2
  enddo ; enddo
  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0  ! Change this to the desired expression.
  enddo ; enddo

  !    Set the surface friction velocity, in units of [Z T-1 ~> m s-1].  ustar
  !  is always positive.
  if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    forces%tau_mag(i,j) = G%mask2dT(i,j) * (CS%gust_const + &
            US%L_to_Z*sqrt(0.5*((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2)) + &
                           0.5*((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2))))
    if (associated(forces%ustar)) &
      forces%ustar(i,j) = G%mask2dT(i,j) * sqrt(forces%tau_mag(i,j) * (1.0/CS%Rho0))
  enddo ; enddo ; endif

end procedure USER_wind_forcing
module procedure USER_buoyancy_forcing
  real :: Temp_restore   ! The temperature that is being restored toward [C ~> degC].
  real :: Salin_restore  ! The salinity that is being restored toward [S ~> ppt]
  real :: density_restore  ! The potential density that is being restored
  real :: rhoXcp         ! The mean density times the heat capacity [Q R C-1 ~> J m-3 degC-1].
  real :: buoy_rest_const  ! A constant relating density anomalies to the
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
    "User forcing routine called without modification." )

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call safe_alloc_ptr(fluxes%evap, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lprec, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%fprec, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lrunoff, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%frunoff, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%vprec, isd, ied, jsd, jed)

    call safe_alloc_ptr(fluxes%sw, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lw, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%latent, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%sens, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of [R Z T-1 ~> kg m-2 s-1]
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%mask2dT(i,j)
      fluxes%lprec(i,j) = 0.0 * G%mask2dT(i,j)

      ! vprec will be set later, if it is needed for salinity restoring.
      fluxes%vprec(i,j) = 0.0

      !   Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
      fluxes%lw(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%latent(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sens(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sw(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  else ! This is the buoyancy only mode.
    do j=js,je ; do i=is,ie
      !   fluxes%buoy is the buoyancy flux into the ocean [L2 T-3 ~> m2 s-3].  A positive
      ! buoyancy flux is of the same sign as heating the ocean.
      fluxes%buoy(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      call safe_alloc_ptr(fluxes%heat_added, isd, ied, jsd, jed)
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%rho_restore * fluxes%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in [C ~> degC]) and
        ! salinity (in [S ~> ppt]) that are being restored toward.
        Temp_restore = 0.0
        Salin_restore = 0.0

        fluxes%heat_added(i,j) = (G%mask2dT(i,j) * (rhoXcp * CS%Flux_const)) * &
            (Temp_restore - sfc_state%SST(i,j))
        fluxes%vprec(i,j) = - (G%mask2dT(i,j) * (CS%rho_restore*CS%Flux_const)) * &
            ((Salin_restore - sfc_state%SSS(i,j)) / (0.5 * (Salin_restore + sfc_state%SSS(i,j))))
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%rho_restore
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [R ~> kg m-3] that is being restored toward.
        density_restore = 1030.0*US%kg_m3_to_R

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          (density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end procedure USER_buoyancy_forcing
module procedure USER_surface_forcing_init
# include "version_variable.h"
  character(len=40)  :: mdl = "user_surface_forcing" ! This module's name.
  if (associated(CS)) then
    call MOM_error(WARNING, "USER_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)

  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2*US%L_to_Z)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z/(86400.0*US%s_to_T))
  endif
  call get_param(param_file, mdl, "RESTORE_FLUX_RHO", CS%rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat "//&
                 "fluxes with RESTORE_SALINITY or RESTORE_TEMPERATURE.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(CS%Flux_const==0.0).or.(.not.CS%restorebuoy))

end procedure USER_surface_forcing_init
end submodule user_surface_forcing_s
