submodule (MESO_surface_forcing) MESO_surface_forcing_s
  implicit none
contains
module procedure MESO_buoyancy_forcing
  real :: density_restore  ! The potential density that is being restored toward [R ~> kg m-3].
  real :: rhoXcp ! The mean density times the heat capacity [Q R C-1 ~> J m-3 degC-1].
  real :: buoy_rest_const  ! A constant relating density anomalies to the
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.

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
    call safe_alloc_ptr(fluxes%heat_content_lprec, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.
  if (CS%restorebuoy .and. first_call) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
    call safe_alloc_ptr(CS%T_Restore, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%S_Restore, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%Heat, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%PmE, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%Solar, isd, ied, jsd, jed)

    call MOM_read_data(trim(CS%inputdir)//trim(CS%SSTrestore_file), "SST", &
             CS%T_Restore(:,:), G%Domain, scale=US%degC_to_C)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%salinityrestore_file), "SAL", &
             CS%S_Restore(:,:), G%Domain, scale=US%ppt_to_S)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%heating_file), "Heat", &
             CS%Heat(:,:), G%Domain, scale=US%W_m2_to_QRZ_T)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%PmE_file), "PmE", &
             CS%PmE(:,:), G%Domain, scale=US%m_to_Z*US%T_to_s)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%Solar_file), "NET_SOL", &
             CS%Solar(:,:), G%Domain, scale=US%W_m2_to_QRZ_T)
    first_call = .false.
  endif

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of [R Z T-1 ~> kg m-2 s-1]
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%mask2dT(i,j)
      fluxes%lprec(i,j) =  CS%PmE(i,j) * CS%Rho0 * G%mask2dT(i,j)

      ! vprec will be set later, if it is needed for salinity restoring.
      fluxes%vprec(i,j) = 0.0

      !   Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
      fluxes%lw(i,j)     = 0.0 * G%mask2dT(i,j)
      fluxes%latent(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sens(i,j)   = CS%Heat(i,j) * G%mask2dT(i,j)
      fluxes%sw(i,j)     = CS%Solar(i,j) * G%mask2dT(i,j)
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
!      call MOM_error(FATAL, "MESO_buoyancy_surface_forcing: " // &
!        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%rho_restore * fluxes%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in degC) and
        ! salinity (in ppt or PSU) that are being restored toward.
        if (G%mask2dT(i,j) > 0.0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%rho_restore * CS%Flux_const) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "MESO_buoyancy_surface_forcing: " // &
        "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%rho_restore
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [R ~> kg m-3] that is being restored toward.
        density_restore = 1030.0 * US%kg_m3_to_R

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                           (density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end procedure MESO_buoyancy_forcing
module procedure MESO_surface_forcing_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MESO_surface_forcing" ! This module's name.
  if (associated(CS)) then
    call MOM_error(WARNING, "MESO_surface_forcing_init called with an associated "// &
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
                 "The background gustiness in the winds.", units="Pa", default=0.0, &
                 scale=US%Pa_to_RLZ_T2)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)

  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z/(86400.0*US%s_to_T))

    call get_param(param_file, mdl, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in "//&
                 "variable TEMP.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to "//&
                 "restore in variable SALT.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%heating_file, &
                 "The file with the non-shortwave heat flux in "//&
                 "variable Heat.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "PRECIP_FILE", CS%PmE_file, &
                 "The file with the net precipiation minus evaporation "//&
                 "in variable PmE.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SHORTWAVE_FILE", CS%Solar_file, &
                 "The file with the shortwave heat flux in "//&
                 "variable NET_SOL.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
    CS%inputdir = slasher(CS%inputdir)
    call get_param(param_file, mdl, "RESTORE_FLUX_RHO", CS%rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat "//&
                 "fluxes with RESTORE_SALINITY or RESTORE_TEMPERATURE.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(CS%Flux_const==0.0).or.(.not.CS%restorebuoy))
  endif

end procedure MESO_surface_forcing_init
end submodule MESO_surface_forcing_s
