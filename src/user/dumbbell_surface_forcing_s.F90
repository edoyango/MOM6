submodule (dumbbell_surface_forcing) dumbbell_surface_forcing_s
  implicit none
contains
module procedure dumbbell_buoyancy_forcing
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  ! Allocate and zero out the forcing arrays, as necessary.
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

      ! Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
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

  if (CS%use_temperature .and. CS%restorebuoy) then
    do j=js,je ; do i=is,ie
      if (CS%forcing_mask(i,j)>0.) then
        fluxes%vprec(i,j) = - (G%mask2dT(i,j) * CS%Flux_const) * &
                ((CS%S_restore(i,j) - sfc_state%SSS(i,j)) /  (0.5 * (CS%S_restore(i,j) + sfc_state%SSS(i,j))))

      endif
    enddo ; enddo
  endif

end procedure dumbbell_buoyancy_forcing
module procedure dumbbell_dynamic_forcing
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  integer :: idays, isecs
  real :: deg_rad  ! A conversion factor from degrees to radians [nondim]
  real :: rdays    ! The elapsed time [days]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  deg_rad = atan(1.0)*4.0/180.

  call get_time(day,isecs,idays)
  rdays = real(idays) + real(isecs)/8.64e4
  ! This could be:  rdays = time_type_to_real(day)/8.64e4

  ! Allocate and zero out the forcing arrays, as necessary.
  call safe_alloc_ptr(fluxes%p_surf, isd, ied, jsd, jed)
  call safe_alloc_ptr(fluxes%p_surf_full, isd, ied, jsd, jed)

  do j=js,je ; do i=is,ie
    fluxes%p_surf(i,j) = CS%forcing_mask(i,j)* CS%slp_amplitude * &
                         G%mask2dT(i,j) * sin(deg_rad*(rdays/CS%slp_period))
    fluxes%p_surf_full(i,j) = CS%forcing_mask(i,j) * CS%slp_amplitude * &
                         G%mask2dT(i,j) * sin(deg_rad*(rdays/CS%slp_period))
  enddo ; enddo

end procedure dumbbell_dynamic_forcing
module procedure dumbbell_surface_forcing_init
  real :: S_surf  ! Initial surface salinity [S ~> ppt]
  real :: S_range ! Range of the initial vertical distribution of salinity [S ~> ppt]
  real :: x       ! Latitude normalized by the domain size [nondim]
  real :: Rho0          ! The density used in the Boussinesq approximation [R ~> kg m-3]
  real :: rho_restore   ! The density that is used to convert piston velocities into salt
  integer :: i, j
  logical :: dbrotate    ! If true, rotate the domain.
# include "version_variable.h"
  character(len=40)  :: mdl = "dumbbell_surface_forcing" ! This module's name.
  if (associated(CS)) then
    call MOM_error(WARNING, "dumbbell_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state variables.", default=.true.)

  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)
  call get_param(param_file, mdl, "RHO_0", Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DUMBBELL_SLP_AMP", CS%slp_amplitude, &
                 "Amplitude of SLP forcing in reservoirs.", &
                 units="Pa", default=10000.0, scale=US%Pa_to_RL2_T2)
  call get_param(param_file, mdl, "DUMBBELL_SLP_PERIOD", CS%slp_period, &
                 "Periodicity of SLP forcing in reservoirs.", &
                 units="days", default=1.0)
  call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.',&
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl,"INITIAL_SSS", S_surf, &
                 "Initial surface salinity", &
                 units="ppt", default=34.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl,"INITIAL_S_RANGE", S_range, &
                 "Initial salinity range (bottom - surface)", &
                 units="ppt", default=2., scale=US%ppt_to_S, do_not_log=.true.)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default=.false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s)
    call get_param(param_file, mdl, "RESTORE_FLUX_RHO", rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat "//&
                 "fluxes with RESTORE_SALINITY or RESTORE_TEMPERATURE.", &
                 units="kg m-3", default=Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(CS%Flux_const==0.0))
    ! Convert FLUXCONST from m day-1 to m s-1 and Flux_const to [R Z T-1 ~> kg m-2 s-1]
    CS%Flux_const = rho_restore * (CS%Flux_const / 86400.0)


    allocate(CS%forcing_mask(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    allocate(CS%S_restore(G%isd:G%ied, G%jsd:G%jed))

    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        ! Compute normalized zonal coordinates (x,y=0 at center of domain)
!       x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
!       y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
        if (dbrotate) then
          ! This is really y in the rotated case
          x = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
        else
          x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
        endif
        CS%forcing_mask(i,j)=0
        CS%S_restore(i,j) = S_surf
        if ((x>0.25)) then
          CS%forcing_mask(i,j) = 1
          CS%S_restore(i,j) = S_surf + S_range
        elseif ((x<-0.25)) then
          CS%forcing_mask(i,j) = 1
          CS%S_restore(i,j) = S_surf - S_range
        endif
      enddo
    enddo
  endif

end procedure dumbbell_surface_forcing_init
end submodule dumbbell_surface_forcing_s
