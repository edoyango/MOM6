submodule (oil_tracer) oil_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_oil_tracer
  character(len=40)  :: mdl = "oil_tracer" ! This module's name.
# include "version_variable.h"
  real, dimension(NTR_MAX) :: oil_decay_days  !< Decay time scale of oil [days]
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=3)   :: name_tag ! String for creating identifying oils
  character(len=48) :: flux_units ! The units for tracer fluxes, here
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [kg m-3]
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_oil_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "OIL_IC_FILE", CS%IC_file, &
                 "The file in which the oil tracer initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/OIL_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "OIL_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, OIL_IC_FILE is in depth space, not layer space", &
                 default=.false.)

  call get_param(param_file, mdl, "OIL_MAY_REINIT", CS%oil_may_reinit, &
                 "If true, oil tracers may go through the initialization "//&
                 "code if they are not found in the restart files. "//&
                 "Otherwise it is a fatal error if the oil tracers are not "//&
                 "found in the restart files of a restarted run.", &
                 default=.false.)
  call get_param(param_file, mdl, "OIL_SOURCE_LONGITUDE", CS%oil_source_longitude, &
                 "The geographic longitude of the oil source.", units="degrees_E", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "OIL_SOURCE_LATITUDE", CS%oil_source_latitude, &
                 "The geographic latitude of the oil source.", units="degrees_N", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "OIL_SOURCE_LAYER", CS%oil_source_k, &
                 "The layer into which the oil is introduced, or a "//&
                 "negative number for a vertically uniform source, "//&
                 "or 0 not to use this tracer.", units="Layer", default=0)
  call get_param(param_file, mdl, "OIL_SOURCE_RATE", CS%oil_source_rate, &
                 "The rate of oil injection.", &
                 units="kg s-1", scale=US%T_to_s, default=1.0)
  call get_param(param_file, mdl, "OIL_DECAY_DAYS", oil_decay_days, &
                 "The decay timescale in days (if positive), or no decay "//&
                 "if 0, or use the temperature dependent decay rate of "//&
                 "Adcroft et al. (GRL, 2010) if negative.", units="days", &
                 default=0.0)
  call get_param(param_file, mdl, "OIL_DATED_START_YEAR", CS%oil_start_year, &
                 "The time at which the oil source starts", units="years", &
                 default=0.0)
  call get_param(param_file, mdl, "OIL_DATED_END_YEAR", CS%oil_end_year, &
                 "The time at which the oil source ends", units="years", &
                 default=1.0e99)

  CS%ntr = 0
  CS%oil_decay_rate(:) = 0.
  do m=1,NTR_MAX
    if (CS%oil_source_k(m)/=0) then
      write(name_tag(1:3),'("_",I2.2)') m
      CS%ntr = CS%ntr + 1
      CS%tr_desc(m) = var_desc("oil"//trim(name_tag), "kg m-3", "Oil Tracer", caller=mdl)
      CS%IC_val(m) = 0.0
      if (oil_decay_days(m) > 0.) then
        CS%oil_decay_rate(m) = 1. / (86400.0*US%s_to_T * oil_decay_days(m))
      elseif (oil_decay_days(m) < 0.) then
        CS%oil_decay_rate(m) = -1.
      endif
    endif
  enddo
  call log_param(param_file, mdl, "OIL_DECAY_RATE", CS%oil_decay_rate(1:CS%ntr), &
                 units="s-1", unscale=US%s_to_T)

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "kg s-1"
  else ; flux_units = "kg m-3 kg s-1" ; endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, caller="register_oil_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=CS%tr_desc(m), &
                         registry_diags=.true., flux_units=flux_units, restart_CS=restart_CS, &
                         mandatory=.not.CS%oil_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_oil_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_oil_tracer = .true.

end procedure register_oil_tracer
module procedure initialize_oil_tracer
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Establish location of source
  do j=G%jsdB+1,G%jed ; do i=G%isdB+1,G%ied
    ! This test for i,j index is specific to a lat/lon (non-rotated grid).
    ! and needs to be generalized to work properly on the tri-polar grid.
    if (CS%oil_source_longitude<G%geoLonBu(I,J) .and. &
        CS%oil_source_longitude>=G%geoLonBu(I-1,J) .and. &
        CS%oil_source_latitude<G%geoLatBu(I,J) .and. &
        CS%oil_source_latitude>=G%geoLatBu(I,J-1) ) then
      CS%oil_source_i=i
      CS%oil_source_j=j
    endif
  enddo ; enddo

  CS%Time => day
  CS%diag => diag

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, caller="initialize_oil_tracer")
    if ((.not.restart) .or. (CS%oil_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then

      if (len_trim(CS%IC_file) > 0) then
  !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%IC_file, G%Domain)) &
          call MOM_error(FATAL, "initialize_oil_tracer: Unable to open "//CS%IC_file)

        if (CS%Z_IC_file) then
          OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name, &
                             G, GV, US, -1e34, 0.0) ! CS%land_val(m))
          if (.not.OK) then
            OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, &
                     trim(name), G, GV, US, -1e34, 0.0) ! CS%land_val(m))
            if (.not.OK) call MOM_error(FATAL,"initialize_oil_tracer: "//&
                    "Unable to read "//trim(name)//" from "//&
                    trim(CS%IC_file)//".")
          endif
        else
          call MOM_read_data(CS%IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
        endif
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) < 0.5) then
            CS%tr(i,j,k,m) = CS%land_val(m)
          else
            CS%tr(i,j,k,m) = CS%IC_val(m)
          endif
        enddo ; enddo ; enddo
      endif
      call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif ! restart
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! Put something here...
  endif

end procedure initialize_oil_tracer
module procedure oil_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: Isecs_per_year = 1.0 / (365.0*86400.0) ! Conversion factor from seconds to year [year s-1]
  real :: vol_scale ! A conversion factor for volumes into m3 [m3 H-1 L-2 ~> 1 or m3 kg-1]
  real :: year      ! Time in fractional years [years]
  real :: h_total   ! A running sum of thicknesses [H ~> m or kg m-2]
  real :: decay_timescale ! Chemical decay timescale for oil [T ~> s]
  real :: ldecay    ! Chemical decay rate of oil [T-1 ~> s-1]
  integer :: i, j, k, is, ie, js, je, nz, m, k_max
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  year = time_type_to_real(CS%Time) * Isecs_per_year

  ! Decay tracer (limit decay rate to 1/dt - just in case)
  do m=2,CS%ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      !CS%tr(i,j,k,m) = CS%tr(i,j,k,m) - dt*CS%oil_decay_rate(m)*CS%tr(i,j,k,m) ! Simple
      !CS%tr(i,j,k,m) = CS%tr(i,j,k,m) - min(dt*CS%oil_decay_rate(m),1.)*CS%tr(i,j,k,m) ! Safer
      if (CS%oil_decay_rate(m)>0.) then
        CS%tr(i,j,k,m) = G%mask2dT(i,j)*max(1. - dt*CS%oil_decay_rate(m),0.)*CS%tr(i,j,k,m) ! Safest
      elseif (CS%oil_decay_rate(m)<0.) then
        decay_timescale = (12.0 * (3.0**(-(tv%T(i,j,k)-20.0*US%degC_to_C)/10.0*US%degC_to_C))) * &
                          (86400.0*US%s_to_T) ! Timescale [T ~> s]
        ldecay = 1. / decay_timescale ! Rate [T-1 ~> s-1]
        CS%tr(i,j,k,m) = G%mask2dT(i,j)*max(1. - dt*ldecay,0.)*CS%tr(i,j,k,m)
      endif
    enddo ; enddo ; enddo
  enddo

  ! Add oil at the source location
  if (year>=CS%oil_start_year .and. year<=CS%oil_end_year .and. &
      CS%oil_source_i>-999 .and. CS%oil_source_j>-999) then
    i = CS%oil_source_i ; j = CS%oil_source_j
    k_max = nz ; h_total = 0.
    vol_scale = GV%H_to_m * US%L_to_m**2
    do k=nz, 2, -1
      h_total = h_total + h_new(i,j,k)
      if (h_total < 10.*GV%m_to_H) k_max=k-1 ! Find bottom most interface that is 10 m above bottom
    enddo
    do m=1,CS%ntr
      k = CS%oil_source_k(m)
      if (k>0) then
        k = min(k,k_max) ! Only insert k or first layer with interface 10 m above bottom
        CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + CS%oil_source_rate*dt / &
                (vol_scale * (h_new(i,j,k)+GV%H_subroundoff) * G%areaT(i,j) )
      elseif (k<0) then
        h_total = GV%H_subroundoff
        do k=1, nz
          h_total = h_total + h_new(i,j,k)
        enddo
        do k=1, nz
          CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + CS%oil_source_rate*dt / (vol_scale * h_total * G%areaT(i,j) )
        enddo
      endif
    enddo
  endif

end procedure oil_tracer_column_physics
module procedure oil_stock
  integer :: m
  oil_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="oil_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  oil_stock = CS%ntr

end procedure oil_stock
module procedure oil_tracer_surface_state
  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/), turns=G%HI%turns)
    enddo
  endif

end procedure oil_tracer_surface_state
module procedure oil_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure oil_tracer_end
end submodule oil_tracer_s
