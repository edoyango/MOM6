submodule (advection_test_tracer) advection_test_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_advection_test_tracer
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "advection_test_tracer" ! This module's name.
  character(len=200) :: inputdir   ! The directory where the input file can be found
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to a tracer array [conc]
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_advection_test_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "ADVECTION_TEST_X_ORIGIN", CS%x_origin, &
        "The x-coordinate of the center of the test-functions.", units=G%x_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_Y_ORIGIN", CS%y_origin, &
        "The y-coordinate of the center of the test-functions.", units=G%y_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_X_WIDTH", CS%x_width, &
        "The x-width of the test-functions.", units=G%x_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_Y_WIDTH", CS%y_width, &
        "The y-width of the test-functions.", units=G%y_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the tracers, or blank to initialize "//&
                 "them internally.", default=" ")

  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/ADVECTION_TEST_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)


  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR), source=0.0)

  do m=1,NTR
    write(name,'("tr",I0)') m
    write(longname,'("Concentration of Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif


    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_advection_test_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_advection_test_tracer = .true.
end procedure register_advection_test_tracer
module procedure initialize_advection_test_tracer
  character(len=16) :: name ! A variable's name in a NetCDF file.
  real :: locx, locy        ! x- and y- positions relative to the center of the tracer patch
  real :: h_neglect         ! A thickness that is so small it is usually lost
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%diag => diag
  CS%ntr = NTR
  do m=1,NTR
    call query_vardesc(CS%tr_desc(m), name=name, &
                       caller="initialize_advection_test_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
      k=1 ! Square wave
      do j=js,je ; do i=is,ie
        if (abs(G%geoLonT(i,j)-CS%x_origin)<0.5*CS%x_width .and. &
            abs(G%geoLatT(i,j)-CS%y_origin)<0.5*CS%y_width) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=2 ! Triangle wave
      do j=js,je ; do i=is,ie
        locx = abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy = abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        CS%tr(i,j,k,m) = max(0.0, 1.0-locx)*max(0.0, 1.0-locy)
      enddo ; enddo
      k=3 ! Cosine bell
      do j=js,je ; do i=is,ie
        locx = min(1.0, abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width) * (acos(0.0)*2.)
        locy = min(1.0, abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width) * (acos(0.0)*2.)
        CS%tr(i,j,k,m) = (1.0+cos(locx))*(1.0+cos(locy))*0.25
      enddo ; enddo
      k=4 ! Cylinder
      do j=js,je ; do i=is,ie
        locx = abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy = abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if ((locx**2) + (locy**2) <= 1.0) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=5 ! Cut cylinder
      do j=js,je ; do i=is,ie
        locx = (G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy = (G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if ((locx**2) + (locy**2) <= 1.0) CS%tr(i,j,k,m) = 1.0
        if (locx>0.0 .and. abs(locy)<0.2) CS%tr(i,j,k,m) = 0.0
      enddo ; enddo

      call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif ! restart
  enddo


end procedure initialize_advection_test_tracer
module procedure advection_test_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

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
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end procedure advection_test_tracer_column_physics
module procedure advection_test_tracer_surface_state
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

end procedure advection_test_tracer_surface_state
module procedure advection_test_stock
  integer :: is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  advection_test_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="advection_test_stock")
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  advection_test_stock = CS%ntr

end procedure advection_test_stock
module procedure advection_test_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure advection_test_tracer_end
end submodule advection_test_tracer_s
