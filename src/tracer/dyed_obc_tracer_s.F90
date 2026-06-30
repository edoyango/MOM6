submodule (dyed_obc_tracer) dyed_obc_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_dyed_obc_tracer
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "dyed_obc_tracer" ! This module's name.
  character(len=200) :: inputdir
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [conc]
  integer :: isd, ied, jsd, jed, nz, m
  integer :: n_dye           ! Number of regionsl dye tracers
  integer :: advect_scheme   ! Advection scheme value for this tracer
  character(len=256) :: mesg ! Advection scheme name for this tracer
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "dyed_obc_register_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NUM_DYED_TRACERS", CS%ntr, &
                 "The number of dyed_obc tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.  "//&
                 "If not present, use NUM_DYE_TRACERS.", default=-1)
  if (CS%ntr == -1) then
    !for backward compatibility
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", CS%ntr, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate boundary segment.", default=0)
    n_dye = 0
  else
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", n_dye, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate region.", default=0, do_not_log=.true.)
  endif
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))

  call get_param(param_file, mdl, "dyed_obc_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the dyed_obc tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/dyed_obc_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif

  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)

  call get_param(param_file, mdl, "DYED_TRACER_ADVECTION_SCHEME", mesg, &
        desc="The horizontal transport scheme for dyed_obc tracers:\n"//&
        trim(TracerAdvectionSchemeDoc)//&
        "\n Set to blank (the default) to use TRACER_ADVECTION_SCHEME.", default="")

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m=1,CS%ntr
    write(name,'("dye_",I2.2)') m+n_dye  !after regional dye tracers
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Get the integer value of the tracer scheme
    call set_tracer_advect_scheme(advect_scheme, mesg)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit, &
                         advect_scheme=advect_scheme)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_dyed_obc_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_dyed_obc_tracer = .true.
end procedure register_dyed_obc_tracer
module procedure initialize_dyed_obc_tracer
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  real :: h_neglect         ! A thickness that is so small it is usually lost
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%Time => day
  CS%diag => diag

  do m=1,CS%ntr
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      if (len_trim(CS%tracer_IC_file) >= 1) then
        !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
          call MOM_error(FATAL, "dyed_obc_initialize_tracer: Unable to open "// &
                          CS%tracer_IC_file)
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_dyed_obc_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 0.0
        enddo ; enddo ; enddo
      endif
      call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif ! restart
  enddo ! Tracer loop

end procedure initialize_dyed_obc_tracer
module procedure dyed_obc_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, m
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
      if (nz > 1) call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%ntr
      if (nz > 1) call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end procedure dyed_obc_tracer_column_physics
module procedure dyed_obc_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)

    deallocate(CS)
  endif
end procedure dyed_obc_tracer_end
end submodule dyed_obc_tracer_s
