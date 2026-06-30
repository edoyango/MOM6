submodule (boundary_impulse_tracer) boundary_impulse_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_boundary_impulse_tracer
  character(len=40)  :: mdl = "boundary_impulse_tracer" ! This module's name.
  character(len=48)  :: var_name ! The variable's name.
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
# include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [CU ~> conc]
  real, pointer :: rem_time_ptr => NULL() ! The ramaining injection time [T ~> s]
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_boundary_impulse_tracer called with an "// &
                           "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "IMPULSE_SOURCE_TIME", CS%remaining_source_time, &
                 "Length of time for the boundary tracer to be injected "//&
                 "into the mixed layer. After this time has elapsed, the "//&
                 "surface becomes a sink for the boundary impulse tracer.", &
                 units="s", default=31536000.0, scale=US%s_to_T)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)
  CS%ntr = NTR_MAX
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  CS%nkml = max(GV%nkml,1)

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    CS%tr_desc(m) = var_desc(trim("boundary_impulse"), "kg kg-1", &
        "Boundary impulse tracer", caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, caller="register_boundary_impulse_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=CS%tr_desc(m), &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_boundary_impulse_tracer")
  enddo
  ! Register remaining source time as a restart field
  rem_time_ptr => CS%remaining_source_time
  call register_restart_field(rem_time_ptr, "bir_remain_time", &
                              .not.CS%tracers_may_reinit, restart_CS, &
                              "Remaining time to apply BIR source", "s", conversion=US%T_to_s)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_boundary_impulse_tracer = .true.

end procedure register_boundary_impulse_tracer
module procedure initialize_boundary_impulse_tracer
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  name = "boundary_impulse"

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, caller="initialize_boundary_impulse_tracer")
    if ((.not.restart) .or. (.not. query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      do k=1,CS%nkml ; do j=jsd,jed ; do i=isd,ied
        CS%tr(i,j,k,m) = 1.0
      enddo ; enddo ; enddo
      call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

end procedure initialize_boundary_impulse_tracer
module procedure boundary_impulse_tracer_column_physics
  integer :: i, j, k, is, ie, js, je, nz, m
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  ! This uses applyTracerBoundaryFluxesInOut, usually in ALE mode
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,1), dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,1), G, GV)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,1), G, GV)
  endif

  ! Set surface conditions
  do m=1,1
    if (CS%remaining_source_time > 0.0) then
      do k=1,CS%nkml ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 1.0
      enddo ; enddo ; enddo
      CS%remaining_source_time = CS%remaining_source_time-dt
    else
      do k=1,CS%nkml ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
    endif

  enddo

end procedure boundary_impulse_tracer_column_physics
module procedure boundary_impulse_stock
  integer :: m
  boundary_impulse_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,1
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="boundary_impulse_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo

  boundary_impulse_stock = CS%ntr

end procedure boundary_impulse_stock
module procedure boundary_impulse_tracer_surface_state
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

end procedure boundary_impulse_tracer_surface_state
module procedure boundary_impulse_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure boundary_impulse_tracer_end
end submodule boundary_impulse_tracer_s
