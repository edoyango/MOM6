submodule (regional_dyes) regional_dyes_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_dye_tracer
  character(len=40)  :: mdl = "regional_dyes" ! This module's name.
  character(len=48)  :: var_name ! The variable's name.
  character(len=48)  :: desc_name ! The variable's descriptor.
  character(len=48)  :: param_name ! The param's name suffix.
# include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers [CU ~> conc]
  integer :: isd, ied, jsd, jed, nz, m
  integer :: advect_scheme   ! Advection scheme value for this tracer
  character(len=256) :: mesg ! Advection scheme name for this tracer
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_dye_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NUM_DYE_TRACERS", CS%ntr, &
                 "The number of dye tracers in this run. Each tracer "//&
                 "should have a separate region.", default=0)
  allocate(CS%dye_source_minlon(CS%ntr), &
           CS%dye_source_maxlon(CS%ntr), &
           CS%dye_source_minlat(CS%ntr), &
           CS%dye_source_maxlat(CS%ntr), &
           CS%dye_source_mindepth(CS%ntr), &
           CS%dye_source_maxdepth(CS%ntr))
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))
  allocate(CS%id_tr_dia_diff(CS%ntr))
  CS%id_tr_dia_diff(:) = -1

  CS%dye_source_minlon(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINLON", CS%dye_source_minlon, &
                 "This is the starting longitude at which we start injecting dyes.", &
                 units="degrees_E", fail_if_missing=.true.)
               ! units=G%x_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%dye_source_minlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINLON ")

  CS%dye_source_maxlon(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXLON", CS%dye_source_maxlon, &
                 "This is the ending longitude at which we finish injecting dyes.", &
                 units="degrees_E", fail_if_missing=.true.)
               ! units=G%x_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%dye_source_maxlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXLON ")

  CS%dye_source_minlat(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINLAT", CS%dye_source_minlat, &
                 "This is the starting latitude at which we start injecting dyes.", &
                 units="degrees_N", fail_if_missing=.true.)
               ! units=G%y_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%dye_source_minlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINLAT ")

  CS%dye_source_maxlat(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXLAT", CS%dye_source_maxlat, &
                 "This is the ending latitude at which we finish injecting dyes.", &
                 units="degrees_N", fail_if_missing=.true.)
               ! units=G%y_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%dye_source_maxlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXLAT ")

  CS%dye_source_mindepth(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINDEPTH", CS%dye_source_mindepth, &
                 "This is the minimum depth at which we inject dyes.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  if (minval(CS%dye_source_mindepth(:)) < -1.e29*US%m_to_Z) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINDEPTH")

  CS%dye_source_maxdepth(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXDEPTH", CS%dye_source_maxdepth, &
                 "This is the maximum depth at which we inject dyes.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  if (minval(CS%dye_source_maxdepth(:)) < -1.e29*US%m_to_Z) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXDEPTH")

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m = 1, CS%ntr
    write(param_name(:),'(A,I3.3,A)') "DYE",m,"_TRACER_ADVECTION_SCHEME"
    call get_param(param_file, mdl, trim(param_name), mesg, &
          desc="The horizontal transport scheme for dye tracer:\n"//&
          trim(TracerAdvectionSchemeDoc)//&
          "\n Set to blank (the default) to use TRACER_ADVECTION_SCHEME.", default="")
    ! Get the integer value of the tracer scheme
    call set_tracer_advect_scheme(advect_scheme, mesg)

    write(var_name(:),'(A,I3.3)') "dye",m
    write(desc_name(:),'(A,I3.3)') "Dye Tracer ",m
    CS%tr_desc(m) = var_desc(trim(var_name), "conc", trim(desc_name), caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_dye_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit,&
                         advect_scheme=advect_scheme)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_dye_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_dye_tracer = .true.
end procedure register_dye_tracer
module procedure initialize_dye_tracer
  character(len=64)  :: var_name, longname
  real    :: dz(SZI_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real    :: z_bot    ! Height of the bottom of the layer relative to the sea surface [Z ~> m]
  real    :: z_center ! Height of the center of the layer relative to the sea surface [Z ~> m]
  integer :: i, j, k, m
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Register vertical flux diagnostic
  do m = 1, CS%ntr
    write(var_name,'(A,I3.3,A)') "dye",m,"_dia_diff"
    write(longname,'(A,I3.3,A)') "Vertical diffusive flux of dye ",m," (positive up)"
    CS%id_tr_dia_diff(m) = register_diag_field('ocean_model', trim(var_name), &
        diag%axesTi, day, trim(longname), 'conc H s-1', conversion=GV%H_to_MKS*US%s_to_T)
  enddo

  ! Establish location of source
  do j=G%jsc,G%jec
    call thickness_to_dz(h, tv, dz, j, G, GV)
    do m=1,CS%ntr ; do i=G%isc,G%iec
      ! A dye is set dependent on the center of the cell being inside the rectangular box.
      if (CS%dye_source_minlon(m) < G%geoLonT(i,j) .and. &
          CS%dye_source_maxlon(m) >= G%geoLonT(i,j) .and. &
          CS%dye_source_minlat(m) < G%geoLatT(i,j) .and. &
          CS%dye_source_maxlat(m) >= G%geoLatT(i,j) .and. &
          G%mask2dT(i,j) > 0.0 ) then
        z_bot = 0.0
        do k = 1, GV%ke
          z_bot = z_bot - dz(i,k)
          z_center = z_bot + 0.5*dz(i,k)
          if ( z_center > -CS%dye_source_maxdepth(m) .and. &
               z_center < -CS%dye_source_mindepth(m) ) then
            CS%tr(i,j,k,m) = 1.0
          endif
        enddo
      endif
    enddo ; enddo
  enddo

end procedure initialize_dye_tracer
module procedure dye_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: vert_flux ! Vertical tracer flux positive upward
  real    :: dz(SZI_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real    :: z_bot    ! Height of the bottom of the layer relative to the sea surface [Z ~> m]
  real    :: z_center ! Height of the center of the layer relative to the sea surface [Z ~> m]
  real    :: Idt      ! Inverse of timestep [T-1 ~> s-1]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  Idt = 1.0 / dt

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ; do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)

      ! Calculate net vertical flux from entrainment
      ! Net flux = upward component - downward component
      ! Upward (from below): eb(k) * tr(k+1), Downward (from above): ea(k+1) * tr(k)
      do K=2,nz ; do j=js,je ; do i=is,ie
        vert_flux(i,j,K) = (eb(i,j,k-1) * CS%tr(i,j,k,m) - ea(i,j,k) * CS%tr(i,j,k-1,m)) * Idt
      enddo ; enddo ; enddo
      do j=js,je ; do i=is,ie ; vert_flux(i,j,1) = 0.0 ; vert_flux(i,j,nz+1) = 0.0 ; enddo ; enddo

      ! Post diagnostic
      if (CS%id_tr_dia_diff(m) > 0) &
        call post_data(CS%id_tr_dia_diff(m), vert_flux, CS%diag)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)

      ! Calculate net vertical flux from entrainment
      ! Net flux = upward component - downward component
      ! Upward (from below): eb(k) * tr(k+1), Downward (from above): ea(k+1) * tr(k)
      do K=2,nz ; do j=js,je ; do i=is,ie
        vert_flux(i,j,K) = (eb(i,j,k-1) * CS%tr(i,j,k,m) - ea(i,j,k) * CS%tr(i,j,k-1,m)) * Idt
      enddo ; enddo ; enddo
      do j=js,je ; do i=is,ie ; vert_flux(i,j,1) = 0.0 ; vert_flux(i,j,nz+1) = 0.0 ; enddo ; enddo

      ! Post diagnostic
      if (CS%id_tr_dia_diff(m) > 0) &
        call post_data(CS%id_tr_dia_diff(m), vert_flux, CS%diag)
    enddo
  endif

  do j=js,je
    call thickness_to_dz(h_new, tv, dz, j, G, GV)
    do m=1,CS%ntr ; do i=is,ie
      ! A dye is set dependent on the center of the cell being inside the rectangular box.
      if (CS%dye_source_minlon(m) < G%geoLonT(i,j) .and. &
          CS%dye_source_maxlon(m) >= G%geoLonT(i,j) .and. &
          CS%dye_source_minlat(m) < G%geoLatT(i,j) .and. &
          CS%dye_source_maxlat(m) >= G%geoLatT(i,j) .and. &
          G%mask2dT(i,j) > 0.0 ) then
        z_bot = 0.0
        do k=1,nz
          z_bot = z_bot - dz(i,k)
          z_center = z_bot + 0.5*dz(i,k)
          if ( z_center > -CS%dye_source_maxdepth(m) .and. &
               z_center < -CS%dye_source_mindepth(m) ) then
            CS%tr(i,j,k,m) = 1.0
          endif
        enddo
      endif
    enddo ; enddo
  enddo

end procedure dye_tracer_column_physics
module procedure dye_stock
  integer :: m
  dye_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="dye_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  dye_stock = CS%ntr

end procedure dye_stock
module procedure dye_tracer_surface_state
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

end procedure dye_tracer_surface_state
module procedure regional_dyes_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure regional_dyes_end
end submodule regional_dyes_s
