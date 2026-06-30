submodule (ideal_age_example) ideal_age_example_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_ideal_age_tracer
# include "version_variable.h"
  character(len=40)  :: mdl = "ideal_age_example" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [years]
  logical :: do_ideal_age, do_vintage, do_ideal_age_dated, do_BL_residence
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_ideal_age_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DO_IDEAL_AGE", do_ideal_age, &
                 "If true, use an ideal age tracer that is set to 0 age "//&
                 "in the boundary layer and ages at unit rate in the interior.", &
                 default=.true.)
  call get_param(param_file, mdl, "DO_IDEAL_VINTAGE", do_vintage, &
                 "If true, use an ideal vintage tracer that is set to an "//&
                 "exponentially increasing value in the boundary layer and "//&
                 "is conserved thereafter.", default=.false.)
  call get_param(param_file, mdl, "DO_IDEAL_AGE_DATED", do_ideal_age_dated, &
                 "If true, use an ideal age tracer that is everywhere 0 "//&
                 "before IDEAL_AGE_DATED_START_YEAR, but the behaves like "//&
                 "the standard ideal age tracer - i.e. is set to 0 age in "//&
                 "the boundary layer and ages at unit rate in the interior.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_BL_RESIDENCE", do_BL_residence, &
                 "If true, use a residence tracer that is set to 0 age "//&
                 "in the interior and ages at unit rate in the boundary layer.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_REAL_BL_DEPTH", CS%use_real_BL_depth, &
                 "If true, the ideal age tracers will use the boundary layer "//&
                 "depth diagnosed from the BL or bulkmixedlayer scheme.", &
                 default=.false.)
  call get_param(param_file, mdl, "AGE_IC_FILE", CS%IC_file, &
                 "The file in which the age-tracer initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/AGE_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "AGE_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, AGE_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)

  CS%ntr = 0
  if (do_ideal_age) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age", "yr", "Ideal Age Tracer", cmor_field_name="agessc", caller=mdl)
    CS%tracer_ages(m) = .true. ; CS%growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
  endif

  if (do_vintage) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("vintage", "yr", "Exponential Vintage Tracer", &
                            caller=mdl)
    CS%tracer_ages(m) = .false. ; CS%growth_rate(m) = 1.0/30.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 1e-20 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mdl, "IDEAL_VINTAGE_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the ideal vintage tracer starts.", &
                 units="years", default=0.0)
  endif

  if (do_ideal_age_dated) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age_dated","yr","Ideal Age Tracer with a Start Date",&
                            caller=mdl)
    CS%tracer_ages(m) = .true. ; CS%growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mdl, "IDEAL_AGE_DATED_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the dated ideal age tracer starts.", &
                 units="years", default=0.0)
  endif

  CS%BL_residence_num = 0
  if (do_BL_residence) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr ; CS%BL_residence_num = CS%ntr
    CS%tr_desc(m) = var_desc("BL_age", "yr", "BL Residence Time Tracer", caller=mdl)
    CS%tracer_ages(m) = .true. ; CS%growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
  endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_ideal_age_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=CS%tr_desc(m), &
                         registry_diags=.true., restart_CS=restart_CS, &
                         mandatory=.not.CS%tracers_may_reinit, &
                         flux_scale=GV%H_to_m)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_ideal_age_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_ideal_age_tracer = .true.
end procedure register_ideal_age_tracer
module procedure initialize_ideal_age_tracer
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  logical :: use_real_BL_depth
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  CS%nkbl = max(GV%nkml,1)

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, &
                       caller="initialize_ideal_age_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then

      if (len_trim(CS%IC_file) > 0) then
  !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%IC_file, G%Domain)) &
          call MOM_error(FATAL, "initialize_ideal_age_tracer: "// &
                                 "Unable to open "//CS%IC_file)

        if (CS%Z_IC_file) then
          OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name,&
                             G, GV, US, -1e34, 0.0) ! CS%land_val(m))
          if (.not.OK) then
            OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, &
                     trim(name), G, GV, US, -1e34, 0.0) ! CS%land_val(m))
            if (.not.OK) call MOM_error(FATAL,"initialize_ideal_age_tracer: "//&
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
  ! Steal from updated DOME in the fullness of time.
  endif

end procedure initialize_ideal_age_tracer
module procedure ideal_age_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G)) :: BL_layers ! Stores number of layers in boundary layer [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: young_val       ! The "young" value for the tracers [years] or other units
  real :: Isecs_per_year  ! The inverse of the amount of time in a year [T-1 ~> s-1]
  real :: year            ! The time in years [years]
  real :: layer_frac      ! The fraction of the current layer that is within the mixed layer [nondim]
  integer :: i, j, k, is, ie, js, je, nz, m, nk
  character(len=255) :: msg
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (CS%use_real_BL_depth .and. .not. present(Hbl)) then
    call MOM_error(FATAL, "Attempting to use real boundary layer depth for ideal age tracers, " &
         // "but no valid boundary layer scheme was found")
  endif

  if (CS%use_real_BL_depth .and. present(Hbl)) then
    call count_BL_layers(G, GV, h_old, Hbl, BL_layers)
  endif

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

  Isecs_per_year = 1.0 / (365.0*86400.0*US%s_to_T)
  !   Set the surface value of tracer 1 to increase exponentially
  ! with a 30 year time scale.
  year = time_to_real(CS%Time, scale=US%s_to_T) * Isecs_per_year

  do m=1,CS%ntr

    if (CS%growth_rate(m) == 0.0) then
      young_val = CS%young_val(m)
    else
      young_val = CS%young_val(m) * &
          exp((year-CS%tracer_start_year(m)) * CS%growth_rate(m))
    endif

    if (m == CS%BL_residence_num) then

      if (CS%use_real_BL_depth) then
        do j=js,je ; do i=is,ie
          nk = floor(BL_layers(i,j))

          do k=1,nk
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo

          k = MIN(nk+1,nz)

          if (G%mask2dT(i,j) > 0.0) then
            layer_frac = BL_layers(i,j)-nk
            CS%tr(i,j,k,m) = layer_frac * (CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt &
                             *Isecs_per_year) + (1.-layer_frac) * young_val
          else
            CS%tr(i,j,k,m) = CS%land_val(m)
          endif


          do k=nk+2,nz
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = young_val
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo
       enddo ; enddo

      else  ! use real BL depth
        do j=js,je ; do i=is,ie
          do k=1,CS%nkbl
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo

          do k=CS%nkbl+1,nz
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = young_val
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo
        enddo ; enddo

     endif ! use real BL depth

    else ! if BL residence tracer

      if (CS%use_real_BL_depth) then
        do j=js,je ; do i=is,ie
          nk = floor(BL_layers(i,j))
          do k=1,nk
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = young_val
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo

          k = MIN(nk+1,nz)
          if (G%mask2dT(i,j) > 0.0) then
            layer_frac = BL_layers(i,j)-nk
            CS%tr(i,j,k,m) = (1.-layer_frac) * (CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt &
                             *Isecs_per_year) + layer_frac * young_val
          else
            CS%tr(i,j,k,m) = CS%land_val(m)
          endif

          do k=nk+2,nz
            if (G%mask2dT(i,j) > 0.0) then
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
            else
              CS%tr(i,j,k,m) = CS%land_val(m)
            endif
          enddo
        enddo ; enddo

      else ! use real BL depth
        do k=1,CS%nkbl ; do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) > 0.0) then
            CS%tr(i,j,k,m) = young_val
          else
            CS%tr(i,j,k,m) = CS%land_val(m)
          endif
        enddo ; enddo ; enddo

        if (CS%tracer_ages(m) .and. (year>=CS%tracer_start_year(m))) then
          !$OMP parallel do default(none) shared(is,ie,js,je,CS,nz,G,dt,Isecs_per_year,m)
          do k=CS%nkbl+1,nz ; do j=js,je ; do i=is,ie
            CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
          enddo ; enddo ; enddo
        endif


      endif ! if use real BL depth
    endif ! if BL residence tracer

  enddo ! loop over all tracers

end procedure ideal_age_tracer_column_physics
module procedure ideal_age_stock
  integer :: m
  ideal_age_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="ideal_age_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  ideal_age_stock = CS%ntr

end procedure ideal_age_stock
module procedure ideal_age_tracer_surface_state
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

end procedure ideal_age_tracer_surface_state
module procedure ideal_age_example_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure ideal_age_example_end
module procedure count_BL_layers
  real :: current_depth  ! Distance from the free surface [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, m, nk
  character(len=255) :: msg
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  BL_layers(:,:) = 0.
  do j=js,je
    do i=is,ie
      current_depth = 0.
      do k=1,nz
        current_depth = current_depth + h(i,j,k)
        if (Hbl(i,j) <= current_depth) then
          BL_layers(i,j) = BL_layers(i,j) + (1.0 - (current_depth - Hbl(i,j)) / h(i,j,k))
          exit
        else
          BL_layers(i,j) = BL_layers(i,j) + 1.0
        endif
      enddo
    enddo
  enddo

end procedure count_BL_layers
end submodule ideal_age_example_s
