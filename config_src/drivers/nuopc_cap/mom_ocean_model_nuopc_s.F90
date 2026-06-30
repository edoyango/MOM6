submodule (MOM_ocean_model_nuopc) MOM_ocean_model_nuopc_s
#include <MOM_memory.h>
  implicit none
contains
module procedure ocean_model_init
  real :: Rho0        ! The Boussinesq ocean density, in kg m-3.
  real :: G_Earth     ! The gravitational acceleration in m s-2.
  real :: HFrz        !< If HFrz > 0 (m), melt potential will be computed.
  logical :: use_melt_pot !< If true, allocate melt_potential array
#include "version_variable.h"
  character(len=40)  :: mdl = "ocean_model_init"  ! This module's name.
  character(len=48)  :: stagger
  integer :: secs, days
  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  logical :: use_temperature
  call callTree_enter("ocean_model_init(), ocean_model_MOM.F90")
  if (associated(OS)) then
    call MOM_error(WARNING, "ocean_model_init called with an associated "// &
                   "ocean_state_type structure. Model is already initialized.")
    return
  endif
  allocate(OS)

  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
  if (.not.OS%is_ocean_pe) return

  call time_interp_external_init

  OS%Time = Time_in
  call initialize_MOM(OS%Time, Time_init, param_file, OS%dirs, OS%MOM_CSp, &
                      Time_in, offline_tracer_mode=OS%offline_tracer_mode, &
                      input_restart_file=input_restart_file, &
                      diag_ptr=OS%diag, count_calls=.true., tracer_flow_CSp=OS%tracer_flow_CSp, &
                      waves_CSp=OS%Waves, ensemble_num=inst_index)
  call get_MOM_state_elements(OS%MOM_CSp, G=OS%grid, GV=OS%GV, US=OS%US, C_p=OS%C_p, &
                              C_p_scaled=OS%fluxes%C_p, use_temp=use_temperature)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "SINGLE_STEPPING_CALL", OS%single_step_call, &
                 "If true, advance the state of MOM with a single step "//&
                 "including both dynamics and thermodynamics.  If false, "//&
                 "the two phases are advanced with separate calls.", default=.true.)
  call get_param(param_file, mdl, "DT", OS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that is actually "//&
                 "used will be an integer fraction of the forcing time-step.", &
                 units="s", scale=OS%US%s_to_T, fail_if_missing=.true.)
  call get_param(param_file, mdl, "DT_THERM", OS%dt_therm, &
                 "The thermodynamic and tracer advection time step. "//&
                 "Ideally DT_THERM should be an integer multiple of DT "//&
                 "and less than the forcing or coupling time-step, unless "//&
                 "THERMO_SPANS_COUPLING is true, in which case DT_THERM "//&
                 "can be an integer multiple of the coupling timestep.  By "//&
                 "default DT_THERM is set to DT.", &
                 units="s", default=OS%US%T_to_s*OS%dt, scale=OS%US%s_to_T)
  call get_param(param_file, "MOM", "THERMO_SPANS_COUPLING", OS%thermo_spans_coupling, &
                 "If true, the MOM will take thermodynamic and tracer "//&
                 "timesteps that can be longer than the coupling timestep. "//&
                 "The actual thermodynamic timestep that is used in this "//&
                 "case is the largest integer multiple of the coupling "//&
                 "timestep that is less than or equal to DT_THERM.", default=.false.)
  call get_param(param_file, mdl, "DIABATIC_FIRST", OS%diabatic_first, &
                 "If true, apply diabatic and thermodynamic processes, "//&
                 "including buoyancy forcing and mass gain or loss, "//&
                 "before stepping the dynamics forward.", default=.false.)

  call get_param(param_file, mdl, "RESTART_CONTROL", OS%Restart_control, &
                 "An integer whose bits encode which restart files are "//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd "//&
                 "(bit 0) for a non-time-stamped file.  A restart file "//&
                 "will be saved at the end of the run segment for any "//&
                 "non-negative value.", default=1)
  call get_param(param_file, mdl, "OCEAN_SURFACE_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the "//&
                 "staggering of the surface velocity field that is "//&
                 "returned to the coupler.  Valid values include "//&
                 "'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; Ocean_sfc%stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ocean_sfc%stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ocean_sfc%stagger = CGRID_NE
  else ; call MOM_error(FATAL,"ocean_model_init: OCEAN_SURFACE_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif

  call get_param(param_file, mdl, "EPS_OMESH",OS%eps_omesh, &
                 "Maximum allowable difference between ESMF mesh and "//&
                 "MOM6 domain coordinates in nuopc cap.", &
                 units="degrees", default=1.e-4)
  call get_param(param_file, mdl, "RESTORE_SALINITY",OS%restore_salinity, &
                 "If true, the coupled driver will add a globally-balanced "//&
                 "fresh-water flux that drives sea-surface salinity "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RESTORE_TEMPERATURE",OS%restore_temp, &
                 "If true, the coupled driver will add a "//&
                 "heat flux that drives sea-surface temperature "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RHO_0", Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "G_EARTH", G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mdl, "ICE_SHELF",  OS%use_ice_shelf, &
                 "If true, enables the ice shelf model.", default=.false.)

  call get_param(param_file, mdl, "ICEBERGS_APPLY_RIGID_BOUNDARY",  OS%icebergs_alter_ocean, &
                 "If true, allows icebergs to change boundary condition felt by ocean", default=.false.)

  OS%press_to_z = 1.0/(Rho0*G_Earth)

  call get_param(param_file, mdl, "HFREEZE", HFrz, &
                 "If HFREEZE > 0, melt potential will be computed. The actual depth "//&
                 "over which melt potential is computed will be min(HFREEZE, OBLD), "//&
                 "where OBLD is the boundary layer depth. If HFREEZE <= 0 (default), "//&
                 "melt potential will not be computed.", units="m", default=-1.0, do_not_log=.true.)

  if (HFrz .gt. 0.0) then
    use_melt_pot=.true.
  else
    use_melt_pot=.false.
  endif

  call get_param(param_file, mdl, "USE_WAVES", OS%Use_Waves, &
       "If true, enables surface wave modules.", default=.false.)
  call get_param(param_file, mdl, "USE_MARBL_TRACERS", OS%use_MARBL, &
                 default=.false., do_not_log=.true.)

  !   Consider using a run-time flag to determine whether to do the diagnostic
  ! vertical integrals, since the related 3-d sums are not negligible in cost.
  call allocate_surface_state(OS%sfc_state, OS%grid, use_temperature, &
                              do_integrals=.true., gas_fields_ocn=gas_fields_ocn, &
                              use_meltpot=use_melt_pot, use_iceshelves=OS%use_ice_shelf, &
                              use_marbl_tracers=OS%use_MARBL)

  call surface_forcing_init(Time_in, OS%grid, OS%US, param_file, OS%diag, &
                            OS%forcing_CSp, OS%restore_salinity, OS%restore_temp, OS%use_waves)

  if (OS%use_ice_shelf)  then
    call initialize_ice_shelf(param_file, OS%grid, OS%Time, OS%ice_shelf_CSp, &
                              OS%diag, Time_init, OS%dirs%output_directory, OS%forces, OS%fluxes)
  endif
  if (OS%icebergs_alter_ocean)  then
    call marine_ice_init(OS%Time, OS%grid, param_file, OS%diag, OS%marine_ice_CSp)
    if (.not. OS%use_ice_shelf) &
      call allocate_forcing_type(OS%grid, OS%fluxes, shelf=.true.)
  endif

  if (OS%Use_Waves) then
    call get_param(param_file, mdl, "WAVE_METHOD", OS%wave_method, default="EMPTY", do_not_log=.true.)
  endif
  call allocate_forcing_type(OS%grid, OS%fluxes, waves=.true., lamult=(trim(OS%wave_method)=="EFACTOR"))

  ! MOM_wave_interface_init is called regardless of the value of USE_WAVES because
  ! it also initializes statistical waves.
  call MOM_wave_interface_init(OS%Time, OS%grid, OS%GV, OS%US, param_file, OS%Waves, OS%diag)

  if (associated(OS%grid%Domain%maskmap)) then
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%diag, maskmap=OS%grid%Domain%maskmap, &
                                      gas_fields_ocn=gas_fields_ocn)
  else
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%diag, gas_fields_ocn=gas_fields_ocn)
  endif

  ! This call can only occur here if the coupler_bc_type variables have been
  ! initialized already using the information from gas_fields_ocn.
  if (present(gas_fields_ocn)) then
    call coupler_type_set_diags(Ocean_sfc%fields, "ocean_sfc", &
                                Ocean_sfc%axes(1:2), Time_in)

    call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

    if (OS%use_ice_shelf .and. allocated(OS%sfc_state%frazil)) &
      call adjust_ice_sheet_frazil(OS%sfc_state, OS%fluxes, OS%Ice_shelf_CSp)

    call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)

  endif

  call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

  if (OS%use_ice_shelf .and. allocated(OS%sfc_state%frazil)) &
    call adjust_ice_sheet_frazil(OS%sfc_state, OS%fluxes, OS%Ice_shelf_CSp)

! get number of processors and PE list for stocasthci physics initialization
  call get_param(param_file, mdl, "DO_SPPT", OS%do_sppt, &
                 "If true, then stochastically perturb the thermodynamic "//&
                 "tendencies of T,S, and h.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  call get_param(param_file, mdl, "PERT_EPBL", OS%pert_epbl, &
                 "If true, then stochastically perturb the kinetic energy "//&
                 "production and dissipation terms.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)

  call close_param_file(param_file)
  call diag_mediator_close_registration(OS%diag)

  if (is_root_pe()) &
    write(stdout,'(/12x,a/)') '======== COMPLETED MOM INITIALIZATION ========'

  call callTree_leave("ocean_model_init(")
end procedure ocean_model_init
module procedure update_ocean_model
  type(time_type) :: Master_time ! This allows step_MOM to temporarily change
  type(time_type) :: Time1       ! The value of the ocean model's time at the
  integer :: index_bnds(4)       ! The computational domain index bounds in the
  real :: weight          ! Flux accumulation weight
  real :: dt_coupling     ! The coupling time step in rescaled seconds [T ~> s].
  integer :: nts          ! The number of baroclinic dynamics time steps
  real :: dt_therm        ! A limited and quantized version of OS%dt_therm [T ~> s]
  real :: dt_dyn          ! The dynamics time step [T ~> s]
  real :: dtdia           ! The diabatic time step [T ~> s]
  real :: t_elapsed_seg   ! The elapsed time in this update segment [T ~> s]
  integer :: n, n_max, n_last_thermo
  type(time_type) :: Time2  ! A temporary time.
  logical :: thermo_does_span_coupling ! If true, thermodynamic forcing spans
  logical :: do_dyn       ! If true, step the ocean dynamics and transport.
  logical :: do_thermo    ! If true, step the ocean thermodynamics.
  logical :: step_thermo           ! If true, take a thermodynamic step.
  integer :: secs, days
  integer :: is, ie, js, je
  call callTree_enter("update_ocean_model(), MOM_ocean_model_nuopc.F90")
  call get_time(Ocean_coupling_time_step, secs, days)
  dt_coupling = OS%US%s_to_T*(86400.0*real(days) + real(secs))

  if (time_start_update /= OS%Time) then
    call MOM_error(WARNING, "update_ocean_model: internal clock does not "//&
                            "agree with time_start_update argument.")
  endif
  if (.not.associated(OS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

  do_dyn = .true. ; if (present(update_dyn)) do_dyn = update_dyn
  do_thermo = .true. ; if (present(update_thermo)) do_thermo = update_thermo
  ! This is benign but not necessary if ocean_model_init_sfc was called or if
  ! OS%sfc_state%tr_fields was spawned in ocean_model_init.  Consider removing it.
  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  ! Translate Ice_ocean_boundary into fluxes.
  call mpp_get_compute_domain(Ocean_sfc%Domain, index_bnds(1), index_bnds(2), &
                              index_bnds(3), index_bnds(4))

  weight = 1.0

  call convert_IOB_to_forces(Ice_ocean_boundary, OS%forces, index_bnds, OS%Time, &
                             OS%grid, OS%US, OS%forcing_CSp)

  if (OS%fluxes%fluxes_used) then

    ! enable_averages() is necessary to post forcing fields to diagnostics
    call enable_averages(dt_coupling, OS%Time + Ocean_coupling_time_step, OS%diag)

    if (do_thermo) &
      call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%fluxes, index_bnds, OS%Time, dt_coupling, &
                               OS%grid, OS%US, OS%forcing_CSp, OS%sfc_state, &
                               OS%restore_salinity, OS%restore_temp)

    ! Add ice shelf fluxes
    if (OS%use_ice_shelf) then
      if (do_thermo) &
        call shelf_calc_flux(OS%sfc_state, OS%fluxes, OS%Time, dt_coupling, OS%Ice_shelf_CSp)
      if (do_dyn) &
        call add_shelf_forces(OS%grid, OS%US, OS%Ice_shelf_CSp, OS%forces)
    endif
    if (OS%icebergs_alter_ocean)  then
      if (do_dyn) &
        call iceberg_forces(OS%grid, OS%forces, OS%use_ice_shelf, &
                            OS%sfc_state, dt_coupling, OS%marine_ice_CSp)
      if (do_thermo) &
        call iceberg_fluxes(OS%grid, OS%US, OS%fluxes, OS%use_ice_shelf, &
                          OS%sfc_state, dt_coupling, OS%marine_ice_CSp)
    endif

    ! Fields that exist in both the forcing and mech_forcing types must be copied.
    call copy_common_forcing_fields(OS%forces, OS%fluxes, OS%grid, skip_pres=.true.)

#ifdef _USE_GENERIC_TRACER
    call enable_averages(dt_coupling, OS%Time + Ocean_coupling_time_step, OS%diag) !Is this needed?
    call MOM_generic_tracer_fluxes_accumulate(OS%fluxes, weight) !here weight=1, just saving the current fluxes
#endif
  else
    OS%flux_tmp%C_p = OS%fluxes%C_p
    if (do_thermo) &
      call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%flux_tmp, index_bnds, OS%Time, dt_coupling, &
                               OS%grid, OS%US, OS%forcing_CSp, OS%sfc_state, OS%restore_salinity, OS%restore_temp)

    if (OS%use_ice_shelf) then
      if (do_thermo) &
        call shelf_calc_flux(OS%sfc_state, OS%flux_tmp, OS%Time, dt_coupling, OS%Ice_shelf_CSp)
      if (do_dyn) &
        call add_shelf_forces(OS%grid, OS%US, OS%Ice_shelf_CSp, OS%forces)
    endif
    if (OS%icebergs_alter_ocean)  then
      if (do_dyn) &
        call iceberg_forces(OS%grid, OS%forces, OS%use_ice_shelf, &
                            OS%sfc_state, dt_coupling, OS%marine_ice_CSp)
      if (do_thermo) &
        call iceberg_fluxes(OS%grid, OS%US, OS%flux_tmp, OS%use_ice_shelf, &
                          OS%sfc_state, dt_coupling, OS%marine_ice_CSp)
    endif

    call forcing_accumulate(OS%flux_tmp, OS%forces, OS%fluxes, OS%grid, weight)
    ! Some of the fields that exist in both the forcing and mech_forcing types
    ! (e.g., ustar) are time-averages must be copied back to the forces type.
    call copy_back_forcing_fields(OS%fluxes, OS%forces, OS%grid)

#ifdef _USE_GENERIC_TRACER
    call MOM_generic_tracer_fluxes_accumulate(OS%flux_tmp, weight) !weight of the current flux in the running average
#endif
  endif
  call set_derived_forcing_fields(OS%forces, OS%fluxes, OS%grid, OS%US, OS%GV%Rho0)
  call set_net_mass_forcing(OS%fluxes, OS%forces, OS%grid, OS%US)

  if (OS%use_waves) then
    if (OS%wave_method /= "EFACTOR") then
      call Update_Surface_Waves(OS%grid, OS%GV, OS%US, OS%time, ocean_coupling_time_step, OS%waves, OS%forces)
    endif
  endif

  if (OS%nstep==0) then
    call finish_MOM_initialization(OS%Time, OS%dirs, OS%MOM_CSp)
  endif

  if (do_thermo) &
    call call_tracer_set_forcing(OS%sfc_state, OS%fluxes, OS%Time, &
                                 real_to_time_type(dt_coupling), OS%grid, OS%US, OS%GV%Rho0, &
                                 OS%tracer_flow_CSp)

  call disable_averaging(OS%diag)
  Master_time = OS%Time ; Time1 = OS%Time

  if (OS%offline_tracer_mode) then
    call step_offline(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp)
  elseif ((.not.do_thermo) .or. (.not.do_dyn)) then
    ! The call sequence is being orchestrated from outside of update_ocean_model.
    call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp, &
                  Waves=OS%Waves, do_dynamics=do_thermo, do_thermodynamics=do_dyn, &
                  reset_therm=Ocn_fluxes_used)
 !### What to do with these?   , start_cycle=(n==1), end_cycle=.false., cycle_length=dt_coupling)

  elseif (OS%single_step_call) then
    call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp, Waves=OS%Waves)
  else
    n_max = 1 ; if (dt_coupling > OS%dt) n_max = ceiling(dt_coupling/OS%dt - 0.001)
    dt_dyn = dt_coupling / real(n_max)
    thermo_does_span_coupling = (OS%thermo_spans_coupling .and. &
                                (OS%dt_therm > 1.5*dt_coupling))

    if (thermo_does_span_coupling) then
      dt_therm = dt_coupling * floor(OS%dt_therm / dt_coupling + 0.001)
      nts = floor(dt_therm/dt_dyn + 0.001)
    else
      nts = MAX(1,MIN(n_max,floor(OS%dt_therm/dt_dyn + 0.001)))
      n_last_thermo = 0
    endif

    Time2 = Time1 ; t_elapsed_seg = 0.0
    do n=1,n_max
      if (OS%diabatic_first) then
        if (thermo_does_span_coupling) call MOM_error(FATAL, &
            "MOM is not yet set up to have restarts that work with "//&
            "THERMO_SPANS_COUPLING and DIABATIC_FIRST.")
        if (modulo(n-1,nts)==0) then
          dtdia = dt_dyn*min(nts,n_max-(n-1))
          call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time2, dtdia, OS%MOM_CSp, &
                        Waves=OS%Waves, do_dynamics=.false., do_thermodynamics=.true., &
                        start_cycle=(n==1), end_cycle=.false., cycle_length=dt_coupling)
        endif

        call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time2, dt_dyn, OS%MOM_CSp, &
                      Waves=OS%Waves, do_dynamics=.true., do_thermodynamics=.false., &
                      start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_coupling)
      else
        call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time2, dt_dyn, OS%MOM_CSp, &
                      Waves=OS%Waves, do_dynamics=.true., do_thermodynamics=.false., &
                      start_cycle=(n==1), end_cycle=.false., cycle_length=dt_coupling)

        step_thermo = .false.
        if (thermo_does_span_coupling) then
          dtdia = dt_therm
          step_thermo = MOM_state_is_synchronized(OS%MOM_CSp, adv_dyn=.true.)
        elseif ((modulo(n,nts)==0) .or. (n==n_max)) then
          dtdia = dt_dyn*(n - n_last_thermo)
          n_last_thermo = n
          step_thermo = .true.
        endif

        if (step_thermo) then
          ! Back up Time2 to the start of the thermodynamic segment.
          Time2 = Time2 - set_time(int(floor(OS%US%T_to_s*(dtdia - dt_dyn) + 0.5)))
          call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time2, dtdia, OS%MOM_CSp, &
                        Waves=OS%Waves, do_dynamics=.false., do_thermodynamics=.true., &
                        start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_coupling)
        endif
      endif

      t_elapsed_seg = t_elapsed_seg + dt_dyn
      Time2 = Time1 + set_time(int(floor(OS%US%T_to_s*t_elapsed_seg + 0.5)))
    enddo
  endif

  OS%Time = Master_time + Ocean_coupling_time_step
  OS%nstep = OS%nstep + 1

  call mech_forcing_diags(OS%forces, dt_coupling, OS%grid, OS%Time, OS%diag, OS%forcing_CSp%handles)

  if (OS%fluxes%fluxes_used) then
    if (cesm_coupled) then
      call forcing_diagnostics(OS%fluxes, OS%sfc_state, OS%grid, OS%US, OS%Time, OS%diag, &
                               OS%forcing_CSp%handles, enthalpy=.true.)
    else
      call forcing_diagnostics(OS%fluxes, OS%sfc_state, OS%grid, OS%US, OS%Time, OS%diag, OS%forcing_CSp%handles)
    endif
  endif

  !only make ice-shelf frazil adjustments if sfc_state%frazil was updated (do_thermo=True)
  if (do_thermo .and. OS%use_ice_shelf .and. allocated(OS%sfc_state%frazil)) &
    call adjust_ice_sheet_frazil(OS%sfc_state, OS%fluxes, OS%Ice_shelf_CSp)

! Translate state into Ocean.
!  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, &
!                                   Ice_ocean_boundary%p, OS%press_to_z)
  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)
  call coupler_type_send_data(Ocean_sfc%fields, OS%Time)

  call callTree_leave("update_ocean_model()")
end procedure update_ocean_model
module procedure ocean_model_restart
  if (.not.MOM_state_is_synchronized(OS%MOM_CSp)) &
      call MOM_error(WARNING, "End of MOM_main reached with inconsistent "//&
         "dynamics and advective times.  Additional restart fields "//&
         "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_restart "//&
      "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
      "restart files can only be created after the buoyancy forcing is applied.")

  if (present(restartname)) then
    call save_MOM_restart(OS%MOM_CSp, OS%dirs%restart_output_dir, OS%Time, &
        OS%grid, GV=OS%GV, filename=restartname, num_rest_files=num_rest_files)
    call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
         OS%dirs%restart_output_dir) ! Is this needed?
    if (OS%use_ice_shelf) then
      call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, &
           OS%dirs%restart_output_dir)
    endif
  else
    if (BTEST(OS%Restart_control,1)) then
      call save_MOM_restart(OS%MOM_CSp, OS%dirs%restart_output_dir, OS%Time, &
          OS%grid, time_stamped=.true., GV=OS%GV)
      call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
           OS%dirs%restart_output_dir, time_stamped=.true.)
      if (OS%use_ice_shelf) then
        call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir, .true.)
      endif
    endif
    if (BTEST(OS%Restart_control,0)) then
      call save_MOM_restart(OS%MOM_CSp, OS%dirs%restart_output_dir, OS%Time, &
          OS%grid, GV=OS%GV)
      call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
           OS%dirs%restart_output_dir)
      if (OS%use_ice_shelf) then
        call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
      endif
    endif
  endif
  if (present(stoch_restartname)) then
    if (OS%do_sppt .OR. OS%pert_epbl) then
      call write_stoch_restart_ocn('RESTART/'//trim(stoch_restartname))
    endif
  endif

end procedure ocean_model_restart
module procedure ocean_model_end
  if (write_restart) call ocean_model_save_restart(Ocean_state, Time)
  call diag_mediator_end(Time, Ocean_state%diag, end_diag_manager=.true.)
  call MOM_end(Ocean_state%MOM_CSp)
  if (Ocean_state%use_ice_shelf) call ice_shelf_end(Ocean_state%Ice_shelf_CSp)
end procedure ocean_model_end
module procedure ocean_model_save_restart
  character(len=200) :: restart_dir
  if (.not.MOM_state_is_synchronized(OS%MOM_CSp)) &
    call MOM_error(WARNING, "ocean_model_save_restart called with inconsistent "//&
         "dynamics and advective times.  Additional restart fields "//&
         "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_save_restart "//&
       "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
       "restart files can only be created after the buoyancy forcing is applied.")

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = OS%dirs%restart_output_dir ; endif

  call save_MOM_restart(OS%MOM_CSp, restart_dir, Time, OS%grid, GV=OS%GV)

  call forcing_save_restart(OS%forcing_CSp, OS%grid, Time, restart_dir)

  if (OS%use_ice_shelf) then
    call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
  endif
end procedure ocean_model_save_restart
module procedure initialize_ocean_public_type
  integer :: xsz, ysz, layout(2)
  integer :: isc, iec, jsc, jec
  call mpp_get_layout(input_domain,layout)
  call mpp_get_global_domain(input_domain, xsize=xsz, ysize=ysz)
  if (PRESENT(maskmap)) then
    call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain, maskmap=maskmap)
  else
    call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain)
  endif
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc, iec, jsc, jec)

  allocate(Ocean_sfc%t_surf (isc:iec,jsc:jec),  &  ! time averaged sst (Kelvin) passed to atmosphere/ice model
           Ocean_sfc%s_surf (isc:iec,jsc:jec),  &  ! time averaged sss (psu) passed to atmosphere/ice models
           Ocean_sfc%u_surf (isc:iec,jsc:jec),  &  ! time averaged u-current (m/sec) passed to atmosphere/ice models
           Ocean_sfc%v_surf (isc:iec,jsc:jec),  &  ! time averaged v-current (m/sec)  passed to atmosphere/ice models
           Ocean_sfc%sea_lev(isc:iec,jsc:jec),  &  ! time averaged thickness of top model grid cell (m) plus
                                                   ! patm/rho0/grav
           Ocean_sfc%frazil (isc:iec,jsc:jec),  &  ! time accumulated frazil (J/m^2) passed to ice model
           Ocean_sfc%melt_potential(isc:iec,jsc:jec), &  ! time accumulated melt potential (J/m^2) passed to ice model
           Ocean_sfc%area   (isc:iec,jsc:jec),  &
           Ocean_sfc%OBLD   (isc:iec,jsc:jec),  &  ! ocean boundary layer depth, in m
           Ocean_sfc%fco2_ocn(isc:iec,jsc:jec), &  ! time averaged co2 flux (kg/m^2/s) passed to atmosphere model
           source=0.0)

  Ocean_sfc%axes    = diag%axesT1%handles !diag axes to be used by coupler tracer flux diagnostics

  if (present(gas_fields_ocn)) then
    call coupler_type_spawn(gas_fields_ocn, Ocean_sfc%fields, (/isc,isc,iec,iec/), &
                              (/jsc,jsc,jec,jec/), suffix = '_ocn', as_needed=.true.)
  endif

end procedure initialize_ocean_public_type
module procedure convert_state_to_ocean_type
  real :: IgR0
  character(len=48)  :: val_str
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  integer :: i, j, i0, j0, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call pass_vector(sfc_state%u, sfc_state%v, G%Domain)

  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)
  if (present(patm)) then
    ! Check that the inidicies in patm are (isc_bnd:iec_bnd,jsc_bnd:jec_bnd).
    if (.not.present(press_to_z)) call MOM_error(FATAL, &
        'convert_state_to_ocean_type: press_to_z must be present if patm is.')
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  if (sfc_state%T_is_conT) then
    ! Convert the surface T from conservative T to potential T.
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = gsw_pt_from_ct(US%S_to_ppt*sfc_state%SSS(i+i0,j+j0), &
                               US%C_to_degC*sfc_state%SST(i+i0,j+j0)) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = US%C_to_degC*sfc_state%SST(i+i0,j+j0) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  endif
  if (sfc_state%S_is_absS) then
    ! Convert the surface S from absolute salinity to practical salinity.
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = gsw_sp_from_sr(US%S_to_ppt*sfc_state%SSS(i+i0,j+j0))
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = US%S_to_ppt*sfc_state%SSS(i+i0,j+j0)
    enddo ; enddo
  endif

  if (present(patm)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = US%Z_to_m * sfc_state%sea_lev(i+i0,j+j0) + patm(i,j) * press_to_z
      Ocean_sfc%area(i,j) = US%L_to_m**2 * G%areaT(i+i0,j+j0)
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = US%Z_to_m * sfc_state%sea_lev(i+i0,j+j0)
      Ocean_sfc%area(i,j) = US%L_to_m**2 * G%areaT(i+i0,j+j0)
    enddo ; enddo
  endif

  if (allocated(sfc_state%frazil)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%frazil(i,j) = US%Q_to_J_kg*US%RZ_to_kg_m2 * sfc_state%frazil(i+i0,j+j0)
    enddo ; enddo
  endif

  if (allocated(sfc_state%melt_potential)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%melt_potential(i,j) = US%Q_to_J_kg*US%RZ_to_kg_m2 * sfc_state%melt_potential(i+i0,j+j0)
    enddo ; enddo
  endif

  if (allocated(sfc_state%Hml)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%OBLD(i,j) = US%Z_to_m * sfc_state%Hml(i+i0,j+j0)
    enddo ; enddo
  endif

  if (allocated(sfc_state%fco2)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%fco2_ocn(i,j) = US%RZ_T_to_kg_m2s * sfc_state%fco2(i+i0,j+j0)
    enddo ; enddo
  endif

  if (Ocean_sfc%stagger == AGRID) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dT(i+i0,j+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I-1+i0,j+j0))
      Ocean_sfc%v_surf(i,j) = G%mask2dT(i+i0,j+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0,J-1+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == BGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dBu(I+i0,J+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I+i0,j+j0+1))
      Ocean_sfc%v_surf(i,j) = G%mask2dBu(I+i0,J+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0+1,J+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == CGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dCu(I+i0,j+j0) * US%L_T_to_m_s * sfc_state%u(I+i0,j+j0)
      Ocean_sfc%v_surf(i,j) = G%mask2dCv(i+i0,J+j0) * US%L_T_to_m_s * sfc_state%v(i+i0,J+j0)
    enddo ; enddo
  else
    write(val_str, '(I8)') Ocean_sfc%stagger
    call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
      "Ocean_sfc%stagger has the unrecognized value of "//trim(val_str))
  endif

  if (coupler_type_initialized(sfc_state%tr_fields)) then
    if (.not.coupler_type_initialized(Ocean_sfc%fields)) then
      call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
               "Ocean_sfc%fields has not been initialized.")
    endif
    call coupler_type_copy_data(sfc_state%tr_fields, Ocean_sfc%fields)
  endif

end procedure convert_state_to_ocean_type
module procedure ocean_model_init_sfc
  integer :: is, ie, js, je
  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

  if (OS%use_ice_shelf .and. allocated(OS%sfc_state%frazil)) &
    call adjust_ice_sheet_frazil(OS%sfc_state, OS%fluxes, OS%Ice_shelf_CSp)

  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)

end procedure ocean_model_init_sfc
module procedure ocean_model_flux_init
  logical :: OS_is_set
  integer :: verbose
  OS_is_set = .false. ; if (present(OS)) OS_is_set = associated(OS)

  ! Use this to control the verbosity of output; consider rethinking this logic later.
  verbose = 5 ; if (OS_is_set) verbose = 3
  if (present(verbosity)) verbose = verbosity

  call call_tracer_flux_init(verbosity=verbose)

end procedure ocean_model_flux_init
module procedure query_ocean_state
  logical :: undo_scaling
  undo_scaling = .false. ; if (present(unscale)) undo_scaling = unscale

  if (present(use_waves)) use_waves = OS%use_waves
  if (present(NumWaveBands)) call query_wave_properties(OS%Waves, NumBands=NumWaveBands)
  if (present(Wavenumbers) .and. undo_scaling) then
    call query_wave_properties(OS%Waves, WaveNumbers=WaveNumbers, US=OS%US)
  elseif (present(Wavenumbers)) then
    call query_wave_properties(OS%Waves, WaveNumbers=WaveNumbers)
  endif
  if (present(wave_method)) wave_method = OS%wave_method
  if (present(use_MARBL)) use_MARBL = OS%use_MARBL

end procedure query_ocean_state
module procedure Ocean_stock_pe
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT,ISTOCK_SALT
  real :: salt
  value = 0.0
  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  select case (index)
    case (ISTOCK_WATER)  ! Return the mass of fresh water in the ocean in kg.
      if (OS%GV%Boussinesq) then
        call get_ocean_stocks(OS%MOM_CSp, mass=value, on_PE_only=.true.)
      else  ! In non-Boussinesq mode, the mass of salt needs to be subtracted.
        call get_ocean_stocks(OS%MOM_CSp, mass=value, salt=salt, on_PE_only=.true.)
        value = value - salt
      endif
    case (ISTOCK_HEAT)  ! Return the heat content of the ocean in J.
      call get_ocean_stocks(OS%MOM_CSp, heat=value, on_PE_only=.true.)
    case (ISTOCK_SALT)  ! Return the mass of the salt in the ocean in kg.
      call get_ocean_stocks(OS%MOM_CSp, salt=value, on_PE_only=.true.)
    case default ; value = 0.0
  end select
  ! If the FMS coupler is changed so that Ocean_stock_PE is only called on
  ! ocean PEs, uncomment the following and eliminate the on_PE_only flags above.
  !  if (.not.is_root_pe()) value = 0.0

end procedure Ocean_stock_pe
module procedure ocean_public_type_chksum
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  root = is_root_pe()

  if (root) write(stdout,*) "BEGIN CHECKSUM(ocean_type):: ", id, timestep
  chks = field_chksum(ocn%t_surf ) ; if (root) write(stdout,100) 'ocean%t_surf   ', chks
  chks = field_chksum(ocn%s_surf ) ; if (root) write(stdout,100) 'ocean%s_surf   ', chks
  chks = field_chksum(ocn%u_surf ) ; if (root) write(stdout,100) 'ocean%u_surf   ', chks
  chks = field_chksum(ocn%v_surf ) ; if (root) write(stdout,100) 'ocean%v_surf   ', chks
  chks = field_chksum(ocn%sea_lev) ; if (root) write(stdout,100) 'ocean%sea_lev  ', chks
  chks = field_chksum(ocn%frazil ) ; if (root) write(stdout,100) 'ocean%frazil   ', chks
  chks = field_chksum(ocn%melt_potential) ; if (root) write(stdout,100) 'ocean%melt_potential   ', chks
  call coupler_type_write_chksums(ocn%fields, stdout, 'ocean%')
100 FORMAT("   CHECKSUM::",A20," = ",Z20)

end procedure ocean_public_type_chksum
module procedure get_ocean_grid
  Gridp => OS%grid
  return
end procedure get_ocean_grid
module procedure get_eps_omesh
  get_eps_omesh = OS%eps_omesh ; return
end procedure get_eps_omesh
end submodule MOM_ocean_model_nuopc_s
