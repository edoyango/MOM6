submodule (MOM) MOM_s
#include <MOM_memory.h>
  implicit none
contains
module procedure step_MOM
  type(ocean_grid_type),   pointer :: G => NULL()  ! pointer to a structure containing
  type(ocean_grid_type),   pointer :: G_in => NULL()  ! Input grid metric
  type(verticalGrid_type), pointer :: GV => NULL() ! Pointer to the vertical grid structure
  type(unit_scale_type),   pointer :: US => NULL() ! Pointer to a structure containing
  integer       :: ntstep  ! number of time steps between diabatic forcing updates
  integer       :: ntastep ! number of time steps between tracer advection updates
  integer       :: n_max  ! number of steps to take in this call
  integer :: halo_sz, dynamics_stencil
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: time_interval   ! time interval covered by this run segment [T ~> s].
  real :: dt              ! baroclinic time step [T ~> s]
  real :: dtdia           ! time step for diabatic processes [T ~> s]
  real :: dt_tr_adv       ! time step for tracer advection [T ~> s]
  real :: dt_therm        ! a limited and quantized version of CS%dt_therm [T ~> s]
  real :: dt_tradv_here   ! a further limited value of dt_tr_adv [T ~> s]
  real :: wt_end, wt_beg  ! Fractional weights of the future pressure at the end
  real :: bbl_time_int    ! The amount of time over which the calculated BBL
  real :: rel_time = 0.0  ! relative time since start of this call [T ~> s].
  logical :: do_advection    ! If true, do tracer advection.
  logical :: do_diabatic     ! If true, do diabatic update.
  logical :: thermo_does_span_coupling ! If true,thermodynamic (diabatic) forcing spans
  logical :: tradv_does_span_coupling  ! If true, tracer advection spans
  logical :: do_dyn     ! If true, dynamics are updated with this call.
  logical :: do_thermo  ! If true, thermodynamics and remapping may be applied with this call.
  logical :: debug_redundant ! If true, check redundant values on PE boundaries when debugging.
  logical :: nonblocking_p_surf_update ! A flag to indicate whether surface properties
  logical :: cycle_start ! If true, do calculations that are only done at the start of
  logical :: cycle_end  ! If true, do calculations and diagnostics that are only done at
  logical :: therm_reset ! If true, reset running sums of thermodynamic quantities.
  real :: cycle_time    ! The length of the coupled time-stepping cycle [T ~> s].
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: &
    U_star      ! The wind friction velocity, calculated using the Boussinesq reference density or
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: &
    ssh         ! sea surface height, which may be based on eta_av [Z ~> m]
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
    dz          ! Vertical distance across layers [Z ~> m]
  real, dimension(:,:,:), pointer :: &
    u => NULL(), & ! u : zonal velocity component [L T-1 ~> m s-1]
    v => NULL(), & ! v : meridional velocity component [L T-1 ~> m s-1]
    h => NULL()    ! h : layer thickness [H ~> m or kg m-2]
  real, dimension(:,:), pointer :: &
    p_surf => NULL() ! A pointer to the ocean surface pressure [R L2 T-2 ~> Pa].
  real :: I_wt_ssh  ! The inverse of the time weights [T-1 ~> s-1]
  type(time_type) :: Time_local, end_time_thermo
  type(time_type) :: Time_end_diag ! End time of a diagnostic segment, as a time type
  type(group_pass_type) :: pass_tau_ustar_psurf
  logical :: showCallTree
  type(mech_forcing), pointer :: forces     ! Mechanical forcing
  type(forcing), pointer :: fluxes          ! Boundary fluxes
  type(surface), pointer :: sfc_state_diag  ! Surface boundary fields
  integer :: turns  ! Number of quarter turns from input to model indexing
  G => CS%G ; G_in => CS%G_in ; GV => CS%GV ; US => CS%US
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  u => CS%u ; v => CS%v ; h => CS%h

  time_interval = time_int_in
  do_dyn = .true. ; if (present(do_dynamics)) do_dyn = do_dynamics
  do_thermo = .true. ; if (present(do_thermodynamics)) do_thermo = do_thermodynamics
  if (.not.(do_dyn .or. do_thermo)) call MOM_error(FATAL,"Step_MOM: "//&
    "Both do_dynamics and do_thermodynamics are false, which makes no sense.")
  cycle_start = .true. ; if (present(start_cycle)) cycle_start = start_cycle
  cycle_end = .true. ; if (present(end_cycle)) cycle_end = end_cycle
  cycle_time = time_interval ; if (present(cycle_length)) cycle_time = cycle_length
  therm_reset = cycle_start ; if (present(reset_therm)) therm_reset = reset_therm

  call cpu_clock_begin(id_clock_ocean)
  call cpu_clock_begin(id_clock_other)

  if (CS%debug) then
    call query_debugging_checks(do_redundant=debug_redundant)
    call MOM_state_chksum("Beginning of step_MOM ", u, v, h, CS%uh, CS%vh, G, GV, US)
  endif

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("step_MOM(), MOM.F90")

  ! Rotate the forces from G_in to G
  if (CS%rotate_index) then
    turns = G%HI%turns
    allocate(forces)
    call allocate_mech_forcing(forces_in, G, forces)
    call rotate_mech_forcing(forces_in, turns, forces)

    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, G, fluxes, turns=turns)
    call rotate_forcing(fluxes_in, fluxes, turns)
  else
    forces => forces_in
    fluxes => fluxes_in
  endif

  ! Homogenize the forces
  if (CS%homogenize_forcings) then
    ! Homogenize all forcing and fluxes fields.
    call homogenize_mech_forcing(forces, G, US, GV%Rho0, CS%update_ustar)
    ! Note the following computes the mean ustar as the mean of ustar rather than
    !  ustar of the mean of tau.
    call homogenize_forcing(fluxes, G, GV, US)
    if (CS%update_ustar) then
      ! These calls corrects the ustar values
      call copy_common_forcing_fields(forces, fluxes, G)
      call set_derived_forcing_fields(forces, fluxes, G, US, GV%Rho0)
    endif
  endif

  ! This will be replaced later with the pressures from forces or fluxes if they are available.
  if (associated(CS%tv%p_surf)) CS%tv%p_surf(:,:) = 0.0

  ! First determine the time step that is consistent with this call and an
  ! integer fraction of time_interval.
  if (do_dyn) then
    n_max = 1
    if (time_interval > CS%dt) n_max = ceiling(time_interval/CS%dt - 0.001)

    dt = time_interval / real(n_max)
    thermo_does_span_coupling = (CS%thermo_spans_coupling .and. &
                                (CS%dt_therm > 1.5*cycle_time))
    tradv_does_span_coupling = (CS%tradv_spans_coupling .and. &
                                (CS%dt_tr_adv > 1.5*cycle_time))
    if (thermo_does_span_coupling) then
      ! Set dt_therm to be an integer multiple of the coupling time step.
      dt_therm = cycle_time * floor(CS%dt_therm / cycle_time + 0.001)
      ntstep = floor(dt_therm/dt + 0.001)
    elseif (.not.do_thermo) then
      dt_therm = CS%dt_therm
      if (present(cycle_length)) dt_therm = min(CS%dt_therm, cycle_length)
      ntstep = 1 ! ntstep is initialized to avoid an error in a secondary logical test,
                 ! but the nonzero value of ntstep does not matter when do_thermo is false.
    else
      ntstep = MAX(1, MIN(n_max, floor(CS%dt_therm/dt + 0.001)))
      dt_therm = dt*ntstep
    endif
    if (tradv_does_span_coupling) then
      ! Set dt_tr_adv to be an integer multiple of the coupling time step.
      dt_tr_adv = cycle_time * floor(CS%dt_tr_adv / cycle_time + 0.001)
      ntastep = floor(dt_tr_adv/dt + 0.001)
    elseif (.not.do_thermo) then
      dt_tr_adv = CS%dt_tr_adv
      if (present(cycle_length)) dt_tr_adv = min(CS%dt_tr_adv, cycle_length)
      ! ntastep is not used.
    else
      ntastep = MAX(1, MIN(n_max, floor(CS%dt_tr_adv/dt + 0.001)))
      dt_tr_adv = dt*ntastep
    endif

    !---------- Initiate group halo pass of the forcing fields
    call cpu_clock_begin(id_clock_pass)
    ! Halo updates for surface pressure need to be completed before calling calc_resoln_function
    ! among other routines if the surface pressure is used in the equation of state.
    nonblocking_p_surf_update = G%nonblocking_updates .and. &
        .not.(associated(CS%tv%p_surf) .and. associated(forces%p_surf) .and. &
              allocated(CS%tv%SpV_avg) .and. associated(CS%tv%T))
    if (.not.associated(forces%taux) .or. .not.associated(forces%tauy)) &
         call MOM_error(FATAL,'step_MOM:forces%taux,tauy not associated')
    call create_group_pass(pass_tau_ustar_psurf, forces%taux, forces%tauy, G%Domain)
    if (associated(forces%ustar)) &
      call create_group_pass(pass_tau_ustar_psurf, forces%ustar, G%Domain)
    if (associated(forces%tau_mag)) &
      call create_group_pass(pass_tau_ustar_psurf, forces%tau_mag, G%Domain)
    if (associated(forces%p_surf)) &
      call create_group_pass(pass_tau_ustar_psurf, forces%p_surf, G%Domain)
    if (nonblocking_p_surf_update) then
      call start_group_pass(pass_tau_ustar_psurf, G%Domain)
    else
      call do_group_pass(pass_tau_ustar_psurf, G%Domain)
    endif
    call cpu_clock_end(id_clock_pass)

    if (associated(forces%p_surf)) p_surf => forces%p_surf
    if (.not.associated(forces%p_surf)) CS%interp_p_surf = .false.
    if (associated(CS%tv%p_surf) .and. associated(forces%p_surf)) then
      do j=jsd,jed ; do i=isd,ied ; CS%tv%p_surf(i,j) = forces%p_surf(i,j) ; enddo ; enddo

      if (allocated(CS%tv%SpV_avg) .and. associated(CS%tv%T)) then
        ! The internal ocean state depends on the surface pressues, so update SpV_avg.
        dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
        call calc_derived_thermo(CS%tv, h, G, GV, US, halo=dynamics_stencil, debug=CS%debug)
      endif
    endif

  else
    ! This step only updates the thermodynamics so setting timesteps is simpler.
    n_max = 1
    if ((time_interval > CS%dt_therm) .and. (CS%dt_therm > 0.0)) &
      n_max = ceiling(time_interval/CS%dt_therm - 0.001)

    dt = time_interval / real(n_max)
    dt_therm = dt ; ntstep = 1

    if (CS%UseWaves .and. associated(fluxes%ustar)) &
      call pass_var(fluxes%ustar, G%Domain, clock=id_clock_pass, halo=1)
    if (CS%UseWaves .and. associated(fluxes%tau_mag)) &
      call pass_var(fluxes%tau_mag, G%Domain, clock=id_clock_pass, halo=1)

    if (associated(fluxes%p_surf)) p_surf => fluxes%p_surf
    if (associated(CS%tv%p_surf) .and. associated(fluxes%p_surf)) then
      do j=js,je ; do i=is,ie ; CS%tv%p_surf(i,j) = fluxes%p_surf(i,j) ; enddo ; enddo
      if (allocated(CS%tv%SpV_avg)) then
        call pass_var(CS%tv%p_surf, G%Domain, clock=id_clock_pass)
        ! The internal ocean state depends on the surface pressues, so update SpV_avg.
        call extract_diabatic_member(CS%diabatic_CSp, diabatic_halo=halo_sz)
        halo_sz = max(halo_sz, 1)
        call calc_derived_thermo(CS%tv, h, G, GV, US, halo=halo_sz, debug=CS%debug)
      endif
    endif
  endif

  if (therm_reset) then
    CS%time_in_thermo_cycle = 0.0
    if (associated(CS%tv%frazil)) then
      CS%tv%frazil(:,:) = 0.0
      CS%tv%frazil_was_reset = .true.
    endif
    if (associated(CS%tv%salt_deficit))  CS%tv%salt_deficit(:,:)  = 0.0
    if (associated(CS%tv%TempxPmE))      CS%tv%TempxPmE(:,:)      = 0.0
    if (associated(CS%tv%internal_heat)) CS%tv%internal_heat(:,:) = 0.0
  endif

  if (cycle_start) then
    CS%time_in_cycle = 0.0
    do j=js,je ; do i=is,ie ; CS%ssh_rint(i,j) = 0.0 ; enddo ; enddo

    if (CS%VarMix%use_variable_mixing) then
      Time_end_diag = Time_start + real_to_time(cycle_time, unscale=US%T_to_s)
      call enable_averages(cycle_time, Time_end_diag, CS%diag)
      call calc_resoln_function(h, CS%tv, G, GV, US, CS%VarMix, CS%MEKE, CS%OBC, dt)
      call calc_depth_function(G, CS%VarMix)
      call disable_averaging(CS%diag)
    endif
  endif
  ! advance the random pattern if stochastic physics is active
  if (CS%stoch_CS%do_sppt .OR. CS%stoch_CS%pert_epbl .OR. CS%stoch_CS%do_skeb) &
    call update_stochastics(CS%stoch_CS)

  if (do_dyn) then
    if (nonblocking_p_surf_update) &
      call complete_group_pass(pass_tau_ustar_psurf, G%Domain, clock=id_clock_pass)

    if (CS%interp_p_surf) then
      if (.not.associated(CS%p_surf_end))   allocate(CS%p_surf_end(isd:ied,jsd:jed))
      if (.not.associated(CS%p_surf_begin)) allocate(CS%p_surf_begin(isd:ied,jsd:jed))
      if (.not.CS%p_surf_prev_set) then
        do j=jsd,jed ; do i=isd,ied
          CS%p_surf_prev(i,j) = forces%p_surf(i,j)
        enddo ; enddo
        CS%p_surf_prev_set = .true.
      endif
    else
      CS%p_surf_end  => forces%p_surf
    endif
    if (CS%UseWaves) then
      ! Update wave information, which is presently kept static over each call to step_mom
      Time_end_diag = Time_start + real_to_time(time_interval, unscale=US%T_to_s)
      call enable_averages(time_interval, Time_end_diag, CS%diag)
      call find_ustar(forces, CS%tv, U_star, G, GV, US, halo=1)
      call thickness_to_dz(h, CS%tv, dz, G, GV, US, halo_size=1)
      call Update_Stokes_Drift(G, GV, US, Waves, dz, U_star, time_interval, do_dyn)
      call disable_averaging(CS%diag)
    endif
  else ! not do_dyn.
    if (CS%UseWaves) then ! Diagnostics are not enabled in this call.
      call find_ustar(fluxes, CS%tv, U_star, G, GV, US, halo=1)
      call thickness_to_dz(h, CS%tv, dz, G, GV, US, halo_size=1)
      call Update_Stokes_Drift(G, GV, US, Waves, dz, U_star, time_interval, do_dyn)
    endif
  endif

  if (CS%debug) then
    if (cycle_start) &
      call MOM_state_chksum("Before steps ", u, v, h, CS%uh, CS%vh, G, GV, US)
    if (cycle_start .and. debug_redundant) &
      call check_redundant("Before steps ", u, v, G, unscale=US%L_T_to_m_s)
    if (do_dyn) call MOM_mech_forcing_chksum("Before steps", forces, G, US, haloshift=0)
    if (do_dyn .and. debug_redundant) &
      call check_redundant("Before steps ", forces%taux, forces%tauy, G, &
                                     unscale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  endif
  call cpu_clock_end(id_clock_other)

  rel_time = 0.0
  do n=1,n_max
    rel_time = rel_time + dt ! The relative time at the end of the step.
    ! Set the universally visible time to the middle of the time step.
    CS%Time = Time_start + real_to_time(rel_time - 0.5*dt, unscale=US%T_to_s)
    ! Set the local time to the end of the time step.
    Time_local = Time_start + real_to_time(rel_time, unscale=US%T_to_s)

    if (showCallTree) call callTree_enter("DT cycles (step_MOM) n=",n)

    ! Update the vertically extensive diagnostic grids so that they are
    ! referenced to the beginning timestep
    call diag_update_remap_grids(CS%diag, update_intensive = .false., update_extensive = .true. )

    !===========================================================================
    ! This is the first place where the diabatic processes and remapping could occur.
    if (CS%diabatic_first .and. (CS%t_dyn_rel_adv==0.0) .and. do_thermo) then ! do thermodynamics.

      if (.not.do_dyn) then
        dtdia = dt
      elseif (thermo_does_span_coupling) then
        dtdia = dt_therm
        if ((fluxes%dt_buoy_accum > 0.0) .and. (dtdia > time_interval) .and. &
            (abs(fluxes%dt_buoy_accum - dtdia) > 1e-6*dtdia)) then
          call MOM_error(FATAL, "step_MOM: Mismatch between long thermodynamic "//&
            "timestep and time over which buoyancy fluxes have been accumulated.")
        endif
        call MOM_error(FATAL, "MOM is not yet set up to have restarts that work "//&
          "with THERMO_SPANS_COUPLING and DIABATIC_FIRST.")
      else
        dtdia = dt*min(ntstep,n_max-(n-1))
      endif

      end_time_thermo = Time_local
      if (dtdia > dt) then
        ! If necessary, temporarily reset CS%Time to the center of the period covered
        ! by the call to step_MOM_thermo, noting that they begin at the same time.
        CS%Time = CS%Time + real_to_time(0.5*(dtdia-dt), unscale=US%T_to_s)
        ! The end-time of the diagnostic interval needs to be set ahead if there
        ! are multiple dynamic time steps worth of thermodynamics applied here.
        end_time_thermo = Time_local + real_to_time(dtdia-dt, unscale=US%T_to_s)
      endif

      ! Apply diabatic forcing, do mixing, and regrid.
      call step_MOM_thermo(CS, G, GV, US, u, v, h, CS%tv, fluxes, dtdia, &
                           end_time_thermo, .true., Waves=Waves)
      if ( CS%use_ALE_algorithm ) &
        call ALE_regridding_and_remapping(CS, G, GV, US, u, v, h, CS%tv, dtdia, Time_local)
      call post_diabatic_halo_updates(CS, G, GV, US, u, v, h, CS%tv)
      CS%time_in_thermo_cycle = CS%time_in_thermo_cycle + dtdia

      ! The diabatic processes are now ahead of the dynamics by dtdia.
      CS%t_dyn_rel_thermo = -dtdia
      if (showCallTree) call callTree_waypoint("finished diabatic_first (step_MOM)")

      if (dtdia > dt) & ! Reset CS%Time to its previous value.
        CS%Time = Time_start + real_to_time(rel_time - 0.5*dt, unscale=US%T_to_s)
    endif ! end of block "(CS%diabatic_first .and. (CS%t_dyn_rel_adv==0.0))"

    if (do_dyn) then
      ! Store pre-dynamics thicknesses for proper diagnostic remapping for transports or
      ! advective tendencies.  If there are more than one dynamics steps per advective
      ! step (i.e DT_THERM > DT), this needs to be stored at the first dynamics call.
      if (.not.CS%preadv_h_stored .and. (CS%t_dyn_rel_adv == 0.)) then
        call diag_copy_diag_to_storage(CS%diag_pre_dyn, h, CS%diag)
        CS%preadv_h_stored = .true.
      endif

      ! The pre-dynamics velocities might be stored for debugging truncations.
      if (associated(CS%u_prev) .and. associated(CS%v_prev)) then
        do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
          CS%u_prev(I,j,k) = u(I,j,k)
        enddo ; enddo ; enddo
        do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
          CS%v_prev(I,j,k) = v(i,J,k)
        enddo ; enddo ; enddo
      endif

      if (CS%interface_filter_dt_bug) then
        dt_tradv_here = dt_therm
        if (do_thermo .and. do_dyn .and. .not.thermo_does_span_coupling) &
          dt_tradv_here = dt*min(ntstep, n_max-n+1)
      else
        dt_tradv_here = dt_tr_adv
        if (do_thermo .and. do_dyn .and. .not.tradv_does_span_coupling) &
          dt_tradv_here = dt*min(ntstep, n_max-n+1)
      endif

      ! Indicate whether the bottom boundary layer properties need to be
      ! recalculated, and if so for how long an interval they are valid.
      bbl_time_int = 0.0
      if (do_thermo) then
        if ((CS%t_dyn_rel_adv == 0.0) .or. (n==1)) &
          bbl_time_int = max(dt, min(dt_therm - CS%t_dyn_rel_adv, dt*(1+n_max-n)) )
      else
        if ((CS%t_dyn_rel_adv == 0.0) .or. ((n==1) .and. cycle_start)) &
          bbl_time_int = min(dt_therm, cycle_time)
      endif

      if (CS%interp_p_surf) then
        wt_end = real(n) / real(n_max)
        wt_beg = real(n-1) / real(n_max)
        do j=jsd,jed ; do i=isd,ied
          CS%p_surf_end(i,j) = wt_end * forces%p_surf(i,j) + &
                          (1.0-wt_end) * CS%p_surf_prev(i,j)
          CS%p_surf_begin(i,j) = wt_beg * forces%p_surf(i,j) + &
                          (1.0-wt_beg) * CS%p_surf_prev(i,j)
        enddo ; enddo
      endif

      call step_MOM_dynamics(forces, CS%p_surf_begin, CS%p_surf_end, dt, &
                             dt_tradv_here, bbl_time_int, CS, &
                             Time_local, Waves=Waves)

      !===========================================================================
      ! This is the start of the tracer advection part of the algorithm.
      if (tradv_does_span_coupling .or. .not.do_thermo) then
        do_advection = ((CS%t_dyn_rel_adv + 0.5*dt > dt_tr_adv) .or. &
                        (CS%t_dyn_rel_thermo + 0.5*dt > dt_therm))
      else
        do_advection = ((MOD(n,ntastep) == 0) .or. (n==n_max))
      endif

      if (do_advection) then ! Do advective transport and lateral tracer mixing.
        call step_MOM_tracer_dyn(CS, G, GV, US, h, Time_local)
        if (CS%diabatic_first .and. abs(CS%t_dyn_rel_thermo) > 1e-6*dt) call MOM_error(FATAL, &
                "step_MOM: Mismatch between the dynamics and diabatic times "//&
                "with DIABATIC_FIRST.")
      endif
    endif ! end of (do_dyn)

    !===========================================================================
    ! This is the second place where the diabatic processes and remapping could occur.
    if (thermo_does_span_coupling .or. .not.do_dyn) then
      do_diabatic = (do_thermo .and. (CS%t_dyn_rel_thermo + 0.5*dt > dt_therm))
    else
      do_diabatic = (do_thermo .and. ((MOD(n,ntstep) == 0) .or. (n==n_max)))
    endif
    if ((CS%t_dyn_rel_adv==0.0) .and. (.not.CS%diabatic_first) .and. do_diabatic) then

      dtdia = CS%t_dyn_rel_thermo
      ! If the MOM6 dynamic and thermodynamic time stepping is being orchestrated
      ! by the coupler, the value of diabatic_first does not matter.
      if ((CS%t_dyn_rel_thermo==0.0) .and. .not.do_dyn) dtdia = dt

      if (CS%thermo_spans_coupling .and. (CS%dt_therm > 1.5*cycle_time) .and. &
          (abs(dt_therm - dtdia) > 1e-6*dt_therm)) then
        call MOM_error(FATAL, "step_MOM: Mismatch between dt_therm and dtdia "//&
                       "before call to diabatic.")
      endif

      ! If necessary, temporarily reset CS%Time to the center of the period covered
      ! by the call to step_MOM_thermo, noting that they end at the same time.
      if (dtdia > dt) &
        CS%Time = CS%Time - real_to_time(0.5*(dtdia-dt), unscale=US%T_to_s)

      ! Apply diabatic forcing, do mixing, and regrid.
      call step_MOM_thermo(CS, G, GV, US, u, v, h, CS%tv, fluxes, dtdia, &
                           Time_local, .false., Waves=Waves)
      if ( CS%use_ALE_algorithm ) &
        call ALE_regridding_and_remapping(CS, G, GV, US, u, v, h, CS%tv, dtdia, Time_local)
      call post_diabatic_halo_updates(CS, G, GV, US, u, v, h, CS%tv)
      CS%time_in_thermo_cycle = CS%time_in_thermo_cycle + dtdia

      if ((CS%t_dyn_rel_thermo==0.0) .and. .not.do_dyn) then
        ! The diabatic processes are now ahead of the dynamics by dtdia.
        CS%t_dyn_rel_thermo = -dtdia
      else ! The diabatic processes and the dynamics are synchronized.
        CS%t_dyn_rel_thermo = 0.0
      endif

      ! Reset CS%Time to its previous value.
      if (dtdia > dt) &
        CS%Time = Time_start + real_to_time(rel_time - 0.5*dt, unscale=US%T_to_s)
    endif

    if (do_dyn) then
      call cpu_clock_begin(id_clock_dynamics)
      ! Determining the time-average sea surface height is part of the algorithm.
      ! This may be eta_av if Boussinesq, or need to be diagnosed if not.
      CS%time_in_cycle = CS%time_in_cycle + dt
      call find_eta(h, CS%tv, G, GV, US, ssh, CS%eta_av_bc, dZref=G%Z_ref)
      do j=js,je ; do i=is,ie
        CS%ssh_rint(i,j) = CS%ssh_rint(i,j) + dt*ssh(i,j)
      enddo ; enddo
      if (CS%IDs%id_ssh_inst > 0) then
        call enable_averages(dt, Time_local, CS%diag)
        call post_data(CS%IDs%id_ssh_inst, ssh, CS%diag)
        call disable_averaging(CS%diag)
      endif
      call cpu_clock_end(id_clock_dynamics)
    endif

    !===========================================================================
    ! Calculate diagnostics at the end of the time step if the state is self-consistent.
    if (MOM_state_is_synchronized(CS)) then
    !### Perhaps this should be if (CS%t_dyn_rel_thermo == 0.0)
      call cpu_clock_begin(id_clock_other) ; call cpu_clock_begin(id_clock_diagnostics)
      ! Diagnostics that require the complete state to be up-to-date can be calculated.

      call enable_averages(CS%t_dyn_rel_diag, Time_local, CS%diag)
      call calculate_diagnostic_fields(u, v, h, CS%uh, CS%vh, CS%tv, CS%ADp,  &
                          CS%CDp, p_surf, CS%t_dyn_rel_diag, CS%diag_pre_sync,&
                          G, GV, US, CS%diagnostics_CSp)
      call post_tracer_diagnostics_at_sync(CS%Tracer_reg, h, CS%diag_pre_sync, CS%diag, G, GV, CS%t_dyn_rel_diag)
      call diag_copy_diag_to_storage(CS%diag_pre_sync, h, CS%diag)
      if (showCallTree) call callTree_waypoint("finished calculate_diagnostic_fields (step_MOM)")
      call disable_averaging(CS%diag)
      CS%t_dyn_rel_diag = 0.0

      call cpu_clock_end(id_clock_diagnostics) ; call cpu_clock_end(id_clock_other)
    endif

    if (do_dyn .and. .not.CS%count_calls) CS%nstep_tot = CS%nstep_tot + 1
    if (showCallTree) call callTree_leave("DT cycles (step_MOM)")

  enddo ! complete the n loop

  if (CS%count_calls .and. cycle_start) CS%nstep_tot = CS%nstep_tot + 1

  call cpu_clock_begin(id_clock_other)

  if (CS%time_in_cycle > 0.0) then
    I_wt_ssh = 1.0/CS%time_in_cycle
    do j=js,je ; do i=is,ie
      ssh(i,j) = CS%ssh_rint(i,j)*I_wt_ssh
      CS%ave_ssh_ibc(i,j) = ssh(i,j)
    enddo ; enddo
    if (associated(CS%HA_CSp)) call HA_accum('ssh', ssh, Time_local, G, CS%HA_CSp)
    if (do_dyn) then
      call adjust_ssh_for_p_atm(CS%tv, G, GV, US, CS%ave_ssh_ibc, forces%p_surf_SSH, &
                                CS%calc_rho_for_sea_lev)
    elseif (do_thermo) then
      call adjust_ssh_for_p_atm(CS%tv, G, GV, US, CS%ave_ssh_ibc, fluxes%p_surf_SSH, &
                                CS%calc_rho_for_sea_lev)
    endif
  endif

  if (do_dyn .and. CS%interp_p_surf) then ; do j=jsd,jed ; do i=isd,ied
    CS%p_surf_prev(i,j) = forces%p_surf(i,j)
  enddo ; enddo ; endif

  if (CS%ensemble_ocean) then
    ! store ensemble vector in odaCS
    call set_prior_tracer(CS%Time, G, GV, CS%h, CS%tv, CS%odaCS)
    ! call DA interface
    call oda(CS%Time,CS%odaCS)
    ! update the time for the next analysis step if needed
    call set_analysis_time(CS%Time,CS%odaCS)
  endif

  if (showCallTree) call callTree_waypoint("calling extract_surface_state (step_MOM)")
  ! NOTE: sfc_state uses input indexing, since it is also used by drivers.
  call extract_surface_state(CS, sfc_state)

  ! Do diagnostics that only occur at the end of a complete forcing step.
  if (cycle_end) then
    if (showCallTree) call callTree_waypoint("Do cycle end diagnostics (step_MOM)")
    if (CS%rotate_index) then
      allocate(sfc_state_diag)
      call rotate_surface_state(sfc_state, sfc_state_diag, G, turns)
    else
      sfc_state_diag => sfc_state
    endif

    call cpu_clock_begin(id_clock_diagnostics)
    if (CS%time_in_cycle > 0.0) then
      call enable_averages(CS%time_in_cycle, Time_local, CS%diag)
      call post_surface_dyn_diags(CS%sfc_IDs, G, CS%diag, sfc_state_diag, ssh)
    endif
    if (CS%time_in_thermo_cycle > 0.0) then
      call enable_averages(CS%time_in_thermo_cycle, Time_local, CS%diag)
      call post_surface_thermo_diags(CS%sfc_IDs, G, GV, US, CS%diag, CS%time_in_thermo_cycle, &
                                     sfc_state_diag, CS%tv, ssh, CS%ave_ssh_ibc)
    endif
    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_diagnostics)
    if (CS%rotate_index) then
      call deallocate_surface_state(sfc_state_diag)
    endif
    if (showCallTree) call callTree_waypoint("Done with end cycle diagnostics (step_MOM)")
  endif

  ! Accumulate the surface fluxes for assessing conservation
  if (do_thermo .and. fluxes%fluxes_used) &
    call accumulate_net_input(fluxes, sfc_state, CS%tv, fluxes%dt_buoy_accum, &
                              G, US, CS%sum_output_CSp)

  if (MOM_state_is_synchronized(CS)) &
    call write_energy(CS%u, CS%v, CS%h, CS%tv, Time_local, CS%nstep_tot, &
                      G, GV, US, CS%sum_output_CSp, CS%tracer_flow_CSp, &
                      dt_forcing=real_to_time(time_interval, unscale=US%T_to_s) )

  call cpu_clock_end(id_clock_other)

  ! De-rotate fluxes and copy back to the input, since they can be changed.
  if (CS%rotate_index) then
    call rotate_forcing(fluxes, fluxes_in, -turns)
    call rotate_mech_forcing(forces, -turns, forces_in)
    call deallocate_mech_forcing(forces)
    deallocate(forces)
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif

  if (showCallTree) call callTree_leave("step_MOM()")
  call cpu_clock_end(id_clock_ocean)

end procedure step_MOM
module procedure step_MOM_dynamics
  type(ocean_grid_type),   pointer :: G => NULL()  ! pointer to a structure containing
  type(verticalGrid_type), pointer :: GV => NULL() ! Pointer to the vertical grid structure
  type(unit_scale_type),   pointer :: US => NULL() ! Pointer to a structure containing
  type(MOM_diag_IDs), pointer :: IDs => NULL() ! A structure with the diagnostic IDs.
  real, dimension(:,:,:), pointer :: &
    u => NULL(), & ! u : zonal velocity component [L T-1 ~> m s-1]
    v => NULL(), & ! v : meridional velocity component [L T-1 ~> m s-1]
    h => NULL()    ! h : layer thickness [H ~> m or kg m-2]
  type(time_type) :: Time_end_diag ! End time of a diagnostic segment, as a time type
  logical :: calc_dtbt  ! Indicates whether the dynamically adjusted
  logical :: showCallTree
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  G => CS%G ; GV => CS%GV ; US => CS%US ; IDs => CS%IDs
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  u => CS%u ; v => CS%v ; h => CS%h
  showCallTree = callTree_showQuery()

  call cpu_clock_begin(id_clock_dynamics)
  call cpu_clock_begin(id_clock_stoch)
  if (CS%use_stochastic_EOS) call MOM_stoch_eos_run(G, u, v, dt, Time_local, CS%stoch_eos_CS)
  call cpu_clock_end(id_clock_stoch)
  call cpu_clock_begin(id_clock_varT)
  if (CS%use_stochastic_EOS) then
    call MOM_calc_varT(G, GV, US, h, CS%tv, CS%stoch_eos_CS, dt)
    if (associated(CS%tv%varT)) call pass_var(CS%tv%varT, G%Domain, clock=id_clock_pass, halo=1)
  endif
  call cpu_clock_end(id_clock_varT)

  if ((CS%t_dyn_rel_adv == 0.0) .and. CS%thickness_diffuse_first .and. &
      (CS%thickness_diffuse .or. CS%interface_filter)) then

    Time_end_diag = Time_local + real_to_time(dt_tr_adv - dt, unscale=US%T_to_s)
    call enable_averages(dt_tr_adv, Time_end_diag, CS%diag)
    if (CS%thickness_diffuse) then
      call cpu_clock_begin(id_clock_thick_diff)
      if (CS%VarMix%use_variable_mixing) &
        call calc_slope_functions(h, CS%tv, dt, G, GV, US, CS%VarMix, OBC=CS%OBC)
      call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dt_tr_adv, G, GV, US, &
                             CS%MEKE, CS%VarMix, CS%CDp, CS%thickness_diffuse_CSp, &
                             CS%stoch_CS, u, v)
      call cpu_clock_end(id_clock_thick_diff)
      call pass_var(h, G%Domain, clock=id_clock_pass, halo=CS%dyn_h_stencil)
      if (showCallTree) call callTree_waypoint("finished thickness_diffuse_first (step_MOM)")
    endif

    if (CS%interface_filter) then
      if (allocated(CS%tv%SpV_avg)) call pass_var(CS%tv%SpV_avg, G%Domain, clock=id_clock_pass)
      CS%tv%valid_SpV_halo = min(G%Domain%nihalo, G%Domain%njhalo)
      call cpu_clock_begin(id_clock_int_filter)
      call interface_filter(h, CS%uhtr, CS%vhtr, CS%tv, dt_tr_adv, G, GV, US, &
                            CS%CDp, CS%interface_filter_CSp)
      call cpu_clock_end(id_clock_int_filter)
      call pass_var(h, G%Domain, clock=id_clock_pass, halo=CS%dyn_h_stencil)
      if (showCallTree) call callTree_waypoint("finished interface_filter_first (step_MOM)")
    endif

    call disable_averaging(CS%diag)
    ! Whenever thickness changes let the diag manager know, target grids
    ! for vertical remapping may need to be regenerated.
    call diag_update_remap_grids(CS%diag)
  endif

  ! Update porous barrier fractional cell metrics
  if (CS%use_porbar) then
    call enable_averages(dt, Time_local, CS%diag)
    call porous_widths_layer(h, CS%tv, G, GV, US, CS%pbv, CS%por_bar_CS)
    call disable_averaging(CS%diag)
    call pass_vector(CS%pbv%por_face_areaU, CS%pbv%por_face_areaV, &
                     G%Domain, direction=To_All+SCALAR_PAIR, clock=id_clock_pass, halo=CS%cont_stencil)
  endif

  ! The bottom boundary layer properties need to be recalculated.
  if (bbl_time_int > 0.0) then
    Time_end_diag = Time_local + real_to_time(bbl_time_int - dt, unscale=US%T_to_s)
    call enable_averages(bbl_time_int, Time_end_diag, CS%diag)
    ! Calculate the BBL properties and store them inside visc (u,h).
    call cpu_clock_begin(id_clock_BBL_visc)
    call set_viscous_BBL(CS%u, CS%v, CS%h, CS%tv, CS%visc, G, GV, US, CS%set_visc_CSp, CS%pbv)
    call cpu_clock_end(id_clock_BBL_visc)
    if (showCallTree) call callTree_wayPoint("done with set_viscous_BBL (step_MOM)")
    call disable_averaging(CS%diag)
  endif

  !OBC segment data update for some fields can be less frequent than others
  if (associated(CS%OBC)) then
    CS%OBC%update_OBC_seg_data = .false.
    if (CS%dt_obc_seg_period == 0.0) CS%OBC%update_OBC_seg_data = .true.
    if (CS%dt_obc_seg_period > 0.0) then
      if (Time_local >= CS%dt_obc_seg_time) then
        CS%OBC%update_OBC_seg_data = .true.
        CS%dt_obc_seg_time = CS%dt_obc_seg_time + CS%dt_obc_seg_interval
      endif
    endif
  endif
  ! if (CS%debug_OBCs .and. associated(CS%OBC)) call chksum_OBC_segments(CS%OBC, G, GV, US, 3)

  if (CS%do_dynamics .and. CS%split) then !--------------------------- start SPLIT
    ! This section uses a split time stepping scheme for the dynamic equations,
    ! basically the stacked shallow water equations with viscosity.

    calc_dtbt = .false.
    if (CS%dtbt_reset_period == 0.0) calc_dtbt = .true.
    if (CS%dtbt_reset_period > 0.0) then
      if (Time_local >= CS%dtbt_reset_time) then  !### Change >= to > here.
        calc_dtbt = .true.
        CS%dtbt_reset_time = CS%dtbt_reset_time + CS%dtbt_reset_interval
      endif
    endif

    if (CS%use_alt_split) then
      call step_MOM_dyn_split_RK2b(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
                  p_surf_begin, p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
                  CS%eta_av_bc, G, GV, US, CS%dyn_split_RK2b_CSp, calc_dtbt, CS%VarMix, &
                  CS%MEKE, CS%thickness_diffuse_CSp, CS%pbv, waves=waves)
    else
      call step_MOM_dyn_split_RK2(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
                  p_surf_begin, p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
                  CS%eta_av_bc, G, GV, US, CS%dyn_split_RK2_CSp, calc_dtbt, CS%VarMix, &
                  CS%MEKE, CS%thickness_diffuse_CSp, CS%pbv, CS%stoch_CS, waves=waves)
    endif
    if (showCallTree) call callTree_waypoint("finished step_MOM_dyn_split (step_MOM)")

  elseif (CS%do_dynamics) then ! ------------------------------------ not SPLIT
    !   This section uses an unsplit stepping scheme for the dynamic
    ! equations; basically the stacked shallow water equations with viscosity.
    ! Because the time step is limited by CFL restrictions on the external
    ! gravity waves, the unsplit is usually much less efficient that the split
    ! approaches. But because of its simplicity, the unsplit method is very
    ! useful for debugging purposes.

    if (CS%use_RK2) then
      call step_MOM_dyn_unsplit_RK2(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
               p_surf_begin, p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
               CS%eta_av_bc, G, GV, US, CS%dyn_unsplit_RK2_CSp, CS%VarMix, CS%MEKE, CS%pbv, &
               CS%stoch_CS)
    else
      call step_MOM_dyn_unsplit(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
               p_surf_begin, p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
               CS%eta_av_bc, G, GV, US, CS%dyn_unsplit_CSp, CS%VarMix, CS%MEKE, CS%pbv, &
               CS%stoch_CS, Waves=Waves)
    endif
    if (showCallTree) call callTree_waypoint("finished step_MOM_dyn_unsplit (step_MOM)")

  endif ! -------------------------------------------------- end SPLIT

  if (CS%use_particles .and. CS%do_dynamics .and. (.not. CS%use_uh_particles)) then
    if (CS%thickness_diffuse_first) call MOM_error(WARNING,"particles_run: "//&
      "Thickness_diffuse_first is true and use_uh_particles is false. "//&
      "This is usually a bad combination.")
    !Run particles using unweighted velocity
    call particles_run(CS%particles, Time_local, CS%u, CS%v, CS%h, &
                       CS%tv, dt, CS%use_uh_particles)
    call particles_to_z_space(CS%particles, h)
  endif

  ! Update the model's current to reflect wind-wave growth
  if (Waves%Stokes_DDT .and. (.not.Waves%Passive_Stokes_DDT)) then
    do J=jsq,jeq ; do i=is,ie
      v(i,J,:) = v(i,J,:) + Waves%ddt_us_y(i,J,:)*dt
    enddo ; enddo
    do j=js,je ; do I=isq,ieq
      u(I,j,:) = u(I,j,:) + Waves%ddt_us_x(I,j,:)*dt
    enddo ; enddo
    call pass_vector(u, v, G%Domain)
  endif
  ! Added an additional output to track Stokes drift time tendency.
  ! It is mostly for debugging, and perhaps doesn't need to hang
  ! around permanently.
  if (Waves%Stokes_DDT .and. (Waves%id_3dstokes_y_from_ddt>0)) then
    do J=jsq,jeq ; do i=is,ie
      Waves%us_y_from_ddt(i,J,:) = Waves%us_y_from_ddt(i,J,:) + Waves%ddt_us_y(i,J,:)*dt
    enddo ; enddo
  endif
  if (Waves%Stokes_DDT .and. (Waves%id_3dstokes_x_from_ddt>0)) then
    do j=js,je ; do I=isq,ieq
      Waves%us_x_from_ddt(I,j,:) = Waves%us_x_from_ddt(I,j,:) + Waves%ddt_us_x(I,j,:)*dt
    enddo ; enddo
  endif


  if ((CS%thickness_diffuse .or. CS%interface_filter) .and. &
      .not.CS%thickness_diffuse_first) then

    if (CS%debug) call hchksum(h,"Pre-thickness_diffuse h", G%HI, haloshift=0, unscale=GV%H_to_MKS)

    if (CS%thickness_diffuse) then
      call cpu_clock_begin(id_clock_thick_diff)
      if (CS%VarMix%use_variable_mixing) &
        call calc_slope_functions(h, CS%tv, dt, G, GV, US, CS%VarMix, OBC=CS%OBC)
      call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dt, G, GV, US, &
                             CS%MEKE, CS%VarMix, CS%CDp, CS%thickness_diffuse_CSp, &
                             CS%stoch_CS, u, v)
      call cpu_clock_end(id_clock_thick_diff)
      call pass_var(h, G%Domain, clock=id_clock_pass, halo=CS%dyn_h_stencil)
      if (CS%debug) call hchksum(h,"Post-thickness_diffuse h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      if (showCallTree) call callTree_waypoint("finished thickness_diffuse (step_MOM)")
    endif

    if (CS%interface_filter) then
      if (allocated(CS%tv%SpV_avg)) call pass_var(CS%tv%SpV_avg, G%Domain, clock=id_clock_pass)
      CS%tv%valid_SpV_halo = min(G%Domain%nihalo, G%Domain%njhalo)
      call cpu_clock_begin(id_clock_int_filter)
      if (CS%interface_filter_dt_bug) then
        call interface_filter(h, CS%uhtr, CS%vhtr, CS%tv, dt_tr_adv, G, GV, US, &
                              CS%CDp, CS%interface_filter_CSp)
      else
        call interface_filter(h, CS%uhtr, CS%vhtr, CS%tv, dt, G, GV, US, &
                              CS%CDp, CS%interface_filter_CSp)
      endif
      call cpu_clock_end(id_clock_int_filter)
      call pass_var(h, G%Domain, clock=id_clock_pass, halo=CS%dyn_h_stencil)
      if (showCallTree) call callTree_waypoint("finished interface_filter (step_MOM)")
    endif
  endif

  ! apply the submesoscale mixed layer restratification parameterization
  if (CS%mixedlayer_restrat) then
    if (CS%debug) then
      call hchksum(h,"Pre-mixedlayer_restrat h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      call uvchksum("Pre-mixedlayer_restrat uhtr", &
                    CS%uhtr, CS%vhtr, G%HI, haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
    endif
    call cpu_clock_begin(id_clock_ml_restrat)
    call mixedlayer_restrat(h, CS%uhtr, CS%vhtr, CS%tv, forces, dt, CS%visc%MLD, CS%visc%h_ML, &
                            CS%visc%sfc_buoy_flx, CS%VarMix, G, GV, US, CS%mixedlayer_restrat_CSp)
    call cpu_clock_end(id_clock_ml_restrat)
    call pass_var(h, G%Domain, clock=id_clock_pass, halo=CS%dyn_h_stencil)
    if (CS%debug) then
      call hchksum(h,"Post-mixedlayer_restrat h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      call uvchksum("Post-mixedlayer_restrat [uv]htr", &
                    CS%uhtr, CS%vhtr, G%HI, haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
    endif
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  if (CS%useMEKE .and. CS%MEKE_in_dynamics) then
    call step_forward_MEKE(CS%MEKE, h, CS%VarMix%SN_u, CS%VarMix%SN_v, &
                           CS%visc, dt, G, GV, US, CS%MEKE_CSp, CS%uhtr, CS%vhtr, &
                           CS%u, CS%v, CS%tv, Time_local)
  endif
  call disable_averaging(CS%diag)

  ! Advance the dynamics time by dt.
  CS%t_dyn_rel_adv = CS%t_dyn_rel_adv + dt

  if (CS%use_particles .and. CS%do_dynamics .and. CS%use_uh_particles .and. &
      CS%uh_particles_bug) then
    ! Run particles using thickness-weighted velocity
    call particles_run(CS%particles, Time_local, CS%uhtr, CS%vhtr, CS%h, &
                       CS%tv, CS%t_dyn_rel_adv, CS%use_uh_particles)
  endif

  CS%n_dyn_steps_in_adv = CS%n_dyn_steps_in_adv + 1
  if (CS%alternate_first_direction) then
    call set_first_direction(G, MODULO(G%first_direction+1,2))
    CS%first_dir_restart = real(G%first_direction)
  elseif (CS%use_particles .and. CS%do_dynamics .and. (.not.CS%use_uh_particles)) then
    call particles_to_k_space(CS%particles, h)
  endif
  CS%t_dyn_rel_thermo = CS%t_dyn_rel_thermo + dt
  if (abs(CS%t_dyn_rel_thermo) < 1e-6*dt) CS%t_dyn_rel_thermo = 0.0
  CS%t_dyn_rel_diag = CS%t_dyn_rel_diag + dt

  call cpu_clock_end(id_clock_dynamics)

  call cpu_clock_begin(id_clock_other) ; call cpu_clock_begin(id_clock_diagnostics)
  call enable_averages(dt, Time_local, CS%diag)
  ! These diagnostics are available after every time dynamics step.
  if (IDs%id_u > 0) call post_data(IDs%id_u, u, CS%diag)
  if (IDs%id_v > 0) call post_data(IDs%id_v, v, CS%diag)
  if (IDs%id_h > 0) call post_data(IDs%id_h, h, CS%diag)
  if (CS%use_stochastic_EOS) call post_stoch_EOS_diags(CS%stoch_eos_CS, CS%tv, CS%diag)
  call disable_averaging(CS%diag)
  call cpu_clock_end(id_clock_diagnostics) ; call cpu_clock_end(id_clock_other)

end procedure step_MOM_dynamics
module procedure step_MOM_tracer_dyn
  type(group_pass_type) :: pass_T_S
  integer :: halo_sz ! The size of a halo where data must be valid.
  logical :: x_first ! If true, advect tracers first in the x-direction, then y.
  logical :: showCallTree
  showCallTree = callTree_showQuery()

  if (CS%debug) then
    call cpu_clock_begin(id_clock_other)
    call hchksum(h,"Pre-advection h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    call uvchksum("Pre-advection uhtr", CS%uhtr, CS%vhtr, G%HI, &
                  haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
    if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Pre-advection T", G%HI, haloshift=1, unscale=US%C_to_degC)
    if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Pre-advection S", G%HI, haloshift=1, unscale=US%S_to_ppt)
    if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, "Pre-advection frazil", G%HI, haloshift=0, &
                                               unscale=US%Q_to_J_kg*US%RZ_to_kg_m2)
    if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, &
                   "Pre-advection salt deficit", G%HI, haloshift=0, unscale=US%S_to_ppt*US%RZ_to_kg_m2)
  ! call MOM_thermo_chksum("Pre-advection ", CS%tv, G, US)
    call cpu_clock_end(id_clock_other)
  endif

  call cpu_clock_begin(id_clock_thermo) ; call cpu_clock_begin(id_clock_tracer)
  call enable_averages(CS%t_dyn_rel_adv, Time_local, CS%diag)

  if (CS%use_particles .and. CS%use_uh_particles .and. (.not. CS%uh_particles_bug)) then
    ! Run particles using thickness-weighted velocity
    call particles_run(CS%particles, Time_local, CS%uhtr, CS%vhtr, CS%h, &
                       CS%tv, CS%t_dyn_rel_adv, CS%use_uh_particles)
  endif


  if (CS%alternate_first_direction) then
    ! This calculation of the value of G%first_direction from the start of the accumulation of
    ! mass transports for use by the tracers is the equivalent to adding 2*n_dyn_steps before
    ! subtracting n_dyn_steps so that the mod will be taken of a non-negative number.
    x_first = (MODULO(G%first_direction+CS%n_dyn_steps_in_adv,2) == 0)
  else
    x_first = (MODULO(G%first_direction,2) == 0)
  endif
  if (CS%debug) call MOM_tracer_chksum("Pre-advect ", CS%tracer_Reg, G)
  call advect_tracer(h, CS%uhtr, CS%vhtr, CS%OBC, CS%t_dyn_rel_adv, G, GV, US, &
                     CS%tracer_adv_CSp, CS%tracer_Reg, x_first_in=x_first)
  if (CS%debug) call MOM_tracer_chksum("Post-advect ", CS%tracer_Reg, G)
  call tracer_hordiff(h, CS%t_dyn_rel_adv, CS%MEKE, CS%VarMix, CS%visc, G, GV, US, &
                      CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
  if (CS%debug) call MOM_tracer_chksum("Post-diffuse ", CS%tracer_Reg, G)
  if (showCallTree) call callTree_waypoint("finished tracer advection/diffusion (step_MOM)")
  if (associated(CS%OBC)) then
    call pass_vector(CS%uhtr, CS%vhtr, G%Domain)
    call update_segment_tracer_reservoirs(G, GV, CS%uhtr, CS%vhtr, h, CS%OBC, &
                     CS%tracer_Reg)
  endif
  call cpu_clock_end(id_clock_tracer) ; call cpu_clock_end(id_clock_thermo)

  call cpu_clock_begin(id_clock_other) ; call cpu_clock_begin(id_clock_diagnostics)
  call post_transport_diagnostics(G, GV, US, CS%uhtr, CS%vhtr, h, CS%transport_IDs, &
           CS%diag_pre_dyn, CS%diag, CS%t_dyn_rel_adv, CS%tracer_reg)
  ! Rebuild the remap grids now that we've posted the fields which rely on thicknesses
  ! from before the dynamics calls
  call diag_update_remap_grids(CS%diag)

  call disable_averaging(CS%diag)
  call cpu_clock_end(id_clock_diagnostics) ; call cpu_clock_end(id_clock_other)

  ! Reset the accumulated transports to 0 and record that the dynamics
  ! and advective times now agree.
  call cpu_clock_begin(id_clock_thermo) ; call cpu_clock_begin(id_clock_tracer)
  CS%uhtr(:,:,:) = 0.0
  CS%vhtr(:,:,:) = 0.0
  CS%n_dyn_steps_in_adv = 0
  CS%t_dyn_rel_adv = 0.0
  call cpu_clock_end(id_clock_tracer) ; call cpu_clock_end(id_clock_thermo)

  if (CS%useMEKE .and. (.not. CS%MEKE_in_dynamics)) then
    call step_forward_MEKE(CS%MEKE, h, CS%VarMix%SN_u, CS%VarMix%SN_v, &
                           CS%visc, CS%t_dyn_rel_adv, G, GV, US, CS%MEKE_CSp, CS%uhtr, CS%vhtr, &
                           CS%u, CS%v, CS%tv, Time_local)
  endif

  if (associated(CS%tv%T)) then
    call extract_diabatic_member(CS%diabatic_CSp, diabatic_halo=halo_sz)
    ! The bottom boundary layer calculation may need halo values of SpV_avg, including the corners.
    if (allocated(CS%tv%SpV_avg)) halo_sz = max(halo_sz, 1)
    if (halo_sz > 0) then
      call create_group_pass(pass_T_S, CS%tv%T, G%Domain, To_All, halo=halo_sz)
      call create_group_pass(pass_T_S, CS%tv%S, G%Domain, To_All, halo=halo_sz)
      call do_group_pass(pass_T_S, G%Domain, clock=id_clock_pass)
    elseif (CS%diabatic_first) then
      ! Temperature and salinity need halo updates because they will be used
      ! in the dynamics before they are changed again.
      call create_group_pass(pass_T_S, CS%tv%T, G%Domain, To_All+Omit_Corners, halo=1)
      call create_group_pass(pass_T_S, CS%tv%S, G%Domain, To_All+Omit_Corners, halo=1)
      call do_group_pass(pass_T_S, G%Domain, clock=id_clock_pass)
      halo_sz = 1
    endif

    ! Update derived thermodynamic quantities.
    if (allocated(CS%tv%SpV_avg)) then
      call calc_derived_thermo(CS%tv, h, G, GV, US, halo=halo_sz, debug=CS%debug)
    endif
  endif

  CS%preadv_h_stored = .false.

end procedure step_MOM_tracer_dyn
module procedure step_MOM_thermo
  logical :: debug_redundant ! If true, check redundant values on PE boundaries when debugging.
  logical :: showCallTree
  type(group_pass_type) :: pass_T_S
  integer :: dynamics_stencil  ! The computational stencil for the calculations
  integer :: halo_sz ! The size of a halo where data must be valid.
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("step_MOM_thermo(), MOM.F90")
  if (CS%debug) call query_debugging_checks(do_redundant=debug_redundant)

  call enable_averages(dtdia, Time_end_thermo, CS%diag)

  if (associated(CS%odaCS)) then
    if (CS%debug) then
      call MOM_thermo_chksum("Pre-oda ", tv, G, US, haloshift=0)
    endif
    call apply_oda_tracer_increments(dtdia, Time_end_thermo, G, GV, tv, h, CS%odaCS)
    if (CS%debug) then
      call MOM_thermo_chksum("Post-oda ", tv, G, US, haloshift=0)
    endif
  endif

  if (associated(fluxes%p_surf) .or. associated(fluxes%p_surf_full)) then
    call extract_diabatic_member(CS%diabatic_CSp, diabatic_halo=halo_sz)
    if (halo_sz > 0) then
      if (associated(fluxes%p_surf_full)) &
        call pass_var(fluxes%p_surf_full, G%Domain, &
                      clock=id_clock_pass, halo=halo_sz, complete=.not.associated(fluxes%p_surf))
      call pass_var(fluxes%p_surf, G%Domain, clock=id_clock_pass, halo=halo_sz, complete=.true.)
    endif
  endif

  if (update_BBL) then
    !   Calculate the BBL properties and store them inside visc (u,h).
    ! This is here so that CS%visc is updated before diabatic() when
    ! DIABATIC_FIRST=True. Otherwise diabatic() is called after the dynamics
    ! and set_viscous_BBL is called as a part of the dynamic stepping.
    call cpu_clock_begin(id_clock_BBL_visc)
    !update porous barrier fractional cell metrics
    if (CS%use_porbar) then
      call porous_widths_interface(h, CS%tv, G, GV, US, CS%pbv, CS%por_bar_CS)
      call pass_vector(CS%pbv%por_layer_widthU, CS%pbv%por_layer_widthV, &
                      G%Domain, direction=To_ALL+SCALAR_PAIR, clock=id_clock_pass, halo=CS%cont_stencil)
    endif
    call set_viscous_BBL(u, v, h, tv, CS%visc, G, GV, US, CS%set_visc_CSp, CS%pbv)
    call cpu_clock_end(id_clock_BBL_visc)
    if (showCallTree) call callTree_wayPoint("done with set_viscous_BBL (step_MOM_thermo)")
  endif

  call cpu_clock_begin(id_clock_thermo)
  if (.not.CS%adiabatic) then
    if (CS%debug) then
      call uvchksum("Pre-diabatic [uv]", u, v, G%HI, haloshift=2, unscale=US%L_T_to_m_s)
      call hchksum(h,"Pre-diabatic h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      call uvchksum("Pre-diabatic [uv]h", CS%uhtr, CS%vhtr, G%HI, &
                    haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
    ! call MOM_state_chksum("Pre-diabatic ", u, v, h, CS%uhtr, CS%vhtr, G, GV, vel_scale=1.0)
      call MOM_thermo_chksum("Pre-diabatic ", tv, G, US, haloshift=0)
      if (debug_redundant) &
        call check_redundant("Pre-diabatic ", u, v, G, unscale=US%L_T_to_m_s)
      call MOM_forcing_chksum("Pre-diabatic", fluxes, G, US, haloshift=0)
    endif

    call cpu_clock_begin(id_clock_diabatic)

    call diabatic(u, v, h, tv, CS%Hml, fluxes, CS%visc, CS%ADp, CS%CDp, dtdia, &
                  Time_end_thermo, G, GV, US, CS%diabatic_CSp, CS%stoch_CS, CS%OBC, Waves)
    fluxes%fluxes_used = .true.

    if (CS%stoch_CS%do_skeb) then
       call apply_skeb(CS%G,CS%GV,CS%stoch_CS,CS%u,CS%v,CS%h,CS%tv,dtdia,Time_end_thermo)
    endif

    if (showCallTree) call callTree_waypoint("finished diabatic (step_MOM_thermo)")

    if (CS%debug) then
      call uvchksum("Post-diabatic u", u, v, G%HI, haloshift=2, unscale=US%L_T_to_m_s)
      call hchksum(h, "Post-diabatic h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      call uvchksum("Post-diabatic [uv]h", CS%uhtr, CS%vhtr, G%HI, &
                    haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
    ! call MOM_state_chksum("Post-diabatic ", u, v, &
    !                       h, CS%uhtr, CS%vhtr, G, GV, haloshift=1)
      if (associated(tv%T)) call hchksum(tv%T, "Post-diabatic T", G%HI, haloshift=1, unscale=US%C_to_degC)
      if (associated(tv%S)) call hchksum(tv%S, "Post-diabatic S", G%HI, haloshift=1, unscale=US%S_to_ppt)
      if (associated(tv%frazil)) call hchksum(tv%frazil, "Post-diabatic frazil", G%HI, haloshift=0, &
                                              unscale=US%Q_to_J_kg*US%RZ_to_kg_m2)
      if (associated(tv%salt_deficit)) call hchksum(tv%salt_deficit, &
                               "Post-diabatic salt deficit", G%HI, haloshift=0, unscale=US%RZ_to_kg_m2)
    ! call MOM_thermo_chksum("Post-diabatic ", tv, G, US)
      if (debug_redundant) &
        call check_redundant("Post-diabatic ", u, v, G, unscale=US%L_T_to_m_s)
    endif
    call disable_averaging(CS%diag)

    call cpu_clock_end(id_clock_diabatic)
  else   ! complement of "if (.not.CS%adiabatic)"

    call cpu_clock_begin(id_clock_adiabatic)
    call adiabatic(h, tv, fluxes, dtdia, G, GV, US, CS%diabatic_CSp)
    fluxes%fluxes_used = .true.
    call cpu_clock_end(id_clock_adiabatic)

    if (associated(tv%T)) then
      dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
      call create_group_pass(pass_T_S, tv%T, G%Domain, To_All+Omit_Corners, halo=dynamics_stencil)
      call create_group_pass(pass_T_S, tv%S, G%Domain, To_All+Omit_Corners, halo=dynamics_stencil)
      call do_group_pass(pass_T_S, G%Domain, clock=id_clock_pass)
      if (CS%debug) then
        if (associated(tv%T)) call hchksum(tv%T, "Post-diabatic T", G%HI, haloshift=1, unscale=US%C_to_degC)
        if (associated(tv%S)) call hchksum(tv%S, "Post-diabatic S", G%HI, haloshift=1, unscale=US%S_to_ppt)
      endif

      ! Update derived thermodynamic quantities.
      if (allocated(tv%SpV_avg)) then
        call calc_derived_thermo(tv, h, G, GV, US, halo=dynamics_stencil, debug=CS%debug)
      endif
    endif

  endif   ! endif for the block "if (.not.CS%adiabatic)"
  call cpu_clock_end(id_clock_thermo)

  call disable_averaging(CS%diag)

! This works in general:
!  if (associated(tv%T)) &
!    call totalTandS(G%HI, h, G%areaT, tv%T, tv%S, "End of step_MOM", US, GV%H_to_mks)
! This works only if there is no rescaling being used:
!  if (associated(tv%T)) &
!    call totalTandS(G%HI, h, G%areaT, tv%T, tv%S, "End of step_MOM")

  if (showCallTree) call callTree_leave("step_MOM_thermo(), MOM.F90")

end procedure step_MOM_thermo
module procedure ALE_regridding_and_remapping
  real :: h_new(SZI_(G),SZJ_(G),SZK_(GV))      ! Layer thicknesses after regridding [H ~> m or kg m-2]
  real :: dzRegrid(SZI_(G),SZJ_(G),SZK_(GV)+1) ! The change in grid interface positions due to regridding,
  real :: h_old_u(SZIB_(G),SZJ_(G),SZK_(GV))   ! Source grid thickness at zonal
  real :: h_old_v(SZI_(G),SZJB_(G),SZK_(GV))   ! Source grid thickness at meridional
  real :: h_new_u(SZIB_(G),SZJ_(G),SZK_(GV))   ! Destination grid thickness at zonal
  real :: h_new_v(SZI_(G),SZJB_(G),SZK_(GV))   ! Destination grid thickness at meridional
  logical :: PCM_cell(SZI_(G),SZJ_(G),SZK_(GV)) ! If true, PCM remapping should be used in a cell.
  logical :: use_ice_shelf ! Needed for selecting the right ALE interface.
  logical :: debug_redundant ! If true, check redundant values on PE boundaries when debugging.
  logical :: showCallTree
  type(group_pass_type) :: pass_T_S_h
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  use_ice_shelf = .false.
  if (associated(CS%frac_shelf_h)) use_ice_shelf = .true.
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("ALE_regridding_and_remapping(), MOM.F90")
  if (CS%debug) call query_debugging_checks(do_redundant=debug_redundant)

  call cpu_clock_begin(id_clock_remap)

  ! Regridding/remapping is done here, at end of thermodynamics time step
  ! (that may comprise several dynamical time steps)
  ! The routine 'ALE_regrid' can be found in 'MOM_ALE.F90'.
  call enable_averages(dtdia, Time_end_thermo, CS%diag)

  call cpu_clock_begin(id_clock_pass)
  if (associated(tv%T)) &
    call create_group_pass(pass_T_S_h, tv%T, G%Domain, To_All+Omit_Corners, halo=1)
  if (associated(tv%S)) &
    call create_group_pass(pass_T_S_h, tv%S, G%Domain, To_All+Omit_Corners, halo=1)
  call create_group_pass(pass_T_S_h, h, G%Domain, To_All+Omit_Corners, halo=1)
  call do_group_pass(pass_T_S_h, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call preAle_tracer_diagnostics(CS%tracer_Reg, G, GV)

  if (CS%use_particles) then
    call particles_to_z_space(CS%particles, h)
  endif

  if (CS%debug) then
    call MOM_state_chksum("Pre-ALE ", u, v, h, CS%uh, CS%vh, G, GV, US, omit_corners=.true.)
    call hchksum(tv%T,"Pre-ALE T", G%HI, haloshift=1, omit_corners=.true., unscale=US%C_to_degC)
    call hchksum(tv%S,"Pre-ALE S", G%HI, haloshift=1, omit_corners=.true., unscale=US%S_to_ppt)
    if (debug_redundant) &
      call check_redundant("Pre-ALE ", u, v, G, unscale=US%L_T_to_m_s)
  endif
  call cpu_clock_begin(id_clock_ALE)

  call pre_ALE_diagnostics(G, GV, US, h, u, v, tv, CS%ALE_CSp)
  call ALE_update_regrid_weights(dtdia, CS%ALE_CSp)
  ! Do any necessary adjustments ot the state prior to remapping.
  call pre_ALE_adjustments(G, GV, US, h, tv, CS%tracer_Reg, CS%ALE_CSp, u, v)
  ! Adjust the target grids for diagnostics, in case there have been thickness adjustments.
  call diag_update_remap_grids(CS%diag)

  if (use_ice_shelf) then
    call ALE_regrid(G, GV, US, h, h_new, dzRegrid, tv, CS%ALE_CSp, CS%frac_shelf_h, PCM_cell)
  else
    call ALE_regrid(G, GV, US, h, h_new, dzRegrid, tv, CS%ALE_CSp, PCM_cell=PCM_cell)
  endif

  if (showCallTree) call callTree_waypoint("new grid generated")
  ! Remap all variables from the old grid h onto the new grid h_new
  call ALE_remap_tracers(CS%ALE_CSp, G, GV, h, h_new, CS%tracer_Reg, showCallTree, dtdia, PCM_cell)

  ! Determine the old and new grid thicknesses at velocity points.
  call ALE_remap_set_h_vel(CS%ALE_CSp, G, GV, h, h_old_u, h_old_v, CS%OBC, debug=showCallTree)
  if (CS%remap_uv_using_old_alg) then
    call ALE_remap_set_h_vel_via_dz(CS%ALE_CSp, G, GV, h_new, h_new_u, h_new_v, CS%OBC, h, dzRegrid, showCallTree)
  else
    call ALE_remap_set_h_vel(CS%ALE_CSp, G, GV, h_new, h_new_u, h_new_v, CS%OBC, debug=showCallTree)
  endif

  ! Remap the velocity components.
  call ALE_remap_velocities(CS%ALE_CSp, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, u, v, showCallTree, &
                            dtdia, allow_preserve_variance=.true.)

  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

  if (CS%remap_aux_vars) then
    if (CS%split .and. CS%use_alt_split) then
      call remap_dyn_split_RK2b_aux_vars(G, GV, CS%dyn_split_RK2b_CSp, h_old_u, h_old_v, &
                                         h_new_u, h_new_v, CS%ALE_CSp)
    elseif (CS%split) then
      call remap_dyn_split_RK2_aux_vars(G, GV, CS%dyn_split_RK2_CSp, h_old_u, h_old_v, h_new_u, h_new_v, CS%ALE_CSp)
    endif

    if (associated(CS%OBC) .or. associated(CS%visc%Kv_shear_Bu)) then
      call pass_var(h, G%Domain, complete=.false.)
      call pass_var(h_new, G%Domain, complete=.true.)
    endif

    if (associated(CS%OBC)) &
      call remap_OBC_fields(G, GV, h, h_new, CS%OBC, PCM_cell=PCM_cell)

    call remap_vertvisc_aux_vars(G, GV, CS%visc, h, h_new, CS%ALE_CSp, CS%OBC)
    if (associated(CS%visc%Kv_shear)) &
      call pass_var(CS%visc%Kv_shear, G%Domain, To_All+Omit_Corners, clock=id_clock_pass, halo=1)
  endif

  ! Replace the old grid with new one.  All remapping must be done by this point in the code.
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (showCallTree) call callTree_waypoint("finished ALE_regrid (ALE_regridding_and_remapping)")
  call cpu_clock_end(id_clock_ALE)

  ! Update derived thermodynamic quantities.
  if (allocated(CS%tv%SpV_avg)) then
    call calc_derived_thermo(CS%tv, CS%h, G, GV, US, halo=1, debug=CS%debug)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.  In non-Boussinesq mode,
  ! calc_derived_thermo needs to be called before diag_update_remap_grids.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

  call postALE_tracer_diagnostics(CS%tracer_Reg, G, GV, CS%diag, dtdia)

  if (CS%debug .and. CS%use_ALE_algorithm) then
    call MOM_state_chksum("Post-ALE ", u, v, h, CS%uh, CS%vh, G, GV, US)
    call hchksum(tv%T, "Post-ALE T", G%HI, haloshift=1, unscale=US%C_to_degC)
    call hchksum(tv%S, "Post-ALE S", G%HI, haloshift=1, unscale=US%S_to_ppt)
    if (debug_redundant) &
      call check_redundant("Post-ALE ", u, v, G, unscale=US%L_T_to_m_s)
  endif
  if (CS%debug) then
    call uvchksum("Post-ALE, Post-diabatic u", u, v, G%HI, haloshift=2, unscale=US%L_T_to_m_s)
    call hchksum(h, "Post-ALE, Post-diabatic h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    call uvchksum("Post-ALE, Post-diabatic [uv]h", CS%uhtr, CS%vhtr, G%HI, &
                  haloshift=0, unscale=GV%H_to_MKS*US%L_to_m**2)
  ! call MOM_state_chksum("Post-diabatic ", u, v, &
  !                       h, CS%uhtr, CS%vhtr, G, GV, haloshift=1)
    if (associated(tv%T)) call hchksum(tv%T, "Post-ALE, Post-diabatic T", G%HI, haloshift=1, unscale=US%C_to_degC)
    if (associated(tv%S)) call hchksum(tv%S, "Post-ALE, Post-diabatic S", G%HI, haloshift=1, unscale=US%S_to_ppt)
    if (associated(tv%frazil)) call hchksum(tv%frazil, "Post-ALE, Post-diabatic frazil", G%HI, haloshift=0, &
                                            unscale=US%Q_to_J_kg*US%RZ_to_kg_m2)
    if (associated(tv%salt_deficit)) call hchksum(tv%salt_deficit, &
                             "Post-ALE, Post-diabatic salt deficit", G%HI, haloshift=0, unscale=US%RZ_to_kg_m2)
  ! call MOM_thermo_chksum("Post-diabatic ", tv, G, US)
    if (debug_redundant) &
      call check_redundant("Post-ALE, Post-diabatic ", u, v, G, unscale=US%L_T_to_m_s)
  endif
  call disable_averaging(CS%diag)

  call cpu_clock_end(id_clock_remap)

  if (showCallTree) call callTree_leave("ALE_regridding_and_remapping(), MOM.F90")

end procedure ALE_regridding_and_remapping
module procedure post_diabatic_halo_updates
  logical :: debug_redundant ! If true, check redundant values on PE boundaries when debugging.
  logical :: showCallTree
  type(group_pass_type) :: pass_uv_T_S_h
  integer :: dynamics_stencil  ! The computational stencil for the calculations
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("post_diabatic_halo_updates, MOM.F90")
  if (CS%debug) call query_debugging_checks(do_redundant=debug_redundant)

  if (CS%use_particles) then
    call particles_to_k_space(CS%particles, h)
  endif

  dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
  call create_group_pass(pass_uv_T_S_h, u, v, G%Domain, halo=dynamics_stencil)
  if (associated(tv%T)) &
    call create_group_pass(pass_uv_T_S_h, tv%T, G%Domain, halo=dynamics_stencil)
  if (associated(tv%S)) &
    call create_group_pass(pass_uv_T_S_h, tv%S, G%Domain, halo=dynamics_stencil)
  call create_group_pass(pass_uv_T_S_h, h, G%Domain, halo=dynamics_stencil)
  call do_group_pass(pass_uv_T_S_h, G%Domain, clock=id_clock_pass)

  if (associated(tv%frazil) .and. (.not.tv%frazil_was_reset) .and. CS%vertex_shear) &
    call pass_var(tv%frazil, G%Domain, halo=1)

  ! Update derived thermodynamic quantities.
  if (allocated(tv%SpV_avg)) then
    call calc_derived_thermo(tv, h, G, GV, US, halo=dynamics_stencil, debug=CS%debug)
  endif
  if (showCallTree) call callTree_leave("post_diabatic_halo_updates, MOM.F90")
end procedure post_diabatic_halo_updates
module procedure step_offline
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
  type(unit_scale_type),      pointer :: US => NULL() ! Pointer to a structure containing
  logical :: first_iter    !< True if this is the first time step_offline has been called in a given interval
  logical :: last_iter     !< True if this is the last time step_tracer is to be called in an offline interval
  logical :: do_vertical   !< If enough time has elapsed, do the diabatic tracer sources/sinks
  logical :: adv_converged !< True if all the horizontal fluxes have been used
  real, allocatable, dimension(:,:,:) :: h_new    ! Layer thicknesses after regridding [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: dzRegrid ! The change in grid interface positions due to regridding,
  real :: dt_offline          ! The offline timestep for advection [T ~> s]
  real :: dt_offline_vertical ! The offline timestep for vertical fluxes and remapping [T ~> s]
  logical :: skip_diffusion
  type(time_type), pointer :: accumulated_time => NULL()
  type(time_type), pointer :: vertical_time => NULL()
  integer :: dynamics_stencil  ! The computational stencil for the calculations
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  real, dimension(:,:,:), pointer :: &
    uhtr => NULL(), &  ! Accumulated zonal thickness fluxes to advect tracers [H L2 ~> m3 or kg]
    vhtr => NULL(), &  ! Accumulated meridional thickness fluxes to advect tracers [H L2 ~> m3 or kg]
    eatr => NULL(), &  ! Layer entrainment rates across the interface above [H ~> m or kg m-2]
    ebtr => NULL(), &  ! Layer entrainment rates across the interface below [H ~> m or kg m-2]
    h_end => NULL()    ! Layer thicknesses at the end of a step [H ~> m or kg m-2]
  type(time_type) :: Time_end    ! End time of a segment, as a time type
  G => CS%G ; GV => CS%GV ; US => CS%US

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed

  call cpu_clock_begin(id_clock_offline_tracer)
  call extract_offline_main(CS%offline_CSp, uhtr, vhtr, eatr, ebtr, h_end, accumulated_time, &
                            vertical_time, dt_offline, dt_offline_vertical, skip_diffusion)
  Time_end = increment_date(Time_start, seconds=floor(US%T_to_s*time_interval+0.001))

  call enable_averages(time_interval, Time_end, CS%diag)

  ! Check to see if this is the first iteration of the offline interval
  first_iter = (accumulated_time == real_to_time(0.0))

  ! Check to see if vertical tracer functions should be done
  do_vertical = (first_iter .or. (accumulated_time >= vertical_time))
  if (do_vertical) vertical_time = accumulated_time + real_to_time(dt_offline_vertical, unscale=US%T_to_s)

  ! Increment the amount of time elapsed since last read and check if it's time to roll around
  accumulated_time = accumulated_time + real_to_time(time_interval, unscale=US%T_to_s)

  last_iter = (accumulated_time >= real_to_time(dt_offline, unscale=US%T_to_s))

  if (CS%use_ALE_algorithm) then
    ! If this is the first iteration in the offline timestep, then we need to read in fields and
    ! perform the main advection.
    if (first_iter) then
      call MOM_mesg("Reading in new offline fields")
      ! Read in new transport and other fields
      ! call update_transport_from_files(G, GV, CS%offline_CSp, h_end, eatr, ebtr, uhtr, vhtr, &
      !     CS%tv%T, CS%tv%S, fluxes, CS%use_ALE_algorithm)
      ! call update_transport_from_arrays(CS%offline_CSp)
      call update_offline_fields(CS%offline_CSp, G, GV, US, CS%h, fluxes, CS%use_ALE_algorithm)

      ! Apply any fluxes into the ocean
      call offline_fw_fluxes_into_ocean(G, GV, CS%offline_CSp, fluxes, CS%h)

      if (.not.CS%diabatic_first) then
        call offline_advection_ale(fluxes, Time_start, time_interval, G, GV, US, CS%offline_CSp, &
                                   id_clock_ALE, CS%h, uhtr, vhtr, converged=adv_converged)

        ! Redistribute any remaining transport
        call offline_redistribute_residual(CS%offline_CSp, G, GV, US, CS%h, uhtr, vhtr, adv_converged)

        ! Perform offline diffusion if requested
        if (.not. skip_diffusion) then
          if (CS%VarMix%use_variable_mixing) then
            call pass_var(CS%h, G%Domain)
            call calc_resoln_function(CS%h, CS%tv, G, GV, US, CS%VarMix, CS%MEKE, CS%OBC, dt_offline)
            call calc_depth_function(G, CS%VarMix)
            call calc_slope_functions(CS%h, CS%tv, dt_offline, G, GV, US, CS%VarMix, OBC=CS%OBC)
          endif
          call tracer_hordiff(CS%h, dt_offline, CS%MEKE, CS%VarMix, CS%visc, G, GV, US, &
                              CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
        endif
      endif
    endif
    ! The functions related to column physics of tracers is performed separately in ALE mode
    if (do_vertical) then
      call offline_diabatic_ale(fluxes, Time_start, Time_end, G, GV, US, CS%offline_CSp, &
                                CS%h, CS%tv, eatr, ebtr)
    endif

    ! Last thing that needs to be done is the final ALE remapping
    if (last_iter) then
      if (CS%diabatic_first) then
        call offline_advection_ale(fluxes, Time_start, time_interval, G, GV, US, CS%offline_CSp, &
                                   id_clock_ALE, CS%h, uhtr, vhtr, converged=adv_converged)

        ! Redistribute any remaining transport and perform the remaining advection
        call offline_redistribute_residual(CS%offline_CSp, G, GV, US, CS%h, uhtr, vhtr, adv_converged)
                ! Perform offline diffusion if requested
        if (.not. skip_diffusion) then
          if (CS%VarMix%use_variable_mixing) then
            call pass_var(CS%h, G%Domain)
            call calc_resoln_function(CS%h, CS%tv, G, GV, US, CS%VarMix, CS%MEKE, CS%OBC, dt_offline)
            call calc_depth_function(G, CS%VarMix)
            call calc_slope_functions(CS%h, CS%tv, dt_offline, G, GV, US, CS%VarMix, OBC=CS%OBC)
          endif
          call tracer_hordiff(CS%h, dt_offline, CS%MEKE, CS%VarMix, CS%visc, G, GV, US, &
              CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
        endif
      endif

      call MOM_mesg("Last iteration of offline interval")

      ! Apply freshwater fluxes out of the ocean
      call offline_fw_fluxes_out_ocean(G, GV, CS%offline_CSp, fluxes, CS%h)
      ! These diagnostic can be used to identify which grid points did not converge within
      ! the specified number of advection sub iterations
      call post_offline_convergence_diags(G, GV, CS%offline_CSp, CS%h, h_end, uhtr, vhtr)

      ! Call ALE one last time to make sure that tracers are remapped onto the layer thicknesses
      ! stored from the forward run
      call cpu_clock_begin(id_clock_ALE)

      ! Do any necessary adjustments ot the state prior to remapping.
      call pre_ALE_adjustments(G, GV, US, h_end, CS%tv, CS%tracer_Reg, CS%ALE_CSp)

      allocate(h_new(isd:ied, jsd:jed, nz), source=0.0)
      allocate(dzRegrid(isd:ied, jsd:jed, nz+1), source=0.0)

      ! Generate the new grid based on the tracer grid at the end of the interval.
      call ALE_regrid(G, GV, US, h_end, h_new, dzRegrid, CS%tv, CS%ALE_CSp)

      ! Remap the tracers from the previous tracer grid onto the new grid.  The thicknesses that
      ! are used are intended to ensure that in the case where transports don't quite conserve,
      ! the offline layer thicknesses do not drift too far away from the online model.
      call ALE_remap_tracers(CS%ALE_CSp, G, GV, CS%h, h_new, CS%tracer_Reg, debug=CS%debug)
      if (allocated(CS%tv%SpV_avg)) CS%tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

      ! Update the tracer grid.
      do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        CS%h(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

      deallocate(h_new, dzRegrid)

      call cpu_clock_end(id_clock_ALE)
      call pass_var(CS%h, G%Domain)
    endif
  else ! NON-ALE MODE...NOT WELL TESTED
    call MOM_error(WARNING, &
        "Offline tracer mode in non-ALE configuration has not been thoroughly tested")
    ! Note that for the layer mode case, the calls to tracer sources and sinks is embedded in
    ! main_offline_advection_layer. Warning: this may not be appropriate for tracers that
    ! exchange with the atmosphere
    if (abs(time_interval - dt_offline) > 1.0e-6*US%s_to_T) then
      call MOM_error(FATAL, &
          "For offline tracer mode in a non-ALE configuration, dt_offline must equal time_interval")
    endif
    call update_offline_fields(CS%offline_CSp, G, GV, US, CS%h, fluxes, CS%use_ALE_algorithm)
    call offline_advection_layer(fluxes, Time_start, time_interval, G, GV, US, CS%offline_CSp, &
                                 CS%h, eatr, ebtr, uhtr, vhtr)
    ! Perform offline diffusion if requested
    if (.not. skip_diffusion) then
      call tracer_hordiff(h_end, dt_offline, CS%MEKE, CS%VarMix, CS%visc, G, GV, US, &
                          CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
    endif

    CS%h = h_end

    call pass_var(CS%tv%T, G%Domain)
    call pass_var(CS%tv%S, G%Domain)
    call pass_var(CS%h, G%Domain)

  endif

  call adjust_ssh_for_p_atm(CS%tv, G, GV, US, CS%ave_ssh_ibc, forces%p_surf_SSH, &
                            CS%calc_rho_for_sea_lev)
  call extract_surface_state(CS, sfc_state)

  call disable_averaging(CS%diag)
  call pass_var(CS%tv%T, G%Domain)
  call pass_var(CS%tv%S, G%Domain)
  call pass_var(CS%h, G%Domain)

  fluxes%fluxes_used = .true.

  ! Update derived thermodynamic quantities.
  if (allocated(CS%tv%SpV_avg)) then
    dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
    call calc_derived_thermo(CS%tv, CS%h, G, GV, US, halo=dynamics_stencil)
  endif

  if (last_iter) then
    accumulated_time = real_to_time(0.0)
  endif

  call cpu_clock_end(id_clock_offline_tracer)

end procedure step_offline
module procedure initialize_MOM
  type(ocean_grid_type),  pointer :: G => NULL()    ! A pointer to the metric grid use for the run
  type(ocean_grid_type),  pointer :: G_in => NULL() ! Pointer to the input grid
  type(hor_index_type),   pointer :: HI => NULL()   ! A hor_index_type for array extents
  type(hor_index_type),   target  :: HI_in          ! HI on the input grid
  type(hor_index_type)            :: HI_in_unmasked ! HI on the unmasked input grid
  type(verticalGrid_type), pointer :: GV => NULL()
  type(dyn_horgrid_type), pointer :: dG => NULL(), test_dG => NULL()
  type(dyn_horgrid_type), pointer :: dG_in => NULL()
  type(dyn_horgrid_type), pointer :: dG_unmasked_in => NULL()
  type(diag_ctrl),        pointer :: diag => NULL()
  type(unit_scale_type),  pointer :: US => NULL()
  type(MOM_restart_CS),   pointer :: restart_CSp => NULL()
  character(len=4), parameter :: vers_num = 'v2.0'
  integer :: turns   ! Number of grid quarter-turns
  logical :: point_calving
  real, allocatable         :: u_in(:,:,:) ! Initial zonal velocities [L T-1 ~> m s-1]
  real, allocatable         :: v_in(:,:,:) ! Initial meridional velocities [L T-1 ~> m s-1]
  real, allocatable         :: h_in(:,:,:) ! Initial layer thicknesses [H ~> m or kg m-2]
  real, allocatable, target :: frac_shelf_in(:,:) ! Initial fraction of the total cell area occupied
  real, allocatable, target :: mass_shelf_in(:,:) ! Initial mass of ice shelf contained within a grid cell
  real, allocatable, target :: T_in(:,:,:) ! Initial temperatures [C ~> degC]
  real, allocatable, target :: S_in(:,:,:) ! Initial salinities [S ~> ppt]
  type(ocean_OBC_type), pointer :: OBC_in => NULL()
  type(sponge_CS), pointer :: sponge_in_CSp => NULL()
  type(ALE_sponge_CS), pointer :: ALE_sponge_in_CSp => NULL()
  type(oda_incupd_CS),pointer :: oda_incupd_in_CSp => NULL()
# include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  real    :: dtbt              ! If negative, this specifies the barotropic timestep as a fraction
  real, allocatable, dimension(:,:)   :: eta ! free surface height or column mass [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: h_new    ! Layer thicknesses after regridding [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: dzRegrid ! The change in grid interface positions due to regridding,
  real, allocatable, dimension(:,:,:) :: h_old_u  ! Source grid thickness at zonal velocity points [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: h_old_v  ! Source grid thickness at meridional velocity
  real, allocatable, dimension(:,:,:) :: h_new_u  ! Destination grid thickness at zonal
  real, allocatable, dimension(:,:,:) :: h_new_v  ! Destination grid thickness at meridional
  logical, allocatable, dimension(:,:,:) :: PCM_cell ! If true, PCM remapping should be used in a cell.
  type(group_pass_type) :: tmp_pass_uv_T_S_h, pass_uv_T_S_h
  real    :: Hmix_z, Hmix_UV_z ! Temporary variables with averaging depths [Z ~> m]
  real    :: HFrz_z            ! Temporary variable with the melt potential depth [Z ~> m]
  real    :: default_val       ! The default value for DTBT_RESET_PERIOD [s]
  logical :: write_geom_files  ! If true, write out the grid geometry files.
  logical :: new_sim           ! If true, this has been determined to be a new simulation
  logical :: use_geothermal    ! If true, apply geothermal heating.
  logical :: use_EOS           ! If true, density calculated from T & S using an equation of state.
  logical :: symmetric         ! If true, use symmetric memory allocation.
  logical :: save_IC           ! If true, save the initial conditions.
  logical :: do_unit_tests     ! If true, call unit tests.
  logical :: fpmix             ! Needed to decide if BLD should be passed to RK2.
  logical :: test_grid_copy = .false.
  logical :: bulkmixedlayer    ! If true, a refined bulk mixed layer scheme is used
  logical :: use_temperature   ! If true, temperature and salinity used as state variables.
  logical :: use_p_surf_in_EOS ! If true, always include the surface pressure contributions
  logical :: use_frazil        ! If true, liquid seawater freezes if temp below freezing,
  logical :: bound_salinity    ! If true, salt is added to keep salinity above
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: use_conT_absS     ! If true, the prognostics T & S are conservative temperature
  logical :: advect_TS         ! If false, then no horizontal advection of temperature
  logical :: use_ice_shelf     ! Needed for ALE
  logical :: global_indexing   ! If true use global horizontal index values instead
  logical :: bathy_at_vel      ! If true, also define bathymetric fields at the
  logical :: calc_dtbt         ! Indicates whether the dynamically adjusted barotropic
  logical :: debug_truncations ! If true, turn on diagnostics useful for debugging truncations.
  integer :: first_direction   ! An integer that indicates which direction is to be
  logical :: enable_bugs       ! If true, the defaults for certain recently added bug-fix flags are
  logical :: non_Bous          ! If true, this run is fully non-Boussinesq
  logical :: Boussinesq        ! If true, this run is fully Boussinesq
  logical :: semi_Boussinesq   ! If true, this run is partially non-Boussinesq
  logical :: use_KPP           ! If true, diabatic is using KPP vertical mixing
  logical :: MLE_use_PBL_MLD   ! If true, use stored boundary layer depths for submesoscale restratification.
  logical :: OBC_reservoir_init_bug
  integer :: nkml, nkbl, verbosity, write_geom, number_of_OBC_segments
  integer :: dynamics_stencil  ! The computational stencil for the calculations
  real :: salin_underflow      ! A tiny value of salinity below which the it is set to 0 [S ~> ppt]
  real :: temp_underflow       ! A tiny magnitude of temperatures below which they are set to 0 [C ~> degC]
  real :: conv2watt            ! A conversion factor from temperature fluxes to heat
  real :: conv2salt            ! A conversion factor for salt fluxes [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  character(len=48) :: S_flux_units
  type(vardesc) :: vd_T, vd_S  ! Structures describing temperature and salinity variables.
  type(time_type)                 :: Start_time
  type(ocean_internal_state)      :: MOM_internal_state
  type(MOM_domain_type), pointer  :: MOM_dom_unmasked => null() ! Unmasked MOM domain instance
  character(len=240) :: geom_file ! Name of the ocean geometry file
  CS%Time => Time

  id_clock_ocean    = cpu_clock_id('Ocean', grain=CLOCK_COMPONENT)
  id_clock_init = cpu_clock_id('Ocean Initialization', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_ocean) ; call cpu_clock_begin(id_clock_init)

  Start_time = Time ; if (present(Time_in)) Start_time = Time_in

  ! Read paths and filenames from namelist and store in "dirs".
  ! Also open the parsed input parameter file(s) and setup param_file.
  call get_MOM_input(param_file, dirs, default_input_filename=input_restart_file, ensemble_num=ensemble_num)

  verbosity = 2 ; call read_param(param_file, "VERBOSITY", verbosity)
  call MOM_set_verbosity(verbosity, .true.)
  call callTree_enter("initialize_MOM(), MOM.F90")

  call find_obsolete_params(param_file)

  ! Determining the internal unit scaling factors for this run.
  call unit_scaling_init(param_file, CS%US)
  US => CS%US

  ! Read relevant parameters and write them to the model log.
  call log_version(param_file, "MOM", version, "", log_to_all=.true., layout=.true., debugging=.true.)
  call get_param(param_file, "MOM", "VERBOSITY", verbosity,  &
                 "Integer controlling level of messaging\n" // &
                 "\t0 = Only FATAL messages\n" // &
                 "\t2 = Only FATAL, WARNING, NOTE [default]\n" // &
                 "\t6 = Above plus call tree messages\n" //&
                 "\t9 = All)", default=2, debuggingParam=.true.)
  call get_param(param_file, "MOM", "DO_UNIT_TESTS", do_unit_tests, &
                 "If True, exercises unit tests at model start up.", &
                 default=.false., debuggingParam=.true.)
  if (do_unit_tests) then
    id_clock_unit_tests = cpu_clock_id('(Ocean unit tests)', grain=CLOCK_MODULE)
    call cpu_clock_begin(id_clock_unit_tests)
    call unit_tests(verbosity)
    call cpu_clock_end(id_clock_unit_tests)
  endif

  call get_param(param_file, "MOM", "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)
  call get_param(param_file, "MOM", "SPLIT_RK2B", CS%use_alt_split, &
                 "If true, use a version of the split explicit time stepping scheme that "//&
                 "exchanges velocities with step_MOM that have the average barotropic phase over "//&
                 "a baroclinic timestep rather than the instantaneous barotropic phase.", &
                 default=.false., do_not_log=.not.CS%split)
  if (CS%split) then
    CS%use_RK2 = .false.
  else
    call get_param(param_file, "MOM", "USE_RK2", CS%use_RK2, &
                 "If true, use RK2 instead of RK3 in the unsplit time stepping.", &
                 default=.false.)
  endif

  ! FPMIX is needed to decide if boundary layer depth should be passed to RK2
  call get_param(param_file, '', "FPMIX", fpmix, &
                 "If true, add non-local momentum flux increments and diffuse down the Eulerian gradient.", &
                 default=.false., do_not_log=.true.)

  if (fpmix .and. .not. CS%split)  then
    call MOM_error(FATAL, "initialize_MOM: "//&
       "FPMIX=True only works when SPLIT=True.")
  endif

  call get_param(param_file, "MOM", "BOUSSINESQ", Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true., do_not_log=.true.)
  call get_param(param_file, "MOM", "SEMI_BOUSSINESQ", semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=.true.)
  non_Bous = .not.(Boussinesq .or. semi_Boussinesq)
  call get_param(param_file, "MOM", "CALC_RHO_FOR_SEA_LEVEL", CS%calc_rho_for_sea_lev, &
                 "If true, the in-situ density is used to calculate the "//&
                 "effective sea level that is returned to the coupler. If false, "//&
                 "the Boussinesq parameter RHO_0 is used.", default=non_Bous)
  call get_param(param_file, "MOM", "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, "MOM", "USE_EOS", use_EOS, &
                 "If true,  density is calculated from temperature and "//&
                 "salinity with an equation of state.  If USE_EOS is "//&
                 "true, ENABLE_THERMODYNAMICS must be true as well.", &
                 default=use_temperature)
  call get_param(param_file, "MOM", "DIABATIC_FIRST", CS%diabatic_first, &
                 "If true, apply diabatic and thermodynamic processes, "//&
                 "including buoyancy forcing and mass gain or loss, "//&
                 "before stepping the dynamics forward.", default=.false.)
  call get_param(param_file, "MOM", "USE_CONTEMP_ABSSAL", use_conT_absS, &
                 "If true, the prognostics T&S are the conservative temperature "//&
                 "and absolute salinity. Care should be taken to convert them "//&
                 "to potential temperature and practical salinity before "//&
                 "exchanging them with the coupler and/or reporting T&S diagnostics.", &
                 default=.false.)
  CS%tv%T_is_conT = use_conT_absS ; CS%tv%S_is_absS = use_conT_absS
  call get_param(param_file, "MOM", "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is true.  "//&
                 "This assumes that KD = 0.0 and that there is no buoyancy forcing, "//&
                 "but makes the model faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, "MOM", "DO_DYNAMICS", CS%do_dynamics, &
                 "If False, skips the dynamics calls that update u & v, as well as "//&
                 "the gravity wave adjustment to h. This may be a fragile feature, "//&
                 "but can be useful during development", default=.true.)
  call get_param(param_file, "MOM", "ADVECT_TS", advect_TS, &
                 "If True, advect temperature and salinity horizontally "//&
                 "If False, T/S are registered for advection. "//&
                 "This is intended only to be used in offline tracer mode "//&
                 "and is by default false in that case.", &
                 do_not_log=.true., default=.true.)
  if (present(offline_tracer_mode)) then ! Only read this parameter in enabled modes
    call get_param(param_file, "MOM", "OFFLINE_TRACER_MODE", CS%offline_tracer_mode, &
                 "If true, barotropic and baroclinic dynamics, thermodynamics "//&
                 "are all bypassed with all the fields necessary to integrate "//&
                 "the tracer advection and diffusion equation are read in from "//&
                 "files stored from a previous integration of the prognostic model. "//&
                 "NOTE: This option only used in the ocean_solo_driver.", default=.false.)
    if (CS%offline_tracer_mode) then
      call get_param(param_file, "MOM", "ADVECT_TS", advect_TS, &
                   "If True, advect temperature and salinity horizontally "//&
                   "If False, T/S are registered for advection. "//&
                   "This is intended only to be used in offline tracer mode, "//&
                   "and is by default false in that case", &
                   default=.false. )
    endif
  endif
  call get_param(param_file, "MOM", "USE_REGRIDDING", CS%use_ALE_algorithm, &
                 "If True, use the ALE algorithm (regridding/remapping). "//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(param_file, "MOM", "REMAP_UV_USING_OLD_ALG", CS%remap_uv_using_old_alg, &
                 "If true, uses the old remapping-via-a-delta-z method for "//&
                 "remapping u and v. If false, uses the new method that remaps "//&
                 "between grids described by an old and new thickness.", &
                 default=.false., do_not_log=.not.CS%use_ALE_algorithm)
  call get_param(param_file, "MOM", "REMAP_AUXILIARY_VARS", CS%remap_aux_vars, &
                 "If true, apply ALE remapping to all of the auxiliary 3-dimensional "//&
                 "variables that are needed to reproduce across restarts, similarly to "//&
                 "what is already being done with the primary state variables.  "//&
                 "The default should be changed to true.", default=.false., &
                 do_not_log=.not.CS%use_ALE_algorithm)
  call get_param(param_file, "MOM", "BULKMIXEDLAYER", bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer "//&
                 "with transitional buffer layers.  Layers 1 through "//&
                 "NKML+NKBL have variable densities. There must be at "//&
                 "least NKML+NKBL+1 layers if BULKMIXEDLAYER is true. "//&
                 "BULKMIXEDLAYER can not be used with USE_REGRIDDING. "//&
                 "The default is influenced by ENABLE_THERMODYNAMICS.", &
                 default=use_temperature .and. .not.CS%use_ALE_algorithm)
  call get_param(param_file, "MOM", "USE_POROUS_BARRIER", CS%use_porbar, &
                 "If true, use porous barrier to constrain the widths "//&
                 "and face areas at the edges of the grid cells. ", &
                 default=.false.)
  call get_param(param_file, "MOM", "BATHYMETRY_AT_VEL", bathy_at_vel, &
                 "If true, there are separate values for the basin depths "//&
                 "at velocity points.  Otherwise the effects of topography "//&
                 "are entirely determined from thickness points.", &
                 default=.false.)
  call get_param(param_file, "MOM", "USE_WAVES", CS%UseWaves, default=.false., &
                 do_not_log=.true.)

  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, "MOM", "DEBUG_TRUNCATIONS", debug_truncations, &
                 "If true, calculate all diagnostics that are useful for "//&
                 "debugging truncations.", default=.false., debuggingParam=.true.)
  call get_param(param_file, "MOM", "OBC_NUMBER_OF_SEGMENTS", number_of_OBC_segments, &
                 default=0, do_not_log=.true.)
  call get_param(param_file, "MOM", "DEBUG_OBCS", CS%debug_OBCs, &
                 "If true, write out verbose debugging data about OBCs.", &
                 default=.false., debuggingParam=.true., do_not_log=(number_of_OBC_segments<=0))
  call get_param(param_file, "MOM", "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 "If true, the defaults for certain recently added bug-fix flags are set to "//&
                 "recreate the bugs so that the code can be moved forward without changing "//&
                 "answers for existing configurations.  The defaults for groups of bug-fix "//&
                 "flags are periodically changed to correct the bugs, at which point this "//&
                 "parameter will no longer be used to set their default.  Setting this to false "//&
                 "means that bugs are only used if they are actively selected, but it also "//&
                 "means that answers may change when code is updated due to newly found bugs.", &
                 default=.true.)

  call get_param(param_file, "MOM", "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that "//&
                 "is actually used will be an integer fraction of the "//&
                 "forcing time-step (DT_FORCING in ocean-only mode or the "//&
                 "coupling timestep in coupled mode.)", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
  call get_param(param_file, "MOM", "DT_THERM", CS%dt_therm, &
                 "The thermodynamic time step. Ideally DT_THERM should be an "//&
                 "integer multiple of DT and of DT_TRACER_ADVECT "//&
                 "and less than the forcing or coupling time-step. However, if "//&
                 "THERMO_SPANS_COUPLING is true, DT_THERM can be an integer multiple "//&
                 "of the coupling timestep. By default DT_THERM is set to DT.", &
                 units="s", scale=US%s_to_T, default=US%T_to_s*CS%dt)
  call get_param(param_file, "MOM", "THERMO_SPANS_COUPLING", CS%thermo_spans_coupling, &
                 "If true, the MOM will take thermodynamic "//&
                 "timesteps that can be longer than the coupling timestep. "//&
                 "The actual thermodynamic timestep that is used in this "//&
                 "case is the largest integer multiple of the coupling "//&
                 "timestep that is less than or equal to DT_THERM.", default=.false.)
  call get_param(param_file, "MOM", "DT_TRACER_ADVECT", CS%dt_tr_adv, &
                 "The tracer advection time step. Ideally DT_TRACER_ADVECT should be an "//&
                 "integer multiple of DT, less than DT_THERM, and less than the forcing "//&
                 "or coupling time-step. However, if TRADV_SPANS_COUPLING is true, "//&
                 "DT_TRACER_ADVECT can be longer than the coupling timestep. By "//&
                 "default DT_TRACER_ADVECT is set to DT_THERM.", &
                 units="s", scale=US%s_to_T, default=US%T_to_s*CS%dt_therm)
  call get_param(param_file, "MOM", "TRADV_SPANS_COUPLING", CS%tradv_spans_coupling, &
                 "If true, the MOM will take tracer advection "//&
                 "timesteps that can be longer than the coupling timestep. "//&
                 "The actual tracer advection timestep that is used in this "//&
                 "case is the largest integer multiple of the coupling "//&
                 "timestep that is less than or equal to DT_TRACER_ADVECT.", &
                 default=CS%thermo_spans_coupling)
  if ( CS%diabatic_first .and. (CS%dt_tr_adv /= CS%dt_therm) ) then
    call MOM_error(FATAL,"MOM: If using DIABATIC_FIRST, DT_TRACER_ADVECT must equal DT_THERM.")
  endif
  call get_param(param_file, "MOM", "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, isopycnal surfaces are diffused with a Laplacian "//&
                 "coefficient of KHTH.", default=.false.)
  call get_param(param_file, "MOM", "APPLY_INTERFACE_FILTER", CS%interface_filter, &
                 "If true, model interface heights are subjected to a grid-scale "//&
                 "dependent spatial smoothing, often with biharmonic filter.", default=.false.)
  call get_param(param_file, "MOM", "THICKNESSDIFFUSE_FIRST", CS%thickness_diffuse_first, &
                 "If true, do thickness diffusion or interface height smoothing before dynamics.  "//&
                 "This is only used if THICKNESSDIFFUSE or APPLY_INTERFACE_FILTER is true.", &
                 default=.false., do_not_log=.not.(CS%thickness_diffuse.or.CS%interface_filter))
  CS%interface_filter_dt_bug = .false.
  if ((.not.CS%thickness_diffuse_first .and. CS%interface_filter) .or. &
      (CS%thickness_diffuse_first .and. (CS%thickness_diffuse .or. CS%interface_filter) &
          .and. (CS%dt_tr_adv /= CS%dt_therm))) then
    call get_param(param_file, "MOM", "INTERFACE_FILTER_DT_BUG", CS%interface_filter_dt_bug, &
                   "If true, uses the wrong time interval in calls to interface_filter "//&
                   "and thickness_diffuse.  Has no effect when THICKNESSDIFFUSE_FIRST is "//&
                   "true and DT_TRACER_ADVECT = DT_THERMO or when THICKNESSDIFFUSE_FIRST "//&
                   "is false and APPLY_INTERFACE_FILTER is false. ", default=.false.)
  endif

  if (bulkmixedlayer) then
    CS%Hmix = -1.0 ; CS%Hmix_UV = -1.0
  else
    call get_param(param_file, "MOM", "HMIX_SFC_PROP", Hmix_z, &
                 "If BULKMIXEDLAYER is false, HMIX_SFC_PROP is the depth "//&
                 "over which to average to find surface properties like "//&
                 "SST and SSS or density (but not surface velocities).", &
                 units="m", default=1.0, scale=US%m_to_Z)
    call get_param(param_file, "MOM", "HMIX_UV_SFC_PROP", Hmix_UV_z, &
                 "If BULKMIXEDLAYER is false, HMIX_UV_SFC_PROP is the depth "//&
                 "over which to average to find surface flow properties, "//&
                 "SSU, SSV. A non-positive value indicates no averaging.", &
                 units="m", default=0.0, scale=US%m_to_Z)
  endif
  call get_param(param_file, "MOM", "HFREEZE", HFrz_z, &
                 "If HFREEZE > 0, melt potential will be computed. The actual depth "//&
                 "over which melt potential is computed will be min(HFREEZE, OBLD), "//&
                 "where OBLD is the boundary layer depth. If HFREEZE <= 0 (default), "//&
                 "melt potential will not be computed.", &
                 units="m", default=-1.0, scale=US%m_to_Z)
  call get_param(param_file, "MOM", "INTERPOLATE_P_SURF", CS%interp_p_surf, &
                 "If true, linearly interpolate the surface pressure "//&
                 "over the coupling time step, using the specified value "//&
                 "at the end of the step.", default=.false.)

  if (CS%split) then
    call get_param(param_file, "MOM", "DTBT", dtbt, units="s or nondim", default=-0.98)
    default_val = US%T_to_s*CS%dt_therm ; if (dtbt > 0.0) default_val = -1.0
    CS%dtbt_reset_period = -1.0
    call get_param(param_file, "MOM", "DTBT_RESET_PERIOD", CS%dtbt_reset_period, &
                 "The period between recalculations of DTBT (if DTBT <= 0). "//&
                 "If DTBT_RESET_PERIOD is negative, DTBT is set based "//&
                 "only on information available at initialization.  If 0, "//&
                 "DTBT will be set every dynamics time step. The default "//&
                 "is set by DT_THERM.  This is only used if SPLIT is true.", &
                 units="s", default=default_val, scale=US%s_to_T, do_not_read=(dtbt > 0.0))
  endif

  ! This is here in case these values are used inappropriately.
  use_frazil = .false. ; bound_salinity = .false. ; use_p_surf_in_EOS = .false.
  CS%tv%P_Ref = 2.0e7*US%Pa_to_RL2_T2
  if (use_temperature) then
    call get_param(param_file, "MOM", "FRAZIL", use_frazil, &
                 "If true, water freezes if it gets too cold, and the "//&
                 "accumulated heat deficit is returned in the "//&
                 "surface state.  FRAZIL is only used if "//&
                 "ENABLE_THERMODYNAMICS is true.", default=.false.)
    call get_param(param_file, "MOM", "DO_GEOTHERMAL", use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
    call get_param(param_file, "MOM", "BOUND_SALINITY", bound_salinity, &
                 "If true, limit salinity to being positive. (The sea-ice "//&
                 "model may ask for more salt than is available and "//&
                 "drive the salinity negative otherwise.)", default=.false.)
    call get_param(param_file, "MOM", "MIN_SALINITY", CS%tv%min_salinity, &
                 "The minimum value of salinity when BOUND_SALINITY=True.", &
                 units="PPT", default=0.0, scale=US%ppt_to_S, do_not_log=.not.bound_salinity)
    call get_param(param_file, "MOM", "SALINITY_UNDERFLOW", salin_underflow, &
                 "A tiny value of salinity below which the it is set to 0.  For reference, "//&
                 "one molecule of salt per square meter of ocean is of order 1e-29 ppt.", &
                 units="PPT", default=0.0, scale=US%ppt_to_S)
    call get_param(param_file, "MOM", "TEMPERATURE_UNDERFLOW", temp_underflow, &
                 "A tiny magnitude of temperatures below which they are set to 0.", &
                 units="degC", default=0.0, scale=US%degC_to_C)
    call get_param(param_file, "MOM", "C_P", CS%tv%C_p, &
                 "The heat capacity of sea water, approximated as a constant. "//&
                 "This is only used if ENABLE_THERMODYNAMICS is true. The default "//&
                 "value is from the TEOS-10 definition of conservative temperature.", &
                 units="J kg-1 K-1", default=3991.86795711963, scale=US%J_kg_to_Q*US%C_to_degC)
    call get_param(param_file, "MOM", "USE_PSURF_IN_EOS", use_p_surf_in_EOS, &
                 "If true, always include the surface pressure contributions "//&
                 "in equation of state calculations.", default=.true.)
  endif
  if (use_EOS) call get_param(param_file, "MOM", "P_REF", CS%tv%P_Ref, &
                 "The pressure that is used for calculating the coordinate "//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS are true.", &
                 units="Pa", default=2.0e7, scale=US%Pa_to_RL2_T2)

  if (bulkmixedlayer) then
    call get_param(param_file, "MOM", "NKML", nkml, &
                 "The number of sublayers within the mixed layer if "//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=2)
    call get_param(param_file, "MOM", "NKBL", nkbl, &
                 "The number of layers that are used as variable density buffer "//&
                 "layers if BULKMIXEDLAYER is true.", units="nondim", default=2)
  endif

  call get_param(param_file, "MOM", "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so "//&
                 "that corresponding points on different processors have "//&
                 "the same index. This does not work with static memory.", &
                 default=.false., layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if (global_indexing) call MOM_error(FATAL, "initialize_MOM: "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#endif
  call get_param(param_file, "MOM", "FIRST_DIRECTION", first_direction, &
                 "An integer that indicates which direction goes first "//&
                 "in parts of the code that use directionally split "//&
                 "updates, with even numbers (or 0) used for x- first "//&
                 "and odd numbers used for y-first.", default=0)
  call get_param(param_file, "MOM", "ALTERNATE_FIRST_DIRECTION", CS%alternate_first_direction, &
                 "If true, after every dynamic timestep alternate whether the x- or y- "//&
                 "direction updates occur first in directionally split parts of the calculation. "//&
                 "If this is true, FIRST_DIRECTION applies at the start of a new run or if "//&
                 "the next first direction can not be found in the restart file.", default=.false.)
  call get_param(param_file, "MOM", "CHECK_BAD_SURFACE_VALS", CS%check_bad_sfc_vals, &
                 "If true, check the surface state for ridiculous values.", &
                 default=.false.)
  if (CS%check_bad_sfc_vals) then
    call get_param(param_file, "MOM", "BAD_VAL_SSH_MAX", CS%bad_val_ssh_max, &
                 "The value of SSH above which a bad value message is "//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="m", default=20.0, scale=US%m_to_Z)
    call get_param(param_file, "MOM", "BAD_VAL_SSS_MAX", CS%bad_val_sss_max, &
                 "The value of SSS above which a bad value message is "//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="PPT", default=45.0, scale=US%ppt_to_S)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MAX", CS%bad_val_sst_max, &
                 "The value of SST above which a bad value message is "//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=45.0, scale=US%degC_to_C)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MIN", CS%bad_val_sst_min, &
                 "The value of SST below which a bad value message is "//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=-2.1, scale=US%degC_to_C)
    call get_param(param_file, "MOM", "BAD_VAL_COLUMN_THICKNESS", CS%bad_val_col_thick, &
                 "The value of column thickness below which a bad value message is "//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="m", default=0.0, scale=US%m_to_Z)
  endif
  call get_param(param_file, "MOM", "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, "MOM", "SURFACE_ANSWER_DATE", CS%answer_date, &
               "The vintage of the expressions for the surface properties.  Values below "//&
               "20190101 recover the answers from the end of 2018, while higher values "//&
               "use updated and more robust forms of the same expressions.", &
               default=default_answer_date, do_not_log=non_Bous)
  if (non_Bous) CS%answer_date = 99991231

  call get_param(param_file, "MOM", "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, write the initial conditions to a file given "//&
                 "by IC_OUTPUT_FILE.", default=.false.)
  call get_param(param_file, "MOM", "IC_OUTPUT_FILE", CS%IC_file, &
                 "The file into which to write the initial conditions.", &
                 default="MOM_IC")
  call get_param(param_file, "MOM", "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files. "//&
                 "If =1, write the geometry and vertical grid files only for "//&
                 "a new simulation. If =2, always write the geometry and "//&
                 "vertical grid files. Other values are invalid.", default=1)
  if (write_geom<0 .or. write_geom>2) call MOM_error(FATAL,"MOM: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  call get_param(param_file, "MOM", "GEOM_FILE", geom_file, &
                 "The file into which to write the ocean geometry.", &
                 default="ocean_geometry")
  call get_param(param_file, "MOM", "USE_DBCLIENT", CS%use_dbclient, &
                 "If true, initialize a client to a remote database that can "//&
                 "be used for online analysis and machine-learning inference.",&
                 default=.false.)

  ! Check for inconsistent parameter settings.
  if (CS%use_ALE_algorithm .and. bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: BULKMIXEDLAYER can not currently be used with the ALE algorithm.")
  if (CS%use_ALE_algorithm .and. .not.use_temperature) call MOM_error(FATAL, &
     "MOM: At this time, USE_EOS should be True when using the ALE algorithm.")
  if (CS%adiabatic .and. use_temperature) call MOM_error(WARNING, &
    "MOM: ADIABATIC and ENABLE_THERMODYNAMICS both defined is usually unwise.")
  if (use_EOS .and. .not.use_temperature) call MOM_error(FATAL, &
    "MOM: ENABLE_THERMODYNAMICS must be defined to use USE_EOS.")
  if (CS%adiabatic .and. bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: ADIABATIC and BULKMIXEDLAYER can not both be defined.")
  if (bulkmixedlayer .and. .not.use_EOS) call MOM_error(FATAL, &
      "initialize_MOM: A bulk mixed layer can only be used with T & S as "//&
      "state variables. Add USE_EOS = True to MOM_input.")

  use_ice_shelf = .false.
  if (present(ice_shelf_CSp)) then
    call get_param(param_file, "MOM", "ICE_SHELF", use_ice_shelf, &
       "If true, enables the ice shelf model.", default=.false.)
  endif

  call get_param(param_file, "MOM", "USE_PARTICLES", CS%use_particles, &
                 "If true, use the particles package.", default=.false.)
  call get_param(param_file, "MOM", "USE_UH_PARTICLES", CS%use_uh_particles, &
                 "If true, use the uh velocity in the particles package.", &
                 default=.false., do_not_log=.not.CS%use_particles)
  call get_param(param_file, "MOM", "UH_PARTICLES_BUG", CS%uh_particles_bug, &
                 "If true, use a bug in which the particles are advected inconsistently"//&
                 "with the dynamics timestep instead of the tracer timestep.", &
                 default=enable_bugs, do_not_log=.not.CS%use_uh_particles)
  CS%ensemble_ocean=.false.
  call get_param(param_file, "MOM", "ENSEMBLE_OCEAN", CS%ensemble_ocean, &
                 "If False, The model is being run in serial mode as a single realization. "//&
                 "If True, The current model realization is part of a larger ensemble "//&
                 "and at the end of step MOM, we will perform a gather of the ensemble "//&
                 "members for statistical evaluation and/or data assimilation.", default=.false.)

  call callTree_waypoint("MOM parameters read (initialize_MOM)")

  call get_param(param_file, "MOM", "HOMOGENIZE_FORCINGS", CS%homogenize_forcings, &
                 "If True, homogenize the forces and fluxes.", default=.false.)
  call get_param(param_file, "MOM", "UPDATE_USTAR",CS%update_ustar, &
                 "If True, update ustar from homogenized tau when using the "//&
                 "HOMOGENIZE_FORCINGS option.  Note that this will not work "//&
                 "with a non-zero gustiness factor.", default=.false., &
                 do_not_log=.not.CS%homogenize_forcings)

  ! Grid rotation test
  call get_param(param_file, "MOM", "ROTATE_INDEX", CS%rotate_index, &
      "Enable rotation of the horizontal indices.", default=.false., &
      debuggingParam=.true.)
  if (CS%rotate_index) then
    ! TODO: Index rotation currently only works when index rotation does not
    !   change the MPI rank of each domain.  Resolving this will require a
    !   modification to FMS PE assignment.
    !   For now, we only permit single-core runs.

    if (num_PEs() /= 1) &
      call MOM_error(FATAL, "Index rotation is only supported on one PE.")

    ! Alternate_first_direction is not permitted with index rotation.
    !   This feature can be added later in the future if needed.
    if (CS%alternate_first_direction) &
      call MOM_error(FATAL, "Alternating_first_direction is not compatible with index rotation.")

    call get_param(param_file, "MOM", "INDEX_TURNS", turns, &
        "Number of counterclockwise quarter-turn index rotations.", &
        default=1, debuggingParam=.true.)
  else
    turns = 0
  endif

  ! Set up the model domain and grids.
#ifdef SYMMETRIC_MEMORY_
  symmetric = .true.
#else
  symmetric = .false.
#endif
  G_in => CS%G_in
#ifdef STATIC_MEMORY_
  call MOM_domains_init(G_in%domain, param_file, symmetric=symmetric, &
                        static_memory=.true., NIHALO=NIHALO_, NJHALO=NJHALO_, &
                        NIGLOBAL=NIGLOBAL_, NJGLOBAL=NJGLOBAL_, NIPROC=NIPROC_, &
                        NJPROC=NJPROC_, US=US, MOM_dom_unmasked=MOM_dom_unmasked)
#else
  call MOM_domains_init(G_in%domain, param_file, symmetric=symmetric, &
                        domain_name="MOM_in", US=US, MOM_dom_unmasked=MOM_dom_unmasked)
#endif

  ! Copy input grid (G_in) domain to active grid G
  ! Swap axes for quarter and 3-quarter turns
  if (CS%rotate_index) then
    allocate(CS%G)
    call clone_MOM_domain(G_in%Domain, CS%G%Domain, turns=turns, domain_name="MOM_rot")
  else
    CS%G => G_in
  endif

  ! TODO: It is unlikely that test_grid_copy and rotate_index would work at the
  !   same time.  It may be possible to enable both but for now we prevent it.
  if (test_grid_copy .and. CS%rotate_index) &
    call MOM_error(FATAL, "Grid cannot be copied during index rotation.")

  if (test_grid_copy) then ; allocate(G)
  else ; G => CS%G ; endif

  call callTree_waypoint("domains initialized (initialize_MOM)")

  call MOM_debugging_init(param_file)
  call diag_mediator_infrastructure_init()
  call MOM_io_init(param_file)

  ! Create HI and dG on the input index map.
  call hor_index_init(G_in%Domain, HI_in, param_file, &
                      local_indexing=.not.global_indexing)
  call create_dyn_horgrid(dG_in, HI_in, bathymetry_at_vel=bathy_at_vel)
  call clone_MOM_domain(G_in%Domain, dG_in%Domain)
  ! Also allocate the input ocean_grid_type type at this point based on the same information.
  call MOM_grid_init(G_in, param_file, US, HI_in, bathymetry_at_vel=bathy_at_vel)

  ! Allocate initialize time-invariant MOM variables.
  call MOM_initialize_fixed(dG_in, US, OBC_in, param_file)

  call get_param(param_file, "MOM", "DT_OBC_SEG_UPDATE_OBGC", CS%dt_obc_seg_period, &
                 "The time between OBC segment data updates for OBGC tracers.  This must be an "//&
                 "integer multiple of DT and DT_THERM.  The default is set to DT.", units="s", &
                 default=US%T_to_s*CS%dt, scale=US%s_to_T, do_not_log=.not.associated(OBC_in))

  ! Copy the grid metrics and bathymetry to the ocean_grid_type
  call copy_dyngrid_to_MOM_grid(dG_in, G_in, US)

  call callTree_waypoint("returned from MOM_initialize_fixed() (initialize_MOM)")

  call verticalGridInit( param_file, CS%GV, US )
  GV => CS%GV

  ! Now that the vertical grid has been initialized, rescale parameters that depend on factors
  ! that are set with the vertical grid to their desired units.  This added rescaling step would
  ! be unnecessary if the vertical grid were initialized earlier in this routine.
  if (.not.bulkmixedlayer) then
    CS%Hmix = (US%Z_to_m * GV%m_to_H) * Hmix_z
    CS%Hmix_UV = (US%Z_to_m * GV%m_to_H) * Hmix_UV_z
  endif
  CS%HFrz = (US%Z_to_m * GV%m_to_H) * HFrz_z

  !   Shift from using the temporary dynamic grid type to using the final (potentially static)
  ! and properly rotated ocean-specific grid type and horizontal index type.
  if (CS%rotate_index) then
    allocate(HI)
    call rotate_hor_index(HI_in, turns, HI)
    ! NOTE: If indices are rotated, then G and G_in must both be initialized separately, and
    ! the dynamic grid must be created to handle the grid rotation. G%domain has already been
    ! initialized above.
    call MOM_grid_init(G, param_file, US, HI, bathymetry_at_vel=bathy_at_vel)
    call create_dyn_horgrid(dG, HI, bathymetry_at_vel=bathy_at_vel)
    call clone_MOM_domain(G%Domain, dG%Domain)
    call rotate_dyn_horgrid(dG_in, dG, US, turns)
    call copy_dyngrid_to_MOM_grid(dG, G, US)

    if (associated(OBC_in)) then
      allocate(CS%OBC)
      call rotate_OBC_config(OBC_in, dG_in, CS%OBC, dG, turns)
    endif

    call destroy_dyn_horgrid(dG)
  else
    ! If not rotated, then G_in and G are the same grid.
    HI => HI_in
    G => G_in
    CS%OBC => OBC_in
  endif
  ! dG_in is retained for now so that it can be used with write_ocean_geometry_file() below.

  if (is_root_PE()) call check_MOM6_scaling_factors(CS%GV, US)

  call callTree_waypoint("grids initialized (initialize_MOM)")

  call MOM_timing_init(CS)

  call tracer_registry_init(param_file, CS%tracer_Reg)

  ! Allocate and initialize space for the primary time-varying MOM variables.
  is   = HI%isc   ; ie   = HI%iec  ; js   = HI%jsc  ; je   = HI%jec ; nz = GV%ke
  isd  = HI%isd   ; ied  = HI%ied  ; jsd  = HI%jsd  ; jed  = HI%jed
  IsdB = HI%IsdB  ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB
  ALLOC_(CS%u(IsdB:IedB,jsd:jed,nz))   ; CS%u(:,:,:) = 0.0
  ALLOC_(CS%v(isd:ied,JsdB:JedB,nz))   ; CS%v(:,:,:) = 0.0
  ALLOC_(CS%h(isd:ied,jsd:jed,nz))     ; CS%h(:,:,:) = GV%Angstrom_H
  ALLOC_(CS%uh(IsdB:IedB,jsd:jed,nz))  ; CS%uh(:,:,:) = 0.0
  ALLOC_(CS%vh(isd:ied,JsdB:JedB,nz))  ; CS%vh(:,:,:) = 0.0
  if (use_temperature) then
    ALLOC_(CS%T(isd:ied,jsd:jed,nz))   ; CS%T(:,:,:) = 0.0
    ALLOC_(CS%S(isd:ied,jsd:jed,nz))   ; CS%S(:,:,:) = 0.0
    CS%tv%T => CS%T ; CS%tv%S => CS%S
    if (CS%tv%T_is_conT) then
      vd_T = var_desc(name="contemp", units="Celsius", longname="Conservative Temperature", &
                      cmor_field_name="bigthetao", cmor_longname="Sea Water Conservative Temperature", &
                      conversion=US%C_to_degC)
    else
      vd_T = var_desc(name="temp", units="degC", longname="Potential Temperature", &
                      cmor_field_name="thetao", cmor_longname="Sea Water Potential Temperature", &
                      conversion=US%C_to_degC)
    endif
    if (CS%tv%S_is_absS) then
      vd_S = var_desc(name="abssalt", units="g kg-1", longname="Absolute Salinity", &
                      cmor_field_name="absso", cmor_longname="Sea Water Absolute Salinity", &
                      conversion=US%S_to_ppt)
    else
      vd_S = var_desc(name="salt", units="psu", longname="Salinity", &
                      cmor_field_name="so", cmor_longname="Sea Water Salinity", &
                      conversion=US%S_to_ppt)
    endif

    if (advect_TS) then
      S_flux_units = get_tr_flux_units(GV, "psu") ! Could change to "kg m-2 s-1"?
      conv2watt    = GV%H_to_kg_m2 * US%Q_to_J_kg*CS%tv%C_p
      if (GV%Boussinesq) then
        conv2salt = US%S_to_ppt*GV%H_to_m ! Could change to US%S_to_ppt*GV%H_to_kg_m2 * 0.001?
      else
        conv2salt = US%S_to_ppt*GV%H_to_kg_m2
      endif
      call register_tracer(CS%tv%T, CS%tracer_Reg, param_file, HI, GV, &
                           tr_desc=vd_T, registry_diags=.true., conc_scale=US%C_to_degC, &
                           flux_nameroot='T', flux_units='W', flux_longname='Heat', &
                           net_surfflux_name='KPP_QminusSW', NLT_budget_name='KPP_NLT_temp_budget', &
                           net_surfflux_longname='Net temperature flux ignoring short-wave, as used by [CVMix] KPP', &
                           flux_scale=conv2watt, convergence_units='W m-2', &
                           convergence_scale=conv2watt, CMOR_tendprefix="opottemp", &
                           diag_form=2, underflow_conc=temp_underflow, Tr_out=CS%tv%tr_T)
      call register_tracer(CS%tv%S, CS%tracer_Reg, param_file, HI, GV, &
                           tr_desc=vd_S, registry_diags=.true., conc_scale=US%S_to_ppt, &
                           flux_nameroot='S', flux_units=S_flux_units, flux_longname='Salt', &
                           net_surfflux_name='KPP_netSalt', NLT_budget_name='KPP_NLT_saln_budget', &
                           flux_scale=conv2salt, convergence_units='kg m-2 s-1', &
                           convergence_scale=0.001*US%S_to_ppt*GV%H_to_kg_m2, CMOR_tendprefix="osalt", &
                           diag_form=2, underflow_conc=salin_underflow, Tr_out=CS%tv%tr_S)
    endif
  endif

  if (use_p_surf_in_EOS) allocate(CS%tv%p_surf(isd:ied,jsd:jed), source=0.0)
  if (use_frazil) then
    allocate(CS%tv%frazil(isd:ied,jsd:jed), source=0.0)
    CS%tv%frazil_was_reset = .true.
  endif
  if (bound_salinity) allocate(CS%tv%salt_deficit(isd:ied,jsd:jed), source=0.0)

  allocate(CS%Hml(isd:ied,jsd:jed), source=0.0)

  if (bulkmixedlayer) then
    GV%nkml = nkml ; GV%nk_rho_varies = nkml + nkbl
  else
    GV%nkml = 0 ; GV%nk_rho_varies = 0
  endif
  if (CS%use_ALE_algorithm) then
    call get_param(param_file, "MOM", "NK_RHO_VARIES", GV%nk_rho_varies, default=0) ! Will default to nz later... -AJA
  endif

  ALLOC_(CS%uhtr(IsdB:IedB,jsd:jed,nz)) ; CS%uhtr(:,:,:) = 0.0
  ALLOC_(CS%vhtr(isd:ied,JsdB:JedB,nz)) ; CS%vhtr(:,:,:) = 0.0
  CS%t_dyn_rel_adv = 0.0 ; CS%t_dyn_rel_thermo = 0.0 ; CS%t_dyn_rel_diag = 0.0
  CS%n_dyn_steps_in_adv = 0

  if (debug_truncations) then
    allocate(CS%u_prev(IsdB:IedB,jsd:jed,nz), source=0.0)
    allocate(CS%v_prev(isd:ied,JsdB:JedB,nz), source=0.0)
    MOM_internal_state%u_prev => CS%u_prev
    MOM_internal_state%v_prev => CS%v_prev
    call safe_alloc_ptr(CS%ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(CS%ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    if (.not.CS%adiabatic) then
      call safe_alloc_ptr(CS%ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
      call safe_alloc_ptr(CS%ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
    endif
  endif

  MOM_internal_state%u => CS%u ; MOM_internal_state%v => CS%v
  MOM_internal_state%h => CS%h
  MOM_internal_state%uh => CS%uh ; MOM_internal_state%vh => CS%vh
  if (use_temperature) then
    MOM_internal_state%T => CS%T ; MOM_internal_state%S => CS%S
  endif

  CS%CDp%uh => CS%uh ; CS%CDp%vh => CS%vh

  if (CS%interp_p_surf) allocate(CS%p_surf_prev(isd:ied,jsd:jed), source=0.0)

  ALLOC_(CS%ssh_rint(isd:ied,jsd:jed)) ; CS%ssh_rint(:,:) = 0.0
  ALLOC_(CS%ave_ssh_ibc(isd:ied,jsd:jed)) ; CS%ave_ssh_ibc(:,:) = 0.0
  ALLOC_(CS%eta_av_bc(isd:ied,jsd:jed)) ; CS%eta_av_bc(:,:) = 0.0 ! -G%Z_ref
  CS%time_in_cycle = 0.0 ; CS%time_in_thermo_cycle = 0.0

  !allocate porous topography variables
  allocate(CS%pbv%por_face_areaU(IsdB:IedB,jsd:jed,nz), source=1.0)
  allocate(CS%pbv%por_face_areaV(isd:ied,JsdB:JedB,nz), source=1.0)
  allocate(CS%pbv%por_layer_widthU(IsdB:IedB,jsd:jed,nz+1), source=1.0)
  allocate(CS%pbv%por_layer_widthV(isd:ied,JsdB:JedB,nz+1), source=1.0)

  ! Use the Wright equation of state by default, unless otherwise specified
  ! Note: this line and the following block ought to be in a separate
  ! initialization routine for tv.
  if (use_EOS) then
    allocate(CS%tv%eqn_of_state)
    call EOS_init(param_file, CS%tv%eqn_of_state, US, use_conT_absS)
  endif
  if (use_temperature) then
    allocate(CS%tv%TempxPmE(isd:ied,jsd:jed), source=0.0)
    if (use_geothermal) then
      allocate(CS%tv%internal_heat(isd:ied,jsd:jed), source=0.0)
    endif
  endif
  call callTree_waypoint("state variables allocated (initialize_MOM)")

  ! Set the fields that are needed for bitwise identical restarting
  ! the time stepping scheme.
  call restart_init(param_file, CS%restart_CS)
  restart_CSp => CS%restart_CS

  call set_restart_fields(GV, US, param_file, CS, restart_CSp)
  if (CS%split .and. CS%use_alt_split) then
    call register_restarts_dyn_split_RK2b(HI, GV, US, param_file, &
             CS%dyn_split_RK2b_CSp, restart_CSp, CS%uh, CS%vh)
  elseif (CS%split) then
    call register_restarts_dyn_split_RK2(HI, GV, US, param_file, &
             CS%dyn_split_RK2_CSp, restart_CSp, CS%uh, CS%vh)
  elseif (CS%use_RK2) then
    call register_restarts_dyn_unsplit_RK2(HI, GV, param_file, &
           CS%dyn_unsplit_RK2_CSp)
  else
    call register_restarts_dyn_unsplit(HI, GV, param_file, &
           CS%dyn_unsplit_CSp)
  endif

  ! This subroutine calls user-specified tracer registration routines.
  ! Additional calls can be added to MOM_tracer_flow_control.F90.
  call call_tracer_register(G, GV, US, param_file, CS%tracer_flow_CSp, &
                            CS%tracer_Reg, restart_CSp)

  call MEKE_alloc_register_restart(HI, US, param_file, CS%MEKE, restart_CSp)
  call set_visc_register_restarts(HI, G, GV, US, param_file, CS%visc, restart_CSp, use_ice_shelf)
  call mixedlayer_restrat_register_restarts(HI, GV, US, param_file, &
           CS%mixedlayer_restrat_CSp, restart_CSp)

  if (associated(CS%OBC)) then
    ! This call initializes the relevant vertical remapping structures.
    call open_boundary_setup_vert(GV, US, CS%OBC)

    ! Set up remaining information about open boundary conditions that is needed for OBCs.
    ! Package specific changes to OBCs occur here.
    call call_OBC_register(G, GV, US, param_file, CS%update_OBC_CSp, CS%OBC, CS%tracer_Reg)

    ! This is the equivalent to 2 calls to register_segment_tracer (per segment), which
    ! could occur with the call to update_OBC_data or after the main initialization.
    if (use_temperature) &
      call register_temp_salt_segments(GV, US, CS%OBC, CS%tracer_Reg, param_file)
    ! This is the equivalent call to register_temp_salt_segments for external tracers with OBC
    call call_tracer_register_obc_segments(GV, param_file, CS%tracer_flow_CSp, CS%tracer_Reg, CS%OBC)

    ! Set up the thickness reservoirs if using them.
    if (CS%OBC%use_h_res) &
      call segment_thickness_reservoir_init(GV, US, CS%OBC, param_file)

    ! This needs the number of tracers and to have called any code that sets whether
    ! reservoirs are used.
    call open_boundary_register_restarts(HI, GV, US, CS%OBC, CS%tracer_Reg, &
                          param_file, restart_CSp, use_temperature)

    ! This call allocates the arrays on the segments for open boundary data, but it must occur
    ! after any calls to call_tracer_register_obc_segments.
    call initialize_segment_data(GV, US, CS%OBC, param_file, turns, use_temperature)

    if (CS%debug_OBCs) call write_OBC_info(CS%OBC, G, GV, US)
  endif

  if (present(waves_CSp)) then
    call waves_register_restarts(waves_CSp, HI, GV, US, param_file, restart_CSp)
  endif

  if (use_temperature) then
    call stoch_EOS_register_restarts(HI, param_file, CS%stoch_eos_CS, restart_CSp)
  endif

  if (.not. CS%adiabatic) then
    call register_diabatic_restarts(G, GV, US, param_file, CS%int_tide_CSp, restart_CSp, CS%diabatic_CSp)
  endif

  call callTree_waypoint("restart registration complete (initialize_MOM)")
  call restart_registry_lock(restart_CSp)

  ! Write out all of the grid data used by this run.
  new_sim = determine_is_new_run(dirs%input_filename, dirs%restart_input_dir, G_in, restart_CSp)
  write_geom_files = ((write_geom==2) .or. ((write_geom==1) .and. new_sim))
  if (write_geom_files) then
    if (associated(MOM_dom_unmasked)) then
      call hor_index_init(MOM_dom_unmasked, HI_in_unmasked, param_file, &
                          local_indexing=.not.global_indexing)
      call create_dyn_horgrid(dG_unmasked_in, HI_in_unmasked, bathymetry_at_vel=bathy_at_vel)
      call clone_MOM_domain(MOM_dom_unmasked, dG_unmasked_in%Domain)
      call MOM_initialize_fixed(dG_unmasked_in, US, OBC_in, param_file)
      call write_ocean_geometry_file(dG_unmasked_in, param_file, dirs%output_directory, US=US, geom_file=geom_file)
      call deallocate_MOM_domain(MOM_dom_unmasked)
      call destroy_dyn_horgrid(dG_unmasked_in)
    else
      call write_ocean_geometry_file(dG_in, param_file, dirs%output_directory, US=US, geom_file=geom_file)
    endif
  endif
  call destroy_dyn_horgrid(dG_in)

  ! Initialize dynamically evolving fields, perhaps from restart files.
  call cpu_clock_begin(id_clock_MOM_init)
  call MOM_initialize_coord(GV, US, param_file, CS%tv, G%max_depth)
  call callTree_waypoint("returned from MOM_initialize_coord() (initialize_MOM)")

  if (CS%use_ALE_algorithm) then
    call ALE_init(param_file, G, GV, US, G%max_depth, CS%ALE_CSp)
    call callTree_waypoint("returned from ALE_init() (initialize_MOM)")
  endif

  ! Set a few remaining fields that are specific to the ocean grid type.
  if (CS%rotate_index) then
    call set_first_direction(G, modulo(first_direction + turns, 2))
  else
    call set_first_direction(G, modulo(first_direction, 2))
  endif
  ! Allocate the auxiliary non-symmetric domain for debugging or I/O purposes.
  if (CS%debug .or. G%symmetric) then
    call clone_MOM_domain(G%Domain, G%Domain_aux, symmetric=.false.)
  else ; G%Domain_aux => G%Domain ; endif
  ! Copy common variables from the vertical grid to the horizontal grid.
  ! Consider removing this later?
  G%ke = GV%ke

  if (use_ice_shelf) then
    point_calving = .false. ; if (present(calve_ice_shelf_bergs)) point_calving = calve_ice_shelf_bergs
  endif

  if (CS%rotate_index) then
    G_in%ke = GV%ke

    ! Allocate the auxiliary non-symmetric domain for debugging or I/O purposes.
    if (CS%debug .or. G_in%symmetric) then
      call clone_MOM_domain(G_in%Domain, G_in%Domain_aux, symmetric=.false.)
    else ; G_in%Domain_aux => G_in%Domain ; endif

    allocate(u_in(G_in%IsdB:G_in%IedB, G_in%jsd:G_in%jed, nz), source=0.0)
    allocate(v_in(G_in%isd:G_in%ied, G_in%JsdB:G_in%JedB, nz), source=0.0)
    allocate(h_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz), source=GV%Angstrom_H)

    if (use_temperature) then
      allocate(T_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz), source=0.0)
      allocate(S_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz), source=0.0)

      CS%tv%T => T_in
      CS%tv%S => S_in

      if (associated(CS%OBC)) then
        ! Log this parameter in MOM_initialize_state
        call get_param(param_file, "MOM", "OBC_RESERVOIR_INIT_BUG", OBC_reservoir_init_bug, &
                   "If true, set the OBC tracer reservoirs at the startup of a new run from the "//&
                   "interior tracer concentrations regardless of properties that may be explicitly "//&
                   "specified for the reservoir concentrations.", default=enable_bugs, do_not_log=.true.)
        if (OBC_reservoir_init_bug .and. (allocated(CS%OBC%tres_x) .or. allocated(CS%OBC%tres_y))) &
          call MOM_error(FATAL, "OBC_RESERVOIR_INIT_BUG can not be set to true with grid rotation.")
      endif
    endif

    if (use_ice_shelf) then
      ! These arrays are not initialized in most solo cases, but are needed
      ! when using an ice shelf. Passing the ice shelf diagnostics CS from MOM
      ! for legacy reasons. The actual ice shelf diag CS is internal to the ice shelf
      call initialize_ice_shelf(param_file, G, Time, ice_shelf_CSp, diag_ptr, &
                                Time_init, dirs%output_directory, calve_ice_shelf_bergs=point_calving)
      allocate(frac_shelf_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed), source=0.0)
      allocate(mass_shelf_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed), source=0.0)
      allocate(CS%frac_shelf_h(isd:ied, jsd:jed), source=0.0)
      allocate(CS%mass_shelf(isd:ied, jsd:jed), source=0.0)
      call ice_shelf_query(ice_shelf_CSp, G, CS%frac_shelf_h, CS%mass_shelf)
      ! MOM_initialize_state is using the  unrotated metric
      call rotate_array(CS%frac_shelf_h, -turns, frac_shelf_in)
      call rotate_array(CS%mass_shelf, -turns, mass_shelf_in)
      call MOM_initialize_state(u_in, v_in, h_in, CS%tv, Time, G_in, GV, US, &
          param_file, dirs, restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
          sponge_in_CSp, ALE_sponge_in_CSp, oda_incupd_in_CSp, OBC_in, Time_in, &
          frac_shelf_h=frac_shelf_in, mass_shelf=mass_shelf_in)
    else
      call MOM_initialize_state(u_in, v_in, h_in, CS%tv, Time, G_in, GV, US, &
          param_file, dirs, restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
          sponge_in_CSp, ALE_sponge_in_CSp, oda_incupd_in_CSp, OBC_in, Time_in)
    endif

    if (use_temperature) then
      CS%tv%T => CS%T
      CS%tv%S => CS%S
    endif

    ! Reset the first direction if it was found in a restart file
    if (CS%first_dir_restart > -1.0) then
      call set_first_direction(G, modulo(NINT(CS%first_dir_restart) + turns, 2))
    else
      CS%first_dir_restart = real(modulo(first_direction, 2))
    endif

    call rotate_initial_state(u_in, v_in, h_in, T_in, S_in, use_temperature, &
        turns, CS%u, CS%v, CS%h, CS%T, CS%S)

    if (associated(sponge_in_CSp)) then
      ! TODO: Implementation and testing of non-ALE sponge rotation
      call MOM_error(FATAL, "Index rotation of non-ALE sponge is not yet implemented.")
    endif

    if (associated(ALE_sponge_in_CSp)) then
      call rotate_ALE_sponge(ALE_sponge_in_CSp, G_in, CS%ALE_sponge_CSp, G, GV, US, turns, param_file)
      call update_ALE_sponge_field(CS%ALE_sponge_CSp, T_in, G, GV, CS%T)
      call update_ALE_sponge_field(CS%ALE_sponge_CSp, S_in, G, GV, CS%S)
    endif

   ! Deallocate the unrotated arrays and types that are no longer needed.
    deallocate(u_in)
    deallocate(v_in)
    deallocate(h_in)
    if (use_temperature) then
      deallocate(T_in)
      deallocate(S_in)
    endif
    if (use_ice_shelf) deallocate(frac_shelf_in, mass_shelf_in)
    if (associated(OBC_in)) call open_boundary_end(OBC_in)

  else  ! The model is being run without grid rotation.  This is true of all production runs.
    if (use_ice_shelf) then
      call initialize_ice_shelf(param_file, G, Time, ice_shelf_CSp, diag_ptr, Time_init, &
                               dirs%output_directory, calve_ice_shelf_bergs=point_calving)
      allocate(CS%frac_shelf_h(isd:ied, jsd:jed), source=0.0)
      allocate(CS%mass_shelf(isd:ied, jsd:jed), source=0.0)
      call ice_shelf_query(ice_shelf_CSp,G,CS%frac_shelf_h, CS%mass_shelf)
      call MOM_initialize_state(CS%u, CS%v, CS%h, CS%tv, Time, G, GV, US, &
          param_file, dirs, restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
          CS%sponge_CSp, CS%ALE_sponge_CSp, CS%oda_incupd_CSp, CS%OBC, Time_in, &
          frac_shelf_h=CS%frac_shelf_h, mass_shelf=CS%mass_shelf, OBC_for_bug=CS%OBC)
    else
      call MOM_initialize_state(CS%u, CS%v, CS%h, CS%tv, Time, G, GV, US, &
          param_file, dirs, restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
          CS%sponge_CSp, CS%ALE_sponge_CSp, CS%oda_incupd_CSp, CS%OBC, Time_in, OBC_for_bug=CS%OBC)
    endif

    ! Reset the first direction if it was found in a restart file.
    if (CS%first_dir_restart > -1.0) then
      call set_first_direction(G, NINT(CS%first_dir_restart))
    else
      CS%first_dir_restart = real(modulo(first_direction, 2))
    endif
  endif

  ! Allocate any derived densities or other equation of state derived fields.
  if (.not.(GV%Boussinesq .or. GV%semi_Boussinesq)) then
    allocate(CS%tv%SpV_avg(isd:ied,jsd:jed,nz), source=0.0)
    CS%tv%valid_SpV_halo = -1  ! This array does not yet have any valid data.
  endif

  if (associated(CS%OBC)) then
    call MOM_initialize_OBCs(CS%h, CS%tv, CS%OBC, Time, G, GV, US, param_file, restart_CSp, CS%tracer_Reg)

    if (use_temperature) then
      call pass_var(CS%tv%T, G%Domain, complete=.false.)
      call pass_var(CS%tv%S, G%Domain, complete=.true.)
    endif
    call calc_derived_thermo(CS%tv, CS%h, G, GV, US)

    ! Call this during initialization to fill boundary arrays from fixed values
    call read_OBC_segment_data(G, GV, US, CS%OBC, CS%tv, CS%h, Time)
    call update_OBC_segment_data(G, GV, US, CS%OBC, CS%h, Time)
    call initialize_OBC_segment_reservoirs(GV, CS%OBC)
  endif

  if (use_ice_shelf .and. CS%debug) then
    call hchksum(CS%frac_shelf_h, "MOM:frac_shelf_h", G%HI, haloshift=0)
    call hchksum(CS%mass_shelf, "MOM:mass_shelf", G%HI, haloshift=0, unscale=US%RZ_to_kg_m2)
  endif

  call cpu_clock_end(id_clock_MOM_init)
  call callTree_waypoint("returned from MOM_initialize_state() (initialize_MOM)")

  ! From this point, there may be pointers being set, so the final grid type
  ! that will persist throughout the run has to be used.

  if (test_grid_copy) then
    !  Copy the data from the temporary grid to the dyn_hor_grid to CS%G.
    call create_dyn_horgrid(test_dG, G%HI)
    call clone_MOM_domain(G%Domain, test_dG%Domain)

    call clone_MOM_domain(G%Domain, CS%G%Domain)
    call MOM_grid_init(CS%G, param_file, US)

    call copy_MOM_grid_to_dyngrid(G, test_dG, US)
    call copy_dyngrid_to_MOM_grid(test_dG, CS%G, US)

    call destroy_dyn_horgrid(test_dG)
    call MOM_grid_end(G) ; deallocate(G)

    G => CS%G
    if (CS%debug .or. CS%G%symmetric) then
      call clone_MOM_domain(CS%G%Domain, CS%G%Domain_aux, symmetric=.false.)
    else ; CS%G%Domain_aux => CS%G%Domain ; endif
    G%ke = GV%ke
  endif

  ! At this point, all user-modified initialization code has been called.  The
  ! remainder of this subroutine is controlled by the parameters that have
  ! have already been set.

  if (ALE_remap_init_conds(CS%ALE_CSp) .and. .not. query_initialized(CS%h,"h",restart_CSp)) then
    ! This block is controlled by the ALE parameter REMAP_AFTER_INITIALIZATION.
    ! \todo This block exists for legacy reasons and we should phase it out of all examples. !###
    if (CS%debug) then
      call uvchksum("Pre ALE adjust init cond [uv]", CS%u, CS%v, G%HI, haloshift=1, unscale=US%L_T_to_m_s)
      call hchksum(CS%h,"Pre ALE adjust init cond h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    endif
    call callTree_waypoint("Calling adjustGridForIntegrity() to remap initial conditions (initialize_MOM)")
    call adjustGridForIntegrity(CS%ALE_CSp, G, GV, CS%h )
    if (allocated(CS%tv%SpV_avg)) call calc_derived_thermo(CS%tv, CS%h, G, GV, US, halo=1)
    call pre_ALE_adjustments(G, GV, US, CS%h, CS%tv, CS%tracer_Reg, CS%ALE_CSp, CS%u, CS%v)

    call callTree_waypoint("Calling ALE_regrid() to remap initial conditions (initialize_MOM)")
    allocate(h_new(isd:ied, jsd:jed, nz), source=0.0)
    allocate(dzRegrid(isd:ied, jsd:jed, nz+1), source=0.0)
    allocate(PCM_cell(isd:ied, jsd:jed, nz), source=.false.)
    allocate(h_old_u(IsdB:IedB, jsd:jed, nz), source=0.0)
    allocate(h_new_u(IsdB:IedB, jsd:jed, nz), source=0.0)
    allocate(h_old_v(isd:ied, JsdB:JedB, nz), source=0.0)
    allocate(h_new_v(isd:ied, JsdB:JedB, nz), source=0.0)
    if (use_ice_shelf) then
      call ALE_regrid(G, GV, US, CS%h, h_new, dzRegrid, CS%tv, CS%ALE_CSp, CS%frac_shelf_h, PCM_cell)
    else
      call ALE_regrid(G, GV, US, CS%h, h_new, dzRegrid, CS%tv, CS%ALE_CSp, PCM_cell=PCM_cell)
    endif

    if (callTree_showQuery()) call callTree_waypoint("new grid generated")
    ! Remap all variables from the old grid h onto the new grid h_new
    call ALE_remap_tracers(CS%ALE_CSp, G, GV, CS%h, h_new, CS%tracer_Reg, CS%debug, PCM_cell=PCM_cell)

    ! Determine the old and new grid thicknesses at velocity points.
    call ALE_remap_set_h_vel(CS%ALE_CSp, G, GV, CS%h, h_old_u, h_old_v, CS%OBC, debug=CS%debug)
    if (CS%remap_uv_using_old_alg) then
      call ALE_remap_set_h_vel_via_dz(CS%ALE_CSp, G, GV, h_new, h_new_u, h_new_v, CS%OBC, CS%h, dzRegrid, CS%debug)
    else
      call ALE_remap_set_h_vel(CS%ALE_CSp, G, GV, h_new, h_new_u, h_new_v, CS%OBC, debug=CS%debug)
    endif

    ! Remap the velocity components.
    call ALE_remap_velocities(CS%ALE_CSp, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, CS%u, CS%v, CS%debug)

    if (allocated(CS%tv%SpV_avg)) CS%tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

    ! Replace the old grid with new one.  All remapping must be done at this point.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      CS%h(i,j,k) = h_new(i,j,k)
    enddo ; enddo ; enddo

    deallocate(h_new, dzRegrid, PCM_cell, h_old_u, h_new_u, h_old_v, h_new_v)

    call cpu_clock_begin(id_clock_pass_init)
    call create_group_pass(tmp_pass_uv_T_S_h, CS%u, CS%v, G%Domain)
    if (use_temperature) then
      call create_group_pass(tmp_pass_uv_T_S_h, CS%tv%T, G%Domain)
      call create_group_pass(tmp_pass_uv_T_S_h, CS%tv%S, G%Domain)
    endif
    call create_group_pass(tmp_pass_uv_T_S_h, CS%h, G%Domain)
    call do_group_pass(tmp_pass_uv_T_S_h, G%Domain)
    call cpu_clock_end(id_clock_pass_init)

    if (CS%debug) then
      call uvchksum("Post ALE adjust init cond [uv]", CS%u, CS%v, G%HI, haloshift=1, unscale=US%L_T_to_m_s)
      call hchksum(CS%h, "Post ALE adjust init cond h", G%HI, haloshift=2, unscale=GV%H_to_MKS)
      if (use_temperature) then
        call hchksum(CS%tv%T, "Post ALE adjust init cond T", G%HI, haloshift=2, unscale=US%C_to_degC)
        call hchksum(CS%tv%S, "Post ALE adjust init cond S", G%HI, haloshift=2, unscale=US%S_to_ppt)
      endif
    endif
  endif
  if ( CS%use_ALE_algorithm ) then
   call ALE_set_extrap_boundaries (param_file, CS%ALE_CSp)
   call callTree_waypoint("returned from ALE_init() (initialize_MOM)")
   call ALE_updateVerticalGridType( CS%ALE_CSp, GV )
  endif
  ! The basic state variables have now been fully initialized, so update their halos and
  ! calculate any derived thermodynmics quantities.

  !--- set up group pass for u,v,T,S and h. pass_uv_T_S_h also is used in step_MOM
  call cpu_clock_begin(id_clock_pass_init)
  dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
  call create_group_pass(pass_uv_T_S_h, CS%u, CS%v, G%Domain, halo=dynamics_stencil)
  if (use_temperature) then
    call create_group_pass(pass_uv_T_S_h, CS%tv%T, G%Domain, halo=dynamics_stencil)
    call create_group_pass(pass_uv_T_S_h, CS%tv%S, G%Domain, halo=dynamics_stencil)
  endif
  call create_group_pass(pass_uv_T_S_h, CS%h, G%Domain, halo=dynamics_stencil)

  call do_group_pass(pass_uv_T_S_h, G%Domain)
  if (associated(CS%tv%p_surf)) call pass_var(CS%tv%p_surf, G%Domain, halo=dynamics_stencil)
  call cpu_clock_end(id_clock_pass_init)

  ! Update derived thermodynamic quantities.
  if (allocated(CS%tv%SpV_avg)) then
    call calc_derived_thermo(CS%tv, CS%h, G, GV, US, halo=dynamics_stencil, debug=CS%debug)
  endif


  diag => CS%diag
  ! Initialize the diag mediator.
  call diag_mediator_init(G, GV, US, GV%ke, param_file, diag, doc_file_dir=dirs%output_directory)
  if (associated(CS%OBC)) then
    call diag_mediator_set_OBC_info(G, CS%OBC%segnum_u, CS%OBC%segnum_v, diag)
  endif
  if (present(diag_ptr)) diag_ptr => CS%diag

  ! Initialize the diagnostics masks for native arrays.
  ! This step has to be done after call to MOM_initialize_state
  ! and before MOM_diagnostics_init
  call diag_masks_set(G, GV%ke, diag)

  ! Set up pointers within diag mediator control structure,
  ! this needs to occur _after_ CS%h etc. have been allocated.
  call diag_set_state_ptrs(CS%h, CS%tv, diag)

  ! This call sets up the diagnostic axes. These are needed,
  ! e.g. to generate the target grids below.
  call set_axes_info(G, GV, US, param_file, diag)

  ! Whenever thickness/T/S changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.  In non-Boussinesq mode,
  ! calc_derived_thermo needs to be called before diag_update_remap_grids.
  call diag_update_remap_grids(diag)

  ! Setup the diagnostic grid storage types
  call diag_grid_storage_init(CS%diag_pre_sync, G, GV, diag)
  call diag_grid_storage_init(CS%diag_pre_dyn, G, GV, diag)

  ! Calculate masks for diagnostics arrays in non-native coordinates
  ! This step has to be done after set_axes_info() because the axes needed
  ! to be configured, and after diag_update_remap_grids() because the grids
  ! must be defined.
  call set_masks_for_axes(G, diag)

  ! Register the volume cell measure (must be one of first diagnostics)
  call register_cell_measure(G, CS%diag, Time)

  call cpu_clock_begin(id_clock_MOM_init)
  ! Diagnose static fields AND associate areas/volumes with axes
  call write_static_fields(G, GV, US, CS%tv, CS%diag)
  call callTree_waypoint("static fields written (initialize_MOM)")

  if (CS%use_ALE_algorithm) then
    call ALE_writeCoordinateFile( CS%ALE_CSp, GV, dirs%output_directory )
    call callTree_waypoint("ALE initialized (initialize_MOM)")
  elseif (write_geom_files) then
    call write_vertgrid_file(GV, US, param_file, dirs%output_directory)
  endif
  call cpu_clock_end(id_clock_MOM_init)

  if (CS%use_dbclient) call database_comms_init(param_file, CS%dbcomms_CS)
  CS%useMEKE = MEKE_init(Time, G, GV, US, param_file, diag, CS%dbcomms_CS, CS%MEKE_CSp, CS%MEKE, &
                         restart_CSp, CS%MEKE_in_dynamics)

  call VarMix_init(Time, G, GV, US, param_file, diag, CS%VarMix)
  call set_visc_init(Time, G, GV, US, param_file, diag, CS%visc, CS%set_visc_CSp, restart_CSp, CS%OBC)
  call thickness_diffuse_init(Time, G, GV, US, param_file, diag, CS%CDp, CS%thickness_diffuse_CSp)
  if (CS%interface_filter) &
    call interface_filter_init(Time, G, GV, US, param_file, diag, CS%CDp, CS%interface_filter_CSp)

  new_sim = is_new_run(restart_CSp)
  if (use_temperature) then
    CS%use_stochastic_EOS = MOM_stoch_eos_init(Time, G, GV, US, param_file, diag, CS%stoch_eos_CS, restart_CSp)
  else
    CS%use_stochastic_EOS = .false.
  endif

  if (CS%use_porbar) &
    call porous_barriers_init(Time, GV, US, param_file, diag, CS%por_bar_CS)

  if (CS%split) then
    allocate(eta(SZI_(G),SZJ_(G)), source=0.0)
    if (CS%use_alt_split) then
      call initialize_dyn_split_RK2b(CS%u, CS%v, CS%h, CS%tv, CS%uh, CS%vh, eta, Time, &
              G, GV, US, param_file, diag, CS%dyn_split_RK2b_CSp, CS%HA_CSp, restart_CSp, &
              CS%dt, CS%ADp, CS%CDp, MOM_internal_state, CS%VarMix, CS%MEKE, &
              CS%thickness_diffuse_CSp, CS%OBC, CS%update_OBC_CSp, CS%ALE_CSp, CS%set_visc_CSp, &
              CS%visc, dirs, CS%ntrunc, CS%pbv, calc_dtbt=calc_dtbt, &
              cont_stencil=CS%cont_stencil, dyn_h_stencil=CS%dyn_h_stencil)
    else
      call initialize_dyn_split_RK2(CS%u, CS%v, CS%h, CS%tv, CS%uh, CS%vh, eta, Time, &
              G, GV, US, param_file, diag, CS%dyn_split_RK2_CSp, CS%HA_CSp, restart_CSp, &
              CS%dt, CS%ADp, CS%CDp, MOM_internal_state, CS%VarMix, CS%MEKE, &
              CS%thickness_diffuse_CSp, CS%OBC, CS%update_OBC_CSp, CS%ALE_CSp, CS%set_visc_CSp, &
              CS%visc, dirs, CS%ntrunc, CS%pbv, calc_dtbt=calc_dtbt, &
              cont_stencil=CS%cont_stencil, dyn_h_stencil=CS%dyn_h_stencil)
    endif
    if (CS%dtbt_reset_period > 0.0) then
      CS%dtbt_reset_interval = real_to_time(CS%dtbt_reset_period, unscale=US%T_to_s)
      ! Set dtbt_reset_time to be the next even multiple of dtbt_reset_interval.
      CS%dtbt_reset_time = Time_init + CS%dtbt_reset_interval * &
                                 ((Time - Time_init) / CS%dtbt_reset_interval)
      if ((CS%dtbt_reset_time > Time) .and. calc_dtbt) then
        ! Back up dtbt_reset_time one interval to force dtbt to be calculated,
        ! because the restart was not aligned with the interval to recalculate
        ! dtbt, and dtbt was not read from a restart file.
        CS%dtbt_reset_time = CS%dtbt_reset_time - CS%dtbt_reset_interval
      endif
    endif
  elseif (CS%use_RK2) then
    call initialize_dyn_unsplit_RK2(CS%u, CS%v, CS%h, CS%tv, Time, G, GV,  &
            US, param_file, diag, CS%dyn_unsplit_RK2_CSp,                  &
            CS%ADp, CS%CDp, MOM_internal_state, CS%OBC,                    &
            CS%update_OBC_CSp, CS%ALE_CSp, CS%set_visc_CSp, CS%visc, dirs, &
            CS%ntrunc, cont_stencil=CS%cont_stencil, dyn_h_stencil=CS%dyn_h_stencil)
  else
    call initialize_dyn_unsplit(CS%u, CS%v, CS%h, CS%tv, Time, G, GV,      &
            US, param_file, diag, CS%dyn_unsplit_CSp,                      &
            CS%ADp, CS%CDp, MOM_internal_state, CS%OBC,                    &
            CS%update_OBC_CSp, CS%ALE_CSp, CS%set_visc_CSp, CS%visc, dirs, &
            CS%ntrunc, cont_stencil=CS%cont_stencil, dyn_h_stencil=CS%dyn_h_stencil)
  endif
  CS%dyn_h_stencil = max(2, CS%dyn_h_stencil)

  !Set OBC segment data update period
  if (associated(CS%OBC) .and. CS%dt_obc_seg_period > 0.0) then
    CS%dt_obc_seg_interval = real_to_time(CS%dt_obc_seg_period, unscale=US%T_to_s)
    CS%dt_obc_seg_time = Time + CS%dt_obc_seg_interval
  endif

  call callTree_waypoint("dynamics initialized (initialize_MOM)")

  CS%mixedlayer_restrat = mixedlayer_restrat_init(Time, G, GV, US, param_file, diag, &
                                                  CS%mixedlayer_restrat_CSp, restart_CSp)

  if (GV%Boussinesq .and. associated(CS%visc%h_ML)) then
    ! This is here to allow for a transition of restart files between model versions.
    call get_param(param_file, "MOM", "MLE_USE_PBL_MLD", MLE_use_PBL_MLD, &
                   default=.false., do_not_log=.true.)
    if (MLE_use_PBL_MLD .and. .not.query_initialized(CS%visc%h_ML, "h_ML", restart_CSp) .and. &
        associated(CS%visc%MLD)) then
      do j=js,je ; do i=is,ie ; CS%visc%h_ML(i,j) = GV%Z_to_H * CS%visc%MLD(i,j) ; enddo ; enddo
    endif
  endif

  if (CS%mixedlayer_restrat) then
    if (.not.(bulkmixedlayer .or. CS%use_ALE_algorithm)) &
      call MOM_error(FATAL, "MOM: MIXEDLAYER_RESTRAT true requires a boundary layer scheme.")
    ! When DIABATIC_FIRST=False and using CS%visc%ML in mixedlayer_restrat we need to update after a restart
    if (.not. CS%diabatic_first .and. associated(CS%visc%MLD)) &
      call pass_var(CS%visc%MLD, G%domain, halo=1)
    if (.not. CS%diabatic_first .and. associated(CS%visc%h_ML)) &
      call pass_var(CS%visc%h_ML, G%domain, halo=1)
  endif

  call MOM_diagnostics_init(MOM_internal_state, CS%ADp, CS%CDp, Time, G, GV, US, &
                            param_file, diag, CS%diagnostics_CSp, CS%tv)
  call diag_copy_diag_to_storage(CS%diag_pre_sync, CS%h, CS%diag)


  if (CS%adiabatic) then
    call adiabatic_driver_init(Time, G, param_file, diag, CS%diabatic_CSp, &
                               CS%tracer_flow_CSp)
  else
    call diabatic_driver_init(Time, G, GV, US, param_file, CS%use_ALE_algorithm, diag, &
                              CS%ADp, CS%CDp, CS%diabatic_CSp, CS%tracer_flow_CSp, &
                              CS%sponge_CSp, CS%ALE_sponge_CSp, CS%oda_incupd_CSp, CS%int_tide_CSp)
  endif

  CS%vertex_shear = kappa_shear_at_vertex(param_file)

  ! GMM, the following is needed to get BLDs into the dynamics module
  if (CS%split .and. fpmix) then
    call init_dyn_split_RK2_diabatic(CS%diabatic_CSp, CS%dyn_split_RK2_CSp)
  endif

  if (associated(CS%sponge_CSp)) &
    call init_sponge_diags(Time, G, GV, US, diag, CS%sponge_CSp)

  if (associated(CS%oda_incupd_CSp)) &
    call init_oda_incupd_diags(Time, G, GV, diag, CS%oda_incupd_CSp, US)

  call tracer_advect_init(Time, G, US, param_file, diag, CS%tracer_adv_CSp)
  call tracer_hor_diff_init(Time, G, GV, US, param_file, diag, CS%tv%eqn_of_state, CS%diabatic_CSp, &
                            CS%tracer_diff_CSp)

  call lock_tracer_registry(CS%tracer_Reg)
  call callTree_waypoint("tracer registry now locked (initialize_MOM)")

  ! now register some diagnostics since the tracer registry is now locked
  call register_surface_diags(Time, G, US, CS%sfc_IDs, CS%diag, CS%tv)
  call register_diags(Time, G, GV, US, CS%IDs, CS%diag)
  call register_transport_diags(Time, G, GV, US, CS%transport_IDs, CS%diag)
  call extract_diabatic_member(CS%diabatic_CSp, use_KPP=use_KPP)
  call register_tracer_diagnostics(CS%tracer_Reg, CS%h, Time, diag, G, GV, US, &
                                   CS%use_ALE_algorithm, use_KPP)
  if (CS%use_ALE_algorithm) then
    call ALE_register_diags(Time, G, GV, US, diag, CS%ALE_CSp)
  endif

  ! Do any necessary halo updates on any auxiliary variables that have been initialized.
  call cpu_clock_begin(id_clock_pass_init)
  if (associated(CS%visc%Kv_shear)) &
    call pass_var(CS%visc%Kv_shear, G%Domain, To_All+Omit_Corners, halo=1)

  if (associated(CS%visc%Kv_slow)) &
    call pass_var(CS%visc%Kv_slow, G%Domain, To_All+Omit_Corners, halo=1)
  call cpu_clock_end(id_clock_pass_init)

  ! This subroutine initializes any tracer packages.
  call tracer_flow_control_init(.not.new_sim, Time, G, GV, US, CS%h, param_file, &
             CS%diag, CS%OBC, CS%tracer_flow_CSp, CS%sponge_CSp, &
             CS%ALE_sponge_CSp, CS%tv)
  if (present(tracer_flow_CSp)) tracer_flow_CSp => CS%tracer_flow_CSp

  if (associated(CS%ALE_sponge_CSp)) &
    call init_ALE_sponge_diags(Time, G, diag, CS%ALE_sponge_CSp, US)

  ! If running in offline tracer mode, initialize the necessary control structure and
  ! parameters
  if (present(offline_tracer_mode)) offline_tracer_mode=CS%offline_tracer_mode

  if (CS%offline_tracer_mode) then
    ! Setup some initial parameterizations and also assign some of the subtypes
    call offline_transport_init(param_file, CS%offline_CSp, CS%diabatic_CSp, G, GV, US)
    call insert_offline_main( CS=CS%offline_CSp, ALE_CSp=CS%ALE_CSp, diabatic_CSp=CS%diabatic_CSp, &
                              diag=CS%diag, OBC=CS%OBC, tracer_adv_CSp=CS%tracer_adv_CSp, &
                              tracer_flow_CSp=CS%tracer_flow_CSp, tracer_Reg=CS%tracer_Reg, &
                              tv=CS%tv, x_before_y=(MODULO(first_direction,2)==0), debug=CS%debug )
    call register_diags_offline_transport(Time, CS%diag, CS%offline_CSp, GV, US)
  endif

  if (associated(CS%OBC)) then
    ! At this point any information related to the tracer reservoirs has either been read from
    ! the restart file or has been specified in the segments.  Initialize the tracer reservoir
    ! values from the segments if they have not been set via the restart file.
    call setup_OBC_tracer_reservoirs(G, GV, CS%OBC, restart_CSp)
    call setup_OBC_thickness_reservoirs(G, GV, CS%OBC, restart_CSp)
    call open_boundary_halo_update(G, CS%OBC)
  endif

  call register_obsolete_diagnostics(param_file, CS%diag)

  if (use_frazil) then
    if (.not.query_initialized(CS%tv%frazil, "frazil", restart_CSp)) then
      CS%tv%frazil(:,:) = 0.0
      call set_initialized(CS%tv%frazil, "frazil", restart_CSp)
    endif
  endif

  if (CS%interp_p_surf) then
    CS%p_surf_prev_set = query_initialized(CS%p_surf_prev, "p_surf_prev", restart_CSp)

    if (CS%p_surf_prev_set) then
      call pass_var(CS%p_surf_prev, G%domain)
    endif
  endif

  if (.not.query_initialized(CS%ave_ssh_ibc, "ave_ssh", restart_CSp)) then
    if (CS%split) then
      call find_eta(CS%h, CS%tv, G, GV, US, CS%ave_ssh_ibc, eta, dZref=G%Z_ref)
    else
      call find_eta(CS%h, CS%tv, G, GV, US, CS%ave_ssh_ibc, dZref=G%Z_ref)
    endif
    call set_initialized(CS%ave_ssh_ibc, "ave_ssh", restart_CSp)
  endif
  if (CS%split) deallocate(eta)

  CS%nstep_tot = 0
  if (present(count_calls)) CS%count_calls = count_calls
  call MOM_sum_output_init(G_in, GV, US, param_file, dirs%output_directory, &
                           CS%ntrunc, Time_init, CS%sum_output_CSp)

  ! Flag whether to save initial conditions in finish_MOM_initialization() or not.
  CS%write_IC = save_IC .and. &
                .not.((dirs%input_filename(1:1) == 'r') .and. &
                      (LEN_TRIM(dirs%input_filename) == 1))

  if (CS%ensemble_ocean) then
    call init_oda(Time, G, GV, US, CS%diag, CS%odaCS)
  endif

  ! initialize stochastic physics
  call stochastics_init(CS%dt_therm, CS%G, CS%GV, CS%stoch_CS, param_file, diag, Time)

  call callTree_leave("initialize_MOM()")
  call cpu_clock_end(id_clock_init) ; call cpu_clock_end(id_clock_ocean)

end procedure initialize_MOM
module procedure finish_MOM_initialization
  type(ocean_grid_type),   pointer :: G => NULL()  ! pointer to a structure containing
  type(verticalGrid_type), pointer :: GV => NULL() ! Pointer to the vertical grid structure
  type(unit_scale_type),   pointer :: US => NULL() ! Pointer to a structure containing
  type(MOM_restart_CS),    pointer :: restart_CSp_tmp => NULL()
  real, allocatable :: z_interface(:,:,:) ! Interface heights [Z ~> m]
  call cpu_clock_begin(id_clock_init)
  call callTree_enter("finish_MOM_initialization()")

  ! Pointers for convenience
  G => CS%G ; GV => CS%GV ; US => CS%US

  if (CS%use_particles) then
    call particles_init(CS%particles, G, CS%Time, CS%dt_therm, CS%u, CS%v, CS%h)
  endif

  ! Write initial conditions
  if (CS%write_IC) then
    allocate(restart_CSp_tmp)
    restart_CSp_tmp = CS%restart_CS
    call restart_registry_lock(restart_CSp_tmp, unlocked=.true.)
    allocate(z_interface(SZI_(G),SZJ_(G),SZK_(GV)+1))
    call find_eta(CS%h, CS%tv, G, GV, US, z_interface, dZref=G%Z_ref)
    call register_restart_field(z_interface, "eta", .true., restart_CSp_tmp, &
                                "Interface heights", "meter", z_grid='i', conversion=US%Z_to_m)
    ! NOTE: write_ic=.true. routes routine to fms2 IO write_initial_conditions interface
    call save_restart(dirs%output_directory, Time, CS%G_in, &
                      restart_CSp_tmp, filename=CS%IC_file, GV=GV, write_ic=.true.)
    deallocate(z_interface)
    deallocate(restart_CSp_tmp)
  endif

  call write_energy(CS%u, CS%v, CS%h, CS%tv, Time, 0, G, GV, US, &
                    CS%sum_output_CSp, CS%tracer_flow_CSp)

  call callTree_leave("finish_MOM_initialization()")
  call cpu_clock_end(id_clock_init)

end procedure finish_MOM_initialization
module procedure register_diags
  character(len=48) :: thickness_units
  thickness_units = get_thickness_units(GV)

  ! Diagnostics of the rapidly varying dynamic state
  IDs%id_u = register_diag_field('ocean_model', 'u_dyn', diag%axesCuL, Time, &
      'Zonal velocity after the dynamics update', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_v = register_diag_field('ocean_model', 'v_dyn', diag%axesCvL, Time, &
      'Meridional velocity after the dynamics update', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_h = register_diag_field('ocean_model', 'h_dyn', diag%axesTL, Time, &
      'Layer Thickness after the dynamics update', thickness_units, conversion=GV%H_to_MKS, &
      v_extensive=.true.)
  IDs%id_ssh_inst = register_diag_field('ocean_model', 'SSH_inst', diag%axesT1, &
      Time, 'Instantaneous Sea Surface Height', 'm', conversion=US%Z_to_m)

end procedure register_diags
module procedure MOM_timing_init
  id_clock_dynamics = cpu_clock_id('Ocean dynamics', grain=CLOCK_SUBCOMPONENT)
  id_clock_thermo   = cpu_clock_id('Ocean thermodynamics and tracers', grain=CLOCK_SUBCOMPONENT)
  id_clock_remap    = cpu_clock_id('Ocean grid generation and remapping', grain=CLOCK_SUBCOMPONENT)
  id_clock_other    = cpu_clock_id('Ocean Other', grain=CLOCK_SUBCOMPONENT)
  id_clock_MOM_end  = cpu_clock_id('Ocean MOM_end', grain=CLOCK_SUBCOMPONENT)
  id_clock_tracer   = cpu_clock_id('(Ocean tracer advection)', grain=CLOCK_MODULE_DRIVER)
  if (.not.CS%adiabatic) then
    id_clock_diabatic = cpu_clock_id('(Ocean diabatic driver)', grain=CLOCK_MODULE_DRIVER)
  else
    id_clock_adiabatic = cpu_clock_id('(Ocean adiabatic driver)', grain=CLOCK_MODULE_DRIVER)
  endif

  id_clock_continuity = cpu_clock_id('(Ocean continuity equation *)', grain=CLOCK_MODULE)
  id_clock_BBL_visc   = cpu_clock_id('(Ocean set BBL viscosity)', grain=CLOCK_MODULE)
  id_clock_pass       = cpu_clock_id('(Ocean message passing *)', grain=CLOCK_MODULE)
  id_clock_MOM_init   = cpu_clock_id('(Ocean MOM_initialize_state)', grain=CLOCK_MODULE)
  id_clock_pass_init  = cpu_clock_id('(Ocean init message passing *)', grain=CLOCK_ROUTINE)
  if (CS%thickness_diffuse) &
    id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion *)', grain=CLOCK_MODULE)
  if (CS%interface_filter) &
    id_clock_int_filter = cpu_clock_id('(Ocean interface height filter *)', grain=CLOCK_MODULE)
 !if (CS%mixedlayer_restrat) &
    id_clock_ml_restrat = cpu_clock_id('(Ocean mixed layer restrat)', grain=CLOCK_MODULE)
  id_clock_diagnostics  = cpu_clock_id('(Ocean collective diagnostics)', grain=CLOCK_MODULE)
  id_clock_Z_diag       = cpu_clock_id('(Ocean Z-space diagnostics)', grain=CLOCK_MODULE)
  id_clock_ALE          = cpu_clock_id('(Ocean ALE)', grain=CLOCK_MODULE)
  if (CS%offline_tracer_mode) then
    id_clock_offline_tracer = cpu_clock_id('Ocean offline tracers', grain=CLOCK_SUBCOMPONENT)
  endif
  id_clock_stoch = cpu_clock_id('(Stochastic EOS)', grain=CLOCK_MODULE)
  id_clock_varT = cpu_clock_id('(SGS Temperature Variance)', grain=CLOCK_MODULE)

  id_clock_save_restart   = cpu_clock_id('(Ocean MOM save_restart)', grain=CLOCK_MODULE)

end procedure MOM_timing_init
module procedure set_restart_fields
  logical :: use_ice_shelf ! Needed to determine whether to add CS%Hml to restarts
  character(len=48) :: thickness_units, flux_units
  type(vardesc) :: u_desc, v_desc
  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

  if (associated(CS%tv%T)) &
    call register_restart_field(CS%tv%T, "Temp", .true., restart_CSp, &
                                "Potential Temperature", "degC", conversion=US%C_to_degC)
  if (associated(CS%tv%S)) &
    call register_restart_field(CS%tv%S, "Salt", .true., restart_CSp, &
                                "Salinity", "PPT", conversion=US%S_to_ppt)

  call register_restart_field(CS%h, "h", .true., restart_CSp, &
                              "Layer Thickness", thickness_units, conversion=GV%H_to_MKS)

  u_desc = var_desc("u", "m s-1", "Zonal velocity", hor_grid='Cu')
  v_desc = var_desc("v", "m s-1", "Meridional velocity", hor_grid='Cv')
  call register_restart_pair(CS%u, CS%v, u_desc, v_desc, .true., restart_CSp, conversion=US%L_T_to_m_s)

  if (associated(CS%tv%frazil)) &
    call register_restart_field(CS%tv%frazil, "frazil", .false., restart_CSp, &
                                "Frazil heat flux into ocean", &
                                "J m-2", conversion=US%Q_to_J_kg*US%RZ_to_kg_m2)

  if (CS%interp_p_surf) then
    call register_restart_field(CS%p_surf_prev, "p_surf_prev", .false., restart_CSp, &
                                "Previous ocean surface pressure", "Pa", conversion=US%RL2_T2_to_Pa)
  endif

  if (associated(CS%tv%p_surf)) &
    call register_restart_field(CS%tv%p_surf, "p_surf_EOS", .false., restart_CSp, &
                                "Ocean surface pressure used in EoS", "Pa", conversion=US%RL2_T2_to_Pa)

  call register_restart_field(CS%ave_ssh_ibc, "ave_ssh", .false., restart_CSp, &
                              "Time average sea surface height", "meter", conversion=US%Z_to_m)

  ! hML is needed when using the ice shelf module
  call get_param(param_file, '', "ICE_SHELF", use_ice_shelf, default=.false., &
                 do_not_log=.true.)
  if (use_ice_shelf .and. associated(CS%Hml)) then
    call register_restart_field(CS%Hml, "hML", .false., restart_CSp, &
                                "Mixed layer thickness", "m", conversion=US%Z_to_m)
  endif

  ! Register scalar unit conversion factors.
  call register_restart_field(CS%first_dir_restart, "First_direction", .false., restart_CSp, &
                              "Indicator of the first direction in split calculations.", "nondim")

end procedure set_restart_fields
module procedure adjust_ssh_for_p_atm
  real :: Rho_conv(SZI_(G))  ! The density used to convert surface pressure to
  real :: IgR0        ! The SSH conversion factor from R L2 T-2 to Z [Z T2 R-1 L-2 ~> m Pa-1].
  logical :: calc_rho
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  EOSdom(:) = EOS_domain(G%HI)
  if (associated(p_atm)) then
    calc_rho = use_EOS .and. associated(tv%eqn_of_state)
    ! Correct the output sea surface height for the contribution from the ice pressure.
    do j=js,je
      if (calc_rho) then
        call calculate_density(tv%T(:,j,1), tv%S(:,j,1), 0.5*p_atm(:,j), Rho_conv, &
                               tv%eqn_of_state, EOSdom)
        do i=is,ie
          IgR0 = 1.0 / (Rho_conv(i) * GV%g_Earth)
          ssh(i,j) = ssh(i,j) + p_atm(i,j) * IgR0
        enddo
      else
        IgR0 = 1.0 / (GV%Rho0 * GV%g_Earth)
        do i=is,ie
          ssh(i,j) = ssh(i,j) + p_atm(i,j) * IgR0
        enddo
      endif
    enddo
  endif

end procedure adjust_ssh_for_p_atm
module procedure extract_surface_state
  real :: hu, hv  ! Thicknesses interpolated to velocity points [H ~> m or kg m-2]
  type(ocean_grid_type),   pointer :: G => NULL() !< pointer to a structure containing
  type(ocean_grid_type),   pointer :: G_in => NULL() !< Input grid metric
  type(verticalGrid_type), pointer :: GV => NULL() !< structure containing vertical grid info
  type(unit_scale_type),   pointer :: US => NULL() !< structure containing various unit conversion factors
  type(surface),           pointer :: sfc_state => NULL()  ! surface state on the model grid
  real, dimension(:,:,:),  pointer :: h => NULL()    !< h : layer thickness [H ~> m or kg m-2]
  real :: depth(SZI_(CS%G))  !< Distance from the surface in depth units [Z ~> m] or [H ~> m or kg m-2]
  real :: depth_ml           !< Depth over which to average to determine mixed
  real :: dh                 !< Thickness of a layer within the mixed layer [Z ~> m] or [H ~> m or kg m-2]
  real :: mass               !< Mass per unit area of a layer [R Z ~> kg m-2]
  real :: I_depth            !< The inverse of depth [Z-1 ~> m-1] or [H-1 ~> m-1 or m2 kg-1]
  real :: missing_depth      !< The portion of depth_ml that can not be found in a column [H ~> m or kg m-2]
  real :: H_rescale          !< A conversion factor from thickness units to the units used in the
  real :: T_freeze(SZI_(CS%G)) !< freezing temperature [C ~> degC]
  real :: pres(SZI_(CS%G))   !< Pressure to use for the freezing temperature calculation [R L2 T-2 ~> Pa]
  real :: delT(SZI_(CS%G))   !< Depth integral of T-T_freeze [H C ~> m degC or degC kg m-2]
  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  integer :: i, j, k, is, ie, js, je, nz, numberOfErrors, ig, jg
  integer :: isd, ied, jsd, jed
  integer :: iscB, iecB, jscB, jecB, isdB, iedB, jsdB, jedB
  logical :: localError
  logical :: use_iceshelves
  character(240) :: msg
  integer :: turns    ! Number of quarter turns
  call callTree_enter("extract_surface_state(), MOM.F90")
  G => CS%G ; G_in => CS%G_in ; GV => CS%GV ; US => CS%US
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  iscB = G%iscB ; iecB = G%iecB ; jscB = G%jscB ; jecB = G%jecB
  isdB = G%isdB ; iedB = G%iedB ; jsdB = G%jsdB ; jedB = G%jedB
  h => CS%h

  use_temperature = associated(CS%tv%T)

  use_iceshelves=.false.
  if (associated(CS%frac_shelf_h)) use_iceshelves = .true.

  turns = 0
  if (CS%rotate_index) &
    turns = G%HI%turns

  if (.not.sfc_state_in%arrays_allocated)  then
    !  Consider using a run-time flag to determine whether to do the vertical
    ! integrals, since the 3-d sums are not negligible in cost.
    call allocate_surface_state(sfc_state_in, G_in, use_temperature, &
          do_integrals=.true., omit_frazil=.not.associated(CS%tv%frazil),&
          use_iceshelves=use_iceshelves)
  endif

  if (CS%rotate_index) then
    allocate(sfc_state)
    call allocate_surface_state(sfc_state, G, use_temperature, &
              do_integrals=.true., omit_frazil=.not.associated(CS%tv%frazil),&
              use_iceshelves=use_iceshelves, sfc_state_in=sfc_state_in, turns=turns)
  else
    sfc_state => sfc_state_in
  endif

  sfc_state%T_is_conT = CS%tv%T_is_conT
  sfc_state%S_is_absS = CS%tv%S_is_absS

  do j=js,je ; do i=is,ie
    sfc_state%sea_lev(i,j) = CS%ave_ssh_ibc(i,j)
  enddo ; enddo

  if (allocated(sfc_state%frazil) .and. associated(CS%tv%frazil)) then ; do j=js,je ; do i=is,ie
    sfc_state%frazil(i,j) = CS%tv%frazil(i,j)
  enddo ; enddo ; endif

  ! copy Hml into sfc_state, so that caps can access it
  do j=js,je ; do i=is,ie
    sfc_state%Hml(i,j) = CS%Hml(i,j)
  enddo ; enddo

  if (CS%Hmix < 0.0) then  ! A bulk mixed layer is in use, so layer 1 has the properties
    if (use_temperature) then ; do j=js,je ; do i=is,ie
      sfc_state%SST(i,j) = CS%tv%T(i,j,1)
      sfc_state%SSS(i,j) = CS%tv%S(i,j,1)
    enddo ; enddo ; endif
    do j=js,je ; do I=is-1,ie
      sfc_state%u(I,j) = CS%u(I,j,1)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      sfc_state%v(i,J) = CS%v(i,J,1)
    enddo ; enddo

  else  ! (CS%Hmix >= 0.0)
    H_rescale = 1.0
    depth_ml = CS%Hmix
    if (CS%answer_date < 20190101) then
      H_rescale = GV%H_to_Z
      depth_ml = GV%H_to_Z*CS%Hmix
    endif
    ! Determine the mean tracer properties of the uppermost depth_ml fluid.

    !$OMP parallel do default(shared) private(depth,dh)
    do j=js,je
      do i=is,ie
        depth(i) = 0.0
        if (use_temperature) then
          sfc_state%SST(i,j) = 0.0 ; sfc_state%SSS(i,j) = 0.0
        else
          sfc_state%sfc_density(i,j) = 0.0
        endif
      enddo

      do k=1,nz ; do i=is,ie
        if (depth(i) + h(i,j,k)*H_rescale < depth_ml) then
          dh = h(i,j,k)*H_rescale
        elseif (depth(i) < depth_ml) then
          dh = depth_ml - depth(i)
        else
          dh = 0.0
        endif
        if (use_temperature) then
          sfc_state%SST(i,j) = sfc_state%SST(i,j) + dh * CS%tv%T(i,j,k)
          sfc_state%SSS(i,j) = sfc_state%SSS(i,j) + dh * CS%tv%S(i,j,k)
        else
          sfc_state%sfc_density(i,j) = sfc_state%sfc_density(i,j) + dh * GV%Rlay(k)
        endif
        depth(i) = depth(i) + dh
      enddo ; enddo
  ! Calculate the average properties of the mixed layer depth.
      do i=is,ie
        if (CS%answer_date < 20190101) then
          if (depth(i) < GV%H_subroundoff*H_rescale) &
              depth(i) = GV%H_subroundoff*H_rescale
          if (use_temperature) then
            sfc_state%SST(i,j) = sfc_state%SST(i,j) / depth(i)
            sfc_state%SSS(i,j) = sfc_state%SSS(i,j) / depth(i)
          else
            sfc_state%sfc_density(i,j) = sfc_state%sfc_density(i,j) / depth(i)
          endif
        else
          if (depth(i) < GV%H_subroundoff*H_rescale) then
            I_depth = 1.0 / (GV%H_subroundoff*H_rescale)
            missing_depth = GV%H_subroundoff*H_rescale - depth(i)
            if (use_temperature) then
              sfc_state%SST(i,j) = (sfc_state%SST(i,j) + missing_depth*CS%tv%T(i,j,1)) * I_depth
              sfc_state%SSS(i,j) = (sfc_state%SSS(i,j) + missing_depth*CS%tv%S(i,j,1)) * I_depth
            else
              sfc_state%sfc_density(i,j) = (sfc_state%sfc_density(i,j) + &
                                            missing_depth*GV%Rlay(1)) * I_depth
            endif
          else
            I_depth = 1.0 / depth(i)
            if (use_temperature) then
              sfc_state%SST(i,j) = sfc_state%SST(i,j) * I_depth
              sfc_state%SSS(i,j) = sfc_state%SSS(i,j) * I_depth
            else
              sfc_state%sfc_density(i,j) = sfc_state%sfc_density(i,j) * I_depth
            endif
          endif
        endif
      enddo
    enddo ! end of j loop

!   Determine the mean velocities in the uppermost depth_ml fluid.
    ! NOTE: Velocity loops start on `[ij]s-1` in order to update halo values
    !       required by the speed diagnostic on the non-symmetric grid.
    !       This assumes that u and v halos have already been updated.
    if (CS%Hmix_UV>0.) then
      depth_ml = CS%Hmix_UV
      if (CS%answer_date < 20190101) depth_ml = GV%H_to_Z*CS%Hmix_UV
      !$OMP parallel do default(shared) private(depth,dh,hv)
      do J=js-1,je
        do i=is,ie
          depth(i) = 0.0
          sfc_state%v(i,J) = 0.0
        enddo
        do k=1,nz ; do i=is,ie
          hv = 0.5 * (h(i,j,k) + h(i,j+1,k)) * H_rescale
          if (depth(i) + hv < depth_ml) then
            dh = hv
          elseif (depth(i) < depth_ml) then
            dh = depth_ml - depth(i)
          else
            dh = 0.0
          endif
          sfc_state%v(i,J) = sfc_state%v(i,J) + dh * CS%v(i,J,k)
          depth(i) = depth(i) + dh
        enddo ; enddo
        ! Calculate the average properties of the mixed layer depth.
        do i=is,ie
          sfc_state%v(i,J) = sfc_state%v(i,J) / max(depth(i), GV%H_subroundoff*H_rescale)
        enddo
      enddo ! end of j loop

      !$OMP parallel do default(shared) private(depth,dh,hu)
      do j=js,je
        do I=is-1,ie
          depth(I) = 0.0
          sfc_state%u(I,j) = 0.0
        enddo
        do k=1,nz ; do I=is-1,ie
          hu = 0.5 * (h(i,j,k) + h(i+1,j,k)) * H_rescale
          if (depth(i) + hu < depth_ml) then
            dh = hu
          elseif (depth(I) < depth_ml) then
            dh = depth_ml - depth(I)
          else
            dh = 0.0
          endif
          sfc_state%u(I,j) = sfc_state%u(I,j) + dh * CS%u(I,j,k)
          depth(I) = depth(I) + dh
        enddo ; enddo
        ! Calculate the average properties of the mixed layer depth.
        do I=is-1,ie
          sfc_state%u(I,j) = sfc_state%u(I,j) / max(depth(I), GV%H_subroundoff*H_rescale)
        enddo
      enddo ! end of j loop
    else ! Hmix_UV<=0.
      do j=js,je ; do I=is-1,ie
        sfc_state%u(I,j) = CS%u(I,j,1)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        sfc_state%v(i,J) = CS%v(i,J,1)
      enddo ; enddo
    endif
  endif  ! (CS%Hmix >= 0.0)


  if (allocated(sfc_state%melt_potential)) then
    !$OMP parallel do default(shared) private(depth_ml, dh, T_freeze, depth, pres, delT)
    do j=js,je
      do i=is,ie
        depth(i) = 0.0
        delT(i) = 0.0
        pres(i) = 0.0
        ! Here it is assumed that p=0 is OK, since HFrz ~ 10 to 20m, but under ice-shelves this
        ! can be a very bad assumption.  ###To fix this, uncomment the following...
        !   pres(i) = p_surface(i) + 0.5*(GV%g_Earth*GV%H_to_RZ)*h(i,j,1)
      enddo

      do k=1,nz
        call calculate_TFreeze(CS%tv%S(is:ie,j,k), pres(is:ie), T_freeze(is:ie), CS%tv%eqn_of_state)
        do i=is,ie
          depth_ml = min(CS%HFrz, CS%visc%h_ML(i,j))
          if (depth(i) + h(i,j,k) < depth_ml) then
            dh = h(i,j,k)
          elseif (depth(i) < depth_ml) then
            dh = depth_ml - depth(i)
          else
            dh = 0.0
          endif

          depth(i) = depth(i) + dh
          delT(i) =  delT(i) + dh * (CS%tv%T(i,j,k) - T_freeze(i))
        enddo
      ! If there is a pressure-dependent freezing point calculation uncomment the following.
      ! if (k<nz) then ; do i=is,ie
      !   pres(i) = pres(i) + 0.5*(GV%g_Earth*GV%H_to_RZ) * (h(i,j,k) + h(i,j,k+1))
      ! enddo ; endif
      enddo

      do i=is,ie
        ! set melt_potential to zero to avoid passing previous values
        sfc_state%melt_potential(i,j) = 0.0

        if (G%mask2dT(i,j)>0.) then
          ! instantaneous melt_potential [Q R Z ~> J m-2]
          sfc_state%melt_potential(i,j) = CS%tv%C_p * GV%H_to_RZ * delT(i)
        endif
      enddo
    enddo ! end of j loop
  endif   ! melt_potential

  if (allocated(sfc_state%taux_shelf) .and. allocated(CS%visc%taux_shelf)) then
    !$OMP parallel do default(shared)
    do j=js,je ; do I=is-1,ie
      sfc_state%taux_shelf(I,j) = CS%visc%taux_shelf(I,j)
    enddo ; enddo
  endif
  if (allocated(sfc_state%tauy_shelf) .and. allocated(CS%visc%tauy_shelf)) then
    !$OMP parallel do default(shared)
    do J=js-1,je ; do i=is,ie
      sfc_state%tauy_shelf(i,J) = CS%visc%tauy_shelf(i,J)
    enddo ; enddo
  endif

  if (allocated(sfc_state%ocean_mass) .and. allocated(sfc_state%ocean_heat) .and. &
      allocated(sfc_state%ocean_salt)) then
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      sfc_state%ocean_mass(i,j) = 0.0
      sfc_state%ocean_heat(i,j) = 0.0 ; sfc_state%ocean_salt(i,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared) private(mass)
    do j=js,je ; do k=1,nz ; do i=is,ie
      mass = GV%H_to_RZ*h(i,j,k)
      sfc_state%ocean_mass(i,j) = sfc_state%ocean_mass(i,j) + mass
      sfc_state%ocean_heat(i,j) = sfc_state%ocean_heat(i,j) + mass * CS%tv%T(i,j,k)
      sfc_state%ocean_salt(i,j) = sfc_state%ocean_salt(i,j) + mass * (1.0e-3*CS%tv%S(i,j,k))
    enddo ; enddo ; enddo
  else
    if (allocated(sfc_state%ocean_mass)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_mass(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nz ; do i=is,ie
        sfc_state%ocean_mass(i,j) = sfc_state%ocean_mass(i,j) + GV%H_to_RZ*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (allocated(sfc_state%ocean_heat)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_heat(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared) private(mass)
      do j=js,je ; do k=1,nz ; do i=is,ie
        mass = GV%H_to_RZ*h(i,j,k)
        sfc_state%ocean_heat(i,j) = sfc_state%ocean_heat(i,j) + mass * CS%tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (allocated(sfc_state%ocean_salt)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_salt(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared) private(mass)
      do j=js,je ; do k=1,nz ; do i=is,ie
        mass = GV%H_to_RZ*h(i,j,k)
        sfc_state%ocean_salt(i,j) = sfc_state%ocean_salt(i,j) + mass * (1.0e-3*CS%tv%S(i,j,k))
      enddo ; enddo ; enddo
    endif
  endif

  if (associated(CS%tracer_flow_CSp)) then
    call call_tracer_surface_state(sfc_state, h, G, GV, US, CS%tracer_flow_CSp)
  endif

  if (CS%check_bad_sfc_vals) then
    numberOfErrors=0 ! count number of errors
    do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j)>0.) then
        localError = sfc_state%sea_lev(i,j) < -G%bathyT(i,j) - G%Z_ref &
                .or. sfc_state%sea_lev(i,j) >=  CS%bad_val_ssh_max + (G%meanSL(i,j) - G%Z_ref) &
                .or. sfc_state%sea_lev(i,j) <= -CS%bad_val_ssh_max + (G%meanSL(i,j) - G%Z_ref) &
                .or. sfc_state%sea_lev(i,j) + G%bathyT(i,j) + G%Z_ref < CS%bad_val_col_thick
        if (use_temperature) localError = localError &
                .or. sfc_state%SSS(i,j)<0.                        &
                .or. sfc_state%SSS(i,j)>=CS%bad_val_sss_max       &
                .or. sfc_state%SST(i,j)< CS%bad_val_sst_min       &
                .or. sfc_state%SST(i,j)>=CS%bad_val_sst_max
        if (localError) then
          numberOfErrors=numberOfErrors+1
          if (numberOfErrors<9) then ! Only report details for the first few errors
            ig = i + G%HI%idg_offset ! Global i-index
            jg = j + G%HI%jdg_offset ! Global j-index
            if (use_temperature) then
              write(msg(1:240),'(2(a,I0,1x),4(a,f8.3,1x),8(a,es11.4,1x))') &
                'Extreme surface sfc_state detected: i=',ig,'j=',jg, &
                'lon=',G%geoLonT(i,j), 'lat=',G%geoLatT(i,j), &
                'x=',G%gridLonT(ig), 'y=',G%gridLatT(jg), &
                'D=',US%Z_to_m*(G%bathyT(i,j)+G%Z_ref), 'SSH=',US%Z_to_m*sfc_state%sea_lev(i,j), &
                'SST=',US%C_to_degC*sfc_state%SST(i,j), 'SSS=',US%S_to_ppt*sfc_state%SSS(i,j), &
                'U-=',US%L_T_to_m_s*sfc_state%u(I-1,j), 'U+=',US%L_T_to_m_s*sfc_state%u(I,j), &
                'V-=',US%L_T_to_m_s*sfc_state%v(i,J-1), 'V+=',US%L_T_to_m_s*sfc_state%v(i,J)
            else
              write(msg(1:240),'(2(a,I0,1x),4(a,f8.3,1x),6(a,es11.4))') &
                'Extreme surface sfc_state detected: i=',ig,'j=',jg, &
                'lon=',G%geoLonT(i,j), 'lat=',G%geoLatT(i,j), &
                'x=',G%gridLonT(ig), 'y=',G%gridLatT(jg), &
                'D=',US%Z_to_m*(G%bathyT(i,j)+G%Z_ref), 'SSH=',US%Z_to_m*sfc_state%sea_lev(i,j), &
                'U-=',US%L_T_to_m_s*sfc_state%u(I-1,j), 'U+=',US%L_T_to_m_s*sfc_state%u(I,j), &
                'V-=',US%L_T_to_m_s*sfc_state%v(i,J-1), 'V+=',US%L_T_to_m_s*sfc_state%v(i,J)
            endif
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          elseif (numberOfErrors==9) then ! Indicate once that there are more errors
            call MOM_error(WARNING, 'There were more unreported extreme events!', all_print=.true.)
          endif ! numberOfErrors
        endif ! localError
      endif ! mask2dT
    enddo ; enddo
    call sum_across_PEs(numberOfErrors)
    if (numberOfErrors>0) then
      write(msg(1:240),'(a,i0,a)') 'There were a total of ',numberOfErrors, &
          ' locations detected with extreme surface values!'
      call MOM_error(FATAL, trim(msg))
    endif
  endif

  if (CS%debug) call MOM_surface_chksum("Post extract_sfc", sfc_state, G, US, haloshift=0, symmetric=.true.)

  ! Rotate sfc_state back onto the input grid, sfc_state_in
  if (CS%rotate_index) then
    call rotate_surface_state(sfc_state, sfc_state_in, G_in, -turns)
    call deallocate_surface_state(sfc_state)
  endif

  call callTree_leave("extract_surface_sfc_state()")
end procedure extract_surface_state
module procedure rotate_initial_state
  call rotate_vector(u_in, v_in, turns, u, v)
  call rotate_array(h_in, turns, h)
  if (use_temperature) then
    call rotate_array(T_in, turns, T)
    call rotate_array(S_in, turns, S)
  endif
end procedure rotate_initial_state
module procedure MOM_state_is_synchronized
  logical :: adv_only
  adv_only = .false. ; if (present(adv_dyn)) adv_only = adv_dyn

  if (adv_only) then
    in_synch = (CS%t_dyn_rel_adv == 0.0)
  else
    in_synch = ((CS%t_dyn_rel_adv == 0.0) .and. (CS%t_dyn_rel_thermo == 0.0))
  endif

end procedure MOM_state_is_synchronized
module procedure get_MOM_state_elements
  if (present(G)) G => CS%G_in
  if (present(GV)) GV => CS%GV
  if (present(US)) US => CS%US
  if (present(C_p)) C_p = CS%US%Q_to_J_kg*US%degC_to_C * CS%tv%C_p
  if (present(C_p_scaled)) C_p_scaled = CS%tv%C_p
  if (present(use_temp)) use_temp = associated(CS%tv%T)
end procedure get_MOM_state_elements
module procedure get_ocean_stocks
  if (present(mass)) &
    mass = global_mass_integral(CS%h, CS%G, CS%GV, on_PE_only=on_PE_only)
  if (present(heat)) &
    heat = CS%US%Q_to_J_kg*CS%US%RZL2_to_kg * CS%tv%C_p * &
           global_mass_integral(CS%h, CS%G, CS%GV, CS%tv%T, on_PE_only=on_PE_only, tmp_scale=CS%US%C_to_degC)
  if (present(salt)) &
    salt = 1.0e-3 * global_mass_integral(CS%h, CS%G, CS%GV, CS%tv%S, on_PE_only=on_PE_only, unscale=CS%US%S_to_ppt)

end procedure get_ocean_stocks
module procedure save_MOM_restart
  logical :: showCallTree
  showCallTree = callTree_showQuery()

  call cpu_clock_begin(id_clock_ocean) ; call cpu_clock_begin(id_clock_save_restart)
  if (showCallTree) call callTree_waypoint("About to call save_restart (step_MOM)")
  call save_restart(directory, time, G, CS%restart_CS, &
      time_stamped=time_stamped, filename=filename, GV=GV, &
      num_rest_files=num_rest_files, write_IC=write_IC)
  if (showCallTree) call callTree_waypoint("Done with call to save_restart (step_MOM)")

  if (CS%use_particles) call particles_save_restart(CS%particles, CS%h, directory, time, time_stamped)
  call cpu_clock_end(id_clock_save_restart) ; call cpu_clock_end(id_clock_ocean)
end procedure save_MOM_restart
module procedure MOM_end
  call cpu_clock_begin(id_clock_ocean) ; call cpu_clock_begin(id_clock_MOM_end)

  call MOM_sum_output_end(CS%sum_output_CSp)

  if (CS%use_ALE_algorithm) call ALE_end(CS%ALE_CSp)

  !deallocate porous topography variables
  deallocate(CS%pbv%por_face_areaU) ; deallocate(CS%pbv%por_face_areaV)
  deallocate(CS%pbv%por_layer_widthU) ; deallocate(CS%pbv%por_layer_widthV)

  ! NOTE: Allocated in PressureForce_FV_Bouss
  if (associated(CS%tv%varT)) deallocate(CS%tv%varT)

  call tracer_advect_end(CS%tracer_adv_CSp)
  call tracer_hor_diff_end(CS%tracer_diff_CSp)
  call tracer_registry_end(CS%tracer_Reg)
  call tracer_flow_control_end(CS%tracer_flow_CSp)

  if (.not. CS%adiabatic) then
    call diabatic_driver_end(CS%diabatic_CSp)
    deallocate(CS%diabatic_CSp)
  endif

  call MOM_diagnostics_end(CS%diagnostics_CSp, CS%ADp, CS%CDp)

  if (CS%offline_tracer_mode) call offline_transport_end(CS%offline_CSp)

  if (CS%split .and. CS%use_alt_split) then
    call end_dyn_split_RK2b(CS%dyn_split_RK2b_CSp)
  elseif (CS%split) then
    call end_dyn_split_RK2(CS%dyn_split_RK2_CSp)
  elseif (CS%use_RK2) then
    call end_dyn_unsplit_RK2(CS%dyn_unsplit_RK2_CSp)
  else
    call end_dyn_unsplit(CS%dyn_unsplit_CSp)
  endif

  if (CS%use_particles) then
    call particles_end(CS%particles, CS%h)
    deallocate(CS%particles)
  endif

  call thickness_diffuse_end(CS%thickness_diffuse_CSp, CS%CDp)
  if (CS%interface_filter) call interface_filter_end(CS%interface_filter_CSp, CS%CDp)
  call VarMix_end(CS%VarMix)
  call set_visc_end(CS%visc, CS%set_visc_CSp)
  call MEKE_end(CS%MEKE)

  if (associated(CS%tv%internal_heat)) deallocate(CS%tv%internal_heat)
  if (associated(CS%tv%TempxPmE)) deallocate(CS%tv%TempxPmE)

  DEALLOC_(CS%ave_ssh_ibc) ; DEALLOC_(CS%ssh_rint) ; DEALLOC_(CS%eta_av_bc)

  ! TODO: debug_truncations deallocation

  DEALLOC_(CS%uhtr) ; DEALLOC_(CS%vhtr)

  if (associated(CS%Hml)) deallocate(CS%Hml)
  if (associated(CS%tv%salt_deficit)) deallocate(CS%tv%salt_deficit)
  if (associated(CS%tv%frazil)) deallocate(CS%tv%frazil)
  if (allocated(CS%tv%SpV_avg)) deallocate(CS%tv%SpV_avg)

  if (associated(CS%tv%T)) then
    DEALLOC_(CS%T) ; CS%tv%T => NULL() ; DEALLOC_(CS%S) ; CS%tv%S => NULL()
  endif

  DEALLOC_(CS%u) ; DEALLOC_(CS%v) ; DEALLOC_(CS%h)
  DEALLOC_(CS%uh) ; DEALLOC_(CS%vh)

  if (associated(CS%update_OBC_CSp)) call OBC_register_end(CS%update_OBC_CSp)
  if (associated(CS%OBC)) call open_boundary_end(CS%OBC)

  call verticalGridEnd(CS%GV)
  call MOM_grid_end(CS%G)

  if (CS%debug .or. CS%G%symmetric) &
    call deallocate_MOM_domain(CS%G%Domain_aux)

  if (CS%rotate_index) &
    call deallocate_MOM_domain(CS%G%Domain)

  ! The MPP domains may be needed by an external coupler, so use `cursory`.
  ! TODO: This may create a domain memory leak, and needs investigation.
  call deallocate_MOM_domain(CS%G_in%domain, cursory=.true.)

  call unit_scaling_end(CS%US)

  call cpu_clock_end(id_clock_MOM_end) ; call cpu_clock_end(id_clock_ocean)

end procedure MOM_end
end submodule MOM_s
