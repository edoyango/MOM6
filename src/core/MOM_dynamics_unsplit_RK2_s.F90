submodule (MOM_dynamics_unsplit_RK2) MOM_dynamics_unsplit_RK2_s
#include <MOM_memory.h>
  implicit none
contains
module procedure step_MOM_dyn_unsplit_RK2
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: h_av ! Averaged layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: hp ! Predicted layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: dz ! Distance between the interfaces around a layer [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: up ! Predicted zonal velocities [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vp ! Predicted meridional velocities [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: ueffA   ! Effective Area of U-Faces [H L ~> m2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: veffA   ! Effective Area of V-Faces [H L ~> m2]
  real, dimension(:,:), pointer :: p_surf => NULL() ! A pointer to the surface pressure [R L2 T-2 ~> Pa]
  real :: dt_pred   ! The time step for the predictor part of the baroclinic time stepping [T ~> s]
  real :: dt_visc   ! The time step for a part of the update due to viscosity [T ~> s]
  logical :: dyn_p_surf
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: cor_stencil  ! Stencil size for Coriolis schemes [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  dt_pred = dt * CS%BE
  cor_stencil = CoriolisAdv_stencil(CS%CoriolisAdv)

  h_av(:,:,:) = 0 ; hp(:,:,:) = 0
  up(:,:,:) = 0
  vp(:,:,:) = 0

  dyn_p_surf = associated(p_surf_begin) .and. associated(p_surf_end)
  if (dyn_p_surf) then
    call safe_alloc_ptr(p_surf,G%isd,G%ied,G%jsd,G%jed) ; p_surf(:,:) = 0.0
  else
    p_surf => forces%p_surf
  endif

! Runge-Kutta second order accurate two step scheme is used to step
! all of the fields except h.  h is stepped separately.

  if (CS%debug) then
    call MOM_state_chksum("Start Predictor ", u_in, v_in, h_in, uh, vh, G, GV, US)
  endif

! diffu = horizontal viscosity terms (u,h)
  call enable_averages(dt,Time_local, CS%diag)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_in, v_in, h_in, uh, vh, CS%diffu, CS%diffv, MEKE, VarMix, &
                            G, GV, US, CS%hor_visc, tv, dt, STOCH=STOCH)
  call cpu_clock_end(id_clock_horvisc)
  call disable_averaging(CS%diag)
  call pass_vector(CS%diffu, CS%diffv, G%Domain, clock=id_clock_pass)

! This continuity step is solely for the Coroilis terms, specifically in the
! denominator of PV and in the mass transport or PV.
! uh = u[n-1]*h[n-1/2]
! hp = h[n-1/2] + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  ! This is a duplicate calculation of the last continuity from the previous step
  ! and could/should be optimized out. -AJA
  call continuity(u_in, v_in, h_in, hp, uh, vh, dt_pred, G, GV, US, CS%continuity_CSp, CS%OBC, pbv)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(hp, G%Domain, halo=cor_stencil, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, halo=cor_stencil, clock=id_clock_pass)
  if (cor_stencil > 2) then
    call pass_vector(u_in, v_in, G%Domain, halo=cor_stencil, clock=id_clock_pass)
  endif

! h_av = (h + hp)/2  (used in PV denominator)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js-cor_stencil,je+cor_stencil ; do i=is-cor_stencil,ie+cor_stencil
    h_av(i,j,k) = (h_in(i,j,k) + hp(i,j,k)) * 0.5
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

! CAu = -(f+zeta)/h_av vh + d/dx KE  (function of u[n-1] and uh[n-1])
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_in, v_in, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, US, CS%CoriolisAdv, pbv)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)  (function of h[n-1/2])
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.5*p_surf_begin(i,j) + 0.5*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_in, tv, CS%PFu, CS%PFv, G, GV, US, CS%PressureForce_CSp, CS%ALE_CSp, CS%ADp, p_surf)
  call cpu_clock_end(id_clock_pres)
  call pass_vector(CS%PFu, CS%PFv, G%Domain, clock=id_clock_pass)
  call pass_vector(CS%CAu, CS%CAv, G%Domain, clock=id_clock_pass)

  if (associated(CS%OBC)) then ; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, US, tv, h_in, CS%update_OBC_CSp, Time_local)
  endif ; endif
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, CS%PFu, CS%PFv)
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, CS%CAu, CS%CAv)
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, CS%diffu, CS%diffv)
  endif

! up+[n-1/2] = u[n-1] + dt_pred * (PFu + CAu)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) + dt_pred * &
                   ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt_pred * &
                   ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) &
    call MOM_accel_chksum("Predictor 1 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G, GV, US)

 ! up[n-1/2] <- up*[n-1/2] + dt/2 d/dz visc d/dz up[n-1/2]
  call cpu_clock_begin(id_clock_vertvisc)
  call enable_averages(dt, Time_local, CS%diag)
  dt_visc = dt ; if (CS%dt_visc_bug) dt_visc = dt_pred
  call set_viscous_ML(u_in, v_in, h_av, tv, forces, visc, dt_visc, G, GV, US, CS%set_visc_CSp)
  call disable_averaging(CS%diag)

  call thickness_to_dz(h_av, tv, dz, G, GV, US, halo_size=1)
  call vertvisc_coef(up, vp, h_av, dz, forces, visc, tv, dt_pred, G, GV, US, CS%vertvisc_CSp, CS%OBC, VarMix)
  call vertvisc(up, vp, h_av, forces, visc, dt_pred, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, US, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call pass_vector(up, vp, G%Domain, clock=id_clock_pass)

! uh = up[n-1/2] * h[n-1/2]
! h_av = h + dt div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h_in, hp, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(hp, G%Domain, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, clock=id_clock_pass)

! h_av <- (h + hp)/2   (centered at n-1/2)
  do k=1,nz ; do j=js-cor_stencil,je+cor_stencil ; do i=is-cor_stencil,ie+cor_stencil
    h_av(i,j,k) = (h_in(i,j,k) + hp(i,j,k)) * 0.5
  enddo ; enddo ; enddo

  if (CS%debug) &
    call MOM_state_chksum("Predictor 1", up, vp, h_av, uh, vh, G, GV, US)

! CAu = -(f+zeta(up))/h_av vh + d/dx KE(up)  (function of up[n-1/2], h[n-1/2])
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(up, vp, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, US, CS%CoriolisAdv, pbv)
  call cpu_clock_end(id_clock_Cor)
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, CS%CAu, CS%CAv)
  endif

! call enable_averages(dt, Time_local, CS%diag)  ?????????????????????/

! up* = u[n] + (1+gamma) * dt * ( PFu + CAu )  Extrapolated for damping
! u*[n+1] = u[n] + dt * ( PFu + CAu )
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) + dt * (1.+CS%begw) * &
            ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
    u_in(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) +  dt * &
            ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt * (1.+CS%begw) * &
            ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
    v_in(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt * &
            ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
  enddo ; enddo ; enddo

! up[n] <- up* + dt d/dz visc d/dz up
! u[n] <- u*[n] + dt d/dz visc d/dz u[n]
  call cpu_clock_begin(id_clock_vertvisc)
  call thickness_to_dz(h_av, tv, dz, G, GV, US, halo_size=1)
  call vertvisc_coef(up, vp, h_av, dz, forces, visc, tv, dt, G, GV, US, CS%vertvisc_CSp, CS%OBC, VarMix)
  call vertvisc(up, vp, h_av, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call vertvisc_coef(u_in, v_in, h_av, dz, forces, visc, tv, dt, G, GV, US, CS%vertvisc_CSp, CS%OBC, VarMix)
  call vertvisc(u_in, v_in, h_av, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp,&
                G, GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call cpu_clock_end(id_clock_vertvisc)
  call pass_vector(up, vp, G%Domain, clock=id_clock_pass)
  call pass_vector(u_in, v_in, G%Domain, clock=id_clock_pass)

! uh = up[n] * h[n]  (up[n] might be extrapolated to damp GWs)
! h[n+1] = h[n] + dt div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h_in, h_in, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(h_in, G%Domain, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, clock=id_clock_pass)

! Accumulate mass flux for tracer transport
  do k=1,nz
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      uhtr(I,j,k) = uhtr(I,j,k) + dt*uh(I,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      vhtr(i,J,k) = vhtr(i,J,k) + dt*vh(i,J,k)
    enddo ; enddo
  enddo

  if (CS%debug) then
    call MOM_state_chksum("Corrector", u_in, v_in, h_in, uh, vh, G, GV, US)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV, US)
  endif

  if (GV%Boussinesq) then
    do j=js,je ; do i=is,ie ; eta_av(i,j) = -GV%Z_to_H*G%bathyT(i,j) ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie ; eta_av(i,j) = 0.0 ; enddo ; enddo
  endif
  do k=1,nz ; do j=js,je ; do i=is,ie
    eta_av(i,j) = eta_av(i,j) + h_av(i,j,k)
  enddo ; enddo ; enddo

  if (dyn_p_surf) deallocate(p_surf)

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%PFu, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%PFv, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%CAu, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%CAv, CS%diag)

!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)

! Calculate effective areas and post data
  if (CS%id_ueffA > 0) then
    ueffA(:,:,:) = 0
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      if (abs(up(I,j,k)) > 0.) ueffA(I,j,k) = uh(I,j,k)/up(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_ueffA, ueffA, CS%diag)
  endif

  if (CS%id_veffA > 0) then
    veffA(:,:,:) = 0
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      if (abs(vp(i,J,k)) > 0.) veffA(i,J,k) = vh(i,J,k)/vp(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_veffA, veffA, CS%diag)
  endif


end procedure step_MOM_dyn_unsplit_RK2
module procedure register_restarts_dyn_unsplit_RK2
  character(len=48) :: thickness_units, flux_units
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

! This is where a control structure that is specific to this module would be allocated.
  if (associated(CS)) then
    call MOM_error(WARNING, "register_restarts_dyn_unsplit_RK2 called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)

  ALLOC_(CS%diffu(IsdB:IedB,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC_(CS%diffv(isd:ied,JsdB:JedB,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC_(CS%CAu(IsdB:IedB,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC_(CS%CAv(isd:ied,JsdB:JedB,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC_(CS%PFu(IsdB:IedB,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC_(CS%PFv(isd:ied,JsdB:JedB,nz)) ; CS%PFv(:,:,:) = 0.0

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

!  No extra restart fields are needed with this time stepping scheme.

end procedure register_restarts_dyn_unsplit_RK2
module procedure initialize_dyn_unsplit_RK2
  character(len=40) :: mdl = "MOM_dynamics_unsplit_RK2" ! This module's name.
  character(len=48) :: flux_units
# include "version_variable.h"
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "initialize_dyn_unsplit_RK2 called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "initialize_dyn_unsplit_RK2 called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "BE", CS%be, &
                 "If SPLIT is true, BE determines the relative weighting "//&
                 "of a 2nd-order Runga-Kutta baroclinic time stepping "//&
                 "scheme (0.5) and a backward Euler scheme (1) that is "//&
                 "used for the Coriolis and inertial terms.  BE may be "//&
                 "from 0.5 to 1, but instability may occur near 0.5. "//&
                 "BE is also applicable if SPLIT is false and USE_RK2 "//&
                 "is true.", units="nondim", default=0.6)
  call get_param(param_file, mdl, "BEGW", CS%begw, &
                 "If SPLIT is true, BEGW is a number from 0 to 1 that "//&
                 "controls the extent to which the treatment of gravity "//&
                 "waves is forward-backward (0) or simulated backward "//&
                 "Euler (1).  0 is almost always used. "//&
                 "If SPLIT is false and USE_RK2 is true, BEGW can be "//&
                 "between 0 and 0.5 to damp gravity waves.", &
                 units="nondim", default=0.0)

  call get_param(param_file, mdl, "UNSPLIT_DT_VISC_BUG", CS%dt_visc_bug, &
                 "If false, use the correct timestep in the viscous terms applied in the first "//&
                 "predictor step with the unsplit time stepping scheme, and in the calculation "//&
                 "of the turbulent mixed layer properties for viscosity with unsplit or "//&
                 "unsplit_RK2.  If true, an older incorrect value is used.", &
                 default=.false.)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "TIDES", CS%use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%use_tides)

  allocate(CS%taux_bot(IsdB:IedB,jsd:jed), source=0.0)
  allocate(CS%tauy_bot(isd:ied,JsdB:JedB), source=0.0)

  MIS%diffu => CS%diffu ; MIS%diffv => CS%diffv
  MIS%PFu => CS%PFu ; MIS%PFv => CS%PFv
  MIS%CAu => CS%CAu ; MIS%CAv => CS%CAv

  CS%ADp => Accel_diag ; CS%CDp => Cont_diag
  Accel_diag%diffu => CS%diffu ; Accel_diag%diffv => CS%diffv
  Accel_diag%PFu => CS%PFu ; Accel_diag%PFv => CS%PFv
  Accel_diag%CAu => CS%CAu ; Accel_diag%CAv => CS%CAv

  call continuity_init(Time, G, GV, US, param_file, diag, CS%continuity_CSp, CS%OBC)
  cont_stencil = continuity_stencil(CS%continuity_CSp)
  call CoriolisAdv_init(Time, G, GV, US, param_file, diag, CS%ADp, CS%CoriolisAdv)
  dyn_h_stencil = max(cont_stencil, CoriolisAdv_stencil(CS%CoriolisAdv))
  if (CS%calculate_SAL) call SAL_init(h, tv, G, GV, US, param_file, CS%SAL_CSp)
  if (CS%use_tides) call tidal_forcing_init(Time, G, US, param_file, CS%tides_CSp)
  call PressureForce_init(Time, G, GV, US, param_file, diag, CS%PressureForce_CSp, CS%ADp, &
                          CS%SAL_CSp, CS%tides_CSp)
  call hor_visc_init(Time, G, GV, US, param_file, diag, CS%hor_visc)
  call vertvisc_init(MIS, Time, G, GV, US, param_file, diag, CS%ADp, dirs, &
                     ntrunc, CS%vertvisc_CSp)
  CS%set_visc_CSp => set_visc

  if (associated(ALE_CSp)) CS%ALE_CSp => ALE_CSp
  if (associated(OBC)) CS%OBC => OBC

  flux_units = get_flux_units(GV)
  CS%id_uh = register_diag_field('ocean_model', 'uh', diag%axesCuL, Time, &
      'Zonal Thickness Flux', flux_units, conversion=GV%H_to_MKS*US%L_to_m**2*US%s_to_T, &
      y_cell_method='sum', v_extensive=.true.)
  CS%id_vh = register_diag_field('ocean_model', 'vh', diag%axesCvL, Time, &
      'Meridional Thickness Flux', flux_units, conversion=GV%H_to_MKS*US%L_to_m**2*US%s_to_T, &
      x_cell_method='sum', v_extensive=.true.)
  CS%id_CAu = register_diag_field('ocean_model', 'CAu', diag%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', diag%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', diag%axesCuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', diag%axesCvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_ueffA = register_diag_field('ocean_model', 'ueffA', diag%axesCuL, Time, &
       'Effective U-Face Area', 'm^2', conversion=GV%H_to_m*US%L_to_m, &
       y_cell_method='sum', v_extensive=.true.)
  CS%id_veffA = register_diag_field('ocean_model', 'veffA', diag%axesCvL, Time, &
       'Effective V-Face Area', 'm^2', conversion=GV%H_to_m*US%L_to_m, &
       x_cell_method='sum', v_extensive=.true.)

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)

end procedure initialize_dyn_unsplit_RK2
module procedure end_dyn_unsplit_RK2
  DEALLOC_(CS%diffu) ; DEALLOC_(CS%diffv)
  DEALLOC_(CS%CAu)   ; DEALLOC_(CS%CAv)
  DEALLOC_(CS%PFu)   ; DEALLOC_(CS%PFv)

  if (CS%calculate_SAL) call SAL_end(CS%SAL_CSp)
  if (CS%use_tides) call tidal_forcing_end(CS%tides_CSp)

  deallocate(CS)
end procedure end_dyn_unsplit_RK2
end submodule MOM_dynamics_unsplit_RK2_s
