submodule (DOME_initialization) DOME_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure DOME_initialize_topography
  real :: min_depth  ! The minimum ocean depth [Z ~> m]
  real :: shelf_depth ! The ocean depth on the shelf in the DOME configuration [Z ~> m]
  real :: slope      ! The bottom slope in the DOME configuration [Z L-1 ~> nondim]
  real :: shelf_edge_lat ! The latitude of the edge of the topographic shelf in the same units as geolat, often [km]
  real :: inflow_lon ! The edge longitude of the DOME inflow in the same units as geolon, often [km]
  real :: inflow_width ! The longitudinal width of the DOME inflow channel in the same units as geolat, often [km]
  real :: km_to_grid_unit ! The conversion factor from km to the units of latitude often 1 [nondim],
# include "version_variable.h"
  character(len=40)  :: mdl = "DOME_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "DOME_initialization: "//&
          "DOME_initialize_topography is only set to work with Cartesian axis units.")
  if (abs(G%grid_unit_to_L*US%L_to_m - 1000.0) < 1.0e-3) then ! The grid latitudes are in km.
    km_to_grid_unit = 1.0
  elseif (abs(G%grid_unit_to_L*US%L_to_m - 1.0) < 1.0e-6) then ! The grid latitudes are in m.
    km_to_grid_unit = 1000.0
  else
    call MOM_error(FATAL, "DOME_initialization: "//&
          "DOME_initialize_topography is not recognizing the value of G%grid_unit_to_L.")
  endif

  call MOM_mesg("  DOME_initialization.F90, DOME_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", default=0.0, units="m", scale=US%m_to_Z)
  call get_param(param_file, mdl, "DOME_TOPOG_SLOPE", slope, &
                 "The slope of the bottom topography in the DOME configuration.", &
                 default=0.01, units="nondim", scale=US%L_to_Z)
  call get_param(param_file, mdl, "DOME_SHELF_DEPTH", shelf_depth, &
                 "The bottom depth in the shelf inflow region in the DOME configuration.", &
                 default=600.0, units="m", scale=US%m_to_Z)
  call get_param(param_file, mdl, "DOME_SHELF_EDGE_LAT", shelf_edge_lat, &
                 "The latitude of the shelf edge in the DOME configuration.", &
                 default=600.0, units="km", scale=km_to_grid_unit)
  call get_param(param_file, mdl, "DOME_INFLOW_LON", inflow_lon, &
                 "The edge longitude of the DOME inflow.", units="km", default=1000.0, scale=km_to_grid_unit)
  call get_param(param_file, mdl, "DOME_INFLOW_WIDTH", inflow_width, &
                 "The longitudinal width of the DOME inflow channel.", &
                 units="km", default=100.0, scale=km_to_grid_unit)

  do j=js,je ; do i=is,ie
    if (G%geoLatT(i,j) < shelf_edge_lat) then
      D(i,j) = min(shelf_depth - slope * (G%geoLatT(i,j)-shelf_edge_lat)*G%grid_unit_to_L, max_depth)
    else
      if ((G%geoLonT(i,j) > inflow_lon) .AND. (G%geoLonT(i,j) < inflow_lon+inflow_width)) then
        D(i,j) = shelf_depth
      else
        D(i,j) = 0.5*min_depth
      endif
    endif

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end procedure DOME_initialize_topography
module procedure DOME_initialize_thickness
  real :: e0(SZK_(GV)+1)    ! The resting interface heights [Z ~> m], usually
  real :: eta1D(SZK_(GV)+1) ! Interface height relative to the sea surface
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (just_read) return ! This subroutine has no run-time parameters.

  call MOM_mesg("  DOME_initialization.F90, DOME_initialize_thickness: setting thickness", 5)

  e0(1)=0.0
  do k=2,nz
    e0(K) = -G%max_depth * (real(k-1)-0.5)/real(nz-1)
  enddo

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least an !
!  Angstrom thick, and 2.  the interfaces are where they should be   !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_Z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

end procedure DOME_initialize_thickness
module procedure DOME_initialize_sponges
  real :: eta(SZI_(G),SZJ_(G),SZK_(GV)+1) ! A temporary array for interface heights [Z ~> m].
  real :: temp(SZI_(G),SZJ_(G),SZK_(GV))  ! A temporary array for other variables [various]
  real :: Idamp(SZI_(G),SZJ_(G))          ! The sponge damping rate [T-1 ~> s-1]
  real :: e_tgt(SZK_(GV)+1) ! Target interface heights [Z ~> m].
  real :: min_depth    ! The minimum depth at which to apply damping [Z ~> m]
  real :: damp_W, damp_E ! Damping rates in the western and eastern sponges [T-1 ~> s-1]
  real :: peak_damping ! The maximum sponge damping rates as the edges [T-1 ~> s-1]
  real :: edge_dist    ! The distance to an edge [L ~> m]
  real :: sponge_width ! The width of the sponges [L ~> m]
  character(len=40)  :: mdl = "DOME_initialize_sponges" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "DOME_initialization: "//&
          "DOME_initialize_sponges is only set to work with Cartesian axis units.")

  !   Set up sponges for the DOME configuration
  call get_param(PF, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(PF, mdl, "DOME_SPONGE_DAMP_RATE", peak_damping, &
                 "The largest damping rate in the DOME sponges.", &
                 default=10.0, units="day-1", scale=1.0/(86400.0*US%s_to_T))
  call get_param(PF, mdl, "DOME_SPONGE_WIDTH", sponge_width, &
                 "The width of the DOME sponges.", &
                 default=200.0, units="km", scale=1.0e3*US%m_to_L)

  ! Here the inverse damping time [T-1 ~> s-1], is set. Set Idamp to 0 wherever
  ! there is no sponge, and the subroutines that are called will automatically
  ! set up the sponges only where Idamp is positive and mask2dT is 1.

  Idamp(:,:) = 0.0
  do j=js,je ; do i=is,ie ; if (depth_tot(i,j) > min_depth) then
    edge_dist = (G%geoLonT(i,j) - G%west_lon) * G%grid_unit_to_L
    if (edge_dist < 0.5*sponge_width) then
      damp_W = peak_damping
    elseif (edge_dist < sponge_width) then
      damp_W = peak_damping * (sponge_width - edge_dist) / (0.5*sponge_width)
    else
      damp_W = 0.0
    endif

    edge_dist = ((G%len_lon + G%west_lon) - G%geoLonT(i,j)) * G%grid_unit_to_L
    if (edge_dist < 0.5*sponge_width) then
      damp_E = peak_damping
    elseif (edge_dist < sponge_width) then
      damp_E = peak_damping * (sponge_width - edge_dist) / (0.5*sponge_width)
    else
      damp_E = 0.0
    endif

    Idamp(i,j) = max(damp_W, damp_E)
  endif ; enddo ; enddo

  e_tgt(1) = 0.0
  do K=2,nz ; e_tgt(K) = -(real(K-1)-0.5)*G%max_depth / real(nz-1) ; enddo
  e_tgt(nz+1) = -G%max_depth
  eta(:,:,:) = 0.0
  do K=1,nz+1 ; do j=js,je ; do i=is,ie
    ! These target interface heights will be rescaled inside of apply_sponge, so
    ! they can be in depth space for Boussinesq or non-Boussinesq models.
    eta(i,j,K) = max(e_tgt(K), GV%Angstrom_Z*(nz+1-K) - depth_tot(i,j))
  enddo ; enddo ; enddo

  !  This call stores the sponge damping rates and target interface heights.
  call initialize_sponge(Idamp, eta, G, PF, CSp, GV)

  !   Now register all of the fields which are damped in the sponge.
  ! By default, momentum is advected vertically within the sponge, but
  ! momentum is typically not damped within the layer-mode sponge.

! At this point, the layer-mode DOME configuration is done. The following are here as a
! template for other configurations.

  !  The remaining calls to set_up_sponge_field can be in any order.
  if ( associated(tv%T) ) then
    temp(:,:,:) = 0.0
    call MOM_error(FATAL, "DOME_initialize_sponges is not set up for use with "//&
                          "temperatures defined.")
    ! This should use the target values of T in temp.
    call set_up_sponge_field(temp, tv%T, G, GV, nz, CSp)
    ! This should use the target values of S in temp.
    call set_up_sponge_field(temp, tv%S, G, GV, nz, CSp)
  endif

end procedure DOME_initialize_sponges
module procedure register_DOME_OBC
  if (OBC%number_of_segments /= 1) then
    call MOM_error(FATAL, 'Error in register_DOME_OBC - DOME should have 1 OBC segment', .true.)
  endif

  ! Store this information for use in setting up the OBC restarts for tracer reservoirs.
  OBC%ntr = tr_Reg%ntr
  if (.not. allocated(OBC%tracer_x_reservoirs_used)) then
    allocate(OBC%tracer_x_reservoirs_used(OBC%ntr))
    allocate(OBC%tracer_y_reservoirs_used(OBC%ntr))
    OBC%tracer_x_reservoirs_used(:) = .false.
    OBC%tracer_y_reservoirs_used(:) = .false.
    OBC%tracer_y_reservoirs_used(1) = .true.
  endif

end procedure register_DOME_OBC
module procedure DOME_set_OBC_data
  real :: T0(SZK_(GV))       ! A profile of target temperatures [C ~> degC]
  real :: S0(SZK_(GV))       ! A profile of target salinities [S ~> ppt]
  real :: pres(SZK_(GV))     ! An array of the reference pressure [R L2 T-2 ~> Pa].
  real :: drho_dT(SZK_(GV))  ! Derivative of density with temperature [R C-1 ~> kg m-3 degC-1].
  real :: drho_dS(SZK_(GV))  ! Derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_guess(SZK_(GV)) ! Potential density at T0 & S0 [R ~> kg m-3].
  real :: S_ref             ! A default value for salinities [S ~> ppt]
  real :: T_light           ! A first guess at the temperature of the lightest layer [C ~> degC]
  real :: tr_0              ! The total integrated inflow transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: tr_k              ! The integrated inflow transport of a layer [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: v_k               ! The velocity of a layer at the edge [L T-1 ~> m s-1]
  real :: yt, yb            ! The log of these variables gives the fractional velocities at the
  real :: rst, rsb          ! The relative position of the top and bottom of a layer [nondim],
  real :: rc                ! The relative position of the center of a layer [nondim]
  real :: lon_im1           ! An extrapolated value for the longitude of the western edge of a
  real :: D_edge            ! The thickness [Z ~> m] of the dense fluid at the
  real :: RLay_range        ! The range of densities [R ~> kg m-3].
  real :: Rlay_Ref          ! The surface layer's target density [R ~> kg m-3].
  real :: f_0               ! The reference value of the Coriolis parameter [T-1 ~> s-1]
  real :: f_inflow          ! The value of the Coriolis parameter used to determine DOME inflow
  real :: g_prime_tot       ! The reduced gravity across all layers [L2 Z-1 T-2 ~> m s-2]
  real :: Def_Rad           ! The deformation radius, based on fluid of thickness D_edge [L ~> m]
  real :: inflow_lon        ! The edge longitude of the DOME inflow in the same units as geolon, often [km]
  real :: I_Def_Rad         ! The inverse of the deformation radius in the same units as G%geoLon [km-1]
  real :: Ri_trans          ! The shear Richardson number in the transition
  real :: km_to_grid_unit   ! The conversion factor from km to the units of latitude often 1 [nondim],
  character(len=32)  :: name ! The name of a tracer field.
  character(len=40)  :: mdl = "DOME_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, m, nz, ntherm, ntr_id
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()
  type(tracer_type), pointer      :: tr_ptr => NULL()
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "DOME_initialization: "//&
          "DOME_initialize_topography is only set to work with Cartesian axis units.")
  if (abs(G%grid_unit_to_L*US%L_to_m - 1000.0) < 1.0e-3) then ! The grid latitudes are in km.
    km_to_grid_unit = 1.0
  elseif (abs(G%grid_unit_to_L*US%L_to_m - 1.0) < 1.0e-6) then ! The grid latitudes are in m.
    km_to_grid_unit = 1000.0
  else
    call MOM_error(FATAL, "DOME_initialization: "//&
          "DOME_initialize_topography is not recognizing the value of G%grid_unit_to_L.")
  endif

  call get_param(PF, mdl, "DOME_INFLOW_THICKNESS", D_edge, &
                 "The thickness of the dense DOME inflow at the inner edge.", &
                 default=300.0, units="m", scale=US%m_to_Z)
  call get_param(PF, mdl, "DOME_INFLOW_RI_TRANS", Ri_trans, &
                 "The shear Richardson number in the transition region of the specified "//&
                 "DOME inflow shear profile.", default=(1.0/3.0), units="nondim")
  call get_param(PF, mdl, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0, scale=US%kg_m3_to_R)
  call get_param(PF, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)
  call get_param(PF, mdl, "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the betaplane option.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(PF, mdl, "DOME_INFLOW_F", f_inflow, &
                 "The value of the Coriolis parameter that is used to determine the DOME "//&
                 "inflow properties.", units="s-1", default=f_0*US%s_to_T, scale=US%T_to_s)
  call get_param(PF, mdl, "DOME_INFLOW_LON", inflow_lon, &
                 "The edge longitude of the DOME inflow.", units="km", default=1000.0, scale=km_to_grid_unit)
  if (associated(tv%S) .or. associated(tv%T)) then
    call get_param(PF, mdl, "S_REF", S_ref, &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(PF, mdl, "DOME_T_LIGHT", T_light, &
                 "A first guess at the temperature of the lightest layer in the DOME test case.", &
                 units="degC", default=25.0, scale=US%degC_to_C)
  endif

  if (.not.associated(OBC)) return

  if (GV%Boussinesq) then
    g_prime_tot = (GV%g_Earth / GV%Rho0) * Rlay_range
    Def_Rad = sqrt(D_edge*g_prime_tot) / abs(f_inflow)
    tr_0 = (-D_edge*sqrt(D_edge*g_prime_tot)*0.5*Def_Rad) * GV%Z_to_H
  else
    g_prime_tot = (GV%g_Earth / (Rlay_Ref + 0.5*Rlay_range)) * Rlay_range
    Def_Rad = sqrt(D_edge*g_prime_tot) / abs(f_inflow)
    tr_0 = (-D_edge*sqrt(D_edge*g_prime_tot)*0.5*Def_Rad) * (Rlay_Ref + 0.5*Rlay_range) * GV%RZ_to_H
  endif

  I_Def_Rad = 1.0 / ((1.0e-3*US%L_to_m*km_to_grid_unit) * Def_Rad)
  ! This is mathematically equivalent to
  ! I_Def_Rad = G%grid_unit_to_L / Def_Rad

  if (OBC%number_of_segments /= 1) then
    call MOM_error(WARNING, 'Error in DOME OBC segment setup', .true.)
    return   !!! Need a better error message here
  endif

  segment => OBC%segment(1)
  if (.not. segment%on_pe) return

  ! Set up space for the OBCs to use for all the tracers.
  ntherm = 0
  if (associated(tv%S)) ntherm = ntherm + 1
  if (associated(tv%T)) ntherm = ntherm + 1
  allocate(segment%field(ntherm+tr_Reg%ntr))

  do k=1,nz
    rst = -1.0
    if (k>1) rst = -1.0 + (real(k-1)-0.5)/real(nz-1)

    rsb = 0.0
    if (k<nz) rsb = -1.0 + (real(k-1)+0.5)/real(nz-1)
    rc = -1.0 + real(k-1)/real(nz-1)

    ! These come from assuming geostrophy and a constant Ri profile.
    yt = (2.0*Ri_trans*rst + Ri_trans + 2.0)/(2.0 - Ri_trans)
    yb = (2.0*Ri_trans*rsb + Ri_trans + 2.0)/(2.0 - Ri_trans)
    tr_k = tr_0 * (2.0/(Ri_trans*(2.0-Ri_trans))) * &
           ((log(yt)+1.0)/yt - (log(yb)+1.0)/yb)
    v_k = -sqrt(D_edge*g_prime_tot)*log((2.0 + Ri_trans*(1.0 + 2.0*rc)) / &
                                        (2.0 - Ri_trans))
    if (k == nz)  tr_k = tr_k + tr_0 * (2.0/(Ri_trans*(2.0+Ri_trans))) * &
                                       log((2.0+Ri_trans)/(2.0-Ri_trans))
    ! New way
    isd = segment%HI%isd ; ied = segment%HI%ied
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
    do J=JsdB,JedB ; do i=isd,ied
      ! Here lon_im1 estimates G%geoLonBu(I-1,J), which may not have been set if
      ! the symmetric memory mode is not being used.
      lon_im1 = 2.0*G%geoLonCv(i,J) - G%geoLonBu(I,J)
      segment%normal_trans(i,J,k) = tr_k * (exp(-2.0*(lon_im1 - inflow_lon) * I_Def_Rad) - &
                                  exp(-2.0*(G%geoLonBu(I,J) - inflow_lon) * I_Def_Rad))
      segment%normal_vel(i,J,k) = v_k * exp(-2.0*(G%geoLonCv(i,J) - inflow_lon) * I_Def_Rad)
    enddo ; enddo
  enddo

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (associated(tv%S)) then
    ! In this example, all S inflows have values given by S_ref.
    name = 'salt'
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
    call register_segment_tracer(tr_ptr, ntr_id, PF, GV, segment, OBC_scalar=S_ref, scale=US%ppt_to_S)
  endif
  if (associated(tv%T)) then
    ! In this example, the T values are set to be consistent with the layer
    ! target density and a salinity of S_ref.  This code is taken from
    ! USER_initialize_temp_sal.
    pres(:) = tv%P_Ref ; S0(:) = S_ref ; T0(1) = T_light
    call calculate_density(T0(1), S0(1), pres(1), rho_guess(1), tv%eqn_of_state)
    call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, tv%eqn_of_state, (/1,1/) )

    do k=1,nz ; T0(k) = T0(1) + (GV%Rlay(k)-rho_guess(1)) / drho_dT(1) ; enddo
    do itt=1,6
      call calculate_density(T0, S0, pres, rho_guess, tv%eqn_of_state)
      call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, tv%eqn_of_state)
      do k=1,nz ; T0(k) = T0(k) + (GV%Rlay(k)-rho_guess(k)) / drho_dT(k) ; enddo
    enddo

    ! Temperature is tracer 1 for the OBCs.
    allocate(segment%field(1)%buffer_src(segment%HI%isd:segment%HI%ied,segment%HI%JsdB:segment%HI%JedB,nz))
    do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
      ! With the revised OBC code, buffer_src uses the same rescaled units as for tracers.
      segment%field(1)%buffer_src(i,j,k) = T0(k)
    enddo ; enddo ; enddo
    name = 'temp'
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
    call register_segment_tracer(tr_ptr, ntr_id, PF, GV, segment, OBC_array=.true., scale=US%degC_to_C)
  endif

  ! Set up dye tracers
  ! First dye - only one with OBC values
  ! This field(ntherm+1) requires tr_D1 to be the first tracer after temperature and salinity.
  allocate(segment%field(ntherm+1)%buffer_src(segment%HI%isd:segment%HI%ied,segment%HI%JsdB:segment%HI%JedB,nz))
  do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed ; do i=segment%HI%isd,segment%HI%ied
    if (k < nz/2) then ; segment%field(ntherm+1)%buffer_src(i,j,k) = 0.0
    else ; segment%field(ntherm+1)%buffer_src(i,j,k) = 1.0 ; endif
  enddo ; enddo ; enddo
  name = 'tr_D1'
  call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
  call register_segment_tracer(tr_ptr, ntr_id, PF, GV, OBC%segment(1), OBC_array=.true.)

  ! All tracers but the first have 0 concentration in their inflows. As 0 is the
  ! default value for the inflow concentrations, the following calls are unnecessary.
  do m=2,tr_Reg%ntr
    write(name,'("tr_D",I0)') m
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
    call register_segment_tracer(tr_ptr, ntr_id, PF, GV, OBC%segment(1), OBC_scalar=0.0)
  enddo

end procedure DOME_set_OBC_data
end submodule DOME_initialization_s
