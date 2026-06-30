submodule (MOM_ALE) MOM_ALE_s
#include <MOM_memory.h>
  implicit none
contains
module procedure ALE_init
  character(len=40) :: mdl = "MOM_ALE" ! This module's name.
  character(len=80) :: string, vel_string ! Temporary strings
  real              :: filter_shallow_depth, filter_deep_depth ! Depth ranges of filtering [H ~> m or kg m-2]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical           :: check_reconstruction
  logical           :: check_remapping
  logical           :: force_bounds_in_subcell
  logical           :: local_logical
  logical           :: remap_boundary_extrap
  logical           :: init_boundary_extrap
  logical           :: om4_remap_via_sub_cells
  type(hybgen_regrid_CS), pointer :: hybgen_regridCS => NULL() ! Control structure for hybgen regridding
  real :: h_neglect, h_neglect_edge ! small thicknesses [H ~> m or kg m-2]
  if (associated(CS)) then
    call MOM_error(WARNING, "ALE_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%show_call_tree = callTree_showQuery()
  if (CS%show_call_tree) call callTree_enter("ALE_init(), MOM_ALE.F90")

  call get_param(param_file, mdl, "REMAP_UV_USING_OLD_ALG", CS%remap_uv_using_old_alg, &
                 "If true, uses the old remapping-via-a-delta-z method for "//&
                 "remapping u and v. If false, uses the new method that remaps "//&
                 "between grids described by an old and new thickness.", &
                 default=.false.)

  ! Initialize and configure regridding
  call ALE_initRegridding(G, GV, US, max_depth, param_file, mdl, CS%regridCS)
  call regridding_preadjust_reqs(CS%regridCS, CS%do_conv_adj, CS%use_hybgen_unmix, &
                                 hybgen_CS=hybgen_regridCS)

  ! Initialize and configure remapping that is orchestrated by ALE.
  call get_param(param_file, mdl, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: \n"//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mdl, "VELOCITY_REMAPPING_SCHEME", vel_string, &
                 "This sets the reconstruction scheme used for vertical remapping "//&
                 "of velocities. By default it is the same as REMAPPING_SCHEME. "//&
                 "It can be one of the following schemes: \n"//&
                 trim(remappingSchemesDoc), default=trim(string))
  call get_param(param_file, mdl, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
                 "If true, cell-by-cell reconstructions are checked for "//&
                 "consistency and if non-monotonicity or an inconsistency is "//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mdl, "FATAL_CHECK_REMAPPING", check_remapping, &
                 "If true, the results of remapping are checked for "//&
                 "conservation and new extrema and if an inconsistency is "//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mdl, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
                 "If true, the values on the intermediate grid used for remapping "//&
                 "are forced to be bounded, which might not be the case due to "//&
                 "round off.", default=.false.)
  call get_param(param_file, mdl, "REMAP_BOUNDARY_EXTRAP", remap_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant", default=.false.)
  call get_param(param_file, mdl, "INIT_BOUNDARY_EXTRAP", init_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant during initialization.  "//&
                 "Defaults to REMAP_BOUNDARY_EXTRAP.", default=remap_boundary_extrap)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "This selects the remapping algorithm used in OM4 that does not use "//&
                 "the full reconstruction for the top- and lower-most sub-layers, but instead "//&
                 "assumes they are always vanished (untrue) and so just uses their edge values. "//&
                 "We recommend setting this option to false.", default=.true.)
  call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H * 1.0e-30 ; h_neglect_edge = GV%m_to_H * 1.0e-10
  else
    h_neglect = GV%kg_m2_to_H * 1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H * 1.0e-10
  endif

  call initialize_remapping( CS%remapCS, string, nk=GV%ke, &
                             boundary_extrapolation=init_boundary_extrap, &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell, &
                             om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                             answer_date=CS%answer_date, &
                             h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping( CS%vel_remapCS, vel_string, nk=GV%ke, &
                             boundary_extrapolation=init_boundary_extrap, &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell, &
                             om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                             answer_date=CS%answer_date, &
                             h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

  call get_param(param_file, mdl, "PARTIAL_CELL_VELOCITY_REMAP", CS%partial_cell_vel_remap, &
                 "If true, use partial cell thicknesses at velocity points that are masked out "//&
                 "where they extend below the shallower of the neighboring bathymetry for "//&
                 "remapping velocity.", default=.false.)

  call get_param(param_file, mdl, "REMAP_AFTER_INITIALIZATION", CS%remap_after_initialization, &
                 "If true, applies regridding and remapping immediately after "//&
                 "initialization so that the state is ALE consistent. This is a "//&
                 "legacy step and should not be needed if the initialization is "//&
                 "consistent with the coordinate mode.", default=.true.)

  call get_param(param_file, mdl, "REGRID_USE_DEPTH_BASED_TIME_FILTER", local_logical, &
                 "If true, always uses depth-based time filtering code that updates the "//&
                 "generated grid using REGRID_TIME_SCALE, REGRID_FILTER_SHALLOW_DEPTH, "//&
                 "REGRID_FILTER_DEEP_DEPTH parameters. Setting to True always uses "//&
                 "filtering but setting to False bypasses calculations when filter times = 0.", &
                 default=.true.)
  call set_regrid_params(CS%regridCS, use_depth_based_time_filter=local_logical)
  call get_param(param_file, mdl, "REGRID_TIME_SCALE", CS%regrid_time_scale, &
                 "The time-scale used in blending between the current (old) grid "//&
                 "and the target (new) grid. A short time-scale favors the target "//&
                 "grid (0. or anything less than DT_THERM) has no memory of the old "//&
                 "grid. A very long time-scale makes the model more Lagrangian.", &
                 units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "REGRID_FILTER_SHALLOW_DEPTH", filter_shallow_depth, &
                 "The depth above which no time-filtering is applied. Above this depth "//&
                 "final grid exactly matches the target (new) grid.", &
                 units="m", default=0., scale=GV%m_to_H)
  call get_param(param_file, mdl, "REGRID_FILTER_DEEP_DEPTH", filter_deep_depth, &
                 "The depth below which full time-filtering is applied with time-scale "//&
                 "REGRID_TIME_SCALE. Between depths REGRID_FILTER_SHALLOW_DEPTH and "//&
                 "REGRID_FILTER_DEEP_DEPTH the filter weights adopt a cubic profile.", &
                 units="m", default=0., scale=GV%m_to_H)
  call set_regrid_params(CS%regridCS, depth_of_time_filter_shallow=filter_shallow_depth, &
                         depth_of_time_filter_deep=filter_deep_depth)
  call get_param(param_file, mdl, "REGRID_USE_OLD_DIRECTION", local_logical, &
                 "If true, the regridding integrates upwards from the bottom for "//&
                 "interface positions, much as the main model does. If false "//&
                 "regridding integrates downward, consistent with the remapping code.", &
                 default=.true., do_not_log=.true.)
  call set_regrid_params(CS%regridCS, integrate_downward_for_e=.not.local_logical)

  call get_param(param_file, mdl, "REMAP_VEL_MASK_BBL_THICK", CS%BBL_h_vel_mask, &
                 "A thickness of a bottom boundary layer below which velocities in thin layers "//&
                 "are zeroed out after remapping, following practice with Hybgen remapping, "//&
                 "or a negative value to avoid such filtering altogether.", &
                 default=-0.001, units="m", scale=GV%m_to_H)
  call get_param(param_file, mdl, "REMAP_VEL_MASK_H_THIN", CS%h_vel_mask, &
                 "A thickness at velocity points below which near-bottom layers are zeroed out "//&
                 "after remapping, following practice with Hybgen remapping, "//&
                 "or a negative value to avoid such filtering altogether.", &
                 default=1.0e-6, units="m", scale=GV%m_to_H, do_not_log=(CS%BBL_h_vel_mask<=0.0))

  if (CS%use_hybgen_unmix) &
      call init_hybgen_unmix(CS%hybgen_unmixCS, GV, US, param_file, hybgen_regridCS)

  call get_param(param_file, mdl, "REMAP_VEL_CONSERVE_KE", CS%conserve_ke, &
                 "If true, a correction is applied to the baroclinic component of velocity "//&
                 "after remapping so that total KE is conserved. KE may not be conserved "//&
                 "when (CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)", &
                 default=.false.)
  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  ! Keep a record of values for subsequent queries
  CS%nk = GV%ke

  if (CS%show_call_tree) call callTree_leave("ALE_init()")
end procedure ALE_init
module procedure ALE_set_extrap_boundaries
  logical :: remap_boundary_extrap
  call get_param(param_file, "MOM_ALE", "REMAP_BOUNDARY_EXTRAP", remap_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant", default=.false.)
  call remapping_set_param(CS%remapCS, boundary_extrapolation=remap_boundary_extrap)
end procedure ALE_set_extrap_boundaries
module procedure ALE_set_OM4_remap_algorithm
  call remapping_set_param(CS%remapCS, om4_remap_via_sub_cells=om4_remap_via_sub_cells )

end procedure ALE_set_OM4_remap_algorithm
module procedure ALE_register_diags
  character(len=48)  :: thickness_units
  CS%diag => diag
  thickness_units = get_thickness_units(GV)

  ! These diagnostics of the state variables before ALE are useful for
  ! debugging the ALE code.
  CS%id_u_preale = register_diag_field('ocean_model', 'u_preale', diag%axesCuL, Time, &
      'Zonal velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_preale = register_diag_field('ocean_model', 'v_preale', diag%axesCvL, Time, &
      'Meridional velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_h_preale = register_diag_field('ocean_model', 'h_preale', diag%axesTL, Time, &
      'Layer Thickness before remapping', thickness_units, conversion=GV%H_to_MKS, &
      v_extensive=.true.)
  CS%id_T_preale = register_diag_field('ocean_model', 'T_preale', diag%axesTL, Time, &
      'Temperature before remapping', 'degC', conversion=US%C_to_degC)
  CS%id_S_preale = register_diag_field('ocean_model', 'S_preale', diag%axesTL, Time, &
      'Salinity before remapping', 'PSU', conversion=US%S_to_ppt)
  CS%id_e_preale = register_diag_field('ocean_model', 'e_preale', diag%axesTi, Time, &
      'Interface Heights before remapping', 'm', conversion=US%Z_to_m)

  CS%id_dzRegrid = register_diag_field('ocean_model', 'dzRegrid', diag%axesTi, Time, &
      'Change in interface height due to ALE regridding', 'm', conversion=GV%H_to_m)
  CS%id_vert_remap_h = register_diag_field('ocean_model', 'vert_remap_h', diag%axestl, Time, &
      'layer thicknesses after ALE regridding and remapping', &
      thickness_units, conversion=GV%H_to_MKS, v_extensive=.true.)
  CS%id_vert_remap_h_tendency = register_diag_field('ocean_model', &
      'vert_remap_h_tendency', diag%axestl, Time, &
      'Layer thicknesses tendency due to ALE regridding and remapping', &
      trim(thickness_units)//" s-1", conversion=GV%H_to_MKS*US%s_to_T, v_extensive=.true.)
  CS%id_remap_delta_integ_u2 = register_diag_field('ocean_model', 'ale_u2', diag%axesCu1, Time, &
      'Rate of change in half rho0 times depth integral of squared zonal '//&
      'velocity by remapping. If REMAP_VEL_CONSERVE_KE is .true. then '//&
      'this measures the change before the KE-conserving correction is applied.', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*GV%H_to_RZ*US%L_to_Z**2)
  CS%id_remap_delta_integ_v2 = register_diag_field('ocean_model', 'ale_v2', diag%axesCv1, Time, &
      'Rate of change in half rho0 times depth integral of squared meridional '//&
      'velocity by remapping. If REMAP_VEL_CONSERVE_KE is .true. then '//&
      'this measures the change before the KE-conserving correction is applied.', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*GV%H_to_RZ*US%L_to_Z**2)

end procedure ALE_register_diags
module procedure adjustGridForIntegrity
  call inflate_vanished_layers_old( CS%regridCS, G, GV, h(:,:,:) )

end procedure adjustGridForIntegrity
module procedure ALE_end
  call end_remapping( CS%remapCS )

  if (CS%use_hybgen_unmix) call end_hybgen_unmix( CS%hybgen_unmixCS )
  call end_regridding( CS%regridCS )

  deallocate(CS)

end procedure ALE_end
module procedure pre_ALE_diagnostics
  real :: eta_preale(SZI_(G),SZJ_(G),SZK_(GV)+1)  ! Interface heights before remapping [Z ~> m]
  if (CS%id_u_preale > 0) call post_data(CS%id_u_preale, u,    CS%diag)
  if (CS%id_v_preale > 0) call post_data(CS%id_v_preale, v,    CS%diag)
  if (CS%id_h_preale > 0) call post_data(CS%id_h_preale, h,    CS%diag)
  if (CS%id_T_preale > 0) call post_data(CS%id_T_preale, tv%T, CS%diag)
  if (CS%id_S_preale > 0) call post_data(CS%id_S_preale, tv%S, CS%diag)
  if (CS%id_e_preale > 0) then
    call find_eta(h, tv, G, GV, US, eta_preale, dZref=G%Z_ref)
    call post_data(CS%id_e_preale, eta_preale, CS%diag)
  endif

end procedure pre_ALE_diagnostics
module procedure pre_ALE_adjustments
  integer :: ntr
  if (CS%do_conv_adj) call convective_adjustment(G, GV, h, tv)

  if (CS%use_hybgen_unmix) then
    ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr
    call hybgen_unmix(G, GV, US, CS%hybgen_unmixCS, tv, Reg, ntr, h)
  endif

end procedure pre_ALE_adjustments
module procedure ALE_regrid
  logical :: showCallTree
  showCallTree = callTree_showQuery()

  if (showCallTree) call callTree_enter("ALE_regrid(), MOM_ALE.F90")

  ! Build the new grid and store it in h_new. The old grid is retained as h.
  ! Both are needed for the subsequent remapping of variables.
  dzRegrid(:,:,:) = 0.0
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, US, h, tv, h_new, dzRegrid, &
                        frac_shelf_h=frac_shelf_h, PCM_cell=PCM_cell)

  if (CS%id_dzRegrid>0) then ; if (query_averaging_enabled(CS%diag)) then
    call post_data(CS%id_dzRegrid, dzRegrid, CS%diag, alt_h=h_new)
  endif ; endif

  if (showCallTree) call callTree_leave("ALE_regrid()")

end procedure ALE_regrid
module procedure ALE_offline_inputs
  integer :: nk, i, j, k, isc, iec, jsc, jec
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV))   :: h_new    ! Layer thicknesses after regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: h_src   ! Source grid thicknesses at velocity points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: h_dest  ! Destination grid thicknesses at velocity points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: temp_vec ! Transports on the destination grid [H L2 ~> m3 or kg]
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = GV%ke
  dzRegrid(:,:,:) = 0.0
  h_new(:,:,:) = 0.0

  if (debug) call MOM_tracer_chkinv("Before ALE_offline_inputs", G, GV, h, Reg%Tr, Reg%ntr)

  ! Build new grid from the Zstar state onto the requested vertical coordinate. The new grid is stored
  ! in h_new. The old grid is h. Both are needed for the subsequent remapping of variables. Convective
  ! adjustment right now is not used because it is unclear what to do with vanished layers
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, US, h, tv, h_new, dzRegrid)
  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_offline_inputs)")

  ! Remap all variables from old grid h onto new grid h_new
  call ALE_remap_tracers(CS, G, GV, h, h_new, Reg, debug=CS%show_call_tree)
  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.
  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_inputs)")

  ! Reintegrate mass transports from Zstar to the offline vertical coordinate
  do j=jsc,jec ; do i=G%iscB,G%iecB
    if (G%mask2dCu(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i+1,j,:))
      call reintegrate_column(nk, h_src, uhtr(I,j,:), nk, h_dest, temp_vec)
      uhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo
  do j=G%jscB,G%jecB ; do i=isc,iec
    if (G%mask2dCv(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i,j+1,:))
      call reintegrate_column(nk, h_src, vhtr(I,j,:), nk, h_dest, temp_vec)
      vhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j)>0.) then
      if (check_column_integrals(nk, h_src, nk, h_dest)) then
        call MOM_error(FATAL, "ALE_offline_inputs: Kd interpolation columns do not match")
      endif
      call interpolate_column(nk, h(i,j,:), Kd(i,j,:), nk, h_new(i,j,:), Kd(i,j,:), .true.)
    endif
  enddo ; enddo

  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%T, h_new, tv%T)
  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%S, h_new, tv%S)

  if (debug) call MOM_tracer_chkinv("After ALE_offline_inputs", G, GV, h_new, Reg%Tr, Reg%ntr)

  ! Copy over the new layer thicknesses
  do k = 1,nk  ; do j = jsc-1,jec+1 ; do i = isc-1,iec+1
    h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

  if (CS%show_call_tree) call callTree_leave("ALE_offline_inputs()")
end procedure ALE_offline_inputs
module procedure ALE_regrid_accelerated
  integer :: i, j, itt, nz
  type(thermo_var_ptrs) :: tv_local ! local/intermediate temp/salt
  type(group_pass_type) :: pass_T_S_h ! group pass if the coordinate has a stencil
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))         :: h_loc  ! A working copy of layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))         :: h_orig ! The original layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: T      ! local temporary temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: S      ! local temporary salinities [S ~> ppt]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV))        :: h_old_u ! Source grid thickness at zonal
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV))        :: h_old_v ! Source grid thickness at meridional
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV))        :: h_new_u ! Destination grid thickness at zonal
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV))        :: h_new_v ! Destination grid thickness at meridional
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzInterface ! Interface height changes within
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzIntTotal  ! Cumulative interface position changes [H ~> m or kg m-2]
  nz = GV%ke

  ! initial total interface displacement due to successive regridding
  if (CS%remap_uv_using_old_alg) &
      dzIntTotal(:,:,:) = 0.

  call create_group_pass(pass_T_S_h, T, G%domain)
  call create_group_pass(pass_T_S_h, S, G%domain)
  call create_group_pass(pass_T_S_h, h_loc, G%domain)

  ! copy original temp/salt and set our local tv_pointers to them
  tv_local = tv
  T(:,:,:) = tv%T(:,:,:)
  S(:,:,:) = tv%S(:,:,:)
  tv_local%T => T
  tv_local%S => S

  ! get local copy of thickness and save original state for remapping
  h_loc(:,:,:) = h(:,:,:)
  h_orig(:,:,:) = h(:,:,:)

  ! Apply timescale to regridding (for e.g. filtered_grid_motion)
  if (present(dt)) &
      call ALE_update_regrid_weights(dt, CS)

  do itt = 1, n_itt

    call do_group_pass(pass_T_S_h, G%domain)

    ! generate new grid
    if (CS%do_conv_adj) call convective_adjustment(G, GV, h_loc, tv_local)

    ! Update the layer specific volumes if necessary
    if (allocated(tv_local%SpV_avg)) call calc_derived_thermo(tv_local, h, G, GV, US, halo=1)

    call regridding_main(CS%remapCS, CS%regridCS, G, GV, US, h_loc, tv_local, h, dzInterface)
    if (CS%remap_uv_using_old_alg) &
        dzIntTotal(:,:,:) = dzIntTotal(:,:,:) + dzInterface(:,:,:)

    ! remap from original grid onto new grid
    do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
      call remapping_core_h(CS%remapCS, nz, h_orig(i,j,:), tv%S(i,j,:), nz, h(i,j,:), &
                            tv_local%S(i,j,:))
      call remapping_core_h(CS%remapCS, nz, h_orig(i,j,:), tv%T(i,j,:), nz, h(i,j,:), &
                            tv_local%T(i,j,:))
    enddo ; enddo

    ! starting grid for next iteration
    h_loc(:,:,:) = h(:,:,:)
  enddo

  ! remap all state variables (including those that weren't needed for regridding)
  call ALE_remap_tracers(CS, G, GV, h_orig, h, Reg)

  call ALE_remap_set_h_vel(CS, G, GV, h_orig, h_old_u, h_old_v, OBC)
  if (CS%remap_uv_using_old_alg) then
    call ALE_remap_set_h_vel_via_dz(CS, G, GV, h, h_new_u, h_new_v, OBC, h_orig, dzIntTotal)
  else
    call ALE_remap_set_h_vel(CS, G, GV, h, h_new_u, h_new_v, OBC)
  endif

  call ALE_remap_velocities(CS, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, u, v)

  ! save total dzregrid for diags if needed?
  if (present(dzRegrid)) dzRegrid(:,:,:) = dzIntTotal(:,:,:)

  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

end procedure ALE_regrid_accelerated
module procedure ALE_remap_tracers
  real :: tr_column(GV%ke)  ! A column of updated tracer concentrations [CU ~> Conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_conc ! The rate of change of concentrations [Conc T-1 ~> Conc s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_cont ! The rate of change of cell-integrated tracer
  real, dimension(SZI_(G),SZJ_(G))          :: work_2d ! The rate of change of column-integrated tracer
  logical :: PCM(GV%ke) ! If true, do PCM remapping from a cell.
  real :: Idt           ! The inverse of the timestep [T-1 ~> s-1]
  real :: h1(GV%ke)     ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h2(GV%ke)     ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  logical :: show_call_tree
  type(tracer_type), pointer :: Tr => NULL()
  integer :: i, j, k, m, nz, ntr
  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug

  if (show_call_tree) call callTree_enter("ALE_remap_tracers(), MOM_ALE.F90")

  nz = GV%ke

  ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr

  if (present(dt)) then
    Idt = 1.0/dt
    work_conc(:,:,:) = 0.0
    work_cont(:,:,:) = 0.0
  endif

  ! Remap all registered tracers, including temperature and salinity.
  if (ntr>0) then
    if (show_call_tree) call callTree_waypoint("remapping tracers (ALE_remap_tracers)")
    !$OMP parallel do default(shared) private(h1,h2,tr_column,Tr,PCM,work_conc,work_cont,work_2d)
    do m=1,ntr ! For each tracer
      Tr => Reg%Tr(m)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec ; if (G%mask2dT(i,j)>0.) then
        ! Build the start and final grids
        h1(:) = h_old(i,j,:)
        h2(:) = h_new(i,j,:)
        if (present(PCM_cell)) then
          PCM(:) = PCM_cell(i,j,:)
          call remapping_core_h(CS%remapCS, nz, h1, Tr%t(i,j,:), nz, h2, tr_column, PCM_cell=PCM)
        else
          call remapping_core_h(CS%remapCS, nz, h1, Tr%t(i,j,:), nz, h2, tr_column)
        endif

        ! Possibly underflow any very tiny tracer concentrations to 0.  Note that this is not conservative!
        if (Tr%conc_underflow > 0.0) then ; do k=1,GV%ke
          if (abs(tr_column(k)) < Tr%conc_underflow) tr_column(k) = 0.0
        enddo ; endif

        ! Intermediate steps for tendency of tracer concentration and tracer content.
        if (present(dt)) then
          if (Tr%id_remap_conc > 0) then
            do k=1,GV%ke
              work_conc(i,j,k) = (tr_column(k) - Tr%t(i,j,k)) * Idt
            enddo
          endif
          if (Tr%id_remap_cont > 0 .or. Tr%id_remap_cont_2d > 0) then
            do k=1,GV%ke
              work_cont(i,j,k) = (tr_column(k)*h2(k) - Tr%t(i,j,k)*h1(k)) * Idt
            enddo
          endif
        endif

        ! update tracer concentration
        Tr%t(i,j,:) = tr_column(:)
      endif ; enddo ; enddo

      ! tendency diagnostics.
      if (present(dt)) then
        if (Tr%id_remap_conc > 0) then
          call post_data(Tr%id_remap_conc, work_conc, CS%diag)
        endif
        if (Tr%id_remap_cont > 0) then
          call post_data(Tr%id_remap_cont, work_cont, CS%diag)
        endif

        if (Tr%id_remap_cont_2d > 0) then
          do j = G%jsc,G%jec ; do i = G%isc,G%iec
            work_2d(i,j) = 0.0
            do k = 1,GV%ke
              work_2d(i,j) = work_2d(i,j) + work_cont(i,j,k)
            enddo
          enddo ; enddo
          call post_data(Tr%id_remap_cont_2d, work_2d, CS%diag)
        endif
      endif
    enddo ! m=1,ntr

  endif  ! endif for ntr > 0


  if (CS%id_vert_remap_h > 0) call post_data(CS%id_vert_remap_h, h_old, CS%diag)
  if ((CS%id_vert_remap_h_tendency > 0) .and. present(dt)) then
    do k = 1, nz ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
      work_cont(i,j,k) = (h_new(i,j,k) - h_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_vert_remap_h_tendency, work_cont, CS%diag)
  endif

  if (show_call_tree) call callTree_leave("ALE_remap_tracers(), MOM_ALE.F90")

end procedure ALE_remap_tracers
module procedure ALE_remap_set_h_vel
  logical :: show_call_tree
  integer :: i, j, k
  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_set_h_vel()")

  ! Build the u- and v-velocity grid thicknesses for remapping.

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_u(I,j,k) = 0.5*(h_new(i,j,k) + h_new(i+1,j,k))
  endif ; enddo ; enddo ; enddo
  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_v(i,J,k) = 0.5*(h_new(i,j,k) + h_new(i,j+1,k))
  endif ; enddo ; enddo ; enddo

  ! Mask out blocked portions of velocity cells.
  if (CS%partial_cell_vel_remap) call ALE_remap_set_h_vel_partial(CS, G, GV, h_new, h_u, h_v)

  ! Take open boundary conditions into account.
  if (associated(OBC)) call ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)

  if (show_call_tree) call callTree_leave("ALE_remap_set_h_vel()")

end procedure ALE_remap_set_h_vel
module procedure ALE_remap_set_h_vel_via_dz
  logical :: show_call_tree
  integer :: i, j, k
  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_set_h_vel()")

  ! Build the u- and v-velocity grid thicknesses for remapping using the old grid and interface movement.

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_u(I,j,k) = max( 0., 0.5*(h_old(i,j,k) + h_old(i+1,j,k)) + &
            0.5 * (( dzInterface(i,j,k) + dzInterface(i+1,j,k) ) - &
                   ( dzInterface(i,j,k+1) + dzInterface(i+1,j,k+1) )) )
  endif ; enddo ; enddo ; enddo

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_v(i,J,k) = max( 0., 0.5*(h_old(i,j,k) + h_old(i,j+1,k)) + &
            0.5 * (( dzInterface(i,j,k) + dzInterface(i,j+1,k) ) - &
                   ( dzInterface(i,j,k+1) + dzInterface(i,j+1,k+1) )) )
  endif ; enddo ; enddo ; enddo

  ! Mask out blocked portions of velocity cells.
  if (CS%partial_cell_vel_remap) call ALE_remap_set_h_vel_partial(CS, G, GV, h_old, h_u, h_v)

  ! Take open boundary conditions into account.
  if (associated(OBC)) call ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)

  if (show_call_tree) call callTree_leave("ALE_remap_set_h_vel()")

end procedure ALE_remap_set_h_vel_via_dz
module procedure ALE_remap_set_h_vel_partial
  real, dimension(SZI_(G),SZJ_(G)) :: h_tot  ! The vertically summed thicknesses [H ~> m or kg m-2]
  real :: h_mask_vel ! A depth below which the thicknesses at a velocity point are masked out [H ~> m or kg m-2]
  integer :: i, j, k
  h_tot(:,:) = 0.0
  do k=1,GV%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    h_tot(i,j) = h_tot(i,j) + h_mask(i,j,k)
  enddo ; enddo ; enddo

  !$OMP parallel do default(shared) private(h_mask_vel)
  do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_mask_vel = min(h_tot(i,j), h_tot(i+1,j))
    call apply_partial_cell_mask(h_u(I,j,:), h_mask_vel)
  endif ; enddo ; enddo

  !$OMP parallel do default(shared) private(h_mask_vel)
  do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_mask_vel = min(h_tot(i,j), h_tot(i,j+1))
    call apply_partial_cell_mask(h_v(i,J,:), h_mask_vel)
  endif ; enddo ; enddo

end procedure ALE_remap_set_h_vel_partial
module procedure ALE_remap_set_h_vel_OBC
  integer :: i, j, k, nz, is_OBC, ie_OBC, js_OBC, je_OBC
  if (.not.associated(OBC)) return

  nz = GV%ke

  ! Take open boundary conditions into account.
  if (OBC%u_E_OBCs_on_PE) then
    js_OBC = max(G%jsc,  OBC%js_u_E_obc) ; je_OBC = min(G%jec,  OBC%je_u_E_obc)
    Is_OBC = max(G%IscB, OBC%Is_u_E_obc) ; Ie_OBC = min(G%IecB, OBC%Ie_u_E_obc)
    !$OMP parallel do default(shared)
    do j=js_OBC,je_OBC ; do I=Is_OBC,Ie_OBC ; if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
      do k=1,nz ; h_u(I,j,k) = h_new(i,j,k) ; enddo
    endif ; enddo ; enddo
  endif
  if (OBC%u_W_OBCs_on_PE) then
    js_OBC = max(G%jsc,  OBC%js_u_W_obc) ; je_OBC = min(G%jec,  OBC%je_u_W_obc)
    Is_OBC = max(G%IscB, OBC%Is_u_W_obc) ; Ie_OBC = min(G%IecB, OBC%Ie_u_W_obc)
    !$OMP parallel do default(shared)
    do j=js_OBC,je_OBC ; do I=Is_OBC,Ie_OBC ; if (OBC%segnum_u(I,j) < 0) then !  OBC_DIRECTION_W
      do k=1,nz ; h_u(I,j,k) = h_new(i+1,j,k) ; enddo
    endif ; enddo ; enddo
  endif

  if (OBC%v_N_OBCs_on_PE) then
    Js_OBC = max(G%JscB, OBC%Js_v_N_obc) ; Je_OBC = min(G%JecB, OBC%Je_v_N_obc)
    is_OBC = max(G%isc,  OBC%is_v_N_obc) ; ie_OBC = min(G%iec,  OBC%ie_v_N_obc)
    !$OMP parallel do default(shared)
    do J=Js_OBC,Je_OBC ; do i=is_OBC,ie_OBC ; if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
      do k=1,nz ; h_v(i,J,k) = h_new(i,j,k) ; enddo
    endif ; enddo ; enddo
  endif
  if (OBC%v_S_OBCs_on_PE) then
    Js_OBC = max(G%JscB, OBC%Js_v_S_obc) ; Je_OBC = min(G%JecB, OBC%Je_v_S_obc)
    is_OBC = max(G%isc,  OBC%is_v_S_obc) ; ie_OBC = min(G%iec,  OBC%ie_v_S_obc)
    !$OMP parallel do default(shared)
    do J=Js_OBC,Je_OBC ; do i=is_OBC,ie_OBC ; if (OBC%segnum_v(i,J) < 0) then !  OBC_DIRECTION_S
      do k=1,nz ; h_v(i,J,k) = h_new(i,j+1,k) ; enddo
    endif ; enddo ; enddo
  endif

end procedure ALE_remap_set_h_vel_OBC
module procedure ALE_remap_velocities
  real :: h_mask_vel ! A depth below which the thicknesses at a velocity point are masked out [H ~> m or kg m-2]
  real :: u_src(GV%ke)  ! A column of u-velocities on the source grid [L T-1 ~> m s-1]
  real :: u_tgt(GV%ke)  ! A column of u-velocities on the target grid [L T-1 ~> m s-1]
  real :: v_src(GV%ke)  ! A column of v-velocities on the source grid [L T-1 ~> m s-1]
  real :: v_tgt(GV%ke)  ! A column of v-velocities on the target grid [L T-1 ~> m s-1]
  real :: h1(GV%ke)     ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h2(GV%ke)     ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: rescale_coef  ! Factor that scales the baroclinic velocity to conserve ke [nondim]
  real :: u_bt, v_bt    ! Depth-averaged velocity components [L T-1 ~> m s-1]
  real :: ke_c_src, ke_c_tgt ! \int [u_c or v_c]^2 dz on src and tgt grids [H L2 T-2 ~> m3 s-2]
  real, dimension(SZIB_(G),SZJ_(G)) :: du2h_tot  ! The rate of change of vertically integrated
  real, dimension(SZI_(G),SZJB_(G)) :: dv2h_tot  ! The rate of change of vertically integrated
  real :: u2h_tot, v2h_tot   ! The vertically integrated u**2 and v**2 [H L2 T-2 ~> m3 s-2 or kg s-2]
  real :: I_dt               ! 1 / dt [T-1 ~> s-1]
  logical :: variance_option ! Contains the value of allow_preserve_variance when present, else false
  logical :: show_call_tree
  integer :: i, j, k, nz
  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_velocities()")

  ! Setup related to KE conservation
  variance_option = .false.
  if (present(allow_preserve_variance)) variance_option=allow_preserve_variance
  if (present(dt)) I_dt = 1.0 / dt

  if (CS%id_remap_delta_integ_u2>0) du2h_tot(:,:) = 0.
  if (CS%id_remap_delta_integ_v2>0) dv2h_tot(:,:) = 0.

  if (((CS%id_remap_delta_integ_u2>0) .or. (CS%id_remap_delta_integ_v2>0)) .and. .not.present(dt))&
      call MOM_error(FATAL, "ALE KE diagnostics requires passing dt into ALE_remap_velocities")

  nz = GV%ke

  ! --- Remap u profiles from the source vertical grid onto the new target grid.

  !$OMP parallel do default(shared) private(h1,h2,u_src,h_mask_vel,u_tgt, &
  !$OMP                                     u_bt,ke_c_src,ke_c_tgt,rescale_coef, &
  !$OMP                                     u2h_tot,v2h_tot)
  do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    ! Make a 1-d copy of the start and final grids and the source velocity
    do k=1,nz
      h1(k) = h_old_u(I,j,k)
      h2(k) = h_new_u(I,j,k)
      u_src(k) = u(I,j,k)
    enddo

    if (CS%id_remap_delta_integ_u2>0) then
      u2h_tot = 0.
      do k=1,nz
        u2h_tot = u2h_tot - h1(k) * (u_src(k)**2)
      enddo
    endif

    call remapping_core_h(CS%vel_remapCS, nz, h1, u_src, nz, h2, u_tgt)

    if (variance_option .and. CS%conserve_ke) then
    ! Conserve ke_u by correcting baroclinic component.
    ! Assumes total depth doesn't change during remap, and
    ! that \int u(z) dz doesn't change during remap.
      ! First get barotropic component
      u_bt = 0.0
      do k=1,nz
        u_bt = u_bt + h2(k) * u_tgt(k) ! Dimensions [H L T-1 ~> m2 s-1 or kg m-1 s-1]
      enddo
      u_bt = u_bt / (sum(h2(1:nz)) + GV%H_subroundoff) ! Dimensions return to [L T-1 ~> m s-1]
      ! Next get baroclinic ke = \int (u-u_bt)^2 from source and target
      ke_c_src = 0.0
      ke_c_tgt = 0.0
      do k=1,nz
        ke_c_src = ke_c_src + h1(k) * (u_src(k) - u_bt)**2
        ke_c_tgt = ke_c_tgt + h2(k) * (u_tgt(k) - u_bt)**2
      enddo
      ! Next rescale baroclinic component on target grid to conserve ke
      ! The values 1.5625 = 1.25**2 and 1.25 below mean that the KE-conserving
      ! correction cannot amplify the baroclinic part of velocity by more
      ! than 25%. This threshold is somewhat arbitrary. It was added to
      ! prevent unstable behavior when the amplification factor is large.
      if (ke_c_src < 1.5625 * ke_c_tgt) then
        rescale_coef = sqrt(ke_c_src / ke_c_tgt)
      else
        rescale_coef = 1.25
      endif
      do k=1,nz
        u_tgt(k) = u_bt + rescale_coef * (u_tgt(k) - u_bt)
      enddo
    endif

    if (CS%id_remap_delta_integ_u2>0) then
      do k=1,nz
        u2h_tot = u2h_tot + h2(k) * (u_tgt(k)**2)
      enddo
      du2h_tot(I,j) = u2h_tot * I_dt
    endif

    if ((CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)) &
        call mask_near_bottom_vel(u_tgt, h2, CS%BBL_h_vel_mask, CS%h_vel_mask, nz)

    ! Copy the column of new velocities back to the 3-d array
    do k=1,nz
      u(I,j,k) = u_tgt(k)
    enddo !k
  endif ; enddo ; enddo

  if (CS%id_remap_delta_integ_u2>0) call post_data(CS%id_remap_delta_integ_u2, du2h_tot, CS%diag)

  if (show_call_tree) call callTree_waypoint("u remapped (ALE_remap_velocities)")


  ! --- Remap v profiles from the source vertical grid onto the new target grid.

  !$OMP parallel do default(shared) private(h1,h2,v_src,h_mask_vel,v_tgt, &
  !$OMP                                     v_bt,ke_c_src,ke_c_tgt,rescale_coef, &
  !$OMP                                     u2h_tot,v2h_tot)
  do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then

    do k=1,nz
      h1(k) = h_old_v(i,J,k)
      h2(k) = h_new_v(i,J,k)
      v_src(k) = v(i,J,k)
    enddo

    if (CS%id_remap_delta_integ_v2>0) then
      v2h_tot = 0.
      do k=1,nz
        v2h_tot = v2h_tot - h1(k) * (v_src(k)**2)
      enddo
    endif

    call remapping_core_h(CS%vel_remapCS, nz, h1, v_src, nz, h2, v_tgt)

    if (variance_option .and. CS%conserve_ke) then
    ! Conserve ke_v by correcting baroclinic component.
    ! Assumes total depth doesn't change during remap, and
    ! that \int v(z) dz doesn't change during remap.
      ! First get barotropic component
      v_bt = 0.0
      do k=1,nz
        v_bt = v_bt + h2(k) * v_tgt(k) ! Dimensions [H L T-1 ~> m2 s-1 or kg m-1 s-1]
      enddo
      v_bt = v_bt / (sum(h2(1:nz)) + GV%H_subroundoff) ! Dimensions return to [L T-1 ~> m s-1]
      ! Next get baroclinic ke = \int (u-u_bt)^2 from source and target
      ke_c_src = 0.0
      ke_c_tgt = 0.0
      do k=1,nz
        ke_c_src = ke_c_src + h1(k) * (v_src(k) - v_bt)**2
        ke_c_tgt = ke_c_tgt + h2(k) * (v_tgt(k) - v_bt)**2
      enddo
      ! Next rescale baroclinic component on target grid to conserve ke
      if (ke_c_src < 1.5625 * ke_c_tgt) then
        rescale_coef = sqrt(ke_c_src / ke_c_tgt)
      else
        rescale_coef = 1.25
      endif
      do k=1,nz
        v_tgt(k) = v_bt + rescale_coef * (v_tgt(k) - v_bt)
      enddo
    endif

    if (CS%id_remap_delta_integ_v2>0) then
      do k=1,nz
        v2h_tot = v2h_tot + h2(k) * (v_tgt(k)**2)
      enddo
      dv2h_tot(I,j) = v2h_tot * I_dt
    endif

    if ((CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)) then
      call mask_near_bottom_vel(v_tgt, h2, CS%BBL_h_vel_mask, CS%h_vel_mask, nz)
    endif

    ! Copy the column of new velocities back to the 3-d array
    do k=1,nz
      v(i,J,k) = v_tgt(k)
    enddo !k
  endif ; enddo ; enddo

  if (CS%id_remap_delta_integ_v2>0) call post_data(CS%id_remap_delta_integ_v2, dv2h_tot, CS%diag)

  if (show_call_tree) call callTree_waypoint("v remapped (ALE_remap_velocities)")
  if (show_call_tree) call callTree_leave("ALE_remap_velocities()")

end procedure ALE_remap_velocities
module procedure ALE_remap_interface_vals
  real :: val_src(GV%ke+1)  ! A column of interface values on the source grid [A]
  real :: val_tgt(GV%ke+1)  ! A column of interface values on the target grid [A]
  real :: h_src(GV%ke)      ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h_tgt(GV%ke)      ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  integer :: i, j, k, nz
  nz = GV%ke

  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (G%mask2dT(i,j)>0.) then
    do k=1,nz
      h_src(k) = h_old(i,j,k)
      h_tgt(k) = h_new(i,j,k)
    enddo

    do K=1,nz+1
      val_src(K) = int_val(i,j,K)
    enddo

    call interpolate_column(nz, h_src, val_src, nz, h_tgt, val_tgt, .false.)

    do K=1,nz+1
      int_val(i,j,K) = val_tgt(K)
    enddo
  endif ; enddo ; enddo

end procedure ALE_remap_interface_vals
module procedure ALE_remap_vertex_vals
  real :: val_src(GV%ke+1)  ! A column of interface values on the source grid [A]
  real :: val_tgt(GV%ke+1)  ! A column of interface values on the target grid [A]
  real :: h_src(GV%ke)      ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h_tgt(GV%ke)      ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: I_mask_sum        ! The inverse of the tracer point masks surrounding a corner [nondim]
  integer :: i, j, k, nz
  nz = GV%ke

  do J=G%JscB,G%JecB ; do I=G%IscB,G%IecB
    if ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) > 0.0 ) then
      I_mask_sum = 1.0 / ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                          (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)))

    do k=1,nz
      h_src(k) = ((G%mask2dT(i,j) * h_old(i,j,k) + G%mask2dT(i+1,j+1) * h_old(i+1,j+1,k)) + &
          (G%mask2dT(i+1,j) * h_old(i+1,j,k) + G%mask2dT(i,j+1) * h_old(i,j+1,k)) ) * I_mask_sum
      h_tgt(k) = ((G%mask2dT(i,j) * h_new(i,j,k) + G%mask2dT(i+1,j+1) * h_new(i+1,j+1,k)) + &
          (G%mask2dT(i+1,j) * h_new(i+1,j,k) + G%mask2dT(i,j+1) * h_new(i,j+1,k)) ) * I_mask_sum
    enddo

    do K=1,nz+1
      val_src(K) = vert_val(I,J,K)
    enddo

    call interpolate_column(nz, h_src, val_src, nz, h_tgt, val_tgt, .false.)

    do K=1,nz+1
      vert_val(I,J,K) = val_tgt(K)
    enddo
  endif ; enddo ; enddo

end procedure ALE_remap_vertex_vals
module procedure apply_partial_cell_mask
  real :: h1_rsum  ! The running sum of h1 [H ~> m or kg m-2]
  integer :: k
  h1_rsum = 0.0
  do k=1,size(h1)
    if (h1(k) > h_mask - h1_rsum) then
      ! This thickness is reduced because it extends below the shallower neighboring bathymetry.
      h1(k) = max(h_mask - h1_rsum, 0.0)
      h1_rsum = h_mask
    else
      h1_rsum = h1_rsum + h1(k)
    endif
  enddo
end procedure apply_partial_cell_mask
module procedure mask_near_bottom_vel
  real :: h_from_bot  ! The distance between the top of a layer and the seafloor [H ~> m or kg m-2]
  integer :: k
  if ((h_BBL < 0.0) .or. (h_thin < 0.0)) return

  h_from_bot = 0.0
  do k=nk,1,-1
    h_from_bot = h_from_bot + h(k)
    if (h_from_bot > h_BBL) return
    ! Set the velocity to zero in thin, near-bottom layers.
    if (h(k) <= h_thin) vel(k) = 0.0
  enddo !k

end procedure mask_near_bottom_vel
module procedure ALE_remap_scalar
  integer :: i, j, k, n_points
  real :: dx(GV%ke+1) ! Change in interface position [H ~> m or kg m-2]
  logical :: ignore_vanished_layers, use_remapping_core_w
  ignore_vanished_layers = .false.
  if (present(all_cells)) ignore_vanished_layers = .not. all_cells
  use_remapping_core_w = .false.
  if (present(old_remap)) use_remapping_core_w = old_remap
  n_points = nk_src

  !$OMP parallel do default(shared) firstprivate(n_points,dx)
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j) > 0.) then
      if (ignore_vanished_layers) then
        n_points = 0
        do k = 1, nk_src
          if (h_src(i,j,k)>0.) n_points = n_points + 1
        enddo
        s_dst(i,j,:) = 0.
      endif
      if (use_remapping_core_w) then
        call dzFromH1H2( n_points, h_src(i,j,1:n_points), GV%ke, h_dst(i,j,:), dx )
        call remapping_core_w(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), &
                              GV%ke, dx, s_dst(i,j,:))
      else
        call remapping_core_h(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), &
                              GV%ke, h_dst(i,j,:), s_dst(i,j,:))
      endif
    else
      s_dst(i,j,:) = 0.
    endif
  enddo ; enddo

end procedure ALE_remap_scalar
module procedure TS_PLM_edge_values
  call ALE_PLM_edge_values( CS, G, GV, h, tv%S, bdry_extrap, S_t, S_b )
  call ALE_PLM_edge_values( CS, G, GV, h, tv%T, bdry_extrap, T_t, T_b )

end procedure TS_PLM_edge_values
module procedure ALE_PLM_edge_values
  integer :: i, j, k
  real :: slp(GV%ke) ! Tracer slope times the cell width [A]
  real :: mslp       ! Monotonized tracer slope times the cell width [A]
  real :: h_neglect  ! Tiny thicknesses used in remapping [H ~> m or kg m-2]
  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30
  endif

  !$OMP parallel do default(shared) private(slp,mslp)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    slp(1) = 0.
    do k = 2, GV%ke-1
      slp(k) = PLM_slope_wa(h(i,j,k-1), h(i,j,k), h(i,j,k+1), h_neglect, &
                            Q(i,j,k-1), Q(i,j,k), Q(i,j,k+1))
    enddo
    slp(GV%ke) = 0.

    do k = 2, GV%ke-1
      mslp = PLM_monotonized_slope(Q(i,j,k-1), Q(i,j,k), Q(i,j,k+1), slp(k-1), slp(k), slp(k+1))
      Q_t(i,j,k) = Q(i,j,k) - 0.5 * mslp
      Q_b(i,j,k) = Q(i,j,k) + 0.5 * mslp
    enddo
    if (bdry_extrap) then
      mslp = - PLM_extrapolate_slope(h(i,j,2), h(i,j,1), h_neglect, Q(i,j,2), Q(i,j,1))
      Q_t(i,j,1) = Q(i,j,1) - 0.5 * mslp
      Q_b(i,j,1) = Q(i,j,1) + 0.5 * mslp
      mslp = PLM_extrapolate_slope(h(i,j,GV%ke-1), h(i,j,GV%ke), h_neglect, &
                                   Q(i,j,GV%ke-1), Q(i,j,GV%ke))
      Q_t(i,j,GV%ke) = Q(i,j,GV%ke) - 0.5 * mslp
      Q_b(i,j,GV%ke) = Q(i,j,GV%ke) + 0.5 * mslp
    else
      Q_t(i,j,1) = Q(i,j,1)
      Q_b(i,j,1) = Q(i,j,1)
      Q_t(i,j,GV%ke) = Q(i,j,GV%ke)
      Q_b(i,j,GV%ke) = Q(i,j,GV%ke)
    endif

  enddo ; enddo

end procedure ALE_PLM_edge_values
module procedure TS_PPM_edge_values
  integer :: i, j, k
  real    :: hTmp(GV%ke) ! A 1-d copy of h [H ~> m or kg m-2]
  real    :: tmp(GV%ke)  ! A 1-d copy of a column of temperature [C ~> degC] or salinity [S ~> ppt]
  real, dimension(CS%nk,2) :: &
      ppol_E            ! Edge value of polynomial in [C ~> degC] or [S ~> ppt]
  real, dimension(CS%nk,3) :: &
      ppol_coefs        ! Coefficients of polynomial, all in [C ~> degC] or [S ~> ppt]
  real :: h_neglect, h_neglect_edge ! Tiny thicknesses [H ~> m or kg m-2]
  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! Determine reconstruction within each column
  !$OMP parallel do default(shared) private(hTmp,tmp,ppol_E,ppol_coefs)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1

    ! Build current grid
    hTmp(:) = h(i,j,:)
    tmp(:) = tv%S(i,j,:)

    ! Reconstruct salinity profile
    ppol_E(:,:) = 0.0
    ppol_coefs(:,:) = 0.0
    call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=h_neglect_edge, &
                                  answer_date=CS%answer_date )
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answer_date=CS%answer_date )
    if (bdry_extrap) &
        call PPM_boundary_extrapolation( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect )

    do k = 1,GV%ke
      S_t(i,j,k) = ppol_E(k,1)
      S_b(i,j,k) = ppol_E(k,2)
    enddo

    ! Reconstruct temperature profile
    ppol_E(:,:) = 0.0
    ppol_coefs(:,:) = 0.0
    tmp(:) = tv%T(i,j,:)
    if (CS%answer_date < 20190101) then
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=1.0e-10*GV%m_to_H, &
                                  answer_date=CS%answer_date )
    else
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=GV%H_subroundoff, &
                                  answer_date=CS%answer_date )
    endif
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answer_date=CS%answer_date )
    if (bdry_extrap) &
        call PPM_boundary_extrapolation(GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect )

    do k = 1,GV%ke
      T_t(i,j,k) = ppol_E(k,1)
      T_b(i,j,k) = ppol_E(k,2)
    enddo

  enddo ; enddo

end procedure TS_PPM_edge_values
module procedure TS_PLM_WLS_edge_values
  integer :: i, j, k
  type(PLM_WLS) :: recon !< A PLM-WLS reconstruction
  call recon%init(GV%ke, h_neglect=GV%H_subroundoff)

  !$OMP parallel do default(shared) firstprivate(recon)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1

    call recon%reconstruct(h(i,j,:), tv%T(i,j,:))
    T_t(i,j,:) = recon%ul(:)
    T_b(i,j,:) = recon%ur(:)

    call recon%reconstruct(h(i,j,:), tv%S(i,j,:))
    S_t(i,j,:) = recon%ul(:)
    S_b(i,j,:) = recon%ur(:)

  enddo ; enddo

  call recon%destroy()

end procedure TS_PLM_WLS_edge_values
module procedure ALE_initRegridding
  character(len=30) :: coord_mode
  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", coord_mode, &
                 "Coordinate mode for vertical regridding. "//&
                 "Choose among the following possibilities: "//&
                 trim(regriddingCoordinateModeDoc), &
                 default=DEFAULT_COORDINATE_MODE, fail_if_missing=.true.)

  call initialize_regridding(regridCS, G, GV, US, max_depth, param_file, mdl, coord_mode, '', '')

end procedure ALE_initRegridding
module procedure ALE_getCoordinate
  real, dimension(CS%nk+1) :: ALE_getCoordinate !< The coordinate positions, in the appropriate units
  ALE_getCoordinate(:) = getCoordinateInterfaces( CS%regridCS, undo_scaling=.true. )

end procedure ALE_getCoordinate
module procedure ALE_getCoordinateUnits
  ALE_getCoordinateUnits = getCoordinateUnits( CS%regridCS )

end procedure ALE_getCoordinateUnits
module procedure ALE_remap_init_conds
  ALE_remap_init_conds = .false.
  if (associated(CS)) ALE_remap_init_conds = CS%remap_after_initialization
end procedure ALE_remap_init_conds
module procedure ALE_update_regrid_weights
  real :: w  ! An implicit weighting estimate [nondim]
  if (associated(CS)) then
    w = 0.0
    if (CS%regrid_time_scale > 0.0) then
      w = CS%regrid_time_scale / (CS%regrid_time_scale + dt)
    endif
    call set_regrid_params(CS%regridCS, old_grid_weight=w)
  endif

end procedure ALE_update_regrid_weights
module procedure ALE_updateVerticalGridType
  integer :: nk
  nk = GV%ke
  GV%sInterface(1:nk+1) = getCoordinateInterfaces( CS%regridCS, undo_scaling=.true. )
  GV%sLayer(1:nk) = 0.5*( GV%sInterface(1:nk) + GV%sInterface(2:nk+1) )
  GV%zAxisUnits = getCoordinateUnits( CS%regridCS )
  GV%zAxisLongName = getCoordinateShortName( CS%regridCS )
  GV%direction = -1 ! Because of ferret in z* mode. Need method to set
                    ! as function of coordinate mode.

end procedure ALE_updateVerticalGridType
module procedure ALE_writeCoordinateFile
  character(len=240) :: filepath
  filepath = trim(directory) // trim("Vertical_coordinate.nc")

  call write_regrid_file(CS%regridCS, GV, filepath)

end procedure ALE_writeCoordinateFile
module procedure ALE_initThicknessToCoord
  real :: scale ! A scaling value for the thicknesses [nondim] or [H Z-1 ~> nondim or kg m-3]
  integer :: i, j
  scale = GV%Z_to_H
  if (present(height_units)) then ; if (height_units) scale = 1.0 ; endif
  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = scale * getStaticThickness( CS%regridCS, 0., max(G%meanSL(i,j)+G%bathyT(i,j), 0.0) )
  enddo ; enddo

end procedure ALE_initThicknessToCoord
end submodule MOM_ALE_s
