submodule (MOM_dyn_horgrid) MOM_dyn_horgrid_s
  implicit none
contains
module procedure create_dyn_horgrid
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg
  if (associated(G)) then
    call MOM_error(WARNING, "create_dyn_horgrid called with an associated horgrid_type.")
  else
    allocate(G)
  endif

  G%HI = HI

  G%isc = HI%isc ; G%iec = HI%iec ; G%jsc = HI%jsc ; G%jec = HI%jec
  G%isd = HI%isd ; G%ied = HI%ied ; G%jsd = HI%jsd ; G%jed = HI%jed
  G%isg = HI%isg ; G%ieg = HI%ieg ; G%jsg = HI%jsg ; G%jeg = HI%jeg

  G%IscB = HI%IscB ; G%IecB = HI%IecB ; G%JscB = HI%JscB ; G%JecB = HI%JecB
  G%IsdB = HI%IsdB ; G%IedB = HI%IedB ; G%JsdB = HI%JsdB ; G%JedB = HI%JedB
  G%IsgB = HI%IsgB ; G%IegB = HI%IegB ; G%JsgB = HI%JsgB ; G%JegB = HI%JegB

  G%idg_offset = HI%idg_offset ; G%jdg_offset = HI%jdg_offset
  G%isd_global = G%isd + HI%idg_offset ; G%jsd_global = G%jsd + HI%jdg_offset
  G%symmetric = HI%symmetric

  G%bathymetry_at_vel = .false.
  if (present(bathymetry_at_vel)) G%bathymetry_at_vel = bathymetry_at_vel

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  allocate(G%dxT(isd:ied,jsd:jed), source=0.0)
  allocate(G%dxCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%dxCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dxBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IdxT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IdxCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IdxCu_OBCmask(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IdxCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IdxBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%dyT(isd:ied,jsd:jed), source=0.0)
  allocate(G%dyCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%dyCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dyBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IdyT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IdyCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IdyCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IdyCv_OBCmask(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IdyBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%areaT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IareaT(isd:ied,jsd:jed), source=0.0)
  allocate(G%areaBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IareaBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%mask2dT(isd:ied,jsd:jed), source=0.0)
  allocate(G%mask2dCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%mask2dCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%mask2dBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%OBCmaskCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%OBCmaskCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLatT(isd:ied,jsd:jed), source=0.0)
  allocate(G%geoLatCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%geoLatCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLatBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%geoLonT(isd:ied,jsd:jed), source=0.0)
  allocate(G%geoLonCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%geoLonCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLonBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%dx_Cv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dy_Cu(IsdB:IedB,jsd:jed), source=0.0)

  allocate(G%areaCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%areaCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IareaCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IareaCv(isd:ied,JsdB:JedB), source=0.0)

  allocate(G%porous_DminU(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%porous_DmaxU(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%porous_DavgU(IsdB:IedB,jsd:jed), source=0.0)

  allocate(G%porous_DminV(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%porous_DmaxV(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%porous_DavgV(isd:ied,JsdB:JedB), source=0.0)

  allocate(G%bathyT(isd:ied, jsd:jed), source=0.0)
  allocate(G%meanSL(isd:ied, jsd:jed), source=0.0)
  allocate(G%CoriolisBu(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(G%Coriolis2Bu(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(G%dF_dx(isd:ied, jsd:jed), source=0.0)
  allocate(G%dF_dy(isd:ied, jsd:jed), source=0.0)

  allocate(G%sin_rot(isd:ied,jsd:jed), source=0.0)
  allocate(G%cos_rot(isd:ied,jsd:jed), source=1.0)

  if (G%bathymetry_at_vel) then
    allocate(G%Dblock_u(IsdB:IedB, jsd:jed), source=0.0)
    allocate(G%Dopen_u(IsdB:IedB, jsd:jed), source=0.0)
    allocate(G%Dblock_v(isd:ied, JsdB:JedB), source=0.0)
    allocate(G%Dopen_v(isd:ied, JsdB:JedB), source=0.0)
  endif

  ! gridLonB and gridLatB are used as edge values in some cases, so they
  ! always need to use symmetric memory allcoations.
  allocate(G%gridLonT(isg:ieg), source=0.0)
  allocate(G%gridLonB(isg-1:ieg), source=0.0)
  allocate(G%gridLatT(jsg:jeg), source=0.0)
  allocate(G%gridLatB(jsg-1:jeg), source=0.0)

end procedure create_dyn_horgrid
module procedure rotate_dyn_horgrid
  call rotate_array(G_in%geoLonT, turns, G%geoLonT)
  call rotate_array(G_in%geoLatT, turns, G%geoLatT)
  call rotate_array_pair(G_in%dxT, G_in%dyT, turns, G%dxT, G%dyT)
  call rotate_array(G_in%areaT, turns, G%areaT)
  call rotate_array(G_in%bathyT, turns, G%bathyT)
  call rotate_array(G_in%meanSL, turns, G%meanSL)

  call rotate_array_pair(G_in%df_dx, G_in%df_dy, turns, G%df_dx, G%df_dy)
  call rotate_array(G_in%sin_rot, turns, G%sin_rot)
  call rotate_array(G_in%cos_rot, turns, G%cos_rot)
  call rotate_array(G_in%mask2dT, turns, G%mask2dT)

  ! Face points
  call rotate_array_pair(G_in%geoLonCu, G_in%geoLonCv, turns, G%geoLonCu, G%geoLonCv)
  call rotate_array_pair(G_in%geoLatCu, G_in%geoLatCv, turns, G%geoLatCu, G%geoLatCv)
  call rotate_array_pair(G_in%dxCu, G_in%dyCv, turns, G%dxCu, G%dyCv)
  call rotate_array_pair(G_in%dxCv, G_in%dyCu, turns, G%dxCv, G%dyCu)
  call rotate_array_pair(G_in%dx_Cv, G_in%dy_Cu, turns, G%dx_Cv, G%dy_Cu)

  call rotate_array_pair(G_in%mask2dCu, G_in%mask2dCv, turns, G%mask2dCu, G%mask2dCv)
  call rotate_array_pair(G_in%OBCmaskCu, G_in%OBCmaskCv, turns, G%OBCmaskCu, G%OBCmaskCv)
  call rotate_array_pair(G_in%areaCu, G_in%areaCv, turns, G%areaCu, G%areaCv)
  call rotate_array_pair(G_in%IareaCu, G_in%IareaCv, turns, G%IareaCu, G%IareaCv)

  call rotate_array_pair(G_in%porous_DminU, G_in%porous_DminV, &
       turns, G%porous_DminU, G%porous_DminV)
  call rotate_array_pair(G_in%porous_DmaxU, G_in%porous_DmaxV, &
       turns, G%porous_DmaxU, G%porous_DmaxV)
  call rotate_array_pair(G_in%porous_DavgU, G_in%porous_DavgV, &
       turns, G%porous_DavgU, G%porous_DavgV)


  ! Vertex point
  call rotate_array(G_in%geoLonBu, turns, G%geoLonBu)
  call rotate_array(G_in%geoLatBu, turns, G%geoLatBu)
  call rotate_array_pair(G_in%dxBu, G_in%dyBu, turns, G%dxBu, G%dyBu)
  call rotate_array(G_in%areaBu, turns, G%areaBu)
  call rotate_array(G_in%CoriolisBu, turns, G%CoriolisBu)
  call rotate_array(G_in%Coriolis2Bu, turns, G%Coriolis2Bu)
  call rotate_array(G_in%mask2dBu, turns, G%mask2dBu)

  ! Topography at the cell faces
  G%bathymetry_at_vel = G_in%bathymetry_at_vel
  if (G%bathymetry_at_vel) then
    call rotate_array_pair(G_in%Dblock_u, G_in%Dblock_v, turns, G%Dblock_u, G%Dblock_v)
    call rotate_array_pair(G_in%Dopen_u, G_in%Dopen_v, turns, G%Dopen_u, G%Dopen_v)
  endif

  ! Nominal grid axes
  ! TODO: We should not assign lat values to the lon axis, and vice versa.
  !   We temporarily copy lat <-> lon since several components still expect
  !   lat and lon sizes to match the first and second dimension sizes.
  !   But we ought to instead leave them unchanged and adjust the references to
  !   these axes.
  if (modulo(turns, 2) /= 0) then
    G%gridLonT(:) = G_in%gridLatT(G_in%jeg:G_in%jsg:-1)
    G%gridLatT(:) = G_in%gridLonT(:)
    G%gridLonB(:) = G_in%gridLatB(G_in%jeg:(G_in%jsg-1):-1)
    G%gridLatB(:) = G_in%gridLonB(:)
  else
    G%gridLonT(:) = G_in%gridLonT(:)
    G%gridLatT(:) = G_in%gridLatT(:)
    G%gridLonB(:) = G_in%gridLonB(:)
    G%gridLatB(:) = G_in%gridLatB(:)
  endif

  G%x_axis_units = G_in%y_axis_units
  G%y_axis_units = G_in%x_axis_units
  G%x_ax_unit_short = G_in%y_ax_unit_short
  G%y_ax_unit_short = G_in%x_ax_unit_short
  G%south_lat = G_in%south_lat
  G%west_lon = G_in%west_lon
  G%len_lat = G_in%len_lat
  G%len_lon = G_in%len_lon

  ! Rotation-invariant fields
  G%grid_unit_to_L = G_in%grid_unit_to_L
  G%areaT_global = G_in%areaT_global
  G%IareaT_global = G_in%IareaT_global
  G%Rad_Earth_L = G_in%Rad_Earth_L
  G%max_depth = G_in%max_depth

  call set_derived_dyn_horgrid(G, US)
end procedure rotate_dyn_horgrid
module procedure rescale_dyn_horgrid_bathymetry
  real :: rescale ! The inverse of m_in_new_units, used in rescaling bathymetry [Z m-1 ~> 1]
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (m_in_new_units == 1.0) return
  if (m_in_new_units < 0.0) &
    call MOM_error(FATAL, "rescale_dyn_horgrid_bathymetry: Negative depth units are not permitted.")
  if (m_in_new_units == 0.0) &
    call MOM_error(FATAL, "rescale_dyn_horgrid_bathymetry: Zero depth units are not permitted.")

  rescale = 1.0 / m_in_new_units
  do j=jsd,jed ; do i=isd,ied
    G%bathyT(i,j) = rescale*G%bathyT(i,j)
    G%meanSL(i,j) = rescale*G%meanSL(i,j)
  enddo ; enddo
  if (G%bathymetry_at_vel) then ; do j=jsd,jed ; do I=IsdB,IedB
    G%Dblock_u(I,j) = rescale*G%Dblock_u(I,j) ; G%Dopen_u(I,j) = rescale*G%Dopen_u(I,j)
  enddo ; enddo ; endif
  if (G%bathymetry_at_vel) then ; do J=JsdB,JedB ; do i=isd,ied
    G%Dblock_v(i,J) = rescale*G%Dblock_v(i,J) ; G%Dopen_v(i,J) = rescale*G%Dopen_v(i,J)
  enddo ; enddo ; endif
  G%max_depth = rescale*G%max_depth

end procedure rescale_dyn_horgrid_bathymetry
module procedure set_derived_dyn_horgrid
  integer :: i, j, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do j=jsd,jed ; do i=isd,ied
    if (G%dxT(i,j) < 0.0) G%dxT(i,j) = 0.0
    if (G%dyT(i,j) < 0.0) G%dyT(i,j) = 0.0
    G%IdxT(i,j) = Adcroft_reciprocal(G%dxT(i,j))
    G%IdyT(i,j) = Adcroft_reciprocal(G%dyT(i,j))
    G%IareaT(i,j) = Adcroft_reciprocal(G%areaT(i,j))
  enddo ; enddo

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dxCu(I,j) < 0.0) G%dxCu(I,j) = 0.0
    if (G%dyCu(I,j) < 0.0) G%dyCu(I,j) = 0.0
    G%IdxCu(I,j) = Adcroft_reciprocal(G%dxCu(I,j))
    G%IdyCu(I,j) = Adcroft_reciprocal(G%dyCu(I,j))
    G%IdxCu_OBCmask(I,j) = G%OBCmaskCu(I,j) * G%IdxCu(I,j) ! This may be reset when the masks are set.
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dxCv(i,J) < 0.0) G%dxCv(i,J) = 0.0
    if (G%dyCv(i,J) < 0.0) G%dyCv(i,J) = 0.0
    G%IdxCv(i,J) = Adcroft_reciprocal(G%dxCv(i,J))
    G%IdyCv(i,J) = Adcroft_reciprocal(G%dyCv(i,J))
    G%IdyCv_OBCmask(i,J) = G%OBCmaskCv(i,J) * G%IdyCv(i,J) ! This may be reset when the masks are set.
  enddo ; enddo

  do J=JsdB,JedB ; do I=IsdB,IedB
    if (G%dxBu(I,J) < 0.0) G%dxBu(I,J) = 0.0
    if (G%dyBu(I,J) < 0.0) G%dyBu(I,J) = 0.0

    G%IdxBu(I,J) = Adcroft_reciprocal(G%dxBu(I,J))
    G%IdyBu(I,J) = Adcroft_reciprocal(G%dyBu(I,J))
    ! areaBu has usually been set to a positive area elsewhere.
    if (G%areaBu(I,J) <= 0.0) G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    G%IareaBu(I,J) =  Adcroft_reciprocal(G%areaBu(I,J))
  enddo ; enddo

end procedure set_derived_dyn_horgrid
module procedure Adcroft_reciprocal
  I_val = 0.0 ; if (val /= 0.0) I_val = 1.0/val
end procedure Adcroft_reciprocal
module procedure destroy_dyn_horgrid
  if (.not.associated(G)) then
    call MOM_error(FATAL, "destroy_dyn_horgrid called with an unassociated horgrid_type.")
  endif

  deallocate(G%dxT)  ; deallocate(G%dxCu)  ; deallocate(G%dxCv)  ; deallocate(G%dxBu)
  deallocate(G%IdxT) ; deallocate(G%IdxCu) ; deallocate(G%IdxCv) ; deallocate(G%IdxBu)

  deallocate(G%dyT)  ; deallocate(G%dyCu)  ; deallocate(G%dyCv)  ; deallocate(G%dyBu)
  deallocate(G%IdyT) ; deallocate(G%IdyCu) ; deallocate(G%IdyCv) ; deallocate(G%IdyBu)

  deallocate(G%areaT)  ; deallocate(G%IareaT)
  deallocate(G%areaBu) ; deallocate(G%IareaBu)
  deallocate(G%areaCu) ; deallocate(G%IareaCu)
  deallocate(G%areaCv) ; deallocate(G%IareaCv)

  deallocate(G%mask2dT)  ; deallocate(G%mask2dCu) ; deallocate(G%OBCmaskCu)
  deallocate(G%mask2dCv) ; deallocate(G%OBCmaskCv) ; deallocate(G%mask2dBu)
  deallocate(G%IdxCu_OBCmask) ; deallocate(G%IdyCv_OBCmask)

  deallocate(G%geoLatT)  ; deallocate(G%geoLatCu)
  deallocate(G%geoLatCv) ; deallocate(G%geoLatBu)
  deallocate(G%geoLonT)  ; deallocate(G%geoLonCu)
  deallocate(G%geoLonCv) ; deallocate(G%geoLonBu)

  deallocate(G%dx_Cv) ; deallocate(G%dy_Cu)

  deallocate(G%porous_DminU) ; deallocate(G%porous_DmaxU) ; deallocate(G%porous_DavgU)
  deallocate(G%porous_DminV) ; deallocate(G%porous_DmaxV) ; deallocate(G%porous_DavgV)

  deallocate(G%bathyT)     ; deallocate(G%meanSL)
  deallocate(G%CoriolisBu) ; deallocate(G%Coriolis2Bu)
  deallocate(G%dF_dx)      ; deallocate(G%dF_dy)
  deallocate(G%sin_rot)    ; deallocate(G%cos_rot)

  if (allocated(G%Dblock_u)) deallocate(G%Dblock_u)
  if (allocated(G%Dopen_u)) deallocate(G%Dopen_u)
  if (allocated(G%Dblock_v)) deallocate(G%Dblock_v)
  if (allocated(G%Dopen_v)) deallocate(G%Dopen_v)

  deallocate(G%gridLonT) ; deallocate(G%gridLatT)
  deallocate(G%gridLonB) ; deallocate(G%gridLatB)

  ! CS%debug is required to validate Domain_aux, so use allocation test
  if (associated(G%Domain_aux)) call deallocate_MOM_domain(G%Domain_aux)

  call deallocate_MOM_domain(G%Domain)

  deallocate(G)

end procedure destroy_dyn_horgrid
end submodule MOM_dyn_horgrid_s
