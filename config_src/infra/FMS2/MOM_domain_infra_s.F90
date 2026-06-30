submodule (MOM_domain_infra) MOM_domain_infra_s
  implicit none
contains
module procedure pass_var_3d
  integer :: dirflag
  logical :: block_til_complete
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position, &
                        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                          complete=block_til_complete, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_3d
module procedure pass_var_2d
  real, allocatable, dimension(:,:) :: tmp
  integer :: pos, i_halo, j_halo
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  integer :: inner, i, j, isfw, iefw, isfe, iefe, jsfs, jefs, jsfn, jefn
  integer :: dirflag
  logical :: block_til_complete
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true. ; if (present(complete)) block_til_complete = complete
  pos = CENTER ; if (present(position)) pos = position

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    ! Store the original values.
    allocate(tmp(size(array,1), size(array,2)))
    tmp(:,:) = array(:,:)
    block_til_complete = .true.
  endif ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position, &
                        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position)
  endif

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    call mpp_get_compute_domain(MOM_dom%mpp_domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(MOM_dom%mpp_domain, isd, ied, jsd, jed)
    ! Convert to local indices for arrays starting at 1.
    isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
    jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
    i_halo = min(inner_halo, isc-1) ; j_halo = min(inner_halo, jsc-1)

    ! Figure out the array index extents of the eastern, western, northern and southern regions to copy.
    if (pos == CENTER) then
      if (size(array,1) == ied) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for CENTER array.") ; endif
      if (size(array,2) == jed) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for CENTER array.") ; endif
    elseif (pos == CORNER) then
      if (size(array,1) == ied) then
        isfw = max(isc - (i_halo+1), 1) ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      elseif (size(array,1) == ied+1) then
        isfw = isc - i_halo ; iefw = isc+1 ; isfe = iec+1 ; iefe = min(iec + 1 + i_halo, ied+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for CORNER array.") ; endif
      if (size(array,2) == jed) then
        jsfs = max(jsc - (j_halo+1), 1) ; jefs = jsc ; jsfn = jec ; jefn = jec + j_halo
      elseif (size(array,2) == jed+1) then
        jsfs = jsc - j_halo ; jefs = jsc+1 ; jsfn = jec+1 ; jefn = min(jec + 1 + j_halo, jed+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for CORNER array.") ; endif
    elseif (pos == NORTH_FACE) then
      if (size(array,1) == ied) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for NORTH_FACE array.") ; endif
      if (size(array,2) == jed) then
        jsfs = max(jsc - (j_halo+1), 1) ; jefs = jsc ; jsfn = jec ; jefn = jec + j_halo
      elseif (size(array,2) == jed+1) then
        jsfs = jsc - j_halo ; jefs = jsc+1 ; jsfn = jec+1 ; jefn = min(jec + 1 + j_halo, jed+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for NORTH_FACE array.") ; endif
    elseif (pos == EAST_FACE) then
      if (size(array,1) == ied) then
        isfw = max(isc - (i_halo+1), 1) ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      elseif (size(array,1) == ied+1) then
        isfw = isc - i_halo ; iefw = isc+1 ; isfe = iec+1 ; iefe = min(iec + 1 + i_halo, ied+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for EAST_FACE array.") ; endif
      if (size(array,2) == jed) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for EAST_FACE array.") ; endif
    else
      call MOM_error(FATAL, "pass_var_2d: Unrecognized position")
    endif

    ! Copy back the stored inner halo points
    do j=jsfs,jefn ; do i=isfw,iefw ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefn ; do i=isfe,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefs ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfn,jefn ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo

    deallocate(tmp)
  endif ; endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_2d
module procedure pass_var_start_2d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_var_start_2d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_var_start_2d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_start_2d
module procedure pass_var_start_3d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_var_start_3d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_var_start_3d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_start_3d
module procedure pass_var_complete_2d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                     flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_complete_2d
module procedure pass_var_complete_3d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                     flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_var_complete_3d
module procedure pass_vector_2d
  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_2d
module procedure fill_vector_symmetric_edges_2d
  integer :: stagger_local
  integer :: dirflag
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  real, allocatable, dimension(:) :: sbuff_x, sbuff_y, wbuff_x, wbuff_y
  logical :: block_til_complete
  if (.not. MOM_dom%symmetric) then
      return
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  if (.not.(stagger_local == CGRID_NE .or. stagger_local == BGRID_NE)) return

  call mpp_get_compute_domain(MOM_dom%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(MOM_dom%mpp_domain, isd, ied, jsd, jed)

  ! Adjust isc, etc., to account for the fact that the input arrays indices all
  ! start at 1 (and are effectively on a SW grid!).
  isc = isc - (isd-1) ; iec = iec - (isd-1)
  jsc = jsc - (jsd-1) ; jec = jec - (jsd-1)
  IscB = isc ; IecB = iec+1 ; JscB = jsc ; JecB = jec+1

  dirflag = To_All ! 60
  if (present(scalar)) then ; if (scalar) dirflag = To_All+SCALAR_PAIR ; endif

  if (stagger_local == CGRID_NE) then
    allocate(wbuff_x(jsc:jec)) ; allocate(sbuff_y(isc:iec))
    wbuff_x(:) = 0.0 ; sbuff_y(:) = 0.0
    call mpp_get_boundary(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          wbufferx=wbuff_x, sbuffery=sbuff_y, &
                          gridtype=CGRID_NE)
    do i=isc,iec
      v_cmpt(i,JscB) = sbuff_y(i)
    enddo
    do j=jsc,jec
      u_cmpt(IscB,j) = wbuff_x(j)
    enddo
    deallocate(wbuff_x) ; deallocate(sbuff_y)
  elseif  (stagger_local == BGRID_NE) then
    allocate(wbuff_x(JscB:JecB)) ; allocate(sbuff_x(IscB:IecB))
    allocate(wbuff_y(JscB:JecB)) ; allocate(sbuff_y(IscB:IecB))
    wbuff_x(:) = 0.0 ; wbuff_y(:) = 0.0 ; sbuff_x(:) = 0.0 ; sbuff_y(:) = 0.0
    call mpp_get_boundary(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          wbufferx=wbuff_x, sbufferx=sbuff_x, &
                          wbuffery=wbuff_y, sbuffery=sbuff_y, &
                          gridtype=BGRID_NE)
    do I=IscB,IecB
      u_cmpt(I,JscB) = sbuff_x(I) ; v_cmpt(I,JscB) = sbuff_y(I)
    enddo
    do J=JscB,JecB
      u_cmpt(IscB,J) = wbuff_x(J) ; v_cmpt(IscB,J) = wbuff_y(J)
    enddo
    deallocate(wbuff_x) ; deallocate(sbuff_x)
    deallocate(wbuff_y) ; deallocate(sbuff_y)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure fill_vector_symmetric_edges_2d
module procedure pass_vector_3d
  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_3d
module procedure pass_vector_start_2d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_vector_start_2d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_vector_start_2d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_start_2d
module procedure pass_vector_start_3d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_vector_start_3d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_vector_start_3d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_start_3d
module procedure pass_vector_complete_2d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
             whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_complete_2d
module procedure pass_vector_complete_3d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure pass_vector_complete_3d
module procedure create_var_group_pass_2d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position, whalo=halo, ehalo=halo, &
                                 shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure create_var_group_pass_2d
module procedure create_var_group_pass_3d
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position, whalo=halo, ehalo=halo, &
                                 shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure create_var_group_pass_3d
module procedure create_vector_group_pass_2d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local, whalo=halo, ehalo=halo, &
            shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure create_vector_group_pass_2d
module procedure create_vector_group_pass_3d
  integer :: stagger_local
  integer :: dirflag
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local, whalo=halo, ehalo=halo, &
            shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure create_vector_group_pass_3d
module procedure do_group_pass
  real :: d_type
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_do_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure do_group_pass
module procedure start_group_pass
  real                                 :: d_type
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_start_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure start_group_pass
module procedure complete_group_pass
  real                                 :: d_type
  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_complete_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end procedure complete_group_pass
module procedure redistribute_array_2d
  logical :: do_complete
  do_complete=.true. ; if (PRESENT(complete)) do_complete = complete

  call mpp_redistribute(Domain1, array1, Domain2, array2, do_complete)

end procedure redistribute_array_2d
module procedure redistribute_array_3d
  logical :: do_complete
  do_complete=.true. ; if (PRESENT(complete)) do_complete = complete

  call mpp_redistribute(Domain1, array1, Domain2, array2, do_complete)

end procedure redistribute_array_3d
module procedure redistribute_array_4d
  logical :: do_complete
  do_complete=.true. ; if (PRESENT(complete)) do_complete = complete

  call mpp_redistribute(Domain1, array1, Domain2, array2, do_complete)

end procedure redistribute_array_4d
module procedure rescale_comp_data_4d
  logical :: unsign_zeros ! If true, convert negative zeros into ordinary signless zeros.
  integer :: is, ie, js, je, i, j, k, m
  unsign_zeros = .false. ; if (present(zero_zeros)) unsign_zeros = zero_zeros

  if ((scale == 1.0) .and. (.not.unsign_zeros)) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  if (scale /= 1.0) &
    array(is:ie,js:je,:,:) = scale*array(is:ie,js:je,:,:)

  if (unsign_zeros) then ! Convert negative zeros into zeros
    do m=1,size(array,4) ; do k=1,size(array,3) ; do j=js,je ; do i=is,ie
      if (array(i,j,k,m) == 0.0) array(i,j,k,m) = 0.0
    enddo ; enddo ; enddo ; enddo
  endif

end procedure rescale_comp_data_4d
module procedure rescale_comp_data_3d
  logical :: unsign_zeros ! If true, convert negative zeros into ordinary signless zeros.
  integer :: is, ie, js, je, i, j, k
  unsign_zeros = .false. ; if (present(zero_zeros)) unsign_zeros = zero_zeros

  if ((scale == 1.0) .and. (.not.unsign_zeros)) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  if (scale /= 1.0) &
    array(is:ie,js:je,:) = scale*array(is:ie,js:je,:)

  if (unsign_zeros) then ! Convert negative zeros into zeros
    do k=1,size(array,3) ; do j=js,je ; do i=is,ie
      if (array(i,j,k) == 0.0) array(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

end procedure rescale_comp_data_3d
module procedure rescale_comp_data_2d
  logical :: unsign_zeros ! If true, convert negative zeros into ordinary signless zeros.
  integer :: is, ie, js, je, i, j
  unsign_zeros = .false. ; if (present(zero_zeros)) unsign_zeros = zero_zeros

  if ((scale == 1.0) .and. (.not.unsign_zeros)) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  if (scale /= 1.0) &
    array(is:ie,js:je) = scale*array(is:ie,js:je)

  if (unsign_zeros) then ! Convert negative zeros into zeros
    do j=js,je ; do i=is,ie
      if (array(i,j) == 0.0) array(i,j) = 0.0
    enddo ; enddo
  endif

end procedure rescale_comp_data_2d
module procedure create_MOM_domain
  integer, dimension(4) :: global_indices ! The lower and upper global i- and j-index bounds
  integer :: X_FLAGS  ! A combination of integers encoding the x-direction grid connectivity.
  integer :: Y_FLAGS  ! A combination of integers encoding the y-direction grid connectivity.
  character(len=200) :: mesg    ! A string for use in error messages
  logical :: mask_table_exists  ! Mask_table is present and the file it points to exists
  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
  endif

  MOM_dom%name = "MOM" ; if (present(domain_name)) MOM_dom%name = trim(domain_name)

  X_FLAGS = 0 ; Y_FLAGS = 0
  if (reentrant(1)) X_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (reentrant(2)) Y_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (tripolar_N) then
    Y_FLAGS = FOLD_NORTH_EDGE
    if (reentrant(2)) call MOM_error(FATAL,"MOM_domains: "// &
      "TRIPOLAR_N and REENTRANT_Y may not be used together.")
  endif

  MOM_dom%nonblocking_updates = .false.
  if (present(nonblocking)) MOM_dom%nonblocking_updates = nonblocking
  MOM_dom%thin_halo_updates = .false.
  if (present(thin_halos)) MOM_dom%thin_halo_updates = thin_halos
  MOM_dom%symmetric = .true. ; if (present(symmetric)) MOM_dom%symmetric = symmetric
  MOM_dom%niglobal = n_global(1) ; MOM_dom%njglobal = n_global(2)
  MOM_dom%nihalo = n_halo(1) ; MOM_dom%njhalo = n_halo(2)

  ! Save the extra data for creating other domains of different resolution that overlay this domain.
  MOM_dom%X_FLAGS = X_FLAGS
  MOM_dom%Y_FLAGS = Y_FLAGS
  MOM_dom%layout(:) = layout(:)

  ! Set up the io_layout, with error handling.
  MOM_dom%io_layout(:) = (/ 1, 1 /)
  if (present(io_layout)) then
    if (io_layout(1) == 0) then
      MOM_dom%io_layout(1) = layout(1)
    elseif (io_layout(1) > 1) then
      MOM_dom%io_layout(1) = io_layout(1)
      if (modulo(layout(1), io_layout(1)) /= 0) then
        write(mesg,'("MOM_domains_init: The i-direction I/O-layout, IO_LAYOUT(1)=",i4, &
              &", does not evenly divide the i-direction layout, NIPROC=,",i4,".")') io_layout(1), layout(1)
        call MOM_error(FATAL, mesg)
      endif
    endif

    if (io_layout(2) == 0) then
      MOM_dom%io_layout(2) = layout(2)
    elseif (io_layout(2) > 1) then
      MOM_dom%io_layout(2) = io_layout(2)
      if (modulo(layout(2), io_layout(2)) /= 0) then
        write(mesg,'("MOM_domains_init: The j-direction I/O-layout, IO_LAYOUT(2)=",i4, &
              &", does not evenly divide the j-direction layout, NJPROC=,",i4,".")') io_layout(2), layout(2)
        call MOM_error(FATAL, mesg)
      endif
    endif
  endif

  if (present(mask_table)) then
    mask_table_exists = file_exists(mask_table)
    if (mask_table_exists) then
      allocate(MOM_dom%maskmap(layout(1), layout(2)))
      call parse_mask_table(mask_table, MOM_dom%maskmap, MOM_dom%name)
    endif
  else
    mask_table_exists = .false.
  endif

  ! Initialize as an unrotated domain
  MOM_dom%turns = 0

  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain)

end procedure create_MOM_domain
module procedure deallocate_MOM_domain
  logical :: invasive  ! If true, deallocate fields associated with the underlying infrastructure
  integer :: n
  invasive = .true. ; if (present(cursory)) invasive = .not.cursory

  if (associated(MOM_domain)) then
    if (associated(MOM_domain%mpp_domain)) then
      if (invasive) call mpp_deallocate_domain(MOM_domain%mpp_domain)
      deallocate(MOM_domain%mpp_domain)
    endif
    if (associated(MOM_domain%mpp_domain_d)) then
      if (invasive) then ; do n=1,size(MOM_domain%mpp_domain_d)
       call mpp_deallocate_domain(MOM_domain%mpp_domain_d(n))
      enddo ; endif
      deallocate(MOM_domain%mpp_domain_d)
    endif
    if (associated(MOM_domain%maskmap)) deallocate(MOM_domain%maskmap)
    deallocate(MOM_domain)
  endif

end procedure deallocate_MOM_domain
module procedure MOM_thread_affinity_set
  MOM_thread_affinity_set = .false.
  !$ call fms_affinity_init()
  !$OMP PARALLEL
  !$OMP   MASTER
  !$        ocean_nthreads = omp_get_num_threads()
  !$OMP   END MASTER
  !$OMP END PARALLEL
  !$ MOM_thread_affinity_set = (ocean_nthreads > 1 )
end procedure MOM_thread_affinity_set
module procedure set_MOM_thread_affinity
end procedure set_MOM_thread_affinity
module procedure get_domain_components_MD
  call mpp_get_domain_components(MOM_dom%mpp_domain, x_domain, y_domain)
end procedure get_domain_components_MD
module procedure get_domain_components_d2D
  call mpp_get_domain_components(domain, x_domain, y_domain)
end procedure get_domain_components_d2D
module procedure clone_MD_to_MD
  integer :: global_indices(4)
  logical :: mask_table_exists
  integer, dimension(:), allocatable :: exni ! The extents of the grid for each i-row of the layout.
  integer, dimension(:), allocatable :: exnj ! The extents of the grid for each j-row of the layout.
  integer :: qturns ! The number of quarter turns, restricted to the range of 0 to 3.
  integer :: i, j, nl1, nl2
  integer :: io_layout_in(2)
  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (present(io_layout)) then
    io_layout_in(:) = io_layout(:)
  else
    io_layout_in(:) = MD_in%io_layout(:)
  endif

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
  endif

! Save the extra data for creating other domains of different resolution that overlay this domain
  MOM_dom%symmetric = MD_in%symmetric
  MOM_dom%nonblocking_updates = MD_in%nonblocking_updates
  MOM_dom%thin_halo_updates = MD_in%thin_halo_updates

  if (modulo(qturns, 2) /= 0) then
    MOM_dom%niglobal = MD_in%njglobal ; MOM_dom%njglobal = MD_in%niglobal
    MOM_dom%nihalo = MD_in%njhalo ; MOM_dom%njhalo = MD_in%nihalo
    call get_layout_extents(MD_in, exnj, exni)

    MOM_dom%X_FLAGS = MD_in%Y_FLAGS ; MOM_dom%Y_FLAGS = MD_in%X_FLAGS
    ! Correct the position of a tripolar grid, assuming that flags are not additive.
    if (modulo(qturns, 4) == 1) then
      if (MD_in%Y_FLAGS == FOLD_NORTH_EDGE) MOM_dom%X_FLAGS = FOLD_EAST_EDGE
      if (MD_in%Y_FLAGS == FOLD_SOUTH_EDGE) MOM_dom%X_FLAGS = FOLD_WEST_EDGE
      if (MD_in%X_FLAGS == FOLD_EAST_EDGE) MOM_dom%Y_FLAGS = FOLD_SOUTH_EDGE
      if (MD_in%X_FLAGS == FOLD_WEST_EDGE) MOM_dom%Y_FLAGS = FOLD_NORTH_EDGE
    elseif (modulo(qturns, 4) == 3) then
      if (MD_in%Y_FLAGS == FOLD_NORTH_EDGE) MOM_dom%X_FLAGS = FOLD_WEST_EDGE
      if (MD_in%Y_FLAGS == FOLD_SOUTH_EDGE) MOM_dom%X_FLAGS = FOLD_EAST_EDGE
      if (MD_in%X_FLAGS == FOLD_EAST_EDGE) MOM_dom%Y_FLAGS = FOLD_NORTH_EDGE
      if (MD_in%X_FLAGS == FOLD_WEST_EDGE) MOM_dom%Y_FLAGS = FOLD_SOUTH_EDGE
    endif

    MOM_dom%layout(:) = MD_in%layout(2:1:-1)
    MOM_dom%io_layout(:) = io_layout_in(2:1:-1)
  else
    MOM_dom%niglobal = MD_in%niglobal ; MOM_dom%njglobal = MD_in%njglobal
    MOM_dom%nihalo = MD_in%nihalo ; MOM_dom%njhalo = MD_in%njhalo
    call get_layout_extents(MD_in, exni, exnj)

    MOM_dom%X_FLAGS = MD_in%X_FLAGS ; MOM_dom%Y_FLAGS = MD_in%Y_FLAGS
    ! Correct the position of a tripolar grid, assuming that flags are not additive.
    if (modulo(qturns, 4) == 2) then
      if (MD_in%Y_FLAGS == FOLD_NORTH_EDGE) MOM_dom%Y_FLAGS = FOLD_SOUTH_EDGE
      if (MD_in%Y_FLAGS == FOLD_SOUTH_EDGE) MOM_dom%Y_FLAGS = FOLD_NORTH_EDGE
      if (MD_in%X_FLAGS == FOLD_EAST_EDGE) MOM_dom%X_FLAGS = FOLD_WEST_EDGE
      if (MD_in%X_FLAGS == FOLD_WEST_EDGE) MOM_dom%X_FLAGS = FOLD_EAST_EDGE
    endif

    MOM_dom%layout(:) = MD_in%layout(:)
    MOM_dom%io_layout(:) = io_layout_in(:)
  endif

  ! Ensure that the points per processor are the same on the source and destination grids.
  select case (qturns)
    case (1) ; call invert(exni)
    case (2) ; call invert(exni) ; call invert(exnj)
    case (3) ; call invert(exnj)
  end select

  if (associated(MD_in%maskmap)) then
    mask_table_exists = .true.
    allocate(MOM_dom%maskmap(MOM_dom%layout(1), MOM_dom%layout(2)))

    nl1 = MOM_dom%layout(1) ; nl2 = MOM_dom%layout(2)
    select case (qturns)
      case (0)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(i, j)
        enddo ; enddo
      case (1)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(j, nl1+1-i)
        enddo ; enddo
      case (2)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(nl1+1-i, nl2+1-j)
        enddo ; enddo
      case (3)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(nl2+1-j, i)
        enddo ; enddo
    end select
  else
    mask_table_exists = .false.
  endif

  ! Optionally enhance the grid resolution.
  if (present(refine)) then ; if (refine > 1) then
    MOM_dom%niglobal = refine*MOM_dom%niglobal ; MOM_dom%njglobal = refine*MOM_dom%njglobal
    MOM_dom%nihalo = refine*MOM_dom%nihalo ; MOM_dom%njhalo = refine*MOM_dom%njhalo
    do i=1,MOM_dom%layout(1) ; exni(i) = refine*exni(i) ; enddo
    do j=1,MOM_dom%layout(2) ; exnj(j) = refine*exnj(j) ; enddo
  endif ; endif

  ! Optionally enhance the grid resolution.
  if (present(extra_halo)) then ; if (extra_halo > 0) then
    MOM_dom%nihalo = MOM_dom%nihalo + extra_halo ; MOM_dom%njhalo = MOM_dom%njhalo + extra_halo
  endif ; endif

  if (present(halo_size) .and. present(min_halo)) call MOM_error(FATAL, &
      "clone_MOM_domain can not have both halo_size and min_halo present.")

  if (present(min_halo)) then
    MOM_dom%nihalo = max(MOM_dom%nihalo, min_halo(1))
    min_halo(1) = MOM_dom%nihalo
    MOM_dom%njhalo = max(MOM_dom%njhalo, min_halo(2))
    min_halo(2) = MOM_dom%njhalo
  endif

  if (present(halo_size)) then
    MOM_dom%nihalo = halo_size ; MOM_dom%njhalo = halo_size
  endif

  if (present(symmetric)) then ; MOM_dom%symmetric = symmetric ; endif

  if (present(domain_name)) then
    MOM_dom%name = trim(domain_name)
  else
    MOM_dom%name = MD_in%name
  endif

  MOM_dom%turns = qturns
  if (qturns /= 0) then
    MOM_dom%domain_in => MD_in
  endif

  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain, xextent=exni, yextent=exnj)

end procedure clone_MD_to_MD
module procedure clone_MD_to_d2D
  integer :: global_indices(4)
  integer :: nihalo, njhalo
  logical :: symmetric_dom, do_coarsen
  character(len=64) :: dom_name
  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for MOM_domain to domain2d")

  if (present(halo_size) .and. present(min_halo)) call MOM_error(FATAL, &
      "clone_MOM_domain can not have both halo_size and min_halo present.")

  do_coarsen = .false. ; if (present(coarsen)) then ; do_coarsen = (coarsen > 1) ; endif

  nihalo = MD_in%nihalo ; njhalo = MD_in%njhalo
  if (do_coarsen) then
    nihalo = int(MD_in%nihalo / coarsen) ; njhalo = int(MD_in%njhalo / coarsen)
  endif

  if (present(min_halo)) then
    nihalo = max(nihalo, min_halo(1))
    njhalo = max(njhalo, min_halo(2))
    min_halo(1) = nihalo ; min_halo(2) = njhalo
  endif
  if (present(halo_size)) then
    nihalo = halo_size ; njhalo = halo_size
  endif

  symmetric_dom = MD_in%symmetric
  if (present(symmetric)) then ; symmetric_dom = symmetric ; endif

  dom_name = MD_in%name
  if (do_coarsen) dom_name = trim(MD_in%name)//"c"
  if (present(domain_name)) dom_name = trim(domain_name)

  global_indices(1:4) = (/ 1, MD_in%niglobal, 1, MD_in%njglobal /)
  if (do_coarsen) then
    global_indices(1:4) = (/ 1, (MD_in%niglobal/coarsen), 1, (MD_in%njglobal/coarsen) /)
  endif

  if (associated(MD_in%maskmap)) then
    call mpp_define_domains( global_indices, MD_in%layout, mpp_domain, &
                xflags=MD_in%X_FLAGS, yflags=MD_in%Y_FLAGS, xhalo=nihalo, yhalo=njhalo, &
                xextent=xextent, yextent=yextent, symmetry=symmetric_dom, name=dom_name, &
                maskmap=MD_in%maskmap )
  else
    call mpp_define_domains( global_indices, MD_in%layout, mpp_domain, &
                xflags=MD_in%X_FLAGS, yflags=MD_in%Y_FLAGS, xhalo=nihalo, yhalo=njhalo, &
                symmetry=symmetric_dom, xextent=xextent, yextent=yextent, name=dom_name)
  endif

  if ((MD_in%io_layout(1) + MD_in%io_layout(2) > 0) .and. &
      (MD_in%layout(1)*MD_in%layout(2) > 1)) then
    call mpp_define_io_domain(mpp_domain, MD_in%io_layout)
  else
    call mpp_define_io_domain(mpp_domain, (/ 1, 1 /) )
  endif

end procedure clone_MD_to_d2D
module procedure get_domain_extent_MD
  integer :: isg_, ieg_, jsg_, jeg_
  integer :: ind_off, idg_off, jdg_off, coarsen_lev
  logical :: local
  local = .true. ; if (present(local_indexing)) local = local_indexing
  ind_off = 0 ; if (present(index_offset)) ind_off = index_offset

  coarsen_lev = 0 ; if (present(coarsen)) coarsen_lev = coarsen

  if (coarsen_lev == 0) then
    call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain%mpp_domain, isg_, ieg_, jsg_, jeg_)
  else
    if (.not.associated(Domain%mpp_domain_d)) call MOM_error(FATAL, &
            "get_domain_extent called with coarsen_lev, but Domain%mpp_domain_d(coarsen_lev) is not associated.")
    call mpp_get_compute_domain(Domain%mpp_domain_d(coarsen_lev), isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain%mpp_domain_d(coarsen_lev), isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain%mpp_domain_d(coarsen_lev), isg_, ieg_, jsg_, jeg_)
  endif

  if (local) then
    ! This code institutes the MOM convention that local array indices start at 1.
    idg_off = isd - 1 ; jdg_off = jsd - 1
    isc = isc - isd + 1 ; iec = iec - isd + 1 ; jsc = jsc - jsd + 1 ; jec = jec - jsd + 1
    ied = ied - isd + 1 ; jed = jed - jsd + 1
    isd = 1 ; jsd = 1
  else
    idg_off = 0 ; jdg_off = 0
  endif
  if (ind_off /= 0) then
    idg_off = idg_off + ind_off ; jdg_off = jdg_off + ind_off
    isc = isc + ind_off ; iec = iec + ind_off
    jsc = jsc + ind_off ; jec = jec + ind_off
    isd = isd + ind_off ; ied = ied + ind_off
    jsd = jsd + ind_off ; jed = jed + ind_off
  endif
  if (present(isg)) isg = isg_
  if (present(ieg)) ieg = ieg_
  if (present(jsg)) jsg = jsg_
  if (present(jeg)) jeg = jeg_
  if (present(idg_offset)) idg_offset = idg_off
  if (present(jdg_offset)) jdg_offset = jdg_off
  if (present(symmetric)) symmetric = Domain%symmetric

end procedure get_domain_extent_MD
module procedure get_domain_extent_d2D
  integer :: isd_, ied_, jsd_, jed_, jsg_, jeg_, isg_, ieg_
  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain, isd_, ied_, jsd_, jed_)

  if (present(isd)) isd = isd_
  if (present(ied)) ied = ied_
  if (present(jsd)) jsd = jsd_
  if (present(jed)) jed = jed_

end procedure get_domain_extent_d2D
module procedure get_simple_array_i_ind
  logical :: sym
  character(len=120) :: mesg, mesg2
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)

  isc = isc-isd+1 ; iec = iec-isd+1 ; ied = ied-isd+1 ; isd = 1
  sym = Domain%symmetric ; if (present(symmetric)) sym = symmetric

  if (size == ied) then ; is = isc ; ie = iec
  elseif (size == 1+iec-isc) then ; is = 1 ; ie = size
  elseif (sym .and. (size == 1+ied)) then ; is = isc ; ie = iec+1
  elseif (sym .and. (size == 2+iec-isc)) then ; is = 1 ; ie = size+1
  else
    write(mesg,'("Unrecognized size ", i6, "in call to get_simple_array_i_ind.  \")') size
    if (sym) then
      write(mesg2,'("Valid sizes are : ", 2i7)') ied, 1+iec-isc
    else
      write(mesg2,'("Valid sizes are : ", 4i7)') ied, 1+iec-isc, 1+ied, 2+iec-isc
    endif
    call MOM_error(FATAL, trim(mesg)//trim(mesg2))
  endif

end procedure get_simple_array_i_ind
module procedure get_simple_array_j_ind
  logical :: sym
  character(len=120) :: mesg, mesg2
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)

  jsc = jsc-jsd+1 ; jec = jec-jsd+1 ; jed = jed-jsd+1 ; jsd = 1
  sym = Domain%symmetric ; if (present(symmetric)) sym = symmetric

  if (size == jed) then ; js = jsc ; je = jec
  elseif (size == 1+jec-jsc) then ; js = 1 ; je = size
  elseif (sym .and. (size == 1+jed)) then ; js = jsc ; je = jec+1
  elseif (sym .and. (size == 2+jec-jsc)) then ; js = 1 ; je = size+1
  else
    write(mesg,'("Unrecognized size ", i6, "in call to get_simple_array_j_ind.  \")') size
    if (sym) then
      write(mesg2,'("Valid sizes are : ", 2i7)') jed, 1+jec-jsc
    else
      write(mesg2,'("Valid sizes are : ", 4i7)') jed, 1+jec-jsc, 1+jed, 2+jec-jsc
    endif
    call MOM_error(FATAL, trim(mesg)//trim(mesg2))
  endif

end procedure get_simple_array_j_ind
module procedure invert
  integer :: i, ni, swap
  ni = size(array)
  do i=1,ni
    swap = array(i)
    array(i) = array(ni+1-i)
    array(ni+1-i) = swap
  enddo
end procedure invert
module procedure get_global_shape
  niglobal = domain%niglobal
  njglobal = domain%njglobal
end procedure get_global_shape
module procedure compute_block_extent
  call mpp_compute_block_extent(isg, ieg, ndivs, ibegin, iend)
end procedure compute_block_extent
module procedure compute_extent
  call mpp_compute_extent(isg, ieg, ndivs, ibegin, iend)
end procedure compute_extent
module procedure broadcast_domain
  call mpp_broadcast_domain(domain)
end procedure broadcast_domain
module procedure global_field
  call mpp_global_field(domain, local, global)
end procedure global_field
module procedure same_domain
  integer :: isc_a, iec_a, jsc_a, jec_a, isc_b, iec_b, jsc_b, jec_b
  integer :: layout_a(2), layout_b(2)
  call mpp_get_layout(domain_a, layout_a)
  call mpp_get_layout(domain_b, layout_b)

  call get_domain_extent(domain_a, isc_a, iec_a, jsc_a, jec_a)
  call get_domain_extent(domain_b, isc_b, iec_b, jsc_b, jec_b)

  same_domain = (layout_a(1) == layout_b(1)) .and. (layout_a(2) == layout_b(2)) .and. &
                (iec_a - isc_a == iec_b - isc_b) .and. (jec_a - jsc_a == jec_b - jsc_b)

end procedure same_domain
module procedure get_layout_extents
  if (allocated(extent_i)) deallocate(extent_i)
  if (allocated(extent_j)) deallocate(extent_j)
  allocate(extent_i(domain%layout(1))) ; extent_i(:) = 0
  allocate(extent_j(domain%layout(2))) ; extent_j(:) = 0
  call mpp_get_domain_extents(domain%mpp_domain, extent_i, extent_j)
end procedure get_layout_extents
module procedure set_domain
end procedure set_domain
module procedure nullify_domain
end procedure nullify_domain
end submodule MOM_domain_infra_s
