submodule (MOM_debugging) MOM_debugging_s
  implicit none
contains
module procedure MOM_debugging_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_debugging" ! This module's name.
  call log_version(param_file, mdl, version, debugging=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_CHKSUMS", debug_chksums, &
                 "If true, checksums are performed on arrays in the "//&
                 "various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", debug_redundant, &
                 "If true, debug redundant data points during calls to "//&
                 "the various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)

  call MOM_checksums_init(param_file)

end procedure MOM_debugging_init
module procedure query_debugging_checks
  if (present(do_debug)) do_debug = debug
  if (present(do_chksums)) do_chksums = debug_chksums
  if (present(do_redundant)) do_redundant = debug_redundant

end procedure query_debugging_checks
module procedure check_redundant_vC3d
  character(len=24) :: mesg_k
  integer :: k
  do k=1,size(u_comp,3)
    write(mesg_k,'(" Layer ",i0," ")') k
    call check_redundant_vC2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction, unscale)
  enddo
end procedure check_redundant_vC3d
module procedure check_redundant_vC2d
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)  ! A nonsymmetric version of u_comp [A ~> a]
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)  ! A nonsymmetric version of v_comp [A ~> a]
  real :: u_resym(G%IsdB:G%IedB,G%jsd:G%jed) ! A reconstructed symmetric version of u_comp [A ~> a]
  real :: v_resym(G%isd:G%ied,G%JsdB:G%JedB) ! A reconstructed symmetric version of v_comp [A ~> a]
  real :: sc  ! A factor that undoes the scaling for the arrays to give consistent output [a A-1 ~> 1]
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  sc  = 1.0 ; if (present(unscale)) sc = unscale

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL, &
    " check_redundant called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction)

  do I=IsdB,IedB ; do j=jsd,jed ; u_resym(I,j) = u_comp(I,j) ; enddo ; enddo
  do i=isd,ied ; do J=JsdB,JedB ; v_resym(i,J) = v_comp(i,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," on pe ",I0)') &
           sc*u_comp(i,j), sc*u_resym(i,j), sc*(u_comp(i,j)-u_resym(i,j)), i, j, pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," x,y = ",2(1pe12.4)," on pe ",I0)') &
           sc*v_comp(i,j), sc*v_resym(i,j), sc*(v_comp(i,j)-v_resym(i,j)), i, j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo

end procedure check_redundant_vC2d
module procedure check_redundant_sB3d
  character(len=24) :: mesg_k
  integer :: k
  do k=1,size(array,3)
    write(mesg_k,'(" Layer ",i0," ")') k
    call check_redundant_sB2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je, unscale)
  enddo
end procedure check_redundant_sB3d
module procedure check_redundant_sB2d
  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)    ! A nonsymmetric version of array [A ~> a]
  real :: a_resym(G%IsdB:G%IedB,G%JsdB:G%JedB) ! A reconstructed symmetric version of array [A ~> a]
  real :: sc  ! A factor that undoes the scaling for the arrays to give consistent output [a A-1 ~> 1]
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  sc = 1.0 ; if (present(unscale)) sc = unscale

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL, &
    " check_redundant called with a non-associated auxiliary domain the grid type.")
  call pass_vector(a_nonsym, a_nonsym, G%Domain_aux, &
                   direction=To_All+Scalar_Pair, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB ; a_resym(I,J) = array(I,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    a_resym(i,j) = a_nonsym(i,j)
  enddo ; enddo
  call pass_vector(a_resym, a_resym, G%Domain, direction=To_All+Scalar_Pair, &
                   stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_resym(i,j) /= array(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," on pe ",I0)') &
           sc*array(i,j), sc*a_resym(i,j), sc*(array(i,j)-a_resym(i,j)), i, j, pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end procedure check_redundant_sB2d
module procedure check_redundant_vB3d
  character(len=24) :: mesg_k
  integer :: k
  do k=1,size(u_comp,3)
    write(mesg_k,'(" Layer ",i0," ")') k
    call check_redundant_vB2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction, unscale)
  enddo
end procedure check_redundant_vB3d
module procedure check_redundant_vB2d
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)    ! A nonsymmetric version of u_comp [A ~> a]
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)    ! A nonsymmetric version of v_comp [A ~> a]
  real :: u_resym(G%IsdB:G%IedB,G%JsdB:G%JedB) ! A reconstructed symmetric version of u_comp [A ~> a]
  real :: v_resym(G%IsdB:G%IedB,G%JsdB:G%JedB) ! A reconstructed symmetric version of v_comp [A ~> a]
  real :: sc  ! A factor that undoes the scaling for the arrays to give consistent output [a A-1 ~> 1]
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  sc = 1.0 ; if (present(unscale)) sc = unscale

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL, &
    " check_redundant called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB
    u_resym(I,J) = u_comp(I,J) ; v_resym(I,J) = v_comp(I,J)
  enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction, stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," on pe ",I0)') &
           sc*u_comp(i,j), sc*u_resym(i,j), sc*(u_comp(i,j)-u_resym(i,j)), i, j, pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo
  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," x,y = ",2(1pe12.4)," on pe ",I0)') &
           sc*v_comp(i,j), sc*v_resym(i,j), sc*(v_comp(i,j)-v_resym(i,j)), i, j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end procedure check_redundant_vB2d
module procedure check_redundant_sT3d
  character(len=24) :: mesg_k
  integer :: k
  do k=1,size(array,3)
    write(mesg_k,'(" Layer ",i0," ")') k
    call check_redundant_sT2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je, unscale)
  enddo
end procedure check_redundant_sT3d
module procedure check_redundant_sT2d
  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)  ! A version of array with halo points updated by message passing [A ~> a]
  real :: sc ! A factor that undoes the scaling for the arrays to give consistent output [a A-1 ~> 1]
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  sc = 1.0 ; if (present(unscale)) sc = unscale

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  call pass_var(a_nonsym, G%Domain)

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_nonsym(i,j) /= array(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," on pe ",I0)') &
           sc*array(i,j), sc*a_nonsym(i,j), sc*(array(i,j)-a_nonsym(i,j)), i, j, pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo

end procedure check_redundant_sT2d
module procedure check_redundant_vT3d
  character(len=24) :: mesg_k
  integer :: k
  do k=1,size(u_comp,3)
    write(mesg_k,'(" Layer ",i0," ")') k
    call check_redundant_vT2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction, unscale)
  enddo
end procedure check_redundant_vT3d
module procedure check_redundant_vT2d
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed) ! A version of u_comp with halo points updated by message passing [A ~> a]
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed) ! A version of v_comp with halo points updated by message passing [A ~> a]
  real :: sc ! A factor that undoes the scaling for the arrays to give consistent output [a A-1 ~> 1]
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  sc = 1.0 ; if (present(unscale)) sc = unscale

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  call pass_vector(u_nonsym, v_nonsym, G%Domain, direction, stagger=AGRID)

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_nonsym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," on pe ",I0)') &
           sc*u_comp(i,j), sc*u_nonsym(i,j), sc*(u_comp(i,j)-u_nonsym(i,j)), i, j, pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_nonsym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",I0,",",I0," x,y = ",2(1pe12.4)," on pe ",I0)') &
           sc*v_comp(i,j), sc*v_nonsym(i,j), sc*(v_comp(i,j)-v_nonsym(i,j)), i, j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo

end procedure check_redundant_vT2d
module procedure chksum_vec_C3d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_C3d
module procedure chksum_vec_C2d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_C2d
module procedure chksum_vec_B3d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos, unscale=unscale)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_B3d
module procedure chksum_vec_B2d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos, symmetric=symmetric, unscale=unscale)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos, symmetric=symmetric, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_B2d
module procedure chksum_vec_A3d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos, unscale=unscale)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_A3d
module procedure chksum_vec_A2d
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos, unscale=unscale)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos, unscale=unscale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair, unscale=unscale)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G, unscale=unscale)
    endif
  endif

end procedure chksum_vec_A2d
module procedure totalStuff
  real  :: tmp_for_sum(HI%isc:HI%iec, HI%jsc:HI%jec)  ! The column integrated amount of stuff in a
  integer :: i, j, k, nz
  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * stuff(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  totalStuff = reproducing_sum(tmp_for_sum, unscale=unscale)

end procedure totalStuff
module procedure totalTandS
  real, save :: totalH = 0.   ! The total ocean volume or mass, saved for the next
  real, save :: totalT = 0.   ! The total volume integrated ocean temperature, saved for the next
  real, save :: totalS = 0.   ! The total volume integrated ocean salinity, saved for the next
  logical, save :: firstCall = .true.
  real :: tmp_for_sum(HI%isc:HI%iec, HI%jsc:HI%jec) ! The volume of each column [H L2 ~> m3 or kg] or [m3] or [kg]
  real :: thisH, delH  ! The total ocean volume and the change from the last call [H L2 ~> m3 or kg] or [m3] or [kg]
  real :: thisT, delT  ! The current total volume integrated temperature and the change from the last
  real :: thisS, delS  ! The current total volume integrated salinity and the change from the last
  real :: H_unscale    ! A constant that translates thickness units to its MKS units (m or kg m-2) based on
  real :: HL2_unscale  ! An overall unscaling factor for cell mass or volume [m3 H-1 L-2 ~> 1] or [kg H-1 L-2 ~> 1]
  real :: T_unscale    ! An overall unscaling factor for cell-integrated temperature [degC m3 C-1 H-1 L-2 ~> 1] or
  real :: S_unscale    ! An overall unscaling factor for cell-integrated salinity [ppt m3 S-1 H-1 L-2 ~> 1] or
  integer :: i, j, k, nz
  H_unscale = 1.0 ; if (present(H_to_mks)) H_unscale = H_to_mks
  if (present(US)) then
    HL2_unscale = US%L_to_m**2 * H_unscale
    T_unscale = US%C_to_degC * HL2_unscale ; S_unscale = US%S_to_ppt * HL2_unscale
  else
    HL2_unscale = H_unscale
    T_unscale = HL2_unscale ; S_unscale = HL2_unscale
  endif

  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  thisH = reproducing_sum(tmp_for_sum, unscale=HL2_unscale)
  thisT = totalStuff(HI, hThick, areaT, temperature, unscale=T_unscale)
  thisS = totalStuff(HI, hThick, areaT, salinity, unscale=S_unscale)

  if (is_root_pe()) then
    if (firstCall) then
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(stdout,*) 'Totals H,T,S:', thisH*HL2_unscale, thisT*T_unscale, thisS*S_unscale, ' ', mesg
      firstCall = .false.
    else
      delH = thisH - totalH
      delT = thisT - totalT
      delS = thisS - totalS
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Tot/del H,T,S:', thisH*HL2_unscale, thisT*T_unscale, thisS*S_unscale, &
                                   delH*HL2_unscale,  delT*T_unscale,  delS*S_unscale, ' ', mesg
    endif
  endif

end procedure totalTandS
module procedure check_column_integral
  real    :: u_sum    ! The vertical sum of the field [arbitrary]
  real    :: error    ! An estimate of the roundoff error in the sum [arbitrary]
  real    :: expected ! The expected vertical sum [arbitrary]
  integer :: k
  u_sum = field(1)
  error = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk
    u_sum = u_sum + field(k)
    error = error + EPSILON(u_sum)*MAX(ABS(u_sum),ABS(field(k)))
  enddo

  ! Assign expected answer to either the optional input or 0
  if (present(known_answer)) then
    expected = known_answer
  else
    expected = 0.
  endif

  ! Compare the column integrals against calculated roundoff error
  if (abs(u_sum-expected) > error) then
    check_column_integral = .true.
  else
    check_column_integral = .false.
  endif

end procedure check_column_integral
module procedure check_column_integrals
  real    :: u1_sum, u2_sum ! The vertical sums of the two fields [arbitrary]
  real    :: error1, error2 ! Estimates of the roundoff errors in the sums [arbitrary]
  real    :: misval         ! The missing value flag, indicating elements that are to be omitted
  integer :: k
  if (present(missing_value)) then
    misval = missing_value
  else
    misval = 0.
  endif

  u1_sum = field_1(1)
  error1 = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk_1
    if (field_1(k) /= misval) then
      u1_sum = u1_sum + field_1(k)
      error1 = error1 + EPSILON(u1_sum)*MAX(ABS(u1_sum),ABS(field_1(k)))
    endif
  enddo

  u2_sum = field_2(1)
  error2 = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk_2
    if (field_2(k) /= misval) then
      u2_sum = u2_sum + field_2(k)
      error2 = error2 + EPSILON(u2_sum)*MAX(ABS(u2_sum),ABS(field_2(k)))
    endif
  enddo

  ! Compare the column integrals against calculated roundoff error
  if (abs(u1_sum-u2_sum) > (error1+error2)) then
    check_column_integrals = .true.
  else
    check_column_integrals = .false.
  endif

end procedure check_column_integrals
end submodule MOM_debugging_s
