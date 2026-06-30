submodule (MOM_PointAccel) MOM_PointAccel_s
#include <MOM_memory.h>
  implicit none
contains
module procedure write_u_accel
  real    :: CFL              ! The local velocity-based CFL number [nondim]
  real    :: Angstrom         ! A negligibly small thickness [H ~> m or kg m-2]
  real    :: du               ! A velocity change [L T-1 ~> m s-1]
  real    :: Inorm(SZK_(GV))  ! The inverse of the normalized velocity change [L T-1 ~> m s-1]
  real    :: e(SZK_(GV)+1)    ! Simple estimates of interface heights based on the sum of thicknesses [m]
  real    :: h_scale          ! A scaling factor for thicknesses [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  real    :: vel_scale        ! A scaling factor for velocities [m T s-1 L-1 ~> 1]
  real    :: uh_scale         ! A scaling factor for transport per unit length [m2 T s-1 L-1 H-1 ~> 1]
  real    :: temp_scale       ! A scaling factor for temperatures [degC C-1 ~> 1]
  real    :: saln_scale       ! A scaling factor for salinities [ppt S-1 ~> 1]
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(GV)+1)
  logical :: prev_avail
  integer :: file
  Angstrom = GV%Angstrom_H + GV%H_subroundoff
  h_scale = GV%H_to_mks ; vel_scale = US%L_T_to_m_s ; uh_scale = h_scale*vel_scale
  temp_scale = US%C_to_degC ; saln_scale = US%S_to_ppt

!  if (.not.associated(CS)) return
  nz = GV%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%u_file == -1) then
      if (len_trim(CS%u_trunc_file) < 1) return
      call open_ASCII_file(CS%u_file, trim(CS%u_trunc_file), action=APPEND_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%u_file == -1) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%u_trunc_file)//'.')
        return
      endif
    endif
    file = CS%u_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

  ! Determine which layers to write out accelerations for.
    do k=1,nz
      if (((max(CS%u_av(I,j,k),um(I,j,k)) >= vel_rpt) .or. &
           (min(CS%u_av(I,j,k),um(I,j,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%u_av(I,j,k), um(I,j,k)) >= vel_rpt) .or. &
           (min(CS%u_av(I,j,k), um(I,j,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1 ; ke = nz ; write(file,'("U: Unable to set ks & ke.")')
    endif
    if (CS%full_column) then
      ks = 1 ; ke = nz
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",I0," ",I0," ",F6.2," U-velocity violation at ",I0,": ",I0,", ",I0, &
        & " (",F7.2," E ",F7.2," N) Layers ",I0," to ",I0,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), I, j, &
        G%geoLonCu(I,j), G%geoLatCu(I,j), ks, ke, US%T_to_s*dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ")', advance='no') (k) ; enddo
    write(file,'(/,"u(m):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*um(I,j,k)) ; enddo
    if (prev_avail) then
      write(file,'(/,"u(mp): ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%u_prev(I,j,k)) ; enddo
    endif
    write(file,'(/,"u(3):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%u_av(I,j,k)) ; enddo

    write(file,'(/,"CFL u: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(um(I,j,k)) * dt * G%dy_Cu(I,j)
      if (um(I,j,k) < 0.0) then ; CFL = CFL * G%IareaT(i+1,j)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ")', advance='no') CFL
    endif ; enddo
    write(file,'(/,"CFL0 u:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    abs(um(I,j,k)) * dt * G%IdxCu(I,j) ; enddo

    if (prev_avail) then
      write(file,'(/,"du:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*(um(I,j,k)-CS%u_prev(I,j,k))) ; enddo
    endif
    write(file,'(/,"CAu:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%CAu(I,j,k)) ; enddo
    write(file,'(/,"PFu:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%PFu(I,j,k)) ; enddo
    write(file,'(/,"diffu: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%diffu(I,j,k)) ; enddo

    if (associated(ADp%gradKEu)) then
      write(file,'(/,"KEu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%gradKEu(I,j,k)) ; enddo
    endif
    if (associated(ADp%rv_x_v)) then
      write(file,'(/,"Coru:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k)) ; enddo
    endif
    if (associated(ADp%du_dt_visc)) then
      write(file,'(/,"ubv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*(um(I,j,k) - dt*ADp%du_dt_visc(I,j,k)) ; enddo
      write(file,'(/,"duv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%du_dt_visc(I,j,k)) ; enddo
    endif
    if (associated(ADp%du_other)) then
      write(file,'(/,"du_other: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*ADp%du_other(I,j,k)) ; enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",ES10.3," ")', advance='no') h_scale*a(I,j,ks)*dt
      do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') (h_scale*a(I,j,K)*dt) ; enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hv(I,j,k) ; enddo
    endif
    if (present(str)) then
      write(file,'(/,"Stress:  ",ES10.3)', advance='no') (uh_scale*GV%RZ_to_H) * (str*dt)
    endif

    if (associated(CS%u_accel_bt)) then
      write(file,'(/,"dubt:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*CS%u_accel_bt(I,j,k)) ; enddo
    endif
    write(file,'(/)')

    write(file,'(/,"h--:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j-1,k)) ; enddo
    write(file,'(/,"h+-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j-1,k)) ; enddo
    write(file,'(/,"h-0:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j,k)) ; enddo
    write(file,'(/,"h+0:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j,k)) ; enddo
    write(file,'(/,"h-+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j+1,k)) ; enddo
    write(file,'(/,"h++:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j+1,k)) ; enddo


    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"e-:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i+1,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i+1,j,k) ; enddo
    write(file,'(/,"e+:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo
    if (associated(CS%T)) then
      write(file,'(/,"T-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j,k) ; enddo
      write(file,'(/,"T+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i+1,j,k) ; enddo
    endif
    if (associated(CS%S)) then
      write(file,'(/,"S-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j,k) ; enddo
      write(file,'(/,"S+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i+1,j,k) ; enddo
    endif

    if (prev_avail) then
      write(file,'(/,"v--:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J-1,k)) ; enddo
      write(file,'(/,"v-+:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J,k)) ; enddo
      write(file,'(/,"v+-:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i+1,J-1,k)) ; enddo
      write(file,'(/,"v++:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i+1,J,k)) ; enddo
    endif

    write(file,'(/,"vh--:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%vh(i,J-1,k)*G%IdxCv(i,J-1)) ; enddo
    write(file,'(/," vhC--:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                        (0.5*CS%v_av(i,j-1,k)*uh_scale*(hin(i,j-1,k) + hin(i,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp--:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (0.5*CS%v_prev(i,j-1,k)*uh_scale*(hin(i,j-1,k) + hin(i,j,k))) ; enddo
    endif

    write(file,'(/,"vh-+:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%vh(i,J,k)*G%IdxCv(i,J)) ; enddo
    write(file,'(/," vhC-+:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                        (0.5*CS%v_av(i,J,k)*uh_scale*(hin(i,j,k) + hin(i,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp-+:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (0.5*CS%v_prev(i,J,k)*uh_scale*(hin(i,j,k) + hin(i,j+1,k))) ; enddo
    endif

    write(file,'(/,"vh+-:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (uh_scale*CDp%vh(i+1,J-1,k)*G%IdxCv(i+1,J-1)) ; enddo
    write(file,'(/," vhC+-:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                    (0.5*CS%v_av(i+1,J-1,k)*uh_scale*(hin(i+1,j-1,k) + hin(i+1,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp+-:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                      (0.5*CS%v_prev(i+1,J-1,k)*uh_scale*(hin(i+1,j-1,k) + hin(i+1,j,k))) ; enddo
    endif

    write(file,'(/,"vh++:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (uh_scale*CDp%vh(i+1,J,k)*G%IdxCv(i+1,J)) ; enddo
    write(file,'(/," vhC++:")', advance='no')
         do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                     (0.5*CS%v_av(i+1,J,k)*uh_scale*(hin(i+1,j,k) + hin(i+1,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp++:")', advance='no')
           do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                       (0.5*CS%v_av(i+1,J,k)*uh_scale*(hin(i+1,j,k) + hin(i+1,j+1,k))) ; enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') US%Z_to_m*(G%bathyT(i,j) + G%Z_ref), US%Z_to_m*(G%bathyT(i+1,j) + G%Z_ref)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        du = um(I,j,k) - CS%u_prev(I,j,k)
        if (abs(du) < 1.0e-6*US%m_s_to_L_T) du = 1.0e-6*US%m_s_to_L_T
        Inorm(k) = 1.0 / du
      enddo

      write(file,'(2/,"Norm:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') (vel_scale / Inorm(k)) ; enddo

      write(file,'(/,"du:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      ((um(I,j,k)-CS%u_prev(I,j,k)) * Inorm(k)) ; enddo

      write(file,'(/,"CAu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%CAu(I,j,k) * Inorm(k)) ; enddo

      write(file,'(/,"PFu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%PFu(I,j,k) * Inorm(k)) ; enddo

      write(file,'(/,"diffu: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%diffu(I,j,k) * Inorm(k)) ; enddo

      if (associated(ADp%gradKEu)) then
        write(file,'(/,"KEu:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%gradKEu(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%rv_x_v)) then
        write(file,'(/,"Coru:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k)) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%du_dt_visc)) then
        write(file,'(/,"duv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%du_dt_visc(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%du_other)) then
        write(file,'(/,"du_other: ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (ADp%du_other(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(CS%u_accel_bt)) then
        write(file,'(/,"dubt:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*CS%u_accel_bt(I,j,k) * Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    flush(file)
  endif

end procedure write_u_accel
module procedure write_v_accel
  real    :: CFL              ! The local velocity-based CFL number [nondim]
  real    :: Angstrom         ! A negligibly small thickness [H ~> m or kg m-2]
  real    :: dv               ! A velocity change [L T-1 ~> m s-1]
  real    :: Inorm(SZK_(GV))  ! The inverse of the normalized velocity change [L T-1 ~> m s-1]
  real    :: e(SZK_(GV)+1)    ! Simple estimates of interface heights based on the sum of thicknesses [m]
  real    :: h_scale          ! A scaling factor for thicknesses [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  real    :: vel_scale        ! A scaling factor for velocities [m T s-1 L-1 ~> 1]
  real    :: uh_scale         ! A scaling factor for transport per unit length [m2 T s-1 L-1 H-1 ~> 1]
  real    :: temp_scale       ! A scaling factor for temperatures [degC C-1 ~> 1]
  real    :: saln_scale       ! A scaling factor for salinities [ppt S-1 ~> 1]
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(GV)+1)
  logical :: prev_avail
  integer :: file
  Angstrom = GV%Angstrom_H + GV%H_subroundoff
  h_scale = GV%H_to_mks ; vel_scale = US%L_T_to_m_s ; uh_scale = h_scale*vel_scale
  temp_scale = US%C_to_degC ; saln_scale = US%S_to_ppt

!  if (.not.associated(CS)) return
  nz = GV%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%v_file == -1) then
      if (len_trim(CS%v_trunc_file) < 1) return
      call open_ASCII_file(CS%v_file, trim(CS%v_trunc_file), action=APPEND_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%v_file == -1) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%v_trunc_file)//'.')
        return
      endif
    endif
    file = CS%v_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

    do k=1,nz
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= vel_rpt) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= vel_rpt) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1 ; ke = nz ; write(file,'("V: Unable to set ks & ke.")')
    endif
    if (CS%full_column) then
      ks = 1 ; ke = nz
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",I0," ",I0," ",F6.2," V-velocity violation at ",I0,": ",I0,", ",I0, &
        & " (",F7.2," E ",F7.2," N) Layers ",I0," to ",I0,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), i, J, &
        G%geoLonCv(i,J), G%geoLatCv(i,J), ks, ke, US%T_to_s*dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ")', advance='no') (k) ; enddo
    write(file,'(/,"v(m):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*vm(i,J,k)) ; enddo

    if (prev_avail) then
      write(file,'(/,"v(mp): ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J,k)) ; enddo
    endif

    write(file,'(/,"v(3):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_av(i,J,k)) ; enddo
    write(file,'(/,"CFL v: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(vm(i,J,k)) * dt * G%dx_Cv(i,J)
      if (vm(i,J,k) < 0.0) then ; CFL = CFL * G%IareaT(i,j+1)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ")', advance='no') CFL
    endif ; enddo
    write(file,'(/,"CFL0 v:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    abs(vm(i,J,k)) * dt * G%IdyCv(i,J) ; enddo

    if (prev_avail) then
      write(file,'(/,"dv:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*(vm(i,J,k)-CS%v_prev(i,J,k))) ; enddo
    endif

    write(file,'(/,"CAv:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%CAv(i,J,k)) ; enddo

    write(file,'(/,"PFv:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%PFv(i,J,k)) ; enddo

    write(file,'(/,"diffv: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%diffv(i,J,k)) ; enddo

    if (associated(ADp%gradKEv)) then
      write(file,'(/,"KEv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%gradKEv(i,J,k)) ; enddo
    endif
    if (associated(ADp%rv_x_u)) then
      write(file,'(/,"Corv:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                 vel_scale*dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k)) ; enddo
    endif
    if (associated(ADp%dv_dt_visc)) then
      write(file,'(/,"vbv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*(vm(i,J,k) - dt*ADp%dv_dt_visc(i,J,k)) ; enddo

      write(file,'(/,"dvv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%dv_dt_visc(i,J,k)) ; enddo
    endif
    if (associated(ADp%dv_other)) then
      write(file,'(/,"dv_other: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*ADp%dv_other(i,J,k)) ; enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",ES10.3," ")', advance='no') h_scale*a(i,J,ks)*dt
      do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') (h_scale*a(i,J,K)*dt) ; enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hv(i,J,k) ; enddo
    endif
    if (present(str)) then
      write(file,'(/,"Stress:  ",ES10.3)', advance='no') (uh_scale*GV%RZ_to_H) * (str*dt)
    endif

    if (associated(CS%v_accel_bt)) then
      write(file,'("dvbt:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*CS%v_accel_bt(i,J,k)) ; enddo
    endif
    write(file,'(/)')

    write(file,'("h--:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i-1,j,k) ; enddo
    write(file,'(/,"h0-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"h+-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i+1,j,k) ; enddo
    write(file,'(/,"h-+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i-1,j+1,k) ; enddo
    write(file,'(/,"h0+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i,j+1,k) ; enddo
    write(file,'(/,"h++:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i+1,j+1,k) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"e-:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j+1) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j+1,k) ; enddo
    write(file,'(/,"e+:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo
    if (associated(CS%T)) then
      write(file,'(/,"T-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j,k) ; enddo
      write(file,'(/,"T+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j+1,k) ; enddo
    endif
    if (associated(CS%S)) then
      write(file,'(/,"S-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j,k) ; enddo
      write(file,'(/,"S+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j+1,k) ; enddo
    endif

    if (prev_avail) then
      write(file,'(/,"u--:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I-1,j,k) ; enddo
      write(file,'(/,"u-+:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I-1,j+1,k) ; enddo
      write(file,'(/,"u+-:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I,j,k) ; enddo
      write(file,'(/,"u++:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I,j+1,k) ; enddo
    endif

    write(file,'(/,"uh--:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I-1,j,k)*G%IdyCu(I-1,j)) ; enddo
    write(file,'(/," uhC--: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I-1,j,k) * uh_scale*0.5*(hin(i-1,j,k) + hin(i,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp--:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I-1,j,k) * uh_scale*0.5*(hin(i-1,j,k) + hin(i,j,k))) ; enddo
    endif

    write(file,'(/,"uh-+:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I-1,j+1,k)*G%IdyCu(I-1,j+1)) ; enddo
    write(file,'(/," uhC-+: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I-1,j+1,k) * uh_scale*0.5*(hin(i-1,j+1,k) + hin(i,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp-+:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I-1,j+1,k) * uh_scale*0.5*(hin(i-1,j+1,k) + hin(i,j+1,k))) ; enddo
    endif

    write(file,'(/,"uh+-:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I,j,k)*G%IdyCu(I,j)) ; enddo
    write(file,'(/," uhC+-: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I,j,k) * uh_scale*0.5*(hin(i,j,k) + hin(i+1,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp+-:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I,j,k) * uh_scale*0.5*(hin(i,j,k) + hin(i+1,j,k))) ; enddo
    endif

    write(file,'(/,"uh++:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I,j+1,k)*G%IdyCu(I,j+1)) ; enddo
    write(file,'(/," uhC++: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I,j+1,k) * uh_scale*0.5*(hin(i,j+1,k) + hin(i+1,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp++:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I,j+1,k) * uh_scale*0.5*(hin(i,j+1,k) + hin(i+1,j+1,k))) ; enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') US%Z_to_m*(G%bathyT(i,j) + G%Z_ref), US%Z_to_m*(G%bathyT(i,j+1) + G%Z_ref)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        dv = vm(i,J,k) - CS%v_prev(i,J,k)
        if (abs(dv) < 1.0e-6*US%m_s_to_L_T) dv = 1.0e-6*US%m_s_to_L_T
        Inorm(k) = 1.0 / dv
      enddo

      write(file,'(2/,"Norm:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') (vel_scale / Inorm(k)) ; enddo
      write(file,'(/,"dv:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      ((vm(i,J,k)-CS%v_prev(i,J,k)) * Inorm(k)) ; enddo
      write(file,'(/,"CAv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%CAv(i,J,k) * Inorm(k)) ; enddo
      write(file,'(/,"PFv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%PFv(i,J,k) * Inorm(k)) ; enddo
      write(file,'(/,"diffv: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%diffv(i,J,k) * Inorm(k)) ; enddo

      if (associated(ADp%gradKEu)) then
        write(file,'(/,"KEv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%gradKEv(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%rv_x_u)) then
        write(file,'(/,"Corv:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k)) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%dv_dt_visc)) then
        write(file,'(/,"dvv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%dv_dt_visc(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%dv_other)) then
        write(file,'(/,"dv_other: ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (ADp%dv_other(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(CS%v_accel_bt)) then
        write(file,'(/,"dvbt:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*CS%v_accel_bt(i,J,k) * Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    flush(file)
  endif

end procedure write_v_accel
module procedure PointAccel_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_PointAccel" ! This module's name.
  if (associated(CS)) return
  allocate(CS)

  CS%diag => diag ; CS%Time => Time

  CS%T => MIS%T ; CS%S => MIS%S
  CS%u_accel_bt => MIS%u_accel_bt ; CS%v_accel_bt => MIS%v_accel_bt
  CS%u_prev => MIS%u_prev ; CS%v_prev => MIS%v_prev
  CS%u_av => MIS%u_av ; if (.not.associated(MIS%u_av)) CS%u_av => MIS%u(:,:,:)
  CS%v_av => MIS%v_av ; if (.not.associated(MIS%v_av)) CS%v_av => MIS%v(:,:,:)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", debugging=.true.)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to zonal velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is "//&
                 "not needed.", default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to meridional velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is "//&
                 "not needed.", default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "MAX_TRUNC_FILE_SIZE_PER_PE", CS%max_writes, &
                 "The maximum number of columns of truncations that any PE "//&
                 "will write out during a run.", default=50, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_FULL_COLUMN", CS%full_column, &
                 "If true, write out the accelerations in all massive layers; otherwise "//&
                 "just document the ones with large velocities.", &
                 default=.false., debuggingParam=.true.)

  if (len_trim(dirs%output_directory) > 0) then
    if (len_trim(CS%u_trunc_file) > 0) &
      CS%u_trunc_file = trim(dirs%output_directory)//trim(CS%u_trunc_file)
    if (len_trim(CS%v_trunc_file) > 0) &
      CS%v_trunc_file = trim(dirs%output_directory)//trim(CS%v_trunc_file)
    call log_param(param_file, mdl, "output_dir/U_TRUNC_FILE", CS%u_trunc_file, debuggingParam=.true.)
    call log_param(param_file, mdl, "output_dir/V_TRUNC_FILE", CS%v_trunc_file, debuggingParam=.true.)
  endif
  CS%u_file = -1 ; CS%v_file = -1 ; CS%cols_written = 0

end procedure PointAccel_init
end submodule MOM_PointAccel_s
