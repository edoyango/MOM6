submodule (MOM_wave_drag) MOM_wave_drag_s
#include <MOM_memory.h>
  implicit none
contains
module procedure wave_drag_init
  character(len=40)  :: mdl = "MOM_wave_drag"          !< This module's name
  character(len=50)  :: filter_name_str                !< List of drag coefficients to be used
  character(len=2),  allocatable, dimension(:) :: filter_names !< Names of drag coefficients
  character(len=80)  :: var_names(4)                   !< Names of variables in wave_drag_file
  character(len=200) :: mesg
  real               :: var_scale                      !< Scaling factors of drag coefficients [nondim]
  integer            :: c
  call get_param(param_file, mdl, "N_FILTERS", CS%nf, &
                 "Number of streaming band-pass filters to be used in the simulation.", &
                 default=0, do_not_log=.true.)
  call get_param(param_file, mdl, "FILTER_NAMES", filter_name_str, &
                 "Names of streaming band-pass filters to be used in the simulation.", &
                 do_not_log=.true.)

  allocate(CS%coef_u(G%IsdB:G%IedB,G%jsd:G%jed,CS%nf)) ; CS%coef_u(:,:,:) = 0.0
  allocate(CS%coef_v(G%isd:G%ied,G%JsdB:G%JedB,CS%nf)) ; CS%coef_v(:,:,:) = 0.0
  allocate(CS%coef_uv(G%IsdB:G%IedB,G%jsd:G%jed,CS%nf)) ; CS%coef_uv(:,:,:) = 0.0
  allocate(CS%coef_vu(G%isd:G%ied,G%JsdB:G%JedB,CS%nf)) ; CS%coef_vu(:,:,:) = 0.0
  allocate(filter_names(CS%nf)) ; read(filter_name_str, *) filter_names

  CS%tensor_drag = .false.

  if (len_trim(wave_drag_file) > 0) then
    do c=1,CS%nf
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_U", &
                     var_names(1), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at u points.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_V", &
                     var_names(2), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at v points.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_UV", &
                     var_names(3), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at u points, corresponding to the off-diagonal "//&
                     "component of the wave drag tensor.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_VU", &
                     var_names(4), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at v points, corresponding to the off-diagonal "//&
                     "component of the wave drag tensor.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_SCALE", &
                     var_scale, "A scaling factor for the drag coefficient of the "//&
                     trim(filter_names(c))//" frequency.", default=1.0, units="nondim")

      if (len_trim(var_names(1))>0 .and. len_trim(var_names(2))>0 .and. var_scale>0.0) then
        call MOM_read_data(wave_drag_file, trim(var_names(1)), CS%coef_u(:,:,c), G%Domain, &
                           position=EAST_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
        call MOM_read_data(wave_drag_file, trim(var_names(2)), CS%coef_v(:,:,c), G%Domain, &
                           position=NORTH_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
        call pass_vector(CS%coef_u(:,:,c), CS%coef_v(:,:,c), G%domain, &
                         direction=To_All+SCALAR_PAIR)

        if (len_trim(var_names(3))>0 .and. len_trim(var_names(4))>0) then
          CS%tensor_drag = .true.

          call MOM_read_data(wave_drag_file, trim(var_names(3)), CS%coef_uv(:,:,c), G%Domain, &
                             position=EAST_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
          call MOM_read_data(wave_drag_file, trim(var_names(4)), CS%coef_vu(:,:,c), G%Domain, &
                             position=NORTH_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
          call pass_vector(CS%coef_uv(:,:,c), CS%coef_vu(:,:,c), G%domain, &
                           direction=To_All+SCALAR_PAIR)
        endif

        write(mesg, *) "MOM_wave_drag: ", trim(filter_names(c)), &
                       " coefficients read from file, scaling factor = ", var_scale
        call MOM_error(NOTE, trim(mesg))
      endif ! (len_trim(var_names(1))+len_trim(var_names(2))>0 .and. var_scale>0.0)
    enddo ! k=1,CS%nf
  endif ! (len_trim(wave_drag_file) > 0)

end procedure wave_drag_init
module procedure wave_drag_calc
  integer :: is, ie, js, je, i, j, c
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  drag_u(:,:) = 0.0 ; drag_v(:,:) = 0.0

  if (CS%tensor_drag) then
    call pass_vector(u(:,:,1:CS%nf), v(:,:,1:CS%nf), G%domain, direction=To_All+SCALAR_PAIR)
    !$OMP do
    do j=js,je ; do I=is-1,ie ; do c=1,CS%nf ; if (G%mask2dCu(I,j) * CS%coef_u(I,j,c) > 0.0) then
      drag_u(I,j) = drag_u(I,j) + (u(I,j,c) * CS%coef_u(I,j,c) + &
                    0.25 * ((v(i+1,J,c) + v(i,J-1,c)) + (v(i,J,c) + v(i+1,J-1,c))) * CS%coef_uv(I,j,c))
    endif ; enddo ; enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie ; do c=1,CS%nf ; if (G%mask2dCv(i,J) * CS%coef_v(i,J,c) > 0.0) then
      drag_v(i,J) = drag_v(i,J) + (v(i,J,c) * CS%coef_v(i,J,c) + &
                    0.25 * ((u(I-1,j,c) + u(I,j+1,c)) + (u(I,j,c) + u(I-1,j+1,c))) * CS%coef_vu(i,J,c))
    endif ; enddo ; enddo ; enddo
  else ! (.not.CS%tensor_drag)
    !$OMP do
    do j=js,je ; do I=is-1,ie ; do c=1,CS%nf ; if (G%mask2dCu(I,j) * CS%coef_u(I,j,c) > 0.0) then
      drag_u(I,j) = drag_u(I,j) + u(I,j,c) * CS%coef_u(I,j,c)
    endif ; enddo ; enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie ; do c=1,CS%nf ; if (G%mask2dCv(i,J) * CS%coef_v(i,J,c) > 0.0) then
      drag_v(i,J) = drag_v(i,J) + v(i,J,c) * CS%coef_v(i,J,c)
    endif ; enddo ; enddo ; enddo
  endif ! (CS%tensor_drag)

end procedure wave_drag_calc
end submodule MOM_wave_drag_s
