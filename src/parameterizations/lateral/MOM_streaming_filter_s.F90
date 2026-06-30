submodule (MOM_streaming_filter) MOM_streaming_filter_s
#include <MOM_memory.h>
  implicit none
contains
module procedure Filt_register
  type(axis_info) :: filter_axis(1)
  real, dimension(:), allocatable :: n_filters         !< Labels of filters [nondim]
  integer :: c
  CS%nf  = nf
  CS%key = key

  select case (trim(grid))
    case ('h')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('u')
      CS%is = HI%IsdB ; CS%ie = HI%IedB ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('v')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%JsdB ; CS%je = HI%JedB
    case default
      call MOM_error(FATAL, "MOM_streaming_filter: horizontal grid not supported")
  end select

  allocate(CS%s1(CS%is:CS%ie, CS%js:CS%je, nf), source=0.0)
  allocate(CS%u1(CS%is:CS%ie, CS%js:CS%je, nf), source=0.0)

  ! Register restarts for s1 and u1
  allocate(n_filters(nf))

  do c=1,nf ; n_filters(c) = c ; enddo

  call set_axis_info(filter_axis(1), "n_filters", "", "number of filters", nf, n_filters, "N", 1)

  call register_restart_field(CS%s1(:,:,:), "Filter_"//trim(key)//"_s1", .false., restart_CS, &
                              longname="Dummy variable for streaming band-pass filter", &
                              hor_grid=trim(grid), z_grid="1", t_grid="s", extra_axes=filter_axis)
  call register_restart_field(CS%u1(:,:,:), "Filter_"//trim(key)//"_u1", .false., restart_CS, &
                              longname="Output of streaming band-pass filter", &
                              hor_grid=trim(grid), z_grid="1", t_grid="s", extra_axes=filter_axis)

end procedure Filt_register
module procedure Filt_init
  character(len=40)  :: mdl = "MOM_streaming_filter"   !< This module's name
  character(len=50)  :: filter_name_str                !< List of filters to be registered
  character(len=200) :: mesg
  integer :: c
  call get_param(param_file, mdl, "FILTER_NAMES", filter_name_str, &
                 "Names of streaming band-pass filters to be used in the simulation.", &
                 fail_if_missing=.true.)
  allocate(CS%filter_names(CS%nf))
  allocate(CS%filter_omega(CS%nf))
  allocate(CS%filter_alpha(CS%nf))
  read(filter_name_str, *) CS%filter_names

  do c=1,CS%nf
    ! If filter_name_str consists of tidal constituents, use tidal frequencies.
    call get_param(param_file, mdl, "FILTER_"//trim(CS%filter_names(c))//"_OMEGA", &
                   CS%filter_omega(c), "Target frequency of the "//trim(CS%filter_names(c))//&
                   " filter. This is used if USE_FILTER is true and "//trim(CS%filter_names(c))//&
                   " is in FILTER_NAMES.", units="rad s-1", scale=US%T_to_s, default=0.0)
    call get_param(param_file, mdl, "FILTER_"//trim(CS%filter_names(c))//"_ALPHA", &
                   CS%filter_alpha(c), "Bandwidth parameter of the "//trim(CS%filter_names(c))//&
                   " filter. Must be positive.", units="nondim", fail_if_missing=.true.)

    if (CS%filter_omega(c)<=0.0) CS%filter_omega(c) = tidal_frequency(trim(CS%filter_names(c)))
    if (CS%filter_alpha(c)<=0.0) call MOM_error(FATAL, "MOM_streaming_filter: bandwidth <= 0")

    write(mesg,*) "MOM_streaming_filter: ", trim(CS%filter_names(c)), &
                  " filter registered, target frequency = ", CS%filter_omega(c), &
                  ", bandwidth = ", CS%filter_alpha(c)
    call MOM_error(NOTE, trim(mesg))
  enddo

  if (query_initialized(CS%s1, "Filter_"//trim(CS%key)//"_s1", restart_CS)) then
    write(mesg,*) "MOM_streaming_filter: Dummy variable for filter ", trim(CS%key), &
                  " found in restart files."
  else
    write(mesg,*) "MOM_streaming_filter: Dummy variable for filter ", trim(CS%key), &
                  " not found in restart files. The filter will spin up from zeros."
  endif
  call MOM_error(NOTE, trim(mesg))

  if (query_initialized(CS%u1, "Filter_"//trim(CS%key)//"_u1", restart_CS)) then
    write(mesg,*) "MOM_streaming_filter: Output of filter ", trim(CS%key), &
                  " found in restart files."
  else
    write(mesg,*) "MOM_streaming_filter: Output of filter ", trim(CS%key), &
                  " not found in restart files. The filter will spin up from zeros."
  endif
  call MOM_error(NOTE, trim(mesg))

end procedure Filt_init
module procedure Filt_accum
  real    :: now, &              !< The current model time [T ~> s]
             dt, &               !< Time step size for the filter equations [T ~> s]
             c1, c2              !< Coefficients for the filter equations [nondim]
  integer :: i, j, k
  now = time_to_real(Time, scale=US%s_to_T)

  ! Initialize CS%old_time at the first time step
  if (CS%old_time<0.0) CS%old_time = now

  ! Timestep the filter equations only if we are in a new time step
  if (CS%old_time<now) then
    dt = now - CS%old_time
    CS%old_time = now

    do k=1,CS%nf
      c1 = CS%filter_omega(k) * dt
      c2 = 1.0 - CS%filter_alpha(k) * c1

      do j=CS%js,CS%je ; do i=CS%is,CS%ie
        CS%s1(i,j,k) =  c1 *  CS%u1(i,j,k) + CS%s1(i,j,k)
        CS%u1(i,j,k) = -c1 * (CS%s1(i,j,k) - CS%filter_alpha(k) * u(i,j)) + c2 * CS%u1(i,j,k)
      enddo ; enddo
    enddo ! k=1,CS%nf
  endif ! (CS%old_time<now)

  u1 => CS%u1

end procedure Filt_accum
end submodule MOM_streaming_filter_s
