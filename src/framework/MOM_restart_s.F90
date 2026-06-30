submodule (MOM_restart) MOM_restart_s
  implicit none
contains
module procedure register_restart_field_as_obsolete
  CS%num_obsolete_vars = CS%num_obsolete_vars+1
  CS%restart_obsolete(CS%num_obsolete_vars)%field_name = field_name
  CS%restart_obsolete(CS%num_obsolete_vars)%replacement_name = replacement_name
end procedure register_restart_field_as_obsolete
module procedure register_restart_field_ptr3d
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr3d")

  CS%var_ptr3d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end procedure register_restart_field_ptr3d
module procedure register_restart_field_ptr4d
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr4d")

  CS%var_ptr4d(CS%novars)%p => f_ptr
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end procedure register_restart_field_ptr4d
module procedure register_restart_field_ptr2d
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr2d")

  CS%var_ptr2d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end procedure register_restart_field_ptr2d
module procedure register_restart_field_ptr1d
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr1d")

  CS%var_ptr1d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end procedure register_restart_field_ptr1d
module procedure register_restart_field_ptr0d
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr0d")

  CS%var_ptr0d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()

end procedure register_restart_field_ptr0d
module procedure register_restart_pair_ptr2d
  real :: a_conv, b_conv  ! Factors to multipy the a- and b-components by before they are written,
  call lock_check(CS, a_desc)
  call set_conversion_pair(a_conv, b_conv, CS%turns, conversion, scalar_pair)

  if (modulo(CS%turns, 2) == 0) then  ! This is the usual case.
    call register_restart_field(a_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(b_ptr, b_desc, mandatory, CS, conversion=b_conv)
  else
    call register_restart_field(b_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(a_ptr, b_desc, mandatory, CS, conversion=b_conv)
  endif
end procedure register_restart_pair_ptr2d
module procedure register_restart_pair_ptr3d
  real :: a_conv, b_conv  ! Factors to multipy the a- and b-components by before they are written,
  call lock_check(CS, a_desc)
  call set_conversion_pair(a_conv, b_conv, CS%turns, conversion, scalar_pair)

  if (modulo(CS%turns, 2) == 0) then  ! This is the usual case.
    call register_restart_field(a_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(b_ptr, b_desc, mandatory, CS, conversion=b_conv)
  else
    call register_restart_field(b_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(a_ptr, b_desc, mandatory, CS, conversion=b_conv)
  endif
end procedure register_restart_pair_ptr3d
module procedure register_restart_pair_ptr4d
  real :: a_conv, b_conv  ! Factors to multipy the a- and b-components by before they are written,
  call lock_check(CS, a_desc)
  call set_conversion_pair(a_conv, b_conv, CS%turns, conversion, scalar_pair)

  if (modulo(CS%turns, 2) == 0) then  ! This is the usual case.
    call register_restart_field(a_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(b_ptr, b_desc, mandatory, CS, conversion=b_conv)
  else
    call register_restart_field(b_ptr, a_desc, mandatory, CS, conversion=a_conv)
    call register_restart_field(a_ptr, b_desc, mandatory, CS, conversion=b_conv)
  endif
end procedure register_restart_pair_ptr4d
module procedure set_conversion_pair
  integer :: q_turns
  logical :: scalars
  u_conv = 1.0 ; v_conv = 1.0
  if (present(conversion)) then
    u_conv = conversion ; v_conv = conversion
  endif

  scalars = .false. ; if (present(scalar_pair)) scalars = scalar_pair
  if (scalars) return

  q_turns = modulo(turns, 4)
  if (q_turns == 1) then
    v_conv = -1.0*v_conv
  elseif (q_turns == 2) then
    u_conv = -1.0*u_conv ; v_conv = -1.0*v_conv
  elseif (q_turns == 3) then
    u_conv = -1.0*u_conv
  endif

end procedure set_conversion_pair
module procedure register_restart_field_4d
  type(vardesc) :: vd
  character(len=32), dimension(:), allocatable :: dim_names
  integer :: n, n_extradims
  if (present(extra_axes)) then
    n_extradims = size(extra_axes)
    allocate(dim_names(n_extradims+2))
    dim_names(1) = ""
    dim_names(2) = ""
    do n=3,n_extradims+2
      dim_names(n) = extra_axes(n-2)%name
    enddo
  endif

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_4d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

  if (present(extra_axes)) then
    vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                  z_grid=z_grid, t_grid=t_grid, dim_names=dim_names, extra_axes=extra_axes)
  else
    vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                  z_grid=z_grid, t_grid=t_grid)
  endif

  call register_restart_field_ptr4d(f_ptr, vd, mandatory, CS, conversion)

end procedure register_restart_field_4d
module procedure register_restart_field_3d
  type(vardesc) :: vd
  character(len=32), dimension(:), allocatable :: dim_names
  integer :: n, n_extradims
  if (present(extra_axes)) then
    n_extradims = size(extra_axes)
    allocate(dim_names(n_extradims+2))
    dim_names(1) = ""
    dim_names(2) = ""
    do n=3,n_extradims+2
      dim_names(n) = extra_axes(n-2)%name
    enddo
  endif

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

  if (present(extra_axes)) then
    vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                  z_grid=z_grid, t_grid=t_grid, dim_names=dim_names, extra_axes=extra_axes)
  else
    vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                  z_grid=z_grid, t_grid=t_grid)
  endif

  call register_restart_field_ptr3d(f_ptr, vd, mandatory, CS, conversion)

end procedure register_restart_field_3d
module procedure register_restart_field_2d
  type(vardesc) :: vd
  character(len=8) :: Zgrid
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_2d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  zgrid = '1' ; if (present(z_grid)) zgrid = z_grid

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=zgrid, t_grid=t_grid)

  call register_restart_field_ptr2d(f_ptr, vd, mandatory, CS, conversion)

end procedure register_restart_field_2d
module procedure register_restart_field_1d
  type(vardesc) :: vd
  character(len=8) :: hgrid
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  hgrid = '1' ; if (present(hor_grid)) hgrid = hor_grid

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hgrid, &
                z_grid=z_grid, t_grid=t_grid)

  call register_restart_field_ptr1d(f_ptr, vd, mandatory, CS, conversion)

end procedure register_restart_field_1d
module procedure register_restart_field_0d
  type(vardesc) :: vd
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_0d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid='1', &
                z_grid='1', t_grid=t_grid)

  call register_restart_field_ptr0d(f_ptr, vd, mandatory, CS, conversion)

end procedure register_restart_field_0d
module procedure query_initialized_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (trim(name) == CS%restart_field(m)%var_name) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if ((n==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE,"MOM_restart: Unknown restart variable "//name// &
                        " queried for initialization.")

  if ((is_root_pe()) .and. query_initialized) &
    call MOM_error(NOTE,"MOM_restart: "//name// &
                         " initialization confirmed by name.")

end procedure query_initialized_name
module procedure query_initialized_0d
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo

end procedure query_initialized_0d
module procedure query_initialized_1d
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo

end procedure query_initialized_1d
module procedure query_initialized_2d
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo

end procedure query_initialized_2d
module procedure query_initialized_3d
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo

end procedure query_initialized_3d
module procedure query_initialized_4d
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo

end procedure query_initialized_4d
module procedure query_initialized_0d_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end procedure query_initialized_0d_name
module procedure query_initialized_1d_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end procedure query_initialized_1d_name
module procedure query_initialized_2d_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end procedure query_initialized_2d_name
module procedure query_initialized_3d_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end procedure query_initialized_3d_name
module procedure query_initialized_4d_name
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end procedure query_initialized_4d_name
module procedure set_initialized_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (trim(name) == trim(CS%restart_field(m)%var_name)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if ((m==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE,"MOM_restart: Unknown restart variable "//name// &
                        " used in set_initialized call.")

end procedure set_initialized_name
module procedure set_initialized_0d_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if (m==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    call set_initialized_name(name, CS)
  endif

end procedure set_initialized_0d_name
module procedure set_initialized_1d_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if (m==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    call set_initialized_name(name, CS)
  endif

end procedure set_initialized_1d_name
module procedure set_initialized_2d_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if (m==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    call set_initialized_name(name, CS)
  endif

end procedure set_initialized_2d_name
module procedure set_initialized_3d_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if (m==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    call set_initialized_name(name, CS)
  endif

end procedure set_initialized_3d_name
module procedure set_initialized_4d_name
  integer :: m
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "set_initialized: Module must be initialized before it is used.")

  do m=1,CS%novars ; if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
    CS%restart_field(m)%initialized = .true. ; exit
  endif ; enddo

  if (m==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    call set_initialized_name(name, CS)
  endif

end procedure set_initialized_4d_name
module procedure only_read_restart_field_4d
  character(len=:), allocatable :: file_path ! The full path to the file with the variable
  logical :: found     ! True if the variable was found.
  logical :: is_global ! True if the variable is in a global file.
  found = find_var_in_restart_files(varname, G, CS, file_path, filename, directory, is_global)

  if (found) then
    call MOM_read_data(file_path, varname, f_ptr, G%domain, timelevel=1, position=position, &
                       scale=scale, global_file=is_global)
  endif
  if (present(success)) success = found

end procedure only_read_restart_field_4d
module procedure only_read_restart_field_3d
  character(len=:), allocatable :: file_path ! The full path to the file with the variable
  logical :: found     ! True if the variable was found.
  logical :: is_global ! True if the variable is in a global file.
  found = find_var_in_restart_files(varname, G, CS, file_path, filename, directory, is_global)

  if (found) then
    call MOM_read_data(file_path, varname, f_ptr, G%domain, timelevel=1, position=position, &
                       scale=scale, global_file=is_global)
  endif
  if (present(success)) success = found

end procedure only_read_restart_field_3d
module procedure only_read_restart_field_2d
  character(len=:), allocatable :: file_path ! The full path to the file with the variable
  logical :: found     ! True if the variable was found.
  logical :: is_global ! True if the variable is in a global file.
  found = find_var_in_restart_files(varname, G, CS, file_path, filename, directory, is_global)

  if (found) then
    call MOM_read_data(file_path, varname, f_ptr, G%domain, timelevel=1, position=position, &
                       scale=scale, global_file=is_global)
  endif
  if (present(success)) success = found

end procedure only_read_restart_field_2d
module procedure only_read_restart_pair_3d
  character(len=:), allocatable :: file_path_a ! The full path to the file with the first variable
  character(len=:), allocatable :: file_path_b ! The full path to the file with the second variable
  integer :: a_pos, b_pos       ! A coded position for the two variables.
  logical :: a_found, b_found   ! True if the variables were found.
  logical :: global_a, global_b ! True if the variables are in global files.
  a_found = find_var_in_restart_files(a_name, G, CS, file_path_a, filename, directory, global_a)
  b_found = find_var_in_restart_files(b_name, G, CS, file_path_b, filename, directory, global_b)

  a_pos = EAST_FACE ; b_pos = NORTH_FACE
  if (present(stagger)) then ; select case (stagger)
    case (AGRID)    ; a_pos = CENTER ; b_pos = CENTER
    case (BGRID_NE) ; a_pos = CORNER ; b_pos = CORNER
    case (CGRID_NE) ; a_pos = EAST_FACE ; b_pos = NORTH_FACE
    case default    ; a_pos = EAST_FACE ; b_pos = NORTH_FACE
  end select ; endif

  if (a_found .and. b_found) then
    call MOM_read_data(file_path_a, a_name, a_ptr, G%domain, timelevel=1, position=a_pos, &
                       scale=scale, global_file=global_b, file_may_be_4d=.true.)
    call MOM_read_data(file_path_b, b_name, b_ptr, G%domain, timelevel=1, position=b_pos, &
                       scale=scale, global_file=global_b, file_may_be_4d=.true.)
  endif
  if (present(success)) success = (a_found .and. b_found)

end procedure only_read_restart_pair_3d
module procedure find_var_in_restart_files
  character(len=240), allocatable, dimension(:) :: file_paths ! The possible file names.
  character(len=:), allocatable :: dir ! The directory to read from.
  character(len=:), allocatable :: fname ! The list of file names.
  logical, allocatable, dimension(:) :: global_file  ! True if the file is global
  integer :: n, num_files
  dir = "./INPUT/" ; if (present(directory)) dir = trim(directory)

  ! Set the default return values.
  found = .false.
  file_path = ""
  if (present(is_global)) is_global = .false.

  fname = 'r'
  if (present(filename)) then
    if (.not.((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F'))) fname = filename
  endif

  num_files = get_num_restart_files(fname, dir, G, CS)
  if (num_files == 0) return
  allocate(file_paths(num_files), global_file(num_files))
  num_files = open_restart_units(fname, dir, G, CS, file_paths=file_paths, global_files=global_file)

  do n=1,num_files ; if (field_exists(file_paths(n), varname, MOM_Domain=G%domain)) then
    found = .true.
    file_path = file_paths(n)
    if (present(is_global)) is_global = global_file(n)
    exit
  endif ; enddo

  deallocate(file_paths, global_file)

end procedure find_var_in_restart_files
module procedure copy_restart_var_3d
  logical :: keep_rotation
  character(len=256) :: size_msg  !< The array sizes
  integer :: m, n
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  keep_rotation = .true. ; if (present(unrotate)) keep_rotation = .not.unrotate

  n = CS%novars+1
  do m=1,CS%novars
    if (trim(name) == CS%restart_field(m)%var_name) then
      if (.not.associated(CS%var_ptr3d(m)%p)) then
        call MOM_error(FATAL, "MOM_restart: copy_restart_var(_3d) "//&
                      "attempted to copy restart variable "//name//" with the wrong rank.")
      elseif (CS%restart_field(m)%initialized) then
        if (CS%turns == 0 .or. keep_rotation) then
          if ( size_mismatch_3d(var, CS%var_ptr3d(m)%p, CS%turns, size_msg) ) &
            call MOM_error(FATAL, "MOM_restart: copy_restart_var(_3d) "//&
                      "attempted to copy restart variable "//name//" with the wrong sizes, "//trim(size_msg))

          var(:,:,:) = CS%var_ptr3d(m)%p(:,:,:)
        else
          call rotate_array(CS%var_ptr3d(m)%p, -CS%turns, var)
        endif
      else
        call MOM_error(NOTE, "MOM_restart: copy_restart_var(_3d) "//&
                      "attempted to copy uninitialized restart variable "//name//".")
      endif
      n = m ; exit
    endif
  enddo
  if ((n==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE, "MOM_restart: copy_restart_var(_3d) "//&
                  "attempted to copy unknown restart variable "//name//".")

end procedure copy_restart_var_3d
module procedure copy_restart_vector_3d
  logical :: keep_rotation, scalars
  character(len=256) :: size_msg  !< The array sizes
  integer :: m, n_u, n_v
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  keep_rotation = .true. ; if (present(unrotate)) keep_rotation = .not.unrotate

  n_u = CS%novars+1 ; n_v = CS%novars+1
  do m=1,CS%novars
    if (trim(u_name) == CS%restart_field(m)%var_name) then
      if (.not.associated(CS%var_ptr3d(m)%p)) then
        call MOM_error(FATAL, "MOM_restart: copy_restart_vector(_3d) "//&
                      "attempted to copy restart variable "//trim(u_name)//" with the wrong rank.")
      elseif (CS%restart_field(m)%initialized) then
        n_u = m
      else
        call MOM_error(NOTE, "MOM_restart: copy_restart_vector(_3d) "//&
                      "attempted to copy uninitialized restart variable "//trim(u_name)//".")
        n_u = -1
      endif
    endif
    if (trim(v_name) == CS%restart_field(m)%var_name) then
      if (.not.associated(CS%var_ptr3d(m)%p)) then
        call MOM_error(FATAL, "MOM_restart: copy_restart_vector(_3d) "//&
                      "attempted to copy restart variable "//trim(v_name)//" with the wrong rank.")
      elseif (CS%restart_field(m)%initialized) then
        n_v = m
      else
        call MOM_error(NOTE, "MOM_restart: copy_restart_vector(_3d) "//&
                      "attempted to copy uninitialized restart variable "//trim(v_name)//".")
        n_v = -1
      endif
    endif
  enddo
  if ((n_u==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE, "MOM_restart: copy_restart_vector(_3d) "//&
                  "attempted to copy unknown restart variable "//trim(u_name)//".")
  if ((n_v==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE, "MOM_restart: copy_restart_vector(_3d) "//&
                  "attempted to copy unknown restart variable "//trim(v_name)//".")

  if ((n_u>0) .and. (n_u<=CS%novars) .and. (n_v>0) .and. (n_v<=CS%novars)) then
    ! Now actually update the vector.
    if ( size_mismatch_3d(u_var, CS%var_ptr3d(n_u)%p, CS%turns, size_msg) ) &
      call MOM_error(FATAL, "MOM_restart: copy_restart_vector(_3d) "//&
                "attempted to copy restart variable "//trim(u_name)//" with the wrong sizes, "//trim(size_msg))
    if ( size_mismatch_3d(v_var, CS%var_ptr3d(n_v)%p, CS%turns, size_msg) ) &
      call MOM_error(FATAL, "MOM_restart: copy_restart_vector(_3d) "//&
                "attempted to copy restart variable "//trim(v_name)//" with the wrong sizes, "//trim(size_msg))

    if (CS%turns == 0 .or. keep_rotation) then
      u_var(:,:,:) = CS%var_ptr3d(n_u)%p(:,:,:)
      v_var(:,:,:) = CS%var_ptr3d(n_v)%p(:,:,:)
    else
      scalars = .false. ; if (present(scalar_pair)) scalars = scalar_pair
      if ((modulo(CS%turns, 2) == 0) .and. scalars) then
        call rotate_array_pair(CS%var_ptr3d(n_u)%p, CS%var_ptr3d(n_v)%p, -CS%turns, u_var, v_var)
      elseif (modulo(CS%turns, 2) == 0) then
        call rotate_vector(CS%var_ptr3d(n_u)%p, CS%var_ptr3d(n_v)%p, -CS%turns, u_var, v_var)
      elseif (scalars) then  ! This is less common
        call rotate_array_pair(CS%var_ptr3d(n_v)%p, CS%var_ptr3d(n_u)%p, -CS%turns, u_var, v_var)
      else
        call rotate_vector(CS%var_ptr3d(n_v)%p, CS%var_ptr3d(n_u)%p, -CS%turns, u_var, v_var)
      endif
    endif
  endif

end procedure copy_restart_vector_3d
module procedure size_mismatch_3d
  if (modulo(turns, 2) == 0) then
    size_mismatch_3d = ( (size(var_a,1) /= size(var_b,1)) .or. &
                         (size(var_a,2) /= size(var_b,2)) .or. &
                         (size(var_a,3) /= size(var_b,3)) )
  else
    size_mismatch_3d = ( (size(var_a,1) /= size(var_b,2)) .or. &
                         (size(var_a,2) /= size(var_b,1)) .or. &
                         (size(var_a,3) /= size(var_b,3)) )
  endif
  write(size_msg, '(3(1x,I0), " vs ", 3(1x,I0))') size(var_a,1), size(var_a,2), size(var_a,3), &
                                                  size(var_b,1), size(var_b,2), size(var_b,3)
end procedure size_mismatch_3d
module procedure save_restart
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
  type(MOM_field) :: fields(CS%max_fields) ! Opaque types containing metadata describing
  character(len=512) :: restartpath     ! The restart file path (dir/file).
  character(len=256) :: restartname     ! The restart file name (no dir).
  character(len=8)   :: suffix          ! A suffix (like _2) that is appended
  integer(kind=int64) :: var_sz, size_in_file ! The size in bytes of each variable
  integer(kind=int64), parameter :: max_file_size = 4294967292_int64 ! The maximum size in bytes for the
  integer :: start_var, next_var        ! The starting variables of the
  type(MOM_infra_file) :: IO_handle     ! The I/O handle of the open fileset
  integer :: m, nz, na
  integer :: num_files                  ! The number of restart files that will be used.
  integer :: seconds, days, year, month, hour, minute
  character(len=8) :: z_grid, t_grid    ! Variable grid info.
  integer :: pos                        ! A coded integer indicating the horizontal staggering of a variable
  real :: conv                          ! Shorthand for the conversion factor [a A-1 ~> 1]
  real :: restart_time                  ! The model time at whic the restart file is being written [days]
  character(len=32) :: filename_appendix = '' ! Appendix to filename for ensemble runs
  integer :: length                     ! The length of a text string.
  character(len=256) :: mesg, var_name
  integer(kind=int64) :: check_val(CS%max_fields,1)
  logical :: verbose
  integer :: isL, ieL, jsL, jeL
  integer :: turns                      ! Number of quarter turns from input to model domain
  integer, parameter :: nmax_extradims = 5
  type(axis_info), dimension(:), allocatable :: extra_axes
  turns = CS%turns

  allocate (extra_axes(nmax_extradims))

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "save_restart: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)
  verbose = (is_root_pe() .and. (MOM_get_verbosity() >= 7))

  ! With parallel read & write, it is possible to disable the following...
  num_files = 0
  next_var = 0
  nz = 1 ; if (present(GV)) nz = GV%ke

  restart_time = time_type_to_real(time) / 86400.0

  restartname = trim(CS%restartfile)
  if (present(filename)) restartname = trim(filename)
  if (PRESENT(time_stamped)) then ; if (time_stamped) then
    call get_date(time, year, month, days, hour, minute, seconds)
    ! Compute the year-day, because I don't like months. - RWH
    do m=1,month-1
      days = days + days_in_month(set_date(year, m, 2, 0, 0, 0))
    enddo
    seconds = seconds + 60*minute + 3600*hour
    if (year <= 9999) then
      write(restartname,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
    elseif (year <= 99999) then
      write(restartname,'("_Y",I5.5,"_D",I3.3,"_S",I5.5)') year, days, seconds
    else
      write(restartname,'("_Y",I10.10,"_D",I3.3,"_S",I5.5)') year, days, seconds
    endif
    restartname = trim(CS%restartfile)//trim(restartname)
  endif ; endif

  ! Determine if there is a filename_appendix (used for ensemble runs).
  call get_filename_appendix(filename_appendix)
  if (len_trim(filename_appendix) > 0) then
    length = len_trim(restartname)
    if (restartname(length-2:length) == '.nc') then
      restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
    else
      restartname = restartname(1:length)  //'.'//trim(filename_appendix)
    endif
  endif

  next_var = 1
  do while (next_var <= CS%novars )
    start_var = next_var
    size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)

    do m=start_var,CS%novars
      call query_vardesc(CS%restart_field(m)%vars, position=pos, &
                         z_grid=z_grid, t_grid=t_grid, caller="save_restart", &
                         extra_axes=extra_axes)

      var_sz = get_variable_byte_size(pos, z_grid, t_grid, G, nz)
      ! factor in size of extra axes, or multiply by 1
      do na=1,nmax_extradims
        var_sz = var_sz*extra_axes(na)%ax_size
      enddo

      if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
        size_in_file = size_in_file + var_sz
      else ; exit
      endif

    enddo
    next_var = m

    restartpath = trim(directory) // trim(restartname)

    write(suffix,'("_",I0)') num_files

    length = len_trim(restartpath)
    if (length < 3) then  ! This case is very uncommon but this test avoids segmentation-faults.
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    elseif (restartpath(length-2:length) == ".nc") then
      if (num_files > 0) restartpath = restartpath(1:length-3)//trim(suffix)//".nc"
    else
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    endif

    do m=start_var,next_var-1
      vars(m-start_var+1) = CS%restart_field(m)%vars
    enddo
    call query_vardesc(vars(1), t_grid=t_grid, position=pos, caller="save_restart")
    t_grid = adjustl(t_grid)
    if (t_grid(1:1) /= 'p') &
      call modify_vardesc(vars(1), t_grid='s', caller="save_restart")

    !Prepare the checksum of the restart fields to be written to restart files
    do m=start_var,next_var-1

      call query_vardesc(vars(m), position=pos, name=var_name, caller="save_restart")
      if (modulo(turns, 2) == 0) then
        call get_checksum_loop_ranges(G, CS, pos, isL, ieL, jsL, jeL)
      else   ! Note that G is always the unrotated grid as it is seen by the driver level.
        call get_checksum_loop_ranges(G, CS, pos, jsL, jeL, isL, ieL)
      endif
      if (verbose) then
        if (pos == CENTER) then
          write(mesg, '(" is in CENTER position, checksum range",4(1x,I0))') isL, ieL, jsL, jeL
        elseif (pos == CORNER) then
          write(mesg, '(" is in CORNER position, checksum range",4(1x,I0))') isL, ieL, jsL, jeL
        elseif (pos == NORTH_FACE) then
          write(mesg, '(" is in NORTH_FACE position, checksum range",4(1x,I0))') isL, ieL, jsL, jeL
        elseif (pos == EAST_FACE) then
          write(mesg, '(" is in EAST_FACE position, checksum range",4(1x,I0))') isL, ieL, jsL, jeL
        else
          write(mesg, '(" is in another position, ",I0,", checksum range",4(1x,I0))') pos, isL, ieL, jsL, jeL
        endif
        call MOM_mesg(trim(var_name)//mesg)
      endif

      conv = CS%restart_field(m)%conv
      if (associated(CS%var_ptr3d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:), turns=-turns, unscale=conv)
      elseif (associated(CS%var_ptr2d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL), turns=-turns, unscale=conv)
      elseif (associated(CS%var_ptr4d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:), turns=-turns, unscale=conv)
      elseif (associated(CS%var_ptr1d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr1d(m)%p(:), unscale=conv)
      elseif (associated(CS%var_ptr0d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr0d(m)%p, pelist=(/PE_here()/), unscale=conv)
      endif
    enddo

    if (CS%parallel_restartfiles) then
      call create_MOM_file(IO_handle, trim(restartpath), vars, next_var-start_var, &
          fields, MULTIPLE, G=G, GV=GV, checksums=check_val, extra_axes=extra_axes)
    else
      call create_MOM_file(IO_handle, trim(restartpath), vars, next_var-start_var, &
          fields, SINGLE_FILE, G=G, GV=GV, checksums=check_val, extra_axes=extra_axes)
    endif

    do m=start_var,next_var-1
      if (associated(CS%var_ptr3d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr3d(m)%p, &
                             restart_time, unscale=CS%restart_field(m)%conv, turns=-turns, &
                             zero_zeros=CS%unsigned_zeros)
      elseif (associated(CS%var_ptr2d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr2d(m)%p, &
                             restart_time, unscale=CS%restart_field(m)%conv, turns=-turns, &
                             zero_zeros=CS%unsigned_zeros)
      elseif (associated(CS%var_ptr4d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr4d(m)%p, &
                             restart_time, unscale=CS%restart_field(m)%conv, turns=-turns, &
                             zero_zeros=CS%unsigned_zeros)
      elseif (associated(CS%var_ptr1d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr1d(m)%p, &
                             restart_time, unscale=CS%restart_field(m)%conv, &
                             zero_zeros=CS%unsigned_zeros)
      elseif (associated(CS%var_ptr0d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr0d(m)%p, &
                             restart_time, unscale=CS%restart_field(m)%conv, &
                             zero_zeros=CS%unsigned_zeros)
      endif
    enddo

    call IO_handle%close()

    num_files = num_files+1

  enddo

  if (present(num_rest_files)) num_rest_files = num_files

end procedure save_restart
module procedure restore_state
  real :: scale  ! A scaling factor for reading a field [A a-1 ~> 1] to convert
  real :: conv   ! The output conversion factor for writing a field [a A-1 ~> 1]
  character(len=512) :: mesg      ! A message for warnings.
  character(len=80) :: varname    ! A variable's name.
  integer :: num_file        ! The number of files (restart files and others
  integer :: i, n, m, missing_fields
  integer :: isL, ieL, jsL, jeL
  integer :: nvar, ntime, pos
  type(MOM_infra_file) :: IO_handles(CS%max_fields) ! The I/O units of all open files.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.
  real    :: t1, t2 ! Two times from the start of different files [days].
  real, allocatable :: time_vals(:)  ! Times from a file extracted with getl_file_times [days]
  type(MOM_field), allocatable :: fields(:)
  logical            :: is_there_a_checksum ! Is there a valid checksum that should be checked.
  integer(kind=int64) :: checksum_file  ! The checksum value recorded in the input file.
  integer(kind=int64) :: checksum_data  ! The checksum value for the data that was read in.
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")

  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! Get NetCDF ids for all of the restart files.
  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_file = open_restart_units('r', directory, G, CS, IO_handles=IO_handles, &
                     file_paths=unit_path, global_files=unit_is_global)
  else
    num_file = open_restart_units(filename, directory, G, CS, IO_handles=IO_handles, &
                     file_paths=unit_path, global_files=unit_is_global)
  endif

  if (num_file == 0) then
    write(mesg,'("Unable to find any restart files specified by  ",A,"  in directory ",A,".")') &
                  trim(filename), trim(directory)
    call MOM_error(FATAL,"MOM_restart: "//mesg)
  endif

  ! Get the time from the first file in the list that has one.
  do n=1,num_file
    call IO_handles(n)%get_file_times(time_vals, ntime)
    if (ntime < 1) cycle

    t1 = time_vals(1)
    deallocate(time_vals)

    day = real_to_time(t1*86400.0)
    exit
  enddo

  if (n>num_file) call MOM_error(WARNING, "MOM_restart: No times found in restart files.")

  ! Check the remaining files for different times and issue a warning
  ! if they differ from the first time.
    do m = n+1,num_file
      call IO_handles(n)%get_file_times(time_vals, ntime)
      if (ntime < 1) cycle

      t2 = time_vals(1)
      deallocate(time_vals)

      if (t1 /= t2 .and. is_root_PE()) then
        write(mesg,'("WARNING: Restart file ",I0," has time ",F10.4,"whereas &
              &simulation is restarted at ",F10.4," (differing by ",F10.4,").")') &
              m, t1, t2, t1-t2
        call MOM_error(WARNING, "MOM_restart: "//mesg)
      endif
    enddo

  ! Read each variable from the first file in which it is found.
  do n=1,num_file
    call IO_handles(n)%get_file_info(nvar=nvar)

    allocate(fields(nvar))
    call IO_handles(n)%get_file_fields(fields(1:nvar))

    do m=1, nvar
      call IO_handles(n)%get_field_atts(fields(m), name=varname)
      do i=1,CS%num_obsolete_vars
        if (adjustl(lowercase(trim(varname))) == adjustl(lowercase(trim(CS%restart_obsolete(i)%field_name)))) then
            call MOM_error(FATAL, "MOM_restart restore_state: Attempting to use obsolete restart field "//&
                           trim(varname)//" - the new corresponding restart field is "//&
                           trim(CS%restart_obsolete(i)%replacement_name))
        endif
      enddo
    enddo

    missing_fields = 0

    do m=1,CS%novars
      if (CS%restart_field(m)%initialized) cycle
      call query_vardesc(CS%restart_field(m)%vars, position=pos, caller="restore_state")
      conv = CS%restart_field(m)%conv
      if (conv == 0.0) then ; scale = 1.0 ; else ; scale = 1.0 / conv ; endif

      if (modulo(CS%turns, 2) == 0) then
        call get_checksum_loop_ranges(G, CS, pos, isL, ieL, jsL, jeL)
      else   ! Note that G is always the unrotated grid as it is used during initialization.
        call get_checksum_loop_ranges(G, CS, pos, jsL, jeL, isL, ieL)
      endif
      do i=1, nvar
        call IO_handles(n)%get_field_atts(fields(i), name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(CS%restart_field(m)%var_name))) then
          checksum_data = -1
          if (CS%checksum_required) then
            call IO_handles(n)%read_field_chksum(fields(i), checksum_file, is_there_a_checksum)
          else
            checksum_file = -1
            is_there_a_checksum = .false. ! Do not need to do data checksumming.
          endif

          if (associated(CS%var_ptr1d(m)%p))  then
            ! Read a 1d array, which should be invariant to domain decomposition.
            call MOM_read_data(unit_path(n), varname, CS%var_ptr1d(m)%p, &
                               timelevel=1, scale=scale, MOM_Domain=G%Domain)
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr1d(m)%p(:), unscale=conv)
          elseif (associated(CS%var_ptr0d(m)%p)) then ! Read a scalar...
            call MOM_read_data(unit_path(n), varname, CS%var_ptr0d(m)%p, &
                               timelevel=1, scale=scale, MOM_Domain=G%Domain)
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr0d(m)%p, pelist=(/PE_here()/), unscale=conv)
          elseif (associated(CS%var_ptr2d(m)%p)) then  ! Read a 2d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                                 G%Domain, timelevel=1, position=pos, scale=scale, turns=CS%turns)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 2-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p,no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL), unscale=conv)
          elseif (associated(CS%var_ptr3d(m)%p)) then  ! Read a 3d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                                 G%Domain, timelevel=1, position=pos, scale=scale, turns=CS%turns)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 3-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:), unscale=conv)
          elseif (associated(CS%var_ptr4d(m)%p)) then  ! Read a 4d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                                 G%Domain, timelevel=1, position=pos, scale=scale, &
                                 global_file=unit_is_global(n), turns=CS%turns)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 4-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:), unscale=conv)
          else
            call MOM_error(FATAL, "MOM_restart restore_state: No pointers set for "//trim(varname))
          endif

          if (is_root_pe() .and. is_there_a_checksum .and. (checksum_file /= checksum_data)) then
             write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field "// trim(varname)//" ",checksum_data,&
                                          " does not match value ", checksum_file, &
                                          " stored in "//trim(unit_path(n)//"." )
             call MOM_error(FATAL, "MOM_restart(restore_state): "//trim(mesg) )
          endif

          CS%restart_field(m)%initialized = .true.
          exit ! Start search for next restart variable.
        endif
      enddo
      if (i>nvar) missing_fields = missing_fields+1
    enddo

    deallocate(fields)
    if (missing_fields == 0) exit
  enddo

  do n=1,num_file
    call IO_handles(n)%close()
  enddo

  ! Check whether any mandatory fields have not been found.
  CS%restart = .true.
  do m=1,CS%novars
    if (.not.(CS%restart_field(m)%initialized)) then
      CS%restart = .false.
      if (CS%restart_field(m)%mand_var) then
        call MOM_error(FATAL,"MOM_restart: Unable to find mandatory variable " &
                       //trim(CS%restart_field(m)%var_name)//" in restart files.")
      endif
    endif
  enddo

  ! Lock the restart registry so that no further variables can be registered.
  CS%locked = .true.

end procedure restore_state
module procedure restart_files_exist
  integer :: num_files
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "restart_files_exist: Module must be initialized before it is used.")

  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_files = get_num_restart_files('r', directory, G, CS)
  else
    num_files = get_num_restart_files(filename, directory, G, CS)
  endif
  restart_files_exist = (num_files > 0)
end procedure restart_files_exist
module procedure determine_is_new_run
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run: Module must be initialized before it is used.")

  if (LEN_TRIM(filename) > 1) then
    CS%new_run = .false.
  elseif (LEN_TRIM(filename) == 0) then
    CS%new_run = .true.
  elseif (filename(1:1) == 'n') then
    CS%new_run = .true.
  elseif (filename(1:1) == 'F') then
    CS%new_run = (get_num_restart_files('r', directory, G, CS) == 0)
  else
    CS%new_run = .false.
  endif

  CS%new_run_set = .true.
  is_new_run = CS%new_run
end procedure determine_is_new_run
module procedure is_new_run
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "is_new_run: Module must be initialized before it is used.")

  if (.not.CS%new_run_set) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run must be called for a restart file before is_new_run.")

  is_new_run = CS%new_run
end procedure is_new_run
module procedure open_restart_units
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
  integer :: num_restart     ! The number of restart files that have already
  integer :: start_char      ! The location of the starting character in the
  integer :: nf              ! The number of files that have been found so far
  integer :: m, length
  logical :: still_looking   ! If true, the code is still looking for automatically named files
  logical :: fexists         ! True if a file has been found
  character(len=32) :: filename_appendix = '' ! Filename appendix for ensemble runs
  character(len=80) :: restartname
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "open_restart_units: Module must be initialized before it is used.")

  ! Get NetCDF ids for all of the restart files.
  num_restart = 0 ; nf = 0 ; start_char = 1
  do while (start_char <= len_trim(filename) )
    do m=start_char,len_trim(filename)
      if (filename(m:m) == ' ') exit
    enddo
    fname = filename(start_char:m-1)
    start_char = m
    do while (start_char <= len_trim(filename))
      if (filename(start_char:start_char) == ' ') then
        start_char = start_char + 1
      else
        exit
      endif
    enddo

    if ((fname(1:1)=='r') .and. ( len_trim(fname) == 1)) then
      still_looking = (num_restart <= 0) ! Avoid going through the file list twice.
      do while (still_looking)
        restartname = trim(CS%restartfile)

        ! Determine if there is a filename_appendix (used for ensemble runs).
        call get_filename_appendix(filename_appendix)
        if (len_trim(filename_appendix) > 0) then
          length = len_trim(restartname)
          if (restartname(length-2:length) == '.nc') then
            restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
          else
            restartname = restartname(1:length)  //'.'//trim(filename_appendix)
          endif
        endif
        filepath = trim(directory) // trim(restartname)

        write(suffix,'("_",I0)') num_restart
        if (num_restart > 0) filepath = trim(filepath) // suffix

        filepath = trim(filepath)//".nc"

        num_restart = num_restart + 1
        ! Look for a global netCDF file.
        inquire(file=filepath, exist=fexists)
        if (fexists) then
          nf = nf + 1
          if (present(IO_handles)) &
            call IO_handles(nf)%open(trim(filepath), READONLY_FILE, &
                MOM_domain=G%Domain, threading=MULTIPLE, fileset=SINGLE_FILE)
          if (present(global_files)) global_files(nf) = .true.
          if (present(file_paths)) file_paths(nf) = filepath
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          fexists = file_exists(filepath, G%Domain)
          if (fexists) then
            nf = nf + 1
            if (present(IO_handles)) &
              call IO_handles(nf)%open(trim(filepath), READONLY_FILE, MOM_domain=G%Domain)
            if (present(global_files)) global_files(nf) = .false.
            if (present(file_paths)) file_paths(nf) = filepath
          endif
        endif

        if (fexists) then
          if (is_root_pe() .and. (present(IO_handles))) &
            call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
        else
          still_looking = .false. ; exit
        endif
      enddo ! while (still_looking) loop
    else
      filepath = trim(directory)//trim(fname)
      inquire(file=filepath, exist=fexists)
      if (.not. fexists) filepath = trim(filepath)//".nc"

      inquire(file=filepath, exist=fexists)
      if (fexists) then
        nf = nf + 1
        if (present(IO_handles)) &
          call IO_handles(nf)%open(trim(filepath), READONLY_FILE, &
              MOM_Domain=G%Domain, threading=MULTIPLE, fileset=SINGLE_FILE)
        if (present(global_files)) global_files(nf) = .true.
        if (present(file_paths)) file_paths(nf) = filepath
        if (is_root_pe() .and. (present(IO_handles))) &
          call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
      else
        if (present(IO_handles)) &
          call MOM_error(WARNING, "MOM_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < len_trim(filename)) loop
  num_files = nf

end procedure open_restart_units
module procedure get_num_restart_files
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_restart " // &
      "get_num_restart_files: Module must be initialized before it is used.")

  ! This call uses open_restart_units without the optional arguments needed to actually
  ! open the files to determine the number of restart files.
  num_files = open_restart_units(filenames, directory, G, CS, file_paths=file_paths)

end procedure get_num_restart_files
module procedure restart_init
  logical :: rotate_index
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_restart"   ! This module's name.
  logical :: all_default   ! If true, all parameters are using their default values.
  if (associated(CS)) then
    call MOM_error(WARNING, "restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  CS%initialized = .true.

  ! Determine whether all paramters are set to their default values.
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, default=100, do_not_log=.true.)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "RESTART_SYMMETRIC_CHECKSUMS", CS%symmetric_checksums, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "RESTART_UNSIGNED_ZEROS", CS%unsigned_zeros, &
                 default=.false., do_not_log=.true.)
  all_default = ((.not.CS%parallel_restartfiles) .and. (CS%max_fields == 100) .and. &
                 (CS%checksum_required) .and. (.not.CS%symmetric_checksums) .and. (.not.CS%unsigned_zeros))
  if (.not.present(restart_root)) then
    call get_param(param_file, mdl, "RESTARTFILE", CS%restartfile, &
                   default="MOM.res", do_not_log=.true.)
    all_default = (all_default .and. (trim(CS%restartfile) == trim("MOM.res")))
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", all_default=all_default)
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 "If true, the IO layout is used to group processors that write to the same "//&
                 "restart file or each processor writes its own (numbered) restart file. "//&
                 "If false, a single restart file is generated combining output from all PEs.", &
                 default=.false.)

  if (present(restart_root)) then
    CS%restartfile = restart_root
    call log_param(param_file, mdl, "RESTARTFILE from argument", CS%restartfile)
  else
    call get_param(param_file, mdl, "RESTARTFILE", CS%restartfile, &
                 "The name-root of the restart file.", default="MOM.res")
  endif
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 "If true, require the restart checksums to match and error out otherwise. "//&
                 "Users may want to avoid this comparison if for example the restarts are "//&
                 "made from a run with a different mask_table than the current run, "//&
                 "in which case the checksums will not match and cause crash.",&
                 default=.true.)
  call get_param(param_file, mdl, "RESTART_SYMMETRIC_CHECKSUMS", CS%symmetric_checksums, &
                 "If true, do the restart checksums on all the edge points for a non-reentrant "//&
                 "grid.  This requires that SYMMETRIC_MEMORY_ is defined at compile time.", &
                 default=.false.)
  call get_param(param_file, mdl, "RESTART_UNSIGNED_ZEROS", CS%unsigned_zeros, &
                 "If true, convert any negative zeros that would be written to the restart file "//&
                 "into ordinary unsigned zeros.  This does not change answers, but it can be "//&
                 "helpful in comparing restart files after grid rotation, for example.", &
                 default=.false.)
  call get_param(param_file, mdl, "REENTRANT_X", CS%reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "REENTRANT_Y", CS%reentrant_y, &
                 "If true, the domain is meridionally reentrant.", default=.false., do_not_log=.true.)

  ! Maybe not the best place to do this?
  call get_param(param_file, mdl, "ROTATE_INDEX", rotate_index, &
      default=.false., do_not_log=.true.)

  CS%turns = 0
  if (rotate_index) then
    call get_param(param_file, mdl, "INDEX_TURNS", CS%turns, &
        default=1, do_not_log=.true.)
  endif

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%restart_obsolete(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

  CS%locked = .false.

end procedure restart_init
module procedure lock_check
  character(len=256) :: var_name  ! A variable name.
  if (CS%locked) then
    if (present(var_desc)) then
      call query_vardesc(var_desc, name=var_name)
      call MOM_error(FATAL, "Attempted to register "//trim(var_name)//" but the restart registry is locked.")
    elseif (present(name)) then
      call MOM_error(FATAL, "Attempted to register "//trim(name)//" but the restart registry is locked.")
    else
      call MOM_error(FATAL, "Attempted to register a variable but the restart registry is locked.")
    endif
  endif

end procedure lock_check
module procedure restart_registry_lock
  CS%locked = .true.
  if (present(unlocked)) CS%locked = .not.unlocked
end procedure restart_registry_lock
module procedure restart_init_end
  if (associated(CS)) then
    CS%locked = .true.

    if (CS%novars == 0) call restart_end(CS)
  endif

end procedure restart_init_end
module procedure restart_end
  if (associated(CS%restart_field)) deallocate(CS%restart_field)
  if (associated(CS%restart_obsolete)) deallocate(CS%restart_obsolete)
  if (associated(CS%var_ptr0d)) deallocate(CS%var_ptr0d)
  if (associated(CS%var_ptr1d)) deallocate(CS%var_ptr1d)
  if (associated(CS%var_ptr2d)) deallocate(CS%var_ptr2d)
  if (associated(CS%var_ptr3d)) deallocate(CS%var_ptr3d)
  if (associated(CS%var_ptr4d)) deallocate(CS%var_ptr4d)
  deallocate(CS)

end procedure restart_end
module procedure restart_error
  character(len=16)  :: num  ! String for error messages
  if (CS%novars > CS%max_fields) then
    write(num,'(I0)') CS%novars
    call MOM_error(FATAL,"MOM_restart: Too many fields registered for " // &
           "restart.  Set MAX_FIELDS to be at least "//trim(num)//" in the MOM input file.")
  else
    call MOM_error(FATAL,"MOM_restart: Unspecified fatal error.")
  endif
end procedure restart_error
module procedure get_checksum_loop_ranges
  isL = G%isc-G%isd+1
  ieL = G%iec-G%isd+1
  jsL = G%jsc-G%jsd+1
  jeL = G%jec-G%jsd+1

  ! Expand range east or south for symmetric arrays
  if (CS%symmetric_checksums) then
    if (.not.G%symmetric) call MOM_error(FATAL, &
        "Setting SYMMETRIC_RESTART_CHECKSUMS to true only works with symmetric memory allocation, "//&
        "which is specified at compile time by defining the cpp macro SYMMETRIC_MEMORY_.")

    if (((pos == EAST_FACE) .or. (pos == CORNER)) .and. (.not.CS%reentrant_x)) then ! For u-, q-points only
      if (G%isc+G%idg_offset == 1) isL = isL - 1 ! Include western edge in checksums only for western PEs
    endif
    if (((pos == NORTH_FACE) .or. (pos == CORNER)) .and. (.not.CS%reentrant_y)) then ! For v-, q-points only
      if (G%jsc+G%jdg_offset == 1) jsL = jsL - 1 ! Include southern edge in checksums only for southern PEs
    endif
  endif

end procedure get_checksum_loop_ranges
module procedure get_variable_byte_size
  integer :: var_periods  ! The number of entries in a time-periodic axis
  character(len=8) :: t_grid_read, t_grid_tmp ! Modified versions of t_grid
  if (pos == 0) then
    var_sz = 8
  else ! This may be an overestimate, as it is based on symmetric-memory corner points.
    var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
  endif

  select case (trim(z_grid))
    case ('L') ; var_sz = var_sz * num_z
    case ('i') ; var_sz = var_sz * (num_z+1)
  end select

  t_grid_tmp = adjustl(t_grid)
  if (t_grid_tmp(1:1) == 'p') then
    if (len_trim(t_grid_tmp(2:8)) > 0) then
      var_periods = -1
      t_grid_read = adjustl(t_grid_tmp(2:8))
      read(t_grid_read,*) var_periods
      if (var_periods > 1) var_sz = var_sz * var_periods
    endif
  endif

end procedure get_variable_byte_size
end submodule MOM_restart_s
