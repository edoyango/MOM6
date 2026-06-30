submodule (numerical_testing_type) numerical_testing_type_s
  implicit none
contains
module procedure test
  logical :: ignore_this_fail
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore

  this%num_tests_checked = this%num_tests_checked + 1
  if (state) then
    if (.not. ignore_this_fail) then
      this%state = .true.
      this%num_tests_failed = this%num_tests_failed + 1
      if (this%num_tests_failed<=100) this%ifailed(this%num_tests_failed) = this%num_tests_checked
      if (this%num_tests_failed == 1) this%label_first_fail = label
      write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" FAILED!'
      write(this%stderr, '(2x,3a)') 'Test "',trim(label),'" FAILED!'
    else
      write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" IGNORED!'
      write(this%stderr, '(2x,3a)') 'Test "',trim(label),'" IGNORED!'
    endif
  elseif (this%verbose) then
    write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" passed'
  endif
  if (this%stop_instantly .and. this%state .and. .not. ignore_this_fail) stop 1
end procedure test
module procedure set
  if (present(verbose)) then
    this%verbose = verbose
  endif
  if (present(stdout)) then
    this%stdout = stdout
  endif
  if (present(stderr)) then
    this%stderr = stderr
  endif
  if (present(stop_instantly)) then
    this%stop_instantly = stop_instantly
  endif
  if (present(ignore_fail)) then
    this%ignore_fail = ignore_fail
  endif
end procedure set
module procedure summarize
  integer :: i
  if (this%state) then
    write(this%stdout,'(a," : ",a,", ",i4," failed of ",i4," tested")') &
         'FAIL', trim(label), this%num_tests_failed, this%num_tests_checked
    write(this%stdout,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,min(100,this%num_tests_failed))
    write(this%stdout,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,min(100,this%num_tests_failed))
    write(this%stderr,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a," : ",a)') trim(label),'FAILED'
  else
    write(this%stdout,'(a," : ",a,", all ",i4," tests passed")') &
         'Pass', trim(label), this%num_tests_checked
  endif
  summarize = this%state
end procedure summarize
module procedure real_scalar
  logical :: this_test, ignore_this_fail
  real :: tolerance, err ! Tolerance and error [A]
  tolerance = 0.0
  if (present(tol)) tolerance = tol
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  if (present(robits)) tolerance = abs(u_true) * float(robits) * epsilon(err)
  if (abs(u_test - u_true) > tolerance) this_test = .true.

  if (this_test) then
    if (ignore_this_fail) then
      if (this%verbose) then
        write(this%stdout,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- IGNORING"
        write(this%stderr,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- IGNORING"
      endif
      this_test = .false.
    else
      write(this%stdout,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- WRONG"
      write(this%stderr,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- WRONG"
    endif
  elseif (this%verbose) then
    write(this%stdout,'(2(a,1p1e24.16,1x),a)') "Calculated value =",u_test,"Correct value =",u_true,label
  endif

  call this%test( this_test, label, ignore=ignore_this_fail ) ! Updates state and counters in this
end procedure real_scalar
module procedure real_arr
  integer :: k
  logical :: this_test, ignore_this_fail
  real :: tolerance, err ! Tolerance and error [A]
  tolerance = 0.0
  if (present(tol)) tolerance = tol
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (present(robits)) tolerance = abs(u_true(k)) * float(robits) * epsilon(err)
    if (abs(u_test(k) - u_true(k)) > tolerance) this_test = .true.
  enddo

  ! If either being verbose, or an error was measured then display results
  if (this_test .or. this%verbose) then
    write(this%stdout,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    if (this_test) write(this%stderr,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    do k = 1, n
      if (present(robits)) tolerance = abs(u_true(k)) * float(robits) * epsilon(err)
      err = u_test(k) - u_true(k)
      if ( ( abs(err) > tolerance .and. ignore_this_fail ) .or. &
           ( abs(err) > 0. .and. abs(err) <= tolerance ) ) then
        write(this%stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- IGNORING'
      elseif (abs(err) > tolerance) then
        write(this%stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- WRONG'
        write(this%stderr,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- WRONG'
      else
        write(this%stdout,'(i4,1p2e24.16)') k, u_test(k), u_true(k)
      endif
    enddo
  endif

  call this%test( this_test, label, ignore=ignore_this_fail ) ! Updates state and counters in this
end procedure real_arr
module procedure int_arr
  integer :: k
  logical :: this_test, ignore_this_fail
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (i_test(k) /= i_true(k)) this_test = .true.
  enddo

  if (this%verbose) then
    write(this%stdout,'(a14," : calculated =",30i3)') label, i_test
    write(this%stdout,'(14x,"      correct =",30i3)') i_true
    if (this_test) then
      if (ignore_this_fail) then
        write(this%stdout,'(3x,a,8x,"error =",30i3)') 'IGNORE --->', i_test(:) - i_true(:)
      else
        write(this%stdout,'(3x,a,8x,"error =",30i3)') ' FAIL  --->', i_test(:) - i_true(:)
      endif
    endif
  endif

  if (ignore_this_fail) this_test = .false.

  if (this_test) then
    write(this%stderr,'(a14," : calculated =",30i3)') label, i_test
    write(this%stderr,'(14x,"      correct =",30i3)') i_true
    write(this%stderr,'("   FAIL --->        error =",30i3)') i_test(:) - i_true(:)
  endif

  call this%test( this_test, label ) ! Updates state and counters in this
end procedure int_arr
module procedure numerical_testing_type_unit_tests
  type(testing) :: tester ! An instance to record tests
  type(testing) :: test ! The instance used for testing (is mutable)
  logical :: tmpflag ! Temporary for return flags
  numerical_testing_type_unit_tests = .false. ! Assume all is well at the outset
  if (verbose) write(test%stdout,*) "  ===== testing_type: numerical_testing_type_unit_tests ====="
  call tester%set( verbose=verbose ) ! Sets the verbosity flag in tester

  call test%set( verbose=verbose ) ! Sets the verbosity flag in test
  call test%set( stderr=6 ) ! Sets stderr (redirect errors for "test" since they are not real)
  call test%set( stdout=6 ) ! Sets stdout
  call test%set( stop_instantly=.false. ) ! Sets stop_instantly
  call test%set( ignore_fail=.false. ) ! Sets ignore_fail

  ! Check that %summary() reports nothing when %state is unset
  ! (note this has to be confirmed visually since everything is in stdout)
  tmpflag = test%summarize("Summary is for a passing state")
  call tester%test(tmpflag, "test%summarize() with no fails")

  ! Check that %test(.false.,...) leaves %state unchanged
  call test%test( .false., "test(F) should pass" )
  call tester%test(test%state, "test%test(F)")

  ! Check that %test(.true.,...,ignore=.true.) leaves %state unchanged
  call test%test( .true., "test(T) should fail but be ignored", ignore=.true. )
  call tester%test(test%state, "test%test(T,ignore)")

  ! Check that %test(.true.,...) sets %state
  call test%test( .true., "test(T) should fail" )
  call tester%test(.not. test%state, "test%test(T,ignore)")
  test%state = .false. ! reset

  ! Check that %real_scalar(a,a,...) leaves %state unchanged
  call test%real_scalar(1., 1., "real_scalar(s,s) should pass", robits=0, tol=0.)
  call tester%test(test%state, "test%real_scalar(s,s)")

  ! Check that %real_scalar(a,b,...,ignore=.true.) leaves %state unchanged
  call test%real_scalar(1., 2., "real_scalar(s,t) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%real_scalar(s,t,ignore)")

  ! Check that %real_scalar(a,a,...) sets %state
  call test%real_scalar(1., 2., "s != t should fail")
  call tester%test(.not. test%state, "test%real_scalar(s,t)")
  test%state = .false. ! reset

  ! Check that %real_arr(a,a,...) leaves %state unchanged
  call test%real_arr(2, (/1.,2./), (/1.,2./), "real_arr(a,a) should pass", robits=0, tol=0.)
  call tester%test(test%state, "test%real_arr(a,a)")

  ! Check that %real_arr(a,b,...,ignore=.true.) leaves %state unchanged
  call test%real_arr(2, (/1.,2./), (/3.,4./), "real_arr(a,b) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%real_arr(a,b,ignore)")

  ! Check that %real_arr(a,b,...) sets %state
  call test%real_arr(2, (/1.,2./), (/3.,4./), "real(a,b) should fail")
  call tester%test(.not. test%state, "test%real_arr(a,b)")
  test%state = .false. ! reset

  ! Check that %int_arr(a,a,...) leaves %state unchanged
  call test%int_arr(2, (/1,2/), (/1,2/), "int_arr(i,i) should pass")
  call tester%test(test%state, "test%int_arr(i,i)")

  ! Check that %int_arr(a,b,...,ignore=.true.) leaves %state unchanged
  call test%int_arr(2, (/1,2/), (/3,4/), "int_arr(i,j) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%int_arr(i,j,ignore)")

  ! Check that %int_arr(a,b,...) sets %state
  call test%int_arr(2, (/1,2/), (/3,4/), "int(arr(i,j) should fail")
  call tester%test(.not. test%state, "test%int_arr(i,j)")
  test%state = .false. ! reset

  ! Check that %summary() reports nothing when %state is set
  ! (note this has to be confirmed visually since everything is in stdout)
  test%state = .true. ! reset to fail for testing %summary()
  tmpflag = test%summarize("This summary should report 4 fails")
  call tester%test(.not. tmpflag, "test%summarize() with fails")

  numerical_testing_type_unit_tests = tester%summarize("numerical_testing_type_unit_tests")

end procedure numerical_testing_type_unit_tests
end submodule numerical_testing_type_s
