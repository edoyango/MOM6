submodule (MOM_string_functions) MOM_string_functions_s
  implicit none
contains
module procedure lowercase
  integer, parameter :: co=iachar('a')-iachar('A') ! case offset
  integer :: k
  lowercase = input_string
  do k=1, len_trim(input_string)
    if (lowercase(k:k) >= 'A' .and. lowercase(k:k) <= 'Z') &
        lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
  enddo
end procedure lowercase
module procedure uppercase
  integer, parameter :: co=iachar('A')-iachar('a') ! case offset
  integer :: k
  uppercase = input_string
  do k=1, len_trim(input_string)
    if (uppercase(k:k) >= 'a' .and. uppercase(k:k) <= 'z') &
        uppercase(k:k) = achar(ichar(uppercase(k:k))+co)
  enddo
end procedure uppercase
module procedure left_int
  character(len=19) :: tmp
  write(tmp(1:19),'(I19)') i
  write(left_int(1:19),'(A)') adjustl(tmp)
end procedure left_int
module procedure left_ints
  character(len=1320) :: tmp
  integer :: j
  write(left_ints(1:1320),'(A)') trim(left_int(i(1)))
  if (size(i)>1) then
    do j=2,size(i)
      tmp=left_ints
      write(left_ints(1:1320),'(A,", ",A)') trim(tmp),trim(left_int(i(j)))
    enddo
  endif
end procedure left_ints
module procedure left_real
  integer :: l, ind
  if ((abs(val) < 1.0e4) .and. (abs(val) >= 1.0e-3)) then
    write(left_real, '(F30.11)') val
    if (.not.isFormattedFloatEqualTo(left_real,val)) then
      write(left_real, '(F30.12)') val
      if (.not.isFormattedFloatEqualTo(left_real,val)) then
        write(left_real, '(F30.13)') val
        if (.not.isFormattedFloatEqualTo(left_real,val)) then
          write(left_real, '(F30.14)') val
          if (.not.isFormattedFloatEqualTo(left_real,val)) then
            write(left_real, '(F30.15)') val
            if (.not.isFormattedFloatEqualTo(left_real,val)) then
              write(left_real, '(F30.16)') val
            endif
          endif
        endif
      endif
    endif
    do
      l = len_trim(left_real)
      if ((l<2) .or. (left_real(l-1:l) == ".0") .or. &
          (left_real(l:l) /= "0")) exit
      left_real(l:l) = " "
    enddo
  elseif (val == 0.) then
    left_real = "0.0"
  else
    if ((abs(val) <= 1.0e-100) .or. (abs(val) >= 1.0e100)) then
      write(left_real(1:32), '(ES24.14E3)') val
      if (.not.isFormattedFloatEqualTo(left_real,val)) &
        write(left_real(1:32), '(ES24.15E3)') val
    else
      write(left_real(1:32), '(ES23.14)') val
      if (.not.isFormattedFloatEqualTo(left_real,val)) &
        write(left_real(1:32), '(ES23.15)') val
    endif
    do
      ind = index(left_real,"0E")
      if (ind == 0) exit
      if (left_real(ind-1:ind-1) == ".") exit
      left_real = left_real(1:ind-1)//left_real(ind+1:)
    enddo
  endif
  left_real = adjustl(left_real)
end procedure left_real
module procedure left_reals
  integer :: j, n, ns
  logical :: doWrite
  character(len=10) :: separator
  n=1 ; doWrite=.true. ; left_reals=''
  if (present(sep)) then
    separator=sep ; ns=len(sep)
  else
    separator=', ' ; ns=2
  endif
  do j=1,size(r)
    doWrite=.true.
    if (j<size(r)) then
      if (r(j)==r(j+1)) then
        n=n+1
        doWrite=.false.
      endif
    endif
    if (doWrite) then
      if (len(left_reals)>0) then ! Write separator if a number has already been written
        left_reals = left_reals // separator(1:ns)
      endif
      if (n>1) then
        left_reals = left_reals // trim(left_int(n)) // "*" // trim(left_real(r(j)))
      else
        left_reals = left_reals // trim(left_real(r(j)))
      endif
      n=1
    endif
  enddo
end procedure left_reals
module procedure isFormattedFloatEqualTo
  real :: scannedVal ! The value extraced from str, in arbitrary units [A]
  isFormattedFloatEqualTo=.false.
  read(str(1:),*,err=987) scannedVal
  if (scannedVal == val) isFormattedFloatEqualTo=.true.
 987 return
end procedure isFormattedFloatEqualTo
module procedure extractWord
  extractWord = extract_word(string, ' ,', n)

end procedure extractWord
module procedure extract_word
  integer :: ns, i, b, e, nw
  logical :: lastCharIsSeperator
  extract_word = ''
  lastCharIsSeperator = .true.
  ns = len_trim(string)
  i = 0 ; b=0 ; e=0 ; nw=0
  do while (i<ns)
    i = i+1
    if (lastCharIsSeperator) then ! search for end of word
      if (verify(string(i:i),separators)==0) then
        continue ! Multiple separators
      else
        lastCharIsSeperator = .false. ! character is beginning of word
        b = i
        continue
      endif
    else ! continue search for end of word
      if (verify(string(i:i),separators)==0) then
        lastCharIsSeperator = .true.
        e = i-1 ! Previous character is end of word
        nw = nw+1
        if (nw==n) then
          extract_word = trim(string(b:e))
          return
        endif
      endif
    endif
  enddo
  if (b<=ns .and. nw==n-1) extract_word = trim(string(b:ns))
end procedure extract_word
module procedure extract_integer
  character(len=20) :: word
  word = extract_word(string, separators, n)

  if (len_trim(word)>0) then
    read(word(1:len_trim(word)),*) extract_integer
  else
    if (present(missing_value)) then
      extract_integer = missing_value
    else
      extract_integer = 0
    endif
  endif

end procedure extract_integer
module procedure extract_real
  character(len=20) :: word
  word = extract_word(string, separators, n)

  if (len_trim(word)>0) then
    read(word(1:len_trim(word)),*) extract_real
  else
    if (present(missing_value)) then
      extract_real = missing_value
    else
      extract_real = 0
    endif
  endif

end procedure extract_real
module procedure remove_spaces
  integer :: ns, i, o
  logical :: lastCharIsSeperator
  lastCharIsSeperator = .true.
  ns = len_trim(string)
  i = 0 ; o = 0
  do while (i<ns)
    i = i+1
    if (string(i:i) /= ' ') then ! Copy character to output string
      o = o + 1
      remove_spaces(o:o) = string(i:i)
    endif
  enddo
  do i = o+1, 120
    remove_spaces(i:i) = ' ' ! Wipe any non-empty characters
  enddo
  remove_spaces = trim(remove_spaces)
end procedure remove_spaces
module procedure string_functions_unit_tests
  integer :: i(5) = (/ -1, 1, 3, 3, 0 /)
  real :: r(8) = (/ 0., 1., -2., 1.3, 3.E-11, 3.E-11, 3.E-11, -5.1E12 /)
  logical :: fail, v
  fail = .false.
  v = verbose
  write(stdout,*) '==== MOM_string_functions: string_functions_unit_tests ==='
  fail = fail .or. localTestS(v,left_int(-1),'-1')
  fail = fail .or. localTestS(v,left_ints(i(:)),'-1, 1, 3, 3, 0')
  fail = fail .or. localTestS(v,left_real(0.),'0.0')
  fail = fail .or. localTestS(v,left_reals(r(:)),'0.0, 1.0, -2.0, 1.3, 3*3.0E-11, -5.1E+12')
  fail = fail .or. localTestS(v,left_reals(r(:),sep=' '),'0.0 1.0 -2.0 1.3 3*3.0E-11 -5.1E+12')
  fail = fail .or. localTestS(v,left_reals(r(:),sep=','),'0.0,1.0,-2.0,1.3,3*3.0E-11,-5.1E+12')
  fail = fail .or. localTestS(v,ints_to_string(i(:),5),'_-0001_0001_0003_0003_0000')
  fail = fail .or. localTestS(v,ints_to_string(i(2:),2),'_0001_0003')
  fail = fail .or. localTestS(v,ints_to_string(i(:)),'_-0001_0001_0003')
  fail = fail .or. localTestS(v,trim_trailing_commas("One, Two, Three, "), "One, Two, Three")
  fail = fail .or. localTestS(v,extractWord("One Two,Three",1),"One")
  fail = fail .or. localTestS(v,extractWord("One Two,Three",2),"Two")
  fail = fail .or. localTestS(v,extractWord("One Two,Three",3),"Three")
  fail = fail .or. localTestS(v,extractWord("One Two,  Three",3),"Three")
  fail = fail .or. localTestS(v,extractWord(" One Two,Three",1),"One")
  fail = fail .or. localTestS(v,extract_word("One,Two,Three",",",3),"Three")
  fail = fail .or. localTestS(v,extract_word("One,Two,Three",",",4),"")
  fail = fail .or. localTestS(v,remove_spaces("1 2 3"),"123")
  fail = fail .or. localTestS(v,remove_spaces(" 1 2 3"),"123")
  fail = fail .or. localTestS(v,remove_spaces("1 2 3 "),"123")
  fail = fail .or. localTestS(v,remove_spaces("123"),"123")
  fail = fail .or. localTestS(v,remove_spaces(" "),"")
  fail = fail .or. localTestS(v,remove_spaces(""),"")
  fail = fail .or. localTestI(v,extract_integer("1","",1),1)
  fail = fail .or. localTestI(v,extract_integer("1,2,3",",",1),1)
  fail = fail .or. localTestI(v,extract_integer("1,2",",",2),2)
  fail = fail .or. localTestI(v,extract_integer("1,2",",",3),0)
  fail = fail .or. localTestI(v,extract_integer("1,2",",",4,4),4)
  fail = fail .or. localTestR(v,extract_real("1.","",1),1.)
  fail = fail .or. localTestR(v,extract_real("1.,2.,3.",",",1),1.)
  fail = fail .or. localTestR(v,extract_real("1.,2.",",",2),2.)
  fail = fail .or. localTestR(v,extract_real("1.,2.",",",3),0.)
  fail = fail .or. localTestR(v,extract_real("1.,2.",",",4,4.),4.)
  if (.not. fail) write(stdout,*) 'Pass'
  string_functions_unit_tests = fail
end procedure string_functions_unit_tests
module procedure localTestS
  localTestS=.false.
  if (trim(str1)/=trim(str2)) localTestS=.true.
  if (localTestS .or. verbose) then
    write(stdout,*) '>'//trim(str1)//'<'
    if (localTestS) then
      write(stdout,*) trim(str1),':',trim(str2), '<-- FAIL'
      write(stderr,*) trim(str1),':',trim(str2), '<-- FAIL'
    endif
  endif
end procedure localTestS
module procedure localTestI
  localTestI=.false.
  if (i1/=i2) localTestI=.true.
  if (localTestI .or. verbose) then
    write(stdout,*) i1,i2
    if (localTestI) then
      write(stdout,*) i1,'!=',i2, '<-- FAIL'
      write(stderr,*) i1,'!=',i2, '<-- FAIL'
    endif
  endif
end procedure localTestI
module procedure localTestR
  localTestR=.false.
  if (r1/=r2) localTestR=.true.
  if (localTestR .or. verbose) then
    write(stdout,*) r1,r2
    if (localTestR) then
      write(stdout,*) r1,'!=',r2, '<-- FAIL'
      write(stderr,*) r1,'!=',r2, '<-- FAIL'
    endif
  endif
end procedure localTestR
module procedure slasher
  if (len_trim(dir) == 0) then
    slasher = "./"
  elseif (dir(len_trim(dir):len_trim(dir)) == '/') then
    slasher = trim(dir)
  else
    slasher = trim(dir)//"/"
  endif
end procedure slasher
module procedure trim_trailing_commas
  out_str = trim(adjustl(in_str))
  if (len_trim(out_str) > 0) then
    if (out_str(len_trim(out_str):len_trim(out_str)) == ",") then
      out_str = out_str(1:len_trim(out_str) - 1)
    endif
    out_str = trim(out_str)
  endif

end procedure trim_trailing_commas
module procedure ints_to_string
  character(len=8) :: i2s_temp
  integer :: i, n_max
  n_max = 3
  if (present(n)) n_max = n

  i2s = ''
  do i=1,min(size(a), n_max)
    if (a(i) < 0) then
      write (i2s_temp, '(I5.4)') a(i)
    else
      write (i2s_temp, '(I4.4)') a(i)
    endif
    i2s = trim(i2s) //'_'// trim(i2s_temp)
  enddo
  i2s = adjustl(i2s)
end procedure ints_to_string
end submodule MOM_string_functions_s
