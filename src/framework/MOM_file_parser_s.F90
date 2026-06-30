submodule (MOM_file_parser) MOM_file_parser_s
  implicit none
contains
module procedure open_param_file
  logical :: file_exists, Netcdf_file, may_check, reopened_file
  integer :: ios, iounit, strlen, i
  character(len=240) :: doc_path
  character(len=5)  :: ensemble_suffix
  type(parameter_block), pointer :: block => NULL()
  may_check = .true. ; if (present(checkable)) may_check = checkable

  ! Check for non-blank filename
  strlen = len_trim(filename)
  if (strlen == 0) then
    call MOM_error(FATAL, "open_param_file: Input file has not been specified.")
  endif

  ! Check that this file has not already been opened
  if (CS%nfiles > 0) then
    reopened_file = .false.

    if (is_root_pe()) then
      inquire(file=trim(filename), number=iounit)
      if (iounit /= -1) then
        do i = 1, CS%nfiles
          if (CS%iounit(i) == iounit) then
            call assert(trim(CS%filename(1)) == trim(filename), &
                "open_param_file: internal inconsistency! "//trim(filename)// &
                " is registered as open but has the wrong unit number!")
            call MOM_error(WARNING, &
                "open_param_file: file "//trim(filename)// &
                " has already been opened. This should NOT happen!"// &
                " Did you specify the same file twice in a namelist?")
            reopened_file = .true.
          endif ! unit numbers
        enddo ! i
      endif
    endif

    if (any_across_PEs(reopened_file)) return
  endif

  ! Check that the file exists to readstdlog
  if (is_root_pe()) then
    inquire(file=trim(filename), exist=file_exists)
    if (.not.file_exists) call MOM_error(FATAL, &
        "open_param_file: Input file '"// trim(filename)//"' does not exist.")
  endif

  Netcdf_file = .false.
  if (strlen > 3) then
    if (filename(strlen-2:strlen) == ".nc") Netcdf_file = .true.
  endif

  if (Netcdf_file) &
    call MOM_error(FATAL,"open_param_file: NetCDF files are not yet supported.")

  if (is_root_pe()) then
    open(newunit=iounit, file=trim(filename), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)
    if (ios /= 0) call MOM_error(FATAL, "open_param_file: Error opening '"//trim(filename)//"'.")
  else
    iounit = 1
  endif

  ! Store/register the unit and details
  i = CS%nfiles + 1
  CS%nfiles = i
  CS%iounit(i) = iounit
  CS%filename(i) = filename
  CS%NetCDF_file(i) = Netcdf_file

  if (associated(CS%blockName)) deallocate(CS%blockName)
  allocate(block) ; block%name = '' ; CS%blockName => block

  call MOM_mesg("open_param_file: "// trim(filename)//" has been opened successfully.", 5)

  call populate_param_data(iounit, filename, CS%param_data(i))
  ! Increment the maximum line length, but always report values in blocks of 4 characters.
  CS%max_line_len = max(CS%max_line_len, 4 + 4*(max_input_line_length(CS, i) - 1) / 4)

  call read_param(CS,"SEND_LOG_TO_STDOUT",CS%log_to_stdout)
  call read_param(CS,"REPORT_UNUSED_PARAMS",CS%report_unused)
  call read_param(CS,"FATAL_UNUSED_PARAMS",CS%unused_params_fatal)
  CS%doc_file = "MOM_parameter_doc"
  if (present(ensemble_num)) then
    ! append instance suffix to doc_file
    write(ensemble_suffix,'(A,I0.4)') '_', ensemble_num
    CS%doc_file = trim(CS%doc_file)//ensemble_suffix
  endif
  if (present(component)) CS%doc_file = trim(component)//"_parameter_doc"
  call read_param(CS,"DOCUMENT_FILE", CS%doc_file)
  if (.not.may_check) then
    CS%report_unused = .false.
    CS%unused_params_fatal = .false.
  endif

  ! Open the log file.
  CS%stdlog = stdlog() ; CS%stdout = stdout()
  CS%log_open = (stdlog() > 0)

  doc_path = CS%doc_file
  if (len_trim(CS%doc_file) > 0) then
    CS%complete_doc = complete_doc_default
    call read_param(CS, "COMPLETE_DOCUMENTATION", CS%complete_doc)
    CS%minimal_doc = minimal_doc_default
    call read_param(CS, "MINIMAL_DOCUMENTATION", CS%minimal_doc)
    if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
      doc_path = trim(slasher(doc_file_dir))//trim(CS%doc_file)
    endif ; endif
  else
    CS%complete_doc = .false.
    CS%minimal_doc = .false.
  endif
  call doc_init(doc_path, CS%doc, minimal=CS%minimal_doc, complete=CS%complete_doc, &
                layout=CS%complete_doc, debugging=CS%complete_doc)

end procedure open_param_file
module procedure close_param_file
  logical :: all_default
  character(len=128) :: docfile_default
  character(len=40)  :: mdl   ! This module's name.
# include "version_variable.h"
  integer :: i, n, num_unused
  if (present(quiet_close)) then ; if (quiet_close) then
    do i = 1, CS%nfiles
      if (is_root_pe()) close(CS%iounit(i))
      call MOM_mesg("close_param_file: "// trim(CS%filename(i))// &
                    " has been closed successfully.", 5)
      CS%iounit(i) = -1
      CS%filename(i) = ''
      CS%NetCDF_file(i) = .false.
      do n=1,CS%param_data(i)%num_lines ; deallocate(CS%param_data(i)%fln(n)%line) ; enddo
      deallocate (CS%param_data(i)%fln)
      deallocate (CS%param_data(i)%line_used)
    enddo
    CS%log_open = .false.
    call doc_end(CS%doc)
    deallocate(CS%doc)
    return
  endif ; endif

  ! Log the parameters for the parser.
  docfile_default = "MOM_parameter_doc"
  if (present(component)) docfile_default = trim(component)//"_parameter_doc"

  all_default = (CS%log_to_stdout .eqv. log_to_stdout_default)
  all_default = all_default .and. (trim(CS%doc_file) == trim(docfile_default))
  if (len_trim(CS%doc_file) > 0) then
    all_default = all_default .and. (CS%complete_doc .eqv. complete_doc_default)
    all_default = all_default .and. (CS%minimal_doc .eqv. minimal_doc_default)
  endif

  mdl = "MOM_file_parser"
  call log_version(CS, mdl, version, "", debugging=.true., log_to_all=.true., all_default=all_default)
  call log_param(CS, mdl, "SEND_LOG_TO_STDOUT", CS%log_to_stdout, &
                 "If true, all log messages are also sent to stdout.", &
                 default=log_to_stdout_default)
  call log_param(CS, mdl, "REPORT_UNUSED_PARAMS", CS%report_unused, &
                 "If true, report any parameter lines that are not used "//&
                 "in the run.", default=report_unused_default, &
                 debuggingParam=.true.)
  call log_param(CS, mdl, "FATAL_UNUSED_PARAMS", CS%unused_params_fatal, &
                 "If true, kill the run if there are any unused "//&
                 "parameters.", default=unused_params_fatal_default, &
                 debuggingParam=.true.)
  call log_param(CS, mdl, "DOCUMENT_FILE", CS%doc_file, &
                 "The basename for files where run-time parameters, their "//&
                 "settings, units and defaults are documented. Blank will "//&
                 "disable all parameter documentation.", default=docfile_default)
  if (len_trim(CS%doc_file) > 0) then
    call log_param(CS, mdl, "COMPLETE_DOCUMENTATION",  CS%complete_doc, &
                  "If true, all run-time parameters are "//&
                  "documented in "//trim(CS%doc_file)//&
                  ".all .", default=complete_doc_default)
    call log_param(CS, mdl, "MINIMAL_DOCUMENTATION", CS%minimal_doc, &
                  "If true, non-default run-time parameters are "//&
                  "documented in "//trim(CS%doc_file)//&
                  ".short .", default=minimal_doc_default)
  endif

  num_unused = 0
  do i = 1, CS%nfiles
    if (is_root_pe() .and. (CS%report_unused .or. &
                            CS%unused_params_fatal)) then
      ! Check for unused lines.
      do n=1,CS%param_data(i)%num_lines
        if (.not.CS%param_data(i)%line_used(n)) then
          num_unused = num_unused + 1
          if (CS%report_unused) &
            call MOM_error(WARNING, "Unused line in "//trim(CS%filename(i))// &
                            " : "//trim(CS%param_data(i)%fln(n)%line))
        endif
      enddo
    endif

    if (is_root_pe()) close(CS%iounit(i))
    call MOM_mesg("close_param_file: "// trim(CS%filename(i))//" has been closed successfully.", 5)
    CS%iounit(i) = -1
    CS%filename(i) = ''
    CS%NetCDF_file(i) = .false.
    do n=1,CS%param_data(i)%num_lines ; deallocate(CS%param_data(i)%fln(n)%line) ; enddo
    deallocate (CS%param_data(i)%fln)
    deallocate (CS%param_data(i)%line_used)
  enddo
  deallocate(CS%blockName)

  if (is_root_pe() .and. (num_unused>0) .and. CS%unused_params_fatal) &
    call MOM_error(FATAL, "Run stopped because of unused parameter lines.")

  CS%log_open = .false.
  call doc_end(CS%doc)
  deallocate(CS%doc)
end procedure close_param_file
module procedure populate_param_data
  character(len=INPUT_STR_LENGTH) :: line
  character(len=1), allocatable, dimension(:) :: char_buf
  integer, allocatable, dimension(:) :: line_len ! The trimmed length of each processed input line
  integer :: n, num_lines, total_chars, ch, rsc, llen, int_buf(2)
  logical :: inMultiLineComment
  if (is_root_pe()) then
    ! rewind the parameter file
    rewind(iounit)

    ! count the number of valid entries in the parameter file
    num_lines = 0
    total_chars = 0
    inMultiLineComment = .false.
    do while(.true.)
      read(iounit, '(a)', end=8) line
      line = replaceTabs(line)
      if (inMultiLineComment) then
        if (closeMultiLineComment(line)) inMultiLineComment=.false.
      else
        if (lastNonCommentNonBlank(line)>0) then
          line = removeComments(line)
          line = simplifyWhiteSpace(line(:len_trim(line)))
          num_lines = num_lines + 1
          total_chars = total_chars + len_trim(line)
        endif
        if (openMultiLineComment(line)) inMultiLineComment=.true.
      endif
    enddo ! while (.true.)
 8  continue ! get here when read() reaches EOF

    if (inMultiLineComment .and. is_root_pe()) &
      call MOM_error(FATAL, 'MOM_file_parser : A C-style multi-line comment '// &
                      '(/* ... */) was not closed before the end of '//trim(filename))


    int_buf(1) = num_lines
    int_buf(2) = total_chars
  endif  ! (is_root_pe())

  ! Broadcast the number of valid entries in parameter file
  call broadcast(int_buf, 2, root_pe())
  num_lines = int_buf(1)
  total_chars = int_buf(2)

  ! Set up the space for storing the actual lines.
  param_data%num_lines = num_lines
  allocate (line_len(num_lines), source=0)
  allocate (char_buf(total_chars), source=" ")

  ! Read the actual lines.
  if (is_root_pe()) then
    ! rewind the parameter file
    rewind(iounit)

    ! Populate param_data%fln%line
    num_lines = 0
    rsc = 0
    do while(.true.)
      read(iounit, '(a)', end=18) line
      line = replaceTabs(line)
      if (inMultiLineComment) then
        if (closeMultiLineComment(line)) inMultiLineComment=.false.
      else
        if (lastNonCommentNonBlank(line)>0) then
          line = removeComments(line)
          if ((len_trim(line) > 1000) .and. is_root_PE()) then
            call MOM_error(WARNING, "MOM_file_parser: Consider using continuation to split up "//&
                                    "the excessivley long parameter input line "//trim(line))
          endif
          line = simplifyWhiteSpace(line(:len_trim(line)))
          num_lines = num_lines + 1
          llen = len_trim(line)
          line_len(num_lines) = llen
          do ch=1,llen ; char_buf(rsc+ch)(1:1) = line(ch:ch) ; enddo
          rsc = rsc + llen
        endif
        if (openMultiLineComment(line)) inMultiLineComment=.true.
      endif
    enddo ! while (.true.)
18  continue ! get here when read() reaches EOF

    call assert(num_lines == param_data%num_lines, &
        'MOM_file_parser: Found different number of valid lines on second ' &
        // 'reading of '//trim(filename))
  endif  ! (is_root_pe())

  ! Broadcast the populated arrays line_len and char_buf
  call broadcast(line_len, num_lines, root_pe())
  call broadcast(char_buf(1:total_chars), 1, root_pe())

  ! Allocate space to hold contents of the parameter file, including the lines in param_data%fln
  allocate(param_data%fln(num_lines))
  allocate(param_data%line_used(num_lines))
  param_data%line_used(:) = .false.
  ! Populate param_data%fln%line with the keyword lines from parameter file
  rsc = 0
  do n=1,num_lines
    line(1:INPUT_STR_LENGTH) = " "
    do ch=1,line_len(n) ; line(ch:ch) = char_buf(rsc+ch)(1:1) ; enddo
    param_data%fln(n)%line = trim(line)
    rsc = rsc + line_len(n)
  enddo

  deallocate(char_buf) ; deallocate(line_len)

end procedure populate_param_data
module procedure openMultiLineComment
  integer :: icom, last
  openMultiLineComment = .false.
  last = lastNonCommentIndex(string)+1
  icom = index(string(last:), "/*")
  if (icom > 0) then
    openMultiLineComment=.true.
    last = last+icom+1
  endif
  icom = index(string(last:), "*/") ; if (icom > 0) openMultiLineComment=.false.
end procedure openMultiLineComment
module procedure closeMultiLineComment
  closeMultiLineComment = .false.
  if (index(string, "*/")>0) closeMultiLineComment=.true.
end procedure closeMultiLineComment
module procedure lastNonCommentIndex
  integer :: icom, last
  last = len_trim(string)
  icom = index(string(:last), "!") ; if (icom > 0) last = icom-1 ! F90 style
  icom = index(string(:last), "//") ; if (icom > 0) last = icom-1 ! C++ style
  icom = index(string(:last), "/*") ; if (icom > 0) last = icom-1 ! C style
  lastNonCommentIndex = last
end procedure lastNonCommentIndex
module procedure lastNonCommentNonBlank
  lastNonCommentNonBlank = len_trim(string(:lastNonCommentIndex(string))) ! Ignore remaining trailing blanks
end procedure lastNonCommentNonBlank
module procedure replaceTabs
  integer :: i
  do i=1, len(string)
    if (string(i:i)==achar(9)) then
      replaceTabs(i:i)=" "
    else
      replaceTabs(i:i)=string(i:i)
    endif
  enddo
end procedure replaceTabs
module procedure removeComments
  integer :: last
  removeComments=repeat(" ",len(string))
  last = lastNonCommentNonBlank(string)
  removeComments(:last)=adjustl(string(:last)) ! Copy only the non-comment part of string
end procedure removeComments
module procedure simplifyWhiteSpace
  integer :: i, j
  logical :: nonBlank = .false., insideString = .false.
  character(len=1) :: quoteChar=" "
  nonBlank  = .false. ; insideString = .false. ! NOTE: For some reason this line is needed??
  i=0
  simplifyWhiteSpace=repeat(" ",len(string)+16)
  do j=1,len_trim(string)
    if (insideString) then ! Do not change formatting inside strings
      i=i+1
      simplifyWhiteSpace(i:i)=string(j:j)
      if (string(j:j)==quoteChar) insideString=.false. ! End of string
    else ! The following is outside of string delimiters
      if (string(j:j)==" " .or. string(j:j)==achar(9)) then ! Space or tab
        if (nonBlank) then ! Only copy a blank if the preceding character was non-blank
          i=i+1
          simplifyWhiteSpace(i:i)=" " ! Not string(j:j) so that tabs are replace by blanks
          nonBlank=.false.
        endif
      elseif (string(j:j)=='"' .or. string(j:j)=="'") then ! Start a sting
        i=i+1
        simplifyWhiteSpace(i:i)=string(j:j)
        insideString=.true.
        quoteChar=string(j:j) ! Keep copy of starting quote
        nonBlank=.true.       ! For exit from string
      elseif (string(j:j)=='=') then
        ! Insert spaces if this character is "=" so that line contains " = "
        if (nonBlank) then
          i=i+1
          simplifyWhiteSpace(i:i)=" "
        endif
        i=i+2
        simplifyWhiteSpace(i-1:i)=string(j:j)//" "
        nonBlank=.false.
      else ! All other characters
        i=i+1
        simplifyWhiteSpace(i:i)=string(j:j)
        nonBlank=.true.
      endif
    endif ! if (insideString)
  enddo ! j
  if (insideString) then ! A missing close quote should be flagged
    if (is_root_pe()) call MOM_error(FATAL, &
      "There is a mismatched quote in the parameter file line: "// &
      trim(string))
  endif
end procedure simplifyWhiteSpace
module procedure read_param_int
  character(len=CS%max_line_len) :: value_string(1)
  logical            :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,err = 1001) value
    if (present(set)) set = .true.
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_int: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_int: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
    if (present(set)) set = .false.
  endif
  return
 1001 call MOM_error(FATAL,'read_param_int: read error for integer variable '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end procedure read_param_int
module procedure read_param_int_array
  character(len=CS%max_line_len) :: value_string(1)
  logical            :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    if (present(set)) set = .true.
    read(value_string(1),*,end=991,err=1002) value
 991 return
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_int_array: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_int_array: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
    if (present(set)) set = .false.
  endif
  return
 1002 call MOM_error(FATAL,'read_param_int_array: read error for integer array '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end procedure read_param_int_array
module procedure read_param_real
  character(len=CS%max_line_len) :: value_string(1)
  logical            :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,err=1003) value
    if (present(scale)) value = scale*value
    if (present(set)) set = .true.
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_real: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_real: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
    if (present(set)) set = .false.
  endif
  return
 1003 call MOM_error(FATAL,'read_param_real: read error for real variable '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end procedure read_param_real
module procedure read_param_real_array
  character(len=CS%max_line_len) :: value_string(1)
  logical                        :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991,err=1004) value
991 continue
    if (present(scale)) value(:) = scale*value(:)
    if (present(set)) set = .true.
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_real_array: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_real_array: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
    if (present(set)) set = .false.
  endif
  return
 1004 call MOM_error(FATAL,'read_param_real_array: read error for real array '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end procedure read_param_real_array
module procedure read_param_char
  character(len=CS%max_line_len) :: value_string(1)
  logical            :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found) then
    value = trim(strip_quotes(value_string(1)))
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL, 'Unable to find variable '//trim(varname)//' in any input files.')
  endif ; endif

  if (present(set)) set = found

end procedure read_param_char
module procedure read_param_char_array
  character(len=CS%max_line_len) :: value_string(1), loc_string
  logical            :: found, defined
  integer            :: i, i_out
  call get_variable_line(CS, varname, found, defined, value_string)
  if (found) then
    loc_string = trim(value_string(1))
    i = index(loc_string,",")
    i_out = 1
    do while(i>0)
      value(i_out) = trim(strip_quotes(loc_string(:i-1)))
      i_out = i_out+1
      loc_string = trim(adjustl(loc_string(i+1:)))
      i = index(loc_string,",")
    enddo
    if (len_trim(loc_string)>0) then
      value(i_out) = trim(strip_quotes(adjustl(loc_string)))
      i_out = i_out+1
    endif
    do i=i_out,SIZE(value) ; value(i) = " " ; enddo
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL, 'Unable to find variable '//trim(varname)//' in any input files.')
  endif ; endif

  if (present(set)) set = found

end procedure read_param_char_array
module procedure read_param_logical
  character(len=CS%max_line_len) :: value_string(1)
  logical            :: found, defined
  call get_variable_line(CS, varname, found, defined, value_string, paramIsLogical=.true.)
  if (found) then
    value = defined
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL, 'Unable to find variable '//trim(varname)//' in any input files.')
  endif ; endif

  if (present(set)) set = found

end procedure read_param_logical
module procedure read_param_time
  character(len=CS%max_line_len) :: value_string(1)
  character(len=240) :: err_msg
  logical            :: found, defined
  real               :: real_time, time_unit
  integer            :: vals(7)
  if (present(date_format)) date_format = .false.

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    ! Determine whether value string should be parsed for a real number
    ! or a date, in either a string format or a comma-delimited list of values.
    if ((INDEX(value_string(1),'-') > 0) .and. &
        (INDEX(value_string(1),'-',back=.true.) > INDEX(value_string(1),'-'))) then
      ! There are two dashes, so this must be a date format.
      value = set_date(value_string(1), err_msg=err_msg)
      if (LEN_TRIM(err_msg) > 0) call MOM_error(FATAL,'read_param_time: '//&
          trim(err_msg)//' in integer list read error for time-type variable '//&
          trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      if (present(date_format)) date_format = .true.
    elseif (INDEX(value_string(1),',') > 0) then
      ! Initialize vals with an invalid date.
      vals(:) = (/ -999, -999, -999, 0, 0, 0, 0 /)
      read(value_string(1), *, end=995, err=1005) vals
      995 continue
      if ((vals(1) < 0) .or. (vals(2) < 0) .or. (vals(3) < 0)) &
        call MOM_error(FATAL,'read_param_time: integer list read error for time-type variable '//&
                       trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      value = set_date(vals(1), vals(2), vals(3), vals(4), vals(5), vals(6), &
                       vals(7), err_msg=err_msg)
      if (LEN_TRIM(err_msg) > 0) call MOM_error(FATAL,'read_param_time: '//&
          trim(err_msg)//' in integer list read error for time-type variable '//&
          trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      if (present(date_format)) date_format = .true.
    else
      time_unit = 1.0 ; if (present(timeunit)) time_unit = timeunit
      read( value_string(1), *) real_time
      value = real_to_time(real_time*time_unit)
    endif
    if (present(set)) set = .true.
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL, 'Unable to find variable '//trim(varname)//' in any input files.')
      else
        call MOM_error(FATAL, 'Variable '//trim(varname)//' found but not set in input files.')
      endif
    endif ; endif
    if (present(set)) set = .false.
  endif
  return

  1005 call MOM_error(FATAL, 'read_param_time: read error for time-type variable '//&
                             trim(varname)// ' parsing "'//trim(value_string(1))//'"')
end procedure read_param_time
module procedure strip_quotes
  integer :: i
  strip_quotes = val_str
  i = index(strip_quotes,ACHAR(34)) ! Double quote
  do while (i>0)
    if (i > 1) then ; strip_quotes = strip_quotes(:i-1)//strip_quotes(i+1:)
    else ; strip_quotes = strip_quotes(2:) ; endif
    i = index(strip_quotes,ACHAR(34)) ! Double quote
  enddo
  i = index(strip_quotes,ACHAR(39)) ! Single quote
  do while (i>0)
    if (i > 1) then ; strip_quotes = strip_quotes(:i-1)//strip_quotes(i+1:)
    else ; strip_quotes = strip_quotes(2:) ; endif
    i = index(strip_quotes,ACHAR(39)) ! Single quote
  enddo
end procedure strip_quotes
module procedure max_input_line_length
  character(len=FILENAME_LENGTH) :: filename
  character :: last_char
  integer :: ipf, ipf_s, ipf_e
  integer :: last, line_len, count, contBufSize
  logical :: continuedLine
  max_len = 0
  ipf_s = 1 ; ipf_e = CS%nfiles
  if (present(pf_num)) then
    if (pf_num > CS%nfiles) return
    ipf_s = pf_num ; ipf_e = pf_num
  endif

  paramfile_loop: do ipf = ipf_s, ipf_e
    filename = CS%filename(ipf)
    contBufSize = 0
    continuedLine = .false.

    ! Scan through each line of the file
    do count = 1, CS%param_data(ipf)%num_lines
      ! line = CS%param_data(ipf)%fln(count)%line
      last = len_trim(CS%param_data(ipf)%fln(count)%line)
      last_char = " "
      if (last > 0) last_char = CS%param_data(ipf)%fln(count)%line(last:last)
      ! Check if line ends in continuation character (either & or \)
      ! Note achar(92) is a backslash
      if (last_char == achar(92) .or. last_char == "&") then
        contBufSize = contBufSize + last - 1
        continuedLine = .true.
        if (count==CS%param_data(ipf)%num_lines .and. is_root_pe()) &
           call MOM_error(FATAL, "MOM_file_parser : the last line of the file ends in a"// &
                 " continuation character but there are no more lines to read. "// &
                 " Line: '"//trim(CS%param_data(ipf)%fln(count)%line(:last))//"'"// &
                 " in file "//trim(filename)//".")
        cycle ! cycle inorder to append the next line of the file
      elseif (continuedLine) then
        ! If we reached this point then this is the end of line continuation
        line_len = contBufSize + last
        contBufSize = 0
        continuedLine = .false.
      else  ! This is a simple line with no continuation.
        line_len = last
      endif
      max_len = max(max_len, line_len)
    enddo ! CS%param_data(ipf)%num_lines
  enddo paramfile_loop

end procedure max_input_line_length
module procedure get_variable_line
  character(len=CS%max_line_len) :: val_str, lname, origLine
  character(len=CS%max_line_len) :: line, continuationBuffer
  character(len=240) :: blockName
  character(len=FILENAME_LENGTH) :: filename
  integer            :: is, id, isd, isu, ise, iso, ipf
  integer            :: last, last1, ival, oval, max_vals, count, contBufSize
  character(len=52)  :: set
  logical            :: found_override, found_equals
  logical            :: found_define, found_undef
  logical            :: force_cycle, defined_in_line, continuedLine
  logical            :: variableKindIsLogical, valueIsSame
  logical            :: inWrongBlock, fullPathParameter
  logical, parameter :: requireNamedClose = .false.
  integer, parameter :: verbose = 1
  set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  continuationBuffer = repeat(" ", CS%max_line_len)
  contBufSize = 0

  variableKindIsLogical=.false.
  if (present(paramIsLogical)) variableKindIsLogical = paramIsLogical

  ! Find the first instance (if any) where the named variable is found, and
  ! return variables indicating whether this variable is defined and the string
  ! that contains the value of this variable.
  found = .false.
  oval = 0 ; ival = 0
  max_vals = SIZE(value_string)
  do is=1,max_vals ; value_string(is) = " " ; enddo

  paramfile_loop: do ipf = 1, CS%nfiles
    filename = CS%filename(ipf)
    continuedLine = .false.
    blockName = ''

    ! Scan through each line of the file
    do count = 1, CS%param_data(ipf)%num_lines
      line = CS%param_data(ipf)%fln(count)%line
      last = len_trim(line)

      last1 = max(1,last)
      ! Check if line ends in continuation character (either & or \)
      ! Note achar(92) is a backslash
      if (line(last1:last1) == achar(92).or.line(last1:last1) == "&") then
        continuationBuffer(contBufSize+1:contBufSize+len_trim(line))=line(:last-1)
        contBufSize=contBufSize + len_trim(line)-1
        continuedLine = .true.
        if (count==CS%param_data(ipf)%num_lines .and. is_root_pe()) &
           call MOM_error(FATAL, "MOM_file_parser : the last line"// &
                 " of the file ends in a continuation character but"// &
                 " there are no more lines to read. "// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        cycle ! cycle inorder to append the next line of the file
      elseif (continuedLine) then
        ! If we reached this point then this is the end of line continuation
        continuationBuffer(contBufSize+1:contBufSize+len_trim(line))=line(:last)
        line = continuationBuffer
        continuationBuffer=repeat(" ",CS%max_line_len) ! Clear for next use
        contBufSize = 0
        continuedLine = .false.
        last = len_trim(line)
      endif

      origLine = trim(line) ! Keep original for error messages

      ! Check for '#override' at start of line
      found_override = .false. ; found_define = .false. ; found_undef = .false.
      iso = index(line(:last), "#override " )! ; if (is > 0) found_override = .true.
      if (iso>1) call MOM_error(FATAL, "MOM_file_parser : #override was found "// &
                 " but was not the first keyword."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      if (iso==1) then
        found_override = .true.
        if (index(line(:last), "#override define ")==1) found_define = .true.
        if (index(line(:last), "#override undef ")==1) found_undef = .true.
        line = trim(adjustl(line(iso+10:last))) ; last = len_trim(line)
      endif

      ! Newer form of parameter block, block%, %block or block%param or
      iso=index(line(:last),'%')
      fullPathParameter = .false.
      if (iso==1) then ! % is first character means this is a close
        if (len_trim(blockName)==0 .and. is_root_pe()) call MOM_error(FATAL, &
            'get_variable_line: An extra close block was encountered. Line="'// &
            trim(line(:last))//'"' )
        if (last>1 .and. trim(blockName)/=trim(line(2:last)) .and. is_root_pe()) &
            call MOM_error(FATAL, 'get_variable_line: A named close for a parameter'// &
            ' block did not match the open block. Line="'//trim(line(:last))//'"' )
        if (last==1 .and. requireNamedClose) & ! line = '%' is a generic (unnamed) close
            call MOM_error(FATAL, 'get_variable_line: A named close for a parameter'// &
            ' block is required but found "%". Block="'//trim(blockName)//'"' )
        blockName = popBlockLevel(blockName)
        call flag_line_as_read(CS%param_data(ipf)%line_used,count)
      elseif (iso==last) then ! This is a new block if % is last character
        blockName = pushBlockLevel(blockName, line(:iso-1))
        call flag_line_as_read(CS%param_data(ipf)%line_used,count)
      else ! This is of the form block%parameter = ... (full path parameter)
        iso=index(line(:last),'%',.true.)
        ! Check that the parameter block names on the line matches the state set by the caller
        if (iso>0 .and. trim(CS%blockName%name)==trim(line(:iso-1))) then
          fullPathParameter = .true.
          line = trim(line(iso+1:last)) ! Strip away the block name for subsequent processing
          last = len_trim(line)
        endif
      endif

      ! We should only interpret this line if this block is the active block
      inWrongBlock = .false.
      if (len_trim(blockName)>0) then ! In a namelist block in file
        if (trim(CS%blockName%name)/=trim(blockName)) inWrongBlock = .true. ! Not in the required block
      endif
      if (len_trim(CS%blockName%name)>0) then ! In a namelist block in the model
        if (trim(CS%blockName%name)/=trim(blockName)) inWrongBlock = .true. ! Not in the required block
      endif

      if (inWrongBlock .and. .not. fullPathParameter) then
        if (index(" "//line(:last+1), " "//trim(varname)//" ")>0) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
               ' found outside of block '//trim(CS%blockName%name)//'%. Ignoring.')
        cycle
      endif

      ! Determine whether this line mentions the named parameter or not
      if (index(" "//line(:last)//" ", " "//trim(varname)//" ") == 0) cycle

      ! Detect keywords
      found_equals = .false.
      isd = index(line(:last), "define" )! ; if (isd > 0) found_define = .true.
      isu = index(line(:last), "undef" )! ; if (isu > 0) found_undef = .true.
      ise = index(line(:last), " = " ) ; if (ise > 1) found_equals = .true.
      if (index(line(:last), "#define ")==1) found_define = .true.
      if (index(line(:last), "#undef ")==1) found_undef = .true.

      ! Check for missing, mutually exclusive or incomplete keywords
      if (.not. (found_define .or. found_undef .or. found_equals)) then
        if (found_override) then
          call MOM_error(FATAL, "MOM_file_parser : override was found " // &
              " without a define or undef." // &
              " Line: '" // trim(line(:last)) // "'" // &
              " in file " // trim(filename) // ".")
        else
          call MOM_error(FATAL, "MOM_file_parser : the parameter name '" // &
              trim(varname) // "' was found without define or undef." // &
              " Line: '" // trim(line(:last)) // "'" // &
              " in file " // trim(filename) // ".")
        endif
      endif

      if (found_equals .and. (found_define .or. found_undef)) &
             call MOM_error(FATAL, &
               "MOM_file_parser : Both 'a=b' and 'undef/define' syntax occur."// &
               " Line: '"//trim(line(:last))//"'"//&
               " in file "//trim(filename)//".")

      ! Interpret the line and collect values, if any
      ! NOTE: At least one of these must be true
      if (found_define) then
        ! Move starting pointer to first letter of defined name.
        is = isd + 5 + scan(line(isd+6:last), set)

        id = scan(line(is:last), ' ')  ! Find space between name and value
        if ( id == 0 ) then
          ! There is no space so the name is simply being defined.
          lname = trim(line(is:last))
          if (trim(lname) /= trim(varname)) cycle
          val_str = " "
        else
          ! There is a string or number after the name.
          lname = trim(line(is:is+id-1))
          if (trim(lname) /= trim(varname)) cycle
          val_str = trim(adjustl(line(is+id:last)))
        endif
        found = .true. ; defined_in_line = .true.
      elseif (found_undef) then
        ! Move starting pointer to first letter of undefined name.
        is = isu + 4 + scan(line(isu+5:last), set)

        id = scan(line(is:last), ' ')  ! Find the first space after the name.
        if (id > 0) last = is + id - 1
        lname = trim(line(is:last))
        if (trim(lname) /= trim(varname)) cycle
        val_str = " "
        found = .true. ; defined_in_line = .false.
      elseif (found_equals) then
        ! Move starting pointer to first letter of defined name.
        is = scan(line(1:ise), set)
        lname = trim(line(is:ise-1))
        if (trim(lname) /= trim(varname)) cycle
        val_str = trim(adjustl(line(ise+3:last)))
        if (variableKindIsLogical) then ! Special handling for logicals
          read(val_str(:len_trim(val_str)),*) defined_in_line
        else
          defined_in_line = .true.
        endif
        found = .true.
      endif

      ! This line has now been used.
      call flag_line_as_read(CS%param_data(ipf)%line_used,count)

      ! Detect inconsistencies
      force_cycle = .false.
      valueIsSame = (trim(val_str) == trim(value_string(max_vals)))
      if (found_override .and. (oval >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line .neqv. defined) .or. .not. valueIsSame) then
            call MOM_error(FATAL,"MOM_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent overrides."// &
                     " Line A: '"//trim(value_string(max_vals))//"'"//&
                     " Line B: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                     " over-ridden more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (.not.found_override .and. (oval > 0)) then
        if (is_root_pe()) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                   " has already been over-ridden."// &
                   " Line: '"//trim(line(:last))//"'"//&
                   " in file "//trim(filename)//" is being ignored.")
        force_cycle = .true.
      endif
      if (.not.found_override .and. (ival >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line .neqv. defined) .or. .not. valueIsSame) then
            call MOM_error(FATAL,"MOM_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent definitions."// &
                     " Line A: '"//trim(value_string(max_vals))//"'"//&
                     " Line B: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                     " occurs more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (force_cycle) cycle

      ! Store new values
      if (found_override) then
        oval = oval + 1
        value_string(oval) = trim(val_str)
        defined = defined_in_line
        if (verbose > 0 .and. ival > 0 .and. is_root_pe() .and. &
            .not. overrideWarningHasBeenIssued(CS%chain, trim(varname)) ) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                 " over-ridden.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      else ! (.not. found_overide)
        ival = ival + 1
        value_string(ival) = trim(val_str)
        defined = defined_in_line

        if (verbose > 1 .and. is_root_pe()) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                 " set.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      endif

    enddo ! CS%param_data(ipf)%num_lines

    if (len_trim(blockName)>0 .and. is_root_pe()) call MOM_error(FATAL, &
      'A namelist/parameter block was not closed. Last open block appears '// &
      'to be "'//trim(blockName)//'".')

  enddo paramfile_loop

end procedure get_variable_line
module procedure flag_line_as_read
  line_used(count) = .true.
end procedure flag_line_as_read
module procedure overrideWarningHasBeenIssued
  type(link_parameter), pointer :: newLink => NULL(), this => NULL()
  overrideWarningHasBeenIssued = .false.
  this => chain
  do while( associated(this) )
    if (trim(varName) == trim(this%name)) then
      overrideWarningHasBeenIssued = .true.
      return
    endif
    this => this%next
  enddo
  allocate(newLink)
  newLink%name = trim(varName)
  newLink%hasIssuedOverrideWarning = .true.
  newLink%next => chain
  chain => newLink
end procedure overrideWarningHasBeenIssued
module procedure log_version_cs
  character(len=240) :: mesg
  mesg = trim(modulename)//": "//trim(version)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) call doc_module(CS%doc, modulename, desc, log_to_all, all_default, layout, debugging)

end procedure log_version_cs
module procedure log_version_plain
  character(len=240) :: mesg
  mesg = trim(modulename)//": "//trim(version)
  if (is_root_pe()) then
    write(stdlog(),'(a)') trim(mesg)
  endif

end procedure log_version_plain
module procedure log_param_int
  character(len=240) :: mesg, myunits
  write(mesg, '("  ",a," ",a,": ",a)') trim(modulename), trim(varname), trim(left_int(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits = " " ; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_int
module procedure log_param_int_array
  character(len=CS%max_line_len+120) :: mesg
  character(len=240) :: myunits
  write(mesg, '("  ",a," ",a,": ",A)') trim(modulename), trim(varname), trim(left_ints(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits = " " ; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, defaults, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_int_array
module procedure log_param_real
  real :: log_val ! The parameter value that is written out
  character(len=240) :: mesg, myunits
  log_val = value ; if (present(unscale)) log_val = unscale * value

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(left_real(log_val))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, log_val, default, &
                   debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_real
module procedure log_param_real_array
  real, dimension(size(value)) :: log_val ! The array of parameter values that is written out
  character(len=:), allocatable :: mesg
  character(len=240) :: myunits
  log_val(:) = value(:) ; if (present(unscale)) log_val(:) = unscale * value(:)

 !write(mesg, '("  ",a," ",a,": ",ES19.12,99(",",ES19.12))') &
 !write(mesg, '("  ",a," ",a,": ",G,99(",",G))') &
 !  trim(modulename), trim(varname), value
  mesg = "  " // trim(modulename) // " " // trim(varname) // ": " // trim(left_reals(log_val))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, log_val, default, defaults, &
                   debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_real_array
module procedure log_param_logical
  character(len=240) :: mesg, myunits
  if (value) then
    write(mesg, '("  ",a," ",a,": True")') trim(modulename), trim(varname)
  else
    write(mesg, '("  ",a," ",a,": False")') trim(modulename), trim(varname)
  endif
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits = "Boolean" ; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_logical
module procedure log_param_char
  character(len=:), allocatable :: mesg
  character(len=240) :: myunits
  mesg = "  " // trim(modulename) // " " // trim(varname) // ": " // trim(value)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits = " " ; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end procedure log_param_char
module procedure log_param_time
  real :: real_time, real_default
  logical :: use_timeunit, date_format
  character(len=240) :: mesg, myunits
  character(len=80) :: date_string, default_string
  integer :: days, secs, ticks
  use_timeunit = .false.
  date_format = .false. ; if (present(log_date)) date_format = log_date

  call get_time(value, secs, days, ticks)

  if (ticks == 0) then
    write(mesg, '("  ",a," ",a," (Time): ",i0,":",i0)') trim(modulename), &
       trim(varname), days, secs
  else
    write(mesg, '("  ",a," ",a," (Time): ",i0,":",i0,":",i0)') trim(modulename), &
       trim(varname), days, secs, ticks
  endif
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) then
    if (present(timeunit)) use_timeunit = (timeunit > 0.0)
    if (date_format) then
      myunits='[date]'

      date_string = convert_date_to_string(value)
      if (present(default)) then
        default_string = convert_date_to_string(default)
        call doc_param(CS%doc, varname, desc, myunits, date_string, &
                       default=default_string, layoutParam=layoutParam, &
                       debuggingParam=debuggingParam, like_default=like_default)
      else
        call doc_param(CS%doc, varname, desc, myunits, date_string, &
                       layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)
      endif
    elseif (use_timeunit) then
      if (present(units)) then
        write(myunits(1:240),'(A)') trim(units)
      else
        if (abs(timeunit-1.0) < 0.01) then ; myunits = "seconds"
        elseif (abs(timeunit-3600.0) < 1.0) then ; myunits = "hours"
        elseif (abs(timeunit-86400.0) < 1.0) then ; myunits = "days"
        elseif (abs(timeunit-3.1e7) < 1.0e6) then ; myunits = "years"
        else ; write(myunits,'(es8.2," sec")') timeunit ; endif
      endif
      real_time = (86400.0/timeunit)*days + secs/timeunit
      if (ticks > 0) real_time = real_time + &
                           real(ticks) / (timeunit*get_ticks_per_second())
      if (present(default)) then
        call get_time(default, secs, days, ticks)
        real_default = (86400.0/timeunit)*days + secs/timeunit
        if (ticks > 0) real_default = real_default + &
                           real(ticks) / (timeunit*get_ticks_per_second())
        call doc_param(CS%doc, varname, desc, myunits, real_time, real_default, like_default=like_default)
      else
        call doc_param(CS%doc, varname, desc, myunits, real_time, like_default=like_default)
      endif
    else
      call doc_param(CS%doc, varname, desc, value, default, units=units, like_default=like_default)
    endif
  endif

end procedure log_param_time
module procedure convert_date_to_string
  character(len=40) :: sub_string
  real    :: real_secs
  integer :: yrs, mons, days, hours, mins, secs, ticks, ticks_per_sec
  call get_date(date, yrs, mons, days, hours, mins, secs, ticks)
  write (date_string, '(i8.4)') yrs
  write (sub_string, '("-", i2.2, "-", I2.2, " ", i2.2, ":", i2.2, ":")') &
         mons, days, hours, mins
  date_string = trim(adjustl(date_string)) // trim(sub_string)
  if (ticks > 0) then
    ticks_per_sec = get_ticks_per_second()
    real_secs = secs + ticks/ticks_per_sec
    if (ticks_per_sec <= 100) then
      write (sub_string, '(F7.3)') real_secs
    else
      write (sub_string, '(F10.6)') real_secs
    endif
  else
    write (sub_string, '(i2.2)') secs
  endif
  date_string = trim(date_string) // trim(adjustl(sub_string))

end procedure convert_date_to_string
module procedure get_param_int
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  integer :: new_name_value  ! The value that is set when the standard name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value = value
      call read_param_int(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_int(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = (value == new_name_value)
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_int(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_int(CS, modulename, varname, value, desc, units, &
                       default, layoutParam, debuggingParam)
  endif

end procedure get_param_int
module procedure get_param_int_array
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  integer :: new_name_value(size(value))  ! The values that are set when the old name is used.
  integer :: m
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (present(defaults)) then
    if (present(default)) call MOM_error(FATAL, &
          "get_param_int_array: Only one of default and defaults can be specified at a time.")
    if (size(defaults) /= size(value)) call MOM_error(FATAL, &
          "get_param_int_array: The size of defaults and value are not the same.")
  endif

  if (do_read) then
    if (present(default)) value(:) = default
    if (present(defaults)) value(:) = defaults(:)

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value(:) = value(:)
      call read_param_int_array(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_int_array(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = .true.
        do m=1,size(value) ; if (value(m) /= new_name_value(m)) same_value = .false. ; enddo
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_int_array(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_int_array(CS, modulename, varname, value, desc, units, &
                             default, defaults, layoutParam, debuggingParam)
  endif

end procedure get_param_int_array
module procedure get_param_real
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  real :: new_name_value  ! The value that is set when the old name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value = value
      call read_param_real(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_real(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = (new_name_used .and. old_name_used .and. (value == new_name_value))
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_real(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_real(CS, modulename, varname, value, desc, units, &
                        default, debuggingParam)
  endif

  if (present(unscaled)) unscaled = value
  if (present(scale)) value = scale*value

end procedure get_param_real
module procedure get_param_real_array
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  real    :: new_name_value(size(value))  ! The values that are set when the standard name is used.
  integer :: m
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (present(defaults)) then
    if (present(default)) call MOM_error(FATAL, &
          "get_param_real_array: Only one of default and defaults can be specified at a time.")
    if (size(defaults) /= size(value)) call MOM_error(FATAL, &
          "get_param_real_array: The size of defaults and value are not the same.")
  endif

  if (do_read) then
    if (present(default)) value(:) = default
    if (present(defaults)) value(:) = defaults(:)

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value(:) = value(:)
      call read_param_real_array(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_real_array(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = .true.
        do m=1,size(value) ; if (value(m) /= new_name_value(m)) same_value = .false. ; enddo
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_real_array(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_real_array(CS, modulename, varname, value, desc, &
                              units, default, defaults, debuggingParam)
  endif

  if (present(unscaled)) unscaled(:) = value(:)
  if (present(scale)) value(:) = scale*value(:)

end procedure get_param_real_array
module procedure get_param_char
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  character(len=:), allocatable :: new_name_value  ! The value that is set when the standard name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value = value
      call read_param_char(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_char(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = (trim(value) == trim(new_name_value))
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_char(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_char(CS, modulename, varname, value, desc, units, &
                        default, layoutParam, debuggingParam)
  endif

end procedure get_param_char
module procedure get_param_char_array
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  integer :: i, m, len_tot, len_val
  character(len=:), allocatable :: cat_val
  character(len=:), allocatable :: new_name_value(:)  ! The value that is set when the standard name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value(:) = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value(:) = value(:)
      call read_param_char_array(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_char_array(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = .true.
        do m=1,size(value) ; if (trim(value(m)) /= trim(new_name_value(m))) same_value = .false. ; enddo
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_char_array(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    cat_val = trim(value(1)) ; len_tot = len_trim(value(1))
    do i=2,size(value)
      len_val = len_trim(value(i))
      if ((len_val > 0) .and. (len_tot + len_val + 2 < 240)) then
        cat_val = trim(cat_val)//ACHAR(34)// ", "//ACHAR(34)//trim(value(i))
        len_tot = len_tot + len_val
      endif
    enddo
    call log_param_char(CS, modulename, varname, cat_val, desc, &
                        units, default)
  endif

end procedure get_param_char_array
module procedure get_param_logical
  logical :: do_read, do_log
  logical :: new_name_used, old_name_used, same_value
  logical :: new_name_value  ! The value that is set when the standard name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value = value
      call read_param_logical(CS, old_name, value, set=old_name_used)
      if (old_name_used) then
        call read_param_logical(CS, varname, new_name_value, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = (value .eqv. new_name_value)
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_logical(CS, varname, value, fail_if_missing)
    endif
  endif

  if (do_log) then
    call log_param_logical(CS, modulename, varname, value, desc, &
                           units, default, layoutParam, debuggingParam)
  endif

end procedure get_param_logical
module procedure get_param_time
  logical :: do_read, do_log, log_date
  logical :: new_name_used, old_name_used, same_value
  type(time_type) :: new_name_value  ! The value that is set when the standard name is used.
  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log
  log_date = .false.

  if (do_read) then
    if (present(default)) value = default

    old_name_used = .false.
    if (present(old_name)) then
      new_name_value = value
      call read_param_time(CS, old_name, value, timeunit, date_format=log_date, set=old_name_used)
      if (old_name_used) then
        call read_param_time(CS, varname, new_name_value, timeunit, date_format=log_date, set=new_name_used)

        ! Issue appropriate warnings or error messages.
        same_value = (value == new_name_value)
        call archaic_param_name_message(varname, old_name, new_name_used, same_value)
      endif
    endif

    if (.not.old_name_used) then ! Old name is either not present or not set.
      call read_param_time(CS, varname, value, timeunit, fail_if_missing, date_format=log_date)
    endif
  endif

  if (do_log) then
    if (present(log_as_date)) log_date = log_as_date
    call log_param_time(CS, modulename, varname, value, desc, units, default, &
                        timeunit, layoutParam=layoutParam, &
                        debuggingParam=debuggingParam, log_date=log_date)
  endif

end procedure get_param_time
module procedure archaic_param_name_message
  if (new_name_used .and. same_value) then
    call MOM_error(WARNING, "The runtime parameter "//trim(varname)//&
                 " is also being set consistently via its older name of "//trim(old_name)//&
                 ".  Please migrate to only using "//trim(varname)//".")
  elseif (new_name_used .and. .not.same_value) then
    call MOM_error(FATAL, "The runtime parameter "//trim(varname)//&
                 " is also being set inconsistently via its older name of "//trim(old_name)//&
                 ".  Only use "//trim(varname)//".")
  else
    call MOM_error(WARNING, "The runtime parameter "//trim(varname)//&
                   " is being set via its soon to be obsolete name of "//trim(old_name)//&
                   ".  Please migrate to using "//trim(varname)//".")
  endif
end procedure archaic_param_name_message
module procedure clearParameterBlock
  type(parameter_block), pointer :: block => NULL()
  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = ''
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'clearParameterBlock: A clear was attempted before allocation.')
  endif
end procedure clearParameterBlock
module procedure openParameterBlock
  type(parameter_block), pointer :: block => NULL()
  logical :: do_log
  do_log = .true.
  if (present(do_not_log)) do_log = .not. do_not_log

  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = pushBlockLevel(block%name,blockName)
    if (do_log) then
      call doc_openBlock(CS%doc, block%name, desc)
      block%log_access = .true.
    else
      block%log_access = .false.
    endif
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'openParameterBlock: A push was attempted before allocation.')
  endif
end procedure openParameterBlock
module procedure closeParameterBlock
  type(parameter_block), pointer :: block => NULL()
  if (associated(CS%blockName)) then
    block => CS%blockName
    if (is_root_pe().and.len_trim(block%name)==0) call MOM_error(FATAL, &
      'closeParameterBlock: A pop was attempted on an empty stack. ("'//&
      trim(block%name)//'")')
    if (block%log_access) call doc_closeBlock(CS%doc, block%name)
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'closeParameterBlock: A pop was attempted before allocation.')
  endif
  block%name = popBlockLevel(block%name)
end procedure closeParameterBlock
module procedure pushBlockLevel
  if (len_trim(oldBlockName)>0) then
    pushBlockLevel=trim(oldBlockName)//'%'//trim(newBlockName)
  else
    pushBlockLevel=trim(newBlockName)
  endif
end procedure pushBlockLevel
module procedure popBlockLevel
  integer :: i
  i = index(trim(oldBlockName), '%', .true.)
  if (i>1) then
    popBlockLevel = trim(oldBlockName(1:i-1))
  elseif (i==0) then
    popBlockLevel = ''
  else ! i==1
    if (is_root_pe()) call MOM_error(FATAL, &
      'popBlockLevel: A pop was attempted leaving an empty block name.')
  endif
end procedure popBlockLevel
end submodule MOM_file_parser_s
