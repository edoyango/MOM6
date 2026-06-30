submodule (MOM_document) MOM_document_s
  implicit none
contains
module procedure doc_param_none
  integer :: numspc
  character(len=mLen) :: mesg
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    numspc = max(1,doc%commentColumn-8-len_trim(varname))
    mesg = "#define "//trim(varname)//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc)
  endif
end procedure doc_param_none
module procedure doc_param_logical
  character(len=mLen) :: mesg
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    if (val) then
      mesg = define_string(doc, varname, STRING_TRUE, units)
    else
      mesg = undef_string(doc, varname, units)
    endif

    equalsDefault = .false.
    if (present(like_default)) equalsDefault = like_default
    if (present(default)) then
      if (val .eqv. default) equalsDefault = .true.
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end procedure doc_param_logical
module procedure doc_param_logical_array
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen) :: valstring
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    if (vals(1)) then ; valstring = STRING_TRUE ; else ; valstring = STRING_FALSE ; endif
    do i=2,min(size(vals),128)
      if (vals(i)) then
        valstring = trim(valstring)//", "//STRING_TRUE
      else
        valstring = trim(valstring)//", "//STRING_FALSE
      endif
    enddo

    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) .neqv. default) equalsDefault = .false. ; enddo
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif
    if (present(like_default)) then ; if (like_default) equalsDefault = .true. ; endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end procedure doc_param_logical_array
module procedure doc_param_int
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = int_string(val)
    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(like_default)) equalsDefault = like_default
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end procedure doc_param_int
module procedure doc_param_int_array
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen)  :: valstring
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = int_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(int_string(vals(i)))
    enddo

    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif
    if (present(defaults)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= defaults(i)) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//trim(int_array_string(defaults))
    endif
    if (present(like_default)) then ; if (like_default) equalsDefault = .true. ; endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif

end procedure doc_param_int_array
module procedure doc_param_real
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = real_string(val)
    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(like_default)) equalsDefault = like_default
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif
end procedure doc_param_real
module procedure doc_param_real_array
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen) :: valstring
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = trim(real_array_string(vals(:)))

    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif
    if (present(defaults)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= defaults(i)) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//trim(real_array_string(defaults))
    endif
    if (present(like_default)) then ; if (like_default) equalsDefault = .true. ; endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif

end procedure doc_param_real_array
module procedure doc_param_char
  character(len=mLen) :: mesg
  logical :: equalsDefault
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = define_string(doc, varname, '"'//trim(val)//'"', units)

    equalsDefault = .false.
    if (present(like_default)) equalsDefault = like_default
    if (present(default)) then
      if (trim(val) == trim(default)) equalsDefault = .true.
      mesg = trim(mesg)//' default = "'//trim(adjustl(default))//'"'
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif

end procedure doc_param_char
module procedure doc_openBlock
  character(len=mLen) :: mesg
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = trim(blockName)//'%'

    if (present(desc)) then
      call writeMessageAndDesc(doc, mesg, desc)
    else
      call writeMessageAndDesc(doc, mesg, '')
    endif
  endif
  doc%blockPrefix = trim(doc%blockPrefix)//trim(blockName)//'%'
end procedure doc_openBlock
module procedure doc_closeBlock
  character(len=mLen) :: mesg
  integer :: i
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = '%'//trim(blockName)

    call writeMessageAndDesc(doc, mesg, '')
  endif
  i = index(trim(doc%blockPrefix), trim(blockName)//'%', .true.)
  if (i>1) then
    doc%blockPrefix = trim(doc%blockPrefix(1:i-1))
  else
    doc%blockPrefix = ''
  endif
end procedure doc_closeBlock
module procedure doc_param_time
  character(len=mLen)              :: mesg          ! The output message
  character(len=doc%commentColumn) :: valstring     ! A string with the formatted value.
  logical                          :: equalsDefault ! True if val = default.
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = time_string(val)
    if (present(units)) then
      mesg = define_string(doc, varname, valstring, units)
    else
      mesg = define_string(doc, varname, valstring, "[days : seconds]")
    endif

    equalsDefault = .false.
    if (present(like_default)) equalsDefault = like_default
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//trim(time_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif

end procedure doc_param_time
module procedure writeMessageAndDesc
  character(len=mLen) :: mesg          ! A full line of a message including indents.
  character(len=mLen) :: mesg_text     ! A line of message text without preliminary indents.
  integer :: start_ind = 1             ! The starting index in the description for the next line.
  integer :: nl_ind, tab_ind, end_ind  ! The indices of new-lines, tabs, and the end of a line.
  integer :: len_text, len_tab, len_nl ! The lengths of the text string, tabs and new-lines.
  integer :: len_cor                   ! The permitted length corrected for tab sizes in a line.
  integer :: len_desc                  ! The non-whitespace length of the description.
  integer :: substr_start              ! The starting index of a substring to search for tabs.
  integer :: indnt, msg_pad            ! Space counts used to format a message.
  logical :: msg_done, reset_msg_pad   ! Logicals used to format messages.
  logical :: all, short, layout, debug ! Flags indicating which files to write into.
  layout = .false. ; if (present(layoutParam)) layout = layoutParam
  debug = .false. ; if (present(debuggingParam)) debug = debuggingParam
  all = doc%complete .and. (doc%unitAll > 0) .and. .not. (layout .or. debug)
  short = doc%minimal .and. (doc%unitShort > 0) .and. .not. (layout .or. debug)
  if (present(valueWasDefault)) short = short .and. (.not. valueWasDefault)

  if (all) write(doc%unitAll, '(a)') trim(vmesg)
  if (short) write(doc%unitShort, '(a)') trim(vmesg)
  if (layout) write(doc%unitLayout, '(a)') trim(vmesg)
  if (debug) write(doc%unitDebugging, '(a)') trim(vmesg)

  if (len_trim(desc) == 0) return

  len_tab = len_trim("_\t_") - 2
  len_nl = len_trim("_\n_") - 2

  indnt = doc%commentColumn ; if (present(indent)) indnt = indent
  len_text = doc%max_line_len - (indnt + 2)
  start_ind = 1 ; msg_pad = 0 ; msg_done = .false.
  do
    if (len_trim(desc(start_ind:)) < 1) exit

    len_cor = len_text - msg_pad

    substr_start = start_ind
    len_desc = len_trim(desc)
    do ! Adjust the available line length for anomalies in the size of tabs, counting \t as 2 spaces.
      if (substr_start >= start_ind+len_cor) exit
      tab_ind = index(desc(substr_start:min(len_desc,start_ind+len_cor)), "\t")
      if (tab_ind == 0) exit
      substr_start = substr_start + tab_ind
      len_cor = len_cor + (len_tab - 2)
    enddo

    nl_ind = index(desc(start_ind:), "\n")
    end_ind = 0
    if ((nl_ind > 0) .and. (len_trim(desc(start_ind:start_ind+nl_ind-2)) > len_cor)) then
      ! This line is too long despite the new-line character.  Look for an earlier space to break.
      end_ind = scan(desc(start_ind:start_ind+len_cor), " ", back=.true.) - 1
      if (end_ind > 0) nl_ind = 0
    elseif ((nl_ind == 0) .and. (len_trim(desc(start_ind:)) > len_cor)) then
      ! This line is too long and does not have a new-line character.  Look for a space to break.
      end_ind = scan(desc(start_ind:start_ind+len_cor), " ", back=.true.) - 1
    endif

    reset_msg_pad = .false.
    if (nl_ind > 0) then
      mesg_text = trim(desc(start_ind:start_ind+nl_ind-2))
      start_ind = start_ind + nl_ind + len_nl - 1
      reset_msg_pad = .true.
    elseif (end_ind > 0) then
      mesg_text = trim(desc(start_ind:start_ind+end_ind))
      start_ind = start_ind + end_ind + 1
      ! Adjust the starting point to move past leading spaces.
      start_ind = start_ind + (len_trim(desc(start_ind:)) - len_trim(adjustl(desc(start_ind:))))
    else
      mesg_text = trim(desc(start_ind:))
      msg_done = .true.
    endif

    do ; tab_ind = index(mesg_text, "\t") ! Replace \t with 2 spaces.
      if (tab_ind == 0) exit
      mesg_text(tab_ind:) = "  "//trim(mesg_text(tab_ind+len_tab:))
    enddo

    mesg = repeat(" ",indnt)//"! "//repeat(" ",msg_pad)//trim(mesg_text)

    if (reset_msg_pad) then
      msg_pad = 0
    elseif (msg_pad == 0) then ! Indent continuation lines.
      msg_pad = len_trim(mesg_text) - len_trim(adjustl(mesg_text))
      ! If already indented, indent an additional 2 spaces.
      if (msg_pad >= 2) msg_pad = msg_pad + 2
    endif

    if (all) write(doc%unitAll, '(a)') trim(mesg)
    if (short) write(doc%unitShort, '(a)') trim(mesg)
    if (layout) write(doc%unitLayout, '(a)') trim(mesg)
    if (debug) write(doc%unitDebugging, '(a)') trim(mesg)

    if (msg_done) exit
  enddo

end procedure writeMessageAndDesc
module procedure time_string
  integer :: secs, days, ticks, ticks_per_sec
  call get_time(Time, secs, days, ticks)

  time_string = trim(adjustl(int_string(days))) // ":" // trim(adjustl(int_string(secs)))
  if (ticks /= 0) then
    ticks_per_sec = get_ticks_per_second()
    time_string = trim(time_string) // ":" // &
                  trim(adjustl(int_string(ticks)))//"/"//trim(adjustl(int_string(ticks_per_sec)))
  endif

end procedure time_string
module procedure real_string
  integer :: len, ind
  if ((abs(val) < 1.0e4) .and. (abs(val) >= 1.0e-3)) then
    write(real_string, '(F30.11)') val
    if (.not.testFormattedFloatIsReal(real_string,val)) then
      write(real_string, '(F30.12)') val
      if (.not.testFormattedFloatIsReal(real_string,val)) then
        write(real_string, '(F30.13)') val
        if (.not.testFormattedFloatIsReal(real_string,val)) then
          write(real_string, '(F30.14)') val
          if (.not.testFormattedFloatIsReal(real_string,val)) then
            write(real_string, '(F30.15)') val
            if (.not.testFormattedFloatIsReal(real_string,val)) then
              write(real_string, '(F30.16)') val
            endif
          endif
        endif
      endif
    endif
    do
      len = len_trim(real_string)
      if ((len<2) .or. (real_string(len-1:len) == ".0") .or. &
          (real_string(len:len) /= "0")) exit
      real_string(len:len) = " "
    enddo
  elseif (val == 0.) then
    real_string = "0.0"
  else
    if ((abs(val) < 1.0e-99) .or. (abs(val) >= 1.0e100)) then
      write(real_string(1:32), '(ES24.14E4)') val
      if (scan(real_string, "eE") == 0) then  ! Fix a bug with a missing E in PGI formatting
        ind = scan(real_string, "-+", back=.true.)
        if (ind > index(real_string, ".") ) &  ! Avoid changing a leading sign.
          real_string = real_string(1:ind-1)//"E"//real_string(ind:)
      endif
      if (.not.testFormattedFloatIsReal(real_string, val)) then
        write(real_string(1:32), '(ES25.15E4)') val
        if (scan(real_string, "eE") == 0) then  ! Fix a bug with a missing E in PGI formatting
          ind = scan(real_string, "-+", back=.true.)
          if (ind > index(real_string, ".") ) &  ! Avoid changing a leading sign.
            real_string = real_string(1:ind-1)//"E"//real_string(ind:)
        endif
      endif
      ! Remove a leading 0 from the exponent, if it is there.
      ind = max(index(real_string, "E+0"), index(real_string, "E-0"))
      if (ind > 0) real_string = real_string(1:ind+1)//real_string(ind+3:)
    else
      write(real_string(1:32), '(ES23.14)') val
      if (.not.testFormattedFloatIsReal(real_string, val)) &
        write(real_string(1:32), '(ES23.15)') val
    endif
    do  ! Remove extra trailing 0s before the exponent.
      ind = index(real_string, "0E")
      if (ind == 0) exit
      if (real_string(ind-1:ind-1) == ".") exit ! Leave at least one digit after the decimal point.
      real_string = real_string(1:ind-1)//real_string(ind+1:)
    enddo
  endif
  real_string = adjustl(real_string)
end procedure real_string
module procedure real_array_string
  integer :: j, n, ns
  logical :: doWrite
  character(len=10) :: separator
  n = 1 ; doWrite = .true. ; real_array_string = ''
  if (present(sep)) then
    separator = sep ; ns = len(sep)
  else
    separator = ', ' ; ns = 2
  endif
  do j=1,size(vals)
    doWrite = .true.
    if (j < size(vals)) then
      if (vals(j) == vals(j+1)) then
        n = n+1
        doWrite = .false.
      endif
    endif
    if (doWrite) then
      if (len(real_array_string) > 0) then ! Write separator if a number has already been written
        real_array_string = real_array_string // separator(1:ns)
      endif
      if (n>1) then
        real_array_string = real_array_string // trim(int_string(n)) // "*" // trim(real_string(vals(j)))
      else
        real_array_string = real_array_string // trim(real_string(vals(j)))
      endif
      n=1
    endif
  enddo
end procedure real_array_string
module procedure int_array_string
  integer :: j, m, n, ns
  logical :: doWrite
  character(len=10) :: separator
  n = 1 ; doWrite = .true. ; int_array_string = ''
  if (present(sep)) then
    separator = sep ; ns = len(sep)
  else
    separator = ', ' ; ns = 2
  endif
  do j=1,size(vals)
    doWrite = .true.
    if (j < size(vals)) then
      if (vals(j) == vals(j+1)) then
        n = n+1
        doWrite = .false.
      endif
    endif
    if (doWrite) then
      if (len(int_array_string) > 0) then ! Write separator if a number has already been written
        int_array_string = int_array_string // separator(1:ns)
      endif
      if (n>1) then
        if (size(vals) > 6) then  ! The n*val syntax is convenient in long lists of integers.
          int_array_string = int_array_string // trim(int_string(n)) // "*" // trim(int_string(vals(j)))
        else  ! For short lists of integers, do not use the n*val syntax as it is less convenient.
          do m=1,n-1
            int_array_string = int_array_string // trim(int_string(vals(j))) // separator(1:ns)
          enddo
          int_array_string = int_array_string // trim(int_string(vals(j)))
        endif
      else
        int_array_string = int_array_string // trim(int_string(vals(j)))
      endif
      n=1
    endif
  enddo
end procedure int_array_string
module procedure testFormattedFloatIsReal
  real :: scannedVal
  read(str(1:),*) scannedVal
  if (scannedVal == val) then
    testFormattedFloatIsReal=.true.
  else
    testFormattedFloatIsReal=.false.
  endif
end procedure testFormattedFloatIsReal
module procedure int_string
  write(int_string, '(i24)') val
  int_string = adjustl(int_string)
end procedure int_string
module procedure logical_string
  write(logical_string, '(l24)') val
  logical_string = adjustl(logical_string)
end procedure logical_string
module procedure define_string
  integer :: numSpaces
  define_string = repeat(" ",mLen) ! Blank everything for safety
  if (doc%defineSyntax) then
    define_string = "#define "//trim(varName)//" "//valString
  else
    define_string = trim(varName)//" = "//valString
  endif
  numSpaces = max(1, doc%commentColumn - len_trim(define_string) )
  define_string = trim(define_string)//repeat(" ",numSpaces)//"!"
  if (len_trim(units) > 0) define_string = trim(define_string)//"   ["//trim(units)//"]"
end procedure define_string
module procedure undef_string
  integer :: numSpaces
  undef_string = repeat(" ",240) ! Blank everything for safety
  undef_string = "#undef "//trim(varName)
  if (doc%defineSyntax) then
    undef_string = "#undef "//trim(varName)
  else
    undef_string = trim(varName)//" = "//STRING_FALSE
  endif
  numSpaces = max(1, doc%commentColumn - len_trim(undef_string) )
  undef_string = trim(undef_string)//repeat(" ",numSpaces)//"!"
  if (len_trim(units) > 0) undef_string = trim(undef_string)//"   ["//trim(units)//"]"
end procedure undef_string
module procedure doc_module
  character(len=mLen) :: mesg
  logical :: repeat_doc
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    ! Add a blank line for delineation
    call writeMessageAndDesc(doc, '', '', valueWasDefault=all_default, &
                             layoutParam=layoutMod, debuggingParam=debuggingMod)
    mesg = "! === module "//trim(modname)//" ==="
    call writeMessageAndDesc(doc, mesg, desc, valueWasDefault=all_default, indent=0, &
                             layoutParam=layoutMod, debuggingParam=debuggingMod)
    if (present(log_to_all)) then ; if (log_to_all) then
      ! Log the module version again if the previous call was intercepted for use to document
      ! a layout or debugging module.
      repeat_doc = .false.
      if (present(layoutMod)) then ; if (layoutMod) repeat_doc = .true. ; endif
      if (present(debuggingMod)) then ; if (debuggingMod) repeat_doc = .true. ; endif
      if (repeat_doc) then
        call writeMessageAndDesc(doc, '', '', valueWasDefault=all_default)
        call writeMessageAndDesc(doc, mesg, desc, valueWasDefault=all_default, indent=0)
      endif
    endif ; endif
  endif
end procedure doc_module
module procedure doc_subroutine
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end procedure doc_subroutine
module procedure doc_function
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end procedure doc_function
module procedure doc_init
  if (.not. associated(doc)) then
    allocate(doc)
  endif

  doc%docFileBase = docFileBase
  if (present(minimal)) doc%minimal = minimal
  if (present(complete)) doc%complete = complete
  if (present(layout)) doc%layout = layout
  if (present(debugging)) doc%debugging = debugging

end procedure doc_init
module procedure open_doc_file
  logical :: opened, new_file
  integer :: ios
  character(len=240) :: fileName
  if (.not. (is_root_pe() .and. associated(doc))) return

  if ((len_trim(doc%docFileBase) > 0) .and. doc%complete .and. (doc%unitAll<0)) then
    new_file = .true. ; if (doc%unitAll /= -1) new_file = .false.
    doc%unitAll = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.all'
    if (new_file) then
      open(doc%unitAll, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitAll, '(a)') &
       '! This file was written by the model and records all non-layout '//&
       'or debugging parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitAll, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitAll, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%minimal .and. (doc%unitShort<0)) then
    new_file = .true. ; if (doc%unitShort /= -1) new_file = .false.
    doc%unitShort = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.short'
    if (new_file) then
      open(doc%unitShort, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitShort, '(a)') &
       '! This file was written by the model and records the non-default parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitShort, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitShort, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%layout .and. (doc%unitLayout<0)) then
    new_file = .true. ; if (doc%unitLayout /= -1) new_file = .false.
    doc%unitLayout = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.layout'
    if (new_file) then
      open(doc%unitLayout, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitLayout, '(a)') &
       '! This file was written by the model and records the layout parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitLayout, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitLayout, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%debugging .and. (doc%unitDebugging<0)) then
    new_file = .true. ; if (doc%unitDebugging /= -1) new_file = .false.
    doc%unitDebugging = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.debugging'
    if (new_file) then
      open(doc%unitDebugging, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitDebugging, '(a)') &
       '! This file was written by the model and records the debugging parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitDebugging, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitDebugging, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

end procedure open_doc_file
module procedure find_unused_unit_number
  logical :: opened
  do find_unused_unit_number=512,42,-1
    inquire( find_unused_unit_number, opened=opened)
    if (.not.opened) exit
  enddo
  if (opened) call MOM_error(FATAL, &
    "doc_init failed to find an unused unit number.")
end procedure find_unused_unit_number
module procedure doc_end
  type(link_msg), pointer :: this => NULL(), next => NULL()
  if (.not.associated(doc)) return

  if (doc%unitAll > 0) then
    close(doc%unitAll)
    doc%unitAll = -2
  endif

  if (doc%unitShort > 0) then
    close(doc%unitShort)
    doc%unitShort = -2
  endif

  if (doc%unitLayout > 0) then
    close(doc%unitLayout)
    doc%unitLayout = -2
  endif

  if (doc%unitDebugging > 0) then
    close(doc%unitDebugging)
    doc%unitDebugging = -2
  endif

  doc%filesAreOpen = .false.

  this => doc%chain_msg
  do while( associated(this) )
    next => this%next
    deallocate(this)
    this => next
  enddo
end procedure doc_end
module procedure mesgHasBeenDocumented
  type(link_msg), pointer :: newLink => NULL(), this => NULL(), last => NULL()
  mesgHasBeenDocumented = .false.

!!if (mesg(1:1) == '!') return ! Ignore commented parameters

  ! Search through list for this parameter
  last => NULL()
  this => doc%chain_msg
  do while( associated(this) )
    if (trim(doc%blockPrefix)//trim(varName) == trim(this%name)) then
      mesgHasBeenDocumented = .true.
      if (trim(mesg) == trim(this%msg)) return
      ! If we fail the above test then cause an error
      if (mesg(1:1) == '!') return ! Do not cause error for commented parameters
      call MOM_error(WARNING, "Previous msg:"//trim(this%msg))
      call MOM_error(WARNING, "New message :"//trim(mesg))
      call MOM_error(WARNING, "Encountered inconsistent documentation line for parameter "&
                     //trim(varName)//"!")
    endif
    last => this
    this => this%next
  enddo

  ! Allocate a new link
  allocate(newLink)
  newLink%name = trim(doc%blockPrefix)//trim(varName)
  newLink%msg = trim(mesg)
  newLink%next => NULL()
  if (.not. associated(doc%chain_msg)) then
    doc%chain_msg => newLink
  else
    if (.not. associated(last)) call MOM_error(FATAL, &
         "Unassociated LINK in mesgHasBeenDocumented: "//trim(mesg))
    last%next => newLink
  endif
end procedure mesgHasBeenDocumented
end submodule MOM_document_s
