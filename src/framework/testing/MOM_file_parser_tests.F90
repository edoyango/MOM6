! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module MOM_file_parser_tests

use posix, only : chmod

use MOM_file_parser, only : param_file_type
use MOM_file_parser, only : open_param_file
use MOM_file_parser, only : close_param_file
use MOM_file_parser, only : read_param
use MOM_file_parser, only : log_param
use MOM_file_parser, only : get_param
use MOM_file_parser, only : log_version
use MOM_file_parser, only : clearParameterBlock
use MOM_file_parser, only : openParameterBlock
use MOM_file_parser, only : closeParameterBlock

use MOM_time_manager, only : time_type
use MOM_time_manager, only : set_date
use MOM_time_manager, only : set_ticks_per_second
use MOM_time_manager, only : set_calendar_type
use MOM_time_manager, only : NOLEAP, NO_CALENDAR

use MOM_error_handler, only : assert
use MOM_error_handler, only : MOM_error
use MOM_error_handler, only : FATAL

use MOM_unit_testing, only : TestSuite
use MOM_unit_testing, only : string
use MOM_unit_testing, only : create_test_file
use MOM_unit_testing, only : delete_test_file

implicit none ; private

public :: run_file_parser_tests

character(len=*), parameter :: param_filename = 'TEST_input'
character(len=*), parameter :: missing_param_filename = 'MISSING_input'
character(len=*), parameter :: netcdf_param_filename = 'TEST_input.nc'

character(len=*), parameter :: sample_param_name = 'SAMPLE_PARAMETER'
character(len=*), parameter :: missing_param_name = 'MISSING_PARAMETER'

character(len=*), parameter :: module_name = "SAMPLE_module"
character(len=*), parameter :: module_version = "SAMPLE_version"
character(len=*), parameter :: module_desc = "Description here"

character(len=9), parameter :: param_docfiles(4) = [ &
  "all      ", &
  "debugging", &
  "layout   ", &
  "short    " &
]


  interface
module subroutine test_open_param_file

end subroutine test_open_param_file
module subroutine test_close_param_file_quiet

end subroutine test_close_param_file_quiet
module subroutine test_open_param_file_component

end subroutine test_open_param_file_component
module subroutine cleanup_open_param_file_component

end subroutine cleanup_open_param_file_component
module subroutine test_open_param_file_docdir
  ! TODO: Make a new directory...?

end subroutine test_open_param_file_docdir
module subroutine test_open_param_file_empty_filename

end subroutine test_open_param_file_empty_filename
module subroutine test_open_param_file_long_name
  !> Store filename in a variable longer than FILENAME_LENGTH

end subroutine test_open_param_file_long_name
module subroutine test_missing_param_file

end subroutine test_missing_param_file
module subroutine test_open_param_file_ioerr
  ! NOTE: Induce an I/O error in open() by making the file unreadable

end subroutine test_open_param_file_ioerr
module subroutine cleanup_open_param_file_ioerr

end subroutine cleanup_open_param_file_ioerr
module subroutine test_open_param_file_netcdf

end subroutine test_open_param_file_netcdf
module subroutine cleanup_open_param_file_netcdf

end subroutine cleanup_open_param_file_netcdf
module subroutine test_open_param_file_checkable

end subroutine test_open_param_file_checkable
module subroutine test_reopen_param_file

end subroutine test_reopen_param_file
module subroutine test_open_param_file_no_doc

end subroutine test_open_param_file_no_doc
module subroutine test_read_param_int

end subroutine test_read_param_int
module subroutine test_read_param_int_missing

end subroutine test_read_param_int_missing
module subroutine test_read_param_int_undefined

end subroutine test_read_param_int_undefined
module subroutine test_read_param_int_type_err

end subroutine test_read_param_int_type_err
module subroutine test_read_param_int_array

end subroutine test_read_param_int_array
module subroutine test_read_param_int_array_missing

end subroutine test_read_param_int_array_missing
module subroutine test_read_param_int_array_undefined

end subroutine test_read_param_int_array_undefined
module subroutine test_read_param_int_array_type_err

end subroutine test_read_param_int_array_type_err
module subroutine test_read_param_real

end subroutine test_read_param_real
module subroutine test_read_param_real_missing

end subroutine test_read_param_real_missing
module subroutine test_read_param_real_undefined

end subroutine test_read_param_real_undefined
module subroutine test_read_param_real_type_err

end subroutine test_read_param_real_type_err
module subroutine test_read_param_real_array

end subroutine test_read_param_real_array
module subroutine test_read_param_real_array_missing

end subroutine test_read_param_real_array_missing
module subroutine test_read_param_real_array_undefined

end subroutine test_read_param_real_array_undefined
module subroutine test_read_param_real_array_type_err

end subroutine test_read_param_real_array_type_err
module subroutine test_read_param_logical

end subroutine test_read_param_logical
module subroutine test_read_param_logical_missing

end subroutine test_read_param_logical_missing
module subroutine test_read_param_char_no_delim

end subroutine test_read_param_char_no_delim
module subroutine test_read_param_char_quote_delim

end subroutine test_read_param_char_quote_delim
module subroutine test_read_param_char_apostrophe_delim

end subroutine test_read_param_char_apostrophe_delim
module subroutine test_read_param_char_missing

end subroutine test_read_param_char_missing
module subroutine test_read_param_char_array

end subroutine test_read_param_char_array
module subroutine test_read_param_char_array_missing

end subroutine test_read_param_char_array_missing
module subroutine test_read_param_time_date

end subroutine test_read_param_time_date
module subroutine test_read_param_time_date_bad_format

end subroutine test_read_param_time_date_bad_format
module subroutine test_read_param_time_tuple

end subroutine test_read_param_time_tuple
module subroutine test_read_param_time_bad_tuple

end subroutine test_read_param_time_bad_tuple
module subroutine test_read_param_time_bad_tuple_values

end subroutine test_read_param_time_bad_tuple_values
module subroutine test_read_param_time_unit

end subroutine test_read_param_time_unit
module subroutine test_read_param_time_missing

end subroutine test_read_param_time_missing
module subroutine test_read_param_time_undefined

end subroutine test_read_param_time_undefined
module subroutine test_read_param_time_type_err

end subroutine test_read_param_time_type_err
module subroutine test_read_param_unused_fatal

end subroutine test_read_param_unused_fatal
module subroutine test_read_param_replace_tabs

end subroutine test_read_param_replace_tabs
module subroutine test_read_param_pad_equals

end subroutine test_read_param_pad_equals
module subroutine test_read_param_multiline_param

end subroutine test_read_param_multiline_param
module subroutine test_read_param_multiline_param_unclosed

end subroutine test_read_param_multiline_param_unclosed
module subroutine test_read_param_multiline_comment


end subroutine test_read_param_multiline_comment
module subroutine test_read_param_multiline_comment_unclosed

end subroutine test_read_param_multiline_comment_unclosed
module subroutine test_read_param_misplaced_quote

end subroutine test_read_param_misplaced_quote
module subroutine test_read_param_define

end subroutine test_read_param_define
module subroutine test_read_param_define_as_flag

end subroutine test_read_param_define_as_flag
module subroutine test_read_param_override

end subroutine test_read_param_override
module subroutine test_read_param_override_misplaced

end subroutine test_read_param_override_misplaced
module subroutine test_read_param_override_twice

end subroutine test_read_param_override_twice
module subroutine test_read_param_override_repeat

end subroutine test_read_param_override_repeat
module subroutine test_read_param_override_warn_chain

end subroutine test_read_param_override_warn_chain
module subroutine test_read_param_assign_after_override

end subroutine test_read_param_assign_after_override
module subroutine test_read_param_override_no_def

end subroutine test_read_param_override_no_def
module subroutine test_read_param_assign_twice

end subroutine test_read_param_assign_twice
module subroutine test_read_param_assign_repeat

end subroutine test_read_param_assign_repeat
module subroutine test_read_param_null_stmt

end subroutine test_read_param_null_stmt
module subroutine test_read_param_assign_in_define

end subroutine test_read_param_assign_in_define
module subroutine test_read_param_block

end subroutine test_read_param_block
module subroutine test_read_param_block_stack

end subroutine test_read_param_block_stack
module subroutine test_read_param_block_inline_stack

end subroutine test_read_param_block_inline_stack
module subroutine test_read_param_block_empty_pop

end subroutine test_read_param_block_empty_pop
module subroutine test_read_param_block_close_unnamed

end subroutine test_read_param_block_close_unnamed
module subroutine test_read_param_block_close_unopened

end subroutine test_read_param_block_close_unopened
module subroutine test_read_param_block_unmatched

end subroutine test_read_param_block_unmatched
module subroutine test_open_unallocated_block

end subroutine test_open_unallocated_block
module subroutine test_close_unallocated_block

end subroutine test_close_unallocated_block
module subroutine test_clear_unallocated_block

end subroutine test_clear_unallocated_block
module subroutine test_read_param_block_outside_block

end subroutine test_read_param_block_outside_block
module subroutine test_log_version_cs

end subroutine test_log_version_cs
module subroutine test_log_version_plain
end subroutine test_log_version_plain
module subroutine test_log_param_int

end subroutine test_log_param_int
module subroutine test_log_param_int_array

end subroutine test_log_param_int_array
module subroutine test_log_param_real

end subroutine test_log_param_real
module subroutine test_log_param_real_array

end subroutine test_log_param_real_array
module subroutine test_log_param_time

end subroutine test_log_param_time
module subroutine test_log_param_time_as_date

end subroutine test_log_param_time_as_date
module subroutine test_log_param_time_as_date_default

end subroutine test_log_param_time_as_date_default
module subroutine test_log_param_time_as_date_tick

end subroutine test_log_param_time_as_date_tick
module subroutine test_log_param_time_with_unit

end subroutine test_log_param_time_with_unit
module subroutine test_log_param_time_with_timeunit

end subroutine test_log_param_time_with_timeunit
module subroutine test_get_param_int

end subroutine test_get_param_int
module subroutine test_get_param_int_no_read_no_log

end subroutine test_get_param_int_no_read_no_log
module subroutine test_get_param_int_array

end subroutine test_get_param_int_array
module subroutine test_get_param_int_array_no_read_no_log

end subroutine test_get_param_int_array_no_read_no_log
module subroutine test_get_param_real

end subroutine test_get_param_real
module subroutine test_get_param_real_no_read_no_log

end subroutine test_get_param_real_no_read_no_log
module subroutine test_get_param_real_array

end subroutine test_get_param_real_array
module subroutine test_get_param_real_array_no_read_no_log

end subroutine test_get_param_real_array_no_read_no_log
module subroutine test_get_param_char

end subroutine test_get_param_char
module subroutine test_get_param_char_no_read_no_log

end subroutine test_get_param_char_no_read_no_log
module subroutine test_get_param_char_array

end subroutine test_get_param_char_array
module subroutine test_get_param_logical

end subroutine test_get_param_logical
module subroutine test_get_param_logical_no_read_no_log

end subroutine test_get_param_logical_no_read_no_log
module subroutine test_get_param_logical_default

end subroutine test_get_param_logical_default
module subroutine test_get_param_time

end subroutine test_get_param_time
module subroutine test_get_param_time_no_read_no_log

end subroutine test_get_param_time_no_read_no_log
module subroutine cleanup_file_parser

end subroutine cleanup_file_parser
module subroutine run_file_parser_tests
  ! testing...

  ! Delete any pre-existing test parameter files
end subroutine run_file_parser_tests
  end interface

end module MOM_file_parser_tests
