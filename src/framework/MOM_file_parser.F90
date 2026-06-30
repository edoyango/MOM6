! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The MOM6 facility to parse input files for runtime parameters
module MOM_file_parser

use MOM_coms, only : root_PE, broadcast
use MOM_coms, only : any_across_PEs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, assert
use MOM_error_handler, only : is_root_pe, stdlog, stdout
use MOM_time_manager, only : get_time, time_type, get_ticks_per_second
use MOM_time_manager, only : set_date, get_date, real_to_time, operator(-), operator(==), set_time
use MOM_document, only : doc_param, doc_module, doc_init, doc_end, doc_type
use MOM_document, only : doc_openBlock, doc_closeBlock
use MOM_string_functions, only : left_int, left_ints, slasher
use MOM_string_functions, only : left_real, left_reals

implicit none ; private

! These are hard-coded limits that are used in the following code.  They should be set
! generously enough not to impose any significant limitations.
integer, parameter, public :: MAX_PARAM_FILES = 5 !< Maximum number of parameter files.
integer, parameter :: INPUT_STR_LENGTH = 1024 !< Maximum line length in parameter file.  Lines that
                                              !! are combined by ending in '\' or '&' can exceed
                                              !! this limit after merging.
integer, parameter :: FILENAME_LENGTH = 200   !< Maximum number of characters in file names.


!>@{ Default values for parameters
logical, parameter :: report_unused_default = .true.
logical, parameter :: unused_params_fatal_default = .false.
logical, parameter :: log_to_stdout_default = .false.
logical, parameter :: complete_doc_default = .true.
logical, parameter :: minimal_doc_default = .true.
!>@}


!> A simple type to allow lines in an array to be allocated with variable sizes.
type, private :: file_line_type ; private
  character(len=:), allocatable :: line !< An allocatable line with content
end type file_line_type

!> The valid lines extracted from an input parameter file without comments
type, private :: file_data_type ; private
  integer :: num_lines = 0 !< The number of lines in this type
  type(file_line_type), allocatable, dimension(:) :: fln !< Lines with the input content.
  logical,                  pointer, dimension(:) :: line_used => NULL() !< If true, the line has been read
end type file_data_type

!> A link in the list of variables that have already had override warnings issued
type, private :: link_parameter ; private
  type(link_parameter), pointer :: next => NULL() !< Facilitates linked list
  character(len=80) :: name                       !< Parameter name
  logical :: hasIssuedOverrideWarning = .false.   !< Has a default value
end type link_parameter

!> Specify the active parameter block
type, private :: parameter_block ; private
  character(len=240) :: name = ''   !< The active parameter block name
  logical :: log_access = .true.
    !< Log the entry and exit of the block (but not its contents)
end type parameter_block

!> A structure that can be parsed to read and document run-time parameters.
type, public :: param_file_type ; private
  integer  :: nfiles = 0            !< The number of open files.
  integer  :: iounit(MAX_PARAM_FILES) !< The unit numbers of open files.
  character(len=FILENAME_LENGTH)  :: filename(MAX_PARAM_FILES) !< The names of the open files.
  logical  :: NetCDF_file(MAX_PARAM_FILES) !< If true, the input file is in NetCDF.
                                    ! This is not yet implemented.
  type(file_data_type) :: param_data(MAX_PARAM_FILES) !< Structures that contain
                                    !! the valid data lines from the parameter
                                    !! files, enabling all subsequent reads of
                                    !! parameter data to occur internally.
  logical  :: report_unused = report_unused_default !< If true, report any
                                    !! parameter lines that are not used in the run.
  logical  :: unused_params_fatal = unused_params_fatal_default  !< If true, kill
                                    !! the run if there are any unused parameters.
  logical  :: log_to_stdout = log_to_stdout_default !< If true, all log
                                    !! messages are also sent to stdout.
  logical  :: log_open = .false.    !< True if the log file has been opened.
  integer  :: max_line_len = 4      !< The maximum number of characters in the lines
                                    !! in any of the files in this param_file_type after
                                    !! any continued lines have been combined.
  integer  :: stdout                !< The unit number from stdout().
  integer  :: stdlog                !< The unit number from stdlog().
  character(len=240) :: doc_file    !< A file where all run-time parameters, their
                                    !! settings and defaults are documented.
  logical  :: complete_doc = complete_doc_default !< If true, document all
                                    !! run-time parameters.
  logical  :: minimal_doc = minimal_doc_default !< If true, document only those
                                    !! run-time parameters that differ from defaults.
  type(doc_type), pointer :: doc => NULL() !< A structure that contains information
                                    !! related to parameter documentation.
  type(link_parameter), pointer :: chain => NULL() !< Facilitates linked list
  type(parameter_block), pointer :: blockName => NULL() !< Name of active parameter block
end type param_file_type

public read_param, open_param_file, close_param_file, log_param, log_version
public doc_param, get_param
public clearParameterBlock, openParameterBlock, closeParameterBlock

!> An overloaded interface to read various types of parameters
interface read_param
  module procedure read_param_int, read_param_real, read_param_logical, &
                   read_param_char, read_param_char_array, read_param_time, &
                   read_param_int_array, read_param_real_array
end interface
!> An overloaded interface to log the values of various types of parameters
interface log_param
  module procedure log_param_int, log_param_real, log_param_logical, &
                   log_param_char, log_param_time, &
                   log_param_int_array, log_param_real_array
end interface
!> An overloaded interface to read and log the values of various types of parameters
interface get_param
  module procedure get_param_int, get_param_real, get_param_logical, &
                   get_param_char, get_param_char_array, get_param_time, &
                   get_param_int_array, get_param_real_array
end interface

!> An overloaded interface to log version information about modules
interface log_version
  module procedure log_version_cs, log_version_plain
end interface


  interface
module subroutine open_param_file(filename, CS, checkable, component, doc_file_dir, ensemble_num)
  character(len=*),           intent(in) :: filename !< An input file name, optionally with the full path
  type(param_file_type),   intent(inout) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  logical,          optional, intent(in) :: checkable   !< If this is false, it disables checks of this
                                         !! file for unused parameters.  The default is True.
  character(len=*), optional, intent(in) :: component   !< If present, this component name is used
                                         !! to generate parameter documentation file names; the default is"MOM"
  character(len=*), optional, intent(in) :: doc_file_dir !< An optional directory in which to write out
                                         !! the documentation files.  The default is effectively './'.
  integer, optional, intent(in)          :: ensemble_num !< ensemble number to be appended to _doc filenames (optional)

  ! Local variables

end subroutine open_param_file
module subroutine close_param_file(CS, quiet_close, component)
  type(param_file_type),   intent(inout) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  logical,          optional, intent(in) :: quiet_close !< if present and true, do not do any
                                         !! logging with this call.
  character(len=*), optional, intent(in) :: component   !< If present, this component name is used
                                         !! to generate parameter documentation file names
  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine close_param_file
module subroutine populate_param_data(iounit, filename, param_data)
  integer,                 intent(in) :: iounit !< The IO unit number that is open for filename
  character(len=*),        intent(in) :: filename !< An input file name, optionally with the full path
  type(file_data_type), intent(inout) :: param_data !< A list of the input lines that set parameters
                                                !! after comments have been stripped out.

  ! Local variables

  ! Find the number of keyword lines in a parameter file
end subroutine populate_param_data
module function openMultiLineComment(string)
  character(len=*), intent(in) :: string  !< The input string to process
  logical                      :: openMultiLineComment

  ! Local variables

end function openMultiLineComment
module function closeMultiLineComment(string)
  character(len=*), intent(in) :: string  !< The input string to process
  logical                      :: closeMultiLineComment
! True if a */ appears on this line
end function closeMultiLineComment
module function lastNonCommentIndex(string)
  character(len=*), intent(in) :: string  !< The input string to process
  integer                      :: lastNonCommentIndex

  ! Local variables

  ! This subroutine is the only place where a comment needs to be defined
end function lastNonCommentIndex
module function lastNonCommentNonBlank(string)
  character(len=*), intent(in) :: string  !< The input string to process
  integer                      :: lastNonCommentNonBlank

end function lastNonCommentNonBlank
module function replaceTabs(string)
  character(len=*), intent(in) :: string  !< The input string to process
  character(len=len(string))   :: replaceTabs


end function replaceTabs
module function removeComments(string)
  character(len=*), intent(in) :: string  !< The input string to process
  character(len=len(string))   :: removeComments


end function removeComments
module function simplifyWhiteSpace(string)
  character(len=*), intent(in) :: string !< A string to modify to simplify white space
  character(len=len(string)+16)   :: simplifyWhiteSpace

  ! Local variables

end function simplifyWhiteSpace
module subroutine read_param_int(CS, varname, value, fail_if_missing, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  integer,             intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.
  ! Local variables

end subroutine read_param_int
module subroutine read_param_int_array(CS, varname, value, fail_if_missing, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  integer, dimension(:),  intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.
  ! Local variables

end subroutine read_param_int_array
module subroutine read_param_real(CS, varname, value, fail_if_missing, scale, set)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),      intent(in) :: varname !< The case-sensitive name of the parameter to read
  real,               intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,     optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  real,        optional, intent(in) :: scale   !< A scaling factor that the parameter is multiplied
                                         !! by before it is returned.
  logical,     optional, intent(out) :: set    !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.

  ! Local variables

end subroutine read_param_real
module subroutine read_param_real_array(CS, varname, value, fail_if_missing, scale, set)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),      intent(in) :: varname !< The case-sensitive name of the parameter to read
  real, dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,     optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  real,        optional, intent(in) :: scale   !< A scaling factor that the parameter is multiplied
                                         !! by before it is returned.
  logical,     optional, intent(out) :: set    !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.

  ! Local variables

end subroutine read_param_real_array
module subroutine read_param_char(CS, varname, value, fail_if_missing, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  character(len=*),    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.
  ! Local variables

end subroutine read_param_char
module subroutine read_param_char_array(CS, varname, value, fail_if_missing, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  character(len=*), dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.

  ! Local variables

end subroutine read_param_char_array
module subroutine read_param_logical(CS, varname, value, fail_if_missing, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  logical,             intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.

  ! Local variables

end subroutine read_param_logical
module subroutine read_param_time(CS, varname, value, timeunit, fail_if_missing, date_format, set)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  type(time_type),     intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  real,         optional, intent(in) :: timeunit !< The number of seconds in a time unit for real-number input.
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: date_format !< If present, this indicates whether this
                                         !! parameter was read in a date format, so that it can
                                         !! later be logged in the same format.
  logical,     optional, intent(out) :: set     !< If present, this indicates whether this parameter
                                         !! has been found and successfully set in the input files.

  ! Local variables

end subroutine read_param_time
module function strip_quotes(val_str)
  character(len=*), intent(in) :: val_str !< The character string to work on
  character(len=len(val_str)) :: strip_quotes
  ! Local variables
end function strip_quotes
module function max_input_line_length(CS, pf_num) result(max_len)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                                !! it is also a structure to parse for run-time parameters
  integer,      optional, intent(in) :: pf_num  !< If present, only work on a single file in the
                                                !! param_file_type, or return 0 if this exceeds the
                                                !! number of files in the param_file_type.
  integer :: max_len !< The maximum number of characters in any input lines after they
                     !! have been combined by any line continuation.

  ! Local variables

end function max_input_line_length
module subroutine get_variable_line(CS, varname, found, defined, value_string, paramIsLogical)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                                !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  logical,               intent(out) :: found   !< If true, this parameter has been found in CS
  logical,               intent(out) :: defined !< If true, this parameter is set (or true) in the CS
  character(len=*),      intent(out) :: value_string(:) !< A string that encodes the new value
  logical, optional,      intent(in) :: paramIsLogical  !< If true, this is a logical parameter
                                                !! that can be simply defined without parsing a value_string.

  ! Local variables
end subroutine get_variable_line
module subroutine flag_line_as_read(line_used, count)
  logical, dimension(:), pointer    :: line_used !< A structure indicating which lines have been read
  integer,               intent(in) :: count !< The parameter on this line number has been read
end subroutine flag_line_as_read
module function overrideWarningHasBeenIssued(chain, varName)
  type(link_parameter), pointer    :: chain   !< The linked list of variables that have already had
                                              !! override warnings issued
  character(len=*),     intent(in) :: varName !< The name of the variable being queried for warnings
  logical                          :: overrideWarningHasBeenIssued
  ! Local variables

end function overrideWarningHasBeenIssued
module subroutine log_version_cs(CS, modulename, version, desc, log_to_all, all_default, layout, debugging)
  type(param_file_type),      intent(in) :: CS         !< File parser type
  character(len=*),           intent(in) :: modulename !< Name of calling module
  character(len=*),           intent(in) :: version    !< Version string of module
  character(len=*), optional, intent(in) :: desc       !< Module description
  logical,          optional, intent(in) :: log_to_all !< If present and true, log this parameter to the
                                                       !! ..._doc.all files, even if this module also has layout
                                                       !! or debugging parameters.
  logical,          optional, intent(in) :: all_default !< If true, all parameters take their default values.
  logical,          optional, intent(in) :: layout     !< If present and true, this module has layout parameters.
  logical,          optional, intent(in) :: debugging  !< If present and true, this module has debugging parameters.
  ! Local variables

end subroutine log_version_cs
module subroutine log_version_plain(modulename, version)
  character(len=*),           intent(in) :: modulename !< Name of calling module
  character(len=*),           intent(in) :: version    !< Version string of module
  ! Local variables

end subroutine log_version_plain
module subroutine log_param_int(CS, modulename, varname, value, desc, units, &
                         default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the module using this parameter
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  integer,                    intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  integer,          optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.


end subroutine log_param_int
module subroutine log_param_int_array(CS, modulename, varname, value, desc, &
                               units, default, defaults, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the module using this parameter
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  integer, dimension(:),      intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  integer,          optional, intent(in) :: default !< The uniform default value of this parameter
  integer,          optional, intent(in) :: defaults(:) !< The element-wise default values of this parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.


end subroutine log_param_int_array
module subroutine log_param_real(CS, modulename, varname, value, desc, units, &
                          default, debuggingParam, like_default, unscale)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  real,                       intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*),           intent(in) :: units   !< The units of this parameter
  real,             optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.
  real,             optional, intent(in) :: unscale   !< A reciprocal scaling factor that the parameter is
                                         !! multiplied by before it is logged


end subroutine log_param_real
module subroutine log_param_real_array(CS, modulename, varname, value, desc, &
                                units, default, defaults, debuggingParam, like_default, unscale)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  real, dimension(:),         intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                             !! present, this parameter is not written to a doc file
  character(len=*),           intent(in) :: units   !< The units of this parameter
  real,             optional, intent(in) :: default !< A uniform default value of the parameter
  real,             optional, intent(in) :: defaults(:) !< The element-wise defaults of the parameter
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.
  real,             optional, intent(in) :: unscale   !< A reciprocal scaling factor that the parameter is
                                         !! multiplied by before it is logged


end subroutine log_param_real_array
module subroutine log_param_logical(CS, modulename, varname, value, desc, &
                             units, default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  logical,                    intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  logical,          optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.


end subroutine log_param_logical
module subroutine log_param_char(CS, modulename, varname, value, desc, units, &
                          default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  character(len=*),           intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  character(len=*), optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.


end subroutine log_param_char
module subroutine log_param_time(CS, modulename, varname, value, desc, units, &
                          default, timeunit, layoutParam, debuggingParam, log_date, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  type(time_type),            intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  type(time_type),  optional, intent(in) :: default !< The default value of the parameter
  real,             optional, intent(in) :: timeunit !< The number of seconds in a time unit for
                                         !! real-number output.
  logical,          optional, intent(in) :: log_date   !< If true, log the time_type in date format.
                                         !! If missing the default is false.
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  ! Local variables

end subroutine log_param_time
module function convert_date_to_string(date) result(date_string)
  type(time_type), intent(in) :: date !< The date to be translated into a string.
  character(len=40) :: date_string    !< A date string in a format like YYYY-MM-DD HH:MM:SS.sss

  ! Local variables

end function convert_date_to_string
module subroutine get_param_int(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               layoutParam, debuggingParam, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  integer,                    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  integer,          optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_int
module subroutine get_param_int_array(CS, modulename, varname, value, desc, units, &
               default, defaults, fail_if_missing, do_not_read, do_not_log, &
               layoutParam, debuggingParam, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  integer, dimension(:),      intent(inout) :: value   !< The value of the parameter that may be reset
                                         !! from the parameter file
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  integer,          optional, intent(in)    :: default !< The uniform default value of this parameter
  integer,          optional, intent(in)    :: defaults(:) !< The element-wise default values of this parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_int_array
module subroutine get_param_real(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               debuggingParam, scale, unscaled, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  real,                       intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*),           intent(in)    :: units   !< The units of this parameter
  real,             optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  real,             optional, intent(in)    :: scale   !< A scaling factor that the parameter is
                                         !! multiplied by before it is returned.
  real,             optional, intent(out)   :: unscaled !< The value of the parameter that would be
                                         !! returned without any multiplication by a scaling factor.
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_real
module subroutine get_param_real_array(CS, modulename, varname, value, desc, units, &
               default, defaults, fail_if_missing, do_not_read, do_not_log, debuggingParam, &
               scale, unscaled, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  real, dimension(:),         intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*),           intent(in)    :: units   !< The units of this parameter
  real,             optional, intent(in)    :: default !< A uniform default value of the parameter
  real,             optional, intent(in)    :: defaults(:) !< The element-wise defaults of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  real,             optional, intent(in)    :: scale   !< A scaling factor that the parameter is
                                         !! multiplied by before it is returned.
  real, dimension(:), optional, intent(out) :: unscaled !< The value of the parameter that would be
                                         !! returned without any multiplication by a scaling factor.
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_real_array
module subroutine get_param_char(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               layoutParam, debuggingParam, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  character(len=*),           intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  character(len=*), optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_char
module subroutine get_param_char_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  character(len=*), dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  character(len=*), optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_char_array
module subroutine get_param_logical(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               layoutParam, debuggingParam, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  logical,                    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  logical,          optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_logical
module subroutine get_param_time(CS, modulename, varname, value, desc, units, &
                          default, fail_if_missing, do_not_read, do_not_log, &
                          timeunit, layoutParam, debuggingParam, &
                          log_as_date, old_name)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  type(time_type),            intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  type(time_type),  optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  real,             optional, intent(in)    :: timeunit !< The number of seconds in a time unit for
                                         !! real-number input to be translated to a time.
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in)    :: log_as_date  !< If true, log the time_type in date
                                         !! format. The default is false.
  character(len=*), optional, intent(in)    :: old_name !< A case-sensitive archaic name of the parameter
                                         !! to read.  Errors or warnings are issued if the old name
                                         !! is being used.

  ! Local variables

end subroutine get_param_time
module subroutine archaic_param_name_message(varname, old_name, new_name_used, same_value)
  character(len=*), intent(in) :: varname  !< The case-sensitive name of the parameter to read
  character(len=*), intent(in) :: old_name !< The case-sensitive archaic name of the parameter
  logical,          intent(in) :: new_name_used  !< True if varname is used in the parameter file.
  logical,          intent(in) :: same_value !< True if varname and old_name give the same values.

end subroutine archaic_param_name_message
module subroutine clearParameterBlock(CS)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters

end subroutine clearParameterBlock
module subroutine openParameterBlock(CS, blockName, desc, do_not_log)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: blockName !< The name of a parameter block being added
  character(len=*), optional, intent(in) :: desc    !< A description of the parameter block being added
  logical, optional, intent(in) :: do_not_log
    !< Log block entry if true.  This only prevents logging of entry to the block, and not the contents.


end subroutine openParameterBlock
module subroutine closeParameterBlock(CS)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters


end subroutine closeParameterBlock
module function pushBlockLevel(oldblockName,newBlockName)
  character(len=*),        intent(in) :: oldBlockName  !< A sequence of hierarchical parameter block names
  character(len=*),        intent(in) :: newBlockName  !< A new block name to add to the end of the sequence
  character(len=len(oldBlockName)+40) :: pushBlockLevel

end function pushBlockLevel
module function popBlockLevel(oldblockName)
  character(len=*),        intent(in) :: oldBlockName !< A sequence of hierarchical parameter block names
  character(len=len(oldBlockName)+40) :: popBlockLevel

end function popBlockLevel
  end interface

end module MOM_file_parser
