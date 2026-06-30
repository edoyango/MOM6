! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The subroutines here provide hooks for document generation functions at
!! various levels of granularity.
module MOM_document

use MOM_time_manager,  only : time_type, operator(==), get_time, get_ticks_per_second
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe

implicit none ; private

public doc_param, doc_subroutine, doc_function, doc_module, doc_init, doc_end
public doc_openBlock, doc_closeBlock

!> Document parameter values
interface doc_param
  module procedure doc_param_none, &
                   doc_param_logical, doc_param_logical_array, &
                   doc_param_int,     doc_param_int_array, &
                   doc_param_real,    doc_param_real_array, &
                   doc_param_char, &
                   doc_param_time
end interface

integer, parameter :: mLen = 1240 !< Length of interface/message strings

!> A structure that controls where the documentation occurs, its veborsity and formatting.
type, public :: doc_type ; private
  integer :: unitAll = -1           !< The open unit number for docFileBase + .all.
  integer :: unitShort = -1         !< The open unit number for docFileBase + .short.
  integer :: unitLayout = -1        !< The open unit number for docFileBase + .layout.
  integer :: unitDebugging  = -1    !< The open unit number for docFileBase + .debugging.
  logical :: filesAreOpen = .false. !< True if any files were successfully opened.
  character(len=mLen) :: docFileBase = '' !< The basename of the files where run-time
                                    !! parameters, settings and defaults are documented.
  logical :: complete = .true.      !< If true, document all parameters.
  logical :: minimal = .true.       !< If true, document non-default parameters.
  logical :: layout = .true.        !< If true, document layout parameters.
  logical :: debugging = .true.     !< If true, document debugging parameters.
  logical :: defineSyntax = .false. !< If true, use '\#def' syntax instead of a=b syntax
  logical :: warnOnConflicts = .false. !< Cause a WARNING error if defaults differ.
  integer :: commentColumn = 32     !< Number of spaces before the comment marker.
  integer :: max_line_len = 112     !< The maximum length of message lines.
  type(link_msg), pointer :: chain_msg => NULL() !< Database of messages
  character(len=240) :: blockPrefix = '' !< The full name of the current block.
end type doc_type

!> A linked list of the parameter documentation messages that have been issued so far.
type :: link_msg ; private
  type(link_msg), pointer :: next => NULL()  !< Facilitates linked list
  character(len=80) :: name                  !< Parameter name
  character(len=620) :: msg                  !< Parameter value and default
end type link_msg

character(len=4), parameter :: STRING_TRUE  = 'True'  !< A string for true logicals
character(len=5), parameter :: STRING_FALSE = 'False' !< A string for false logicals


  interface
module subroutine doc_param_none(doc, varname, desc, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varname !< The name of the parameter being documented
  character(len=*), intent(in) :: desc    !< A description of the parameter being documented
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
! This subroutine handles parameter documentation with no value.

end subroutine doc_param_none
module subroutine doc_param_logical(doc, varname, desc, units, val, default, &
                             layoutParam, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  logical,           intent(in) :: val     !< The value of this parameter
  logical, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for logicals.

end subroutine doc_param_logical
module subroutine doc_param_logical_array(doc, varname, desc, units, vals, default, &
                                   layoutParam, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  logical,           intent(in) :: vals(:) !< The array of values to record
  logical, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for arrays of logicals.

end subroutine doc_param_logical_array
module subroutine doc_param_int(doc, varname, desc, units, val, default, &
                         layoutParam, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  integer,           intent(in) :: val     !< The value of this parameter
  integer, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for integers.

end subroutine doc_param_int
module subroutine doc_param_int_array(doc, varname, desc, units, vals, default, defaults, &
                               layoutParam, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  integer,           intent(in) :: vals(:) !< The array of values to record
  integer, optional, intent(in) :: default !< The uniform default value of this parameter
  integer, optional, intent(in) :: defaults(:) !< The element-wise default values of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for arrays of integers.

end subroutine doc_param_int_array
module subroutine doc_param_real(doc, varname, desc, units, val, default, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  real,              intent(in) :: val     !< The value of this parameter
  real,    optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for reals.

end subroutine doc_param_real
module subroutine doc_param_real_array(doc, varname, desc, units, vals, default, defaults, &
                                debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  real,              intent(in) :: vals(:) !< The array of values to record
  real,    optional, intent(in) :: default !< A uniform default value of this parameter
  real,    optional, intent(in) :: defaults(:) !< The element-wise default values of this parameter
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for arrays of reals.

end subroutine doc_param_real_array
module subroutine doc_param_char(doc, varname, desc, units, val, default, &
                          layoutParam, debuggingParam, like_default)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  character(len=*),  intent(in) :: val     !< The value of the parameter
  character(len=*), &
           optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical, optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                           !! it has the default value, even if there is no default.
! This subroutine handles parameter documentation for character strings.

end subroutine doc_param_char
module subroutine doc_openBlock(doc, blockName, desc)
  type(doc_type),   pointer    :: doc       !< A pointer to a structure that controls where the
                                            !! documentation occurs and its formatting
  character(len=*), intent(in) :: blockName !< The name of the parameter block being opened
  character(len=*), optional, intent(in) :: desc !< A description of the parameter block being opened
! This subroutine handles documentation for opening a parameter block.

end subroutine doc_openBlock
module subroutine doc_closeBlock(doc, blockName)
  type(doc_type),   pointer    :: doc !< A pointer to a structure that controls where the
                                      !! documentation occurs and its formatting
  character(len=*), intent(in) :: blockName !< The name of the parameter block being closed
! This subroutine handles documentation for closing a parameter block.

end subroutine doc_closeBlock
module subroutine doc_param_time(doc, varname, desc, val, default, units, debuggingParam, like_default)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varname !< The name of the parameter being documented
  character(len=*), intent(in) :: desc    !< A description of the parameter being documented
  type(time_type),  intent(in) :: val     !< The value of the parameter
  type(time_type),  optional, intent(in) :: default !< The default value of this parameter
  character(len=*), optional, intent(in) :: units   !< The units of the parameter being documented
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as though
                                             !! it has the default value, even if there is no default.

  ! Local varables

end subroutine doc_param_time
module subroutine writeMessageAndDesc(doc, vmesg, desc, valueWasDefault, indent, &
                               layoutParam, debuggingParam)
  type(doc_type),    intent(in) :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: vmesg   !< A message with the parameter name, units, and default value.
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  logical, optional, intent(in) :: valueWasDefault !< If true, this parameter has its default value
  integer, optional, intent(in) :: indent      !< An amount by which to indent this message
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.

  ! Local variables

end subroutine writeMessageAndDesc
module function time_string(time)
  type(time_type), intent(in) :: time !< The time type being translated
  character(len=40) :: time_string

  ! Local variables

end function time_string
module function real_string(val)
  real, intent(in)  :: val !< The value being written into a string
  character(len=32) :: real_string
! This function returns a string with a real formatted like '(G)'

end function real_string
module function real_array_string(vals, sep)
  character(len=:) ,allocatable :: real_array_string !< The output string listing vals
  real,      intent(in)  :: vals(:) !< The array of values to record
  character(len=*), &
    optional, intent(in) :: sep     !< The separator between successive values,
                                    !! by default it is ', '.
! Returns a character string of a comma-separated, compact formatted, reals
! e.g. "1., 2., 5*3., 5.E2"
  ! Local variables
end function real_array_string
module function int_array_string(vals, sep)
  character(len=:), allocatable :: int_array_string !< The output string listing vals
  integer,          intent(in)  :: vals(:) !< The array of values to record
  character(len=*), &
           optional, intent(in) :: sep     !< The separator between successive values,
                                           !! by default it is ', '.

  ! Local variables
end function int_array_string
module function testFormattedFloatIsReal(str, val)
  character(len=*), intent(in) :: str !< The string that match val
  real,             intent(in) :: val !< The value being tested
  logical                      :: testFormattedFloatIsReal
  ! Local variables

end function testFormattedFloatIsReal
module function int_string(val)
  integer, intent(in)  :: val !< The value being written into a string
  character(len=24)    :: int_string
! This function returns a string with an integer formatted like '(I)'
end function int_string
module function logical_string(val)
  logical, intent(in)  :: val !< The value being written into a string
  character(len=24)    :: logical_string
! This function returns a string with an logical formatted like '(L)'
end function logical_string
module function define_string(doc, varName, valString, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varName !< The name of the parameter being documented
  character(len=*), intent(in) :: valString !< A string containing the value of the parameter
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
  character(len=mLen) :: define_string
! This function returns a string for formatted parameter assignment
end function define_string
module function undef_string(doc, varName, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varName !< The name of the parameter being documented
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
  character(len=mLen) :: undef_string
! This function returns a string for formatted false logicals
end function undef_string
module subroutine doc_module(doc, modname, desc, log_to_all, all_default, layoutMod, debuggingMod)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: desc    !< A description of the module being documented
  logical, optional, intent(in) :: log_to_all !< If present and true, log this parameter to the
                                          !! ..._doc.all files, even if this module also has layout
                                          !! or debugging parameters.
  logical, optional, intent(in) :: all_default  !< If true, all parameters take their default values.
  logical, optional, intent(in) :: layoutMod    !< If present and true, this module has layout parameters.
  logical, optional, intent(in) :: debuggingMod !< If present and true, this module has debugging parameters.

  ! This subroutine handles the module documentation

end subroutine doc_module
module subroutine doc_subroutine(doc, modname, subname, desc)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: subname !< The name of the subroutine being documented
  character(len=*), intent(in) :: desc    !< A description of the subroutine being documented
! This subroutine handles the subroutine documentation
end subroutine doc_subroutine
module subroutine doc_function(doc, modname, fnname, desc)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: fnname  !< The name of the function being documented
  character(len=*), intent(in) :: desc    !< A description of the function being documented
! This subroutine handles the function documentation
end subroutine doc_function
module subroutine doc_init(docFileBase, doc, minimal, complete, layout, debugging)
  character(len=*),  intent(in)  :: docFileBase !< The base file name for this set of parameters,
                                             !! for example MOM_parameter_doc
  type(doc_type),    pointer     :: doc      !< A pointer to a structure that controls where the
                                             !! documentation occurs and its formatting
  logical, optional, intent(in)  :: minimal  !< If present and true, write out the files (.short) documenting
                                             !! those parameters that do not take on their default values.
  logical, optional, intent(in)  :: complete !< If present and true, write out the (.all) files documenting all
                                             !! parameters
  logical, optional, intent(in)  :: layout   !< If present and true, write out the (.layout) files documenting
                                             !! the layout parameters
  logical, optional, intent(in)  :: debugging !< If present and true, write out the (.debugging) files documenting
                                             !! the debugging parameters

end subroutine doc_init
module subroutine open_doc_file(doc)
  type(doc_type), pointer :: doc !< A pointer to a structure that controls where the
                                 !! documentation occurs and its formatting


end subroutine open_doc_file
module function find_unused_unit_number()
! Find an unused unit number.
! Returns >0 if found. FATAL if not.
  integer :: find_unused_unit_number
end function find_unused_unit_number
module subroutine doc_end(doc)
  type(doc_type), pointer :: doc !< A pointer to a structure that controls where the
                                 !! documentation occurs and its formatting

end subroutine doc_end
module function mesgHasBeenDocumented(doc,varName,mesg)
  type(doc_type),   pointer     :: doc  !< A pointer to a structure that controls where the
                                        !! documentation occurs and its formatting
  character(len=*), intent(in)  :: varName !< The name of the parameter being documented
  character(len=*), intent(in)  :: mesg !< A message with parameter values, defaults, and descriptions
                                        !! to compare with the message that was written previously
  logical                       :: mesgHasBeenDocumented
! Returns true if documentation has already been written

end function mesgHasBeenDocumented
  end interface

end module MOM_document
