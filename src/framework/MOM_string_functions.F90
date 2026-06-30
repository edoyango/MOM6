! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Handy functions for manipulating strings
module MOM_string_functions

use iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public lowercase, uppercase
public left_int, left_ints
public left_real, left_reals
public string_functions_unit_tests
public extractWord
public extract_word
public extract_integer
public extract_real
public remove_spaces
public slasher
public trim_trailing_commas
public ints_to_string


  interface
module function lowercase(input_string)
  character(len=*),     intent(in) :: input_string !< The string to modify
  character(len=len(input_string)) :: lowercase !< The modified output string
!   This function returns a string in which all uppercase letters have been
! replaced by their lowercase counterparts.  It is loosely based on the
! lowercase function in mpp_util.F90.

end function lowercase
module function uppercase(input_string)
  character(len=*),     intent(in) :: input_string !< The string to modify
  character(len=len(input_string)) :: uppercase !< The modified output string
!   This function returns a string in which all lowercase letters have been
! replaced by their uppercase counterparts.  It is loosely based on the
! uppercase function in mpp_util.F90.

end function uppercase
module function left_int(i)
  integer, intent(in) :: i !< The integer to convert to a string
  character(len=19) :: left_int !< The output string

end function left_int
module function left_ints(i)
  integer, intent(in) :: i(:) !< The array of integers to convert to a string
  character(len=1320) :: left_ints !< The output string

end function left_ints
module function left_real(val)
  real, intent(in)  :: val !< The real variable to convert to a string, in arbitrary units [A]
  character(len=32) :: left_real !< The output string


end function left_real
module function left_reals(r,sep)
  real, intent(in) :: r(:) !< The array of real variables to convert to a string, in arbitrary units [A]
  character(len=*), optional, intent(in) :: sep !< The separator between
                                    !! successive values, by default it is ', '.
  character(len=:), allocatable :: left_reals !< The output string


end function left_reals
module function isFormattedFloatEqualTo(str, val)
  character(len=*), intent(in) :: str !< The string to parse
  real,             intent(in) :: val !< The real value to compare with, in arbitrary units [A]
  logical                      :: isFormattedFloatEqualTo
  ! Local variables

end function isFormattedFloatEqualTo
character(len=120) module function extractWord(string, n)
  character(len=*),   intent(in) :: string !< The string to scan
  integer,            intent(in) :: n      !< Number of word to extract

end function extractWord
character(len=120) module function extract_word(string, separators, n)
  character(len=*),   intent(in) :: string     !< String to scan
  character(len=*),   intent(in) :: separators !< Characters to use for delineation
  integer,            intent(in) :: n          !< Number of word to extract
  ! Local variables
end function extract_word
integer module function extract_integer(string, separators, n, missing_value)
  character(len=*),   intent(in) :: string     !< String to scan
  character(len=*),   intent(in) :: separators !< Characters to use for delineation
  integer,            intent(in) :: n          !< Number of word to extract
  integer, optional,  intent(in) :: missing_value !< Value to assign if word is missing
  ! Local variables

end function extract_integer
real module function extract_real(string, separators, n, missing_value)
  character(len=*), intent(in) :: string     !< String to scan
  character(len=*), intent(in) :: separators !< Characters to use for delineation
  integer,          intent(in) :: n          !< Number of word to extract
  real, optional,   intent(in) :: missing_value !< Value to assign if word is missing, in arbitrary units [A]
  ! Local variables

end function extract_real
character(len=120) module function remove_spaces(string)
  character(len=*),   intent(in) :: string     !< String to scan
  ! Local variables
end function remove_spaces
logical module function string_functions_unit_tests(verbose)
  ! Arguments
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  ! This is an array of real test values, in arbitrary units [A]
end function string_functions_unit_tests
logical module function localTestS(verbose,str1,str2)
  logical, intent(in) :: verbose !< If true, write results to stdout
  character(len=*), intent(in) :: str1 !< String
  character(len=*), intent(in) :: str2 !< String
end function localTestS
logical module function localTestI(verbose,i1,i2)
  logical, intent(in) :: verbose !< If true, write results to stdout
  integer, intent(in) :: i1 !< Integer
  integer, intent(in) :: i2 !< Integer
end function localTestI
logical module function localTestR(verbose,r1,r2)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: r1 !< The first value to compare, in arbitrary units [A]
  real, intent(in) :: r2 !< The first value to compare, in arbitrary units [A]
end function localTestR
module function slasher(dir)
  character(len=*), intent(in) :: dir !< A directory to be terminated with a "/"
                                      !! or changed to "./" if it is blank.
  character(len=len(dir)+2) :: slasher

end function slasher
module function trim_trailing_commas(in_str) result(out_str)
  character(len=*), intent(in) :: in_str  !< A string that is to be left adjusted and have
                                          !! its trailing commas and white space removed.
  character(len=len(in_str))   :: out_str !< A left-adjusted version of in_str with
                                          !! trailing commas and white space removed

end function trim_trailing_commas
module function ints_to_string(a, n) result(i2s)
  integer, dimension(:), intent(in) :: a !< The array of integers to translate
  integer, optional    , intent(in) :: n !< The number of elements to translate, by default the lesser
                                         !! of 3 or all of the integers
  character(len=5*size(a)+1) :: i2s !< The returned underscore delimited string of integers


end function ints_to_string
  end interface

end module MOM_string_functions
