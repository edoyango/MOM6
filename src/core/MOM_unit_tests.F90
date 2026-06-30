! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Invokes unit tests in all modules that have them
module MOM_unit_tests

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform,            only : symmetric_sum_unit_tests
use MOM_diag_buffers,               only : diag_buffer_unit_tests_2d, diag_buffer_unit_tests_3d
use MOM_error_handler,              only : MOM_error, FATAL, is_root_pe
use MOM_hor_bnd_diffusion,          only : near_boundary_unit_tests
use MOM_intrinsic_functions,        only : intrinsic_functions_unit_tests
use MOM_mixed_layer_restrat,        only : mixedlayer_restrat_unit_tests
use MOM_neutral_diffusion,          only : neutral_diffusion_unit_tests
use MOM_random,                     only : random_unit_tests
use MOM_remapping,                  only : remapping_unit_tests
use MOM_string_functions,           only : string_functions_unit_tests
use MOM_CFC_cap,                    only : CFC_cap_unit_tests
use MOM_EOS,                        only : EOS_unit_tests

implicit none ; private

public unit_tests


  interface
module subroutine unit_tests(verbosity)
  ! Arguments
  integer, intent(in) :: verbosity !< The verbosity level
  ! Local variables

end subroutine unit_tests
  end interface

end module MOM_unit_tests
