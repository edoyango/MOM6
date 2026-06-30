submodule (MOM_unit_tests) MOM_unit_tests_s
  implicit none
contains
module procedure unit_tests
  logical :: verbose
  verbose = verbosity>=5

  if (is_root_pe()) then ! The following need only be tested on 1 PE
    if (string_functions_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: string_functions_unit_tests FAILED")
    if (symmetric_sum_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: symmetric_sum_unit_tests FAILED")
    if (EOS_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: EOS_unit_tests FAILED")
    if (remapping_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: remapping_unit_tests FAILED")
    if (intrinsic_functions_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: intrinsic_functions_unit_tests FAILED")
    if (neutral_diffusion_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: neutralDiffusionUnitTests FAILED")
    if (random_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: random_unit_tests FAILED")
    if (near_boundary_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: near_boundary_unit_tests FAILED")
    if (CFC_cap_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: CFC_cap_unit_tests FAILED")
    if (mixedlayer_restrat_unit_tests(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: mixedlayer_restrat_unit_tests FAILED")
    if (diag_buffer_unit_tests_2d(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: diag_buffer_unit_tests_2d FAILED")
    if (diag_buffer_unit_tests_3d(verbose)) call MOM_error(FATAL, &
       "MOM_unit_tests: diag_buffer_unit_tests_3d FAILED")
  endif

end procedure unit_tests
end submodule MOM_unit_tests_s
