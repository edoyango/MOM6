submodule (MOM_get_input) MOM_get_input_s
  implicit none
contains
module procedure get_MOM_input
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=240) :: &
    parameter_filename(npf), & ! List of files containing parameters.
    output_directory,        & ! Directory to use to write the model output.
    restart_input_dir,       & ! Directory for reading restart and input files.
    restart_output_dir         ! Directory into which to write restart files.
  character(len=2048) :: &
    input_filename             ! A string that indicates the input files or how
  character(len=240) :: output_dir
  integer :: unit, io, ierr, valid_param_files
  type(stat_buf) :: buf
  namelist /MOM_input_nml/ output_directory, input_filename, parameter_filename, &
                           restart_input_dir, restart_output_dir
  parameter_filename(:) = ' '
  output_directory = ' '
  restart_input_dir = ' '
  restart_output_dir = ' '
  input_filename  = ' '
  if (present(default_input_filename)) input_filename = trim(default_input_filename)

  ! Open namelist
  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call MOM_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ! Read namelist parameters
  ! NOTE: Every rank is reading MOM_input_nml
  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=MOM_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'MOM_input_nml')
  enddo
10 call close_file(unit)

  ! Store parameters in container
  if (present(dirs)) then
    if (present(ensemble_num)) then
      dirs%output_directory = slasher(ensembler(output_directory,ensemble_num))
      dirs%restart_output_dir = slasher(ensembler(restart_output_dir,ensemble_num))
      dirs%restart_input_dir = slasher(ensembler(restart_input_dir,ensemble_num))
      dirs%input_filename = ensembler(input_filename,ensemble_num)
    else
      dirs%output_directory = slasher(ensembler(output_directory))
      dirs%restart_output_dir = slasher(ensembler(restart_output_dir))
      dirs%restart_input_dir = slasher(ensembler(restart_input_dir))
      dirs%input_filename = ensembler(input_filename)
    endif

    ! Create the RESTART directory if absent
    if (is_root_PE()) then
      if (stat(trim(dirs%restart_output_dir), buf) == -1) then
        ierr = mkdir(trim(dirs%restart_output_dir), int(o'700'))
        if (ierr == -1) &
          call MOM_error(FATAL, 'Restart directory could not be created.')
      endif
    endif
  endif

  ! Open run-time parameter file(s)
  if (present(param_file)) then
    output_dir = slasher(ensembler(output_directory))
    valid_param_files = 0
    do io = 1, npf
      if (len_trim(trim(parameter_filename(io))) > 0) then
        if (present(ensemble_num)) then
          call open_param_file(ensembler(parameter_filename(io),ensemble_num), param_file, &
               check_params, doc_file_dir=output_dir, ensemble_num=ensemble_num)
        else
          call open_param_file(ensembler(parameter_filename(io)), param_file, &
               check_params, doc_file_dir=output_dir)
        endif
        valid_param_files = valid_param_files + 1
      endif
    enddo
    if (valid_param_files == 0) call MOM_error(FATAL, "There must be at "//&
         "least 1 valid entry in input_filename in MOM_input_nml in input.nml.")
  endif

end procedure get_MOM_input
end submodule MOM_get_input_s
