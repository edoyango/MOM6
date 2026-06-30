submodule (marbl_interface) marbl_interface_s
  implicit none
contains
    module procedure put_setting
        call MOM_error(FATAL, error_msg)
    end procedure put_setting
    module procedure get_setting
      log_out = .false.
      call MOM_error(FATAL, error_msg)
    end procedure get_setting
    module procedure init
        call MOM_error(FATAL, error_msg)
    end procedure init
    module procedure compute_totChl
      call MOM_error(FATAL, error_msg)

    end procedure compute_totChl
    module procedure surface_flux_compute
        call MOM_error(FATAL, error_msg)

    end procedure surface_flux_compute
    module procedure interior_tendency_compute
        call MOM_error(FATAL, error_msg)

    end procedure interior_tendency_compute
    module procedure add_output_for_GCM
        output_id = 0
        field_source = ""

    end procedure add_output_for_GCM
    module procedure shutdown
        call MOM_error(FATAL, error_msg)

    end procedure shutdown
end submodule marbl_interface_s
