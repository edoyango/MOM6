submodule (MOM_tracer_advect_schemes) MOM_tracer_advect_schemes_s
  implicit none
contains
module procedure set_tracer_advect_scheme
  select case (trim(advect_scheme_name))
    case ("")
      scheme_value = -1
    case ("PLM")
      scheme_value = ADVECT_PLM
    case ("PPM:H3")
      scheme_value = ADVECT_PPMH3
    case ("PPM")
      scheme_value = ADVECT_PPM
    case default
      call MOM_error(FATAL, "set_tracer_advect_scheme: "//&
           "Unknown TRACER_ADVECTION_SCHEME = "//trim(advect_scheme_name))
  end select
end procedure set_tracer_advect_scheme
end submodule MOM_tracer_advect_schemes_s
