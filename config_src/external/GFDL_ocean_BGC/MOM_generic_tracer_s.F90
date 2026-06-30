submodule (MOM_generic_tracer) MOM_generic_tracer_s
#include <MOM_memory.h>
#ifdef NO_F2000
#define _ALLOCATED associated
#else
#define _ALLOCATED allocated
#endif
  implicit none
contains
module procedure register_MOM_generic_tracer
  register_MOM_generic_tracer = .false.

  call MOM_error(FATAL, "register_MOM_generic_tracer should not be called with the stub code "// &
         "in MOM6/config_src/external, as it does nothing.  Recompile using the full MOM_generic_tracer package.")

end procedure register_MOM_generic_tracer
module procedure register_MOM_generic_tracer_segments
end procedure register_MOM_generic_tracer_segments
module procedure initialize_MOM_generic_tracer
end procedure initialize_MOM_generic_tracer
module procedure MOM_generic_tracer_column_physics
end procedure MOM_generic_tracer_column_physics
module procedure MOM_generic_tracer_stock
  integer :: m
  MOM_generic_tracer_stock = 0

  ! These should never be used, but they are set to avoid compile-time warnings
  do m=1,size(names) ; names(m) = "" ; enddo
  do m=1,size(units) ; units(m) = "" ; enddo
  do m=1,size(stocks) ; stocks(m) = real_to_EFP(0.0) ; enddo

end procedure MOM_generic_tracer_stock
module procedure MOM_generic_tracer_min_max
  integer :: m
  MOM_generic_tracer_min_max = 0

  ! These should never be used, but they are set to avoid compile-time warnings.  Note that the minimum values
  ! are delibarately set to be larger than the maximum values.
  got_minmax(:) = .false.
  gmax(:) = -huge(gmax)
  gmin(:) = huge(gmin)
  do m=1,size(names) ; names(m) = "" ; enddo
  do m=1,size(units) ; units(m) = "" ; enddo
  if (present(xgmin)) xgmin(:) = 0.0
  if (present(ygmin)) ygmin(:) = 0.0
  if (present(zgmin)) zgmin(:) = 0.0
  if (present(xgmax)) xgmax(:) = 0.0
  if (present(ygmax)) ygmax(:) = 0.0
  if (present(zgmax)) zgmax(:) = 0.0

end procedure MOM_generic_tracer_min_max
module procedure MOM_generic_tracer_surface_state
end procedure MOM_generic_tracer_surface_state
module procedure MOM_generic_flux_init
end procedure MOM_generic_flux_init
module procedure MOM_generic_tracer_fluxes_accumulate
end procedure MOM_generic_tracer_fluxes_accumulate
module procedure MOM_generic_tracer_get
  real, dimension(:,:,:),   pointer :: array_ptr  ! The tracer in the generic tracer structures, in
  character(len=128), parameter :: sub_name = 'MOM_generic_tracer_get'
  array(:,:,:) = huge(array)

end procedure MOM_generic_tracer_get
module procedure end_MOM_generic_tracer
end procedure end_MOM_generic_tracer
end submodule MOM_generic_tracer_s
