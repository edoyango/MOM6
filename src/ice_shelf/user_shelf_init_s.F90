submodule (user_shelf_init) user_shelf_init_s
#include <MOM_memory.h>
  implicit none
contains
module procedure USER_initialize_shelf_mass
  character(len=40) :: mdl = "USER_initialize_shelf_mass" ! This subroutine's name.
  if (.not.associated(CS)) allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  if (CS%first_call) call write_user_log(param_file)
  CS%first_call = .false.
  call get_param(param_file, mdl, "RHO_0", CS%Rho_ocean, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "SHELF_MAX_DRAFT", CS%max_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "SHELF_MIN_DRAFT", CS%min_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "FLAT_SHELF_WIDTH", CS%flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_SLOPE_SCALE", CS%shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_EDGE_POS_0", CS%pos_shelf_edge_0, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mdl, "SHELF_SPEED", CS%shelf_speed, &
                 units="axis_units day-1", default=0.0)

  call USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, set_time(0,0), new_sim)

end procedure USER_initialize_shelf_mass
module procedure USER_init_ice_thickness
  real, dimension(SZI_(G),SZJ_(G)) :: mass_shelf ! The ice shelf mass per unit area averaged
  type(user_ice_shelf_CS), pointer :: CS => NULL()
  call USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, US, CS, param_file, .true.)

end procedure USER_init_ice_thickness
module procedure USER_update_shelf_mass
  real :: c1        ! The inverse of the range over which the shelf slopes [km-1]
  real :: edge_pos  ! The time-evolving position the ice shelf edge [km]
  real :: slope_pos ! The time-evolving position of the start of the ice shelf slope [km]
  integer :: i, j
  edge_pos = CS%pos_shelf_edge_0 + CS%shelf_speed*(time_type_to_real(Time) / 86400.0)

  slope_pos = edge_pos - CS%flat_shelf_width
  c1 = 0.0 ; if (CS%shelf_slope_scale > 0.0) c1 = 1.0 / CS%shelf_slope_scale


  do j=G%jsd,G%jed

    if (((j+G%jdg_offset) <= G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+G%jdg_offset) >= G%domain%njhalo+1)) then

      do i=G%isc,G%iec

!    if (((i+G%idg_offset) <= G%domain%niglobal+G%domain%nihalo) .AND. &
!           ((i+G%idg_offset) >= G%domain%nihalo+1)) then

        if ((j >= G%jsc) .and. (j <= G%jec)) then
          if (new_sim) then ; if (G%geoLonCu(i-1,j) >= edge_pos) then
            ! Everything past the edge is open ocean.
            mass_shelf(i,j) = 0.0
            area_shelf_h(i,j) = 0.0
            hmask (i,j) = 0.0
            h_shelf (i,j) = 0.0
          else
            if (G%geoLonCu(i,j) > edge_pos) then
              area_shelf_h(i,j) = G%areaT(i,j) * (edge_pos - G%geoLonCu(i-1,j)) / &
                                  (G%geoLonCu(i,j) - G%geoLonCu(i-1,j))
              hmask (i,j) = 2.0
            else
              area_shelf_h(i,j) = G%areaT(i,j)
              hmask (i,j) = 1.0
            endif

            if (G%geoLonT(i,j) > slope_pos) then
              h_shelf (i,j) = CS%min_draft
              mass_shelf(i,j) = CS%Rho_ocean * CS%min_draft
            else
              mass_shelf(i,j) = CS%Rho_ocean * (CS%min_draft + &
                     (CS%max_draft - CS%min_draft) * &
                     min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
              h_shelf(i,j) = (CS%min_draft + &
                     (CS%max_draft - CS%min_draft) * &
                     min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
            endif
          endif ; endif
        endif

        if ((i+G%idg_offset) == G%domain%nihalo+1) then
          hmask(i-1,j) = 3.0
        endif

      enddo
    endif
  enddo

end procedure USER_update_shelf_mass
module procedure write_user_log
  character(len=128) :: version = '$Id: user_shelf_init.F90,v 1.1.2.7 2012/06/19 22:15:52 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: MOM_ogrp $'
  character(len=40)  :: mdl = "user_shelf_init" ! This module's name.
  call log_version(param_file, mdl, version, tagname)

end procedure write_user_log
end submodule user_shelf_init_s
