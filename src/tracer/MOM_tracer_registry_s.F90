submodule (MOM_tracer_registry) MOM_tracer_registry_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_tracer
  logical :: mand
  type(tracer_type), pointer :: Tr=>NULL()
  character(len=256) :: mesg    ! Message for error messages.
  if (.not. associated(Reg)) call tracer_registry_init(param_file, Reg)

  if (Reg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I0," to allow for &
        &all the tracers being registered via register_tracer.")') Reg%ntr+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  Reg%ntr = Reg%ntr + 1

  Tr => Reg%Tr(Reg%ntr)
  if (present(Tr_out)) Tr_out => Reg%Tr(Reg%ntr)

  if (present(name)) then
    Tr%name = name
    Tr%longname = name ; if (present(longname)) Tr%longname = longname
    Tr%units = "Conc" ; if (present(units)) Tr%units = units

    Tr%cmor_name = ""
    if (present(cmor_name)) Tr%cmor_name = cmor_name

    Tr%cmor_units = Tr%units
    if (present(cmor_units)) Tr%cmor_units = cmor_units

    Tr%cmor_longname = ""
    if (present(cmor_longname)) Tr%cmor_longname = cmor_longname

    if (present(tr_desc)) call MOM_error(WARNING, "MOM register_tracer: "//&
      "It is a bad idea to use both name and tr_desc when registring "//trim(name))
  elseif (present(tr_desc)) then
    call query_vardesc(tr_desc, name=Tr%name, units=Tr%units, &
                       longname=Tr%longname, cmor_field_name=Tr%cmor_name, &
                       cmor_longname=Tr%cmor_longname, caller="register_tracer")
    Tr%cmor_units = Tr%units
  else
    call MOM_error(FATAL,"MOM register_tracer: Either name or "//&
                   "tr_desc must be present when registering a tracer.")
  endif

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(Tr%name)//&
      " with a locked tracer registry.")

  Tr%conc_scale = 1.0
  if (present(conc_scale)) Tr%conc_scale = conc_scale

  Tr%conc_underflow = 0.0
  if (present(underflow_conc)) Tr%conc_underflow = underflow_conc

  Tr%flux_nameroot = Tr%name
  if (present(flux_nameroot)) then
    if (len_trim(flux_nameroot) > 0) Tr%flux_nameroot = flux_nameroot
  endif

  Tr%flux_longname = Tr%longname
  if (present(flux_longname)) then
    if (len_trim(flux_longname) > 0) Tr%flux_longname = flux_longname
  endif

  Tr%net_surfflux_name = "KPP_net"//trim(Tr%name)
  if (present(net_surfflux_name)) then
    Tr%net_surfflux_name = net_surfflux_name
  endif

  Tr%NLT_budget_name = 'KPP_NLT_'//trim(Tr%flux_nameroot)//'_budget'
  if (present(NLT_budget_name)) then
    Tr%NLT_budget_name = NLT_budget_name
  endif

  Tr%net_surfflux_longname = 'Effective net surface '//trim(lowercase(Tr%flux_longname))//&
                             ' flux, as used by [CVMix] KPP'
  if (present(net_surfflux_longname)) then
    Tr%net_surfflux_longname = net_surfflux_longname
  endif

  Tr%flux_units = ""
  if (present(flux_units)) Tr%flux_units = flux_units

  Tr%flux_scale = GV%H_to_MKS*Tr%conc_scale
  if (present(flux_scale)) Tr%flux_scale = flux_scale

  Tr%conv_units = ""
  if (present(convergence_units)) Tr%conv_units = convergence_units

  Tr%cmor_tendprefix = ""
  if (present(cmor_tendprefix)) Tr%cmor_tendprefix = cmor_tendprefix

  Tr%conv_scale = GV%H_to_MKS*Tr%conc_scale
  if (present(convergence_scale)) then
    Tr%conv_scale = convergence_scale
  elseif (present(flux_scale)) then
    Tr%conv_scale = flux_scale
  endif

  Tr%diag_form = 1
  if (present(diag_form)) Tr%diag_form = diag_form

  Tr%advect_scheme = -1
  if (present(advect_scheme)) Tr%advect_scheme = advect_scheme

  Tr%t => tr_ptr

  if (present(registry_diags)) Tr%registry_diags = registry_diags

  if (present(ad_x)) then ; if (associated(ad_x)) Tr%ad_x => ad_x ; endif
  if (present(ad_y)) then ; if (associated(ad_y)) Tr%ad_y => ad_y ; endif
  if (present(df_x)) then ; if (associated(df_x)) Tr%df_x => df_x ; endif
  if (present(df_y)) then ; if (associated(df_y)) Tr%df_y => df_y ; endif
! if (present(OBC_inflow)) Tr%OBC_inflow_conc = OBC_inflow
! if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) Tr%OBC_in_u => OBC_in_u ; endif
! if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) Tr%OBC_in_v => OBC_in_v ; endif
  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Tr%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Tr%ad2d_y => ad_2d_y ; endif
  if (present(df_2d_x)) then ; if (associated(df_2d_x)) Tr%df2d_x => df_2d_x ; endif
  if (present(df_2d_y)) then ; if (associated(df_2d_y)) Tr%df2d_y => df_2d_y ; endif

  if (present(advection_xy)) then
    if (associated(advection_xy)) Tr%advection_xy => advection_xy
  endif

  if (present(restart_CS)) then
    ! Register this tracer to be read from and written to restart files.
    mand = .true. ; if (present(mandatory)) mand = mandatory

    call register_restart_field(tr_ptr, Tr%name, mand, restart_CS, &
                                longname=Tr%longname, units=Tr%units, conversion=conc_scale)
  endif
end procedure register_tracer
module procedure lock_tracer_registry
  if (.not. associated(Reg)) call MOM_error(WARNING, &
    "lock_tracer_registry called with an unassociated registry.")

  Reg%locked = .True.

end procedure lock_tracer_registry
module procedure register_tracer_diagnostics
  character(len=24)  :: name     ! A variable's name in a NetCDF file.
  character(len=24)  :: shortnm  ! A shortened version of a variable's name for
  character(len=72)  :: longname ! The long name of that tracer variable.
  character(len=72)  :: flux_longname ! The tracer name in the long names of fluxes.
  character(len=48)  :: units    ! The dimensions of the tracer.
  character(len=48)  :: flux_units ! The units for fluxes, either
  character(len=48)  :: conv_units ! The units for flux convergences, either
  character(len=48)  :: unit2    ! The dimensions of the tracer squared
  character(len=72)  :: cmorname ! The CMOR name of this tracer.
  character(len=120) :: cmor_longname ! The CMOR long name of that variable.
  character(len=120) :: var_lname      ! A temporary longname for a diagnostic.
  character(len=120) :: cmor_var_lname ! The temporary CMOR long name for a diagnostic
  real :: conversion ! Temporary term while we address a bug [conc m CU-1 H-1 ~> 1] or [conc kg m-2 CU-1 H-1 ~> 1]
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: i, j, k, is, ie, js, je, nz, m, m2, nTr_in
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not. associated(Reg)) call MOM_error(FATAL, "register_tracer_diagnostics: "//&
       "register_tracer must be called before register_tracer_diagnostics")

  nTr_in = Reg%ntr

  do m=1,nTr_in ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
!    call query_vardesc(Tr%vd, name, units=units, longname=longname, &
!                       cmor_field_name=cmorname, cmor_longname=cmor_longname, &
!                       caller="register_tracer_diagnostics")
    name = Tr%name ; units=adjustl(Tr%units) ; longname = Tr%longname
    cmorname = Tr%cmor_name ; cmor_longname = Tr%cmor_longname
    shortnm = Tr%flux_nameroot
    flux_longname = Tr%flux_longname
    if (len_trim(cmor_longname) == 0) cmor_longname = longname

    if (len_trim(Tr%flux_units) > 0) then ; flux_units = Tr%flux_units
    elseif (GV%Boussinesq) then ; flux_units = trim(units)//" m3 s-1"
    else ; flux_units = trim(units)//" kg s-1" ; endif

    if (len_trim(Tr%conv_units) > 0) then ; conv_units = Tr%conv_units
    elseif (GV%Boussinesq) then ; conv_units = trim(units)//" m s-1"
    else ; conv_units = trim(units)//" kg m-2 s-1" ; endif

    if (len_trim(cmorname) == 0) then
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
          Time, trim(longname), trim(units), conversion=Tr%conc_scale)
    else
      Tr%id_tr = register_diag_field("ocean_model", trim(name), diag%axesTL, &
          Time, trim(longname), trim(units), conversion=Tr%conc_scale, &
          cmor_field_name=cmorname, cmor_long_name=cmor_longname, &
          cmor_units=Tr%cmor_units, cmor_standard_name=cmor_long_std(cmor_longname))
    endif
    Tr%id_tr_post_horzn = register_diag_field("ocean_model", &
        trim(name)//"_post_horzn", diag%axesTL, Time, &
        trim(longname)//" after horizontal transport (advection/diffusion) has occurred", &
        trim(units), conversion=Tr%conc_scale)
    if (Tr%diag_form == 1) then
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, trim(flux_longname)//" advective zonal flux" , &
          trim(flux_units), v_extensive=.true., y_cell_method='sum', &
          conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T)
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, trim(flux_longname)//" advective meridional flux" , &
          trim(flux_units), v_extensive=.true., x_cell_method='sum', &
          conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T)
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_dfx", &
          diag%axesCuL, Time, trim(flux_longname)//" diffusive zonal flux" , &
          trim(flux_units), v_extensive=.true., y_cell_method='sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_dfy", &
          diag%axesCvL, Time, trim(flux_longname)//" diffusive meridional flux" , &
          trim(flux_units), v_extensive=.true., x_cell_method='sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_hbd_dfx = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffx", &
          diag%axesCuL, Time, trim(flux_longname)//" diffusive zonal flux " //&
          "from the horizontal boundary diffusion scheme", trim(flux_units), v_extensive=.true., &
          y_cell_method='sum', conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
      Tr%id_hbd_dfy = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffy", &
          diag%axesCvL, Time, trim(flux_longname)//" diffusive meridional " //&
          "flux from the horizontal boundary diffusion scheme", &
          trim(flux_units), v_extensive=.true., x_cell_method='sum', &
          conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T)
    else
      Tr%id_adx = register_diag_field("ocean_model", trim(shortnm)//"_adx", &
          diag%axesCuL, Time, "Advective (by residual mean) Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, y_cell_method='sum')
      Tr%id_ady = register_diag_field("ocean_model", trim(shortnm)//"_ady", &
          diag%axesCvL, Time, "Advective (by residual mean) Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, x_cell_method='sum')
      Tr%id_dfx = register_diag_field("ocean_model", trim(shortnm)//"_diffx", &
          diag%axesCuL, Time, "Diffusive Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          y_cell_method='sum')
      Tr%id_dfy = register_diag_field("ocean_model", trim(shortnm)//"_diffy", &
          diag%axesCvL, Time, "Diffusive Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          x_cell_method='sum')
      Tr%id_hbd_dfx = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffx", &
          diag%axesCuL, Time, &
          "Horizontal Boundary Diffusive Zonal Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          y_cell_method='sum')
      Tr%id_hbd_dfy = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffy", &
          diag%axesCvL, Time, &
          "Horizontal Boundary Diffusive Meridional Flux of "//trim(flux_longname), &
          flux_units, v_extensive=.true., conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
          x_cell_method='sum')
    endif
    Tr%id_zint = register_diag_field("ocean_model", trim(shortnm)//"_zint", &
        diag%axesT1, Time, &
        "Thickness-weighted integral of " // trim(longname), &
        trim(units) // " m")
    Tr%id_zint_100m = register_diag_field("ocean_model", trim(shortnm)//"_zint_100m", &
        diag%axesT1, Time, &
        "Thickness-weighted integral of "// trim(longname) // " over top 100m", &
        trim(units) // " m")
    Tr%id_surf = register_diag_field("ocean_model", trim(shortnm)//"_SURF", &
        diag%axesT1, Time, "Surface values of "// trim(longname), trim(units))
    if (Tr%id_adx > 0) call safe_alloc_ptr(Tr%ad_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_ady > 0) call safe_alloc_ptr(Tr%ad_y,isd,ied,JsdB,JedB,nz)
    if (Tr%id_dfx > 0) call safe_alloc_ptr(Tr%df_x,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_dfy > 0) call safe_alloc_ptr(Tr%df_y,isd,ied,JsdB,JedB,nz)
    if (Tr%id_hbd_dfx > 0) call safe_alloc_ptr(Tr%hbd_dfx,IsdB,IedB,jsd,jed,nz)
    if (Tr%id_hbd_dfy > 0) call safe_alloc_ptr(Tr%hbd_dfy,isd,ied,JsdB,JedB,nz)

    Tr%id_adx_2d = register_diag_field("ocean_model", trim(shortnm)//"_adx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Advective Zonal Flux of "//trim(flux_longname), &
        flux_units, conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, y_cell_method='sum')
    Tr%id_ady_2d = register_diag_field("ocean_model", trim(shortnm)//"_ady_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Advective Meridional Flux of "//trim(flux_longname), &
        flux_units, conversion=Tr%flux_scale*(US%L_to_m**2)*US%s_to_T, x_cell_method='sum')
    Tr%id_dfx_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffx_2d", &
        diag%axesCu1, Time, &
        "Vertically Integrated Diffusive Zonal Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        y_cell_method='sum')
    Tr%id_dfy_2d = register_diag_field("ocean_model", trim(shortnm)//"_diffy_2d", &
        diag%axesCv1, Time, &
        "Vertically Integrated Diffusive Meridional Flux of "//trim(flux_longname), &
        flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        x_cell_method='sum')
    Tr%id_hbd_dfx_2d = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffx_2d", &
        diag%axesCu1, Time, "Vertically-integrated zonal diffusive flux from the horizontal boundary diffusion "//&
        "scheme for "//trim(flux_longname), flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        y_cell_method='sum')
    Tr%id_hbd_dfy_2d = register_diag_field("ocean_model", trim(shortnm)//"_hbd_diffy_2d", &
        diag%axesCv1, Time, "Vertically-integrated meridional diffusive flux from the horizontal boundary diffusion "//&
        "scheme for "//trim(flux_longname), flux_units, conversion=(US%L_to_m**2)*Tr%flux_scale*US%s_to_T, &
        x_cell_method='sum')

    if (Tr%id_adx_2d > 0) call safe_alloc_ptr(Tr%ad2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_ady_2d > 0) call safe_alloc_ptr(Tr%ad2d_y,isd,ied,JsdB,JedB)
    if (Tr%id_dfx_2d > 0) call safe_alloc_ptr(Tr%df2d_x,IsdB,IedB,jsd,jed)
    if (Tr%id_dfy_2d > 0) call safe_alloc_ptr(Tr%df2d_y,isd,ied,JsdB,JedB)
    if (Tr%id_hbd_dfx_2d > 0) call safe_alloc_ptr(Tr%hbd_dfx_2d,IsdB,IedB,jsd,jed)
    if (Tr%id_hbd_dfy_2d > 0) call safe_alloc_ptr(Tr%hbd_dfy_2d,isd,ied,JsdB,JedB)

    Tr%id_adv_xy = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy", &
        diag%axesTL, Time, &
        'Horizontal convergence of residual mean advective fluxes of '//&
        trim(lowercase(flux_longname)), &
        conv_units, v_extensive=.true., conversion=Tr%conv_scale*US%s_to_T)
    Tr%id_adv_xy_2d = register_diag_field('ocean_model', trim(shortnm)//"_advection_xy_2d", &
        diag%axesT1, Time, &
        'Vertical sum of horizontal convergence of residual mean advective fluxes of '//&
        trim(lowercase(flux_longname)), conv_units, conversion=Tr%conv_scale*US%s_to_T)
    if ((Tr%id_adv_xy > 0) .or. (Tr%id_adv_xy_2d > 0)) &
      call safe_alloc_ptr(Tr%advection_xy,isd,ied,jsd,jed,nz)

    Tr%id_tendency = register_diag_field('ocean_model', trim(shortnm)//'_tendency', &
        diag%axesTL, Time, &
        'Net time tendency for '//trim(lowercase(longname)), &
        trim(units)//' s-1', conversion=Tr%conc_scale*US%s_to_T)

    if (Tr%id_tendency > 0) then
      call safe_alloc_ptr(Tr%t_prev,isd,ied,jsd,jed,nz)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Tr%t_prev(i,j,k) = Tr%t(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! Neutral/Horizontal diffusion convergence tendencies
    if (Tr%diag_form == 1) then
      Tr%id_dfxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency', &
          diag%axesTL, Time, "Neutral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, v_extensive=.true.)

      Tr%id_dfxy_cont_2d = register_diag_field("ocean_model", &
          trim(shortnm)//'_dfxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated neutral diffusion tracer content "//&
          "tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T)

      Tr%id_hbdxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_hbdxy_cont_tendency', &
          diag%axesTL, Time, "Horizontal boundary diffusion tracer content tendency for "//&
          trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, v_extensive=.true.)

      Tr%id_hbdxy_cont_2d = register_diag_field("ocean_model", &
          trim(shortnm)//'_hbdxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated horizontal boundary diffusion tracer content "//&
          "tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T)
    else
      cmor_var_lname = 'Tendency of '//trim(lowercase(cmor_longname))//' expressed as '//&
          trim(lowercase(flux_longname))//&
          ' content due to parameterized mesoscale neutral diffusion'
      Tr%id_dfxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_cont_tendency', &
          diag%axesTL, Time, "Neutral diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, v_extensive=.true., &
          cmor_field_name=trim(Tr%cmor_tendprefix)//'pmdiff', &
          cmor_long_name=trim(cmor_var_lname), &
          cmor_standard_name=trim(cmor_long_std(cmor_var_lname)))

      cmor_var_lname = 'Tendency of '//trim(lowercase(cmor_longname))//' expressed as '//&
                       trim(lowercase(flux_longname))//&
                       ' content due to parameterized mesoscale neutral diffusion'
      Tr%id_dfxy_cont_2d = register_diag_field("ocean_model", &
          trim(shortnm)//'_dfxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated neutral diffusion tracer "//&
          "content tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          cmor_field_name=trim(Tr%cmor_tendprefix)//'pmdiff_2d', &
          cmor_long_name=trim(cmor_var_lname), &
          cmor_standard_name=trim(cmor_long_std(cmor_var_lname)))

      Tr%id_hbdxy_cont = register_diag_field("ocean_model", trim(shortnm)//'_hbdxy_cont_tendency', &
          diag%axesTL, Time, &
          "Horizontal boundary diffusion tracer content tendency for "//trim(shortnm), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, v_extensive=.true.)

      Tr%id_hbdxy_cont_2d = register_diag_field("ocean_model", &
          trim(shortnm)//'_hbdxy_cont_tendency_2d', &
          diag%axesT1, Time, "Depth integrated horizontal boundary diffusion of tracer "//&
          "content tendency for "//trim(shortnm), conv_units, conversion=Tr%conv_scale*US%s_to_T)
    endif
    Tr%id_dfxy_conc = register_diag_field("ocean_model", trim(shortnm)//'_dfxy_conc_tendency', &
        diag%axesTL, Time, "Neutral diffusion tracer concentration tendency for "//trim(shortnm), &
        trim(units)//' s-1', conversion=Tr%conc_scale*US%s_to_T)

    Tr%id_hbdxy_conc = register_diag_field("ocean_model", trim(shortnm)//'_hbdxy_conc_tendency', &
        diag%axesTL, Time, &
        "Horizontal diffusion tracer concentration tendency for "//trim(shortnm), &
        trim(units)//' s-1', conversion=Tr%conc_scale*US%s_to_T)

    var_lname = "Net time tendency for "//lowercase(flux_longname)
    if (len_trim(Tr%cmor_tendprefix) == 0) then
      Tr%id_trxh_tendency = register_diag_field('ocean_model', trim(shortnm)//'h_tendency', &
          diag%axesTL, Time, var_lname, conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          v_extensive=.true.)
      Tr%id_trxh_tendency_2d = register_diag_field('ocean_model', trim(shortnm)//'h_tendency_2d', &
          diag%axesT1, Time, "Vertical sum of "//trim(lowercase(var_lname)), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T)
    else
      cmor_var_lname = "Tendency of "//trim(cmor_longname)//" Expressed as "//&
                        trim(flux_longname)//" Content"
      Tr%id_trxh_tendency = register_diag_field('ocean_model', trim(shortnm)//'h_tendency', &
          diag%axesTL, Time, var_lname, conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          cmor_field_name=trim(Tr%cmor_tendprefix)//"tend", &
          cmor_standard_name=cmor_long_std(cmor_var_lname), cmor_long_name=cmor_var_lname, &
          v_extensive=.true.)
      cmor_var_lname = trim(cmor_var_lname)//" Vertical Sum"
      Tr%id_trxh_tendency_2d = register_diag_field('ocean_model', trim(shortnm)//'h_tendency_2d', &
          diag%axesT1, Time, "Vertical sum of "//trim(lowercase(var_lname)), &
          conv_units, conversion=Tr%conv_scale*US%s_to_T, &
          cmor_field_name=trim(Tr%cmor_tendprefix)//"tend_2d", &
          cmor_standard_name=cmor_long_std(cmor_var_lname), cmor_long_name=cmor_var_lname)
    endif
    if ((Tr%id_trxh_tendency > 0) .or. (Tr%id_trxh_tendency_2d > 0)) then
      call safe_alloc_ptr(Tr%Trxh_prev,isd,ied,jsd,jed,nz)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Tr%Trxh_prev(i,j,k) = Tr%t(i,j,k) * h(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! Vertical regridding/remapping tendencies
    if (use_ALE .and. Tr%remap_tr) then
      var_lname = "Vertical remapping tracer concentration tendency for "//trim(Reg%Tr(m)%name)
      Tr%id_remap_conc= register_diag_field('ocean_model', &
          trim(Tr%flux_nameroot)//'_tendency_vert_remap', diag%axesTL, Time, var_lname, &
          trim(units)//' s-1', conversion=Tr%conc_scale*US%s_to_T)

      var_lname = "Vertical remapping tracer content tendency for "//trim(Reg%Tr(m)%flux_longname)
      Tr%id_remap_cont = register_diag_field('ocean_model', &
          trim(Tr%flux_nameroot)//'h_tendency_vert_remap', &
          diag%axesTL, Time, var_lname, conv_units, v_extensive=.true., conversion=Tr%conv_scale*US%s_to_T)

      var_lname = "Vertical sum of vertical remapping tracer content tendency for "//&
                  trim(Reg%Tr(m)%flux_longname)
      Tr%id_remap_cont_2d = register_diag_field('ocean_model', &
          trim(Tr%flux_nameroot)//'h_tendency_vert_remap_2d', &
          diag%axesT1, Time, var_lname, conv_units, conversion=Tr%conv_scale*US%s_to_T)

    endif

    if (use_ALE .and. (Reg%ntr<MAX_FIELDS_) .and. Tr%remap_tr) then
      unit2 = trim(units)//"2"
      if (index(units(1:len_trim(units))," ") > 0) unit2 = "("//trim(units)//")2"
      Tr%id_tr_vardec = register_diag_field('ocean_model', trim(shortnm)//"_vardec", diag%axesTL, &
          Time, "ALE variance decay for "//lowercase(longname), &
          trim(unit2)//" s-1", conversion=Tr%conc_scale**2*US%s_to_T)
      if (Tr%id_tr_vardec > 0) then
        ! Set up a new tracer for this tracer squared
        m2 = Reg%ntr+1
        Tr%ind_tr_squared = m2
        call safe_alloc_ptr(Reg%Tr(m2)%t,isd,ied,jsd,jed,nz) ; Reg%Tr(m2)%t(:,:,:) = 0.0
        Reg%Tr(m2)%name = trim(shortnm)//"2"
        Reg%Tr(m2)%longname = "Squared "//trim(longname)
        Reg%Tr(m2)%units = unit2
        Reg%Tr(m2)%registry_diags = .false.
        Reg%Tr(m2)%ind_tr_squared = -1
        ! Augment the total number of tracers, including the squared tracers.
        Reg%ntr = Reg%ntr + 1
      endif
    endif

    ! KPP nonlocal term diagnostics
    if (use_KPP) then
      Tr%id_net_surfflux = register_diag_field('ocean_model', Tr%net_surfflux_name, diag%axesT1, Time, &
          Tr%net_surfflux_longname, trim(units)//' m s-1', conversion=Tr%conc_scale*GV%H_to_m*US%s_to_T)
      Tr%id_NLT_tendency = register_diag_field('ocean_model', "KPP_NLT_d"//trim(shortnm)//"dt", &
          diag%axesTL, Time, &
          trim(longname)//' tendency due to non-local transport of '//trim(lowercase(flux_longname))//&
          ', as calculated by [CVMix] KPP', trim(units)//' s-1', conversion=Tr%conc_scale*US%s_to_T)
      if (Tr%conv_scale == 0.001*GV%H_to_kg_m2) then
        conversion = GV%H_to_kg_m2
      else
        conversion = Tr%conv_scale
      endif
      ! We actually want conversion=Tr%conv_scale for all tracers, but introducing the local variable
      ! 'conversion' and setting it to GV%H_to_kg_m2 instead of 0.001*GV%H_to_kg_m2 for salt tracers
      ! keeps changes introduced by this refactoring limited to round-off level; as it turns out,
      ! there is a bug in the code and the NLT budget term for salinity is off by a factor of 10^3
      ! so introducing the 0.001 here will fix that bug.
      Tr%id_NLT_budget = register_diag_field('ocean_model', Tr%NLT_budget_name, &
          diag%axesTL, Time, &
          trim(flux_longname)//&
          ' content change due to non-local transport, as calculated by [CVMix] KPP', &
          conv_units, conversion=conversion*US%s_to_T, v_extensive=.true.)
    endif

  endif ; enddo

end procedure register_tracer_diagnostics
module procedure preALE_tracer_diagnostics
  integer :: i, j, k, is, ie, js, je, nz, m, m2
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do m=1,Reg%ntr ; if (Reg%Tr(m)%ind_tr_squared > 0) then
    m2 = Reg%Tr(m)%ind_tr_squared
  ! Update squared quantities
    do k=1,nz ; do j=js,je ; do i=is,ie
      Reg%Tr(m2)%T(i,j,k) = Reg%Tr(m)%T(i,j,k)**2
    enddo ; enddo ; enddo
  endif ; enddo

end procedure preALE_tracer_diagnostics
module procedure postALE_tracer_diagnostics
  real    :: work(SZI_(G),SZJ_(G),SZK_(GV)) ! Variance decay [CU2 T-1 ~> conc2 s-1]
  real    :: Idt ! The inverse of the time step [T-1 ~> s-1]
  integer :: i, j, k, is, ie, js, je, nz, m, m2
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! The "if" is to avoid NaNs if the diagnostic is called for a zero length interval
  Idt = 0.0 ; if (dt /= 0.0) Idt = 1.0 / dt

  do m=1,Reg%ntr ; if (Reg%Tr(m)%id_tr_vardec > 0) then
    m2 = Reg%Tr(m)%ind_tr_squared
    if (m2 < 1) call MOM_error(FATAL, "Bad value of Tr%ind_tr_squared for "//trim(Reg%Tr(m)%name))
  ! Update squared quantities
    do k=1,nz ; do j=js,je ; do i=is,ie
      work(i,j,k) = (Reg%Tr(m2)%T(i,j,k) - Reg%Tr(m)%T(i,j,k)**2) * Idt
    enddo ; enddo ; enddo
    call post_data(Reg%Tr(m)%id_tr_vardec, work, diag)
  endif ; enddo

end procedure postALE_tracer_diagnostics
module procedure post_tracer_diagnostics_at_sync
  real    :: work3d(SZI_(G),SZJ_(G),SZK_(GV)) ! The time tendency of a diagnostic [CU T-1 ~> conc s-1]
  real    :: work2d(SZI_(G),SZJ_(G)) ! The vertically integrated time tendency of a diagnostic
  real    :: Idt ! The inverse of the time step [T-1 ~> s-1]
  type(tracer_type), pointer :: Tr=>NULL()
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  Idt = 0. ; if (dt/=0.) Idt = 1.0 / dt ! The "if" is in case the diagnostic is called for a zero length interval

  ! Tendency diagnostics need to be posted on the grid from the last call to this routine
  call diag_save_grids(diag)
  call diag_copy_storage_to_diag(diag, diag_prev)
  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_tr > 0) call post_data(Tr%id_tr, Tr%t, diag)
    if (Tr%id_tendency > 0) then
      work3d(:,:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work3d(i,j,k)    = (Tr%t(i,j,k) - Tr%t_prev(i,j,k))*Idt
        tr%t_prev(i,j,k) =  Tr%t(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_tendency, work3d, diag, alt_h=diag_prev%h_state)
    endif
    if ((Tr%id_trxh_tendency > 0) .or. (Tr%id_trxh_tendency_2d > 0)) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work3d(i,j,k)     = (Tr%t(i,j,k)*h(i,j,k) - Tr%Trxh_prev(i,j,k)) * Idt
        Tr%Trxh_prev(i,j,k) =  Tr%t(i,j,k) * h(i,j,k)
      enddo ; enddo ; enddo
      if (Tr%id_trxh_tendency > 0) call post_data(Tr%id_trxh_tendency, work3d, diag, &
                                                  alt_h=diag_prev%h_state)
      if (Tr%id_trxh_tendency_2d > 0) then
        work2d(:,:) = 0.0
        do k=1,nz ; do j=js,je ; do i=is,ie
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo ; enddo ; enddo
        call post_data(Tr%id_trxh_tendency_2d, work2d, diag)
      endif
    endif
  endif ; enddo
  call diag_restore_grids(diag)

end procedure post_tracer_diagnostics_at_sync
module procedure post_tracer_transport_diagnostics
  integer :: i, j, k, is, ie, js, je, nz, m, khi
  real    :: work2d(SZI_(G),SZJ_(G))      ! The vertically integrated convergence of lateral advective
  real    :: frac_under_100m(SZI_(G),SZJ_(G),SZK_(GV)) ! weights used to compute 100m vertical integrals [nondim]
  real    :: ztop(SZI_(G),SZJ_(G)) ! position of the top interface [H ~> m or kg m-2]
  real    :: zbot(SZI_(G),SZJ_(G)) ! position of the bottom interface [H ~> m or kg m-2]
  type(tracer_type), pointer :: Tr=>NULL()
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! If any tracers are posting 100m vertical integrals, compute weights
  frac_under_100m(:,:,:) = 0.0
  ! khi will be the largest layer index corresponding where ztop < 100m and ztop >= 100m
  ! in any column (we can reduce computation of 100m integrals by only looping through khi
  ! rather than GV%ke)
  khi = 0
  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_zint_100m > 0) then
      zbot(:,:) = 0.0
      do k=1, nz
        do j=js,je ; do i=is,ie
          ztop(i,j) = zbot(i,j)
          zbot(i,j) = ztop(i,j) + h_diag(i,j,k)*GV%H_to_m
          if (zbot(i,j) <= 100.0) then
            frac_under_100m(i,j,k) = 1.0
          elseif (ztop(i,j) < 100.0) then
            frac_under_100m(i,j,k) = (100.0 - ztop(i,j)) / (zbot(i,j) - ztop(i,j))
          else
            frac_under_100m(i,j,k) = 0.0
          endif
          ! frac_under_100m(i,j,k) = max(0, min(1.0, (100.0 - ztop(i,j)) / (zbot(i,j) - ztop(i,j))))
        enddo ; enddo
        if (any(frac_under_100m(:,:,k) > 0)) khi = k
      enddo
      exit
    endif
  endif ; enddo

  do m=1,Reg%ntr ; if (Reg%Tr(m)%registry_diags) then
    Tr => Reg%Tr(m)
    if (Tr%id_tr_post_horzn> 0) call post_data(Tr%id_tr_post_horzn, Tr%t, diag)
    if (Tr%id_adx > 0) call post_data(Tr%id_adx, Tr%ad_x, diag, alt_h=h_diag)
    if (Tr%id_ady > 0) call post_data(Tr%id_ady, Tr%ad_y, diag, alt_h=h_diag)
    if (Tr%id_dfx > 0) call post_data(Tr%id_dfx, Tr%df_x, diag, alt_h=h_diag)
    if (Tr%id_dfy > 0) call post_data(Tr%id_dfy, Tr%df_y, diag, alt_h=h_diag)
    if (Tr%id_adx_2d > 0) call post_data(Tr%id_adx_2d, Tr%ad2d_x, diag)
    if (Tr%id_ady_2d > 0) call post_data(Tr%id_ady_2d, Tr%ad2d_y, diag)
    if (Tr%id_dfx_2d > 0) call post_data(Tr%id_dfx_2d, Tr%df2d_x, diag)
    if (Tr%id_dfy_2d > 0) call post_data(Tr%id_dfy_2d, Tr%df2d_y, diag)
    if (Tr%id_adv_xy > 0) call post_data(Tr%id_adv_xy, Tr%advection_xy, diag, alt_h=h_diag)
    if (Tr%id_adv_xy_2d > 0) then
      work2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work2d(i,j) = work2d(i,j) + Tr%advection_xy(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_adv_xy_2d, work2d, diag)
    endif

    ! A few diagnostics introduce with MARBL driver
    ! Compute full-depth vertical integral
    if (Tr%id_zint > 0) then
      work2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work2d(i,j) = work2d(i,j) + (h_diag(i,j,k)*GV%H_to_m)*tr%t(i,j,k)
      enddo ; enddo ; enddo
      call post_data(Tr%id_zint, work2d, diag)
    endif

    ! Compute 100m vertical integral
    if (Tr%id_zint_100m > 0) then
      work2d(:,:) = 0.0
      do k=1,khi ; do j=js,je ; do i=is,ie
        work2d(i,j) = work2d(i,j) + frac_under_100m(i,j,k)*((h_diag(i,j,k)*GV%H_to_m)*tr%t(i,j,k))
      enddo ; enddo ; enddo
      call post_data(Tr%id_zint_100m, work2d, diag)
    endif

    ! Surface values of tracers
    if (Tr%id_SURF > 0) call post_data(Tr%id_SURF, Tr%t(:,:,1), diag)
  endif ; enddo

end procedure post_tracer_transport_diagnostics
module procedure tracer_array_chksum
  integer :: m
  do m=1,ntr
    call hchksum(Tr(m)%t, mesg//trim(Tr(m)%name), G%HI, unscale=Tr(m)%conc_scale)
  enddo

end procedure tracer_array_chksum
module procedure tracer_Reg_chksum
  integer :: m
  if (.not.associated(Reg)) return

  do m=1,Reg%ntr
    call hchksum(Reg%Tr(m)%t, mesg//trim(Reg%Tr(m)%name), G%HI, unscale=Reg%Tr(m)%conc_scale)
  enddo

end procedure tracer_Reg_chksum
module procedure tracer_array_chkinv
  real :: vol_scale ! The dimensional scaling factor to convert volumes to m3 [m3 H-1 L-2 ~> 1] or cell
  real :: tr_inv(SZI_(G),SZJ_(G),SZK_(GV)) ! Volumetric or mass-based tracer inventory in
  real :: total_inv ! The total amount of tracer [conc m3] or [conc kg]
  integer :: is, ie, js, je, nz
  integer :: i, j, k, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  vol_scale = GV%H_to_MKS*G%US%L_to_m**2
  do m=1,ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_inv(i,j,k) = Tr(m)%conc_scale*Tr(m)%t(i,j,k) * &
                      (vol_scale * h(i,j,k) * G%areaT(i,j)*G%mask2dT(i,j))
    enddo ; enddo ; enddo
    total_inv = reproducing_sum(tr_inv, is+(1-G%isd), ie+(1-G%isd), js+(1-G%jsd), je+(1-G%jsd))
    if (is_root_pe()) write(0,'(A,1X,A5,1X,ES25.16,1X,A)') &
                      "h-point: inventory", Tr(m)%name, total_inv, mesg
  enddo

end procedure tracer_array_chkinv
module procedure tracer_Reg_chkinv
  real :: vol_scale ! The dimensional scaling factor to convert volumes to m3 [m3 H-1 L-2 ~> 1] or cell
  real :: tr_inv(SZI_(G),SZJ_(G),SZK_(GV)) ! Volumetric or mass-based tracer inventory in
  real :: total_inv ! The total amount of tracer [conc m3] or [conc kg]
  integer :: is, ie, js, je, nz
  integer :: i, j, k, m
  if (.not.associated(Reg)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  vol_scale = GV%H_to_MKS*G%US%L_to_m**2
  do m=1,Reg%ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_inv(i,j,k) = Reg%Tr(m)%conc_scale*Reg%Tr(m)%t(i,j,k) * &
                      (vol_scale * h(i,j,k) * G%areaT(i,j)*G%mask2dT(i,j))
    enddo ; enddo ; enddo
    total_inv = reproducing_sum(tr_inv, is+(1-G%isd), ie+(1-G%isd), js+(1-G%jsd), je+(1-G%jsd))
    if (is_root_pe()) write(0,'(A,1X,A5,1X,ES25.16,1X,A)') &
                      "h-point: inventory", Reg%Tr(m)%name, total_inv, mesg
  enddo

end procedure tracer_Reg_chkinv
module procedure tracer_name_lookup
  do n=1,Reg%ntr
    if (lowercase(Reg%Tr(n)%name) == lowercase(name)) then
      tr_ptr => Reg%Tr(n)
      return
    endif
  enddo

  call MOM_error(FATAL,"MOM cannot find registered tracer: "//name)

end procedure tracer_name_lookup
module procedure tracer_registry_init
  integer, save :: init_calls = 0
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.
  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", all_default=.true.)

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("tracer_registry_init called ",I0, &
      &" times with different registry pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_tracer "//mesg)
  endif

end procedure tracer_registry_init
module procedure tracer_registry_end
  if (associated(Reg)) deallocate(Reg)
end procedure tracer_registry_end
end submodule MOM_tracer_registry_s
