submodule (MOM_EOS_base_type) MOM_EOS_base_type_s
  implicit none
contains
  module procedure a_density_fn
    if (present(rho_ref)) then
      a_density_fn = this%density_anomaly_elem(T, S, pressure, rho_ref)
    else
      a_density_fn = this%density_elem(T, S, pressure)
    endif

  end procedure a_density_fn
  module procedure a_calculate_density_scalar
    if (present(rho_ref)) then
      rho = this%density_anomaly_elem(T, S, pressure, rho_ref)
    else
      rho = this%density_elem(T, S, pressure)
    endif

  end procedure a_calculate_density_scalar
  module procedure a_calculate_density_array
    integer :: js, je
    js = start
    je = start+npts-1

    if (present(rho_ref)) then
      rho(js:je) = this%density_anomaly_elem(T(js:je), S(js:je), pressure(js:je), rho_ref)
    else
      rho(js:je) = this%density_elem(T(js:je), S(js:je), pressure(js:je))
    endif

  end procedure a_calculate_density_array
  module procedure a_spec_vol_fn
    if (present(spv_ref)) then
      a_spec_vol_fn = this%spec_vol_anomaly_elem(T, S, pressure, spv_ref)
    else
      a_spec_vol_fn = this%spec_vol_elem(T, S, pressure)
    endif

  end procedure a_spec_vol_fn
  module procedure a_calculate_spec_vol_scalar
    if (present(spv_ref)) then
      specvol = this%spec_vol_anomaly_elem(T, S, pressure, spv_ref)
    else
      specvol = this%spec_vol_elem(T, S, pressure)
    endif

  end procedure a_calculate_spec_vol_scalar
  module procedure a_calculate_spec_vol_array
    integer :: js, je
    js = start
    je = start+npts-1

    if (present(spv_ref)) then
      specvol(js:je) = this%spec_vol_anomaly_elem(T(js:je), S(js:je), pressure(js:je), spv_ref)
    else
      specvol(js:je) = this%spec_vol_elem(T(js:je), S(js:je), pressure(js:je) )
    endif

  end procedure a_calculate_spec_vol_array
  module procedure a_calculate_density_derivs_scalar
    call this%calculate_density_derivs_elem(T, S, P, drho_dt, drho_ds)

  end procedure a_calculate_density_derivs_scalar
  module procedure a_calculate_density_derivs_array
    integer :: js, je
    js = start
    je = start+npts-1

    call this%calculate_density_derivs_elem(T(js:je), S(js:je), pressure(js:je), drho_dt(js:je), drho_ds(js:je))

  end procedure a_calculate_density_derivs_array
  module procedure a_calculate_density_second_derivs_scalar
    call this%calculate_density_second_derivs_elem(T, S, pressure, &
                      drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)

  end procedure a_calculate_density_second_derivs_scalar
  module procedure a_calculate_density_second_derivs_array
    integer :: js, je
    js = start
    je = start+npts-1

    call this%calculate_density_second_derivs_elem(T(js:je), S(js:je), pressure(js:je), &
                              drho_ds_ds(js:je), drho_ds_dt(js:je), drho_dt_dt(js:je), &
                              drho_ds_dp(js:je), drho_dt_dp(js:je))

  end procedure a_calculate_density_second_derivs_array
  module procedure a_calculate_specvol_derivs_array
    integer :: js, je
    js = start
    je = start+npts-1

    call this%calculate_specvol_derivs_elem(T(js:je), S(js:je), pressure(js:je), &
                                            dSV_dT(js:je), dSV_dS(js:je))

  end procedure a_calculate_specvol_derivs_array
  module procedure a_calculate_compress_array
    integer :: js, je
    js = start
    je = start+npts-1

    call this%calculate_compress_elem(T(js:je), S(js:je), pressure(js:je), &
                                      rho(js:je), drho_dp(js:je))

  end procedure a_calculate_compress_array
end submodule MOM_EOS_base_type_s
