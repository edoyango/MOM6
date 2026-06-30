submodule (MOM_temperature_convert) MOM_temperature_convert_s
  implicit none
contains
module procedure poTemp_to_consTemp
  real :: x2 ! Absolute salinity normalized by a plausible salinity range [nondim]
  real :: x  ! Square root of normalized absolute salinity [nondim]
  x2 = max(I_S0 * Sa, 0.0)
  x = sqrt(x2)

  Tc = H00 + (T*(H01 + T*(H02 +  T*(H03 +  T*(H04  + T*(H05 + T*(H06 + T* H07)))))) &
                    + x2*(H20 + (T*(H21 +  T*(H22  + T*(H23 + T*(H24 + T*(H25 + T*H26))))) &
                              +  x*(H30 + (T*(H31  + T*(H32 + T*(H33 + T* H34))) &
                                         + x*(H40 + (T*(H41 + T*(H42 + T*(H43 + T*(H44 + T*H45)))) &
                                                  +  x*(H50 + x*(H60 + x* H70)) )) )) )) )

end procedure poTemp_to_consTemp
module procedure dTc_dTp
  real :: x2 ! Absolute salinity normalized by a plausible salinity range [nondim]
  real :: x  ! Square root of normalized absolute salinity [nondim]
  x2 = max(I_S0 * Sa, 0.0)
  x = sqrt(x2)

  dTc_dTp = (     H01 + T*(2.*H02 + T*(3.*H03 + T*(4.*H04 + T*(5.*H05 + T*(6.*H06 + T*(7.*H07)))))) ) &
      + x2*(     (H21 + T*(2.*H22 + T*(3.*H23 + T*(4.*H24 + T*(5.*H25 + T*(6.*H26)))))) &
         +  x*(  (H31 + T*(2.*H32 + T*(3.*H33 + T*(4.*H34)))) &
             + x*(H41 + T*(2.*H42 + T*(3.*H43 + T*(4.*H44 + T*(5.*H45))))) ) )

end procedure dTc_dTp
module procedure consTemp_to_poTemp
  real :: Tp_num    ! The numerator  of a simple expression for potential temperature [degC]
  real :: I_Tp_den  ! The inverse of the denominator of a simple expression for potential temperature [nondim]
  real :: Tc_diff   ! The difference between an estimate of conservative temperature and its target [degC]
  real :: Tp_old    ! A previous estimate of the potential tempearture [degC]
  real :: dTp_dTc   ! The partial derivative of potential temperature with conservative temperature [nondim]
  real, parameter :: TPN00 = -1.446013646344788e-2               ! Simple fit numerator constant [degC]
  real, parameter :: TPN10 = -3.305308995852924e-3*Sprac_Sref    ! Simple fit numerator Sa coef. [degC ppt-1]
  real, parameter :: TPN20 =  1.062415929128982e-4*Sprac_Sref**2 ! Simple fit numerator Sa**2 coef. [degC ppt-2]
  real, parameter :: TPN01 =  9.477566673794488e-1               ! Simple fit numerator Tc coef. [nondim]
  real, parameter :: TPN11 =  2.166591947736613e-3*Sprac_Sref    ! Simple fit numerator Sa * Tc coef. [ppt-1]
  real, parameter :: TPN02 =  3.828842955039902e-3               ! Simple fit numerator Tc**2 coef. [degC-1]
  real, parameter :: TPD10 =  6.506097115635800e-4*Sprac_Sref    ! Simple fit denominator Sa coef. [ppt-1]
  real, parameter :: TPD01 =  3.830289486850898e-3               ! Simple fit denominator Tc coef. [degC-1]
  real, parameter :: TPD02 =  1.247811760368034e-6               ! Simple fit denominator Tc**2 coef. [degC-2]
  Tp_num = TPN00 + (Sa*(TPN10 + TPN20*Sa) + Tc*(TPN01 + (TPN11*Sa + TPN02*Tc)))
  I_Tp_den = 1.0 / (1.0 + (TPD10*Sa + Tc*(TPD01 + TPD02*Tc)))
  Tp = Tp_num*I_Tp_den
  dTp_dTc = ((TPN01 + (TPN11*Sa + 2.*TPN02*Tc)) - (TPD01 + 2.*TPD02*Tc)*Tp)*I_Tp_den

  ! Start the 1.5 iterations through the modified Newton-Raphson iterative method, which is also known
  ! as the Newton-McDougall method.  In this case 1.5 iterations converge to 64-bit machine precision
  ! for oceanographically relevant temperatures and salinities.

  Tc_diff = poTemp_to_consTemp(Tp, Sa) - Tc
  Tp_old = Tp
  Tp = Tp_old - Tc_diff*dTp_dTc

  dTp_dTc = 1.0 / dTc_dTp(0.5*(Tp + Tp_old), Sa)

  Tp = Tp_old - Tc_diff*dTp_dTc
  Tc_diff = poTemp_to_consTemp(Tp, Sa) - Tc
  Tp_old = Tp

  Tp = Tp_old - Tc_diff*dTp_dTc

end procedure consTemp_to_poTemp
end submodule MOM_temperature_convert_s
