submodule (MOM_CoriolisAdv) MOM_CoriolisAdv_s
#include <MOM_memory.h>
  implicit none
contains
module procedure CorAdCalc
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    q, &        ! Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    qS, &       ! Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    Ih_q, &     ! The inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    h_q, &      ! The thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    Area_q      ! The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    a, b, c, d  ! a, b, c, & d are combinations of the potential vorticities
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Area_h, &   ! The ocean area at h points [L2 ~> m2].  Area_h is used to find the
                ! average thickness in the denominator of q.  0 for land points.
    KE          ! Kinetic energy per unit mass [L2 T-2 ~> m2 s-2], KE = (u^2 + v^2)/2.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    hArea_u, &  ! The cell area weighted thickness interpolated to u points
                ! times the effective areas [H L2 ~> m3 or kg].
    KEx, &      ! The zonal gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEx = d/dx KE.
    uh_center   ! Transport based on arithmetic mean h at u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    hArea_v, &  ! The cell area weighted thickness interpolated to v points
                ! times the effective areas [H L2 ~> m3 or kg].
    KEy, &      ! The meridional gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEy = d/dy KE.
    vh_center   ! Transport based on arithmetic mean h at v-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    uh_min, uh_max, &   ! The smallest and largest estimates of the zonal volume fluxes through
                        ! the faces (i.e. u*h*dy) [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh_min, vh_max, &   ! The smallest and largest estimates of the meridional volume fluxes through
                        ! the faces (i.e. v*h*dx) [H L2 T-1 ~> m3 s-1 or kg s-1]
    ep_u, ep_v  ! Additional pseudo-Coriolis terms in the Arakawa and Lamb
                ! discretization [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! Contributions to the circulation around q-points [L2 T-1 ~> m2 s-1]
    dvSdx, duSdy, & ! idem. for Stokes drift [L2 T-1 ~> m2 s-1]
    rel_vort, & ! Relative vorticity at q-points [T-1 ~> s-1].
    abs_vort, & ! Absolute vorticity at q-points [T-1 ~> s-1].
    stk_vort, & ! Stokes vorticity at q-points [T-1 ~> s-1].
    q2          ! Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    PV, &       ! A diagnostic array of the potential vorticities [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    RV          ! A diagnostic array of the relative vorticities [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: CAuS ! Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: CAvS ! Stokes contribution to CAv [L T-2 ~> m s-2]
  real :: fv1, fv2, fv3, fv4   ! (f+rv)*v at the 4 points surrounding a u points[L T-2 ~> m s-2]
  real :: fu1, fu2, fu3, fu4   ! -(f+rv)*u at the 4 points surrounding a v point [L T-2 ~> m s-2]
  real :: max_fv, max_fu       ! The maximum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real :: min_fv, min_fu       ! The minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]

  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  real, parameter :: C1_24 = 1.0 / 24.0 ! C1_24 = 1/24 [nondim]
  real :: max_Ihq, min_Ihq       ! The maximum and minimum of the nearby Ihq [H-1 ~> m-1 or m2 kg-1].
  real :: hArea_q                ! The sum of area times thickness of the cells
                                 ! surrounding a q point [H L2 ~> m3 or kg].
  real :: vol_neglect            ! A volume so small that is expected to be
                                 ! lost in roundoff [H L2 ~> m3 or kg].
  real :: area_neglect           ! An area so small that is expected to be
                                 ! lost in roundoff [L2 ~> m2].
  real :: temp1, temp2           ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: eps_vel                ! A tiny, positive velocity [L T-1 ~> m s-1].

  real :: uhc, vhc               ! Centered estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uhm, vhm               ! The input estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: c1, c2, c3, slope      ! Nondimensional parameters for the Coriolis limiter scheme [nondim]

  real :: Fe_m2         ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_lin       ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_m1        ! The ratio of the maximum neighboring inverse thickness
                        ! to the minimum inverse thickness minus 1 [nondim]. rat_m1 >= 0.
  real :: AL_wt         ! The relative weight of the Arakawa & Lamb scheme to the
                        ! Arakawa & Hsu scheme [nondim], between 0 and 1.
  real :: Sad_wt        ! The relative weight of the Sadourny energy scheme to
                        ! the other two with the ARAKAWA_LAMB_BLEND scheme [nondim],
                        ! between 0 and 1.

  real :: Heff1, Heff2  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: Heff3, Heff4  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: h_tiny        ! A very small thickness [H ~> m or kg m-2].
  real :: UHeff, VHeff  ! More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: QUHeff,QVHeff ! More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
  integer :: i, j, k, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: Is_q, Ie_q, Js_q, Je_q  ! The scheme-dependent range of values at which vorticity is set.
  logical :: Stokes_VF
  real :: u_v, v_u      ! u_v is the u velocity at v point, v_u is the v velocity at u point [L T-1 ~> m s-1]
  real :: q_v, q_u      ! PV at the u and v points [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  integer :: seventh_order, fifth_order, third_order ! Order of accuracy for the WENO calculations
  real :: u_q8(8) ! Eight-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q6(6) ! Six-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q4(4) ! Four-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q8(8) ! Eight-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q6(6) ! Six-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q4(4) ! Four-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  integer :: stencil    ! Stencil size of WENO scheme

! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
!   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_CoriolisAdv: Module must be initialized before it is used.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke
  vol_neglect = GV%H_subroundoff * (1e-4 * US%m_to_L)**2
  area_neglect = (1e-4 * US%m_to_L)**2
  eps_vel = 1.0e-10*US%m_s_to_L_T
  h_tiny = GV%Angstrom_H  ! Perhaps this should be set to h_neglect instead.

  stencil = CoriolisAdv_stencil(CS)

  if ((CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) .or. (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) .or. &
      (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO)) then
    Is_q = is - stencil ; Ie_q = ie + stencil - 1 ; Js_q = js - stencil ; Je_q = je + stencil - 1
  else
    Is_q = G%IscB - 1 ; Ie_q = G%IecB + 1 ; Js_q = G%JscB - 1 ; Je_q = G%JecB + 1
  endif

  !$OMP parallel do default(private) shared(Is_q,Ie_q,Js_q,Je_q,G,Area_h)
  do j=Js_q,Je_q+1 ; do I=Is_q,Ie_q+1
    Area_h(i,j) = G%mask2dT(i,j) * G%areaT(i,j)
  enddo ; enddo
  if (associated(OBC)) then ; do n=1,OBC%number_of_segments
    if (.not. OBC%segment(n)%on_pe) cycle
    I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
    if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
      do i = max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          Area_h(i,j+1) = Area_h(i,j)
        else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
          Area_h(i,j) = Area_h(i,j+1)
        endif
      enddo
    elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
      do j = max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          Area_h(i+1,j) = Area_h(i,j)
        else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
          Area_h(i,j) = Area_h(i+1,j)
        endif
      enddo
    endif
  enddo ; endif
  !$OMP parallel do default(private) shared(Is_q,Ie_q,Js_q,Je_q,G,Area_h,Area_q)
  do J=Js_q,Je_q ; do I=Is_q,Ie_q
    Area_q(i,j) = (Area_h(i,j) + Area_h(i+1,j+1)) + &
                  (Area_h(i+1,j) + Area_h(i,j+1))
  enddo ; enddo

  Stokes_VF = .false.
  if (present(Waves)) then ; if (associated(Waves)) then
    Stokes_VF = Waves%Stokes_VF
  endif ; endif

  !$OMP parallel do default(private) shared(u,v,h,uh,vh,CAu,CAv,G,GV,CS,AD,Area_h,Area_q,&
  !$OMP                        RV,PV,is,ie,js,je,Isq,Ieq,Jsq,Jeq,Is_q,Ie_q,Js_q,Je_q,nz,vol_neglect,&
  !$OMP                        h_tiny,OBC,eps_vel,area_neglect,pbv,Stokes_VF,stencil)
  do k=1,nz

    ! Here the second order accurate layer potential vorticities, q,
    ! are calculated.  hq is  second order accurate in space.  Relative
    ! vorticity is second order accurate everywhere with free slip b.c.s,
    ! but only first order accurate at boundaries with no slip b.c.s.
    ! First calculate the contributions to the circulation around the q-point.
    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do J=Js_q,Je_q ; do I=Is_q,Ie_q
          dvSdx(I,J) = (-Waves%us_y(i+1,J,k)*G%dyCv(i+1,J)) - &
                       (-Waves%us_y(i,J,k)*G%dyCv(i,J))
          duSdy(I,J) = (-Waves%us_x(I,j+1,k)*G%dxCu(I,j+1)) - &
                       (-Waves%us_x(I,j,k)*G%dxCu(I,j))
        enddo ; enddo
      endif
      if (.not. Waves%Passive_Stokes_VF) then
        do J=Js_q,Je_q ; do I=Is_q,Ie_q
          dvdx(I,J) = ((v(i+1,J,k)-Waves%us_y(i+1,J,k))*G%dyCv(i+1,J)) - &
                      ((v(i,J,k)-Waves%us_y(i,J,k))*G%dyCv(i,J))
          dudy(I,J) = ((u(I,j+1,k)-Waves%us_x(I,j+1,k))*G%dxCu(I,j+1)) - &
                      ((u(I,j,k)-Waves%us_x(I,j,k))*G%dxCu(I,j))
        enddo ; enddo
      else
        do J=Js_q,Je_q ; do I=Is_q,Ie_q
          dvdx(I,J) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
          dudy(I,J) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
        enddo ; enddo
      endif
    else
      do J=Js_q,Je_q ; do I=Is_q,Ie_q
        dvdx(I,J) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
        dudy(I,J) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
      enddo ; enddo
    endif
    do J=Js_q,Je_q ; do i=Is_q,Ie_q+1
      hArea_v(i,J) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i,j+1) * h(i,j+1,k)))
    enddo ; enddo
    do j=Js_q,Je_q+1 ; do I=Is_q,Ie_q
      hArea_u(I,j) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i+1,j) * h(i+1,j,k)))
    enddo ; enddo

    if (CS%Coriolis_En_Dis) then
      do j=Jsq,Jeq+1 ; do I=is-1,ie
        uh_center(I,j) = 0.5 * ((G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k)) * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo
      do J=js-1,je ; do i=Isq,Ieq+1
        vh_center(i,J) = 0.5 * ((G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k)) * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo
    endif

    ! Adjust circulation components to relative vorticity and thickness projected onto
    ! velocity points on open boundaries.
    if (associated(OBC)) then ; do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
        select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dvdx(I,J) = 0. ; dudy(I,J) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dudy(I,J) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*G%dxCu(I,j)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J) = 2.0*(u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dxCu(I,j+1)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j)*G%dyBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j+1)*G%dyBu(I,J)
              endif
            enddo
        end select

        ! Project thicknesses across OBC points with a no-gradient condition.
        do i = max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j,k)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
          endif
        enddo

        if (CS%Coriolis_En_Dis) then
          do i = max(Isq,OBC%segment(n)%HI%isd), min(Ieq+1,OBC%segment(n)%HI%ied)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              vh_center(i,J) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              vh_center(i,J) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j+1,k)
            endif
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
        select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J) = 0. ; dudy(I,J) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*G%dyCv(i,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J) = 2.0*(v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dyCv(i+1,J)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i,J)*G%dxBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i+1,J)*G%dxBu(I,J)
              endif
            enddo
        end select

        ! Project thicknesses across OBC points with a no-gradient condition.
        do j = max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            hArea_u(I,j) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i,j,k)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            hArea_u(I,j) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i+1,j,k)
          endif
        enddo
        if (CS%Coriolis_En_Dis) then
          do j = max(Jsq,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              uh_center(I,j) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              uh_center(I,j) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i+1,j,k)
            endif
          enddo
        endif
      endif
    enddo ; endif

    if (associated(OBC)) then ; do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      ! Now project thicknesses across cell-corner points in the OBCs.  The two
      ! projections have to occur in sequence and can not be combined easily.
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
        do I = max(Is_q,OBC%segment(n)%HI%IsdB), min(Ie_q,OBC%segment(n)%HI%IedB)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            if (Area_h(i,j) + Area_h(i+1,j) > 0.0) then
              hArea_u(I,j+1) = hArea_u(I,j) * ((Area_h(i,j+1) + Area_h(i+1,j+1)) / &
                                               (Area_h(i,j) + Area_h(i+1,j)))
            else ; hArea_u(I,j+1) = 0.0 ; endif
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            if (Area_h(i,j+1) + Area_h(i+1,j+1) > 0.0) then
              hArea_u(I,j) = hArea_u(I,j+1) * ((Area_h(i,j) + Area_h(i+1,j)) / &
                                               (Area_h(i,j+1) + Area_h(i+1,j+1)))
            else ; hArea_u(I,j) = 0.0 ; endif
          endif
        enddo
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
        do J = max(Js_q,OBC%segment(n)%HI%JsdB), min(Je_q,OBC%segment(n)%HI%JedB)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            if (Area_h(i,j) + Area_h(i,j+1) > 0.0) then
              hArea_v(i+1,J) = hArea_v(i,J) * ((Area_h(i+1,j) + Area_h(i+1,j+1)) / &
                                               (Area_h(i,j) + Area_h(i,j+1)))
            else ; hArea_v(i+1,J) = 0.0 ; endif
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
            if (Area_h(i+1,j) + Area_h(i+1,j+1) > 0.0) then
              hArea_v(i,J) = hArea_v(i+1,J) * ((Area_h(i,j) + Area_h(i,j+1)) / &
                                               (Area_h(i+1,j) + Area_h(i+1,j+1)))
            else ; hArea_v(i,J) = 0.0 ; endif
          endif
        enddo
      endif
    enddo ; endif

    if (CS%no_slip) then
      do J=Js_q,Je_q ; do I=Is_q,Ie_q
        rel_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvdx(I,J) - dudy(I,J)) * G%IareaBu(I,J)
      enddo ; enddo
      if (Stokes_VF) then
        if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
          do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
            stk_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J) - duSdy(I,J)) * G%IareaBu(I,J)
          enddo ; enddo
        endif
      endif
    else
      do J=Js_q,Je_q ; do I=Is_q,Ie_q
        rel_vort(I,J) = G%mask2dBu(I,J) * (dvdx(I,J) - dudy(I,J)) * G%IareaBu(I,J)
      enddo ; enddo
      if (Stokes_VF) then
        if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
          do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
            stk_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J) - duSdy(I,J)) * G%IareaBu(I,J)
          enddo ; enddo
        endif
      endif
    endif

    do J=Js_q,Je_q ; do I=Is_q,Ie_q
      abs_vort(I,J) = G%CoriolisBu(I,J) + rel_vort(I,J)
    enddo ; enddo

    do J=Js_q,Je_q ; do I=Is_q,Ie_q
      hArea_q = (hArea_u(I,j) + hArea_u(I,j+1)) + (hArea_v(i,J) + hArea_v(i+1,J))
      Ih_q(I,J) = Area_q(I,J) / (hArea_q + vol_neglect)
      h_q(I,J) = hArea_q / max(Area_q(I,J), area_neglect)
      q(I,J) = abs_vort(I,J) * Ih_q(I,J)
    enddo ; enddo

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          qS(I,J) = stk_vort(I,J) * Ih_q(I,J)
        enddo ; enddo
      endif
    endif

    if (CS%id_rv > 0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        RV(I,J,k) = rel_vort(I,J)
      enddo ; enddo
    endif

    if (CS%id_PV > 0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        PV(I,J,k) = q(I,J)
      enddo ; enddo
    endif

    if (associated(AD%rv_x_v) .or. associated(AD%rv_x_u)) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        q2(I,J) = rel_vort(I,J) * Ih_q(I,J)
      enddo ; enddo
    endif

    !   a, b, c, and d are combinations of neighboring potential
    ! vorticities which form the Arakawa and Hsu vorticity advection
    ! scheme.  All are defined at u grid points.

    if (CS%Coriolis_Scheme == ARAKAWA_HSU90) then
      do j=Jsq,Jeq+1
        do I=is-1,Ieq
          a(I,j) = (q(I,J) + (q(I+1,J) + q(I,J-1))) * C1_12
          d(I,j) = ((q(I,J) + q(I+1,J-1)) + q(I,J-1)) * C1_12
        enddo
        do I=Isq,Ieq
          b(I,j) = (q(I,J) + (q(I-1,J) + q(I,J-1))) * C1_12
          c(I,j) = ((q(I,J) + q(I-1,J-1)) + q(I,J-1)) * C1_12
        enddo
      enddo
    elseif (CS%Coriolis_Scheme == ARAKAWA_LAMB81) then
      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        a(I-1,j) = (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        d(I-1,j) = ((q(I,j) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        b(I,j) =   ((q(I,J) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        c(I,j) =   (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        ep_u(i,j) = ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == AL_BLEND) then
      Fe_m2 = CS%F_eff_max_blend - 2.0
      rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend, 1.0e-16)

      ! This allows the code to always give Sadourny Energy
      if (CS%F_eff_max_blend <= 2.0) then ; Fe_m2 = -1. ; rat_lin = -1.0 ; endif

      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        min_Ihq = MIN(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        max_Ihq = MAX(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        rat_m1 = 1.0e15
        if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
        ! The weights used here are designed to keep the effective Coriolis
        ! acceleration from any one point on its neighbors within a factor
        ! of F_eff_max.  The minimum permitted value is 2 (the factor for
        ! Sadourny's energy conserving scheme).

        ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
        if (rat_m1 <= Fe_m2) then ; AL_wt = 1.0
        elseif (rat_m1 < 1.5*Fe_m2) then ; AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0
        else ; AL_wt = 0.0 ; endif

        ! Determine the relative weights of Sadourny Energy vs. the other two.
        if (rat_m1 <= 1.5*Fe_m2) then ; Sad_wt = 0.0
        elseif (rat_m1 <= rat_lin) then
          Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1
        elseif (rat_m1 < 2.0*rat_lin) then
          Sad_wt = 1.0 - (CS%wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin)
        else ; Sad_wt = 1.0 ; endif

        a(I-1,j) = Sad_wt * 0.25 * q(I-1,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J) + AL_wt*q(I,J-1)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        d(I-1,j) = Sad_wt * 0.25 * q(I-1,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J-1) + AL_wt*q(I,J)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        b(I,j) =   Sad_wt * 0.25 * q(I,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J) + AL_wt*q(I-1,J-1)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        c(I,j) =   Sad_wt * 0.25 * q(I,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J-1) + AL_wt*q(I-1,J)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        ep_u(i,j) = AL_wt  * ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = AL_wt * (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    endif

    if (CS%Coriolis_En_Dis) then
    !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
      c1 = 1.0-1.5*0.5 ; c2 = 1.0-0.5 ; c3 = 2.0 ; slope = 0.5

      do j=Jsq,Jeq+1 ; do I=is-1,ie
        uhc = uh_center(I,j)
        uhm = uh(I,j,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dy_Cu(I,j) == 0.0) uhc = uhm

        if (abs(uhc) < 0.1*abs(uhm)) then
          uhm = 10.0*uhc
        elseif (abs(uhc) > c1*abs(uhm)) then
          if (abs(uhc) < c2*abs(uhm)) then ; uhc = (3.0*uhc+(1.0-c2*3.0)*uhm)
          elseif (abs(uhc) <= c3*abs(uhm)) then ; uhc = uhm
          else ; uhc = slope*uhc+(1.0-c3*slope)*uhm
          endif
        endif

        if (uhc > uhm) then
          uh_min(I,j) = uhm ; uh_max(I,j) = uhc
        else
          uh_max(I,j) = uhm ; uh_min(I,j) = uhc
        endif
      enddo ; enddo
      do J=js-1,je ; do i=Isq,Ieq+1
        vhc = vh_center(i,J)
        vhm = vh(i,J,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dx_Cv(i,J) == 0.0) vhc = vhm

        if (abs(vhc) < 0.1*abs(vhm)) then
          vhm = 10.0*vhc
        elseif (abs(vhc) > c1*abs(vhm)) then
          if (abs(vhc) < c2*abs(vhm)) then ; vhc = (3.0*vhc+(1.0-c2*3.0)*vhm)
          elseif (abs(vhc) <= c3*abs(vhm)) then ; vhc = vhm
          else ; vhc = slope*vhc+(1.0-c3*slope)*vhm
          endif
        endif

        if (vhc > vhm) then
          vh_min(i,J) = vhm ; vh_max(i,J) = vhc
        else
          vh_max(i,J) = vhm ; vh_min(i,J) = vhc
        endif
      enddo ; enddo
    endif

    ! Calculate KE and the gradient of KE
    call gradKE(u(:,:,k), v(:,:,k), h(:,:,k), KE, KEx, KEy, G, GV, US, CS)

    ! Calculate the tendencies of zonal velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAu =  q * vh - d(KE)/dx.
    if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (CS%Coriolis_En_Dis) then
        ! Energy dissipating biased scheme, Hallberg 200x
        do j=js,je ; do I=Isq,Ieq
          if (q(I,J)*u(I,j,k) == 0.0) then
            temp1 = q(I,J) * ( (vh_max(i,j)+vh_max(i+1,j)) &
                             + (vh_min(i,j)+vh_min(i+1,j)) )*0.5
          elseif (q(I,J)*u(I,j,k) < 0.0) then
            temp1 = q(I,J) * (vh_max(i,j)+vh_max(i+1,j))
          else
            temp1 = q(I,J) * (vh_min(i,j)+vh_min(i+1,j))
          endif
          if (q(I,J-1)*u(I,j,k) == 0.0) then
            temp2 = q(I,J-1) * ( (vh_max(i,j-1)+vh_max(i+1,j-1)) &
                               + (vh_min(i,j-1)+vh_min(i+1,j-1)) )*0.5
          elseif (q(I,J-1)*u(I,j,k) < 0.0) then
            temp2 = q(I,J-1) * (vh_max(i,j-1)+vh_max(i+1,j-1))
          else
            temp2 = q(I,J-1) * (vh_min(i,j-1)+vh_min(i+1,j-1))
          endif
          CAu(I,j,k) = 0.25 * G%IdxCu(I,j) * (temp1 + temp2)
        enddo ; enddo
      else
        ! Energy conserving scheme, Sadourny 1975
        do j=js,je ; do I=Isq,Ieq
          CAu(I,j,k) = 0.25 * &
            ((q(I,J) * (vh(i+1,J,k) + vh(i,J,k))) + &
             (q(I,J-1) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
        enddo ; enddo
      endif
    elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        CAu(I,j,k) = 0.125 * (G%IdxCu(I,j) * (q(I,J) + q(I,J-1))) * &
                     ((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      enddo ; enddo
    elseif ((CS%Coriolis_Scheme == ARAKAWA_HSU90) .or. &
            (CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
            (CS%Coriolis_Scheme == AL_BLEND)) then
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
      do j=js,je ; do I=Isq,Ieq
        CAu(I,j,k) = (((a(I,j) * vh(i+1,J,k)) +  (c(I,j) * vh(i,J-1,k)))  + &
                      ((b(I,j) * vh(i,J,k)) +  (d(I,j) * vh(i+1,J-1,k)))) * G%IdxCu(I,j)
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == ROBUST_ENSTRO) then
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_v that should be returned by the
      !       continuity solver. AJA
      do j=js,je ; do I=Isq,Ieq
        Heff1 = abs(vh(i,J,k) * G%IdxCv(i,J)) / (eps_vel+abs(v(i,J,k)))
        Heff1 = max(Heff1, min(h(i,j,k),h(i,j+1,k)))
        Heff1 = min(Heff1, max(h(i,j,k),h(i,j+1,k)))
        Heff2 = abs(vh(i,J-1,k) * G%IdxCv(i,J-1)) / (eps_vel+abs(v(i,J-1,k)))
        Heff2 = max(Heff2, min(h(i,j-1,k),h(i,j,k)))
        Heff2 = min(Heff2, max(h(i,j-1,k),h(i,j,k)))
        Heff3 = abs(vh(i+1,J,k) * G%IdxCv(i+1,J)) / (eps_vel+abs(v(i+1,J,k)))
        Heff3 = max(Heff3, min(h(i+1,j,k),h(i+1,j+1,k)))
        Heff3 = min(Heff3, max(h(i+1,j,k),h(i+1,j+1,k)))
        Heff4 = abs(vh(i+1,J-1,k) * G%IdxCv(i+1,J-1)) / (eps_vel+abs(v(i+1,J-1,k)))
        Heff4 = max(Heff4, min(h(i+1,j-1,k),h(i+1,j,k)))
        Heff4 = min(Heff4, max(h(i+1,j-1,k),h(i+1,j,k)))
        if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
          CAu(I,j,k) = 0.5*(abs_vort(I,J)+abs_vort(I,J-1)) * &
                       ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) ) /  &
                       (h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) * G%IdxCu(I,j)
        elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
          VHeff = ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) )
          QVHeff = 0.5*( ((abs_vort(I,J)+abs_vort(I,J-1))*VHeff) &
                       - ((abs_vort(I,J)-abs_vort(I,J-1))*abs(VHeff)) )
          CAu(I,j,k) = (QVHeff / ( h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) ) * G%IdxCu(I,j)
        endif
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
        ! check whether there is masked land points in the stencil
        third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                       G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

        fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3)
        seventh_order = fifth_order * G%mask2dCu(I,j-4) * G%mask2dCu(I,j+4)


        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
          ! all values are valid, we use seventh order reconstruction
          u_q8(:) = (u(I,j-4:j+3,k) + u(I,j-3:j+4,k)) * 0.5
          call weno_seven_h_weight_reconstruction(abs_vort(I,J-4:J+3), &
                                         h_q(I,J-4:J+3), &
                                         u_q8, &
                                         GV%H_subroundoff, v_u, q_u, cs%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)

        elseif (fifth_order == 1) then
          ! all values are valid, we use fifth order reconstruction
          u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
          call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2), &
                                        h_q(I,J-3:J+2), &
                                        u_q6, &
                                        GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)

        elseif (third_order == 1) then
          ! only the middle values are valid, we use third order reconstruction
          u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
          call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1), &
                                         h_q(I,J-2:J+1), &
                                         u_q4, &
                                         GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)
        else ! Upwind first order
          if (v_u>0.) then
              q_u = q(I,J-1)
          else
              q_u = q(I,J)
          endif
          CAu(I,j,k) = (q_u * v_u)

        endif
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
        third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                       G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

        fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3)

        if (fifth_order == 1) then
          ! all values are valid, we use fifth order reconstruction
          u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
          call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2), &
                                        h_q(I,J-3:J+2), &
                                        u_q6, &
                                        GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)

        elseif (third_order == 1) then
          ! only the middle values are valid, we use third order reconstruction
          u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
          call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1), &
                                         h_q(I,J-2:J+1), &
                                         u_q4, &
                                         GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)

        else ! Upwind first order
          if (v_u>0.) then
              q_u = q(I,J-1)
          else
              q_u = q(I,J)
          endif
          CAu(I,j,k) = (q_u * v_u)
        endif
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
        third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                       G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))


        if (third_order == 1) then
          ! only the middle values are valid, we use third order reconstruction
          u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
          call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1), &
                                         h_q(I,J-2:J+1), &
                                         u_q4, &
                                         GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
          CAu(I,j,k) = (q_u * v_u)

        else ! Upwind first order
          if (v_u>0.) then
              q_u = q(I,J-1)
          else
              q_u = q(I,J)
          endif
          CAu(I,j,k) = (q_u * v_u)
        endif
      enddo ; enddo
    endif
    ! Add in the additional terms with Arakawa & Lamb.
    if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
        (CS%Coriolis_Scheme == AL_BLEND)) then ; do j=js,je ; do I=Isq,Ieq
      CAu(I,j,k) = CAu(I,j,k) + &
            ((ep_u(i,j)*uh(I-1,j,k)) - (ep_u(i+1,j)*uh(I+1,j,k))) * G%IdxCu(I,j)
    enddo ; enddo ; endif

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        ! Computing the diagnostic Stokes contribution to CAu
        do j=js,je ; do I=Isq,Ieq
          CAuS(I,j,k) = 0.25 * &
                ((qS(I,J) * (vh(i+1,J,k) + vh(i,J,k))) + &
                 (qS(I,J-1) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
        enddo ; enddo
      endif
    endif

    if (CS%bound_Coriolis) then
      do j=js,je ; do I=Isq,Ieq
        fv1 = abs_vort(I,J) * v(i+1,J,k)
        fv2 = abs_vort(I,J) * v(i,J,k)
        fv3 = abs_vort(I,J-1) * v(i+1,J-1,k)
        fv4 = abs_vort(I,J-1) * v(i,J-1,k)

        max_fv = max(fv1, fv2, fv3, fv4)
        min_fv = min(fv1, fv2, fv3, fv4)

        CAu(I,j,k) = min(CAu(I,j,k), max_fv)
        CAu(I,j,k) = max(CAu(I,j,k), min_fv)
      enddo ; enddo
    endif

    ! Term - d(KE)/dx.
    do j=js,je ; do I=Isq,Ieq
      CAu(I,j,k) = CAu(I,j,k) - KEx(I,j)
    enddo ; enddo

    if (associated(AD%gradKEu)) then
      do j=js,je ; do I=Isq,Ieq
        AD%gradKEu(I,j,k) = -KEx(I,j)
      enddo ; enddo
    endif

    ! Calculate the tendencies of meridional velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAv = - q * uh - d(KE)/dy.
    if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (CS%Coriolis_En_Dis) then
        ! Energy dissipating biased scheme, Hallberg 200x
        do J=Jsq,Jeq ; do i=is,ie
          if (q(I-1,J)*v(i,J,k) == 0.0) then
            temp1 = q(I-1,J) * ( (uh_max(i-1,j)+uh_max(i-1,j+1)) &
                               + (uh_min(i-1,j)+uh_min(i-1,j+1)) )*0.5
          elseif (q(I-1,J)*v(i,J,k) > 0.0) then
            temp1 = q(I-1,J) * (uh_max(i-1,j)+uh_max(i-1,j+1))
          else
            temp1 = q(I-1,J) * (uh_min(i-1,j)+uh_min(i-1,j+1))
          endif
          if (q(I,J)*v(i,J,k) == 0.0) then
            temp2 = q(I,J) * ( (uh_max(i,j)+uh_max(i,j+1)) &
                             + (uh_min(i,j)+uh_min(i,j+1)) )*0.5
          elseif (q(I,J)*v(i,J,k) > 0.0) then
            temp2 = q(I,J) * (uh_max(i,j)+uh_max(i,j+1))
          else
            temp2 = q(I,J) * (uh_min(i,j)+uh_min(i,j+1))
          endif
          CAv(i,J,k) = -0.25 * G%IdyCv(i,J) * (temp1 + temp2)
        enddo ; enddo
      else
        ! Energy conserving scheme, Sadourny 1975
        do J=Jsq,Jeq ; do i=is,ie
          CAv(i,J,k) = - 0.25* &
              ((q(I-1,J)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
               (q(I,J)*(uh(I,j,k) + uh(I,j+1,k)))) * G%IdyCv(i,J)
        enddo ; enddo
      endif
    elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        CAv(i,J,k) = -0.125 * (G%IdyCv(i,J) * (q(I-1,J) + q(I,J))) * &
                     ((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))
      enddo ; enddo
    elseif ((CS%Coriolis_Scheme == ARAKAWA_HSU90) .or. &
            (CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
            (CS%Coriolis_Scheme == AL_BLEND)) then
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
      do J=Jsq,Jeq ; do i=is,ie
        CAv(i,J,k) = - (((a(I-1,j)   * uh(I-1,j,k)) + &
                         (c(I,j+1)   * uh(I,j+1,k)))  &
                      + ((b(I,j)     * uh(I,j,k)) +   &
                         (d(I-1,j+1) * uh(I-1,j+1,k)))) * G%IdyCv(i,J)
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == ROBUST_ENSTRO) then
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_u that should be returned by the
      !       continuity solver. AJA
      do J=Jsq,Jeq ; do i=is,ie
        Heff1 = abs(uh(I,j,k) * G%IdyCu(I,j)) / (eps_vel+abs(u(I,j,k)))
        Heff1 = max(Heff1, min(h(i,j,k),h(i+1,j,k)))
        Heff1 = min(Heff1, max(h(i,j,k),h(i+1,j,k)))
        Heff2 = abs(uh(I-1,j,k) * G%IdyCu(I-1,j)) / (eps_vel+abs(u(I-1,j,k)))
        Heff2 = max(Heff2, min(h(i-1,j,k),h(i,j,k)))
        Heff2 = min(Heff2, max(h(i-1,j,k),h(i,j,k)))
        Heff3 = abs(uh(I,j+1,k) * G%IdyCu(I,j+1)) / (eps_vel+abs(u(I,j+1,k)))
        Heff3 = max(Heff3, min(h(i,j+1,k),h(i+1,j+1,k)))
        Heff3 = min(Heff3, max(h(i,j+1,k),h(i+1,j+1,k)))
        Heff4 = abs(uh(I-1,j+1,k) * G%IdyCu(I-1,j+1)) / (eps_vel+abs(u(I-1,j+1,k)))
        Heff4 = max(Heff4, min(h(i-1,j+1,k),h(i,j+1,k)))
        Heff4 = min(Heff4, max(h(i-1,j+1,k),h(i,j+1,k)))
        if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
          CAv(i,J,k) = - 0.5*(abs_vort(I,J)+abs_vort(I-1,J)) * &
                         ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                          (uh(I-1,j  ,k)+uh(I  ,j+1,k)) ) /    &
                      (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
        elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
          UHeff = ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                   (uh(I-1,j  ,k)+uh(I  ,j+1,k)) )
          QUHeff = 0.5*( ((abs_vort(I,J)+abs_vort(I-1,J))*UHeff) &
                       - ((abs_vort(I,J)-abs_vort(I-1,J))*abs(UHeff)) )
          CAv(i,J,k) = - QUHeff / &
                       (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
        endif
      enddo ; enddo
    ! Calculate the tendencies of meridional velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAv = - q * uh - d(KE)/dy.
    elseif (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

        ! check whether there is any masked land values within the stencils
        third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                       G%mask2dCv(i+2,J))
        fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J)
        seventh_order = fifth_order * G%mask2dCv(i-4,J) * G%mask2dCv(i+4,J)



        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
          v_q8(:) = (v(i-4:i+3,J,k) + v(i-3:i+4,J,k)) * 0.5
          ! all values are valid, we use seventh order reconstruction
          call weno_seven_h_weight_reconstruction(abs_vort(I-4:I+3,J), &
                                         h_q(I-4:I+3,J), &
                                         v_q8, &
                                         GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)

        elseif (fifth_order == 1) then
          v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
          ! all values are valid, we use fifth order reconstruction
          call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J), &
                                        h_q(I-3:I+2,J), &
                                        v_q6, &
                                        GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)

        elseif (third_order == 1) then
          v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
          call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J), &
                                                 h_q(I-2:I+1,J), &
                                                 v_q4, &
                                                 GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)
        else ! Upwind first order!
          if (u_v>0.) then
              q_v = q(I-1,J)
          else
              q_v = q(I,J)
          endif
          CAv(i,J,k) = - (q_v * u_v)
        endif

      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

        third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                       G%mask2dCv(i+2,J))
        fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J)


        ! compute the masking to make sure that inland values are not used
        if (fifth_order == 1) then
          v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
          ! all values are valid, we use fifth order reconstruction
          call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J), &
                                        h_q(I-3:I+2,J), &
                                        v_q6, &
                                        GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)

        elseif (third_order == 1) then
          v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
          call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J), &
                                                 h_q(I-2:I+1,J), &
                                                 v_q4, &
                                                 GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)

        else
          if (u_v>0.) then
              q_v = q(I-1,J)
          else
              q_v = q(I,J)
          endif
          CAv(i,J,k) = - (q_v * u_v)
        endif

      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

        third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                       G%mask2dCv(i+2,J))


        ! compute the masking to make sure that inland values are not used
        if (third_order == 1) then
          v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
          call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J), &
                                                 h_q(I-2:I+1,J), &
                                                 v_q4, &
                                                 GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
          CAv(i,J,k) = - (q_v * u_v)

        else
          if (u_v>0.) then
              q_v = q(I-1,J)
          else
              q_v = q(I,J)
          endif
          CAv(i,J,k) = - (q_v * u_v)
        endif

      enddo ; enddo
    endif
    ! Add in the additonal terms with Arakawa & Lamb.
    if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
        (CS%Coriolis_Scheme == AL_BLEND)) then ; do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) + &
            ((ep_v(i,j)*vh(i,J-1,k)) - (ep_v(i,j+1)*vh(i,J+1,k))) * G%IdyCv(i,J)
    enddo ; enddo ; endif

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        ! Computing the diagnostic Stokes contribution to CAv
        do J=Jsq,Jeq ; do i=is,ie
          CAvS(i,J,k) = 0.25 * &
                ((qS(I,J) * (uh(I,j+1,k) + uh(I,j,k))) + &
                 (qS(I-1,J) * (uh(I-1,j,k) + uh(I-1,j+1,k)))) * G%IdyCv(i,J)
        enddo ; enddo
      endif
    endif

    if (CS%bound_Coriolis) then
      do J=Jsq,Jeq ; do i=is,ie
        fu1 = -abs_vort(I,J) * u(I,j+1,k)
        fu2 = -abs_vort(I,J) * u(I,j,k)
        fu3 = -abs_vort(I-1,J) * u(I-1,j+1,k)
        fu4 = -abs_vort(I-1,J) * u(I-1,j,k)

        max_fu = max(fu1, fu2, fu3, fu4)
        min_fu = min(fu1, fu2, fu3, fu4)

        CAv(I,j,k) = min(CAv(I,j,k), max_fu)
        CAv(I,j,k) = max(CAv(I,j,k), min_fu)
      enddo ; enddo
    endif

    ! Term - d(KE)/dy.
    do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) - KEy(i,J)
    enddo ; enddo
    if (associated(AD%gradKEv)) then
      do J=Jsq,Jeq ; do i=is,ie
        AD%gradKEv(i,J,k) = -KEy(i,J)
      enddo ; enddo
    endif

    if (associated(AD%rv_x_u) .or. associated(AD%rv_x_v)) then
      ! Calculate the Coriolis-like acceleration due to relative vorticity.
      if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
        if (associated(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = - 0.25* &
              ((q2(I-1,j)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
               (q2(I,j)*(uh(I,j,k) + uh(I,j+1,k)))) * G%IdyCv(i,J)
          enddo ; enddo
        endif

        if (associated(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = 0.25 * &
              ((q2(I,j) * (vh(i+1,J,k) + vh(i,J,k))) + &
               (q2(I,j-1) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
          enddo ; enddo
        endif
      else
        if (associated(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = -G%IdyCv(i,J) * C1_12 * &
              (((((q2(I,J) + q2(I-1,J-1)) + q2(I-1,J)) * uh(I-1,j,k)) + &
                (((q2(I-1,J) + q2(I,J+1)) + q2(I,J)) * uh(I,j+1,k))) + &
               ((((q2(I-1,J) + q2(I,J-1)) + q2(I,J)) * uh(I,j,k))+ &
                (((q2(I,J) + q2(I-1,J+1)) + q2(I-1,J)) * uh(I-1,j+1,k))))
          enddo ; enddo
        endif

        if (associated(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = G%IdxCu(I,j) * C1_12 * &
              (((((q2(I+1,J) + q2(I,J-1)) + q2(I,J)) * vh(i+1,J,k)) + &
                (((q2(I-1,J-1) + q2(I,J)) + q2(I,J-1)) * vh(i,J-1,k))) + &
               ((((q2(I-1,J) + q2(I,J-1)) + q2(I,J)) * vh(i,J,k)) + &
                (((q2(I+1,J-1) + q2(I,J)) + q2(I,J-1)) * vh(i+1,J-1,k))))
          enddo ; enddo
        endif
      endif
    endif

  enddo ! k-loop.

  ! Here the various Coriolis-related derived quantities are offered for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_rv > 0) call post_data(CS%id_rv, RV, CS%diag)
    if (CS%id_PV > 0) call post_data(CS%id_PV, PV, CS%diag)
    if (CS%id_gKEu>0) call post_data(CS%id_gKEu, AD%gradKEu, CS%diag)
    if (CS%id_gKEv>0) call post_data(CS%id_gKEv, AD%gradKEv, CS%diag)
    if (CS%id_rvxu > 0) call post_data(CS%id_rvxu, AD%rv_x_u, CS%diag)
    if (CS%id_rvxv > 0) call post_data(CS%id_rvxv, AD%rv_x_v, CS%diag)
    if (Stokes_VF) then
      if (CS%id_CAuS > 0) call post_data(CS%id_CAuS, CAuS, CS%diag)
      if (CS%id_CAvS > 0) call post_data(CS%id_CAvS, CAvS, CS%diag)
    endif

    ! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_gKEu etc. are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    ! if (CS%id_hf_gKEu > 0) call post_product_u(CS%id_hf_gKEu, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_gKEv > 0) call post_product_v(CS%id_hf_gKEv, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    ! if (CS%id_hf_rvxv > 0) call post_product_u(CS%id_hf_rvxv, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_rvxu > 0) call post_product_v(CS%id_hf_rvxu, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_hf_gKEu_2d > 0) call post_product_sum_u(CS%id_hf_gKEu_2d, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_gKEv_2d > 0) call post_product_sum_v(CS%id_hf_gKEv_2d, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    if (CS%id_intz_gKEu_2d > 0) call post_product_sum_u(CS%id_intz_gKEu_2d, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_gKEv_2d > 0) call post_product_sum_v(CS%id_intz_gKEv_2d, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_hf_rvxv_2d > 0) call post_product_sum_u(CS%id_hf_rvxv_2d, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_rvxu_2d > 0) call post_product_sum_v(CS%id_hf_rvxu_2d, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_h_gKEu > 0) call post_product_u(CS%id_h_gKEu, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_gKEv > 0) call post_product_v(CS%id_h_gKEv, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)
    if (CS%id_h_rvxv > 0) call post_product_u(CS%id_h_rvxv, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_rvxu > 0) call post_product_v(CS%id_h_rvxu, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_intz_rvxv_2d > 0) call post_product_sum_u(CS%id_intz_rvxv_2d, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_rvxu_2d > 0) call post_product_sum_v(CS%id_intz_rvxu_2d, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)
  endif

end procedure CorAdCalc
module procedure gradKE
  real :: um, up, vm, vp         ! Temporary variables [L T-1 ~> m s-1].
  real :: um2, up2, vm2, vp2     ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: um2a, up2a, vm2a, vp2a ! Temporary variables [L4 T-2 ~> m4 s-2].
  real :: third_order_u, third_order_v  ! Product of mask values to determine the boundary
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  real, parameter     :: C1_12 = 1.0/12.0   ! The ratio of 1/12 [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB


  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
  if (CS%KE_Scheme == KE_ARAKAWA) then
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      KE(i,j) = ( ( (G%areaCu( I ,j)*(u( I ,j)*u( I ,j))) + &
                    (G%areaCu(I-1,j)*(u(I-1,j)*u(I-1,j))) ) + &
                  ( (G%areaCv(i, J )*(v(i, J )*v(i, J ))) + &
                    (G%areaCv(i,J-1)*(v(i,J-1)*v(i,J-1))) ) )*0.25*G%IareaT(i,j)
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_SIMPLE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme which does not take into account any geometric factors
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j) + ABS( u(I-1,j) ) ) ; up2 = up*up
      um = 0.5*( u( I ,j) - ABS( u( I ,j) ) ) ; um2 = um*um
      vp = 0.5*( v(i,J-1) + ABS( v(i,J-1) ) ) ; vp2 = vp*vp
      vm = 0.5*( v(i, J ) - ABS( v(i, J ) ) ) ; vm2 = vm*vm
      KE(i,j) = ( max(up2,um2) + max(vp2,vm2) ) *0.5
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j) + ABS( u(I-1,j) ) ) ; up2a = up*up*G%areaCu(I-1,j)
      um = 0.5*( u( I ,j) - ABS( u( I ,j) ) ) ; um2a = um*um*G%areaCu( I ,j)
      vp = 0.5*( v(i,J-1) + ABS( v(i,J-1) ) ) ; vp2a = vp*vp*G%areaCv(i,J-1)
      vm = 0.5*( v(i, J ) - ABS( v(i, J ) ) ) ; vm2a = vm*vm*G%areaCv(i, J )
      KE(i,j) = ( max(um2a,up2a) + max(vm2a,vp2a) )*0.5*G%IareaT(i,j)
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_UP3) then
    ! The following discretization of KE is based on the one-dimensional third-order
    ! upwind scheme which does not take horizontal grid factors into account
    if (CS%KE_use_limiter) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        ! compute the masking to make sure that inland values are not used
        third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j)* &
                       G%mask2dCu(I,j) * G%mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j) + u(I,j)) - (u(I-2,j) + u(I+1,j))) * C1_12
          call UP3_Koren_limiter_reconstruction(u(I-2:I+1,j), up, um)
        else
          up = (u(I-1,j) + u(I,j))*0.5
          if (up>0.) then
            um = u(I-1,j)
          elseif (up<0.) then
            um = u(I,j)
          else
            um = up
          endif
        endif

        third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1)* &
                       G%mask2dCv(i,J) * G%mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1) + v(i,J)) - (v(i,J-2) + v(i,J+1))) * C1_12
          call UP3_Koren_limiter_reconstruction(v(i,J-2:J+1), vp, vm)
        else
          vp = (v(i,J-1) + v(i,J))*0.5
          if (vp>0.) then
            vm = v(i,J-1)
          elseif (vp<0.) then
            vm = v(i,J)
          else
            vm = vp
          endif
        endif

        KE(i,j) = ( (um*um) + (vm*vm) )*0.5
      enddo ; enddo
    else
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        ! compute the masking to make sure that inland values are not used
        third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j)* &
                       G%mask2dCu(I,j) * G%mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j) + u(I,j)) - (u(I-2,j) + u(I+1,j))) * C1_12
          call UP3_reconstruction(u(I-2:I+1,j), up, um)
        else
          up = (u(I-1,j) + u(I,j))*0.5
          if (up>0.) then
            um = u(I-1,j)
          elseif (up<0.) then
            um = u(I,j)
          else
            um = up
          endif
        endif

        third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1)* &
                       G%mask2dCv(i,J) * G%mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1) + v(i,J)) - (v(i,J-2) + v(i,J+1))) * C1_12
          call UP3_reconstruction(v(i,J-2:J+1), vp, vm)
        else
          vp = (v(i,J-1) + v(i,J))*0.5
          if (vp>0.) then
            vm = v(i,J-1)
          elseif (vp<0.) then
            vm = v(i,J)
          else
            vm = vp
          endif
        endif

        KE(i,j) = ( (um*um) + (vm*vm) )*0.5
      enddo ; enddo
    endif
  endif

  ! Term - d(KE)/dx.
  do j=js,je ; do I=Isq,Ieq
    KEx(I,j) = (KE(i+1,j) - KE(i,j)) * G%IdxCu_OBCmask(I,j)
  enddo ; enddo

  ! Term - d(KE)/dy.
  do J=Jsq,Jeq ; do i=is,ie
    KEy(i,J) = (KE(i,j+1) - KE(i,j)) * G%IdyCv_OBCmask(i,J)
  enddo ; enddo

end procedure gradKE
module procedure UP3_reconstruction
  real, parameter :: C1_6 = 1.0/6.0       ! The ratio of 1/6 [nondim]
  if (u>0.) then
    qr = ((2.*q4(3) + 5.*q4(2)) - q4(1)) * C1_6
  else
    qr = ((2.*q4(2) + 5.*q4(3)) - q4(4)) * C1_6
  endif

end procedure UP3_reconstruction
module procedure UP3_Koren_limiter_reconstruction
  real                :: theta            ! Ratio of gradient [nondim]
  real                :: psi              ! Limiter function [nondim]
  real, parameter     :: C1_3 = 1.0/3.0   ! The ratio of 1/3 [nondim]
  real, parameter     :: C1_6 = 1.0/6.0   ! The ratio of 1/6 [nondim]
  if (u>0.) then
    if (q4(3) == q4(2)) then
      qr = q4(2)
    else
      theta = (q4(2) - q4(1))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta)) ! limiter introduced by Koren (1993)
      qr = q4(2) + psi*(q4(3) - q4(2))
    endif
  else
    if (q4(3) == q4(2)) then
      qr = q4(3)
    else
      theta = (q4(4) - q4(3))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta))
      qr = q4(3) + psi*(q4(2) - q4(3))
    endif
  endif

end procedure UP3_Koren_limiter_reconstruction
module procedure fac_fn
  fac = 1.0e40 ; if (abs(b) > 1.0e-20*tau) fac = (1 + tau / b)**2

end procedure fac_fn
module procedure weno_three_h_weight_reconstruction
    real :: vr                            ! Reconstruction of hq [A ~> a]
    real :: hr                            ! Reconstruction of h [L ~> m]
    real :: c0, c1                        ! Intermediate reconstruction of q [A ~> a]
    real :: d0, d1                        ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1                        ! Smoothness indicator [A ~> a]
    real :: tau                           ! Difference of smoothness indicator [A ~> a]
    real :: w0, w1                        ! Weights [nondim]
    real :: s                             ! Temporary variables [nondim]
    real, parameter :: C2_3 = 2.0/3.0     ! The ratio of 2/3 [nondim]
    real, parameter :: C1_3 = 1.0/3.0     ! The ratio of 1/3 [nondim]
    if (u>0.) then
      call weno_three_reconstruction_0(q4(2:3), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(1:2), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(2:3), d0)
      call weno_three_reconstruction_1(h4(1:2), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      endif
    else
      call weno_three_reconstruction_0(q4(3:2:-1), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(4:3:-1), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(3:2:-1), d0)
      call weno_three_reconstruction_1(h4(4:3:-1), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      endif
    endif

    tau = abs(b0-b1)
    w0  = C2_3 * fac_fn(tau, b0)
    w1  = C1_3 * fac_fn(tau, b1)

    s = 1. / (w0 + w1)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s

    vr = (w0 * c0) + (w1 * c1)
    hr = (w0 * d0) + (w1 * d1)
!    vr = min(max(q4(3), q4(2)), vr) ; vr = max(min(q4(3), q4(2)), vr) !Impose a monotonicity limiter
    hr = min(max(h4(3), h4(2)), hr) ; hr = max(min(h4(3), h4(2)), hr) ! A monotonicity limiter

    qr = vr / max(hr, h_tiny)

end procedure weno_three_h_weight_reconstruction
module procedure weno_three_weight
    w0 = (q2(1) - q2(2))**2

end procedure weno_three_weight
module procedure weno_three_reconstruction_0
    w0 = (q2(1) + q2(2)) * 0.5

end procedure weno_three_reconstruction_0
module procedure weno_three_reconstruction_1
    w0 = (- q2(1) + 3 * q2(2)) * 0.5

end procedure weno_three_reconstruction_1
module procedure weno_five_h_weight_reconstruction
    real :: vr                                    ! Reconstruction of hq [A ~> a]
    real :: hr                                    ! Reconstruction of h [L ~> m]
    real :: c0, c1, c2                            ! Intermediate reconstruction of hq[A ~> a]
    real :: d0, d1, d2                            ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1, b2                            ! Smoothness indicator [A ~> a]
    real :: tau                                   ! Difference of smoothness indicators [A ~> a]
    real :: w0, w1, w2                            ! Weights [nondim]
    real :: s                                     ! Temporary variables [nondim]
    real, parameter :: C3_10 = 3.0/10.0           ! The ratio of 3/10 [nondim]
    real, parameter :: C3_5 = 3.0/5.0             ! The ratio of 3/5 [nondim]
    real, parameter :: C1_10 = 1.0/10.0           ! The ratio of 1/10 [nondim]
    if (u>0.) then
      call weno_five_reconstruction_0(q6(3:5), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(2:4), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(1:3), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(3:5), d0)
      call weno_five_reconstruction_1(h6(2:4), d1)
      call weno_five_reconstruction_2(h6(1:3), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(3:5), b0)   ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(2:4), b1)   ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(1:3), b2)   ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(3:5), b0)
        call weno_five_weight_1(q6(2:4), b1)
        call weno_five_weight_2(q6(1:3), b2)
      endif
    else
      call weno_five_reconstruction_0(q6(4:2:-1), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(5:3:-1), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(6:4:-1), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(4:2:-1), d0)
      call weno_five_reconstruction_1(h6(5:3:-1), d1)
      call weno_five_reconstruction_2(h6(6:4:-1), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(4:2:-1), b0) ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(5:3:-1), b1) ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(6:4:-1), b2) ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(4:2:-1), b0)
        call weno_five_weight_1(q6(5:3:-1), b1)
        call weno_five_weight_2(q6(6:4:-1), b2)
      endif
    endif

    tau = abs(b0 - b2)
    w0  = C3_10 * fac_fn(tau, b0)
    w1  = C3_5  * fac_fn(tau, b1)
    w2  = C1_10 * fac_fn(tau, b2)

    s = 1. / ((w0 + w1) + w2)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s
    w2 = w2 * s

    vr = ((w0 * c0) + (w1 * c1)) + (w2 * c2)
    hr = ((w0 * d0) + (w1 * d1)) + (w2 * d2)
!    vr = min(max(q6(3), q6(4)), vr) ; vr = max(min(q6(3), q6(4)), vr) !Impose a monotonicity limiter
    hr = min(max(h6(3), h6(4)), hr) ; hr = max(min(h6(3), h6(4)), hr) !Impose a monotonicity limiter

    qr = vr / max(hr, h_tiny)

end procedure weno_five_h_weight_reconstruction
module procedure weno_five_weight_0
  w0 = (q3(1) * ((10 * q3(1) - 31 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 19 * q3(3))) + 4 * (q3(3) * q3(3)))

end procedure weno_five_weight_0
module procedure weno_five_weight_1
  w1 = (q3(1) * ((4 * q3(1) - 13 * q3(2)) + 5 * q3(3))) + &
       ((q3(2) * (13 * q3(2) - 13 * q3(3))) + 4 * (q3(3) * q3(3)))

end procedure weno_five_weight_1
module procedure weno_five_weight_2
  w2 = (q3(1) * ((4 * q3(1) - 19 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 31 * q3(3))) + 10 * (q3(3) * q3(3)))

end procedure weno_five_weight_2
module procedure weno_five_reconstruction_0
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]
  p0 = ((2*q3(1) + 5*q3(2)) - q3(3)) * C1_6

end procedure weno_five_reconstruction_0
module procedure weno_five_reconstruction_1
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]
  p1 = ((-q3(1) + 5*q3(2)) + 2*q3(3)) * C1_6

end procedure weno_five_reconstruction_1
module procedure weno_five_reconstruction_2
  real, parameter :: C1_6 = 1.0/6.0  ! One sixth [nondim]
  p2 = ((2*q3(1) - 7*q3(2)) + 11*q3(3)) * C1_6

end procedure weno_five_reconstruction_2
module procedure weno_seven_h_weight_reconstruction
  real :: vr                  ! Reconstruction of hq [A ~> a]
  real :: hr                  ! Reconstruction of h [L ~> m]
  real :: c0, c1, c2, c3      ! Intermediate reconstruction of hq [A ~> a]
  real :: d0, d1, d2, d3      ! Intermediate reconstruction of h [L ~> m]
  real :: b0, b1, b2, b3      ! Smoothness indicator [A ~> a]
  real :: tau                 ! Difference of smoothness indicators [A ~> a]
  real :: w0, w1, w2, w3      ! Weights [nondim]
  real :: s                   ! Temporary variables [nondim]
  real, parameter :: C4_35 = 4.0/35.0 ! The ratio of 4/35 [nondim]
  real, parameter :: C18_35 = 18.0/35.0 ! The ratio of 18/35 [nondim]
  real, parameter :: C12_35 = 12.0/35.0 ! The ratio of 12/35 [nondim]
  real, parameter :: C1_35 = 1.0/35.0   ! The ratio of 1/35 [nondim]
  if (u>0.) then
    call weno_seven_reconstruction_0(q8(4:7), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(3:6), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(2:5), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(1:4), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(4:7), d0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(h8(3:6), d1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(h8(2:5), d2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(h8(1:4), d3) ! Reconstruction in the first upwind stencil
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(4:7), b0)       ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(3:6), b1)       ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(2:5), b2)       ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(1:4), b3)       ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(4:7), b0)
      call weno_seven_weight_1(q8(3:6), b1)
      call weno_seven_weight_2(q8(2:5), b2)
      call weno_seven_weight_3(q8(1:4), b3)
    endif
  else
    call weno_seven_reconstruction_0(q8(5:2:-1), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(6:3:-1), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(7:4:-1), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(8:5:-1), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(5:2:-1), d0)
    call weno_seven_reconstruction_1(h8(6:3:-1), d1)
    call weno_seven_reconstruction_2(h8(7:4:-1), d2)
    call weno_seven_reconstruction_3(h8(8:5:-1), d3)
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(5:2:-1), b0)    ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(6:3:-1), b1)    ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(7:4:-1), b2)    ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(8:5:-1), b3)    ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(5:2:-1), b0)
      call weno_seven_weight_1(q8(6:3:-1), b1)
      call weno_seven_weight_2(q8(7:4:-1), b2)
      call weno_seven_weight_3(q8(8:5:-1), b3)
    endif
  endif

  tau = abs((b0 - b3) + 3 * (b1 - b2))
  w0  = C4_35  * fac_fn(tau, b0)
  w1  = C18_35 * fac_fn(tau, b1)
  w2  = C12_35 * fac_fn(tau, b2)
  w3  = C1_35  * fac_fn(tau, b3)

  s = 1. / ((w0 + w1) + (w2 + w3))
  w0 = w0 * s   ! Weights of the stencils
  w1 = w1 * s
  w2 = w2 * s
  w3 = w3 * s

  vr = ((w0 * c0) + (w1 * c1)) + ((w2 * c2) + (w3 * c3))
  hr = ((w0 * d0) + (w1 * d1)) + ((w2 * d2) + (w3 * d3))

!  vr = min(max(q4, q5), vr) ; vr = max(min(q4, q5), vr)
  hr = min(max(h8(4), h8(5)), hr) ; hr = max(min(h8(4), h8(5)), hr) ! Impose a monotonicity limiter

  qr = vr / max(hr, h_tiny)

end procedure weno_seven_h_weight_reconstruction
module procedure weno_seven_weight_0
  w0 = ((q4(1) * ((2.107 * q4(1) - 9.402 * q4(2)) + (7.042 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((11.003 * q4(2) - 17.246 * q4(3)) + 4.642 * q4(4)))) + &
     ((q4(3) * (7.043 * q4(3) - 3.882 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end procedure weno_seven_weight_0
module procedure weno_seven_weight_1
  w1 = ((q4(1) * ((0.547 * q4(1) - 2.522 * q4(2)) + (1.922 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((3.443 * q4(2) - 5.966 * q4(3)) + 1.602 * q4(4)))) + &
     ((q4(3) * (2.843 * q4(3) - 1.642 * q4(4))) + 0.267 * (q4(4) * q4(4)))

end procedure weno_seven_weight_1
module procedure weno_seven_weight_2
  w2 = ((q4(1) * ((0.267 * q4(1) - 1.642 * q4(2)) + (1.602 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((2.843 * q4(2) - 5.966 * q4(3)) + 1.922 * q4(4)))) + &
     ((q4(3) * (3.443 * q4(3) - 2.522 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end procedure weno_seven_weight_2
module procedure weno_seven_weight_3
  w3 = ((q4(1) * ((0.547  * q4(1) - 3.882 * q4(2)) + (4.642 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((7.043 * q4(2) - 17.246 * q4(3)) + 7.042 * q4(4)))) + &
     ((q4(3) * (11.003 * q4(3) - 9.402 * q4(4))) + 2.107 * (q4(4) * q4(4)))

end procedure weno_seven_weight_3
module procedure weno_seven_reconstruction_0
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]
  p0 = (((6 * q4(1) + 26 * q4(2)) - 10 * q4(3)) + 2 * q4(4)) * C1_24

end procedure weno_seven_reconstruction_0
module procedure weno_seven_reconstruction_1
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]
  p1 = (14 * (q4(2) + q4(3)) - 2 * (q4(1) + q4(4))) * C1_24

end procedure weno_seven_reconstruction_1
module procedure weno_seven_reconstruction_2
  real, parameter :: C1_24 = 1.0/24.0   ! One twenty fourth [nondim]
  p2 = (((2 * q4(1) - 10 * q4(2)) + 26 * q4(3)) + 6 * q4(4)) * C1_24

end procedure weno_seven_reconstruction_2
module procedure weno_seven_reconstruction_3
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]
  p3 = (((-6 * q4(1) + 26 * q4(2)) - 46 * q4(3)) + 50 * q4(4)) * C1_24

end procedure weno_seven_reconstruction_3
module procedure CoriolisAdv_stencil
  stencil = 2
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) stencil = 4
  if (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) stencil = 3

end procedure CoriolisAdv_stencil
module procedure CoriolisAdv_init
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_CoriolisAdv" ! This module's name.
  character(len=20)  :: tmpstr
  character(len=400) :: mesg
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mdl, "CORIOLIS_EN_DIS", CS%Coriolis_En_Dis, &
                 "If true, two estimates of the thickness fluxes are used "//&
                 "to estimate the Coriolis term, and the one that "//&
                 "dissipates energy relative to the other one is used.", &
                 default=.false.)

  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
  call get_param(param_file, mdl, "CORIOLIS_SCHEME", tmpstr, &
                 "CORIOLIS_SCHEME selects the discretization for the "//&
                 "Coriolis terms. Valid values are: \n"//&
                 "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n"//&
                 "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n"//&
                 "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n"//&
                 "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n"//&
                 "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n"//&
                 "\t                      Arakawa & Hsu and Sadourny energy \n"//&
                 "\t WENOVI5TH_PV_ENSTRO   - 5th-order WENO PV enstrophy \n"//&
                 "\t WENOVI3RD_PV_ENSTRO   - 3rd-order WENO PV enstrophy \n"//&
                 "\t WENOVI7TH_PV_ENSTRO  - 7th-order WENO PV enstrophy \n", &
                 default=SADOURNY75_ENERGY_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (SADOURNY75_ENERGY_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENERGY
    case (ARAKAWA_HSU_STRING)
      CS%Coriolis_Scheme = ARAKAWA_HSU90
    case (SADOURNY75_ENSTRO_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENSTRO
    case (ARAKAWA_LAMB_STRING)
      CS%Coriolis_Scheme = ARAKAWA_LAMB81
    case (AL_BLEND_STRING)
      CS%Coriolis_Scheme = AL_BLEND
    case (ROBUST_ENSTRO_STRING)
      CS%Coriolis_Scheme = ROBUST_ENSTRO
      CS%Coriolis_En_Dis = .false.
    case (WENOVI7TH_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi7th_PV_ENSTRO
    case (WENOVI5TH_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi5th_PV_ENSTRO
    case (WENOVI3RD_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi3rd_PV_ENSTRO
    case default
      call MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: Unrecognized setting "// &
            "#define CORIOLIS_SCHEME "//trim(tmpstr)//" found in input file.")
  end select

  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO .or. &
          CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO .or. CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
    call get_param(param_file, mdl, "WENO_VELOCITY_SMOOTH", CS%weno_velocity_smooth, &
            "If true, use velocity to compute weighting for WENO. ", &
                  default=.false.)
  endif

  if (CS%Coriolis_Scheme == AL_BLEND) then
    call get_param(param_file, mdl, "CORIOLIS_BLEND_WT_LIN", CS%wt_lin_blend, &
                 "A weighting value for the ratio of inverse thicknesses, "//&
                 "beyond which the blending between Sadourny Energy and "//&
                 "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME "//&
                 "is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.", &
                 units="nondim", default=0.125)
    call get_param(param_file, mdl, "CORIOLIS_BLEND_F_EFF_MAX", CS%F_eff_max_blend, &
                 "The factor by which the maximum effective Coriolis "//&
                 "acceleration from any point can be increased when "//&
                 "blending different discretizations with the "//&
                 "ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be "//&
                 "greater than 2.0 (the max value for Sadourny energy).", &
                 units="nondim", default=4.0)
    CS%wt_lin_blend = min(1.0, max(CS%wt_lin_blend,1e-16))
    if (CS%F_eff_max_blend < 2.0) call MOM_error(WARNING, "CoriolisAdv_init: "//&
           "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
  endif

  mesg = "If true, the Coriolis terms at u-points are bounded by "//&
         "the four estimates of (f+rv)v from the four neighboring "//&
         "v-points, and similarly at v-points."
  if (CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) then
    mesg = trim(mesg)//"  This option is "//&
                 "always effectively false with CORIOLIS_EN_DIS defined and "//&
                 "CORIOLIS_SCHEME set to "//trim(SADOURNY75_ENERGY_STRING)//"."
  else
    mesg = trim(mesg)//"  This option would "//&
                 "have no effect on the SADOURNY Coriolis scheme if it "//&
                 "were possible to use centered difference thickness fluxes."
  endif
  call get_param(param_file, mdl, "BOUND_CORIOLIS", CS%bound_Coriolis, mesg, &
                 default=.false.)
  if ((CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) .or. &
      (CS%Coriolis_Scheme == ROBUST_ENSTRO)) CS%bound_Coriolis = .false.

  ! Set KE_Scheme (selects discretization of KE)
  call get_param(param_file, mdl, "KE_SCHEME", tmpstr, &
                 "KE_SCHEME selects the discretization for acceleration "//&
                 "due to the kinetic energy gradient. Valid values are: \n"//&
                 "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV, KE_UP3", &
                 default=KE_ARAKAWA_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (KE_ARAKAWA_STRING); CS%KE_Scheme = KE_ARAKAWA
    case (KE_SIMPLE_GUDONOV_STRING); CS%KE_Scheme = KE_SIMPLE_GUDONOV
    case (KE_GUDONOV_STRING); CS%KE_Scheme = KE_GUDONOV
    case (KE_UP3_STRING); CS%KE_Scheme = KE_UP3
    case default
      call MOM_mesg('CoriolisAdv_init: KE_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
               "#define KE_SCHEME "//trim(tmpstr)//" in input file is invalid.")
  end select

  if (CS%KE_Scheme == KE_UP3) then
    call get_param(param_file, mdl, "KE_USE_LIMITER", CS%KE_use_limiter, &
            "If true, use Koren limiter for KE_UP3 scheme", &
                  default=.True.)
  endif

  ! Set PV_Adv_Scheme (selects discretization of PV advection)
  call get_param(param_file, mdl, "PV_ADV_SCHEME", tmpstr, &
                 "PV_ADV_SCHEME selects the discretization for PV "//&
                 "advection. Valid values are: \n"//&
                 "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n"//&
                 "\t PV_ADV_UPWIND1  - upwind, first order", &
                 default=PV_ADV_CENTERED_STRING)
  select case (uppercase(tmpstr))
    case (PV_ADV_CENTERED_STRING)
      CS%PV_Adv_Scheme = PV_ADV_CENTERED
    case (PV_ADV_UPWIND1_STRING)
      CS%PV_Adv_Scheme = PV_ADV_UPWIND1
    case default
      call MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
                     "#DEFINE PV_ADV_SCHEME in input file is invalid.")
  end select

  CS%id_rv = register_diag_field('ocean_model', 'RV', diag%axesBL, Time, &
     'Relative Vorticity', 's-1', conversion=US%s_to_T)

  CS%id_PV = register_diag_field('ocean_model', 'PV', diag%axesBL, Time, &
     'Potential Vorticity', 'm-1 s-1', conversion=GV%m_to_H*US%s_to_T)

  CS%id_gKEu = register_diag_field('ocean_model', 'gKEu', diag%axesCuL, Time, &
     'Zonal Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_gKEv = register_diag_field('ocean_model', 'gKEv', diag%axesCvL, Time, &
     'Meridional Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxu = register_diag_field('ocean_model', 'rvxu', diag%axesCvL, Time, &
     'Meridional Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxv = register_diag_field('ocean_model', 'rvxv', diag%axesCuL, Time, &
     'Zonal Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_CAuS = register_diag_field('ocean_model', 'CAu_Stokes', diag%axesCuL, Time, &
     'Zonal Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  CS%id_CAvS = register_diag_field('ocean_model', 'CAv_Stokes', diag%axesCvL, Time, &
     'Meridional Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  !CS%id_hf_gKEu = register_diag_field('ocean_model', 'hf_gKEu', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEu_2d = register_diag_field('ocean_model', 'hf_gKEu_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_gKEv = register_diag_field('ocean_model', 'hf_gKEv', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEv_2d = register_diag_field('ocean_model', 'hf_gKEv_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_gKEu = register_diag_field('ocean_model', 'h_gKEu', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEu_2d = register_diag_field('ocean_model', 'intz_gKEu_2d', diag%axesCu1, Time, &
     'Depth-integral of Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_gKEv = register_diag_field('ocean_model', 'h_gKEv', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEv_2d = register_diag_field('ocean_model', 'intz_gKEv_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  !CS%id_hf_rvxu = register_diag_field('ocean_model', 'hf_rvxu', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxu_2d = register_diag_field('ocean_model', 'hf_rvxu_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_rvxv = register_diag_field('ocean_model', 'hf_rvxv', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxv_2d = register_diag_field('ocean_model', 'hf_rvxv_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_rvxu = register_diag_field('ocean_model', 'h_rvxu', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxu_2d = register_diag_field('ocean_model', 'intz_rvxu_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_rvxv = register_diag_field('ocean_model', 'h_rvxv', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxv_2d = register_diag_field('ocean_model', 'intz_rvxv_2d', diag%axesCu1, Time, &
     'Depth-integral of Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  ! Allocate memory needed for the diagnostics that have been enabled.
  if ((CS%id_gKEu > 0) .or. (CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. &
      (CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_gKEv > 0) .or. (CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. &
      (CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEv, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxu > 0) .or. (CS%id_hf_rvxu_2d > 0) .or. &
    ! (CS%id_hf_rvxu > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_u, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxv > 0) .or. (CS%id_hf_rvxv_2d > 0) .or. &
    ! (CS%id_hf_rvxv > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_v, IsdB, IedB, jsd, jed, nz)
  endif

  if ((CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. (CS%id_hf_rvxu > 0) .or. &
      (CS%id_hf_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_v, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. (CS%id_hf_rvxv > 0) .or. &
      (CS%id_hf_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_u, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hv, isd, ied, JsdB, JedB, nz)
  endif

end procedure CoriolisAdv_init
module procedure CoriolisAdv_end
end procedure CoriolisAdv_end
end submodule MOM_CoriolisAdv_s
