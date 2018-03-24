module GeomCone

use GeomStruct
use GeomFonct
use Constantes
use Incident_mod
! ***************************************************************************
contains
! ---------------------------------------------------------------------------
! Calcul du rayon 
! ---------------------------------------------------------------------------
subroutine Computation_radius(v,geom,radius) 
    implicit none
    real(rp),intent(in) :: v
    !f2py integer*1, dimension(1000) :: geom
    type(cone),intent(in) :: geom
    real(rp),intent(out) :: radius
!   local
    real(rp) :: r1,r2
    real(rp) :: denom

    r1 = geom%r1
    r2 = geom%r2
    denom = (geom%vmax - geom%vmin)
    if(abs(denom).gt.Epsilon)then
        radius = r2*(v-geom%vmin) + r1*(geom%vmax-v)
        radius = radius/denom
    else
        print*,'Error : division by 0' 
        stop
    endif

end subroutine Computation_radius
!
! ---------------------------------------------------------------------------
! Calcul coordonnées cartésienne à partir des coordonnées paramétrique (u,v)
! ---------------------------------------------------------------------------
subroutine Computation_eq_cone3D(u,v,geom,eq_cone3D)
    implicit none
    real(rp),intent(in) :: u,v
    !f2py integer*1, dimension(1000) :: geom
    type(cone),intent(in) :: geom
    real(rp),dimension(3),intent(out) :: eq_cone3D
!   local
    real(rp),dimension(3) :: x
    real(rp) :: r

    call Computation_radius(v,geom,r)
    x = [r*cos(u),r*sin(u),v]
    call loc2cart(x,geom%repere,eq_cone3D)

end subroutine Computation_eq_cone3D
!
! ----------------------------------------------------------------------------
! Calcul coordonnees parametrique d'un point
! ----------------------------------------------------------------------------
subroutine cart2cone(P,geom,param)
    implicit none
    !f2py integer*1, dimension(1000) :: P
    type(point),intent(in) :: P
    !f2py integer*1, dimension(1000) :: geom
    type(cone),intent(in) :: geom
    real(rp),dimension(2),intent(out) :: param
!   local
    !f2py integer*1, dimension(1000) :: Q
    type(point) :: Q
    real(rp) :: r, x
!
    call cart2loc(P%coord,geom%repere,Q%coord) ! Previously only Q
    call Computation_radius(Q%coord(3),geom,r)
    x = Q%coord(1)/r
!    
    if(abs(x-1).lt.Epsilon)then
      param(1) = 0.
    elseif(abs(x+1).lt.Epsilon)then
      param(1) = PI
    else
      param(1) = dacos(x)
    endif
!
    if (Q%coord(2).lt.0) then
      param(1) = -param(1)
    endif
    param(1) = mod(param(1),TWOPI)
!
    param(2) = Q%coord(3)
!    
end subroutine cart2cone
!
! -----------------------------------------------------------------------------
! Gradient selon u
! -----------------------------------------------------------------------------
subroutine Gu_cone(geom,u,v,x)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cone),intent(in) :: geom
  real(rp),intent(in) :: u, v
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp),dimension(3) :: y
  real(rp),dimension(3) :: e1,e2,e3
  real(rp) :: r  
  call Computation_radius(v,geom,r)
  e1 = geom%repere%e1
  e2 = geom%repere%e2
  e3 = geom%repere%e3
  y(1) = -r*sin(u)
  y(2) = r*cos(u)
  y(3) = 0.
  x = y(1)*e1+y(2)*e2+y(3)*e3
end subroutine Gu_cone
!
! ------------------------------------------------------------------------------
! Gradient selon v
! ------------------------------------------------------------------------------
subroutine Gv_cone(geom,u,v,x)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cone),intent(in) :: geom
  real(rp),intent(in) :: u, v
  real(rp),dimension(3) :: x
! local 
  real(rp) :: drdv
  real(rp),dimension(3) :: y
  real(rp),dimension(3) :: e1,e2,e3
  e1 = geom%repere%e1
  e2 = geom%repere%e2
  e3 = geom%repere%e3
  drdv = geom%r2 - geom%r1
  y(1) = drdv*cos(u)
  y(2) = drdv*sin(u)
  y(3) = 1.
  x = y(1)*e1 + y(2)*e2 + y(3)*e3
end subroutine Gv_cone
!
! ------------------------------------------------------------------------------
! Calcul systeme differentiel
! ------------------------------------------------------------------------------
subroutine derive_cone(P,geom,t,derive)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: geom
  type(cone),intent(in) :: geom
  real(rp),intent(in) :: t
  real(rp),dimension(2),intent(out) :: derive
! local
  real(rp),dimension(2) :: U
  real(rp),dimension(3) :: Su,Sv
  real(rp),dimension(3) :: N2
  real(rp) :: E,F,H
  real(rp) :: SvN,SuN
  real(rp) :: denom
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
!  
  rep = geom%repere
  call cart2cone(P,geom,U)
  call Gu_cone(geom,U(1),U(2),Su)
  call Gv_cone(geom,U(1),U(2),Sv)
  
  call CGEta0(P%coord,t,N2)
  !N2(1) = !N2(1)
  !N2(2) = -N2(2)
  N2(3) = -1.
!  
  E = dot_product(Su,Su)
  F = dot_product(Su,Sv)
  H = dot_product(Sv,Sv)
!
  SvN = Sv(1)*N2(1)+Sv(2)*N2(2)+Sv(3)*N2(3)
  SuN = Su(1)*N2(1)+Su(2)*N2(2)+Su(3)*N2(3)
!  
  denom = E*SvN*SvN-2*F*SvN*SuN+H*SuN*SuN
!  
  if (denom.ge.0)then
    denom = sqrt(denom)
  else
    print*,'warning : error denom neg.'
  endif
!
  derive(1) =  SvN/denom
  derive(2) = -SuN/denom   
!  
end subroutine derive_cone 
!
! ---------------------------------------------------------------------------
! Traitement des points proche du bord du domaine
! ---------------------------------------------------------------------------
subroutine check_border_cone(param,geom,y,yout,is_border,ierror)
    implicit none
    real(rp),dimension(2),intent(in) :: param
    !f2py integer*1, dimension(1000) :: geom
    type(cone),intent(in) :: geom
    real(rp),dimension(2),intent(in) :: y
    real(rp),dimension(2),intent(inout) :: yout
    logical,intent(inout) :: is_border
    integer,intent(inout) :: ierror
!   local
    logical :: inside_u,inside_v
    real(rp) :: umin,umax,vmin,vmax
    real(rp) :: denom, s
    real(rp),dimension(2) :: A,B,C,D
    real(rp),dimension(2) :: P1,P2

    ierror = 0
    is_border = .false.
    inside_u = .false.
    inside_v = .false.

    umin = 0
    if(Symmetry)then
        umax= pi
    else
        umax = 2*pi
    endif
    vmin = geom%vmin
    vmax = geom%vmax

    A = [umin,vmin]
    B = [umax,vmin]
    C = [umax,vmax]
    D = [umin,vmax]
  
    P1 = [y(1),y(2)]
    P2 = [param(1),param(2)]
!
    if(param(1).gt.umin-Epsilon .and. param(1).lt.umax+Epsilon)then
        yout(1) = param(1)
        inside_u = .true.
        if(abs(param(2)-vmax).lt.Epsilon)then
            yout(2) = vmax
            is_border = .true.
        elseif(abs(param(2)-vmin).lt.Epsilon)then
            yout(2) = vmin
            is_border = .true.
        endif
    endif

    if(param(2).gt.vmin-Epsilon .and. param(2).lt.vmax+Epsilon)then
        yout(2) = param(2)
        inside_v = .true.
        if(abs(param(1)-umax).lt.Epsilon)then
            yout(1) = umax  
            is_border = .true.
        elseif(abs(param(1)-umin).lt.Epsilon)then
            yout(1) = umin
            is_border = .true.
        endif
    endif

    if(.not.(inside_u .and. inside_v))then

    if(y(1).gt.umin-Epsilon .and. y(1).lt.umax+Epsilon)then
        if(abs(y(2)-vmax).lt.Epsilon .or. abs(y(2)-vmin).lt.Epsilon)then
            yout(1:2) = y(1:2)
            is_border = .true.
        endif
    endif
    
    if(.not.is_border)then
      if(y(2).gt.vmin-Epsilon .and. y(2).lt.vmax+Epsilon)then
        if(abs(y(1)-umax).lt.Epsilon .or. abs(y(1)-umin).lt.Epsilon)then
            yout(1:2) = y(1:2)
            is_border = .true.
        endif
      endif
    endif


    if(.not.is_border)then
      call segment_intersect_2D(P1,P2,A,B,is_border)
      if(is_border)then
        yout(2) = vmin
        denom = P1(2)-P2(2)
        if(abs(denom).gt.Epsilon)then
          s = (P1(2)-vmin)/denom
          yout(1) = s*P2(1)+(1._rp-s)*P1(1)
        else
          ierror = 100
          goto 9999
        endif
      endif
    endif 
  
    if(.not.is_border)then
      call segment_intersect_2D(P1,P2,B,C,is_border)
      if(is_border)then
        yout(1) = umax
        denom = P1(1)-P2(1)
        if(abs(denom).gt.Epsilon)then
          s = (P1(1)-umax)/denom
          yout(2) = s*P2(2)+(1._rp-s)*P1(2)
        else
          ierror = 110
          goto 9999
        endif
      endif
    endif
     
    if(.not.is_border)then  
      call segment_intersect_2D(P1,P2,C,D,is_border)
      if(is_border)then
        yout(2) = vmax
        denom = P1(2)-P2(2)
        if(abs(denom).gt.Epsilon)then
          s = (P1(2)-vmax)/denom
          yout(1) = s*P2(1)+(1._rp-s)*P1(1)
        else
          ierror = 120
          goto 9999
        endif
      endif
    endif
      
    if(.not.is_border)then
      call segment_intersect_2D(P1,P2,D,A,is_border)
      if(is_border)then
        yout(1) = umin
        denom = P1(1)-P2(1)
        if(abs(denom).gt.Epsilon)then
          s = (P1(1)-umin)/denom
          yout(2) = s*P2(2)+(1._rp-s)*P1(2)
        else
          ierror = 130
          goto 9999
        endif
      endif
    endif       

    endif
  
9999 continue
  if(ierror/=0)then 
    write(*,99),ierror
  endif
99 format('** error #',i3,' : pb. traitement point proche paroi (cone)') 
!  
end subroutine check_border_cone 

subroutine assign_cone(geom,rep,C,normale,r,angle,h,nface,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(inout) :: geom
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: normale
    real(rp),intent(in) :: r,angle,h
    integer,intent(inout) :: nface
    integer,intent(inout) :: ierror
!   local   
    integer :: k
    !f2py integer*1, dimension(1000) :: cone
    type(cone) :: cone

    ierror = 0
    
    nface = nface+1

    call adapt_repere(rep,cone%repere,normale,C)
    cone%long = h
    cone%r1 = r
    if(abs(angle-0.5*pi).lt.Epsilon)then
        ierror = 100
        goto 9999
    else
        cone%r2 = r - h*tan(angle)
    endif
    cone%vmin = 0._rp
    cone%vmax = h
    cone%index = nface
    k = geom%ncone + 1
    geom%cone(k) = cone
    geom%ncone = k    
!
9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
99 format('** error #',i3,' : cannont define cone')
!
end subroutine assign_cone
!
! *************************************************************************
end module GeomCone   