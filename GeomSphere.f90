module GeomSphere

use GeomStruct
use GeomFonct
use Constantes
use Incident_mod

contains

subroutine Computation_eq_sph3D(u,v,geom,eq_sph3D)
    implicit none
    real(rp),intent(in) :: u,v
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),dimension(3),intent(out) :: eq_sph3D
    
    real(rp),dimension(3) :: x
    real(rp) :: r
    
    r = geom%radius
    x = [r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v)]
    call loc2cart(x,geom%repere,eq_sph3D)
    
end subroutine Computation_eq_sph3D

subroutine cart2sph_f(P,geom,param)
    implicit none
    !f2py integer*1, dimension(1000) :: P
    type(point),intent(in) :: P
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),dimension(2),intent(out) :: param
!   local
    integer :: ierror
    real(rp) :: zr
    !f2py integer*1, dimension(1000) :: Q
    type(point) :: Q
!
    ierror = 0
!
    call cart2loc(P%coord,geom%repere,Q%coord)
!
    if(Q%coord(1).gt.Epsilon)then
        param(1) = atan(Q%coord(2)/Q%coord(1))
    elseif(Q%coord(1) .lt. -Epsilon)then
        param(1) = atan(Q%coord(2)/Q%coord(1))+PI
    else
        param(1) = 0.5_rp*PI*sign(1._rp,Q%coord(2))
    endif
    param(1) = mod(param(1),TWOPI)
!
    zr = Q%coord(3)/geom%radius
    if(abs(zr).gt.1)then
        zr = sign(1._rp,zr)
    endif
    param(2) = asin(zr)
!
9999 continue
  if (ierror .ne. 0)then
    write(*,*) '** error : calcul coordonne spherique'  
  endif
!
end subroutine cart2sph_f

subroutine proj2cart(M,geom,x)
    implicit none
    real(rp),dimension(2),intent(in) :: M
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),dimension(3),intent(out) :: x
!   local
    real(rp) :: r,l2
    r = geom%radius
    l2 = M(1)*M(1)+M(2)*M(2)
    x(1:2) = M(1:2)
    x(3) = sqrt(r*r - l2)
    call loc2cart(x,geom%repere,x)
end subroutine proj2cart

subroutine Gu_sph(geom,u,v,x)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
!   local
    real(rp) :: r
    real(rp),dimension(3) :: e1,e2,e3
    r = geom%radius
    e1 = geom%repere%e1
    e2 = geom%repere%e2
    e3 = geom%repere%e3
    x(1) = -r*cos(v)*sin(u)
    x(2) = r*cos(v)*cos(u)
    x(3) = 0._rp
    x = x(1)*e1 + x(2)*e2 + x(3)*e3
end subroutine Gu_sph

subroutine Gv_sph(geom,u,v,x)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
!   local
    real(rp) :: r
    real(rp),dimension(3) :: e1,e2,e3
    r = geom%radius
    e1 = geom%repere%e1
    e2 = geom%repere%e2
    e3 = geom%repere%e3
    x(1) = -r*sin(v)*cos(u)
    x(2) = -r*sin(v)*sin(u)
    x(3) = r*cos(v)
    x = x(1)*e1 + x(2)*e2 + x(3)*e3
end subroutine Gv_sph

subroutine derive_sphere(P,geom,t,derive)
    implicit none
    !f2py integer*1, dimension(1000) :: P
    type(point),intent(in) :: P
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),intent(in) :: t
    real(rp),dimension(2),intent(out) :: derive
!   local
    real(rp),dimension(2) :: U
    real(rp),dimension(3) :: Su,Sv
    real(rp),dimension(3) :: N2
    real(rp) :: E,F,H
    real(rp) :: SvN,SuN
    real(rp) :: denom

    call cart2sph_f(P,geom,U)
    call Gu_sph(geom,U(1),U(2),Su)
    call Gv_sph(geom,U(1),U(2),Sv)

    call CGEta0(P%coord,t,N2)
    N2(3) = -1._rp ! ##a verifier

    E = dot_product(Su,Su)
    F = dot_product(Su,Sv)
    H = dot_product(Sv,Sv)

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
end subroutine derive_sphere

subroutine check_border_sphere(param,geom,y,yout,is_border,ierror)
    implicit none
    real(rp),dimension(2),intent(in) :: param
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
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

    umin = geom%umin
    umax = geom%umax
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
99 format('** error #',i3,' : pb. traitement point proche paroi (sphere)') 

end subroutine check_border_sphere

subroutine assign_sphere(geom,rep,C,normale,r,phi1,phi2,nface,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(inout) :: geom
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: normale
    real(rp),intent(in) :: r,phi1,phi2
    integer,intent(inout) :: nface
    integer,intent(inout) :: ierror
!   local
    integer :: k
    !f2py integer*1, dimension(1000) :: sphere
    type(sphere) :: sphere
!
    ierror = 0
!
    nface = nface+1
!
    call adapt_repere(rep,sphere%repere,normale,C)
    sphere%index = nface
    sphere%radius = r
    sphere%umin = 0._rp
    if(Symmetry)then
        sphere%umax = pi
    else
        sphere%umax = 2*pi
    endif
    sphere%vmin = phi1
    sphere%vmax = phi2
    k = geom%nsphere+1
    geom%sphere(k) = sphere
    geom%nsphere = k
!
end subroutine assign_sphere

end module GeomSphere