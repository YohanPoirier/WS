module GeomDisque
!
use GeomStruct
use GeomFonct
use Incident_mod
!
contains
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! - calcul coordonnees cartesienne depuis u,v
! ----------------------------------------------------------------------------
subroutine eq_dis3D(u,v,geom,x)
  implicit none
  real(rp),intent(in) :: u,v
  !f2py integer*1, dimension(1000) :: geom
  type(cylindrical_surface2),intent(in) :: geom
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp) :: r,d  
  !f2py integer*1, dimension(1000) :: N,v_translat
  type(vector) :: N,v_translat
  !f2py integer*1, dimension(1000) :: C
  type(point)  :: C
  real(rp) :: gamma,phi
!  
  r = geom%radius
  N = geom%normal
  C = geom%centre
  call assign_vector_coord(v_translat,C%coord)
  gamma = geom%gamma
  phi = geom%phi
!  
  d = sqrt(C%coord(1)*C%coord(1)+C%coord(2)*C%coord(2)+C%coord(3)*C%coord(3))
  x = (/v*cos(u),v*sin(u),C%coord(3)/)
  if (abs(gamma).gt.Epsilon) call mat_rotation(x,gamma,2,[0._rp,0._rp,0._rp],x)
  if (abs(phi).gt.Epsilon) call mat_rotation(x,phi,3,[0._rp,0._rp,0._rp],x)
!  
end subroutine eq_dis3D 
!
! ---------------------------------------------------------------------------
! Calcul coordonnées paramétrique d'un point
! ---------------------------------------------------------------------------
subroutine cart2dis(P,geom,param)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: geom
  type(cylindrical_surface2),intent(in) :: geom
  real(rp),dimension(2),intent(out) :: param
! local  
  !f2py integer*1, dimension(1000) :: Q,C
  type(point) :: Q,C
  real(rp) :: inv_r,r,x
  integer :: ierror
!  
  r = geom%radius
  C = geom%centre
  call mat_rotation(P%coord,-geom%phi,3,C%coord,Q%coord)
  call mat_rotation(Q%coord,-geom%gamma,2,C%coord,Q%coord)
  param(2) = sqrt(Q%coord(1)*Q%coord(1)+Q%coord(2)*Q%coord(2))
  inv_r = 1/param(2)
!  
  x = Q%coord(1)*inv_r
  if(abs(x+1.).lt.0.001)then
    param(1) = PI
  elseif(abs(x-1.).lt.0.001)then
    param(1) = 0._rp
  elseif(abs(x).lt.1)then
    param(1) = dacos(x)
  else
    ierror = 1
    goto 9999
  endif
  if(Q%coord(2).lt.0)then
    param(1) = -param(1)
  endif
!   
9999 continue
  if(ierror .ne.0)then
    write(*,*),'*** error : calcul coordonnee cylindrique'
  endif
!
end subroutine cart2dis
! 
! -----------------------------------------------------------------------------
! Calcul du gradient selon u
! -----------------------------------------------------------------------------
subroutine Gu_dis(geom,u,v,x)
 implicit none
 !f2py integer*1, dimension(1000) :: geom
  type(cylindrical_surface2),intent(in) :: geom
  real(rp),intent(in) :: u
  real(rp),intent(in) :: v
  real(rp),dimension(3),intent(out) :: x  
  real(rp) :: r
  !f2py integer*1, dimension(1000) :: O3
  type(point) :: O3
  r = geom%radius
  O3 = geom%centre
  x(1) = -v*sin(u)
  x(2) = v*cos(u)
  x(3) = 0.
  call mat_rotation2(x,geom%gamma,geom%phi,[0._rp,0._rp,0._rp],x)
end subroutine Gu_dis  
!
! ------------------------------------------------------------------------------
! Calcul du gradient selon v
! ------------------------------------------------------------------------------
subroutine Gv_dis(geom,u,v,x)
 implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cylindrical_surface2),intent(in) :: geom
  real(rp),intent(in) :: u
  real(rp),intent(in) :: v
  real(rp),dimension(3),intent(out) :: x  
  real(rp) :: r
  !f2py integer*1, dimension(1000) :: O3
  type(point) :: O3
  r = geom%radius
  O3 = geom%centre
  x(1) = cos(u)
  x(2) = sin(u)
  x(3) = 0.
  call mat_rotation2(x,geom%gamma,geom%phi,[0._rp,0._rp,0._rp],x)
end subroutine Gv_dis 
!
! --------------------------------------------------------------------------
! Calcul système différentiel
! --------------------------------------------------------------------------
subroutine derive_disque(P,disque,t,derive)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: disque
  type(disque2),intent(in) :: disque
  real(rp),intent(in) :: t
  real(rp),dimension(2),intent(out) :: derive
! local
  real(rp),dimension(3) :: Su,Sv
  real(rp),dimension(3) :: N2
  real(rp) :: E,F,H
  real(rp) :: SvN,SuN
  real(rp) :: denom
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
!  
  rep = disque%repere
  Su = rep%e1
  Sv = rep%e2
  
  call CGEta0(P%coord,t,N2)
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
  
end subroutine derive_disque 
! 
! ------------------------------------------------------------------------
! Traitement point proche du bord du domaine
! ------------------------------------------------------------------------
subroutine check_border_disque(param,disque,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  !f2py integer*1, dimension(1000) :: disque
  type(disque2),intent(in) :: disque
  real(rp),dimension(2),intent(in) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror
! local
  real(rp) :: r2,r2max
  real(rp) :: a,b,c,d
  real(rp) :: a1,b1,c1
  real(rp) :: det,s1,s2,s
  
  ierror = 0
  is_border = .false.
  r2max = disque%r2max
  
  r2 = dot_product(param,param)
  if(abs(r2-r2max).lt.Epsilon)then
    yout = y
    is_border = .true.
    goto 100
  elseif(r2 .gt. r2max-Epsilon)then  
    is_border = .true.
  endif
    
  if(is_border)then
    a = param(1) - y(1)
    b = y(1)
    c = param(2) - y(2)
    d = y(2)
    a1 = a*a+c*c
    b1 = 2*a*b+2*c*d
    c1 = b*b+d*d-r2max
    det = b1*b1-4.*a1*c1
    s1 = 0.5*(-b1+sqrt(det))/a1
    s2 = 0.5*(-b1-sqrt(det))/a1 
    if(s1.ge.0 .and. s1.le.1) then
      s = s1
    elseif(s2.ge.0.and.s2.le.1)then
      s = s2
    else
      ierror = 100
      goto 9999
    endif
    yout(1) = s*param(1)+(1-s)*y(1)
    yout(2) = s*param(2)+(1-s)*y(2)  
  endif
  
100 continue
9999 continue

  if(ierror.ne.0)then
    write(*,99),ierror
  endif
99 format('** error #',i3,' : check_border_disque')    
  
end subroutine check_border_disque 

subroutine assign_disque(geom,rep,C,normale,r,nface,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(inout) :: geom
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: normale
    real(rp),intent(in) :: r
    integer,intent(inout) :: nface
    integer,intent(inout) :: ierror
!   local
    integer :: k
    !f2py integer*1, dimension(1000) :: disque
    type(disque2) :: disque
!
    ierror = 0
!    
    nface = nface+1    
!
    call adapt_repere(rep,disque%repere,normale,C)
    disque%r2max = r*r
    disque%index = nface
    k = geom%ndisque + 1
    geom%disque(k) = disque
    geom%ndisque = k    
!
end subroutine assign_disque
! -------------------------------------------------------------------------------
end module GeomDisque