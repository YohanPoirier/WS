module GeomPlan
!
use GeomStruct
use GeomFonct
use Incident_mod
!
contains
! -----------------------------------------------------------------------------
!
! -----------------------------------------------------------------------------
! Calcul systeme differentiel
! -----------------------------------------------------------------------------
subroutine derive_plan(P,plan,t,derive)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: plan
  type(Gplan),intent(in) :: plan
  real(rp),intent(in) :: t
  real(rp),dimension(2) :: derive
! local
  real(rp),dimension(3) :: Su,Sv
  real(rp),dimension(3) :: N2
  real(rp) :: E,F,H
  real(rp) :: SvN,SuN
  real(rp) :: denom
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
!  
  rep = plan%repere
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
  
end subroutine derive_plan
!
! -----------------------------------------------------------------------
! Traitement point proche du bord du domaine
! -----------------------------------------------------------------------
subroutine check_border_plan(param,plan,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  !f2py integer*1, dimension(1000) :: plan
  type(Gplan),intent(in) :: plan
  real(rp),dimension(2) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror
! local
  real(rp) :: umin,umax
  real(rp) :: vmin,vmax
  real(rp) :: denom,s
  real(rp),dimension(2) :: A,B,C,D
  real(rp),dimension(2) :: P1,P2
  
  ierror = 0
  is_border = .false.
  
  umin = plan%umin
  umax = plan%umax
  vmin = plan%vmin
  vmax = plan%vmax
  
  A = [umin,vmin]
  B = [umax,vmin]
  C = [umax,vmax]
  D = [umin,vmax]
  
  P1 = [y(1),y(2)]
  P2 = [param(1),param(2)]
  
  if(param(1).lt.umin .or. param(1).gt.umax .or. param(2).lt.vmin .or. param(2).gt.vmax)then
    if((abs(y(1)-umin).lt.Epsilon .or. abs(y(1)-umax).lt.Epsilon).or.&
    &  (abs(y(2)-vmin).lt.Epsilon .or. abs(y(2)-vmax).lt.Epsilon))then
      yout = y
      is_border = .true.
    endif
  endif
  
  if(.not.is_border)then
    if((abs(param(1)-umin).lt.Epsilon .or. abs(param(1)-umax).lt.Epsilon).or.&
    &  (abs(param(2)-vmin).lt.Epsilon .or. abs(param(2)-vmax).lt.Epsilon))then
      yout = param
      is_border = .true.
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
  
9999 continue
  if(ierror/=0)then 
    write(*,99),ierror
  endif
99 format('** error #',i3,' : pb. traitement point proche paroi (plan)')    

end subroutine check_border_plan
! ----------------------------------------------------------------------
end module GeomPlan
