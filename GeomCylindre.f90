module GeomCylindre

use GeomStruct
use GeomFonct
use Constantes
use Incident_mod

contains
! -----------------------------------------------------------------------------
! - Calcul coordonnées cartésienne à partir des coordonnées paramétriques
! (u,v) compris entre [0,1]
! -----------------------------------------------------------------------------
subroutine Computation_eq_cyl3D(u,v,geom,eq_cyl3D)
  implicit none
  real(rp),intent(in) :: u,v
  !f2py integer*1, dimension(1000) :: geom
  type(cylindre_2),intent(in) ::geom
  real(rp),dimension(3),intent(out) :: eq_cyl3D
! Local
  real(rp),dimension(3) :: e1,e2,e3
  real(rp),dimension(3) :: x
  real(rp) :: r
!  
  r = geom%rayon
  e1 = geom%repere%e1
  e2 = geom%repere%e2
  e3 = geom%repere%e3
  x = [r*cos(u),r*sin(u),v]
  call loc2cart(x,geom%repere,eq_cyl3D)
!  
end subroutine Computation_eq_cyl3D
!
! --------------------------------------------------------------------------
! - Calcul coordonnees paramétriques d'un point
! --------------------------------------------------------------------------
subroutine cart2cyl(P,cyl,param)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: cyl
  type(cylindre_2),intent(in) :: cyl
  real(rp),dimension(2),intent(out) :: param
  real(rp) :: inv_r,x
  integer :: ierror
! local
!f2py integer*1, dimension(1000) :: Q
  type(point) :: Q
  real(rp) :: L
!  
  ierror = 0
!  
  inv_r = 1/cyl%rayon
  L = cyl%long
!  
  call cart2loc(P%coord,cyl%repere,Q%coord)
  x = Q%coord(1)*inv_r
  if(abs(x).gt.1) then
    x = sign(1._rp,x)
  endif
  if (abs(x+1.).lt.Epsilon) then
    param(1) = PI
  elseif(abs(x-1.).lt.Epsilon) then
    param(1) = 0._rp
  !elseif(abs(x).lt.1)then
  !if(abs(x-1).lt.Epsgeom)then
  !  param(1) = 0.
  !elseif(abs(x+1).lt.Epsgeom)then
  !  param(1) = PI
  else
    param(1) = dacos(x)
  endif
  !else
 !   ierror = 1
  !  goto 9999
  !endif
  if (Q%coord(2).lt.-Epsilon) then
    param(1) = -param(1)
  endif
!
  param(2) = Q%coord(3)
!  
9999 continue
  if (ierror .ne. 0)then
    write(*,*) '** error : calcul coordonne cylindrique'  
  endif
!    
end subroutine cart2cyl
!
! --------------------------------------------------------------------------
! - Gradient selon u
! ----------------------------------------------------------------------------
subroutine Gu_cyl(geom,u,v,x)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cylindre_2),intent(in) :: geom
  real(rp),intent(in) :: u
  real(rp),intent(in) :: v
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp),dimension(3) :: y
  real(rp),dimension(3) :: e1,e2,e3
  real(rp) :: r  
  r = geom%rayon
  e1 = geom%repere%e1
  e2 = geom%repere%e2
  e3 = geom%repere%e3
  y(1) = -r*sin(u)
  y(2) = r*cos(u)
  y(3) = 0.
  x = y(1)*e1+y(2)*e2+y(3)*e3
end subroutine Gu_cyl
!
! ---------------------------------------------------------------------
! Gradient selon v
!------------------------------------------------------------------------
subroutine Gv_cyl(geom,u,v,x)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cylindre_2),intent(in) :: geom
  real(rp),intent(in) :: u
  real(rp),intent(in) :: v
  real(rp),dimension(3),intent(out) :: x
! local 
  real(rp),dimension(3) :: y
  real(rp),dimension(3) :: e3
  e3 = geom%repere%e3
  y(1) = 0.
  y(2) = 0.
  y(3) = 1.
  x = e3
end subroutine Gv_cyl
!
!
! ---------------------------------------------------------------------------
! Calcul système différentiel
! ---------------------------------------------------------------------------
!
subroutine derive_cylindre(P,cyl,t,derive)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: cyl
  type(cylindre_2),intent(in) :: cyl
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
  rep = cyl%repere
  call cart2cyl(P,cyl,U)
  call Gu_cyl(cyl,U(1),U(2),Su)
  call Gv_cyl(cyl,U(1),U(2),Sv)
  
  call CGEta0(P%coord,t,N2)
  N2(3) = -1. ! ## a verifier
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
end subroutine derive_cylindre 
!
!
! ----------------------------------------------------------------------------
! Traitement point proche du bord du domaine
! ----------------------------------------------------------------------------
!
subroutine check_border_cylindre(param,cyl,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  !f2py integer*1, dimension(1000) :: cyl
  type(cylindre_2),intent(in) :: cyl
  real(rp),dimension(2),intent(in) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror
! local
  real(rp) :: vmin,vmax
!   
  ierror = 0
  is_border = .false. 
  vmax = cyl%vmax
  vmin = cyl%vmin
!  
  if(param(2).gt.(vmax-Epsilon) .or. param(2).lt.(vmin+Epsilon))then
    is_border = .true.
  endif
!  
  if(is_border)then
    if (abs(param(2)-vmax).lt.Epsilon) then 
       yout(1) = param(1)
       yout(2) = vmax
    elseif(abs(param(2)-vmin).lt.Epsilon) then
      yout(1) = param(1)
      yout(2) = vmin
    elseif(param(2).lt.vmin .and. abs(y(2)-vmin).gt.Epsilon) then
      yout(2) = vmin
      yout(1) = y(1)+(param(1)-y(1))*(param(2)-y(2))/(vmin-y(2))
    elseif(param(2).lt.vmin .and. abs(y(2)-vmin).le.Epsilon) then
      yout(2) = vmin
      yout(1) = y(1)
    elseif(param(2).gt.vmax .and. abs(vmax-y(2)).gt.Epsilon) then
      yout(2) = vmax
      yout(1) = y(1)+(param(1)-y(1))*(param(2)-y(2))/(vmax-y(2))
    elseif(param(2).gt.vmax .and. abs(vmax-y(2)).le.Epsilon) then
      yout(2) = vmax
      yout(1) = y(1)
    endif
  endif
!  
end subroutine check_border_cylindre  

subroutine assign_cylindre(geom,rep,C,normale,r,L,nface,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(inout) :: geom
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: normale
    real(rp),intent(in) :: r,L
    integer,intent(inout) :: nface
    integer,intent(inout) :: ierror
!   local
    integer :: k
    !f2py integer*1, dimension(1000) :: cylindre
    type(cylindre_2) :: cylindre
!
    ierror = 0
!    
    nface = nface+1
!
    call adapt_repere(rep,cylindre%repere,normale,C)
    cylindre%long = L
    cylindre%rayon = r
    cylindre%vmin = -0.5_rp*L
    cylindre%vmax = 0.5_rp*L
    cylindre%index = nface
    k = geom%ncylindre + 1
    geom%cylindre(k) = cylindre
    geom%ncylindre = k 
!
end subroutine assign_cylindre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local base for cylinder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LocalBase_Cyl(Normale,Puvw)
    implicit none
    real(rp),dimension(:),intent(in) :: Normale
    real(rp),dimension(:,:),intent(inout) :: Puvw
!   local
    real(rp),dimension(3),parameter :: X = [1._rp,0._rp,0._rp], Y=[0._rp,1._rp,0._rp], Z = [0._rp,0._rp,1._rp]

    if(abs(dot_product(Normale,Z)).lt.10.*Epsilon2)then
        Puvw(1:3,2) = Z
        Puvw(1:3,3) = Normale
        call Computation_vect_product(Z,Normale,Puvw(1:3,1) )
    else
        Puvw(1:3,1) = -X
        Puvw(1:3,2) = Y
        Puvw(1:3,3) = Normale
    endif

end subroutine LocalBase_Cyl

! ------------------------------------------------------------------------------
end module GeomCylindre