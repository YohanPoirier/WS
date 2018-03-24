module GeomWigley

use GeomStruct
use GeomFonct
use Constantes
use Incident_mod
use FonctionsCommunes

implicit none

contains

! -----------------------------------------------------------------------------
!   EQ_Wigleyparam3D : Cartesian coordinate from parametric coordinate
! -----------------------------------------------------------------------------
subroutine eq_Wigley3D(u,v,Wigley,val)
    
    real(rp),intent(in) :: u,v
    !f2py integer*1, dimension(1000)    :: Wigley
    type(Twigley),intent(in) :: Wigley
    real(rp),dimension(3),intent(out) :: val

    real(rp) :: yhull
    real(rp),dimension(3) :: x
    
    !x = 0.5_RP*[Wigley%L*u,Wigley%B*Hull_function(u,v),Wigley%D*v]
    call Hull_function(2._RP*u/Wigley%L,v/Wigley%D, yhull)
    x = [u,0.5_RP*Wigley%B*yhull,v]
    call loc2cart(x,Wigley%repere, val)
    !val = loc2cart(x,Wigley%repere)
!
end subroutine eq_Wigley3D

! -----------------------------------------------------------------------------
!   CART2Wigleyparam : Parametric coordinate of a given point
! -----------------------------------------------------------------------------
subroutine cart2Wigley(P,Wigley,param)
    
    !f2py integer*1, dimension(1000)    :: P
    type(point),intent(in) :: P
    !f2py integer*1, dimension(1000)    :: Wigley
    type(Twigley),intent(in) :: Wigley
    real(rp),dimension(2),intent(out) :: param
    
    !f2py integer*1, dimension(1000)    :: Q
    type(point) :: Q
    
    call cart2loc(P%coord,Wigley%repere, Q%coord)
    !Q = cart2loc(P%coord,Wigley%repere)
    !param = 2._RP*[Q%coord(1)/Wigley%L, Q%coord(3)/Wigley%D]
    param = [Q%coord(1), Q%coord(3)]
    
end subroutine cart2Wigley

! -----------------------------------------------------------------------------
!   GU_Wigleyparam : Gradient along u 
! -----------------------------------------------------------------------------
subroutine Gu_Wigley(Wigley,u,v,x)
    
    !f2py integer*1, dimension(1000)    :: Wigley
    type(Twigley),intent(in) :: Wigley
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
    
    real(rp),dimension(3) :: y
    
    !y(1) = 0.5_RP*Wigley%L
    !y(2) = 0.5_RP*Wigley%B*DHullDx(u,v)
    y(1) = 1._RP
    y(2) = Wigley%B/Wigley%L*DHullDx(2._RP*u/Wigley%L,v/Wigley%D)
    y(3) = 0._rp
    x = y(1)*Wigley%repere%e1 + y(2)*Wigley%repere%e2 + y(3)*Wigley%repere%e3
    
end subroutine Gu_Wigley

! ------------------------------------------------------------------------------
!   GV_Wigleyparam : Gradient along v
! ------------------------------------------------------------------------------
subroutine Gv_Wigley(Wigley,u,v,x)
    
    !f2py integer*1, dimension(1000)    :: Wigley
    type(Twigley),intent(in) :: Wigley
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
    
    real(rp),dimension(3) :: y
    
    y(1) = 0._RP
    !y(2) = 0.5_RP*Wigley%B*DHullDz(u,v)
    !y(3) = 0.5_RP*Wigley%D
    y(2) = 0.5_RP*Wigley%B/Wigley%D*DHullDz(2._RP*u/Wigley%L,v/Wigley%D)
    y(3) = 1._RP
    x = y(1)*Wigley%repere%e1 + y(2)*Wigley%repere%e2 + y(3)*Wigley%repere%e3
!
end subroutine Gv_Wigley

! -------------------------------------------------------------------------------
!   DERIVE_Wigleyparam : compute differential equation
! -------------------------------------------------------------------------------
subroutine derive_Wigley(P,Wigley,t,derive)
  
!f2py integer*1, dimension(1000)    :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000)    :: Wigley
  type(Twigley),intent(in) :: Wigley
  real(rp),intent(in) :: t
  real(rp),dimension(2),intent(out) :: derive
    
  real(rp),dimension(2) :: U
  real(rp),dimension(3) :: Su,Sv
  real(rp),dimension(3) :: N2
  real(rp) :: E,F,H
  real(rp) :: SvN,SuN
  real(rp) :: denom
  !f2py integer*1, dimension(1000)    :: rep
  type(repere3d) :: rep

  rep = Wigley%repere
  call cart2Wigley(P,Wigley,U)
  call Gu_Wigley(Wigley,U(1),U(2),Su)
  call Gv_Wigley(Wigley,U(1),U(2),Sv)

  call CGEta0(P%coord,t,N2)
  N2(3) = -1.

  E = dot_product(Su,Su)
  F = dot_product(Su,Sv)
  H = dot_product(Sv,Sv)

  SvN = Sv(1)*N2(1)+Sv(2)*N2(2)+Sv(3)*N2(3)
  SuN = Su(1)*N2(1)+Su(2)*N2(2)+Su(3)*N2(3)

  denom = E*SvN*SvN-2*F*SvN*SuN+H*SuN*SuN

  if (denom.ge.0)then
    denom = sqrt(denom)
  else
    print*,'warning : error denom neg.'
  endif
  
  derive(1) =  SvN/denom
  derive(2) = -SuN/denom  

end subroutine derive_Wigley

subroutine normal_Wigley(u, v, Wigley, Normal)

real(rp), intent(in) :: u,v
!f2py integer*1, dimension(1000)    :: Wigley
type(Twigley), intent(in) :: Wigley
real(rp), dimension(3), intent(out) :: Normal

real(rp) :: temp
real(rp), dimension(3) :: Du, Dv

call Gu_Wigley(Wigley,u,v,Du)
call Gv_Wigley(Wigley,u,v,Dv)
! Base avec y et z inversee
temp = Du(2)
Du(2) = Du(3)
Du(3) = temp
temp = Dv(2)
Dv(2) = Dv(3)
Dv(3) = temp
!
call Computation_vect_product(Du,Dv,Normal)
Normal = Normal/norm2(Normal)

end subroutine normal_Wigley
! -------------------------------------------------------------------------
!   DRDZ_AXI : 
! -------------------------------------------------------------------------
!function drdz_axi(v,n,Rn,Zn) result (drdv)
!    implicit none
!    real(rp),intent(in) :: v
!    real(rp),dimension(n),intent(in) :: Rn, Zn
!    integer,intent(in) :: n
!    real(rp) :: drdv
!!   local
!    integer :: j
!    integer,dimension(2) :: iP
!    real(rp) :: r1,z1,r2,z2      
!!
!    iP = find_voisin(Zn,v)
!    r1 = Rn(iP(1))
!    z1 = Zn(iP(1))
!    r2 = Rn(iP(2))
!    z2 = Zn(iP(2))    
!!
!    if(iP(1)/=iP(2))then    
!        drdv = (r2-r1)/(z2-z1)
!    elseif(iP(1)==1)then
!        drdv = (Rn(2)-r1)/(Zn(2)-z1)
!    elseif(iP(1)==n)then
!        drdv = (r1-Rn(n-1))/(z1-Zn(n-1))
!    else
!        j = iP(1)
!        drdv = (Rn(j+1)-r1)/(Zn(j+1)-z1)
!        drdv = drdv + (r1-Rn(j-1))/(z1-Zn(j-1))
!        drdv = 0.5_rp*drdv
!    endif    
!!
!end function drdz_axi
!
!
!function fz_wigley(M,wigley) result (z)
!! Parameters
!    real(rp),dimension(2),intent(in) :: M
!    type(Twigley) :: wigley
!    real(rp) :: z
!! Locals
!    real(rp) :: r
!! Begin
!    r = radius_axi_param(M(2),n,Rn,Zn)
!    if(abs(abs(M(1))-r).lt.Epsilon)then
!        z = 0._rp
!    elseif(abs(M(1)).lt.r)then
!        z = sqrt(r*r-M(1)*M(1))
!    else
!        !print*,'Warning : P_opt out of bound in fz'
!        !print*,' M : ',[M(1),M(2)],' r : ',r
!        z = 0._rp
!    endif
!! End
!end function fz_wigley
!
! -------------------------------------------------------------------------
!   MATRIX_RIEMANN : compute the rieman matrix 
! -------------------------------------------------------------------------
subroutine matrix_riemann_wigley(M,Wigley,mat,istat)
    
    real(rp),dimension(2),intent(in) :: M
    !f2py integer*1, dimension(1000)    :: Wigley
    type(TWigley), intent(in) :: Wigley
    real(rp),dimension(2,2),intent(inout) :: mat
    integer, intent(out) :: istat
    
    real(rp),dimension(3) :: Su, Sv
    real(rp) :: E, F, G, Y
    
    istat = 0
    call Hull_function(2._RP*M(1)/Wigley%L,M(2)/Wigley%D, Y)
    Y = 0.5_RP*Y
    call Gu_wigley(Wigley, M(1), M(2), Su)
    call Gv_wigley(Wigley, M(1), M(2), Sv)
    E = dot_product(Su,Su)
    F = dot_product(Su,Sv)
    G = dot_product(Sv,Sv)
    mat(1:2,1) = [E,F]
    mat(1:2,2) = [F,G]
    
end subroutine matrix_riemann_wigley

subroutine norm_wigley(A,B,Wigley,val)
    
    real(rp),dimension(2),intent(in) :: A,B
    !f2py integer*1, dimension(1000)    :: Wigley
    type(TWigley), optional :: Wigley
    real(rp),intent(out) :: val
    
    integer :: istat
    real(rp),dimension(2,2) :: matR
    
    call matrix_riemann_wigley(0.5_rp*(A+B),Wigley,matR,istat)
    call norm_riemann(matR,B-A, val)
    !val = norm_riemann(matR,B-A)   
    
end subroutine norm_wigley
!
!subroutine LocalBase_axi(Normale,Puvw,geom)
!    implicit none
!    real(rp),dimension(:),intent(in) :: Normale
!    real(rp),dimension(:,:),intent(inout) :: Puvw
!    type(axisym),intent(in) :: geom
!!   local
!    real(rp) :: denom
!    real(rp),dimension(3) :: e1,e2,e3,ev    
!!
!    e1 = geom%repere%e1
!    e2 = geom%repere%e2
!    e3 = geom%repere%e3
!!    
!    ev = vect_product(e3,Normale)
!    denom = norm2(ev)
!!
!    if(denom.gt.Epsilon2)then
!        Puvw(1:3,1) = ev/denom
!    else
!        Puvw(1:3,1) = e1
!    endif
!!
!    Puvw(1:3,3) = Normale
!    Puvw(1:3,2) = vect_product(Puvw(1:3,3),Puvw(1:3,1))    
!!
!end subroutine LocalBase_axi
!
subroutine test_prox_wigley(A,B,Q,Wigley,delta,bool)
    
    real(rp), dimension(2), intent(in) :: A, B, Q
    !f2py integer*1, dimension(1000)    :: Wigley
    type(TWigley), intent(in) :: Wigley
    real(rp), intent(in) :: delta
    logical, intent(out) :: bool
    
    integer :: flag
    real(rp) :: hA_2D, hB_2D, dist, hr2, ps2, ps1, zM, zP
    real(rp),dimension(2) :: v_AQ, v_BQ, V_AB, N_2D
    real(rp),dimension(2) :: M, P
    real(rp), dimension(3) :: V_3D
    real(rp),dimension(2,2) :: matR, matN
    
    v_AQ = Q-A
    v_BQ = Q-B

    hA_2D = norm2(v_AQ)
    hB_2D = norm2(v_BQ)

    if (hA_2D .gt. delta .and. hB_2D .gt. delta) then
        bool = .false.
    else
        M = 0.5_rp*(A+B)
        call matrix_riemann_wigley(M,Wigley,matR,flag)
        matN(1,1:2) = -matR(2,1:2)
        matN(2,1:2) = matR(1,1:2)
        V_AB = B(1:2)-A(1:2)
        dist = norm2(V_AB)
        V_AB = V_AB/dist
        N_2D(1:2) = matmul(matN,V_AB)
        call norm_riemann(matR,N_2D,hr2)
        call ps_riemann(matR,V_AB,v_AQ,ps1)
        call ps_riemann(matR,N_2D,v_AQ,ps2)
        !hr2 = norm_riemann(matR,N_2D)
        !ps1 = ps_riemann(matR,V_AB,v_AQ)
        !ps2 = ps_riemann(matR,N_2D,v_AQ)
        ps2 = ps2*hr2
        if(flag==1)then
            P = M+ps2*N_2D
            call Hull_function(2._RP*M(1)/Wigley%L,M(2)/Wigley%D, zM)
            zM = 0.5_RP*zM
            call Hull_function(2._RP*P(1)/Wigley%L,P(2)/Wigley%D, zP)
            zP = 0.5_RP*zP
            V_3D(1:2) = P-M
            V_3D(3) = zP-zM
            ps2 = norm2(V_3D)
        endif
        bool = ps1.gt.-Epsilon .and. ps1.lt.dist+Epsilon .and. abs(ps2).lt.0.1*delta
    endif
    
end subroutine test_prox_wigley    

! -------------------------------------------------------------------------
!   ASSIGN_AXISYM : create geom2
! -------------------------------------------------------------------------
subroutine WigleyGeom(rep,Dimensions,nface,nline,ierror,fgeom)
    
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep              ! Frame
    real(rp), dimension(3)              :: Dimensions       ! Dimensions
    integer,intent(inout)               :: nface,nline      ! Number of faces and lines in the geometry.
    integer,intent(inout)               :: ierror           ! Error flag
    !f2py integer*1, dimension(1000)    :: fgeom
    type(type_geom),intent(inout)       :: fgeom            ! Geometry
    
    integer                             :: iline            ! Line number
    integer                             :: j                ! Loop parameter
    !f2py integer*1, dimension(1000)    :: P1,P2,P3,P4
    type(point)                         :: P1, P2, P3, P4   ! Points
    !f2py integer*1, dimension(1000)    :: Wigley
    type(Twigley)                       :: wigley           ! Wigley structure
    
    ! This subroutine creates the geometry of a Wigley hull.   
    
    ierror = 0
    iline = 0
    
    nface = nface+1
    
    wigley%index = nface
    wigley%repere = rep
    
    wigley%L = Dimensions(1)
    wigley%B = Dimensions(2)
    wigley%D = Dimensions(3)
    
    wigley%umin = -0.5_RP*wigley%L
    wigley%umax = 0.5_RP*wigley%L
    wigley%vmin = -wigley%D
    wigley%vmax = wigley%D
    
    if(Symmetry)then
        call init_geom(fgeom,rep,[12,1,2],[1,4,4],3,ierror)
    else
        call init_geom(fgeom,rep,[12,1,2],[2,4,4],3,ierror)
    end if
    
    fgeom%nwigley = fgeom%nwigley+1
    fgeom%wigley(fgeom%nwigley) = wigley
    
    if(not(Symmetry))then
        nface = nface+1
        wigley%index = nface
        fgeom%nwigley = fgeom%nwigley+1
        wigley%repere%e1 = -wigley%repere%e1
        wigley%repere%e2 = -wigley%repere%e2
        fgeom%wigley(fgeom%nwigley) = wigley
    end if
    
    P1 = -0.5_RP*wigley%L*rep%e1 + wigley%D*rep%e3
    P2 = P1%coord - 2._rp*wigley%D*rep%e3
    P3 = P2%coord + wigley%L*rep%e1
    P4 = P3%coord + 2._rp*wigley%D*rep%e3
    
    P1%nface = 2
    P1%face(1:2) = [fgeom%wigley(1)%index,-1]
    P2%nface = 2
    P2%face(1:2) = [fgeom%wigley(1)%index,-1]
    P3%nface = 2
    P3%face(1:2) = [fgeom%wigley(1)%index,-1]
    P4%nface = 2
    P4%face(1:2) = [fgeom%wigley(1)%index,-1]
    
    fgeom%npoint = 4
    fgeom%point(1:4) = [P1,P2,P3,P4]
    
    fgeom%narete = 4
    call arete_from_point(P1,P2,fgeom%arete(1))
    call arete_from_point(P2,P3,fgeom%arete(2))
    call arete_from_point(P3,P4,fgeom%arete(3))
    call arete_from_point(P4,P1,fgeom%arete(4))
    
    iline = nline
    do j=1,4
        fgeom%arete(j)%index = iline+j
    end do
    iline = iline+4
    
    fgeom%arete(1)%face = [fgeom%wigley(1)%index,-1]
    fgeom%arete(2)%face = [fgeom%wigley(1)%index,-1]
    fgeom%arete(3)%face = [fgeom%wigley(1)%index,-1]

    fgeom%point(1)%nedge = 2
    fgeom%point(1)%edge(1:2) = [1,4] + nline
    fgeom%point(2)%nedge = 2
    fgeom%point(2)%edge(1:2) = [1,2] + nline
    fgeom%point(3)%nedge = 2
    fgeom%point(3)%edge(1:2) = [2,3] + nline
    fgeom%point(4)%nedge = 2
    fgeom%point(4)%edge(1:2) = [3,4] + nline
    
    nline = iline
    
    9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
    99 format('** error #',i3,' : cannot generate geometrie axisym')
    
end subroutine WigleyGeom

end module GeomWigley