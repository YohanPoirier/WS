module GeomAxiSym

use GeomStruct
use GeomFonct
use Constantes
use Incident_mod
use FonctionsCommunes

implicit none

contains

! ------------------------------------------------------------------------
!   READ_AXISYM : Read the input file
! ------------------------------------------------------------------------
subroutine read_axisym(filename,geom)
    
    character(len=50),intent(in)        :: filename ! File name.
    !f2py integer*1, dimension(1000)    :: geom
    type(axisym),intent(inout)          :: geom     ! Axisym geometry.
    
    integer                             :: j,n      ! Loop parameters.
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d)                      :: rep      ! Frame.
    real(rp),dimension(3)               :: C,vect   ! Vectors.
    integer,dimension(:),allocatable    :: order    ! ?
    
    ! This subroutine reads the axisym input file.
    
    print*,"read_axisym info:"
    print*,"First line: origine of the local frame."
    print*,"Second line: vertical unit vector coordinates."
    print*,"Third line: number of points."
    print*,"Then, coordinates of the points (radius,height)"
    
    rep%origine = [0._rp,0._rp,0._rp]
    rep%e1 = [1._rp,0._rp,0._rp]
    rep%e2 = [0._rp,1._rp,0._rp]
    rep%e3 = [0._rp,0._rp,1._rp]    
    
    ! Opening
    open(unit=ioaxisym,file=filename)
    
    ! Center
    read(ioaxisym,*) C
    
    ! Direction vector
    read(ioaxisym,*) vect    
    call adapt_repere(rep,geom%repere,vect,C)
    
    read(ioaxisym,*) n
    geom%npoint = n

    allocate(geom%P(n,2))
    allocate(order(n))
    
    ! Coordinates
    do j=1,n
        read(ioaxisym,*) geom%P(j,1:2)
    enddo
    
    call quick_sort(geom%P(:,2),order)
    call sort_order(geom%P(:,1),order)    
    
    ! Closing
    close(ioaxisym)

    deallocate(order)

end subroutine read_axisym

! ----------------------------------------------------------------------------
!   RADIUS_AXI : Compute radius from geom def.
! ----------------------------------------------------------------------------
subroutine radius_axi(geom,v,val)
    
    !f2py integer*1, dimension(1000) :: geom
    type(axisym),intent(in) :: geom
    real(rp),intent(in) :: v
    real(rp),intent(out) :: val
    
    integer,dimension(2) :: iP
    real(rp) :: r1,z1,r2,z2
    
    ! This subroutine computes the radius at the height v in using a linear interpolation between two neighbours.
    
    ! Neighbours.
    call find_voisin(geom%P(:,2),v,iP)
    r1 = geom%P(iP(1),1)
    z1 = geom%P(iP(1),2)
    r2 = geom%P(iP(2),1)
    z2 = geom%P(iP(2),2)
    
    ! Computation of the radius at v.
    if(iP(1)/=iP(2))then
        val = (r2*(v-z1)+r1*(z2-v))/(z2-z1) ! Linear interpolation between (r1,z1) and (r2,z2) to find r.
    else
        val = r1
    endif
    
end subroutine radius_axi

subroutine radius_axi_param(v,n,Rn,Zn,val)
    implicit none
    real(rp),intent(in) :: v
    integer,intent(in) :: n
    real(rp),dimension(:),intent(in) :: Rn, Zn
    real(rp),intent(out) :: val
!   local
    integer,dimension(2) :: iP
    real(rp) :: r1,z1,r2,z2
!
    call find_voisin(Zn(:),v,iP)
    r1 = Rn(iP(1))
    z1 = Zn(iP(1))
    r2 = Rn(iP(2))
    z2 = Zn(iP(2))
!
    if(iP(1)/=iP(2))then
        val = (r2*(v-z1)+r1*(z2-v))/(z2-z1)
    else
        val = r1
    endif
!
end subroutine radius_axi_param

! -----------------------------------------------------------------------------
!   EQ_AXI3D : Cartesian coordinate from parametric coordinate
! -----------------------------------------------------------------------------
subroutine eq_axi3D(u,v,geom,val)
    
    real(rp),intent(in)                 :: u,v  ! Parameter coordinates
    !f2py integer*1, dimension(1000)    :: geom
    type(axisym),intent(in)             :: geom ! Geometry
    real(rp),dimension(3),intent(out)   :: val  ! Cartesian coordinates
    
    real(rp),dimension(3)               :: x    ! Cartesian coordinates
    real(rp)                            :: r    ! Radius
   
    ! This subroutine computes the cartesian coordinates from the parametric ones.
    
    call radius_axi(geom,v,r)
    x = [r*cos(u),r*sin(u),v]
    call loc2cart(x,geom%repere,val)
    
end subroutine eq_axi3D

! -----------------------------------------------------------------------------
!   CART2AXI : Parametric coordinate of a given point
! -----------------------------------------------------------------------------
subroutine cart2axi(P,geom,param,ierr)
    
    !f2py integer*1, dimension(1000)    :: P
    type(point),intent(in)              :: P    ! Point.
    !f2py integer*1, dimension(1000)    :: geom 
    type(axisym),intent(in)             :: geom ! Geometry.
    real(rp),dimension(2),intent(out)   :: param! Axisym coordinates (theta, z).
    integer,optional                    :: ierr ! Error flag.
    
    !f2py integer*1, dimension(1000)    :: Q
    type(point)                         :: Q    ! Point.
    real(rp)                            :: r,x  ! Radius and x-coordinate.
    
    integer,dimension(2) :: iP
    real(rp) :: r1,z1,r2,z2
    
    ! This subroutine transforms the Cartesian coordinates into the axisym coordinates.
    
    
    call cart2loc(P%coord,geom%repere,Q%coord) ! Previously only Q
    call radius_axi(geom,Q%coord(3),r)
    x = Q%coord(1)/r
    
    ! abs(x) can be higher than 1 because of the linear interpolation in the subroutine radius_axi.
    if(x.gt.1._RP)then
        x = 1._RP
    else if(x.lt.(-1._RP))then
        x = -1._RP
    end if
    
    if(abs(x-1).lt.Epsilon)then ! x = 1
        param(1) = 0._rp
    elseif(abs(x+1).lt.Epsilon)then ! x = -1
        param(1) = PI
    else ! -1 < x < 1
        param(1) = dacos(x)
    endif
    
    if(Q%coord(2).lt.0)then
        param(1) = -param(1)
    endif
    param(1) = mod(param(1),TWOPI)
    
    param(2) = Q%coord(3)
    
    9999 continue    
    
end subroutine cart2axi

! -----------------------------------------------------------------------------
!   GU_AXI : Gradient along u 
! -----------------------------------------------------------------------------
subroutine Gu_axi(geom,u,v,x)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(axisym),intent(in) :: geom
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
!   local
    real(rp),dimension(3) :: y
    real(rp) :: r
!
    call radius_axi(geom,v,r) 
    y(1) = -r*sin(u)
    y(2) = r*cos(u)
    y(3) = 0._rp
    x = y(1)*geom%repere%e1 + y(2)*geom%repere%e2 + y(3)*geom%repere%e3
!
end subroutine Gu_axi

! ------------------------------------------------------------------------------
!   GV_AXI : Gradient along v
! ------------------------------------------------------------------------------
subroutine Gv_axi(geom,u,v,x)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(axisym),intent(in) :: geom
    real(rp),intent(in) :: u,v
    real(rp),dimension(3),intent(out) :: x
!   local
    integer :: j,n
    integer,dimension(2) :: iP
    real(rp),dimension(3) :: y
    real(rp) :: r1,r2,z1,z2,drdv
!
    call find_voisin(geom%P(:,2),v,iP)
    r1 = geom%P(iP(1),1)
    z1 = geom%P(iP(1),2)
    r2 = geom%P(iP(2),1)
    z2 = geom%P(iP(2),2)
    n = geom%npoint
!
    if(iP(1)/=iP(2))then    
        drdv = (r2-r1)/(z2-z1)
    elseif(iP(1)==1)then
        drdv = (geom%P(2,1)-r1)/(geom%P(2,2)-z1)
    elseif(iP(1)==n)then
        drdv = (r1-geom%P(n-1,1))/(z1-geom%P(n-1,2))
    else
        j = iP(1)
        drdv = (geom%P(j+1,1)-r1)/(geom%P(j+1,2)-z1)
        drdv = drdv + (r1-geom%P(j-1,1))/(z1-geom%P(j-1,2))
        drdv = 0.5_rp*drdv
    endif
!    
    y(1) = drdv*cos(u)
    y(2) = drdv*sin(u)
    y(3) = 1._rp
    x = y(1)*geom%repere%e1 + y(2)*geom%repere%e2 + y(3)*geom%repere%e3
!
end subroutine Gv_axi

! -------------------------------------------------------------------------------
!   DERIVE_AXI : compute differential equation
! -------------------------------------------------------------------------------
subroutine derive_axi(P,geom,t,derive)
  implicit none
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(in) :: P
  !f2py integer*1, dimension(1000) :: geom
  type(axisym),intent(in) :: geom
  real(rp),intent(in) :: t
  real(rp),dimension(2) :: derive
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
  call cart2axi(P,geom,U)
  call Gu_axi(geom,U(1),U(2),Su)
  call Gv_axi(geom,U(1),U(2),Sv)
!  
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
! 
end subroutine derive_axi

! -------------------------------------------------------------------------
!   DRDZ_AXI : 
! -------------------------------------------------------------------------
subroutine drdz_axi(v,n,Rn,Zn,drdv)
    implicit none
    real(rp),intent(in) :: v
    real(rp),dimension(n),intent(in) :: Rn, Zn
    integer,intent(in) :: n
    real(rp),intent(out) :: drdv
!   local
    integer :: j
    integer,dimension(2) :: iP
    real(rp) :: r1,z1,r2,z2      
!
    call find_voisin(Zn,v,iP)
    r1 = Rn(iP(1))
    z1 = Zn(iP(1))
    r2 = Rn(iP(2))
    z2 = Zn(iP(2))    
!
    if(iP(1)/=iP(2))then    
        drdv = (r2-r1)/(z2-z1)
    elseif(iP(1)==1)then
        drdv = (Rn(2)-r1)/(Zn(2)-z1)
    elseif(iP(1)==n)then
        drdv = (r1-Rn(n-1))/(z1-Zn(n-1))
    else
        j = iP(1)
        drdv = (Rn(j+1)-r1)/(Zn(j+1)-z1)
        drdv = drdv + (r1-Rn(j-1))/(z1-Zn(j-1))
        drdv = 0.5_rp*drdv
    endif    
!
end subroutine drdz_axi

! -------------------------------------------------------------------------
!   BASELOCAL_AXISYM : compute the local base assosciated with the geom.
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
!   NORMAL_AXISYM : compute the normale
! -------------------------------------------------------------------------
subroutine normal_axisym(P,n,Rn,Zn,vect)
    implicit none
    real(rp),dimension(3),intent(in) :: P
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: Rn, Zn
    real(rp),dimension(3),intent(out) :: vect
!   local
    real(rp) :: r, drdz, a, d
    real(rp),dimension(3) :: Gu, Gv

    r = sqrt(P(1)*P(1)+P(3)*P(3))

    if(abs(r).lt.Epsilon)then
        vect = [0._rp,-1._rp,0._rp]
    else
        call drdz_axi(P(2),n,Rn,Zn,drdz)
        a = drdz/r

        Gu = [-P(3),0._rp,P(1)]
        Gv = [a*P(1),1._rp,a*P(3)]

        call Computation_vect_product(Gv,Gu,vect)
        d = norm2(vect)
        if(d.gt.Epsilon)then
            vect = vect / norm2(vect)
        else
            vect = [-999.,-999.,-999.]
        endif
    endif

end subroutine normal_axisym
!
subroutine fz_axisym(M,n,Rn,Zn,z)
    implicit none
    real(rp),dimension(2),intent(in) :: M
    integer,intent(in) :: n
    real(rp),dimension(n) :: Rn, Zn
    real(rp),intent(out) :: z
    real(rp) :: r
    call radius_axi_param(M(2),n,Rn,Zn,r)
    if(abs(abs(M(1))-r).lt.Epsilon)then
        z = 0._rp
    elseif(abs(M(1)).lt.r)then
        z = sqrt(r*r-M(1)*M(1))
    else
        !print*,'Warning : P_opt out of bound in fz'
        !print*,' M : ',[M(1),M(2)],' r : ',r
        z = 0._rp
    endif
end subroutine fz_axisym

! -------------------------------------------------------------------------
!   MATRIX_RIEMANN : compute the rieman matrix 
! -------------------------------------------------------------------------
subroutine matrix_riemann_axisym(M,n,Rn,Zn,mat,flag)
    implicit none
    real(rp),dimension(2),intent(in) :: M
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: Rn, Zn
    real(rp),dimension(2,2),intent(inout) :: mat
    integer,intent(inout) :: flag
    
    real(rp),dimension(3) :: Su, Sv
    real(rp) :: E, F, G
    real(rp) :: r, drdz, denom
    
    flag = 0
    
    call radius_axi_param(M(2),n,Rn,Zn,r)
    call drdz_axi(M(2),n,Rn,Zn,drdz)
    
    if(abs(abs(M(1)*M(1))-r*r).lt.Epsilon)then
        mat(1:2,1) = [1._rp,0._rp]
        mat(1:2,2) = [0._rp,1._rp]
        flag = 1
    elseif(abs(M(1)*M(1)).lt.r*r)then
        denom = sqrt(r*r - M(1)*M(1))
        Su = [1._rp, 0._rp, -M(1)/denom]
        Sv = [0._rp, 1._rp, r*drdz/denom]
        E = dot_product(Su,Su)
        F = dot_product(Su,Sv)
        G = dot_product(Sv,Sv)
        mat(1:2,1) = [E,F]
        mat(1:2,2) = [F,G]
    else
        mat(1:2,1) = [1._rp,0._rp]
        mat(1:2,2) = [0._rp,1._rp]
        flag = 2
    endif
    
    if(idebug>0 .and. flag>0)then
        write(*,10) flag
    endif
10 format('matrix riemann flag = ',i3)
     
end subroutine matrix_riemann_axisym

subroutine norm_axisym(A,B,n,Rn,Zn,val)
    implicit none
    real(rp),dimension(2),intent(in) :: A, B
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: Rn, Zn
    real(rp),intent(out) :: val
!   local
    integer :: flag_A, flag_B
    real(rp),dimension(2,2) :: matR_A, matR_B
    real(rp),dimension(2) :: V_AB
    real(rp),dimension(3) :: V_3D
    real(rp) :: hA, hB, zA, zB

    !if(idebug>0) print*,'norm2 axi...'
!
    call matrix_riemann_axisym(A,n,Rn,Zn,matR_A,flag_A)
    call matrix_riemann_axisym(B,n,Rn,Zn,matR_B,flag_B)
!
    V_AB = B-A
!   
    if(flag_A==0 .and. flag_B == 0)then
        call norm_riemann(matR_A,V_AB,hA)
        call norm_riemann(matR_B,V_AB,hB)
        val = 0.5_rp*(hA+hB)
    else
        call fz_axisym(A,n,Rn,Zn,zA)
        call fz_axisym(B,n,Rn,Zn,zB)
        V_3D(1:2) = V_AB
        V_3D(3) = zB-zA
        val = norm2(V_3D)
    endif
!
end subroutine norm_axisym

subroutine norm_axisym2(A,B,n,Rn,Zn,val)
    implicit none
    real(rp),dimension(2),intent(in) :: A,B
    integer,intent(in) :: n
    real(rp),dimension(n),intent(in) :: Rn,Zn
    real(rp),intent(out) :: val
!   local
    integer :: flag
    real(rp),dimension(2,2) :: matR
    real(rp),dimension(2) :: M,V_AB

    !if(idebug>0) print*,'norm2 ...'
!
    M = 0.5_rp*(A+B)
    V_AB = B-A
!
    call matrix_riemann_axisym(M,n,Rn,Zn,matR,flag)
    call norm_riemann(matR,V_AB,val)   
!
end subroutine norm_axisym2

subroutine LocalBase_axi(Normale,Puvw,axi)
    implicit none
    real(rp),dimension(3),intent(in) :: Normale
    real(rp),dimension(3,3),intent(inout) :: Puvw
    !f2py integer*1, dimension(1000) :: axi
    type(axisym) :: axi
    
!   local
    real(rp) :: denom
    real(rp),dimension(3) :: e1,e2,e3,ev    
!
    e1 = axi%repere%e1
    e2 = axi%repere%e2
    e3 = axi%repere%e3
!    
    call Computation_vect_product(e3,Normale,ev)
    denom = norm2(ev)
!
    if(denom.gt.Epsilon2)then
        Puvw(1:3,1) = ev/denom
    else
        Puvw(1:3,1) = e1
    endif
!
    Puvw(1:3,3) = Normale
    call Computation_vect_product(Puvw(1:3,3),Puvw(1:3,1),Puvw(1:3,2))    
!
end subroutine LocalBase_axi

subroutine test_prox_axisym(A,B,Q,param,delta,bool)
    implicit none
    real(rp),dimension(2),intent(in) :: A, B, Q
    real(rp),dimension(*),intent(in) :: param
    real(rp),intent(in) :: delta
    logical,intent(out) :: bool
!   local
    integer :: n, flag
    real(rp) :: hA_2D, hB_2D, dist, hr2, ps2, ps1, zM, zP
    real(rp),dimension(2) :: v_AQ, v_BQ, V_AB, N_2D
    real(rp),dimension(3) :: V_3D
    real(rp),dimension(2) :: M, P
    real(rp),dimension(2,2) :: matR, matN

    !if(idebug>0) print*,'test prox ...'

    n = int(param(1))

    v_AQ = Q-A
    v_BQ = Q-B

    hA_2D = norm2(v_AQ)
    hB_2D = norm2(v_BQ)

    if(hA_2D .gt. delta .and. hB_2D .gt. delta)then
      bool = .false.
    else

      M = 0.5_rp*(A+B)
      call matrix_riemann_axisym(M,n,param(2:n+1),param(n+2:2*n+1),matR,flag)
      matN(1,1:2) = -matR(2,1:2)
      matN(2,1:2) = matR(1,1:2)

      V_AB = B(1:2)-A(1:2)
      dist = norm2(V_AB)
      V_AB = V_AB/dist
      N_2D(1:2) = matmul(matN,V_AB)

      !hr1 = norm_riemann(matR,V_AB)
      call norm_riemann(matR,N_2D,hr2)

      call ps_riemann(matR,V_AB,v_AQ,ps1)
      call ps_riemann(matR,N_2D,v_AQ,ps2)
      !ps1 = ps1*hr1
      ps2 = ps2*hr2

      if(flag==1)then
        P = M+ps2*N_2D
        call fz_axisym(M,n,param(2:n+1),param(n+2:2*n+1),zM)
        call fz_axisym(P,n,param(2:n+1),param(n+2:2*n+1),zP)
        V_3D(1:2) = P-M
        V_3D(3) = zP-zM
        ps2 = norm2(V_3D)
      endif

      bool = ps1.gt.-Epsilon .and. ps1.lt.dist+Epsilon .and. abs(ps2).lt.0.1*delta
    
    endif

end subroutine test_prox_axisym    
      
! -------------------------------------------------------------------------
!   ASSIGN_AXISYM : create geom
! -------------------------------------------------------------------------
subroutine wavestar1(filename,rep,nface,nline,ierror,fgeom)
    implicit none
    character(len=50),intent(in) :: filename
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    integer,intent(inout) :: nface,nline
    integer,intent(inout) :: ierror
    !f2py integer*1, dimension(1000) :: fgeom
    type(type_geom),intent(out) :: fgeom
!   local
    integer :: n,np,j,np0,narete,i1,i2,nedge,nf
    !f2py integer*1, dimension(1000) :: axisym
    type(axisym) :: axisym
    !f2py integer*1, dimension(1000) :: rep1
    type(repere3d) :: rep1
    !f2py integer*1, dimension(1000) :: P
    type(point) :: P
    integer,dimension(2) :: face

    ierror = 0

    nface = nface+1

    call read_axisym(filename,axisym)
    axisym%index = nface

    axisym%vmin = minval(axisym%P(:,2))
    axisym%vmax = maxval(axisym%P(:,2))

    axisym%umin = 0._rp
    if(Symmetry)then
        axisym%umax = PI
    else
        axisym%umax = TWOPI
    endif

    n = axisym%npoint
    np0 = 0
    do j=1,n
        if(abs(axisym%P(j,2)).lt.Epsilon)then
            np0 = np0+1
        endif
    enddo

    call init_geom(fgeom,rep,[10,1,2],[1,2*n-np0,2*(n-1)],3,ierror)

    fgeom%naxisym = fgeom%naxisym+1
    fgeom%axisym(fgeom%naxisym) = axisym

    rep1 = axisym%repere

    np = fgeom%npoint
    narete = fgeom%narete

    face = [axisym%index,-1]
    if(Symmetry)then
        nf = 2
    else
        nf = 1
    endif

    P%nface = nf
    P%face(1:nf) = face(1:nf)

    P%nedge = 0
    P%edge(:) = 0

    P%coord = rep1%origine + axisym%P(1,1)*rep1%e1 + axisym%P(1,2)*rep1%e3
    np = np+1
    fgeom%point(np) = P
    i1 = np ; i2 = np
    if(abs(P%coord(1)).gt.Epsilon)then
        np = np+1
        P%coord(1) = -P%coord(1)
        fgeom%point(np) = P
        i2 = np
    endif

    do j=2,n
!
        P%coord = rep1%origine + axisym%P(j,1)*rep1%e1 + axisym%P(j,2)*rep1%e3
        np=np+1
        fgeom%point(np) = P
!
        narete = narete+1
        call arete_from_point(fgeom%point(i1),fgeom%point(np),fgeom%arete(narete))
        nline = nline+1
        fgeom%arete(narete)%index = nline
        fgeom%arete(narete)%nface = nf
        fgeom%arete(narete)%face(1:nf) = face(1:nf)        
!
        nedge = fgeom%point(i1)%nedge
        nedge = nedge+1
        fgeom%point(i1)%edge(nedge) = nline
        fgeom%point(i1)%nedge = nedge
        nedge = fgeom%point(np)%nedge
        nedge = nedge+1
        fgeom%point(np)%edge(nedge) = nline
        fgeom%point(np)%nedge = nedge
!
        i1 = np
!
        if(abs(P%coord(1)).gt.Epsilon)then
            np=np+1
            P%coord(1) = -P%coord(1)
            fgeom%point(np) = P
        endif
!
        narete = narete+1
        call arete_from_point(fgeom%point(np),fgeom%point(i2),fgeom%arete(narete))
        nline = nline+1
        fgeom%arete(narete)%index = nline
        fgeom%arete(narete)%nface = nf
        fgeom%arete(narete)%face(1:nf) = face(1:nf)
!          
        nedge = fgeom%point(i2)%nedge
        nedge = nedge+1
        fgeom%point(i2)%edge(nedge) = nline
        fgeom%point(i2)%nedge = nedge
        nedge = fgeom%point(np)%nedge
        nedge = nedge+1
        fgeom%point(np)%edge(nedge) = nline
        fgeom%point(np)%nedge = nedge            
 !         
        i2 = np

    enddo
!
    fgeom%npoint = np
    fgeom%narete = narete
!
end subroutine wavestar1


! -------------------------------------------------------------------------
!   ASSIGN_AXISYM : create geom2
! -------------------------------------------------------------------------
subroutine wavestar2(filename,rep,nface,nline,ierror,fgeom)
    
    character(len=50),intent(in)        :: filename
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep
    integer,intent(inout)               :: nface,nline
    integer,intent(inout)               :: ierror
    !f2py integer*1, dimension(1000)    :: fgeom
    type(type_geom),intent(out)         :: fgeom
    
    !f2py integer*1, dimension(1000)    :: rep1
    type(repere3d)                      :: rep1
    integer                             :: index_1, index_2, npolyline
    integer                             :: n, np0, np, j, nf
    integer,dimension(2)                :: face
    !f2py integer*1, dimension(1000)    :: P
    type(point)                         :: P
    real(rp),dimension(3)               :: M
    !f2py integer*1, dimension(1000)    :: axisym
    type(axisym)                        :: axisym
    !f2py integer*1, dimension(1000)    :: polyline
    type(GPolyline)                     :: polyline
    
    ! This subroutine creates the geometry of an axisymetric body.
    
    ierror = 0
    
    nface = nface+1
    
    ! Read and creation of the axisym geometry
    call read_axisym(filename,axisym)
    axisym%index = nface
    
    rep1%origine = axisym%repere%origine   
    rep1%e1 = axisym%repere%e1
    rep1%e2 = axisym%repere%e3
    rep1%e3 = -axisym%repere%e2
    
    axisym%vmin = minval(axisym%P(:,2))
    axisym%vmax = maxval(axisym%P(:,2))
    
    axisym%umin = 0._rp
    axisym%umax = PI
    
    n = axisym%npoint
    np0 = 0
    do j=1,n
        if(abs(axisym%P(j,2)).lt.Epsilon)then
            np0 = np0+1
        endif
    enddo
    
    if(Symmetry)then
        call init_geom(fgeom,rep,[10,11,1],[1,2,4-np0],3,ierror) ! 1 axisym, 2 polyline, 4-np0 points
    else ! No symmetry
        call init_geom(fgeom,rep,[10,11,1],[2,2,4-np0],3,ierror) ! 1 axisym, 2 polyline, 4-np0 points
    endif
    
    fgeom%naxisym = fgeom%naxisym+1
    fgeom%axisym(fgeom%naxisym) = axisym
    
    if(not(Symmetry))then
        nface = nface+1
        axisym%index = nface
        fgeom%naxisym = fgeom%naxisym+1
        axisym%repere%e1 = -axisym%repere%e1
        axisym%repere%e2 = -axisym%repere%e2
        fgeom%axisym(fgeom%naxisym) = axisym
    endif     
    
    ! Creation of polylines
    
    nf = 2    

    if(Symmetry)then
        face = [fgeom%axisym(1)%index,-1]
    else
        face = [fgeom%axisym(1)%index,fgeom%axisym(2)%index]
    endif
    
    ! If one wave start with a cylindrical domain, nline = 2
    nline = nline+1 ! nline = 3 in the same conditions
    polyline%index = nline
    index_1 = nline
    polyline%repere = rep1
    polyline%n = n  
    polyline%nface = nf
    polyline%face(1:nf) = face(1:nf)  
    
    allocate(polyline%P(n,3))
    polyline%P(:,3) = 0._rp
    polyline%P(:,1:2) = axisym%P(:,:)
    
    npolyline = fgeom%npolyline
    fgeom%polyline(npolyline+1) = polyline ! 1st polyline
    
    nline = nline+1
    polyline%P(:,1) = -polyline%P(:,1)
    polyline%index = nline ! nline = 4 in the same conditions
    index_2 = nline
    fgeom%polyline(npolyline+2) = polyline ! 2nd polyline
    
    fgeom%npolyline = npolyline+2
    
    ! Creation of corresponding points

    np = fgeom%npoint
    P%nface = nf
    P%face(1:nf) = face(1:nf)
    
    ! First point
    M = polyline%P(1,:)
    np = np+1
    call loc2cart(M,rep1,P%coord)
    P%nedge = 1
    P%edge(1) = index_1
    fgeom%point(np) = P
    if(abs(M(1)).lt.Epsilon)then
        fgeom%point(np)%nedge = 2
        fgeom%point(np)%edge(2) = index_2
    else
        np = np+1
        M(1) = -M(1)
        call loc2cart(M,rep1,P%coord)
        P%nedge = 1
        P%edge(1) = index_2
        fgeom%point(np) = P
    endif

    ! Last point
    M = polyline%P(n,:)
    np = np+1
    call loc2cart(M,rep1,P%coord)
    P%nedge = 1
    P%edge(1) = index_1
    fgeom%point(np) = P
    if(abs(M(1)).lt.Epsilon)then
        fgeom%point(np)%nedge = 2
        fgeom%point(np)%edge(2) = index_2
    else
        np = np+1
        M(1) = -M(1)
        call loc2cart(M,rep1,P%coord)
        P%nedge = 1
        P%edge(1) = index_2
        fgeom%point(np) = P
    endif

    fgeom%npoint = np

9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
99 format('** error #',i3,' : cannot generate geometrie axisym')

    deallocate(polyline%P)

end subroutine wavestar2

end module GeomAxiSym