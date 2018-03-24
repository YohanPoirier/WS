module GeomGen
use GeomStruct
use GeomFonct
use GeomSphere
use GeomCone
use GeomDisque
use GeomCylindre
use GeomPlan
implicit none 

contains

! ************************************************************
! Generateurs d'objets
! ************************************************************

! ------------------------------------------------------------
! Cylindre
! ------------------------------------------------------------
subroutine cylindre_from_dist(rep,L,r,nface,nline,isens,iflag,fgeom)
  
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(inout)        :: rep          ! Frame.
    real(rp),intent(in)                 :: L,r          ! Length and radius of the cylinder.
    !f2py integer*1, dimension(1000)    :: fgeom
    type(type_geom),intent(out)         :: fgeom        ! Geometry of the cylinder.
    integer,intent(inout)               :: nface,nline  ! Number of faces and lines.
    integer,intent(in)                  :: isens        ! Sens (clockwise or not).
    integer,intent(in)                  :: iflag        ! iflag = 0 for the domain, 1 for the bodies
    
    !f2py integer*1, dimension(1000)    :: P1,P2,P3,P4
    type(point)                         :: P1,P2,P3,P4  ! Points.
    !f2py integer*1, dimension(1000)    :: cyl
    type(cylindre_2)                    :: cyl          ! Body of the cylinder.
    !f2py integer*1, dimension(1000)    :: disque
    type(disque2)                       :: disque       ! Bottom and top Disces.
    !f2py integer*1, dimension(1000)    :: cercle
    type(GCercle)                       :: cercle       ! Bottom and top circles.
    integer                             :: iface,iline  ! Face and line number.
    integer                             :: k,narete,j   ! Loop parameter.
    real(rp),dimension(3)               :: pos          ! Position of the origin of the frame.
    integer                             :: ierror       ! Error flag.
    
    ! This subroutine creates the geometry of a cylinder.    
        
    if(Symmetry)then
        call init_geom(fgeom,rep,[4,5,3,9,1,2],[1,2,2,6,4,4],6,ierror)
        ! 1 cylinder
        ! 2 disc
        ! 2 circles
        ! 6 connections
        ! 4 points
        ! 4 aretes
    else
        select case(iflag)
        case(0) ! Domain
            call init_geom(fgeom,rep,[4,5,3,9],[1,2,2,2],4,ierror)
            ! 1 cylinder
            ! 2 disc
            ! 2 circles
            ! 2 connections
        case(1) ! Bodies
            call init_geom(fgeom,rep,[4,5,3,9,1,2],[1,2,2,6,4,2],6,ierror)
            ! 1 cylinder
            ! 2 disc
            ! 2 circles
            ! 6 connections
            ! 4 points to catch the intersection curve
            ! 2 edges to catch the intersection curve
        end select
    endif
    
    pos = rep%origine
    
    iface = nface
    iline = nline
    
    ! Top of the cylinder.
    iface = iface + 1 ! iface = 1 for the free surface in case of a cylindrical tank.
    
    disque%repere%e1 = isens*rep%e1
    disque%repere%e2 = rep%e2
    disque%repere%e3 = isens*rep%e3
    disque%repere%origine = pos + L/2.*rep%e3
    disque%r2max = r*r
    disque%index = iface
    k = fgeom%ndisque
    fgeom%disque(k+1) = disque
    fgeom%ndisque = k+1
    
    ! Bottom of the cylinder.
    iface = iface + 1 ! iface = 2 for the sea bed in case of a cylindrical tank.
    
    disque%repere%e1 = -isens*rep%e1
    disque%repere%e2 = rep%e2
    disque%repere%e3 = -isens*rep%e3
    disque%repere%origine = pos - L/2.*rep%e3
    disque%r2max = r*r
    disque%index = iface
    k = fgeom%ndisque
    fgeom%disque(k+1) = disque
    fgeom%ndisque = k+1  
    
    ! Body of the cylinder.
    iface = iface + 1 ! iface = 3 for the tank in case of a cylindrical tank.
    
    cyl%repere%e1 = rep%e1
    cyl%repere%e2 = rep%e2
    cyl%repere%e3 = isens*rep%e3
    cyl%repere%origine = pos
    cyl%long = L
    cyl%rayon = r  
    cyl%vmin = -L/2.
    cyl%vmax = L/2.
    cyl%index = iface
    k=fgeom%ncylindre
    fgeom%cylindre(k+1) = cyl
    fgeom%ncylindre = k+1  
    
    iline = iline+1
    
    cercle%repere%e1 = isens*rep%e1
    cercle%repere%e2 = rep%e2
    cercle%repere%e3 = isens*rep%e3
    cercle%repere%origine = pos + L/2.*rep%e3
    cercle%rayon = r
    cercle%index = iline
    k = fgeom%ncercle
    fgeom%cercle(k+1) = cercle
    fgeom%ncercle = k+1
    
    iline = iline+1

    cercle%repere%e1 = isens*rep%e1
    cercle%repere%e2 = rep%e2
    cercle%repere%e3 = isens*rep%e3
    cercle%repere%origine = pos - L/2.*rep%e3
    cercle%rayon = r
    cercle%index = iline
    k = fgeom%ncercle
    fgeom%cercle(k+1) = cercle
    fgeom%ncercle = k+1
    
    fgeom%nconnect = 2
    if(abs(fgeom%disque(2)%repere%origine(3)+Ldom(3)).lt.Epsilon .and. iflag==1 .and. cuve_ferme)then
        fgeom%connect(2,:) = [fgeom%cylindre(1)%index,2]
    else
        fgeom%connect(2,:) = [fgeom%cylindre(1)%index,fgeom%disque(2)%index]
    endif
    fgeom%connect(1,:) = [fgeom%cylindre(1)%index,fgeom%disque(1)%index]
    fgeom%cercle(1)%face = [fgeom%cylindre(1)%index,fgeom%disque(1)%index]    
    fgeom%cercle(2)%face = fgeom%connect(2,1:2)

    if(Symmetry)then
        
        P1 = -r*rep%e1 - 0.5_rp*L*rep%e3
        P2 = P1%coord + 2._rp*r*rep%e1
        P3 = P1%coord + L*rep%e3
        P4 = P3%coord + 2._rp*r*rep%e1

        P1%nface = 3
        P1%face(1:3) = [fgeom%disque(2)%index,fgeom%cylindre(1)%index,-1]
        P2%nface = 3
        P2%face(1:3) = [fgeom%disque(2)%index,fgeom%cylindre(1)%index,-1]
        P3%nface = 3
        P3%face(1:3) = [fgeom%disque(1)%index,fgeom%cylindre(1)%index,-1]
        P4%nface = 3
        P4%face(1:3) = [fgeom%disque(1)%index,fgeom%cylindre(1)%index,-1]

        fgeom%npoint = 4
        fgeom%point(1:4) = [P1,P2,P3,P4]

        narete = fgeom%narete
        call arete_from_point(P1,P3,fgeom%arete(1))
        call arete_from_point(P2,P4,fgeom%arete(2))
        call arete_from_point(P3,P4,fgeom%arete(3))
        call arete_from_point(P1,P2,fgeom%arete(4))
        
        do j=1,4
            fgeom%arete(j)%index = iline+j
        enddo
        iline = iline+4

        if(.not.bottom_sym .or. iflag==1)then
            fgeom%narete = narete+4
        else
            fgeom%narete = narete+3
        endif
    
        fgeom%nconnect = iline
        fgeom%connect(3,1:2) = [fgeom%cylindre(1)%index,-1]
        fgeom%arete(1)%face = [fgeom%cylindre(1)%index,-1]
        fgeom%connect(4,1:2) = [fgeom%cylindre(1)%index,-1]
        fgeom%arete(2)%face = [fgeom%cylindre(1)%index,-1]
        fgeom%connect(5,1:2) = [fgeom%disque(1)%index,-1]
        fgeom%arete(3)%face = [fgeom%disque(1)%index,-1]
        fgeom%connect(6,1:2) = [fgeom%disque(2)%index,-1]
        fgeom%arete(4)%face = [fgeom%disque(2)%index,-1]

        fgeom%point(1)%nedge = 3
        fgeom%point(1)%edge(1:3) = [3,2,6]+nline
        fgeom%point(2)%nedge = 3
        fgeom%point(2)%edge(1:3) = [2,4,6]+nline
        fgeom%point(3)%nedge = 3
        fgeom%point(3)%edge(1:3) = [3,1,5]+nline
        fgeom%point(4)%nedge = 3
        fgeom%point(4)%edge(1:3) = [1,4,5]+nline
    else ! No symmetry
        if(iflag == 1)then ! Geometry
            
            ! Points to catch the intersection curve
            !P1 = -r*rep%e1 - 0.5_rp*L*rep%e3
            P1 = r*rep%e2 - 0.5_rp*L*rep%e3
            P3 = P1%coord + L*rep%e3
            
            P1%nface = 3
            P1%face(1:3) = [fgeom%disque(2)%index,fgeom%cylindre(1)%index,-1]
            
            P3%nface = 3
            P3%face(1:3) = [fgeom%disque(1)%index,fgeom%cylindre(1)%index,-1]
            
            fgeom%npoint = 2
            fgeom%point(1:2) = [P1,P3]
            
            ! Edge to catch the intersection curve
            narete = fgeom%narete
            fgeom%narete = narete+1
            call arete_from_point(P1,P3,fgeom%arete(1))
            
            fgeom%arete(1)%index = iline+1
            iline = iline+1
            
            fgeom%nconnect = iline
            fgeom%connect(3,1:2) = [fgeom%cylindre(1)%index,-1]
            fgeom%arete(1)%face = [fgeom%cylindre(1)%index,-1]
            
            fgeom%point(1)%nedge = 3
            fgeom%point(1)%edge(1:3) = [3,2,6]+nline
            fgeom%point(2)%nedge = 3
            fgeom%point(2)%edge(1:3) = [2,4,6]+nline
            
            ! New ---
            narete = fgeom%narete
            fgeom%narete = narete+1
            
            P2 = -r*rep%e2 - 0.5_rp*L*rep%e3
            P4 = P2%coord + L*rep%e3
            
            P2%nface = 3
            P2%face(1:3) = [fgeom%disque(2)%index,fgeom%cylindre(1)%index,-1]
            
            P4%nface = 3
            P4%face(1:3) = [fgeom%disque(1)%index,fgeom%cylindre(1)%index,-1]
            
            fgeom%npoint = 4
            fgeom%point(1:4) = [P1,P3,P2,P4]
            
            call arete_from_point(P2,P4,fgeom%arete(2))
            
            fgeom%arete(2)%index = iline+1
            iline = iline+1
            
            fgeom%nconnect = iline
            fgeom%connect(3,1:2) = [fgeom%cylindre(1)%index,-1]
            fgeom%arete(2)%face = [fgeom%cylindre(1)%index,-1]
            
            fgeom%point(3)%nedge = 3
            fgeom%point(3)%edge(1:3) = [3,2,6]+nline
            fgeom%point(4)%nedge = 3
            fgeom%point(4)%edge(1:3) = [2,4,6]+nline
            
            
            fgeom%cercle(1)%P1 = r*rep%e1 - 0.5_rp*L*rep%e3
            fgeom%cercle(1)%P1%nface = 2
            fgeom%cercle(1)%P1%edge(1:2) = [fgeom%cylindre(1)%index,fgeom%disque(1)%index]
            
            fgeom%cercle(1)%P2 = -r*rep%e1 - 0.5_rp*L*rep%e3
            fgeom%cercle(1)%P2%nface = 2
            fgeom%cercle(1)%P2%edge(1:2) = [fgeom%cylindre(1)%index,fgeom%disque(1)%index]
            
            fgeom%cercle(2)%P1 = r*rep%e1 + 0.5_rp*L*rep%e3
            fgeom%cercle(2)%P1%nface = 2
            fgeom%cercle(2)%P1%edge(1:2) = [fgeom%cylindre(1)%index,fgeom%disque(2)%index]
            
            fgeom%cercle(2)%P2 = -r*rep%e1 + 0.5_rp*L*rep%e3
            fgeom%cercle(2)%P2%nface = 2
            fgeom%cercle(2)%P2%edge(1:2) = [fgeom%cylindre(1)%index,fgeom%disque(2)%index]
            
            ! ---
            
        end if
    endif

    nface = iface
    nline = iline

end subroutine cylindre_from_dist

! ---------------------------------------------------------
! Disque
! ---------------------------------------------------------

subroutine disque_from_dist(rep,r,nface,isens,fgeom)
    
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),intent(in) :: r
    !f2py integer*1, dimension(1000) :: fgeom
    type(type_geom),intent(out) :: fgeom
    integer,intent(inout) :: nface
    integer,intent(in) :: isens
    
    !f2py integer*1, dimension(1000) :: disque
    type(disque2) :: disque
    !f2py integer*1, dimension(1000) :: cercle
    type(GCercle) :: cercle
    integer :: iface
    integer :: k
    real(rp),dimension(3) :: pos
    
    ! This subroutine creates the geometry of a disc.
    
    allocate(fgeom%disque(1))
    allocate(fgeom%cercle(1))
  
    fgeom%repere = rep
    fgeom%npoint = 0  
    fgeom%narete = 0
    fgeom%ncercle = 0
    fgeom%ncylindre = 0
    fgeom%ndisque = 0
    fgeom%nplan = 0
    fgeom%nconnect = 0
    
    pos = rep%origine
    
    iface = nface
    
    iface = iface+1
    
    disque%repere%e1 = isens*rep%e1
    disque%repere%e2 = rep%e2
    disque%repere%e3 = isens*rep%e3
    disque%repere%origine = pos
    disque%r2max = r*r
    disque%index = iface
    k = fgeom%ndisque
    fgeom%disque(k+1) = disque
    fgeom%ndisque = k+1        
    
    cercle%repere%e1 = isens*rep%e1
    cercle%repere%e2 = rep%e2
    cercle%repere%e3 = isens*rep%e3
    cercle%repere%origine = pos
    cercle%rayon = r
    cercle%index = 1
    k = fgeom%ncercle
    fgeom%cercle(k+1) = cercle
    fgeom%ncercle = k+1
  
    nface = iface
    allocate(fgeom%cmd_dis(1))
    fgeom%cmd_dis(1) = 1   
  
end subroutine  disque_from_dist

!----------------------------------------------------------------
! Cube
! ---------------------------------------------------------------

subroutine cube_from_length(rep,L,nface,fgeom)
    
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),intent(in) :: L
    integer,intent(inout) :: nface
    !f2py integer*1, dimension(1000) :: fgeom
    type(type_geom),intent(out) :: fgeom
    ! local
    !f2py integer*1, dimension(1000) :: O
    type(point) :: O
    !f2py integer*1, dimension(1000) :: e1,e2,e3
    type(vector):: e1,e2,e3
    !
    !f2py integer*1, dimension(1000) :: P1,P2,P3,P4,P5,P6,P7,P8
    type(point) :: P1,P2,P3,P4,P5,P6,P7,P8
    !f2py integer*1, dimension(1000) :: lplan
    type(Gplan),dimension(6) :: lplan
    !f2py integer*1, dimension(1000) :: larete
    type(Garete),dimension(12) :: larete
    integer :: nplan,narete
    integer :: j 
    
    ! This subroutine creates the geometry of a cube.
    
    allocate(fgeom%plan(6))
    allocate(fgeom%arete(12))
    allocate(fgeom%point(8))
    allocate(fgeom%connect(12,2))

    O = rep%origine
    call assign_vector_coord(e1,rep%e1)
    call assign_vector_coord(e2,rep%e2)
    call assign_vector_coord(e3,rep%e3)
    
    fgeom%repere = rep
    fgeom%npoint = 0
    fgeom%narete = 0
    fgeom%ncercle = 0
    fgeom%ncylindre = 0
    fgeom%ndisque = 0
    fgeom%nplan = 0
    fgeom%nconnect = 0
    
    P1 = O%coord - L/2.*(e1%coord+e2%coord+e3%coord)
    P2 = P1%coord + L*e1%coord
    P3 = P2%coord + L*e2%coord
    P4 = P3%coord - L*e1%coord
    P5 = P1%coord + L*e3%coord
    P6 = P5%coord + L*e1%coord
    P7 = P6%coord + L*e2%coord
    P8 = P7%coord - L*e1%coord
    
    call plan_from_point(P1,P4,P3,P2,lplan(1))
    lplan(1)%index = nface+1
    call plan_from_point(P5,P6,P7,P8,lplan(2))
    lplan(2)%index = nface+2
    call plan_from_point(P1,P2,P6,P5,lplan(3))
    lplan(3)%index = nface+3
    call plan_from_point(P2,P3,P7,P6,lplan(4))
    lplan(4)%index = nface+4
    call plan_from_point(P3,P4,P8,P7,lplan(5))
    lplan(5)%index = nface+5
    call plan_from_point(P4,P1,P5,P8,lplan(6))
    lplan(6)%index = nface+6
    nface = nface+6
    
    nplan = fgeom%nplan
    do j=1,6
        fgeom%plan(nplan+j) = lplan(j)
    enddo
    fgeom%nplan = nplan+6
    
    P1%nface=3
    P1%face(1:3) = [lplan(1)%index,lplan(3)%index,lplan(6)%index]
    P2%nface=3
    P2%face(1:3) = [lplan(1)%index,lplan(3)%index,lplan(4)%index]
    P3%nface=3
    P3%face(1:3) = [lplan(1)%index,lplan(4)%index,lplan(5)%index]
    P4%nface=3
    P4%face(1:3) = [lplan(1)%index,lplan(5)%index,lplan(6)%index]
    P5%nface=3
    P5%face(1:3) = [lplan(2)%index,lplan(3)%index,lplan(6)%index]
    P6%nface=3
    P6%face(1:3) = [lplan(2)%index,lplan(3)%index,lplan(4)%index]
    P7%nface=3
    P7%face(1:3) = [lplan(2)%index,lplan(4)%index,lplan(5)%index]
    P8%nface=3
    P8%face(1:3) = [lplan(2)%index,lplan(5)%index,lplan(6)%index]  
    
    fgeom%point(1) = P1
    fgeom%point(2) = P2
    fgeom%point(3) = P3
    fgeom%point(4) = P4
    fgeom%point(5) = P5
    fgeom%point(6) = P6
    fgeom%point(7) = P7
    fgeom%point(8) = P8
    fgeom%npoint = 8

    call arete_from_point(P1,P2,larete(1))
    call arete_from_point(P2,P3,larete(2))
    call arete_from_point(P3,P4,larete(3))
    call arete_from_point(P4,P1,larete(4))
    call arete_from_point(P1,P5,larete(5))
    call arete_from_point(P2,P6,larete(6))
    call arete_from_point(P3,P7,larete(7))
    call arete_from_point(P4,P8,larete(8))
    call arete_from_point(P5,P6,larete(9))
    call arete_from_point(P6,P7,larete(10))
    call arete_from_point(P7,P8,larete(11))
    call arete_from_point(P8,P5,larete(12))
  
    
    narete = fgeom%narete
    do j=1,12
        fgeom%arete(narete+j) = larete(j)
        fgeom%arete(narete+j)%index = narete+j
    enddo
    fgeom%narete = narete+12    
    
    fgeom%nconnect = 12
    fgeom%connect(1,:)  = [lplan(1)%index,lplan(3)%index]
    fgeom%connect(2,:)  = [lplan(1)%index,lplan(4)%index]
    fgeom%connect(3,:)  = [lplan(1)%index,lplan(5)%index]
    fgeom%connect(4,:)  = [lplan(1)%index,lplan(6)%index]
    fgeom%connect(5,:)  = [lplan(3)%index,lplan(6)%index]
    fgeom%connect(6,:)  = [lplan(3)%index,lplan(4)%index]
    fgeom%connect(7,:)  = [lplan(4)%index,lplan(5)%index]
    fgeom%connect(8,:)  = [lplan(5)%index,lplan(6)%index]
    fgeom%connect(9,:)  = [lplan(2)%index,lplan(3)%index]
    fgeom%connect(10,:) = [lplan(2)%index,lplan(4)%index]
    fgeom%connect(11,:) = [lplan(2)%index,lplan(5)%index]
    fgeom%connect(12,:) = [lplan(2)%index,lplan(6)%index]
    
    allocate(fgeom%cmd_pl(6))
    fgeom%cmd_pl = [0,0,0,0,0,0] 
  
end subroutine cube_from_length 

! -----------------------------------------------------------------
! Rectangle
! -----------------------------------------------------------------

subroutine rectangle_from_length(rep,LL,L,H,nface,is_close,fgeom)
    
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),intent(in) :: LL,L,H
    integer,intent(inout) :: nface
    logical,intent(in) :: is_close
    !f2py integer*1, dimension(1000) :: fgeom
    type(type_geom),intent(out) :: fgeom
    
    !f2py integer*1, dimension(1000) :: O
    type(point) :: O
    !f2py integer*1, dimension(1000) :: e1,e2,e3
    type(vector) :: e1,e2,e3
    !f2py integer*1, dimension(1000) :: P1,P2,P3,P4,P5,P6,P7,P8
    type(point) :: P1,P2,P3,P4,P5,P6,P7,P8
    !f2py integer*1, dimension(1000) :: lplan
    type(Gplan),dimension(6) :: lplan
    !f2py integer*1, dimension(1000) :: larete
    type(Garete),dimension(12) :: larete
    integer :: nplan,narete
    integer :: j,np,nl
    
    ! This subroutine creates the geometry of a rectangle.

    if(.not.is_close)then
        np = 1 ; nl = 4
        allocate(fgeom%plan(1))
        allocate(fgeom%arete(4))
        allocate(fgeom%point(4))
        allocate(fgeom%connect(4,2))
        allocate(fgeom%cmd_pl(1))
        elseif(is_close)then
        np = 6 ; nl = 12
        allocate(fgeom%plan(6))
        allocate(fgeom%arete(12))
        allocate(fgeom%point(12))
        allocate(fgeom%connect(12,2))
        allocate(fgeom%cmd_pl(6))
    endif
    
    O = rep%origine
    call assign_vector_coord(e1,rep%e1)
    call assign_vector_coord(e2,rep%e2)
    call assign_vector_coord(e3,rep%e3)
    
    fgeom%repere = rep
    fgeom%npoint = 0
    fgeom%narete = 0
    fgeom%ncercle = 0
    fgeom%ncylindre = 0
    fgeom%ndisque = 0
    fgeom%nplan = 0
    fgeom%nconnect = 0
    
    P1 =  - 0.5*LL*e1%coord -0.5*L*e2%coord - 0.5*H*e3%coord
    P2 = P1%coord + LL*e1%coord  
    P3 = P2%coord + L*e2%coord
    P4 = P3%coord - LL*e1%coord
    P5 = P1%coord + H*e3%coord
    P6 = P5%coord + LL*e1%coord
    P7 = P6%coord + L*e2%coord
    P8 = P7%coord - LL*e1%coord
  
    if(Symmetry)then
        P1%coord(2) = 0._rp
        P2%coord(2) = 0._rp
        P5%coord(2) = 0._rp
        P6%coord(2) = 0._rp
    endif
  
    call plan_from_point(P5,P6,P7,P8,lplan(1))
    lplan(1)%index = nface+1
    lplan(1)%repere%origine = lplan(1)%repere%origine + O%coord
    lplan(1)%repere%e1 = -lplan(1)%repere%e1
    lplan(1)%repere%e3 = -lplan(1)%repere%e3
    nface = nface+1
    
    if(is_close)then
        call plan_from_point(P1,P2,P3,P4,lplan(2))
        lplan(2)%index = nface+1
        lplan(2)%repere%origine = lplan(2)%repere%origine + O%coord  
    
        call plan_from_point(P1,P5,P6,P2,lplan(3))
        lplan(3)%index = nface+2
        lplan(3)%repere%origine = lplan(3)%repere%origine + O%coord
    
        call plan_from_point(P2,P6,P7,P3,lplan(4))
        lplan(4)%index = nface+3
        lplan(4)%repere%origine = lplan(4)%repere%origine + O%coord
    
        call plan_from_point(P3,P7,P8,P4,lplan(5))
        lplan(5)%index = nface+4
        lplan(5)%repere%origine = lplan(5)%repere%origine + O%coord
    
        call plan_from_point(P4,P8,P5,P1,lplan(6))
        lplan(6)%index = nface+5
        lplan(6)%repere%origine = lplan(6)%repere%origine + O%coord
        nface = nface+5
    endif
    
    nplan = fgeom%nplan
    do j=1,np
        fgeom%plan(nplan+j) = lplan(j)
    enddo
    fgeom%nplan = nplan+np
    
    P5%nface = 1
    P5%face(1) = lplan(1)%index
    P6%nface = 1
    P6%face(1) = lplan(1)%index
    P7%nface = 1
    P7%face(1) = lplan(1)%index
    P8%nface = 1
    P8%face(1) = lplan(1)%index
  
    if(is_close)then 
        P1%nface=3
        P1%face(1:3) = [lplan(2)%index,lplan(3)%index,lplan(6)%index]
        P2%nface=3
        P2%face(1:3) = [lplan(2)%index,lplan(3)%index,lplan(4)%index]
        P3%nface=3
        P3%face(1:3) = [lplan(2)%index,lplan(4)%index,lplan(5)%index]
        P4%nface=3
        P4%face(1:3) = [lplan(2)%index,lplan(5)%index,lplan(6)%index]
        P5%nface=3
        P5%face(1:3) = [lplan(1)%index,lplan(3)%index,lplan(6)%index]
        P6%nface=3
        P6%face(1:3) = [lplan(1)%index,lplan(3)%index,lplan(4)%index]
        P7%nface=3
        P7%face(1:3) = [lplan(1)%index,lplan(4)%index,lplan(5)%index]
        P8%nface=3
        P8%face(1:3) = [lplan(1)%index,lplan(5)%index,lplan(6)%index]  
    endif
    
    fgeom%point(1) = P5
    fgeom%point(2) = P6
    fgeom%point(3) = P7
    fgeom%point(4) = P8
    fgeom%npoint = 4
    
    if(is_close)then
        fgeom%point(5) = P1
        fgeom%point(6) = P2
        fgeom%point(7) = P3
        fgeom%point(8) = P4
        fgeom%npoint = 8
    endif
    
    call arete_from_point(P5,P6,larete(1))
    call arete_from_point(P6,P7,larete(2))
    call arete_from_point(P7,P8,larete(3))
    call arete_from_point(P8,P5,larete(4))
    
    if(is_close)then
        call arete_from_point(P1,P2,larete(5))
        call arete_from_point(P2,P3,larete(6))
        call arete_from_point(P3,P4,larete(7))
        call arete_from_point(P4,P1,larete(8))
        call arete_from_point(P1,P5,larete(9))
        call arete_from_point(P2,P6,larete(10))
        call arete_from_point(P3,P7,larete(11))
        call arete_from_point(P4,P8,larete(12))
    endif
    
    narete = fgeom%narete
    do j=1,nl
        fgeom%arete(narete+j) = larete(j)
        fgeom%arete(narete+j)%index = narete+j
    enddo
    fgeom%narete = narete+nl    
    
    fgeom%nconnect = nl
    fgeom%connect(1,:)  = [lplan(1)%index,lplan(1)%index]
    fgeom%connect(2,:)  = [lplan(1)%index,lplan(1)%index]
    fgeom%connect(3,:)  = [lplan(1)%index,lplan(1)%index]
    fgeom%connect(4,:)  = [lplan(1)%index,lplan(1)%index]
    if(is_close)then
        fgeom%connect(1,:)  = [lplan(1)%index,lplan(3)%index]
        fgeom%connect(2,:)  = [lplan(1)%index,lplan(4)%index]
        fgeom%connect(3,:)  = [lplan(1)%index,lplan(5)%index]
        fgeom%connect(4,:)  = [lplan(1)%index,lplan(6)%index]
        fgeom%connect(5,:)  = [lplan(3)%index,lplan(6)%index]
        fgeom%connect(6,:)  = [lplan(3)%index,lplan(4)%index]
        fgeom%connect(7,:)  = [lplan(4)%index,lplan(5)%index]
        fgeom%connect(8,:)  = [lplan(5)%index,lplan(6)%index]
        fgeom%connect(9,:)  = [lplan(2)%index,lplan(3)%index]
        fgeom%connect(10,:) = [lplan(2)%index,lplan(4)%index]
        fgeom%connect(11,:) = [lplan(2)%index,lplan(5)%index]
        fgeom%connect(12,:) = [lplan(2)%index,lplan(6)%index]
    endif
    
    fgeom%cmd_pl = 0 
        
end subroutine rectangle_from_length

subroutine geom_wavestar0(rep,nface,nline,geom)
    implicit none
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    integer,intent(inout) :: nface,nline
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(out) :: geom
!   local
    !f2py integer*1, dimension(1000) :: P1,P2,P3,P4,P5,P6
    type(point) :: P1,P2,P3,P4,P5,P6
    real(rp),dimension(3) :: C, pos
    real(rp),dimension(3) :: normale
    real(rp) :: r,h,angle,phi1,phi2
    real(rp) :: scale
    integer :: narete,j,nline0
    integer :: ierror

    ierror = 0

    nline0 = nline

    scale = 20._rp

    if(symmetry)then
        call init_geom(geom,rep,[8,7,5,3,9,1,2],[1,2,1,4,9,6,5],7,ierror)
    else
        call init_geom(geom,rep,[8,7,5,3,9],[1,2,1,3,3],5,ierror)
    endif

    pos = scale*[0._rp,0._rp,-0.055_rp]

    C = scale*[0._rp,0._rp,-0.055_rp]-pos
    normale = [0._rp,0._rp,-1._rp]
    r = scale*0.125_rp
    phi2 = 0.5_rp*pi
    phi1 = 0.1221434839_rp
    call assign_sphere(geom,rep,C,normale,r,phi1,phi2,nface,ierror)

    C = scale*[0._rp,0._rp,-0.00975_rp]-pos
    normale = [0._rp,0._rp,-1._rp]
    r = scale*0.1317_rp
    angle = 0.12549443288_rp
    h = scale*0.06048_rp
    call assign_cone(geom,rep,C,normale,r,angle,h,nface,ierror)

    C = scale*[0._rp,0._rp,-0.00975_rp]-pos
    normale = [0._rp,0._rp,1._rp]
    r = scale*0.1317_rp
    angle = 1.3843975277
    h = scale*0.00975_rp
    call assign_cone(geom,rep,C,normale,r,angle,h,nface,ierror)

    C = scale*[0._rp,0._rp,0._rp]-pos
    normale = [0._rp,0._rp,1._rp]
    r = scale*0.080_rp
    call assign_disque(geom,rep,C,normale,r,nface,ierror)

    normale = [0._rp,0._rp,-1._rp]

    C = scale*[0._rp,0._rp,0._rp]-pos
    r = scale*0.08_rp
    call assign_cercle(geom,rep,C,normale,r,nline,ierror)
    
    C = scale*[0._rp,0._rp,-0.00975_rp]-pos
    r = scale*0.1317_rp
    call assign_cercle(geom,rep,C,normale,r,nline,ierror)
    
    C = scale*[0._rp,0._rp,-0.07023_rp]-pos
    r = scale*0.12407_rp
    call assign_cercle(geom,rep,C,normale,r,nline,ierror)
    
    geom%nconnect = nline
    geom%connect(1,1:2) = [geom%disque(1)%index,geom%cone(2)%index]
    geom%cercle(1)%face = [geom%disque(1)%index,geom%cone(2)%index]
    geom%connect(2,1:2) = [geom%cone(1)%index,geom%cone(2)%index]
    geom%cercle(2)%face = [geom%cone(1)%index,geom%cone(2)%index]
    geom%connect(3,1:2) = [geom%cone(1)%index,geom%sphere(1)%index] 
    geom%cercle(3)%face = [geom%cone(1)%index,geom%sphere(1)%index] 

    if(symmetry)then

        P1 = scale*[-0.12407_rp,0._rp,-0.07023_rp]-pos
        P2 = scale*[0.12407_rp,0._rp,-0.07023_rp]-pos
        P3 = scale*[-0.1317_rp,0._rp,-0.00975_rp]-pos
        P4 = scale*[0.1317_rp,0._rp,-0.00975_rp]-pos
        P5 = scale*[-0.08_rp,0._rp,0._rp]-pos
        P6 = scale*[0.08_rp,0._rp,0._rp]-pos

        P1%nface = 3
        P1%face(1:3) = [geom%sphere(1)%index,geom%cone(1)%index,-1]
        P1%bf = 2
        P2%nface = 3
        P2%face(1:3) = [geom%sphere(1)%index,geom%cone(1)%index,-1]
        P2%bf = 2
        P3%nface = 3
        P3%face(1:3) = [geom%cone(1)%index,geom%cone(2)%index,-1]
        P3%bf = 2
        P4%nface = 3
        P4%face(1:3) = [geom%cone(1)%index,geom%cone(2)%index,-1]
        P4%bf = 2
        P5%nface = 3
        P5%face(1:3) = [geom%cone(2)%index,geom%disque(1)%index,-1]
        P5%bf = 2
        P6%nface = 3
        P6%face(1:3) = [geom%cone(2)%index,geom%disque(1)%index,-1] 
        P6%bf = 2

        C = scale*[0._rp,0._rp,-0.055_rp]-pos
        r = scale*0.125_rp
        normale = [0._rp,-1._rp,0._rp]
        call assign_cercle(geom,rep,C,normale,r,nline,ierror,P1,P2)

        geom%npoint = 6
        geom%point(1:6) = [P1,P2,P3,P4,P5,P6]

        narete = geom%narete
        call arete_from_point(P1,P3,geom%arete(narete+1))
        call arete_from_point(P4,P2,geom%arete(narete+2))
        call arete_from_point(P3,P5,geom%arete(narete+3))
        call arete_from_point(P4,P6,geom%arete(narete+4))
        call arete_from_point(P5,P6,geom%arete(narete+5))
        do j=1,5
            geom%arete(narete+j)%index = nline+j
        enddo
        geom%narete = narete+5

        geom%nconnect = nline+6
        geom%connect(4,1:2) = [geom%sphere(1)%index,-1]
        geom%cercle(4)%face = [geom%sphere(1)%index,-1]
        geom%connect(5,1:2) = [geom%cone(1)%index,-1]
        geom%arete(1)%face =  [geom%cone(1)%index,-1]
        geom%connect(6,1:2) = [geom%cone(1)%index,-1]
        geom%arete(2)%face = [geom%cone(1)%index,-1]
        geom%connect(7,1:2) = [geom%cone(2)%index,-1]
        geom%arete(3)%face = [geom%cone(2)%index,-1]
        geom%connect(8,1:2) = [geom%cone(2)%index,-1]
        geom%arete(4)%face = [geom%cone(2)%index,-1]
        geom%connect(9,1:2) = [geom%disque(1)%index,-1]
        geom%arete(5)%face = [geom%disque(1)%index,-1]

        geom%point(1)%nedge = 3
        geom%point(1)%edge(1:3) = [3,4,5]+nline0
        geom%point(2)%nedge = 3
        geom%point(2)%edge(1:3) = [3,4,6]+nline0
        geom%point(3)%nedge = 3
        geom%point(3)%edge(1:3) = [2,5,7]+nline0
        geom%point(4)%nedge = 3
        geom%point(4)%edge(1:3) = [2,6,8]+nline0
        geom%point(5)%nedge = 3
        geom%point(5)%edge(1:3) = [1,7,9]+nline0 
        geom%point(6)%nedge = 3
        geom%point(6)%edge(1:3) = [1,8,9]+nline0

    endif

end subroutine geom_wavestar0
 
subroutine Sphere_from_radius(rep,radius,nface,nline,geom)
    
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep          ! Frame.
    real(rp),intent(in)                 :: radius       ! Radius of the sphere.
    integer,intent(inout)               :: nface,nline  ! Number of faces and lines.
    !f2py integer*1, dimension(1000)    :: geom
    type(type_geom),intent(out)         :: geom         ! Geometry.
    
    real(rp),dimension(3)               :: C,normale    ! Center of the sphere.
    real(rp)                            :: phi1,phi2    ! Angles of the sphere.
    integer                             :: ierror       ! Error flag.

    ! This subroutine creates the geometry of a sphere in case of non symmetry.

    ierror = 0
        
    if(symmetry)then
        print*,"Please use a axisym case for a sphere with a symmetry plan."
        call exit()
    else
        call init_geom(geom,rep,[8],[1],1,ierror) ! 8 = sphere.
    endif
    
    C = [0._RP,0._RP,0._RP]
    normale = [0._rp,0._rp,-1._rp]
    phi1 = -0.5_rp**pi
    phi2 = 0.5_rp*pi
    call assign_sphere(geom,rep,C,normale,radius,phi1,phi2,nface,ierror)

end subroutine Sphere_from_radius

! -----------------------------------------------------------------
! Mise a jour de la position de l'objet
! -----------------------------------------------------------------

subroutine update_geom(rep0,fgeom,angle,pos)

    !f2py integer*1, dimension(1000)    :: rep0
    type(repere3d),intent(in)           :: rep0     ! Inertial frame.
    !f2py integer*1, dimension(1000)    :: fgeom
    type(type_geom),intent(inout)       :: fgeom    ! Geometry.
    real(rp),dimension(3),intent(in)    :: angle    ! Angular position (rad).
    real(rp),dimension(3),intent(in)    :: pos      ! Linear position.

    !f2py integer*1, dimension(1000)    :: rep,repj
    type(repere3d)                      :: rep,repj ! Local frames.
    integer                             :: n        ! Number of sub-geometries.
    integer                             :: j        ! Loop parameter.
    real(rp),dimension(3)               :: M1,M2    ! Points.

	! This subroutine updates the linear and angular position of the geometries.

    !  Passage coordonnees locales.
   rep = fgeom%repere
   n = fgeom%ncercle
   do j=1,n
     repj = fgeom%cercle(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%cercle(j)%repere = repj
   enddo
   n = fgeom%narete
   do j=1,n
     M1 = fgeom%arete(j)%P1%coord
     M2 = fgeom%arete(j)%P2%coord
     call cart2loc(M1,rep,fgeom%arete(j)%P1%coord)
     call cart2loc(M2,rep,fgeom%arete(j)%P2%coord)     
   enddo   
   n = fgeom%ncylindre
   do j=1,n
     repj = fgeom%cylindre(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%cylindre(j)%repere = repj
   enddo
   n = fgeom%ndisque
   do j=1,n
     repj = fgeom%disque(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%disque(j)%repere = repj
   enddo
   n = fgeom%nplan
   do j=1,n
     repj =fgeom%plan(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%plan(j)%repere = repj
   enddo
   n = fgeom%ncone
   do j=1,n
     repj = fgeom%cone(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%cone(j)%repere = repj
   enddo
   n = fgeom%nsphere
   do j=1,n
     repj = fgeom%sphere(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%sphere(j)%repere = repj
   enddo
   n = fgeom%naxisym
   do j=1,n
     repj = fgeom%axisym(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%axisym(j)%repere = repj
   enddo
   n = fgeom%npolyline
   do j=1,n
     repj = fgeom%polyline(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%polyline(j)%repere = repj
   enddo
   n = fgeom%nwigley
   do j=1,n
     repj = fgeom%wigley(j)%repere
     call cart2loc(repj%origine,rep,repj%origine)
     call cart2loc(repj%e1+rep%origine,rep,repj%e1)
     call cart2loc(repj%e2+rep%origine,rep,repj%e2)
     call cart2loc(repj%e3+rep%origine,rep,repj%e3)
     fgeom%wigley(j)%repere = repj
   enddo
   
    ! Mise a jour position repere principal.
    fgeom%repere%origine = pos
    fgeom%repere%e1 = rep0%e1
    fgeom%repere%e2 = rep0%e2
    fgeom%repere%e3 = rep0%e3
    if(abs(angle(1)).gt.Epsilon)then
        call mat_rotation(fgeom%repere%e1,angle(1),1,[0._rp,0._rp,0._rp],fgeom%repere%e1)
        call mat_rotation(fgeom%repere%e3,angle(1),1,[0._rp,0._rp,0._rp],fgeom%repere%e3)
        call mat_rotation(fgeom%repere%e2,angle(1),1,[0._rp,0._rp,0._rp],fgeom%repere%e2)
    end if
    if(abs(angle(2)).gt.Epsilon)then
        call mat_rotation(fgeom%repere%e1,angle(2),2,[0._rp,0._rp,0._rp],fgeom%repere%e1)
        call mat_rotation(fgeom%repere%e2,angle(2),2,[0._rp,0._rp,0._rp],fgeom%repere%e2)
        call mat_rotation(fgeom%repere%e3,angle(2),2,[0._rp,0._rp,0._rp],fgeom%repere%e3)
    end if
    if(abs(angle(3)).gt.Epsilon)then
        call mat_rotation(fgeom%repere%e1,angle(3),3,[0._rp,0._rp,0._rp],fgeom%repere%e1)
        call mat_rotation(fgeom%repere%e2,angle(3),3,[0._rp,0._rp,0._rp],fgeom%repere%e2)
        call mat_rotation(fgeom%repere%e3,angle(3),3,[0._rp,0._rp,0._rp],fgeom%repere%e3)
    end if
        
    ! Mise a jour position repere surface.
    rep = fgeom%repere    
    n = fgeom%ncercle
    do j=1,n
        repj = fgeom%cercle(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%cercle(j)%repere = repj
        
        M1 = fgeom%cercle(j)%P1%coord
        M2 = fgeom%cercle(j)%P2%coord
        call loc2cart(M1,rep,fgeom%cercle(j)%P1%coord)
        call loc2cart(M2,rep,fgeom%cercle(j)%P2%coord)
    end do
    n = fgeom%narete
    do j=1,n        
        M1 = fgeom%arete(j)%P1%coord
        M2 = fgeom%arete(j)%P2%coord
        call loc2cart(M1,rep,fgeom%arete(j)%P1%coord)
        call loc2cart(M2,rep,fgeom%arete(j)%P2%coord)
    end do
    n = fgeom%ncylindre
    do j=1,n
        repj = fgeom%cylindre(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%cylindre(j)%repere = repj        
    end do
    n = fgeom%ndisque
    do j=1,n
        repj = fgeom%disque(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%disque(j)%repere = repj
    end do
    n = fgeom%nplan
    do j=1,n
        repj =fgeom%plan(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%plan(j)%repere = repj
    end do   
    n = fgeom%ncone
    do j=1,n
        repj = fgeom%cone(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%cone(j)%repere = repj
    end do
    n = fgeom%nsphere
    do j=1,n
        repj = fgeom%sphere(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%sphere(j)%repere = repj
    end do
    n = fgeom%naxisym
    do j=1,n
        repj = fgeom%axisym(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%axisym(j)%repere = repj
    end do
    n = fgeom%npolyline
    do j=1,n
        repj = fgeom%polyline(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%polyline(j)%repere = repj
    end do
    n = fgeom%nwigley
    do j=1,n
        repj = fgeom%wigley(j)%repere
        call loc2cart(repj%origine,rep,repj%origine)
        call loc2cart(repj%e1,rep,repj%e1)
        repj%e1 = repj%e1-pos
        call loc2cart(repj%e2,rep,repj%e2)
        repj%e2 = repj%e2-pos
        call loc2cart(repj%e3,rep,repj%e3)
        repj%e3 = repj%e3-pos
        fgeom%wigley(j)%repere = repj
    end do
    
end subroutine update_geom   

subroutine update_position(P_new,P_old,C,vect,angle,ierror)
    
    real(rp),dimension(3),intent(inout) :: P_new
    real(rp),dimension(3),intent(in) :: P_old
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: vect,angle
    integer,intent(inout) :: ierror
    
    integer :: k
    real(rp),dimension(3) :: P, Ptemp
    
    ierror = 0
    
    P = P_old(1:3)
    do k=1,3
        Ptemp = P-C
        Ptemp(k) = 0._RP
        if(norm2(Ptemp).gt.Epsgeom)then
            call mat_rotation(P,angle(k),k,C,P)
        endif
    end do

    P_new(1:3) = P + vect

9999 continue
    if(ierror/=0)then
        write(99,*) ierror
    endif
99 format('** error #',i3,' : pb. update position.')

end subroutine update_position

end module GeomGen
