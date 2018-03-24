module GenMaillage
use Parameters
use Constantes
use FonctionsCommunes
use Structuresdonnees
use Incident_mod
use GeomMesh
use Preplot
implicit none

!integer, parameter :: decalage = 1            ! Mesh à décalage
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                         Création d'un maillage en coordonnées cartésiennes                       !
!                                 LETOURNEL Lucas    03/13                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bodygen(Mesh,fgeom_vect,InputData)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(inout)      :: Mesh                                 ! Mesh
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                           ! Geometries
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                            ! Input data

    integer                             :: Nface
    integer                             :: j                                    ! Loop parameter
    real(rp), dimension(3)              :: Orig, P1, P2, P3, P4, P5, P6, P7, P8 ! Main points
    real(rp)                            :: iabs, nx, ny, nz, Nzero, iabsz       !
    real(rp), dimension(3,18)           :: Facette                              !
    real(rp), dimension(4,6)            :: n                                    ! Nodes number table
    integer, dimension(6)               :: typeFrontiere                        !
    
    ! This subroutine creates a mesh from a cartesian mesh of the tank.
    
    ! Number of nodes in each direction.
    Nzero=0._RP
    Ldom(4) = Ldom(5)
    nx=float(int(Ldom(1)/Ldom(4)))+1 ! Number of nodes in x direction.
    ny=float(int(Ldom(2)/Ldom(4)))+1 ! Number of nodes in y direction.
    nz=float(int(Ldom(3)/Ldom(4)))+1 ! Number of nodes in z direction.
    print*, "Rectangular mesh of the domain"
    Write(6,'("Length : ", F5.1," , Width : ", F5.1," , Depth : ", F5.1)') Ldom(1:3)
    Write(6,'("nx : ", F5.1," , ny : ", F5.1," , nz : ", F5.1)') nx,ny,nz
    
    ! Definition of the origin (on the free surface) 
    if (Symmetry) then
        Orig=[-Ldom(1)/2._RP,0._RP,0._RP]
    else
        Orig=[-Ldom(1)/2._RP,-Ldom(2)/2._RP,0._RP]
    end if
    
    iabs = 1.25_RP
    iabsz = 2._RP
    
    if (not(cuve_ferme)) nz = 1
    
    ! Points to define the geometry of the rectangular tank.
    P1=[0._RP,0._RP,Nzero]
    P2=[Ldom(1),0._RP,Nzero]
    P3=[Ldom(1),Ldom(2),Nzero]
    P4=[0._RP,Ldom(2),Nzero] 
    P5=[0._RP,0._RP,-Ldom(3)]
    P6=[Ldom(1),0._RP,-Ldom(3)]
    P7=[Ldom(1),Ldom(2),-Ldom(3)]
    P8=[0._RP,Ldom(2),-Ldom(3)]
    
    ! Main panel of the rectangular tank.
    Facette(:,1:3) = reshape([P2,P3,P1],(/3,3/))
    n(:,1)=[Ldom(4),Ldom(4),iabs,1.0_RP]
    typeFrontiere(1) = 0

    Facette(:,4:6) = reshape([P3,P7,P4],(/3,3/))
    n(:,2)=[Ldom(4),Ldom(4),iabsz,1.0_RP]
    typeFrontiere(2) = 1

    Facette(:,7:9) = reshape([P1,P4,P5],(/3,3/))
    n(:,3)=[Ldom(4),Ldom(4),iabs,iabsz]
    typeFrontiere(3) = 1

    Facette(:,10:12) = reshape([P2,P6,P3],(/3,3/))
    n(:,4)=[Ldom(4),Ldom(4),iabsz,iabs]
    typeFrontiere(4) = 1
    
    ! Case of the symmetry.
    if (Symmetry) then
        if (Bottom_Sym) then
            Nface = 4
        else
            Nface = 5
            Facette(:,13:15) = reshape([P8,P7,P5],(/3,3/))
            n(:,5)=[nx,ny,iabs,1._RP]
            typeFrontiere(5) = 1
        end if    
    else
        Facette(:,13:15) = reshape([P2,P1,P6],(/3,3/))
        n(:,5)=[nx,nz,iabs,iabsz]
        typeFrontiere(5) = 1
        if (Bottom_Sym) then
            Nface = 5
        else
            Nface = 6
            Facette(:,16:18) = reshape([P8,P7,P5],(/3,3/))
            n(:,6)=[nx,ny,iabs,1._RP]
            typeFrontiere(6) = 1
        end if  
    end if
    
    ! Mesh parameters
    Mesh%Nnoeud=0
    Mesh%Nfacette=0
    Mesh%DimTank = Ldom
    Mesh%TypeM = 1
    Mesh%Origine = Orig
    
    call mesh_double( Facette(:,1:3), Mesh, n(:,1), typeFrontiere(1), Orig, 0, InputData%Position(1,1,1))
    Mesh%FS%IndFS = [1, 1, Mesh%Nnoeud, Mesh%Nfacette]

    ! Mesh of the tank
    Mesh%NBody = Mesh%NBody + 1 ! At that moment Mesh%NBody = 1, the tank mesh is the first one.
    Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1] ! Tank nodes in the total mesh
    Mesh%Body(Mesh%NBody)%CMD = .false. ! By default: not moved (CMD = corps maillé et déformé)
    if (DeformMesh) Mesh%Body(Mesh%NBody)%CMD(2) = .true. ! If moving mesh: CMD = true
    Mesh%Body(Mesh%NBody)%is_tank = .true. ! Presence of the tank

Mesh%Body(Mesh%NBody)%DimBody = Ldom(1:3)
    j = 2
    call mesh_double( Facette(:,3*(j-1)+1:3*j), Mesh, n(:,j), typeFrontiere(j), Orig, j, InputData%Position(1,1,1))
    do j=3,Nface
        call mesh_double( Facette(:,3*(j-1)+1:3*j), Mesh, n(:,j), typeFrontiere(j), Orig, j)
    end do
    Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette] ! Tank panels in the total mesh
    
    ! Mesh of the body
    call MeshBody(Mesh,InputData)
    
    ! Case of the closed tank
    if (cuve_ferme) then
        Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3)
        Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4)
    else
        Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3) - Mesh%Body(1)%IndBody(3) + 1 ! The nodes of the tank are taken off.
        Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4) - Mesh%Body(1)%IndBody(4) + 1 ! The panels of the tank are taken off.
    end if
    
    ! Geometrical properties of the mesh
    call GeomInit(Mesh, fgeom_vect, 0._RP,InputData, .false.)

    print*, 'Mesh Rectangulaire cree : Nnoeud = ', Mesh%Nnoeud, 'Nfacette = ', Mesh%Nfacette
    
end subroutine bodygen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                         Création d'un maillage en coordonnées cartésiennes                       !
!                                 LETOURNEL Lucas    03/13                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mesh_double( Facette, Mesh , n, typeFrontiere, Orig, Npanneau, Pcorps)

    ! Paramètres
    real(rp), dimension(3,3), intent(in) :: Facette
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage), intent(inout) :: Mesh
    !real(rp), dimension(3,PointMax), intent(inout) :: P
    !real(rp), dimension(4,PointMax), intent(inout) :: S
    real(rp), dimension(4), intent(in) :: n
    integer, intent(in) :: typeFrontiere
    real(rp), dimension(3), intent(in), optional :: Orig
    integer, intent(in) :: NPanneau
    real(rp), intent(in), optional :: Pcorps
    !Variables Locales
    integer :: j, k, nx, ny, i1, i2, i3, i4, Nnoeud, Nfacette, NnTemp, NfTemp
    real(rp) :: AB, AC, ix, iy, dx, dy, AP, BP
    real(rp), allocatable :: L(:), LL(:), tempLL(:)
    real(rp), dimension(3) :: A, B, C, x, y

    dx=n(1) ! nombre de point de discrétisation dans la direction x
    dy=n(2) ! nombre de point de discrétisation dans la direction y
    ix=n(3) ! coefficient permettant de distendre le maillage
    iy=n(4) ! coefficient permettant de distendre le maillage
    A=Facette(:,1)
    B=Facette(:,2)
    C=Facette(:,3)
    AB=norm2(B-A)
    AC=norm2(C-A)
    x = (B-A)/AB
    y = (C-A)/AC
    nx = AB/dx
    ny = AC/dy
    Allocate(L(2*nx), LL(2*ny+1))
    L=0
    LL=0
    ! X direction
    j=2
    L(j) = dx
    do while(L(j).lt.AB)
        j=j+1
        L(j) = dx*(j-1)**ix
    end do
    L(j) = AB
    if (L(j)-L(j-1).lt.L(j-1)-L(j-2)) then
        j=j-1
        L(j) = AB
    end if
    nx = j-1
    ! Y direction
    if (present(Pcorps)) then
    AP = (A(1) + Orig(1)) - Pcorps 
    BP = AC - AP

    j=2
    LL(j) = dy
    do while(LL(j).lt.AP)
        j=j+1
        LL(j) = dy*(j-1)**iy
    end do
    LL(j) = AP
    if (LL(j)-LL(j-1).lt.LL(j-1)-LL(j-2)) then
        j=j-1
        LL(j) = AP
    end if
    ny = j
    allocate(tempLL(ny))
    tempLL(1:ny) = LL(1:ny)
    LL = 0._RP
    do j=1,ny-1
        LL(j+1) = AP-tempLL(ny-j)
    end do
    deallocate(tempLL)

    j=j+1
    LL(j) = AP + dy
    do while(LL(j).lt.AC)
        j=j+1
        LL(j) = AP + dy*(j-ny)**iy
    end do
    LL(j) = AC
    if (LL(j)-LL(j-1).lt.LL(j-1)-LL(j-2)) then
        j=j-1
        LL(j) = AC
    end if
    ny = j-1
    else
        j=2
        LL(j) = dy
        do while(LL(j).lt.AC)
            j=j+1
            LL(j) = dy*(j-1)**iy
        end do
        LL(j) = AC
        if (LL(j)-LL(j-1).lt.LL(j-1)-LL(j-2)) then
            j=j-1
            LL(j) = AC
        end if
        ny = j-1    
    end if
    
    Nnoeud=Mesh%Nnoeud
    NnTemp=Nnoeud
    do j=1,nx+1
        do k=1,ny+1 
            Nnoeud=Nnoeud+1 
            Mesh%Tnoeud(Nnoeud)%Pnoeud=A+Orig+L(j)*x+LL(k)*y
            Mesh%Tnoeud(Nnoeud)%TypeNoeud = TypeFrontiere
        end do
    end do

    Nfacette=Mesh%Nfacette-1
    NfTemp=Nfacette
    do j=1,nx
        do k=1,ny
            Nfacette=Nfacette+2
            i1= Mesh%Nnoeud+(k+(j-1)*(ny+1))
            i2=Mesh%Nnoeud+(k+(j-1)*(ny+1)+1)
            i3=Mesh%Nnoeud+(k+j*(ny+1)+1)
            i4=Mesh%Nnoeud+(k+j*(ny+1))
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1,i2,i3]
            Mesh%Tfacette(Nfacette)%TypeFrontiere=typeFrontiere
            Mesh%Tfacette(Nfacette+1)%Tnoeud=[i1,i3,i4]
            Mesh%Tfacette(Nfacette+1)%TypeFrontiere=typeFrontiere       
        end do
    end do 
    Mesh%Tnoeud(NnTemp+1:Nnoeud)%Npanneau = Npanneau
    Mesh%Tfacette(NfTemp+2:Nfacette+1)%Npanneau = Npanneau
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette = NFacette+1

end subroutine mesh_double

subroutine MeshBody(Mesh, InputData, Origine)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage)                                 :: Mesh                                         ! Mesh
    real(rp), dimension(3), optional, intent(in)    :: Origine                                      ! Origine of the inertial frame
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData                                    ! Input data
    
    integer                                         :: j, jr, jz, Ntheta, Nz, modu, itemp, Nrtemp,nc! Parameters
    real(rp)                                        :: dz, Rayon
    real(rp), dimension(:), allocatable             :: RTab, Ztab
    real(rp), dimension(3)                          :: Orig, OrigineBody
    real(rp), dimension(3,3)                        :: P1, P2, P3, Pt
    
    ! This subroutine creates the mesh of the floaters.
    
    ! Origine definition
    Orig=0._RP
    if (present(Origine)) Orig=Origine
    
    do nc=1,NBodies
                
        ! Origine definition
        OrigineBody = InputData%Position(1:3,1,nc) + Orig ! Actually Orig is always equal to [0,0,0].
        if (not(InputData%free_body(nc)).and.not(lineaireBody)) then
            do j=1,3
                if(InputData%dll_dir(j,nc))then
                    OrigineBody(j) = OrigineBody(j) + InputData%ACorps(nc) ! ???
                end if
            end do
        end if
    
        ! Mesh parameters
        Mesh%NBody = Mesh%NBody + 1 ! At that moment Mesh%NBody = 2 for the first loop step
        Mesh%Body(Mesh%NBody)%GBody(1:3) = OrigineBody ! PositionG (A changer PYW)
        Mesh%Body(Mesh%NBody)%MBody = 0._RP
        Mesh%Body(Mesh%NBody)%VBody = 0._RP
        Mesh%Body(Mesh%NBody)%ABody = 0._RP
        Mesh%Body(Mesh%NBody)%DimBody(1:2) = InputData%Lgeom(1:2,nc)
        Mesh%Body(Mesh%NBody)%CMD = .false.
        Mesh%Body(Mesh%NBody)%CMD(1) = .true. ! Body
        Mesh%Body(Mesh%NBody)%is_tank = .false.
        Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1] ! Body nodes in the total mesh
    
        ! Mesh of the body
        if (InputData%igtype(nc).eq.0) then ! Sphere
            
            print*,"This igtype value is not used anymore. Pb with Lgeom!"
            pause
            
            Rayon = InputData%Lgeom(1,nc) ! Radius
            if (.true.) then ! ! Maillage Non Polaire
                call mesh_sphere_unstruct(Mesh, Rayon, OrigineBody,InputData,nc)
            else ! Maillage Polaire : à modifier (LL)
                if (Symmetry) then
                    call mesh_spheresym(Mesh, Rayon, OrigineBody, InputData%NphiSphere(1))
                else
                    call mesh_sphere(Mesh, Rayon, OrigineBody, 2*InputData%NphiSphere(1))
                end if
            end if
        
        elseif (InputData%igtype(nc).eq.2) then ! Cylinder
        
            Rayon = InputData%Lgeom(2,nc)
            Nrtemp = InputData%NphiSphere(nc)
            Ntheta = Nth
            allocate(Rtab(Nrtemp))
            Rtab(1:Nrtemp) = (/( Rayon/float(Nrtemp-1)*(jr-1),jr=1,Nrtemp)/)
        
            ! Upper disc
            OrigineBody = InputData%Position(1:3,1,nc) + Orig + 0.5_RP*InputData%Lgeom(2,nc)*[0,0,1]
            modu = 0 !mod(modu+NphiSphere+1,2) ! On commence par le centre, d'ou le NphiSphere
            if (Symmetry) then
                call mesh_cylhorsym(Mesh, InputData%NphiSphere(nc), Ntheta, Rtab, 0._RP, OrigineBody, .true., modu)
            else
                call mesh_cylhor(Mesh, InputData%NphiSphere(nc), Ntheta, Rtab, 0._RP, OrigineBody, .true., modu)
            end if
        
            ! Body of the cylinder
            dz = Pi/InputData%NphiSphere(nc); Nz = InputData%Lgeom(2,nc)/dz; dz = -InputData%Lgeom(2,nc)/float(Nz); Nz = Nz+1
            if (mod(Nz,2).eq.0) then
                dz = -InputData%Lgeom(2,nc)/float(Nz)
                Nz=Nz+1
            end if
            allocate(Ztab(Nz))
            Ztab(1:Nz) = (/(dz*(jz-1),jz=1,Nz)/)
            modu = mod(modu+InputData%NphiSphere(nc)+1,2)
            if (Symmetry) then
                call mesh_cylvertsym(Mesh, Nz, InputData%NphiSphere(nc), Ztab, Rayon, OrigineBody, .false., modu)
            else
                call mesh_cylvert(Mesh, Nz, InputData%NphiSphere(nc), Ztab, Rayon, OrigineBody, .false., modu)
            end if
        
            ! Lower disc
            modu = mod(modu+Nz+InputData%NphiSphere(nc),2) ! On commence par le centre, d'ou le NphiSphere
            itemp = Mesh%Nnoeud+1
            OrigineBody = InputData%Position(1:3,1,nc) + Orig - 0.5_RP*InputData%Lgeom(2,nc)*[0._RP,0._RP,1._RP]
            if (Symmetry) then
                call mesh_cylhorsym(Mesh, InputData%NphiSphere(nc), Ntheta, Rtab, 0._RP, OrigineBody, .true., modu)
            else
                call mesh_cylhor(Mesh, InputData%NphiSphere(nc), Ntheta, Rtab, 0._RP, OrigineBody, .true., modu)
            end if
        
            ! Updating the position of the mesh.
            P1 = reshape([-1._RP,0._RP,0._RP,0._RP,1._RP,0._RP,0._RP,0._RP,-1._RP],(/3,3/))
            do j=itemp,Mesh%Nnoeud
                Mesh%Tnoeud(j)%Pnoeud = OrigineBody + matmul(P1,Mesh%Tnoeud(j)%Pnoeud-OrigineBody)
            end do
            deallocate(Rtab,Ztab)
        
        elseif (InputData%igtype(nc).eq.1) then ! Cube
            call Mesh_Body_Cube(Mesh,InputData)
        elseif(InputData%igtype(nc).eq.3) then
            call extract_mesh2bis(Mesh,InputData,nc)
        end if
        
        Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = 1
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = 1
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 2
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 2
        
        ! Updating the final position of the mesh.
        if (norm2(InputData%Position(1:3,2,nc)).gt.Epsilon) then
            P1 = reshape([1._RP,0._RP,0._RP,0._RP,cos(InputData%Position(1,2,nc)),-sin(InputData%Position(1,2,nc)),0._RP,sin(InputData%Position(1,2,nc)),cos(InputData%Position(1,2,nc))],(/3,3/))
            P2 = reshape([cos(InputData%Position(2,2,nc)),0._RP,-sin(InputData%Position(2,2,nc)),0._RP,1._RP,0._RP,sin(InputData%Position(2,2,nc)),0._RP,cos(InputData%Position(2,2,nc))],(/3,3/))
            P3 = reshape([cos(InputData%Position(3,2,nc)),sin(InputData%Position(3,2,nc)),0._RP,-sin(InputData%Position(3,2,nc)),cos(InputData%Position(3,2,nc)),0._RP,0._RP,0._RP,1._RP],(/3,3/))
            Pt = matmul(P1,P2)
            Pt = matmul(Pt,P3)
            do j=Mesh%Body(Mesh%NBody)%IndBody(1),Mesh%Body(Mesh%NBody)%IndBody(3)
                Mesh%Tnoeud(j)%Pnoeud = Mesh%Body(Mesh%NBody)%GBody(1:3) + matmul(Pt,Mesh%Tnoeud(j)%Pnoeud-Mesh%Body(Mesh%NBody)%GBody(1:3))
            end do
        end if
        
        ! New initialization for the next body
        OrigineBody = 0._RP
        
    end do
    
end subroutine MeshBody

subroutine Mesh_Body_Cube(Mesh,InputData)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh         ! Mesh
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data
    
    integer :: Nface
    integer :: j
    real(rp), dimension(3):: Orig, P1, P2, P3, P4, P5, P6, P7, P8
    real(rp) :: iabs, nx, ny, nz, Nzero, iabsz, Lcube, dcube
    real(rp), dimension(3,18) :: Facette
    real(rp), dimension(4,6) :: n
    logical, dimension(2,6) :: cote2
    integer, dimension(6) :: typeFrontiere
    
    ! This subroutine creates the mesh of a cube.
    
    Nzero=0._RP
    Lcube = InputData%Lgeom(1,1)
    dcube = InputData%Lgeom(2,1)
    nx=dcube
    ny=nx
    nz=nx
    if (Symmetry) then
        Orig=InputData%Position(1:3,1,1) + [-Lcube/2._RP,0._RP,Lcube/2._RP] + InputData%ACorps(1)*[0._RP,0._RP,-1._RP]
    else
        Orig=InputData%Position(1:3,1,1) + [-Lcube/2._RP,-Lcube/2._RP,Lcube/2._RP] + InputData%ACorps(1)*[0._RP,0._RP,-1._RP]
    end if
    iabs = 0.75_RP
    iabsz = 1._RP
    P1=[0._RP,0._RP,Nzero]
    P2=[Lcube,0._RP,Nzero]
    P3=[Lcube,Lcube,Nzero]
    P4=[0._RP,Lcube,Nzero]
    P5=[0._RP,0._RP,-Lcube]
    P6=[Lcube,0._RP,-Lcube]
    P7=[Lcube,Lcube,-Lcube]
    P8=[0._RP,Lcube,-Lcube]
    if (Symmetry) then
        P3(2) = P3(2)/2._RP
        P4(2) = P4(2)/2._RP
        P7(2) = P7(2)/2._RP
        P8(2) = P8(2)/2._RP
        ny = float(int(ny/2))
    end if

    Facette(:,1:3) = reshape([P1,P2,P4],(/3,3/))
    n(:,1)=[nx,ny,iabs,iabs]
    cote2(:,1) = [.true.,.true.]
    typeFrontiere(1) = 1

    Facette(:,4:6) = reshape([P4,P3,P8],(/3,3/))
    n(:,2)=[nx,nz,iabs,iabs]
    cote2(:,2) = [.true.,.true.]
    typeFrontiere(2) = 1

    Facette(:,7:9) = reshape([P1,P4,P5],(/3,3/))
    n(:,3)=[ny,nz,iabs,iabs]
    cote2(:,3) = [.true.,.true.]
    typeFrontiere(3) = 1

    Facette(:,10:12) = reshape([P6,P7,P2],(/3,3/))
    n(:,4)=[ny,nz,iabs,iabs]
    cote2(:,4) = [.true.,.true.]
    typeFrontiere(4) = 1

    if (Symmetry) then
        cote2(:,1) = [.true.,.false.]
        cote2(:,3) = [.false.,.true.]
        cote2(:,4) = [.false.,.true.]
        Nface = 5
        Facette(:,13:15) = reshape([P6,P5,P7],(/3,3/))
        n(:,5)=[nx,ny,iabs,iabs]
        cote2(:,5) = [.true.,.false.]
        typeFrontiere(5) = 1
    else
        Facette(:,13:15) = reshape([P2,P1,P6],(/3,3/))
        n(:,5)=[nx,nz,iabs,iabsz]
        cote2(:,1) = [.true.,.true.]
        typeFrontiere(5) = 1
        Nface = 6
        Facette(:,16:18) = reshape([P8,P7,P5],(/3,3/))
        n(:,6)=[nx,ny,iabs,iabs]
        cote2(:,1) = [.true.,.true.]
        typeFrontiere(6) = 1 
    end if
    do j=1,Nface
        call mesh_doubleSpe( Facette(:,3*(j-1)+1:3*j), Mesh, n(:,j), cote2(:,j), typeFrontiere(j), Orig, j)
    end do
    
end subroutine Mesh_Body_Cube

subroutine mesh_doubleSpe( Facette, Mesh , n, cote2, typeFrontiere, Orig, Npanneau )

    ! Paramètres
    real(rp), dimension(3,3), intent(in) :: Facette
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage), intent(inout) :: Mesh
    !real(rp), dimension(3,PointMax), intent(inout) :: P
    !real(rp), dimension(4,PointMax), intent(inout) :: S
    real(rp), dimension(4), intent(in) :: n
    logical, dimension(2) :: cote2
    integer, intent(in) :: typeFrontiere
    real(rp), dimension(3), intent(in), optional :: Orig
    integer, intent(in) :: NPanneau
    !Variables Locales
    integer :: j, k, nx, ny, i1, i2, i3, i4, Nnoeud, Nfacette, nab, NnTemp, NfTemp
    real(rp) :: AB, AC, ix, iy, dx, dy, Lr
    real(rp), allocatable :: L(:), LL(:), Ltemp(:)
    real(rp), dimension(3) :: A, B, C, x, y

    nx=n(1) ! nombre de point de discrétisation dans la direction x
    ny=n(2) ! nombre de point de discrétisation dans la direction y
    ix=n(3) ! coefficient permettant de distendre le maillage
    iy=n(4) ! coefficient permettant de distendre le maillage
    A=Facette(:,1)
    B=Facette(:,2)
    C=Facette(:,3)
    AB=norm2(B-A)
    AC=norm2(C-A)
    x = (B-A)/AB
    y = (C-A)/AC
    dx = AB/real(nx)
    dy = AC/real(ny)
    Allocate(L(2*nx), LL(2*ny))
    L=0
    LL=0
    ! Selon x
    if (abs(ix-1._RP).lt.Epsilon2) then
        do j=0,nx
            L(j+1) = AB*(j/real(nx))
        end do
    else
        Lr = AB*1._RP/3._RP
        nab=(AB-Lr)/dx
        do j=0,nab
            L(j+1) = (AB-Lr)*(j/real(nab))
        end do
        j=nab
        do while(L(j+1).lt.AB)
            j=j+1
            if (abs((j-nab)**ix-1/10._RP).lt.Epsilon) then
                L(j+1) = (AB-Lr) + dx/10._RP
            else
                L(j+1) = (AB-Lr) + dx*(j-nab)**ix
            end if
        end do
        L(j) = AB
        nx=j
        Allocate(Ltemp(nx))
        if (cote2(1)) then
            Ltemp(1:nx) = L(1:nx)/2._RP
            L = 0._RP
            if ( modulo(nx,2).lt.Epsilon) then
                L(nx/2+1:nx:1) = AB/2._RP + Ltemp(2:nx:2)
                L(1:nx/2:1) = AB/2._RP - Ltemp(nx:2:-2)            
            else
                L(nx/2+1:nx:1) = AB/2._RP + Ltemp(1:nx:2)
                L(1:nx/2:1) = AB/2._RP - Ltemp(nx:3:-2)
            end if
        end if
        deallocate(Ltemp)
        nx=nx-1
    end if
    ! Selon y
    if (abs(iy-1._RP).lt.Epsilon2) then
        do j=0,ny
            LL(j+1) = AC*(j/real(ny))
        end do
    else
        Lr = AC*1._RP/3._RP
        nab=(AC-Lr)/dy
        do j=0,nab
            LL(j+1) = (AC-Lr)*(j/real(nab))
        end do
        j=nab
        do while(LL(j+1).lt.AC)
            j=j+1
            if (abs((j-nab)**iy-1/10._RP).lt.Epsilon) then
                LL(j+1) = (AC-Lr) + dy/10._RP
            else
                LL(j+1) = (AC-Lr) + dy*(j-nab)**iy
            end if
        end do
        LL(j) = AC
        ny=j
        Allocate(Ltemp(ny))
        if (cote2(2)) then
            Ltemp(1:ny) = LL(1:ny)/2._RP
            LL = 0._RP
            if ( modulo(ny,2).lt.Epsilon) then
                LL(ny/2+1:ny:1) = AC/2._RP + Ltemp(2:ny:2)
                LL(1:ny/2:1) = AC/2._RP - Ltemp(ny:2:-2)            
            else
                LL(ny/2+1:ny:1) = AC/2._RP + Ltemp(1:ny:2)
                LL(1:ny/2:1) = AC/2._RP - Ltemp(ny:3:-2)
            end if
        end if
        deallocate(Ltemp)
        ny=ny-1
    end if

    Nnoeud=Mesh%Nnoeud
    NnTemp=Nnoeud
    do j=1,nx+1
        do k=1,ny+1 
            Nnoeud=Nnoeud+1 
            Mesh%Tnoeud(Nnoeud)%Pnoeud=A+Orig+L(j)*x+LL(k)*y
            Mesh%Tnoeud(Nnoeud)%TypeNoeud = TypeFrontiere
        end do
    end do

    Nfacette=Mesh%Nfacette-1
    NfTemp=Nfacette
    do j=1,nx
        do k=1,ny
            Nfacette=Nfacette+2
            i1= Mesh%Nnoeud+(k+(j-1)*(ny+1))
            i2=Mesh%Nnoeud+(k+(j-1)*(ny+1)+1)
            i3=Mesh%Nnoeud+(k+j*(ny+1)+1)
            i4=Mesh%Nnoeud+(k+j*(ny+1))
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1,i2,i3]
            Mesh%Tfacette(Nfacette)%TypeFrontiere=typeFrontiere
            Mesh%Tfacette(Nfacette+1)%Tnoeud=[i1,i3,i4]
            Mesh%Tfacette(Nfacette+1)%TypeFrontiere=typeFrontiere       
        end do
    end do 
    Mesh%Tnoeud(NnTemp+1:Nnoeud)%Npanneau = Npanneau
    Mesh%Tfacette(NfTemp+2:Nfacette+1)%Npanneau = Npanneau
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette = NFacette+1

end subroutine mesh_doubleSPe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                         Création d'un maillage en coordonnées cylindriques                       !
!                                 LETOURNEL Lucas    03/13                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mesh_cyl( Mesh, fgeom_vect,InputData, Origine)

    
    !f2py integer*1, dimension(1000)                        :: Mesh
    type(TMaillage), intent(inout)                          :: Mesh                                         ! Mesh
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                                           ! Geometries
    !f2py integer*1, dimension(1000)                        :: InputData
    type(InputDataStruct),intent(inout)                        :: InputData                                    ! Input data
    real(rp), dimension(3), optional, intent(in)            :: Origine                                      ! Origine
    
    integer                                                 :: j, jr, jz, Ntheta, Nz, modu, TypeN, Nrtemp   !
    real(rp)                                                :: dr, dtheta, dz, Test                         !
    real(rp), dimension(PointMax)                           :: Ztab,RtabSL                                  !
    real(rp), dimension(3)                                  :: Orig                                         ! Origine
    real(rp)                                                :: Cramp_1,Cramp_2,Cramp_3,Cramp_4              ! Ramps

    ! This subroutine creates a mesh from a cylindrical mesh of the tank.
    
    ! Origin
    Orig=0
    if (present(Origine)) Orig=Origine
    
    ! Mesh parameters
    Mesh%TypeM = 0          ! Circular Mesh
    Mesh%DimTank = LDom
    Mesh%Origine = Orig     ! Origine du maillage

    modu = 0
    if (Symmetry) then
        Ntheta = Nth
    else
        Ntheta = 2*(Nth-1)
    end if
    Mesh%Nnoeud = 0
    Mesh%Nfacette = 0

    ! Mesh of the Free Surface (upper disc)
    if (.true.) then    
        dTheta = 2._RP*PI/(Ntheta)
        if (Ldom(2)/Nr - Ldom(5) .lt.Epsilon) then
            dr = Ldom(2)/Nr
        else
            dr = Ldom(5)
        end if
        RtabSL(1) = 0
        RtabSL(2) = dr
    !RtabSL(1) = Lgeom(1)
    !RtabSL(2) = Lgeom(1) + dr
        j = 2
        Nrtemp = 2
        do while (RtabSL(j).lt.Mesh%DimTank(2))
            j=j+1
            call Computation_Cramp(1._RP*float(j),1._RP*float(Nrtemp+5),1._RP,Cramp_1)
            call Computation_Cramp(1._RP*float(j),1._RP*float(Nrtemp+5),1._RP,Cramp_2)
            Test = Cramp_1*0.35_RP*dTheta*RtabSL(j-1)/2._RP + (1-Cramp_2)*dr
            if (Test.lt.abs(RtabSL(j-1)-RtabSL(j-2))) then ! Permet d'avoir des mailles de taille régulière au début de la plage absorbante
                RtabSL(j) = RtabSL(j-1) + abs(RtabSL(j-1)-RtabSL(j-2))
            else
                RtabSL(j) = RtabSL(j-1) + Test
            end if
        end do
        Nrtemp=j
        Mesh%DimTank(2) = RtabSL(Nrtemp)
    
        do while (RtabSL(j).lt.Mesh%DimTank(1))
            j=j+1
            call Computation_Cramp(1._RP*float(j),1._RP*float(Nrtemp+5),1._RP,Cramp_3)
            call Computation_Cramp(1._RP*float(j),1._RP*float(Nrtemp+5),1._RP,Cramp_4)
            Test = Cramp_3*0.75_RP*dTheta*(RtabSL(j-1)-1._RP*Mesh%DimTank(2))/2._RP + (1-Cramp_4)*dr
            if (Test.lt.abs(RtabSL(j-1)-RtabSL(j-2))) then ! Permet d'avoir des mailles de taille régulière au début de la plage absorbante
                RtabSL(j) = RtabSL(j-1) + abs(RtabSL(j-1)-RtabSL(j-2))
            else
                RtabSL(j) = RtabSL(j-1) + Test
            end if
        end do
        Nrtemp=j
        Mesh%DimTank(1) = RtabSL(Nrtemp)
        Mesh%FS%DimFS(1) = RtabSL(Nrtemp)
        Mesh%FS%IndFS(1:2) = [1, 1]
        if (Symmetry) then
            call mesh_cylhorsym(Mesh, Nrtemp, Ntheta, RtabSL, 0._RP, Orig, .true., modu)
        else
            call mesh_cylhor(Mesh, Nrtemp, Ntheta, RtabSL, 0._RP, Orig, .true., modu)
        end if
    else
        Mesh%FS%DimFS = Mesh%DimTank
        Mesh%FS%IndFS(1:2) = [1, 1]
        call mesh_FS_unstruc(Mesh, Mesh%DimTank(1), Orig)
    end if

    Mesh%FS%IndFS(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%typeNoeud = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%typeFrontiere = 0
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Npanneau = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%Npanneau = 0
    
    ! Body of the cylinder
    if (cuve_ferme) then
        dz = RtabSL(Nrtemp)-RtabSL(Nrtemp-1)
        Nz = Mesh%DimTank(3)/dz
        TypeN = 1
    else
        Nz = 1
        TypeN = 10
    end if
    dz = -Mesh%DimTank(3)/float(Nz)
    Nz = Nz+1
    if (Nz.eq.1) then
        Nz = 2
        dz = -Mesh%DimTank(3)
    end if
    Ztab(1:Nz) = (/(dz*(jz-1),jz=1,Nz)/)
    
    if(Symmetry)then
        modu =  mod(Nrtemp+1,2)
    else
        modu =  mod(Nrtemp,2)
    end if
    Mesh%NBody = Mesh%NBody + 1 ! At that moment Mesh%NBody = 1
    Mesh%Body(Mesh%NBody)%GBody(1:3) = Ldom(1)*100._RP*[0._RP,1._RP,0._RP]
    Mesh%Body(Mesh%NBody)%MBody = 0._RP
    Mesh%Body(Mesh%NBody)%VBody = 0._RP
    Mesh%Body(Mesh%NBody)%ABody = 0._RP
    Mesh%Body(Mesh%NBody)%DimBody = Ldom(1:3)
    Mesh%Body(Mesh%NBody)%CMD = .false. ! Tank
    
    if(DeformMesh) Mesh%Body(Mesh%NBody)%CMD(2) = .true.
    Mesh%Body(Mesh%NBody)%is_tank = .true.
    Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
    
    if (Symmetry) then
        call mesh_cylvertsym(Mesh, Nz, Ntheta, Ztab, Mesh%DimTank(1), Orig, .true., modu)
    else
        call mesh_cylvert(Mesh, Nz, Ntheta, Ztab, Mesh%DimTank(1), Orig, .true., modu)
    end if
    
    Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = TypeN
    Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = TypeN
    Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 1 
    Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 1
        
    ! Mesh of the bottom
    if (bottom_sym) then
        ! Nothing
    else
        if (cuve_ferme) then
            Nrtemp=Nrtemp+1
        else
            Nrtemp=2
            RtabSL(1:Nrtemp)=(/( Mesh%DimTank(1)/float(Nrtemp-1)*(jr-1),jr=1,Nrtemp)/)
        end if
        modu = mod(modu+Nrtemp+Nz,2) ! On commence par le centre, d'ou le Nrtemp
        
        Mesh%NBody = Mesh%NBody + 1 ! At that moment Mesh%NBody = 2
        
        Mesh%Body(Mesh%NBody)%GBody(1:3) = -1000._RP*[0._RP,0._RP,1._RP] 
        Mesh%Body(Mesh%NBody)%MBody = 0._RP
        Mesh%Body(Mesh%NBody)%VBody = 0._RP
        Mesh%Body(Mesh%NBody)%ABody = 0._RP
        Mesh%Body(Mesh%NBody)%DimBody = LDom(1:3)
        Mesh%Body(Mesh%NBody)%CMD = .false.
        Mesh%Body(Mesh%NBody)%is_tank = .true.
        Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
        if (Symmetry) then
            call mesh_cylhorsym(Mesh, Nrtemp, Ntheta, RtabSL, -Mesh%DimTank(3), Orig, .false., modu)    
        else
            call mesh_cylhor(Mesh, Nrtemp, Ntheta, RtabSL, -Mesh%DimTank(3), Orig, .false., modu)
        end if
        Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = TypeN
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = TypeN
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 3
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 3
    end if
    
    ! Mesh of the body
    call MeshBody(Mesh,InputData, Origine)
    
    ! Case of the closed tank
    if (cuve_ferme) then
        Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3)
        Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4)
    else
        if (is_body.and.Int_Body.ne.0) then
            Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3) - Mesh%Body(1)%IndBody(3) + 1 ! The nodes of the tank are taken off.
            Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4) - Mesh%Body(1)%IndBody(4) + 1 ! The panels of the tank are taken off4)
        else
            Mesh%Nsys = Mesh%FS%IndFS(3)
            Mesh%Nfsys = Mesh%FS%IndFS(4)
        end if
    end if
    
    ! Geometrical properties of the mesh
    call GeomInit(Mesh, fgeom_vect, 0._RP,InputData, .false.)
    
end subroutine mesh_cyl

subroutine mesh_cylvert(Mesh, Nz, Ntheta, Ztab, R, Orig, Normale, u)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    integer :: Nz, Ntheta
    real(rp), intent(in) :: Ztab(Nz), R, Orig(3)
    logical, intent(in) :: Normale ! définit le sens de la normale
    integer, intent(in), optional :: u
    
    integer :: jz, jtheta, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu
    real(rp) :: X, Y, Z, dtheta
    
    ! This subroutine creates the mesh of the body of a cylinder without symmetry.
    
    dtheta = 2._RP*PI/float(Ntheta)
    Nnoeud=Mesh%Nnoeud
    Nntemp=Nnoeud
    Nfacette=Mesh%Nfacette
    modu = 0
    if (present(u)) modu = u
    do jz=1,Nz
        do jtheta=1,Ntheta 
            Nnoeud=Nnoeud+1
            X = R*COS((jtheta-1)*dtheta+0.5_rp*dtheta*MOD(decalage*jz+modu,2))
            Y = R*SIN((jtheta-1)*dtheta+0.5_rp*dtheta*MOD(decalage*jz+modu,2))
            Z = Ztab(jz)
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z]+Orig
        end do
    end do
    do jz=1,Nz-1
        do jtheta=1,Ntheta-1
            i1 = NnTemp + jtheta+(jz-1)*Ntheta
            if (Normale) then
                i2 = NnTemp + jtheta+jz*Ntheta
                i3 = NnTemp + jtheta+1+(jz-1)*Ntheta
            else
                i2 = NnTemp + jtheta+1+(jz-1)*Ntheta
                i3 = NnTemp + jtheta+jz*Ntheta
            end if
            i4 = NnTemp + jtheta+1+jz*Ntheta        
	        if (mod(decalage*jz+modu,2).eq.1) then
	            Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]
		        Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
	        else
	            Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
		        Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	        end if
        end do
        i1=NnTemp+(jz*(Ntheta))
        if (Normale) then
            i2=NnTemp+((jz+1)*(Ntheta))
            i3=NnTemp+((jz-1)*(Ntheta)+1)
        else
            i2=NnTemp+((jz-1)*(Ntheta)+1)
            i3=NnTemp+((jz+1)*(Ntheta))
        end if    
        i4=NnTemp+((jz)*(Ntheta)+1)
        if (mod(decalage*jz+modu,2).eq.1) then
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]
	        Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
        else
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
	        Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
        end if
    end do
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette = Nfacette
    
end subroutine mesh_cylvert

subroutine mesh_cylhor(Mesh, Nr, Ntheta, Rtab, Z, Orig, Normale, u)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    integer :: Nr, Ntheta
    real(rp), intent(in) :: RTab(Nr), Z, Orig(3)
    logical, intent(in) :: Normale ! définit le sens de la normale
    integer, intent(in), optional :: u
    
    integer :: jr, jtheta, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu, modv
    real(rp) :: X, Y, dtheta, ThetaTab(Ntheta)
    
    ! This subroutine creates a mesh of a disc with symmetry.
    
    dtheta = 2._RP*PI/float(Ntheta)
    ThetaTab(1:Ntheta) = (/(jtheta*dtheta,jtheta=0,Ntheta-1)/)
    Nnoeud=Mesh%Nnoeud
    Nntemp=Nnoeud
    NFacette=Mesh%NFacette
    
    modu = 0
    if (present(u)) modu=u
    modv = modu
    if (RTab(1).lt.Epsilon2) then
        Nnoeud=Nnoeud+1
        Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, Z]+Orig
        do jtheta=1,Ntheta-1
            Nfacette=Nfacette+1
            if (Normale) then
                Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jtheta+2, jtheta+1]
            else
                Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jtheta+1, jtheta+2]
            end if
        end do
        Nfacette=Nfacette+1
        if (Normale) then
            Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, 2, Ntheta+1]
        else
            Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, Ntheta+1, 2]
        end if
    end if
    NnTemp = Nnoeud
    do jr=1,Nr
        if (RTab(jr).gt.Epsilon2) then
            do jtheta=1,Ntheta
	            Nnoeud=Nnoeud+1
	            X = RTab(jr)*COS(ThetaTab(jtheta)+0.5_rp*dtheta*MOD(decalage*jr+modu,2))
		        Y = RTab(jr)*SIN(ThetaTab(jtheta)+0.5_rp*dtheta*MOD(decalage*jr+modu,2))
	            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
            end do
        end if
    end do
    if (RTab(1).lt.Epsilon2) then
        Nr=Nr-1
        modv = modu+1
    end if
    do jr=1,Nr-1
        do jtheta=1,Ntheta-1
            i1 = NnTemp + jtheta+(jr-1)*Ntheta
            if (Normale) then
	            i2 = NnTemp + jtheta+jr*Ntheta
	            i3 = NnTemp + jtheta+1+(jr-1)*Ntheta
	        else
	            i2 = NnTemp + jtheta+1+(jr-1)*Ntheta
	            i3 = NnTemp + jtheta+jr*Ntheta
	        end if
	        i4 = NnTemp + jtheta+1+jr*Ntheta
	        if (MOD(decalage*jr+modv,2).EQ.1) THEN
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]
		        Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
	        ELSE
	            Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
			    Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	        end if
        end do
        i1 = NnTemp + jr*Ntheta
        if (Normale) then
            i2 = NnTemp + (jr+1)*Ntheta
            i3 = NnTemp + (jr-1)*Ntheta+1
        else
            i2 = NnTemp + (jr-1)*Ntheta+1
            i3 = NnTemp + (jr+1)*Ntheta
        end if
        i4 = NnTemp + jr*Ntheta+1
        if (MOD(decalage*jr+modv,2).EQ.1) THEN
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]	
	        Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
        ELSE
	        Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
	        Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
        end if
    end do
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette= Nfacette
    
end subroutine mesh_cylhor

subroutine mesh_sphere(Mesh, Rsphere, OrigineSphere, NthS, demi)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    real(rp), intent(in) :: Rsphere
    real(rp), dimension(3) :: OrigineSphere
    integer, intent(in) :: NthS
    logical, optional :: demi
    
    integer :: jpsi, jphi, Nphi, Npsi, Nfacette, Nnoeud, Nntemp, i1, i2, i3, i4
    real(rp) :: dphi, dPsi, X, Y, Z
    real(rp), allocatable :: phiTab(:), PsiTab(:)
    
    ! This subroutine creates a polar mesh of a sphere without symmetry.
    
    Nphi = NthS
    Npsi = Nphi/2
    dphi=2._RP*PI/(Nphi)
    dPsi = dphi
    
    if (present(demi)) then
        if (demi) then
            dPsi = dPsi/2._RP
            Npsi = Npsi/2._RP
        end if
    end if
    allocate(phiTab(1:Nphi), PsiTab(1:Npsi+1))
    phiTab(1:Nphi) = (/(jphi*dphi,jphi=0,Nphi-1)/)
    PsiTab(1:Npsi+1) = (/(Pi/2._RP-jpsi*dPsi,jpsi=0,Npsi)/)
    Nnoeud = Mesh%Nnoeud
    Nfacette = Mesh%Nfacette
    Nntemp = Nnoeud
    do jpsi=1,Npsi+1
        if (jpsi.eq.1) then
            Nnoeud=Nnoeud+1
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, RSphere] + OrigineSphere
            do jphi=1,Nphi-1
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jphi+1, jphi+2]
            end do
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, Nphi+1, 2]
         elseif (jpsi.ne.1 .and. jpsi.ne.(Npsi+1)) then
            do jphi=1,Nphi
                Nnoeud=Nnoeud+1
                X = Rsphere*cos(phiTab(jphi)+0.5_rp*dphi*mod(0*jpsi,2))*cos(PsiTab(jpsi))
	            Y = Rsphere*sin(phiTab(jphi)+0.5_rp*dphi*mod(0*jpsi,2))*cos(PsiTab(jpsi))
                Z = Rsphere*sin(PsiTab(jpsi))
                Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + OrigineSphere
           end do       
        elseif (jpsi.eq.Npsi+1) then
            Nnoeud=Nnoeud+1
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, -RSphere] + OrigineSphere
            do jphi=1,Nphi-1
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=Nntemp + 1 + (Npsi-2)*Nphi + [jphi+1, jphi, Nphi+1]
            end do
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=Nntemp + 1 + (Npsi-2)*Nphi + [1, Nphi, Nphi+1]
        end if
    end do
    do jpsi=2,Npsi-1  
       do jphi=1,Nphi-1  
            i1 = NnTemp + 1 + (jpsi-2)*Nphi + jphi
            i2 = NnTemp + 1 + (jpsi-2)*Nphi + jphi + 1
            i3 = NnTemp + 1 + (jpsi-1)*Nphi + jphi
            i4 = NnTemp + 1 + (jpsi-1)*Nphi + jphi + 1
            if (mod(0*jpsi,2).eq.1) then
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]	
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
            else
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]	
                Nfacette=Nfacette+1
                Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
            end if
        end do
        i1 = NnTemp + 1 + (jpsi-1)*Nphi
        i2 = NnTemp + 1 + (jpsi-2)*Nphi + 1
        i3 = NnTemp + 1 + jpsi*Nphi
        i4 = NnTemp + 1 + (jpsi-1)*Nphi + 1 
        if (mod(0*jpsi,2).eq.1) then
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]	
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
        else
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]	
            Nfacette=Nfacette+1
            Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
        end if
    end do
    Mesh%Nnoeud=Nnoeud
    Mesh%Nfacette=Nfacette
    
end subroutine mesh_sphere

subroutine mesh_cylvertsym(Mesh, Nz, Ntheta, Ztab, R, Orig, Normale, u)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    integer :: Nz, Ntheta
    real(rp), intent(in) :: Ztab(Nz), R, Orig(3)
    logical, intent(in) :: Normale ! définit le sens de la normale
    integer, intent(in), optional :: u
    
    integer :: jz, jtheta, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu
    real(rp) :: X, Y, Z, dtheta
    
    ! This subroutine creates the mesh of the body of a cylinder with symmetry.
    
    dtheta = PI/float(Ntheta-1)
    Nnoeud=Mesh%Nnoeud
    Nntemp=Nnoeud
    Nfacette=Mesh%Nfacette
    modu = 0
    if (present(u)) modu = u
    do jz=1,Nz
        if (decalage.eq.0 .or. MOD(jz+modu,2).eq.0) then
            do jtheta=1,Ntheta
                Nnoeud=Nnoeud+1
                X = R*COS((jtheta-1)*dtheta)
                Y = R*SIN((jtheta-1)*dtheta)
                Z = Ztab(jz)
                Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
            end do
        else
            Nnoeud=Nnoeud+1
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [R, 0._RP, Ztab(jz)] + Orig
            do jtheta=1,Ntheta-1
                Nnoeud=Nnoeud+1
                X = R*COS((jtheta-1)*dtheta+0.5_rp*dtheta*MOD(decalage*jz+modu,2))
                Y = R*SIN((jtheta-1)*dtheta+0.5_rp*dtheta*MOD(decalage*jz+modu,2))
                Z = Ztab(jz)
                Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
            end do
            Nnoeud=Nnoeud+1
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [-R, 0._RP, Ztab(jz)] + Orig
        end if
    end do
    do jz=1,Nz-1
        if (decalage.eq.0) then
            do jtheta=1,Ntheta-1
                i1 = NnTemp + jtheta+(jz-1)*Ntheta
                if (Normale) then
                    i2 = NnTemp + jtheta+jz*Ntheta
                    i3 = NnTemp + jtheta+1+(jz-1)*Ntheta
                else
                    i2 = NnTemp + jtheta+1+(jz-1)*Ntheta
                    i3 = NnTemp + jtheta+jz*Ntheta
                end if
                i4 = NnTemp + jtheta+1+jz*Ntheta        
	            if (mod(0*jz+modu,2).eq.1) then
	                Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]
		            Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
	            else
	                Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
		            Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	            end if
            end do
        else
            if (mod(decalage*jz+modu,2).EQ.0) then
                i1 = NnTemp + (jz-1)*(Ntheta+1) + 1 - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
                i2 = NnTemp + jz*(Ntheta+1) - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
                i3 = NnTemp + jz*(Ntheta+1) + 1 - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
            else
                i1 = NnTemp + jz*(Ntheta+1) - jz/2
                i2 = NnTemp + jz*(Ntheta+1) - 1 - jz/2
                i3 = NnTemp + (jz+1)*(Ntheta+1) - 1 - jz/2
            end if
            Nfacette=Nfacette+1        
            if (Normale) then
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
            else
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i2, i3]
            end if
            do jtheta=1,Ntheta-1
                i1 = NnTemp + jtheta + (jz-1)*(Ntheta+1) - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
                i2 = NnTemp + jtheta + jz*(Ntheta+1) - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
                i3 = NnTemp + jtheta + jz*(Ntheta+1) + 1 - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
                i4 = NnTemp + jtheta + (jz-1)*(Ntheta+1) + 1 - jz/2 + MOD(modu+1,2)*mod(decalage*(jz+1)+modu,2)
	            if (Normale) then
                    Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
		            Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i4, i3, i2]
	            else
	                Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i2, i4]
			        Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	            end if
            end do
        end if
    end do
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette= Nfacette
    
end subroutine mesh_cylvertsym

subroutine mesh_cylhorsym(Mesh, Nra, Ntheta, Rtab, Z, Orig, Normale, u)
    ! Parameters
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    integer, intent(in) :: Nra, Ntheta
    real(rp), intent(in) :: RTab(Nra), Z, Orig(3)
    logical, intent(in) :: Normale ! définit le sens de la normale
    integer, intent(in), optional :: u
    ! Local
    integer :: jr, jtheta, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu, modv, Nrtemp
    real(rp) :: X, Y, dtheta, ThetaTab(Ntheta)
    
    ! This subroutine creates the mesh of a disc without symmetry.
    
    dtheta = PI/float(Ntheta-1)
    ThetaTab(1:Ntheta) = (/(jtheta*dtheta,jtheta=0,Ntheta-1)/)
    Nnoeud=Mesh%Nnoeud
    Nntemp=Nnoeud
    NFacette=Mesh%NFacette
    modu = 0
    if (present(u)) modu=u
    modv = modu
    if (RTab(1).lt.Epsilon2) then
        Nnoeud=Nnoeud+1
        Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, Z]+Orig
        do jtheta=1,Ntheta-1+mod(modu,2)
            Nfacette=Nfacette+1
            if (Normale) then
                Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jtheta+2, jtheta+1]
            else
                Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jtheta+1, jtheta+2]
            end if
        end do
    end if
    NnTemp = Nnoeud
    do jr=1,Nra
        if (RTab(jr).gt.Epsilon2) then
            if (decalage.eq.0 .or. MOD(jr+modu,2).eq.0) then
                do jtheta=1,Ntheta
	                Nnoeud=Nnoeud+1
	                X = RTab(jr)*COS(ThetaTab(jtheta))
		            Y = RTab(jr)*SIN(ThetaTab(jtheta))
	                Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
                end do
            else
                Nnoeud=Nnoeud+1
                Mesh%Tnoeud(Nnoeud)%Pnoeud = [RTab(jr), 0._RP, Z] + Orig
                do jtheta=1,Ntheta-1
	                Nnoeud=Nnoeud+1
	                X = RTab(jr)*COS(ThetaTab(jtheta)+0.5_rp*dtheta*MOD(jr+modu,2))
		            Y = RTab(jr)*SIN(ThetaTab(jtheta)+0.5_rp*dtheta*MOD(jr+modu,2))
	                Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
                end do
                Nnoeud=Nnoeud+1
                Mesh%Tnoeud(Nnoeud)%Pnoeud = [-RTab(jr), 0._RP, Z] + Orig
            end if
        end if
    end do
    if (RTab(1).lt.Epsilon2) then
        Nrtemp=Nra-1
        modv = modu+1
    end if
    do jr=1,Nrtemp-1
        if (decalage.eq.0) then
            do jtheta=1,Ntheta-1
                i1 = NnTemp + jtheta+(jr-1)*Ntheta
                if (Normale) then
	                i2 = NnTemp + jtheta+jr*Ntheta
	                i3 = NnTemp + jtheta+1+(jr-1)*Ntheta
	            else
	                i2 = NnTemp + jtheta+1+(jr-1)*Ntheta
	                i3 = NnTemp + jtheta+jr*Ntheta
	            end if
	            i4 = NnTemp + jtheta+1+jr*Ntheta
	            if (mod(0*jr+modv,2).EQ.1) then
                    Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i4]
		            Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
	            else
	                Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
			        Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	            end if
            end do
        else
            if (mod(decalage*jr+modv,2).EQ.0) then
                i1 = NnTemp + (jr-1)*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i2 = NnTemp + jr*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i3 = NnTemp + jr*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
            else
                i1 = NnTemp + jr*(Ntheta+1) - jr/2
                i2 = NnTemp + jr*(Ntheta+1) - 1 - jr/2
                i3 = NnTemp + (jr+1)*(Ntheta+1) - 1 - jr/2
            end if
            Nfacette=Nfacette+1        
            if (Normale) then
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i3, i2]
            else
                Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i2, i3]
            end if
            do jtheta=1,Ntheta-1
                i1 = NnTemp + jtheta + (jr-1)*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i2 = NnTemp + jtheta + jr*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i3 = NnTemp + jtheta + jr*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i4 = NnTemp + jtheta + (jr-1)*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
	            if (Normale) then
                    Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i4, i2]
		            Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i4, i3, i2]
	            else
	                Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i1, i2, i4]
			        Nfacette=Nfacette+1
                    Mesh%Tfacette(Nfacette)%Tnoeud=[i2, i3, i4]
	            end if
            end do
        end if
    end do
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette= Nfacette
end subroutine mesh_cylhorsym

subroutine mesh_spheresym(Mesh, Rsphere, OrigineSphere, NthS, demi)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    real(rp), intent(in) :: Rsphere
    real(rp), dimension(3) :: OrigineSphere
    integer, intent(in) :: NthS
    logical, optional :: demi ! demi sphere
    
    integer :: jpsi, jphi, Nphi, Npsi, Nfacette, Nnoeud, Nntemp, i1, i2, i3, i4
    real(rp) :: dphi, dPsi, X, Y, Z
    real(rp), allocatable :: phiTab(:), PsiTab(:)
    
    ! This subroutine creates a polar mesh of a sphere with symmetry.
    
    Nphi = NthS
    Npsi = Nphi
    dphi=PI/(Nphi-1)
    dPsi = dphi
    if (present(demi)) then
        if (demi) then
            dPsi = dPsi/2._RP
            !Npsi = Npsi/2._RP
        end if
    end if
    allocate(phiTab(1:Nphi), PsiTab(1:Npsi))
    phiTab(1:Nphi) = (/(jphi*dphi,jphi=0,Nphi-1)/)
    PsiTab(1:Npsi) = (/(Pi/2._RP-jpsi*dPsi,jpsi=0,Npsi-1)/)
    Nnoeud = Mesh%Nnoeud
    Nfacette = Mesh%Nfacette
    Nntemp = Nnoeud
    ! jpsi = 1
    Nnoeud=Nnoeud+1
    Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, RSphere] + OrigineSphere
    do jphi=1,Nphi-1
        Nfacette=Nfacette+1
        Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + [1, jphi+1, jphi+2]
    end do
    do jpsi=2,Npsi-1
        do jphi=1,Nphi
            Nnoeud=Nnoeud+1
            X = Rsphere*cos(phiTab(jphi))*cos(PsiTab(jpsi))
	        Y = Rsphere*sin(phiTab(jphi))*cos(PsiTab(jpsi))
            Z = Rsphere*sin(PsiTab(jpsi))
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + OrigineSphere        
        end do
    end do
    if (demi) then ! jpsi = Npsi
        do jphi=1,Nphi
            Nnoeud=Nnoeud+1
            X = Rsphere*cos(phiTab(jphi))*cos(PsiTab(Npsi))
	        Y = Rsphere*sin(phiTab(jphi))*cos(PsiTab(Npsi))
            Z = Rsphere*sin(PsiTab(Npsi))
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + OrigineSphere        
        end do
    else
        Nnoeud=Nnoeud+1
        Mesh%Tnoeud(Nnoeud)%Pnoeud = [0._RP, 0._RP, -RSphere] + OrigineSphere
        do jphi=1,Nphi-1
            Nfacette = Nfacette + 1
            Mesh%Tfacette(Nfacette)%Tnoeud = Nntemp + 1 + (Npsi-3)*Nphi + [jphi+1, jphi, Nphi+1]
        end do
    end if
    do jpsi=2,Npsi-2
       do jphi=1,Nphi-1  
            i1 = NnTemp + 1 + (jpsi-2)*Nphi + jphi
            i2 = NnTemp + 1 + (jpsi-2)*Nphi + jphi + 1
            i3 = NnTemp + 1 + (jpsi-1)*Nphi + jphi
            i4 = NnTemp + 1 + (jpsi-1)*Nphi + jphi + 1
            if (mod(0*jpsi,2).eq.1) then
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i3, i4]	
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i4, i2]
            else
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i3, i2]	
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i2, i3, i4]
            end if
        end do
    end do
    if (demi) then
        jpsi=Npsi-1
        do jphi=1,Nphi-1  
            i1 = NnTemp + 1 + (jpsi-2)*Nphi + jphi
            i2 = NnTemp + 1 + (jpsi-2)*Nphi + jphi + 1
            i3 = NnTemp + 1 + (jpsi-1)*Nphi + jphi
            i4 = NnTemp + 1 + (jpsi-1)*Nphi + jphi + 1
            if (mod(0*jpsi,2).eq.1) then
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i3, i4]	
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i4, i2]
            else
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i1, i3, i2]	
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = [i2, i3, i4]
            end if
        end do
    end if
    Mesh%Nnoeud = Nnoeud
    Mesh%Nfacette = Nfacette
    
end subroutine mesh_spheresym

subroutine mesh_sphere_unstruct(Mesh, Rsphere, OrigineSphere,InputData,NumBody)
    
    type Tedge
    logical                 :: Done     ! edge which middle point have been determined
    integer, dimension(3)   :: Tnoeud   ! 2 vertex, plus middle index
    integer, dimension(2)   :: Next
    integer                 :: Nfacette
    integer, dimension(2)   :: Tfacette
    real(rp)                :: Length
    real(rp), dimension(3)  :: Normale
    end type Tedge
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    real(rp), intent(in)                :: Rsphere      ! Radius of the sphere
    real(rp), dimension(3)              :: OrigineSphere! Origin of the sphere
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data
    integer,intent(in)                  :: NumBody      ! Body number
    
    logical :: condition
    logical, dimension(PointMax) :: cond
    integer :: j, k, r, s, i1, Nnoeud, Nfacette, Nftemp, Nft, Nedge, Nfe, Np
    integer, dimension(2,3) :: Tabe
    real(rp) :: Dmini, Eps ! minimum size of an edge
    real(rp), dimension(3) :: M, P
    real(rp), allocatable :: theta(:)
    integer, dimension(PointMax,3) :: TabF ! Tableau des arêtes pour chaque facette
    !f2py integer*1, dimension(1000) :: Edge,EdgeTemp
    type(Tedge), dimension(2*PointMax) :: Edge, EdgeTemp
    !f2py integer*1, dimension(1000) :: Facette,FacetteTemp
    type(TFacette), dimension(PointMax) :: Facette, FacetteTemp
    
    ! This subroutine creates a non-polar mesh of a sphere.
    
    Dmini = InputData%Lgeom(2,NumBody) !Mesh%Body(Int_Body)%DimBody(2) !2._RP
    
    Nnoeud = Mesh%Nnoeud
    Nfacette = Mesh%Nfacette
    
    Mesh%Nnoeud = Mesh%Nnoeud + 1
    Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineSphere + [0._RP,RSphere,0._RP]
    Nft = 0
    Nedge = 0
    Np = 4
    allocate(theta(Np))
    theta = (/(2._Rp*Pi/float(Np)*(j-1),j=1,Np)/)
    do j=1,Np
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineSphere + RSphere*[cos(theta(j)),0._RP,-sin(theta(j))]
    end do
    deallocate(theta)
    
    if (Symmetry) then    
        do k=1,Np-1
            Nft = Nft + 1
            Facette(Nft)%Tnoeud = [1,k+1,k+2] + Nnoeud
            Nedge = Nedge + 1
            Edge(Nedge)%Tnoeud(1:2) = [1,k+2] + Nnoeud
            Edge(Nedge)%Nfacette = 2
            Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [Nft, Nft+1]
            Nedge = Nedge + 1
            Edge(Nedge)%Tnoeud(1:2) = [k+1,k+2] + Nnoeud
            Edge(Nedge)%Nfacette = 1
            Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = Nft
        end do
        Nft = Nft + 1
        Facette(Nft)%Tnoeud = [1,2,Np+1] + Nnoeud
        Nedge = Nedge + 1
        Edge(Nedge)%Tnoeud(1:2) = [1,2] + Nnoeud
        Edge(Nedge)%Nfacette = 2
        Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [1,Nft]
        Nedge = Nedge + 1
        Edge(Nedge)%Tnoeud(1:2) = [Np+1,2] + Nnoeud
        Edge(Nedge)%Nfacette = 1
        Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = Nft
    else
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineSphere + [0._RP,-RSphere,0._RP]
        do k=1,Np-1
            Nft = Nft + 1
            Facette(Nft)%Tnoeud = [1,k+1,k+2] + Nnoeud
            Nft = Nft + 1
            Facette(Nft)%Tnoeud = [Np+2,k+1,k+2] + Nnoeud
            Nedge = Nedge + 1
            Edge(Nedge)%Tnoeud(1:2) = [1,k+2] + Nnoeud
            Edge(Nedge)%Nfacette = 2
            Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [Nft-1, Nft+1]
            Nedge = Nedge + 1
            Edge(Nedge)%Tnoeud(1:2) = [Np+2,k+2] + Nnoeud
            Edge(Nedge)%Nfacette = 2
            Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [Nft, Nft+2]
            Nedge = Nedge + 1
            Edge(Nedge)%Tnoeud(1:2) = [k+1,k+2] + Nnoeud
            Edge(Nedge)%Nfacette = 2
            Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [Nft-1, Nft]
        end do
        Nft = Nft + 1
        Facette(Nft)%Tnoeud = [1,2,Np+1] + Nnoeud
        Nft = Nft + 1
        Facette(Nft)%Tnoeud = [Np+2,2,Np+1] + Nnoeud
        Nedge = Nedge + 1
        Edge(Nedge)%Tnoeud(1:2) = [1,2] + Nnoeud
        Edge(Nedge)%Nfacette = 2
        Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [1,Nft-1]
        Nedge = Nedge + 1
        Edge(Nedge)%Tnoeud(1:2) = [Np+2,2] + Nnoeud
        Edge(Nedge)%Nfacette = 2
        Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [2,Nft]
        Nedge = Nedge + 1
        Edge(Nedge)%Tnoeud(1:2) = [Np+1,2] + Nnoeud
        Edge(Nedge)%Nfacette = 2
        Edge(Nedge)%Tfacette(1:Edge(Nedge)%Nfacette) = [Nft-1, Nft]
    end if
    
    call GeomSphere(Mesh, Facette, Edge, TabF, Nft, Nedge)
    Nftemp = Nft
    condition = .true.
    cond = .true.
    do while (condition)
        condition = .false.
        Nft = 0
        Nfe = 0
        
        ! Boucle sur les facettes
        do j=1,Nftemp 
            if (cond(j)) then
                cond(j) = .false.
                ! Boucle sur les arêtes: on cherche les centres des arêtes et on découpe l'arête originale en deux.
                do k=1,3 
                    if (.not.Edge(TabF(j,k))%Done) then
                        Eps = 1._RP
                        P = 0.5_RP*(Mesh%Tnoeud(Edge(TabF(j,k))%Tnoeud(1))%Pnoeud+Mesh%Tnoeud(Edge(TabF(j,k))%Tnoeud(2))%Pnoeud)
                        M = P
                        do while (abs(Eps).gt.Epsilon)
                            M = M + Eps*Edge(TabF(j,k))%Normale
                            Eps = (Rsphere-norm2(M-OrigineSphere))*sign(1._RP,dot_product(M-OrigineSphere,P-OrigineSphere))
                        end do
                        Mesh%Nnoeud = Mesh%Nnoeud + 1                        
                        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = M
                        Edge(TabF(j,k))%Tnoeud(3) = Mesh%Nnoeud
                        Nfe = Nfe + 2
                        EdgeTemp(Nfe-1)%Tnoeud(1:2) = [Edge(TabF(j,k))%Tnoeud(1),Mesh%Nnoeud]
                        EdgeTemp(Nfe)%Tnoeud(1:2) = [Edge(TabF(j,k))%Tnoeud(2),Mesh%Nnoeud]
                        Edge(TabF(j,k))%Next = [Nfe-1, Nfe]
                        Edge(TabF(j,k))%Done = .true.
                    end if
                end do
                ! On relie les centres des arêtes entre eux
                Nfe = Nfe + 3
                EdgeTemp(Nfe-2)%Tnoeud(1:2) = [Edge(TabF(j,1))%Tnoeud(3),Edge(TabF(j,2))%Tnoeud(3)]
                EdgeTemp(Nfe-1)%Tnoeud(1:2) = [Edge(TabF(j,2))%Tnoeud(3),Edge(TabF(j,3))%Tnoeud(3)]
                EdgeTemp(Nfe  )%Tnoeud(1:2) = [Edge(TabF(j,3))%Tnoeud(3),Edge(TabF(j,1))%Tnoeud(3)]
                Nft = Nft+1
                FacetteTemp(Nft)%Tnoeud(1:3) = [Edge(TabF(j,1))%Tnoeud(3),Edge(TabF(j,2))%Tnoeud(3),Edge(TabF(j,3))%Tnoeud(3)]
                EdgeTemp(Nfe-2:Nfe)%Nfacette = 2
                EdgeTemp(Nfe-2:Nfe)%Tfacette(1) = Nft
                
                ! On connecte les 6 arêtes extérieures avec les 3 intérieures pour former 3 nouvelles facettes
                ! boucle sur les nouvelles arêtes centrales
                Tabe = reshape([1,2,2,3,3,1],(/2,3/))
                do k=1,3
                    i1 = 0
                    do r=1,2
                        do s=1,2
                            if (Edge(TabF(j,Tabe(1,k)))%Tnoeud(r).eq.Edge(TabF(j,Tabe(2,k)))%Tnoeud(s)) then
                                i1 = Edge(TabF(j,k))%Tnoeud(r)
                                cycle
                            end if
                        end do
                    end do
                    Nft = Nft+1
                    FacetteTemp(Nft)%Tnoeud(1:3) = [EdgeTemp(Nfe-3+k)%Tnoeud(1:2), i1]
                    EdgeTemp(Nfe-3+k)%Tfacette(2) = Nft
                    do r=1,2
                        if (EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tnoeud(1).eq.i1) then
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette = EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette + 1
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tfacette(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette) = Nft
                        elseif(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tnoeud(1).eq.i1) then
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette = EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette + 1
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tfacette(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette) = Nft
                        else
                            print*, 'erreur connection arete, Next = ', Edge(TabF(j,Tabe(r,k)))%Next, EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tnoeud(1:2)-Nnoeud, EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tnoeud(1:2)-Nnoeud, 'i1 = ', i1-Nnoeud
                            pause
                        end if
                    end do
                end do
            else
                Nft = Nft + 1
                FacetteTemp(Nft) = Facette(j)
            end if
        end do
        Nedge = Nfe
        Edge(1:Nedge) = EdgeTemp(1:Nedge)
        EdgeTemp(:)%Nfacette = 0
        Nftemp = Nft
        Facette(1:Nftemp) = FacetteTemp(1:Nftemp)
        call GeomSphere(Mesh, Facette, Edge, TabF, Nft, Nedge)
        do j=1,Nftemp
            do k=1,3
                if (Edge(TabF(j,k))%Length.gt.Dmini) then 
                    cond(j) = .true.
                    condition = .true.
                    cycle
                end if
            end do
        end do
    end do
    Mesh%Nfacette = Mesh%Nfacette + Nft
    do j=1,Nft
        Mesh%Tfacette(Nfacette+j)%Tnoeud(1:3) = Facette(j)%Tnoeud(1:3)
    end do
contains

    subroutine GeomSphere(Mesh, Facette, Edge, TabF, Nft, Nedge)
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    !f2py integer*1, dimension(1000) :: Edge
    type(Tedge), dimension(2*PointMax) :: Edge
    !f2py integer*1, dimension(1000) :: Facette
    type(TFacette), dimension(PointMax) :: Facette
    integer, dimension(PointMax) :: NabF
    integer, dimension(PointMax,3) :: TabF ! Tableau des arêtes pour chaque facette
    integer, intent(in) :: Nft, Nedge
    
    integer :: j, k, jtemp
    real(rp), dimension(3) :: Nor
    real(rp), dimension(3,3) :: Pg
    
    do j=1,Nft

        Facette(j)%Pnoeud=reshape([Mesh%Tnoeud(Facette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Facette(j)%Tnoeud(2))%Pnoeud,&
                                        & Mesh%Tnoeud(Facette(j)%Tnoeud(3))%Pnoeud],(/3,3/))    
        ! Calcul du centre de gravité de chaque facette, de son aire, de sa normale et de la distance maximale du centre aux sommets    
        Pg=Facette(j)%Pnoeud
        Facette(j)%Gfacette=sum(Pg,2)/3._RP
        call Computation_vect_product(Pg(:,2)-Pg(:,1), Pg(:,3)-Pg(:,1),Nor)
        if (norm2(Nor).lt.Epsilon) then
            print*, 'pb Normale sphere unstructured, Facette : ', j
            print*, Nor
            print*, Pg
            print*, Facette(j)%Tnoeud
            pause
        else
            if (dot_product(Nor,Facette(j)%Gfacette-OrigineSphere).lt.Epsilon) then
                Nor = - Nor
                jtemp = Facette(j)%Tnoeud(2)
                Facette(j)%Tnoeud(2) = Facette(j)%Tnoeud(3)
                Facette(j)%Tnoeud(3) = jtemp
                Facette(j)%Pnoeud=reshape([Mesh%Tnoeud(Facette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Facette(j)%Tnoeud(2))%Pnoeud,&
                                            & Mesh%Tnoeud(Facette(j)%Tnoeud(3))%Pnoeud],(/3,3/))
            end if
            Facette(j)%Normale = Nor/norm2(Nor)
        end if
    end do

    Nabf = 0
    do j=1,Nedge
        Edge(j)%Normale = 0._RP
        do k=1,Edge(j)%Nfacette
            Nabf(Edge(j)%Tfacette(k)) = Nabf(Edge(j)%Tfacette(k)) + 1
            Tabf(Edge(j)%Tfacette(k), Nabf(Edge(j)%Tfacette(k))) = j
            Edge(j)%Normale = Edge(j)%Normale + 0.5_RP*(Facette(Edge(j)%Tfacette(k))%Normale)
        end do
        if (Symmetry .and. Edge(j)%Nfacette.eq.1) Edge(j)%Normale(2) = 0._RP
        Edge(j)%Length = norm2(Mesh%Tnoeud(Edge(j)%Tnoeud(1))%Pnoeud-Mesh%Tnoeud(Edge(j)%Tnoeud(2))%Pnoeud)
        Edge(j)%Done=.false.
    end do
    
end subroutine GeomSphere

end subroutine mesh_sphere_unstruct

subroutine SphereHull(Mesh, Rwh, Twh, profile,InputData,size_Rwh_Twh)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage) :: Mesh
    integer :: size_Rwh_Twh
    real(rp), dimension(size_Rwh_Twh), intent(out)  :: Rwh, Twh
    logical, intent(in) :: profile ! gives only the profile
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                            ! Input data

    logical :: demi ! demi sphere
    integer :: jpsi, jphi, Nphi, Npsi, Nfacette, Nnoeud, Nntemp, i1, i2, i3, i4, NthS
    real(rp) :: dphi, dPsi, X, Y, Z, Rsphere
    real(rp), dimension(3) :: OrigineSphere
    real(rp), allocatable :: phiTab(:), PsiTab(:)

    demi = .true.
    print*,"SphereHull: Lgeom and other physical parameters are from the 1st body."
    Rsphere = InputData%Lgeom(1,1)
    OrigineSphere = InputData%Position(1:3,1,1)
    Nths = InputData%NphiSphere(1)

    Nphi = NthS
    Npsi = Nphi
    dphi=PI/(Nphi-1)
    dPsi = dphi
    !if (present(demi)) then
        if (demi) then
            dPsi = dPsi/2._RP
        end if
    !end if
    allocate(phiTab(1:Nphi), PsiTab(1:Npsi))
    phiTab(1:Nphi) = (/(jphi*dphi,jphi=0,Nphi-1)/)
    PsiTab(1:Npsi) = (/(Pi/2._RP-jpsi*dPsi,jpsi=0,Npsi-1)/)

    if (not(profile)) then
        Mesh%NBody = Mesh%NBody + 1
        Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
        Mesh%Body(Mesh%NBody)%GBody(1:3) = InputData%Position(1:3,1,1)
        Mesh%Body(Mesh%NBody)%MBody = 0._RP
        Mesh%Body(Mesh%NBody)%VBody = 0._RP
        Mesh%Body(Mesh%NBody)%ABody = 0._RP
        Mesh%Body(Mesh%NBody)%DimBody(1:3) = InputData%Lgeom(1:3,1)
        Mesh%Body(Mesh%NBody)%is_tank = .false.
        Nnoeud = Mesh%Nnoeud
        Nfacette = Mesh%Nfacette
        Nntemp = Nnoeud
        ! jpsi = 1
        Nnoeud=Nnoeud+1
        Mesh.Tnoeud(Nnoeud).Pnoeud = [0._RP, 0._RP, -RSphere] + OrigineSphere
        do jphi=1,Nphi-1
            Nfacette=Nfacette+1
            Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + [1, jphi+1, jphi+2]
        end do
        do jpsi=2,Npsi-1
            do jphi=1,Nphi
                Nnoeud=Nnoeud+1
                X = Rsphere*cos(phiTab(jphi))*cos(PsiTab(jpsi))
	            Y = Rsphere*sin(phiTab(jphi))*cos(PsiTab(jpsi))
                Z = -Rsphere*sin(PsiTab(jpsi))
                Mesh.Tnoeud(Nnoeud).Pnoeud = [X, Y, Z] + OrigineSphere        
            end do
        end do
        if (demi) then ! jpsi = Npsi
            do jphi=1,Nphi
                Nnoeud=Nnoeud+1
                X = Rsphere*cos(phiTab(jphi))*cos(PsiTab(Npsi))
	            Y = Rsphere*sin(phiTab(jphi))*cos(PsiTab(Npsi))
                Z = -Rsphere*sin(PsiTab(Npsi))
                Mesh.Tnoeud(Nnoeud).Pnoeud = [X, Y, Z] + OrigineSphere        
            end do
        else
            Nnoeud=Nnoeud+1
            Mesh.Tnoeud(Nnoeud).Pnoeud = [0._RP, 0._RP, RSphere] + OrigineSphere
            do jphi=1,Nphi-1
                Nfacette = Nfacette + 1
                Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + 1 + (Npsi-3)*Nphi + [jphi+1, jphi, Nphi+1]
            end do
        end if
        do jpsi=2,Npsi-2
           do jphi=1,Nphi-1  
                i1 = NnTemp + 1 + (jpsi-2)*Nphi + jphi
                i2 = NnTemp + 1 + (jpsi-2)*Nphi + jphi + 1
                i3 = NnTemp + 1 + (jpsi-1)*Nphi + jphi
                i4 = NnTemp + 1 + (jpsi-1)*Nphi + jphi + 1
                if (mod(0*jpsi,2).eq.1) then
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i3, i4]	
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i4, i2]
                else
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i3, i2]	
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i2, i3, i4]
                end if
            end do
        end do
        if (demi) then
            jpsi=Npsi-1
            do jphi=1,Nphi-1  
                i1 = NnTemp + 1 + (jpsi-2)*Nphi + jphi
                i2 = NnTemp + 1 + (jpsi-2)*Nphi + jphi + 1
                i3 = NnTemp + 1 + (jpsi-1)*Nphi + jphi
                i4 = NnTemp + 1 + (jpsi-1)*Nphi + jphi + 1
                if (mod(0*jpsi,2).eq.1) then
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i3, i4]	
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i4, i2]
                else
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i1, i3, i2]	
                    Nfacette = Nfacette + 1
                    Mesh.Tfacette(Nfacette).Tnoeud = [i2, i3, i4]
                end if
            end do
        end if
        Mesh%Nnoeud = Nnoeud
        Mesh%Nfacette = Nfacette
        Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = 1
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = 1
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 2 !6
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 2 !6
    else
        ! Hull profil output
        Rwh = RSphere
        Twh = phiTab
    end if
    
end subroutine SphereHull

subroutine WigleyHull(Mesh, Rwh, Twh, profile,InputData,size_Rwh_Thw)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage) :: Mesh
    integer,intent(in) :: size_Rwh_Thw
    real(rp), dimension(size_Rwh_Thw), intent(out)  :: Rwh, Twh
    logical, intent(in) :: profile
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData    ! Input data

    integer :: j, k, j1, j2, j3, j4, Nnoeud, NFacette
    integer :: Nx, Nz
    real(RP) :: L, B, D, y
    real(rp), allocatable :: x(:), z(:)

    print*,"WigleyHull: Lgeom and other physical parameters are from the 1st body."
    
    L = InputData%Lgeom(1,1)
    B = InputData%Lgeom(2,1)
    D = InputData%Lgeom(3,1)

    Nx = 0.5*L/Mesh%DimTank(5)
    Nz = 2
    if (Nx*D/L*2 .ge. 2) Nz = Nx*D/L*2
    Nz = Nx/3

    allocate(x(-Nx:Nx),z(0:Nz))
    x = (/(j/float(Nx),j=-Nx,Nx)/)
    z = -(/(j/float(Nz),j=0,Nz)/)

    if (not(profile)) then
        Mesh%NBody = Mesh%NBody + 1
        Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
        Mesh%Body(Mesh%NBody)%GBody(1:3) = InputData%Position(1:3,1,1)
        Mesh%Body(Mesh%NBody)%MBody = 0._RP
        Mesh%Body(Mesh%NBody)%VBody = 0._RP
        Mesh%Body(Mesh%NBody)%ABody = 0._RP
        Mesh%Body(Mesh%NBody)%DimBody(1:3) = InputData%Lgeom(1:3,1)
        Mesh%Body(Mesh%NBody)%is_tank = .false.
        ! Vertex creation
        Nnoeud = Mesh%Nnoeud
        do k=0,Nz
            do j=-Nx,Nx
                Mesh%Nnoeud = Mesh%Nnoeud + 1
                call Hull_function(x(j),z(k), y)
                Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [0.5_RP*L*x(j),0.5_RP*B*y,d*z(k)]
            end do
        end do
        if (not(symmetry)) then
            do k=0,Nz
                do j=-Nx,Nx
                    Mesh%Nnoeud = Mesh%Nnoeud + 1
                    call Hull_function(x(j),z(k), y)
                    Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [0.5_RP*L*x(j),-0.5_RP*B*y,d*z(k)]
                end do
            end do
        end if
        ! Connectivity
        Nfacette = Mesh%Nfacette
        do j=1,Nx
            do k=1,Nz
                j1 = (k-1)*(2*Nx+1) + j
                j2 = (k-1)*(2*Nx+1) + j+1
                j3 =     k*(2*Nx+1) + j
                j4 =     k*(2*Nx+1) + j+1
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + [j1,j2,j3]
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + [j4,j3,j2]
            end do
        end do
        do j=Nx+1,2*Nx
            do k=1,Nz
                j1 = (k-1)*(2*Nx+1) + j
                j2 = (k-1)*(2*Nx+1) + j+1
                j3 =     k*(2*Nx+1) + j
                j4 =     k*(2*Nx+1) + j+1
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + [j1,j2,j4]
                Nfacette = Nfacette + 1
                Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + [j4,j3,j1]
            end do
        end do
        if (not(symmetry)) then
            do j=1,2*Nx
                do k=1,Nz
                    j1 = (k-1)*(2*Nx+1) + j
                    j2 = (k-1)*(2*Nx+1) + j+1
                    j3 =     k*(2*Nx+1) + j
                    j4 =     k*(2*Nx+1) + j+1
                    Nfacette = Nfacette + 1
                    Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + (2*Nx+1)*(Nz+1) + [j1,j2,j3]
                    Nfacette = Nfacette + 1
                    Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + (2*Nx+1)*(Nz+1) + [j4,j3,j2]
                end do
            end do
            do j=1,2*Nx
                do k=1,Nz
                    j1 = (k-1)*(2*Nx+1) + j
                    j2 = (k-1)*(2*Nx+1) + j+1
                    j3 =     k*(2*Nx+1) + j
                    j4 =     k*(2*Nx+1) + j+1
                    Nfacette = Nfacette + 1
                    Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + (2*Nx+1)*(Nz+1) + [j1,j2,j4]
                    Nfacette = Nfacette + 1
                    Mesh%Tfacette(Nfacette)%Tnoeud = Nnoeud + (2*Nx+1)*(Nz+1) + [j4,j3,j1]
                end do
            end do
        end if
        Mesh%Nfacette = Nfacette
        Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = 1
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = 1
        Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 2 !6
        Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 2 !6
    else
        ! Intersection profil output for FS meshing
        k=0
        do j=-Nx,Nx
            k = k + 1
            call Hull_function(x(j),z(0), y)
            Rwh(k) = sqrt((0.5_RP*L*x(j))**2 + (0.5_RP*B*y)**2)
            !Twh(k) = asin(y/Rwh(k))
            Twh(k) = acos(0.5_RP*L*x(j)/Rwh(k))
            !if (y.lt.0) TwH(k) = -Twh(k)
            !if (j .eq. -Nx+1 ) print*, 0.5_RP*L*x(j), 0.5_RP*B*y, Rwh(k), Twh(k)
        end do
        if (not(symmetry)) then
            do j=-Nx+1,Nx-1
                k = k + 1
                call Hull_function(x(j),z(0), y)
                y = -y
                Rwh(k) = sqrt((0.5_RP*L*x(j))**2 + (0.5_RP*B*y)**2)
                !Twh(k) = asin(y/Rwh(k))
                Twh(k) = acos(0.5_RP*L*x(j)/Rwh(k))
                if (y.le.0) then
                    TwH(k) = Twh(k) - PI
                !elseif (y.eq.0 .and. x.gt.0) then
                !    TwH(k) = Twh(k) + PI
                end if
                !if (j .eq. -Nx+1 ) print*, 0.5_RP*L*x(j), 0.5_RP*B*y, Rwh(k), Twh(k)
            end do
        end if
    end if
    deallocate(x, z)
    
end subroutine WigleyHull

subroutine MeshFS_WigleyHull(Mesh,InputData,fgeom_vect,Origine)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage), intent(inout) :: Mesh
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                            ! Input data
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                           ! Geometries
    real(rp), dimension(3), optional, intent(in) :: Origine
    
    logical :: wigley
    integer :: jt, jr, jz, Ntheta, Nr0, Nz, typeN, modu
    real(rp) :: q, dr, Rtemp, Xtemp, Ytemp, dz, Rint
    real(rp), dimension(3) :: Orig, Profil, Pext, vect
    real(rp), allocatable :: Rwh(:), Twh(:), Rtab(:,:), Ttab(:,:), theta(:), Ztab(:), RtabB(:)
    integer :: size_Rwh_Thw,Nx
    real(rp) :: L
    
    wigley = .true.
    Orig=0
    if (present(Origine)) Orig=Origine
    ! Mesh parameters
    Mesh%TypeM = 0                          ! Circular Mesh
    Mesh%DimTank = LDom
    Mesh%Origine = Orig                       ! Origine du maillage

    Mesh%Nnoeud = 0
    Mesh%Nfacette = 0
    if (wigley) then
        ! Wigley Hull mesh
        L = InputData%Lgeom(1,1)
        Nx = 0.5*L/Mesh%DimTank(5)
        if(symmetry)then
            size_Rwh_Thw = 2*Nx+1
            allocate(Rwh(1:2*Nx+1), Twh(1:2*Nx+1))
        else
            size_Rwh_Thw = 4*Nx
            allocate(Rwh(1:4*Nx), Twh(1:4*Nx))
        end if
        call WigleyHull(Mesh, Rwh, Twh, .true.,InputData,size_Rwh_Thw)
        if (symmetry) Twh = PI - Twh
    else
        ! Sphere Mesh
        allocate(Rwh(1:InputData%NphiSphere(1)),Twh(1:InputData%NphiSphere(1)))
        call SphereHull(Mesh, Rwh, Twh, .true.,InputData,InputData%NphiSphere(1))
    end if

    if (allocated(Rwh)) then
        Ntheta = size(Rwh)
    else
        print*, 'MeshFS_WigleyHull : Rwh not allocated, go check it !'
        return
    end if
    !!! Mesh of the Free Surface
    Nr0 = 1.5*Nr
    !Nr = Nr0/2
    dr = Mesh%DimTank(5)
    Rint = Mesh%DimTank(2)-Rwh(1)
    !Nr0 = Rint/dr
    !dr = Rint/float(Nr0-1)
    !if (not(symmetry)) Ntheta = Ntheta-1
    allocate(Rtab(Nr0+Nr,Ntheta), Ttab(Nr0+Nr,Ntheta), theta(Ntheta))
    if (symmetry) then
        theta = (/( PI/float(Ntheta-1)*(jt-1),jt=1,Ntheta )/)
    else
        theta = (/( PI - 2._RP*PI/float(Ntheta)*(jt-1),jt=1,Ntheta )/)
    end if
    q=1._RP+dr
    do jt=1,Ntheta
        Profil = Rwh(jt)*[cos(Twh(jt)), sin(Twh(jt)), 0._RP]
        Pext = Mesh%DimTank(1)*[cos(theta(jt)), sin(theta(jt)), 0._RP]
        ! Détermination du vecteur directeur
        vect(:) = (Pext-Profil) / norm2(Pext-Profil)
        !
        Rtemp = 0._RP
        Xtemp = Profil(1)
        Ytemp = Profil(2)
        Rtab(1,jt) = sqrt(Xtemp**2 + Ytemp**2)
        Ttab(1,jt) = acos(Xtemp/Rtab(1,jt))
        Ttab(1,jt) = sign(Ttab(1,jt),Twh(jt))
        !do jr=2,Nr0
        !    Rtemp = Rtemp + dr
        !    Xtemp = Profil(1) + Rtemp*vect(1)
        !    Ytemp = Profil(2) + Rtemp*vect(2)
        !    Rtab(jr,jt) = sqrt(Xtemp**2 + Ytemp**2)
        !    Ttab(jr,jt) = acos(Xtemp/Rtab(jr,jt))
        !end do
        dr=Mesh%DimTank(5)
        call Serie_Geom(Nr0-1,Rint,dr,q)
        do jr=2,Nr0
            Rtemp = Rtemp + dr*q**(jr-1)
            Xtemp = Profil(1) + Rtemp*vect(1)
            Ytemp = Profil(2) + Rtemp*vect(2)
            Rtab(jr,jt) = sqrt(Xtemp**2 + Ytemp**2)
            Ttab(jr,jt) = acos(Xtemp/Rtab(jr,jt))
            Ttab(jr,jt) = sign(Ttab(jr,jt),Twh(jt))
        end do
        dr = dr*q**(jr-2)
        call Serie_Geom(Nr-1,norm2(Pext-Profil)-Rint,dr,q)
        do jr=Nr0+1,Nr0+Nr-1
            Rtemp = Rtemp + dr*q**(jr-Nr0-1)
            Xtemp = Profil(1) + Rtemp*vect(1)
            Ytemp = Profil(2) + Rtemp*vect(2)
            Rtab(jr,jt) = sqrt(Xtemp**2 + Ytemp**2)
            Ttab(jr,jt) = acos(Xtemp/Rtab(jr,jt))
            Ttab(jr,jt) = sign(Ttab(jr,jt),Twh(jt))
        end do
        Xtemp = Pext(1)
        Ytemp = Pext(2)
        Rtab(Nr0+Nr,jt) = sqrt(Xtemp**2 + Ytemp**2)
        Ttab(Nr0+Nr,jt) = acos(Xtemp/Rtab(Nr0+Nr,jt))
        Ttab(Nr0+Nr,jt) = sign(Ttab(Nr0+Nr,jt),Twh(jt))
    end do
    Mesh%FS%IndFS(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
    if (symmetry) then
        call mesh_cylhorsym_WH(Mesh, Rtab, Ttab, 0._RP, Orig, .true.,Nr0+Nr,Ntheta)
    else
        call mesh_cylhor_WH(Mesh, Rtab, Ttab, 0._RP, Orig, .true.,Nr0+Nr,Ntheta)
    end if
    Mesh%FS%IndFS(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%typeNoeud = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%typeFrontiere = 0
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Npanneau = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%Npanneau = 0

    !!! Mesh of the vertical wall
    Mesh.NBody = Mesh.NBody + 1
    Mesh%Body(Mesh%NBody)%GBody(1:3) = 1000._RP*[0._RP,1._RP,0._RP] !,0._RP,0._RP,0._RP]
    Mesh%Body(Mesh%NBody)%MBody = 0._RP !Mesh%Body(Mesh%NBody)%MBody(1:3)
    Mesh%Body(Mesh%NBody)%VBody = 0._RP
    Mesh%Body(Mesh%NBody)%ABody = 0._RP
    Mesh%Body(Mesh%NBody)%DimBody = LDom(1:3)
    Mesh%Body(Mesh%NBody)%CMD = .false.
    if(DeformMesh) Mesh%Body(Mesh%NBody)%CMD(2) = .true.
    Mesh%Body(Mesh%NBody)%is_tank = .true.
    Mesh%Body(Mesh%NBody)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]

    if (cuve_ferme) then
        dz = Rtab(Nr0+Nr,1)-Rtab(Nr0+Nr-1,1)
        Nz = Mesh%DimTank(3)/dz
        TypeN = 1
    else
        Nz = 1
        TypeN = 10
    end if
    dz = -Mesh%DimTank(3)/float(Nz)
    Nz = Nz+1
    if (Nz.eq.1) then
        Nz = 2
        dz = -Mesh%DimTank(3)
    end if
    allocate(Ztab(Nz))
    Ztab(1:Nz) = (/(dz*(jz-1),jz=1,Nz)/)
    modu =  mod(Nr0+Nr,2)

    if (Symmetry) then
        call mesh_cylvertsym(Mesh, Nz, Ntheta, Ztab, Mesh%DimTank(1), Orig, .true., modu)
    else
        call mesh_cylvert(Mesh, Nz, Ntheta, Ztab, Mesh%DimTank(1), Orig, .true., modu)
    end if
    ! Mesh of the bottom
    if (not(bottom_sym)) then
        if (cuve_ferme) then
            Nr=Nr+1
        else
            Nr=2
            RtabB(1:Nr)=(/( Mesh%DimTank(1)/float(Nr-1)*(jr-1),jr=1,Nr)/)
        end if
        modu = mod(modu+Nr+Nz,2) ! On commence par le centre, d'ou le Nr
        if (Symmetry) then
            call mesh_cylhorsym(Mesh, Nr, Ntheta, RtabB, -Mesh%DimTank(3), Orig, .false., modu)    
        else
            call mesh_cylhor(Mesh, Nr, Ntheta, RtabB, -Mesh%DimTank(3), Orig, .false., modu)
        end if
    end if
    Mesh%Body(Mesh%NBody)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%typeNoeud = TypeN
    Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%typeFrontiere = TypeN
    Mesh%Tnoeud(Mesh%Body(Mesh%NBody)%IndBody(1):Mesh%Body(Mesh%NBody)%IndBody(3))%Npanneau = 1
    Mesh%Tfacette(Mesh%Body(Mesh%NBody)%IndBody(2):Mesh%Body(Mesh%NBody)%IndBody(4))%Npanneau = 1

    ! Deallocation
    deallocate(Rwh,Twh,Rtab,Ttab,theta, Ztab)

    if (wigley) then
        ! Wigley Hull mesh
        L = InputData%Lgeom(1,1)
        Nx = 0.5*L/Mesh%DimTank(5)
        if(symmetry)then
            size_Rwh_Thw = 2*Nx+1
            allocate(Rwh(1:2*Nx+1), Twh(1:2*Nx+1))
        else
            size_Rwh_Thw = 4*Nx
            allocate(Rwh(1:4*Nx), Twh(1:4*Nx))
        end if
        call WigleyHull(Mesh, Rwh, Twh, .false.,InputData,size_Rwh_Thw)
    else
        ! Sphere Mesh
        allocate(Rwh(1:InputData%NphiSphere(1)),Twh(1:InputData%NphiSphere(1)))
        call SphereHull(Mesh, Rwh, Twh, .false.,InputData,InputData%NphiSphere(1))
    end if
    !
    if (cuve_ferme) then
        Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3)
        Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4)
    else
        if (is_body.and.Int_Body.ne.0) then
            Mesh%Nsys = Mesh%Body(Int_Body)%IndBody(3)
            Mesh%Nfsys = Mesh%Body(Int_Body)%IndBody(4)
        else
            Mesh%Nsys = Mesh%FS%IndFS(3)
            Mesh%Nfsys = Mesh%FS%IndFS(4)
        end if
    end if
    call GeomInit(Mesh, fgeom_vect, 0._RP,InputData, .false.)
    
end subroutine MeshFS_WigleyHull


subroutine mesh_cylhor_WH(Mesh, Rtab, Ttab, Z, Orig, Normale,size_X,size_Y)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage) :: Mesh
    integer :: size_X,size_Y
    real(rp), dimension(size_X,size_Y), intent(in) :: Rtab, Ttab
    real(rp), intent(in) :: Z, Orig(3)
    logical, intent(in) :: Normale ! définit le sens de la normale
    
    integer :: Nr, Ntheta
    integer :: jr, jtheta, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu, modv
    real(rp) :: X, Y, dtheta
    
    Nr = size(Rtab,1)
    Ntheta = size(Rtab,2)
    if (Nr.ne.size(Ttab,1) .or. Ntheta.ne.size(Ttab,2)) print*, 'mesh_cylhor_WH : Rtab et Ttab de taille différente --> go check it!' 

    Nnoeud=Mesh.Nnoeud
    Nntemp=Nnoeud
    NFacette=Mesh.NFacette
    modu = 1
    modv = modu
    !if (RTab(1).lt.Epsilon2) then
    !    Nnoeud=Nnoeud+1
    !    Mesh.Tnoeud(Nnoeud).Pnoeud = [0._RP, 0._RP, Z]+Orig
    !    do jtheta=1,Ntheta-1
    !        Nfacette=Nfacette+1
    !        if (Normale) then
    !            Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + [1, jtheta+2, jtheta+1]
    !        else
    !            Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + [1, jtheta+1, jtheta+2]
    !        end if
    !    end do
    !    Nfacette=Nfacette+1
    !    if (Normale) then
    !        Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + [1, 2, Ntheta+1]
    !    else
    !        Mesh.Tfacette(Nfacette).Tnoeud = Nntemp + [1, Ntheta+1, 2]
    !    end if
    !end if
    NnTemp = Nnoeud
    do jr=1,Nr
        !if (RTab(jr).gt.Epsilon2) then
            do jtheta=1,Ntheta
	            Nnoeud=Nnoeud+1
                if (jtheta.eq.Ntheta) then
                    dtheta = (Ttab(jr,jtheta)-Ttab(jr,jtheta-1))*0.5_RP*MOD(jr+modu,2)
                else
                    dtheta = (Ttab(jr,jtheta+1)-Ttab(jr,jtheta))*0.5_RP*MOD(jr+modu,2)
                end if
	            X = RTab(jr,jtheta)*COS(Ttab(jr,jtheta)+0.5_rp*dtheta*MOD(decalage*jr+modu,2))
		        Y = RTab(jr,jtheta)*SIN(Ttab(jr,jtheta)+0.5_rp*dtheta*MOD(decalage*jr+modu,2))
	            Mesh.Tnoeud(Nnoeud).Pnoeud = [X, Y, Z] + Orig
            end do
        !end if
    end do
    !if (RTab(1).lt.Epsilon2) then
    !    Nr=Nr-1
    !    modv = modu+1
    !end if
    do jr=1,Nr-1
        do jtheta=1,Ntheta-1
            i1 = NnTemp + jtheta+(jr-1)*Ntheta
            if (Normale) then
	            i2 = NnTemp + jtheta+jr*Ntheta
	            i3 = NnTemp + jtheta+1+(jr-1)*Ntheta
	        else
	            i2 = NnTemp + jtheta+1+(jr-1)*Ntheta
	            i3 = NnTemp + jtheta+jr*Ntheta
	        end if
	        i4 = NnTemp + jtheta+1+jr*Ntheta
	        if (MOD(decalage*jr+modv,2).EQ.1) THEN
                Nfacette=Nfacette+1
                Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i4]
		        Nfacette=Nfacette+1
                Mesh.Tfacette(Nfacette).Tnoeud=[i1, i4, i2]
	        ELSE
	            Nfacette=Nfacette+1
                Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i2]
			    Nfacette=Nfacette+1
                Mesh.Tfacette(Nfacette).Tnoeud=[i2, i3, i4]
	        end if
        end do
        i1 = NnTemp + jr*Ntheta
        if (Normale) then
            i2 = NnTemp + (jr+1)*Ntheta
            i3 = NnTemp + (jr-1)*Ntheta+1
        else
            i2 = NnTemp + (jr-1)*Ntheta+1
            i3 = NnTemp + (jr+1)*Ntheta
        end if
        i4 = NnTemp + jr*Ntheta+1
        if (MOD(decalage*jr+modv,2).EQ.1) THEN
            Nfacette=Nfacette+1
            Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i4]	
	        Nfacette=Nfacette+1
            Mesh.Tfacette(Nfacette).Tnoeud=[i1, i4, i2]
        ELSE
	        Nfacette=Nfacette+1
            Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i2]
	        Nfacette=Nfacette+1
            Mesh.Tfacette(Nfacette).Tnoeud=[i2, i3, i4]
        end if
    end do
    Mesh.Nnoeud = Nnoeud
    Mesh.Nfacette= Nfacette
    
end subroutine mesh_cylhor_WH

subroutine mesh_cylhorsym_WH(Mesh, Rtab, Ttab, Z, Orig, Normale,size_X,size_Y)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage) :: Mesh
    integer :: size_X,size_Y
    real(rp), dimension(size_X,size_Y), intent(in) :: Rtab, Ttab
    real(rp), intent(in) :: Z
    real(rp), dimension(3), intent(in) :: Orig
    logical, intent(in) :: Normale ! définit le sens de la normale
    
    integer :: jr, jt, Nnoeud, Nfacette, Nntemp, i1, i2, i3, i4, modu, modv, Nr, Ntheta
    real(rp) :: X, Y, dtheta
    
    Nr = size(Rtab,1)
    Ntheta = size(Rtab,2)
    if (Nr.ne.size(Ttab,1) .or. Ntheta.ne.size(Ttab,2)) print*, 'mesh_cylhorsym_WH : Rtab et Ttab de taille différente --> go check it!' 

    Nnoeud=Mesh%Nnoeud
    Nntemp=Nnoeud
    NFacette=Mesh%NFacette
    modu = 1
    modv = modu
    ! Vertex definition
    do jr=1,Nr
        if (decalage.eq.0 .or. MOD(jr+modu,2).eq.0) then
            do jt=1,Ntheta
	            Nnoeud=Nnoeud+1
	            X = RTab(jr,jt)*COS(Ttab(jr,jt))
		        Y = RTab(jr,jt)*SIN(Ttab(jr,jt))
	            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
                !if (jt.eq.1 .and. jr.eq.2) print*, X, Y, RTab(jr,jt), TTab(jr,jt)
            end do
        else
            Nnoeud=Nnoeud+1
            Mesh.Tnoeud(Nnoeud).Pnoeud = [RTab(jr,1), 0._RP, Z] + Orig
            do jt=1,Ntheta-1
	            Nnoeud=Nnoeud+1
                dtheta = (Ttab(jr,jt+1)-Ttab(jr,jt))*0.5_RP*MOD(jr+modu,2)
	            X = RTab(jr,jt)*COS(Ttab(jr,jt)+dtheta)
		        Y = RTab(jr,jt)*SIN(Ttab(jr,jt)+dtheta)
	            Mesh%Tnoeud(Nnoeud)%Pnoeud = [X, Y, Z] + Orig
            end do
            Nnoeud=Nnoeud+1
            Mesh%Tnoeud(Nnoeud)%Pnoeud = [-RTab(jr,Ntheta), 0._RP, Z] + Orig
        end if
    end do
    ! Facettes
    do jr=1,Nr-1
        if (decalage.eq.0) then
           ! do jt=1,Ntheta-1! Ntheta/2
           !     i1 = NnTemp + jt+(jr-1)*Ntheta
           !     if (Normale) then
	          !      i2 = NnTemp + jt+jr*Ntheta
	          !      i3 = NnTemp + jt+1+(jr-1)*Ntheta
	          !  else
	          !      i2 = NnTemp + jt+1+(jr-1)*Ntheta
	          !      i3 = NnTemp + jt+jr*Ntheta
	          !  end if
	          !  i4 = NnTemp + jt+1+jr*Ntheta
	          !  if (mod(0*jr+modv,2).EQ.1) then
           !         Nfacette=Nfacette+1
           !         Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i4]
		         !   Nfacette=Nfacette+1
           !         Mesh.Tfacette(Nfacette).Tnoeud=[i1, i4, i2]
	          !  else
	          !      Nfacette=Nfacette+1
           !         Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i2]
			        !Nfacette=Nfacette+1
           !         Mesh.Tfacette(Nfacette).Tnoeud=[i2, i3, i4]
	          !  end if
           ! end do
            do jt=1,Ntheta-1!Ntheta/2+1,Ntheta-1
                i1 = NnTemp + jt+(jr-1)*Ntheta
                if (Normale) then
	                i2 = NnTemp + jt+jr*Ntheta
	                i3 = NnTemp + jt+1+(jr-1)*Ntheta
	            else
	                i2 = NnTemp + jt+1+(jr-1)*Ntheta
	                i3 = NnTemp + jt+jr*Ntheta
	            end if
	            i4 = NnTemp + jt+1+jr*Ntheta
	            if (mod(0*jr+modv,2).EQ.1) then
	                Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i2]
			        Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i2, i3, i4]
	            else
                    Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i4]
		            Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i1, i4, i2]
	            end if
            end do
        else
            if (mod(decalage*jr+modv,2).EQ.0) then
                i1 = NnTemp + (jr-1)*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i2 = NnTemp + jr*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i3 = NnTemp + jr*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
            else
                i1 = NnTemp + jr*(Ntheta+1) - jr/2
                i2 = NnTemp + jr*(Ntheta+1) - 1 - jr/2
                i3 = NnTemp + (jr+1)*(Ntheta+1) - 1 - jr/2
            end if
            Nfacette=Nfacette+1        
            if (Normale) then
                Mesh.Tfacette(Nfacette).Tnoeud=[i1, i3, i2]
            else
                Mesh.Tfacette(Nfacette).Tnoeud=[i1, i2, i3]
            end if
            do jt=1,Ntheta-1
                i1 = NnTemp + jt + (jr-1)*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i2 = NnTemp + jt + jr*(Ntheta+1) - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i3 = NnTemp + jt + jr*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
                i4 = NnTemp + jt + (jr-1)*(Ntheta+1) + 1 - jr/2 + MOD(modv+1,2)*mod(decalage*(jr+1)+modv,2)
	            if (Normale) then
                    Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i1, i4, i2]
		            Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i4, i3, i2]
	            else
	                Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i1, i2, i4]
			        Nfacette=Nfacette+1
                    Mesh.Tfacette(Nfacette).Tnoeud=[i2, i3, i4]
	            end if
            end do
        end if
    end do
    Mesh.Nnoeud = Nnoeud
    Mesh.Nfacette= Nfacette
    
end subroutine mesh_cylhorsym_WH

subroutine mesh_FS_unstruc(Mesh, RFS, OrigineFS)
    
    type Tedge
    logical                 :: Done     ! edge which middle point have been determined
    integer, dimension(3)   :: Tnoeud   ! 2 vertex, plus middle index
    integer, dimension(2)   :: Next
    integer                 :: Nfacette
    integer, dimension(2)   :: Tfacette
    real(rp)                :: Length
    real(rp), dimension(3)  :: Normale
    end type Tedge
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    real(rp), intent(in) :: RFS
    real(rp), dimension(3) :: OrigineFS
    
    logical :: condition
    logical, dimension(PointMax) :: cond
    integer :: j, k, r, s, i1, Nnoeud, Nfacette, Nftemp, Nft, Nedge, Nfe
    integer, dimension(2,3) :: Tabe
    real(rp) :: Dmini, Eps ! minimum size of an edge
    real(rp), dimension(3) :: M, P
    integer, dimension(PointMax,3) :: TabF ! Tableau des arêtes pour chaque facette
    !f2py integer*1, dimension(1000) :: Edge,EdgeTemp
    type(Tedge), dimension(2*PointMax) :: Edge, EdgeTemp
    !f2py integer*1, dimension(1000) :: Facette,FacetteTemp
    type(TFacette), dimension(PointMax) :: Facette, FacetteTemp
    Dmini = Mesh%FS%DimFS(5)
    
    ! This subroutine creates an unstructured mesh of the upper disc of a cylindrical tank.
    
    Nnoeud = Mesh%Nnoeud
    Nfacette = Mesh%Nfacette
    
    ! Initialisation du triangle
    Nft = 0
    Nedge = 0
    Mesh%Nnoeud = Mesh%Nnoeud + 1
    Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineFS + [RFS,0._RP,0._RP]
    if (symmetry) then
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineFS + [-RFS,0._RP,0._RP]
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineFS + [0._RP,RFS,0._RP]
    else
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineFS + 0.5_RP*sqrt(2._RP)*RFS*[-1._RP,1._RP,0._RP]
        Mesh%Nnoeud = Mesh%Nnoeud + 1
        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = OrigineFS + 0.5_RP*sqrt(2._RP)*RFS*[-1._RP,-1._RP,0._RP]
    end if
    Nft = Nft + 1
    Facette(Nft)%Tnoeud = [1,2,3] + Nnoeud
    Nedge = Nedge + 1
    Edge(Nedge)%Tnoeud(1:2) = [1,2] + Nnoeud
    Edge(Nedge)%Nfacette = 1
    Edge(Nedge)%Tfacette(1) = 1
    Nedge = Nedge + 1
    Edge(Nedge)%Tnoeud(1:2) = [1,3] + Nnoeud
    Edge(Nedge)%Nfacette = 1
    Edge(Nedge)%Tfacette(1) = 1
    Nedge = Nedge + 1
    Edge(Nedge)%Tnoeud(1:2) = [3,2] + Nnoeud
    Edge(Nedge)%Nfacette = 1
    Edge(Nedge)%Tfacette(1) = 1
    call GeomFS(Mesh, Facette, Edge, TabF, Nft, Nedge)
    
    ! Initialisation de la boucle
    Nftemp = Nft
    condition = .true.
    cond = .true.
    do while (condition)
        condition = .false.
        Nft = 0
        Nfe = 0
        do j=1,Nftemp ! Boucle sur les facettes
            if (cond(j)) then
                cond(j) = .false.
                do k=1,3 ! Boucle sur les arêtes: on cherche les centres des arêtes et on découpe l'arête originale en deux.
                    if (.not.Edge(TabF(j,k))%Done) then
                        Eps = 1._RP
                        P = 0.5_RP*(Mesh%Tnoeud(Edge(TabF(j,k))%Tnoeud(1))%Pnoeud+Mesh%Tnoeud(Edge(TabF(j,k))%Tnoeud(2))%Pnoeud)
                        M = P
                        Mesh%Nnoeud = Mesh%Nnoeud + 1
                        Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = M
                        Edge(TabF(j,k))%Tnoeud(3) = Mesh%Nnoeud
                        Nfe = Nfe + 2
                        EdgeTemp(Nfe-1)%Tnoeud(1:2) = [Edge(TabF(j,k))%Tnoeud(1),Mesh%Nnoeud]
                        EdgeTemp(Nfe)%Tnoeud(1:2) = [Edge(TabF(j,k))%Tnoeud(2),Mesh%Nnoeud]
                        Edge(TabF(j,k))%Next = [Nfe-1, Nfe]
                        Edge(TabF(j,k))%Done = .true.
                    end if
                end do
                ! On relie les centres des arêtes entre eux
                Nfe = Nfe + 3
                EdgeTemp(Nfe-2)%Tnoeud(1:2) = [Edge(TabF(j,1))%Tnoeud(3),Edge(TabF(j,2))%Tnoeud(3)]
                EdgeTemp(Nfe-1)%Tnoeud(1:2) = [Edge(TabF(j,2))%Tnoeud(3),Edge(TabF(j,3))%Tnoeud(3)]
                EdgeTemp(Nfe  )%Tnoeud(1:2) = [Edge(TabF(j,3))%Tnoeud(3),Edge(TabF(j,1))%Tnoeud(3)]
                Nft = Nft+1
                FacetteTemp(Nft)%Tnoeud(1:3) = [Edge(TabF(j,1))%Tnoeud(3),Edge(TabF(j,2))%Tnoeud(3),Edge(TabF(j,3))%Tnoeud(3)]
                EdgeTemp(Nfe-2:Nfe)%Nfacette = 2
                EdgeTemp(Nfe-2:Nfe)%Tfacette(1) = Nft
                ! On connecte les 6 arêtes extérieures avec les 3 intérieures pour former 3 nouvelles facettes
                ! boucle sur les nouvelles arêtes centrales
                Tabe = reshape([1,2,2,3,3,1],(/2,3/))
                do k=1,3
                    i1 = 0
                    do r=1,2
                        do s=1,2
                            if (Edge(TabF(j,Tabe(1,k)))%Tnoeud(r).eq.Edge(TabF(j,Tabe(2,k)))%Tnoeud(s)) then
                                i1 = Edge(TabF(j,k))%Tnoeud(r)
                                cycle
                            end if
                        end do
                    end do
                    Nft = Nft+1
                    FacetteTemp(Nft)%Tnoeud(1:3) = [EdgeTemp(Nfe-3+k)%Tnoeud(1:2), i1]
                    EdgeTemp(Nfe-3+k)%Tfacette(2) = Nft
                    do r=1,2
                        if (EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tnoeud(1).eq.i1) then
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette = EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette + 1
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tfacette(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Nfacette) = Nft
                        elseif(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tnoeud(1).eq.i1) then
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette = EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette + 1
                            EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tfacette(EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Nfacette) = Nft
                        else
                            print*, 'erreur connection arete, Next = ', Edge(TabF(j,Tabe(r,k)))%Next, EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(1))%Tnoeud(1:2)-Nnoeud, EdgeTemp(Edge(TabF(j,Tabe(r,k)))%Next(2))%Tnoeud(1:2)-Nnoeud, 'i1 = ', i1-Nnoeud
                            pause
                        end if
                    end do
                end do
            else
                Nft = Nft + 1
                FacetteTemp(Nft) = Facette(j)
            end if
        end do
        Nedge = Nfe
        Edge(1:Nedge) = EdgeTemp(1:Nedge)
        EdgeTemp(:)%Nfacette = 0
        Nftemp = Nft
        Facette(1:Nftemp) = FacetteTemp(1:Nftemp)
        call GeomFS(Mesh, Facette, Edge, TabF, Nft, Nedge)
        do j=1,Nftemp
            do k=1,3
                if (Edge(TabF(j,k))%Length.gt.Dmini) then 
                    cond(j) = .true.
                    condition = .true.
                    cycle
                end if
            end do
        end do
    end do
    Mesh%Nfacette = Mesh%Nfacette + Nft
    do j=1,Nft
        Mesh%Tfacette(Nfacette+j)%Tnoeud(1:3) = Facette(j)%Tnoeud(1:3)
    end do

    contains
    subroutine GeomFS(Mesh, Facette, Edge, TabF, Nft, Nedge)
    ! Parameters
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    !f2py integer*1, dimension(1000) :: Edge
    type(Tedge), dimension(2*PointMax) :: Edge
    !f2py integer*1, dimension(1000) :: Facette
    type(TFacette), dimension(PointMax) :: Facette
    integer, dimension(PointMax) :: NabF
    integer, dimension(PointMax,3) :: TabF ! Tableau des arêtes pour chaque facette
    integer, intent(in) :: Nft, Nedge
    ! Variables
    integer :: j, k, jtemp
    real(rp), dimension(3) :: Nor
    real(rp), dimension(3,3) :: Pg
    ! Begin
    do j=1,Nft
        Facette(j)%Pnoeud=reshape([Mesh%Tnoeud(Facette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Facette(j)%Tnoeud(2))%Pnoeud,&
                                        & Mesh%Tnoeud(Facette(j)%Tnoeud(3))%Pnoeud],(/3,3/))    
        ! Calcul du centre de gravité de chaque facette, de son aire, de sa normale et de la distance maximale du centre aux sommets    
        Pg=Facette(j)%Pnoeud
        Facette(j)%Gfacette=sum(Pg,2)/3._RP
        call Computation_vect_product(Pg(:,2)-Pg(:,1), Pg(:,3)-Pg(:,1),Nor)
        if (norm2(Nor).lt.Epsilon) then
            print*, 'pb Normale sphere unstructured, Facette : ', j
            print*, Nor
            print*, Pg
            print*, Facette(j)%Tnoeud
            pause
        else
            if (Nor(3).gt.0) then
                Nor = - Nor
                jtemp = Facette(j)%Tnoeud(2)
                Facette(j)%Tnoeud(2) = Facette(j)%Tnoeud(3)
                Facette(j)%Tnoeud(3) = jtemp
                Facette(j)%Pnoeud=reshape([Mesh%Tnoeud(Facette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Facette(j)%Tnoeud(2))%Pnoeud,&
                                            & Mesh%Tnoeud(Facette(j)%Tnoeud(3))%Pnoeud],(/3,3/))
            end if
            Facette(j)%Normale = Nor/norm2(Nor)
        end if
    end do

    Nabf = 0
    do j=1,Nedge
        Edge(j)%Normale = 0._RP
        do k=1,Edge(j)%Nfacette
            Nabf(Edge(j)%Tfacette(k)) = Nabf(Edge(j)%Tfacette(k)) + 1
            Tabf(Edge(j)%Tfacette(k), Nabf(Edge(j)%Tfacette(k))) = j
            Edge(j)%Normale = Edge(j)%Normale + 0.5_RP*(Facette(Edge(j)%Tfacette(k))%Normale)
        end do
        if (Symmetry .and. Edge(j)%Nfacette.eq.1) Edge(j)%Normale(2) = 0._RP
        Edge(j)%Length = norm2(Mesh%Tnoeud(Edge(j)%Tnoeud(1))%Pnoeud-Mesh%Tnoeud(Edge(j)%Tnoeud(2))%Pnoeud)
        Edge(j)%Done=.false.
    end do
    
    end subroutine geomFS

end subroutine mesh_FS_unstruc

subroutine Mesh_cart_Wigley(Mesh, InputData,fgeom_vect,Origine)
    
    !f2py integer*1, dimension(1000):: Mesh
    type(TMaillage), intent(inout) :: Mesh
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                            ! Input data
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                           ! Geometries
    real(rp), dimension(3), optional, intent(in) :: Origine
    
    integer :: j, k, j1, j2, j3, j4
    integer :: Nx, Nxe, Nxa, Nxtotal, Ny, Nz, Nnoeud
    real(rp) :: dx, dxa, dy, dz, q, maxValue, X, Y, Lxe, Crampe1
    real(rp), dimension(:), allocatable :: Xtab, Ytab, Ztab, Xtemp, Rwh, Twh
    integer :: size_Rwh_Thw,Nx2
    real(rp) :: L

    print*,"Mesh_cart_Wigley: Lgeom and other physical parameters are from the 1st body."
    
    ! Free-surface mesh
    Mesh%DimTank = LDom
    ! Ldom(1) : length of the domain after the hull with no mesh elongation
    ! Ldom(2) : width of the domain
    ! Ldom(3) : depth of the domain
    ! Ldom(4) : length of the domain with mesh elongation, upstream and downstream (damping beach)
    ! Ldom(5) : mesh discretization on the FS and body
    Mesh%Origine = Origine

    dx = Mesh%DimTank(5)
    ! Discretization for the mesh elongation
    Lxe = Mesh%DimTank(4) ! Lxe = 0.25_RP*Lgeom(1)
    Nxe = 0.75_RP*Lxe/dx ! Condition to vary to get mesh elongation (0.5< <=1, 1 --> no elongation)
    allocate(Xtemp(Nxe+1))
    q = 1.5_RP
    call Serie_Geom(Nxe-1,Lxe,dx,q)
    Xtemp(1) = dx
    do j=2,Nxe
        Xtemp(j) = Xtemp(j-1) + dx*q**(j-1)
    end do
    maxValue = Xtemp(Nxe) - Xtemp(Nxe-1)

    ! Downstream discretization (no mesh elongation)
    dxa = maxValue
    Nxa = Mesh%DimTank(1)/dxa
    dxa = Mesh%DimTank(1)/Nxa
    ! Hull discretization
    Nx = InputData%Lgeom(1,1)/dx
    if (mod(Nx,2).eq.1) Nx = Nx - 1
    dx = InputData%Lgeom(1,1)/Nx
    allocate(Xtab(Nxa+2*Nxe+Nx+1))
    Xtab(1:Nxa) = (/(dxa*(j-1),j=1,Nxa)/) - 0.5*InputData%Lgeom(1,1) - Lxe - Mesh%DimTank(1)
    Xtab(Nxa+1:Nxa+Nxe) = -Xtemp(Nxe:1:-1) - 0.5*InputData%Lgeom(1,1)
    Xtab(Nxa+Nxe+1:Nxa+Nxe+Nx) = (/(dx*(j-1),j=1,Nx)/) - 0.5*InputData%Lgeom(1,1)
    Xtab(Nxa+Nxe+Nx+1) = 0.5*InputData%Lgeom(1,1)
    ! Upstream discretization
    Xtab(Nxa+Nxe+Nx+2:Nxa+Nxe+Nx+Nxe+1) = Xtemp(1:Nxe) + 0.5*InputData%Lgeom(1,1)
    Nxtotal = Nxa+2*Nxe+Nx+1
    Mesh%FS%DimFS(1) = Xtab(Nxtotal) - Xtab(1)

    ! Transversal discretization
    dy = Mesh%DimTank(5)
    Ny = 0.75_RP*Mesh%DimTank(2)/dy ! Condition to vary to get mesh elongation (0.5< <=1, 1 --> no elongation)
    allocate(Ytab(Ny+1))
    call Serie_Geom(Ny-1,Mesh%DimTank(2),dy,q)
    Ytab(1) = 0
    Ytab(2) = dy
    do j=2,Ny
        Ytab(j+1) = Ytab(j) + dy*q**(j-1)
    end do
    Ny = Ny + 1
    Mesh%FS%DimFS(2) = Mesh%DimTank(2)
    maxValue = max(maxValue,Ytab(Ny) - Ytab(Ny-1))

    ! Creation of the points
    Mesh%Nnoeud = 0
    Mesh%Nfacette = 0
    do j=1,Ny
        do k=1,Nxtotal
            Mesh%Nnoeud = Mesh%Nnoeud + 1
            X = Xtab(k) + Origine(1)
            call Hull_function(-1._RP,0._RP, Y)
            Y = 0.5_RP*InputData%Lgeom(2,1)*Y + Origine(2)
            if (X.lt.-0.5_RP*InputData%Lgeom(1,1) .or. X .gt. 0.5_RP*InputData%Lgeom(1,1)) then
                call Computation_Crampe(Ytab(j),Ytab(1),Ytab(Ny),Crampe1)
                Y = Ytab(j) + (1._RP-Crampe1)*Y
            else
                call Hull_function(X/(0.5_RP*InputData%Lgeom(1,1)),0._RP, Y)
                Y = 0.5_RP*InputData%Lgeom(2,1)*Y + Origine(2)
                call Computation_Crampe(Ytab(j),Ytab(1),Ytab(Ny),Crampe1)
                Y = Ytab(j) + (1._RP-Crampe1)*Y
            end if
            Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [X,Y,0._RP]
        end do
    end do
    ! Connectivity
    do k=1,Ny-1
        do j=1,Nxtotal-1
            j1 = (k-1)*Nxtotal + j
            j2 = (k-1)*Nxtotal + j+1
            j3 =     k*Nxtotal + j
            j4 =     k*Nxtotal + j+1
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = [j1,j3,j2]
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = [j4,j2,j3]
        end do
    end do
    Mesh%FS%IndFS(1:4) = [1, 1, Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%typeNoeud = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%typeFrontiere = 0
    Mesh%Tnoeud(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Npanneau = 0
    Mesh%Tfacette(Mesh%FS%IndFS(2):Mesh%FS%IndFS(4))%Npanneau = 0

    ! Walls Mesh
    Mesh%NBody = 1
    Mesh%Body(1)%IndBody(1:2) = [Mesh%Nnoeud+1, Mesh%Nfacette+1]
    Mesh%Body(Mesh%NBody)%GBody(1:3) = 1000._RP*[0._RP,1._RP,0._RP]
    Mesh%Body(Mesh%NBody)%MBody = 0._RP
    Mesh%Body(Mesh%NBody)%VBody = 0._RP
    Mesh%Body(Mesh%NBody)%ABody = 0._RP
    Mesh%Body(Mesh%NBody)%DimBody = LDom(1:3)
    Mesh%Body(Mesh%NBody)%CMD = .false.
    if(DeformMesh) Mesh%Body(Mesh%NBody)%CMD(2) = .true.
    Mesh%Body(Mesh%NBody)%is_tank = .true.
    dz = maxValue
    dz = dy ! Comment if you want larger elements on the wall
    if (0.5*Mesh%DimTank(3)/dz.lt.2) then
        Nz = 2
        allocate(Ztab(2))
        Ztab = [0._RP,-Mesh%DimTank(3)]
    else
        Nz = 0.5*Mesh%DimTank(3)/dz
        allocate(Ztab(Nz+1))
        call Serie_Geom(Nz-1,Mesh%DimTank(3),dz,q)
        Ztab(1) = 0._RP
        Ztab(2) = -dz
        do j=2,Nz
            Ztab(j+1) = Ztab(j) - dz*q**(j-1)
        end do
        Nz = Nz + 1
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Upstream Wall (x<0)
    Nnoeud = Mesh%Nnoeud
    do j=1,Ny
        do k=1,Nz
            Mesh%Nnoeud = Mesh%Nnoeud + 1
            Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [- 0.5*InputData%Lgeom(1,1) - Lxe - Ldom(1), Ytab(j), Ztab(k)]
        end do
    end do
    ! Connectivity
    do k=1,Ny-1
        do j=1,Nz-1
            j1 = (k-1)*Nz + j
            j2 = (k-1)*Nz + j+1
            j3 =     k*Nz + j
            j4 =     k*Nz + j+1
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j1,j2,j3]
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j4,j3,j2]
        end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lateral Wall (y>0)
    Nnoeud = Mesh%Nnoeud
    do j=1,Nxtotal
        do k=1,Nz
            Mesh%Nnoeud = Mesh%Nnoeud + 1
            Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [Xtab(j), Ldom(2), Ztab(k)]
        end do
    end do
    ! Connectivity
    do k=1,Nxtotal-1
        do j=1,Nz-1
            j1 = (k-1)*Nz + j
            j2 = (k-1)*Nz + j+1
            j3 =     k*Nz + j
            j4 =     k*Nz + j+1
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j1,j2,j3]
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j4,j3,j2]
        end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Downstream Wall (x>0)
    Nnoeud = Mesh%Nnoeud
    do j=1,Ny
        do k=1,Nz
            Mesh%Nnoeud = Mesh%Nnoeud + 1
            Mesh%Tnoeud(Mesh%Nnoeud)%Pnoeud = [0.5_RP*InputData%Lgeom(1,1) + Lxe, Ytab(j), Ztab(k)]
        end do
    end do
    ! Connectivity
    do k=1,Ny-1
        do j=1,Nz-1
            j1 = (k-1)*Nz + j
            j2 = (k-1)*Nz + j+1
            j3 =     k*Nz + j
            j4 =     k*Nz + j+1
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j1,j3,j2]
            Mesh%Nfacette = Mesh%Nfacette + 1
            Mesh%Tfacette(Mesh%Nfacette)%Tnoeud = Nnoeud + [j4,j2,j3]
        end do
    end do
    Mesh%Body(1)%IndBody(3:4) = [Mesh%Nnoeud, Mesh%Nfacette]
    Mesh%Tnoeud(Mesh%Body(1)%IndBody(1):Mesh%Body(1)%IndBody(3))%typeNoeud = 1
    Mesh%Tfacette(Mesh%Body(1)%IndBody(2):Mesh%Body(1)%IndBody(4))%typeFrontiere = 1
    Mesh%Tnoeud(Mesh%Body(1)%IndBody(1):Mesh%Body(1)%IndBody(3))%Npanneau = 1
    Mesh%Tfacette(Mesh%Body(1)%IndBody(2):Mesh%Body(1)%IndBody(4))%Npanneau = 1

    ! Wigley Hull mesh
    L = InputData%Lgeom(1,1)
    Nx2 = 0.5*L/Mesh%DimTank(5)
    if(symmetry)then
        size_Rwh_Thw = 2*Nx2+1
        allocate(Rwh(1:2*Nx2+1), Twh(1:2*Nx2+1))
    else
        size_Rwh_Thw = 4*Nx2
        allocate(Rwh(1:4*Nx2), Twh(1:4*Nx2))
    end if
    call WigleyHull(Mesh, Rwh, Twh, .false.,InputData,size_Rwh_Thw)
    
    if (cuve_ferme) then
        Mesh%Nsys = Mesh%Body(Mesh%NBody)%IndBody(3)
        Mesh%Nfsys = Mesh%Body(Mesh%NBody)%IndBody(4)
    else
        if (is_body.and.Int_Body.ne.0) then
            Mesh%Nsys = Mesh%Body(Int_Body)%IndBody(3)
            Mesh%Nfsys = Mesh%Body(Int_Body)%IndBody(4)
        else
            Mesh%Nsys = Mesh%FS%IndFS(3)
            Mesh%Nfsys = Mesh%FS%IndFS(4)
        end if
    end if
    call GeomInit(Mesh, fgeom_vect, 0._RP,InputData, .false.)
    if (allocated(Xtab)) deallocate(Xtab)
    if (allocated(Ytab)) deallocate(Ytab)
    if (allocated(Ztab)) deallocate(Ztab)
    if (allocated(Rwh)) deallocate(Rwh)
    if (allocated(Twh)) deallocate(Twh)
    
end subroutine Mesh_cart_Wigley
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!   Reconstruction d'un maillage à partir d'un maillage initial et de la surface libre incidente   !
!                                 LETOURNEL Lucas    03/13                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine Remaillage(t, Mesh0, Mesh, InputData, ierror)
!    use BodyMotion_mod
!
!    !Paramètres
!    !character(len=50) :: filename
!    !f2py integer*1, dimension(1000) :: Mesh0,Mesh
!    type(TMaillage) :: Mesh0, Mesh
!    real(rp) :: t
!    !f2py integer*1, dimension(1000)    :: InputData
!    type(InputDataStruct),intent(inout) :: InputData                            ! Input data
!    integer, intent(inout) :: ierror
!
!    ! Variables Locales
!    integer :: j, k, lf, nc !, Interf
!    !f2py integer*1, dimension(1000) :: MeshTemp
!    type(TMaillage) :: MeshTemp
!    real(rp) :: Eta0,dlx
!    integer :: i1,i2,i3
!    real(rp), dimension(3) :: vect_product_1
!    real(rp) :: Crampe_1,Crampe_2
!
!    ierror = 0
!    dlx = 0._RP
!
!    if (.not.DeformMesh) then
!        call CopyMaillage(Mesh,Mesh0)
!        call GeomInit(Mesh, fgeom_vect, t, InputData, .false.)
!    else
!        call CopyMaillage(Mesh, Mesh0)
!    
!        do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
!            call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
!            Mesh%Tnoeud(j)%Pnoeud(3) = Eta0
!        end do
!        do nc=1,Mesh%Nbody
!            if (nc.ne.Int_Body) then
!                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
!                    call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
!                    call Computation_Crampe(-Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Mesh%Body(nc)%DimBody(3),Crampe_1)
!                    if (Crampe_1.lt.0.1_RP) then
!                        Mesh%Tnoeud(j)%Pnoeud(3) = Mesh%Tnoeud(j)%Pnoeud(3) + Eta0
!                    else
!                        call Computation_Crampe(-Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Mesh%Body(nc)%DimBody(3),Crampe_2)
!                        Mesh%Tnoeud(j)%Pnoeud(3) = Mesh%Tnoeud(j)%Pnoeud(3) + Eta0*(1._RP-Crampe_2)
!                    end if
!                end do
!            end if
!        end do
!        if (is_body.and.free_body) then
!            Mesh%Body(Int_Body)%GBody(1:3) = Mesh0%Body(Int_Body)%GBody(1:3) + Mesh%Body(Int_Body)%VBody(1:3)*dt
!            do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%body(Int_Body)%IndBody(3)
!                call Computation_vect_product(Mesh%Body(nc)%MBody(4:6),Mesh0%Tnoeud(j)%Pnoeud(1:3)-Mesh0%Body(nc)%GBody(1:3),vect_product_1)
!                Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud + Mesh%Body(nc)%MBody(1:3) + vect_product_1
!            end do
!        end if
!
!        call GeomInit(Mesh, t, .false.)   
!
!    end if
!
!    9999 continue
!      if(ierror/=0)then 
!        write(*,99),ierror
!      endif
!    99 format('** error #',i3,' : pb. remaillage')
!
!      select case(ierror)
!      case(200)
!        i1 = MeshTemp%Tfacette(k)%Tnoeud(1)
!        i2 = MeshTemp%Tfacette(k)%Tnoeud(2)
!        i3 = MeshTemp%Tfacette(k)%Tnoeud(3)
!        print*,'** ifacette = ',k,'  Noeud_1 : ',i1,'  Noeud_2 : ',i2,'  Noeud_3 : ',i3
!        print*,'** lf = ',lf
!      case default
!      end select
!
!    contains
!
!    subroutine Computation_Interf(t,Mesh,indice,Interf)
!
!    !Paramètre
!    !f2py integer*1, dimension(1000) :: Mesh
!    type(TMaillage),intent(in) :: Mesh
!    integer,intent(in) :: indice
!    real(rp),intent(in) :: t
!    integer,intent(out) :: Interf
!
!    !Variables locales
!    character (len=50) :: filecoeff
!    integer :: j, k
!    real(rp) Eta0
!    real(rp), dimension(3) :: M
!    Interf=0
!    do j=1,Mesh%Tnoeud(indice)%Nfacette
!        do k=1,3
!            M = Mesh%Tnoeud(Mesh%Tfacette(Mesh%Tnoeud(indice)%Tfacette(j,1))%Tnoeud(k))%Pnoeud
!            call CEta0( M, t, Eta0)
!            if ((Eta0-M(3)).gt.Epsilon2) then ! Valeur de 0.01 permet de ne pas considérer des points trop proches et d'éviter les triangles plats
!                Interf=1
!            end if
!        end do
!    end do
!    end subroutine Computation_Interf
!
!end subroutine Remaillage


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                               Extraction d'un maillage depuis un fichier                         !
!                                       LETOURNEL Lucas 14/06/12                                   !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extract_maillage(filemesh,fgeom_vect,Mesh,InputData)
    !!!!! Problème :
    !   Extraire le maillage d'un fichier texte, contenant :
    !       o Le nom du maillage sur la première ligne (100 caractère maxi)
    !       o Le nombre de points et de facettes
    !       o Les dimensions de la cuve: (optionnel)
    !           + L
    !           + LL
    !           + Profondeur
    !       o Indice puis coordonnés du point 
    !       o Indice de la facette puis indice des sommets puis type de frontière (0--> SL, autre --> SM)
    !Ex:
    !mailllarge
    !Nombre de Points =    558	 Nombre de Facettes =    864	 
    !dimensions de la cuve : 
    !L = 5.000000
    !l = 0.500000
    !h = 2.000000
    !1	0.000000E+00	0.000000E+00	0.000000E+00	
    !2	0.000000E+00	2.500000E-01	0.000000E+00	
    !3	0.000000E+00	5.000000E-01	0.000000E+00
    !...
    !1	5	2	1	0	
    !2	1	4	5	0	
    !3	6	3	2	0	
    !4	2	5	6	0
    !!!!!

    character(len=50)                   :: filemesh
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                       ! Geometries.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                        ! Input data.

    integer                             ::  j, nf, Nnoeud, Nfacette, ios
    integer                             :: NPFS, NPBD, NFFS
    integer                             :: NFBD
    !real(rp), dimension(3)             :: te
    real(rp), dimension(11)             :: lignep
    integer, dimension(3)               :: lignef
    character(len=100)                  :: truc
    
    ! This subroutine extracts a mesh from filemesh.
    
    print*,"Extract_maillage: this subroutine is outdated."
    call exit()
    print*,"Extract_maillage: this subroutine is not MB."
    pause
    
    ! Création d'un objet Mesh
    call NewMaillage(Mesh,Pointmax,NBodies+1) ! +1 for the tank.

    ! Ouverture du fichier de maillage
    open(10,file=filemesh, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier de maillage"

    ! Lecture des informations du maillage
    read(10,fmt='(A100)') truc
    read(10,fmt='(A100)') truc
    read(10,fmt='(A100)') truc

    read(truc(9:20),*) Nnoeud
    read(truc(25:37),*) Nfacette
 
    Mesh%Nnoeud=Nnoeud
    Mesh%Nfacette=Nfacette

    NPFS = 0
    NPBD = 0
    NFFS = 0
    NFBD = 0

    ! Lecture des noeuds du maillage
    do j=1,Mesh%Nnoeud
        read(10,*) lignep
        !np=lignep(1)
        Mesh%Tnoeud(j)%Pnoeud=[lignep(1), lignep(2), lignep(3)] 
        Mesh%Tnoeud(j)%Normale = [lignep(4),lignep(5),lignep(6)]   
        Mesh%Tnoeud(j)%typeNoeud = int(lignep(7))
        Mesh%Tnoeud(j)%Npanneau = int(lignep(8))
        Mesh%Tnoeud(j)%Damping = lignep(9)
        Mesh%Tnoeud(j)%NVoisin = int(lignep(10))
        Mesh%Tnoeud(j)%Ndouble = int(lignep(11))
        if(int(lignep(8)).eq.0)then
            NPFS = NPFS + 1
        elseif(int(lignep(8)).eq.-9)then
            NPBD = NPBD+ 1
        else
            print*,'Erreur Npanneau'
        endif    
    end do
    ! Lecture des facettes du maillage
    do j=1,Mesh%Nfacette
        read(10,*) lignef
        nf=lignef(1)
        ! Indices des sommets de la facette
        Mesh%Tfacette(j)%Tnoeud=[lignef(1), lignef(2), lignef(3)]
        Mesh%Tfacette(j)%Npanneau = Mesh%Tnoeud(lignef(1))%Npanneau
        Mesh%Tfacette(j)%typeFrontiere = Mesh%Tnoeud(lignef(1))%typeNoeud
        if(Mesh%Tfacette(j)%Npanneau == 0)then
            NFFS = NFFS + 1
        elseif(Mesh%Tfacette(j)%Npanneau .eq.-9)then
            NFBD = NFBD + 1
        endif
    end do
    close(10)

    Mesh%FS%IndFS(1:4) = [1,1,NPFS,NFFS]
    Mesh%Body(Int_Body)%IndBody(1:4) = [NPFS+1,NFFS+1,NPFS+NPBD,NFFS+NFBD]

    call geomInit(Mesh, fgeom_vect, 0._RP,InputData)
    print*,"Extract_maillage: t0 is fixed to 0."

    Nnoeud = Mesh%Nnoeud
    Nfacette = Mesh%Nfacette

    print*, '   Extraction du maillage ok'

end subroutine extract_maillage

! ----------------------------------------------------------------------
! Update_corres_mesh
! ----------------------------------------------------------------------

subroutine update_corres_mesh(tab1,tab2,n1,tab_out,ierror)
  implicit none
  integer,dimension(:),intent(in) :: tab1,tab2
  integer,intent(in) :: n1
  integer,dimension(:),intent(inout) :: tab_out
  integer,intent(inout) :: ierror
! local
  integer :: j,j2
  
  ierror = 0
  tab_out(:) = 0
  
  do j=1,n1
    if(tab1(j).gt.0 .and. tab2(j).gt.0)then
      j2 = tab2(j)
      tab_out(j2) = tab1(j)  
    endif
  enddo

end subroutine update_corres_mesh
	
subroutine extract_mesh2bis(Mesh,InputData,NumBody)

    ! Nécessité d'avoir les panneau les uns à la suite des autres -> Eventuellement obsolète (pas testé)
    ! Géométrie convexe imposée
    ! Necessity to have the different pannels one after the other. -> maybe not true (not tested)
    ! Convex geometry imposed

    ! File format :
    ! int  int              (unused, symmetry)
    ! int  real real real   (NodeIndex, NodePosition)
    ! ...  ...  ...  ...
    ! int  real real real
    ! 0    0.   0.   0.
    ! int  int  int  int    (node index of each quadrilateral face)
    ! ...  ...  ...  ...
    ! int  int  int  int
    ! 0    0    0    0
    
    ! Paramètres
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage),intent(inout)           :: Mesh         ! Mesh
    !f2py integer*1, dimension(1000)        :: InputData
    type(InputDataStruct),intent(in)        :: InputData    ! Input data
    integer,intent(in)                      :: NumBody      ! Body number
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Meshtemp

    ! Loops index
    integer ::  j, k,l, m, ind
    ! File opening
    integer :: ios
    ! Symmetry (mesh input)
    logical :: sym
    ! Node Index (while reading)
    integer :: NodeInd
    ! Counting integer (in loops)
    integer :: compt, comptf, comptNpanneau
    integer ::  NNoeudCorps,NFacetteCorps
    ! Reading
    integer, dimension(2) :: FstLine, IndBody   ! First line
    real(rp), dimension(3) :: NodeLine          ! Nodes
    integer, dimension(4) :: FLine              ! Faces
    ! Vector norm2
    real(rp) :: norme
    ! Tables : repartition nodes and faces by pannels
    integer, dimension(100,FacetteMax) :: FPan
    integer, dimension(100,FacetteMax) :: NPan
    ! Sizes of precedent fields of each pannels
    integer, dimension(100):: FpanCompt,NPanCompt
    ! Normals for differnt pannels
    real(rp),dimension(3)::Normalej, Normalek, vp
    ! Vectors : PiM where m is in a pannel and (Pi )i=1,3 are the node of a face of a different pannel
    real(rp),dimension(3,3) :: v
    ! Node tables of faces
    integer,dimension(3) :: TNk,TNj
    ! Idem but whith 2 faces
    integer,dimension(2,3) :: vnoeud
    ! 3* Center of 2 faces (sum of positions)
    real(rp),dimension(2,3) :: CtreFacette
    ! Table pannel i and j are coplanar <=> coplaniare(i,j) = 1
    integer,dimension(100,100) :: coplanaire
    ! Copy(j) = true <=> Node i has already been copied
    integer,dimension(Pointmax) :: copy
    ! Continue loop or not
    logical :: cont
    real(rp),dimension(3) :: vect_product_1,vect_product_2,vect_product_3,vect_product_4,vect_product_5
    real(rp),dimension(3) :: HorizontalPosition

    ! This subroutine extracts the mesh of the file nemoh.dat.
    
    print*,"extract_mesh2bis: this subroutine is not MB."
    
    ! Horizontal position
    HorizontalPosition = InputData%Position(1:3,1,NumBody)
    HorizontalPosition(3) = 0._RP
    print*,"extract_mesh2bis: The vertical position is the one of Nemoh.dat."
    
    !Use of a temporary Mesh structure
    call NewMaillage(MeshTemp,PointMax,NBodies+1) ! +1 for the tank.
    
    MeshTemp=Mesh
    
    ! Open mesh file: nemoh.dat
    open(10,file=filenemoh, iostat=ios)
    if (ios/=0) stop "Erreur à l'ouverture du fichier de maillage"
    
    ! Lecture premiere ligne (symetrie)
    ! Reading first line (symmetry)
    read(10,*) FstLine
    if(FstLine(2).eq.1) then
        sym = .true.
        if (.not. Symmetry) then
            print*, 'Please use a complete BodyMesh in input or allow symmetrys'
        end if
    else
        sym = .false.
    end if
    
    ! Reading Nodes -> in MeshTemp%Tnoeud
    ! Lecture des noeuds stockés dans tableau Meshtemp%Tnoeud
    IndBody=Mesh%Body(Mesh%NBody)%IndBody(1:2)
    compt = IndBody(1)-1
    do
        read(10,*) NodeInd, NodeLine(1:3)
        ! Stop reading nodes when arriving at the zero line
        if(NodeInd.eq.0) then
            NNoeudCorps = compt - IndBody(1)+1
            MeshTemp%Nnoeud = compt
            MeshTemp%Body(Mesh%NBody)%IndBody(3) = compt
            exit
        endif
        compt=compt+1
        MeshTemp%TNoeud(compt)%Pnoeud=NodeLine(1:3) + HorizontalPosition
        MeshTemp%TNoeud(compt)%TypeNoeud=1
        MeshTemp%Tnoeud(compt)%Nfacette=0
        MeshTemp%Tnoeud(compt)%Npanneau=0
    end do

    ! Reading faces -> in MeshTemp%Tfacette
    ! Lecture des facettes stockées dans MeshTemp%TFacette
    comptNpanneau=0
    compt=IndBody(2)-1
    FPancompt=0
    NPancompt=0
    Npan=0
    Fpan =0
    do
        !Lecture fichier d'entrée : quadrangle
        read(10,*) FLine
        ! Stop reading faces when arriving at the zero line
        if(FLine(1).eq.0) then
            MeshTemp%Nfacette= compt
            NFacetteCorps = compt-IndBody(2)+1
            MeshTemp%Body(Mesh%NBody)%IndBody(4) = compt
            exit
        end if
        cont=.true.
        ! Not taking into account the nodes that have y<0 if symmetry is on
        if(sym) then
            do j=1,4
                if(MeshTemp%Tnoeud(Fline(j)+IndBody(1)-1)%Pnoeud(2).lt.-epsilon) then
                    cont = .false.
                endif
            enddo
        endif
       if(cont) then
        !Séparation quadrangle en deux triangles
        !Prise en compte premier triangle (noeuds 1,2,3)
        !Vérification facette non plate
            call Computation_vect_product(MeshTemp%Tnoeud(Fline(3)+IndBody(1)-1)%Pnoeud - MeshTemp%Tnoeud(Fline(1)+IndBody(1)-1)%Pnoeud , &
                                  MeshTemp%Tnoeud(Fline(3)+IndBody(1)-1)%Pnoeud - MeshTemp%Tnoeud(Fline(2)+IndBody(1)-1)%Pnoeud,vect_product_1)
            norme = norm2(vect_product_1)
            if(norme.gt.epsilon ) then
                ! Parameters of the face
                compt=compt+1
                MeshTemp%Tfacette(compt)%Tnoeud=FLine(1:3)+IndBody(1)-1
                MeshTemp%Tfacette(compt)%TypeFrontiere = 1
                MeshTemp%Tfacette(compt)%Npanneau = 0
                do j=1,3
                    !Parameters of the nodes
                    MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette=MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette+1
                    MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Tfacette(MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette,1)=compt

                    ! Pannel repartition for faces and nodes

                    !Répartition par panneau
                    !Si un des noeuds possède déjà un numéro de panneau...
                    if(MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau.ne.0) then
                        !...la facette est sur le même panneau
                        MeshTemp%Tfacette(compt)%Npanneau =  MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau
                    endif
                end do
                if(MeshTemp%Tfacette(compt)%Npanneau.eq.0) then !Sinon (aucun des noeuds n'a déjà été visité)
                    !On attribue à la facette le numéro de panneau suivant
                    comptNpanneau=comptNpanneau+1
                    MeshTemp%Tfacette(compt)%Npanneau = comptNPanneau
                endif
                do j=1,3
                    if( MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau.ne.MeshTemp%Tfacette(compt)%Npanneau) then
                        Npancompt(MeshTemp%Tfacette(compt)%Npanneau) = NPanCompt(MeshTemp%Tfacette(compt)%Npanneau) + 1
                        Npan(MeshTemp%Tfacette(compt)%Npanneau,NPanCompt(MeshTemp%Tfacette(compt)%Npanneau))=MeshTemp%Tfacette(compt)%Tnoeud(j)
                        MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau = MeshTemp%Tfacette(compt)%Npanneau
                    endif
                enddo
                Fpancompt(MeshTemp%Tfacette(compt)%Npanneau) = Fpancompt(MeshTemp%Tfacette(compt)%Npanneau) +1
                Fpan(MeshTemp%Tfacette(compt)%Npanneau,Fpancompt(MeshTemp%Tfacette(compt)%Npanneau)) = compt
                ! End of pannel repartition
            else
                print*, 'Facette plate vérifier maillage'
                print*, norme
                print*, Fline(1),FLine(2),FLine(3)
            end if

            !Idem nodes 4, 3 et 1
            call Computation_vect_product(MeshTemp%Tnoeud(Fline(4)+IndBody(1)-1)%Pnoeud &
                                  - MeshTemp%Tnoeud(Fline(3)+IndBody(1)-1)%Pnoeud , &
                                    MeshTemp%Tnoeud(Fline(1)+IndBody(1)-1)%Pnoeud   &
                                  - MeshTemp%Tnoeud(Fline(3)+IndBody(1)-1)%Pnoeud,vect_product_2)
            norme = norm2(vect_product_2)
            if(norme.gt.epsilon ) then
                compt=compt+1
                MeshTemp%Tfacette(compt)%Tnoeud(1)=FLine(3)+IndBody(1)-1
                MeshTemp%Tfacette(compt)%Tnoeud(2)=FLine(4)+IndBody(1)-1
                MeshTemp%Tfacette(compt)%Tnoeud(3)=FLine(1)+IndBody(1)-1
                do j=1,3
                    MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette=MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette+1
                    MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Tfacette(MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Nfacette,1)=compt
                    if(MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau.ne.0) then
                        MeshTemp%Tfacette(compt)%Npanneau =  MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau
                    endif
                end do
                if(MeshTemp%Tfacette(compt)%Npanneau.eq.0) then
                    comptNpanneau=comptNpanneau+1
                    MeshTemp%Tfacette(compt)%Npanneau = comptNPanneau
                    Fpan(comptNpanneau,1)=compt
                    Fpancompt(comptNPanneau)=1
                endif
                do j=1,3
                    if( MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau.ne.MeshTemp%Tfacette(compt)%Npanneau) then
                        Npancompt(MeshTemp%Tfacette(compt)%Npanneau) = NPanCompt(MeshTemp%Tfacette(compt)%Npanneau) + 1
                        Npan(MeshTemp%Tfacette(compt)%Npanneau,NPanCompt(MeshTemp%Tfacette(compt)%Npanneau))=MeshTemp%Tfacette(compt)%Tnoeud(j)
                        MeshTemp%Tnoeud(MeshTemp%Tfacette(compt)%Tnoeud(j))%Npanneau = MeshTemp%Tfacette(compt)%Npanneau
                    endif
                enddo
                Fpancompt(MeshTemp%Tfacette(compt)%Npanneau) = Fpancompt(MeshTemp%Tfacette(compt)%Npanneau) +1
                Fpan(MeshTemp%Tfacette(compt)%Npanneau,Fpancompt(MeshTemp%Tfacette(compt)%Npanneau)) = compt
            else
                print*, 'Facette plate vérifier maillage'
                print*, 'Noeuds'
                print*,'NNoeudCorps',NNoeudCorps
                print*, Fline(1),FLine(3),FLine(4)
                print*, Fline(1)+IndBody(1)-1,FLine(3)+IndBody(1)-1,FLine(4)+IndBody(1)-1
                print*,MeshTemp%Tnoeud(Fline(4)+IndBody(1)-1)%Pnoeud
                print*,MeshTemp%Tnoeud(Fline(3)+IndBody(1)-1)%Pnoeud
                print*,MeshTemp%Tnoeud(Fline(1)+IndBody(1)-1)%Pnoeud
            end if
        endif
    end do
    
    !Test coplanarité : panneau i et j coplanaire ssi normales aux facettes colinéaires et
    coplanaire=0
    do j=1,comptNpanneau-1
        TNj=MeshTemp%Tfacette(Fpan(j,1))%Tnoeud
        call Computation_vect_product(MeshTemp%Tnoeud(TNj(2))%Pnoeud-MeshTemp%Tnoeud(TNj(1))%Pnoeud,MeshTemp%Tnoeud(TNj(3))%Pnoeud-MeshTemp%Tnoeud(TNj(1))%Pnoeud,Normalej)
        do k=j+1,comptNpanneau
            TNk=MeshTemp%Tfacette(Fpan(k,1))%Tnoeud
            call Computation_vect_product(MeshTemp%Tnoeud(TNk(2))%Pnoeud-MeshTemp%Tnoeud(TNk(1))%Pnoeud,MeshTemp%Tnoeud(TNk(3))%Pnoeud-MeshTemp%Tnoeud(TNk(1))%Pnoeud,Normalek)
            call Computation_vect_product(Normalej, Normalek,vp)
            if(norm2(vp).lt.epsilon)then
                vp = MeshTemp%Tnoeud(TNj(1))%Pnoeud-MeshTemp%Tnoeud(TNk(1))%Pnoeud
                if(abs(dot_product(vp, Normalej)).lt.epsilon) then
                    print*,'Panneau ',j,' et ',k,' coplanaires'
                    coplanaire(j,k)=1
                    coplanaire(k,j)=1
                endif
            endif
        enddo
    enddo
    
    ! Dealing with non-conformities :
    ! Common nodes for different pannels
    ! Node in an edge of another pannel
    copy=0
    compt=Mesh%Body(Mesh%NBody)%IndBody(1)-1
    comptf= Mesh%Body(Mesh%NBody)%IndBody(2)-1
    do j=1,comptNPanneau
        !Recherche des exceptions
            do k=1,comptNpanneau
                l=0
                do while(l.lt.Fpancompt(j))
                l=l+1
                do m=1,NPancompt(k)  ! Parcours sur les noeuds de k
                    cont=.true.
                    ! Si noeud commun on les attributs du noeud de j sont copiés dans celui de k
                    ! Le noeud du panneau j ne sera pas recopié par la suite
                            
                    do ind=1,3
                        v(ind,:) = MeshTemp%Tnoeud(MeshTemp%Tfacette(FPan(j,l))%Tnoeud(ind))%Pnoeud-MeshTemp%Tnoeud(Npan(k,m))%Pnoeud
                        if(norm2(v(ind,:)).lt.epsilon) then !Si noeud commun on les attributs du noeud de j sont copiés dans celui de k
                            if(k.gt.j)then
                                if(coplanaire(k,j).eq.1)then
!                                    print*, 'Pnoeud',MeshTemp%Tnoeud(Npan(k,m))%Pnoeud
!                                    print*, 'FPan : ',FPan(j,l), ' NPan : ',NPan(k,m)
                                MeshTemp%Tfacette(FPan(j,l))%Tnoeud(ind)= NPan(k,m)
                                MeshTemp%Tnoeud(NPan(k,m))%Nfacette=MeshTemp%Tnoeud(NPan(k,m))%Nfacette+1
                                MeshTemp%Tnoeud(Npan(k,m))%Tfacette(MeshTemp%Tnoeud(NPan(k,m))%Nfacette,1)=Fpan(j,l)
                                endif
                            endif
                            cont=.false.
                        endif
                    enddo
                            
                    ! Sinon on teste si le noeud m du panneau k est situé sur une arrete de la facette l du panneau j
                    ! Dans ce cas la facette est divisée en deux facettes, toujours orientées dans le sens direct
                    if(cont) then
                        call Computation_vect_product(v(1,:),v(2,:),vect_product_3)
                        if(norm2(vect_product_3).lt.epsilon) then
                            if(dot_product(v(1,:),v(2,:)).lt.epsilon) then
                                if(coplanaire(j,k).eq.0) then
                                    MeshTemp%Nnoeud = MeshTemp%Nnoeud +1
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Pnoeud = MeshTemp%Tnoeud(NPan(k,m))%Pnoeud
                                    ind = MeshTemp%Nnoeud
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Npanneau = j

                                else
                                    ind = NPan(k,m)
                                endif
                                MeshTemp%Nfacette=MeshTemp%Nfacette+1
                                MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(1) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1)
                                MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(2) = ind
                                MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(3) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3)
                                MeshTemp%Tfacette(MeshTemp%Nfacette)%Npanneau = j
                                MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1) = ind
                                MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2)
                                MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3)
                                Fpancompt(j)=Fpancompt(j)+1
                                Fpan(j,Fpancompt(j))=MeshTemp%Nfacette
                                cont=.false.
                            endif
                        endif
                    endif
                    if(cont) then
                        call Computation_vect_product(v(1,:),v(3,:),vect_product_4)
                        if(norm2(vect_product_4).lt.epsilon) then
                            if(dot_product(v(1,:),v(3,:)).lt.epsilon) then
                                if(coplanaire(j,k).eq.0) then
                                    MeshTemp%Nnoeud = MeshTemp%Nnoeud +1
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Pnoeud = MeshTemp%Tnoeud(NPan(k,m))%Pnoeud
                                    ind = MeshTemp%Nnoeud
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Npanneau = j
                                else
                                    ind = NPan(k,m)
                                endif
                            MeshTemp%Nfacette=MeshTemp%Nfacette+1
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(1) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1)
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(2) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2)
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(3) = ind
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Npanneau = j
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1) = ind
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2)
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3)
                            Fpancompt(j)=Fpancompt(j)+1
                            Fpan(j,Fpancompt(j))=MeshTemp%Nfacette
                            cont=.false.
                            endif
                        endif
                    endif
                    if(cont) then
                        call Computation_vect_product(v(2,:),v(3,:),vect_product_5)
                        if(norm2(vect_product_5).lt.epsilon) then
                            if(dot_product(v(2,:),v(3,:)).lt.epsilon) then
                                if(coplanaire(j,k).eq.0) then
                                    MeshTemp%Nnoeud = MeshTemp%Nnoeud +1
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Pnoeud = MeshTemp%Tnoeud(NPan(k,m))%Pnoeud
                                    ind = MeshTemp%Nnoeud
                                    MeshTemp%Tnoeud(MeshTemp%Nnoeud)%Npanneau = j
                                else
                                    ind = NPan(k,m)
                                endif
                            MeshTemp%Nfacette=MeshTemp%Nfacette+1
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(1) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1)
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(2) = ind
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Tnoeud(3) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3)
                            MeshTemp%Tfacette(MeshTemp%Nfacette)%Npanneau = j
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(1)
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2) = MeshTemp%Tfacette(FPan(j,l))%Tnoeud(2)
                            MeshTemp%Tfacette(FPan(j,l))%Tnoeud(3) = ind
                            Fpancompt(j)=Fpancompt(j)+1
                            Fpan(j,Fpancompt(j))=MeshTemp%Nfacette
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
    
    ! Other exceptions : faces written twice
    ! Elimitnation des facettes doubles
    do j=1,comptNpanneau
        do k =j, comptNpanneau
            l=1
            do while(l.le.Fpancompt(j))
                m=1
                do while(m.le.Fpancompt(k))
                    if(k.ne.j .or.(m.gt.l)) then
                    vnoeud(1,:) = MeshTemp%Tfacette(Fpan(j,l))%Tnoeud
                    vnoeud(2,:) = MeshTemp%Tfacette(Fpan(k,m))%Tnoeud
                    CtreFacette(1,:) = 0
                    CtreFacette(2,:) = 0
                    do ind=1,3
                        CtreFacette(1,:) = Ctrefacette(1,:) + MeshTemp%Tnoeud(vnoeud(1,ind))%Pnoeud
                        CtreFacette(2,:) = Ctrefacette(2,:) + MeshTemp%Tnoeud(vnoeud(2,ind))%Pnoeud
                    enddo
                    if(norm2(CtreFacette(1,:)-Ctrefacette(2,:)).lt.epsilon) then
                        Fpan(k,m:Fpancompt(k)-1)=Fpan(k, m+1:Fpancompt(k))
                        Fpancompt(k)=Fpancompt(k)-1
                        print*, 'Facettes identiques'
                    else
                        m=m+1
                    endif
                    else
                        m = m+1
                    endif
                enddo
                l=l+1
           enddo
       enddo
    enddo
    
    ! Copy of MeshTemp in Mesh
    ! Loop on the faces -> Common Nodes are not copied
        
    ! Copie du maillage MeshTemp dans Mesh à partir des facettes
    ! Parcours des facettes
    do j=1,comptNpanneau
        do k=1,Fpancompt(j)
            comptf=comptf+1
            Mesh%Tfacette(comptf)= MeshTemp%Tfacette(Fpan(j,k))
            do l=1,3
                if(copy(MeshTemp%Tfacette(Fpan(j,k))%Tnoeud(l)).eq.0) then
                    !Si le sommet n'a pas encore été recopié
                    compt=compt+1
                    !On le note comme recopié dans le noeud compt
                    copy(MeshTemp%Tfacette(Fpan(j,k))%Tnoeud(l))=compt
                    !La facette comptf  à comme sommet ce noeud compt
                    Mesh%Tfacette(comptf)%Tnoeud(l)=compt
                    ! Le noeud est recopié dans le noeud compt
                    Mesh%Tnoeud(compt)=MeshTemp%Tnoeud(MeshTemp%Tfacette(Fpan(j,k))%Tnoeud(l))
                    Mesh%Tnoeud(compt)%Nfacette=1
                    Mesh%Tnoeud(compt)%Tfacette(1,1)= comptf
                else ! Si il a déjà été recopié, onlui affect la facette comptf et vice versa
                    ind = copy(MeshTemp%Tfacette(Fpan(j,k))%Tnoeud(l))
                    Mesh%Tfacette(comptf)%Tnoeud(l) = ind
                    Mesh%Tnoeud(ind)%Nfacette=Mesh%Tnoeud(ind)%Nfacette+1
                    Mesh%Tnoeud(ind)%Tfacette(Mesh%Tnoeud(ind)%Nfacette,1)=comptf
                endif
            enddo
        enddo
    enddo
    !Mesh%Body(Mesh%NBody)%IndBody(3)=compt
    Mesh%NNoeud=compt
    !Mesh%Body(Mesh%NBody)%IndBody(4)=comptf
    Mesh%Nfacette=comptf
    Mesh%Body(Mesh%NBody)%IndBody(3)=MeshTemp%Body(Mesh%NBody)%IndBody(3)
    Mesh%Body(Mesh%NBody)%IndBody(4)=MeshTemp%Body(Mesh%NBody)%IndBody(4)
    print*, 'IndBody', Mesh%Body(Mesh%NBody)%IndBody
    print*, 'fin extract mesh'
    
    call DelMaillage(MeshTemp)
    close(10)
    
end subroutine extract_mesh2bis


end module GenMaillage