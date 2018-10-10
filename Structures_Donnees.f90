module Structuresdonnees
use Constantes
use iso_c_binding
implicit none

! type liste d'ecoulements (parareal)

type T_liste_ecoulements
    type(TEcoulement), dimension(:,:), allocatable :: G ! Liste d'ecoulements obtenus pas la methode grossiere
    type(TEcoulement), dimension(:), allocatable :: F ! Liste d'ecoulements obtenus pas la methode fine
    type(TEcoulement), dimension(:,:), allocatable :: lambda ! Liste d'ecoulements correspondants aux differents lambda
end type T_liste_ecoulements

type T_liste_maillages
    type(TMaillage), dimension(:,:), allocatable :: G ! Liste des maillages obtenus pas la methode grossiere
    type(TMaillage), dimension(:), allocatable :: F ! Liste des maillages obtenus pas la methode fine
end type T_liste_maillages


type T_liste_bodies 
    type(TBody), dimension(:,:,:), allocatable :: G ! Liste de corps obtenus pas la methode grossiere
    type(TBody), dimension(:,:), allocatable :: F ! Liste de corps obtenus pas la methode fine
    type(TBody), dimension(:,:,:), allocatable :: lambda ! Liste de corps correspondants aux differents lambda
end type T_liste_bodies


!	type noeud
type TNoeud
    real(rp),dimension(3)		:: Pnoeud			! Position des noeuds dans le repère global.
    integer						:: Nfacette			! Nombre de facettes dont fait parti le noeud.
    integer, dimension(Ncnxmx,2):: Tfacette			! Table des facettes connectées au noeud, numéro de sommet du noeud et angle au sommet.
    integer						:: TypeNoeud      	! type de noeud (0 SL, 1 SM, 10 intersection SL/SM, 11 intersection SM/SM).
    real(rp), dimension(3)	    :: Normale          ! Normale au noeud.
    real(rp), dimension(3,3,2)  :: Plocal           ! Base Locale (1) et sa transposée(2) ((1:3,1,1) = u, (1:3,2,1) = v, (1:3,3,1) = normal).
    real(rp), dimension(5)      :: DLocal           ! Dérivées locales selon les vecteurs tangents de la base locale, premières (1 et 2) et secondes (3 et 4) et delta.
    real(rp),dimension(Ncnxmx)  :: Angle            ! Angle au sommet pour chaque facette auquel est relié le noeud.
    real(rp)                    :: Aire             ! Aire associée au noeud.
    integer                     :: Ndouble          ! Nombre de noeuds doubles.
    integer, dimension(3)       :: double           ! Table des noeuds doubles.
    integer                     :: TypeDouble       ! = 0 if the normal is (0,0,-1) = 1 otherwise. -> useless ?
    integer, dimension(2)       :: Nvoisin          ! Nombre de Voisins (1er : N facettes Voisines, 2e : N Points Voisins(1st order and 2nd order)).
    integer                     :: Ordre            ! Ordre des BSplines à utiliser.
    real(rp)                    :: Rvoisin          ! Critère de distance poure la définition du voisinage.
    integer, dimension(NPVMax,2):: Tvoisin          ! Table des noeuds voisins NVoisin(2)x(indice, ordre de voisinage).
    real(rp)                    :: Damping          ! Coefficient d'amortissement du noeud (0< <1).
    integer                     :: NPanneau         ! Type de structure (0 = free surface, 1 = tank, 2 and more = floaters).
    logical                     :: mobility         ! Mobility of the node in the local coordinate.
    logical                     :: control_point    ! Definition des points de controle.
    real(rp), dimension(3)      :: Velocity         ! Velocity of the node of the mesh.
    integer                     :: indmesh          ! indice of the corresponding node in mesh.
end type TNoeud

!	type Facette
type TFacette
    integer, dimension(3)		:: Tnoeud			! Indexes of the 3 nodes/vertexes of the panel
    real(rp), dimension(3,3)    :: PNoeud           ! Positions des sommets de la facette
    real(rp), dimension(3)      :: Gfacette         ! Centre de gravité de la facette
    real(rp)                    :: Rmax             ! Distance maximale du centre de gravité aux sommets
    real(rp)                    :: Aire             ! Aire de la facette
    real(rp), dimension(3)		:: Normale			! Normale à la facette
    real(rp), dimension(3,3)    :: ds               ! Matrice pour le Gradient Surfacique de la facette
    real(rp), dimension(3,3)    :: dsT              ! Matrice transposée pour le Gradient Surfacique de la facette
    integer						:: TypeFrontiere	! type de frontière du domaine (0. Surface libre, 1. Solide)
    integer                     :: NPanneau         ! Type de structure (0 = free surface, 1 = tank, 2 and more = floaters)
end type TFacette

type TBody
    logical                             :: is_tank  ! Body is an external boundary (numerical tank wall).
    integer, dimension(4)               :: IndBody  ! Indices of the 1st and last Nodes/Elements in the mesh : (1st Node, 1st Element, last Node, last Element).
    logical, dimension(2)               :: CMD      ! Moving or fixed body.
    real(rp), dimension(3)              :: GBody    ! Gravity center with respect to the inertial frame.
    real(rp), dimension(3)              :: BBody    ! Buoyancy Center.
    real(rp), dimension(6)              :: CSolv    ! Position (linear, angular) where the motion equation is solved (Center of gravity for classical motion, fix point for fix motion).
    real(rp), dimension(6,6)            :: IBody    ! Inertia.
    real(rp), dimension(3)              :: DimBody  ! Dimensions (ex : length/width/draft).
    real(rp), dimension(6)              :: MBody    ! Linear and rotational displacement of the body (MBody = VBody*dt, CSolv = CSolv_old + VBody*dt).
    real(rp), dimension(6)              :: VBody    ! Linear and rotational velocity of the body (time-differentiation of CSolv).
    real(rp), dimension(6)              :: ABody    ! Linear and rotational acceleration of the body (time-differentiation of VBody).
    real(rp), dimension(12)             :: FBody    ! Used in Energy.f90.
    real(rp), dimension(3,2)            :: Pmline   ! Fixation Point of the mooring line.
    real(rp)                            :: Aire     ! Surface.
    real(rp)                            :: Volume   ! Volume.
    real(rp)                            :: Mass     ! Mass.
    real(rp),dimension(:,:),allocatable :: CK       ! CK coefficients at A if iFixPoint, at G otherwise.
    real(rp),dimension(:,:),allocatable :: CT       ! CT coefficients at A if iFixPoint, at G otherwise.
    real(rp),dimension(:),allocatable   :: Q        ! Q coefficients.
    real(rp),dimension(:),allocatable   :: Th       ! Th coefficients.
    logical                             :: Active   ! Active body in the WSC simulation (immerged of piercing body) (True) or not (False).
end type TBody

! Free surface parameters
type TFS
    integer, dimension(4) :: IndFS                  ! Indices de 1er et dernier Noeud/Facette dans le maillage : (1er Noeud, 1ere Facette, dernier Noeud, dernière Facette)
    real(rp), dimension(5) :: DimFS                 ! Dimensions (ex : Radius of the mesh/mean damping radius/depth)
    logical, dimension(2) :: CMD
    real(rp)               :: Aire                 ! Area of the Free Surface (introduced for the filtering)
end type TFS

type TMaillage
	real(rp), dimension(3)                      :: Origine  ! Origine du maillage.
	integer                                     :: Nnoeud   ! Nombre de noeuds total du maillage (corps, surface libre et cuve).
	type(TNoeud) ,dimension(:), allocatable     :: Tnoeud	! Table des noeuds du maillage.
	integer                                     :: Nfacette ! Nombre de facettes du maillage.
	type(TFacette), dimension(:), allocatable   :: Tfacette	! Table des noeuds sommets des facettes.
	integer                                     :: NBody    ! Total number of bodies (tank included).
	type(TBody), dimension(:), allocatable      :: Body     ! Table des corps du maillage.
	type(TFS)                                   :: FS       ! Free-surface parameters.
	integer                                     :: Nsys     ! Nombre d'inconnues dans BVP (différent de Nnoeud si la cuve est fixe par exemple).
	integer                                     :: Nfsys    ! Nombre de facettes liées aux inconnues (juste pour visu).
	integer                                     :: TypeM    ! Type of the mesh : 0 --> circular, 1 --> rectangular.
	real(rp), dimension(5)                      :: DimTank  ! Dimensions of the tank ( Length/diameter, width/diameter, depth, damping dimension parameter, any other parameter).    
end type TMaillage

type TPotentiel
	real(rp) :: incident		! Valeur de la composante incidente du potentiel
	real(rp) :: perturbation    ! Valeur de la composante de perturbation
end type TPotentiel

type TEcoulement
	type(TPotentiel), dimension(:),     allocatable :: Phi          ! Velocity potential.
	type(TPotentiel), dimension(:),     allocatable :: DPhiDn       ! Normal derivative of the velocity potential.
	type(TPotentiel), dimension(:,:),   allocatable :: GPhi         ! Gradient of the velocity potential.
	type(TPotentiel), dimension(:,:),   allocatable :: GPhi2        ! dérivées de Phi dans la base locale (dPhi/ds1, dPhi/ds2, d²Phi/ds1², d²Phi/ds2², d²Phi/ds1dn, d²Phi/ds2dn)
	type(TPotentiel), dimension(:,:,:), allocatable :: GradGrad     ! test
	type(TPotentiel), dimension(:,:),   allocatable :: DGPhiDz      ! 
	type(TPotentiel), dimension(:),     allocatable :: DPhiDt       ! Partial time derivative of the velocity potenial (d_rond/d_rondt) -> use in the 2nd BEM problem.
	type(TPotentiel), dimension(:),     allocatable :: DpPhiDt      ! Particular time derivative of the velocity potential -> use to time step the velocity potential on the free surface.
	type(TPotentiel), dimension(:),     allocatable :: DDPhiDnDt    
	type(TPotentiel), dimension(:),     allocatable :: DGradPhiSqDn   ! Dérivée normale de grad de phi carré
	type(TPotentiel), dimension(:,:),   allocatable :: GsDPhiDn       ! Gradient surfacique de la dérivée normale dans la base locale
	type(TPotentiel), dimension(:),     allocatable :: Eta
	type(TPotentiel), dimension(:,:),   allocatable :: GEta
	type(TPotentiel), dimension(:),     allocatable :: DEtaDt
	type(TPotentiel), dimension(:),     allocatable :: Pression
	type(TPotentiel), dimension(:,:),   allocatable :: Fpression
end type TEcoulement

!	type HouleRF
type THouleRF
    integer                     :: index    ! Index, fixed equal to 1.
    integer                     :: N        ! Nombre de coefficients
    real(rp), dimension(0:50)   :: An, Bn   ! Coefficients de Houle R&F
    real(rp)                    :: Landa    ! Longueur d'onde
    real(rp)                    :: Hcc      ! Hauteur crête à creux
    real(rp)                    :: k        ! Nombre d'onde
    real(rp)                    :: T        ! Période de la houle
    real(rp)                    :: C        ! Célérité
    real(rp)                    :: Cp       ! Vitesse Particulaire
    real(rp)                    :: Cc       ! Vitesse du Courant
    real(rp)                    :: d        ! Profondeur
end type THouleRF

type THoule
    integer :: Htype
    integer :: Nhoule
    real(rp):: depth
    logical :: infinite_depth
    real(rp), allocatable :: w(:), dir(:), Aphi(:), konde(:), lambda(:)
end type THoule

! type Input data
type InputDataStruct
    
    ! This type is used to manage the input data of the bodies between their reading and the final structure Mesh.
    
    ! Fundamental parameters
    logical,dimension(:),allocatable            :: free_body        ! Free motion: 1, fixed body or forced motion: 0.
    logical,dimension(:),allocatable            :: MovingBody       ! FreeBody: 1, forced body: 0.
    logical,dimension(:),allocatable            :: DeformBody       ! Remeshing of the body: true, no: false. Created in InputDataStruct and not TMaillage because when it is initialized, TMaillage does not exit.
    integer,dimension(:),allocatable            :: is_piercing      ! Immerged body = 0, piercing = 1, above the free surface = 2. Created in InputDataStruct and not TMaillage because when it is initialized, TMaillage does not exit.
    
    ! Geometrical properties
    logical,dimension(:),allocatable            :: is_Mesh_file     ! Mesh file has to be read (T) or not (F).
    character(len=50),dimension(:),allocatable  :: Mesh_file        ! Mesh file if Mesh_type = 3.
    integer,dimension(:),allocatable            :: igtype           ! Type of geometry (1 = cube, 2 = cylinder, 3 = Nemoh.dat if mesh strategy of LL, 5 = axisym).
    real(rp), dimension(:,:),allocatable        :: Lgeom            ! Length (1) and radius (2) of the body.
    character(len=50),dimension(:),allocatable  :: file_axisym      ! Axisym file.
    
    ! Mesh properties
    integer,dimension(:),allocatable            :: NphiSphere       ! Parameter in case of using the mesh strategy of LL.
    real(rp),dimension(:),allocatable           :: dx2              ! Length of the panels on the body.
        
    ! Physical properties
    real(rp),dimension(:,:,:),allocatable       :: Position         ! Position (linear and angular (rad)) of the orgine of the body frame from the inertial frame (OeOj).
    real(rp),dimension(:,:),allocatable         :: PositionG        ! Center of gravity in the body frame (OjGj).
    real(rp),dimension(:,:),allocatable         :: GOt0             ! Position of the center of the mesh - Center of gravity at t = 0 s = Fixed position in the body frame of the center of the mesh wrt the center of gravity . Useful to update the frames of the geometries the subroutine Remesh.
    real(rp),dimension(:),allocatable           :: Mass             ! Mass of the body.
    logical,dimension(:),allocatable            :: is_Inertia_file  ! Inertia file has to be read (T) or not (F).
    character(len=50),dimension(:),allocatable  :: file_inertia     ! Inertia file.
    logical,dimension(:),allocatable            :: is_IBodyA        ! Flag to know if the inertia matrix if file_inertia is given in Gj or A.
        
    ! Wave probes
    real(rp),dimension(:,:),allocatable         :: PositionWP       ! Positions of the wave probes.
    
    ! External loads
    real(rp),dimension(:),allocatable           :: Raideur          ! Spring stiffness.
    real(rp),dimension(:),allocatable           :: Lressort         ! Spring length.
    real(rp),dimension(:,:),allocatable         :: Pressort         ! Spring Position.
    real(rp),dimension(:),allocatable           :: MuPTO            ! Damping coefficient.
    real(rp),dimension(:),allocatable           :: Cd_Morison       ! Morison coefficient.
    
    ! Degrees of freedom
    integer,dimension(:),allocatable            :: ndll             ! Number of degrees of freedom.
    logical,dimension(:,:),allocatable          :: dll_dir          ! Degree of freedom vector.
    
    ! Forced motion
    real(rp),dimension(:),allocatable           :: Vcst             ! Constant velocity
    real(rp),dimension(:),allocatable           :: ACorps           ! Motion amplitude.
    real(rp),dimension(:),allocatable           :: WCorps           ! Wave frequency of the motion.
    real(rp),dimension(:),allocatable           :: PhiCorps         ! Phase of the motion.
    real(rp),dimension(:),allocatable           :: LambdACorps      ! Wave length of the motion of the body.
    
end type InputDataStruct

contains
    
subroutine NewMaillage(Mesh, N, nCorps)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(out)        :: Mesh     ! Mesh.
    integer, intent(in)                 :: N        ! Number of nodes.
    integer,intent(in)                  :: nCorps   ! Number of bodies (tank included).
    
    ! This subroutine creates a new structure TMaillage
        
    allocate( Mesh%Tnoeud(N), Mesh%Tfacette(2*N), Mesh%Body(nCorps), STAT=ierr) ! Why 2N for the number of panels? - NBodies+1 = Number of body + the tank.
        
end subroutine NewMaillage

subroutine DelMaillage(Mesh)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh ! Mesh
    
    ! This subroutine deallocates the structure TMaillage
    
    deallocate( Mesh%Tnoeud, Mesh%Tfacette, Mesh%Body)
    
end subroutine DelMaillage

subroutine CopyMaillage(Maillage2, Maillage1)
    
    !f2py integer*1, dimension(1000)    :: Maillage1
    !f2py integer*1, dimension(1000)    :: Maillage2
    type(TMaillage)                     :: Maillage1,Maillage2  ! Meshes.
    
    integer                             :: k                    ! Loop parameter.
    
    ! This subroutine copies a structure TMaillage (1) into a new structure TMaillage (2).
 
    Maillage2%Origine = Maillage1%Origine
    Maillage2%Nnoeud = Maillage1%Nnoeud
    Maillage2%Nfacette = Maillage1%Nfacette
    Maillage2%NBody = Maillage1%NBody
    Maillage2%NSys = Maillage1%NSys
    Maillage2%Nfsys = Maillage1%Nfsys
    Maillage2%TypeM = Maillage1%TypeM
    Maillage2%DimTank = Maillage1%DimTank

    do k=1,Maillage1%Nnoeud
        Maillage2%Tnoeud(k) = Maillage1%Tnoeud(k)
    end do

    do k=1,Maillage1%Nfacette
        Maillage2%Tfacette(k) = Maillage1%Tfacette(k)
    end do

    do k=1,Maillage1%NBody
        Maillage2%Body(k) = Maillage1%Body(k)
    enddo
    Maillage2%FS = Maillage1%FS
    

    return
    
end subroutine CopyMaillage

!	Initialisation d un potentiel ecoulement par un reel
subroutine IniPotentiel( Potentiel, Valeur)
    implicit none
    !f2py integer*1, dimension(1000) :: Potentiel
    type(TPotentiel) :: Potentiel
    real(rp) :: Valeur
    Potentiel%incident=Valeur
    Potentiel%perturbation=Valeur
    return
end subroutine IniPotentiel

!	Copie d une variable de type Potentiel1 dans une variable Potentiel2
subroutine CopyPotentiel( Potentiel2, Potentiel1)
    implicit none
    !f2py integer*2, dimension(2000) :: Potentiel1,Potentiel2
    type(TPotentiel) :: Potentiel1,Potentiel2
    Potentiel2%incident=Potentiel1%incident
    Potentiel2%perturbation=Potentiel1%perturbation
    return
end subroutine CopyPotentiel

!	Construction d une nouvelle variable de type ecoulement
subroutine NewEcoulement( Ecoulement, N)
    implicit none
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement) :: Ecoulement
    integer :: N
    allocate( Ecoulement%GradGrad(3,3,N)) !, Ecoulement%DEtaDn(N))
    allocate( Ecoulement%Phi(N), Ecoulement%DPhiDn(N), Ecoulement%DDPhiDnDt(N), Ecoulement%GPhi(3,N), Ecoulement%GPhi2(6,N), Ecoulement%DGPhiDz(3,N), Ecoulement%DPhiDt(N), Ecoulement%DpPhiDt(N), Ecoulement%DGradPhiSqDn(N), Ecoulement%GsDPhiDn(3,N), STAT=ierr)
    allocate( Ecoulement%Eta(N), Ecoulement%GEta(3,N), Ecoulement%DEtaDt(N), Ecoulement%Pression(N), Ecoulement%Fpression(3,N), STAT=ierr)
return
end subroutine NewEcoulement

subroutine DelEcoulement(Ecoulement)
    
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement)                   :: Ecoulement   ! Flow parameters.
    
    ! This subroutine deletes an object TEcoulement.
    
    deallocate( Ecoulement%GradGrad) !, Ecoulement%DEtaDn)
    deallocate( Ecoulement%Phi, Ecoulement%DPhiDn, Ecoulement%DDPhiDnDt, Ecoulement%GPhi, Ecoulement%GPhi2, Ecoulement%DGPhiDz, Ecoulement%DPhiDt, Ecoulement%DpPhiDt, Ecoulement%DGradPhiSqDn, Ecoulement%GsDPhiDn)
    deallocate( Ecoulement%Eta, Ecoulement%GEta, Ecoulement%DEtaDt, Ecoulement%Pression, Ecoulement%Fpression)
    return
    
end subroutine DelEcoulement

!	Initialisation d un ecoulement par un reel
subroutine IniEcoulement( Ecoulement, N, Valeur)
    implicit none
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement) :: Ecoulement
    real(rp) :: Valeur
    integer :: N
    !f2py integer*1, dimension(1000) ::  PPotentiel	
    type(TPotentiel) :: PPotentiel	
    integer :: j,k,r
    call IniPotentiel( PPotentiel, Valeur)
    do j=1,N
	    call CopyPotentiel( Ecoulement%Phi(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%DPhiDn(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%DDPhiDnDt(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%DGradPhiSqDn(j), PPotentiel)
	    do k=1,3
		    call CopyPotentiel( Ecoulement%GPhi(k,j), PPotentiel)
		    call CopyPotentiel( Ecoulement%DGPhiDz(k,j), PPotentiel)
		    call CopyPotentiel( Ecoulement%GsDPhiDn(k,j), PPotentiel)
		    do r=1,3
		        call CopyPotentiel( Ecoulement%GradGrad(r,k,j), PPotentiel)
		    end do
	    end do
	    do k=1,4
	        call CopyPotentiel( Ecoulement%GPhi2(k,j), PPotentiel)
	    end do
	    call CopyPotentiel( Ecoulement%DPhiDt(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%DpPhiDt(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%Eta(j), PPotentiel)
	    do k=1,3
		    call CopyPotentiel( Ecoulement%GEta(k,j), PPotentiel)
		    call CopyPotentiel( Ecoulement%Fpression(k,j), PPotentiel)
	    end do	
	    call CopyPotentiel( Ecoulement%DEtaDt(j), PPotentiel)
	    call CopyPotentiel( Ecoulement%Pression(j), PPotentiel)
    end do
    return
end subroutine IniEcoulement

subroutine CopyEcoulement( Ecoulement2, Ecoulement1, N)
    
    !f2py integer*2, dimension(2000)    :: Ecoulement1,Ecoulement2
    type(TEcoulement)                   :: Ecoulement1, Ecoulement2 ! Flow parameters.
    integer                             :: N                        ! Number of nodes in the mesh.
    
    integer                             :: j,k,r                    ! Loop parameters.
    
    ! This subroutine copies the object Ecoulement1 into Ecoulement2.
    
    do j = 1,N
	    call CopyPotentiel( Ecoulement2%Phi(j), Ecoulement1%Phi(j))
	    call CopyPotentiel( Ecoulement2%DPhiDn(j), Ecoulement1%DPhiDn(j))
        call CopyPotentiel( Ecoulement2%DDPhiDnDt(j), Ecoulement1%DDPhiDnDt(j))
	    call CopyPotentiel( Ecoulement2%DGradPhiSqDn(j), Ecoulement1%DGradPhiSqDn(j))
	    do k = 1,3
		    call CopyPotentiel( Ecoulement2%GPhi(k,j), Ecoulement1%GPhi(k,j))
		    call CopyPotentiel( Ecoulement2%DGPhiDz(k,j), Ecoulement1%DGPhiDz(k,j))
		    call CopyPotentiel( Ecoulement2%GsDPhiDn(k,j), Ecoulement1%GsDPhiDn(k,j))
		    do r = 1,3
		        call CopyPotentiel( Ecoulement2%GradGrad(r,k,j), Ecoulement1%GradGrad(r,k,j))
		    end do
	    end do
	    do k = 1,4
	        call CopyPotentiel( Ecoulement2%GPhi2(k,j), Ecoulement1%GPhi2(k,j))
	    end do
	    call CopyPotentiel( Ecoulement2%DPhiDt(j), Ecoulement1%DPhiDt(j))
	    call CopyPotentiel( Ecoulement2%DpPhiDt(j), Ecoulement1%DpPhiDt(j))
	    call CopyPotentiel( Ecoulement2%Eta(j), Ecoulement1%Eta(j))
	    do k = 1,3
		    call CopyPotentiel( Ecoulement2%GEta(k,j), Ecoulement1%GEta(k,j))
		    call CopyPotentiel( Ecoulement2%Fpression(k,j), Ecoulement1%Fpression(k,j))
	    end do	
	    call CopyPotentiel( Ecoulement2%DEtaDt(j), Ecoulement1%DEtaDt(j))
	    call CopyPotentiel( Ecoulement2%Pression(j), Ecoulement1%Pression(j))
    end do
    return
    
end subroutine CopyEcoulement

subroutine extract_ecoulement(Ecoulement,t,fileflow,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(inout) :: Ecoulement
    real(rp),intent(inout) :: t
    character(len=50),intent(in) :: fileflow
    integer,intent(inout) :: ierror
    ! local
    character(len=100) :: poub,lignep
    real(rp),dimension(:),allocatable :: tab
    integer :: Nnoeud,Nfacette
    integer :: j
    
    print*,"Extract_ecoulement: this subroutine is outdated."
    call exit()
    
    ierror = 0
    
    open(unit=122,file=fileflow,iostat=ierror)
    if(ierror/=0)then 
        ierror = 100
        goto 9999
    endif
  
    read(122,'(a)') poub
    read(122,'(a)') poub
    read(122,'(a)') lignep  
  
    read(lignep(11:15),*) t
    read(lignep(30:41),*) Nnoeud
    read(lignep(46:58),*) Nfacette
  
    print*,'t = ',t,' ; Nnoeud = ',Nnoeud,' ; Nfacette = ',Nfacette
  
    allocate(tab(19))
  
    do j=1,Nnoeud
        read(122,*) tab(1:19)
        Ecoulement%Phi(j)%perturbation = tab(5)
        Ecoulement%GPhi(1,j)%perturbation = tab(6)
        Ecoulement%GPhi(2,j)%perturbation = tab(7)
        Ecoulement%GPhi(3,j)%perturbation = tab(8)
        Ecoulement%DPhiDt(j)%perturbation = tab(14)
        Ecoulement%Eta(j)%perturbation = tab(4)
        Ecoulement%GEta(1,j)%perturbation = tab(16)
        Ecoulement%GEta(2,j)%perturbation = tab(17)
        Ecoulement%GEta(3,j)%perturbation = tab(18)
        Ecoulement%DEtaDt(j)%perturbation = tab(19)
    enddo
  

    9999 continue
    if(ierror/=0)then
    write(*,99),ierror
    endif
    99 format('** error #',i3,' : extraction ecoulement impossible')

    select case(ierror)
    case(100)
    print*,'** pb. ouverture fichier ecoulement : '!,filename
    end select
  
    close(122)
  
    if(allocated(tab)) deallocate(tab)

end subroutine extract_ecoulement  

subroutine New_InputData(InputData,n,p)
    
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData    ! Input data.
    integer,intent(in)                  :: n            ! Number of bodies.
    integer,intent(in)                  :: p            ! Number of wave probes.

    ! This subroutine initialized the structure InputData.
    
    ! Fundamental parameters
    allocate(InputData%free_body(n))
    allocate(InputData%MovingBody(n))
    allocate(InputData%DeformBody(n))
    allocate(InputData%is_piercing(n))
    
    ! Geometrical properties
    allocate(InputData%is_Mesh_file(n))
    allocate(InputData%Mesh_file(n))
    allocate(InputData%igtype(n))
    allocate(InputData%Lgeom(3,n))
    allocate(InputData%file_axisym(n))
    
    ! Mesh properties
    allocate(InputData%NphiSphere(n))
    allocate(InputData%dx2(n))

    ! Physical properties
    allocate(InputData%Position(3,2,n))
    allocate(InputData%PositionG(3,n))
    allocate(InputData%GOt0(3,n))
    allocate(InputData%Mass(n))
    allocate(InputData%is_Inertia_file(n))
    allocate(InputData%file_inertia(n))
    allocate(InputData%is_IBodyA(n))
        
    ! Wave probes
    allocate(InputData%PositionWP(3,p))
    
    ! External loads
    allocate(InputData%Raideur(n))
    allocate(InputData%Lressort(n))
    allocate(InputData%Pressort(3,n))
    allocate(InputData%MuPTO(n))
    allocate(InputData%Cd_Morison(n))

    ! Degrees of freedom
    allocate(InputData%ndll(n))
    allocate(InputData%dll_dir(6,n))
    
    ! Forced motion
    allocate(InputData%Vcst(n))
    allocate(InputData%ACorps(n))
    allocate(InputData%WCorps(n))
    allocate(InputData%PhiCorps(n))
    allocate(InputData%LambdACorps(n))

end subroutine New_InputData

subroutine Initialization_InputData(InputData)
    
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData    ! Input data.
    
    ! This subroutine initialized the structure InputData.
    
    ! Fundamental parameters
    InputData%free_body = 0
    InputData%MovingBody = 0
    InputData%DeformBody = 0
    InputData%is_piercing = 0
    
    ! Geometrical properties
    InputData%igtype = 0._RP
    InputData%Lgeom = 0._RP
        
    ! Mesh properties
    InputData%NphiSphere = 0._RP
    InputData%dx2 = 0._RP
        
    ! Physical properties
    InputData%Position = 0._RP
    InputData%PositionG = 0._RP
    InputData%GOt0 = 0._RP
    InputData%Mass = 0._RP
    InputData%is_Inertia_file = 0._RP
    InputData%is_IBodyA = 0._RP
        
    ! External loads
    InputData%Raideur = 0._RP
    InputData%Lressort = 0._RP
    InputData%Pressort = 0._RP
    InputData%MuPTO = 0._RP
    InputData%Cd_Morison = 0._RP
    
    ! Degrees of freedom
    InputData%ndll = 0._RP
    InputData%dll_dir = 0._RP
    
    ! Forced motion
    InputData%Vcst = 0._RP
    InputData%ACorps = 0._RP
    InputData%WCorps = 0._RP
    InputData%PhiCorps = 0._RP
    InputData%LambdACorps = 0._RP
    
end subroutine Initialization_InputData

subroutine Delete_InputData(InputData)
    
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData    ! Input data.
    
    ! This subroutine deallocates the structure InputData.
    
    ! Fundamental parameters
    deallocate(InputData%free_body)
    deallocate(InputData%MovingBody)
    deallocate(InputData%DeformBody)
    deallocate(InputData%is_piercing)
    
    ! Geometrical properties
    deallocate(InputData%is_Mesh_file)
    deallocate(InputData%Mesh_file)
    deallocate(InputData%igtype)
    deallocate(InputData%Lgeom)
    deallocate(InputData%file_axisym)
    
    ! Mesh properties
    deallocate(InputData%NphiSphere)
    deallocate(InputData%dx2)
    
    ! Physical properties
    deallocate(InputData%Position)
    deallocate(InputData%PositionG)
    deallocate(InputData%GOt0)
    deallocate(InputData%Mass)
    deallocate(InputData%is_Inertia_file)
    deallocate(InputData%is_IBodyA)
    deallocate(InputData%file_inertia)
        
    ! Wave probes
    deallocate(InputData%PositionWP)
    
    ! External loads
    deallocate(InputData%Raideur)
    deallocate(InputData%Lressort)
    deallocate(InputData%Pressort)
    deallocate(InputData%MuPTO)
    deallocate(InputData%Cd_Morison)
    
    ! Degrees of freedom
    deallocate(InputData%ndll)
    deallocate(InputData%dll_dir)
    
    ! Forced motion
    deallocate(InputData%Vcst)
    deallocate(InputData%ACorps)
    deallocate(InputData%WCorps)
    deallocate(InputData%PhiCorps)
    deallocate(InputData%LambdACorps)
    
end subroutine Delete_InputData

end module Structuresdonnees
