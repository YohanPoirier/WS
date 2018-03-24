module GeomStruct
use Constantes
use Parameters
use FonctionsCommunes
use iso_c_binding

implicit none
 
! ----------- Geometry Structure ------------------

! Repere dans l'espace 3D
type repere3d
    integer :: index
    real(rp),dimension(3) :: origine
    real(rp),dimension(3) :: e1
    real(rp),dimension(3) :: e2
    real(rp),dimension(3) :: e3
    real(rp) :: phi              ! rotation selon e3
    real(rp) :: gamma            ! rotation selon e2
    real(rp) :: psi              ! rotation selon e1
end type repere3d

! Point dans l'espace 3D (remarque bf et flag à modifier par la suite)
type point
  character(len=8)  :: name
  real(rp),dimension(3) :: coord
  integer :: bf                     ! Type de frontiere
  integer :: nface                  ! Nombre de faces auxquelles appartient le point
  integer,dimension(nfvmax) :: face ! Liste des faces
  integer :: nedge
  integer,dimension(nfvmax) :: edge
  integer :: flag                   ! Identification d'une liste d'appartenance du point
end type point

! Operateurs 
interface assignment (=)
  module procedure assign_point
end interface

! Vecteur
type vector
  character(len=8)  :: name
  real(rp) :: coord(3)
  real(rp) :: length
end type vector

! Arete
type GArete
  integer :: index
  type(point) :: P1
  type(point) :: P2
  integer :: nface
  integer,dimension(2) :: face
end type GArete

! Cercle
type GCercle
  integer :: index
  type(repere3d) :: repere
  real(rp) :: rayon,theta
  type(point) :: P1,P2
  integer,dimension(2) :: face
end type GCercle 

! GPolyline
type GPolyline
  integer                               :: index    ! Index
  type(repere3d)                        :: repere   ! Frame
  integer                               :: n        ! Number of points
  real(rp),dimension(:,:),allocatable   :: P        ! Table of points
  integer                               :: nface    ! Number of faces
  integer,dimension(2)                  :: face     ! Faces
end type GPolyline

! Disque (déf. 1)
type cylindrical_surface2
  integer :: index
  type(point) :: centre
  type(vector) :: normal
  real(rp) :: radius
  real(rp) :: gamma
  real(rp) :: phi
end type cylindrical_surface2

! Disque (déf. 2) 
type disque2
  integer :: index
  type(repere3d) :: repere
  real(rp) :: r2max
end type disque2

! Plan
type Gplan
  integer :: index
  type(repere3d) :: repere
  real(rp) :: umin
  real(rp) :: umax
  real(rp) :: vmin
  real(rp) :: vmax
end type Gplan

! Surface cylindrique
type cylindre
  integer :: index
  type(cylindrical_surface2) :: base
  real(rp) :: length
  real(rp) :: phi
  real(rp) :: gamma
  real(rp) :: vmin,vmax
  type(point) :: G
end type cylindre

type cylindre_2
  integer :: index
  type(repere3d) :: repere
  real(rp) :: long
  real(rp) :: rayon
  real(rp) :: vmin
  real(rp) :: vmax
end type cylindre_2

! Surface spherique
type sphere
    integer :: index
    type(repere3d) :: repere
    real(rp) :: radius 
    real(rp) :: umin,umax
    real(rp) :: vmin,vmax
end type sphere

! Surface Conique
type cone
  integer :: index
  type(repere3d) :: repere
  real(rp) :: long
  real(rp) :: r1, r2
  real(rp) :: vmin, vmax
end type cone

! Surface axiSym
type axisym
    integer :: index
    type(repere3d) :: repere
    real(rp) :: umin,umax
    real(rp) :: vmin,vmax
    integer  :: npoint
    real(rp),dimension(:,:),allocatable :: P
end type axisym

! /* Surface Wigley
type Twigley
    integer :: index
    type(repere3d) :: repere
    real(rp) :: L, B, D
    real(rp) :: umin,umax
    real(rp) :: vmin,vmax
end type
!
! */ Objet définit par la liste de ses points, arêtes et faces
type type_geom
  integer :: index
  type(repere3d) :: repere
  integer :: npoint
  type(point),dimension(:),allocatable :: point
  integer :: narete
  type(GArete),dimension(:),allocatable :: arete
  integer,dimension(:),allocatable :: cmd_ar
  integer :: ncercle
  type(GCercle),dimension(:),allocatable :: cercle
  integer,dimension(:),allocatable :: cmd_cer
  integer :: ncylindre
  type(cylindre_2),dimension(:),allocatable :: cylindre
  integer,dimension(:),allocatable :: cmd_cyl
  integer :: ndisque
  type(disque2),dimension(:),allocatable :: disque
  integer,dimension(:),allocatable :: cmd_dis
  integer :: nplan
  type(Gplan),dimension(:),allocatable :: plan
  integer,dimension(:),allocatable :: cmd_pl
  integer :: ncone
  type(cone),dimension(:),allocatable :: cone
  integer, dimension(:),allocatable :: cmd_cone
  integer :: nsphere
  type(sphere),dimension(:),allocatable :: sphere
  integer,dimension(:),allocatable :: cmd_sph
  integer :: naxisym
  type(axisym),dimension(:),allocatable :: axisym
  integer,dimension(:),allocatable :: cmd_axi
  integer :: nwigley
  type(Twigley),dimension(:),allocatable :: wigley
  integer,dimension(:),allocatable :: cmd_wigley
  integer :: npolyline
  type(GPolyline),dimension(:),allocatable :: polyline
  integer :: nconnect
  integer,dimension(:,:),allocatable :: connect
end type type_geom

type type_GeomVect
    
    ! This type is a vector of geometry. There are NBodies gemetries.
    
    type(type_geom),dimension(:),allocatable    :: geom         ! Gemetry
    integer,dimension(:),allocatable            :: nface_vect   ! Number of faces for each body
    logical,dimension(:),allocatable            :: Active       ! Active geometry in the WSC simulation (immerged of piercing body) (True) or not (False).
    
end type type_GeomVect

! Type auxiliaires
type domaine
  integer                               :: dim      ! Number of points
  type(point),dimension(:),allocatable  :: liste    ! Points
end type domaine

! **************************************************************************
! Definition operateurs
! **************************************************************************
! Type point
! --------------------------------------------------------------------------

contains

subroutine assign_point(this,a)
    
    !f2py integer*1, dimension(1000)    :: this
    type(point),intent(inout)           :: this ! Point
    real(rp),dimension(3),intent(in)    :: a    ! Vector
    
    integer                             :: j    ! Loop parameter
    
    ! This subroutine initializes a point from a vector.
    
    this%name = ''
    do j=1,3
        this%coord(j)=a(j)
    enddo    
    this%bf = -9
    this%nedge = 0
    this%edge = -9
    
end subroutine assign_point

subroutine add_point(a,b,c)
  !f2py integer*2, dimension(1000) :: a,b
  type(point),intent(in) :: a,b
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)+b%coord(j)
  enddo
end subroutine add_point

subroutine less_point(a,b,c)
  !f2py integer*2, dimension(1000) :: a,b
  type(point),intent(in)  :: a,b
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)-b%coord(j)
  enddo
end subroutine less_point

subroutine add_point_vector(a,b,c)
  !f2py integer*1, dimension(1000) :: a
  type(point),intent(in) :: a
  !f2py integer*1, dimension(1000) :: b
  type(vector),intent(in) :: b
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)+b%coord(j)
  enddo
end subroutine add_point_vector

subroutine add_vector_point(a,b,c)
  !f2py integer*1, dimension(1000) :: b
  type(point),intent(in) :: b
  !f2py integer*1, dimension(1000) :: a
  type(vector),intent(in) :: a
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)+b%coord(j)
  enddo
end subroutine add_vector_point 

subroutine add_array_array(a,b,c)
  !f2py integer*1, dimension(1000) :: a
  real(rp),dimension(3),intent(in) :: a
  !f2py integer*1, dimension(1000) :: b
  real(rp),dimension(3),intent(in) :: b
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a(j)+b(j)
  enddo
end subroutine add_array_array

subroutine less_point_vector(a,b,c)
  !f2py integer*1, dimension(1000) :: a
  type(point),intent(in) :: a
  !f2py integer*1, dimension(1000) :: b
  type(vector),intent(in) :: b
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)-b%coord(j)
  enddo
end subroutine less_point_vector

subroutine less_vector_point(a,b,c)
  !f2py integer*1, dimension(1000) :: b
  type(point),intent(in) :: b
  !f2py integer*1, dimension(1000) :: a
  type(vector),intent(in) :: a
  !f2py integer*1, dimension(1000) :: c
  type(point),intent(out) :: c
  integer :: j
  do j=1,3
    c%coord(j)=a%coord(j)-b%coord(j)
  enddo
end subroutine less_vector_point

subroutine equal_point(a,b,bool)
  !f2py integer*1, dimension(1000) :: a
  type(point),intent(in) :: a
  !f2py integer*1, dimension(1000) :: b
  type(point),intent(in) :: b
  logical,intent(out) :: bool
  bool = abs(a%coord(1)-b%coord(1)).le.Epsgeom
  bool = bool .and. abs(a%coord(2)-b%coord(2)).le.Epsgeom
  bool = bool .and. abs(a%coord(3)-b%coord(3)).le.Epsgeom
end subroutine equal_point
    
! --------------------------------------------------------------------------
! Type vector
! --------------------------------------------------------------------------
subroutine assign_vector_coord(this,a)
    
    !f2py integer*1, dimension(1000)    :: this
    type(vector),intent(out)            :: this     ! Vector
    real(rp),dimension(3),intent(in)    :: a        ! Table of size 3
    real(rp)                            :: x,y,z    ! Coordinates
  
    ! This subroutine initializes a structure vector from a vector.
    
    ! Definition
    x = a(1)
    y = a(2)
    z = a(3)
    this%name = ''
    this%coord(1) = x
    this%coord(2) = y
    this%coord(3) = z
    
    ! Length
    this%length = sqrt(x*x+y*y+z*z)
    if(abs(this%length).lt.Epsilon)then
        print*,"assign_vector_coord: the length of the vector is nul."
    end if
    
end subroutine assign_vector_coord

subroutine dot_vector_scalar(V,a,Vout)
  !f2py integer*1, dimension(1000) :: V
  type(vector),intent(in)  :: V
  real(rp),intent(in)  :: a
  !f2py integer*1, dimension(1000) :: Vout
  type(vector),intent(out) :: Vout
  integer :: j
  do j=1,3
    Vout%coord(j)=V%coord(j)*a
  enddo
  ! Length
    Vout%length = sqrt(Vout%coord(1)*Vout%coord(1)+Vout%coord(2)*Vout%coord(2)+Vout%coord(3)*Vout%coord(3))
    if(abs(Vout%length).lt.Epsilon)then
        print*,"dot_vector_scalar: the length of the vector is nul."
    end if
    
end subroutine dot_vector_scalar

subroutine dot_scalar_vector(a,V,Vout)
  !f2py integer*1, dimension(1000) :: V
  type(vector),intent(in)  :: V
  !f2py integer*1, dimension(1000) :: a
  real(rp),intent(in)  :: a
  !f2py integer*1, dimension(1000) :: Vout
  type(vector),intent(out) :: Vout
  integer :: j
  do j=1,3
    Vout%coord(j)=V%coord(j)*a
  enddo
  ! Length
    Vout%length = sqrt(Vout%coord(1)*Vout%coord(1)+Vout%coord(2)*Vout%coord(2)+Vout%coord(3)*Vout%coord(3))
    if(abs(Vout%length).lt.Epsilon)then
        print*,"dot_scalar_vector: the length of the vector is nul."
    end if
end subroutine dot_scalar_vector

subroutine dot_vector_vector(V1,V2,Vout)
  !f2py integer*2, dimension(1000) :: V1,V2
  type(vector),intent(in)  :: V1,V2
  !f2py integer*1, dimension(1000) :: vout
  type(vector),intent(out) :: Vout
  integer :: j
  do j=1,3
    Vout%coord(j)=V1%coord(j)*V2%coord(j)
  enddo
  ! Length
    Vout%length = sqrt(Vout%coord(1)*Vout%coord(1)+Vout%coord(2)*Vout%coord(2)+Vout%coord(3)*Vout%coord(3))
    if(abs(Vout%length).lt.Epsilon)then
        print*,"dot_vector_vector: the length of the vector is nul."
    end if
end subroutine dot_vector_vector

subroutine add_vector(V1,V2,Vout)
  !f2py integer*2, dimension(1000) :: V1,V2
  type(vector),intent(in)  :: V1,V2
  !f2py integer*1, dimension(1000) :: Vout
  type(vector) :: Vout
  integer :: j
  do j=1,3
    Vout%coord(j)=V1%coord(j)+V2%coord(j)
  enddo
  ! Length
    Vout%length = sqrt(Vout%coord(1)*Vout%coord(1)+Vout%coord(2)*Vout%coord(2)+Vout%coord(3)*Vout%coord(3))
    if(abs(Vout%length).lt.Epsilon)then
        print*,"add_vector: the length of the vector is nul."
    end if
end subroutine add_vector

subroutine less_vector(V1,V2,Vout)
  !f2py integer*2, dimension(1000) :: V1,V2
  type(vector),intent(in)  :: V1,V2
  !f2py integer*1, dimension(1000) :: vout
  type(vector) :: Vout
  integer :: j
  do j=1,3
    Vout%coord(j)=V1%coord(j)-V2%coord(j)
  enddo
  ! Length
    Vout%length = sqrt(Vout%coord(1)*Vout%coord(1)+Vout%coord(2)*Vout%coord(2)+Vout%coord(3)*Vout%coord(3))
    if(abs(Vout%length).lt.Epsilon)then
        print*,"less_vector: the length of the vector is nul."
    end if
end subroutine less_vector  

subroutine less_array_array(a,b,V)
  real(rp),dimension(3),intent(in) :: a
  real(rp),dimension(3),intent(in) :: b
  !f2py integer*1, dimension(1000) :: V
  type(vector) :: V
  integer :: j
  do j=1,3
    V%coord(j)=a(j)-b(j)
  enddo
  ! Length
    V%length = sqrt(V%coord(1)*V%coord(1)+V%coord(2)*V%coord(2)+V%coord(3)*V%coord(3))
    if(abs(V%length).lt.Epsilon)then
        print*,"less_array_array: the length of the vector is nul."
    end if
end subroutine less_array_array

subroutine Computation_vect_product_vect(A,B,vect_product)
    
    real(rp), dimension(3), intent(in)  :: A,B          ! Two vectors.
    !f2py integer*1, dimension(1000)    :: vect_product
    type(vector),intent(out)            :: vect_product ! Cross product of the two vectors.
    
    ! This subroutine computes the cross product.
    
    vect_product%coord(1)=A(2)*B(3)-A(3)*B(2)
    vect_product%coord(2)=A(3)*B(1)-A(1)*B(3)
    vect_product%coord(3)=A(1)*B(2)-A(2)*B(1)
    
    ! Length
    vect_product%length = sqrt(vect_product%coord(1)*vect_product%coord(1)+vect_product%coord(2)*vect_product%coord(2)+vect_product%coord(3)*vect_product%coord(3))
    if(abs(vect_product%length).lt.Epsilon)then
        print*,"Computation_vect_product_vect: the length of the vector is nul."
    end if
    
end subroutine Computation_vect_product_vect

subroutine assign_repere(index,O,e1,e2,e3,phi,gamma,geom)
    
    integer,intent(in) :: index
    real(rp),dimension(3),intent(in) :: O
    real(rp),dimension(3),intent(in) :: e1,e2,e3
    real(rp),intent(in) :: phi,gamma
    !f2py integer*1, dimension(1000) :: geom
    type(repere3d),intent(out) :: geom
    
    ! This subroutine initialized a frame.
    
    geom%index = index
    geom%origine = O
    geom%e1 = e1
    geom%e2 = e2
    geom%e3 = e3
    geom%phi = phi
    geom%gamma = gamma
    
end subroutine assign_repere

subroutine plan_from_point(P1,P2,P3,P4,plan)
  implicit none
  !f2py integer*4, dimension(1000) :: P1,P2,P3,P4
  type(point),intent(in) :: P1,P2,P3,P4
  !f2py integer*1, dimension(1000) :: plan
  type(Gplan) :: plan
  !f2py integer*3, dimension(1000) :: e1,e2,e3
  type(vector) :: e1,e2,e3
  !f2py integer*1, dimension(1000) :: O
  type(point)  :: O
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep

  !e1 = P2%coord - P1%coord
  call less_array_array(P2%coord,P1%coord,e1)
  
  !e2 = P4%coord - P1%coord
  call less_array_array(P4%coord,P1%coord,e2)
  
  call Computation_vect_product_vect(e1%coord,e2%coord,e3)
  call add_array_array(0.5*P1%coord,0.5*P3%coord,O)
  
  rep%origine = O%coord
  rep%e1 = e1%coord/e1%length
  rep%e2 = e2%coord/e2%length
  rep%e3 = e3%coord/e3%length
  
  plan%repere = rep
  plan%umax = e1%length/2.
  plan%umin = -plan%umax
  plan%vmax = e2%length/2.
  plan%vmin = -plan%vmax  
  
end subroutine plan_from_point  

subroutine arete_from_point(P1,P2,arete)
  implicit none
  !f2py integer*2, dimension(1000) :: P1,P2
  type(point),intent(in) :: P1,P2
  !f2py integer*1, dimension(1000) :: arete
  type(GArete),intent(out) :: arete
  arete%P1 = P1
  arete%P2 = P2
  arete%nface = 0
  arete%face = -9
end subroutine arete_from_point  

subroutine GCercle_from_point(P1,P2,P3,cercle,ierror)
  implicit none
  real(rp),dimension(3),intent(in) :: P1,P2,P3
  !f2py integer*1, dimension(1000) :: cercle
  type(GCercle),intent(inout) :: cercle
  integer,intent(inout) :: ierror
  !f2py integer*1, dimension(1000) :: repere
  type(repere3d) :: repere
  real(rp),dimension(3) :: V1,V2,M1,M2
  real(rp),dimension(3) :: C
  real(rp),dimension(3) :: e1,e2,e3
  real(rp) :: r,norm2
  real(rp),dimension(3,3) :: mat
  real(rp),dimension(3) :: vect
  integer,dimension(3) :: ipiv

  V1 = P3-P1
  M1 = 0.5*(P3+P1)
  V2 = P2-P1
  M2 = 0.5*(P2+P1)
  
  call Computation_vect_product(V1,V2,e3)
  norm2 = sqrt(e3(1)*e3(1)+e3(2)*e3(2)+e3(3)*e3(3))
  e3 = e3/norm2
  
  mat(1,:) = V1
  mat(2,:) = V2
  mat(3,:) = e3  
  vect = [dot_product(M1,V1),dot_product(M2,V2),dot_product(e3,P1)]
  
  call dgesv(3,1,mat,3,ipiv,vect,3,ierror)
  
  C = vect
  
  e1 = P1-C
  r = sqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
  e1 = e1/r

  call Computation_vect_product(e3,e1,e2)
!  
  repere%origine = C
  repere%e1 = e1
  repere%e2 = e2
  repere%e3 = e3
!  
  cercle%index = -9
  cercle%repere = repere
  cercle%rayon = r
  
9999 continue
  if(ierror/=0)then
    write(*,99),ierror
  endif
99 format('** error #',i3,' : creation du cercle impossible')  
   
end subroutine GCercle_from_point  
 
subroutine GCercle_from_2PointCenter(P1,P2,C,cercle,ierror)
    
    real(rp),dimension(3),intent(in)    :: P1,P2,C
    !f2py integer*1, dimension(1000)    :: cercle
    type(GCercle),intent(inout)         :: cercle
    integer,intent(inout)               :: ierror
    
    real(rp)                            :: r1,r3
    !f2py integer*1, dimension(1000)    :: repere
    type(repere3d)                      :: repere
    real(rp),dimension(3)               :: V1,V2,V3    

    ierror = 0

    V1 = P1-C
    V2 = P2-C
    call Computation_vect_product(V1,V2,V3)

    repere%origine = C

    r1 = norm2(V1)
    r3 = norm2(V3)
    
    if(abs(r1).lt.Epsilon .or. abs(r3).lt.Epsilon)then
        ierror = 100
        goto 9999
    else
        repere%e1 = V1/r1
        repere%e3 = V3/r3
        call Computation_vect_product(repere%e3,repere%e1,repere%e2)
    endif

    cercle%index = -9
    cercle%repere = repere
    cercle%rayon = r1
    if(symmetry)then
        cercle%theta = PI
    else
        cercle%theta = 2*PI
    end if

9999 continue
    if(ierror/=0)then
      write(*,99),ierror
    endif
99 format('** error #',i3,' : not able to create the cercle from given point')

end subroutine GCercle_from_2PointCenter

subroutine init_geom(geom,rep,id,ndim,ntab,ierror)
    
    !f2py integer*1, dimension(1000)    :: geom
    type(type_geom),intent(inout)       :: geom     ! Geometry.
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep      ! Frame.
    integer,dimension(*),intent(in)     :: id       ! Kind of structures.
    integer,dimension(*),intent(in)     :: ndim     ! Number of each structure.
    integer,intent(in)                  :: ntab     ! Number of kind of different structures.
    integer,intent(inout)               :: ierror   ! Error flag.
    
    integer :: j,idim
    
    ! This subroutine allocates and initializes the geometry.
    
    ierror = 0

    ! Nullify all dimension.

    geom%npoint = 0
    geom%narete = 0
    geom%ncercle = 0
    geom%ncylindre = 0
    geom%ndisque = 0
    geom%nplan = 0
    geom%ncone = 0
    geom%nsphere = 0
    geom%nwigley = 0
    geom%nconnect = 0

    ! Specify repere.
    
    geom%repere = rep

    ! Allocate array.

    do j = 1,ntab
      idim = ndim(j)
      if(idim.gt.0)then
        select case (id(j))
        case(1)
            allocate(geom%point(idim))
        case(2)
            allocate(geom%arete(idim))
        case(3)
            allocate(geom%cercle(idim))
        case(4)
            allocate(geom%cylindre(idim))
            allocate(geom%cmd_cyl(idim))
            geom%cmd_cyl(:) = 0
        case(5)
            allocate(geom%disque(idim))
            allocate(geom%cmd_dis(idim))
            geom%cmd_dis(:) = 0
        case(6)
            allocate(geom%plan(idim))
            allocate(geom%cmd_pl(idim))
            geom%cmd_pl(:) = 0
        case(7)
            allocate(geom%cone(idim))
            allocate(geom%cmd_cone(idim))
            geom%cmd_cone(:) = 0
        case(8)
            allocate(geom%sphere(idim))
            allocate(geom%cmd_sph(idim))
            geom%cmd_sph(:) = 0
        case(9)
            allocate(geom%connect(idim,2))
        case(10)
            allocate(geom%axisym(idim))
            allocate(geom%cmd_axi(idim))
            geom%cmd_axi(:) = 0
        case(11)
            allocate(geom%polyline(idim))
        case(12)
            allocate(geom%wigley(idim))
            allocate(geom%cmd_wigley(idim))
        case default
            ierror = 100
            goto 9999
        end select
      else
        ierror = 200
        goto 9999
      endif
    enddo

9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
99 format('** error #',i3,' : wrong definition of geometry.')

end subroutine init_geom

subroutine init_GeomVect(VectGeom)
    
    implicit none
    !f2py integer*1, dimension(1000)    :: Vectgeom
    type(type_GeomVect),intent(inout)   :: Vectgeom ! Geometries.
        
    ! This subroutine allocates and initializes the structure type_VectGeom.
    
    allocate(VectGeom%geom(NBodies))
    allocate(VectGeom%nface_vect(NBodies+1)) ! + 1 for the tank.
    allocate(VectGeom%Active(NBodies))
    VectGeom%nface_vect = 0._RP
    
end subroutine init_GeomVect

subroutine deallocate_geom(geom,ierror)
    implicit none
    !f2py integer*1, dimension(1000)    :: geom
    type(type_geom),intent(inout)       :: geom     ! Geometry
    integer,intent(inout)               :: ierror   ! Error flag
    
    ! This subroutine deallocates the structure type_geom.
    
    ierror = 0

    if(allocated(geom%point))     deallocate(geom%point)
    if(allocated(geom%arete))     deallocate(geom%arete)
    if(allocated(geom%cercle))    deallocate(geom%cercle)
    if(allocated(geom%cylindre))  deallocate(geom%cylindre)
    if(allocated(geom%disque))    deallocate(geom%disque)
    if(allocated(geom%plan))      deallocate(geom%plan)
    if(allocated(geom%cone))      deallocate(geom%cone)
    if(allocated(geom%sphere))    deallocate(geom%sphere)
    if(allocated(geom%wigley))    deallocate(geom%wigley)

    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    end if
    99  format('** error #',i3,' : deallocation geom')    
    
end subroutine deallocate_geom  

subroutine deallocate_GeomVect(VectGeom,ierror)
    
    implicit none
    !f2py integer*1, dimension(1000)    :: Vectgeom
    type(type_GeomVect),intent(inout)   :: Vectgeom     ! Geometries
    integer,intent(inout)               :: ierror       ! Error flag
    
    integer                             :: j            ! Loop parameter
    
    ! This subroutine deallocates the structure type_VectGeom.
    
    do j = 1,NBodies
        call deallocate_geom(VectGeom%geom(j),ierror)
    end do
    deallocate(VectGeom%nface_vect)
    
end subroutine deallocate_GeomVect

end module GeomStruct
