module MeshStruct
use Constantes
use GeomStruct
use iso_c_binding
implicit none

type MVertex
    integer                             :: index        ! Index dans MGrid
    real(rp),dimension(3)               :: coord        ! Coordinates of the point
    integer                             :: nv_arrete    ! Nombre d'arrete auxquelles le sommet est connecte
    integer,dimension(nvmax)            :: va           ! index des arretes auxquelles le sommet est connecté dans MGrid
    integer                             :: nv_face      ! Nombre de face du maillage auxquelles le sommet est connecte 
    integer,dimension(nvmax)            :: vf           ! index des faces du maillage auxquelles le sommet est connecté dans MGrid
    integer                             :: nv_point     ! Nombre de points du maillage auxquelles le sommet est connecte 
    integer,dimension(nvmax)            :: vp           ! index des points du maillage auxquelles le sommet est connecté dans MGrid
    integer                             :: nface        ! Nombre de face de la geometrie auxquelles le sommet est connecte 
    integer,dimension(nvmax)            :: face         ! index des faces de la geometrie auxquelles le sommet est connecté
    integer                             :: nfaceAxisym  ! Number of face before using the subroutine merge_face which merges the face of the two parts of the axisym geometry.
    integer,dimension(:),allocatable    :: faceAxisym   ! Faces before using the subroutine merge_face which merges the face of the two parts of the axisym geometry.
    integer                             :: nedge        ! Nombre des arretes de la geometrie auxquelles le sommet est connecte 
    integer,dimension(nvmax)            :: edge         ! index des arretes de la geometrie auxquelles le sommet est connecté
    integer                             :: bf           ! Boundary face/flag
end type MVertex

interface assignment (=)
    module procedure assign_mvertex_point,assign_mvertex_coord
end interface

type MEdge
    integer                 :: index    ! Index dans MGrid
    integer,dimension(2)    :: iP       ! Index des points dans MGrid
    integer,dimension(2)    :: iT       ! Index des triangles auquel l'arrete est connecte dans MGrid
    real(rp),dimension(3)   :: P1       ! 1er point
    real(rp),dimension(3)   :: P2       ! 2eme point
    real(rp),dimension(3)   :: dir      ! Vecteur tangent
    real(rp)                :: length   ! Longueur
    integer                 :: nface    ! Numero de la face (surface) ou est MEdge.
    integer,dimension(nvmax):: face     ! Surfaces auxquelles l'arete est reliée. (nvmax <= 2)
    integer,dimension(nvmax):: isens    ! Changement de sens des noeuds
    integer                 :: bf       ! Boundary face/flag
end type MEdge

type MElement
    integer :: index              ! Index dans Mgrid
    integer,dimension(3) :: iP    ! Index des points
    integer,dimension(3) :: iA    ! Index des arretes qui composent l'element dans MGrid
    integer,dimension(3) :: iT    ! Index des elements adjacent dans MGrid
    real(rp),dimension(3) :: P1   ! 1er point
    real(rp),dimension(3) :: P2   ! 2eme point
    real(rp),dimension(3) :: P3   ! 3eme point
    real(rp) :: area              ! Aire
    integer :: face               ! Numero de la face dans la geometrie
    integer :: typeFrontiere      ! 1 = frontiere mat, 0 = FS, 2 = external boundaries ?
end type MElement

type MGrid
    type(MVertex) ,dimension(:),allocatable :: point  ! Points
    type(MEdge)   ,dimension(:),allocatable :: arrete ! Edges
    type(MElement),dimension(:),allocatable :: tri    ! Triangles
end type MGrid
!
! **/ Maillon liste chainee
!
type maillon_chaine
  integer :: index
  integer :: id
  type(maillon_chaine),pointer :: prec,suiv
end type maillon_chaine
!
type chaine_point
  type(point) :: val
  type(chaine_point),pointer :: prec,suiv
end type chaine_point

type chaine_point_pt
  type(chaine_point),pointer :: pt
end type chaine_point_pt
!
! **/ Maillon pile de points (non present dans MGrid)
!
type near_point_cell
  type(MVertex) :: val
  real(rp) :: dist
  type(near_point_cell),pointer :: suiv
end type near_point_cell
!
type pile_element
  integer :: val
  real(rp) :: poid
  integer :: bf
  type(pile_element),pointer :: suiv
end type pile_element


interface add_element
  module procedure add_element_point,add_element_arrete,add_element_tri
end interface

interface depile
  module procedure depile_near_point_cell,depile_pile_element
end interface


contains 

! */ Surcharge des operateurs

! **/ Creation d'un type Vertex à partir d'un type point

subroutine assign_mvertex_point(this,a)
    !f2py integer*1, dimension(1000)    :: this
    type(MVertex),intent(out)           :: this
    !f2py integer*1, dimension(1000)    :: a
    type(point),intent(in)              :: a
    
    integer                             :: j,nface
    
    ! This subroutine initializes a Mvertex structure from a point.
    
    this%index = -9
    do j=1,3
        this%coord(j) = a%coord(j)
    enddo
    this%nv_arrete = 0
    this%va = -9
    nface = a%nface
    this%nface = nface
    this%face(1:nface) = a%face(1:nface)
    this%bf = -9
    this%nv_point = 0
    this%vp = -9
    this%nv_face = 0
    this%vf = -9
    this%nedge = a%nedge
    this%edge(1:a%nedge) = a%edge(1:a%nedge)
            
end subroutine assign_mvertex_point

subroutine assign_mvertex_coord(this,a)

    !f2py integer*1, dimension(1000)    :: this
    type(MVertex),intent(out)           :: this ! Vertex
    !f2py integer*1, dimension(1000)    :: a
    real(rp),dimension(3)               :: a    ! Vector
    
    integer                             :: j    ! Loop parameter
    
    ! This subroutine initializes a Mvertex structure from a vector.
    
    this%index = -9
    do j=1,3
        this%coord(j) = a(j)
    enddo
    this%nv_arrete = 0
    this%va = -9
    this%nv_face = 0
    this%vf = -9
    this%bf = -9
    this%nv_point = 0
    this%vp = -9
    this%nedge = 0
    this%edge = -9
    
end subroutine assign_mvertex_coord

! Creation d'une arrete a partir de deux vertex

subroutine edge_from_vertex(P1,P2,edge)
    
    !f2py integer*1, dimension(1000)    :: P1,P2
    type(MVertex),intent(in),target     :: P1,P2
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge),intent(out)             :: edge
    
    ! This subroutine initializes an edge from two vertexes.
    
    type(vector) :: V,V_tmp
    edge%index = -9
    edge%iP = [P1%index,P2%index]
    edge%iT = -9
    edge%P1 = P1%coord
    edge%P2 = P2%coord
    !V = P2%coord-P1%coord
    call less_array_array(P2%coord,P1%coord,V_tmp)
    call assign_vector_coord(V,V_tmp%coord)
  
    edge%dir = V%coord / V%length
    edge%length = V%length
    edge%bf = -9
    edge%nface = 0
    
end subroutine edge_from_vertex

subroutine triangle_from_edge(E1,E2,E3,tri)
  !f2py integer*1, dimension(1000) :: E1,E2,E3
  type(MEdge),intent(in) :: E1,E2,E3
  !f2py integer*1, dimension(1000) :: tri
  type(MElement),intent(out) :: tri
  ! local
  integer :: index1,index2,index3
  index1 = E3%iP(1)
  index2 = E3%iP(2)
  if(E2%iP(1).ne.index1 .and. E2%iP(1).ne.index2)then
    index3 = E2%iP(1)
  elseif(E2%iP(2).ne.index1 .and. E2%iP(2).ne.index2) then
    index3 = E2%iP(2)
  endif
  tri%index = -9
  tri%iP = [index1,index2,index3]
  tri%iA = [E1%index,E2%index,E3%index]
  tri%iT = -9
  tri%P1 = -9
  tri%P2 = -9
  tri%P3 = -9
  tri%iT = [max(E1%iT(1),E1%iT(2)),max(E2%iT(1),E2%iT(2)),max(E3%iT(1),E3%iT(2))]
  tri%area = -9.
end subroutine triangle_from_edge
  

! Inegalite entre arrete
subroutine ineg_edge(E1,E2,res)
  !f2py integer*1, dimension(1000) :: E1,E2
  type(MEdge),intent(in) :: E1,E2
  logical,intent(out) :: res
  ! local
  logical :: bool1,bool2
  bool1 = E1%iP(1).eq.E2%iP(1) .and. E1%iP(2).eq.E2%iP(2)
  bool2 = E1%iP(1).eq.E2%iP(2) .and. E1%iP(2).eq.E2%iP(1)
  res = .not.(bool1 .or. bool2)
end subroutine ineg_edge

subroutine init_mesh(mesh,nb_point,nb_arete,nb_tri)
    
    !f2py integer*1, dimension(1000)    :: mesh                     
    type(MGrid),intent(inout)           :: mesh                     ! Mgrid.
    integer,intent(inout)               :: nb_point,nb_arete,nb_tri ! Number of points, edges and panels.
        
    ! This subroutine initialized a MGrid structure.
    
    if(allocated(mesh%point)) deallocate(mesh%point)
    if(allocated(mesh%arrete)) deallocate(mesh%arrete)
    if(allocated(mesh%tri)) deallocate(mesh%tri)
    
    nb_point = 0
    nb_arete = 0
    nb_tri   = 0
    
    allocate(mesh%point(10*PointMax))
    allocate(mesh%arrete(10*PointMax))
    allocate(mesh%tri(10*FacetteMax))
  
end subroutine init_mesh

subroutine CopyMGrid(Mesh1,Mesh2,nb_point,nb_arete,nb_tri)
    
    !f2py integer*1, dimension(1000)    :: Mesh1
    !f2py integer*1, dimension(1000)    :: Mesh2
    type(MGrid)                         :: Mesh1,Mesh2              ! Meshes.
    integer,intent(inout)               :: nb_point,nb_arete,nb_tri ! Number of points, edges and panels.
    
    integer                             :: k                        ! Loop parameter.
    
    ! This subroutine copies a structure MGrid (1) into a new structure MGrid (2).
    
    do k=1,nb_point
        mesh2%point(k) = mesh1%point(k)
    end do
    do k=1,nb_arete
        mesh2%arrete(k) = mesh1%arrete(k)
    end do
    do k=1,nb_tri
        mesh2%tri(k) = mesh1%tri(k)
    enddo
    
end subroutine CopyMGrid

subroutine DelMGrid(mesh)
    
    !f2py integer*1, dimension(1000)    :: mesh   
    type(MGrid),intent(inout)   :: mesh ! Mgrid
    
    ! This subroutine deallocates the structure MGrid.
    
    if(allocated(mesh%point)) deallocate(mesh%point)
    if(allocated(mesh%arrete)) deallocate(mesh%arrete)
    if(allocated(mesh%tri)) deallocate(mesh%tri)
    
end subroutine DelMgrid

! */ Fonction deallocation chaine
subroutine dealoc_chain(struct)
  !f2py integer*1, dimension(1000) :: struct
  type(maillon_chaine),pointer :: struct
  !f2py integer*1, dimension(1000) :: ptr
  type(maillon_chaine),pointer :: ptr
  do while (associated(struct%suiv))
    ptr => struct%suiv
    nullify(ptr%prec)
    deallocate(struct)
    struct => ptr
  enddo
  deallocate(struct)
  nullify(struct)
end subroutine 

! Fonction ajout element

! Ajout element a la suite d'un element liste chainee

subroutine add_element_suiv(struct,elem)
  !f2py integer*1, dimension(1000) :: struct
  type(maillon_chaine),pointer :: struct
  integer,intent(in) :: elem
  !f2py integer*1, dimension(1000) :: nouv,suiv
  type(maillon_chaine),pointer :: nouv,suiv
  if (.not.associated(struct)) then
    allocate(struct)
    struct%id = elem
    nullify(struct%prec)
    nullify(struct%suiv)
  elseif (.not.associated(struct%suiv)) then
    allocate(nouv)
    nouv%id = elem
    nouv%prec => struct
    struct%suiv => nouv
    nullify(nouv%suiv)
  elseif (associated(struct%suiv)) then
    suiv => struct%suiv
    allocate(nouv)    
    nouv%id=elem
    nouv%suiv => suiv
    nouv%prec => struct
    struct%suiv => nouv
    suiv%prec   => nouv
  endif
end subroutine add_element_suiv

subroutine add_point_suiv(struct,elem)
  !f2py integer*1, dimension(1000) :: struct
  type(chaine_point),pointer :: struct
  !f2py integer*1, dimension(1000) :: elem
  type(point),intent(in) :: elem
! local
  type(chaine_point),pointer :: nouv,suiv
  if(.not.associated(struct))then
    allocate(struct)
    struct%val = elem
    nullify(struct%suiv)
    nullify(struct%prec)
  elseif(.not.associated(struct%suiv))then
    allocate(nouv)
    nouv%val = elem
    nouv%prec => struct
    struct%suiv => nouv
    nullify(nouv%suiv)
  elseif(associated(struct%suiv))then
    suiv => struct%suiv
    allocate(nouv)
    nouv%val = elem
    nouv%suiv => suiv
    nouv%prec => struct
    struct%suiv => nouv
    suiv%prec => nouv
  endif
end subroutine add_point_suiv

! **/ Ajout element en amont dans la liste chainee

subroutine add_element_prec(struct,elem)
  !f2py integer*1, dimension(1000) :: struct
  type(maillon_chaine),pointer :: struct
  integer,intent(in) :: elem
  !f2py integer*1, dimension(1000) :: nouv,prec
  type(maillon_chaine),pointer :: nouv,prec
  if (.not.associated(struct)) then
    allocate(struct)
    struct%id = elem
    nullify(struct%prec)
    nullify(struct%suiv)
  elseif (.not.associated(struct%prec)) then
    allocate(nouv)
    nouv%id = elem
    nouv%suiv => struct
    struct%prec => nouv
    nullify(nouv%prec)
  elseif (associated(struct%prec)) then
    prec => struct%prec
    allocate(nouv)    
    nouv%id=elem
    nouv%prec => prec
    nouv%suiv => struct
    struct%prec => nouv
    prec%suiv   => nouv
  endif
end subroutine add_element_prec

subroutine add_point_prec(struct,elem)
  !f2py integer*1, dimension(1000) :: struct
  type(chaine_point),pointer :: struct
  !f2py integer*1, dimension(1000) :: elem
  type(point),intent(in) :: elem
! local
  !f2py integer*1, dimension(1000) :: nouv,prec
  type(chaine_point),pointer :: nouv,prec
  if (.not.associated(struct)) then
    allocate(struct)
    struct%val = elem
    nullify(struct%prec)
    nullify(struct%suiv)
  elseif (.not.associated(struct%prec)) then
    allocate(nouv)
    nouv%val = elem
    nouv%suiv => struct
    struct%prec => nouv
    nullify(nouv%prec)
  elseif (associated(struct%prec)) then
    prec => struct%prec
    allocate(nouv)    
    nouv%val = elem
    nouv%prec => prec
    nouv%suiv => struct
    struct%prec => nouv
    prec%suiv   => nouv
  endif
end subroutine add_point_prec

subroutine point_in_chaine(struct,elem,is_pres,ds,Pout,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: struct
  type(chaine_point),intent(inout),pointer :: struct
  !f2py integer*1, dimension(1000) :: elem
  type(point),intent(in) :: elem
  logical,intent(inout) :: is_pres
  real(rp),intent(in) :: ds
  !f2py integer*1, dimension(1000) :: Pout
  type(point),intent(inout) :: Pout
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: ptr0,ptr
  type(chaine_point),pointer :: ptr0,ptr
  !f2py integer*1, dimension(1000) :: P
  type(point) :: P
 
  ierror = 0

  is_pres = .false.
  ptr0 => struct
  !is_pres = ptr0%val==elem
  is_pres = norm2(ptr0%val%coord-elem%coord)<ds
  if(is_pres)then
    ptr0%val%nedge = elem%nedge
    ptr0%val%edge = elem%edge
    ptr0%val%bf = elem%bf
  endif
  if(associated(ptr0%suiv) .and. .not.is_pres)then
    ptr => ptr0%suiv
    do while(associated(ptr) .and. .not.is_pres .and. .not.associated(ptr,ptr0)) ! la derniere condition evite de boucler
      P = ptr%val
      !is_pres = all(abs(P%coord - elem%coord)<ds)
      is_pres = norm2(P%coord-elem%coord)<ds
      if(is_pres)then
        ptr%val%nedge = elem%nedge
        ptr%val%edge = elem%edge
      endif
      ptr=> ptr%suiv
    enddo
  endif
    
end subroutine point_in_chaine  
!
! **/ Ajout element de type Vertex dans la pile et ordonne selon distance
!
subroutine add_vertex_pile(struct,elem,dist)
  !f2py integer*1, dimension(1000) :: struct
  type(near_point_cell),pointer :: struct
  !f2py integer*1, dimension(1000) :: elem
  type(MVertex),intent(in) :: elem
  !f2py integer*1, dimension(1000) :: dist
  real(rp),intent(in) :: dist
  ! local
  !f2py integer*1, dimension(1000) :: new,aux,prev
  type(near_point_cell),pointer :: new,aux,prev
  logical :: element_added
  integer :: ierror
!  
  ierror = 0
  allocate(new)
  nullify(new%suiv) ! Probably a problem here with the nullify of the pointer. 
  new%val=elem
  new%dist=dist
  nullify(prev)
!  
  if (.not.associated(struct)) then
    struct => new
  else
    element_added = .false.
    aux => struct
    nullify(prev)
    do while (associated(aux) .and. .not.element_added)
      if (aux%dist < new%dist-Epsilon) then
        prev => aux
        aux  => aux%suiv
      else
        new%suiv => aux
        if(.not.associated(prev))then
          struct => new
        else
          prev%suiv => new
        endif
        element_added = .true.
      endif
    end do
! ajout fin de pile
    if (associated(prev) .and. .not.element_added) then
      prev%suiv => new
    elseif(.not.element_added)then
      ierror = 1
      goto 9999
    endif
  endif
  
9999 continue
  if(ierror /= 0) then
    print*,"** error #",ierror," : impossible d'ajouter l'element à la pile"
    deallocate(new)
  endif  

end subroutine add_vertex_pile

subroutine depile_near_point_cell(struct,ierror)
  !f2py integer*1, dimension(1000) :: struct
  type(near_point_cell),pointer :: struct
  integer :: ierror
  !f2py integer*1, dimension(1000) :: ptr
  type(near_point_cell),pointer :: ptr
  ierror = 0
  if(.not.associated(struct)) then
    ierror = 1
    goto 9999
  else
    ptr => struct%suiv
    deallocate(struct)
    struct => ptr
  endif
9999 continue
  if (ierror.ne.0) then  
    print*,'** error #',ierror,' : depile'
  endif
end subroutine depile_near_point_cell

subroutine ajout_fin_pile(pile,elem,dist)
  implicit none
  !f2py integer*1, dimension(1000) :: pile
  type(near_point_cell),pointer :: pile
  !f2py integer*1, dimension(1000) :: elem
  type(MVertex),intent(in) :: elem
  !f2py integer*1, dimension(1000) :: dist
  real(rp),intent(in) :: dist
! local
  !f2py integer*1, dimension(1000) :: ptr,new,prev
  type(near_point_cell),pointer :: ptr,new,prev
  integer :: ierror
!
  ierror = 0
  allocate(new)
  nullify(new%suiv)
  new%val = elem
  new%dist = dist
!
  nullify(prev)
  ptr => pile
  do while (associated(ptr))
    prev => ptr
    ptr => ptr%suiv
  enddo
  if (.not.associated(prev)) then
    pile => new
  else
    prev%suiv => new
  endif

end subroutine ajout_fin_pile

subroutine add_element_pile(struct,elem,poid,bf)
    
    !f2py integer*1, dimension(1000)    :: struct
    type(pile_element),pointer          :: struct
    integer,intent(in)                  :: elem
    real(rp),intent(in)                 :: poid
    integer,intent(in)                  :: bf
        
    !f2py integer*1, dimension(1000)    :: new,aux,prev
    type(pile_element),pointer          :: new,aux,prev
    logical                             :: element_added
    integer                             :: ierror           ! Error flag.
        
    ierror = 0
    allocate(new)
    nullify(new%suiv)
    new%val = elem
    new%poid = poid
    new%bf = bf
    nullify(prev)
    
    if(not(associated(struct)))then
        struct => new
    else
        element_added = .false.
        aux => struct 
        nullify(prev)
        do while (associated(aux) .and. .not.element_added)
            if(aux%poid < new%poid - Epsilon) then ! .or. &
                !&  (new%bf==0 .and. aux%bf==1) ) then
                prev => aux
                aux => aux%suiv
            else
                new%suiv => aux
                if (.not.associated(prev))then
                    struct => new
                else
                    prev%suiv => new
                end if
                element_added = .true.
            end if
        end do
        if(associated(prev) .and. .not.element_added)then
            prev%suiv => new
        elseif(.not.element_added) then
            ierror = 1
            goto 9999
        end if
    end if
        
    9999 continue
    if(ierror /= 0) then
        print*,"** error #",ierror," : impossible d'ajouter l'element à la pile"
        deallocate(new)
    end if  

end subroutine add_element_pile

subroutine depile_pile_element(struct,ierror)
  !f2py integer*1, dimension(1000) :: struct
  type(pile_element),pointer :: struct
  integer :: ierror
  !f2py integer*1, dimension(1000) :: ptr
  type(pile_element),pointer :: ptr
  ierror = 0
  if(.not.associated(struct)) then
    ierror = 1
    goto 9999
  else
    ptr => struct%suiv
    deallocate(struct)
    struct => ptr
  endif
9999 continue
  if (ierror.ne.0) then  
    print*,'** error #',ierror,' : depile'
  endif
end subroutine depile_pile_element

subroutine remove_element_pile(struct,id,ierror)
  !f2py integer*1, dimension(1000) :: struct
  type(pile_element),pointer :: struct
  integer,intent(in) :: id
  integer :: ierror
  ! local
  !f2py integer*1, dimension(1000) :: aux1,aux2
  type(pile_element),pointer :: aux1,aux2
  logical :: is_found
  
  if(.not.associated(struct)) then
    ierror = 1
    goto 9999
  endif
  
  is_found = .false.
  aux1 => struct
  aux2 => struct%suiv
  
  if(struct%val == id)then
    call depile(struct,ierror)
    is_found = .true.
  endif
  
  do while (.not.is_found .and. associated(aux2))
    if (aux2%val == id) then
      aux1%suiv => aux2%suiv
      deallocate(aux2)
      is_found = .true.
    else
      aux1 => aux2
      aux2 => aux2%suiv
    endif
  enddo
  
9999 continue
  if(ierror.ne.0) then
    print*,'** error #',ierror," : impossible d'enlever l'element à la pile"
  endif  

end subroutine remove_element_pile  

subroutine libere_pile(struct,ierror)
  !f2py integer*1, dimension(1000) :: struct
  type(pile_element),pointer :: struct
  integer :: ierror
  !f2py integer*1, dimension(1000) :: ptr
  type(pile_element),pointer :: ptr
  ierror = 0
  do while(associated(struct))
    ptr => struct%suiv
    deallocate(struct)
    struct => ptr
  enddo
end subroutine libere_pile

subroutine add_element_point(struct,nb_point,Vertex)
    
    !f2py integer*1, dimension(1000)    :: struct
    type(MVertex),dimension(*)          :: struct   ! MGrid
    integer,intent(inout)               :: nb_point ! Number of the point
    !f2py integer*1, dimension(1000)    :: Vertex    
    type(MVertex),intent(in)            :: Vertex    ! MVertex
    
    ! This subroutine adds the point into the structure MGrid.
    
    nb_point = nb_point+1
    struct(nb_point) = Vertex
    struct(nb_point)%index = nb_point
    
end subroutine add_element_point

subroutine add_element_arrete(mesh,nb_arrete,arrete)
    
    !f2py integer*1, dimension(1000):: mesh
    type(MGrid)                     :: mesh                 ! MGrid
    integer,intent(inout)           :: nb_arrete            ! Number of the edge
    !f2py integer*1, dimension(1000):: arrete
    type(MEdge),intent(in)          :: arrete               ! Edge
    
    integer                         :: i1,i2                !
    integer                         :: nv_arrete,nv_point   ! Number of edge and point
    
    ! This subroutine adds the edge into the structure MGrid.
    
    ! Ajout element liste arrete dynamique
    nb_arrete = nb_arrete+1
    mesh%arrete(nb_arrete) = arrete
    mesh%arrete(nb_arrete)%index = nb_arrete
    
    ! Mise à jour données sur les noeuds
    i1 = arrete%iP(1)
    i2 = arrete%iP(2)
    nv_arrete = mesh%point(i1)%nv_arrete
    mesh%point(i1)%nv_arrete = nv_arrete+1
    mesh%point(i1)%va(nv_arrete+1) = nb_arrete
    nv_arrete = mesh%point(i2)%nv_arrete
    mesh%point(i2)%nv_arrete = nv_arrete+1
    mesh%point(i2)%va(nv_arrete+1) = nb_arrete 
    
    nv_point = mesh%point(i1)%nv_point
    mesh%point(i1)%nv_point = nv_point+1
    mesh%point(i1)%vp(nv_point+1) = i2
    nv_point = mesh%point(i2)%nv_point
    mesh%point(i2)%nv_point = nv_point+1
    mesh%point(i2)%vp(nv_point+1) = i1 
  
end subroutine add_element_arrete

subroutine add_element_tri(struct,nb_tri,tri)
  
    !f2py integer*1, dimension(1000)    :: struct
    type(MElement),dimension(*)         :: struct   ! MGrid
    integer,intent(inout)               :: nb_tri   ! Number of panel
    !f2py integer*1, dimension(1000)    :: tri
    type(MElement),intent(in)           :: tri      ! Panel
    
    ! This subroutine adds the panel/triangle into the structure MGrid.
    
    nb_tri = nb_tri+1
    struct(nb_tri) = tri
    struct(nb_tri)%index = nb_tri
    
end subroutine add_element_tri

end module MeshStruct
