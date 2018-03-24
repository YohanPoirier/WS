module MAdvanceFront
use Constantes
use GeomStruct
use MeshFonct
implicit none

contains

subroutine advancefront(mesh,nb_point,nb_arrete,nb_tri,xgrid,ygrid,matdref,iface,&
&                       imetric,param,size_param,nx,ny,ierr,InputData,tin,ioMeshevol,fgeom)

    !f2py integer*1, dimension(1000)            :: mesh
    type(MGrid),intent(inout)                   :: mesh                     ! MGrid.
    integer,intent(inout)                       :: nb_point,nb_arrete,nb_tri! Number of points, edges and triangles for the present surface.
    integer,intent(inout)                       :: ierr                     ! Error flag.
    integer,intent(in)                          :: nx,ny                    ! Size of xgrid,ygrid ad matdref
    real(rp),dimension(nx),intent(in)           :: xgrid                    ! Cartesian grid.
    real(rp),dimension(ny),intent(in)           :: ygrid                    ! Cartesian grid.
    real(rp),dimension(nx,ny,3),intent(in)      :: matdref                  ! Reference length for each cell of the Cartesian grid.
    integer,intent(in)                          :: iface                    ! Number of the face.
    integer,intent(in)                          :: imetric                  ! Metric (0 for a cylinder, a disc and a plan, 1 for a sphere, 2 for a cone, 3 for an axisym and 4 for a Wigley hull).
    integer,intent(in)                          :: size_param               ! Size of param.
    real(rp),dimension(size_param),intent(in)   :: param                    ! 8 param for a cone, 2*axisym%npoint+1 for a axisym, 1 for a cylinder, a disc, a sphere, a plan and a Wigley hull.
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(in)            :: InputData                ! Input data.
    real(rp),intent(inout)                      :: tin                      ! Plotting time for Advance_front.dat.
    integer,intent(in)                          :: ioMeshevol               ! Plotting or not Advance_front.dat or Advance_front_Remesh.dat.
    !f2py integer*1, dimension(1000)            :: fgeom
    type(type_geom),intent(in)                  :: fgeom                    ! Geometry.
            
    !f2py integer*1, dimension(1000)            :: front_edge,front_vertex
    type(pile_element),pointer                  :: front_edge,front_vertex  ! Edges of the mesh front (front_edge) and nodex of the mesh front (front_vertex).
    !f2py integer*1, dimension(1000)            :: ptr
    type(pile_element),pointer                  :: ptr                      ! Pointer.
    !f2py integer*1, dimension(1000)            :: new_point
    type(MVertex)                               :: new_point                ! New point.
    integer                                     :: j                        ! Loop parameter.
            
    ! This subroutine uses the advance front method.
    
    ! Creation frontiere fermee.
    nullify(front_edge)
    nullify(front_vertex)
    
    do j = 1,nb_arrete
        if(mesh%arrete(j)%bf /= 0)then
            call add_element_pile(front_edge,j,mesh%arrete(j)%length,1)
        endif
    enddo
    
    do j = 1,nb_point
        if(mesh%point(j)%bf /= 0)then
            call add_element_pile(front_vertex,j,1._RP,1)
        end if
    end do
    
    if (iprint == 1)then 
        print*,''
        print*,'Update front...'
        print*,''
    endif
    
    ! Advance front method
    j = 1
    do while(associated(front_edge).and. j<=30000) ! 30000 is subjective (exit condition).
        
        tin = tin + 1._RP ! 1._RP is purely artificial.
                
        call optimalpoint(mesh,nb_point,nb_arrete,front_edge,front_vertex,new_point,&
        &    xgrid,ygrid,matdref,imetric,param,size_param,nx,ny,ierr,InputData,fgeom)
        if (ierr .ne. 0) then
            ierr = 401
            goto 9999
        endif
        
        call updatefront(mesh,nb_point,nb_arrete,nb_tri,front_edge,front_vertex,new_point,iface,ierr)
        if (ierr .ne. 0) then
            ierr = 402
            goto 9000
        endif
        
        if ((iwevol .and. ioMeshevol.eq.73) .or. (iwevolRemesh .and. ioMeshevol.eq.733)) then
            call write_tec(ioMeshevol,mesh,tin,nb_point,nb_tri,ierr) 
            if (ierr .ne. 0) then
                ierr = 403
                goto 9999
            endif
        endif
        j = j + 1
    enddo
        
    ! Affichage des piles 
    
    9000 continue
    
    9999 continue
             
    if(ierr/=0)then
        print*,"j = ",j
        write(*,99),ierr
    endif
    99 format('** error #',i3,' : AdvanceFront')

    if (iwfront == 1) then
        open(unit = 2,file="front_edge.dat")
    endif
    
    if (iwfront == 1) then
        open(unit = 4,file="front_vertex.dat")
    endif

    if (iwfront == 1) then
        print*,''
        print*,'Front Edge:'
        ptr => front_edge
        do while(associated(ptr))
            write(2,'(3i10)'),ptr%val,mesh%arrete(ptr%val)%iP(1),mesh%arrete(ptr%val)%iP(2)
            ptr => ptr%suiv
        enddo
        print*,''
        print*,'Front Vertex:'
        ptr => front_vertex
        do while(associated(ptr))
            write(4,41),ptr%val
            ptr => ptr%suiv
        enddo
    endif  
    41 format('ind = ',I4)
    
    ! Affichage triangles
    if (iprint == 1) then
        print*,'' 
        print*,'Nombre de triangle formé : ',nb_tri
        do j=1,nb_tri
            write(*,60),mesh%tri(j)%index,mesh%tri(j)%iP
        enddo
    endif
    60 format(' ind = ',i5,' P1 : ',i5,' P2 : ',i5,' P3 : ',i5)

    if (iwfront == 1) then 
        close(2)
        close(4)
    endif
    nullify(front_edge)
    nullify(front_vertex) 

end subroutine advancefront

end module MAdvanceFront 
