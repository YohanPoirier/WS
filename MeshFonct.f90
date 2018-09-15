module MeshFonct
use Constantes
use Structuresdonnees
use FonctionsCommunes
use GeomFonct
use GeomDisque
use GeomCylindre
use GeomPlan
use GeomSphere
use GeomCone
use GeomAxiSym
use GeomWigley
use MeshStruct
use Houle_RF
implicit none

! **************************************************************************
! Mesh Fonction
! **************************************************************************
! --------------------------------------------------------------------------
! Fonction de manipulation de liste chainee
! --------------------------------------------------------------------------
! Fonction de l'interface
! --------------------------------------------------------------------------
contains 

subroutine element_in_pile(struct,elem,res)
  !f2py integer*1, dimension(1000) :: struct
  type(pile_element),pointer :: struct
  integer,intent(in) :: elem
  logical,intent(out) :: res
  ! local
  !f2py integer*1, dimension(1000) :: ptr
  type(pile_element),pointer :: ptr
  res = .false.
  ptr => struct
  do while(associated(ptr) .and. .not.res)
    res = ptr%val==elem
    ptr => ptr%suiv
  enddo
end subroutine element_in_pile

subroutine mesh_segment(vertex1,vertex2,dl1,dl2,mesh,nb_point,nb_arete,bf,&
&                        lface,nface,Method_interpolation)
    
    !f2py integer*1, dimension(1000)    :: vertex1,vertex2
    type(MVertex),intent(in)            :: vertex1, vertex2
    real(rp),intent(in)                 :: dl1, dl2
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                     ! MGrid
    integer,intent(inout)               :: nb_point, nb_arete
    integer,intent(in)                  :: bf
    integer,dimension(2),intent(in)     :: lface
    integer,intent(in)                  :: nface
    integer,intent(in)                  :: Method_interpolation     ! = 0: Linear, = 1: parabolique.

    !f2py integer*1, dimension(1000)    :: P1,P2,P
    type(point)                         :: P1, P2, P
    !f2py integer*1, dimension(1000)    :: V
    type(vector)                        :: V
    !f2py integer*1, dimension(1000)    :: vt1,vt2
    type(MVertex)                       :: vt1, vt2
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge)                         :: edge
    real(rp)                            :: poid1, poid2, rn
    integer                             :: j, N, j1, j2, jsens
    real(rp),dimension(500)             :: RTab
    real(rp),dimension(3)               :: Tmp_Vect                 ! Temporary vector to create a new point.
    
    ! This subroutine creates the mesh for a segment.
    
    call assign_point_mvertex(P1,vertex1)
    call assign_point_mvertex(P2,vertex2)
    
    call assign_vector_coord(V,P2%coord - P1%coord)

    poid1 = dl1/V%length
    poid2 = dl2/V%length

    rn = 0
    N = 0
    do while (rn < 1.)
        N = N+1
        Rtab(N) = rn
        if(Method_interpolation == 0)then
            rn = rn + rn*poid2 + (1.-rn)*poid1
        elseif(Method_interpolation == 1)then
            rn = rn + (2.*rn-1)*(2.*rn-1)*poid1 + 4.*rn*(1.-rn)*poid2
        endif
    enddo
    
    if(N>=1)then
        if(abs(rn-1).lt.Epsilon)then
            N = N+1
            Rtab(N)=1.
        elseif(1.-RTab(N).lt.0.5*poid2)then
            Rtab(N) = 1.
        else
            N = N+1
            Rtab(N)=1.
        endif    
    endif
    
    vt1 = vertex1    
    
    j1 = 2
    j2 = N-1
    jsens = 1
    
    do j=j1,j2,jsens
        
        ! Point
        call Computation_eq_line(P1,P2,RTab(j),Tmp_Vect)
        call assign_point(P,Tmp_Vect) ! Vector -> Point
        
        call assign_mvertex_point(vt2,P) ! Point -> Mvertex
        vt2%bf = bf   
        vt2%nface = nface
        vt2%face(1:nface) = lface(1:nface)
        if(vt2%nedge .eq. 0)then
            vt2%edge = 0._RP
        end if
        
        call add_element(mesh%point,nb_point,vt2)        
        vt2%index = nb_point
        
        ! Edge
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = bf
        edge%nface=nface
        edge%face(1:nface) = lface(1:nface)
        call add_element(mesh,nb_arete,edge)
   
        vt1 = vt2         
    enddo
    
    vt2 = vertex2
    call edge_from_vertex(vt1,vt2,edge)
    edge%bf = bf
    edge%nface = nface
    edge%face(1:nface) = lface(1:nface)
    call add_element(mesh,nb_arete,edge)

end subroutine mesh_segment

subroutine mesh_segment_sym(vertex1,vertex2,dl1,dl2,mesh,nb_point,nb_arete,bf,&
&                        lface,nface)
    
    !f2py integer*1, dimension(1000) :: vertex1,vertex2
    type(MVertex),intent(in) :: vertex1, vertex2
    real(rp),intent(in) :: dl1, dl2
    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(inout) :: mesh
    integer,intent(inout) :: nb_point, nb_arete
    integer,intent(in) :: bf
    integer,dimension(2),intent(in) :: lface
    integer,intent(in) :: nface
    ! local
    !f2py integer*1, dimension(1000) :: P1,P2,P
    type(point) :: P1, P2, P
    !f2py integer*1, dimension(1000) :: V
    type(vector) :: V
    !f2py integer*1, dimension(1000) :: vt1,vt2
    type(MVertex) :: vt1, vt2
    !f2py integer*1, dimension(1000) :: edge
    type(MEdge) :: edge
    real(rp) :: poid1, poid2, rn
    integer :: j, N, j1, j2, jsens
    real(rp),dimension(500) :: RTab
  
    call assign_point_mvertex(P1,vertex1)
    call assign_point_mvertex(P2,vertex2)
    call assign_vector_coord(V,P2%coord - P1%coord)
 
    poid1 = 2.*dl1/V%length
    poid2 = 2.*dl2/V%length

    rn = 0
    N = 0
    do while (rn < 0.5)
        N = N+1
        Rtab(N) = rn
        rn = rn + rn*poid2 + (0.5-rn)*poid1
    enddo
    do while (rn < 1.)
        N = N+1
        Rtab(N) = rn
        rn = rn + (rn-0.5)*poid1 + (1.-rn)*poid2
    enddo
    if(N>=1)then
        if(abs(rn-1).lt.Epsilon)then
            N = N+1
            Rtab(N)=1.
        elseif(1.-RTab(N).lt.0.5*poid1)then
            Rtab(N) = 1.
        else
            N = N+1
            Rtab(N)=1.
        endif    
    endif
    
    vt1 = vertex1
    
    j1 = 2
    j2 = N-1
    jsens = 1
     
    do j=j1,j2,jsens
        
        ! Point
        call Computation_eq_line(P1,P2,RTab(j),P%coord)
        vt2 = P
        vt2%bf = bf   
        vt2%nface=nface
        vt2%face(1:nface) = lface(1:nface)
        if(vt2%nedge .eq. 0)then
            vt2%edge = 0._RP
        end if
        
        call add_element(mesh%point,nb_point,vt2)
        
        ! Edge
        vt2%index = nb_point
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = bf
        edge%nface=nface
        edge%face(1:nface) = lface(1:nface)
        call add_element(mesh,nb_arete,edge)
    
        vt1 = vt2
            
    enddo
      
    vt2 = vertex2
    call edge_from_vertex(vt1,vt2,edge)
    edge%bf = bf
    edge%nface = nface
    edge%face(1:nface) = lface(1:nface)
    call add_element(mesh,nb_arete,edge)
     
end subroutine mesh_segment_sym

subroutine mesh_periodic_boundary(vertex_1,vertex_2,dx1,dx2,mesh,nb_point,nb_arrete,&
&          ind_pt,ind_arr,a,b,i1,i2,iface)
  implicit none
  !f2py integer*1, dimension(1000) :: vertex_1,vertex_2
  type(MVertex),intent(in) :: vertex_1,vertex_2
  real(rp),intent(in) :: dx1,dx2
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  integer,intent(inout) :: nb_point,nb_arrete
  integer,dimension(*),intent(inout) :: ind_pt,ind_arr
  real(rp),intent(in) :: a,b
  integer,intent(inout) :: i1,i2
  integer,intent(in) :: iface
! local
  !f2py integer*1, dimension(1000) :: V
  type(vector) :: V
  !f2py integer*1, dimension(1000) :: P1,P2,P
  type(point) :: P1,P2,P
  !f2py integer*1, dimension(1000) :: vt1,vt2
  type(MVertex) :: vt1,vt2
  !f2py integer*1, dimension(1000) :: vt1p,vt2p
  type(MVertex) :: vt1p,vt2p
  !f2py integer*1, dimension(1000) :: edge
  type(MEdge) :: edge
  integer :: j,N,ind
  real(rp),dimension(100) :: RTab
  real(rp) :: poid1,poid2,rn,T
!  
  P1 = vertex_1
  P2 = vertex_2
  !V = [P2%coord(1)-P1%coord(1),P2%coord(2)-P1%coord(2),0.d0]
  call assign_vector_coord(V,[P2%coord(1)-P1%coord(1),P2%coord(2)-P1%coord(2),0.d0])
  
  poid1 = dx1/V%length
  poid2 = dx2/V%length
  rn = poid1
  N = 0
  do while (rn < 1.)
    N = N+1
    Rtab(N) = rn
    rn = rn + rn*poid2 + (1.-rn)*poid1
  enddo
  if(N>=1)then
    if(abs(rn-1).lt.Epsilon)then
      N = N+1
      Rtab(N)=1.
    elseif(1.-RTab(N).lt.poid2)then
      Rtab(N) = 1.
    else
      N = N+1
      Rtab(N)=1.
    endif    
  endif
    
  !N = int(V%length/dx)

!
  vt1 = vertex_1
  T = a*vt1%coord(2)+b
  vt1p = [vt1%coord(1)+T,vt1%coord(2),vt1%coord(3)]
  vt1p%bf = -12
  call add_element(mesh%point,nb_point,vt1p)
  vt1p%index = nb_point
  i1 = nb_point
  ind_pt(nb_point) = vt1%index
!  
  do j=1,N-1
  
    !P = eq_line(P1,P2,dble(j)/dble(N))
    call Computation_eq_line(P1,P2,RTab(j),P%coord)
    vt2 = P
    vt2%bf = -11    
    vt2%nface=1
    vt2%face(1) = iface
    call add_element(mesh%point,nb_point,vt2)
    vt2%index = nb_point
    call edge_from_vertex(vt1,vt2,edge)
    edge%bf = -11
    edge%nface=1
    edge%face(1)=iface
    call add_element(mesh,nb_arrete,edge)
    ind = nb_arrete
    vt1 = vt2
    
    T = a*vt2%coord(2)+b
    vt2p = [vt2%coord(1)+T,vt2%coord(2),vt2%coord(3)]
    vt2p%bf = -12
    call add_element(mesh%point,nb_point,vt2p)
    vt2p%index = nb_point
    ind_pt(nb_point) = vt2%index
    call edge_from_vertex(vt2p,vt1p,edge) ! sens inverse de l'arrete precedente
    edge%bf = -12
    call add_element(mesh,nb_arrete,edge)
    ind_arr(nb_arrete) = ind
    vt1p = vt2p
        
  enddo
  
  vt2 = vertex_2
  call edge_from_vertex(vt1,vt2,edge)
  edge%bf = -11
  edge%nface=1
  edge%face(1)=iface
  call add_element(mesh,nb_arrete,edge)
  ind = nb_arrete
  

  T = a*vt2%coord(2)+b
  vt2p = [vt2%coord(1)+T,vt2%coord(2),vt2%coord(3)]
  vt2p%bf = -12
  call add_element(mesh%point,nb_point,vt2p)
  vt2p%index = nb_point
  i2 = nb_point
  ind_pt(nb_point) = vt2%index
  call edge_from_vertex(vt2p,vt1p,edge)
  edge%bf = -12
  call add_element(mesh,nb_arrete,edge)
  ind_arr(nb_arrete) = ind
  
end subroutine mesh_periodic_boundary

subroutine mesh_polyline(tab,ntab,is_close,dx,mesh,nb_point,nb_arrete,bf,iajout,ierror)
    
    !f2py integer*1, dimension(1000)    :: tab
    type(point),dimension(*),intent(in) :: tab                  ! Table of the pointers toward the intersection points.
    integer,intent(in)                  :: ntab                 ! Number of intersection curve.
    logical,intent(in)                  :: is_close             ! True if the intersection curve is closed, false otherwise.
    real(rp),intent(in)                 :: dx                   ! Panel size.
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                 ! MGrid.
    integer,intent(inout)               :: nb_point,nb_arrete   ! Number of points and edges.
    integer,intent(in)                  :: bf                   ! Boudary flag.
    integer,intent(in)                  :: iajout               ! = 1 in the subroutine mesh_inter and 0 in mesh1D.
    integer,intent(inout)               :: ierror               ! Error flag.
    
    !f2py integer*1, dimension(1000)    :: P1,P2,P
    type(point)                         :: P1,P2,P              ! Points.
    !f2py integer*1, dimension(1000)    :: vertex_1,vertex_2
    type(MVertex)                       :: vertex_1,vertex_2    ! Vertexes.
    !f2py integer*1, dimension(1000)    :: V
    type(vector)                        :: V                    ! Vector.
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge)                         :: edge                 ! Edge.
    integer                             :: j,k,N,ip1            ! Loop parameters, ip1 = index of the first point.
    integer                             :: nface                ! Number of the face.
    integer,dimension(nvmax)            :: iface                ! ?
    
    ! This subroutine creates the mesh intersection (from tab) curve when creating the mesh or remeshing.
    
    ip1 = -9
    
    ! At least two points on the intersection curve.
    if(ntab.ne.1)then
        
        do j=1,ntab-1
            
            P1 = tab(j)
            P2 = tab(j+1)
                        
            call assign_vector_coord(V,[P2%coord(1)-P1%coord(1),P2%coord(2)-P1%coord(2),P2%coord(3)-P1%coord(3)]) ! Coord -> Vector
            N = int(V%length/dx) 
            
            if (j==1) then
                ! Points to Mvertex.
                call assign_mvertex_point(vertex_1,P1)
            
                ! bf
                vertex_1%bf = bf
            
                ! Add the vertex or not.
                if(iajout==1)then ! The vertex is added to MGrid.
                    call add_element(mesh%point,nb_point,vertex_1)
                    vertex_1%index = nb_point
                else
                    vertex_1%index = P1%flag
                endif
            
                ip1 = nb_point ! Index of the 1st point of the intersection curve
            else
                vertex_1 = vertex_2 ! Previous vertex
            endif 
        
            call common_int(P1%face,P1%nface,P2%face,P2%nface,iface,nface,ierror)  
        
            ! Second discretization of the polyline if the first one with tab is not enought.
            do k = 1,N-1
            
                ! Point
                call Computation_eq_line(P1,P2,dble(k)/dble(N),P%coord) ! Information about face and edge from P1 and P2 in P ???
                        
                vertex_2 = P
                print*,"Mesh_polyline: maybe a problem with and vertex2%edge because of the operator overloading."
            
                vertex_2%bf = bf
                call add_element(mesh%point,nb_point,vertex_2)
                vertex_2%index = nb_point
            
                ! Edge
                call edge_from_vertex(vertex_1,vertex_2,edge)
                edge%bf = bf
                edge%nface = nface
                edge%face(1:nface) = iface(1:nface)
                call add_element(mesh,nb_arrete,edge)
                vertex_1 = vertex_2
            enddo
            
            ! Last point if the intersection curve is closed.
            if(j==ntab-1 .and. is_close) then
                vertex_2 = mesh%point(ip1)      
            else
            
                call assign_mvertex_point(vertex_2,P2)
            
                ! if %nedge = 0, then %edge is not updated by assign_mvertex_point so the previous values are not modified.
                if(P2%nedge .eq. 0)then
                    vertex_2%edge = 0._RP
                end if
                
                ! bf
                vertex_2%bf = bf
                
                ! Add the vertex or not.
                if(iajout==1)then ! The vertex is added to MGrid.
                    call add_element(mesh%point,nb_point,vertex_2)
                    vertex_2%index = nb_point
                else
                    vertex_2%index = P2%flag
                endif      
            endif
        
            call edge_from_vertex(vertex_1,vertex_2,edge)
            edge%bf = bf
            edge%nface = nface
            edge%face(1:nface) = iface(1:nface)
            call add_element(mesh,nb_arrete,edge)
        
        end do 
        
    else if(ntab.eq.1)then ! The intersection curve is only one point.
        
        ! Unique point.
        P1 = tab(1)
        
        ! Points to Mvertex.
        call assign_mvertex_point(vertex_1,P1)
        
        ! bf
        vertex_1%bf = bf
        
        ! Add the vertex or not.
        if(iajout==1)then ! The vertex is added to MGrid.
            call add_element(mesh%point,nb_point,vertex_1)
            vertex_1%index = nb_point            
        endif
        
    end if
    
end subroutine mesh_polyline

subroutine mesh_polyline2(vertex_start,vertex_end,np,P_loc,rep,dx,mesh,nb_point,nb_arrete,bf,iedge,ierror)
    
    !f2py integer*1, dimension(1000)    :: vertex_start,vertex_end
    type(MVertex),intent(in)            :: vertex_start, vertex_end
    integer,intent(in)                  :: np
    !f2py integer*1, dimension(1000)    :: P_loc
    real(rp),dimension(np,3),intent(in) :: P_loc
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep
    real(rp),intent(in)                 :: dx
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh
    integer,intent(inout)               :: nb_point, nb_arrete
    integer,intent(in)                  :: bf
    integer,intent(in)                  :: iedge
    integer,intent(inout)               :: ierror
    
    !f2py integer*1, dimension(1000)    :: vertex1,vertex2
    type(MVertex)                       :: vertex1, vertex2
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge)                         :: edge
    integer                             :: j, nface, ip1, ip2
    integer,dimension(nvmax)            :: face
    real(rp)                            :: l_end, l, dl, lf, l_tot
    real(rp),dimension(3)               :: M, M1, M2, P1, P2, Vdir
    real(rp),dimension(np,3)            :: P_cart
    integer,dimension(2)                :: ips, ipe
    logical                             :: loop
    
    ! This subroutine creates the mesh of the polylines which are under the sea level in mesh1D.

    ierror = 0
    
    call cart2loc(vertex_start%coord,rep,M1)
    call cart2loc(vertex_end%coord,rep,M2)
    
    
    
    call find_voisin(P_loc(:,2),M1(2),ips)
    call find_voisin(P_loc(:,2),M2(2),ipe)
    
    if(ips(1)==ips(2))then
        if(ips(1)<np)then
            ips(2) = ips(1)+1
        else
            goto 9999
        endif
    endif

    ip1 = ips(1)
    ip2 = ips(2)

    if(ipe(1)==ipe(2))then
        if(ipe(1)>1)then
            ipe(1) = ipe(1)-1
        else
            goto 9999
        endif
    endif
    
    call common_int(vertex_start%face,vertex_start%nface,&
&                   vertex_end%face,vertex_end%nface,face,nface)

    do j=1,np
        call loc2cart(P_loc(j,1:3),rep,P_cart(j,1:3))
    enddo
    
    P1 = P_cart(ips(1),1:3)
    P2 = P_cart(ips(2),1:3)
    Vdir = P2-P1
    Vdir = Vdir / norm2(Vdir)
    
    vertex1 = vertex_start
    vertex1%nedge = 1
    vertex1%edge(1) = iedge
    vertex2 = vertex1
    
    M = vertex_start%coord
    lf = norm2(P2-M)
    l = 0._rp
    l_end = norm2(vertex_end%coord-P_cart(ipe(1),:))
    l_tot = lf + l_end 
    do j=ip2,ipe(1)-1
        l_tot = l_tot + norm2(P_cart(j+1,:)-P_cart(j,:))
    enddo
    dl = min(lf,dx-l)

    l_tot = l_tot-dl
    
    loop = .true.
    do while (loop)
        
        M = M + dl*Vdir
        l=l+dl
        lf = lf-dl
        
        if(abs(lf).lt.Epsilon .and. ip1<np .and. ip2<np)then
            M = P2
            ip1 = ip1+1
            ip2 = ip2+1
            P2 = P_cart(ip2,:)
            P1 = P_cart(ip1,:)
            Vdir = P2-P1
            lf = norm2(P2-P1)
            Vdir = Vdir/lf
        endif
         
        if(abs(l-dx).lt.Epsilon)then
            vertex2%coord = M
            
            
            call add_element(mesh%point,nb_point,vertex2)
            vertex2%index = nb_point
            call edge_from_vertex(vertex1,vertex2,edge)
            edge%bf = bf
            edge%nface = nface
            edge%face(1:nface) = face(1:nface)
            call add_element(mesh,nb_arrete,edge)
            vertex1 = vertex2
            l = 0._rp
        endif
        
        dl = min(lf,dx-l)
        
        l_tot = l_tot-dl
        if(l<0.5_rp*dx .and. l_tot<dx .or. l>0.5_rp*dx .and. l_tot<=l)then
            loop = .false.
        endif
        
    enddo
    
    ! Added last edge
    call edge_from_vertex(vertex1,vertex_end,edge)
    edge%bf = bf
    edge%nface = nface
    edge%face(1:nface) = face(1:nface)  
    call add_element(mesh,nb_arrete,edge)
    
    9999 continue
        if(ierror/=0)then
            write(*,99) ierror
        endif
    99 format('** error #',i3,' : cannot mesh polyline')
    
end subroutine mesh_polyline2 

subroutine mesh_arc_cercle(vertex_1,vertex_2,dx,mesh,nb_point,&
&          nb_arrete,bf,cercle,lface,nface,isens,ierr,igeom,HorizontalCyl)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                 ! Mesh.
    integer,intent(inout)               :: nb_point,nb_arrete   ! Number of points and edges.
    !f2py integer*1, dimension(1000)    :: vertex_1,vertex_2
    type(MVertex),intent(in)            :: vertex_1,vertex_2    ! Vertexes.
    !f2py integer*1, dimension(1000)    :: cercle
    type(GCercle),intent(in)            :: cercle               ! Circle.
    real(rp),intent(in)                 :: dx                   ! dx = dx2 for the floaters, dx1 for the domain.
    integer,intent(in)                  :: bf                   ! boundary flag.
    integer,intent(in)                  :: isens                ! Sens of the path of the nodes (actually isens is fixed to 1).
    integer,dimension(*),intent(in)     :: lface                !
    integer,intent(in)                  :: nface                ! Number of faces.
    integer,intent(inout)               :: ierr                 ! Error flags.
    integer,intent(in)                  :: igeom                ! igeom = 0 for the domain, 1 for the floaters.
    logical,intent(in)                  :: HorizontalCyl        ! Horizontal cylinder modification.
    
    real(rp)                            :: L,angle,dtheta       !
    real(rp),dimension(3)               :: y1,y2,yl             ! Points.
    !f2py integer*1, dimension(1000)    :: vt1,vt2
    type(MVertex)                       :: vt1,vt2              ! Vertexes.
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge)                         :: edge                 ! Edge.
    integer                             :: N,j                  !
    real(rp),dimension(3)               :: Tmp_Vect             ! Temporal vector.
    
    ! This subroutine creates the mesh of a part of a circle.
    
    ierr = 0
    
    call cart2sph(vertex_1%coord,cercle%repere,y1)
    call cart2sph(vertex_2%coord,cercle%repere,y2)
    
    angle = y2(2)-y1(2)    
    if(isens == -1) then
        angle = -2.*PI + angle
    endif
    
    ! Only for the floaters
    if(igeom.eq.1 .and. HorizontalCyl)then 
        print*,"mesh_arc_cercle: definition of angle modified for the horizontal cylinders."
        angle = TWOPI - abs(angle) ! Horizontal cylinder
    end if
      
    L = abs(angle) * cercle%rayon
    N = int(L/dx)
    
    if(N.eq.0)then
        N = 1
    end if
    
    if(N.eq.0)then
        print*,"mesh_arc_cercle: N is equal to 0."
        print*,"angle = ",angle
        print*,"L = ",L
        ierr = 100
        go to 10
    else
        dtheta = angle/dble(N)
    end if
    
    vt1 = vertex_1
    yl = y1
    do j=1,N-1
        
        ! Point
        yl(2) = y1(2)+j*dtheta      
        call sph2cart(yl,cercle%repere,Tmp_Vect)
        call assign_mvertex_coord(vt2,Tmp_Vect) ! Vect -> MVertex
        vt2%bf = bf
        vt2%nface = nface
        vt2%face(1:nface) = lface(1:nface)
        if(vt2%nedge .eq. 0)then
            vt2%edge = 0._RP
        end if
        
        call add_element(mesh%point,nb_point,vt2)
        
        ! Edge
        vt2%index = nb_point
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = bf
        edge%nface = nface
        edge%face(1:nface) = lface(1:nface)
        call add_element(mesh,nb_arrete,edge)
        vt1 = vt2
    enddo
  
    vt2 = vertex_2
    call edge_from_vertex(vt1,vt2,edge)
    edge%bf = bf
    edge%nface = nface
    edge%face(1:nface) = lface(1:nface)
    call add_element(mesh,nb_arrete,edge)
    
    10    continue
    
end subroutine mesh_arc_cercle
  
subroutine mesh_cercle(mesh,nb_point,nb_arrete,dx,&
&          cercle,lface,nface,isens,bf,ierr)
    
    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(inout) :: mesh
    !f2py integer*1, dimension(1000) :: nb_point,nb_arrete
    integer,intent(inout) :: nb_point,nb_arrete
    real(rp),intent(in) :: dx
    !f2py integer*1, dimension(1000) :: cercle
    type(GCercle),intent(in) :: cercle
    integer,dimension(*),intent(in) :: lface
    integer,intent(in) :: nface
    integer,intent(in) :: isens
    integer,intent(in) :: bf
    integer,intent(inout) :: ierr
    
    !f2py integer*1, dimension(1000) :: vt1,vt2,vt0
    type(Mvertex) :: vt1,vt2,vt0
    !f2py integer*1, dimension(1000) :: edge
    type(MEdge) :: edge
    real(rp),dimension(3) :: yl,y0
    real(rp) :: L,dtheta
    integer :: N,j,Nmax
    real(rp),dimension(3) :: Tmp_Vect
    
    ! This subroutine creates the mesh of a circle.
    
    ierr = 0
    vt1%face = 0
    vt2%face = 0
    vt0%face = 0
    edge%face = 0

    y0 = [cercle%rayon,0.d0,PI/2.]

    if(Symmetry)then
        L = PI*cercle%rayon
        N = int(L/dx)
        dtheta = PI/dble(N)
        Nmax = N
    else
        L = 2.*PI*cercle%rayon
        N = int(L/dx)
        dtheta = 2.*PI/dble(N)
        Nmax = N-1
    endif
    if(isens==-1)then
        dtheta = -dtheta
    endif

    call sph2cart(y0,cercle%repere,vt0%coord)
    vt0%nface = nface
    vt0%face(1:nface) = lface(1:nface)
    vt0%bf = bf
    if(Symmetry) then 
        vt0%bf = bf
        vt0%nface = vt0%nface+1
        vt0%face(vt0%nface)=-1
    endif
    call add_element(mesh%point,nb_point,vt0)
    vt0%index = nb_point

    vt1 = vt0
    yl = y0 
    do j=1,Nmax
        
        ! Point
        yl(2) = y0(2)+j*dtheta        
        call sph2cart(yl,cercle%repere,Tmp_Vect)
        call assign_mvertex_coord(vt2,Tmp_Vect) ! Vect -> MVertex
        vt2%bf = bf
        vt2%nface = nface
        vt2%face(1:nface) = lface(1:nface)
        if(vt2%nedge .eq. 0)then
            vt2%edge = 0._RP
        end if
        
        call add_element(mesh%point,nb_point,vt2)
        
        ! Edge
        vt2%index = nb_point
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = bf
        edge%nface = nface
        edge%face(1:nface) = lface(1:nface)
        call add_element(mesh,nb_arrete,edge)
        vt1 = vt2
    enddo
    
    ! Closing of the mesh
    if(.not.Symmetry)then
        vt2 = vt0
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = bf
        edge%nface = nface
        edge%face(1:nface) = lface(1:nface)
        call add_element(mesh,nb_arrete,edge)
    else
        mesh%point(nb_point)%bf = bf
        mesh%point(nb_point)%nface = mesh%point(nb_point)%nface+1
        mesh%point(nb_point)%face(mesh%point(nb_point)%nface) = -1
    endif

end subroutine mesh_cercle  

subroutine extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
&          ind_pt,ind_arr,iface,ierror)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(in)              :: mesh                 ! Mgrid
    integer,intent(in)                  :: nb_arete             ! Number of edges
    !f2py integer*1, dimension(1000)    :: mesh2
    type(MGrid),intent(inout)           :: mesh2                ! New Mgrid
    integer,intent(inout)               :: nb_point2,nb_arete2  ! New number of points and edges
    integer,dimension(*),intent(inout)  :: ind_pt,ind_arr       ! 
    integer,intent(in)                  :: iface                ! Index of the face
    integer,intent(inout)               :: ierror               ! Error flag
        
    !f2py integer*1, dimension(1000)    :: edge
    type(MEdge)                         :: edge                 ! Edge
    !f2py integer*1, dimension(1000)    :: vertex
    type(MVertex),dimension(2)          :: vertex               ! Vertex
    integer,dimension(100)              :: aux
    integer                             :: j,jj,k               ! Loop parameters
    integer                             :: ind,nface,naux,ip
    integer                             :: bf                   ! Boundary flag
    logical                             :: bool 
    
    ! This subroutine extracts the points of mesh (from mesh1D) which match with iface.
    
    if(idebug/=0) print*," MP_extract_point_mesh .."
  
    ind_pt(1:-1) = -9
    ind_arr(1:-1) = -9
    
    do j=1,nb_arete
        edge = mesh%arrete(j)
        bf = edge%bf
        ind = edge%index
        nface = edge%nface
        call common_int(edge%face(1:nface),nface,[iface],1,aux,naux,ierror)
        bool = naux==1
        if(bool)then ! Face in common
            do jj=1,2
                ip = edge%iP(jj)
                k = 1
                do while(ip /= ind_pt(k) .and. k <= nb_point2)
                    k=k+1
                enddo
                if(k <= nb_point2) then
                    edge%iP(jj) = k
                    vertex(jj) = mesh2%point(k)
                else
                    vertex(jj) = mesh%point(ip)
                    if(vertex(jj)%bf.le.0 .and. vertex(jj)%bf.ge.-8)then
                        if(vertex(jj)%bf.eq.-4)then
                            vertex(jj)%bf=4
                        else
                            vertex(jj)%bf=1
                        endif
                    else
                        vertex(jj)%bf=bf
                    endif
                    call add_element(mesh2%point,nb_point2,vertex(jj))
                    vertex(jj)%index = nb_point2
                    ind_pt(nb_point2) = ip
                endif
            enddo
            
            call edge_from_vertex(vertex(1),vertex(2),edge)
            
            if(bf.le.0 .and. bf.ge.-8)then
                if(edge%bf.eq.-4)then
                    edge%bf = 4
                else
                    edge%bf = 1
                endif        
            else
                edge%bf = bf
            endif
            
            call add_element(mesh2,nb_arete2,edge)
            ind_arr(nb_arete2) = ind
        endif
    enddo 
    
end subroutine extract_point_mesh

subroutine change_edge(mesh,edge,vertex_1,vertex_2,bf)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  !f2py integer*1, dimension(1000) :: edge
  type(MEdge),intent(inout) :: edge
  !f2py integer*1, dimension(1000) :: vertex_1,vertex_2
  type(MVertex),intent(inout) :: vertex_1,vertex_2
  integer,intent(in) :: bf
! local
  integer :: ind,ind2,ip1,ip2
  integer :: ivt1,ivt2,n
  !f2py integer*1, dimension(1000) :: V
  type(vector) :: V
  integer :: jj,k,k2
  !f2py integer*1, dimension(1000) :: P
  type(MVertex) :: P
  
  
  ind = edge%index
  ip1 = edge%iP(1)
  ip2 = edge%iP(2)
  ivt1 = vertex_1%index
  ivt2 = vertex_2%index 
  
! Actualisation arrete  
  edge%iP(1) = ivt1
  edge%iP(2) = ivt2
  edge%P1 = vertex_1%coord
  edge%P2 = vertex_2%coord
  
  call assign_vector_coord(V,vertex_2%coord - vertex_1%coord)
  edge%length = V%length
  edge%dir = V%coord/V%length
  edge%bf = bf 

! Actualisation noeuds  
  do jj=1,2
    P = mesh%point(ip1)
    n = P%nv_arrete
    P%nv_arrete = n-1
    k2 = 0
    do k=1,n
      ind2 = P%va(k)
      if(ind2 /= ind) then
        k2 = k2+1
        P%va(k2) = ind2
      endif
    enddo
    P%va(k2+1) = -9
  
    n = P%nv_point
    P%nv_point = n-1
    k2 = 0
    do k=1,n
      ind2 = P%vp(k)
      if(jj== 1 .and. ind2 /= ip2 .or. jj==2 .and. ind2 /= ip1 ) then
        k2=k2+1
        P%vp(k2)=ind2      
      endif
    enddo
    P%vp(k2+1)=-9
  enddo

! Nouveaux noeuds
  n = vertex_1%nv_arrete
  vertex_1%nv_arrete=n+1
  vertex_1%va(n+1)=ind
  
  n = vertex_1%nv_point
  vertex_1%nv_point = n+1
  vertex_1%vp(n+1)=ivt2
  mesh%point(ivt1) = vertex_1
  
  n = vertex_2%nv_arrete
  vertex_2%nv_arrete=n+1
  vertex_2%va(n+1)=ind
  
  n = vertex_2%nv_point
  vertex_2%nv_point = n+1
  vertex_2%vp(n+1)= ivt1
  mesh%point(ivt2) = vertex_2
!  
end subroutine change_edge  

subroutine get_index(mesh,nb_arrete,vertex_1,vertex_2,ind)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(in) :: mesh
  integer,intent(in) :: nb_arrete
  !f2py integer*1, dimension(1000) :: vertex_1,vertex_2
  type(MVertex),intent(in) :: vertex_1,vertex_2
  integer,intent(out) :: ind
! local
  !f2py integer*1, dimension(1000) :: edge
  type(MEdge) :: edge
  integer :: ivt1,ivt2
  integer :: j
  logical :: bool
  integer :: ierror
  
  ind = -9
  ierror = 0
  bool = .false.
  ivt1 = vertex_1%index
  ivt2 = vertex_2%index
  
  j=1
  do while (j <= nb_arrete .and. .not.bool)
    edge = mesh%arrete(j)
    bool = edge%iP(1)==ivt1 .and. edge%iP(2)==ivt2
    bool = bool .or. edge%iP(2)==ivt1 .and. edge%iP(1)==ivt2
    j = j+1
  enddo
  
  if(bool)then
    ind = j-1
  else
    ierror = 1
    goto 9999
  endif
  
9999 continue
  if (ierror /= 0) then
    print*,'** error : get_index'
    print*,'**       : vertex_1 : ',vertex_1%index
    print*,'**       : vertex_2 : ',vertex_2%index
  endif
  
end subroutine get_index

subroutine optimalpoint(mesh,nb_point,nb_arrete,front_edge,front_vertex,new_point,xgrid,ygrid,matdref,imetric,param,size_param,nxx,nyy,ierr,InputData,fgeom)
    
    !f2py integer*1, dimension(1000)                    :: mesh
    type(MGrid)                                         :: mesh                                                             ! MGrid.
    integer,intent(in)                                  :: nxx,nyy                                                          ! Size of xgrid,ygrid and matdref.
    integer,intent(inout)                               :: nb_point,nb_arrete                                               ! Number of points and edges for the present surface.
    !f2py integer*1, dimension(1000)                    :: front_edge
    type(pile_element),pointer                          :: front_edge                                                       ! Edges of the mesh front.
    !f2py integer*1, dimension(1000)                    :: front_vertex
    type(pile_element),pointer                          :: front_vertex                                                     ! Nodes of the mesh front.
    !f2py integer*1, dimension(1000)                    :: new_point    
    type(MVertex)                                       :: new_point                                                        ! New point.
    real(rp),dimension(nxx),intent(in)                  :: xgrid                                                            ! Cartesian grid.
    real(rp),dimension(nyy),intent(in)                  :: ygrid                                                            ! Cartesian grid.
    real(rp),dimension(nxx,nyy,3),intent(in)            :: matdref                                                          ! Reference length for each cell of the Cartesian grid.
    integer,intent(in)                                  :: imetric                                                          ! Metric (0 for a cylinder, a disc and a plan, 1 for a sphere, 2 for a cone, 3 for an axisym and 4 for a Wigley hull).
    integer,intent(in)                                  :: size_param                                                       ! Size of param.
    real(rp),dimension(size_param),intent(in)           :: param                                                            ! 8 param for a cone, 2*axisym%npoint+1 for a axisym, 1 for a cylinder, a disc, a sphere, a plan and a Wigley hull.
    integer,intent(inout)                               :: ierr                                                             ! Error flag.
    !f2py integer*1, dimension(1000)                    :: InputData
    type(InputDataStruct),intent(in)                    :: InputData                                                        ! Input data. 
    !f2py integer*1, dimension(1000)                    :: fgeom
    type(type_geom),intent(in)                          :: fgeom                                                            ! Geometry.
        
    !f2py integer*1, dimension(1000)                    :: active_edge
    type(MEdge)                                         :: active_edge                                                      ! Current adge of the mesh front.
    !f2py integer*1, dimension(1000)                    :: Q
    type(MVertex)                                       :: Q                                                                ! Possible point to test.
    real(rp),dimension(3)                               :: P_opt,P1,P2,P                                                    ! Optimal point (P_opt), endless points of the edges in the intersection test (P1 and P2), point of the Cartesian grid (P).
    real(rp),dimension(3)                               :: A,B,M,M0                                                         ! A and B are the endpoints of the current edge of the mesh front, M and M0 are the middle of the of the current edge AB.
    real(rp),dimension(3)                               :: dir,vect                                                         ! dir: surely uesless, vect is used in the computation of the metric.
    real(rp),dimension(3)                               :: dir2                                                             ! Direction used in the computation of the metric.
    real(rp),dimension(2)                               :: N_AB                                                             ! Normal vector to the current edge of the mesh front.
    !f2py integer*1, dimension(1000)                    :: near_point,ptr_np
    type(near_point_cell),pointer                       :: near_point,ptr_np                                                ! Pointer to store all the possible nodes (near_point) and a useless pointer (ptr_np).
    !f2py integer*1, dimension(1000)                    :: ptr_fv
    type(pile_element),pointer                          :: ptr_fv                                                           ! Pointer for seaching the nearby point.
    !f2py integer*1, dimension(1000)                    :: ptr_ed
    type(pile_element),pointer                          :: ptr_ed                                                           ! Pointer to check all the edges of the mesh front in the intersection test.
    real(rp)                                            :: delta_1,delta_2                                                  ! Delta_1 and delta_2.
    real(rp)                                            :: d_ref                                                            ! Delta_ref.
    real(rp)                                            :: d_ref1                                                           ! Delta_ref1
    real(rp)                                            :: d,dd,dist,ps_2D, ps_3D                                           ! Distances for searching the new points and dot product results.
    real(rp)                                            :: dist2,dist3                                                      ! Distances for searching the new points.
    real(rp)                                            :: angle_maxA, angle_maxB                                           ! Limit angles.
    integer                                             :: j,k,i1,i2,j1,j2,n                                                ! Loop parameters.
    integer                                             :: index_arr,index_ed,index_A,index_B,index_Q                       ! Indexes of edges.
    integer                                             :: index_1,index_2                                                  ! Indexes of the edge to check in the intersection test.
    integer                                             :: ierror                                                           ! Error flag.
    logical                                             :: is_tri,is_inter1,is_inter2,is_inter,is_prox,is_intnode,is_angle  ! Flags to know is a triangle is admissible (false) or not (true).
    integer                                             :: count,it                                                         !
    integer,dimension(2)                                :: ivoisin,jvoisin                                                  !
    real(rp),dimension(2,2)                             :: mat2                                                             !
    integer                                             :: nx,ny                                                            ! Number of points in x and y in the Cartesian grid.
    real(rp)                                            :: xmin,xmax,ymin,ymax                                              ! Coordinates of the Cartesian grid.
    logical                                             :: bool                                                             ! Boolean.
    real(rp),dimension(3)                               :: NormaleA, NormaleB, Normale, N_3D                                ! Normal vectors.
    real(rp),dimension(3)                               :: vect_product_1                                                   ! Cross product result.
    logical                                             :: PcloserM                                                         ! Searching for an optimal point closer to M (dist3 -> 0.5*dist3) according to the equation (4.12) of CC.
    integer                                             :: NumPoints                                                        ! Number of points which can be used to create a new panel.
    real(rp)                                            :: iter
    
    integer,parameter                                   :: ndiv = 2                                                         ! Number of division (useless).
    
    ! This subroutine searches the optimal point, its neighboors, then creates a new panel and finally tests it.
    
    iter = 0
    
    ! Initialisation flag
    is_intnode = .false.
    ierror = 0
    count = 0
    new_point%coord = -9.
    new_point%index = -9
    PcloserM = .false.
    NumPoints = 0
    
    ! Cartesian grid parameters
    nx = size(xgrid)
    xmin = xgrid(1)
    xmax = xgrid(nx)
    
    ny = size(ygrid)
    ymin = ygrid(1)
    ymax = ygrid(ny)  
  
    ! Arrete du front courant
    index_arr = front_edge%val
    active_edge = mesh%arrete(index_arr)  
    
    ! Points qui delimitent l'arrete
    index_A = active_edge%iP(1)
    index_B = active_edge%iP(2)
    A = mesh%point(index_A)%coord
    B = mesh%point(index_B)%coord
    
    ! Normal to the edge
    N_AB = [A(2)-B(2),B(1)-A(1)] ! Another normal would be [B(2)-A(2),A(1)-B(1)]
    N_AB = N_AB / norm2(N_AB)
    
    3 continue
    
    ! Centre of the edge (M in the PhD of CC (p 67))
    M = 0.5_rp*(A+B)
    M0 = M

    if(imetric==3)then ! imetric = 3 for the axisym geometry
        n = int(param(1))
        call normal_axisym(A,n,param(2:n+1),param(n+2:2*n+1),NormaleA)
        call normal_axisym(B,n,param(2:n+1),param(n+2:2*n+1),NormaleB)
        Normale = 0.5_rp*(NormaleA+NormaleB)
        call Computation_vect_product(B-A,-Normale,N_3D)
        N_3D = N_3D/norm2(N_3D)
    endif
    
    ! Calcul distance de reference  
    P = M ! Centre of the edge
    do it = 1,2
        
        ! Closest cells of the Cartesian grid in the x direction.
        call find_voisin(xgrid,P(1),ivoisin)
        i1 = ivoisin(1) ; i2 = ivoisin(2)
        
        ! Closest cells of the Cartesian grid in the y direction.
        call find_voisin(ygrid,P(2),jvoisin)
        j1 = jvoisin(1) ; j2 = jvoisin(2)
        
        do j=1,2
            do k=1,2
                mat2(j,k) = matdref(ivoisin(j),jvoisin(k),1)
            enddo
        end do
        
        if(it==2) d_ref1 = d_ref
        
        ! Bilinear interpolation of the reference size.
        call int_bilin_cart(P(1:2),d_ref,xgrid(i1),xgrid(i2),ygrid(j1),ygrid(j2),&
        &    mat2(1,1),mat2(2,1),mat2(1,2),mat2(2,2),ierror)   
        if(ierror/=0)then
            ierror = 100
            goto 9999
        endif
    
        !## FIXME
        if(imetric==3) d_ref = InputData%dx2(1) ! Axisym geometry
        !##
    
        if(it==2) d_ref = 0.5*(d_ref1+d_ref)
        
        d_ref = min(d_ref,1.5*active_edge%length)

        if (0.55*active_edge%length < d_ref .and. 2.*active_edge%length > d_ref) then
            delta_1 = d_ref
        elseif (0.55*active_edge%length < d_ref) then
            delta_1 = 0.55*d_ref
        elseif (2.*active_edge%length > d_ref) then
            delta_1 = 2.*d_ref
        else
            ierror = 2            
            goto 9999
        endif
        
        ! Equation (4.9) of CC
        if(not(PcloserM))then
            dist3 = 0.5*sqrt(4.*delta_1*delta_1-d_ref*d_ref)
        elseif(PcloserM)then
            dist3 = 0.5*0.5*sqrt(4.*delta_1*delta_1-d_ref*d_ref)
        end iF
        
        ! Searching a new point in the disc of radius dist3.
        call estimate_position(A,B,dist3,imetric,param,P,ierror,fgeom) ! P is modified here.
        
        delta_2 = delta_1

        if(imetric/=2)then  
            if(P(1).lt.xmin-Epsilon .or. P(1).gt.xmax+Epsilon .or.&
        &  P(2).lt.ymin-Epsilon .or. P(2).gt.ymax+Epsilon) then
                goto 10
            endif
        endif
        
    enddo
    
    10 continue
    P_opt = P
    
    100 continue
    
    ! Liste des points potentiellement admissibles
    nullify(near_point)
    
    ptr_fv => front_vertex
    do while (associated(ptr_fv))
        index_Q = ptr_fv%val
        Q = mesh%point(index_Q)
                
        if(index_Q /= index_A .and. index_Q /= index_B) then
    
            d = norm2(P_opt(1:2)-Q%coord(1:2))
            ps_2D = (Q%coord(1)-M(1))*N_AB(1) + (Q%coord(2)-M(2))*N_AB(2)
            bool = d < NormOverDref*delta_2 .and. ps_2D > Epsilon ! PYW
            
            if(imetric == 3 .and. bool)then ! Axisym geometry                
                call norm_metric(Q%coord,P_opt,imetric,param,d)
                ps_3D = dot_product(Q%coord-M,N_3D)
                bool = d < NormOverDref*delta_2 .and. ps_3D > Epsilon ! bool = d < 0.8*delta_2 .and. ps_3D > Epsilon
            endif
            
            if(bool)then                
                call add_vertex_pile(near_point,Q,d)
                NumPoints = NumPoints + 1
            endif
    
        endif
        ptr_fv => ptr_fv%suiv
    end do
    
    ! Test if the angle between adjacent fronts is below 75 degrees, if yes, creation of a new triangle. 
    if(imetric==3)then ! Axisym geometry
        call adjacent_front_edge(mesh,nb_arrete,front_edge,near_point,angle_maxA,angle_maxB,ierror)
    endif
    
    ! Ajout point la fin de la liste de points
    ptr_np => near_point
    do while (associated(ptr_np)) 
        ptr_np => ptr_np%suiv
    enddo
    Q%index = nb_point+1
    Q%coord = P_opt
            
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Tests on the new points
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Test Proximite (boucle sur les aretes du front)
    is_prox = .false.
    index_ed = 1 !ptr_ed => front_edge
    
    ! FIXME : mesure des normes selon la metrique associee
    do while (index_ed .le. nb_arrete .and. .not.is_prox) !do while (associated(ptr_ed) .and. .not.is_prox)
        
        !index_ed = ptr_ed%val
        index_1 = mesh%arrete(index_ed)%iP(1)
        index_2 = mesh%arrete(index_ed)%iP(2)
        P1 = mesh%point(index_1)%coord
        P2 = mesh%point(index_2)%coord
        
        if(imetric==2)then
            call test_prox_cone(P1(1:2),P2(1:2),Q%coord(1:2),param,delta_1,is_prox)
        elseif(imetric==3)then
            call test_prox_axisym(P1(1:2),P2(1:2),Q%coord(1:2),param,delta_1,is_prox)
        else
            !dist = sqrt((P2(1)-P1(1))**2+(P2(2)-P1(2))**2)
            dist = norm2(P2-P1)
            !dir2 = [P2(1)-P1(1),P2(2)-P1(2)]/dist         !dir2 = mesh%arrete(index_ed)%dir(1:2) 
            dir2 = (P2-P1)/dist
            !vect = [Q%coord(1)-P1(1),Q%coord(2)-P1(2),0._rp]
            vect = Q%coord-P1
            !dist2 = dot_product(vect(1:2),dir2)                  ! P1Q dot P1P2
            dist2 = dot_product(vect,dir2)
            !dist3 = dot_product(vect,[-dir2(2),dir2(1),0._rp])
            call Computation_vect_product(dir2,vect,vect_product_1)
            dist3  = norm2(vect_product_1)
            is_prox = dist2 .le. dist+Epsilon .and. dist2 .gt. -Epsilon .and. abs(dist3) .lt. 0.1*delta_1
        endif
        index_ed = index_ed+1 !ptr_ed => ptr_ed%suiv
        
    enddo
    
    ! Ajout P_opt if close enough
    if(.not.is_prox)then
        call ajout_fin_pile(near_point,Q,990._rp)
    endif
    
    ! This part is surely useless.
    do j=1,ndiv
        dd = (ndiv-j)*dist3 / ndiv
        Q%coord = [M(1) + dd*dir(1),M(2)+dd*dir(2),0._RP] ! dir is not initialized, nothing is done with Q.
    enddo
    
    ! Tests triangle admissible
    is_tri = .false.
    do while (.not.is_tri .and. associated(near_point)) ! Boucle sur les aretes du front
        
        is_angle = .false.
        is_inter = .false.
        is_inter1 = .false.
        is_inter2 = .false.
        
        ! Point to test
        Q = near_point%val
        
        ! Tests intersection for each node of the mesh front.
        ptr_ed => front_edge
        do while (associated(ptr_ed) .and. .not.is_inter)
            
            if (ptr_ed%val /= index_arr) then
                index_ed = ptr_ed%val
                index_1 = mesh%arrete(index_ed)%iP(1)
                index_2 = mesh%arrete(index_ed)%iP(2)
                P1 = mesh%point(index_1)%coord
                P2 = mesh%point(index_2)%coord
                
                ! Test of the intersection between QA and P1P2.
                if(imetric==3)then
                    if(.not.((index_A==index_1 .and. Q%index==index_2) .or. (index_A==index_2 .and. Q%index==index_1)))then
                        call segment_intersect_2D_2(Q%coord(1:2),A(1:2),P1(1:2),P2(1:2),is_inter1)
                    endif
                else
                    call segment_intersect_2D(Q%coord(1:2),A(1:2),P1(1:2),P2(1:2),is_inter1)
                endif
                
                ! Test of the intersection between QB and P1P2.
                if(imetric==3)then
                    if(.not.((index_B==index_1 .and. Q%index==index_2) .or. (index_A==index_2 .and. Q%index==index_1)))then
                        call segment_intersect_2D_2(Q%coord(1:2),B(1:2),P1(1:2),P2(1:2),is_inter2)
                    endif
                else
                    call segment_intersect_2D(Q%coord(1:2),B(1:2),P1(1:2),P2(1:2),is_inter2)
                endif
                
                is_inter = is_inter1 .or. is_inter2
                
            endif
            ptr_ed => ptr_ed%suiv !index_ed = index_ed+1
            
        end do
                
        ! Q is not one of the extremities of the edge.
        if( (abs(Q%coord(1)-mesh%point(index_A)%coord(1)).gt.Epsilon*10 .and. abs(Q%coord(2)-mesh%point(index_A)%coord(2)).gt.Epsilon*10) .or.  (abs(Q%coord(1)-mesh%point(index_B)%coord(1)).gt.Epsilon*10 .and. abs(Q%coord(2)-mesh%point(index_B)%coord(2)).gt.Epsilon*10) ) then ! PYW - Q != A or B
            
            if(.not.is_inter)then
                call internal_node2D(mesh,mesh%point(index_A),mesh%point(index_B),Q,front_vertex,is_intnode,ierror)
            endif
            
            if(.not. (is_inter .or. is_intnode) .and. imetric==3 .and. index_A == 103 .and. index_B == 147)then
                call limit_angle(mesh%point(index_A)%coord,mesh%point(index_B)%coord,Q%coord,angle_maxA,angle_maxB,is_angle,ierror)
            else
                is_angle = .false.
            endif
        end if
        
        if (is_inter .or. is_intnode .or. is_angle) then
            !call depile(near_point,ierror)
            !if(ierror/=0)then
            !    ierror = 3
            !    goto 9999
            !endif
            near_point => near_point%suiv ! PYW
        else ! if (is_inter .or. is_intnode .or. is_angle) = .false.
            is_tri = .true.
        endif
         
    enddo
    
    if (is_tri) then ! A new panel is created.
        new_point%coord = Q%coord
        new_point%index = Q%index        
    else ! A new panel cannot be created. Modification of delta_2 to try again.
        if (count.lt.5)then
            ! Wider search area of the neighboors of Popt.
            delta_2 = (2**(count+1))*delta_1
            count = count + 1
            if(idebug>=1) write(*,30),index_A,index_B      
            goto 100
        else
            ! Popt is searched closer to M (eq 4.12 of CC).
            if(not(PcloserM))then
                PcloserM = .true.
                count = 0
                goto 3
            elseif(PcloserM)then
                ierror = 5
                goto 9999
            end if
        endif
    endif
    
    ! Affichage coord. du point
    if (idebug > 1) then
        write(*,50)index_A,index_B,new_point%index,new_point%coord
    endif
    50 format("[",i5,",",i5,'] -> ',i5,' : ',3f8.2)
    30 format("** warning : boucle pour la recherche point optimal [A:",i3,",B:",i3,']')
    
    9999 continue
    if(ierror /= 0) then
        print*,""
        print*,"** error #",ierror,": optimal point research"
        print*,"** index_A = ",index_A
        print*,"** index_B = ",index_B
        print*,"** P   = ",P
        print*,"** index_pt= ",new_point%index
        print*,"** P_opt   = ",P_opt
        print*,"** delta_1 = ",delta_1
        print*,"** delta_2 = ",delta_2
        if(abs(delta_1).ge.Epsilon)then
            print*,"** delta_2/delta_1 = ",delta_2/delta_1
        end if
        print*,""
        ierr = ierror
    endif      
    
    select case(ierror)
    case(100)
        print*,'** i1 = ',i1,' i2 = ',i2,' j1 = ',j1,' j2 = ',j2
        print*,'** x1 = ',xgrid(i1),' x2 = ',xgrid(i2),' y1 = ',ygrid(j1),' y2 = ',ygrid(j2)
    case default
    end select
    
    nullify(ptr_ed,ptr_fv,near_point)
    
end subroutine optimalpoint  

subroutine updatefront(mesh,nb_point,nb_arrete,nb_tri,front_edge,front_vertex,new_point,iface,ierr)

    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(inout) :: mesh
    integer,intent(inout) :: nb_point,nb_arrete,nb_tri
    !f2py integer*1, dimension(1000) :: front_edge,front_vertex
    type(pile_element),pointer :: front_edge,front_vertex
    !f2py integer*1, dimension(1000) :: new_point
    type(MVertex),intent(inout) :: new_point
    integer,intent(in) :: iface
    integer,intent(inout) :: ierr
    ! local
    integer :: index_pt,index_A,index_B,index_arr
    integer :: iA
    integer :: nv_arrete,j
    !f2py integer*1, dimension(1000) :: edge_1,edge_2
    type(MEdge) :: edge_1,edge_2
    !f2py integer*1, dimension(1000) :: tri
    type(MElement) :: tri
    logical :: is_pres,res
    integer,dimension(2) :: ipt
  
    ! This subroutine updates the mesh front.
  
    index_pt = new_point%index
    index_arr = front_edge%val
    call depile(front_edge,ierr)
    index_A = mesh%arrete(index_arr)%iP(1)
    index_B = mesh%arrete(index_arr)%iP(2)
    ipt = -9
    if(index_pt > nb_point) then
        ! Ajout nouveau point au maillage
        new_point%nface=1
        new_point%face(1)=iface
        call add_element(mesh%point,nb_point,new_point)
        call add_element_pile(front_vertex,index_pt,1._rp,0)
        
        call edge_from_vertex(mesh%point(index_A),mesh%point(index_pt),edge_1)
        call edge_from_vertex(mesh%point(index_pt),mesh%point(index_B),edge_2)
        
        edge_1%bf = 0
        edge_1%nface=1
        edge_1%face(1)=iface
        call add_element(mesh,nb_arrete,edge_1)
        edge_1%index = nb_arrete
        edge_2%bf = 0
        edge_2%nface = 1
        edge_2%face(1) = iface
        call add_element(mesh,nb_arrete,edge_2)
        edge_2%index = nb_arrete
        
        call add_element_pile(front_edge,edge_1%index,edge_1%length,0)
        call add_element_pile(front_edge,edge_2%index,edge_2%length,0) 
    else
        ! Utilisation d'un point appartenant au front
        nv_arrete = mesh%point(index_pt)%nv_arrete
        !
        call edge_from_vertex(mesh%point(index_A),mesh%point(index_pt),edge_1)
        is_pres = .false.
        j = 1
        ! Recherche si l'arrete cr\E9\E9e appartient dej\E0 au front
        do while(.not.is_pres .and. j<=nv_arrete)
            iA = mesh%point(index_pt)%va(j)
            call ineg_edge(edge_1, mesh%arrete(iA),res)
            is_pres = .not.res
            j = j+1
        enddo
        if(.not.is_pres) then
            ! Creation nouvelle arrete avec un point du front
            edge_1%bf = 0
            edge_1%nface = 1
            edge_1%face(1) = iface
            call add_element(mesh,nb_arrete,edge_1)
            edge_1%index = nb_arrete
            call add_element_pile(front_edge,edge_1%index,edge_1%length,0)
        else
            call remove_element_pile(front_edge,iA,ierr)
        endif
        
        nv_arrete = mesh%point(index_pt)%nv_arrete
        call edge_from_vertex(mesh%point(index_pt),mesh%point(index_B),edge_2)
        is_pres = .false.
        ipt = -9
        j = 1
        ! Recherche si l'arrete appartient dej\E0 au front
        do while(.not.is_pres .and. j<=nv_arrete)
            iA = mesh%point(index_pt)%va(j)
            call ineg_edge(edge_2,mesh%arrete(iA),res)
            is_pres = .not.res
            j = j+1
        enddo
        if(.not.is_pres) then
            ! Creation nouvelle arrete avec un point appartenant au front
            edge_2%bf = 0
            edge_2%nface = 1
            edge_2%face(1) = iface
            call add_element(mesh,nb_arrete,edge_2)
            edge_2%index = nb_arrete
            call add_element_pile(front_edge,edge_2%index,edge_2%length,0)
        else
            call remove_element_pile(front_edge,iA,ierr)
        endif
    end if
    
    call triangle_from_edge(edge_1,edge_2,mesh%arrete(index_arr),tri)
    tri%face = iface
    call add_element(mesh%tri,nb_tri,tri)
    tri%index = nb_tri
    
    ! Update front edge
    
    ! Point A
    nv_arrete = mesh%point(index_A)%nv_arrete
    is_pres = .false.
    j = 1
    do while(j<=nv_arrete .and. .not.is_pres)
        call element_in_pile(front_edge,mesh%point(index_A)%va(j),is_pres)
        j=j+1
    enddo
    if(.not.is_pres)then
        call remove_element_pile(front_vertex,index_A,ierr)
    endif
    
    ! Point B
    nv_arrete = mesh%point(index_B)%nv_arrete
    is_pres = .false.
    j = 1
    do while(j<=nv_arrete .and. .not.is_pres)
        call element_in_pile(front_edge,mesh%point(index_B)%va(j),is_pres)
        j=j+1
    enddo
    if(.not.is_pres)then
        call remove_element_pile(front_vertex,index_B,ierr)
    endif
    
    ! Point Optimal
    nv_arrete = mesh%point(index_pt)%nv_arrete
    is_pres = .false.
    j = 1
    do while(j<=nv_arrete .and. .not.is_pres)
        call element_in_pile(front_edge,mesh%point(index_pt)%va(j),is_pres)
        j=j+1
    enddo
    if(.not.is_pres)then
        call remove_element_pile(front_vertex,index_pt,ierr)
    endif  
    
end subroutine updatefront

subroutine updatemesh(mesh1,nb_point1,nb_arrete1,nb_tri1,mesh2,nb_point2,&
&                     nb_arrete2,nb_tri2,ind_pt,ind_arr,typeFrontiere,iface,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh1
  type(MGrid),intent(inout) :: mesh1
  integer,intent(inout) :: nb_point1,nb_arrete1,nb_tri1
  !f2py integer*1, dimension(1000) :: mesh2
  type(MGrid),intent(in) :: mesh2
  integer,intent(in) :: nb_point2,nb_arrete2,nb_tri2
  integer,dimension(*),intent(inout) :: ind_pt,ind_arr
  integer,intent(in) :: typeFrontiere
  integer,intent(in) :: iface
  integer,intent(inout) :: ierror
! local
  integer :: j,ind
  integer :: ip1,ip2,ip3
  !f2py integer*1, dimension(1000) :: vertex
  type(MVertex) :: vertex
  !f2py integer*1, dimension(1000) :: edge
  type(MEdge)   :: edge
  !f2py integer*1, dimension(1000) :: tri
  type(MElement) :: tri
  !f2py integer*1, dimension(1000) :: V
  type(vector) :: V
  
  do j=1,nb_point2
    vertex = mesh2%point(j)
    vertex%nv_arrete = 0
    vertex%va = -9
    vertex%nv_point = 0
    vertex%vp = -9
    vertex%bf = mesh2%point(j)%bf
    if(vertex%bf .le. 0 .and. vertex%bf.ge.-8 .or. vertex%bf==-11)then
      if(vertex%bf==-11) vertex%bf = 0
      !vertex%nface = 1
      !vertex%face(1) = iface
      call add_element(mesh1%point,nb_point1,vertex)
      ind_pt(j) = nb_point1
    elseif(vertex%bf == -12)then ! actualisation correspondant point virtuel avec mesh1
      ind = ind_pt(j)
      ind_pt(j) = ind_pt(ind)    
    elseif(vertex%bf == -10) then
      ind = ind_pt(j)
      mesh1%point(ind)%bf = 0 
    endif
  enddo    
  
  do j=1,nb_arrete2
    edge = mesh2%arrete(j)
    if(edge%bf.le.0 .and. edge%bf.ge.-8 .or. edge%bf==-11)then
      ip1 = ind_pt(edge%iP(1))
      ip2 = ind_pt(edge%iP(2))
      edge%iP(1) = ip1
      edge%iP(2) = ip2
      edge%P1 = mesh1%point(ip1)%coord
      edge%P2 = mesh1%point(ip2)%coord
      !V = edge%P2 - edge%P1
      call assign_vector_coord(V,edge%P2 - edge%P1)
      edge%dir = V%coord/V%length
      edge%length = V%length
      if (edge%bf == -8 .or. edge%bf==-11) edge%bf = 0
      !edge%nface = 1
      !edge%face(1) = iface
      call add_element(mesh1,nb_arrete1,edge)
!     Actualisation tableau correspondance \E0 rajouter si utilis\E9   
    elseif(edge%bf==-12) then
      ind = ind_arr(j)
      ind_arr(j) = ind_arr(ind)
    elseif(edge%bf == -10) then
      ind = ind_arr(j)
      mesh1%arrete(ind)%bf = 0
    endif
  enddo
  
  do j=1,nb_tri2
    tri = mesh2%tri(j) 
    ip1 = ind_pt(tri%iP(1))
    ip2 = ind_pt(tri%iP(2))
    ip3 = ind_pt(tri%iP(3))
    tri%iP(1) = ip1
    tri%iP(2) = ip2
    tri%iP(3) = ip3
    tri%P1 = mesh1%point(ip1)%coord
    tri%P2 = mesh1%point(ip2)%coord
    tri%P3 = mesh1%point(ip3)%coord
    tri%typeFrontiere = typeFrontiere
    !tri%face = iface
    call add_element(mesh1%tri,nb_tri1,tri)
  enddo
    
end subroutine updatemesh
  
subroutine int_bilin_cart(P,val,xmin,xmax,ymin,ymax,val11,val21,val12,val22,ierror)
        
    real(rp),dimension(2),intent(in)    :: P                        ! Point where the value is computed
    real(rp),intent(inout)              :: val                      ! Output value
    real(rp),intent(in)                 :: xmin,xmax,ymin,ymax      ! Boundaries of the cell of the Cartesian grid
    real(rp),intent(in)                 :: val11,val12,val21,val22  ! Values on the boundaries of the cell of the Cartesian grid
    integer,intent(inout)               :: ierror                   ! Error flag
    
    real(rp)                            :: alpha,beta               ! Parameters for the bilinear interpolation
    
    ! This subroutine does a bilinear interpolation on a cell of the Cartesian grid.
    
    ierror = 0
  
    if(abs(xmax-xmin).gt.Epsilon)then
        alpha = (P(1)-xmin)/(xmax-xmin)
    else
        alpha = 0._rp
    endif
  
    if(abs(ymax-ymin).gt.Epsilon)then
        beta = (P(2)-ymin)/(ymax-ymin)
    else
        beta = 0._rp
    endif
  
    val = (1._rp-alpha)*(1._rp-beta)*val11 + alpha*(1._rp-beta)*val21 +&
    &     alpha*beta*val22 + (1._rp-alpha)*beta*val12
  
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error #',i3,' : interpolation bilineaire impossible') 

    select case(ierror)
    case(100)
        print*,'** xmax = ',xmax,' xmin = ',xmin
    case(200)
        print*,'** ymax = ',ymax,' ymin = ',ymin
    case default
    end select
  
end subroutine int_bilin_cart

subroutine write_tec(iunit,mesh,t,nb_point,nb_tri,ierror)
    
    integer,intent(in)                  :: iunit                ! Output reference
    !f2py integer*1, dimension(1000)    :: mesh                 
    type(MGrid),intent(in)              :: mesh                 ! MGrid
    real(rp),intent(in)                 :: t                    ! Current time
    integer,intent(in)                  :: nb_point,nb_tri      ! Number of points and triangles of the surface
    integer,intent(inout)               :: ierror               ! Error flag
      
    integer                             :: j                    ! Loop parameter
  
    ! This subroutine writes the mesh during the advance front.
  
    if(nb_tri.gt.0)then
        write(iunit,10),t,nb_point,nb_tri,t
10      format('Zone T= "',f16.1,' seconds", N =',i16,', E = ',i16,&
        &         ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME =',f16.1)
    end if
    
    do j =1,nb_point
        write(iunit,'(3f24.5,i16,i5,i5)'),mesh%point(j)%coord,j,mesh%point(j)%bf,mesh%point(j)%nface
    enddo
    do j=1,nb_tri
        write(iunit,'(3i16)'),mesh%tri(j)%iP
    enddo
  
    9999 continue
    if(ierror.ne.0)then
        print*,'** error #',ierror,' : ecriture du fichier .tec impossible'  
    endif
  
end subroutine write_tec

subroutine write_Cartesian_grid(xgrid,ygrid,nx,ny,matdref,t,index)

    integer,intent(in)                      :: nx,ny    ! Size of xgrid and ygrid.
    real(rp),dimension(nx),intent(in)       :: xgrid    ! x-positions of the Cartesian grid.
    real(rp),dimension(ny),intent(in)       :: ygrid    ! y-positions of the Cartesian grid.
    real(rp),dimension(nx,ny,3),intent(in)  :: matdref  ! Reference length for each cell of the Cartesian grid.
    real(rp),intent(in)                     :: t        ! Current time.
    integer,intent(in)                      :: index    ! Index to identify the part of the mesh.
    
    integer                                 :: j,p,q    ! Loop parameters.
    character(len=20)                       :: num      ! Solution time.
    integer                                 :: ierror   ! Error flag.
    character(len=200)                      :: fileGrid ! Temporary file name.
      
    ! This subroutine writes the Cartesian grid into a Tecplot file.
        
    if(iwCartGrid)then
        
        ! Current time.
        write(num, '( f0.5 )') t
        
        ! File name.
        write(fileGrid,'("Cartesian_grid_",i2,".dat")'),index
        
        ! Opening.
        open(ioCartGrid,file=fileGrid, iostat=ierror)
        write(ioCartGrid,fmt='(50a)') 'Title = "Cartesian grid"'
        write(ioCartGrid,fmt='(150a)') 'VARIABLES = "X","Y","Number","MatDref"'
        write(ioCartGrid,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', nx*ny, ', E=', nx*(ny-1)+ny*(nx-1),& ! There are nx*ny points and nx*(ny-1)+ny*(nx-1) segments.
        &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
        
        ! Filling points.
        j = 1
        do p = 1,nx
            do q = 1,ny
                write(ioCartGrid,'(2E,I)'),xgrid(p),ygrid(q),j,matdref(p,q,1)
                j = j + 1
            end do
        end do
        
        ! Filling vertical segments.
        do p = 1,nx
            do q = 1,ny-1
                write(ioCartGrid,'(4I)'),q + (p-1)*ny,q + 1 + (p-1)*ny,q + 1 + (p-1)*ny ! Two times because they are triangles
            end do
        end do
        
        ! Filling horizontal segments.
        do q = 1,ny
            do p = 1,nx-1
                write(ioCartGrid,'(4I)'),q + (p-1)*ny,q + nx + (p-1)*ny,q + nx + (p-1)*ny ! Two times because they are triangles
            end do
        end do
        
        ! Closing.
        close(ioCartGrid)
        
    end if
    
end subroutine write_Cartesian_grid

subroutine write_debug(mesh,nb_point,nb_arrete,nb_tri,t,&
                     & filename1,filename2,filename3,ierr,imetric,param)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                             ! Mgrid
    integer,intent(in)                  :: nb_point,nb_arrete,nb_tri        ! Number of points, edges and triangles.
    real(rp),intent(in)                 :: t                                ! Current time
    character(*)                        :: filename1,filename2,filename3    ! Output file names
    integer,optional                    :: imetric                          !
    real(rp),dimension(*),optional      :: param                            ! Parameters
    integer,intent(inout)               :: ierr                             ! Error flag
    
    integer                             :: k                                ! Loop parameter
    real(rp),dimension(3)               :: M                                ! Vector
    
    ! This subroutine writes the output files point.dat, arrete.dat and mesh.dat from the MGrid structure.
    
    if(present(imetric))then
        if(imetric>0)then
            do k=1,nb_point
                M = mesh%point(k)%coord
                M(3) = fz(M(1:2),imetric,param)
                mesh%point(k)%coord = M
            enddo
        endif
    endif
    
    ! Point.dat
    print*,'** Ecriture fichiers : ',filename2,' ...'
    open(unit=40,file=filename2)
    write(40,*) 'Title = "Noeud du maillage"'
    write(40,'(100a)') 'VARIABLES = "X","Y","Z","bf","nface","face1","face2","face3","nedge","edge1","edge2","edge3"'
    do k=1,nb_point
        write(40,*),mesh%point(k)%coord,mesh%point(k)%bf,mesh%point(k)%nface,&
        & mesh%point(k)%face(1:3),mesh%point(k)%nedge,mesh%point(k)%edge(1:3)
    enddo
    close(40)
    
    ! Arrete.dat
    print*,'**                   : ', filename3,' ...'
    open(unit=50,file=filename3)
    write(50,*) 'Title = "Table des aretes"'
    write(50,*)'VARIABLES = "ind","P1","P2","bf","nface","face1","face2"'
    do k=1,nb_arrete
        write(50,*),k,mesh%arrete(k)%iP(1:2),mesh%arrete(k)%bf,&
        & mesh%arrete(k)%nface,mesh%arrete(k)%face(1),mesh%arrete(k)%face(2)
    end do
    10 format(' ind = ',i5,' : [',i6,i6,']',i4,3i5)   
    close(50) 
    
    ! Mesh.dat
    print*,' **                 : ',filename1,' ...'
    open(unit=60,file=filename1)   
    write(60,*),'TITLE = "Mesh2D"'
    write(60,*),'VARIABLES = "X","Y","Z","LABEL","bf","nface"'
    call write_tec(60,mesh,t,nb_point,nb_tri,ierr)
    close(60)
   
 end subroutine write_debug
 
subroutine local_coord_cyl(mesh,nb_point,nb_arete,nb_tri,fgeom,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  integer,intent(in) :: nb_point,nb_arete,nb_tri
  !f2py integer*1, dimension(1000) :: fgeom
  type(cylindre_2),intent(in) :: fgeom
  integer,intent(inout) :: ierror
! local
  real(rp),dimension(2) :: y
  !f2py integer*1, dimension(1000) :: P
  type(point) :: P
  integer :: j
!    
  if(idebug/=0) print*,"MP_local_coord_cyl"

  ierror = 0
!  
  do j=1,nb_point
    P = mesh%point(j)%coord
    call cart2cyl(P,fgeom,y)
    mesh%point(j)%coord = [fgeom%rayon*y(1),y(2),0._rp]
  end do
!  
end subroutine local_coord_cyl

subroutine local_coord_cone(mesh,nb_point,nb_arete,nb_tri,fgeom,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  integer,intent(in) :: nb_point,nb_arete,nb_tri
  !f2py integer*1, dimension(1000) :: fgeom
  type(cone),intent(in) :: fgeom
  integer,intent(inout) :: ierror
! local
  real(rp),dimension(2) :: y
  real(rp),dimension(2) :: M
  !f2py integer*1, dimension(1000) :: P
  type(point) :: P
  integer :: j
  real(rp) :: r,a,b,z
 
  ierror = 0
  a = (fgeom%r2 - fgeom%r1)/fgeom%long
  b = (fgeom%r1*fgeom%vmax - fgeom%r2*fgeom%vmin)/fgeom%long

  do j=1,nb_point
    P = mesh%point(j)%coord
    call cart2cone(P,fgeom,y)
    r = a*(y(2)-fgeom%vmin) + fgeom%r1
    M(1:2) = [r*y(1),y(2)]
    !y(1) = y(1)-PI
    call fz_cone(M,a,b,z)
    mesh%point(j)%coord = [r*y(1),y(2),z]
  enddo
!
end subroutine local_coord_cone

subroutine local_coord_rep(mesh,nb_point,nb_arete,nb_tri,rep,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  integer,intent(in) :: nb_point,nb_arete,nb_tri
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d),intent(in) :: rep
  integer,intent(inout) :: ierror
! local
  integer :: j
  real(rp),dimension(3) :: M
!  
  ierror = 0
!  
  do j=1,nb_point
    M = mesh%point(j)%coord
    call cart2loc(M,rep,mesh%point(j)%coord)
  enddo
! 
end subroutine local_coord_rep

subroutine orient_hole_mesh(mesh,nb_point,nb_arete,np0,na0,ierror)
    
    !f2py integer*1, dimension(1000)    :: mesh                 
    type(MGrid),intent(inout)           :: mesh                 ! MGrid.
    integer,intent(in)                  :: nb_point,nb_arete    ! Number of points and edges after mesh_inter.
    integer,intent(in)                  :: np0,na0              ! Number of points and edges in mesh0D.
    integer,intent(inout)               :: ierror               ! Error flag.
    
    integer                             :: j                    ! Loop parameter.
    real(rp)                            :: ps                   ! Dot product result.
    real(rp),dimension(2)               :: yg,ym,y,yn,v         ! 2D vectors.
    integer                             :: i1,i2                ! Indexes.
    real(rp),dimension(2)               :: P1,P2                ! Points.
    real(rp)                            :: xmin,xmax,ymin,ymax  
    
    ! This subroutine orientates the edges of the intersection curve for the advance front method.
    
    ierror = 0

    ! Calcul coordonnees barycentriques (centre of gravity of the projection of the intersection curve on the plane (x,y))
    if(.false.)then
        yg = [0._rp,0._rp]
        do j=np0+1,nb_point
            y = mesh%point(j)%coord(1:2)
            yg = yg + y
        enddo
        yg = yg/dble(nb_point-np0)
    else
        xmin = 999. ; xmax = -999.
        ymin = 999. ; ymax = -999.
        do j=np0+1,nb_point            
            y = mesh%point(j)%coord(1:2)                        
            if (y(1).lt.xmin)then
                xmin = y(1)
            endif
            if(y(1).gt.xmax)then
                xmax = y(1)
            endif
            if(y(2).lt.ymin)then
                ymin = y(2)
            endif
            if(y(2).gt.ymax)then
                ymax = y(2)
            endif
        enddo
        yg(1) = 0.5_rp*(xmin+xmax)
        yg(2) = 0.5_rp*(ymin+ymax)
    endif
    
    do j=na0+1,nb_arete
        i1 = mesh%arrete(j)%iP(1)
        i2 = mesh%arrete(j)%iP(2)
        P1 = mesh%point(i1)%coord(1:2)
        P2 = mesh%point(i2)%coord(1:2)
                
        ! Normal to the edge
        yn = [P1(2)-P2(2),P2(1)-P1(1)]
        yn = yn/norm2(yn) ! Unit vector
        
        ! Coordinate of the centre of the edge : M
        ym = 0.5*(P1+P2)
        
        ! Vector GM
        v = ym-yg 
        
        ! GM.n
        ps = dot_product(v,yn) ! GM.n
                
        ! Orientation de l'arete    
        if(abs(ps).lt.Epsilon)then
            ierror = 1
            goto 9999
        elseif(ps > 0._rp)then ! Interversion of the end points of the edge.
            mesh%arrete(j)%iP(1) = mesh%arrete(j)%iP(2)
            mesh%arrete(j)%iP(2) = i1
        endif
    enddo   
  
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error #',i3,' : orient_hole_mesh: Orientation hole impossible')  
  
end subroutine orient_hole_mesh

subroutine orient_hole_mesh2(mesh,nb_point,nb_arete,np0,na0,ierror)
  
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                     ! MGrid.
    integer,intent(in)                  :: nb_point,nb_arete        ! Number of points and edges after mesh_inter.
    integer,intent(in)                  :: np0,na0                  ! Number of points and edges in mesh0D.
    integer,intent(inout)               :: ierror                   ! Error flag.
    
    integer                             :: j                        ! Loop parameter.
    real(rp)                            :: ps                       ! Dot product result.
    real(rp),dimension(2)               :: yg,ym,y,yn,v
    integer                             :: i1,i2                    ! Indexes
    real(rp),dimension(2)               :: P1,P2                    ! Points.
    real(rp)                            :: xmin,xmax,ymin,ymax
    integer                             :: ncom1,ncom2,tab(nfvmax)
  
    ! This subroutine orientates the edges of the intersection curve for the advance front method.
  
    ierror = 0

    ! Calcul coordonnees barycentriques (centre of gravity of the projection of the intersection curve on the plane (x,y))
    if(.false.)then
        yg = [0._rp,0._rp]
        do j=np0+1,nb_point
            y = mesh%point(j)%coord(1:2)
            yg = yg + y
        enddo
        yg = yg/dble(nb_point-np0)
    else
        xmin = 999. ; xmax = -999.
        ymin = 999. ; ymax = -999.
        do j=np0+1,nb_point
            y = mesh%point(j)%coord(1:2)        
            if (y(1).lt.xmin)then
                xmin = y(1)
            endif
            if(y(1).gt.xmax)then
                xmax = y(1)
            endif
            if(y(2).lt.ymin)then
                ymin = y(2)
            endif
            if(y(2).gt.ymax)then
                ymax = y(2)
            endif
        enddo
        yg(1) = 0.5_rp*(xmin+xmax)
        yg(2) = 0.5_rp*(ymin+ymax)
    endif
    
    do j=na0+1,nb_arete
        i1 = mesh%arrete(j)%iP(1)
        i2 = mesh%arrete(j)%iP(2)
        P1 = mesh%point(i1)%coord(1:2)
        P2 = mesh%point(i2)%coord(1:2)

        !call common_int(mesh%point(i1)%face,mesh%point(i1)%nface,[1],1,tab,ncom1)
        !call common_int(mesh%point(i2)%face,mesh%point(i2)%nface,[1],1,tab,ncom2)
        !is_FS = ncom1==1 .and. ncom2 ==1
        
        ! Normal to the edge
        yn = [P1(2)-P2(2),P2(1)-P1(1)]    
        yn = yn/norm2(yn) ! Unit vector
        
        ! Coordinate of the centre of the edge : M   
        ym = 0.5*(P1+P2)
        
        ! vector GM
        v = ym-yg
        
        ! GM.n
        ps = dot_product(v,yn)
        
        !   Orienation de l'arete
        if(abs(ps).lt.Epsilon)then
            ierror = 1
            goto 9999
        !elseif(ncom1==1 .and. ncom2==1 .and. P2(1).gt.P1(1))then
        !  mesh%arrete(j)%iP(1) = mesh%arrete(j)%iP(2)
        !  mesh%arrete(j)%iP(2) = i1
        !elseif((ncom1/=1 .or. ncom2/=1 ).and. ps > 0._rp)then
        elseif(ps>0._rp)then
            mesh%arrete(j)%iP(1) = mesh%arrete(j)%iP(2)
            mesh%arrete(j)%iP(2) = i1
        endif
    enddo   
  
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error #',i3,' : orient_hole_mesh2: Orientation hole impossible')  
  
end subroutine orient_hole_mesh2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Change type Mesh
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine change_type_mesh(mesh,nb_point,nb_tri,t,Maillage,nface,fgeom_vect,ierror,InputData)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                         ! MGrid
    integer,intent(in)                  :: nb_point,nb_tri              ! Number of points and triangles
    real(rp),intent(in)                 :: t                            ! Current time
    !f2py integer*1, dimension(1000)    :: maillage
    type(TMaillage),intent(inout)       :: maillage                     ! Mesh
    integer,intent(in)                  :: nface                        ! Number of faces
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                   ! Geometries of the bodies
    integer,intent(inout)               :: ierror                       ! Error flag
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                    ! Input data
    
    integer                             :: j, k, n, iface, Nnoeud,jj    !
    integer                             :: id_sl, np_sl, itemp, id_nc   !
    integer                             :: index_sl, ind, nc            !
    integer                             :: i1, i2, i3                   !
    integer                             :: typeFrontiere                !
    integer                             :: NCorps                       ! NBodies + 1 (tank)
    integer,dimension(nface)            :: ICorps                       !
    integer,dimension(:),allocatable    :: id_noeud, id_tri             !
    integer,dimension(:,:),allocatable  :: table_corres                 !
    integer,dimension(3)                :: iP                           !
    real(rp)                            :: L, LL, prof                  !
    !f2py integer*1, dimension(1000)    :: TNtemp
    type(TNoeud)                        :: TNtemp                       !
    !f2py integer*1, dimension(1000)    :: TFtemp
    type(TFacette)                      :: TFtemp                       !
    !f2py integer*1, dimension(1000)    :: vertex
    type(MVertex)                       :: vertex                       !
    logical                             :: bool                         ! Boolean used because fgeom_vect%Active is not defined for the tank.
        
    ! This subroutine creates the structure TMaillage from MGrid.
    
    ierror = 0
    index_sl = HouleRF%index    
    
    ! NCorps
    NCorps = NBodies + 1 ! + 1 = Tank
        
    ! Dimension de la cuve

    if(idtype==1)then
        L = Ldom(1)
        LL = Ldom(2)
    elseif(idtype==2)then
        L = 2.*Ldom(5)
        LL = L
    endif
    prof = Ldom(3)
    maillage%DimTank(1:4) = [L, LL, prof, Labs]
    
    Maillage%NBody = NCorps
    
    ! Definition of ICorps
    ICorps = 0._RP
    jj = 1
    do j = 1,NBodies+1 ! +1 for the tank.
        if(j.gt.1)then
            bool = fgeom_vect%Active(j-1)
        else
            bool = .false.
        end if
        if(j.eq.1 .or. bool)then! Active geometry
            do k = 1,fgeom_vect%nface_vect(j) ! nface = Somme des nface_vect
                if(j == 1)then ! Free surface and tank
                    if (k == 1)then
                        ICorps(jj) = 0 ! Free surface
                    else
                        ICorps(jj) = 1 ! Tank
                    end if
                else ! Bodies
                    ICorps(jj) = j
                end if
                jj = jj + 1
            end do
        else
            do k = 1,fgeom_vect%nface_vect(j) ! nface = Somme des nface_vect
                ICorps(jj) = -1 ! No body
                jj = jj + 1
            end do
        end if
    end do
    
    allocate(table_corres(PointMax,nface))

    ! Change face number for non symetrical case
    do jj = 1,NBodies
        if(not(Symmetry) .and. InputData%igtype(jj).eq.5)then
            if(fgeom_vect%Active(jj))then
                if (.true.) call merge_face(Mesh,nb_point,nb_tri,fgeom_vect%geom(jj)%axisym(1)%index,fgeom_vect%geom(jj)%axisym(2)%index,ierror)
                if(ierror/=0) goto 9999
            end if
        endif
    end do
    
    ! ---------------------------------------------------------------------------------
    ! Operation sur les noeuds
    ! ---------------------------------------------------------------------------------
    
    Nnoeud = 0
    np_sl = 0
    
    ! Remplissage -------------------------------------------------------------------
    do j=1,nb_point
        vertex = mesh%point(j)
        n = vertex%nface
        do k=1,n
            iface = vertex%face(k)
            if(bottom_sym .and. iface==  2) goto 1000
            if(Symmetry   .and. iface== -1) goto 1000
            if(iface.ne.-1)then ! Probably iface can be equal to -1 if the definition of the structure arete: fgeom%arete(1)%face = [fgeom%cylindre(1)%index,-1]
                Nnoeud = Nnoeud+1
                if(Nnoeud.gt.Pointmax)then
                    print*,"change_type_mesh: The maximum number of nodes is reached. Please reduce it!"
                    pause
                end if
                maillage%Tnoeud(Nnoeud)%Pnoeud(1:3) = vertex%coord(1:3)
                maillage%Tnoeud(Nnoeud)%indmesh = j
                table_corres(j,iface) = Nnoeud
                if(iface==index_sl)then
                    maillage%Tnoeud(Nnoeud)%typeNoeud = 0
                    maillage%Tnoeud(Nnoeud)%Npanneau = 0
                    id_sl = Nnoeud
                    np_sl = np_sl+1
                else
                    ind = ICorps(iface)
                    maillage%Tnoeud(Nnoeud)%typeNoeud = 1
                    maillage%Tnoeud(Nnoeud)%Npanneau = ind
                endif
            end if
1000 continue
        enddo
    enddo
    maillage%Nnoeud = Nnoeud
    
    ! Tri
    allocate(id_noeud(Nnoeud), id_tri(Nnoeud))
    do j=1,Nnoeud
        id_noeud(j) = j
    enddo
    
    nc = 0
    j = 1
    maillage%FS%IndFS(1) = j ! Free surface
    do while(j .lt. id_sl)
        if(maillage%Tnoeud(j)%Npanneau > nc)then
            TNtemp = maillage%Tnoeud(j)
            maillage%Tnoeud(j) = maillage%Tnoeud(id_sl)
            maillage%Tnoeud(id_sl) = TNtemp
            itemp = id_noeud(j)
            id_noeud(j) = id_noeud(id_sl)
            id_noeud(id_sl) = itemp
            do while(maillage%Tnoeud(id_sl)%Npanneau .ne. nc)
                id_sl = id_sl-1
            enddo
        endif
        j = j+1
    enddo 
    maillage%FS%IndFS(3) = id_sl
    
    j = id_sl+1
    do nc = 1,NCorps ! Bodies and tank
        if(nc.gt.1)then
            bool = fgeom_vect%Active(nc-1)
        else
            bool = .false.
        end if
        if(nc.eq.1 .or. bool)then
            maillage%Body(nc)%IndBody(1) = j
            id_nc = Nnoeud
            do while(maillage%Tnoeud(id_nc)%Npanneau .ne. nc)
                id_nc = id_nc-1
            enddo
            do while(j .lt. id_nc)
                if(maillage%Tnoeud(j)%Npanneau > nc)then
                    TNtemp = maillage%Tnoeud(j)
                    maillage%Tnoeud(j) = maillage%Tnoeud(id_nc)
                    maillage%Tnoeud(id_nc) = TNtemp
                    itemp = id_noeud(j)
                    id_noeud(j) = id_noeud(id_nc)
                    id_noeud(id_nc) = itemp
                    do while(maillage%Tnoeud(id_nc)%Npanneau .ne. nc)
                        id_nc = id_nc-1
                    enddo
                endif
                j = j +1
            enddo
            maillage%Body(nc)%IndBody(3) = id_nc                             
            j = id_nc+1
        end if
    enddo
    
    ! -------------------------------------------------------------------------------------------
    ! Inversion du tableau id_noeud
    ! -------------------------------------------------------------------------------------------
    do j = 1,Nnoeud
        id_tri(id_noeud(j)) = j
    enddo    
    
    ! --------------------------------------------------------------------------------------------
    ! Operation sur les facettes
    ! --------------------------------------------------------------------------------------------    
    ! Remplissage
    do j=1,nb_tri
        typeFrontiere = mesh%tri(j)%typeFrontiere
        iface = mesh%tri(j)%face
        iP(1:3) = mesh%tri(j)%iP(1:3)
        i1 = table_corres(iP(1),iface)
        i2 = table_corres(iP(2),iface)
        i3 = table_corres(iP(3),iface)
        maillage%Tfacette(j)%Tnoeud(1:3) = [id_tri(i1),id_tri(i2),id_tri(i3)]
        maillage%Tfacette(j)%typeFrontiere = typeFrontiere
        if(iface == index_sl)then
            maillage%Tfacette(j)%Npanneau = 0
            maillage%Tfacette(j)%typeFrontiere = 0
            id_sl = j
        else
            ind = ICorps(iface) 
            maillage%Tfacette(j)%Npanneau = ind
            maillage%Tfacette(j)%typeFrontiere = 1
        endif
    enddo
    maillage%Nfacette = nb_tri
    
    ! Tri
    nc = 0
    j = 1
    maillage%FS%IndFS(2) = j ! Surface libre
    do while(j .lt. id_sl)
        if(maillage%Tfacette(j)%Npanneau > nc)then
            TFtemp = maillage%Tfacette(j)
            maillage%Tfacette(j) = maillage%Tfacette(id_sl)
            maillage%Tfacette(id_sl) = TFtemp
            do while(maillage%Tfacette(id_sl)%Npanneau .ne. nc)
                id_sl = id_sl-1
            enddo
        endif
    j = j+1
    enddo 
    maillage%FS%IndFS(4) = j
    j = j+1
     
    print*,""
    print*,maillage%FS%IndFS(4)
    print*,maillage%TFacette(maillage%FS%IndFS(4))%NPanneau
    print*,""
    if(maillage%TFacette(maillage%FS%IndFS(4))%NPanneau.eq.1)then
         pause
    end if
    
    do nc = 1,NCorps ! Corps
        if(nc.gt.1)then
            bool = fgeom_vect%Active(nc-1)
        else
            bool = .false.
        end if        
        if(nc.eq.1 .or. bool)then
            maillage%Body(nc)%IndBody(2) = j                            
            id_nc = nb_tri
            do while(maillage%Tfacette(id_nc)%Npanneau .ne. nc)
                id_nc = id_nc-1
            enddo
            do while(j .lt. id_nc)
                if(maillage%Tfacette(j)%Npanneau > nc)then
                    TFtemp = maillage%Tfacette(j)
                    maillage%Tfacette(j) = maillage%Tfacette(id_nc)
                    maillage%Tfacette(id_nc) = TFtemp
                    do while(maillage%Tfacette(id_nc)%Npanneau .ne. nc)
                        id_nc = id_nc-1
                    enddo
                endif
                j = j +1
            enddo
            maillage%Body(nc)%IndBody(4) = j                             
            j = j+1
        end if
    enddo       
    
    ! -----------------------------------------------------------------------------------------------
    ! Mise a jour indice tableau
    ! --------------------------------------------------------------------------------------------
    
    if(NCorps.gt.0)then
        maillage%Nsys = maillage%FS%IndFS(3)
        maillage%Nfsys = maillage%FS%IndFS(4)
        
        do nc = 1,NCorps ! Corps
            if(nc.gt.1)then
                bool = fgeom_vect%Active(nc-1)
            else
                bool = .true.
            end if
            if(nc.eq.1 .or. bool)then
                maillage%Nsys = maillage%Nsys + maillage%Body(nc)%IndBody(3) - maillage%Body(nc)%IndBody(1) + 1
                maillage%Nfsys = maillage%Nfsys + maillage%Body(nc)%IndBody(4) - maillage%Body(nc)%IndBody(2) + 1
            end if
        end do
        !maillage%Nsys = maillage%Body(NCorps)%IndBody(3)
        !maillage%Nfsys = maillage%Body(NCorps)%IndBody(4)
    else
        maillage%Nsys = maillage%FS%IndFS(3)
        maillage%Nfsys = maillage%FS%IndFS(4)
    endif
    
    ! Updating Mesh%Active
    do nc = 1,maillage%NBody
        if(nc.eq.1)then ! Tank
            maillage%Body(nc)%Active = .true.
        else
            maillage%Body(nc)%Active = fgeom_vect%Active(nc-1)
        end if
    end do
    
    9999 continue
        if(ierror/=0)then
            write(*,99) ierror
        endif
    99 format('** error #',i3,' : error change type mesh')
    
    if(allocated(id_noeud)) deallocate(id_noeud)
    if(allocated(id_tri))   deallocate(id_tri)
    if(allocated(table_corres)) deallocate(table_corres)

end subroutine change_type_mesh

subroutine Partial_initialization(Maillage,InputData)
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage), intent(inout)      :: Maillage     ! Mesh.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data.
    
    integer                             :: nc,jj        ! Loop parameters.
    

    ! This subroutine initializes partially the physical characteristics of the bodies before calling GeomInit for the first time.
    
    ! A better solution could be found. A second initialization is done in the subroutine Initialisation called just before beginning the temporal loop.
    
    ! Tank
    maillage%Body(1)%CMD = [.false., .false.]
    
    ! Bodies
    jj = 1
    do nc = Int_Body,maillage%NBody
        maillage%Body(jj)%CMD = [.true., .false.]
        maillage%Body(nc)%GBody(1:3) = InputData%Position(1:3,1,jj)
        if(iFixPoint)then
            maillage%Body(nc)%CSolv(1:3) = FixPointPos(1:3)
        else
            Maillage%Body(nc)%CSolv(1:3) = Maillage%Body(nc)%GBody(1:3)
        endif
        maillage%Body(nc)%CSolv(4:6) = InputData%Position(1:3,2,jj)
    
        maillage%Body(nc)%DimBody(1) = InputData%Lgeom(2,jj)
        
        jj = jj + 1
    end do
    
    do nc = Int_Body,maillage%NBody
        maillage%Body(nc)%IBody = 0._rp
        maillage%Body(nc)%VBody = 0._rp
        maillage%Body(nc)%ABody = 0._rp
    end do

end subroutine Partial_initialization

subroutine insert_mesh_global(mesh,nb_point,nb_tri,Maillage,nface,fgeom_vect,ierror,InputData,LocalRemeshFS)
    
    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(inout)           :: mesh                                 ! MGrid.
    integer,intent(in)                  :: nb_point, nb_tri                     !
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(inout)       :: Maillage                             ! Mesh.
    integer,intent(in)                  :: nface                                ! Number of faces.
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                           ! Geometries of the bodies.
    integer,intent(inout)               :: ierror                               ! Error flag.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                            ! Input data.
    logical,intent(in)                  :: LocalRemeshFS                        ! = true: remeshing the free surface, = false: no remeshing.
        
    integer                             :: j, k, n, jk, ind,jj,nc               !
    integer                             :: index_sl, Nnoeud, iface, ind_face0   !
    integer                             :: typeFrontiere, i1, i2, i3            !
    integer                             :: NCorps                               !
    integer                             :: np_sl                                !
    integer                             :: iP(3)                                !
    integer,dimension(nface)            :: ICorps                               !
    integer,dimension(:,:),allocatable  :: table_corres                         !
    integer,dimension(:),allocatable    :: id_noeud, id_tri                     !
    !f2py integer*1, dimension(1000)    :: vertex
    type(MVertex)                       :: vertex                               !
    integer                             :: id_nc,id_sl,itemp                    !
    !f2py integer*1, dimension(1000)    :: TNtemp
    type(TNoeud)                        :: TNtemp                               !
    !f2py integer*1, dimension(1000)    :: TFtemp
    type(TFacette)                      :: TFtemp                               !
    
    ! This subroutine recreates the structure TMaillage from MGrid and in adding the new floater meshes into the free surface mesh.
        
    ierror = 0
    index_sl = HouleRF%index
    
    ! NCorps
    NCorps = NBodies + 1 ! + 1 = Tank
            
    ! Definition of ICorps
    ICorps = 0._RP
    jj = 1
    do j = 1,NBodies+1
        do k = 1,fgeom_vect%nface_vect(j) ! nface = Somme des nface_vect
            if(j == 1)then ! Free surface and tank
                if (k == 1)then
                    ICorps(jj) = 0 ! Free surface
                else
                    ICorps(jj) = 1 ! Tank
                end if
            else ! Bodies
                ICorps(jj) = j
            end if
            jj = jj + 1
        end do
    end do
    allocate(table_corres(PointMax,nface))
    
    ! Change face number for non symetrical case
    do jj = 1,NBodies
        if(not(Symmetry) .and. InputData%igtype(jj).eq.5)then
            call merge_face(Mesh,nb_point,nb_tri,fgeom_vect%geom(jj)%axisym(1)%index,fgeom_vect%geom(jj)%axisym(2)%index,ierror)
            if(ierror/=0) goto 9999
        endif
    end do
    
    ! ---------------------------------------------------------------------------------
    ! Operation sur les noeuds
    ! ---------------------------------------------------------------------------------
    
    if(LocalRemeshFS)then ! Free surface remeshing
        Nnoeud = 0
    else
        Nnoeud = Maillage%Body(int_body)%IndBody(1)-1 ! Total number of nodes for the free surface and the tank.
    end if
    np_sl = 0
    
    ! Remplissage -------------------------------------------------------------------
    do j=1,nb_point
        vertex = mesh%point(j)
        n = vertex%nface
        do k=1,n
            iface = vertex%face(k)
                                 
            if(bottom_sym .and. iface ==  2) goto 1000
            if(Symmetry   .and. iface == -1) goto 1000
            
            if(iface.ne.-1)then ! Probably iface can be equal to -1 because the definition of the structure arete.
                
                if(not(LocalRemeshFS))then ! No free surface remeshing
                    if(ICorps(iface).lt.2) goto 1000
                end if
                
                Nnoeud = Nnoeud+1 ! Adding the nodes of the bodies.
                if(Nnoeud.gt.Pointmax)then
                    print*,"insert_mesh_global: The maximum number of nodes is reached in the remeshing."
                end if
                maillage%Tnoeud(Nnoeud)%Pnoeud(1:3) = vertex%coord(1:3)
                if(LocalRemeshFS)then
                    maillage%Tnoeud(Nnoeud)%indmesh = j
                end if
                table_corres(j,iface) = Nnoeud
                if(LocalRemeshFS)then ! Free surface remeshing
                    if(iface==index_sl)then
                        maillage%Tnoeud(Nnoeud)%typeNoeud = 0
                        maillage%Tnoeud(Nnoeud)%Npanneau = 0
                        id_sl = Nnoeud
                        np_sl = np_sl+1
                    else
                        ind = ICorps(iface)
                        maillage%Tnoeud(Nnoeud)%typeNoeud = 1
                        maillage%Tnoeud(Nnoeud)%Npanneau = ind
                    endif
                else ! No free surface remeshing
                    ind = ICorps(iface)
                    maillage%Tnoeud(Nnoeud)%typeNoeud = 1
                    maillage%Tnoeud(Nnoeud)%Npanneau = ind
                end if
                
            end if
            1000 continue
        enddo
    enddo
        
    ! Total number of nodes
    maillage%Nnoeud = Nnoeud
        
    ! Tri
    allocate(id_noeud(Nnoeud), id_tri(Nnoeud))
    do j=1,Nnoeud
        id_noeud(j) = j
    enddo
    
    if(not(LocalRemeshFS))then ! No free surface remeshing
        id_sl = maillage%FS%IndFS(3)
    end if
    
    nc = 0
    j = 1
    if(LocalRemeshFS)then
        maillage%FS%IndFS(1) = j
    end if
    do while(j .lt. id_sl) ! Free surface
        if(maillage%Tnoeud(j)%Npanneau > nc)then
            if(LocalRemeshFS)then
                TNtemp = maillage%Tnoeud(j)
                maillage%Tnoeud(j) = maillage%Tnoeud(id_sl)
                maillage%Tnoeud(id_sl) = TNtemp
            end if
            itemp = id_noeud(j)
            id_noeud(j) = id_noeud(id_sl)
            id_noeud(id_sl) = itemp
            do while(maillage%Tnoeud(id_sl)%Npanneau .ne. nc)
                id_sl = id_sl-1
            enddo
        endif
        j = j+1
    enddo 
    if(LocalRemeshFS)then
        maillage%FS%IndFS(3) = id_sl
    end if
        
    j = id_sl+1 ! First node of the tank
    do nc = 1,NCorps ! Tank and bodies
        maillage%Body(nc)%IndBody(1) = j
        id_nc = Nnoeud
        do while(maillage%Tnoeud(id_nc)%Npanneau .ne. nc)
            id_nc = id_nc-1
        enddo
        do while(j .lt. id_nc)
            if(maillage%Tnoeud(j)%Npanneau > nc)then
                TNtemp = maillage%Tnoeud(j)
                maillage%Tnoeud(j) = maillage%Tnoeud(id_nc)
                maillage%Tnoeud(id_nc) = TNtemp
                itemp = id_noeud(j)
                id_noeud(j) = id_noeud(id_nc)
                id_noeud(id_nc) = itemp
                do while(maillage%Tnoeud(id_nc)%Npanneau .ne. nc)
                    id_nc = id_nc-1
                enddo
            endif
            j = j + 1
        enddo
        maillage%Body(nc)%IndBody(3) = id_nc                             
        j = id_nc+1
    enddo
    
    ! -------------------------------------------------------------------------------------------
    ! Inversion du tableau id_noeud
    ! -------------------------------------------------------------------------------------------
    do j = 1,Nnoeud
        id_tri(id_noeud(j)) = j
    enddo 
    
    ! --------------------------------------------------------------------------------------------
    ! Operation sur les facettes
    ! --------------------------------------------------------------------------------------------    
    ind_face0 = Maillage%Body(int_body)%IndBody(2)-1
    
    ! Remplissage ------------------------------------------------------------------------------    
    do j=1,nb_tri
        if(LocalRemeshFS)then
            jk = j
        else
            jk = ind_face0 + j
        end if        
        typeFrontiere = mesh%tri(j)%typeFrontiere
        iface = mesh%tri(j)%face
        iP(1:3) = mesh%tri(j)%iP(1:3)
        i1 = table_corres(iP(1),iface)
        i2 = table_corres(iP(2),iface)
        i3 = table_corres(iP(3),iface)
        maillage%Tfacette(jk)%Tnoeud(1:3) = [id_tri(i1),id_tri(i2),id_tri(i3)] ! PYW
        maillage%Tfacette(jk)%typeFrontiere = typeFrontiere
        if(LocalRemeshFS)then
            if(iface == index_sl)then
                maillage%Tfacette(j)%Npanneau = 0
                maillage%Tfacette(j)%typeFrontiere = 0
                id_sl = j
            else
                ind = ICorps(iface) 
                maillage%Tfacette(j)%Npanneau = ind
                maillage%Tfacette(j)%typeFrontiere = 1
            endif
        else
            ind = ICorps(iface) 
            maillage%Tfacette(jk)%Npanneau = ind
            maillage%Tfacette(jk)%typeFrontiere = 1
        end if
    enddo
    if(LocalRemeshFS)then
        maillage%Nfacette = nb_tri
    else
        maillage%Nfacette = nb_tri + ind_face0 ! Number of panels of the bodies + number of panels of tank and the free surface.
    end if
        
    ! Tri
    if(LocalRemeshFS)then
        nc = 0
        j = 1
        maillage%FS%IndFS(2) = j ! Surface libre
        do while(j .lt. id_sl)            
            if(maillage%Tfacette(j)%Npanneau > nc)then
                TFtemp = maillage%Tfacette(j)
                maillage%Tfacette(j) = maillage%Tfacette(id_sl)
                maillage%Tfacette(id_sl) = TFtemp
                do while(maillage%Tfacette(id_sl)%Npanneau .ne. nc)
                    id_sl = id_sl-1
                enddo
            endif
            j = j+1
        enddo 
        maillage%FS%IndFS(4) = j
        j = j + 1
        
    else
        j = maillage%FS%IndFS(4)+1
    end if
    
    print*,""
    print*,maillage%FS%IndFS(4)
    print*,maillage%TFacette(maillage%FS%IndFS(4))%NPanneau
    print*,""
    if(maillage%TFacette(maillage%FS%IndFS(4))%NPanneau.eq.1)then
        pause
    end if
    
    do nc = 1,NCorps ! Bodies and tank
                
        maillage%Body(nc)%IndBody(2) = j   
        if(LocalRemeshFS)then
            id_nc = nb_tri
        else
            id_nc = maillage%Nfacette
        end if
        
        do while(maillage%Tfacette(id_nc)%Npanneau .ne. nc)
            id_nc = id_nc-1
        enddo
                
        do while(j .lt. id_nc)
            if(maillage%Tfacette(j)%Npanneau > nc)then
                TFtemp = maillage%Tfacette(j)
                maillage%Tfacette(j) = maillage%Tfacette(id_nc)
                maillage%Tfacette(id_nc) = TFtemp
                do while(maillage%Tfacette(id_nc)%Npanneau .ne. nc)
                    id_nc = id_nc-1
                enddo
            endif
            j = j + 1
        enddo
        
        maillage%Body(nc)%IndBody(4) = j  
        j = j + 1
    enddo   
        
    ! Mise a jour taille systeme 
    if(NCorps.gt.0)then
        maillage%Nsys = maillage%Body(NCorps)%IndBody(3)
        maillage%Nfsys = maillage%Body(NCorps)%IndBody(4)
    else
        maillage%Nsys = maillage%FS%IndFS(3)
        maillage%Nfsys = maillage%FS%IndFS(4)
    endif
    
    9999 continue
        if(ierror/=0)then
            write(*,99) ierror
        endif
    99 format('** error #',i3,' : error merge mesh')
    
    ! Deallocation  
    if(allocated(id_noeud)) deallocate(id_noeud)
    if(allocated(id_tri))   deallocate(id_tri)
    if(allocated(table_corres)) deallocate(table_corres)
    
end subroutine insert_mesh_global

subroutine create_dref_grid(mesh,nb_arete,xmin,xmax,ymin,ymax,&
&          xgrid,ygrid,matdref,nx,ny,iface,imeth,dx_out,ierror,InputData)
    
    !$ use OMP_LIB
    
    !f2py integer*1, dimension(1000)            :: mesh
    type(MGrid),intent(in)                      :: mesh                     ! Mesh.
    integer,intent(in)                          :: nb_arete                 ! Number of edges.
    real(rp),intent(in)                         :: xmin,xmax                ! Boundaries according to x in the Cartesian grid.
    real(rp),intent(in)                         :: ymin,ymax                ! Boundaries according to y in the Cartesian grid.
    integer,intent(in)                          :: nx,ny                    ! Number of points in the Cartesian grid.
    real(rp),dimension(nx),intent(inout)        :: xgrid                    ! Cartesian grid according to x.
    real(rp),dimension(ny),intent(inout)        :: ygrid                    ! Cartesian grid according to y.
    real(rp),dimension(nx,ny,3),intent(inout)   :: matdref                  ! Reference length for each cell of the Cartesian grid.
    integer,intent(in)                          :: iface                    ! Number of the face.
    integer,intent(in)                          :: imeth                    ! idref.
    real(rp),intent(in)                         :: dx_out                   ! Panel discretization.
    integer,intent(inout)                       :: ierror                   ! Error flag.
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(in)            :: InputData                ! Input data.

    integer                                     :: i1,i2,ind                ! Edge parameters.
    !f2py integer*1, dimension(1000)            :: edge
    type(MEdge)                                 :: edge                     ! Edge.
    integer                                     :: j,k,ii,jj,p              ! Loop parameters.
    real(rp)                                    :: dx,dy                    ! Spatial step.
    real(rp),dimension(2)                       :: M                        ! Point.
    real(rp),dimension(3)                       :: P1,P2                    ! Points
    integer                                     :: ndim                     ! nx*ny.
    integer,dimension(:),allocatable            :: ipiv                     ! Pivot indices
    real(rp),dimension(:,:),allocatable         :: A                        ! Matrix A of the linear system.
    real(rp),dimension(:),allocatable           :: B,Sol                    ! Vector B and solution of the linear system.
    real(rp)                                    :: precond, lnorm, a2, b2   ! Parameters of the linear approximation.
    real(rp),dimension(2)                       :: V                        ! Normal to an edge.
    integer                                     :: nNonZero                 ! Number the non-zero elements of A.
    integer,allocatable,dimension(:)            :: ja                       ! Contains column indices of the sparse matrix A (CSR3).
    integer,allocatable,dimension(:)            :: ia                       ! ia(i) points to the first column index of row i in the array ja (CSR3).
    real(rp),allocatable,dimension(:)           :: aa                       ! Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja.
    
    real(rp) :: debut,fin
    
    ! This subroutine creates the Cartesian grid and fills matdref.
    
    matdref = 0._rp
    ndim = nx*ny ! Only used with idref = 1 (Boundary problem).
    ierror = 0
    
    ! Creating of the Cartesian grid.
    dx = (xmax-xmin)/dble(nx-1)
    dy = (ymax-ymin)/dble(ny-1)
    do j = 1,nx
        xgrid(j) = xmin + (j-1)*dx
    end do
    do j = 1,ny
        ygrid(j) = ymin + (j-1)*dy
    end do
    
    ! Reference length for each panel of the Cartesian grid.
    select case (imeth) ! idref
    
    case(0)
        ! --------------------------------------------------------------
        !   Sparse Laplace equation : idref = 1
        ! --------------------------------------------------------------
        
        ! Assign boundary values.
        do j = 1,nb_arete
            edge = mesh%arrete(j) ! For the free surface, edges of the intersection curves and the boundary of the tank.
            i1 = edge%iP(1)
            i2 = edge%iP(2)
            P1 = mesh%point(i1)%coord
            P2 = mesh%point(i2)%coord            
            M = 0.5*(P1(1:2)+P2(1:2)) ! Center of the edge.
            call find_pos(xgrid,M(1),ii)
            call find_pos(ygrid,M(2),jj)
            matdref(ii,jj,1) = matdref(ii,jj,1) + edge%length ! 0 if the cell of the Cartesian grid does not intersect any part of the intersection curves or the boudary of the tank (case of the free surface).
            matdref(ii,jj,2) = matdref(ii,jj,2) + 1._rp
            
            ! Boundary layer of mesh with the same reference length as the floaters.
            if(d_bl.gt.Epsilon .and. edge%bf == 10)then ! On the intersection curves.
                
                ! Normal vector
                V = [P1(2)-P2(2),P2(1)-P1(1)]/edge%length
                M = M + d_bl*V
                if(M(1).le.xmax .and. M(1).ge.xmin .and. M(2).le.ymax .and. M(2).ge.ymin)then
                    call find_pos(xgrid,M(1),ii)
                    call find_pos(ygrid,M(2),jj)
                    matdref(ii,jj,1) = matdref(ii,jj,1) + edge%length ! edge%length = dx.
                    matdref(ii,jj,2) = matdref(ii,jj,2) + 1._rp
                end if
            end if
        end do
    
        where(matdref(:,:,2).gt.0)   
            matdref(:,:,1) = matdref(:,:,1)/matdref(:,:,2) ! Either 0 ou 1 in matdref.
        end where
    
        ! nNonZero.
        nNonZero = 0
        do k = 1,ny
            do j = 1,nx
                if(matdref(j,k,2)==0 .and. j.gt.1 .and. j.lt.nx .and.&
                    &  k.gt.1 .and. k.lt.ny) then
                    nNonZero = nNonZero + 5 ! 5 terms of the discretization of the Laplacian equation.
                elseif(matdref(j,k,2)==0)then ! Dans le cas ou le bord de la grille de ref n'est pas affecte de taille (cas solide).
                    nNonZero = nNonZero + 1 ! 1 term of the discretization of the Laplacian equation.
                else ! matdref(j,k,2)==1                    
                    nNonZero = nNonZero + 1 ! 1 term of the discretization of the Laplacian equation.
                end if
            end do
        end do
        
        ! Sol.
        allocate(Sol(ndim))
    
        ! B.
        allocate(B(ndim))
    
        ! aa.
        allocate(aa(nNonZero))
    
        ! ia.
        allocate(ia(ndim+1))
        ia(ndim+1) = nNonZero + 1 ! Number of non-zero elements in A plus one.
    
        ! ja.
        allocate(ja(nNonZero))
        ind = 0
        p = 1
        precond = -0.5_rp/(dx*dx+dy*dy)
        do j = 1,nx
            do k = 1,ny
                ind = ind + 1
                ia(ind) = p ! The first coefficient in the row ind is the p-th in aa.
                if(matdref(j,k,2)==0 .and. j.gt.1 .and. j.lt.nx .and.&
                    &  k.gt.1 .and. k.lt.ny) then
                    B(ind) = 0._rp
                    ! In ja, coeffcients must be in increasing order per row.
                    ja(p) = ind-ny 
                    ja(p+1) = ind-1
                    ja(p+2) = ind
                    ja(p+3) = ind+1
                    ja(p+4) = ind+ny
                    aa(p) = dx*dx*precond*1._rp   ! dref(i,j-1)
                    aa(p+1) = dy*dy*precond*1._rp ! dref(i-1,j)
                    aa(p+2) = 1._rp               ! dref(i,j)
                    aa(p+3) = dy*dy*precond*1._rp ! dref(i+1,j)
                    aa(p+4) = dx*dx*precond*1._rp ! dref(i,j+1)
                    p = p + 5
                elseif(matdref(j,k,2)==0)then ! Dans le cas ou le bord de la grille de ref n'est pas affecte de taille (cas solide).
                    B(ind) = dx_out
                    ja(p) = ind
                    aa(p) = 1._RP
                    p = p + 1
                else ! Reference size (dx1 on the tank, dx3 for the intersection curves, dx2 for the floaters).
                    B(ind) = matdref(j,k,1)
                    ja(p) = ind  
                    aa(p) = 1._RP
                    p = p + 1
                end if
            end do
        end do
    
        ! Parsido solver.
        call pardiso_solver(ja,ia,aa,nNonZero,ndim,ndim + 1,B,Sol,ierror)
        
        ! Distribution of the solution.
        ind = 0
        do j = 1,nx
            do k = 1,ny
                ind = ind+1
                if(matdref(j,k,2)==0)then
                    matdref(j,k,1) = Sol(ind)
                end if
            end do
        end do    
        
        ! Deallocating.
        if(allocated(ia))   deallocate(ia)
        if(allocated(ja))   deallocate(ja)
        if(allocated(aa))   deallocate(aa)
        if(allocated(B))    deallocate(B)
        if(allocated(Sol))  deallocate(Sol)
        
    case(1)
        ! --------------------------------------------------------------
        !   Full Laplace equation : idref = 1
        ! --------------------------------------------------------------
        
        ! Assign boundary values.
        do j = 1,nb_arete
            edge = mesh%arrete(j) ! For the free surface, edges of the intersection curves and the boundary of the tank.
            i1 = edge%iP(1)
            i2 = edge%iP(2)
            P1 = mesh%point(i1)%coord
            P2 = mesh%point(i2)%coord            
            M = 0.5*(P1(1:2)+P2(1:2)) ! Center of the edge.
            call find_pos(xgrid,M(1),ii)
            call find_pos(ygrid,M(2),jj)
            matdref(ii,jj,1) = matdref(ii,jj,1) + edge%length ! 0 if the cell of the Cartesian grid does not intersect any part of the intersection curves or the boudary of the tank (case of the free surface).
            matdref(ii,jj,2) = matdref(ii,jj,2) + 1._rp
            
            ! Boundary layer of mesh with the same reference length as the floaters.
            if(d_bl.gt.Epsilon .and. edge%bf == 10)then ! On the intersection curves.
                V = [P1(2)-P2(2),P2(1)-P1(1)]/edge%length
                M = M + d_bl*V
                if(M(1).le.xmax .and. M(1).ge.xmin .and. M(2).le.ymax .and. M(2).ge.ymin)then
                    call find_pos(xgrid,M(1),ii)
                    call find_pos(ygrid,M(2),jj)
                    matdref(ii,jj,1) = matdref(ii,jj,1) + edge%length ! edge%length = dx.
                    matdref(ii,jj,2) = matdref(ii,jj,2) + 1._rp
                end if
            end if
        end do
        
        where(matdref(:,:,2).gt.0)   
            matdref(:,:,1) = matdref(:,:,1)/matdref(:,:,2) ! Either 0 ou 1 in matdref.
        end where
        
        ! Allocate temporary table
        allocate(B(ndim))
        allocate(A(ndim,ndim))
        allocate(ipiv(ndim))
        allocate(Sol(ndim))
        
        B(1:ndim) = 0._rp
        A(1:ndim,1:ndim) = 0._rp
        
        precond = -0.5_rp/(dx*dx+dy*dy) 
        
        ! Assign Linear system.
        ind = 0
        do j = 1,nx
            do k = 1,ny
                ind = ind + 1
                if(matdref(j,k,2)==0 .and. j.gt.1 .and. j.lt.nx .and.&
                    &  k.gt.1 .and. k.lt.ny) then
                    B(ind) = 0._rp
                    A(ind,ind) = 1._rp                  ! dref(i,j)
                    A(ind,ind-1)  = dy*dy*precond*1._rp ! dref(i-1,j)
                    A(ind,ind+1)  = dy*dy*precond*1._rp ! dref(i+1,j)
                    A(ind,ind-ny) = dx*dx*precond*1._rp ! dref(i,j-1)
                    A(ind,ind+ny) = dx*dx*precond*1._rp ! dref(i,j+1)
                elseif(matdref(j,k,2)==0)then ! Dans le cas ou le bord de la grille de ref n'est pas affecte de taille (cas solide).
                    B(ind) = dx_out
                    A(ind,ind) = 1._rp
                else
                    B(ind) = matdref(j,k,1) ! Reference size (dx1 on the tank, dx3 for the intersection curves, dx2 for the floaters).
                    A(ind,ind) = 1._rp
                end if
            end do
        end do
        
        ! Solver.      
        call LU(A, B, Sol, ndim)
        
        ! Distribution of the solution.
        ind = 0
        do j = 1,nx
            do k = 1,ny
                ind = ind+1
                if(matdref(j,k,2)==0)then
                    matdref(j,k,1) = Sol(ind)
                end if
            end do
        end do    
        
        ! Deallocating.
        if(allocated(A))    deallocate(A) 
        if(allocated(ipiv)) deallocate(ipiv)
        if(allocated(B))    deallocate(B)
        if(allocated(Sol))  deallocate(Sol)
        
    case(2)
        ! ------------------------------------------------------------
        !   User define linear evolution : idref = 2
        ! ------------------------------------------------------------
        
        if(NBodies.gt.1)then
            print*,"WARNING: pay attention in using idref = 2 with several bodies. Linear evolution of the mesh is probably not appropriate!"
            print*,"WARNING: mesh parameters of the 1st body are used in the case of a linear evolution."
        end if
    
        if(iface.eq.1 .and. is_body)then ! Free surface.
            
            ! Rectangular tank.
            if(idtype.ne.2 .and. idtype.ne.4)then
                print*,"create_dref_grid shoud be adapted to rectangular tank!" ! Ldom(5) only works with circular tank. Lgeom(2) only works with circular floater.
                pause
            end if
            
            ! Linear approximation coefficients.
            if(is_immerged)then ! No intersection curve.
                a2 = (dx1-dx2Domain)/(Ldom(5)-InputData%Lgeom(2,1)) ! Ldom(5) : width of the tank, InputData%Lgeom(2,1) : width of the first floater.
                b2 = dx2Domain
			else ! It exists an intersection curve.
                a2 = (dx1-InputData%dx2(1))/(Ldom(5)-InputData%Lgeom(2,1))
                b2 = InputData%dx2(1)
			end if
            
            ! Initialization of matdref.
            do j = 1,nx
                do k = 1,ny
                    M = [xgrid(j),ygrid(k)]
                    lnorm = sqrt((M(1)-InputData%Position(1,1,1))**2 + (M(2)-InputData%Position(2,1,1))**2)
                    matdref(j,k,1) = max(a2*(lnorm-InputData%Lgeom(2,1))+b2,b2) ! Linear approximation
                    matdref(j,k,2) = 1
                end do
            end do
        elseif(iface.eq.2 .or. iface.eq.3)then
            do j = 1,nx
                do k=1,ny
                    matdref(j,k,1) = dx1
                    matdref(j,k,2) = 1
                end do
            end do  
        elseif(iface.ge.4)then
            do j=1,nx
                do k=1,ny
                    matdref(j,k,1) = InputData%dx2(1)
                    matdref(j,k,2) = 1
                end do
            end do
        endif
        
    end select
    
end subroutine create_dref_grid

subroutine init_cylindrical_mesh(mesh,nb_point,nb_arete,nb_tri,geom,iface,dx,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: mesh
  type(MGrid),intent(inout) :: mesh
  integer,intent(inout) :: nb_point,nb_arete,nb_tri
  !f2py integer*1, dimension(1000) :: geom
  type(disque2),intent(in) :: geom
  integer,intent(in) :: iface
  real(rp),intent(in) :: dx
  integer,intent(inout) :: ierror
! variables locales  
  real(rp) :: Orig(3), X, Y
  integer  :: np0, na0, Ntheta, jtheta, na1, j
  !f2py integer*1, dimension(1000) :: vertex,vertex_0
  type(MVertex) :: vertex,vertex_0
  !f2py integer*1, dimension(1000) :: edge
  type(MEdge)   :: edge
  !f2py integer*1, dimension(1000) :: tri
  type(MElement):: tri
  
  ierror = 0
  
  np0 = nb_point
  na0 = nb_arete
  Orig = geom%repere%origine
  
  vertex = [0._rp,0._rp,0._rp]
  vertex%nface = 1
  vertex%face(1) = iface
  vertex%bf = 0
  call add_element(mesh%point,nb_point,vertex)
  
  Ntheta = int(2*pi)
  do jtheta=1,Ntheta
    X = dx*cos(float(jtheta)-1.)
    Y = dx*sin(float(jtheta)-1.)
    vertex = [X, Y, 0._rp]
    vertex%nface = 1
    vertex%face=iface
    vertex%bf = -1
    call add_element(mesh%point,nb_point,vertex)
  enddo
  
  vertex_0 = mesh%point(np0+1)
  do j = np0+2,nb_point-1
    call edge_from_vertex(vertex_0,mesh%point(j),edge)
    edge%bf = 0
    call add_element(mesh,nb_arete,edge)
    
    call edge_from_vertex(mesh%point(j+1),mesh%point(j),edge)
    edge%bf = -1
    call add_element(mesh,nb_arete,edge)
  enddo
  call edge_from_vertex(vertex_0,mesh%point(nb_point),edge)
  edge%bf = 0
  call add_element(mesh,nb_arete,edge)
  call edge_from_vertex(mesh%point(np0+2),mesh%point(nb_point),edge)
  edge%bf = -1
  call add_element(mesh,nb_arete,edge)
  
  do jtheta = 1,Ntheta-1
    na1 = na0 + 2*(jtheta-1) + 1
    call triangle_from_edge(mesh%arrete(na1+2),mesh%arrete(na1+1),mesh%arrete(na1),tri)
    tri%face = iface
    call add_element(mesh%tri,nb_tri,tri)
  enddo
  na1 = na0 + 2*(Ntheta-1) + 1
  call triangle_from_edge(mesh%arrete(na0+1),mesh%arrete(na1+1),mesh%arrete(na1),tri)
  tri%face = iface
  call add_element(mesh%tri,nb_tri,tri)    
  
end subroutine init_cylindrical_mesh

subroutine smooth_corner(mesh,nb_point,ierror,InputData)
    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(inout) :: mesh
    integer,intent(inout) :: nb_point
    integer,intent(inout) :: ierror 
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data

    integer :: j,k
    integer :: index_sl, n
    real(rp) :: Rmax, Zmax, r, r1, z1, r2
    !f2py integer*1, dimension(1000) :: vertex
    type(MVertex) :: vertex
    logical :: bool
    real(rp),dimension(3) :: M, M2
    !real(rp),parameter :: Ra = 4.

    ierror = 0
    index_sl = HouleRF%index

    Rmax = InputData%Lgeom(2,1)-Ra
    Zmax = 0.5*InputData%Lgeom(1,1)-Ra
    print*,"Warning: smooth_corner uses the geometrical properties of the first floater."
    do j=1,nb_point
        vertex = mesh%point(j)
        n = vertex%nface
        bool = .false.
        do k=1,n
            bool = bool .or. vertex%face(k)==index_sl
        enddo
        if(.not.bool .and. sqrt(vertex%coord(1)**2 + vertex%coord(2)**2).lt.Ldom(5)-LAbs)then
            M = vertex%coord(1:3) - InputData%Position(1:3,1,1)
            r = sqrt(M(1)*M(1)+M(2)*M(2))
            if(r>=Rmax .and. abs(M(3))>=Zmax)then
                r1 = r-Rmax
                z1 = abs(M(3))-Zmax
                r2 = Ra/sqrt(1._rp+z1*z1/(r1*r1))
                M2(3) = Ra/sqrt(1._rp+r1*r1/(z1*z1))+Zmax
                M2(1) = (r2+Rmax)/r*M(1)
                M2(2) = (r2+Rmax)/r*M(2)
                vertex%coord(1:2) = M2(1:2) + InputData%Position(1:2,1,1)
                vertex%coord(3)   = InputData%Position(3,1,1) + sign(1.,M(3))*M2(3)
                mesh%point(j) = vertex
            endif
        endif
    enddo

end subroutine smooth_corner

function fz(M,imetric,param) result(res)
    implicit none
    real(rp),dimension(2),intent(in) :: M
    integer,intent(in) :: imetric
    real(rp),dimension(*),intent(in) :: param
    real(rp) :: res
    integer :: n

    select case (imetric)
      case(1)
        call fz_sphere(M,param(1),res)
      case(2)
        call fz_cone(M,param(1),param(2),res)
      case(3)
        n = int(param(1))
        call fz_axisym(M,n,param(2:n+1),param(n+2:2*n+1),res)
      case(4)
        call Hull_function(M(1),M(2),res)
        res = 0.5_RP*res
      case default 
        res = 0._rp
    endselect

end function fz

subroutine fz_sphere(M,r,z)
    implicit none
    real(rp),dimension(2),intent(in)    :: M
    real(rp),intent(in)                 :: r
    real(rp),intent(out)                :: z
    
    if(r*r-(M(1)*M(1)+M(2)*M(2)) .lt. Epsilon)then
        z = 0._RP
    else
        z = sqrt(r*r-(M(1)*M(1)+M(2)*M(2)))
    end if
    
end subroutine fz_sphere

subroutine fz_cone(M,a,b,z)
    implicit none
    real(rp),dimension(2),intent(in) :: M
    real(rp),intent(in) :: a,b
    real(rp),intent(out) :: z
    z = a*M(2)+b
end subroutine fz_cone

subroutine fdir_cone(A,B,M,param,dir)
    implicit none   
    real(rp),dimension(3),intent(in) :: A,B
    real(rp),dimension(3),intent(in) :: M
    real(rp),dimension(*),intent(in) :: param
    real(rp),dimension(3),intent(out) :: dir
!   local
    real(rp),dimension(2) :: P1,P2
    real(rp) :: rv,dist
    real(rp),dimension(2) :: vect,sm,sol
    real(rp),dimension(2,2) :: mat
    integer :: ierror
    integer :: saving_A

    rv = param(1)*M(2)+param(2)
    P1 = [M(1)*param(5)/rv,param(7)]
    P2 = [M(1)*param(6)/rv,param(8)]
    vect = P2-P1
    dist = sqrt(dot_product(vect,vect)+(param(5)-param(6))**2)
    mat(:,1) = [1._rp,0._rp]
    mat(:,2) = vect/dist    
    sm = B(1:2)-A(1:2)
    sm = sm/norm2(sm)

    Sol = 0._rp
    saving_A = 1
    call LU(mat,sm,Sol,2,ierror,saving_A)

    dir(1:2) = -Sol(2)*mat(:,1) + Sol(1)*mat(:,2)  
    dir(3) = 0._rp  

end subroutine fdir_cone

subroutine norm_metric(A,B,imetric,param,val,fgeom)
    
    real(rp),dimension(3),intent(in) :: A,B
    integer,intent(in) :: imetric
    real(rp),dimension(*),intent(in) :: param
    real(rp),intent(out) :: val
    !f2py integer*1, dimension(1000):: fgeom
    type(type_geom), optional :: fgeom
    
    integer :: n

    select case (imetric)
        case(2)
            call norm_cone(A(1:2),B(1:2),param,val)
        case(3)
            n = int(param(1))
            call norm_axisym2(A(1:2),B(1:2),n,param(2:n+1),param(n+2:2*n+1),val)
        case(4)
            call norm_wigley(A(1:2),B(1:2),fgeom%wigley(1),val)
        case default 
            val = norm2(B(1:3) - A(1:3))
    end select
    
end subroutine norm_metric

subroutine norm_cone(A,B,param,val)
    implicit none
    real(rp),dimension(2),intent(in) :: A,B
    real(rp),dimension(*),intent(in) :: param
    real(rp),intent(out) :: val
!   local
    real(rp),dimension(3) :: M,N
    real(rp) :: rA,thetaA,rB,thetaB

    rA = param(1)*A(2)+param(2)
    thetaA = A(1)/rA
    rB = param(1)*B(2)+param(2)
    thetaB = B(1)/rB

    M = [rA*cos(thetaA),rA*sin(thetaA),A(2)]
    N = [rB*cos(thetaB),rB*sin(thetaB),B(2)]

    val = norm2(M-N)

    if(abs(thetaA-thetaB).gt.PI)then
        val = 999.
    endif
    
end subroutine norm_cone

subroutine test_prox_cone(A,B,Q,param,delta,bool)
    implicit none
    real(rp),dimension(2),intent(in) :: A,B,Q
    real(rp),dimension(*),intent(in) :: param
    real(rp),intent(in) :: delta
    logical,intent(out) :: bool
!   local
    real(rp),dimension(3) :: M,N,P,pv1
    real(rp) :: rA,rB,rQ
    real(rp) :: thetaA,thetaB,thetaQ    
    real(rp) :: dist,ps1,h
!    
    rA = param(1)*A(2)+param(2)
    thetaA = A(1)/rA
    rB = param(1)*B(2)+param(2)
    thetaB = B(1)/rB
    rQ = param(1)*Q(2)+param(2)
    thetaQ = Q(1)/rQ
!
    M = [rA*cos(thetaA),rA*sin(thetaA),A(2)]
    N = [rB*cos(thetaB),rB*sin(thetaB),B(2)]
    P = [rQ*cos(thetaQ),rQ*sin(thetaQ),Q(2)]
!
    dist = norm2(N-M)
    ps1 = dot_product(P-M,N-M)
    call Computation_vect_product(M-N,P-M,pv1)
    h = norm2(pv1)/dist
!
    !bool = ps1 .le. dist+Epsilon .and. dist .gt. -Epsilon .and. abs(h) .le. 0.1*delta 
    bool = abs(h).le.0.1*delta

    M = [A(1),A(2),0._rp]
    N = [B(1),B(2),0._rp]
    P = [Q(1),Q(2),0._rp]
    dist = norm2(N-M)
    ps1 = dot_product(P-M,N-M)
    call Computation_vect_product(M-N,P-M,pv1)
    h = norm2(pv1)/dist
    bool = bool .or. norm2(pv1).le.0.001
!
end subroutine test_prox_cone

subroutine adjacent_front_edge(mesh,nb_arete,front_edge,near_point,&
&                              angle_maxA,angle_maxB,ierror)
    
    !f2py integer*1, dimension(1000)                            :: mesh
    type(MGrid),intent(in)                                      :: mesh
    integer,intent(in)                                          :: nb_arete
    !f2py integer*1, dimension(1000)                            :: front_edge
    type(pile_element),pointer                                  :: front_edge
    !f2py integer*1, dimension(1000)                            :: near_point
    type(near_point_cell),pointer                               :: near_point
    real(rp),intent(inout)                                      :: angle_maxA, angle_maxB
    integer,intent(inout)                                       :: ierror
    
    integer                                                     :: index_arr, ind, j
    integer                                                     :: index_A, index_B
    !f2py integer*1, dimension(1000)                            :: active_edge,local_edge
    type(MEdge)                                                 :: active_edge, local_edge
    real(rp),dimension(3)                                       :: A, B, C
    real(rp)                                                    :: theta, cosA, cosB, vectA, vectB
    real(rp)                                                    :: d1, d2, d3
    !f2py integer*1, dimension(1000)                            :: ptr_ed
    type(pile_element) ,pointer                                 :: ptr_ed
    integer                                                     :: na, nb
    integer,dimension(10)                                       :: A_adj, B_adj
    !f2py integer*1, dimension(1000)                            :: ptr,local_pile
    type(near_point_cell),pointer                               :: ptr, local_pile
    real(rp),parameter                                          :: theta_max = 1.3089969389957472_rp ! 75 deg in radians.
    
    ! This subroutine tests if two adjacent fronts can form a new triangle.
    ! This algorithm is based on the part 3.4.3.2 of the paper of Tristano ("Advancing front surface mesh generation in parametric space using a Riemannian surface definition).
    
    ierror = 0
    na = 0 ; nb = 0
    
    index_arr = front_edge%val
    active_edge = mesh%arrete(index_arr)
    
    index_A = active_edge%iP(1)
    index_B = active_edge%iP(2)
    
    ! Loop on the front edge to find adjacent vertex
    ptr_ed => front_edge
    do while(associated(ptr_ed))
        ind = ptr_ed%val
        local_edge = mesh%arrete(ind)
        if(ind/=index_arr)then
            if(local_edge%iP(1)==index_A)then
                na = na+1
                A_adj(na) = local_edge%iP(2)
            elseif(local_edge%iP(2)==index_A)then
                na = na+1
                A_adj(na) = local_edge%iP(1)
            elseif(local_edge%iP(1)==index_B)then
                nb = nb+1
                B_adj(nb) = local_edge%iP(2)
            elseif(local_edge%iP(2)==index_B)then
                nb = nb+1
                B_adj(nb) = local_edge%iP(1)
            endif
        endif
        ptr_ed => ptr_ed%suiv
    enddo
    
    ! Test the location of adjacent vertex (angle < 75°).
    
    A = mesh%point(index_A)%coord(1:3)
    B = mesh%point(index_B)%coord(1:3)
    
    d1 = active_edge%length

    angle_maxA = PI
    angle_maxB = PI

    nullify(local_pile)
    do j=1,na
        ind = A_adj(j)
        C = mesh%point(ind)%coord
        d2 = norm2(C-A)
        d3 = norm2(C-B)
        cosA = (d2*d2+d1*d1-d3*d3)/(2._rp*d2*d1)
        vectA = (B(1)-A(1))*(C(2)-A(2))-(B(2)-A(2))*(C(1)-A(1))
        if(abs(cosA).lt.1._rp)then
            theta = acos(cosA) ! Equation 8 of Tristano.
            angle_maxA = min(theta,angle_maxA)           
        else
            theta = 999.
        endif
        if(theta .lt. theta_max .and. vectA .gt. Epsilon)then
            call add_vertex_pile(local_pile,mesh%point(ind),theta)
        endif        
    enddo
    
    do j=1,nb
        ind = B_adj(j)
        C = mesh%point(ind)%coord
        d2 = norm2(C-B)
        d3 = norm2(C-A)
        cosB = (d2*d2+d1*d1-d3*d3)/(2._rp*d2*d1)
        vectB = (B(1)-A(1))*(C(2)-B(2))-(B(2)-A(2))*(C(1)-B(1))
        if(abs(cosB).lt.1._rp)then
            theta = acos(cosB)
            angle_maxB = min(theta,angle_maxB)
        else
            theta = 999.
        endif
        if(theta .lt. theta_max .and. vectB .gt. Epsilon)then
            call add_vertex_pile(local_pile,mesh%point(ind),theta)
        endif        
    enddo
    
    ! Add vertex to the top of near point pile
    ptr => local_pile
    if(associated(ptr))then
        do while (associated(ptr%suiv))
            ptr => ptr%suiv
        enddo
        ptr%suiv => near_point
        near_point => local_pile 
    endif

    9999 continue
        if(ierror/=0)then
            write(*,99) ierror
        endif
    99 format('** error #',i3,' : cannot find adjacent front vertex')   

    nullify(ptr,local_pile)         
    
end subroutine adjacent_front_edge

subroutine internal_node2D(mesh,P_A,P_B,P_C,front_vertex,bool,ierror)
    
    !f2py integer*1, dimension(1000) :: mesh
    type(MGrid),intent(in) :: mesh
    !f2py integer*1, dimension(1000) :: P_A,P_B,P_C
    type(MVertex),intent(in) :: P_A, P_B, P_C
    !f2py integer*1, dimension(1000) :: front_vertex
    type(pile_element),pointer :: front_vertex
    logical,intent(inout) :: bool
    integer,intent(inout) :: ierror
    
    integer :: ind
    real(rp),dimension(2) :: A, B, C, Q
    !f2py integer*1, dimension(1000) :: ptr
    type(pile_element),pointer :: ptr
    real(rp),parameter :: max_tol = -0.0154321012_rp
    real(rp) :: a1,b1,c1,a2,b2,c2,detT,inv_detT
    real(rp),dimension(3) :: Sol
    
    ierror = 0
    
    A = P_A%coord(1:2)
    B = P_B%coord(1:2)
    C = P_C%coord(1:2)
    
    bool = .false.
    ptr => front_vertex
    
    detT = (B(2)-C(2))*(A(1)-C(1))+(C(1)-B(1))*(A(2)-C(2))
    if(abs(detT).gt.Epsilon)then
        inv_detT = 1._rp/detT
    else
        print*,""
        print*,"A: ",A(1)," ",A(2)
        print*,"B: ",B(1)," ",B(2)
        print*,"C: ",C(1)," ",C(2)
        print*,"detT: ",detT
        print*,"internal_node2D: cannot verify internal node criteria"
        print*,""
        goto 9999
    endif
     
    a1 = inv_detT*(B(2)-C(2))
    b1 = inv_detT*(C(1)-B(1))
    c1 = -a1*C(1)-b1*C(2)
    a2 = inv_detT*(C(2)-A(2))
    b2 = inv_detT*(A(1)-C(1))
    c2 = -a2*C(1)-b2*C(2)
    
    do while (associated(ptr) .and. .not. bool)
        ind = ptr%val
        Q = mesh%point(ind)%coord(1:2)   
        if(ind/=P_A%index .and. ind/=P_B%index .and. ind/=P_C%index)then
            Sol(1) = a1*Q(1)+b1*Q(2)+c1
            Sol(2) = a2*Q(1)+b2*Q(2)+c2
            Sol(3) = 1-Sol(1)-Sol(2)
            bool = all(Sol>max_tol)
        endif       
        ptr=>ptr%suiv
    enddo
    
    9999 continue
    
end subroutine internal_node2D

subroutine limit_angle(A,B,Q,angle_maxA,angle_maxB,is_angle,ierror)
    implicit none
    real(rp),dimension(3),intent(in) :: A, B, Q
    real(rp),intent(in) :: angle_maxA, angle_maxB
    logical,intent(inout) :: is_angle
    integer,intent(inout) :: ierror
!   local
    real(rp) :: d1, d2 ,d3
    real(rp) :: cosA, cosB, thetaA, thetaB
!
    ierror = 0
    is_angle = .false.
!    
    d1 = norm2(Q-B)
    d2 = norm2(Q-A)
    d3 = norm2(B-A)
!    
    cosA = (d2*d2+d3*d3-d1*d1)/(2._rp*d2*d3)
    if(abs(cosA).lt.1._rp)then
        thetaA = acos(cosA)
    else
        thetaA = 999.
    endif
!
    cosB = (d1*d1+d3*d3-d2*d2)/(2._rp*d1*d3)
    if(abs(cosB).lt.1._rp)then
        thetaB = acos(cosB)
    else
        thetaB = 999.
    endif
!    
    is_angle = thetaA .gt. angle_maxA+Epsilon .or. thetaB .gt. angle_maxB+Epsilon
!
9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
99 format('** error #',i3,' : cannot verify limit angle')   
!    
end subroutine limit_angle

! ----------------------------------------------------------------------
!   ESTIMATE POSITION : 
! ----------------------------------------------------------------------
subroutine estimate_position(A,B,h_3D,imetric,param,P,ierror,fgeom)
    
    real(rp),dimension(3),intent(in)        :: A, B
    real(rp),intent(inout)                  :: h_3D
    integer,intent(in)                      :: imetric
    real(rp),dimension(*),intent(in)        :: param
    real(rp),dimension(3),intent(inout)     :: P
    integer,intent(inout)                   :: ierror
    !f2py integer*1, dimension(1000)        :: fgeom
    type(type_geom),intent(in)              :: fgeom                    ! Geometry.
    
    integer                                 :: istat
    real(rp)                                :: hr, x_start, h_2D, yhull
    real(rp),dimension(2)                   :: N_AB, V_AB, N_2D
    real(rp),dimension(3)                   :: M0, M
    real(rp),dimension(2,4)                 :: RK_N2D
    real(rp),dimension(4)                   :: RK_hr
    logical                                 :: high_order
    real(rp)                                :: h_3D_modif               ! h_3d modified in case of high curvature surface (ex at the boundary between two half spheres).
    
    ! This subroutine searches a new point in the disc of radius dist3.
    
    ierror = 0
    
    N_AB = [A(2)-B(2),B(1)-A(1)]
    N_AB = N_AB / norm2(N_AB)
    
    V_AB = B(1:2)-A(1:2)
    V_AB = V_AB / norm2(V_AB)
    
    M0 = 0.5_rp*(A+B)
    
    call fdir2D(M0,A,B,imetric,param,N_2D,hr,high_order,istat,ierror,fgeom)
    if(ierror/=0) goto 9999
        
    if(high_order .and. istat==0)then
        
        ! Integration with Runge-Kutta if high order needed
        RK_N2D(:,1) = N_2D
        RK_hr(1) = hr
        
        M(1:2) = M0(1:2) + 0.5_rp*h_3D/hr*N_2D
        call fdir2D(M,A,B,imetric,param,N_2D,hr,high_order,istat,ierror,fgeom)
        if(ierror/=0) goto 9999
        RK_N2D(:,2) = N_2D
        RK_hr(2) = hr
        
        if(istat==0)then
            M(1:2) = M0(1:2) + 0.5_rp*h_3D/hr*N_2D
            call fdir2D(M,A,B,imetric,param,N_2D,hr,high_order,istat,ierror,fgeom)
            if(ierror/=0) goto 9999
            RK_N2D(:,3) = N_2D
            RK_hr(3) = hr
        endif
        
        if(istat==0)then
            M(1:2) = M0(1:2) + h_3D/hr*N_2D
            call fdir2D(M,A,B,imetric,param,N_2D,hr,high_order,istat,ierror,fgeom)
            if(ierror/=0) goto 9999
            RK_N2D(:,4) = N_2D
            RK_hr(4) = hr
        endif
        
        if(istat==0)then
            P(1:2) = M0(1:2) + h_3D/6._rp*(RK_N2D(:,1)/RK_hr(1)+2.*RK_N2D(:,2)/RK_hr(2)+2.*RK_N2D(:,3)/RK_hr(3)+RK_N2D(:,4)/RK_hr(4))
        else
            P(1:2) = M0(1:2) + h_3D*RK_N2D(:,1)/RK_hr(1)
        endif
        
    elseif(istat==1 .or. istat==2)then
        
        if(idebug>0) write(*,*) 'dicho...'
        x_start = h_3D
        call dicho_fz(fz,M0,N_2D,x_start,h_3D,Epsilon,imetric,param,h_2D)
        
        if(h_2D .lt. -900.)then
            ierror = 100
            goto 9999
        endif
        P(1:2) = M0(1:2) + h_2D*N_2D
        
    else
        ! Simple first order integration in other cases  
        P(1:2) = M0(1:2) + h_3D/hr*N_2D
        
    endif       

    if (imetric.eq.4) then
        call Hull_function(2._RP*P(1)/fgeom%Wigley(1)%L, P(2)/fgeom%Wigley(1)%D,yhull)
        P(3) = 0.5_RP*yhull*fgeom%Wigley(1)%B
    else
        P(3) = fz(P(1:2),imetric,param)
    end if
9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
99 format('** error #',i3,' : cannot determine front direction')        
    
end subroutine estimate_position

subroutine dicho_fz(fz,M0,dir,s0,valf,tol,imetric,param,s)
    implicit none
    real(rp) :: fz
    real(rp),dimension(2),intent(in) :: M0
    real(rp),dimension(2),intent(in) :: dir
    real(rp),intent(in) :: s0,valf,tol
    integer,intent(in) :: imetric
    real(rp),dimension(*),intent(in) :: param
    real(rp),intent(out) :: s
!   local
    integer :: n
    real(rp),dimension(2) :: M
    real(rp) :: f0,fi,s1,s2,df
    integer,parameter :: itmax = 8000

    f0 = fz(M0,imetric,param)
  
    M = M0+s0*dir
    fi = fz(M,imetric,param)

    df = fi-f0

    n = 0 ; s1 = 0._rp ; s2 = s0
    
    do while(abs(df-valf).gt.tol .and. n.le.itmax)
        s = 0.5_rp*(s1+s2)
        M = M0+s*dir
        fi = fz(M,imetric,param)
        df = fi-f0
        if(df > valf)then
            s2 = s
        else
            s1 = s
        endif
        n = n+1
    enddo
    
    if(n==itmax)then
        write(*,*) 'warning : end of loop obtained during dichotomy'
    endif

end subroutine dicho_fz

! ---------------------------------------------------------------
!   FDIR2D
! ---------------------------------------------------------------
subroutine fdir2D(M,A,B,imetric,param,N_2D,hr,high_order,istat,ierror,fgeom)
    
    real(rp),dimension(3),intent(in) :: M, A, B
    integer,intent(in) :: imetric
    real(rp),dimension(*),intent(in) :: param
    real(rp),dimension(2),intent(inout) :: N_2D
    real(rp),intent(inout) :: hr
    logical,intent(inout) :: high_order
    integer,intent(inout) :: istat
    integer,intent(inout) :: ierror
    !f2py integer*1, dimension(1000) :: fgeom
    type(type_geom), optional :: fgeom
    
    integer :: n
    real(rp),dimension(3) :: V3, N_3D
    real(rp),dimension(2) :: V2
    real(rp),dimension(2,2) :: matR, matN
    real(rp),dimension(3) :: vect_product_1
    
    ierror = 0
    istat = 0
    high_order = .false.
    hr = 1._rp
         
    select case (imetric)
    
    case(1)
        
        V3 = B-A
        V3 = V3 / norm2(V3)
        call Computation_vect_product(V3,-M,vect_product_1)
        N_3D =  vect_product_1/norm2(M)
        N_2D = N_3D(1:2)
        
    case(2)
        
        call fdir_cone(A,B,M,param,N_3D)
        N_2D = N_3D(1:2)
        
    case(3)
        
        n = int(param(1))
        call matrix_riemann_axisym(M(1:2),n,param(2:n+1),param(n+2:2*n+1),matR,istat)
        matN(1,1:2) = -matR(2,1:2)
        matN(2,1:2) = matR(1,1:2)
        V2 = B(1:2)-A(1:2)
        V2 = V2 / norm2(V2)
        N_2D = matmul(matN,V2)
        call norm_riemann(matR,N_2D,hr)
        if(abs(hr).gt.90000)then
            ierror = 100
            goto 9999
        endif     

        high_order = matR(1,1).gt.dS2max .or. matR(2,2).gt.dS2max
        
    case(4)
        call matrix_riemann_Wigley(M(1:2),fgeom%wigley(1),matR,istat)
        matN(1,1:2) = -matR(2,1:2)
        matN(2,1:2) = matR(1,1:2)
        V2 = B(1:2)-A(1:2)
        V2 = V2 / norm2(V2)
        N_2D = matmul(matN,V2)
        call norm_riemann(matR,N_2D,hr)
        if(abs(hr).gt.90000)then
            ierror = 100
            goto 9999
        endif     
        
        high_order = matR(1,1).gt.dS2max .or. matR(2,2).gt.dS2max
    case default
    
        N_2D = [A(2)-B(2),B(1)-A(1)]
        N_2D = N_2D / norm2(N_2D)
        
    end select
    
9999 continue
    if(ierror/=0)then
        write(*,99) ierror
    endif
99 format('** error #',i3,' : cannont estimate direction 2D')
!
end subroutine fdir2D

subroutine merge_face(Mesh,nb_point,nb_tri,iface1,iface2,ierror)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(MGrid),intent(inout)           :: Mesh             ! Mgrid
    integer,intent(in)                  :: nb_point,nb_tri  ! Number of points and triangles
    integer,intent(in)                  :: iface1,iface2    ! Index of the first and second axisym parts (there are always two parts in a axisym geometry)
    integer,intent(inout)               :: ierror           ! Error flag
    
    !f2py integer*1, dimension(1000)    :: vertex
    type(MVertex)                       :: vertex           ! Vertex
    integer                             :: j,k,n,ncom       ! Loop parameters
    integer,dimension(nfvmax)           :: fcom             ! Table
    
    ! This subroutine merges the number of faces and the faces themselves of the points of the second part of a axisym geometry in the MGrid structure.
    
    ierror = 0

    do j=1,nb_point
        vertex = mesh%point(j)
        n = vertex%nface
        call common_int(vertex%face,n,[iface1,iface2],2,fcom,ncom)
        if(ierror/=0) goto 9999
        
        ! New definition of nface and face for vertex
        if(ncom==1)then ! All the points of the two parts of the axisym geometry.
            if(fcom(1)==iface2)then ! Only the points of the second parts.
                do k=1,n
                    if(vertex%face(k)==iface2)then
                        vertex%face(k) = iface1 ! They are now linked to the first part.
                    endif
                enddo
            endif
        elseif(ncom==2)then ! Points linked to the two parts of the axisym geometry (closing points)
            do k=1,n
                if(vertex%face(k)==iface2)then
                    vertex%nface = vertex%nface-1 ! They have one less face.
                    if(k>1 .and. k<n)then ! Between the first and the last face.
                        vertex%face(1:k-1) = vertex%face(1:k-1)
                        vertex%face(k:n-1) = vertex%face(k+1:n)
                    elseif(k==1)then ! First face
                        vertex%face(1:n-1) = vertex%face(2:n)
                    elseif(k==n)then ! Last face
                        vertex%face(n) = -9
                    endif
                 endif
            enddo
        endif
        
        ! Copy of vertex into mesh
        if(ncom>0)then
            ! nface
            mesh%point(j)%nfaceAxisym = mesh%point(j)%nface
            mesh%point(j)%nface = vertex%nface
            
            ! face
            if(not(allocated(mesh%point(j)%faceAxisym))) allocate(mesh%point(j)%faceAxisym(nvmax))
            mesh%point(j)%faceAxisym(:) = 0._RP
            do k = 1,mesh%point(j)%nfaceAxisym
                mesh%point(j)%faceAxisym(k) = mesh%point(j)%face(k)
            end do
            mesh%point(j)%face(:) = vertex%face(:)
        else
            mesh%point(j)%nfaceAxisym = 0
        endif
        
    enddo  

    do j=1,nb_tri
        if(mesh%tri(j)%face==iface2)then
            mesh%tri(j)%face = iface1
        endif
    enddo

    9999 continue
        if(ierror/=0)then
            write(*,99) ierror
        endif
    99 format('** error #',i3,' : error merge face')    

end subroutine merge_face

subroutine demerge_face(Mesh,nb_point,nb_tri,iface1,iface2,ierror)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(MGrid),intent(inout)           :: Mesh             ! Mgrid
    integer,intent(in)                  :: nb_point,nb_tri  ! Number of points and triangles
    integer,intent(in)                  :: iface1,iface2    ! Index of the first and second axisym parts (there are always two parts in a axisym geometry)
    integer,intent(inout)               :: ierror           ! Error flag
    
    integer                             :: j,k              ! Loop parameters
        
    ! This subroutine separates again the face of the two parts of the axisym geometries.
    
    ierror = 0
    do j=1,nb_point
        if(mesh%point(j)%nfaceAxisym .ne. 0)then ! Former second part of an axisym geometry.
            mesh%point(j)%nface = mesh%point(j)%nfaceAxisym
            do k = 1,mesh%point(j)%nfaceAxisym
                mesh%point(j)%face(k) = mesh%point(j)%faceAxisym(k)
            end do
        end if
    enddo  
    
end subroutine demerge_face


end module MeshFonct 
  
