module MeshModule
use Constantes
use FonctionsCommunes
use PrePlot
use MAdvanceFront
use GeomMesh
use Structuresdonnees
implicit none

! ------------------------------------------------------------------------
!               Mesh module
! ------------------------------------------------------------------------

contains

  subroutine mesh0D(mesh,nb_point,fgeom,t,iflag,ierror)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh     ! MGrid.
    integer,intent(inout)               :: nb_point
    !f2py integer*1, dimension(1000)    :: fgeom
    type(type_geom),intent(in)          :: fgeom    ! Geometry.
    real(rp),intent(in)                 :: t        ! Current time.
    integer,intent(in)                  :: iflag    ! Intersection curves present or not.
    integer,intent(inout)               :: ierror   ! Error flag.
    
    integer                             :: npoint   ! Number of points in the geometry.
    integer                             :: j        ! Loop parameter.
    !f2py integer*1, dimension(1000)    :: P
    type(point)                         :: P        ! Point.
    real(rp)                            :: eta      ! Wave elevation.
    !f2py integer*1, dimension(1000)    :: vertex
    type(MVertex)                       :: vertex   ! Vertex.
    integer                             :: naux     !
    integer,dimension(1)                :: aux      !
    
    ! This subroutine creates mesh0D (points under or at the sea level).
    
    ierror = 0
    ! iflag = 0 when creating the mesh, 1 in case of remeshing.
    
    npoint = fgeom%npoint

    
    do j=1,npoint
        P = fgeom%point(j)

        ! Local to cartesian coordinates.
        call loc2cart(P%coord,fgeom%repere,P%coord)
        
        ! Wave elevation at the points of the geometry.
        call CEta0(P%coord,t,eta)        
        

        ! Points on the free surface (z = 0).
        call common_int(P%face,P%nface,[HouleRF%index],1,aux,naux,ierror) ! HouleRF%index is fixed equal to 1.
        
        if((P%coord(3)<eta .or. abs(P%coord(3)-eta).lt.Epsilon .or. naux==1).and.&
            & (iflag==0 .or. (iflag==1 .and. naux==0)) )then ! If the points is under the free surace or if the points is at the free surface and if the mesh is creating or if it is a remeshing, just consider the points under the sea level.
            
            ! Points to Mvertex.
            call assign_mvertex_point(vertex,P) 
            
            ! Updating the vertical position for the points on the sea level.
            if(naux==1) vertex%coord(3) = eta 
            
            ! Updating bf
            vertex%bf = 1 ! bf = 1: free surface.
            if(abs(P%coord(2)).lt.Epsilon .and. Symmetry)then
                vertex%bf = 2 ! bf = 2: points on the plan of symmetry y = 0.
            endif
            
            ! Adding the vertex to the MGrid structure.
            call add_element(mesh%point,nb_point,vertex)
            
        endif
    enddo
    
  end subroutine mesh0D 

! *********************************************************************
! Mesh intersection SL/SM
! *********************************************************************
  subroutine mesh_inter(mesh,nb_point,nb_arete,dx,ierror,tab,n_tab)
    
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh             ! Mgrid
    integer,intent(inout)               :: nb_point,nb_arete! Number of points and edges
    real(rp),intent(in)                 :: dx               ! Size of the panel
    integer,intent(inout)               :: ierror           ! Error flag
    !f2py integer*1, dimension(1000)    :: tab
    type(chaine_point_pt),dimension(:)  :: tab              ! Table of the pointers toward the intersection points
    integer                             :: n_tab            ! Number of intersection curve
    
    !f2py integer*1, dimension(1000)    :: domaine2
    type(domaine)                       :: domaine2         ! Geometry of the domain (local variable)
    !f2py integer*1, dimension(1000)    :: P
    type(point)                         :: P                ! Point
    !f2py integer*1, dimension(1000)    :: ptr0,ptr1,ptr2
    type(chaine_point),pointer          :: ptr0,ptr1,ptr2   ! Chain of points
    integer                             :: j                ! Parameter loop
    integer                             :: ndim,na0,np0     !
    logical                             :: is_close         ! True if the intersection curve is closed, false otherwise
    integer                             :: nPts             ! Number of points on the intersection curves.
    
    ! This subroutine computes the intersection line between the free surface and the body.
    
    ierror = 0   
   
    ! Stockage a part des points appartenant a la meme courbe d'intersection 
    if(n_tab>0 .and. .not.allocated(domaine2%liste))then
        nPts = 1000*NBodies
        allocate(domaine2%liste(nPts))
    endif     

    
    do j=1,n_tab ! Number of intersection curves
        
        ! Getting back the intersection points from tab.
        ndim = 1
        P = tab(j)%pt%val ! First point of the intersection curve

        domaine2%liste(ndim) = P
        ptr0 => tab(j)%pt
        ptr1 => ptr0
        ptr2 => ptr0%suiv
        
        do while(associated(ptr2) .and. .not.associated(ptr2,ptr0)) ! Following points.
            P = ptr2%val
            ndim = ndim+1
            domaine2%liste(ndim) = P ! Intersection points
                        
            if(.not.associated(ptr1,ptr2%suiv))then
                ptr1 => ptr2
                ptr2 => ptr2%suiv
            elseif(.not.associated(ptr1,ptr2%prec))then
                ptr1 => ptr2
                ptr2 => ptr2%prec
            else
                ierror = 103
                goto 9999
            endif
        enddo
        
        ! Closing the intersection curve or not.
        is_close = .false.
        if(domaine2%liste(1)%bf/=2 .and. domaine2%liste(ndim)%bf/=2)then
            ndim=ndim+1
            domaine2%liste(ndim) = domaine2%liste(1)    
            is_close = .true.
        endif
        
        ! Number of points for the intersection curve.
        domaine2%dim = ndim
        if(ndim.eq.1)then
            print*,""
            print*,"Error Mesh_inter: there is only one point on the free surface. This configuration is not good for the mes generator. It would be better not to be in this configuration."
            pause
        end if
        
        ! Mesh of the intersection curve (polyline)
        np0 = nb_point 
        na0 = nb_arete
        call mesh_polyline(domaine2%liste,domaine2%dim,is_close,2.*dx,mesh,nb_point,nb_arete,10,1,ierror) ! 10 = bf, 2.*dx can be changed !!!
        
        if(ierror/=0)then
            ierror = 102
            goto 9999
        endif
     
        ! Orientation des aretes pour l'avance de front     
        call orient_hole_mesh(mesh,nb_point,nb_arete,np0,na0,ierror)
        if(ierror/=0)then
            ierror = 202
            goto 9999
        endif
     
    enddo 
    
    ! Deallocating
    if(allocated(domaine2%liste)) deallocate(domaine2%liste)
    
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error#',i3," : maillage courbe d'intersection impossible.")

  end subroutine mesh_inter 
  
! *********************************************************************
!  Mesh 1D 
! *********************************************************************
  subroutine mesh1D(mesh,nb_point,nb_arete,fgeom,igeom,dx,t,imode,ierror,InputData,NumBody,tab,n_tab,n_tab_lines)
    
    !f2py integer*1, dimension(1000)                    :: mesh
    type(MGrid),intent(inout)                           :: mesh                                             ! Mgrid.
    integer,intent(inout)                               :: nb_point,nb_arete                                ! Number of points and edges.
    !f2py integer*1, dimension(1000)                    :: fgeom
    type(type_geom),intent(in)                          :: fgeom                                            ! Geometry.
    real(rp),intent(in)                                 :: dx                                               ! dx = dx2 for the floaters, dx1 for the domain.
    real(rp),intent(in)                                 :: t                                                ! Current time.
    integer,intent(in)                                  :: imode, igeom                                     ! igeom = 0 for the domain, 1 for the floaters.
    integer,intent(inout)                               :: ierror                                           ! Error flag.
    !f2py integer*1, dimension(1000)                    :: InputData
    type(InputDataStruct),intent(in)                    :: InputData                                        ! Input data.
    integer                                             :: NumBody                                          ! Body number.
    !f2py integer*1, dimension(1000)                    :: tab
    type(chaine_point_pt),dimension(:),optional         :: tab                                              ! Table of the pointers toward the intersection points.
    integer,intent(in),optional                         :: n_tab,n_tab_lines                                ! Number of intersection curves and lines.
    
    integer                                             :: j,k,k1                                           ! Loop parameters.
    integer                                             :: narete,ncom,nface,naux,ncercle,N,ncom0,ncom_sym  !
    integer                                             :: isens                                            ! Sens of the path of the nodes (actually isens is fixed to 1).
    !f2py integer*1, dimension(1000)                    :: arete
    type(GArete)                                        :: arete                                            ! Edge.
    !f2py integer*1, dimension(1000)                    :: cercle
    type(GCercle)                                       :: cercle                                           ! Circle.
    !f2py integer*1, dimension(1000)                    :: polyline
    type(GPolyline)                                     :: polyline                                         ! Polyline.
    !f2py integer*1, dimension(1000)                    :: domaine2
    type(domaine)                                       :: domaine2                                         ! Domain.
    integer,dimension(2)                                :: connectj                                         !
    integer,dimension(nfvmax)                           :: aux2, aux2_sym                                   !
    integer,dimension(1000)                             :: aux                                              !
    integer                                             :: naux_HorizontalCyl                               !
    integer,dimension(1000)                             :: aux_HorizontalCyl                                !
    !f2py integer*1, dimension(1000)                    :: P,Q,Pcart,Ploc
    type(point)                                         :: P,Q,Pcart,Ploc                                   ! Points.
    real(rp)                                            :: nom,denom,s                                      !
    !f2py integer*1, dimension(1000)                    :: vertex
    type(MVertex)                                       :: vertex                                           ! Mvertex.
    !f2py integer*1, dimension(1000)                    :: pile1
    type(pile_element),pointer                          :: pile1                                            !
    real(rp),dimension(3)                               :: y1,y2,yl                                         !
    real(rp),dimension(3)                               :: M                                                ! Vector.
    real(rp)                                            :: angle,dtheta,L,eta, dxl, dxl4                    !
    integer                                             :: ind1,ind2,na0,np0                                !
    !f2py integer*1, dimension(1000)                    :: tab2
    type(point),dimension(10)                           :: tab2                                             !
    !f2py integer*1, dimension(1000)                    :: tab3
    type(MVertex),dimension(10)                         :: tab3                                             !
    integer                                             :: ntab2,nsurface,nb_point0,ntab3,nline             !
    logical                                             :: HorizontalCyl                                    ! Horizontal cylinder (true) or not (false). Necessary to find the intersection points in common between the circles and the body of the cylinder. BRICOLAGE !!!
    
    ! This subroutine creates the mesh1D (segments, polylines, etc.).
    
    ierror = 0
   
    nsurface = fgeom%nplan + fgeom%ndisque + fgeom%ncylindre + fgeom%ncone + fgeom%nsphere  + fgeom%naxisym
    nline = fgeom%ncercle + fgeom%narete + fgeom%npolyline
    
    ! --------------------------------------------------------------------
    ! Mesh contour surface libre seule   
    ! --------------------------------------------------------------------
    
    if(nsurface==1 .and. present(tab) .and. n_tab.ne.0 .and. nline==0)then ! Only for the free surface when nsurface = 1 (rarely).
                
        ntab2 = 0
        nb_point0 = nb_point
        if(fgeom%nplan==1)then ! Plan
            do k=1,nb_point
                if(mesh%point(k)%nface==1 .and. mesh%point(k)%face(1)==HouleRF%index)then
                    ntab2 = ntab2+1
                    vertex = mesh%point(k) ! Operator overloading problem ?
                    vertex%coord(3) = 0._rp
                    nface = vertex%nface
                    tab2(ntab2) = vertex%coord
                    tab2(ntab2)%nface = nface
                    tab2(ntab2)%face(1:nface) = vertex%face(1:nface)
                    tab2(ntab2)%bf = vertex%bf
                    tab2(ntab2)%flag = k
                endif
            enddo
            ntab2 = ntab2+1
            tab2(ntab2) = tab2(1)
            if(ntab2>=2)then
                call mesh_polyline(tab2,ntab2,.false.,dx1,mesh,nb_point,nb_arete,1,0,ierror) ! 1 = bf
            else
                ierror = 104
                goto 9999
            endif
        elseif(fgeom%ndisque==1)then ! Disc
            call mesh_cercle(mesh,nb_point,nb_arete,dx1,fgeom%cercle(1),&
            &   [fgeom%disque(1)%index,0],1,1,-1,ierror)
        endif
    
        do k=nb_point0+1,nb_point
            P = mesh%point(k)%coord
            call CEta0(P%coord,t,eta)
            mesh%point(k)%coord(3) = eta
        enddo
    
    elseif(nsurface>1 .or. nline>0)then
        
        ! ---------------------------------------------------------------------
        ! Mesh des aretes immergées
        ! ---------------------------------------------------------------------
        nullify(pile1)
        narete = fgeom%narete
        do j=1,narete
               
            arete = fgeom%arete(j)
            connectj = arete%face 
            if(imode==1)then
                call common_int(connectj,2,[HouleRF%index],1,aux2,ncom,ierror)
                if(ncom==1) goto 1000
            endif
            call common_int(connectj,2,[HouleRF%index,-1],2,aux2_sym,ncom_sym,ierror)
            
            if(is_body)then
                dxl = InputData%dx2(NumBody)
            else
                dxl = dx2Domain
            end if
                        
            ! Recherche des points appartenant à l'arête.
            k1 = 0
            do k = 1,nb_point
                if(mesh%point(k)%nedge>0)then                    
                    call common_int(mesh%point(k)%edge,mesh%point(k)%nedge,[arete%index],1,aux2,ncom,ierror)
                    if(ncom==1)then
                        k1 = k1 + 1   ! Number of points in common (mesh and edge from the geometry)
                        aux(k1) = k ! Reference of points in common
                    end if
                end if
            end do  
                        
            naux = k1
                        
            ! Repositionnement des points selon abscisse.
            denom = norm2(arete%P2%coord - arete%P1%coord) ! Length of the edge
            do k = 1,naux
                Pcart = mesh%point(aux(k)) ! Point in common in cartesian
                call cart2loc(Pcart%coord,fgeom%repere,Ploc%coord) ! Point in common in local coordinates           
                nom = norm2(Ploc%coord - arete%P1%coord)
                s = nom/denom
                call add_element_pile(pile1,aux(k),s,k)
            end do
                        
            ! Maillage des segments. 
            ntab3 = 0
            do while (associated(pile1))
                ntab3 = ntab3+1
                if(ntab3.gt.10)then
                    print*,"Mesh1D: ntab3 is greather than 10 for the edge: ",j
                end if
                tab3(ntab3) = mesh%point(pile1%val)
                call depile(pile1,ierror)
            enddo
            nb_point0 = nb_point
            
            if(ntab3>=2)then
                if(ntab3==4)then
                    call mesh_segment(tab3(1),tab3(2),dx1,dxl,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                    call mesh_segment(tab3(3),tab3(4),dxl,dx1,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                elseif(ntab3==3)then
                    call mesh_segment(tab3(1),tab3(2),dx1,dxl,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                    call mesh_segment(tab3(2),tab3(3),dxl,dx1,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                elseif(ntab3 == 6)then
                    call mesh_segment(tab3(1),tab3(2),dx1,dxl,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                    call mesh_segment(tab3(3),tab3(4),dxl,10*dxl,mesh,nb_point,nb_arete,1,connectj,2,1) ! dx1 != dxl, 10 to have less point between the two floaters
                    print*,"Between two floaters, the panel size is 10*dx2."
                    call mesh_segment(tab3(5),tab3(6),dxl,dx1,mesh,nb_point,nb_arete,1,connectj,2,0) ! dx1 != dxl
                elseif(ncom_sym==2)then
                    call mesh_segment_sym(tab3(1),tab3(2),dx,dxl,mesh,nb_point,nb_arete,1,connectj,2) ! dx1 != dxl
                else
                    call mesh_segment(tab3(1),tab3(2),dx,dx,mesh,nb_point,nb_arete,1,connectj,2,0)
                endif
            endif
            
            call common_int(connectj,2,[HouleRF%index],1,aux2,ncom,ierror)
            if(ncom==1)then
                do k=nb_point0+1,nb_point
                P = mesh%point(k)%coord
                call CEta0(P%coord,t,eta)
                mesh%point(k)%coord(3) = eta
                enddo
            endif

            1000 continue    
 
        end do
        
        ! -------------------------------------------------------------------------
        ! Mesh des cercles immergées ou de la surface libre
        ! -------------------------------------------------------------------------
        
        ncercle = fgeom%ncercle
        
        if(igeom==0)then
            dxl4 = dx4
        else
            dxl4 = dx3
        endif
        
        do j = 1,ncercle
            
            HorizontalCyl = .false.
            cercle = fgeom%cercle(j)
            connectj = cercle%face 
            call common_int(connectj,2,[HouleRF%index],1,aux2,ncom0,ierror)
            if(imode==1 .and. ncom0==1) goto 2000
            np0 = nb_point
                        
            ! Recherche des points appartenant à l'arête     
            k1 = 0
            do k = 1,nb_point                
                if(mesh%point(k)%nedge>0)then                    
                    call common_int(mesh%point(k)%edge,mesh%point(k)%nedge,[cercle%index],1,aux2,ncom,ierror)
                    if(ncom==1)then
                        k1 = k1 + 1
                        aux(k1) = k
                    end if
                end if   
            end do
            naux = k1
            
            ! Horizontal cylinders - Only for the floaters - BRICOLAGE !!!
            if(igeom.eq.1 .and. fgeom%ncylindre.ne.0)then
                if(present(n_tab_lines))then
                    if(n_tab_lines.ge.(3*fgeom%ncylindre))then ! If the cylinder is horizontal, there are several parts in the intersection curves (at least 3, only 1 for a vertical cylinder) - BRICOLAGE !!!
                        ! Horizontal cylinder on the (x0z) plan.
                        print*,"Mesh1D: Modification for the horizontal cylinders."
                        k1 = 0
                        ! Searching the endpoints of the circles in connection with the body of the cylinder on the intersection curve. Why? Because the intersection algorithm does not suceed to do it itself. I don't know why (PYW).
                        ! Thus naux = 2 and in aux there are the two endpoints of the circles.
                        do k = 1,nb_point
                            if(mesh%point(k)%nface>0)then              
                                call common_int(mesh%point(k)%face,mesh%point(k)%nface,[fgeom%disque(j)%index],1,aux2,ncom,ierror) ! Disque%index and not cercle%index, %face and not %edge. Problem to get the end points back. BRICOLAGE !!!
                                if(ncom==1)then
                                    if(k1.eq.0)then
                                        k1 = 1
                                        aux_HorizontalCyl(1) = k
                                        aux_HorizontalCyl(2) = k
                                    else
                                        ! Only work in the cylinder is in the plan (x0z).
                                        ! Could be done in local coordinates to be more consistent.
                                        if(mesh%point(k)%coord(2) .le. mesh%point(aux_HorizontalCyl(1))%coord(2))then ! End points with min y.
                                            aux_HorizontalCyl(1) = k
                                            k1 = 2
                                        end if
                                        if(mesh%point(k)%coord(2) .ge. mesh%point(aux_HorizontalCyl(2))%coord(2))then ! End points with max y.
                                            aux_HorizontalCyl(2) = k
                                            k1 = 2
                                        end if
                                    end if
                                end if
                            end if   
                        end do
                        naux_HorizontalCyl = 2
                        ! In case there is no such point, that means the cylinder is considered as vertical.
                        if(aux_HorizontalCyl(1).ne.0 .and. aux_HorizontalCyl(2).ne.0)then 
                            if(abs(mesh%point(aux_HorizontalCyl(1))%coord(2)).gt.Epsilon .and. abs(mesh%point(aux_HorizontalCyl(2))%coord(2)).gt.Epsilon)then ! If mesh%point(aux_HorizontalCyl(1:2))%coord = (0,0,0) : no intersection -> full circle.
                                naux = naux_HorizontalCyl
                                aux(1) = aux_HorizontalCyl(1)
                                aux(2) = aux_HorizontalCyl(2)
                                HorizontalCyl = .true.
                                print*,"Mesh1D: aux and naux modified for the horizontal cylinders."
                            end if
                        end if
                    end if
                end if
            end if
                        
            ! Mesh des arcs de cercle
            do k = 1,naux-1,2 ! From 1 to naux-1 with a step of 2.
                ind1 = aux(k)
                ind2 = aux(k+1)
                isens = 1                
                P = mesh%point(ind1)%coord ! Only the coordinates are initialized, P and Q are a point structure.
                Q = mesh%point(ind2)%coord
                call cart2sph(P%coord,cercle%repere,y1)
                call cart2sph(Q%coord,cercle%repere,y2)   
                angle = y2(2)-y1(2)
                if(abs(angle).lt.Epsilon)then
                    angle = y2(3)-y1(3)
                end if
                L = abs(angle)*cercle%rayon
                N = int(L/dxl4)
                if(N.ne.0)then
                    dtheta = angle/dble(N)
                else
                    ! The points are too close.
                    print*,""
                    print*,"naux = ",naux
                    if(present(n_tab_lines))then
                        print*,"n_tab_lines = ",n_tab_lines
                    end if
                    print*,"P = ",P%coord
                    print*,"Q = ",Q%coord
                    print*,"y1 = ",y1
                    print*,"y2 = ",y2
                    print*,"angle = ",angle
                    print*,"L = ",L
                    print*,"dxl4 = ",dxl4
                    print*,"N = ",N
                    print*,"mesh1D: N is equal to 0."
                end if
                          
                yl = y1
                yl(2) = y1(2) + dtheta                
                call sph2cart(yl,cercle%repere,M)
                call CEta0(M,t,eta)
                
                call mesh_arc_cercle(mesh%point(ind1),mesh%point(ind2),dxl4,mesh,nb_point,&
                &                    nb_arete,-1,cercle,connectj,2,isens,ierror,igeom,HorizontalCyl) ! -1 = bf
                
                if(ierror /= 0)then
                    go to 9999
                end if
                
            end do
                        
            ! Sinon maillage cercle entier si le centre de celui-ci se trouve sous la surface libre ou si c'est la surface libre.
            if(naux==0)then
                np0 = nb_point+1
                na0 = nb_arete+1
                P = cercle%repere%origine
                
                ! In case of the tank, CEta0 cannot be used because if the phase is not good, the wave elevation is below 0 at P. eta must igeom be above 0 in case of wave.
                call CEta0(P%coord,t,eta)
                                
                ! Mesh of the circle.
                if(P%coord(3).le.eta+Epsilon .or. igeom.eq.0)then ! Verification if the circle is above or below the free surface. In case of igeom = 0 (tank), all the circles must be mesh.
                    call mesh_cercle(mesh,nb_point,nb_arete,dxl4,cercle,connectj,2,1,1,ierror)
                end if
                
            end if
                        
            if(ncom0==1)then
                do k=np0+1,nb_point
                    M = mesh%point(k)%coord
                    call CEta0(M,t,mesh%point(k)%coord(3))
                end do
            end if
                        
        2000 continue     
     
        end do
    end if
    
    ! ----------------------------------------------------------------------
    ! Mesh of polylines 
    ! ----------------------------------------------------------------------
    
    nullify(pile1)
    n = fgeom%npolyline
    
    do j = 1,n
                
        polyline = fgeom%polyline(j)
        if(is_body)then
            !dxl = InputData%dx2(NumBody)
            dxl = 0.8*InputData%dx2(NumBody)
            print*,"Mesh1D: dxl = 0.8*dx2."
        else
            dxl = dx2Domain
        end if
                
        ! Find point on the polyline.
        k1 = 0
        do k = 1,nb_point
            if(mesh%point(j)%nedge>0)then
                call common_int(mesh%point(k)%edge,mesh%point(k)%nedge,[polyline%index],1,aux2,ncom)
                if(ncom==1)then
                    k1 = k1+1
                    aux(k1) = k
                end if
            end if
        end do
        naux = k1
                        
        ! Positionning node along z coordinate.
        do k = 1,naux
            P = mesh%point(aux(k))
		    
            call cart2loc(P%coord,polyline%repere,P%coord) ! It is not a good thing to have the same structure as both input and output.
            s = P%coord(2)
            
            call add_element_pile(pile1,aux(k),s,k)
            
        end do
                
        ntab3 = 0  
        do while (associated(pile1))
            ntab3 = ntab3 + 1
            if(ntab3.gt.10)then
                print*,"Mesh1D: ntab3 is greather than 10 for the polyline: ",j
            end if
            tab3(ntab3) = mesh%point(pile1%val)
            call depile(pile1,ierror)
        end do
                
        if(ntab3==2)then
            call mesh_polyline2(tab3(1),tab3(2),polyline%n,polyline%P,polyline%repere,dxl,&
    &                               mesh,nb_point,nb_arete,1,polyline%index,ierror) ! bf  = 1            
        else if(ntab3 == 0)then
            ! No intersection curve.
            exit
        else if(ntab3 == 1)then            
            ! Only one point of the intersection curve.
            exit
        else
            print*,'Mesh1D, polylines: case not implemented.'
            print*,"ntab3 = ",ntab3
            pause
            ierror = 100
            goto 9999
        end if
		
    end do
    
    ! ----------------------------------------------------------------------   
    ! Sortie Maillage 1D
    ! ----------------------------------------------------------------------- 

    9999 continue

    if(allocated(domaine2%liste)) deallocate(domaine2%liste)

    if(ierror /= 0)then 
        write(*,99),ierror
    endif  
    99 format('** error #',i3,': mesh1d')        
       
end subroutine mesh1D       

! ***********************************************************************
! Maillage 2D
! ***********************************************************************
subroutine mesh2d(mesh,nb_point,nb_arete,nb_tri,fgeom,t,dx,ierror,InputData,tin,ioMeshevol,CartesianGrid)
    
    !f2py integer*1, dimension(1000)        :: mesh
    type(MGrid),intent(inout)               :: mesh                                                             ! MGrid structure.
    integer,intent(inout)                   :: nb_point,nb_arete,nb_tri                                         ! Number of points, edges and triangles.
    !f2py integer*1, dimension(1000)        :: fgeom
    type(type_geom),intent(in)              :: fgeom                                                            ! Geometry.
    real(rp),intent(in)                     :: t                                                                ! Current time.
    real(rp),intent(in)                     :: dx                                                               ! Panel discretization.
    integer,intent(inout)                   :: ierror                                                           ! Error flag.
    !f2py integer*1, dimension(1000)        :: InputData
    type(InputDataStruct),intent(in)        :: InputData                                                        ! Input data.
    real(rp),intent(inout)                  :: tin                                                              ! Plotting time for Advance_front.dat.
    integer,intent(in)                      :: ioMeshevol                                                       ! Plotting or not Advance_front.dat or Advance_front_Remesh.dat.
    logical,optional                        :: CartesianGrid                                                    ! Plotting CartesianGrid.dat if 1st first of the domain.
        
    !f2py integer*1, dimension(1000)        :: mesh2
    type(MGrid)                             :: mesh2                                                            ! Temporary MGrid structure of the geometry.
    integer                                 :: nb_point2,nb_arete2,nb_tri2                                      ! Number of points, edges and triangles of the mesh of the geometry.
    !f2py integer*1, dimension(1000)        :: cyl
    type(cylindre_2)                        :: cyl                                                              ! Cylinder.
    !f2py integer*1, dimension(1000)        :: disque1
    type(disque2)                           :: disque1                                                          ! Disc.
    !f2py integer*1, dimension(1000)        :: plan
    type(Gplan)                             :: plan                                                             ! Plan.
    !f2py integer*1, dimension(1000)        :: cone
    type(cone)                              :: cone                                                             ! Cone.
    !f2py integer*1, dimension(1000)        :: sphere
    type(sphere)                            :: sphere                                                           ! Sphere.
    !f2py integer*1, dimension(1000)        :: axisym
    type(axisym)                            :: axisym                                                           ! Axisym geoemtry.
    !f2py integer*1, dimension(1000)        :: wigley
    type(TWigley)                           :: wigley                                                           ! Wigley hull.
    real(rp)                                :: x1,x2,d,poids,r,periode,eta,inv_r                                !
    real(rp)                                :: xmin,xmax,ymin,ymax,rd                                           ! Parameters for the Cartesian grid.
    integer                                 :: j,jj,k                                                           ! Loop parameters.
    integer                                 :: n1,n2,i1,i2,i5,i6,ind,ip                                         !
    integer                                 :: iface,bf                                                         ! Number of the face and boundary flag.
    integer,dimension(:),allocatable        :: ind_pt,ind_arr
    !f2py integer*1, dimension(1000)        :: pile1,ptr1,ptr2
    type(pile_element),pointer              :: pile1,ptr1,ptr2                                                  ! Pointers.
    !f2py integer*1, dimension(1000)        :: vertex_11,vertex_12,vertex_21,vertex_22,vertex_p21,vertex_p11
    type(MVertex)                           :: vertex_11,vertex_12,vertex_21,vertex_22,vertex_p21,vertex_p11    !    
    !f2py integer*1, dimension(1000)        :: edge
    type(MEdge)                             :: edge                                                             ! Edge.
    real(rp),dimension(3)                   :: M1,M2,M                                                          ! 3D vectors.
    !f2py integer*1, dimension(1000)        :: V
    type(vector)                            :: V                                                                ! Vector.
    !f2py integer*1, dimension(1000)        :: P
    type(point)                             :: P                                                                ! Points.
    !f2py integer*1, dimension(1000)        :: rep,prej
    type(repere3d)                          :: rep,repj                                                         ! Frames.
    logical                                 :: bool                                                             ! Boolean.
    real(rp),dimension(:,:,:),allocatable   :: matdref                                                          ! Reference length for each cell of the Cartesian grid.
    real(rp),dimension(:),allocatable       :: xgrid,ygrid                                                      ! Cartesian grid.
    integer                                 :: nx,ny,npoint                                                     ! Number of points in the Cartesian grid.
    real(rp)                                :: rmax,a,a2                                                        !
    real(rp),dimension(:),allocatable       :: param                                                            ! 8 param for a cone, 2*axisym%npoint+1 for a axisym, 1 for a cylinder, a disc, a sphere, a plan and a Wigley hull.
         
    ! This subroutine creates the mesh2D (cylinders, surfaces, etc.).
    
    ierror = 0
    allocate(ind_pt(40000))
    allocate(ind_arr(80000))
    
    ! ------------------------------------------------------------------------
    ! Maillage des surfaces cylindriques
    ! -------------------------------------------------------------------------
    
    n2 = fgeom%ncylindre
    do j = 1,n2
        if(fgeom%cmd_cyl(j) == 1)then
                                
            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)
            cyl = fgeom%cylindre(j)
            iface = cyl%index
            bf = -9
            
            ! In case of a cylindrical tank (idtype.eq.2 .or. idtype.eq.4) in closed domain (cuve_ferme), the external boundaries (iface=3) are meshed, not otherwise.
            if(iface.ne.3 .or. ((idtype.eq.2 .or. idtype.eq.4) .and. iface.eq.3 .and. cuve_ferme))then
                
                if(idebug>0)then
                    print*,'Maillage surface cylindrique : ',iface
                endif
            
                ! Recuperation des points appartenant à la surface.
                call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
                &    ind_pt,ind_arr,iface,ierror)
                        
                call local_coord_cyl(mesh2,nb_point2,nb_arete2,nb_tri2,cyl,ierror)
            
                ! Maillage frontière virtuelle périodique dans l'espace des paramètres.
                k=0
                nullify(pile1)
                do jj = 1,nb_arete2
                    i1 = mesh2%arrete(jj)%iP(1)
                    i2 = mesh2%arrete(jj)%iP(2)
                    x1 = mesh2%point(i1)%coord(1)
                    x2 = mesh2%point(i2)%coord(1)
                    d = x2-x1
                    if(abs(d) .ge. PI*cyl%rayon .and. not(Symmetry))then
                        k = k + 1
                        poids = mesh2%point(i1)%coord(2)
                        if(x1<x2)then
                            call add_element_pile(pile1,i1,poids,i2)
                        else
                            call add_element_pile(pile1,i2,poids,i1)
                        endif
                    endif
                end do
            
                do while(associated(pile1))
                    ptr1 => pile1
                    if(associated(pile1%suiv)) then
                        ptr2 => pile1%suiv
                    else
                        goto 100
                    endif  
                    vertex_11 = mesh2%point(ptr1%val)
                    vertex_12 = mesh2%point(ptr1%bf)
                    vertex_21 = mesh2%point(ptr2%val)
                    vertex_22 = mesh2%point(ptr2%bf)
    
                    r = cyl%rayon
                    periode = 2.*PI*r                
                    call mesh_periodic_boundary(vertex_21,vertex_11,dx,dx,mesh2,nb_point2,nb_arete2,&
                    &                             ind_pt,ind_arr,0._rp,periode,i5,i6,iface)
      
                    vertex_p21 = mesh2%point(i5)
                    vertex_p11 = mesh2%point(i6)
      
                    call get_index(mesh2,nb_arete2,vertex_11,vertex_12,ind)
                    call change_edge(mesh2,mesh2%arrete(ind),&
                    &                vertex_12,vertex_p11,bf)
     
                    call get_index(mesh2,nb_arete2,vertex_21,vertex_22,ind)
                    call change_edge(mesh2,mesh2%arrete(ind),&
                    &                vertex_p21,vertex_22,bf)
      
                100 continue
    
                    call depile_pile_element(pile1,ierror)
                    call depile_pile_element(pile1,ierror)
         
                end do    
                        
                ! Orientation des aretes
                ymin = cyl%vmin
                ymax = cyl%vmax
                xmin = -PI*cyl%rayon
                if(Symmetry) xmin = 0._rp
                xmax = PI*cyl%rayon
                inv_r = 1.d0/cyl%rayon
            
                do jj = 1,nb_arete2
                    edge = mesh2%arrete(jj)
                    M1 = mesh2%point(edge%iP(1))%coord
                    M2 = mesh2%point(edge%iP(2))%coord
                    M = 0.5_RP*(M1+M2)
                    call assign_vector_coord(V,[M1(2)-M2(2),M2(1)-M1(1),0._RP])
                    P = M + 0.4_RP*V%coord
                    call Computation_eq_cyl3D(P%coord(1)*inv_r,P%coord(2),cyl,P%coord)
                    call CEta0(P%coord,t,eta)
      
                    bool=.false.
      
                    if(abs(M1(2)-ymin).lt.Epsilon .and. abs(M2(2)-ymin).lt.Epsilon .and. M2(1) < M1(1))then
                        bool = .true.
                    elseif(abs(M1(2)-ymax).lt.Epsilon .and. abs(M2(2)-ymax).lt.Epsilon .and. M1(1) < M2(1))then
                        bool = .true.
                    elseif(abs(M1(1)-xmin).lt.Epsilon .and. abs(M2(1)-xmin).lt.Epsilon .and. M2(2)>M1(2))then
                        bool = .true.
                    elseif(abs(M1(1)-xmax).lt.Epsilon .and. abs(M2(1)-xmax).lt.Epsilon .and. M2(2)<M1(2))then
                        bool = .true.
                    elseif(P%coord(3) > eta .and. edge%bf /= -12 .and. edge%bf /= -11)then
                        bool = .true.  
                    endif
                    if(bool)then
                        ip = mesh2%arrete(jj)%iP(1)
                        mesh2%arrete(jj)%iP(1) = mesh2%arrete(jj)%iP(2)
                        mesh2%arrete(jj)%iP(2) = ip
                        mesh2%arrete(jj)%dir = -mesh2%arrete(jj)%dir
                    end if
                end do   
                
                ! Creation grille des tailles de reference
                if(nb_point2.gt.0)then
                    if(xmax-xmin .le. 2.*hrefx)then
                        nx = 2
                    else
                        nx = int((xmax-xmin)/hrefx)
                    endif
                        if(ymax-ymin .le. 2.*hrefy)then
                        ny = 2
                    else
                        ny = int((ymax-ymin)/hrefy)
                    endif
                    
                    ! Protection contre nombre de grille null
                    allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                    
                    call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                    &    matdref,nx,ny,iface,idref,dx,ierror,InputData)  
                    
                    if(ierror/=0)then
                        ierror = 203
                        goto 9999
                    endif  
                endif
                
                ! Maillage 2D de la surface
                if(allocated(xgrid) .and. allocated(ygrid))then
                    call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,iface,&
                    &                     0,[1._rp],1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                    if(ierror/=0)then
                        ierror = 101
                        goto 9999
                    endif
                end if
                
                inv_r = 1.d0/cyl%rayon
                do k=1,nb_point2
                    M = mesh2%point(k)%coord
                    call Computation_eq_cyl3D(M(1)*inv_r,M(2),cyl,mesh2%point(k)%coord)
                enddo
   
                !   Mise à jour maillage   
                call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
                &               ind_pt,ind_arr,1,iface,ierror)
                if(ierror/=0)then
                    ierror = 102
                    goto 9999
                endif
    
                if(allocated(xgrid)) deallocate(xgrid)
                if(allocated(ygrid)) deallocate(ygrid)
                if(allocated(matdref)) deallocate(matdref)
            end if
        end if   
    end do
    
    ! -------------------------------------------------------------------------
    ! Maillage des cones
    ! -------------------------------------------------------------------------
    
    n1 = n2
    n2 = fgeom%ncone
    
    do j=1,n2
   
        if(fgeom%cmd_cone(j)==1)then

            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)
            cone = fgeom%cone(j)
            iface = cone%index
            bf = -9.
  
            if(idebug>0)then
                print*,'Maillage surface conique : ',iface
            endif

            call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
            &        ind_pt,ind_arr,iface,ierror)

            call local_coord_cone(mesh2,nb_point2,nb_arete2,nb_tri2,cone,ierror)

            a = (cone%r2-cone%r1)/cone%long

            k=0
            nullify(pile1)
            do jj=1,nb_arete2
                i1 = mesh2%arrete(jj)%iP(1)
                i2 = mesh2%arrete(jj)%iP(2)
                x1 = mesh2%point(i1)%coord(1)
                x2 = mesh2%point(i2)%coord(1)
                d = x2-x1
                if(abs(d).ge.PI*cone%r1 .and. not(Symmetry))then
                    k=k+1
                    poids = mesh2%point(i1)%coord(2)
                    if(x1<x2)then
                        call add_element_pile(pile1,i1,poids,i2)
                    else
                        call add_element_pile(pile1,i2,poids,i1)
                    endif
                endif
            enddo
    
            do while(associated(pile1))
                ptr1 => pile1
                if(associated(pile1%suiv))then
                    ptr2 => pile1%suiv
                else
                    goto 110
                endif
                vertex_11 = mesh2%point(ptr1%val)
                vertex_12 = mesh2%point(ptr1%bf)
                vertex_21 = mesh2%point(ptr2%val)
                vertex_22 = mesh2%point(ptr2%bf)

                periode = 2.*PI*cone%r1
                a2 = 2.*PI*a
                ! FIME : rd
                rd = sqrt((PI**2+1)*(cone%r2-cone%r1)**2+cone%long**2)
                rd = rd / sqrt(PI**2*(cone%r2-cone%r1)**2+cone%long**2)
                call mesh_periodic_boundary(vertex_21,vertex_11,dx/rd,dx/rd,mesh2,nb_point2,nb_arete2,&
            &                                   ind_pt,ind_arr,a2,periode,i5,i6,iface)

                vertex_p21 = mesh2%point(i5)
                vertex_p11 = mesh2%point(i6)

                call get_index(mesh2,nb_arete2,vertex_11,vertex_12,ind)
                call change_edge(mesh2,mesh2%arrete(ind),&
            &                        vertex_12,vertex_p11,bf)

                call get_index(mesh2,nb_arete2,vertex_21,vertex_22,ind)
                call change_edge(mesh2,mesh2%arrete(ind),&
            &                        vertex_p21,vertex_22,bf)

            110 continue

                call depile_pile_element(pile1,ierror)
                call depile_pile_element(pile1,ierror)
                enddo

                ! Orientation des aretes
                ymin = cone%vmin
                ymax = cone%vmax
                xmin = -PI*max(cone%r1,cone%r2)
                if(Symmetry) xmin = 0._rp
                xmax = abs(xmin)
                do jj=1,nb_arete2
                    edge = mesh2%arrete(jj)
                    M1 = mesh2%point(edge%iP(1))%coord
                    M2 = mesh2%point(edge%iP(2))%coord
                    M = 0.5_rp*(M1+M2)
                    call assign_vector_coord(V,[M1(2)-M2(2),M2(1)-M1(1),0._rp])
                    P = M + 0.4*V%coord
                    inv_r = 1._rp/(a*P%coord(2)+cone%r1)
                    call Computation_eq_cone3D(P%coord(1)*inv_r,P%coord(2),cone,P%coord)
                    call CEta0(P%coord,t,eta)

                    bool = .false.

                    if(abs(M1(2)-ymin).lt.Epsilon .and. abs(M2(2)-ymin).lt.Epsilon .and. M2(1) < M1(1))then
                        bool = .true.
                    elseif(abs(M1(2)-ymax).lt.Epsilon .and. abs(M2(2)-ymax).lt.Epsilon .and. M1(1) < M2(1))then
                        bool = .true.
                    elseif(P%coord(3) > eta .and. edge%bf /= -12 .and. edge%bf /= -11)then
                        bool = .true.  
                    endif
                    if(bool)then
                        ip = mesh2%arrete(jj)%iP(1)
                        mesh2%arrete(jj)%iP(1) = mesh2%arrete(jj)%iP(2)
                        mesh2%arrete(jj)%iP(2) = ip
                        mesh2%arrete(jj)%dir = -mesh2%arrete(jj)%dir
                    endif
                enddo

            ! Creation grille des tailles de reference
            if(nb_point2.gt.0)then
                if(xmax-xmin .le. 2.*hrefx)then
                    nx = 2
                else
                    nx = int((xmax-xmin)/hrefx)
                endif
                    if(ymax-ymin .le. 2.*hrefy)then
                    ny = 2
                else
                    ny = int((ymax-ymin)/hrefy)
                endif
                ! Protection contre nombre de grille null
                allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                &    matdref,nx,ny,iface,idref,dx,ierror,InputData)  
                if(ierror/=0)then
                    ierror = 203
                    goto 9999
                endif  
            endif    
    
            ! Maillage 2D de la surface   
            allocate(param(8))
            param(1) = (cone%r2-cone%r1)/cone%long
            param(2) = (cone%r1*cone%vmax-cone%r2*cone%vmin)/cone%long
            param(3) = (cone%r2-cone%r1)/sqrt((cone%r2-cone%r1)**2+cone%long)
            param(4) = (cone%vmax-cone%vmin)/sqrt((cone%r2-cone%r1)**2+cone%long)
            param(5) = cone%r1
            param(6) = cone%r2
            param(7) = cone%vmin
            param(8) = cone%vmax
            if(allocated(xgrid) .and. allocated(ygrid))then
                call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
                &                     iface,2,param(1:8),8,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                if(ierror/=0)then
                    ierror = 101
                    goto 9999
                endif
            end if
    

            do k=1,nb_point2
                M = mesh2%point(k)%coord
                inv_r = 1._rp/(a*M(2)+cone%r1)
                call Computation_eq_cone3D(M(1)*inv_r,M(2),cone,mesh2%point(k)%coord)
            enddo

            !   Mise à jour maillage   
            call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
            &               ind_pt,ind_arr,1,iface,ierror)
            if(ierror/=0)then
                ierror = 102
                goto 9999
            endif
    
            if(allocated(xgrid)) deallocate(xgrid)
            if(allocated(ygrid)) deallocate(ygrid)
            if(allocated(matdref)) deallocate(matdref)
            if(allocated(param)) deallocate(param)
      
        endif   
           
    enddo

    ! -------------------------------------------------------------------------  
    ! Maillage des disques
    ! -------------------------------------------------------------------------
    
    n1 = n2
    n2 = fgeom%ndisque
    
    do j=1,n2
                
        if(fgeom%cmd_dis(j)==1 .or.fgeom%cmd_dis(j)==3 )then

  			if(idebug>0)then
                print*,'Maillage disque :',fgeom%disque(j)%index
            endif
  
            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)  
            disque1 = fgeom%disque(j)
            iface = disque1%index
            rep%origine = disque1%repere%origine
            rep%e1 = disque1%repere%e1
            rep%e2 = disque1%repere%e2
            rep%e3 = disque1%repere%e3
            rmax = sqrt(disque1%r2max)
            xmin = -rmax ; xmax = rmax
            ymin = -rmax ; ymax = rmax
            if(Symmetry) ymin = 0._rp
            bf = -1

            ! Recuperation des points appartenant à la surface    
            call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
            &    ind_pt,ind_arr,iface,ierror)
    
            call local_coord_rep(mesh2,nb_point2,nb_arete2,nb_tri2,fgeom%disque(j)%repere,ierror)
            
            ! If there is no body piercing the free surface.
            if( (not(is_body) .and. (iface==1 .or. iface==2) .and. .not. Symmetry) .or. (not(Symmetry) .and. is_immerged .and. iface.eq.1) )then
                call init_cylindrical_mesh(mesh2,nb_point2,nb_arete2,nb_tri2,fgeom%disque(j),iface,dx2Domain,ierror)
            endif
       
            ! Orientation des aretes   
            if(iface/=1)then
                call orient_hole_mesh2(mesh2,nb_point2,nb_arete2,0,0,ierror)
            endif
 
            do jj=1,nb_arete2
                edge = mesh2%arrete(jj)
                M1 = mesh2%point(edge%iP(1))%coord
                M2 = mesh2%point(edge%iP(2))%coord
                M = 0.5*(M1+M2)
                call assign_vector_coord(V,[M1(2)-M2(2),M2(1)-M1(1),0.d0])
                P = M+0.4*V%coord
                call loc2cart(P%coord,fgeom%disque(j)%repere,P%coord)
                call CEta0(P%coord,t,eta)
      
                bool = .false.
                if(P%coord(3)>eta .and. edge%bf/=-12 .and. edge%bf/=-11 .and. iface.gt.1)then
                    bool = .true.
                endif

                if(Symmetry)then
                    if(abs(M1(2)-ymin).lt.Epsilon .and. abs(M2(2)-ymin).lt.Epsilon .and. M1(1)>M2(1))then
                        bool = .true.
                    endif
                endif
                if (bool) then
                    ip = mesh2%arrete(jj)%iP(1)
                    mesh2%arrete(jj)%iP(1) = mesh2%arrete(jj)%iP(2)
                    mesh2%arrete(jj)%iP(2) = ip
                    mesh2%arrete(jj)%dir = -mesh2%arrete(jj)%dir
                endif
            enddo     

            ! Creation grille des tailles de reference
            if(nb_point2.gt.0)then
                if(xmax-xmin.le.2.*hrefx)then
                    nx = 2
                else
                    nx = int((xmax-xmin)/hrefx)
                end if
                if(ymax-ymin.le.2.*hrefy)then
                    ny = 2
                else
                    ny = int((ymax-ymin)/hrefy)
                end if
                
                ! Modification of the free surface in case of a cylindrical tank.
                if(iface.eq.1 .and. (idtype.eq.2 .or. idtype.eq.4))then
                    if(xmax-xmin.le.2.*hrefxFS)then
                        nx = 2
                    else
                        nx = int((xmax-xmin)/hrefxFS)
                    end if
                    if(ymax-ymin.le.2.*hrefyFS)then
                        ny = 2
                    else
                        ny = int((ymax-ymin)/hrefyFS)
                    end if
                end if
                
                ! Protection contre nombre de maille null
                allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                &    matdref,nx,ny,iface,idref,dx,ierror,InputData)
                
                ! Cartesian_grid.dat
                ! If the free surface is rectangular, should be copied elsewhere.
                if(present(CartesianGrid) .and. iface.eq.1 .and. (idtype.eq.2 .or. idtype.eq.4))then ! iface = 1 for the FS, can be different if the free surface is in a rectangular domain.
                    call write_Cartesian_grid(xgrid,ygrid,nx,ny,matdref,t,iface) ! Cartesian_grid.dat
                end if
                
                if(ierror/=0)then
                    ierror = 203
                    goto 9999
                endif
            endif  

    
            ! Maillage de la surface 2D 
            if(allocated(xgrid) .and. allocated(ygrid))then

                call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
                &                     iface,0,[1._rp],1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                
                if(ierror/=0)then
                    ierror = 101
                    goto 9999
                endif
            end if
    
            do k=1,nb_point2
                M = mesh2%point(k)%coord
                call loc2cart(M,fgeom%disque(j)%repere,mesh2%point(k)%coord)
            enddo
    
            if(iface.eq.1)then
                do k=1,nb_point2
                    M = mesh2%point(k)%coord
                    call CEta0(M,t,eta)
                    mesh2%point(k)%coord(3)=eta
                enddo
            endif 


            ! Mise à jour maillage    
            call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
            &               ind_pt,ind_arr,1,iface,ierror)
            if(ierror/=0)then
                ierror = 102
                goto 9999
            endif
    
            if(allocated(xgrid)) deallocate(xgrid)
            if(allocated(ygrid)) deallocate(ygrid)
            if(allocated(matdref)) deallocate(matdref)
    
        end if
        
    end do
  
    ! ------------------------------------------------------------------------
    ! Maillage des spheres
    ! ------------------------------------------------------------------------
    n1 = n2
    n2 = fgeom%nsphere
    
    do j=1,n2

        if(fgeom%cmd_sph(j)==1)then
            
            if(idebug>0)then
                print*,'Maillage sphere : ',fgeom%sphere(j)%index
            end if
            
            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)
            sphere = fgeom%sphere(j)
            iface = sphere%index
            
            rep = sphere%repere
            rmax = sphere%radius
            xmin = -rmax ; xmax = rmax
            ymin = -rmax ; ymax = rmax
            if(symmetry) ymin = 0._rp
            bf = -1
            
            ! Recuperation des points appartenant à la surface
            call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
            &        ind_pt,ind_arr,iface,ierror)
            
            ! Projection dans le plan xOy
            call local_coord_rep(mesh2,nb_point2,nb_arete2,nb_tri2,rep,ierror)
            
            ! Orientation des arretes
            do jj=1,nb_arete2
                edge = mesh2%arrete(jj) 
                ! FIXME : gestion automatique de l'orientation des arretes
                M1 = mesh2%point(edge%iP(1))%coord
                M2 = mesh2%point(edge%iP(2))%coord
                M = 0.5_rp*(M1+M2)
                call assign_vector_coord(V,[M1(2)-M2(2),M2(1)-M1(1),0._rp])        
        
                bool = dot_product(M(1:2),V%coord(1:2)) > 0.
                if(Symmetry)then
                    if(abs(M1(2)-ymin).lt.Epsilon .and. abs(M2(2)-ymin).lt.Epsilon .and. M1(1)>M2(1))then
                        bool = .true.
                    endif
                endif
        
                if(bool)then
                    ip = mesh2%arrete(jj)%iP(1)
                    mesh2%arrete(jj)%iP(1) = mesh2%arrete(jj)%iP(2)
                    mesh2%arrete(jj)%iP(2) = ip
                    mesh2%arrete(jj)%dir = -mesh2%arrete(jj)%dir
                endif
            enddo

            ! Creation grille des tailles de reference
            !if(nb_point2.gt.0)then
                if(xmax-xmin.le.2.*hrefx)then
                    nx = 2
                else
                    nx = int((xmax-xmin)/hrefx)
                endif
                if(ymax-ymin.le.2.*hrefy)then
                    ny = 2
                else
                    ny = int((ymax-ymin)/hrefy)
                endif
        
                ! Protection contre nombre de maille null
                allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                &    matdref,nx,ny,iface,idref,dx,ierror,InputData)  
                if(ierror/=0)then
                    ierror = 203
                    goto 9999
                endif
            !endif  
    
            ! Maillage de la surface 2D
            if(allocated(xgrid) .and. allocated(ygrid))then
                call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
                    &                     iface,1,[sphere%radius],1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                if(ierror/=0)then
                    ierror = 101
                    goto 9999
                endif
            end if

            do k=1,nb_point2
                M = mesh2%point(k)%coord
                call proj2cart(M(1:2),sphere,mesh2%point(k)%coord)
            enddo
    
            if(iface.eq.1)then
                do k=1,nb_point2
                    M = mesh2%point(k)%coord
                    call CEta0(M,t,eta)
                    mesh2%point(k)%coord(3)=eta
                enddo
            endif 

            ! Mise à jour maillage    
            call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
            &               ind_pt,ind_arr,1,iface,ierror)
            if(ierror/=0)then
                ierror = 102
                goto 9999
            endif
    
            if(allocated(xgrid)) deallocate(xgrid)
            if(allocated(ygrid)) deallocate(ygrid)
            if(allocated(matdref)) deallocate(matdref)
            
        endif  
        
    enddo
    
    ! ------------------------------------------------------------------------  
    ! Maillage des plans
    ! ------------------------------------------------------------------------
    
    n1 = n2
    n2 = fgeom%nplan
    
    do j=1,n2
        
        if(fgeom%cmd_pl(j)/=0)then
            
            if(idebug>0)then
                print*,'Maillage plan :',j
            endif
            
            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)  
            plan = fgeom%plan(j)
            iface = plan%index
            bf = -9
            
            call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
            &    ind_pt,ind_arr,iface,ierror)
            
            call local_coord_rep(mesh2,nb_point2,nb_arete2,nb_tri2,plan%repere,ierror) 
            
            if(iface==HouleRF%index)then
                do k=1,nb_point2
                    mesh2%point(k)%coord(3)=0._rp
                end do
            endif
            
            ! Orientation des aretes
            ymin = plan%vmin
            ymax = plan%vmax
            xmin = plan%umin
            xmax = plan%umax
            
            if(idebug>0)then
                print*,'xmin = ',xmin,'; xmax = ',xmax,'; ymin = ',ymin,'; ymax = ',ymax
            endif
            
            do jj=1,nb_arete2
                edge = mesh2%arrete(jj)
                M1 = mesh2%point(edge%iP(1))%coord
                M2 = mesh2%point(edge%iP(2))%coord
                M = 0.5*(M1+M2)
                call assign_vector_coord(V,[M1(2)-M2(2),M2(1)-M1(1),0.d0])
                P = M+0.4*V%coord
                call loc2cart(P%coord,plan%repere,P%coord)
                call CEta0(P%coord,t,eta)
      
                bool = .false.
                
                if(abs(M1(2)-ymin).lt.Epsilon .and. abs(M2(2)-ymin).lt.Epsilon .and. M2(1) < M1(1))then
                    bool = .true.
                elseif(abs(M1(2)-ymax).lt.Epsilon .and. abs(M2(2)-ymax).lt.Epsilon .and. M2(1) > M1(1))then
                    bool = .true.
                elseif(abs(M1(1)-xmin).lt.Epsilon .and. abs(M2(1)-xmin).lt.Epsilon .and. M2(2) > M1(2))then
                    bool = .true.
                elseif(abs(M1(1)-xmax).lt.Epsilon .and. abs(M2(1)-xmax).lt.Epsilon .and. M2(2) < M1(2))then
                    bool = .true.
                elseif(P%coord(3) > eta)then 
                    bool = .true.
                endif
                if(bool)then
                    ip = mesh2%arrete(jj)%iP(1)
                    mesh2%arrete(jj)%iP(1) = mesh2%arrete(jj)%iP(2)
                    mesh2%arrete(jj)%iP(2)=ip
                endif
            enddo
            
            ! Creation grille des tailles de reference
            if(nb_point2.gt.0)then
                if(xmax-xmin.le.2.*hrefx)then
                    nx = 2
                else
                    nx = int((xmax-xmin)/hrefx)
                endif
                if(ymax-ymin.le.2.*hrefy)then
                    ny = 2
                else
                    ny = int((ymax-ymin)/hrefy)
                endif
                allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                &    matdref,nx,ny,iface,idref,dx,ierror,InputData)  
                if(ierror/=0)then
                    ierror = 203
                    goto 9999
                endif  
            endif         
            
            ! Maillage 2D de la surface 
            if(allocated(xgrid) .and. allocated(ygrid))then
                call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
                &                     iface,0,[1._rp],1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                if(ierror/=0)then
                    ierror = 101
                    goto 9999
                endif
            end if
            
            do k=1,nb_point2
                M = mesh2%point(k)%coord
                call loc2cart(M,plan%repere,mesh2%point(k)%coord)
            enddo
            
            if(plan%index==HouleRF%index)then          
                do k=1,nb_point2
                    P = mesh2%point(k)%coord
                    call CEta0(P%coord,t,eta)
                    mesh2%point(k)%coord(3) = eta
                enddo
            endif
            
            ! Mise a jour maillage 
            call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
            &               ind_pt,ind_arr,1,iface,ierror)
            if(ierror/=0)then
                ierror = 102
                goto 9999
            endif
            
            if(allocated(xgrid)) deallocate(xgrid)
            if(allocated(ygrid)) deallocate(ygrid)
            if(allocated(matdref)) deallocate(matdref)
            
        endif
        
    enddo  
    
    ! ------------------------------------------------------------------------
    ! Mesh of the axisym surfaces
    ! ------------------------------------------------------------------------
    n1 = n2
    n2 = fgeom%naxisym
    do j=1,n2
        if(fgeom%cmd_axi(j)==1)then
            
            if(idebug>0)then
                print*,'Maillage axisym : ',fgeom%axisym(j)%index
            endif
            call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)
            axisym = fgeom%axisym(j)
            iface = axisym%index
            
            rep = axisym%repere
            rep%origine = axisym%repere%origine
            rep%e1 = axisym%repere%e1
            rep%e2 = axisym%repere%e2
            rep%e3 = axisym%repere%e3
            rmax = maxval(abs(axisym%P(:,1)))
            xmin = -rmax ; xmax = rmax
            ymin = axisym%vmin ; ymax = axisym%vmax
            bf = -1
            
            ! Recuperation des points appartenant à la surface
            call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
            &        ind_pt,ind_arr,iface,ierror)
            
            ! Projection dans le plan xOy
            repj = axisym%repere
            repj%origine = axisym%repere%origine
            repj%e1 = axisym%repere%e1
            repj%e2 = axisym%repere%e3
            repj%e3 = axisym%repere%e2
            
            call local_coord_rep(mesh2,nb_point2,nb_arete2,nb_tri2,repj,ierror)
            
            ! Orientation of the edges
            call orient_hole_mesh2(mesh2,nb_point2,nb_arete2,0,0,ierror)
            
            ! Creation grille des tailles de reference
            if(nb_point2.gt.0)then
                if(xmax-xmin.le.2.*hrefx)then
                    nx = 2
                else
                    nx = int((xmax-xmin)/hrefx)
                endif
                if(ymax-ymin.le.2.*hrefy)then
                    ny = 2
                else
                    ny = int((ymax-ymin)/hrefy)
                endif
                
                ! Protection contre nombre de maille null
                allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
                call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
                &    matdref,nx,ny,iface,idref,dx,ierror,InputData)
                
                if(present(CartesianGrid))then
                    call write_Cartesian_grid(xgrid,ygrid,nx,ny,matdref,t,iface) ! Cartesian_grid.dat
                end if
                
                if(ierror/=0)then
                    ierror = 203
                    goto 9999
                endif
            endif  
            
            ! Maillage de la surface 2D
            npoint = axisym%npoint
            allocate(param(2*npoint+1))
            param(1) = npoint
            param(2:npoint+1) = axisym%P(1:npoint,1)
            param(npoint+2:2*npoint+1) = axisym%P(1:npoint,2)
            
            if(allocated(xgrid) .and. allocated(ygrid))then
                call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
                &                     iface,3,param,2*npoint+1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
                if(ierror/=0)then
                    ierror = 101
                    goto 9999
                endif
            end if
            
            repj = axisym%repere
            repj%origine = axisym%repere%origine
            repj%e1 = axisym%repere%e1
            repj%e2 = axisym%repere%e3
            repj%e3 = axisym%repere%e2
            do k=1,nb_point2
                M = mesh2%point(k)%coord
                call loc2cart(M,repj,mesh2%point(k)%coord)
            enddo
            
            if(iface.eq.1)then
                do k=1,nb_point2
                    M = mesh2%point(k)%coord
                    call CEta0(M,t,eta)
                    mesh2%point(k)%coord(3)=eta
                enddo
            endif 
            
            !   Mise à jour maillage    
            call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
            &               ind_pt,ind_arr,1,iface,ierror)
            
            if(ierror/=0)then
                ierror = 102
                goto 9999
            endif
                        
            if(allocated(xgrid)) deallocate(xgrid)
            if(allocated(ygrid)) deallocate(ygrid)
            if(allocated(matdref)) deallocate(matdref)
            if(allocated(param)) deallocate(param)
                        
        endif  

    enddo
	
	! ------------------------------------------------------------------------
    ! Mesh of Wigley surfaces
    ! ------------------------------------------------------------------------
    n1 = n2
    n2 = fgeom%nwigley

    do j=1,n2

        if(fgeom%cmd_wigley(j)==1)then

        if(idebug>0)then
            print*,'Maillage wigley : ',fgeom%wigley(j)%index
        endif

        call init_mesh(mesh2,nb_point2,nb_arete2,nb_tri2)
        wigley = fgeom%wigley(j)
        iface = wigley%index

        rep = wigley%repere
        xmin = wigley%umin ; xmax = wigley%umax
        ymin = wigley%vmin ; ymax = wigley%vmax
        bf = -1

        ! Recuperation des points appartenant à la surface
        call extract_point_mesh(mesh,nb_arete,mesh2,nb_point2,nb_arete2,&
        &        ind_pt,ind_arr,iface,ierror)

        ! Projection dans le plan xOy
        repj%origine = rep%origine
        repj%e1 = rep%e1
        repj%e2 = rep%e3
        repj%e3 = rep%e2
        call local_coord_rep(mesh2,nb_point2,nb_arete2,nb_tri2,repj,ierror)

        call orient_hole_mesh2(mesh2,nb_point2,nb_arete2,0,0,ierror)

        ! Creation grille des tailles de reference
        if(nb_point2.gt.0)then
            if(xmax-xmin.le.2.*hrefx)then
                nx = 2
            else
                nx = int((xmax-xmin)/hrefx)
            endif
            if(ymax-ymin.le.2.*hrefy)then
                ny = 2
            else
                ny = int((ymax-ymin)/hrefy)
            endif
            
            allocate(xgrid(nx),ygrid(ny),matdref(nx,ny,3))
            call create_dref_grid(mesh2,nb_arete2,xmin,xmax,ymin,ymax,xgrid,ygrid,&
            &    matdref,nx,ny,iface,idref,dx,ierror,InputData)  
            if(ierror/=0)then
            ierror = 203
            goto 9999
            endif
        endif  
        
        !  Maillage de la surface 2D
        
        if(allocated(xgrid) .and. allocated(ygrid))then
            call advancefront(mesh2,nb_point2,nb_arete2,nb_tri2,xgrid,ygrid,matdref,&
            &                     iface,4,[1._RP],1,nx,ny,ierror,InputData,tin,ioMeshevol,fgeom)
            if(ierror/=0)then
                ierror = 101
                goto 9999
            endif
        end if

        do k=1,nb_point2
            M = mesh2%point(k)%coord
            call loc2cart(M,repj,mesh2%point(k)%coord)
        enddo
    
        if(iface.eq.1)then
            do k=1,nb_point2
            M = mesh2%point(k)%coord
            call CEta0(M,t,eta)
            mesh2%point(k)%coord(3)=eta
            enddo
        endif 

        ! Mise à jour maillage    
        call updatemesh(mesh,nb_point,nb_arete,nb_tri,mesh2,nb_point2,nb_arete2,nb_tri2,&
        &               ind_pt,ind_arr,1,iface,ierror)
        if(ierror/=0)then
            ierror = 102
            goto 9999
        endif
    
        if(allocated(xgrid)) deallocate(xgrid)
        if(allocated(ygrid)) deallocate(ygrid)
        if(allocated(matdref)) deallocate(matdref)
       
        endif  
    
    enddo

    ! --------------------------------------------------------------   

    9999 continue
    
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error #',i3,' : mesh2D') 
    
    if(allocated(mesh2%point)) deallocate(mesh2%point)
    if(allocated(mesh2%arrete)) deallocate(mesh2%arrete)
    if(allocated(mesh2%tri)) deallocate(mesh2%tri)
    if(allocated(ind_pt)) deallocate(ind_pt)
    if(allocated(ind_arr)) deallocate(ind_arr)
    
end subroutine mesh2d

subroutine compute_mesh(maillage,t,fgeom_vect,nface,dx,mesh,nb_point,&
&                       nb_arete,nb_tri,ierror,InputData,fdomaine,tab2,n_tab2,n_tab)
    
    !f2py integer*1, dimension(1000)                        :: maillage
    type(TMaillage),intent(inout)                           :: maillage                 ! Mesh.
    real(rp),intent(in)                                     :: t                        ! Current time.
    !f2py integer*1, dimension(1000)                        :: fgeom_vect
    type(type_GeomVect),intent(in)                          :: fgeom_vect               ! Geometries.
    integer,intent(in)                                      :: nface                    ! Number of faces.
    real(rp),intent(in)                                     :: dx                       ! dx = dx1.
    !f2py integer*1, dimension(1000)                        :: mesh
    type(MGrid),intent(inout)                               :: mesh                     ! MGrid.
    integer,intent(inout)                                   :: nb_point,nb_arete,nb_tri ! Number of points, edges and triangles in the mesh.
    integer,intent(inout)                                   :: ierror                   ! Error floag.
    !f2py integer*1, dimension(1000)                        :: InputData
    type(InputDataStruct),intent(inout)                     :: InputData                ! Input data.
    !f2py integer*1, dimension(1000)                        :: fdomaine
    type(type_geom),intent(inout),optional                  :: fdomaine                 ! Geometry of the domain.
    !f2py integer*1, dimension(1000)                        :: tab2
    type(chaine_point_pt),dimension(:),intent(in),optional  :: tab2                     ! Table of the pointers toward the intersection points.
    integer,intent(inout),optional                          :: n_tab2,n_tab             ! Number of intersection curves and lines.
    
    integer                                                 :: j, iflag,jj,kk              ! Loop parameters.
    type(point)                                             :: P                        ! Point.
    real(rp)                                                :: eta                      ! Wave elevation.
    integer                                                 :: ios                      ! Output parameter.
    real(rp)                                                :: tin                      ! Plotting time for Advance_front.dat.
        
    ! This subroutine creates the mesh with the mesh strategy of CC.
    
    


    
    

    
    print*,"compute_mesh: Optimal point is searched with dx2 of the first floater."
    ierror = 0
    
    ! Following the mesh evolution with Tecplot.
    if (iwevol) then
        open(unit=ioevol,file="Advance_front.dat", iostat=ios)
        write(unit=ioevol,fmt='(150a)') 'VARIABLES = "X","Y","Z","Number","bf","nface"'
    end if
    
    if(.not.present(tab2) .and. n_tab2.ne.0) then
        iflag = 1
    else
        iflag = 0
    endif
    

    ! Wave elevation of each point
    if(iflag==1)then  
        do j=1,nb_point
            P = mesh%point(j)%coord
            call CEta0(P%coord,t,eta)
            mesh%point(j)%coord(3) = eta
        enddo
    endif  
  
    ! Maillage 0D : ajout point de la geometrie dans le maillage
    if(iprint==1) print*,'** mesh0D ...'
    
    ! Tank
    call mesh0D(mesh,nb_point,fdomaine,t,iflag,ierror)
    
    if(ierror/=0)then
        ierror = 220
        goto 9999
    endif
    
    ! Bodies
    if(is_body)then
        do jj = 1,NBodies
            if(fgeom_vect%Active(jj))then
                call mesh0D(mesh,nb_point,fgeom_vect%geom(jj),t,iflag,ierror)
                if(ierror/=0)then
                    ierror = 221
                    goto 9999    
                endif
            end if
        end do
    endif
    
    if(idebug>0)then
        if(is_body)then
            call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                &                'Mesh_0D.dat','point_0D.dat','arrete_0D.dat',ierror) 
        end if
    endif
    
    ! Maillage 1D de la courbe d'intersection
    if(iflag==0)then 
        if(is_body)then
            call mesh_inter(mesh,nb_point,nb_arete,InputData%dx2(1),ierror,tab2,n_tab2)
        end if
        if(ierror/=0)then
            ierror = 231
            goto 9999
        endif
    endif
        
    if(idebug>0)then
        if(is_body)then
            call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                &                'Mesh_intersection.dat','point_intersection.dat','arrete_intersection.dat',ierror) 
        end if
    endif
        
    
    print*, "nb_point :   ", nb_point, nb_arete, nb_tri
    print*, "()()()()()()()()()()()"
    
    
    
    ! Maillage 1D : maillage des aretes de la geometrie et de la courbe d'intersection.
    if(iprint==1) print*,'** mesh1D ...'
    
    if(present(tab2) .and. n_tab2.ne.0)then
        ! Bodies.
        do jj = 1,NBodies
            if(fgeom_vect%Active(jj))then
                call mesh1D(mesh,nb_point,nb_arete,fgeom_vect%geom(jj),1,InputData%dx2(jj),t,iflag,ierror,InputData,jj,tab2(:),n_tab2,n_tab)
            end if
        end do
        
        if(idebug>0)then
            call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                &                'Mesh1Bodies.dat','point1Bodies.dat','arrete1Bodies.dat',ierror) 
        endif
        
        ! Tank.
        call mesh1D(mesh,nb_point,nb_arete,fdomaine,0,dx,t,iflag,ierror,InputData,1,tab2(:),n_tab2)
        
    else
        ! Bodies.
        if(is_body)then
            do jj = 1,NBodies
                if(fgeom_vect%Active(jj))then
                    call mesh1D(mesh,nb_point,nb_arete,fgeom_vect%geom(jj),1,InputData%dx2(jj),t,iflag,ierror,InputData,jj)                    
                end if
            end do
        end if
        
        if(idebug>0)then
            if(is_body)then
                call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                    &                'Mesh1Bodies.dat','point1Bodies.dat','arrete1Bodies.dat',ierror) 
            end if
        end if
        
        ! Tank.
        call mesh1D(mesh,nb_point,nb_arete,fdomaine,0,dx,t,iflag,ierror,InputData,1)
                
    end if
    if(ierror/=0)then
        ierror = 230
        goto 9999
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh1.dat','point1.dat','arrete1.dat',ierror) 
    end if
        
    ! Maillage 2D  : maillage des faces du domaine.
    if(iprint==1) print*,'** mesh2D ...'
    tin = 0._RP
    ! Tank and free surface.
    call mesh2D(mesh,nb_point,nb_arete,nb_tri,fdomaine,t,dx,ierror,InputData,tin,ioevol,.true.) ! .true. means to write the Cartesian grid for the free surface in an output file if iwCartGrid = T in *.in.
    if(ierror/=0)then
        ierror = 240
        goto 9999
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh2domaine.dat','point2domaine.dat','arrete2domaine.dat',ierror) 
    endif
    
    ! Bodies.
    if(is_body)then
        do jj = 1,NBodies
            if(fgeom_vect%Active(jj))then
                call mesh2D(mesh,nb_point,nb_arete,nb_tri,fgeom_vect%geom(jj),t,InputData%dx2(jj),ierror,InputData,tin,ioevol,.true.) ! .true. means to write the Cartesian grid.
                if(ierror/=0)then
                    ierror = 241
                    goto 9999
                endif
            end if
        end do
    end if
    
    ! Closing Advance_front.dat
    if (iwevol) then
        close(ioevol)
    end if
        
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh2.dat','point2.dat','arrete2.dat',ierror) 
    endif
    
    if(Ra.gt.Epsilon)then
        if(iprint>0) print*,"smooth corner ..."
        call smooth_corner(mesh,nb_point,ierror,InputData)
        if(ierror/=0)then
            ierror = 242
            goto 9999
        endif
    endif
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh.dat','point.dat','arrete.dat',ierror) 
    endif
    
    1000 continue
    
    ! Transfert type fichier maillage 
    if(iprint==1) print *,'** transfert type donnees maillage'
    call change_type_mesh(mesh,nb_point,nb_tri,t,maillage,nface,fgeom_vect,ierror,InputData)
    
    ! Partial initialization of the physical characteristics of the bodies.
    call Partial_initialization(Maillage,InputData)
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh_change_type_mesh.dat','point_change_type_mesh.dat','arrete_change_type_mesh.dat',ierror) 
    endif
    
    ! Geometry initialization
    call GeomInit(maillage,fgeom_vect,t,InputData,.false.)
    
    if(ierror/=0)then
        ierror = 320
        goto 9999
    endif
    
    if(iprint==1) then
        print*," Fin Compute mesh : nb_point = ",maillage%Nnoeud,' nb_tri = ',maillage%Nfacette
    endif
    
    9999 continue
    if(ierror/=0)then
    write(*,99),ierror
    endif
    99 format('** error #',i3,' : pb. generation maillage')  
    
    if(ierror/=0 .and. idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh.dat','point.dat','arrete.dat',ierror) 
    endif
    
end subroutine compute_mesh

subroutine extract_inter(Maillage,mesh,nb_point,nb_arete,mesh_old,ierror,InputData)
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(in)          :: Maillage                                                         ! Mesh
    !f2py integer*1, dimension(1000)    :: mesh
    type(MGrid),intent(inout)           :: mesh                                                             ! New MGrid structure
    integer,intent(inout)               :: nb_point,nb_arete                                                ! Number of points and edges
    !f2py integer*1, dimension(1000)    :: mesh_old
    type(MGrid),intent(in)              :: mesh_old                                                         ! Former MGrid structure
    integer,intent(inout)               :: ierror                                                           ! Error flag
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                                                        ! InputData
    
    integer                             :: j, k, jk, jj, ind, nface, PreviousnFace, nfaceOther,nc,NumBody
    integer,dimension(nvmax)            :: face, PreviousFace, faceOther
    integer                             :: int_sl0, int_sl1, Ndouble, np0,npa
    !f2py integer*1, dimension(1000)    :: vertex,vt1,vt2,vt1Other,vt2Other
    type(MVertex)                       :: vertex, vt1, vt2, vt1Other, vt2Other
    !f2py integer*1, dimension(1000)    :: edge, edgeOther
    type(MEdge)                         :: edge, edgeOther
    logical                             :: bool
    logical                             :: AxiSymFSparts                                                    ! The intersection curve of an axisym floater is made of two parts, consequently the subroutine must see one body for two free surface parts.
    integer,dimension(NBodies,2)        :: nb_point_tab                                                     ! Saving the number of points for each floater. (:,1) : 1st point, (:2) last point.
    
    ! This subroutines stores the points and the edges of the intersection curve between the floaters and the free surface of the TMaillage into MGrid. Then the subroutine orientates the edges for the advance front method. 
    
    ierror = 0
    
    int_sl0 = Maillage%FS%IndFS(1)
    int_sl1 = Maillage%FS%IndFS(3)
    
    np0 = nb_point+1
        
    nb_point_tab = 0._RP
    nb_point_tab(1,1) = nb_point+1 ! 1st point of the 1st floater
    
    npa = nb_arete
    
    ! Looking for the nodes at the intersection between the free surface and the piercing bodies
    do nc = 1,NBodies
        
        NumBody = nc + 1 ! 0 is FS, 1 is tank.
        
        ! Loop over the nodes of the free surface
        do j = int_sl0,int_sl1
            
            Ndouble = Maillage%Tnoeud(j)%Ndouble
            
            ! Nodes of the free surface in intersection with the floaters
            bool = .false.
            do k = 1,Ndouble
                jk = Maillage%Tnoeud(j)%double(k)
                bool = bool .or. Maillage%Tnoeud(jk)%Npanneau .eq. NumBody ! Intersection FS - Floaters
            enddo
        
            if(bool)then
                ind = Maillage%Tnoeud(j)%indmesh
                vertex = mesh_old%point(ind)
                vertex%coord = Maillage%Tnoeud(j)%Pnoeud(1:3)
                
                vertex%bf = 10 ! Boundary flag of the flotaer
                if(Symmetry .and. abs(vertex%coord(2).lt.Epsilon2))then
                    vertex%bf = 2 ! Boundary flag of the free surface on the symmetry plan
                endif
                                
                call add_element(mesh%point,nb_point,vertex) ! add_element_point
            endif
            
        enddo
                
        ! Last point of the current body
        nb_point_tab(nc,2) = nb_point
        
        ! 1st point of the next body
        if(nc.ne.NBodies)then ! Not for the last body
            nb_point_tab(nc+1,1) = nb_point+1
        end if
        
    end do
    
    ! FIXME : cas axisym: works only when the floater is on the plan of symmetry (even if there is no symmetry).
    if(not(Symmetry) .and. InputData%igtype(1).eq.5)then ! 1 and not Int_Body or NumBody !!!
        do j=np0,nb_point
            vt1 = mesh%point(j)
            if(abs(vt1%coord(2)).lt.Epsgeom)then ! Bad condition !!!
                mesh%point(j)%nface = vt1%nface+1
                mesh%point(j)%face(3) = 5
            elseif(vt1%coord(2)<0)then ! Negative y
                mesh%point(j)%face(2) = 5
            endif
        enddo
    endif
    
    ! On suppose que les points sont rangés dans l'ordre.
    jj = 0
    nc = 1
    AxiSymFSparts = .false.
    do j=np0,nb_point-1 ! -1 because of vt2
        vt1 = mesh%point(j)
        vt2 = mesh%point(j+1)
        
        ! Keep the previous nface
        if(jj.ge.1)then
            PreviousnFace = nface
            PreviousFace(1:nface) = face(1:nface)
        end if
        
        ! nface in common
        call common_int(vt1%face,vt1%nface,vt2%face,vt2%nface,face,nface,ierror)
        
        ! In case of several bodies, an edge would go from one body to another one, the only face in common is 1 (free surface). This edge is deleted.
        if(nface.eq.1 .and. face(nface).eq.1) cycle ! Only the free surface in common
        
        ! face(1) or face(2) can be equal to -1 in case of symmetry (for a cylinder for instance)
        if(nface.eq.2 .and. (face(1).eq.-1 .or. face(2).eq.-1)) cycle ! In case of Symmetry
        
        ! PreviousFace .ne. face means a change of floater, consequently the mesh should be closed (linked the last point to the first one).
        if((PreviousFace(1).ne.face(1) .or. PreviousFace(2).ne.face(2)) .and. jj.ge.1)then
            
            if(.not.Symmetry)then
                
                if(InputData%igtype(nc).ne.5 .or. (InputData%igtype(nc).eq.5 .and. AxisymFSparts))then
                    
                    ! Creating a new edge to close the intersection curve for this floater
                    vt1Other = mesh%point(j-1) ! Last point of the previous body
                    vt2Other = mesh%point(nb_point_tab(nc,1)) ! 1st point of the previous body
                    call common_int(vt1Other%face,vt1Other%nface,vt2Other%face,vt2Other%nface,faceOther,nfaceOther,ierror)
                    call edge_from_vertex(vt1Other,vt2Other,edgeOther)
                    edgeOther%bf = 1
                    edgeOther%nface = nfaceOther
                    edgeOther%face(1:nfaceOther) = faceOther(1:nfaceOther)
                    call add_element(mesh,nb_arete,edgeOther) ! add_element_arrete
                end if
                
                ! Axisym case
                if(InputData%igtype(nc).eq.5)then
                    if(.not.AxisymFSparts)then
                        AxisymFSparts = .true.
                        goto 3 ! No updating of nc, npa and no call to orient_hole_mesh.
                    else
                        AxisymFSparts = .false.
                    end if
                end if
                
            end if
            
            call orient_hole_mesh(mesh,nb_point_tab(nc,2),nb_arete,nb_point_tab(nc,1)-1,npa,ierror)
            
            ! Updating npa and nc
            npa = nb_arete
            nc = nc + 1
            
            3 continue
            
        end if
        
        call edge_from_vertex(vt1,vt2,edge) 
        edge%bf = 1
        edge%nface = nface
        edge%face(1:nface) = face(1:nface)
        call add_element(mesh,nb_arete,edge) ! add_element_arrete
        jj = jj + 1
    end do
    
    ! Closing of the mesh for the last floater (linked the last point to the first one).
    if(.not.Symmetry)then
        vt1 = mesh%point(nb_point) ! Last point of the last body
        vt2 = mesh%point(nb_point_tab(nc,1)) ! 1st point of the last body
        call common_int(vt1%face,vt1%nface,vt2%face,vt2%nface,face,nface,ierror)
        call edge_from_vertex(vt1,vt2,edge)
        edge%bf = 1
        edge%nface = nface
        edge%face(1:nface) = face(1:nface)
        call add_element(mesh,nb_arete,edge) ! add_element_arrete
    end if
    
    call orient_hole_mesh(mesh,nb_point,nb_arete,nb_point_tab(nc,1)-1,npa,ierror)
    
end subroutine extract_inter

subroutine Interpolation_FS(OldMesh,NewMesh,Ecoulement,t,LocalRemeshFS,ierror)
    
    !f2py integer*1, dimension(1000)    :: OldMesh,NewMesh
    type(TMaillage),intent(inout)       :: OldMesh,NewMesh  ! Old and new meshs.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement       ! Old and new flow parameters.
    real(rp),intent(in)                 :: t                ! Current time.
    logical,intent(in)                  :: LocalRemeshFS    ! = true: remeshing the free surface, = false: no remeshing.
    integer,intent(inout)               :: ierror           ! Error flag.
    
    integer                             :: j,k              ! Loop parameters.
    integer,allocatable,dimension(:)    :: ClosestPoints    ! Closest point of the old mesh for every point of the new mesh.
    real(rp)                            :: dist,PjPk        ! Distances between two points.
    integer                             :: Na, Nvoisin      ! Spline ordre and number of neighbours.
    real(rp),allocatable                :: Pvoisin(:,:)     ! Position of the neighbours.
    real(rp),allocatable                :: A(:,:), B(:,:)   ! Matrix A and vector B for the B-spline problem.
    integer                             :: info             !
    character(len=1)                    :: trans            !
    integer                             :: lda, nrhs        ! Size of the matrix A and number of rhs.
    integer, allocatable                :: ipiv(:)          !
    type(TEcoulement)                   :: EcoulementTmp    ! Temporary flow parameters.
    
    ! This subroutine initializes Phi_p and Eta_p on the free surface of the new mesh from the old one.
    
    ! Allocation
    allocate(ClosestPoints(NewMesh%FS%IndFS(3)-NewMesh%FS%IndFS(1)+1))
    ClosestPoints = 0
    
    ! EcoulementTmp
    call NewEcoulement(EcoulementTmp, NewMesh%Nnoeud)
    
    ! Initialisation de l'écoulement à 0
    call IniEcoulement(EcoulementTmp, NewMesh%Nnoeud, 0._RP)
    
    ! Interpolation if FS remeshing.
    if(LocalRemeshFS)then
    
        do j = NewMesh%FS%IndFS(1),NewMesh%FS%IndFS(3) ! New mesh
                
            ! Searching the closest points of NewMesh%Tnoeud(j)%Pnoeud.
            dist = 999._RP
            do k = OldMesh%FS%IndFS(1), OldMesh%FS%IndFS(3) ! Old mesh
                PjPk = norm2(NewMesh%Tnoeud(j)%Pnoeud-OldMesh%Tnoeud(k)%Pnoeud)
                if(PjPk .le. dist)then
                    dist = PjPk
                    ClosestPoints(j) = k
                end if
            end do
        
            if(is_BS)then
            
                ! Spline ordre at the new mesh point (may we could choose the closest point of the old mesh). Actually it is the same, %Ordre is given for all points of the FS. 
                Na = Nordre(NewMesh%Tnoeud(j)%Ordre)
            
                ! Number of neighbours of the closest point.
                NVoisin = OldMesh%Tnoeud(ClosestPoints(j))%NVoisin(2)-1 ! -1 because the point itself is counted.
            
                ! Allocations
                allocate(Pvoisin(3,NVoisin+1), B(NVoisin+Na,2), A(NVoisin+Na,NVoisin+Na))
            
                ! Pvoisin and B
                do k = 1,NVoisin+1 
                
                    ! Neighbours of the points in the NEW mesh.
                    Pvoisin(1:3,k) = OldMesh%Tnoeud(abs(OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(k,1)))%Pnoeud ! k = 1 is OldMesh%Tnoeud(ClosestPoints(j)) itself.
                
                    if(OldMesh%Tnoeud(abs(OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(k,1)))%NPanneau.ne.0)then
                        print*,"Interpolation_FS: Point is not on the FS:"
                        print*,j
                        pause
                    end if
                
                    if (OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(k,1).lt.0) Pvoisin(2,k) = - Pvoisin(2,k)
                    B(k,1) = Ecoulement%Phi(abs(OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(k,1)))%perturbation
                    B(k,2) = Ecoulement%Eta(abs(OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(k,1)))%perturbation
                
                end do
                B(Nvoisin+2:Nvoisin+Na,1:2) = 0._RP
            
                ! A
                call SplineMatrix(Nvoisin, NewMesh%Tnoeud(j)%Ordre, Pvoisin, A) ! If Na is given by the closest point of the old mesh: NewMesh%Tnoeud(j)%Ordre --> OldMesh%Tnoeud(ClosestPoints(j))%Ordre.
            
                ! Inversion of the linear system
                lda = Nvoisin+Na
                allocate(ipiv(lda))
                ipiv = 0
                trans = 'n'
                nrhs = 2 ! Number of rhs.
                call dgetrf(lda,lda,A,lda,ipiv,info)
                if (info.ne.0) then
                    print*,'WARNING Interpolation: info =', info, lda, j, OldMesh%Tnoeud(ClosestPoints(j))%NVoisin(2), Na
                    print*,OldMesh%Tnoeud(ClosestPoints(j))%TVoisin(1:OldMesh%Tnoeud(ClosestPoints(j))%NVoisin(2),1)
                    pause
                end if
                call dgetrs(trans,lda,nrhs,A,lda,ipiv,B,lda,info) ! The B-spline coefficients are in B.
                deallocate(ipiv)
            
                ! Phi_p
                call SplineF(Nvoisin, NewMesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,1), NewMesh%Tnoeud(j)%Pnoeud, Pvoisin, EcoulementTmp%Phi(j)%perturbation) ! NewMesh%Tnoeud(j)%Pnoeud is where SplineF is computed. If Na is given by the closest point of the old mesh: NewMesh%Tnoeud(j)%Ordre --> OldMesh%Tnoeud(ClosestPoints(j))%Ordre.
            
                ! Eta_p
                call SplineF(Nvoisin, NewMesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,2), NewMesh%Tnoeud(j)%Pnoeud, Pvoisin, EcoulementTmp%Eta(j)%perturbation) ! NewMesh%Tnoeud(j)%Pnoeud is where SplineF is computed. If Na is given by the closest point of the old mesh: NewMesh%Tnoeud(j)%Ordre --> OldMesh%Tnoeud(ClosestPoints(j))%Ordre.
                        
                ! Deallocations
                if(allocated(Pvoisin)) deallocate(Pvoisin)
                if(allocated(B)) deallocate(B)
                if(allocated(A)) deallocate(A)
            
            else
                print*,"Interpolation_FS: interpolation algorithm only works with B-splines."
                call exit()
            end if
        
        end do
        
    else ! No FS remeshing: necessity to update the size of Ecoulement.
        
        do j = NewMesh%FS%IndFS(1),NewMesh%FS%IndFS(3) ! New mesh
            
            ! Phi_p.
            EcoulementTmp%Phi(j)%perturbation = Ecoulement%Phi(j)%perturbation ! Same j because structures are store in the following way: (FS, tank, bodies) and only tank and bodies change in size.
            
            ! Eta_p.
            EcoulementTmp%Eta(j)%perturbation = Ecoulement%Eta(j)%perturbation ! Same j because structures are store in the following way: (FS, tank, bodies) and only tank and bodies change in size.
            
        end do
            
    end if
    
    ! Copy Ecoulement
    call DelEcoulement(Ecoulement)
    call NewEcoulement(Ecoulement, NewMesh%Nnoeud)
    call IniEcoulement(Ecoulement, NewMesh%Nnoeud, 0._RP) ! %incident and %perturbation are 0.
    call CopyEcoulement(Ecoulement, EcoulementTmp, NewMesh%Nnoeud)
    call DelEcoulement(EcoulementTmp)
    
    ! Computing the incident flow. If the nodes change their location, the incident quantity are not updated automatically, hence the call to Incident.
    call Incident(t, NewMesh, Ecoulement)
    
    ! Deallocation
    if(allocated(ClosestPoints)) deallocate(ClosestPoints)
    
end subroutine Interpolation_FS

subroutine compute_mesh_fgeom(Maillage,Ecoulement,t,fgeom_vect,fdomaine,nface,mesh_old,nb_point_old,nb_tri_old,ierror,InputData,nRemesh,LocalRemeshFS,tab2,n_tab2,n_tab)
        
    !f2py integer*1, dimension(1000)                        :: Maillage
    type(TMaillage),intent(inout)                           :: Maillage                     ! Mesh.
    !f2py integer*1, dimension(1000)                        :: Ecoulement
    type(TEcoulement),intent(inout)                         :: Ecoulement                   ! Old and new flow parameters.
    real(rp),intent(inout)                                  :: t                            ! Current time.
    !f2py integer*1, dimension(1000)                        :: fgeom_vect
    type(type_GeomVect),intent(in)                          :: fgeom_vect                   ! Geometry of the floaters.
    !f2py integer*1, dimension(1000)                        :: fdomaine
    type(type_geom),intent(in)                              :: fdomaine                     ! Geometry of the domain.
    integer,intent(in)                                      :: nface                        ! Total number of faces.
    !f2py integer*1, dimension(1000)                        :: mesh_old
    type(MGrid),intent(inout)                               :: mesh_old                     ! MGrid.
    integer                                                 :: nb_point_old, nb_tri_old     ! Number of points and triangles in the mesh_old.
    integer,intent(inout)                                   :: ierror                       ! Error flag.
    !f2py integer*1, dimension(1000)                        :: InputData
    type(InputDataStruct),intent(inout)                     :: InputData                    ! Input data.
    integer,intent(in)                                      :: nRemesh                      ! Number of remeshing (use for the number of Advance_Front_Remesh.dat).
    logical,intent(inout)                                   :: LocalRemeshFS                ! = true: remeshing the free surface, = false: no remeshing.
    !f2py integer*1, dimension(1000)                        :: tab2
    type(chaine_point_pt),dimension(:),intent(in),optional  :: tab2                         ! Table of the pointers toward the intersection points.
    integer,intent(inout),optional                          :: n_tab2,n_tab                 ! Number of intersection curves and lines.
            
    !f2py integer*1, dimension(1000)                        :: MaillageTemp
    type(TMaillage)                                         :: MaillageTemp                 ! Temporary mesh.
    !f2py integer*1, dimension(1000)                        :: mesh
    type(MGrid)                                             :: mesh                         ! MGrid.
    integer                                                 :: nb_point, nb_arete, nb_tri   ! Number of points, edges and triangles in the mesh.
    integer                                                 :: iflag                        ! Flag to know if there are given inersection points or not (tab2, n_tab2).
    integer                                                 :: jj                           ! Loop parameter.
    integer                                                 :: ios                          ! Output parameter.
    real(rp)                                                :: tin                          ! Printing time for Advance_front_Remesh.dat.
    character(len=200)                                      :: fileRemesh                   ! Temporary file name for the remeshing.
    character (len=50)                                      :: filemaillRemesh              ! Name of the output remesh file.
            
    ! This subroutine creates a new mesh of the floaters.
    
    print*,""
    print*,"Remeshing number ",nRemesh
    print*,""
    
    ! LocalRemeshFS.
    if(not(RemeshFS))then
        ! If RemeshFS = false: no remeshing of the free surface anyway.
        LocalRemeshFS = .false.
        
        ! If RemeshFS = true: remeshing the FS mesh if one of its panels needs it or if a body goes across the free surface or if a forced remeshing is aksed.
    else
        LocalRemeshFS = .true.
        ! LocalRemeshFS was a first trial to remesh the FS without remeshing the floaters. It does not work right now.
        ! It is necessary to create a new subroutine to insert the mesh of the FS in keeping the mesh of the floaters.
        ! That is to say, creating a "new" insert_mesh_global (just add the mesh of the floaters, not the FS) -> insert_meshFS_global. TO DO.
        ! So far if the FS is remeshed without the floater meshes, that leads to numerical problem at the intersection curves.
    end if
        
    if(LocalRemeshFS)then
        if(not(present(tab2)) .and. n_tab2.ne.0) then
            iflag = 1
        else
            iflag = 0
        endif
    else
        iflag = 0
    end if
    
    ! Following the mesh evolution with Tecplot.
    if (iwevolRemesh) then
        write(fileRemesh,'("Advance_front_Remesh_",i3,".dat")'),nRemesh
        open(unit=ioevolRemesh,file=fileRemesh,iostat=ios)
        write(unit=ioevolRemesh,fmt='(150a)') 'VARIABLES = "X","Y","Z","Number","bf","nface"'
    end if
    tin = 0._RP
    
    ! Correction of mesh_old in case of axisym geometries
    if(not(LocalRemeshFS))then
        do jj = 1,NBodies
            if(not(Symmetry) .and. InputData%igtype(jj).eq.5)then
                call demerge_face(Mesh_old,nb_point_old,nb_tri_old,fgeom_vect%geom(jj)%axisym(1)%index,fgeom_vect%geom(jj)%axisym(2)%index,ierror)
                if(ierror/=0) goto 9999
            endif
        end do
    end if
    
    ! Nouveau maillage local 
    call init_mesh(mesh,nb_point,nb_arete,nb_tri) 
    
    ! Mesh0D of the tank
    if(LocalRemeshFS)then
        call mesh0D(mesh,nb_point,fdomaine,t,iflag,ierror)
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh0D_Remesh_Domain.dat','point0D_Remesh_Domain.dat','arrete0D_Remesh_Domain.dat',ierror) 
    end if
    
    ! Mesh0D of the bodies
    if(is_body)then
        do jj = 1,NBodies
            if(idebug>0)then
                print*,"Mesh 0D : ",jj
            end if
            if(fgeom_vect%Active(jj))then
                call mesh0D(mesh,nb_point,fgeom_vect%geom(jj),t,iflag,ierror)
            end if
            if(ierror/=0) goto 9999
        end do
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh0D_Remesh.dat','point0D_Remesh.dat','arrete0D_Remesh.dat',ierror) 
    end if
        
    ! Recupération des points de l'intersection
    if(is_body)then
        if(LocalRemeshFS)then
            if(iflag==0)then 
                call mesh_inter(mesh,nb_point,nb_arete,InputData%dx2(1),ierror,tab2,n_tab2)
                
                if(idebug>0)then
                    call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                        &                'Mesh_Remesh_Intersection.dat','point_Remesh_Intersection.dat','arrete_Remesh_Intersection.dat',ierror) 
                end if                
                
            endif
        else
            call extract_inter(Maillage,mesh,nb_point,nb_arete,mesh_old,ierror,InputData)
            
            if(idebug>0)then
                call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
                    &                'Mesh_Remesh_Extract_inter.dat','point_Remesh_Extract_inter.dat','arrete_Remesh_Extract_inter.dat',ierror) 
            end if
            
        end if
        if(ierror/=0) goto 9999
    end if
        
    if(ierror/=0) goto 9999     
    
    ! Maillage des aretes de la geometrie
    if(is_body)then
        do jj = 1,NBodies
            if(idebug>0)then
                print*,"Mesh 1D : ",jj
            end if
            
            ! Bodies
            if(fgeom_vect%Active(jj))then
                if(LocalRemeshFS)then
                    call mesh1D(mesh,nb_point,nb_arete,fgeom_vect%geom(jj),1,InputData%dx2(jj),t,iflag,ierror,InputData,jj,tab2(:),n_tab2,n_tab)
                else
                    call mesh1D(mesh,nb_point,nb_arete,fgeom_vect%geom(jj),1,InputData%dx2(jj),t,iflag,ierror,InputData,jj)
                end if
            end if
            if(ierror/=0) goto 9999
        end do
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh1D_Remesh_Bodies.dat','point1D_Remesh_Bodies.dat','arrete1D_Remesh_Bodies.dat',ierror) 
    end if
        
    ! Tank
    if(LocalRemeshFS)then
        call mesh1D(mesh,nb_point,nb_arete,fdomaine,0,dx1,t,iflag,ierror,InputData,1,tab2(:),n_tab2)
    end if
    if(ierror/=0) goto 9999
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh1D_Remesh.dat','point1D_Remesh.dat','arrete1D_Remesh.dat',ierror) 
    end if
    
    ! Maillage 3D
    
    ! Tank and free surface
    if(LocalRemeshFS)then
        call mesh2D(mesh,nb_point,nb_arete,nb_tri,fdomaine,t,dx1,ierror,InputData,tin,ioevolRemesh,.true.) ! .true. means to write the Cartesian grid for the free surface in an output file if iwCartGrid = T in *.in.
    end if
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh2D_Remesh_Domain.dat','point2D_Remesh_Domain.dat','arrete2D_Remesh_Domain.dat',ierror) 
    endif
    if(ierror/=0) goto 9999
    
    if(is_body)then
        do jj = 1,NBodies
            if(idebug>0)then
                print*,"Mesh 2D : ",jj
            end if
            ! Bodies
            if(fgeom_vect%Active(jj))then
                call mesh2D(mesh,nb_point,nb_arete,nb_tri,fgeom_vect%geom(jj),t,InputData%dx2(jj),ierror,InputData,tin,ioevolRemesh)
            end if
            if(ierror/=0) goto 9999
        end do
    end if
    
    ! Closing Advance_front_Remesh.dat
    if(iwevolRemesh)then
        close(ioevolRemesh)
    end if
    
    ! Change face number for non symetrical case
    if(not(LocalRemeshFS))then
        do jj = 1,NBodies
            if(not(Symmetry) .and. InputData%igtype(jj).eq.5)then
                call merge_face(Mesh,nb_point_old,nb_tri_old,fgeom_vect%geom(jj)%axisym(1)%index,fgeom_vect%geom(jj)%axisym(2)%index,ierror)
                if(ierror/=0) goto 9999
            endif
        end do
    end if
    
    if(idebug>0)then
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
            &                'Mesh2D_Remesh.dat','point2D_Remesh.dat','arrete2D_Remesh.dat',ierror) 
    endif
    
    if(ierror/=0) goto 9999
    
    ! Maillage global
    call NewMaillage(MaillageTemp, PointMax,NBodies+1) ! +1 for the tank.
    call CopyMaillage(MaillageTemp, Maillage) 
    
    ! Mesh -> MaillageTemp
    if(LocalRemeshFS)then
        call change_type_mesh(mesh,nb_point,nb_tri,t,MaillageTemp,nface,fgeom_vect,ierror,InputData)
    else
        call insert_mesh_global(mesh,nb_point,nb_tri,MaillageTemp,nface,fgeom_vect,ierror,InputData,.false.) ! LocalRemeshFS = .false. is it true ?
    end if
    
    ! Initialization of the characterisitics of the mesh.
    call GeomInit(MaillageTemp,fgeom_vect,t,InputData,.false.,ierror)
    if(ierror/=0) goto 9999
    
    ! Creation of the new structure Ecoulement (good size) (even without FS remeshing!) and interpolation of Phi_p and Eta_p between the old mesh and the new one.
    call Interpolation_FS(Maillage,MaillageTemp,Ecoulement,t,LocalRemeshFS,ierror)

    ! Old TMaillage -> New TMaillage   
    call DelMaillage(Maillage)
    call NewMaillage(Maillage,MaillageTemp%Nnoeud,NBodies+1) ! +1 for the tank.
    call CopyMaillage(Maillage,MaillageTemp)
    call DelMaillage(MaillageTemp)
    
    ! Writting the mesh in an output file
    filemaillRemesh = 'Maillage_Remesh.dat'
    call PlotMaill(filemaillRemesh, Maillage)
    
    9999 continue
        
    if(ierror/=0)then
        
        ! Writting the structure MGrid in an output file.
        call write_debug(mesh,nb_point,nb_arete,nb_tri,t,&
        &                'Err_mesh_Remesh.dat','Err_point_Remesh.dat','Err_arrete_Remesh.dat',ierror)
        
        ! Writting the structure TMaillage in an output file.
        filemaillRemesh = 'MaillageTemp_Remesh.dat'
        call PlotMaill(filemaillRemesh, MaillageTemp)
        
    endif
    
    if(allocated(mesh%point)) deallocate(mesh%point)
    if(allocated(mesh%arrete)) deallocate(mesh%arrete)
    if(allocated(mesh%tri)) deallocate(mesh%tri)
    if(ierror/=0)then 
        write(*,99),ierror
    endif
    99 format('** error #',i3,' compute_mesh_fgeom')  
        
end subroutine compute_mesh_fgeom

end module MeshModule

