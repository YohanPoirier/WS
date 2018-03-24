module MeshGen
use Parameters
use MeshStruct
use MeshModule
use PrePlot
use GenMaillage
implicit none

contains

subroutine Generation_Mesh(maillage,fdomaine,fgeom_vect,nface,mesh,nb_point,nb_tri,ierror,InputData,get_State,tab2,n_tab2,n_tab)
    
    !f2py integer*1, dimension(1000)                        :: maillage
    type(TMaillage),intent(inout)                           :: maillage                 ! Total mesh (domain + floater).
    !f2py integer*1, dimension(1000)                        :: fdomaine
    type(type_geom),intent(inout)                           :: fdomaine                 ! Geometry of the domain.
    !f2py integer*1, dimension(1000)                        :: fgeom_vect
    type(type_GeomVect),intent(inout)                       :: fgeom_vect               ! Geometry of the bodies.
    integer,intent(in)                                      :: nface                    ! Number of faces in both the floater and the domain.
    !f2py integer*1, dimension(1000)                        :: mesh
    type(MGrid),intent(out)                                 :: mesh                     ! Transitional mesh.
    integer,intent(out)                                     :: nb_point,nb_tri          ! Number of points and triangles in the Mgrid.
    integer                                                 :: ierror                   ! Error flag.
    !f2py integer*1, dimension(1000)                        :: InputData
    type(InputDataStruct),intent(inout)                     :: InputData                ! Input data.
    logical,intent(in)                                      :: get_State                ! True if the state input file is present, false otherwise.
    !f2py integer*1, dimension(1000)                        :: tab2
    type(chaine_point_pt),dimension(:),intent(in),optional  :: tab2                     ! Table of intersection points.
    integer,intent(inout),optional                          :: n_tab2,n_tab             ! Number of intersection curves and lines.
    
    integer                                                 :: nb_arete                 ! Number of edges in the Mgrid.
    integer                                                 :: iarg                     ! Argument about the status of the third input file.
    character (len=50)                                      :: filemaill                ! Name of the output mesh file.
    real(rp), dimension(3)                                  :: Origine                  ! Origine of the inertial frame.
    logical                                                 :: get_mesh                 ! Flag to known if there is a mesh input file.
    
    ! This subroutine generates mesh of the domain, the free surface and the floater.
    
    ! Initialization of the final mesh structure
    call NewMaillage(maillage, PointMax,NBodies+1) ! +1 for the tank.
    
    ! Testing if there is a mesh input file
    !call get_command_argument(3,filemaill,status=iarg) ! This third input is not a mesh anymore (State.txt).
    
    ! Flag to know if there is a mesh input file or not
    !if(iarg==0)then
    !    get_mesh = .true.
    !else
    !    get_mesh = .false.
    !endif
    get_mesh = .false. ! This command is totaly desactivated.
    
    ! Reading of the mesh input file
    if(get_mesh)then
        call Extract_Maillage(filemaill,fgeom_vect,maillage,InputData)
        print*,' Utilisation du maillage initial : ',filemaill,' [',maillage%Nnoeud,' noeuds ; ',&
        &      maillage%Nfacette,' facettes]'
    endif
  
    print*,"Compute mesh ..."
    
    ! Mesh generation of Lucas
    if(Mesh_type==1 .and. .not.get_mesh)then
        Origine = [0._rp, 0._rp, 0._rp]
        if(idtype.eq.1)then ! Domain = Rectangular
            print*,"Generation_Mesh: igtype of the fist body."
            if (InputData%igtype(1).eq.6) then ! Specific mesh generation for a Wigley Hull
                call Mesh_cart_Wigley(maillage, InputData,fgeom_vect,Origine)
            else
                call bodygen(maillage,fgeom_vect,InputData)
            end if
        elseif(idtype.eq.2)then ! Domain = Cylinder
            print*,"Generation_Mesh: igtype of the fist body."
            if (InputData%igtype(1).eq.6) then ! Specific mesh generation for a Wigley Hull
                call MeshFS_WigleyHull(maillage,InputData,fgeom_vect,Origine)
            else
                call mesh_cyl(maillage, fgeom_vect,InputData,Origine)
            end if
        elseif(idtype.eq.0)then ! Domain = mesh input file
            filemaill = 'Maillage.dat'
            call extract_maillage(filemaill,fgeom_vect,maillage,InputData)
        else
            ierror = 200
            goto 9999
        endif
        
    ! Mesh generation of Camille
    elseif(Mesh_type==2 .and. .not.get_mesh)then
        
        ! Initilization of the transitional mesh.
        call init_mesh(mesh,nb_point,nb_arete,nb_tri)
        if(ierror/=0)then
            ierror = 210
            goto 9999
        endif
        
        ! Creation nouveau maillage.
        call NewMaillage(maillage,Pointmax,NBodies+1) ! +1 for the tank.
        
        ! Generation of the final mesh.
        call compute_mesh(maillage,t0,fgeom_vect,nface,dx1,mesh,&
        &                 nb_point,nb_arete,nb_tri,ierror,InputData,fdomaine,tab2,n_tab2,n_tab)
        if(ierror/=0)then
            ierror = 211
            goto 9999
        end if
        
    endif

    print*,'Nnoeud = ',maillage%Nnoeud
    Print*,'Nfacette = ',maillage%Nfacette
    
    if(Mesh_type .eq. 1)then ! Method of Lucas
        if(allocated(mesh%point)) deallocate(mesh%point)
        if(allocated(mesh%arrete)) deallocate(mesh%arrete)
        if(allocated(mesh%tri)) deallocate(mesh%tri)
        call deallocate_geom(fdomaine,ierror)
        !call deallocate_GeomVect(fgeom_vect,ierror)
    endif 

    ! Initial mesh output file
    filemaill = 'Maillage.dat'
    call PlotMaill(filemaill, maillage) 
    
    ! STL format
    call ExportMesh_STL(Maillage,fgeom_vect,.true.)
    
    ! Checking of the surface area of the floater
    call verif_surf(maillage,fgeom_vect)
    
    ! Checking of the normal vector length
    call verif_normale(maillage,fgeom_vect)
    
    print*,"End of Generation Mesh"
    
    ! Error
    9999 continue     
    if(ierror/=0)then
        write(*,99),ierror
        stop
    endif
    99 format('** error #',i3,' : Generation_Mesh')
    
end subroutine Generation_Mesh

subroutine verif_surf(Maillage,fgeom_vect)

    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(in)          :: Maillage     ! Total mesh (domain + floater)
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect   ! Geometry of the floaters.
    
    integer                             :: j,k,n1,n2,nc ! Parameters
    real(rp)                            :: Dn, S1, S2   ! Normal, surfaces, etc
    real(rp),dimension(3)               :: FH0, DS, Xf  ! Hydrostatic loads, etc
    
    ! This subroutine checks the quality of mesh in computing of the area of the body.
    
    do nc = Int_Body,Maillage%Nbody
        
        if(fgeom_vect%Active(nc-1))then
            n1 = Maillage%Body(nc)%IndBody(1)
            n2 = Maillage%Body(nc)%IndBody(3)
        
            S1 = 0._rp
            FH0 = 0._rp
            do j=n1,n2
                Dn = 0._rp 
                do k=1,Maillage%Tnoeud(j)%Nfacette
                    Dn = Dn + Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Aire*Maillage%Tnoeud(j)%Angle(k)
                enddo
                DS = Dn*Maillage%Tnoeud(j)%Normale
                DS = Maillage%Tnoeud(j)%Aire*Maillage%Tnoeud(j)%Normale
                S1 = S1 + Dn
                FH0 = FH0 + Maillage%Tnoeud(j)%Pnoeud(3)*DS
            enddo
        
            Xf = 0._rp
            do j=n1,n2
                Dn = 0._rp 
                do k=1,Maillage%Tnoeud(j)%Nfacette
                    Dn = Dn + Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Aire*Maillage%Tnoeud(j)%Angle(k)
                enddo
                Xf = Xf + Maillage%Tnoeud(j)%Pnoeud(1:3)*Dn
            enddo
        
            Xf = Xf/S1
        
            n1 = Maillage%Body(nc)%IndBody(2)
            n2 = Maillage%Body(nc)%IndBody(4)    
        
            S2 = 0._rp
            do j=n1,n2
                S2 = S2 + Maillage%Tfacette(j)%Aire
            enddo
        
            if(Symmetry)then
                S1 = 2._rp*S1
                S2 = 2._rp*S2
                FH0 = 2._rp*FH0 
                FH0(2) = 0._rp
            endif
        
            write(*,*) '-------------------------------------------'
            write(*,*) 'Body: ',nc
            write(*,*) 'Compute surface from point : S1 = ',S1
            write(*,*) 'Compute surface from face  : S2 = ',S2    
            write(*,*) 'Volume = ', FH0(3)
            write(*,*) 'FH0(1) = ', FH0(1)*ro
            write(*,*) 'FH0(2) = ', FH0(2)*ro
            write(*,*) 'FH0(3) = ', FH0(3)*ro
            write(*,'(a5,3f8.4)') 'Xf = ',Xf
        
        end if
        
    end do
    write(*,*) '-------------------------------------------'

end subroutine verif_surf

subroutine verif_normale(Maillage,fgeom_vect )
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(in)          :: Maillage                 ! Total mesh (domain + floater)
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect               ! Geometry of the floaters.
    
    integer                             :: j,n1, n2, nc             ! Parameters
    real(rp)                            :: dmin, dmax, dmean, val   ! length of the normal vectors
    
    ! This subroutine checks the quality of mesh in computing of the size of the normal vectors.
    
    do nc = Int_Body,Maillage%NBody
        
        if(fgeom_vect%Active(nc-1))then
            n1 = Maillage%Body(nc)%IndBody(1)
            n2 = Maillage%Body(nc)%IndBody(3)
        
            dmax = -999.
            dmin = 999.
        
            dmean = 0._rp
            do j=n1,n2
                val = norm2(Maillage%Tnoeud(j)%Normale)
                dmean = dmean + val
                if(val.gt.dmax) dmax = val
                if(val.lt.dmin) dmin = val
                if (dot_product(Maillage%Tnoeud(j)%Pnoeud-Maillage%Body(nc)%GBody(1:3),Maillage%Tnoeud(j)%Normale).lt.Epsilon) then
                    print*, 'Erreur Normale Point :',j, ' de normale : ', Maillage%Tnoeud(j)%Normale
                end if
            enddo
            dmean = dmean / (n2-n1+1)
        
            write(*,*) 'Body: ',nc
            write(*,*) 'Normalisation des vecteurs normaux :'
            write(*,*) ' dmean = ',dmean
            write(*,*) ' dmin  = ',dmin
            write(*,*) ' dmax  = ',dmax
            write(*,*) '-------------------------------------------'
        end if
        
    end do

end subroutine verif_normale

end module MeshGen

    