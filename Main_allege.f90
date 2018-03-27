   
Program Main
    use Exec_mod
    use GeomDef
    use MeshGen
    use CoeffInfluence
    use CoeffInfluence_CL
    use Structuresdonnees

    
    character (len=50)                      :: fileparam,filegeom                           ! Input files
    type(InputDataStruct)                   :: InputData                                    ! Input data.
    logical                                 :: get_State                                    ! True if the state input file is present, false otherwise.

    type(type_GeomVect)                     :: fgeom_vect                                   ! Geometry of the floaters (not the domain).
    type(type_geom)                         :: fdomaine                                     ! Geometry of the domain.
    integer                                 :: nface,nline                                  ! Number of faces and lines in both the floater and the domain geometries.
    type(chaine_point_pt),dimension(100)    :: tab2                                         ! Table of intersection points ! Why 100? (PYW).
    integer                                 :: n_tab2,n_tab                                 ! Number of intersection curves and lines.
    type(repere3d)                          :: rep0                                         ! Inertial frame.
    integer                                 :: ierror                                       ! Error flag.

    type(TMaillage)                         :: Mesh                                         ! Final mesh (named Maillage in Main.f90). 
    type(MGrid)                             :: Grid                                         ! Transitional mesh (named Mesh in Main.f90).
    integer                                 :: nb_point,nb_tri                              ! Number of points and triangles in the Mgrid.
 
    real(rp),allocatable                    :: CD(:,:), CS(:,:)                             ! Influence coefficients
    real(rp),allocatable                    :: CD2(:,:), CS2(:,:)                             ! Influence coefficients
    
    real(rp), allocatable :: L_X(:,:),  L_n(:,:), L_GS(:,:,:), L_G(:,:), L_R_max(:)
    integer :: Nnodes, i_f_i, i_f_f, N_sym
    integer, allocatable :: L_T(:,:)

    

    integer :: k,j
      
    fileparam = 'ws.in'
    filegeom = 'test.geom'

    get_State = .false.


    ! Lecture input
    call Execution(fileparam,filegeom,InputData,get_State)

    ! Generation Geometry
    call Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)

    ! Generation maillage
    call Generation_Mesh(Mesh,fdomaine,fgeom_vect,nface,Grid,nb_point,nb_tri,ierror,InputData,get_State,tab2,n_tab2,n_tab)
   
    
    
    ! CI calculation (cudaLIKE)
    Nnodes = Mesh%Nnoeud
    Nfacettes = Mesh%Nfacette
    allocate(CS2(Nnodes,Nnodes), CD2(Nnodes,Nnodes))
    CD2 = 0._RP ; CS2 = 0._RP
    
    allocate(L_X(Nnodes,3))
    allocate(L_T(Nfacettes,3), L_n(Nfacettes,3), L_GS(Nfacettes,3,3), L_G(Nfacettes,3), L_R_max(Nfacettes))
    do k = 1, Nnodes
        L_X(k, :) = Mesh%Tnoeud(k)%Pnoeud
    end do
    
    
    do k = 1, Nfacettes
        L_T(k,:) = Mesh%Tfacette(k)%Tnoeud
        L_n(k,:) = Mesh%Tfacette(k)%normale
        L_GS(k,:,:) = Mesh%Tfacette(k)%ds
        L_G(k,:) = Mesh%Tfacette(k)%Gfacette
        L_R_max(k) = Mesh%Tfacette(k)%Rmax
    end do
    
    i_f_i = 1
    i_f_f = Nfacettes
    N_sym = 1

    call  mat_CI_kernel(L_X, L_T, L_n, L_GS, L_G, L_R_max, CD2, CS2, Nnodes,Nfacettes, i_f_i, i_f_f, N_sym, Ldom(3))
    call angle_solide_kernel(CD2, Nnodes)
    
    ! CI calculation
    Nnodes = Mesh%Nnoeud
    allocate(CS(Nnodes,Nnodes), CD(Nnodes,Nnodes))
    CD = 0._RP ; CS = 0._RP
    call CoeffInfl(Mesh, CD, CS,Nnodes)
    
    
    print *, "CS"
    do k = 1, 5
        do j = 1, 5
            print*, k,j,CS(k,j) -CS2(k,j)
        end do
    end do
    
    print*, "CD"
    do k = 1, 5
        do j = 1, 5
            print*, k,j,CD(k,j), CD2(k,j)
        end do
    end do
    
    pause
    
end program