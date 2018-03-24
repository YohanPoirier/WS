   
Program Main
    use Exec_mod
    use GeomDef
    use MeshGen
    use CoeffInfluence
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
      
    fileparam = 'ws.in'
    filegeom = 'test.geom'

    get_State = .false.


    ! Lecture input
    call Execution(fileparam,filegeom,InputData,get_State)

    ! Generation Geometry
    call Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)

    ! Generation maillage
    call Generation_Mesh(Mesh,fdomaine,fgeom_vect,nface,Grid,nb_point,nb_tri,ierror,InputData,get_State,tab2,n_tab2,n_tab)
   
    ! CI calculation
    Nnodes = Mesh%Nnoeud
    allocate(CS(Nnodes,Nnodes), CD(Nnodes,Nnodes))
    CD = 0._RP ; CS = 0._RP
    call CoeffInfl(Mesh, CD, CS,Nnodes)
    
end program