include 'mkl_pardiso.f90'
Program Main
use Exec_mod
use GeomDef
use MeshGen
use BoucleTemp
!use Unit_Test
implicit none


! Variables
type(TMaillage)                         :: Maillage             ! Final mesh.
type(TMaillage)                         :: Maillage_State       ! Mesh read in the input file (WARNING: only Tnoeud and Tfacette are filled).
type(TEcoulement)                       :: Ecoulement_State     ! Flow parameters read in the state input file.
integer                                 :: nface                ! Number of faces in both the geometries of thefloaters and the domain.
integer                                 :: ierror               ! Error flag.
type(type_geom)                         :: fdomaine             ! Geometry of the domain.
type(type_GeomVect)                     :: fgeom_vect           ! Geometries of the floaters.
type(chaine_point_pt),dimension(100)    :: tab2                 ! Table of intersection points ! Why 100? (PYW).
integer                                 :: n_tab2,n_tab         ! Number of intersection curves and lines.
type(repere3d)                          :: rep0                 ! Inertial frame.
type(MGrid)                             :: mesh                 ! MGrid.
integer                                 :: nb_point,nb_tri      ! Number of points and triangles in the Mgrid.
character(len=50)                       :: fileparam,filegeom   ! *.ws and *.geom input files.
character(len=50)                       :: fileState            ! State input file.
type(InputDataStruct)                   :: InputData            ! Input data.
integer                                 :: iargState            ! = 0 (state input file present), otherwise no third input file.
logical                                 :: get_State            ! True if the state input file is present, false otherwise.
integer                                 :: jt0                  ! Initial time step parameter.
real(rp)                                :: t_tmp                ! Back-up of t0 in case of input state file.

character(len=150)                      :: fileMesh
integer                                 :: j

! This is the main function of the weak-scatterer code.

!if (is_Licence) call Licence

! *.ws
call get_command_argument(1,fileparam)

! *. geom
call get_command_argument(2,filegeom)

! State.txt
call get_command_argument(3,fileState,status=iargState) ! This fourth input is not a mesh anymore.
if (iargState.eq.0) then
    get_State = .true.
else
    get_State = .false.
endif

! Reading the input files *.ws and *.geom.
call Execution(fileparam,filegeom,InputData,get_State)

!fileMesh = 'Sphere_WSC_clipped.dat'
!call Extract_Mesh_Body(filemesh,Maillage,fgeom_vect,InputData)
!ierror = 0
!call GeomInit(Maillage, fgeom_vect, 0._RP, InputData,.false., ierror)
!
!fileMesh = 'Sphere_WSC_before_deformation.dat'
!call PlotMaill(fileMesh,Maillage)
!
!call MeshVelBody_Vertex(Maillage,0._RP,InputData,0._RP)
!
!fileMesh = 'Bonjour.dat'
!call PlotMaill(fileMesh,Maillage)
!
!call GeomInit(Maillage, fgeom_vect, 0._RP, InputData,.false., ierror)
!
!fileMesh = 'Sphere_WSC_after_deformation.dat'
!call PlotMaill(fileMesh,Maillage)
!
!call exit()

! Updating InputData in case of state input file.
if(get_State)then
    t_tmp = t0 ! Back-up t0 of *.in.
    call read_State_InputData(fileState,InputData,t0,jt0)
    print*,""
    print*,"Simulation starts at t = ",t0," s"
    print*,""
end if

! Generation of the geometries.
call Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)

if(ierror/=0)then
    print*,""
    print*,"Main: Error in the generation of the geometrie."
    go to 10
end if

! Generation of the mesh.
call Generation_Mesh(Maillage,fdomaine,fgeom_vect,nface,mesh,nb_point,nb_tri,ierror,InputData,get_state,tab2,n_tab2,n_tab)
if(get_State)then
    t0 = t_tmp ! t0 is used in the generation of the mesh (explicit wave elevation).
end if

if(ierror/=0)then
    print*,""
    print*,"Main: Error in the generation of the mesh."
    go to 10
end if

! Unit-Tests.
!call UnitTests(Maillage, fgeom, nface, rep0, mesh)

! Temporal loop.
call BoucleTemporelle_RK4(Maillage, fgeom_vect, fdomaine, nface, mesh,nb_point,nb_tri,InputData,n_tab2,n_tab,get_State,fileState,jt0)

10 continue

! Deallocating.
if(allocated(mesh%point)) deallocate(mesh%point)
if(allocated(mesh%arrete)) deallocate(mesh%arrete)
if(allocated(mesh%tri)) deallocate(mesh%tri)

call deallocate_geom(fdomaine,ierror)
if(Mesh_type .eq. 1)then
    call deallocate_GeomVect(fgeom_vect,ierror)
end if
call Delete_InputData(InputData)

end Program Main
