module Unit_Test
use Test_Solver
use Test_Gradient
use Test_Lissage
!use Test_MeshDeformation
use GeomDef
implicit none
contains
    
subroutine UnitTests(Mesh, fgeom, nface, rep0, Grid)
! ************************************************************************
! Declaration Variables
! ************************************************************************
!f2py integer*1, dimension(1000)                    :: Mesh
type(TMaillage), intent(inout)                      :: Mesh
!f2py integer*1, dimension(1000)                    :: fgeom
type(type_geom), intent(inout)                      :: fgeom                                ! Geometry of the floater (not the domain)
integer, intent(in)                                 :: nface                                ! Number of faces in both the floater and the domain geometries
!f2py integer*1, dimension(1000)                    :: rep0
type(repere3d), intent(in)                          :: rep0                                 ! Inertial frame
!f2py integer*1, dimension(1000)                    :: Grid
type(MGrid), intent(inout)                          :: Grid 

logical :: t_Spline, t_Solver, t_Lissage, t_MeshDeformation

! This subroutine calls the unit tests.

t_Solver = .false.
if (t_Solver) call TestSolver(Mesh)

t_Spline = .false.
if (t_Spline) call TestGradient(Mesh)

t_Lissage = .true. .and. Cd_Morison.ne.0
if (t_Lissage) call TestLissage(Mesh)

t_MeshDeformation = .false.
if (t_MeshDeformation .and. Mesh_type.eq.2) call TestMeshDeformation(Mesh,fgeom,nface,rep0,Grid)


end subroutine UnitTests

end module Unit_Test