program Main
use Constantes
use Structuresdonnees
use Parameters
use GenMaillage
use GeomMesh
use PrePlot
use Exec

implicit none
! ************************************************************************
! Declaration Variables
! ************************************************************************
! Variables
character (len=50) :: filemaill, ligne, filemaill2
integer :: j, ios, ierror
real(rp) :: periode, time_begin, time_end, nPeriod
real(rp), allocatable :: t(:)
real(rp), dimension(3) :: Origine
type(TMaillage) :: MeshInit
! Begin
call Execution
igtype=3
! Mesh
call NewMaillage(MeshInit, PointMax)
Origine = [0._RP, 0._RP, 0._RP]
call mesh_cyl(MeshInit, Origine)

! Write Mesh File
filemaill = 'Mesh.dat'
call PlotMaill(fileMaill,MeshInit)
end program Main
