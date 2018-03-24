program Main
use GenMaillage
use GeomMesh
use PrePlot
use BoucleTemp
use BoucleTemp_lineaire
use Exec
use BVP
use Validation_Solver
use Validation_Gradient
implicit none
! ************************************************************************
! Declaration Variables
! ************************************************************************
! Variables
character (len=50) :: filemaill, ligne
integer :: j, ios, ierror
real(rp) :: periode, time_begin, time_end, nPeriod
real(rp), allocatable :: t(:)
real(rp), dimension(3) :: Origine
type(TMaillage) :: MeshInit
! Begin
call Execution
! Mesh 
call NewMaillage(MeshInit, PointMax)
if (idtype.eq.1) then
    ! Rectangular Mesh
    call bodygen(MeshInit)
elseif (idtype.eq.2) then
    ! Cylindrical Mesh 
    Origine = [0._RP, 0._RP, 0._RP]
    call mesh_cyl(MeshInit, Origine)
elseif (idtype.eq.0) then
    ! Extraction from a file
    filemaill = 'Maillage.dat'
    call extract_maillage(filemaill,MeshInit)
end if
! Write Mesh File
filemaill = 'MeshInit.dat'
call PlotMaill(filemaill, MeshInit)
t_Solver = .false.
t_Spline = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Test Solveur!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (t_Solver) call Test_Solver(MeshInit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Test B-Splines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (t_Spline) call Test_Spline(MeshInit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! */ Advance in Time
!1 periode = 2._RP*Pi/w(1) ! Period
!nPeriod = 8._RP ! Number of periods 
!dt = nPeriod*periode/nt ! Time Step
allocate(t(nt+1))
!go to 999
print*, 'Periode =', periode, 'dt = ', dt, 'Temps de simulation = ', dt*nt
t=(/((j-1)*dt,j=1,nt+1)/)
call BoucleTemporelle(nt+1, t, MeshInit, ierror)
!call BoucleTemporelle_lineaire(nt+1, t, MeshInit, ierror)
999 deallocate(t)
! End
!end if
    end program Main
