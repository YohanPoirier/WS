subroutine TestGenMaillage
!program Main
use Constantes
use Structuresdonnees
use GenMaillage
use PrePlot
use BoucleTemp
use Incident_mod
use Airy
use Houle_RF
implicit none

! Variables
character (len=50) :: filemaill
integer :: j
integer, parameter :: nt=300
real(rp) :: t0, h
real(rp), dimension(nt) :: t
type(TMaillage) :: MaillageInit, MaillageT
type(TEcoulement) :: EcoulementDT
type(THouleRF) :: HouleRF
type(type_geom) :: fgeom
type(type_geom) :: domaine1
real(rp) :: prof,dx
integer :: nface
type(repere3d) :: rep0



h=0.01
t0=-h+pi/(w(1)*2._RP)
t=(/(t0+j*h,j=1,nt)/)
call NewMaillage(MaillageInit, PointMax)
!call GenMaillage(t(1), MaillageInit)
call bodygen(MaillageInit)
!call mesh_cyl(MaillageInit, 0)
filemaill='Mesh%dat'
call PlotMaill(filemaill, MaillageInit)

!call Extract_WaveRF(HouleRF)
!
!print*, HouleRF.Landa, HouleRF.Hcc, HouleRF.k, HouleRF.T, HouleRF.C, HouleRF.Cp, HouleRF.Cc, HouleRF.N
!call NewMaillage(MaillageT, PointMax)
!call Remaillage(t(1), MaillageInit, MaillageT)
!!call Remaillage(0._RP, MaillageInit, MaillageT)
!filemaill='MaillageT.dat'
!call PlotMaill(filemaill, MaillageT)

call BoucleTemporelle(nt, t, MaillageInit,HouleRF,fgeom,domaine1,nface,rep0,dx)

!end program Main
end subroutine TestGenMaillage