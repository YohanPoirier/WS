module Validation_Rotation
use Constantes
use Parameters
use StructuresDonnees
use FonctionsCommunes
use BoucleTemp
use GenMaillage

implicit none

contains

subroutine rotation(Mesh)
    type(TMaillage) :: Mesh, Mesh2
    character(len=50) :: filemaill, num
    type(TEcoulement) :: Ecoulement
    integer :: j

    !call Initialisation(Ecoulement, Mesh, 0._Rp)
    Htype=0
    filemaill = 'RotationMaillage_'//trim(num)//'.dat'
    Mesh%Body(Int_Body)%MBody(1:3)=[0._RP,0._RP,0._RP]
    Mesh%Body(Int_Body)%MBody(4:6)=Position(:,2)
    call NewMaillage(Mesh2,Mesh%Nnoeud)
    Mesh2=Mesh
    open(unit=98, file='RotMesh_'//filename)
    do j=1,100
        Mesh%Body(Int_Body)%CMD(1)=.true.
        Mesh%Body(Int_Body)%Gbody(1:3,1)=Mesh%Body(Mesh%NBody)%Gbody(1:3,1)+[0._RP,0._RP,0._RP]
        Mesh%Body(Int_Body)%Mbody(4:6)=[0._RP,0._RP+2*PI*j*0.01_RP,0._RP]
        call remesh(Mesh,Mesh2, 0._RP, 0._RP)
        !call PlotPropHoule(j*1.0_RP, Mesh, Ecoulement, 0.001_Rp)
        call plotMaillage(Mesh, j*0.001_Rp)
    end do
    close(98)
end subroutine rotation

subroutine plotMaillage(mesh,t)
 type(TMaillage), intent(in) :: Mesh
  !type(TEcoulement), intent(inout), optional :: Ecoulement
  real(rp), intent(in) :: t
  !real(rp), dimension(3,Mesh%Nnoeud), intent(in), optional :: Fpre0, Fpre1, Fpre2
! Variables Locales
  integer :: j
  character(len=10) :: num
!
  write( num, '( f0.4 )' ) t
!

  write(98,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette,&
  &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
  do j=1,Mesh%Nnoeud
     write(98,'(19E)') Mesh%Tnoeud(j)%Pnoeud(1:2), Mesh%Tnoeud(j)%Pnoeud(3)
  end do
!
  do j=1,Mesh%Nfacette
    write(98,'(3I)') Mesh%Tfacette(j)%Tnoeud
  end do
end subroutine plotMaillage
end module Validation_Rotation
