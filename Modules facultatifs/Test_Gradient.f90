module Test_Gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                  !
!!                        Validation du calcul des gradients sur la SL                              !
!!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! On impose une solution analytique sur la perturbation. On calcule les dérivées globales, locales 
!! sur la surface libre et le corps. On compare à la solution analytique imposée.
!!
!!23/05/2013 Lucas LETOURNEL
!!
use Constantes
use Parameters
use Incident_mod
use Spline
use FonctionsCommunes

implicit none

contains

subroutine TestGradient(Mesh)
implicit none
! Parameters
type(TMaillage) :: Mesh
! Variables
real(rp) :: time_begin, time_end
type(TEcoulement) :: Ecoulement
! Begin
call NewEcoulement(Ecoulement, Mesh%Nnoeud)
call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)
! Initialisation of the solution
call InitGradient(Ecoulement, Mesh)
! Derivatives Calculation
call cpu_time(time_begin)
call Gradient(Mesh, Ecoulement)
call cpu_time(time_end)
print*, time_end - time_begin
! Plot of the solution
call PlotFSGradient(Mesh, Ecoulement)
call PlotBodyGradient(Mesh, Ecoulement)
call DelEcoulement(Ecoulement)
! End
end subroutine TestGradient

subroutine InitGradient(Ecoulement, Mesh)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
! Locals
integer :: j, typeInput
real(rp) :: t, aph, bph, cph, expa, expb, expc
real(rp), dimension(3) :: M, GPhiLocal
real(rp), dimension(3,3) :: Puvw, Pt
type(THoule) :: HouleTemp, HouleObj
! Begin
typeInput = 0
if (typeInput.eq.0) then
    ! Airy or RF Wave solution
    t = 0._RP
    ! Definition of the temporary wave for the test
    call TempWave(HouleObj)
    call copyHoule(HouleObj, HouleTemp)
    ! Calculation of the incident components
    call Incident(t, Mesh, Ecoulement)
    ! Copy back the original wave
    call CopyBackHoule(HouleTemp)
    deallocate(HouleObj%w, HouleObj%dir, HouleObj%Aphi, HouleObj%konde, HouleObj%lambda)
elseif (typeInput.eq.1) then
    ! vertical oscillating sphere in infinite domain
    do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
        Mesh%Body(Int_Body)%VBody(1:3) = [0._RP, 0._RP, 1._RP]
        Ecoulement%Phi(j)%Incident = -0.5*Mesh%Body(Int_Body)%DimBody(1)*dot_product(Mesh%Body(Int_Body)%VBody(1:3),Mesh%TNoeud(j)%Normale)
        Ecoulement%DPhiDn(j)%Incident = dot_product(Mesh%Body(Int_Body)%VBody(1:3), Mesh%Tnoeud(j)%Normale)
        if (.true.) then
            Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
            Pt(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
            ! 1er Ordre
            GPhiLocal(1) = -0.5_RP*dot_product(Mesh%Body(Int_Body)%VBody(1:3), Puvw(1:3,1))
            GPhiLocal(2) = 0._RP
            GPhiLocal(3) = Ecoulement%DPhiDn(j)%Incident    
            Ecoulement%GPhi(1:3,j)%Incident = matmul(Puvw,GPhiLocal)
            Ecoulement%GPhi2(1:2,j)%Incident = GPhiLocal(1:2)    
            ! 2nd Ordre
            Ecoulement%GPhi2(3:4,j)%Incident =  0.5_RP*Mesh%TNoeud(j)%DLocal(3:4)*Ecoulement%DPhiDn(j)%Incident
        else    
            M = Mesh%Tnoeud(j)%Pnoeud(1:3) - Mesh%Body(Int_Body)%GBody(1:3)
            Ecoulement%GPhi(1,j)%Incident = 1.5_RP*Mesh%Body(Int_Body)%VBody(3)*M(1)*M(3)/Mesh%Body(Int_Body)%DimBody(1)**2
            Ecoulement%GPhi(2,j)%Incident = 1.5_RP*Mesh%Body(Int_Body)%VBody(3)*M(2)*M(3)/Mesh%Body(Int_Body)%DimBody(1)**2
            Ecoulement%GPhi(3,j)%Incident =-0.5_RP*Mesh%Body(Int_Body)%VBody(3)*(1._RP-(M(3)/Mesh%Body(Int_Body)%DimBody(1))**2)
            Ecoulement%DPhiDn(j)%Incident = dot_product(Ecoulement%GPhi(1:3,j)%Incident, Mesh%Tnoeud(j)%Normale)
            Ecoulement%GradGrad(1:3,1:3,j)%Incident = 0._RP
        end if
    end do    
elseif(typeInput.eq.2) then
    ! Potentiel Phi = aph.x^(expa)+bph.y^(expb)+cph.z^(expc)
    aph = 1._RP
    bph = -2._RP
    cph = 0._RP
    expa = 4._RP
    expb = 4._RP
    expc = 4._RP
    do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
        M = Mesh%Tnoeud(j)%Pnoeud
        Ecoulement%Phi(j)%Incident = aph*M(1)**expa + bph*M(2)**expb + cph*M(3)**expc
        Ecoulement%GPhi(1:3,j)%Incident = [expa*aph*M(1)**(expa-1), expb*bph*M(2)**(expb-1), expc*cph*M(3)**(expc-1)]
        Ecoulement%DPhiDn(j)%Incident = dot_product(Ecoulement%GPhi(1:3,j)%Incident, Mesh%Tnoeud(j)%Normale)
        GPhiLocal = matmul(Ecoulement%GPhi(1:3,j)%Incident,Mesh%Tnoeud(j)%Plocal(1:3,1:3,1))
        Ecoulement%GPhi2(1:2,j)%Incident = GPhiLocal(1:2)
        Ecoulement%GradGrad(1:3,1:3,j)%Incident = reshape([ expa*(expa-1)*aph*M(1)**(expa-2), 0._RP, 0._Rp, & 
                                                            0._RP, expb*(expb-1)*bph*M(2)**(expb-2), 0._RP,&
                                                            0._RP, 0._RP, expc*(expc-1)*cph*M(3)**(expc-2)],(/3,3/))
        Ecoulement%GPhi2(3,j)%incident = dot_product(matmul(Ecoulement%GradGrad(1:3,1:3,j)%Incident,Mesh%Tnoeud(j)%Plocal(1:3,1,1)),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
        Ecoulement%GPhi2(4,j)%incident = dot_product(matmul(Ecoulement%GradGrad(1:3,1:3,j)%Incident,Mesh%Tnoeud(j)%Plocal(1:3,2,1)),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
        Ecoulement%GPhi2(3:4,j)%incident = Ecoulement%GPhi2(3:4,j)%incident - Mesh%TNoeud(j)%DLocal(3:4)*Ecoulement%DPhiDn(j)%incident
    end do
end if
Ecoulement%Eta(1:Mesh%Nnoeud)%Perturbation = Ecoulement%Eta(1:Mesh%Nnoeud)%Incident
Ecoulement%Phi(1:Mesh%Nnoeud)%Perturbation = Ecoulement%Phi(1:Mesh%Nnoeud)%Incident
Ecoulement%DPhiDn(1:Mesh%Nnoeud)%Perturbation = Ecoulement%DPhiDn(1:Mesh%Nnoeud)%Incident
! End
end subroutine InitGradient

subroutine PlotFSGradient(Mesh, Ecoulement)
implicit none
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
! Variables
integer :: j, i1, i2, ios
real(rp), dimension(1:Mesh%FS%IndFS(3)) :: EpsGEta, EpsGPhi
real(rp), dimension(4) :: Maxi
! Begin
i1 = Mesh%FS%IndFS(3)
i2 = Mesh%FS%IndFS(4)
Maxi = 0._RP
do j=1,i1
    Ecoulement%GPhi(1:3,j)%incident = Ecoulement%GPhi(1:3,j)%incident - Ecoulement%DPhiDn(j)%incident*Mesh%Tnoeud(j)%Normale
    Maxi(1) = max(Maxi(1),norm2(Ecoulement%GEta(1:3,j)%incident))
    Maxi(2) = max(Maxi(2),norm2(Ecoulement%GPhi(1:3,j)%incident))
    if (Maxi(1) .lt. Epsilon) Maxi(1) = 1._RP
    if (Maxi(2) .lt. Epsilon) Maxi(2) = 1._RP
end do
do j=1,i1
    Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation - Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
    EpsGEta(j) = 100._RP*(norm2(Ecoulement%GEta(1:3,j)%perturbation - Ecoulement%GEta(1:3,j)%incident))/Maxi(1)
    EpsGPhi(j) = 100._RP*(norm2(Ecoulement%GPhi(1:3,j)%perturbation - Ecoulement%GPhi(1:3,j)%incident))/Maxi(2)
    Maxi(3) = max(Maxi(3),EpsGEta(j))
    Maxi(4) = max(Maxi(4),EpsGPhi(j))
end do
print*, Maxi(3:4)
! Plot
open(unit=3,file='TestGradient_FS_'//filename, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
write(3,fmt='(50a)') 'Title= "Gradients sur la SL"'
write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","EpsGEta","GEtax","GEtay","GEtax0","GEtay0","Eta","EpsGPhi","GPhix","GPhiy","GPhiz","GPhix0","GPhiy0","GPhiz0","Phi","Nvoisin" '
write(3,fmt='(a,i,a,i,a)') 'Zone N =', i1, ', E=', i2, ' , ET=TRIANGLE, F=FEPOINT' 
do j=1,i1
    write(3,'(17E,I)') Mesh%Tnoeud(j)%Pnoeud, EpsGEta(j), Ecoulement%GEta(1:2,j)%perturbation, Ecoulement%GEta(1:2,j)%incident, Ecoulement%Eta(j)%incident&
                                           &, EpsGPhi(j), Ecoulement%GPhi(1:3,j)%perturbation, Ecoulement%GPhi(1:3,j)%incident, Ecoulement%Phi(j)%incident, Mesh%Tnoeud(j)%Nvoisin(2)
end do
do j=1,i2
    write(3,'(3I)') Mesh%Tfacette(j)%TNoeud
end do 
close(unit=3)
! End
end subroutine PlotFSGradient

subroutine PlotBodyGradient(Mesh, Ecoulement)
implicit none
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp), dimension(1:Mesh%Nnoeud,1:2) :: EpsDPhi, EpsD2Phi
real(rp), dimension(1:4,1:Mesh%Nnoeud) :: GPhi
! Variables
integer :: j, k, i1, i2, ios
real(rp), dimension(4) :: Maxi
real(rp), dimension(3) :: Grad
real(rp), dimension(3,3) :: Puvw, Pt
! Begin
i1 = Mesh%Body(Int_Body)%IndBody(1)
i2 = Mesh%Body(Int_Body)%IndBody(3)
Maxi = 0._RP
do j=i1,i2
    if (.false.) then
        Puvw = Mesh%TNoeud(j)%Plocal(1:3,1:3,1)
        Pt = Mesh%TNoeud(j)%Plocal(1:3,1:3,2)
        !! Dérivées Premières Locales
        Grad = matmul(Pt,Ecoulement%GPhi(:,j)%incident)
        GPhi(1:2,j) = Grad(1:2)
        !! Dérivées Secondes Locales
        GPhi(3,j) = dot_product(matmul(Ecoulement%GradGrad(:,:,j)%incident,Puvw(1:3,1)),Puvw(1:3,1)) - 0*Mesh%TNoeud(j)%Dlocal(3)*dot_product(Ecoulement%GPhi(:,j)%incident,Mesh%TNoeud(j)%Normale)
        GPhi(4,j) = dot_product(matmul(Ecoulement%GradGrad(:,:,j)%incident,Puvw(1:3,2)),Puvw(1:3,2)) - 0*Mesh%TNoeud(j)%Dlocal(4)*dot_product(Ecoulement%GPhi(:,j)%incident,Mesh%TNoeud(j)%Normale)
    else
        GPhi(1:4,j) = Ecoulement%GPhi2(1:4,j)%incident
    end if
    do k=1,4
        Maxi(k) = Max(Maxi(k),abs(GPhi(k,j)))
    end do
end do

do j=i1,i2
    !! Dérivées Premières Locales
    EpsDPhi(j,1) = 100._RP*abs(Ecoulement%GPhi2(1,j)%perturbation - GPhi(1,j))/Maxi(1)
    EpsDPhi(j,2) = 100._RP*abs(Ecoulement%GPhi2(2,j)%perturbation - GPhi(2,j))/Maxi(2)
    !! Dérivées Secondes Locales
    EpsD2Phi(j,1) = 100._RP*abs(Ecoulement%GPhi2(3,j)%perturbation - GPhi(3,j))/Maxi(3)
    EpsD2Phi(j,2) = 100._RP*abs(Ecoulement%GPhi2(4,j)%perturbation - GPhi(4,j))/Maxi(4)
end do
! Plot
open(unit=3,file='TestGradient_Body_Local_'//filename, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
write(3,fmt='(50a)') 'Title= "Remaillage de la cuve"'
write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","EpsDPhiDs1","EpsDPhiDs2","EpsD2PhiDs1","EpsD2PhiDs2","0DPhiDs1","0DPhiDs2","0D2PhiDs1","0D2PhiDs2","DPhiDs1","DPhiDs2","D2PhiDs1","D2PhiDs2","Phi","DPhiDn","Nvoisin" '
write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(Int_Body)%IndBody(3)-Mesh%Body(Int_Body)%IndBody(1)+1, ', E=', Mesh%Body(Int_Body)%IndBody(4)-Mesh%Body(Int_Body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
    write(3,'(17E,I)') Mesh%Tnoeud(j)%Pnoeud, EpsDPhi(j,1:2), EpsD2Phi(j,1:2), GPhi(1:4,j), Ecoulement%GPhi2(1:4,j)%perturbation, Ecoulement%Phi(j)%incident, Ecoulement%DPhiDn(j)%incident, Mesh%Tnoeud(j)%Nvoisin(2)
end do
do j=Mesh%Body(Int_Body)%IndBody(2),Mesh%Body(Int_Body)%IndBody(4)
    write(3,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(Int_Body)%IndBody(1)-1)*[1,1,1]
end do 
close(unit=3)


i1 = Mesh%Body(Int_Body)%IndBody(1)
i2 = Mesh%Body(Int_Body)%IndBody(3)
Maxi = 0._RP
do j=i1,i2
    Ecoulement%GPhi(1:3,j)%incident = Ecoulement%GPhi(1:3,j)%incident - Ecoulement%DPhiDn(j)%incident*Mesh%Tnoeud(j)%Normale
    Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation-Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
    Maxi(1) = max(Maxi(1),norm2(Ecoulement%GPhi(1:3,j)%incident))
    if (Maxi(1) .lt. Epsilon) Maxi(1) = 1._RP
end do
do j=i1,i2
    EpsDPhi(j,1) = 100._RP*(norm2(Ecoulement%GPhi(1:3,j)%perturbation - Ecoulement%GPhi(1:3,j)%incident))/Maxi(1)
    !Maxi(4) = max(Maxi(4),EpsGPhi(j))
end do
open(unit=3,file='TestGradient_Body_Global_'//filename, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
write(3,fmt='(50a)') 'Title= "Remaillage de la cuve"'
write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","EpsGPhi","0DPhiDx","0DPhiDy","0DPhiDz","DPhiDx","DPhiDy","DPhiDz","Phi","DPhiDn","Nvoisin" '
write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(Int_Body)%IndBody(3)-Mesh%Body(Int_Body)%IndBody(1)+1, ', E=', Mesh%Body(Int_Body)%IndBody(4)-Mesh%Body(Int_Body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
do j=i1,i2
    write(3,'(12E,I)') Mesh%Tnoeud(j)%Pnoeud, EpsDPhi(j,1), Ecoulement%GPhi(1:3,j)%incident, Ecoulement%GPhi(1:3,j)%perturbation, Ecoulement%Phi(j)%incident, Ecoulement%DPhiDn(j)%incident, Mesh%Tnoeud(j)%Nvoisin(2)
end do
do j=Mesh%Body(Int_Body)%IndBody(2),Mesh%Body(Int_Body)%IndBody(4)
    write(3,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(Int_Body)%IndBody(1)-1)*[1,1,1]
end do 
close(unit=3)

open(unit=3,file='TestGradient_BaseLocale_'//filename, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
write(3,fmt='(50a)') 'Title= "Base Locale"'
write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","ux","uy","uz","vx","vy","vz","nx","ny","nz"'
write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(Int_Body)%IndBody(3)-Mesh%Body(Int_Body)%IndBody(1)+1, ', E=', Mesh%Body(Int_Body)%IndBody(4)-Mesh%Body(Int_Body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
    write(3,'(16E,I)') Mesh%Tnoeud(j)%Pnoeud, Mesh%Tnoeud(j)%Plocal(1:3,1,1), Mesh%Tnoeud(j)%Plocal(1:3,2,1), Mesh%Tnoeud(j)%Plocal(1:3,3,1)
end do
do j=Mesh%Body(Int_Body)%IndBody(2),Mesh%Body(Int_Body)%IndBody(4)
    write(3,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(Int_Body)%IndBody(1)-1)*[1,1,1]
end do 
close(unit=3)

!print*, 'Test Intégration sur le corps de GPhi2'
!Test=0._RP
!do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
!    Test = Test + ro*0.5_RP*dot_product(Ecoulement%GPhi(1:3,j)%perturbation,Ecoulement%GPhi(1:3,j)%perturbation)*Mesh%Tnoeud(j)%Aire*Mesh%Tnoeud(j)%Normale
!end do
!print*, Test
!
!Test=0._RP
!do j=Mesh%Body(Int_Body)%IndBody(2),Mesh%Body(Int_Body)%IndBody(4)
!    do k=1,3
!        PhiFacette(k)=Ecoulement%Phi(Mesh%Tfacette(j)%Tnoeud(k))%perturbation
!    end do
!    Grad = matmul(Mesh%Tfacette(j)%ds,PhiFacette)
!    Test = Test + ro*0.5_RP*dot_product(Grad,Grad)*Mesh%Tfacette(j)%Aire*Mesh%Tfacette(j)%Normale
!end do
!print*, Test
!deallocate(test)
! End
end subroutine PlotBodyGradient

subroutine TempWave(HouleObj)
! Parameter
Type(THoule) :: HouleObj
! Begin
HouleObj%Htype = 1;      HouleObj%nhoule = 1
allocate(HouleObj%w(HouleObj%nhoule), HouleObj%dir(HouleObj%nhoule), HouleObj%Aphi(HouleObj%nhoule), HouleObj%konde(HouleObj%nhoule), HouleObj%lambda(HouleObj%nhoule))
HouleObj%w(1) = 1._RP;     HouleObj%dir(1) = 90._RP/pi*0._RP;     HouleObj%Aphi(1) = 0.01_RP
call find_konde(HouleObj%konde,HouleObj%w,profondeur,HouleObj%NHoule)
HouleObj%lambda(1:NHoule) = 2._rp*Pi/HouleObj%konde(1:NHoule)
! End
end subroutine TempWave

end module Test_Gradient
