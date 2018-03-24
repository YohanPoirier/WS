module Validation_IntTemp
use Parameters
use Constantes
use FonctionsCommunes
use Structuresdonnees
implicit none
!real(rp) :: Lambda = Ldom(1)/2._RP
!real(rp) :: Klambda = 2._RP*pi/Lambda
real(rp) :: A_batteur = 0.05_RP
!real(rp) :: Omeg = sqrt(g*Klambda*tanh(Klambda*Ldom(3)))   
    contains
    
    
subroutine Init_OndeStat(Mesh, Ecoulement, t)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
! Variables
integer :: j
real(rp) :: Amplitude, Lambda, Klambda, Omeg
! Begin
Lambda = Mesh%DimTank(1)
Klambda = 2*pi/Lambda
Amplitude = 0.1*Lambda
Omeg = sqrt(g*Klambda*tanh(Klambda*Mesh%DimTank(3)))
do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
    !Ecoulement%Eta(j)%Perturbation = Amplitude*cos(2*Pi/Lambda*norm(Mesh%Tnoeud(j)%Pnoeud))
    Ecoulement%Eta(j)%Perturbation = Amplitude*cos(2*Pi/Lambda*Mesh%Tnoeud(j)%Pnoeud(1))*cos(Omeg*t)
end do
! End
end subroutine Init_OndeStat

subroutine Plot_OndeStat(Mesh, Ecoulement, t)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
! Variables
character(len=10) :: num
integer :: j
real(rp) :: Maxi
real(rp), dimension(Mesh%Nnoeud) :: EpsOS
type(TEcoulement) :: EcoulementA
! Begin
Write( num, '( f0.4 )' ) t
open(27,file='OndeStat_'//filename, position='append')
write(27,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)

call NewEcoulement(EcoulementA, Mesh%Nnoeud)
call IniEcoulement(EcoulementA, Mesh%Nnoeud, 0._RP)
call Init_OndeStat(Mesh, EcoulementA, t)
maxi = 0._RP
EpsOS = 0._RP
do j=1,Mesh%Nnoeud
    maxi = max(maxi,EcoulementA%Eta(j)%Perturbation)
end do
if (abs(maxi).lt.Epsilon) Maxi = 1._RP
do j=1,Mesh%Nnoeud
    EpsOS(j) = 100._RP*abs(Ecoulement%Eta(j)%Perturbation - EcoulementA%Eta(j)%Perturbation)/maxi
    write(27,'(6E)') Mesh%Tnoeud(j)%Pnoeud(1:2), Mesh%Tnoeud(j)%Pnoeud(3)+Ecoulement%Eta(j)%Perturbation, EpsOS(j), Ecoulement%Eta(j)%Perturbation, EcoulementA%Eta(j)%Perturbation
end do
do j=1,Mesh%Nfacette
    write(27,'(3I)') Mesh%Tfacette(j)%TNoeud
end do
close(27)
call DelEcoulement(EcoulementA)
! End
end subroutine Plot_OndeStat


subroutine Batteur(Mesh, Ecoulement, t)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
! Variables
integer :: j, nc
real(rp) :: Amplitude, Lambda, Klambda, Omeg
! Begin
Lambda = Mesh%DimTank(1)/5._RP
Klambda = 2*pi/Lambda
Amplitude = A_batteur*Lambda
Omeg = sqrt(g*Klambda*tanh(Klambda*Mesh%DimTank(3)))

do nc=1,Mesh%NBody
    if (nc.ne.Int_body) then
        do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
            if (abs(Mesh%Tnoeud(j)%Normale(1)-1._RP).lt.Epsilon) then
                Ecoulement%DPhiDn(j)%perturbation = Amplitude*sin(Omeg*t)
            end if             
        end do
    end if
end do
!End
end subroutine Batteur

subroutine EnergieVol_Batteur(Mesh, EcoulementTemp, EcoulementT, Integre, DPhiDt, t, dt)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: EcoulementT, EcoulementTemp
real(rp) :: t, dt
real(rp), dimension(:) :: DPhiDt
real(rp), dimension(2) :: Integre
! Variables
integer :: j, nc
! Begin
Integre = 0._RP
do nc=1,Mesh%Nbody
    do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
        !if (abs(dt).lt.Epsilon) then
        !    DPhiDt = EcoulementT%DPhiDt(j)%perturbation
        !    !EcoulementTemp%DPhiDt(j)%perturbation = EcoulementT%DPhiDt(j)%perturbation
        !else
        !    DPhiDt(j) = (EcoulementTemp%Phi(j)%perturbation - EcoulementT%Phi(j)%perturbation)/dt
        !    !EcoulementTemp%DPhiDt(j)%perturbation = (EcoulementTemp%Phi(j)%perturbation - EcoulementT%Phi(j)%perturbation)/dt
        !end if
        !Integre(2) = Integre(2) - ro*Mesh%Tnoeud(j)%Aire*EcoulementTemp%DPhiDn(j)%perturbation*DPhiDt(j)
        Integre(1) = Integre(1) - ro*Mesh%Tnoeud(j)%Aire*EcoulementTemp%DPhiDn(j)%perturbation*EcoulementTemp%DPhiDt(j)%perturbation
    end do
end do
call Vol_Batteur(Integre(2), Mesh, t)
! End
end subroutine EnergieVol_Batteur

subroutine InitBVP_Batteur(Mesh,Ecoulement,t)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
! Locals
integer :: j, nc
real(rp) :: Amplitude, Lambda, Klambda, Omeg
! Begin
Lambda = Mesh%DimTank(1)/5._RP
Klambda = 2*pi/Lambda
Amplitude = A_batteur*Lambda
Omeg = sqrt(g*Klambda*tanh(Klambda*Mesh%DimTank(3)))
!Amplitude = 1
!Omeg = w(1)
do nc=1,Mesh%Nbody
    do j=Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)
        if (abs(Mesh%Tnoeud(j)%Normale(1)-1._RP).lt.Epsilon) then
            Ecoulement%DDPhiDnDt(j)%perturbation = Amplitude*Omeg*cos(Omeg*t)
        else
            Ecoulement%DDPhiDnDt(j)%perturbation = 0._RP
        end if
    end do
end do
! End
end subroutine InitBVP_Batteur

subroutine Vol_Batteur(Vol, Mesh, t)
! Parameters
real(rp) :: Vol, t
type(TMaillage) :: Mesh
! Locals
real(rp) :: Amplitude, Lambda, Klambda, Omeg
! Begin
Lambda = Mesh%DimTank(1)/5._RP
Klambda = 2*pi/Lambda
Amplitude = A_batteur*Lambda
Omeg = sqrt(g*Klambda*tanh(Klambda*Mesh%DimTank(3)))
!Amplitude = 1
!Omeg = w(1)
Vol = Ldom(2)*LDom(3)*Amplitude*sin(Omeg*t)
if (symmetry) Vol =2._RP*Vol
! End
end subroutine Vol_Batteur

subroutine PrePlot_Batteur
! Parameters
! Locals
integer :: ios
! Begin
open(unit=101,file='Domm_'//filename, iostat=ios)
if (ios/=0) stop "Erreur � l'ouverture du fichier de comparaison Lin�aire/Non-Lin�aire"
write(101,fmt='(50a)') 'Title = "Comparaison Lin�aire/Non-Lin�aire"'
write(101,fmt='(50a)') 'VARIABLES = "X","<greek>e</greek><sub><greek>h</greek></sub>","<greek>h</greek><sub>Num</sub>","<greek>h</greek><sub>Ana</sub>"'
close(101)
open(unit=27, file='Batteur_'//filename)
write(27,fmt='(50a)') 'Title = "Batteur 3D"'
write(27,fmt='(50a)') 'VARIABLES = "X","Y","Z","<greek>e</greek><sub><greek>h</greek></sub>","<greek>h</greek><sub>Num</sub>","<greek>h</greek><sub>Ana</sub>"'
close(27)
! End
end subroutine PrePlot_Batteur

subroutine Plot_Batteur(Mesh, Ecoulement, Nt, t, jt)
! Parameters
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
integer :: jt, Nt
real(rp), dimension(:) :: t
! Variables
character(len=10) :: num
integer :: j
real(rp) :: Maxi
real(rp), dimension(Mesh%Nnoeud) :: EpsOS
real(rp), dimension(PointMax) :: EtaDomm
! Begin
Write( num, '( f0.4 )' ) t(jt)*w(1)/(2._RP*pi)
open(101,file='Domm_'//filename, position='append')
write(101,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", SOLUTIONTIME = '//trim(num) !, N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
call Dommermuth88(nt, t, jt, Mesh, EtaDomm)
EpsOS = 0._RP
maxi = 0._RP
do j=1,Mesh%Nnoeud
    maxi = max(maxi,EtaDomm(j))
end do
if (abs(maxi).lt.Epsilon) Maxi = 1._RP
do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
    EpsOS(j) = 100._RP*abs(Ecoulement%Eta(j)%Perturbation - EtaDomm(j))/maxi
    if (abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2) then
        write(101,'(4E)') Mesh%Tnoeud(j)%Pnoeud(1), EpsOS(j), Ecoulement%Eta(j)%Perturbation, EtaDomm(j) !, ECc(NErr), EPp(NErr)
    end if
end do
close(101)

open(27,file='Batteur_'//filename, position='append')
write(27,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
do j=1,Mesh%Nnoeud
    write(27,'(6E)') Mesh%Tnoeud(j)%Pnoeud(1:2), Mesh%Tnoeud(j)%Pnoeud(3)+Ecoulement%Eta(j)%Perturbation, EpsOS(j), Ecoulement%Eta(j)%Perturbation, EtaDomm(j) !, EcoulementA%Eta(j)%Perturbation
end do
do j=1,Mesh%Nfacette
    write(27,'(3I)') Mesh%Tfacette(j)%TNoeud
end do
close(27)
!
!call DelEcoulement(EcoulementA)
! End
end subroutine Plot_Batteur

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Dommermuth88(nt, t, jt, Mesh, EtaDomm)
! Parametres
integer, intent(in) :: nt, jt
real(rp), dimension(nt), intent(in) :: t
type(Tmaillage), intent(in) :: Mesh
real(rp), dimension(PointMax), intent(out) :: EtaDomm
! Variables Locales
integer :: j,k,r
integer, parameter :: Nn=1000
real(rp) :: L, B, Xmod
real(rp), dimension(Nn) :: A, Z
! D�but
L=Mesh%DimTank(1)
!print*, 'A REPRENDRE !!! '  !! Dimensions L et LL
call PreDomm(nt, t, jt, L, A, B, Z, Nn)
EtaDomm = B/L
do j=Mesh%Fs%IndFS(1),Mesh%Fs%IndFS(3)
    !if (abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2) then
        do k=1,Nn
            Xmod = Mesh%Tnoeud(j)%Pnoeud(1) + L/2._RP
            EtaDomm(j) = EtaDomm(j) + 4/L*Z(k)*cos(k*pi/L*Xmod)*A(k)
        end do        
    !end if
end do
! Fin
end subroutine Dommermuth88

subroutine PreDomm(nt, t, jt, L,  A, B, Z, Nn)
! Parametres
integer, intent(in) :: nt, jt, Nn
real(rp), dimension(nt), intent(in) :: t
real(rp), intent(in) :: L
real(rp), dimension(Nn), intent(out) :: A, Z
real(rp), intent(out) :: B
! Variables Locales
integer j,k,m
real(rp) :: kn, wn, lm
real(rp), dimension(jt) :: U, V
real(rp) :: Amplitude, Lambda, Klambda, Omeg
! Begin
Lambda = LDom(1)/5._RP
Klambda = 2*pi/Lambda
Amplitude = 0.05*Lambda
Omeg = sqrt(g*Klambda*tanh(Klambda*LDom(3)))
!Amplitude = 1
!Omeg = w(1)

do j=1,jt
    !U(j) = sin(2._RP*sqrt(g/Ldom(3))*t(j))
    U(j) = Amplitude*sin(Omeg*t(j))
end do
B=0
!do j=1,jt-1
!    ! Int�gration par la m�thode des trap�zes
!    B=B+(U(j+1)+U(j))*(t(j+1)-t(j))/2
!end do
do j=1,jt-2,2
    ! Int�gration par la m�thode de Simpson
    B = B + (U(j) + 4*U(j+1) + U(j+2))*(t(j+2)-t(j))/6
end do
A=0
Z=0
do k=1,Nn
    kn = k*pi/L
    do m=0,Nn
        lm = (m+0.5_RP)*pi
        Z(k) = Z(k) + 1/(kn**2 + lm**2)
    end do
    wn = sqrt(g*kn*tanh(kn*LDom(3)))
    do j=1,jt
        V(j) = U(j)*cos(wn*(t(jt)-t(j)))
    end do
!    do j=1,jt-1
!        ! Int�gration par la m�thode des trap�zes
!        A(k)=A(k)+(V(j+1)+V(j))*(t(j+1)-t(j))/2
!    end do   
    do j=1,jt-2,2
        ! Int�gration par la m�thode de Simpson
        A(k) = A(k) + (V(j) + 4*V(j+1) + V(j+2))*(t(j+2)-t(j))/6
    end do
end do
! Fin
end subroutine PreDomm

end module Validation_IntTemp