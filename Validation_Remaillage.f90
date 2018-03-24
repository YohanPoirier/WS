subroutine Validation_Remaillage
use Parameters
use Constantes
use GenMaillage
use Structuresdonnees
use FonctionsCommunes
use Incident_mod
use BoucleTemp
use PrePlot
implicit none
! Variables Locales
character(len=50) :: fileout
type(TMaillage) :: Maillage0, Maillage, MaillageDT2
type(TEcoulement) :: Ecoulement0, Ecoulement, EcoulementTemp, EcoulementDT2
integer j, ios
real( rp) :: t, h, Phi0
real(rp), dimension(3) :: S1,S2
real(rp), allocatable :: RKp1(:), RKe1(:)
real(rp), allocatable :: CD(:,:), CS(:,:)
t=9
h=0.1

call Extract_Maillage(Maillage0)
!call NewEcoulement(Ecoulement0, Maillage0.Nnoeud)        
!call Incident(t, Ecoulement0, Maillage0)
call NewMaillage(Maillage, PointMax)
call Remaillage(t, Maillage0, Maillage)

fileout = 'RemaillageT_'//filename
call PlotMaill(fileout, Maillage)
    
call NewEcoulement(Ecoulement, Maillage%Nnoeud)        
call Incident(t, Maillage, Ecoulement)       
     
!   Solution Analytique des champs potentiel et vitesse normale
!    S1=[2.5, 0.25, 1]
!    S2=[2.5, 0.25,-3]
!    call PhiSS(Maillage, S1, S2, Ecoulement)

!   Initialisation de l'écoulement à partir du champ analytique
!call Ini(Maillage, Ecoulement)

allocate(CS(Maillage%Nnoeud,Maillage%Nnoeud), CD(Maillage%Nnoeud,Maillage%Nnoeud))
CS=-1
CD=-1
!   Resolution du Probleme aux Limites
print*, 'Resolution du probleme aux limites'
call solBVP(Ecoulement, Maillage, filename, CD, CS)
!   Récupération des données à afficher    
print*, 'Recuperation des sorties'
call outplot(filename, t, Maillage, Ecoulement)

do j=1,Maillage%Nnoeud
    Ecoulement%Eta(j)%perturbation = Ecoulement%Eta(j)%incident    
end do

!Calcul des Dérivées Temporelles et Spatiales de Phi(t) et Eta(t)
call Derive(t, Ecoulement, Maillage)

allocate(RKp1(Maillage%Nnoeud), RKe1(Maillage%Nnoeud))
!Calcul des pentes dérivées pour la première passe.
do j=1,Maillage%Nnoeud
    if (Maillage%Tnoeud(j)%TypeNoeud.eq.0) then
        RKp1(j)=Ecoulement%DPhiDt(j)%perturbation
        RKe1(j)=Ecoulement%DEtaDt(j)%perturbation
    end if
end do

!Calcul de Phi(t+h/2) et Eta(t+h/2) sur le MaillageT, à t
call NewEcoulement(EcoulementTemp,Maillage%Nnoeud)
call CopyEcoulement(EcoulementTemp, Ecoulement, Maillage%Nnoeud)

do j=1,Maillage%Nnoeud
    if (Maillage%Tnoeud(j)%TypeNoeud.eq.0) then
        EcoulementTemp.Phi(j)%perturbation = Ecoulement%Phi(j)%perturbation + RKp1(j)*h/2
        EcoulementTemp.Eta(j)%perturbation = Ecoulement%Eta(j)%perturbation + RKe1(j)*h/2
    end if
end do
call outplot(filename, t, Maillage, EcoulementTemp)
!Interpolation sur le MaillageDT2, à t+h/2
call NewMaillage(MaillageDT2,PointMax)
call Remaillage(t, Maillage0, MaillageDT2)

fileout = 'RemaillageDT2_'//filename
call PlotMaill(fileout, MaillageDT2)

call NewEcoulement(EcoulementDT2, MaillageDT2.Nnoeud)

call Interpolation(Maillage, EcoulementTemp, MaillageDT2, EcoulementDT2)

do j=1,MaillageDT2.Nnoeud
    if (MaillageDT2%Tnoeud(j)%TypeNoeud .ne. 0) then
        call CPhi0(MaillageDT2%Tnoeud(j)%Pnoeud,t+h/2, Phi0)
        EcoulementDT2.DPhiDn(j)%perturbation =  Phi0
    end if
end do
fileout = 'RK1_'//filename
call outplot(fileout, t, MaillageDT2, EcoulementDT2)

print*, '   RK 2eme passe'
!Résolution du Problème aux limites au temps t+h/2
CD=-1
CS=-1
call solBVP( EcoulementDT2, MaillageDT2, filename, CD, CS)

call Derive(t, EcoulementDT2, MaillageDT2)

fileout = 'RK2_'//filename
call outplot(fileout, t, MaillageDT2, EcoulementDT2)

end subroutine Validation_Remaillage