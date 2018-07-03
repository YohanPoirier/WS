!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                                              Validation du Solveur                               !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   On choisit une géométrie simple: une cuve parapipédique, un champ engendré par une singularité 
! extérieur au domaine. Discrétisation à panneaux constants
! Le potentiel sur la surface libre et la vitesse normale sur les parois engendrées par cette singularité
!constituent les entrées du solveur. En sortie, les solutions numériques: potentiel sur les parois et
!vitesse normale sur la surface libre, sont comparées aux valeurs analytiques exactes.
!
!
!27/05/2013 Lucas Letournel
!
!


!program Main
!use Constantes
!use Structuresdonnees
!use GenMaillage
!use PrePlot
!use Incident_mod
!use BVP
!use BVP_PnCst
!implicit none
!! Variables
!character (len=50) :: filename, filemaill
!logical :: oplot
!integer :: j, ios
!integer, parameter :: nt=3000
!type(TMaillage) :: Maillage
!type(TEcoulement) :: Ecoulement
!type(THouleRF) :: HouleRF
!real(rp), allocatable :: CD(:,:), CS(:,:)
!
!real(rp) :: L, a
!real(rp), dimension(101) :: x
!
!oplot=.false.
!filename='TestSolveur.dat'
!call NewMaillage(Maillage, PointMax)
!!call GenMaill(t(1), MaillageInit)
!call bodygen(filename, Maillage)
!!call mesh_cyl(Maillage)
!!filemaill='maill2.dat'
!!call Extract_Maillage(filemaill, Maillage)
!filemaill='Maillage.dat'
!call PlotMaill(filemaill, Maillage)
!
!allocate(CD(Maillage%Nfacette,Maillage%Nfacette),CS(Maillage%Nfacette,Maillage%Nfacette))
!call NewEcoulement(Ecoulement, Maillage%Nfacette)
!!Initialisation de l'écoulement
!call InitValSolveur_PnCst(Maillage, Ecoulement)
!CD(1,1)=-1
!call solBVP_PnCst( Ecoulement, Maillage, filename, CD, CS, 0._RP)
!call outplot_PnCSt(filename, Maillage, Ecoulement)
!
!pause
!end program Main

subroutine InitValSolveur_PnCst(Maillage, Ecoulement)
use Constantes
use StructuresDonnees
use FonctionsCommunes
implicit none
! Paramètres
real(rp) :: M1P, M2P
real(rp), dimension(3) :: M, S1, S2
Type(Tmaillage) :: Maillage
Type(TEcoulement) :: Ecoulement
! Variables Locales
integer :: j
! Début
S1=[0._RP, 0._RP, 1._RP]
S2=[0._RP, 0._RP, Maillage%Profondeur-1._RP]
do j=1,Maillage%Nfacette
    M = Maillage%Tfacette(j)%Gfacette
    M1P=1/norm2(M - S1)
    M2P=1/norm2(M - S2)
    Ecoulement%Phi(j)%incident = M1P + M2P        
    Ecoulement%DPhiDn(j)%incident = (M1P**3)*dot_product(S1-M,Maillage%Tfacette(j)%Normale) &
                                 & + (M2P**3)*dot_product(S2-M,Maillage%Tfacette(j)%Normale)
    if (Maillage%Tfacette(j)%TypeFrontiere == 0) then
        Ecoulement%Phi(j)%perturbation = Ecoulement%Phi(j)%incident
        Ecoulement%Eta(j)%perturbation = 0   
    else 
        Ecoulement%Eta(j)%perturbation = 0 
        Ecoulement%DPhiDn(j)%perturbation = Ecoulement%DPhiDn(j)%incident
    end if
end do
! Fin
end subroutine InitValSolveur_PnCst

module BVP_PnCst
use Constantes
use FonctionsCommunes
use StructuresDonnees
use BVP
implicit none
contains

subroutine solBVP_PnCst(Ecoulement, Maillage, filename, CD, CS, t, Option)
!!!!!Problème :
!   Résoudre le système des équations intégrales par
!       o Calcul ou récupération des coefficients d'influence
!       o Définition du système linéaire
!       o Résolution du système linéaire
!       o Redistribution des solutions
!!!!!
implicit none
!   Paramètres:
type(TEcoulement) :: Ecoulement
type(TMaillage), intent(in) ::  Maillage
character(len=50), intent(in) :: filename
real(rp), dimension(Maillage%Nfacette,Maillage%Nfacette):: CD, CS
real(rp), intent(in), optional :: t
logical, intent(in), optional :: Option !Pour le calcul de DPhiDt qur le corps
!   Variables Locales:         
logical:: CCI, Opt
integer :: j,k,r
real(rp) :: Solv, time_begin, time_end
real(rp),allocatable, dimension(:,:) ::  A
real(rp),allocatable, dimension(:) :: B, Sol      
 ! Début
allocate(A(Maillage%Nfacette,Maillage%Nfacette), B(Maillage%Nfacette), Sol(Maillage%Nfacette))    

call CoeffInfl_PnCst(Maillage, CD, CS)   

call SystLin_PnCst(CD, CS, Ecoulement, Maillage, A, B, filename, Option)
call cpu_time(time_begin)
Solv=0
if (Solv==0) then 
    call GMRES(Maillage, A, B, Sol, Maillage%Nfacette)           
else
    call LU(Maillage, A, B, Sol)  
end if  
call cpu_time(time_end)
print*, 'Le temps de la resolution du SL a ete de ', time_end - time_begin, ' secondes'
open(unit=11,file='temps.dat')
    write(11, '(a,e)') 'Temps de calcul', time_end - time_begin
    write(11, '(a,i)') 'Nombre de Noeud =', Maillage%Nnoeud
close(unit=11)
call postSL_PnCst(Sol, Maillage, Ecoulement, Option)
deallocate(A, B, Sol) !CS,CD,
! Fin
end subroutine solBVP_PnCst


subroutine CoeffInfl_PnCst(Maillage,CD,CS)
!!!!! Problème :
!   On calcule les coefficients d'influence pour un maillage donné. Pour chaque point du maillage, on
!   calcule l'influence de chaque facette. Par collocation on redistribue dans les matrices de coefficients.
!   Selon le critère de distance du point de contrôle à la facette, on utilise la solution exacte ou approchée.
!!!!! 
implicit none
!Paramètres
type(TMaillage), intent(in) :: Maillage
real(rp), dimension(Maillage%Nfacette,Maillage%Nfacette), intent(out)   :: CD,CS
!Locales
real(rp)  :: Ssigma, Smu
real(rp), dimension(3)    :: Isigma, Imu, Ar
real(rp), dimension(3)      :: Css, Cdd
integer :: j,k, r, test, testp
real(rp) ::  Rho, Cr, Z
real(rp), dimension(3) ::   MG
real(rp), dimension(3,3) :: ds
CD(:,:)=0
CS(:,:)=0
do j=1,Maillage%Nfacette
    do k=1,Maillage%Nfacette
!        if (j.ne.k) then
            MG=Maillage%Tfacette(j)%Gfacette - Maillage%Tfacette(k)%Gfacette
            Rho=sqrt(MG(1)**2 + MG(2)**2 + MG(3)**2 ) 
            Cr=Rho/Maillage%Tfacette(k)%Rmax
            if (Cr<=18) then !
    !       Calcul Exact
                call Guevel(Maillage%Tfacette(j)%Gfacette, Maillage%Tfacette(k),Ssigma,Smu)
            else
    !       Calcul Asymptotique 
                Ssigma=Maillage%Tfacette(k)%Aire/Rho
                Z= dot_product(MG, Maillage%Tfacette(k)%Normale)
                Smu=Z*Maillage%Tfacette(k)%Aire/(Rho**3)
            end if
            CS(j,k) = Ssigma
            CD(j,k) = Smu
!        end if
    end do
end do
do j=1,Maillage%Nnoeud
    CD(j,j)=0   ! Evite de reprendre une valeur de l'angle solide éventuellement déjà calculé.
    CD(j,j)= -sum(CD(j,:))
end do
end subroutine CoeffInfl_PnCst

subroutine SystLin_PnCst(CD, CS, Ecoulement, Maillage, A, B, filename, Option)
!!!!! Problème :
!   Initialiser le système linéaire à partir du type de frontière
!       o Surface Libre: On cherche la vitesse normale --> A = CS et B = B + CD*Phi
!       o Surfaces Matérielles: On cherche le potentiel --> A= -CD et B = B - CS*DPhiDn
!!!!!
implicit none
!   Paramètres
character(len=50), optional, intent(in) :: filename
type(TMaillage), intent(in)     :: Maillage
type(TEcoulement) , intent(in)  :: Ecoulement
real(rp), dimension(:,:), intent(in)   :: CD, CS
real(rp), dimension(:,:), intent(out)   :: A
real(rp), dimension(:), intent(out)   ::  B
logical, intent(in), optional :: Option
!   Variables Locales
logical :: Opt
integer :: j, k
character(len=50) :: filesyst, fileInfl
real(rp) :: Diag
!   Début
do j=1,Maillage%Nfacette
    B(j)=0._RP
    do k=1,Maillage%Nfacette
        if (Maillage%Tfacette(k)%TypeFrontiere.eq.0) then
            A(j,k) = CS(j,k)                
            B(j) = B(j) + CD(j,k)*Ecoulement%Phi(k)%perturbation
        else 
            A(j,k) = - CD(j,k)
            B(j) = B(j) - CS(j,k)*Ecoulement%DPhiDn(k)%perturbation
        end if
    end do
end do   

!    Preconditionnement de la matrice du système linéaire
do j=1,Maillage%Nfacette
	Diag=A(j,j)
	if (abs(Diag).GT.Epsilon) then
		B(j)=B(j)/Diag
		A(j,:)=A(j,:)/Diag
	else
		print*, ' Le coefficient ',j,' de la diagonale de A est nul !'
    end if
end do

if (.false.) then
    fileInfl='Matrix_'//filename
    open(3,file=fileInfl)
    write(3,'(50a)') 'Title= "Matrice du Système Linéaire"'
    write(3,'(a,i,a,i)') 'Nombre de Noeuds =', Maillage%Nnoeud, ', Nombre de Facettes=', Maillage%Nfacette

    do j=1,Maillage%Nfacette
        do k=1,Maillage%Nfacette
            write(3,*), A(j,k)
        end do
    end do
    do j=1,Maillage%Nfacette
        write(3,*), B(j) 
    end do   
end if
close(3)   
    
end subroutine SystLin_PnCst


subroutine PostSL_PnCst(X, Maillage, Ecoulement, Option)
!!!!! Problème :
!   Travail inverse de la fonction Ini, on réorganise les solutions à leur place dans Ecoulement   
!       o Surface Libre : DPhiDn = solution
!       o Surface Matérielle : Phi = solution   
!!!!!
implicit none
!   Paramètres
type(TMaillage), intent(in) ::  Maillage
Real(RP), dimension(:), intent(in) :: X
type(TEcoulement), intent(inout)   :: Ecoulement
logical, intent(in), optional :: Option
!   Variables Locales
integer :: j, k
logical :: Opt
! Début
do j=1,Maillage%Nfacette
    if (Maillage%Tfacette(j)%TypeFrontiere.eq.0) then
        Ecoulement%DPhiDn(j)%perturbation = X(j)
    else
        Ecoulement%Phi(j)%perturbation = X(j)
    end if            
end do
!   Fin
return
end subroutine PostSL_PnCst

subroutine outplot_PnCst(filename, Maillage, Ecoulement)
!!!!! Problème :
!   Récupérer les données à afficher:
!       o Ecoulement : solution du problème aux limites 
!       o EcoulementSS : solution analytique déterminée dans PhiSS
!!!!!
implicit none
character(len=50) :: filename,filemaillage,fileout,filemax
!	Maillage de la géométrie
type(TMaillage)			:: Maillage
!
!	Tables potentiel et de son gradient normal aux noeuds du maillage
type(TEcoulement)		:: Ecoulement
! Variables Locales
integer :: j,k,l,iMin,Nresult
integer, dimension(1) :: MaxLo, MinLo
real(rp) :: Min, MinErr, MaxErr, MaxAn, MinAn
real(rp):: IntN, IntA, Eps, IntE
!real(rp), dimension(Maillage%Nnoeud,6) :: results, resultss
!real(rp), dimension(Maillage%Nnoeud) :: pert
real(rp), dimension(:,:), allocatable :: results, resultss
real(rp), dimension(:), allocatable :: pert

allocate(results(Maillage%Nfacette,7), resultss(Maillage%Nfacette,7), pert(Maillage%Nfacette))
!	Initialisation des tableaux
results=0
resultss=0	
!   Récupération des données dans Ecoulement
l=1
if (Maillage%typeM.eq.0) then
    do j=1,Maillage%Nfacette
        if (Maillage%Tfacette(j)%Npanneau.eq.0 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2) then
            if (Maillage%Tfacette(j)%Gfacette(1).lt.0) then
                results(l,1)=-Maillage%LL - Maillage%Profondeur + Maillage%Tfacette(j)%Gfacette(1) ! compris entre 0 et Maillage%L	        
            else
                results(l,1)= Maillage%LL + Maillage%Profondeur + Maillage%Tfacette(j)%Gfacette(1)
            end if
            results(l,2)=Ecoulement%DPhiDn(j)%perturbation
            results(l,3)=Ecoulement%DPhiDn(j)%incident
            results(l,4)=Ecoulement%Phi(j)%perturbation
            results(l,5)=Ecoulement%Phi(j)%incident
            results(l,6)=Maillage%Tfacette(j)%Npanneau
            results(l,7)=results(l,3)-results(l,2)
            l=l+1
        elseif (Maillage%Tfacette(j)%Npanneau.eq.1 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2) then
            if (Maillage%Tfacette(j)%Gfacette(1).lt.0) then
                results(l,1) = -Maillage%LL - Maillage%Profondeur - Maillage%Tfacette(j)%Gfacette(3) ! compris entre 0 et Maillage%L	        
            else
                results(l,1) =  Maillage%LL + Maillage%Profondeur + Maillage%Tfacette(j)%Gfacette(3)
            end if
            results(l,2)=Ecoulement%Phi(j)%perturbation
            results(l,3)=Ecoulement%Phi(j)%incident
            results(l,4)=Ecoulement%DPhiDn(j)%perturbation
            results(l,5)=Ecoulement%DPhiDn(j)%incident
            results(l,6)=Maillage%Tfacette(j)%Npanneau
            results(l,7)=results(l,3)-results(l,2)
            l=l+1
        elseif (Maillage%Tfacette(j)%Npanneau.eq.2 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2) then
            results(l,1)=Maillage%Tfacette(j)%Gfacette(1) ! compris entre -Maillage%Profondeur et 0
            results(l,2)=Ecoulement%Phi(j)%perturbation
            results(l,3)=Ecoulement%Phi(j)%incident
            results(l,4)=Ecoulement%DPhiDn(j)%perturbation
            results(l,5)=Ecoulement%DPhiDn(j)%incident
            results(l,6)=Maillage%Tfacette(j)%Npanneau
            results(l,7)=results(l,3)-results(l,2)
            l=l+1
        end if
    end do
else
    do j=1,Maillage%Nfacette ! On teste tous les noeuds
    !        Le noeud appartient à la surface libre et au plan de symétrie
    !	    if (Maillage%Tnoeud(j)%TypeNoeud.eq.10 .and. abs(Maillage%Tnoeud(j)%Pnoeud(2)-0*Maillage%LL/2).lt.Epsilon2 .and. abs(Maillage%Tnoeud(j)%Pnoeud(1)).gt.Epsilon2 .and. Maillage%Tnoeud(j)%Pnoeud(1).lt.Maillage%L ) then !.and. Maillage%Tnoeud(j)%Pnoeud(3)==0
        if (Maillage%Tfacette(j)%Npanneau.eq.1 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2 ) then !.and. Maillage%Tnoeud(j)%Pnoeud(3)==0
            results(l,1)=Maillage%L/2 + Maillage%Tfacette(j)%Gfacette(1) ! compris entre -L/2 et L/2	        
            results(l,2)=Ecoulement%DPhiDn(j)%perturbation
            results(l,3)=Ecoulement%DPhiDn(j)%incident
            results(l,4)=Ecoulement%Phi(j)%perturbation
            results(l,5)=Ecoulement%Phi(j)%incident
            results(l,6)=Maillage%Tfacette(j)%Npanneau
            results(l,7)=results(l,3)-results(l,2)
            l=l+1
    !	    Le noeud appartient à la paroi latérale gauche et au plan de symétrie    
!        else if (Maillage%Tfacette(j)%Npanneau.eq.5 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2) then 
!            results(l,1)=Maillage%Tfacette(j)%Gfacette(3) ! compris entre -Maillage%Profondeur et 0
!            results(l,2)=Ecoulement%Phi(j)%perturbation
!            results(l,3)=Ecoulement%Phi(j)%incident
!            results(l,4)=Ecoulement%DPhiDn(j)%perturbation
!            results(l,5)=Ecoulement%DPhiDn(j)%incident
!            results(l,6)=Maillage%Tfacette(j)%Npanneau
!            results(l,7)=results(l,3)-results(l,2)
!            l=l+1
!    !	    Le noeud appartient à la paroi latérale droite et au plan de symétrie    
!        else if (Maillage%Tfacette(j)%Npanneau.eq.6 .and. abs(Maillage%Tfacette(j)%Gfacette(2)).lt.Epsilon2) then 
!            results(l,1)=Maillage%L - Maillage%Tfacette(j)%Gfacette(3) ! compris entre Maillage%L+1 et Maillage%L+1+Maillage%Profondeur
!            results(l,2)=Ecoulement%Phi(j)%perturbation
!            results(l,3)=Ecoulement%Phi(j)%incident
!            results(l,4)=Ecoulement%DPhiDn(j)%perturbation
!            results(l,5)=Ecoulement%DPhiDn(j)%incident
!            results(l,6)=Maillage%Tfacette(j)%Npanneau
!            results(l,7)=results(l,3)-results(l,2)
!            l=l+1
        endif
    end do
end if
Nresult=l-1
!   Réorganisation des résultats selon l'"abscisse" du noeud considéré (results(k,1))
do j=1,Nresult
    Min=100
    ! Recherche du minimum dans la première colonne 
    do k=1,Nresult
        if (results(k,1).LE.Min) then
            Min=results(k,1)
            iMin=k
        end if
    end do  
    ! Réagencement des lignes selon le minimum précédemment déterminé
    resultss(j,:)=results(iMin,:)
    results(iMin,1)=1000
end do
    
!   Création du fichier de sortie des données pour Tecplot
fileout='outEcoulement_'//filename
open(3,file=fileout)
write(3,fmt='(50a)') 'Title= "Maillage de la cuve"'
write(3,fmt='(50a)') 'VARIABLES ="X","Solution Numérique","Solution Analytique","CL Num","CL Ana","typeNoeud","Erreur Relative"'	 
do j=1,Nresult
    write(3,'(7E)') resultss(j,:)
!	    write(3,'(6E)') resultss(j,1), resultss(j,2), resultss(j,3), resultss(j,4), resultss(j,5), resultss(j,6)
end do
close(3)    

!   Calcul de l'erreur relative entre les deux solutions intégrées 
IntN=0
IntA=0 
IntE=0   
!do j1=,Nresult-1
!    if (resultss(j+1,6).eq.1) then
!        ! Intégration par la méthode des trapèzes
!        IntN=IntN+abs(resultss(j+1,2)+resultss(j,2))*(resultss(j+1,1)-resultss(j,1))/2
!        IntA=IntA+abs(resultss(j+1,3)+resultss(j,3))*(resultss(j+1,1)-resultss(j,1))/2
!        IntE=IntE+abs(resultss(j+1,2)+resultss(j,2)-resultss(j+1,3)-resultss(j,3))*(resultss(j+1,1)-resultss(j,1))/2
!    end if
!end do
do j=1,Nresult-2
    if (resultss(j+1,6).eq.1) then
        ! Intégration par la méthode de Simpson
!        IntN=IntN+abs(resultss(j+1,2)+resultss(j,2))*(resultss(j+1,1)-resultss(j,1))/2
        IntA=IntA+abs(resultss(j,3)+4*resultss(j+1,3)+resultss(j+2,3))*(resultss(j+2,1)-resultss(j,1))/6
        IntE=IntE+abs(resultss(j,2)+4*resultss(j+1,2)+resultss(j+2,2)-(resultss(j,3)+4*resultss(j+1,3)+resultss(j+2,3)))*(resultss(j+2,1)-resultss(j,1))/6
    end if
end do
!    print*, IntN, IntA
Eps=100*IntE/IntA
print*, 'Eps = ', Eps
print*, 'IntE = ', IntE
print*, 'IntA = ', IntA
!   Création du fichier de sortie pour l'erreur relative
filemax='Eps_'//filename
open(unit=11,file=filemax)
    write(11, fmt='(50a)') filename
    write(11, '(a,e)') 'IntN = ', IntE, 'IntA = ', IntA, 'Eps = ', Eps
    write(11, '(a,i)') 'Nombre de Noeud =', Maillage%Nnoeud
close(unit=11)

deallocate(results, resultss, pert)
end subroutine outplot_PnCst

end module BVP_PnCst