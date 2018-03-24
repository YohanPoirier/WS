!-------------------------------------------------------------
!
!	Remaillage des surfaces � partir des surfaces des cuves et de la SL incidente au temps t
!
!--------------------------------------------------------------
!
SUBROUTINE Remaillage(t, Maillage0, Maillage)

use Constantes
use StructuresDonnees
use FonctionsCommunes
!use Incident
use Houle_RF
implicit none

!Param�tres
!character(len=50) :: filename
type(TMaillage) :: Maillage0, Maillage
real(rp) :: t

! Variables Locales
character (len=50) :: filecoeff
integer :: j, k, r, li, lf, tf, mn, Tpresf, Interf
type(TMaillage) :: MaillageTemp
real(rp), dimension(3) :: M, Nor, ds1, ds2, A, B, GEta0
real(rp) :: Eta, delta
real(rp), dimension(3,3) :: Pg
real(rp) :: LongueurDOnde, AmplitudeCC, NombreDOnde, PeriodeH, VitessePhase, VitesseParticulaire, VitesseCourant, Direction, Profondeur
integer :: N
real(rp), dimension(NMax) :: An, Bn

call NewMaillage(MaillageTemp,PointMax) ! Cr�ation d'un maillage temporaire

MaillageTemp.Longueur=Maillage0.Longueur !
MaillageTemp.Largeur=Maillage0.Largeur   ! R�cup�ration des caract�ristiques du maillage
MaillageTemp.Hauteur=Maillage0.Hauteur   !
Direction = 0
Profondeur = Maillage0.Hauteur
filecoeff='waverf.cof'
call extract_WaveRF(filecoeff, profondeur, An, Bn, LongueurDOnde, AmplitudeCC, NombreDOnde, PeriodeH, VitessePhase, VitesseParticulaire, VitesseCourant, N)


li=0
do j=1,Maillage0.Nnoeud
    M=Maillage0%Tnoeud(j)%Pnoeud
    Eta=EtaIncidentRF(M, t, An, N, VitessePhase, NombreDOnde, Direction)
    if (Maillage0%Tnoeud(j)%typeNoeud==0) then  ! Si le point appartient � la Surface Libre, on le replace dessus : z=eta
        li=li+1
        MaillageTemp%Tnoeud(li)=Maillage0%Tnoeud(j)
        MaillageTemp%Tnoeud(li)%Pnoeud(3)=Eta
        do k=1,MaillageTemp%Tnoeud(li)%Nfacette ! Ajout des facettes connect�es � ce point au tableau des facettes et mise � jour des coordonn�es du point dans les facettes
!             Test(MaillageTemp%Tnoeud(li)%Tfacette(k,1))=1            
            MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%typeFrontiere=Maillage0%Tfacette(Maillage0%Tnoeud(j)%Tfacette(k,1))%typeFrontiere
            tf=MaillageTemp%Tnoeud(li)%Tfacette(k,2)
            MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Tnoeud(tf)=li
            MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Pnoeud(:,tf)=MaillageTemp%Tnoeud(li)%Pnoeud
        end do
    else ! Le point n'appartient pas � la surface libre
        if ((M(3).lt.Eta).and.(abs(M(3)-Eta).gt.0.05)) then ! Le point est au dessous de la SL et relativement loin (distance sup�rieure � 0.05)
            li=li+1
            MaillageTemp%Tnoeud(li)=Maillage0%Tnoeud(j)
            ! Modification de l'indice dans le tableau des noeuds dans les facettes correspondantes
            do k=1,MaillageTemp%Tnoeud(li)%Nfacette
                MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%typeFrontiere=Maillage0%Tfacette(Maillage0%Tnoeud(j)%Tfacette(k,1))%typeFrontiere
                tf = MaillageTemp%Tnoeud(li)%Tfacette(k,2)
                MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Tnoeud(tf)=li
                MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Pnoeud(:,tf)=MaillageTemp%Tnoeud(li)%Pnoeud
            end do
        else
            Tpresf=Interf(t, Maillage0, j); ! On v�rifie si la surface libre et les facettes connect�es au point d'indice j ont une intersection non nulle
            if ((Tpresf==1)) then ! Si c'est le cas, 
                M(3)=Eta        
                do r=1,li ! Test sur les points d�j� trait�s : Fusion des points 
                    if ((norm(M-MaillageTemp%Tnoeud(r)%Pnoeud).lt.0.05).and.(MaillageTemp%Tnoeud(r)%typeNoeud==Maillage0%Tnoeud(j)%typeNoeud).and.&
                    &(norm(Maillage0%Tfacette(MaillageTemp%Tnoeud(r)%Tfacette(1,1))%Normale-Maillage0%Tfacette(Maillage0%Tnoeud(j)%Tfacette(1,1))%Normale).lt.Epsilon)) then
                        do k=1,Maillage0%Tnoeud(j)%Nfacette
    !                        if (1.gt.0) then
                            MaillageTemp%Tfacette(Maillage0%Tnoeud(j)%Tfacette(k,1))%Tnoeud(Maillage0%Tnoeud(j)%Tfacette(k,2))=r
                            MaillageTemp%Tfacette(Maillage0%Tnoeud(j)%Tfacette(k,1))%Pnoeud(:,Maillage0%Tnoeud(j)%Tfacette(k,2))=M
    !                        end if
                        end do
                        go to 42
                    end if
                end do        
                ! Ne correspond � aucun point d�j� trait�, on l'ajoute � la suite
                li=li+1
                MaillageTemp%Tnoeud(li)=Maillage0%Tnoeud(j)
                MaillageTemp%Tnoeud(li)%Pnoeud(3)=Eta
                do k=1,MaillageTemp%Tnoeud(li)%Nfacette
    !             Test(MaillageTemp%Tnoeud(li)%Tfacette(k,1))=1
                    MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%typeFrontiere=Maillage0%Tfacette(Maillage0%Tnoeud(j)%Tfacette(k,1))%typeFrontiere
                    tf=MaillageTemp%Tnoeud(li)%Tfacette(k,2)
                    MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Tnoeud(tf)=li
                    MaillageTemp%Tfacette(MaillageTemp%Tnoeud(li)%Tfacette(k,1))%Pnoeud(:,tf)=MaillageTemp%Tnoeud(li)%Pnoeud
                end do
         42   end if
        end if
    end if
end do

MaillageTemp.Nnoeud=li
do j=1,MaillageTemp.Nnoeud
!    do k=j+1,MaillageTemp.Nnoeud
!        if ((norm(MaillageTemp%Tnoeud(j)%Pnoeud-MaillageTemp%Tnoeud(k)%Pnoeud).lt.Epsilon).and.(MaillageTemp%Tnoeud(j)%typeNoeud==MaillageTemp%Tnoeud(k)%typeNoeud)) then
!            print*, 'Points redondants: j = ', j, ', k = ', k            
!        end if
!    end do

    !    Calcul de la normale au point d'apr�s la variation locale (gradient surfacique) de la d�form�e de surface libre incidente
    if (MaillageTemp%Tnoeud(j)%typeNoeud.eq.0 ) then
        call GradEtaIncidentRF(M, t, An, N, VitessePhase, NombreDOnde, Profondeur, Direction, GEta0)
        A=[1._RP,0._RP,GEta0(1)]
        B=[0._RP,1._RP,GEta0(2)]    
        MaillageTemp%Tnoeud(j)%Normale=vect_product(A,B)/norm(vect_product(A,B))
    end if    
    !   R�cup�ration des points doubles et triples
    MaillageTemp%Tnoeud(j)%NDouble=0
    if (j.gt.1) then
        do k=1,j-1
            if (norm(MaillageTemp%Tnoeud(j)%Pnoeud-MaillageTemp%Tnoeud(k)%Pnoeud).lt.Epsilon) then
                MaillageTemp%Tnoeud(j)%NDouble=MaillageTemp%Tnoeud(j)%NDouble+1
                MaillageTemp%Tnoeud(k)%NDouble=MaillageTemp%Tnoeud(k)%NDouble+1
                MaillageTemp%Tnoeud(j)%Double(MaillageTemp%Tnoeud(j)%NDouble)=k
                MaillageTemp%Tnoeud(k)%Double(MaillageTemp%Tnoeud(k)%NDouble)=j
            end if
        end do
    end if
end do

MaillageTemp%Nfacette=Maillage0%Nfacette
Maillage.Nnoeud=MaillageTemp.Nnoeud
Maillage.Longueur=MaillageTemp.Longueur
Maillage.Largeur=MaillageTemp.Largeur
Maillage.Hauteur=MaillageTemp.Hauteur
do j=1,Maillage.Nnoeud
    Maillage%Tnoeud(j)=MaillageTemp%Tnoeud(j)
end do
! Traitement des facettes, pour supprimer les facettes de points non trait�s
lf=0
do j=1,MaillageTemp%Nfacette
    if (size(MaillageTemp%Tfacette(j)%Tnoeud(:),1).eq.3) then
!        print*, MaillageTemp%Tfacette(j)%Pnoeud(:,1), MaillageTemp%Tfacette(j)%Pnoeud(:,2), MaillageTemp%Tfacette(j)%Pnoeud(:,3)
        if(MaillageTemp%Tfacette(j)%Tnoeud(1).gt.0 .and. MaillageTemp%Tfacette(j)%Tnoeud(2).gt.0 .and. MaillageTemp%Tfacette(j)%Tnoeud(3).gt.0 .and.&
        ! V�rification que les points de la facette ne sont pas � la m�me position
        & norm(MaillageTemp%Tfacette(j)%Pnoeud(:,1)-MaillageTemp%Tfacette(j)%Pnoeud(:,2)).gt.0 .and. norm(MaillageTemp%Tfacette(j)%Pnoeud(:,1)-MaillageTemp%Tfacette(j)%Pnoeud(:,3)).gt.0 .and. norm(MaillageTemp%Tfacette(j)%Pnoeud(:,2)-MaillageTemp%Tfacette(j)%Pnoeud(:,3)).gt.0 ) then
            lf=lf+1
            Maillage%Tfacette(lf)=MaillageTemp%Tfacette(j)
            do k=1,3 ! Mise � jour de l'indice de la facette dans le tableau des facettes connect�es aux points
                if (MaillageTemp%Tfacette(j)%Tnoeud(k).lt.0) then
                    print*, 'erreur facettes : jfacette = ', lf, 'k = ', MaillageTemp%Tfacette(j)%Tnoeud, MaillageTemp%Tfacette(j)%Pnoeud
                    pause
                else
                    do mn=1,MaillageTemp%Tnoeud(Maillage%Tfacette(lf)%Tnoeud(k))%Nfacette
                        if (Maillage%Tnoeud(Maillage%Tfacette(lf)%Tnoeud(k))%Tfacette(mn,1)==j) then
                            Maillage%Tnoeud(Maillage%Tfacette(lf)%Tnoeud(k))%Tfacette(mn,1)=lf
                            Maillage%Tnoeud(Maillage%Tfacette(lf)%Tnoeud(k))%Tfacette(mn,2)=k
                        end if
                    end do
                end if
            end do
        end if
    end if
end do
Maillage%Nfacette=lf

do j=1,Maillage%Nfacette   
    if ((Maillage%Tfacette(j)%Tnoeud(1).eq.Maillage%Tfacette(j)%Tnoeud(2)).or.(Maillage%Tfacette(j)%Tnoeud(1).eq.Maillage%Tfacette(j)%Tnoeud(3)).or.(Maillage%Tfacette(j)%Tnoeud(2).eq.Maillage%Tfacette(j)%Tnoeud(3))) then
        print*, 'Probl�me de connectivit� sur la facette ', j
    else
        ! Calcul de la normale � la facette
        Nor=vect_product(Maillage%Tfacette(j)%Pnoeud(:,2)-Maillage%Tfacette(j)%Pnoeud(:,1),&
                           & Maillage%Tfacette(j)%Pnoeud(:,3)-Maillage%Tfacette(j)%Pnoeud(:,1))    
        if (norm(Nor)<0.0001) then
            print *, '        probl�me de calcul de la normale pour la facette', j, ' de coordonn�es', Maillage%Tfacette(j)%Pnoeud 
        else
            Maillage%Tfacette(j)%Normale=Nor/norm(Nor)
        end if
        ! Calcul du centre de gravit� de chaque facette, de son aire et de la distance maximale du centre aux sommets
        Pg=Maillage%Tfacette(j)%Pnoeud
        Maillage%Tfacette(j)%Gfacette=sum(Maillage%Tfacette(j)%Pnoeud,2)/3
        Maillage%Tfacette(j)%Aire=0.5*norm(vect_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1)))
        Maillage%Tfacette(j).Rmax=max(norm(Pg(:,1)-Maillage%Tfacette(j)%Gfacette),&
                                      &norm(Pg(:,2)-Maillage%Tfacette(j)%Gfacette),&
                                      &norm(Pg(:,3)-Maillage%Tfacette(j)%Gfacette))
        
       ! Calcul du gradient constant surfacique de la facette (ne d�pend que de la g�om�trie de la facette)    
        delta=(norm(Pg(:,2)-Pg(:,1))**2)*norm(Pg(:,3)-Pg(:,1))**2-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))**2
        ds1=(Pg(:,2)-Pg(:,1))*norm(Pg(:,3)-Pg(:,1))**2 - dot_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1))*(Pg(:,3)-Pg(:,1))
        ds2=(Pg(:,3)-Pg(:,1))*norm(Pg(:,2)-Pg(:,1))**2 - dot_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1))*(Pg(:,2)-Pg(:,1))
        Maillage%Tfacette(j)%ds=reshape([-(ds1+ds2),ds1,ds2]/delta,(/3,3/));    
    end if
end do



!! Calcul temporaire de la normale en un point d'apr�s les normales des facettes associ�es
!do j=1,Maillage.Nnoeud
!    Maillage%Tnoeud(j)%NormaleR(:) = 0
!    do k=1,Maillage%Tnoeud(j)%Nfacette
!        Maillage%Tnoeud(j)%NormaleR = Maillage%Tnoeud(j)%NormaleR + Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale
!    end do
!    Maillage%Tnoeud(j)%NormaleR = Maillage%Tnoeud(j)%NormaleR/norm(Maillage%Tnoeud(j)%NormaleR)
!end do

!do j=1,Maillage.Nnoeud
!    Maillage%Tnoeud(j)%NDouble=0
!    if (j.gt.1) then
!        do k=1,j-1
!            if (norm(Maillage%Tnoeud(j)%Pnoeud-Maillage%Tnoeud(k)%Pnoeud).lt.Epsilon) then
!                Maillage%Tnoeud(j)%NDouble=Maillage%Tnoeud(j)%NDouble+1
!                Maillage%Tnoeud(k)%NDouble=Maillage%Tnoeud(k)%NDouble+1
!                Maillage%Tnoeud(j)%Double(Maillage%Tnoeud(j)%NDouble)=k
!                Maillage%Tnoeud(k)%Double(Maillage%Tnoeud(k)%NDouble)=j
!            end if
!        end do
!    end if
!end do
!
!do j=1,Maillage.Nnoeud
!    if (Maillage%Tnoeud(j)%NDouble.eq.1) Maillage%NDouble=Maillage%NDouble+1
!    if (Maillage%Tnoeud(j)%NDouble.eq.2) Maillage.Ntriple=Maillage.Ntriple+1
!end do

!Maillage.NPSy=Maillage.Nnoeud-Maillage%NDouble-Maillage.Ntriple

!print*, Maillage.Nnoeud, Maillage%NDouble, Maillage.Ntriple
!print*, 'Np = ', Maillage.NPSy

end subroutine Remaillage


integer function Interf(t,Maillage,indice)

use Constantes
use StructuresDonnees
use FonctionsCommunes
use Incident

!Param�tre
type(TMaillage) :: Maillage
integer indice
real(rp) t

!Variables locales
integer :: j, k
real(rp) Eta
real(rp), dimension(3) :: M

Interf=0
do j=1,Maillage%Tnoeud(indice)%Nfacette
    do k=1,3
        M = Maillage%Tnoeud(Maillage%Tfacette(Maillage%Tnoeud(indice)%Tfacette(j,1))%Tnoeud(k))%Pnoeud
        Eta=eta0(M,t)
        if ((Eta-M(3)).gt.Epsilon2) then ! Valeur de 0.01 permet de ne pas consid�rer des points trop proches et d'�viter les triangles plats
            Interf=1
        end if
    end do
end do


end function
