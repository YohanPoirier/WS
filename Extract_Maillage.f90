!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                               Extraction d'un maillage depuis un fichier                         !
!                                           Lucas 14/06/12                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine extract_maillage(filename, Maillage)
!!!!! Probl�me :
!   Extraire le maillage d'un fichier texte, contenant :
!       o Le nom du maillage sur la premi�re ligne (100 caract�re maxi)
!       o Le nombre de points et de facettes
!       o Les dimensions de la cuve: (optionnel)
!           + Longueur
!           + Largeur
!           + Hauteur
!       o Indice puis coordonn�s du point 
!       o Indice de la facette puis indice des sommets puis type de fronti�re (0--> SL, autre --> SM)
!Ex:
!mailllarge
!Nombre de Points =    558	 Nombre de Facettes =    864	 
!Dimensions de la cuve : 
!L = 5.000000
!l = 0.500000
!h = 2.000000
!1	0.000000E+00	0.000000E+00	0.000000E+00	
!2	0.000000E+00	2.500000E-01	0.000000E+00	
!3	0.000000E+00	5.000000E-01	0.000000E+00
!...
!1	5	2	1	0	
!2	1	4	5	0	
!3	6	3	2	0	
!4	2	5	6	0
!!!!!

use Constantes
use StructuresDonnees
use FonctionsCommunes

implicit none

!Param�tres
character(len=50) :: filename
type(TMaillage) :: Maillage

!Variables Locales
integer ::  j, k, r, m, np, nf, t, Nnoeud, Nfacette, ios
real(rp) :: delta
real(rp), dimension(3) :: nor, ds1, ds2, te
real(rp), dimension(4) :: lignep
integer, dimension(5)  :: lignef
!integer, dimension(3)  :: te
real(rp), dimension(3,3):: Pg
character(len=100) :: truc
character(len=19) :: muche
character(len=23) :: machin,forma

! Cr�ation d'un objet Maillage
call NewMaillage(Maillage,Pointmax)

! Ouverture du fichier de maillage
open(10,file=filename, iostat=ios)
if (ios/=0) stop "Erreur � l'ouverture du fichier de maillage"

! Lecture des informations du maillage
read(10,fmt='(A100)') truc
print *, truc
!if (filename=='cuvereduite.dat') then
!    forma='(A19,I3,A23,I4)'
!else
!    forma='(A19,I4,A23,I4)'
!end if
forma='(A19,I6,A23,I6)'
read(10,fmt=forma) muche, Nnoeud, machin, Nfacette
!60 format(A100,I4,A100,I4)
print*,  muche, Nnoeud, machin, Nfacette
Maillage.Nnoeud=Nnoeud
Maillage%Nfacette=Nfacette
read(10,'(a)') truc
read(10,fmt='(A3,f8.7)') muche, Maillage.Longueur
read(10,fmt='(A3,F8.7)') muche, Maillage.Largeur
read(10,fmt='(A3,F8.7)') muche, Maillage.Hauteur
print*, 'L = ', Maillage.Longueur, ', l = ', Maillage.Largeur, ', h = ', Maillage.Hauteur

! Lecture des noeuds du maillage
do j=1,Maillage.Nnoeud
    read(10,*) lignep
    np=lignep(1)
    Maillage%Tnoeud(j)%Pnoeud=[lignep(2), lignep(3), lignep(4)]    
end do

! Lecture des facettes du maillage
do j=1,Maillage%Nfacette
    read(10,*) lignef
    nf=lignef(1)
    ! Indices des sommets de la facette
    Maillage%Tfacette(j)%Tnoeud=[lignef(2), lignef(3), lignef(4)]
    ! Type de fronti�re de la facette
    Maillage%Tfacette(j)%typeFrontiere=lignef(5)
    ! R�cup�ration des coordonn�es des sommets de la facette
    Maillage%Tfacette(j)%Pnoeud=reshape([Maillage%Tnoeud(lignef(2))%Pnoeud,Maillage%Tnoeud(lignef(3))%Pnoeud,&
                                    &Maillage%Tnoeud(lignef(4))%Pnoeud],(/3,3/))
    Nor=vect_product(Maillage%Tfacette(j)%Pnoeud(:,2)-Maillage%Tfacette(j)%Pnoeud(:,1),&
                       & Maillage%Tfacette(j)%Pnoeud(:,3)-Maillage%Tfacette(j)%Pnoeud(:,1))
    ! Calcul de la normale � la facette
    if (norm(Nor)<0.0001) print *, '        probleme de calcul de la normale pour la facette', j, ' de coordonn�es', Maillage%Tfacette(j)%Pnoeud
      
    Maillage%Tfacette(j)%Normale=Nor/norm(Nor)
    
    ! Pour chaque sommet de la facette, on renvoie l'indice de la facette et la position du noeud (1, 2 ou 3) dans le tableau de d�nombrement des facettes dont le noeud est un sommet 
    do r=1,3
        k=1
        do while (Maillage%Tnoeud(lignef(r+1))%Tfacette(k,1)>0.01)
        k=k+1
        end do
        Maillage%Tnoeud(lignef(r+1))%Tfacette(k,:)=[nf,r]
    end do
        
    ! Calcul du centre de gravit� de chaque facette, de son aire et de la distance maximale du centre aux sommets
    Pg=Maillage%Tfacette(j)%Pnoeud
    Maillage%Tfacette(j)%Gfacette=sum(Maillage%Tfacette(j)%Pnoeud,2)/3
    Maillage%Tfacette(j)%Aire=0.5*norm(vect_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1)))
    Maillage%Tfacette(j).Rmax=max(norm(Pg(:,1)-Maillage%Tfacette(j)%Gfacette),&
                                  &norm(Pg(:,2)-Maillage%Tfacette(j)%Gfacette),&
                                  &norm(Pg(:,3)-Maillage%Tfacette(j)%Gfacette))
    
   ! Calcul du gradient constant surfacique de la facette (ne d�pend que de la g�om�trie de la facette)    
    delta=(norm(Pg(:,2)-Pg(:,1))**2)*norm(Pg(:,3)-Pg(:,1))**2-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))**2
    ds1=(norm(Pg(:,3)-Pg(:,1))**2*(Pg(:,2)-Pg(:,1))-(dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1))))*(Pg(:,3)-Pg(:,1)))
    ds2=(-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))*(Pg(:,2)-Pg(:,1))+norm(Pg(:,2)-Pg(:,1))**2*(Pg(:,3)-Pg(:,1)))
    Maillage%Tfacette(j)%ds=reshape([-(ds1+ds2),ds1,ds2]/delta,(/3,3/))    
end do

! Op�rations sur les noeuds
do j=1,Maillage.Nnoeud
    ! D�nombrement du nombre de facettes dont le noeud est un sommet
    m=1
    do while (Maillage%Tnoeud(j)%Tfacette(m,1)>0.001)
    m=m+1
    end do
    Maillage%Tnoeud(j)%Nfacette=m-1
    
    ! R�cup�ration du type de fronti�re et de la normale associ�s au noeud: 
    ! La cr�ation de points doubles aux intersections de parois nous permet d'avoir des points de coordonn�es g�om�triques identiques, 
    ! mais appartenant � des parois diff�rentes --> leurs normales et types de fronti�re sont alors diff�rentes
    ! La normale et le type de fronti�re sont par cons�quent ceux de la paroi associ�e. 
    ! Une paroi est d�finie comme une ensemble de facette ayant une normale et un type de fronti�re identique.
    
    
    !! Pour les maillages � noeuds doubles :
    Maillage%Tnoeud(j)%typeNoeud = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(1,1))%typeFrontiere
    Maillage%Tnoeud(j)%Normale = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(1,1))%Normale    
    do k=1,Maillage%Tnoeud(j)%Nfacette
        if (Maillage%Tnoeud(j)%Normale(1) /= Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(1)&
     & .and.Maillage%Tnoeud(j)%Normale(2) /= Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(2)&
     & .and.Maillage%Tnoeud(j)%Normale(3) /= Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(3)) print*, 'probl�me de Normale pour le noeud', j
        if (Maillage%Tnoeud(j)%typeNoeud /= Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%typeFrontiere) print*, 'probl�me de Type de Frontiere pour le noeud', j
    end do

    Maillage%Tnoeud(j)%NDouble=0
    if (j.gt.1) then
        do k=1,j-1
            if (norm(Maillage%Tnoeud(j)%Pnoeud-Maillage%Tnoeud(k)%Pnoeud).lt.Epsilon) then
                Maillage%Tnoeud(j)%NDouble=Maillage%Tnoeud(j)%NDouble+1
                Maillage%Tnoeud(k)%NDouble=Maillage%Tnoeud(k)%NDouble+1
                Maillage%Tnoeud(j)%Double(Maillage%Tnoeud(j)%NDouble)=k
                Maillage%Tnoeud(k)%Double(Maillage%Tnoeud(k)%NDouble)=j
            end if
        end do
    end if

!    !! Pour les maillages sans noeuds doubles :
!    Maillage%Tnoeud(j)%typeNoeud = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(1,1))%typeFrontiere
!    Maillage%Tnoeud(j)%Normale = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(1,1))%Normale
!    
!    
!    do k=2,Maillage%Tnoeud(j)%Nfacette
!        if (Maillage%Tnoeud(j)%typeNoeud.ne.Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%typeFrontiere) then
!            Maillage%Tnoeud(j)%typeNoeud=10
!            Maillage%Tnoeud(j)%Normale2 = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale
!            exit
!!        elseif (norm(vect_product(Maillage%Tnoeud(j)%Normale,Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale)).gt.Epsilon2) then
!!            Maillage%Tnoeud(j)%typeNoeud=11
!!            Maillage%Tnoeud(j)%Normale2 = Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale
!!            exit
!        end if
!    end do
!    Nor = 0
!    do k=1,Maillage%Tnoeud(j)%Nfacette
!        if (Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%typeFrontiere.ne.0) then
!            Nor = Nor + Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale*Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Aire
!        end if
!    end do
!    if (Maillage%Tnoeud(j)%typeNoeud.eq.1) then
!        Maillage%Tnoeud(j)%Normale = Nor/norm(Nor)  
!        Maillage%Tnoeud(j)%Normale2 = 0 
!    elseif (Maillage%Tnoeud(j)%typeNoeud.eq.10) then
!        Maillage%Tnoeud(j)%Normale2 = Nor/norm(Nor)
!    else  
!        Maillage%Tnoeud(j)%Normale2 = 0 
!    end if
    
!    if (Maillage%Tnoeud(j)%typeNoeud.eq.11) then
!        te=[0,0,0]    
!        do k=2,Maillage%Tnoeud(j)%Nfacette           
!            if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(1))>0.001) then
!                te(1)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(1);
!            end if
!            if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(2))>=0.001) then
!                te(2)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(2);
!            end if
!            if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(3))>0.001) then
!                te(3)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(3);
!            endif
!        end do
!        Maillage%Tnoeud(j)%Normale=te/norm(te)
!    end if 
    
!    Maillage%Tnoeud(j)%typeNoeud=10
!    Maillage%Tnoeud(j)%Normale=[0,0,0]    
!       
!    
!    do k=1,Maillage%Tnoeud(j)%Nfacette
!        Maillage%Tnoeud(j)%typeNoeud=min(Maillage%Tnoeud(j)%typeNoeud,Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%typeFrontiere)
!        
!        if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(1))>0.001) then
!            te(1)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(1);
!        end if
!        if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(2))>=0.001) then
!            te(2)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(2);
!        end if
!        if (abs(Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(3))>0.001) then
!            te(3)=Maillage%Tfacette(Maillage%Tnoeud(j)%Tfacette(k,1))%Normale(3);
!        endif     
!    end do
!    Maillage%Tnoeud(j)%Normale=te/norm(te)
end do
close(10)
print*, '   Extraction du maillage ok'

end subroutine extract_maillage








subroutine Geom(Maillage)  !! n'est plus utile, tout est dans la fonction Extract_Maillage !!
use Constantes
use StructuresDonnees
use FonctionsCommunes
implicit none 

! Param�tres
Type(TMaillage) :: Maillage

! Locales
real(rp) :: delta
real(rp), dimension(3)  ::Nor, ds1, ds2
real(rp), dimension(3,3):: Pg

Maillage%Tfacette(1)%Pnoeud=reshape([Maillage%Tnoeud(1)%Pnoeud,Maillage%Tnoeud(2)%Pnoeud,&
                                    &Maillage%Tnoeud(3)%Pnoeud],(/3,3/))
                                    
Pg=Maillage%Tfacette(1)%Pnoeud

Nor=vect_product(Pg(:,2)-Pg(:,1), Pg(:,3)-Pg(:,1))
Maillage%Tfacette(1)%Normale=Nor/norm(Nor)

Maillage%Tfacette(1)%Gfacette=sum(Pg,2)/3.
Maillage%Tfacette(1)%Aire=0.5*norm(vect_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1)))
Maillage%Tfacette(1).Rmax=max(norm(Pg(:,1)-Maillage%Tfacette(1)%Gfacette),&
                                  &norm(Pg(:,2)-Maillage%Tfacette(1)%Gfacette),&
                                  &norm(Pg(:,3)-Maillage%Tfacette(1)%Gfacette))
   !! Calcul du gradient constant surfacique de la facette (ne d�pend que de la g�om�trie de la facette)
    
delta=(norm(Pg(:,2)-Pg(:,1))**2)*norm(Pg(:,3)-Pg(:,1))**2-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))**2;
ds1=(norm(Pg(:,3)-Pg(:,1))**2*(Pg(:,2)-Pg(:,1))-(dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1))))*(Pg(:,3)-Pg(:,1)));
ds2=(-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))*(Pg(:,2)-Pg(:,1))+norm(Pg(:,2)-Pg(:,1))**2*(Pg(:,3)-Pg(:,1)));
Maillage%Tfacette(1)%ds=reshape([-(ds1+ds2),ds1,ds2]/delta,(/3,3/));


end subroutine Geom
