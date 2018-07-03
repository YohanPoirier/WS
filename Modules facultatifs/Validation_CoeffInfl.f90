!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                            Validation du calcul des Coefficients d'Influence                     !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   On choisit une facette triangulaire de référence: isocèle rectangle de côté 1. On choisit un point
! influencé M.

!program main
subroutine validation_CoeffInfl(Maillage)

use Constantes
use GenMaillage
use Structuresdonnees
use FonctionsCommunes

! Paramètres!
character(len=50)   :: fileInfl
type(TMaillage) :: Maillage
!real(rp), dimension(3)  :: pointM
!integer :: N,j
real(rp),allocatable   ::  CD(:,:),CS(:,:) 

call extract_maillage(Maillage)

allocate(CD(Maillage%Nnoeud,Maillage%Nnoeud),CS(Maillage%Nnoeud,Maillage%Nnoeud))

!call Infl(Maillage,CD,CS)

!call CoeffInfl(Maillage,CD,CS)

fileInfl='CoeffInfl_'//filename
open(4,file=fileInfl)
write(4,'(50a)') 'Title= "Coefficients d`Influence"'
write(4,'(a,i,a,i)') 'Nombre de Noeuds =', Maillage%Nnoeud, ', Nombre de Facettes=', Maillage%Nfacette
 
    do j=1,Maillage%Nnoeud
        do k=1,Maillage%Nnoeud
            write(4,*), CD(j,k), CS(j,k)
        end do
    end do
close(4)

end subroutine validation_CoeffInfl
!end program main

subroutine Infl(Maillage,CD,CS) !! N'est plus utile --> Voir la fonction CoeffInfl !!
use Constantes
use Structuresdonnees
use FonctionsCommunes
implicit none

!Paramètres
type(TMaillage) :: Maillage

!Résultats
!integer, dimension(Maillage%Nnoeud*Maillage%Nfacette,3) :: test
real(rp), dimension(Maillage%Nnoeud*Maillage%Nfacette)  :: Ssigma, Smu
real(rp), dimension(Maillage%Nnoeud*Maillage%Nfacette,3)    :: Isigma, Imu, Ar
!real(rp), dimension(3,3,Maillage%Nnoeud*Maillage%Nfacette)  :: Jsigma
real(rp), dimension(3)      :: Css, Cdd
!real(rp),allocatable   ::  CD(:,:),CS(:,:),CDl(:,:),CSl(:,:) 
real(rp), dimension(Maillage%Nnoeud*Maillage%Nfacette,Maillage%Nnoeud)   :: CDl,CSl
real(rp), dimension(Maillage%Nnoeud,Maillage%Nnoeud)   :: CD,CS

!Locales
    integer :: j,k,l,t
    real(rp) ::  Rho, Cr, Z
    real(rp), dimension(3) ::   MG
    real(rp), dimension(3,3) :: ds
    
!allocate(CD(Maillage%Nnoeud,Maillage%Nnoeud),CS(Maillage%Nnoeud,Maillage%Nnoeud))
!allocate(CDl(Maillage%Nnoeud*Maillage%Nfacette,Maillage%Nnoeud),CSl(Maillage%Nnoeud*Maillage%Nfacette,Maillage%Nnoeud))
CDl(:,:)=0
CSl(:,:)=0
CD(:,:)=0
CS(:,:)=0
l=1
do j=1,Maillage%Nnoeud
    print*, 100*j/Maillage%Nnoeud, '.'
    do k=1,Maillage%Nfacette    
        !call Iints(Maillage%Tnoeud(j), Maillage%Tfacette(k), Isigma(l,:), Imu(l,:))
        MG=Maillage%Tnoeud(j)%Pnoeud - Maillage%Tfacette(k)%Gfacette

! Détermination de la distance du point de contrôle à la facette --> Critère de calcul
        Rho=norm2(MG)
        Cr=Rho/Maillage%Tfacette(k)%Rmax
        if (Cr<=8) then
            !call Sints(Maillage%Tnoeud(j), Maillage%Tfacette(k),Ssigma(l),Smu(l))
        else
            Ssigma(l)=Maillage%Tfacette(k)%Aire/Rho
            Z= dot_product(MG, Maillage%Tfacette(k)%Normale)
            Smu(l)=Z*Maillage%Tfacette(k)%Aire/(Rho**3)
        end if
        ds=transpose(Maillage%Tfacette(k)%ds)        
        Ar(l,:)=[1,1,1]/3. + matmul(ds,MG)
        Css=Ssigma(l)*Ar(l,:) + matmul(ds,Isigma(l,:))
        Cdd=Smu(l)*Ar(l,:) + matmul(ds,Imu(l,:))
        
        CS(j,Maillage%Tfacette(k)%Tnoeud(1))=CS(j,Maillage%Tfacette(k)%Tnoeud(1))+Css(1)
        CS(j,Maillage%Tfacette(k)%Tnoeud(2))=CS(j,Maillage%Tfacette(k)%Tnoeud(2))+Css(2)
        CS(j,Maillage%Tfacette(k)%Tnoeud(3))=CS(j,Maillage%Tfacette(k)%Tnoeud(3))+Css(3)
        CD(j,Maillage%Tfacette(k)%Tnoeud(1))=CD(j,Maillage%Tfacette(k)%Tnoeud(1))+Cdd(1)
        CD(j,Maillage%Tfacette(k)%Tnoeud(2))=CD(j,Maillage%Tfacette(k)%Tnoeud(2))+Cdd(2)
        CD(j,Maillage%Tfacette(k)%Tnoeud(3))=CD(j,Maillage%Tfacette(k)%Tnoeud(3))+Cdd(3)
        
        
        CSl(l,Maillage%Tfacette(k)%Tnoeud(1))=Css(1)
        CSl(l,Maillage%Tfacette(k)%Tnoeud(2))=Css(2)
        CSl(l,Maillage%Tfacette(k)%Tnoeud(3))=Css(3)
        CDl(l,Maillage%Tfacette(k)%Tnoeud(1))=Cdd(1)
        CDl(l,Maillage%Tfacette(k)%Tnoeud(2))=Cdd(2)
        CDl(l,Maillage%Tfacette(k)%Tnoeud(3))=Cdd(3)        
        l=l+1
    end do
end do

open(4,file='CoeffInfl.dat')
write(4,'(50a)') 'Title= "Coefficients d`Influence"'
write(4,'(a,i,a,i)') 'Nombre de Noeuds =', Maillage%Nnoeud, ', Nombre de Facettes=', Maillage%Nfacette
 
    do j=1,Maillage%Nnoeud
        do k=1,Maillage%Nnoeud
            write(4,*), CD(j,k), CS(j,k)
        end do
    end do

do j=1,Maillage%Nnoeud
    write(4,fmt='(a,26F)') 'CD = ', CD(j,:)
end do
do j=1,Maillage%Nnoeud
    write(4,fmt='(a,26F)') 'CS = ', CS(j,:)
end do



l=1
do j=1,Maillage%Nnoeud
    do k=1,Maillage%Nfacette
        write(4,fmt='(a,26F,2i)') 'CDl = ', CDl(l,:),j,k
        l=l+1
    enddo
end do

l=1
do j=1,Maillage%Nnoeud
    do k=1,Maillage%Nfacette
    write(4,fmt='(a,26F,2i)') 'CSl = ', CSl(l,:),j,k
        l=l+1
    enddo
end do
!close(4)





!open(44,file='valid_coeffInfl.dat')
    print*, 'test'
    do j=1,Maillage%Nnoeud*Maillage%Nfacette
        write(4,fmt='(a,3E)') 'Isigma = ', Isigma(j,:)
    end do
    do j=1,Maillage%Nnoeud*Maillage%Nfacette   
        write(4,fmt='(a,3E)') 'Imu = ', Imu(j,:)
    end do
    do j=1,Maillage%Nnoeud*Maillage%Nfacette  
        write(4,fmt='(a,E)') 'Ssigma = ', Ssigma(j)
    end do
    do j=1,Maillage%Nnoeud*Maillage%Nfacette   
        write(4,fmt='(a,E)') 'Smu = ', Smu(j)
    end do
    do j=1,Maillage%Nnoeud*Maillage%Nfacette   
        write(4,fmt='(a,3E)') 'Ar = ', Ar(j,:)
    end do
!    do j=1,Maillage%Nnoeud*Maillage%Nfacette
!        do k=1,3   
!            write(4,fmt='(a,3E)') 'Jsigma = ', Jsigma(:,k,j)
!        end do
!    end do
!    l=1
!    do j=1,Maillage%Nnoeud
!        do k=1,Maillage%Nfacette   
!            write(4,fmt='(a,5i)') 'test = ', test(l,:), j, k
!            l=l+1
!        end do
!    end do
close(4)



end subroutine Infl


        subroutine CoeffInfl0(Maillage,CD,CS)   !!! Pas fini?!
 
use Constantes
use Structuresdonnees
use FonctionsCommunes
implicit none

!Paramètres
type(TMaillage) :: Maillage

!Résultats
real(rp)  :: Ssigma, Smu
real(rp), dimension(3)    :: Isigma, Imu, Ar,M
real(rp), dimension(3,4)    :: Pm
real(rp), dimension(3)      :: Css, Cdd
real(rp), dimension(Maillage%Nfacette,Maillage%Nfacette)   :: CD,CS

!Locales
    integer :: j,k
    real(rp) ::  Rho, Cr, Z
    real(rp), dimension(3) ::   MG
    real(rp), dimension(3,3) :: ds

CD(:,:)=0
CS(:,:)=0

do j=1,Maillage%Nnoeud
    !print*, 100*j/Maillage%Nnoeud, '.'
    do k=1,Maillage%Nfacette   
    
        if (j==k) then
            CD(j,k)=0;
            M=Maillage%Tfacette(j)%Gfacette
            Pm(:,1:3)=Maillage%Tfacette(k)%Pnoeud
            Pm(:,4)=Maillage%Tfacette(k)%Pnoeud(1,:)
        end if
    
     
        !call Iints(Maillage%Tnoeud(j), Maillage%Tfacette(k), Isigma, Imu)
        MG=Maillage%Tnoeud(j)%Pnoeud - Maillage%Tfacette(k)%Gfacette

! Détermination de la distance du point de contrôle à la facette --> Critère de calcul
        Rho=norm2(MG)
        Cr=Rho/Maillage%Tfacette(k)%Rmax
        if (Cr<=8) then
            !call Sints(Maillage%Tnoeud(j), Maillage%Tfacette(k),Ssigma,Smu)
        else
            Ssigma=Maillage%Tfacette(k)%Aire/Rho
            Z= dot_product(MG, Maillage%Tfacette(k)%Normale)
            Smu=Z*Maillage%Tfacette(k)%Aire/(Rho**3)
        end if
        ds=transpose(Maillage%Tfacette(k)%ds)        
        Ar=[1,1,1]/3. + matmul(ds,MG)
        Css=Ssigma*Ar + matmul(ds,Isigma)
        Cdd=Smu*Ar + matmul(ds,Imu)
        
        CS(j,Maillage%Tfacette(k)%Tnoeud(1))=CS(j,Maillage%Tfacette(k)%Tnoeud(1))+Css(1)
        CS(j,Maillage%Tfacette(k)%Tnoeud(2))=CS(j,Maillage%Tfacette(k)%Tnoeud(2))+Css(2)
        CS(j,Maillage%Tfacette(k)%Tnoeud(3))=CS(j,Maillage%Tfacette(k)%Tnoeud(3))+Css(3)
        CD(j,Maillage%Tfacette(k)%Tnoeud(1))=CD(j,Maillage%Tfacette(k)%Tnoeud(1))+Cdd(1)
        CD(j,Maillage%Tfacette(k)%Tnoeud(2))=CD(j,Maillage%Tfacette(k)%Tnoeud(2))+Cdd(2)
        CD(j,Maillage%Tfacette(k)%Tnoeud(3))=CD(j,Maillage%Tfacette(k)%Tnoeud(3))+Cdd(3)
    end do
end do

end subroutine CoeffInfl0