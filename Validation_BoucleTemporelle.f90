!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                  !
!!                        Validation de la Boucle Temporelle                                        !
!!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    On considère une cuve fermée contenant une surface libre. On lâche une onde stationnaire dans la cuve sous la forme d'une Houle d'Airy. 
!!    On propage la houle à l'intérieur de la cuve et on observe le caractère stationnaire. 
!!
!!
!!
!!program main
!subroutine Validation_BoucleTemporelle(filename)
!
!use Constantes
!use Structuresdonnees
!use FonctionsCommunes
!use GenMaillage
!use Incident_mod
!use BoucleTemp
!use PrePlot
!
!implicit none
!
!character(len=50) :: filename !, fileout, filemaill
!type(TMaillage) :: MaillageInit !, MaillageT, MaillageDT
!type(TEcoulement) :: EcoulementDT !EcoulementInit, EcoulementT, 
!integer :: j !, k, l, ios
!real(rp) :: t0 !, time_begin, tmoy, intvelocity
!real(rp) :: h=0.1
!!real(rp), dimension(3) :: GPhi0, M
!real(rp), dimension(100) :: t !, time_end
!!real(rp), allocatable :: CDT(:,:), CST(:,:), CDDT(:,:), CSDT(:,:)
!logical :: oplot
!!print*, 'konde = ', konde
!t0=-h+pi/(w(1)*2._RP)
!t=(/(t0+j*h,j=1,100)/)
!oplot=.false.
!!   Récupération du maillage depuis le fichier maillage
!call Extract_Maillage(filename, MaillageInit)
!
!call BoucleTemporelle(filename, 100, t, MaillageInit, EcoulementDT, oplot)
!
!!call NewMaillage(MaillageT, PointMax)
!!call Remaillage(0._RP, MaillageInit, MaillageT)
!!filemaill = 'RemaillageT_'//filename
!!call PlotMaill(filemaill, MaillageT)
!!
!!call NewEcoulement(EcoulementT, MaillageT.Nnoeud)
!
!!call CCL0( EcoulementT, MaillageT)
!!
!!do j=1,MaillageT.Nnoeud
!!!CL pour onde stationnaire
!!    M = MaillageT%Tnoeud(j)%Pnoeud
!!    EcoulementT.Phi(j)%incident=Phi0(M,t(1))
!!    call CGPhi0(M,t(1),GPhi0)
!!    EcoulementT.DPhiDn(j)%incident=-dot_product(GPhi0,MaillageT%Tnoeud(j)%Normale)
!!    
!!!CL pour batteur
!!!    EcoulementT.Phi(j)%incident=0
!!!    if (abs(MaillageT%Tnoeud(j)%Pnoeud(1)-0).lt.Epsilon2) then
!!!        EcoulementT.DPhiDn(j)%incident=intvelocity(t(1))
!!!    else
!!!        EcoulementT.DPhiDn(j)%incident=0
!!!    end if
!!    
!!    if (MaillageT%Tnoeud(j)%TypeNoeud.eq.0) then
!!        EcoulementT.Phi(j)%perturbation=EcoulementT.Phi(j)%incident
!!    elseif (MaillageT%Tnoeud(j)%TypeNoeud.eq.1) then
!!        EcoulementT.DPhiDn(j)%perturbation=0 !EcoulementT.DPhiDn(j)%incident
!!    elseif (MaillageT%Tnoeud(j)%TypeNoeud.eq.10) then
!!        EcoulementT.Phi(j)%perturbation=EcoulementT.Phi(j)%incident
!!        EcoulementT.DPhiDn2(j)%incident = 0 !-dot_product(GPhi0,Maillage%Tnoeud(j).Normale2)
!!        EcoulementT.DPhiDn2(j)%perturbation=EcoulementT.DPhiDn2(j)%incident
!!   end if
!!    EcoulementT.Eta(j)%perturbation = EcoulementT.Eta(j)%incident
!!end do
!
!!fileout='output_'//filename
!!open(unit=33,file=fileout, iostat=ios)
!!if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
!!write(33,fmt='(50a)') 'Title = "Maillage de la cuve"'
!!write(33,fmt='(50a)') 'VARIABLES = "X","Y","Eta"'
!!
!!call cpu_time(time_begin)
!!allocate(CST(MaillageT.Nnoeud,MaillageT.Nnoeud), CDT(MaillageT.Nnoeud,MaillageT.Nnoeud))
!!if (1>0) then
!!    CDT=-1
!!    do j=1,size(t)
!!        print*, '   '
!!        print*, 'Pas de temps : t = ', t(j), ' , indice : ', j
!!        oplot=.false.
!!        call BoucleTemporelle(filename, t(j), h, MaillageT, MaillageDT, MaillageInit, EcoulementT, EcoulementDT, CDT, CST, CDDT, CSDT, oplot)
!!        print*, '   --> Pas de temps : t = ', t(j), ' acheve'    
!!    !    call outplot(filename, MaillageDT, EcoulementDT)
!!        call PlotTemp(fileout, t(j), MaillageDT, EcoulementDT)    
!!
!!        call DelEcoulement(EcoulementT)
!!        call NewEcoulement(EcoulementT, MaillageDT.Nnoeud)
!!        call CopyEcoulement(EcoulementT, EcoulementDT, MaillageDT.Nnoeud)
!!        call DelEcoulement(EcoulementDT)
!!        
!!        deallocate(CDT, CST)
!!        allocate(CST(MaillageDT.Nnoeud,MaillageDT.Nnoeud), CDT(MaillageDT.Nnoeud,MaillageDT.Nnoeud))
!!        CDT=CDDT
!!        CST=CSDT
!!        deallocate(CDDT, CSDT)
!!        
!!        call DelMaillage(MaillageT)
!!        call NewMaillage(MaillageT, PointMax)
!!        call CopyMaillage(MaillageT, MaillageDT)
!!        call DelMaillage(MaillageDT)
!!        
!!        call cpu_time(time_end(j))    
!!        if (j==1) then
!!        tmoy = time_end(j) - time_begin
!!        print*, 'Le temps de l''operation a ete de ', time_end(j) - time_begin, ' secondes'
!!        else
!!        tmoy = (tmoy*(j-1) + time_end(j) - time_end(j-1))/j
!!        print*, 'Le temps de l''operation a ete de ', time_end(j) - time_end(j-1), ' secondes'
!!        end if
!!        print*, 'fin prevu dans ', tmoy*(size(t)-j)/60._RP, ' minutes'
!!    end do
!!else
!!    CDT=-1
!!    do j=1,10 !size(t)
!!        print*, '   '
!!        print*, 'Pas de temps : t = ', t(j), ' , indice : ', j
!!        oplot=.false.
!!        call BoucleTempLin(t(j), h, MaillageT, EcoulementT, EcoulementDT, filename, CDT, CST, oplot)
!!        print*, '   --> Pas de temps : t = ', t(j), ' acheve'
!!        call PlotTemp(fileout, t(j), MaillageT, EcoulementDT)    
!!
!!        call DelEcoulement(EcoulementT)
!!        call NewEcoulement(EcoulementT, MaillageT.Nnoeud)
!!        call CopyEcoulement(EcoulementT, EcoulementDT, MaillageT.Nnoeud)
!!        call DelEcoulement(EcoulementDT)
!!
!!        call cpu_time(time_end(j))    
!!        if (j==1) then
!!        tmoy = time_end(j) - time_begin
!!        print*, 'Le temps de l''operation a ete de ', time_end(j) - time_begin, ' secondes'
!!        else
!!        tmoy = (tmoy*(j-1) + time_end(j) - time_end(j-1))/j
!!        print*, 'Le temps de l''operation a ete de ', time_end(j) - time_end(j-1), ' secondes'
!!        end if
!!        print*, 'fin prevu dans ', tmoy*(size(t)-j)/60._RP, ' minutes'
!!    end do
!!end if
!!close(33)
!
!end subroutine Validation_BoucleTemporelle
!
!!function intvelocity(t)
!!use Constantes
!!implicit none
!!real(rp) :: t
!!real(rp) :: intvelocity
!!real(rp) :: Amplitude, Omega
!!
!!Amplitude = 0
!!Omega = 1/pi
!!intvelocity=Amplitude*sin(2*pi*omega*t)
!!
!!end function
!
