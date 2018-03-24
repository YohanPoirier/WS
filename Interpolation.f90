!!-------------------------------------------------------------
!!
!!	Interpolation de la perturbation sur la surface au noeud 
!!	du nouveau maillage
!!
!!--------------------------------------------------------------
!!
!	subroutine Interpolation( Maillage0, Maillage, Ecoulement0, Ecoulement,t  )
!!
!	use Constantes
!	use Structuresdonnees
!	use FonctionsCommunes
!!
!	implicit none
!	
!	!Paramètres
!	type(TMaillage) :: Maillage0, Maillage
!	type(TEcoulement) :: Ecoulement0, Ecoulement 	
!	real(rp) :: t
!	
!	!Variables Locales
!	integer:: j, k
!	real(rp) :: d
!	real(rp), dimension(3) :: M0, M1, Phifacette
!	real(rp), dimension(3,3) :: ds
!	
!	d=100
!	
!	do j=1,Maillage.Nnoeud
!	    do k=1,Maillage0.Nnoeud
!	        M0 = Maillage0.Tnoeud(k).Pnoeud
!	        M1 = Maillage.Tnoeud(j).Pnoeud
! 	        d=norm(M1-M0)
!	        if (d<minD) then
!	            minD=d
!	            minI=k	        
!	        end if
!	    end do
!	    
!	    M0 = Maillage0.Tnoeud(minI).Pnoeud
!	    M1 = Maillage.Tnoeud(j).Pnoeud
!	    do k=1,Maillage0.Tnoeud(minI).Nfacette
!	        if (Ecoulement0.GPhi(minI)==0) then
!	            ds=transpose(Maillage0.Tfacette(Maillage.Tnoeud(minI).Tfacette(k,1)).ds)
!	            Phifacette=[Ecoulement0.Phi(Maillage0.Tfacette(Maillage.Tnoeud(minI).Tnoeud(1)).perturbation,&&
!	            Ecoulement0.Phi(Maillage0.Tfacette(Maillage.Tnoeud(minI).Tnoeud(2)).perturbation,&&
!	            Ecoulement0.Phi(Maillage0.Tfacette(Maillage.Tnoeud(minI).Tnoeud(3)).perturbation,]
!	            Ecoulement0.GPhi(minI)=matmul(ds,Phifacette)
!	        end if
!	        Phi(k)=Ecoulement0.Phi(minI).perturbation + dot_product(Ecoulement0.GPhi(minI),M1-M0)
!	        Eta(k)=Ecoulement0.Eta(minI).perturbation + dot_product(Ecoulement0.GEta(minI),M1-M0)
!	    end do
!	    Ecoulement.Phi(j).perturbation=sum(Phi)/Maillage0.Tnoeud(minI).Nfacette
!	    Ecoulement.Eta(j).perturbation=sum(Eta)/Maillage0.Tnoeud(minI).Nfacette
!	    if (abs(Ecoulement.Eta(j).perturbation-eta0(M1,t)) .gt. 0.1) then
!	        print*, 'composante de perturbation de Eta a t = ', t, ' supérieure de 10. de Eta incident'
!	    end if
!	end do
!	
!	
!	end subroutine Interpolation