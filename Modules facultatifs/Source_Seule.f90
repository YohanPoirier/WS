!module SourceSeule
!    use Constantes
!    use Structuresdonnees
!    use FonctionsCommunes
!    implicit none
!    
!    contains
!    
!    subroutine PhiSS(Maillage,S1,S2,Ecoulement)
!    ! donne la Solution Analytique du potentiel et de la vitesse normale engendrés par 2 sources singulières
!    ! disposées aux positions S1 et S2
!    use Constantes
!    use Structuresdonnees
!    use FonctionsCommunes
!    !Paramètres
!    type(TMaillage) :: Maillage
!    real(rp), dimension(3) :: S1,S2 ! Positions des deux sources
!    !Résultats
!    type(TEcoulement) :: Ecoulement
!    !Locales
!    real :: M1P,M2P
!    integer :: j    
!    
!    call NewEcoulement(Ecoulement,Maillage.Nnoeud)
!    
!    do j=1,Maillage.Nnoeud
!        M1P=1/norm(Maillage.Tnoeud(j).Pnoeud - S1)
!        M2P=0 !1/norm(Maillage.Tnoeud(j).Pnoeud - S2)
!        Ecoulement.Phi(j).Incident = M1P + M2P
!        
!        Ecoulement.DPhiDn(j).Incident=(M1P**3)*dot_product(S1-Maillage.Tnoeud(j).Pnoeud,Maillage.Tnoeud(j).Normale)&
!                                   & + (M2P**3)*dot_product(S2-Maillage.Tnoeud(j).Pnoeud,Maillage.Tnoeud(j).Normale)
!    end do  
!    end subroutine
!    
!    
!    
!end module SourceSeule