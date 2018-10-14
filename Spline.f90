module Spline
use Constantes
use FonctionsCommunes
use Incident_mod
implicit none

contains

subroutine Gradient(Mesh, Ecoulement)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh         ! Mesh
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement   ! Flow parameters

    
    
    ! This subroutine computes the gradients (potential and wave elevation) on the free surface, on the floaters and at the intersections (because of the double nodes).

    call GradientFS(Mesh, Ecoulement)
    call GradientBody(Mesh, Ecoulement)
    call GradientIntersection(Mesh, Ecoulement)

end subroutine Gradient



subroutine GradientFS_modif(Mesh, Ecoulement)
	
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement                   ! Flow parameters.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh                         ! Mesh.
    
    integer                             :: j, k, r                      ! Loop parameters.
    integer                             :: Na, Nvoisin                  ! Spline ordre and number of neighbours.
    real(rp)                            :: Ponderation
    real(rp), dimension(2)              :: DF1, DF2
    real(rp), dimension(2,2)            :: TenseurMetrique
    real(rp), dimension(3)              :: DFDu, DFDv, Phifacette
    real(rp), allocatable               :: Pvoisin(:,:), A(:,:), B(:,:)
    integer                             :: info
    character(len=1)                    :: trans
    integer                             :: lda, nrhs
    integer, allocatable                :: ipiv(:) 
    

    
    ! This subroutine computes the gradients (potential and wave elevation) of the velocity potential on the free surface.
    
    do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        
        if (.false.) then
            call gradspline( j, Ecoulement, Mesh, .true.)
            Ecoulement%GPhi(:,j)%perturbation = Ecoulement%GPhi(:,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
        elseif(is_BS_GS) then
                

            Na = Nordre(Mesh%Tnoeud(j)%Ordre)
            NVoisin = Mesh%Tnoeud(j)%NVoisin(2)-1
            allocate(Pvoisin(3,NVoisin+1), B(NVoisin+Na,2), A(NVoisin+Na,NVoisin+Na))
            
            do k = 1,NVoisin+1
                Pvoisin(1:3,k) = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%Pnoeud
                if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) Pvoisin(2,k) = - Pvoisin(2,k)
                B(k,1) = Ecoulement%Phi(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%perturbation
                B(k,2) = Ecoulement%Eta(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%perturbation
            end do
            B(Nvoisin+2:Nvoisin+Na,1:2) = 0._RP

            ! On construit A.
            call SplineMatrix(Nvoisin, Mesh%Tnoeud(j)%Ordre, Pvoisin, A)
            

            ! Résolution du Système Linéaire.
            lda = Nvoisin+Na
            allocate( ipiv(lda) )
            ipiv=0
            trans='n'
            nrhs=2 ! Nombre de second membre.
            call dgetrf(lda,lda,A,lda,ipiv,info)
            if (info.ne.0) then
                print*, 'Warning GradientFS 69 : info =', info, lda, j, Mesh%Tnoeud(j)%NVoisin(2), Na
                print*, Mesh%Tnoeud(j)%TVoisin(1:Mesh%Tnoeud(j)%NVoisin(2),1)
                pause
            end if
            call dgetrs(trans,lda,nrhs,A,lda,ipiv,B,lda,info)
            deallocate(ipiv)
            
 
            ! Calcul des Dérivées Premières de la géométrie locale.
            call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,1), Pvoisin(1:3,1), Pvoisin, DF1)
            call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,2), Pvoisin(1:3,1), Pvoisin, DF2)
            
            ! Retour à la base globale.
            DfDu = [1._RP,0._RP,Mesh%Tnoeud(j)%DLocal(1)]
            DfDv = [0._RP,1._RP,Mesh%Tnoeud(j)%DLocal(2)]
            TenseurMetrique = reshape([1+Mesh%Tnoeud(j)%DLocal(2)**2, -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2),&
                         -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2), 1+Mesh%Tnoeud(j)%DLocal(1)**2]/Mesh%Tnoeud(j)%DLocal(5),(/2,2/))
            
            ! Calcul du gradient surfacique de la déformée de surface libre.
            Ecoulement%GEta(1:2,j)%perturbation = DF2
            Ecoulement%GEta(3,j)%perturbation = 0._RP
            
            ! Calcul du gradient surfacique du potentiel.
            Ecoulement%GPhi(1:3,j)%perturbation = DF1(1)*(TenseurMetrique(1,1)*DfDu + TenseurMetrique(1,2)*DfDv) + DF1(2)*(TenseurMetrique(2,1)*DfDu + TenseurMetrique(2,2)*DfDv)
            
            ! Ajout de la composante normale.
            Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
            deallocate(Pvoisin, A, B)
    
        else
            
            ! Linear Panel approximation.
            Ecoulement%GEta(1:3,j)%perturbation = 0._RP
            Ecoulement%GPhi(1:3,j)%perturbation = 0._RP
            Ponderation = 0._RP
            
            do k = 1,Mesh%Tnoeud(j)%Nfacette
                Ponderation = Ponderation + Mesh%Tnoeud(j)%Angle(k)
                do r = 1,3
                    Phifacette(r)=Ecoulement%Eta(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation
                end do
                Ecoulement%GEta(1:3,j)%perturbation = Ecoulement%GEta(1:3,j)%perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k)
                do r = 1,3
                    Phifacette(r)=Ecoulement%Phi(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation
                end do
                Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k)
            end do
            
            Ecoulement%GEta(1:3,j)%perturbation = Ecoulement%GEta(1:3,j)%perturbation/Ponderation
            Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation/Ponderation
            
            if (Symmetry.and.abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon) then
                !Ecoulement%GEta(1:3,j)%perturbation = 2._RP*Ecoulement%GEta(1:3,j)%perturbation
                Ecoulement%GEta(2,j)%perturbation = 0._RP
                !Ecoulement%GPhi(1:3,j)%perturbation = 2._RP*Ecoulement%GPhi(1:3,j)%perturbation
                Ecoulement%GPhi(2,j)%perturbation = 0._RP
            end if
            
            Ecoulement%GPhi(:,j)%perturbation = Ecoulement%GPhi(:,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
        end if
        
    end do
    
end subroutine GradientFS_modif





subroutine GradientFS(Mesh, Ecoulement)
	
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement                   ! Flow parameters.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh                         ! Mesh.
    
    integer                             :: j, k, r                      ! Loop parameters.
    integer                             :: Na, Nvoisin                  ! Spline ordre and number of neighbours.
    real(rp)                            :: Ponderation
    real(rp), dimension(2)              :: DF1, DF2
    real(rp), dimension(2,2)            :: TenseurMetrique
    real(rp), dimension(3)              :: DFDu, DFDv, Phifacette
    real(rp), allocatable               :: Pvoisin(:,:), A(:,:), B(:,:)
    integer                             :: info
    character(len=1)                    :: trans
    integer                             :: lda, nrhs
    integer, allocatable                :: ipiv(:) 
    

    
    ! This subroutine computes the gradients (potential and wave elevation) of the velocity potential on the free surface.
    
    

    
    do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        
        if (.false.) then
            call gradspline( j, Ecoulement, Mesh, .true.)
            Ecoulement%GPhi(:,j)%perturbation = Ecoulement%GPhi(:,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
        elseif(is_BS_GS) then
                

            Na = Nordre(Mesh%Tnoeud(j)%Ordre)
            NVoisin = Mesh%Tnoeud(j)%NVoisin(2)-1
            allocate(Pvoisin(3,NVoisin+1), B(NVoisin+Na,2), A(NVoisin+Na,NVoisin+Na))
            
            do k = 1,NVoisin+1
                Pvoisin(1:3,k) = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%Pnoeud
                if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) Pvoisin(2,k) = - Pvoisin(2,k)
                B(k,1) = Ecoulement%Phi(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%perturbation
                B(k,2) = Ecoulement%Eta(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%perturbation
            end do
            B(Nvoisin+2:Nvoisin+Na,1:2) = 0._RP

            ! On construit A.
            call SplineMatrix(Nvoisin, Mesh%Tnoeud(j)%Ordre, Pvoisin, A)
            

            ! Résolution du Système Linéaire.
            lda = Nvoisin+Na
            allocate( ipiv(lda) )
            ipiv=0
            trans='n'
            nrhs=2 ! Nombre de second membre.
            call dgetrf(lda,lda,A,lda,ipiv,info)
            if (info.ne.0) then
                print*, 'Warning GradientFS 69 : info =', info, lda, j, Mesh%Tnoeud(j)%NVoisin(2), Na
                print*, Mesh%Tnoeud(j)%TVoisin(1:Mesh%Tnoeud(j)%NVoisin(2),1)
                pause
            end if
            call dgetrs(trans,lda,nrhs,A,lda,ipiv,B,lda,info)
            deallocate(ipiv)
            
 
            ! Calcul des Dérivées Premières de la géométrie locale.
            call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,1), Pvoisin(1:3,1), Pvoisin, DF1)
            call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,2), Pvoisin(1:3,1), Pvoisin, DF2)
            
            ! Retour à la base globale.
            DfDu = [1._RP,0._RP,Mesh%Tnoeud(j)%DLocal(1)]
            DfDv = [0._RP,1._RP,Mesh%Tnoeud(j)%DLocal(2)]
            TenseurMetrique = reshape([1+Mesh%Tnoeud(j)%DLocal(2)**2, -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2),&
                         -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2), 1+Mesh%Tnoeud(j)%DLocal(1)**2]/Mesh%Tnoeud(j)%DLocal(5),(/2,2/))
            
            ! Calcul du gradient surfacique de la déformée de surface libre.
            Ecoulement%GEta(1:2,j)%perturbation = DF2
            Ecoulement%GEta(3,j)%perturbation = 0._RP
            
            ! Calcul du gradient surfacique du potentiel.
            Ecoulement%GPhi(1:3,j)%perturbation = DF1(1)*(TenseurMetrique(1,1)*DfDu + TenseurMetrique(1,2)*DfDv) + DF1(2)*(TenseurMetrique(2,1)*DfDu + TenseurMetrique(2,2)*DfDv)
            
            ! Ajout de la composante normale.
            Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
            deallocate(Pvoisin, A, B)
    
        else
            ! Linear Panel approximation.
            Ecoulement%GEta(1:3,j)%perturbation = 0._RP
            Ecoulement%GPhi(1:3,j)%perturbation = 0._RP
            Ponderation = 0._RP
            
            do k = 1,Mesh%Tnoeud(j)%Nfacette
                Ponderation = Ponderation + Mesh%Tnoeud(j)%Angle(k)
                do r = 1,3
                    Phifacette(r)=Ecoulement%Eta(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation
                end do
                Ecoulement%GEta(1:3,j)%perturbation = Ecoulement%GEta(1:3,j)%perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k)
                do r = 1,3
                    Phifacette(r)=Ecoulement%Phi(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation
                end do
                Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k)
            end do
            
            Ecoulement%GEta(1:3,j)%perturbation = Ecoulement%GEta(1:3,j)%perturbation/Ponderation
            Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation/Ponderation
            
            if (Symmetry.and.abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon) then
                !Ecoulement%GEta(1:3,j)%perturbation = 2._RP*Ecoulement%GEta(1:3,j)%perturbation
                Ecoulement%GEta(2,j)%perturbation = 0._RP
                !Ecoulement%GPhi(1:3,j)%perturbation = 2._RP*Ecoulement%GPhi(1:3,j)%perturbation
                Ecoulement%GPhi(2,j)%perturbation = 0._RP
            end if
            
            Ecoulement%GPhi(:,j)%perturbation = Ecoulement%GPhi(:,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
            
        end if
        
    end do
    
end subroutine GradientFS

subroutine GradientBody(Mesh, Ecoulement)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(in)          :: Mesh                                         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement                                   ! Flow parameters.
        
    integer                             :: j, k ,r, s, nc                               ! Loop parameters.
    real(rp), dimension(3)              :: Phifacette
    integer                             :: Na, Nvoisin, i1, i2
    real(rp)                            :: Ponderation
    real(rp), dimension(2)              :: DF1,DF2
    real(rp), dimension(3)              :: GPhiLocal
    real(rp), dimension(2,2)            :: TenseurMetrique
    real(rp), dimension(3)              :: Pk, M, DFDu, DFDv
    real(rp), dimension(3,3)            :: Puvw, Pt, GGPhi
    real(rp), allocatable               :: Pvoisin(:,:), A(:,:), B(:,:), B1(:), Beta(:)
    integer, dimension(:), allocatable  :: ipiv
    character(len=1)                    :: trans
    integer                             :: nrhs,info
    
    ! This subroutine computes the gradient of the velocity potential on the surface of the bodies.
    
    do nc = Int_Body,Mesh%NBody
        
        if(Mesh%Body(nc)%Active)then
        
            i1 = Mesh%Body(nc)%IndBody(1)
            i2 = Mesh%Body(nc)%IndBody(3)
            if (is_BS_GS) then
            
                ! B-Splines approximation.
                do j=i1,i2
                    
                    ! Calcul du gradient du potentiel.
                    Na = Nordre(Mesh%Tnoeud(j)%Ordre)
                    NVoisin = Mesh%Tnoeud(j)%Nvoisin(2)-1
                    allocate(Pvoisin(3,NVoisin+1), B1(NVoisin+Na), A(NVoisin+Na,NVoisin+Na),Beta(NVoisin+Na))
                    Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                    Pt(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
                    
                    ! Changement de variables.
                    do k=1,NVoisin+1
                        Pk = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%Tvoisin(k,1)))%Pnoeud
                        B1(k) = Ecoulement%Phi(abs(Mesh%Tnoeud(j)%Tvoisin(k,1)))%Perturbation
                        if (symmetry .and. Mesh%Tnoeud(j)%Tvoisin(k,1).lt.0) Pk(2)=-Pk(2)
                        M = Mesh%Tnoeud(j)%Pnoeud
                        Pvoisin(1:3,k) = matmul(Pt,Pk-M)
                    end do
                
                    ! On construit B.
                    B1(Nvoisin+2:Nvoisin+Na) = 0._RP
                
                    ! On construit A.
                    call SplineMatrix(Nvoisin, Mesh%Tnoeud(j)%Ordre, Pvoisin, A)
                    call LU(A, B1, Beta, Nvoisin+Na)
                
                    ! Calcul des Dérivées Premières.
                    call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, Beta(1:Nvoisin+Na), Pvoisin(1:3,1), Pvoisin, DF1)
                
                    ! Dérivées locales surfaciques du potentiel dans la base locale.
                    Ecoulement%GPhi2(1:2,j)%perturbation = DF1
                    DfDu=[1._RP,0._RP,Mesh%Tnoeud(j)%DLocal(1)]
                    DfDv=[0._RP,1._RP,Mesh%Tnoeud(j)%DLocal(2)]
                    TenseurMetrique=reshape([1+Mesh%Tnoeud(j)%DLocal(2)**2, -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2),&
                                 -Mesh%Tnoeud(j)%DLocal(1)*Mesh%Tnoeud(j)%DLocal(2), 1+Mesh%Tnoeud(j)%DLocal(1)**2]/Mesh%Tnoeud(j)%DLocal(5),(/2,2/))
                    GPhiLocal = DF1(1)*(TenseurMetrique(1,1)*DfDu + TenseurMetrique(1,2)*DfDv) + &
                                                            & DF1(2)*(TenseurMetrique(2,1)*DfDu + TenseurMetrique(2,2)*DfDv)
                
                    ! Retour à la base globale.
                    Ecoulement%GPhi(1:3,j)%perturbation = matmul(Puvw,GPhiLocal)
                
                    ! Ajout de la composante normale.
                    Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
                
                    deallocate(Pvoisin, A, B1, Beta)
                end do
            
                do j = i1,i2
                    
                    ! Gradient du gradient du potentiel par dérivation des dérivées d'ordre 1.
                    Na = Nordre(Mesh%Tnoeud(j)%Ordre)
                    NVoisin = Mesh%Tnoeud(j)%Nvoisin(2)-1
                    allocate(Pvoisin(3,NVoisin+1), A(NVoisin+Na,NVoisin+Na), B(NVoisin+Na,3), ipiv(NVoisin+Na),Beta(NVoisin+Na))
                    Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                    Pt(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
                
                    ! Changement de variables.
                    B = 0._RP
                    do k = 1,NVoisin+1
                        Pk = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%Tvoisin(k,1)))%Pnoeud
                        if (Mesh%Tnoeud(j)%Tvoisin(k,1).lt.0) Pk(2)=-Pk(2)
                        M = Mesh%Tnoeud(j)%Pnoeud
                        Pvoisin(1:3,k) = matmul(Pk - M, Puvw)                        
                        GPhiLocal = Ecoulement%GPhi(1:3,abs(Mesh%Tnoeud(j)%Tvoisin(k,1)))%Perturbation
                        if (Mesh%Tnoeud(j)%Tvoisin(k,1).lt.0) GPhiLocal(2) = - GPhiLocal(2)
                        B(k,1:3) = matmul(GPhiLocal, Puvw)
                    end do
                
                    ! On construit A.
                    call SplineMatrix(Nvoisin, Mesh%Tnoeud(j)%Ordre, Pvoisin, A)
                    ipiv = 0._RP
                    trans='n'
                    nrhs=3 ! Nombre de second membre.
                
                    ! Décomposition LU de A.
                    call dgetrf(Nvoisin+Na,Nvoisin+Na,A,Nvoisin+Na, ipiv,info)
                
                    ! Résolution du système à partir de la décomposition LU.
                    call dgetrs(trans,Nvoisin+Na,nrhs,A,Nvoisin+Na,ipiv,B,Nvoisin+Na,info) ! The solution is in B.
                
                    ! Calcul des dérivées partielles secondes.
                    call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,1), Pvoisin(1:3,1), Pvoisin, DF1)
                    call SplineDF(Nvoisin, Mesh%Tnoeud(j)%Ordre, B(1:Nvoisin+Na,2), Pvoisin(1:3,1), Pvoisin, DF2)
                    Ecoulement%GPhi2(3:4,j)%perturbation = [DF1(1),DF2(2)]
                    Ecoulement%GPhi2(3:4,j)%perturbation = Ecoulement%GPhi2(3:4,j)%perturbation - Mesh%TNoeud(j)%DLocal(3:4)*Ecoulement%DPhiDn(j)%Perturbation        
                
                    deallocate(Pvoisin, A, B, Beta, ipiv)
                end do
            else
                ! Linear Panel approximation.
                do j=i1,i2
                    
                    ! Surface Gradient of Phi.
                    Ecoulement%GPhi(1:3,j)%perturbation = 0._RP
                    Ponderation = 0._RP
                    do k=1,Mesh%Tnoeud(j)%Nfacette
                        Ponderation = Ponderation + Mesh%Tnoeud(j)%Angle(k)
                        do r=1,3
                            Phifacette(r)=Ecoulement%Phi(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation
                        end do
                        Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k) !*0.5_RP
                    end do
                    Ecoulement%GPhi(1:3,j)%perturbation = Ecoulement%GPhi(1:3,j)%perturbation/Ponderation
                    if (symmetry.and.abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon) then
                        Ecoulement%GPhi(2,j)%perturbation = 0._RP
                    end if
                    Ecoulement%GPhi(:,j)%perturbation = Ecoulement%GPhi(:,j)%perturbation + Ecoulement%DPhiDn(j)%perturbation*Mesh%Tnoeud(j)%Normale
                    Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                    Pt(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
                    GPhiLocal = matmul(Pt,Ecoulement%GPhi(:,j)%perturbation)
                    Ecoulement%GPhi2(1:2,j)%perturbation = GPhiLocal(1:2)
                end do
                do j = i1,i2
                    ! Surface Gradient of GPhi.
                    do s=1,3
                        Ecoulement%GradGrad(1:3,s,j)%Perturbation = 0._RP
                        Ponderation = 0._RP
                        do k = 1,Mesh%Tnoeud(j)%Nfacette
                            Ponderation = Ponderation + Mesh%Tnoeud(j)%Angle(k)
                            do r = 1,3
                                Phifacette(r) = Ecoulement%Gphi(s,Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation - Ecoulement%DphiDn(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%perturbation*Mesh%TNoeud(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%Normale(s)
                            end do
                            Ecoulement%GradGrad(1:3,s,j)%Perturbation = Ecoulement%GradGrad(1:3,s,j)%Perturbation + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Phifacette)*Mesh%Tnoeud(j)%Angle(k) !*0.5_RP
                        end do
                        Ecoulement%GradGrad(1:3,s,j)%Perturbation = Ecoulement%GradGrad(1:3,s,j)%Perturbation/Ponderation
                    end do
                    GGPhi = Ecoulement%GradGrad(1:3,1:3,j)%perturbation
                    Ecoulement%GPhi2(3,j)%perturbation = dot_product(matmul(GGPhi,Mesh%Tnoeud(j)%Plocal(1:3,1,1)),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                    Ecoulement%GPhi2(4,j)%perturbation = dot_product(matmul(GGPhi,Mesh%Tnoeud(j)%Plocal(1:3,2,1)),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                end do
            end if
            
        end if
        
    end do
    
end subroutine GradientBody

subroutine GradientIntersection(Mesh, Ecoulement)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(in)         :: Mesh                     ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement), intent(inout)    :: Ecoulement               ! Flow parameters.
    
    integer                             :: j, k2, k3                ! Loop parameters.
    real(rp), dimension(3)              :: N1, N2, N3, matB, Sol    ! Normals, vector B and the solution of the linear system.
    real(rp), dimension(3,3)            :: matA                     ! Matrix A of the linear system.
    logical                             :: bool                     ! True if the normals of the twin nodes are almost equal, false otherwise.
        
    ! This subroutine computes the gradient of the velocity potential on the intersection lines (nothing abour GEta because it is not consistent on the floater panels).
    
    ! Traitement particulier du gradient à l'interface SL/SM.
    do j = 1,Mesh%Nnoeud
        
        if(Mesh%Tnoeud(j)%Ndouble.ne.0)then ! Ajouter conditions pour ne pas doubler le calcul.
            
            bool = .false. ! In case the normals of the twin nodes are equal.
            
            N1 = Mesh%Tnoeud(j)%Normale
            matB(1) = Ecoulement%DPhiDn(j)%perturbation ! GradPhi.N = DPhiDn for the first twin node.
            
            if(Mesh%Tnoeud(j)%Ndouble.ge.1)then
                k2 = Mesh%Tnoeud(j)%Double(1)
                N2 = Mesh%Tnoeud(k2)%Normale
                matB(2) = Ecoulement%DPhiDn(k2)%perturbation ! GradPhi.N = DPhiDn for the second twin node.
                
                if(Mesh%Tnoeud(j)%Ndouble.eq.2)then
                    k3 = Mesh%Tnoeud(j)%Double(2)
                    N3 = Mesh%Tnoeud(k3)%Normale
                    matB(3) = Ecoulement%DPhiDn(k3)%perturbation
                else
                    call Computation_vect_product(N1,N2,N3)
                    N3 = N3 / norm2(N3) ! Unit norm.
                    matB(3) = dot_product(Ecoulement%GPhi(:,j)%perturbation,N3) ! GradPhi.N = GradPhi.N for the two twin nodes.
                end if
                
                ! If N1 and N2 are almost equal (ex: beginning of a hydrodynamic impact with a sphere in still water, the first panels of the sphere are almost horizontal, like the panels of the free surface. The normals at the intersection are the same.
                ! Thus the linear system has two identical equations: GradPhi.N = DPhiDn along the same normals (floater and free surface). The other components are not well defined.
                if(dot_product(N1,N2).gt.DPNormals)then
                    
                    bool = .true. ! Normals almost equal.
                    
                    ! Second condition.
                    N2 = N3
                    matB(2) = dot_product(Ecoulement%GPhi(:,j)%perturbation,N2) ! GradPhi.N = GradPhi.N for the two twin nodes.
                    
                    ! Thid condition.
                    call Computation_vect_product(N1,N2,N3)
                    N3 = N3 / norm2(N3) ! Unit norm.
                    matB(3) = dot_product(Ecoulement%GPhi(:,j)%perturbation,N3) ! GradPhi.N = GradPhi.N for the two twin nodes.
                    
                end if
                    
            end if
            
            ! Systeme lineaire.
            matA(1,1:3) = N1
            matA(2,1:3) = N2
            matA(3,1:3) = N3
            
            ! Resolution.
            if(iprint>0) print*,"GradientIntersection: Resolution systeme linéaire intersection..."
            call LU(matA,matB,Sol,3)    
            if(iprint>0) print*,"GradientIntersection: done"
            
            ! Affectation.
            Ecoulement%GPhi(:,j)%perturbation = Sol(1:3)
            if(not(bool))then
                if (abs(dot_product(N1,Sol) - Ecoulement%DPhiDn(j)%perturbation).gt. Epsilon) then
                    print*, 'GradientIntersection: FS/Body Interface :', dot_product(N1,Sol), Ecoulement%DPhiDn(j)%perturbation
                end if
                if (abs(dot_product(N2,Sol) - Ecoulement%DPhiDn(k2)%perturbation).gt. Epsilon) then
                    print*, 'GradientIntersection: FS/Body Interface :', dot_product(N2,Sol), Ecoulement%DPhiDn(k2)%perturbation
                end if
            end if
            
        end if 
        
    end do
    
end subroutine GradientIntersection

subroutine SplineMatrix(Nvoisin, Ordre, Pvoisin, A)
    !!!!!
    ! Construction of the linear matrix for the determination of the B-Spline coefficients, from the position of the Nvoisin points Pvoisin.
    !!!!
    
    integer, intent(in) :: Nvoisin, Ordre
    real(rp), dimension(3,Nvoisin+1),intent(in) :: Pvoisin
    real(rp),dimension(:,:),intent(inout) :: A
    
    integer :: j,k, Na
    real(rp) :: Norme
    real(rp), dimension(3) :: M1
    
    ! Size of A
    Na = Nordre(Ordre)
    
    ! Construction of A
    do k = 1,NVoisin+1
        do j = 1,Nvoisin+1
            Norme = norm2([Pvoisin(1,j)-Pvoisin(1,k),Pvoisin(2,j)-Pvoisin(2,k),0._RP])
            if(Ordre.eq.0)then
                if(j.eq.k)then
                    A(j,k) = 0._RP
                else
                    A(j,k) = (Norme**2)*Log(Norme)
                end if
            else
                A(j,k) = Norme**(2*Ordre-1)
            end if
        end do
    end do
    
    A(Nvoisin+2:Nvoisin+Na,Nvoisin+2:Nvoisin+Na) = 0._RP
    
    do k = 1,Nvoisin+1
        M1 = Pvoisin(1:3,k)	
        A(k,Nvoisin+2) = 1
        A(k,Nvoisin+3) = M1(1)
        A(k,Nvoisin+4) = M1(2)    	
        A(Nvoisin+2,k) = 1
        A(Nvoisin+3,k) = M1(1)
        A(Nvoisin+4,k) = M1(2)
        if(Ordre.gt.0)then
            A(k,Nvoisin+5) = M1(1)**2
            A(k,Nvoisin+6) = M1(1) * M1(2)
            A(k,Nvoisin+7) = M1(2)**2
            A(Nvoisin+5,k) = M1(1)**2
            A(Nvoisin+6,k) = M1(1) * M1(2)
            A(Nvoisin+7,k) = M1(2)**2
            if(Ordre.eq.4)then
                A(k,Nvoisin+8) = M1(1)**3
                A(k,Nvoisin+9) = M1(2)*M1(1)**2
                A(k,Nvoisin+10) = M1(1)*M1(2)**2
                A(k,Nvoisin+11) = M1(2)**3
                A(Nvoisin+8,k) = M1(1)**3
                A(Nvoisin+9,k) = M1(2)*M1(1)**2
                A(Nvoisin+10,k) = M1(1)*M1(2)**2
                A(Nvoisin+11,k) = M1(2)**3
            end if
        end if	
    end do
    
end subroutine SplineMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                                                       !
!                                               SPLINEDF                                                                                                !
!                                                                                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SplineDF(Nvoisin, Ordre, Alpha, Pref, Pvoisin, DF)
    
    integer, intent(in)                             :: Nvoisin, Ordre   ! Number of neighbours and ordre of the B-splines.
    real(rp), intent(in)                            :: Alpha(:)         ! B-spline coefficients.
    real(rp), dimension(3), intent(in)              :: Pref             ! Position of computation.
    real(rp), dimension(3, Nvoisin+1), intent(in)   :: Pvoisin          ! Neighbours.
    real(rp), dimension(2), intent(out)             :: DF               ! First differentiation, (1) = d/dx, (2) = d/dy.
    
    integer                                         :: j                ! Loop parameter.
    real(rp)                                        :: Norme            ! Norm.
    
    ! This subroutine computes the first differentation of F using B-splines from (Eq 4.37) of LL.
        
    ! Differentiation of the polynom.
    DF(1) = Alpha(Nvoisin+3)
    DF(2) = Alpha(Nvoisin+4)
    if (Ordre.gt.0) then
        DF(1) = DF(1) + 2._RP*Alpha(Nvoisin+5)*Pref(1) + Alpha(Nvoisin+6)*Pref(2)
        DF(2) = DF(2) + 2._RP*Alpha(Nvoisin+7)*Pref(2) + Alpha(Nvoisin+6)*Pref(1)
        if (Ordre.eq.4) then
            DF(1) = DF(1) + 3._RP*Alpha(Nvoisin+8)*Pref(1)**2 + 2._RP*Alpha(Nvoisin+9)*Pref(2)*Pref(1) + Alpha(Nvoisin+10)*Pref(2)**2
            DF(2) = DF(2) + Alpha(Nvoisin+9)*Pref(1)**2 + 2._RP*Alpha(Nvoisin+10)*Pref(1)*Pref(2) + 3._RP*Alpha(Nvoisin+11)*Pref(2)**2
        end if
    end if
    
    ! Differentiation of the sum.
    do j=2,NVoisin+1
        Norme = norm2([Pvoisin(1,j)-Pref(1),Pvoisin(2,j)-Pref(2),0._RP]) ! ||u - ui||
        if (Ordre.eq.0) then
            DF(1) = DF(1) + Alpha(j)*(Pref(1)-Pvoisin(1,j))*(2._RP*Log(Norme) + 1._RP)
            DF(2) = DF(2) + Alpha(j)*(Pref(2)-Pvoisin(2,j))*(2._RP*Log(Norme) + 1._RP)
        else
            DF(1) = DF(1) + (2._RP*Ordre-1._RP)*Alpha(j)*(Pref(1)-Pvoisin(1,j))*Norme**(2*Ordre-3)
            DF(2) = DF(2) + (2._RP*Ordre-1._RP)*Alpha(j)*(Pref(2)-Pvoisin(2,j))*Norme**(2*Ordre-3)
        end if
    end do
    
end subroutine SplineDF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                                                       !
!                                                 SPLINED2F                                                                                             !
!                                                                                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SplineD2F(Nvoisin, Ordre, Alpha, Pref, Pvoisin, D2F)
    !!!!!
    ! Approximation of the second derivatives DF on the point Pref, from the B-Spline coefficients Alpha.
    !!!!
    ! Parameters
    integer, intent(in) :: Nvoisin, Ordre
    real(rp), intent(in) :: Alpha(:)
    real(rp), dimension(3), intent(in) :: Pref
    real(rp), dimension(3, Nvoisin+1), intent(in) :: Pvoisin
    real(rp), dimension(3), intent(out) :: D2F ! Derivatives d²/du², d²/dv² and d²/dudv
    ! Locals
    integer :: j
    real(rp) :: Norme, inv_Norme
    ! Begin
    D2F = 0._RP
    if (Ordre.gt.0) then
        D2F(1) = D2F(1) + 2._RP*Alpha(Nvoisin+5)
        D2F(2) = D2F(2) + 2._RP*Alpha(Nvoisin+7)
        D2F(3) = D2F(3) + Alpha(Nvoisin+6)    
        if (Ordre.eq.4) then
            D2F(1) = D2F(1) + 6._RP*Alpha(Nvoisin+8)*Pref(1) + 2._RP*Alpha(Nvoisin+9)*Pref(2)
            D2F(2) = D2F(2) + 2._RP*Alpha(Nvoisin+10)*Pref(1) + 6._RP*Alpha(Nvoisin+11)*Pref(2)
            D2F(3) = D2F(3) + Alpha(Nvoisin+6) + 2._RP*Alpha(Nvoisin+9)*Pref(1) + 2._RP*Alpha(Nvoisin+10)*Pref(2)
        end if
    end if
    do j=2,NVoisin+1
        Norme = norm2([Pvoisin(1,j)-Pref(1),Pvoisin(2,j)-Pref(2),0._RP])
        inv_Norme = 1._RP/Norme**2
        if (Ordre.eq.0) then
            D2F(1) = D2F(1) + Alpha(j)*( 2._RP*Log(Norme) + 1._RP + inv_Norme*(Pref(1)-Pvoisin(1,j))**2)
            D2F(2) = D2F(2) + Alpha(j)*( 2._RP*Log(Norme) + 1._RP + inv_Norme*(Pref(2)-Pvoisin(2,j))**2)
            D2F(3) = D2F(3) + Alpha(j)*(Pref(1)-Pvoisin(1,j))*(Pref(2)-Pvoisin(2,j))*inv_Norme
        else
            D2F(1) = D2F(1) + (2._RP*Ordre-1._RP)*Alpha(j)*Norme**(2*Ordre-3) + (2._RP*Ordre-1._RP)*(2._RP*Ordre-3._RP)*0.5_RP*Alpha(j)*(Pref(1)-Pvoisin(1,j))*Norme**(2*Ordre-5)
            D2F(2) = D2F(2) + (2._RP*Ordre-1._RP)*Alpha(j)*Norme**(2*Ordre-3) + (2._RP*Ordre-1._RP)*(2._RP*Ordre-3._RP)*0.5_RP*Alpha(j)*(Pref(2)-Pvoisin(2,j))*Norme**(2*Ordre-5)
            D2F(3) = D2F(3) + (2._RP*Ordre-1._RP)*(2._RP*Ordre-3._RP)*Alpha(j)*(Pref(1)-Pvoisin(1,j))*(Pref(2)-Pvoisin(2,j))*Norme**(2*Ordre-5)
        end if
    end do
    ! End
end subroutine SplineD2F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                                                       !
!                                                 SPLINEF                                                                                               !
!                                                                                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SplineF(Nvoisin, Ordre, Alpha, Pref, Pvoisin, F)
    !!!!!
    ! Approximation of the value F on the point Pref, from the B-Spline coefficients Alpha.
    !!!!
    
    integer, intent(in)                             :: Nvoisin, Ordre
    real(rp), intent(in)                            :: Alpha(:)
    real(rp), dimension(3), intent(in)              :: Pref
    real(rp), dimension(3, Nvoisin+1), intent(in)   :: Pvoisin
    real(rp), intent(out)                           :: F
    
    integer                                         :: j
    real(rp)                                        :: Norme
    
    ! This subroutine computes F using B-splines from (Eq 4.38) or (4.44) of LL. (?)
    
    ! Polynom.
    F = Alpha(Nvoisin+2)+Alpha(Nvoisin+3)*Pref(1)+Alpha(Nvoisin+4)*Pref(2)
    if (Ordre.gt.0) then
        F = F + Alpha(Nvoisin+5)*Pref(1)**2 + Alpha(Nvoisin+6)*Pref(2)*Pref(1) + Alpha(Nvoisin+7)*Pref(2)**2
        if (Ordre.eq.4) then
            F = F + Alpha(Nvoisin+8)*Pref(1)**3 + Alpha(Nvoisin+9)*Pref(2)*Pref(1)**2 + Alpha(Nvoisin+10)*Pref(1)*Pref(2)**2 + Alpha(Nvoisin+11)*Pref(2)**3
        end if
    end if
    
    ! Sum
    do j = 2,NVoisin+1
        Norme = norm2([Pvoisin(1,j)-Pref(1),Pvoisin(2,j)-Pref(2),0._RP])
        if (Ordre.eq.0) then
            F = F + Alpha(j)*(Norme**2)*Log(Norme)
        else
            F = F + Alpha(j)*Norme**(2*Ordre-1)
        end if
    end do
    
end subroutine SplineF

subroutine gradA(indice,Mesh,A,size_A)
    !! Paramètres 
    integer, intent(in) :: indice
    integer,intent(inout) :: size_A
    !real(rp),dimension(:,:), intent(inout) :: A
    real(rp),dimension(size_A,size_A),intent(inout) :: A
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh
    
    ! Variables Locales
    integer :: j, k, NVoisin, Ordre
    real(rp), allocatable :: Voisin(:)
    real(rp), dimension(3) :: vecprojete, M1, M2
    real(rp) :: Norme
    ! Début
    NVoisin = Mesh.Tnoeud(indice).Nvoisin(2)-1		
    Ordre = Mesh.Tnoeud(indice).Ordre 
    allocate(Voisin(0:NVoisin))
    Voisin(0:NVoisin) = Mesh.Tnoeud(indice).TVoisin(1:NVoisin+1,1)
    ! On construit A
    do j=1,NVoisin+1
        M1 = Mesh.Tnoeud(abs(Voisin(j-1))).Pnoeud
        if (Voisin(j-1).lt.0) M1(2) = -M1(2)
	    do k=1,Nvoisin+1
	        M2 = Mesh.Tnoeud(abs(Voisin(k-1))).Pnoeud
            if (Voisin(k-1).lt.0) M2(2) = -M2(2)
	        vecprojete = M2 - M1
		    vecprojete(3) = 0._RP			
		    Norme = norm2(vecprojete)
            if (k.eq.j) then 
	            A(k,j) = 0
            else
                if (Norme.lt.epsilon) print*, Voisin(j-1), Voisin(k-1)
                if (Ordre==0) then
	                A(k,j)=(Norme**2) *Log(Norme)
	            elseif (Ordre==3) then
	                A(k,j)= Norme**5
	            else
	                print*, 'erreur gradA dans Ordre des Splines : Ordre = ', Ordre
	            end if        
            end if
	    end do
    end do
    do j=Nvoisin+2,Nvoisin+Ordre+4
	    do k=Nvoisin+2,Nvoisin+Ordre+4	
	    A(k,j) = 0
	    end do
    end do
    if (Ordre==0) then
        do j=1,Nvoisin+1
            M1 = Mesh.Tnoeud(abs(Voisin(j-1))).Pnoeud
            if (Voisin(j-1).lt.0) M1(2) = -M1(2)
	        A(j,Nvoisin+2) = 1
	        A(j,Nvoisin+3) = M1(1)
	        A(j,Nvoisin+4) = M1(2)
	        A(Nvoisin+2,j) = 1
	        A(Nvoisin+3,j) = M1(1)
	        A(Nvoisin+4,j) = M1(2)
        end do
    elseif (Ordre==3) then
        do j=1,Nvoisin+1
	        M1 = Mesh.Tnoeud(abs(Voisin(j-1))).Pnoeud
            if (Voisin(j-1).lt.0) M1(2) = -M1(2)	
	        A(j,Nvoisin+2) = 1
	        A(j,Nvoisin+3) = M1(1)
	        A(j,Nvoisin+4) = M1(2)
	        A(j,Nvoisin+5) = M1(1) * M1(1)
	        A(j,Nvoisin+6) = M1(1) * M1(2)
	        A(j,Nvoisin+7) = M1(2) * M1(2)	
	        A(Nvoisin+2,j) = 1
	        A(Nvoisin+3,j) = M1(1)
	        A(Nvoisin+4,j) = M1(2)
	        A(Nvoisin+5,j) = M1(1) * M1(1)
	        A(Nvoisin+6,j) = M1(1) * M1(2)
	        A(Nvoisin+7,j) = M1(2) * M1(2)	
        end do
    end if
    deallocate(Voisin)
    ! Fin
end subroutine gradA

!-------------------------------------------------------------
!
!	Calcul des Coefficients de Spline, selon l'ordre de la Spline voulu
!                       
!-------------------------------------------------------------
subroutine CoeffSpline(indice, Ecoulement, Mesh, Alpha, size_Alpha,Beta,size_Beta,Opt)
! Paramètres
integer, intent(in) :: indice
!f2py integer*1, dimension(1000) :: Ecoulement
type(TEcoulement), intent(in) :: Ecoulement
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
integer,intent(in) :: size_Alpha
real(rp), intent(out) :: Alpha(size_Alpha)
integer,intent(in) :: size_Beta
real(rp),intent(out) :: Beta(size_Beta)
logical, optional, intent(in) :: Opt
! Variables
integer :: j, Nvoisin, Ordre
integer :: size_A
real(rp),dimension(:,:),allocatable :: B,A
real(rp),dimension(:),allocatable :: Voisin
!!!! LU Solver variables:
character(len=1) :: trans
integer, allocatable :: ipiv(:)
integer :: info, lda, ldb, nrhs
!Début
NVoisin = Mesh.Tnoeud(indice).Nvoisin(2)-1		
Ordre = Mesh.Tnoeud(indice).Ordre 
size_A = Nvoisin+Ordre+4
allocate(Voisin(0:NVoisin), B(NVoisin+Ordre+4,2), A(size_A,size_A))!, Alpha(0:Nvoisin+Ordre+3))
allocate(ipiv(Nvoisin+Ordre+4))
Voisin(0:NVoisin) = Mesh.Tnoeud(indice).TVoisin(1:NVoisin+1,1)
!On construit B
do j=1,Nvoisin+1
    B(j,1) = Ecoulement%Eta(abs(Voisin(j-1)))%perturbation
end do
nrhs=1 ! nombre de second membre
if (present(Opt)) then
    if (Opt) then
        !allocate(Beta(0:Nvoisin+Ordre+3))
        do j=1,Nvoisin+1
            B(j,2) = Ecoulement%Phi(abs(Voisin(j-1)))%perturbation
        end do
        nrhs = 2
    end if
end if
B(Nvoisin+2:Nvoisin+Ordre+4,:) = 0._RP

call gradA(indice, Mesh, A,size_A)

! Résolution du Système Linéaire
lda=NVoisin+Ordre+4 ! leading dimension de A
ldb=lda ! leading dimension de B
ipiv=0 
trans='n'
! Décomposition LU de A
if(iprint>0) print*,"Decomposition LU de A pour gradA..."
call dgetrf(lda,lda,A,lda, ipiv, info) 
if (info.ne.0) print*, '       info =', info, indice
if(info.ne.0)then
  print*," indice  = ",indice," ordre = ",Ordre," Nvoisin = ",Nvoisin
endif
! Résolution du système à partir de la décomposition Lu
call dgetrs(trans,NVoisin+Ordre+4,nrhs,A,lda,ipiv,B,ldb, info)
if (info.ne.0) print*, '       info =', info, indice

Alpha(1:NVoisin+Ordre+4) = B(1:NVoisin+Ordre+4,1)
Beta(1:NVoisin+Ordre+4) = B(1:NVoisin+Ordre+4,2)

!!print*,"done"
deallocate(Voisin, B, A, ipiv)
!! Fin
end subroutine CoeffSpline
!-------------------------------------------------------------
!
!	Gradients surfaciques de Phi et Eta, par approximation par BSplines
!                       
!-------------------------------------------------------------
subroutine gradspline(indice, Ecoulement, Mesh, Opt)
! Paramètres
integer, intent(in) :: indice
!f2py integer*1, dimension(1000) :: Ecoulement
type(TEcoulement), intent(inout):: Ecoulement
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
! Permet de calculer le gradient de Eta plus en option celui de Phi: Opt = 0 --> calcul que de grad(Eta), Opt = 1 --> calcul des deux
logical, optional, intent(in) :: Opt 
! Variables Locales
integer :: Nvoisin, Ordre
real(rp), dimension(2,2) :: TenseurMetrique
real(rp), dimension(3) :: DfDu, DfDv
real(rp) :: delta
real(rp), dimension(2) :: DF
real(rp), dimension(:), allocatable :: Alpha, Beta
integer ::  size_Alpha, size_Beta
! Début
NVoisin = Mesh%Tnoeud(indice)%NVoisin(2)-1
Ordre = Mesh%Tnoeud(indice)%Ordre
size_Alpha = NVoisin+Ordre+4
size_Beta = NVoisin+Ordre+4
!allocate(Voisin(0:NVoisin))
allocate(Alpha(size_Alpha), Beta(size_Beta))

!Voisin(0:NVoisin) = Mesh.Tnoeud(indice).TVoisin(1:NVoisin+1,1)
! Calcul des Coefficients des BSplines
call CoeffSpline(indice, Ecoulement, Mesh, Alpha, size_Alpha, Beta,size_Beta,Opt)

! Calcul du gradient surfacique de Eta
if (Ordre==0) then    
    call GradSplinePM(indice, Mesh, Alpha, DF)
elseif (Ordre==3) then
    call GradSplinePP3(indice, Mesh, Alpha, DF)
else
    print*, 'erreur GradSpline dans Ordre des Splines : Ordre = ', Ordre, 'test'
end if

Ecoulement%GEta(1:2,indice)%perturbation = DF
Ecoulement%GEta(3,indice)%perturbation = 0
! Calcul du gradient surfacique de Phi
if (present(Opt)) then
    if (Opt) then
        delta=(1+Ecoulement%GEta(1,indice)%incident**2)*(1+Ecoulement%GEta(2,indice)%incident**2)-(Ecoulement%GEta(1,indice)%incident*Ecoulement%GEta(2,indice)%incident)**2
        DfDu=[1._RP,0._RP,Ecoulement%GEta(1,indice)%incident]
        DfDv=[0._RP,1._RP,Ecoulement%GEta(2,indice)%incident]
        ! Calcul temporaire de la normale aux noeuds
!        Mesh%Tnoeud(indice).NormaleS=vect_product(DfDu,DfDv)/norm2(vect_product(DfDu,DfDv))
        TenseurMetrique=reshape([1+Ecoulement%GEta(2,indice)%incident**2, -Ecoulement%GEta(1,indice)%incident*Ecoulement%GEta(2,indice)%incident,&
                     -Ecoulement%GEta(1,indice)%incident*Ecoulement%GEta(2,indice)%incident, 1+Ecoulement%GEta(1,indice)%incident**2]/delta,(/2,2/))
        if (Ordre==0) then
            call GradSplinePM(indice, Mesh, Beta, DF)
        elseif (Ordre==3) then
            call GradSplinePP3(indice, Mesh, Beta, DF) 
        end if
        Ecoulement%GPhi(:,indice)%perturbation = DF(1)*(TenseurMetrique(1,1)*DfDu + TenseurMetrique(1,2)*DfDv) + &
                                                & DF(2)*(TenseurMetrique(2,1)*DfDu + TenseurMetrique(2,2)*DfDv)
    end if  
end if
!deallocate(Voisin)
deallocate(Alpha,Beta)
! Fin
end subroutine gradspline

!-------------------------------------------------------------
!
!	Gradient surfacique pour BSpline Plaque Mince
!
!-------------------------------------------------------------
subroutine GradSplinePM(indice, Mesh, Alpha, DF, F)
! Paramètres 
integer, intent(in) :: indice
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
real(rp), dimension(0:), intent(in) :: Alpha
real(rp), dimension(2), intent(out) :: DF
real(rp), intent(out), optional :: F
! Variables Locales
integer :: j, Nvoisin
real(rp) :: Norme
real(rp), dimension(3) :: Vecprojete, M1
integer, allocatable :: Voisin(:)
! Début
NVoisin = Mesh%Tnoeud(indice)%NVoisin(2)-1
allocate(Voisin(0:NVoisin))
Voisin(0:NVoisin) = Mesh%Tnoeud(indice)%TVoisin(1:NVoisin+1,1)
DF(1)=Alpha(Nvoisin+2)
DF(2)=Alpha(Nvoisin+3)
if (present(F)) F = Alpha(Nvoisin+1) + Alpha(Nvoisin+2)*Mesh%Tnoeud(indice)%Pnoeud(1) + Alpha(Nvoisin+3)*Mesh%Tnoeud(indice)%Pnoeud(2)
do j=1,Nvoisin
    M1 = Mesh%Tnoeud(abs(Voisin(j)))%Pnoeud
    if (Voisin(j).lt.0) M1(2) = -M1(2)
    vecprojete = M1 - Mesh%Tnoeud(indice)%Pnoeud
	vecprojete(3) = 0._RP
	Norme = norm2(vecprojete)
!	if (Norme.lt.Epsilon) print*, indice, Voisin(j-1), j
!	Norme = sqrt(vecprojete(1)**2+vecprojete(2)**2)
    DF(1) = DF(1) + Alpha(j)*(Mesh%Tnoeud(indice)%Pnoeud(1)-M1(1))*(2*log(Norme)+1)
    DF(2) = DF(2) + Alpha(j)*(Mesh%Tnoeud(indice)%Pnoeud(2)-M1(2))*(2*log(Norme)+1)	    
    if (present(F)) F = F  + Alpha(j)*log(Norme)*Norme**2
end do
deallocate(Voisin)
! Fin
end subroutine GradSplinePM
!-------------------------------------------------------------
!
!	Gradient surfacique pour BSpline Pseudo Polynomiale d'ordre 3
!
!-------------------------------------------------------------
subroutine GradSplinePP3(indice, Mesh, Alpha, DF, F)
! Paramètres 
integer, intent(in) :: indice
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage), intent(in) :: Mesh
real(rp), dimension(0:), intent(in) :: Alpha
real(rp), dimension(2), intent(out) :: DF
real(rp), intent(out), optional :: F
! Variables Locales
integer :: j, Nvoisin
real(rp) :: Norme
real(rp), dimension(3) :: Vecprojete, M1
integer, allocatable :: Voisin(:)
NVoisin = Mesh%Tnoeud(indice)%NVoisin(2)-1
allocate(Voisin(0:NVoisin))
Voisin(0:NVoisin) = Mesh%Tnoeud(indice)%TVoisin(1:NVoisin+1,1)

DF(1)=Alpha(Nvoisin+2)+2*Alpha(Nvoisin+4)*Mesh%Tnoeud(indice)%Pnoeud(1)+Alpha(Nvoisin+5)*Mesh%Tnoeud(indice)%Pnoeud(2)
DF(2)=Alpha(Nvoisin+3)+2*Alpha(Nvoisin+6)*Mesh%Tnoeud(indice)%Pnoeud(2)+Alpha(Nvoisin+5)*Mesh%Tnoeud(indice)%Pnoeud(1)
if (present(F)) F = Alpha(Nvoisin+1) + Alpha(Nvoisin+2)*Mesh%Tnoeud(indice)%Pnoeud(1) + Alpha(Nvoisin+3)*Mesh%Tnoeud(indice)%Pnoeud(2)&
& + Alpha(Nvoisin+4)*Mesh%Tnoeud(indice)%Pnoeud(1)**2 + Alpha(Nvoisin+6)*Mesh%Tnoeud(indice)%Pnoeud(2)**2&
& + Alpha(Nvoisin+5)*Mesh%Tnoeud(indice)%Pnoeud(1)*Mesh%Tnoeud(indice)%Pnoeud(2)
do j=1,Nvoisin
    M1 = Mesh%Tnoeud(abs(Voisin(j)))%Pnoeud
    if (Voisin(j).lt.0) M1(2) = -M1(2)
    vecprojete = M1 - Mesh%Tnoeud(indice)%Pnoeud
	vecprojete(3) = 0._RP
	Norme = norm2(vecprojete)
    DF(1) = DF(1) + 5*Alpha(j)*(Mesh%Tnoeud(indice)%Pnoeud(1)-M1(1))*Norme**3
    DF(2) = DF(2) + 5*Alpha(j)*(Mesh%Tnoeud(indice)%Pnoeud(2)-M1(2))*Norme**3
    if (present(F)) F = F  + Alpha(j)*Norme**5
end do
deallocate(Voisin)
! Fin
end subroutine GradSplinePP3

end module Spline
