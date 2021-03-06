!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                   !
!                                   Calcul des Coefficients d'Influence                             !
!                                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module CoeffInfluence
use FonctionsCommunes
!$ use OMP_LIB
implicit none

type TVariableChange
    real(rp), dimension(3,3) :: Puvw
    real(rp)                 :: Jacobien, Delta
end type TVariableChange

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                           Calcul Partiel des Coefficients d'Influence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CoeffInfl(Mesh, CD, CS, Nnodes, bornes)
    !!!!! Probl�me :
    !   On ne calcule les coefficients d'influence que d'une partie du maillage sur une autre partie du Maillage%
    !   Ces parties sont d�limit�es en noeuds et facettes par les bornes des indices:
    !       o borne(1,1) : indice du 1er noeud du maillage concern�, 
    !       o borne(1,2) : indice du dernier noeud du maillage concern�
    !       o borne(2,1) : indice de la 1ere facette du maillage concern�,
    !       o borne(2,2) : indice de la derniere facette du maillage concern�
    !!!!! 
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage), intent(in)                         :: Mesh                                     ! Mesh.
    integer,intent(in)                                  :: Nnodes                                   ! Number of nodes in the mesh.
    real(rp), dimension(Nnodes,Nnodes), intent(inout)   :: CD,CS                                    ! Influence coefficient matrices.
    integer, dimension(2,2), intent(in), optional       :: bornes                                   ! Boundaries for the computation of the influence coefficients (nodes = rows and panels != columns !!!).
                
    integer                                             :: j,k, r, jtemp, NDouble , js,ns,js1,ns1   ! Loop parameters, number of twin nodes.
    integer, dimension(2,2)                             :: borne                                    ! Boundaries.
    real(rp)                                            :: Ssigma, Smu                              ! Integrals.
    real(rp), dimension(3)                              :: Isigma, Imu, Ar, TDouble                 ! Integrals.
    real(rp), dimension(3)                              :: Css, Cdd                                 ! Influence coefficient.
    real(rp)                                            :: Cr, ZG, invRho                           ! Parameters.
    real(rp), dimension(3)                              :: MG, M, M0                                ! Points.
    !f2py integer*1, dimension(1000)                    :: Facette
    type(TFacette)                                      :: Facette                                  ! Panel.
    real(rp)                                            :: Crmax                                    ! Distance criteria for the asymptotic development.
    real(rp),parameter                                  :: inv3 = 0.3333333333333333_rp             ! 1/3.
       

    real(rp) :: ostart, oend ! Pour chronometerer en openmp
    
    !$ ostart = omp_get_wtime()

    ! This subroutine was previously used. But it still works.
    ! This subroutine computes the influence coefficients from each node on each panel between the boundaries bornes.
    
    if(idebug>1) print*,'CoeffInfl'
        
    ! Boundaries for the nodes.
    borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    if(present(bornes)) borne = bornes
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    endif
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    endif


    
    !$OMP PARALLEL DO NUM_THREADS(NThreads) SHARED(Mesh,borne,CD,CS,ns,ns1,Ldom) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC,(borne(1,2)-borne(1,1)+1)/Nthreads)
    do j = borne(1,1),borne(1,2) ! 1st node, last node.

        
        NDouble = Mesh%Tnoeud(j)%Ndouble ! Number of twin nodes.
        TDouble = Mesh%Tnoeud(j)%double ! Table of twin nodes.
               
        if(NDouble.eq.0 .or. (NDouble.ne.0.and.j.lt.TDouble(1)) .or. present(bornes))then
            
            M0 = Mesh%Tnoeud(j)%Pnoeud
            
            do k = borne(2,1),borne(2,2) ! 1st panel, last panel.
                
            
            
                  
                    if (j ==1 .and. (k==317 .or. k==318 .or. k==379)) then
                
                        print_debug = 1
                                
                    else
                        print_debug = 0
                                
                          
                    end if
                            
                            
                            
                M = M0
                Facette = Mesh%Tfacette(k) 
                Crmax = 64._rp*Facette%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.
                Css = 0._rp ; Cdd = 0._rp
                
                
                ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
                do js1 = 1,-ns1,-2
                                                        
                    M(2) = js1*M(2)
                    
                    ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
                    do js = 1,-ns,-2
         
                        M(3) = js*M(3)-(1._rp-js)*Ldom(3)
                        MG = M-Facette%Gfacette                   
                        
                        ! Distance criteria to use or not an symptotic development.
                        Cr = dot_product(MG,MG) ! Cr = GM^2
                        
                        

                        
                        
                        
                        
                        if(Cr.lt.Crmax)then ! Exact computation.

                            ! Exact computation.
                           
                            call Iints(M, Facette, Isigma, Imu) ! Line integrations.
                            call Sints(M, Facette, Ssigma, Smu) ! Surface integrations.
                            
                            

                            
                            ! Prise en compte du param�trage de la facette.
                            Ar = [inv3,inv3,inv3] + matmul(Facette%dsT,MG)
                            Css = Css + Ssigma*Ar + matmul(Facette%dsT,Isigma)
                            Cdd = Cdd + Smu*Ar + matmul(Facette%dsT,Imu)
                            

                            
                            
                        else ! Asymptotic development.
                        
                            invRho = 1._RP/norm2(MG) ! invRho = 1/||GM||
                            ZG = dot_product(MG,Facette%Normale) ! ZG = GM.n
                            
                            ! I_sigma (Eq 3.53 of LL).
                            Css = Css + Facette%Aire*invRho*inv3
                            
                            ! I_mu (Eq 3.54 of LL).
                            Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
 
                        end if
                        
                    end do   
                end do


                
                ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).                
                CS(j,Facette%Tnoeud(1)) = CS(j,Facette%Tnoeud(1)) + Css(1)
                CS(j,Facette%Tnoeud(2)) = CS(j,Facette%Tnoeud(2)) + Css(2)
                CS(j,Facette%Tnoeud(3)) = CS(j,Facette%Tnoeud(3)) + Css(3)
                
                CD(j,Facette%Tnoeud(1)) = CD(j,Facette%Tnoeud(1)) + Cdd(1)
                CD(j,Facette%Tnoeud(2)) = CD(j,Facette%Tnoeud(2)) + Cdd(2)
                CD(j,Facette%Tnoeud(3)) = CD(j,Facette%Tnoeud(3)) + Cdd(3)
                
                ! Traitement des points doubles pour �viter redondance dans le calcul des Coeff d'Infl.
                if(NDouble.ne.0.and.j.lt.TDouble(1).and..not.present(bornes))then
                    do r = 1,Mesh%TNoeud(j)%Ndouble
                        jtemp = TDouble(r)
                        CS(jtemp,Facette%Tnoeud(1)) = CS(jtemp,Facette%Tnoeud(1)) + Css(1)
                        CS(jtemp,Facette%Tnoeud(2)) = CS(jtemp,Facette%Tnoeud(2)) + Css(2)
                        CS(jtemp,Facette%Tnoeud(3)) = CS(jtemp,Facette%Tnoeud(3)) + Css(3)
                        
                        CD(jtemp,Facette%Tnoeud(1)) = CD(jtemp,Facette%Tnoeud(1)) + Cdd(1)
                        CD(jtemp,Facette%Tnoeud(2)) = CD(jtemp,Facette%Tnoeud(2)) + Cdd(2)
                        CD(jtemp,Facette%Tnoeud(3)) = CD(jtemp,Facette%Tnoeud(3)) + Cdd(3)
                        
                    end do
                end if   
            end do            
        end if
    end do
    !$OMP END PARALLEL DO

    
    !$ oend = omp_get_wtime()
    
    !write(*,*) "Temps du calcul des coefficients d'influence : ", oend - ostart
    
    ! Calcul de l'Angle Solide
    if(not(present(bornes)))then
        call Angle_solid(CD,Mesh%Nnoeud)
    end if
    
end subroutine CoeffInfl

subroutine CoeffInfl_v2(Mesh, CD, CS, Nnodes, bornes)
    !!!!! Probl�me :
    !   On ne calcule les coefficients d'influence que d'une partie du maillage sur une autre partie du Maillage%
    !   Ces parties sont d�limit�es en noeuds et facettes par les bornes des indices:
    !       o borne(1,1) : indice du 1er noeud du maillage concern�, 
    !       o borne(1,2) : indice du dernier noeud du maillage concern�
    !       o borne(2,1) : indice de la 1ere facette du maillage concern�,
    !       o borne(2,2) : indice de la derniere facette du maillage concern�
    !!!!! 
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage), intent(in)                         :: Mesh                                     ! Mesh.
    integer,intent(in)                                  :: Nnodes                                   ! Number of nodes in the mesh.
    real(rp), dimension(Nnodes,Nnodes), intent(inout)   :: CD,CS                                    ! Influence coefficient matrices.
    integer, dimension(2,2), intent(in), optional       :: bornes                                   ! Boundaries for the computation of the influence coefficients (nodes and panels).
            
    integer                                             :: j,k, r, jtemp, NDouble,p,js,ns,js1,ns1   ! Loop parameters, number of twin nodes.
    integer, dimension(2,2)                             :: borne                                    ! Boundaries.
    real(rp)                                            :: Ssigma, Smu                              ! Integrals.
    real(rp), dimension(3)                              :: Isigma, Imu, Ar, TDouble                 ! Integrals.
    real(rp), dimension(3)                              :: Css, Cdd                                 ! Influence coefficient.
    real(rp)                                            :: Cr, ZG, invRho                           ! Parameters.
    real(rp), dimension(3)                              :: MG, M, M0                                ! Points.
    !f2py integer*1, dimension(1000)                    :: Facette
    type(TFacette)                                      :: Facette                                  ! Panel.
    real(rp)                                            :: Crmax                                    ! Distance criteria for the asymptotic development.
    real(rp),parameter                                  :: inv3 = 0.3333333333333333_rp             ! 1/3.
    real(rp),dimension(:),allocatable                   :: col1_CD,col2_CD,col3_CD                  ! Three columns per panel to store the coefficients of CD.
    real(rp),dimension(:),allocatable                   :: col1_CS,col2_CS,col3_CS                  ! Three columns per panel to store the coefficients of CS.
    real(rp),dimension(:),allocatable                   :: diag                                     ! Vector to store the diagonal coefficients of CD (solid angles).
    real(rp), dimension(:,:),allocatable                :: CD_tmp,CS_tmp                            ! Temporary influence coefficient matrices.
            
    ! This subroutine computes the influence coefficients from each panel on each node between the boundaries bornes from .
    
    if(idebug>1) print*,'CoeffInfl'
        
    ! Boundaries.
    borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    if(present(bornes)) borne = bornes
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    endif
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    endif
    
    !if(not(present(bornes)))then
    !    allocate(diag(Mesh%Nnoeud))
    !    diag = 0._RP
    !end if
    
    do k = borne(2,1),borne(2,2) ! 1st panel, last panel.
        
        Facette = Mesh%Tfacette(k)
        Crmax = 64._rp*Mesh%Tfacette(k)%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.
                
        !allocate(col1_CD(borne(1,2)-borne(1,1)+1),col2_CD(borne(1,2)-borne(1,1)+1),col3_CD(borne(1,2)-borne(1,1)+1))
        !allocate(col1_CS(borne(1,2)-borne(1,1)+1),col2_CS(borne(1,2)-borne(1,1)+1),col3_CS(borne(1,2)-borne(1,1)+1))
        !
        !col1_CD = 0._RP
        !col2_CD = 0._RP
        !col3_CD = 0._RP
        !
        !col1_CS = 0._RP
        !col2_CS = 0._RP
        !col3_CS = 0._RP
        
        allocate(CD_tmp(Nnodes,Nnodes),CS_tmp(Nnodes,Nnodes))
        CD_tmp = 0._RP
        CS_tmp = 0._RP
        
        do j = borne(1,1),borne(1,2) ! 1st node, last node.
            
            NDouble = Mesh%Tnoeud(j)%Ndouble ! Number of twin nodes.
            TDouble = Mesh%Tnoeud(j)%double ! Table of twin nodes.
            M = Mesh%Tnoeud(j)%Pnoeud
            
            if(NDouble.eq.0 .or. (NDouble.ne.0.and.j.lt.TDouble(1)) .or. present(bornes))then
                                
                Css = 0._rp ; Cdd = 0._rp
                
                ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
                do js1 = 1,-ns1,-2
                    
                    M(2) = js1*M(2)
                    
                    ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
                    do js = 1,-ns,-2    
                
                        M(3) = js*M(3)-(1._rp-js)*Ldom(3)    
                        MG = M - Facette%Gfacette
                        
                        ! Distance criteria to use or not an symptotic development.
                        Cr = dot_product(MG,MG) ! Cr = GM^2
                        
                        if(Cr.lt.Crmax)then ! Exact computation.
                            
                            ! Exact computation.
                            call Iints(M, Facette , Isigma, Imu) ! Line integrations.
                            call Sints(M, Facette , Ssigma, Smu) ! Surface integrations.
                            
                            ! Prise en compte du param�trage de la facette.
                            Ar = [inv3,inv3,inv3] + matmul(Facette%dsT,MG)
                            Css = Css + Ssigma*Ar + matmul(Facette%dsT,Isigma)
                            Cdd = Cdd + Smu*Ar + matmul(Facette%dsT,Imu)
                            
                        else ! Asymptotic development.
                            
                            invRho = 1._RP/norm2(MG) ! invRho = 1/||GM||
                            ZG = dot_product(MG,Facette%Normale) ! ZG = GM.n
                            
                            ! I_sigma (Eq 3.53 of LL).
                            Css = Css + Facette%Aire*invRho*inv3
                            
                            ! I_mu (Eq 3.54 of LL).
                            Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
                        
                        end if
                                 
                    end do
                    
                end do
                
                ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).
                !CS(j,Mesh%Tfacette(k)%Tnoeud(1)) = CS(j,Mesh%Tfacette(k)%Tnoeud(1)) + Css(1)
                !CS(j,Mesh%Tfacette(k)%Tnoeud(2)) = CS(j,Mesh%Tfacette(k)%Tnoeud(2)) + Css(2)
                !CS(j,Mesh%Tfacette(k)%Tnoeud(3)) = CS(j,Mesh%Tfacette(k)%Tnoeud(3)) + Css(3)
                !
                !CD(j,Mesh%Tfacette(k)%Tnoeud(1)) = CD(j,Mesh%Tfacette(k)%Tnoeud(1)) + Cdd(1)
                !CD(j,Mesh%Tfacette(k)%Tnoeud(2)) = CD(j,Mesh%Tfacette(k)%Tnoeud(2)) + Cdd(2)
                !CD(j,Mesh%Tfacette(k)%Tnoeud(3)) = CD(j,Mesh%Tfacette(k)%Tnoeud(3)) + Cdd(3)
                
                !p = j - borne(1,1) + 1
                !col1_CS(p) = col1_CS(p) + Css(1)
                !col2_CS(p) = col2_CS(p) + Css(2)
                !col3_CS(p) = col3_CS(p) + Css(3)
                !
                !col1_CD(p) = col1_CD(p) + Cdd(1)
                !col2_CD(p) = col2_CD(p) + Cdd(2)
                !col3_CD(p) = col3_CD(p) + Cdd(3)
                
                CS(j,Mesh%Tfacette(k)%Tnoeud(1)) = CS_tmp(j,Mesh%Tfacette(k)%Tnoeud(1)) + Css(1)
                CS(j,Mesh%Tfacette(k)%Tnoeud(2)) = CS_tmp(j,Mesh%Tfacette(k)%Tnoeud(2)) + Css(2)
                CS(j,Mesh%Tfacette(k)%Tnoeud(3)) = CS_tmp(j,Mesh%Tfacette(k)%Tnoeud(3)) + Css(3)
                
                CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(1)) = CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(1)) + Cdd(1)
                CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(2)) = CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(2)) + Cdd(2)
                CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(3)) = CD_tmp(j,Mesh%Tfacette(k)%Tnoeud(3)) + Cdd(3)
                
                ! Traitement des points doubles pour �viter redondance dans le calcul des Coeff d'Infl.
                if(NDouble.ne.0.and.j.lt.TDouble(1).and..not.present(bornes))then
                    do r = 1,NDouble
                        
                        jtemp = TDouble(r)
                        CS(jtemp,Facette%Tnoeud(1)) = CS(jtemp,Facette%Tnoeud(1)) + Css(1)
                        CS(jtemp,Facette%Tnoeud(2)) = CS(jtemp,Facette%Tnoeud(2)) + Css(2)
                        CS(jtemp,Facette%Tnoeud(3)) = CS(jtemp,Facette%Tnoeud(3)) + Css(3)
                        
                        CD(jtemp,Facette%Tnoeud(1)) = CD(jtemp,Facette%Tnoeud(1)) + Cdd(1)
                        CD(jtemp,Facette%Tnoeud(2)) = CD(jtemp,Facette%Tnoeud(2)) + Cdd(2)
                        CD(jtemp,Facette%Tnoeud(3)) = CD(jtemp,Facette%Tnoeud(3)) + Cdd(3)
                        
                        !col1_CS(jtemp - borne(1,1) + 1) = col1_CS(jtemp - borne(1,1) + 1) + Css(1)
                        !col2_CS(jtemp - borne(1,1) + 1) = col2_CS(jtemp - borne(1,1) + 1) + Css(2)
                        !col3_CS(jtemp - borne(1,1) + 1) = col3_CS(jtemp - borne(1,1) + 1) + Css(3)
                        !
                        !col1_CD(jtemp - borne(1,1) + 1) = col1_CD(jtemp - borne(1,1) + 1) + Cdd(1)
                        !col2_CD(jtemp - borne(1,1) + 1) = col2_CD(jtemp - borne(1,1) + 1) + Cdd(2)
                        !col3_CD(jtemp - borne(1,1) + 1) = col3_CD(jtemp - borne(1,1) + 1) + Cdd(3)
                        
                    end do
                end if
                
            end if
            
        end do
        
        ! Solid angle.
        !if(not(present(bornes)))then
        !    
        !    ! Sum of the three colums to compute the partial solid angle.
        !    diag = diag - col1_CD - col2_CD - col3_CD ! - borne(1,1) + 1 not used because borne(1,1) = 1 in this case.
        !    
        !    ! The diagonal coefficient is substracts to ensure CD(j,j) = 0.
        !    diag(Facette%Tnoeud(1)) = diag(Facette%Tnoeud(1)) + Col1_CD(Facette%Tnoeud(1))
        !    diag(Facette%Tnoeud(2)) = diag(Facette%Tnoeud(2)) + Col2_CD(Facette%Tnoeud(2))
        !    diag(Facette%Tnoeud(3)) = diag(Facette%Tnoeud(3)) + Col3_CD(Facette%Tnoeud(3))
        !                
        !end if        
        
        !CD(borne(1,1):borne(1,2),Facette%Tnoeud(1)) = CD(borne(1,1):borne(1,2),Facette%Tnoeud(1)) + col1_CD
        !CD(borne(1,1):borne(1,2),Facette%Tnoeud(2)) = CD(borne(1,1):borne(1,2),Facette%Tnoeud(2)) + col2_CD
        !CD(borne(1,1):borne(1,2),Facette%Tnoeud(3)) = CD(borne(1,1):borne(1,2),Facette%Tnoeud(3)) + col3_CD
        !
        !CS(borne(1,1):borne(1,2),Facette%Tnoeud(1)) = CS(borne(1,1):borne(1,2),Facette%Tnoeud(1)) + col1_CS
        !CS(borne(1,1):borne(1,2),Facette%Tnoeud(2)) = CS(borne(1,1):borne(1,2),Facette%Tnoeud(2)) + col2_CS
        !CS(borne(1,1):borne(1,2),Facette%Tnoeud(3)) = CS(borne(1,1):borne(1,2),Facette%Tnoeud(3)) + col3_CS
        
        !deallocate(col1_CD,col2_CD,col3_CD)
        !deallocate(col1_CS,col2_CS,col3_CS)
                
    end do
            
    ! Calcul de l'Angle Solide.
    if(not(present(bornes)))then
        do j = 1,Mesh%Nnoeud
            CD(j,j) = 0._RP
            CD(j,j) = -sum(CD(j,1:Mesh%Nnoeud))
            !CD(j,j) = diag(j) ! - borne(1,1) + 1 not used because borne(1,1) = 1 in this case.
        end do
    end if
    
    !if(not(present(bornes))) deallocate(diag)
        
end subroutine CoeffInfl_v2

! ------------------------------------------------------------------------
! CoeffInfl_Line : Calcul des coefficients d'influence pour une ligne donn�e
! ------------------------------------------------------------------------
subroutine CoeffInfl_Line(Mesh, CS, CD, size_CDCS, ind, j0, bornes)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage),intent(in)                      :: Mesh                             ! Mesh.
    integer,intent(in)                              :: size_CDCS                        ! Number of nodes in the row.
    real(rp),dimension(size_CDCS),intent(inout)     :: CS, CD                           ! Influence coefficient VECTORS (not matrix!).
    integer,intent(in)                              :: ind                              ! Node from which the influence coefficients are computed for each panel of a part of a mesh (bornes).
    integer,intent(in)                              :: j0                               ! borne(1,1) : first node of the part of the mesh.
    integer, dimension(2,2), intent(in), optional   :: bornes                           ! Boundaries for the computation of the influence coefficients (nodes and panels).
    
    integer, dimension(2,2)                         :: borne                            ! Boundaries.
    integer                                         :: k, js, ns, js1, ns1, j1, j2, j3  ! Loop parameters.
    real(rp)                                        :: Cr, Crmax, invRho, ZG            ! Parameters.
    real(rp)                                        :: Ssigma, Smu                      ! Integrals.
    real(rp), dimension(3)                          :: Css, Cdd                         ! Influence coefficients
    real(rp), dimension(3)                          :: M0, M, MG                        ! Points.
    real(rp), dimension(3)                          :: Isigma, Imu, Ar                  ! Integrals
    !f2py integer*1, dimension(1000)                :: Facette
    type(TFacette)                                  :: Facette                          ! Panel.
    real(rp),parameter                              :: inv3 = 0.3333333333333333_rp     ! 1/3.
    
    ! This subroutine computes the influence coefficients for a node (influence of each panel of a part of the mesh).
    
    if(idebug>1) print*,"CoeffInl_Lin"
    
    ! Boundaries.
    if (present(bornes))then
        borne = bornes
    else
        borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    end if
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    end if
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    end if
    
    ! The following algorithm matches with the subroutine CoeffInfl but only for one node (influence of one node on some panels).
    
    ! Node from which the influence coefficients are computed.
    M0 = Mesh%Tnoeud(ind)%Pnoeud
    
    do k = borne(2,1),borne(2,2) ! Panels.
        
        M = M0
        Facette = Mesh%Tfacette(k)
        Crmax = 64._rp*Facette%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.
        Css = 0._rp ; Cdd = 0._rp
        
        ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
        do js1 = 1,-ns1,-2
            
            M(2) = js1*M(2)

            
            ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
            do js = 1,-ns,-2
                
                M(3) = js*M(3)-(1._rp-js)*Ldom(3)
                MG = M-Facette%Gfacette
                
                ! Distance criteria to use or not an symptotic development.
                Cr = dot_product(MG,MG)
                
                if(Cr.lt.Crmax)then ! Exact computation.
                    
                    ! Exact computation.
                    call Iints(M, Facette, Isigma, Imu)
                    call Sints(M, Facette, Ssigma, Smu)
                    
                    ! Prise en compte du param�trage de la facette.
                    Ar = [inv3,inv3,inv3] + matmul(Facette%dsT ,MG)
                    Css = Css + Ssigma*Ar + matmul(Facette%dsT ,Isigma)
                    Cdd = Cdd + Smu*Ar + matmul(Facette%dsT ,Imu)
                    
                else ! Asymptotic development.
                    
                    invRho=1._RP/norm2(MG) ! invRho = 1/||GM||
                    ZG = dot_product(MG, Facette%Normale) ! ZG = GM.n
                    
                    ! I_sigma (Eq 3.53 of LL).
                    Css = Css + Facette%Aire*invRho*inv3
                    
                    ! I_mu (Eq 3.54 of LL).
                    Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
                    
                end if
                
            end do
            
        end do
        
        ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).
        j1 = Facette%Tnoeud(1) - j0 + 1
        j2 = Facette%Tnoeud(2) - j0 + 1
        j3 = Facette%Tnoeud(3) - j0 + 1
        
        CS(j1) = CS(j1) + Css(1)
        CS(j2) = CS(j2) + Css(2)
        CS(j3) = CS(j3) + Css(3)
        
        CD(j1) = CD(j1) + Cdd(1)
        CD(j2) = CD(j2) + Cdd(2)
        CD(j3) = CD(j3) + Cdd(3)
        
    end do      

end subroutine CoeffInfl_Line

subroutine CoeffInfl_Line_Full(Mesh, Nnodes, CS, CD, LTab, Border,j0,bornes)
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage),intent(in)                          :: Mesh                             ! Mesh.
    integer, intent(in)                                 :: Nnodes                           ! Number of nodes in the linear system.
    real(rp), dimension(Nnodes,Nnodes), intent(inout)   :: CD,CS                            ! Influence coefficient matrices.
    logical, dimension(Nnodes), intent(in)              :: Ltab                             ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    integer, dimension(2,2), intent(in), optional       :: bornes                           ! Boundaries for the computation of the influence coefficients (nodes and panels).
    integer,intent(in)                                  :: j0                               ! borne(1,1) : first node of the part of the mesh.
    integer,dimension(1),intent(in)                     :: Border                           ! Boundaries of the modification of CD and CS.
    
    integer, dimension(2,2)                             :: borne                            ! Boundaries.
    integer                                             :: k, js, ns, js1, ns1,j1,j2,j3,jk  ! Loop parameters.
    real(rp)                                            :: Cr, Crmax, invRho, ZG            ! Parameters.
    real(rp)                                            :: Ssigma, Smu                      ! Integrals.
    real(rp), dimension(3)                              :: Css, Cdd                         ! Influence coefficients
    real(rp), dimension(3)                              :: M0, M, MG                        ! Points.
    real(rp), dimension(3)                              :: Isigma, Imu, Ar                  ! Integrals
    !f2py integer*1, dimension(1000)                    :: Facette
    type(TFacette)                                      :: Facette                          ! Panel.
    real(rp),parameter                                  :: inv3 = 0.3333333333333333_rp     ! 1/3.
    real(rp),dimension(:),allocatable                   :: CS_tmp, CD_tmp                   ! Influence coefficient VECTORS (not matrix!).
        
    ! This subroutine computes the influence coefficients for a node (influence of each panel of a part of the mesh).
    
    if(idebug>1) print*,"CoeffInl_Lin"
    
    ! Boundaries.
    if (present(bornes))then
        borne = bornes
    else
        borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    end if
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    end if
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    end if
    
    ! The following algorithm matches with the subroutine CoeffInfl but only for one node (influence of one node on some panels).
    
    ! Initialization of the lines of CD and CS.
    do jk = borne(1,1),borne(1,2) ! Nodes.
        if(not(LTab(jk)))then
            CS(jk,Border(1):Border(2)) = 0._rp
            CD(jk,Border(1):Border(2)) = 0._rp
        end if
    end do
    
    !$OMP PARALLEL DO NUM_THREADS(NThreads) SHARED(Mesh,borne,CD,CS,LTab,ns,ns1,Ldom,Border,j0) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC,(borne(1,2)-borne(1,1)+1)/Nthreads)
    do jk = borne(1,1),borne(1,2) ! Nodes.
                
        if(not(LTab(jk)))then ! If the node is taken into account (LTab = false).
            
            allocate(CS_tmp(Border(2)-Border(1)+1),CD_tmp(Border(2)-Border(1)+1))
            CS_tmp = 0._RP
            CD_tmp = 0._RP
            
            ! Node from which the influence coefficients are computed.
            M0 = Mesh%Tnoeud(jk)%Pnoeud
            
            do k = borne(2,1),borne(2,2) ! Panels.
        
                M = M0
                Facette = Mesh%Tfacette(k)
                Crmax = 64._rp*Facette%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.
                Css = 0._rp ; Cdd = 0._rp
                
        
                ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
                do js1 = 1,-ns1,-2
            
                    M(2) = js1*M(2)
                    
                
            
                    ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
                    do js = 1,-ns,-2
                
                        M(3) = js*M(3)-(1._rp-js)*Ldom(3)
                        MG = M-Facette%Gfacette
                        
                        ! Distance criteria to use or not an symptotic development.
                        Cr = dot_product(MG,MG)
                        
                        if(Cr.lt.Crmax)then ! Exact computation.
                            
                            ! Exact computation.
                            call Iints(M, Facette, Isigma, Imu)
                            call Sints(M, Facette, Ssigma, Smu)
                            
                            ! Prise en compte du param�trage de la facette.
                            Ar = [inv3,inv3,inv3] + matmul(Facette%dsT ,MG)
                            Css = Css + Ssigma*Ar + matmul(Facette%dsT ,Isigma)
                            Cdd = Cdd + Smu*Ar + matmul(Facette%dsT ,Imu)
                            
                        else ! Asymptotic development.
                            
                            invRho = 1._RP/norm2(MG) ! invRho = 1/||GM||
                            ZG = dot_product(MG, Facette%Normale) ! ZG = GM.n
                            
                            ! I_sigma (Eq 3.53 of LL).
                            Css = Css + Facette%Aire*invRho*inv3
                            
                            ! I_mu (Eq 3.54 of LL).
                            Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
                            
                        end if
                        
                    end do
                    
                end do
                
                ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).                
                j1 = Facette%Tnoeud(1) - j0 + 1
                j2 = Facette%Tnoeud(2) - j0 + 1
                j3 = Facette%Tnoeud(3) - j0 + 1
                
                CS_tmp(j1) = CS_tmp(j1) + Css(1)
                CS_tmp(j2) = CS_tmp(j2) + Css(2)
                CS_tmp(j3) = CS_tmp(j3) + Css(3)
                
                CD_tmp(j1) = CD_tmp(j1) + Cdd(1)
                CD_tmp(j2) = CD_tmp(j2) + Cdd(2)
                CD_tmp(j3) = CD_tmp(j3) + Cdd(3)
                
            end do
            
            CS(jk,Border(1):Border(2)) = CS(jk,Border(1):Border(2)) + CS_tmp(:)
            CD(jk,Border(1):Border(2)) = CD(jk,Border(1):Border(2)) + CD_tmp(:)
            
            deallocate(CS_tmp,CD_tmp)
            
        end if
        
    end do
    !$OMP END PARALLEL DO
    
end subroutine CoeffInfl_Line_Full

! -----------------------------------------------------------------------
! CoeffInfl_Col : Calcul des coefficients d'influence pour des colonnes donn�es
! -----------------------------------------------------------------------
subroutine CoeffInfl_Col(Mesh, CD1, CS1, BPoint1, CD2, CS2, BPoint2, CD3, CS3, BPoint3, size_CDCS,Nnodes,ind, Ltab, bornes)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage),intent(in)                      :: Mesh                                 ! Mesh.
    integer,intent(in)                              :: size_CDCS                            ! Number of nodes in the column.
    real(rp), dimension(size_CDCS), intent(inout)   :: CD1, CD2, CD3, CS1, CS2, CS3         ! Influence coefficient vectors (three columns for three nodes per panel).
    integer, intent(in)                             :: ind                                  ! Panel from which the influence coefficients are computed for each node of a part of a mesh (bornes).
    integer, intent(in)                             :: Nnodes                               ! Number of nodes in the linear system.
    logical, dimension(Nnodes), intent(in)          :: Ltab                                 ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    integer, dimension(2,2), optional               :: bornes                               ! Boundaries for the computation of the influence coefficients (nodes and panels).
    logical, intent(in)                             :: BPoint1, BPoint2, BPoint3            ! Neighbours (nodes) of the nodes of LTab.
    
    integer, dimension(2,2)                         :: borne                                ! Boundaries.
    integer                                         :: j, js, ns, js1, ns1, j1              ! Loop parameters.
    real(rp)                                        :: Crmax, Cr, Ssigma, Smu, invRho, ZG   ! Parameters and integrals.
    real(rp), dimension(3)                          :: Css, Cdd, M, MG, Ar, Isigma, Imu     ! Influence coefficients and integrals.
    !f2py integer*1, dimension(1000)                :: Facette
    type(TFacette)                                  :: Facette                              ! Panel.
    real(rp),parameter                              :: inv3 = 0.3333333333333333_rp         ! 1/3.
    
    ! This subroutine computes the influence coefficients for a panel (influence of each node of a part of the mesh).
    
    ! Boundaries.
    if(present(bornes))then
        borne = bornes
    else
        borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    end if
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    end if
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    end if
        
    ! Panel from which the influence coefficients are computed.
    Facette = Mesh%Tfacette(ind)
    Crmax = 64._rp*Facette%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.

    do j = borne(1,1),borne(1,2) ! Nodes.
        
        M = Mesh%Tnoeud(j)%Pnoeud
        
        ! The node requieres to compute its influence coefficients again.
        if(Ltab(j))then
            
            Css = 0._rp ; Cdd = 0._rp
            
            ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
            do js1 = 1,-ns1,-2
                
                M(2) = js1*M(2)
                
                ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
                do js = 1,-ns,-2
                    
                    M(3) = js*M(3)-(1._rp-js)*Ldom(3)
                    MG = M-Facette%Gfacette
                    
                    ! Distance criteria to use or not an symptotic development.
                    Cr = dot_product(MG,MG)

                    if(Cr.lt.Crmax)then ! Exact computation.
                        
                        ! Exact computation.
                        call Iints(M, Facette, Isigma, Imu)
                        call Sints(M, Facette, Ssigma, Smu)
                        
                        ! Prise en compte du param�trage de la facette.
                        Ar = [inv3,inv3,inv3] + matmul(Facette%dsT ,MG)
                        Css = Css + Ssigma*Ar + matmul(Facette%dsT ,Isigma)
                        Cdd = Cdd + Smu*Ar + matmul(Facette%dsT ,Imu)
                        
                    else! Asymptotic development.
                        
                        invRho = 1._RP/norm2(MG) ! invRho = 1/||GM||
                        ZG = dot_product(MG, Facette%Normale)
                        
                        ! I_sigma (Eq 3.53 of LL).
                        Css = Css + Facette%Aire*invRho*inv3
                        
                        ! I_mu (Eq 3.54 of LL).
                        Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
                        
                    end if
                    
                end do
                
            end do
            
            ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).
            j1 = j - borne(1,1) + 1
            
            if(BPoint1)then
                CS1(j1) = CS1(j1) + Css(1)
                CD1(j1) = CD1(j1) + Cdd(1)
            end if
            
            if(BPoint2)then
                CS2(j1) = CS2(j1) + Css(2)
                CD2(j1) = CD2(j1) + Cdd(2)
            end if
            
            if(BPoint3)then
                CS3(j1) = CS3(j1) + Css(3)        
                CD3(j1) = CD3(j1) + Cdd(3)
            end if
            
        end if
        
    end do

end subroutine CoeffInfl_Col

subroutine CoeffInfl_Col_Full(Mesh, CD, CS, Nnodes, Nfacettes,LTab,BPoint,BFace,bornes)    
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage), intent(in)                         :: Mesh                                 ! Mesh.
    integer,intent(in)                                  :: Nnodes                               ! Number of nodes in the mesh.
    integer,intent(in)                                  :: Nfacettes                            ! Number of panels in the mesh.
    real(rp), dimension(Nnodes,Nnodes), intent(inout)   :: CD,CS                                ! Influence coefficient matrices.
    integer, dimension(2,2), intent(in), optional       :: bornes                               ! Boundaries for the computation of the influence coefficients (nodes = rows and panels != columns !!!).
    logical, dimension(Nnodes), intent(in)              :: Ltab                                 ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    logical, dimension(Nnodes)                          :: BPoint                               ! Neighbours (nodes) of the nodes of LTab.
    logical,dimension(Nfacettes)                        :: BFace                                ! Panels linked to the nodes of LTab and their neighbours (panels).
    
    integer, dimension(2,2)                             :: borne                                ! Boundaries.
    integer                                             :: j, js, ns, js1, ns1, j1,jk           ! Loop parameters.
    real(rp)                                            :: Crmax, Cr, Ssigma, Smu, invRho, ZG   ! Parameters and integrals.
    real(rp), dimension(3)                              :: Css, Cdd, M, MG, Ar, Isigma, Imu     ! Influence coefficients and integrals.
    !f2py integer*1, dimension(1000)                    :: Facette
    type(TFacette)                                      :: Facette                              ! Panel.
    real(rp),parameter                                  :: inv3 = 0.3333333333333333_rp         ! 1/3.
    real(rp), dimension(:),allocatable                  :: CD1, CD2, CD3, CS1, CS2, CS3         ! Influence coefficient vectors (three columns for three nodes per panel).
    logical                                             :: BPoint1, BPoint2, BPoint3            ! Neighbours (nodes) of the nodes of LTab.
    integer                                             :: j11,j22,j33                          ! Indexes of the nodes of a panel.
    
    ! This subroutine computes the influence coefficients for a panel (influence of each node of a part of the mesh).
    
    ! Boundaries.
    if(present(bornes))then
        borne = bornes
    else
        borne = reshape((/1,1,Mesh%Nnoeud,Mesh%Nfacette/), (/2,2/))
    end if
    
    ! Bottom symmetry.
    if(bottom_sym)then
        ns = 1
    else
        ns = 0
    end if
    
    ! Symmetry (xOz).
    if(Symmetry)then
        ns1 = 1
    else
        ns1 = 0
    end if
    
    !$OMP PARALLEL DO NUM_THREADS(NThreads) SHARED(Mesh,borne,CD,CS,BFace,Ltab,BPoint,Ldom,ns,ns1) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC,(borne(2,2)-borne(2,1)+1)/Nthreads)
    do jk = borne(2,1),borne(2,2)
        
        if(BFace(jk))then
            
            allocate(CD1(borne(1,2)-borne(1,1)+1),CD2(borne(1,2)-borne(1,1)+1),CD3(borne(1,2)-borne(1,1)+1))
            allocate(CS1(borne(1,2)-borne(1,1)+1),CS2(borne(1,2)-borne(1,1)+1),CS3(borne(1,2)-borne(1,1)+1))
            
            CD1 = 0._RP
            CD2 = 0._RP
            CD3 = 0._RP
            
            CS1 = 0._RP
            CS2 = 0._RP
            CS3 = 0._RP
            
            j11 =  Mesh%Tfacette(jk)%Tnoeud(1)
            j22 =  Mesh%Tfacette(jk)%Tnoeud(2)
            j33 =  Mesh%Tfacette(jk)%Tnoeud(3)
            
            BPoint1 = BPoint(j11)
            BPoint2 = BPoint(j22)
            BPoint3 = BPoint(j33)
                
            ! Panel from which the influence coefficients are computed.
            Facette = Mesh%Tfacette(jk)
            Crmax = 64._rp*Facette%Rmax**2 ! Crmax = (8*Rmax)^2, the value of 8 is justified in 3.4 of LL.
            
            do j = borne(1,1),borne(1,2) ! Nodes.
                
                M = Mesh%Tnoeud(j)%Pnoeud
                
                ! The node requieres to compute its influence coefficients again.
                if(Ltab(j))then
                    
                    Css = 0._rp ; Cdd = 0._rp
                    
                    ! If (xOz) symmetry, computation for M(x,y,z) and M(x,-y,z).
                    do js1 = 1,-ns1,-2
                        
                        M(2) = js1*M(2)
                        
                        ! If bottom symmetry, computation for M(x,y,z) and M(x,y,-z-2h).
                        do js = 1,-ns,-2
                            
                            M(3) = js*M(3)-(1._rp-js)*Ldom(3)
                            MG = M-Facette%Gfacette
                            
                            ! Distance criteria to use or not an symptotic development.
                            Cr = dot_product(MG,MG)

                            if(Cr.lt.Crmax)then ! Exact computation.
                                
                                ! Exact computation.
                                call Iints(M, Facette, Isigma, Imu)
                                call Sints(M, Facette, Ssigma, Smu)
                                
                                ! Prise en compte du param�trage de la facette.
                                Ar = [inv3,inv3,inv3] + matmul(Facette%dsT,MG)
                                Css = Css + Ssigma*Ar + matmul(Facette%dsT,Isigma)
                                Cdd = Cdd + Smu*Ar + matmul(Facette%dsT,Imu)
                                
                            else! Asymptotic development.
                                
                                invRho = 1._RP/norm2(MG) ! invRho = 1/||GM||
                                ZG = dot_product(MG, Facette%Normale)
                                
                                ! I_sigma (Eq 3.53 of LL).
                                Css = Css + Facette%Aire*invRho*inv3
                                
                                ! I_mu (Eq 3.54 of LL).
                                Cdd = Cdd + Facette%Aire*ZG*invRho/Cr*inv3
                                
                            end if
                            
                        end do
                        
                    end do
                    
                    ! Redistribution par colocation (computations are done for a panel but CD/CS but are wanted for the 3 nodes of the panel).
                    j1 = j - borne(1,1) + 1
                    
                    if(BPoint1)then
                        CS1(j1) = CS1(j1) + Css(1)
                        CD1(j1) = CD1(j1) + Cdd(1)
                    end if
                    
                    if(BPoint2)then
                        CS2(j1) = CS2(j1) + Css(2)
                        CD2(j1) = CD2(j1) + Cdd(2)
                    end if
                    
                    if(BPoint3)then
                        CS3(j1) = CS3(j1) + Css(3)        
                        CD3(j1) = CD3(j1) + Cdd(3)
                    end if
                    
                end if
                
            end do
            
            CD(borne(1,1):borne(1,2),j11) = CD(borne(1,1):borne(1,2),j11) + CD1(1:borne(1,2)-borne(1,1)+1)
            CS(borne(1,1):borne(1,2),j11) = CS(borne(1,1):borne(1,2),j11) + CS1(1:borne(1,2)-borne(1,1)+1)
            
            CD(borne(1,1):borne(1,2),j22) = CD(borne(1,1):borne(1,2),j22) + CD2(1:borne(1,2)-borne(1,1)+1)
            CS(borne(1,1):borne(1,2),j22) = CS(borne(1,1):borne(1,2),j22) + CS2(1:borne(1,2)-borne(1,1)+1)
            
            CD(borne(1,1):borne(1,2),j33) = CD(borne(1,1):borne(1,2),j33) + CD3(1:borne(1,2)-borne(1,1)+1)
            CS(borne(1,1):borne(1,2),j33) = CS(borne(1,1):borne(1,2),j33) + CS3(1:borne(1,2)-borne(1,1)+1)
            
            deallocate(CD1,CD2,CD3,CS1,CS2,CS3)
            
        end if
        
    end do
    !$OMP END PARALLEL DO
    
end subroutine CoeffInfl_Col_Full

subroutine ZoneInfl(Mesh, LTab, borne, BPoint, BFace,Nnodes,Nfacettes)

    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage),intent(in)              :: Mesh                 ! Mesh.
    integer,intent(in)                      :: Nnodes,Nfacettes     ! Number of nodes and panel.
    logical, dimension(Nnodes), intent(in)  :: LTab                 ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    integer, dimension(2,2), intent(in)     :: borne                ! Boundaries of LTab.
    logical, dimension(Nnodes)              :: BPoint               ! Neighbours (nodes) of the nodes of LTab.
    logical,dimension(Nfacettes)            :: BFace                ! Panels linked to the nodes of LTab and their neighbours (panels).
    
    integer                                 :: j, k ,jk, nf, np     ! Loop parameters.
    
    ! This subroutine fills the tables BPoint and BFace.
    
    do j = borne(1,1),borne(1,2) ! All the nodes.
        if(.not.LTab(j))then ! Node taken into account.
            
            nf = Mesh%Tnoeud(j)%Nfacette ! Number of panels linked to the node j.
            np = Mesh%Tnoeud(j)%NVoisin(2) ! Number of neigbours (nodes) of the node j.
            
            ! All the panels linked to the node are taken into account.
            do k = 1,nf
                jk = Mesh%Tnoeud(j)%TFacette(k,1)
                BFace(jk) = .true.
            end do
            
            ! All the neighboors of the node are taken into account.
            k = 1
            do while(Mesh%Tnoeud(j)%TVoisin(k,2).le.1 .and. k.le.np) ! TVoisin(k,2).le.1 : ordre de voisinage <= 1.
                jk = Mesh%Tnoeud(j)%TVoisin(k,1)
                BPoint(abs(jk)) = .true.
                k = k + 1
            end do
            
        end if
    end do
    
    ! The panels where the nodes taken into account are, are taken into account as well.
    do j = borne(1,1),borne(1,2) ! All the nodes.
        if(BPoint(j))then
            nf = Mesh%Tnoeud(j)%Nfacette
            do k = 1,nf
                jk = Mesh%Tnoeud(j)%TFacette(k,1)
                BFace(jk) = .true.
            end do
        end if
    end do
    
end subroutine ZoneInfl

subroutine DeplNoeud(Mesh, LTab, borne,Nnodes)

    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage),intent(in)                  :: Mesh     ! Mesh.
    integer,intent(in)                          :: Nnodes   ! Number of nodes in the mesh.
    logical, dimension(Nnodes), intent(inout)   :: LTab     ! Table to know if a node moves enough to compute its influence coefficients again (= false) or not (= true).
    integer, dimension(2,2), intent(in)         :: borne    ! Boundaries of LTab.

    integer                                     :: j, k, jk ! Loop parameters.

    ! This subroutine fills the table LTab, usefull to know which influence coefficients have to be computed again or not.
    
    ! LTab(borne(1,1):borne(1,2)) = .true. ! Useless, already done in solBVP.
    
    do j = borne(1,1),borne(1,2)
        
        ! If the node moves enought, it is taken into account in the computation of the influence computation.
        if(norm2(Mesh%Tnoeud(j)%Velocity).gt.Epsilon) Ltab(j) = .false.
        
        ! Case of the twin nodes.
        do k = 1,Mesh%Tnoeud(j)%Ndouble
            jk = Mesh%Tnoeud(j)%double(k)
            LTab(j) = LTab(j) .and. Mesh%Tnoeud(jk)%typeNoeud .ne. 0
        enddo    
    enddo

end subroutine DeplNoeud



subroutine Iints(M, Facette, Isigma, Imu)
    !!!!! Probl�me :
    !   Calculer Isigma et Imu d'apr�s les donn�es g�om�triques d'une noeud et d'une facette
    !   Solution analytique d�velopp�e dans le rapport bibliographique
    !   Int�grales de contour d�compos�es sur chaque segment de la facette, puis somm�es vectoriellement
    !!!!!
    
    real(rp),dimension(3),intent(in)    :: M                                                        ! Point.
    !f2py integer*1, dimension(1000)    :: Facette
    type(TFacette),intent(in)           :: Facette                                                  ! Panel.
    real(rp), dimension(3),intent(out)  :: Isigma, Imu                                              ! Isigma and Imu.
    
    integer                             :: j
    real(rp)                            :: K1, q1, q2, AM, BM, AB, sgn, d1, d2, ABM, invK1, invAB   !
    real(rp), dimension(3)              :: dsT,vect_product                                         !
    real(rp), dimension(3)              :: A, B                                                     ! End points of one segment.
    real(rp), dimension(3,4)            :: Pt                                                       ! Positions of the vertexes of the panel.
    real(rp), dimension(3,3)            :: Jsigma, Jmu                                              ! Integral Isigma and Imu.
    
    integer :: i_n, i_f
    
    ! This subroutine computes Isigma and Imu.
    
    Pt(:,1:3) = Facette%Pnoeud
    Pt(:,4) = Pt(:,1) ! permet de revenir au point de d�part
    
    ! Line integration over the three segments of the panel.
    do j=1,3
        
        A = Pt(:,j)        
        B = Pt(:,j+1)
        AM = norm2(M-A)
        AB = norm2(B-A)
        BM = norm2(M-B)
        call Computation_vect_product(Facette%Normale,B-A,dsT)
        
        if (abs(AM).lt.Epsilon.or.abs(BM).lt.Epsilon)then ! M = A ou M = B
            Jsigma(:,j) = -0.5_RP*AB*dsT
            Jmu(:,j) = 0._RP
        else

            
            
            
            if (abs(abs(dot_product(B-A,M-A))-(AB*AM)) .lt. Epsilon)then !((K1<1.0E-8_RP)) then ! M est sur la droite AB en dehors du segment [AB]
                sgn = -dot_product(B-A,M-A)/AM
                Jsigma(:,j) = -abs(AM+0.5*sgn)*dsT
                Jmu(:,j) = 0._RP
            else ! Cas R�gulier
                invAB = 1._RP/AB
                ABM = dot_product(B-A,M-A)*invAB
                K1 = AM**2 - (ABM)**2
                invK1 = 1._RP/sqrt(K1)
                q1 = -ABM*invK1
                q2 = (AB - ABM)*invK1
                d1 = asinh(q1)
                d2 = asinh(q2)
                
     
                
                
                
                Jsigma(:,j) = -0.5_RP*K1*( d2 - d1 + q2*sqrt(1._RP+q2*q2) - q1*sqrt(1._RP+q1*q1) )*dsT*invAB
                call Computation_vect_product(A-M,B-A,vect_product)
                dsT = vect_product*invAB
                Jmu(:,j) = dsT*(d1-d2)
                
                
            end if
        end if
    end do
    
    Isigma = sum(Jsigma,2)
    Imu = sum(Jmu,2)
    
end subroutine  Iints
    
subroutine Sints(M, Facette, Ssigma, Smu)
    !!!!! Probl�me :
    !   Calcul des int�grales de surface Ssigma et Smu par les expressions de Guevel
    !   Traitement du cas singulier, M est un des sommets de la facette
    !
    !!!!!
    
    real(rp),dimension(3),intent(in)    :: M            ! Point.
    !f2py integer*1, dimension(1000)    :: Facette
    type(TFacette)                      :: Facette      ! Panel.
    real(rp)                            :: Ssigma, Smu  ! Ssigma and Smu.

    integer                             :: j, t         ! Loop parameters.
    real(rp)                            :: a, b         !
    !f2py integer*1, dimension(1000)    :: VC
    type(TVariableChange)               :: VC           !
    real(rp), dimension(3,3)            :: Pt           !
    real(rp), dimension(3,6)            :: Pinit        !
    
    ! Rep�rage point de contr�le sur la facette.
    t = 0
    do j = 1,3
        if(dot_product(M-Facette%Pnoeud(:,j),M-Facette%Pnoeud(:,j)).lt.1.0E-05)then
            t = j ! Control point not on the panel.
        else
            t = t ! t = 0, control point on the panel.
        end if
    end do
    
    if (t.ne.0) then !Point de contr�le sur la facette, t=n� du sommet
        
        !  Control point on the panel, that means at a vertex of the panel.
        
        Pinit(1:3,1:3) = Facette%Pnoeud
        Pinit(1:3,4:6) = Facette%Pnoeud
        Pt = Pinit(:,t:t+2)
        call VarChange(Pt, VC) ! Changement de variables locales
        a = (norm2(VC%Puvw(:,2))**2)/sqrt(VC%Delta)
        b = (dot_product(VC%Puvw(:,1),VC%Puvw(:,2)))/sqrt(VC%Delta)
        Ssigma = (VC%Jacobien)/norm2(VC%Puvw(:,2))*(asinh(a+b)-asinh(-a+b)) ! Eq 3.28 of LL.
        Smu = 0 ! No singurality with Smu because it is the solid angle.
        
    else
        
        ! Control point not on the panel.
        
        call Guevel(M,Facette, Ssigma, Smu) ! Guevel's expression for Ssigma but nor for Smu even if both of them are computed in the subroutine Guevel.
        
    endif
    
end subroutine Sints
    
subroutine Guevel(M, Facette, Ssigma, Smu)
    !!!!! Probl�me :
    !   Calcul de Ssigma et Smu par la solution exacte analytique donn�e par Guevel.
    !   Les singularit�s ont �t� trait�es s�par�ment dans la subroutine Sints
    !   Il ne reste donc que le cas o� M n'est pas sur la facette.
    !!!!!
    
    real(rp), dimension(3)              :: M                                ! Point.
    !f2py integer*1, dimension(1000)    :: Facette
    type(TFacette)                      :: Facette                          ! Panel.
    real(rp)                            :: Ssigma, Smu                      ! Ssigma and Smu when the point is not on the panel.
    
    integer                             :: j                                ! Loop parameter.
    real(rp)                            :: Z, Z1                            !
    real(rp), dimension(3)              :: J1,J2,d,Y,D1,Dt,N1,vect_product  !
    real(rp), dimension(4)              :: R                                !
    real(rp), dimension(3,3)            :: P,test                           !
    
    Z = dot_product(M-Facette%Gfacette,Facette%Normale)
    R(1) = norm2(M-Facette%Pnoeud(:,1))
    
    ! Loop over the 3 nodes of the panel.
    do j = 1,3
        
        if(j==3)then
            P(j,:) = Facette%Pnoeud(:,1)-Facette%Pnoeud(:,3)
            R(j+1) = norm2(M-Facette%Pnoeud(:,1))
        else
            P(j,:) = Facette%Pnoeud(:,j+1)-Facette%Pnoeud(:,j)
            R(j+1) = norm2(M-Facette%Pnoeud(:,j+1))
        endif
        
        d(j) = norm2(P(j,:))
        call Computation_vect_product(Facette%Normale,P(j,:),test(j,:))
        call Computation_vect_product(Facette%Normale,P(j,:),vect_product)
        Y(j) = dot_product(M-Facette%PNoeud(:,j),vect_product)/d(j)
        N1(j) = R(j+1)+R(j)+d(j)
        D1(j) = R(j+1)+R(j)-d(j)
        if (D1(j).lt.Epsilon) print*, 'Erreur Guevel : facette de coordonn�es : ', Facette%Tnoeud, Facette%Pnoeud, 'Noeud', M
        Dt(j) = (R(j+1)+R(j))**2-d(j)**2+2*abs(Z)*(R(j+1)+R(j))
        if (D1(j).lt.Epsilon) then
            print*,M
            print*,Facette%Gfacette
            print*,Facette%Normale
        end if
        
        J1(j) = Y(j)*log(N1(j)/D1(j))
        J2(j) = 2*atan(2*Y(j)*d(j)/Dt(j))
        
    end do
    
    Ssigma = sum(J1-abs(Z)*J2) ! Eq 3.18 of LL.
    Z1 = 1
    Smu = sum(J2)*sign(Z1,Z) ! Eq 3.43 of LL.
    
end subroutine Guevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       Fonction VarChange                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine VarChange(PNoeud, VariableChange)
!!!!! Probl�me :
!   Effectue le changement de variable correspondant au param�trage de la facette d�fini dans le rapport Biblio
!   Renvoie le tableau VarChange constitu� de
!       o la matrice de changement de base
!       o le jacobien de la transformation
!       o delta= norm2(i)�.norm2(j)�-<i,j>�
!!!!!
! Param�tres

real(rp), dimension(3,3) :: PNoeud
!f2py integer*1, dimension(1000) :: VariableChange
type(TVariableChange) :: VariableChange
! Variables Locales
real(rp), dimension(3,3) :: Puvw
! D�but            
Puvw(:,1) = 0.5*(PNoeud(:,2)-2*PNoeud(:,1)+PNoeud(:,3))
Puvw(:,2) = 0.5*(PNoeud(:,3)-PNoeud(:,2))
call Computation_vect_product(Puvw(:,1), Puvw(:,2),Puvw(:,3))
VariableChange%Puvw = Puvw
VariableChange%Jacobien = norm2(Puvw(:,3))
VariableChange%Delta = ((norm2(Puvw(:,1)))**2)*((norm2(Puvw(:,2)))**2) - dot_product(Puvw(:,1),Puvw(:,2))**2
! Fin
end subroutine VarChange

subroutine Angle_solid(CD,Nnoeud)
    
    real(rp),dimension(Nnoeud,Nnoeud),intent(inout) :: CD       ! Influence coefficient matrix.
    integer,intent(in)                              :: Nnoeud   ! Number of nodes.
    
    integer                                         :: j        ! Loop parameter.
    
    ! This subroutine computes the solid angles of each node.
    
    !$OMP PARALLEL DO NUM_THREADS(NThreads) SHARED(Nnoeud,CD) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC,Nnoeud/Nthreads)
    do j = 1,Nnoeud
        CD(j,j) = 0._rp
        CD(j,j) = -sum(CD(j,1:Nnoeud))
    end do
    !$OMP END PARALLEL DO

end subroutine Angle_solid

end module CoeffInfluence
