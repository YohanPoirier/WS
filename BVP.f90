module BVP
use FonctionsCommunes
use CoeffInfluence

implicit none

contains

!
!
! Dans cette subroutine on calcule les matrices A des deux systemes et onles inverse
subroutine initialisation_lineaire(Ecoulement, Mesh, Nnodes, CD, CS)

    !f2py integer*1, dimension(1000)        :: Ecoulement
    type(TEcoulement)                       :: Ecoulement               ! Flow parameters.
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage)                         :: Mesh                     ! Mesh.

    integer                                 :: Nnodes
    real(rp), dimension(Nnodes,Nnodes)      :: CD, CS                   ! Influence coefficient matrices.
  
    logical                                 :: Option                   ! present(Option) = false for the BVP on Phi, present(Option) = true and Option = true for the BVP on DPhiDt.
  
    real(rp),allocatable, dimension(:,:)    :: A                        ! A matrix of the linear system.

    real(rp) :: det
    integer :: N,i,j
    

    
    ! Allocation.
    N = Mesh%Nsys
    allocate(A(N,N))
    allocate(Ainv1(N,N)) 
    allocate(cond1(N)) 
    
    !Calculation of the coeffiient of influence
    CD = 0._RP ; CS = 0._RP
    call CoeffInfl(Mesh, CD, CS,N)
    
    ! Building of A, B and initialization of Sol.
    option = .false. !BVP sur phi
    call Calcul_A(CD, CS, Ecoulement, Mesh, A, cond1, Nnodes, N)
    
    ! Inversion of A
    call inv_matrice(A,Ainv1, N)

    ! Deallocating.
    deallocate(A)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                   !
!                                   Solveur du Problème aux Limites                                 !
!                        (Yohan : les coefficients CD et CS sont deja calcules)                     !
!                                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solBVP_cuda(Ecoulement, Mesh, CD, CS,Nnodes,time,boolRemesh, t, Option)
    !!!!!Problème :
    !   Résoudre le système des équations intégrales par
    !       o Calcul ou récupération des coefficients d'influence
    !       o Définition du système linéaire
    !       o Résolution du système linéaire
    !       o Redistribution des solutions
    !!!!!

    !f2py integer*1, dimension(1000)        :: Ecoulement
    type(TEcoulement)                       :: Ecoulement               ! Flow parameters.
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage), intent(in)             :: Mesh                     ! Mesh.
    integer,intent(in)                      :: Nnodes                   ! Number of nodes in the mesh.
    real(rp), dimension(Nnodes,Nnodes)      :: CD, CS                   ! Influence coefficient matrices.
    logical,intent(in),optional             :: boolRemesh               ! Boolean to know if the mesh of the bodies or/and the free surface was regenerated (true) or not (false).
    
    real(rp),dimension(20) :: time
    
    real(rp), intent(in), optional          :: t                        ! Parameter to know if the influence coefficients have to be computed or not (t = -1 : yes, t > 0 : no).
    logical, intent(in), optional           :: Option                   ! present(Option) = false for the BVP on Phi, present(Option) = true and Option = true for the BVP on DPhiDt.
    
    logical                                 :: Opt                      ! = false : BVP on Phi, = true : BVP on DPhiDt.
    integer                                 :: j,k, N, jk, j1, j2, j3   ! Parameters.
    integer, dimension(2,2)                 :: borne                    ! Boundaries for the computation of the influence coefficients (nodes and panels).
    real(rp)                                :: CIPartiel                ! Parameter to know if the partial computation of the influence coefficient has to be done or not.
    real(rp),allocatable, dimension(:,:)    :: A                        ! A matrix of the linear system.
    real(rp),allocatable, dimension(:)      :: B, Sol                   ! B matrix and the solution of the linear system.
    logical, dimension(:),allocatable       :: LTab                     ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    logical, dimension(:),allocatable       :: BPoint                   ! Neighbours (nodes) of the nodes of LTab.
    logical, dimension(:),allocatable       :: BFace                    ! Panels linked to the nodes of LTab and their neighbours (panels).
    integer                                 :: ierror                   ! Error flag.
    logical                                 :: bool                     ! Computation of the CI (true) or not (false) when it was already done after the remeshing (only for jk = 1).
    integer,dimension(2)                    :: Border                   ! Boundaries of the modification of CD and CS in case of calling CoeffInfl_Line_Full.
        
    ! This subroutine computes the influence coefficients, creates the linear system and solves it.

    ! CD and CS have already a value, consequently, in case of partial computation of the influence coefficients, only coefficients will be updates. The other ones will keep their value.
    
    ierror = 0
    
    ! Allocation.
    N = Mesh%Nsys
    allocate(A(N,N), B(N), Sol(N))

    ! Initialization of Opt.
    ! Opt = false for the BVP on Phi and true for the BVP on DPhiDt (forced motion).
    Opt = .false.
    if(present(Option))then
        ! Case of the second call of solBVP for a forced mesh.
        if(Option) Opt = .true.
    end if
        
    if(not(Opt))then

        call cpu_time(time(9))

    end if
    
    ! Building of A, B and initialization of Sol.
    call SystLin(CD, CS, Ecoulement, Mesh, A, B, Sol,Nnodes,N, Option)
    
    if(not(Opt))then

        call cpu_time(time(10))

    end if
    
    ! Solving of the linear system.
    if(iprint>0) print*,"Resolution système linéaire ..."
    if (Solv==0) then
        call GMRES(A, B, Sol, N, ierror)
    else  
        call LU(A, B, Sol, N)
    end if
    
    if(not(Opt))then

        call cpu_time(time(11))

    end if
    
    if(ierror/=0) goto 9999
    
    ! Distribution of the solution.
    if(iprint>0) print*,"Mise à jour champ ecoulement ..." 
    call postSL(Sol,N, Mesh, Ecoulement, Option)
    
    9999 continue
        if(ierror/=0)then
            write(*,90),ierror
        end if
    90 format('error #',i3,' : solBVP: pb. in solving the Boundary Value Problem.')

    ! Deallocating.
    deallocate(A, B, Sol)

end subroutine solBVP_cuda






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                   !
!                                   Solveur du Problème aux Limites                                 !
!                                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solBVP(Ecoulement, Mesh, CD, CS, Nnodes,time,boolRemesh, t, Option)
    !!!!!Problème :
    !   Résoudre le système des équations intégrales par
    !       o Calcul ou récupération des coefficients d'influence
    !       o Définition du système linéaire
    !       o Résolution du système linéaire
    !       o Redistribution des solutions
    !!!!!

    !f2py integer*1, dimension(1000)        :: Ecoulement
    type(TEcoulement)                       :: Ecoulement               ! Flow parameters.
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage), intent(in)             :: Mesh                     ! Mesh.
    integer,intent(in)                      :: Nnodes                   ! Number of nodes in the mesh.
    real(rp), dimension(Nnodes,Nnodes)      :: CD, CS                   ! Influence coefficient matrices.
    logical,intent(in),optional             :: boolRemesh               ! Boolean to know if the mesh of the bodies or/and the free surface was regenerated (true) or not (false).
    
    real(rp),dimension(20) :: time
    
    real(rp), intent(in), optional          :: t                        ! Parameter to know if the influence coefficients have to be computed or not (t = -1 : yes, t > 0 : no).
    logical, intent(in), optional           :: Option                   ! present(Option) = false for the BVP on Phi, present(Option) = true and Option = true for the BVP on DPhiDt.
    
    logical                                 :: Opt                      ! = false : BVP on Phi, = true : BVP on DPhiDt.
    integer                                 :: j,k, N, jk, j1, j2, j3   ! Parameters.
    integer, dimension(2,2)                 :: borne                    ! Boundaries for the computation of the influence coefficients (nodes and panels).
    real(rp)                                :: CIPartiel                ! Parameter to know if the partial computation of the influence coefficient has to be done or not.
    real(rp),allocatable, dimension(:,:)    :: A                        ! A matrix of the linear system.
    real(rp),allocatable, dimension(:)      :: B, Sol                   ! B matrix and the solution of the linear system.
    logical, dimension(:),allocatable       :: LTab                     ! Table to know if a node moved enough to compute its influence coefficients again (= false) or not (= true).
    logical, dimension(:),allocatable       :: BPoint                   ! Neighbours (nodes) of the nodes of LTab.
    logical, dimension(:),allocatable       :: BFace                    ! Panels linked to the nodes of LTab and their neighbours (panels).
    integer                                 :: ierror                   ! Error flag.
    logical                                 :: bool                     ! Computation of the CI (true) or not (false) when it was already done after the remeshing (only for jk = 1).
    integer,dimension(2)                    :: Border                   ! Boundaries of the modification of CD and CS in case of calling CoeffInfl_Line_Full.
        
    integer :: i_test, compt ! a supprimer
    real(rp) :: ta,tb,tc,td,te

    
    
    call CPU_TIME(ta)
    
    ! This subroutine computes the influence coefficients, creates the linear system and solves it.

    ! CD and CS have already a value, consequently, in case of partial computation of the influence coefficients, only coefficients will be updates. The other ones will keep their value.
    
    ierror = 0
    
    ! Allocation.
    N = Mesh%Nsys
    allocate(A(N,N), B(N), Sol(N))

    ! Initialization of Opt.
    ! Opt = false for the BVP on Phi and true for the BVP on DPhiDt (forced motion).
    Opt = .false.
    if(present(Option))then
        ! Case of the second call of solBVP for a forced mesh.
        if(Option) Opt = .true.
    end if
    
    ! Initialization of CIPartiel.
    ! For the second BVP, t = ti > 0, therefore CIPartiel > 0.
    CIPartiel = 1._RP
    if(present(t))then
        CIPartiel = t 
    end if
    
    ! Initialization of bool (not computation of the influence coefficient matrices when boolRemesh = true (in case of remeshing when jk = 1).
    bool = .true.
    if(present(boolRemesh))then
        if(boolRemesh)then ! Remeshing was done.
            bool = .false.
        end if
    end if
    
    
        
    ! Computation of the influence coefficients
    if(CCI .and. DeformMesh .and. bool)then ! Partial computation of the influence coefficients, the mesh is moving and no remeshing.

        
        if(CIPartiel.lt.0)then ! If the partial computation of the influence coefficients is activated.
            allocate(LTab(Mesh%Nnoeud))
            Ltab = .true.
            allocate(BPoint(Mesh%Nnoeud))
            BPoint = .false.
            allocate(BFace(Mesh%Nfacette))
            BFace = .false.
    
            borne = reshape((/1, 1, Mesh%Nsys, Mesh%Nfsys/), (/2,2/)) !  All the nodes of the system.
            
            ! Filling LTab.
            call DeplNoeud(Mesh, LTab, borne,Mesh%Nnoeud) ! LTab(j) = false : computation in using CoeffInfl_Line.
            
            

            ! Filling of BPoint and BFace.
            call ZoneInfl(Mesh, Ltab, borne, BPoint, BFace,Mesh%Nnoeud,Mesh%Nfacette) ! BPoint(j) or BFace(j) = true : computation in using CoeffInfl_Col.

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !           Free surface panels <-> Free surface panels
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! If the free surface is meshed.
            if(Mesh%FS%CMD(2))then 
                
                ! Boundaries for the free surface nodes.
                borne = reshape(Mesh%FS%IndFS(1:4), (/2,2/))
                
                ! Rows.
                if(NThreads.eq.1)then ! Sequantial.
                    do jk = borne(1,1),borne(1,2) ! Nodes.
                        if(not(LTab(jk)))then ! If the node is taken into account (LTab = false)
                        
                            ! Initialization of the lines of CD and CS for the nodes of the free surface.
                            CS(jk,borne(1,1):borne(1,2)) = 0._rp
                            CD(jk,borne(1,1):borne(1,2)) = 0._rp
                        
                            ! Computation of the line which matches with the node.
                            call CoeffInfl_Line(Mesh, CS(jk,borne(1,1):borne(1,2)), CD(jk,borne(1,1):borne(1,2)),&
                            &                                     borne(1,2)-borne(1,1)+1,jk, borne(1,1), borne)
                        end if
                    end do
                else ! Parallel.
                    Border(1) = borne(1,1)
                    Border(2) = borne(1,2)
                    call CoeffInfl_Line_Full(Mesh,Nnodes, CS, CD,LTab, Border,borne(1,1),borne)
                end if
                
                ! Columns
                do jk = borne(1,1),borne(1,2) ! Nodes.
                    if(BPoint(jk))then
                        ! Initialization of the columns of CD and CS for the nodes of the free surface.
                        where(LTab(borne(1,1):borne(1,2))) CS(borne(1,1):borne(1,2),jk) = 0._rp ! where(LTab(borne(1,1):borne(1,2))) : when LTab(j) = T do CS(:,jk) = 0.
                        where(LTab(borne(1,1):borne(1,2))) CD(borne(1,1):borne(1,2),jk) = 0._rp ! where(LTab(borne(1,1):borne(1,2))) : when LTab(j) = T do CD(:,jk) = 0.
                    end if
                end do
                
                !if(NThreads.eq.1)then ! Sequantial.
                do jk = borne(2,1),borne(2,2) ! Panels.
                    if(BFace(jk))then
                        j1 = Mesh%Tfacette(jk)%Tnoeud(1)
                        j2 = Mesh%Tfacette(jk)%Tnoeud(2)
                        j3 = Mesh%Tfacette(jk)%Tnoeud(3)
                        
                        ! Computation of the column which matches with the node.
                        call CoeffInfl_Col(Mesh,&
                        &                    CD(borne(1,1):borne(1,2),j1), CS(borne(1,1):borne(1,2),j1), BPoint(j1),&
                        &                    CD(borne(1,1):borne(1,2),j2), CS(borne(1,1):borne(1,2),j2), BPoint(j2),&
                        &                    CD(borne(1,1):borne(1,2),j3), CS(borne(1,1):borne(1,2),j3), BPoint(j3),&
                        &                    borne(1,2)-borne(1,1)+1,N,jk, Ltab, borne)
                    end if
                end do
                !else ! Parallel.
                    !call CoeffInfl_Col_Full(Mesh, CD, CS, Nnodes,Mesh%Nfacette, LTab,BPoint,BFace,borne)
                !end if
                
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !              Free surface panels <-> Body panels
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! All the nodes of both the floater and the free surface meshes are taken into account here.
            
            do j = 1,Mesh%NBody
                
                if(Mesh%Body(j)%Active)then
                    
                    ! If both the free surface and the floaters are meshed.    
                    if(Mesh%FS%CMD(2) .or. Mesh%Body(j)%CMD(2) .or. Mesh%Body(j)%CMD(1))then
                        
                        ! Influence of the nodes of the floater on the nodes of the free surface.
                        borne = reshape((/Mesh%Body(j)%IndBody(1), Mesh%FS%IndFS(2), Mesh%Body(j)%IndBody(3), Mesh%FS%IndFS(4)/), (/2,2/))
                        CD(borne(1,1):borne(1,2) , Mesh%FS%IndFS(1):Mesh%FS%IndFS(3)) = 0._RP
                        CS(borne(1,1):borne(1,2) , Mesh%FS%IndFS(1):Mesh%FS%IndFS(3)) = 0._RP
                        call CoeffInfl(Mesh, CD, CS,N, borne)
                        
                        ! Influence of the nodes of the free surface on the nodes of the floater.
                        borne = reshape((/Mesh%FS%IndFS(1), Mesh%Body(j)%IndBody(2), Mesh%FS%IndFS(3), Mesh%Body(j)%IndBody(4)/), (/2,2/))
                        CD(borne(1,1):borne(1,2), Mesh%Body(j)%IndBody(1):Mesh%Body(j)%IndBody(3)) = 0._RP
                        CS(borne(1,1):borne(1,2), Mesh%Body(j)%IndBody(1):Mesh%Body(j)%IndBody(3)) = 0._RP
                        call CoeffInfl(Mesh, CD, CS,N, borne)
                        
                    endif
                
                end if
                
            end do
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !                  Body panels <-> Body panels
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            do j = 1,Mesh%NBody
                do k = 1,Mesh%NBody
                                        
                    if(Mesh%Body(j)%Active .and. Mesh%Body(k)%Active)then
                        
                        ! If the floaters are meshed.
                        if(Mesh%Body(k)%CMD(2) .or. Mesh%Body(j)%CMD(2) .or.&
                            &          j.ne.k .and. (Mesh%Body(k)%CMD(1) .or. Mesh%Body(j)%CMD(1)) )then
                            
                            ! Boundaries for the floater nodes.
                            borne = reshape((/Mesh%Body(j)%IndBody(1), Mesh%Body(k)%IndBody(2), Mesh%Body(j)%IndBody(3), Mesh%Body(k)%IndBody(4)/), (/2,2/))
                        
                            ! If the relative velocity of the two floaters is great enough : computation of all the new floater influence coefficients.
                            if(not(norm2(Mesh%Body(j)%VBody(1:3)-Mesh%Body(k)%VBody(1:3)).lt.Epsilon))then
                                CD(borne(1,1):borne(1,2) , Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)) = 0._RP
                                CS(borne(1,1):borne(1,2) , Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)) = 0._RP
                                call CoeffInfl(Mesh, CD, CS,N, borne)
                                
                            else ! If the relative velocity of the two nodes is NOT great enough.
                                
                                ! In case of one floater : always this option.
                                
                                ! Rows.
                                if(NThreads.eq.1)then ! Sequantial.
                                    do jk = borne(1,1),borne(1,2)
                                        if(.not. LTab(jk)) then
                                            CS(jk , Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)) = 0._rp
                                            CD(jk , Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)) = 0._rp
                                            call CoeffInfl_Line(Mesh, CS(jk,Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)) ,&
                                            &                                                 CD(jk,Mesh%Body(k)%IndBody(1):Mesh%Body(k)%IndBody(3)), &
                                            &                                                 Mesh%Body(k)%IndBody(3)-Mesh%Body(k)%IndBody(1)+1,jk, Mesh%Body(k)%IndBody(1), borne)
                                        end if
                                    end do
                                else ! Parallel.
                                    Border(1) = Mesh%Body(k)%IndBody(1)
                                    Border(2) = Mesh%Body(k)%IndBody(3)
                                    call CoeffInfl_Line_Full(Mesh, Nnodes,CS, CD, LTab, Border, Mesh%Body(k)%IndBody(1),borne)
                                end if
                                
                                ! Columns.
                                do jk = Mesh%Body(k)%IndBody(1),Mesh%Body(k)%IndBody(3)
                                    if(BPoint(jk))then
                                        where (LTab(borne(1,1):borne(1,2))) CS(borne(1,1):borne(1,2),jk) = 0._rp
                                        where (LTab(borne(1,1):borne(1,2))) CD(borne(1,1):borne(1,2),jk) = 0._rp
                                    end if
                                end do          
                                
                                !if(NThreads.eq.1)then ! Sequantial.
                                do jk = borne(2,1),borne(2,2)
                                    if(BFace(jk))then
                                            j1 =  Mesh%Tfacette(jk)%Tnoeud(1)
                                            j2 =  Mesh%Tfacette(jk)%Tnoeud(2)
                                            j3 =  Mesh%Tfacette(jk)%Tnoeud(3)
                                            
                                            call CoeffInfl_Col(Mesh,&
                                        &                               CD(borne(1,1):borne(1,2),j1), CS(borne(1,1):borne(1,2),j1), BPoint(j1),&
                                        &                               CD(borne(1,1):borne(1,2),j2), CS(borne(1,1):borne(1,2),j2), BPoint(j2),&
                                        &                               CD(borne(1,1):borne(1,2),j3), CS(borne(1,1):borne(1,2),j3), BPoint(j3),&
                                        &                               borne(1,2)-borne(1,1)+1,N,jk, Ltab, borne)    
                                    end if
                                end do
                                !else ! Parallel.
                                    !call CoeffInfl_Col_Full(Mesh, CD, CS, Nnodes,Mesh%Nfacette, LTab,BPoint,BFace,borne) 
                                !end if
                                
                            end if
                        end if
                    end if
                end do
            end do
            
            ! Solid angle.
            call Angle_solid(CD, Mesh%Nnoeud)
            
            ! Deallocating.
            deallocate(Ltab)
            deallocate(BPoint)
            deallocate(BFace)
            
        end if
    elseif(DeformMesh .and. bool)then ! No partial computation of the influence coefficients, the mesh is moving and no remeshing.
        
        if(.not. Opt)then ! If Opt = false (every case but forced motion when solBVP is called for the second time).
            CD = 0._RP ; CS = 0._RP
            call CoeffInfl(Mesh, CD, CS,N)
            

        end if
        
    end if
        
    if(not(Opt))then
#ifdef _OPENMP
    !$ time(9) = omp_get_wtime()
#else
    call cpu_time(time(9))
#endif
    end if
    
    
    
    
    call CPU_TIME(tb)
    
    ! Building of A, B and initialization of Sol.
        
    if (grossier) then
        call calcul_B(CD, CS, Ecoulement, Mesh, B, cond1, Nnodes,N, opt)
    else
        call SystLin(CD, CS, Ecoulement, Mesh, A, B, Sol,Nnodes,N, Option)
    end if
    
    call CPU_TIME(tc)
    

    if(not(Opt))then
#ifdef _OPENMP
    !$ time(10) = omp_get_wtime()
#else
    call cpu_time(time(10))
#endif
    end if
    

    ! Solving of the linear system.
    if(iprint>0) print*,"Resolution système linéaire ..."
      
    if (grossier) then  
        sol = matmul(Ainv1,B)
    else
        if (Solv==0) then
            call GMRES(A, B, Sol, N, ierror)
        else  
            call LU(A, B, Sol, N)
        end if
    end if

    
    call CPU_TIME(td)
    if(not(Opt))then
#ifdef _OPENMP
    !$ time(11) = omp_get_wtime()
#else
    call cpu_time(time(11))
#endif
    end if
    
    if(ierror/=0) goto 9999
    
    ! Distribution of the solution.
    if(iprint>0) print*,"Mise à jour champ ecoulement ..." 
    call postSL(Sol,N, Mesh, Ecoulement, Option)
    
    9999 continue
        if(ierror/=0)then
            write(*,90),ierror
        end if
    90 format('error #',i3,' : solBVP: pb. in solving the Boundary Value Problem.')

    
    call CPU_TIME(te)
    
   ! write(1111,*) "Youpiiiiii : ",tb-ta,tc-tb,td-tc,te-td
    ! Deallocating.
    deallocate(A, B, Sol)

end subroutine solBVP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine Calcul_A(CD, CS, Ecoulement, Mesh, A, cond, Nnodes,Nsys)
    !!!!! Problème :
    !   Initialiser le système linéaire à partir du type de frontière
    !       o Surface Libre: On cherche la vitesse normale --> A = CS et B = B + CD*Phi
    !       o Surfaces Matérielles: On cherche le potentiel --> A = -CD et B = B - CS*DPhiDn
    !!!!!
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement) , intent(in)                  :: Ecoulement                       ! Flow parameters.
    integer,intent(in)                              :: Nnodes,Nsys                      ! Number of nodes and unknowns in the linear system.
    real(rp), dimension(Nnodes,Nnodes), intent(in)  :: CD, CS                           ! Influence coefficients.
    real(rp), dimension(Nsys,Nsys), intent(out)     :: A                                ! Matrix A of the linear system.

    integer                                         :: j, k                             ! Loop parameters.
    integer, dimension(3,2)                         :: N                                ! Boundaries for the nodes and the panels.
    real(rp), dimension(Nnodes)                     :: cond                             ! Inverse of the diagonal coefficients.
    
    ! This subroutine builds A

        
    ! Only the first body can be above the free surface, otherwise a bug will appear.
    do j = Int_Body,Mesh%NBody
        if(not(Mesh%Body(j)%Active) .and. j.ne.Mesh%NBody)then
            print*,""
            print*,"Only the last floater can be above the free surface. Otherwise a bug will appear in the 1st BVP problem. See SystLin in BVP.f90."
            pause
        end if
    end do
    
    ! Boundaries for the nodes and the panels.
    if(cuve_ferme)then
        if(Mesh%Body(Mesh%NBody)%Active)then
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is not above the free surface.
        else
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody-1)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is above the free surface.
        end if
    else
        N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(Int_Body)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) 
    end if

    ! Building A
    A(N(1,1):N(1,2),N(2,1):N(2,2)) = CS(N(1,1):N(1,2),N(2,1):N(2,2)) ! CS(:,FS)
    A(N(1,1):N(1,2),N(3,1):N(3,2)) = -CD(N(1,1):N(1,2),N(3,1):N(3,2)) ! -CD(:,Ext) -CD(:,B0) ... -CD(:,Bn)
  

    do j = 1,N(1,2)
	    if (abs(A(j,j)).GT.Epsilon) then
	        cond(j) = 1._RP/A(j,j)
	    else
		    print*,"Preconditioning of the 1st BEM problem: the coefficient ",j," of the diagonal of A is nul!"
            pause 
        end if
    end do
    
        
    do j = 1,N(1,2)
        do k = 1,N(1,2)
            A(k,j) = A(k,j)*cond(k)
        end do
    end do
        
end subroutine Calcul_A


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine Calcul_B(CD, CS, Ecoulement, Mesh, B, cond, Nnodes,Nsys, opt)
    !!!!! Problème :
    !   Initialiser le système linéaire à partir du type de frontière
    !       o Surface Libre: On cherche la vitesse normale --> A = CS et B = B + CD*Phi
    !       o Surfaces Matérielles: On cherche le potentiel --> A = -CD et B = B - CS*DPhiDn
    !!!!!
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement) , intent(in)                  :: Ecoulement                       ! Flow parameters.
    integer,intent(in)                              :: Nnodes,Nsys                      ! Number of nodes and unknowns in the linear system.
    real(rp), dimension(Nnodes,Nnodes), intent(in)  :: CD, CS                           ! Influence coefficients.
    real(rp), dimension(Nsys),intent(out)           :: B                                ! Matrix A of the linear system.

    integer                                         :: j, k                             ! Loop parameters.
    integer, dimension(3,2)                         :: N                                ! Boundaries for the nodes and the panels.
    real(rp), dimension(Nsys)                       :: cond                             ! Inverse of the diagonal coefficients.
    real(rp), allocatable                           :: Phi(:), DphiDn(:)                ! Phi, DPhiDn and inverse
   logical, intent(in)                              :: Opt                              ! Opt = false for the BVP on Phi and true for the BVP on DPhiDt (forced motion).
    
    ! This subroutine builds B

        
    ! Only the first body can be above the free surface, otherwise a bug will appear.
    do j = Int_Body,Mesh%NBody
        if(not(Mesh%Body(j)%Active) .and. j.ne.Mesh%NBody)then
            print*,""
            print*,"Only the last floater can be above the free surface. Otherwise a bug will appear in the 1st BVP problem. See SystLin in BVP.f90."
            pause
        end if
    end do
    
    ! Boundaries for the nodes and the panels.
    if(cuve_ferme)then
        if(Mesh%Body(Mesh%NBody)%Active)then
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is not above the free surface.
        else
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody-1)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is above the free surface.
        end if
    else
        N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(Int_Body)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) 
    end if

    ! Building B
    allocate(Phi(N(1,2)),DPhiDn(N(1,2)))
    
    if(opt)then
        ! BVP on DPhiDt.
        Phi(N(2,1):N(2,2)) = Ecoulement%DPhiDt(N(2,1):N(2,2))%perturbation
        DPhiDn(N(3,1):N(3,2)) = Ecoulement%DDPhiDnDt(N(3,1):N(3,2))%perturbation
    else
        ! BVP on Phi.
        Phi(N(2,1):N(2,2)) = Ecoulement%Phi(N(2,1):N(2,2))%perturbation
        DPhiDn(N(3,1):N(3,2)) = Ecoulement%DPhiDn(N(3,1):N(3,2))%perturbation
    end if

    Phi(N(2,1):N(2,2)) = Ecoulement%Phi(N(2,1):N(2,2))%perturbation
    DPhiDn(N(3,1):N(3,2)) = Ecoulement%DPhiDn(N(3,1):N(3,2))%perturbation


    B = 0._RP
    do k = N(2,1),N(2,2) ! Free surface
        do j = 1,N(1,2)
            B(j) = B(j) + CD(j,k)*Phi(k) ! CD(:,FS)*Phi(FS)
        end do
    end do
        
    do k = N(3,1),N(3,2) ! Bodies
        do j = 1,N(1,2)
            B(j) = B(j) - CS(j,k)*DPhiDn(k) ! -CS(:,Ext)*Phi_n(Ext) -CS(:,B0)*Phi_n(B0) ... -CS(:,Bn)*Phi_n(Bn)
        end do    
    end do

    deallocate(Phi,DPhiDn)
    
    B = B*cond
        
end subroutine Calcul_B


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine SystLin(CD, CS, Ecoulement, Mesh, A, B, Sol,Nnodes,Nsys, Option)
    !!!!! Problème :
    !   Initialiser le système linéaire à partir du type de frontière
    !       o Surface Libre: On cherche la vitesse normale --> A = CS et B = B + CD*Phi
    !       o Surfaces Matérielles: On cherche le potentiel --> A = -CD et B = B - CS*DPhiDn
    !!!!!
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement) , intent(in)                  :: Ecoulement                       ! Flow parameters.
    integer,intent(in)                              :: Nnodes,Nsys                      ! Number of nodes and unknowns in the linear system.
    real(rp), dimension(Nnodes,Nnodes), intent(in)  :: CD, CS                           ! Influence coefficients.
    real(rp), dimension(Nsys,Nsys), intent(out)     :: A                                ! Matrix A of the linear system.
    real(rp), dimension(Nsys), intent(out)          :: B, Sol                           ! Vector B of the linear system.
    logical, intent(in), optional                   :: Option                           ! present(Option) = false for the BVP on Phi, present(Option) = true and Option = true for the BVP on DPhiDt.
    
    logical                                         :: Opt                              ! = false : BVP on Phi, = true : BVP on DPhiDt.
    integer                                         :: j, k                             ! Loop parameters.
    integer, dimension(3,2)                         :: N                                ! Boundaries for the nodes and the panels.
    real(rp), allocatable                           :: Phi(:), DphiDn(:), inv_diag(:)   ! Phi, DPhiDn and inverse of the diagonal coefficients.
    
    ! This subroutine builds A, B and initializes the solution.
    
    ! Initialization of Opt.
    ! Opt = false for the BVP on Phi and true for the BVP on DPhiDt (forced motion).
    Opt = .false.
    if(present(Option))then
        ! Case of the second call of solBVP for a forced mesh.
        if(Option) Opt = Option
    end if
        
    ! Only the first body can be above the free surface, otherwise a bug will appear.
    do j = Int_Body,Mesh%NBody
        if(not(Mesh%Body(j)%Active) .and. j.ne.Mesh%NBody)then
            print*,""
            print*,"Only the last floater can be above the free surface. Otherwise a bug will appear in the 1st BVP problem. See SystLin in BVP.f90."
            pause
        end if
    end do
    
    ! Boundaries for the nodes and the panels.
    if(cuve_ferme)then
        if(Mesh%Body(Mesh%NBody)%Active)then
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is not above the free surface.
        else
            N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody-1)%IndBody(3)],(/3,2/)) ! Body(Mesh%NBody) is above the free surface.
        end if
    else
        N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(Int_Body)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/)) 
    end if

    ! Building A
    A(N(1,1):N(1,2),N(2,1):N(2,2)) = CS(N(1,1):N(1,2),N(2,1):N(2,2)) ! CS(:,FS)
    A(N(1,1):N(1,2),N(3,1):N(3,2)) = -CD(N(1,1):N(1,2),N(3,1):N(3,2)) ! -CD(:,Ext) -CD(:,B0) ... -CD(:,Bn)
    
    ! Building B
    allocate(Phi(N(1,2)),DPhiDn(N(1,2)))
    
    ! Opt = false : BVP on Phi, = true : BVP on DPhiDt.
    if(Opt)then
        ! BVP on DPhiDt.
        do j=N(2,1),N(2,2)
            !Phi(j) = Ecoulement%DPhiDt(j)%perturbation - dot_product(Mesh%Tnoeud(j)%Velocity,Ecoulement%GPhi(:,j)%perturbation)
            Phi(j) = Ecoulement%DPhiDt(j)%perturbation
        enddo
        DPhiDn(N(3,1):N(3,2)) = Ecoulement%DDPhiDnDt(N(3,1):N(3,2))%perturbation
    else
        ! BVP on Phi.
        Phi(N(2,1):N(2,2)) = Ecoulement%Phi(N(2,1):N(2,2))%perturbation
        DPhiDn(N(3,1):N(3,2)) = Ecoulement%DPhiDn(N(3,1):N(3,2))%perturbation
    end if
    
    if (.false.) then
        do j = 1,N(1,2)
            B(j) = dot_product(CD(j,N(2,1):N(2,2)),Phi(N(2,1):N(2,2))) - dot_product(CS(j,N(3,1):N(3,2)),DPhiDn(N(3,1):N(3,2)))
        end do
    else
        B = 0._RP
        do k = N(2,1),N(2,2) ! Free surface
            do j = 1,N(1,2)
                B(j) = B(j) + CD(j,k)*Phi(k) ! CD(:,FS)*Phi(FS)
            end do
        end do
        
        do k = N(3,1),N(3,2) ! Bodies
            do j = 1,N(1,2)
                B(j) = B(j) - CS(j,k)*DPhiDn(k) ! -CS(:,Ext)*Phi_n(Ext) -CS(:,B0)*Phi_n(B0) ... -CS(:,Bn)*Phi_n(Bn)
            end do    
        end do
    end if
    deallocate(Phi,DPhiDn)
    
    ! Building the solution.
    if(Opt)then
        ! BVP on DPhiDt.
        Sol(N(2,1):N(2,2)) = Ecoulement%DDPhiDnDt(N(2,1):N(2,2))%perturbation
        Sol(N(3,1):N(3,2)) = Ecoulement%DPhiDt(N(3,1):N(3,2))%perturbation
    else
        ! BVP on Phi.
        Sol(N(2,1):N(2,2)) = Ecoulement%DPhiDn(N(2,1):N(2,2))%perturbation ! Phi_n(FS)
        Sol(N(3,1):N(3,2)) = Ecoulement%Phi(N(3,1):N(3,2))%perturbation ! Phi(Ext) Phi(B0) ... Phi(Bn)
    end if
    
    ! Preconditionnement de la matrice du système linéaire.
    allocate(inv_diag(N(1,2)))
    do j = 1,N(1,2)
	    if (abs(A(j,j)).GT.Epsilon) then
	        inv_diag(j) = 1._RP/A(j,j)
		    B(j) = B(j)*inv_diag(j)
	    else
		    print*,"Preconditioning of the 1st BEM problem: the coefficient ",j," of the diagonal of A is nul!"
            pause 
        end if
    end do
    
    do j = 1,N(1,2)
        do k = 1,N(1,2)
            A(k,j) = A(k,j)*inv_diag(k)
        end do
    end do
    
    deallocate(inv_diag)
        
end subroutine SystLin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                               !
!                                     Post-Process de la Résolution                                             !
!                                                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine PostSL(X,size_X, Mesh, Ecoulement, Option)
    !!!!! Problème :
    !   Travail inverse de la fonction Ini, on réorganise les solutions à leur place dans Ecoulement   
    !       o Surface Libre : DPhiDn = solution
    !       o Surface Matérielle : Phi = solution   
    !!!!!
    
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage), intent(in)             :: Mesh         ! Mesh.
    integer,intent(in)                      :: size_X       ! Size of the solution.
    Real(RP), dimension(size_X), intent(in) :: X            ! Solution.
    !f2py integer*1, dimension(1000)        :: Ecoulement
    type(TEcoulement), intent(inout)        :: Ecoulement   ! Flow parameters.
    logical, intent(in), optional           :: Option       ! present(Option) = false for the BVP on Phi, present(Option) = true and Option = true for the BVP on DPhiDt.
    
    integer                                 :: j, k         ! Loop parameters.
    integer, dimension(3,2)                 :: N            ! Boundaries for the nodes and the panels.
    logical                                 :: Opt          ! = false : BVP on Phi, = true : BVP on DPhiDt.
    
    ! This subroutine distributes the solution of the linear system (1st BVP or 2nd one in case of forced motion).
    
    ! Initialization of Opt.
    ! Opt = false for the BVP on Phi and true for the BVP on DPhiDt (forced motion).
    Opt = .false.
    if(present(Option))then
        ! Case of the second call of solBVP for a forced mesh.
        if (Option) Opt = .true.
    end if
    
    ! Boundaries for the nodes and the panels.
    if(cuve_ferme)then
        N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(1)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Mesh%NBody)%IndBody(3)],(/3,2/))
    else
        N = reshape([1,Mesh%FS%IndFS(1),Mesh%Body(Int_Body)%IndBody(1),Mesh%Nsys,Mesh%FS%IndFS(3),Mesh%Body(Int_Body)%IndBody(3)],(/3,2/))
    end if
    
    ! Distribution of the solution.
    if(Opt)then
        Ecoulement%DDPhiDnDt(N(2,1):N(2,2))%Perturbation = X(N(2,1):N(2,2)) 
        do j = N(3,1),N(3,2)
            Ecoulement%DPhiDt(j)%perturbation = X(j) ! + dot_product(Mesh%Tnoeud(j)%Velocity,Ecoulement%GPhi(:,j)%perturbation)
        enddo
    else
        Ecoulement%DPhiDn(N(2,1):N(2,2))%Perturbation = X(N(2,1):N(2,2))
        Ecoulement%Phi(N(3,1):N(3,2))%Perturbation = X(N(3,1):N(3,2))
    end if
    
    ! Continuity of Phi of DPhiDt for the twin nodes.
    if (cuve_ferme) then
        do j = 1,Mesh%Nsys
            if (Mesh%Tnoeud(j)%TypeNoeud.eq.1 .and. Mesh%Tnoeud(j)%Ndouble.ne.0) then
                do k = 1,Mesh%Tnoeud(j)%Ndouble
                    if (Mesh%Tnoeud(Mesh%Tnoeud(j)%double(k))%TypeNoeud.ne.0 .and. j.lt.Mesh%Tnoeud(j)%double(k) ) then
                        if (abs(Ecoulement%Phi(j)%perturbation-Ecoulement%Phi(Mesh%Tnoeud(j)%double(k))%perturbation).gt.Epsilon2) then
                            print*,'PostSL: Phi(',j,') = ', Ecoulement%Phi(j)%perturbation, ' , Phi(',Mesh%Tnoeud(j)%double(k),') = ', Ecoulement%Phi(Mesh%Tnoeud(j)%double(k))%perturbation
                        elseif (abs(Ecoulement%Phi(j)%perturbation-Ecoulement%Phi(Mesh%Tnoeud(j)%double(k))%perturbation).gt.Epsilon) then
                            print*,'PostSL: Ajustement Phi au point ', j, abs(Ecoulement%Phi(j)%perturbation- Ecoulement%Phi(Mesh%Tnoeud(j)%double(k))%perturbation)
                        end if
                    end if
                end do
            end if
        end do
    end if
    
    return
    
end subroutine PostSL

end module BVP
    
    
    
    
    
