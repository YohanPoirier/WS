module BodyMotion_mod
use Constantes
use Parameters
use FonctionsCommunes
use Spline
use PrePlot
use Exec_mod
!$ use OMP_LIB
implicit none
logical, parameter :: forcage_hydrostatique_lineaire = .false.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                            FreeBodyMotion                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FreeBodyMotion(Mesh, Ecoulement, CD, CS,t,Nnodes,InputData,time)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage),intent(inout)                   :: Mesh                     ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement               ! Flow parameters (phi, etc.)
    integer,intent(in)                              :: Nnodes                   ! Number of nodes in the mesh (size of CD and CS)
    real(rp), dimension(Nnodes,Nnodes), intent(in)  :: CD, CS                   ! Influence coefficients
    real(rp), intent(in)                            :: t                        ! Current time (of the RK4 step)
    !f2py integer*1, dimension(1000)                :: InputData                
    type(InputDataStruct),intent(in)                :: InputData                ! Input data
    
    real(rp),dimension(20) :: time
    
    integer                                         :: Ndof, N, Nb,size_NBody   ! Number of dof, size of A, number of nodes in the mesh of the floaters and size of the structure NBody
    real(rp), dimension(:,:), allocatable           :: A                        ! A matrix of the linear system
    real(rp), dimension(:), allocatable             :: B, Sol                   ! B and X matrices of the linear system
    integer                                         :: jj,nc                    ! Loop parameter
    integer                                         :: ierror                   ! Error flag
    
    ! This subroutine solves the coupling between the hydrodynamics and the dynamics;
    
    ierror = 0
    Ndof = 0
    Nb = 0
    jj = 1
    do nc = Int_Body,Mesh%Nbody
        if (Mesh%Body(nc)%CMD(1)) then
            Nb = Nb + Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1 ! Total number of nodes for all the bodies (not the tank).
            Ndof = Ndof + InputData%Ndll(jj)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
    N = Mesh%Nsys + Nb + Ndof ! Size of the linear system.
    size_NBody = Mesh%NBody
    allocate(A(N,N), B(N), Sol(N))
    
    ! Computation of q'%perturbation
    call AccConvec(Mesh, Ecoulement, t)
            
    ! Computation of the intertia matrix of the body
    call Inertia(Mesh,InputData,.false.)
    
    ! Initialisation of the linear system
    call SystLin_FreeBody(CD, CS, Ecoulement, Mesh, A, B, Sol, t, size_NBody,Nnodes,N,InputData)
    
    !$ time(15) = omp_get_wtime()
    
    ! Resolution of the linear system
    
    
    if (bool_coarse) then  
        sol = matmul(Mesh%Ainv2,B)
    else
        if (Solv.eq.0) then
            call GMRES(A, B, Sol, N,ierror)
        else  
            call LU(A, B, Sol, N, ierror)
        end if
    end if
    
    if(ierror/=0) goto 9999
    
    !$ time(16) = omp_get_wtime()
    
    ! Distribution of the solutions
    call postSL_FreeBody(Sol, N,Mesh, Ecoulement,N,InputData)
    
    deallocate(A, B, Sol)

    9999 continue
        if(ierror/=0)then
            write(*,90),ierror
        endif
    90 format('error #',i3,' : pb. in resolution of free motion.')

end subroutine FreeBodyMotion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                            SystLin_FreeBody                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SystLin_FreeBody(CD, CS, Ecoulement, Mesh, A, B, Sol, t,size_NBody, Nnodes, N,InputData)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(inout)                  :: Mesh                 ! Total mesh (domain and floater).
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement           ! Flow parameters (phi, etc.).
    integer,intent(in)                              :: Nnodes               ! Number of nodes in the mesh.
    integer,intent(in)                              :: N                    ! Number of unknowns in the BVP problem.
    real(rp), dimension(Nnodes,Nnodes), intent(in)  :: CD, CS               ! Influence coefficients.
    real(rp), dimension(N,N), intent(inout)         :: A                    ! Matrix A of the linear system.
    real(rp), dimension(N), intent(inout)           :: B                    ! Matrix B of the linear system and the initialization of the solution of it.
    real(rp),dimension(N),intent(inout)             :: Sol
    real(rp),intent(in)                             :: t                    ! Current time.
    integer,intent(inout)                           :: size_NBody           ! Size of NBody.
    !f2py integer*1, dimension(1000)                :: InputData                
    type(InputDataStruct),intent(in)                :: InputData            ! Input data.
    
    integer                                         :: Nsys,Nb              ! Number of unknowns in the mesh, nodes in the mesh of the floaters and the number of degree of freedom.
    integer, dimension(3)                           :: Nsl                  ! Index of the nodes for the free surface.
    integer, dimension(size_NBody,3)                :: NBody                ! Index of the nodes for the floater.
    real(rp)                                        :: TimeRamp             ! Ramp parameter.
    real(rp), allocatable                           :: Fext(:,:)            ! External forces.
    real(rp), allocatable                           :: S(:,:,:),DSDt(:,:,:) ! Transformation matrices between the Cardan frames and the inertial frame and its time-differentiation.
    real(rp), allocatable                           :: M(:,:,:)             ! Mass matrices.
    real(rp), allocatable                           :: InertiaTerms(:,:)    ! Inertia terms in the dynamical momentum expression.
    real(rp),dimension(3,3,2)                       :: Tob                  ! eR0.
    

    
    ! This subroutine creates the linear system of the coupling between the hydrodynamics and the dynamics.
        
    ! Time ramp
    call Ramp_Management(t,TimeRamp)
    
    ! Number of nodes and index
    call Index_Number_Nodes(Mesh,Nsl,NBody,Nsys,Nb,size_NBody)
    
    ! Allocation of CK, CT, Q, Th, S, DSDt, Fext, M, InertiaTerms
    call Init_CK_CT_Q_Th(Mesh,InputData)
    allocate(S(3,3,NBodies))
    allocate(DSDt(3,3,NBodies))
    allocate(Fext(6,NBodies))
    allocate(M(6,6,NBodies))
    allocate(Inertiaterms(3,NBodies))
    
    ! Computation of S
    call get_S(Mesh,S,NBodies)
    
    ! Computation of DSDt
    call get_DSDt(Mesh,DSDt,NBodies)
    
    ! Computation of CK
    call get_CK(Mesh)
    
    ! Computation of CT and updating of CK
    call get_CT(Mesh,S,NBodies)
        
    ! Computation of the inertia terms used in the dynamic momentum
    call get_Inertia_terms(Mesh,M,S,DSDt,InertiaTerms,InputData,NBodies,Tob)
        
    ! Computation of Q
    call get_Q(Mesh,Ecoulement,TimeRamp)
    
    ! External forces
    call ExternalForces(Mesh, Ecoulement, Fext, t, 1,InputData,NBodies) ! In a pure WSC simulation, the weight is always computed.
    
    ! Computation of Th
    call get_Th(Mesh,Ecoulement,Fext,InertiaTerms,TimeRamp,InputData,NBodies)
    
    ! Building of A
    call Building_A(Mesh,A,CS,CD,M,Nsl,NBody,size_NBody,N,Nsys,Nb,Nnodes,InputData,NBodies)
        
    
    ! Building of B
    call Building_B(Mesh,Ecoulement,B,CD,Nsys,Nb,Nsl,Nnodes,N,InputData,NBody,size_NBody)
        
    
    ! Initialization of the solution with the one at the last time step
    call Initialization_solution(Mesh,Ecoulement,Nsys,Nb,Nbody,size_NBody,Sol,N,InputData)
    
    ! Preconditioning of the linear system
    call Preconditioning(A,B,N)
    
    ! Deallocating of the useless coefficients
    call deallocate_CK_CT_Q_Th(Mesh)
    deallocate(S,DSDt,Fext,M,Inertiaterms)
    
end subroutine SystLin_FreeBody

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                            ExternalForces                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ExternalForces(Mesh, Ecoulement, Fext, t, WeightOrNot,InputData,size_Fext)

    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                     ! Total mesh (domain and floater).
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement), intent(in)                   :: Ecoulement               ! Flow parameters (phi, etc.).
    integer,intent(in)                              :: size_Fext                ! NBodies.
    real(rp), dimension(6,size_Fext), intent(inout) :: Fext                     ! External forces.
    real(rp), intent(in)                            :: t                        ! Current time.
    integer, intent(in)                             :: WeightOrNot              ! = 1 if the weight is computed, = 0 otherwise (useful for the coupling between InWave and the WSC code).
    !f2py integer*1, dimension(1000)                :: InputData                
    type(InputDataStruct),intent(in)                :: InputData                ! Input data.
    
    integer :: j
    real(rp)                                        :: TimeRamp, Surface, ds    ! Time ramp, surface, ds.
    real(rp), dimension(3)                          :: AG, Vhoule               ! AG and the velocity of the waves.
    real(rp), dimension(6)                          :: Fadd                     ! Extra loads.
    integer                                         :: jj,nc
    
    ! This subroutine computes the external forces.
    
    Fext = 0._RP
    
    ! Time ramp.
    call Ramp_Management(t,TimeRamp)
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
        Fadd = 0._RP
        
        ! Weight.

        if (WeightOrNot == 1)then
            if (hydrostatique .and. not(lineaireBody .or. forcage_hydrostatique_lineaire)) then
                Fext(1:3,jj) = [0._rp, 0._rp,-Mesh%Body(nc)%Mass*g] ! Force
       
                !Fext(1:3,jj) = Fext(1:3,jj) + [0._rp, 0._rp, 1.772184543_rp] !  a suppr (1111)

            end if
                        
            
            if(iFixPoint)then
                ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                AG = Mesh%Body(nc)%GBody(1:3) - FixPointPos
            else
                ! Motion equation solved at G (CSolv = (GBody,angles)).
                AG = 0._RP
            end if
            call Computation_vect_product(AG,Fext(1:3,jj),Fext(4:6,jj)) ! Moment
        end if
        
        ! PTO
        do j = 1,3 ! Force.
            if(InputData%dll_dir(j,jj))then
                Fadd(j) = - InputData%MuPTO(jj)*Mesh%Body(nc)%VBody(j)
                if((Mesh%Body(nc)%CSolv(j) - InputData%Pressort(j,jj)).ge.InputData%Lressort(jj))then ! No compression. (Not the case at the origin of the CETO test case).
                    Fadd(j) = Fadd(j) - InputData%Raideur(jj)*(Mesh%Body(nc)%CSolv(j) - InputData%Pressort(j,jj) - InputData%Lressort(jj))
                end if                
            else    
                Fadd(j) = 0._rp
            end if
        end do
        

        do j = 4,6 ! Moment.
            if(InputData%dll_dir(j,jj))then
                Fadd(j) = -InputData%MuPTO(jj)*Mesh%Body(nc)%VBody(j)
                Fadd(j) = Fadd(j) - InputData%Raideur(jj)*Mesh%Body(nc)%CSolv(j)                                
            else    
                Fadd(j) = 0._rp
            end if
        end do
                
        Fext(1:6,jj) = Fext(1:6,jj) + Fadd
        
        ! Morison load.
        if (.true.) then
            Surface = 2._RP/3._RP*PI*Mesh%Body(nc)%DimBody(jj)**2
            if (symmetry) Surface = 0.5_RP*Surface
            Vhoule = 0._RP
            Fadd(1:3) = -0.5*ro*InputData%Cd_Morison(jj)*Surface*norm2(Mesh%Body(nc)%VBody(1:3)-Vhoule)*(Mesh%Body(nc)%VBody(1:3)-Vhoule)
            Fadd(4:6) = 0._RP
        else
            Fadd(1:3) = 0._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                ds = Mesh%Tnoeud(j)%Aire*dot_product(Mesh%Tnoeud(j)%Normale,Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)
                Fadd(1:3) = Fadd(1:3) + ds*norm2(Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)*(Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)
            end do
            Fadd(1:3) = 0.5_RP*ro*InputData%Cd_Morison(jj)*Fadd(1:3)
            Fadd(4:6) = 0._RP
        end if
        

        
        ! Total external loads.
        Fext(1:6,jj) = Fext(1:6,jj) + Fadd

        jj = jj + 1
        
        
        
    end do

end subroutine ExternalForces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Calcul terme convection condition sur le corps                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AccConvec(Mesh, Ecoulement, t) 
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                                         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement)                   :: Ecoulement                                   ! Flow parameters.
    real(rp), optional                  :: t                                            ! Current time.
    
    integer                             :: j, nc,jj                                     ! Loop parameters.
    integer                             :: Int_Body0,Int_Body1                          ! Body indexes.
    real(rp)                            :: DPhiDn, TimeRamp                             ! DPhiDn and the time ramp.
    real(rp), dimension(2)              :: Curv                                         ! Curvature.
    real(rp), dimension(3)              :: V, ThetaP, GPhi, GPhi2, OmegaMG              ! Vectors.
    real(rp), dimension(3)              :: AM, Corr, Ctre, Omega, Vrot                  ! Vectors.
    real(rp), dimension(3,3)            :: Pt, Puvw                                     ! Local basis.
    real(rp), allocatable               :: S(:,:,:), DSDt(:,:,:)                        ! Transformation matrices between the Cardan frames and the inertial frame and its time-differentiation.
    real(rp), dimension(3)              :: Vect_product_1,Vect_product_2,Vect_product_3 ! Cross products.
    
    ! This subroutine computes the body condition for BVP on Phi_t.

    allocate(S(3,3,NBodies))
    allocate(DSDt(3,3,NBodies))
    
    if(present(t))then
        call Ramp_Management(t,TimeRamp)
    else    
        TimeRamp = 1._rp
    endif
    
    ! S
    call get_S(Mesh,S,NBodies)
        
    ! DSDt
    call get_DSDt(Mesh,DSDt,NBodies)
    
    jj = 1    
    do nc = 1,Mesh%Nbody
        
        Int_Body0 = Mesh%Body(nc)%IndBody(1)
        Int_Body1 = Mesh%Body(nc)%IndBody(3)
                
        if (Mesh%Body(nc)%CMD(1)) then
            
            if(iFixPoint)then
                ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                Ctre = FixPointPos
            else
                ! Motion equation solved at G (CSolv = (GBody,angles)).
                Ctre = Mesh%Body(nc)%CSolv(1:3)
            end if
            
            ! Angular velocity
            Vrot(1:3) = Mesh%Body(nc)%VBody(4:6)
        
            ! Omega
            Omega = matmul(S(1:3,1:3,jj),Vrot)
        
            do j = Int_Body0,Int_Body1
                
                ! Correction term for rotational motion.
                AM = Mesh%Tnoeud(j)%Pnoeud(1:3) - Ctre ! AM = AM if iFixPoint, = GM if not.
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    call Computation_vect_product(Omega,AM,Vect_product_1) ! Vect_product_1 = Omega x AM or GM
                    call Computation_vect_product(Omega, Vect_product_1,Vect_product_2) ! Vect_product_2 = Omega x (Omega x AM)
                    call Computation_vect_product(matmul(DSDt(1:3,1:3,jj),VRot),AM,Vect_product_3) ! Vect_product_3 = (DSDt*Theta) x AM
                    corr = Vect_product_3 + Vect_product_2 ! (DSDt*Theta) x AM + Omega x (Omega x AM)
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    Vect_product_1 = 0._RP
                    Vect_product_2 = 0._RP
                    call Computation_vect_product(matmul(DSDt(1:3,1:3,jj),VRot),AM,Vect_product_3) ! Vect_product_3 = (DSDt*Theta) x GM
                    corr = Vect_product_3 ! (DSDt*Theta) x GM
                end if

                Puvw = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                Pt = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
                
                ! Velocity.
                call Computation_vect_product(Omega,AM,OmegaMG) ! OmegaMG = Omega x AM or GM
                V = Mesh%Body(nc)%VBody(1:3) + OmegaMG + [1.0_RP,0._RP,0._RP]*Forward_Velocity*TimeRamp
                V = matmul(V,Puvw) ! Projection of V in the local basis.
                
                ThetaP = matmul(Omega,Puvw)
                
                ! FIXME
                Curv = Mesh%Tnoeud(j)%Dlocal(3:4) ! = 1/R
            
                ! Perturbation
                GPhi = TimeRamp*matmul(Pt, Ecoulement%GPhi(1:3,j)%incident)
                GPhi(1:2) = GPhi(1:2) + Ecoulement%GPhi2(1:2,j)%perturbation
                GPhi2(1:2) = TimeRamp*Ecoulement%Gphi2(3:4,j)%incident
                GPhi2(1:2) = Gphi2(1:2) + Ecoulement%GPhi2(3:4,j)%perturbation         
                DPhiDn = Ecoulement%DPhiDn(j)%perturbation + TimeRamp*Ecoulement%DPhiDn(j)%incident
            
                ! Computation of q%perturbation
                Ecoulement%DGradPhiSqDn(j)%perturbation = (2._rp*V(1)-GPhi(1))*ThetaP(2) - (2._rp*V(2)-GPhi(2))*ThetaP(1) &
    &                + (GPhi(1)-V(1))*Curv(1)*V(1) + (GPhi(2)-V(2))*Curv(2)*V(2) + (GPhi2(1)+GPhi2(2))*V(3) + (Curv(1)+Curv(2))*DPhiDn*V(3)
            
                ! Computation of q'%perturbation
                Ecoulement%DGradPhiSqDn(j)%perturbation = Ecoulement%DGradPhiSqDn(j)%perturbation + dot_product(Corr,Mesh%Tnoeud(j)%Normale(1:3))
                
            end do
            jj = jj + 1
        else
            if(Mesh%Body(nc)%Active)then ! Not for the body above the FS.
                Ecoulement%DGradPhiSqDn(Int_Body0:Int_Body1)%incident = 0._RP
                Ecoulement%DGradPhiSqDn(Int_Body0:Int_Body1)%perturbation = 0._RP
            end if
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
            
        end if
    end do
    
    ! Deallocating.
    deallocate(S)
    deallocate(DSDt)
    
end subroutine AccConvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     Post-Process de la Résolution                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PostSL_FreeBody(X, size_X,Mesh, Ecoulement,size_NBody,InputData)
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(inout)              :: Mesh                     ! Mesh.
    integer,intent(in)                          :: size_X
    real(RP), dimension(size_X), intent(in)     :: X                        ! The solution of the linear system.
    !f2py integer*1, dimension(1000)            :: Ecoulement
    type(TEcoulement), intent(inout)            :: Ecoulement               ! Flow parameters.
    integer,intent(in)                          :: size_NBody
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(in)            :: InputData                ! Input data.

    integer                                     :: k, j,jj                  ! Loop parameters.
    integer                                     :: Nc, Nsys, Ndof, Nb, Nt   ! Parameters.
    integer, dimension(size_NBody,3)            :: NBody                    ! Number of bodies.
    
    ! This subroutine distributes the solution of the linear system.
        
    Nsys = Mesh%Nsys
    Nb = 0
    do nc=1,Mesh%NBody
        Nbody(nc,1:2) = [Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)]
        NBody(nc,3) = Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1
        if (Mesh%Body(nc)%CMD(1)) then 
            Nb = Nb + NBody(nc,3)
        end if
    end do
    
    ! Phi(FS)_tn
           
    Ecoulement%DDPhiDnDt(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Perturbation = X(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))
    
    Nt = 0
    jj = 1
    Ndof = 0
    
    do nc = 1,Mesh%NBody

        ! Phi(B0)_t and Phi(Ext)_t
        if(Mesh%Body(nc)%Active)then
            Ecoulement%DPhiDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation = X(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))
        end if

        if (Mesh%Body(nc)%CMD(1)) then
            Ecoulement%DDPhiDnDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation = X(Nsys+Nt+1:Nsys+Nt+NBody(nc,3))
            k=0
            do j = 1,6
                if(InputData%dll_dir(j,jj))then
                    k = k+1
                    ! Acceleration of the floater
                    Mesh%Body(nc)%ABody(j) = X(Nsys+Nb+Ndof+k)
                else
                    Mesh%Body(nc)%ABody(j) = 0._rp
                end if
            enddo
            
            Ndof = Ndof + InputData%ndll(jj)
            Nt = Nt + NBody(nc,3)
            jj = jj + 1
        else
            Mesh%Body(nc)%ABody = 0._RP
            
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
            
        end if
    end do
    
end subroutine PostSL_FreeBody

subroutine PostSL_FreeBody_MB(X, size_X,Mesh, Ecoulement,size_NBody)
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(inout)              :: Mesh             ! Mesh
    integer,intent(in)                          :: size_X           ! Size of the solution
    real(RP), dimension(size_X), intent(in)     :: X                ! The solution of the linear system.
    !f2py integer*1, dimension(1000)            :: Ecoulement
    type(TEcoulement), intent(inout)            :: Ecoulement       ! Flow
    integer,intent(in)                          :: size_NBody       ! Size of NBody
    
    integer                                     :: Nc, Nsys, Nb, Nt ! Parameters
    integer, dimension(size_NBody,3)            :: NBody            ! Number of bodies
    
    ! This subroutine distributes the solution of the linear system for the WSC components.
    
    Nsys = Mesh%Nsys
    Nb = 0
    do nc = 1,Mesh%NBody
        Nbody(nc,1:2) = [Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)]
        NBody(nc,3) = Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1
        if (Mesh%Body(nc)%CMD(1)) then 
            Nb = Nb + NBody(nc,3)
        end if
    end do
    
    ! Phi(FS)_tn
    Ecoulement%DDPhiDnDt(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Perturbation = X(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))
    
    Nt = 0
    do nc = 1,Mesh%NBody
        ! Phi(B0)_t and Phi(Ext)_t
        if(Mesh%Body(nc)%Active)then
            Ecoulement%DPhiDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation = X(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))
        end if
        
        if (Mesh%Body(nc)%CMD(1)) then
            Ecoulement%DDPhiDnDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation = X(Nsys+Nt+1:Nsys+Nt+NBody(nc,3))
            Nt = Nt + NBody(nc,3)
        end if
    end do

end subroutine PostSL_FreeBody_MB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                            ForceBodyMotion                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ForceBodyMotion(Mesh, Ecoulement, t, InputData)

    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh                 ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement           ! Flow parameters.
    real(rp), intent(in)                :: t                    ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData            ! Input data.
    
    integer                             :: j,nc,jj              ! Loop parameters.
    integer                             :: Int_Body0, Int_Body1 ! Body node indexes.
    real(rp),dimension(3)               :: M,AM,AMvN            ! Vectors.
    real(rp),dimension(6)               :: Acc                  ! Acceleration.
    real(rp)                            :: TimeRamp             ! Time ramp.
    real(rp), allocatable               :: S(:,:,:)             ! Transformation matrices between the Cardan frames and the inertial frame and its time-differentiation.
    
    ! This subroutine computes the body condition on Phi_t for a forced motion.
    
    ! Calculation of the advection term in the body condition of the BVP on DPhiDt
    call AccConvec(Mesh, Ecoulement, t)
    
    ! Time ramp
    call Ramp_Management(t,TimeRamp)
    
    ! S
    allocate(S(3,3,NBodies))
    call get_S(Mesh,S,NBodies)
    
    ! Conditions limites sur le corps (DDPhiDnDt)
    jj = 1
    do nc = Int_Body,Mesh%NBody
        if(Mesh%Body(nc)%Active)then
            
            ! Accelerations.
            Acc = Mesh%Body(nc)%ABody(1:6)
            do j = 1,6
                if(InputData%dll_dir(j,jj))then ! Forced dof.
                    if(abs(InputData%ACorps(jj)).gt.Epsilon)then ! Sinusodial velocity.
                        Mesh%Body(nc)%ABody(j) = -InputData%ACorps(jj)*InputData%WCorps(jj)**2*sin(InputData%WCorps(jj)*t + InputData%PhiCorps(jj))
                    elseif(abs(InputData%Vcst(jj)).gt.Epsilon)then ! Constant velocity.                
                        Mesh%Body(nc)%ABody(j) = 0._RP
                    end if
                else ! Blocked dof.
                    Mesh%Body(nc)%ABody(j) = 0._RP
                endif
            end do
            
            ! Time-differentiation of the body boundary condition.
            Int_Body0 = Mesh%Body(nc)%IndBody(1)
            Int_Body1 = Mesh%Body(nc)%IndBody(3)
            do j = Int_Body0,Int_Body1
                
                ! GM x n or AM x n
                M = Mesh%Tnoeud(j)%Pnoeud
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    AM = M - FixPointPos
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    AM = M - Mesh%Body(nc)%CSolv(1:3)
                end if
                call Computation_vect_product(AM,Mesh%Tnoeud(j)%Normale,AMvN) ! AMvN = AM x n if iFixPoint, = GM x n if not.
                
                ! Time-differentiation of the body condition.
                Ecoulement%DDPhiDnDt(j)%perturbation =  dot_product(Acc(1:3), Mesh%Tnoeud(j)%Normale) + dot_product(matmul(S(1:3,1:3,jj),Acc(4:6)),AMvN) - &
        &           TimeRamp*Ecoulement%DDPhiDnDt(j)%incident + Ecoulement%DGradPhiSqDn(j)%perturbation
            end do
            
            jj = jj + 1
            
        end if
    end do
    
    ! Deallocating.
    deallocate(S)
    
end subroutine ForceBodyMotion




!Yohan : comme PlotForces mais avec le fichier où on écrit en argument

subroutine PlotForces2(Mesh, Ecoulement, t, nc, io)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement)                   :: Ecoulement                       ! Flow parameters.
 
    integer                             :: nc ! Indice du corps
    integer                             :: io ! Indice du fichier ou on écrit                           
    real(rp)                            :: t
    real(rp)                            :: F_i
    real(rp), dimension(3)              :: Force
    real(rp), dimension(3)              :: normale
    real(rp)                            :: aire
    real(rp)                            :: DPhiDt, GPhi2, TimeRamp
    real(rp), dimension(3)              :: GPhip, GPhi0
    
    integer :: i,j,k
    
    call Ramp_Management(t,TimeRamp)
    
    
    Force = 0.
    
    do i = Mesh%Body(nc)%IndBody(2), Mesh%Body(nc)%IndBody(4)
        
        normale = Mesh%TFacette(i)%Normale
        aire = Mesh%TFacette(i)%Aire
        

        do j = 1, 3
        
            k = Mesh%TFacette(i)%Tnoeud(j)

            GPhip = Ecoulement%GPhi(:,k)%perturbation
            GPhi0 = TimeRamp*Ecoulement%GPhi(:,k)%incident
            DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
            GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip) - dot_product(GPhip,GPhip)
                
    
            Force = Force +  aire*normale*(DPhiDt + 0.5_rp*GPhi2 + g*Mesh%TNoeud(k)%PNoeud(3))

        end do
        
    end do
    
    Force = Force*ro*(1._rp/3._rp) ! 1/3 car précédemment la valeur sur la facette correspond à la moyenne des valeurs aux sommets
    
    write(io, '(4E)') t, Force  
    
    
end subroutine PlotForces2




subroutine PlotForces(Mesh, Ecoulement, t, InputData)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement)                   :: Ecoulement                       ! Flow parameters.
    real(rp), intent(in)                :: t                                ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                        ! Input data.
    
    integer                             :: j,nc,jj                          ! Loop parameter.
    real(rp)                            :: DPhiDt, GPhi2, TimeRamp,Surface
    real(rp), dimension(3)              :: GPhip, GPhi0, Ds, M
    real(rp), dimension(3)              :: Vect_product_1
    real(rp), dimension(3)              :: Fhs,Mhs,Fhydro,Mhydro            ! Hydrodynamic and hydrostatic loads.
    real(rp), dimension(6)              :: Weight,PTO,Morison               ! Waight, PTO and Morison loads.
    real(rp), dimension(3)              :: AG, Vhoule                       ! AG and the velocity of the waves.
    
    ! This subroutine writes the weight, hydrostatic, hydrodynamic, PTO and Morison loads in the output file.
    
    ! Time ramp
    call Ramp_Management(t,TimeRamp)
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
    

    
        ! Initialization
        Fhs = 0._RP
        Mhs = 0._RP
        Fhydro = 0._RP
        Mhydro = 0._RP
        Weight = 0._RP
        PTO = 0._RP
        Morison = 0._RP
        
        ! Hydrostatic and hydrodynamic loads.
        if(Mesh%Body(nc)%Active)then
            
            do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                
                M = Mesh%Tnoeud(j)%Pnoeud
                GPhip = Ecoulement%GPhi(:,j)%perturbation
                GPhi0 = TimeRamp*Ecoulement%GPhi(:,j)%incident
                DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
                GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip) - dot_product(GPhip,GPhip)
                
                ! Control surface.
                Ds = Mesh%Tnoeud(j)%Normale*Mesh%Tnoeud(j)%Aire
                
                ! Hydrostatic force.
                Fhs  = Fhs + g*M(3)*Ds
                
                ! Moment des forces hydrostatiques, calculé au point de résolution de l'équation de mouvement CSolv.
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    AG = M(1:3) - FixPointPos
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    AG = M(1:3) - Mesh%Body(nc)%CSolv(1:3)
                end if
                call Computation_vect_product(M(1:3)-Mesh%Body(nc)%CSolv(1:3), Ds,Vect_product_1)
                
                Mhs = Mhs + g*M(3)*Vect_product_1
                
                ! Hydrodynamic Force.
                Fhydro  = Fhydro + (DPhiDt + 0.5_rp*GPhi2 + g*M(3))*Ds
                
                ! Moment des forces hydrodynamiques, calculé au point de résolution de l'équation de mouvement CSolv.
                Mhydro = Mhydro + (DPhiDt + 0.5_rp*GPhi2 + g*M(3))*Vect_product_1
                
            end do
            
            Fhs  = ro*Fhs
            Mhs = ro*Mhs
            Fhydro = ro*Fhydro
            Mhydro = ro*Mhydro
            
            ! Symmetry
            if(Symmetry)then
                Fhs  = 2._rp*Fhs
                Mhs = 2._rp*Mhs
                Fhydro = 2._rp*Fhydro
                Mhydro = 2._rp*Mhydro
                Fhs(2)  = 0._rp
                Mhs(1) = 0._rp
                Mhs(3) = 0._rp
                Fhydro(2) = 0._rp
                Mhydro(1) = 0._rp
                Mhydro(3)  = 0._RP
            endif
        end if
        
        Weight = 0._RP
        PTO = 0._RP
        Morison = 0._RP
        
        ! Weight
        Weight(1:3) = [0._rp, 0._rp,-Mesh%Body(nc)%Mass*g] ! Force
        if(iFixPoint)then
            ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
            AG = Mesh%Body(nc)%GBody(1:3) - FixPointPos
        else
            ! Motion equation solved at G (CSolv = (GBody,angles)).
            AG = 0._RP
        end if
        call Computation_vect_product(AG,Weight(1:3),Weight(4:6)) ! Moment
                
        ! PTO
        do j=1,3 ! Force
            if(InputData%dll_dir(j,jj))then
                PTO(j) = - InputData%MuPTO(jj)*Mesh%Body(nc)%VBody(j)
                PTO(j) = PTO(j) - InputData%Raideur(jj)*(Mesh%Body(nc)%CSolv(j) - InputData%Pressort(j,jj) - InputData%Lressort(jj))
            else    
                PTO(j) = 0._rp
            endif
        enddo
        
        do j=4,6 ! Moment
            if(InputData%dll_dir(j,jj))then
                PTO(j) = - InputData%MuPTO(jj)*Mesh%Body(nc)%VBody(j)
                PTO(j) = PTO(j) - InputData%Raideur(jj)*Mesh%Body(nc)%CSolv(j)
            else    
                PTO(j) = 0._rp
            endif
        enddo
            
        ! Morison
        if (.true.) then
            Surface = 2._RP/3._RP*PI*Mesh%Body(nc)%DimBody(1)**2
            if (symmetry) Surface = 0.5_RP*Surface
            Vhoule = HouleRF%C*[1._RP,0._RP,0._RP] !suivant direction de propagation de la houle
            Morison(1:3) = -0.5*ro*InputData%Cd_Morison(jj)*Surface*norm2(Mesh%Body(nc)%VBody(1:3)-Vhoule)*(Mesh%Body(nc)%VBody(1:3)-Vhoule)
            Morison(4:6) = 0._RP
        else
            Morison(1:3) = 0._RP
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                ds = Mesh%Tnoeud(j)%Aire*dot_product(Mesh%Tnoeud(j)%Normale,Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)
                Morison(1:3) = Morison(1:3) + ds*norm2(Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)*(Mesh%Body(nc)%VBody(1:3)-Ecoulement%GPhi(1:3,j)%Incident)
            end do
            Morison(1:3) = 0.5_RP*ro*InputData%Cd_Morison(jj)*Morison(1:3)
            Morison(4:6) = 0._RP
        end if
        
        write(ioLoads+nc, '(37E)')  t, Weight(1:6),Fhs(1:3),Mhs(1:3),Fhydro(1:3),Mhydro(1:3),PTO(1:6),Morison(1:6),Weight(1:3)+Fhydro(1:3)+PTO(1:3)+Morison(1:3),Weight(4:6)+Mhydro(1:3)+PTO(4:6)+Morison(4:6)
        
        jj = jj + 1
        
    end do
            
end subroutine PlotForces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       Body Motion                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BodyMotion(Mesh,Mesh0,h)
    
    !f2py integer*1, dimension(1000)    :: Mesh,Mesh0
    type(TMaillage),intent(inout)       :: Mesh,Mesh0   ! Mesh and the copy of it
    real(rp), intent(in)                :: h            ! Time step

    integer                             :: nc           ! Loop parameter
    
    ! This subroutine updates the position of the floater at the end of a RK4 step.
    
    do nc = 1,Mesh%NBody
        if (Mesh%Body(nc)%CMD(1)) then
            Mesh%Body(nc)%MBody = Mesh%Body(nc)%VBody*h
            Mesh%Body(nc)%CSolv = Mesh0%Body(nc)%CSolv + Mesh%Body(nc)%MBody
            
            ! GBody.
            if(not(iFixPoint))then
                ! Motion equation solved at G (CSolv = (GBody,angles)).
                Mesh%Body(nc)%GBody = Mesh%Body(nc)%CSolv(1:3)
            end if
            ! If iFixPoint: GBody is updated in update_position when the subroutine Remesh is called.
            
        else !  Tank
            Mesh%Body(nc)%MBody = 0._RP            
        end if
    end do
    
end subroutine BodyMotion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       Body Velocity                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BodyVelocity(Mesh0,Mesh,t,h,InputData)
    !f2py integer*1, dimension(1000)    :: Mesh,Mesh0
    type(TMaillage)                     :: Mesh,Mesh0   ! Mesh and the copy of it
    real(rp), intent(in)                :: t, h         ! Current time and time step
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! InputData
    
    integer                             :: j, nc,jj     ! Loop parameters
    
    ! This subroutine updates the velocity of the floater at the end of a RK4 step.
    
    ! Tank
    Mesh%Body(1)%VBody = 0._RP
    Mesh%Body(1)%ABody = 0._RP
    
    ! Bodies
    jj = 1
    do nc = Int_Body,Mesh%NBody
        if(InputData%Free_Body(jj))then ! Free motion
            Mesh%Body(nc)%VBody = Mesh0%Body(nc)%VBody + Mesh%Body(nc)%ABody*h
        else
            do j=1,6
                if(InputData%dll_dir(j,jj))then ! Forced motion
                    if(abs(InputData%ACorps(jj)).gt.Epsilon)then
                        Mesh%Body(nc)%VBody(j) = InputData%ACorps(jj)*InputData%WCorps(jj)*cos(InputData%WCorps(jj)*t + InputData%PhiCorps(jj))
                        Mesh%Body(nc)%ABody(j) = -InputData%ACorps(jj)*InputData%WCorps(jj)**2*sin(InputData%WCorps(jj)*t + InputData%PhiCorps(jj))
                    elseif(abs(InputData%Vcst(jj)).gt.Epsilon)then                            
                        Mesh%Body(nc)%VBody(j) = InputData%Vcst(jj)
                        Mesh%Body(nc)%ABody(j) = 0._RP
                    end if
                else ! Fixed
                    Mesh%Body(nc)%VBody(j) = 0._RP
                    Mesh%Body(nc)%ABody(j) = 0._RP
                endif
            end do
        end if 
        jj = jj + 1
    end do
    
end subroutine BodyVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                           Inertia                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Inertia(Mesh,InputData,ConstantInertia)
    ! Calcul de la matrice d'inertie, par rapport au centre de gravité du corps, 
    ! à l'aide des grandeurs géométriques du corps. Déplacée au centre de résolution dans RotationCoeff
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh             ! Mesh
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData        ! InputData
    logical,intent(in)                  :: ConstantInertia  ! Inertia is constant and already computes (read_Inertia).
    
    integer                             :: j,jj,nc,jjj      ! Loop parameters
    real(RP)                            :: dm               ! Elementary mass
    real(RP), dimension(3)              :: M                ! Point
    integer                             :: ierror           ! Error flag
    
        
    ! This subroutine computes the inertia matrix at the geometrical center.
    
    ! Zeroing the total mass matrix of the tank.
    Mesh%Body(1)%IBody = 0._RP
    
    ! Floaters
    jjj = 1
    do nc = Int_Body,Mesh%NBody
        
        if(Mesh%Body(nc)%Active)then
        
            if(InputData%is_Inertia_File(jjj))then
            
                ! Reading the Inertia file when creating the mesh, not during the temporal loop.
                if(ConstantInertia)then
                    ! Inertia file
                
                    ! Zeroing
                    Mesh%Body(nc)%IBody = 0._RP
                
                    ! Mass and Inertia
                    call read_Inertia(InputData%file_inertia(jjj),Mesh,nc,ierror)
                
                    ! Consistency of the input data
                    if(abs(Mesh%Body(nc)%IBody(1,1)-Mesh%Body(nc)%Mass).gt.Epsilon)then
                        print*,"The mass of the body given in *.geom is not the same as the mass given in the inertia file for the body: ",jjj
                        print*,"*.geom: ",Mesh%Body(nc)%Mass
                        print*,"Inertia file: ",Mesh%Body(nc)%IBody(1,1)
                    end if
                
                end if
            
            else
                ! Total mass matrix from Numerical computations.
            
                ! Zeroing
                Mesh%Body(nc)%IBody = 0._RP
            
                ! Mass
                Mesh%Body(nc)%IBody(1,1) = Mesh%Body(nc)%Mass
                Mesh%Body(nc)%IBody(2,2) = Mesh%Body(nc)%IBody(1,1)
                Mesh%Body(nc)%IBody(3,3) = Mesh%Body(nc)%IBody(1,1)
            
                ! Inertia
                dm = ro*Mesh%Body(nc)%Volume/Mesh%Body(nc)%Aire
                do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    M = Mesh%Tnoeud(j)%Pnoeud-Mesh%Body(nc)%GBody(1:3) ! With respect to the center of gravity
                    Mesh%Body(nc)%IBody(4,4) = Mesh%Body(nc)%IBody(4,4) + dm*Mesh%Tnoeud(j)%Aire*(M(2)**2+M(3)**2)
                    Mesh%Body(nc)%IBody(5,5) = Mesh%Body(nc)%IBody(5,5) + dm*Mesh%Tnoeud(j)%Aire*(M(1)**2+M(3)**2)
                    Mesh%Body(nc)%IBody(6,6) = Mesh%Body(nc)%IBody(6,6) + dm*Mesh%Tnoeud(j)%Aire*(M(2)**2+M(1)**2)
                    Mesh%Body(nc)%IBody(4,5) = Mesh%Body(nc)%IBody(4,5) - dm*Mesh%Tnoeud(j)%Aire*M(1)*M(2)
                    Mesh%Body(nc)%IBody(4,6) = Mesh%Body(nc)%IBody(4,6) - dm*Mesh%Tnoeud(j)%Aire*M(1)*M(3)
                    Mesh%Body(nc)%IBody(5,6) = Mesh%Body(nc)%IBody(5,6) - dm*Mesh%Tnoeud(j)%Aire*M(2)*M(3)
                end do
                Mesh%Body(nc)%IBody(5,4) = Mesh%Body(nc)%IBody(4,5)
                Mesh%Body(nc)%IBody(6,4) = Mesh%Body(nc)%IBody(4,6)
                Mesh%Body(nc)%IBody(6,5) = Mesh%Body(nc)%IBody(5,6)
            
                ! Symmetry
                if (symmetry) then ! The symmetry is done because Mesh%Body(nc)%Volume was not updated by Symmetry.
                    do j = 4,6
                        do jj = 4,6
                            Mesh%Body(nc)%IBody(j,jj) = Mesh%Body(nc)%IBody(j,jj)*2._RP
                        end do
                    end do
                end if
            
            end if
        
            ! Writing the Inertia output file.
            call Write_Inertia(Mesh,nc)
        
        end if
        
        jjj = jjj + 1
        
    end do

end subroutine Inertia

subroutine Finite_differences(ti,Mesh,Ecoulement,EcoulementDF,Phi,size_Phi,h,InputData)
    
    real(rp),intent(in)                         :: ti                        ! Current time
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(in)                 :: Mesh                     ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)            :: Ecoulement, EcoulementDF
    type(TEcoulement),intent(inout)             :: Ecoulement, EcoulementDF ! Flow parameters (phi, etc.) and the copy of them
    integer,intent(in)                          :: size_Phi                 ! Size of Phi
    real(rp),dimension(size_Phi),intent(inout)  :: Phi                      ! Phi for the finite difference method
    real(rp),intent(in)                         :: h                        ! Time step
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(in)            :: InputData                ! Input data
    
    integer                                     :: j,k                      ! Loop parameters
    real(rp)                                    :: TimeRamp                 ! Ramp parameter
    integer                                     :: Int_body0,Int_body1      ! Index of the floater nodes
    
    ! This subroutine solves the dynamics using a finite difference method on DPhiDt.
    
    Int_body0 = Mesh%Body(Int_Body)%IndBody(1)
    Int_Body1 = Mesh%Body(Int_Body)%IndBody(3)
    
    ! Time ramp   
    call Ramp_Management(ti,TimeRamp)
    
    ! Local Derivatives on the body
    call GradientBody(Mesh,EcoulementDF)
                
    ! Calcul de DPhiDt par différences finies
    k = 1
    do j = Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
        ! Material Derivative following the body's motion
        EcoulementDF%DpPhiDt(j)%Perturbation = (Ecoulement%Phi(j)%Perturbation - Phi(k))/(h*2._RP)
        
        ! Partial derivative
        EcoulementDF%DPhiDt(j)%Perturbation = EcoulementDF%DpPhiDt(j)%perturbation &
        & - TimeRamp*dot_product(Mesh%Tnoeud(j)%Velocity + [1.0_RP,0._RP,0._RP]*Forward_Velocity,EcoulementDF%GPhi(1:3,j)%Perturbation)
        
        k = k + 1
    end do
    
    call PlotForces(Mesh,EcoulementDF,ti,InputData)
    
    k = 1
    do j = Int_Body0,Int_Body1
        Phi(k) = EcoulementDF%Phi(j)%perturbation
        k = k + 1
    end do
    
    call CopyEcoulement(EcoulementDF, Ecoulement, Mesh%Nnoeud)
    
end subroutine Finite_differences

subroutine Index_Number_Nodes(Mesh,Nsl,NBody,Nsys,Nb,size_NBody)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh         ! Total mesh (domain and floater)
    integer,intent(inout)                           :: Nsys,Nb      ! Number of nodes in the mesh, in the mesh of the floaters and the number of degree of freedom
    integer, dimension(3),intent(out)               :: Nsl          ! Index of the nodes for the free surface
    integer,intent(inout)                           :: size_NBody   ! Size of the structure NBody
    integer, dimension(size_NBody,3),intent(inout)  :: NBody        ! Index of the nodes for the floater
    
    integer                                         :: nc           ! Loop parameter
    
    ! This subroutine computes the number of nodes for the whole system and for the floater, the number of degree of freedom, the index of the nodes for the free surface and the floater and the size of the linear system.
    
    ! Number of unknows in the BVP problem
    Nsys = Mesh%Nsys
    
    ! Initialization
    Nb = 0
    
    ! Number of degree of freedom and number of nodes in the mesh of the floater
    do nc = 1,Mesh%NBody
        ! Index of the nodes of the floater
        NBody(nc,1:3) = [Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3), Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1]
        if (Mesh%Body(nc)%CMD(1)) then
            Nb = Nb + NBody(nc,3)
        end if
    end do
        
    ! Index of the nodes for the free surface
    Nsl = [Mesh%FS%IndFS(1), Mesh%FS%IndFS(3), Mesh%FS%IndFS(3) - Mesh%FS%IndFS(1) + 1]
    
end subroutine Index_Number_Nodes

subroutine get_S(Mesh,S,size_S)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh     ! Total mesh (domain and floater)
    integer,intent(in)                              :: size_S   ! NBodies
    real(rp), dimension(3,3,size_S),intent(inout)   :: S        ! Transformation matrix between the Cardan frames and the inertial frame
    
    integer                                         :: j,nc,jj  ! Loop parameter
    real(rp), dimension(3)                          :: theta    ! Angular position of the floater
    real(rp), dimension(3,2)                        :: Trig2    ! Cos and sin of theta

    ! This subroutine computes the matrices S (transformation matrix between the Cardan frames and the inertial frame).
    
    ! Initialization
    S = 0._RP
    
    ! Computation
    jj = 1
    do nc = Int_Body,Mesh%NBody
        
        ! Angular position
        theta(1:3) = Mesh%Body(nc)%Csolv(4:6)
    
        ! Cos and sin 
        do j=1,3
            Trig2(j,1) = cos(theta(j))
            Trig2(j,2) = sin(theta(j))
        end do
    
        ! S
        S(1:3,1:3,jj) = reshape([Trig2(2,1)*Trig2(3,1),Trig2(2,1)*Trig2(3,2),-Trig2(2,2),&
        -Trig2(3,2),Trig2(3,1),0._Rp,0._Rp,0._Rp,1._Rp],(/3,3/))
        
        jj = jj + 1
    end do
    
end subroutine get_S

subroutine get_DSDt(Mesh,DSDt,size_DSDt)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh             ! Total mesh (domain and floater).
    integer,intent(in)                              :: size_DSDt        ! NBodies.
    real(rp), dimension(3,3,size_DSDt),intent(out)  :: DSDt             ! Transformation matrix between the Cardan frames and the inertial frame.
    
    integer                                         :: nc,jj            ! Loop parameters.
    real(rp)                                        :: phi,theta,psi    ! Anglular position.
    real(rp)                                        :: phip,thetap,psip ! Anglular velocity.
    
    ! This subroutine computes the time-differentiation of the transformation matrices S.
    
    ! Initialization
    DSDt = 0._RP
    
    ! Computation
    jj = 1
    do nc = Int_Body,Mesh%NBody
        
        ! Cardan angles
        phi = Mesh%Body(nc)%CSolv(4)
        theta = Mesh%Body(nc)%CSolv(5)
        psi = Mesh%Body(nc)%CSolv(6)
        
        ! Angular velocities
        phip = Mesh%Body(nc)%VBody(4)
        thetap = Mesh%Body(nc)%VBody(5)
        psip = Mesh%Body(nc)%VBody(6)
        
        ! DSDt        
        DSDt(1,1,jj) = -thetap*sin(theta)*cos(psi) - psip*cos(theta)*sin(psi)
        DSDt(1,2,jj) = -psip*cos(psi)
        DSDt(1,3,jj) = 0._RP
        
        DSDt(2,1,jj) = -thetap*sin(theta)*sin(psi) + psip*cos(theta)*cos(psi)
        DSDt(2,2,jj) = -psip*sin(psi)
        DSDt(2,3,jj) = 0._RP
        
        DSDt(3,1,jj) = -thetap*cos(theta)
        DSDt(3,2,jj) = 0._RP
        DSDt(3,3,jj) = 0._RP
        
        jj = jj + 1
                        
    end do

end subroutine get_DSDt

subroutine get_CK(Mesh)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(inout)      :: Mesh     ! Total mesh (domain and floater).
    
    integer                             :: j        ! Loop parameter.
    integer                             :: nc       ! Loop parameter.
    integer                             :: k,p      ! Loop parameters.
    real(rp),dimension(:,:),allocatable :: Table    ! Table to save the transpose of CK(:,4:6) in order to the memory correctly.

    ! This subroutine computes the coefficients CK.
    
    ! Computation.
    do nc = 1,Mesh%NBody ! j = 0 should be after this loop !!! (PYW)
        if (Mesh%Body(nc)%CMD(1))then
            
            ! Allocating.
            allocate(Table(3,Mesh%Body(nc)%IndBody(3)-Mesh%Body(nc)%IndBody(1)+1))
            
            ! Translation.
            do p = 1,3
                j = 0
                do k = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    j = j + 1 !  Because k does not start at 1.
                    Mesh%Body(nc)%CK(j,p) = Mesh%Tnoeud(k)%Normale(p)
                end do
            end do
            
            ! Table (use to respect the column major in Fortran).
            j = 0
            do k = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                j = j + 1 ! Because k does not start at 1.
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    call Computation_vect_product(Mesh%Tnoeud(k)%Pnoeud - FixPointPos, Mesh%Tnoeud(k)%Normale(1:3),Table(1:3,j))
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    call Computation_vect_product(Mesh%Tnoeud(k)%Pnoeud - Mesh%Body(nc)%CSolv(1:3), Mesh%Tnoeud(k)%Normale(1:3),Table(1:3,j))
                end if
            end do
            
            ! Rotation.
            do p = 4,6
                j = 0
                do k = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    j = j + 1 ! Because k does not start at 1.
                    Mesh%Body(nc)%CK(j,p) = Table(p-3,j)
                end do
            end do
                        
            ! Deallocating.
            deallocate(Table)
            
        end if
    end do
    
end subroutine get_CK

subroutine get_CT(Mesh,S,size_S)
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage), intent(inout)              :: Mesh     ! Total mesh (domain and floater).
    integer,intent(in)                          :: size_S   ! NBodies.
    real(rp),dimension(3,3,size_S),intent(in)   :: S        ! Transformation matrix between the Cardan frames and the inertial frame.
    
    integer                                     :: j,jj     ! Loop parameters.
    integer                                     :: nc       ! Loop parameter.
    integer                                     :: k        ! Loop parameter.
    
    ! This subroutine computes the coefficients CT and updated CK.
    
    ! Computation
    jj = 1
    do nc=1,Mesh%NBody ! j = 0 should be after this loop !!! (PYW)
        if (Mesh%Body(nc)%CMD(1)) then
            j = 0
            do k=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                j = j+1 !  Because k does not start at 1.
                Mesh%Body(nc)%CT(1:6,j) = Mesh%Body(nc)%CK(j,1:6)*Mesh%Tnoeud(k)%Aire*ro
            end do
            
            ! In case of symmetry, the hydrodynamic loads must be counted twice.
            if(Symmetry)then
                Mesh%Body(nc)%CT = 2.*Mesh%Body(nc)%CT
                Mesh%Body(nc)%CT(2,:) = 0._rp ! Because the symmetry is around the x axis.
            endif
            
            ! Updating the CK coefficients (done here because previous CK coefficients were used to defined CT coefficients.
            Mesh%Body(nc)%CK(1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1,4:6) = matmul(Mesh%Body(nc)%CK(1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1,4:6),S(1:3,1:3,jj))
            
            jj = jj + 1
        
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
            
        end if
    end do
    
end subroutine get_CT

subroutine get_Inertia_terms(Mesh,M,S,DSDt,InertiaTerms,InputData,size_S,Tob)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                                 ! Total mesh (domain and floater)
    integer,intent(in)                              :: size_S                               ! NBodies
    real(rp), dimension(6,6,size_S),intent(inout)   :: M                                    ! Mass matrix used in the linear system
    real(rp),dimension(3,3,size_S),intent(in)       :: S                                    ! Transformation matrix between the Cardan frames and the inertial frame
    real(rp), dimension(3,3,size_S),intent(in)      :: DSDt                                 ! Transformation matrix between the Cardan frames and the inertial frame
    real(rp), dimension(3,size_S), intent(out)      :: Inertiaterms                         ! Inertia terms in the dynamical momentum expression
    real(rp),dimension(3,3,2),intent(out)           :: Tob                                  ! eR0
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData                            ! Input data
    
    real(rp), dimension(3,2)                        :: Trig2                                ! Cos and sin of theta
    real(rp), dimension(3)                          :: Omega                                ! Angular velocity in the inertia frame, Omega = S*Theta
    real(rp), dimension(3)                          :: AG                                   ! AG
    real(rp), dimension(3)                          :: TmpVect,TmpTmpVect,Vect_product      ! m*vect(AG,vect(omega,vect(omega,AG)))
    real(rp), dimension(6)                          :: MBody, VBody                         ! Position and velocity of the floater
    real(rp), dimension(3,3,2)                      :: DTobDt                           ! Transformation matrix and its time-differentiation
    real(rp), dimension(3,3)                        :: Topsi, Tpsitheta, Tthetab            ! Elementary transformation matrices
    real(rp), dimension(3,3)                        :: DTopsiDt, DTpsithetaDt, DTthetabDt   ! Time-differentiation of the elementary transformation matrices
    real(rp), dimension(3,3)                        :: IBody,Io, DIoDt,IBodyA               ! Sevral inertia matrices
    real(rp),dimension(3,3)                         :: A                                    ! A = S(AG)
    real(rp),dimension(3,3)                         :: SOmega                               ! OmegaMat = S(omega)
    integer                                         :: j,nc,jj                              ! Loop parameters
    
    ! This subroutine computes the mass matrix and the inertia terms of the dynamical momentum expression.
    
    M = 0._RP
    Inertiaterms = 0._RP
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
        
        ! Mass.
        M(1:3,1:3,jj) = Mesh%Body(nc)%IBody(1:3,1:3)
        
        ! Linear and angular position.
        MBody(1:6) = Mesh%Body(nc)%Csolv(1:6)
        
        ! Linear and angular velocity.
        VBody(1:6) = Mesh%Body(nc)%VBody(1:6)
        
        ! Angular velocity in the inertia frame.
        Omega = matmul(S(1:3,1:3,jj), VBody(4:6))
        
        ! Inertia in the body frame.
        IBody(1:3,1:3) = Mesh%Body(nc)%IBody(4:6,4:6)
                
        ! AG.
        if(iFixPoint .and. not(InputData%is_IBodyA(jj))) then
            AG(1:3) = Mesh%Body(nc)%GBody(1:3) - FixPointPos
        else
            AG(1:3) = 0._Rp
        endif
        
        ! Matrix A such as A*x = vect(AG,x).
        A = reshape([0._Rp, AG(3), -AG(2), -AG(3),0._Rp, AG(1), AG(2), -AG(1), 0._Rp],[3,3])
        
        ! Matrix SOmega such as SOmega*X = vect(Omega,X).
        SOmega = reshape([0._Rp, Omega(3), -Omega(2), -Omega(3),0._Rp, Omega(1), Omega(2), -Omega(1), 0._Rp],[3,3])
        
        ! Cos and sin.
        do j = 1,3
            Trig2(j,1) = cos(MBody(j+3))
            Trig2(j,2) = sin(MBody(j+3))
        end do
        
        ! Rotational Matrices.
        Topsi = reshape([Trig2(3,1),Trig2(3,2),0._Rp,-Trig2(3,2),Trig2(3,1),0._Rp,0._Rp,0._Rp, 1._Rp],[3,3])
        Tpsitheta = reshape([Trig2(2,1),0._Rp,-Trig2(2,2),0._Rp,1._Rp,0._Rp,Trig2(2,2),0._Rp,Trig2(2,1)],[3,3])
        Tthetab = reshape([1._Rp, 0._Rp,0._Rp, 0._Rp,Trig2(1,1),Trig2(1,2),0._Rp,-Trig2(1,2),Trig2(1,1)],[3,3])
        Tob(1:3,1:3,1) = matmul(Topsi(1:3,1:3),matmul(Tpsitheta(1:3,1:3),Tthetab(1:3,1:3)))
        Tob(1:3,1:3,2) = transpose(Tob(1:3,1:3,1))
        
        ! Inertia matrix in the inertia frame.
        Io(1:3,1:3) = matmul(Tob(1:3,1:3,1),matmul(IBody(1:3,1:3),Tob(1:3,1:3,2)))
                
        ! Inertia matrix in the inertia frame around the point of computation.
        IBodyA = Io - M(1,1,jj)*matmul(A,A)
        
        ! Time-differentiation of the rotational matrices.
        DTopsiDt = reshape([-VBody(6)*Trig2(3,2),VBody(6)*Trig2(3,1),0._Rp,-VBody(6)*Trig2(3,1),-VBody(6)*Trig2(3,2),0._Rp, 0._Rp,0._Rp, 0._Rp],[3,3])
        
        DTpsithetaDt = reshape([-VBody(5)*Trig2(2,2),0._Rp,-VBody(5)*Trig2(2,1),0._Rp,0._Rp,0._Rp,VBody(5)*Trig2(2,1),0._Rp,-VBody(5)*Trig2(2,2)],(/3,3/))
        
        DTthetabDt(1:3,1:3) = reshape([0._Rp,0._Rp,0._Rp,0._Rp,-VBody(4)*Trig2(1,2), VBody(4)*Trig2(1,1),0._Rp,-VBody(4)*Trig2(1,1), -VBody(4)*Trig2(1,2)],(/3,3/))
        
        DTobDt(1:3,1:3,1) = matmul(matmul(DTopsiDt(1:3,1:3),Tpsitheta(1:3,1:3)),Tthetab(1:3,1:3)) + &
                          matmul(matmul(Topsi(1:3,1:3),DTpsithetaDt(1:3,1:3)),Tthetab(1:3,1:3)) + &
                          matmul(matmul(Topsi(1:3,1:3),Tpsitheta(1:3,1:3)),DTthetabDt(1:3,1:3))
        DTobDt(1:3,1:3,2) = transpose(DTobDt(1:3,1:3,1))
        
        ! Inertia terms.
        if(iFixPoint)then
            
            ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
            
            ! TmpVect = m*vect(AG,vect(omega,vect(omega,AG)))
            call Computation_vect_product(Omega, AG,TmpTmpVect)
            call Computation_vect_product(Omega, TmpTmpVect,TmpVect)
            call Computation_vect_product(AG,TmpVect,Vect_product)
            TmpVect = M(1,1,jj)*Vect_product
        
            ! DIoDt
            DIoDt = matmul(matmul(DTobDt(1:3,1:3,1),IBody),Tob(1:3,1:3,2))
            DIoDt = DIoDt + matmul(matmul(Tob(1:3,1:3,1),IBody),DTobDt(1:3,1:3,2))
            
            InertiaTerms(1:3,jj) = matmul( matmul(IBodyA,DSDt(1:3,1:3,jj)) + matmul(DIoDt,S(1:3,1:3,jj)) ,VBody(4:6)) + TmpVect ! CC ! DIoDt = 0 (PYW)
            !InertiaTerms(1:3,jj) = matmul(matmul(IBodyA,DSDt(1:3,1:3,jj)),VBody(4:6)) + matmul(SOmega,matmul(Io,Omega)) + TmpVect ! PYW
        else
            
            ! Motion equation solved at G (CSolv = (GBody,angles)).
            Vect_product = 0._RP
            call Computation_vect_product(Omega, matmul(IBodyA,Omega),Vect_product) ! Vect_product = Omega x (IBodyA.Omega)
            InertiaTerms(1:3,jj) = matmul(matmul(IBodyA,DSDt(1:3,1:3,jj)),VBody(4:6)) + Vect_product
            
        end if
                
        ! Inertia
        M(4:6,4:6,jj) = matmul(IBodyA,S(1:3,1:3,jj)) ! Pay attention: I_A*S
        
        jj = jj + 1
        
    end do
    
end subroutine get_Inertia_terms

subroutine get_Q(Mesh,Ecoulement,TimeRamp)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(inout)      :: Mesh         ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(in)        :: Ecoulement   ! Flow parameters (phi, etc.)
    real(rp),intent(in)                 :: TimeRamp     ! Ramp parameter
        
    integer                             :: j,k,nc       ! Loop parameter
    
    ! This subroutine computes the Q coefficients.
    
    ! Computation
    do nc=1,Mesh%NBody
        if (Mesh%Body(nc)%CMD(1)) then
            j = 0
            do k=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                j = j + 1
                Mesh%Body(nc)%Q(j) =  - TimeRamp*Ecoulement%DDPhiDnDt(k)%incident + Ecoulement%DGradPhiSqDn(k)%perturbation
            end do
        end if
    end do    
    
end subroutine get_Q

subroutine get_Th(Mesh,Ecoulement,Fext,InertiaTerms,TimeRamp,InputData,size_Fext)

    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(inout)                  :: Mesh             ! Total mesh (domain and floater).
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement       ! Flow parameters (phi, etc.).
    integer,intent(in)                              :: size_Fext        ! NBodies.
    real(rp), dimension(6,size_Fext), intent(in)    :: Fext             ! External forces.
    real(rp), dimension(3,size_Fext), intent(in)    :: Inertiaterms     ! Inertia terms in the dynamical momentum expression.
    real(rp),intent(in)                             :: TimeRamp         ! Ramp parameter.
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData        ! InputData.
    
    integer                                         :: j,k,nc,idir,jj   ! Loop parameter.
    real(rp),allocatable                            :: SM(:)            ! Gradient and hydrostatics before the integration in the Bernoulli equation.
    real(rp)                                        :: GPhi2            ! Square of the grandient.
    real(rp), dimension(3)                          :: GPhi0            ! Gradient of the incident flow.
    


    ! This subroutine computes the coefficients Th.

    ! Computation
    jj = 1
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then
            j = 0
            
            Mesh%Body(nc)%Th = 0._RP
            
            ! Allocation
            allocate(SM(Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1))
            SM = 0._RP

            
            do k = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                j = j + 1
                
                ! Computation of the gradient term of the Bernoulli equation
                GPhi0 = TimeRamp*Ecoulement%GPhi(:,k)%incident
                GPhi2 = 0.5_rp*dot_product(GPhi0,GPhi0) + dot_product(GPhi0,Ecoulement%GPhi(:,k)%perturbation)
                SM(j) = TimeRamp*Ecoulement%DPhiDt(k)%incident + GPhi2
                                
                ! Hydrostatics
                if (hydrostatique) then ! Probably a problem with the hydrostatics (PYW)

                    if (lineaireBody .or. forcage_hydrostatique_lineaire) then 
                        SM(j) = SM(j) + g*Mesh%Body(nc)%CSolv(3)
                    else
                        SM(j) = SM(j) + g*Mesh%Tnoeud(k)%Pnoeud(3)
                    end if                 
                end if
            end do
            
            
            k = 0
            do idir = 1,6
                if(InputData%dll_dir(idir,jj))then
                    k = k + 1
                    ! Hydrodynamics
                    Mesh%Body(nc)%Th(k) = dot_product(Mesh%Body(nc)%CT(idir,1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1),SM(1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1))

                    
                    ! External forces
                    Mesh%Body(nc)%Th(k) = Mesh%Body(nc)%Th(k) + Fext(idir,jj)
                                    
                    
                    ! Inertia terms for the rotation
                    if(idir.gt.3)then
                        Mesh%Body(nc)%Th(k) = Mesh%Body(nc)%Th(k) - InertiaTerms(idir-3,jj)
                    end if
                                        
                end if
            enddo
            

             
            ! Deallocating
            deallocate(SM)

            
            jj = jj + 1
            
        ! The slip condition on DDPhiDnDt is obvious if the floater does not move (or tank)
        else 
            
            if(Mesh%Body(nc)%Active)then
                do k = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    Ecoulement%DDPhiDnDt(k)%perturbation = 0 
                end do
            end if
            
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
            
        end if
        
    end do
    
end subroutine get_Th

subroutine Building_A(Mesh,A,CS,CD,M,Nsl,NBody,size_NBody,N,Nsys,Nb,Nnodes,InputData,size_M)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                         ! Total mesh (domain and floater).
    real(rp),dimension(N,N),intent(inout)           :: A                            ! Matrix A of the linear system.
    real(rp), dimension(Nnodes,Nnodes),intent(in)   :: CD, CS                       ! Influence coefficients.
    integer,intent(in)                              :: size_M                       ! NBodies.
    real(rp), dimension(6,6,size_M),intent(in)      :: M                            ! Mass matrix used in the linear system.
    integer,intent(in)                              :: N,Nsys,Nb,Nnodes             ! Number of nodes in the linear system, unknowns in the mesh, in the mesh of the floaters and the number of degree of freedom, number of nodes in the mesh.
    integer, dimension(3),intent(in)                :: Nsl                          ! Index of the nodes for the free surface.
    integer,intent(in)                              :: size_NBody                   ! Size of the structure NBody.
    integer, dimension(size_NBody,3),intent(in)     :: NBody                        ! Index of the nodes for the floater.
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData                    ! InputData.
    
    integer                                         :: Nt,j,k,nc,idir,idir2,jj,Ndof ! Loop parameters.
    integer                                         :: p,q                          ! Loop parameters.
    
    ! This subroutine builds the matrix A of the linear system.
    
    ! Initialization
    A = 0._RP
    
    ! Integral equation
    ! CS(FS)
    A(1:Nsys,Nsl(1):Nsl(2)) = CS(1:Nsys,Nsl(1):Nsl(2))
    
    Nt = 0
    Ndof = 0
    jj = 1
    do nc = 1,Mesh%NBody
                
        if (not(cuve_ferme).and.Mesh%Body(nc)%is_tank) cycle ! No tank so Body(1) is not the tank.
        
        ! -CD(Ext) CD(B1) ... CD(Bn)
        if(Mesh%Body(nc)%Active)then
            A(1:Nsys,Nbody(nc,1):NBody(nc,2)) = -CD(1:Nsys,Nbody(nc,1):NBody(nc,2)) ! With the tank
        end if
        
        if (Mesh%Body(nc)%CMD(1)) then ! Not the tank.
            
            ! CS(B1) ... CS(Bn)
            A(1:Nsys,Nsys+Nt+1:Nsys+Nt+NBody(nc,3)) = CS(1:Nsys,Nbody(nc,1):NBody(nc,2))
            
            ! Motion equation (mass matrix)
            j = 0
            do idir = 1,6
                if (InputData%dll_dir(idir,jj)) then
                    j = j + 1
                    k = 0
                    do idir2 = 1,6
                        if (InputData%dll_dir(idir2,jj)) then
                            k = k+1
                            ! M
                            A(Nsys+Nb+Ndof+j,Nsys+Nb+Ndof+k) = M(idir,idir2,jj)
                        end if
                    end do
                end if
            end do
            
            Nt = Nt + NBody(nc,3)
            Ndof = Ndof + InputData%ndll(jj)
            jj = jj + 1
            
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
    
    ! Motion equation
    jj = 1
    Ndof = 0
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then
            ! This algorithm with lots of loop parameters is made to respect the column major in Fortran.
            q = 0
            do p = Nbody(nc,1),Nbody(nc,2)
                q = q + 1
                k = 0
                do idir = 1,6
                    if(InputData%dll_dir(idir,jj))then
                        k = k + 1
                        ! -CT0
                        A(Nsys+Nb+Ndof+k,p) = -Mesh%Body(nc)%CT(idir,q)
                    end if
                end do
            end do
            Ndof = Ndof + InputData%ndll(jj)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        endif
    enddo
    
    ! Slip condition
    jj = 1
    Nt = 0
    Ndof = 0
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then
            k = 0
            do idir = 1,6
                if(InputData%dll_dir(idir,jj)) then
                    k = k+1
                    ! -CK0
                    A(Nsys+Nt+1:Nsys+Nt+NBody(nc,3),Nsys+Nb+Ndof+k) = -Mesh%Body(nc)%CK(1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1,idir)
                endif
            end do
            Nt = Nt + NBody(nc,3)
            Ndof = Ndof + InputData%ndll(jj)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
    
    Nt = 0
    jj = 1
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then
            do j = Nsys+Nt+1,Nsys+Nt+Nbody(nc,3)
                ! I_Nb
                A(j,j) = 1._rp
            end do
            Nt = Nt + NBody(nc,3)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
    
    

end subroutine Building_A

subroutine Building_B(Mesh,Ecoulement,B,CD,Nsys,Nb,Nsl,Nnodes,N,InputData,NBody,size_NBody)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage),intent(in)                      :: Mesh                 ! Mesh
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(in)                    :: Ecoulement           ! Flow parameters (phi, etc.)
    integer,intent(in)                              :: Nnodes,N             ! Number of nodes in the mesh and number of unknowns in the linear system
    real(rp),dimension(N),intent(inout)             :: B                    ! Matrix B of the linear system
    real(rp), dimension(Nnodes,Nnodes),intent(in)   :: CD                   ! Influence coefficients
    integer,intent(in)                              :: Nsys,Nb              ! Number of nodes in the linear system, in the mesh of the floater and the number of degree of freedom
    integer, dimension(3),intent(in)                :: Nsl                  ! Index of the nodes for the free surface
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData            ! InputData
    integer,intent(in)                              :: size_NBody           ! Size of the structure NBody
    integer, dimension(size_NBody,3),intent(in)     :: NBody                ! Index of the nodes for the floater
    
    integer                                         :: j,k,nc,jj,Nt,Ndof    ! Loop parameters
    
    ! This subroutine builds the matrix B of the linear system.
    
    ! Initialization
    B = 0._RP
    

    ! Integral equation    
    do k = Nsl(1),Nsl(2)  
        do j = 1,Nsys
            B(j) = B(j) + CD(j,k)*Ecoulement%DPhiDt(k)%perturbation
        end do    
    end do
    
        
    ! Slip condition
    Nt = 0
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then
            B(Nsys+Nt+1:Nsys+Nt+NBody(nc,3)) = Mesh%Body(nc)%Q(1:Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1)
                        
            Nt = Nt + NBody(nc,3)
        end if
    end do
    
    ! Motion equation
    jj = 1
    Ndof = 0
    do nc = 1,Mesh%NBody
        if(Mesh%Body(nc)%CMD(1))then            
            B(Nsys+Nb+Ndof+1:Nsys+Nb+Ndof+InputData%ndll(jj)) = Mesh%Body(nc)%Th(1:InputData%ndll(jj))
            
            Ndof = Ndof + InputData%ndll(jj)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
    
end subroutine Building_B



subroutine Initialization_solution(Mesh,Ecoulement,Nsys,Nb,Nbody,size_NBody,Sol,N,InputData)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh                 ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement           ! Flow parameters (phi, etc.)
    integer,intent(in)                              :: Nsys,Nb              ! Number of nodes in the mesh and in the mesh of the floater
    integer,intent(in)                              :: size_NBody           ! Size of the structure NBody
    integer, dimension(size_NBody,3),intent(in)     :: NBody                ! Index of the nodes for the floater
    integer,intent(in)                              :: N                    ! Number of unknowns in the linear system
    real(rp), dimension(N), intent(out)             :: Sol                  ! Solution for the linear system
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData            ! InputData
    
    integer                                         :: nc,Nt,k,idir,jj,Ndof ! Loop parameters
    
    ! size_NBody = Mesh%NBody
    
    ! This subroutine initialized the solution of the GMRES solver.
    
    ! Phi(FS)_th    
    Sol(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3)) = Ecoulement%DDPhiDnDt(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Perturbation
    
    Nt = 0
    Ndof = 0
    jj = 1
    do nc = 1,Mesh%NBody
        
        ! Phi(B0)_t and Phi(Ext)_t
        if(Mesh%Body(nc)%Active)then
            Sol(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3)) = Ecoulement%DPhiDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation
        end if
        
        if (Mesh%Body(nc)%CMD(1)) then
            ! Phi(B0)_tn
            Sol(Nsys+Nt+1:Nsys+Nt+NBody(nc,3)) = Ecoulement%DDPhiDnDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation
            
            ! Acceleration
            k=0
            do idir=1,6
                if(InputData%dll_dir(idir,jj)) then
                    k=k+1
                    ! Acc
                    Sol(Nsys+Nb+Ndof+k) = Mesh%Body(nc)%ABody(idir)
                endif
            enddo
            Nt = Nt + NBody(nc,3)
            Ndof = Ndof + InputData%ndll(jj)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do

end subroutine Initialization_solution

subroutine Initialization_solution_MB(Mesh,Ecoulement,Nsys,Nbody,size_NBody,Sol,N)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh             ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement       ! Flow parameters (phi, etc.)
    integer,intent(in)                              :: Nsys             ! Number of nodes in the mesh and in the mesh of the floater, number of dof of the floater
    integer,intent(in)                              :: size_NBody       ! Size of the structure NBody
    integer, dimension(size_NBody,3),intent(in)     :: NBody            ! Index of the nodes for the floater
    integer,intent(in)                              :: N                ! Number of unknowns in the linear system
    real(rp), dimension(N), intent(out)             :: Sol              ! Solution for the linear system
    
    integer                                         :: nc,Nt,jj         ! Loop parameters
    
    ! size_NBody = Mesh%NBody
    
    ! This subroutine initialized the solution of the GMRES solver for the WSC components
        
    ! Phi(FS)_th    
    Sol(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3)) = Ecoulement%DDPhiDnDt(Mesh%FS%IndFS(1):Mesh%FS%IndFS(3))%Perturbation
    
    Nt = 0
    jj = 1
    do nc = 1,Mesh%NBody
        
        ! Phi(Bj)_t and Phi(Ext)_t
        if(Mesh%Body(nc)%Active)then
            Sol(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3)) = Ecoulement%DPhiDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation
        end if
        
        if (Mesh%Body(nc)%CMD(1)) then
            ! Phi(Bj)_tn
            Sol(Nsys+Nt+1:Nsys+Nt+NBody(nc,3)) = Ecoulement%DDPhiDnDt(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Perturbation
            
            Nt = Nt + NBody(nc,3)
            jj = jj + 1
        else
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do

end subroutine Initialization_solution_MB

subroutine Preconditioning(A,B,N)
    
    real(rp),dimension(N,N),intent(inout)   :: A        ! Matrix A of the linear system.
    real(rp),dimension(N),intent(inout)     :: B        ! Matrix B of the linear system.
    integer,intent(in)                      :: N        ! Number of nodes in the linear system.
    
    integer                                 :: j,p      ! Loop parameters.
    real(rp),dimension(:),allocatable       :: inv_diag ! Inversion of the diagonal coefficients.
    
    ! This subroutine preconditiones the linear system of the 2nd BVP.
    
    ! Allocating.
    allocate(inv_diag(N))
    inv_diag = 0._RP
    
    ! Back-up the inverse of the diagonal coefficients.
    do j = 1,N
	    if(abs(A(j,j)).GT.Epsilon)then
            inv_diag(j) = 1._RP/A(j,j)
            B(j) = B(j)*inv_diag(j)            
        else
		    print*, 'Preconditioning of the 2nd BEM problem: the coefficient ',j,' of the diagonal of A is nul!'
            pause
        end if
    end do
     
    ! Preconditioning of A.
    do j = 1,N
        do p = 1,N
		    A(p,j) = A(p,j)*inv_diag(p)
        end do
    end do
    
    ! Deallocating.
    deallocate(inv_diag)
    
end subroutine Preconditioning

subroutine Init_CK_CT_Q_Th(Mesh,InputData)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh         ! Mesh
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input Data
    
    integer                             :: nc,jj        ! Loop parameters
    
    ! This subroutine allocates and initializes CK, CT, Q and Th for each body.
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
        allocate(Mesh%Body(nc)%CK(Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1,6))
        Mesh%Body(nc)%CK = 0._RP
        
        allocate(Mesh%Body(nc)%CT(6,Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1))
        Mesh%Body(nc)%CT = 0._RP
        
        allocate(Mesh%Body(nc)%Q(Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1))
        Mesh%Body(nc)%Q = 0._RP
        
        allocate(Mesh%Body(nc)%Th(InputData%ndll(jj)))
        Mesh%Body(nc)%Th = 0._RP
        
        jj = jj + 1
    end do

end subroutine Init_CK_CT_Q_Th

subroutine deallocate_CK_CT_Q_Th(Mesh)

    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh     ! Mesh
        
    integer                             :: nc       ! Loop parameter
    
    ! This subroutine deallocates CK, CT, Q and Th for each body.
    
    do nc = Int_Body,Mesh%NBody
        deallocate(Mesh%Body(nc)%CK)
        deallocate(Mesh%Body(nc)%CT)
        deallocate(Mesh%Body(nc)%Q)
        deallocate(Mesh%Body(nc)%Th)
    end do
    
end subroutine deallocate_CK_CT_Q_Th

end module BodyMotion_mod
