module BoucleTemp
use Constantes
use Parameters
use Structuresdonnees
use FonctionsCommunes
use GenMaillage
use GeomMesh
use Incident_mod
use BVP
use Spline
use PrePlot
use GeomStruct
use SolvNum
use MeshModule
use BodyMotion_mod
use EnergyVolume
use rbf_interp
use GeomGen
use Remesh_mod
!use Test_Solver
implicit none

contains

subroutine BoucleTemporelle_RK4(Mesh, fgeom_vect, fdomaine, nface, Grid, nb_point,nb_tri, InputData,n_tab2,n_tab,get_State,fileState,jt0) 
    
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage), intent(inout)          :: Mesh                                         ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)        :: fgeom_vect
    type(type_GeomVect), intent(inout)      :: fgeom_vect                                   ! Geometry of the floaters (not the domain)
    !f2py integer*1, dimension(1000)        :: fdomaine
    type(type_geom),intent(in)              :: fdomaine                                     ! Geometry of the domain
    integer,intent(in)                      :: nface                                        ! Number of faces in both the floater and the domain geometries.
    !f2py integer*1, dimension(1000)        :: Grid
    type(MGrid),intent(inout)               :: Grid                                         ! Transitional mesh
    integer,intent(in)                      :: nb_point,nb_tri                              ! Number of points and triangles in the Mgrid
    !f2py integer*1, dimension(1000)        :: InputData
    type(InputDataStruct),intent(inout)     :: InputData                                    ! Input data
    integer,intent(inout)                   :: n_tab2,n_tab                                 ! Number of intersection curves and lines.
    logical,intent(in)                      :: get_State                                    ! True if the state input file is present, false otherwise.
    character(len=50),intent(in)            :: fileState                                    ! State input file.
    integer,intent(in)                      :: jt0                                          ! Initial time step parameter.
    
    logical                                 :: boolRemesh                                   ! Boolean to know if the mesh of the bodies or/and the free surface was regenerated (true) or not (false).
    logical                                 :: boolRemeshFS                                 ! = true: remeshing the free surface, = false: no remeshing.
    integer                                 :: jk, jt, nc, ierror,jj                        ! Parameters
    real(rp)                                :: ti, h, time_begin, rCI, tmoy                 ! Time parameters
    real(rp),allocatable                    :: RKABody(:,:,:), RKVBody(:,:,:)               ! RK4 data for the floater
    integer                                 :: Nnodes                                       ! Number of nodes in the mesh
    integer                                 :: Nnodes_FS                                    ! Number of nodes in the mesh of the free surface
    real(rp), dimension(nt)                 :: time_end                                     ! Final time of the RK4 loop
    real(rp), allocatable                   :: CD(:,:), CS(:,:), RK(:,:,:), RKVel(:,:,:)    ! Influence coefficients, RK4 data for the flow and the mesh
    real(rp),dimension(:),allocatable       :: Phi                                          ! Phi for the finite difference method
    !f2py integer*1, dimension(1000)        :: Mesh0
    type(TMaillage)                         :: Mesh0                                        ! Local copy of Mesh in the RK4 loop
    !f2py integer*1, dimension(1000)        :: Ecoulement, Ecoulement0, EcoulementDF
    type(TEcoulement)                       :: Ecoulement, Ecoulement0, EcoulementDF        ! Flow parameters (phi, etc.) and the copy of them
    real(rp), parameter                     :: denom = 1._RP/6._RP                          ! Denominateur of the RK4 mean
    real(rp)                                :: periode                                      ! Wave period
    real(rp), dimension(:),allocatable      :: t                                            ! Time vector
    integer                                 :: NumBodyNodes                                 ! Index of the floater nodes of the current floater mesh
    integer                                 :: size_Phi                                     ! Size of Phi
    integer,dimension(:,:),allocatable      :: IndBody_former                               ! (1st node, last node)xNbodies: Index of the floater nodes of the former mesh when new mesh generation.
    integer                                 :: nRemesh                                      ! Number of remeshing (use for the number of Advance_Front_Remesh.dat)
    logical                                 :: ForcedRemesh                                 ! Force the remeshing if necessary (crossing the free surface).
    logical									:: CrossingFS									! A body crossed the free surface (True) or not (False).
    !f2py integer*1, dimension(1000)        :: Mesh_State
    type(TMaillage)                         :: Mesh_State                                   ! Mesh read in the state input file (WARNING: only Tnoeud and Tfacette are (partially) filled).
    !f2py integer*1, dimension(1000)        :: Ecoulement_State
    type(TEcoulement)                       :: Ecoulement_State                             ! Flow parameters read in the state input file.
    integer                                 :: jt_init                                      ! First time step parameter.
    real(rp)                                :: Starting_time                                ! Starting time saved in case of calling a state input file.
    integer                                 :: jFiltering                                   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
    integer                                 :: ierrorRemesh                                 ! = 0 when the remeshing worked, 1 otherwise.
    
    real(rp),dimension(20) :: time
    real(rp),dimension(20) :: time_delta
    real(rp)  :: t_total
        
    ! This subroutine is the temporal loop. It updates the state vector at each time step.
    
    ! State vector = (Eta_p, Phi_p on the FS, Node_positions, floater_position, floater_velocity).
    
    ! Only a generation of mesh at t0.
    if(nt.lt.1)then
        go to 25
    end if
    
    ! Errors.
    ierror = 0
    ierrorRemesh = 1
    
    ! ForcedRemesh and CrossingFS.
    ForcedRemesh = .false.
    CrossingFS = .false.
    
    ! nRemesh
    nRemesh = 0
    
    ! jFiltering
    jFiltering = 0
    
    ! Nnodes and Nnodes_FS
    Nnodes = Mesh%Nnoeud
    Nnodes_FS = Mesh%FS%IndFS(3) - Mesh%FS%IndFS(1) + 1
    
    ! Allocation
    allocate(IndBody_former(2,NBodies))
    IndBody_former = 0._RP
    NumBodyNodes = 0
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
        if(Mesh%Body(nc)%Active)then
            IndBody_former(1,jj) = Mesh%Body(nc)%IndBody(1)
            IndBody_former(2,jj) = Mesh%Body(nc)%IndBody(3)
            NumBodyNodes = NumBodyNodes + IndBody_former(2,jj) - IndBody_former(1,jj) + 1
            jj = jj + 1
        end if
    end do
    
    allocate(t(nt))
    allocate(CS(Nnodes,Nnodes), CD(Nnodes,Nnodes))
    allocate(RK(Nnodes_FS,4,2))
	allocate(RKVBody(6,4,Nbodies),RKABody(6,4,Nbodies))
    if(DeformMesh) allocate(RKVel(3,Nnodes,4)) ! RKVel != RK
    if(FiniteDifference_ForceCalcul)then
        size_Phi = NumBodyNodes
        allocate(Phi(NumBodyNodes))
        Phi = 0._RP
    end if
    
    ! Time vector
    call Time_vector(t,nt,periode,InputData)
    
    ! Initialization of the output files
    call PrePlots(InputData,Mesh%NBody,get_State)
    
    ! Initialization of Mesh, Mesh0, Ecoulement, Ecoulement0 and eventually EcoulementDF for the updating at the end of each RK4 step
    call Initialization_Mesh_Ecoulement(Mesh,Mesh0,Ecoulement,Ecoulement0,t(1),EcoulementDF,InputData)
    
    ! Updating Ecoulement in case of state input file.
    if(get_State)then
        
        ! Reading the state input file.
        call read_State(fileState,Mesh_State,Ecoulement_State,Starting_time,jFiltering)
        
        ! Updating nliss in case of crossing the free surface.
        if(jFiltering.ne.0)then
            nliss = 1
		end if
        
        ! Interpolation of Phi_p and Eta_p from the Mesh_State.
        
        
        call Interpolation_FS(Mesh_State,Mesh,Ecoulement_State,t(1),.true.,ierror) ! True because there is a FS remeshing, therefore an interpolation.

        ! Initilization of Ecoulement from Ecoulement_State.
        call CopyEcoulement(Ecoulement, Ecoulement_State, Mesh%Nnoeud)
        
        ! Updating the velocity of the floaters.
        do nc = Int_Body,Mesh%NBody
            Mesh%Body(nc)%VBody = Mesh_State%Body(nc)%VBody
        end do
        
        ! Deleting.
        call DelEcoulement(Ecoulement_State)
        call DelMaillage(Mesh_State)
        
        ! Initialization of jt_init.
        jt_init = jt0
        
    else
        
        ! Initialization of jt_init.
        jt_init = 1
        
    end if
    
    ! First calculation of the influence coefficient to know how long the simulation will lasts.
    call Time_evaluation(Mesh,CD,CS,time_begin,time_end,Nnodes,nt,get_state,jt_init,Starting_time)
    
    time = 0._RP
    time_delta = 0._RP
    t_total = 0._RP
    
    ! Temporal loop
    do jt = jt_init,nt
        
        ! RK4 steps
        do jk = 1,4
            
            ! Computation of the time step of the RK4 step and the current time
            call Time_step_Current_time(t,jk,jt,ti,h)
            
            ! Influence coefficient parameter (rCI = -1 : Computation of CD and CS for this RK4 step, if rCI > 0 : no computation)
            call rCI_manager(jk,rCI)
#ifdef _OPENMP
    !$ time(1) = omp_get_wtime()
#else
    call cpu_time(time(1))
#endif
                        
            ! Application of a Gaussian filter on the free surface for Eta
            if(jk.eq.1 .and. jt.ne.1)then
                if(nliss.ne.0)then
                    if(mod(jt+1,nliss).eq.0)then !.and. t(jt).gt.T2
                        call LissageGaussien(Mesh, Ecoulement, ierror)
                    endif
                end if
            end if
#ifdef _OPENMP
    !$ time(2) = omp_get_wtime()
#else
    call cpu_time(time(2))
#endif
            
            !----------------------------------------------------------------------------------
            !           Computation of the time differentiation of the state vector
            !----------------------------------------------------------------------------------
            
            ! Computing the incident flow            
            call Incident(ti, Mesh, Ecoulement)
            
#ifdef _OPENMP
    !$ time(3) = omp_get_wtime()
#else
    call cpu_time(time(3))
#endif
            
            ! Test solveur BVP
            !call InitValidationSolveur(Mesh, Ecoulement)
            !if(jk.eq.1 .and. jt.ne.1)then
            !    call PlotBVP(Mesh, Ecoulement, ti,.true.) ! True = Test_BVP.
            !end if
            
            ! Writting the info in the output files.
            if(jk.eq.1)then
                if(nout.ne.0)then
                    if(mod(jt,nout).eq.0)then
                        call Plots(ti, Mesh, Ecoulement,InputData) ! Should be called after the four RK steps, but this definition matches with LL and CC. (PYW)
                    end if
                end if
            end if
#ifdef _OPENMP
    !$ time(4) = omp_get_wtime()
#else
    call cpu_time(time(4))
#endif

            ! Actually, the new time step starts here.
            
            ! New mesh if necessary.
            if(jk.eq.1 .and. jt.ne.1)then
                
                ! Generation of a new mesh. 
         
                call Regeneration_Mesh(Mesh,Ecoulement,ti,boolRemesh,boolRemeshFS,fgeom_vect,fdomaine,nface,Grid,nb_point,nb_tri,IndBody_former,NBodies,Nnodes,Nnodes_FS,ierror,InputData,nRemesh,n_tab2,n_tab,ForcedRemesh,CrossingFS)
                
                
                
                ! Updating nliss in case of crossing the free surface.
                call Updating_nliss(CrossingFS,jFiltering)
                
                ! Deallocation and new allocation.
                if(Mesh_type.eq.2)then ! Mesh strategy of CC.
                    
                    if(boolRemesh)then
                        
                        if(ierror == 0)then ! Remeshing was sucessful.
                            
                            ierrorRemesh = 0
                            
                            ! Deallocating
                            if(allocated(CS)) deallocate(CS)
                            if(allocated(CD)) deallocate(CD)
                            if(allocated(RK)) deallocate(RK)
                            if(allocated(RKVel)) deallocate(RKVel)
                            
                            allocate(CS(Nnodes,Nnodes), CD(Nnodes,Nnodes))
                            CD = 0._RP; CS = 0._RP
                            
                            ! New initialization of CD and CS if the remeshing was valid.
                            call CoeffInfl(Mesh, CD, CS, Nnodes)
                            !print*,"CoeffInfl is not called after Regeneration_Mesh."
                            
                            if (DeformMesh) allocate(RKVel(3,Nnodes,4))
                            allocate(RK(Nnodes_FS,4,2))
                            
                            ! State vector.
                            if(iwState)   call Write_State(Mesh,Ecoulement,ti,jt,time_begin,jFiltering) ! Only when the remeshing was sucessful to be able to regenerate the mesh.
                            
                        else
                            boolRemesh = .false.
                            boolRemeshFS = .false.
                            ierrorRemesh = 1
                        end if
                        
                    end if
                    
                else ! Mesh strategy of LL.
                    
                    ! State vector.
                    if(iwState)   call Write_State(Mesh,Ecoulement,ti,jt,time_begin,jFiltering) ! At each time step.
                    
                end if
                
                ! Zeroing the quantities for the new points and creates a new Mesh0.
                call Zeroing_new_points(Mesh,Mesh0,Ecoulement,Ecoulement0,boolRemesh,boolRemeshFS,IndBody_former,NBodies,ierror)
                
            elseif(jk.ne.1)then
                
                ! boolRemesh.
                boolRemesh = .false.
                
            end if
            
#ifdef _OPENMP
    !$ time(5) = omp_get_wtime()
#else
    call cpu_time(time(5))
#endif
                        
            ! Reinitialisation de Mesh0 et Ecoulement0 avant la resolution de la passe RK1.
            if (jk.eq.1) then
                call CopyMaillage(Mesh0, Mesh)
                call CopyEcoulement(Ecoulement0, Ecoulement, Mesh%Nnoeud)
            end if
            
#ifdef _OPENMP
    !$ time(6) = omp_get_wtime()
#else
    call cpu_time(time(6))
#endif
            ! Body condition
            call BodyCondition(ti, Mesh, Ecoulement)
            
#ifdef _OPENMP
    !$ time(7) = omp_get_wtime()
#else
    call cpu_time(time(7))
#endif
            
            ! Calcul des vitesses des noeuds du maillage (updating the velocity of the nodes)
            if (DeformMesh) call MeshVel(Mesh, ti, 1._RP,InputData) ! PYW: MeshVel uses a time step equal to 1.0 instead of h: validated with the floating cylinder test case.
            !if (DeformMesh .and. (ierror == 0)) call MeshVel(Mesh, ti, 1._RP,InputData) ! Deformation only when the remeshing worked (for horizontal cylinder case).
            !print*,"Condition on MeshVel modified!"
#ifdef _OPENMP
    !$ time(8) = omp_get_wtime()
#else
    call cpu_time(time(8))
#endif
            ! Computation of the influence coefficients, Phi and DPhiDn: boundary element solver (1st BVP).
            call solBVP(Ecoulement, Mesh, CD, CS, Nnodes,time,boolRemesh,rCI) ! rCI = -1 execept for frozen RK.
#ifdef _OPENMP
    !$ time(12) = omp_get_wtime()
#else
    call cpu_time(time(12))
#endif
            ! Computation of the velocity potential and FS elevation local derivatives on the FS and Bodies (potential only).
            call Gradient(Mesh, Ecoulement)
#ifdef _OPENMP
    !$ time(13) = omp_get_wtime()
#else
    call cpu_time(time(13))
#endif
            
            ! Computation of DPhiDt and DEtaDt and spatial differentiations
            call Derive(ti, Ecoulement, Mesh)
#ifdef _OPENMP
    !$ time(14) = omp_get_wtime()
#else
    call cpu_time(time(14))
#endif
            ! Computation of the Hydrodynamic loads and motions
            if(FreeBodies)then
                
                call FreeBodyMotion(Mesh, Ecoulement, CD, CS, ti,Nnodes,InputData,time)
                
            else ! not(FreeBodies)
                if(FiniteDifference_ForceCalcul .and. jt.ne.1)then
					
					print*,"No MB in these subroutines!"
                    
                    if(jk.eq.1)then ! Only for the first RK4 step
                        call Finite_differences(t(jt-1),Mesh,Ecoulement,EcoulementDF,Phi,size_Phi,h,InputData)
                    end if
                    
                else ! not(FiniteDifference_ForceCalcul)
                    
                    call ForceBodyMotion(Mesh, Ecoulement, ti,InputData)
                    
                    ! Computation of the DPhiDt and DDPhiDnDt: 2nd BVP.
                    call solBVP(Ecoulement, Mesh, CD, CS,Nnodes, time, .false., ti, .true.) ! rCI -> ti therefore, rCI > 0. Moreover Option = true, hence no new computation of the influence coefficients.
                                        
                end if
            end if
                  
#ifdef _OPENMP
    !$ time(17) = omp_get_wtime()
#else
    call cpu_time(time(17))
#endif
            
            ! Writing hydrodynamic loads in the output files
            if (jk.eq.1 .and. (FreeBodies .or. not(FiniteDifference_ForceCalcul)) .and. iwLoads) call PlotForces(Mesh,Ecoulement,ti,InputData)
            
#ifdef _OPENMP
    !$ time(18) = omp_get_wtime()
#else
    call cpu_time(time(18))
#endif
            
            ! Creation of the internal state vector for this RK step
            call RK_manager(Mesh,Ecoulement,RK,RKVBody,RKABody,RKVel,jk,Nnodes,Nnodes_FS,InputData,NBodies)
            
            !----------------------------------------------------------------------------------
            !               Preparation of the state vector for the next RK4 step
            !----------------------------------------------------------------------------------
            
#ifdef _OPENMP
    !$ time(19) = omp_get_wtime()
#else
    call cpu_time(time(19))
#endif
            if(jk.ne.4)then ! No updating for the 4th RK4 step.
                
                ! Computation of the time step of the RK4 step and the current time.
                call Time_step_Current_time(t,jk+1,jt,ti,h) ! jk+1 -> for the next RK4 step
                
                ! Position of the floater.
                call BodyMotion(Mesh,Mesh0,h) ! Could use RKVBody(1:6,jk)
                
                ! Mise a jour de la vitesse du corps au temps ti.
                call BodyVelocity(Mesh0, Mesh, ti, h,InputData) ! Could use RKABody(1:6,jk)
                
                ! Phi_p, Eta_p.
                call Free_surface_RK_step(Mesh,Ecoulement,Ecoulement0,RK,h,jk,Nnodes_FS)
                
                ! Mise a jour des geometries et de la position des noeuds du maillage a partir de la vitesse des noeuds au temps ti (deformation).

                if (DeformMesh) call Remesh(Mesh, Mesh0, ti,InputData, h, fgeom_vect) ! Could use RKVel.
                              
            end if
            
#ifdef _OPENMP
    !$ time(20) = omp_get_wtime()
#else
    call cpu_time(time(20))
#endif
            
            ! Time.
            do jj = 1,19
                t_total = t_total + time(jj+1) - time(jj)
                time_delta(jj) = time_delta(jj) + time(jj+1) - time(jj)
            end do
        
        end do ! End of RK4 stepss
        
        ! Time-stepping
        call Time_stepping(Mesh,Mesh0,Ecoulement,Ecoulement0,fgeom_vect,RK,RKVBody,RKABody,RKVel,denom,Nnodes,Nnodes_FS,ti,dt,InputData,NBodies)
        
        ! Writting the time info in the command window
        call Writting_time_info(t,time_end,time_begin,tmoy,jt,get_State,jt_init)
                
    end do ! End of temporal loop
    
    print*,""
    print*,"t_total = ",t_total,t_total/t_total*100
    print*,""
    print*,"Lissage = ",time_delta(1),time_delta(1)/t_total*100
    print*,"Incident = ",time_delta(2),time_delta(2)/t_total*100
    print*,"Plots = ",time_delta(3),time_delta(3)/t_total*100
    print*,"Remeshing = ",time_delta(4),time_delta(4)/t_total*100
    print*,"Mesh0/Ecoulement0 = ",time_delta(5),time_delta(5)/t_total*100
    print*,"BodyCondition = ",time_delta(6),time_delta(6)/t_total*100
    print*,"MeshVel = ",time_delta(7),time_delta(7)/t_total*100
    print*,"Influence coefficients = ",time_delta(8),time_delta(8)/t_total*100
    print*,"A, B, X - 1st BVP = ",time_delta(9),time_delta(9)/t_total*100
    print*,"Inversion - 1st BVP = ",time_delta(10),time_delta(10)/t_total*100
    print*,"Distribution - 1st BVP = ",time_delta(11),time_delta(11)/t_total*100
    print*,"Gradient = ",time_delta(12),time_delta(12)/t_total*100
    print*,"Derive = ",time_delta(13),time_delta(13)/t_total*100
    print*,"A, B, X - 2nd BVP = ",time_delta(14),time_delta(14)/t_total*100
    print*,"Inversion - 2nd BVP = ",time_delta(15),time_delta(15)/t_total*100
    print*,"Distribution - 2nd BVP = ",time_delta(16),time_delta(16)/t_total*100
    print*,"PlotForces = ",time_delta(17),time_delta(17)/t_total*100
    print*,"RK_manager = ",time_delta(18),time_delta(18)/t_total*100
    print*,"Remesh = ",time_delta(19),time_delta(19)/t_total*100
    print*,""    
    
    ! Last plot at t(nt) + dt, otherwise the computations after the last plots are never saved.
    call Plots(t(nt)+dt, Mesh, Ecoulement,InputData)
    
    ! State vector.
    if(iwState)   call Write_State(Mesh,Ecoulement,t(nt)+dt,jt,time_begin,jFiltering) ! At each time step.
    
    ! Closing of the output files
    call Closing(nt,time_end,time_begin,Mesh%NBody)
    
    ! Deallocating
    if(allocated(Phi))          deallocate(Phi)
    if(allocated(CD))           deallocate(CD)
    if(allocated(CS))           deallocate(CS)
    if(allocated(RK))           deallocate(RK)
    if(allocated(RKVel))        deallocate(RKVel)
    if(allocated(t))            deallocate(t)
    if(allocated(RKABody))      deallocate(RKABody)
    if(allocated(RKVBody))      deallocate(RKVBody)
    if(allocated(IndBody_former)) deallocate(IndBody_former)
    
    25 continue
    
end subroutine BoucleTemporelle_RK4

!-------------------------------------------------------------
!
!	Initialisation de l'écoulement pour l'intégration en temps
!
!--------------------------------------------------------------
!
subroutine BodyCondition(t, Mesh, Ecoulement)
    
    real(rp)                            :: t                        ! Current time.
    !f2py integer*1, dimension(1000)    :: Mesh 
    Type(Tmaillage)                     :: Mesh                     ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    Type(TEcoulement)                   :: Ecoulement               ! Flow parameters.
    
    integer                             :: j, k,nc,jj               ! Loop parameters.
    real(rp)                            :: TimeRamp                 ! Ramp.
    real(rp),dimension(3)               :: OmegaGM,vect_product_1 ! Cross product results.
    
    ! This subroutine computes the body condition on the body mesh.
    
    ! Ramp.
    call Ramp_Management(t,TimeRamp)
    
    ! Bodies.
    if(is_body)then
        jj = 1
        do nc = Int_Body,Mesh%NBody
            if(Mesh%Body(nc)%Active)then
                
                OmegaGM = 0._RP
                do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    
                    if(iFixPoint)then
                        ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-FixPointPos,OmegaGM) ! OmegaGM = VBody(4:6) x (OM - OA).
                        Ecoulement%DPhiDn(j)%perturbation = TimeRamp *(- Ecoulement%DPhiDn(j)%incident + dot_product(OmegaGM + [1.0_RP,0._RP,0._RP]*Forward_Velocity,Mesh%Tnoeud(j)%Normale))
                    else
                        ! Motion equation solved at G (CSolv = (GBody,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-Mesh%Body(nc)%CSolv(1:3),OmegaGM) ! OmegaGM = VBody(4:6) x (OM - OG).
                        Ecoulement%DPhiDn(j)%perturbation = TimeRamp *(- Ecoulement%DPhiDn(j)%incident + dot_product(Mesh%Body(nc)%VBody(1:3) + OmegaGM + [1.0_RP,0._RP,0._RP]*Forward_Velocity , Mesh%Tnoeud(j)%Normale))
                    end if
                    
                end do
                
            end if
            jj = jj + 1
        end do
    endif
    
    ! Tank.
    if (cuve_ferme) then
        do k = 1,Mesh%NBody
            if (k.lt.Int_body) then ! Tank
                do j = Mesh%Body(k)%IndBody(1),Mesh%Body(k)%IndBody(3)
                    Ecoulement%DPhiDn(j)%perturbation = 0._RP
                end do
            end if
        end do
    end if
    
end subroutine BodyCondition

!-------------------------------------------------------------
!
!	Dérivations spatiales de Phi et Eta
!   Dérivations temporelles de phi et Eta grâce aux CSL
!
!-------------------------------------------------------------

subroutine Derive(t, Ecoulement, Mesh)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh                                         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement),intent(inout)     :: Ecoulement                                   ! Flow parameters.
    real(rp),intent(in)                 :: t                                            ! Current time.
    
    integer                             :: j
    real(rp)                            :: Eta, Phi 
    real(rp), dimension(3)              :: VEl, M, GEta0, GEta, GPhi, GPhi0, DGPhi0Dz
    real(rp)                            :: D2Phi0DzDt, damp1, damp2, TimeRamp
    
    ! This subroutine computes the spatial and temporal differentiation of Phi and Eta in using the boundary conditions.
    
    ! Ramp.
    call Ramp_Management(t,TimeRamp)
    
    ! DEtaDt and DpPhiDt.
    if(lineaireFS) then ! Free surface mesh at z = 0
        do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            if (cuve_ferme) then
	            damp1 = 1._rp
	            damp2 = Mesh%Tnoeud(j)%Damping
            else
	            damp1 = Mesh%Tnoeud(j)%Damping
	            damp2 = 0._RP
            end if

            Gphi = Ecoulement%GPhi(:,j)%perturbation ! Ecoulement%GPhi(:,j)%incident = 0
            GEta = Ecoulement%GEta(:,j)%perturbation ! Ecoulement%GEta(:,j)%incident = 0
            Vel = Mesh%Tnoeud(j)%Velocity + [1.0_RP,0._RP,0._RP]*Forward_Velocity*TimeRamp

            Ecoulement%DPhiDt(j)%perturbation = -g*Ecoulement%Eta(j)%perturbation 
	        Ecoulement%DpPhiDt(j)%perturbation  = Ecoulement%DPhiDt(j)%perturbation + dot_product(Vel,GPhi)*(1._RP-damp2)
	        Ecoulement%DpPhiDt(j)%perturbation = damp1*Ecoulement%DpPhiDt(j)%perturbation - damp2*Ecoulement%Phi(j)%perturbation
			Ecoulement%DPhiDt(j)%perturbation  = damp1*Ecoulement%DPhiDt(j)%perturbation - damp2*Ecoulement%Phi(j)%perturbation
	
	        Ecoulement%DEtaDt(j)%perturbation = -Ecoulement%DPhiDn(j)%perturbation + dot_product(Vel,GEta)*(1._RP-damp2)
			Ecoulement%DEtaDt(j)%perturbation = damp1*Ecoulement%DEtaDt(j)%perturbation - damp2*Ecoulement%Eta(j)%perturbation 
    	end do

	else
    	do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)

	        if (cuve_ferme) then
		        damp1 = 1._rp
		        damp2 = Mesh%Tnoeud(j)%Damping
	        else
		        damp1 = Mesh%Tnoeud(j)%Damping
		        damp2 = 0._RP
	        end if

	        DGPhi0Dz = Ecoulement%DGPhiDz(:,j)%incident
			GPhi0 = Ecoulement%GPhi(:,j)%incident
			GEta0 = Ecoulement%GEta(:,j)%incident
			Gphi = Ecoulement%GPhi(:,j)%perturbation
			GEta = Ecoulement%GEta(:,j)%perturbation			
			Eta = Ecoulement%Eta(j)%perturbation
			Phi = Ecoulement%Phi(j)%perturbation
            M = Mesh%Tnoeud(j)%Pnoeud            
			call CD2Phi0DzDt(M,t,D2Phi0DzDt)
	        Vel = Mesh%Tnoeud(j)%Velocity + [1.0_RP,0._RP,0._RP]*Forward_Velocity*TimeRamp
			
			! Potential time derivative calculation
            ! Partial Derivative
            Ecoulement%DPhiDt(j)%perturbation = -g*Eta - dot_product(Gphi0,Gphi) - Eta*(dot_product(DGPhi0Dz,Gphi0) + D2Phi0DzDt)
            !Ecoulement%DPhiDt(j)%perturbation = - Ecoulement%DPhiDt(j)%incident - g*Ecoulement%Eta(j)%incident - 0.5_RP*dot_product(Gphi0,Gphi0) &
		            !                              -g*Eta - dot_product(Gphi0,Gphi) - Eta*(dot_product(DGPhi0Dz,Gphi0) + D2Phi0DzDt)
            
            !if (j == Mesh%FS%IndFS(1)) print*,'Phi :', Ecoulement%DPhiDt(j)%perturbation, -g*Eta, -dot_product(Gphi0,Gphi), - Eta*(dot_product(DGPhi0Dz,Gphi0) + D2Phi0DzDt)
            ! Material Derivative           
	        Ecoulement%DpPhiDt(j)%perturbation = Ecoulement%DPhiDt(j)%perturbation + dot_product(Vel,Gphi)
            
            ! Damping zone coefficients addition
	        Ecoulement%DPhiDt(j)%perturbation  = damp1*Ecoulement%DPhiDt(j)%perturbation - damp2*Phi
	        Ecoulement%DpPhiDt(j)%perturbation = damp1*Ecoulement%DpPhiDt(j)%perturbation - damp2*Phi 
            
            ! Wave elevation time derivative calculations
            ! Partial Derivative
	        Ecoulement%DEtaDt(j)%perturbation = Gphi(3) - dot_product(Gphi0,GEta) - dot_product(Gphi,Geta0) - Eta*(dot_product(DGPhi0Dz,GEta0)-DGphi0Dz(3))
            !Ecoulement%DEtaDt(j)%perturbation = - Ecoulement%DEtaDt(j)%incident + Gphi0(3) - dot_product(Gphi0,Geta0) &
            !                                    + Gphi(3) - dot_product(Gphi0,GEta) - dot_product(Gphi,Geta0) - Eta*(dot_product(DGPhi0Dz,GEta0)-DGphi0Dz(3))
            !if (j == Mesh%FS%IndFS(1)) print*, 'Eta :', Ecoulement%DEtaDt(j)%perturbation, + Gphi(3), - dot_product(Gphi0,GEta)- dot_product(Gphi,Geta0),- Eta*(dot_product(DGPhi0Dz,GEta0) - DGphi0Dz(3)),  + dot_product(Vel(1:2),GEta(1:2))
            
            ! Material Derivative 
            Ecoulement%DEtaDt(j)%perturbation = Ecoulement%DEtaDt(j)%perturbation + dot_product(Vel(1:2),GEta(1:2))
            
            ! Damping zone coefficients addition
	        Ecoulement%DEtaDt(j)%perturbation = damp1*Ecoulement%DEtaDt(j)%perturbation - damp2*Eta 
            
        end do
    end if

end subroutine Derive

subroutine Initialisation(Ecoulement, Mesh, t,InputData)
    
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement)                   :: Ecoulement                   ! Flow parameters.
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                         ! Mesh.
    real(rp), intent(in)                :: t                            ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                    ! Input data.
    
    character(len=50)                   :: fileflow                     ! Flow input file.
    integer                             :: iarg, ierror, j,nc,jj        ! Loop parameters.
    logical                             :: get_flow                     ! True if the flow input file is present, false otherwise.
    real(rp)                            :: tstart                       ! Starting time.
    real(rp),dimension(3,3)             :: Topsi, Tpsitheta, Tthetab    ! Elementary transformation matrices.
    real(rp),dimension(3,3)             :: Tob                          ! eR0.
    real(rp),dimension(3,2)             :: Trig2                        ! Cos and sin of theta.
    
    ! This subroutine initialized some physical parameters of Mesh.
    
    !call get_command_argument(4, fileflow, status=iarg) ! This fourth input is not a mesh anymore (State.txt).
    !if (iarg.eq.0) then
    !    get_flow = .true.
    !else
    !    get_flow = .false.
    !endif
    get_flow = .false. ! This command is totaly desactivated.
    

    
    if(get_flow)then
        call extract_ecoulement(Ecoulement, tstart, fileflow, ierror)
        print*,' Initialisation flow at t = ',tstart,' : ',fileflow
    else
        ! Initialisation des vitesses des points du maillage à 0
        do j=1,Mesh%Nnoeud
            Mesh%Tnoeud(j)%Velocity = 0._rp
        end do
        
        ! Initialisation des position, vitesse et accélération du corps
        jj = 1
        do nc = Int_Body,Mesh%NBody

            ! Positions.
            if(iFixPoint)then
                ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                Mesh%Body(nc)%GBody(1:3) = InputData%PositionG(1:3,jj)
                Mesh%Body(nc)%CSolv(1:3) = FixPointPos(1:3)
            else
                ! Motion equation solved at G (CSolv = (GBody,angles)). PositionG is given in the body frame -> shoule move to the inertial frame.
                
                ! Cos and sin.
                do j = 1,3
                    Trig2(j,1) = cos(Mesh%Body(nc)%CSolv(j+3))
                    Trig2(j,2) = sin(Mesh%Body(nc)%CSolv(j+3))
                end do
                    
                ! Rotational Matrices
                Topsi = reshape([Trig2(3,1),Trig2(3,2),0._Rp,-Trig2(3,2),Trig2(3,1),0._Rp,0._Rp,0._Rp, 1._Rp],[3,3])
                Tpsitheta = reshape([Trig2(2,1),0._Rp,-Trig2(2,2),0._Rp,1._Rp,0._Rp,Trig2(2,2),0._Rp,Trig2(2,1)],[3,3])
                Tthetab = reshape([1._Rp, 0._Rp,0._Rp, 0._Rp,Trig2(1,1),Trig2(1,2),0._Rp,-Trig2(1,2),Trig2(1,1)],[3,3])
                Tob(1:3,1:3) = matmul(Topsi(1:3,1:3),matmul(Tpsitheta(1:3,1:3),Tthetab(1:3,1:3)))
                
                Mesh%Body(nc)%GBody(1:3) = matmul(Tob(1:3,1:3),InputData%PositionG(1:3,jj))
                Mesh%Body(nc)%CSolv(1:3) = Mesh%Body(nc)%GBody(1:3) ! Motion equation is solved at G except in case of a fix point motion.
            end if
            Mesh%Body(nc)%CSolv(4:6) = InputData%Position(1:3,2,jj)
            Mesh%Body(nc)%DimBody(1) = InputData%Lgeom(2,jj)
            Mesh%Body(nc)%MBody = 0._rp
            
            ! Velocities and accelerations.
            if(InputData%free_body(jj))then
                Mesh%Body(nc)%VBody = 0._rp
                Mesh%Body(nc)%ABody = 0._rp
            else if(abs(InputData%ACorps(jj)).gt.Epsilon .or. abs(InputData%Vcst(jj)).gt.Epsilon)then ! Forced motion.
                do j=1,6
                    if(InputData%dll_dir(j,jj))then ! Forced dof.
                        if(abs(InputData%ACorps(jj)).gt.Epsilon)then ! Sinusoidal velocity.
                            Mesh%Body(nc)%VBody(j) = InputData%ACorps(jj)*InputData%WCorps(jj)*cos(InputData%WCorps(jj)*t0 + InputData%PhiCorps(jj)) ! xp(t) = A*w*cos(w*t + phi) therefore xp(0) = A*w*cos(w*t0 + phi).
                            Mesh%Body(nc)%ABody = -InputData%ACorps(jj)*InputData%WCorps(jj)**2*sin(InputData%WCorps(jj)*t0 + InputData%PhiCorps(jj)) ! xpp(t) = -A*w^2*sin(w*t + phi) therefore xpp(0) = -A*w^2*sin(w*t0 + phi)
                        elseif(abs(InputData%Vcst(jj)).gt.Epsilon)then ! Constant velocity.
                            Mesh%Body(nc)%VBody(j) = InputData%Vcst(jj)
                            Mesh%Body(nc)%ABody = 0._rp
                        end if
                    else ! Blocked dof.
                        Mesh%Body(nc)%VBody(j) = 0._RP
                        Mesh%Body(nc)%ABody = 0._rp
                    end if
                end do
            end if
            jj = jj + 1
            
        end do
        
        ! Repositionnement de la surface libre sur la position de la houle incidente
        if(not(lineaireFS) .and. Mesh_type.eq.1) then
            call Remesh(Mesh, Mesh, t, InputData,dt)
        endif
        
        ! Initialisation de l'écoulement à 0
        call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)
        
        ! Initialisation des composantes de houle incidente
        call Incident(t, Mesh, Ecoulement)
        
        ! Imposition de la condition limite sur le corps
        call BodyCondition(t, Mesh, Ecoulement)
        
        ! Inertia matrices
        call Inertia(Mesh,InputData,.true.)
        
    endif
    
end subroutine Initialisation

subroutine Writting_time_info(t,time_end,time_begin,tmoy,jt,get_State,jt_init)
    
    real(rp),intent(in),dimension(:)        :: t            ! Time vector
    real(rp),dimension(:),intent(inout)     :: time_end     ! Final time of the RK4 loop
    real(rp),intent(in)                     :: time_begin   ! Time parameters
    real(rp),intent(inout)                  :: tmoy         ! Mean time
    integer,intent(in)                      :: jt           ! Time loop parameter
    logical,intent(in)                      :: get_State    ! True if the state input file is present, false otherwise.
    integer,intent(in)                      :: jt_init      ! First time step parameter.
    
    ! This subroutine displays the time information at each time step on the command window.
    
    ! Information are written for the NEXT time step

    ! Fortran -> Preprocessor -> Preprocessor source file -> Yes
#ifdef _OPENMP
    !$ time_end(jt) = omp_get_wtime()
#else
    call cpu_time(time_end(jt))
#endif
        
    10 format(' Tcpu(',f8.4,') = ',f8.2,'s. Completion = ',i3,'%')
    20 format(' End in ',i3,' h, ',f4.1,' min')
           
    if ( (not(get_State) .and. jt.eq.1) .or. (get_State .and. jt.eq.jt_init) )then
        tmoy = 2*(time_end(jt) - time_begin)/3._RP
        write(*,10) t(jt+1), time_end(jt)- time_begin, int(dble(jt+1)/dble(nt)*100)
    else if(jt == nt)then
        go to 30 ! No time info.
    else
        tmoy = (tmoy*jt + time_end(jt) - time_end(jt-1))/(jt+1)
        write(*,10) t(jt+1), time_end(jt) - time_end(jt-1), int(dble(jt+1)/dble(nt)*100)
    end if
    
    write(*,20) int(tmoy*(size(t)-(jt+1))/3600), ( tmoy*(size(t)-(jt+1))/3600._RP - int(tmoy*(size(t)-(jt+1))/3600))*60._RP

    30    continue
    
end subroutine Writting_time_info

subroutine Closing(nt,time_end,time_begin,NBody)
    
    integer,intent(in)                  :: nt           ! Number of time steps.
    real(rp),dimension(:),intent(in)    :: time_end     ! Final time of the RK4 loop.
    real(rp),intent(in)                 :: time_begin   ! Time of the first time step.
    integer,intent(in)                  :: NBody        ! Mesh%NBody.
    
    ! This subroutine closes the output files.

    ! Closing of the output files.
    print*,'Temps total du calcul : ', time_end(nt)-time_begin, ' secondes'
    call close_output(time_end(nt)-time_begin,NBody)
    print*,""
    print*,"Alors, heureux ?"
    
end subroutine Closing

subroutine Time_evaluation(Mesh,CD,CS,time_begin,time_end,Nnodes,nt,get_State,jt_init,Starting_time)
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage),intent(in)                          :: Mesh         ! Total mesh (domain and floater)
    integer,intent(in)                                  :: Nnodes       ! Number of nodes in the mesh
    real(rp), dimension(Nnodes,Nnodes),intent(inout)    :: CD, CS       ! Influence coefficients
    real(rp),intent(inout)                              :: time_begin   ! Time of the first time step
    integer,intent(in)                                  :: nt           ! Number of time step
    real(rp),dimension(nt),intent(inout)                :: time_end     ! Final time of the RK4 loop
    logical,intent(in)                                  :: get_State    ! True if the state input file is present, false otherwise.
    integer,intent(in)                                  :: jt_init      ! First time step parameter.
    real(rp),intent(in)                                 :: Starting_time! Starting time saved in case of calling a state input file.
    
    ! This subroutine evaluates the time necessary to compute once the influence coefficients.
    
    if(get_State)then
        time_begin = Starting_time
    else
    
#ifdef _OPENMP
    !$ time_begin = omp_get_wtime()
#else
    call cpu_time(time_begin)
#endif
        
    end if
    
    ! Only time where CD and CS are initialized to 0 except after a remeshing.
    ! Initilization of CD and CS.
    CD = 0._RP ; CS = 0._RP
    call CoeffInfl(Mesh, CD, CS,Nnodes)      
    !print*,"CoeffInfl is not called in Time_evaluation."
    
    if(get_State)then
#ifdef _OPENMP
    !$ time_end(jt_init) = omp_get_wtime()
#else
    call cpu_time(time_end(jt_init))
#endif
        print*, 'Temps CI :', time_end(jt_init)-time_begin
    else
#ifdef _OPENMP
    !$ time_end(1) = omp_get_wtime()
#else
    call cpu_time(time_end(1))
#endif
        
        print*, 'Temps CI :', time_end(1)-time_begin        
    end if
	
#ifdef _OPENMP
    !$ time_begin = omp_get_wtime()
#else
    call cpu_time(time_begin)
#endif
            
end subroutine Time_evaluation

subroutine Initialization_Mesh_Ecoulement(Mesh,Mesh0,Ecoulement,Ecoulement0,t1,EcoulementDF,InputData)
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage),intent(in)                  :: Mesh                     ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)            :: Mesh0
    type(TMaillage),intent(out)                 :: Mesh0                    ! Local copy of Mesh in the RK4 loop
    !f2py integer*1, dimension(1000)            :: Ecoulement,Ecoulement0
    type(TEcoulement),intent(out)               :: Ecoulement,Ecoulement0   ! Flow parameters (phi, etc.) and the copy of them
    real(rp)                                    :: t1                       ! t(1).
    !f2py integer*1, dimension(1000)            :: EcoulementDF
    type(TEcoulement),optional,intent(out)      :: EcoulementDF             ! Flow parameters (phi, etc.) in case of finite differences
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(inout)         :: InputData                ! Input data
    
    ! This subroutines initilizes the copy of the mesh (Mesh0) and creates the flow parameters (Ecoulement) and a copy of it (Ecoulement0).
    
    ! Initialisation of the copy of the mesh
    call NewMaillage(Mesh0, Mesh%Nnoeud, NBodies+1) ! +1 for the tank.
    call CopyMaillage(Mesh0, Mesh)
    
    ! Initialisation of the copy of the flow parameters
    call NewEcoulement(Ecoulement, Mesh%Nnoeud)
    call NewEcoulement(Ecoulement0, Mesh%Nnoeud)
    if(FiniteDifference_ForceCalcul)then
        call NewEcoulement(EcoulementDF, PointMax)
        call IniEcoulement(EcoulementDF, Mesh%Nnoeud, 0._RP)
    endif
    call Initialisation(Ecoulement, Mesh,t1,InputData)
    call CopyEcoulement(Ecoulement0, Ecoulement, Mesh%Nnoeud)
    call CopyMaillage(Mesh0, Mesh)

end subroutine Initialization_Mesh_Ecoulement

subroutine Time_step_Current_time(t,jk,jt,ti,h)
    
    real(rp),dimension(:),intent(in)    :: t    ! Time vector
    integer,intent(in)                  :: jk   ! RK4 loop parameter
    integer,intent(in)                  :: jt   ! Temporal loop parameter
    real(rp),intent(out)                :: ti   ! Current time (of the RK4 step)
    real(rp),intent(out)                :: h    ! Time step
    
    ! This subroutine computes the time step h and the current time of the RK4 step.
    
    if(jk.eq.1)then
        h = 0.0_RP
    elseif(jk.eq.2)then
        h = 0.5_RP*dt
    elseif(jk.eq.3)then
        h = 0.5_RP*dt
    elseif(jk.eq.4)then
        h = dt
    end if
    
    ti = t(jt) + h

end subroutine Time_step_Current_time

subroutine Time_vector(t,nt,periode,InputData)
    
    integer,intent(in)                      :: nt           ! Number of time steps.
    real(rp),dimension(nt),intent(inout)    :: t            ! Time vector.
    !f2py integer*1, dimension(1000)        :: InputData
    type(InputDataStruct),intent(in)        :: InputData    ! Input data.
    
    real(rp)                                :: periode      ! Wave period.
    real(rp)                                :: lambda_wave  ! Wave length.
    integer                                 :: j            ! Parameter.

    ! This subroutine initializes the time vector.
    
    ! Free surface parameters.
    if (Htype.eq.0)then 
        if(is_body)then
            if(abs(InputData%wcorps(1)).gt.Epsilon)then ! Oscillating motion if there is a body.
                periode = 2._RP*Pi/InputData%wcorps(1)
                lambda_wave = 0._RP
                print*,"No wave length defined"
            else
                periode = 0._RP
                lambda_wave = 0._RP
                print*,"No wave length/period defined."
            end if
        else
            periode = 0._RP
            lambda_wave = 0._RP
            print*,"No wave length/period defined."
        end if
    end if
    if(Htype.eq.1)then 
        periode = 2._RP*Pi/w(1)
        lambda_wave = lambda(1)
    end if
    if(Htype.eq.2)then
        periode = HouleRF%T
        lambda_wave = 0._RP
    end if
    
    ! Time vector.
    t = (/((j-1)*dt+t0,j=1,nt)/)
    print*,'Wave period =', periode,"s"
    print*,"Wave length =",lambda_wave,"m"
    print*,'Time step = ', dt,"s"
    print*,'Simulation time = ', dt*nt,"s"

end subroutine Time_vector

subroutine rCI_manager(jk,rCI)

    integer,intent(in)  :: jk   ! RK4 loop parameter.
    real(rp),intent(out):: RCI  ! Influence coefficient parameter.

    ! This subroutines handles the influence coefficient parameter.

    ! Computation or not of the influence coefficient.
    rCI = -1._RP
    
    ! No influence coefficients calculation for RK4 step 2,3 and 4 (one computation per time step).
    if(Rk_fige .and. jk.ne.1) rCI = 1._RP ! jk.ne.1 : only the first step computes the influence coeffients.

end subroutine rCI_manager

subroutine RK_manager(Mesh,Ecoulement,RK,RKVBody,RKABody,RKVel,jk,Nnodes,Nnodes_FS,InputData,size_RKVBody)
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage), intent(in)                         :: Mesh                 ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                    :: Ecoulement
    type(TEcoulement),intent(in)                        :: Ecoulement           ! Flow parameters (phi, etc.)
    integer,intent(in)                                  :: Nnodes               ! Number of nodes in the mesh
    integer,intent(in)                                  :: Nnodes_FS            ! Number of nodes in the mesh of the free surface
    real(rp), dimension(Nnodes_FS,4,2),intent(inout)    :: RK                   ! RK4 data for the flow
    real(rp),dimension(3,Nnodes,4),intent(inout)        :: RKVel                ! RK4 data for the mesh
    integer,intent(in)                                  :: size_RKVBody         ! NBodies
    real(rp), dimension(6,4,size_RKVBody),intent(inout) :: RKABody, RKVBody     ! RK4 data for the floater
    integer,intent(in)                                  :: jk                   ! RK4 loop parameter
    !f2py integer*1, dimension(1000)                    :: InputData
    type(InputDataStruct),intent(in)                    :: InputData            ! InputData
    
    integer                                             :: Int_sl0,Int_sl1      ! Index of the free surface nodes
    integer                                             :: j,nc,jj              ! Parameter
        
    ! This subroutine handles the data for the time-stepping at the end of the four RK steps.
    
    ! Index of the floater and free surface nodes.
    Int_sl0 = Mesh%FS%IndFS(1)
    Int_sl1 = Mesh%FS%IndFS(3)
    
    ! Internal state vector: DpPhi_pdt and DEta_pdt (particular derivative) on the free surface.
    RK(Int_sl0:Int_sl1,jk,1) = Ecoulement%DpPhiDt(Int_sl0:Int_sl1)%perturbation
    RK(Int_sl0:Int_sl1,jk,2) = Ecoulement%DEtaDt(Int_sl0:Int_sl1)%perturbation
    
    ! Internal state vector: velocity and acceleration of the floater.
    jj = 1
    do nc = Int_Body,Mesh%NBody
        if(InputData%free_body(jj) .or. abs(InputData%ACorps(jj)).gt.Epsilon .or. abs(InputData%Vcst(jj)).gt.Epsilon)then
            RKVBody(1:6,jk,jj) = Mesh%Body(nc)%VBody(1:6)            
            RKABody(1:6,jk,jj) = Mesh%Body(nc)%ABody(1:6)
        end if
        jj = jj + 1
    end do
       
    ! Internal state vector: velocity of the nodes of the mesh.
    if (DeformMesh) then
        do j=1,Mesh%Nnoeud
            RKVel(1:3,j,jk) = Mesh%Tnoeud(j)%Velocity(1:3)
        enddo
    end if

end subroutine RK_manager

subroutine Time_stepping(Mesh,Mesh0,Ecoulement,Ecoulement0,fgeom_vect,RK,RKVBody,RKABody,RKVel,denom,Nnodes,Nnodes_FS,t,h,InputData,size_RKVBody)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(inout)                  :: Mesh                 ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                :: Mesh0
    type(TMaillage), intent(in)                     :: Mesh0                ! Local copy of Mesh in the RK4 loop
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement           ! Flow parameters (phi, etc.)
    !f2py integer*1, dimension(1000)                :: Ecoulement0
    type(TEcoulement),intent(in)                    :: Ecoulement0          ! Local copy of the flow parameter in the RK4 loop
    !f2py integer*1, dimension(1000)                :: fgeom_vect
    type(type_GeomVect), intent(inout)              :: fgeom_vect           ! Geometry of the floaters (not the domain)
    integer,intent(in)                              :: Nnodes               ! Number of nodes in the mesh
    integer,intent(in)                              :: Nnodes_FS            ! Number of nodes in the mesh of the free surface
    real(rp), dimension(Nnodes_FS,4,2),intent(in)   :: RK                   ! RK4 data for the flow
    real(rp), dimension(3,Nnodes,4),intent(in)      :: RKVel                ! RK4 data for the mesh
    integer,intent(in)                              :: size_RKVBody         ! NBodies
    real(rp), dimension(6,4,size_RKVBody),intent(in):: RKABody, RKVBody     ! RK4 data for the floater
    real(rp),intent(in)                             :: denom                ! 1/6
    real(rp),intent(in)                             :: t                    ! Current time
    real(rp),intent(in)                             :: h                    ! time step
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(inout)             :: InputData            ! InputData
    
    integer                                         :: Int_sl0,Int_sl1      ! Index of the free surface nodes
    integer                                         :: j,nc,jj              ! Parameter
    
    ! This subroutine time-steps the state vector.
    
    ! Index of the free surface nodes
    Int_sl0 = Mesh%FS%IndFS(1)
    Int_sl1 = Mesh%FS%IndFS(3)
    
    ! Time-stepping of DpPhi_pDt and DEta_pDt (particular derivative) on the free surface.
    Ecoulement%Phi(Int_sl0:Int_sl1)%perturbation = Ecoulement0%Phi(Int_sl0:Int_sl1)%perturbation + denom*dt*(RK(Int_sl0:Int_sl1,1,1)+2._RP*RK(Int_sl0:Int_sl1,2,1)+2._RP*RK(Int_sl0:Int_sl1,3,1)+RK(Int_sl0:Int_sl1,4,1))
    
    ! Time-stepping of Eta_p and DpEtaDt
    Ecoulement%Eta(Int_sl0:Int_sl1)%perturbation = Ecoulement0%Eta(Int_sl0:Int_sl1)%perturbation + denom*dt*(RK(Int_sl0:Int_sl1,1,2)+2._RP*RK(Int_sl0:Int_sl1,2,2)+2._RP*RK(Int_sl0:Int_sl1,3,2)+RK(Int_sl0:Int_sl1,4,2))
        
    ! Time-stepping of the position of the nodes
    if (DeformMesh) then
        do j=1,Mesh%Nnoeud
            Mesh%Tnoeud(j)%Pnoeud(1:3) = Mesh0%Tnoeud(j)%Pnoeud(1:3) + denom*dt*(RKVel(1:3,j,1) + 2._rp*RKVel(1:3,j,2) + 2._rp*RKVel(1:3,j,3) + RKVel(1:3,j,4))
        end do
    end if
        
    ! Time-step of the position, velocity and acceleration of the floater
    jj = 1
    do nc = int_Body,Mesh%NBody
        if (InputData%free_body(jj) .or. abs(InputData%ACorps(jj)).gt.Epsilon .or. abs(InputData%Vcst(jj)).gt.Epsilon)then
            Mesh%Body(nc)%MBody = denom*dt*(RKVBody(:,1,jj) + 2._RP*RKVBody(:,2,jj) + 2._RP*RKVBody(:,3,jj) + RKVBody(:,4,jj))
            Mesh%Body(nc)%CSolv = Mesh0%Body(nc)%CSolv + denom*dt*(RKVBody(:,1,jj) + 2._RP*RKVBody(:,2,jj) + 2._RP*RKVBody(:,3,jj) + RKVBody(:,4,jj))
            ! GBody.
            if(not(iFixPoint))then
                ! Motion equation solved at G (CSolv = (GBody,angles)).
                Mesh%Body(nc)%GBody = Mesh%Body(nc)%CSolv(1:3)
            end if
            ! If iFixPoint: GBody is updated in update_position when the subroutine Remesh is called.
            Mesh%Body(nc)%VBody = Mesh0%Body(nc)%VBody + denom*dt*(RKABody(:,1,jj) + 2._rp*RKABody(:,2,jj) + 2._rp*RKAbody(:,3,jj) + RKABody(:,4,jj))
            Mesh%Body(nc)%ABody = denom*(RKABody(:,1,jj) + 2._rp*RKABody(:,2,jj) + 2._rp*RKAbody(:,3,jj) + RKABody(:,4,jj))
        end if
        jj = jj + 1
    end do
    
    ! Time-stepping of the velocity of the nodes
    if (DeformMesh) then
        do j=1,Mesh%Nnoeud
            Mesh%Tnoeud(j)%Velocity(1:3) = denom*(RKVel(1:3,j,1) + 2._rp*RKVel(1:3,j,2) + 2._rp*RKVel(1:3,j,3) + RKVel(1:3,j,4))
        end do
    end if
    
    ! Mise a jour de la position des noeuds du maillage a partir de la vitesse des noeuds au temps ti (deformation)
    if (DeformMesh) call Remesh(Mesh, Mesh0, t, InputData,h, fgeom_vect) ! Could use RKVel 
    
end subroutine Time_stepping

subroutine Free_surface_RK_step(Mesh,Ecoulement,Ecoulement0,RK,h,jk,Nnodes_FS)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage), intent(in)                     :: Mesh             ! Total mesh (domain and floater)
    !f2py integer*1, dimension(1000)                :: Ecoulement
    type(TEcoulement),intent(inout)                 :: Ecoulement       ! Flow parameters (phi, etc.)
    !f2py integer*1, dimension(1000)                :: Ecoulement0
    type(TEcoulement),intent(in)                    :: Ecoulement0      ! Local copy of the flow parameter in the RK4 loop
    integer,intent(in)                              :: Nnodes_FS        ! Number of nodes in the mesh of the free surface
    real(rp), dimension(Nnodes_FS,4,2),intent(in)   :: RK               ! RK4 data for the flow
    real(rp),intent(in)                             :: h                ! Time step
    integer,intent(in)                              :: jk               ! RK4 loop parameter
    
    integer                                         :: Int_sl0,Int_sl1  ! Index of the free surface nodes
    
    ! This subroutine updates Phi_p and Eta_p for the next RK step.
    
    ! Index of the free surface nodes
    Int_sl0 = Mesh%FS%IndFS(1)
    Int_sl1 = Mesh%FS%IndFS(3)
    
    ! Phi_p
    Ecoulement%Phi(Int_sl0:Int_sl1)%perturbation = Ecoulement0%Phi(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,jk,1)*h
    
    ! Eta_p
    Ecoulement%Eta(Int_sl0:Int_sl1)%perturbation = Ecoulement0%Eta(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,jk,2)*h
    
end subroutine Free_surface_RK_step

subroutine Updating_nliss(CrossingFS,jFiltering)
	
    logical,intent(in)	    :: CrossingFS   ! A body crossed the free surface (True) or not (False).
    integer,intent(inout)   :: jFiltering   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
    
    ! This subroutine updates nliss in case of crossing the free surface.
    
    ! New nliss in case of crossing the FS.
    if(CrossingFS)then
        print*,"nliss was equal to ",nliss
        nliss = 1
        print*,"It is now equal to ",nliss
        jFiltering = 1
	end if
    
    if(jFiltering.ne.0)then
        jFiltering = jFiltering + 1
	end if
    
    ! Old nliss after 200 time steps.
    if(jFiltering.eq.200)then ! 200 is subjective.
        jFiltering = 0
        print*,"nliss was equal to ",nliss
        nliss = nliss_input
        print*,"It is now equal to ",nliss
	end if
    
end subroutine Updating_nliss

end module BoucleTemp
