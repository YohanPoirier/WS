!module Test_MeshDeformation
!    use Constantes
!    use Structuresdonnees
!    use FonctionsCommunes
!    use GeomDef
!    use Remesh_mod
!    use BodyMotion_mod
!    use MeshGen
!    use BoucleTemp
!    implicit none
!    contains
!    
!    subroutine TimeLoop_MeshDeformation(Mesh,fgeom_vect,nface,rep0,Grid,InputData)
!        
!    !f2py integer*1, dimension(1000)        :: Mesh
!    type(TMaillage), intent(inout)          :: Mesh                                         ! Total mesh (domain and floater)
!    !f2py integer*1, dimension(1000)        :: fgeom_vect
!    type(type_geom), intent(inout)          :: fgeom_vect                                   ! Geometry of the floater (not the domain)
!    integer,intent(in)                      :: nface                                        ! Number of faces in both the floater and the domain geometries
!    !f2py integer*1, dimension(1000)        :: rep0
!    type(repere3d),intent(in)               :: rep0                                         ! Inertial frame
!    !f2py integer*1, dimension(1000)        :: Grid
!    type(MGrid),intent(inout)               :: Grid                                         ! Transitional mesh
!    !f2py integer*1, dimension(1000)        :: InputData
!    type(InputDataStruct),intent(inout)     :: InputData                                    ! Input data
!    
!    logical                                 :: bgenerate                                    ! Boolean to know if the mesh needs to be generated or not
!    integer                                 :: j,jj,nc, jk, jt, ierror, Nnodes, Nnodes_FS   ! Parameters
!    real(rp)                                :: ti, h, time_begin, tmoy               
!    real(rp), dimension(nt)                 :: time_end                                     ! Final time of the RK4 loop
!    real(rp), allocatable                   :: RKVel(:,:,:)                                 ! Influence coefficients, RK4 data for the flow and the mesh
!    !f2py integer*1, dimension(1000)        :: Mesh0
!    type(TMaillage)                         :: Mesh0                                        ! Local copy of Mesh in the RK4 loop
!    real(rp), parameter                     :: denom = 1._RP/6._RP                          ! Denominateur of the RK4 mean
!    real(rp)                                :: periode                                      ! Wave period
!    real(rp), dimension(:),allocatable      :: t                                            ! Time vector
!    integer                                 :: NumBodyNodes                                 ! Index of the floater nodes of the current floater mesh
!    integer                                 :: size_Phi                                     ! Size of Phi
!    integer,dimension(:,:),allocatable      :: IndBody_former                               ! (1st node, last node)xNbodies: Index of the floater nodes of the former mesh when new mesh generation
!    
!    ! This subroutine is the temporal loop. It updates the state vector at each time step.
!    ! State vector = (Eta_p, Phi_p, Node_positions, floater_position, floater_velocity)
!    
!    ! Nnodes and Nnodes_FS
!    Nnodes = Mesh%Nnoeud
!    Nnodes_FS = Mesh%FS%IndFS(3) - Mesh%FS%IndFS(1) + 1
!    
!    ! Allocation
!    allocate(IndBody_former(2,NBodies))
!    IndBody_former = 0._RP
!    NumBodyNodes = 0
!    
!    jj = 1
!    do nc = Int_Body,Mesh%NBody
!        IndBody_former(1,jj) = Mesh%Body(nc)%IndBody(1)
!        IndBody_former(2,jj) = Mesh%Body(nc)%IndBody(3)
!        NumBodyNodes = NumBodyNodes + IndBody_former(2,jj) - IndBody_former(1,jj) + 1
!        jj = jj + 1
!    end do
!    
!    allocate(t(nt))
!    allocate(RKVel(3,Nnodes,4))
!    
!    ! Time vector
!    call Time_vector(t,nt,periode)
!    call cpu_time(time_begin)
!
!    ! Initialization of the output files
!    call PrePlotMeshDeformation()
!    
!    ! Initialisation of the copy of the mesh
!    call NewMaillage(Mesh0, Mesh%Nnoeud,NBodies+1) ! +1 for the tank.
!    
!    ! Initialisation des vitesses des points du maillage à 0
!    do j=1,Mesh%Nnoeud
!        Mesh%Tnoeud(j)%Velocity = 0._rp
!    end do
!    
!    ! Initialisation des position, vitesse et accélération du corps
!    jj = 1
!    do nc = Int_Body,Mesh%NBody
!        Mesh%Body(nc)%GBody(1:3) = InputData%PositionG(1:3,jj)
!        Mesh%Body(nc)%CSolv(1:3) = InputData%Position(1:3,1,jj)
!        Mesh%Body(nc)%CSolv(4:6) = InputData%Position(1:3,2,jj)
!
!        Mesh%Body(nc)%MBody = 0._rp
!        Mesh%Body(nc)%VBody = 0._rp
!        Mesh%Body(nc)%ABody = 0._rp
!        jj = jj + 1
!    end do
!    
!    call CopyMaillage(Mesh0, Mesh)
!    
!    ! Temporal loop
!    do jt=1,nt
!        
!        ! RK4 steps
!        do jk=1,4 
!            
!            ! Computation of the time step of the RK4 step and the current time
!            call Time_step_Current_time(t,jk,jt,ti,h)
!
!            ! New mesh if necessary
!            if (jk.eq.1 .and. jt.ne.1 ) then !.and. .false.
!                ! New mesh
!                call Regeneration_Mesh(Mesh,Ecoulement,ti,bgenerate,fgeom_vect,nface,rep0,Grid,IndBody_former,NBodies,Nnodes,Nnodes_FS,ierror,InputData)
!                  
!                ! Deallocation and new allocation
!                if (Mesh_type.eq.2) then
!                    if (bgenerate) then
!                        if(ierror == 0) then
!                            if(allocated(RKVel)) deallocate(RKVel)
!                            allocate(RKVel(3,Nnodes,4))
!                            call DelMaillage(Mesh0)
!                            call NewMaillage(Mesh0,Mesh%Nnoeud,NBodies) ! +1 for the tank.
!                        end if
!                    end if
!                end if
!            end if
!            
!            !----------------------------------------------------------------------------------
!            !           Computation of the time differentiation of the state vector
!            !----------------------------------------------------------------------------------
!            
!            ! Reinitialisation de Mesh0 et Ecoulement0 avant la resolution de la passe RK1 and writting the info in the output files.
!            if (jk.eq.1) then
!                call CopyMaillage(Mesh0, Mesh)
!            end if
!
!            ! Calcul des vitesses des noeuds du maillage (updating the velocity of the nodes)
!            if (DeformMesh) call MeshVel(Mesh, ti, 1._RP,InpuData) ! PYW: MeshVel uses a time step equal to 1.0 instead of h: validated with the floating cylinder test case.
!
!            ! Internal state vector: velocity of the nodes of the mesh
!            if (DeformMesh) then
!                do j=1,Mesh%Nnoeud
!                    RKVel(1:3,j,jk) = Mesh%Tnoeud(j)%Velocity(1:3)
!                enddo
!            end if
!
!            !----------------------------------------------------------------------------------
!            !               Preparation of the state vector for the next RK4 step
!            !----------------------------------------------------------------------------------
!            
!            if(jk.ne.4)then ! No updating for the 4th RK4 step.
!                
!                ! Computation of the time step of the RK4 step and the current time
!                call Time_step_Current_time(t,jk+1,jt,ti,h) ! jk+1 -> for the next RK4 step
!                
!                ! Position of the floater
!                call BodyMotion(Mesh,Mesh0,h) ! Could use RKVBody(1:6,jk)
!                
!                ! Mise a jour de la vitesse du corps au temps ti
!                call BodyVelocity(Mesh0, Mesh, ti, h,InputData) ! Could use RKABody(1:6,jk)
!                
!                ! Mise a jour de la position des noeuds du maillage a partir de la vitesse des noeuds au temps ti (deformation)
!                if (DeformMesh) call Remesh(Mesh, Mesh0, ti, h, fgeom_vect) ! Could use RKVel
!
!            end if
!            
!        end do ! End of RK4 steps
!    
!        ! Time-stepping of the velocity of the nodes
!        if (DeformMesh) then
!            do j=1,Mesh%Nnoeud
!                Mesh%Tnoeud(j)%Velocity(1:3) = denom*(RKVel(1:3,j,1) + 2._rp*RKVel(1:3,j,2) + 2._rp*RKVel(1:3,j,3) + RKVel(1:3,j,4))
!            end do
!        end if
!            
!        ! Mise a jour de la position des noeuds du maillage a partir de la vitesse des noeuds au temps ti (deformation)
!        if (DeformMesh) call Remesh(Mesh, Mesh0, ti, h, fgeom_vect) ! Could use RKVel
!        
!        call PlotMeshDeformation(ti, Mesh)
!    
!        ! Writting the time info in the command window
!        call Writting_time_info(t,time_end,time_begin,tmoy,jt)
!        
!    end do ! End of temporal loop
!    
!    ! Deallocating
!    if(allocated(RKVel))         deallocate(RKVel)
!    if(allocated(t))             deallocate(t)
!    
!    end subroutine TimeLoop_MeshDeformation
!    
!    subroutine TestMeshDeformation(Mesh,fgeom,nface,rep0,Grid)
!        !f2py integer*1, dimension(1000)                    :: Mesh
!        type(TMaillage), intent(inout)                      :: Mesh
!        !f2py integer*1, dimension(1000)                    :: fgeom
!        type(type_geom), intent(inout)                      :: fgeom                                ! Geometry of the floater (not the domain)
!        integer,intent(in)                                  :: nface                                ! Number of faces in both the floater and the domain geometries
!        !f2py integer*1, dimension(1000)                    :: rep0
!        type(repere3d), intent(in)                           :: rep0                                 ! Inertial frame
!        !f2py integer*1, dimension(1000)                    :: Grid
!        type(MGrid), intent(inout)                           :: Grid 
!        ! Subroutine to test the deformation of a floating body and the pierced FS meshes
!        logical, dimension(5) :: LogicalTemp
!        logical, dimension(6) :: ddl_dirTemp
!        real(rp) :: Tcorps
!        real(rp), dimension(2) :: Ttemp
!        type(TMaillage) :: MeshInit
!        ! 
!        call NewMaillage(MeshInit, Mesh%Nnoeud,NBodies+1) ! +1 for the tank.
!        MeshInit = Mesh
!        !call CopyMaillage(Mesh,MeshInit,Mesh%Nnoeud)
!        Ttemp = [ACorps,wcorps]
!        LogicalTemp = [DeformBody, DeformFS, Mesh%Body(Int_Body)%CMD(1), Free_Body, DeformMesh]
!        Free_Body = .false.
!        DeformBody = .true.
!        DeformFS   = .true.
!        DeformMesh = .true.
!        
!        Mesh%Body(Int_Body)%CMD(1) = .true.
!        ddl_dirTemp = dll_dir
!        dll_dir(3) = .true.
!    
!        ACorps = 0.125_RP*Mesh%Body(Int_Body)%DimBody(3)
!        wcorps = 1._RP
!        Tcorps = 2._RP*Pi/wcorps
!        print*, 'Test_MeshDeformation : ACorps = ', ACorps, ', Tcorps = ', Tcorps 
!        
!        call TimeLoop_MeshDeformation(Mesh,fgeom,nface,rep0,Grid)
!        
!        close(33)
!        call DelMaillage(Mesh)
!        call NewMaillage(Mesh,MeshInit%Nnoeud,NBodies+1) ! +1 for the tank.
!        Mesh = MeshInit
!        !call CopyMaillage(MeshInit,Mesh,MeshInit%Nnoeud)
!        call DelMaillage(MeshInit)
!        !
!        ACorps = Ttemp(1)
!        wcorps = Ttemp(2)
!        DeformBody = LogicalTemp(1)
!        DeformFS   = LogicalTemp(2)
!        Mesh%Body(Int_Body)%CMD(1) = LogicalTemp(3)
!        Free_Body = LogicalTemp(4)
!        DeformMesh = LogicalTemp(5)
!        dll_dir = ddl_dirTemp
!        !
!    end subroutine TestMeshDeformation
!    
!    subroutine PrePlotMeshDeformation()
!        integer :: ios
!        open(unit=33,file='Test_MeshDeformation'//filename, iostat=ios)
!        if (ios/=0) stop "Erreur à l'ouverture du fichier de sortie principal"
!        write(33,fmt='(50a)') 'Title = "Mesh de la cuve"'
!        write(33,fmt='(50a)') 'VARIABLES = "X","Y","Z"'
!    end subroutine PreplotMeshDeformation
!    
!    subroutine PlotMeshDeformation(t, Mesh)
!        ! Paramètres
!        real(rp), intent(in) :: t
!        !f2py integer*1, dimension(1000) :: Mesh
!        type(TMaillage), intent(in) :: Mesh
!        ! Variables Locales
!        integer :: j
!        character(len=10) :: num
!        !
!        write( num, '( f0.4 )' ) t
!        ! 
!        write(33,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Nnoeud, ', E=', Mesh%Nfacette,&
!        &  ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
!        do j=1,Mesh%Nnoeud
!            write(33,'(3E)') Mesh%Tnoeud(j)%Pnoeud(1:3)
!        end do
!        !
!        do j=1,Mesh%Nfacette
!        write(33,'(3I)') Mesh%Tfacette(j)%Tnoeud
!        end do
!    end subroutine PlotMeshDeformation
!    
!    subroutine CheckMeshAxisym(Mesh)
!    type(TMaillage) :: Mesh
!    integer :: j
!    real(rp) :: Radius
!    print*, '**********************************'
!    do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
!        Radius = norm2(Mesh%Body(Int_Body)%GBody(1:3)-Mesh%Tnoeud(j)%Pnoeud)
!        if (abs(Radius - Mesh%Body(Int_Body)%DimBody(1)).gt.Epsilon2) then
!            print*, j, Radius
!            print*, Mesh%Tnoeud(j)%Pnoeud
!        end if
!        
!    end do
!    end subroutine CheckMeshAxisym
!end module Test_MeshDeformation
