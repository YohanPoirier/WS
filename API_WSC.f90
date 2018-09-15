module API_WSC
    use Exec_mod
    use GeomDef
    use MeshGen
    use BoucleTemp
    use Parameters
    
    implicit none
    
    ! Parareal
    
    !f2py integer*1, dimension(1000)        :: Mesh_ref 
    type(TMaillage)                         :: Mesh_ref                               ! Maillage de reference (lineaire) pour comparer les différentes itérations
    !f2py integer*1, dimension(1000)        :: L_Ecoulements 
    type(T_liste_Ecoulements) :: L_Ecoulements  ! Liste d'ecoulements
    !f2py integer*1, dimension(1000)        :: L_Bodies
    type(T_liste_Bodies) :: L_Bodies  ! Liste des corps
    
    ! Variables of Main (Maillage -> Mesh ; Mesh -> Grid)
    !f2py integer*1, dimension(1000)        :: Mesh
    type(TMaillage)                         :: Mesh                                         ! Final mesh (named Maillage in Main.f90).
    integer                                 :: nface,nline                                  ! Number of faces and lines in both the floater and the domain geometries.
    integer                                 :: ierror                                       ! Error flag.
    !f2py integer*1, dimension(1000)        :: fgeom_vect
    type(type_GeomVect)                     :: fgeom_vect                                   ! Geometry of the floaters (not the domain).
    !f2py integer*1, dimension(1000)        :: fdomaine
    type(type_geom)                         :: fdomaine                                     ! Geometry of the domain.
    !f2py integer*1, dimension(1000)        :: tab2
    type(chaine_point_pt),dimension(100)    :: tab2                                         ! Table of intersection points ! Why 100? (PYW).
    integer                                 :: n_tab2,n_tab                                 ! Number of intersection curves and lines.
    !f2py integer*1, dimension(1000)        :: rep0
    type(repere3d)                          :: rep0                                         ! Inertial frame.
    !f2py integer*1, dimension(1000)        :: Grid
    type(MGrid)                             :: Grid                                         ! Transitional mesh (named Mesh in Main.f90).
    integer                                 :: nb_point,nb_tri                              ! Number of points and triangles in the Mgrid.
    integer                                 :: jt0                                          ! Initial time step parameter.
    real(rp)                                :: t_tmp                                        ! Back-up of t0 in case of input state file.
    integer                                 :: jt_init                                      ! First time step parameter.
    real(rp)                                :: Starting_time                                ! Starting time saved in case of calling a state input file.
    logical                                 :: get_State                                    ! True if the state input file is present, false otherwise.
    character(len=50)                       :: fileState                                    ! State input file.
    
    ! Variables of BoucleTemporelle_RK4
    logical                                 :: boolRemesh                                   ! Boolean to know if the mesh of the bodies or the free surface was regenerated (true) or not (false).
    logical                                 :: boolRemeshFS                                 ! = true: remeshing the free surface, = false: no remeshing.
    integer                                 :: j, nc,jj                                     ! Parameters.
    real(rp)                                :: ti, h, time_begin, rCI, tmoy, TimeRamp       ! Time parameters.
    real(rp),allocatable                    :: RKABody(:,:,:), RKVBody(:,:,:)               ! RK4 data for the floater integer.
    integer                                 :: Nnodes                                       ! Number of nodes in the mesh.
    integer                                 :: Nnodes_FS                                    ! Number of nodes in the mesh of the free surface.
    real(rp),dimension(:),allocatable       :: time_end                                     ! Final time of the RK4 loop.
    real(rp),allocatable                    :: CD(:,:), CS(:,:), RK(:,:,:), RKVel(:,:,:)    ! Influence coefficients, RK4 data for the flow and the mesh.
    real(rp),dimension(:),allocatable       :: Phi                                          ! Phi for the finite difference method.
    !f2py integer*1, dimension(1000)        :: Mesh0
    type(TMaillage)                         :: Mesh0                                        ! Local copy of Mesh in the RK4 loop.
    !f2py integer*1, dimension(1000)        :: Ecoulement, Ecoulement0, EcoulementDF
    type(TEcoulement)                       :: Ecoulement, Ecoulement0, EcoulementDF        ! Flow parameters (phi, etc.) and the copy of them.
    real(rp), parameter                     :: denom = 1._RP/6._RP                          ! Denominateur of the RK4 mean.
    real(rp)                                :: periode                                      ! Wave period.
    real(rp), dimension(:),allocatable      :: t                                            ! Time vector.
    integer,dimension(:,:),allocatable      :: IndBody_former                               ! (1st node, last node)xNbodies: Index of the floater nodes of the former mesh when new mesh generation.
    integer                                 :: size_Phi                                     ! Size of Phi.
    integer                                 :: NumBodyNodes                                 ! Index of the floater nodes of the current floater mesh.
    !f2py integer*1, dimension(1000)        :: Mesh_State
    type(TMaillage)                         :: Mesh_State                                   ! Mesh read in the state input file (WARNING: only Tnoeud and Tfacette are (partially) filled).
    !f2py integer*1, dimension(1000)        :: Ecoulement_State
    type(TEcoulement)                       :: Ecoulement_State                             ! Flow parameters read in the state input file.
    integer                                 :: jFiltering                                   ! Counter to wait for 200 time steps before using the original parameters in case of a body crossing the free surface.
    
    ! Variable of Regeneration_Mesh
    integer                                 :: nRemesh                                      ! Number of remeshing (use for the number of Advance_Front_Remesh.dat).
    logical                                 :: ForcedRemesh                                 ! Force the remeshing if necessary (crossing the free surface).
    logical                                 :: CrossingFS                                   ! A body crossed the free surface (True) or not (False).
    
    ! Variables of FreeBodyMotion
    integer                                 :: Ndof, N, Nb,size_NBody                       ! Number of dof, size of A, number of nodes in the mesh of the floater and size of the structure NBody.
    real(rp), dimension(:,:), allocatable   :: A                                            ! A matrix of the linear system.
    real(rp), dimension(:), allocatable     :: B, Sol                                       ! B and X matrices of the linear system.
    
    ! Variables of SystLin_FreeBody
    integer                                 :: Nsys                                         ! Number of unknowns in the mesh, nodes in the mesh of the floater and the number of degree of freedom.
    integer, dimension(3)                   :: Nsl                                          ! Index of the nodes for the free surface.
    integer, dimension(:,:), allocatable    :: NBody                                        ! Index of the nodes for the floater.
    real(rp), allocatable                   :: Fext(:,:)                                    ! External forces.
    real(rp), allocatable                   :: S(:,:,:),DSDt(:,:,:)                         ! Transformation matrices between the Cardan frames and the inertial frame and its time-differentiation.
    real(rp), allocatable                   :: M(:,:,:)                                     ! Mass matrices.
    real(rp), allocatable                   :: InertiaTerms(:,:)                            ! Inertia terms in the dynamical momentum expression.
    real(rp),dimension(3,3,2)               :: Tob                                          ! eR0.
    
    ! Variables of get_Inertia
    real(rp),dimension(6,6)                 :: Inertia_mat                                  ! Inertia matrix of the floater.
    
    ! Variable of get_MBody
    real(rp),dimension(6)                   :: MBody                                        ! Displacement of the floater.
    
    ! Variable of get_CSolv
    real(rp),dimension(6)                   :: CSolv                                        ! Position of the floater.
    
    ! Variable of get_VBody
    real(rp),dimension(6)                   :: VBody                                        ! Velocity of the floater.
    
    ! Variable of get_ABody
    real(rp),dimension(6)                   :: ABody                                        ! Acceleration of the floater.
    
    ! Variable of get_Gj
    real(rp),dimension(3)                   :: Gj                                           ! Center of gravity of the floater.
    
    ! Variable of Execution
    !f2py integer*1, dimension(1000)        :: InputData
    type(InputDataStruct)                   :: InputData                                    ! Input data.
    
    ! Variable of get_CMD
    integer                                 :: CMD                                          ! Mesh%CMD(1).
    
    ! Variable of get_CT_Body
    real(rp),dimension(:,:),allocatable     :: CT                                           ! CT coefficients.
    
    ! Variable of get_CK_Body
    real(rp),dimension(:,:),allocatable     :: CK                                           ! CK coefficients.
    
    ! Variable of get_Q_Body
    real(rp),dimension(:),allocatable       :: Q                                            ! Q coefficients.
    
    ! Variable of get_Th_Body
    real(rp),dimension(:),allocatable       :: Th                                           ! Th coefficients.
    
    ! Variable of get_Fhydro_Body
    real(rp),dimension(:),allocatable       :: Fhydro                                       ! Total hydrodynamic loads (DPhiDt%Perturbation included).
    
    ! Variable of get_DPhiDt_Body
    real(rp),dimension(:),allocatable       :: DPhiDt                                       ! DPhiDt.
    
    ! Variable of get_Active_Body
    logical                                 :: Active                                       ! Mesh%Body%Active
    
    ! Variable of get_InputData
    integer,dimension(:),allocatable        :: ndll                                         ! Number of DOF for each body.
    integer,dimension(:,:),allocatable      :: dll_dir                                      ! DOF for each body.
    
    ! Variable of compute_eJ0
    real(rp),dimension(6,6)                 :: eJ0                                          ! eJ0.
    
    ! Variable of compute_eJ0p
    real(rp),dimension(6,6)                 :: eJ0p                                         ! eJ0p.
    
    ! Variable of CDPhi
    real(rp),dimension(:),allocatable       :: CDPhi                                        ! CD(:FS)¨*Phi_t(FS).
    
    ! Variable of Initialization_solution_MB
    real(rp),dimension(:),allocatable       :: Sol_WSC                                      ! Solution for the GMRES solver for the WSC components.
    
    ! Variables of API_Geometry_Domain
    real(rp),dimension(3)                   :: O                                            ! Origine of the inertial frame.
    !f2py integer*1, dimension(1000)        :: VX,VY,VZ
    type(vector)                            :: VX,VY,VZ                                     ! Unit vector of the basis of the inertial frame.
    integer                                 :: nface_old                                    ! Number of faces for the previous body (in order to compute nface_vect.
    logical                                 :: AboveFS                                      ! All parts of the geometry are above the free surface (True) or not (False).
    
    ! Variables of API_Mesh0d
    integer                                 :: ios                                          ! Output parameter.
    integer                                 :: iflag                                        ! Flag.
    !f2py integer*1, dimension(1000)        :: P
    type(point)                             :: P                                            ! Point.
    integer                                 :: nb_arete                                     ! Number of edges in Mgrid.
    
    real(rp),dimension(20) :: time2 ! To be deleted
    

    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                       Subroutines of Main.f90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    ! Yohan ----------------------------------------
    
    
    subroutine print_file(text, io)
        character (len=50),intent(in) :: text
        integer :: io
        
        write(io,*) text
        
    end subroutine print_file
        
        
        
        
    subroutine import_temp(temp_py)
        real(rp), dimension(4,5000) :: temp_py
        !f2py intent(out) :: temp_py
        
        temp_py = temp_fort
        
    end subroutine import_temp
    
    
    ! Parareal
    
    subroutine API_open_file_debug()
        open(unit=1111,file="debug_file.txt",iostat=ios)
        
    end subroutine API_open_file_debug
    
    subroutine API_close_file_debug()
        close(unit=1111)
        
    end subroutine API_close_file_debug
        
    
    subroutine API_parareal_init(fileparam, filegeom, N_iterations, N_ordis)
        
        character (len=50),intent(in)           :: fileparam,filegeom   ! Input files
        integer :: N_iterations, N_ordis
        integer :: j,k
        
        ! Lecture des fichiers d'entree
        call API_Execution(fileparam,filegeom)
        close(ioIntersection)
        
        lineaireFS = .True.
        lineaireBody = .True.
        
        ! Construction d'un maillage qui servira de reference dans toute la suite
        
     
        call Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)
        
        call Generation_Mesh(Mesh_ref,fdomaine,fgeom_vect,nface,Grid,nb_point,nb_tri,ierror,InputData,get_State,tab2,n_tab2,n_tab)
    
        
        ! Initialisation de la liste des ecoulements
        allocate(L_ecoulements%G(N_iterations+1, N_ordis))
        allocate(L_ecoulements%lambda(N_iterations+1, N_ordis))
        allocate(L_ecoulements%F(N_ordis))
        

        
        do j = 1, N_iterations+1
            do k = 1, N_ordis
                call NewEcoulement(L_ecoulements%G(j,k),Mesh_ref%Nnoeud)
                call NewEcoulement(L_ecoulements%lambda(j,k),Mesh_ref%Nnoeud)
            end do
        end do
        
        do k = 1, N_ordis
            call NewEcoulement(L_ecoulements%F(k),Mesh_ref%Nnoeud)
        end do

        ! Initialisation de la liste des corps
        allocate(L_bodies%G(N_iterations+1, N_ordis, Mesh_ref%Nbody))
        allocate(L_bodies%lambda(N_iterations+1, N_ordis, Mesh_ref%Nbody))
        allocate(L_bodies%F(N_ordis, Mesh_ref%Nbody))
 
    end subroutine API_parareal_init
    
    
    
    
    subroutine API_parareal_save_G(i_iter, i_ordi)

        integer :: i_iter, i_ordi
        
        ! Interpolation vers le maillage de reference
        call Interpolation_FS(Mesh, Mesh_ref,Ecoulement,ti,.true.,ierror)! Interpolation vers le maillage de reference
        write(1111,*) "Save G"
        write(1111,*) Mesh%Nnoeud
        write(1111,*) Mesh_ref%Nnoeud
        
       ! Copie de l'ecoulement dans la liste        
        call CopyEcoulement( L_ecoulements%G(i_iter+1,i_ordi+1), Ecoulement, Mesh_ref%Nnoeud)

        ! Copie des corps dans la liste
        do j = 1, Mesh_ref%Nbody
            L_bodies%G(i_iter+1,i_ordi+1,j) = Mesh%Body(j)
        end do
        
        
    end subroutine API_parareal_save_G
    
    
    
    subroutine API_parareal_save_F(i_ordi)
        
        integer :: i_ordi
        integer :: j
        
        ! Interpolation vers le maillage de reference
        call Interpolation_FS(Mesh, Mesh_ref,Ecoulement,ti,.true.,ierror)
        write(1111,*) "Save F"
        write(1111,*) Mesh%Nnoeud
        write(1111,*) Mesh_ref%Nnoeud
        
       ! Copie de l'ecoulement dans la liste        
        call CopyEcoulement( L_ecoulements%F(i_ordi+1), Ecoulement, Mesh_ref%Nnoeud)
        
        ! Copie des corps dans la liste
        do j = 1, Mesh_ref%Nbody
            L_bodies%F(i_ordi+1,j) = Mesh%Body(j)
        end do
        
        
    end subroutine API_parareal_save_F
    
    
    
    
    subroutine API_parareal_calcul_lambda(i_iter, i_ordi, jt, filestate_out)
    
        integer :: i_iter, i_ordi
        character(len=50),intent(in) :: filestate_out
        integer, intent(in) ::jt
        
        ! Ecoulement
        call DelEcoulement(Ecoulement)
        call NewEcoulement(Ecoulement, Mesh_ref%Nnoeud)
        call IniEcoulement(Ecoulement, Mesh_ref%Nnoeud, 0._RP)
        
        ! Phi_p, Eta_p.
        do j = Mesh_ref%FS%IndFS(1),Mesh_ref%FS%IndFS(3)
            Ecoulement%Phi(j)%perturbation = L_ecoulements%F(i_ordi+1)%Phi(j)%perturbation &
            + L_ecoulements%G(i_iter+2, i_ordi+1)%Phi(j)%perturbation &
            - L_ecoulements%G(i_iter+1, i_ordi+1)%Phi(j)%perturbation
            
            Ecoulement%Eta(j)%perturbation = L_ecoulements%F(i_ordi+1)%Eta(j)%perturbation &
            + L_ecoulements%G(i_iter+2, i_ordi+1)%Eta(j)%perturbation &
            - L_ecoulements%G(i_iter+1, i_ordi+1)%Eta(j)%perturbation
        end do
      
        ! Corps
        do nc = 1,Mesh_ref%NBody
            Mesh_ref%Body(nc)%MBody = L_bodies%F(i_ordi+1, nc)%MBody + L_bodies%G(i_iter+2, i_ordi+1, nc)%MBody -L_bodies%G(i_iter+1, i_ordi+1, nc)%MBody
            Mesh_ref%Body(nc)%CSolv = L_bodies%F(i_ordi+1, nc)%CSolv + L_bodies%G(i_iter+2, i_ordi+1, nc)%CSolv -L_bodies%G(i_iter+1, i_ordi+1, nc)%CSolv
            Mesh_ref%Body(nc)%GBody = L_bodies%F(i_ordi+1, nc)%GBody + L_bodies%G(i_iter+2, i_ordi+1, nc)%GBody -L_bodies%G(i_iter+1, i_ordi+1, nc)%GBody
            Mesh_ref%Body(nc)%VBody = L_bodies%F(i_ordi+1, nc)%VBody + L_bodies%G(i_iter+2, i_ordi+1, nc)%VBody -L_bodies%G(i_iter+1, i_ordi+1, nc)%VBody
        end do
        
        call Write_State2(Mesh_ref,Ecoulement,ti,jt,Starting_time,jFiltering, filestate_out)

    end subroutine API_parareal_calcul_lambda
    
    
    
    
    
    subroutine API_Write_State_with_interpolation(filename, jt)
    
        character(len=50),intent(in) :: filename
        integer, intent(in) ::jt

        
        call Interpolation_FS(Mesh, Mesh_ref,Ecoulement,ti,.true.,ierror)
        call Write_State2(Mesh_ref,Ecoulement,ti,jt,Starting_time,jFiltering, filename)
        
    end subroutine
    
    
    
    
    subroutine API_open_file_WP(filename, io)
    
        character (len=50),intent(in)       :: filename         ! Fichier d'ecriture
        integer, intent(in)                 :: io               ! Indice du fichier ou on ecrit

        open(unit=io,file=filename,iostat=ios)
        if(ios/=0) stop "Erreur lors de l'ouverture du fichier Wave_elevation"
        write(io,'(50a)') 'Title = "Wave elevation"'
        write(io,'(50a)') 'VARIABLES = "t","Eta_Incident","Eta_Perturbation","Eta"'
    
    end subroutine
    
    
    
    
    
    subroutine API_write_WP(io, nc)
    
        integer, intent(in)                 :: nc               ! Indice de la probe
        integer, intent(in)                 :: io               ! Indice du fichier ou on ecrit
    
        call PlotWaveElevation2(ti,InputData,Mesh,Ecoulement, io, nc)
    end subroutine
    
    
    
    
    subroutine API_close_file_WP(io)
    
        integer, intent(in)                 :: io               ! Indice du fichier ou on ecrit
    
        close(unit = io)

    end subroutine
    
    
    
    
    
    subroutine import_syst(Nn, A_py, B_py)
        
        integer :: Nn
        real(rp), dimension(Nn,Nn) :: A_py
        real(rp), dimension(Nn) :: B_py
          
        !f2py intent(in) :: Nn
        !f2py intent(out) :: A_py, B_py
        
        A_py = A
        B_py = B
                       
    end subroutine
    
    subroutine import_indices(ind)

        integer, dimension(6) :: ind

        !f2py intent(out) :: ind
    
                
        if(cuve_ferme)then
            if(Mesh%Body(Mesh%NBody)%Active)then
                ind = [1,Mesh%Nsys, Mesh%FS%IndFS(1),Mesh%FS%IndFS(3),Mesh%Body(1)%IndBody(1),Mesh%Body(Mesh%NBody)%IndBody(3)] ! Body(Mesh%NBody) is not above the free surface.
            else
                ind = [1,Mesh%Nsys, Mesh%FS%IndFS(1),Mesh%FS%IndFS(3),Mesh%Body(1)%IndBody(1),Mesh%Body(Mesh%NBody-1)%IndBody(3)] ! Body(Mesh%NBody) is above the free surface.
            end if
        else
            ind = [1,Mesh%Nsys, Mesh%FS%IndFS(1),Mesh%FS%IndFS(3),Mesh%Body(1)%IndBody(1),Mesh%Body(Mesh%NBody)%IndBody(3)]
        end if
    
    
    end subroutine 
    
    
    subroutine import_ecoulement(N, ind, X_connu, X_inconnu_0)
        
        integer :: N
        integer, dimension(6) :: ind
        real(rp), dimension(N) :: X_connu, X_inconnu_0
        
        !f2py intent(in) :: N, ind
        !f2py intent(out) :: X_connu, X_inconnu_0
        
        X_connu(ind(3):ind(4)) = Ecoulement%Phi(ind(3):ind(4))%perturbation
        X_connu(ind(5):ind(6)) = Ecoulement%DPhiDn(ind(5):ind(6))%perturbation
    
        X_inconnu_0(ind(3):ind(4)) = Ecoulement%DPhiDn(ind(3):ind(4))%perturbation
        X_inconnu_0(ind(5):ind(6)) = Ecoulement%Phi(ind(5):ind(6))%perturbation
        
    end subroutine
    
    
    subroutine export_ecoulement(N, ind, sol)
        integer :: N
        integer, dimension(6) :: ind
        real(rp), dimension(N) :: sol
        
        !f2py intent(in) :: N, ind, sol
        
        Ecoulement%DPhiDn(ind(3):ind(4))%perturbation = sol(ind(3):ind(4))
        Ecoulement%Phi(ind(5):ind(6))%perturbation = sol(ind(5):ind(6))
        
    end subroutine
    
    
    subroutine import_Ldom(Ldom_py)
        real(rp), dimension(5) :: Ldom_py
        !f2py intent(out) :: Ldom_py
        
        Ldom_py = Ldom
        
    end subroutine
    
    subroutine import_Mesh_dim(Nf,Nn,N_body)
        integer :: Nf, Nn, N_body
        !f2py intent(out) :: Nf, Nn, N_body
        
        Nf = Mesh%Nfacette
        Nn = Mesh%Nnoeud
        N_body = Mesh%Nbody
         
    end subroutine
    
    
    subroutine import_Mesh(Nf, Nn, L_P, L_ds, L_T, L_G, L_N, L_Rmax, L_type)
        
        integer :: Nf, Nn
        real(rp), dimension(Nf) :: L_Rmax
        integer, dimension(Nf,3) :: L_T
        real(rp), dimension(Nf,3) :: L_G, L_N
        real(rp), dimension(Nf,3,3) :: L_ds
        real(rp), dimension(Nn,3) :: L_P
        integer, dimension(Nn) :: L_type


        integer :: i
        
        !f2py intent(in) :: Nf, Nn
        !f2py intent(out) :: L_P, L_ds, L_T, L_G, L_N, L_Rmax, L_type
        
        do i = 1, Nn
            L_P(i,:) = Mesh%Tnoeud(i)%Pnoeud
            L_type(i) = Mesh%Tnoeud(i)%TypeNoeud
        end do
        
        do i = 1, Nf
            L_ds(i,:,:) = Mesh%Tfacette(i)%ds
            L_T(i,:) = Mesh%Tfacette(i)%Tnoeud - 1
            L_G(i,:) = Mesh%Tfacette(i)%Gfacette
            L_N(i,:) = Mesh%Tfacette(i)%Normale
            L_Rmax(i) = Mesh%Tfacette(i)%Rmax
        end do
               
    end subroutine
    
    
    
       
    subroutine import_Mesh_vois(Nn, Nmax, L_vois, L_vois_N) ! A supprimer
        
        integer :: Nn, Nmax
        integer, dimension(Nn,Nmax) :: L_vois
        integer, dimension(Nn) :: L_vois_N
               
        integer :: i,j
        
        !f2py intent(in) :: Nn, Nmax
        !f2py intent(out) :: L_vois, L_vois_N
        
        do i = 1, Nn
            L_vois_N(i) = Mesh%Tnoeud(i)%Nfacette 
            do j = 1, Mesh%Tnoeud(i)%Nfacette
                L_vois(i,j) = Mesh%Tnoeud(i)%Tfacette(j,1) - 1
            end do
        end do
                       
    end subroutine
    
    subroutine export_CI(Nn, CD_py, CS_py)
        
        integer :: Nn
        real(rp), dimension(Nn,Nn) :: CD_py, CS_py
           
        integer :: i
        
        !f2py intent(in) :: Nn, CD_py, CS_py

        CD = CD_py
        CS = CS_py
        
    end subroutine
    
    

    
    subroutine import_CI(Nn, CD_py, CS_py)
        
        integer :: Nn
        real(rp), dimension(Nn,Nn) :: CD_py, CS_py
          
        !f2py intent(in) :: Nn
        !f2py intent(out) :: CD_py, CS_py
        
        CD_py = CD
        CS_py = CS
                       
    end subroutine
    ! --------------------------------------------------------
    
    
    subroutine API_Execution(fileparam,filegeom,get_State_in,fileState_in)
        
        character (len=50),intent(in)           :: fileparam,filegeom   ! Input files
        logical,intent(in),optional             :: get_State_in         ! True if the state input file is present, false otherwise.
        character(len=50),intent(in),optional   :: fileState_in         ! State input file.
                        
        ! This subroutine wraps the subroutine Execution.

        
        ! get_State.
        if(present(get_State_in))then
            get_State = get_State_in
            fileState = fileState_in
        else ! Not present.
            get_State = .false.
        end if
        
        if (allocated(InputData%free_body)) then
            call Delete_InputData(InputData)     
        end if
        
        ! Execution.
        call Execution(fileparam,filegeom,InputData,get_State)
        

        
        ! Updating InputData in case of state input file.

        if(get_State)then
            call read_State_InputData(fileState,InputData,t0,jt0)
            print*,"Simulation starts at t = ",t0," s"
            ti = t0
            t_tmp = t0 ! Back-up t0 of *.in.
        end if


        
    end subroutine API_Execution
    
    subroutine API_Geometry()
    
        ! This subroutine wraps the subroutine Generation_Geometry.
        call Generation_Geometry(fgeom_vect,fdomaine,nface,tab2,n_tab2,rep0,InputData,ierror,n_tab)
    
    end subroutine API_Geometry
    
    subroutine API_Mesh()
                
        ! This subroutine wraps the subroutine Generation_Mesh.
        
        call Generation_Mesh(Mesh,fdomaine,fgeom_vect,nface,Grid,nb_point,nb_tri,ierror,InputData,get_State,tab2,n_tab2,n_tab)
        !if(get_State)then
        !    t0 = t_tmp ! t0 is used in the generation of the mesh (explicit wave elevation).
        !end if
    
    end subroutine API_Mesh
    
    subroutine API_BoucleTemporelle_RK4()
                
        ! This subroutine wraps the subroutine BoucleTemporelle_RK4.
        
        call BoucleTemporelle_RK4(Mesh, fgeom_vect, fdomaine, nface, Grid,nb_point,nb_tri,InputData,n_tab2,n_tab,get_State,fileState,jt0)
        
    end subroutine API_BoucleTemporelle_RK4
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Subroutines to create a developed temporal loop in Python
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine Pre_Temporal_loop()
            
        ! This subroutine calls all the subroutines which are used before the beginning of the temporal loop.
        
        ! Error.
        ierror = 0
        
        ! ForcedRemesh and CrossingFS.
        ForcedRemesh = .false.
        CrossingFS = .false.
        
        ! nRemesh.
        nRemesh = 0
        
        ! jFiltering.
        jFiltering = 0
        
        ! Nnodes and Nnodes_FS.
        Nnodes = Mesh%Nnoeud
        Nnodes_FS = Mesh%FS%IndFS(3) - Mesh%FS%IndFS(1) + 1
        
        ! Allocation.
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
        allocate(time_end(nt))
        allocate(CS(Nnodes,Nnodes), CD(Nnodes,Nnodes))
        allocate(RK(Nnodes_FS,4,2))
        allocate(RKVBody(6,4,Nbodies),RKABody(6,4,Nbodies))
        if (DeformMesh) allocate(RKVel(3,Nnodes,4)) ! RKVel != RK
        if (FiniteDifference_ForceCalcul) then
            size_Phi = NumBodyNodes
            allocate(Phi(NumBodyNodes))
            Phi = 0._RP
        end if
        
        ! Time vector
        call Time_vector(t,nt,periode,InputData)
        
        ! Initialization of the output files
        call PrePlots(InputData,Mesh%NBody,get_State)
        
        ! Initialization of Mesh0, Ecoulement, Ecoulement0 and eventually EcoulementDF for the updating at the end of each RK4 step
        call Initialization_Mesh_Ecoulement(Mesh,Mesh0,Ecoulement,Ecoulement0,t(1),EcoulementDF,InputData)
        
        ! Updating Ecoulement in case of state input file.
        

        if(get_State)then
            if(get_State)then
				
                ! Reading the state input file.
                call read_State(fileState,Mesh_State,Ecoulement_State,Starting_time,jFiltering)
                
                ! Updating nliss in case of crossing the free surface.
                if(jFiltering.ne.0)then
                    nliss = 1
		        end if
                
                ! Interpolation of Phi_p and Eta_p from the Mesh_State.
                call Interpolation_FS(Mesh_State,Mesh,Ecoulement_State,t(1),.true.,ierror) ! True because there is a FS remeshing, therefore an interpolation.
                write(1111,*) "Pretemporal"
                write(1111,*) Mesh_State%Nnoeud
                write(1111,*) Mesh%Nnoeud
                
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
                            
        end if
        
        ! First calculation of the influence coefficient to know how long the simulation will lasts.
        call Time_evaluation(Mesh,CD,CS,time_begin,time_end,Nnodes,nt,get_state,jt_init,Starting_time)
        
    end subroutine Pre_Temporal_loop
    
    subroutine API_Time_step_Current_time(jt,jk)
    
        integer,intent(in) :: jt,jk ! Temporal loop and RK4 parameters
        
        ! This subroutine wraps the subroutine Time_step_Current_time.
        
        call Time_step_Current_time(t,jk,jt,ti,h)
    
    end subroutine API_Time_step_Current_time
    
    subroutine API_rCI_manager(jk)
    
        integer,intent(in) :: jk ! RK4 parameter
        
        ! This subroutine wraps the subroutine rCI_manager.
        
        call rCI_manager(jk,rCI)
    
    end subroutine API_rCI_manager
    
    subroutine Initialization_Mesh0_Ecoulement0()
    
        ! This subroutine re-initializes Mesh0 and Ecoulement0.
        
        call CopyMaillage(Mesh0, Mesh)
        call CopyEcoulement(Ecoulement0, Ecoulement, Mesh%Nnoeud)
    
    end subroutine Initialization_Mesh0_Ecoulement0
    
    subroutine API_Plots()
        
        ! This subroutine wraps the subroutine Plots.
        
        call Plots(ti, Mesh, Ecoulement,InputData)
        
    end subroutine API_Plots
    
    subroutine API_Incident()
    
        ! This subroutine wraps the subroutine Incident.
        
        call Incident(ti, Mesh, Ecoulement)
    
    end subroutine API_Incident
    
    subroutine API_BodyCondition()
    
        ! This subroutine wraps the subroutine BodyCondition.
    
        call BodyCondition(ti, Mesh, Ecoulement)
        
    end subroutine API_BodyCondition
    
    subroutine API_MeshVel()
    
        ! This subroutine wraps the subroutine MeshVel.
        
        call MeshVel(Mesh, ti, 1._RP,InputData)
        
    end subroutine API_MeshVel
    
    subroutine API_solBVP(bool_CUDA)
    
        logical :: bool_CUDA
        ! This subroutine wraps the subroutine solBVP.
        
        time2 = 0._RP ! To be deleted.
        if (bool_CUDA) then
            call solBVP2( Ecoulement, Mesh, CD, CS, Nnodes,time2,boolRemesh, rCI)
        else
            !allocate(A(Mesh%Nsys,Mesh%Nsys), B(Mesh%Nsys))
            !call solBVP3( Ecoulement, Mesh, CD, CS,A,B, Nnodes,time2,boolRemesh, rCI)
            call solBVP( Ecoulement, Mesh, CD, CS, Nnodes,time2,boolRemesh, rCI)
        end if
        
        
    end subroutine API_solBVP
    
    subroutine API_Derive()
    
        ! This subroutine wraps the subroutine Derive.
        
        call Derive(ti, Ecoulement, Mesh)
        
    end subroutine API_Derive
    
    subroutine API_FreeBodyMotion()
    
        ! This subroutine wraps the subroutine FreeBodyMotion.
        
        time2 = 0._RP ! To be deleted.
        call FreeBodyMotion(Mesh, Ecoulement, CD, CS, ti,Nnodes,InputData,time2)
        
    end subroutine API_FreeBodyMotion
    
    subroutine API_Finite_differences(jt)
                    
        integer,intent(in) :: jt ! Time loop parameter.
    
        ! This subroutine wraps the subroutine Finite_differences.
        
        call Finite_differences(t(jt-1),Mesh,Ecoulement,EcoulementDF,Phi,size_Phi,h,InputData)
    
    end subroutine API_Finite_differences
    
    subroutine ForceBodyMotion_solBVP()

        ! This subroutine computes the hydrodynamic loads if the floater is still and without finite difference approximation.
        
        time2 = 0._RP !  To be deleted
        call ForceBodyMotion(Mesh, Ecoulement, ti,InputData)

        call solBVP(Ecoulement, Mesh, CD, CS,Nnodes, time2, .false., ti, .true.) ! .false. should be deleted.

    end subroutine ForceBodyMotion_solBVP
    
    subroutine API_PlotForces()
    
        ! This subroutine wraps the subroutine PlotForces.
        
        call PlotForces(Mesh,Ecoulement,ti,InputData)
        
    end subroutine API_PlotForces
     
    subroutine API_RK_manager(jk)
        
        integer,intent(in) :: jk ! RK4 parameter
    
        ! This subroutine wraps the subroutine RK_manager.
        
        call RK_manager(Mesh,Ecoulement,RK,RKVBody,RKABody,RKVel,jk,Nnodes,Nnodes_FS,InputData,NBodies)
        
    end subroutine API_RK_manager
    
    subroutine API_BodyMotion()
    
        ! This subroutine wraps the subroutine BodyMotion.
        
        call BodyMotion(Mesh,Mesh0,h)
        
    end subroutine API_BodyMotion
    
    subroutine API_BodyVelocity()
        
        ! This subroutine wraps the subroutine BodyVelocity.
        
        call BodyVelocity(Mesh0, Mesh, ti, h,InputData)
        
    end subroutine API_BodyVelocity
    
    subroutine API_Free_surface_RK_step(jk)
    
        integer,intent(in) :: jk ! RK4 parameter
    
        ! This subroutine wraps the subroutine Free_surface_RK_step.
        
        call Free_surface_RK_step(Mesh,Ecoulement,Ecoulement0,RK,h,jk,Nnodes_FS)
        
    end subroutine API_Free_surface_RK_step
    
    subroutine API_Remesh()
    
        ! This subroutine wraps the subroutine Remesh.
        
        call Remesh(Mesh, Mesh0, ti, InputData,h, fgeom_vect)
        
    end subroutine API_Remesh
    
    subroutine API_Regeneration_Mesh()
    
        ! This subroutine wraps the subroutine Regeneration_Mesh.
        
        call Regeneration_Mesh(Mesh,Ecoulement,ti,boolRemesh,boolRemeshFS,fgeom_vect,fdomaine,nface,Grid,nb_point,nb_tri,IndBody_former,NBodies,Nnodes,Nnodes_FS,ierror,InputData,nRemesh,n_tab2,n_tab,ForcedRemesh,CrossingFS)
        
    end subroutine API_Regeneration_Mesh
    
    subroutine New_allocation(jt)
        
        integer,intent(in) :: jt ! Time loop parameter.
        
        ! This subroutine deallocates the former variables and allocates the new variables if case of a new mesh.
        
        ! Updating nliss in case of crossing the free surface.
        call Updating_nliss(CrossingFS,jFiltering)
        
        if(Mesh_type.eq.2)then ! Mesh strategy of CC.
            
            if(boolRemesh)then
                
                if(ierror == 0)then ! Remeshing was sucessful.
                    
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
                end if
                                
            end if
            
        else ! Mesh strategy of LL.
                    
            ! State vector.
            if(iwState)   call Write_State(Mesh,Ecoulement,ti,jt,time_begin,jFiltering) ! At each time step.
                    
        end if
        
    end subroutine New_allocation
    
    subroutine API_Write_State(filename, jt)
    
        character(len=50),intent(in) :: filename
        integer, intent(in) ::jt
        
        print*,""
        print*,""
        print*, "Nombre de noeuds", Mesh%Nnoeud
        print*,""
        call Write_State2(Mesh,Ecoulement,ti,jt,Starting_time,jFiltering, filename)
        
    end subroutine
    
    
    subroutine API_Zeroing_new_points()
    
        ! This subroutine wraps the subroutine Zeroing_new_points.
        
        call Zeroing_new_points(Mesh,Mesh0,Ecoulement,Ecoulement0,boolRemesh,boolRemeshFS,IndBody_former,NBodies,ierror)
        
    end subroutine API_Zeroing_new_points
    
    subroutine API_Time_stepping()
    
        ! This subroutine wraps the subroutine Time_stepping.
        
        call Time_stepping(Mesh,Mesh0,Ecoulement,Ecoulement0,fgeom_vect,RK,RKVBody,RKABody,RKVel,denom,Nnodes,Nnodes_FS,ti,dt,InputData,NBodies)
        
    end subroutine API_Time_stepping
    
    subroutine API_Writting_time_info(jt)
    
        integer,intent(in)          :: jt           ! Temporal loop parameter
                
        ! This subroutine wraps the subroutine Writting_time_info.
        
        call Writting_time_info(t,time_end,time_begin,tmoy,jt,get_State,jt_init)
        
    end subroutine API_Writting_time_info
    
    subroutine API_Lissage(jt)
    
        integer,intent(in) :: jt ! Temporal loop parameter
        
        ! This subroutine wraps the subroutine LissageGaussien.
        
        if (nliss.ne.0) then
            if(mod(jt+1,nliss).eq.0)then
                call LissageGaussien(Mesh,Ecoulement,ierror)
                if(ierror/=0)then
                    ierror = 200
                    print*,"API_Lissage: Erreur dans le lissage: division par 0."
                endif
            endif
        end if
    
    end subroutine API_Lissage
    
    subroutine API_Closing()
    
        ! This subroutine wraps the subroutine Closing.
        
        call Closing(nt,time_end,time_begin,Mesh%NBody)
        
    end subroutine API_Closing
    
    subroutine Deallocating_Temporal_Loop()
    
        ! This subroutine deallocates all the variables of the temporal loop.

        if(allocated(Phi))           deallocate(Phi)
        if(allocated(CD))            deallocate(CD)
        if(allocated(CS))            deallocate(CS)
        if(allocated(RK))            deallocate(RK)
        if(allocated(RKVel))         deallocate(RKVel)
        if(allocated(t))             deallocate(t)
        if(allocated(RKABody))      deallocate(RKABody)
        if(allocated(RKVBody))      deallocate(RKVBody)
        if(allocated(IndBody_former))deallocate(IndBody_former)
        if(allocated(time_end))      deallocate(time_end)

    end subroutine Deallocating_Temporal_Loop
    
    ! Subroutines of FreeBodyMotion
    
    subroutine Allocation_FreeBodyMotion()
        
        ! This subroutine allocates A,B and Sol of the subroutine FreeBodyMotion.
        
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
        N = Mesh%Nsys + Nb + Ndof
        size_NBody = Mesh%NBody
        allocate(A(N,N), B(N), Sol(N))
        
    end subroutine Allocation_FreeBodyMotion
    
    subroutine API_Gradient()
    
        ! This subroutine wraps the subroutine Gradient.
        
        call Gradient(Mesh,Ecoulement)
    
    end subroutine API_Gradient
    
    subroutine API_AccConvec()
    
        ! This subroutine wraps the subroutine AccConvec.
        
        call AccConvec(Mesh, Ecoulement, ti) ! t -> ti
        
    end subroutine API_AccConvec
    
    subroutine API_Inertia()
    
        ! This subroutine wraps the subroutine Inertia.
        
        call Inertia(Mesh,InputData,.false.)
        
    end subroutine API_Inertia
    
    subroutine API_SystLin_FreeBody()
    
        ! This subroutine wraps the subroutine SystLin_FreeBody.
        
        call SystLin_FreeBody(CD, CS, Ecoulement, Mesh, A, B, Sol, ti, size_NBody,Nnodes, N,InputData)
        
    end subroutine API_SystLin_FreeBody
    
    subroutine Solving_SystLin()
    
        ! This subroutine solved the linear system.
        
        if (Solv.eq.0) then
            call GMRES(A, B, Sol, N,ierror)
        else  
            call LU(A, B, Sol, N, ierror)
        end if
        
        if(ierror/=0) goto 9999
        
        9999 continue
        if(ierror/=0)then
            write(*,90),ierror
        endif
        90      format('error #',i3,' : pb. in resolution of free motion.')
        
    end subroutine Solving_SystLin
    
    subroutine API_postSL_FreeBody()
    
        ! This subroutine wraps the subroutine postSL_FreeBody.
        
        call postSL_FreeBody(Sol, N,Mesh, Ecoulement,size_NBody,InputData)
    
    end subroutine API_postSL_FreeBody
    
    subroutine API_postSL_FreeBody_MB(Nsystem)
            
        integer,intent(in)  :: Nsystem  ! Nsystem
        
        ! This subroutine wraps the subroutine postSL_FreeBody.
        
        call postSL_FreeBody_MB(Sol, Nsystem,Mesh, Ecoulement,size_NBody)
    
    end subroutine API_postSL_FreeBody_MB
    
    subroutine API_PlotBVPPhit_Deallocation()
    
        ! This subroutine wraps the subroutine PlotBVPPhit and deallocates A, B and Sol.
    
        !if(iwphit)then
        !    call PlotBVPPhit(ti, Mesh, Ecoulement)
        !endif
    
        deallocate(A, B, Sol)
        
    end subroutine API_PlotBVPPhit_Deallocation
    
    ! Subroutines of SystLin_FreeBody
    
    subroutine Computation_M_CK_CT_Q_Th()
    
        ! This subroutine computes the matrices M, CK, CT, Q and Th.
        
        ! Time ramp
        call Ramp_Management(ti,TimeRamp)

        ! Number of nodes and index
        allocate(NBody(size_NBody,3))
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

        ! Computation of mass matrix and the inertia terms used in the dynamic momentum
        call get_Inertia_terms(Mesh,M,S,DSDt,InertiaTerms,InputData,NBodies,Tob)

        ! Computation of Q
        call get_Q(Mesh,Ecoulement,TimeRamp)

        ! External forces
        call ExternalForces(Mesh, Ecoulement, Fext, ti, 1,InputData,NBodies) ! In the coupling between InWave and the WSC code, the weight is computed in InWave that's why WeightOrNot = 0.

        ! Computation of Th
        call get_Th(Mesh,Ecoulement,Fext,InertiaTerms,TimeRamp,InputData,NBodies)

    end subroutine Computation_M_CK_CT_Q_Th

    subroutine Computation_M_CK_CT_Q_Th_MB()
        
        ! This subroutine computes the matrices M, CK, CT, Q and Th with the multibody theory.
        
        ! Time ramp
        call Ramp_Management(ti,TimeRamp)
        
        ! Number of nodes and index
        allocate(NBody(size_NBody,3))
        call Index_Number_Nodes(Mesh,Nsl,NBody,Nsys,Nb,size_NBody)
        
        ! Allocation of CK, CT, Q, Th, S, DSDt, Fext, M, InertiaTerms
        call Init_CK_CT_Q_Th(Mesh,InputData)
        allocate(S(3,3,NBodies))
        allocate(DSDt(3,3,NBodies))
        allocate(Fext(6,NBodies))
        allocate(Inertiaterms(3,NBodies))
        allocate(M(6,6,NBodies))
        
        ! Computation of S
        call get_S(Mesh,S,NBodies)
        
        ! Computation of DSDt
        call get_DSDt(Mesh,DSDt,NBodies)
        
        ! Computation of CK
        call get_CK(Mesh)
        
        ! Computation of CT (updating of CK).
        call get_CT(Mesh,S,NBodies)
          
        ! Inertia terms
        call get_Inertia_terms(Mesh,M,S,DSDt,InertiaTerms,InputData,NBodies,Tob)
        InertiaTerms = 0._RP
        
        ! Computation of Q
        call get_Q(Mesh,Ecoulement,TimeRamp)
        
        ! External forces
        call ExternalForces(Mesh, Ecoulement, Fext, ti, 0,InputData,NBodies) ! In the coupling between InWave and the WSC code, the weight is computed in InWave that's why WeightOrNot = 0.
                
        ! Computation of Th
        call get_Th(Mesh,Ecoulement,Fext,InertiaTerms,TimeRamp,InputData,NBodies)

    end subroutine Computation_M_CK_CT_Q_Th_MB
    
    subroutine API_get_S()
        
        ! Thue subroutine computes the matrix S.
        
        ! Allocation
        allocate(S(3,3,NBodies))
        
        ! Computation of S
        call get_S(Mesh,S,NBodies)
        
    end subroutine API_get_S
    
    subroutine API_get_DSDt()
        
        ! Thue subroutine computes the matrix DSDt.
        
        ! Allocation
        allocate(DSDt(3,3,NBodies))
        
        ! Computation of DSDt
        call get_DSDt(Mesh,DSDt,NBodies)
        
    end subroutine API_get_DSDt
    
    subroutine Deallocating_S()
    
        ! This subroutine deallocates the matrix S.
        
        if(allocated(S)) deallocate(S)
        
    end subroutine Deallocating_S
    
    subroutine Deallocating_DSDt()
    
        ! This subroutine deallocates the matrix DSDt.
        
        if(allocated(DSDt)) deallocate(DSDt)
        
    end subroutine Deallocating_DSDt
    
    subroutine API_Building_A
    
        ! This subroutine wraps the subroutine Building_A.
        
        call Building_A(Mesh,A,CS,CD,M,Nsl,NBody,size_NBody,N,Nsys,Nb,Nnodes,InputData,NBodies)
            
    end subroutine API_Building_A
    
    subroutine API_Building_B
    
        ! This subroutine wraps the subroutine Building_B.
    
        call Building_B(Mesh,Ecoulement,B,CD,Nsys,Nb,Nsl,Nnodes,N,InputData,NBody,size_NBody)
        
    end subroutine API_Building_B
    
    subroutine API_Initialization_solution()
    
        ! This subroutine wraps the subroutine Initialization_solution.
        
        call Initialization_solution(Mesh,Ecoulement,Nsys,Nb,Nbody,size_NBody,Sol,N,InputData)
    
    end subroutine API_Initialization_solution
    
    subroutine API_Preconditioning()
    
        ! This subroutine wraps the subroutine Preconditioning.
        
        ! Preconditioning of the linear system
        call Preconditioning(A,B,N)
    
    end subroutine API_Preconditioning
    
    subroutine API_Deallocating()
        
        ! This subroutine deallocates all the structures used in SystLin_FreeBody.
        
        call deallocate_CK_CT_Q_Th(Mesh)
        deallocate(S,DSDt,Fext,Inertiaterms,NBody)
        if(allocated(M)) deallocate(M)
        
    end subroutine API_Deallocating
    
    subroutine get_MBody(NumBody)
        
        integer,intent(in) :: NumBody  ! Body number
        
        ! This subroutine gives the displacement of the floater.
        
        if(NumBody.ge.Int_Body)then
            MBody = Mesh%Body(NumBody)%MBody(1:6)
        else
            print*,"NumBody in get_MBody is wrong!"
        end if
        
    end subroutine get_MBody
    
    subroutine get_CSolv(NumBody)
        
        integer,intent(in) :: NumBody  ! Body number
        
        ! This subroutine gives the position of the floater.
        
        if(NumBody.ge.Int_Body)then
            CSolv = Mesh%Body(NumBody)%CSolv(1:6)
        else
            print*,"NumBody in get_CSolv is wrong!"
        end if
        
    end subroutine get_CSolv
    
    subroutine get_VBody(NumBody)
        
        integer,intent(in) :: NumBody  ! Body number
        
        ! This subroutine gives the velocity of the floater.
        
        if(NumBody.ge.Int_Body)then
            VBody = Mesh%Body(NumBody)%VBody(1:6)
        else
            print*,"NumBody in get_VBody is wrong!"
        end if
        
    end subroutine get_VBody
    
    subroutine get_ABody(NumBody)
        
        integer,intent(in) :: NumBody  ! Body number
        
        ! This subroutine gives the acceleration of the floater.
        
        if(NumBody.ge.Int_Body)then
            ABody = Mesh%Body(NumBody)%ABody(1:6)
        else
            print*,"NumBody in get_ABody is wrong!"
        end if
        
    end subroutine get_ABody
    
    subroutine get_Nsys_Nb_Ndof()
                
        integer             :: nc       ! Loop parameter
        
        ! This subroutine gives the number of nodes in the mesh and in the meshes of each floater
        
        Nsys = Mesh%Nsys
        Nb = 0
        do nc = 1,Mesh%Nbody
            if (Mesh%Body(nc)%CMD(1)) then ! Bodies
                Nb = Nb + Mesh%Body(nc)%IndBody(3) - Mesh%Body(nc)%IndBody(1) + 1 ! Total number of nodes for all the bodies (not the tank).
            end if
        end do
        if(is_body)then
            Ndof = InputData%Ndll(1) ! Only for the first floater
        else
            Ndof = 0
        end if
        size_NBody = Mesh%NBody
        
    end subroutine get_Nsys_Nb_Ndof
    
    subroutine get_CMD(NumBody)
            
        integer,intent(in) :: NumBody  ! Body number
        
        ! This subroutine gives a value to CMD.
        
        CMD = Mesh%Body(NumBody)%CMD(1)
    
    end subroutine get_CMD
    
    subroutine get_CT_Body(NumBody)
        
        integer,intent(in)  :: Numbody   ! Body number
        
        ! This subroutine gives the CT coefficients of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
            allocate(CT(6,Mesh%Body(NumBody)%IndBody(3) - Mesh%Body(NumBody)%IndBody(1) + 1))
            CT = 0._RP
            CT = Mesh%Body(NumBody)%CT
        else
            print*,"NumBody in get_CT_Body is wrong!"
        end if
        
    end subroutine get_CT_Body
    
    subroutine Deallocate_CT
        
        ! This subroutine deallocates CT.
        
        deallocate(CT)
        
    end subroutine Deallocate_CT
    
    subroutine get_CK_Body(NumBody)
        
        integer,intent(in)  :: Numbody   ! Body number
        
        ! This subroutine gives the CK coefficients of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
            allocate(CK(Mesh%Body(NumBody)%IndBody(3) - Mesh%Body(NumBody)%IndBody(1) + 1,6))
            CK = 0._RP
            CK = Mesh%Body(NumBody)%CK
        else
            print*,"NumBody in get_CK_Body is wrong!"
        end if
        
    end subroutine get_CK_Body
    
    subroutine Deallocate_CK
        
        ! This subroutine deallocates CK.
        
        deallocate(CK)
        
    end subroutine Deallocate_CK
    
    subroutine get_Q_Body(NumBody)
        
        integer,intent(in)  :: Numbody   ! Body number
        
        ! This subroutine gives the Q coefficients of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
            allocate(Q(Mesh%Body(Numbody)%IndBody(3) - Mesh%Body(Numbody)%IndBody(1) + 1))
            Q = 0._RP
            Q = Mesh%Body(NumBody)%Q
        else
            print*,"NumBody in get_Q_Body is wrong!"
        end if
        
    end subroutine get_Q_Body
    
    subroutine Deallocate_Q
        
        ! This subroutine deallocates Q.
        
        deallocate(Q)
        
    end subroutine Deallocate_Q
    
    subroutine get_Th_Body(NumBody)
        
        integer,intent(in)  :: Numbody  ! Body number
        
        integer             :: k,j      ! Loop parameters
        
        ! This subroutine gives the Th coefficients of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
            allocate(Th(6))                        
            Th = 0._RP
            k = 1
            do j = 1,6
                if(InputData%dll_dir(j,NumBody-1))then
                    Th(j) = Mesh%Body(NumBody)%Th(k)
                    k = k + 1
                end if
            end do
        else
            print*,"NumBody in get_Th_Body is wrong!"
        end if
        
    end subroutine get_Th_Body
    
    subroutine Deallocate_Th
        
        ! This subroutine deallocates Th.
        
        deallocate(Th)
        
    end subroutine Deallocate_Th
    
    subroutine get_Fhydro_Body(NumBody,Method)
        
        integer,intent(in)      :: Numbody              ! Body number.
        integer,intent(in)      :: Method               ! Method of computation of Fhydro (1: direct, 2: with CT and Th).
        
        real(rp), dimension(3)  :: GPhip, GPhi0, Ds, M  !
        real(rp)                :: DPhiDt, GPhi2        ! DPhiDt and the square of the gradient of Phi.
        real(rp), dimension(3)  :: Vect_product_1       ! Vector.
        integer                 :: j,idir,k             ! Loop parameters.
                
        ! This subroutine gives the total hydrodynamic loads of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
                        
            allocate(Fhydro(6))
            Fhydro = 0._RP
            
            if(Method.eq.1)then ! Direct Method
                
                ! Time ramp
                call Ramp_Management(ti,TimeRamp)
                
                do j = Mesh%Body(NumBody)%IndBody(1),Mesh%Body(NumBody)%IndBody(3)
                
                    M = Mesh%Tnoeud(j)%Pnoeud
                    GPhip = Ecoulement%GPhi(:,j)%perturbation
                    GPhi0 = TimeRamp*Ecoulement%GPhi(:,j)%incident
                    DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
                    GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip) - dot_product(GPhip,GPhip)
                                        
                    ! Control surface
                    Ds = Mesh%Tnoeud(j)%Normale*Mesh%Tnoeud(j)%Aire
                    
                    ! Hydrodynamic Force
                    if(hydrostatique) then ! Included the hydrodstatics.
                        Fhydro(1:3) = Fhydro(1:3) + (DPhiDt + 0.5_rp*GPhi2 + g*M(3))*Ds
                    else ! Without the hydrostatics.
                        Fhydro(1:3) = Fhydro(1:3) + (DPhiDt + 0.5_rp*GPhi2)*Ds
                    end if
                                
                    ! Moment des forces hydrodynamiques, calculé au point de résolution de l'équation de mouvement CSolv
                    print*,"get_Fhydro_Body: maybe a modification should be done with GBody here."
                    call Computation_vect_product(M(1:3)-Mesh%Body(NumBody)%CSolv(1:3), Ds,Vect_product_1)
                    if(hydrostatique) then ! Included the hydrodstatics.
                        Fhydro(4:6) = Fhydro(4:6) + (DPhiDt + 0.5_rp*GPhi2 + g*M(3))*Vect_product_1
                    else ! Without the hydrostatics.
                        Fhydro(4:6) = Fhydro(4:6) + (DPhiDt + 0.5_rp*GPhi2)*Vect_product_1
                    end if
                    
                end do
                
                Fhydro = ro*Fhydro
                
                ! Symmetry
                if(Symmetry)then
                    Fhydro = 2._rp*Fhydro
                    Fhydro(2) = 0._rp
                    Fhydro(4) = 0._rp
                    Fhydro(6)  = 0._RP
                endif
                
            else if(Method.eq.2)then ! Computation from CT
                
                ! Ramp_Management is not called because if CT is used, that involves Computation_M_CK_CT_Q_Th or Computation_M_CK_CT_Q_Th_MB was run previsouly where Ramp_Management is called.
                
                k = 0
                do idir = 1,6
                    if(InputData%dll_dir(idir,NumBody-1))then
                        ! CT
                        Fhydro(idir) = Fhydro(idir) + dot_product(Mesh%Body(NumBody)%CT(idir,1:Mesh%Body(NumBody)%IndBody(3) - Mesh%Body(NumBody)%IndBody(1) + 1),Ecoulement%DPhiDt(Mesh%Body(NumBody)%IndBody(1):Mesh%Body(NumBody)%IndBody(3))%perturbation)
                        
                        ! Th
                        k = k + 1
                        Fhydro(idir) = Fhydro(idir) + Mesh%Body(NumBody)%Th(k) ! DPhiDt%Incident is in Th.
                        
                    end if
                end do
            else
                print*,"Method is not equal to 1 or 2. Fhydro = 0."
            end if
            
        else
            print*,"NumBody in get_Fhydro_Body is wrong!"
        end if
        
    end subroutine get_Fhydro_Body
    
    subroutine Deallocate_Fhydro
        
        ! This subroutine deallocates Fhydro.
        
        deallocate(Fhydro)
        
    end subroutine Deallocate_Fhydro
    
    subroutine get_DPhiDt_Body(NumBody)
        
        integer,intent(in)  :: Numbody  ! Body number
        
        integer             :: k,j      ! Loop parameters
        
        ! This subroutine gives DPhiDt_perturbation of the body NumBody.
        
        if(NumBody.ge.Int_Body)then
            allocate(DPhiDt(Mesh%Body(NumBody)%IndBody(3) - Mesh%Body(NumBody)%IndBody(1) + 1)) 
            k = 1
            do j = Mesh%Body(NumBody)%IndBody(1),Mesh%Body(NumBody)%IndBody(3)
                DPhiDt(k) = Ecoulement%DPhiDt(j)%perturbation
                k = k + 1
            end do
        else
            print*,"NumBody in get_DPhiDt_Body is wrong!"
        end if
        
    end subroutine get_DPhiDt_Body
    
    subroutine Deallocate_DPhiDt
        
        ! This subroutine deallocates DPhiDt.
        
        deallocate(DPhiDt)
        
    end subroutine Deallocate_DPhiDt
    
    subroutine get_InputData()
        
        ! This subroutine allows to have an access to the data of the structure InputData.
        
        ! ndll
        allocate(ndll(NBodies))
        ndll = InputData%ndll
        
        ! dll_dir
        allocate(dll_dir(6,NBodies))
        dll_dir = InputData%dll_dir
        
    end subroutine get_InputData
    
    subroutine Deallocate_InputData()
    
        ! This subroutine deallocates the variables allocated by get_InputData.
        
        deallocate(ndll)
        deallocate(dll_dir)
    
    end subroutine Deallocate_InputData
    
    subroutine compute_eJ0(NumBody)
        
        integer,intent(in)      :: NumBody          ! Body number
        
        real(rp),dimension(3,3) :: eOmega0,eR0      ! eOmega0 and eR0
        real(rp)                :: phi,theta,psi    ! Angles
        
        ! This subroutine computes the matrix eJ0 for the body NumBody.
        
        eJ0 = 0._RP
        eOmega0 = 0._RP
        eR0 = 0._RP
        
        phi = Mesh%Body(NumBody)%CSolv(4)
        theta = Mesh%Body(NumBody)%CSolv(5)
        psi = Mesh%Body(NumBody)%CSolv(6)
        
        ! eR0
        eR0(1, 1) = cos(psi)*cos(theta); eR0(1, 2) = -sin(psi)*cos(phi) + cos(psi)*sin(theta)*sin(phi); eR0(1, 3) = sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta)
	    eR0(2, 1) = sin(psi)*cos(theta); eR0(2, 2) = cos(psi)*cos(phi) + sin(phi)*sin(theta)*sin(psi); eR0(2, 3) = -cos(psi)*sin(phi) + sin(psi)*cos(phi)*sin(theta)
	    eR0(3, 1) = -sin(theta); eR0(3, 2) = cos(theta)*sin(phi); eR0(3, 3) = cos(theta)*cos(phi)
        
        ! eOmega0
        eOmega0(1, 1) = 1; eOmega0(1, 2) = sin(phi)*tan(theta); eOmega0(1, 3) = cos(phi)*tan(theta)
	    eOmega0(2, 1) = 0; eOmega0(2, 2) = cos(phi); eOmega0(2, 3) = -sin(phi)
	    eOmega0(3, 1) = 0; eOmega0(3, 2) = sin(phi) / cos(theta); eOmega0(3, 3) = cos(phi) / cos(theta)
        
        ! eJ0
        eJ0(1:3,1:3) = eR0(1:3,1:3)
        eJ0(4:6,4:6) = eOmega0(1:3,1:3)
        
    end subroutine compute_eJ0
    
    subroutine compute_eJ0p(NumBody)
        
        integer,intent(in)      :: NumBody          ! Body number
        
        real(rp),dimension(3,3) :: eOmega0p,eR0p    ! eOmega0 and eR0
        real(rp)                :: phi,theta,psi    ! Angles
        real(rp)                :: cphi,ctheta,cpsi ! Cosinus of the angles
        real(rp)                :: sphi,stheta,spsi ! Sinus of the angles
        real(rp)                :: ttheta           ! tangent of theta
        real(rp)                :: phip,thetap,psip ! Time differentiaiton of the angles
        
        ! This subroutine computes the matrix eJ0 for the body NumBody.
        
        print*,"The coefficients of this matrix are not verified. There are probably some mistakes."
        
        eJ0p = 0._RP
        eOmega0p = 0._RP
        eR0p = 0._RP
        
        phi = Mesh%Body(NumBody)%CSolv(4)
        theta = Mesh%Body(NumBody)%CSolv(5)
        psi = Mesh%Body(NumBody)%CSolv(6)
        
        cphi = cos(phi); ctheta = cos(theta); cpsi = cos(psi)
	    sphi = sin(phi); stheta = sin(theta); spsi = sin(psi)
        ttheta = tan(theta)
        
        phip = Mesh%Body(NumBody)%VBody(4)
        thetap = Mesh%Body(NumBody)%VBody(5)
        psip = Mesh%Body(NumBody)%VBody(6)
        
        ! eR0p
        eR0p(1, 1) = psip*spsi*ctheta - thetap*stheta*cpsi;
	    eR0p(1, 2) = -psip*cpsi*cphi + phip*spsi*cphi - psip*spsi*stheta*sphi + thetap*cpsi*ctheta*sphi + phip*cpsi*stheta*cphi;
	    eR0p(1, 3) = psip*cpsi*sphi + phip*spsi*cphi - psip*spsi*cphi*stheta - phip*cpsi*sphi*stheta + thetap*cpsi*cphi*ctheta;

	    eR0p(2, 1) = psip*cpsi*ctheta - thetap*spsi*stheta;
	    eR0p(2, 2) = -psip*spsi*cphi - phip*cpsi*sphi + phip*cphi*stheta*spsi + thetap*sphi*ctheta*spsi + psip*sphi*stheta*cpsi;
	    eR0p(2, 3) = psip*spsi*sphi - phip*cpsi*cphi + psip*cpsi*cphi*stheta - phip*spsi*sphi*stheta + thetap*spsi*cphi*ctheta;

	    eR0p(3, 1) = -thetap*ctheta;
	    eR0p(3, 2) = -thetap*stheta*sphi + phip*ctheta*cphi;
	    eR0p(3, 3) = -thetap*stheta*cphi - phip*ctheta*sphi;
        
        ! eOmega0p
        eOmega0p(1, 1) = 0; 
	    eOmega0p(1, 2) = psip*cphi*ttheta + (thetap*sphi / (ctheta*ctheta));
	    eOmega0p(1, 3) = -phip*sphi*ttheta + (thetap*cphi / (ctheta*ctheta));

	    eOmega0p(2, 1) = 0; 
	    eOmega0p(2, 2) = -phip*sphi;
	    eOmega0p(2, 3) = -phip*cphi;

	    eOmega0p(3, 1) = 0; 
	    eOmega0p(3, 2) = (phip*cphi*ctheta + thetap*stheta*sphi) / (ctheta*ctheta);
	    eOmega0p(3, 3) = (-phip*sphi*ctheta + thetap*stheta*cphi) / (ctheta*ctheta);
        
        ! eJ0p
        eJ0p(1:3,1:3) = eR0p(1:3,1:3)
        eJ0p(4:6,4:6) = eOmega0p(1:3,1:3)
        
    end subroutine compute_eJ0p
    
    subroutine compute_CDPhi
        
        integer :: j,k  ! Loop parameters
        
        ! This subroutine computes the quantity CD(:,FS)*Phi_t(FS).
        
        allocate(CDPhi(Nsys))
        CDPhi = 0._RP
        
        do j = 1,Nsys
            CDPhi(j) = 0._RP
            do k = Nsl(1),Nsl(2)
                CDPhi(j) = CDPhi(j) + CD(j,k)*Ecoulement%DPhiDt(k)%perturbation
            end do
        end do
        
    end subroutine compute_CDPhi
    
    subroutine deallocate_CDPhi()
    
        ! This subroutine deallocates CDPhi.
        
        deallocate(CDPhi)
    
    end subroutine deallocate_CDPhi
    
    subroutine API_Initialization_solution_MB
    
        ! This subroutine wraps the subroutine Initialization_solution_MB.
        
        allocate(Sol_WSC(Nsys + Nb))
        Sol_WSC = 0._RP
        call Initialization_solution_MB(Mesh,Ecoulement,Nsys,Nbody,size_NBody,Sol_WSC,N)
    
    end subroutine API_Initialization_solution_MB
    
    subroutine deallocate_Sol_WSC
    
        ! This subroutine deallocates Sol_WSC.
        
        deallocate(Sol_WSC)
    
    end subroutine deallocate_Sol_WSC
    
    subroutine allocating_A_B_Sol(Nsystem)
    
        integer,intent(in) :: Nsystem    ! Nsystem
        
        ! This subroutine allocates A, B and Sol.
        
        allocate(A(Nsystem,Nsystem),B(Nsystem),Sol(Nsystem))
    
    end subroutine allocating_A_B_Sol
        
    subroutine API_GMRES(nsystem)
        
        integer,intent(in) :: Nsystem    ! Nsystem
        
        ! This subroutine calls the solver GMRES.
        
        if (Solv.eq.0) then
            call GMRES(A, B, Sol, Nsystem,ierror)
        else  
            call LU(A, B, Sol, Nsystem, ierror)
        end if
        
        if(ierror/=0) goto 9999
        
        9999 continue
        if(ierror/=0)then
            write(*,90),ierror
        endif
        90      format('error #',i3,' : pb. in resolution of free motion.')
        
    end subroutine API_GMRES
    
    subroutine deallocate_A_B_Sol()
        
        ! This subroutine deallocate A, B and Sol.
        
        deallocate(A)
        deallocate(B)
        deallocate(Sol)
    
    end subroutine deallocate_A_B_Sol
    
    subroutine get_Inertia(NumBody)
        
        integer,intent(in) :: NumBody   ! Body number
        
        ! This subroutine gets the inertia matrix of the body NumBody.
        
        Inertia_mat = Mesh%Body(NumBody)%IBody
        
    end subroutine get_Inertia
    
    subroutine get_Gj(NumBody)

         integer,intent(in) :: NumBody   ! Body number
    
        ! This subroutine gives center of gravity of the floater.
    
        Gj = Mesh%Body(NumBody)%GBody(1:3)

    end subroutine get_Gj
    
    subroutine RK_Body(ABody,NumBody)
        
        real(rp),dimension(6),intent(in)    :: ABody    ! Acceleration of the floater.
        integer,intent(in)                  :: NumBody  ! Body number.
                
        ! This subroutine fills RKVBody and RKABody for the body NumBody.
        
        if(InputData%free_body(NumBody-1))then
            Mesh%Body(NumBody)%ABody = ABody
        end if
        
    end subroutine RK_Body
    
    subroutine get_Active_Body(NumBody)
        
        integer,intent(in) :: NumBody   ! Body number
        
        ! This subroutine gets Mesh%Body(NumBody)%Active.
        
        Active = Mesh%Body(NumBody)%Active
        
    end subroutine get_Active_Body
    
    subroutine API_Extract_Mesh_Body(fileMesh)
        
        character(len=150),intent(in) :: filemesh     ! Tecplot format input file.
        
        ! This subroutine wraps the subroutine Extract_Mesh_Body.
          
        call Extract_Mesh_Body(filemesh,Mesh,fgeom_vect,InputData)
        
    end subroutine API_Extract_Mesh_Body
    
    subroutine API_GeomInit()
    
        ! This subroutine wraps the subroutine GeomInit.
        
        ierror = 0
        call GeomInit(Mesh, fgeom_vect, 0._RP, InputData,.false., ierror)
        
    end subroutine API_GeomInit
    
    subroutine API_PlotMaill(fileMesh)
        
        character(len=150),intent(in) :: filemesh     ! Tecplot format input file.
        
        ! This subroutine wraps the subroutine PlotMaill.
        
        call PlotMaill(fileMesh,Mesh)
    
    end subroutine API_PlotMaill
    
    subroutine API_MeshVelBody_Vertex(tin,Alpha)
        
        real(rp),intent(in) :: tin      ! Current time.
        real(rp),intent(in) :: Alpha    ! k = L**Alpha.
    
        ! This subroutine wraps the subroutine MeshVelBody_Vertex.
        
        call MeshVelBody_Vertex(Mesh,tin,InputData,Alpha)
    
    end subroutine API_MeshVelBody_Vertex
    
    subroutine API_Geometry_Domain()
    
        ! This subroutine creates the geometry of the domain.
        
        ! Allocating of fgeom_vect
        call init_GeomVect(fgeom_vect)
                
        ! Creating of the inertial frame
        O   = [0._rp,0._rp,0._rp]
        VX%coord  = [1._rp,0._rp,0._rp] ! Previously only VX
        VY%coord  = [0._rp,1._rp,0._rp] ! Previously only VY
        VZ%coord  = [0._rp,0._rp,1._rp] ! Previously only VZ
        call assign_repere(1,O,VX%coord,VY%coord,VZ%coord,0._rp,0._rp,rep0)
        
        ! Initalization
        nface = 0 ; nline = 0
        HouleRF%index = 1       ! HouleRF%Index fixed equal to 1.
        
        ! Creating of the geomtry of the domain
        call create_geom(fdomaine,nface,nline,0,InputData,0) ! 0 is useless.
        fgeom_vect%nface_vect(1) = nface
        
        ! Intersection curves
        ierror = 0 ! Initialization of the error flag
        n_tab2 = 0 ! No body, no intersection curve.
        
        ! Flag to know if at least one floater pierces the free surface.
        !is_immerged = is_body .and. n_tab2==0
        is_immerged = .false.
        
        ! Update fgeom_vect%Active.
        AboveFS = .false.
        do jj = 1,NBodies           
            fgeom_vect%Active(jj) = .true. ! Body immerged or piercing.
        end do
    
    end subroutine API_Geometry_Domain
    
    subroutine API_Geometry_Body(BodyInt)
        
        integer :: BodyInt ! Index of the body (between 1 to NBodies).
    
        ! This subroutine creates the geometry of a body according to the method of CC.
        
        ! Creating of the geometry of the floaters
        if(is_body)then
            if(not(InputData%is_Mesh_File(BodyInt)) .and. Mesh_Type.eq.3)then
                nface_old = nface ! Saving the previous value of nface.
                call create_geom(fgeom_vect%geom(BodyInt),nface,nline,1,InputData,BodyInt)
                fgeom_vect%nface_vect(BodyInt+1) = nface - nface_old ! Number of faces for this geometry.
                call update_geom(rep0,fgeom_vect%geom(BodyInt),InputData%Position(:,2,BodyInt),InputData%Position(:,1,BodyInt))
            end if
        end if
        
    end subroutine API_Geometry_Body
    
    subroutine API_Mesh0d()
    
        ! This subroutine runs compute_mesh until mesh_inter (not run).
                
        ! Initilization of the transitional mesh.
        call init_mesh(Grid,nb_point,nb_arete,nb_tri)
        if(ierror/=0)then
            ierror = 210
            goto 9999
        endif
        
        ! Creation nouveau maillage.
        call NewMaillage(Mesh,Pointmax,NBodies+1) ! +1 for the tank.
        
        print*,"compute_mesh: Optimal point is searched with dx2 of the first floater."
        ierror = 0
    
        ! Following the mesh evolution with Tecplot.
        if (iwevol) then
            open(unit=ioevol,file="Advance_front.dat", iostat=ios)
            write(unit=ioevol,fmt='(150a)') 'VARIABLES = "X","Y","Z","Number","bf","nface"'
        end if
    
        iflag = 0
        
        ! Maillage 0D : ajout point de la geometrie dans le maillage
        if(iprint==1) print*,'** mesh0D ...'
    
        ! Tank
        call mesh0D(Grid,nb_point,fdomaine,t0,iflag,ierror) ! t = t0 !!!
    
        if(ierror/=0)then
            ierror = 220
            goto 9999
        end if
    
        ! Bodies
        if(is_body)then
            do jj = 1,NBodies
                if(fgeom_vect%Active(jj))then
                    call mesh0D(Grid,nb_point,fgeom_vect%geom(jj),t0,iflag,ierror) ! t = t0 !!!
                    if(ierror/=0)then
                        ierror = 221
                        goto 9999    
                    end if
                end if
            end do
        end if
    
        if(idebug>0)then
            if(is_body)then
                call write_debug(Grid,nb_point,nb_arete,nb_tri,t0,&
                    &                'Mesh_0D.dat','point_0D.dat','arrete_0D.dat',ierror) ! t = t0 !!!
            end if
        end if
        
        9999 continue
        if(ierror/=0)then
            write(*,99),ierror
        end if
        99 format('** error #',i3,' : pb. generation maillage')  
    
        if(ierror/=0 .and. idebug>0)then
            call write_debug(Grid,nb_point,nb_arete,nb_tri,t0,&
                &                'Mesh.dat','point.dat','arrete.dat',ierror) ! t = t0 !!!
        end if
        
    end subroutine API_Mesh0d
    
end module API_WSC