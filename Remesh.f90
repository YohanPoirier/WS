module Remesh_mod
use Constantes
use Structuresdonnees
use Parameters
use FonctionsCommunes
use Incident_mod
use GeomMesh
use GeomAxiSym
use GeomStruct
use GeomDef
use SolvNum
use PrePlot
use rbf_interp
use MeshModule
implicit none
contains

subroutine Remesh(Mesh, Mesh0, t, InputData,dtl, fgeom_vect, bcorrect)
    ! ----------------------------------------------------------------------------------
    ! Mise à jour de la position des noeuds du maillage
    ! - Elevation de la surface libre
    ! - Points doubles des bords extérieurs
    ! - Deplacement points du flotteur 
    ! Rq : La mise à jour de la surface libre ne prend pas en compte l'intersection
    ! et les déplacement dans le plan xOy
    ! - Déplacement des noeuds de la surface libre en fonction de leur vitesse dans le plan
    ! xOy. Puis replacement des points sur la surface libre (plus précis)
    ! - Déplacement des noeuds du flotteur en fonction de la vitesse de déplacement
    ! calculé précédement
    ! ----------------------------------------------------------------------------------
    
    !f2py integer*1, dimension(1000)            :: Mesh
    type(TMaillage),intent(inout)               :: Mesh                             ! Mesh.
    !f2py integer*1, dimension(1000)            :: Mesh0
    type(TMaillage),intent(in)                  :: Mesh0                            ! Copy at the first RK4 step of Mesh.
    real(rp),intent(in)                         :: t                                ! Current time.
    !f2py integer*1, dimension(1000)            :: InputData
    type(InputDataStruct),intent(inout)         :: InputData                        ! Input data.
    real(rp),intent(in),optional                :: dtl                              ! Time step of the RK4 step.
    !f2py integer*1, dimension(1000)            :: fgeom_vect
    type(type_GeomVect),intent(inout),optional  :: fgeom_vect                       ! Geometries.
    logical,optional                            :: bcorrect                         ! 
    
    !f2py integer*1, dimension(1000)            :: rep0
    type(repere3d)                              :: rep0                             ! Frame.
    character(len=50)                           :: filemaill, num                   ! Characters.
    integer                                     :: nc, j, k, j2, jn, jtemp, it,jj   ! Loop parameters.
    real(rp)                                    :: Eta0, h                          ! Free surface elevation.
    real(rp),dimension(3)                       :: M, N2                            ! Points.
    real(rp),dimension(3)                       :: O, VX, VY, VZ                    ! Vectors.
    logical                                     :: bool, boolc                      ! Booleans.
    integer,parameter                           :: itmax = 20                       ! Maximum number of iterations.
    integer                                     :: ierror                           ! Error flag.
    real(rp),dimension(3,3)                     :: Topsi, Tpsitheta, Tthetab        ! Elementary transformation matrices.
    real(rp),dimension(3,3)                     :: Tob                              ! eR0.
    real(rp),dimension(3,2)                     :: Trig2                            ! Cos and sin of theta.
    real(rp),dimension(3)                       :: GOt                              ! GOt0 in the Cartesian frame.
    
    ! This subroutine updates the position of the nodes and the geometries.
    
    ierror = 0
    
    O   = [0._rp,0._rp,0._rp]
    VX  = [1._rp,0._rp,0._rp]
    VY  = [0._rp,1._rp,0._rp]
    VZ  = [0._rp,0._rp,1._rp]
    call assign_repere(1,O,VX,VY,VZ,0._rp,0._rp,rep0)
    
    ! Test presence dtl
    if(present(dtl))then
        h = dtl
    else
        h = dt
    endif
    
    if(present(bcorrect))then
        boolc = bcorrect
    else
        boolc = .false.
    endif
    
    ! Modification of the FS.
    if(DeformFS)then
        do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            if(Mesh%Tnoeud(j)%Ndouble .eq. 0)then
                ! Pas de noeuds doubles : coeur de la surface libre
                Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud + Mesh%Tnoeud(j)%Velocity(1:3)*h
                call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
                Mesh%Tnoeud(j)%Pnoeud(3) = Eta0
            else
                ! Noeuds doubles : intersection surface libre / parois du domaine et surface libre / corps
                Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud 
                call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
                Mesh%Tnoeud(j)%Pnoeud(3) = Eta0
                do jn = 1,Mesh%Tnoeud(j)%Ndouble
                    j2 = Mesh%Tnoeud(j)%double(jn)
                    Mesh%Tnoeud(j2)%Pnoeud = Mesh%Tnoeud(j)%Pnoeud
                end do
            end if
        end do
    end if
    
    ! Modification of the Body.
    if(is_body)then
        jj = 1
        do nc = 1,Mesh%NBody
            if (nc.eq.1 .or. (nc.ge.Int_Body .and. InputData%DeformBody(jj))) then ! DeformBody is not defined for the tank.
                if (Mesh%Body(nc)%CMD(1)) then ! The tank is excluded anyway.
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                        Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud + Mesh%Tnoeud(j)%Velocity*h
                    end do
                    ! Ecrasement de la position des noeuds a l'intersection SL/corps, precedemment calculé dans ! Modification of the FS
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                        do jn = 1,Mesh%Tnoeud(j)%Ndouble
                            j2 = Mesh%Tnoeud(j)%double(jn)
                            Mesh%Tnoeud(j2)%Pnoeud = Mesh%Tnoeud(j)%Pnoeud
                        end do
                    end do
                    jj = jj + 1
                else
                    if(nc.ge.Int_Body)then 
                       jj = jj + 1
                    end if
                end if
            end if
        end do
    end if
    
    ! Repositionnement des points sur la surface libre.
    if(boolc)then
        do nc = Int_Body,Mesh%NBody
            do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                bool = .false.
                do k = 1,Mesh%Tnoeud(j)%Ndouble
                    jtemp = Mesh%Tnoeud(j)%double(k)
                    bool = bool .or. Mesh%Tnoeud(jtemp)%typeNoeud == 0
                end do
                if(bool)then
                    M = Mesh%Tnoeud(j)%Pnoeud
                    N2 = Mesh%Tnoeud(j)%Plocal(1:3,2,1)
                    call CEta0(M,t,Eta0)
                    it=1
                    do while (abs(M(3)-Eta0).gt.Epsilon .and. it.lt.itmax)
                        M = M + (Eta0-M(3))/N2(3)*N2
                        call CEta0(M,t,Eta0)
                        it=it+1
                    end do
                    Mesh%Tnoeud(j)%Pnoeud(1:3) = M(1:3)
                    do k = 1,Mesh%Tnoeud(j)%Ndouble
                        jtemp = Mesh%Tnoeud(j)%double(k)
                        Mesh%Tnoeud(jtemp)%Pnoeud = M(1:3)
                    end do
                end if
            end do
        end do
    end if
    
    ! Modification of the geometry of the cylinders in order to compute the intersection curve
    call Active_point_arete(fgeom_vect,InputData)
    
    jj = 1
    do nc = Int_Body,Mesh%NBody
        if(is_body .and. InputData%DeformBody(jj))then
            if(iFixPoint)then
                call update_position(Mesh%Body(nc)%GBody(1:3),Mesh0%Body(nc)%GBody(1:3),&
        &                     Mesh0%Body(nc)%CSolv(1:3),Mesh%Body(nc)%MBody(1:3),Mesh%Body(nc)%MBody(4:6),ierror)
            endif
            if (Mesh_type==2) then
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    call update_geom(rep0,fgeom_vect%geom(jj),Mesh%Body(nc)%CSolv(4:6),FixPointPos)
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    
                    GOt = 0._RP
                    
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
                    
                    ! GO at t in the Cartesian frame (O = center of the mesh = Position1).
                    GOt = matmul(Tob(1:3,1:3),InputData%GOt0(1:3,nc-1))
                                        
                    ! Updating the geometries from the motion of the floater.
                    call update_geom(rep0,fgeom_vect%geom(jj),Mesh%Body(nc)%CSolv(4:6),Mesh%Body(nc)%CSolv(1:3) + GOt) ! OeOm = OeG(t) + GOm = Position in the inertial frame of the center of gravity at t - Position of the center of the mesh wrt CoG (fixed during the simulation) = Position of the center of the mesh at t. nc-1 because there is not the tank data in InputData.
                end if
            end if
        endif
        jj = jj + 1
    end do
        
    ! Correction of the geometry of the cylinders
    call Correction_Geom(fgeom_vect,InputData)
    
    ! Calcul caractéristique geometrie aux noeuds
    call GeomInit(Mesh, fgeom_vect,t,InputData,.true.)
    
    if(idebug.gt.0) then
        Write( num, '( f0.4 )' ) t
        filemaill = 'Maillage_'//trim(num)//'.dat'
        call PlotMaill(filemaill, Mesh)
    end if
    
end subroutine Remesh  

subroutine MeshVel(Maillage,t,dtl,InputData)
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(inout)       :: Maillage     ! Mesh
    real(rp),intent(in)                 :: t            ! Current time
    real(rp),intent(in)                 :: dtl
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data
    
    integer                             :: j            ! Loop parameters
    
    ! This subroutine computes the velocity of the nodes.
    
    ! Calcul Vitesse des noeuds du flotteur
    call MeshVelBody(Maillage,t,dtl,InputData)
    
    ! Calcul Vitesse des noeuds de la surface libre
    if (DeformFS) then
        call MeshVelFS(Maillage,t,dtl)
    else
        do j = Maillage%FS%IndFS(1),Maillage%FS%IndFS(3)
            Maillage%Tnoeud(j)%Velocity = 0._RP
        end do
    end if
    
end subroutine MeshVel

subroutine MeshVelBody(Mesh,t,dtl,InputData)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                         ! Mesh.
    real(rp)                            :: t, dtl                       ! Current time and time step (actually 1._RP).
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                    ! Input data.
    
    integer                             :: j, nc, j2, jn                ! Loop parameters.
    real(rp), dimension(3)              :: OmegaMG                      ! Omega vect GM.
    real(rp),dimension(3,3)             :: A                            ! Matrix A of the linear system.
    real(rp),dimension(3)               :: B, Sol, V                    ! Vector B and the solution of the linear system.
    real(rp),dimension(3)               :: GPhi0                        ! GPhi of the incident waves.
    logical                             :: bool                         ! Boolean to know if there are twin nodes with the free surface mesh.
    integer                             :: ierror                       ! Error flag.
    integer                             :: Nk, k,jj                     ! Loop parameters.
    integer,parameter                   :: nPlocal = 1000               ! Maximum number of the local basis to save. Should be 1000*Mesh%NBody.
    integer,dimension(nPlocal)          :: indPlocal                    ! Index of the local basis to save.
    real(rp),dimension(3,3,nPlocal)     :: Plocal                       ! Local basis saved.
    real(rp),dimension(3)               :: vect_product_1,vect_product_2! Cross products.
            
    ! This subroutine computes the velocity of the nodes of the bodies based on a segment spring analogy.
    
    ierror = 0
    
    if (is_immerged) then ! All bodies are immerged, no intersection curve floater/free surface. Can appear even if the mesh strategy of CC is used.
        
        jj = 1
        do nc=1,Mesh%NBody            
            if (Mesh%Body(nc)%CMD(1)) then ! Body can be deformed.
                do j = Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                    if(iFixPoint)then
                        ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-FixPointPos,vect_product_1)
                    else
                        ! Motion equation solved at G (CSolv = (GBody,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-Mesh%Body(nc)%CSolv(1:3),vect_product_1)
                    end if
                    Mesh%Tnoeud(j)%Velocity = Mesh%Body(nc)%VBody(1:3) + vect_product_1
                end do
            else ! No deformation (tank).
                if(Mesh%Body(nc)%Active)then
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                        Mesh%Tnoeud(j)%Velocity = 0._rp
                    end do
                end if
            end if
        end do
        
    else ! At least one body pierces the free surface.
        
        ! Definition of the boundary conditions for every node before using a segment spring analogy.
        
        ! ---------------------------------------------------------------------------------------
        !   Sauvegarde base locale
        ! ---------------------------------------------------------------------------------------
        Nk = 0
        jj = 1
        do nc = Int_Body,Mesh%NBody            
            if(Mesh%Body(nc)%CMD(1))then ! No tank
                do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    if(Mesh%Tnoeud(j)%NDouble.gt.0 .or. (Symmetry .or. InputData%igtype(jj).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
                        Nk = Nk + 1
                        if(Nk.eq.nPlocal)then
                            print*,"MeshVelBody: nPlocal has reached its maximum value: ",nPlocal
                        end if
                        Plocal(1:3,1:3,Nk) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                        indPlocal(Nk) = j
                    end if
                end do
            end if
            jj = jj + 1
        end do
        
        ! ----------------------------------------------------------------------------------------
        ! Vitesse des noeuds du maillage
        ! ----------------------------------------------------------------------------------------
        jj = 1
        do nc = 1,Mesh%NBody
            if (Mesh%Body(nc)%CMD(1)) then ! No tank
                
                ! Linear velocity.
                if(iFixPoint)then
                    ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                    V = 0._RP
                else
                    ! Motion equation solved at G (CSolv = (GBody,angles)).
                    V = Mesh%Body(nc)%VBody(1:3)
                end if
                                
                do j = Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                    A = 0._RP
                    B = 0._RP
                    
                    ! Omega x GM
                    if(iFixPoint)then
                        ! Motion equation solved at A fixed in the inertial frame (Csolv = (FixPointPos,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3) - FixPointPos,OmegaMG) ! It is written MG but it is AM!
                    else
                        ! Motion equation solved at G (CSolv = (GBody,angles)).
                        call Computation_vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3) - Mesh%Body(nc)%CSolv(1:3),OmegaMG) ! It is written MG but it is GM!
                    end if
                    
                    bool = .false.
                    do jn = 1,Mesh%Tnoeud(j)%Ndouble
                        j2 = Mesh%Tnoeud(j)%double(jn)
                        bool = bool .or. Mesh%Tnoeud(j2)%typeNoeud == 0 ! Mesh%Tnoeud(j2)%typeNoeud == 0 -> Free surface
                    end do
                    
                    if(bool)then ! Noeuds Doubles avec surface libre.
                        B(1) = dot_product(V + OmegaMG, Mesh%Tnoeud(j)%Normale)
                        
                        j2 = Mesh%Tnoeud(j)%double(1)
                        if (Mesh%Tnoeud(j2)%typeNoeud /= 0) then
                            j2 = Mesh%Tnoeud(j)%double(2)
                            if(Mesh%Tnoeud(j2)%typeNoeud .ne. 0) print*,'MeshVelBody: error j2'
                        end if
                        call CGPhi0(Mesh%Tnoeud(j)%Pnoeud,t,GPhi0)
                        
                        B(2) = dot_product(GPhi0,Mesh%Tnoeud(j2)%Normale) ! Mesh%Tnoeud(j2)%Normale = Normal_FS
                                         
                        A(1,1:3) = Mesh%Tnoeud(j)%Normale ! Mesh%Tnoeud(j)%Normale = Normal_SM
                        A(2,1:3) = Mesh%Tnoeud(j2)%Normale ! Mesh%Tnoeud(j2)%Normale = Normal_FS
                        
                        if((symmetry .or. InputData%igtype(jj).eq.5).and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then ! Plane y = 0 in case of symmetry.
                            A(3,1:3) = [0._rp, 1._rp, 0._rp] 
                            B(3) = 0._rp
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = A(3,1:3)
                            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
                            call Computation_vect_product(A(1,1:3),A(3,1:3),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
                        elseif(Mesh%Tnoeud(j)%Ndouble .gt.1)then ! Ndouble > 1
                            if(j2.eq.Mesh%Tnoeud(j)%double(1))then
                                jn = Mesh%Tnoeud(j)%double(2)
                                A(3,1:3) = Mesh%Tnoeud(jn)%Normale
                            else
                                jn = Mesh%Tnoeud(j)%double(1)
                                A(3,1:3) = Mesh%Tnoeud(jn)%Normale
                            end if
                            B(3) = dot_product(V+OmegaMG,A(3,1:3))
                            call Computation_vect_product(A(1,1:3),A(2,1:3),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
                            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
                            call Computation_vect_product(A(1,1:3),Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
                        else ! Ndouble = 1
                            call Computation_vect_product(A(1,1:3),A(2,1:3),A(3,1:3))
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = A(3,1:3)
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
                            call Computation_vect_product(A(1,1:3),A(3,1:3),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
                            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
                            B(3) = 0._rp
                        end if
                        
                        call LU(A,B,Sol,3)
                        
                        Mesh%Tnoeud(j)%Velocity = Sol(1:3)
                        
                        do jn = 1,Mesh%Tnoeud(j)%Ndouble
                            j2 = Mesh%Tnoeud(j)%double(jn)
                            Mesh%Tnoeud(j2)%Velocity = Mesh%Tnoeud(j)%Velocity ! Twin nodes Floater - Free surface
                        end do  
                        
                    else ! Other nodes (not the nodes on the intersection curve).
                        
                        ! In case the norm of the normal is null.
                        if(norm2(Mesh%Tnoeud(j)%Normale) .lt. Epsilon)then
                            go to 10
                        end if
                        
                        ! Modification BaseLocale pour les noeuds de bords.                     
                        if ((Symmetry.or.InputData%igtype(jj).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2) then
                            ! Normal
                            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = Mesh%Tnoeud(j)%Normale
                            ! v
                            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = [0._rp, 1._rp, 0._rp]
                            ! u
                            call Computation_vect_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(j)%Plocal(1:3,3,1),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
                            
                        elseif(Mesh%Tnoeud(j)%Ndouble .ge. 1)then
                            
                            ! Nodes at the intersection of more than two interfaces (on the low circle of a cylinder for example)
                            j2 = Mesh%Tnoeud(j)%double(1)
                            
                            ! Local basis
                            ! Normal_SM
                            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = Mesh%Tnoeud(j)%Normale ! Mesh%Tnoeud(j)%Normale = Normal_SM and Mesh%Tnoeud(j2)%Normale = Normal_SL
                            ! u = Normal_SM vect Normal_SL
                            call Computation_vect_product(Mesh%Tnoeud(j)%Normale,Mesh%Tnoeud(j2)%Normale,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u is tangent at the intersection between the two surfaces (defined by the normal vectors)
                            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit vector
                            ! v = Normal_SM vect u
                            call Computation_vect_product(Mesh%Tnoeud(j)%Plocal(1:3,3,1),Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit vector
                            
                        end if
                                                
                        ! dw/dt
                        Mesh%Tnoeud(j)%Velocity = dot_product(V + OmegaMG,Mesh%Tnoeud(j)%Plocal(1:3,3,1))*Mesh%Tnoeud(j)%Plocal(1:3,3,1) ! ((xp + (Omega vect GM)).Normal)Normal
                        
                        if(Mesh%Tnoeud(j)%Ndouble.eq.1)then
                            ! Nodes at the intersection of TWO interfaces.
                            
                            ! dv/dt
                            Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + dot_product(V+OmegaMG,Mesh%Tnoeud(j)%Plocal(1:3,2,1))*Mesh%Tnoeud(j)%Plocal(1:3,2,1) ! ((xp + (Omega vect GM)).v)v
                                                            
                        elseif(Mesh%Tnoeud(j)%Ndouble.gt.1)then
                            ! Nodes at the intersection of THREE or more interfaces.
                            
                            ! dv/dt
                            Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + dot_product(V+OmegaMG,Mesh%Tnoeud(j)%Plocal(1:3,2,1))*Mesh%Tnoeud(j)%Plocal(1:3,2,1) ! ((xp + (Omega vect GM)).v)v
                            ! du/dt
                            Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + dot_product(V+OmegaMG,Mesh%Tnoeud(j)%Plocal(1:3,1,1))*Mesh%Tnoeud(j)%Plocal(1:3,1,1) ! ((xp + (Omega vect GM)).u)u
                            
                        end if
                        
                        if (Symmetry.or.InputData%igtype(jj).eq.5) then
                            if(abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2 .and. Mesh%Tnoeud(j)%Ndouble.eq.1)then ! For the nodes on the plan of symmetry.
                                Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + dot_product(V+OmegaMG,Mesh%Tnoeud(j)%Plocal(1:3,1,1))*Mesh%TNoeud(j)%Plocal(1:3,1,1)
                            end if
                        end if
                        
                        ! Updating twin nodes
                        do jn = 1,Mesh%Tnoeud(j)%Ndouble
                            j2 = Mesh%Tnoeud(j)%double(jn)
                            Mesh%Tnoeud(j2)%Velocity = Mesh%Tnoeud(j)%Velocity ! Twin nodes Floater - Floater
                        end do
                    end if
                                        
                end do
                
                ! Deformation of the mesh of the floater with the spring analogy method
                call SysLin_SpringBody(Mesh,dtl,nc,ierror,InputData,jj)
                
                jj = jj + 1
            else ! Tank and bodies which are above the free surface.
                
                if(Mesh%Body(nc)%Active)then
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                        Mesh%Tnoeud(j)%Velocity = 0._rp
                    end do
                end if
                
                ! If a floater is fixed, CMD = F but jj must be increased of 1.
                if(nc.ge.Int_Body)then 
                   jj = jj + 1
                end if
            end if
        end do
        
        10  continue
        
        ! ----------------------------------------------------------------------------
        ! Reaffectation base locale
        ! ----------------------------------------------------------------------------
        do j = 1,Nk
            k = indPlocal(j)
            Mesh%Tnoeud(k)%Plocal(1:3,1:3,1) = Plocal(1:3,1:3,j)
        end do
        
    end if
    
end subroutine MeshVelBody
    
subroutine SysLin_SpringBody(Mesh,dtl,nc,ierror,InputData,NumBody)
    ! ---------------------------------------------------------
    ! Optimisation du deplacement des noeuds sur la surface immergé
    ! - Calcul des vitesses de glissement en fonction des vitesses 
    !   normale aux parois
    ! ----------------------------------------------------------
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh                             ! Mesh.
    real(rp),intent(in)                 :: dtl                              ! 1._rp.
    integer,intent(in)                  :: nc                               ! Number of the body in Mesh.
    integer, intent(inout)              :: ierror                           ! Error flag.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                        ! Input data.
    integer,intent(in)                  :: NumBody                          ! Number of the body in InputData.

    integer                             :: j, j1, j2, ndim, NBD, NBD1, NBD2 ! Loop parameters, body indexes.
    integer                             :: k, k1, k2, jk                    ! Loop parameters.
    integer                             :: icase                            ! = 0 for single node, more  otherwise.
    real(rp),dimension(:,:),allocatable :: mat
    real(rp),dimension(:),allocatable   :: sm, Sol
    real(rp),dimension(3)               :: V, TDir, Vj, Vk
    real(rp)                            :: aj1,ak1,aj2,ak2,TRaideur
        
    ! This subroutine deforms the mesh of the floating bodies based on a segment spring analogy.
    
    NBD1 = Mesh%Body(nc)%IndBody(1)
    NBD2 = Mesh%Body(nc)%IndBody(3)
    NBD = NBD2-NBD1+1    
    ndim = 2*NBD
    
    allocate(mat(1:ndim,1:ndim))
    allocate(sm(1:ndim))
    allocate(Sol(1:ndim))

    mat(1:ndim,1:ndim) = 0._rp
    sm(1:ndim) = 0._rp
    Sol(1:ndim) = 0._rp
    
    ! Creation of the system.
    do j = NBD1,NBD2
                
        j1 = 2*(j-NBD1)+1 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        j2 = 2*(j-NBD1+1) ! 2* because Fi.ui = 0 and Fi.vi = 0.
        icase = Mesh%Tnoeud(j)%Ndouble
    
        if((Symmetry .or. InputData%igtype(NumBody).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
            icase = icase+1
        elseif(Symmetry .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2 .and. Mesh%Tnoeud(j)%Ndouble.eq.1)then ! Twin nodes on the symmetry plan.
            icase = icase+1
        endif

        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+Ldom(3)).lt.Epsilon)then
            icase = 2
        endif
        
        select case (icase)
        case(0) ! Simple node.
            k = 2
            do while (Mesh%Tnoeud(j)%Tvoisin(k,2).eq.1) ! Loop over the neighbours.
                
                jk = Mesh%Tnoeud(j)%TVoisin(k,1) ! 1st order of neighbourhood.
                k1 = 2*(jk-NBD1)+1
                k2 = 2*(jk-NBD1+1)
                
                ! Stiffness.
                if(Mesh%Tnoeud(jk)%Ndouble.eq.0)then
                    TRaideur = 1./norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**2 ! CC
                else
                    TRaideur = 1./norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**3 ! PYW
                end if
                                
                TDir = (Mesh%Tnoeud(j)%Pnoeud-Mesh%Tnoeud(jk)%Pnoeud)/norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud) ! TDir = nij = (pi-pj)/||pi-pj||.
                
                ! if(Mesh%Tnoeud(jk)%Ndouble.gt.0) TRaideur = 4.*TRaideur ! Pourquoi ?
                
                aj1 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),TDir) ! aj1 = ui.nij
                aj2 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),TDir) ! aj2 = vi.nij
                ak1 = dot_product(Mesh%Tnoeud(jk)%Plocal(1:3,1,1),TDir) ! ak1 = uj.nij
                ak2 = dot_product(Mesh%Tnoeud(jk)%Plocal(1:3,2,1),TDir) ! ak2 = vj.nij
                
                ! Fi.ui = 0
                mat(j1,j1) = mat(j1,j1) + TRaideur*aj1*aj1 ! kij*(ui.nij)**2
                mat(j1,j2) = mat(j1,j2) + TRaideur*aj1*aj2 ! kij*(ui.nij)*(vi.nij)
                mat(j1,k1) = -TRaideur*ak1*aj1
                mat(j1,k2) = -TRaideur*ak2*aj1
                
                ! Fi.vi = 0
                mat(j2,j2) = mat(j2,j2) + TRaideur*aj2*aj2 ! kij*(vi.nij)**2
                mat(j2,j1) = mat(j2,j1) + TRaideur*aj2*aj1 ! kij*(vi.nij)*(ui.nij)
                mat(j2,k1) = -TRaideur*ak1*aj2
                mat(j2,k2) = -TRaideur*ak2*aj2
                
                Vj = dot_product(Mesh%Tnoeud(j)%Velocity,Mesh%Tnoeud(j)%Plocal(1:3,3,1))*Mesh%Tnoeud(j)%Plocal(1:3,3,1) ! At least dw/dt*nij
                Vk = dot_product(Mesh%Tnoeud(jk)%Velocity,Mesh%Tnoeud(jk)%Plocal(1:3,3,1))*Mesh%Tnoeud(jk)%Plocal(1:3,3,1) ! At least dw/dt*nij
                sm(j1) = sm(j1) - TRaideur*dot_product(Vj*dtl,TDir)*aj1 &
        &                           + TRaideur*dot_product(Vk*dtl,TDir)*aj1
                sm(j2) = sm(j2) - TRaideur*dot_product(Vj*dtl,TDir)*aj2 &
        &                           + TRaideur*dot_product(Vk*dtl,TDir)*aj2
                k = k + 1
            end do
            
        case(1) ! Twin nodes
            k = 2
            do while (Mesh%Tnoeud(j)%Tvoisin(k,2).eq.1) 
                
                jk = Mesh%Tnoeud(j)%TVoisin(k,1)
                if(jk.gt.0 .and. Mesh%Tnoeud(jk)%NPanneau.ge.Int_Body)then ! In case a neighbour is on the FS.
                    
                    k1 = 2*(jk-NBD1)+1
                    k2 = 2*(jk-NBD1+1)
                    
                    ! Stiffness.
                    ! TRaideur = 1./norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**2 ! CC
                    TRaideur = 1./norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**3 ! PYW
                    
                    TDir = (Mesh%Tnoeud(j)%Pnoeud-Mesh%Tnoeud(jk)%Pnoeud)/norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud) ! Fij = kij*(dj-di)
                    
                    ! TRaideur = 4.*TRaideur ! Pourquoi ?
                    
                    aj1 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),TDir) ! ui.nij
                    aj2 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),TDir) ! vi.nij
                    ak1 = dot_product(Mesh%Tnoeud(jk)%Plocal(1:3,1,1),TDir) ! uj.nij
                    ak2 = dot_product(Mesh%Tnoeud(jk)%Plocal(1:3,2,1),TDir) ! vj.nij
                    mat(j1,j1) = mat(j1,j1) + TRaideur*aj1*aj1
                    mat(j1,j2) = mat(j1,j2) + TRaideur*aj1*aj2
                    mat(j1,k1) = -TRaideur*ak1*aj1
                    mat(j1,k2) = -TRaideur*ak2*aj1
                    Vj = dot_product(Mesh%Tnoeud(j)%Velocity,Mesh%Tnoeud(j)%Plocal(1:3,3,1))*Mesh%Tnoeud(j)%Plocal(1:3,3,1) ! Normal
                    Vk = dot_product(Mesh%Tnoeud(jk)%Velocity,Mesh%Tnoeud(jk)%Plocal(1:3,3,1))*Mesh%Tnoeud(jk)%Plocal(1:3,3,1) ! Normal
                    sm(j1) = sm(j1) - TRaideur*dot_product(Vj*dtl,TDir)*aj1 &
        &                           + TRaideur*dot_product(Vk*dtl,TDir)*aj1      
                endif
                k = k + 1
                
            enddo
            
            mat(j2,j2) = 1._rp
            sm(j2) = dot_product(Mesh%Tnoeud(j)%Velocity(1:3)*dtl,Mesh%Tnoeud(j)%Plocal(1:3,2,1))
        case(2)
            mat(j1,j1) = 1._rp
            mat(j2,j2) = 1._rp
            sm(j1) = dot_product(Mesh%Tnoeud(j)%Velocity(1:3)*dtl,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u
            sm(j2) = dot_product(Mesh%Tnoeud(j)%Velocity(1:3)*dtl,Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! v
        end select
        
    end do
    
    ! Inversion of the system.
    call LU(mat,sm,Sol,ndim,ierror)
    
    if(ierror/=0)then
        ierror = 100 ; goto 9999
    endif
    
    ! Distribution of the solution.
    do j = NBD1,NBD2
        j1 = 2*(j-NBD1)+1
        j2 = 2*(j-NBD1+1)
    
        icase = Mesh%Tnoeud(j)%Ndouble
    
        if((Symmetry .or. InputData%igtype(NumBody).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
            icase = icase+1
        elseif(Symmetry .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2 .and. Mesh%Tnoeud(j)%Ndouble.eq.1)then
            icase = icase + 1
        end if
        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+Ldom(3)).lt.Epsilon)then
            icase = 2
        end if
    
        if(icase == 0)then
            V = Sol(j1)/dtl * Mesh%Tnoeud(j)%Plocal(1:3,1,1) + Sol(j2)/dtl * Mesh%Tnoeud(j)%Plocal(1:3,2,1)
            Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + V
        elseif(icase == 1)then ! Twin nodes.
            V = Sol(j1)/dtl * Mesh%Tnoeud(j)%Plocal(1:3,1,1) ! Sol*u
            Mesh%Tnoeud(j)%Velocity = Mesh%Tnoeud(j)%Velocity + V
        endif
                
        ! Same velocity for the twin nodes.
        do k = 1,Mesh%Tnoeud(j)%Ndouble
            jk = Mesh%Tnoeud(j)%double(k)
            if(Mesh%Tnoeud(jk)%typeNoeud .eq. 0)then
                Mesh%Tnoeud(jk)%Velocity = Mesh%Tnoeud(j)%Velocity ! Twin nodes Floater - Free surface
            endif
        end do
    end do
    
    ! Deallocating.
    if(allocated(mat)) deallocate(mat)
    if(allocated(sm)) deallocate(sm)
    if(allocated(Sol)) deallocate(Sol)
    
    9999 continue
    if(ierror/=0)then
        write(*,90),ierror
    endif
    90 format('error #',i3,' : SUB_SystLin_SpringBody')

end subroutine SysLin_SpringBody

subroutine MeshVelBody_Vertex(Mesh,t,InputData,Alpha)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage)                     :: Mesh                         ! Mesh.
    real(rp)                            :: t                            ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                    ! Input data.
    real(rp),intent(in)                 :: Alpha                        ! Stiffness coefficient.
    
    integer                             :: j, nc, j2, jn                ! Loop parameters.
    real(rp), dimension(3)              :: OmegaMG                      ! Omega vect GM.
    real(rp),dimension(3,3)             :: A                            ! Matrix A of the linear system.
    real(rp),dimension(3)               :: B, Sol, V                    ! Vector B and the solution of the linear system.
    real(rp),dimension(3)               :: GPhi0                        ! GPhi of the incident waves.
    logical                             :: bool                         ! Boolean to know if there are twin nodes with the free surface mesh.
    integer                             :: ierror                       ! Error flag.
    integer                             :: Nk, k,jj                     ! Loop parameters.
    integer,parameter                   :: nPlocal = 1000               ! Maximum number of the local basis to save. Should be 1000*Mesh%NBody.
    integer,dimension(nPlocal)          :: indPlocal                    ! Index of the local basis to save.
    real(rp),dimension(3,3,nPlocal)     :: Plocal                       ! Local basis saved.
    real(rp),dimension(3)               :: vect_product_1,vect_product_2! Cross products.
            
    ! This subroutine computes the velocity of the nodes of the bodies based on a vertex spring analogy.
    
    ierror = 0
                
    ! Definition of the boundary conditions for every node before using a vertex spring analogy.
        
    ! ---------------------------------------------------------------------------------------
    !   Sauvegarde base locale
    ! ---------------------------------------------------------------------------------------
    Nk = 0
    jj = 1
    do nc = Int_Body,Mesh%NBody            
        if(Mesh%Body(nc)%CMD(1))then ! No tank
            do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                
                ! Twin nodes on the free surface.
                !call CEta0(Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Eta)
                !if(abs(Mesh%Tnoeud(j)%Pnoeud(3)-Eta).lt.1.0E-5)then ! +1 is hard coded here.
                if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+1._RP).lt.1.0E-5)then ! +1 is hard coded here. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    Nk = Nk + 1
                    if(Nk.eq.nPlocal)then
                        print*,"MeshVelBody: nPlocal has reached its maximum value: ",nPlocal
                    end if
                    Plocal(1:3,1:3,Nk) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                    indPlocal(Nk) = j
                    
                    ! n
                    Mesh%Tnoeud(j)%Plocal(1:3,3,1) = Mesh%Tnoeud(j)%Normale ! nSM
                    Mesh%Tnoeud(j)%Plocal(1:3,2,1) = [0._RP,0._RP,-1._RP] ! v = nSL
                    
                    ! u
                    call Computation_vect_product(Mesh%Tnoeud(j)%Normale,[0._RP,0._RP,-1._RP] ,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u = nSM x nSL
                    Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                    
                    ! v
                    call Computation_vect_product(Mesh%Tnoeud(j)%Normale,Mesh%Tnoeud(j)%Plocal(1:3,1,1) ,Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! v = nSM x (nSM x nSL)
                    Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                    
                end if
                
                !if(Mesh%Tnoeud(j)%NDouble.gt.0 .or. (Symmetry .or. InputData%igtype(jj).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
                !    Nk = Nk + 1
                !    if(Nk.eq.nPlocal)then
                !        print*,"MeshVelBody: nPlocal has reached its maximum value: ",nPlocal
                !    end if
                !    Plocal(1:3,1:3,Nk) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                !    indPlocal(Nk) = j
                !end if
            end do
        end if
        jj = jj + 1
    end do
        
    ! ----------------------------------------------------------------------------------------
    ! Vitesse des noeuds du maillage
    ! ----------------------------------------------------------------------------------------
    jj = 1
    do nc = 1,Mesh%NBody
        if (Mesh%Body(nc)%CMD(1)) then ! No tank
            !                    
            !do j = Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
            !    A = 0._RP
            !    
            !    bool = .false.
            !    do jn = 1,Mesh%Tnoeud(j)%Ndouble
            !        j2 = Mesh%Tnoeud(j)%double(jn)
            !        bool = bool .or. Mesh%Tnoeud(j2)%typeNoeud == 0 ! Mesh%Tnoeud(j2)%typeNoeud == 0 -> Free surface
            !    end do
            !    
            !    ! Twin nodes.
            !    if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+1).lt.1.0E-5)then ! +1 is hard coded here.
            !        bool = .true.
            !    end if
            !    
            !    if(bool)then ! Noeuds Doubles avec surface libre.
            !        
            !        !j2 = Mesh%Tnoeud(j)%double(1)
            !        !if (Mesh%Tnoeud(j2)%typeNoeud /= 0) then
            !        !    j2 = Mesh%Tnoeud(j)%double(2)
            !        !    if(Mesh%Tnoeud(j2)%typeNoeud .ne. 0) print*,'MeshVelBody: error j2'
            !        !end if
            !        
            !        call CGPhi0(Mesh%Tnoeud(j)%Pnoeud,t,GPhi0)
            !        
            !        A(1,1:3) = Mesh%Tnoeud(j)%Normale ! Mesh%Tnoeud(j)%Normale = Normal_SM
            !        !A(2,1:3) = Mesh%Tnoeud(j2)%Normale ! Mesh%Tnoeud(j2)%Normale = Normal_FS
            !        A(2,1:3) = [0._RP,0._RP,-1._RP] ! Mesh%Tnoeud(j2)%Normale = Normal_FS
            !            
            !        if((symmetry .or. InputData%igtype(jj).eq.5).and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
            !            A(3,1:3) = [0._rp, 1._rp, 0._rp] 
            !            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = A(3,1:3)
            !            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
            !            call Computation_vect_product(A(1,1:3),A(3,1:3),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
            !            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
            !        !elseif(Mesh%Tnoeud(j)%Ndouble .gt.1)then ! Ndouble > 1
            !        !    if(j2.eq.Mesh%Tnoeud(j)%double(1))then
            !        !        jn = Mesh%Tnoeud(j)%double(2)
            !        !        A(3,1:3) = Mesh%Tnoeud(jn)%Normale
            !        !    else
            !        !        jn = Mesh%Tnoeud(j)%double(1)
            !        !        A(3,1:3) = Mesh%Tnoeud(jn)%Normale
            !        !    end if
            !        !    call Computation_vect_product(A(1,1:3),A(2,1:3),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
            !        !    Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
            !        !    Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
            !        !    call Computation_vect_product(A(1,1:3),Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
            !        !    Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
            !        else ! Ndouble = 1
            !            call Computation_vect_product(A(1,1:3),A(2,1:3),A(3,1:3))
            !            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = A(3,1:3)
            !            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
            !            call Computation_vect_product(A(1,1:3),A(3,1:3),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
            !            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit norm.
            !            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = A(1,1:3)
            !        end if
            !            
            !    else ! Other nodes (not the nodes on the intersection curve).
            !            
            !        ! In case the norm of the normal is null.
            !        if(norm2(Mesh%Tnoeud(j)%Normale) .lt. Epsilon)then
            !            go to 10
            !        end if
            !            
            !        ! Modification BaseLocale pour les noeuds de bords.                     
            !        if ((Symmetry.or.InputData%igtype(jj).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2) then
            !            ! Normal
            !            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = Mesh%Tnoeud(j)%Normale
            !            ! v
            !            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = [0._rp, 1._rp, 0._rp]
            !            ! u
            !            call Computation_vect_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(j)%Plocal(1:3,3,1),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
            !            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit norm.
            !                
            !        elseif(Mesh%Tnoeud(j)%Ndouble .ge. 1)then
            !                
            !            ! Nodes at the intersection of more than two interfaces (on the low circle of a cylinder for example)
            !            !j2 = Mesh%Tnoeud(j)%double(1)
            !                
            !            ! Local basis
            !            ! Normal_SM
            !            Mesh%Tnoeud(j)%Plocal(1:3,3,1) = Mesh%Tnoeud(j)%Normale ! Mesh%Tnoeud(j)%Normale = Normal_SM and Mesh%Tnoeud(j2)%Normale = Normal_SL
            !            ! u = Normal_SM vect Normal_SL
            !            !call Computation_vect_product(Mesh%Tnoeud(j)%Normale,Mesh%Tnoeud(j2)%Normale,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u is tangent at the intersection between the two surfaces (defined by the normal vectors)
            !            call Computation_vect_product(Mesh%Tnoeud(j)%Normale,[0._RP,0._RP,-1._RP],Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u is tangent at the intersection between the two surfaces (defined by the normal vectors)
            !            Mesh%Tnoeud(j)%Plocal(1:3,1,1) = Mesh%Tnoeud(j)%Plocal(1:3,1,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! Unit vector
            !            ! v = Normal_SM vect u
            !            call Computation_vect_product(Mesh%Tnoeud(j)%Plocal(1:3,3,1),Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
            !            Mesh%Tnoeud(j)%Plocal(1:3,2,1) = Mesh%Tnoeud(j)%Plocal(1:3,2,1)/norm2(Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! Unit vector
            !                
            !        end if
            !            
            !    end if
            !    
            !end do
            
            ! Deformation of the mesh of the floater with the spring analogy method
            call SysLin_SpringBody_vertex(Mesh,nc,ierror,InputData,jj,Alpha)
                
            jj = jj + 1
        else ! Tank and bodies which are above the free surface.
            ! If a floater is fixed, CMD = F but jj must be increased of 1.
            if(nc.ge.Int_Body)then 
                jj = jj + 1
            end if
        end if
    end do
        
    10  continue
        
    ! ----------------------------------------------------------------------------
    ! Reaffectation base locale
    ! ----------------------------------------------------------------------------
    do j = 1,Nk
        k = indPlocal(j)
        Mesh%Tnoeud(k)%Plocal(1:3,1:3,1) = Plocal(1:3,1:3,j)
    end do
    
end subroutine MeshVelBody_Vertex

subroutine SysLin_SpringBody_Vertex(Mesh,nc,ierror,InputData,NumBody,Alpha)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh                             ! Mesh.
    integer,intent(in)                  :: nc                               ! Number of the body in Mesh.
    integer, intent(inout)              :: ierror                           ! Error flag.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                        ! Input data.
    integer,intent(in)                  :: NumBody                          ! Number of the body in InputData.
    real(rp),intent(in)                 :: Alpha                            ! Stiffness coefficient.

    integer                             :: j, j1, j2,j3,ndim, NBD,NBD1, NBD2 ! Loop parameters, body indexes.
    integer                             :: k, k1, k2, k3, jk                ! Loop parameters.
    integer                             :: icase                            ! = 0 for single node, more  otherwise.
    real(rp),dimension(:,:),allocatable :: mat
    real(rp),dimension(:),allocatable   :: sm, Sol
    real(rp),dimension(3)               :: X, TDir, xj, xk
    real(rp)                            :: aj1,ak1,aj2,aj3,ak2,ak3,TRaideur
    real(rp)                            :: bj1,bj2,bj3
    real(rp)                            :: eta                              ! Wave elevation.
        
    ! This subroutine deforms the mesh of the floating bodies based on a vertex spring analogy.
    
    NBD1 = Mesh%Body(nc)%IndBody(1)
    NBD2 = Mesh%Body(nc)%IndBody(3)
    NBD = NBD2-NBD1+1    
    ndim = 3*NBD
    
    allocate(mat(1:ndim,1:ndim))
    allocate(sm(1:ndim))
    allocate(Sol(1:ndim))

    mat(1:ndim,1:ndim) = 0._rp
    sm(1:ndim) = 0._rp
    Sol(1:ndim) = 0._rp
    
    ! Creation of the system.
    do j = NBD1,NBD2
        
        j1 = 3*(j-NBD1) + 1 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        j2 = 3*(j-NBD1) + 2 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        j3 = 3*(j-NBD1) + 3 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        
        icase = Mesh%Tnoeud(j)%Ndouble
    
        if((Symmetry .or. InputData%igtype(NumBody).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
            icase = icase + 1
        elseif(Symmetry .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2 .and. Mesh%Tnoeud(j)%Ndouble.eq.1)then ! Twin nodes on the symmetry plan.
            icase = icase + 1
        end if
        
        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+Ldom(3)).lt.Epsilon)then
            icase = 2
        end if
        
        ! Twin nodes on the free surface. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call CEta0(Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Eta)
        !if(abs(Mesh%Tnoeud(j)%Pnoeud(3)-Eta).lt.1.0E-5)then ! +1 is hard coded here.
        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+1._RP).lt.1.0E-5)then ! +1 is hard coded here.
            icase = 1
        end if
        
        select case (icase)
        case(0) ! Simple node.
            
            k = 2
            do while (Mesh%Tnoeud(j)%Tvoisin(k,2).eq.1) ! Loop over the 1st order neighbours.
                
                jk = Mesh%Tnoeud(j)%TVoisin(k,1) ! 1st order of neighbourhood.
                k1 = 3*(jk-NBD1) + 1
                k2 = 3*(jk-NBD1) + 2
                k3 = 3*(jk-NBD1) + 3
                
                ! Stiffness.
                TRaideur = norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**Alpha
                                
                aj1 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Plocal(1:3,1,1)) ! ui.uj
                aj2 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Plocal(1:3,2,1)) ! ui.vj
                aj3 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Normale) ! ui.nSMj
                
                ak1 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(jk)%Plocal(1:3,1,1)) ! vi.uj
                ak2 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(jk)%Plocal(1:3,2,1)) ! vi.vj
                ak3 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(jk)%Normale) ! vi.nSMj
                
                ! Fi.ui = 0
                mat(j1,j1) = mat(j1,j1) - TRaideur
                mat(j1,k1) = TRaideur*aj1
                mat(j1,k2) = TRaideur*aj2
                mat(j1,k3) = TRaideur*aj3
                
                ! Fi.vi = 0
                mat(j2,j2) = mat(j2,j2) - TRaideur
                mat(j2,k1) = TRaideur*ak1
                mat(j2,k2) = TRaideur*ak2
                mat(j2,k3) = TRaideur*ak3
                
                k = k + 1
            end do
            
            ! (xi-xiold).nSMi = 0
            mat(j3,j3) = 1._RP
            sm(j3) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale)
            
        case(1) ! Twin nodes
            
            k = 2
            do while (Mesh%Tnoeud(j)%Tvoisin(k,2).eq.1) ! Loop over the 1st order neighbours.
                
                jk = Mesh%Tnoeud(j)%TVoisin(k,1) ! 1st order of neighbourhood.
                k1 = 3*(jk-NBD1) + 1
                k2 = 3*(jk-NBD1) + 2
                k3 = 3*(jk-NBD1) + 3
                
                ! Stiffness.
                TRaideur = norm2(Mesh%Tnoeud(j)%Pnoeud - Mesh%Tnoeud(jk)%Pnoeud)**Alpha
                                
                aj1 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Plocal(1:3,1,1)) ! ui.uj
                aj2 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Plocal(1:3,2,1)) ! ui.vj
                aj3 = dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(jk)%Normale) ! ui.nSMj
                
                ! Fi.ui = 0
                mat(j1,j1) = mat(j1,j1) - TRaideur
                mat(j1,k1) = TRaideur*aj1
                mat(j1,k2) = TRaideur*aj2
                mat(j1,k3) = TRaideur*aj3
                
                k = k + 1
            end do
            
            !! (xi-xiold).ui = 0
            !mat(j1,j1) = 1._rp
            !sm(j1) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u
            
            ! (xi-xiold).nSL = 0
            mat(j2,j2) = 1._RP
            sm(j2) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                        
            ! (xi-xiold).nSMi = 0
            mat(j3,j3) = 1._RP
            sm(j3) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale)
            
        case(2)
            
            ! (xi-xiold).ui = 0
            mat(j1,j1) = 1._rp
            sm(j1) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,1,1)) ! u
            
            ! (xi-xiold).vi = 0
            mat(j2,j2) = 1._rp
            sm(j2) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,2,1)) ! v
            
            ! (xi-xiold).nSMi = 0
            mat(j3,j3) = 1._rp            
            sm(j3) = dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale) ! nSM
            
        end select
        
    end do
    
    ! Inversion of the system.
    call LU(mat,sm,Sol,ndim,ierror)
    
    if(ierror/=0)then
        ierror = 100 ; goto 9999
    endif
    
    ! Distribution of the solution.
    do j = NBD1,NBD2
        j1 = 3*(j-NBD1) + 1 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        j2 = 3*(j-NBD1) + 2 ! 2* because Fi.ui = 0 and Fi.vi = 0.
        j3 = 3*(j-NBD1) + 3 ! 2* because Fi.ui = 0 and Fi.vi = 0.
    
        icase = Mesh%Tnoeud(j)%Ndouble
        
        if((Symmetry .or. InputData%igtype(NumBody).eq.5) .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2)then
            icase = icase + 1
        elseif(Symmetry .and. abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon2 .and. Mesh%Tnoeud(j)%Ndouble.eq.1)then
            icase = icase + 1
        end if
        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+Ldom(3)).lt.Epsilon)then
            icase = 2
        end if
        
        ! Twin nodes. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+1).lt.1.0E-5 .or. Mesh%Tnoeud(j)%Pnoeud(3).le.-1.01)then ! +1 is hard coded here.
        if(abs(Mesh%Tnoeud(j)%Pnoeud(3)+1._RP).lt.1.0E-5)then ! +1 is hard coded here.
            icase = 1
        end if
        
        if(icase == 0)then
            X = Sol(j1) * Mesh%Tnoeud(j)%Plocal(1:3,1,1) + Sol(j2) * Mesh%Tnoeud(j)%Plocal(1:3,2,1) + Sol(j3) * Mesh%Tnoeud(j)%Normale
            Mesh%Tnoeud(j)%Pnoeud = X
        elseif(icase == 1)then ! Twin nodes.
            X = Sol(j1) * Mesh%Tnoeud(j)%Plocal(1:3,1,1) + Sol(j2) * Mesh%Tnoeud(j)%Plocal(1:3,2,1) + Sol(j3) * Mesh%Tnoeud(j)%Normale
            
            !print*,Sol(j3)-dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale)
            !print*,Sol(j2) - dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,2,1))
            !print*,Sol(j1) - dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,1,1))
            !print*,Mesh%Tnoeud(j)%Pnoeud(3)
            
            !print*,dot_product(Mesh%Tnoeud(j)%Plocal(1:3,1,1),Mesh%Tnoeud(j)%Plocal(1:3,2,1)),dot_product(Mesh%Tnoeud(j)%Plocal(1:3,2,1),Mesh%Tnoeud(j)%Plocal(1:3,3,1)),dot_product(Mesh%Tnoeud(j)%Plocal(1:3,3,1),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
            
            !X = Sol(j1) * dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,1,1))*Mesh%Tnoeud(j)%Plocal(1:3,1,1) + dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Plocal(1:3,2,1))*Mesh%Tnoeud(j)%Plocal(1:3,2,1) + dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale)*Mesh%Tnoeud(j)%Normale
            
            Mesh%Tnoeud(j)%Pnoeud = X
            !Mesh%Tnoeud(j)%Pnoeud = Mesh%Tnoeud(j)%Pnoeud
            
        endif
        
        ! Same position for the twin nodes.
        !do k = 1,Mesh%Tnoeud(j)%Ndouble
        !    jk = Mesh%Tnoeud(j)%double(k)
        !    if(Mesh%Tnoeud(jk)%typeNoeud .eq. 0)then
        !        Mesh%Tnoeud(jk)%Pnoeud = Mesh%Tnoeud(j)%Pnoeud ! Twin nodes Floater - Free surface
        !    endif
        !end do
    end do
    
    ! Deallocating.
    if(allocated(mat)) deallocate(mat)
    if(allocated(sm)) deallocate(sm)
    if(allocated(Sol)) deallocate(Sol)
    
    9999 continue
    if(ierror/=0)then
        write(*,90),ierror
    endif
    90 format('error #',i3,' : SUB_SystLin_SpringBody')

end subroutine SysLin_SpringBody_Vertex

subroutine MeshVelFS(Maillage,t,dtl)
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(inout)       :: Maillage     ! Mesh
    real(rp),intent(in)                 :: t            ! Current time
    real(rp),intent(in)                 :: dtl
    
    integer :: NFS1,NFS2,nd,ni,j,k,jk, jkn
    real(rp),dimension(:,:),allocatable :: xd,xi,fi
    real(rp),dimension(:),allocatable :: tab_corres,tab_data
    real(rp),dimension(:),allocatable :: dlx,dly,w
    real(rp),dimension(3,3) :: A
    real(rp),dimension(3) :: B, Sol, M, GPhi0
    real(rp),dimension(:,:),allocatable :: mat  ! A in the problem of the coefficients of the radial basis approximation.
    integer                             :: ndim ! Dimension of the linear system
    integer                             :: nCoord
        
    ! This subroutine computes the velocity of the nodes of the free surface.
    
    NFS1 = Maillage%FS%IndFS(1)
    NFS2 = Maillage%FS%IndFS(3)
    
    ! Number of points of controle (nd) and mobile points (ni).
    nd = 0 ; ni = 0
    do j = NFS1,NFS2
        if(Maillage%Tnoeud(j)%control_point)then ! Controle points : boundary points (free surface and intersection points between the FS and the floates).
            nd = nd+1
        elseif(Maillage%Tnoeud(j)%mobility)then ! Free surface points inside the surface delimited by the controle points + points of the floaters.
            ni = ni+1
        end if
    end do
    
    ! Dimension of the linear system.
    if(Htype.eq.0)then
        nCoord = 2 ! 2 because the 2 coordinates (x,y) are used to compute the interpolation (the free surface is on z = 0 so the vertical coordinate is useless).
    else 
        nCoord = 3 ! 3 because the 3 coordinates (x,y,z) are used to compute the interpolation.
    end if
    ndim = nd + 1 + nCoord 
    
    ! Allocating.
    allocate(xd(nCoord,nd)) 
    allocate(xi(nCoord,ni))
    allocate(fi(2,ni)) ! 2 because the interpolation is only applied for the displacement in x and y.
    allocate(tab_data(nd)) ! Link between the number of controle points between [1,NFS] and [1,ni].
    allocate(tab_corres(ni)) ! Link between the number of mobile points between [1,NFS] and [1,ni].
    
    ! Initialisation des tableaux a zero.
    fi(1:2,1:ni) = 0._rp
    
    ! Recherche des points de controle.
    nd = 0 ; ni = 0
    do j = NFS1,NFS2
        if(Maillage%Tnoeud(j)%control_point)then ! Controle points : boundary points (free surface and intersection points between the FS and the floates).
            nd = nd+1
            do k = 1,nCoord
                xd(k,nd) = Maillage%Tnoeud(j)%Pnoeud(k)
            end do
            tab_data(nd) = j
        elseif(Maillage%Tnoeud(j)%mobility)then ! Free surface points inside the surface delimited by the controle points + points of the floaters.
            ni = ni+1
            do k = 1,nCoord
                xi(k,ni) = Maillage%Tnoeud(j)%Pnoeud(k)
            end do
            tab_corres(ni) = j
        else ! Free surface points outside the surface delimited by the controle points.
            Maillage%Tnoeud(j)%Velocity(1:2) = 0._rp
            call CDEta0Dt(Maillage%Tnoeud(j)%Pnoeud(1:3),t,Maillage%Tnoeud(j)%Velocity(3)) ! Only vertical velocity.
        end if
    end do
    
    ! RBF method.
    if(not(is_immerged))then ! At least one body is piercing.
        allocate(w(ndim))
        allocate(dlx(nd),dly(nd))
        allocate(mat(ndim,ndim))
        dlx(1:nd) = 0._rp
        dly(1:nd) = 0._rp
        mat = 0._RP
        
        ! Definition des vitesses de deplacement sur la surface libre.
         do j = 1,nd
            k = tab_data(j)
            dlx(j) = Maillage%Tnoeud(k)%Velocity(1)*dtl
            dly(j) = Maillage%Tnoeud(k)%Velocity(2)*dtl
        end do
        
        ! Calcul des deplacements des noeuds de la surface libre.
        
        call rbf_computing_A(nCoord,nd,xd,phi_log,mat,ndim)
        
        call rbf_weight(nCoord, nd, dlx(1:nd), w,mat,ndim) ! Alpha and beta for the x-velocities.
        call rbf_funct(nCoord, nd, xd, phi_log, w, ni, xi, fi(1,1:ni))
        
        call rbf_weight(nCoord, nd, dly(1:nd), w,mat,ndim) ! Alpha and beta for the y-velocities.
        call rbf_funct(nCoord, nd, xd, phi_log, w, ni, xi, fi(2,1:ni))
        
        ! Mise a jour de la vitesse des noeuds en fonction des deplacements.
        A(1,1:3) = [1._rp,0._rp,0._rp]
        A(2,1:3) = [0._rp,1._rp,0._rp]
        do j = 1,ni
            k = tab_corres(j)       
            A(3,1:3) = Maillage%Tnoeud(k)%Normale
            B(1:2) = fi(1:2,j)/dtl
            if(Symmetry .and. abs(Maillage%Tnoeud(k)%Pnoeud(2)).lt.Epsilon2)then
                B(2) = 0._rp
            end if
            M = Maillage%Tnoeud(k)%Pnoeud
            call CGPhi0(M,t,GPhi0)
            B(3) = dot_product(GPhi0,A(3,1:3))
            call LU(A,B,Sol,3)
            Maillage%Tnoeud(k)%Velocity(1:3) = Sol            
        end do
        
        ! Mise a jour de la vitesse des noeuds de controle.
        do j = 1,nd
            k = tab_data(j)
            A(3,1:3) = Maillage%Tnoeud(k)%Normale
            B(1:2) = [dlx(j),dly(j)]
            M = Maillage%Tnoeud(k)%Pnoeud
            if(Maillage%Tnoeud(k)%Ndouble.ne.0) then ! Twin nodes.
                do jk = 1,Maillage%Tnoeud(k)%Ndouble
                    jkn = Maillage%Tnoeud(k)%double(jk)
                    ! Ne pas traiter les noeuds à l'intersection avec le corps.
                    if (Maillage%Tnoeud(jkn)%NPanneau.lt.Int_Body) then ! Int_Body is for the first floater.
                        call CGPhi0(M,t,GPhi0)
                        B(3) = dot_product(GPhi0,A(3,1:3))
                        call LU(A,B,Sol,3)
                        Maillage%Tnoeud(k)%Velocity(1:3) = [0._rp,0._rp,Sol(3)]
                        jkn = Maillage%Tnoeud(k)%double(jk)
                        Maillage%Tnoeud(jkn)%Velocity = Maillage%Tnoeud(k)%Velocity
                    end if
                end do
            else ! Ndouble = 0.
                call CGPhi0(M,t,GPhi0)
                B(3) = dot_product(GPhi0,A(3,1:3))
                call LU(A,B,Sol,3)
                Maillage%Tnoeud(k)%Velocity(1:3) = [0._rp,0._rp,Sol(3)]
            end if
        end do
    end if
    
    ! Deallocation des tableaux.
    if(allocated(xd)) deallocate(xd)
    if(allocated(xi)) deallocate(xi)
    if(allocated(fi)) deallocate(fi)
    if(allocated(tab_corres)) deallocate(tab_corres)
    if(allocated(tab_data)) deallocate(tab_data)
    if(allocated(dlx)) deallocate(dlx)
    if(allocated(dly)) deallocate(dly)
    if(allocated(w)) deallocate(w)
    if(allocated(mat)) deallocate(mat)
    
end subroutine MeshVelFS

subroutine adjust_waterline2(Mesh,fgeom_vect,CrossingFS,t,ierror,InputData,tab2,n_tab2,n_tab)
    
    !f2py integer*1, dimension(1000)                :: Mesh
    type(TMaillage),intent(inout)                   :: Mesh                 ! Mesh
    !f2py integer*1, dimension(1000)                :: fgeom_vect
    type(type_GeomVect), intent(inout)              :: fgeom_vect           ! Geometry of the floater (not the domain)
    logical,intent(inout)                           :: CrossingFS           ! A body crossed the free surface (True) or not (False)
    real(rp),intent(in)                             :: t                    ! Current time
    integer,intent(inout)                           :: ierror               ! Error flag
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(inout)             :: InputData            ! Input data
    !f2py integer*1, dimension(1000)                :: tab2
    type(chaine_point_pt),dimension(100),optional   :: tab2                 ! Table of the pointers toward the new intersection points
    integer,optional                                :: n_tab2,n_tab         ! New number of intersection curves and lines.
    
    !f2py integer*1, dimension(1000)                :: Noeud        
    type(point)                                     :: Noeud
    integer                                         :: NBD1, NBD2
    integer                                         :: j,k,jk,n,jtemp,jn    ! Loop parameters.
    logical                                         :: bool
    real(rp)                                        :: dist, dist_min
    real(rp),dimension(3)                           :: M, M0
    integer,dimension(300)                          :: tab                  ! Table of the nodes of the intersection curve in the mesh.
    !f2py integer*1, dimension(1000)                :: ptr0,ptr1
    type(chaine_point),pointer                      :: ptr0,ptr1
    integer                                         :: nc                   ! Loop parameter.
    integer                                         :: n_tab2_old           ! Old number of intersection curves.
    logical                                         :: badjust              ! Boolean to know some bodies present in the structure Tmaillage pierce the free surface.
        
    ! This subroutine updates the position of the intersection curve.
    
    ! Initialization
    ierror = 0
    CrossingFS = .false.
    badjust = .false.
    n_tab2_old = n_tab2
        
    ! This following algorithm searches the nodes of the structure TMaillage which are twin with both a floater and the free surfarce.
    ! If a floater goes across the free surface, it does not see that.
    n = 0
    do nc = Int_Body,Mesh%NBody
        if(Mesh%Body(nc)%Active)then
            NBD1 = Mesh%Body(nc)%IndBody(1)
            NBD2 = Mesh%Body(nc)%IndBody(3)
            do j = NBD1,NBD2
                bool = .false.
                do k = 1,Mesh%Tnoeud(j)%Ndouble
                    jtemp = Mesh%Tnoeud(j)%double(k)
                    bool = bool .or. Mesh%Tnoeud(jtemp)%typeNoeud == 0 ! Twin nodes on the free surface
                end do
                if(bool)then
                    n = n+1
                    tab(n) = j
                end if
            end do
        end if
    end do
    
    ! Computation of the intersection curve.
    if(is_body)then
        call compute_intersection(t,fgeom_vect,tab2,n_tab2,5,ierror,InputData,inputData%dx2(1),n_tab)
    end if
    
    ! Either twin nodes Floater-FS are detected (with old mesh) or a new number of intersection curves is found (futur mesh).
    if(n>0 .or. n_tab2.ne.n_tab2_old)then
        
        ! Some bodies present in the structure Tmaillage pierce the free surface.
        badjust = .true.
        
        ! Flag to know if all the floaters are immerged or not.
        is_immerged = is_body .and. n_tab2==0 ! Il faut dire que les corps sont en dessous de la SF et non au dessus aussi !
        
    end if
    
    if(n_tab2.ne.n_tab2_old)then
        CrossingFS = .true.
    end if
    
    if(ierror.ne.0)then
        go to 9999
	end if
    
    1   continue
        
    ! Updating the nodes on the intersection curve, even if there is no remeshing after calling adjust_waterline2.
    if(not(badjust))then
        do j=1,n
            jk = tab(j)
            M0 = Mesh%Tnoeud(jk)%Pnoeud
            M = M0
            dist_min = 999.
            do k=1,n_tab2
                ptr0 => tab2(k)%pt
                ptr1 => ptr0
                Noeud = ptr1%val
                dist = norm2(Noeud%coord-M0)
                if(dist<dist_min)then
                    M = Noeud%coord
                    dist_min = dist
                endif
                ptr1 => ptr1%suiv
                do while(associated(ptr1).and. .not.associated(ptr1,ptr0))
                    Noeud = ptr1%val
                    dist = norm2(Noeud%coord-M0)
                    if(dist<dist_min)then
                        M = Noeud%coord
                        dist_min = dist
                    endif
                    ptr1 => ptr1%suiv
                enddo
            enddo
            Mesh%Tnoeud(jk)%Pnoeud  = M
            do k=1,Mesh%Tnoeud(jk)%Ndouble
                jn = Mesh%Tnoeud(jk)%double(k)
                Mesh%Tnoeud(jn)%Pnoeud = M
            enddo
        enddo
    end if
    
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('** error #',i3,' : cannot adjust waterline2.')

end subroutine adjust_waterline2

subroutine CheckMesh3D(Maillage,ti,boolBodies,boolFS,InputData)
    
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(in)          :: Maillage                     ! Mesh.
    real(rp),intent(in)                 :: ti                           ! Current time.
    logical,intent(inout)               :: boolBodies                   ! bool = false: no remeshing of the bodies, bool = true: remeshing of the bodies.
    logical,intent(inout)               :: boolFS                       ! bool = false: no remeshing of the free surface, bool = true: remeshing of the free surface.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData                    ! Input data.
    
    integer                             :: j,NF1,NFT,NFB1,NFBT          ! Node parameters.
    integer,dimension(3)                :: id                           ! Indexes of the 3 nodes/vertexes of a panel.
    real(rp),dimension(3,3)             :: PNoeud                       ! Coordinates of the three vertexes of a panel.
    real(rp),dimension(:),allocatable   :: fsize                        ! Function about the size of a panel (between 0 and 1).
    real(rp),dimension(:),allocatable   :: fshape                       ! Function about the shape of a pane (between 0 and 1).
    real(rp)                            :: area                         ! Area of a panel.
    real(rp),dimension(3)               :: L1,L2,L3                     ! Edge vectors of a panel.
    real(rp)                            :: d1, d2, d3                   ! Size of each edge of a panel.
    real(rp),parameter                  :: c_area = 0.25_rp*sqrt(3._rp) ! Area of a equilateral area of which the edge size is dx2 / dx2^2.
    integer                             :: nc                           ! Loop parameter.
    real(rp)                            :: dx                           ! Reference length of a panel.
    
    ! This subroutine checks the quality of the mesh.
    
    boolBodies = .false.
    boolFS = .false.
    NF1 = 1
    NFT = Maillage%Nfacette
    
    ! Allocating.
    allocate(fsize(NFT))
    allocate(fshape(NFT))
    
    ! Initialization.
    fsize(:) = -1._RP
    fshape(:) = -1._RP
    
    do j = NF1,NFT
        
        id = Maillage%Tfacette(j)%Tnoeud(1:3)
        PNoeud(1:3,1) = Maillage%Tnoeud(id(1))%PNoeud ! 1st vertex of the panel
        PNoeud(1:3,2) = Maillage%Tnoeud(id(2))%Pnoeud ! 2nd vertex of the panel
        PNoeud(1:3,3) = Maillage%Tnoeud(id(3))%Pnoeud ! 3rd vertex of the panel
        area = Maillage%Tfacette(j)%Aire
        
        ! Edge.
        L1 = PNoeud(1:3,3) - PNoeud(1:3,2) ! 1st edge of the panel
        L2 = PNoeud(1:3,1) - PNoeud(1:3,3) ! 2nd edge of the panel
        L3 = PNoeud(1:3,2) - PNoeud(1:3,1) ! 3rd edge of the panel
        d1 = norm2(L1) ! Size of the 1st edge
        d2 = norm2(L2) ! Size of the 2nd edge
        d3 = norm2(L3) ! Size of the 3rd edge
                    
        ! Metric shape.
	    fshape(j) = d3*d3+d2*d2+dot_product(L2,L3)
        fshape(j) = (2._rp*sqrt(3._rp)*area)/fshape(j)
        
        if(Maillage%Tfacette(j)%NPanneau .ge. Int_Body)then ! Bodies
            
            ! Reference length for each panel.
            dx = InputData%dx2(Maillage%Tfacette(j)%NPanneau-1) ! -1 because NPanneau starts at 2 for the bodies and %dx2 starts at 1.
            
            ! Metric fsize.
            fsize(j) = area/(c_area*dx*dx) ! Area of the panel / Area of a reference panel
            
        end if
                
    end do
    
    ! Free surface.
    do j = Maillage%FS%IndFS(2),Maillage%FS%IndFS(4)
        
        ! fshape.
        if (fshape(j).lt.0.25_rp) then
            id = Maillage%Tfacette(j)%Tnoeud(1:3)
            print*,""
            print*,"Free surface"
            print*, 'Shape :', j, fshape(j)
            print*,"Node 1:"
            print*, Maillage%Tnoeud(id(1))%PNoeud
            print*,"Node 2:"
            print*, Maillage%Tnoeud(id(2))%PNoeud
            print*,"Node 3:"
            print*, Maillage%Tnoeud(id(3))%PNoeud
            boolFS = .true.
            exit
        end if
        
        ! fsize.
        ! No fsize for the free surface because there is no reference size for each panel.
        
    end do
    
    ! Bodies.
    do nc = Int_Body,Maillage%NBody
        if(Maillage%Body(nc)%Active)then
            NFB1 = Maillage%Body(nc)%IndBody(2)
            NFBT = Maillage%Body(nc)%IndBody(4)
            do j = NFB1,NFBT
            
                ! fshape.
                if (fshape(j).lt.0.25_rp) then
                    id = Maillage%Tfacette(j)%Tnoeud(1:3)
                    print*,""
                    print*,"Body",nc
                    print*,'Shape :', j, fshape(j)
                    print*,"Node 1:"
                    print*, Maillage%Tnoeud(id(1))%PNoeud
                    print*,"Node 2:"
                    print*, Maillage%Tnoeud(id(2))%PNoeud
                    print*,"Node 3:"
                    print*, Maillage%Tnoeud(id(3))%PNoeud
                    boolBodies = .true.
                    exit
                end if
            
                ! fsize.
                if (fsize(j).lt.0.25_rp .or. fsize(j).gt.4._rp)then
                    id = Maillage%Tfacette(j)%Tnoeud(1:3)
                    print*,"Body",nc
                    print*,'Size :', j, fsize(j)
                    print*,"Node 1:"
                    print*, Maillage%Tnoeud(id(1))%PNoeud
                    print*,"Node 2:"
                    print*, Maillage%Tnoeud(id(2))%PNoeud
                    print*,"Node 3:"
                    print*, Maillage%Tnoeud(id(3))%PNoeud
                    boolBodies = .true.
                    exit
                end if
            
            end do
        end if
    end do
        
    if(iwsize)then
        write(iosize,'(4f8.2)') ti,minval(fshape(NFB1:NFBT)),minval(fsize(NFB1:NFBT)),maxval(fsize(NFB1:NFBT))
    endif
    
    if(allocated(fsize)) deallocate(fsize)
    if(allocated(fshape)) deallocate(fshape)
    
end subroutine CheckMesh3D

subroutine InitMesh(Mesh)
! Parameter
!f2py integer*1, dimension(1000) :: Mesh
type(TMaillage) :: Mesh
! Locals
integer :: j, nc
real(rp) :: theta
real(rp), dimension(3) :: InitRot, A, M
real(rp), dimension(3,3) :: Stheta
! Begin
!do nc = 2,Mesh%NBody
print*,"InitMesh: No MB in this subroutine."
pause
nc = Int_Body
    InitRot(1:3) = Mesh%Body(nc)%CSolv(4:6)
    theta = Mesh%Body(nc)%CSolv(5)
    Stheta = reshape([cos(theta), 0._RP, sin(theta), 0._RP, 1._RP, 0._RP, -sin(theta), 0._RP, cos(theta)],(/3,3/))
    A = Mesh%Body(nc)%GBody(1:3)
    if (iFixPoint) A = FixPointPos
    do j = Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)
        M = Mesh%Tnoeud(j)%Pnoeud - A
        Mesh%Tnoeud(j)%Pnoeud = A + matmul(Stheta,M)
    end do
    if (iFixPoint) then
        M = Mesh%Body(nc)%GBody(1:3) - A
        Mesh%Body(nc)%GBody(1:3) = A + matmul(Stheta,M)
    end if
!end do

! End
end subroutine InitMesh

subroutine Regeneration_Mesh(Mesh,Ecoulement,ti,boolRemesh,boolRemeshFS,fgeom_vect,fdomaine,nface,Grid,nb_point, nb_tri,IndBody_former,size_IndBody_former,Nnodes,Nnodes_FS,ierror,InputData,nRemesh,n_tab2,n_tab,ForcedRemesh,CrossingFS)

    !f2py integer*1, dimension(1000)                        :: Mesh
    type(TMaillage),intent(inout)                           :: Mesh                 ! Total mesh (domain and floater).
    !f2py integer*1, dimension(1000)                        :: Ecoulement
    type(TEcoulement),intent(inout)                         :: Ecoulement           ! Old and new flow parameters.
    real(rp),intent(inout)                                  :: ti                   ! Current time (of the RK4 step).
    logical,intent(inout)                                   :: boolRemesh           ! Boolean to know if the mesh of the bodies or/and the free surface was regenerated (True) or not (False).
    logical,intent(inout)                                   :: boolRemeshFS         ! = true: remeshing the free surface, = false: no remeshing.
    !f2py integer*1, dimension(1000)                        :: fgeom_vect
    type(type_GeomVect),intent(inout)                       :: fgeom_vect           ! Geometry of the floater (not the domain).
    !f2py integer*1, dimension(1000)                        :: fdomaine
    type(type_geom),intent(in)                              :: fdomaine             ! Geometry of the domain.
    integer,intent(in)                                      :: nface                ! Number of faces in both the floater and the domain geometries.
    !f2py integer*1, dimension(1000)                        :: Grid
    type(MGrid),intent(inout)                               :: Grid                 ! Transitional mesh.
    integer                                                 :: nb_point, nb_tri     ! Number of points and triangles in the Grid.
    integer,intent(in)                                      :: size_IndBody_former  ! Index of the floater nodes of the former mesh.
    integer,dimension(2,size_IndBody_former),intent(inout)  :: IndBody_former       ! (1st node, last node)xNbodies: Index of the floater nodes of the former mesh when new mesh generation.
    integer,intent(inout)                                   :: Nnodes               ! Number of nodes in the mesh.
    integer,intent(inout)                                   :: Nnodes_FS            ! Number of nodes in the mesh of the free surface.
    integer,intent(inout)                                   :: ierror               ! Error flag.
    !f2py integer*1, dimension(1000)                        :: InputData
    type(InputDataStruct),intent(inout)                     :: InputData            ! Input data.
    integer,intent(inout)                                   :: nRemesh              ! Number of remeshing (use for the number of Advance_Front_Remesh.dat).
    integer,intent(inout)                                   :: n_tab2,n_tab         ! Number of intersection curves and lines.
    logical,intent(inout)                                   :: ForcedRemesh         ! Force the remeshing if necessary (crossing the free surface).
    logical,intent(inout)                                   :: CrossingFS           ! A body crossed the free surface (True) or not (False).
    
    logical                                                 :: MeshQualityBodies    ! = false: no remeshing of the bodies, = true: remeshing of the bodies.
    logical                                                 :: MeshQualityFS        ! = false: no remeshing of the free surface, = true: remeshing of the free surface.
    integer                                                 :: nc,jj                ! Loop parameters.
    character(len=50)                                       :: filemaill            ! Output mesh file.
    type(chaine_point_pt),dimension(100)                    :: tab2                 ! Table of the pointers toward the intersection points.
    logical                                                 :: AboveFS              ! All parts of the geometry are above the free surface (True) or not (False).
    logical,dimension(Nbodies+1)                            :: Active_old           ! Former fgeom%Active in case of the remeshing does not work.
            
    ! If the deformation of the mesh is too important, this function will generate a new mesh from the geometry.
    
    filemaill = 'mesh_'//filename
    
    if(Mesh_type.eq.2)then ! Mesh strategy of CC.
        
        ! Updating IndBody_former
        jj = 1
        do nc = Int_Body,Mesh%NBody
            IndBody_former(1,jj) = Mesh%Body(nc)%IndBody(1)
            IndBody_former(2,jj) = Mesh%Body(nc)%IndBody(3)
            jj = jj + 1
        end do
        
        ! Checking the quality of the mesh
        call CheckMesh3D(Mesh,ti,MeshQualityBodies,MeshQualityFS,InputData) ! bool = false: no remeshing, bool = true: remeshing
        
        ! Printing the mesh
        call PlotMaill(filemaill,Mesh)
        
        ! Computation of the intersection curves.
        call adjust_waterline2(Mesh,fgeom_vect,CrossingFS,ti,ierror,InputData,tab2,n_tab2,n_tab)
        if(ierror.ne.0)then
            go to 2 ! Remeshing is useless.
	    end if
        
        ! When there is a crossing of the free surface, to know which body is concerned.
        do nc = 1,NBodies
            
            ! Is a body above the free surface?
            if(RemeshFS)then
                call isGeom_aboveFS(fgeom_vect%geom(nc),ti,AboveFS)
            else
                AboveFS = .false.
            end if
            
            Active_old(nc) = fgeom_vect%Active(nc)
            
            if(not(AboveFS))then
                fgeom_vect%Active(nc) = .true.  ! Body immerged or piercing.
            else
                fgeom_vect%Active(nc) = .false. ! Body above the free surface.
            end if
        end do
        
        ! Flag to known if the free surface mesh needs to be remeshed or not.
        DeformFS =  not(lineaireFS) .or. lineaireFS .and. not(lineaireBody) .and.  not(is_immerged)
        
        ! Remeshing (True) or not (False).
        boolRemesh = MeshQualityBodies .or. MeshQualityFS .or. CrossingFS .or. ForcedRemesh ! Either the quality of the mesh is not good enough (MeshQuality = True) or a body goes across the free surface (CrossingFS = True) or the remeshing is forced (ForcedRemesh).
        
        ! boolRemeshFS.
        boolRemeshFS = MeshQualityFS .or. CrossingFS .or. ForcedRemesh ! Remeshing the FS mesh if one of its panels needs it or if a body goes across the free surface or if a forced remeshing is aksed.
        
        ! Forced remeshing of the free surface (all the mesh actually).
        if(ForcedRemeshFS)then
            boolRemesh = .true.
            boolRemeshFS = .true.
        end if
                
        ! Remeshing if necessary.
        if(boolRemesh)then
            
            print*,""
            print*,"Remeshing wanted:"
            print*,"Bodies mesh quality: ",MeshQualityBodies
            print*,"Free surface mesh quality: ",MeshQualityFS
            print*,"Crossing free surface: ",CrossingFS
            print*,"Forced remeshing: ",ForcedRemesh
            do nc = 1,NBodies
                print*,"Body: ",nc+1,"active: ",fgeom_vect%Active(nc)
            end do
            print*,"Number of intersection curves: ",n_tab2
            print*,"is_immerged: ",is_immerged
            print*,"DeformFS: ",DeformFS
            print*,"FS Remehing asked: ",boolRemeshFS
            print*,""
            
            nRemesh = nRemesh + 1
            
            ! Information on the old number of nodes for each body.
            print*,""
            print*,"Free surface:"
            print*,"Nodes: ",Mesh%FS%IndFS(1), Mesh%FS%IndFS(3)
            print*,"Panels: ",Mesh%FS%IndFS(2), Mesh%FS%IndFS(4)
            print*,""
            print*,"Tank:"
            print*,"Nodes: ",Mesh%Body(1)%IndBody(1), Mesh%Body(1)%IndBody(3)
            print*,"Panels: ",Mesh%Body(1)%IndBody(2), Mesh%Body(1)%IndBody(4)
            do nc = Int_Body,Mesh%NBody
                if(Mesh%Body(nc)%Active)then
                    print*,""
                    print*,"Body:",nc
                    print*,"Nodes: ",Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)
                    print*,"Panels: ",Mesh%Body(nc)%IndBody(2), Mesh%Body(nc)%IndBody(4)
                end if
            end do
            print*,""
            
            1 continue
            
            ! Remeshing.
            if(RemeshFS .or. (not(RemeshFS) .and. is_body))then ! No remeshing at all if no body and no remeshing of the FS.
                call compute_mesh_fgeom(Mesh,Ecoulement,ti,fgeom_vect,fdomaine,nface,Grid,nb_point,nb_tri,ierror,InputData,nRemesh,boolRemeshFS,tab2,n_tab2,n_tab)
            end if         
            
            ! Monitoring.
            if(iwmonitor)then
                call PlotMonitor(ti,MeshQualityBodies,MeshQualityFS,CrossingFS,ForcedRemesh,boolRemeshFS,boolRemesh,ierror)
            end if
            
            ! Information on the new number of nodes for each body.
            print*,""
            print*,"Free surface:"
            print*,"Nodes: ",Mesh%FS%IndFS(1), Mesh%FS%IndFS(3)
            print*,"Panels: ",Mesh%FS%IndFS(2), Mesh%FS%IndFS(4)
            print*,""
            print*,"Tank:"
            print*,"Nodes: ",Mesh%Body(1)%IndBody(1), Mesh%Body(1)%IndBody(3)
            print*,"Panels: ",Mesh%Body(1)%IndBody(2), Mesh%Body(1)%IndBody(4)
            do nc = Int_Body,Mesh%NBody
                if(Mesh%Body(nc)%Active)then
                    print*,""
                    print*,"Body:",nc
                    print*,"Nodes: ",Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)
                    print*,"Panels: ",Mesh%Body(nc)%IndBody(2), Mesh%Body(nc)%IndBody(4)
                end if
            end do
            print*,""
            
            ! Updating the number of nodes
            Nnodes = Mesh%Nnoeud
            Nnodes_FS = Mesh%FS%IndFS(3) - Mesh%FS%IndFS(1) + 1
            
            if(ierror/=0) then
                
                print*,"Regeneration_Mesh: error during the remeshing process, former mesh reused."
                print*,""
                
                ! Keep the previous value of fgeom_vect%Active.
                do nc = 1,NBodies
                    fgeom_vect%Active(nc) = Active_old(nc)
                end do
                
                ! Updating the geometrical properties
                call GeomInit(Mesh, fgeom_vect,ti, InputData,.true.)
                
            else ! No error
                
                if(ForcedRemesh)then
                    ! Forced remshing worked, the condition on ForcedRemesh can be deleted.
                    ForcedRemesh = .false.
                end if
                
            end if
            
            ! If there is a crossing of the free surface, a remeshing is done at the next time step as well.
            if(CrossingFS)then
                ForcedRemesh = .true.
            end if
        
        else ! No remeshing.
            
            ! Monitoring.
            if(iwmonitor)then
                call PlotMonitor(ti,.false.,.false.,.false.,.false.,.false.,.false.,0)
            end if
            
        end if
        
    else ! Mesh stategy of LL.
        
        ! Updating the geometrical properties
        call GeomInit(Mesh, fgeom_vect,ti, InputData,.true.)
        
        ! boolRemesh.
        boolRemesh = .false.
        
        ! boolRemeshFS.
        boolRemeshFS = .false.
        
        ! Monitoring.
        if(iwmonitor)then
            call PlotMonitor(ti,.false.,.false.,.false.,.false.,.false.,.false.,0)
        end if
        
    end if 
    
    2 continue
    
end subroutine Regeneration_Mesh

subroutine Zeroing_new_points(Mesh,Mesh0,Ecoulement,Ecoulement0,bgenerate,boolRemeshFS,IndBody_former,size_IndBody_former,ierror)
    
    !f2py integer*1, dimension(1000)                    :: Mesh
    type(TMaillage),intent(inout)                       :: Mesh                     ! Total mesh (domain and floater).
    !f2py integer*1, dimension(1000)                    :: Mesh0
    type(TMaillage),intent(inout)                       :: Mesh0                    ! Local mesh (domain and floater).
    !f2py integer*1, dimension(1000)                    :: Ecoulement,Ecoulement0
    type(TEcoulement)                                   :: Ecoulement,Ecoulement0   ! Flow parameters (phi, etc.).
    logical,intent(in)                                  :: bgenerate                ! Boolean to know if the mesh needs to be generated (true) or not (false).
    logical,intent(in)                                  :: boolRemeshFS             ! = true: remeshing the free surface, = false: no remeshing.
    integer,intent(in)                                  :: size_IndBody_former      ! Index of the floater nodes of the former mesh.
    integer,dimension(2,size_IndBody_former),intent(in) :: IndBody_former           ! (1st node, last node)xNbodies: Index of the floater nodes of the former mesh when new mesh generation.
    integer,intent(inout)                               :: ierror                   ! Error flag.

    integer                                             :: j,nc,jj                  ! Loop parameters.
    
    ! This subroutine zeros the physical quantities for the new points created by the subroutine Regeneration_mesh and creates a new Mesh0.

    if (Mesh_type.eq.2) then
        if (bgenerate) then
            if(ierror == 0) then
                
                ! New Mesh0.
                call DelMaillage(Mesh0)
                call NewMaillage(Mesh0,Mesh%Nnoeud,NBodies+1) ! +1 for the tank.
                
                ! New Ecoulement0.
                call DelEcoulement(Ecoulement0)
                call NewEcoulement(Ecoulement0, Mesh%Nnoeud)
                call IniEcoulement(Ecoulement0, Mesh%Nnoeud, 0._RP)
                
                ! Mise à zero du potentiel des nouveaux points sur le corps.
                jj = 1
                do nc = Int_Body,Mesh%NBody
                    if(Mesh%Body(nc)%Active)then
                        do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                            Ecoulement%Phi(j)%Perturbation = 0._RP
                            Ecoulement%DPhiDt(j)%Perturbation = 0._RP
                            Ecoulement%DDPhiDnDt(j)%Perturbation = 0._RP
                        end do
                    end if
                    jj = jj + 1
                end do
                
            end if
        end if    
    end if

end subroutine Zeroing_new_points

end module Remesh_mod