module GeomMesh
use FonctionsCommunes
use Incident_mod
use GeomStruct
use GeomAxiSym
use GeomCylindre
use Spline
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                       Calcul des caractéristiques géométriques d'un maillage                     !
!                                 LETOURNEL Lucas    06/12                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GeomInit(Mesh, fgeom_vect, t, InputData,bool, ierror)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(inout)      :: Mesh                                                 ! Mesh.
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                                           ! Geometries
    real(rp), intent(in)                :: t                                                    ! Current time.
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData                                            ! Input data.
    logical,intent(in),optional         :: bool                                                 ! Boolean to define same_topo.
    integer, optional                   :: ierror                                               ! Error flag.
    
    logical                             :: error                                                ! Error flag.
    integer                             :: j, k, jk, r, n1, i1, jtemp, nc,jj                    ! Loop parameters.
    integer                             :: ind_body0, ind_body1                                 ! 1st and last nodes for a body.
    real(rp)                            :: delta, Alpha, InvPi                                  !
    real(rp)                            :: M0, vol                                              ! Mass and volume.
    real(rp), dimension(3)              :: Nor, ds1, ds2, M, GEta0,VectResult_1,VectResult_2    !
    real(rp), dimension(3,3)            :: PG                                                   !
    real(rp), dimension(3,6)            :: Pt                                                   !
    real(rp)                            :: Crampe_1,Crampe_2,Crampe_3,Crampe_4,Crampe_5,Crampe_6! Ramps
    logical                             :: same_topo, isPC                                      !
    
    ! This subroutine computes the geometrical properties of the mesh.
    
    error = .false.
    
    ! same_topo = false when creating the mesh or remeshing;
    !           = true if no remeshing.
    
    if(present(bool))then
        same_topo = bool
    else
        same_topo = .false.
    endif
    
    if(.not.same_topo)then
        Mesh%Tnoeud(:)%Nfacette = 0
        Mesh%Tnoeud(:)%NDouble = 0
        Mesh%Tnoeud(:)%TypeDouble = 1
    endif
    InvPi = 1._RP/Pi
    
    ! is_piercing and %Active for the floaters
    if(not(same_topo))then ! %Active can only change when there is a remeshing.
        call Updating_is_piercing(Mesh,fgeom_vect,Inputdata) ! Use %NPanneau so must be called after the definition of %NPanneau (here-above).
    end if
    
    ! Tank
    Mesh%Body(1)%Active = .true.
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                   Free surface
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SL Facettes
    ! Pour chaque sommet de la facette, on renvoie l'indice de la facette et la position du noeud (1, 2 ou 3) dans le tableau de dénombrement des facettes dont le noeud est un sommet 
    if(.not.same_topo)then
        if(idebug > 0) print*," GeomInit : recherche facettes appartenant au noeud"
        do j=Mesh%FS%IndFS(2),Mesh%FS%IndFS(4)
            do r=1,3
                Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette = Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette + 1
                Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Tfacette(Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette,:)=[j,r]
            end do    
        enddo
    endif
    
    if(idebug > 0) print*," GeomInit : calcul G, aire, normale, Rmax"
    
    Mesh%FS%Aire = 0._RP
    do j=Mesh%FS%IndFS(2),Mesh%FS%IndFS(4)

        Mesh%Tfacette(j)%Pnoeud=reshape([Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(2))%Pnoeud,&
                                        & Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(3))%Pnoeud],(/3,3/))
                                    
        ! Calcul du centre de gravité de chaque facette, de son aire, de sa normale et de la distance maximale du centre aux sommets    
        Pg=Mesh%Tfacette(j)%Pnoeud
        Mesh%Tfacette(j)%Gfacette=sum(Pg,2)/3._RP
        call Computation_vect_product(Pg(:,2)-Pg(:,1), Pg(:,3)-Pg(:,1),Nor)
        if (norm2(Nor)<Epsilon) then
            print *, '        probleme de normale pour la facette', j, ' de coordonnées', Pg
            error = .true.
            go to 9999
        else
            Mesh%Tfacette(j)%Normale=Nor/norm2(Nor)
        end if
        call Computation_vect_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1),VectResult_1)
        Mesh.Tfacette(j).Aire=0.5_RP*norm2(VectResult_1)    
        if (Mesh.Tfacette(j).Aire.lt.Epsilon) then
            print *, '        probleme Aire : facette plate', j, ' de coordonnées', Pg
            error = .true.
            go to 9999
        end if
        
        ! Maximal distance between the gravity center and the apexes.
        Mesh%Tfacette(j)%Rmax = max(norm2(Pg(:,1)-Mesh%Tfacette(j)%Gfacette),&
                                      &norm2(Pg(:,2)-Mesh%Tfacette(j)%Gfacette),&
                                      &norm2(Pg(:,3)-Mesh%Tfacette(j)%Gfacette))  
                                  
       ! Calcul du gradient constant surfacique de la facette (ne dépend que de la géométrie de la facette)    
        delta = (norm2(Pg(:,2)-Pg(:,1))**2) * norm2(Pg(:,3)-Pg(:,1))**2 - dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))**2 ! ||P1P2||^2*||P1P3||^2 - (P1P2.P1P3)^2
        ds1 = (norm2(Pg(:,3)-Pg(:,1))**2*(Pg(:,2)-Pg(:,1))-(dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1))))*(Pg(:,3)-Pg(:,1))) ! ||P1P3||^2.P1P2 - (P1P2.P1P3).P1P3
        ds2 = (-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))*(Pg(:,2)-Pg(:,1))+norm2(Pg(:,2)-Pg(:,1))**2*(Pg(:,3)-Pg(:,1))) ! -(P1P2.P1P3).P1P2 + ||P1P2||^2.P1P3
        Mesh%Tfacette(j)%ds = reshape([-(ds1+ds2),ds1,ds2]/delta,(/3,3/))
        Mesh%Tfacette(j)%dsT = transpose(Mesh%Tfacette(j)%ds)
        
    end do
    
    ! SL Noeuds
    if(idebug>0) print*," GeomInit : calcul angle sommet, normale, zone abs"
    
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        
        ! Calcul de l'angle au sommet
        Mesh%Tnoeud(j)%Angle = 0._RP
        do k=1,Mesh%Tnoeud(j)%Nfacette
            n1=Mesh%Tnoeud(j)%Tfacette(k,1)
            i1=Mesh%Tnoeud(j)%Tfacette(k,2)
            Pt(:,1:3)=Mesh%Tfacette(n1)%Pnoeud
            Pt(:,4:6)=Mesh%Tfacette(n1)%Pnoeud
            Alpha=acos(dot_product(Pt(:,i1+1)-Pt(:,i1),Pt(:,i1+2)-Pt(:,i1))/(norm2(Pt(:,i1+1)-Pt(:,i1))*norm2(Pt(:,i1+2)-Pt(:,i1))))
            Mesh%Tnoeud(j)%Angle(k) = Alpha
        end do
        Mesh%Tnoeud(j)%Angle = Mesh%Tnoeud(j)%Angle*InvPi

	    Mesh%Tnoeud(j)%Aire = 0._rp
        do k=1,Mesh%Tnoeud(j)%Nfacette
		    Mesh%Tnoeud(j)%Aire = Mesh%Tnoeud(j)%Aire + Mesh%TFacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire * Mesh%Tnoeud(j)%Angle(k)
	    enddo

        ! Calcul de la Normale aux Noeuds
        if (Htype.eq.0) then
            Mesh%Tnoeud(j)%Normale = [0._RP,0._RP,-1._RP]
            Mesh%Tnoeud(j)%DLocal(1:2) = [0._RP,0._RP]
            Mesh%Tnoeud(j)%DLocal(5) =  1._RP
        else
            ! Calcul de la normale au point sur la SL d'après la variation locale (gradient surfacique) de la déformée de surface libre incidente
            M = Mesh%Tnoeud(j)%Pnoeud
            call CGEta0(M, t, GEta0)  
            Mesh%Tnoeud(j)%Normale = [GEta0(1),GEta0(2),-1._RP]/norm2([GEta0(1),GEta0(2),-1._RP]) !-vect_product(A,B)/norm2(vect_product(A,B))
            Mesh%Tnoeud(j)%DLocal(1:2) = GEta0(1:2)
            Mesh%Tnoeud(j)%DLocal(5) = GEta0(1)**2 + GEta0(2)**2 + 1._RP            
        end if
        
        ! Vérification du type du noeud par rapport au type des facettes liées au noeud
        if (idebug.gt.0) then
            do k=1,Mesh%Tnoeud(j)%Nfacette
                if ((Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%TypeFrontiere.ne.0)) then
                    print*, 'problème de type de Frontiere pour le noeud', j
                    error = .true.
                    go to 9999
                end if
            end do
        end if
        
        ! Définition de la zone de Damping
        Mesh%Tnoeud(j)%Damping = 0._RP
        M = Mesh%Tnoeud(j)%Pnoeud
        M(3) = 0
        M0 = Ldom(5)-LAbs
        
        if(abs(LAbs).le. 2*Epsilon2)then
            Mesh%Tnoeud(j)%Damping = 0._rp
        elseif (Mesh_type.eq.1) then
            if (idtype.eq.1) then
                print*,"Want to use Lgeom but for which body?"
                print*,"Next lines are commented."
                pause
                !if (igtype.eq.6) then !Wigley Hull damping condition
                    !call Computation_Crampe(M(1), 0.5_RP*Lgeom(1)+Ldom(4)-Labs,0.5_RP*Lgeom(1)+Ldom(4),Crampe_1)
                    !call Computation_Crampe(-M(1),0.5_RP*Lgeom(1)+Ldom(4)+Ldom(1)-Labs,0.5_RP*Lgeom(1)+Ldom(4)+Ldom(1),Crampe_2)
                    !Mesh%Tnoeud(j)%Damping = Crampe_1 + Crampe_2
                    !call Computation_Crampe(M(2),Ldom(2)-Labs,Ldom(2),Crampe_3)
                    !Mesh%Tnoeud(j)%Damping = Mesh%Tnoeud(j)%Damping*(1-Crampe_3) + Crampe_3
                !else
                    call Computation_Crampe(M(1), Mesh%DimTank(1)/2._RP-Labs,Mesh%DimTank(1)/2._RP,Crampe_1)
                    call Computation_Crampe(-M(1),Mesh%DimTank(1)/2._RP-Labs,Mesh%DimTank(1)/2._RP,Crampe_2)
                    Mesh%Tnoeud(j)%Damping = Crampe_1 + Crampe_2
                    call Computation_Crampe(M(2), Mesh%DimTank(2)-Labs,Mesh%DimTank(2),Crampe_3)
                    Mesh%Tnoeud(j)%Damping = Mesh%Tnoeud(j)%Damping*(1-Crampe_3) + Crampe_3
                !end if
            elseif(idtype.eq.2) then
                call Computation_Crampe(norm2(M),Mesh%DimTank(1)-Labs,Mesh%DimTank(1),Crampe_4)
                Mesh.Tnoeud(j)%Damping = Crampe_4
            end if
            if (.not.cuve_ferme) Mesh%Tnoeud(j)%Damping = 1._RP - Mesh%Tnoeud(j)%Damping
        elseif (Mesh_type.eq.2 .and. (idtype==1 .or. idtype.eq.3)) then
            if(M(1).ge.M0)then
                call Computation_Crampe(M(1),M0,0.5*Ldom(1),Crampe_5)
                Mesh%Tnoeud(j)%Damping = Crampe_5
                Mesh%Tnoeud(j)%mobility = .false.
            else
                Mesh%Tnoeud(j)%Damping = 0._rp
                Mesh%Tnoeud(j)%mobility = .true.
            endif
        elseif(Mesh_type.eq.2 .and. (idtype==2 .or. idtype.eq.4))then
            if(norm2(M).ge.M0)then
                call Computation_Crampe(norm2(M),M0,Ldom(5),Crampe_6)
                Mesh%Tnoeud(j)%Damping = Crampe_6
                Mesh%Tnoeud(j)%mobility = .false.
            else
                Mesh%Tnoeud(j)%Damping = 0._rp
                Mesh%Tnoeud(j)%mobility = .true.
            endif
        else
            print*, ' Error : Mesh_type = ', Mesh_type
            error = .true.
            go to 9999
        end if
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Bodies (tank included)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Pour chaque sommet de la facette, on renvoie l'indice de la facette et la position du noeud (1, 2 ou 3) dans le tableau de dénombrement des facettes dont le noeud est un sommet 
    if(.not.same_topo)then
        if(idebug>0) print*," GeomInit : recherche facette appartenant au noeud (flotteur)"
        do nc=1,Mesh%NBody
            if(Mesh%Body(nc)%Active)then
                Mesh%Body(nc)%Aire = 0._RP
                Mesh%Body(nc)%Volume = 0._RP
                ! Facettes
                do j=Mesh%Body(nc)%IndBody(2),Mesh%Body(nc)%IndBody(4)        
                    do r=1,3
                        Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette = Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette + 1
                        Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Tfacette(Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(r))%Nfacette,:)=[j,r]
                    end do   
                enddo
            end if
        enddo
    endif
    if(idebug > 0) print*," GeomInit : calcul G, aire, normale, Rmax"
    
    ! Panels   
    do nc=1,Mesh%NBody
        
        if(Mesh%Body(nc)%Active)then            
            Mesh%Body(nc)%BBody(1:3) = 0._RP
            Mesh%Body(nc)%Aire = 0._RP
            Mesh%Body(nc)%Volume = 0._RP
                        
            do j=Mesh%Body(nc)%IndBody(2),Mesh%Body(nc)%IndBody(4)
	            Mesh%Tfacette(j)%Pnoeud=reshape([Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(2))%Pnoeud,&
                                                & Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(3))%Pnoeud],(/3,3/)) 
                ! Calcul du centre de gravité de chaque facette, de son aire, de sa normale et de la distance maximale du centre aux sommets    
                Pg=Mesh%Tfacette(j)%Pnoeud ! Coordinates of the three points of the panel.
                Mesh%Tfacette(j)%Gfacette=sum(Pg,2)/3._RP ! sum(Pg,2) = Sum of the components (sum of xj, sum of yj, sum of zj) Eq 3.2 of the PhD of LL.
                call Computation_vect_product(Pg(:,2)-Pg(:,1), Pg(:,3)-Pg(:,1),Nor)
            
                ! Calcul de la normale à la facette
                if (norm2(Nor)/(norm2(Pg(:,2)-Pg(:,1))*norm2(Pg(:,3)-Pg(:,1))) .lt. Epsilon) then
                    print *, '        probleme de normale pour la facette', j, ' de coordonnées', Pg
                    error = .true.
                    go to 9999
                else
                    ! Test de la direction de la normale sur le corps
                    if ((nc.ge.Int_Body).and.dot_product(Nor, Mesh%Tfacette(j)%Gfacette-Mesh%Body(nc)%GBody(1:3)).lt.Epsilon) then ! Bodies except the tank
                        Nor = - Nor
                        jtemp = Mesh%Tfacette(j)%Tnoeud(2)
                        Mesh%Tfacette(j)%Tnoeud(2) = Mesh%Tfacette(j)%Tnoeud(3)
                        Mesh%Tfacette(j)%Tnoeud(3) = jtemp
                        Mesh%Tfacette(j)%Pnoeud=reshape([Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(1))%Pnoeud, Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(2))%Pnoeud,&
                                                & Mesh%Tnoeud(Mesh%Tfacette(j)%Tnoeud(3))%Pnoeud],(/3,3/))    
                        Pg=Mesh%Tfacette(j)%Pnoeud
                    end if          
                    Mesh%Tfacette(j)%Normale=Nor/norm2(Nor)
                end if
            
                ! Aire
                call Computation_vect_product(Pg(:,2)-Pg(:,1),Pg(:,3)-Pg(:,1),VectResult_2)
                Mesh%Tfacette(j)%Aire=0.5_RP*norm2(VectResult_2)    
                if (Mesh%Tfacette(j)%Aire/(norm2(Pg(:,2)-Pg(:,1))*norm2(Pg(:,3)-Pg(:,1))) .lt. Epsilon) then
                    print *, '        probleme Aire : facette plate', j, ' de coordonnées', Pg
                    error = .true.
                    go to 9999
                end if
                Mesh%Tfacette(j)%Rmax=max(norm2(Pg(:,1)-Mesh%Tfacette(j)%Gfacette),&
                                              &norm2(Pg(:,2)-Mesh%Tfacette(j)%Gfacette),&
                                              &norm2(Pg(:,3)-Mesh%Tfacette(j)%Gfacette))    
               ! Calcul du gradient constant surfacique de la facette (ne dépend que de la géométrie de la facette)    
                delta=(norm2(Pg(:,2)-Pg(:,1))**2)*norm2(Pg(:,3)-Pg(:,1))**2-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))**2
                ds1=(norm2(Pg(:,3)-Pg(:,1))**2*(Pg(:,2)-Pg(:,1))-(dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1))))*(Pg(:,3)-Pg(:,1)))
                ds2=(-dot_product((Pg(:,2)-Pg(:,1)),(Pg(:,3)-Pg(:,1)))*(Pg(:,2)-Pg(:,1))+norm2(Pg(:,2)-Pg(:,1))**2*(Pg(:,3)-Pg(:,1)))
                Mesh%Tfacette(j)%ds = reshape([-(ds1+ds2),ds1,ds2]/delta,(/3,3/))
                Mesh%Tfacette(j)%dsT=transpose(Mesh%Tfacette(j)%ds)
            end do
            
            ! Noeuds
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                
                ! Calcul de l'angle au sommet
                Mesh%Tnoeud(j)%Angle = 0._RP
                do k=1,Mesh%Tnoeud(j)%Nfacette
                    n1=Mesh%Tnoeud(j)%Tfacette(k,1)
                    i1=Mesh%Tnoeud(j)%Tfacette(k,2)
                    Pt(:,1:3)=Mesh%Tfacette(n1)%Pnoeud
                    Pt(:,4:6)=Mesh%Tfacette(n1)%Pnoeud
                    Alpha=acos(dot_product(Pt(:,i1+1)-Pt(:,i1),Pt(:,i1+2)-Pt(:,i1))/(norm2(Pt(:,i1+1)-Pt(:,i1))*norm2(Pt(:,i1+2)-Pt(:,i1))))
                    Mesh%Tnoeud(j)%Angle(k) = Alpha
                end do
                Mesh%Tnoeud(j)%Angle = Mesh%Tnoeud(j)%Angle*InvPi
            
                ! Calcul de la Normale et de l'Aire aux Noeuds
                Mesh%Tnoeud(j)%Normale = 0._RP
                Mesh%Tnoeud(j)%Aire = 0._RP
                do k=1,Mesh%Tnoeud(j)%Nfacette
                    Mesh%Tnoeud(j)%Normale = Mesh%Tnoeud(j)%Normale + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Normale*Mesh%Tnoeud(j)%Angle(k)
                    Mesh%Tnoeud(j)%Aire = Mesh%Tnoeud(j)%Aire + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*Mesh%Tnoeud(j)%Angle(k)
                end do
            
                if (Symmetry.and.abs(Mesh%Tnoeud(j)%Pnoeud(2)).lt.Epsilon) then
                    Mesh%Tnoeud(j)%Normale = 2._RP*Mesh%Tnoeud(j)%Normale
                    Mesh%Tnoeud(j)%Normale(2) = 0._RP
                end if
            
                if (norm2(Mesh%Tnoeud(j)%Normale).lt.Epsilon) then
                    print*, 'GeomInit: the length of the normale vector on the mesh of a floater is null. Floater: ',nc,', node: ',j,", Nfacette: ",Mesh%Tnoeud(j)%Nfacette,", Npanneau: ", Mesh%Tnoeud(j)%Npanneau,", area: ",Mesh%Tnoeud(j)%Aire 
                    print*,"Pnoeud:",Mesh%Tnoeud(j)%Pnoeud
                    print*,"Normal:",Mesh%Tnoeud(j)%Normale
                    Mesh%Tnoeud(j)%Normale = 0
                    error = .true.
                    go to 9999
                    do k=1,Mesh%Tnoeud(j)%Nfacette
                        print*, Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Normale, Mesh%Tnoeud(j)%Angle(k)
                    end do
                else
                    Mesh%Tnoeud(j)%Normale = Mesh%Tnoeud(j)%Normale/norm2(Mesh%Tnoeud(j)%Normale)
                end if
            
                Mesh%Body(nc)%BBody(1:3) = Mesh%Body(nc)%BBody(1:3) + Mesh%Tnoeud(j)%Pnoeud(1:3)*Mesh%Tnoeud(j)%Aire
                Mesh%Body(nc)%Aire = Mesh%Body(nc)%Aire + Mesh%Tnoeud(j)%Aire
                Mesh%Body(nc)%Volume = Mesh%Body(nc)%Volume + dot_product(Mesh%Tnoeud(j)%Pnoeud,Mesh%Tnoeud(j)%Normale)*Mesh%Tnoeud(j)%Aire
        
                ! Détermination des Noeuds Doubles
                if (j.gt.1 .and. .not.same_topo) then
                    do k=1,j-1
                        if (norm2(Mesh%Tnoeud(j)%Pnoeud-Mesh%Tnoeud(k)%Pnoeud).lt.Epsilon) then
                            Mesh%Tnoeud(j)%Ndouble=Mesh%Tnoeud(j)%Ndouble+1
                            Mesh%Tnoeud(k)%Ndouble=Mesh%Tnoeud(k)%Ndouble+1
                            Mesh%Tnoeud(j)%double(Mesh%Tnoeud(j)%Ndouble)=k
                            Mesh%Tnoeud(k)%double(Mesh%Tnoeud(k)%Ndouble)=j
                        
                            ! Initialization Typedouble
                            if(abs(Mesh%Tnoeud(j)%Normale(1)).lt.Epsilon .and. abs(Mesh%Tnoeud(j)%Normale(2)).lt.Epsilon .and. abs(Mesh%Tnoeud(j)%Normale(3)+1._RP).lt.Epsilon)then ! Necessity for MeshVelBody (Remesh)
                                Mesh%Tnoeud(j)%Typedouble = 0 ! The normal is equal to (0,0,-1) -> the node is linked to the free surfarce or the low disc of a floating cylinder (Normal_SL).
                            end if
                            if(abs(Mesh%Tnoeud(k)%Normale(1)).lt.Epsilon .and. abs(Mesh%Tnoeud(k)%Normale(2)).lt.Epsilon .and. abs(Mesh%Tnoeud(k)%Normale(3)+1._RP).lt.Epsilon)then ! Necessity for MeshVelBody (Remesh)
                                Mesh%Tnoeud(k)%Typedouble = 0 ! The normal is not equal to (0,0,-1) -> the node is linked to the body of the floater (Normal_SM).
                            end if
                        
                            ! Checking the logical choise of Typedouble
                            if(Mesh%Tnoeud(j)%Typedouble.eq.0 .and. Mesh%Tnoeud(k)%Typedouble.eq.0)then
                                print*,"Normals of the twin nodes ",j," and ",k," are equal to (0,0,-1)"
                            end if
                        
                        end if
                    end do
                    if (Mesh%Tnoeud(j)%Ndouble .gt.2) then
                        print*, 'Probleme Noeuds Doubles', j, Mesh%Tnoeud(j)%Ndouble, Mesh%Tnoeud(j)%double
                        print*, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(1))%Pnoeud, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(1))%TypeNoeud
                        print*, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(2))%Pnoeud, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(2))%TypeNoeud
                        print*, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(3))%Pnoeud, Mesh%Tnoeud(Mesh%Tnoeud(j)%double(3))%TypeNoeud
                    end if
                end if        
            
                ! Vérification du type du noeud par rapport au type des facettes liées au noeud
                if (idebug.gt.0) then
                    do k=1,Mesh%Tnoeud(j)%Nfacette
                        if (Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%TypeFrontiere .ne. Mesh%Tnoeud(j)%TypeNoeud) then
                            print*, 'problème de type de Frontiere pour le noeud', j
                            error = .true.
                            go to 9999
                        end if
                    end do
                end if
                Mesh%Tnoeud(j)%Damping = 0._RP
            end do
            
            if(Mesh%Body(nc)%Aire.ne.0)then
                Mesh%Body(nc)%BBody(1:3) = Mesh%Body(nc)%BBody(1:3)/Mesh%Body(nc)%Aire
                if (Symmetry) Mesh%Body(nc)%BBody(2) = 0._RP
            end if
            
            Mesh%Body(nc)%Volume = Mesh%Body(nc)%Volume/3._RP
            if(Htype.eq.0 .and. abs(Mesh%Body(nc)%VBody(3).lt.Epsilon))then
                ind_body0 = Mesh%Body(nc)%IndBody(1)
                ind_body1 = Mesh%Body(nc)%IndBody(3)
                Mesh%Tnoeud(ind_body0:ind_body1)%mobility = .false.
            endif
            
        end if
    end do
        
    ! Spline order
    ! Bodies
    if (is_Body) then
        do nc=Int_body,Mesh%NBody
            if(Mesh%Body(nc)%Active)then
                Mesh%Tnoeud(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Ordre = SplineOrder
            end if
        end do
    end if
        
    ! Tank
    if (cuve_ferme) then
        do nc=1,Mesh%NBody
            if (nc.lt.Int_Body) then ! Tank except bodies
                Mesh%Tnoeud(Mesh%Body(nc)%IndBody(1):Mesh%Body(nc)%IndBody(3))%Ordre = 0
            end if
        end do
    end if
        
    ! Free surface
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        if(Mesh%Tnoeud(j)%Ndouble.gt.0 .and. norm2(Mesh%Tnoeud(j)%Pnoeud).lt.Ldom(5)-LAbs)then
            Mesh%Tnoeud(j)%Ordre = 3
            !Mesh%Tnoeud(j)%Ordre = SplineOrder
        else
            if(RemeshFS)then ! As soon as the FS remeshing can be activated, the spline order is fixed during the simulation.
                Mesh%Tnoeud(j)%Ordre = SplineOrder
            else
                Mesh%Tnoeud(j)%Ordre = 0
            end if
            !Mesh%Tnoeud(j)%Ordre = SplineOrder
            !Mesh%Tnoeud(j)%Ordre = 0
        endif
    enddo
        
    Mesh%Tnoeud(1:Mesh%Nnoeud)%RVoisin = RSpline ! Useless
        
    if(.not.same_topo)then
        if(idebug>0) print*," GeomInit : recherche des noeuds voisins"
        call CVoisin(Mesh,fgeom_vect)
    endif
        
    ! Recherche des points de controle
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        if(.not.Mesh%Tnoeud(j)%mobility)then ! Controle points on the line (r = R - Labs)
            isPC = .false.
            k = 1
            do while (k .le.Mesh%Tnoeud(j)%Nvoisin(2) .and. Mesh%Tnoeud(j)%Tvoisin(k,2).lt.2)
                jk = abs(Mesh%Tnoeud(j)%Tvoisin(k,1))
                isPC = isPC .or. Mesh%Tnoeud(jk)%mobility==.true.
                k = k +1
            enddo
            Mesh%Tnoeud(j)%control_point = isPC
        !if(abs(sqrt(Mesh%Tnoeud(j)%Pnoeud(1)**2 + Mesh%Tnoeud(j)%Pnoeud(2)**2)-Ldom(5)).lt.Epsilon .and. Mesh%Tnoeud(j)%Ndouble.ge.1)then ! Tank
        !    Mesh%Tnoeud(j)%control_point = .true.
        elseif(Mesh%Tnoeud(j)%mobility .and. Mesh%Tnoeud(j)%Ndouble.ge.1)then ! Intersection Floater - Free surface
            Mesh%Tnoeud(j)%control_point = .true.
        else
            Mesh%Tnoeud(j)%control_point = .false.
        endif
    enddo
        
    ! Calcul de la masse (à passer dans la structure du maillage par la suite)
    jj = 1
    do nc=Int_Body,Mesh%NBody
        if(Mesh%Body(nc)%Active)then
            if(is_body .and. InputData%Mass(jj).lt.Epsilon) then ! Negative mass in input data.
                if(idebug>0) print*," GeomInit : calcul de la masse du flotteur"
                vol = 0._RP
                Mesh%Body(nc)%Mass = 0._RP
                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    vol = vol + Mesh%Tnoeud(j)%Pnoeud(3)*Mesh%Tnoeud(j)%Normale(3)*Mesh%Tnoeud(j)%Aire
                end do
                InputData%Mass(jj) = ro*vol
                Mesh%Body(nc)%Mass = ro*vol
            
                if(Symmetry)then
                    InputData%Mass(jj) = 2._rp*InputData%Mass(jj)
                    Mesh%Body(nc)%Mass = 2._rp*Mesh%Body(nc)%Mass
                    vol = 2._rp*vol
                end if
                print*, 'Mass Body ',nc,' : ', Mesh%Body(nc)%Mass
            elseif(is_body .and. InputData%Mass(jj).ge.Epsilon) then ! Strictly positive mass in input data.
                Mesh%Body(nc)%Mass = InputData%Mass(jj)
            endif
        end if
        jj = jj + 1
    end do
        
    10 format('mass_body = ',f16.4,' ; volume = ',f16.4,' ; Tarea = ',f16.4,' ; Narea = ',f16.4)
    20     format('erreur position noeud j = ',i3,' coord = [',3f8.4,'] distance to COG : ',f16.8)
    
    ! NPanneau - Free surface
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        Mesh%Tnoeud(j)%NPanneau = 0
    end do
        
    ! NPanneau - Bodies and CMD
    Mesh%FS%CMD = Htype.ne.0 .and. not(lineaireFS)
    jj = 1
    do nc=1,Mesh%NBody
        if(nc.ge.Int_Body)then ! Floaters
            if(Mesh%Body(nc)%Active)then
                if(Htype.eq.0)then ! Still water
                    Mesh%Body(nc)%CMD = [InputData%free_body(jj) .or. abs(InputData%ACorps(jj)).gt.Epsilon .or. abs(InputData%Vcst(jj)).gt.Epsilon,.true.]
                else ! If the free surface is not still, the remeshing of the floaters is done.
                    Mesh%Body(nc)%CMD = [.true.,.true.]
                end if
                Mesh%Body(nc)%is_tank = .false.
                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    Mesh%Tnoeud(j)%NPanneau = nc
                end do
            end if
            jj = jj + 1
        else ! Tank
            Mesh%Body(nc)%CMD(1) = .false.
            Mesh%Body(nc)%CMD(2) = Htype.ne.0 .and. not(lineaireFS)
            Mesh%Body(nc)%is_tank = .true.
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                Mesh%Tnoeud(j)%NPanneau = 1
            end do
        end if
    end do
        
    ! Local basis (use Mesh%Active so need to be call after Updating_is_piercing)
    if (is_Body) then
        if(idebug>0) print*," GeomInit : base local"
        call BLocal(Mesh)
    end if
        
    ! Flag to known if the floater mesh needs to be remeshed or not
    do nc = 1,NBodies
        if(Mesh%Body(nc)%Active)then
            InputData%DeformBody(nc)  =   not(lineaireBody)   .or.    not(lineaireFS).and. InputData%is_piercing(nc) == 1 ! is_piercing == 1: piercing body
            
        end if
    end do
        
    ! Flag to known if the toal mesh needs to be deform or not
    DeformMesh = .false.
    do nc=1,NBodies
        if(Mesh%Body(nc)%Active)then
            DeformMesh = DeformMesh .or. InputData%DeformBody(nc)
            
        end if
    end do
    DeformMesh = DeformFS .or. DeformMesh
        
    9999 if (error) then 
             print*, 'GeomInit stopped', error
             if (present(ierror)) ierror = error
         end if
         
end subroutine GeomInit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                       Détermination des points voisins pour tout le maillage                     !
!                                 LETOURNEL Lucas    06/13                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CVoisin(Mesh,fgeom_vect)
    !!!!! Problème :
    ! Détermination des points voisins pour tout le maillage pour l'approximation par BSplines.
    !   On relève les points voisins pour chaque points du maillage selon l'ordre de voisinage.
    !   Pour un point donné, le voisinage d'ordre 0 correspond à lui-même. Le voisinage d'ordre 1
    !   correspond à tous les points ayant une facette commune avec lui. Le voisinage d'ordre 2
    !   est étendu logiquement à tous les points ayant une facette commune avec le voisinage d'ordre 1.
    !!!!!
    
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage),intent(inout)       :: Mesh                                                                 ! Mesh
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect                                                           ! Geometries
    
    logical                             :: condition, cond
    logical, dimension(Mesh%Nnoeud)     :: condOrdre
    integer                             :: j, k, r, indice, jface, ordre, jtemp, modificateur,potar_voisin,nc
    integer, dimension(Ncnxmx*NPVMax)   :: jvoisin
    real(rp), dimension(3)              :: Pvo
    
    ! This subroutine gives the neighbours of each node at different order.
    
    ordre = 1
    do indice=1,Mesh%Nsys
        
        ! The first point of TVoisin of Tnoeud(indice) is the point Tnoeud(indice) itself.
        Mesh%Tnoeud(indice)%NVoisin(2) = 1
        Mesh%Tnoeud(indice)%TVoisin(1,1) = indice
        Mesh%Tnoeud(indice)%TVoisin(1,2) = ordre
        jvoisin = 0
        ! Nodes of the panel which are in common with Tnoeud(indice).
        do k=1,Mesh%Tnoeud(indice)%Nfacette
            jface = Mesh%Tnoeud(indice)%Tfacette(k,1)
            jvoisin(3*(k-1)+1) = Mesh%Tfacette(jface)%Tnoeud(1)
            jvoisin(3*(k-1)+2) = Mesh%Tfacette(jface)%Tnoeud(2)
            jvoisin(3*(k-1)+3) = Mesh%Tfacette(jface)%Tnoeud(3)
        end do
        
        ! On teste les indices des points récupérés pour vérifier qu'on ne les a pas déjà.
        do k=1,3*Mesh%Tnoeud(indice)%Nfacette
            condition = .true.
            do r=1,Mesh%Tnoeud(indice)%NVoisin(2)
                if (jvoisin(k).eq.0 .or. jvoisin(k).eq.Mesh%Tnoeud(indice)%TVoisin(r,1)) then
                    condition = .false.
                    exit
                end if
            end do
            if (condition) then
                Mesh%Tnoeud(indice)%NVoisin(2) = Mesh%Tnoeud(indice)%NVoisin(2) + 1				
                Mesh%Tnoeud(indice)%TVoisin(Mesh%Tnoeud(indice)%NVoisin(2),1) = jvoisin(k)
                Mesh%Tnoeud(indice)%TVoisin(Mesh%Tnoeud(indice)%NVoisin(2),2) = ordre
            end if
        end do
            !if (.false.) then
                if (Symmetry .and. Mesh.Tnoeud(indice).Pnoeud(2).lt.Epsilon .and. Mesh.Tnoeud(indice).Normale(2).lt.Epsilon) then
                    jvoisin = -jvoisin
                    do k=1,3*Mesh%Tnoeud(indice)%Nfacette
                        ! On teste les indices des points récupérés pour vérifier qu'on ne les a pas déjà.
                        condition = .true.
                        do r=1,Mesh%Tnoeud(indice)%NVoisin(2)
                            if (jvoisin(k).eq.0 .or. jvoisin(k).eq.Mesh%Tnoeud(indice)%TVoisin(r,1) &
                            .or. Mesh%Tnoeud(-jvoisin(k))%Pnoeud(2).lt.Epsilon) then
                                condition = .false.
                                exit
                            end if
                        end do
                        if (condition) then
                            Mesh%Tnoeud(indice)%NVoisin(2) = Mesh%Tnoeud(indice)%NVoisin(2) + 1				
                            Mesh%Tnoeud(indice)%TVoisin(Mesh%Tnoeud(indice)%NVoisin(2),1) = jvoisin(k)
                            Mesh%Tnoeud(indice)%TVoisin(Mesh%Tnoeud(indice)%NVoisin(2),2) = ordre
                        end if
                    end do
                end if
            !end if
    end do
        
    ! Free surface
    cond = .true.
    condOrdre = .true.
    ordre = 1
    do while(cond)
        ordre = ordre + 1
        cond = .false.
        do indice = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
            if(condOrdre(indice)) then

                Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(indice)%Tvoisin(2,1)))%Pnoeud
                if (Mesh%Tnoeud(indice)%Tvoisin(2,1).lt.Epsilon) Pvo(2) = -Pvo(2)
                Mesh%Tnoeud(indice)%Rvoisin = norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo)
                do j=3,Mesh%Tnoeud(indice)%NVoisin(2)
                    Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(indice)%Tvoisin(j,1)))%Pnoeud
                    if (Mesh%Tnoeud(indice)%Tvoisin(j,1).lt.Epsilon) Pvo(2) = -Pvo(2)
                    Mesh%Tnoeud(indice)%Rvoisin = Mesh%Tnoeud(indice)%Rvoisin + norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo) !min(Mesh%Tnoeud(indice)%Rvoisin,norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo))
                end do
                Mesh%Tnoeud(indice)%Rvoisin = 1._RP*Mesh%Tnoeud(indice)%Rvoisin/float(Mesh%Tnoeud(indice)%Nvoisin(2))

    21          condOrdre(indice) = .false.
                potar_voisin = 20 ! Permet de regler manuellement le nombre de noeuds voisins
                do j=2,Mesh%Tnoeud(indice)%Nvoisin(2)
                    jtemp=Mesh%Tnoeud(indice)%Tvoisin(j,1)
                    modificateur = sign(1,jtemp)
                    jtemp = abs(jtemp)
                    do k=2,Mesh%Tnoeud(jtemp)%Nvoisin(2)
                        condition = .true.
                        do r=1,Mesh%Tnoeud(indice)%NVoisin(2)
                            Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(jtemp)%Tvoisin(k,1)))%Pnoeud
                            if ((modificateur.lt.Epsilon.and.Mesh%Tnoeud(jtemp)%Tvoisin(k,1).gt.Epsilon).or.(modificateur.gt.Epsilon.and.Mesh%Tnoeud(jtemp)%Tvoisin(k,1).lt.Epsilon)) Pvo(2) = - Pvo(2)
                            if (Mesh%Tnoeud(jtemp)%Tvoisin(k,1)*modificateur .eq. Mesh%Tnoeud(indice)%Tvoisin(r,1) .or. Mesh%Tnoeud(jtemp)%Tvoisin(k,2).gt.1 .or.&
                                (modificateur.lt.0 .and. Mesh%Tnoeud(abs(Mesh%Tnoeud(jtemp)%Tvoisin(k,1)))%Pnoeud(2).lt.Epsilon)) then ! .or. & !
                                 !norm2(Pvo-Mesh%Tnoeud(indice)%Pnoeud).gt.(0*RSpline+Mesh%Tnoeud(indice)%Rvoisin)) then ! Permet de prendre les points voisins dans une sphère de rayon donné
                                condition = .false.
                                exit
                            end if
                        end do
                        if (condition) then
                            Mesh%Tnoeud(indice)%Nvoisin(2) = Mesh%Tnoeud(indice)%Nvoisin(2) + 1
                            Mesh%Tnoeud(indice)%Tvoisin(Mesh%Tnoeud(indice)%Nvoisin(2),1) = Mesh%Tnoeud(jtemp)%Tvoisin(k,1)*modificateur
                            Mesh%Tnoeud(indice)%Tvoisin(Mesh%Tnoeud(indice)%Nvoisin(2),2) = ordre
                            condOrdre(indice) = .true.
                        end if
                    end do
                end do
            end if
    22  end do
        if (t_spline) then
            j = Mesh%FS%IndFS(3) + 3
            open(unit=3,file='Voisin.dat')
            write(3,fmt='(50a)') 'Title= "Voisinage"'
            write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z"'
            do k=1,Mesh%Tnoeud(j)%NVoisin(2)
                if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) then
                    Pvo = Mesh%Tnoeud(-Mesh%Tnoeud(j)%TVoisin(k,1))%Pnoeud
                    Pvo(2) = - Pvo(2)
                    write(3,'(13E,I)') Pvo
                else
                    write(3,'(13E,I)') Mesh%Tnoeud(Mesh%Tnoeud(j)%TVoisin(k,1))%Pnoeud
                end if      
            end do
           close(unit=3)
        end if
    end do
        
    ! Bodies
    do nc = Int_Body,Mesh%NBody
        
        if(fgeom_vect%Active(nc-1))then
            cond = .true.
            condOrdre = .true.
            ordre = 1
            do while(cond)
                ordre = ordre + 1
                cond = .false.
                do indice = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    
                    if(condOrdre(indice)) then
                        
                        Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(indice)%Tvoisin(2,1)))%Pnoeud                        
                        if (Mesh%Tnoeud(indice)%Tvoisin(2,1).lt.Epsilon) Pvo(2) = -Pvo(2)
                        Mesh%Tnoeud(indice)%Rvoisin = norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo)
                        do j=3,Mesh%Tnoeud(indice)%NVoisin(2)
                            Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(indice)%Tvoisin(j,1)))%Pnoeud
                            if (Mesh%Tnoeud(indice)%Tvoisin(j,1).lt.Epsilon) Pvo(2) = -Pvo(2)
                            Mesh%Tnoeud(indice)%Rvoisin = Mesh%Tnoeud(indice)%Rvoisin + norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo) !min(Mesh%Tnoeud(indice)%Rvoisin,norm2(Mesh%Tnoeud(indice)%Pnoeud-Pvo))
                        end do
                        Mesh%Tnoeud(indice)%Rvoisin = 1._RP*Mesh%Tnoeud(indice)%Rvoisin/float(Mesh%Tnoeud(indice)%Nvoisin(2))
                        condOrdre(indice) = .false.
                        potar_voisin = 20 ! Permet de regler manuellement le nombre de noeuds voisins
                        do j=2,Mesh%Tnoeud(indice)%Nvoisin(2)
                            jtemp=Mesh%Tnoeud(indice)%Tvoisin(j,1)
                            modificateur = sign(1,jtemp)
                            jtemp = abs(jtemp)
                            do k=2,Mesh%Tnoeud(jtemp)%Nvoisin(2)
                                condition = .true.
                                do r=1,Mesh.Tnoeud(indice).Nvoisin(2)
                                    Pvo = Mesh%Tnoeud(abs(Mesh%Tnoeud(jtemp)%Tvoisin(k,1)))%Pnoeud
                                    if ((modificateur.lt.Epsilon.and.Mesh%Tnoeud(jtemp)%Tvoisin(k,1).gt.Epsilon).or.(modificateur.gt.Epsilon.and.Mesh%Tnoeud(jtemp)%Tvoisin(k,1).lt.Epsilon)) Pvo(2) = - Pvo(2)
                                    if (Mesh%Tnoeud(jtemp)%Tvoisin(k,1)*modificateur .eq. Mesh%Tnoeud(indice)%Tvoisin(r,1) .or. Mesh%Tnoeud(jtemp)%Tvoisin(k,2).gt.1 .or.&
                                        (modificateur.lt.0 .and. Mesh%Tnoeud(abs(Mesh%Tnoeud(jtemp)%Tvoisin(k,1)))%Pnoeud(2).lt.Epsilon)) then ! .or. & !
                                            !norm2(Pvo-Mesh%Tnoeud(indice)%Pnoeud).gt.(0*RSpline+Mesh%Tnoeud(indice)%Rvoisin)) then ! Permet de prendre les points voisins dans une sphère de rayon donné
                                        condition = .false.
                                        exit
                                    end if
                                end do
                                if (condition) then
                                    Mesh%Tnoeud(indice)%Nvoisin(2) = Mesh%Tnoeud(indice)%Nvoisin(2) + 1
                                    Mesh%Tnoeud(indice)%Tvoisin(Mesh%Tnoeud(indice)%Nvoisin(2),1) = Mesh%Tnoeud(jtemp)%Tvoisin(k,1)*modificateur
                                    Mesh%Tnoeud(indice)%Tvoisin(Mesh%Tnoeud(indice)%Nvoisin(2),2) = ordre
                                    condOrdre(indice) = .true.
                                end if
                            end do
                        end do
                                                
                    end if
                    
                end do
        
                if (t_spline) then
                    j = Mesh%FS%IndFS(3) + 3
                    open(unit=3,file='Voisin.dat')
                    write(3,fmt='(50a)') 'Title= "Voisinage"'
                    write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z"'
                    do k=1,Mesh%Tnoeud(j)%NVoisin(2)
                        if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) then
                            Pvo = Mesh%Tnoeud(-Mesh%Tnoeud(j)%TVoisin(k,1))%Pnoeud
                            Pvo(2) = - Pvo(2)
                            write(3,'(13E,I)') Pvo
                        else
                            write(3,'(13E,I)') Mesh%Tnoeud(Mesh%Tnoeud(j)%TVoisin(k,1))%Pnoeud
                        end if      
                    end do
                   close(unit=3)
                end if
            
            end do
        end if
    end do
    
    ! Ajout des points voisins des points doubles ayant la même normale (?)
    do j=1,Mesh%Nnoeud
        if (Mesh%Tnoeud(j)%NDouble.ne.0) then
            do k=1,Mesh%Tnoeud(j)%NDouble
                if (norm2(Mesh%Tnoeud(j)%Normale-Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%Normale).lt.Epsilon) then
                if (j.lt.Mesh%Tnoeud(j)%Double(k)) then
                    jtemp = Mesh%Tnoeud(j)%NVoisin(2)
                    do r=2,Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%NVoisin(2)
                        if (Mesh%Tnoeud(abs(Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%TVoisin(r,1)))%NDouble.eq.0) then
                            Mesh%Tnoeud(j)%NVoisin(2) = Mesh%Tnoeud(j)%NVoisin(2) + 1
                            Mesh%Tnoeud(j)%TVoisin(Mesh%Tnoeud(j)%NVoisin(2),1) = Mesh%TNoeud(Mesh%Tnoeud(j)%Double(k))%TVoisin(r,1)
                        end if
                    end do
                    do r=2,jtemp
                        if (Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%TVoisin(r,1)))%NDouble.eq.0) then
                            Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%NVoisin(2) = Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%NVoisin(2) + 1
                            Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%TVoisin(Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%NVoisin(2),1) = Mesh%TNoeud(j)%TVoisin(r,1)
                        end if
                    end do
                end if
                end if
            end do
        end if
    end do
        
end subroutine CVoisin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                  !
!                       Base Locale et Courbures associées aux vecteurs tangents                   !
!                                 LETOURNEL Lucas    01/14                                         !
!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BLocal(Mesh)
    !!!!! Problème :
    ! Détermination d'une base locale en chaque point du maillage et des courbures associées aux vecteurs
    ! tangents.
    !!!!!
    
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage) :: Mesh
    
    integer :: j, k , r, s, Ordre, Na, Nvoisin, nc, ios
    real(rp), dimension(2) :: DF
    real(rp), dimension(3) :: M, Normale
    real(rp), dimension(3,3) :: Puvw, Pt, Gnor
    real(rp), allocatable :: Pvoisin(:,:), Pk(:,:), B(:), Alpha(:), A(:,:), B1(:,:), Courbure(:,:), EpsNormale(:,:)
    character(len=1) :: trans
    integer :: info, nrhs
    integer, allocatable :: ipiv(:)
    real(rp), dimension(2) :: DF1, DF2
    real(rp),dimension(3),parameter :: E3 = [0._rp, 0._rp, 1._rp]
    
    
    if(idebug>0) print*,"** BLocal"
    
    ! Initialization
    ios = 0
    
    if (is_BS) then
        do nc=Int_Body,Mesh%NBody
            if(Mesh%Body(nc)%Active)then
                allocate(Courbure(2,Mesh%Body(nc)%IndBody(3)-Mesh%Body(nc)%IndBody(1)+1),stat=ios)
                if(ios/=0) goto 9999
                allocate(EpsNormale(2,Mesh%Body(nc)%IndBody(3)-Mesh%Body(nc)%IndBody(1)+1),stat=ios)
                if(ios/=0) goto 9999
                
                if(idebug>0) print*," Calcul de la normale par BSplines"
                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    
                    ! Recalcul de la normale à l'aide des BSplines
                    Ordre = Mesh%Tnoeud(j)%Ordre
                    Na = Nordre(Ordre)
                    NVoisin = Mesh%Tnoeud(j)%NVoisin(2)-1
                    
                    ! If not enough neighbours with the polynomial splines, use of the thin plane ones.
                    if(Nvoisin.le.6 .and. Mesh%Tnoeud(j)%Ordre.gt.0)then
                        !print*,""
                        !print*,"BLocal: not enough neighboors points"
                        !print*,"Pnoeud = ",Mesh%Tnoeud(j)%Pnoeud
                        !print*,"SplineOrder = ",Mesh%Tnoeud(j)%Ordre
                        !print*,"Na = ",Na
                        !print*,"NVoisin = ",NVoisin
                        !print*,"It would be necessary to reduce the size of the panel for the floater ",nc
                        !print*,"%Order of the node ",j," of the body ",nc," is ",0
                        Mesh%Tnoeud(j)%Ordre = 0 ! Thin plate slines.
                        Ordre = Mesh%Tnoeud(j)%Ordre
                        Na = Nordre(Ordre)
                        !print*,"SplineOrder = ",Ordre
                        !print*,"Na = ",Na
                        !print*,"NVoisin = ",NVoisin
                    end if
                    
                    allocate(Pvoisin(3,NVoisin+1), Pk(3,NVoisin+1), B(NVoisin+Na), A(NVoisin+Na,NVoisin+Na), Alpha(1:NVoisin+Na))
                    
                    ! Base Locale à partir de la normale au point, précédemment calculée
                    M = Mesh%Tnoeud(j)%Pnoeud
                    call LocalBase(M, Mesh%Tnoeud(j)%Normale, Puvw)
                    Pt = transpose(Puvw)
                    
                    ! Changement de variables
                    do k=1,NVoisin+1
                        Pk(1:3,k) = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%Pnoeud
                        if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) Pk(2,k)=-Pk(2,k)
                        Pvoisin(1:3,k) = matmul(Pt,Pk(1:3,k) - M)                        
                    end do
                    
                    ! On construit B
                    B(1:Nvoisin+1) = Pvoisin(3,1:Nvoisin+1)
                    B(Nvoisin+2:Nvoisin+Na) = 0._RP
                    
                    ! On construit A
                    call SplineMatrix(Nvoisin, Ordre, Pvoisin, A)
                    
                    ! Résolution du Système Linéaire
                    Alpha = 0._RP
                    call LU(A, B, Alpha, Nvoisin+Na)
                    
                    ! Calcul des Dérivées Premières de la géométrie locale
                    call SplineDF(Nvoisin, Ordre, Alpha, Pvoisin(1:3,1), Pvoisin, DF)
                    Normale = [-DF(1),-DF(2),1._RP]/norm2([-DF(1),-DF(2),1._RP])
                    Normale = matmul(Puvw,Normale)
                    
                    if(Mesh%Tnoeud(j)%Ndouble.ne.0)then
                        Normale = Mesh%Tnoeud(j)%Normale
                    else
                        Mesh%Tnoeud(j)%Normale = Normale
                    endif
                    
                    ! Base Locale à partir de la nouvelle normale
                    call LocalBase(M, Normale, Puvw)
                    Pt = transpose(Puvw)
                    
                    ! Changement de variables
                    do k=1,NVoisin+1
                        Pvoisin(1:3,k) = matmul(Pt,Pk(1:3,k) - M)
                    end do
                    
                    ! On construit B
                    B(1:Nvoisin+1) = Pvoisin(3,1:Nvoisin+1)
                    B(Nvoisin+2:Nvoisin+Na) = 0._RP
                    
                    ! On construit A
                    call SplineMatrix(Nvoisin, Ordre, Pvoisin, A)
                    
                    ! Résolution du Système Linéaire
                    Alpha = 0._RP
                    call LU(A, B, Alpha, Nvoisin+Na)
                    
                    ! Dérivées partielles spatiales 1ères
                    call SplineDF(Nvoisin, Ordre, Alpha, Pvoisin(1:3,1), Pvoisin, DF)
                    Mesh%Tnoeud(j)%DLocal(1:2) = abs(DF(1:2))
                    
                    ! Dérivées partielles spatiales 2ndes
                    Mesh%Tnoeud(j)%DLocal(5) = DF(1)**2 + DF(2)**2 + 1._RP
                    
                    ! Local base and its transpose
                    Mesh%Tnoeud(j).PLocal(1:3,1:3,1) = Puvw
                    Mesh%Tnoeud(j).PLocal(1:3,1:3,2) = Pt
                    
                    ! Deallocating
                    deallocate(Pvoisin, Pk, A, B, Alpha)
                end do
    
                if(idebug>0) print*," derivee locale"
                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    M = Mesh%Tnoeud(j)%Pnoeud
                    Ordre = Mesh%Tnoeud(j)%Ordre
                    Na = Nordre(Ordre)
                    NVoisin = Mesh%Tnoeud(j)%Nvoisin(2)-1
                    allocate(Pk(3,NVoisin+1), Pvoisin(3,NVoisin+1), B1(NVoisin+Na,2), A(NVoisin+Na,NVoisin+Na), ipiv(NVoisin+Na))
                    Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
                    Pt(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,2)
            
                    ! Changement de variables
                    do k=1,NVoisin+1
                        Pk(1:3,k) = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%TVoisin(k,1)))%Pnoeud
                        if (Mesh%Tnoeud(j)%TVoisin(k,1).lt.0) Pk(2,k)=-Pk(2,k)
                        Pvoisin(1:3,k) = matmul(Pt,Pk(1:3,k) - M)
                    end do
            
                    ! On construit A
                    call SplineMatrix(Nvoisin, Ordre, Pvoisin, A)
            
                    ! Passage des Normales des points voisins dans la base locale
                    do k=1,NVoisin+1
                        Normale = Mesh%Tnoeud(abs(Mesh%Tnoeud(j)%Tvoisin(k,1)))%Normale
                        if (Mesh%Tnoeud(j)%Tvoisin(k,1).lt.0) Normale(2)=-Normale(2)
                        Normale = matmul(Pt,Normale) ! tableau des normales dans la base locale
                        B1(k,1:2) = Normale(1:2)
                    end do
                    B1(Nvoisin+2:Nvoisin+Na,1:2) = 0._RP
                    ipiv=0 
                    trans='n'
                    nrhs=2 ! nombre de second membre
            
                    ! Décomposition LU de A
                    call dgetrf(Nvoisin+Na,Nvoisin+Na,A,Nvoisin+Na,ipiv,info)
                    if(idebug>0 .and. info.ne.0) print*,' BLocal, décomposition LU:      info =', info
            
                    ! Résolution du système à partir de la décomposition Lu
                    call dgetrs(trans,Nvoisin+Na,nrhs,A,Nvoisin+Na,ipiv,B1,Nvoisin+Na,info)
                    if(idebug>0 .and. info.ne.0) print*,' BLocal, décomposition LU:      info =', info
                    DF1 = 0._RP
                    DF2 = 0._RP
                    call SplineDF(Nvoisin, Ordre, B1(:,1), Pvoisin(1:3,1), Pvoisin, DF1)
                    call SplineDF(Nvoisin, Ordre, B1(:,2), Pvoisin(1:3,1), Pvoisin, DF2)
                    Mesh%Tnoeud(j)%DLocal(3:4) = [DF1(1),DF2(2)]
                    deallocate(Pk,Pvoisin, A, B1, ipiv)
                end do
            
                ! Deallocating
                if(allocated(Courbure)) deallocate(Courbure)
                if(allocated(EpsNormale)) deallocate(EpsNormale)
                
            end if
        end do
    else
        do nc = Int_Body,Mesh%NBody
            if(Mesh%Body(nc)%Active)then
                do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                    M = Mesh%Tnoeud(j)%Pnoeud
                    call LocalBase(M, Mesh%Tnoeud(j)%Normale, Puvw)
                    Mesh%Tnoeud(j)%PLocal(1:3,1:3,1) = Puvw
                    Mesh%Tnoeud(j)%PLocal(1:3,1:3,2) = transpose(Puvw)
                    ! Mauvais calcul à reprendre...
                    do s=1,2
                        Gnor = 0._RP
                        do k=1,Mesh%Tnoeud(j)%Nfacette
                            do r=1,3
                                Normale(r) = Mesh%Tnoeud(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Tnoeud(r))%Normale(s)
                            end do
                            Gnor(1:3,s) = Gnor(1:3,s) + matmul(Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%ds,Normale)*0.5_RP*Mesh%Tnoeud(j)%Angle(k)
                        end do
                    end do
                    Mesh%Tnoeud(j)%DLocal(3) = dot_product(matmul(Gnor,Mesh%Tnoeud(j)%Plocal(1:3,1,1)),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                    Mesh%Tnoeud(j)%DLocal(4) = dot_product(matmul(Gnor,Mesh%Tnoeud(j)%Plocal(1:3,2,1)),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                    Mesh%Tnoeud(j)%DLocal(3:4) = 1._RP/Mesh%Body(nc)%DimBody(1) !! Sphere uniquement...
                    Mesh%Tnoeud(j)%DLocal(1:2) = 0._RP
                    Mesh%Tnoeud(j)%DLocal(5) =  1._RP
                end do
            end if
        end do
    endif
    
    ! Previously not commented
    
    !if (t_spline) then
    !    print*,"BLocal: No MB in this part of the routine."
    !    open(unit=3,file='DLocal.dat', iostat=ios)
    !    if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
    !    write(3,fmt='(50a)') 'Title= "Remaillage de la cuve"'
    !    write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","Du","Dv","D2u","D2v","Ux","Uy","Uz","Vx","Vy","Vz","NVoisin"' !
    !    write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(1)%IndBody(3)-Mesh%Body(1)%IndBody(1)+1, ', E=', Mesh%Body(1)%IndBody(4)-Mesh%Body(1)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
    !    do j=Mesh%Body(1)%IndBody(1),Mesh%Body(1)%IndBody(3)
    !        write(3,'(13E,I)') Mesh%Tnoeud(j)%Pnoeud,  Mesh%Tnoeud(j)%DLocal(1:4), Mesh%Tnoeud(j)%Plocal(1:3,1,1), Mesh%Tnoeud(j)%Plocal(1:3,2,1), Mesh%Tnoeud(j)%NVoisin(2)
    !    end do
    !    do j=Mesh%Body(1)%IndBody(2),Mesh%Body(1)%IndBody(4)
    !        write(3,'(3I)') Mesh%Tfacette(j)%Tnoeud - (Mesh%Body(1)%IndBody(1)-1)*[1,1,1]
    !    end do 
    !    close(unit=3)
    !
    !    do j=Mesh%Body(1)%IndBody(1),Mesh%Body(1)%IndBody(3)
    !        NormaleT = (Mesh%Tnoeud(j)%Pnoeud-Mesh%Body(1)%GBody(1:3))/Mesh%Body(1)%DimBody(1)
    !        EpsNormale(2,j) = 100._RP*norm2(Mesh%Tnoeud(j)%Normale - NormaleT)/norm2(NormaleT)
    !        Courbure(1,j)=100._RP*abs(Mesh%Body(1)%DimBody(1) - abs(1._RP/Mesh%Tnoeud(j)%DLocal(3)))/Mesh%Body(1)%DimBody(1)
    !        Courbure(2,j)=100._RP*abs(Mesh%Body(1)%DimBody(1) - abs(1._RP/Mesh%Tnoeud(j)%DLocal(4)))/Mesh%Body(1)%DimBody(1)
    !    end do
    !    open(unit=3,file='Courbure.dat', iostat=ios)
    !    if (ios/=0) stop "Erreur à l'ouverture du fichier de gradient"
    !    write(3,fmt='(50a)') 'Title= "Remaillage de la cuve"'
    !    write(3,fmt='(50a)') 'VARIABLES = "X","Y","Z","Epsk1","Epsk2","EpsL","EpsBS","NVoisin"'
    !    write(3,fmt='(a,i,a,i,a)') 'Zone N =', Mesh%Body(1)%IndBody(3)-Mesh%Body(1)%IndBody(1)+1, ', E=', Mesh%Body(1)%IndBody(4)-Mesh%Body(1)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT' 
    !    do j=Mesh%Body(1)%IndBody(1),Mesh%Body(1)%IndBody(3)
    !        write(3,'(7E,I)') Mesh%Tnoeud(j)%Pnoeud, Courbure(1:2,j), EpsNormale(1:2,j), Mesh%Tnoeud(j)%NVoisin(2)
    !    end do
    !    do j=Mesh%Body(1)%IndBody(2),Mesh%Body(1)%IndBody(4)
    !        write(3,'(3I)') Mesh%Tfacette(j)%Tnoeud - (Mesh%Body(1)%IndBody(1)-1)*[1,1,1]
    !    end do 
    !    close(unit=3)
    !    deallocate(Courbure, EpsNormale)
    !end if
    
    9999 continue
    
    if(ios/=0)then
        print*,"BLocal: allocating error."
    endif
        
    if(idebug>0) print*," ** Fin BLocal"
    
end subroutine BLocal 

subroutine LocalBase(M, Normale, Puvw)
    implicit none
    real(rp),dimension(3),intent(in) :: M
    real(rp),dimension(3),intent(in) :: Normale
    real(rp),dimension(3,3),intent(inout) :: Puvw
    
    call LocalBase_gen(M, Normale, Puvw)

end subroutine LocalBase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                      Local Base                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LocalBase_gen(M, Normal, Puvw)
    
    real(rp), dimension(3), intent(in)      :: M, Normal    ! Point and normal to the panel
    real(rp), dimension(3,3), intent(out)   :: Puvw         ! Local base
    
    integer                                 :: k            ! Loop parameter
    real(rp), dimension(3)                  :: M1           ! Point
    
    ! This subroutine computes the local base from the normal.
    
    ! Normal 
    Puvw(1:3,3) = Normal(1:3)
    
    ! u
    M1 = [random(1),random(2),random(3)]
    if (norm2(M).lt.Epsilon) then
        M1 = M1/norm2(M1)
    else
        M1 = M1/norm2(M1)*norm2(M)
    end if
    k = 1
    do while (1._RP-abs(dot_product((M-M1),Puvw(1:3,3))/norm2((M-M1))).lt.Epsilon)
        k=k+10
        M1 = [random(k),random(k+1),random(k+2)]
        M1 = M1/norm2(M1)*norm2(M)
    end do
    Puvw(1:3,1) = (M-M1) - dot_product((M-M1),Puvw(1:3,3))*Puvw(1:3,3)
    Puvw(1:3,1) = Puvw(1:3,1)/norm2(Puvw(1:3,1))
    
    ! v
    call Computation_vect_product(Puvw(1:3,3),Puvw(1:3,1),Puvw(1:3,2))
    
end subroutine LocalBase_gen

subroutine Updating_is_piercing(Maillage,fgeom_vect,InputData)
        
    !f2py integer*1, dimension(1000)    :: Maillage
    type(TMaillage),intent(inout)       :: Maillage         ! Mesh.
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(in)      :: fgeom_vect       ! Geometries
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(inout) :: InputData        ! InputData.
    
    integer                             :: nc,j,jk,k        ! Loop parameters.
    integer                             :: NumBody          ! Number of a body.
    integer                             :: int_sl0, int_sl1 ! 1st and last node of the free surface mesh.
    integer                             :: Ndouble          ! Number of twin nodes.
    logical                             :: Bool             ! Boolean to know if a node is on the intersection curves or not.
            
    ! This subroutine updates the values of Mesh%ody%Active and is_piercing for every body.
    
    int_sl0 = Maillage%FS%IndFS(1)
    int_sl1 = Maillage%FS%IndFS(3)
    
    ! Looking for the nodes at the intersection between the free surface and the piercing bodies
    do nc = 1,NBodies
        
        NumBody = nc + 1 ! 0 is FS, 1 is tank.
        Bool = .false.
                
        ! Loop over the nodes of the free surface
        do j = int_sl0,int_sl1
            
            Ndouble = Maillage%Tnoeud(j)%Ndouble
            
            ! Nodes of the free surface in intersection with the floaters
            do k = 1,Ndouble
                jk = Maillage%Tnoeud(j)%double(k)
                Bool = Bool .or. Maillage%Tnoeud(jk)%Npanneau .eq. NumBody ! Intersection FS - Floaters
            enddo
            
            ! If the body pierces the free surface, moving to 10.
            if(Bool)then
                go to 10
            end if
            
        end do
        
        10 continue
        
        if(fgeom_vect%Active(nc))then
            Maillage%Body(NumBody)%Active = .true.
            ! Bool = True if a body number NumBody is piercing the FS, False otherwise.
            if(Bool)then
                InputData%is_piercing(nc) = 1 ! Piercing body.
            else
                InputData%is_piercing(nc) = 0 ! Immerged body.
            end if
        else
            Maillage%Body(NumBody)%Active = .false.
            InputData%is_piercing(nc) = 2 ! Body above the free surface.
            Maillage%Body(NumBody)%CMD = [.false.,.false.] ! No deformation.
        end if
    end do
    
end subroutine Updating_is_piercing

end module GeomMesh
