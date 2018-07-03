module Validation_DPhiDn
    use Constantes
    use Parameters
    use Structuresdonnees
    use CoeffInfluence
    use FonctionsCommunes
    use BVP
    use Spline
    use BodyMotion_mod
    implicit none
    contains
    
    subroutine BT_DPhiDn(nt, t, Mesh)
    ! Parameters
    integer, intent(in) :: nt
    real(rp), dimension(nt), intent(in) :: t
    type(TMaillage) :: Mesh
    ! Locals
    integer :: jt
    real(rp) :: CCIP, tmoy, time_begin
    real(rp), dimension(nt) :: time_end
    real(rp), allocatable :: CD(:,:), CS(:,:)
    type(TEcoulement) :: Ecoulement
    ! Begin
    open(unit=38, file='VerifDPhiDt_'//filename)
    write(38,fmt='(50a)') 'Title= "Verification du potentiel sur le corps"'
    write(38,fmt='(50a)') 'VARIABLES = "X","Y","Z","DPhiDtNum","DPhiDtAna","ErreurDPhiDt","DDPhiDnDtNum","DDPhiDnDtAna","ErreurDDPhiDnDt","QNum","Qana","ErreurQ","Acc""Vel"'
    close(38)
    call NewEcoulement(Ecoulement, Mesh%Nnoeud)
    call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)      
    allocate(CD(Mesh%Nnoeud,Mesh%Nnoeud), CS(Mesh%Nnoeud,Mesh%Nnoeud))
    call cpu_time(time_begin)
    CD=0._RP
    CS=0._RP
    call CoeffInfl(Mesh, CD, CS)
    call cpu_time(time_end(1))
    print*, 'Temps CI :', time_end(1)-time_begin
    call cpu_time(time_begin)
    CCIP = 1._RP
    call BodyVel(Mesh,0._RP,0._RP)
    call BodyNode(Mesh)
    call Remesh_DDPhiDnDt(Mesh, Mesh, t(1))
    do jt=1,nt
        if (not(lineaire)) CCIP = -t(jt)
        ! Initialisation BVP
        call InitBVP_DPhiDt(Mesh, Ecoulement)
        ! Résolution BVP
        call SolBVP(Ecoulement, Mesh, CD, CS, CCIP, .true.)
        ! Verification Solution BVP
        call Verif_DPhiDt(Mesh, Ecoulement, t(jt))
        ! Calcul position, vitesse, acceleration du centre du corps
        call BodyVel(Mesh,t(jt+1),t(jt+1)-t(jt))
        ! Calcul position, vitesse, acceleration des noeuds du corps
        call BodyNode(Mesh)
        if (not(lineaire)) call Remesh_DDPhiDnDt(Mesh, Mesh, t(jt+1))
        call cpu_time(time_end(jt))    
        if (jt==1) then
            tmoy = 2*(time_end(jt) - time_begin)/3._RP
            print*, 'Le temps de l''operation a ete de ', time_end(jt) - time_begin, ' secondes'
        else
            tmoy = (tmoy*(jt-1) + time_end(jt) - time_end(jt-1))/jt
            print*, 'Le temps de l''operation a ete de ', time_end(jt) - time_end(jt-1), ' secondes'
        end if
        print*, 'Fin prevue dans ', tmoy*(size(t)-jt)/60._RP, ' minutes'
    end do
    
    deallocate(CD,CS)
    call DelEcoulement(Ecoulement)
    ! End
    end subroutine BT_DPhiDn
    
    subroutine Verif_DPhiDt(Mesh, Ecoulement, t)
    ! Parameters
    type(TMaillage) :: Mesh
    type(TEcoulement) :: Ecoulement
    real(rp) :: t
    ! Locals
    character(len=10) :: num
    integer :: j
    real(rp) :: DPhiDtN, DPhiDtA, erreurDPhiDt
    real(rp) :: DDPhiDNDtN, DDPhiDNDtA, erreurDDPhiDNDt
    real(rp) :: QN, QA, erreurQ
    ! Begin
    open(unit=38, file='VerifDPhiDt_'//filename, position="append")
    Write( num, '( f0.4 )' ) t
    write(38,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Body(Int_Body)%IndBody(3)+(-Mesh%Body(Int_Body)%IndBody(1)+1), ', E=', Mesh%Body(Int_Body)%IndBody(4)+(-Mesh%Body(Int_Body)%IndBody(2)+1), ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    !allocate(Phi(Mesh%Body(Int_Body)%IndBody(1):Mesh%Body(Int_Body)%IndBody(3))
    do j=Mesh%Body(Int_Body)%IndBody(1), Mesh%Body(Int_Body)%IndBody(3)
        if (.true.) then
            ! DPhiDt
            DPhiDtN = Ecoulement%DPhiDt(j)%Perturbation
            DPhiDtA = Ecoulement%DPhiDt(j)%Incident
            erreurDPhiDt = 100._RP*abs(DPhiDtN - DPhiDtA)
            ! DDPhiDnDt
            DDPhiDnDtN = Ecoulement%DDPhiDnDt(j)%Perturbation
            DDPhiDnDtA = Ecoulement%DDPhiDnDt(j)%Incident
            erreurDDPhiDNDt = 100._RP*abs(DDPhiDNDtN - DDPhiDNDtA)
            ! Q
            QN = DDPhiDnDtN - dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale)
            QA = DDPhiDnDtA - dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale)
            erreurQ = 100._RP*abs(QN - QA)
        else
            ! Phi
            DPhiDtN = Ecoulement%Phi(j)%Perturbation
            DPhiDtA = Ecoulement%Phi(j)%Incident
            erreurDPhiDt = 100._RP*abs(DPhiDtN - DPhiDtA)
            ! DPhiDn
            DDPhiDnDtN = Ecoulement%DPhiDn(j)%Perturbation
            DDPhiDnDtA = Ecoulement%DPhiDn(j)%Incident
            erreurDDPhiDNDt = 100._RP*abs(DDPhiDnDtN - DDPhiDnDtA)
            ! Q
            QN = DDPhiDnDtN - dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale)
            QA = DDPhiDnDtA - dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale)
            erreurQ = 100._RP*abs(QN - QA)
        end if
        ! Output
        write(38,'(20E)') Mesh%Tnoeud(j)%Pnoeud(1:3), DPhiDtN, DPhiDtA, erreurDPhiDt, DDPhiDnDtN, DDPhiDnDtA, erreurDDPhiDNDt, QN, QA, erreurQ &
        , dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale), dot_product(Mesh%Tnoeud(j)%Velocity, Mesh%Tnoeud(j)%Normale)
    end do
    do j=Mesh%Body(1)%IndBody(2),Mesh%Body(1)%IndBody(4)
        write(38,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(1)%IndBody(1)-1)*[1,1,1]
    end do 
    close(38)
    ! End
    end subroutine Verif_DPhiDt

    subroutine InitBVP_DPhiDt(Mesh,Ecoulement)
    ! Parameters
    type(TMaillage) :: Mesh
    type(TEcoulement) :: Ecoulement
    ! Locals
    integer :: j, nc
    real(rp), dimension(3) :: GPhiLocal
    real(rp), dimension(3,3) :: Puvw
    ! Begin
    do nc=1,Mesh%Nbody
        do j=Mesh%Body(nc)%IndBody(1), Mesh%Body(nc)%IndBody(3)
            Ecoulement%Phi(j)%incident = -0.5*Mesh%Body(nc)%DimBody(1)*dot_product(Mesh%Body(nc)%VBody(1:3),Mesh%TNoeud(j)%Normale)
            Ecoulement%DPhiDn(j)%incident = dot_product(Mesh%Tnoeud(j)%Velocity,Mesh%Tnoeud(j)%Normale)
            Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
            GPhiLocal(1) = -0.5_RP*dot_product(Mesh%Tnoeud(j)%Velocity, Puvw(1:3,1)) ! dot_product(s1,z)
            GPhiLocal(2) = 0._RP
            GPhiLocal(3) = Ecoulement%DPhiDn(j)%incident
            Ecoulement%GPhi(1:3,j)%incident = matmul(Puvw,GPhiLocal)
            Ecoulement%DGradPhiSqDn(j)%incident = 3._RP/Mesh%Body(nc)%DimBody(1)*dot_product(Mesh%TNoeud(j)%Velocity,Mesh%TNoeud(j)%Normale)**2 &
                                                    - 1.5_RP/Mesh%Body(nc)%DimBody(1)*dot_product(Mesh%TNoeud(j)%Velocity,Puvw(1:3,1))**2
            Ecoulement%DDPhiDnDt(j)%incident = dot_product(Mesh%Tnoeud(j)%Acceleration,Mesh%Tnoeud(j)%Normale) + 0*Ecoulement%DGradPhiSqDn(j)%incident
            Ecoulement%DPhiDt(j)%incident = -0.5*Mesh%Body(nc)%DimBody(1)*dot_product(Mesh%Body(nc)%ABody(1:3), Mesh%TNoeud(j)%Normale) - 0*dot_product(Mesh%Body(nc)%VBody(1:3), Ecoulement%GPhi(:,j)%incident)
            if (.false.) then
                Ecoulement%Phi(j)%perturbation = Ecoulement%Phi(j)%incident
            elseif(.false.) then
                Ecoulement%DPhiDn(j)%perturbation = Ecoulement%DPhiDn(j)%incident
            elseif(.false.) then
                Ecoulement%DPhiDt(j)%perturbation = Ecoulement%DPhiDt(j)%incident
            elseif(.true.) then
                Ecoulement%DDPhiDnDt(j)%perturbation = Ecoulement%DDPhiDnDt(j)%incident
            end if
        end do
    end do
    ! End
    end subroutine InitBVP_DPhiDt
    
    subroutine Remesh_DDPhiDnDt(Mesh, Mesh0, t)
    ! Parameters
    type(TMaillage) :: Mesh0, Mesh
    real(rp) :: t
    ! Locals
    integer :: nc, j
    !Begin
    if (is_body) then
        do nc=Int_Body,Int_Body
            do j=Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud + Mesh%Body(nc)%MBody(1:3) !+ vect_product(Mesh%Body(nc)%MBody(4:6),Mesh0%Tnoeud(j)%Pnoeud(1:3)-Mesh0%Body(nc)%GBody(1:3))
            end do
        end do
    end if
    ! End
    end subroutine Remesh_DDPhiDnDt
    
    
    !subroutine VerifPhi(Mesh, Ecoulement, t)
    !! Parameters
    !type(TMaillage) :: Mesh
    !type(TEcoulement) :: Ecoulement
    !real(rp) :: t
    !! Variables
    !character(len=10) :: num
    !integer :: j
    !real(rp) :: cost, erreurPhi, Phi
    !!real(rp), allocatable :: 
    !! Begin
    !open(unit=38, file='VerifPhi_'//filename, position="append")
    !Write( num, '( f0.4 )' ) t
    !write(38,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Body(Int_Body)%IndBody(3)-Mesh%Body(Int_Body)%IndBody(1)+1, ', E=', Mesh%Body(Int_Body)%IndBody(4)-Mesh%Body(Int_Body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    !!allocate(Phi(Mesh%Body(Int_Body)%IndBody(1):Mesh%Body(Int_Body)%IndBody(3))
    !do j=Mesh%Body(Int_Body)%IndBody(1), Mesh%Body(Int_Body)%IndBody(3)
    !    !cost = (Mesh%Body(Int_Body)%GBody(3)-Mesh%TNoeud(j)%Pnoeud(3)) / norm2(Mesh%Body(Int_Body)%GBody-Mesh%TNoeud(j)%Pnoeud)
    !    !Phi = 0.5_RP*Mesh%Body(Int_Body)%DimBody(1)*Mesh%Body(Int_Body)%VBody(3)*cost
    !    Phi = -0.5*Mesh%Body(Int_Body)%DimBody(1)*dot_product(Mesh%Body(Int_Body)%VBody(1:3),Mesh%TNoeud(j)%Normale)
    !    erreurPhi = 100._RP*abs(Ecoulement%Phi(j)%Perturbation - Phi)
    !    !if (abs(Phi).lt.Epsilon) then
    !    !    erreurPhi = abs(Ecoulement%Phi(j)%Perturbation)
    !    !else
    !    !    erreurPhi = 100._RP*abs(Ecoulement%Phi(j)%Perturbation - Phi)/abs(Phi)
    !    !end if
    !    write(38,'(6E)') Mesh%Tnoeud(j)%Pnoeud(1:3), Ecoulement%Phi(j)%Perturbation, Phi, erreurPhi
    !    !Ecoulement%Phi(j)%Perturbation = Phi
    !end do
    !do j=Mesh%Body(1)%IndBody(2),Mesh%Body(1)%IndBody(4)
    !    write(38,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(1)%IndBody(1)-1)*[1,1,1]
    !end do 
    !close(38)
    !!deallocate(Phi)
    !! End
    !end subroutine VerifPhi
    !
    !subroutine VerifDPhiDt(Mesh, Ecoulement, t)
    !! Parameters
    !type(TMaillage) :: Mesh
    !type(TEcoulement) :: Ecoulement
    !real(rp) :: t
    !! Variables
    !character(len=10) :: num
    !integer :: j
    !real(rp) :: cost, erreurDPhiDt, DPhiDt, erreurDDPhiDNDt, DDPhiDNDt, erreurQ, Q
    !!real(rp), allocatable :: 
    !! Begin
    !open(unit=38, file='VerifDPhiDt_'//filename, position="append")
    !Write( num, '( f0.4 )' ) t
    !write(38,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =', Mesh%Body(Int_Body)%IndBody(3)-Mesh%Body(Int_Body)%IndBody(1)+1, ', E=', Mesh%Body(Int_Body)%IndBody(4)-Mesh%Body(Int_Body)%IndBody(2)+1, ' , ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    !!allocate(Phi(Mesh%Body(Int_Body)%IndBody(1):Mesh%Body(Int_Body)%IndBody(3))
    !do j=Mesh%Body(Int_Body)%IndBody(1), Mesh%Body(Int_Body)%IndBody(3)
    !    !DPhiDt = 0.5*Mesh%Body(Int_Body)%ABody(3)*(Mesh%Body(Int_Body)%GBody(3) - Mesh%TNoeud(j)%Pnoeud(3)) - dot_product(Mesh%Body(Int_Body)%VBody, Ecoulement%GPhi(:,j)%Perturbation)
    !    DPhiDt = -0.5*Mesh%Body(Int_Body)%DimBody(1)*dot_product(Mesh%Body(Int_Body)%ABody(1:3), Mesh%TNoeud(j)%Normale) - 0*dot_product(Mesh%Body(Int_Body)%VBody(1:3), Ecoulement%GPhi(:,j)%Perturbation)    
    !    !DPhiDt = -0.5*dot_product(Mesh%Body(Int_Body)%ABody(1:3), Mesh%TNoeud(j)%Normale) - 0.5_RP*dot_product(Mesh%Body(Int_Body)%VBody(1:3), Mesh%Tnoeud(j)%Plocal(1:3,1,1))**2 - dot_product(Mesh%Body(Int_Body)%VBody(1:3), Mesh%TNoeud(j)%Normale)**2    
    !    erreurDPhiDt = 100._RP*abs(Ecoulement%DPhiDt(j)%Perturbation - DPhiDt)
    !    DDPhiDnDt = Ecoulement%DDPhiDnDt(j)%Incident
    !    erreurDDPhiDNDt = 100._RP*abs(Ecoulement%DDPhiDnDt(j)%Perturbation - DDPhiDNDt)
    !
    !    Q = DDPhiDnDt - dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale)
    !    erreurQ = 100._RP*abs(Ecoulement%DGradPhiSqDn(j)%perturbation - Q)
    !    write(38,'(20E)') Mesh%Tnoeud(j)%Pnoeud(1:3), Ecoulement%DPhiDt(j)%Perturbation, DPhiDt, erreurDPhiDt, &
    !    Ecoulement%DDPhiDnDt(j)%Perturbation, DDPhiDnDt, erreurDDPhiDNDt, &
    !    Ecoulement%DGradPhiSqDn(j)%perturbation, Q, erreurQ &
    !    , dot_product(Mesh%Tnoeud(j)%Acceleration, Mesh%Tnoeud(j)%Normale), dot_product(Mesh%Tnoeud(j)%Velocity, Mesh%Tnoeud(j)%Normale)
    !    !Ecoulement%DPhiDt(j)%Perturbation = DPhiDt
    !end do
    !do j=Mesh%Body(1)%IndBody(2),Mesh%Body(1)%IndBody(4)
    !    write(38,'(3I)') Mesh%Tfacette(j)%TNoeud - (Mesh%Body(1)%IndBody(1)-1)*[1,1,1]
    !end do 
    !close(38)
    !!deallocate(Phi)
    !! End
    !end subroutine VerifDPhiDt
    
end module Validation_DPhiDn
    