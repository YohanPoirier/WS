module sys_lineaire
    use FonctionsCommunes
    implicit none

contains


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
    
    print*, Mesh%Nnoeud, Mesh%Nfacette
    print("B BVP :" ,N)
    pause
    
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

end module