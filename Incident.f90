module Incident_mod
use parameters
use Constantes
use StructuresDonnees
use Airy
use Houle_RF
implicit none
contains

subroutine Incident(t, Mesh, Ecoulement)
    
    real(rp), intent(in)                :: t                                        ! Current time
    !f2py integer*1, dimension(1000)    :: Mesh
    type(TMaillage), intent(in)         :: Mesh                                     ! Mesh
    !f2py integer*1, dimension(1000)    :: Ecoulement
    type(TEcoulement), intent(inout)    :: Ecoulement                               ! Flow parameters
    
    integer                             :: j, nc                                    ! Loop parameters
    real(rp) , dimension(3)	            :: M, GPhi, GPhi2, GEta, DGPhiDz, GdPhidt   ! 
    real(rp), dimension(3,3)            :: GradGrad                                 !
    real(rp)                            :: EtaAiry,DEtaDtAiry,PhiAiry,DPhiDtAiry    ! Airy quantities
    real(rp)                            :: EtaRF,DEtaDtRF,PhiRF,DPhiDtRF            ! RF quantities
    
    ! This subroutine initializes the incident flow parameters.
    
    if(idebug>0) then
      print*,' CCL : t = ',t,' Htype = ',Htype,' Ind_sl0 = ',Mesh%FS%IndFS(1),&
      & ' Ind_sl1 = ',Mesh%FS%IndFS(3)
    endif
    
    select case(Htype)
        case(0)
            ! Still water
            
            ! Eta
            Ecoulement%Eta(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DEtaDt(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%GEta(1:3,1:Mesh%Nsys)%incident = 0._RP
            
            ! Phi
            Ecoulement%Phi(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DPhiDt(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%GPhi(1:3,1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%GPhi2(1:4,1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%GsDPhiDn(1:3,1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DGPhiDz(1:3,1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DPhiDn(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DGradPhiSqDn(1:Mesh%Nsys)%incident = 0._RP
            Ecoulement%DDPhiDnDt(1:Mesh%Nsys)%incident = 0._RP
            
        case(1)
            ! Airy
            
            ! Free surface
            do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
                
                M = Mesh%Tnoeud(j)%Pnoeud
                                
                ! Eta
                call Computation_EtaAiry(M,t,EtaAiry)
                call Computation_DEtaDtAiry(M,t,DEtaDtAiry)
                Ecoulement%Eta(j)%incident = EtaAiry
                Ecoulement%DEtaDt(j)%incident = DEtaDtAiry
                call GEtaAiry(M, t, GEta)
                Ecoulement%GEta(:,j)%incident = GEta
                
                ! Phi
                call Computation_PhiAiry(M,t,profondeur,PhiAiry)
                call Computation_DPhiDtAiry(M,t,profondeur,DPhiDtAiry)
                Ecoulement%Phi(j)%incident = PhiAiry
                Ecoulement%DPhiDt(j)%incident = DPhiDtAiry                
                call GphiAiry(M, t, profondeur, Gphi)
                Ecoulement%GPhi(:,j)%incident = Gphi
                Ecoulement%DPhiDn(j)%incident = dot_product(GPhi,Mesh%Tnoeud(j)%Normale)
                call DGPhiDtAiry(M, t, profondeur, GDPhiDt)
                Ecoulement%DDPhiDnDt(j)%incident = dot_product(GDPhiDt,Mesh%Tnoeud(j)%Normale)
                call Gphi2Airy(M, t, profondeur, Gphi2, GradGrad)
                Ecoulement%GPhi2(:,j)%incident = Gphi2
                Ecoulement%GradGrad(:,:,j)%incident = GradGrad
                call DGPhiDzAiry(M, t, profondeur, DGPhiDz)
                Ecoulement%DGPhiDz(:,j)%incident = DGPhiDz
            end do
            
            ! Tank and bodies
            do nc = 1,Mesh%NBody
                if(Mesh%Body(nc)%Active)then
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                        
                        M = Mesh%Tnoeud(j)%Pnoeud
                    
                        ! Phi
                        call Computation_PhiAiry(M,t,profondeur,PhiAiry)
                        call Computation_DPhiDtAiry(M,t, profondeur,DPhiDtAiry)
                        Ecoulement%Phi(j)%incident = PhiAiry
                        Ecoulement%DPhiDt(j)%incident = DPhiDtAiry                        
                        call GphiAiry(M, t, profondeur, Gphi)
                        Ecoulement%GPhi(:,j)%incident = Gphi
                        Ecoulement%DPhiDn(j)%incident = dot_product(GPhi,Mesh%Tnoeud(j)%Normale)
                        call DGPhiDtAiry(M, t, profondeur, GDPhiDt)
                        Ecoulement%DDPhiDnDt(j)%incident = dot_product(GDPhiDt,Mesh%Tnoeud(j)%Normale)
                        call Gphi2Airy(M, t, profondeur, Gphi2, GradGrad)
                        Ecoulement%GradGrad(:,:,j)%incident = GradGrad
                        Ecoulement%GPhi2(1,j)%incident =dot_product(Gphi, Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                        Ecoulement%GPhi2(2,j)%incident =dot_product(Gphi, Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                        Ecoulement%GPhi2(3,j)%incident = dot_product(matmul(GradGrad,Mesh%Tnoeud(j)%Plocal(1:3,1,1)),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                        Ecoulement%GPhi2(4,j)%incident = dot_product(matmul(GradGrad,Mesh%Tnoeud(j)%Plocal(1:3,2,1)),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                        Ecoulement%GPhi2(3:4,j)%incident = Ecoulement%GPhi2(3:4,j)%incident - Mesh%TNoeud(j)%DLocal(3:4)*Ecoulement%DPhiDn(j)%incident
                        Ecoulement%DGPhiDz(:,j)%incident = GradGrad(1:3,3)
                    end do
                end if
            end do
            
        case(2)
            ! RF
            
            ! Free surface
            do j = Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
                
                M = Mesh%Tnoeud(j)%Pnoeud
                
                ! Eta
                call Computation_EtaRF(M, t,EtaRF)
                call Computation_DEtaDtRF(M, t,DEtaDtRF)
                Ecoulement%Eta(j)%incident = EtaRF
                Ecoulement%DEtaDt(j)%incident = DEtaDtRF
                call GEtaRF(M, t, GEta)
                Ecoulement%GEta(:,j)%incident = GEta
                
                ! Phi
                call Computation_PhiRF(M,t,PhiRF)
                call Computation_DPhiDtRF(M, t,DPhiDtRF)
                Ecoulement%Phi(j)%incident = PhiRF
                Ecoulement%DPhiDt(j)%incident = DPhiDtRF
                call GphiRF(M, t, Gphi)
                Ecoulement%GPhi(:,j)%incident = Gphi
                Ecoulement%DPhiDn(j)%incident = dot_product(GPhi,Mesh%Tnoeud(j)%Normale)
                call DGPhiDtRF(M,t,GDPhiDt)
                Ecoulement%DDPhiDnDt(j)%incident = dot_product(GDPhiDt,Mesh%Tnoeud(j)%Normale)
                call Gphi2RF(M, t, Gphi2, GradGrad)
                Ecoulement%GPhi2(:,j)%incident = Gphi2
                Ecoulement%GradGrad(:,:,j)%incident = GradGrad
                call DGPhiDzRF(M, t, DGPhiDz)
                Ecoulement%DGPhiDz(:,j)%incident = DGPhiDz
            end do
            
            ! Tank and bodies
            do nc = 1,Mesh%NBody
                if(Mesh%Body(nc)%Active)then
                    do j = Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                        
                        M = Mesh%Tnoeud(j)%Pnoeud
                    
                        ! Phi
                        call Computation_PhiRF(M,t,PhiRF)
                        call Computation_DPhiDtRF(M, t,DPhiDtRF)
                        Ecoulement%Phi(j)%incident = PhiRF
                        Ecoulement%DPhiDt(j)%incident = DPhiDtRF
                        call GphiRF(M, t, Gphi)
                        Ecoulement%GPhi(:,j)%incident = Gphi
                        Ecoulement%DPhiDn(j)%incident = dot_product(GPhi,Mesh%Tnoeud(j)%Normale)
                        call DGPhiDtRF(M,t,GDPhiDt)
                        Ecoulement%DDPhiDnDt(j)%incident = dot_product(GDPhiDt,Mesh%Tnoeud(j)%Normale)
                        call Gphi2RF(M, t, Gphi2, GradGrad)!
                        Ecoulement%GradGrad(:,:,j)%incident = GradGrad
                        Ecoulement%GPhi2(1,j)%incident = dot_product(Gphi, Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                        Ecoulement%GPhi2(2,j)%incident = dot_product(Gphi, Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                        Ecoulement%GPhi2(3,j)%incident = dot_product(matmul(GradGrad,Mesh%Tnoeud(j)%Plocal(1:3,1,1)),Mesh%Tnoeud(j)%Plocal(1:3,1,1))
                        Ecoulement%GPhi2(4,j)%incident = dot_product(matmul(GradGrad,Mesh%Tnoeud(j)%Plocal(1:3,2,1)),Mesh%Tnoeud(j)%Plocal(1:3,2,1))
                        Ecoulement%GPhi2(3:4,j)%incident = Ecoulement%GPhi2(3:4,j)%incident - Mesh%TNoeud(j)%DLocal(3:4)*Ecoulement%DPhiDn(j)%incident
                        call DGPhiDzRF(M, t, DGPhiDz)
                        Ecoulement%DGPhiDz(:,j)%incident = GradGrad(1:3,3)
                    end do
                end if
            end do
    end select
    
end subroutine Incident

subroutine CEta0(M, t, Eta0)
    
    real(rp)                :: Eta0 ! Wave elevation
    real(rp), dimension(3)  :: M    ! Point of computation
    real(rp)                :: t    ! Current time
    
    ! This subroutine computes the wave elevation of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            Eta0 = 0._RP
        case(1)
            ! Airy
            call Computation_EtaAiry(M,t,Eta0)
        case(2)
            ! RF
            call Computation_EtaRF(M, t,Eta0)
    end select

    if(lineaireFS)then
        Eta0 = 0._rp
    endif
    
end subroutine CEta0

subroutine CDEta0Dt(M, t, DEta0Dt)

    real(rp)                :: DEta0Dt  ! Time differentiation of the wave elevation
    real(rp), dimension(3)  :: M        ! Point of computation
    real(rp)                :: t        ! Current time
    
    ! This subroutine computes the time differentiation of the wave elevation of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            DEta0Dt = 0._RP
        case(1)
            ! Airy
            call Computation_DEtaDtAiry(M,t,DEta0Dt)
        case(2)
            ! RF
            call Computation_DEtaDtRF(M, t, DEta0Dt)
    end select
    
end subroutine CDEta0Dt

subroutine CGEta0(M,t,GEta0)
    ! Paramètres :
    real(rp), dimension(3)  :: M        ! Point of computation
    real(rp)                :: t        ! Current time
    real(rp), dimension(3)  :: GEta0    ! Gradient of the wave elevation
    
    ! This subroutine computes the gradient of the wave elevation of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            GEta0 = 0._RP
        case(1)
            ! Airy
            call GEtaAiry(M, t, GEta0)
        case(2)
            ! RF
            call GEtaRF(M, t, GEta0)
    end select

    if(lineaireFS)then
        GEta0 = 0._rp
    endif
    
end subroutine CGETa0

subroutine CGPhi0(M, t, GPhi0)
    
    real(rp), dimension(3)  :: M        ! Point
    real(rp)                :: t        ! Current time
    real(rp), dimension(3)  :: GPhi0    ! Gradient of Phi_incident
        
    ! This subroutine computes the gradient of the potential of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            GPhi0 = [0._rp, 0._rp, 0._rp]
        case(1)
            ! Airy
            call GPhiAiry(M, t, profondeur, GPhi0)
        case(2)
            ! RF
            call GPhiRF(M, t, GPhi0)
    end select

    if(lineaireFS)then
        GPhi0 = 0._rp
    endif
    
end subroutine CGPhi0

subroutine CGDPhi0Dt(M, t, GdPhi0dt)

    real(rp), dimension(3)  :: M        ! Point of computation
    real(rp)                :: t        ! Current time
    real(rp), dimension(3)  :: GdPhi0dt ! Gradient of the time differentiation of the potential
    
    ! This subroutine computes the gradient of the time differentiation of the potential of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            GdPhi0dt = [0._rp, 0._rp, 0._rp]
        case(1)
            ! Airy
            call DGPhiDtAiry(M, t, profondeur, GdPhi0dt)
        case(2)
            ! RF
            call DGPhiDtRF(M, t, GdPhi0dt)
    end select

end subroutine CGDPhi0Dt

subroutine CD2Phi0DzDt(M, t, D2Phi0DzDt)
    
    real(rp), dimension(3)  :: M            ! Point of computation
    real(rp)                :: t            ! Current time
    real(rp)                :: D2Phi0DzDt   ! Vertical differentiation of the time differentiation of the potential
    
    real(rp), dimension(3)  :: GdPhi0dt     ! Gradient of the time differentiation of the potential
    
    ! This subroutine computes the vertical differentiation of the time differentiation of the potential of the incident waves.
    
    select case(Htype)
        case(0)
            ! Still water
            GdPhi0dt = [0._rp, 0._rp, 0._rp]
        case(1)
            ! Airy
            call DGPhiDtAiry(M, t, profondeur, GdPhi0dt)
        case(2)
            ! RF
            call DGPhiDtRF(M, t, GdPhi0dt)
    end select
    D2Phi0DzDt = GdPhi0dt(3)
    
end subroutine CD2Phi0DzDt

end module Incident_mod