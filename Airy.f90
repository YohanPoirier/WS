module Airy
use parameters
use Constantes
use Structuresdonnees
use FonctionsCommunes   
implicit none

contains

!--------------------------------------------------
!
!	Calcul de la deformee de surface libre associee 
!	au champ incident au point M et à l'instant t
!
!--------------------------------------------------
subroutine Computation_EtaAiry(M,t,EtaAiry)
	    
    real(rp),intent(out)                :: EtaAiry  ! Wave elevation.
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    real(rp),intent(in)                 :: t        ! Current time.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: wbar     ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: phase    ! j*(k*x-Omega*t).
        
    ! This subroutine computes the wave elevation for Airy waves.
    
    EtaAiry = 0._RP
    do j=1,Nhoule
	    wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	    phase = i*(konde(j)*wbar-w(j)*t + PhiWave(j))
	    EtaAiry = EtaAiry+real(Aphi(j)*EXP(phase))
    end do
    
end subroutine Computation_EtaAiry

!--------------------------------------------------
!
!	Calcul de la derivee temporelle de la
!	deformee de surface libre associee 
!	au champ incident au point (X,Y) et à l'instant t
!
!--------------------------------------------------
subroutine Computation_DEtaDtAiry(M,t,DEtaDtAiry)
    
    real(rp)                :: DEtaDtAiry   ! Time-differentiation of the wave elevation.
    real(rp),dimension(3)   :: M            ! Point of computation.
    real(rp)                :: t            ! Current time.
    
    integer                 :: j            ! Loop parameter.
    real(rp)                :: wbar         ! x*cos(beta) + y*sin(beta).
    COMPLEX                 :: phase        ! j*(k*x-Omega*t).
    
    ! This subroutine computes the time-differentiation of the wave elevation for Airy waves.
    
    DEtaDtAiry=0._RP
	do j=1,Nhoule
		wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
		phase = i*(konde(j)*wbar-w(j)*t + PhiWave(j))
		DEtaDtAiry = DEtaDtAiry+real(-Aphi(j)*i*w(j)*EXP(phase))
	end do	
    
end subroutine Computation_DEtaDtAiry

!--------------------------------------------------
!
!	Calcul du gradient de la deformee de surface libre 
!	associee au champ incident au point (X,Y) 
!	et à l'instant t
!
!--------------------------------------------------
subroutine GEtaAiry(M, t, Geta)
    
    real(rp),dimension(3)   :: Geta     ! Gradient of the wave elevation.
    real(rp),dimension(3)   :: M        ! Point of computation.
    real(rp)                :: t        ! Current time.
    
    integer                 :: j        ! Loop parameter.
    real(rp)                :: wbar     ! x*cos(beta) + y*sin(beta).
    COMPLEX                 :: CG,phase !
            
    ! This subroutine computes the gradient of the wave elevation for Airy waves.
    
    Geta(1:3)=0._RP
	do j=1,Nhoule
		wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
		phase  =i*(konde(j)*wbar-w(j)*t + PhiWave(j))
		CG = i*konde(j)*Aphi(j)*EXP(phase)
		Geta(1) = Geta(1)+real(cos(DirectionAiry(j))*CG)
		Geta(2) = Geta(2)+real(sin(DirectionAiry(j))*CG)
	end do
    
end subroutine GEtaAiry
!
!--------------------------------------------------
!
!	Calcul du champ incident au point M à l'instant t
!
!--------------------------------------------------
subroutine Computation_PhiAiry(M,t,Depth,PhiAiry)
    
    real(rp),intent(out)                :: PhiAiry  ! Potential of the incident Airy waves.
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    real(rp),intent(in)                 :: t, Depth ! Current time and water depth.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: wbar     ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: phase    ! Phase.
    
    ! This subroutine computes the potential of the incident Airy waves.
    
    PhiAiry=0._RP
    if(infinite_depth)then
        do j=1,Nhoule
	        wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	        phase = konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
	        PhiAiry = PhiAiry+real(-i*Aphi(j)*g/w(j)*EXP(phase))
        enddo
    else
        do j=1,Nhoule
	        wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	        phase  = i*(konde(j)*wbar-w(j)*t + PhiWave(j))
	        PhiAiry = PhiAiry+real(-i*Aphi(j)*g/w(j)*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)*exp(phase))
        enddo
    endif
    
end subroutine Computation_PhiAiry
!
!--------------------------------------------------
!
!	Calcul de la derivee temporelle du
!	champ incident au point M et à l'instant t
!
!--------------------------------------------------
subroutine Computation_DPhiDtAiry(M,t,Depth,DPhiDtAiry)
    
    real(rp),intent(out)                :: DPhiDtAiry   ! Time-differentiation of the potential of the incident Airy waves.
    real(rp),dimension(3),intent(in)    :: M            ! Point of computation.
    real(rp),intent(in)                 :: t, Depth     ! Current time and water depth.
    
    integer                             :: j            ! Loop parameter.
    real(rp)                            :: wbar         ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: phase        ! Phase.
        
    ! This subroutine computes the time-differentiation of the potential of the incident Airy waves.
    
    DPhiDtAiry=0._RP
    if(infinite_depth)then
        do j=1,Nhoule
            wbar=M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
            phase=konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
            dPhiDTAiry=dPhiDTAiry+real(-Aphi(j)*g*EXP(phase))
        enddo
    else
        do j=1,Nhoule
	        wbar=M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	        phase = i*(konde(j)*wbar-w(j)*t + PhiWave(j))
            dPhidtAiry=dPhidtAiry+real(-Aphi(j)*g*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)*exp(phase))
        enddo
    endif
end subroutine Computation_DPhiDtAiry

!--------------------------------------------------
!
!	Calcul du gradient du champ incident au point M
!	et à l'instant t
!
!--------------------------------------------------
subroutine GphiAiry(M, t, Depth, Gphi0)
    
    real(rp),dimension(3), intent(out)  :: Gphi0        ! Gradient of the potential of the incident Airy waves.
    real(rp),dimension(3), intent(in)   :: M            ! Point of computation.
    real(rp), intent(in)                :: t, Depth     ! Current time and water depth.
    
    integer                             :: j            ! Loop parameter.
    real(rp)                            :: wbar         ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: CG,phase,CG2 !
    
    ! This subroutine computes the gradient of the potential of the incident Airy waves.
    
    Gphi0(:) = 0._RP
    if(infinite_depth)then	
        do j = 1,Nhoule
	        wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	        phase = konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
	        CG = i*konde(j)*Aphi(j)*g/w(j)*EXP(phase)	
	        Gphi0(1) = Gphi0(1) + real(-i*cos(DirectionAiry(j))*CG)
	        Gphi0(2) = Gphi0(2) + real(-i*sin(DirectionAiry(j))*CG)
	        Gphi0(3) = Gphi0(3) + real(-CG)
        end do
    else
        do j = 1,Nhoule
	        wbar = M(1)*cos(DirectionAiry(j)) + M(2)*sin(DirectionAiry(j))
	        phase = i*(konde(j)*wbar-w(j)*t + PhiWave(j))
            CG  = i*konde(j)*Aphi(j)*g/w(j)*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)*exp(phase)
            CG2 = i*konde(j)*Aphi(j)*g/w(j)*sinh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)*exp(phase)
	        Gphi0(1) = Gphi0(1) + real(-i*cos(DirectionAiry(j))*CG)
	        Gphi0(2) = Gphi0(2) + real(-i*sin(DirectionAiry(j))*CG)
	        Gphi0(3) = GPhi0(3) + real(-CG2)
        end do
    end if
    
end subroutine GphiAiry
!
!--------------------------------------------------
!
!	Calcul des dérivées secondes du champ incident
!   au point M et à l'instant t
!
!--------------------------------------------------
subroutine GPhi2Airy(M, t, Depth, Gphi2, GradGrad)
    
    real(rp), dimension(3), intent(in)              :: M                ! Point of computation.
    real(rp), intent(in)                            :: t, Depth         ! Current time and water depth.
    real(rp), dimension(3), intent(out)             :: Gphi2            ! Square gradient.
    real(rp), dimension(3,3), intent(out), optional :: GradGrad         ! Gradient of the gradient?
    
    integer                                         :: j                ! Loop parameter.
    real(rp)                                        :: wbar             ! x*cos(beta) + y*sin(beta).
    COMPLEX                                         :: CG, CG2, phase   !
    
    ! This subroutine computes the second order spatial differentiation of the potential of the incident Airy waves.
    
    Gphi2(:)=0._RP
    GradGrad=0._RP
    if(infinite_depth)then	
        do j=1,Nhoule
	        wbar=M(1)*cos(DirectionAiry(j)) + M(2)*sin(DirectionAiry(j))
	        phase=konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
	        CG=i*konde(j)**2*Aphi(j)*g/w(j)*EXP(phase)
	        Gphi2(1) = Gphi2(1) + real((cos(DirectionAiry(j)))**2*CG)
	        Gphi2(2) = Gphi2(2) + real((sin(DirectionAiry(j)))**2*CG)
	        Gphi2(3) = Gphi2(3) + real(-CG)
	        GradGrad(1,2) = GradGrad(1,2) + real(cos(DirectionAiry(j))*sin(DirectionAiry(j))*CG)
	        GradGrad(1,3) = GradGrad(1,3) + real(-i*cos(DirectionAiry(j))*CG)
	        GradGrad(2,3) = GradGrad(2,3) + real(-i*sin(DirectionAiry(j))*CG)
        end do
    else
        do j=1,Nhoule
	        wbar=M(1)*cos(DirectionAiry(j)) + M(2)*sin(DirectionAiry(j))
            phase=i*(konde(j)*wbar-w(j)*t + PhiWave(j))
            CG = i*konde(j)*konde(j)*Aphi(j)*g/w(j)*EXP(phase)*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)
            CG2 =  konde(j)*konde(j)*APhi(j)*g/w(j)*EXP(phase)*sinh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)
	        Gphi2(1) = Gphi2(1) + real((cos(DirectionAiry(j)))**2*CG)
	        Gphi2(2) = Gphi2(2) + real((sin(DirectionAiry(j)))**2*CG)
	        Gphi2(3) = Gphi2(3) + real(-CG)
	        GradGrad(1,2) = GradGrad(1,2) + real(cos(DirectionAiry(j))*sin(DirectionAiry(j))*CG)
	        GradGrad(1,3) = GradGrad(1,3) + real(cos(DirectionAiry(j))*CG2)
	        GradGrad(2,3) = GradGrad(2,3) + real(sin(DirectionAiry(j))*CG2)
        end do
    end if
    
    GradGrad(1,1) = Gphi2(1)
    GradGrad(2,2) = Gphi2(2)
    GradGrad(3,3) = Gphi2(3)
    GradGrad(2,1) = GradGrad(1,2)
    GradGrad(3,2) = GradGrad(2,3)
    GradGrad(3,1) = GradGrad(1,3)
    
end subroutine GPhi2Airy

!--------------------------------------------------
!
!	Calcul du gradient de la derivee verticale 
!	du champ incident au point M 
!	et à l'instant t
!
!--------------------------------------------------
subroutine DGPhiDzAiry(M, t, Depth, DGPhiDz)
    
    real(rp), dimension(3), intent(in)  :: M            ! Point of Computation.
    real(rp), intent(in)                :: t, Depth     ! Current time and water depth.
    real(rp), dimension(3), intent(out) :: DGPhiDz      ! Gradient of the vertical differentiation of the potential of the incident Airy waves.
    
    integer                             :: j            ! Loop parameter.
    real(rp)                            :: wbar         ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: CG,CG2,phase
    
    ! This subroutine computes the grdient of the vertical differentiation of the potential of the incident Airy waves.
    
    DGPhiDz(:)=0._RP
    if(infinite_depth)then
        do j=1,Nhoule
	    wbar=M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	    phase=konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
	    CG=-i*konde(j)*Aphi(j)*g/w(j)*EXP(phase)
	    DGPhiDz(1)=DGPhiDz(1)+real(i*cos(DirectionAiry(j))*CG)
	    DGPhiDz(2)=DGPhiDz(2)+real(i*sin(DirectionAiry(j))*CG)
	    DGPhiDz(3)=DGPhiDz(3)+real(CG)
        enddo
    else
        do j=1,Nhoule
	    wbar=M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	    phase=i*(konde(j)*wbar-w(j)*t + PhiWave(j))
        CG=i*konde(j)*konde(j)/w(j)*Aphi(j)*g*EXP(phase)*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)
        CG2=i*konde(j)*konde(j)*APhi(j)*g/w(j)*EXP(phase)*sinh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)
        DGPhiDz(1) = DGPhiDz(1)+real(-i*cos(DirectionAiry(j))*CG2)
        DGPhiDz(2) = DGPhiDz(2)+real(-i*sin(DirectionAiry(j))*CG2)
	    DGPhiDz(3)=DGPhiDz(3)+real(-CG)
        enddo
    endif
    
end subroutine DGPhiDzAiry
!
!--------------------------------------------------
!
!	Calcul du gradient de la derivee temporelle du
!	champ incident au point M et à l'instant t
!
!--------------------------------------------------
subroutine DGPhiDtAiry(M, t, Depth, GDPhiDt)
    
    real(rp), dimension(3), intent(in)  :: M                ! Point of computation.
    real(rp), intent(in)                :: t, Depth         ! Current time and water depth.
    real(rp), dimension(3), intent(out) :: GDPhiDt          ! Gradient of the time differentiation of the potential of the incident Airy waves.
    
    integer                             :: j                ! Loop parameter.
    real(rp)                            :: wbar             ! x*cos(beta) + y*sin(beta).
    COMPLEX                             :: CG, CG2, phase   !
        
    ! This subroutine computes the gradient of the time differentiation of the potential of the incident Airy waves.
    
    GDPhiDt(:) = 0._RP
    if(infinite_depth)then
        do j=1,Nhoule
	        wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
	        phase = konde(j)*(M(3)+i*wbar)-i*w(j)*t+i*PhiWave(j)
	        CG = w(j)*w(j)*Aphi(j)*EXP(phase)
	        GDPhiDt(1) = GDPhiDt(1)+real(-i*cos(DirectionAiry(j))*CG)
	        GDPhiDt(2) = GDPhiDt(2)+real(-i*sin(DirectionAiry(j))*CG)
	        GDPhiDt(3) = GDPhiDt(3)+real(-CG)
        end do
    else
        do j=1,Nhoule
            wbar = M(1)*cos(DirectionAiry(j))+M(2)*sin(DirectionAiry(j))
            phase = i*(konde(j)*wbar-w(j)*t + PhiWave(j))            
            CG = -konde(j)*g*Aphi(j)*EXP(phase)*cosh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)            
            CG2 = -konde(j)*g*APhi(j)*EXP(phase)*sinh(konde(j)*(M(3)+Depth))/cosh(konde(j)*Depth)
            GDPhiDt(1) = GDPhiDt(1) + real(i*cos(DirectionAiry(j))*CG)
            GDPhiDt(2) = GDPhiDt(2) + real(i*sin(DirectionAiry(j))*CG)
            GDPhiDt(3) = GDPhiDt(3) + real(CG2)
        end do
    end if
    
end subroutine DGPhiDtAiry

end module Airy