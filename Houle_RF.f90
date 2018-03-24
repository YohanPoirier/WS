module Houle_RF
use Parameters
use Constantes
use Structuresdonnees
use iso_c_binding
implicit none

! Include components of the Rienecker & Fenton wave model

!f2py integer*1, dimension(1000) :: HouleRF
type(THouleRF) :: HouleRF

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_EtaRF(M, t,EtaRF)
    
    real(rp),intent(out)                :: EtaRF    ! Wave elevation.
    real(rp),intent(in)                 :: t        ! Current time.
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y     ! x and y coordinates.
    
    ! This subroutine computes the wave elevation for RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    EtaRF = HouleRF%An(0)/2
    do j=1,HouleRF%N
        EtaRF = EtaRF + HouleRF%An(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)
    end do
    
end subroutine Computation_EtaRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_DEtaDtRF(M, t,DEtaDtRF)
    
    real(rp),intent(out)                :: DEtaDtRF ! Time-differentiation of the wave elevation.
    real(rp),intent(in)                 :: t        ! Current time.
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y     ! x and y coordinates.
    
    ! This subroutine computes the time-differentiation of the wave elevation for RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    DEtaDtRF = 0._RP
    do j=1,HouleRF%N
        DEtaDtRF = DEtaDtRF + j*HouleRF%k*HouleRF%C*HouleRF%An(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)
    end do
    
end subroutine Computation_DEtaDtRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GEtaRF(M, t, GEta)
    
    real(rp), dimension(3), intent(in)  :: M        ! Point of computation.
    real(rp), intent(in)                :: t        ! Current time.
    real(rp), dimension(3), intent(out) :: GEta     ! Gradient of the wave elevation.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the gradient of the wave elevation for RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    GEta(1:3) = 0._RP
    do j=1,HouleRF%N
        GEta(1) = GEta(1) - j*HouleRF%k*cos(DirectionRF)*HouleRF%An(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)
        GEta(2) = GEta(2) - j*HouleRF%k*sin(DirectionRF)*HouleRF%An(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)
    end do
    
end subroutine GEtaRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_PhiRF(M, t,PhiRF)
    
    real(rp),intent(out)                :: PhiRF    ! Potential of the incident RF waves.
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    real(rp),intent(in)                 :: t        ! Current time.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    PhiRF = (HouleRF%Bn(0) + HouleRF%C)*(x+y)
    do j=1,HouleRF%N
        PhiRF = PhiRF + HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
end subroutine Computation_PhiRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_DPhiDtRF(M, t,DPhiDtRF)
    
    real(rp),intent(out)                :: DPhiDtRF     ! Time-differentiation of the potential of the incident RF waves.
    real(rp), dimension(3),intent(in)   :: M            ! Point of computation.
    real(rp),intent(in)                 :: t            ! Current time.
    
    integer                             :: j            ! Loop parameter.
    real(rp)                            :: x, y, z      ! x, y and z coordinates.
    
    ! This subroutine computes the time-differentiation of the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    DPhiDtRF = 0._RP
    do j=1,HouleRF%N
        DPhiDtRF = DPhiDtRF - j*HouleRF%k*HouleRF%C*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
    
end subroutine Computation_DPhiDtRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GPhiRF(M, t, GPhi)
    
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    real(rp), intent(in)                :: t        ! Current time.
    real(rp), dimension(3), intent(out) :: GPhi     ! Gradient of the potential of the incident RF waves.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the gradient of the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    GPhi(1) = (HouleRF%Bn(0) + HouleRF%C)*cos(DirectionRF)
    GPhi(2) = (HouleRF%Bn(0) + HouleRF%C)*sin(DirectionRF)
    GPhi(3) = 0._RP
    do j=1,HouleRF%N
        GPhi(1) = GPhi(1) + j*HouleRF%k*cos(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        GPhi(2) = GPhi(2) + j*HouleRF%k*sin(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        GPhi(3) = GPhi(3) + j*HouleRF%k*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
    
end subroutine GPhiRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GPhi2RF(M, t, GPhi2, GradGrad)    !!! Problème sur GradGrad !!! A verifier absolument !!!
    
    real(rp), dimension(3),intent(in)       :: M        ! Point of computation.
    real(rp), intent(in)                    :: t        ! Current time.
    real(rp), dimension(3), intent(out)     :: GPhi2    ! Square gradient.
    real(rp), dimension(3,3), intent(out)   :: GradGrad ! Gradient of the gradient?
    
    integer                                 :: j        ! Loop parameter.
    real(rp)                                :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the second order spatial differentiation of the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    GPhi2(1:3) = 0._RP
    GradGrad(1:3,1:3) = 0._RP
    do j=1,HouleRF%N
        GPhi2(1) = GPhi2(1) - (j*HouleRF%k*cos(DirectionRF))**2*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        GPhi2(2) = GPhi2(2) - (j*HouleRF%k*sin(DirectionRF))**2*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        GPhi2(3) = GPhi2(3) + (j*HouleRF%k)**2*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        GradGrad(1,2) = GradGrad(1,2) - (j*HouleRF%k)**2*cos(DirectionRF)*sin(DirectionRF)*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
	    GradGrad(1,3) = GradGrad(1,3) + (j*HouleRF%k)**2*cos(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
	    GradGrad(2,3) = GradGrad(2,3) + (j*HouleRF%k)**2*sin(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
    GradGrad(1,1) = Gphi2(1)
    GradGrad(2,2) = Gphi2(2)
    GradGrad(3,3) = Gphi2(3)
    GradGrad(2,1) = GradGrad(1,2)
    GradGrad(3,2) = GradGrad(2,3)
    GradGrad(3,1) = GradGrad(1,3)
    
end subroutine GPhi2RF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DGPhiDzRF(M, t, DGPhiDz)
    
    real(rp), dimension(3),intent(in)   :: M        ! Point of Computation.
    real(rp),intent(in)                 :: t        ! Current time.
    real(rp), dimension(3), intent(out) :: DGPhiDz  ! Gradient of the vertical differentiation of the potential of the incident RF waves.
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the grdient of the vertical differentiation of the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    DGPhiDz(1:3) = 0._RP
    do j=1,HouleRF%N
        DGPhiDz(1) = DGPhiDz(1) + (j*HouleRF%k)**2*cos(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        DGPhiDz(2) = DGPhiDz(2) + (j*HouleRF%k)**2*sin(DirectionRF)*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        DGPhiDz(3) = DGPhiDz(3) + (j*HouleRF%k)**2*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
    
end subroutine DGPhiDzRF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DGPhiDtRF(M, t, DGPhiDt)
    
    real(rp), dimension(3),intent(in)   :: M        ! Point of computation.
    real(rp),intent(in)                 :: t        ! Current time.
    real(rp), dimension(3), intent(out) :: DGPhiDt  ! Gradient of the time differentiation of the potential of the incident RF waves
    
    integer                             :: j        ! Loop parameter.
    real(rp)                            :: x, y, z  ! x, y and z coordinates.
    
    ! This subroutine computes the gradient of the time differentiation of the potential of the incident RF waves.
    
    x=M(1)*cos(DirectionRF)
    y=M(2)*sin(DirectionRF)
    z=M(3)
    DGPhiDt = 0._RP
    do j=1,HouleRF%N
        DGPhiDt(1) = DGPhiDt(1) + (j*HouleRF%k)**2*HouleRF%C*cos(DirectionRF)*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        DGPhiDt(2) = DGPhiDt(2) + (j*HouleRF%k)**2*HouleRF%C*sin(DirectionRF)*HouleRF%Bn(j)*sin(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*cosh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
        DGPhiDt(3) = DGPhiDt(3) - (j*HouleRF%k)**2*HouleRF%C*HouleRF%Bn(j)*cos(j*HouleRF%k*(x+y - HouleRF%C*t) + PhiWaveRF)*sinh(j*HouleRF%k*(z+HouleRF%d))/cosh(j*HouleRF%k*HouleRF%d)
    end do
    
end subroutine DGPhiDtRF

end module Houle_RF