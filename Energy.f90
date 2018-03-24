module EnergyVolume
use Constantes
use Structuresdonnees
use FonctionsCommunes
use Airy
use Incident_mod
use BodyMotion_mod
implicit none
contains

! ---------------------------------------------------------------------------
!   Energy : interface entre Energy Pseudo Linéaire et Energy Weak-Scatterer
! ---------------------------------------------------------------------------

!subroutine Energy(Mesh, Ecoulement, t, EtSL0, EtSL, iwrite, EtBb)
!    implicit none
!    !f2py integer*1, dimension(1000) :: Mesh
!    type(TMaillage),intent(in) :: Mesh  
!    !f2py integer*1, dimension(1000) :: Ecoulement
!    type(TEcoulement),intent(in) :: Ecoulement
!    real(rp),intent(in) :: t
!    real(rp),intent(in) :: EtSL0
!    real(rp),intent(inout) :: EtSL
!    logical,intent(in) :: iwrite
!    real(rp),intent(in) :: EtBb
!!
!    if(Htype.eq.0 .and. .not.Tcase.eq.4)then
!        call Energy_PseudoLin(Mesh, Ecoulement, t, EtSL0, EtSL, iwrite, EtBb)
!    else
!        call Energy_WeakScatterer(Mesh, Ecoulement, t, EtSL0, EtSL, iwrite, EtBb)
!    endif
!!
!end subroutine Energy

! --------------------------------------------------------------------------
! Energy_PseudoLin : Calcul de l'energie totale du système avec une approche
!                    pseudo-linéaire
! --------------------------------------------------------------------------

subroutine Energy_PseudoLin(Mesh, Ecoulement, t, EtSL0, EtSL, iwrite, EtBb)
    implicit none
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh  
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
    real(rp),intent(in) :: t
    real(rp),intent(in) :: EtSL0
    real(rp),intent(inout) :: EtSL
    logical,intent(in) :: iwrite
    real(rp),intent(in) :: EtBb
!   local
    integer :: j
    integer :: ind_sl0, ind_sl1, ind_body0, ind_body1
    real(rp) :: Phi, DPhiDn, Dn
    real(rp) :: Eta
    real(rp) :: EcSL, EpSL, EcSM 
    real(rp),dimension(3) :: M
 
! -- Surface libre
   
    ind_sl0 = Mesh%FS%IndFS(1)
    ind_sl1 = Mesh%FS%IndFS(3)
    
    EcSL = 0._rp
    EpSL = 0._rp
    do j=ind_sl0,ind_sl1   
      M = Mesh%Tnoeud(j)%Pnoeud     
      Dn = Mesh%Tnoeud(j)%Aire
      DPhiDn = Ecoulement%DPhiDn(j)%perturbation
      Phi = Ecoulement%Phi(j)%perturbation
      Eta = Ecoulement%Eta(j)%perturbation
      EcSL = EcSL - Phi*DPhiDn*Dn
      EpSL = EpSL + Eta**2*Dn
    enddo

    EcSL = 0.5_rp*ro*EcSL 
    EpSL = 0.5_rp*g*ro*EpSL

! -- Surface du flotteur

    ind_body0 = Mesh%Body(Int_Body)%IndBody(1)
    ind_body1 = Mesh%Body(Int_Body)%IndBody(3)

    EcSM = 0._rp
    do j=ind_body0,ind_body1
      M = Mesh%Tnoeud(j)%Pnoeud
      Dn = Mesh%Tnoeud(j)%Aire
      DPhiDn = Ecoulement%DPhiDn(j)%perturbation
      Phi = Ecoulement%Phi(j)%perturbation
      EcSM = EcSM - Phi*DPhiDn*Dn
    enddo

    EcSM = 0.5_rp*ro*EcSM
    EtSL = EcSL + EcSM

! -- Symetrie

    if(Symmetry)then
      EtSL = 2._rp*EtSL
    endif

! -- Energie totale 

    EtSL = EtSL - EtSL0

! -- Ecriture dans fichier

    if(iwrite) write(ioenergy,'(12E)') t, EtSL, EcSL, EpSL, EcSM, EtBb, Mesh%Body(Int_Body)%VBody(1:3), Mesh%Body(Int_Body)%GBody(1:3)

end subroutine Energy_PseudoLin

! ----------------------------------------------------------------------
! Energy_WeakScatterer : Calcul de l'energie totale du systeme avec une 
!                        approche Weak-Scatterer
! ----------------------------------------------------------------------

subroutine Energy_WeakScatterer(Mesh, Ecoulement, t,EtSL0,EtSL,iwrite,EtBb)
! Parameters
  !f2py integer*1, dimension(1000) :: Mesh
  type(TMaillage) :: Mesh
  !f2py integer*1, dimension(1000) :: Ecoulement
  type(TEcoulement) :: Ecoulement
  real(rp), intent(in) :: t
  real(rp), intent(in) :: EtSL0
  real(rp), intent(inout) :: EtSL
  logical, intent(in) :: iwrite
  real(rp), optional :: EtBb
! Locals
  integer :: j, nc
  real(rp) :: EcSL, EpSL, EcB, EcFl, EpFl
  real(rp) :: Dn, Phi,DPhiDn
  real(rp),dimension(3) :: M
  real(rp),parameter :: inv3 = 1._rp/3._rp
 

  EcSL = 0._RP
  EcFl = 0._RP
  EpSL = 0._RP
  EpFl = 0._rp

! -- Surface libre

  do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
    M = Mesh%Tnoeud(j)%Pnoeud
    Dn = Mesh%Tnoeud(j)%Aire
    DPhiDn = Ecoulement%DPhiDn(j)%perturbation + Ecoulement%DPhiDn(j)%incident
    Phi    = Ecoulement%Phi(j)%perturbation + Ecoulement%Phi(j)%incident
    EcSL = EcSL - 0.5_RP*ro*DPhiDn*Phi*Dn
    EpSL = EpSL - 0.5_rp*ro*g*M(3)**2*Dn*Mesh%Tnoeud(j)%Normale(3)     
  end do

! --- Surface du corps

  do nc=Int_Body,Int_Body
    do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
      M = Mesh%Tnoeud(j)%Pnoeud
      Dn = Mesh%Tnoeud(j)%Aire
      DPhiDn = Ecoulement%DPhiDn(j)%perturbation + Ecoulement%DPhiDn(j)%incident
      Phi    = Ecoulement%Phi(j)%perturbation + Ecoulement%Phi(j)%incident
      EcFl = EcFl - 0.5_RP*ro*Phi*DPhiDn*Dn
    end do
  end do

! -- Symetrie

  if(Symmetry)then
    EcSL = 2._rp*EcSL
    EpSL = 2._rp*EpSL
    EcFL = 2._rp*EcFL
    EpFl = 2._rp*EpFL
  endif

! -- Energie totale

  EcB = 0.5_rp*Mesh%Body(Int_Body)%IBody(1,1)*dot_product(Mesh%Body(nc)%VBody(1:3),Mesh%Body(nc)%VBody(1:3)) ! a modifier pour le cas general
  EtSL = EcSL + EpSL + EcFl - EtSL0

! --- Ecriture sur sortie fichier

  if(iwrite) write(ioenergy,'(12E)') t, EtSL, EcSL, EpSL, EcFl, EtBb, Mesh%Body(Int_Body)%VBody(1:3), Mesh%Body(Int_Body)%GBody(1:3)

end subroutine Energy_WeakScatterer

! ----------------------------------------------------------------------------
! FPuissance : interface entre FPuissance Pseudo Linéaire et FPuissance
!              WeakScatterer
! ----------------------------------------------------------------------------
subroutine FPuissance(Mesh, Ecoulement, t, Puissance)
    implicit none
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(inout) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
    real(rp),intent(in) :: t
    real(rp),intent(inout) :: Puissance
!   
    if(Htype.eq.0 .and. Tcase.eq.3)then
        call FPuissance_PseudoLin(Mesh, Ecoulement, Puissance)
    else
        call FPuissance_WeakScatterer(Mesh, Ecoulement, t, Puissance)
    endif
!
end subroutine FPuissance

! ----------------------------------------------------------------------------
! FPuissance_WeakScatterer : calcul de la variation d'énergie dans le domaine
!                            avec une approche WeakScatterer
! ----------------------------------------------------------------------------
subroutine FPuissance_WeakScatterer(Mesh, Ecoulement, t, Puissance)
    implicit none
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(inout) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in)  :: Ecoulement
    real(rp),intent(in) :: t
    real(rp),intent(inout) :: Puissance
!   local
    integer :: j, k
    integer :: Ind_sl0, Ind_sl1
    integer :: Ind_body0, Ind_body1
    integer :: Ind_sf0, Ind_sf1
    real(rp) :: W_WSC, W_SL0, W_Sigma
    real(rp) :: DPhiDt, DPhi0Dn, DPhipDn, GPhi2, Dn, D2Phi0DzDt
    real(rp) :: Etap, TimeRamp
    real(rp), dimension(3) :: GPhi0, GPhip, DS, FHydro, M, DGPhi0Dz, Fr
    real(rp), dimension(3) :: FHydro0, FHydro1, FHydro2
    real(rp), parameter :: inv3 = 1._rp/3._rp

    Ind_sl0 = Mesh%FS%IndFS(1)
    Ind_sl1 = Mesh%FS%IndFS(3)  
       
    call Ramp_Management(t,TimeRamp)
    
    W_WSC = 0._rp
    W_SL0 = 0._rp
    do j=Ind_sl0,Ind_sl1
        M = Mesh%Tnoeud(j)%Pnoeud
        GPhi0 = TimeRamp*Ecoulement%GPhi(:,j)%incident
        GPhip = Ecoulement%GPhi(:,j)%perturbation
        Etap = Ecoulement%Eta(j)%perturbation
        DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
        !DPhiDt = Ecoulement%DPhiDt(j)%perturbation
        DPhiDt = DPhiDt - dot_product(Mesh%Tnoeud(j)%Velocity,GPhip) ! Lagrange to Euler
        DPhi0Dn = TimeRamp*Ecoulement%DPhiDn(j)%incident
        DPhipDn = Ecoulement%DPhiDn(j)%perturbation
        call CD2Phi0DzDt(M,t,D2Phi0DzDt)
        DGPhi0Dz = TimeRamp*Ecoulement%DGPhiDz(:,j)%incident
        Dn = 0._RP
        do k=1,Mesh%Tnoeud(j)%Nfacette
            Dn = Dn + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*Mesh%Tnoeud(j)%Angle(k)
            !Dn = Dn + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*inv3
        end do
!        
        W_WSC = W_WSC + Etap*(TimeRamp*D2Phi0DzDt + dot_product(DGPhi0Dz(:),GPhi0) + g)*DPhi0Dn*Dn
        W_WSC = W_WSC - 0.5_rp*dot_product(GPhip,GPhip)*DPhi0Dn*Dn
        W_SL0 = W_SL0 - DPhiDt*DPhipDn*Dn
!
    enddo
    W_WSC = ro*W_WSC
    W_SL0 = ro*W_SL0

    Ind_body0 = Mesh%Body(Int_Body)%IndBody(1)
    Ind_body1 = Mesh%Body(Int_Body)%IndBody(3)
    Ind_sf0   = Mesh%Body(Int_Body)%IndBody(2)
    Ind_sf1   = Mesh%Body(Int_Body)%IndBody(4)

    W_Sigma = 0._rp
    FHydro = 0._rp
    FHydro0 = 0._rp
    FHydro1 = 0._rp
    FHydro2 = 0._rp
    do j=Ind_body0,Ind_Body1
        M = Mesh%Tnoeud(j)%Pnoeud
        GPhi0 = TimeRamp*Ecoulement%GPhi(:,j)%incident
        GPhip = Ecoulement%GPhi(:,j)%perturbation
        DPhiDt = TimeRamp*Ecoulement%DPhiDt(j)%incident + Ecoulement%DPhiDt(j)%perturbation
        !DPhiDt = DPhiDt - dot_product(Mesh%Tnoeud(j)%Velocity,GPhip)
        DPhi0Dn = TimeRamp*Ecoulement%DPhiDn(j)%incident
        GPhi2 = dot_product(GPhi0+GPhip,GPhi0+GPhip)
        Dn = 0._RP
        do k=1,Mesh%Tnoeud(j)%Nfacette
            Dn = Dn + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*Mesh%Tnoeud(j)%Angle(k)
        end do
        DS = Mesh%Tnoeud(j)%Normale*Dn

        FHydro  = FHydro + (DPhiDt + 0.5_rp*GPhi2 + g*M(3))*DS
        FHydro0 = FHydro0 + g*M(3) * DS
        FHydro1 = FHydro1 + DPhiDt * DS
        FHydro2 = FHydro2 + 0.5_rp*GPhi2*DS
    enddo

    FHydro = ro*FHydro
    FHydro0 = ro*FHydro0
    FHydro1 = ro*FHydro1
    FHydro2 = ro*FHydro2
    W_Sigma = ro*W_Sigma

    if(Symmetry)then
        FHydro = 2._rp*FHydro
        FHydro0 = 2._rp*FHydro0
        FHydro1 = 2._rp*FHydro1
        FHydro2 = 2._rp*FHydro2
        FHydro(2) = 0._rp        
        W_Sigma = 2._rp*W_Sigma
        W_WSC = 2._rp*W_WSC
        W_SL0 = 2._rp*W_SL0
    endif

    Fr(1:3) = FHydro(1:3)+[0._rp, 0._rp, -g*Mesh%Body(Int_Body)%Mass]
    Puissance = W_WSC + W_SL0 + W_Sigma - dot_product(Fr,Mesh%Body(Int_Body)%VBody(1:3))
    !Puissance = W_Sigma
    !Puissance = - dot_product(FHydro(1:3),Mesh%Body(Int_Body)%VBody(1:3))

    Mesh%Body(Int_Body)%FBody(1:12) = [FHydro(1:3),FHydro0(1:3),FHydro1(1:3),FHydro2(1:3)]
        

end subroutine FPuissance_WeakScatterer

! ----------------------------------------------------------------------------
! FPuissance_PseudoLin : calcul des variation d'énergie dans le domaine avec 
!                        une approche pseudo linéaire
! ----------------------------------------------------------------------------
subroutine FPuissance_PseudoLin(Mesh, Ecoulement, Puissance)
    implicit none
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(inout) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
    real(rp),intent(inout) :: Puissance
!   local
    integer :: j, k
    integer :: ind_sl0, ind_sl1,ind_body0, ind_body1
    real(rp) :: WEC_SL, WPP_SL, W_SM
    real(rp) :: DPhiDt, DPhiDn, Dn
    real(rp) :: Eta, DEtaDt
    real(rp),dimension(3) :: M, DS, FHydro

    ind_sl0 = Mesh%FS%IndFS(1)
    ind_sl1 = Mesh%FS%IndFS(3)

    WEC_SL = 0._rp
    WPP_SL = 0._rp
    do j=ind_sl0,ind_sl1
        M = Mesh%Tnoeud(j)%Pnoeud   
        !DPhi0Dt = DPhiDtAiry(M,t)
        !DPhiDt = DPhi0Dt + Ecoulement%DPhiDt(j)%perturbation
        DPhiDt = Ecoulement%DPhiDt(j)%perturbation
        !call GPhiAiry(M,t,GPhi0)
        !DPhi0Dn = dot_product(GPhi0,Mesh%Tnoeud(j)%Normale)
        !DPhiDn = DPhi0Dn + Ecoulement%DPhiDn(j)%perturbation
        DPhiDn = Ecoulement%DPhiDn(j)%perturbation
        !Eta0 = EtaAiry(M,t)
        !Eta = Eta0 + Ecoulement%Eta(j)%perturbation
        Eta = Ecoulement%Eta(j)%perturbation
        !DEta0Dt = DEtaDtAiry(M,t)
        !DEtaDt = DEta0Dt + Ecoulement%DEtaDt(j)%perturbation
        DEtaDt = Ecoulement%DEtaDt(j)%perturbation
        Dn = 0._RP
        do k=1,Mesh%Tnoeud(j)%Nfacette
            Dn = Dn + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*Mesh%Tnoeud(j)%Angle(k)
        end do
        WEC_SL = WEC_SL - DPhiDt*DPhiDn*Dn  
        WPP_SL = WPP_SL + g*Eta*DEtaDt*Dn
    enddo

    ind_body0 = Mesh%Body(Int_Body)%IndBody(1)
    ind_body1 = Mesh%Body(Int_Body)%IndBody(3)

    W_SM = 0._rp
    FHydro = 0._rp
    do j=ind_body0,ind_body1
        M = Mesh%Tnoeud(j)%Pnoeud   
        !DPhi0Dt = DPhiDtAiry(M,t)
        !DPhiDt = DPhi0Dt + Ecoulement%DPhiDt(j)%perturbation
        DPhiDt = Ecoulement%DPhiDt(j)%perturbation
        !call GPhiAiry(M,t,GPhi0)
        !DPhi0Dn = dot_product(GPhi0,Mesh%Tnoeud(j)%Normale)
        !DPhiDn = DPhi0Dn + Ecoulement%DPhiDn(j)%perturbation
        DPhiDn = Ecoulement%DPhiDn(j)%perturbation
        Dn = 0._RP
        do k=1,Mesh%Tnoeud(j)%Nfacette
            Dn = Dn + Mesh%Tfacette(Mesh%Tnoeud(j)%Tfacette(k,1))%Aire*Mesh%Tnoeud(j)%Angle(k)
        end do
        W_SM = W_SM - DPhiDt*DPhiDn*Dn
        if(sqrt(M(1)*M(1)+M(2)*M(2)).lt.Ldom(5)-LAbs)then
            DS = Mesh%Tnoeud(j)%Normale*Dn
            FHydro = FHydro + DPhiDt*DS
        endif
    enddo

    FHydro = ro*FHydro
    Puissance = ro*(WEC_SL+W_SM)

    if(symmetry)then
        FHydro = 2._rp*FHydro 
        FHydro(2) = 0._rp
        Puissance = 2._rp*Puissance
    endif

    Mesh%Body(Int_Body)%FBody(1:12) = [FHydro(1:3),0._rp,0._rp,0._rp,FHydro(1:3),0._rp,0._rp,0._rp]        

end subroutine FPuissance_PseudoLin

subroutine mass_conservation(Mesh,Ecoulement,t,ierror)
    implicit none
    !f2py integer*1, dimension(1000) :: Mesh
    type(TMaillage),intent(in) :: Mesh
    !f2py integer*1, dimension(1000) :: Ecoulement
    type(TEcoulement),intent(in) :: Ecoulement
    real(rp),intent(in) :: t
    integer,intent(inout) :: ierror
!   local
    integer :: j, k
    integer,dimension(3) :: TNoeud
    real(rp) :: h_moy, mass, eta
    real(rp),parameter :: inv3 = 0.33333333333333_rp        

    ierror = 0

    mass = 0._rp
    do j=Mesh%FS%IndFS(2),Mesh%FS%IndFS(4)
        h_moy = Ldom(3)
        Tnoeud = Mesh%Tfacette(j)%Tnoeud
        do k=1,3
            eta = Ecoulement%Eta(Tnoeud(k))%perturbation + Ecoulement%Eta(Tnoeud(k))%incident
            h_moy = h_moy + eta*inv3
        enddo
        mass = mass + h_moy*Mesh%Tfacette(j)%Aire
    enddo
    
    write(iomass,'(f8.4,f16.8)') t,mass    

end subroutine mass_conservation

end module EnergyVolume
