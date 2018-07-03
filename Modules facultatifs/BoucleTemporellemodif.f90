module BoucleTemp
use Constantes
use Parameters
use Structuresdonnees
use FonctionsCommunes
use GenMaillage
use GeomMesh
use Incident
use BVP
use Spline
use PrePlot
use GeomStruct
use SolvNum
use MeshModule
use BodyMotion
use EnergyVolume
use rbf_interp
use GeomGen
!use Energy
implicit none

contains

subroutine BoucleTemporelle(nt, t, MeshT,ierror,fgeom,fdomaine,nface,rep0,dx,&
&  mesh,nb_point,nb_arete,nb_tri)

implicit none

integer, intent(in) :: nt
real(rp), dimension(nt), intent(in) :: t
type(TMaillage),intent(inout) :: MeshT
integer, optional :: ierror
type(type_geom),intent(inout), optional :: fgeom
type(type_geom),intent(inout), optional :: fdomaine
integer,intent(in), optional :: nface
type(repere3d),intent(in), optional :: rep0
real(rp),intent(in), optional :: dx
type(MGrid),intent(inout), optional :: mesh
integer,intent(inout), optional :: nb_point,nb_arete,nb_tri

! Locales
character(len=50) :: filetemp, fileout, filemaill, num
integer :: j, k, r, s, jt, ios, Int_sl0, Int_sl1, Int_body0, Int_body1
integer :: Nnoeud
integer, dimension(2,2) :: borne
real(rp) :: h, time_begin, tmoy, rCI 
real(rp) :: Vcorps, ti, h2, inv_R, EtB, EtSL0, EtSL, Ethm1
real(rp), dimension(MeshT%NBody, 6) :: VBody
real(rp), dimension(3) :: Dn, M, GDPhiDt, Vg, Acc
real(rp), dimension(6,4) :: RKABody, RKVBody
real(rp), dimension(12,4) :: RKFBody
real(rp) :: Puissance(4), P_WSC, P_Bound
real(rp), dimension(:,:,:),allocatable :: RKVel
real(rp), dimension(nt) :: time_end, FpreNBVP, FpreFBVP, FpreNPP, FpreFPP
real(rp), allocatable :: CD(:,:), CS(:,:), RK(:,:,:), RKDPhiDt(:,:)
real(rp), allocatable :: DPhiDtCorps(:), PhiCorps(:,:), GradPhiCorps(:,:,:), PhiCorpsTemp(:)
real(rp), allocatable :: FPression(:,:), FPression0(:,:),  FPression1(:,:), FPression2(:,:)
real(rp), allocatable :: taux0(:), area0(:), fshape0(:), TMetric(:,:)
type(TMaillage) :: MeshTemp
type(TEcoulement) :: EcoulementT,  EcoulementTemp
logical :: bgenerate

! Parametres
real(rp), parameter ::   denom = 1._RP/6._RP
logical,parameter :: RK_fige = .false. 

! --- Début
call PrePlotPropHoule
call PrePlotPression2('Pression_PP_'//filename)  
if(iwts)then
    open(unit=iots, file='Ts.dat')
    write(iots,fmt='(50a)') 'Title = "Ts"'
    write(iots,fmt='(50a)') 'VARIABLES = "t","Ts(1)","CTdSM"'
    close(iots)
endif
open(unit=22, file='ABody.dat')
write(22,fmt='(50a)') 'Title = "ABody"'
write(22,fmt='(50a)') 'VARIABLES = "t","ABody(1)","ABody(3)","ABody(5)"'
close(22)
open(unit=23, file='GBody.dat')
write(23,fmt='(50a)') 'Title = "GBody"'
write(23,fmt='(50a)') 'VARIABLES = "t","Gx","Gy","Gz"'
close(23)

filemaill = 'mesh_'//filename
ierror = 0

  allocate(CS(MeshT.Nnoeud,MeshT.Nnoeud), CD(MeshT.Nnoeud,MeshT.Nnoeud))
  allocate(RK(MeshT%Nnoeud,4,2))
  allocate(PhiCorps(MeshT%Nnoeud,nt), GradPhiCorps(3,MeshT%Nnoeud,nt))
  allocate(RKVel(3,MeshT%Nnoeud,4))
  allocate(PhiCorpsTemp(MeshT%NNoeud))
  allocate(RKDPhiDt(MeshT%Nnoeud,4))
 ! 
  !allocate(CS(PointMax,PointMax), CD(PointMax,PointMax))
  !allocate(RK(PointMax,4,2))
  !allocate(PhiCorps(PointMax,nt), GradPhiCorps(3,PointMax,nt))
  !allocate(RKVel(3,PointMax,4))
  !allocate(PhiCorpsTemp(PointMax))
  !allocate(RKDPhiDt(PointMax,4))

! 
  call NewMaillage(MeshTemp,MeshT%Nnoeud)
!
  !call NewEcoulement(EcoulementT, MeshT%Nnoeud) 
  !call NewEcoulement(EcoulementTemp, MeshT%Nnoeud)
  call NewEcoulement(EcoulementT, PointMax)
  call NewEcoulement(EcoulementTemp, PointMax)
  call Initialisation(EcoulementT, MeshT, t(1))

  if(iwenergy)then
     EtB = 0._rp
     EtSL0 = 0._rp
     Ethm1 = 0._rp
     call Energy(MeshT,EcoulementT,t(1),h,EtSL0,EtSL,Ethm1,.false.,EtB)
     EtSL0 = EtSL
     Ethm1 = 0._rp
  endif

  PhiCorpsTemp(1:MeshT%Nnoeud) = EcoulementT%Phi(1:MeshT%Nnoeud)%perturbation

  !PhiT(1:MeshT%Nnoeud) = 0._rp

! -- Vitesse initial du flotteur et noeuds du maillage
  VBody(:,:) = 0._rp
  call MeshVel(MeshT,EcoulementT,t(1),dt,VBody)
  Nnoeud = MeshT%Nnoeud  
  call CopyMaillage(MeshTemp,MeshT)

! --------------------------------------------------------------------
!!   CONTROLE MAILLAGE
!    if(icheck)then
!!       Allocation dynamique
!        allocate(taux0(MeshTemp%Nfacette))
!        taux0(1:MeshTemp%Nfacette) = 1._rp
!        allocate(area0(MeshTemp%Nfacette))
!        allocate(fshape0(MeshTemp%Nfacette))
!        allocate(TMetric(MeshTemp%NNoeud,2))
!
!!       Qualite Maillage
!        call CheckMesh3D(MeshTemp,taux0,area0,fshape0,t(1),TMetric)
!        open(unit=iometric,file='metric.dat')
!        call PlotMetric(t0,MeshTemp,TMetric)
!
!!       Conservation Geometrie
!        open(unit=ioint,file="intersection_quality.dat")
!        open(unit=ioqual,file="immerged_surf_quality.dat")
!        call CheckGeom(MeshTemp,fgeom,t(1),ioint,ioqual)
!    endif
! -----------------------------------------------------------------------
!
! Calculation of the Influence Coefficients
  call cpu_time(time_begin)

  CD = 0._RP ; CS = 0._RP
  call CoeffInfl(MeshTemp, CD, CS)

  call cpu_time(time_end(1))
  print*, 'Temps CI :', time_end(1)-time_begin
  call cpu_time(time_begin)

  inv_R = 1._rp/Lgeom(2)

! Boucle en temps
  do jt=1,nt

    h=t(2)-t(1)
    h2 = 0.5_rp*h
    ti = t(jt)

! ## CC Forçage déplacement initial pour mouvement libre
    !if(ti.lt.9._rp)then
    !    free_body=.false.
    !else
    !    free_body=.true.
    !endif
! ##

    RK = 0._RP
    
    Nnoeud = MeshT%Nnoeud

    VBody(Int_Body,1:6) = MeshT%Body(Int_Body)%VBody(1:6)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       RK 1ere passe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Stockage des coeff RK pour les vitesses
    do j=1,Nnoeud
      RKVel(1:3,j,1) = MeshT%Tnoeud(j)%Velocity(1:3)
    enddo
    RKVBody(1:3,1) = MeshT%Body(Int_Body)%VBody(1:3)
    RKABody(1:3,1) = MeshTemp%Body(Int_Body)%ABody(1:3)

! -- Gestion du recalcul partiel des coeff d'influence
    rCI = -1._RP
    !if (RK_fige) then
    !    ! Modification du maillage
    !    call Remesh(Mesh, Mesh0, ti+h)
    !    rCI = -1._RP
    !end if
    !call CopyMaillage(Mesh0,Mesh)

! --- Résolution du Problème aux limites au temps t
    call solBVP( EcoulementT, MeshT, CD, CS, rCI)
    if(oplot)then
        Write( num, '( f0.4 )' ) ti
        fileout = 'RK1_'//trim(num)//filename
        call PlotBVP(fileout, MeshT, EcoulementT)
    endif

    if(iwenergy .and. jt.gt.1) call Energy(MeshT,EcoulementT,ti+h,h,EtSL0,EtSL,Ethm1,.false.,EtB)
    
!   Calcul des Dérivées Temporelles et Spatiales de Phi(t) et Eta(t)
    call Derive_2(ti, EcoulementT, MeshT)      
    if (oplot) call PlotEcoul(fileout, MeshT, EcoulementT)

!   Stockage des coeff RK pour Phi et Eta
    Int_sl0=MeshT%FS%IndFS(1)
    Int_sl1=MeshT%FS%IndFS(3)
    RK(Int_sl0:Int_sl1,1,1)=EcoulementT%DPhiDt(Int_sl0:Int_sl1)%perturbation
    RK(Int_sl0:Int_sl1,1,2)=EcoulementT%DEtaDt(Int_sl0:Int_sl1)%perturbation

!   Calcul des pentes dérivées pour la première passe, puis de Phi(t+h/2) et Eta(t+h/2) sur le MaillageT, à t
    call CopyEcoulement(EcoulementTemp, EcoulementT, MeshT%Nnoeud)

!   Calcul des efforts hydrodynamiques 
    if(free_body)then
        call FreeBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti)
    else
        call ForceBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti)
        call solBVP(EcoulementTemp, MeshTemp, CD, CS, ti, .true.)
    endif

!   Puissance
    call FPuissance(MeshT, EcoulementTemp, t(jt), Puissance(1))
    RKFBody(1:12,1) = MeshT%Body(Int_Body)%FBody(1:12)
    Int_body0 = MeshT%Body(Int_Body)%IndBody(1)
    Int_Body1 = MeshT%Body(Int_Body)%IndBody(3)
    do j=Int_Body0,Int_Body1
        RKDPhiDt(j,1) = EcoulementTemp%DPhiDt(j)%Perturbation !+ dot_product(MeshT%Tnoeud(j)%Velocity,EcoulementTemp%GPhi(1:3,j)%perturbation)
    enddo

!   Mise a jour Ecoulement t+h/2
    EcoulementTemp%Phi(Int_sl0:Int_sl1)%perturbation = EcoulementT%Phi(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,1,1)*h2
    EcoulementTemp%Eta(Int_sl0:Int_sl1)%perturbation = EcoulementT%Eta(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,1,2)*h2
    
!   Remaillage
    if (not(RK_fige) .and. moving_mesh) call Remesh(MeshTemp, MeshT, ti+h2,h2,fgeom%repere%e3)

!   Vitesse de déplacement du maillage
    call MeshVel(MeshTemp,EcoulementTemp,ti+h2,h2,VBody(:,:))
    do j=1,Nnoeud
        RKVel(1:3,j,2) = MeshTemp%Tnoeud(j)%Velocity(1:3)
    enddo
    RKVBody(1:3,2) = MeshTemp%Body(Int_Body)%VBody(1:3)
    RKABody(1:3,2) = MeshTemp%Body(Int_Body)%ABody(1:3)   

! --- Initialisation next step
    call CCL(ti+h2, MeshTemp, EcoulementTemp)
    call BodyCondition(ti+h2, MeshTemp, EcoulementTemp)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !       RK 2eme passe
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rCI = -1._RP
    if (RK_fige)  rCI = 1._RP

! --- Résolution du Problème aux limites au temps t+h/2
    call solBVP( EcoulementTemp, MeshTemp, CD, CS, rCI)
    if(oplot)then
        fileout = 'RK2_'//trim(num)//filename
        call PlotBVP(fileout, MeshTemp, EcoulementTemp)
    endif

! --- Calcul des dérivées Temporelles et Spatiales de Phi(t) et Eta(t)
    call Derive_2(ti+h2, EcoulementTemp, MeshTemp)
    if (oplot) call PlotEcoul(fileout, MeshTemp, EcoulementTemp)
    RK(Int_sl0:Int_sl1,2,1) = EcoulementTemp%DPhiDt(Int_sl0:Int_sl1)%perturbation
    RK(Int_sl0:Int_sl1,2,2) = EcoulementTemp%DEtaDt(Int_sl0:Int_sl1)%perturbation

! -- Force hydrodynamique
    if(free_body)then
        call FreeBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h2)
    else
        call ForceBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h2)
        call solBVP(EcoulementTemp, MeshTemp, CD, CS, ti+h2, .true.)
    endif

!   Puissance  
    call FPuissance(MeshTemp, EcoulementTemp, ti+h2, Puissance(2))
    RKFBody(1:12,2) = MeshTemp%Body(Int_Body)%FBody(1:12)
    do j=Int_Body0,Int_Body1
        RKDPhiDt(j,2) = EcoulementTemp%DPhiDt(j)%perturbation + dot_product(MeshTemp%Tnoeud(j)%Velocity,EcoulementTemp%GPhi(1:3,j)%perturbation)
    enddo

!   Mise a jour Ecoulement a t+h/2
    EcoulementTemp%Phi(Int_sl0:Int_sl1)%perturbation = EcoulementT%Phi(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,2,1)*h2
    EcoulementTemp%Eta(Int_sl0:Int_sl1)%perturbation = EcoulementT%Eta(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,2,2)*h2

!   Body Velocity
    if (not(RK_fige) .and. moving_mesh) call Remesh(MeshTemp, MeshT, ti+h2, h2, fgeom%repere%e3)
    call MeshVel(MeshTemp,EcoulementTemp,ti+h2,h2,VBody(:,:))
    do j=1,Nnoeud
        RKVel(1:3,j,3) = MeshTemp%Tnoeud(j)%Velocity(1:3)
    enddo
    RKVBody(1:3,3) = MeshTemp%Body(Int_Body)%VBody(1:3)
    RKABody(1:3,3) = MeshTemp%Body(Int_Body)%ABody(1:3)

! --- Initialisation next step
    call BodyCondition(ti+h2, MeshTemp, EcoulementTemp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       RK 3eme passe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rCI = -1._RP
    if (not(RK_fige).and.is_body.and.Free_Body)  rCI = -1._RP

! --- Résolution du Problème aux Limites au temps t+h/2
    call solBVP( EcoulementTemp, MeshTemp, CD, CS, rCI)
    if(oplot)then
        fileout = 'RK3_'//trim(num)//filename
        call PlotBVP(fileout, MeshTemp, EcoulementTemp)   
    endif

! --- Calcul des Dérivées Temporelles et Spatiales de Phi(t) et Eta(t)
    call Derive_2(ti+h2, EcoulementTemp, MeshTemp)
    if (oplot) call PlotEcoul(fileout, MeshTemp, EcoulementTemp)
    RK(Int_sl0:Int_sl1,3,1)=EcoulementTemp%DPhiDt(Int_sl0:Int_sl1)%perturbation
    RK(Int_sl0:Int_sl1,3,2)=EcoulementTemp%DEtaDt(Int_sl0:Int_sl1)%perturbation

!   Calcul des efforts hydrodynamiques 
    if(free_body)then
        call FreeBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h2)
    else
        call ForceBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h2)
        call solBVP(EcoulementTemp, MeshTemp, CD, CS, ti+h2, .true.)
    endif

!   Puissance
    call FPuissance(MeshTemp, EcoulementTemp, ti+h2, Puissance(3))
    RKFBody(1:12,3) = MeshTemp%Body(Int_Body)%FBody(1:12)
    do j=Int_Body0,Int_Body1
        RKDPhiDt(j,3) = EcoulementTemp%DPhiDt(j)%perturbation + dot_product(MeshTemp%Tnoeud(j)%Velocity,EcoulementTemp%GPhi(1:3,j)%perturbation)
    enddo

!   Mise a jour Ecoulement a t+h
    EcoulementTemp%Phi(Int_sl0:Int_sl1)%perturbation = EcoulementT%Phi(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,3,1)*h
    EcoulementTemp%Eta(Int_sl0:Int_sl1)%perturbation = EcoulementT%Eta(Int_sl0:Int_sl1)%Perturbation + RK(Int_sl0:Int_sl1,3,2)*h

! --- Body Velocity
    if (not(RK_fige) .and. moving_mesh) call Remesh(MeshTemp, MeshT, ti+h, h, fgeom%repere%e3)
    call MeshVel(MeshTemp,EcoulementTemp,ti+h,h,VBody(:,:))
    do j=1,Nnoeud
        RKVel(1:3,j,4) = MeshTemp%Tnoeud(j)%Velocity(1:3)
    enddo
    RKVBody(1:3,4) = MeshTemp%Body(Int_Body)%VBody(1:3)
    RKABody(1:3,4) = MeshTemp%Body(Int_Body)%ABody(1:3)
    
! --- Initialisation next step
    call CCL(ti+h, MeshTemp, EcoulementTemp)
    call BodyCondition(ti+h, MeshTemp, EcoulementTemp)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       RK 4eme passe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rCI = -1._RP
    if (RK_fige)  rCI = 1._RP

! --- Résolution du Problème aux limites au temps t+h
    call solBVP(EcoulementTemp, MeshTemp, CD, CS, rCI)
    if(oplot)then
        fileout = 'RK4_'//trim(num)//filename
        call PlotBVP(fileout, MeshTemp, EcoulementTemp)
    endif
   
! --- Calcul des dérivées Temporelles et Spatiales de Phi(t) et Eta(t)
    call Derive_2(ti+h, EcoulementTemp, MeshTemp)
    if (oplot) call PlotEcoul(fileout, MeshTemp, EcoulementTemp)
    RK(Int_sl0:Int_sl1,4,1)=EcoulementTemp%DPhiDt(Int_sl0:Int_sl1)%perturbation
    RK(Int_sl0:Int_sl1,4,2)=EcoulementTemp%DEtaDt(Int_sl0:Int_sl1)%perturbation

!   Calcul des efforts hydrodynamiques 
    if(free_body)then
        call FreeBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h)
    else
        call ForceBodyMotion(MeshTemp, EcoulementTemp, CD, CS, ti+h)
        call solBVP(EcoulementTemp, MeshTemp, CD, CS, ti+h, .true.)
    endif

!   Puissance
    call FPuissance(MeshTemp, EcoulementTemp, ti+h, Puissance(4))
    RKFBody(1:12,4) = MeshTemp%Body(Int_Body)%FBody(1:12)
    do j=Int_Body0,Int_Body1
        RKDPhiDt(j,4) = EcoulementTemp%DPhiDt(j)%perturbation + dot_product(MeshTemp%Tnoeud(j)%Velocity,EcoulementTemp%GPhi(1:3,j)%perturbation)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Actualisation fin RK4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Calcul de Phi et Eta au temps t+h        
    EcoulementTemp%Phi(Int_sl0:Int_sl1)%perturbation = EcoulementT%Phi(Int_sl0:Int_sl1)%Perturbation&
                                                     + (RK(Int_sl0:Int_sl1,1,1)+2._RP*RK(Int_sl0:Int_sl1,2,1)+2._RP*RK(Int_sl0:Int_sl1,3,1)+RK(Int_sl0:Int_sl1,4,1))*h*denom
    EcoulementTemp%Eta(Int_sl0:Int_sl1)%perturbation = EcoulementT%Eta(Int_sl0:Int_sl1)%Perturbation&
                                                     + (RK(Int_sl0:Int_sl1,1,2)+2._RP*RK(Int_sl0:Int_sl1,2,2)+2._RP*RK(Int_sl0:Int_sl1,3,2)+RK(Int_sl0:Int_sl1,4,2))*h*denom
    EcoulementT.DPhiDt(Int_sl0:Int_sl1).perturbation = RK(Int_sl0:Int_sl1,1,1)
    EcoulementT.DEtaDt(Int_sl0:Int_sl1).perturbation = RK(Int_sl0:Int_sl1,1,2)
    EcoulementTemp.DPhiDt(Int_sl0:Int_sl1).perturbation = (RK(Int_sl0:Int_sl1,1,1)+2._RP*RK(Int_sl0:Int_sl1,2,1)+2._RP*RK(Int_sl0:Int_sl1,3,1)+RK(Int_sl0:Int_sl1,4,1))*denom
    EcoulementTemp.DEtaDt(Int_sl0:Int_sl1).perturbation = (RK(Int_sl0:Int_sl1,1,2)+2._RP*RK(Int_sl0:Int_sl1,2,2)+2._RP*RK(Int_sl0:Int_sl1,3,2)+RK(Int_sl0:Int_sl1,4,2))*denom

! --- Mise a jour de DPhiDt sur le corps pour verification formulation green sur DPhiDt
    EcoulementTemp%DPhiDt(Int_Body0:Int_Body1)%Perturbation = (RKDPhiDt(Int_Body0:Int_Body1,1) + 2._rp*RKDPhiDt(Int_Body0:Int_Body1,2) + 2._rp*RKDPhiDt(Int_Body0:Int_Body1,3) +&
&                                                              RKDPhiDt(Int_Body0:Int_Body1,4))*denom


! --- Body Velocity
    !MeshTemp%Body(Int_Body)%ABody = RKABody(:,1) ! a modifier car champ a t non a t+dt
    MeshTemp%Body(Int_Body)%ABody = (RKABody(:,1) + 2._RP*RKABody(:,2) + 2._RP*RKABody(:,3) + RKABody(:,4))*denom
    if(free_body)then
      !MeshTemp%Body(Int_Body)%VBody = MeshT%Body(Int_Body)%VBody + (RKABody(:,1) + 2._RP*RKABody(:,2) + 2._RP*RKABody(:,3) + RKABody(:,4))*denom*h
      MeshTemp%Body(Int_Body)%VBody = (RKVBody(:,1) + 2._rp*RKVBody(:,2) + 2._rp*RKVbody(:,3) + RKVBody(:,4))*denom
    else
      !MeshTemp%Body(Int_Body)%VBody(1:3) = RKVBody(1:3,4)
      MeshTemp%Body(Int_Body)%VBody = (RKVBody(:,1) + 2._rp*RKVBody(:,2) + 2._rp*RKVbody(:,3) + RKVBody(:,4))*denom
    endif
    do j=1,Nnoeud
        MeshTemp%Tnoeud(j)%Velocity(1:3) = (RKVel(1:3,j,1) + 2._rp*RKVel(1:3,j,2) + 2._rp*RKVel(1:3,j,3) + RKVel(1:3,j,4))*denom
        !MeshT%Tnoeud(j)%Velocity(1:3) = MeshTemp%Tnoeud(j)%Velocity
    enddo
    MeshT%Body(Int_Body)%VBody = MeshTemp%Body(Int_Body)%VBody

! --- Force hydro
    !MeshTemp%Body(Int_Body)%FBody = (RKFBody(:,1) + 2._rp*RKFBody(:,2) + 2._rp*RKFBody(:,3) + RKFBody(:,4))*denom
    MeshTemp%Body(Int_Body)%FBody = RKFBody(:,1)
    !EtB = EtB - (Puissance(1) + 2._rp*Puissance(2) + 2._rp*Puissance(3) + Puissance(4))*h*denom
    EtB = EtB + (Puissance(1) + 2._rp*Puissance(2) + 2._rp*Puissance(3) + Puissance(4))*h*denom
    if(iwpress) write(iopress,'(13E)') ti, MeshTemp%Body(Int_Body)%FBody(1:12)

    if (not(RK_fige) .and. moving_mesh) then
        call Remesh(MeshTemp, MeshT, ti+h,h,fgeom%repere%e3,.false.) 
        call CheckMesh3D_2(MeshTemp,ti+h,bgenerate)
        if(bgenerate)then
            call update_geom(rep0,fgeom,[0._rp,0._rp,0._rp],MeshTemp%Body(Int_Body)%GBody(1:3))
            Position(1:3,1) = MeshTemp%Body(Int_Body)%GBody(1:3)
            call compute_mesh_fgeom(MeshTemp,ti+h,fgeom,Ldom(3),nface,rep0,dx2,ierror)
            if(ierror/=0) goto 9999
            Nnoeud = MeshTemp%Nnoeud

            if(allocated(CS)) deallocate(CS)
            if(allocated(CD)) deallocate(CD)
            if(allocated(RK)) deallocate(RK)
            if(allocated(PhiCorps)) deallocate(PhiCorps)
            if(allocated(RKVel)) deallocate(RKVel)
            if(allocated(GradPhiCorps)) deallocate(GradPhiCorps)
            if(allocated(PhiCorpsTemp)) deallocate(PhiCorpsTemp)
            if(allocated(RKDPhiDt)) deallocate(RKDPhiDt)

            allocate(CS(Nnoeud,Nnoeud), CD(Nnoeud,Nnoeud))
            allocate(RK(Nnoeud,4,2))
            allocate(PhiCorps(Nnoeud,nt), GradPhiCorps(3,Nnoeud,nt))
            allocate(RKVel(3,Nnoeud,4))
            allocate(PhiCorpsTemp(Nnoeud))
            allocate(RKDPhiDt(Nnoeud,4))

        endif
    endif   
    !if (is_body) call BodyVel(MeshTemp,ti+h,h,VBody)

! --- Initialisation next step
    call CCL(ti+h, MeshTemp, EcoulementTemp)
    call BodyCondition(ti+h, MeshTemp, EcoulementTemp)
! CC
    call solBVP(EcoulementTemp, MeshTemp, CD, CS)
    if (is_body) call MeshVel(MeshTemp,EcoulementTemp,ti+h,h,VBody(:,:))

    !call FHydro_ordre1(MeshT,EcoulementT,ti,MeshTemp,EcoulementTemp,ti+h)
    !Vg = 0.5_rp*(MeshT%Body(Int_Body)%VBody(1:3) + MeshTemp%Body(Int_Body)%VBody(1:3))
    !Vg = MeshT%Body(Int_Body)%VBody(1:3)
    !EtB = EtB - dot_product(MeshT%Body(Int_Body)%FBody(1:3),Vg)*h
    !if(iwpress) write(iopress,'(13E)') ti, MeshT%Body(Int_Body)%FBody(1:12)

! CC
    !if (is_body .and. not(free_body)) then
        !Int_Body0 = MeshTemp%Body(Int_Body)%IndBody(1)
        !Int_Body1 = MeshTemp%Body(Int_Body)%IndBody(3)
        !PhiCorps(Int_Body0:Int_Body1,jt) = EcoulementT.Phi(Int_Body0:Int_Body1).perturbation
        !do j=Int_Body0,Int_Body1
        !    GradPhiCorps(:,j,jt) = EcoulementT%GPhi(:,j)%perturbation
        !end do
    !else
        !call FHydro(MeshT, EcoulementT, EcoulementTemp, t(jt),dt)
    !    call Energy(Mesh, EcoulementT, t(jt))
    !end if

! --- Verification mise a jour Phi sur le corps
    !do j=Int_Body0,Int_Body1
    !    PhiCorpsTemp(j) = PhiCorpsTemp(j) + h*EcoulementTemp%DPhiDt(j)%perturbation
    !enddo
    !call PlotPhiCorps(MeshTemp,EcoulementTemp,PhiCorpsTemp,ti+h)

    call CopyEcoulement(EcoulementT,EcoulementTemp,MeshT%Nnoeud)
    call CopyMaillage(MeshT,MeshTemp)

    if(iwenergy) then
      call Energy(MeshT,EcoulementT,ti+h,h,EtSL0,EtSL,Ethm1,.true.,EtB)
      if(ti+h.gt.T2-Epsilon .and. abs(ti+h-T2).lt.dt)then
        EtB = EtSL
      endif
    endif

 ! --- Lissage des solutions EcoulementT
    !if(.not.Htype.eq.0 .and. mod(jt+1,nliss).eq.0)then
     if(mod(jt+1,nliss).eq.0)then
      call lissage(MeshT,EcoulementT,ierror)
      if(ierror/=0)then
        ierror = 200
        goto 9999
      endif
    endif

    call PlotPropHoule2(ti+h, MeshT, EcoulementT, h)

    if(iwmesh) call PlotMaill(filemaill,MeshT,ti+h,iomesh)

    if(iwpress3D) call PlotPression3D(ti+h,MeshT, EcoulementT)

    if(iwmbody) call PlotMBody(ti+h,MeshTemp)

    call IniEcoulement(EcoulementTemp, MeshTemp%Nnoeud, 0._RP)
! ----------------------------------------------------------------
!!   Controle Maillage
!    if(icheck)then
!        call CheckGeom(MeshT,fgeom,ti+h,ioint,ioqual)
!        call CheckMesh3D(MeshT,taux0,area0,fshape0,ti+h,TMetric)
!        call PlotMetric(ti+h,MeshT,TMetric)
!    endif 
! -----------------------------------------------------------------
!
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
print*, 'Temps total du calcul : ', time_end(nt)-time_begin, ' secondes'

call close_output

! Post process pression
!if (is_body .and. (.not.free_body)) then
    !allocate(FpressionF(3,MaillageT.Nfacette),FpressionF0(3,MaillageT.Nfacette),FpressionF1(3,MaillageT.Nfacette),FpressionF2(3,MaillageT.Nfacette))
!    allocate(Fpression(3,MeshT.Nnoeud),Fpression0(3,MeshT.Nnoeud),Fpression1(3,MeshT.Nnoeud),Fpression2(3,MeshT.Nnoeud))
!    allocate(DPhiDtCorps(MeshT%Body(Int_Body)%IndBody(1):MeshT%Body(Int_Body)%IndBody(3)))
!    call PrePlotPression2('Pression_PP_'//filename)    
!    do jt=1,nt-1
!        ! Somme sur les noeuds
!        !Vcorps = -Acorps*wcorps*sin(wcorps*t(jt))
!        do j=MeshT%Body(Int_Body)%IndBody(1),MeshT%Body(Int_Body)%IndBody(3)
!            ! Calcul DphiDt
!            if (jt.eq.1) then
!                DphiDtCorps(j) = (PhiCorps(j,jt+1)-PhiCorps(j,jt))/(t(jt+1)-t(jt))
!            else
!                DphiDtCorps(j) = (PhiCorps(j,jt+1)-PhiCorps(j,jt-1))/(t(jt+1)-t(jt-1))
!           end if
!            Dn = 0._RP
!            do k=1,MeshT.Tnoeud(j).Nfacette
!                Dn = Dn + MeshT.Tfacette(MeshT.Tnoeud(j).Tfacette(k,1)).Aire*MeshT.Tnoeud(j).Angle(k)*MeshT.Tnoeud(j).Normale
!            end do
!            Fpression0(:,j) = ro*g*MeshT.Tnoeud(j).Pnoeud(3)*Dn
!            Fpression1(:,j) = ro*(DphiDtCorps(j) - dot_product(MeshT%Tnoeud(j)%Velocity,GradPhiCorps(:,j,jt)))*Dn
!            Fpression2(:,j) = ro*0.5_RP*dot_product(GradPhiCorps(:,j,jt),GradPhiCorps(:,j,jt))*Dn
!            Fpression(:,j) = Fpression0(:,j) + Fpression1(:,j) + Fpression2(:,j)
!            FpreNPP(jt) = FpreNPP(jt) + Fpression1(3,j) + Fpression2(3,j)
!        end do
!        call plotPression2(t(jt), MeshT, Fpression, Fpression0, Fpression1, Fpression2)
!    end do
!    close(58)
!    deallocate(Fpression, Fpression0, Fpression1, Fpression2, DPhiDtCorps)
!end if

filetemp='Temps_'//filename
open(unit=17,file=filetemp, iostat=ios)
if (ios/=0) stop "Erreur à l'ouverture du fichier Temps"
write(17,fmt='(a,e,a)')'Temps total du calcul : ', time_end(nt)-time_begin, ' secondes'
close(17)

9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
99 format("error #",i3)

    if(ierror/=0)then
      Write( num, '( f0.4 )' ) ti+h
      !filemaill = 'Maillage_'//trim(num)//'.dat'
      filemaill = 'MeshTemp_last.dat'
      call PlotMaill(filemaill, MeshTemp)
    endif

if(allocated(CD))  deallocate(CD)
if(allocated(CS))  deallocate(CS)
if(allocated(RK))  deallocate(RK)
if(allocated(GradPhiCorps))  deallocate(GradPhiCorps)
if(allocated(PhiCorps))      deallocate(PhiCorps)
if(allocated(RKVel)) deallocate(RKVel)
if(allocated(RKDPhiDt)) deallocate(RKDPhiDt)
if(allocated(PhiCorpsTemp)) deallocate(PhiCorpsTemp)
if(allocated(taux0)) deallocate(taux0)
if(allocated(area0)) deallocate(area0)
if(allocated(fshape0)) deallocate(fshape0)
if(allocated(TMetric)) deallocate(TMetric)

if(icheck)then
    close(iometric)
    close(ioint)
    close(ioqual)
endif
!
end subroutine BoucleTemporelle

!-------------------------------------------------------------
!
!	Initialisation de l'écoulement pour l'intégration en temps
!
!--------------------------------------------------------------
!
subroutine BodyCondition(t, Mesh, Ecoulement)
! Parameters
real(rp) :: t
Type(Tmaillage) :: Mesh
Type(TEcoulement) :: Ecoulement
! Locals
integer :: j, k
real(rp) :: TimeRamp, Vn
real(rp),dimension(3) :: OmegaMG
real(rp), dimension(3,3) :: Puvw
! Begin
!if(is_body)then
    !if (Free_Body) then
    !    do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%body(Int_Body)%IndBody(3)
    !        Mesh%Tnoeud(j)%Velocity = Mesh%Body(nc)%VBody(1:3) + vect_product(Mesh%Body(nc)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-Mesh%Body(nc)%GBody(1:3))
    !    end do
    !end if
    if (t.ge.T2) then
        TimeRamp = 1._RP
    else
        TimeRamp = CRampe(t,T1,T2)
    end if
    do j=Mesh%Body(Int_Body)%IndBody(1),Mesh%Body(Int_Body)%IndBody(3)
        ! Vitesse Normale
        !OmegaMG = vect_product(Mesh%Body(Int_Body)%VBody(4:6),Mesh%Tnoeud(j)%Pnoeud(1:3)-Mesh%Body(Int_Body)%GBody(1:3))
        !Vn = dot_product(Mesh%Body(Int_Body)%VBody(1:3) + OmegaMG,Mesh%Tnoeud(j)%Normale)
        !if(norm(Mesh%Tnoeud(j)%Pnoeud).gt.Ldom(5)-LAbs .and. abs(Mesh%Tnoeud(j)%Normale(3)).lt.Epsilon2)then
        if(norm(Mesh%Tnoeud(j)%Pnoeud).gt.Ldom(5)-LAbs)then
        !if(norm(Mesh%Tnoeud(j)%Pnoeud).gt.Ldom(5)-Epsilon2)then
            Ecoulement%DPhiDn(j)%perturbation = 0._rp
        else
            Ecoulement%DPhiDn(j)%perturbation = TimeRamp*intvelocity(t,Mesh%Tnoeud(j),Mesh%Body(Int_Body)%VBody,Mesh%Body(Int_Body)%GBody)&
&                                             - TimeRamp*Ecoulement%DPhiDn(j)%incident
        endif
        !if(Mesh%Tnoeud(j)%Pnoeud(1) .lt. 0. .and. abs(Mesh%Tnoeud(j)%Normale(2)).lt.Epsilon)then
        !    Ecoulement%DPhiDn(j)%perturbation = intvelocity(t,Mesh%Tnoeud(j),Mesh%Body(Int_Body)%VBody,Mesh%Body(Int_Body)%GBody)&
!&                                             - TimeRamp*Ecoulement%DPhiDn(j)%incident
        !else
        !    Ecoulement%DPhiDn(j)%perturbation = 0._rp
        !endif
        !Ecoulement%DPhiDn(j)%perturbation =  dot_product(Mesh%Tnoeud(j)%Velocity,Mesh%Tnoeud(j)%Normale) - Ecoulement%DPhiDn(j)%incident
        ! Dérivées partielles de la vitesse normale (dans la base locale)
        Puvw(1:3,1:3) = Mesh%Tnoeud(j)%Plocal(1:3,1:3,1)
        Ecoulement%GsDPhiDn(1,j)%perturbation = -dot_product(Mesh%Body(Int_Body)%VBody(4:6),Puvw(1:3,2)) &
                                                - dot_product(Mesh%Tnoeud(j)%Velocity,Puvw(1:3,1))*Mesh%Tnoeud(j)%DLocal(3)
        Ecoulement%GsDPhiDn(2,j)%perturbation =  dot_product(Mesh%Body(Int_Body)%VBody(4:6),Puvw(1:3,1)) &
                                                - dot_product(Mesh%Tnoeud(j)%Velocity,Puvw(1:3,2))*Mesh%Tnoeud(j)%DLocal(4)
    end do
!endif
if (cuve_ferme) then
    do k=1,Mesh%NBody
        if (k.ne.Int_body) then
            do j=Mesh%Body(k)%IndBody(1),Mesh%Body(k)%IndBody(3)
                Ecoulement%DPhiDn(j)%perturbation = 0._RP
            end do
        end if
    end do
end if
! End
end subroutine BodyCondition

!-------------------------------------------------------------
!
!	Dérivations spatiales de Phi et Eta
!   Dérivations temporelles de phi et Eta grâce aux CSL
!
!-------------------------------------------------------------
subroutine Derive(t, Ecoulement, Mesh, matv)
! Paramètres
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
real(rp),dimension(:,:),intent(in), optional :: matv
! Variables Locales
integer :: j, k, r, s
real(rp) :: Rayon, RayonTot, mu, Eta, Phi, Cs, RDamping 
real(rp), dimension(3) :: M, OM, GEta0, GEta, GPhi, GPhi0, DGPhi0Dz
real(rp), dimension(3,3) :: Puvw, Pt
! Début
if (Htype.eq.0) then
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        Ecoulement.DPhiDt(j).perturbation = - g*Ecoulement%Eta(j)%Perturbation*Mesh.Tnoeud(j).Damping
        Ecoulement.DpPhiDt(j).perturbation = Ecoulement.DPhiDt(j).perturbation
        Ecoulement.DEtaDt(j).perturbation = -Ecoulement%DPhiDn(j)%Perturbation*Mesh.Tnoeud(j).Damping
    end do
else
    call CCL(t,Mesh, Ecoulement)
    call GradientFS(Ecoulement,Mesh,t)
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        DGPhi0Dz = Ecoulement.DGPhiDz(:,j).incident
		GPhi0 = Ecoulement.GPhi(:,j).incident
		GEta0 = Ecoulement.GEta(:,j).incident
		Gphi=Ecoulement.GPhi(:,j).perturbation
		GEta=Ecoulement.GEta(:,j).perturbation			
		Eta = Ecoulement.Eta(j).perturbation
		Phi = Ecoulement.Phi(j).perturbation
		!! Derivees Temporelles
        Ecoulement.DPhiDt(j).perturbation = (- Ecoulement.DPhiDt(j).incident - g*Ecoulement.Eta(j).incident - 0.5_RP*dot_product(Gphi0,Gphi0) &
		                                    -g*Eta - dot_product(Gphi0,Gphi) + Gphi0(3)*Gphi(3) + 0.5_RP*Eta*dot_product(DGPhi0Dz,Gphi0))*Mesh.Tnoeud(j).Damping
		Ecoulement.DpPhiDt(j).perturbation = (Ecoulement.DPhiDt(j).perturbation - GPhi0(3)*(GPhi0(3)+Gphi(3)))*Mesh.Tnoeud(j).Damping
		Ecoulement.DEtaDt(j).perturbation = (- Ecoulement.DEtaDt(j).incident + Gphi0(3) - dot_product(Gphi0,Geta0) + Gphi(3) - dot_product(Gphi0,GEta)&
		                                     - dot_product(Gphi,Geta0)-Eta*dot_product(DGPhi0Dz,GEta0))*Mesh.Tnoeud(j).Damping
    end do
end if
! Fin
end subroutine Derive

subroutine Derive_2(t, Ecoulement, Mesh)
! Paramètres
type(TMaillage) :: Mesh
type(TEcoulement) :: Ecoulement
real(rp) :: t
! Variables Locales
integer :: j, k, r, s
real(rp) :: Rayon, RayonTot, mu, Eta, Phi, Cs, RDamping 
real(rp), dimension(3) :: M, OM, GEta0, GEta, GPhi, GPhi0, DGPhi0Dz !, GDEta0Dt, GDEta0Dz, GDPhi0Dt, GDPhiDz, GDPhiDt
real(rp),dimension(3) :: Vel
real(rp), dimension(3,3) :: Puvw, Pt
real(rp) :: D2Phi0DzDt, damp1, damp2
integer,parameter :: tdamp = 2
!real(rp), dimension(:,:) :: matv ! A mettre dans le type Tnoeud
! Début

if (Htype.eq.0) then
    call GradientFS(Ecoulement,Mesh,t)
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        select case(tdamp)
	    case(1)
		    damp1 = Mesh%Tnoeud(j)%Damping
		    damp2 = 0._rp
	    case(2) 
		    damp1 = 1._rp
		    damp2 = Mesh%Tnoeud(j)%Damping  
        end select
        Gphi=Ecoulement.GPhi(:,j).perturbation
        GEta=Ecoulement.GEta(:,j).perturbation
        Vel = Mesh%Tnoeud(j)%Velocity
        Ecoulement.DPhiDt(j).perturbation = - g*Ecoulement%Eta(j)%Perturbation + dot_product(Vel,GPhi)
        !Ecoulement.DPhiDt(j).perturbation = - g*Ecoulement%Eta(j)%Perturbation - 0.5_rp*dot_product(GPhi,GPhi) + dot_product(Vel,GPhi)
		Ecoulement%DPhiDt(j)%perturbation = damp1*Ecoulement%DPhiDt(j)%perturbation - damp2*Ecoulement%Phi(j)%perturbation
        Ecoulement.DpPhiDt(j).perturbation = Ecoulement.DPhiDt(j).perturbation
        Ecoulement.DEtaDt(j).perturbation = -Ecoulement%DPhiDn(j)%perturbation + dot_product(Vel,GEta)
        !Ecoulement.DEtaDt(j).perturbation = -Ecoulement%DPhiDn(j)%perturbation - dot_product(GPhi,GEta) + dot_product(Vel,GEta)
		Ecoulement%DEtaDt(j)%perturbation = damp1*Ecoulement%DEtaDt(j)%perturbation - damp2*Ecoulement%Eta(j)%perturbation
    end do
else
    call CCL(t,Mesh, Ecoulement)
    call GradientFS(Ecoulement,Mesh,t)
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        select case(tdamp)
	    case(1)
		    damp1 = Mesh%Tnoeud(j)%Damping
		    damp2 = 0._rp
	    case(2) 
		    damp1 = 1._rp
		    damp2 = Mesh%Tnoeud(j)%Damping  
        end select
        DGPhi0Dz = Ecoulement.DGPhiDz(:,j).incident
		GPhi0 = Ecoulement.GPhi(:,j).incident
		GEta0 = Ecoulement.GEta(:,j).incident
		Gphi=Ecoulement.GPhi(:,j).perturbation
		GEta=Ecoulement.GEta(:,j).perturbation			
		Eta = Ecoulement.Eta(j).perturbation
		Phi = Ecoulement.Phi(j).perturbation
		call CD2Phi0DzDt(M,t,D2Phi0DzDt)
        Vel = Mesh%Tnoeud(j)%Velocity
		!! Derivees Temporelles
        Ecoulement.DpPhiDt(j).perturbation = - Ecoulement.DPhiDt(j).incident - g*Ecoulement.Eta(j).incident - 0.5_RP*dot_product(Gphi0,Gphi0) &
		                                    -g*Eta - dot_product(Gphi0,Gphi) - Eta*dot_product(DGPhi0Dz,Gphi0) - Eta*D2Phi0DzDt
		Ecoulement%DPhiDt(j)%perturbation = Ecoulement%DpPhiDt(j)%perturbation + dot_product(Vel,Gphi)
		!Ecoulement%DpPhiDt(j)%perturbation = Ecoulement%DPhiDt(j)%perturbation - dot_product(Vel,Gphi)
		Ecoulement%DPhiDt(j)%perturbation = damp1*Ecoulement%DPhiDt(j)%perturbation - damp2*Phi
		Ecoulement%DpPhiDt(j)%perturbation = damp1*Ecoulement%DpPhiDt(j)%perturbation - damp2*Phi 

		Ecoulement.DEtaDt(j).perturbation = - Ecoulement.DEtaDt(j).incident + Gphi0(3) - dot_product(Gphi0,Geta0) + Gphi(3) - dot_product(Gphi0,GEta)&
		                                     - dot_product(Gphi,Geta0)-Eta*(dot_product(DGPhi0Dz,GEta0)-DGphi0Dz(3))
        Ecoulement%DEtaDt(j)%perturbation = Ecoulement%DEtaDt(j)%perturbation + dot_product(Vel(1:2),GEta(1:2))
		Ecoulement%DEtaDt(j)%perturbation = damp1*Ecoulement%DEtaDt(j)%perturbation - damp2*Eta 
    end do
end if

! Fin
end subroutine Derive_2


!!-------------------------------------------------------------
!!
!!	Interpolation de la perturbation sur la surface au noeud 
!!	du nouveau maillage
!!
!!--------------------------------------------------------------
!!
!subroutine Interpolation( Mesh0, Ecoulement0, Mesh, Ecoulement )	
!! Paramètres
!type(TMaillage) :: Mesh0, Mesh
!type(TEcoulement) :: Ecoulement0, Ecoulement
!! Variables Locales
!integer:: j, k
!real(rp) :: d2
!real(rp), dimension(3) :: M0, M1
!! Début	
!do j=Mesh0%FS%IndFS(1),Mesh0%FS%IndFS(3)
!    Ecoulement%Eta(j)%perturbation = Ecoulement0%Eta(j)%perturbation
!    Ecoulement%DPhiDt(j)%perturbation = Ecoulement0%DPhiDt(j)%perturbation
!    Ecoulement%DEtaDt(j)%perturbation = Ecoulement0%DEtaDt(j)%perturbation
!    Ecoulement%DPhiDn(j)%perturbation = Ecoulement0%DPhiDn(j)%perturbation        
!end do
!do j=1,Mesh%Nsys
!    Ecoulement%Phi(j)%perturbation = Ecoulement0%Phi(j)%perturbation
!end do
!! Fin		
!end subroutine Interpolation

subroutine Initialisation(Ecoulement, Mesh, t)
! Parameters
type(TEcoulement) :: Ecoulement
type(TMaillage) :: Mesh
real(rp) :: t
! Locals
character(len=50) :: fileflow
integer :: iarg, ierror
logical :: get_flow
real(rp) :: t0
! Begin
call get_command_argument(4, fileflow, status=iarg)
if(iarg.eq.0)then
  get_flow = .true.
else
  get_flow = .false.
endif
if(get_flow)then
  call extract_ecoulement(Ecoulement, t0, fileflow, ierror)
  print*,' Initialisation flow à t = ',t0,' : ',fileflow
else
  call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)
  call CCL(t, Mesh, Ecoulement)
  call BodyCondition(t, Mesh, Ecoulement)
endif
! End
end subroutine Initialisation

subroutine RemeshLL(Mesh, Mesh0, t)
! Parameters
type(TMaillage) :: Mesh, Mesh0
real(rp) :: t, dt
! Locals
character(len=50) :: filemaill, num
integer :: nc, j, k
real(rp) :: Eta0
! Begin
if (Htype.ne.0) then
    ! Modification of the FS
    do j=Mesh%FS%IndFS(1),Mesh%FS%IndFS(3)
        call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
        Mesh%Tnoeud(j)%Pnoeud(3) = Eta0
    end do
    ! Modification of the tank walls
    do nc=1,Mesh%Nbody
        if (nc.ne.Int_Body) then
            do j=Mesh%Body(nc)%IndBody(1),Mesh%Body(nc)%IndBody(3)
                call CEta0(Mesh%Tnoeud(j)%Pnoeud, t, Eta0)
                !if (Crampe(-Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Mesh%Body(nc)%DimBody(3)).lt.0.1_RP) then
                    !Mesh%Tnoeud(j)%Pnoeud(3) = Mesh%Tnoeud(j)%Pnoeud(3) + Eta0
                !else
                !    Mesh%Tnoeud(j)%Pnoeud(3) = Mesh%Tnoeud(j)%Pnoeud(3) + Eta0*(1._RP-Crampe(-Mesh%Tnoeud(j)%Pnoeud(3),0._RP,Mesh%Body(nc)%DimBody(3)))
                !end if
                if (Mesh%Tnoeud(j)%NDouble.ne.0) then
                    do k=1,Mesh%Tnoeud(j)%NDouble
                        if (Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%TypeNoeud.eq.0) then
                            Mesh%Tnoeud(j)%Pnoeud = Mesh%Tnoeud(Mesh%Tnoeud(j)%Double(k))%Pnoeud
                            exit
                        end if                        
                    end do
                end if
            end do
        end if
    end do
end if
! Modification of the Body
if (is_body) then
    do nc=1,Mesh%NBody
        if (Mesh%Body(nc)%CMD(1)) then
            Mesh%Body(nc)%GBody = Mesh0%Body(nc)%GBody + Mesh%Body(nc)%MBody(1:3)
            do j=Mesh%Body(nc)%IndBody(1),Mesh%body(nc)%IndBody(3)
                Mesh%Tnoeud(j)%Pnoeud = Mesh0%Tnoeud(j)%Pnoeud + Mesh%Body(nc)%MBody(1:3) + vect_product(Mesh%Body(nc)%MBody(4:6),Mesh0%Tnoeud(j)%Pnoeud(1:3)-Mesh0%Body(nc)%GBody(1:3))
            end do
        end if
    end do
end if
!if (is_Body) then
!    if (free_body) then
!        call BodyMot(Mesh, dt, t)
!    else
!        call BodyAcc(Mesh, t, dt)
!    end if
!    call BodyVel(Mesh, t)
!end if
call GeomInit(Mesh, t)
if(idebug.gt.0) then
    Write( num, '( f0.4 )' ) t
    filemaill = 'Maillage_'//trim(num)//'.dat'
    call PlotMaill(filemaill, Mesh)
endif
! End
end subroutine RemeshLL

subroutine lissage(Maillage,Ecoulement,ierror)
  implicit none
  type(TMaillage),intent(in) :: Maillage
  type(TEcoulement),intent(inout) :: Ecoulement
  integer,intent(inout) :: ierror
! local
  integer :: j,k,kj
  integer :: NVoisin,Int_sl0,Int_sl1
  integer,dimension(NPVMax,2) :: Tvoisin
  real(rp) :: temp,temp2,lj,poids
  real(rp),dimension(3) :: M,Mj
  real(rp),parameter :: cl = 0.02  
  
  ierror = 0
  
  Int_sl0 = Maillage%FS%IndFS(1)
  Int_sl1 = Maillage%FS%IndFS(3)
  do j=Int_sl0,Int_sl1
    temp = 0._rp
    temp2 = 0._rp
    poids = 0._rp
    NVoisin = Maillage%Tnoeud(j)%NVoisin(2)-1
    M = Maillage%Tnoeud(j)%Pnoeud
    Tvoisin(1:NVoisin,1) = Maillage%Tnoeud(j)%Tvoisin(2:Nvoisin+1,1)
    Tvoisin(1:NVoisin,2) = Maillage%Tnoeud(j)%Tvoisin(2:Nvoisin+1,2)
    do k=1,NVoisin
      kj = abs(Tvoisin(k,1))
      if(Tvoisin(k,2).lt.2)then
        Mj = Maillage%Tnoeud(kj)%Pnoeud
        !if(Tvoisin(k,1).lt.0) Mj(2) = -Mj(2)
        lj = dot_product(Mj-M,Mj-M)
        temp = temp+lj*Ecoulement%Eta(kj)%perturbation
        temp2 = temp2+lj*Ecoulement%Phi(kj)%perturbation
        poids = poids+lj
      endif
    enddo
    if(poids.gt.Epsilon)then
      temp = temp/poids
      temp2 = temp2/poids
    else
      ierror = 100
      goto 9999
    endif
    Ecoulement%Eta(j)%perturbation = (1._rp-cl)*Ecoulement%Eta(j)%perturbation + cl*temp
    Ecoulement%Phi(j)%perturbation = (1._rp-cl)*Ecoulement%Phi(j)%perturbation + cl*temp2
  enddo
      
9999 continue
  if(ierror/=0)then
    if(ierror==100)then
      print*,'error #',ierror,' : division par 0'
      print*,'j = ',j,' Nvoisin = ',Nvoisin
      print*,'Tvoisin(:,1) = ',Tvoisin(:,1)
      print*,'Tvoisin(:,2) = ',Tvoisin(:,2)
    else
      write(*,99),ierror
    endif
  endif
99 format('error #',i3,' : division par 0')     
  
  
end subroutine lissage

end module BoucleTemp
