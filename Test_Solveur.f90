module Test_Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                  !
!!                                              Validation du Solveur                               !
!!                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   On choisit une géométrie simple: une cuve parapipédique, un champ engendré par une singularité 
!! extérieur au domaine.
!! Le potentiel sur la surface libre et la vitesse normale sur les parois engendrées par cette singularité
!!constituent les entrées du solveur. En sortie, les solutions numériques: potentiel sur les parois et
!!vitesse normale sur la surface libre, sont comparées aux valeurs analytiques exactes.
!!
!!
!!22/05/2013 Lucas Letournel
!!
!!
use Parameters
use StructuresDonnees
use BVP
!use Preplot
implicit none

contains

!subroutine TestSolver(Mesh)
!implicit none
!!Parameters
!type(TMaillage) :: Mesh 
!! Variables
!type(TEcoulement) :: Ecoulement
!real(rp) :: time_begin, time_end
!real(rp), allocatable :: CD(:,:), CS(:,:)
!! Begin
!call NewEcoulement(Ecoulement, Mesh%Nnoeud)
!! Initialisation de l'écoulement
!call IniEcoulement(Ecoulement, Mesh%Nnoeud, 0._RP)
!! Init Solution
!call InitValidationSolveur(Mesh, Ecoulement)
!!call InitBoxProblem(Mesh, Ecoulement)
!! Infl Coeff Calculation
!call cpu_time(time_begin)
!allocate(CD(Mesh%Nnoeud,Mesh%Nnoeud), CS(Mesh%Nnoeud,Mesh%Nnoeud))
!CD=0._RP
!CS=0._RP
!call CoeffInfl(Mesh, CD, CS, Mesh%Nnoeud)
!call cpu_time(time_end)
!print*, 'CoeffInfl', time_end-time_begin
!! Solving
!call cpu_time(time_begin)
!call SolBVP(Ecoulement, Mesh, CD, CS, Mesh%Nnoeud)
!call cpu_time(time_end)
!print*, 'SolBVP', time_end-time_begin
!! Plot
!call PlotBVP(Mesh, Ecoulement, 0._RP,.true.,'TestSolver_'//filename)
!call DelEcoulement(Ecoulement)
!deallocate(CD, CS)
!! End
!end subroutine TestSolver

subroutine InitValidationSolveur(Mesh, Ecoulement)
    
    !f2py integer*1, dimension(1000)    :: Mesh
    Type(Tmaillage)                     :: Mesh         ! Mesh.
    !f2py integer*1, dimension(1000)    :: Ecoulement
    Type(TEcoulement)                   :: Ecoulement   ! Flow parameters.
    
    real(rp)                            :: M1P, M2P     ! Norms.
    integer                             :: j            ! Loop parameter.
    real(rp)                            :: Amp          ! Amplitude of the source.
    real(rp), dimension(3)              :: M, S1, S2    ! Vectors.
    
    ! This subroutine initializes Ecoulement with a Rankine source located above the free surface.
    
    ! Amplitude.
    Amp = 1 
    S1 = 2*[0._RP, 0._RP, 1._RP] + 0*Mesh%Origine
    if (Bottom_Sym) then
        S2 = S1
        S2(3) = - 2._rp*Ldom(3) - S2(3)
    end if
    do j=1,Mesh%Nnoeud
        M = Mesh%Tnoeud(j)%Pnoeud
        M1P=1._RP/norm2(M - S1)
        
        ! Incident
        if (Bottom_Sym) then
            M2P=1._RP/norm2(M - S2)
            Ecoulement%Phi(j)%Incident = Amp*( M1P + M2P )
            Ecoulement%DPhiDn(j)%Incident = Amp*(M1P**3)*dot_product(S1-M,Mesh%Tnoeud(j)%Normale) + Amp*(M2P**3)*dot_product(S2-M,Mesh%Tnoeud(j)%Normale)
        else
            Ecoulement%Phi(j)%Incident = Amp*M1P
            Ecoulement%DPhiDn(j)%Incident = Amp*(M1P**3)*dot_product(S1-M,Mesh%Tnoeud(j)%Normale)
        end if
        
        ! Perturbation
        if (Mesh%Tnoeud(j)%typeNoeud == 0) then ! Free surface
            if(Cuve_ferme) then
                Ecoulement%Phi(j)%perturbation = Ecoulement%Phi(j)%incident
            else
                Ecoulement%Phi(j)%perturbation = Ecoulement%Phi(j)%incident*Mesh%Tnoeud(j)%Damping
            end if 
        else ! Floaters
            Ecoulement%DPhiDn(j)%perturbation = Ecoulement%DPhiDn(j)%incident
        end if
    end do
    
end subroutine InitValidationSolveur

!subroutine InitBoxProblem(Mesh, Ecoulement)
!implicit none
!! Parameters
!type(TMaillage) :: Mesh
!type(TEcoulement) :: Ecoulement
!! Local
!integer :: j
!! Begin
!if (cuve_ferme.and..false.) then
!    do j=1,Mesh%Nnoeud
!        if (Mesh%Tnoeud(j)%TypeNoeud .eq. 0) then
!            Ecoulement%Phi(j)%perturbation = Mesh%Tnoeud(j)%Pnoeud(1)
!        else
!            Ecoulement%DPhiDn(j)%perturbation = Mesh%Tnoeud(j)%Normale(1)
!        end if
!    end do
!else
!    do j=1,Mesh%Nsys
!        Ecoulement%Phi(j)%Incident = Mesh%Tnoeud(j)%Pnoeud(3)
!        Ecoulement%DPhiDn(j)%Incident = Mesh%Tnoeud(j)%Normale(3)
!        if (Mesh%Tnoeud(j)%TypeNoeud.eq.0) then
!            Ecoulement%Phi(j)%perturbation = Ecoulement%Phi(j)%Incident
!        else
!            Ecoulement%DPhiDn(j)%perturbation = Ecoulement%DPhiDn(j)%Incident
!        end if
!    end do
!end if
!
!! End
!    end subroutine InitBoxProblem
    
end module Test_Solver