module Parameters

use Constantes

integer :: aff ! a suppr (utilise pour debugue

real(rp), dimension(4,5000) :: temp_fort
! =========================================================================
! -------------------------------------------------------------------------
!   GENERAL PARAMETERS
! -------------------------------------------------------------------------
!
! -- Project parameters ---------------------------------------------------
character(len=50)                   :: filename
logical, parameter                  :: is_Licence = .false.
! -- Wave parameters ------------------------------------------------------
integer                             :: Htype                        ! Type de Houle incidente (0 : Nulle; 1: Airy ; 2: RF).
real(rp)                            :: directionRF                  ! Direction of the RF wave.
real(rp)                            :: PhiWaveRF                    ! Phase of the RF wave.
real(rp)                            :: profondeur
logical                             :: infinite_depth 
! -- Rienecker & Fenton wave parameters --
character(len=50)                   :: filecoeff
integer                             :: NMAX
! -- Airy wave parameters ----------------
integer                             :: Nhoulemx
integer                             :: Nhoule		                                ! Nombre de composantes de la houle (max 10).
real(rp), dimension(10)             :: w, konde, lambda, DirectionAiry, Aphi,PhiWave! Wave frequency, wave number, wave length, wave direction, wave amplitude and phase of each Airy wave (there are lots of "wave" in this sentence, maybe it is too vague...).
real(rp)                            :: Trig
! -- Mesh Parameters ------------------------------------------------------
real(rp)                            :: Prof                         ! Profondeur
logical                             :: Symmetry                     ! Symétrie
logical                             :: bottom_sym                   ! Symetrisation du fond
logical                             :: cuve_ferme                   ! Domaine fermé = 1, Domaine ouvert = 0
integer                             :: Int_body                     ! Indice du corps dans le maillage
integer                             :: Mesh_Type                    ! Methode generation du maillage (1:Struct. , 2:Frontale)
integer                             :: decalage
integer                             :: Nth                          ! Nombre de pas de discrétisation en theta pour la SL
integer                             :: Nr                           ! Nombre de pas de discrétisation en rayon pour la SL
real(rp)                            :: dx1, dx3, dx4                ! Paramètres de discrétisation (longueur de maillage) for the external boudary of the free surface (dx1) and the bottom of the floaters (dx3). dx4 is not used.
real(rp)                            :: dx2Domain                    ! Panel size in the center of the tank if there is no body.
real(rp)                            :: d_bl                         ! Taille de couche limite autour des flotteurs.
real(rp)                            :: LAbs                         ! Longueur de la plage absorbante
real(rp)                            :: hrefx                        ! Size according to x of the Cartesian grid except on the free surface.
real(rp)                            :: hrefy                        ! Size according to y of the Cartesian grid except on the free surface.
real(rp)                            :: hrefxFS                      ! Size according to x of the Cartesian grid on the free surface.
real(rp)                            :: hrefyFS                      ! Size according to y of the Cartesian grid on the free surface.
real(rp)                            :: dsi                          ! Pas de calcul pour la recherche de la courbe d'intersection SL/SM
real(rp)                            :: NormOverDref                 ! Coefficient in front of delta_ref in the equation (4.11) of the PhD thesis of CC (Camille used 0.8 except for the axisym geometry: 1.0).
real(rp)                            :: Ra                           ! Rayon d'arondi des angles de bords (version experimentale) 
logical                             :: DeformMesh                   ! True = Mesh can be deformed, false otherwise.
real(rp),   parameter               :: dS2max = 9.
! -- Numerical parameters --------------------------------------------------
logical                             :: CCI                          ! Partial computation of the influence coefficients (T/F).
logical                             :: is_BS                        ! true--> Bspline, false --> Linear discretization.
logical                             :: is_BS_GS                     ! true--> Bspline, false --> Linear discretization (pour le calcul du gradient surfacique)


real(rp)                            :: DPNormals                    ! Dot product of the FS normal and the body normal at the intersection curve. If nFS.nBody >= DPNormals (usually around 1), the linear system to update the gradient is ill-conditionned (cf the subroutine GradientIntersection).
integer                             :: SplineOrder                  ! Orders of the B-Spline(0, 3 or 4) and order of the neighboring points.
integer                             :: Theory                       ! Potential Theory used (0: Fully Linear (FS and Body),
                                                                    !                        1: Body Exact (linear FS, exact Body)
                                                                    !                        2: Weak-Scatterer (weakly NL FS, exact Body)
                                                                    !                        3: Fully NL (à venir ?! :D).
logical                             :: lineaireFS                   ! Linear FS condition, FS position on z=0.
logical                             :: lineaireBody                 ! Body position fixed on mean position.
logical                             :: RemeshFS                     ! Remeshing of the FS (T or F).
logical                             :: ForcedRemeshFS               ! Remeshing of all the mesh at each time step (T or F).
logical                             :: DeformFS                     ! Deformation of the FS.
logical                             :: DeformBody                   ! Deformation of the body.
logical                             :: RK_fige                      ! Maillage figé dans le schéma RK4.
logical                             :: is_latching                  ! Latching method.
real(rp)                            :: t0                           ! Décalage par rapport à l'instant initial de la houle.
real(rp)                            :: dt                           ! Pas de temps pour RK4.                
integer                             :: nt                           ! Nombre de pas de temps.
real(rp)                            :: T1, T2                       ! Parametre de rampe pour la condition limite.
real(rp)                            :: Tout                         ! Instant a partir duquel on stocke le resultat.
integer                             :: nout                         ! Nombre de pas de temps entre chaque ecriture.
integer                             :: Solv                         ! Schema numerique utilise pour la resolution du probleme aux limite (0: GMRES, 1: LU).
real(rp)                            :: TolGMRES                     ! Critere de convergence de la methode numerique.
integer                             :: nliss_input                  ! Number of time steps between two calls of the subroutine LissageGaussien, value in the *.ws input file.
integer                             :: nliss                        ! Number of time steps between two calls of the subroutine LissageGaussien, value used in the temporal looped. nliss_input and nliss are different in case of a body crossing the free surface.
logical                             :: hydrostatique                ! Computation of the hydrostics (weight and Archimede): True or not: False. 
integer                             :: nWP                          ! Number of wave probes.
real(rp)                            :: ro                           ! Water density.
real(rp)                            :: g                            ! Acceleration of gravity.
integer                             :: NBodies                      ! Number of bodies.
real(rp)                            :: eps
logical, parameter                  :: FiniteDifference_ForceCalcul= .false. ! Force calculation using finite differences in forced motion.
logical                             :: FreeBodies                   ! = 1 if at least one body is free to move, 0 if all InputData%free_body = false
integer                             :: NThreads                     ! Number of threads for OpenMP.
! -- Body Motion parameters ------------------------------------------------
integer                             :: Tcase                        ! Choix de CL pour déplacement, vitesse et accélération du corps dans Incident      Useless !!!


! -----PARAREAL---------

logical                                     :: bool_coarse
logical                                     :: bool_coarse_init


! ==========================================================================
! --------------------------------------------------------------------------
!   OUTPUT PARAMETERS 
! --------------------------------------------------------------------------

integer                             :: iinfodiv = 0                 ! Pour enelever plein de sorties
integer                             :: iprint                       ! Sortie Ecran (O/N)

integer                             :: idebug                       ! Information de debugage (O/N)
integer                             :: iwbound                      ! Sortie geometrie des frontières (O/N)
integer                             :: iwfront                      ! Sortie tableau du front de la méthode frontale (O/N)

integer,parameter                   :: ioMaillage = 3               ! Maillage.dat
integer,parameter                   :: ioMeshBody = 300             ! Mesh_body_.dat
integer,parameter                   :: ioMeshBodSTL = 400           ! Mesh_body_.stl
integer,parameter                   :: ioTmpBodyMesh = 10           ! Temporary body mesh file.

logical                             :: iwevol                       ! Sortie des étapes de création du maillage (O/N)                    
integer,parameter                   :: ioevol = 73                  ! If you change the number 73, please change it in AdvanceFront as well.

logical                             :: iwevolRemesh                 ! Sortie des étapes de création du maillage lors du remaillage (O/N) 
integer,parameter                   :: ioevolRemesh = 733           ! If you change the number 733, please change it in AdvanceFront as well.

logical                             :: iwCartGrid                   ! Tecplot file for the Cartesian grid
integer,parameter                   :: ioCartGrid = 74

integer                             :: iwmesh1
logical                             :: oplot                        ! Plot ou pas?

integer,parameter                   :: iomesh = 37                  ! Sortie maillage a chaque iteration
logical,parameter                   :: iwmesh = .false.

integer,parameter                   :: iometric = 40                ! Sortie mesure qualite maillage
logical,parameter                   :: iwmetric= .false.

integer,parameter                   :: iopress = 48, iopress2 = 488 ! Sortie forces hydrodynamiques
logical,parameter                   :: iwpress = .true.

integer,parameter                   :: ioenergy = 49                ! Sortie de l'energie
logical                             :: iwenergy

integer,parameter                   :: iophi = 50                   ! Sortie potentiel 3D
logical,parameter                   :: iwphi = .false.

integer,parameter                   :: ioBVP = 53                   ! Sortie solveur BVP
logical,parameter                   :: iwBVP = .false.

integer,parameter                   :: ioTest_BVP = 533             ! Sortie solveur Test_BVP
logical,parameter                   :: iwTest_BVP = .true.

integer,parameter                   :: ioDerive = 54                ! Sortie Dérivées SL
logical,parameter                   :: iwDerive = .false.

integer,parameter                   :: iodphidt = 50, iodphidt2 = 550! Sortie DPhiDt 3D
logical,parameter                   :: iwdphidt = .false.

logical                             :: icheck                       ! Sortie controle precision maillage
integer,parameter                   :: ioint = 51, ioqual = 52

integer,parameter                   :: iopress3D = 24               ! Sortie Pression 3D
logical,parameter                   :: iwpress3D = .false.

integer,parameter                   :: ioabody = 22                 ! Sortie Acceleration (cas mouvement libre)
logical,parameter                   :: iwabody = .false.

integer,parameter                   :: iots = 21                    ! Sortie TS (cas mouvement libre)
logical,parameter                   :: iwts = .false.

integer,parameter                   :: iombody = 3000               ! Sortie Mouvement du Corps
logical                             :: iwmbody 

integer,parameter                   :: iophit = 27                  ! Sortie terme du problème aux limites sur phit
logical,parameter                   :: iwphit = .false.

logical                             :: iwhoule

integer,parameter                   :: ioaxisym = 28
character(len=50)                   :: file_axisym 

integer,parameter                   :: iosize = 29
logical,parameter                   :: iwsize = .false.

integer,parameter                   :: iomonitor = 30
logical                             :: iwmonitor 

integer,parameter                   :: iomass = 38
logical,parameter                   :: iwmass = .false.

integer,parameter                   :: iowave = 3900
logical,parameter                   :: iwwave = .true.              ! Wave probes

integer,parameter                   :: ioLoads = 4000
logical,parameter                   :: iwLoads = .true.             ! Loads (weight, hydrostatics, hydrodynamics, PTO and Morison).

integer,parameter                   :: ioparam = 122                ! *.param input file.
integer,parameter                   :: iogeom = 123                 ! *.geom input file.
integer,parameter                   :: ioinertia = 124              ! Inertia file to read.
integer,parameter                   :: ioWriteInertia = 125         ! inertia file to write (when there is nothing to read!).

integer,parameter                   :: ioIntersection = 41          ! Intersection curves.
logical,parameter                   :: iwIntersection = .true.

integer,parameter                   :: ioIntersectionDebug = 42     ! Intersection curves debug.
logical,parameter                   :: iwIntersectionDebug = .true.

integer,parameter                   :: ioState = 1789               ! State vector output file. Why 1789? Because the input state file is a Revolution!
logical,parameter                   :: iwState = .true.

logical                             :: iwmatv                       ! Stockage de la vitesse de deplacement du maillage dans un fichier.

logical                             :: t_spline                     ! Flag for the test of B-Spline gradient calculations.
logical                             :: t_Solver                     ! Flag for the test of the resolution of the BVP solver.

! ========================================================================
! ------------------------------------------------------------------------
!   FLOATING BODY PARAMETERS 
! ------------------------------------------------------------------------
logical                             :: is_body                      ! If there is a body: 1, otherwise: 0.
! -- Body geometry -------------------------------------------------------
logical                             :: iFixPoint
real(rp), dimension(3)              :: FixPointPos(1:3) 
! -- Domain geometry -------------------------------------------------------
integer                             :: idtype                       ! Type de geometrie pour le domaine (1:rectangulaire, 2: circulaire)
real(rp), dimension(5)              :: Ldom                         ! Longueur, largeur et profondeur du domaine, puis discrétisation pour le domaine de type cylindrique(dlong,larg,dprof,dxgeom,drayon)

real(rp)                            :: Forward_Velocity             ! Forward velocity for EACH body (all the body msut have the same forward velocity).

! ========================================================================
! ------------------------------------------------------------------------
! SUPPLEMENTARY PARAMETERS
! ------------------------------------------------------------------------

character(len=50)                   :: filenemoh = 'nemoh.dat'
integer                             :: idref                        ! = 0 Sparse Laplacian to create the mesh = 1 : Full Laplacian , = 2 : linear approximation.
logical                             :: is_immerged                  ! = true no body is piercing (no intersection curve floater/free surface), false otherwise (at least one body pierces the free surface). 

! Spline parameters
integer, dimension(0:4)             :: Nordre                       ! Spline type and order (cf Exec.f90). Pay attention, this vector starts at 0 and not 1.
real(rp), parameter                 :: RSpline = 2                  ! Useless

end module Parameters
