	module Constantes
    
    ! Bibliotheque du compilateur utilisee
	use ifport
    
    ! Definition de symboles pour les types reels (rp) et complexes (CP)
    
    ! Les reels, simple ou double precision
	integer, PARAMETER :: SP = KIND(1.0)
	integer, PARAMETER :: DP = KIND(1.0D0)
    
    ! Les complexes, simple ou double precision
	integer, PARAMETER :: SPC = KIND(COMPLEX(1.0,1.0,4))
	integer, PARAMETER :: DPC = KIND(COMPLEX(1.0D0,1.0D0,8))
    
    ! Les types courants
	integer, PARAMETER :: RP = DP
	integer, PARAMETER :: CP = DPC
    
    ! Le zero machine
    !real(rp),PARAMETER :: Epsilon = 1.0E-5  ! CC
	!real(RP),PARAMETER :: Epsilon2 = 1.0E-2 ! CC
    !real(rp),parameter :: Epsgeom = 1.0E-3  ! CC
	!real(rp),PARAMETER :: Epsilon = 1.0E-10
	!real(RP),PARAMETER :: Epsilon2 = 1.0E-10
    !real(rp),parameter :: Epsgeom = 1.0E-10
    real(rp),PARAMETER :: Epsilon = 1.0E-7
	real(RP),PARAMETER :: Epsilon2 = 1.0E-7
    real(rp),parameter :: Epsgeom = 1.0E-7
    
    
    
    integer :: print_debug = 0 ! A sUPPRIMER
    
    ! Erreur allocation
	integer :: ierr
    
    ! Les constantes mathematiques usuelles i, pi, 2pi, pi/2, 1./pi, racine de 2.
	complex(CP), PARAMETER :: i     = CMPLX(0.0_rp, 1.0_rp, KIND=CP)
	real(rp), PARAMETER    :: PI    = 3.141592653589793238462643383279502888419_rp ! pi
	real(rp), PARAMETER    :: PIO2  = 1.570796326794896619231321691639751442098_rp ! Pi/2
	real(rp), PARAMETER    :: TWOPI = 6.283185307179586476925286766559005768394_rp ! 2*pi
	real(rp), PARAMETER    :: SQ2   = 1.414213562373095048801688724209698078569_rp ! pi^2
	real(rp), PARAMETER	   :: SPI	= 1._RP/PI ! 1/pi
    
    ! Une tolerance sur les distances physiques des noeuds
	real(rp),PARAMETER :: Tol = 1.0E-02_RP
    
    ! Une valeur de dissipation numerique
	real(rp),PARAMETER :: Nnu = 0.1_RP              ! ???
    
	integer, PARAMETER :: Ncnxmx = 100            	! Nombre maximal de facettes connectees a un sommet
	integer, parameter :: NPVMax = 500           	! Nombre maximal de Points voisins à un point (BSPline)
	integer, parameter :: NFCPMax = 5000        	! Nombre maximal de Facettes du champ proches (Coeff Infl)
!	integer, parameter :: NFVMax = Ncnxmx       	! Nombre maximal de Facettes voisines à une facette (geom)
	integer, PARAMETER :: PointMax = 11000     		! Nombre maximal de Points du maillage
	integer, PARAMETER :: FacetteMax = 2*PointMax   ! Nombre maximal de Facettes du maillage
	
	integer, parameter :: nvmax = 30                ! Nombre max de point voisin à un point donné
    integer, parameter :: nfvmax = 5                ! Nombre max de face contenant le point
    
    integer, parameter :: NHMAX = 200               ! ???
    
    ! Dichotomy constants used in the computation of the intersection curves.
    real(rp), parameter :: Residu_Dichotomie = 0.000001 ! Residual of the dichotomy method.
    integer, parameter  :: itmax_Dichotomie = 800       ! Maximum number of loops with the dichotomy method.
    
	end module Constantes	
