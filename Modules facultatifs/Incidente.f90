!module MRampe
!use Constantes
!implicit none
!
!contains
!
!FUNCTION rampe (t)
!
!	REAL(RP) :: periode 
!	REAL(RP) :: rampe
!	REAL(RP) :: t
!
!	periode = 8.
!
!	rampe = ( 3 / (periode ** 4 ) ) * (t**4)  - ( 8/(periode**3) ) * ( t **3  ) + ( 6/(periode**2) ) * ( t **2  )
!
!END FUNCTION
!
!
!FUNCTION drampedt (t)
!
!	REAL(RP) :: periode = 8.
!	REAL(RP) :: drampedt
!	REAL(RP) :: t
!
!	drampedt = (12 / (periode ** 4 ) ) * (t**3)  - ( 24/(periode**3) ) * ( t**2  ) + ( 12 /(periode**2) ) * ( t )
!
!END FUNCTION
!
!
!
!FUNCTION d2rampedt2 (t)
!
!
!	REAL(RP) :: periode = 8.
!	REAL(RP) :: d2rampedt2
!	REAL(RP) :: t
!
!	d2rampedt2 = ( 36 / (periode ** 4 ) ) * (t**2 )  - ( 48/(periode**3) ) * (t) + ( 12/ (periode**2 ) ) 
!
!END FUNCTION
!
!end module MRampe




module IIncident

use Constantes
use Structuresdonnees
use DonneesEnvironnement
use MRampe

implicit none
    
contains

subroutine CCL(t, Ecoulement, Maillage)

! Temps courant
real(rp) :: t
type(TMaillage)			:: Maillage
type(TEcoulement)		:: Ecoulement

! Locales
integer :: k
real(rp) :: Scalaire
real(rp) , dimension(3)	:: M, GPhi0, GEta0

	
do k=1,Maillage.Nnoeud
    M=Maillage.Tnoeud(k).Pnoeud
    Ecoulement.Phi(k).incident=Phi0(M,t)
!    Ecoulement.DPhiDt(k).incident=dPhi0dt(M,t)
    call CGPhi0(M,t,GPhi0)
    Ecoulement.GPhi(:,k).incident=GPhi0
!    Ecoulement.DPhiDn(k).incident=dot_product(GPhi0,Maillage.Tnoeud(k).Normale)
!    Ecoulement.GPhi(:,k).incident=Ecoulement.GPhi(:,k).incident-Ecoulement.DPhiDn(k).incident
    
    Ecoulement.Eta(k).incident=Eta0(M,t)
!    Ecoulement.DEtaDt(k).incident=dEta0dt(M,t)
    call CGEta0(M,t,GEta0)
    Ecoulement.GEta(:,k).incident=GEta0
!    Ecoulement.DEtaDn(k).incident=dot_product(GEta0,Maillage.Tnoeud(k).Normale)
    

    Ecoulement.Phi(k).perturbation=0
    Ecoulement.DPhiDn(k).perturbation=0
    Ecoulement.GPhi(:,k).perturbation=0
    
    Ecoulement.Eta(k).perturbation=0
    Ecoulement.DEtaDn(k).perturbation=0
    Ecoulement.GEta(:,k).perturbation=0

end do

end subroutine CCL


!real(rp) function Eta0(M,t)
!! Paramètres
!!real(rp), intent(out) :: Eta0
!real(rp), dimension(3), intent(in) :: M
!real(rp), intent(in) :: t
!
!! Variables Locales
!
!Eta0 = Aphi(1)*cos(konde(1)*(M(1)*cos(dir(1))+M(2)*sin(dir(1)))-w(1)*t)
!
!end function Eta0
!
!subroutine CGEta0(M,t,Geta0)
!! Paramètres
!real(rp), dimension(3), intent(out) :: GEta0
!real(rp), dimension(3), intent(in) :: M
!real(rp), intent(in):: t
!! Variables Locales
!real(rp):: Comm
!
!Comm=-Aphi(1)*konde(1)*sin(konde(1)*(M(1)*cos(dir(1))+M(2)*sin(dir(1)))-w(1)*t)
!GEta0(1)=cos(dir(1))*Comm
!GEta0(2)=sin(dir(1))*Comm
!Geta0(3)=0
!
!end subroutine CGEta0
!
!real(rp) function Phi0(M,t)
!! Paramètres
!!real(rp), intent(out) :: Phi0
!real(rp), dimension(3), intent(in) :: M
!real(rp), intent(in) :: t
!
!Phi0 = Aphi(1)*g/w(1)*sin(konde(1)*(M(1)*cos(dir(1))+M(2)*sin(dir(1)))-w(1)*t)
!
!end function Phi0
!
!subroutine CGPhi0(M,t,GPhi0)
!! Paramètres
!real(rp), dimension(3), intent(out) :: GPhi0
!real(rp), dimension(3), intent(in) :: M
!real(rp), intent(in) :: t
!
!! Variables Locales
!real(rp) :: Comm
!real(rp), dimension(3) :: Geta0
!
!call CGEta0(M,t,GEta0)
!
!Comm=konde(1)*Aphi(1)*g/w(1)*exp(konde(1)*M(3))*cos(konde(1)*(M(1)*cos(dir(1))+M(2)*sin(dir(1)))-w(1)*t)
!GPhi0(3)= konde(1)*Aphi(1)*g/w(1)*exp(konde(1)*M(3))*sin(konde(1)*(M(1)*cos(dir(1))+M(2)*sin(dir(1)))-w(1)*t)
!GPhi0(1)=GEta0(1)*GPhi0(3)+cos(dir(1))*Comm
!GPhi0(2)=GEta0(2)*GPhi0(3)+sin(dir(1))*Comm
!
!end subroutine CGPhi0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------
!
!	Calcul de la deformee de surface libre associee 
!	au champ incident au point M et à l'instant t
!
!--------------------------------------------------

	FUNCTION eta0(M,t)

	USE DonneesEnvironnement
	USE Constantes

	IMPLICIT NONE

!	Hauteur de surface libre
	REAL(RP) :: eta0
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar
	COMPLEX :: phase
	REAL(RP):: periode

!	Pour rampe
	REAL(RP):: dPhi0dt 
	REAL(RP) , DIMENSION ( 3 ) :: M0 ! Point à z =0

	
	
	eta0=0.

periode = 8.

	IF ( t .LT. periode ) THEN

		M0(1) = M(1)
		M0(2) = M(2)
		M0(3) = 0.

		eta0 = -  dPhi0dt ( M0 , t ) / g
	ELSE
 
		DO j=1,Nhoule
			wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
			phase=i*(-konde(j)*wbar+w(j)*t)
			eta0=eta0+REAL(Aphi(j)*EXP(phase))
		END DO

!
	!Rampe



	END IF

	END FUNCTION eta0

!--------------------------------------------------
!
!	Calcul de la derivee temporelle de la
!	deformee de surface libre associee 
!	au champ incident au point (X,Y) et à l'instant t
!
!--------------------------------------------------
	FUNCTION deta0dt(M,t)


!	Hauteur de surface libre
	REAL(RP) :: deta0dt
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode , eta0 , Phiprime , Phi , Phiseconde
	COMPLEX :: phase

		REAL(RP) :: rampe, drampedt , d2rampedt2

	deta0dt=0.
	
	periode = 8.



	!Rampe

	IF ( t .LT. periode ) THEN
			
			wbar=M(1)*COS(dir(1))+M(2)*SIN(dir(1))
			phase=i*(-konde(1)*wbar+w(1)*t)
			

			Phi=REAL(i*Aphi(1)*g/w(1)*EXP(phase))
			Phiprime =REAL(- 1 * g * Aphi(1)*EXP(phase))
			Phiseconde =REAL(- g * Aphi(1)*i*w(1)*EXP(phase))


	deta0dt = - Phiseconde * rampe ( t ) / g	- 2 * drampedt ( t ) * Phiprime / g -  d2rampedt2 (t ) * Phi / g 
	ELSE

		DO j=1,Nhoule
			wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
			phase=i*(-konde(j)*wbar+w(j)*t)
			deta0dt=deta0dt+REAL(Aphi(j)*i*w(j)*EXP(phase))
		END DO




	
	END IF



	END FUNCTION deta0dt

!--------------------------------------------------
!
!	Calcul du gradient de la deformee de surface libre 
!	associee au champ incident au point (X,Y) 
!	et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGeta0(M,t,Geta0)

	USE DonneesEnvironnement
	USE Constantes

	IMPLICIT NONE

!	Gradient 
	REAL(RP),DIMENSION(3) :: Geta0
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode
	COMPLEX :: CG,phase

!	Pour la rampe

	REAL(RP),DIMENSION(3) :: GdPhi0dt
	REAL(RP),DIMENSION(3) :: M0

	DO j=1,3
		Geta0(j)=0.
	END DO

	periode = 8.


		!Rampe

	IF ( t .LT. periode ) THEN
		
		M0(1) = M(1)
		M0(2) = M(2)
		M0(3) = 0.
		CALL CGdphi0dt(M0,t,Gdphi0dt)
		Geta0 = -  Gdphi0dt  / g	
		
		
	ELSE 

		DO j=1,Nhoule
			wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
			phase=i*(-konde(j)*wbar+w(j)*t)
			CG=i*-konde(j)*Aphi(j)*EXP(phase)
			Geta0(1)=Geta0(1)+REAL(COS(dir(j))*CG)
			Geta0(2)=Geta0(2)+REAL(SIN(dir(j))*CG)
		END DO


	

	END IF


	END SUBROUTINE CGeta0

!--------------------------------------------------
!
!	Calcul du gradient de la derivee temporelle de la 
!	deformee de surface libre 
!	associee au champ incident au point M
!	et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGdeta0dt(M,t,Gdeta0dt)


!	Gradient 
	REAL(RP),DIMENSION(3) :: Gdeta0dt
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode 
	COMPLEX :: CG,phase
	REAL ( RP) , DIMENSION ( 3 ) :: Geta0
	
!	Pour la rampe
	REAL(RP) :: rampe, drampedt , d2rampedt2
	REAL ( RP) , DIMENSION ( 3 ) :: GPhi , Gphiprime , Gphiseconde
	
	
	DO j=1,3
		Gdeta0dt(j)=0.
	END DO


!rampe
periode = 8.
IF ( t .LT. periode ) THEN

		wbar=M(1)*COS(dir(1))+M(2)*SIN(dir(1))
		phase=konde(1)*(-i*wbar)+i*w(1)*t
		CG=i*w(1)*Aphi(1)*EXP(phase)
		Gphi(1)=REAL(-i*COS(dir(1))*CG)
		Gphi(2)=REAL(-i*SIN(dir(1))*CG)
		Gphi(3)=REAL(CG)

		wbar=M(1)*COS(dir(1))+M(2)*SIN(dir(1))
		phase=konde(1)*(-i*wbar)+i*w(1)*t
		CG=-w(1)*w(1)*Aphi(1)*EXP(phase)
		Gphiprime(1)=REAL(-i*COS(dir(1))*CG)
		Gphiprime(2)=REAL(-i*SIN(dir(1))*CG)
		Gphiprime(3)=REAL(CG)

		wbar=M(1)*COS(dir(1))+M(2)*SIN(dir(1))
		phase=konde(1)*(-i*wbar)+i*w(1)*t
		CG=-i*w(1)*w(1)*w(1)*Aphi(1)*EXP(phase)
		Gphiseconde(1)=REAL(-i*COS(dir(1))*CG)
		Gphiseconde(2)=REAL(-i*SIN(dir(1))*CG)
		Gphiseconde(3)=REAL(CG)




	Gdeta0dt (1) = - GPhiseconde(1) * rampe ( t ) / g	- 2 * drampedt ( t ) * GPhiprime(1) / g -  d2rampedt2 ( t ) * GPhi (1) / g 
	Gdeta0dt (2) = - GPhiseconde(2) * rampe ( t ) / g	- 2 * drampedt ( t ) * GPhiprime(2) / g -  d2rampedt2  ( t )* GPhi (2) / g 
	Gdeta0dt (3) = 0.



ELSE

	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=i*(-konde(j)*wbar+w(j)*t)
		CG=w(j)*konde(j)*Aphi(j)*EXP(phase)
		Gdeta0dt(1)=Gdeta0dt(1)+REAL(COS(dir(j))*CG)
		Gdeta0dt(2)=Gdeta0dt(2)+REAL(SIN(dir(j))*CG)
	END DO

END IF

	END SUBROUTINE CGdeta0dt

!--------------------------------------------------
!
!	Calcul du gradient de la derivee verticale de la 
!	deformee de surface libre 
!	associee au champ incident au point M 
!	et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGdeta0dz(M,t,Gdeta0dz)


!	Gradient 
	REAL(RP),DIMENSION(3) :: Gdeta0dz
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar
	COMPLEX :: CG,phase
	

	DO j=1,3
		Gdeta0dz(j)=0.
	END DO

	END SUBROUTINE CGdeta0dz

!--------------------------------------------------
!
!	Calcul du champ incident au point M à l'instant t
!
!--------------------------------------------------
!
	FUNCTION Phi0(M,t)
!
	USE DonneesEnvironnement
	USE Constantes
!
	IMPLICIT NONE
!
!	Potentiel du champ incident
	REAL(RP) :: Phi0
!
!	Point de Calcul
	REAL(RP), DIMENSION(3) :: M
!
!	Instant de Calcul
	REAL(RP) :: t
!
!	Locales
	INTEGER :: j,k
	REAL(RP) :: wbar , periode
	COMPLEX :: phase

!	Pour la rampe
	REAL(RP) :: rampe
	
!
	Phi0=0._RP

	periode = 8.

	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=konde(j)*(M(3)-i*wbar)+i*w(j)*t
		Phi0=Phi0+REAL(i*Aphi(j)*g/w(j)*EXP(phase))
	END DO
!
	!Rampe

	IF ( t .LT. periode ) THEN

		Phi0 = rampe ( t ) * Phi0

	END IF

	END FUNCTION Phi0

!--------------------------------------------------
!
!	Calcul de la derivee temporelle du
!	champ incident au point M et à l'instant t
!
!--------------------------------------------------

	FUNCTION dPhi0dt(M,t)


!	Hauteur de surface libre
	REAL(RP) :: dPhi0dt
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode , phi0
	COMPLEX :: phase
	REAL(RP) :: rampe, drampedt , d2rampedt2
	
	dPhi0dt=0.

	periode = 8.

	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=konde(j)*(M(3)-i*wbar)+i*w(j)*t
		dPhi0dt=dPhi0dt+REAL(-Aphi(j)*g*EXP(phase))
	END DO

		!Rampe

	IF ( t .LT. periode ) THEN

		dphi0dt = rampe (t )* dphi0dt + drampedt ( t ) * phi0 ( M, t )

	END IF

	END FUNCTION dPhi0dt

!--------------------------------------------------
!
!	Calcul du gradient du champ incident au point M
!	et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGphi0(M,t,Gphi0)

	USE DonneesEnvironnement
	USE Constantes

	IMPLICIT NONE

!	Gradient 
	REAL(RP),DIMENSION(3) :: Gphi0
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode
	COMPLEX :: CG,phase 
	REAL(RP) :: rampe

	DO j=1,3
		Gphi0(j)=0.
	END DO

	periode = 8.

	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=konde(j)*(M(3)-i*wbar)+i*w(j)*t
		CG=i*w(j)*Aphi(j)*EXP(phase)
		Gphi0(1)=Gphi0(1)+REAL(-i*COS(dir(j))*CG)
		Gphi0(2)=Gphi0(2)+REAL(-i*SIN(dir(j))*CG)
		Gphi0(3)=Gphi0(3)+REAL(CG)
	END DO
	!Rampe

	IF ( t .LT. periode ) THEN

		GPhi0 = rampe (t )* GPhi0

	END IF

	END SUBROUTINE CGphi0

!--------------------------------------------------
!
!	Calcul du gradient de la derivee verticale 
!	du champ incident au point M 
!	et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGdPhi0dz(M,t,GdPhi0dz)


!	Gradient 
	REAL(RP),DIMENSION(3) :: GdPhi0dz
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode
	COMPLEX :: CG,phase
	REAL(RP) :: rampe

	DO j=1,3
		GdPhi0dz(j)=0.
	END DO

	periode = 8.
	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=konde(j)*(M(3)-i*wbar)+i*w(j)*t
		CG=konde(j)*i*w(j)*Aphi(j)*EXP(phase)
		Gdphi0dz(1)=Gdphi0dz(1)+REAL(-i*COS(dir(j))*CG)
		Gdphi0dz(2)=Gdphi0dz(2)+REAL(-i*SIN(dir(j))*CG)
		Gdphi0dz(3)=Gdphi0dz(3)+REAL(CG)
	END DO

		!Rampe

	IF ( t .LT. periode ) THEN

		GdPhi0dz = rampe (t ) * GdPhi0dz

	END IF

	END SUBROUTINE CGdPhi0dz

!--------------------------------------------------
!
!	Calcul du gradient de la derivee temporelle du
!	champ incident au point M et à l'instant t
!
!--------------------------------------------------

	SUBROUTINE CGdphi0dt(M,t,Gdphi0dt)


!	Gradient 
	REAL(RP),DIMENSION(3) :: Gdphi0dt
!	Point de Calcul
	REAL(RP),DIMENSION(3) :: M
!	Instant de Calcul
	REAL(RP) :: t
!	Indices
	INTEGER :: j,k
!	Locales
	REAL(RP) :: wbar , periode
	COMPLEX :: CG,phase 
	REAL(RP), DIMENSION (3) :: GPhi0
		REAL(RP) :: rampe, drampedt , d2rampedt2

	DO j=1,3
		Gdphi0dt(j)=0.
	END DO


	periode = 8.

	DO j=1,Nhoule
		wbar=M(1)*COS(dir(j))+M(2)*SIN(dir(j))
		phase=konde(j)*(M(3)-i*wbar)+i*w(j)*t
		CG=-w(j)*w(j)*Aphi(j)*EXP(phase)
		Gdphi0dt(1)=Gdphi0dt(1)+REAL(-i*COS(dir(j))*CG)
		Gdphi0dt(2)=Gdphi0dt(2)+REAL(-i*SIN(dir(j))*CG)
		Gdphi0dt(3)=Gdphi0dt(3)+REAL(CG)
	END DO



		!Rampe

	IF ( t .LT. periode ) THEN
	

		wbar=M(1)*COS(dir(1))+M(2)*SIN(dir(1))
		phase=konde(1)*(M(3)-i*wbar)+i*w(1)*t
		CG=i*w(1)*Aphi(1)*EXP(phase)
		Gphi0(1)=REAL(-i*COS(dir(1))*CG)
		Gphi0(2)=REAL(-i*SIN(dir(1))*CG)
		Gphi0(3)=REAL(CG)


		Gdphi0dt = rampe (t)* Gdphi0dt + drampedt(t) * Gphi0 

	END IF


	END SUBROUTINE CGdphi0dt

end module Iincident