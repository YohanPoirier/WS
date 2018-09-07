module SolvNum
use Constantes
use GeomStruct
use GeomFonct
use GeomDisque
use GeomCylindre
use GeomPlan
use GeomSphere
use GeomCone
use GeomAxiSym
use GeomWigley
implicit none

interface intersection_3D_2
  module procedure intersection_3D_cylindre,&
&                  intersection_3D_disque,&
&                  intersection_3D_plan,&
&                  intersection_3D_sphere,&
&                  intersection_3D_cone,&
&                  intersection_3D_axisym,&
&                  intersection_3D_Wigley
end interface intersection_3D_2

interface point_inter
  module procedure point_inter_cyl,&
&                  point_inter_disque,&
&                  point_inter_plan,&
&                  point_inter_sph,&
&                  point_inter_cone,&
&                  point_inter_Wigley
end interface

contains

subroutine intersect_line_surf(P1,P2,Pout,t)
    
    !f2py integer*1, dimension(1000)    :: P1,P2
    type(point),intent(in)              :: P1,P2        ! End points.
    !f2py integer*1, dimension(1000)    :: Pout         
    type(point),intent(inout)           :: Pout         ! Result.
    real(rp),intent(in)                 :: t            ! Current time.
    
    !f2py integer*1, dimension(1000)    :: P
    type(point)                         :: P            ! Point.
    integer                             :: n            ! Loop parameter.
    real(rp)                            :: s,s1,s2      ! Curvilinear coordinates.
    real(rp)                            :: z1,z2,d      ! Vertical distances.
    real(rp)                            :: deltaz,eta   ! Delta z and the wave elevation.
    
    ! This subroutine computes the intersection between a line and a surface in using a dichomoty method.
    
    ! Verification
    call CEta0(P1%coord,t,z1)
    call CEta0(P2%coord,t,z2)
    z1 = P1%coord(3)-z1
    z2 = P2%coord(3)-z2
    
    d = min(abs(z1),abs(z2))
    
    if (idebug>0 .and. (z1.gt.0 .and. z2.lt.0 .or. z1.lt.0 .and. z2.gt.0)) then
        print*,"Line intersect the free surf."
    elseif(iprint ==1) then
        print*,"warning : line doesn't intersect free surf"
        print*,"        : P1 = ",P1%coord," ; z1 = ",z1
        print*,"        : P2 = ",P2%coord," ; z2 = ",z2  
    endif
    
    ! Initialisation variables
    n = 0 ; s1 = 0 ; s2 = 1.
    
    ! Boucle dichotomie
    do while (d.gt.Residu_Dichotomie .and. n.le.itmax_Dichotomie) ! itmax_Dichotomie = 800 (cf Constantes.f90).
        n = n+1
        s = (s1+s2)/2.
        call Computation_eq_line(P1,P2,s,P%coord) ! Previously only P
        call CEta0(P%coord,t,eta)
        deltaz = P%coord(3)-eta
    
        if (deltaz.gt.0 .and. z2 > 0 .or. deltaz <= 0 .and. z2 <= 0.) then
            s2 = s
        else
            s1 = s
        endif
        d = abs(deltaz)
    enddo
    
    ! Result
    if (d.le.Residu_Dichotomie) then
        if(n.gt.0)then
            ! The dichotomy loop was used.
            Pout = P
        elseif(n.eq.0)then
            ! The dichotomy loop was not used.
            if(abs(z1).le.abs(z2))then
                Pout = P1
            else ! (abs(z1).gt.abs(z2))
                Pout = P2
            end if
        end if
    else
        print*,""
        print*,'intersect_line_surf : solver cannot find P0'
        print*,"P1 = ",P1
        print*,"P2 = ",P2
        print*,"d = ",d
        print*,"n = ",n
        if(n.ne.0)then
            print*,"Pout = ",Pout
        end if
        print*,"t = ",t
        pause
    endif
    
end subroutine intersect_line_surf

subroutine intersect_curve_surf(P1,P2,cercle,Pout,t)
    
    real(rp),dimension(3),intent(in)    :: P1,P2        ! End points (cart).
    !f2py integer*1, dimension(1000)    :: cercle
    type(GCercle),intent(in)            :: cercle       ! Circle.
    real(rp),dimension(3),intent(inout) :: Pout         ! Result.
    real(rp),intent(in)                 :: t            ! Current time.
    
    real(rp),dimension(3)               :: P
    integer                             :: n
    real(rp)                            :: s,s1,s2
    real(rp)                            :: z1,z2,d,eta
    real(rp)                            :: deltaz
    
    real(rp) :: theta1,theta2,theta
    real(rp),dimension(3)               :: Q1,Q2
    
    ! This subroutine computes the intersection between a circle and a surface in using a dichomoty method.
    
    call CEta0(P1,t,z1)
    call CEta0(P2,t,z2)
    z1 = P1(3)-z1
    z2 = P2(3)-z2
    
    d = min(abs(z1),abs(z2))
    
    if (idebug>0 .and. (z1.gt.0 .and. z2.lt.0 .or. z1.lt.0 .and. z2.gt.0)) then
        print*,"Line intersect the free surf."
    elseif(iprint ==1) then
        print*,"warning : line doesn't intersect free surf"
        print*,"        : P1 = ",P1," ; z1 = ",z1
        print*,"        : P2 = ",P2," ; z2 = ",z2  
    endif
    
    ! Initialisation variables
    n = 0 ; s1 = 0._rp ; s2 = 1._rp
    
    ! Boucle dichotomie
    do while (d.gt.Residu_Dichotomie .and. n.le.itmax_Dichotomie) ! itmax_Dichotomie = 800 (cf Constantes.f90).
        n = n+1
        s = (s1+s2)/2.
        call eq_curve(cercle,P1,P2,s,P) ! seule modif par rapport a intersect_line_surf
        call CEta0(P,t,eta)
        deltaz = P(3)-eta
        if (deltaz.gt.0 .and. z2 > 0 .or. deltaz <= 0 .and. z2 <= 0.) then
            s2 = s
        else
            s1 = s
        endif
        d = abs(deltaz)
    enddo
    
    ! Result
    if (d.le.Residu_Dichotomie) then
        if(n.gt.0)then
            ! The dichotomy loop was used.
            Pout = P
        elseif(n.eq.0)then
            ! The dichotomy loop was not used.
            if(abs(z1).le.abs(z2))then
                Pout = P1
            else ! (abs(z1).gt.abs(z2))
                Pout = P2
            end if
        end if
    else
        print*,""
        print*,'intersect_curve_surf : solver cannot find P0'
        print*,"P1 = ",P1
        print*,"P2 = ",P2
        print*,"d = ",d
        print*,"n = ",n
        if(n.ne.0)then
            print*,"Pout = ",Pout
        end if
        print*,"t = ",t
        pause
    endif
    
end subroutine intersect_curve_surf

subroutine create_tabcurve_intersect(P1,P2,t,cercle,dx,tab_point,np,iface)
  implicit none
  !f2py integer*1, dimension(1000) :: P1,P2
  type(point),intent(inout) :: P1,P2
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000) :: cercle
  type(GCercle),intent(in) :: cercle
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
  integer,intent(in) :: iface
! local  
  !f2py integer*1, dimension(1000) :: A,B,M,Q
  type(point) :: A,B,M,Q
  real(rp),dimension(3) :: P0
  real(rp) :: eta,eta1,eta2,L
  real(rp) :: r,angle
  integer :: k,N
  integer :: ierror
   
  ! This subroutine creates the table of the intersection points for circles.
  
 ! np = 0
  ierror = 0
  
  r = cercle%rayon
  M = P1
  Q = P2
  !theta1 = cart2cer(P1%coord,cercle)
  !theta2 = cart2cer(P2%coord,cercle)
  call loc2cart(P1%coord,cercle%repere,M%coord)
  call loc2cart(P2%coord,cercle%repere,Q%coord)
  
  !if(abs(theta2-theta1).lt.Epsilon)then
  !  theta2 = TWOPI
  !elseif(theta2 < theta1)then
  !  theta2 = theta2+2.*PI
  !endif  
  
  !angle = theta2-theta1
  angle = cercle%theta
  L = angle*r
  N = int(L/dx)
  
  call CEta0(M%coord,t,eta1)
  call CEta0(Q%coord,t,eta2)
  
  if(abs(M%coord(3)-eta1) .lt. Epsilon) then
    M%bf = 2
    np=np+1
    tab_point(np)=M%coord
    tab_point(np)%nface = 2
    tab_point(np)%face(1:2) = [HouleRF%index,iface]
    tab_point(np)%bf = 2
  elseif(M%coord(3) .gt. eta1) then
    M%bf = 0
  else
    M%bf = 1
  endif
  
  if(abs(Q%coord(3)-eta2) .lt. Epsilon) then
    Q%bf = 2
    np=np+1
    tab_point(np)=Q%coord
    tab_point(np)%nface = 2
    tab_point(np)%face(1:2) = [HouleRF%index,iface]
    tab_point(np)%bf = 2
  elseif(Q%coord(3) .gt. eta2) then
    Q%bf = 0
  else
    Q%bf = 1
  endif
  
  A = M
  do k=1,N-1
    call eq_curve(cercle,M%coord,Q%coord,dble(k)/dble(N),B%coord) ! Previously only B
    
    call CEta0(B%coord,t,eta)
    if(abs(B%coord(3)-eta).lt.Epsilon .or. B%coord(3).ge.eta .and. A%coord(3).le.eta .or. B%coord(3).le.eta .and. A%coord(3).ge.eta)then
    !if(abs(B%coord(3)-eta).lt.Epsilon)then
      B%bf = 2
      np = np+1
      tab_point(np) = B%coord
      tab_point(np)%nface = 2
      tab_point(np)%face(1:2) = [HouleRF%index,iface]
    elseif(B%coord(3).gt.eta)then
      B%bf = 0
    else
      B%bf = 1
    endif
    if (A%bf /= B%bf .and. A%bf/=2 .and. B%bf/=2) then
      call intersect_curve_surf(A%coord,B%coord,cercle,P0,t)
      np = np+1
      tab_point(np) = P0
      tab_point(np)%nface = 2
      tab_point(np)%face(1:2) = [HouleRF%index,iface]
    endif
    A = B
  enddo
!
  B = Q
  if(A%bf/=B%bf .and. A%bf/=2 .and. B%bf/=2) then
    call intersect_curve_surf(A%coord,B%coord,cercle,P0,t)
    np = np+1
    tab_point(np) = P0
    tab_point(np)%nface = 2
    tab_point(np)%face(1:2) = [HouleRF%index,iface]
  endif
  
  do k =1,np-1
  ! supprimer deux points d'intersection trop proches (<2*dx)
  enddo  
  
9999 continue  
  if(ierror/=0)then
    write(*,100),ierror
100 format('error #',i3," : Creation de la table des points d'intersection")
  endif    
!  
end subroutine create_tabcurve_intersect

subroutine create_tab_intersect(P1,P2,t,dx,tab_point,np,iface)
    
    !f2py integer*1, dimension(1000)        :: P1,P2
    type(point),intent(inout)               :: P1,P2            ! End points.
    real(rp),intent(in)                     :: t,dx             ! Current time and panel discretization.
    !f2py integer*1, dimension(1000)        :: tab_point
    type(point),dimension(*),intent(inout)  :: tab_point        ! List of points of the intersection curves.
    integer,intent(inout)                   :: np               ! Number of points in tab_point.
    integer,intent(in)                      :: iface            ! Number of the face.
    
    !f2py integer*1, dimension(1000)        :: A,B,P0
    type(point)                             :: A,B,P0           ! Points.
    real(rp)                                :: eta,eta1,eta2,L  ! Wave elevations and length of the line segment P1P2.
    integer                                 :: k,N              ! Loop parameters.
    real(rp),dimension(3)                   :: V                ! Vector.
    integer                                 :: ierror           ! Error flag.
   
    ! This subroutine creates the table of the intersection points for lines.
    
    np = 0
    ierror = 0
    
    ! Parameters for the dichotomy method.
    V = P2%coord - P1%coord ! P1P2
    L = sqrt(V(1)*V(1)+V(2)*V(2)+V(3)*V(3)) ! ||P1P1||
    N = int(L/dx)
    
    ! Wave elevations
    call CEta0(P1%coord,t,eta1)
    call CEta0(P2%coord,t,eta2)
    
    ! Intersection with P1
    if(abs(P1%coord(3)-eta1) .lt. Epsilon) then
        P1%bf = 2
        np=np+1
        tab_point(np)=P1%coord
        tab_point(np)%nface = 2
        tab_point(np)%face(1:2) = [HouleRF%index,iface]
    elseif(P1%coord(3) .gt. eta1) then
        P1%bf = 0
    else
        P1%bf = 1
    endif
        
    ! Intersection with P2
    if(abs(P2%coord(3)-eta2) .lt. Epsilon) then
        P2%bf = 2
        np=np+1
        tab_point(np)=P2%coord
        tab_point(np)%nface = 2
        tab_point(np)%face(1:2) = [HouleRF%index,iface]
    elseif(P2%coord(3) .gt. eta2) then
        P2%bf = 0
    else
        P2%bf = 1
    endif
    
    ! Dichotomy method between A and B
    A = P1
    do k=1,N-1
        call Computation_eq_line(P1,P2,dble(k)/dble(N),B%coord) ! Previously only B
        call CEta0(B%coord,t,eta)
        
        ! Searching the point which intersects P1P2 and the free surface.
        !if(abs(B%coord(3)-eta).lt.Epsilon .or. B%coord(3).ge.eta .and. A%coord(3).le.eta .or. B%coord(3).le.eta .and. A%coord(3).ge.eta)then
        if(abs(B%coord(3)-eta).lt.Epsilon)then
            B%bf = 2
            np = np+1
            tab_point(np) = B%coord
            tab_point(np)%nface = 2
            tab_point(np)%face(1:2) = [HouleRF%index,iface]
        elseif(B%coord(3).gt.eta)then
            B%bf = 0
        else
            B%bf = 1
        endif
        if (A%bf /= B%bf .and. A%bf/=2 .and. B%bf/=2) then
            call intersect_line_surf(A,B,P0,t)
            np = np+1
            tab_point(np) = P0%coord   
            tab_point(np)%nface = 2
            tab_point(np)%face(1:2) = [HouleRF%index,iface]
        endif
        A = B
        if(np.ne.0)then
            go to 9999 ! There is already an intersection point, it is useless to seach another one.
        end if
    enddo
    
    B = P2
    if(A%bf/=B%bf .and. A%bf/=2 .and. B%bf/=2) then
        call intersect_line_surf(A,B,P0,t)
        np = np+1
        
        tab_point(np) = P0%coord
        tab_point(np)%nface = 2
        tab_point(np)%face(1:2) = [HouleRF%index,iface]
    endif
    
    do k =1,np-1
        ! supprimer deux points d'intersection trop proches (<2*dx)
        ! Done in update_point_inter? (PYW)
    enddo  
    
    9999 continue  
    if(ierror/=0)then
        write(*,100),ierror
        100 format('error #',i3," : Creation de la table des points d'intersection")
    endif    
    
end subroutine create_tab_intersect

subroutine point_inter_cyl(cyl,t,dx,tab_point,np)
  implicit none
  !f2py integer*1, dimension(1000) :: cyl
  type(cylindre_2),intent(in) :: cyl
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
! local
  !f2py integer*1, dimension(1000) :: cercle
  type(GCercle) :: cercle
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
  !f2py integer*1, dimension(1000) :: P1,P2,P3,P4,P5,P6
  type(point) :: P1,P2,P3,P4,P5,P6
  integer :: iface
  real(rp) :: r
  real(rp) :: vmin,vmax
  integer :: ierror
  
  np = 0
  ierror = 0  

  iface = cyl%index
  rep = cyl%repere
  r = cyl%rayon
  vmin = cyl%vmin
  vmax = cyl%vmax
  
  if(abs(rep%e1(3)).lt.Epsilon)then
    P1 = [-r,0._rp,vmin]
    P2 = [-r,0._rp,vmax]
    P3 = [r,0._rp,vmin]
    P4 = [r,0._rp,vmax]
  elseif(rep%e1(3) < 0.)then
    P1 = [r,0._rp,vmin]
    P2 = [r,0._rp,vmax]
    P3 = [-r,0._rp,vmin]
    P4 = [-r,0._rp,vmax]
  elseif(rep%e1(3) > 0.)then
    P1 = [-r,0._rp,vmin]
    P2 = [-r,0._rp,vmax]
    P3 = [r,0._rp,vmin]
    P4 = [r,0._rp,vmax]
  endif
  call loc2cart(P1%coord,rep,P1%coord)
  call loc2cart(P2%coord,rep,P2%coord)
  call loc2cart(P3%coord,rep,P3%coord)
  call loc2cart(P4%coord,rep,P4%coord) 
  
  P1%nface = 1
  P1%face(1) = iface
  P2%nface = 1
  P2%face(1) = iface
  P3%nface = 1
  P3%face(1) = iface
  P4%nface = 1
  P4%face(1) = iface
  
  call create_tab_intersect(P1,P2,t,dx,tab_point,np,iface)
  
  if(np==0)then
    call create_tab_intersect(P3,P4,t,dx,tab_point,np,iface)
  endif   
  
  if(np==0)then
    if(abs(rep%e1(3)).lt.Epsilon)then
      P5 = [0._rp,-r,vmin]
      P6 = [0._rp,r,vmin]
    elseif(rep%e1(3) > 0.)then
      P5 = [0._rp,r,vmin]
      P6 = [0._rp,-r,vmin]
    elseif(rep%e1(3) > 0.)then
      P5 = [0._rp,r,vmin]
      P6 = [0._rp,-r,vmin]
    endif
    call loc2cart(P5%coord,rep,P5%coord)
    call loc2cart(P6%coord,rep,P6%coord)
    call GCercle_from_point(P1%coord,P3%coord,P5%coord,cercle,ierror)
    call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
    call GCercle_from_point(P1%coord,P3%coord,P6%coord,cercle,ierror)
    call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
  endif
  
  if(np==0)then
    if(abs(rep%e1(3)).lt.Epsilon)then
      P5 = [0._rp,-r,vmax]
      P6 = [0._rp,r,vmax]
    elseif(rep%e1(3) > 0.)then
      P5 = [0._rp,r,vmax]
      P6 = [0._rp,-r,vmax]
    elseif(rep%e1(3) > 0.)then
      P5 = [0._rp,r,vmax]
      P6 = [0._rp,-r,vmax]
    endif
    call loc2cart(P5%coord,rep,P5%coord)
    call loc2cart(P6%coord,rep,P6%coord)
    call GCercle_from_point(P2%coord,P4%coord,P5%coord,cercle,ierror)
    call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
    call GCercle_from_point(P4%coord,P2%coord,P6%coord,cercle,ierror)
    call create_tabcurve_intersect(P4,P2,t,cercle,dx,tab_point,np,iface)
  endif
  
end subroutine point_inter_cyl

subroutine point_inter_cone(geom,t,dx,tab_point,np)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cone),intent(in) :: geom
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
! local
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
  !f2py integer*1, dimension(1000) :: P1,P2,P3,P4
  type(point) :: P1,P2,P3,P4
  integer :: iface
  real(rp) :: r1,r2
  real(rp) :: vmin,vmax
  integer :: ierror
  
  np = 0
  ierror = 0  

  iface = geom%index
  rep = geom%repere
  r1 = geom%r1
  r2 = geom%r2
  vmin = geom%vmin
  vmax = geom%vmax
  
  if(abs(rep%e1(3)).lt.Epsilon)then
    P1 = [-r1,0._rp,vmin]
    P2 = [-r2,0._rp,vmax]
    P3 = [r1,0._rp,vmin]
    P4 = [r2,0._rp,vmax]
  elseif(rep%e1(3) < 0.)then
    P1 = [r1,0._rp,vmin]
    P2 = [r2,0._rp,vmax]
    P3 = [-r1,0._rp,vmin]
    P4 = [-r2,0._rp,vmax]
  elseif(rep%e1(3) > 0.)then
    P1 = [-r1,0._rp,vmin]
    P2 = [-r2,0._rp,vmax]
    P3 = [r1,0._rp,vmin]
    P4 = [r2,0._rp,vmax]
  endif
  call loc2cart(P1%coord,rep, P1%coord)
  call loc2cart(P2%coord,rep,P2%coord)
  call loc2cart(P3%coord,rep,P3%coord)
  call loc2cart(P4%coord,rep,P4%coord) 
  
  P1%nface = 1
  P1%face(1) = iface
  P2%nface = 1
  P2%face(1) = iface
  P3%nface = 1
  P3%face(1) = iface
  P4%nface = 1
  P4%face(1) = iface
  
  call create_tab_intersect(P1,P2,t,dx,tab_point,np,iface)
  
  if(np==0)then
    call create_tab_intersect(P3,P4,t,dx,tab_point,np,iface)
  endif   
  
! ## FIXME : calcul separe
!  if(np==0)then
!    if(abs(rep%e1(3)).lt.Epsilon)then
!      P5 = [0._rp,-r1,vmin]
!      P6 = [0._rp,r1,vmin]
!    elseif(rep%e1(3) > 0.)then
!      P5 = [0._rp,r1,vmin]
!      P6 = [0._rp,-r1,vmin]
!    elseif(rep%e1(3) > 0.)then
!      P5 = [0._rp,r1,vmin]
!      P6 = [0._rp,-r1,vmin]
!    endif
!    P5%coord = loc2cart(P5%coord,rep)
!    P6%coord = loc2cart(P6%coord,rep)
!    call GCercle_from_point(P1%coord,P3%coord,P5%coord,cercle,ierror)
!    call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
!    call GCercle_from_point(P1%coord,P3%coord,P6%coord,cercle,ierror)
!    call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
!  endif
!  
!  if(np==0)then
!    if(abs(rep%e1(3)).lt.Epsilon)then
!      P5 = [0._rp,-r2,vmax]
!      P6 = [0._rp,r2,vmax]
!    elseif(rep%e1(3) > 0.)then
!      P5 = [0._rp,r2,vmax]
!      P6 = [0._rp,-r2,vmax]
!    elseif(rep%e1(3) > 0.)then
!      P5 = [0._rp,r2,vmax]
!      P6 = [0._rp,-r2,vmax]
!    endif
!    P5%coord = loc2cart(P5%coord,rep)
!    P6%coord = loc2cart(P6%coord,rep)
!    call GCercle_from_point(P2%coord,P4%coord,P5%coord,cercle,ierror)
!    call create_tabcurve_intersect(P2,P4,t,cercle,dx,tab_point,np,iface)
!    call GCercle_from_point(P4%coord,P2%coord,P6%coord,cercle,ierror)
!    call create_tabcurve_intersect(P4,P2,t,cercle,dx,tab_point,np,iface)
!  endif
  
end subroutine point_inter_cone

subroutine point_inter_disque(disque,t,dx,tab_point,np)
  implicit none
  !f2py integer*1, dimension(1000) :: disque
  type(disque2),intent(in) :: disque
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
! local
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
  !f2py integer*1, dimension(1000) :: P1,P2,P3,P4
  type(point) :: P1,P2,P3,P4
  integer :: iface
  real(rp) :: r,r2max
  real(rp) :: vmin,vmax
  
  iface = disque%index
  rep = disque%repere
  r2max = disque%r2max
  r = sqrt(r2max)
  vmin = 0._RP
  vmax = 0._RP
  
  if(abs(rep%e1(3)).lt.Epsilon)then
    P1 = [-r,0._rp,vmin]
    P2 = [-r,0._rp,vmax]
    P3 = [r,0._rp,vmin]
    P4 = [r,0._rp,vmax]
  elseif(rep%e1(3) < 0.)then
    P1 = [r,0._rp,vmin]
    P2 = [r,0._rp,vmax]
    P3 = [-r,0._rp,vmin]
    P4 = [-r,0._rp,vmax]
  elseif(rep%e1(3) > 0.)then
    P1 = [-r,0._rp,vmin]
    P2 = [-r,0._rp,vmax]
    P3 = [r,0._rp,vmin]
    P4 = [r,0._rp,vmax]
  endif
  call loc2cart(P1%coord,rep,P1%coord)
  call loc2cart(P2%coord,rep,P2%coord)
  call loc2cart(P3%coord,rep,P3%coord)
  call loc2cart(P4%coord,rep,P4%coord) 
  
  P1%nface = 1
  P1%face(1) = iface
  P2%nface = 1
  P2%face(1) = iface
  P3%nface = 1
  P3%face(1) = iface
  P4%nface = 1
  P4%face(1) = iface
  
  call create_tab_intersect(P1,P3,t,dx,tab_point,np,iface) 
  
end subroutine point_inter_disque  
  
subroutine point_inter_plan(plan,t,dx,tab_point,np)  
  implicit none
  !f2py integer*1, dimension(1000) :: plan
  type(Gplan),intent(in) :: plan
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
! local
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d) :: rep
  !f2py integer*1, dimension(1000) :: P1,P2,P3,P4
  type(point) :: P1,P2,P3,P4
  integer :: iface
  real(rp) :: umin,umax
  real(rp) :: vmin,vmax
  !f2py integer*1, dimension(1000) :: tab_point2
  type(point),dimension(50) :: tab_point2
  integer :: np2,j
  
  np = 0
  
  iface = plan%index
  rep = plan%repere
  umin = plan%umin
  umax = plan%umax
  vmin = plan%vmin
  vmax = plan%vmax
  
  P1 = [umin,vmin,0._rp]
  P2 = [umax,vmin,0._rp]
  P3 = [umax,vmax,0._rp]
  P4 = [umin,vmax,0._rp]
  
  call loc2cart(P1%coord,rep,P1%coord)
  call loc2cart(P2%coord,rep,P2%coord)
  call loc2cart(P3%coord,rep,P3%coord)
  call loc2cart(P4%coord,rep,P4%coord)
  
  P1%nface = 1
  P1%face(1) = iface
  P2%nface = 1
  P2%face(1) = iface
  P3%nface = 1
  P3%face(1) = iface
  P4%nface = 1
  P4%face(1) = iface
  
  call create_tab_intersect(P1,P2,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P2,P3,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P3,P4,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P4,P1,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif    
!
end subroutine point_inter_plan

subroutine point_inter_Wigley(Wigley,t,dx,tab_point,np)  
  
    !f2py integer*1, dimension(1000)    :: Wigley
  type(Twigley),intent(in) :: Wigley
  real(rp),intent(in) :: t,dx
  !f2py integer*1, dimension(1000)    :: tab_point
  type(point),dimension(*),intent(inout) :: tab_point
  integer,intent(inout) :: np
  
  !f2py integer*1, dimension(1000)    :: rep
  type(repere3d) :: rep
  !f2py integer*1, dimension(1000)    :: P1,P2,P3,P4
  type(point) :: P1,P2,P3,P4
  integer :: iface
  real(rp) :: umin,umax
  real(rp) :: vmin,vmax
  !f2py integer*1, dimension(1000)    :: tab_point2
  type(point),dimension(50) :: tab_point2
  integer :: np2,j
  
  np = 0
  
  iface = Wigley%index
  rep = Wigley%repere
  umin = Wigley%umin
  umax = Wigley%umax
  vmin = Wigley%vmin
  vmax = Wigley%vmax
  
  P1 = [umin,vmin,0._rp]
  P2 = [umax,vmin,0._rp]
  P3 = [umax,vmax,0._rp]
  P4 = [umin,vmax,0._rp]
  
  call loc2cart(P1%coord,rep,P1%coord)
  call loc2cart(P2%coord,rep,P2%coord)
  call loc2cart(P3%coord,rep,P3%coord)
  call loc2cart(P4%coord,rep,P4%coord)
  
  P1%nface = 1
  P1%face(1) = iface
  P2%nface = 1
  P2%face(1) = iface
  P3%nface = 1
  P3%face(1) = iface
  P4%nface = 1
  P4%face(1) = iface
  
  call create_tab_intersect(P1,P2,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P2,P3,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P3,P4,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif
  
  call create_tab_intersect(P4,P1,t,dx,tab_point2,np2,iface)
  if(np2>0)then
    do j=1,np2
      tab_point(np+j) = tab_point2(j)
    enddo
    np=np+np2
  endif    
!
end subroutine point_inter_Wigley

subroutine point_inter_sph(geom,t,dx,tab_point,np)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(sphere),intent(in) :: geom
    real(rp),intent(in) :: t,dx
    !f2py integer*1, dimension(1000) :: tab_point
    type(point),dimension(*),intent(inout) :: tab_point
    integer,intent(inout) :: np
!   local
    !f2py integer*1, dimension(1000) :: cercle
    type(GCercle) :: cercle
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d) :: rep
    !f2py integer*1, dimension(1000) :: P1,P2,P3,P4
    type(point) :: P1,P2,P3,P4
    integer :: iface
    real(rp) :: r
    real(rp) :: vmin,vmax
    real(rp),dimension(3) :: C
    integer :: ierror
    logical :: is_truncated
    
    np = 0
    ierror = 0
    
    iface = geom%index
    rep = geom%repere
    C = rep%origine
    r = geom%radius
    vmin = geom%vmin
    vmax = geom%vmax
    
    call Computation_eq_sph3D(0._rp,vmin,geom,P1%coord) ! Previously only P1
    call Computation_eq_sph3D(pi,vmin,geom,P2%coord) ! Previously only P2
    call Computation_eq_sph3D(0._rp,vmax,geom,P3%coord) ! Previously only P3
    if(vmax.lt.0.5_rp*pi-Epsilon)then
        call computation_eq_sph3D(pi,vmax,geom,P4%coord) ! Previously only P4
        is_truncated = .true.
    else
        P4 = [-999._rp,-999._rp,-999._rp]
        is_truncated = .false.
    endif
    
    P1%nface = 1
    P1%face(1) = iface
    P2%nface = 1
    P2%face(1) = iface
    P3%nface = 1
    P3%face(1) = iface
    P4%nface = 1
    P4%face(1) = iface
    
    if(is_truncated)then
        call GCercle_from_2PointCenter(P1%coord,P3%coord,C,cercle,ierror)
        call create_tabcurve_intersect(P1,P3,t,cercle,dx,tab_point,np,iface)
        call GCercle_from_2PointCenter(P2%coord,P4%coord,C,cercle,ierror)
        call create_tabcurve_intersect(P2,P4,t,cercle,dx,tab_point,np,iface)
    else    
        call GCercle_from_2PointCenter(P1%coord,P2%coord,C,cercle,ierror)
        call create_tabcurve_intersect(P1,P2,t,cercle,dx,tab_point,np,iface)
    endif
    
    if(np==0)then
        ! FIXME : complete for non-symetric test case
    endif    
    
end subroutine point_inter_sph

subroutine intersection_3D_cylindre(cyl,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
  
    !f2py integer*1, dimension(1000) :: cyl
    type(cylindre_2),intent(in) :: cyl
    !f2py integer*1, dimension(1000) :: tab
    type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
    integer,intent(inout) :: n_tab
    !f2py integer*1, dimension(1000) :: list_point
    type(point),dimension(:),intent(in) :: list_point
    integer,intent(in) :: nlist
    real(rp),intent(in) :: t
    real(rp),intent(in) :: dx
    integer,intent(inout) :: nflag
    integer,intent(inout) :: ierror
    
    !f2py integer*1, dimension(1000) :: P0,Pn,Pnp1,Q,P
    type(point) :: P0,Pn,Pnp1,Q,P
    !f2py integer*1, dimension(1000) :: tab_point
    type(point),dimension(100) :: tab_point
    integer :: j,k
    integer :: np,iflag,iface,n_bord
    integer :: isens
    real(rp) :: dl
    real(rp) :: h,hh,h6
    !f2py integer*1, dimension(1000) :: V,V0
    type(vector) :: V,V0
    logical :: is_loop,is_fin,is_border,is_border0,is_pres
    integer,parameter :: n=2
    real(rp),dimension(n) :: y,yl,yout
    real(rp),dimension(2) :: diff1,diff2,diff3
    integer,parameter :: jmax =10000
    integer,dimension(10) :: aux
    integer :: naux
    !f2py integer*3, dimension(3000) :: ptr1,ptr2,ptr
    type(chaine_point),pointer :: ptr1,ptr2,ptr
    !f2py integer*1, dimension(1000) :: tab2
    type(chaine_point_pt),dimension(100,2) :: tab2
    integer :: n_tab2
    real(rp),dimension(:,:),allocatable :: tbound
    integer :: nbound
    
    ! This subroutine computes the intersection of a cylinder with the free surface.
    
    np = 0 
    n_tab2 = 0 
    !ds = 0.01
    
    if(Symmetry)then
        allocate(tbound(2,2))
        tbound(1,1:2) = [0._rp,PI]
        tbound(2,1:2) = [cyl%vmin,cyl%vmax]
        nbound = 2
    else
        allocate(tbound(1,2))
        tbound(1,1:2) = [cyl%vmin,cyl%vmax]
        nbound = 1
    endif
    
    iface = cyl%index
    ! Recherche point appartenant a la surface  
    do j=1,nlist
        P = list_point(j)
        call common_int(P%face,P%nface,[iface],1,aux,naux,ierror)
        if(naux==1)then
            np=np+1
            tab_point(np) = P
        endif
    enddo
    
    do j=1,n_tab
        do k=1,2
            P = tab(j,k)%pt%val
            call common_int(P%face,P%nface,[iface],1,aux,naux,ierror)
            if(naux==1)then
            np = np+1
            tab_point(np)=P
            endif
        enddo
    enddo
    
    iflag = nflag
    
    ! Debut boucle calcul courbes d'intersection  
    do k=1,np
      
    nullify(ptr1)
    nullify(ptr2)
    
    P0 = tab_point(k)   
        
    do j=1,n_tab2
        ptr => tab2(j,1)%pt
        call point_in_chaine(ptr,P0,is_pres,0.4*dx,P,ierror)
        if(is_pres)then
            goto 100    ! Test point deja present => fin boucle k
        endif
    enddo
    
    n_tab2 = n_tab2+1
    iflag = iflag+1
    P0%flag = iflag
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
       
    j=1
    dl = 0._rp
    n_bord = 0
    isens = 1
        
    Pn = P0
    call cart2cyl(Pn,cyl,y)
        
    call check_border_ND(y,tbound,nbound,y,yout,is_border0,ierror)
    if(is_border0)then
        ptr1%val%bf = 1
        if(Symmetry .and. abs(ptr1%val%coord(2)).lt.Epsilon)then
            ptr1%val%bf = 2
        endif
        n_bord = 1
    endif
    
    do while (.not.is_fin)
      
        h = dble(isens)*dsi
        hh=0.5*h
        h6=h/6.
         
        ! RK4 1ere passe
        call derive_cylindre(Pn,cyl,t,diff1)
        do j=1,n
        yl(j)=y(j)+hh*diff1(j)
        enddo
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)

  
        ! RK4 2eme passe
        if (.not.is_border) then
            call Computation_eq_cyl3D(yl(1),yl(2),cyl,Q%coord) ! Previously only Q
            call derive_cylindre(Q,cyl,t,diff2)
            do j=1,n
                yl(j)=y(j)+hh*diff2(j)
            enddo
            call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
        endif  

        ! RK4 3eme passe  
        if(.not.is_border) then
            call Computation_eq_cyl3D(yl(1),yl(2),cyl,Q%coord) ! Previously only Q
            call derive_cylindre(Q,cyl,t,diff3)
            do j=1,n
                yl(j)=y(j)+h*diff3(j)
                diff3(j) = diff2(j)+diff3(j)
            enddo
            call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
        endif
  
        ! RK4 4eme passe 
        if(.not.is_border)then
            call Computation_eq_cyl3D(yl(1),yl(2),cyl,Q%coord)
            call derive_cylindre(Q,cyl,t,diff2)
            do j=1,n
                yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
            enddo
            call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
        endif 

        if(.not.is_border)then
            yout = yl
        endif
  
        call Computation_eq_cyl3D(yout(1),yout(2),cyl,Pnp1%coord)
        
        ! Norm between the previous point and the current point.
        call assign_vector_coord(V,Pnp1%coord - Pn%coord)
        
        ! Norm between P0 and the current point to initialized is_loop.
        call assign_vector_coord(V0,Pnp1%coord - P0%coord)
        
        ! Ajout point intersection     
        dl = dl+V%length
        if(dl >= dx .and. V0%length.ge.dx .and. .not.is_border)then ! Point is not too close from the the previsous one (dl>=dx) and P0 (V0%length.ge.dx) and it is not on the border.
            Pnp1%nface = 2
            Pnp1%face(1:2) = [HouleRF%index,iface]
            Pnp1%flag = iflag
            Pnp1%bf = 0
            Pnp1%nedge = 0
            Pnp1%edge = 0                
            if(isens > 0)then
                call add_point_suiv(ptr2,Pnp1)
                ptr2 => ptr2%suiv
            else
                call add_point_prec(ptr1,Pnp1)
                ptr1 => ptr1%prec
            endif
            dl = 0._rp
        endif  
        
        ! Ajout point frontiere      
        if(is_border)then
            Pnp1%nface=2
            Pnp1%face(1:2) = [HouleRF%index,iface]
            Pnp1%flag = iflag
            
            if(abs(Pnp1%coord(2)).lt.Epsilon .and. Symmetry)then
                Pnp1%bf = 2
                Pnp1%nface = 3
                Pnp1%face(3) = -1
            else
                Pnp1%bf = 1        
            endif
            if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
                if(isens > 0)then
                    ptr2%val=Pnp1
                else
                    ptr1%val=Pnp1
                endif
            n_bord = n_bord+1
            elseif(abs(dl).gt.Epsilon)then
                if(isens > 0)then
                    call add_point_suiv(ptr2,Pnp1)
                    ptr2 => ptr2%suiv
                else
                    call add_point_prec(ptr1,Pnp1)
                    ptr1 => ptr1%prec
                endif
                n_bord = n_bord+1
            endif
            ! Mise a jour isens 
            if(isens == 1 .and. n_bord < 2)then
                isens = -1
                Pnp1 = P0
                call cart2cyl(P0,cyl,yout)
            else
                is_fin = .true.
            endif
        endif
        
        ! Ajout point pour fermer la boucle 
        is_loop = V0%length < 0.9*dsi .and. j>1 .and. isens==1     
        if(is_loop)then
            if(dl < 0.5*dx)then
                if(isens > 0)then
                ptr2 => ptr2%prec
                deallocate(ptr2%suiv)
                ptr2%suiv => ptr1
                ptr1%prec => ptr2
                else
                ptr1 => ptr1%suiv
                deallocate(ptr1%prec)
                ptr1%prec => ptr2
                ptr2%suiv => ptr1
                endif
            else
                ptr2%suiv => ptr1
                ptr1%prec => ptr2
            endif
        endif
      
        is_fin = is_fin .or. is_loop
      
        j=j+1
        Pn = Pnp1
        y = yout
     
    enddo   
    tab2(n_tab2,1)%pt => ptr1
    tab2(n_tab2,2)%pt => ptr2
    
    100 continue
    
    enddo 
  
    nflag = iflag   
  
    do j=1,n_tab2
        do k=1,2
            tab(n_tab+j,k)%pt => tab2(j,k)%pt
        enddo
    enddo  
    n_tab = n_tab+n_tab2
      
    if(allocated(tbound)) deallocate(tbound)
  
end subroutine intersection_3D_cylindre

subroutine intersection_3D_disque(disque,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: disque
  type(disque2),intent(in) :: disque
  !f2py integer*1, dimension(1000) :: tab
  type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
  integer,intent(inout) :: n_tab
  !f2py integer*1, dimension(1000) :: list_point
  type(point),dimension(:),intent(in) :: list_point
  integer,intent(in) :: nlist
  real(rp),intent(in) :: t
  real(rp),intent(in) :: dx
  integer,intent(inout) :: nflag
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: P0,Pn,Pnp1,Q,P
  type(point) :: P0,Pn,Pnp1,Q,P
  !f2py integer*1, dimension(1000) :: tab_point,tab_point2
  type(point),dimension(100) :: tab_point,tab_point2
  integer :: j,jj,k
  integer :: np,np2,iflag,iface,n_bord
  integer :: isens
  real(rp) :: ds,dl
  real(rp) :: h,hh,h6
  !f2py integer*1, dimension(1000) :: V,V0
  type(vector) :: V,V0
  logical :: is_loop,is_fin,is_border,is_border0,is_pres
  integer,parameter :: n=2
  real(rp),dimension(n) :: y,yl,yout
  real(rp),dimension(3) :: x
  real(rp),dimension(2) :: diff1,diff2,diff3 
  real(rp) :: r2,r2max
  integer,parameter :: jmax = 5000 
  integer,dimension(10) :: aux
  integer :: naux
  !f2py integer*1, dimension(1000) :: ptr1,ptr2,ptr
  type(chaine_point),pointer :: ptr1,ptr2,ptr
  !f2py integer*1, dimension(1000) :: tab2
  type(chaine_point_pt),dimension(100,2) :: tab2
  integer :: n_tab2
    
  np = 0  
  np2 = 0
  n_tab2 = 0
  ds = 0.01

  iface = disque%index
  r2max = disque%r2max
! Recherche point appartenant \E0 la surface  
  do j=1,nlist
    P = list_point(j)
    call common_int(P%face,P%nface,[iface],1,aux,naux)
    if(naux==1)then
        np=np+1
        tab_point(np) = P
    endif
  enddo

  do j=1,n_tab
    do k=1,2
      P = tab(j,k)%pt%val
      x = P%coord
      call cart2loc(x,disque%repere,x)
      r2 = x(1)*x(1)+x(2)*x(2)
      if(abs(x(3)).lt.0.001 .and. abs(r2-r2max).lt.0.001)then
        np=np+1
        tab_point(np) = P
      endif  
    enddo
  enddo
! Recherche nouveaux points d'intersection    
  call point_inter(disque,t,dx,tab_point2,np2)
  
  do j=1,np2
    tab_point(np+j) = tab_point2(j)
  enddo
  np=np+np2
  
  iflag = nflag
  
  do k=1,np
  
    nullify(ptr1)
    nullify(ptr2)
  
    P0 = tab_point(k)
    do j=1,n_tab2
      ptr => tab2(j,1)%pt
      call point_in_chaine(ptr,P0,is_pres,0.6*dx,P,ierror)
      if(is_pres)then
        goto 100
      endif
    enddo
    
    n_tab2 = n_tab2+1
    
    iflag = iflag+1
    P0%flag = iflag
    P0%nface = 2
    P0%face(1:2) = [HouleRF%index,iface]
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
    
    isens = 1   
    jj=1
    dl = 0._rp
    n_bord = 0
    
    Pn = P0
    call cart2loc(Pn%coord,disque%repere,x)
    y = x(1:2)
    
    call check_border_disque(y,disque,y,yout,is_border0,ierror)
    if(is_border0)then
      ptr1%val%bf=1
      n_bord = 1
    endif
    
    do while (.not.is_fin)
      h = dble(isens)*dsi
      hh=0.5*h
      h6=h/6.
    
! 1ere passe
      call derive_disque(Pn,disque,t,diff1)
      do j=1,n
        yl(j)=y(j)+hh*diff1(j)
      enddo
      call check_border_disque(yl,disque,y,yout,is_border,ierror)
  
! 2eme passe
      if (.not.is_border) then
        call loc2cart([yl(1),yl(2),0._rp],disque%repere,Q%coord) ! Previously only Q
        call derive_disque(Q,disque,t,diff2)
        do j=1,n
          yl(j)=y(j)+hh*diff2(j)
        enddo
        call check_border_disque(yl,disque,y,yout,is_border,ierror)
      endif  

! 3eme passe  
      if(.not.is_border) then
        call loc2cart([yl(1),yl(2),0._rp],disque%repere,Q%coord) ! Previously only Q
        call derive_disque(Q,disque,t,diff3)
        do j=1,n
          yl(j)=y(j)+h*diff3(j)
          diff3(j) = diff2(j)+diff3(j)
        enddo
        call check_border_disque(yl,disque,y,yout,is_border,ierror)
      endif
  
! 4eme passe  
      if(.not.is_border)then
        call loc2cart([yl(1),yl(2),0._rp],disque%repere,Q%coord) ! Previously only Q
        call derive_disque(Q,disque,t,diff2)
        do j=1,n
          yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
        enddo
        call check_border_disque(yl,disque,y,yout,is_border,ierror)
      endif 

      if(.not.is_border)then
        yout = yl
      endif
      
      call loc2cart([yout(1),yout(2),0._rp],disque%repere,Pnp1%coord)
      
      call assign_vector_coord(V,Pnp1%coord - Pn%coord)
      call assign_vector_coord(V0,Pnp1%coord - P0%coord)
      
      dl = dl+V%length
      if(dl >= dx .and. .not.is_border)then
       !n_tab = n_tab+1
        Pnp1%nface = 2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 0
        Pnp1%nedge = 0
        Pnp1%edge = 0
        
        if(isens > 0)then
          call add_point_suiv(ptr2,Pnp1)
          ptr2 => ptr2%suiv
        else
          call add_point_prec(ptr1,Pnp1)
          ptr1 => ptr1%prec
        endif
        dl = 0._rp
      endif
      
      if(is_border)then
        Pnp1%nface=2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 1
        
        if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
          if(isens > 0)then
            ptr2%val=Pnp1
          else
            ptr1%val=Pnp1
          endif
          n_bord = n_bord+1
        elseif(abs(dl).gt.Epsilon)then
          if(isens > 0)then
            call add_point_suiv(ptr2,Pnp1)
            ptr2 => ptr2%suiv
          else
            call add_point_prec(ptr1,Pnp1)
            ptr1 => ptr1%prec
          endif
          n_bord = n_bord+1
        endif
        
!       Mise a jour isens      
        if(isens == 1 .and. n_bord < 2)then
          isens = -1
          Pnp1 = P0
          call cart2loc(P0%coord,disque%repere,x)
          yout = x(1:2)
        else
          is_fin = .true.
        endif
      endif
      
      is_loop = V0%length < 0.9*dsi .and. jj>1 .and. isens==1 
      if(is_loop)then
        if(dl < 0.5*dx)then
          if(isens > 0)then
            ptr2 => ptr2%prec
            deallocate(ptr2%suiv)
            ptr2%suiv => ptr1
            ptr1%prec => ptr2
          else
            ptr1 => ptr1%suiv
            deallocate(ptr1%prec)
            ptr1%prec => ptr2
            ptr2%suiv => ptr1
          endif
        else
          ptr1%prec => ptr2
          ptr2%suiv => ptr1
        endif
      endif
      
      is_fin = is_fin .or. is_loop
      
      jj=jj+1
      Pn = Pnp1
      y = yout
     
    enddo    
        
     tab2(n_tab2,1)%pt => ptr1
     tab2(n_tab2,2)%pt => ptr2
!
100 continue
!
  enddo    
  
  nflag = iflag
  
  do j=1,n_tab2
    do k=1,2
      tab(n_tab+j,k)%pt => tab2(j,k)%pt
    enddo
  enddo
  n_tab = n_tab + n_tab2
!  
end subroutine intersection_3D_disque

subroutine intersection_3D_plan(plan,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: plan
  type(Gplan),intent(in) :: plan
  !f2py integer*1, dimension(1000) :: tab
  type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
  integer,intent(inout) :: n_tab
  !f2py integer*1, dimension(1000) :: list_point
  type(point),dimension(:),intent(in) :: list_point
  integer,intent(in) :: nlist
  real(rp),intent(in) :: t
  real(rp),intent(in) :: dx
  integer,intent(inout) :: nflag
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: P0,Pn,Pnp1,Q,P
  type(point) :: P0,Pn,Pnp1,Q,P
  !f2py integer*1, dimension(1000) :: tab_point,tab_point2
  type(point),dimension(100) :: tab_point,tab_point2
  integer :: j,k
  integer :: np,np2,iflag,iface,n_bord
  integer :: isens
  real(rp) :: ds,dl
  real(rp) :: h,hh,h6
  !f2py integer*1, dimension(1000) :: V,V0
  type(vector) :: V,V0
  logical :: is_loop,is_fin,is_border,is_border0,is_pres
  integer,parameter :: n=2
  real(rp),dimension(n) :: y,yl,yout
  real(rp),dimension(3) :: x
  real(rp),dimension(2) :: diff1,diff2,diff3 
  real(rp) :: umin,umax,vmin,vmax
  integer,parameter :: jmax = 5000 
  integer,dimension(10) :: aux
  integer :: naux
  !f2py integer*1, dimension(1000) :: ptr1,ptr2,ptr
  type(chaine_point),pointer :: ptr1,ptr2,ptr
  !f2py integer*1, dimension(1000) :: tab2
  type(chaine_point_pt),dimension(100,2) :: tab2
  integer :: n_tab2
  
  np = 0  
  np2 = 0
  n_tab2 = 0
  ds = 0.01

  iface = plan%index
  umin = plan%umin
  umax = plan%umax
  vmin = plan%vmin
  vmax = plan%vmax
!
! Recherche point appartenant \E0 la surface  
  do j=1,nlist
    P = list_point(j)
    call common_int(P%face,P%nface,[iface],1,aux,naux)
    if(naux==1)then
        np=np+1
        tab_point(np) = P
    endif
  enddo

  do j=1,n_tab
    do k=1,2
      P = tab(j,k)%pt%val
      x = P%coord
      call cart2loc(x,plan%repere,x)
      if((abs(x(1)-umin).lt.0.0001d0 .or. abs(x(1)-umax).lt.0.0001d0 .or.&
      &  abs(x(2)-vmin).lt.0.0001d0 .or. abs(x(2)-vmax).lt.0.0001d0).and.&
      &  abs(x(3)).lt.0.0001d0)then
        np=np+1
        tab_point(np) = P
      endif  
    enddo
  enddo
! Recherche nouveaux points d'intersection    
  call point_inter(plan,t,dx,tab_point2,np2)
  
  do j=1,np2
    tab_point(np+j) = tab_point2(j)
  enddo
  np=np+np2
  
  iflag = nflag
  
  do k=1,np
  
    nullify(ptr1)
    nullify(ptr2)
  
    P0 = tab_point(k)
    do j=1,n_tab2
      ptr => tab2(j,1)%pt
      call point_in_chaine(ptr,P0,is_pres,0.6*dx,P,ierror)
      if(is_pres)then
        goto 100
      endif
    enddo
    
    n_tab2 = n_tab2+1
    iflag = iflag+1
    P0%flag = iflag
    P0%nface = 2
    P0%face(1:2) = [HouleRF%index,iface]
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
    
    isens = 1   
    j=1
    dl = 0._rp
    n_bord = 0
    
    Pn = P0
    call cart2loc(Pn%coord,plan%repere,x)
    y = x(1:2)
    
    call check_border_plan(y,plan,y,yout,is_border0,ierror)
    if(is_border0)then
      ptr1%val%bf=1
      n_bord = 1
    endif
    
    do while (.not.is_fin)
      
      h = dble(isens)*dsi
      hh=0.5*h
      h6=h/6.
    
! 1ere passe
      call derive_plan(Pn,plan,t,diff1)
      do j=1,n
        yl(j)=y(j)+hh*diff1(j)
      enddo
      call check_border_plan(yl,plan,y,yout,is_border,ierror)
  
! 2eme passe
      if (.not.is_border) then
        call loc2cart([yl(1),yl(2),0._rp],plan%repere,Q%coord) ! Previously only Q
        call derive_plan(Q,plan,t,diff2)
        do j=1,n
          yl(j)=y(j)+hh*diff2(j)
        enddo
        call check_border_plan(yl,plan,y,yout,is_border,ierror)
      endif  

! 3eme passe  
      if(.not.is_border) then
        call loc2cart([yl(1),yl(2),0._rp],plan%repere,Q%coord) ! Previously only Q
        call derive_plan(Q,plan,t,diff3)
        do j=1,n
          yl(j)=y(j)+h*diff3(j)
          diff3(j) = diff2(j)+diff3(j)
        enddo
        call check_border_plan(yl,plan,y,yout,is_border,ierror)
      endif
  
! 4eme passe  
      if(.not.is_border)then
        call loc2cart([yl(1),yl(2),0._rp],plan%repere,Q%coord) ! Previously only Q
        call derive_plan(Q,plan,t,diff2)
        do j=1,n
          yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
        enddo
        call check_border_plan(yl,plan,y,yout,is_border,ierror)
      endif 

      if(.not.is_border)then
        yout = yl
      endif
      
      call loc2cart([yout(1),yout(2),0._rp],plan%repere,Pnp1%coord)
       
      call assign_vector_coord(V,Pnp1%coord - Pn%coord)
      call assign_vector_coord(V0,Pnp1%coord - P0%coord)
     
      dl = dl+V%length
      if(dl >= dx .and. .not.is_border)then
       !n_tab = n_tab+1
        Pnp1%nface = 2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 0
        if(isens > 0)then
          call add_point_suiv(ptr2,Pnp1)
          ptr2 => ptr2%suiv
        else
          call add_point_prec(ptr1,Pnp1)
          ptr1 => ptr1%prec
        endif
        dl = 0._rp
      endif       
      
      if(is_border)then
        Pnp1%nface=2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 1        
        if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
          if(isens > 0)then
            ptr2%val=Pnp1
          else
            ptr1%val=Pnp1
          endif
          n_bord = n_bord+1
        elseif(abs(dl).gt.Epsilon)then
          if(isens > 0)then
            call add_point_suiv(ptr2,Pnp1)
            ptr2 => ptr2%suiv
          else
            call add_point_prec(ptr1,Pnp1)
            ptr1 => ptr1%prec
          endif
          n_bord = n_bord+1
        endif
!       Mise \E0 jour isens      
        if(isens == 1 .and. n_bord < 2)then
          isens = -1
          Pnp1 = P0
          call cart2loc(P0%coord,plan%repere,x)
          yout = x(1:2)
        else
          is_fin = .true.
        endif
      endif
      
      is_loop = V0%length < 0.9*dsi .and. j>1 .and. isens==1 
      if(is_loop)then
        if(dl < 0.5*dx)then
          if(isens > 0)then
            ptr2 => ptr2%prec
            deallocate(ptr2%suiv)
            ptr2%suiv => ptr1
            ptr1%prec => ptr2
          else
            ptr1 => ptr1%suiv
            deallocate(ptr1%prec)
            ptr1%prec => ptr2
            ptr2%suiv => ptr1
          endif
        else
          ptr1%prec => ptr2
          ptr2%suiv => ptr1
        endif
      endif
      
      is_fin = is_fin .or. is_loop
      
      j=j+1
      Pn = Pnp1
      y = yout
     
    enddo        
!
     tab2(n_tab2,1)%pt => ptr1
     tab2(n_tab2,2)%pt => ptr2
!
100 continue
!
  enddo    
  
  nflag = iflag
  
  do j=1,n_tab2
    do k=1,2
      tab(n_tab+j,k)%pt => tab2(j,k)%pt
    enddo
  enddo
  n_tab = n_tab + n_tab2
!  
end subroutine intersection_3D_plan

subroutine intersection_3D_sphere(geom,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(sphere),intent(in) :: geom
  !f2py integer*1, dimension(1000) :: tab
  type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
  integer,intent(inout) :: n_tab
  !f2py integer*1, dimension(1000) :: list_point
  type(point),dimension(:),intent(in) :: list_point
  integer,intent(in) :: nlist
  real(rp),intent(in) :: t
  real(rp),intent(in) :: dx
  integer,intent(inout) :: nflag
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: P0,Pn,Pnp1,Q,P
  type(point) :: P0,Pn,Pnp1,Q,P
  !f2py integer*1, dimension(1000) :: tab_point,tab_point2
  type(point),dimension(100) :: tab_point,tab_point2
  integer :: j,k,jn
  integer :: np,np2,iflag,iface,n_bord
  integer :: isens
  real(rp) :: dl,ln,l0
  real(rp) :: h,hh,h6
  !f2py integer*1, dimension(1000) :: V
  type(vector) :: V
  logical :: is_loop,is_fin,is_border,is_border0,is_pres
  integer,parameter :: n=2
  real(rp),dimension(n) :: y,yl,yout,y0
  real(rp),dimension(2) :: diff1,diff2,diff3
  integer,parameter :: jmax =10000
  integer,dimension(10) :: aux
  integer :: naux
  !f2py integer*1, dimension(1000) :: ptr1,ptr2,ptr
  type(chaine_point),pointer :: ptr1,ptr2,ptr
  !f2py integer*1, dimension(1000) :: tab2
  type(chaine_point_pt),dimension(100,2) :: tab2
  integer :: n_tab2
  real(rp),dimension(:,:),allocatable :: tbound
  integer :: nbound
  
  np = 0 
  n_tab2 = 0 
  !ds = 0.01
  
  if(Symmetry)then
    allocate(tbound(2,2))
    tbound(1,1:2) = [0._rp,PI]
    tbound(2,1:2) = [geom%vmin,geom%vmax]
    nbound = 2
  else
    allocate(tbound(1,2))
    tbound(1,1:2) = [geom%vmin,geom%vmax]
    nbound = 1
  endif

  
  iface = geom%index
! Recherche point appartenant \E0 la surface  
  do j=1,nlist
    P = list_point(j)
    call common_int(P%face,P%nface,[iface],1,aux,naux)
    if(naux==1)then
        np=np+1
        tab_point(np) = P
    endif
  enddo
  
    if(Symmetry)then
    do j=1,np
      if(abs(tab_point(j)%coord(2)).lt.Epsilon)then
        tab_point(j)%bf = 2
      endif
    enddo
  endif

!  do j=1,n_tab
!    do k=1,2
!      P = tab(j,k)%pt%val
!      call common_int(P%face,P%nface,[iface],1,aux,naux)
!      if(naux==1)then
!        np = np+1
!        tab_point(np)=P
!      endif
!    enddo
!  enddo
! Recherche nouveaux points d'intersection  
  if(.not.Symmetry)then
    call point_inter_sph(geom,t,dx,tab_point2,np2)
  endif
  
  do j=1,np2
    P = tab_point2(j)
    if(Symmetry .and. abs(P%coord(2)).lt.Epsilon)then
      P%bf = 2
      P%nface = P%nface+1
      P%face(P%nface) = -1
    endif
    tab_point(np+j)=P
  enddo
  np=np+np2
  
  iflag = nflag
! Debut boucle calcul courbes d'intersection  
  do k=1,np
  
    nullify(ptr1)
    nullify(ptr2)
  
    P0 = tab_point(k)    
    do j=1,n_tab2
      ptr => tab2(j,1)%pt
      call point_in_chaine(ptr,P0,is_pres,0.6*dx,P,ierror)
      if(is_pres)then
        goto 100                         ! Test point deja present => fin boucle k
      endif
    enddo
    
    n_tab2 = n_tab2+1
    iflag = iflag+1
    P0%flag = iflag
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
       
    jn=1
    dl = 0._rp
    n_bord = 0
    isens = 1
        
    Pn = P0
    call cart2sph_f(Pn,geom,y)
    y0 = y
    
    call check_border_ND(y,tbound,nbound,y,yout,is_border0,ierror)
    if(is_border0)then
      ptr1%val%bf = 1
      if(Symmetry .and. abs(ptr1%val%coord(2)).lt.Epsilon)then
        ptr1%val%bf = 2
      endif
      n_bord = 1
    endif
    
    do while (.not.is_fin)
      
      h = dble(isens)*dsi
      hh=0.5*h
      h6=h/6.
         
! RK4 1ere passe
      call derive_sphere(Pn,geom,t,diff1)
      do j=1,n
        yl(j)=y(j)+hh*diff1(j)
      enddo
      yl(1) = mod(yl(1),TWOPI)
      call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)

  
! RK4 2eme passe
      if (.not.is_border) then
        call Computation_eq_sph3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
        call derive_sphere(Q,geom,t,diff2)
        do j=1,n
          yl(j)=y(j)+hh*diff2(j)
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif  

! RK4 3eme passe  
      if(.not.is_border) then
        call Computation_eq_sph3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
        call derive_sphere(Q,geom,t,diff3)
        do j=1,n
          yl(j)=y(j)+h*diff3(j)
          diff3(j) = diff2(j)+diff3(j)
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif
  
! RK4 4eme passe  
      if(.not.is_border)then
        call Computation_eq_sph3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
        call derive_sphere(Q,geom,t,diff2)
        do j=1,n
          yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif 

      if(.not.is_border)then
        yout = yl
      endif
  
      call Computation_eq_sph3D(yout(1),yout(2),geom,Pnp1%coord)
      
      call assign_vector_coord(V,Pnp1%coord - Pn%coord)
      call assign_vector_coord(V,Pnp1%coord - P0%coord)

      ln = geom%radius*sqrt((mod(yout(1)-y(1),TWOPI))**2 + (yout(2)-y(2))**2)
      l0 = geom%radius*sqrt((mod(yout(1)-y0(1),TWOPI))**2 + (yout(2)-y0(2))**2)
      
! Ajout point intersection     
      !dl = dl+V%length
      dl = dl+ln
      if(dl >= dx .and. .not.is_border)then
        !n_tab = n_tab+1
        Pnp1%nface = 2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 0
        Pnp1%nedge = 0
        Pnp1%edge = 0
        if(isens > 0)then
          call add_point_suiv(ptr2,Pnp1)
          ptr2 => ptr2%suiv
        else
          call add_point_prec(ptr1,Pnp1)
          ptr1 => ptr1%prec
        endif
        dl = 0._rp
      endif       
! Ajout point fronti\E8re      
      if(is_border)then
        Pnp1%nface=2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        if(abs(Pnp1%coord(2)).lt.Epsilon .and. Symmetry)then
          Pnp1%bf = 2
          Pnp1%nface = 3
          Pnp1%face(3) = -1
        else
          Pnp1%bf = 1        
        endif
        if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
          if(isens > 0)then
            ptr2%val=Pnp1
          else
            ptr1%val=Pnp1
          endif
          n_bord = n_bord+1
        elseif(abs(dl).gt.Epsilon)then
          if(isens > 0)then
            call add_point_suiv(ptr2,Pnp1)
            ptr2 => ptr2%suiv
          else
            call add_point_prec(ptr1,Pnp1)
            ptr1 => ptr1%prec
          endif
          n_bord = n_bord+1
        endif
!       Mise \E0 jour isens 
        if(isens == 1 .and. n_bord < 2)then
          isens = -1
          Pnp1 = P0
          !yout = cart2sph_f(P0,geom)
          yout = y0
        else
          is_fin = .true.
        endif
      endif
! Ajout point pour fermer la boucle 
      !is_loop = V0%length < 0.9*dsi .and. j>1 .and. isens==1      
      is_loop = l0 < 0.9*dsi .and. j>1 .and. isens==1
      if(is_loop)then
        if(dl < 0.5*dx)then
          if(isens > 0)then
            ptr2 => ptr2%prec
            deallocate(ptr2%suiv)
            ptr2%suiv => ptr1
            ptr1%prec => ptr2
          else
            ptr1 => ptr1%suiv
            deallocate(ptr1%prec)
            ptr1%prec => ptr2
            ptr2%suiv => ptr1
          endif
        else
          ptr2%suiv => ptr1
          ptr1%prec => ptr2
        endif
      endif
      
      is_fin = is_fin .or. is_loop .or. jn>10000
      
      jn=jn+1
      Pn = Pnp1
      y = yout
     
    enddo   
!    
!    tab2(n_tab2,1)%val  =  ptr1%val
!    tab2(n_tab2,1)%suiv => ptr1%suiv
!    tab2(n_tab2,1)%prec => ptr1%prec
!    ptr => ptr1%suiv
!    ptr%prec => tab2(n_tab2,1)
!    tab2(n_tab2,2)%val  =  ptr2%val
!    tab2(n_tab2,2)%suiv => ptr2%suiv
!    tab2(n_tab2,2)%prec => ptr2%prec 
!    ptr => ptr2%prec
!    ptr%suiv => tab2(n_tab2,2)
     tab2(n_tab2,1)%pt => ptr1
     tab2(n_tab2,2)%pt => ptr2

  if(jn>10000)then
    ierror=100
    goto 9999
  endif  
!
100 continue
!    
  enddo 
  
  nflag = iflag   
  
  do j=1,n_tab2
    do k=1,2
     tab(n_tab+j,k)%pt => tab2(j,k)%pt
    enddo
  enddo  
  n_tab = n_tab+n_tab2
  !deallocate(tab2)

9999 continue
  if(ierror/=0)then
    write(*,99),ierror
  endif
99 format('error #',i3,' : pb in computation of intersection.')
  
  if(allocated(tbound)) deallocate(tbound)
  
end subroutine intersection_3D_sphere

subroutine intersection_3D_cone(geom,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: geom
  type(cone),intent(in) :: geom
  !f2py integer*1, dimension(1000) :: tab
  type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
  integer,intent(inout) :: n_tab
  !f2py integer*1, dimension(1000) :: list_point
  type(point),dimension(:),intent(in) :: list_point
  integer,intent(in) :: nlist
  real(rp),intent(in) :: t
  real(rp),intent(in) :: dx
  integer,intent(inout) :: nflag
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: P0,Pn,Pnp1,Q,P
  type(point) :: P0,Pn,Pnp1,Q,P
  !f2py integer*1, dimension(1000) :: tab_point,tab_point2
  type(point),dimension(100) :: tab_point,tab_point2
  integer :: j,k
  integer :: np,np2,iflag,iface,n_bord
  integer :: isens
  real(rp) :: dl
  real(rp) :: h,hh,h6
  !f2py integer*1, dimension(1000) :: V,V0
  type(vector) :: V,V0
  logical :: is_loop,is_fin,is_border,is_border0,is_pres
  integer,parameter :: n=2
  real(rp),dimension(n) :: y,yl,yout,y0
  real(rp),dimension(2) :: diff1,diff2,diff3
  integer,parameter :: jmax =10000
  integer,dimension(10) :: aux
  integer :: naux
  !f2py integer*1, dimension(1000) :: ptr1,ptr2,ptr
  type(chaine_point),pointer :: ptr1,ptr2,ptr
  !f2py integer*1, dimension(1000) :: tab2
  type(chaine_point_pt),dimension(100,2) :: tab2
  integer :: n_tab2
  real(rp),dimension(:,:),allocatable :: tbound
  integer :: nbound
  
  np = 0 ; np2 = 0
  n_tab2 = 0 
  !ds = 0.01
  
  if(Symmetry)then
    allocate(tbound(2,2))
    tbound(1,1:2) = [0._rp,PI]
    tbound(2,1:2) = [geom%vmin,geom%vmax]
    nbound = 2
  else
    allocate(tbound(1,2))
    tbound(1,1:2) = [geom%vmin,geom%vmax]
    nbound = 1
  endif
  
  iface = geom%index
! Recherche point appartenant \E0 la surface 
  do j=1,nlist
    P = list_point(j)
    call common_int(P%face,P%nface,[iface],1,aux,naux)
    if(naux==1)then
        np=np+1
        tab_point(np) = P
    endif
  enddo

  if(Symmetry)then
    do j=1,np
      if(abs(tab_point(j)%coord(2)).lt.Epsilon)then
        tab_point(j)%bf = 2
      endif
    enddo
  endif
 
  do j=1,n_tab
    do k=1,2
      P = tab(j,k)%pt%val
      call common_int(P%face,P%nface,[iface],1,aux,naux)
      if(naux==1)then
        np = np+1
        tab_point(np)=P
      endif
    enddo
  enddo
! Recherche nouveaux points d'intersection  
  if(.not.Symmetry)then
    call point_inter_cone(geom,t,dx,tab_point2,np2)
  endif
  
  do j=1,np2
    P = tab_point2(j)
    if(Symmetry .and. abs(P%coord(2)).lt.Epsilon)then
      P%bf = 2
      P%nface = P%nface+1
      P%face(P%nface) = -1
    endif
    tab_point(np+j)=P
  enddo
  np=np+np2
  
  iflag = nflag
! Debut boucle calcul courbes d'intersection  
  do k=1,np
  
    nullify(ptr1)
    nullify(ptr2)
  
    P0 = tab_point(k)    
    do j=1,n_tab2
      ptr => tab2(j,1)%pt
      call point_in_chaine(ptr,P0,is_pres,0.6*dx,P,ierror)
      if(is_pres)then
        goto 100                         ! Test point deja present => fin boucle k
      endif
    enddo
    
    n_tab2 = n_tab2+1
    iflag = iflag+1
    P0%flag = iflag
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
       
    j=1
    dl = 0._rp
    n_bord = 0
    isens = 1
        
    Pn = P0
    call cart2cone(Pn,geom,y)
    y0 = y
    
    call check_border_ND(y,tbound,nbound,y,yout,is_border0,ierror)
    if(is_border0)then
      ptr1%val%bf = 1
      if(Symmetry .and. abs(ptr1%val%coord(2)).lt.Epsilon)then
        ptr1%val%bf = 2
      endif
      n_bord = 1
    endif
    
    do while (.not.is_fin)
      
      h = dble(isens)*dsi
      hh=0.5*h
      h6=h/6.
         
! RK4 1ere passe
      call derive_cone(Pn,geom,t,diff1)
      do j=1,n
        yl(j)=y(j)+hh*diff1(j)
      enddo
      yl(1) = mod(yl(1),TWOPI)
      call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)

  
! RK4 2eme passe
      if (.not.is_border) then
        call Computation_eq_cone3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q 
        call derive_cone(Q,geom,t,diff2)
        do j=1,n
          yl(j)=y(j)+hh*diff2(j)
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif  

! RK4 3eme passe  
      if(.not.is_border) then
        call Computation_eq_cone3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
        call derive_cone(Q,geom,t,diff3)
        do j=1,n
          yl(j)=y(j)+h*diff3(j)
          diff3(j) = diff2(j)+diff3(j)
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif
  
! RK4 4eme passe  
      if(.not.is_border)then
        call Computation_eq_cone3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
        call derive_cone(Q,geom,t,diff2)
        do j=1,n
          yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
        enddo
        yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif 

      if(.not.is_border)then
        yout = yl
      endif
  
      call Computation_eq_cone3D(yout(1),yout(2),geom,Pnp1%coord)
      
      call assign_vector_coord(V,Pnp1%coord - Pn%coord)
      call assign_vector_coord(V0,Pnp1%coord - P0%coord)
      
! Ajout point intersection     
      dl = dl+V%length
      if(dl >= dx .and. .not.is_border)then
        !n_tab = n_tab+1
        Pnp1%nface = 2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 0
        Pnp1%nedge = 0
        Pnp1%edge = 0
        if(isens > 0)then
          call add_point_suiv(ptr2,Pnp1)
          ptr2 => ptr2%suiv
        else
          call add_point_prec(ptr1,Pnp1)
          ptr1 => ptr1%prec
        endif
        dl = 0._rp
      endif       
! Ajout point fronti\E8re      
      if(is_border)then
        !Pnp1%nface=2
        !Pnp1%face(1:2) = [HouleRF%index,iface]
        !Pnp1%flag = iflag
        !if(abs(Pnp1%coord(2)).lt.Epsilon .and. Symmetry)then
        !  Pnp1%bf = 2
        !  Pnp1%nface = 3
        !  Pnp1%face(3) = -1
        !else
        !  Pnp1%bf = 1        
        !endif
        call update_tab_point(tab_point,k,np,Pnp1,dx/2.,ierror)
        if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
          if(isens > 0)then
            ptr2%val=Pnp1
          else
            ptr1%val=Pnp1
          endif
          n_bord = n_bord+1
        elseif(abs(dl).gt.Epsilon)then
          if(isens > 0)then
            call add_point_suiv(ptr2,Pnp1)
            ptr2 => ptr2%suiv
          else
            call add_point_prec(ptr1,Pnp1)
            ptr1 => ptr1%prec
          endif
          n_bord = n_bord+1
        endif
!       Mise \E0 jour isens 
        if(isens == 1 .and. n_bord < 2)then
          isens = -1
          Pnp1 = P0%coord
          call cart2cone(P0,geom,yout)
        else
          is_fin = .true.
        endif
      endif
! Ajout point pour fermer la boucle 
      is_loop = V0%length < 0.9*dsi .and. j>1 .and. isens==1      
      if(is_loop)then
        if(dl < 0.5*dx)then
          if(isens > 0)then
            ptr2 => ptr2%prec
            deallocate(ptr2%suiv)
            ptr2%suiv => ptr1
            ptr1%prec => ptr2
          else
            ptr1 => ptr1%suiv
            deallocate(ptr1%prec)
            ptr1%prec => ptr2
            ptr2%suiv => ptr1
          endif
        else
          ptr2%suiv => ptr1
          ptr1%prec => ptr2
        endif
      endif
      
      is_fin = is_fin .or. is_loop
      
      j=j+1
      Pn = Pnp1
      y = yout
     
    enddo   
!
     tab2(n_tab2,1)%pt => ptr1
     tab2(n_tab2,2)%pt => ptr2
!
100 continue
!    
  enddo 
  
  nflag = iflag   
  
  do j=1,n_tab2
    do k=1,2
     tab(n_tab+j,k)%pt => tab2(j,k)%pt
    enddo
  enddo  
  n_tab = n_tab+n_tab2
  
  if(allocated(tbound)) deallocate(tbound)
  
end subroutine intersection_3D_cone

subroutine intersection_3D_axisym(geom,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
    
    !f2py integer*1, dimension(1000)                    :: geom
    type(axisym),intent(in)                             :: geom
    !f2py integer*1, dimension(1000)                    :: tab
    type(chaine_point_pt),dimension(:,:),intent(inout)  :: tab
    integer,intent(inout)                               :: n_tab
    !f2py integer*1, dimension(1000)                    :: list_point
    type(point),dimension(:),intent(in)                 :: list_point
    integer,intent(in)                                  :: nlist
    real(rp),intent(in)                                 :: t
    real(rp),intent(in)                                 :: dx
    integer,intent(inout)                               :: nflag
    integer,intent(inout)                               :: ierror
    
    !f2py integer*1, dimension(1000)                    :: P0,Pn,Pnp1,Q,P
    type(point)                                         :: P0,Pn,Pnp1,Q,P
    !f2py integer*1, dimension(1000)                    :: tab_point
    type(point),dimension(100)                          :: tab_point
    integer                                             :: j,k,jn
    integer                                             :: np,iflag,iface,n_bord
    integer                                             :: isens
    real(rp)                                            :: dl,ln,l0,r
    real(rp)                                            :: h,hh,h6
    logical                                             :: is_loop,is_fin,is_border,is_border0,is_pres
    integer,parameter                                   :: n=2
    real(rp),dimension(n)                               :: y,yl,yout,y0
    real(rp),dimension(2)                               :: diff1,diff2,diff3
    integer,parameter                                   :: jmax =10000
    integer,dimension(10)                               :: aux
    integer                                             :: naux
    !f2py integer*1, dimension(1000)                    :: ptr1,ptr2,ptr
    type(chaine_point),pointer                          :: ptr1,ptr2,ptr
    !f2py integer*1, dimension(1000)                    :: tab2
    type(chaine_point_pt),dimension(100,2)              :: tab2
    integer                                             :: n_tab2
    real(rp),dimension(:,:),allocatable                 :: tbound
    integer                                             :: nbound
        
    ! This subroutine computes the intersection curve for an axisym geometry.
    
    np = 0 
    n_tab2 = 0 
        
    allocate(tbound(2,2))
    tbound(1,1:2) = [geom%umin,geom%umax]
    tbound(2,1:2) = [geom%vmin,geom%vmax]
    nbound = 2
    
    iface = geom%index
    
    ! Recherche point appartenant a la surface  
    do j=1,nlist
        P = list_point(j)
        call common_int(P%face,P%nface,[iface],1,aux,naux,ierror)
        if(naux==1)then
            np=np+1
            tab_point(np) = P
        endif
    enddo
    
    if(Symmetry)then
        do j=1,np
            if(abs(tab_point(j)%coord(2)).lt.Epsilon)then
                tab_point(j)%bf = 2
            endif
        enddo
    endif
    
    ! Debut boucle calcul courbes d'intersection
    iflag = nflag  
    do k=1,np
                
        nullify(ptr1)
        nullify(ptr2)
                
        P0 = tab_point(k)    
        do j=1,n_tab2
            ptr => tab2(j,1)%pt
            call point_in_chaine(ptr,P0,is_pres,0.4*dx,P,ierror)
            if(is_pres)then
                goto 100 ! Test point deja present => fin boucle k
            endif
        enddo
        
        n_tab2 = n_tab2+1
        iflag = iflag+1
        P0%flag = iflag
        call add_point_suiv(ptr1,P0)
        ptr2 => ptr1
        
        is_loop = .false.
        is_fin = .false. 
        
        jn=1
        dl = 0._rp
        n_bord = 0
        isens = 1
        
        Pn = P0
        call cart2axi(Pn,geom,y,ierror)
        if(ierror/=0)then
            ierror=101
            goto 9999
        endif  
        y0 = y
    
        call check_border_ND(y,tbound,nbound,y,yout,is_border0,ierror)
        if(is_border0)then
            ptr1%val%bf = 1
            if(Symmetry .and. abs(ptr1%val%coord(2)).lt.Epsilon)then
                ptr1%val%bf = 2
            endif
            n_bord = 1
        endif
                
        do while (.not.is_fin)
                        
            h = dble(isens)*dsi
            hh = 0.5*h
            h6 = h/6.
                        
            ! RK4 1ere passe
            call derive_axi(Pn,geom,t,diff1)
            do j=1,n
                yl(j)=y(j)+hh*diff1(j)
            enddo
            yl(1) = mod(yl(1),TWOPI)
            call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
                        
            ! RK4 2eme passe
            if (.not.is_border) then
                call eq_axi3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
                call derive_axi(Q,geom,t,diff2)
                do j=1,n
                    yl(j)=y(j)+hh*diff2(j)
                enddo
                yl(1) = mod(yl(1),TWOPI)
                call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
            endif  
                        
            ! RK4 3eme passe  
            if(.not.is_border) then
                call eq_axi3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
                call derive_axi(Q,geom,t,diff3)
                do j=1,n
                    yl(j)=y(j)+h*diff3(j)
                    diff3(j) = diff2(j)+diff3(j)
                enddo
                yl(1) = mod(yl(1),TWOPI)
                call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
            endif
                        
            ! RK4 4eme passe  
            if(.not.is_border)then
                call eq_axi3D(yl(1),yl(2),geom,Q%coord) ! Previously only Q
                call derive_axi(Q,geom,t,diff2)
                do j=1,n
                    yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
                enddo
                yl(1) = mod(yl(1),TWOPI)
                call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
            endif 
                        
            if(.not.is_border)then
                yout = yl
            endif
  
            call eq_axi3D(yout(1),yout(2),geom,Pnp1%coord)
            
            call radius_axi(geom,0.5*(yout(2)+y(2)),r)
            ln = r*sqrt((mod(yout(1)-y(1),TWOPI))**2 + (yout(2)-y(2))**2)
            l0 = r*sqrt((mod(yout(1)-y0(1),TWOPI))**2 + (yout(2)-y0(2))**2)
      
            ! Ajout point intersection     
            dl = dl+ln
            if(dl >= dx .and. .not.is_border)then
                Pnp1%nface = 2
                Pnp1%face(1:2) = [HouleRF%index,iface]
                Pnp1%flag = iflag
                Pnp1%bf = 0
                Pnp1%nedge = 0
                Pnp1%edge = 0
                if(isens > 0)then
                    call add_point_suiv(ptr2,Pnp1)
                    ptr2 => ptr2%suiv
                else
                    call add_point_prec(ptr1,Pnp1)
                    ptr1 => ptr1%prec
                endif
                dl = 0._rp
            endif
        
            ! Ajout point frontiere      
            if(is_border)then
                Pnp1%nface=2
                Pnp1%face(1:2) = [HouleRF%index,iface]
                Pnp1%flag = iflag
                if(abs(Pnp1%coord(2)).lt.Epsilon .and. Symmetry)then
                    Pnp1%bf = 2
                    Pnp1%nface = 3
                    Pnp1%face(3) = -1
                else
                    Pnp1%bf = 1        
                endif
        
                if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
                    if(isens > 0)then
                        ptr2%val=Pnp1
                    else
                        ptr1%val=Pnp1
                    endif
                    n_bord = n_bord+1
                elseif(abs(dl).gt.Epsilon)then
                    if(isens > 0)then
                        call add_point_suiv(ptr2,Pnp1)
                        ptr2 => ptr2%suiv
                    else
                        call add_point_prec(ptr1,Pnp1)
                        ptr1 => ptr1%prec
                    endif
                    n_bord = n_bord+1
                endif
        
                ! Mise a jour isens.
                if(isens == 1 .and. n_bord < 2)then
                    isens = -1
                    Pnp1 = P0
                    yout = y0
                else
                    is_fin = .true.
                endif
            endif
            
            ! Ajout point pour fermer la boucle    
            is_loop = l0 < 0.9*dsi .and. j>1 .and. isens==1
            if(is_loop)then
                if(dl < 0.5*dx)then
                    if(isens > 0)then
                        ptr2 => ptr2%prec
                        deallocate(ptr2%suiv)
                        ptr2%suiv => ptr1
                        ptr1%prec => ptr2
                    else
                        ptr1 => ptr1%suiv
                        deallocate(ptr1%prec)
                        ptr1%prec => ptr2
                        ptr2%suiv => ptr1
                    endif
                else
                    ptr2%suiv => ptr1
                    ptr1%prec => ptr2
                endif
            endif
            
            200 continue
            
            is_fin = is_fin .or. is_loop .or. jn>10000
      
            jn = jn+1
            Pn = Pnp1
            y = yout
     
        end do
        
        tab2(n_tab2,1)%pt => ptr1
        tab2(n_tab2,2)%pt => ptr2
        
        if(jn>10000)then
            ierror=100
            goto 9999
        endif  
    
        100 continue
        
    end do 
  
    nflag = iflag   
  
    do j = 1,n_tab2
        do k = 1,2            
            tab(n_tab+j,k)%pt => tab2(j,k)%pt
        end do
    end do  
    n_tab = n_tab + n_tab2

    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99 format('error #',i3,' : pb in computation of intersection.')
  
    if(allocated(tbound)) deallocate(tbound)
  
end subroutine intersection_3D_axisym

subroutine intersection_3D_wigley(geom,tab,n_tab,list_point,nlist,t,dx,nflag,ierror)
    
!f2py integer*1, dimension(1000)    :: geom
  type(TWigley),intent(in) :: geom
  !f2py integer*1, dimension(1000)    :: tab
  type(chaine_point_pt),dimension(:,:),intent(inout) :: tab
  integer,intent(inout) :: n_tab
  !f2py integer*1, dimension(1000)    :: list_point
  type(point),dimension(:),intent(in) :: list_point
  integer,intent(in) :: nlist
  real(rp),intent(in) :: t
  real(rp),intent(in) :: dx
  integer,intent(inout) :: nflag
  integer,intent(inout) :: ierror

  !f2py integer*1, dimension(1000)    :: P0,Pn,Pnp1,Q,P
  type(point) :: P0,Pn,Pnp1,Q,P
  !f2py integer*1, dimension(1000)    :: tab_point
  type(point),dimension(100) :: tab_point
  integer :: j,k,jn
  integer :: np,iflag,iface,n_bord
  integer :: isens
  real(rp) :: dl,ln,l0
  real(rp) :: h,hh,h6
  !f2py integer*1, dimension(1000)    :: V,V0
  type(vector) :: V,V0
  logical :: is_loop,is_fin,is_border,is_border0,is_pres
  integer,parameter :: n=2
  real(rp),dimension(n) :: y,yl,yout,y0
  real(rp),dimension(2) :: diff1,diff2,diff3,diff4
  integer,parameter :: jmax =10000
  integer,dimension(10) :: aux
  integer :: naux
  !f2py integer*1, dimension(1000)    :: ptr1,ptr2,ptr
  type(chaine_point),pointer :: ptr1,ptr2,ptr
  !f2py integer*1, dimension(1000)    :: tab2
  type(chaine_point_pt),dimension(100,2) :: tab2
  integer :: n_tab2
  real(rp),dimension(:,:),allocatable :: tbound
  integer :: nbound
!  
  np = 0 
  n_tab2 = 0 
!  
  allocate(tbound(2,2))
  tbound(1,1:2) = [geom%umin,geom%umax]
  tbound(2,1:2) = [geom%vmin,geom%vmax]
  nbound = 2

!  
  iface = geom%index
! Recherche point appartenant \E0 la surface  
  do j=1,nlist
    P = list_point(j)
    call common_int(P%face,P%nface,[iface],1,aux,naux)
    if(naux==1)then
        np=np+1
        tab_point(np) = P
    endif
  enddo
  
  if(Symmetry)then
    do j=1,np
      if(abs(tab_point(j)%coord(2)).lt.Epsilon)then
        tab_point(j)%bf = 2
      endif
    enddo
  endif
 
  iflag = nflag
! Debut boucle calcul courbes d'intersection  
  do k=1,np
  
    nullify(ptr1)
    nullify(ptr2)
  
    P0 = tab_point(k)    
    do j=1,n_tab2
      ptr => tab2(j,1)%pt
      call point_in_chaine(ptr,P0,is_pres,0.4*dx,P,ierror)
      if(is_pres)then
        goto 100                         ! Test point deja present => fin boucle k
      endif
    enddo
    
    n_tab2 = n_tab2+1
    iflag = iflag+1
    P0%flag = iflag
    call add_point_suiv(ptr1,P0)
    ptr2 => ptr1
    
    is_loop = .false.
    is_fin = .false. 
       
    jn=1
    dl = 0._rp
    n_bord = 0
    isens = 1
        
    Pn = P0
    call cart2Wigley(Pn,geom,y)
    y0 = y
    
    call check_border_ND(y,tbound,nbound,y,yout,is_border0,ierror)
    if(is_border0)then
      ptr1%val%bf = 1
      if(Symmetry .and. abs(ptr1%val%coord(2)).lt.Epsilon)then
        ptr1%val%bf = 2
      endif
      n_bord = 1
    endif
    
    do while (.not.is_fin)
      
      h = dble(isens)*dsi
      hh=0.5*h
      h6=h/6.
         
! RK4 1ere passe
      call derive_Wigley(Pn,geom,t,diff1)
      !diff1 = derive_Wigley(Pn,geom,t)
      do j=1,n
        yl(j)=y(j)+hh*diff1(j)
      enddo
      !yl(1) = mod(yl(1),TWOPI)
      call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)

  
! RK4 2eme passe
      if (.not.is_border) then
        call eq_Wigley3D(yl(1),yl(2),geom,Q%coord)
        call derive_Wigley(Q,geom,t,diff2)
        !Q = eq_Wigley3D(yl(1),yl(2),geom)
        !diff2 = derive2(Q,geom,t)
        do j=1,n
          yl(j)=y(j)+hh*diff2(j)
        enddo
        !yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif  

! RK4 3eme passe  
      if(.not.is_border) then
        call eq_Wigley3D(yl(1),yl(2),geom,Q%coord)
        call derive_Wigley(Q,geom,t,diff3)
        !Q = eq_Wigley3D(yl(1),yl(2),geom)
        !diff3 = derive2(Q,geom,t)
        do j=1,n
          yl(j)=y(j)+h*diff3(j)
          diff3(j) = diff2(j)+diff3(j)
        enddo
        !yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif
  
! RK4 4eme passe  
      if(.not.is_border)then
        call eq_Wigley3D(yl(1),yl(2),geom,Q%coord)
        call derive_Wigley(Q,geom,t,diff4)
        !Q = eq_Wigley3D(yl(1),yl(2),geom)
        !diff4 = derive2(Q,geom,t)
        do j=1,n
          !yl(j)=y(j)+h6*(diff1(j)+diff2(j)+2.*diff3(j))
          yl(j)=y(j)+h6*(diff1(j)+2._RP*diff2(j)+2._RP*diff3(j)+diff4(j))
        enddo
        !yl(1) = mod(yl(1),TWOPI)
        call check_border_ND(yl,tbound,nbound,y,yout,is_border,ierror)
      endif 

      if(.not.is_border)then
        yout = yl
      endif
  
      call eq_Wigley3D(yout(1),yout(2),geom,Pnp1%coord)
      
      call assign_vector_coord(V,Pnp1%coord - Pn%coord)
      call assign_vector_coord(V0,Pnp1%coord - P0%coord)

      !r = radius_axi(geom,0.5*(yout(2)+y(2)))
      !ln = r*sqrt((mod(yout(1)-y(1),TWOPI))**2 + (yout(2)-y(2))**2)
      !l0 = r*sqrt((mod(yout(1)-y0(1),TWOPI))**2 + (yout(2)-y0(2))**2)
      ln = V%length
      l0 = V0%length
      
! Ajout point intersection     
      dl = dl+ln
      if(dl >= dx .and. .not.is_border)then
        Pnp1%nface = 2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        Pnp1%bf = 0
        Pnp1%nedge = 0
        Pnp1%edge = 0
        if(isens > 0)then
          call add_point_suiv(ptr2,Pnp1)
          ptr2 => ptr2%suiv
        else
          call add_point_prec(ptr1,Pnp1)
          ptr1 => ptr1%prec
        endif
        dl = 0._rp
      endif       
! Ajout point fronti\E8re      
      if(is_border)then
        Pnp1%nface=2
        Pnp1%face(1:2) = [HouleRF%index,iface]
        Pnp1%flag = iflag
        if(abs(Pnp1%coord(2)).lt.Epsilon .and. Symmetry)then
          Pnp1%bf = 2
          Pnp1%nface = 3
          Pnp1%face(3) = -1
        else
          Pnp1%bf = 1        
        endif
        if(dl<dx/2. .and. abs(dl).gt.Epsilon)then
          if(isens > 0)then
            ptr2%val=Pnp1
          else
            ptr1%val=Pnp1
          endif
          n_bord = n_bord+1
        elseif(abs(dl).gt.Epsilon)then
          if(isens > 0)then
            call add_point_suiv(ptr2,Pnp1)
            ptr2 => ptr2%suiv
          else
            call add_point_prec(ptr1,Pnp1)
            ptr1 => ptr1%prec
          endif
          n_bord = n_bord+1
        endif
!       Mise \E0 jour isens 
        if(isens == 1 .and. n_bord < 2)then
          isens = -1
          Pnp1 = P0
          yout = y0
        else
          is_fin = .true.
        endif
      endif
! Ajout point pour fermer la boucle    
      is_loop = l0 < 0.9*dsi .and. j>1 .and. isens==1
      if(is_loop)then
        if(dl < 0.5*dx)then
          if(isens > 0)then
            ptr2 => ptr2%prec
            deallocate(ptr2%suiv)
            ptr2%suiv => ptr1
            ptr1%prec => ptr2
          else
            ptr1 => ptr1%suiv
            deallocate(ptr1%prec)
            ptr1%prec => ptr2
            ptr2%suiv => ptr1
          endif
        else
          ptr2%suiv => ptr1
          ptr1%prec => ptr2
        endif
      endif
      
      is_fin = is_fin .or. is_loop .or. jn>10000
      
      jn=jn+1
      Pn = Pnp1
      y = yout
     
    enddo   
!    
     tab2(n_tab2,1)%pt => ptr1
     tab2(n_tab2,2)%pt => ptr2

  if(jn>10000)then
    ierror=100
    goto 9999
  endif  
!
100 continue
!    
  enddo 
  
  nflag = iflag   
  
  do j=1,n_tab2
    do k=1,2
     tab(n_tab+j,k)%pt => tab2(j,k)%pt
    enddo
  enddo  
  n_tab = n_tab+n_tab2

9999 continue
  if(ierror/=0)then
    write(*,99),ierror
  endif
99 format('error #',i3,' : pb in comutation of intersection.')
  
  if(allocated(tbound)) deallocate(tbound)
  
end subroutine intersection_3D_wigley

subroutine update_tab_point(tab_point,k,np,P,dlmax,ierror)
  implicit none
  !f2py integer*1, dimension(1000) :: tab_point
  type(point),dimension(:),intent(inout) :: tab_point
  integer,intent(in) :: k
  integer,intent(inout) :: np
  !f2py integer*1, dimension(1000) :: P
  type(point),intent(inout) :: P
  real(rp),intent(in) :: dlmax
  integer,intent(inout) :: ierror
! local
  !f2py integer*1, dimension(1000) :: Q
  type(point) :: Q
  logical :: bool
  integer :: j
  
  !print*,'Begin update_point_inter...'
  !print*,' k = ',k
  !print*,' np = ',np
  !print*,'dlmax = ',dlmax
  !print*,'P = ',P%coord
  
  ierror= 0
  bool = .false.
  
  Q = tab_point(k)
  if(norm2(P%coord-Q%coord).lt.dlmax)then
    P = Q
    bool = .true.
  endif
  
  j = k+1
  do while(.not.bool .and. j<=np)
    Q = tab_point(j)
    if(norm2(P%coord-Q%coord).lt.dlmax)then
      P = Q
      tab_point(j:np-1) = tab_point(j+1:np)
      np=np-1
      bool = .true.
      !print*,'Correspondance found j = ',j
      !print*,'Q = ',Q%coord
    endif
    j = j+1
  enddo 
  
end subroutine update_tab_point

subroutine update_point_inter(tab,n_tab,tab2,n_tab2,ierror,InputData,NumBody)
  
    !f2py integer*1, dimension(1000)                    :: tab
    type(chaine_point_pt),dimension(:,:),intent(inout)  :: tab
    integer,intent(in)                                  :: n_tab
    !f2py integer*1, dimension(1000)                    :: tab2
    type(chaine_point_pt),dimension(:),intent(inout)    :: tab2
    integer,intent(inout)                               :: n_tab2
    integer,intent(inout)                               :: ierror
    !f2py integer*1, dimension(1000)                    :: InputData
    type(InputDataStruct),intent(in)                    :: InputData    ! Input data
    integer,intent(in)                                  :: NumBody      ! Body number = 1
    
    integer                                             :: j,k,j1,k1,n
    integer                                             :: iflag,iflag1,iflag2,icolor
    logical                                             :: is_pres
    !f2py integer*1, dimension(1000)                    :: P,Q
    type(point)                                         :: P,Q
    integer,dimension(50)                               :: aux
    integer                                             :: naux
    !f2py integer*1, dimension(1000)                    :: ptr2,ptr1,nouv,ptr
    type(chaine_point),pointer                          :: ptr2,ptr1,nouv,ptr
    integer,dimension(100,2)                            :: color_tab
    
    ierror = 0
    n_tab2 = 0
    iflag = 0
	
    do j = 1,n_tab
        ptr1 => tab(j,1)%pt
        ptr2 => tab(j,2)%pt
        if(associated(ptr2%suiv,ptr1).and.associated(ptr1%prec,ptr2) .or. ptr1%val%bf==2 .and. ptr2%val%bf==2)then
            iflag = iflag+1
            color_tab(j,1) = iflag
            color_tab(j,2) = iflag
        else
            do k = 1,2
                P = tab(j,k)%pt%val
                j1 = j
                is_pres = .false.
                do while (j1 < n_tab .and. .not. is_pres)
                    k1 = 0
                    j1=j1+1
                    do while (k1 < 2 .and. .not. is_pres)
                        k1=k1+1         
                        Q = tab(j1,k1)%pt%val
                        is_pres = norm2(P%coord-Q%coord).le.0.5_RP*InputData%dx2(NumBody)                                                
                    end do
                end do
                if(is_pres)then
                    iflag = iflag+1
                    color_tab(j,k)=iflag
                    color_tab(j1,k1)=iflag
                    call entire_int(P%face,P%nface,Q%face,Q%nface,aux,naux,ierror)
                    P%nface = naux
                    P%face(1:naux) = aux(1:naux)
                    P%bf = 1
                    tab(j,k)%pt%val = P
                    if(.not.associated(tab(j1,k1)%pt%prec) .and. .not.associated(tab(j,k)%pt%prec))then
                        if(not(associated(tab(j1,k1)%pt%suiv)))then
                            ierror = 1492 ! Year of the discover of a new continent: the pointers in Fortran!
                            goto 9999
						end if
                        tab(j,k)%pt%prec => tab(j1,k1)%pt%suiv
                        tab(j1,k1)%pt%suiv%prec => tab(j,k)%pt
                    elseif(.not.associated(tab(j1,k1)%pt%prec) .and. .not.associated(tab(j,k)%pt%suiv))then
                        tab(j,k)%pt%suiv => tab(j1,k1)%pt%suiv
                        tab(j1,k1)%pt%suiv%prec => tab(j,k)%pt
                    elseif(.not.associated(tab(j1,k1)%pt%suiv) .and. .not.associated(tab(j,k)%pt%prec))then
                        tab(j,k)%pt%prec => tab(j1,k1)%pt%prec
                        tab(j1,k1)%pt%prec%suiv => tab(j,k)%pt
                    elseif(.not.associated(tab(j1,k1)%pt%suiv) .and. .not.associated(tab(j,k)%pt%suiv))then
                        tab(j,k)%pt%suiv => tab(j1,k1)%pt%prec
                        tab(j1,k1)%pt%prec%suiv => tab(j,k)%pt
                    else
                        ierror = 100
                        goto 9999
                    end if
                end if
            end do
        end if
    end do
    
    do j=1,n_tab
        iflag1 = color_tab(j,1)
        iflag2 = color_tab(j,2)
        if(iflag1/=iflag2)then
            where (color_tab==iflag2)
                color_tab = iflag1
            end where
        end if
    end do
    
    do icolor=1,iflag
        n=0
        do j=1,n_tab
            do k=1,2
                is_pres = color_tab(j,k)==icolor
                if(is_pres)then
                    n=n+1
                    if(n==1)then
                        n_tab2 = n_tab2+1
                        tab2(n_tab2)%pt => tab(j,k)%pt
                    elseif(n>1 .and. associated(tab(j,k)%pt%suiv) .and. associated(tab(j,k)%pt%prec))then
                        allocate(nouv)
                        nouv%val = tab(j,k)%pt%val
                        ptr1 => tab(j,k)%pt%prec
                        ptr2 => tab(j,k)%pt%suiv
                        ptr  => tab(j,k)%pt
                        nouv%suiv => ptr2
                        nouv%prec => ptr1
                        if(associated(ptr1%suiv,ptr))then
                            ptr1%suiv => nouv
                        elseif(associated(ptr1%prec,ptr))then
                            ptr1%prec => nouv
                        else
                            ierror = 130
                            goto 9999
                        endif
                        if(associated(ptr2%suiv,ptr))then
                            ptr2%suiv => nouv
                        elseif(associated(ptr2%prec,ptr))then
                            ptr2%prec => nouv
                        else
                            ierror = 131
                            goto 9999
                        endif
                    elseif(n>1 .and. Symmetry .and. tab(j,k)%pt%val%bf==2 .and. k==1)then
                        allocate(nouv)
                        nouv%val  = tab2(n_tab2)%pt%val
                        ptr1 => tab2(n_tab2)%pt%prec
                        ptr2 => tab2(n_tab2)%pt%suiv
                        ptr  => tab2(n_tab2)%pt
                        nouv%suiv => ptr2
                        nouv%prec => ptr1
                        if(associated(ptr1))then
                            ptr1%suiv => nouv
                        endif
                        if(associated(ptr2))then
                            ptr2%prec => nouv
                        endif
                        tab2(n_tab2)%pt => tab(j,k)%pt
                    elseif(n>1 .and. Symmetry .and. tab(j,k)%pt%val%bf==2 .and. k==2)then
                        allocate(nouv)
                        nouv%val = tab(j,k)%pt%val
                        ptr1 => tab(j,k)%pt%prec
                        ptr2 => tab(j,k)%pt%suiv
                        nouv%suiv => ptr2
                        nouv%prec => ptr1
                        if(associated(ptr1))then
                            ptr1%suiv => nouv
                        end if
                        if(associated(ptr2))then
                            ptr2%prec => nouv
                        end if   
                    end if
                end if
            end do
        end do
    end do
    
    if(n_tab2<1)then
        ierror = 120
        goto 9999
    end if
    
    9999 continue
    if(ierror/=0)then
        print*,'** error #',ierror,' : update_point_inter'
    end if  
    
end subroutine update_point_inter

subroutine compute_intersection(t,fgeom_vect,tab2,n_tab2,ilink,ierror,InputData,dx2,n_tab)
    
    real(rp),intent(in)                             :: t                    ! Current time.
    !f2py integer*1, dimension(1000)                :: fgeom_vect
    type(type_GeomVect),intent(inout)               :: fgeom_vect           ! Geometries.
    !f2py integer*1, dimension(1000)                :: tab2
    type(chaine_point_pt),dimension(:),intent(inout):: tab2                 ! Table of the pointers toward the intersection points.
    integer,intent(inout)                           :: n_tab2               ! Number of intersection curves.
    integer,intent(in)                              :: ilink
    integer,intent(inout)                           :: ierror               ! Error flag.
    !f2py integer*1, dimension(1000)                :: InputData
    type(InputDataStruct),intent(in)                :: InputData            ! Input data.
    real(rp),intent(in)                             :: dx2                  ! dx2.
    integer,intent(inout),optional                  :: n_tab                ! Number of intersection lines.
            
    !f2py integer*1, dimension(1000)                :: arete
    type(Garete)                                    :: arete                ! Edge
    !f2py integer*1, dimension(1000)                :: cercle
    type(GCercle)                                   :: cercle               ! Circle
    !f2py integer*1, dimension(1000)                :: polyline
    type(GPolyline)                                 :: polyline             ! Polyline
    integer                                         :: nflag,nedge
    !f2py integer*1, dimension(1000)                :: tab
    type(chaine_point_pt),dimension(100,2)          :: tab  
    integer                                         :: n,jj,k,nl,np,nk,jk
    !f2py integer*1, dimension(1000)                :: tab_loc,tab_point
    type(point),dimension(100)                      :: tab_loc,tab_point
    !f2py integer*1, dimension(1000)                :: rep
    type(repere3d)                                  :: rep
    !f2py integer*1, dimension(1000)                :: P1,P2
    type(point)                                     :: P1,P2
    !f2py integer*1, dimension(1000)                :: ptr
    type(chaine_point),pointer                      :: ptr
    integer                                         :: jjj                  ! Loop parameter
        
    ! This subroutine computes the intersection line of the body.
        
    ierror = 0  
    
    nflag = 0
    n_tab = 0
    np = 0
    
    ! Modification of the geometry of the cylinders in order to compute the intersection curve
    call Active_point_arete(fgeom_vect,InputData)
    
    ! Lines: 1D
    if(is_body)then
        do jjj = 1,NBodies
            
            ! Edges
            n = fgeom_vect%geom(jjj)%narete
            do jj=1,n
                arete = fgeom_vect%geom(jjj)%arete(jj)
                call create_tab_intersect(arete%P1,arete%P2,t,dx2,tab_loc,nl,arete%face(1))
                do k=1,nl
                    call add_point_face(tab_loc(k),[arete%face(2)],1)
                    tab_loc(k)%nedge = 1
                    tab_loc(k)%edge(1) = arete%index
                enddo
                tab_point(np+1:np+nl) = tab_loc(1:nl)
                np = np+nl
            enddo
                        
            ! Circles
            !n = fgeom_vect%geom(jjj)%ncercle
            !do jj=1,n
            !    cercle = fgeom_vect%geom(jjj)%cercle(jj)
            !    rep = cercle%repere
            !    P1 = cercle%P1
            !    P2 = cercle%P2
            !    nl = 0
            !    call create_tabcurve_intersect(P1,P2,t,cercle,dx2,tab_loc,nl,cercle%face(1))
            !    do k=1,nl
            !        call add_point_face(tab_loc(k),[cercle%face(2)],1)
            !        tab_loc(k)%nedge = 1
            !        tab_loc(k)%edge(1) = cercle%index
            !    enddo
            !    tab_point(np+1:np+nl) = tab_loc(1:nl)
            !    np=np+nl
            !enddo
            
            ! Polylines
            n = fgeom_vect%geom(jjj)%npolyline
            do jj = 1,n
                polyline = fgeom_vect%geom(jjj)%polyline(jj)   
                rep = polyline%repere
                nk = polyline%n-1
                
                do jk = 1,nk
                    
                    call loc2cart(polyline%P(jk,:),rep,P1%coord)
                    call loc2cart(polyline%P(jk+1,:),rep,P2%coord)
                    nl = 0
                    
                    call create_tab_intersect(P1,P2,t,dx2,tab_loc,nl,polyline%face(1))
                                        
                    do k = 1,nl                        
                        call add_point_face(tab_loc(k),[polyline%face(2)],1)
                        tab_loc(k)%nedge = 1
                        tab_loc(k)%edge(1) = polyline%index                        
                    end do
                    tab_point(np+1:np+nl) = tab_loc(1:nl)
                    np = np + nl
                end do
            end do
        end do
    end if
        
    ! Correction of the geometry of the cylinders
    call Correction_Geom(fgeom_vect,InputData)
    
    if(Symmetry)then
        do jj=1,np
            if(abs(tab_point(jj)%coord(2)).lt.Epsilon)then
                nedge = tab_point(jj)%nedge 
                tab_point(jj)%edge(nedge+1) = ilink
                tab_point(jj)%nedge = nedge+1
            endif
        enddo
    endif
    
    ! Surface: 2D
    if(is_body)then
        do jjj = 1,NBodies
            
            ! Cylinders
            n = fgeom_vect%geom(jjj)%ncylindre
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%cylindre(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)  
                if(ierror/=0)then
                    ierror = 100
                    goto 9999
                endif
            enddo
                        
            ! Discs
            n = fgeom_vect%geom(jjj)%ndisque
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%disque(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
                if(ierror/=0)then
                    ierror = 200
                    goto 9999
                endif
            enddo
                        
            ! Plans
            n = fgeom_vect%geom(jjj)%nplan
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%plan(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
                if(ierror/=0)then
                    ierror = 300
                    goto 9999
                endif
            enddo
            
            ! Cones    
            n = fgeom_vect%geom(jjj)%ncone
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%cone(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
                if(ierror/=0)then
                    ierror = 400
                    goto 9999
                endif
            enddo
            
            ! Spheres
            n = fgeom_vect%geom(jjj)%nsphere
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%sphere(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
                if(ierror/=0)then
                    ierror = 500
                    goto 9999
                endif
            enddo
                        
            ! Axisyms
            n = fgeom_vect%geom(jjj)%naxisym
            do jj=1,n
                call intersection_3D_2(fgeom_vect%geom(jjj)%axisym(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
                if(ierror/=0)then
                    ierror = 600
                    goto 9999   
                endif
            enddo
			
			! Wigley hull
			n = fgeom_vect%geom(jjj)%nwigley
		    do jj=1,n
		      	call intersection_3D_2(fgeom_vect%geom(jjj)%wigley(jj),tab,n_tab,tab_point,np,t,dx2,nflag,ierror)
		      	if(ierror/=0)then
		        	ierror = 700
		        	goto 9999   
		      	endif
		    enddo
            
        end do
    endif
    
    if(n_tab.gt.0)then
    
        call update_point_inter(tab,n_tab,tab2,n_tab2,ierror,InputData,1) ! 1 = NumBody for dx2 (better if changed)

        if(ierror/=0)then
            ierror = 700
            goto 9999
        endif
        
        ! y coordinate of the end points is 0 in case of Symmetry.
        if(Symmetry)then
            do jj=1,n_tab2
                ptr => tab2(jj)%pt
                do while(associated(ptr))
                    if(ptr%val%bf == 2) ptr%val%coord(2) = 0._rp
                    ptr => ptr%suiv
                enddo
            enddo
        endif
        
        call write_intersection(tab2(:),n_tab2,t,ierror)
        if(ierror/=0)then
            ierror = 800
            goto 9999
        endif        
    else
        n_tab2 = 0
    endif
    
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
        call write_intersection_debug(tab(:,:),n_tab,ierror)
    endif
    99  format('** error #',i3,' : pb. calcul intersection')  
        
end subroutine compute_intersection

subroutine Active_point_arete(fgeom_vect,InputData)
    
    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(inout)   :: fgeom_vect   ! Geometry of the floaters
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data
    
    integer                             :: jj           ! Loop parameter
    
    ! This subroutine gives the good number of points and edges of the body geometries which were done to obtain a good intersection curve.
    
    if(is_body)then
        do jj = 1,NBodies
            if(InputData%igtype(jj)==2)then ! Cylinder
                if(not(Symmetry))then
                    fgeom_vect%geom(jj)%narete = 2
                    fgeom_vect%geom(jj)%npoint = 0
                end if
            end if
        end do
    end if

end subroutine Active_point_arete

subroutine Correction_Geom(fgeom_vect,InputData)

    !f2py integer*1, dimension(1000)    :: fgeom_vect
    type(type_GeomVect),intent(inout)   :: fgeom_vect   ! Geometry of the floaters
    !f2py integer*1, dimension(1000)    :: InputData
    type(InputDataStruct),intent(in)    :: InputData    ! Input data
    
    integer                             :: jj           ! Loop parameter
    
    ! This subroutine zeroes the number of points and edges of the body geometries which were done to obtain a good intersection curve.
    
    if(is_body)then
        do jj = 1,NBodies
            if(InputData%igtype(jj) == 2)then ! Cylinder
                if(not(Symmetry))then
                    fgeom_vect%geom(jj)%narete = 0
                    fgeom_vect%geom(jj)%npoint = 0
                end if
            end if
        end do
    end if
    
end subroutine Correction_Geom

end module SolvNum  
