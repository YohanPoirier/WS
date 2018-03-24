module GeomFonct
use Constantes
use GeomStruct
use MeshStruct
implicit none

interface assignment (=)
  module procedure assign_point_mvertex
end interface

contains

subroutine assign_point_mvertex(this,a)
    
    !f2py integer*1, dimension(1000)    :: this
    type(Point),intent(out)             :: this
    !f2py integer*1, dimension(1000)    :: a
    type(Mvertex),intent(in)            :: a
    
    integer                             :: j,nface
    
    ! This subroutine initializes a point from a MVertex structure.
    
    do j=1,3
        this%coord(j) = a%coord(j)
    enddo
    nface = a%nface
    this%nface = nface
    this%face(1:nface) = a%face(1:nface)
    this%bf = -9
    this%nedge = a%nedge
    this%edge(1:a%nedge) = a%edge(1:a%nedge)
            
end subroutine assign_point_mvertex

subroutine assign_base(C,V,r,base)
  !f2py integer*1, dimension(1000) :: C
  type(point),intent(in) :: C
  !f2py integer*1, dimension(1000) :: V
  type(vector),intent(in) :: V
  real(rp),intent(in) :: r
  !f2py integer*1, dimension(1000) :: base
  type(cylindrical_surface2),intent(out) :: base
  !base%name   = ''
  base%centre = C
  base%normal = V
  base%radius = r
end subroutine assign_base
! *************************************************************************
! Fonction de changement de repere
! *************************************************************************

subroutine mat_rotation2(y,gamma,phi,O3,yout)
  implicit none
  real(rp),dimension(3),intent(in) :: y
  real(rp),intent(in) :: gamma,phi
  real(rp),dimension(3),intent(out) :: yout
  real(rp),dimension(3) :: O3
  real(rp),dimension(3,3) :: R_gamma,R_phi
  yout=y
  R_gamma = reshape((/cos(gamma),0._rp,sin(gamma),0._rp,1._rp,0._rp,-sin(gamma),0._rp,cos(gamma)/),(/3,3/))
  R_phi = reshape((/cos(phi),sin(phi),0._rp,-sin(phi),cos(phi),0._rp,0._rp,0._rp,1._rp/),(/3,3/))
  yout = matmul(R_gamma,yout)-matmul(R_gamma,O3)+O3
  yout = matmul(R_phi,yout)-matmul(R_phi,O3)+O3
end subroutine mat_rotation2

subroutine mat_rotation(y,angle,dir,O3,yout)
  
    real(rp),dimension(3),intent(in)    :: y        ! Vector before the rotation.
    integer,intent(in)                  :: dir      ! 1 = rotation around x, 2 = rotation around y, 3 = rotation around z.
    real(rp),intent(in)                 :: angle    ! Angle.
    real(rp),dimension(3)               :: O3       ! Origin.
    real(rp),dimension(3)               :: yout     ! Vector after the rotation.
    
    real(rp),dimension(3,3)             :: R_angle  ! Rotation matrix.
    
    ! Rotation matrix.
    select case (dir) 
    case(1) ! Phi
        R_angle = reshape((/1._rp,0._rp,0._rp,0._rp,cos(angle),sin(angle),0._rp,-sin(angle),cos(angle)/),(/3,3/))
    case (2) ! Theta
        R_angle = reshape((/cos(angle),0._rp,-sin(angle),0._rp,1._rp,0._rp,sin(angle),0._rp,cos(angle)/),(/3,3/))
    case (3) ! Psi
        R_angle = reshape((/cos(angle),sin(angle),0._rp,-sin(angle),cos(angle),0._rp,0._rp,0._rp,1._rp/),(/3,3/))
    end select
    
    ! Rotation of the vector.
    yout = matmul(R_angle,y)-matmul(R_angle,O3)+O3
    
 end subroutine mat_rotation
 

! ******************************************************************
! Fonctions generales
! ******************************************************************
!
! Calcul coordonnées locales
!
subroutine cart2loc(P,rep,param)
    
    real(rp),dimension(3),intent(in)    :: P        ! Cartesian coordinate vector.
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep      ! Frame.
    real(rp),dimension(3),intent(out)   :: param    ! Local coordinate vector.

    real(rp),dimension(3)               :: O,Q      ! Points.
    real(rp),dimension(3)               :: e1,e2,e3 ! Unit vectors.
    
    ! This subroutine transforms a Cartesian coordinate vector into a local coordinate one.
    
    ! Frame.
    O = rep%origine
    e1 = rep%e1
    e2 = rep%e2
    e3 = rep%e3
    
    ! Local coordinates.
    Q = P-O
    param(1) = dot_product(Q,e1)
    param(2) = dot_product(Q,e2)
    param(3) = dot_product(Q,e3)

end subroutine cart2loc

! Calcul coordonnées dans le repère global

subroutine loc2cart(param,rep,x)
    
    real(rp),dimension(3),intent(in)    :: param    ! Local coordinate vector.
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep      ! Frame.
    real(rp),dimension(3),intent(out)   :: x        ! Cartesian coordinate vector.

    real(rp),dimension(3)               :: O        ! Origin.
    real(rp),dimension(3)               :: e1,e2,e3 ! Unit vector.
    
    ! This subroutine transforms the local coordinates into the Cartesian coordinates.
    
    ! Frame.
    O = rep%origine
    e1 = rep%e1
    e2 = rep%e2
    e3 = rep%e3
    
    ! Cartesian coordinates.
    x = O + param(1)*e1 + param(2)*e2 + param(3)*e3
    
end subroutine loc2cart  

! Calcul coordonnées sphérique d'un point dans un repère donné

subroutine cart2sph(P,rep,param)
    
    real(rp),dimension(3),intent(in)    :: P                ! Cartesian coordinate vector.
    !f2py integer*1, dimension(1000)    :: rep
    type(repere3d),intent(in)           :: rep              ! Frame.
    real(rp),dimension(3),intent(out)   :: param            ! Spheric coordinates.
    
    real(rp),dimension(3)               :: x                ! Point.
    real(rp)                            :: inv_r,theta,phi  !
    
    ! This subroutine transforms a Cartesian coordinate vector into a spheric coordinate one.
    
    call cart2loc(P,rep,x)
    
    param(1) = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)) ! Radius.
    
    if(param(1) > Epsilon)then
        inv_r = 1./param(1)
        
        ! In case where |x(1)*inv_r| is higher than 1 (numerical problem).
        if(x(1)*inv_r.gt.1_RP)then
            theta = PI
        else if(x(1)*inv_r.lt.(-1._RP))then
            theta = 0._RP
        else
            theta = acos(x(1)*inv_r)
            if(x(2) < -Epsilon)then
                theta = 2.*PI-theta
            endif
        end if
        
        ! In case where |x(3)*inv_r| is higher than 1 (numerical problem).
        if(x(3)*inv_r.gt.1_RP)then
            phi = PI
        else if(x(3)*inv_r.lt.(-1._RP))then
            phi = 0._RP
        else
            phi = acos(x(3)*inv_r)
        end if
        
    else
        theta = 0._rp
        phi   = 0._rp
    endif
  
    param(2) = theta
    param(3) = phi 

end subroutine cart2sph

! Calcul coordoonée global

subroutine sph2cart(P,rep,x)
  implicit none
  real(rp),dimension(3),intent(in) :: P
  !f2py integer*1, dimension(1000) :: rep
  type(repere3d),intent(in) :: rep
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp) :: r,theta,phi  
  
  r = P(1)
  theta = P(2)
  phi = P(3)
  x(1) = r*cos(theta)
  x(2) = r*sin(theta)
  x(3) = r*cos(phi)
  call loc2cart(x,rep,x) 

end subroutine sph2cart


! ******************************************************************
! Fonctions générales
! ******************************************************************
!
subroutine Computation_eq_line(P1,P2,s,eq_line)
    
    !f2py integer*1, dimension(1000)    :: P1,P2
    type(point),intent(in)              :: P1,P2    ! Boundary points
    real(rp),intent(in)                 :: s        ! Curvilinear abscissa
    real(rp),dimension(3),intent(out)   :: eq_line  ! Vector
    
    ! This subroutine computes the equation of a straight line.
    
    eq_line(1) = s*P2%coord(1)+(1._rp-s)*P1%coord(1)
    eq_line(2) = s*P2%coord(2)+(1._rp-s)*P1%coord(2)
    eq_line(3) = s*P2%coord(3)+(1._rp-s)*P1%coord(3)
    
end subroutine Computation_eq_line

subroutine cart2cer(P,cercle,param)
    
    real(rp),dimension(3),intent(in)    :: P        ! Cartesian coordinate vector.
    !f2py integer*1, dimension(1000)    :: cercle
    type(GCercle),intent(in)            :: cercle   ! Circle.
    real(rp),intent(out)                :: param    ! Polar coordinate vector.
    ! local
    real(rp)                            :: inv_r,x  !
    real(rp),dimension(3)               :: Q        ! Point.
    integer                             :: ierror   ! Error flag.
    
    ! This subroutine transforms a Cartesian coordinate vector into a polar coordinate one.
    
    ierror = 0
  
    inv_r = 1/cercle%rayon
    call cart2loc(P,cercle%repere,Q)
    x = Q(1)*inv_r
  
    ! In case where |x| is higher than 1 (numerical problem).
    if(x.gt.1._RP)then
        x = 1._RP
    else if(x.lt.(-1._RP))then
        x = -1._RP
    end if
  
    if (abs(x+1.).lt.0.001) then
        param = PI
    elseif(abs(x-1.).lt.0.001) then
        param = 0._rp
    elseif(abs(x).lt.1)then
        param = dacos(x)
    else
        ierror = 1
        goto 9999
    endif
    if (Q(2).lt.0) then
        param = 2.*PI-param
    endif
  
    9999 continue
    if (ierror .ne. 0)then
        print*,""
        print*,"P (cart) = ",P
        print*,"Q (loc) = ",Q
        print*,"x = Q(1)/r = ",x
        write(*,*) '** error : calcul coordonne cercle'  
        call exit()
    endif  
    
end subroutine cart2cer  

subroutine cer2cart(param,cercle,x)
  implicit none
  real(rp),intent(in) :: param
  !f2py integer*1, dimension(1000) :: cercle
  type(GCercle),intent(in) :: cercle
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp) :: r
!  
  r = cercle%rayon
  x = [r*cos(param),r*sin(param),0._rp]
  call loc2cart(x,cercle%repere,x)
!   
end subroutine cer2cart 

subroutine eq_curve(cercle,P1,P2,s,x)
  implicit none
  !f2py integer*1, dimension(1000) :: cercle
  type(GCercle),intent(in) :: cercle
  real(rp),dimension(3),intent(in) :: P1,P2
  real(rp),intent(in) :: s
  real(rp),dimension(3),intent(out) :: x
! local
  real(rp) :: theta1,theta2,theta
  integer :: ierror

  ierror = 0

  call cart2cer(P1,cercle,theta1)
  call cart2cer(P2,cercle,theta2)
  
  if(abs(theta2-theta1).lt.Epsilon)then
    theta2 = TWOPI
  elseif(theta2 < theta1)then
    theta2 = theta2+2.*PI
  endif  
  
  theta = s*theta2 + (1._rp-s)*theta1
  
  call cer2cart(theta,cercle,x)
  
  if(ierror/=0)then 
     print*,'** error : calcul coordonnee du cercle'
  endif
  
end subroutine eq_curve

subroutine segment_intersect_2D(A,B,C,D,res)

    real(rp),dimension(2),intent(in)    :: A,B,C,D          ! Endpoints of the two line segments.
    logical,intent(out)                 :: res              ! = true if intersection, false if no intersection between the two line segments.
      
    real(rp),dimension(2)               :: r,s,k            ! AB, CD, AC.
    real(rp)                            :: denom,num1,num2  ! r vect s, k vect s, k vect r.
    real(rp)                            :: t,u              ! k.r/(r.r), (k+s).r/(r.r).
    
    ! This subroutine detects if two line segments intersect or not.
    
    ! The following results come from the web site:
    ! https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    
    ! The two lines are AB and CD.
    
    res = .false.
    r = [B(1)-A(1),B(2)-A(2)] ! r = AB
    s = [D(1)-C(1),D(2)-C(2)] ! s = CD
    k = [C(1)-A(1),C(2)-A(2)] ! k = AC
    denom = r(1)*s(2)-r(2)*s(1) ! r vect s
    
    ! Four cases to study.
    if (abs(denom) > Epsilon) then ! r vect s != 0
        ! t = (k vect s)/(r vect s)
        num1 = k(1)*s(2)-k(2)*s(1) 
        t = num1/denom
        
        ! u = (k vect r)/(r vect s)
        num2 = k(1)*r(2)-k(2)*r(1)
        u = num2/denom
        
        ! If 0<=t<=1 and 0<=u<=1, the two line segments meet at the point A + t*r = C + u*s.
        if (t.gt.Epsilon .and. t.lt.(1._rp-Epsilon) .and. &
        &   u.gt.Epsilon .and. u.lt.(1._rp-Epsilon)) then 
            res = .true.
        else
            ! The two line segments are not parallel but not intersect.
            res = .false.
        endif
    else ! r vect s = 0
        num2 = k(1)*r(2)-k(2)*r(1) ! k vect r
        denom = r(1)*r(1)+r(2)*r(2) ! r.r
        if (abs(num2) < Epsilon)then ! r vect s = 0 and k vect r = 0
            ! AB and CD are colinear.
            
            ! Expression of the endpoints of CD in terms of the equation of AB.
            ! t = k.r/(r.r)
            t = k(1)*r(1)+k(2)*r(2) 
            t = t/denom
            ! u = (k+s).r/(r.r) = t + (s.r)/(r.r)
            u = (k(1)+s(1))*r(1)+(k(2)+s(2))*r(2)
            u = u/denom
            
            ! If t and u intersects the interval [0,1], then AB and CD are colinear and overlapping, otherwse they are colinear and disjoint.
            res = t.gt.Epsilon .and. t.lt.(1._rp-Epsilon) ! t > 0 and t < 1.
            res = res .or. u.gt.Epsilon .and. u.lt.(1._rp-Epsilon) ! u > 0 and u < 1.
            res = res .or. t.lt.Epsilon .and. u.gt.(1._rp-Epsilon) ! t < 0 and u > 1.
            res = res .or. u.lt.Epsilon .and. t.gt.(1._rp-Epsilon) ! u < 0 and t > 1.
            res = res .and. .not.(abs(t).lt.Epsilon .and. abs(1._rp-u).lt.Epsilon .or. abs(1._rp-t).lt.Epsilon .and. abs(u).lt.Epsilon) ! t != 0 and u != 1 and t != 1 and u != 0.
        else
            ! AB and CD are parallel and non-intersecting.
            res = .false.
        endif
    endif
      
end subroutine segment_intersect_2D

subroutine segment_intersect_2D_2(A,B,C,D,res)
  real(rp),dimension(2),intent(in) :: A,B,C,D
  logical :: res
  !local
  integer :: signU, signT
  real(rp),dimension(2) :: r,s,k,p
  real(rp),dimension(2) :: r1,s1,k1,p1
  real(rp) :: denom,num1,num2,num2_1
  real(rp) :: t,u
  real(rp),parameter :: smax = 0.5
  logical :: int_u, int_t
  
  res = .false.
  r = [B(1)-A(1),B(2)-A(2)]
  s = [D(1)-C(1),D(2)-C(2)]
  if(norm2(r).gt.Epsilon) r1 = r/norm2(r)
  if(norm2(s).gt.Epsilon) s1 = s/norm2(s)
  denom = r(1)*s(2)-r(2)*s(1)
  k = [C(1)-A(1),C(2)-A(2)]
  p = [D(1)-A(1),D(2)-A(2)]
  if(norm2(k).gt.Epsilon) k1 = k/norm2(k)
  if(norm2(p).gt.Epsilon) p1 = p/norm2(p)
  if (abs(denom) > Epsilon) then
    num1 = k(1)*s(2)-k(2)*s(1)
    num2 = k(1)*r(2)-k(2)*r(1)
    t = num1/denom
    u = num2/denom
    num2_1 = k1(1)*r1(2)-k1(2)*r1(1)
    int_u = u.gt.Epsilon .and. u.lt.(1._rp-Epsilon)
    int_t = t.gt.Epsilon .and. t.lt.(1._rp-Epsilon)
    if(int_t .and. int_u)then
        res = .true.
    elseif((abs(t).lt.Epsilon .or. abs(t-1._rp).lt.Epsilon) .and. int_u)then
        res = .true.
    elseif((abs(u).lt.Epsilon .or. abs(u-1._rp).lt.Epsilon) .and. int_t)then
        res = .true.
    endif
  else
    num2 = k(1)*r(2)-k(2)*r(1)
    denom = r(1)*r(1)+r(2)*r(2)
    if (abs(num2) < Epsilon)then
      t = k(1)*r(1)+k(2)*r(2)
      t = t/denom
      u = (k(1)+s(1))*r(1)+(k(2)+s(2))*r(2)
      u = u/denom
!
      signU = int(sign(1._rp,u))    
      signT = int(sign(1._rp,t))
!
      if(t.gt.Epsilon .and. t.lt.1._rp-Epsilon)then
        res = .true.
      elseif(u.gt.Epsilon .and. u.lt.1._rp-Epsilon)then  
        res = .true.
      elseif(t.lt.Epsilon .and. u.gt.(1._rp-Epsilon))then
        res = .true.
      elseif(u.lt.Epsilon .and. t.gt.(1._rp-Epsilon))then
        res = .true.
      endif
!      res = t.gt.Epsilon .and. t.lt.(1._rp-Epsilon)
!      res = res .or. u.gt.Epsilon .and. u.lt.(1._rp-Epsilon)
!      res = res .or. t.lt.Epsilon .and. u.gt.(1._rp-Epsilon)
!      res = res .or. u.lt.Epsilon .and. t.gt.(1._rp-Epsilon)
      res = res .and. .not.(abs(t).lt.Epsilon .and. abs(1._rp-u).lt.Epsilon .or. abs(1._rp-t).lt.Epsilon .and. abs(u).lt.Epsilon)
!      res = res .and. .not.(abs(abs(t)-smax).lt.Epsilon .and. abs(abs(u)-smax).lt.Epsilon)
    endif
  endif
end subroutine segment_intersect_2D_2

subroutine write_intersection(tab,n_tab,t,ierror)
    
    !f2py integer*1, dimension(1000)                :: tab
    type(chaine_point_pt),dimension(:),intent(in)   :: tab              ! Table of the pointers toward the intersection points.
    integer,intent(inout)                           :: n_tab            ! Number of intersection curves.
    real(rp),intent(in)                             :: t                ! Current time.
    integer,intent(inout)                           :: ierror           ! Error flag.
    
    integer                                         :: j,k              ! Loop parameters.
    !f2py integer*1, dimension(1000)                :: ptr1,ptr2,ptr0
    type(chaine_point),pointer                      :: ptr1,ptr2,ptr0   ! Pointers.
    !f2py integer*1, dimension(1000)                :: P
    type(point)                                     :: P                ! Point
    character(len=20)                               :: num              ! Solution time
    real(rp),dimension(:,:),allocatable             :: TPoints          ! Table of points
    integer,dimension(:,:),allocatable              :: TEdges           ! Table of edges
    integer,parameter                               :: nPts = 1000      ! Number of points
    
    ! This subroutine writes the output file Intersections_curves.dat.
    
    ! Allocating
    allocate(TPoints(nPts*NBodies,3))
    allocate(TEdges(nPts*NBodies,6))
    
    ! Points
    k = 1
    do j=1,n_tab
        P = tab(j)%pt%val        
        TPoints(k,1:3) = P%coord
        TEdges(k,1) = P%nedge
        TEdges(k,2:6) = P%edge(1:5)
        k = k + 1
        ptr0 => tab(j)%pt
        ptr1 => ptr0
        ptr2 => ptr0%suiv
        do while(associated(ptr2) .and. .not. associated(ptr2,ptr0))
            
            ! One point
            P = ptr2%val            
            TPoints(k,1:3) = P%coord
            TEdges(k,1) = P%nedge
            TEdges(k,2:6) = P%edge(1:5)
            k = k + 1
            if(k.eq.nPts)then
                print*, "write_intersection: they are more than 1000 points on the intersections curves. Please change nPts."
                go to 9999
            end if
                        
            ! Next point
            if(.not.associated(ptr1,ptr2%suiv))then
                ptr1 => ptr2
                ptr2 => ptr2%suiv
            elseif(.not.associated(ptr1,ptr2%prec))then
                ptr1 => ptr2
                ptr2 => ptr2%prec
            else
                ierror = 100
                goto 9999
            endif
            
        enddo
    enddo
    
    9999 continue
    
    k = k - 1 ! To balance the last k = k + 1.
        
    ! Writting
    write( num, '( f0.4 )' ) t
    if(k.gt.1)then
        write(unit=ioIntersection,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =',k, ', E =',k-1,', ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    else
        write(unit=ioIntersection,fmt='(a,i,a,i,a)') 'Zone T = "'//trim(num)//'seconds", N =',k, ', E =',k,', ET=TRIANGLE, F=FEPOINT, STRANDID = 1, SOLUTIONTIME = '//trim(num)
    end if
    
    ! Nodes
    do j = 1,k
        write(ioIntersection,'(3f16.5,i5,5i3)') TPoints(j,1:3),TEdges(j,1),TEdges(j,2:6)
    end do
    ! Panels
    if(k.gt.1)then
        do j = 1,k-1
            write(ioIntersection,'(i5,i5,i5)') j,j+1,j+1
        end do
    else
        write(ioIntersection,'(i5,i5,i5)') 1,1,1
    end if
     
    ! Deallocating
    deallocate(TPoints)
    deallocate(TEdges)
    
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99  format('** error #',i3,' : ecriture fichier intersection')

end subroutine write_intersection

subroutine write_intersection_debug(tab,n_tab,ierror)
  
    !f2py integer*1, dimension(1000)                :: tab
    type(chaine_point_pt),dimension(:,:),intent(in) :: tab              ! Table of the pointers toward the intersection points.
    integer,intent(inout)                           :: n_tab            ! Number of intersection curves.
    integer,intent(inout)                           :: ierror           ! Error flag.

    integer                                         :: j,iter           ! Loop parameters.
    !f2py integer*1, dimension(1000)                :: ptr1,ptr2,ptr0
    type(chaine_point),pointer                      :: ptr1,ptr2,ptr0   ! Pointers.
    type(point)                                     :: P                ! Point.
    
    ! This subroutine writes the output file Intersections_curves_Debug.dat.
    
    ! Opening
    open(unit=ioIntersectionDebug,file='Intersection_curves_Debug.dat')
    write(ioIntersectionDebug,'(50a)') 'Title = "Debug - Intersection curves between the floaters and the free surface"'
    write(ioIntersectionDebug,'(50a)') 'VARIABLES = "X","Y","Z","n_tab","Nedge","E1","E2","E3","nface","Face 1","Face 2","Face 3"'
    
    ! Points
    do j=1,n_tab
        P = tab(j,1)%pt%val
        write(ioIntersectionDebug,'(3f16.5,9i3)') P%coord,j,P%nedge,P%edge(1:3),P%nface,P%face(1:3) !,eta,dz,P%nface,P%flag
        ptr0 => tab(j,1)%pt
        ptr1 => ptr0
        ptr2 => ptr0%suiv
        iter = 1
        do while(associated(ptr2).and. .not. associated(ptr2,ptr0) .and. iter.le.1000) ! 1000 is subjectif (exit condition).
            P = ptr2%val
            write(ioIntersectionDebug,'(3f16.5,9i3)') P%coord,j,P%nedge,P%edge(1:3),P%nface,P%face(1:3) !,eta,dz,P%nface,P%flag
            ptr1 => ptr2
            ptr2 => ptr2%suiv
            iter = iter + 1
        enddo
    enddo
    
    9999 continue
    if(ierror/=0)then
        write(*,99),ierror
    endif
    99  format('** error #',i3,' : ecriture fichier intersection debug')  
    
    ! Closing
    close(ioIntersectionDebug)
    
end subroutine write_intersection_debug

! -------------------------------------------------------------------------
! Check Border 1D
! -------------------------------------------------------------------------
! Domaine limite par vmin et vmax
subroutine check_border_1D(param,tbound,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  real(rp),dimension(1,2),intent(in) :: tbound
  real(rp),dimension(2),intent(in) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror

  ierror = 0
  is_border = .false.  
  
  if(param(2).gt.(tbound(1,2)-Epsilon) .or. param(2).lt.(tbound(1,1)+Epsilon))then
    is_border = .true.
  endif
!  
  if(is_border)then
    if (abs(param(2)-tbound(1,2)).lt.Epsilon) then 
       yout(1) = param(1)
       yout(2) = param(2)
    elseif(abs(param(2)-tbound(1,1)).lt.Epsilon) then
      yout(1) = param(1)
      yout(2) = param(2)
    elseif(param(2).lt.tbound(1,1) .and. abs(y(2)-tbound(1,1)).gt.Epsilon) then
      yout(2) = tbound(1,1)
      yout(1) = y(1)+(param(1)-y(1))*(param(2)-y(2))/(tbound(1,1)-y(2))
    elseif(param(2).lt.tbound(1,1) .and. abs(y(2)-tbound(1,1)).le.Epsilon) then
      yout(2) = y(2)
      yout(1) = y(1)
    elseif(param(2).gt.tbound(1,2) .and. abs(tbound(1,2)-y(2)).gt.Epsilon) then
      yout(2) = tbound(1,2)
      yout(1) = y(1)+(param(1)-y(1))*(param(2)-y(2))/(tbound(1,2)-y(2))
    elseif(param(2).gt.tbound(1,2) .and. abs(tbound(1,2)-y(2)).le.Epsilon) then
      yout(2) = y(2)
      yout(1) = y(1)
    endif
  endif  
!  
end subroutine check_border_1D  

! -------------------------------------------------------------------
! Check Border 2D (domaine rectangle)
! -------------------------------------------------------------------
subroutine check_border_2D(param,tbound,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  real(rp),dimension(2,2),intent(in) :: tbound
  real(rp),dimension(2) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror
! local
  real(rp) :: denom,s
  real(rp),dimension(2) :: A,B,C,D,P1,P2
  
  ierror = 0
  is_border = .false.
    
  A = [tbound(1,1),tbound(2,1)]
  B = [tbound(1,2),tbound(2,1)]
  C = [tbound(1,2),tbound(2,2)]
  D = [tbound(1,1),tbound(2,2)]
  
  P1 = [y(1),y(2)]
  P2 = [param(1),param(2)]
  
  if(param(1).lt.tbound(1,1) .or. param(1).gt.tbound(1,2) .or. param(2).lt.tbound(2,1) .or. param(2).gt.tbound(2,2))then
    if((abs(y(1)-tbound(1,1)).lt.Epsilon .or. abs(y(1)-tbound(1,2)).lt.Epsilon).or.&
    &  (abs(y(2)-tbound(2,1)).lt.Epsilon .or. abs(y(2)-tbound(2,2)).lt.Epsilon))then
      yout = y
      is_border = .true.
    endif
  endif
  
  if(.not.is_border)then
    if((abs(param(1)-tbound(1,1)).lt.Epsilon .or. abs(param(1)-tbound(1,2)).lt.Epsilon).or.&
    &  (abs(param(2)-tbound(2,1)).lt.Epsilon .or. abs(param(2)-tbound(2,2)).lt.Epsilon))then
      yout = param
      is_border = .true.
    endif
  endif
  
  if(.not.is_border)then
    call segment_intersect_2D(P1,P2,A,B,is_border)
    if(is_border)then
      yout(2) = tbound(2,1)
      denom = P1(2)-P2(2)
      if(abs(denom).gt.Epsilon)then
        s = (P1(2)-tbound(2,1))/denom
        yout(1) = s*P2(1)+(1._rp-s)*P1(1)
      else
        ierror = 100
        goto 9999
      endif
    endif
  endif 
  
  if(.not.is_border)then
    call segment_intersect_2D(P1,P2,B,C,is_border)
    if(is_border)then
      yout(1) = tbound(1,2)
      denom = P1(1)-P2(1)
      if(abs(denom).gt.Epsilon)then
        s = (P1(1)-tbound(1,2))/denom
        yout(2) = s*P2(2)+(1._rp-s)*P1(2)
      else
        ierror = 110
        goto 9999
      endif
    endif
  endif
  
  if(.not.is_border)then  
    call segment_intersect_2D(P1,P2,C,D,is_border)
    if(is_border)then
      yout(2) = tbound(2,2)
      denom = P1(2)-P2(2)
      if(abs(denom).gt.Epsilon)then
        s = (P1(2)-tbound(2,2))/denom
        yout(1) = s*P2(1)+(1._rp-s)*P1(1)
      else
        ierror = 120
        goto 9999
      endif
    endif
  endif
  
  if(.not.is_border)then
    call segment_intersect_2D(P1,P2,D,A,is_border)
    if(is_border)then
      yout(1) = tbound(1,1)
      denom = P1(1)-P2(1)
      if(abs(denom).gt.Epsilon)then
        s = (P1(1)-tbound(1,1))/denom
        yout(2) = s*P2(2)+(1._rp-s)*P1(2)
      else
        ierror = 130
        goto 9999
      endif
    endif
  endif
  
9999 continue
  if(ierror/=0)then 
    write(*,99),ierror
  endif
99 format('** error #',i3,' : pb. traitement point proche paroi (plan)')     
  
end subroutine check_border_2D

subroutine check_border_ND(param,tbound,nbound,y,yout,is_border,ierror)
  implicit none
  real(rp),dimension(2),intent(in) :: param
  real(rp),dimension(:,:),intent(in) :: tbound
  integer,intent(in) :: nbound
  real(rp),dimension(2) :: y
  real(rp),dimension(2),intent(inout) :: yout
  logical,intent(inout) :: is_border
  integer,intent(inout) :: ierror
!
  select case(nbound)
  case(1)
    call check_border_1D(param,tbound,y,yout,is_border,ierror)
  case(2)
    call check_border_2D(param,tbound,y,yout,is_border,ierror)
  case default
    ierror = 100
    goto 9999
  end select
  
9999 continue
  if(ierror/=0)then
    write(*,99),ierror
  endif
99 format("** error #",i3," : choix dimension incorrect")      

end subroutine check_border_ND

subroutine adapt_repere(rep0,rep1,N,C)
    implicit none
    !f2py integer*1, dimension(1000) :: rep0
    type(repere3d),intent(in) :: rep0
    !f2py integer*1, dimension(1000) :: rep1
    type(repere3d),intent(inout) :: rep1
    real(rp),dimension(3),intent(in) :: N
    real(rp),dimension(3),intent(in) :: C
!   local
    real(rp) :: R,inv_R
    real(rp),dimension(3,3) :: A

    A = 0._rp

    R = sqrt(N(1)*N(1)+N(2)*N(2))
    if(abs(R).gt.Epsilon)then
        inv_R = 1._rp/R
        A(1:3,1) = [N(1)*N(3)*inv_R , -N(2)*inv_R , N(1)]
        A(1:3,2) = [N(2)*N(3)*inv_R ,  N(1)*inv_R , N(2)]
        A(1:3,3) = [-R , 0._rp, N(3)]
    else
        A(1,1) = sign(1._rp,N(3))
        A(2,2) = 1._rp
        A(3,3) = sign(1._rp,N(3))
    endif

    rep1%e3 = N
    rep1%e2 = matmul(A,rep0%e2)
    rep1%e1 = matmul(A,rep0%e1)

    rep1%origine = C

end subroutine adapt_repere

subroutine assign_cercle(geom,rep,C,normale,r,nline,ierror,P1,P2)
    implicit none
    !f2py integer*1, dimension(1000) :: geom
    type(type_geom),intent(inout) :: geom
    !f2py integer*1, dimension(1000) :: rep
    type(repere3d),intent(in) :: rep
    real(rp),dimension(3),intent(in) :: C
    real(rp),dimension(3),intent(in) :: normale
    real(rp),intent(in) :: r
    integer,intent(inout) :: nline
    integer,intent(inout) :: ierror
    !f2py integer*1, dimension(1000) :: P1,P2
    type(point),optional :: P1,P2
!   local
    integer :: k
    !f2py integer*1, dimension(1000) :: cercle
    type(GCercle) :: cercle
    real(rp),dimension(3) :: E1,E2,E3
    !f2py integer*1, dimension(1000) :: M,N
    type(point) :: M,N

    ierror = 0

    nline = nline+1
   
    if(present(P1))then
        cercle%repere%origine = C
        E3 = normale
        E1 = (P1%coord-C)
        E1 = E1/norm2(E1)
        call Computation_vect_product(E3,E1,E2)
        E2 = E2/norm2(E2)
        cercle%repere%e1 = E1
        cercle%repere%e2 = E2
        cercle%repere%e3 = E3
    elseif(Symmetry .and. abs(C(2)).lt.Epsilon .and. abs(normale(2)).lt.Epsilon)then
        cercle%repere%origine = C
        E3 = normale
        E2 = [0._rp,1._rp,0._rp]
        call Computation_vect_product(E2,E3,E1)
        E1 = E1/norm2(E1)
        cercle%repere%e1 = E1
        cercle%repere%e2 = E2
        cercle%repere%e3 = E3
    else
        call adapt_repere(rep,cercle%repere,normale,C)
    endif

    cercle%rayon = r
    cercle%index = nline

    if(present(P1) .and. present(P2))then
        call cart2loc(P1%coord,cercle%repere,M%coord)
        call cart2loc(P2%coord,cercle%repere,N%coord)
        cercle%P1 = M
        cercle%P2 = N
        call cart2cer(N%coord,cercle,cercle%theta)
        !if(abs(N%coord(1)).gt.Epsilon)then
        !    theta = atan(N%coord(2)/N%coord(1))
        !else
        !    theta = sign(0.5_rp*PI,N%coord(2))
        !endif
        !if(M%coord(1).lt.0) theta = theta+PI
        !if(theta.lt.0)      theta = theta+TWOPI
        !cercle%theta = theta
    elseif(Symmetry .and. abs(C(2)).lt.Epsilon .and. abs(normale(2)).lt.Epsilon)then
        cercle%P1 = r*[1._rp,0._rp,0._rp]
        cercle%P2 = -r*[1._rp,0._rp,0._rp]
        cercle%theta = PI
    else
        cercle%P1 = r*[1._rp,0._rp,0._rp]
        cercle%P2 = cercle%P1
        cercle%theta = TWOPI
    endif


    k = geom%ncercle+1
    geom%cercle(k) = cercle
    geom%ncercle = k

end subroutine assign_cercle

subroutine add_point_face(P,iface,nelem)
    implicit none
    !f2py integer*1, dimension(1000) :: P
    type(point),intent(inout) :: P
    integer,dimension(*),intent(in) :: iface
    integer,intent(in) :: nelem
!   local
    integer :: nface
    nface = P%nface
    P%face(nface+1:nface+nelem) = iface(1:nelem)
    P%nface = nface+nelem
end subroutine add_point_face

end module GeomFonct
