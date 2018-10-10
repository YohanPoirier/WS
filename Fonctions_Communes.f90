module FonctionsCommunes
!!!!! Contient
!       o vect_product(A,B) :       Produit vectoriel --> à remplacer par une fonction de LAPACK à terme
!       o Crampe(x,L,a) :           Fonction Rampe
!       o GMRES :                   Solveur iteratif
!       o LU :                      Solveur direct par décomposition LU
!!!!!
use Constantes
use Parameters
use Structuresdonnees
use iso_c_binding
implicit none

contains




! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
subroutine inv_matrice(A, Ainv,N)

    integer, intent(in)                   :: N
    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: Ainv

    real(dp), dimension(:), allocatable   :: work  ! work array for LAPACK
    integer, dimension(:), allocatable    :: ipiv   ! pivot indices
    integer :: info

    
    allocate(work(N), ipiv(N))


    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A


    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
    
        write(1111,*) 'Matrix is numerically singular!'
        stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
        write(1111,*) 'Matrix inv failed...!'
        stop 'Matrix inversion failed!'
    end if
    
    
    deallocate(work, ipiv)
    
end subroutine inv_matrice




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                    Produit Vectoriel                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_vect_product(A,B,vect_product)
    !!!!! Problème
    !
    ! Calcul le produit vectoriel entre deux vecteurs de dimensions 3
    !
    !!!!!
    ! Paramètres
    real(rp), dimension(3), intent(in) :: A,B
    real(rp), dimension(3),intent(out) :: vect_product
    vect_product(1)=A(2)*B(3)-A(3)*B(2)
    vect_product(2)=A(3)*B(1)-A(1)*B(3)
    vect_product(3)=A(1)*B(2)-A(2)*B(1)

end subroutine Computation_vect_product 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   function Crampe                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_Crampe2(x,x0,x1,Crampe,precision)
    !!!!! Problème : Fonction "rampe"
    ! Fonction rampe basée sur la fonction tanh, centrée sur xm et étirée d'une longueur L
    !            <-L->   
    !     1 +         _____
    !       |        /
    !       |       /
    !       |      /
    !       |     /
    !       |____/
    !     0 +----+-+-+------
    !          x0 xm x1
    !!!!!
    ! Paramètres :
    real(rp),intent(in) :: x, x0,x1
    real(rp),intent(out) :: Crampe
    real(rp), optional :: precision
    ! Variables Locales :
    real(rp) :: t, xm, L, prec

    prec = 0.95_RP
    if (present(precision)) prec = precision
    L = (x1-x0)/(2._RP*atanh(prec))
    xm = 0.5_RP*(x0+x1)
    if(abs(L).gt.Epsilon)then
        t=(x-xm)/L
        !t=(x-L)**a
        Crampe = (1+tanh(t))*0.5_RP
    else
        Crampe = 1._rp
    endif

end subroutine Computation_Crampe2

subroutine Computation_Crampe(x,T1,T2,val)
      real(rp),intent(in) :: x
      real(rp),intent(in) :: T1,T2
      real(rp),intent(out) :: val

      if (x-T1 .le.Epsilon) then
      val = 0._RP
  elseif (x-T2 .ge.Epsilon) then
      val = 1._RP
  else
       val = 0.5_RP*(1._RP-cos(pi*(x-T1)/(T2-T1)))
  end if

end subroutine Computation_Crampe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   function Cramp                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Computation_Cramp(x,L,a,Cramp)
    !!!!! Problème : Fonction "rampe"
    ! Fonction rampe basée sur la fonction tanh, centrée sur L et étirée d'un facteur a
    !            <-a->   
    !     1 +         _____
    !       |        /
    !       |       /
    !       |      /
    !       |     /
    !       |____/
    !     0 +------+--------
    !              L 
    !!!!!
    ! Paramètres :
    real(rp),intent(in) :: x, L, a
    real(rp),intent(out) :: Cramp
    ! Variables Locales :
    real(rp) :: t

    t=x**a-L**a
    Cramp = (1+tanh(t))*0.5_RP

end subroutine Computation_Cramp

subroutine Ramp_Management(t,TimeRamp)

    real(rp),intent(in) :: t        ! Current time.
    real(rp),intent(out):: TimeRamp ! Ramp parameter.

    ! This subroutine calls the subroutine CRampe if necessary.

    if (t.ge.T2) then
        TimeRamp = 1._RP
    else
        call Computation_CRampe(t,T1,T2,TimeRamp)
    end if

end subroutine Ramp_Management

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                       GMRES                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GMRES(A, B, Sol, Nn, ierror)
    
    real(rp), dimension(:,:), intent(in)    :: A                                    ! A matrix.
    real(rp), dimension(:), intent(in)      :: B                                    ! Rhs.
    real(rP), dimension(:), intent(out)     :: Sol                                  ! Solution.
    integer, intent(in)                     :: Nn                                   ! Dimension of the linear system.
    integer,intent(inout)                   :: ierror                               ! Error flag.
    
    real(rp), dimension(Nn)                 :: x                                    !
    real(rp), allocatable                   :: tmp(:)                               !
    integer                                 :: RCI_request, itercount, niter = 301  ! 
    integer, dimension(128)                 :: ipar                                 !
    real(rp), dimension(128)                :: dpar                                 !
    
    ! This subroutine solves the linear system A*Sol = B from a GMRES algorithm for a full A matrix.
    
    ierror = 0
    allocate(tmp(((2*niter+1)*Nn+niter*niter+9)/2 + 1))
    
    ! Initialisation de GMRES
    call dfgmres_init(Nn, x, B, RCI_request, ipar, dpar, tmp) 
    ipar(1) = Nn
    ipar(5) = niter ! Nombre maxi d'itérations
    ipar(12) = 1  ! test sur la norm2 du vecteur (?): 0 --> test à faire par l'utilisateur via RCI_request=4, 1 --> test auto par dfgmres 
    x = Sol ! Initialisation de la solution : B --> Second Membre; Sol --> solution du pas de temps précédent.
    
    ! Vérification des paramètres de dfgmres
    call dfgmres_check(Nn, x, B, RCI_request, ipar, dpar, tmp) 
    if (RCI_request/=0) then 
        print*, ' Pb lors de l''init de gmres '
        ierror = 100
        goto 9999
    endif
    
    1   call dfgmres(Nn, x, B, RCI_request, ipar, dpar, tmp)
    if (RCI_request.eq.0) then
        go to 2 ! sortie de la boucle
    elseif (RCI_request.eq.1) then ! conditionnement de la solution avec A
        tmp(ipar(23):ipar(23)+Nn-1)=matmul(A,tmp(ipar(22):ipar(22)+Nn-1))   
        go to 1
    elseif (RCI_request.eq.2) then
        if (dpar(5).gt.TolGMRES) then ! test sur la précision
            go to 1 ! reste dans la boucle
        else
            go to 2 ! sortie de la boucle
        endif
    end if
    
    2  call dfgmres_get(Nn, x, B, RCI_request, ipar, dpar, tmp, itercount) ! récupération de la solution
    if (itercount.eq.niter)then
        print*, 'GMRES : no convergence', itercount
        ierror = 200
        goto 9999
    endif
    
    Sol = x
    
    deallocate(tmp)
    9999 continue
    if(ierror/=0)then
        print*,'error GMRES'
    endif
    
end subroutine GMRES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       Solveur LU                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LU(A, B, Sol, N, ierror,Saving_A)
    
    real(rp), dimension(:,:), intent(inout)     :: A                    ! Matrix A : by default it not the same in input and in output.
    real(rp), dimension(:), intent(inout)       :: B                    ! Vector B.
    real(rP), dimension(:), intent(out)         :: Sol                  ! Solution.
    integer, intent(in)                         :: N                    ! Size of A, B and Sol.
    integer,optional                            :: ierror               ! Error flag.
    integer,optional                            :: Saving_A             ! If Saving_A is present and = 1 the matrix A in output is the same as in input, OTHERWISE IT IS NOT THE CASE !!!
    
    character(len=1)                            :: trans
    integer                                     :: info, lda, ldb, nrhs ! Leading dimensions, number of rhs.
    integer, dimension(N)                       :: ipiv                 ! Pivot indices.
    real(rp), dimension(N)                      :: B1                   ! Local copy of B.
    real(rp), dimension(:,:),allocatable        :: A1                   ! Local copy of A.
    integer                                     :: ierr                 ! Error flag.
    
    ! This subroutine solves the linear system A*Sol = B from a LU factorization for a full A matrix.
    
    ierr = 0
    
    ! Saving A
    if(present(Saving_A) .and. Saving_A.eq.1)then
        allocate(A1(N,N))
        A1 = A
    end if
    
    ! Saving B
    B1 = B 
    
    lda = N ! Leading dimension de A
    ldb = lda ! Leading dimension de B
    if (RP==SP) then ! Test sur le type de réel: simple ou double
        ipiv = 0 
        trans = 'n'
        nrhs = 1 ! Nombre de second membre
    
        ! LU factorization of A.
        call sgetrf(N,N,A,lda, ipiv,info)
        if(info.ne.0) then
            ierr = 100 ;  goto 9999
        endif

        ! Résolution du système à partir de la décomposition Lu.
        call sgetrs(trans,N,nrhs,A,lda,ipiv,B,ldb,info)
        if(info.ne.0) then
            ierr = 200 ; goto 9999
        endif
        
    else
        ipiv = 0 
        trans = 'n'
        nrhs = 1 ! Nombre de second membre.
        
        ! LU factorization of A
        call dgetrf(N,N,A,lda, ipiv,info)
        if(info.ne.0) then
            ierr = 300 ; goto 9999 
        endif
        
        ! Résolution du système à partir de la décomposition Lu.
        call dgetrs(trans,N,nrhs,A,lda,ipiv,B,ldb,info)
        if(info.ne.0) then
            ierr = 400 ; goto 9999
        endif
        
    end if 
    
    Sol = B
    B = B1
    if(present(Saving_A) .and. Saving_A.eq.1)then
        A = A1
        deallocate(A1)
    end if
    
    9999 continue

    if(info/=0)then
        print*, '       info =', info, N
    endif

    if(ierr/=0)then
        write(*,90),ierr
        
        ! Error in dgetrf.
        if(ierr.eq.300)then
            if(info.gt.0)then
                print*,"A(",info,",",info,") is 0."
                print*,A(info,info)
            else
                print*,"A(",info,",",info,") has an illegal value."
                print*,A(info,info)
            end if
        end if
        
        ! Error in dgetrs.
        if(ierr.eq.400)then
            if(info.gt.0)then
                print*,"A(",info,",",info,") has an illegal value."
                print*,A(info,info)
            end if
        end if
            
        pause
    endif
    90 format('error #',i3,' : SUB_LU')
    if(present(ierror))then
        ierror = ierr
    endif
    
end subroutine LU

subroutine InvMat(A,N,invA,ierror)
    
    real(rp), dimension(:,:), intent(inout) :: A            ! Matrix A to inverse
    integer, intent(in)                     :: N            ! Size of A.
    real(rp), dimension(N,N), intent(inout) :: invA         ! Inverse of A.
    integer,optional                        :: ierror       ! Error flag.
    
    character(len=1)                        :: trans
    integer                                 :: info, lda    ! Leading dimensions, number of rhs.
    integer, dimension(N)                   :: ipiv         ! Pivot indices.
    real(rp), dimension(N)                  :: work         ! Work array for LAPACK.
    integer                                 :: ierr         ! Error flag.
    
    ! This subroutine compute the inverse of the fulle matrix A from a LU factorization.
    
    ierr = 0
    
    ipiv = 0 
    trans = 'n'
    lda = N ! Leading dimension de A
    
    ! Copy of A in invA.
    invA = A
    
    ! LU factorization of A
    call dgetrf(N,N,invA,lda, ipiv,info)
    if(info.ne.0) then
        ierr = 300 ; goto 9999 
    endif
    
    ! Résolution du système à partir de la décomposition Lu
    call dgetri(N,invA,lda,ipiv,work,N,info)
    if(info.ne.0) then
        ierr = 400 ; goto 9999
    endif
    
    9999 continue

    if(info/=0)then
        print*, '       info =', info, N
    endif

    if(ierr/=0)then
        write(*,90),ierr
        pause
    endif
    90 format('error #',i3,' : InVMat')
    if(present(ierror))then
        ierror = ierr
    endif
    
end subroutine InVMat

subroutine pardiso_solver(ja,ia,aa,nNonZero,ndim,size_ia,B,Sol,ierror)
    
    integer,intent(in)                          :: nNonZero ! Number the non-zero elements of A.
    integer,intent(in)                          :: ndim     ! Size of the linear system.
    integer,intent(in)                          :: size_ia  ! ndim + 1.
    integer,dimension(nNonZero),intent(in)      :: ja       ! Contains column indices of the sparse matrix A (CSR3).
    integer,dimension(size_ia),intent(in)       :: ia       ! ia(i) points to the first column index of row i in the array ja (CSR3).
    real(rp),dimension(nNonZero),intent(in)     :: aa       ! Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja.
    real(rp),dimension(ndim),intent(inout)      :: B        ! Rhs of the linear system.
    real(rp),dimension(ndim),intent(inout)      :: Sol      ! Solution of the linear system.
    integer,intent(inout)                       :: ierror   ! Error flag.
    
    INTEGER*8                                   :: pt(64)       ! Handle to internal data structure.
    integer                                     :: phase    ! Controls the execution of the Pardiso solver.
    integer                                     :: maxfct   ! Maximum number of factors with identical sparsity structure that must be kept in memory at the same time (?).
    integer                                     :: mnum     ! Indicates the actual matrix for the solution phase (?).
    integer                                     :: mtype    ! Defines the matrix type, which influences the pivoting method.
    integer                                     :: nrhs     ! Number of right-hand sides that need to be solved for.
    integer,dimension(64)                       :: iparm    ! Parameters of the subroutine Paradiso.
    integer                                     :: msglvl   ! Message level information for Pardiso.
    integer                                     :: idum     ! Dummy integer.
    real(rp)                                    :: ddum     ! Dummy real.
    integer                                     :: j        ! Loop parameter.
    
    ! This subroutine solves the linear system A*Sol = B from a LU factorization for a sparse A matrix.
    
    ! Initialization of the internal solver memory. This is only necessary for the FIRST call of the PARDISO solver.
    do j = 1,64
       pt(j) =  0
    end do 
    
    ! Manually set thread number.
    call mkl_set_dynamic(0) ! Enables MKL to dynamically change the number of OpenMP threads
    call mkl_set_num_threads(Nthreads) ! Specifies the number of OpenMP threads to use.
    
    ! Parameters for Pardiso.
    maxfct = 1
    mnum = 1
    mtype = 11 ! Real and nonsymmetric matrix.
    nrhs = 1
    msglvl = 0
    ddum = 0._RP
    idum = 0
    
    iparm = 0
    iparm(1) = 1 ! No solver default.
    iparm(2) = 3 ! fill-in reordering from METIS (?).
    ! iparm(3) must be 0.
    iparm(4) = 0 ! No iterative-direct algorithm (?).
    iparm(5) = 0 ! No user fill-in reducing permutation.
    iparm(6) = 0 ! =0 solution on the first n compoments of x otherwise (=1) solution stores in the rhs B.
    ! iparm(7) is used as output.
    iparm(8) = 1 ! Numbers of iterative refinement steps.
    ! iparm(9) must be 0.
    iparm(10) = 13 ! Perturbe the pivot elements with 1E-13.
    iparm(11) = 1 ! Use nonsymmetric permutation and scaling MPS. Default for nonsymmetric matrices.
    iparm(12) = 0 ! Solve a linear system AX=B.
    iparm(13) = 1 ! Maximum weighted matching algorithm is switched-off (default for nonsymmetric)
    iparm(18) = -1 ! Output: number of nonzeros in the factor LU.
    iparm(19) = -1 ! Output: Mflops for LU factorization.
    iparm(27) = 1 ! check the integer arrays ia and ja.
    
    ! Analysis.
    phase = 11
    call pardiso(pt, maxfct, mnum, mtype, phase, ndim, aa, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
    if(ierror.ne.0)then
        print*,"create_dref_grid: Pardiso error in the analysis: ",ierror
    end if
        
    ! Factorization.
    phase = 22
    call pardiso(pt, maxfct, mnum, mtype, phase, ndim, aa, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
    if(ierror.ne.0)then
        print*,"create_dref_grid: Pardiso error in the factorization: ",ierror
    end if
        
    ! Back substitution and iterative refinement.
    phase = 33
    call pardiso(pt, maxfct, mnum, mtype, phase, ndim, aa, ia, ja,idum, nrhs, iparm, msglvl, B, Sol, ierror)
    
    if(ierror.ne.0)then
        print*,"create_dref_grid: Pardiso error in the back substitution: ",ierror
    end if
    
    ! Release the memory.
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, ndim, ddum, idum, idum,idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
        
    if(ierror.ne.0)then
        print*,"pardiso_solver: ndim = ",ndim
        print*,"pardiso_solver: nNonZero = ",nNonZero
        if(ierror.eq.-1) print*,"pardiso_solver: input inconsistent."
        if(ierror.eq.-2)then 
            print*,"pardiso_solver: not enough memory."
            call exit() ! Stop the executable and run it again from State.dat.
        end if
        if(ierror.eq.-3) print*,"pardiso_solver: reordering problem."
        if(ierror.eq.-4) print*,"pardiso_solver: zero pivot, numerical factorization or iterative refinement problem."
        if(ierror.eq.-5) print*,"pardiso_solver: unclassified (internal) error."
        if(ierror.eq.-6) print*,"pardiso_solver: reordering failed (matrix types 11 and 13 only)."
        if(ierror.eq.-7) print*,"pardiso_solver: diagonal matrix is singular."
        if(ierror.eq.-8) print*,"pardiso_solver: 32-bit integer overflow problem."
        if(ierror.eq.-9) print*,"pardiso_solver: not enough memory for OOC."
        if(ierror.eq.-10) print*,"pardiso_solver: error opening OOC files."
        if(ierror.eq.-11) print*,"pardiso_solver: read/write error with OOC files."
        if(ierror.eq.-12) print*,"pardiso_solver: (pardiso_64 only) pardiso_64 called from 32-bit library."
    end if
    
end subroutine pardiso_solver

subroutine find_konde(konde,w,profondeur,Nhoule)
! Parameters
real(rp), dimension(:) :: konde, w 
real(rp) :: profondeur
integer :: Nhoule
! Variables
integer :: j, jn
real(rp) :: Err
real(rp), dimension(1000) :: ktemp
! Begin
do jn=1,Nhoule
    j = 1
    ktemp(j) = w(jn)**2/g
    j=2
    ktemp(j) = w(jn)**2/g/tanh(ktemp(j-1)*profondeur)
    Err = abs(ktemp(j)-ktemp(j-1))/abs(ktemp(j))
    do while (Err.gt.Epsilon)
        j=j+1
        ktemp(j) = w(jn)**2/g/tanh(ktemp(j-1)*profondeur)
        Err = abs(ktemp(j)-ktemp(j-1))/abs(ktemp(j))
    end do
    konde(jn) = ktemp(j)
end do
! End
end subroutine find_konde

subroutine Serie_Geom(n,Sn,a,q)
! Determine la raison q, d'une suite geometrique, de terme initial a, de somme partielle Sn, pour n termes :
! Sn = a.Sum_{k=0,n} q^k = a.(1-q^(n+1))/(1-q)
! Parameters
integer, intent(in) :: n
real(rp), intent(in) :: Sn, a
real(rp), intent(out) :: q
! Locals
integer :: j, jmax
real(rp) :: S1, S2
! Begin
jmax = 10000
if (abs(Sn).lt.Epsilon) then
    print*, 'Serie_Geom : Sn nul'
    q = 0._RP
    return
end if

if (abs(q).lt.Epsilon) q = a
55 S1 = Sn - a*(1._RP-q**(n+1))/(1._RP-q)
S2 = a*((n+1)*(1-q)*q**n + (1-q)**(n+1))/(1-q)**2
j = 1
do while (not(abs(S1).lt.Epsilon .or. j.gt.jmax))
    j = j+1
    q = q - S1/S2
    S1 = Sn - a*(1._RP-q**(n+1))/(1._RP-q)
    S2 = a*((n+1)*(1-q)*q**n + (1-q)**(n+1))/(1-q)**2
end do
if (j.gt.jmax) print*, 'Error Serie_Geom : ', j, q, S1
if (q.lt.0) then
    q = -q
    go to 55
end if
if (.false.) then
    S1 = a*(1._RP-q**(n+1))/(1._RP-q)
    print*, q, S1, Sn, j
end if
! End
end subroutine Serie_Geom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     Test Vecteurs Entiers                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine common_int(a,na,b,nb,ires,nres,ierror)
    
    integer,dimension(*),intent(in)     :: a,b      ! Vectors to compare
    integer,intent(in)                  :: na,nb    ! Length of a and b
    integer,dimension(*),intent(inout)  :: ires     ! Coordinates in common
    integer,intent(inout)               :: nres     ! Number of coordinates in common
    integer,intent(inout),optional      :: ierror   ! Error flag

    integer                             :: j,k      ! Loop parameters
        
    ! This subroutine tests if the integer vectors a and b have some coordinates in common.
    
    nres = 0
    
    do j=1,na
        do k=1,nb 
            if(a(j)==b(k))then
            nres = nres+1
            ires(nres) = a(j)
            endif
        enddo
    enddo

end subroutine common_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine entire_int(a,na,b,nb,ires,nres,ierror)
  implicit none
  integer,dimension(*),intent(in) :: a,b
  integer,intent(in) :: na,nb
  integer,dimension(*),intent(inout) :: ires
  integer,intent(inout) :: nres
  integer,intent(inout) :: ierror
! local
  integer :: j,n
  integer,dimension(1) :: aux
  integer :: naux
!  
  n=max(na,nb)
  if(n.le.0)then
    ierror = 100
    goto 9999
  endif
!  
  if(na.ge.nb)then
    ires(1:na) = a(1:na)
    nres = na
    do j=1,nb
      call common_int(ires,na,b(j),1,aux,naux)
      if(naux==0)then
        nres = nres+1
        ires(nres) = b(j)
      endif
    enddo
  elseif(na.lt.nb)then
    ires(1:nb) = b(1:nb)
    nres = nb
    do j=1,na
      call common_int(ires,nb,a(j),1,aux,naux)
      if(naux==0)then
        nres = nres+1
        ires(nres) = a(j)
      endif
    enddo
  endif
!  
9999 continue
  if(ierror/=0)then
    write(*,99),ierror
  endif  
99 format('** error #',i3,' : entire_int')  
end subroutine entire_int   

subroutine find_voisin(grid,val,res)
  
    real(rp),dimension(:),intent(in)    :: grid
    real(rp),intent(in)                 :: val
    integer,dimension(2),intent(out)    :: res
    
    integer                             :: n,n1,n2,nmax
    logical                             :: is_find
    integer                             :: p

    n=0
    is_find = .false.
 
    n1 = 1
    n2 = size(grid)
    nmax = 2*n2
  
    do while(.not. is_find .and. n<=nmax)
        n=n+1
        
        if((n2-n1)==1)then 
            res(1) = n1
            res(2) = n2
            is_find = .true.
        else 
            p = n1+(n2-n1)/2 ! integer !
            if(abs(val-grid(p))<=Epsilon)then
                res = [p,p]
                is_find = .true.
            elseif(val > grid(p))then
                n1 = p
            elseif(val < grid(p))then
                n2 = p
            endif
        endif
        
    enddo
    
end subroutine find_voisin

subroutine find_pos(grid,val,res)
  implicit none
  real(rp),dimension(:),intent(in) :: grid
  real(rp),intent(in) :: val
  integer,intent(out) :: res
! local
  integer,dimension(2) :: ivoisin
  real(rp) :: x1,x2
  
  call find_voisin(grid,val,ivoisin)
  
  x1 = abs(val-grid(ivoisin(1)))
  x2 = abs(val-grid(ivoisin(2)))
  
  if(x1 .lt. x2)then
    res = ivoisin(1)
  else
    res = ivoisin(2)
  endif  
  
end subroutine find_pos

RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
REAL(RP), DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort

subroutine sort_order(list,order)
    implicit none
    real(rp),dimension(:),intent(inout) :: list
    integer,dimension(:),intent(in) :: order
!   local
    integer :: j
    real(rp),dimension(:),allocatable :: temp

    allocate(temp(size(list)))
    
    temp(:) = list(:)
    
    do j=1,size(list)
        list(j) = temp(order(j))
    enddo

    deallocate(temp)

end subroutine sort_order

subroutine norm_riemann(mat,vect,val)
    
    real(rp),dimension(2,2),intent(in) :: mat
    real(rp),dimension(2),intent(in) :: vect
    real(rp),intent(out) :: val
    
    real(rp),dimension(2) :: vect_temp
    
    vect_temp = matmul(mat,vect)
    val = dot_product(vect,vect_temp)
    if(abs(val).lt.Epsilon)then
        val = 0._rp
    elseif(val.gt.Epsilon)then
        val = sqrt(val)
    else
        print*,'** error : cannot calculate norm2 riemann'
        val = -999999
    endif
    
end subroutine norm_riemann
    
subroutine LissageGaussien(Maillage,Ecoulement,ierror)
    
    !f2py integer*1, dimension(1000):: Maillage
    type(TMaillage),intent(in)      :: Maillage                                 ! Mesh
    !f2py integer*1, dimension(1000):: Ecoulement
    type(TEcoulement),intent(inout) :: Ecoulement                               ! Flow parameters (phi, etc.)
    integer,intent(inout)           :: ierror                                   ! Error flag.
    
    integer                         :: j, k, kj, NVoisin                        ! Loop parameters.
    integer,dimension(NPVMax,2)     :: Tvoisin                                  ! Neighbour point.
    real(rp)                        :: lj, poids, cl, invsqrt2Pi, sigma, dmax   ! 
    real(rp), dimension(3)          :: M, Mj                                    ! Points.
    real(rp), allocatable           :: Eta(:)                                   ! Updated free surface elevation.
    
    ! This subroutine smoothes the free surface elevation.
    
    ierror = 0
    invsqrt2Pi = 1._RP/sqrt(2*PI)
    allocate(Eta(Maillage%FS%IndFS(1):Maillage%FS%IndFS(3)))
    Eta = 0._rp
    do j=Maillage%FS%IndFS(1),Maillage%FS%IndFS(3)
        sigma = 0.5_RP*sqrt(Maillage%Tnoeud(j)%Aire)*(1._RP + Maillage%Tnoeud(j)%Damping)
        poids = 0._rp
        M = Maillage%Tnoeud(j)%Pnoeud
        NVoisin = Maillage%Tnoeud(j)%NVoisin(2)
        Tvoisin(1:NVoisin,1) = Maillage%Tnoeud(j)%Tvoisin(1:Nvoisin,1)
        dmax = 0._RP
        do k=1,NVoisin
            kj = Tvoisin(k,1)
            Mj = Maillage%Tnoeud(abs(kj))%Pnoeud
            if (kj.lt.0) Mj(2) = -Mj(2)
            lj = exp(-dot_product(Mj-M,Mj-M)/(2._RP*sigma**2))*invsqrt2Pi/sigma
            Eta(j) = Eta(j) + Ecoulement%Eta(abs(kj))%perturbation * lj
            poids = poids + lj
            dmax = max(dmax,norm2(Mj-M))
        end do
        
        if (dmax.lt.sigma) print*, 'Warning Filter : dmax = ', dmax, '< sigma = ', sigma
        
        if(poids.gt.Epsilon)then
            Eta(j) = Eta(j)/poids
        else
            ierror = 100
            goto 9999
        endif
        
    end do
    
    ! Applying the filtering to the free surface elevation.
    do j=Maillage%FS%IndFS(1),Maillage%FS%IndFS(3)
        Ecoulement%Eta(j)%perturbation = Eta(j)
    end do
      
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
    deallocate(Eta)
    
end subroutine LissageGaussien

subroutine CopyHoule(HouleObj, HouleInit)
    
    !f2py integer*1, dimension(1000):: HouleObj
    type(THoule), intent(in)        :: HouleObj     ! Wave structure.
    !f2py integer*1, dimension(1000):: HouleInit
    type(THoule), intent(out)       :: HouleInit    ! Wave structure.
    
    integer                         :: j            ! Loop parameter.
    
    ! This subroutine copies a structure THoule into another one.
    
    ! Stockage des paramètres de houle originaux
    HouleInit%Htype = Htype
    HouleInit%Nhoule = nhoule
    HouleInit%depth = profondeur
    HouleInit%infinite_depth = infinite_depth
    allocate(HouleInit%w(nhoule), HouleInit%dir(nhoule) , HouleInit%Aphi(nhoule), HouleInit%konde(nhoule), HouleInit%lambda(nhoule))
    do j=1,nhoule
        HouleInit%w(j) = w(j)
        HouleInit%dir(j) = DirectionAiry(j)
        HouleInit%Aphi(j) = Aphi(j)
        HouleInit%konde(j) = konde(j)
        HouleInit%lambda(j) = lambda(j)
    end do
    
    ! Modification des paramètres de houle 
    Htype = HouleObj%Htype
    nhoule = HouleObj%Nhoule
    do j=1,nhoule
        w(j) = HouleObj%w(j)
        DirectionAiry(j) = HouleObj%dir(j)
        Aphi(j) = HouleObj%Aphi(j)
        konde(j) = HouleObj%konde(j)
        lambda(j) = HouleObj%lambda(j)
    end do
    
end subroutine CopyHoule

subroutine CopyBackHoule(HouleTemp)
    
    !f2py integer*1, dimension(1000):: HouleTemp
    type(THoule)                    :: HouleTemp    ! Wave structure.
    
    integer                         :: j            ! Loop parameter.
    
    ! This subroutine copies back a wave structure.
    
    Htype = HouleTemp%Htype
    Nhoule = HouleTemp%nhoule
    profondeur = HouleTemp%depth
    infinite_depth = HouleTemp%infinite_depth
    do j=1,nhoule
        w(j) = HouleTemp%w(j)
        DirectionAiry(j) = HouleTemp%dir(j)
        Aphi(j) = HouleTemp%Aphi(j)
        konde(j) = HouleTemp%konde(j)
        lambda(j) = HouleTemp%lambda(j)
    end do
    deallocate(HouleTemp%w, HouleTemp%dir , HouleTemp%Aphi, HouleTemp%konde, HouleTemp%lambda)

end subroutine CopyBackHoule

subroutine ps_riemann(mat,v1,v2,val)
    implicit none
    real(rp),dimension(2,2),intent(in) :: mat
    real(rp),dimension(2),intent(in) :: v1, v2
    real(rp), intent(out) :: val
!   local
    real(rp),dimension(2) :: vtemp
!
    vtemp = matmul(mat,v1)
    val = dot_product(v2,vtemp)    
!
end subroutine ps_riemann
    
subroutine Hull_function(x,z,y)
! Parameters
real(rp), intent(in) :: x, z
real(rp), intent(out) :: y
! Locals
integer :: type_Hull
real(rp) :: ztemp, alpha
! Begin
type_Hull = 1
ztemp = z
if (z.gt.Epsilon) ztemp = 0._RP ! On prolonge la carène au dessus de la surface libre par un plat bord.
if (type_Hull .eq. 1) then ! Modified Wigley Hull
    alpha = 0._RP
    y = (1-ztemp**2)*(1-x**2)*(1_RP + 0.2_RP*x**2) + alpha*ztemp**2*(1-ztemp**8)*(1-x**2)**4
elseif(type_Hull .eq. 2) then ! Sphere 
    if (abs(1._RP-x**2-ztemp**2).gt.Epsilon) then
        y = sqrt(1._RP-x**2-ztemp**2)
    else
        y = 0._RP
    end if
elseif(type_Hull .eq. 3) then ! Wigley Hull (these Delhommeau) 
    y = (1-ztemp**2)*(1-x**2)
end if
! End
return
end subroutine Hull_function

function DHullDx(x,z)
! Parameters
real(rp) :: DHullDx
real(rp), intent(in) :: x, z
! Locals
integer :: type_Hull
real(rp) :: y, ztemp, alpha
! Begin
type_Hull = 1
ztemp = z
if (z.gt.Epsilon) ztemp = 0._RP ! On prolonge la carène au dessus de la surface libre par un plat bord.
if (type_Hull .eq. 1) then ! Modified Wigley Hull
    alpha = 0._RP
    y = (1._RP-ztemp**2)*(-2._RP*x*(1._RP+0.2_RP*x**2)+0.4_RP*x*(1._RP-x**2)) - alpha*6._RP*x*ztemp**2*(1._RP-ztemp**8)*(1._RP-x**2)**3
elseif(type_Hull .eq. 2) then ! Sphere 
    if (abs(1._RP-x**2-ztemp**2).gt.Epsilon) then
        y = -x/sqrt(1._RP-x**2-ztemp**2)
    else
        y = 0._RP
    end if
elseif(type_Hull .eq. 3) then ! Wigley Hull (these Delhommeau) 
    y = -2_RP*x*(1-ztemp**2)
end if
DHullDx = y
! End
return
end function DHullDx

function DHullDz(x,z)
! Parameters
real(rp) :: DHullDz
real(rp), intent(in) :: x, z
! Locals
integer :: type_Hull
real(rp) :: y, ztemp, alpha
! Begin
type_Hull = 1
ztemp = z
if (z.gt.Epsilon) ztemp = 0._RP ! On prolonge la carène au dessus de la surface libre par un plat bord.
if (type_Hull .eq. 1) then ! Modified Wigley Hull
    alpha = 0._RP
    y = -2._RP*(1._RP-x**2)*((1._RP+0.2_RP*x**2) + alpha*4._RP*ztemp**8*(1._RP-x**2)**3 - alpha*(1._RP-ztemp**8)*(1._RP-x**2)**3)
elseif(type_Hull .eq. 2) then ! Sphere 
    if (abs(1._RP-x**2-ztemp**2).gt.Epsilon) then
        y = -ztemp/sqrt(1._RP-x**2-ztemp**2)
    else
        y = 0._RP
    end if
elseif(type_Hull .eq. 3) then ! Wigley Hull (these Delhommeau) 
    y = -2_RP*ztemp*(1-x**2)
end if
DHullDz = y
! End
return
end function DHullDz




end module FonctionsCommunes