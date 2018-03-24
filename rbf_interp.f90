module rbf_interp
use Constantes
use Parameters
use FonctionsCommunes

implicit none

contains 

subroutine phi_log(n,r,val)
    
    integer,intent(in)                  :: n    ! Size of r.
    real(rp),dimension(n),intent(in)    :: r    ! Norm.
    real(rp),dimension(n),intent(inout) :: val  ! Results.
    
    integer                             :: j    ! Loop parameter.
    
    ! This subroutine computes the radial basis function.
    
    do j = 1,n
        if(r(j).gt.Epsilon)then
            val(j) = r(j)**2 * log(r(j)) ! log(x) in Fortran 90 -> ln(x)
        else
            val(j) = 0._rp
        endif
    end do
    
end subroutine phi_log

subroutine rbf_funct(m,nd,xd,phi,w,ni,xi,fi)
    
    integer,intent(in)                      :: m        ! Number of coordinates in P : 2.
    integer,intent(in)                      :: nd       ! Number of controle points.
    real(rp),dimension(m,nd),intent(in)     :: xd       ! x and y positions of the controle points.
    external                                phi         ! Radial basis function.
    real(rp),dimension(:),intent(in)        :: w        ! Alpha and beta.
    integer,intent(in)                      :: ni       ! Number of mobile points.
    real(rp),dimension(m,ni),intent(in)     :: xi       ! x and y positions of the mobile points.
    real(rp),dimension(ni),intent(inout)    :: fi       ! Displacement of the mobile points with respect to x OR y.
    
    integer                                 :: j,k,kk   ! Loop parameters.
    real(rp),dimension(nd)                  :: r        ! Norms.
    real(rp),dimension(:),allocatable       :: val      ! Values of Phi_log.
    integer                                 :: ndim     ! Size of the linear system.
  
    ndim = nd + m + 1 ! m = 3
    allocate(val(ndim))

    do j = 1,ni ! Mobile points
        
        r = 0._RP
        
        do k = 1,nd ! Controle points
            do kk = 1,m
                r(k) = r(k) + (xi(kk,j)-xd(kk,k))**2
            end do
            r(k) = sqrt(r(k))
        end do
    
        call phi(nd,r,val(1:nd))
        val(nd+1) = 1
        val(nd+2:ndim) = xi(1:m,j)
    
        fi(j) = dot_product(val,w)
    
    end do

    deallocate(val)

end subroutine rbf_funct

subroutine rbf_computing_A(m,nd,xd,phi,mat,size_mat)
    
    integer,intent(in)                                  :: m        ! Number of coordinates in P : 2.
    integer,intent(in)                                  :: nd       ! Number of controle points.
    real(rp),dimension(m,nd),intent(in)                 :: xd       ! x and y positions of the controle points.
    external                                            phi         ! Radial basis function.
    integer,intent(in)                                  :: size_mat ! Size of A.
    real(rp),dimension(size_mat,size_mat),intent(inout) :: mat      ! A.
        
    integer                                             :: j,p,k,kk ! Loop parameters.
    real(rp),dimension(:),allocatable                   :: val      ! Values of Phi_log.
    real(rp),dimension(nd)                              :: r        ! Norms.
    integer                                             :: ndim     ! Size of the linear system.
        
    ! This subroutine computes the matrix A to then find alpha and beta in rbf_weight.
    
    ndim = nd+1+m ! m = 3
        
    allocate(val(nd))
    
    ! A.
    do j = 1,nd ! Controle points
        
        r = 0._RP
        
        ! Norm.
        do k = 1,nd
            do kk = 1,m
                r(k) = r(k) + (xd(kk,j)-xd(kk,k))**2
            end do
            r(k) = sqrt(r(k))
        end do   

        call phi(nd,r,val)

        mat(j,1:nd) = val(1:nd) ! M
        mat(j,nd+1) = 1 ! P(1)
        mat(j,nd+2:ndim) = xd(1:m,j) ! P(2:4)
    end do

    mat(nd+1,1:nd) = 1._rp ! Pt(1)
    
    do p = 1,nd
        do j = 1,m
            mat(nd+1+j,p) = xd(j,p) ! Pt(2:4)
        end do
    end do
    
    deallocate(val)
    
end subroutine rbf_computing_A

subroutine rbf_weight(m,nd,fd,w,mat,size_mat)
    
    integer,intent(in)                                  :: m        ! Number of coordinates in P : 2
    integer,intent(in)                                  :: nd       ! Number of controle points
    real(rp),dimension(nd),intent(in)                   :: fd       ! Displacement of the controle points with respect to x OR y
    real(rp),dimension(:),intent(inout)                 :: w        ! Alpha and beta
    integer,intent(in)                                  :: size_mat ! Size of A
    real(rp),dimension(size_mat,size_mat),intent(inout) :: mat      ! A
    
    integer                                             :: ndim     ! Loop parameters
    real(rp),dimension(:),allocatable                   :: rhs      ! X and B
    integer                                             :: saving_A ! Local saving of A in the subroutine LU
    integer                                             :: ierror   ! Error flag
    
    ! This subroutine computes the the coefficients alpha and beta of the radial basis function and the polynomial.
    
    ndim = nd+1+m ! m = 3
    
    allocate(rhs(ndim))
    
    rhs = 0._RP
    
    ! B
    rhs(1:nd) = fd(1:nd)
    rhs(nd+1:) = 0._rp
    
    ! AX = B
    saving_A = 1 ! The matrix mat will be reused, therefore it must be the same one in both input and output of the subroutine.
    call LU(mat,rhs,w,ndim,ierror,saving_A)
    
    deallocate(rhs)
    
end subroutine rbf_weight

end module rbf_interp
