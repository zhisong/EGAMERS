! Using Newton's method to solve the eigen value

module eigen

  ! use Newton's method (inversed vector iteration) to iteratively
  ! find the eigenvalue and eigenvector around a given value

  use matrix_module
  use mhd_matrix
  use trap_matrix
  use radial_grid, only : nelement
  use interfaces
  implicit none

  public

  type(tmatrix), public :: tm
  type(matrix), public :: mat1, mat2
  complex, public :: lambda
  complex, allocatable, dimension(:), public :: v
contains

  subroutine newton_init(lambdain)
    ! initialize the newton iteration
    ! PLEASE DON'T forget to initialize the trap matrix seperately
    implicit none
    
    complex, intent(in) :: lambdain
    integer :: i1

    lambda = lambdain
    write(*,*) 'START ITERATION: OMEGA', real(sqrt(lambda)), 'GAMMA', imag(sqrt(lambda))
    write(*,*)
    call getmat1(mat1)
    call getmat2(mat2)

    allocate(v(2*nelement - 2))
    do i1 = 1, 2*nelement - 2
       v(i1) = 1. / real(2*nelement - 2)
    end do
  end subroutine newton_init
    
  subroutine newton_cleanup()
    ! clean up the allocated objects

    call matrix_destroy(mat1)
    call matrix_destroy(mat2)
    if (allocated(v)) deallocate(v)

  end subroutine newton_cleanup

  subroutine newton_step()
    ! newton iteration
    use paras_num, only : erreig, nmaxit
    implicit none

    type(matrix) :: mat, matprime, mat1, matprime1
    complex :: lambda1, res1, res2
    real :: dlambda
    complex, allocatable, dimension(:) :: vwork, tpdotv, tdotv
    integer, allocatable, dimension(:) :: pmat
    integer :: i1, i2, i3, nele, info
    
    nele = 2*nelement - 2

    ! create and copy to temporary work space
    allocate(vwork(nele))
    allocate(tpdotv(nele))
    allocate(tdotv(nele))
    allocate(pmat(nele))

    do i2 = 1, nmaxit

       if (i2 .eq. 2) then
          mat = mat1
          matprime = matprime1
       else if (i2 .eq. 1) then
          call getmat(lambda, mat, matprime)
          mat1 = mat
          matprime1 = matprime
       else
          call getmat(lambda, mat, matprime)
       end if
       ! solve the linear system T(n) * u(n+1) = T'(n) * v(n)
       ! by PLU decomposition of mat

       ! Step 0, store T(n) * v(n)
       call czgemv('N', nele, nele, (1.,0.), mat%data, nele, v, 1, &
            (0.,0.), tdotv, 1)

       ! Step 1, compute tpdotv = T'(n) * v(n)
       call czgemv('N', nele, nele, (1.,0.), matprime%data, nele, v, 1, &
            (0.,0.), tpdotv, 1)
       ! copy to vwork
       call czcopy(nele, tpdotv, 1, vwork, 1)

       ! Step 2, solve T(n) * u(n+1) = vwork
       call czgetrf(nele, nele, mat%data, nele, pmat, info)
       if (info .ne. 0) then
          write(*,*) 'ERROR IN CZGETRF, INFO =', info
       end if
       call czgetrs('N', nele, 1, mat%data, nele, pmat, vwork, nele, info)
       if (info .ne. 0) then
          write(*,*) 'ERROR IN CZGETRS, INFO =', info
       end if

       ! Step 3, calculate xt(n) T(n) x(n) and xt(n) T'(n) x(n)
       ! lambda(n+1) = lambda -  xt(n) T(n) x(n) / xt(n) T'(n) x(n)
       res1 = czdotu(nele, v, 1, tpdotv, 1) ! xt(n) T'(n) x(n)
       res2 = czdotu(nele, v, 1, tdotv , 1) ! xt(n) T(n) x(n)
       lambda1 = lambda - res2 / res1

       ! Step 4, v(n+1) = u(n+1) / ||u(n+1)||
       res1 = czdotc(nele, vwork, 1, vwork, 1) ! u(n+1)**H * u(n+1)
       call czscal(nele, 1./sqrt(res1), vwork, 1) ! v(n+1) = u(n+1) / ||u(n+1)||

       dlambda = abs(lambda1 - lambda) / abs(lambda1)

       ! update the eigenvalue and eigenvector
       if (i2 .gt. 1) then
          ! discard the first two update of lambda
          lambda = lambda1

          write(*,*) 'ITERATION', i2-1, ' DLAMBDA =', dlambda
          write(*,*) 'OMEGA', real(sqrt(lambda)), 'GAMMA', imag(sqrt(lambda))
          write(*,*)
       end if
       call czcopy(nele, vwork, 1, v, 1) 
       

       ! if precision is reached
       if (dlambda .le. erreig) then
          exit
       end if

    end do
    call matrix_destroy(mat)
    call matrix_destroy(matprime)
    call matrix_destroy(mat1)
    call matrix_destroy(matprime1)
    deallocate(vwork)
    deallocate(tdotv)
    deallocate(tpdotv)
    deallocate(pmat)

  end subroutine newton_step

  subroutine getmat(lambdain, mat, matprime)
    ! get T and T'
    ! INPUT:  lambda
    ! OUTPUT: T(lambda), T'(lambda)
    use paras_num, only : dlambdapctg
    implicit none

    complex, intent(in) :: lambdain
    type(matrix), intent(out) :: mat, matprime

    type(matrix) :: mat3trap, mat3trap1
    complex :: omega, omega1, dlambda

    dlambda = real(dlambdapctg * lambdain)
    omega  = sqrt(lambdain)
    omega1 = sqrt(lambda + dlambda)

    call getmat3trap(tm, omega , mat3trap )
    call getmat3trap(tm, omega1, mat3trap1)

    mat = mat2 + mat3trap + (-lambdain) * mat1
    matprime =  (1./dlambda) * mat3trap1 + (-1./dlambda) * mat3trap &
         + (-1.,0.) * mat1

!!$    mat = mat2 + (-lambdain) * mat1
!!$    matprime =  (-1.,0.) * mat1
    
    call matrix_destroy(mat3trap)
    call matrix_destroy(mat3trap1)

  end subroutine getmat

end module eigen
