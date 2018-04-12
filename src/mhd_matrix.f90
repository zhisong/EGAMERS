! module to fill in MHD part of the matrix

module mhd_matrix

  implicit none

  ! functions to calculate the value of integration
  ! IMPORTANT: core1 and core2 contains physics
  private :: intcore1, intcore2, getmhdint1, getmhdint2

  ! get the matrix for MHD part
  public :: getmat1, getmat2
  
contains

!  INTEGRATION CORE IN MATRIX 1 AND 2
!  (excluding the finite elements)
    real function intcore1(r)
    ! the integration core in matrix 1
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r

    intcore1 = r  * ni(r) / B0**2 * a**2

  end function intcore1
    
  real function intcore2(r)
    ! the integration core in matrix 2
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r

    intcore2 = r  * ni(r) / B0**2 * a**2 * omega_gam(r)**2 

  end function intcore2

!  ///////////////////////////////

  subroutine getmhdint1(ipos, nstep, m)
    ! calculate the integral S Hm(r) Hn(r) * intcore1(r) dr
    ! INPUT : the index of first and second element and number of steps
    ! OUTPUT: a 4*4 matrix containing the integral (2*2 Hermite elements)

    use hermite
    use radial_grid
    use matrix_module
    implicit none

    integer, intent(in) :: ipos, nstep
    complex, dimension(4,4), intent(out) :: m

    real :: x, xstart, xend, dxstep, inttmp
    integer :: i1, i2

    do i1 = 1, 4
       do i2 = 1, 4
          m(i1,i2) = 0
       end do
    end do

    if ((ipos .le. 0) .or. (ipos .ge. nelement)) return

    xstart = rgrid(ipos)
    xend = rgrid(ipos + 1)
    dxstep = (xend - xstart) / real(nstep)

    do i1 = 1, nstep
       x = xstart + real(i1 - 1) * dxstep + dxstep / 2.
       inttmp = intcore1(x)
       m(1,1) = m(1,1) + c1(x, xstart, -1., xend) * inttmp &
            * c1(x, xstart, -1., xend) * dxstep
       m(1,2) = m(1,2) + c1(x, xstart, -1., xend) * inttmp &
            * c2(x, xstart, -1., xend) * dxstep
       m(1,3) = m(1,3) + c1(x, xstart, -1., xend) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(1,4) = m(1,4) + c1(x, xstart, -1., xend) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(2,2) = m(2,2) + c2(x, xstart, -1., xend) * inttmp &
            * c2(x, xstart, -1., xend) * dxstep
       m(2,3) = m(2,3) + c2(x, xstart, -1., xend) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(2,4) = m(2,4) + c2(x, xstart, -1., xend) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(3,3) = m(3,3) + c1(x, xend, xstart, 2.) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(3,4) = m(3,4) + c1(x, xend, xstart, 2.) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(4,4) = m(4,4) + c2(x, xend, xstart, 2.) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
    end do

    do i1 = 1, 4
       do i2 = i1+1, 4
          m(i2,i1) = m(i1,i2)
       end do
    end do

  end subroutine getmhdint1

  subroutine getmhdint2(ipos, nstep, m)
    ! calculate the integral S Hm(r) Hn(r) * intcore2(r) dr
    ! INPUT : the index of first and second element and number of steps
    ! OUTPUT: a 4*4 matrix containing the integral (2*2 Hermite elements)

    use hermite
    use matrix_module
    use radial_grid
    implicit none

    integer, intent(in) :: ipos, nstep
    complex, dimension(4,4), intent(out) :: m

    real :: x, xstart, xend, dxstep, inttmp
    integer :: i1, i2

    do i1 = 1, 4
       do i2 = 1, 4
          m(i1,i2) = 0
       end do
    end do

    if ((ipos .le. 0) .or. (ipos .ge. nelement)) return

    xstart = rgrid(ipos)
    xend = rgrid(ipos + 1)
    dxstep = (xend - xstart) / real(nstep)

    do i1 = 1, nstep
       x = xstart + real(i1 - 1) * dxstep + dxstep / 2.
       inttmp = intcore2(x)
       m(1,1) = m(1,1) + c1(x, xstart, -1., xend) * inttmp &
            * c1(x, xstart, -1., xend) * dxstep
       m(1,2) = m(1,2) + c1(x, xstart, -1., xend) * inttmp &
            * c2(x, xstart, -1., xend) * dxstep
       m(1,3) = m(1,3) + c1(x, xstart, -1., xend) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(1,4) = m(1,4) + c1(x, xstart, -1., xend) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(2,2) = m(2,2) + c2(x, xstart, -1., xend) * inttmp &
            * c2(x, xstart, -1., xend) * dxstep
       m(2,3) = m(2,3) + c2(x, xstart, -1., xend) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(2,4) = m(2,4) + c2(x, xstart, -1., xend) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(3,3) = m(3,3) + c1(x, xend, xstart, 2.) * inttmp &
            * c1(x, xend, xstart, 2.) * dxstep
       m(3,4) = m(3,4) + c1(x, xend, xstart, 2.) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
       m(4,4) = m(4,4) + c2(x, xend, xstart, 2.) * inttmp &
            * c2(x, xend, xstart, 2.) * dxstep
    end do

    do i1 = 1, 4
       do i2 = i1+1, 4
          m(i2,i1) = m(i1,i2)
       end do
    end do

  end subroutine getmhdint2

  subroutine getmat1(mat1)
    ! get matrix1, dimension(2*nelement, 2*nelement)
    use paras_phy
    use paras_num, only : nmhdintegral
    use profile
    use radial_grid
    use matrix_module
    implicit none

    type(matrix), intent(out) :: mat1

    real :: dr1, rmid1, dr2, rmid2
    integer :: i1, i2, i3, ipos, nintegral
    complex, dimension(4,4) :: m
    
    nintegral = nmhdintegral

    if (allocated(mat1%data)) call matrix_destroy(mat1)
    ! matrix size = (n radial grid points) * (n Hermite elements) - 2
    ! - 2 is because boundary condition Er(0) = Er(1) is applied
    ! (c1 is absent in the first and last elements)
    call matrix_init(mat1, 2*nelement-2, 2*nelement-2) 

    ! clear all the data
    do i1 = 1, 2 * nelement-2
       do i2 = 1, 2 * nelement-2  
          mat1%data(i1,i2) = 0.
       end do
    end do

    ! for the first grid point, c1 is absent due to the boundary condition
    call getmhdint1(1, nintegral, m)
    do i2 = 1, 3
       do i3 = 1, 3
          mat1%data(i2, i3) = m(i2+1, i3+1)
       end do
    end do

    do i1 = 2, nelement-2
       call getmhdint1(i1, nintegral, m)
       ipos = i1 * 2 - 2
       do i2 = 1, 4
          do i3 = 1, 4
             mat1%data(ipos+i2-1, ipos+i3-1) = mat1%data(ipos+i2-1, ipos+i3-1) + m(i2, i3)
          end do
       end do
    end do
    
    ! for the last grid point, c1 is absent due to the boundary condition
    call getmhdint1(nelement-1, nintegral, m)
    do i2 = 1, 3
       m(3,i2) = m(4,i2)
       m(i2,3) = m(i2,4)
    end do
    ipos = (nelement - 1) * 2 - 2
    do i2 = 1, 3
       do i3 = 1, 3
          mat1%data(ipos+i2-1, ipos+i3-1) = mat1%data(ipos+i2-1, ipos+i3-1) + m(i2, i3)
       end do
    end do

  end subroutine getmat1

  subroutine getmat2(mat2)
    ! get matrix1, dimension(2*nelement, 2*nelement)
    use paras_phy
    use paras_num, only : nmhdintegral
    use profile
    use radial_grid
    use matrix_module
    implicit none

    type(matrix), intent(out) :: mat2

    real :: dr1, rmid1, dr2, rmid2
    integer :: i1, i2, i3, ipos, nintegral
    complex, dimension(4,4) :: m

    nintegral = nmhdintegral
    
    if (allocated(mat2%data)) call matrix_destroy(mat2)
    call matrix_init(mat2, 2*nelement-2, 2*nelement-2)
    
    ! for the first grid point, c1 is absent due to the boundary condition
    call getmhdint2(1, nintegral, m)
    do i2 = 1, 3
       do i3 = 1, 3
          mat2%data(i2, i3) = m(i2+1, i3+1)
       end do
    end do

    do i1 = 2, nelement-2
       call getmhdint2(i1, nintegral, m)
       ipos = i1 * 2 - 2
       do i2 = 1, 4
          do i3 = 1, 4
             mat2%data(ipos+i2-1, ipos+i3-1) = mat2%data(ipos+i2-1, ipos+i3-1) + m(i2, i3)
          end do
       end do
    end do
    
    ! for the last grid point, c1 is absent due to the boundary condition
    call getmhdint2(nelement-1, nintegral, m)
    do i2 = 1, 3
       m(3,i2) = m(4,i2)
       m(i2,3) = m(i2,4)
    end do
    ipos = (nelement - 1) * 2 - 2
    do i2 = 1, 3
       do i3 = 1, 3
          mat2%data(ipos+i2-1, ipos+i3-1) = mat2%data(ipos+i2-1, ipos+i3-1) + m(i2, i3)
       end do
    end do

  end subroutine getmat2



end module mhd_matrix

