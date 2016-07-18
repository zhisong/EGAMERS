! calculating the orbit integral of 1/pi* H(r)*dr/d zeta * sin(p zeta) 

module orbit_integral

  implicit none
  private
  public :: orbit_int

contains

  subroutine orbit_int(r, theta, norbit, np, vpm, imin, imax)
    ! Inputs : orbitdata(norbit, r theta), number of finite element, number of p
    ! Output : the value of the integral (dimension 2*nelement-2, np),
    !          the index of the most inside/outside elements

    use paras_num
    use paras_phy
    ! sintable: pre-calculated sin function to speed up the calculation
    use sintable
    use hermite
    use radial_grid

    implicit none
    integer, intent(in) :: norbit, np
    real, dimension(norbit), intent(in) :: r, theta
    real, dimension(2 * nelement - 2, np), intent(out) :: vpm
    integer, intent(out) :: imin, imax

    real, dimension(2* nmaxelement, npmax) :: work
    real :: dzeta, dr, drdzeta
    real :: maxr, minr, x0, x1, x2
    real :: tmpint
    integer :: maxrpos, minrpos
    integer :: i1, i2, i3

    ! finite elements dx
    dr = 1. / real(nelement+1)
    dzeta = 2 * pi / real(norbit)

    ! find max and min of r
    maxr = r(1)
    minr = r(floor(norbit / 2.) + 1)
    if (maxr .lt. minr) then
       tmpint = maxr
       maxr = minr
       minr = tmpint
    end if
    
    do i1 = 2, nelement
       if (minr .lt. rgrid(i1)) exit
    end do
    minrpos = i1 - 1

    do i1 = nelement-1, 1, -1
       if (maxr .gt. rgrid(i1)) exit
    end do
    maxrpos = i1 + 1

    if (maxrpos .ge. nelement + 1) then
       maxrpos = nelement
    end if
    if (minrpos .le. 0) then
       minrpos = 1
    end if

    ! clear data
    do i1 = 1, nelement
       do i2 = 1, np
          work(i1*2-1, i2) = 0.
          work(i1*2  , i2) = 0.
       end do
    end do
    ! orbit integral
    do i1 = 1, norbit-1
       drdzeta = (r(i1+1) - r(i1))
       do i2 = 1, np
          tmpint = drdzeta * sinpzeta(i1, i2)
          do i3 = minrpos, maxrpos
             x0 = rgrid(i3)
             if (i3 .eq. 1) then
                x1 = -0.1
             else
                x1 = rgrid(i3 - 1)
             end if
             if (i3 .eq. nelement) then
                x2 = 1.1
             else
                x2 = rgrid(i3+1)
             end if
             if ((i3 .ne. 1) .and. (i3 .ne. nelement)) then
                ! ignore the c1 element for first and last grid points
                work(i3*2-1, i2) = work(i3*2-1,i2) &
                     + c1((r(i1)+r(i1+1))/2.,x0,x1,x2) * tmpint 
             end if
             work(i3*2  , i2) = work(i3*2  ,i2) &
                  + c2((r(i1)+r(i1+1))/2.,x0,x1,x2) * tmpint 
          end do
       end do
    end do
    ! last data point
    i1 = norbit
    drdzeta = (r(1) - r(i1))
    do i2 = 1, np
       tmpint = drdzeta * sinpzeta(i1, i2)
       do i3 = minrpos, maxrpos
          x0 = rgrid(i3)
             if (i3 .eq. 1) then
                x1 = -0.1
             else
                x1 = rgrid(i3 - 1)
             end if
             if (i3 .eq. nelement) then
                x2 = 1.1
             else
                x2 = rgrid(i3+1)
             end if
             if ((i3 .ne. 1) .and. (i3 .ne. nelement)) then
                work(i3*2-1, i2) = work(i3*2-1,i2) + c1((r(i1)+r(1))/2.,x0,x1,x2) * tmpint
             end if
             work(i3*2,i2) = work(i3*2,i2) + c2((r(i1)+r(1))/2.,x0,x1,x2) * tmpint
       end do
    end do
    
    ! transfering to vpm and adding the 1/pi factor
    do i2 = 1, np
       ! special treatment for first and last element
       ! drop the c1 element
       vpm(1, i2) = work(2, i2) / pi
       vpm(2*nelement-2, i2) = work(2*nelement, i2) / pi
       do i3 = 2, nelement-1
          vpm(i3*2-2,i2) = work(i3*2-1,i2) / pi
          vpm(i3*2-1,i2) = work(i3*2  ,i2) / pi
       end do
    end do

    ! the most inside element index
    if (minrpos .eq. 1) then
       imin = 1
    else
       imin = 2 * minrpos - 2
    end if
    
    ! the most outside element index
    if (maxrpos .eq. 1) then
       imax = 1
    else if (maxrpos .eq. nelement) then
       imax = 2 * nelement - 2
    else
       imax = 2 * maxrpos - 1
    end if
    
  end subroutine orbit_int

end module orbit_integral
