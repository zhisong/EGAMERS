! subroutines to get the solution of a cubit function

module cubic

  ! FUNCTION
  ! cubicvalue  - real x   , y = d + c x + b x^2 + a x^3
  ! cubicvaluec - complex x,            "
  !
  ! SUBROUTINE
  ! cubicroot   - real inputs, find real roots for 0 == d + c x + b x^2 + x^3
  ! cubicrootc  - complex inputs, find complex roots for          "
  ! getcoeffs   - given f(x) and f'(x) at x1 and x2 to solve for the coeffs
  
  implicit none

contains
  
  real function cubicvalue(a, b, c, d, x)
    ! get the value of cubic equation y = d + c x + b x^2 + a x^3
    ! all coefficients are real
    implicit none

    real, intent(in) :: a, b, c, d, x
    
    cubicvalue = d + x * (c + x * (b  + x * a))

  end function cubicvalue

  complex function cubicvaluec(a, b, c, d, x)
    ! get the value of cubic equation y = d + c x + b x^2 + a x^3
    ! all coefficients are real, but x is complex
    implicit none

    real, intent(in) :: a, b, c, d
    complex, intent(in) :: x
    
    cubicvaluec = d + x * (c + x * (b  + x * a))
    
  end function cubicvaluec
    
  
  subroutine cubicroot(b, c, d, n, x)
    ! get solution for the cubic equation 0 == d + c x + b x^2 + x^3
    ! all coefficients and roots are real
    ! return n - number of roots, x - roots
    use paras_phy, only : pi
    implicit none
    
    real, intent(in) :: b, c, d
    integer, intent(out) :: n
    real, dimension(3), intent(out) :: x

    real :: q, r, theta, sqrtq, ca, cb, tmp
    integer :: i1, i2
    q = (b**2 - 3. * c) / 9.
    r = (2.* b**3 - 9.* c * b + 27.* d) / 54.
    if (r**2 .le. q**3) then
       ! there are three real roots
       n = 3
       sqrtq = sqrt(q)
       theta = acos(r / sqrtq**3)
       x(1) = -b/3. - 2. * sqrtq * cos(theta / 3.)
       x(2) = -b/3. - 2. * sqrtq * cos((theta + 2.*pi) / 3.)
       x(3) = -b/3. - 2. * sqrtq * cos((theta - 2.*pi) / 3.)
       ! sort ascendingly : bubble
       do i1 = 1, 2
          do i2 = 1, i1
             if (x(i1) .gt.x(i1+1)) then
                tmp = x(i1)
                x(i1) = x(i1+1)
                x(i1+1) = tmp
             end if
          end do
       end do
    else
       ! there is only one real root
       n = 1
       tmp = sqrt(r**2 - q**3)
       if (r .lt. 0.) tmp = - tmp
       ca = - (r + tmp)**(1./3.)
       if (ca .eq. 0.) then
          cb = 0.
       else
          cb = q / ca
       end if
       x(1) = (ca + cb) - b / 3.
    end if
  end subroutine cubicroot

  
  subroutine cubicrootc(b, c, d, n, x)
    ! get solution for the cubic equation 0 == c + b x + a x^2 + x^3
    ! return n - number of roots, x - roots
    ! a, b, c and x are complex
    use paras_phy, only : sqrt3o2i
    implicit none
    
    complex, intent(in) :: b, c, d
    integer, intent(out) :: n
    complex, dimension(3), intent(out) :: x
    
    complex :: q, r, ca, cb, tmp, b3, capcb, camcb
    integer :: i1, i2
    q = (b**2 - 3. * c) / 9.
    r = (2.* b**3 - 9.* b * c + 27.* d) / 54.

    tmp = sqrt(r**2 - q**3)
    
    if (real(conjg(r) * tmp) .le. 0) tmp = -tmp

    ca = - (r + tmp)**(1./3.)
    if (abs(ca) .eq. 0.) then
       cb = 0.
    else
       cb = q / ca
    end if

    n = 3
    capcb = ca + cb
    camcb = ca - cb
    b3 = b / 3.
    x(1) =  capcb - b3
    x(2) = -0.5 * capcb - b3 + sqrt3o2i * camcb
    x(3) = -0.5 * capcb - b3 - sqrt3o2i * camcb

  end subroutine cubicrootc

  subroutine getcoeffs(x2, y1, y1p, y2, y2p, c)
    ! f(x) = y = c(1) x^3 + c(2) x^2 + c(3) x + c(4)
    ! give y1 = f(0), y1p = f'(0), y2 = f(x2), y2p = f'(x2)
    ! solve for the coefficients

    implicit none
    
    real, intent(in) :: y1, y1p, x2, y2, y2p
    real, dimension(4), intent(out) :: c

    real :: xtmp, r1, r2

    if (x2 .eq. 0) then 
       write(*,*) 'getcoeffs: x2 cannot be zero'
    end if

    c(4) = y1
    c(3) = y1p

    r1 = ((y2 - c(4)) / x2 - c(3)) / x2
    r2 = (y2p - c(3)) / x2
    
    c(2) = (3.*r1 - r2)
    c(1) = (r1 - c(2)) / x2

  end subroutine getcoeffs

end module cubic
