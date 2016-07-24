! deal with Landau type of integral

module landau_integral

  use paras_phy, only : pi
  implicit none

contains

  complex function landauint(f, g, t1, b, a)
    ! evaluate the value of Landau type integral
    !   a  f(x)
    !  S  ------ dx  = ?
    !   b g(x)-t
    ! t is a complex defined on the lower half of the complex plane.
    ! f(x) and g(x) are both third order polynomial functions.
    ! INPUT: f is the coefficients of a third order polynomial
    !            y = f(1) x^3 + f(2) x^2 + f(3) x + f(4)
    !        g is the coefficients of a third order polynomial 
    !             y = g(1) x^3 + g(2) x^2 + g(3) x + g(4)
    !        t is the ploe, a and b are upper and lower limits of the integral
    
    use cubic
    implicit none

    real, dimension(4), intent(in) :: f, g
    complex, intent(in) :: t1
    real, intent(in) :: a, b

    complex, dimension(3) :: gx
    integer :: n, nrange
    complex :: ans, t
    real :: angle, dev, xreal
    integer, dimension(3) :: sgn
    
    t = t1
    ! if t gives damping, move it to the real axis
    ! since the analytic continuation does not make sense for splines
    if (imag(t) .gt. 0.) t = real(t)

    call cubicrootc(cmplx(g(2)/g(1)), cmplx(g(3)/g(1)),&
         (g(4) - t)/g(1), n, gx)
    landauint = f(1) / g(1) * (a - b)

    if (imag(t) .eq. 0) then
       ! the pole on real axis, need to treat the Landau integral
       sgn(1) = 1
       sgn(2) = 1
       sgn(3) = 1
       nrange = 0
       xreal = real(gx(1))
       if ((real(b) .le. xreal) .and. (real(a) .ge. xreal)) then
          dev = (3.* g(1) * xreal + 2.* g(2)) * xreal + g(3)
          if (dev .lt. 0.) then
             sgn(1) = -1
          end if
          nrange = nrange + 1
       end if

       xreal = real(gx(2))
       if ((real(b) .le. xreal) .and. (real(a) .ge. xreal)) then
          dev = (3.* g(1) * xreal + 2.* g(2)) * xreal + g(3)
          if (dev .lt. 0.) then
             sgn(2) = -1
          end if
          nrange = nrange + 1
       end if

       xreal = real(gx(3))
       if ((real(b) .le. xreal) .and. (real(a) .ge. xreal)) then
          dev = (3.* g(1) * xreal + 2.* g(2)) * xreal + g(3)
          if (dev .lt. 0.) then
             sgn(3) = -1
          end if
          nrange = nrange + 1
       end if

       if (nrange .gt. 0) then
          if (nrange .eq. 2) then
             if (sgn(1) * sgn(2) * sgn(3) .gt. 0) then
                ! fake roots
                sgn(1) = 0
                sgn(2) = 0
                sgn(3) = 0
             end if
          else if (nrange .eq. 3) then
             if (sgn(1) * sgn(2) .gt. 0) then
                ! fake roots
                sgn(1) = 0
                sgn(2) = 0
             else if (sgn(2) * sgn(3) .gt. 0) then
                sgn(2) = 0
                sgn(3) = 0
             end if
          end if
       end if
    else
       ! the pole on the lower plane, doesn't matter
       sgn(1) = 0
       sgn(2) = 0
       sgn(3) = 0
    end if

    ans = cubicvaluec(f(1), f(2), f(3), f(4), gx(1)) * getlog(gx(1), a, b, sgn(1))
    ans =  ans / (gx(1) - gx(2)) / (gx(1) - gx(3)) / g(1)
    landauint = landauint + ans
    
    ans = cubicvaluec(f(1), f(2), f(3), f(4), gx(2)) * getlog(gx(2), a, b, sgn(2))
    ans = ans / (gx(2) - gx(1)) / (gx(2) - gx(3)) / g(1)
    landauint = landauint + ans

    ans = cubicvaluec(f(1), f(2), f(3), f(4), gx(3)) * getlog(gx(3), a, b, sgn(3))
    ans = ans / (gx(3) - gx(1)) / (gx(3) - gx(2)) / g(1)
    landauint = landauint + ans

  end function landauint

  complex function getlog(t, a, b, sgn)
    ! use to determine value of the log function in a Landau integral
    use paras_phy, only : pi
    implicit none

    complex, intent(in) :: t
    real, intent(in) :: a, b
    integer, intent(in) :: sgn

    real :: angle
    
    if (imag(t) .eq. 0) then 
       if ((b .lt. real(t)) .and. (a .ge. real(t))) then
          if (sgn .gt. 0) then
             angle = -pi
          else if (sgn .lt. 0) then
             angle = pi
          else
             angle = 0.
          end if
       else
          angle = 0.
       end if
    else if (imag(t) .ge. 0) then
       if ((b .lt. real(t)) .and. (a .ge. real(t)) .and. (sgn .gt. 0)) then
          angle = atan2(imag((a-t)/(b-t)), real((a-t)/(b-t))) - 2 * pi
       else
          angle = atan2(imag((a-t)/(b-t)), real((a-t)/(b-t)))
       end if
    else
       if ((b .lt. real(t)) .and. (a .ge. real(t)) .and. (sgn .lt. 0)) then
          angle = atan2(imag((a-t)/(b-t)), real((a-t)/(b-t))) + 2 * pi
       else
          angle = atan2(imag((a-t)/(b-t)), real((a-t)/(b-t)))
       end if
    end if
    getlog = log(abs(a-t)/abs(b-t)) + angle * (0.,1.)

  end function getlog

end module landau_integral
    

