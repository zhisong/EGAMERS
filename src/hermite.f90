! finite elements used in discretizing the radial mode

module hermite

  implicit none
  private
  public :: c1, c2
  
contains

  real function c1(x, x0, x1, x2)
    ! cubic Hermite element base function 1
    ! range: x1 to x2, mid node x0
    implicit none
    real, intent(in) :: x, x0, x1, x2

    if ((x.ge.x1) .and. (x.le.x0)) then
       c1 = ((x-x1)/(x0-x1))**2 * (3. - 2*(x-x1)/(x0-x1))
    else if ((x.gt.x0) .and. (x.le.x2)) then
       c1 = ((x-x2)/(x2-x0))**2 * (3 - 2*(x2-x)/(x2-x0))
    else
       c1 = 0
    end if
  end function c1

  real function c2(x, x0, x1, x2)
    ! cubic Hermite element base function 2
    ! range: x1 to x2, mid node x0
    implicit none
    real, intent(in) :: x, x0, x1, x2

    if ((x.ge.x1) .and. (x.le.x0)) then
       c2 = ((x-x1)/(x0-x1))**2 * (x-x0)
    else if ((x.gt.x0) .and. (x.le.x2)) then
       c2 = ((x-x2)/(x2-x0))**2 * (x-x0)
    else
       c2 = 0
    end if
  end function c2

end module hermite
