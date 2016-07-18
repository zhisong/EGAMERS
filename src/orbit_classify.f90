! orbit classification module
module orbit_classify

  implicit none

  public :: tpbound, trapedge, traplost, ctpassinglost, copassingedge, costagedge
  private :: fbry, dfbrydr, dfbrydy, d2fbrydrdy

contains

  real function tpbound(mub0, pphi, itype)
    ! trapped/passing boundary finder
    ! itype: 1 (TYPE I), 2 (TYPE II), -1 (FAILED)

    use paras_phy
    use paras_num
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    integer, intent(out) :: itype

    real :: r1, r2, psi11, det, vpar
    real :: y1, y2
    real :: costheta, bb, cc, tmp
    integer :: flagnotfound, i, j
    real :: a11, a12, a21, a22, f1, f2
    real :: b11, b12, b21, b22, ddr, ddy

    if (pphi > 0) then
       ! no trapped particles for P_phi > 0
       tpbound = 0.
       return
    endif

    ! initial guess : v_\parallel = 0
    r1 = psitor(-pphi / ei)
    flagnotfound = 0
    i = 0

    ! treat the T/P boundary as TYPE I first (turning point on mid-plane at HFS)
    ! iterative method
    do while ((flagnotfound .eq. 0) .and. (i .le. maxittpbound))
       det = B0**2 * ei**2 * r1**2 * a**2 - 4. * mi * mub0 * q(r1)**2 * (1 - eps * r1)
       if  (det .lt. 0) then
          ! TYPE I search failed
          exit
       end if
       psi11 =  -pphi/ei - B0 / 2.  * r1 * a * (R0 - r1 * a)  / q(r1)
       psi11 = psi11 + sqrt(det) * (R0 - r1 * a) / 2. / q(r1) / ei
       if (psi11 .lt. 0.) then
          ! TYPE I search failed
          exit
       end if
       r2 = psitor(psi11)
       if ( abs((r2 - r1)) / r2 .le. errtpbound) then
          ! the t/p boundary is found
          flagnotfound = 1
          itype = 1
          y1 = -1.
       endif
       r1 = r2
       i = i + 1
    end do

    if (flagnotfound .eq. 0) then
       ! t/p boundary not found, try TYPE II
       do j = 1, 3
          r1 = psitor(-pphi / ei)
          if (j .eq. 1) then
             y1 = -0.5
          else if (j .eq. 2) then
             y1 = 0.
          else
             y1 = -1.
          endif

          flagnotfound = 0
          i = 0
          ! Newton's method to solve dtheta/dt=F(r,y)==0 and dF/dr==0
          do while ((flagnotfound .eq. 0) .and. (i .le. maxittpbound))
             ! matrix elements
             a11 = dfbrydr(r1, y1, mub0, pphi)
             a12 = dfbrydy(r1, y1, mub0, pphi)
             a21 = (dfbrydr(r1+ddx, y1, mub0, pphi) &
                  - dfbrydr(r1, y1, mub0, pphi)) / ddx
             a22 = d2fbrydrdy(r1, y1, mub0, pphi)
             ! inverse matrix
             det = a11*a22 - a12*a21
             if (det .eq. 0.) then
                ! matrix is ill behaved, return error
                write(*,*) 'ill'
                exit
             end if
             b11 =  a22 / det
             b12 = -a12 / det
             b21 = -a21 / det
             b22 =  a11 / det

             ! calculate the function value
             f1 = fbry(r1, y1, mub0, pphi)
             f2 = dfbrydr(r1, y1, mub0, pphi)
             ! calculate the shift
             ddr = f1 * b11 + f2 * b12
             ddy = f1 * b21 + f2 * b22

             r1 = r1 - ddr
             y1 = y1 - ddy

             if ((abs(ddr/r1) .le. errtpbound) .and. (abs(ddy/y1) .le. errtpbound)) then
                ! TYPE II boundary found
                if ((r1 .le. rlost) .and. (r1 .ge. 0.) .and. &
                     (y1 .le. 1.) .and. (y1 .ge. -1.)) then
                   itype = 2
                   flagnotfound = 1
                   psi11 = psi(r1)
                else
                   exit
                end if
             end if
             i = i + 1
          end do
          if (flagnotfound .gt. 0) exit
       end do
    end if

    if (flagnotfound .gt. 0) then
       ! t/p boundary found
       vpar = (pphi + ei * psi11) / mi / R0 / (1 + eps * r1 * y1)
       tpbound = 0.5 * mi * vpar**2 + mub0 * (1 - eps * r1 * y1)
    else
       tpbound = -1.
       itype = -1
    end if

  end function tpbound

  real function trapedge(mub0, pphi, istat)
    ! find the lowest energy limit for trapped particles
    ! istat : success(1), error(-1)
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    integer, intent(out) :: istat

    real :: r1
    
    if ((pphi .le. -ei * psi(1.)) .or. (pphi .ge. 0.)) then
       ! no trapped particles for P_phi > 0 and P_pphi < -e psi(1)
       trapedge = -1.
       istat = 0
       return
    end if

    trapedge = mub0 * (1. - eps * psitor(-pphi/ei))
    istat = 1
  end function trapedge

  real function costagedge(mub0, pphi, istat)
    ! find the lowest energy limit for co-passing stagnation paticles
    ! istat : success(1), lost(2), none(0), error(-1)
    use paras_phy
    use paras_num
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    integer, intent(out) :: istat

    real :: r1, r2, psi11, det, vpar, f1, ddr
    integer :: flagnotfound, i

    if ((pphi .le. -ei * psi(1.))) then
       ! no  particles  P_pphi < -e psi(1)
       costagedge = -1.
       istat = 0
       return
    end if

    ! initial guess : v_\parallel = 0
    r1 = psitor(-pphi / ei)
    r1 = 1.
    flagnotfound = 0
    i = 0

    ! Newton's methods
    do while ((flagnotfound .eq. 0) .and. (i .le. maxittpbound))
       f1 = fbry(r1, 1., mub0, pphi)
       det = dfbrydr(r1, 1., mub0, pphi)
       if (det .eq. 0.) then
          ! matrix singluar, search failed
          exit
       end if
       ddr = f1 / det
       r2 = r1 - ddr
       if (abs(ddr/r2) .le. errtpbound) then
          ! the boundary is found
          flagnotfound = 1
       end if
       r1 = r2
       i = i + 1
    end do

    if (flagnotfound .eq. 0) then
       ! boundary not found
       costagedge = -1.
       istat = -1
       return
    end if
    
    if (r2 .gt. 1) then
       istat = 2
    else
       istat = 1
    end if

    psi11 = psi(r2)
    vpar = (pphi + ei * psi11) / mi / R0 / (1 + eps * r2)
    costagedge = 0.5 * mi * vpar**2 + mub0 * (1 - eps * r2) 
    ! slightly raise the energy to prevent ill-behaviour in numerics
    costagedge = costagedge * (1 + 1.e-12)
  end function costagedge

  real function traplost(mub0, pphi)
    ! find the trapped-confined/trapped-lost boundary

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    real :: eeout

    eeout = mub0 * (1 - eps) + (pphi + ei * psi(1.))**2/mi/2./R0**2/(1.+eps)**2
    traplost = eeout
  end function traplost

  real function ctpassinglost(mub0, pphi)
    ! find the ct-passing-confined/lost boundary

    use paras_phy
    use profile
    implicit none

    real :: mub0, pphi, ee
    real :: eeout

    eeout = mub0 * (1 + eps) + (pphi + ei * psi(1.))**2/mi/2./R0**2/(1.-eps)**2
    ctpassinglost = eeout

  end function ctpassinglost

  real function copassingedge(mub0, pphi)
    ! find the co-passing edge

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    real :: eeout

    eeout = mub0 + (pphi)**2/mi/2./R0**2
    copassingedge = eeout

  end function copassingedge

  real function fbry(r, y, mub0, pphi)
    ! function used to find TYPE II t/p boundary
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r, y, mub0, pphi
    real :: Rbar

    Rbar = (1. + eps * r * y)
    fbry = r * a * ei * B0 * R0 * Rbar * (pphi + ei * psi(r))
    fbry = fbry - q(r) * ((pphi + ei * psi(r))**2  + mi * mub0 * R0**2 * Rbar**3) * y

  end function fbry

  real function dfbrydr(r, y, mub0, pphi)
    ! function used to find TYPE II t/p boundary
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r, y, mub0, pphi
    real :: Rbar

    Rbar = (1. + eps * r * y)
    dfbrydr = B0 * ei * R0 * (pphi + ei * psi(r)) 
    dfbrydr = dfbrydr + B0**2 * ei**2 * R0 * r**2 * a**2 * Rbar / q(r)
    dfbrydr = dfbrydr - mi * y * 3.* mub0 * R0 * y * Rbar**2 * q(r)
    dfbrydr = dfbrydr -  y * (mi * mub0 * R0**2 * Rbar**3 + (pphi + ei*psi(r))**2) * dqdr(r) / a

  end function dfbrydr

  real function dfbrydy(r, y, mub0, pphi)
    ! function used to find TYPE II t/p boundary
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r, y, mub0, pphi
    real :: Rbar

    Rbar = (1. + eps * r * y)
    dfbrydy = -q(r) * ((pphi + ei * psi(r))**2  + mi * mub0 * R0**2 * Rbar**3)
    dfbrydy = dfbrydy  -q(r) * 3. * (mi * mub0 * R0**2 * Rbar**2) * y * eps * r
    dfbrydy = dfbrydy + r * a * ei * B0 * R0 * eps * r * (pphi + ei * psi(r))
  end function dfbrydy

  real function d2fbrydrdy(r, y, mub0, pphi)
    ! function used to find TYPE II t/p boundary
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r, y, mub0, pphi
    real :: Rbar

    Rbar = (1. + eps * r * y)
    d2fbrydrdy = B0**2 * ei**2 * r**3 * a**3 / q(r)
    d2fbrydrdy = d2fbrydrdy - 6. * mi * mub0 * R0 * y * Rbar**2 * q(r)
    d2fbrydrdy = d2fbrydrdy - 6. * mi * mub0 * R0 * y**2 * Rbar * q(r) * r * eps
    d2fbrydrdy = d2fbrydrdy - mi * (mub0 * R0**2 * Rbar**3 + (pphi + ei*psi(r))**2 / mi) * dqdr(r) / a
    d2fbrydrdy = d2fbrydrdy - 3. * mi * y * (mub0 * R0**2 * Rbar**2 * r * eps) * dqdr(r) / a

  end function d2fbrydrdy

end module orbit_classify
