! orbit classification module
module orbit_classify

  implicit none

  ! functions
  ! (public)
  ! tpbound      - the t/p boundary energy, given mub0 and pphi
  ! stagedge     - the energy edge for stagnation particles, given mub0 and pphi
  ! stagedgepphi - the pphi edge for stagnation particles, given mub0 and energy
  ! traplost     - lost energy for trapped particles, given mub0 and pphi
  ! coplost      - lost energy for co-passing particles, given mub0 and pphi
  ! coplostpphi  - lost pphi for co-passing, given mub0 and energy
  ! ctplost      - lost energy for ct-passing particles, given mub0 and pphi
  ! ctplostpphi  - lost pphi for ct-passing particles, given mub0 and pphi
  ! eeonaxis     - energy of particles through the axis, given mub0 and pphi
  ! losttpbound  - the pphi where cop lost boundary hits the tpbound
!!$  ! type12tpbound- find the pphi where TYPE I t/p boundary turns into TYPE II
  !
  ! (private)
  ! rvtdr, drvtdr- functions used in stagedgepphi
  ! fbry, dfbrydr, dfbrydy, df2brydrdy - functions used in stagedge and tpbound
  
  public :: tpbound, stagedge, traplost, coplost, ctplost, eeonaxis, &
       stagedgepphi, coplostpphi, ctplostpphi, losttpbound
  
  private

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

    itype = 0

    if (pphi .gt. 0.) then
       ! no trapped particles for P_phi > 0
       tpbound = 0.
       return
    else if (pphi .eq. 0.) then
       ! special treatment
       tpbound = mub0
       itype = 2
       return
    end if

    ! initial guess : v_\parallel = 0
    r1 = psitor(-pphi / ei)
    flagnotfound = 0
    i = 0
    itype = 0
   
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
       do j = 0, 3
          r1 = psitor(-pphi / ei)
          r1 = r1 / 2.
          if (j .eq. 0) then
             y1 = -0.75
          else if (j .eq. 1) then
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

             r1 = min(max(r1 - ddr, 0.0), 1.0)
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

!!$  real function type12tpbound(mub0, istat)
!!$    ! find the pphi value where type I tpbound turns into type II
!!$    use paras_phy, only : ei
!!$    use paras_num, only : errtpbound, maxittpbound
!!$    use profile, only : psi1
!!$    implicit none
!!$
!!$    real, intent(in) :: mub0
!!$    integer, intent(out) :: istat
!!$
!!$    integer :: i1, typemid, ifound
!!$    real :: pphileft, pphiright, pphimid, tpboundtmp
!!$
!!$    pphileft = - ei * psi1
!!$    pphiright = 0.
!!$    
!!$    i1 = 1
!!$    ifound = 0
!!$    do while ((ifound .le. 0) .and. (i1 .le. maxittpbound))
!!$       ! binary search
!!$       pphimid = (pphileft + pphiright) / 2.
!!$       tpboundtmp = tpbound(mub0, pphimid, typemid)
!!$       write(*,*) pphimid / ei / psi1, typemid
!!$       if (typemid .le. 0) then
!!$          write(*,*) 'tpbound return error in type12tpbound'
!!$          i1 = maxittpbound
!!$          exit
!!$       end if
!!$
!!$       if (abs((pphimid - pphileft)/pphimid) .le. errtpbound) then
!!$          ! precision reached
!!$          ifound = 1
!!$       end if
!!$       
!!$       if (typemid .eq. 1) then
!!$          pphileft = pphimid
!!$       else
!!$          pphiright = pphimid
!!$       end if
!!$          
!!$       i1 = i1 + 1
!!$       
!!$    end do
!!$
!!$    if (ifound .le. 0) then
!!$       write(*,*) 'type12tpbound failed'
!!$       istat = -1
!!$       type12tpbound = -100
!!$    else
!!$       type12tpbound = pphimid
!!$    end if
!!$    
!!$  end function type12tpbound
  
  real function stagedge(mub0, pphi, vsign, istat)
    ! find the edge for stagnation paticles (fixed pphi, output energy)
    ! vsign : =1 co-passing/trapped, =-1 ct-passing 
    ! istat : success(1), none(0), error(-1)
    use paras_phy
    use paras_num
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    integer, intent(in) :: vsign
    integer, intent(out) :: istat

    real :: r1, r2, psi11, det, vpar, f1, ddr, side
    integer :: flagnotfound, i

    if ((pphi .le. -ei * psi(1.)) .and. (vsign .eq. 1)) then
       ! no  particles  P_pphi < -e psi(1) if vsign == 1
       stagedge = -1.
       istat = 0
       return
    else if ((pphi .ge. 0.) .and. (vsign .eq. -1)) then
       ! no  particles  P_pphi > 0 if vsign == -1
       stagedge = -1.
       istat = 0
       return
    end if
    if (vsign .eq. 1) then
       side = 1.
    else if (vsign .eq. -1) then
       side = -1.
    else
       ! vsign must be 1 or -1
       write(*,*) 'In stagedge : vsign must be 1 or -1'
       stagedge = -1.
       istat = 0
       return
    end if

    ! initial guess : v_\parallel = 0
    !    r1 = psitor(-pphi / ei)
    if (vsign .eq. -1) then
       r1 = 0.
    else
       r1 = 1.
    end if
    flagnotfound = 0
    i = 0
    ! Newton's methods
    do while ((flagnotfound .eq. 0) .and. (i .le. maxittpbound))
       f1 = fbry(r1, side, mub0, pphi)
       det = dfbrydr(r1, side, mub0, pphi)
       if (det .eq. 0.) then
          ! matrix singluar, search failed
          exit
       end if
       ddr = f1 / det
       r2 = r1 - ddr

       if ((r2 .gt. 1.) .or. (r2 .le. 0.)) then
          ! search failed
          stagedge = -1
          istat = -1
          return
       end if
    
       if (abs(ddr/r2) .le. errtpbound) then
          ! the boundary is found
          flagnotfound = 1
       end if
       r1 = r2
       i = i + 1
    end do
    if (flagnotfound .eq. 0) then
       ! boundary not found
       stagedge = -1.
       istat = -1
       return
    end if

    psi11 = psi(r2)
    vpar = (pphi + ei * psi11) / mi / R0 / (1 + side * eps * r2)
    stagedge = 0.5 * mi * vpar**2 + mub0 * (1 - side * eps * r2)
    istat = 1
    
    ! slightly modify the energy to prevent ill-behaviour in numerics
    stagedge = stagedge * (1 + side * 1e-10)
  end function stagedge

  real function stagedgepphi(mub0, ee, vsign, istat)
    ! find the edge for stagnation paticles (fixed energy, output pphi)
    ! vsign : =1 co-passing/trapped, =-1 ct-passing 
    ! istat : success(1), none(0), error(-1)
    use paras_phy
    use paras_num
    use profile
    implicit none

    real, intent(in) :: mub0, ee
    integer, intent(in) :: vsign
    integer, intent(out) :: istat

    real :: r1, r2, det, vpar, f1, ddr, side, mub, Rbar
    integer :: flagnotfound, i

    if ((ee .le. mub0 * (1 - eps)) .and. (vsign .eq. 1)) then
       ! no  particles E < mub0 (1-a/R0) if vsign == 1
       stagedgepphi = -1.
       istat = 0
       return
    else if ((ee .le. mub0) .and. (vsign .eq. -1)) then
       ! no  particles E < mub0 if vsign == -1
       stagedgepphi = -1.
       istat = 0
       return
    end if
    if ((vsign .ne. 1) .and. (vsign .ne. -1)) then
       ! vsign must be 1 or -1
       write(*,*) 'In stagedgephi : vsign must be 1 or -1'
       stagedgepphi = -1.
       istat = 0
       return
    end if

    ! initial guess : r = 0 or 1
    if (vsign .eq. 1) then
       r1 = 1.
    else
       r1 = 0.
    end if
    flagnotfound = 0
    i = 0
    side = float(vsign)
    ! Newton's methods
    do while ((flagnotfound .eq. 0) .and. (i .le. maxittpbound))
       f1 = rvt(r1, mub0, ee, vsign)
       det = drvtdr(r1, mub0, ee, vsign)
       if (det .eq. 0.) then
          ! matrix singluar, search failed
          write(*,*) 'singular matrix'
          exit
       end if
       ddr = f1 / det
       r2 = r1 - ddr
       
       if ((r2 .gt. 1.) .or. (r2 .le. 0.)) then
          ! search failed
          stagedgepphi = -1
          istat = -1
          return
       end if
    
       if (abs(ddr/r2) .le. errtpbound) then
          ! the boundary is found
          flagnotfound = 1
       end if
       r1 = r2
       i = i + 1
    end do
    if (flagnotfound .eq. 0) then
       ! boundary not found
       stagedgepphi = -1.
       istat = -1
       return
    end if

    mub = mub0 * (1 - r2 * eps * side)
    Rbar = 1 + r2 * eps * side
    vpar = sqrt((2*ee - 2*mub) / mi)
    stagedgepphi = - ei * psi(r2) + mi * vpar * Rbar * R0 * side
    istat = 1

    ! slightly modify the energy to prevent ill-behaviour in numerics
    if (stagedgepphi .le. 0) then
       stagedgepphi = stagedgepphi * (1 + side * 1e-10)
    else
       stagedgepphi = stagedgepphi * (1 - side * 1e-10)
    end if
  end function stagedgepphi
  
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

  real function coplost(mub0, pphi)
    ! find the cop-confined/lost boundary
    ! which is the same as traplost, so just call it
    implicit none
    
    real, intent(in) :: mub0, pphi
    
    coplost = traplost(mub0, pphi)

  end function coplost

  real function coplostpphi(mub0, ee)
    ! find the cop lost pphi boundary for given mub0 and ee
    
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, ee
    real :: pphiout
    
    pphiout = 2. * (ee - (mub0 * (1-eps))) * mi * R0**2 * (1+eps)**2
    pphiout = sqrt(pphiout) - ei * psi1
    coplostpphi = pphiout

  end function coplostpphi
    
  real function ctplost(mub0, pphi)
    ! find the ct-passing-confined/lost energy boundary for given pphi and mub0 

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    real :: eeout

    eeout = mub0 * (1 + eps) + (pphi + ei * psi1)**2/mi/2./R0**2/(1.-eps)**2
    ctplost = eeout

  end function ctplost

  real function ctplostpphi(mub0, ee)
    ! find the ctp lost boundary pphi for given energy and mub0
    
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, ee
    real :: pphiout

    pphiout = 2.* (ee - mub0 * (1+eps)) * mi * R0**2 * (1-eps)**2
    pphiout = - sqrt(pphiout) - ei * psi1
    ctplostpphi = pphiout

  end function ctplostpphi

    real function losttpbound(mub0, istat)
    ! find the pphi value where tpbound crosses traplost
    use paras_phy, only : ei
    use paras_num, only : errtpbound, maxittpbound
    use profile, only : psi1
    implicit none

    real, intent(in) :: mub0
    integer, intent(out) :: istat

    integer :: i1, typemid, ifound
    real :: pphileft, pphiright, pphimid, tpboundtmp, losttmp

    pphileft = - ei * psi1
    pphiright = 0.
    
    i1 = 1
    ifound = 0
    do while ((ifound .le. 0) .and. (i1 .le. maxittpbound))
       ! binary search
       pphimid = (pphileft + pphiright) / 2.
       tpboundtmp = tpbound(mub0, pphimid, typemid)
       losttmp = traplost(mub0, pphimid)
 
       if (typemid .le. 0) then
          write(*,*) 'tpbound return error in type12tpbound'
          i1 = maxittpbound
          exit
       end if

       if (abs((pphimid - pphileft)/pphimid) .le. errtpbound) then
          ! precision reached
          ifound = 1
       end if
       
       if (tpboundtmp .ge. losttmp) then
          pphileft = pphimid
       else
          pphiright = pphimid
       end if
          
       i1 = i1 + 1
       
    end do

    if (ifound .le. 0) then
       write(*,*) 'losttpbound failed'
       istat = -1
       losttpbound = -100
    else
       istat = 1
       losttpbound = pphimid
    end if
    
  end function losttpbound
  
  real function eeonaxis(mub0, pphi)
    ! find the co-passing edge

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi
    real :: eeout

    eeout = mub0 + (pphi)**2/mi/2./R0**2
    eeonaxis = eeout

  end function eeonaxis

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

  real function rvt(r, mub0, ee, vsign)
    ! R * r (dtheta/dt) given r, mub0, ee and sign of vpar
    use paras_phy
    use profile
    implicit none
    
    real, intent(in) :: r, mub0, ee
    integer, intent(in) :: vsign
    real :: mub, Rbar, vpar, side
    
    side = float(vsign)
    mub = mub0 * (1 - r * eps * side)
    Rbar = 1 + r * eps * side
    vpar = sqrt((2*ee - 2*mub) / mi)

    rvt =  mi * vpar**2 + mub0 * Rbar
    rvt = - rvt / ei / B0 * side  + side * a * r * vpar / q(r)
  end function rvt

  real function drvtdr(r, mub0, ee, vsign)
    ! r derivative of rvt
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: r, mub0, ee
    integer, intent(in) :: vsign
    real :: mub, Rbar, vpar, side
    
    side = float(vsign)
    mub = mub0 * (1 - r * eps * side)
    Rbar = 1 + r * eps * side
    vpar = sqrt((2*ee - 2*mub) / mi)
    
    drvtdr = -3. * mub0 / B0 / ei / R0 
    drvtdr = drvtdr + a * side * vpar * (1 - r * dqdr(r) / q(r)) / q(r)
    drvtdr = drvtdr + a * r * eps * mub0 / mi / vpar / q(r)
    
  end function drvtdr
    
end module orbit_classify
