! calculate the drift orbit

module orbit

  implicit none
  ! subroutine to get orbit data
  public :: getorbit
  ! private subroutines to find the start point and perform Runge-Kutta method
  private :: findrstart, getee, getprime, rkstep, rkstepcat
  ! dr/dt and r d theta/dt in drift orbit calculation
  private :: drdt, rdthetadt

contains

  subroutine getorbit(ee, mub0, pphi, norbit, vsignin, orbit_dt, r, theta, period, istatus)
    ! Inputs : energy, mu B0, initial r, sign of initial v_parallel, time step (in unit of cycotron frequency) 
    ! Outputs: orbit r and theta, orbit period, orbit status
    ! istatus =  lost(2)/success(1)/failed(<=-100)

    use paras_num
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi, orbit_dt
    integer, intent(in) :: norbit, vsignin
    real, intent(out) :: period
    integer, intent(out) :: istatus
    real, dimension(norbit), intent(out) :: r, theta

    real, dimension(nmaxts) :: rtemp, ttemp ! temp orbit data
    real, dimension(100) :: rsmall, tsmall
    integer :: vsign
    real :: rstart, tstart
    real :: newr, newtheta, dt, dtlast, dttemp, yout, ylast
    real :: ttarget, curtime, timediff
    integer :: i, j, ipos, flagtrapped


    ! find the start point
    rstart = findrstart(ee, mub0, pphi, vsignin)

    ! check if the orbit exists
    if ((rstart .le. -rlost) .or. (ee .le. mub0 * (1-eps*rstart))) then
       ! energy must be greater than mub0 * (1 - r/R0)
       istatus = -100
       write(*,*) 'ERROR: orbit does not exist at (e,mub0,r) = ', &
            ee/eunit, mub0/eunit,rstart
       return
    end if

    if (rstart .lt. 0.) then
       ! start from the high field side
       tstart = pi
       rstart = abs(rstart)
    else
       ! start from the low field side
       tstart = 0.
    end if

    ! if the start point is outside the plasma, return "lost"
    if (rstart .gt. rlost) then
       istatus = 2
       return 
    end if

    ! time step = orbit_dt * cycotron period
    dt = orbit_dt * Tc
    ! start point : high field side on mid-plane
    rtemp(1) = rstart
    ttemp(1) = tstart
    vsign = vsignin
    flagtrapped = 0

    ! Runge-Kutta method to solve the ODE, maximum nmaxts steps
    i = 1
    do while (i .lt. nmaxts)

       if (rtemp(i) .le. rcartesian) then
          ! if too close to the magnetic axis, switch to Cartesian (slower)
          call rkstepcat(ee, mub0, pphi, dt, rtemp(i), ttemp(i),&
               rtemp(i+1), ttemp(i+1), yout)
       else
          ! otherwise use cylindrical (faster)
          call rkstep(ee, mub0, pphi, dt, rtemp(i), ttemp(i),&
               rtemp(i+1), ttemp(i+1))
       end if
!!$       write(20,*) rtemp(i+1)*cos(ttemp(i+1)), rtemp(i+1)*sin(ttemp(i+1))

       if (rtemp(i+1) .le. 0.) rtemp(i+1) = 1e-10
       ! check if the orbit is lost
       if (rtemp(i+1).gt.rlost) then
          istatus = 2
          write(*,*) 'ORBIT LOST', rstart, rtemp(i+1)
          return
       end if

       ! check if the orbit is closed (half a orbit)
       if ((ttemp(i) * ttemp(i+1)).lt.0) then
          ! orbit is closed
          ttarget = 0.
          istatus = 1
          ! terminate loop
          exit
       else if ((ttemp(i) -  pi) * (ttemp(i+1) -  pi) .lt. 0) then
          ! orbit is closed
          ttarget =  pi
          istatus = 1
          ! terminate loop
          exit
       else if ((ttemp(i) +  pi) * (ttemp(i+1) +  pi) .lt. 0) then
          ! orbit is closed
          ttarget = - pi
          istatus = 1
          ! terminate loop
          exit
       end if

       ! proceed to next time step
       i = i + 1
    end do

    if (i.eq.nmaxts) then
       ! the orbit is not closed after maximum time steps, return error
       istatus = -200
       return
    end if

    ! calculate poloidal orbit frequency
    ! find time step size at the latest step until orbit perfectly closed
    dtlast = dt
    j = 0
    if (rtemp(i) .le. rcartesian) then
       ylast = rtemp(i)*sin(ttemp(i))
       yout  = rtemp(i+1) * sin(ttemp(i+1))
       do while (abs(yout) .gt. errtheta)          
          dtlast = dtlast /  (yout - ylast) * (-ylast)
          call rkstepcat(ee, mub0, pphi, dtlast, rtemp(i), ttemp(i),&
               rtemp(i+1), ttemp(i+1), yout)
          j = j + 1
          ! give up if exceeding the limit
          if (j .ge. nmaxorbitclose) exit
       end do
    else
       do while (abs(ttemp(i+1)-ttarget) .gt. errtheta) 
          dtlast = dtlast / (ttemp(i+1)-ttemp(i))*(ttarget-ttemp(i))
          call rkstep(ee, mub0, pphi, dtlast, rtemp(i), ttemp(i),&
               rtemp(i+1), ttemp(i+1))
!!$          write(*,*) rtemp(i+1)-rstart, ttemp(i+1)-ttarget
          j = j + 1
          ! give up if exceeding the limit
          if (j .ge. nmaxorbitclose) exit
       end do
    end if
    period = (real(i-1) * dt + dtlast)
    ! put orbit data on equal-distant t/T grid (action-angle)
    do j = 1, norbit / 2 + 1
       ! new grid point
       curtime = real(j-1)/real(norbit / 2) * period
       ! find where 'curtime' is on the old grid with time step dt
       ipos = floor(curtime / dt) + 1
       ! distance from the closet old grid point
       timediff = curtime - real(ipos - 1) * dt
       if (ipos .eq. i) then
          ! if 'curtime' is between the last two time steps
          dttemp = dtlast
       else
          dttemp = dt
       end if
       ! linear approximation between data points
       r(j) = rtemp(ipos) + (rtemp(ipos+1) - rtemp(ipos)) / dttemp * timediff
       theta(j) = ttemp(ipos) + (ttemp(ipos+1) - ttemp(ipos)) / dttemp * timediff
    end do
    ! extend it to another half of the orbit
    do j = norbit / 2 + 2, norbit
       r(j) = r(norbit - j + 2)
       theta(j) = 2. * ttarget - theta(norbit - j + 2)
    end do
    ! we have only calculated half a orbit
    ! the actual period is twice as large
    period = period * 2
  end subroutine getorbit

  subroutine rkstep(ee, mub0, pphi, dt, rin, thetain, rout, thetaout)
    ! Runge-Kutta step, cylindrical coordinate (faster)
    use paras_phy
    implicit none

    real, intent(in) :: ee, mub0, pphi, rin, thetain, dt
    real, intent(out) :: rout, thetaout

    real :: k1r, k2r, k3r, k4r
    real :: k1t, k2t, k3t, k4t
    real :: newr, newtheta
    
    k1r = drdt(ee, mub0, rin, thetain)
    k1t = rdthetadt(ee, mub0, pphi, rin, thetain) / rin

    newr = rin + k1r * dt / 2.
    newtheta = thetain + k1t * dt / 2.
    k2r = drdt(ee, mub0, newr, newtheta)
    k2t = rdthetadt(ee, mub0, pphi, newr, newtheta) / newr

    newr = rin + k2r * dt / 2.
    newtheta = thetain + k2t * dt / 2.
    k3r = drdt(ee, mub0, newr, newtheta)
    k3t = rdthetadt(ee, mub0, pphi, newr, newtheta) / newr

    newr = rin + k3r * dt
    newtheta = thetain + k3t * dt
    k4r = drdt(ee, mub0, newr, newtheta)
    k4t =  rdthetadt(ee, mub0, pphi, newr, newtheta) / newr

    ! evolve the current time step and proceed
    rout     = rin     + dt/6.*(k1r + 2.*k2r + 2.*k3r + k4r)
    thetaout = thetain + dt/6.*(k1t + 2.*k2t + 2.*k3t + k4t)

  end subroutine rkstep

  subroutine rkstepcat(ee, mub0, pphi, dt, rin, thetain, rout, thetaout, yout)
    ! Runge-Kutta step, Cartesian coordinate (slower)
    use paras_phy
    implicit none

    real, intent(in) :: ee, mub0, pphi, rin, thetain, dt
    real, intent(out) :: rout, thetaout, yout

    real :: k1x, k2x, k3x, k4x
    real :: k1y, k2y, k3y, k4y
    real :: xin, yin, newx, newy, dr, rdtheta, newr, newtheta, thetaexp
    
    xin = rin * cos(thetain)
    yin = rin * sin(thetain)
    
    if (xin .ge. 0.) then
       thetaexp = 0.
    else
       thetaexp = pi
    end if

    dr = drdt(ee, mub0, rin, thetain)
    rdtheta = rdthetadt(ee, mub0, pphi, rin, thetain)
    k1x = dr * cos(thetain) - sin(thetain)*rdtheta
    k1y = dr * sin(thetain) + cos(thetain)*rdtheta

    newx = xin + k1x * dt / 2.
    newy = yin + k1y * dt / 2.
    newr = sqrt(newx**2 + newy**2)
    if (newr .le. 1e-10) then
       newtheta = thetaexp
    else
       newtheta = atan2(newy, newx)
    end if

    dr = drdt(ee, mub0, newr, newtheta)
    rdtheta = rdthetadt(ee, mub0, pphi, newr, newtheta)
    k2x = dr * cos(newtheta) - sin(newtheta)*rdtheta
    k2y = dr * sin(newtheta) + cos(newtheta)*rdtheta

    newx = xin + k2x * dt / 2.
    newy = yin + k2y * dt / 2.
    newr = sqrt(newx**2 + newy**2)
    newtheta = atan2(newy, newx)
    if (newr .le. 1e-10) then
       newtheta = thetaexp
    else
       newtheta = atan2(newy, newx)
    end if

    dr = drdt(ee, mub0, newr, newtheta)
    rdtheta = rdthetadt(ee, mub0, pphi, newr, newtheta)
    k3x = dr * cos(newtheta) - sin(newtheta)*rdtheta
    k3y = dr * sin(newtheta) + cos(newtheta)*rdtheta

    newx = xin + k3x * dt
    newy = yin + k3y * dt
    newr = sqrt(newx**2 + newy**2)
    newtheta = atan2(newy, newx)
    if (newr .le. 1e-10) then
       newtheta = thetaexp
    else
       newtheta = atan2(newy, newx)
    end if

    dr = drdt(ee, mub0, newr, newtheta)
    rdtheta = rdthetadt(ee, mub0, pphi, newr, newtheta)
    k4x = dr * cos(newtheta) - sin(newtheta)*rdtheta
    k4y = dr * sin(newtheta) + cos(newtheta)*rdtheta

    ! evolve the current time step and proceed
    newx = xin + dt/6.*(k1x + 2.*k2x + 2.*k3x + k4x)
    newy = yin + dt/6.*(k1y + 2.*k2y + 2.*k3y + k4y)

    rout = sqrt(newx**2 + newy**2)
    thetaout = atan2(newy, newx)
    if (rout .le. 1e-10) then
       thetaout = thetaexp
    else
       thetaout = atan2(newy, newx)
    end if

    yout = newy

  end subroutine rkstepcat

  real function drdt(ee, mub0, r, theta)
    ! dr/dt of the ODE RHS

    use paras_phy
    implicit none

    real, intent(in) :: ee, mub0, r, theta
    real :: mub

    mub = mub0 - mub0 * eps * r * cos(theta)

    drdt = - (2. * (ee-mub) / (1 + eps * r * cos(theta)) + mub0) * sin(theta) * driftc 
  end function drdt

  real function rdthetadt(ee, mub0, pphi, r, theta)
    ! dtheta/dt of the ODE RHS

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, r, theta, pphi
    real ::  Rbar

    real :: ee1

    real :: mub, eemmub

    Rbar = 1 + eps * r * cos(theta)
    mub = mub0 - mub0 * eps * r * cos(theta)
    eemmub = ee - mub

    rdthetadt =  - (2.* eemmub / Rbar  + mub0) * cos(theta) * driftc 
    rdthetadt = rdthetadt + r * (pphi + ei * psi(r)) / mi / R0**2 / q(r) / Rbar**2

  end function rdthetadt

  real function findrstart(ee, mub0, pphi, vsign)

    use paras_phy
    use paras_num
    use profile
    use orbit_classify, only : copassingedge
    implicit none

    real, intent(in) :: ee, mub0, pphi
    integer, intent(in) :: vsign

    real :: r1, r2, rmid, ee1, ee2, eetemp
    real :: r3, r4, f3, f4, fmid
    integer :: i, ifound

    ! default, start searching from the low field side (r>0)
    if (pphi .ge. 0) then
       ! must be co-passing, search from 0 to 1 
       r1 = 0
       r2 = rlost
    else if (vsign .gt. 0) then
       ! must be co-passing / trapped, search from r(pphi) to 1 
       r1 = psitor(- pphi / ei)
       r2 = rlost
    else 
       if (pphi .ge. -psi1 * ei) then
          ! must be ct-passing / trapped, search from 0 to r(pphi)
          if (ee .ge. copassingedge(mub0, pphi)) then
             ! start searching from the high field side
             r1 =  -psitor(- pphi / ei)
             r2 = 0.
          else
             ! start searching from the low field side
             r1 = 0.
             r2 = psitor(-pphi / ei)
          end if
       else
          ! must be ct-passing / trapped, search from 0 to 1
          if (ee .ge. copassingedge(mub0, pphi)) then
             ! start searching from the high field side
             r1 = -rlost
             r2 =  0.
          else
             ! start searching from the low field side
             r1 = 0.
             r2 = rlost
          end if
       end if
    end if

    i = 0
    ifound = 0

    if (r1 .le. 0.) then
       ! high field side, use a Newton's method
       r3 = 0.
       eetemp = getee(mub0, pphi, r3)
       do while ((ifound .eq. 0) .and. (i .le. maxitrstart))
          f3 = getprime(mub0, pphi, r3)
          r4 = r3 - (eetemp - ee) / f3
          if (r4 .le. r1) then
             r4 = r1 / 2.
          end if

          r3 = r4
          eetemp = getee(mub0, pphi, r3)

          if (abs(eetemp - ee) / ee .le. errrstartmin) then
             ifound = 1
             exit
          end if

          i = i + 1
       end do

       rmid = r3
       
    else
       ee1 = getee(mub0, pphi, r1)
       ee2 = getee(mub0, pphi, r2)

       if ((ee1 - ee) * (ee2 - ee) .ge. 0) then
          ! need to find another search start point
          ! binary search for the min, or until getee < 0
          r3 = r1
          r4 = r2
          f3 = getprime(mub0, pphi, r3)
          f4 = getprime(mub0, pphi, r4)

          do while ((ifound .eq. 0) .and. (i .le. maxitrstart))
             rmid = (r3 + r4) / 2.
             eetemp = getee(mub0,pphi,rmid)
             if ((eetemp-ee) * (ee2-ee) .lt. 0) then
                ! a new start point is found
                ifound = 1
                r1 = rmid
                exit 
             end if

             fmid = getprime(mub0, pphi, rmid)
             if (fmid * f3 .le. 0) then
                r4 = rmid
                f4 = fmid
             else
                r3 = rmid
                f3 = fmid
             end if

             if (abs(fmid / ee) .le. errrstart) then
                ifound = 2
             end if
             i = i + 1
          end do

          ! check if the minimum is the solution
          eetemp = getee(mub0,pphi,rmid)
          if (abs(eetemp - ee) / ee .le. errrstartmin) then
             ifound = 1
          else
             if (ifound .eq. 1) then
                ifound = 0
                i = 0
                ee1 = getee(mub0,pphi,r1)
             else
                ! not found
                ifound = 0
                i = maxitrstart
             end if
          end if
       else
          i = 0
          ifound = 0
       end if


       ! binary search

       do while ((ifound .eq. 0) .and. (i .le. maxitrstart))
          rmid = (r1 + r2) / 2.
          eetemp = getee(mub0, pphi, rmid)
          if ((eetemp - ee) * (ee1 - ee) .le. 0) then
             r2 = rmid
             ee2 = eetemp
          else
             r1 = rmid
             ee1 = eetemp
          end if

          if (abs(eetemp - ee)/ee .le. errrstart) then
             ifound = 1
             exit
          end if
          i = i + 1
       end do

    end if
    
    if (ifound .gt. 0) then
       findrstart = rmid
    else
       ! not found
       findrstart = -2.
    end if

  end function findrstart

  real function getee(mub0, pphi, rtry)
    ! function used in findrstart
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi, rtry
    integer :: vsign

    real :: mub, vpar

    mub = mub0 * (1. - eps * rtry)
    vpar= (pphi + ei * psi(abs(rtry))) / mi / R0 / (1. + eps * rtry)
    getee = 0.5 * mi * vpar**2 + mub

  end function getee

  real function getprime(mub0, pphi, rtry)
    ! function used in findrstart, the derivative of getee with respect to r
    use paras_phy
    use profile
    implicit none

    real, intent(in) :: mub0, pphi, rtry

    real :: rbar

    rbar = (1 + eps * rtry)

    getprime = 2. * (pphi + ei * psi(abs(rtry))) * rtry * ei * B0 &
         / q(abs(rtry)) / rbar**2 
    getprime = getprime -  eps * (pphi + ei * psi(abs(rtry)))**2 &
         / rbar**3
    getprime = getprime / mi / R0**2 - mub0 * eps

  end function getprime

end module orbit
