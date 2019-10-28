! calculate m3 for counter-passing particles - mostly copies trap module

module ctp_matrix
!
!  TYPE:
!  ctpmatrix           - a bundle of ctpgrid with different mub0
!
!  SUBROUTINE:
!  (public)
!  ctpmatrix_init      - initialize ctpmatrix and its ctpgrids, must call first
!  ctpmatrix_calculate - calculate the period and finite element weight for each
!                      ctpgrid, must call before getmat3
!  ctpmatrix_destroy   - deallocate everything, must call before program ends
!  getmat3trap         - get M3 for given omega
!  pressurectp         - get fast ion energy over total volume for pressure calcs
!  pressurectp2        - redid pressure calc for debugging, higher by x4




  use mpi
  use paras_phy
  use distribution_fun
  use radial_grid, only : nelement, nele
  use ctp_grid
  use matrix_module
  use landau_integral
  implicit none

  private

  type, public :: ctpmatrix
    ! type containing bundle of ctp grid

    ! contains mub0 grid
    type(ctpgrid), allocatable, dimension(:), public :: grid
    ! number of mub0 grid points
    integer, public :: ngrid
    ! grid mub0
    real, allocatable, dimension(:) :: mub0_table

  end type ctpmatrix

  ! workload allocation for parallel computing
  type(workload), public :: lwork

  public :: ctpmatrix_init, ctpmatrix_calculate, ctpmatrix_destroy, getmat3ctp, pressurectp2 !, slowingdownnorm , ctpmatrix_bcast

contains

  subroutine ctpmatrix_init(this, ngrid, mub0start, mub0end, lambda0, dlambda, npphin, neen, np, ibroadcast, ierr)
    ! initiate grid bundle
    ! INPUT : number of mub0 grid points, mub0 grid start/end points,
    !         number of pphi grid, number of energy grid, 
    !         number of harmonics
    !         ibroadcast - LOGICAL, if the periods and Vp need to be broadcasted to all nodes (not currently used)


    use sintable, only : sinpzeta
    use radial_grid
    implicit none

    type(ctpmatrix) :: this
    integer, intent(in) :: ngrid, npphin, neen, np
    real, intent(in) :: mub0start, mub0end, lambda0, dlambda
    logical, intent(in) :: ibroadcast
    integer, intent(out) :: ierr

    real :: dmub0, mub0
    integer :: i1

    ierr = 0

    if (ngrid .lt. 2) then
      write(*,*) 'ERROR in ctpmatrix_init : ngrid must >= 3'
      ierr = 1
      return
    end if

    if ((mub0start .le. 0.) .or. (mub0end .le. 0) &
        .or. (mub0start .ge. mub0end)) then
      write(*,*) 'ERROR in ctpmatrix_init : illegal mub0 start or end'
      write(*,*) mub0start, mub0end
      ierr = 1
      return
    end if

    if (npphin .le. 4) then
      write(*,*) 'ERROR in ctpmatrix_init : npphi must > 4'
      ierr = 1
      return
    end if

    if (neen .le. 4) then
      write(*,*) 'ERROR in ctpmatrix_init : neen, neeb must > 4'
      ierr = 1
      return
    end if
   
    if (np .lt. 1) then
      write(*,*) 'ERROR in ctpmatrix_init : np must >= 1'
      ierr = 1
      return
    end if

    if (.not. allocated(sinpzeta)) then
      write(*,*) 'ERROR in ctpmatrix_init, sin table must be initialized first'
      ierr = 1
      return
    end if

    if ((.not. allocated(rgrid)) .or. (nelement .le. 0)) then
      write(*,*) 'ERROR in ctpmatrix_init, radial grid must be initialized first'
      ierr = 1
      return
    endif

    call mpi_allocate_work(lwork, ngrid)

    this%ngrid = ngrid
    dmub0 = (mub0end - mub0start) / real(ngrid - 1)

    allocate(this%grid(ngrid))
    allocate(this%mub0_table(ngrid))

    do i1 = 1, ngrid
      mub0 = mub0start + dmub0 * real(i1-1)
      this%mub0_table(i1) = mub0

      if (mpi_is_my_work(lwork, i1)) then
        call ctpgrid_init(this%grid(i1), mub0, lambda0, dlambda, npphin, neen, np)
      end if
    end do

  end subroutine ctpmatrix_init

  ! subroutine ctpmatrix_bcast(this)
  !   ! broadcast the matrics to all nodes
  !   implicit none

  !   type (ctpmatrix) :: this
  !   integer :: i1, i2 ,i3, i4, whose

  !   do i1 = 1, this%ngrid
  !     whose = mpi_whose_work(lwork, i1)
  !     if (this%grid(i1)%npphin .le. 0) cycle ! no grid then skip
  !     ! broadcast the n grids
  !     do i2 = 1, this%grid(i1)%npphin
  !       call mpi_bcast_spline(this%grid(i1)%periodn(i1), whose)
  !       do i3 = 1, nele
  !         do i4 = 1, this%grid(i1)%np
  !            call mpi_bcast_spline(this%grid(i1)%vpmgridn(i3,i4,i1), whose)
  !         end do
  !       end do
  !     end do

  !   end do
  ! end subroutine ctpmatrix_bcast

  subroutine ctpmatrix_calculate(this)
    ! calculate orbit periods and element weights
    use paras_phy, only : eunit
    implicit none

    type(ctpmatrix), intent(in) :: this
    integer :: i1

    do i1 = 1, this%ngrid
      if (mpi_is_my_work(lwork, i1)) then
        write (*,100) mpi_get_rank(), this%grid(i1)%mub0/eunit/1000., i1, this%ngrid
        call ctpgrid_calculate(this%grid(i1))
      end if
    end do

    100 format('on cpu', I2, ': computing orbits for mu B0 =', F8.1, 'keV, ', I4, ' of', I4)

  end subroutine ctpmatrix_calculate

  subroutine ctpmatrix_destroy(this, ibroadcast)
    ! destroy all ctp grids
    ! ibroadcast - LOGICAL, if to broadcast to all nodes (not used currently)
    implicit none

    type(ctpmatrix) :: this
    logical, intent(in) :: ibroadcast
    integer :: i1


    if (allocated(this%mub0_table)) deallocate(this%mub0_table)
    
    do i1 = 1, this%ngrid
       if (mpi_is_my_work(lwork, i1)) then
          call ctpgrid_destroy(this%grid(i1))
       end if
    end do
  end subroutine ctpmatrix_destroy

  subroutine getmat3ctp(this, omega, mat3)
    ! get the value of matrix 3 for the given frequency
    ! INPUT:  frequency
    ! OUTPUT: mat3
    use profile
    implicit none

    type(ctpmatrix), intent(in) :: this
    complex, intent(in) :: omega
    type(matrix), intent(out) :: mat3

    integer :: i1, i2, m, n, p
    type(matrix) :: tmpmat, tmpmat2
    real :: dmub0, dpphi, s
    complex :: tmp

    call matrix_init(mat3, 2*nelement-2, 2*nelement-2)
    call matrix_init(tmpmat, 2*nelement-2, 2*nelement-2)
    call matrix_init(tmpmat2, 2*nelement-2, 2*nelement-2)
    do i1 = 1, this%ngrid
      ! mub0 level

      ! only do the allocated work
      if (.not. mpi_is_my_work(lwork, i1)) cycle
      
      ! clear the tmp matrix
      call matrix_clear(tmpmat2)

      do i2 = 1, this%grid(i1)%npphin
         ! pphi level
         ! clear the tmp matrix
         call matrix_clear(tmpmat)
         
         do m = 1, 2 * nelement - 2
            do n = m, 2 * nelement - 2
               do p = 1, this%grid(i1)%np
                  if (this%grid(i1)%periodn(i2)%d(1) .eq. 0.) then
                     ! fqc constant for given mub0 and pphi, get int will fail
                     tmpmat%data(n, m) = tmpmat%data(n, m) &
                          + getintnormal(this, i1, i2,  omega, p, m, n) &
                          + getintnormal(this, i1, i2, -omega, p, m, n)
                  else
                     tmpmat%data(n, m) = tmpmat%data(n, m) &
                          + getint(this, i1, i2,  omega, p, m, n) &
                          + getintnormal(this, i1, i2, -omega, p, m, n)
                  end if
               end do
               tmpmat%data(m, n) = tmpmat%data(n, m)
            end do
         end do
         
         s = this%grid(i1)%s(i2) 
         dpphi = abs(this%grid(i1)%ds * sgridds(s))
         !write(*,*) s, dpphi
         tmpmat2%data(:,:) = tmpmat2%data(:,:)+ cmplx(dpphi) * tmpmat%data(:,:)
      end do

      if (i1 .eq. 1) then
         dmub0 = this%mub0_table(2) - this%mub0_table(1)
         call matrix_clear(mat3)
      else if (i1 .eq. this%ngrid) then
         dmub0 = this%mub0_table(i1) - this%mub0_table(i1-1)
      else
         dmub0 = this%mub0_table(i1+1) - this%mub0_table(i1-1)
      end if
      mat3 = cmplx(0.5*dmub0) * tmpmat2 + mat3
    end do

   mat3 = cmplx(- nf_ratio * a**2 * ei * pi**2 / mi**2 / mib / R0 /B0) * mat3

   !sum over all cpus
   call matrix_clear(tmpmat2)
   call mpi_sum_matrix(mat3, tmpmat2, 0)
   ! transfer the data back to mat3
   mat3 = tmpmat2
   !mat3%data = real(mat3%data)
   
   call matrix_destroy(tmpmat)
   call matrix_destroy(tmpmat2)

  end subroutine getmat3ctp

  ! subroutine pressurectp(this, ph, pth)
  !   use paras_phy
  !   use profile, only : ti0
  !   use distribution_fun, only : nf_ratio
  !   implicit none

  !   type(ctpmatrix), intent(in) :: this
  !   real, intent(out) :: ph, pth

  !   real :: dmub0, mub00, mub0f, result
  !   integer :: i1


  !   ! pressure = (2pi)^3/em^2 * Int{ ee*f0(ee,mub0,pphi)*(1/wb)*dE*du*dpphi }
  !   ! Use trapezoidal rule:
  !   !  b
  !   ! S  f(x)dx ~ [ Sum[ f(xk) ] + ( f(x0)+f(xf) )/2 ] * dx
  !   !  a
  !   ! Use Fubini's theorem to do this 3 times for a 3D integral over (E,mub0,pphi)

  !   ! mub0 grid spacing
  !   dmub0 = (this%mub0_table(this%ngrid) - this%mub0_table(1))/real(this%ngrid-1)

  !   ! First and last mub0 values
  !   mub00 = getpphiint(this, 1)
  !   mub0f = getpphiint(this, this%ngrid)
  !   result = (mub00 + mub0f)/2.

  !   ! Loop over other mub0 values
  !   do i1 = 2, this%ngrid-1
  !     result = result + getpphiint(this, i1)
  !   end do

  !   ! p_fast = nf_ratio * (2pi)^3 / e m^2 * Sum(...)
  !   ph = nf_ratio * (2.**3. * pi**3. / mi**2. / ei) * dmub0 * result / B0
  !   ! p_bulk = volume * ti0, ti0 in Joules
  !   pth = (pi * a**2. * 2.*pi* R0) * ti0

  ! end subroutine pressurectp

  subroutine pressurectp2(this, ph, pth)
    use paras_phy
    use profile, only : ti0, te0, psi1
    use distribution_fun, only : nf_ratio, f0
    implicit none

    type(ctpmatrix), intent(in) :: this
    real, intent(out) :: ph, pth

    real :: result, mub0, ee, pphi, dmub0, dpphi, dee, w2
    integer :: i1,i2,i3, nee


    dmub0 = (this%mub0_table(this%ngrid) - this%mub0_table(1))/real(this%ngrid-1)

    result = 0
    do i1 = 1, this%ngrid
      mub0 = this%mub0_table(i1)
      ! if ((i1 .eq. 1) .or. (i1 .eq. this%ngrid)) then 
      !   w = dmub0/2.
      ! else
      !   w = dmub0
      ! end if

      dpphi = this%grid(i1)%ds * ei * psi1 ! pphi grid spacing

      do i2 = 1, this%grid(i1)%npphin ! Number of pphi grid points in current grid
        pphi = this%grid(i1)%pphigrid(i2)
        ! if ((i2 .eq. 1) .or. (i2 .eq. this%ngrid)) then 
        !   w = w*dpphi/2.
        ! else
        !   w = w*dpphi
        ! end if

        nee = this%grid(i1)%periodn(i2)%n ! Number of energy points
        dee = (this%grid(i1)%periodn(i2)%x(nee) - this%grid(i1)%periodn(i2)%x(1)) / (nee-1) ! Energy spacing

        do i3 = 1, nee
          ee = this%grid(i1)%periodn(i2)%x(i3)
          ! if ((i3 .eq. 1) .or. (i3 .eq. nee)) then 
          !   w = w*dee/2.
          ! else
          !   w = w*dee
          ! end if

          w2 = dmub0*dpphi*dee
          if ((i1 .eq. 1) .or. (i1 .eq. this%ngrid)) then
            w2 = w2/2.
          end if
          if ((i2 .eq. 1) .or. (i2 .eq. this%grid(i1)%npphin)) then
            w2 = w2/2.
          end if
          if ((i3 .eq. 1) .or. (i3 .eq. nee)) then
            w2 = w2/2.
          end if


          result = result + w2 * f0(ee,mub0,pphi) * ee * this%grid(i1)%periodn(i2)%y(i3) / 2. / pi
        end do
      end do

    end do

    ! p_fast = nf_ratio * (2pi)^3 / e m^2 * result
    ph = nf_ratio * (2.**3. * pi**3. / mi**2. / ei) * result / B0
    ! p_bulk = volume * ti0, ti0 in Joules
    pth = (pi * a**2. * 2.*pi* R0) * (ti0 + te0) / 2. ! Divided by 2 because of T(r) profile

  end subroutine pressurectp2


  ! real function slowingdownnorm(this)
  !   ! attempt to normalize f0 using grid. 
  !   use paras_phy
  !   use profile
  !   use distribution_fun, only : nf_ratio, f0
  !   implicit none

  !   type(ctpmatrix), intent(in) :: this

  !   real :: result, mub0, ee, dmub0, dpphi, dee, w2, eedenom, pitch, J, eec
  !   integer :: i1,i2,i3, nee

  !   eec = te0 * (3*pi**(1./2.) * zf * (2.014*mp)**(3./2.) /4./(me**(1./2.))/(2.014*mp))**(2./3.)

  !   dmub0 = (this%mub0_table(this%ngrid) - this%mub0_table(1))/real(this%ngrid-1)

  !   result = 0
  !   do i1 = 1, this%ngrid
  !     mub0 = this%mub0_table(i1)

  !     do i2 = 1, this%grid(i1)%npphin ! Number of pphi grid points in current grid
  !       nee = this%grid(i1)%periodn(i2)%n ! Number of energy points
  !       dee = (this%grid(i1)%periodn(i2)%x(nee) - this%grid(i1)%periodn(i2)%x(1)) / (nee-1) ! Energy spacing

  !       do i3 = 1, nee
  !         ee = this%grid(i1)%periodn(i2)%x(i3)

  !         w2 = dmub0 *dee ! * this%grid(i1)%ds*ei*psi1
  !         if ((i1 .eq. 1) .or. (i1 .eq. this%ngrid)) then
  !           w2 = w2/2.
  !         end if
  !         ! if ((i2 .eq. 1) .or. (i2 .eq. this%grid(i1)%npphin)) then
  !         !   w2 = w2/2.
  !         ! end if
  !         if ((i3 .eq. 1) .or. (i3 .eq. nee)) then
  !           w2 = w2/2.
  !         end if

  !         eedenom = ee**(3./2.) + eec**(3./2.)
  !         pitch = (mub0/ee - lambda0)**2.
  !         J = pi*B0/(mp**(3./2.)) * sqrt(2./(ee-mub0))
    
  !         result = result + w2 *J* exp(-pitch / dlambda**2.) / eedenom
  !       end do
  !     end do
  !   end do

  !   slowingdownnorm = result/B0
  ! end function slowingdownnorm



! ************ internal subroutines *************
  complex function getint(this, imub0, ipphi, omega, p, m, n)
    ! get the energy Landau integral for a particular mub0, pphi and p, 
    ! for given omega
    ! INPUT : mub0 index, ipphi index, complex frequency, harmonic index p,
    !         finite element indexes m and n
    ! OUTPUT: the value of the integral
    
    use cubic, only : getcoeffs
    use distribution_fun
    implicit none

    type(ctpmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n
    complex, intent(in) :: omega

    complex :: periodp, results, tmpres
    real :: x1, x2, y1, y2, y1p, y2p
    integer :: i1, i2, iup, ipphib
    real, dimension(4) :: c, pc
    real :: mub0, pphi, ee
    
    periodp = 2 * pi * real(p) / omega

    results = (0., 0.)
    ! upper energy index limit of n grid 
    iup = this%grid(imub0)%periodn(ipphi)%n - 1
    
    ! the value on the first grid point
    call getnumeratorn(this, imub0, ipphi, m, n, p, 1, y2, y2p)

    ! integral for the n grids
    do i1 = 1, iup
      x2 =  this%grid(imub0)%periodn(ipphi)%x(i1+1) &
            - this%grid(imub0)%periodn(ipphi)%x(i1)
      ! get value from the last step
      y1 = y2
      y1p = y2p
      call getnumeratorn(this, imub0, ipphi, m, n, p, i1+1, y2, y2p)
      if ((y1 .eq. 0.) .and. (y1p .eq. 0.)) then
          if ((y2 .eq. 0.) .and. (y2p .eq. 0.)) then 
            ! all zero, no need to continue
            cycle
          end if
          c(1) = 0.
          c(2) = 0.
          c(3) = y2p
          c(4) = y2 - c(3) * x2
          x1 = -c(4) / c(3)
          if (x1 .ge. x2) cycle
          if (x1 .le. 0.) x1 = 0.
      else if ((y2 .eq. 0.) .and. (y2p .eq. 0.)) then
          c(1) = 0.
          c(2) = 0.
          c(3) = y1p
          c(4) = y1
          x1 = -c(4) / c(3)
          if (x1 .ge. x2) x2 = x1
          if (x1 .le. 0.) cycle
          x1 = 0.
      else
          x1 = 0.
          call getcoeffs(x2, y1, y1p, y2, y2p, c)
      end if
      pc(1) = this%grid(imub0)%periodn(ipphi)%d(i1)
      pc(2) = this%grid(imub0)%periodn(ipphi)%c(i1)
      pc(3) = this%grid(imub0)%periodn(ipphi)%b(i1)
      pc(4) = this%grid(imub0)%periodn(ipphi)%a(i1)
      tmpres = landauint(c, pc, periodp, x1, x2)

      results = results + tmpres
    end do
    getint = results

  end function getint

  complex function getintnormal(this, imub0, ipphi, omega, p, m, n)
    ! get energy integral when no poles are presented
    ! INPUT:  mub0 index, pphi index, complex frequency, harmonic p, 
    !         finite element indexes m, n
    
    implicit none
    
    type(ctpmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n
    complex, intent(in) :: omega

    complex :: period, periodp, results, tmpres
    real :: dee, y2, y2p
    integer :: i1, i2, iup, ipphib
    
    periodp = 2 * pi * real(p) / omega

    results = (0., 0.)
    ! upper energy index limit of n grid 
    iup = this%grid(imub0)%periodn(ipphi)%n

    do i1 = 1, iup
      if (i1 .eq. 1) then
          dee = (this%grid(imub0)%periodn(ipphi)%x(2) - &
              this%grid(imub0)%periodn(ipphi)%x(1)) * 0.5
      else if (i1 .eq. iup) then
          dee = (this%grid(imub0)%periodn(ipphi)%x(i1) - &
              this%grid(imub0)%periodn(ipphi)%x(i1-1)) * 0.5
      else
          dee = (this%grid(imub0)%periodn(ipphi)%x(i1+1) - &
              this%grid(imub0)%periodn(ipphi)%x(i1-1)) * 0.5
      end if
      call getnumeratorn(this, imub0, ipphi, m, n, p, i1, y2, y2p)
      period = this%grid(imub0)%periodn(ipphi)%y(i1)
      results = results + y2 / (period - periodp) * dee
    end do
    getintnormal = results

  end function getintnormal

  subroutine getnumeratorn(this, imub0, ipphi, m, n, p, ipos, y, yp) 
    ! CONTAINS PHYSICS
    ! get the value of the numerator and its derivative in the energy integral
    ! for n grids
    ! INPUT:  mub0 index, pphi index, finite element indexes m and n,
    !         p the harmonic, the energy n grid index ipos
    ! OUTPUT: the numerator y and its derivative yp

    implicit none
    
    type(ctpmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n, ipos
    real, intent(out) :: y, yp

    real :: mub0, pphi, ee, vpm, vpn, vpmp, vpnp, dee, df0

    mub0 = this%grid(imub0)%mub0
    pphi = this%grid(imub0)%pphigrid(ipphi)
    ee   = this%grid(imub0)%periodn(ipphi)%x(ipos)

    vpn = this%grid(imub0)%vpmgridn(n, p, ipphi)%y(ipos)
    vpm = this%grid(imub0)%vpmgridn(m, p, ipphi)%y(ipos)
    
    if (ipos .eq. this%grid(imub0)%periodn(ipphi)%n) then
      dee = ee - this%grid(imub0)%periodn(ipphi)%x(ipos-1)
      vpnp = ((3. * this%grid(imub0)%vpmgridn(n, p, ipphi)%d(ipos-1) * dee) &
            + 2. * this%grid(imub0)%vpmgridn(n, p, ipphi)%c(ipos-1)) * dee &
            + this%grid(imub0)%vpmgridn(n, p, ipphi)%b(ipos-1)
      vpmp = ((3. * this%grid(imub0)%vpmgridn(m, p, ipphi)%d(ipos-1) * dee) &
            + 2. * this%grid(imub0)%vpmgridn(m, p, ipphi)%c(ipos-1)) * dee &
            + this%grid(imub0)%vpmgridn(m, p, ipphi)%b(ipos-1)
      
    else
      vpnp = this%grid(imub0)%vpmgridn(n, p, ipphi)%b(ipos)
      vpmp = this%grid(imub0)%vpmgridn(m, p, ipphi)%b(ipos)
    end if
      
    df0 = df0de(ee, mub0, pphi)
    y = vpm * vpn
    ! chain rule for derivative
    yp = vpnp * vpm + vpmp * vpn
    yp = yp * df0 +  y * d2f0de2(ee, mub0, pphi) 
    y = y * df0

  end subroutine getnumeratorn

  ! real function getpphiint(this, imub0)
  !   ! Used in calculation of fast ion pressure - 2nd layer of trapezoidal rule
  !   ! Calculate integral over pphi given a specific value of mub0
  !   use paras_phy, only : ei
  !   use profile, only : psi1
  !   implicit none

  !   type(ctpmatrix), intent(in) :: this
  !   integer, intent(in) :: imub0

  !   real :: pphi0, pphif, result
  !   integer :: npphi, i1

  !   npphi = this%grid(imub0)%npphin ! Number of pphi grid points in current grid

  !   pphi0 = geteeint(this, imub0, 1)
  !   pphif = geteeint(this, imub0, npphi)
  !   result = (pphi0 + pphif)/2.

  !   do i1 = 2, npphi-1
  !     result = result + geteeint(this, imub0, i1)
  !   end do

  !   getpphiint  = result * this%grid(imub0)%ds * ei * psi1 ! pphi grid spacing

  ! end function getpphiint

  ! real function geteeint(this, imub0, ipphi)
  !   ! Used in calculation of fast ion pressure - 3rd layer of trapezoidal rule
  !   ! Calculate integral over ee given a specific value of mub0 and pphi
  !   use distribution_fun, only : f0
  !   use paras_phy, only : pi
  !   implicit none

  !   type(ctpmatrix), intent(in) :: this
  !   integer, intent(in) :: imub0, ipphi

  !   real :: mub0, pphi, dee, ee0, t0, eef, tf, result, ee, t
  !   integer :: nee, i1
  !   type(ctpgrid) :: ctpg

  !   mub0 = this%mub0_table(imub0) ! Current mub0 point
  !   ctpg = this%grid(imub0) ! Current ctp grid
  !   pphi = this%grid(imub0)%pphigrid(ipphi) ! Current pphi value

  !   nee = ctpg%periodn(ipphi)%n ! Number of energy points
  !   dee = (ctpg%periodn(ipphi)%x(nee) - ctpg%periodn(ipphi)%x(1)) / (nee-1) ! Energy spacing

  !   ee0 = ctpg%periodn(ipphi)%x(1)
  !   t0 = ctpg%periodn(ipphi)%y(1)

  !   eef = ctpg%periodn(ipphi)%x(nee)
  !   tf = ctpg%periodn(ipphi)%y(nee)

  !   result = (f0(ee0,mub0,pphi)*ee0*t0 /2./pi + f0(eef,mub0,pphi)*eef*tf /2./pi) / 2.

  !   do i1 = 2, nee-1
  !     ee = ctpg%periodn(ipphi)%x(i1) ! Current energy value
  !     t = ctpg%periodn(ipphi)%y(i1) ! Current period value

  !     result = result + f0(ee,mub0,pphi)*ee*t / 2. / pi
  !   end do
  !   geteeint = result*dee

  ! end function geteeint


end module ctp_matrix