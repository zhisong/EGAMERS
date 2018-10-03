! calculate m3 for trap particles

module trap_matrix
!
!  TYPE:
!  tmatrix           - a bundle of tgrid with different mub0
!
!  SUBROUTINE:
!  (public)
!  tmatrix_init      - initialize tmatrix and its tgrids, must call first
!  tmatrix_calculate - calculate the period and finite element weight for each
!                      tgrid, must call before getmat3
!  tmatrix_destroy   - deallocate everything, must call before program ends
!  getmat3trap       - get M3 for given omega
!  (private)
!  getint            - get energy integral
!

  use mpi
  use paras_phy
  use distribution_fun
  use radial_grid, only : nelement, nele
  use trap_grid
  use matrix_module
  use landau_integral
  implicit none

  private
  
  type, public :: tmatrix
     ! a type contains the bundle of trap grid

     ! contains the mub0 grid
     type(tgrid), allocatable, dimension(:), public :: grid
     ! number of mub0 grid
     integer, public :: ngrid
     ! grid mub0 (keep a record because of MPI)
     real, allocatable, dimension(:) :: mub0_table
     
  end type tmatrix

  ! workload allocation for parallel computing
  type(workload) :: lwork

  public :: tmatrix_init, tmatrix_calculate, tmatrix_destroy, &
            tmatrix_bcast, getmat3trap, getint, getintnormal
contains
  
  subroutine getmat3trap(this, omega, mat3)
    ! get the value of matrix 3 for the given frequency
    ! INPUT:  frequency
    ! OUTPUT: mat3
    use profile
    implicit none
    
    type(tmatrix), intent(in) :: this
    complex, intent(in) :: omega
    type(matrix), intent(out) :: mat3

    integer :: i1, i2, m, n, p
    type(matrix) :: tmpmat, tmpmat2
    real :: dmub0, dpphi
    complex :: tmp

    call matrix_init(mat3,   2*nelement-2, 2*nelement - 2)
    call matrix_init(tmpmat, 2*nelement-2, 2*nelement - 2)
    call matrix_init(tmpmat2,2*nelement-2, 2*nelement - 2)
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
          
          do m = this%grid(i1)%ielementmin(i2), this%grid(i1)%ielementmax(i2)
             do n = m, this%grid(i1)%ielementmax(i2)
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
          
          if (i2 .eq. 1) then
             dpphi = this%grid(i1)%pphigrid(2) - this%grid(i1)%pphigrid(1)
             tmpmat2 = cmplx(-0.5*dpphi) * tmpmat
          else if (i2 .eq. this%grid(i1)%npphin) then
             dpphi = this%grid(i1)%pphigrid(this%grid(i1)%npphin) & 
                  - this%grid(i1)%pphigrid(this%grid(i1)%npphin-1)
             tmpmat2 = cmplx(-0.5*dpphi) * tmpmat + tmpmat2
          else
             dpphi = this%grid(i1)%pphigrid(i2+1) - this%grid(i1)%pphigrid(i2-1)
             tmpmat2 = cmplx(-0.5 * dpphi) * tmpmat + tmpmat2
          end if
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
    
    call matrix_destroy(tmpmat)
    call matrix_destroy(tmpmat2)

  end subroutine getmat3trap

  subroutine tmatrix_init(this, ngrid, mub0start, mub0end, npphin, neen, neeb, eeendb, np, ibroadcast, ierr)
    ! initiate the grid bundle, must call first
    ! INPUT : number of mub0 grid point, start and end points
    !         mub0, number of pphi grid, number of n energy grid,
    !         number of b energy grid, b grid upper limit,
    !         number of harmonics
    !         ibroadcast - LOGICAL, if the periods and Vp need to be broadcasted to all nodes

    use sintable, only : sinpzeta
    use radial_grid
    implicit none
    
    type(tmatrix) :: this
    integer, intent(in) :: ngrid, npphin, neen, neeb, np
    real, intent(in) :: mub0start, mub0end, eeendb
    logical, intent(in) :: ibroadcast
    integer, intent(out) :: ierr

    real :: dmub0, mub0
    integer :: i1
    
    ierr = 0
    ! check if the parameters are illegal
    if (ngrid .lt. 2) then
       write(*,*) 'ERROR in tmatrix_init : ngrid must >= 3'
       ierr = 1
       return
    end if

    if ((mub0start .le. 0.) .or. (mub0end .le. 0) &
         .or. (mub0start .ge. mub0end)) then
       write(*,*) 'ERROR in tmatrix_init : illegal mub0 start or end'
       write(*,*) mub0start, mub0end
       ierr = 1
       return
    end if

    if (npphin .le. 4) then
       write(*,*) 'ERROR in tmatrix_init : npphi must > 4'
       ierr = 1
       return
    end if

    if ((neen .le. 4) .or. (neeb .le. 4)) then
       write(*,*) 'ERROR in tmatrix_init : neen, neeb must > 4'
       ierr = 1
       return
    end if
    
    if (np .lt. 1) then
       write(*,*) 'ERROR in tmatrix_init : np must >= 1'
       ierr = 1
       return
    end if
    
    if (eeendb .lt. 4.) then
       write(*,*) 'ERROR in tmatrix_init : eeendb must >= 4.'
       ierr = 1
       return
    end if

    if (.not. allocated(sinpzeta)) then
       write(*,*) 'ERROR in tmatrix_init, sin table must be initialized first'
       ierr = 1
       return
    end if

    if ((.not. allocated(rgrid)) .or. (nelement .le. 0)) then
       write(*,*) 'ERROR in tmatrix_init, radial grid must be initialized first'
       ierr = 1
       return
    endif
    
    ! allocate the workload for parallel computing
    ! parallelizing over muB0 grid, total grid points : ngrid
    call mpi_allocate_work(lwork, ngrid)
    
    this%ngrid = ngrid
    dmub0 = (mub0end - mub0start) / real(ngrid - 1)

    allocate(this%grid(ngrid))
    allocate(this%mub0_table(ngrid))
    do i1 = 1, ngrid
       mub0 = mub0start + dmub0 * real(i1-1)
       ! fill in the mub0 table
       this%mub0_table(i1) = mub0
       
       ! parallelizing over muB0 grid
       ! only allocate the grid if the current cpu needs to
       if (mpi_is_my_work(lwork, i1) .or. ibroadcast) then
          call tgrid_init(this%grid(i1), mub0, npphin, neen, neeb, eeendb, np)
       endif
       
    end do

  end subroutine tmatrix_init
  
  subroutine tmatrix_bcast(this)
    ! broadcast the matrics to all nodes
    implicit none

    type (tmatrix) :: this
    integer :: i1, i2 ,i3, i4, whose

    do i1 = 1, this%ngrid
      whose = mpi_whose_work(lwork, i1)
      if (this%grid(i1)%npphin .le. 0) cycle ! no grid then skip
      ! broadcast the n grids
      do i2 = 1, this%grid(i1)%npphin
        call mpi_bcast_spline(this%grid(i1)%periodn(i1), whose)
        do i3 = 1, nele
          do i4 = 1, this%grid(i1)%np
             call mpi_bcast_spline(this%grid(i1)%vpmgridn(i3,i4,i1), whose)
          end do
        end do
      end do

      ! broadcast the p grids
      do i2 = 1, this%grid(i1)%npphib
        call mpi_bcast_spline(this%grid(i1)%periodb(i1), whose)
        do i3 = 1, nele
          do i4 = 1, this%grid(i1)%np
             call mpi_bcast_spline(this%grid(i1)%vpmgridb(i3,i4,i1), whose)
          end do
        end do
      end do

      ! broadcast the max-min elements
      call mpi_bcast_array(this%grid(i1)%ielementmaxn, whose)
      call mpi_bcast_array(this%grid(i1)%ielementminn, whose)
      call mpi_bcast_array(this%grid(i1)%ielementmaxb, whose)
      call mpi_bcast_array(this%grid(i1)%ielementminb, whose)
      call mpi_bcast_array(this%grid(i1)%ielementmax, whose)
      call mpi_bcast_array(this%grid(i1)%ielementmin, whose)

    end do
  end subroutine tmatrix_bcast

  subroutine tmatrix_destroy(this, ibroadcast)
    ! destroy all tgrids, must call before program ends
    ! ibroadcast - LOGICAL, if to broadcast to all nodes
    implicit none
    
    type(tmatrix) :: this
    logical, intent(in) :: ibroadcast
    integer :: i1

    if (allocated(this%mub0_table)) deallocate(this%mub0_table)
    
    do i1 = 1, this%ngrid
       if (mpi_is_my_work(lwork, i1) .or. ibroadcast) then
          call tgrid_destroy(this%grid(i1))
       endif
    end do
  end subroutine tmatrix_destroy

  subroutine tmatrix_calculate(this)
    ! run the orbit period and element weight calculations
    ! must call before getmat3
    use paras_phy, only : eunit
    implicit none

    type(tmatrix), intent(in) :: this
    integer :: i1
    
    do i1 = 1, this%ngrid
       ! only perfrom the work when the current cpu needs to
       if (mpi_is_my_work(lwork, i1)) then
          write(*,100) mpi_get_rank(), this%grid(i1)%mub0/eunit/1000., i1, this%ngrid
          call tgrid_calculate(this%grid(i1))
       end if
    end do

100 format('on cpu', I2, ': computing orbits for mu B0 =', F8.1, 'keV, ', I4, ' of', I4)
  end subroutine tmatrix_calculate

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

    type(tmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n
    complex, intent(in) :: omega

    complex :: periodp, results, tmpres
    real :: x1, x2, y1, y2, y1p, y2p
    integer :: i1, i2, iup, ipphib
    logical :: ltype1
    real, dimension(4) :: c, pc
    real :: mub0, pphi, ee
    
    periodp = 2 * pi * real(p) / omega

    results = (0., 0.)
    ! upper energy index limit of n grid 
    iup = this%grid(imub0)%periodn(ipphi)%n - 1
    ! A TYPE I t/p boundary?
    ltype1 = istype1(this%grid(imub0), ipphi)
    ! ignore the last grid point if there is a TYPE I t/p boundary
    if (ltype1) iup = iup - 1
    
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

    ! integrate over b grid if there is any
    if (ltype1) then
       results = 0.
       ipphib = indexn2b(this%grid(imub0), ipphi) 
       ! the value on the first grid point
       call getnumeratorb(this, imub0, ipphi, m, n, p, 1, y2, y2p)
       iup = this%grid(imub0)%neeb - 1

       do i1 = 1, iup
          x2 =  this%grid(imub0)%periodb(ipphib)%x(i1+1) &
               - this%grid(imub0)%periodb(ipphib)%x(i1)
          ! get value from the last step
          y1 = y2
          y1p = y2p
          call getnumeratorb(this, imub0, ipphi, m, n, p, i1+1, y2, y2p)
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
             if (x1 .le. 0) x1 = 0.
          else if ((y2 .eq. 0.) .and. (y2p .eq. 0.)) then
             c(1) = 0.
             c(2) = 0.
             c(3) = y1p
             c(4) = y1
             x1 = -c(4) / c(3)
             if (x1 .le. 0.) cycle
             if (x1 .le. x2) x2 = x1
             x1 = 0.
          else
             x1 = 0.
             call getcoeffs(x2, y1, y1p, y2, y2p, c)
          end if
          pc(1) = this%grid(imub0)%periodb(ipphib)%d(i1)
          pc(2) = this%grid(imub0)%periodb(ipphib)%c(i1)
          pc(3) = this%grid(imub0)%periodb(ipphib)%b(i1)
          pc(4) = this%grid(imub0)%periodb(ipphib)%a(i1)
   
          tmpres = landauint(c, pc, periodp, 0., x2)
          results = results + tmpres

       end do
       getint = getint + results
    end if
  end function getint

  complex function getintnormal(this, imub0, ipphi, omega, p, m, n)
    ! get energy integral when no poles are presented
    ! INPUT:  mub0 index, pphi index, complex frequency, harmonic p, 
    !         finite element indexes m, n
    
    implicit none
    
    type(tmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n
    complex, intent(in) :: omega

    complex :: period, periodp, results, tmpres
    real :: dee, y2, y2p
    integer :: i1, i2, iup, ipphib
    logical :: ltype1
    
    periodp = 2 * pi * real(p) / omega

    results = (0., 0.)
    ! upper energy index limit of n grid 
    iup = this%grid(imub0)%periodn(ipphi)%n
    ! A TYPE I t/p boundary?
    ltype1 = istype1(this%grid(imub0), ipphi)
    ! ignore the last grid point if there is a TYPE I t/p boundary
    if (ltype1) iup = iup - 1

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

    ! integrate b grids if there is any
    if (ltype1) then
       results = 0.
       ipphib = indexn2b(this%grid(imub0), ipphi)
       do i1 = 1, this%grid(imub0)%neeb
          if (i1 .eq. 1) then
             dee = (this%grid(imub0)%periodb(ipphib)%x(2) - &
                  this%grid(imub0)%periodb(ipphib)%x(1)) * 0.5
          else if (i1 .eq. iup) then
             dee = (this%grid(imub0)%periodb(ipphib)%x(i1) - &
                  this%grid(imub0)%periodb(ipphib)%x(i1-1)) * 0.5
          else
             dee = (this%grid(imub0)%periodb(ipphib)%x(i1+1) - &
                  this%grid(imub0)%periodb(ipphib)%x(i1-1)) * 0.5
          end if
          call getnumeratorb(this, imub0, ipphi, m, n, p, i1, y2, y2p)
          period = this%grid(imub0)%periodb(ipphib)%y(i1)
          results = results + y2 / (period - periodp) * dee
       end do
       getintnormal = getintnormal + results

    end if
 
  end function getintnormal

  subroutine getnumeratorn(this, imub0, ipphi, m, n, p, ipos, y, yp) 
    ! CONTAINS PHYSICS
    ! get the value of the numerator and its derivative in the energy integral
    ! for n grids
    ! INPUT:  mub0 index, pphi index, finite element indexes m and n,
    !         p the harmonic, the energy n grid index ipos
    ! OUTPUT: the numerator y and its derivative yp

    implicit none
    
    type(tmatrix), intent(in) :: this
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

  subroutine getnumeratorb(this, imub0, ipphi, m, n, p, ipos, y, yp) 
    ! CONTAINS PHYSICS
    ! get the value of the numerator and its derivative in the energy integral
    ! for b grids near the trap/passing boundary
    ! INPUT:  mub0 index, pphi index, finite element indexes m and n,
    !         p the harmonic, the energy n grid index ipos
    ! OUTPUT: the numerator y and its derivative yp

    implicit none
    
    type(tmatrix), intent(in) :: this
    integer, intent(in) :: imub0, ipphi, p, m, n, ipos
    real, intent(out) :: y, yp

    real :: mub0, pphi, ee, eelog, vpm, vpn, vpmp, vpnp, dee, df0, ac
    integer :: ipphib

    ipphib = indexn2b(this%grid(imub0), ipphi)

    mub0 = this%grid(imub0)%mub0
    pphi = this%grid(imub0)%pphigrid(ipphi)
    eelog= this%grid(imub0)%periodb(ipphib)%x(ipos)

    vpn = this%grid(imub0)%vpmgridb(n, p, ipphib)%y(ipos)
    vpm = this%grid(imub0)%vpmgridb(m, p, ipphib)%y(ipos)
    
    if (ipos .eq. this%grid(imub0)%neeb) then
       dee = eelog - this%grid(imub0)%periodb(ipphib)%x(ipos-1)
       vpnp = ((3. * this%grid(imub0)%vpmgridb(n, p, ipphib)%d(ipos-1) * dee) &
            + 2. * this%grid(imub0)%vpmgridb(n, p, ipphib)%c(ipos-1)) * dee &
            + this%grid(imub0)%vpmgridb(n, p, ipphib)%b(ipos-1)
       vpmp = ((3. * this%grid(imub0)%vpmgridb(m, p, ipphib)%d(ipos-1) * dee) &
            + 2. * this%grid(imub0)%vpmgridb(m, p, ipphib)%c(ipos-1)) * dee &
            + this%grid(imub0)%vpmgridb(m, p, ipphib)%b(ipos-1)
       
    else
       vpnp = this%grid(imub0)%vpmgridb(n, p, ipphib)%b(ipos)
       vpmp = this%grid(imub0)%vpmgridb(m, p, ipphib)%b(ipos)
    end if
       
    ee = eelogtoee(this%grid(imub0), eelog, ipphi)
    df0 = df0de(ee, mub0, pphi)
    ac = mub0 * exp(-eelog)
    y = vpm * vpn 
    ! chain rule for derivative
    yp = vpnp * vpm + vpmp * vpn

    yp = (yp - y) * df0  +  y * d2f0de2(ee, mub0, pphi) * ac
    yp = yp * ac
    y = y * df0 * ac

  end subroutine getnumeratorb

end module trap_matrix
