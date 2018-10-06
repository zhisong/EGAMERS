! a grid for the trapped particles for a single mub0

module trap_grid
!
! "n" grids are Normal equal-distant energy grids 
!     from trap edge to t/p boundary or lost boundary
! "b" grid are special grids for TYPE I t/p Boundary
!     exponentially close to the t/p boundary
!
! TYPE:
! tgrid      - contains all the grids and splines
!
! SUBROUTINES:
! (public)
! tgrid_init      - allocate arrays and splines, fill energy grid
!                   must call first
! tgrid_destroy   - deallocate everything, must call before program ends
! tgrid_calcuate  - calcuate orbit period and finite element weight
!                   build splines, must call before spline data are used
! getperiod       - get orbit period and derivative for a given energy and pphi
! getperoidb      - get orbit period when close to t/p boundary
! findperiod      - given a orbit period, find all corresponding energy  
! (private)
! getperiodb1     - same, without checking if the input parameters are illegal
!
! FUNCTIONS:
! (public)
! getvpm          - get finite element weight for a given energy and pphi
! getvpmb         - get finite element weight when close to t/p boundary
! sgrid           - define the pphi grid (equal-distant or others)
! istype1         - logical, if the pphi index has a TYPE 1 t/p boundary
! indexn2b        - integer, convert type n pphi index to type b
! indexb2n        - integer, convert type b pphi index to type n
! eetoeelog       - real, convert energy to log scale close to t/p boundary
! eelogtoee       - real, reverse of eetoeelog
! (private)
! getvpmb1        - same, without checking if the input parameters are illegal
! gettgridtype    - integer, check if the given pphi and ee are out of bound

  use radial_grid, only : nele
  use spline_module
  implicit none

  ! type of trapped particle grid
  type, public :: tgrid

     ! finite element weight for n grid and for t/p boundary
     ! indexes : element label, pth resonance, pphi grid
     type(spline), dimension(:,:,:), allocatable, public :: vpmgridn, vpmgridb
     ! the most inside and most outside element indexes for a certain pphi
     integer, dimension(:), allocatable, public :: ielementmin, ielementmax
     integer, dimension(:,:), allocatable, public :: ielementminn, ielementmaxn, ielementminb, ielementmaxb

     ! orbit period grid for n and b
     type(spline), dimension(:), allocatable, public :: periodn, periodb

     ! number of pphi grid for n and b
     integer, public :: npphin, npphib
     ! start pphi index of b grid
     integer, public :: ibstart
     ! the lost boundary is lower than t/p boundary at
     integer, public :: iloststart
     
     ! pphi grid
     real, dimension(:), allocatable, public :: pphigrid
     ! pphi grid start and end
     real, public :: pphistart, pphiend
     integer, public :: ipphistart, ipphiend

     ! number of energy grid for n and b grid
     integer, public :: neen, neeb
     ! energy of n/b grid switch
     real, dimension(:), allocatable, public :: enbswitch

     ! store grid data for the lowest energy and t/p boundary
     real, dimension(:), allocatable, public :: etrapedge, etpbound, elost

     ! mub0
     real, public :: mub0

     ! harmonic count
     integer, public :: np


  end type tgrid

  private
  public :: tgrid_init, tgrid_destroy, tgrid_calculate, sgrid, getperiod, &
       getperiodb, findperiod, getvpm, getvpmb, istype1, indexn2b, indexb2n, &
       eetoeelog, eelogtoee

contains

  subroutine findperiod(this, period, ipphi, nfound, eelist, eetype)
    ! get the corresponding energy for a given orbit period and pphi grid index
    ! input: target period, pphi index
    ! output: number of energy points found and the energy itself
    !         eetype : 0 - on "n" grid, 1 - log scale on "b" grid
    use paras_num
    implicit none
    
    type(tgrid), intent(in) :: this
    real, intent(in) :: period
    integer, intent(in) :: ipphi
    integer, intent(out) :: nfound
    real, dimension(nmaxroot), intent(out) :: eelist
    integer, dimension(nmaxroot), intent(out) :: eetype

    real :: ki, tmp, prdmax, eelogmax
    integer :: nfound2, i1, ipos
    real, dimension(nmaxroot) :: eelist2
    real, dimension(3) :: abltg

    if ((ipphi .lt. 1) .and. (ipphi .gt. this%npphin)) then
       ! the pphi index is illegal
       nfound = -1
       return
    end if

    ! find in n grid first
    call spline_find(this%periodn(ipphi), period, 1, this%neen-1, nfound, eelist)
    do i1 = 1, nfound
       ! set type tag to "n"
       eetype(i1) = 0
    end do
    
    ! check if a TYPE I boundary exists for the pphi index
    if (istype1(this, ipphi)) then
       ipos = indexn2b(this, ipphi)
       call spline_find(this%periodb(ipos), period, 1, this%neeb, nfound2, eelist2)
       if (nfound2 .ge. 1) then
          ! merge the result
          do i1 = 1, nfound2
             ! set type tag to "b"
             eetype(i1 + nfound) = 1
             eelist(i1 + nfound) = eelist2(i1)
          end do
          nfound = nfound + nfound2
       end if
       prdmax = this%periodb(ipos)%y(this%neeb)
       if (period .gt. prdmax) then
          ! if the required period is still larger than the last b grid
          ! extrapolation needed
          nfound = nfound + 1
          ! set type tag to "b"
          eetype(nfound) = 1
          eelogmax =  this%periodb(ipos)%x(this%neeb)
          tmp = spline_interp(this%periodb(ipos), eelogmax, abltg)
          ki= abltg(1)
          eelist(nfound) = eelogmax +  (period-prdmax) / ki
     end if
  end if
  end subroutine findperiod

  subroutine getperiod(this, ee, ipphi, prd, dprd)
    ! get period by giving the pphi grid and the energy
    ! input: energy and pphi grid number
    ! return period, d period / d E
    implicit none

    type(tgrid), intent(in) :: this
    real, intent(in) :: ee
    integer, intent(in) :: ipphi
    real, intent(out) :: prd
    real, dimension(3) :: dprd

    integer :: ierr, itype1, ipos
    real :: eelog, ki, tmp
    
    ! check if the energy and pphi index are out of range
    ! see if the special treatment close to TYPE I t/p boundary is needed 
    
    itype1 = gettgridtype(this, ee, ipphi, ierr)

    ! if the required data point is out of range
    if (ierr .ge. 1) then
      !write(*,*) ipphi
       prd = -1.
       dprd(:) = -1.
       return
    end if

    if (itype1 .ge. 1) then
       ! between the last grid point and t/p boundary
       eelog = eetoeelog(this, ee, ipphi)
       call getperiodb1(this, eelog, ipphi, prd, dprd)
    else
       prd = spline_interp(this%periodn(ipphi), ee, dprd)
    end if
  end subroutine getperiod

  subroutine getperiodb(this, eelog, ipphi, prd, dprd)
    ! get period and its derivative when close to TYPE I t/p boundary
    ! with a check of energy range, will call getperiodb1
    ! energy input in log scale : E = Etpbound - mub0 exp(-eelog)
    ! ipphi : index of pphi grid
    ! output : prd, dprd
    implicit none

    type(tgrid), intent(in) :: this
    real, intent(in) :: eelog
    integer, intent(in) :: ipphi
    real, intent(out) :: prd
    real, dimension(3), intent(out) :: dprd
    
    integer :: ierr, itype1, i1
    real :: ee
    
    ee = eelogtoee(this, eelog, ipphi)
    itype1 = gettgridtype(this, ee, ipphi, ierr)
    if ((ierr .ge. 1) .and. (itype1 .le. 0)) then
       ! out of bound or no special treatment is needed
       prd = -1.
       do i1 = 1, 3
          dprd(i1) = -1.
       end do
    else
       call getperiodb1(this, eelog, ipphi, prd, dprd)
    end if
    
  end subroutine getperiodb
  
  real function getvpm(this, ee, ipphi, p, m)
    ! get the finite element weight matrix vpm
    ! INPUT : energy, pphi index, harmonic p, finite element index m
    implicit none
    
    type(tgrid), intent(in) :: this
    real, intent(in) :: ee
    integer, intent(in) :: ipphi, p, m
    
    real :: eelog
    integer :: itype1, ierr
    real, dimension(3) :: abltg

    getvpm = 0.

    ! check if the inputs are illegal, and to use n or b grid
    itype1 = gettgridtype(this, ee, ipphi, ierr)
    if ((p .le. 0) .and. (p .ge. this%np)) then
       ierr = 1
    end if
    if ((m .le. 0) .or. (m .gt. nele)) ierr = 1
    if (ierr .ge. 1) return

    ! save some computation time
    if (m < this%ielementmin(ipphi) .or. m > this%ielementmax(ipphi)) return

    if (itype1 .ge. 1) then
       eelog = eetoeelog(this, ee, ipphi)
       getvpm = getvpmb1(this, eelog, ipphi, p, m)
    else
       getvpm = spline_interp(this%vpmgridn(m,p,ipphi), &
            ee, abltg)
    end if

  end function getvpm

  real function getvpmb(this, eelog, ipphi, p, m)
    ! get the finite element weight matrix vpm
    ! INPUT : energy, pphi index e, harmonic p, finite element index m
    use radial_grid
    implicit none
    
    type(tgrid), intent(in) :: this
    real, intent(in) :: eelog
    integer, intent(in) :: ipphi, p, m

    real :: ee
    integer :: itype1, ierr

    getvpmb = 0.

    ! check if the inputs are illegal, and to use n or b grid
    ee = eelogtoee(this, eelog, ipphi)
    itype1 = gettgridtype(this, ee, ipphi, ierr)
    if ((p .le. 0) .and. (p .ge. this%np)) then
       ierr = 1
    end if
    if ((m .le. 0) .or. (m .gt. nele)) ierr = 1
    if ((ierr .ge. 1) .or. (itype1 .le. 0)) return

    getvpmb = getvpmb1(this, eelog, ipphi, p, m)

  end function getvpmb

  subroutine tgrid_init(this, mub0, npphin, neen, neeb, eeendb, np)
    ! initiate the object, create the grid
    ! INPUT : mub0, number of pphi grid, number of n energy grid,
    !         number of b energy grid, b grid upper limit,
    !         number of harmonics
    use paras_num, only : mindiff_el
    use radial_grid
    use orbit_classify
    implicit none
    
    type(tgrid) :: this
    real, intent(in):: mub0, eeendb
    integer, intent(in) :: npphin, neen, neeb, np
    
    real :: ds, traplosttmp, trapedgetmp, dee, eestartb
    integer :: i1, i2, i3, i4, i5, ipos, negrid
    integer :: istat
    
    this%mub0 = mub0
    this%np = np
    this%neen = neen
    this%neeb = neeb

    ! size of the grid, in unit of s flux label
    ds = 1. / real(npphin + 1)
    i1 = 1
    i2 = npphin
    trapedgetmp = stagedge(mub0, sgrid(ds * real(i1)), 1, istat)
    ! determine the range of pphi
    ! right boundary
    do while (istat .le. 0)
       i1 = i1 + 1
       trapedgetmp = stagedge(mub0, sgrid(ds * real(i1)), 1, istat)
       if (i1 .ge. npphin) exit
    end do
    
    ! left boundary, crossing of lost boundary and trap edge
    trapedgetmp = stagedge(mub0, sgrid(ds * real(i2)), 1, istat)
    do while ((trapedgetmp .ge. traplost(mub0, sgrid(ds * real(i2)))&
         - mub0*mindiff_el) .or. (istat .ne. 1))
       ! mindiff_el : the least difference of trapedge and traplost
       i2 = i2 - 1
       trapedgetmp = stagedge(mub0,  sgrid(ds * real(i2)), 1, istat)
       if (i2 .le. 1) exit
    end do

    this%npphin = (i2 - i1) + 1
    this%pphistart =  sgrid(ds * real(i1))
    this%pphiend = sgrid(ds * real(i2))
    this%ipphistart = i1
    this%ipphiend = i2
    ! if grid number is greater than 1
    if (this%npphin .ge. 1) then
       ! allocate grids
       
       allocate(this%pphigrid(this%npphin))
       allocate(this%etpbound(this%npphin))
       allocate(this%etrapedge(this%npphin))
       allocate(this%elost(this%npphin))
       this%ibstart = 0
       this%iloststart = i2
       do i3 = 1, this%npphin
          ipos = i1 + i3 - 1
          this%pphigrid(i3) =  sgrid(ds * real(ipos))

          ! calculate the t/p boundary
          this%etpbound(i3) = tpbound(mub0, this%pphigrid(i3), istat)

          if (istat .eq. 2) then
             ! TYPE II t/p boundary found
             ! no special treatment for the t/p boundary
             this%ibstart = ipos
          end if
          ! calculate the lowest energy of trapped particles

          this%etrapedge(i3) = stagedge(mub0, this%pphigrid(i3), 1, istat)
          if (istat .ne. 1) then
             write(*,*) 'error when calculating trap edge'
          end if

          ! calculate the lost boundary
          traplosttmp = traplost(mub0, this%pphigrid(i3))
          this%elost(i3) = traplosttmp

          if (traplosttmp .le. this%etpbound(i3)) then
             ! the lost boundary is lower than the t/p boundary
             if (this%iloststart .gt. ipos) then
                this%iloststart = ipos
             end if
          end if
       end do
       this%ibstart = this%ibstart + 1
       this%npphib = -this%ibstart + this%iloststart

       ! allocate period spline grid
       allocate(this%periodn(this%npphin))
       allocate(this%periodb(this%npphib))

       ! allocate n/b switch
       allocate(this%enbswitch(this%npphib))

       ! allocate orbit integral spline grid
       allocate(this%vpmgridn(nele, np, this%npphin))
       allocate(this%vpmgridb(nele, np, this%npphib))
       ! allocate the array stores most inside/outside element indexes
       allocate(this%ielementmin(this%npphin))
       allocate(this%ielementmax(this%npphin))
       allocate(this%ielementminn(this%neen, this%npphin))
       allocate(this%ielementmaxn(this%neen, this%npphin))
       allocate(this%ielementminb(this%neeb, this%npphib))
       allocate(this%ielementmaxb(this%neeb, this%npphib))

       this%np = np
       ! fill in energy grid
       do i3 = 1, this%npphin
          ! initiate the spline object
          ! energy grid size
          dee = (this%etpbound(i3) - this%etrapedge(i3)) / real(neen - 1)
          traplosttmp = this%elost(i3)
          if (traplosttmp .le. this%etpbound(i3)) then
             ! if trapped lost boundary is lower than t/p boundary
             ! grid from trap edge to trap lost boundary
             ! reduce the number of grid points
             negrid = floor((traplosttmp - this%etrapedge(i3))/dee)  + 1
             ! min 3 grid to allow spline
             if (negrid .lt. 3) negrid = 3
             ! adjust grid size to let two ends lay on the gridS
             dee = (traplosttmp - this%etrapedge(i3)) / float(negrid - 1)
             call spline_init(this%periodn(i3), negrid)
          else
             ! grid from trap edge to t/p boundary
             call spline_init(this%periodn(i3), neen)
             negrid = neen
          end if

          do i4 = 1, negrid
             ! equal-distant energy grid
             this%periodn(i3)%x(i4) = this%etrapedge(i3) + real(i4 - 1) * dee
          end do
       end do


       ! special grid for TYPE I t/p boundary
        ! fill in the energy grid
       if (this%npphib .gt. 0) then
          do i3 = 1, this%npphib
             ! initiate the spline object
             call spline_init(this%periodb(i3), neeb)
             ipos = indexb2n(this, i3)
             ! start from the second last grid point of the n grid
             eestartb = eetoeelog(this, this%periodn(ipos)%x(neen-1), ipos)
             ! fill in the n/b switch
             this%enbswitch(i3) = this%periodn(ipos)%x(neen-1)
             ! grid size in log scale
             dee = (eeendb - eestartb) / (neeb - 1)
             do i4 = 1, neeb
                ! equal-distant in x (ee = eebound (1-e^-x))
                this%periodb(i3)%x(i4) = eestartb + dee * real(i4 - 1)
             end do
             ! feed the second data point of b grid to n grid
             ! (to allow better interpolation)
             this%periodn(ipos)%x(this%neen) = eelogtoee(this, this%periodb(i3)%x(2), ipos)
          end do
       end if

      ! copy the same energy grid to grid of vpm
       do i1 = 1, this%npphin
          do i4 = 1, this%np
             do i3 = 1, nele
                call spline_init(this%vpmgridn(i3,i4,i1), this%periodn(i1)%n)
                do i2 = 1, this%periodn(i1)%n
                   this%vpmgridn(i3,i4,i1)%x(i2) = this%periodn(i1)%x(i2)
                end do
             end do
          end do
       end do

       do i1 = 1, this%npphib
          do i4 = 1, this%np
             do i3 = 1, nele
                call spline_init(this%vpmgridb(i3,i4,i1), neeb)
                do i2 = 1, this%periodb(i1)%n
                   this%vpmgridb(i3,i4,i1)%x(i2) = this%periodb(i1)%x(i2)
                end do
             end do
          end do
       end do

    else
       write(*,*) 'no trapped particle at (mub0) =', this%mub0
    end if
  end subroutine tgrid_init

  subroutine tgrid_destroy(this)
    ! clean all the allocatables
    use radial_grid
    implicit none
    type(tgrid) :: this

    integer :: i1, i2, i3

    if (allocated(this%vpmgridn)) then
       ! destroy the spline objects one by one
       do i1 = 1, nele
          do i2 = 1, this%np
             do i3 = 1, this%npphin
                call spline_destroy(this%vpmgridn(i1, i2, i3))
             end do
          end do
       end do
       deallocate(this%vpmgridn)
    end if

    if (allocated(this%vpmgridb)) then
       ! destroy the spline objects one by one
       do i1 = 1, nele
          do i2 = 1, this%np
             do i3 = 1, this%npphib
                call spline_destroy(this%vpmgridb(i1, i2, i3))
             end do
          end do
       end do
       deallocate(this%vpmgridb)
    end if

    if (allocated(this%periodn)) then 
       ! destroy the spline objects
       do i1 = 1, this%npphin
          call spline_destroy(this%periodn(i1))
       end do
       deallocate(this%periodn)
    end if

    if (allocated(this%periodb)) then 
       ! destroy the spline objects
       do i1 = 1, this%npphib
          call spline_destroy(this%periodb(i1))
       end do
       deallocate(this%periodb)
    end if

    if (allocated(this%ielementmin )) deallocate(this%ielementmin )
    if (allocated(this%ielementmax )) deallocate(this%ielementmax )
    if (allocated(this%ielementminn)) deallocate(this%ielementminn)
    if (allocated(this%ielementmaxn)) deallocate(this%ielementmaxn)
    if (allocated(this%ielementminb)) deallocate(this%ielementminb)
    if (allocated(this%ielementmaxb)) deallocate(this%ielementmaxb)
    if (allocated(this%enbswitch)) deallocate(this%enbswitch)
    if (allocated(this%pphigrid)) deallocate(this%pphigrid)
    if (allocated(this%etpbound)) deallocate(this%etpbound)
    if (allocated(this%etrapedge)) deallocate(this%etrapedge)

  end subroutine tgrid_destroy

  subroutine tgrid_calculate(this)
    ! fill in the data
    use paras_phy, only : eunit, ei
    use profile, only : psi1
    use paras_num
    use radial_grid
    use orbit
    use orbit_integral
    use orbit_classify
    implicit none

    type(tgrid) :: this

    integer :: i1, i2, i3, i4
    integer :: istat, ipos, imin, imax, istart, iend
    real :: lostbound, ee, dtorbit
    real, dimension(norbitintsample) :: r, theta
    real, allocatable, dimension(:, :) :: work

    
    ! allocate the work space
    allocate(work(nele, this%np))
    ! calculate orbit period and orbit integral on all normal grid points
    do i1 = 1, this%npphin
       do i2 = 1, this%periodn(i1)%n
          call getorbit(this%periodn(i1)%x(i2), this%mub0,this%pphigrid(i1),&
               norbitintsample, 1, dtorbitn, r, theta, &
               this%periodn(i1)%y(i2), istat)
          if (istat .ne. 1) then
             this%periodn(i1)%y(i2) = -1.
             write(*,*) istat
             write(*,*) 'orbit does not exist on n grid at (muB0, E, pphi/e psi1) = '&
                  ,this%mub0/1000./eunit, this%periodn(i1)%x(i2)/1000./eunit, &
                  this%pphigrid(i1) / ei/psi1
             write(*,*) 'orbit status', istat
          else 
             call orbit_int(r, theta, norbitintsample,&
                  this%np, work, imin, imax)
             
             this%ielementminn(i2, i1) = imin
             this%ielementmaxn(i2, i1) = imax
             ! copy to spline objects
             do i3 = 1, nele
                do i4 = 1, this%np
                   this%vpmgridn(i3,i4,i1)%y(i2) = work(i3, i4)
                end do
             end do
          end if
       end do
    end do


    ! calculate orbit period and integral on t/p boundary particles
    do i1 = 1, this%npphib
       ipos = indexb2n(this, i1)
       ! first two energy points : copy from the calculated value
       this%periodb(i1)%y(1) = this%periodn(ipos)%y(this%neen-1)
       this%periodb(i1)%y(2) = this%periodn(ipos)%y(this%neen)
       this%ielementminb(1, i1) = this%ielementminn(this%neen-1, ipos)
       this%ielementmaxb(1, i1) = this%ielementmaxn(this%neen-1, ipos)
       this%ielementminb(2, i1) = this%ielementminn(this%neen  , ipos)
       this%ielementmaxb(2, i1) = this%ielementmaxn(this%neen  , ipos)

       do i3 = 1, nele
          do i4 = 1, this%np
             this%vpmgridb(i3,i4,i1)%y(1) = this%vpmgridn(i3,i4,ipos)%y(this%neen-1)
             this%vpmgridb(i3,i4,i1)%y(2) = this%vpmgridn(i3,i4,ipos)%y(this%neen)
          end do
       end do

       do i2 = 3, this%neeb
          ee = eelogtoee(this, this%periodb(i1)%x(i2), ipos)
          call getorbit(ee, this%mub0, this%pphigrid(ipos), norbitintsample, 1, dtorbitb, r, theta, this%periodb(i1)%y(i2), istat)
          if (istat .ne. 1) then
             this%periodb(i1)%y(i2) = -1.
             write(*,*) 'orbit does not exist on b grid at (muB0,E,ipphi) = ',&
this%mub0/eunit, ee/eunit, ipos
             write(*,*) 'orbit status', istat
          else
             call orbit_int(r, theta, norbitintsample, this%np, work, imin, imax)
             ! fresh the most inside/outside indexes
             this%ielementminb(i2, i1) = imin
             this%ielementmaxb(i2, i1) = imax
             do i3 = 1, nele
                do i4 = 1, this%np
                   this%vpmgridb(i3,i4,i1)%y(i2) = work(i3,i4)
                end do
             end do
          end if
       end do
    end do
   
    ! calculate all the splines
    do i1 = 1, this%npphin
       call spline_build(this%periodn(i1), 0., 0., 2, 1, this%periodn(i1)%n)
       ! equal-distant grid
       this%periodn(i1)%iequaldistant = .true.
       do i3 = 1, nele
          istart = this%periodn(i1)%n
          iend = this%periodn(i1)%n
          call findstartendn(this, i3, i1, istart, iend)
          do i4 = 1, this%np
             call spline_build(this%vpmgridn(i3,i4,i1), 0., 0., 2, istart, iend)
             this%vpmgridn(i3,i4,i1)%iequaldistant = .true.
          end do
       end do
       this%ielementmax(i1)=getmax(this%periodn(i1)%n, this%ielementmaxn(1, i1))
       this%ielementmin(i1)=getmin(this%periodn(i1)%n, this%ielementminn(1, i1))
    end do
    
    do i1 = 1, this%npphib
       ipos = indexb2n(this, i1)
       call spline_build(this%periodb(i1), 0., 0., 2, 1, this%neeb)
       this%periodb(i1)%iequaldistant = .true.
       do i3 = 1, nele
          call findstartendb(this, i3, i1, istart, iend)
          do i4 = 1, this%np
             call spline_build(this%vpmgridb(i3,i4,i1), 0., 0., 2, istart, iend)
             this%vpmgridb(i3,i4,i1)%iequaldistant = .true.
          end do
       end do
       iend   = getmax(this%neeb, this%ielementmaxb(1, i1))
       istart = getmin(this%neeb, this%ielementminb(1, i1))
       if (iend   .gt. this%ielementmax(ipos)) this%ielementmax(ipos) = iend
       if (istart .lt. this%ielementmin(ipos)) this%ielementmin(ipos) = istart
    end do

    if (allocated(work)) deallocate(work)
  end subroutine tgrid_calculate


  logical function istype1(this, ipphi)
    ! check if a TYPE I t/p boundary exists for the given pphi index
    implicit none
    type(tgrid), intent(in) :: this
    integer, intent(in) :: ipphi
    
    istype1 = (ipphi .ge. this%ibstart) .and. (ipphi .lt. this%iloststart)
  end function istype1

!********* INTERNAL FUNCTIONS ************

  subroutine getperiodb1(this, eelog, ipphi, prd, dprd)
    ! get period and its derivative when close to TYPE II t/p boundary
    ! without a check of energy range
    ! energy input in log scale : E = Etpbound - mub0 exp(-eelog)
    ! ipphi : index of pphi grid
    ! output : prd, dprd
    implicit none

    type(tgrid), intent(in) :: this
    real, intent(in) :: eelog
    integer, intent(in) :: ipphi
    real, intent(out) :: prd
    real, dimension(3) :: dprd

    integer :: ipos
    real :: tmp, ki
    ipos = indexn2b(this, ipphi)
         
    if (eelog .ge. this%periodb(ipos)%x(this%neeb)) then
       ! closer to t/p boundary than the last data point
       ! need extrapolation
       tmp = spline_interp(this%periodb(ipos), this%periodb(ipos)%x(this%neeb), dprd)
       ki = dprd(1)
       prd = ki * (eelog - this%periodb(ipos)%x(this%neeb)) + this%periodb(ipos)%y(this%neeb)
       dprd(2) = 0.
       dprd(3) = 0.
    else
       prd = spline_interp(this%periodb(ipos), eelog, dprd)
    end if

  end subroutine getperiodb1

  real function getvpmb1(this, eelog, ipphi, p, m)
    ! get the finite element weight matrix vpm
    ! without checking if the inputs are illegal
    ! INPUT : energy log, pphi index, harmonic p, finite element index m
    ! OUTPUT: vpm
    implicit none
    
    type(tgrid), intent(in) :: this
    real, intent(in) :: eelog
    integer, intent(in) :: ipphi, p, m
    
    integer :: ipos
    real, dimension(3) :: abltg
    real :: eeloglast, prdlast, prd
    real, dimension(3) :: dprd

    ipos = indexn2b(this, ipphi)
    eeloglast = this%periodb(ipos)%x(this%neeb)
    if (eelog .gt. eeloglast) then
       ! extrapolation is needed
       prdlast = this%periodb(ipos)%y(this%neeb)
       call getperiodb1(this, eelog, ipphi, prd, dprd)
       getvpmb1 = spline_interp(this%vpmgridb(m, p,ipos), &
            eeloglast, abltg)
       getvpmb1 = getvpmb1 * prdlast / prd
    else
       getvpmb1 = spline_interp(this%vpmgridb(m, p, ipos), &
            eelog, abltg)
    end if

  end function getvpmb1


  integer function gettgridtype(this, ee, ipphi, ierr)
    ! check if special treatment near t/p bound is needed (return 1)
    ! check if the combination of energy and pphi index is out of range
    ! ierr = 1 : out of range
    implicit none

    type(tgrid), intent(in) :: this
    real, intent(in) :: ee
    integer, intent(in) :: ipphi
    integer, intent(out) :: ierr
    
    ierr = 0
    gettgridtype = 0
    if ((ipphi .le. 0) .or. (ipphi .gt. this%npphin)) then
       ! ipphi index out of bound
       ierr = 1
       return
    end if
    if (ee .lt. this%periodn(ipphi)%x(1)) then
       ! energy lower than trap edge, return error
       !write(*,*) 'type 1 error', this%periodn(ipphi)%x(1), ee
       ierr = 1
       return
    end if
    if (istype1(this, ipphi)) then
       ! if a TYPE I t/p boundary is found for this pphi index
       if (ee .ge. this%etpbound(ipphi)) then
          ! energy higher than t/p boundary, return error
          !write(*,*) 'type 2 error',this%etpbound(ipphi), ee
          ierr = 1
          return
       end if
       if (ee .ge. this%periodn(ipphi)%x(this%periodn(ipphi)%n-1)) then
          ! sufficiently close to the t/p boundary
          gettgridtype = 1
       end if
    else
       if (ee .ge. this%periodn(ipphi)%x(this%periodn(ipphi)%n)) then
          ! energy higher than the upper limit, return error
          !write(*,*) 'type 3 error',this%etpbound(ipphi), ee
          ierr = 1
          return
       end if
    end if
  end function gettgridtype

  real function eelogtoee(this, eelog, ipphi)
    ! ee = etpbound - mub0 * exp(-eelog)
    ! without a check of data range out of bound or not
    implicit none
    type(tgrid), intent(in) :: this
    real, intent(in) :: eelog
    integer, intent(in) :: ipphi

    eelogtoee = this%etpbound(ipphi) - this%mub0 * exp(-eelog)

  end function eelogtoee

  real function eetoeelog(this, ee, ipphi)
    ! ee = etpbound - mub0 * exp(-eelog)
    ! without a check of data range out of bound or not
    implicit none
    type(tgrid), intent(in) :: this
    real, intent(in) :: ee
    integer, intent(in) :: ipphi

    eetoeelog = - log((this%etpbound(ipphi) - ee) / this%mub0)

  end function eetoeelog

  integer function indexn2b(this, ipphin)
    ! convert n pphi index to b pphi index
    implicit none
    type(tgrid), intent(in) :: this
    integer, intent(in) :: ipphin
    
    indexn2b = ipphin - this%ibstart + 1
  end function indexn2b

  integer function indexb2n(this, ipphib)
    ! convert b pphi index to n pphi index
    implicit none
    type(tgrid), intent(in) :: this
    integer, intent(in) :: ipphib
    
    indexb2n = ipphib + this%ibstart - 1
  end function indexb2n

  real function sgrid(s)
    ! define the pphi grid
    use paras_phy, only : ei
    use profile, only : psi1
    implicit none
    real, intent(in) :: s

    sgrid = -ei * psi1 * s
  end function sgrid
  
  subroutine findstartendn(this, ielement, ipphi, istart, iend)
    ! find the start and end energy grid index for element i
    

    type(tgrid), intent(in) :: this
    integer, intent(in) :: ielement, ipphi
    integer, intent(out):: istart, iend

    integer :: i1
    do i1 = 1, this%periodn(ipphi)%n
       if ((this%ielementminn(i1,ipphi) .le. ielement) .and. &
            (this%ielementmaxn(i1,ipphi) .ge. ielement)) then
          istart = i1
          exit
       end if
    end do

    do i1 = this%periodn(ipphi)%n, 1, -1
       if ((this%ielementminn(i1,ipphi) .le. ielement) .and. &
            (this%ielementmaxn(i1,ipphi) .ge. ielement)) then
          iend = i1
          exit
       end if
    end do
  end subroutine findstartendn

  subroutine findstartendb(this, ielement, ipphi, istart, iend)
    ! find the start and end energy grid index for element i
    

    type(tgrid), intent(in) :: this
    integer, intent(in) :: ielement, ipphi
    integer, intent(out):: istart, iend

    integer :: i1
    
    do i1 = 1, this%neeb
       istart = i1
       if ((this%ielementminb(i1,ipphi) .le. ielement) .and. &
            (this%ielementmaxb(i1,ipphi) .ge. ielement)) then
          exit
       end if
    end do

    do i1 = this%neeb, 1, -1
       iend = i1
       if ((this%ielementminb(i1,ipphi) .le. ielement) .and. &
            (this%ielementmaxb(i1,ipphi) .ge. ielement)) then
          exit
       end if
    end do
  end subroutine findstartendb

  integer function getmin(n, x)
    ! get the min  of a given array x
    implicit none
    
    integer, intent(in) :: n
    integer x(*)

    integer :: i1
    
    getmin = x(1) 
    do i1 = 2, n
       if (x(i1) .lt. getmin) getmin = x(i1)
    end do

  end function getmin

  integer function getmax(n, x)
    ! get the max of a given array x
    implicit none
    
    integer, intent(in) :: n
    integer x(*)

    integer :: i1
    
    getmax = x(1) 
    do i1 = 2, n
       if (x(i1) .gt. getmax) getmax = x(i1)
    end do

  end function getmax
  
end module trap_grid
     
