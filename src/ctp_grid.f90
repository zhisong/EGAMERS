  ! a grid for counter-passing particles for a single mub0
  
module ctp_grid
  ! "n" grids are Normal equal-distant energy grids 
  !     from trap edge to t/p boundary or lost boundary
  ! "b" grid are special grids for TYPE I t/p Boundary
  !     exponentially close to the t/p boundary
  !
  ! TYPE:
  ! ctpgrid           -  contains all the grids and splines
  !
  ! SUBROUTINES:
  ! (public)
  ! ctpgrid_init      - allocate arrays and splines, fill energy grid
  !                     must call first

  use radial_grid, only : nele
  use spline_module
  implicit none

  ! type of ctp particle grid
  ! (which is exactly the same as trapped grid "tgrid")
  type, public :: ctpgrid

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
     real, dimension(:), allocatable, public :: estagedge, etpbound, elost

     ! mub0
     real, public :: mub0

     ! harmonic count
     integer, public :: np

     ! enabled t/p boundary
     integer, public :: itpbound

  end type ctpgrid

  private
  public :: ctpgrid_init !, ctpgrid_destroy, ctpgrid_calculate, ctpsgrid

contains

  subroutine ctpgrid_init(this, mub0, eehigh, eelow, npphin, neen, neeb, eeendb, np, itpbound)
    ! initiate the object, create the grid
    ! INPUT : mub0, the upper limit of energy, the lower limit of energy
    !         number of pphi grid, number of n energy grid,
    !         number of b energy grid, b grid upper limit,
    !         number of harmonics, treat the t/p boundary(1) or not(0)
    use paras_phy, only : ei, eunit
    use paras_num, only : mindiff_el
    use profile, only : psi1
    use radial_grid
    use orbit_classify
    implicit none

    type(ctpgrid) :: this
    real, intent(in) :: mub0, eehigh, eelow, eeendb
    integer, intent(in) :: npphin, neen, neeb, np, itpbound

    real :: pphimin, pphimax, ds, tpboundtmp, stagedgetmp, losttmp
    real :: dee, eemin, eemax
    integer :: istat, i1, i2, i3, i4, ipos
    
    this%mub0 = mub0
    this%np = np
    this%neen = neen
    this%neeb = neeb

    ! t/p boundary turn on?
    if (itpbound .le. 0) then
       if (eelow .le. ctplost(mub0, -ei * psi1)) then
          ! specified energy lower then tpbound, tpbound forced to turn on
          this%itpbound = 1
          write(*,*) 'In ctpgrid_init : lowest grid energy from the namelist is lower than the t/p boundary, t/p boundary function forced to enable'
       else
          this%itpbound = 1
       end if
    else
       this%itpbound = 1
    end if

    ! the min |pphi|, 0 or stag of the lowest energy particle
    if (this%itpbound .eq. 0.) then
       pphimin = 0.
    else
       pphimin = stagedgepphi(mub0, eelow, -1, istat)
       if (istat .ne. 1) then
          write(*,*) 'Finding stag edge failed in ctpgrid_init, mub0 = ', &
               mub0 / eunit / 1000., ' keV'
          pphimin = 0.
       end if
    end if

    ! the max |pphi|, lost boundary of the highest energy particle
    pphimax = ctplostpphi(mub0, eehigh)

    ! Step size in pphi, ignore the two end points 
    ds = 1. / real(npphin - 1)
    i1 = 1
    i2 = npphin
    ! right boundary
    ! to find the start pphi index if t/p boundary is considered
    if (this%itpbound .ge. 1) then
       tpboundtmp = tpbound(mub0, sgrid(ds * real(i1), pphimin, pphimax), istat)
       do while (istat .ne. 1)
          ! TYPE II t/p boundary means no ctp particle for given pphi
          ! move outward until a TYPE I t/p boundary is reached
          i1 = i1 + 1
          tpboundtmp = tpbound(mub0, sgrid(ds * real(i1), pphimin, pphimax),&
               istat)
          if (i1 .ge. npphin) exit
       end do
    end if

    this%npphin = (i2 - i1) + 1
    this%pphistart = sgrid(ds * real(i1), pphimin, pphimax)
    this%pphiend = sgrid(ds * real(i1), pphimin, pphimax)
    this%ipphistart = i1
    this%ipphiend = i2
    ! if grid number is greater than 1
    if (this%npphin .ge. 1) then
       ! allocate grids
       allocate(this%pphigrid(this%npphin))
       allocate(this%estagedge(this%npphin))
       allocate(this%elost(this%npphin))
       if (this%itpbound .ge. 1) then
          ! if t/p bound is turned on, allocate the array to store t/p bound
          allocate(this%etpbound(this%npphin))
          ! the start pphi index of t/p bound
          this%ibstart = 1
       else
          ! t/p bound turned off, move the index out of the range
          this%ibstart = i2 + 1
       end if
       this%iloststart = i2
       
       do i3 = 1, this%npphin
          ipos = i1 + i3 - 1
          ! fill in the value of pphi to the grid
          this%pphigrid(i3) = sgrid(ds * real(ipos), pphimin, pphimax)

          ! the stag edge
          this%estagedge(i3) = stagedge(mub0, this%pphigrid(i3), 1 , istat)
          if (istat .ne. 1) write(*,*) 'err in calculating stag edge'
          
          ! if pphi < -ei * psi1
          if (this%pphigrid(i3) .le. -ei * psi1) then
             ! calcuate the lost boundary
             this%elost(i3) = ctplost(mub0, this%pphigrid(i3))
             if (i3 .lt. this%iloststart) then
                ! determine the pphi lost index
                this%iloststart = i3
             end if
          else
             ! calculate the t/p boundary if used
             if (this%itpbound .ge. 1) then
                this%etpbound(i3) = tpbound(mub0, this%pphigrid(i3), istat)
                if (istat .ne. 1) then
                   ! find tpbound failed, return error
                   write(*,*) 'error in finding tpbound in ctpgrid_init'
                   write(*,*) 'pphi/ei psi1 = ', this%pphigrid(i3) / ei / psi1
                   write(*,*) 'mub0         = ', mub0 / eunit/ 1000., ' keV'
                end if
             end if
          end if
       end do
       if (this%itpbound .ge. 1) then
          ! number of b grid
          this%npphib = -this%ibstart + this%iloststart
       end if

       ! allocate period spline grid
       allocate(this%periodn(this%npphin))
       if (this%itpbound .ge. 1) allocate(this%periodb(this%npphib))

       ! allocate n/b switch
       if (this%itpbound .ge. 1) allocate(this%enbswitch(this%npphib))

       ! allocate orbit integral spline grid
       allocate(this%vpmgridn(nele, np, this%npphin))
       if (this%itpbound .ge. 1) allocate(this%vpmgridb(nele, np, this%npphib))
       ! allocate the array stores most inside/outside element indexes
       allocate(this%ielementmax(this%npphin))
       allocate(this%ielementmin(this%npphin))
       allocate(this%ielementmaxn(this%neen, this%npphin))
       allocate(this%ielementminn(this%neen, this%npphin))
       if (this%itpbound .ge. 1) allocate(this%ielementmaxb(this%neeb, this%npphin))
       if (this%itpbound .ge. 1) allocate(this%ielementminb(this%neeb, this%npphin))

       this%np = np
       ! fill in energy grid
       do i3 = 1, this%npphin
          ! initiate the spline object
          call spline_init(this%periodn(i3), neen)
          ! determine the lowest energy of this pphi grid point
          if (this%itpbound .ge. 1) then
             ! need special treatment for t/p bound
             if (i3 .ge. this%iloststart) then
                ! the min energy is the lost bound energy
                eemin = ctplost(mub0, this%pphigrid(i3))
             else
                ! the min energy is the t/p bound energy
                eemin = tpbound(mub0, this%pphigrid(i3), istat)
                this%etpbound(i3) = eemin
                if (istat .ne. 1) then
                   ! find tpbound failed, return error
                   write(*,*) 'error in finding tpbound in ctpgrid_init'
                   write(*,*) 'pphi/ei psi1 = ', this%pphigrid(i3) / ei / psi1
                   write(*,*) 'mub0         = ', mub0 / eunit/ 1000., ' keV'
                end if
             end if
          else
             eemin = eelow
             ! move the min enenrgy to the lost bound if eemin is lower
             if (this%elost(i3) .ge. eemin) eemin = this%elost(i3)
          end if

          ! determine the highest energy of the pphi grid point
          eemax = eehigh
          if (this%estagedge(i3) .le. eemax) eemax = this%estagedge(i3)

          dee = (eemax - eemin) / float(neen - 1)

          do i4 = 1, neen
             ! equal-distant energy grid
             this%periodn(i3)%x(i4) = eemin + float(i4 - 1) * dee
          end do
       end do

       ! special grid for t/p boundary if tpbound is enabled
       
    end if
  end subroutine ctpgrid_init

  real function sgrid(s, pphimin, pphimax)
    ! define the ctp pphi grid
    implicit none

    real, intent(in) :: s, pphimin, pphimax

    sgrid = s * (pphimax - pphimin) + pphimin

  end function sgrid
  
end module ctp_grid


