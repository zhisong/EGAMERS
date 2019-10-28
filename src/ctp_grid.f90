! Creating grids and calculating orbit data for counter passing particles
! Mostly copies trapped module, main differences are:
! - Only equidistant grids are used, no exponential grids
! - Calculation of grid boundaries: bounded above and below by pitch angle 
!   parameters (and also by NBI beam energy)

module ctp_grid
    use spline_module
    use radial_grid, only : nele
    implicit none
  
    type, public :: ctpgrid
  
        ! finite element weight for grid
        ! indexes : element label, pth resonance, pphi grid
  
        type(spline), dimension(:,:,:), allocatable, public :: vpmgridn
  
        ! orbit period grid (spline)
        type(spline), dimension(:), allocatable, public :: periodn

        ! store value of df0/dE to show weight of each point in grid when
        ! using grid to calculate EGAM data
        real, dimension(:,:), allocatable, public :: ctpweight
        real, dimension(:,:), allocatable, public :: ctpgam


  
        ! number of pphi grid points
        integer, public :: npphin
  
        ! pphi grid
        real, dimension(:), allocatable, public :: pphigrid, s
        ! pphi grid start/end points
        real, public :: pphi_start, pphi_end ! Lowest/highest pphi values
        integer, public :: i_pphi_start, i_pphi_end ! Index (centered around -1)
  
        ! Number of energy grid points
        integer, public :: neen
  
        ! Store grid data for ctp lost & on axis boundaries
        real, dimension(:), allocatable, public :: ctp_lost_E, on_axis_E
  
        ! Store lowest value of ctp energy
        real, dimension(:), allocatable, public :: ectp_min 
  
        ! mub0
        real, public :: mub0
  
        ! harmonic count
        integer, public :: np
  
        ! Step size
        real, public :: ds
  
    end type ctpgrid
  
    public :: ctpgrid_init, ctpgrid_destroy, ctpgrid_calculate
  
  contains
  
    subroutine ctpgrid_init(this, mub0, lambda0, dlambda, npphin, neen, np)
        ! initiate the object, create the grid
        ! INPUT: mub0, pitch angle, pitch angle width, number of pphi grid, number of n energy grid, number of harmonics
        use orbit_classify
        use radial_grid
        use paras_phy, only : eunit
        implicit none
    
        type(ctpgrid) :: this
        real, intent(in) :: mub0, lambda0, dlambda
        integer, intent(in) :: npphin, neen, np
  
        real :: ds, ctp_temp, ctp_lost_temp, on_axis_temp, min_ctp_E, max_ctp_E, dee, pphi_pos, offset, ctplost_min, ctplost_max
        integer :: i1, i2, i3, i4, negrid
  
        this%mub0 = mub0
        this%np = np
        this%neen = neen
  
        ! size of grid    
        ds = 1. / real(npphin + 1)  ! Size of grid, units of s flux label
        this%ds = ds
        i1 = 0 ! Counters for finding range of pphi
        i2 = 0

        ctplost_min = ctplost(mub0,sgrid(1.)) ! Bottom of ctplost parabola, minimum energy can't go below this
        ctplost_max = eeonaxis(mub0, sgrid(1.)) ! Energy where eeonaxis parabola crosses the pphi/ei/psi1 = -1 line
        ! Min/max energies to consider, depend on pitch angle. Within 2*delta

        min_ctp_E = mub0 / (lambda0 + 2 * dlambda)
        ! If min_ctp_E goes below allowed minimum, set it to be slightly above threshold
        if (min_ctp_E .le. ctplost_min) then 
            min_ctp_E = ctplost_min * 1.01 
        end if

        ! If denominator for upper bound goes below 0, or max energy goes above
        ! upper bound found earlier, just take 0.99 * upper bound
        if ((lambda0 - 2 * dlambda) .le. 0) then
            max_ctp_E = ctplost_max*0.99 
        else if ((mub0 / (lambda0 - 2 * dlambda)) .ge. ctplost_max) then
            max_ctp_E = ctplost_max*0.99 
        else
            max_ctp_E = mub0 / (lambda0 - 2 * dlambda)
        end if
  
        ! Find minimum pphi (left value) - where ctplosti E intersects max E
        ctp_temp = ctplost(mub0, sgrid(1.)) ! Start at minimum of ctplost parabola
  
        do while (ctp_temp .lt. max_ctp_E) ! Keep searching until ctplost intersects maxE
            i1 = i1 + 1
            ctp_temp = ctplost(mub0, sgrid(1. + ds * real(i1)))
        end do
  
        ! Find max pphi (right value) - where onaxis E intersects min E
       ctp_temp = eeonaxis(mub0, sgrid(1.)) 
       do while(ctp_temp .gt. min_ctp_E)
          i2 = i2 + 1 
          ctp_temp = eeonaxis(mub0, sgrid(1. - ds * real(i2)))
      end do
  
      ! pphi range, don't include final index, it is the first "out of bounds" 
      this%i_pphi_start = i1 - 1! Left pphi index = # ds steps to the left of -1 
      this%i_pphi_end = i2 - 1 ! Right pphi index = # ds steps to the right of -1
  
      this%pphi_start = sgrid(1. + ds * real(i1 - 1)) ! Left pphi value
      this%pphi_end = sgrid(1. - ds * real(i2 - 1)) ! Right pphi value
  
      this%npphin = (i2-1) + (i1-1) + 1 ! Number of pphi grid points
  
      ! if grid number > 1
      if (this%npphin .ge. 1) then
          ! Allocate grids
          allocate(this%pphigrid(this%npphin)) ! pphi grid
          allocate(this%ctp_lost_E(this%npphin)) ! store E of ctp lost boundary
          allocate(this%on_axis_E(this%npphin)) ! store E of on-axis boundary
          allocate(this%s(this%npphin)) ! Store pphi/ei/psi1
  
          do i3 = 1, this%npphin
              ! Get pphi/ei/psi1 val, starting from 1 + ds*(i1-1), until 1 - ds*(i2-1)
              pphi_pos = 1. - ds * real((i3-1) - (i1-1))
              ! Calculate pphi 
              this%pphigrid(i3) = sgrid(pphi_pos)
              ! Store pphi/ei/psi1
              this%s(i3) = pphi_pos
  
              ! Calculate ctp lost boundary
              this%ctp_lost_E(i3) = ctplost(mub0, this%pphigrid(i3))
  
              ! Calculate on axis boundary
              this%on_axis_E(i3) = eeonaxis(mub0, this%pphigrid(i3)) 
          end do
  
          ! allocate min ctp energy storing array
          allocate(this%ectp_min(this%npphin))
          ! allocate period spline grid
          allocate(this%periodn(this%npphin))
          ! allocate orbit integral spline grid
          allocate(this%vpmgridn(nele, np, this%npphin))
          ! allocate grid point df0/dE weights grid
          allocate(this%ctpweight(this%neen, this%npphin))
          allocate(this%ctpgam(this%neen, this%npphin))

  
          ! fill in energy grid
          do i3 = 1, this%npphin
              ctp_lost_temp = this%ctp_lost_E(i3)
              on_axis_temp = this%on_axis_E(i3)
              pphi_pos = -1. + ds * real((i3-1) - (i1-1))
              ! energy grid step size
              dee = (max_ctp_E - min_ctp_E) / real(neen - 1)
              ! this%ectp_min(i3) = min_ctp_E
  
              ! If for pphi/ei/psi1 < -1, the ctp lost energy is greater than the
              ! minimum energy, make grid from ctp_lost_E to max E & reduce
              ! number of grid points
  
              if ((pphi_pos .lt. -1.) .AND. (ctp_lost_temp .gt. min_ctp_E)) then
                 negrid = floor((max_ctp_E - ctp_lost_temp)/dee) + 1
  
                 if (negrid .lt. 5) negrid = 5
                 dee = (max_ctp_E - ctp_lost_temp) / float(negrid - 1) 
                 call spline_init(this%periodn(i3), negrid)
  
                 ! Lowest ctp energy is ctp_lost, not min_ctp_E
                 this%ectp_min(i3) = ctp_lost_temp
  
              else if (on_axis_temp .lt. max_ctp_E) then
                 negrid = floor((on_axis_temp - min_ctp_E)/dee) + 1
                 this%ectp_min(i3) = min_ctp_E
  
                 if (negrid .lt. 5) negrid = 5
                 offset = (on_axis_temp - min_ctp_E) * 0.01
                 dee = (on_axis_temp - min_ctp_E) / float(negrid - 1) - offset / float(negrid - 1) ! Offset a bit so it doesn't lie exactly on-axis
                 call spline_init(this%periodn(i3), negrid)
  
              else
                 call spline_init(this%periodn(i3), neen)
                 this%ectp_min(i3) = min_ctp_E
                 negrid = neen
              end if
  
              do i4 = 1, negrid
                 ! equidistant energy grid
                 this%periodn(i3)%x(i4) = this%ectp_min(i3) + real(i4 - 1)*dee
              end do
          end do
  
          ! copy grid to vpm grid
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
  
      else
          write(*,*) 'no counter passing particle at (mub0) =', this%mub0
      end if
  
    end subroutine ctpgrid_init
  
    subroutine ctpgrid_destroy(this)
      ! clean all the allocatables
      use radial_grid
      implicit none
      type(ctpgrid) :: this
  
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
  
      ! Deallocate periodn splines
      if (allocated(this%periodn)) then 
          ! destroy the spline objects
          do i1 = 1, this%npphin
             call spline_destroy(this%periodn(i1))
          end do
          deallocate(this%periodn)
       end if
  
      if (allocated(this%pphigrid)) deallocate(this%pphigrid)
      if (allocated(this%ctp_lost_E)) deallocate(this%ctp_lost_E)
      if (allocated(this%on_axis_E)) deallocate(this%on_axis_E)
      if (allocated(this%ectp_min)) deallocate(this%ectp_min)
      if (allocated(this%ctpweight)) deallocate(this%ctpweight)
      if (allocated(this%ctpgam)) deallocate(this%ctpgam)      
      if (allocated(this%vpmgridn)) deallocate(this%vpmgridn)
      if (allocated(this%s)) deallocate(this%s)
  
    end subroutine ctpgrid_destroy
    
    subroutine ctpgrid_calculate(this)
      ! fill in the data
      use distribution_fun, only : df0de
      use paras_phy, only : eunit, ei, pi
      use profile, only : psi1, omega_gam
      use paras_num
      use radial_grid
      use orbit
      use orbit_integral
      use orbit_classify
      implicit none
  
      type(ctpgrid) :: this
  
      integer :: i1, i2, i3, i4
      integer :: istat, imin, imax, istart, iend
      real, dimension(norbitintsample) :: r, theta
      real, allocatable, dimension(:, :) :: work
  
      ! allocate work space
      allocate(work(nele, this%np))
      ! calculate orbit period
      do i1 = 1, this%npphin
          do i2 = 1, this%periodn(i1)%n
              call getorbit(this%periodn(i1)%x(i2), this%mub0, this%pphigrid(i1),&
              norbitintsample, -1, dtorbitn, r, theta,&
              this%periodn(i1)%y(i2), istat)
  
              if (istat .ne. 1) then ! if orbit not successful
                  this%periodn(i1)%y(i2) = -1.
                  write(*,*) istat
                  write(*,*) 'orbit does not exist on n grid at (muB0, E, pphi/e psi1) = '&
                       ,this%mub0/1000./eunit, this%periodn(i1)%x(i2)/1000./eunit, &
                       this%pphigrid(i1) / ei/psi1
                  write(*,*) 'orbit status', istat
              ! orbit integral, don't worry about ielement
              else
                  call orbit_int(r, theta, norbitintsample, &
                  this%np, work, imin, imax)

                  ! get df0/dE weight
                  this%ctpweight(i2,i1) = df0de(this%periodn(i1)%x(i2), this%mub0, this%pphigrid(i1))
                  this%ctpgam(i2,i1) = omega_gam(findrstart(this%periodn(i1)%x(i2), this%mub0, this%pphigrid(i1), -1))

                  ! copy to spline objects
                  do i3 = 1, nele
                     do i4 = 1, this%np
                        this%vpmgridn(i3,i4,i1)%y(i2) = work(i3, i4)
                     end do
                  end do                  
              end if
          end do
      end do
  
      ! Calculate the splines
      do i1 = 1, this%npphin
          call spline_build(this%periodn(i1), 0., 0., 2, 1, this%periodn(i1)%n) 
          this%periodn(i1)%iequaldistant = .true.
          do i3 = 1, nele
              istart = 1
              iend = this%periodn(i1)%n
              ! call findstartendn(this, i3, i1, istart, iend) ! Don't worry, don't need (only for efficiency)
              do i4 = 1, this%np
                 call spline_build(this%vpmgridn(i3,i4,i1), 0., 0., 2, istart, iend)
                 this%vpmgridn(i3,i4,i1)%iequaldistant = .true.
              end do
          end do
      end do
  
      if (allocated(work)) deallocate(work)
  
    end subroutine ctpgrid_calculate
  
    !********* INTERNAL FUNCTIONS ************
    real function sgrid(s)
        ! Get pphi value given pphi in units of e*psi1 - equidistant
        use paras_phy, only : ei
        use profile, only : psi1
        implicit none
        real, intent(in) :: s 
        
        sgrid = -ei * s * psi1
    
    end function sgrid
  
    real function sgridds(s)
      ! define the derivative of ds
      use paras_phy, only : ei
      use profile, only : psi1
      implicit none
      real, intent(in) :: s
  
        sgridds = -ei * psi1
  
    end function sgridds
  
  end module ctp_grid
  