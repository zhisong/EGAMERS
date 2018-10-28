! manage the PIC particles
module particles

use mpi

  implicit none 

  private
  integer, parameter :: NSTATE = 4
  integer, parameter :: NGRID_ID = 3

  type, public :: particle_container
    ! state of the particles
    integer, public :: n
    ! number of total particles on a mub0 grid
    integer, public :: n_mub0
    ! start and end indices of particles
    integer, public :: istart, iend
    ! state(x,i), x - 1: deltaf, 2: energy, 3: theta
    real, allocatable, dimension(:,:) :: state
  end type particle_container

  type, public :: particle_container_aux
    ! some additional properties of the particles
    real, public :: mub0
    integer, public :: n
    ! the weight of the particle
    real, allocatable, dimension(:) :: weight
    ! where the particle is in the grid
    ! grid_id(x,i), x - 1: ipphi index, 2: n/p grid
    integer, allocatable, dimension(:,:) :: grid_id
    ! whether the particles are active 
    integer, allocatable, dimension(:) :: active
  end type particle_container_aux

  public :: particles_init, particles_aux_init, particles_destroy, particles_aux_destroy, &
            particles_push_trap, particles_loading_trap

contains

  subroutine particles_push_trap(dgc, dv, gc, gc_aux, efield, tg)
    ! get the particles push from the field
    ! INPUTS:
    ! gc - TYPE(particle_container), the state of the particles
    ! gc_aux - TYPE(particle_aux_container), the aux of the particles
    ! efield - TYPE(field_vector), the field of current time step
    ! tg - TYPE(tgrid), the trapped particle matrix with orbit informations
    ! RETURNS:
    ! dgc - TYPE(particle_container), the push on particle states, will be used in R-K
    ! dv  - real(nele), the particle contribution to fluid

    use paras_phy
    use distribution_fun, only : df0de, nf_ratio
    use trap_grid
    use field, only : field_vector
    implicit none
    type(particle_container) :: dgc, gc
    type(particle_container_aux) :: gc_aux
    type(tgrid) :: tg
    type(field_vector) :: efield
    real, dimension(:), intent(out) :: dv

    integer :: i1, m, p, nele, ipphi
    real :: period, omegab, ee, theta, deedt, df0, factor
    real, dimension(3) :: dperiod
    real, allocatable, dimension(:,:) :: vpm
    real, allocatable, dimension(:) :: sinvalue, vwork

    nele = efield%nele

    ! allocate the working vector and matrix
    allocate(vpm(nele,tg%np))
    allocate(sinvalue(tg%np))
    allocate(vwork(nele))

    dv(:) = 0
    vwork(:) = 0

    do i1 = 1, gc%n
      ! ignore inactive particles (out of bound particles)
      if (gc_aux%active(i1) <= 0) cycle

      ! shorthand
      ee = gc%state(i1, 1)
      theta = gc%state(i1, 2)
      ipphi = gc_aux%grid_id(i1, 1)

      ! the period and omegab
      call getperiod(tg, ee, ipphi, period, dperiod)
      ! check if out of the bound
      if (period < 0.0) then
        ! if yes, set inactive
        gc_aux%active(i1) = 0
        cycle
      end if
      omegab = 2. * pi / period
      

      do p = 1, tg%np
        do m = 1, nele
          vpm(m, p) = getvpm(tg, ee, ipphi, p, m)
        end do
        sinvalue(p) = SIN(p * theta)
      end do
      
      vwork(:) = MATMUL(vpm(:,:), sinvalue(:))
      deedt = ei * omegab * SUM(vwork(:) * efield%eta(:))
      df0 = df0de(ee, gc_aux%mub0, tg%pphigrid(ipphi))
      
      ! set dgc
      ! de/dt
      dgc%state(i1, 1) = deedt
      dgc%state(i1, 2) = omegab
      dgc%state(i1, 3) = - deedt * df0

      ! set dv
      factor = nf_ratio * a**2  / mi**2 / B0 / mib / R0 * omegab
      dv(:) = dv(:) + vwork(:) * gc%state(i1, 3) * gc_aux%weight(i1) * factor
      
    end do

  end subroutine particles_push_trap

  subroutine particles_loading_trap(gc, gc_aux, tg, dmub0, tpbound_ratio)
    ! set the initial state of the particles
    ! INPUTS:
    ! gc - TYPE(particle_container), the state of the particles
    ! gc_aux - TYPE(particle_container_aux), the aux of the particles
    ! tg - TYPE(tgrid), the computed trap grid
    ! dmub0 - REAL, the step of the mub0 grid
    ! tpbound_ratio : REAL, optional, the ratio of particles allocated to the near t/p bound special grid, default = 0.2
    use paras_phy!, only : pi
    use paras_num, only : TP_RATIO_DEFAULT
    use spline_module
    use random
    use trap_grid
    implicit none

    type(particle_container) :: gc
    type(particle_container_aux) :: gc_aux
    type(tgrid) :: tg
    real, intent(in) :: dmub0
    real, intent(in), optional :: tpbound_ratio

    integer, dimension(:), allocatable :: neachpphin, neachpphib
    integer :: npartn, npartb
    integer :: i1, i2, icounter, itotal, ntotalegrid, iup, ipos, iistart, iiend
    real :: tp_ratio, ee, dee, weight, s
    real :: dpphi, period, omegab, eelog
    real, dimension(3) :: dperiod

    if (PRESENT(tpbound_ratio)) then
      if (tpbound_ratio .ge. 0.0 .and. tpbound_ratio .le. 1.0) then
        tp_ratio = tpbound_ratio
      else
        tp_ratio = TP_RATIO_DEFAULT
      end if
    else
      tp_ratio = TP_RATIO_DEFAULT
    endif

    ! assign mub0 to particle_container_aux
    gc_aux%mub0 = tg%mub0
    ! clear the particle counter
    icounter = 0
    ! clear the total counter
    itotal = 0

    ! the number of particles on n/p grid
    if (tg%npphib > 0 ) then
      npartb = FLOOR(float(gc%n_mub0) * tp_ratio)
      npartn = gc%n_mub0 - npartb
    else
      npartn = gc%n_mub0
      npartb = 0
    end if
    ! if (mpi_is_master()) then
    !   write(*,*) gc%n_mub0, gc%n, gc%istart, gc%iend
    ! end if
    ! now collect (approximatly) the number of the e-pphi space
    ! this should approximately equal to the area of the phase space
    ntotalegrid = 0
    do i1 = 1, tg%npphin
      ntotalegrid = ntotalegrid + tg%periodn(i1)%n - 1
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tg, i1)) ntotalegrid = ntotalegrid - 1
    end do

    ! now calculate the number of particles on each pphi (n)grid
    allocate(neachpphin(tg%npphin))
    do i1 = 1, tg%npphin - 1
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tg, i1)) then
        neachpphin(i1) = FLOOR(float(npartn) / float(ntotalegrid) * float(tg%periodn(i1)%n - 2))
      else
        neachpphin(i1) = FLOOR(float(npartn) / float(ntotalegrid) * float(tg%periodn(i1)%n - 1))
      end if
    enddo
    neachpphin(tg%npphin) = npartn - SUM(neachpphin(1:tg%npphin-1))
    
    ! load the particles on ngrid
    do i1 = 1, tg%npphin
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tg, i1)) then
        iup = tg%periodn(i1)%n - 1
      else
        iup = tg%periodn(i1)%n
      end if
      ! equdistant in energy, let's ignore the two ends
      dee = (tg%periodn(i1)%x(iup) - tg%periodn(i1)%x(1)) / float(neachpphin(i1) + 1)
      ! calculate dpphi
      s = tg%s(i1) 
      dpphi = abs(tg%ds * sgridds(s))

      ! load particles
      do i2 = 1, neachpphin(i1)
        ! increase the counter
        itotal = itotal + 1
        ! if not in the range, skip
        if (itotal .lt. gc%istart .or. itotal .gt. gc%iend ) cycle
        icounter = icounter + 1
        ! load the initial energy
        ee = tg%periodn(i1)%x(1) + dee * float(i2)
        gc%state(icounter, 1) = ee
        ! fill in the ipphi
        gc_aux%grid_id(icounter, 1) = i1
        ! set to be on n grid
        gc_aux%grid_id(icounter, 2) = 0
        ! set to be active
        gc_aux%active(icounter) = 1
        ! set the weight of the particle
        call getperiod(tg, ee, i1, period, dperiod)
        weight =  (dmub0) * dpphi * dee * period
        gc_aux%weight(icounter) = weight 

      end do
    end do

    ! //// load the b grid ////
    if (tg%npphib > 0 ) then
      allocate(neachpphib(tg%npphib))
      
      ! now calculate the number of particles on each pphi (b)grid
      neachpphib(1:tg%npphib - 1) = FLOOR(float(npartb) / float(tg%npphib))
      neachpphib(tg%npphib) = npartb - SUM(neachpphib(1:tg%npphib - 1))
      do i1 = 1, tg%npphib
        iup = tg%periodb(i1)%n
        ! equdistant in log energy, let's ignore the two ends
        dee = (tg%periodb(i1)%x(iup) - tg%periodb(i1)%x(1)) / float(neachpphib(i1) + 1)
        ! calculate dpphi
        s = tg%s(i1) 
        dpphi = abs(tg%ds * sgridds(s))
        ! load particles
        do i2 = 1, neachpphib(i1)
          ! increase the counter
          itotal = itotal + 1
          ! if not in the range, skip
          if (itotal .lt. gc%istart .or. itotal .gt. gc%iend ) cycle
          icounter = icounter + 1
          ! n index
          ipos = indexb2n(tg, i1)
          ! load the initial energy
          eelog = tg%periodb(i1)%x(1) + dee * float(i2)
          ee = eelogtoee(tg, eelog, ipos)
          gc%state(icounter, 1) = ee
          ! fill in the ipphi
          gc_aux%grid_id(icounter, 1) = ipos
          ! set to be on b grid
          gc_aux%grid_id(icounter, 2) = 1
          ! set to be active
          gc_aux%active(icounter) = 1
          ! set the weight of the particle
          call getperiod(tg, ee, ipos, period, dperiod)
          weight =  (dmub0) * dpphi * period * (gc_aux%mub0 * exp(-eelog)) * dee 
          gc_aux%weight(icounter) = weight
        end do
      end do
    end if
    ! randomize the starting angle
    call get_rand(gc%state(:, 2)) 
    gc%state(:, 2) = gc%state(:, 2) * 2.0 * pi
    ! inital delta f is zero
    gc%state(:, 3) = 0.0

    if (ALLOCATED(neachpphin)) deallocate(neachpphin)
    if (tg%npphib > 0 ) then
      if (ALLOCATED(neachpphib)) deallocate(neachpphib)
    end if

  end subroutine particles_loading_trap


  subroutine particles_init(this, n, istart, iend)
    ! allocate the memory space for particles
    ! INPUT:
    ! this - type(particle_container)
    ! n    - INTEGER, the number of particles in the container
    ! istart, iend : REAl, optional, only fill in particles with indices between istart and iend

    implicit none

    type(particle_container) :: this
    integer, intent(in) :: n    
    integer, intent(in), optional :: istart, iend

    if (PRESENT(istart)) then
      this%istart = istart
    else
      this%istart = 1
    end if

    if (PRESENT(iend)) then
      this%iend = iend
    else
      this%iend = n
    end if

    ! check 
    if (this%istart .le. 0 .or. this%iend .gt. n) stop 'istart<0 or iend>n'

    this%n_mub0 = n
    this%n = this%iend - this%istart + 1

    allocate(this%state(this%n, NSTATE))

  end subroutine particles_init

  subroutine particles_aux_init(this, this_container)
    ! allocate the memory space for particle aux
    ! INPUT:
    ! this - type(particle_container_aux)
    ! this_container - type(particle_container)
    implicit none

    type(particle_container_aux) :: this
    type(particle_container) :: this_container
    
    this%n = this_container%n
    allocate(this%grid_id(this%n, NGRID_ID))
    allocate(this%weight(this%n))
    allocate(this%active(this%n))

  end subroutine particles_aux_init

  subroutine particles_destroy(this)
    ! deallocate the memory space for particles
    ! INPUT:
    ! this - type(particle_container)
    implicit none
    type(particle_container) :: this

    if (allocated(this%state)) deallocate(this%state)

  end subroutine particles_destroy

  subroutine particles_aux_destroy(this)
    ! deallocate the memory space for particle aux
    ! INPUT:
    ! this - type(particle_container_aux)
    implicit none
    type(particle_container_aux) :: this

    if (allocated(this%grid_id)) deallocate(this%grid_id)
    if (allocated(this%weight)) deallocate(this%weight)
    if (allocated(this%active)) deallocate(this%active)

  end subroutine particles_aux_destroy

  !subroutine push(gc, dgc, field)
    ! push all the particles according to the equation of motion
end module particles