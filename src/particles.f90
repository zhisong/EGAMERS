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
    ! state(x,i), x - 1: deltaf, 2: energy, 3: theta
    real, allocatable, dimension(:,:) :: state
  end type particle_container

  type, public :: particle_container_aux
    ! some additional properties of the particles
    real, public :: mub0
    integer, public :: imub0
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

  subroutine particles_push_trap(dgc, dv, gc, gc_aux, efield, tm, imub0)
    ! get the particles push from the field
    ! INPUTS:
    ! gc - TYPE(particle_container), the state of the particles
    ! gc_aux - TYPE(particle_aux_container), the aux of the particles
    ! efield - TYPE(field_vector), the field of current time step
    ! tm - TYPE(tmatrix), the trapped particle matrix with orbit informations
    ! imub0 - INTEGER, the index of mub0 in the trap matrix
    ! RETURNS:
    ! dgc - TYPE(particle_container), the push on particle states, will be used in R-K
    ! dv  - real(nele), the particle contribution to fluid

    use paras_phy
    use distribution_fun, only : df0de, nf_ratio
    use trap_grid
    use trap_matrix, only : tmatrix
    use field, only : field_vector
    implicit none
    type(particle_container) :: dgc, gc
    type(particle_container_aux) :: gc_aux
    type(tmatrix) :: tm
    type(field_vector) :: efield
    integer, intent(in) :: imub0
    real, dimension(:), intent(out) :: dv

    integer :: i1, m, p, nele, ipphi
    real :: period, omegab, ee, theta, deedt, df0, factor
    real, dimension(3) :: dperiod
    real, allocatable, dimension(:,:) :: vpm
    real, allocatable, dimension(:) :: sinvalue, vwork

    nele = efield%nele

    ! allocate the working vector and matrix
    allocate(vpm(nele,tm%grid(imub0)%np))
    allocate(sinvalue(tm%grid(imub0)%np))
    allocate(vwork(nele))

    dv(:) = 0
    vwork(:) = 0

    do i1 = 1, gc%n
      ! ignore inactive particles (out of bound particles)
      if (gc_aux%active(i1) <= 0) cycle

      ! shorthand
      ee = gc%state(1, i1)
      theta = gc%state(2, i1)
      ipphi = gc_aux%grid_id(1, i1)

      ! the period and omegab
      call getperiod(tm%grid(imub0), ee, ipphi, period, dperiod)
      ! check if out of the bound
      if (period < 0.0) then
        ! if yes, set inactive
        gc_aux%active(i1) = 0
        cycle
      end if
      omegab = 2. * pi / period
      

      do p = 1, tm%grid(imub0)%np
        do m = 1, nele
          vpm(m, p) = getvpm(tm%grid(imub0), ee, ipphi, p, m)
        end do
        sinvalue(p) = SIN(p * theta)
      end do
      
      vwork(:) = MATMUL(vpm(:,:), sinvalue(:))
      deedt = ei * omegab * SUM(vwork(:) * efield%eta(:))
      df0 = df0de(ee, gc_aux%mub0, tm%grid(imub0)%pphigrid(ipphi))
      
      ! set dgc
      ! de/dt
      dgc%state(1, i1) = deedt
      dgc%state(2, i1) = omegab
      dgc%state(3, i1) = - deedt * df0

      ! set dv
      factor = nf_ratio * a**2  / mi**2 / B0 / mib / R0 * omegab
      dv(:) = dv(:) + vwork(:) * gc%state(3, i1) * gc_aux%weight(i1) * factor
      
    end do

  end subroutine particles_push_trap

  subroutine particles_loading_trap(gc, gc_aux, tm, imub0, tpbound_ratio)
    ! set the initial state of the particles
    ! INPUTS:
    ! gc - TYPE(particle_container), the state of the particles
    ! gc_aux - TYPE(particle_container_aux), the aux of the particles
    ! tm - TYPE(tmatrix), the computed trap matrix
    ! imub0 - INTEGER, the index of mub0 in the trap matrix
    ! tpbound_ratio : REAL, optional, the ratio of particles allocated to the near t/p bound special grid, default = 0.2
    use paras_phy!, only : pi
    use paras_num, only : TP_RATIO_DEFAULT
    use spline_module
    use random
    use trap_grid
    use trap_matrix, only : tmatrix
    implicit none

    type(particle_container) :: gc
    type(particle_container_aux) :: gc_aux
    type(tmatrix) :: tm
    integer, intent(in) :: imub0
    real, intent(in), optional :: tpbound_ratio

    integer, dimension(:), allocatable :: neachpphin, neachpphib
    integer :: npartn, npartb
    integer :: i1, i2, icounter, ntotalegrid, iup
    real :: tp_ratio, ee, dee, weight
    real :: dmub0, dpphi, period, omegab
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
    gc_aux%mub0 = tm%mub0_table(imub0)
    gc_aux%imub0 = imub0
    ! calculate the grid size
    dmub0 = abs(tm%mub0_table(2) - tm%mub0_table(1))
    dpphi = abs(tm%grid(imub0)%pphigrid(2) -  tm%grid(imub0)%pphigrid(1))
    ! clear the particle counter
    icounter = 0

    ! the number of particles on n/p grid
    if (tm%grid(imub0)%npphib > 0 ) then
      npartb = FLOOR(float(gc%n) * tp_ratio)
      npartn = gc%n - npartb
    else
      npartn = gc%n
      npartb = 0
    end if

    ! now collect (approximatly) the number of the e-pphi space
    ! this should approximately equal to the area of the phase space
    ntotalegrid = 0
    do i1 = 1, tm%grid(imub0)%npphin
      ntotalegrid = ntotalegrid + tm%grid(imub0)%periodn(i1)%n - 1
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tm%grid(imub0), i1)) ntotalegrid = ntotalegrid - 1
    end do

    ! now calculate the number of particles on each pphi (n)grid
    allocate(neachpphin(tm%grid(imub0)%npphin))
    do i1 = 1, tm%grid(imub0)%npphin - 1
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tm%grid(imub0), i1)) then
        neachpphin(i1) = FLOOR(float(npartn) / float(ntotalegrid) * float(tm%grid(imub0)%periodn(i1)%n - 2))
      else
        neachpphin(i1) = FLOOR(float(npartn) / float(ntotalegrid) * float(tm%grid(imub0)%periodn(i1)%n - 1))
      end if
    enddo
    neachpphin(tm%grid(imub0)%npphin) = npartn - SUM(neachpphin(1:tm%grid(imub0)%npphin-1))
    
    ! load the particles on ngrid
    do i1 = 1, tm%grid(imub0)%npphin
      ! don't count the last grid point for the n grid if it has a Type I t/p bound
      if (istype1(tm%grid(imub0), i1)) then
        iup = tm%grid(imub0)%periodn(i1)%n - 1
      else
        iup = tm%grid(imub0)%periodn(i1)%n
      end if
      ! equdistant in energy, let's ignore the two ends
      dee = (tm%grid(imub0)%periodn(i1)%x(iup) - tm%grid(imub0)%periodn(i1)%x(1)) / float(neachpphin(i1) + 1)

      ! load particles
      do i2 = 1, neachpphin(i1)
        ! increase the counter
        icounter = icounter + 1
        ! load the initial energy
        ee = tm%grid(imub0)%periodn(i1)%x(1) + dee * float(i2)
        gc%state(1,icounter) = ee
        ! fill in the ipphi
        gc_aux%grid_id(1, icounter) = i1
        ! set to be on n grid
        gc_aux%grid_id(2, icounter) = 0
        ! set to be active
        gc_aux%active(icounter) = 1
        ! set the weight of the particle
        call getperiod(tm%grid(imub0), ee, i1, period, dperiod)
        omegab = 2.0 * pi / period
        weight =  (dmub0) * dpphi * dee * period
        gc_aux%weight(icounter) = weight 

      end do
    end do

    ! //// load the b grid ////
    if (tm%grid(imub0)%npphib > 0 ) then
      allocate(neachpphib(tm%grid(imub0)%npphib))
      
      ! now calculate the number of particles on each pphi (b)grid
      neachpphib(1:tm%grid(imub0)%npphib - 1) = FLOOR(float(npartb) / float(tm%grid(imub0)%npphib))
      neachpphib(tm%grid(imub0)%npphib) = npartb - SUM(neachpphib(1:tm%grid(imub0)%npphib - 1))
      if (imub0 == 12) write(*,*) icounter, gc%n, tm%grid(imub0)%npphib
      do i1 = 1, tm%grid(imub0)%npphib
        iup = tm%grid(imub0)%periodb(i1)%n
        ! equdistant in log energy, let's ignore the two ends
        dee = (tm%grid(imub0)%periodb(i1)%x(iup) - tm%grid(imub0)%periodb(i1)%x(1)) / float(neachpphib(i1) - 1)
        ! load particles
        do i2 = 1, neachpphib(i1)
          ! increase the counter
          icounter = icounter + 1
          ! load the initial energy
          gc%state(1,icounter) = eelogtoee(tm%grid(imub0), tm%grid(imub0)%periodb(i1)%x(1) + dee * float(i2), i1)
          ! fill in the ipphi
          gc_aux%grid_id(1, icounter) = indexb2n(tm%grid(imub0), i1)
          ! set to be on b grid
          gc_aux%grid_id(2, icounter) = 1
          ! set to be active
          gc_aux%active(icounter) = 1
          ! set the weight of the particle
          weight = 0.0
          gc_aux%weight(icounter) = weight
        end do
      end do
    end if
    ! randomize the starting angle
    call get_rand(gc%state(2,:)) 
    gc%state(2,:) = gc%state(2,:) * 2.0 * pi

    if (ALLOCATED(neachpphin)) deallocate(neachpphin)
    if (tm%grid(imub0)%npphib > 0 ) then
      if (ALLOCATED(neachpphib)) deallocate(neachpphib)
    end if

  end subroutine particles_loading_trap


  subroutine particles_init(this, n)
    ! allocate the memory space for particles
    ! INPUT:
    ! this - type(particle_container)
    ! n    - INTEGER, the number of particles in the container
    implicit none

    type(particle_container) :: this
    integer, intent(in) :: n

    this%n = n
    allocate(this%state(NSTATE, n))

  end subroutine particles_init

  subroutine particles_aux_init(this, n)
    ! allocate the memory space for particle aux
    ! INPUT:
    ! this - type(particle_container_aux)
    ! n    - INTEGER, the number of particles in the container
    implicit none

    type(particle_container_aux) :: this
    integer, intent(in) :: n

    this%n = n
    allocate(this%grid_id(NGRID_ID, n))
    allocate(this%weight(n))
    allocate(this%active(n))

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