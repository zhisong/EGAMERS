! the PIC module, evolve the particles and the field
module pic

use mpi
use eigen, only : tm
use field
use particles
use trap_matrix
use radial_grid, only : nele

implicit none

  type(field_vector), private :: dotefield, efieldsum, efieldorg
  type(field_vector), public :: efield

  type(particle_container), allocatable, dimension(:), private :: dotgca, gcasum, gcaorg
  type(particle_container), allocatable, dimension(:), public :: gca
  type(particle_container_aux), allocatable, dimension(:), public :: gca_aux

  real, public :: t = 0.0      ! this is the time
  integer, public :: nmarkers   ! number of markers on each cpu

! //// Read from the namelist ////
  integer, public :: nparticles = 100000   ! total number of particles (sum over all cpus)
  integer, public :: ksteps = 10000        ! total number of steps
  integer, public :: dt_adjust = 0         ! ==1:adjust dt according to the orbit frequencies (1/5 of maximum orbit fqc) ==0:no adjust
  real, public    :: dt = 1e-6             ! time step (in seconds)
  real, public    :: initampl = 1e-10      ! initial amplitute (max of random)
  real, public    :: initampldt = 1e-5     ! initial time derivative of amplitute (max of random)
  integer, public :: nsnapfield = 100      ! field output inteval in steps
  integer, public :: nsnappart  = 1000     ! particle output inteval in steps
! ////////////////////////////////

  public :: pic_init, pic_step, pic_destroy

contains

  subroutine pic_init()
    ! initialize the pic simulation
    implicit none
    integer :: i1

    ! set time to zero
    t = 0

    ! initialize MHD matrix
    if (mpi_is_master()) call field_matrix_init()

    ! initialize the fields
    call field_vector_init(efield, .true., initampl, initampldt )
    if (mpi_is_master()) then
      call field_vector_init(efieldsum)
      call field_vector_init(efieldorg)
      call field_vector_init(dotefield)
    end if

    ! initialize the particles
    ! the number of particles on each mub0 grid
    nmarkers = CEILING(float(nparticles) / float(tm%ngrid))
    nparticles = nmarkers * tm%ngrid

    ! allocate gca_aux, gca, dotgca, gcasum and gcaorg
    allocate(gca_aux(tm%ngrid))
    allocate(gca(tm%ngrid))
    allocate(dotgca(tm%ngrid))
    allocate(gcasum(tm%ngrid))
    allocate(gcaorg(tm%ngrid))

    do i1 = 1, tm%ngrid
      if (.not. mpi_is_my_work(lwork, i1)) cycle
      ! initialize all the particle containers
      call particles_init(gca(i1), tm,  nmarkers)
      call particles_init(dotgca(i1), tm, nmarkers)
      call particles_init(gcasum(i1), tm, nmarkers)
      call particles_init(gcaorg(i1), tm, nmarkers)
      call particles_aux_init(gca_aux(i1), nmarkers)
    end do

  end subroutine pic_init
  
  subroutine pic_destroy()
    ! destroy all the variables in use
    implicit none
    integer :: i1
    
    if (mpi_is_master()) call field_matrix_destroy()

    call field_vector_destroy(efield)
    if (mpi_is_master()) then
      call field_vector_destroy(efieldsum)
      call field_vector_destroy(efieldorg)
      call field_vector_destroy(dotefield)
    end if

    do i1 = 1, tm%ngrid
      if (.not. mpi_is_my_work(lwork, i1)) cycle
      ! initialize all the particle containers
      call particles_destroy(gca(i1))
      call particles_destroy(dotgca(i1))
      call particles_destroy(gcasum(i1))
      call particles_destroy(gcaorg(i1))
      call particles_aux_destroy(gca_aux(i1))
    end do

  end subroutine pic_destroy

  subroutine pic_step()
    ! advance in pic step
    ! INPUTS:
    ! tm - TYPE(tmatrix), the trapped matrix
    implicit none
    
    integer :: rkstep, i1
    real :: c0, c1
    real, dimension(nele) :: vfast

    vfast(:) = 0.0

    do rkstep = 1, 4

      vfast(:) = 0.0

      ! evolve the particle
      do i1 = 1, tm%ngrid

        if (.not. mpi_is_my_work(lwork, i1)) cycle

        !call paticles_push(dgca, gca)

        if (rkstep.eq.1) then
          c0 = 1.0/6.0
          c1 = 0.5
      
          gcasum(i1)%state(:,:) = gca(i1)%state(:,:)
          gcaorg(i1)%state(:,:) = gca(i1)%state(:,:)

        else if(rkstep.eq.2)then
          c0 = 1.0/3.0
          c1 = 0.5
        else if(rkstep.eq.3)then
          c0 = 1.0/3.0
          c1 = 1.0
        else if(rkstep.eq.4)then
          c0 = 1.0/6.0
        end if
        
        if(rkstep.ne.4)then

          gcasum(i1)%state(:,:) = gcasum(i1)%state(:,:) + dt * c0 * dotgca(i1)%state(:,:)
          gca(i1)%state(:,:) = gcaorg(i1)%state(:,:) + dt * c1 * dotgca(i1)%state(:,:)
          
        else
          
          gca(i1)%state(:,:) = gcasum(i1)%state(:,:) + dt * c0 * dotgca(i1)%state(:,:)
          
        end if

      end do

      ! evolve the field at master node
      if (mpi_is_master()) then

        call field_evolve(dotefield, efield, vfast)
        

        if (rkstep.eq.1) then
          c0 = 1.0/6.0
          c1 = 0.5
      
          efieldsum%lambda(:) = efield%lambda(:)
          efieldsum%eta(:) = efield%eta(:)
          efieldorg%lambda(:) = efield%lambda(:)
          efieldorg%eta(:) = efield%eta(:)

        else if(rkstep.eq.2)then
          c0 = 1.0/3.0
          c1 = 0.5
        else if(rkstep.eq.3)then
          c0 = 1.0/3.0
          c1 = 1.0
        else if(rkstep.eq.4)then
          c0 = 1.0/6.0
        end if
        
        if(rkstep.ne.4)then

          efieldsum%lambda(:) = efieldsum%lambda(:) + dt * c0 * dotefield%lambda(:)
          efieldsum%eta(:) = efieldsum%eta(:) + dt * c0 * dotefield%eta(:)

          efield%lambda(:) = efieldorg%lambda(:) + dt * c1 * dotefield%lambda(:)
          efield%eta(:) = efieldorg%eta(:) + dt * c1 * dotefield%eta(:)
        else
          
          efield%lambda(:) = efieldsum%lambda(:) + dt * c0 * dotefield%lambda(:)
          efield%eta(:) = efieldsum%eta(:) + dt * c0 * dotefield%eta(:)
          
        end if

      end if

      

      call field_vector_bcast(efield)
    enddo

    t = t + dt

  end subroutine pic_step

end module pic