! manage the PIC particles
module particles

use mpi
use trap_matrix

  implicit none 

  private
  integer, parameter :: NSTATE = 4
  integer, parameter :: NGRID_ID = 3

  type, public :: particle_container
    real, public :: mub0
    integer, public :: imub0
    integer, public :: n
    ! state of the particles
    ! state(x,i), x - 1: deltaf, 2: energy, 3: theta
    real, allocatable, dimension(:,:) :: state
    ! where the particle is in the grid
    ! grid_id(x,i), x - 1: ipphi index, 2: n/p grid
    integer, allocatable, dimension(:,:) :: grid_id
    ! whether the particles are active 
    logical, allocatable, dimension(:) :: active
  end type particle_container

contains

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
    allocate(this%grid_id(NGRID_ID, n))
    allocate(this%active(n))

  end subroutine particles_init

  subroutine particles_destory(this)
    ! allocate the memory space for particles
    ! INPUT:
    ! this - type(particle_container)
    ! n    - INTEGER, the number of particles in the container
    implicit none
    type(particle_container) :: this

    if (allocated(this%state)) deallocate(this%state)
    if (allocated(this%grid_id)) deallocate(this%grid_id)
    if (allocated(this%active)) deallocate(this%active)

  end subroutine particles_destory

  !subroutine push(gca, dgca, field)
    ! push all the particles according to the equation of motion
end module particles