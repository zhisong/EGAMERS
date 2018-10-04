! manage the PIC particles
module particles

use mpi
use trap_matrix

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
            particles_push

contains

  subroutine particles_push(dgca, gca, gca_aux, tm)
    ! get the particles push from the field
    ! INPUTS:
    ! gca - TYPE(particle_container), the state of the particles
    ! gca_aus - TYPE(particle_aux_container), the aux of the particles
    ! tm - TYPE(tmatrix), the trapped particle matrix with orbit informations
    ! RETURNS:
    ! dgca - TYPE(particle_container), the push on particle states, will be used in R-K
    implicit none
    type(particle_container) :: dgca, gca
    type(particle_container_aux) :: gca_aux
    type(tmatrix) :: tm

  end subroutine particles_push

  subroutine particles_init(this, mat3t, n)
    ! allocate the memory space for particles
    ! INPUT:
    ! this - type(particle_container)
    ! mat3t- the computed trap matrix
    ! n    - INTEGER, the number of particles in the container
    implicit none

    type(particle_container) :: this
    type(tmatrix) :: mat3t
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

  !subroutine push(gca, dgca, field)
    ! push all the particles according to the equation of motion
end module particles