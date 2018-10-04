! the PIC module, evolve the particles and the field
module pic

use mpi
use field
use particles
use trap_matrix
use radial_grid, only : nele

implicit none

  type(field_vector), private :: dotefield, efieldsum, efieldorg

  type(field_vector), public :: efield

  real, public :: t = 0.0  ! this is the time

  integer, public :: nparticles = 100000   ! total number of particles (sum over all cpus)
  integer, public :: ksteps = 10000        ! total number of steps
  integer, public :: dt_adjust = 0         ! ==1:adjust dt according to the orbit frequencies (1/5 of maximum orbit fqc) ==0:no adjust
  real, public    :: dt = 1e-6             ! time step (in seconds)
  real, public    :: initampl = 1e-10      ! initial amplitute (max of random)
  real, public    :: initampldt = 1e-5     ! initial time derivative of amplitute (max of random)
  integer, public :: nsnapfield = 100      ! field output inteval in steps
  integer, public :: nsnappart  = 1000     ! particle output inteval in steps

  public :: pic_init, pic_step, pic_destroy

contains

  subroutine pic_init()
    implicit none
    
    if (mpi_is_master()) call field_matrix_init()

    call field_vector_init(efield, .true., initampl, initampldt )
    if (mpi_is_master()) then
      call field_vector_init(efieldsum)
      call field_vector_init(efieldorg)
      call field_vector_init(dotefield)
    end if

    t = 0

  end subroutine pic_init
  
  subroutine pic_destroy()
    implicit none
    
    if (mpi_is_master()) call field_matrix_destroy()

    call field_vector_destroy(efield)
    if (mpi_is_master()) then
      call field_vector_destroy(efieldsum)
      call field_vector_destroy(efieldorg)
      call field_vector_destroy(dotefield)
    end if

  end subroutine pic_destroy

  subroutine pic_step()
    implicit none
    
    integer :: rkstep
    real :: c0, c1
    real, dimension(nele) :: vfast

    vfast(:) = 0.0

    do rkstep = 1, 4

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