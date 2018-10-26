! test particle tracker
module test_particles

use mpi
use field, only : field_vector
use trap_grid
use field
use particles

  implicit none

  private
  integer, parameter :: NSTATE = 4
  integer, parameter :: NGRID_ID = 3

  integer :: istart, iend

  type(field_vector) :: efield_test
  type(particle_container) :: dotgc_test, gcsum_test, gcorg_test
  type(particle_container), public :: gc_test
  type(particle_container_aux), public :: gc_test_aux
  type(tgrid), public :: tg_test

! Namelist variables
  integer, public :: nparticles_test = 100000   ! total number of particles (sum over all cpus)
  integer, public :: ksteps_test = 10000        ! total number of steps
  real,    public :: dt_test = 1e-6             ! time step (in seconds)
  integer, public :: nscreen_test = 1000        ! screen output inteval in steps
  integer, public :: nsnappart_test  = 10       ! particle output inteval in steps (will overwrite, only use for hot start)
  integer, public :: mub0_test = 300            ! the slice of mub0 where test particles are on, in keV
! //////////////////

  real, dimension(:,:), allocatable, public :: lambda_t, eta_t
  real, dimension(:), allocatable, public :: t_list

  real, public :: t_test

  public :: test_pic_init, test_pic_step, test_pic_destroy
contains

  subroutine test_pic_init()
  ! initialize the test particle simulation
    implicit none
    integer :: nmarkers, nmarkers_each, myid, ncpus

    ! get number of cpus
    ncpus = mpi_get_ncpus()
    myid  = mpi_get_rank()

    ! set time to zero
    t_test = 0

    ! initialize the fields
    call field_vector_init(efield_test)

    nmarkers_each = CEILING(float(nparticles_test) / float(ncpus))
    nmarkers = nmarkers_each * ncpus
    istart = nmarkers_each * myid + 1
    iend   = nmarkers_each * (myid + 1)

    call particles_init(gc_test, nmarkers, istart, iend)
    call particles_init(dotgc_test, nmarkers, istart, iend)
    call particles_init(gcsum_test, nmarkers, istart, iend)
    call particles_init(gcorg_test, nmarkers, istart, iend)

    call particles_aux_init(gc_test_aux, gc_test)

    ! load particles
    call particles_loading_trap(gc_test, gc_test_aux, tg_test, 0.0, 0.2)

  end subroutine test_pic_init

  subroutine test_pic_destroy
  ! destroy the objects
    call field_vector_destroy(efield_test)
    call particles_destroy(gc_test)
    call particles_destroy(dotgc_test)
    call particles_destroy(gcsum_test)
    call particles_destroy(gcorg_test)
    call particles_aux_destroy(gc_test_aux)
  end subroutine test_pic_destroy

  subroutine test_pic_step()
  ! evolve the test particles for one step
    use radial_grid, only : nele
    integer :: rkstep, i1
    real :: c0, c1, trk
    real, dimension(nele) :: vtemp

    call fill_efield(t_test)

    do rkstep = 1, 4

      ! evolve the particle
      call particles_push_trap(dotgc_test, vtemp, gc_test, gc_test_aux, efield_test, tg_test)

      if (rkstep.eq.1) then
        c0 = 1.0/6.0
        c1 = 0.5
      
        gcsum_test%state(:,:) = gc_test%state(:,:)
        gcorg_test%state(:,:) = gc_test%state(:,:)

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
        gcsum_test%state(:,:) = gcsum_test%state(:,:) + dt_test * c0 * dotgc_test%state(:,:)
        gc_test%state(:,:) = gcorg_test%state(:,:) + dt_test * c1 * dotgc_test%state(:,:)
        call fill_efield(t_test + dt_test * c1)
      else
        gc_test%state(:,:) = gcsum_test%state(:,:) + dt_test * c0 * dotgc_test%state(:,:)
      end if
    end do
  
    t_test = t_test + dt_test
  
  end subroutine test_pic_step


  subroutine fill_efield(time)
  ! fill in the field_vector using data read from file
    implicit none
    real, intent(in) :: time

    real :: dt_data
    integer :: ipos

    if (mpi_is_master()) then
    ! field is only stored at master node
      dt_data = t_list(2) - t_list(1)
      ipos = CEILING(time / dt_data) + 1
      if (time .eq. 0.0) then
        efield_test%lambda(:) = lambda_t(:, 1)
        efield_test%eta(:) = eta_t(:, 1)
      else if (ipos .gt. SIZE(t_list)) then
        efield_test%lambda(:) = 0.0
        efield_test%eta(:) = 0.0 
      else
      ! linear interp
        efield_test%lambda(:) = lambda_t(:, ipos-1) + (time - t_list(ipos-1)) / dt_data &
            * (lambda_t(:, ipos) - lambda_t(:, ipos-1))
        efield_test%eta(:) = eta_t(:, ipos-1) + (time - t_list(ipos-1)) / dt_data &
            * (eta_t(:, ipos) - eta_t(:, ipos-1))
      end if
    end if

    ! broadcast
    call field_vector_bcast(efield_test) 
  end subroutine fill_efield

end module test_particles
