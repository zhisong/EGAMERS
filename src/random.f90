! random number generator
module random
  use mpi, only : mpi_get_ncpus, mpi_get_rank
  implicit none
  
  logical, private :: initialized = .false.

  interface get_rand
    module procedure &
      get_rand_real_scalar, &
      get_rand_real_1d
  end interface get_rand

  public :: get_rand

contains

  subroutine get_rand_real_scalar(arr)
    implicit none
    real :: arr
    real :: temp
    integer :: i1, rank, ncpus
    rank = mpi_get_rank()
    ncpus = mpi_get_ncpus()
    
    call check_rand_init()
    
    do i1 = 0, rank
      call random_number(temp)
    end do

    arr = temp

    do i1 = rank + 1, ncpus-1
      call random_number(temp)
    end do

  end subroutine get_rand_real_scalar

  subroutine get_rand_real_1d(arr)
    implicit none
    real, dimension(:) :: arr
    real, allocatable, dimension(:) :: temp
    integer :: i1, rank, ncpus
    rank = mpi_get_rank()
    ncpus = mpi_get_ncpus()
    
    call check_rand_init()
    
    allocate(temp(SIZE(arr)))

    do i1 = 0, rank
      call random_number(temp)
    end do

    arr(:) = temp(:)

    do i1 = rank + 1, ncpus-1
      call random_number(temp)
    end do

    if (allocated(temp)) deallocate(temp)

  end subroutine get_rand_real_1d

  subroutine check_rand_init()
    implicit none
    if (.not. initialized) then
      call RANDOM_SEED()
      initialized = .true.
    end if
  end subroutine check_rand_init

end module random