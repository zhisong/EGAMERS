! Define the data type (which is data type specific)
#ifdef MPI
#define MPI_REAL_TYPE MPI_REAL8
#define MPI_COMPLEX_TYPE MPI_COMPLEX16
#define MPI_INTEGER_TYPE MPI_INTEGER
#define MPI_LOGICAL_TYPE MPI_LOGICAL
#endif

module mpi

  implicit none

#ifdef MPI
  include 'mpif.h'
#endif
  
  private

  type, public :: workload
     integer, public :: work_table(1:1024)
     integer, public :: n_work
  end type workload
  
  integer :: rank, nrank
  integer :: ierror, tag !, status(MPI_STATUS_SIZE)

  public :: mpi_start, mpi_end, mpi_is_master, mpi_get_ncpus, mpi_get_rank, &
       mpi_sync, mpi_sum_matrix, mpi_allocate_work, mpi_is_my_work, mpi_whose_work, &
       mpi_bcast_spline, mpi_sum_vector

  interface mpi_bcast_scalar
     module procedure &
          mpi_bcast_scalar_integer, &
          mpi_bcast_scalar_real, &
          mpi_bcast_scalar_complex, &
          mpi_bcast_scalar_logical
  end interface mpi_bcast_scalar

  interface mpi_bcast_array
     module procedure &
          mpi_bcast_array_real_1d, &
          mpi_bcast_array_real_2d, &
          mpi_bcast_array_integer_1d, &
          mpi_bcast_array_integer_2d
  end interface mpi_bcast_array
  
  public :: mpi_bcast_scalar, mpi_bcast_array
  
contains

  ! Wrapper, initializing MPI
  subroutine mpi_start()

    implicit none

#ifdef MPI
    
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    
    if (ierror > 0) then
       write(*,*) 'MPI initializing error'
       stop
    end if
#else
    nrank = 1
    rank = 0
#endif
    
  end subroutine mpi_start

  ! Wrapper, finializing MPI
  subroutine mpi_end()

    implicit none
#ifdef MPI
    call MPI_FINALIZE(ierror)
#endif
  end subroutine mpi_end
  
  ! Wrapper for MPI_
  subroutine mpi_sync()
    integer :: ierr
#ifdef MPI    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr > 0 .and. rank==0) write(*,*) 'MPI sync failed'
#endif
  end subroutine mpi_sync
    
  ! Wrapper, if the current node is the master node (rank==0)
  function mpi_is_master()

    logical :: mpi_is_master

    if (rank == 0) then
       mpi_is_master = .true.
    else
       mpi_is_master = .false.
    endif

  end function mpi_is_master

  ! Wrapper, get the number of cpus
  pure function mpi_get_ncpus()
    integer :: mpi_get_ncpus
    mpi_get_ncpus = nrank
  end function mpi_get_ncpus

  ! Wrapper, get the rank
  pure function mpi_get_rank()
    integer :: mpi_get_rank
    mpi_get_rank = rank
  end function mpi_get_rank
  
  !**********************************
  ! Interfaces for broadcasting
  !**********************************

  ! Wrapper, broadcasting complex scalar
  subroutine mpi_bcast_scalar_complex(data_in, root)

    complex, intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr
#ifdef MPI    
    call MPI_BCAST(data_in, 1, MPI_COMPLEX_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_scalar_complex

  
  ! Wrapper, broadcasting real scalar
  subroutine mpi_bcast_scalar_real(data_in, root)

    real, intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr
#ifdef MPI    
    call MPI_BCAST(data_in, 1, MPI_REAL_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_scalar_real

  ! Wrapper, broadcasting integer scalar
  subroutine mpi_bcast_scalar_integer(data_in, root)

    integer, intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr
#ifdef MPI    
    call MPI_BCAST(data_in, 1, MPI_INTEGER_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_scalar_integer

  ! Wrapper, broadcasting integer scalar
  subroutine mpi_bcast_scalar_logical(data_in, root)

    logical, intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr
#ifdef MPI    
    call MPI_BCAST(data_in, 1, MPI_LOGICAL_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_scalar_logical
  
  ! Wrapper, broadcasting real array
  subroutine mpi_bcast_array_real_1d(data_in, root)

    real, dimension(:), intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr, nelements

    nelements = size(data_in)
#ifdef MPI    
    call MPI_BCAST(data_in, nelements, MPI_REAL_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_array_real_1d

  subroutine mpi_bcast_array_real_2d(data_in, root)

    real, dimension(:,:), intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr, nelements

    nelements = size(data_in)
#ifdef MPI    
    call MPI_BCAST(data_in, nelements, MPI_REAL_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_array_real_2d

    ! Wrapper, broadcasting integer array
  subroutine mpi_bcast_array_integer_1d(data_in, root)

    integer, dimension(:), intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr, nelements

    nelements = size(data_in)
#ifdef MPI    
    call MPI_BCAST(data_in, nelements, MPI_INTEGER_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_array_integer_1d

  subroutine mpi_bcast_array_integer_2d(data_in, root)

    integer, dimension(:,:), intent(inout) :: data_in
    integer, intent(in) :: root
    integer :: ierr, nelements

    nelements = size(data_in)
#ifdef MPI    
    call MPI_BCAST(data_in, nelements, MPI_INTEGER_TYPE, root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (ierr > 0) write(*,*) 'MPI broadcasting error in cpu', rank
#endif    
  end subroutine mpi_bcast_array_integer_2d

  subroutine mpi_bcast_spline(data_in, root)
    ! broadcast the spline
    use spline_module

    type(spline) :: data_in
    integer, intent(in) :: root
    
    call mpi_bcast_array_real_1d(data_in%a, root)
    call mpi_bcast_array_real_1d(data_in%b, root)
    call mpi_bcast_array_real_1d(data_in%c, root)
    call mpi_bcast_array_real_1d(data_in%d, root)
    call mpi_bcast_array_real_1d(data_in%x, root)
    call mpi_bcast_array_real_1d(data_in%y, root)
    call mpi_bcast_scalar_integer(data_in%n, root)
    call mpi_bcast_scalar_logical(data_in%iequaldistant, root)

  end subroutine mpi_bcast_spline

  !**********************************
  ! MPI reduce wrapper
  !**********************************

  subroutine mpi_sum_vector(vec_in, vec_out, root)

    use matrix_module
    implicit none
   
    real, dimension(:), intent(inout) :: vec_in, vec_out
    integer, intent(in) :: root
    integer :: nelement, ierr

    nelement = size(vec_in)
 
#ifdef MPI
    call MPI_reduce(vec_in, vec_out, nelement, MPI_REAL_TYPE, &
         MPI_SUM, root, MPI_COMM_WORLD, ierr)
    if (ierr>0 .and. mpi_is_master()) write(*,*) 'MPI_reduce failed'
#else
    vec_out(:) = vec_in(:)
#endif
  end subroutine mpi_sum_vector

  ! summation over matrix
  subroutine mpi_sum_matrix(mat_in, mat_out, root)

    use matrix_module
    implicit none
   
    type(matrix), intent(inout) :: mat_in, mat_out
    integer, intent(in) :: root
    integer :: nelement, ierr

    nelement = size(mat_in%data)
 
#ifdef MPI
    call MPI_reduce(mat_in%data, mat_out%data, nelement, MPI_COMPLEX_TYPE, &
         MPI_SUM, root, MPI_COMM_WORLD, ierr)
    if (ierr>0 .and. mpi_is_master()) write(*,*) 'MPI_reduce failed'
#else
    mat_out = mat_in
#endif
  end subroutine mpi_sum_matrix
      
  !**********************************
  ! Work allocation section
  !**********************************
  
  ! allocate n_work works to the number of cpus
  subroutine mpi_allocate_work(lwork, n_work)

    implicit none
    
    integer, intent(in) :: n_work
    type(workload) :: lwork

    integer :: icpu, ii
    logical :: ascending

    ascending = .true.
    icpu = 0
    lwork%n_work = n_work
    
    do ii = 1, n_work
       lwork%work_table(ii) = icpu
       if (ascending) then
          icpu = icpu + 1
       else
          icpu = icpu - 1
       end if

       if (icpu==nrank) then
          icpu = nrank - 1
          ascending = .not. ascending
       elseif (icpu==-1) then
          icpu = 0
          ascending = .not. ascending
       endif

    end do
    
  end subroutine mpi_allocate_work

  ! query the worktable lwork to see if the current cpu should do the work with id=workid 
  function mpi_is_my_work(lwork, workid)

    implicit none

    type(workload), intent(in) :: lwork
    integer,intent(in) :: workid

    logical :: mpi_is_my_work
    
    if (lwork%work_table(workid) == rank) then
       mpi_is_my_work = .true.
    else
       mpi_is_my_work = .false.
    end if

    if (workid > lwork%n_work .or. workid < 1) mpi_is_my_work = .false.

  end function mpi_is_my_work

    ! query the worktable lwork to see if the current cpu should do the work with id=workid 
  function mpi_whose_work(lwork, workid)

    implicit none

    type(workload), intent(in) :: lwork
    integer, intent(in) :: workid
    integer :: mpi_whose_work

    if (workid > lwork%n_work .or. workid < 1) then
      mpi_whose_work = -1
    else
      mpi_whose_work = lwork%work_table(workid)
    end if

  end function mpi_whose_work
  
end module
