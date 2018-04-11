! a matrix type for further calculation

module matrix_module

  use interfaces

  public

  type, public :: matrix

     complex, dimension(:,:), allocatable, public :: data
     integer, public :: ncol, nrow

   contains

     final :: matrix_destroy

  end type matrix
  
  interface operator(+)
     ! matrix adding
     module procedure matadd
  end interface operator(+)

  interface operator(*)
     ! matrix times scalar
     module procedure matmultiply
  end interface operator(*)

  interface assignment(=)
     ! matrix give value to another matrix
     module procedure matgivevalue
  end interface assignment(=)

contains

  subroutine matgivevalue(m1, m2)

    implicit none

    type(matrix), intent(out) :: m1
    type(matrix), intent(in) :: m2

    if (allocated(m2%data)) then
       
       call matrix_init(m1, m2%nrow, m2%ncol)
       call czcopy(m2%nrow*m2%ncol, m2%data, 1, m1%data, 1)
    end if
  end subroutine matgivevalue
    
  function matadd(m1, m2) result(m3)

    implicit none
    type(matrix), intent(in) :: m1, m2
    type(matrix) :: m3

    integer :: i1, i2

    if ((m1%ncol .eq. m2%ncol) .and. (m1%nrow .eq. m2%nrow) .and. &
         allocated(m1%data) .and. allocated(m2%data)) then

       call matrix_init(m3, m1%nrow, m1%ncol)
       call czcopy(m1%nrow*m1%ncol, m1%data, 1, m3%data, 1)
       call czaxpy(m1%nrow*m1%ncol, (1.,0.), m2%data, 1, m3%data, 1)

    else
       write(*,*) 'matadd : error'
    end if
  end function matadd


  function matmultiply(s1, m1) result(m3)

    implicit none
    type(matrix), intent(in) :: m1
    complex, intent(in) :: s1
    type(matrix) :: m3

    integer :: i1, i2

    if (allocated(m1%data)) then

       call matrix_init(m3, m1%nrow, m1%ncol)
       call czcopy(m1%nrow * m1%ncol, m1%data, 1, m3%data, 1)
       call czscal(m1%nrow * m3%ncol, s1, m3%data, 1)

    else
       write(*,*) 'matmultiply : error'
    end if

  end function matmultiply

  subroutine matrix_init(this, nrow, ncol)
    ! allocate a matrix

    implicit none
    type(matrix) :: this
    integer, intent(in) :: nrow, ncol
    integer :: i1, i2
    logical :: need_allocate

    if (allocated(this%data)) then
       if (this%nrow==nrow .and. this%ncol==ncol) then
          ! do nothing
          need_allocate = .false.
       else
          call matrix_destroy(this)
          need_allocate = .true.
       end if
    else
       need_allocate = .true.
    end if

    if (need_allocate) then
       if ((nrow .gt. 0) .and. (ncol .gt. 0)) then
          allocate(this%data(nrow, ncol))
          this%ncol = ncol
          this%nrow = nrow

          call matrix_clear(this)
       else
          write(*,*) 'MATRIX: nrow and ncol must be greater than 0'
       end if
    end if

  end subroutine matrix_init

  subroutine matrix_destroy(this)
    ! destroy the matrix
    
    implicit none
    type(matrix) :: this

    if (allocated(this%data)) deallocate(this%data)
  end subroutine matrix_destroy

  subroutine matrix_clear(this)
    ! clear the matrix
    implicit none
    type(matrix) :: this
    integer :: i1,i2

    if (.not. allocated(this%data)) return

    this%data(:,:) = 0.
       
  end subroutine matrix_clear

end module matrix_module
