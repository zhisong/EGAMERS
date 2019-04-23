! the evolution of the field

module field

use mpi
use matrix_module
use radial_grid, only : nele

implicit none

private
type(matrix) :: mat1, mat2, mat1_qr ! MHD matrics M, N and QR decompose of M
integer, allocatable, dimension(:) :: pmat ! QR decompose indices of M
complex, allocatable, dimension(:) :: vwork, tdotv ! temp working space

real, public :: gamma_d = 0.0 ! the "ad-hoc" damping rate
type, public :: field_vector

  real, allocatable, dimension(:) :: lambda ! the finite element coefficients 
  real, allocatable, dimension(:) :: eta    ! the time derivatives of lambda
  integer :: nele ! number of elements, same as nele in module radial_grid

end type field_vector

public :: field_vector_init, field_vector_destroy, field_vector_bcast, &
          field_evolve, field_matrix_init, field_matrix_destroy, field_energy

contains

  subroutine field_vector_init(evector, irandom, ampl, ampldt)
    ! initialize the field vector
    use random
    implicit none

    type(field_vector) :: evector
    logical, intent(in), optional :: irandom
    real, intent(in), optional :: ampl, ampldt

    logical :: iirandom
    real :: aampl, aampldt

    if (.not. PRESENT(irandom)) then
        iirandom = .false.
    else
        iirandom = irandom
    end if

    if (.not. PRESENT(ampl)) then
        aampl = 0.0
    else
        aampl = ampl
    end if

    if (.not. PRESENT(ampldt)) then
        aampldt = 0.0
    else
        aampldt = ampldt
    end if

    evector%nele = nele
    
    allocate(evector%lambda(nele))
    allocate(evector%eta(nele))

    ! if we need to randomize the initial field
    if (iirandom) then
      call get_rand(evector%lambda)
      call get_rand(evector%eta)

      evector%lambda(:) = (evector%lambda(:) - 0.5) * 2.0 * aampl
      evector%eta(:) = (evector%eta(:) - 0.5) * 2.0 * aampldt

    end if

  end subroutine

  subroutine field_vector_destroy(this)
    ! destroy the field vector
    implicit none
    type(field_vector) :: this

    if (ALLOCATED(this%lambda)) deallocate(this%lambda)
    if (ALLOCATED(this%eta)) deallocate(this%eta)

  end subroutine field_vector_destroy

  subroutine field_vector_bcast(this)
    implicit none
    type(field_vector) :: this

    call mpi_bcast_array(this%lambda, 0)
    call mpi_bcast_array(this%eta, 0)
  end subroutine field_vector_bcast

  subroutine field_matrix_init()
    ! initialize the field mhd matrices
    use interfaces, only : czgetrf
    use mhd_matrix
    implicit none
    
    integer :: info

    if (mpi_is_master()) then
      ! compute the MHD matrices
      call getmat1(mat1)
      call getmat2(mat2)
      mat1_qr = mat1
      
      ! allocate the working space
      allocate(pmat(nele))
      allocate(vwork(nele))
      !allocate(tdotv(nele))

      ! QR decompose mat1 for later use
      call czgetrf(nele, nele, mat1_qr%data, nele, pmat, info)
      if (info .ne. 0) then
        write(*,*) 'ERROR IN CZGETRF, INFO =', info
      end if
    end if

  end subroutine

  subroutine field_matrix_destroy()
    ! destroy the field matrices
    implicit none

    if (mpi_is_master()) then
      call matrix_destroy(mat1)
      call matrix_destroy(mat2)
      call matrix_destroy(mat1_qr)
      if (allocated(pmat)) deallocate(pmat)
      if (allocated(vwork)) deallocate(vwork)
      !if (allocated(tdotv)) deallocate(tdotv)
    end if

  end subroutine field_matrix_destroy

  subroutine field_evolve(dotevector, evector, vfast)
    ! calcuate the right hand side of Runge-Kutta for field variables
    ! INPUTS:
    ! evector - TYPE(field_vector), the field at last time step
    ! vfast - REAL,dimenstion(nele), contribution from the fast ions
    ! RETURNS:
    ! dotevector - TYPE(field_vector), the right hand side of R-K
    use interfaces, only : czgetrs
    implicit none
    
    real, dimension(:), intent(in) :: vfast
    type(field_vector) :: evector
    type(field_vector) :: dotevector

    integer :: info

    if (.not. mpi_is_master()) return
    
    ! compute the right hand side of M dot(eta) = - (N lambda + vfast)
    vwork = - MATMUL(mat2%data, evector%lambda) - vfast
    ! inverse the M matrix to get dot(eta) = -M^(-1) (N lambda + vfast)
    call czgetrs('N', nele, 1, mat1_qr%data, nele, pmat, vwork, nele, info)
    if (info .ne. 0) then
       write(*,*) 'ERROR IN CZGETRS, INFO =', info
    end if 
    
    ! compute dot(eta) and dot(lambda)
    dotevector%lambda(:) = evector%eta(:) 
    dotevector%eta(:) = real(vwork(:)) - 2.0 * gamma_d * evector%eta(:)

  end subroutine field_evolve

  real function field_energy(evector)
    ! computes the field energy
    implicit none

    type(field_vector) :: evector
    real :: kinetic_energy, potential_energy
    
    kinetic_energy = SUM(evector%eta * MATMUL(mat1%data, evector%eta))
    potential_energy = SUM(evector%lambda * MATMUL(mat2%data, evector%lambda))
    
  end function field_energy

end module field