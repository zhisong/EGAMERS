! the PIC diagnostics module
module diagnostics

  private

  public :: diag_active_trap_particles, &    ! the number of active trapped particles
            diag_field_energy                ! the total MHD energy

contains

  integer function diag_active_trap_particles()
    ! return the number of active particles in the system
    use trap_matrix, only : lwork
    use eigen, only : tm
    use pic, only : gc_aux
    use mpi
    implicit none
    
    integer :: i1, n_active, n_tmp

    n_active = 0
    n_tmp = 0
    do i1 = 1, tm%ngrid
      if (.not. mpi_is_my_work(lwork, i1)) cycle
      n_tmp = n_tmp + SUM(gc_aux(i1)%active(:))
    enddo

    call mpi_sum_scalar(n_tmp, n_active, 0)
    diag_active_trap_particles = n_active
    
  end function diag_active_trap_particles


  real function diag_field_energy()
    ! computes the field energy
    use field
    use pic, only : efield
    implicit none

    real :: kinetic_energy, potential_energy

    kinetic_energy = SUM(efield%eta * MATMUL(REAL(mat1%data), efield%eta))
    potential_energy = SUM(efield%lambda * MATMUL(REAL(mat2%data), efield%lambda))
    diag_field_energy = kinetic_energy + potential_energy
    
  end function diag_field_energy

end module diagnostics