! pre-calculated sin function value at equal-distance zeta grid points
! to speed up the orbit integration

module sintable
  
  use paras_num
  implicit none

  integer :: ntable
  real, allocatable, dimension(:,:), public :: sinpzeta
  
  public
  
  contains

    subroutine sintable_init(norbit, np)
      ! initialize sin table
      use paras_phy
      implicit none
      integer, intent(in) :: norbit, np
      integer :: i, p
      real :: zeta, dzeta
      
      ! destroy the data if exist
      call sintable_destroy()
      allocate(sinpzeta(norbit, np))
      dzeta = 2. * pi / real(norbit)
      do i = 1, norbit
         do p = 1, np
            zeta =  dzeta * real(i-1) + dzeta / 2.
            sinpzeta(i, p) = sin(real(p) * zeta)
         end do
      end do
    end subroutine sintable_init
      
    subroutine sintable_destroy()
      ! destroy the sin table
      if (allocated(sinpzeta)) deallocate(sinpzeta)
    end subroutine sintable_destroy
  end module sintable
