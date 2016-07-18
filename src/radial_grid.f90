! generate radial grid points of Hermite elements

module radial_grid

  use paras_num
  use matrix_module
  implicit none

  public
  real, allocatable, dimension(:), public  :: rgrid
  integer, public :: nelement, nele
  
  contains

    subroutine rgrid_init(nelementin, ntype, xr1, sig1, xr2, sig2)
      ! initialize radial grids
      ! INPUTS : number of elements, radial grid accumulation point and width
      !          type : 1 (equal-distance), 2 (gaussian)
      !          only type 1 is implemented currently

      implicit none

      integer, intent(in) :: nelementin, ntype
      real, intent(in) :: xr1, sig1, xr2, sig2

      real :: drstep
      integer :: i1
      
      if ((nelementin .le. 2) .or. (nelementin .gt. nmaxelement)) then
         nelement = 30  ! Set to default
         write(*,*) 'Number of radial grid point must be > 2 and <', nmaxelement
         write(*,*) 'Set to default (30)'
      else
         nelement = nelementin
         nele = 2 * nelement - 2
      end if
      
      ! destroy if already allocated
      call rgrid_destroy()

      allocate(rgrid(nelementin))

      if (ntype .eq. 1) then
         drstep = 1. / real(nelement-1)
         do i1 = 1, nelement
            rgrid(i1) = drstep * real(i1 - 1)
         end do
      end if
    end subroutine rgrid_init
   
    subroutine rgrid_destroy()
      ! destroy radial grid
      implicit none
      if (allocated(rgrid)) deallocate(rgrid)
    end subroutine rgrid_destroy

  end module radial_grid
