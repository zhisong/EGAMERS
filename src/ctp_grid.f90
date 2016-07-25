  ! a grid for counter-passing particles for a single mub0
  
module ctp_grid
  ! "n" grids are Normal equal-distant energy grids 
  !     from trap edge to t/p boundary or lost boundary
  ! "b" grid are special grids for TYPE I t/p Boundary
  !     exponentially close to the t/p boundary
  !
  ! TYPE:
  ! ctpgrid           -  contains all the grids and splines
  !
  ! SUBROUTINES:
  ! (public)
  ! ctpgrid_init      - allocate arrays and splines, fill energy grid
  !                     must call first

  use radial_grid, only : nele
  use spline_module
  implicit none

  ! type of ctp particle grid
  ! (which is exactly the same as trapped grid "tgrid")
  type, public :: ctpgrid

     ! finite element weight for n grid and for t/p boundary
     ! indexes : element label, pth resonance, pphi grid
     type(spline), dimension(:,:,:), allocatable, public :: vpmgridn, vpmgridb
     ! the most inside and most outside element indexes for a certain pphi
     integer, dimension(:), allocatable, public :: ielementmin, ielementmax
     integer, dimension(:,:), allocatable, public :: ielementminn, ielementmaxn, ielementminb, ielementmaxb

     ! orbit period grid for n and b
     type(spline), dimension(:), allocatable, public :: periodn, periodb

     ! number of pphi grid for n and b
     integer, public :: npphin, npphib
     ! start pphi index of b grid
     integer, public :: ibstart
     ! the lost boundary is lower than t/p boundary at
     integer, public :: iloststart
     
     ! pphi grid
     real, dimension(:), allocatable, public :: pphigrid
     ! pphi grid start and end
     real, public :: pphistart, pphiend
     integer, public :: ipphistart, ipphiend

     ! number of energy grid for n and b grid
     integer, public :: neen, neeb
     ! energy of n/b grid switch
     real, dimension(:), allocatable, public :: enbswitch

     ! store grid data for the lowest energy and t/p boundary
     real, dimension(:), allocatable, public :: etrapedge, etpbound, elost

     ! mub0
     real, public :: mub0

     ! harmonic count
     integer, public :: np

  end type ctpgrid

  private
  public :: ctpgrid_init !, ctpgrid_destroy, ctpgrid_calculate, ctpsgrid

contains

  subroutine ctpgrid_init(this, mub0, eehigh, eelow, npphin, neen, neeb, eeendb, np, itpbound)
    ! initiate the object, create the grid
    ! INPUT : mub0, the upper limit of energy, the lower limit of energy
    !         number of pphi grid, number of n energy grid,
    !         number of b energy grid, b grid upper limit,
    !         number of harmonics, treat the t/p boundary(1) or not(0)
    use paras_num, only : mindiff_el
    use radial_grid
    use orbit_classify
    implicit none

    type(ctpgrid) :: this
    real, intent(in) :: mub0, eehigh, eelow, eeendb
    integer, intent(in) :: npphin, neen, neeb, np, itpbound

    this%mub0 = mub0
    this%np = np
    this%neen = neen
    this%neeb = neeb

  end subroutine ctpgrid_init

end module ctp_grid
