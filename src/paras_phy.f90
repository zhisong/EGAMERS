! Basic parameters

module paras_phy

  implicit none

  ! ////// PHYSICS CONSTANTS //////
  ! pi
  real, parameter :: pi = 3.141592653589
  ! proton mass
  real, parameter :: mp = 1.66e-27
  ! electron mass
  real, parameter :: me = 9.109e-31
  ! unit charge
  real, parameter :: eunit = 1.60e-19

  ! ////// FAST IONS //////
  ! fast ion mass number
  real :: af = 1
  ! fast ion charge number
  real :: zf = 1
  ! fast ion mass
  real :: mi
  ! fast ion charge
  real :: ei

  ! ////// BULK IONS //////
  ! bulk ion mass number
  real :: ai  = 2.
  ! bulk ion mass
  real :: mib

  ! ////// GEOMETRY //////
  ! major radius (m)
  real :: R0 = 3.0
  ! minor radius (m)
  real :: a = 1.0
  ! magnetic field strengh (T)eex
  real :: B0 = 3.5
  ! inverse aspect ratio
  real :: eps

  ! ////// OTHERS //////
  ! cycotron period (s)
  real :: Tc
  ! constants to accelerate calculations
  real :: or0sqrtmi
  real :: driftc
  complex :: sqrt3o2i

contains
  
  subroutine paras_phy_init()
    ! initialize the physics parameters and constants
    implicit none

    mi = af * mp
    ei = zf * eunit
    mib = ai * mp
    eps = a / R0
    Tc = 2 * pi * mi / ei / B0
    or0sqrtmi = 1. / R0 / sqrt(mi)
    driftc = 1. / ( ei * B0 * R0 * a)
    sqrt3o2i = 3**0.5 / 2 * (0., 1.)
  end subroutine paras_phy_init

end module paras_phy
