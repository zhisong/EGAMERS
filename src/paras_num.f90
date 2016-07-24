! Numerical parameters

module paras_num

  implicit none

  ! drift orbit options
  ! maximum time steps in solving the drift orbit
  integer, parameter :: nmaxts = 200000
  ! maximum iteration step to find the close point of the orbit
  integer, parameter :: nmaxorbitclose = 10
  ! time step in solving drift orbit, for normal particles  
  ! (in unit of cycotron peroid)
  real            :: dtorbitn = 0.5
  ! time step in solving drift orbit, for near t/p boundary
  real            :: dtorbitb = 0.2
  ! switch to Cartesian cooridinate if  r is smaller than
  real, parameter :: rcartesian = 5e-3
  ! precision of the poloidal orbit frequency calculation (in theta)
  real, parameter :: errtheta = 1e-5
  ! number of sampling points in orbit ingetral calculation (must be odd)
  integer         :: norbitintsample = 2000
  ! the orbit is considered lost if r is greater than 
  real, parameter :: rlost = 1.0001

  ! precision of finding rstart
  real, parameter :: errrstart = 1e-12
  ! precision of finding rstart for trap edge
  real, parameter :: errrstartmin = 1e-4
  ! max iteration of finding rstart
  integer, parameter :: maxitrstart = 400
  ! precision in finding trapped/passing boundary
  real, parameter :: errtpbound = 1e-12
  ! max iteration number in finding t/p boundary

  real, parameter :: maxittpbound = 400
  ! maximum number of orbit harmonics
  integer, parameter :: npmax = 5
  ! maximum number of finite elements
  integer, parameter :: nmaxelement = 100
  ! step size for derivative calculation
  real, parameter :: ddx = 0.0001
  ! the min difference between trapedge and traplost energy to 
  ! take into account, in unit of mub0
  real, parameter :: mindiff_el = 1e-5

  ! maximum grid points in a cubic spline interpolation
  integer, parameter :: nmaxspline = 1000
  ! maximum energy root can be found for a given period
  integer, parameter :: nmaxroot = 10

  ! number of grid points in integrals of MHD matrixs
  integer, parameter :: nmhdintegral = 10

  ! the displacement used to calculate the derivative of T(lambda)
  real, parameter :: dlambdapctg = 0.0001

  ! the number of radial grid to spline the q profile
  integer, parameter :: nqspline = 200
  
  ! maximum of eigenvalue to calculate in one run
  integer, parameter:: neigenmax = 10
  ! the tolerance in finding eigenvalue
  real              :: erreig = 1e-4
  ! maximum iteration in finding eigenvalue
  integer           :: nmaxit = 20

end module paras_num
