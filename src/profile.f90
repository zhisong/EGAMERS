! Specify bulk profiles here
module profile

  use paras_phy
  use spline_module
  implicit none
  
  public

  ! on axis electron temperature
  real :: te0 = 3.  !keV
  ! on axis bulk ion temperature
  real :: ti0 = 1.  !keV

  ! profile type : te, ti, ni, q
  integer :: ite_type = 1
  integer :: iti_type = 1
  integer :: ini_type = 1
  integer ::  iq_type = 1
  !      i**_type == 1, polynomial (5th order)
  !               == 2, Gaussian (not used for q profile)
  !               == 3, input on grid points (not implemented)

  ! polynomial te profile, normalized to on axis (5th order)
  real, dimension(6) :: tepoly = (/1.,0.,0.,0.,0.,0./)
  ! polynomial ti profile, normalized to on axis (5th order)
  real, dimension(6) :: tipoly = (/1.,0.,0.,0.,0.,0./)
  ! polynomial ni profile, normalized to on axis (5th order)
  real, dimension(6) :: nipoly = (/1.,0.,0.,0.,0.,0./)
  ! polynomial q profile
  real, dimension(6) :: qpoly = (/2.,0.,0.,0.,0.,0./)

  ! Gaussian te, ti, ni profile, profile width
  real :: te_deltar = 0.5
  real :: ti_deltar = 0.5
  real :: ni_deltar = 0.5 
  
  real :: psi1 = 0.25 * log(2.) * 3.5 * 1**2

  type(spline) :: psir
  
  contains

    subroutine profile_init()
      ! must call after reading namelist
      use paras_num, only : nqspline
      implicit none

      real, dimension(nqspline) :: qr, r2
      real :: r, dr2
      integer :: i1, i2

      te0 = te0 * eunit * 1000.
      ti0 = ti0 * eunit * 1000.

      ! interval
      dr2 = 1 / float(nqspline - 1)
      
      do i1 = 1, nqspline
         ! use r**2 instead of r as the variable
         r2(i1) = dr2 * float(i1 - 1)
         ! r = sqrt(r**2)
         r = sqrt(r2(i1))

         if (iq_type .eq. 1) then
            ! q given by spline
            qr(i1) = qpoly(6)
            do i2 = 2, 6
               qr(i1) = qr(i1) * r + qpoly(7 - i2)
            end do
         else
            write(*,*) 'iq_type invalid'
            stop
         end if

      end do

      ! initialize the spline
      call spline_init(psir, nqspline)
      psir%iequaldistant = .true.
      
      psir%x(1) = 0.
      psir%y(1) = 0.
      
      do i1 = 2, nqspline
         psir%x(i1) = r2(i1)
         psir%y(i1) = (1./qr(i1-1) + 1./qr(i1)) / 4. * dr2 * B0 * a**2 &
              + psir%y(i1-1)
      end do

      call spline_build(psir, 0., 0., 2, 1, psir%n)

      psi1 = psir%y(psir%n)
         
    end subroutine profile_init

    real function q(r)
      ! the q profile
      implicit none

      real, intent(in) :: r

      real :: r2, dpsidr2

      r2 = r**2

      dpsidr2 = spline_interpd1(psir, r2)
      q = 1. / 2. / dpsidr2 * B0 * a**2
      
    end function q

    real function dqdr(r)
      ! radial derivative of the q profile
      implicit none
      real, intent(in) :: r
      real, dimension(3) :: abltg

      real :: r2, dpsidr2, dum, qq, d2psi

      r2 = r**2

      dum = spline_interp(psir, r2, abltg)
      dpsidr2 = abltg(1)
      d2psi = abltg(2)
      qq = 1. / 2. / dpsidr2 * B0 * a**2

      dqdr = - 4.* r * qq**2 * d2psi  / (B0 * a**2)
      
    end function dqdr

    real function psi(r)
      ! the poloidal flux profile as a function of r
      implicit none
      real, intent(in) :: r
      
      psi = spline_interp1(psir, r**2)
      
    end function psi

    real function psitor(psiin)
      ! r as a function of the poloidal flux function psi r(psi)
      use paras_num, only : nmaxroot
      implicit none
      real, intent(in) ::  psiin
      real, dimension(nmaxroot) :: xroots
      integer :: nroots
      
      if (psiin .eq. 0) then
         psitor = 0.
         return
      else if (psiin .eq. psi1) then
         psitor = 1.
         return
      end if

      call spline_find(psir, psiin, 1, psir%n, nroots, xroots)

      if (nroots .ne. 1) then
         write(*,*) 'error in psitor, found',nroots, 'root(s)'
         stop
      end if
      psitor = sqrt(xroots(1))
      
    end function psitor
     
    real function ni(r)
      ! bulk ion density as a function of r, normalized to on-axis value
      implicit none
      real, intent(in) :: r
      integer :: i1
      
      if (ini_type .eq. 1) then
         ni = nipoly(6)
         do i1 = 2, 6
            ni = ni * r + nipoly(7 - i1)
         end do
      else if (ini_type .eq. 2) then
         ni = exp( - r**2 / ni_deltar**2)
      else
         ni = 1.
      end if
    end function ni

    real function ti(r)
      ! bulk ion temperature as a function of r, normalized to on-axis value
      implicit none
      real, intent(in) :: r
      integer :: i1
      
      if (iti_type .eq. 1) then
         ti = tipoly(6)
         do i1 = 2, 6
            ti = ti * r + tipoly(7 - i1)
         end do
      else if (iti_type .eq. 2) then
         ti = exp( - r**2 / ti_deltar**2)
      else
         ti = 1.
      end if
      
    end function ti

    real function te(r)
      ! electron temperature as a function of r, normalized to on-axis value
      implicit none
      real, intent(in) :: r
      integer :: i1
      
      if (ite_type .eq. 1) then
         te = tepoly(6)
         do i1 = 2, 6
            te = te * r + tepoly(7 - i1)
         end do
      else if (ite_type .eq. 2) then
         te = exp( - r**2 / te_deltar**2)
      else
         te = 1.
      end if

    end function te

    real function omega_gam(r)
      ! thermal GAM frequency as a function of r
      implicit none
      real, intent(in) :: r

      real :: omg2, tau_e
      
      tau_e = te0 / ti0 * te(r) / ti(r)
      
      omg2 = 2. * ti0 * ti(r) / mib / R0**2 * (7./4. + tau_e)
      omg2 = omg2 * (1. + 2. * (4.*tau_e**2 + 16.*tau_e + 23)&
           / (7. + 4.*tau_e)**2 / q(r)**2)
 
      omega_gam = sqrt(omg2)
    end function omega_gam
  end module profile
