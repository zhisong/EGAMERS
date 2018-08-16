! specify the distribution function

module distribution_fun

  use paras_phy
  use profile
  implicit none

  public
  
  ! on axis fast ion density to bulk ion density
  real :: nf_ratio = 0.02

  ! perpendicular and parallel temperature on axis
  real :: tper0 = 400  ! keV
  real :: tpar0 = 100  ! keV
  
  ! perpendicular and parallel temperture radial width
  ! T_pper = T_pper0 exp(pphi / ei / dpphitper / psi(1.))
  real :: dpphi_tper = 0.5
  real :: dpphi_tpar = 0.5
  
  ! density width in pphi
  real :: dpphi_nf = 0.2
  ! density peak location
  real :: r_peak = 0.0
  ! the pphi for this peak
  real :: pphi0

  ! ICRH resonance major radius and field strength
  real :: Rres = 2.8

contains

  subroutine distribution_fun_init()

    implicit none

    pphi0 = - ei * psi(r_peak)

  end subroutine distribution_fun_init
  
  real function f0(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres
    
    Bres = B0 * R0 / Rres
    nf = exp((pphi - pphi0) / ei / dpphi_nf / psi1)
    tper = tper0 * 1000. * eunit * exp(pphi / ei / dpphi_tper / psi1)
    tpar = tpar0 * 1000. * eunit * exp(pphi / ei / dpphi_tpar / psi1)

    cons = 2. * (mi / 2. / pi / tper)**(3./2.)
    mubres = mub0 / B0 * Bres
   
    f0 = cons * nf * exp(-mubres/tper - abs(ee-mubres)/tpar)
  end function f0

  
  real function df0de(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres
    
    Bres = B0 * R0 / Rres
    nf = exp((pphi-pphi0) / ei / dpphi_nf / psi1)
    tper = tper0 * 1000. * eunit * exp(pphi / ei / dpphi_tper / psi1)
    tpar = tpar0 * 1000. * eunit * exp(pphi / ei / dpphi_tpar / psi1)

    cons = 2. * (mi / 2. / pi / tper)**(3./2.)
    mubres = mub0 / B0 * Bres
    
    df0de = - 1. / tpar * cons * nf * exp(-mubres/tper - abs(ee-mubres)/tpar)
    if (mubres .ge. ee) then
       df0de = - df0de
    end if
  end function df0de

  real function d2f0de2(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres
    
    Bres = B0 * R0 / Rres
    nf = exp((pphi-pphi0) / ei / dpphi_nf / psi1)
    tper = tper0 * 1000. * eunit * exp(pphi / ei / dpphi_tper / psi1)
    tpar = tpar0 * 1000. * eunit * exp(pphi / ei / dpphi_tpar / psi1)

    cons = 2. * (mi / 2. / pi / tper)**(3./2.)
    mubres = mub0 / B0 * Bres
    
    d2f0de2 = 1. / tpar**2 * cons * nf * exp(-mubres/tper - abs(ee-mubres)/tpar)
    
  end function d2f0de2

end module distribution_fun
