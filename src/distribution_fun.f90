! specify the distribution function

module distribution_fun

  use paras_phy
  use profile
  implicit none

  public

  ! fast particle energy/pitch distribution type
  ! 1 - bi-Maxwellian, 2 - slowing down 
  integer :: ifast_pdf_type = 1
  
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

  ! Slowing down dist. parameters
  real    :: lambda0 = 0.5                   ! pitch angle
  real    :: dlambda = 0.1                   ! pitch angle width delta_lambda
  real    :: dpsi    = 0.5
  real    :: eec                             ! critical beam-ion energy
  real    :: ctp_norm = 5. * (10. ** (-41.)) ! Normailization factor for slowing down dist function

contains

  subroutine distribution_fun_init()

    implicit none
    
    real :: c_beam ! Constant to find critical beam-ion energy

    pphi0 = - ei * psi(r_peak)

    if (ifast_pdf_type .eq. 2) then
      c_beam = 3*pi**(1./2.) * zf * (2.014*mp)**(3./2.) /4./(me**(1./2.))/(2.014*mp)
      eec = te0 * c_beam**(2./3.)
    end if

  end subroutine distribution_fun_init
  
  real function f0(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres, eedenom, pitch, exppsi

    ! slowing down distribution - ifast_pdf_type .eq. 2
    if (ifast_pdf_type .eq. 2) then
      ! Slowing down distribution: f = 1/(E^(3/2) + Ec^(3/2)) * exp(- (lamba - lambda0)^2 / dlambda^2)
      eedenom = ee**(3./2.) + eec**(3./2.)
      pitch = (mub0/ee - lambda0)**2.
      exppsi = pphi / ei / psi1 / dpsi

      f0 = ctp_norm * exp(-pitch / dlambda**2.) * exp(exppsi) / eedenom
    else 
      ! Bi-Maxwellian
      Bres = B0 * R0 / Rres
      nf = exp((pphi - pphi0) / ei / dpphi_nf / psi1)
      tper = tper0 * 1000. * eunit * exp(pphi / ei / dpphi_tper / psi1)
      tpar = tpar0 * 1000. * eunit * exp(pphi / ei / dpphi_tpar / psi1)

      cons = 2. * (mi / 2. / pi / tper)**(3./2.)
      mubres = mub0 / B0 * Bres

      f0 = cons * nf * exp(-mubres/tper - abs(ee-mubres)/tpar)
    end if
  end function f0

  
  real function df0de(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres, denom, term1, term2, term3, pitch, exppsi
    
    ! slowing down distribution
    if (ifast_pdf_type .eq. 2) then
      ! d/dE{ 1/(E^3/2 + Ec^3/2) * exp[-(mub0/E - lambda0)^2 / dlambda^2]}

      denom = 2. * ee**3. *(ee**(3./2.) + eec**(3./2.))**2. * dlambda**2.

      term1 = 4. * mub0**2. *(ee**(3./2.) + eec**(3./2.))
      term2 = -3. * ee**(7./2.) * dlambda**2.
      term3 = -4. * mub0 * ee * (ee**(3./2.) + eec**(3./2.)) * lambda0

      pitch = (mub0/ee - lambda0)**2.
      exppsi = pphi / ei / psi1 / dpsi

      df0de =  ctp_norm * exp(-pitch / dlambda**2.)* exp(exppsi) * (term1 + term2 + term3) / denom

    else
      ! Bi-Maxwellian

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

    end if

    ! print *, "ee, mub0 ", ee/1000./eunit, mub0/1000/eunit
    ! print *, "f0, df0de, d2f0f2e ", f0(ee,mub0,pphi), df0de, d2f0de2(ee, mub0, pphi)

  end function df0de

  real function d2f0de2(ee, mub0, pphi)

    use paras_phy
    use profile
    implicit none

    real, intent(in) :: ee, mub0, pphi

    real :: nf, tper, tpar, cons, mubres, Bres, eecterm, denom, term1, term2,&
            term3, term4, term5, term6, term7, term8, term9, pitch, exppsi
    
    ! slowing down distribution
    if (ifast_pdf_type .eq. 2) then
      ! slowing down second derivative
      eecterm = ee**(3./2.) + eec**(3./2.)
      denom = 4. * ee**6. * eecterm**3. * dlambda**4.

      term1 = 16.*mub0**4. * (ee**3. + eec**3.)
      term2 = 32.*mub0**4.*ee**(3./2.)*eec**(3./2.)
      term3 = -48.*mub0**2. *ee**5. * dlambda**2.
      term4 = -72.*mub0**2. *ee**(7./2.)*eec**(3./2.) *dlambda**2.
      term5 = -24.*mub0**2. *ee**2. *eec**3. *dlambda**2.
      term6 = 15.* ee**7. * dlambda**4.
      term7 = -3.*ee**(11./2.) * eec**(3./2.) * dlambda**4.
      term8 = 8.*mub0*ee*eecterm*lambda0* (-4.*eecterm*mub0**2. + ee**2. * dlambda**2. * (5.*ee**(3./2.) + 2.*eec**(3./2.)))
      term9 = 16.*ee**2. *eecterm**2. * mub0**2. * lambda0**2.

      pitch = (mub0/ee - lambda0)**2.
      exppsi = pphi / ei / psi1 / dpsi

      d2f0de2 = ctp_norm * exp(-pitch / dlambda**2.)* exp(exppsi) * (term1+term2+term3+term4+term5+term6+term7+term8+term9) / denom

    else
      ! Bi-Maxwellian
      Bres = B0 * R0 / Rres
      nf = exp((pphi-pphi0) / ei / dpphi_nf / psi1)
      tper = tper0 * 1000. * eunit * exp(pphi / ei / dpphi_tper / psi1)
      tpar = tpar0 * 1000. * eunit * exp(pphi / ei / dpphi_tpar / psi1)

      cons = 2. * (mi / 2. / pi / tper)**(3./2.)
      mubres = mub0 / B0 * Bres
      
      d2f0de2 = 1. / tpar**2 * cons * nf * exp(-mubres/tper - abs(ee-mubres)/tpar)
    end if
  end function d2f0de2

end module distribution_fun
