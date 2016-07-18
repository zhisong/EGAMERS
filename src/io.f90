! Input and output

module io

  use paras_num, only : neigenmax
  implicit none

  ! file units
  integer, parameter :: ionamelist = 10  ! namelist file
  integer, parameter :: iofield    = 20  ! electric field file
  integer, parameter :: iofqc      = 21  ! continuum and global mode fqc file
  integer, parameter :: iomap      = 22
  integer, parameter :: ioorbit    = 23  ! orbit file

  ! number of radial points put in the output
  integer, parameter :: nfieldoutput = 1000

  ! Input quantities

  ! //////// NAMELIST RUNS ////////
  ! WHAT to do in the code
  ! MODE of running
  !  imode = 1   run as an eigenvalue solver and solve for E_r and omega
  !        = 2   calculate the map of bounce frequency for a given mub0,
  !                  and the resonance lines for a given frequency
  !        = 3   only plot the orbit for a given (E, mub0, Pphi)
  integer :: imode = 1

  ! the follows are used when imode = 1 (eigenvalue solver)
  ! number of trial eigenvalues
  integer :: neigen = 1
  ! list of trial of eigenvalues
  complex omegain(neigenmax)

  ! the follows are used when imode = 2 (frequency map mode)
  real    :: mub0_in = 100

  ! the follows are used when imode = 3 (orbit mode)
  ! the mub0 to calculate particle frequency/orbit in keV (imode = 2,3)
  !real    :: mub0_in = 100
  ! the energy to calculate particle orbit in keV (imode = 3)
  real    :: energy_in = 100
  ! the pphi to calculate particle orbit in psi1 (imode = 3)
  real    :: pphi_in = -0.5
  ! co(1)/counter(-1) passing 
  integer :: vsign_in = 1
  ! plot over or below the t/p boundary
  !   ienergy_tpbound = 0  use energy_in as the energy
  !                   = -1 energy = ee_tpbound - mub0 * exp(-energy_log_in)
  !                   = 1  energy = ee_tpbound + mub0 * exp(-energy_log_in)
  integer :: ienergy_tpbound = 0
  real :: energy_log_in = 25.

  ! //////// NAMELIST NUMS ////////
  
  ! enable the trapped particles or not
  integer :: ienable_trap    = 1
  ! enable the co-passing particle or not
  ! (currently not implemented, please leave to 0)
  integer :: ienable_cop     = 0
  ! enable the counter-passing particles or not
  ! (currently not implemented, please leave to 0)
  integer :: ienable_ctp     = 0

  ! enable the special treatment for trapped/passing boundary or not
  ! THE SPECIAL TREATMENT IS ENFORCED IF TRAPPED PARTICLES ARE ENABLED
  integer :: ienable_tpbound = 1

  ! the number of orbit harmonics for trapped particles (>=1)
  integer :: np_trap         = 1
  ! the number of orbit harmonics for co-passing particles (>=1)
  integer :: np_cop          = 1
  ! the number of orbit harmonics for counter-passing particles (>=1)
  integer :: np_ctp          = 1
  
  ! SOME OTHER namelist members in module paras_num

  ! number of sampling points on the orbit to calculate Vpm (or)
  ! number of sampling points if imode == 3
  ! integer :: norbitintsample = 2000
  
  ! time step in solving drift orbit, in unit of cycotron peroid
  ! for normal grid
  ! real    :: dtorbitn = 0.5
  ! for near t/p boundary
  ! real    :: dtorbitb = 0.2

  ! the tolerance in finding eigenvalue
  ! real    :: erreig = 1e-4
  ! maximum iteration in finding eigenvalue
  ! integer :: nmaxit = 20

  ! //////// NAMELIST GRID ////////
  
  ! number of radial grid points
  integer :: nradial_grid = 30
  
  ! number of mub0 grid points for trapped particles
  integer :: ngtrap_mub0 = 10
  ! mub0 grid starts and ends at (in keV)
  real    :: trap_mub0start = 50.
  real    :: trap_mub0end   = 1000.

  ! number of normal energy grid points for trapped particles
  integer :: ngtrap_energyn = 50
  ! number of special energy grid points for trapped particles
  integer :: ngtrap_energyb = 20
  ! special energy grid ends at E = Etpbound - mub0 exp(-trapebend)
  real    :: trap_ebend = 20.

  ! number of Pphi grid points for trapped particles
  integer :: ngtrap_pphi = 50

  ! //////// NAMELIST PHYS ////////
  
  ! namelist memebers in module paras_phy
  ! real    :: af = 1.     ! fast ion mass number (times of proton mass)
  ! real    :: zf = 1.     ! fast ion charge number (times of electron charge)
  ! real    :: ai = 2.     ! bulk ion mass number (times of proton mass)
  ! real    :: R0 = 3.0    ! major radius (meter)
  ! real    :: a  = 1.0    ! minor radius (meter)
  ! real    :: B0 = 3.5    ! on axis field strength (Tesla)
  
  ! //////// NAMELIST PROF ////////
  ! namelist members in module prof
  ! real    :: te0 = 3.    ! on axis T_e in keV
  ! real    :: ti0 = 1.    ! on axis T_i in keV
  !
  ! integer :: ite_type = 1 ! T_e profile input type
  ! integer :: iti_type = 1 ! T_i profile input type
  ! integer :: ini_type = 2 ! n_i profile input type
  ! integer :: iq_type  = 1 ! q   profile input type
  !      i**_type == 1, polynomial (5th order)
  !               == 2, Gaussian (not used for q profile)
  !               == 3, input on grid points (not implemented)
  !
  ! real, dimension(6) :: tepoly, tipoly, nipoly, qpoly
  !      polynomial coefficients (tepoly, tipoly, nipoly: normalized to on axis)
  !      e.g.  T_e = te0 * (tepoly(0) + tepoly(1) r + ... + tepoly(6) r**5)
  !              q = qpoly(0) + qpoly(1) r + ... + qpoly(6) r**5
  ! real :: ni_deltar, te_deltar, ti_deltar = 0.5
  !      Gaussian parameter
  !      e.g.  T_e = te0 * exp(-r**2 / te_deltar**2)
      
  ! //////// NAMELIST FAST ///////
  
  ! namelist members in module distribution_fun
  ! integer :: ifast_pdf_type = 1! energy/pitch distribution type
  !        =  1,  bi-Maxwellian ~ exp[-mu Bres/Tper -|mu Bres - E|/Tpar]
  !        =  2,  Slowing down, ~ 1/(E^3/2+ Ec^3/2) 
  !                                * exp[ -(muB0/E - lambda0)^2/delta_lambda^2]
  !        =  3,  Bump-on-tail, ~
  ! integer :: ifast_nf_type   = 1 ! density distribution function type
  ! integer :: ifast_tper_type = 1 ! perpendicular temperature type if bi-Max
  ! integer :: ifast_tpar_type = 1 ! parallel temperature type if bi-Max
  !        =  1,  Gaussian, e.g. ~ exp[-Pphi/(dpphi_nf * Pphi1)]
  !        =  2,  Polynomial
  !
  ! real    :: nf_ratio = 0.02     ! fast ion to bulk ion density ratio on axis 
  ! real    :: dpphi_nf = 0.2      ! density width if ifast_nf_type==1
  ! real    :: tper0 = 400.        ! on axis perpendicular temp in keV (bi-Max)
  ! real    :: tpar0 = 100.        ! on axis parallel temp in keV (bi-Max)
  ! real    :: dpphi_tper = 0.5    ! tper width (bi-Max & Gaussian)
  ! real    :: dpphi_tpar = 0.5    ! tpar width (bi-Max & Gaussian)
contains

! ////// INPUT //////

  subroutine readnamelist()
    use paras_phy
    use paras_num
    use profile
    use distribution_fun
    implicit none

    ! what to do in the code
    namelist /RUNS/ imode, neigen, omegain, energy_in, mub0_in, &
         pphi_in, vsign_in, ienergy_tpbound, energy_log_in
    ! numerical settings
    namelist /NUMS/ ienable_trap, ienable_cop, ienable_ctp, &
         ienable_tpbound, np_trap, np_cop, np_ctp, norbitintsample, &
         dtorbitn, dtorbitb, erreig, nmaxit
    ! grid settings
    namelist /GRID/ nradial_grid, ngtrap_mub0, trap_mub0start, trap_mub0end,&
         ngtrap_energyn, ngtrap_energyb, trap_ebend, ngtrap_pphi
    ! physics settings
    namelist /PHYS/ af, zf, ai, R0, a, B0
    ! plasma equilibrium profiles
    namelist /PROF/ te0, ti0, ite_type, iti_type, ini_type, iq_type, &
         tepoly, tipoly, nipoly, qpoly, te_deltar, ti_deltar, ni_deltar
    ! fast particle distribution function
    namelist /FAST/ nf_ratio, dpphi_nf, dpphi_tper, dpphi_tpar, &
         tper0, tpar0, Rres

    open(UNIT=ionamelist, FILE='namelist.in', ACTION='READ')

    read(ionamelist, nml = RUNS)
    read(ionamelist, nml = NUMS)
    read(ionamelist, nml = GRID)
    read(ionamelist, nml = PHYS)
    read(ionamelist, nml = PROF)
    read(ionamelist, nml = FAST)

    close(ionamelist)

  end subroutine readnamelist

! ////// OUTPUT //////  

  subroutine plotcontinuum()
    ! plot the bulk continuum and omega of the global mode
    ! FORMAT :
    !           n_global_mode  n_r
    !           omg_global_mode_1_re     im
    !           omg_global_mode_2_re     im
    !              .....
    !            r       omg_continuum     
    !  (e.g.)   0.1        300000         
    !           0.2        290000        
    !           0.3        280000       
    !           ...         ...          

    use profile, only : omega_gam
    use eigen, only : lambda
    implicit none

    real :: dx, x
    complex :: omega
    integer :: i1
    
    open(UNIT=iofqc, FILE='fqc.out', ACTION='WRITE')

    omega = sqrt(lambda)

    write(iofqc,*) neigen, nfieldoutput
    do i1 = 1, neigen
       write(iofqc,*) real(omega), imag(omega)
    end do

    dx = 1. / real(nfieldoutput - 1)
    do i1 = 1, nfieldoutput
       x = dx * real(i1 - 1)
       write(iofqc,*) x, omega_gam(x)
    end do

    close(UNIT=iofqc)
    
  end subroutine plotcontinuum

  subroutine printfield(v)
    ! write to file the radial field structure
    ! FORMAT:    n_global_mode n_r
    !            r       real(E_r)     Imag(E_r)
    use radial_grid
    use hermite
    implicit none

    complex, dimension(:), intent(in) :: v
    complex, dimension(:), allocatable :: er
    real :: dx, x
    complex :: value
    integer :: i1, i2

    open(UNIT=iofield, FILE='field.out', ACTION='WRITE')
    
    write(iofield,*) neigen, nfieldoutput

    allocate(er(nfieldoutput))
    dx = 1. / real(nfieldoutput - 1)
    
    do i1 = 1, nfieldoutput
       x = dx * real(i1 - 1)
       ! first grid point
       er(i1) = v(1) * c2(x, 0., -1., rgrid(1))
       ! last grid point
       er(i1) = er(i1) + v(nele) * c2(x, 1., rgrid(nelement-1), 1.2) 

       ! adding contribution from other grids
       do i2 = 2, nelement - 1
          er(i1) = er(i1) + v(i2*2-2) * c1(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
          er(i1) = er(i1) + v(i2*2-1) * c2(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
       end do

       write(iofield, *) x, real(er(i1)), imag(er(i1))
    end do

    close(UNIT=iofield)
  end subroutine printfield

  subroutine write_map_header()
    ! write the header of the output file if imode = 2
    implicit none
    
    open(UNIT=iomap, FILE='map.out', ACTION='WRITE')
    write(iomap, *) ienable_trap, ienable_cop, ienable_ctp

  end subroutine write_map_header

  subroutine close_map()
    ! close the output file if imode = 2
    implicit none
    close(UNIT=iomap)
  end subroutine close_map

  subroutine plot_tgrid_map(this)
    ! write to output file the frequency map if imode = 2
    use trap_grid
    use paras_phy, only : ei, pi, eunit
    use profile, only : psi1
    implicit none

    type(tgrid), intent(in) :: this
    integer :: i1, i2

    write(iomap,*) this%npphin, this%neen
    
    do i1 = 1, this%npphin
       do i2 = 1, this%neen
          write(iomap,*) this%pphigrid(i1)/ei/psi1, &
               this%periodn(i1)%x(i2)/eunit/1000., 2*pi/this%periodn(i1)%y(i2)
       end do
    end do

  end subroutine plot_tgrid_map
    
  subroutine printorbit(norbit, rdata, thetadata, peroiddata)
    ! write to file the orbit if imode = 3
    ! format : omega     (line 1)
    !          n         (line 2)   
    !          r, theta  (line 3 - line n+3)
    use paras_phy, only : pi
    implicit none

    integer, intent(in) :: norbit
    real, dimension(:), intent(in) :: rdata, thetadata
    real, intent(in) :: peroiddata

    real :: omegadata
    integer :: i1

    open(UNIT=ioorbit, FILE='orbit.out', ACTION='WRITE')

    if (peroiddata .eq. 0.) then
       omegadata = 0.
    else
       omegadata = 2. * pi / peroiddata
    end if

    write(ioorbit,*) omegadata
    write(ioorbit,*) norbit

    do i1 = 1, norbit
       write(ioorbit,*) rdata(i1), thetadata(i1)
    end do

    close(UNIT=ioorbit)

  end subroutine printorbit

  subroutine printorbittype(ee, mub0, pphi, vsign)
    ! classify the given orbit
    ! give information about t/p boundary and lost boundary
    use paras_phy, only : ei, eunit
    use orbit_classify
    use profile, only : psi1
    implicit none
    
    real, intent(in) :: ee, mub0, pphi
    integer, intent(in) :: vsign

    real :: pphiratio, eetpbound, eetrapedge, eetraplost
    integer :: otype, itype, ierr

    pphiratio = pphi / psi1 / ei

    ierr = 0
    if (pphiratio .le. -1.) then
       ! must be ctp
       otype = -1
       if (vsign .ne. -1) then
          ! error
          otype = -100
       end if
    else if (pphiratio .ge. 1.) then
       ! must be cop
       if (vsign .ne. 1) then
          ! error
          otype = -100
       end if
    else
       ! can be trap, ctp or cop
       eetpbound = tpbound(mub0, pphi, itype)
       if (itype .eq. -1) then
          write(*,*) 'Error finding t/p bound'
          otype = -100
       else
          if (ee .le. eetpbound) then
             ! trapped
             otype = 0
          else if (vsign .eq. -1) then
             ! ctp
             otype = -1
          else if (vsign .eq. 1) then
             ! cop
             otype = 1
          else
             ! error
             otype = -100
          end if
       end if
    end if

    if (otype .eq. 0) then
       eetrapedge = trapedge(mub0, pphi, ierr)
       eetraplost = traplost(mub0, pphi)
       if (ierr .ne. 1) then
          otype = - 100
       end if
    end if
   
    if (otype .eq. 0) then
       if (ee .lt. eetrapedge) then
          write(*,*) 'Energy too low for trap orbit'
          otype = -100
       else if (ee .gt. eetraplost) then
          write(*,*) 'Trapped-lost'
          otype = -100
       else
          write(*,*) 'Trapped'
       end if
          
       write(*,*) 'Min energy    :', eetrapedge/eunit/1000., 'keV'
       write(*,*) 'Trapped lost  :', eetraplost/eunit/1000., 'keV'
       if (itype .eq. 1) then
          write(*,*) 'T/P bound (I) :', eetpbound/eunit/1000., 'keV'
       else
          write(*,*) 'T/P bound (II):', eetpbound/eunit/1000., 'keV'
       end if
    else if (otype .eq. 1) then
       write(*,*) 'Co-passing'
    else if (otype .eq. -1) then
       write(*,*) 'Counter-passing'
    end if
  end subroutine printorbittype

end module io
