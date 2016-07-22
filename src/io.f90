! Input and output

module io

  use paras_num, only : neigenmax
  use nl
  implicit none

  integer, parameter :: iofield    = 20  ! electric field file
  integer, parameter :: iofqc      = 21  ! continuum and global mode fqc file
  integer, parameter :: iomap      = 22
  integer, parameter :: ioorbit    = 23  ! orbit file
  
contains
  
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

    real :: pphiratio, eetpbound, eeedge, eelost
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
       ! stag edge and lost boundary for trapped particles
       eeedge = stagedge(mub0, pphi, 1, ierr)
       eelost = traplost(mub0, pphi)
    else if (otype .eq. 1) then
       ! stag edge and lost boundary for co-passing particles
       eeedge = stagedge(mub0, pphi, 1, ierr)
       eelost = coplost(mub0, pphi)
    else if (otype .eq. -1) then
       ! stag edge and lost boundary for ct-passing particles
       eeedge = stagedge(mub0, pphi, -1, ierr)
       eelost = ctplost(mub0, pphi)
    end if

    if (ierr .ne. 1) then
       otype = - 100
    end if

    if (otype .eq. 0) then
       if (ee .lt. eeedge) then
          write(*,*) 'Energy too low for trap orbit'
       else if (ee .gt. eelost) then
          write(*,*) 'Trapped-lost'
       else
          write(*,*) 'Trapped'
       end if
    else if (otype .eq. 1) then
       if (ee .lt. eeedge) then
          write(*,*) 'Energy too low for co-passing'
       else if (ee .gt. eelost) then
          write(*,*) 'Co-passing lost'
       else
          write(*,*) 'Co-passing'
       end if
    else if (otype .eq. -1) then
       if (ee .ge. eeedge) then
          write(*,*) 'Energy too high for ct-passing'
       else if (ee .lt. eelost) then
          write(*,*) 'Ct-passing lost'
       else
          write(*,*) 'Counter-passing'
       end if
    end if

    if (otype .ne. -100) then
       write(*,*) 'Stag edge     :', eeedge/eunit/1000., 'keV'
       write(*,*) 'Lost          :', eelost/eunit/1000., 'keV'
       if (itype .eq. 1) then
          write(*,*) 'T/P bound (I) :', eetpbound/eunit/1000., 'keV'
       else if (itype .eq. 2) then
          write(*,*) 'T/P bound (II):', eetpbound/eunit/1000., 'keV'
       end if
    end if
  end subroutine printorbittype

end module io
