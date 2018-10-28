! Input and output

module io

! use netCDF for output if defined
#ifdef NC
  use netcdf
#endif

  use paras_num, only : neigenmax
  use nl

  implicit none

  integer, parameter :: iofield    = 20  ! electric field file
  integer, parameter :: iofqc      = 21  ! continuum and global mode fqc file
  integer, parameter :: iomap      = 22
  integer, parameter :: ioorbit    = 23  ! orbit file

#ifdef NC
  integer, parameter :: NFIELD_DIMS = 2
  integer, private :: ncid_field, ncid_part
  integer, private :: nr_dimid, rec_dimid, nr_varid, rec_varid, eta_varid, lambda_varid
  integer, private :: irec, NR
#else
  integer, parameter :: iosnapfield= 24  ! field snapshot file
  integer, parameter :: iosnappart = 25  ! particle snapshot file
#endif

contains

! ////// OUTPUT //////

  subroutine io_snapshot_init()
    ! initialize the snapshot files
    use radial_grid

#ifdef NC
    character (len = *), parameter :: FIELD_FILE = "snapshot.field.nc"
    character (len = *), parameter :: PART_FILE = "snapshot.part.nc"
    character (len = *), parameter :: NR_NAME = "radial_element"
    character (len = *), parameter :: REC_NAME = "time"
    character (len = *), parameter :: ETA_NAME = "eta"
    character (len = *), parameter :: LAMBDA_NAME = "lambda"
    character (len = *), parameter :: LFIELDOUTPUT_NAME = "lfieldoutput"

    ! For the field, we are writing a 2D data, in grid of nele and t
    integer :: dimids(NFIELD_DIMS)

    integer :: i1
    real, dimension(:), allocatable :: r
    real :: dx

    if (lfieldoutput .eq. 0) then
      NR = nele
    elseif (lfieldoutput .eq. 1) then
      NR = nfieldoutput
    else
      stop "lfield output invalid"
    end if

    ! Create the file.
    call check( nf90_create(FIELD_FILE, nf90_clobber, ncid_field) )
    call check( nf90_create(PART_FILE, nf90_clobber, ncid_part) )

    ! Define the dimensions. The record dimension is defined to have
    ! unlimited length - it can grow as needed.
    call check( nf90_def_dim(ncid_field, NR_NAME, NR, nr_dimid) )
    call check( nf90_def_dim(ncid_field, REC_NAME, NF90_UNLIMITED, rec_dimid) )

    ! Define variable
    call check( nf90_def_var(ncid_field, REC_NAME, NF90_REAL8, rec_dimid, rec_varid) )

    dimids = (/ nr_dimid, rec_dimid /)

    call check( nf90_def_var(ncid_field, LAMBDA_NAME, NF90_REAL8, dimids, lambda_varid) )
    call check( nf90_def_var(ncid_field, ETA_NAME, NF90_REAL8, dimids, eta_varid) )
    
    call check( nf90_def_var(ncid_field, NR_NAME, NF90_REAL8, nr_dimid, nr_varid) )
    call check( nf90_put_att(ncid_field, nr_varid, LFIELDOUTPUT_NAME, lfieldoutput) )
    call check( nf90_put_att(ncid_field, nr_varid, "igrid_type", igrid_type) )
    call check( nf90_put_att(ncid_field, nr_varid, "XR1", xr1) )
    call check( nf90_put_att(ncid_field, nr_varid, "XR2", xr2) )
    call check( nf90_put_att(ncid_field, nr_varid, "SIG1", sig1) )
    call check( nf90_put_att(ncid_field, nr_varid, "SIG2", sig2) )

      
    call check( nf90_enddef(ncid_field) )

    ! put in the nr variable
    allocate(r(NR))
    if (lfieldoutput .eq. 0) then   
      r(1) = 0.0
      r(NR) = 1.0
      do i1 = 2, nelement-1
        r(i1*2-2) = rgrid(i1)
        r(i1*2-1) = rgrid(i1)
      enddo
    else
      dx = 1. / real(nfieldoutput - 1)
      do i1 = 1, 100
        r(i1) = dx * float(i1-1)
      enddo
    endif

    call check( nf90_put_var(ncid_field, nr_varid, r) )
    if (ALLOCATED(r)) deallocate(r)

    ! record counter = 0
    irec = 0

#else
    if (lfieldoutput .lt. 0 .or. lfieldoutput .gt. 1) then
      stop "lfield output invalid"
    endif
    open(UNIT=iosnapfield, FILE='snapshot.field.out',ACTION='WRITE',FORM = 'unformatted')
    open(UNIT=iosnappart, FILE='snapshot.part.out',ACTION='WRITE',FORM = 'unformatted')

    ! write the header
    write(iosnapfield) lfieldoutput, nele, nfieldoutput
#endif

  end subroutine io_snapshot_init

  subroutine io_snapshot_destroy()
    ! initialize the snapshot files
#ifdef NC
    call check( nf90_close(ncid_field))
    call check( nf90_close(ncid_part))
#else
    close(UNIT=iosnapfield)
    close(UNIT=iosnappart)
#endif
  end subroutine

  subroutine io_snapshot_field(ev)
    ! write the field to file
    use field
    use pic, only : t
    use radial_grid
    use hermite
    implicit none
    type(field_vector) :: ev

    real, dimension(:), allocatable :: er
    real :: dx, x
    complex :: value
    integer :: i1, i2

#ifdef NC
    integer :: start(NFIELD_DIMS), counts(NFIELD_DIMS)
  
    ! increase the record counter
    irec = irec + 1

    start = (/1, irec/)
    counts = (/NR, 1/)

    ! first write time variable
    call check( nf90_put_var(ncid_field, rec_varid, (/t/), start = (/irec/), &
                              count = (/1/) ) ) 
#endif

    if (lfieldoutput .eq. 0) then
#ifdef NC
      call check( nf90_put_var(ncid_field, lambda_varid, ev%lambda, start = start, &
                              count = counts) )      
      call check( nf90_put_var(ncid_field, eta_varid, ev%eta, start = start, &
                              count = counts) ) 
#else
      write(iosnapfield) ev%lambda(1:nele)
      write(iosnapfield) ev%eta(1:nele)
#endif
    else
      allocate(er(nfieldoutput))
      er(:) = 0.0
      dx = 1. / real(nfieldoutput - 1)

      do i1 = 1, nfieldoutput
        x = dx * real(i1 - 1)
        ! first grid point
        er(i1) = ev%lambda(1) * c2(x, 0., -1., rgrid(1))
        ! last grid point
        er(i1) = er(i1) + ev%lambda(nele) * c2(x, 1., rgrid(nelement-1), 1.2)

        ! adding contribution from other grids
        do i2 = 2, nelement - 1
            er(i1) = er(i1) + ev%lambda(i2*2-2) * c1(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
            er(i1) = er(i1) + ev%lambda(i2*2-1) * c2(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
        end do
      end do
#ifdef NC
      call check( nf90_put_var(ncid_field, lambda_varid, er, start = start, &
                              count = counts) )
#else
      write(iosnapfield) er
#endif

      er(:) = 0.0
      dx = 1. / real(nfieldoutput - 1)

      do i1 = 1, nfieldoutput
        x = dx * real(i1 - 1)
        ! first grid point
        er(i1) = ev%eta(1) * c2(x, 0., -1., rgrid(1))
        ! last grid point
        er(i1) = er(i1) + ev%eta(nele) * c2(x, 1., rgrid(nelement-1), 1.2)

        ! adding contribution from other grids
        do i2 = 2, nelement - 1
            er(i1) = er(i1) + ev%eta(i2*2-2) * c1(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
            er(i1) = er(i1) + ev%eta(i2*2-1) * c2(x, rgrid(i2), rgrid(i2-1), rgrid(i2+1))
        end do
      end do
#ifdef NC
      call check( nf90_put_var(ncid_field, eta_varid, er, start = start, &
                              count = counts) )
#else
      write(iosnapfield) er
#endif
      if(ALLOCATED(er)) deallocate(er)

    end if

#ifdef NC
    call check( nf90_sync(ncid_field) )
#else
    FLUSH(iosnapfield)
#endif

  end subroutine io_snapshot_field

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
    integer :: i1, i2, ndee
    real :: dee, ee

    write(iomap,*) this%npphin, this%neen

    do i1 = 1, this%npphin
       do i2 = 1, this%periodn(i1)%n
          write(iomap,*) this%pphigrid(i1)/ei/psi1, &
               this%periodn(i1)%x(i2)/eunit/1000., 2*pi/this%periodn(i1)%y(i2)
       end do
       if (this%periodn(i1)%n .ne. this%periodn(1)%n) then
          ndee = this%periodn(1)%n - this%periodn(i1)%n
          do i2 = 1, ndee
             write(iomap,*) this%pphigrid(i1)/ei/psi1, &
                  this%periodn(i1)%x(this%periodn(i1)%n) /eunit/1000.,  2*pi/this%periodn(i1)%y(this%periodn(i1)%n)
          end do
       end if
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

! ///// Subroutine only for netcdf
#ifdef NC
    subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped, NETCDF error"
    end if
  end subroutine check
#endif

end module io
