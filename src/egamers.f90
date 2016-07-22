! Energetic Geodesic Acoustic ModE Radial Structure code -EGAMERS
! Copyright Zhisong Qu 2016

program EGAMERS

  use paras_phy, only : eunit, ei, paras_phy_init
  use paras_num, only : norbitintsample, dtorbitn
!!$  use profile, only : psi1
  use profile
  use orbit, only : getorbit
  use nl
  use io
  use sintable
  use radial_grid
  use trap_matrix
  use trap_grid
  use eigen
  use orbit_classify, only : tpbound
  implicit none

  integer, parameter :: nmax = 2000
  type(matrix) :: mat3
  type(tgrid) :: tg1
  integer :: ierr, i1, i2

  real, dimension(:), allocatable :: rdata, thetadata
  real :: perioddata, ee, mub0, pphi, eetpbound

  write(*,*)
  write(*,*) '******** EGAM RADIAL STRUCTURE CODE (EGAMERS) ********'
  write(*,*) '               Version May 2016'
  write(*,*)
  
  call readnamelist()
  call paras_phy_init()
  call profile_init()
  
  if (imode .eq. 1) then
     ! imode == 1
     ! Run as an eigenvalue solver
     call rgrid_init(nradial_grid, 1, 0., 0., 0., 0.)
     call sintable_init(nmax, np_trap)

     call tmatrix_init(tm, ngtrap_mub0, trap_mub0start*1000.*eunit, &
          trap_mub0end*1000.*eunit, ngtrap_pphi, ngtrap_energyn, &
          ngtrap_energyb, trap_ebend, np_trap, ierr)
     call tmatrix_calculate(tm)
     
     call newton_init(omegain(1)**2)
     call newton_step()

     call printfield(v)
     call plotcontinuum()

     call newton_cleanup

     call tmatrix_destroy(tm)
     call rgrid_destroy()
     call sintable_destroy()

  else if (imode .eq. 2) then
     ! imode == 2
     ! Calculate the bounce frequency map for a given muB0

     write(*,*) 'Frequency map mode'
     write(*,*) 'muB0 = ', mub0_in, 'keV'
     mub0 = mub0_in * 1000. * eunit

     call rgrid_init(nradial_grid, 1, 0., 0., 0., 0.)
     call sintable_init(nmax, np_trap)

     ! write header of the output file
     call write_map_header()

     ! calculate frequency for trap grid
     call tgrid_init(tg1, mub0, ngtrap_pphi, ngtrap_energyn, &
          ngtrap_energyb, trap_ebend, np_trap)
     call tgrid_calculate(tg1)
     ! output trap grid frequency
     call plot_tgrid_map(tg1)
     call tgrid_destroy(tg1)

     call close_map()
     call rgrid_destroy()
     call sintable_destroy()

  else if (imode .eq. 3) then
     ! imode == 3
     ! Calculate the orbit and its bounce frequency
     ! for a given (E, muB0, Pphi)

     write(*,*) 'Orbit mode'

     mub0 = mub0_in * 1000. * eunit
     pphi = pphi_in * psi1 * ei

     if (ienergy_tpbound .eq. 0) then
        ee = energy_in * 1000. * eunit
     else if (ienergy_tpbound .eq. -1 .or. ienergy_tpbound .eq. 1) then
        if (pphi_in .ge. -1. .and. pphi_in .le. 0.) then
           write(*,*) 'Orbit close to t/p boundary'
           eetpbound = tpbound(mub0, pphi, ierr)
           
           if (ierr .le. 0) then
              write(*,*) 'Finding t/p boundary failed'
              stop
           end if
        else
           write(*,*) 'pphi_in must be between -1 and 0'
           stop
        end if
        if (ienergy_tpbound .eq. -1) then
           ee = eetpbound - mub0 * exp(-energy_log_in)
        else
           ee = eetpbound + mub0 * exp(-energy_log_in)
        end if
     else
        write(*,*) 'ienergy_bound must be -1, 0 or 1'
        stop
     end if

     write(*,*)
     write(*,*) 'Energy       ', ee/1000./eunit, 'keV'
     write(*,*) 'muB0         ', mub0_in,   'keV'
     write(*,*) 'Pphi/(e psi1)', pphi_in
     write(*,*)
    
     allocate(rdata(norbitintsample))
     allocate(thetadata(norbitintsample))

     call printorbittype(ee, mub0, pphi, vsign_in)

     call getorbit(ee, mub0, pphi, norbitintsample, vsign_in, dtorbitn, &
          rdata, thetadata, perioddata, ierr)
     
     if (ierr .eq. 1) then
        call printorbit(norbitintsample, rdata, thetadata, perioddata)
        write(*,*)
        write(*,*) 'Omega bounce (rad/s)', 2. * pi / perioddata
        write(*,*) 'Orbit frequency (Hz)', 1. / perioddata 
     else if (ierr .eq. 2) then
        write(*,*) 'Orbit lost'
     else
        write(*,*) 'Orbit does not exist'
     end if

     deallocate(rdata)
     deallocate(thetadata)

  else
     ! wrong input
     write(*,*) 'Only imode = 1~3 are supported currently'
     write(*,*) 'Please refer to io.f90 for more information'
     write(*,*)
  end if

  write(*,*) '******** CODE FINISHED ********'

end program EGAMERS
