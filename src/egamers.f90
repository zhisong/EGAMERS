! Energetic Geodesic Acoustic ModE Radial Structure code -EGAMERS
! Copyright Zhisong Qu 2016

program EGAMERS

  use mpi
  use paras_phy, only : eunit, ei, paras_phy_init
  use paras_num, only : norbitintsample, dtorbitn
!!$  use profile, only : psi1
  use profile
  use distribution_fun, only: distribution_fun_init
  use orbit, only : getorbit
  use nl
  use io
  use sintable
  use radial_grid
  use trap_matrix
  use trap_grid
  use continuum
  use eigen
  use pic
  use diagnostics
  use test_particles
  use orbit_classify, only : tpbound
  implicit none
  
  integer, parameter :: nmax = 2000
  type(matrix) :: mat3
  type(tgrid) :: tg1
  integer :: ierr, i1, i2, n_active, nt, nr

  real, dimension(:), allocatable :: rdata, thetadata
  real :: perioddata, ee, mub0, pphi, eetpbound

  call mpi_start()
  
  if (mpi_is_master()) then
    write(*,*)
    write(*,*) '******** EGAM RADIAL STRUCTURE CODE (EGAMERS) ********'
    write(*,*) '               Version May 2016'
#ifdef MPI
    write(*,*) '             MPI Version March 2018'
    write(*,*) '            PIC Version October 2018'
    write(*,*) '******************************************************'
    write(*,*) '   Number of cores used : ', mpi_get_ncpus()
#else
    write(*,*) '              Serialized version'
    write(*,*) '******************************************************'
#endif
  endif
  
  call readnamelist()
  call paras_phy_init()
  call profile_init()
  call distribution_fun_init()
  
  ! we want the title to show completely first
  call mpi_sync()

  if (imode .eq. 0) then
    ! imode == 0
    ! Run PIC simulation
    if (mpi_is_master()) write(*,*) 'PIC simulation mode'

    if (MOD(ngtrap_mub0, mpi_get_ncpus()) .ne. 0 ) then
      if (mpi_is_master()) then 
        write(*,*)
        write(*,*) '######################################################'
        write(*,*) 'WARNING: For good load balance, ngtrap_mub0 needs to be divisible by ncpus!'
        write(*,*) '######################################################'
        write(*,*)
      end if
    end if

    call rgrid_init(nradial_grid, igrid_type, xr1, sig1, xr2, sig2)
    call sintable_init(nmax, np_trap)

    call tmatrix_init(tm, ngtrap_mub0, trap_mub0start*1000.*eunit, &
                      trap_mub0end*1000.*eunit, ngtrap_pphi, ngtrap_energyn, &
                      ngtrap_energyb, trap_ebend, np_trap, ipphi_eqdistant, .false., ierr)
    call tmatrix_calculate(tm)

    ! we'd better finish the previous steps before starting simulation
    call mpi_sync()
    
    if (mpi_is_master()) write(*,*) 'Initializing PIC simulation...'

    ! now we need to initialize the pic simulation
    call pic_init()
    if (mpi_is_master()) call io_snapshot_init()
    call mpi_sync()

    ! initial output
    if (mpi_is_master()) then
      write(*,*) 'Starting PIC simulation...'
      write(*,*) 'ksteps = ', 0
    end if
    if (mpi_is_master() .and. nsnapfield > 0) call io_snapshot_field(efield)
    !if (mpi_is_master() .and. nsnappart > 0) call io_snapshot_part()

    ! run the simulation for ksteps
    do i1 = 1, ksteps
      ! don't run for the last time step
      call pic_step()
      ! mpi barrier after each time step
      call mpi_sync()
      ! output informations at particular time steps
      ! the number of the still active trapped particles
      if (MOD(i1, nscreen).eq.0 .or. i1.eq.ksteps) n_active = diag_active_trap_particles()
      if (mpi_is_master()) then
        ! screen output for every nscreen steps
        if (MOD(i1, nscreen).eq.0 .or. i1.eq.ksteps) then
          write(*,1000) i1, t*1000.0, n_active
1000 format('steps = ', i8, ', t=', e15.3,'ms, active markers ', i8)
        endif
        ! output the field for every nsnapfield steps
        if (MOD(i1, nsnapfield).eq.0 .or. i1.eq.ksteps) call io_snapshot_field(efield)
        ! output particles for every nsnappart steps (will overwrite)
        !if (MOD((i1), nsnappart).eq.0 .or. i1.eq.ksteps-1) call io_snapshot_part()
      end if
    end do

    if (mpi_is_master()) call io_snapshot_destroy()
    call pic_destroy()
    call tmatrix_destroy(tm, .false.)
    call rgrid_destroy()
    call sintable_destroy()

  else if (imode .eq. 1) then
    ! imode == 1
    ! Run as an eigenvalue solver
    if (mpi_is_master()) write(*,*) 'Eigenvalue solver mode'

    call rgrid_init(nradial_grid, igrid_type, xr1, sig1, xr2, sig2)
    call sintable_init(nmax, np_trap)
    
    call tmatrix_init(tm, ngtrap_mub0, trap_mub0start*1000.*eunit, &
                      trap_mub0end*1000.*eunit, ngtrap_pphi, ngtrap_energyn, &
                      ngtrap_energyb, trap_ebend, np_trap, ipphi_eqdistant, .false., ierr)
    call tmatrix_calculate(tm)

    ! we'd better finish the previous step before iterating frequency
    call mpi_sync()

    ! calculate the continuum first
    if (mpi_is_master()) write(*,*) 'Computing EGAM continuum frequency'
    call continuum_init(ngtrap_pphi)
    call continuum_compute(omegain(1)**2)

    ! calculate the global mode next
    if (mpi_is_master()) write(*,*) 'Computing EGAM global frequency'
    call newton_init(omegain(1)**2)
    call newton_step()

    if (mpi_is_master()) then
      call printfield(v)
      call plotcontinuum()
    endif
     
    call newton_cleanup

    call continuum_destroy()

    call tmatrix_destroy(tm, .false.)
    call rgrid_destroy()
    call sintable_destroy()

  else if (imode .eq. 2) then
    ! imode == 2
    ! Calculate the bounce frequency map for a given muB0

    ! we don't need mpi for frequency map since it is quite fast

    if (mpi_is_master()) then
     
      write(*,*) 'Frequency map mode'
      write(*,*) 'muB0 = ', mub0_in, 'keV'
      mub0 = mub0_in * 1000. * eunit

      call rgrid_init(nradial_grid, igrid_type, xr1, sig1, xr2, sig2)
      call sintable_init(nmax, np_trap)

      ! write header of the output file
      call write_map_header()

      ! calculate frequency for trap grid
      call tgrid_init(tg1, mub0, ngtrap_pphi, ngtrap_energyn, &
                      ngtrap_energyb, trap_ebend, np_trap, ipphi_eqdistant)
      call tgrid_calculate(tg1)
      ! output trap grid frequency
      call plot_tgrid_map(tg1)
      call tgrid_destroy(tg1)

      call close_map()
      call rgrid_destroy()
      call sintable_destroy()
    endif

  else if (imode .eq. 3) then
    ! imode == 3
    ! Calculate the orbit and its bounce frequency
    ! for a given (E, muB0, Pphi)

    ! we don't need mpi for running a single orbit

    if (mpi_is_master()) then
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
    endif
  
  else if (imode .eq. 4) then
  ! test particle mode
#ifndef NC
    if (mpi_is_master()) then
      write(*,*) "Currently we only support field snapshot in NETCDF format. Please recompile with NETCDF" 
    end if
    stop "NETCDF needed for test particle mode"
#endif 

    if (mpi_is_master()) then
      write(*,*) 'Test particle mode'
      write(*,*) 'mub0 = ', mub0_test, 'keV'
      write(*,'(A26)', ADVANCE='no') 'Reading field snapshot...'
    end if
      
    call rgrid_init(nradial_grid, igrid_type, xr1, sig1, xr2, sig2)    

    if (mpi_is_master()) then
      call io_read_field_init(nt, nr, rgrid)
      if (nr .ne. nradial_grid*2-2) then 
        stop "radial grid from namelist and field data does not match"
      end if

      ! allocate the space
      allocate(t_list(nt))
      allocate(lambda_t(nr, nt))
      allocate(eta_t(nr, nt))

      do i1 = 1, nt
        call io_read_field(t_list(i1), lambda_t(:,i1), eta_t(:,i1), i1)
      end do
 
      call io_read_field_destroy()

      write(*,*) 'Successful'
      write(*,*) 'Containing: ', nt, 'steps, total ', t_list(nt)*1000, 'ms'
      write(*,*) 'Calculating orbits...'
    end if

    ! now bcast the grid read from file to all cpus
    call mpi_bcast_array(rgrid, 0)

    call mpi_sync()

    call sintable_init(nmax, np_trap)

    mub0 = mub0_test * 1000. * eunit    
    ! calculate frequency for trap grid
    call tgrid_init(tg_test, mub0, ngtrap_pphi, ngtrap_energyn, &
                    ngtrap_energyb, trap_ebend, np_trap, ipphi_eqdistant)
    call tgrid_calculate(tg_test)

    call mpi_sync()

    if (mpi_is_master()) write(*,*) 'Initialing Particles...'

    call test_pic_init()
    call io_snapshot_test_particles_init()
    
    if (mpi_is_master()) write(*,*) 'steps =', 0

    i2 = -1
    do i1 = 1, ksteps_test
      call test_pic_step()
      if (mpi_is_master() .and. MOD(i1, nscreen_test).eq.0) then
        write(*,*) 'steps =', i1
      end if
      call mpi_sync()
      if (i1 .ge. ksteps_snapstart) i2 = i2 + 1
      if (i2 .ge. 0 .and. (MOD(i2, nsnappart_test).eq.0 .or. i1.eq.ksteps)) call io_snapshot_test_particles()
      call mpi_sync()
    end do

    ! cleaning up
    call io_snapshot_test_particles_destroy()
    call test_pic_destroy()
    call tgrid_destroy(tg_test)

    call rgrid_destroy()
    call sintable_destroy()
    
    if (mpi_is_master()) then
      if (ALLOCATED(t_list)) deallocate(t_list)
      if (ALLOCATED(lambda_t)) deallocate(lambda_t)
      if (ALLOCATED(eta_t)) deallocate(eta_t)
    end if

  else
    ! wrong input
    if (mpi_is_master()) then
      write(*,*) 'Only imode = 1~4 are supported currently'
      write(*,*) 'Please refer to io.f90 for more information'
      write(*,*)
    endif
  end if

  ! MPI END
  call mpi_end()
  if (mpi_is_master()) then
    write(*,*) '******** CODE FINISHED ********'
  end if
end program EGAMERS
