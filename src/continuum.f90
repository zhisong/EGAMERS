! module for calculating the EGAM continuum using local approximation
module continuum

  ! computed egam continuum frequency and growth rate (imag part) on radial grid
  complex, allocatable, dimension(:) :: egam_continuum
  ! the radial grid location
  real, allocatable, dimension(:) :: egam_continuum_r
  ! the cooresponding ipphi for the radial grid
  integer, allocatable, dimension(:) :: egam_continuum_ipphi
  
  ! number of radial data points
  integer :: egam_continuum_nr

contains
  subroutine continuum_init(npphi)
    use paras_phy
    use profile
    use trap_grid, only : sgrid
    implicit none

    integer, intent(in) :: npphi

    real :: ds, s
    integer :: i1

    egam_continuum_nr = npphi
    allocate(egam_continuum(npphi))
    allocate(egam_continuum_r(npphi))
    allocate(egam_continuum_ipphi(npphi))

    ds = 1.0 / real(1 + npphi)
    egam_continuum = 0.0
    do i1 = 1, npphi
      egam_continuum_ipphi(i1) = i1
      egam_continuum_r(i1) = psitor(-sgrid(real(i1) * ds) / ei)
    end do

  end subroutine continuum_init

  subroutine continuum_destroy()
    implicit none

    if (allocated(egam_continuum)) deallocate(egam_continuum)
    if (allocated(egam_continuum_r)) deallocate(egam_continuum_r)
    if (allocated(egam_continuum_ipphi)) deallocate(egam_continuum_ipphi)

  end subroutine continuum_destroy

  subroutine continuum_compute(lambdain)
    use paras_num
    use mpi
    use profile
    use trap_grid
    use trap_matrix
    use eigen, only : tm
    use mhd_matrix

    implicit none
    
    complex, intent(in) :: lambdain

    integer :: i1, i2, ipphi
    real :: r, dlambda
    complex :: local1, local2, local3trap, local3trap1, local, localprime, omega, omega1
    complex :: lambda, lambda1, lambda_move, lambda_guess

    do i1 = 1, egam_continuum_nr
      ipphi = egam_continuum_ipphi(i1)
      r = egam_continuum_r(i1)
      lambda_guess = lambdain
      lambda = lambda_guess
      dlambda = 1.e20

      call getlocal1(r, local1)
      call getlocal2(r, local2)

      do i2 = 1, nmaxit

        ! we need to sync before starting
        ! (which is actually trivial because of bcast)
        call mpi_sync()

        ! first, we need to broadcast lambda
        call mpi_bcast_scalar(lambda, 0)
       
        ! then broadcast the stop condition
        call mpi_bcast_scalar(dlambda, 0)

        ! if precision is reached
        if (dlambda .le. erreig) then
          exit
        end if

        lambda_move = real(dlambdapctg * lambda)
        omega  = sqrt(lambda)
        omega1 = sqrt(lambda + lambda_move)
        !write(*,*) 'before getlocal3trap is ok'
        call getlocal3trap(tm, ipphi, r, omega, local3trap)
        call getlocal3trap(tm, ipphi, r, omega1, local3trap1)
        !write(*,*) 'after getlocal3trap'
        ! we only do it in master cpu
        if (mpi_is_master()) then
          local = local2 + local3trap + (-lambda) * local1
          localprime =  (1./lambda_move) * local3trap1 + (-1./lambda_move) * local3trap &
            + (-1.,0.) * local1
        
          lambda1 = lambda - local / localprime
          dlambda = abs(lambda1 - lambda) / abs(lambda1)
          lambda = lambda1
          !write(*,*) r,i2,sqrt(lambda), local3trap
        endif

      end do
      egam_continuum(i1) = sqrt(lambda)

    end do


  end subroutine continuum_compute


end module continuum