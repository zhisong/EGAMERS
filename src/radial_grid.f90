! generate radial grid points of Hermite elements

module radial_grid

  use paras_num
  use matrix_module
  implicit none

  public
  real, allocatable, dimension(:), public  :: rgrid
  integer, public :: nelement, nele
  
  contains

    subroutine rgrid_init(nelementin, ntype, xr1, sig1, xr2, sig2)
      ! initialize radial grids
      ! INPUTS : number of elements, radial grid accumulation point and width
      !          type : 1 (equal-distance), 2 (gaussian)
      !          only type 1 is implemented currently
      
      implicit none
      
      integer, intent(in) :: nelementin, ntype
      real, intent(in) :: xr1, sig1, xr2, sig2

      real :: drstep
      integer :: i1
      
      if ((nelementin .le. 2) .or. (nelementin .gt. nmaxelement)) then
         nelement = 30  ! Set to default
         write(*,*) 'Number of radial grid point must be > 2 and <', nmaxelement
         write(*,*) 'Set to default (30)'
      else
         nelement = nelementin
         nele = 2 * nelement - 2
      end if
      
      ! destroy if already allocated
      call rgrid_destroy()

      allocate(rgrid(nelementin))

      if (ntype .eq. 1) then
         drstep = 1. / real(nelement-1)
         do i1 = 1, nelement
            rgrid(i1) = drstep * real(i1 - 1)
         end do
      else if (ntype .eq. 2) then
        call meshac(xr1, xr2, 0.3, sig1, sig2, 1.0, 0.0, 1.0, rgrid, nelement, nelement-1)
      end if
    end subroutine rgrid_init
   
    subroutine rgrid_destroy()
      ! destroy radial grid
      implicit none
      if (allocated(rgrid)) deallocate(rgrid)
    end subroutine rgrid_destroy


  SUBROUTINE MESHAC(XR1,XR2,BGF,SIG1,SIG2,FACT,RS0,RSA,SGRID,NG,NGINT)
!-----------------------------------------------------------------------
! MESH ACCUMULATION IN XR1 AND XR2
!-----------------------------------------------------------------------
    REAL, INTENT(IN) :: XR1,XR2,BGF,SIG1,SIG2,FACT,RS0,RSA
    INTEGER, INTENT(IN) :: NG, NGINT
    REAL, DIMENSION(:) :: SGRID

    INTEGER :: JINT, J, I, J1INT
    REAL :: ZS, ZSUM, ZDS, ZF, ZNORM, ZI, ZDS0, ZWL, ZFD, ZSUM1
!     CONDITIONS

    IF(BGF .EQ.0.) GO TO 60
    IF(XR1 .EQ.0.) GO TO 60
    IF(XR2 .EQ.0.) GO TO 60
    IF(SIG1.EQ.0.) GO TO 60
    IF(SIG2.EQ.0.) GO TO 60
    IF(FACT.EQ.0.) GO TO 60
!------------------------------------------- EVALUATION OF NORM
    JINT = 100*NGINT+1
    ZS   = RS0
    ZSUM = 0.0
    ZDS  = (RSA - RS0) / FLOAT(JINT-1)
    DO 10 J=1,JINT
        ZF   = FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT)
        ZSUM = ZSUM + ZF * ZDS
        ZS   = ZS + ZDS
10 CONTINUE
    ZNORM   = (RSA-RS0) / (ZSUM)

    J1INT = 100
    JINT  = J1INT * NGINT
    ZFD   = (RSA - RS0) / FLOAT(NGINT)
    ZS    = RS0
    ZSUM  = 0.0
    ZF    = FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT)
    ZF    = ZF * ZNORM
    ZDS0  = (RSA - RS0) * ZF / FLOAT(JINT)
    I = 2
20 CONTINUE
    ZI    = FLOAT(I-1) * ZFD 
    ZDS   = ZDS0 / ZF
    ZS    = ZS + ZDS
    ZSUM1 = ZSUM
    ZF    = FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT)
    ZF    = ZF * ZNORM
    ZSUM  = ZSUM + ZF * ZDS
    IF(ZI.GT.ZSUM) GO TO 20
    ZWL   = (ZI-ZSUM1)/(ZSUM - ZSUM1)
    SGRID(I) = ZS - ZDS * (1.0 - ZWL)
    IF(SGRID(I).LT.SGRID(I-1)) GOTO 40
    I = I + 1
    IF(I.GT.NGINT) GO TO 30
    GO TO 20
30 CONTINUE
    SGRID(1) = RS0
    SGRID(NGINT+1) = RSA
    SGRID(NGINT) = 0.5 * (SGRID(NGINT-1) + SGRID(NGINT+1))
    IF(SGRID(NGINT-1).LT.RSA) GO TO 50
    WRITE(*,31)
40 WRITE(*,41)
    STOP
50 CONTINUE
60 CONTINUE

    RETURN

31 FORMAT('0',' ERR. IN S.R. MESHAC2 : SGRID(NGINT) GT. RSA ')
41 FORMAT('0',' ERR. IN S.R. MESHAC2 : SGRID(I) .GT. SGRID(I+1) ')
51 FORMAT('1',61('*'),' *',20X,'MESHAC GRID POINTS',21X,'*',' ',61('*'))
53 FORMAT(15F8.4)
    END

    FUNCTION FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT)
!-----------------------------------------------------------------------
!     BGF + (1 - BGF) * (GAUSS1 + FACT * GAUSS2) / FACT
!-----------------------------------------------------------------------
    REAL, INTENT(IN) :: ZS, BGF, XR1, XR2, SIG1, SIG2, FACT
    REAL :: FGAUS
    REAL :: ZNORM1, ZNORM2, ZEX1, ZEX2, F1, F2

    ZNORM1 = 0.39894 / SIG1
    ZNORM2 = 0.39894 / SIG2
    ZEX1   = -0.5 * (ZS - XR1)**2 / SIG1**2
    ZEX2   = -0.5 * (ZS - XR2)**2 / SIG2**2
    F1     = ZNORM1 * EXP(ZEX1)
    F2     = ZNORM2 * EXP(ZEX2)
    FGAUS  = BGF + (1.0 - BGF) * (F1 + FACT * F2) / FACT
    RETURN
    END

end module radial_grid
