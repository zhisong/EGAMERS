! The cubic spline interpolating type and subroutines


module spline_module
!
! TYPE:
! spline - store pre-computed cubic spline interpolation of given x and y
!
! SUBROUTINES:
! (public)
! spline_init    - allocate arrays for type(spline), must call first
! spline_build   - build the spline AFTER x and y are filled in
!                  (an object interface for type spline to call build_spline)
! spline_find    - find a list of x for a given y value in the spline
! spline_findc   - find a list of x for a given y value in the spline, y is complex
! spline_destroy - deallocate the arrays, must call before program ends
! (private)
! build_spline   - do the real job to build spline AFTER x and y are filled in
!                  (copied from MISHKA code, original in GERMEN)
!
! FUNCTIONS:
! (public)
! spline_interp  - real, intepolation, can only be used after spline_init and 
!                  spline_build are called
!                  (an object interface for type spline to call spwert)
!                  y(x), y'(x), y''(x) and y'''(x)  
! spline_interp1 - real, intepolation, y(x) only
! spline_interpd1- real, intepolation, y'(x) only  
! spline_interpc - complex, intepolation
!
! (private)
! spwert         - real, do the real job to intepolate
!                  (copied from MISHKA code, original in GERMEN)
! spwert1        - real, y(x)
! spwertd1       - real, y'(x)
  use cubic
  implicit none

  private
  public :: spline_init, spline_destroy, spline_interp, spline_interpc,&
       spline_build, spline_find, spline_interp1, spline_interpd1

  type, public :: spline
     ! the type contains pre-computed cubic splines
     
     ! spline data
     real, allocatable, dimension(:), public :: a, b, c, d
     ! original data
     real, allocatable, dimension(:), public :: x, y
     
     integer, public :: n

     logical :: iequaldistant

  end type spline

contains

  subroutine spline_init(this ,n)
    ! allocate the arrays and an interface to build_spline
    ! INPUT : n, number of data points
    implicit none
    
    type(spline) :: this
    integer, intent(in) :: n

!!$    if (allocated(this%x)) call spline_destroy(this)
    this%n = n
    allocate(this%x(n))
    allocate(this%y(n))
    allocate(this%a(n))
    allocate(this%b(n))
    allocate(this%c(n))
    allocate(this%d(n))
    this%iequaldistant = .false.
    
  end subroutine spline_init

  subroutine spline_destroy(this)
    ! deallocate the spline object
    implicit none
    
    type(spline) :: this
    
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%y)) deallocate(this%y)
    if (allocated(this%a)) deallocate(this%a)
    if (allocated(this%b)) deallocate(this%b)
    if (allocated(this%c)) deallocate(this%c)
    if (allocated(this%d)) deallocate(this%d)

  end subroutine spline_destroy

  subroutine spline_build(this, alfa, beta, typ, istart, iend)
    ! build the spline objec from x and y
    ! interface to call BUILD SPLNE
    ! INPUT : typ - end condition
    !   0 - not-a-knot
    !   1 - specify 1st derivative in alfa and beta
    !   2 -    "    2nd    "
    !   3 -    "    3rd    "
    type(spline) :: this
    integer, intent(in) :: typ, istart, iend
    real, intent(in) :: alfa, beta
    
    integer :: n, i1, istart1, iend1
    n = iend - istart + 1
    if (n .ge. 3) then
       istart1 = istart
       iend1 = iend
       call BUILD_SPLINE(n, this%x(istart), this%y(istart), alfa, beta, typ,&
            this%a(istart), this%b(istart), this%c(istart), this%d(istart))
    else if(n .ge. 1) then
       istart1 = istart - 1
       iend1 = iend + 1
       n = n + 2
       if (istart1 .lt. 1) then
          istart1 = 1
          n = n - 1
       end if
       if (iend1 .gt. this%n) then
          iend1 = this%n
          n = n  - 1
       end if
       if (n .ge. 3) then
          call BUILD_SPLINE(n, this%x(istart1), this%y(istart1), alfa, beta,&
 typ, this%a(istart1), this%b(istart1), this%c(istart1), this%d(istart1))
       else
          iend1 = istart - 1
       end if
    else
       istart1 = 1
       iend1 = 1
    end if
       
       
    do i1 = 1, istart1 - 1
       this%a(i1) = 0.
       this%b(i1) = 0.
       this%c(i1) = 0.
       this%d(i1) = 0.
    end do

    do i1 = iend1+1, this%n
       this%a(i1) = 0.
       this%b(i1) = 0.
       this%c(i1) = 0.
       this%d(i1) = 0.
    end do

  end subroutine spline_build

  subroutine spline_find(this, yin, istart, iend, nfound, xlist) 
    ! find a list of x for a given y
    ! INPUT  : target y, start node, end node
    ! OUTPUT : number of roots found, a list of roots, where they are

    use paras_num
    implicit none

    type(spline), intent(in) :: this
    real, intent(in) :: yin
    integer, intent(in) :: istart, iend
    integer, intent(out) :: nfound
    real, dimension(nmaxroot), intent(out) :: xlist

    integer :: i1, i2, nc
    real :: a, b, c, p1, p2
    real :: xl, xr, tmp
    real, dimension(3) :: xc

    nfound = 0

    if ((istart .le. 0) .or. (istart .ge. iend) .or. (iend .gt. this%n)) then
       ! illegal interval arguments
       write(*,*) 'illegal interval arguments in spline_find', istart, iend
       return
    end if

    do i1 = istart, iend-1
       ! left and right limit
       xl = this%x(i1)
       xr = this%x(i1+1)
       ! the first derivatives on the nodes
       p1 = this%b(i1)
       if (i1 .lt. this%n-1) then
          p2 = this%b(i1+1)
       else
          tmp = spline_interp(this, xr, xc)
          p2 = xc(1)
       end if

       if (((this%y(i1) - yin) * (this%y(i1+1) - yin) .le. 0.) &
            .or. (p1 * p2 .le. 0.)) then
          ! 1. y(i1) and y(i1+1) are on different sides of yin
          !    there must be at least one root inside
          ! or
          ! 2. the first derivative changes sign
          !    there is an extreme inside, possible double roots
          !
          ! solve for x from a cubic equation, call the solver
          ! spline : y = a + b (x-xl) + c (x-xl)^2 + d (x-xl)^3
          !          0 <= x - xl <= xr - xl.
          ! need to convert to input of the solver
          a = this%c(i1) / this%d(i1)
          b = this%b(i1) / this%d(i1)
          c = (this%a(i1) - yin) / this%d(i1)
          call cubicroot(a, b, c, nc, xc)
          ! check if xleft <= x <= xright
          ! if yes, add to the output list
          do i2 = 1, nc
             if ((xc(i2) .ge. 0.) .and. (xc(i2) .lt. xr - xl)) then
                nfound = nfound + 1
                xlist(nfound) = xc(i2) + xl
             else if ((i1 .eq. iend-1) .and. (xc(i2) .eq. xr - xl)) then
                ! special treatment for the last node
                nfound = nfound + 1
                xlist(nfound) = xc(i2) + xl
             end if
          end do
       end if
    end do
  end subroutine spline_find


!!$  subroutine spline_findc(this, yin, istart, iend, nfound, xlist, iposlist) 
!!$    ! find a list of x for a given y, y is complex
!!$    ! INPUT  : target y, start node, end node
!!$    ! OUTPUT : number of roots found, a list of roots, where they are
!!$
!!$    use paras_num
!!$    implicit none
!!$
!!$    type(spline), intent(in) :: this
!!$    complex, intent(in) :: yin
!!$    integer, intent(in) :: istart, iend
!!$    integer, intent(out) :: nfound
!!$    complex, dimension(nmaxroot), intent(out) :: xlist
!!$    integer, dimension(nmaxroot), intent(out) :: iposlist
!!$
!!$    real :: distance
!!$    real, dimension(nmaxroot) :: xlistreal
!!$    complex, dimension(3) :: abltg
!!$    integer :: i1
!!$
!!$    allocate(distance(n))
!!$    do i1 = istart, iend
!!$       distance(i1) = (imag(yin))**2 + (real(yin)-y(istart))**2
!!$    end do
!!$
!!$    ! check 
!!$    deallocate(distance)
!!$
!!$  end subroutine spline_findc

  real function spline_interp(this, xin, abltg)
    ! interpolation to get y(xin)
    ! an interface to call SPWERT
    ! abltg(1 to 3) first, second and third order derivatives
    implicit none
    
    type(spline), intent(in) :: this
    real, intent(in) :: xin
    real, dimension(3), intent(out):: abltg
    !if (xin < this%x(1) .or. xin > this%x(this%n)) write(*,*) "ERROR!!!", xin, this%x(1), this%x(this%n)
    spline_interp = SPWERT(this%n, xin, this%a, this%b, this%c, this%d, this%x, abltg, this%iequaldistant)
    
  end function spline_interp

  real function spline_interp1(this, xin)
    ! interpolation to get y(xin)
    ! an interface to call SPWERT1
    implicit none
    
    type(spline), intent(in) :: this
    real, intent(in) :: xin
    
    spline_interp1 = SPWERT1(this%n, xin, this%a, this%b, this%c, this%d, this%x, this%iequaldistant)

  end function spline_interp1

  real function spline_interpd1(this, xin)
    ! interpolation to get y'(xin)
    ! an interface to call SPWERTD1
    implicit none
    
    type(spline), intent(in) :: this
    real, intent(in) :: xin
    
    spline_interpd1 = SPWERTD1(this%n, xin, this%a, this%b, this%c, this%d, this%x, this%iequaldistant)

  end function spline_interpd1
  
  complex function spline_interpc(this, xin, abltg)
    ! interpolation to get y(xin), xin is complex
    ! an interface to call SPWERT
    ! abltg(1 to 3) first, second and third order derivatives
    implicit none
    
    type(spline), intent(in) :: this
    complex, intent(in) :: xin
    complex, dimension(3), intent(out):: abltg
    
    real :: xreal, yreal
    complex :: ximag
    real, dimension(3) :: abltgreal
    
    xreal = real(xin)
    ximag = xin - xreal
    yreal = SPWERT(this%n, xreal, this%a, this%b, this%c, this%d, this%x, abltgreal, this%iequaldistant)
    spline_interpc = yreal + abltgreal(1) * ximag &
         + 1./2.*abltgreal(2) * ximag**2 + 1./6.*abltgreal(3) * ximag**3
    abltg(1) = abltgreal(1) + abltgreal(2) * ximag &
         + 1./2.*abltgreal(3)
    abltg(2) = abltgreal(2) + abltgreal(3) * ximag
    abltg(3) = abltgreal(3)

  end function spline_interpc

 ! ************** INTERNAL FUNCTIONS ****************

  SUBROUTINE BUILD_SPLINE(N,X,Y,ALFA,BETA,TYP,A,B,C,D)
    ! use the subroutine from MISHKA code
    ! A.B. Mikhailovskii etal. Plasma Physics Reports, 21, 844(1997)
!!$C-----------------------------------------------------------------------
!!$C     INPUT:
!!$C
!!$C     N     ANZAHL DER KNOTEN
!!$C     X     ARRAY DER X-WERTE
!!$C     Y     ARRAY DER Y-WERTE
!!$C     ALFA  RANDBEDINGUNG IN X(1)
!!$C     BETA        "       IN X(N)
!!$C     TYP   =  0  NOT-A-KNOT SPLINE
!!$C              1  ALFA, BETA 1. ABLEITUNGEN VORGEGEBEN
!!$C              2    "    "   2.     "           "
!!$C              3    "    "   3.     "           "
!!$C
!!$C     BEMERKUNG: MIT TYP = 2 UND ALFA = BETA = 0 ERHAELT MAN
!!$C           EINEN NATUERLICHEN SPLINE
!!$C
!!$C     OUTPUT:
!!$C
!!$C     A, B, C, D     ARRAYS DER SPLINEKOEFFIZIENTEN
!!$C       S = A(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2+ D(I)*(X-X(I))**3
!!$C
!!$C     BEI ANWENDUNGSFEHLERN WIRD DAS PROGRAMM MIT ENTSPRECHENDER
!!$C     FEHLERMELDUNG ABGEBROCHEN
!!$C-----------------------------------------------------------------------

    USE INTERFACES
    USE PARAS_NUM

    IMPLICIT NONE
    INTEGER  :: N, TYP
    REAL X(*), Y(*), A(*), B(*), C(*), D(*)
    REAL :: ALFA, BETA
    INTEGER ::  I, IERR
    REAL, DIMENSION(NMAXSPLINE) :: H

    IF((TYP.LT.0).OR.(TYP.GT.3)) THEN
       WRITE(*,*) 'FEHLER IN ROUTINE SPLINE: FALSCHER TYP'
       STOP
    ENDIF

    IF((N.LT.3).OR.(N.GT.NMAXSPLINE)) THEN
       WRITE(*,*) 'FEHLER IN ROUTINE  SPLINE: N < 3 ODER N > NMAXSPLINE'
       STOP
    ENDIF


!!$C     BERECHNE DIFFERENZ AUFEINENDERFOLGENDER X-WERTE UND
!!$C     UNTERSUCHE MONOTONIE

    DO  I = 1, N-1
       H(I) = X(I+1)- X(I)
       IF(H(I).LE.0.0) THEN
          WRITE(*,*) 'MONOTONIEFEHLER IN SPLINE: X(I-1) >= X(I)',I, N, X(I), X(I+1)
          STOP
       ENDIF
    END DO

!!$C     AUFSTELLEN DES GLEICHUNGSSYSTEMS

    DO I = 1, N-2
       A(I) = 3.0 * ((Y(I+2)-Y(I+1)) / H(I+1) - (Y(I+1)-Y(I)) / H(I))
       B(I) = H(I)
       C(I) = H(I+1)
       D(I) = 2.0 * (H(I) + H(I+1))
    END DO

!!$C     BERUECKSICHTIGEN DER RANDBEDINGUNGEN
!!$C
!!$C     NOT-A-KNOT

    IF(TYP.EQ.0) THEN
       A(1)   = A(1) * H(2) / (H(1) + H(2))
       A(N-2) = A(N-2) * H(N-2) / (H(N-1) + H(N-2))
       D(1)   = D(1) - H(1)
       D(N-2) = D(N-2) - H(N-1)
       C(1)   = C(1) - H(1)
       B(N-2) = B(N-2) - H(N-1)
    ENDIF
!!$C
!!$C     1. ABLEITUNG VORGEGEBEN
!!$C
    IF(TYP.EQ.1) THEN
       A(1)   = A(1) - 1.5 * ((Y(2)-Y(1)) / H(1) - ALFA)
       A(N-2) = A(N-2) - 1.5 * (BETA - (Y(N)-Y(N-1)) / H(N-1))
       D(1)   = D(1) - 0.5 * H(1)
       D(N-2) = D(N-2) - 0.5 * H(N-1)
    ENDIF
!!$C
!!$C     2. ABLEITUNG VORGEGEBEN
!!$C
    IF(TYP.EQ.2) THEN
       A(1)   = A(1) - 0.5 * ALFA * H(1)
       A(N-2) = A(N-2) - 0.5 * BETA * H(N-1)
    ENDIF
!!$C
!!$C     3. ABLEITUNG VORGEGEBEN
!!$C
    IF(TYP.EQ.3 ) THEN
       A(1)   = A(1) + 0.5 * ALFA * H(1) * H(1)
       A(N-2) = A(N-2) - 0.5 * BETA * H(N-1)* H(N-1)
       D(1)   = D(1) + H(1)
       D(N-2) = D(N-2) + H(N-1)
    ENDIF
!!$C
!!$C     BERECHNUNG DER KOEFFIZIENTEN
!!$C

    CALL SGTSL(N-2,B,D,C,A,IERR)

    IF(IERR.NE.0) THEN
       WRITE(*,*) 'ERROR IN SGTSL: MATRIX SINGULAR'
       STOP
    ENDIF
!!$C
!!$C     UEBERSCHREIBEN DES LOESUNGSVEKTORS
!!$C
    CALL SDCOPY(N-2,A,1,C(2),1)
!!$C
!!$C     IN ABHAENGIGKEIT VON DEN RANDBEDINGUNGEN WIRD DER 1. UND
!!$C     DER LETZTE WERT VON C KORRIGIERT
!!$C
    IF(TYP.EQ.0) THEN
       C(1) = C(2) + H(1) * (C(2)-C(3)) / H(2)
       C(N) = C(N-1) + H(N-1) * (C(N-1)-C(N-2)) / H(N-2)
    ENDIF

    IF(TYP.EQ.1) THEN
       C(1) = 1.5*((Y(2)-Y(1)) / H(1) - ALFA) / H(1) - 0.5 * C(2)
       C(N) = -1.5*((Y(N)-Y(N-1)) / H(N-1)-BETA) / H(N-1)-0.5*C(N-1)
    ENDIF

    IF(TYP.EQ.2) THEN
       C(1) = 0.5 * ALFA
       C(N) = 0.5 * BETA
    ENDIF

    IF(TYP.EQ.3) THEN
       C(1) = C(2) - 0.5 * ALFA * H(1)
       C(N) = C(N-1) + 0.5 * BETA * H(N-1)
    ENDIF

    CALL SDCOPY(N,Y,1,A,1)

    DO I = 1, N-1
       B(I) = (A(I+1)-A(I)) / H(I) - H(I) * (C(I+1)+2.0 * C(I)) / 3.0
       D(I) = (C(I+1)-C(I)) / (3.0 * H(I))
    END DO

    B(N) = (3.0 * D(N-1) * H(N-1) + 2.0 * C(N-1)) * H(N-1) + B(N-1)

    RETURN

  END SUBROUTINE BUILD_SPLINE

  REAL FUNCTION SPWERT(N,XWERT,A,B,C,D,X,ABLTG,IEQ)
    ! FROM MISHKA CODE
!!$C-----------------------------------------------------------------------
!!$C     INPUT:
!!$C
!!$C     N           ANZAHL DER KNOTENPUNKTE
!!$C     XWERT       STELLE AN DER FUNKTIONSWERTE BERECHNET WERDEN
!!$C     A, B, C, D  ARRAYS DER SPLINEKOEFFIZIENTEN (AUS SPLINE)
!!$C     X           ARRAY DER KNOTENPUNKTE
!!$C
!!$C     OUTPUT:
!!$C
!!$C     SPWERT   FUNKTIONSWERT AN DER STELLE XWERT
!!$C     ABLTG(I) WERT DER I-TEN ABLEITUNG BEI XWERT
!!$C-----------------------------------------------------------------------
!!$C
    IMPLICIT NONE
    INTEGER,INTENT(IN) ::  N
    REAL, INTENT(IN) :: XWERT
    REAL, DIMENSION(N), INTENT(IN) :: A, B, C, D, X 
    REAL, DIMENSION(3), INTENT(OUT) :: ABLTG
    LOGICAL, INTENT(IN) :: IEQ
    INTEGER ::   I, K, M
    REAL :: XX, DX
!!$C
!!$C     SUCHE PASSENDES INTERVALL (BINAERE SUCHE)
!!$C
    
    IF (IEQ) THEN
       ! equal-distant
       DX = X(2) - X(1)
       I = FLOOR((XWERT - X(1)) / DX) + 1
    ELSE
       I = 1
       K = N

10     M = (I+K) / 2

       IF(M.NE.I) THEN
          IF(XWERT.LT.X(M)) THEN
             K = M
          ELSE
             I = M
          ENDIF
          GOTO 10
       ENDIF
    END IF

    XX = XWERT - X(I)


    ABLTG(1) = (3.0 * D(I) * XX + 2.0 * C(I)) * XX + B(I)
    ABLTG(2) = 6.0 * D(I) * XX + 2.0 * C(I)
    ABLTG(3) = 6.0 * D(I)

    SPWERT = ((D(I)*XX + C(I))*XX + B(I))*XX + A(I)

    RETURN
  END FUNCTION SPWERT

  REAL FUNCTION SPWERT1(N,XWERT,A,B,C,D,X,IEQ)
    ! FROM MISHKA CODE
!!$C-----------------------------------------------------------------------
!!$C     INPUT:
!!$C
!!$C     N           ANZAHL DER KNOTENPUNKTE
!!$C     XWERT       STELLE AN DER FUNKTIONSWERTE BERECHNET WERDEN
!!$C     A, B, C, D  ARRAYS DER SPLINEKOEFFIZIENTEN (AUS SPLINE)
!!$C     X           ARRAY DER KNOTENPUNKTE
!!$C
!!$C     OUTPUT:
!!$C
!!$C     SPWERT   FUNKTIONSWERT AN DER STELLE XWERT
!!$C-----------------------------------------------------------------------
!!$C
    IMPLICIT NONE
    INTEGER,INTENT(IN) ::  N
    REAL, INTENT(IN) :: XWERT
    REAL, DIMENSION(N), INTENT(IN) :: A, B, C, D, X 
    LOGICAL, INTENT(IN) :: IEQ
    INTEGER ::   I, K, M
    REAL :: XX, DX
!!$C
!!$C     SUCHE PASSENDES INTERVALL (BINAERE SUCHE)
!!$C
    
    IF (IEQ) THEN
       ! equal-distant
       DX = X(2) - X(1)
       I = FLOOR(XWERT / DX) + 1
    ELSE
       I = 1
       K = N

10     M = (I+K) / 2

       IF(M.NE.I) THEN
          IF(XWERT.LT.X(M)) THEN
             K = M
          ELSE
             I = M
          ENDIF
          GOTO 10
       ENDIF
    END IF

    XX = XWERT - X(I)

    SPWERT1 = ((D(I)*XX + C(I))*XX + B(I))*XX + A(I)

    RETURN
  END FUNCTION SPWERT1

  REAL FUNCTION SPWERTD1(N,XWERT,A,B,C,D,X,IEQ)
    ! FROM MISHKA CODE
!!$C-----------------------------------------------------------------------
!!$C     INPUT:
!!$C
!!$C     N           ANZAHL DER KNOTENPUNKTE
!!$C     XWERT       STELLE AN DER FUNKTIONSWERTE BERECHNET WERDEN
!!$C     A, B, C, D  ARRAYS DER SPLINEKOEFFIZIENTEN (AUS SPLINE)
!!$C     X           ARRAY DER KNOTENPUNKTE
!!$C
!!$C     OUTPUT:
!!$C
!!$C     SPWERT   FUNKTIONSWERT AN DER STELLE XWERT
!!$C-----------------------------------------------------------------------
!!$C
    IMPLICIT NONE
    INTEGER,INTENT(IN) ::  N
    REAL, INTENT(IN) :: XWERT
    REAL, DIMENSION(N), INTENT(IN) :: A, B, C, D, X 
    LOGICAL, INTENT(IN) :: IEQ
    INTEGER ::   I, K, M
    REAL :: XX, DX
!!$C
!!$C     SUCHE PASSENDES INTERVALL (BINAERE SUCHE)
!!$C
    
    IF (IEQ) THEN
       ! equal-distant
       DX = X(2) - X(1)
       I = FLOOR(XWERT / DX) + 1
    ELSE
       I = 1
       K = N

10     M = (I+K) / 2

       IF(M.NE.I) THEN
          IF(XWERT.LT.X(M)) THEN
             K = M
          ELSE
             I = M
          ENDIF
          GOTO 10
       ENDIF
    END IF

    XX = XWERT - X(I)

    SPWERTD1 = (3.0 * D(I) * XX + 2.0 * C(I)) * XX + B(I)

    RETURN
    
  END FUNCTION SPWERTD1
  
end module spline_module
