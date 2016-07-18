! Contain interfaces of BLAS and LAPACK subroutines
!
! CURRENT SETTINGS : compile using real 8, call zxxxx, dxxxx
! MAY need to change when using a different compiler or 
! to compile on a different machine
!
! TO SWITCH : simply between zxxxx and cxxxx, dxxxx and sxxxx
! in each subroutine

module interfaces

  implicit none

  public

  interface

     ! interface for external functions
     subroutine sgtsl(n,c,d,e,b,info)
       ! sgtsl from LINPACK
       integer, intent(in) :: n
       integer, intent(out) :: info
       real c(n), d(n), e(n), b(n)
     end subroutine sgtsl
  end interface
  
contains

  SUBROUTINE SDCOPY(N,SX,INCX,SY,INCY)
    ! DCOPY from BLAS
    INTEGER, INTENT(IN) :: INCX,INCY,N
    REAL SX(*), SY(*)

    CALL DCOPY(N, SX, INCX, SY, INCY)
  END SUBROUTINE SDCOPY

  SUBROUTINE CZSCAL(N,ZA,ZX,INCX)
    ! C(Z)SCAL from BLAS
    ! X = ZA * X
    complex, intent(in) ::  ZA
    integer, intent(in) :: INCX,N
    complex ZX(*)

    CALL ZSCAL(N, ZA, ZX, INCX)

  END SUBROUTINE CZSCAL

  SUBROUTINE CZCOPY(N,SX,INCX,SY,INCY)
    ! C(Z)COPY from BLAS
    ! Y = X
    integer, intent(in) :: INCX,INCY,N
    complex SX(*), SY(*)

    CALL ZCOPY(N, SX, INCX, SY, INCY)
  END SUBROUTINE CZCOPY

  SUBROUTINE CZAXPY(N,ZA,ZX,INCX,ZY,INCY)
    ! C(Z)AXPY from BLAS
    ! Y = ZA * X + Y
    complex, intent(in) ::  ZA
    integer, intent(in) ::  INCX,INCY,N
    COMPLEX ZX(*),ZY(*)

    CALL ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
  END SUBROUTINE CZAXPY

  SUBROUTINE CZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    ! C(Z)GEMV FROM BLAS
    ! Y = ALPHA * A * X + BETA * Y
      COMPLEX, INTENT(IN) ::  ALPHA,BETA
      INTEGER, INTENT(IN) ::  INCX,INCY,LDA,M,N
      CHARACTER, INTENT(IN) :: TRANS

      COMPLEX A(LDA,*),X(*),Y(*)

      CALL ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  END SUBROUTINE CZGEMV

  SUBROUTINE CZGETRF( M, N, A, LDA, IPIV, INFO )
    ! C(Z)GETRF from LAPACK
    ! PLU FACTORIZATION

    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO
    INTEGER            IPIV( * )
    COMPLEX            A( LDA, * )

    CALL ZGETRF( M, N, A, LDA, IPIV, INFO )
  END SUBROUTINE CZGETRF

  SUBROUTINE CZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    ! C(Z)GETRS from LAPACK
    ! SOLVING THE LINEAR SYSTEM A * X = B USING LU FACTORIZATION

     CHARACTER, INTENT(IN) :: TRANS
     INTEGER, INTENT(IN) ::  LDA, LDB, N, NRHS
     INTEGER, INTENT(OUT) :: INFO
     INTEGER            IPIV( * )
     COMPLEX            A( LDA, * ), B( LDB, * )

     CALL ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  END SUBROUTINE CZGETRS
  
  COMPLEX FUNCTION CZDOTU(N,CX,INCX,CY,INCY)
    ! C(Z)DOTU FROM BLAS
    ! COPYED DIRECTLY FROM SOURCE FILE
    ! RETURN X**T * Y
    INTEGER, INTENT(IN) :: INCX,INCY,N
    COMPLEX CX(*),CY(*)

      COMPLEX CTEMP
      INTEGER I,IX,IY
      CTEMP = (0.0,0.0)
      CZDOTU = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CTEMP = CTEMP + CX(I)*CY(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CX(IX)*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CZDOTU = CTEMP

  END FUNCTION CZDOTU

  COMPLEX FUNCTION CZDOTC(N,CX,INCX,CY,INCY)
    ! C(Z)DOTU FROM BLAS
    ! COPY DIRECTLY FROM THE SOURCE FILE
    ! RETURN X**H * Y
    INTEGER, INTENT(IN) :: INCX,INCY,N
    COMPLEX CX(*),CY(*)
    COMPLEX ZDOTC

      COMPLEX CTEMP
      INTEGER I,IX,IY

      INTRINSIC CONJG
      CTEMP = (0.0,0.0)
      CZDOTC = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(I))*CY(I)
         END DO
      ELSE

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CZDOTC = CTEMP
      RETURN
  END FUNCTION CZDOTC
  
end module interfaces
