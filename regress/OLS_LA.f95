MODULE OLS

CONTAINS

SUBROUTINE OLS_LA(Y, X, nobs, KVars, beta, sigma_sq, var_beta, e)
USE LA_PRECISION, ONLY:    WP => SP
USE F95_LAPACK, ONLY:      LA_POSV, LA_SYSV
IMPLICIT NONE
REAL, INTENT(IN), DIMENSION(nobs):: Y
REAL, INTENT(IN), DIMENSION(nobs, KVars):: X
INTEGER, INTENT(IN):: nobs, KVars
REAL, INTENT(OUT), DIMENSION(KVars):: beta
REAL, INTENT(OUT):: sigma_sq
REAL, INTENT(OUT), DIMENSION(KVars, KVars):: var_beta
REAL, INTENT(OUT), DIMENSION(nobs), optional:: e

REAL, DIMENSION(KVars, KVars):: XTX, XX, XTXinverse
REAL:: XTY(KVars)
INTEGER:: i, j
!*****************************************************************
! Calculate X'X and X'Y by using DOT_PRODUCT which is one of     *
! Fortran intrinsic function.                                    *
!*****************************************************************
  do i = 1, KVars
    do j = 1, KVars
      XTX(i,j) = DOT_PRODUCT(X(:,i),X(:,j))
    enddo
    XTY(i) = DOT_PRODUCT(X(:,i),Y)
  enddo

!*****************************************************************
! Solve (X'X)*beta = X'Y for beta using LA_POSV.                 *
!*****************************************************************
  XX = XTX
  CALL LA_POSV(XTX,XTY)
  beta = XTY

!*****************************************************************
! Calculate OLS residuals, OLS sigma squared, and variance       *
! matrix for beta.                                               *
!*****************************************************************
 e = Y - MATMUL(X, beta)
 sigma_sq = DOT_PRODUCT(e,e)/REAL(nobs-KVars)

! Invert XTX.
  XTXinverse = 0.0
  do i = 1, KVars
    XTXinverse(i,i) = 1.0
  end do
  CALL LA_POSV(XX, XTXinverse)

! Multiply the inverse by sigma-squared to obtain
!  the variance matrix for beta
   var_beta = sigma_sq*XTXinverse

END SUBROUTINE OLS_LA

END MODULE OLS
