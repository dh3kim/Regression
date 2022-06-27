MODULE Regress

CONTAINS

SUBROUTINE OLS_LU(Y, X, nobs, KVars, beta, sigma_sq, var_beta, e)
IMPLICIT NONE
REAL, INTENT(IN), DIMENSION(nobs):: Y
REAL, INTENT(IN), DIMENSION(nobs, KVars):: X
INTEGER, INTENT(IN):: nobs, KVars
REAL, INTENT(OUT), DIMENSION(KVars):: beta
REAL, INTENT(OUT):: sigma_sq
REAL, INTENT(OUT), DIMENSION(KVars, KVars):: var_beta
REAL, INTENT(OUT), DIMENSION(nobs), optional:: e 

REAL:: XTX(KVars, KVars), XTY(KVars)
REAL:: identity_beta(KVars, KVars)
INTEGER:: index_beta(KVars), ninterchange, returncode

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
! Solve (X'X)*beta = X'Y for beta using LU decomposition.        *
! Use subroutines LUDCMP and LUBSKSB. The CALL statement is      *
! used to call a subroutine.                                     *
!*****************************************************************
  CALL LUDCMP(XTX,KVars,index_beta,ninterchange,returncode)
 
  if (returncode == 0) then
    CALL LUBKSB(XTX,KVars,index_beta,XTY)
    beta = XTY
  else
    print *, "problem with OLS beta ", returncode
  endif

!*****************************************************************
! Calculate OLS residuals, OLS sigma squared, and variance       *
! matrix for beta.                                               *
!*****************************************************************
 e = Y - MATMUL(X, beta)
 sigma_sq = DOT_PRODUCT(e,e)/REAL(nobs-KVars)

! Invert XTX. Note that XTX is a LU decomposition of XTX here.  
   identity_beta = 0.0
   do i = 1, KVars
     identity_beta(i,i) = 1.0
    enddo

  do i = 1, KVars
   call LUBKSB(XTX,KVars,index_beta,identity_beta(:,i))
  end do

! Multiply the inverse by sigma-squared to obtain 
!  the variance matrix for beta
   var_beta = sigma_sq*identity_beta

END SUBROUTINE OLS_LU

END MODULE Regress

