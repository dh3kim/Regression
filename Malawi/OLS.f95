PROGRAM OLS
! A line with exclamation mark (!) is just a comment line. 
! Start with a statement IMPLICIT NONE and declare all variables
! used in the program. Note that Fortran is case-insensitive.
!##################################################################
! Variable Declaration.                                           # 
!##################################################################
IMPLICIT none
INTEGER, PARAMETER ::            nobs = 376  
INTEGER, PARAMETER ::            KVars = 2   
REAL, DIMENSION(nobs) ::         Y
REAL, DIMENSION(nobs,KVars) ::   X
REAL, DIMENSION(KVars) ::        XTY, beta
REAL, DIMENSION(KVars,KVars) ::  XTX, identity_beta, var_beta
REAL, DIMENSION(nobs) ::         e
REAL ::                          sigma_sq
! Used in LU factorization 
INTEGER ::                       returncode, ninterchange
INTEGER ::                       index_beta(KVars)
INTEGER ::                       i, j, ios

!*****************************************************************
! Read data using OPEN and READ statements                       * 
!*****************************************************************
! Open the YXData.raw file which is in the current
! directory (./) and assign unit number 8 to the file
! for reading. The quotation martk(") or apostrope(')
! is used to indicate a character string. 
 OPEN (Unit = 8, File = "./YXData.raw", ACTION="READ")

! Read the file in unit 8 and close the unit.  
 DO i = 1, nobs
  READ (Unit = 8, Fmt = *, IOSTAT=ios) Y(i), X(i,2:KVars) 
 END DO
 CLOSE (Unit = 8)

! print the ios variable on the default output unit (screen) 
! with default format(*). The ios is an integer variable which 
! is set to zero if the statement is executed successfully.  
 PRINT *, "IO status is ", ios

!*****************************************************************
! Create a OLSResiduals.raw file for saving OLS residuals later. *
!*****************************************************************
 OPEN (Unit = 10, File="./OLSResiduals.raw", ACTION="WRITE")

!*****************************************************************
! Add a ones column to the X matrix and create an identity       * 
! matrix to be used in finding the varaince of beta              *
!*****************************************************************
 X(:,1) = 1.0

 identity_beta = 0.0
 do i = 1, KVars
  identity_beta(i,i) = 1.0
 enddo

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
  do i = 1, KVars
   call LUBKSB(XTX,KVars,index_beta,identity_beta(:,i))
  end do

! Multiply the inverse by sigma-squared to obtain 
!  the variance matrix for beta
   var_beta = sigma_sq*identity_beta

!*****************************************************************
! print the results on the screen and save OLS residuals         * 
! to the file in unit 10. The ampersand (&) is just line         *
! contiuation marker.                                            *
!*****************************************************************
 print *, "beta, standard deviation, t-stat."
 do i = 1, KVars
  print *, beta(i), SQRT(var_beta(i,i)), &
                      beta(i)/SQRT(var_beta(i,i))
 enddo
 print *, "sigma is", SQRT(sigma_sq)

 DO i = 1, nobs
  WRITE(10,*) e(i)
 END DO

END PROGRAM OLS



Subroutine LUDCMP(A,N,INDX,D,CODE)
IMPLICIT none
INTEGER, INTENT(IN)::   N
REAL, INTENT(INOUT)::   A(N,N)
REAL ::                 VV(N)
INTEGER, INTENT(OUT) :: CODE, D, INDX(N)
REAL, PARAMETER ::      TINY=1.5E-16
REAL ::                 AMAX, DUM, SUM0
INTEGER ::              I, J, K, IMAX
!  ***************************************************************
!  * This is a rewrite by Jean-Pierre Moreau of F77 code from    *
!  * Numerical Recipes. I should check it against the original.  *
!  * There is F90 code in the newer edition of the book, but it  *
!  * differs in a couple of ways that need looking at.           *
!  *                                                             *
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************

 D=1; CODE=0

 DO I=1, N
   AMAX=0.0
   DO J=1, N
     IF ( ABS(A(I,J)) > AMAX ) AMAX = ABS(A(I,J))
   END DO ! j loop
   IF(AMAX < TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.0 / AMAX
 END DO ! i loop

 DO J=1, N
   DO I=1, J-1
     SUM0 = A(I,J)
     DO K=1, I-1
       SUM0 = SUM0 - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM0
   END DO ! i loop
   AMAX = 0.0
   DO I=J, N
     SUM0 = A(I,J)
     DO K=1, J-1
       SUM0 = SUM0 - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM0
     DUM = VV(I)*ABS(SUM0)
     IF(DUM >= AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J /= IMAX) THEN
     DO K=1, N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF( ABS(A(J,J)) < TINY ) A(J,J) = TINY

   IF(J /= N) THEN
     DUM = 1.0 / A(J,J)
     DO I=J+1, N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

RETURN
END subroutine LUDCMP


Subroutine LUBKSB(A,N,INDX,B)
IMPLICIT none
INTEGER, INTENT(IN)::        N
INTEGER, INTENT(IN)::        INDX(N)
REAL, INTENT(IN)::           A(N,N)
REAL, INTENT(INOUT)::        B(N)
REAL::                       SUM0
INTEGER::                    I, II, J, LL
!  *********************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is        *
!  * input, not as the matrix A but rather as its LU decomposition,    *
!  * determined by the routine LUDCMP. INDX is input as the permuta-   *
!  * tion vector returned by LUDCMP. B is input as the right-hand      *
!  * side vector B, and returns with the solution vector X. A, N and   *
!  * INDX are not modified by this routine and can be used for suc-    *
!  * cessive calls with different right-hand sides. This routine is    *
!  * also efficient for plain matrix inversion.                        *
!  *********************************************************************

 II = 0

 DO I=1, N
   LL = INDX(I)
   SUM0 = B(LL)
   B(LL) = B(I)
   IF(II /= 0) THEN
     DO J=II, I-1
       SUM0 = SUM0 - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM0 /= 0.0) THEN
     II = I
   END IF
   B(I) = SUM0
 END DO ! i loop

 DO I=N, 1, -1
   SUM0 = B(I)
   IF(I < N) THEN
     DO J=I+1, N
       SUM0 = SUM0 - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM0 / A(I,I)
 END DO ! i loop

RETURN
END subroutine LUBKSB


