MODULE DataParameters
! A module can be used to declare variables used both in the 
! main program and in its subroutines. 
IMPLICIT none
INTEGER, PARAMETER ::            nobs = 376  
INTEGER, PARAMETER ::            n_neighbors = 13 
INTEGER ::                       id(nobs)
INTEGER ::                       weight_index(n_neighbors,nobs)

REAL, DIMENSION(nobs) ::         Y, WY
REAL, DIMENSION(nobs) ::         YTilde, etilde
REAL, DIMENSION(nobs,nobs) ::    W, ITilde

INTEGER, PARAMETER ::            KVars = 2    
REAL, DIMENSION(nobs,KVars) ::   X, WX, XTilde
REAL, DIMENSION(KVars) ::        XTY
REAL, DIMENSION(KVars,KVars) ::  XTX, identity_beta
REAL, DIMENSION(KVars) ::        beta
REAL, DIMENSION(KVars,KVars) ::  var_beta
REAL ::                          sigma_sq, rho
REAL ::                          logdet
REAL ::                          constants, likelihood
REAL ::                          rho_max, rho_limit
REAL ::                          likelihood_OLS, maxlikelihood
REAL ::                          step1, step2, step

! Used in LU factorization 
INTEGER ::                       returncode, ninterchange
INTEGER ::                       index_beta(Kvars)
INTEGER ::                       index_weight(nobs)
END MODULE DataParameters



MODULE LU
! A module can contain functions and subroutines which are called
!  in the program.

CONTAINS

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
 Subroutine LUDCMP(A,N,INDX,D,CODE)
 IMPLICIT none
 INTEGER, INTENT(IN)::   N
 REAL, INTENT(INOUT)::   A(N,N)
 REAL ::                 VV(N)
 INTEGER, INTENT(OUT) :: CODE, D, INDX(N)
 REAL, PARAMETER ::      TINY=1.5E-16
 REAL ::                 AMAX, DUM, SUM0
 INTEGER ::              I, J, K, IMAX

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
 Subroutine LUBKSB(A,N,INDX,B)
 IMPLICIT none
 INTEGER, INTENT(IN)::        N
 INTEGER, INTENT(IN)::        INDX(N)
 REAL, INTENT(IN)::           A(N,N)
 REAL, INTENT(INOUT)::        B(N)

 REAL::                       SUM0
 INTEGER::                    I, II, J, LL

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
END MODULE LU



PROGRAM SEMLU
! This is a main program. Use the above modules. 
USE DataParameters
USE LU
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER, ALLOCATABLE, DIMENSION(:) ::      iwt
CHARACTER(50) ::                           yx_input
CHARACTER(50) ::                           wt_input
!**********************************************************
! This version assumes a SYMMETRIC BINARY weight matrix   *
! It is easily edited for a row-standardized version of   *
! such a weight matrix.                                   *
!                                                         *
! The key to the coding below is the iwt vector. This is  *
! a vector of indices, such that for each observation i,  *
! iwt holds the indices of all neighbors j whose spatial  *
! weights are non-zero. There are at most n_neighbors of  *
! of these, as caculated in the weight matrix program.    *
! The number of neighbors varies depending on i, and in   *
! the weight matrix data for i (the i-th column of the    *
! weight_index matrix) positive indices are associated    *
! with actual neighbors and the remainder of the column   *
! is padded with zeroes. This is why we use the COUNT()   *
! and PACK() commands below, to isolate actual neighbors. *
!                                                         *
! We do not require i to have any neighbors at all.       *
!**********************************************************

!*****************************************************************
! Read data using Fortan open and read function.                 *
!*****************************************************************
 yx_input = "./YXData.raw"
 wt_input = "./weight_index.raw"

 OPEN (Unit = 8, File = yx_input, ACTION="READ")
 DO i = 1, nobs
  READ (Unit = 8, Fmt = *, IOSTAT=ios) Y(i), X(i,2:KVars)
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for Y, X data is ",ios

 OPEN (Unit = 8, File = wt_input, ACTION="READ")
 DO i = 1, nobs
  READ (8, *, IOSTAT=ios) id(i), weight_index(:,i)
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for weight data is ",ios

!*****************************************************************
! Execute subroutine setup                                       *
!*****************************************************************
 call setup

!*****************************************************************
! Calculate weight matrix W as well as WY and WX.                *
!*****************************************************************
! The semicolon (;) is the delimiter. 
 W = 0.0 ; WY = 0.0 ; WX = 0.0

 do i = 1, nobs
  if ( COUNT(weight_index(:,i) > 0) > 0 ) then

    ALLOCATE ( iwt( COUNT(weight_index(:,i) > 0)) )
    iwt = PACK(weight_index(:,i), weight_index(:,i) > 0)

    W(iwt,i) = 1.0         
    WY(i) = sum( Y(iwt) ) 
    do j = 1,KVars
     WX(i,j) = sum( X(iwt,j) )
    end do
  
    DEALLOCATE (iwt)
 
  end if
 end do

!*****************************************************************
! Get OLS (same as ML with rho = 0) results                      *
!*****************************************************************
 call ols

!*****************************************************************
! Begin SEM grid search over rho using the concentrated likelihood
!*****************************************************************
 step1 = 0.025   ! crude grid step
 step2 = 0.0025   ! more refined grid step
 maxlikelihood = likelihood_OLS
 rho_max = 0.0

 rho = step1     ! assuming rho is positive, start grid here

 step = step1
 rho_limit = 1.0

 CALL SEARCH 
 print *, "The best initial rho value is ", rho_max, &
                                  " giving ", maxlikelihood
 print *, "versus the OLS likelihood ", likelihood_OLS

 rho = rho_max - step1 + epsilon(1.0)
 step = step2
 rho_limit = rho_max + step1

 CALL SEARCH
 print *, "Final estimate of rho is ", rho_max, &
                                   "giving ", maxlikelihood

 rho = rho_max

!*****************************************************************
! Get GLS results                                                *
!*****************************************************************
 call gls

STOP
END PROGRAM SEMLU



SUBROUTINE SETUP
USE DataParameters
IMPLICIT none
INTEGER ::          i

! Add a ones column to the X matrix and create an identity matrix
!  to be used in finding the variance of beta
 X(:,1) = 1.0

 identity_beta = 0.0
 do i = 1, KVars
  identity_beta(i,i) = 1.0
 enddo

! Create constants to aid in comparing log-likelihood values
constants = -(real(nobs)/2.0)* log(2.0*3.1415927) - real(nobs)/2.0

RETURN
END SUBROUTINE SETUP



SUBROUTINE OLS
USE DataParameters
USE LU
IMPLICIT none
INTEGER::          i, j 
!*****************************************************************
! Solve for OLS estimator of beta, using the LU decomposition of * 
! X'X and right-hand side vector X'Y.                            *
!*****************************************************************
 do i = 1, KVars
   do j = 1, KVars
     XTX(i,j) = dot_product(X(:,i),X(:,j))
   enddo
   XTY(i) = dot_product(X(:,i),Y)
 enddo

 call ludcmp(XTX,KVars,index_beta,ninterchange,returncode)
 
 if (returncode == 0) then
  call lubksb(XTX,KVars,index_beta,XTY)
  beta = XTY
  print *, "OLS beta is ", beta
 else
  print *, "problem with OLS beta ",returncode
 endif

 sigma_sq = dot_product(Y - matmul(X,beta), Y - matmul(X,beta))/real(nobs)
 print *, "OLS sigma squared is ", sigma_sq

 likelihood_OLS = constants - (real(nobs)/2.0)*log(sigma_sq)

RETURN
END SUBROUTINE OLS



SUBROUTINE SEARCH
USE DataParameters
USE LU
IMPLICIT none
INTEGER ::           i, j

OUTER: DO WHILE (rho < rho_limit)

  YTilde = Y - rho*WY
  XTilde = X - rho*WX
  ITilde = 0.0
  do i=1, nobs; ITilde(i,i) = 1.0; enddo
  WHERE (W /= 0)                   ! Exploit sparse, binary W
   ITilde = ITilde - rho
  END WHERE

!*****************************************************************
! Get the log of the Jacobian for this value of rho. The log of  *
! the Jacobian is the sum of logs of absolute values of diagonal *
! elements of U, the upper triangular matrix of (I-rhoW)         *          
!*****************************************************************  
  call ludcmp(ITilde,nobs,index_weight,ninterchange,returncode)
  if (returncode == 0) then
   logdet = 0.0
   do i = 1, nobs
     logdet = logdet + log(abs(ITilde(i,i)))
   enddo
  else
   print *, "problem with Jacobian for rho equal to ", rho
  endif

!*****************************************************************
! Solve for beta, given rho,  using the LU decomposition of      *
! XTilde'XTilde and XTilde'YTilde. Use the columns of XTilde     *
! and dot products for greater efficiency                        *
!*****************************************************************
  do i = 1, KVars
    do j = 1, KVars
      XTX(i,j) = dot_product(XTilde(:,i),XTilde(:,j))
    enddo
    XTY(i) = dot_product(XTilde(:,i),YTilde)
  enddo

  call ludcmp(XTX,KVars,index_beta,ninterchange,returncode)
 
  if (returncode == 0) then
   call lubksb(XTX,KVars,index_beta,XTY)
   beta = XTY
  else
   print *, "problem with betatilde ",returncode
  endif

  etilde = YTilde - matmul(XTilde,beta)
  sigma_sq = dot_product(etilde,etilde)/real(nobs)

  likelihood = constants -(real(nobs)/2.0)*log(sigma_sq) + logdet

  print *, "rho", rho, "logdet", logdet
 
  if (likelihood > maxlikelihood) then
   maxlikelihood = likelihood
   rho_max = rho
  endif

  rho = rho + step
ENDDO OUTER

RETURN
END SUBROUTINE SEARCH


SUBROUTINE GLS
USE DataParameters
USE LU
IMPLICIT none
INTEGER::          i, j 
!*****************************************************************
! Solve for GLS estimator of beta, using final rho value         *
!*****************************************************************
 YTilde = Y - rho*WY
 XTilde = X - rho*WX

 do i = 1, KVars
  do j = 1, KVars
    XTX(i,j) = dot_product(XTilde(:,i),XTilde(:,j))
  enddo
  XTY(i) = dot_product(XTilde(:,i),YTilde)
 enddo

 call ludcmp(XTX,KVars,index_beta,ninterchange,returncode)
 
 if (returncode == 0) then
  call lubksb(XTX,KVars,index_beta,XTY)
  beta = XTY
  print *, "Final beta is ", beta
 else
  print *, "Problem with final beta ",returncode
 endif

 etilde = YTilde - matmul(XTilde,beta)
 sigma_sq = dot_product(etilde,etilde)/real(nobs)
 print *, "Final sigma squared is ", sigma_sq

! Invert XTX 
 do i = 1, KVars
  call lubksb(XTX,KVars,index_beta,identity_beta(:,i))
 enddo

! Multiply the inverse by sigma-squared to obtain 
!  the variance matrix for beta
 var_beta = sigma_sq*identity_beta

! print the estimation results. 
 print *, "beta, stand. dev., t-statistic"
 do i = 1, KVars
   print *, beta(i), sqrt(var_beta(i,i)), &
                          beta(i)/sqrt(var_beta(i,i))
 enddo
 
 
RETURN
END SUBROUTINE GLS



