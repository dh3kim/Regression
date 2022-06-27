!********************************************************************
! Bayeisan Estimation of Linear Regression Model                    *
! Using independent Normal-Gamma Prior                              *
!******************************************************************** 
! Donghwan Kim                                                      *
!   Department of Economics                                         *
!   Stony Brook University                                          *
!********************************************************************

MODULE DataParameters
IMPLICIT none
 INTEGER, PARAMETER          :: nobs = 546 ! number of observations
 INTEGER, PARAMETER          :: KVars = 5  ! number of explanatory vars 
 REAL, DIMENSION(nobs)       :: Y
 REAL, DIMENSION(nobs,KVars) :: X

! Used in LU factorization 
 INTEGER :: returncode, ninterchange
 INTEGER :: index_beta(KVars)

! Used in Cholesky Decompostion
 INTEGER :: check
 
! Used to check CPU time
 REAL :: start, finish

! Used in Random Sample from gamma dist.
 LOGICAL :: First

!Bayesian
 INTEGER, PARAMETER :: ndraw = 100000
 INTEGER, PARAMETER :: nburn = 50000


! Used to save samples
 REAL, DIMENSION(ndraw-nburn, KVars) :: betasave
 REAL, DIMENSION(ndraw-nburn)        :: sigma_save

! Sampling beta
 REAL, DIMENSION(KVars)        ::  XTY
 REAL, DIMENSION(KVars,KVars)  ::  XTX
 REAL, DIMENSION(KVars)        ::  beta, beta_part
 REAL, DIMENSION(KVars)        ::  beta0, beta1
 REAL, DIMENSION(KVars, KVars) ::  M0, M1
 REAL, DIMENSION(KVars, KVars) ::  identity_M1, inverseM1
 REAL, DIMENSION(KVars)        :: norm_rnd

! Sampling sigma_sq
 REAL, DIMENSION(nobs) :: epsilon
 REAL                  :: sigma_sq
 REAL                  :: S0, S1, v0, v1

! Posterior mean and std
 REAL, DIMENSION(KVars) :: sum_beta, postmean_beta
 REAL, DIMENSION(KVars) :: sum_sq_beta, postvar_beta
 REAL                   :: sum_sigma, postmean_sigma
 REAL                   :: sum_sq_sigma, postvar_sigma

END MODULE DataParameters



MODULE LU
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



MODULE Cholesky
CONTAINS
 Subroutine chodec(n, a, rc)
!**********************************************************************
! It is modified from a code by Jean-Pierre Moreau                    *
!======================================================================
!*                                                                    *
!*  chodec decomposes the symmetric positive definite matrix mat.     *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      n        integer (n > 0)                                      *
!*               Dimension of matrix a                                *
!*      a        REAL matrix (n,n)                                    *
!*               Matrix of left hand coefficients                     *
!*                                                                    *
!*   Output parameter:                                                *
!*   ================                                                 *
!*      a        REAL matix (n,n)                                     *
!*               Cholesky decomposition in the lower triangle         *
!*                                                                    *
!*   Return value rc:                                                 *
!*   ===============                                                  *
!*      = 0      all ok                                               *
!*      = 1      n < 1                                                *
!*      = 2      Matrix not  positive definite                        *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Functions in use:   dsqrt (square root in double precision)      *
!*   ================                                                 *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Constants in use:   EPSQUAD (small number)                       *
!*   ================                                                 *
!*                                                                    *
!======================================================================

IMPLICIT none
INTEGER,INTENT(IN) :: n
INTEGER, INTENT(OUT) :: rc
REAL, DIMENSION(n,n),INTENT(INOUT) ::a(0:n-1, 0:n-1)
REAL  :: EPSQUAD, sum
INTEGER :: j, i, k

  EPSQUAD = 1.d-12

  if (n < 1) then
    rc=1                           ! n < 1  error
    return
  end if
  if (a(0,0) < EPSQUAD) then       ! matrix a not positive definite
    rc=2
    return
  end if
  a(0,0) = sqrt(a(0,0))
  do j = 1, n-1
    a(j,0) = a(j,0) / a(0,0)
  end do

  do i = 1, n-1
    sum = a(i,i)
    do j = 0, i-1
      sum = sum - a(i,j)*a(i,j)
    end do

    if (sum < EPSQUAD) then         ! matrix a not positive definite
      rc=2
      return
    end if
    a(i,i) = sqrt(sum)
    do j = i+1, n-1
      sum = a(j,i)
      do k = 0, i-1
        sum = sum - a(i,k) * a(j,k)
      end do
      a(j,i) = sum / a(i,i)
    end do
  end do

do i = 0, n-2
  do j = i+1, n-1
   a(i, j) = 0.0
  end do
 end do


  rc = 0  ! all Ok
  return
END SUBROUTINE chodec
END MODULE Cholesky



MODULE Random
!***************************************************************
! This module comes from Alan Miller Website                   *
!*************************************************************** 
REAL, PRIVATE :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, &
                  vsmall = TINY(1.0), vlarge = HUGE(1.0)

CONTAINS

FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal


FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s -1.0/3.0
  c = 1.0/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (1.0 + c*x)**3.0
    IF (v > 0.0) EXIT

  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1.0 - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(1.0 - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1

END MODULE Random



PROGRAM BayRegress
!USE DataPaths
USE DataParameters
USE LU
USE Cholesky
USE Random 
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER, ALLOCATABLE, DIMENSION(:) ::      iwt
INTEGER ::                                 iter
!***************************************************************
! Bayeisan Estimation of Linear Regression Model
! Using independent Normal-Gamma Prior
!***************************************************************
! Model: Y = Xbeta + epsilon. epsilon ~ N(0, sigma_sq*I)
! Prior: beta ~ N(beta0, M0^(-1))
!        sigma_sq ~ iG(v0/2, S0/2)
!***************************************************************
 
 OPEN (Unit=8, File="./HPRICE.TXT", ACTION="READ")
 OPEN (Unit=10, File="./Result", ACTION="WRITE")
 OPEN (Unit=12, File="./Samples", ACTION="WRITE")

 ! read Data
 ! Data HPRICE.TXT comes from Gary Koop (2003)'s book           

  DO i = 1, nobs
   READ (Unit = 8, Fmt = *, IOSTAT=ios)  Y(i), X(i, 2:KVars)
  END DO

  CLOSE (Unit = 8)
  print *, "IO status is ",ios

!****************************************************************
! Gibbs Sampling
!****************************************************************
 ! hyperparamters in prior
   CALL SETUP

 ! Storage for Samples
   betasave = 0.0
   sigma_save = 0.0

 ! initial values and quantities used in Gibbs Sampling
   CALL OLS 
   beta = 0.0
   sigma_sq = 100.0
   
 ! quantities used in Gibbs but not changed in iteration
   do i=1, KVars
    do j=1, KVars
      XTX(i,j)=dot_product(X(:,i), X(:,j))
    end do
      XTY(i)=dot_product(X(:,i),Y)
   end do

   v1 = v0 + nobs

 ! Gibbs Sampler
   iter = 1

   print *, "Gibbs Sampler runs"

   CALL CPU_TIME(start)
   CALL RANDOM_SEED
    First=.true.

GIBBS: DO WHILE ( iter <= ndraw)
 
 !*******************************
 ! Sample a vector of beta from
 ! full conditional distribution
 ! beta ~ N(beta1, M1^(-1))
 !*******************************
  M1 = M0 + (1.0/sigma_sq)* XTX
  beta_part = matmul(M0, beta0) +(1.0/sigma_sq)* XTY
  
 ! invert M1 using LU Decomposition
   call ludcmp(M1, KVars,  index_beta, ninterchange, returncode)
   inverseM1 = identity_M1
   do i=1, KVars
    call lubksb(M1, KVars, index_beta, inverseM1(:,i))
   end do

 ! mean of beta
   beta1 = matmul(inverseM1, beta_part)

 ! Sample from N(0,1)
   do i = 1, KVars
     norm_rnd(i) = random_normal()
   end do
 
 ! generate beta
   call chodec(KVars,inverseM1,check)
   if (check == 2) then
    print *, "inverseM1 not a positive definite"
    pause
   end if

   beta = matmul(inverseM1,norm_rnd) +  beta1

 !************************************
 ! Sample sigma_sq
 ! sigma_sq ~ iG(v1/2, S1/2)
 !***********************************
    epsilon = Y - matmul(X,beta)
    S1 = dot_product(epsilon, epsilon) + S0  
   
    ! Sample sigma_sq 
    ! if (iter ==1) then 
    !     First=.true. 
    ! end if
     sigma_sq = S1/(2.0 * random_gamma1(v1/2.0, First))
    ! First = .false.

 ! store samples
   if (iter > nburn) then
     betasave(iter-nburn, 1:KVars) = beta
     sigma_save(iter-nburn) = SQRT(sigma_sq) 
   end if

   iter = iter + 1
END DO GIBBS

   CALL CPU_TIME(finish)
   print *, "end of Gibbs Sampler"
   print *, "time", finish-start, "sec."

  !*************************************************
  ! calculate posterior mean and stadard deviation * 
  !*************************************************
    sum_beta = 0.0
    sum_sigma = 0.0

    do i = 1, ndraw-nburn
     do j = 1, KVars
       sum_beta(j) = sum_beta(j) + betasave(i,j)
     end do
       sum_sigma = sum_sigma + sigma_save(i)
    end do

    do i = 1, KVars
       postmean_beta(i) = sum_beta(i)/real(ndraw-nburn)
    end do
     postmean_sigma = sum_sigma/real(ndraw-nburn)

     sum_sq_beta = 0.0
     sum_sq_sigma = 0.0

    do i = 1, ndraw-nburn
      do j = 1, KVars
        sum_sq_beta(j) = sum_sq_beta(j) + (betasave(i,j) - postmean_beta(j))**2
      end do
        sum_sq_sigma = sum_sq_sigma +&
                              (sigma_save(i) - postmean_sigma)**2
    end do

    do i = 1, KVars
       postvar_beta(i) = sum_sq_beta(i)/real(ndraw-nburn-1)
    end do
       postvar_sigma = sum_sq_sigma/real(ndraw-nburn-1)

   ! print
      print *, " posterior mean and standard deviation of beta"
      do i = 1, KVars
      print *, postmean_beta(i), SQRT(postvar_beta(i))
      end do
      print *, "posterior mean and std of sigma"
      print *, postmean_sigma, SQRT(postvar_sigma)

    ! save samples
     do i = 1, ndraw-nburn
     write(12, *) betasave(i,1:KVars), sigma_save(i)
     end do
END PROGRAM BayRegress



SUBROUTINE SETUP
USE DataParameters
IMPLICIT none
INTEGER ::          i
INTEGER ::   Prior

! add constant terms in matrix of independent variables
  X(:,1) = 1.0

! identity matrix
  identity_M1 = 0.0
  do i = 1, KVars
    identity_M1(i,i) = 1.0
  end do

Prior = 1 

if (Prior == 1) then
! hyperparameters in prior distributions  
! beta
  beta0(1) = 0.0
  beta0(2) = 10.0
  beta0(3) = 5000.0
  beta0(4) = 10000.0
  beta0(5) = 10000.0

  M0 = 0.0
  M0(1,1) = 1/10000.0**2
  M0(2,2) = 1/25.0
  M0(3,3) = 1/2500.0**2
  M0(4,4) = 1/5000.0**2
  M0(5,5) = 1/5000.0**2

! Simga_sq
  v0 = 2.5 
  S0 = 25000.0
 
else if (Prior == 2) then
! prior
 do i =1, KVars
 beta0(i) = 0.0
 end do

 M0 = 0.0
 do i=1, KVars
   M0(i,i) = 0.00000000001
 end do

 v0 = 5.0
 S0 = 500.0

end if

RETURN
END SUBROUTINE SETUP



SUBROUTINE OLS
USE DataParameters
USE LU
IMPLICIT none
INTEGER:: i, j

!************************************************************
! Classical OLS estimator                                   *
! Just to compare it with bayesian estimator                *
!************************************************************
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
   print *, "OLS beta is "
    do i = 1, KVars
     print *, beta(i)
    end do
  else
   print *, "problem with OLS beta ", returncode
  endif

  sigma_sq = dot_product(Y - matmul(X,beta), Y - matmul(X,beta))/real(nobs)
  print *, "OLS sigma squared is ", sigma_sq
  print *, "OLS sigma is", SQRT(sigma_sq)

RETURN
END SUBROUTINE OLS

