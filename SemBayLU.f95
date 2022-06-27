MODULE DataPaths
IMPLICIT none
CHARACTER(100)::                 path = "./"
CHARACTER(200)::                 input, wt_input, output
END MODULE DataPaths

MODULE DataParameters
IMPLICIT none
 INTEGER, PARAMETER ::           nobs = 376
 INTEGER, PARAMETER ::           KVars = 2
 INTEGER, PARAMETER ::           n_neighbors = 13 ! Maximum number of neighbor
 INTEGER, DIMENSION(nobs) ::     id
 REAL, DIMENSION(nobs) ::        Y
 REAL, DIMENSION(nobs,KVars) ::  X
 INTEGER ::                      weight_index(n_neighbors,nobs)

! Used in LU factorization
 INTEGER ::                      returncode, ninterchange
 INTEGER ::                      index_weight(nobs), index_beta(KVars)

! Used in Bayesian Estimation
 INTEGER, PARAMETER ::   ndraw = 10000   ! number of draw
 INTEGER, PARAMETER ::   nburn = 5000    ! number of burn-in
 INTEGER, PARAMETER ::   npar = KVars+2 ! number of parameters 

! Used for saving samples
 REAL, DIMENSION(ndraw-nburn, npar) :: AllDraws 

! Used in sampling beta
 REAL, DIMENSION(nobs,nobs) ::     W
 REAL, DIMENSION(nobs) ::          WY
 REAL, DIMENSION(nobs,KVars) ::    WX

 REAL, DIMENSION(nobs) ::          YTilde
 REAL, DIMENSION(nobs,KVars) ::    XTilde

 REAL, DIMENSION(KVars) ::         XTY
 REAL, DIMENSION(KVars,KVars)  ::  XTX
 REAL, DIMENSION(KVars) ::         beta, beta_part
 REAL, DIMENSION(KVars) ::         beta0, beta1
 REAL, DIMENSION(KVars, KVars) ::  M0, M1
 REAL, DIMENSION(KVars, KVars) ::  identity_M1, inverseM1
 REAL, DIMENSION(KVars) ::         norm_rnd

! Used in sampling sigma_sq
 REAL, DIMENSION(nobs) ::          epsilon
 REAL ::                           sigma_sq
 REAL ::                           S0, S1, v0, v1

! Used in sampling rho
 REAL ::                           rho
 REAL ::                           rho_min, rho_max
 REAL ::                           rhostar, rhotemp
 INTEGER ::                        accept, acc
 REAL ::                           lnp, lnpstar
 REAL ::                           ratio, rnd, cc
 REAL, DIMENSION(ndraw) ::         acc_rate
 REAL ::                           prob
 REAL, DIMENSION(nobs,nobs) ::     ITilde
 REAL ::                           lnconditionalrho

! Used in posterior analysis 
 REAL, DIMENSION(npar) ::          PostMeans, PostVars 
 REAL, DIMENSION(npar) ::          NSE, RNE, CD 
 REAL ::                           frac1, frac2

END MODULE DataParameters


Module Precision
 Integer, Parameter :: sp = Kind(0.0E0)
 Integer, Parameter :: dp = Kind(0.0D0)
 Integer, Parameter :: wp = sp
End Module Precision


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
! It is modified from code by Jean-Pierre Moreau                      *
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
!**************************************************************************
!                                                                         *
! Collection of random number generators                                  *
!                                                                         *
!**************************************************************************

REAL, PRIVATE :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, &
                  vsmall = TINY(1.0), vlarge = HUGE(1.0)

CONTAINS

FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
REAL  :: fn_val

! Local variables

REAL            :: u, sum
REAL, SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
  fn_val = u*sln
END IF

RETURN
END FUNCTION rnorm


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
    x = rnorm()
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

PROGRAM BAYESIANSEM
USE DataPaths
USE DataParameters
USE LU
USE Cholesky
USE Random 
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER, ALLOCATABLE, DIMENSION(:) ::      iwt
INTEGER ::                                 iter
INTEGER ::                                 check
REAL ::                                    start, finish
LOGICAL ::                                 First=.true.
!**********************************************************
! Bayesian estimation of Spatial Error Model              *
!  with Metropolis-Within-Gibbs Algorithm                 *
!**********************************************************
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

!**********************************************************
! Import data:                                            *
!**********************************************************
!(1) YXData.raw: Y and X variables.                           
!(2) weight_index.raw: Spatial unit's and its neighbor's IDs   

 input = adjustl(trim(path))//"YXData.raw"
 wt_input = adjustl(trim(path))//"weight_index.raw"

 OPEN (Unit=8, File=input, ACTION="READ")

! Note that X is read in row-by-row, which is not efficient
 DO i = 1, nobs
    READ (Unit = 8, Fmt = *, IOSTAT=ios) Y(i), X(i,2:KVars)
 ENDDO
 CLOSE (Unit = 8)
 print *, "IO status is ",ios

 OPEN (Unit = 8, File = wt_input, ACTION="READ")
 DO i = 1, nobs
   READ (8, *, IOSTAT=ios) id(i), weight_index(:,i)
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for weight data is ",ios

!****************************************************************
! Create files to save results                                  *
!****************************************************************
  output = adjustl(trim(path))//"EstimationResults"

  OPEN (Unit=10, File=output, ACTION="WRITE")
  OPEN (Unit=11, File="./AllDraws.raw", ACTION="WRITE")

! setup 
CALL SETUP

!****************************************************************
! Calculate weight matrix W as well as WY and WX.               *
! COUNT, PACK intrinsic functions are used.                     *              
!****************************************************************
  W = 0.0; WY = 0.0 ; WX = 0.0

  do i = 1, nobs
    if ( COUNT(weight_index(:,i) > 0) > 0 ) then
      allocate ( iwt( COUNT(weight_index(:,i) > 0)) )
      iwt = PACK(weight_index(:,i), weight_index(:,i) > 0)

      W(iwt,i) = 1.0         
      WY(i) = sum( Y(iwt) )  
      do j = 1,KVars
         WX(i,j) = sum( X(iwt,j) )
      end do
  
      deallocate (iwt)
     end if
   end do

!***********************************************************************
! MCMC Simulation: Metropolis-Within-Gibbs Algorithm                   *
!***********************************************************************
! Storage for Samples
   AllDraws = 0.0

! Set initial values
  call ols
  rho = 0.6

 print *, " "
 print *, "   ************************************************   "
 print *, "   * The number of draws is ", ndraw
 print *, "   * The average runtime is 397 sec. per 100 draws"
 print *, "   ************************************************   "
 print *, "Sampling starts ..."

 CALL CPU_TIME(start)
 CALL RANDOM_SEED

 iter = 1
 DO WHILE ( iter <= ndraw)
 
!**********************************************************************
!(1) Sample a vector of beta ~ N(beta1, M1^(-1))                      *
!**********************************************************************
  XTilde=X-rho*WX
  YTilde=Y-rho*WY
  do i=1, KVars
    do j=1, KVars
     XTX(i,j)=dot_product(XTilde(:,i), XTilde(:,j))
    end do
  XTY(i)=dot_product(XTilde(:,i),YTilde)
  end do

  M1 = M0 + (1.0/sigma_sq) * XTX
  beta_part = matmul(M0, beta0) + (1.0/sigma_sq) * XTY

  !invert M1
  call ludcmp(M1, KVars, index_beta, ninterchange, returncode)
  inverseM1 = identity_M1
  do i=1, KVars
   call lubksb(M1, KVars, index_beta, inverseM1(:,i))
  end do

  beta1 = matmul(inverseM1, beta_part) 

! generate random sample of beta
  call chodec(KVars,inverseM1,check)
  do i=1, KVars
    norm_rnd(i) = rnorm()
  end do

   beta = matmul(inverseM1, norm_rnd) +  beta1

!******************************************************************** 
!(2) Sample sigma2 ~ iG(v1/2, s1/2)                                 *
!********************************************************************
  v1 = nobs + v0 
  epsilon = YTilde - matmul(XTilde,beta)
  S1 = dot_product(epsilon, epsilon) + S0  
  sigma_sq = S1/(2.0 * random_gamma1(v1/2.0, First))

!*******************************************************************
!(3) Sample rho                                                    *
! Metropolis-Hastings algorithm                                    *
!*******************************************************************
  ! Generate rhostar
   rhostar = rho + cc * rnorm()
   accept = 0
   do while (accept == 0)
       if (( rhostar > rho_min).and. (rhostar < rho_max)) then
          accept = 1
       else
          rhostar = rho + cc *  rnorm()
       end if
   end do

  ! Calculate ratio
     rhotemp = rho
     call WEIGHT
     lnp = lnconditionalrho
     
     rhotemp = rhostar
     call WEIGHT
     lnpstar = lnconditionalrho

  !accept rhostar with prob. of ratio
     call random_number(rnd)

     ratio= exp(lnpstar - lnp)
     if (ratio < 1.0) then
         prob = ratio
       else
         prob = 1.0
       end if

     if (rnd < prob) then
        rho = rhostar
        acc = acc + 1
     end if
     acc_rate(iter) = real(acc)/real(iter)

   ! update cc based on std of rho draws
      if (acc_rate(iter) < 0.4) then
         cc = cc/1.1
      else if (acc_rate(iter) > 0.6) then
         cc = cc*1.1
      end if

! save samples
 if (iter > nburn) then
   AllDraws(iter-nburn, 1:KVars) = beta
   AllDraws(iter-nburn, KVars+1) = sqrt(sigma_sq)
   AllDraws(iter-nburn, npar) = rho
 end if

 iter = iter + 1
 END DO

   CALL CPU_TIME(finish)
   print *, "end of sampling: time", finish-start, "sec."

!***********************************************************************
! End of MCMC Simulation                                               *
!***********************************************************************

! Geweke's Convergence Diagnostic test
 frac1 = 0.1
 frac2 = 0.5
 CALL Geweke_Diagnostic(AllDraws, ndraw-nburn, npar, frac1, frac2, CD)

 write(10,*) ""
 write(10,*) "Geweke's Convergence Diagnostic test"
 write(10,*) CD


! Estimation Results 
CALL MC_Statistics(AllDraws, ndraw-nburn, npar, PostMeans, PostVars, NSE, RNE)

 write(10,*) ""
 write(10,*) "Post. Means, Post. Stand. Dev., NSE, RNE"
 do i = 1, npar 
  write(10,"(2X, 4F13.5)") PostMeans(i), sqrt(PostVars(i)), NSE(i), RNE(i)
 end do

      print *, "acceptance rate", iter-1, acc_rate(iter-1)

! Save draws for further analysis. 
  do i = 1, ndraw-nburn
    write(11, *) AllDraws(i,:) 
  end do

END PROGRAM BAYESIANSEM



SUBROUTINE SETUP
USE DataParameters
IMPLICIT none
INTEGER ::          i

! Add a ones column to the X matrix and create an identity matrix
!  to be used in finding the variance of beta
  X(:,1) = 1.0

  identity_M1 = 0.0
  do i = 1, KVars
    identity_M1(i,i) = 1.0
  enddo

! Set values of hyperparameters in prior distributions
! beta ~ N(beta0, M0^(-1))
  do i =1, KVars
    beta0(i) = 0.0
  end do

  M0 = 0.0
  do i=1, KVars
    M0(i,i) = 0.001
  end do

! sigma2 ~ iG(v0/2, s0/2)
  v0 = 0.0
  S0 = 0.0

! maximum and minimum values of rho
  rho_min = -1.0
  rho_max = 1.0

! Used in Metropolis-Hastings algorithm
  cc = 0.2
  acc = 0

RETURN
END SUBROUTINE SETUP



SUBROUTINE OLS
USE DataParameters
USE LU
IMPLICIT none
INTEGER::          i, j
!*****************************************************************
! Solve for OLS estimator of beta, using the LU decomposition of
! X'X and right-hand side vector X'Y. Make use of the columns of
! the X matrix and dot-products for efficiency.
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

RETURN
END SUBROUTINE OLS


SUBROUTINE WEIGHT 
USE DataParameters
USE LU
IMPLICIT none
INTEGER ::           i, j
REAL ::              logdet
!*****************************************************************
! Given a value of rho (rhotemp), it calculates the value of log *
! of the full conditional density for rho, which is used to      *
! calculate acceptance probability alpha                         *
!*****************************************************************

  YTilde = Y - rhotemp*WY
  XTilde = X - rhotemp*WX
  ITilde = 0.0
  do i=1, nobs; ITilde(i,i) = 1.0; enddo
  WHERE (W /= 0)                   
   ITilde = ITilde - rhotemp
  END WHERE

!*****************************************************************
! Get the log of the Jacobian for this value of rho              *
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
! Get the log of the full conditional of rho                     * 
!*****************************************************************
 epsilon = YTilde - matmul(XTilde, beta)
 
 lnconditionalrho = logdet- dot_product(epsilon, epsilon)/(2.0*sigma_sq)  

RETURN
END SUBROUTINE WEIGHT



SUBROUTINE Geweke_Diagnostic(Draws, ndraws, npars, frac1, frac2, CD )
USE Precision
IMPLICIT none
INTEGER, INTENT(IN)::               ndraws, npars
REAL(wp), INTENT(IN)::              frac1, frac2
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN):: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT):: CD
REAL(wp), DIMENSION(npars)::        PostMeans_a, PostMeans_b
REAL(wp), DIMENSION(npars)::        PostVars_a, PostVars_b
REAL(wp), DIMENSION(npars)::        NSE_a, NSE_b
REAL(wp), DIMENSION(npars)::        RNE_a, RNE_b
INTEGER::                           ndraws_a, ndraws_b
INTEGER::                           start_a, start_b
INTEGER::                           i
!********************************************************************
! Using MCMC draws, it returns Geweke (1992)'s convergence          * 
! diagonstic, CD.                                                   *
! Arguments                                                         *
!  Draws: (ndraws * npars) matrix of MCMC draws                     *
!  ndraws: the number of draws                                      *
!  npars: the number of parameters                                  *
!  frac1: first fraction of draws used for CD                       *
!  frac2: second fraction of draws used for CD                      * 
!  CD: convergence diagonstic statistic                             *
!********************************************************************

  ndraws_b = floor(frac2*real(ndraws))
  start_b =  ndraws - ndraws_b + 1
  CALL MC_Statistics(Draws(start_b:ndraws,:), ndraws_b, npars, &
                   PostMeans_b, PostVars_b, NSE_b, RNE_b)

  ndraws_a = floor(frac1*real(ndraws))
  start_a = 1 
  CALL MC_Statistics(Draws(start_a:ndraws_a,:), ndraws_a, npars, &
                   PostMeans_a, PostVars_a, NSE_a, RNE_a)

  CD = (PostMeans_a - PostMeans_b)/(NSE_a + NSE_b)

END SUBROUTINE



SUBROUTINE MC_Statistics(Draws, ndraws, npars, PostMeans, PostVars, NSE, RNE)
USE Precision
IMPLICIT none
INTEGER, INTENT(IN) ::                            ndraws, npars
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN) :: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT) ::        PostMeans, PostVars
REAL(wp), DIMENSION(npars), INTENT(OUT)::         NSE, RNE
REAL(wp), ALLOCATABLE, DIMENSION(:,:)::           Covars
REAL(wp), DIMENSION(npars) ::                     SG 
INTEGER ::                                        lags
INTEGER ::                                        i, j
!***************************************************************************
! Using MCMC draws, it returns posterior mean, variance, NSE (numerical    *
! Standard error, and RNE (relative numerical error).                      *
! Arguments                                                                *
!  Draws: (ndraws * npars) matrix of draws                                 *
!  ndraws: the number of draws                                             *
!  npars: the number of parameters                                         *
!  PostMeans: posterior means of the parameters                            *
!  PostVars: posterior variances of the parameters                         *
!  NSE: numerical standard error                                           *
!  RNE: relative numerical error                                           *
!***************************************************************************

! posterior means 
  PostMeans = SUM(Draws,1)/real(ndraws)

! posterior variances
  do i = 1, npars
    PostVars(i) = SUM( (Draws(:,i)**2) )/real(ndraws) - PostMeans(i)**2
  end do

! variance of posterior mean is SG/ndraws
! numerical standard error (nse) is sqrt(SG/ndraws)
! relative numerical error (rne) is PostVar/SG

  lags = floor( ndraws**(0.25_wp) ) + 1

  Allocate( Covars(lags, npars) )

   do i = 1, npars
    do j = 1, lags 
      Covars(j, i) = &
        dot_product( Draws(1:ndraws-j, i) - PostMeans(i), &
                     Draws(j+1:ndraws, i) - PostMeans(i) )/ &
        real(ndraws)
    end do
   end do  

   do i = 1, npars
    SG(i) = PostVars(i) + ( 2.0_wp*dot_product( (/(i,i=lags,1,-1)/), &
                       Covars(:,i)) )/ real(lags+1) 
   end do

   RNE = PostVars/SG 
   NSE = sqrt( SG/real(ndraws) ) 

  Deallocate ( Covars )

END SUBROUTINE 

