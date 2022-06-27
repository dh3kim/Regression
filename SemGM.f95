! Note:
! This program is originally written with LAPACK (Linear Algebra PACKage)
! but edited for eco522 class. This version uses functions in the LU
! module (see below).  
MODULE DataParameters
IMPLICIT none
 INTEGER, PARAMETER          ::  nobs = 376 
 INTEGER, PARAMETER          ::  KVars = 2
 INTEGER, PARAMETER          ::  n_neighbors = 13 ! Maximum number of neighbors
 INTEGER, DIMENSION(nobs)    ::  id
 REAL, DIMENSION(nobs)       ::  Y
 REAL, DIMENSION(nobs,KVars) ::  X
 INTEGER                     ::  weight_index(n_neighbors,nobs)

 ! Used in OLS
 REAL, DIMENSION(KVars)             ::  XTY
 REAL, DIMENSION(KVars, KVars)      ::  XTX, identity_beta
 REAL, DIMENSION(KVars)             ::  beta
 REAL, DIMENSION(nobs)              ::  e  ! OLS residuals
 REAL                               ::  sigma_sq 

 ! Used in GM
 INTEGER, PARAMETER                 ::  n_moments = 3
 INTEGER, PARAMETER                 ::  n_param = 2
 INTEGER                            ::  num_neighbor
 REAL, DIMENSION(nobs)              ::  We, WWe
 REAL                               ::  trWTW

 REAL, DIMENSION(n_moments, n_moments)  ::  Gn
 REAL, DIMENSION(n_moments)             ::  small_gn

 REAL, DIMENSION(n_param)               :: delta, delta0
 REAL, DIMENSION(n_moments)             :: Error, SE

 REAL, DIMENSION(n_param)               ::  g
 REAL, DIMENSION(n_param, n_param)      ::  H
 REAL, DIMENSION(n_param, n_param)      ::  Hinverse
 REAL, DIMENSION(n_param, n_param)      ::  identity_H
 REAL, DIMENSION(n_param)               ::  ddelta

 REAL                                   ::  rho

 ! Used in GLS
 REAL, DIMENSION(nobs, nobs)    ::   W
 REAL, DIMENSION(nobs)          ::   WY
 REAL, DIMENSION(nobs, KVars)   ::   WX

 REAL, DIMENSION(nobs)              ::  YTilde
 REAL, DIMENSION(nobs,KVars)        ::  XTilde

 REAL, DIMENSION(nobs)              ::  etilde
 REAL, DIMENSION(KVars,KVars)       ::  var_beta

! Used in LU factorization 
INTEGER ::                       returncode, ninterchange
INTEGER ::                       index_beta(Kvars), index_H(n_param)

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



PROGRAM SEMGM
USE DataParameters
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER ::                                 iter
INTEGER, ALLOCATABLE, DIMENSION(:) ::      iwt
REAL ::                                    SSE, newSSE
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

 CALL SETUP

!*****************************************************************
! Get OLS  results                                               *
!*****************************************************************
 CALL OLS


!*****************************************************************
! Calculate weight matrix W as well as WY and WX.                *
!*****************************************************************
! GM: Calculate We, WWe, and tr(WTW) used in sample moments.     *
!     tr(WTW) method used only in symmetric binary matrix.       *
!*****************************************************************
 W = 0.0; WY = 0.0; WX = 0.0
 We = 0.0; WWe = 0.0; trWTW = 0.0

 DO i = 1, nobs
   num_neighbor =  COUNT(weight_index(:,i) > 0 )

   if ( num_neighbor > 0 ) then
      allocate ( iwt( num_neighbor ) )
      iwt = PACK(weight_index(:,i), weight_index(:,i) > 0)

      W(iwt,i) = 1.0         ! exploit symmetry of weights
      WY(i) = sum( Y(iwt) )  ! exploit binary nature of weights
      do j = 1,KVars
        WX(i,j) = sum( X(iwt,j) )
      end do
 
      We(i) = sum ( e(iwt) )
      WWe(i) = sum (We(iwt) ) 
      ! trace(WTW) when a symmetrix binary weight matrix
      trWTW = trWTW + num_neighbor 

      deallocate (iwt)
   end if
 END DO 

!****************************************************
! call elements of sample moments                   *
!****************************************************
 CALL MOMENTS


!****************************************************
! Find GM estimators of rho and sigma_sq            *
! Nonlinear Regression: Newton's method             *
!****************************************************
! initial guesses
  delta(1) = 0.3       ! initial value of rho
  delta(2) = sigma_sq

  delta0 = delta

  call SSEfunc(SSE)

 iter = 15
 DO i = 1, iter
   print *, 'iteration', i
   call direction 

   delta0 = delta + ddelta
   call SSEfunc(newSSE)
   if (newSSE.lt.SSE) then
     delta = delta + ddelta
     print *, 'newSSE = ', newSSE 
     SSE = newSSE
   else
     print *, 'adjusting step length'
     call step(SSE, newSSE)
     delta = delta + ddelta
     SSE = newSSE
   endif

  END DO

  print *, 'GM estimator of rho', delta(1)
  print *, 'GM estimator of sigma_sq', delta(2)

!***************************************************************
! Get GLS                                                      *
!***************************************************************
 rho = delta(1)
 CALL GLS
 
END PROGRAM SEMGM



SUBROUTINE SETUP
USE DataParameters
IMPLICIT none
INTEGER ::          i

! Add a ones column to the X matrix and create an identity matrix
!  to be used in finding the variance of beta
 X(:,1) = 1.0

! Create an identity matrix to be used in 
! inverse calculation with LU decomposition
identity_beta = 0.0
 do i = 1, KVars
  identity_beta(i,i) = 1.0
 enddo
 
 identity_H = 0.0
 do i = 1, n_param
   identity_H(i,i) = 1.0
 end do



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
 
  e = Y - matmul(X, beta)
  sigma_sq = dot_product(e, e)/real(nobs)
  print *, "OLS sigma squared is ", sigma_sq

RETURN
END SUBROUTINE OLS




SUBROUTINE MOMENTS
USE DataParameters
IMPLICIT none

! Elements of samples moments
Gn = 0.0
Gn(1,1) = (2.0/real(nobs))*dot_product(e, We)
Gn(1,2) = -(1.0/real(nobs))*dot_product(We, We)
Gn(1,3) = 1.0
Gn(2,1) = (2.0/real(nobs))*dot_product(We, WWe)
Gn(2,2) = -(1.0/real(nobs))*dot_product(WWe, WWe)
Gn(2,3) = (1.0/real(nobs))*trWTW
Gn(3,1) = (1.0/real(nobs))*(dot_product(e, WWe) + dot_product(We, We))
Gn(3,2) = -(1.0/real(nobs))*dot_product(We, WWe)

small_gn = 0.0
small_gn(1) = (1.0/real(nobs))*dot_product(e, e)
small_gn(2) = (1.0/real(nobs))*dot_product(We, We)
small_gn(3) = (1.0/real(nobs))*dot_product(e, We)

RETURN
END SUBROUTINE MOMENTS 



SUBROUTINE SSEfunc(SSE)
USE DataParameters 
IMPLICIT none
INTEGER ::  i
REAL, INTENT(OUT) ::     SSE

! Calculate SSE (Sum of Squared Error)
 do i = 1, n_moments 
    Error(i) =Gn(i,1)*delta0(1) + Gn(i,2)*delta0(1)**2.0 +&
                          Gn(i,3)*delta0(2)-small_gn(i)
    SE(i) = Error(i)**2.0
 end do

 SSE = sum(SE)

RETURN
END SUBROUTINE SSEfunc



SUBROUTINE direction 
USE DataParameters
USE LU
IMPLICIT none
INTEGER ::  i

 ! gradient g
   g = 0.0
   do i = 1, n_moments 
     g(1) = g(1) + 2.0*Error(i)*(Gn(i,1)+2.0*delta0(1)*Gn(i,2))
     g(2) = g(2) + 2.0*Error(i)*Gn(i,3)
   end do

 ! hessian H
   H = 0.0
   do i = 1, n_moments 
     H(1,1) = H(1,1) + 2.0*(Gn(i,1)+2.0*delta0(1)*Gn(i,2))**2.0 &
                 + 4.0*Error(i)*Gn(i,2)
     H(1,2) = H(1,2) + 2.0*(Gn(i,1)+2.0*delta0(1)*Gn(i,2))*Gn(i,3)
     H(2,1) = H(2,1) + 2.0*(Gn(i,1)+2.0*delta0(1)*Gn(i,2))*Gn(i,3)
     H(2,2) = H(2,2) + 2.0*Gn(i,3)**2.0
   end do
  
 ! inverse of H 
  Hinverse = identity_H
  
  call ludcmp(H,n_param,index_H,ninterchange,returncode)
 
  if (returncode == 0) then
    do i = 1, n_param
      call lubksb(H,n_param,index_H,Hinverse(:,i))
    end do
  else
    print *, "Problem with inverse of Hessian",returncode
  endif
  
  ddelta =-matmul(Hinverse, g)

RETURN
END SUBROUTINE direction 



SUBROUTINE step(SSE, newSSE)
USE DataParameters
IMPLICIT none
REAL, INTENT(IN)    :: SSE
REAL, INTENT(INOUT) :: newSSE

DO WHILE(newSSE < SSE-1d-8)
   print *, 'SSE and last newSSE', SSE, newSSE

   ddelta = 0.8d0*ddelta
   delta0 = delta + ddelta
   call SSEfunc(newSSE)
END DO


RETURN
END SUBROUTINE step



SUBROUTINE GLS
USE DataParameters
USE LU
IMPLICIT none
INTEGER::          i, j
REAL, DIMENSION(KVars,KVars):: XX
!*****************************************************************
! Solve for GLS estimator of beta, using final rho value
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
 else
   print *, "Problem with final beta ",returncode
 endif
  

  etilde = YTilde - matmul(XTilde,beta)
  sigma_sq = dot_product(etilde,etilde)/real(nobs)
  print *, "ML estimator of sigma_sq is ", sigma_sq

  ! Invert XTX 
  do i = 1, KVars
     call lubksb(XTX,KVars,index_beta,identity_beta(:,i))
  enddo

  ! Multiply the inverse by sigma-squared to obtain
  !  the variance matrix for beta
  var_beta = sigma_sq*identity_beta


  print *, "beta, stand. dev., t-statistic"
  do i = 1, KVars
   print *, beta(i), sqrt(var_beta(i,i)), beta(i)/sqrt(var_beta(i,i))
  enddo


RETURN
END SUBROUTINE GLS

