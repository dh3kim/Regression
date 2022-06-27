!*************************************************************************
! Spatial Error Model (SEM): Regression Model with Spatially Correlated  *
! errors.                                                                *    
! Kelejian and Prucha (1999, IER)'s GM Approach                          *
!*************************************************************************
! Donghwan Kim                                                           *
!   Department of Economics                                              *
!   Stony Brook University                                               *
!*************************************************************************

MODULE DataParameters
 USE LA_PRECISION, ONLY:    WP => SP
 USE F95_LAPACK, ONLY:      LA_POSV, LA_SYSV
IMPLICIT none
 INTEGER, PARAMETER          ::  nobs = 376 
 INTEGER, PARAMETER          ::  KVars = 2
 INTEGER, PARAMETER          ::  n_neighbors = 13 ! Maximum number of neighbors
 INTEGER, DIMENSION(nobs)    ::  id
 REAL, DIMENSION(nobs)       ::  Y
 REAL, DIMENSION(nobs,KVars) ::  X
 INTEGER                     ::  weight_index(n_neighbors,nobs)

 ! Used in OLS
 REAL(wp), DIMENSION(KVars)         ::  XTY
 REAL(wp), DIMENSION(KVars, KVars)  ::  XTX
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
 REAL(wp), DIMENSION(n_param, n_param)  ::  H
 REAL(wp), DIMENSION(n_param, n_param)  ::  Hinverse
 REAL(wp), DIMENSION(n_param, n_param)  ::  identity_H
 REAL, DIMENSION(n_param)               ::  ddelta

 REAL                                   ::  rho

 ! Used in GLS
 REAL, DIMENSION(nobs, nobs)    ::   W
 REAL, DIMENSION(nobs)          ::   WY
 REAL, DIMENSION(nobs, KVars)   ::   WX

 REAL, DIMENSION(nobs)              ::  YTilde
 REAL, DIMENSION(nobs,KVars)        ::  XTilde

 REAL, DIMENSION(nobs)              ::  etilde
 REAL(wp), DIMENSION(KVars, KVars)  ::  XTXinverse 
 REAL, DIMENSION(KVars,KVars)       ::  var_beta

END MODULE DataParameters

PROGRAM SEMGM
USE DataParameters
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER ::                                 iter
INTEGER, ALLOCATABLE, DIMENSION(:) ::      iwt
REAL ::                                    SSE, newSSE

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

! read data with weight_index
 OPEN (Unit = 8, File ="./SEMTestFile.raw", ACTION="READ")
 do i = 1, nobs
   READ (Unit = 8, Fmt = *, IOSTAT=ios) id(i), Y(i), X(i,2:KVars), &
                                      weight_index(:,i)
 end do
 CLOSE (Unit = 8)
 print *, "IO status for Y, X data is ",ios


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
  delta(1) = 0.1       ! initial value of rho
  delta(2) = sigma_sq

  delta0 = delta

  call SSEfunc(SSE)

 iter = 12 

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

! Create an identity matrix to be used in GM
 identity_H = 0.0_wp
 do i = 1, n_param
   identity_H(i,i) = 1.0_wp
 end do

! XTX inverse
 XTXinverse = 0.0_wp
 do i = 1, KVars
   XTXinverse(i,i) = 1.0_wp
 end do

RETURN
END SUBROUTINE SETUP



SUBROUTINE OLS
USE DataParameters
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

  call la_posv(XTX,XTY)
  beta = XTY
  print *, "OLS estimators for beta is", beta
  
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
  
  Hinverse = identity_H
  call LA_SYSV(H, Hinverse)

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
IMPLICIT none
INTEGER::          i, j
REAL(wp), DIMENSION(KVars,KVars):: XX
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

  XX = XTX

  call LA_POSV(XTX,XTY)
  beta = XTY
  print *, "Final beta is ", beta

  etilde = YTilde - matmul(XTilde,beta)
  sigma_sq = dot_product(etilde,etilde)/real(nobs)
  print *, "Final sigma squared is ", sigma_sq

  ! Invert XTX
  call LA_POSV(XX, XTXinverse)

  ! Multiply the inverse by sigma-squared to obtain
  !  the variance matrix for beta
  var_beta = sigma_sq*XTXinverse

  do i = 1, KVars
   print *, beta(i), sqrt(var_beta(i,i)), beta(i)/sqrt(var_beta(i,i))
  enddo


RETURN
END SUBROUTINE GLS

