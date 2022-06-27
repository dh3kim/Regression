PROGRAM Malawi 
! A line with exclamation mark (!) is just a comment line.
! Start with a statement IMPLICIT NONE and declare all variables
! used in the program. Note that Fortran is case-insensitive.
!##################################################################
! Variable Declaration.                                           #
!##################################################################
USE OLS 
IMPLICIT none
INTEGER, PARAMETER ::            nobs = 376
INTEGER, PARAMETER ::            KVars = 2
REAL, DIMENSION(nobs) ::         Y
REAL, DIMENSION(nobs,KVars) ::   X
REAL, DIMENSION(KVars) ::        beta
REAL, DIMENSION(KVars,KVars) ::  var_beta
REAL, DIMENSION(nobs) ::         e
REAL ::                          sigma_sq
INTEGER ::                       i, ios

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

!*****************************************************************
! Create a OLSResiduals.raw file for saving OLS residuals later. *
!*****************************************************************
 OPEN (Unit = 10, File="./OLSResiduals.raw", ACTION="WRITE")

!*****************************************************************
! Add a ones column to the X matrix and call OLS_LU              *
!*****************************************************************
 X(:,1) = 1.0

 CALL OLS_LA(Y, X, nobs, KVars, beta, sigma_sq, var_beta, e)

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

END PROGRAM Malawi 

