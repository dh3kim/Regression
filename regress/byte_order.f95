MODULE BYTEORDER

CONTAINS

SUBROUTINE ConvertByteOrder_int(X)
IMPLICIT none
INTEGER, INTENT(INOUT) :: X
INTEGER :: Xtemp
BYTE :: Xbyte(4)
! NOTE that dummy-arguments can not be an
! argument of the EQUIVALENCE statement.
EQUIVALENCE(Xtemp, Xbyte)

Xtemp = X
CALL ByteOrderSwitch(Xbyte)
X = Xtemp

RETURN
END SUBROUTINE ConvertByteOrder_int


SUBROUTINE ConvertByteOrder_real(X)
REAL(8), INTENT(INOUT) :: X
REAL(8) :: Xtemp
BYTE :: Xbyte(8)
EQUIVALENCE(Xtemp, Xbyte)

Xtemp = X
CALL ByteOrderSwitch(Xbyte)
X = Xtemp

RETURN
END SUBROUTINE ConvertByteOrder_real


Subroutine ByteOrderSwitch(X)
Implicit none
BYTE, INTENT(INOUT) :: X(:)
BYTE:: temp
INTEGER :: i, n
n = size(X)

i = 1
do while (i <= n/2)

  temp = X(i)
  X(i) = X(n-i+1)
  X(n-i+1) = temp

  i = i + 1
end do

RETURN
END SUBROUTINE ByteOrderSwitch

subroutine endian(litend)
implicit none
! checks if this is a little endian machine
! returns litend=.true. if it is, litend=.false. if not
! Source: http://www.atmos.washington.edu/~salathe/osx_unix/endian.html
!***********************************************************************
! Description by dokim:
! In the code below, j is a 1-byte (8-bits) integer array of
! dimension 2 and i is a 2-bytes (16-bits) integer of dimension 1.
! By the "equivalence" statement, i and j share same memory spaces
! starting at same memory address (CPU assigns a memory address to
! each byte). By assinging 1 to i, the integer 1 (in decimal) is
! is stored in the variable i's memory address as binary. Note that,
! trivially (?), integer 1 (in decimal) is 0000000000000001 (in binary
! with 16-bits) or 00000001 (in binary with 8-bits).

! If machine is big endian, 1 is stored as 00000000 00000001 at
! i's memory addresses where left 1-byte (8-bits) is the first
! starting (lowest) momory address of i. However, if machine is
! is little endian, 1 is stored as 00000001 00000000.

! So if j(1) is 00000001 (in binary) or 1 (in decimal), the machine
! is little endian. if j(1) is 00000000 (in binary), the machine is
! big endian.
!***********************************************************************

      integer*1 j(2)
      integer*2 i
      equivalence (i,j)
      logical litend

      i = 1
      if (j(1).eq.1) then
         litend = .true.
      else
         litend = .false.
      end if

RETURN
end subroutine endian

END MODULE BYTEORDER
