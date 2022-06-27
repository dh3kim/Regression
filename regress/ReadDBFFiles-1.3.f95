PROGRAM ReadDBFFiles
IMPLICIT none
CHARACTER(80) :: fname
LOGICAL       :: LittleEndian
BYTE          :: Version 
CHARACTER     :: InYear, InMonth, InDay
INTEGER       :: Year, Month, Day
CHARACTER     :: LNRecords(4)
INTEGER       :: NRecords
CHARACTER     :: LLenHeader(2), LLenEachRecord(2)
INTEGER       :: LenHeader, LenEachRecord
INTEGER       :: BytesForFieldDesc
INTEGER       :: NFields
BYTE          :: Ignorable(19)
BYTE          :: FieldNameA(12)
CHARACTER     :: FieldName(12)
BYTE          :: FieldTypeA
!CHARACTER     :: FieldType
INTEGER       :: Address
CHARACTER     :: SFieldLength, SDecimalCount
BYTE          :: FieldLength, DecimalCount
INTEGER       :: IFieldLength
INTEGER       :: Sum_FieldLength
BYTE          :: Unused(13), Unused_lasthead

CHARACTER(12), DIMENSION(:), ALLOCATABLE:: FieldsName, FieldsType
BYTE, DIMENSION(:), ALLOCATABLE:: FieldsLength
BYTE, DIMENSION(:), ALLOCATABLE:: Record

BYTE, DIMENSION(:), ALLOCATABLE:: Record
BYTE, DIMENSION(:), ALLOCATABLE:: Attribute

INTEGER:: Start

INTEGER   :: ios, i, j
!*****************************************************************
! Let us open a dbf file of shape file (the .dbf file            *
! is one of files consisting ESRI's shapefile)                   *
! NOTE: the .shp file is in the form of binary                   *
!*****************************************************************
fname ='./mwi_pov.dbf'
!fname = 'C:/data/mwi_pov.dbf'  ! Absoft version
open(3,file=fname,status='old',form='binary',iostat=ios)
if(ios /= 0) then
    write(*,*) 'Error opening file: ', fname
    stop
endif

!open(unit=6, file=OutName, action="write")
!*****************************************************************
!*****************************************************************
call endian(LittleEndian)

! Read the Main File Header.
! NOTE: the first 32 bytes are the Main Filder Header.

! (1) the first 1 bytes is version. 
   read(3) Version
! (2) Year, Month and Day (last updated) are stored in 1 bytes
!     with unsigned integer. Fortran 95 does not provide unsigned.
!     Read in with character and convert it into its corresponding
!     ASCII code (0 - 256), which becomes unsigned integer.
   read(3) InYear, InMonth, InDay 
   Year = ICHAR(InYear) + 1900 
   Month = ICHAR(InMonth)
   Day = ICHAR(InDay)
! (3) the next 4 bytes is number of records (with unsigned integer,
!      little endian). In big endian machine, switch byte order
!      and then read it with unsigned.
   read(3) LNRecords
   if (.not. LittleEndian) then
     CALL ByteOrderSwitch4(LNRecords)
   end if
  
   NRecords = 0
   do i = 1, 4
      NRecords = NRecords * 256 + ichar ( LNRecords(i) )
   end do

! (4) next 2 bytes are length of header structure and 
!     length of each record.
   read(3) LLenHeader, LLenEachRecord
   if (.not. LittleEndian) then
     CALL ByteOrderSwitch2(LLenHeader)
     CALL ByteOrderSwitch2(LLenEachRecord)
   end if

   LenHeader = 0; LenEachRecord = 0
   do i = 1, 2
     LenHeader = LenHeader * 256 + ichar ( LLenHeader(i) )
     LenEachRecord = LenEachRecord * 256 + ichar (LLenEachRecord(i))
  end do

  BytesForFieldDesc = LenHeader - 33
  NFields = BytesForFieldDesc/32  
   
!  (5) Next 20 bytes.
   read(3) Ignorable

print *,"   "
print *,"   Version: :", Version
print *,"   YYYYMMDD: ", Year, Month, Day
print *,"   NRecords: ", NRecords
print *,"   Length of Header: ", LenHeader
print *,"   Lenght of Each Record: ", LenEachRecord
print *,"   NFields: ", NFields

Sum_FieldLength = 0

! Field information.
! 32 bytes for each field.

ALLOCATE (FieldsName(NFields), FieldsType(NFields))
ALLOCATE (FieldsLength (NFields))

! FieldName and FieldType are stored in ASCII.
j = 1
Do While (BytesForFieldDesc > 0 )
   read(3) FieldNameA, FieldTypeA, Address, SFieldLength, SDecimalCount
   do i = 1, 11
     FieldsName(i)(j:j) = ACHAR( FieldNameA(i) )
   end do
   FieldsType(j) = ACHAR(FieldTypeA)
   FieldsLength(j) = ICHAR(SFieldLength)
   Sum_FieldLength = Sum_FieldLength + FieldsLength(j)
   DecimalCount = ICHAR(SDecimalCount)
   read(3) Unused

! For FieldType, see
!http://www.clicketyclick.dk/databases/xbase/format/data_types.html
!FieldName(:) = FieldsName(:)(j:j) 
print *,"   "
print *,"   *********************************"
print *,"   ", j
print *,"   FieldName: ", FieldsName(:)(j:j)
print *,"   FieldType: ", FieldsType(j)
print *,"   Address: ", Address
print *,"   FieldLength: ", FieldsLength(j)
print *,"   DeciamlCount: ", DecimalCount

j=j+1
BytesForFieldDesc = BytesForFieldDesc - 32

END DO

read(3) Unused_lasthead
    print *, Unused_lasthead
    print *, "Sum of Field lengths are ", Sum_Fieldlength

ALLOCATE (Record(Sum_Fieldlength))

i=1
DO WHILE (i <= NRecords)
 CALL FSEEK(3, 1, 1)
 !Data (contents of fields) are stored in ASCII.
 READ(3) Record(1:Sum_Fieldlength)
 j=1
 Start = 1
 DO WHILE (j <= NFields)

   ALLOCATE (Attribute (FieldsLength(j)) )
   Attribute = Record(Start:Start+FieldsLength(j)-1)

   if (FieldsType(j)=="C" .or. FieldsType(j)=="L") then
   
   ! [coding required here]
 
   end if

   ! [coding required here]

   if (i <= 2 .and. j ==1) then
     write(*,*) FieldsType(j)
     write (*,*) FieldsName(:)(j:j) 
     write (*,*) Attribute
   end if    

   DEALLOCATE (Attribute)
   Start = Start + FieldsLength(j)
   j = j + 1
   END DO

i=i+1 
END DO

DEALLOCATE (FieldsName, FieldsLength)
DEALLOCATE (Record)

END PROGRAM ReadDBFFiles 


SUBROUTINE ByteOrderSwitch2(ByteName)
IMPLICIT none
BYTE :: ByteName(2), revByteName(2)

   revByteName(1:1) = ByteName(2:2)
   revByteName(2:2) = ByteName(1:1)

   bytename = revbytename

RETURN
END SUBROUTINE ByteOrderSwitch2


SUBROUTINE ByteOrderSwitch4(ByteName)
IMPLICIT none
BYTE :: ByteName(4), revByteName(4)

   revByteName(1:1) = ByteName(4:4)
   revByteName(2:2) = ByteName(3:3)
   revByteName(3:3) = ByteName(2:2)
   revByteName(4:4) = ByteName(1:1)

   bytename = revbytename

RETURN
END SUBROUTINE ByteOrderSwitch4

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
