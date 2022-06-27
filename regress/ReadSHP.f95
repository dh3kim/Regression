MODULE ReadShapefile

CONTAINS

SUBROUTINE ReadSHP(fname, RecordsInfo, BoundaryPoints)
USE BYTEORDER
IMPLICIT none
CHARACTER(80), INTENT(IN)    :: fname
CHARACTER(80)                :: Out_Info_Name, Out_Coord_Name, Out_Check_Name
LOGICAL                      :: LittleEndian
INTEGER                      :: FileCode, FileLength, Version, ShapeType
INTEGER, DIMENSION(5)        :: UNUSED
CHARACTER(10)                :: ShapeTypeName
REAL(8)                      :: Xmin, Ymin, Xmax, Ymax
REAL(8)                      :: Zmin, Zmax, Mmin, Mmax                   
INTEGER                      :: NumBytes, StartBytes, Record
INTEGER                      :: RecordNumber, ContentLength
INTEGER                      :: NextRecord, ShapeTypeInRecord
REAL(8)                      :: XminInRecord, YminInRecord
REAL(8)                      :: XmaxInRecord, YmaxInRecord
INTEGER                      :: NumParts, NumPoints

INTEGER, ALLOCATABLE         :: PartStart(:), PartCode(:)
REAL(8)                      :: X_Points, Y_Points
INTEGER:: IndexFileLength, NumRecords
INTEGER, ALLOCATABLE, INTENT(OUT):: RecordsInfo(:,:)
REAL(8), ALLOCATABLE, INTENT(OUT):: BoundaryPoints(:,:)

! A derived type to save X and Y points which are unknown.
TYPE :: ID_X_Y_TYPE
  INTEGER :: ID, Part
  REAL(8) :: X, Y
  TYPE  (ID_X_Y_TYPE), POINTER :: p
END TYPE
TYPE (ID_X_Y_TYPE), POINTER :: head, tail, ptr

CHARACTER(80)                :: fmt   
INTEGER                      :: ios, i, istat
!***************************************************************** 
! Read a shapefile (.shp), binary file written in combination    *
! of little and big endian. Refer to the document below.         *
! http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf     *
!*****************************************************************

!*****************************************************************
! Before reading a shapefile, let us check if the machine is     *
! little endian or big endian. See the endian subroutine below   *
! about details on how to check it.                              *
!***************************************************************** 
call endian(LittleEndian)

!*****************************************************************
! Let us open a shape file (the .shp file                        *
! is one of files consisting ESRI's shapefile)                   *
! NOTE: the .shp file is in the form of binary                   *
!*****************************************************************
open(unit=3,file=trim(adjustl(fname))//'.shp',status='old',form='binary',iostat=ios)
if(ios /= 0) then
    write(*,*) 'Error opening file: ', fname
    stop
endif

open(unit=5,file=trim(adjustl(fname))//'.shx',status='old',form='binary',iostat=ios)
if(ios /= 0) then
    write(*,*) 'Error opening file: ', fname
    stop
endif

Out_Info_Name = './'//trim(adjustl(fname))//'_info.txt'
Out_Coord_Name = './'//trim(adjustl(fname))//'_coord.txt'
Out_Check_Name = './'//trim(adjustl(fname))//'_check.txt'
!OutName = 'C:/data/mwi_pov_out.txt'
open(unit=6, file=Out_Coord_Name, action="write")
open(unit=9, file=Out_Check_Name, action="write")
open(unit=12, file=Out_Info_Name, action="write")
!*****************************************************************
!*****************************************************************
! Read the Main File Header.             
! NOTE: the first 100 bytes are the Main Filder Header.

!*****************************************************************
! The first 28 bytes are stored as big endian. If the machine    *
! SUN SOLARIS (big endian), just read them as it is. If the      *
! machine is Intel (little endian), read them in byte with       *
! variable prefix B (meaning Big endian) and then switch byte    *
! order. By EQUIVALENCE, the original variable without the       *
! can be used.                                                   * 
!*****************************************************************
! (1) the first 24 bytes in which the first 4 bytes are 
!     FileCode (whose value is 9994) and the remaining bytes
!     are not used.
   read(3) FileCode, UNUSED
   if (LittleEndian) CALL ConvertByteOrder_int(FileCode)
 
! (2) The next 4 bytes are File Length. 
!   NOTE: The value for file length is the total length of the
!         file in 16-bit words (including the fifty 16-bit words,
!         100 bytes, that make up the header).
   read(3) FileLength
   if (LittleEndian) CALL ConvertByteOrder_int(FileLength)

!*****************************************************************
! The remaining 72 bytes are stored as little endian and thus    *
!*****************************************************************
! (3) The next 8 bytes in which the first 4 bytes are Version
!     (whose  value of 1000) and the next 4 bytes are Shape
!     Type.
   read(3) Version, ShapeType
   if ( .not. LittleEndian) then
     call ConvertByteOrder_int(Version)
     call ConvertByteOrder_int(ShapeType)
   end if


! The values for shape type are as follows:
SELECT CASE (ShapeType)
case(0)
  ShapeTypeName = "Null Shape"
case(1)
  ShapeTypeName = "Point"
case(3)
  ShapeTypeName = "PolyLine"
case(5)
  ShapeTypeName = "Polygon"
case(8)
  ShapeTypeName = "MultiPoint"
case(11)
  ShapeTypeName = "PointZ"
case(13)
  ShapeTypeName = "PolyLineZ"
case(15)
  ShapeTypeName = "PolygonZ"
case(18)
  ShapeTypeName = "MultiPointZ"
case(21)
  ShapeTypeName = "PointM"
case(23)
  ShapeTypeName = "PolyLineM"
case(25)
  ShapeTypeName = "PolygonM"
case(28)
  ShapeTypeName = "MultiPointM"
case(31)
  ShapeTypeName = "MultiPatch"
END SELECT

! (4) The next 64 bytes in which each 8 bytes represent 
!     bounding box (Xmin, Ymin, Xmax, Ymax, Zmin, Zmax,
!     Mmin, Mmax).
   read(3) Xmin, Ymin, Xmax, Ymax, Zmin, Zmax, Mmin, Mmax
   if (.not. LittleEndian) then
      call ConvertByteOrder_real(Xmin)
      call ConvertByteOrder_real(Ymin)
      call ConvertByteOrder_real(Xmax)
      call ConvertByteOrder_real(Ymax)
      call ConvertByteOrder_real(Zmin)
      call ConvertByteOrder_real(Zmax)
      call ConvertByteOrder_real(Mmin)
      call ConvertByteOrder_real(Mmax)
   end if

write(9,*) " "
write(9,*) "   ********************************************"
write(9,*) "   The Main File Header                       "
write(9,*) "   ********************************************"
write(9,*) "   FileCode:     ", FileCode
write(9,*) "   File Length:  ", FileLength, "(in 16-bits words)"
write(9,*) "   Version:      ", Version
write(9,*) "   Shape Type:   ", ShapeType, "(",  &
                                adjustl(trim(ShapeTypeName)), &
                                  ")" 
write(9,*) "   Xmin:         ", Xmin
write(9,*) "   Ymin:         ", Ymin
write(9,*) "   Xmax:         ", Xmax
write(9,*) "   Ymax:         ", Ymax
write(9,*) "   ********************************************"

!* End of the Main File Header Reading 
!***************************************************************
!***************************************************************


! Check the number of total bytes.
NumBytes = (FileLength * 16)/8


!***************************************************************
!***************************************************************
! In the case of polygon (ShapeType = 5)
! Read each record by loop

!*******************************************************************
!The file length stored in the index file header is the total length
!of the index file in 16-bit words (the fifty 16-bit words of the
! header plus 4 times the number of records.
!********************************************************************
CALL FSEEK(5, 24, 1)
read(5) IndexFileLength
if (LittleEndian) then
   call ConvertByteOrder_int(IndexFileLength)
end if
NumRecords = (IndexFileLength*2 - 100)/8

ALLOCATE (RecordsInfo(NumRecords,2))


StartBytes = 100
Record = 1

do while (StartBytes < NumBytes)

!**************************************************************** 
! (1) Each Record Header (big endian)                           *
! The first 8 bytes in which the first 4 bytes are record number*
! and the next 4 bytes are content length.                      *
! NOTE: The content length for a record is the length of the    *
!       record contents section measured in 16-bit words.       *
!       Each record, therefore, contributes (4 + content length)*
!       16-bit words toward the total length of the file, as    *
!       stored at Byte 24 in the file header.                   *
!****************************************************************
   read(3) RecordNumber, ContentLength
   if (LittleEndian) then
       call ConvertByteOrder_int(RecordNumber)
       call ConvertByteOrder_int(ContentLength)
   end if
   RecordsInfo(Record, 1) = RecordNumber

! Lengh of this record (in Bytes) 
NextRecord = ((ContentLength + 4) * 16)/8

!******************************************************************
! (2) Record Contents (little endian)                             *
!   (2-1) the first 4 bytes are shape type.                       *
!******************************************************************
   read(3) ShapeTypeInRecord
   if (.not. LittleEndian) call ConvertByteOrder_int(ShapeTypeInRecord)

!******************************************************************
! (2-2) Bounding Box (Xmin, Ymin, Xmax, Ymax)                     *
!****************************************************************** 
   read(3) XminInRecord, YminInRecord, XmaxInRecord, YmaxInRecord
   if (.not. LittleEndian) then 
      call ConvertByteOrder_real(XminInRecord)
      call ConvertByteOrder_real(YminInRecord)
      call ConvertByteOrder_real(XmaxInRecord)
      call ConvertByteOrder_real(YmaxInRecord)
   end if

!*******************************************************************
! (2-3) NumParts and NumPoints                                     *
!*******************************************************************
   read(3) NumParts, NumPoints
   if (.not. LittleEndian) then
       call ConvertByteOrder_int(NumParts)
       call ConvertByteOrder_int(NumPoints)
   end if
   RecordsInfo(Record,2) = NumPoints

!*******************************************************************
! (2-4) Parts:                                                     *
!*******************************************************************
ALLOCATE (PartStart(NumParts), PartCode (NumPoints), STAT = istat)
 DO i = 1, NumParts
   read(3) PartStart(i)
    if (.not. LittleEndian)  call ConvertByteOrder_int(PartStart(i))
 END DO

DO i = 1, NumParts 
  if (i < NumParts) then
      PartCode(PartStart(i)+1: PartStart(i+1)) = i
  else if (i == NumParts) then 
      PartCode(PartStart(i)+1: NumPoints) = i
  end if
END DO 

!*******************************************************************
! (2-5) Points                                                     *
!     The first 8 bytes are X coordinates of the first point, the  *
!     next 8 bytes are Y coordinates of the first point, the next  *
!     8 bytes are X coordinates of the sencond point, and so on.   *
!*******************************************************************

i = 1

do while (i <= NumPoints)
      read(3) X_Points, Y_Points
      if (.not. LittleEndian) then
         call ConvertByteOrder_real(X_Points)
         call ConvertByteOrder_real(Y_Points)
      end if
   
   ! Linked list to save   
   IF (.NOT. ASSOCIATED(head) ) THEN
      ALLOCATE (Head, STAT = istat)
      tail => head
      NULLIFY (tail%p)
      tail%ID = RecordNumber
      tail%Part = PartCode(i)
      tail%X = X_Points
      tail%Y = Y_Points
   ELSE
      ALLOCATE (tail%p, STAT = istat)
      tail => tail%p
      NULLIFY (tail%p)
      tail%ID = RecordNumber
      tail%Part = PartCode(i)
      tail%X = X_Points
      tail%Y = Y_Points
   END IF
   
i = i + 1
end do

  write(9,*) "  "
  write(9,*) "   **************************************************"
  write(9,*) "    Record Number is ", RecordNumber   
  write(9,*) "   **************************************************"          
  write(9,*) "    Content Length (in 16-bit words) is ", ContentLength
  write(9,*) "    ShapeType is", ShapeTypeInRecord
  write(9,*) "    Xmin is ", XminInRecord
  write(9,*) "    Ymin is ", YminInRecord
  write(9,*) "    Xmax is ", XmaxInRecord
  write(9,*) "    Ymax is ", YmaxInRecord
  write(9,*) "    Number of Parts is ", NumParts, " ,Starting at ", PartStart
  write(9,*) "    Number of Points is ", NumPoints

DEALLOCATE(PartStart, PartCode)

StartBytes = StartBytes + NextRecord
Record = Record + 1

end do

ALLOCATE (BoundaryPoints (SUM(RecordsInfo(:,2)),2))

! Write out the data.
! Note that list-directed output (write (10,*)) is unreliable
! in that it is a compiler-specific decision how and when to
! break lines of output. Use formatted WRITE below.
fmt = '(2I3, 2F12.6)'
ptr => head
i = 1
Output: DO
  IF (.NOT. ASSOCIATED(ptr) ) EXIT
  WRITE(6,fmt) ptr%ID, ptr%Part, ptr%X, ptr%Y
  BoundaryPoints(i,1) = ptr%X
  BoundaryPoints(i,2) = ptr%Y
  ptr => ptr%p
  i = i+1
END DO Output


!DEALLOCATE (RecordsInfo)

RETURN
END SUBROUTINE ReadSHP


END MODULE ReadShapefile
