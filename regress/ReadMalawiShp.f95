PROGRAM ReadMalawiShp
USE ReadShapefile
IMPLICIT none
CHARACTER(80)                :: fname
INTEGER, ALLOCATABLE:: NumberBoundaryPoints(:,:)
REAL(8), ALLOCATABLE:: BoundaryPoints(:,:)
CHARACTER(80)                :: fmt   
INTEGER                      :: ios, i, istat
!***************************************************************** 
! Read a shapefile (.shp), binary file written in combination    *
! of little and big endian. Refer to the document below.         *
! http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf     *
!*****************************************************************

fname ='./mwi_pov'

CALL ReadSHP(fname, NumberBoundaryPoints, BoundaryPoints)

print *, size(NumberBoundaryPoints,1)
print *, NumberBoundaryPoints(1,:)

do i = 1, 10
  print *, BoundaryPoints(i,:)
end do

DEALLOCATE(NumberBoundaryPoints, BoundaryPoints)
END PROGRAM ReadMalawiShp
