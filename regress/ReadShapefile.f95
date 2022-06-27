P
:q
Q!
OGRAM ReadShapeFile
IMPLICIT none
CHARACTER(80)                :: fname, Out_Coord_Name, Out_Check_Name
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

