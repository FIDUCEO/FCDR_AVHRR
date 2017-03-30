MODULE module_convert
 
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil
  use type
 

#ifdef USE_HDF5
  USE HDF5
#endif

#ifdef USE_GZFILE
  USE gzfile
#endif

#ifdef USE_IEEE
  USE ieee_arithmetic
#endif 

implicit none
PUBLIC

CONTAINS
 SUBROUTINE Convert_Int1_String( Pos, Length, Array, String )

    INTEGER, INTENT(IN) :: Pos
    INTEGER, INTENT(IN) :: Length
    INTEGER(GbcsInt1), INTENT(IN) :: Array(:)
    CHARACTER(LEN=*), INTENT(OUT) :: String

    INTEGER :: I

    DO I=1,Length
       String(I:I) = CHAR(Array(Pos+I))
    END DO

  END SUBROUTINE Convert_Int1_String

  ! Swap Byte order
  INTEGER(GbcsInt2) FUNCTION Header2Short( PosIn, Data )

    INTEGER, INTENT(IN) :: PosIn
    INTEGER(GbcsInt1), INTENT(IN) :: Data(:)

    INTEGER :: Pos
    INTEGER(GbcsInt1) :: DataArray(2)
    INTEGER(GbcsInt2) :: Data2 = 0

    IF( Little_Endian )THEN

       Pos = PosIn+1
       DataArray(1) = Data(Pos+1)
       DataArray(2) = Data(Pos)

       header2Short = TRANSFER(DataArray,Data2)

    ELSE
       
       Pos = PosIn+1
       DataArray(1) = Data(Pos)
       DataArray(2) = Data(Pos+1)

       header2Short = TRANSFER(DataArray,Data2)


    ENDIF

  END FUNCTION Header2Short

  ! Swap Byte order
  INTEGER FUNCTION Header2UShort( PosIn, Data )

    INTEGER, INTENT(IN) :: PosIn
    INTEGER(GbcsInt1), INTENT(IN) :: Data(:)

    INTEGER :: Pos
    INTEGER(GbcsInt1) :: DataArray(2)
    INTEGER(GbcsInt2) :: Data2 = 0

    IF( Little_Endian )THEN
       Pos = PosIn+1
       DataArray(1) = Data(Pos+1)
       DataArray(2) = Data(Pos)
       
       header2UShort = TRANSFER(DataArray,Data2)
    ELSE
       Pos = PosIn+1
       DataArray(1) = Data(Pos)
       DataArray(2) = Data(Pos+1)
       
       header2UShort = TRANSFER(DataArray,Data2)
    ENDIF

    IF( 0 .gt. header2UShort )THEN
       ! remove - sign
       header2ushort = header2ushort + 65536
    ENDIF

  END FUNCTION Header2UShort

  ! Swap Byte order
  INTEGER FUNCTION Header2Int( PosIn, Data )

    INTEGER, INTENT(IN) :: PosIn
    INTEGER(GbcsInt1), INTENT(IN) :: Data(:)

    INTEGER :: Pos
    INTEGER(GbcsInt1) :: DataArray(4)
    INTEGER(GbcsInt4) :: Data4 = 0

    IF( Little_Endian )THEN
       Pos = PosIn+1
       DataArray(1) = Data(Pos+3)
       DataArray(2) = Data(Pos+2)
       DataArray(3) = Data(Pos+1)
       DataArray(4) = Data(Pos)

       Header2Int = TRANSFER(DataArray,Data4)
    ELSE
       Pos = PosIn+1
       DataArray(1) = Data(Pos)
       DataArray(2) = Data(Pos+1)
       DataArray(3) = Data(Pos+2)
       DataArray(4) = Data(Pos+3)

       Header2Int = TRANSFER(DataArray,Data4)
    ENDIF

  END FUNCTION Header2Int

  
  
  
  

   


 

END MODULE module_convert

