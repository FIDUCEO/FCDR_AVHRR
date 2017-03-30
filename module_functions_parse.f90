MODULE module_functions_parse
 
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil
  use type
  use module_convert
 

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
! Get data type from value
  INTEGER FUNCTION Parse_Data_Type( Value )RESULT(type)

    INTEGER(GbcsInt1), INTENT(IN) :: Value

    type = IBITS(Value,4,3)

  END FUNCTION Parse_Data_Type

  ! Routine to get the date from POD guide data
  SUBROUTINE Parse_Date_Time( Header, Pos, year, dayno, hours, UTC_msecs )

    INTEGER(GbcsInt1), INTENT(IN) :: Header(:)
    INTEGER, INTENT(IN) :: Pos
    INTEGER, INTENT(OUT) :: year
    INTEGER, INTENT(OUT) :: dayno
    REAL, INTENT(OUT) :: hours
    INTEGER, INTENT(OUT) :: UTC_msecs

    INTEGER(GbcsInt2) :: DayNumber
    INTEGER(GbcsInt4) :: msecs
    INTEGER(GbcsInt1) :: ByteArray2(2)
    INTEGER(GbcsInt1) :: ByteArray4(4)

    ! First 7 bits are year
    year = IBITS(Header(Pos),1,7)
    if( year .lt. 70 )then
       year = year+2000
    else
       year = year+1900
    endif

    ! Next 9 bits are day number - get this by copying Header(1:2) to a 
    ! 2 byte integer and set to zero the top 7 bits
    IF( Little_Endian )THEN
       ! Byte swap
       ByteArray2(2) = Header(Pos)
       ByteArray2(1) = Header(Pos+1)
    ELSE
       ByteArray2 = Header(Pos:Pos+1)
    ENDIF

    ! Clear year entry
    DayNumber = TRANSFER(ByteArray2,DayNumber)
    DayNumber = IBCLR(DayNumber,15)
    DayNumber = IBCLR(DayNumber,14)
    DayNumber = IBCLR(DayNumber,13)
    DayNumber = IBCLR(DayNumber,12)
    DayNumber = IBCLR(DayNumber,11)
    DayNumber = IBCLR(DayNumber,10)
    DayNumber = IBCLR(DayNumber,9)

    DayNo = DayNumber

    ! Get last 4 bytes - time in mseconds
    IF( Little_Endian )THEN
       ! Byte swap
       ByteArray4(4) = Header(Pos+2)
       ByteArray4(3) = Header(Pos+3)
       ByteArray4(2) = Header(Pos+4)
       ByteArray4(1) = Header(Pos+5)
    ELSE
       ByteArray4 = Header(Pos+2:Pos+5)
    ENDIF
    msecs = TRANSFER(ByteArray4,msecs)
    UTC_msecs = msecs
    Hours = msecs*1e-3/3600. ! same as V2 and V3

  END SUBROUTINE Parse_Date_Time

  ! Routine to get the day and time from HRPT data i.e minor_frame(9:12)
  SUBROUTINE Parse_HRPT_Day_Time( frame, dayno, UTC_msecs )

    INTEGER(GbcsInt2), INTENT(IN) :: frame(4)
    INTEGER, INTENT(OUT) :: dayno
    INTEGER, INTENT(OUT) :: UTC_msecs

    INTEGER :: msecs(3)

    
    dayno = IBITS(frame(1),1,9)
    msecs(1) = IBITS(frame(2),0,7)
    msecs(2) = IBITS(frame(3),0,10)
    msecs(3) = IBITS(frame(4),0,10)
    UTC_msecs = msecs(1)*2**20 + msecs(2)*2**10 + msecs(3)

  END SUBROUTINE Parse_HRPT_Day_Time

  REAL(GbcsDble) FUNCTION Get_Date_Time( year, dayno, hours, yearstart )

    INTEGER, INTENT(IN) :: year
    INTEGER, INTENT(IN) :: dayno
    REAL, INTENT(IN) :: hours
    INTEGER, INTENT(IN) :: yearstart

    INTEGER :: i
    INTEGER :: totday
    
    totday = 0
    do i=yearstart,year-1
       IF( (0 .eq. MOD(i,4) .and. 0 .ne. MOD(i,100)) .or. &
            (0 .eq. MOD(i,400)) )then
          totday = totday + 366
       ELSE
          totday = totday + 365
       ENDIF
    end do
    totday = totday + dayno-1

    get_Date_time = (totday + hours/24.d0)*86400.d0 ! in seconds from start

  END FUNCTION Get_Date_Time
   
  

   


 

END MODULE module_functions_parse

