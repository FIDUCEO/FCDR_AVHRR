MODULE module_functions
 
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
SUBROUTINE greatCircleDistance(lon1, lat1, lon2, lat2, d)
    REAL, INTENT(in) :: lon1, lat1, lon2, lat2
    REAL, INTENT(out) :: d
    REAL :: r
    REAL(GbcsDble) :: arg
    
    r = 6400.0
    arg = SIN(D_Deg2Rad * DBLE(lat1)) * SIN(D_Deg2Rad * DBLE(lat2)) &
         + COS(D_Deg2Rad * DBLE(lat1)) * COS(D_Deg2Rad * DBLE(lat2)) * COS(D_Deg2Rad * (DBLE(lon2) - DBLE(lon1)))
    IF (arg.LT.-1.0d0) THEN
       arg = -1.0d0
    ELSE IF (arg.GT.1.0d0) THEN
       arg = 1.0d0
    END IF
    d = r * ACOS(arg)
    
  END SUBROUTINE greatCircleDistance

 
  ! Check I/O Status
  SUBROUTINE Check_IOS( Ios, String1, String2, String3, AddString, String4 )

    INTEGER, INTENT(IN) :: Ios
    CHARACTER(LEN=*), INTENT(IN) :: String1
    CHARACTER(LEN=*), INTENT(IN) :: String2
    CHARACTER(LEN=*), INTENT(IN) :: String3
    LOGICAL, INTENT(IN) :: AddString
    CHARACTER(LEN=*), INTENT(IN) :: String4

    IF( 0 /= Ios )THEN
       CALL Out_ERROR( String1, String2, String3, AddString, String4 )
    ENDIF

  END SUBROUTINE Check_IOS

  ! Check I/O Status
  SUBROUTINE Out_ERROR( String1, String2, String3, AddString, String4 )

    CHARACTER(LEN=*), INTENT(IN) :: String1
    CHARACTER(LEN=*), INTENT(IN) :: String2
    CHARACTER(LEN=*), INTENT(IN) :: String3
    LOGICAL, INTENT(IN) :: AddString
    CHARACTER(LEN=*), INTENT(IN) :: String4

    IF( AddString )THEN
       write(*,*)' ERROR: '//TRIM(String1)//' '//TRIM(String4)
    ELSE
       write(*,*)' ERROR: '//TRIM(String1)
    ENDIF
    write(*,*)'         SUBROUTINE '//TRIM(String2)
    write(*,*)'         CODE - '//TRIM(String3)
    STOP 1

  END SUBROUTINE Out_ERROR

  ! Check I/O Status
  SUBROUTINE Out_WARNING( String1, String2, String3, AddString, String4 )

    CHARACTER(LEN=*), INTENT(IN) :: String1
    CHARACTER(LEN=*), INTENT(IN) :: String2
    CHARACTER(LEN=*), INTENT(IN) :: String3
    LOGICAL, INTENT(IN) :: AddString
    CHARACTER(LEN=*), INTENT(IN) :: String4

    IF( AddString )THEN
       write(*,*)' WARNING: '//TRIM(String1)//' '//TRIM(String4)
    ELSE
       write(*,*)' WARNING: '//TRIM(String1)
    ENDIF
    write(*,*)'         SUBROUTINE '//TRIM(String2)
    write(*,*)'         CODE - '//TRIM(String3)

  END SUBROUTINE Out_WARNING
 
  FUNCTION Big_Endian()

    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER( GbcsInt2 ) :: Source = 1_GbcsInt2

    ! ------------
    ! The function
    ! ------------

    LOGICAL :: Big_Endian

    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC TRANSFER, ICHAR

    ! ----------------------------------
    ! Initialise result to little-endian
    ! ----------------------------------

    Big_Endian = .FALSE.


    ! ------------------------------------------------------------
    ! Test for "endian-ness".
    !
    ! TRANSFER( source, 'a' ) returns a result with the physical
    !   representation of the number 1, i.e. an integer, but
    !   interpreted as a character (the type of 'a' - a character,
    !   not the value, is what is important).
    !
    ! IACHAR returns the position of a character in the ASCII
    !   collating sequence associated with the kind type parameter
    !   of the character.
    ! ------------------------------------------------------------

    IF ( IACHAR( TRANSFER( Source, 'a' ) ) == 0 ) Big_Endian = .TRUE.

  END FUNCTION Big_Endian

!------------------------------------------------------------------------------
!F+
! NAME:
!       AVHRR_Date
!
! PURPOSE:
!       Convert an AVHRR year,daynum,msec triplet to a Gbcs DateTime structure
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       date = AVHRR_Date(year, daynum, msec)
!
! INPUT ARGUMENTS:
!       year   - integer
!       daynum - integer
!       msec   - integer
!
! RESULT:
!       date
!
! CALLS:
!       Day_Number_To_Date
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! CREATION HISTORY:
!       Written by:   Owen Embury 03/10/2011
!                     IAES, University of Edinburgh
!F-
!------------------------------------------------------------------------------
  ELEMENTAL FUNCTION AVHRR_Date(year, daynum, msec) RESULT(date)
    USE GbcsDateTime
    IMPLICIT NONE
    
    ! ---------
    ! Arguments
    ! ---------
    INTEGER, INTENT(IN) :: year
    INTEGER, INTENT(IN) :: daynum
    INTEGER, INTENT(IN) :: msec

    ! ---------
    ! Function result
    ! ---------
    TYPE(DateTime) :: date

    IF (year.EQ.NAN_I .OR. daynum.EQ.NAN_I .OR. msec.EQ.NAN_I) THEN
       date%Year = NAN_I
       date%Month = NAN_I
       date%Day = NAN_I
       date%utc_offset = NAN_I
       date%Hour = NAN_I
       date%Minute = NAN_I
       date%Seconds = NAN_I
       date%Sec1000 = NAN_I
    ELSE 
       date%Year = year
       CALL Day_Number_To_Date(year, daynum, date%Month, date%Day)
       date%utc_offset = 0
       date%Hour    = msec/3600000
       date%Minute  = MODULO(msec/60000, 60)
       date%Seconds = MODULO(msec/1000,  60)
       date%Sec1000 = MODULO(msec, 1000)
    END IF

  END FUNCTION AVHRR_Date

!------------------------------------------------------------------------------
!F+
! NAME:
!       Invert_AVHRR_Date
!
! PURPOSE:
!       Convert a Gbcs DateTime structure to an AVHRR year,daynum,msec triplet
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL Invert_AVHRR_Date(date, year, daynum, msec)
!
! INPUT ARGUMENTS:
!       date - GBCS DateTime
!
! OUTPUT ARGUMENTS:
!       year   - integer
!       daynum - integer
!       msec   - integer
!
! CALLS:
!       Day_Of_Year_int
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! CREATION HISTORY:
!       Written by:   Mark Filipiak 06/10/2011
!                     IAES, University of Edinburgh
!F-
!------------------------------------------------------------------------------
  SUBROUTINE Invert_AVHRR_Date(date, year, daynum, msec)
    USE GbcsDateTime
    IMPLICIT NONE
    
    ! ---------
    ! Arguments
    ! ---------
    TYPE(DateTime), INTENT(IN) :: date
    INTEGER, INTENT(OUT) :: year
    INTEGER, INTENT(OUT) :: daynum
    INTEGER, INTENT(OUT) :: msec

    IF (date%Year.EQ.NAN_I .OR. date%Month.EQ.NAN_I .OR. date%Day.EQ.NAN_I &
         .OR. date%utc_offset.EQ.NAN_I &
         .OR. date%Hour.EQ.NAN_I .OR. date%Minute.EQ.NAN_I &
         .OR. date%Seconds.EQ.NAN_I .OR. date%Sec1000.EQ.NAN_I ) THEN
       year = NAN_I
       daynum = NAN_I
       msec = NAN_I
    ELSE
       year = date%Year
       daynum = Day_Of_Year_int(date%Year, date%Month, date%Day)
       ! AVHRR times are UTC: date%utc_offset == 0
       ! This subroutine must be an exact inverse of AVHRR_Date
       msec = date%Hour*3600000 + date%Minute*60000 &
            + date%Seconds*1000 + date%Sec1000
    END IF

  END SUBROUTINE Invert_AVHRR_Date

  
  
  
  

   


 

END MODULE module_functions

