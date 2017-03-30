MODULE module_read_data_header 
 
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil
  use type
  use module_convert
  use module_functions
  use module_functions_parse

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

!############### HEADER v3 ##############################################
 ! Read data header
  SUBROUTINE Read_Data_Headerv3( Unit, Coefs, DataType, sizeOfRecords, &
       AVHRR_No, nScans )

#ifdef USE_GZFILE
    INTEGER(c_int), INTENT(IN) :: Unit
#else
    INTEGER, INTENT(IN) :: Unit
#endif
    TYPE(AVHRR_Instrument_Coefs), INTENT(OUT) :: Coefs
    INTEGER, INTENT(OUT) :: DataType
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER, INTENT(OUT) :: AVHRR_No
    INTEGER, INTENT(OUT) :: nScans

    INTEGER :: version
    INTEGER :: intVal

    INTEGER :: STAT

    INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:) :: Header

    INTEGER :: earthLocationBitField

    INTEGER :: I
    INTEGER :: nHeaders

    ALLOCATE(Header(SizeOfRecords),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating Header','Read_Data_Headerv3',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    ! Read header
#ifdef USE_GZFILE
    STAT = gzread(Unit,Header)
    IF(STAT.le.0)THEN
       CALL Check_IOS(1,'Reading Header','Read_Data_Headerv3',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#else
    READ(Unit,IOSTAT=STAT)Header
    CALL Check_IOS(STAT,'Reading Header','Read_Data_Headerv3',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif
    version = Header2ushort( 4, Header )
    if( 3 .gt. version )then
       call out_ERROR('Version Number incorrect',&
            'Read_Data_HeaderV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    endif
    intVal = Header2ushort( 76, Header )
    SELECT CASE(intVal)
    CASE(1)
       print *,'Data is in LAC format'
       DataType = AVHRR_LAC
    CASE(2)
       print *,'Data is in GAC format'
       DataType = AVHRR_GAC
    CASE(13)
       print *,'Data is in FRAC format'
       DataType = AVHRR_LAC
    CASE DEFAULT
       print *,intval
       call out_ERROR('Data not GAC/LAC/FRAC',&
            'Read_Data_HeaderV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END SELECT

    ! older AVHRRs are based on physical records, 2 logical records per physical record
    ! so can have an extra (fill) scan at the end of the file if nScans is odd
    ! need this to keep in step with CLAVRx (which I have also modified), used for cloud mask
    nScans = Header2ushort( 128, Header )

    start_time = AVHRR_Date(header2ushort(84,header), header2ushort(86,header), header2int(88,header))
    end_time   = AVHRR_Date(header2ushort(96,header), header2ushort(98,header), header2int(100,header))

!    ! This is not very informative, it doesn't mean that the attitude is bad
!    earthLocationBitField = header2ushort( 338, Header )
!    IF(.NOT.(BTEST(earthLocationBitField,2) &
!         .OR.BTEST(earthLocationBitField,0)))THEN
!       CALL out_WARNING('No attitude correction: geolocation and land masking may be inaccurate',&
!            'Read_Data_HeaderV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
!    END IF
    
    ! Values for PRT conversion */
    Coefs%prtTempCoefs1(1) = header2short( 200, Header )/1e2
    Coefs%prtTempCoefs2(1) = header2short( 202, Header )/1e5
    Coefs%prtTempCoefs3(1) = header2short( 204, Header )/1e8
    Coefs%prtTempCoefs4(1) = header2short( 206, Header )/1e11
    Coefs%prtTempCoefs5(1) = header2short( 208, Header )/1e14
    Coefs%prtTempCoefs6(1) = header2short( 210, Header )/1e17
    Coefs%prtTempCoefs1(2) = header2short( 212, Header )/1e2
    Coefs%prtTempCoefs2(2) = header2short( 214, Header )/1e5
    Coefs%prtTempCoefs3(2) = header2short( 216, Header )/1e8
    Coefs%prtTempCoefs4(2) = header2short( 218, Header )/1e11
    Coefs%prtTempCoefs5(2) = header2short( 220, Header )/1e14
    Coefs%prtTempCoefs6(2) = header2short( 222, Header )/1e17
    Coefs%prtTempCoefs1(3) = header2short( 224, Header )/1e2
    Coefs%prtTempCoefs2(3) = header2short( 226, Header )/1e5
    Coefs%prtTempCoefs3(3) = header2short( 228, Header )/1e8
    Coefs%prtTempCoefs4(3) = header2short( 230, Header )/1e11
    Coefs%prtTempCoefs5(3) = header2short( 232, Header )/1e14
    Coefs%prtTempCoefs6(3) = header2short( 234, Header )/1e17
    Coefs%prtTempCoefs1(4) = header2short( 236, Header )/1e2
    Coefs%prtTempCoefs2(4) = header2short( 238, Header )/1e5
    Coefs%prtTempCoefs3(4) = header2short( 240, Header )/1e8
    Coefs%prtTempCoefs4(4) = header2short( 242, Header )/1e11
    Coefs%prtTempCoefs5(4) = header2short( 244, Header )/1e14
    Coefs%prtTempCoefs6(4) = header2short( 246, Header )/1e17

    ! Get radiance -> temperature coefficients 
    Coefs%nuC(1)  = header2int(280,Header)/1e2
    Coefs%Aval(1) = header2int(284,Header)/1e5
    Coefs%Bval(1) = header2int(288,Header)/1e6
    Coefs%nuC(2)  = header2int(292,Header)/1e3
    Coefs%Aval(2) = header2int(296,Header)/1e5
    Coefs%Bval(2) = header2int(300,Header)/1e6
    Coefs%nuC(3)  = header2int(304,Header)/1e3
    Coefs%Aval(3) = header2int(308,Header)/1e5
    Coefs%Bval(3) = header2int(312,Header)/1e6

!    print *,'Calibration coefficients (SRF) : '
!    print *,Coefs%nuC(1:3)
!    print *,Coefs%Aval(1:3)
!    print *,Coefs%Bval(1:3)
!    STOP 1

    Coefs%patchCoef(1) = header2int(424,Header)/1e6
    Coefs%patchCoef(2) = header2int(428,Header)/1e6
    Coefs%patchCoef(3) = header2int(432,Header)/1e7
    Coefs%patchCoef(4) = header2int(436,Header)/1e8
    Coefs%patchCoef(5) = header2int(440,Header)/1e9
    Coefs%patchCoef(6) = header2int(444,Header)/1e10

    Coefs%patchCoefExt(1) = header2int(448,Header)/1e6
    Coefs%patchCoefExt(2) = header2int(452,Header)/1e6
    Coefs%patchCoefExt(3) = header2int(456,Header)/1e7
    Coefs%patchCoefExt(4) = header2int(460,Header)/1e8
    Coefs%patchCoefExt(5) = header2int(464,Header)/1e9
    Coefs%patchCoefExt(6) = header2int(468,Header)/1e10
    
    ! Temperature coefficients 
    ! Radiator */
    Coefs%temperatureCoefs1(1) = header2int(496,Header)/1e6
    Coefs%temperatureCoefs1(2) = header2int(500,Header)/1e6
    Coefs%temperatureCoefs1(3) = header2int(504,Header)/1e7
    Coefs%temperatureCoefs1(4) = header2int(508,Header)/1e8
    Coefs%temperatureCoefs1(5) = header2int(512,Header)/1e9
    Coefs%temperatureCoefs1(6) = header2int(516,Header)/1e10    
    ! Electronics 
    Coefs%temperatureCoefs2(1) = header2int(688,Header)/1e6
    Coefs%temperatureCoefs2(2) = header2int(692,Header)/1e6
    Coefs%temperatureCoefs2(3) = header2int(696,Header)/1e7
    Coefs%temperatureCoefs2(4) = header2int(700,Header)/1e8
    Coefs%temperatureCoefs2(5) = header2int(704,Header)/1e9
    Coefs%temperatureCoefs2(6) = header2int(708,Header)/1e10
    ! Cooler 
    Coefs%temperatureCoefs3(1) = header2int(712,Header)/1e6
    Coefs%temperatureCoefs3(2) = header2int(716,Header)/1e6
    Coefs%temperatureCoefs3(3) = header2int(720,Header)/1e7
    Coefs%temperatureCoefs3(4) = header2int(724,Header)/1e8
    Coefs%temperatureCoefs3(5) = header2int(728,Header)/1e9
    Coefs%temperatureCoefs3(6) = header2int(732,Header)/1e10
    ! Baseplate 
    Coefs%temperatureCoefs4(1) = header2int(736,Header)/1e6
    Coefs%temperatureCoefs4(2) = header2int(740,Header)/1e6
    Coefs%temperatureCoefs4(3) = header2int(744,Header)/1e7
    Coefs%temperatureCoefs4(4) = header2int(748,Header)/1e8
    Coefs%temperatureCoefs4(5) = header2int(752,Header)/1e9
    Coefs%temperatureCoefs4(6) = header2int(756,Header)/1e10
    ! Motor 
    Coefs%temperatureCoefs5(1) = header2int(760,Header)/1e6
    Coefs%temperatureCoefs5(2) = header2int(764,Header)/1e6
    Coefs%temperatureCoefs5(3) = header2int(768,Header)/1e7
    Coefs%temperatureCoefs5(4) = header2int(772,Header)/1e8
    Coefs%temperatureCoefs5(5) = header2int(776,Header)/1e9
    Coefs%temperatureCoefs5(6) = header2int(780,Header)/1e10
    ! A/D 
    Coefs%temperatureCoefs6(1) = header2int(784,Header)/1e6
    Coefs%temperatureCoefs6(2) = header2int(788,Header)/1e6
    Coefs%temperatureCoefs6(3) = header2int(792,Header)/1e7
    Coefs%temperatureCoefs6(4) = header2int(796,Header)/1e8
    Coefs%temperatureCoefs6(5) = header2int(800,Header)/1e9
    Coefs%temperatureCoefs6(6) = header2int(804,Header)/1e10
    ! Patch 
    Coefs%temperatureCoefs7(1) = header2int(424,Header)/1e6
    Coefs%temperatureCoefs7(2) = header2int(428,Header)/1e6
    Coefs%temperatureCoefs7(3) = header2int(432,Header)/1e7
    Coefs%temperatureCoefs7(4) = header2int(436,Header)/1e8
    Coefs%temperatureCoefs7(5) = header2int(440,Header)/1e9
    Coefs%temperatureCoefs7(6) = header2int(444,Header)/1e10
  
    ! Motor current coefficients
    Coefs%motorCurrentCoefs(1) = header2int(640,Header)/1e6
    Coefs%motorCurrentCoefs(2) = header2int(644,Header)/1e6
    Coefs%motorCurrentCoefs(3) = header2int(648,Header)/1e7
    Coefs%motorCurrentCoefs(4) = header2int(652,Header)/1e8
    Coefs%motorCurrentCoefs(5) = header2int(656,Header)/1e9
    Coefs%motorCurrentCoefs(6) = header2int(660,Header)/1e10

    intVal = Header2ushort( 72, Header )
    SELECT CASE(intVal)
    CASE(4)
       AVHRR_No = 15
       write(*,*)'Dataset is from AVHRR-15'
!       ! Constant values for PRT conversion */
!       Coefs%prtTempCoefs1(1) = 276.60157
!       Coefs%prtTempCoefs2(1) = 0.051045
!       Coefs%prtTempCoefs3(1) = 1.36328E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.62531
!       Coefs%prtTempCoefs2(2) = 0.050909
!       Coefs%prtTempCoefs3(2) = 1.47266E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.67413
!       Coefs%prtTempCoefs2(3) = 0.050907
!       Coefs%prtTempCoefs3(3) = 1.47656E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.59258
!       Coefs%prtTempCoefs2(4) = 0.050966
!       Coefs%prtTempCoefs3(4) = 1.47656E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -4.50
       Coefs%Nspace(3) = -3.61
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 4.76
       Coefs%nonLinearCoefs4(2) = -0.0932
       Coefs%nonLinearCoefs4(3) = 0.0004524
       Coefs%nonLinearCoefs5(1) = 3.83
       Coefs%nonLinearCoefs5(2) = -0.0659
       Coefs%nonLinearCoefs5(3) = 0.0002811
    CASE(2)
       AVHRR_No = 16
       write(*,*)'Dataset is from AVHRR-16'
!       Coefs%prtTempCoefs1(1) = 276.355
!       Coefs%prtTempCoefs2(1) = 5.562E-02
!       Coefs%prtTempCoefs3(1) = -1.590E-05
!       Coefs%prtTempCoefs4(1) = 2.486E-08
!       Coefs%prtTempCoefs5(1) = -1.199E-11
!       Coefs%prtTempCoefs1(2) = 276.142
!       Coefs%prtTempCoefs2(2) = 5.605E-02
!       Coefs%prtTempCoefs3(2) = -1.707E-05
!       Coefs%prtTempCoefs4(2) = 2.595E-08
!       Coefs%prtTempCoefs5(2) = -1.199E-11
!       Coefs%prtTempCoefs1(3) = 275.996
!       Coefs%prtTempCoefs2(3) = 5.486E-02
!       Coefs%prtTempCoefs3(3) = -1.223E-05
!       Coefs%prtTempCoefs4(3) = 1.862E-08
!       Coefs%prtTempCoefs5(3) = -0.853E-11
!       Coefs%prtTempCoefs1(4) = 276.132
!       Coefs%prtTempCoefs2(4) = 5.494E-02
!       Coefs%prtTempCoefs3(4) = -1.344E-05
!       Coefs%prtTempCoefs4(4) = 2.112E-08
!       Coefs%prtTempCoefs5(4) = -1.001E-11
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -2.467
       Coefs%Nspace(3) = -2.009
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 2.96
       Coefs%nonLinearCoefs4(2) = -0.05411
       Coefs%nonLinearCoefs4(3) = 0.00024532
       Coefs%nonLinearCoefs5(1) = 2.25
       Coefs%nonLinearCoefs5(2) = -0.03665
       Coefs%nonLinearCoefs5(3) = 0.00014854
    CASE(6)
       AVHRR_No = 17
       write(*,*)'Dataset is from AVHRR-17'
!       Coefs%prtTempCoefs1(1) = 276.628
!       Coefs%prtTempCoefs2(1) = 0.05098
!       Coefs%prtTempCoefs3(1) = 1.371E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.538
!       Coefs%prtTempCoefs2(2) = 0.05098
!       Coefs%prtTempCoefs3(2) = 1.371E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.761
!       Coefs%prtTempCoefs2(3) = 0.05097
!       Coefs%prtTempCoefs3(3) = 1.369E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.660
!       Coefs%prtTempCoefs2(4) = 0.05100
!       Coefs%prtTempCoefs3(4) = 1.348E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -8.55
       Coefs%Nspace(3) = -3.97
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 8.22
       Coefs%nonLinearCoefs4(2) = -0.15795
       Coefs%nonLinearCoefs4(3) = 0.00075579
       Coefs%nonLinearCoefs5(1) = 4.31
       Coefs%nonLinearCoefs5(2) = -0.07318
       Coefs%nonLinearCoefs5(3) = 0.00030976
    CASE(7)
       AVHRR_No = 18
       write(*,*)'Dataset is from AVHRR-18'
!       Coefs%prtTempCoefs1(1) = 276.601
!       Coefs%prtTempCoefs2(1) = 0.05090
!       Coefs%prtTempCoefs3(1) = 1.657E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.683
!       Coefs%prtTempCoefs2(2) = 0.05101
!       Coefs%prtTempCoefs3(2) = 1.482E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.565
!       Coefs%prtTempCoefs2(3) = 0.05117
!       Coefs%prtTempCoefs3(3) = 1.313E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.615
!       Coefs%prtTempCoefs2(4) = 0.05103
!       Coefs%prtTempCoefs3(4) = 1.484E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -5.33
       Coefs%Nspace(3) = -2.22
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 5.82
       Coefs%nonLinearCoefs4(2) = -0.11069
       Coefs%nonLinearCoefs4(3) = 0.00052337
       Coefs%nonLinearCoefs5(1) = 2.67
       Coefs%nonLinearCoefs5(2) = -0.04360
       Coefs%nonLinearCoefs5(3) = 0.00017715
    CASE(8)
       AVHRR_No = 19
       write(*,*)'Dataset is from AVHRR-19'
!       Coefs%prtTempCoefs1(1) = 276.6067
!       Coefs%prtTempCoefs2(1) = 0.051111
!       Coefs%prtTempCoefs3(1) = 1.405783E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.6119
!       Coefs%prtTempCoefs2(2) = 0.051090
!       Coefs%prtTempCoefs3(2) = 1.496037E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.6311
!       Coefs%prtTempCoefs2(3) = 0.051033
!       Coefs%prtTempCoefs3(3) = 1.496990E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.6268
!       Coefs%prtTempCoefs2(4) = 0.051058
!       Coefs%prtTempCoefs3(4) = 1.493110E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -5.49
       Coefs%Nspace(3) = -3.39
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 5.70
       Coefs%nonLinearCoefs4(2) = -0.11187
       Coefs%nonLinearCoefs4(3) = 0.00054668
       Coefs%nonLinearCoefs5(1) = 3.58
       Coefs%nonLinearCoefs5(2) = -0.05991
       Coefs%nonLinearCoefs5(3) = 0.00024985
    CASE(12)
       AVHRR_No = -1 ! Metop-A
       write(*,*)'Dataset is from Metop-A'
!       Coefs%prtTempCoefs1(1) = 276.6194
!       Coefs%prtTempCoefs2(1) = 0.050919
!       Coefs%prtTempCoefs3(1) = 1.470892E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.6511
!       Coefs%prtTempCoefs2(2) = 0.050892
!       Coefs%prtTempCoefs3(2) = 1.489000E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.6597
!       Coefs%prtTempCoefs2(3) = 0.05845
!       Coefs%prtTempCoefs3(3) = 1.520646E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.3685
!       Coefs%prtTempCoefs2(4) = 0.050992
!       Coefs%prtTempCoefs3(4) = 1.482390E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -4.98
       Coefs%Nspace(3) = -3.40
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 5.44
       Coefs%nonLinearCoefs4(2) = -0.10152
       Coefs%nonLinearCoefs4(3) = 0.00046964
       Coefs%nonLinearCoefs5(1) = 3.84
       Coefs%nonLinearCoefs5(2) = -0.06249
       Coefs%nonLinearCoefs5(3) = 0.00025239
!*************************
! modif du 8 avril 2016
! on ajoute CASE(-2) pour Metop-B 
    CASE(11)
      AVHRR_No = -2 ! Metop-B
      write(*,*)'Dataset is from Metop-B'
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -4.75
       Coefs%Nspace(3) = -4.39
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 4.85
       Coefs%nonLinearCoefs4(2) = -0.0096771
       Coefs%nonLinearCoefs4(3) = 0.00048091
       Coefs%nonLinearCoefs5(1) = 4.36
       Coefs%nonLinearCoefs5(2) = -0.07663650
       Coefs%nonLinearCoefs5(3) = 0.00033524
!*************************
    CASE DEFAULT
       call out_ERROR('Invalid Satellite ID',&
            'Read_Data_HeaderV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')       
    END SELECT

    ! parameters for reflectance to radiance conversion for the visible channels
    ! see NOAA KLM User's Guide Section 7.1.1.1
    ! (http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c7/sec7-1.htm)
    ! Eqs 7.1.1.1-2 to 7.1.1.1-6
    ! I_lambda = reflectance / 100 * F / pi / w
    ! I_lambda is the radiance W m-2 sr-1 um-1
    ! reflectance is in percent
    ! F is the band-integrated extraterrestrial solar irradiance at normal incidence
    !   at the top of the atmosphere at mean Earth-Sun distance W m-2
    ! w is the band equivalent width um

    Coefs%F(1)  = header2int(256,Header)/1e1
    Coefs%w(1) = header2int(260,Header)/1e3
    Coefs%F(2)  = header2int(264,Header)/1e1
    Coefs%w(2) = header2int(268,Header)/1e3
    Coefs%F(3)  = header2int(272,Header)/1e1
    Coefs%w(3) = header2int(276,Header)/1e3

    ! If extra headers are there - skip them
    nHeaders = header2ushort(14,Header)
    DO I=1,nHeaders-1
       ! Read extra header
#ifdef USE_GZFILE
       STAT = gzread(Unit,Header)
       IF(STAT.LE.0)THEN
          CALL Check_IOS(1,'Reading Extra Headers','Read_Data_Headerv3',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       END IF
#else
       READ(Unit,IOSTAT=STAT)Header
       CALL Check_IOS(STAT,'Reading Header','Read_Data_Headerv3',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif
    END DO

    ! Free header
    IF(ALLOCATED(Header))DEALLOCATE(Header)
    
  END SUBROUTINE Read_Data_Headerv3
!############################################################

!############### HEADER v2 ##############################################
  SUBROUTINE Read_Data_Headerv2( Unit, Coefs, DataType, sizeOfRecords, &
       AVHRR_No, nScans )

#ifdef USE_GZFILE
    INTEGER(c_int), INTENT(IN) :: Unit
#else
    INTEGER, INTENT(IN) :: Unit
#endif
    TYPE(AVHRR_Instrument_Coefs), INTENT(OUT) :: Coefs
    INTEGER, INTENT(OUT) :: DataType
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER, INTENT(OUT) :: AVHRR_No
    INTEGER, INTENT(OUT) :: nScans

    INTEGER(GbcsInt2) :: version
    INTEGER(GbcsInt2) :: intVal

    INTEGER :: STAT

    INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:) :: Header

    ALLOCATE(Header(SizeOfRecords),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating Header','Read_Data_Headerv2',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    ! Read header
#ifdef USE_GZFILE
    STAT = gzread(Unit,Header)
    IF(STAT.LE.0)THEN
       CALL Check_IOS(1,'Reading Header','Read_Data_Headerv2',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#else
    READ(Unit,IOSTAT=STAT)Header
    CALL Check_IOS(STAT,'Reading Header','Read_Data_Headerv2',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

    version = Header2Short( 4, Header )
    ! Some NOAA-15 data has version 1 which is (according to the KLM guide)
    ! impossible - so assume it's version 2
    if( 2 .ne. version .and. 1 .ne. version )then
       call out_ERROR('Version Number incorrect',&
            'Read_Data_HeaderV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    endif
    intVal = Header2Short( 76, Header )
    SELECT CASE(INT(intVal))
    CASE(1)
       DataType = AVHRR_LAC
    CASE(2)
       DataType = AVHRR_GAC
    CASE DEFAULT
       call out_ERROR('Data not GAC/LAC/HRPT',&
            'Read_Data_HeaderV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END SELECT
    
    ! older AVHRRs are based on physical records, 2 logical records per physical record
    ! so can have an extra (fill) scan at the end of the file if nScans is odd
    ! need this to keep in step with CLAVRx (which I have also modified), used for cloud mask
    nScans = Header2ushort( 128, Header )

    start_time = AVHRR_Date(header2ushort(84,header), header2ushort(86,header), header2int(88,header))
    end_time   = AVHRR_Date(header2ushort(96,header), header2ushort(98,header), header2int(100,header))

    ! Values for PRT conversion */
    Coefs%prtTempCoefs1(1) = header2short( 200, Header )/1e2
    Coefs%prtTempCoefs2(1) = header2short( 202, Header )/1e5
    Coefs%prtTempCoefs3(1) = header2short( 204, Header )/1e8
    Coefs%prtTempCoefs4(1) = header2short( 206, Header )/1e11
    Coefs%prtTempCoefs5(1) = header2short( 208, Header )/1e14
    Coefs%prtTempCoefs6(1) = header2short( 210, Header )/1e17
    Coefs%prtTempCoefs1(2) = header2short( 212, Header )/1e2
    Coefs%prtTempCoefs2(2) = header2short( 214, Header )/1e5
    Coefs%prtTempCoefs3(2) = header2short( 216, Header )/1e8
    Coefs%prtTempCoefs4(2) = header2short( 218, Header )/1e11
    Coefs%prtTempCoefs5(2) = header2short( 220, Header )/1e14
    Coefs%prtTempCoefs6(2) = header2short( 222, Header )/1e17
    Coefs%prtTempCoefs1(3) = header2short( 224, Header )/1e2
    Coefs%prtTempCoefs2(3) = header2short( 226, Header )/1e5
    Coefs%prtTempCoefs3(3) = header2short( 228, Header )/1e8
    Coefs%prtTempCoefs4(3) = header2short( 230, Header )/1e11
    Coefs%prtTempCoefs5(3) = header2short( 232, Header )/1e14
    Coefs%prtTempCoefs6(3) = header2short( 234, Header )/1e17
    Coefs%prtTempCoefs1(4) = header2short( 236, Header )/1e2
    Coefs%prtTempCoefs2(4) = header2short( 238, Header )/1e5
    Coefs%prtTempCoefs3(4) = header2short( 240, Header )/1e8
    Coefs%prtTempCoefs4(4) = header2short( 242, Header )/1e11
    Coefs%prtTempCoefs5(4) = header2short( 244, Header )/1e14
    Coefs%prtTempCoefs6(4) = header2short( 246, Header )/1e17

    ! Get Radiance -> temperature coefficients
    Coefs%nuC(1)  = header2int(280,Header)/1e2
    Coefs%Aval(1) = header2int(284,Header)/1e5
    Coefs%Bval(1) = header2int(288,Header)/1e6
    Coefs%nuC(2)  = header2int(292,Header)/1e3
    Coefs%Aval(2) = header2int(296,Header)/1e5
    Coefs%Bval(2) = header2int(300,Header)/1e6
    Coefs%nuC(3)  = header2int(304,Header)/1e3
    Coefs%Aval(3) = header2int(308,Header)/1e5
    Coefs%Bval(3) = header2int(312,Header)/1e6

    Coefs%patchCoef(1) = header2short(424,Header)/1e2
    Coefs%patchCoef(2) = header2short(426,Header)/1e2
    Coefs%patchCoef(3) = header2short(428,Header)/1e2
    Coefs%patchCoef(4) = header2short(430,Header)/1e2
    Coefs%patchCoef(5) = header2short(442,Header)/1e2
    Coefs%patchCoef(6) = 0.

    Coefs%patchCoefExt(1) = header2short(436,Header)/1e2
    Coefs%patchCoefExt(2) = header2short(438,Header)/1e2
    Coefs%patchCoefExt(3) = header2short(440,Header)/1e2
    Coefs%patchCoefExt(4) = header2short(442,Header)/1e2
    Coefs%patchCoefExt(5) = header2short(444,Header)/1e2
    Coefs%patchCoefExt(6) = 0.

    ! Temperature coefficients 
    ! Radiator 
    Coefs%temperatureCoefs1(1) = header2short(460,Header)/1e2
    Coefs%temperatureCoefs1(2) = header2short(462,Header)/1e2
    Coefs%temperatureCoefs1(3) = header2short(464,Header)/1e2
    Coefs%temperatureCoefs1(4) = header2short(466,Header)/1e2
    Coefs%temperatureCoefs1(5) = header2short(468,Header)/1e2
    Coefs%temperatureCoefs1(6) = 0.
    ! Electronics 
    Coefs%temperatureCoefs2(1) =  header2short(556,Header)/1e2
    Coefs%temperatureCoefs2(2) =  header2short(558,Header)/1e2
    Coefs%temperatureCoefs2(3) =  header2short(560,Header)/1e2
    Coefs%temperatureCoefs2(4) =  header2short(562,Header)/1e2
    Coefs%temperatureCoefs2(5) =  header2short(564,Header)/1e2
    Coefs%temperatureCoefs2(6) = 0.
    ! Cooler 
    Coefs%temperatureCoefs3(1) =  header2short(568,Header)/1e2
    Coefs%temperatureCoefs3(2) =  header2short(570,Header)/1e2
    Coefs%temperatureCoefs3(3) =  header2short(572,Header)/1e2
    Coefs%temperatureCoefs3(4) =  header2short(574,Header)/1e2
    Coefs%temperatureCoefs3(5) =  header2short(576,Header)/1e2
    Coefs%temperatureCoefs3(6) = 0.
    ! Baseplate 
    Coefs%temperatureCoefs4(1) =  header2short(580,Header)/1e2
    Coefs%temperatureCoefs4(2) =  header2short(582,Header)/1e2
    Coefs%temperatureCoefs4(3) =  header2short(584,Header)/1e2
    Coefs%temperatureCoefs4(4) =  header2short(586,Header)/1e2
    Coefs%temperatureCoefs4(5) =  header2short(588,Header)/1e2
    Coefs%temperatureCoefs4(6) = 0.
    ! Motor 
    Coefs%temperatureCoefs5(1) =  header2short(592,Header)/1e2
    Coefs%temperatureCoefs5(2) =  header2short(594,Header)/1e2
    Coefs%temperatureCoefs5(3) =  header2short(596,Header)/1e2
    Coefs%temperatureCoefs5(4) =  header2short(598,Header)/1e2
    Coefs%temperatureCoefs5(5) =  header2short(600,Header)/1e2
    Coefs%temperatureCoefs5(6) = 0.
    ! A/D 
    Coefs%temperatureCoefs6(1) =  header2short(604,Header)/1e2
    Coefs%temperatureCoefs6(2) =  header2short(606,Header)/1e2
    Coefs%temperatureCoefs6(3) =  header2short(608,Header)/1e2
    Coefs%temperatureCoefs6(4) =  header2short(610,Header)/1e2
    Coefs%temperatureCoefs6(5) =  header2short(612,Header)/1e2
    Coefs%temperatureCoefs6(6) = 0.
    ! Patch 
    Coefs%temperatureCoefs7(1) =  header2short(424,Header)/1e2
    Coefs%temperatureCoefs7(2) =  header2short(426,Header)/1e2
    Coefs%temperatureCoefs7(3) =  header2short(428,Header)/1e2
    Coefs%temperatureCoefs7(4) =  header2short(430,Header)/1e2
    Coefs%temperatureCoefs7(5) =  header2short(432,Header)/1e2
    Coefs%temperatureCoefs7(6) = 0.

    ! Motor current coefficients
    Coefs%motorCurrentCoefs(1) = header2short(532,Header)/1e2
    Coefs%motorCurrentCoefs(2) = header2short(534,Header)/1e2
    Coefs%motorCurrentCoefs(3) = header2short(536,Header)/1e2
    Coefs%motorCurrentCoefs(4) = header2short(538,Header)/1e2
    Coefs%motorCurrentCoefs(5) = header2short(540,Header)/1e2
    Coefs%motorCurrentCoefs(6) = 0.

    ! Get parameters not stored in file
    intVal = header2short(72,Header)
    SELECT CASE(INT(intVal))
    CASE(4)
       AVHRR_No = 15
       write(*,*)'Dataset is from AVHRR-15'
       ! Constant values for PRT conversion 
!       Coefs%prtTempCoefs1(1) = 276.60157
!       Coefs%prtTempCoefs2(1) = 0.051045
!       Coefs%prtTempCoefs3(1) = 1.36328E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.62531
!       Coefs%prtTempCoefs2(2) = 0.050909
!       Coefs%prtTempCoefs3(2) = 1.47266E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.67413
!       Coefs%prtTempCoefs2(3) = 0.050907
!       Coefs%prtTempCoefs3(3) = 1.47656E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.59258
!       Coefs%prtTempCoefs2(4) = 0.050966
!       Coefs%prtTempCoefs3(4) = 1.47656E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -4.50
       Coefs%Nspace(3) = -3.61
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 4.76
       Coefs%nonLinearCoefs4(2) = -0.0932
       Coefs%nonLinearCoefs4(3) = 0.0004524
       Coefs%nonLinearCoefs5(1) = 3.83
       Coefs%nonLinearCoefs5(2) = -0.0659
       Coefs%nonLinearCoefs5(3) = 0.0002811
    CASE(2)
       AVHRR_No = 16
       write(*,*)'Dataset is from AVHRR-16'
!       Coefs%prtTempCoefs1(1) = 276.355
!       Coefs%prtTempCoefs2(1) = 5.562E-02
!       Coefs%prtTempCoefs3(1) = -1.590E-05
!       Coefs%prtTempCoefs4(1) = 2.486E-08
!       Coefs%prtTempCoefs5(1) = -1.199E-11
!       Coefs%prtTempCoefs1(2) = 276.142
!       Coefs%prtTempCoefs2(2) = 5.605E-02
!       Coefs%prtTempCoefs3(2) = -1.707E-05
!       Coefs%prtTempCoefs4(2) = 2.595E-08
!       Coefs%prtTempCoefs5(2) = -1.199E-11
!       Coefs%prtTempCoefs1(3) = 275.996
!       Coefs%prtTempCoefs2(3) = 5.486E-02
!       Coefs%prtTempCoefs3(3) = -1.223E-05
!       Coefs%prtTempCoefs4(3) = 1.862E-08
!       Coefs%prtTempCoefs5(3) = -0.853E-11
!       Coefs%prtTempCoefs1(4) = 276.132
!       Coefs%prtTempCoefs2(4) = 5.494E-02
!       Coefs%prtTempCoefs3(4) = -1.344E-05
!       Coefs%prtTempCoefs4(4) = 2.112E-08
!       Coefs%prtTempCoefs5(4) = -1.001E-11
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -2.467
       Coefs%Nspace(3) = -2.009
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 2.96
       Coefs%nonLinearCoefs4(2) = -0.05411
       Coefs%nonLinearCoefs4(3) = 0.00024532
       Coefs%nonLinearCoefs5(1) = 2.25
       Coefs%nonLinearCoefs5(2) = -0.03665
       Coefs%nonLinearCoefs5(3) = 0.00014854
    CASE(6)
       AVHRR_No = 17
       write(*,*)'Dataset is from AVHRR-17'
!       Coefs%prtTempCoefs1(1) = 276.628
!       Coefs%prtTempCoefs2(1) = 0.05098
!       Coefs%prtTempCoefs3(1) = 1.371E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.538
!       Coefs%prtTempCoefs2(2) = 0.05098
!       Coefs%prtTempCoefs3(2) = 1.371E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.761
!       Coefs%prtTempCoefs2(3) = 0.05097
!       Coefs%prtTempCoefs3(3) = 1.369E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.660
!       Coefs%prtTempCoefs2(4) = 0.05100
!       Coefs%prtTempCoefs3(4) = 1.348E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -8.55
       Coefs%Nspace(3) = -3.97
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 8.22
       Coefs%nonLinearCoefs4(2) = -0.15795
       Coefs%nonLinearCoefs4(3) = 0.00075579
       Coefs%nonLinearCoefs5(1) = 4.31
       Coefs%nonLinearCoefs5(2) = -0.07318
       Coefs%nonLinearCoefs5(3) = 0.00030976
    CASE(7)
       AVHRR_No = 18
       write(*,*)'Dataset is from AVHRR-18'
!       Coefs%prtTempCoefs1(1) = 276.601
!       Coefs%prtTempCoefs2(1) = 0.05090
!       Coefs%prtTempCoefs3(1) = 1.657E-06
!       Coefs%prtTempCoefs4(1) = 0.
!       Coefs%prtTempCoefs5(1) = 0.
!       Coefs%prtTempCoefs1(2) = 276.683
!       Coefs%prtTempCoefs2(2) = 0.05101
!       Coefs%prtTempCoefs3(2) = 1.482E-06
!       Coefs%prtTempCoefs4(2) = 0.
!       Coefs%prtTempCoefs5(2) = 0.
!       Coefs%prtTempCoefs1(3) = 276.565
!       Coefs%prtTempCoefs2(3) = 0.05117
!       Coefs%prtTempCoefs3(3) = 1.313E-06
!       Coefs%prtTempCoefs4(3) = 0.
!       Coefs%prtTempCoefs5(3) = 0.
!       Coefs%prtTempCoefs1(4) = 276.615
!       Coefs%prtTempCoefs2(4) = 0.05103
!       Coefs%prtTempCoefs3(4) = 1.484E-06
!       Coefs%prtTempCoefs4(4) = 0.
!       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -5.53
       Coefs%Nspace(3) = -2.22
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 5.82
       Coefs%nonLinearCoefs4(2) = -0.11069
       Coefs%nonLinearCoefs4(3) = 0.00052337
       Coefs%nonLinearCoefs5(1) = 2.67
       Coefs%nonLinearCoefs5(2) = -0.04360
       Coefs%nonLinearCoefs5(3) = 0.00017715
    CASE DEFAULT
       call out_ERROR('Invalid Satellite ID',&
            'Read_Data_HeaderV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')       
    END SELECT

    ! parameters for reflectance to radiance conversion for the visible channels
    ! see NOAA KLM User's Guide Section 7.1
    ! (http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c7/sec7-1.htm)
    ! Eqs 7.1.1.1-2 to 7.1.1.1-6
    ! I_lambda = reflectance / 100 * F / pi / w
    ! I_lambda is the radiance W m-2 sr-1 um-1
    ! reflectance is in percent
    ! F is the band-integrated extraterrestrial solar irradiance at normal incidence
    !   at the top of the atmosphere at mean Earth-Sun distance W m-2
    ! w is the band equivalent width um

    Coefs%F(1)  = header2int(256,Header)/1e1
    Coefs%w(1) = header2int(260,Header)/1e3
    Coefs%F(2)  = header2int(264,Header)/1e1
    Coefs%w(2) = header2int(268,Header)/1e3
    Coefs%F(3)  = header2int(272,Header)/1e1
    Coefs%w(3) = header2int(276,Header)/1e3

    ! Free memory
    IF(ALLOCATED(Header))DEALLOCATE(Header)

  END SUBROUTINE Read_Data_Headerv2
!############### HEADER v3 ##############################################

!############### HEADER v1 ##############################################
  SUBROUTINE Read_Data_Headerv1( Unit, Coefs, DataType, sizeOfRecords, &
       AVHRR_No, have_12micron, nScans )

#ifdef USE_GZFILE
    INTEGER(c_int), INTENT(IN) :: Unit
#else
    INTEGER, INTENT(IN) :: Unit
#endif
    TYPE(AVHRR_Instrument_Coefs), INTENT(OUT) :: Coefs
    INTEGER, INTENT(OUT) :: DataType
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER, INTENT(OUT) :: AVHRR_No
    LOGICAL, INTENT(OUT) :: have_12micron
    INTEGER, INTENT(OUT) :: nScans

    INTEGER(GbcsInt2) :: version
    INTEGER(GbcsInt2) :: intVal
    INTEGER :: year
    INTEGER :: dayno
    REAL :: hours
    INTEGER :: UTC_msecs

    INTEGER :: STAT

    INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:) :: Header

    ALLOCATE(Header(SizeOfRecords),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating Header','Read_Data_Headerv2',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    ! Read header
#ifdef USE_GZFILE
    STAT = gzread(Unit,Header)
    IF(STAT.LE.0)THEN
       CALL Check_IOS(1,'Reading Header','Read_Data_HeaderV1',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#else
    READ(Unit,IOSTAT=STAT)Header
    CALL Check_IOS(STAT,'Reading Header','Read_Data_Headerv1',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif
    intVal = Parse_Data_Type( Header(2) )
    SELECT CASE(INT(intVal))
    CASE(1,3)
       DataType = AVHRR_LAC
    CASE(2)
       DataType = AVHRR_GAC
    CASE DEFAULT
       call out_ERROR('Data not GAC/LAC/HRPT',&
            'Read_Data_HeaderV1','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END SELECT

    ! older AVHRRs are based on physical records, 2 logical records per physical record
    ! so can have an extra (fill) scan at the end of the file if nScans is odd
    ! need this to keep in step with CLAVRx (which I have also modified), used for cloud mask
    nScans = Header2ushort( 8, Header )

    CALL Parse_Date_Time( header, 3, year, dayno, hours, UTC_msecs )
    start_time = AVHRR_Date(year, dayno, UTC_msecs)
    CALL Parse_Date_Time( header, 11, year, dayno, hours, UTC_msecs )
    end_time = AVHRR_Date(year, dayno, UTC_msecs)

    ! Get parameters not stored in file

    ! Note:
    ! parameters for reflectance to radiance conversion for the visible channels
    ! see NOAA Polar Orbiter Data User's Guide Section 3.3.2
    ! (http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c3/sec3-3.htm)
    ! Eq 3.3.2-5 and Table 3.3.2-2
    ! 
    ! Satellite         W1      F1      W2      F2
    ! TIROS-N         0.325   443.3   0.303   313.5
    ! NOAA-6          0.109   179.0   0.223   233.7
    ! NOAA-7          0.108   177.5   0.249   261.9
    ! NOAA-8          0.113   183.4   0.230   242.8
    ! NOAA-9          0.117   191.3   0.239   251.8
    ! NOAA-10         0.108   178.8   0.222   231.5
    ! NOAA-11         0.113   184.1   0.229   241.1
    ! NOAA-12         0.124   200.1   0.219   229.9
    ! NOAA-13         0.121   194.09  0.243   249.42
    ! NOAA-14         0.136   221.42  0.245   252.29
    ! 
    ! I_lambda = reflectance / 100 * F / pi / w
    ! I_lambda is the radiance W m-2 sr-1 um-1
    ! reflectance is in percent
    ! F is the band-integrated extraterrestrial solar irradiance at normal incidence
    !   at the top of the atmosphere at mean Earth-Sun distance W m-2
    ! w is the band equivalent width um

    SELECT CASE(INT(Header(1)))
    CASE(1)
       CALL Parse_Date_Time( Header, 3, year, dayno, hours, UTC_msecs )
       IF( 1987 .lt. year )THEN
          AVHRR_No = 11
          write(*,*)'Dataset is from NOAA-11'
          ! Constant values for PRT conversion 
          Coefs%prtTempCoefs1(1) = 276.597
          Coefs%prtTempCoefs2(1) = 0.051275
          Coefs%prtTempCoefs3(1) = 1.363E-6
          Coefs%prtTempCoefs4(1) = 0.
          Coefs%prtTempCoefs5(1) = 0.       
          Coefs%prtTempCoefs1(2) = 276.597  
          Coefs%prtTempCoefs2(2) = 0.051275 
          Coefs%prtTempCoefs3(2) = 1.363E-6 
          Coefs%prtTempCoefs4(2) = 0.       
          Coefs%prtTempCoefs5(2) = 0.       
          Coefs%prtTempCoefs1(3) = 276.597  
          Coefs%prtTempCoefs2(3) = 0.051275 
          Coefs%prtTempCoefs3(3) = 1.363E-6 
          Coefs%prtTempCoefs4(3) = 0.       
          Coefs%prtTempCoefs5(3) = 0.       
          Coefs%prtTempCoefs1(4) = 276.597  
          Coefs%prtTempCoefs2(4) = 0.051275 
          Coefs%prtTempCoefs3(4) = 1.363E-6 
          Coefs%prtTempCoefs4(4) = 0.       
          Coefs%prtTempCoefs5(4) = 0.       
          ! Calibration data 
          Coefs%Nspace(1) = 0.
          Coefs%Nspace(2) = 0.
          Coefs%Nspace(3) = 0. 
          Coefs%nonLinearCoefs3(1) = 0.
          Coefs%nonLinearCoefs3(2) = 0.
          Coefs%nonLinearCoefs3(3) = 0.
          Coefs%nonLinearCoefs4(1) = 0.
          Coefs%nonLinearCoefs4(2) = 0.
          Coefs%nonLinearCoefs4(3) = 0.
          Coefs%nonLinearCoefs5(1) = 0.
          Coefs%nonLinearCoefs5(2) = 0.
          Coefs%nonLinearCoefs5(3) = 0.
          ! Radiance->temp conversion
          Coefs%nuC_numTemps = 3
          Coefs%nuC_minT = (/180.,225.,275.,0./)
          Coefs%nuC_maxT = (/225.,275.,320.,0./)
          Coefs%nuC_Array = RESHAPE(SOURCE=(/2663.50,2668.15,2671.40,0.,&
               926.81,927.36,927.36,0.,841.40,841.81,842.20,0./),&
               SHAPE=(/4,3/))
          ! Non-linear correction lookup table
          Coefs%number_Scene_Temps = 14
          Coefs%number_Target_Temps = 3
          Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
               (/320,315,310,305,295,285,275,265,255,245,235,225,215,205/)
          Coefs%Target_Temp(1:Coefs%number_Target_Temps) = &
               (/10,15,20/)
          Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
               1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
               (/4.29,3.50,2.85,2.23,1.05,0.24,-0.45,-1.06,-1.41,&
               -1.70,-1.87,-1.90,-1.82,-1.54,3.71,2.98,2.33,1.73,&
               0.68,-0.21,-0.79,-1.37,-1.72,-1.96,-2.10,-2.14,&
               -2.02,-1.76,3.25,2.55,1.91,1.32,0.22,-0.67,-1.15,&
               -1.66,-2.03,-2.22,-2.28,-2.36,-2.20,-1.98/),&
               SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
          Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
               1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
               (/1.43,1.23,1.05,0.85,0.43,0.07,-0.19,-0.37,&
               -0.60,-0.72,-0.84,-0.94,-1.12,-1.15,1.26,1.03,&
               0.84,0.64,0.28,-0.07,-0.34,-0.51,-0.77,-0.90,&
               -1.02,-1.06,-1.24,-1.27,1.12,0.89,0.70,0.47,&
               0.09,-0.23,-0.47,-0.60,-0.78,-0.92,-1.00,&
               -1.16,-1.16,-1.23/),&
               SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
          have_12micron = .TRUE.
          ! Reflectance to radiance conversion
          Coefs%F(1)  = 184.1
          Coefs%w(1) = 0.113
          Coefs%F(2)  = 241.1
          Coefs%w(2) = 0.229
          Coefs%F(3)  = 0.0 ! no channel 3A
          Coefs%w(3) = 1.0
       ELSE
          AVHRR_No = 1
          write(*,*)'Dataset is from TIROS-N'
          ! Constant values for PRT conversion 
          Coefs%prtTempCoefs1(1) = 277.73
          Coefs%prtTempCoefs2(1) = 0.047752
          Coefs%prtTempCoefs3(1) = 8.29E-6
          Coefs%prtTempCoefs4(1) = 0.
          Coefs%prtTempCoefs5(1) = 0.
          Coefs%prtTempCoefs1(2) = 277.41
          Coefs%prtTempCoefs2(2) = 0.046637
          Coefs%prtTempCoefs3(2) = 11.01E-6
          Coefs%prtTempCoefs4(2) = 0.
          Coefs%prtTempCoefs5(2) = 0.
          Coefs%prtTempCoefs1(3) = 277.14
          Coefs%prtTempCoefs2(3) = 0.045188
          Coefs%prtTempCoefs3(3) = 14.77E-6
          Coefs%prtTempCoefs4(3) = 0.
          Coefs%prtTempCoefs5(3) = 0.
          Coefs%prtTempCoefs1(4) = 277.42
          Coefs%prtTempCoefs2(4) = 0.046387
          Coefs%prtTempCoefs3(4) = 10.59E-6
          Coefs%prtTempCoefs4(4) = 0.
          Coefs%prtTempCoefs5(4) = 0.
          ! Calibration data 
          Coefs%Nspace(1) = 0.
          Coefs%Nspace(2) = -1.151
          Coefs%Nspace(3) = -1.151 
          Coefs%nonLinearCoefs3(1) = 0.
          Coefs%nonLinearCoefs3(2) = 0.
          Coefs%nonLinearCoefs3(3) = 0.
          Coefs%nonLinearCoefs4(1) = 0.
          Coefs%nonLinearCoefs4(2) = 0.
          Coefs%nonLinearCoefs4(3) = 0.
          Coefs%nonLinearCoefs5(1) = 0.
          Coefs%nonLinearCoefs5(2) = 0.
          Coefs%nonLinearCoefs5(3) = 0.
          ! Radiance->temp conversion
          Coefs%nuC_numTemps = 3
          Coefs%nuC_minT = (/180.,225.,275.,0./)
          Coefs%nuC_maxT = (/225.,275.,320.,0./)
          Coefs%nuC_Array = RESHAPE(SOURCE=(/2631.81,2631.81,2631.81,0.,&
               911.13,911.54,912.01,0.,0.,0.,0.,0./),SHAPE=(/4,3/))
          ! Non-linear correction lookup table
          Coefs%number_Scene_Temps = 9
          Coefs%number_Target_Temps = 1
          Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
               (/304.9,294.9,285.0,275.1,264.9,255.1,234.9,224.9,204.9/)
          Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
               1:Coefs%number_Target_Temps) = &
               RESHAPE(SOURCE=&
               (/1.25,0.98,0.0,-0.03,-0.08,-0.1,-0.75,-0.95,-1.67/),&
               SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
          have_12micron = .FALSE.
          ! Reflectance to radiance conversion
          Coefs%F(1)  = 443.3
          Coefs%w(1) = 0.325
          Coefs%F(2)  = 313.5
          Coefs%w(2) = 0.303
          Coefs%F(3)  = 0.0 ! no channel 3A
          Coefs%w(3) = 1.0
       ENDIF
    CASE(2)
       IF( 1990 .lt. year )THEN
          AVHRR_No = 13
          write(*,*)'Dataset is from NOAA-13'
          Coefs%prtTempCoefs1(1) = 276.597
          Coefs%prtTempCoefs2(1) = 0.051275
          Coefs%prtTempCoefs3(1) = 1.363E-6
          Coefs%prtTempCoefs4(1) = 0.
          Coefs%prtTempCoefs5(1) = 0.
          Coefs%prtTempCoefs1(2) = 276.597
          Coefs%prtTempCoefs2(2) = 0.051275
          Coefs%prtTempCoefs3(2) = 1.363E-6
          Coefs%prtTempCoefs4(2) = 0.
          Coefs%prtTempCoefs5(2) = 0.       
          Coefs%prtTempCoefs1(3) = 276.597  
          Coefs%prtTempCoefs2(3) = 0.051275 
          Coefs%prtTempCoefs3(3) = 1.363E-6 
          Coefs%prtTempCoefs4(3) = 0.       
          Coefs%prtTempCoefs5(3) = 0.       
          Coefs%prtTempCoefs1(4) = 276.597  
          Coefs%prtTempCoefs2(4) = 0.051275 
          Coefs%prtTempCoefs3(4) = 1.363E-6 
          Coefs%prtTempCoefs4(4) = 0.       
          Coefs%prtTempCoefs5(4) = 0.       
          ! Calibration data - 'new' calibration non-linear proceedure
          Coefs%Nspace(1) = 0.
          Coefs%Nspace(2) = -5.31
          Coefs%Nspace(3) = -3.28
          Coefs%nonLinearCoefs3(1) = 0.
          Coefs%nonLinearCoefs3(2) = 0.
          Coefs%nonLinearCoefs3(3) = 0.
          Coefs%nonLinearCoefs4(1) = 5.01
          Coefs%nonLinearCoefs4(2) = 0.91159
          Coefs%nonLinearCoefs4(3) = 0.0003820
          Coefs%nonLinearCoefs5(1) = 3.24
          Coefs%nonLinearCoefs5(2) = 0.94784
          Coefs%nonLinearCoefs5(3) = 0.0002057
          ! Radiance->temp conversion
          Coefs%nuC_numTemps = 4
          Coefs%nuC_minT = (/180.,230.,270.,290./)
          Coefs%nuC_maxT = (/225.,270.,310.,330./)
          Coefs%nuC_Array = RESHAPE(SOURCE=&
               (/2636.124,2640.147,2643.153,2644.382,&
               924.0114,924.5165,924.9732,925.2164,&
               836.1164,836.4339,836.7651,836.9520/),SHAPE=(/4,3/))
          ! Ignore lookup table since has newer calibration coefs listed
          Coefs%number_Scene_Temps = 0
          Coefs%number_Target_Temps = 0
          have_12micron = .TRUE.
          ! Reflectance to radiance conversion
          Coefs%F(1)  = 194.09
          Coefs%w(1) = 0.121
          Coefs%F(2)  = 249.42
          Coefs%w(2) = 0.243
          Coefs%F(3)  = 0.0 ! no channel 3A
          Coefs%w(3) = 1.0
       ELSE
          AVHRR_No = 6
          write(*,*)'Dataset is from NOAA-6'
          Coefs%prtTempCoefs1(1) = 278.3863101
          Coefs%prtTempCoefs2(1) = 3.1620496E-2
          Coefs%prtTempCoefs3(1) = 4.9133034E-5
          Coefs%prtTempCoefs4(1) = 0.
          Coefs%prtTempCoefs5(1) = 0.
          Coefs%prtTempCoefs1(2) = 278.2352898
          Coefs%prtTempCoefs2(2) = 3.0384174E-2
          Coefs%prtTempCoefs3(2) = 5.1555794E-5 
          Coefs%prtTempCoefs4(2) = 0.
          Coefs%prtTempCoefs5(2) = 0.
          Coefs%prtTempCoefs1(3) = 278.0691316
          Coefs%prtTempCoefs2(3) = 4.2528998E-2
          Coefs%prtTempCoefs3(3) = 1.6065146E-5
          Coefs%prtTempCoefs4(3) = 0.
          Coefs%prtTempCoefs5(3) = 0.
          Coefs%prtTempCoefs1(4) = 277.6288372
          Coefs%prtTempCoefs2(4) = 4.0905453E-2
          Coefs%prtTempCoefs3(4) = 1.9771519E-5
          Coefs%prtTempCoefs4(4) = 0.
          Coefs%prtTempCoefs5(4) = 0.
          ! Calibration data 
          Coefs%Nspace(1) = 0.
          Coefs%Nspace(2) = -2.18222
          Coefs%Nspace(3) = -2.18222
          Coefs%nonLinearCoefs3(1) = 0.
          Coefs%nonLinearCoefs3(2) = 0.
          Coefs%nonLinearCoefs3(3) = 0.
          Coefs%nonLinearCoefs4(1) = 0.
          Coefs%nonLinearCoefs4(2) = 0.
          Coefs%nonLinearCoefs4(3) = 0.
          Coefs%nonLinearCoefs5(1) = 0.
          Coefs%nonLinearCoefs5(2) = 0.
          Coefs%nonLinearCoefs5(3) = 0.
          ! Radiance->temp conversion
          Coefs%nuC_numTemps = 3
          Coefs%nuC_minT = (/180.,225.,275.,0./)
          Coefs%nuC_maxT = (/225.,275.,320.,0./)
          Coefs%nuC_Array = RESHAPE(SOURCE=(/2649.90,2653.90,2658.05,0.,&
               910.72,911.41,912.14,0.,0.,0.,0.,0./),SHAPE=(/4,3/))
          ! Non-linear correction lookup table
          Coefs%number_Scene_Temps = 13
          Coefs%number_Target_Temps = 1
          Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
               (/315,305,295,285,275,255,245,235,225,215,205,195,185/)
          Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
               1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
               (/0.8,0.5,0.3,0.0,-0.4,-0.8,-1.4,-1.4,-2.0,-2.0,-2.8,&
               -2.6,-2.0/),&
               SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
          have_12micron = .FALSE.
          ! Reflectance to radiance conversion
          Coefs%F(1)  = 179.0
          Coefs%w(1) = 0.109
          Coefs%F(2)  = 233.7
          Coefs%w(2) = 0.223
          Coefs%F(3)  = 0.0 ! no channel 3A
          Coefs%w(3) = 1.0
       ENDIF
    CASE(4)
       AVHRR_No = 7
       write(*,*)'Dataset is from NOAA-7'
       Coefs%prtTempCoefs1(1) = 277.099
       Coefs%prtTempCoefs2(1) = 5.048E-2
       Coefs%prtTempCoefs3(1) = 2.823E-6
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.
       Coefs%prtTempCoefs1(2) = 276.734
       Coefs%prtTempCoefs2(2) = 5.069E-2
       Coefs%prtTempCoefs3(2) = 2.493E-6
       Coefs%prtTempCoefs4(2) = 0.
       Coefs%prtTempCoefs5(2) = 0.
       Coefs%prtTempCoefs1(3) = 276.876
       Coefs%prtTempCoefs2(3) = 5.148E-2
       Coefs%prtTempCoefs3(3) = 1.040E-6
       Coefs%prtTempCoefs4(3) = 0.
       Coefs%prtTempCoefs5(3) = 0.
       Coefs%prtTempCoefs1(4) = 276.160
       Coefs%prtTempCoefs2(4) = 5.128E-2
       Coefs%prtTempCoefs3(4) = 1.414E-6
       Coefs%prtTempCoefs4(4) = 0.
       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -1.176
       Coefs%Nspace(3) = -1.346
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 0.
       Coefs%nonLinearCoefs4(2) = 0.
       Coefs%nonLinearCoefs4(3) = 0.
       Coefs%nonLinearCoefs5(1) = 0.
       Coefs%nonLinearCoefs5(2) = 0.
       Coefs%nonLinearCoefs5(3) = 0.
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 3
       Coefs%nuC_minT = (/180.,225.,275.,0./)
       Coefs%nuC_maxT = (/225.,275.,320.,0./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/0.,2670.3,2671.9,0.,&
            926.20,926.80,927.22,0.,840.100,840.500,840.872,0./),&
            SHAPE=(/4,3/))
       ! Non-linear correction lookup table
       Coefs%number_Scene_Temps = 9
       Coefs%number_Target_Temps = 1
       Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
            (/315,305,295,285,275,255,235,225,205/)
       Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/1.66,1.05,0.49,0.0,-0.38,-0.66,-0.73,-0.61,-0.19/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/1.08,0.64,0.31,0.0,-0.22,-0.55,-0.86,-0.71,-0.86/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       have_12micron = .TRUE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 177.5
       Coefs%w(1) = 0.108
       Coefs%F(2)  = 261.9
       Coefs%w(2) = 0.249
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE(6)
       AVHRR_No = 8
       write(*,*)'Dataset is from NOAA-8'
       Coefs%prtTempCoefs1(1) = 276.585
       Coefs%prtTempCoefs2(1) = 0.05136
       Coefs%prtTempCoefs3(1) = -9.99E-8
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.
       Coefs%prtTempCoefs1(2) = 276.605
       Coefs%prtTempCoefs2(2) = 0.05122
       Coefs%prtTempCoefs3(2) = 6.86E-8
       Coefs%prtTempCoefs4(2) = 0.
       Coefs%prtTempCoefs5(2) = 0.
       Coefs%prtTempCoefs1(3) = 276.591
       Coefs%prtTempCoefs2(3) = 0.05133
       Coefs%prtTempCoefs3(3) = -1.381E-7
       Coefs%prtTempCoefs4(3) = 0.
       Coefs%prtTempCoefs5(3) = 0.
       Coefs%prtTempCoefs1(4) = 276.592
       Coefs%prtTempCoefs2(4) = 0.05133
       Coefs%prtTempCoefs3(4) = -1.489E-7
       Coefs%prtTempCoefs4(4) = 0.
       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -2.784
       Coefs%Nspace(3) = -2.784
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 0.
       Coefs%nonLinearCoefs4(2) = 0.
       Coefs%nonLinearCoefs4(3) = 0.
       Coefs%nonLinearCoefs5(1) = 0.
       Coefs%nonLinearCoefs5(2) = 0.
       Coefs%nonLinearCoefs5(3) = 0.
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 3
       Coefs%nuC_minT = (/180.,225.,275.,0./)
       Coefs%nuC_maxT = (/225.,275.,320.,0./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/2631.52,2636.05,2639.18,0.,&
            913.360,913.865,914.305,0.,0.,0.,0.,0./),SHAPE=(/4,3/))
       ! Non-linear correction lookup table
       Coefs%number_Scene_Temps = 9
       Coefs%number_Target_Temps = 1
       Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
            (/315,305,295,285,275,255,235,225,205/)
       Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/0.8,0.3,-0.1,-0.3,-0.4,-0.4,0.2,0.7,2.2/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       have_12micron = .FALSE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 183.4
       Coefs%w(1) = 0.113
       Coefs%F(2)  = 242.8
       Coefs%w(2) = 0.230
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE(7)
       AVHRR_No = 9
       write(*,*)'Dataset is from NOAA-9'
       Coefs%prtTempCoefs1(1) = 277.018
       Coefs%prtTempCoefs2(1) = 0.05128
       Coefs%prtTempCoefs3(1) = 0.
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.
       Coefs%prtTempCoefs1(2) = 276.750
       Coefs%prtTempCoefs2(2) = 0.05128
       Coefs%prtTempCoefs3(2) = 0.
       Coefs%prtTempCoefs4(2) = 0.
       Coefs%prtTempCoefs5(2) = 0.
       Coefs%prtTempCoefs1(3) = 276.862
       Coefs%prtTempCoefs2(3) = 0.05128
       Coefs%prtTempCoefs3(3) = 0.
       Coefs%prtTempCoefs4(3) = 0.
       Coefs%prtTempCoefs5(3) = 0.
       Coefs%prtTempCoefs1(4) = 276.546
       Coefs%prtTempCoefs2(4) = 0.05128
       Coefs%prtTempCoefs3(4) = 0.
       Coefs%prtTempCoefs4(4) = 0.
       Coefs%prtTempCoefs5(4) = 0.
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = -3.384
       Coefs%Nspace(3) = -2.313 
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 0.
       Coefs%nonLinearCoefs4(2) = 0.
       Coefs%nonLinearCoefs4(3) = 0.
       Coefs%nonLinearCoefs5(1) = 0.
       Coefs%nonLinearCoefs5(2) = 0.
       Coefs%nonLinearCoefs5(3) = 0.
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 3
       Coefs%nuC_minT = (/180.,225.,275.,0./)
       Coefs%nuC_maxT = (/225.,275.,320.,0./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/2670.93,2674.81,2678.11,0.,&
            928.50,929.02,929.46,0.,844.41,844.80,845.19,0./),&
            SHAPE=(/4,3/))
       ! Non-linear correction lookup table
       Coefs%number_Scene_Temps = 9
       Coefs%number_Target_Temps = 1
       Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
            (/315,305,295,285,275,255,235,225,205/)
       Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/1.8,0.9,0.2,-0.4,-0.9,-1.4,-1.6,-1.5,-1.0/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/1.0,0.6,0.2,-0.1,-0.5,-0.8,-1.0,-1.3,-1.4/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       have_12micron = .TRUE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 191.3
       Coefs%w(1) = 0.117
       Coefs%F(2)  = 251.8
       Coefs%w(2) = 0.239
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE(8)
       AVHRR_No = 10
       write(*,*)'Dataset is from NOAA-10'
       Coefs%prtTempCoefs1(1) = 276.659
       Coefs%prtTempCoefs2(1) = 0.051275
       Coefs%prtTempCoefs3(1) = 1.363E-6
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.       
       Coefs%prtTempCoefs1(2) = 276.659  
       Coefs%prtTempCoefs2(2) = 0.051275 
       Coefs%prtTempCoefs3(2) = 1.363E-6 
       Coefs%prtTempCoefs4(2) = 0.       
       Coefs%prtTempCoefs5(2) = 0.       
       Coefs%prtTempCoefs1(3) = 276.659  
       Coefs%prtTempCoefs2(3) = 0.051275 
       Coefs%prtTempCoefs3(3) = 1.363E-6 
       Coefs%prtTempCoefs4(3) = 0.       
       Coefs%prtTempCoefs5(3) = 0.       
       Coefs%prtTempCoefs1(4) = 276.659  
       Coefs%prtTempCoefs2(4) = 0.051275 
       Coefs%prtTempCoefs3(4) = 1.363E-6 
       Coefs%prtTempCoefs4(4) = 0.       
       Coefs%prtTempCoefs5(4) = 0.       
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = 0.
       Coefs%Nspace(3) = 0.
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 0.
       Coefs%nonLinearCoefs4(2) = 0.
       Coefs%nonLinearCoefs4(3) = 0.
       Coefs%nonLinearCoefs5(1) = 0.
       Coefs%nonLinearCoefs5(2) = 0.
       Coefs%nonLinearCoefs5(3) = 0.
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 3
       Coefs%nuC_minT = (/180.,225.,275.,0./)
       Coefs%nuC_maxT = (/225.,275.,320.,0./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/2658.53,2657.60,2660.76,0.,&
            908.73,909.18,909.58,0.,0.,0.,0.,0./),SHAPE=(/4,3/))
       ! Non-linear correction lookup table
       Coefs%number_Scene_Temps = 13
       Coefs%number_Target_Temps = 3
       Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
            (/320,315,305,295,285,275,265,255,245,235,225,215,205/)
       Coefs%Target_Temp(1:Coefs%number_Target_Temps) = &
            (/10,15,20/)
       Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/3.50,2.93,2.93,1.12,0.20,-0.46,-0.76,-1.33,-1.74,-1.79,-2.22,&
            -2.58,-2.47,2.83,2.19,1.34,0.57,-0.15,-0.53,-0.93,-1.49,-2.09,&
            -2.20,-2.51,-2.65,-2.88,2.54,1.97,1.11,0.12,-0.38,-1.08,-1.37,&
            -1.77,-2.26,-2.53,-2.53,-2.80,-3.27/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       have_12micron = .FALSE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 178.8
       Coefs%w(1) = 0.108
       Coefs%F(2)  = 231.5
       Coefs%w(2) = 0.222
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE(5)
       AVHRR_No = 12
       write(*,*)'Dataset is from NOAA-12'
       Coefs%prtTempCoefs1(1) = 276.597
       Coefs%prtTempCoefs2(1) = 0.051275
       Coefs%prtTempCoefs3(1) = 1.363E-6
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.       
       Coefs%prtTempCoefs1(2) = 276.597  
       Coefs%prtTempCoefs2(2) = 0.051275 
       Coefs%prtTempCoefs3(2) = 1.363E-6 
       Coefs%prtTempCoefs4(2) = 0.       
       Coefs%prtTempCoefs5(2) = 0.       
       Coefs%prtTempCoefs1(3) = 276.597  
       Coefs%prtTempCoefs2(3) = 0.051275 
       Coefs%prtTempCoefs3(3) = 1.363E-6 
       Coefs%prtTempCoefs4(3) = 0.       
       Coefs%prtTempCoefs5(3) = 0.       
       Coefs%prtTempCoefs1(4) = 276.597  
       Coefs%prtTempCoefs2(4) = 0.051275 
       Coefs%prtTempCoefs3(4) = 1.363E-6 
       Coefs%prtTempCoefs4(4) = 0.       
       Coefs%prtTempCoefs5(4) = 0.       
       ! Calibration data 
       Coefs%Nspace(1) = 0.
       Coefs%Nspace(2) = 0.
       Coefs%Nspace(3) = 0.
       Coefs%nonLinearCoefs3(1) = 0.
       Coefs%nonLinearCoefs3(2) = 0.
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 0.
       Coefs%nonLinearCoefs4(2) = 0.
       Coefs%nonLinearCoefs4(3) = 0.
       Coefs%nonLinearCoefs5(1) = 0.
       Coefs%nonLinearCoefs5(2) = 0.
       Coefs%nonLinearCoefs5(3) = 0.
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 4
       Coefs%nuC_minT = (/180.,230.,270.,290./)
       Coefs%nuC_maxT = (/225.,270.,310.,330./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/2632.713,2636.669,2639.61,2640.817,&
            920.0158,920.5504,921.0291,921.2741,&
            836.6847,837.0251,837.3641,837.5612/),SHAPE=(/4,3/))
       ! Look up tables for BT non-linear correction
       Coefs%number_Scene_Temps = 14
       Coefs%number_Target_Temps = 4
       Coefs%Scene_Temp(1:Coefs%number_Scene_Temps) = &
            (/320,315,310,305,295,285,275,265,255,245,235,225,215,205/)
       Coefs%Target_Temp(1:Coefs%number_Target_Temps) = &
            (/10,15,20,25/)
       Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/3.21,2.58,2.04,1.6,0.8,0.16,-0.41,-0.071,&
            -1.04,-1.18,-1.05,-1.33,-1.24,-1.58,2.88,&
            2.39,1.94,1.42,0.53,-0.23,-0.84,-0.97,&
            -1.2,-1.4,-1.59,-1.65,-1.65,-1.8,2.27,&
            1.72,1.28,0.8,0.13,-0.52,-1.05,-1.19,&
            -1.53,-1.58,-1.51,-1.58,-1.49,-1.31,&
            1.91,1.43,0.98,0.52,-0.16,-0.7,-1.19,-1.32,&
            -1.59,-1.62,-1.63,-1.67,-1.53,-1.33/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
            1:Coefs%number_Target_Temps) = RESHAPE(SOURCE=&
            (/0.8,0.8,0.8,0.73,0.37,0.08,-0.21,-0.37,&
            -0.47,-0.63,-0.88,-1.01,-1.15,-1.17,0.8,&
            0.8,0.73,0.61,0.18,-0.08,-0.31,-0.41,-0.53,&
            -0.76,-0.94,-1.1,-1.19,-1.16,0.8,0.73,0.61,&
            0.37,0.08,-0.21,-0.37,-0.47,-0.63,-0.88,-1.01,&
            -1.15,-1.17,-1.19,0.73,0.61,0.37,0.18,-0.08,&
            -0.31,-0.41,-0.53,-0.76,-0.94,-1.1,-1.19,-1.16,&
            -1.23/),&
            SHAPE=(/Coefs%number_Scene_Temps,Coefs%number_Target_Temps/))
       have_12micron = .TRUE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 200.1
       Coefs%w(1) = 0.124
       Coefs%F(2)  = 229.9
       Coefs%w(2) = 0.219
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE(3)
       AVHRR_No = 14
       write(*,*)'Dataset is from NOAA-14'
       Coefs%prtTempCoefs1(1) = 276.597
       Coefs%prtTempCoefs2(1) = 0.051275
       Coefs%prtTempCoefs3(1) = 1.363E-6
       Coefs%prtTempCoefs4(1) = 0.
       Coefs%prtTempCoefs5(1) = 0.       
       Coefs%prtTempCoefs1(2) = 276.597  
       Coefs%prtTempCoefs2(2) = 0.051275 
       Coefs%prtTempCoefs3(2) = 1.363E-6 
       Coefs%prtTempCoefs4(2) = 0.       
       Coefs%prtTempCoefs5(2) = 0.       
       Coefs%prtTempCoefs1(3) = 276.597  
       Coefs%prtTempCoefs2(3) = 0.051275 
       Coefs%prtTempCoefs3(3) = 1.363E-6 
       Coefs%prtTempCoefs4(3) = 0.       
       Coefs%prtTempCoefs5(3) = 0.       
       Coefs%prtTempCoefs1(4) = 276.597  
       Coefs%prtTempCoefs2(4) = 0.051275 
       Coefs%prtTempCoefs3(4) = 1.363E-6 
       Coefs%prtTempCoefs4(4) = 0.       
       Coefs%prtTempCoefs5(4) = 0.       
       ! Calibration data 
       Coefs%Nspace(1) = 0.0069
       Coefs%Nspace(2) = -4.05
       Coefs%Nspace(3) = -2.29
       Coefs%nonLinearCoefs3(1) = -0.0031
       Coefs%nonLinearCoefs3(2) = 1.00359
       Coefs%nonLinearCoefs3(3) = 0.
       Coefs%nonLinearCoefs4(1) = 3.72
       Coefs%nonLinearCoefs4(2) = 0.92378
       Coefs%nonLinearCoefs4(3) = 0.0003822
       Coefs%nonLinearCoefs5(1) = 2.00
       Coefs%nonLinearCoefs5(2) = 0.96194
       Coefs%nonLinearCoefs5(3) = 0.0001742
       ! Radiance->temp conversion
       Coefs%nuC_numTemps = 4
       Coefs%nuC_minT = (/180.,230.,270.,290./)
       Coefs%nuC_maxT = (/225.,270.,310.,330./)
       Coefs%nuC_Array = RESHAPE(SOURCE=(/2638.652,2642.807,2645.899,2647.169,&
            928.2603,928.8284,929.3323,929.5878,&
            834.4496,834.8066,835.1647,835.374/),SHAPE=(/4,3/))
       Coefs%number_Scene_Temps = 0
       Coefs%number_Target_Temps = 0
       have_12micron = .TRUE.
       ! Reflectance to radiance conversion
       Coefs%F(1)  = 221.42
       Coefs%w(1) = 0.136
       Coefs%F(2)  = 252.29
       Coefs%w(2) = 0.245
       Coefs%F(3)  = 0.0 ! no channel 3A
       Coefs%w(3) = 1.0
    CASE DEFAULT
       call out_ERROR('Invalid Satellite ID',&
            'Read_Data_HeaderV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')       
    END SELECT

    ! pre-KLM l1b files are based on physical records (on the original
    ! tapes), with 2 logical records per physical record.  The header and
    ! each scan are one logical record long.  So the file structure is
    ! 
    ! 'physical' record logical record scanline pos scanline num 
    !          1               1            -            -        header
    !          1               2            -            -        fill (=last scan)
    !          2               3            1            1        scan
    !          .
    !          .
    !          .
    !          n             2n-1        nScans      >=nScans     scan
    !          n              2n            -            -        fill
    ! (for odd nScans)
    ! or
    ! 
    !          n             2n-1       nScans-1    >=nScans-1    scan
    !          n              2n         nScans      >=nScans     scan
    ! (for even nScans)
    ! 
    ! KLM l1b files are based on physical records perhaps but have 1 logical
    ! record per physical record.  The header and each scan are one logical
    ! record long.  So the file structure is
    ! 
    ! 'physical' record logical record scanline pos scanline num 
    !          1               1            -            -        header
    !          2               2            1            1        scan
    !          .
    !          .
    !          .
    !          n               n         nScans      >=nScans     scan
    ! 
    ! scanline num can be greater than scanline pos is there are data gaps
    ! (It might be better to output the data with the data gaps filled with
    ! bad data---for the future).

    ! So read the record after the header to get in the right place

#ifdef USE_GZFILE
    STAT = gzread(Unit,Header)
    IF(STAT.LE.0)THEN
       CALL Check_IOS(1,'Reading Header','Read_Data_HeaderV1',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#else
    READ(Unit,IOSTAT=STAT)Header
    CALL Check_IOS(STAT,'Reading Header','Read_Data_Headerv1',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

    ! Free memory
    DEALLOCATE(Header)
 
  END SUBROUTINE Read_Data_Headerv1
!############### HEADER v3 ##############################################   
 
 

   


 

END MODULE module_read_data_header

