MODULE NOAA_LoadAVHRRLevel1B
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsPath
  USE GbcsErrorHandler
  USE GbcsStringUtil
  USE type

  use module_filtering_thresholds
  use module_filtering_thresholds_prt
  use module_average_per_scanline
  use module_average_per_orbite
  use module_urict
  use module_sensitivities      
  use module_radiance_uncertainties
#ifdef USE_HDF5
  USE HDF5
#endif
#ifdef USE_GZFILE
  USE gzfile
#endif
#ifdef USE_IEEE
  USE ieee_arithmetic
#endif

  IMPLICIT NONE
  ! Interfaces
  INTERFACE Reallocate_Arrays
     MODULE PROCEDURE Reallocate_Arrays_Double1D
     MODULE PROCEDURE Reallocate_Arrays_Real1D
     MODULE PROCEDURE Reallocate_Arrays_Int1D
     MODULE PROCEDURE Reallocate_Arrays_Int_2D
     MODULE PROCEDURE Reallocate_Arrays_Real2D
     MODULE PROCEDURE Reallocate_Arrays_Int2_2D
     MODULE PROCEDURE Reallocate_Arrays_Logical1D
  END INTERFACE
  INTERFACE Reallocate_Final_Arrays
     MODULE PROCEDURE Reallocate_Final_Arrays_Double1D
     MODULE PROCEDURE Reallocate_Final_Arrays_Real1D
     MODULE PROCEDURE Reallocate_Final_Arrays_Int1D
     MODULE PROCEDURE Reallocate_Final_Arrays_Int_2D
     MODULE PROCEDURE Reallocate_Final_Arrays_Real2D
     MODULE PROCEDURE Reallocate_Final_Arrays_Int2_2D
     MODULE PROCEDURE Reallocate_Final_Arrays_Logical1D
  END INTERFACE
!!$  INTERFACE Get_Footprint_Angle
!!$     MODULE PROCEDURE Get_Footprint_Angle_Real
!!$     MODULE PROCEDURE Get_Footprint_Angle_Dble
!!$  END INTERFACE

  ! AVHRR data type
  !INTEGER, PARAMETER :: AVHRR_GAC = 1
  !INTEGER, PARAMETER :: AVHRR_LAC = 2
  !INTEGER, PARAMETER :: AVHRR_HRPT = 3
  !INTEGER, PARAMETER :: AVHRR_FRAC = 4

  ! Define Coefficients Structure
  !INTEGER, PARAMETER :: MAX_NUC_TEMP = 4
  !INTEGER, PARAMETER :: NFILTERS = 3
  !INTEGER, PARAMETER :: NCALIB_COEF = 3
  !INTEGER, PARAMETER :: NPRT = 4
  !INTEGER, PARAMETER :: NPRT_IND = 3
  !INTEGER, PARAMETER :: NTEMP_COEFS = 6
  !INTEGER, PARAMETER :: NBBODY_SAMPLES = 10
  !INTEGER, PARAMETER :: MEMORY_ALLOC_STEP = 20000
  !INTEGER, PARAMETER :: NPIXEL_PRT_SMOOTH = 100
  !INTEGER, PARAMETER :: NPIXEL_CNTS_SMOOTH = 27
  !INTEGER, PARAMETER :: NLOOKUP_TABLE_SCENE = 14
  !INTEGER, PARAMETER :: NLOOKUP_TABLE_TARGET = 4

  !INTEGER, PARAMETER :: scanLineStep = 20 ! approx 65 km
  !INTEGER, PARAMETER :: scanElemStep = 3 ! approx 65 km at scan edge
  ! always skip line 1 since it often has location errors (though
  ! many of these were perhaps from the 'two logical records per
  ! physical record': reduced bad tie-point orbits from 4% to 2% (for M2)
  ! often have location errors at end of orbit (N18)
  !INTEGER, PARAMETER :: scanLineStart = 2
  !INTEGER, PARAMETER :: scanElemStart = 1
  ! the maximum timing correction magnitude is 1.75 s in the data that is in
  ! clavr-x (which is the same as the data from UMiami) RESTRICTED to the
  ! dates when there is AVHRR data available.  So expect the shift to be
  ! less than 1.75 * 6 * 1.1 = time/s * lines/s * km/line = 11.50 km, say 20km
  ! as a check.
  !REAL, PARAMETER    :: maxTimingErrorDistance = 20.0

  ! Constants
  !REAL, PARAMETER :: C1 = 1.1910427E-5
  !REAL, PARAMETER :: C2 = 1.4387752
  !REAL, PARAMETER :: eta_ict = 0.985140
  !double precision, parameter :: M_PI = 3.14159265358979323846d0

  !INTEGER, PARAMETER, DIMENSION(13) :: dayNumber = &
  !     (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
  !INTEGER, PARAMETER, DIMENSION(13) :: dayNumberLeap = &
  !     (/0,31,60,91,121,152,182,213,244,274,305,335,366/)

  LOGICAL :: print_missed_matching = .FALSE.
  INTEGER, PARAMETER :: VISCHAN_OPER_YEAR_MIN=2004

  !TYPE AVHRR_Instrument_Coefs
     ! nuC for temperature -> radiance conversion
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: nuC
     ! A for temperature -> radiance conversion
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: Aval
     ! B for temperature -> radiance conversion
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: Bval
     ! F for reflectance -> radiance conversion
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: F
     ! w for reflectance -> radiance conversion
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: w
     ! Older style radiance->temperarure coefficients (pre KLM)
     ! nuC for temperature -> radiance conversion
     ! Number of entries
  !   INTEGER(GbcsInt1) :: nuC_numTemps
     ! minimum temperature range
  !   REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP) :: nuC_minT
     ! maximum temperature range
  !   REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP) :: nuC_maxT
     ! nuC for temperature -> radiance conversion (old way)
  !   REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP,NFILTERS) :: nuC_Array
     ! Space radiance for counts -> radiance
  !   REAL(GbcsReal), DIMENSION(NFILTERS) :: Nspace
     ! detector non-linear coeficients filter 3B
  !   REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs3
     ! detector non-linear coeficients filter 3B
  !   REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs4
     ! detector non-linear coeficients filter 3B
  !   REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs5
     ! Counts -> temperature coefs (0) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs1
     ! Counts -> temperature coefs (1) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs2
     ! Counts -> temperature coefs (2) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs3
     ! Counts -> temperature coefs (3) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs4
     ! Counts -> temperature coefs (4) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs5
     ! Counts -> temperature coefs (4) for NPRTs
  !   REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs6
     ! Counts -> temperature patch
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: patchCoef
     ! Counts -> temperature extended patch
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: patchCoefExt
     ! Counts -> temperature (Radiator) */
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs1
     ! Counts -> temperature (Electronics)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs2
     ! Counts -> temperature (Cooler)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs3
     ! Counts -> temperature (Baseplate)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs4
     ! Counts -> temperature (Motor)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs5
     ! Counts -> temperature (A/D convertor)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs6
     ! Counts -> temperature (Patch)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs7
     ! Counts->Volts->Current (Motor)
  !   REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: motorCurrentCoefs
     ! Lookup table for non-linear correction (old method)
  !   INTEGER :: number_Target_Temps
  !   INTEGER :: number_Scene_Temps
  !   REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE) :: Scene_Temp
  !   REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_TARGET) :: Target_Temp
  !   REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE,NLOOKUP_TABLE_TARGET) :: &
  !        Correction_Factor4
  !   REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE,NLOOKUP_TABLE_TARGET) :: &
  !        Correction_Factor5
  !END TYPE AVHRR_Instrument_Coefs

  !TYPE AVHRR_Data
  !   LOGICAL :: isGAC
  !   LOGICAL :: dataFilled
  !   INTEGER :: AVHRR_No
  !   INTEGER :: nelem
  !   INTEGER :: arraySize
  !   LOGICAL :: filter3a
  !   INTEGER :: start_valid
  !   INTEGER :: stop_valid
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: scanLineNumber
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: badTime
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: badNavigation
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: badCalibration
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: transition3A3B
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Lon
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Lat
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: satZA
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: solZA
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: relAz
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts1
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts2
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts3
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts4
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts5
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array1
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array2
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array3A
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array3B
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array4
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array5
  !   LOGICAL :: walton_there
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array3B_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array4_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: array5_error
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: year
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: month
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: day
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: dayNo
  !   REAL, ALLOCATABLE, DIMENSION(:) :: hours
  !   INTEGER, ALLOCATABLE, DIMENSION(:) :: UTC_msecs
  !   REAL(GbcsDble), ALLOCATABLE, DIMENSION(:) :: time
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt1
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt2
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt1Counts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt2Counts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt3Counts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: prt4Counts
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: prt1CountsAll
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: prt2CountsAll
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: prt3CountsAll
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: prt4CountsAll
  !   REAL, ALLOCATABLE, DIMENSION(:) :: bb3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: bb4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: bb5
  !   REAL, ALLOCATABLE, DIMENSION(:) :: sp3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: sp4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: sp5
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter3
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter4
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter5
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter1
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter2
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter3a
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter3
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter4
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter5
  !   REAL, ALLOCATABLE, DIMENSION(:) :: patch
  !   REAL, ALLOCATABLE, DIMENSION(:) :: patchExtended
  !   REAL, ALLOCATABLE, DIMENSION(:) :: Radiator
  !   REAL, ALLOCATABLE, DIMENSION(:) :: Cooler
  !   REAL, ALLOCATABLE, DIMENSION(:) :: a_d_conv
  !   REAL, ALLOCATABLE, DIMENSION(:) :: motor
  !   REAL, ALLOCATABLE, DIMENSION(:) :: electronics
  !   REAL, ALLOCATABLE, DIMENSION(:) :: baseplate
  !   REAL, ALLOCATABLE, DIMENSION(:) :: motorCurrent
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib1
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib1_2
  !   REAL, ALLOCATABLE, DIMENSION(:) :: calib1_intercept
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib2
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib2_2
  !   REAL, ALLOCATABLE, DIMENSION(:) :: calib2_intercept
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3A
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3A_2
  !   REAL, ALLOCATABLE, DIMENSION(:) :: calib3A_intercept
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib4
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: calib5
  !   LOGICAL :: Clavr_There
  !   INTEGER(GbcsInt2), ALLOCATABLE, DIMENSION(:,:) :: clavr_mask
  !   LOGICAL :: Clavrx_There
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: clavrx_mask ! HDF5 wants integer to read into, even if the
                                                         ! file is byte
  !   INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: clavrx_prb
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_3B
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_4
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_5
  !   LOGICAL :: newCalibration_There
  !   REAL :: orbital_temperature
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib3
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib4
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib5
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt1
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt2
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt1Cnts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt2Cnts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt3Cnts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt4Cnts
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB5
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp3
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp4
  !   REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp5
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: Interpolated
  !   LOGICAL :: solar_contamination_failure
  !   LOGICAl, ALLOCATABLE, DIMENSION(:) :: solar_contamination_3B
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: solar_contamination_4
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: solar_contamination_5
  !   LOGICAL, ALLOCATABLE, DIMENSION(:) :: moon_contamination
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array1
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array2
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3A
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3B
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array4
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array5
  !   REAL, DIMENSION(8) :: new_cal_coefs3
  !   REAL, DIMENSION(8) :: new_cal_coefs4
  !   REAL, DIMENSION(8) :: new_cal_coefs5
  !   REAL, DIMENSION(3) :: gain_stdev
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array1_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array2_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3A_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3B_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array4_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array5_error
  !   REAL, ALLOCATABLE, DIMENSION(:,:) :: noise_cnts
  !   REAL :: earthshine_eta
  !   REAL :: poly_coefs3(3)
  !   REAL :: poly_coefs4(3)
  !   REAL :: poly_coefs5(3)
  !   REAL :: gain_maxdev(3)
  !   LOGICAL :: satelliteAlt_There
  !   REAL, ALLOCATABLE, DIMENSION(:) :: satelliteAltitude
  !END TYPE AVHRR_Data

  !TYPE AVHRR_Bad_Data
  !   LOGICAL :: allow_bad_top
  !   LOGICAL :: allow_bad_nav
  !   LOGICAL :: use_bad_nav
  !   LOGICAL :: allow_bad_time
  !   LOGICAL :: allow_bad_calib
  !END TYPE AVHRR_Bad_Data

  !TYPE AVHRR_Radiative_Coefs
  !   REAL :: nuc(NFILTERS)
  !   REAL :: aVal(NFILTERS)
  !   REAL :: bVal(NFILTERS)
  !END TYPE AVHRR_Radiative_Coefs

  !TYPE AVHRR_Land_Mask
  !   INTEGER :: nelem
  !   INTEGER :: arraySize
  !   INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:,:) :: mask
  !END TYPE AVHRR_Land_Mask

  TYPE Walton_Struct
     INTEGER :: nstored
     REAL, ALLOCATABLE :: sm_bb(:,:)
     REAL, ALLOCATABLE :: sm_sp(:,:)
     REAL, ALLOCATABLE :: sm_prt(:,:)
     REAL :: nuc(3)
     REAL :: aval(3)
     REAL :: bval(3)
     REAL :: Nsp(3)
     REAL :: c(3,3)
  END type Walton_Struct

  ! PRT rotations
  !INTEGER, PARAMETER :: prtRotate_GAC(4) = (/3,1,4,2/)
  !INTEGER, PARAMETER :: prtRotate_LAC(4) = (/1,2,3,4/)

  ! Some global variables
  !LOGICAL, SAVE :: AVHRR_Setup = .FALSE.
  !INTEGER :: ndata
  !INTEGER :: ndata_counts
  !INTEGER :: prtRotate(4)
  !INTEGER :: printOnce(10)
  !INTEGER, PARAMETER :: NDATA_GAC = 409
  !INTEGER, PARAMETER :: NDATA_LAC = 2048
  !INTEGER, PARAMETER :: NDATA_MAX = 2048
  !INTEGER :: prtNumber
  !REAL :: storedPrtTemp(4)
  !REAL :: prtCountsStore(4)
  !REAL :: prtCountsPrevStore(4)
  !REAL :: storedPrtCounts(4)
  !REAL :: prtCountsStoreAll(3,4)
  !REAL :: storedPrtCountsAll(3,4)
  !LOGICAL :: scanAngleSetup = .FALSE.
  !REAL :: scanAngle(NDATA_MAX)
  !REAL :: scanAngleIn(51)
  !INTEGER :: time_yearstart
  !TYPE(AVHRR_Data), TARGET, SAVE :: outDataStore
  !TYPE(DateTime), SAVE :: start_time, end_time
  !LOGICAL :: badDay
  !LOGICAL :: validSolZaDecimal
  !LOGICAL :: validClockDriftInfo
  !LOGICAL :: useL1BClock
  !LOGICAL :: validSBBC
  !REAL, ALLOCATABLE, DIMENSION(:,:) :: clavrLonAnchor, clavrLatAnchor

  !LOGICAL :: Little_Endian = .TRUE.

  PRIVATE
#ifdef USE_GZFILE
  PUBLIC :: Get_Land_Mask
#endif
  PUBLIC :: Get_Tie_Points
  PUBLIC :: Load_Imagery
  PUBLIC :: AVHRR_Data
  PUBLIC :: AVHRR_Land_Mask
  PUBLIC :: AVHRR_Instrument_Coefs
  PUBLIC :: AVHRR_Bad_Data
  PUBLIC :: DeAllocate_outData
!!$  PUBLIC :: AVHRR_Get_Intermediate_Pos
!!$  PUBLIC :: Within_AVHRR_Pixel
  PUBLIC :: AVHRR_Radiative_Coefs
  PUBLIC :: AVHRR_Date
  PUBLIC :: InterpolateLonLat
  PUBLIC :: InterpolateData

  ! Interfaces
  PUBLIC Reallocate_Arrays_Double1D
  PUBLIC Reallocate_Arrays_Real1D
  PUBLIC Reallocate_Arrays_Int1D
  PUBLIC Reallocate_Arrays_Int_2D
  PUBLIC Reallocate_Arrays_Real2D
  PUBLIC Reallocate_Arrays_Int2_2D
  PUBLIC Reallocate_Arrays_Logical1D
  PUBLIC Reallocate_Final_Arrays_Double1D
  PUBLIC Reallocate_Final_Arrays_Real1D
  PUBLIC Reallocate_Final_Arrays_Int1D
  PUBLIC Reallocate_Final_Arrays_Int_2D
  PUBLIC Reallocate_Final_Arrays_Real2D
  PUBLIC Reallocate_Final_Arrays_Int2_2D
  PUBLIC Reallocate_Final_Arrays_Logical1D

  INTEGER, PARAMETER :: ARCHIVE_HEADER_SIZE = 512
  INTEGER, PARAMETER :: ARCHIVE_HEADER_SIZE_POD = 122

  ! -- Module name for error messages
  CHARACTER( * ), PRIVATE, PARAMETER :: MODULE_NAME = &
    'AVHRR/NOAA_LoadAVHRRLevel1B.f90'

CONTAINS

  SUBROUTINE Load_Imagery( IMG, outputData, bad_data, use_new_calibration, &
       calc_new_calibration, use_2000_date, copy_to_output, output_radiances, &
       use_old_37_calibration, getTiePoints, getLandMask, getLocation, &
       use_night_only, correct_solar_contamination, &
       correct_solar_contamination_all, correct_earthshine, no_landmask, &
       moon_events, ignore_solar_contamination, avhrr_rad_coefs, &
       clavrNavPathname, clavrCldPathname, clavrPrbPathname, &
       printAllErrors, doNotCopy, correct_solar_simple, new_vis_cal,&
       calibration_coef_overlap,nedt_scale, out_instr_coefs,&
       use_old_solar,use_iasi_calibration,use_walton)

    !--------------
    !Parameters
    !--------------

    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Load_Imagery'


    TYPE(Imagery), INTENT(INOUT) :: IMG
    TYPE(AVHRR_Data), target, OPTIONAL :: outputData
    TYPE(AVHRR_Bad_Data), OPTIONAL :: bad_data
    LOGICAL, OPTIONAL :: use_new_calibration
    LOGICAL, OPTIONAL :: calc_new_calibration
    LOGICAL, OPTIONAL :: use_2000_date
    LOGICAL, OPTIONAL :: copy_to_output
    LOGICAL, OPTIONAL :: output_radiances
    LOGICAL, OPTIONAL :: use_old_37_calibration
    ! IMG does now have time info, but retain this code for now
    LOGICAL, OPTIONAL :: getTiePoints ! this needs to use outputData since
                                      ! IMG does not have time info
    LOGICAL, OPTIONAL :: getLandMask  ! this needs to use outputData since
                                      ! IMG does not have time info
    LOGICAL, OPTIONAL :: getLocation  ! this needs to use outputData since
                                      ! IMG does not have time info
    LOGICAL, OPTIONAL :: use_night_only
    LOGICAL, OPTIONAL :: correct_solar_contamination
    LOGICAL, OPTIONAL :: correct_solar_contamination_all
    LOGICAL, OPTIONAL :: correct_earthshine
    LOGICAL, OPTIONAL :: no_landmask
    LOGICAL, OPTIONAL :: moon_events
    LOGICAL, OPTIONAL :: ignore_solar_contamination
    TYPE(AVHRR_Radiative_Coefs), OPTIONAL :: avhrr_rad_coefs
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: clavrNavPathname
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: clavrCldPathname
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: clavrPrbPathname
    LOGICAL, OPTIONAL :: printAllErrors
    LOGICAL, OPTIONAL :: doNotCopy
    LOGICAL, OPTIONAL :: correct_solar_simple
    LOGICAL, OPTIONAL :: new_vis_cal
    REAL, OPTIONAL :: calibration_coef_overlap
    REAL, OPTIONAL :: nedt_scale
    TYPE(AVHRR_Instrument_Coefs), INTENT(out), OPTIONAL :: out_instr_coefs
    LOGICAL, INTENT(IN), OPTIONAL :: use_old_solar
    LOGICAL, INTENT(IN), OPTIONAL :: use_iasi_calibration
    LOGICAL, INTENT(IN), OPTIONAL :: use_walton

    ! Local variables

#ifdef USE_GZFILE
    INTEGER(c_int) :: Unit ! for gzfile
#else
    INTEGER :: Unit
#endif
    INTEGER :: Ios
    CHARACTER(LEN=512) :: Filename
    INTEGER :: VersionNo
    INTEGER :: sizeOfRecords
    TYPE(AVHRR_Instrument_Coefs) :: Coefs
    INTEGER :: DataType
    INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:) :: record
    INTEGER :: nScanLines
    TYPE(AVHRR_Bad_Data) :: bad_data_in
    INTEGER :: scanLinePos
    INTEGER :: start_valid
    INTEGER :: stop_valid
    INTEGER :: AVHRR_No
    LOGICAL :: use_new_calibration_in
    LOGICAL :: calc_new_calibration_in
    LOGICAL :: use_old_37_calibration_in
    TYPE(AVHRR_Data), POINTER :: outData => NULL()
    TYPE(AVHRR_Data), POINTER :: DataForFiltering 
    CHARACTER(LEN=512) :: Command
    LOGICAL :: copy_there
    LOGICAL :: out_radiances
    LOGICAL :: have_12micron
    LOGICAL :: load_landmask
    INTEGER :: nScans
    INTEGER :: headerNScans
    INTEGER :: scan
    LOGICAL :: badDayS,badDayE
    LOGICAL :: validSolZaDecimalS,validSolZaDecimalE
    LOGICAL :: validClockDriftInfoS,validClockDriftInfoE
    LOGICAL :: validClockCorrectionS,validClockCorrectionE
    LOGICAL :: validSBBCS,validSBBCE
    LOGICAL :: do_copy
    LOGICAL :: corr_solar_simple
    LOGICAL :: new_vis_calibration
    REAL :: cal_coef_overlap
    REAL :: nedt_scl
    LOGICAL :: use_old_solar_contam
    LOGICAL :: use_iasi_cal
    LOGICAL :: use_walton_in
    
    REAL :: coefs1(8,2)
    REAL :: coefs2(8,2)
    REAL :: coefs3(8,2)
    INTEGER :: ncoefs
    REAL :: coefs_frac1
    REAL :: coefs_frac2
    REAL :: coef_window
    REAL(GbcsDble) :: coefwindow
    LOGICAL :: use_iasi
    INTEGER :: j
    REAL  :: real_fillvalue=9.96921e+36
    INTEGER  :: integer_fillvalue=-32767
    ! WIndow to interpolated case in seconds
    coefwindow = coef_window*30.*86400.

    ! Is this a little endian or big endian machine (it matters...)
    Little_Endian = .NOT.Big_Endian()

    ! New file - initialise some things
    scanLinePos = 0
    sizeOfRecords = 0
    ! top,nav,time,cal,cal3B,cal4,cal5,prtbad,prtmarg
    printOnce = (/0,0,0,0,0,0,0,0,0,0/)
    prtNumber = -1
    storedPrtTemp = (/0.,0.,0.,0./)
    prtCountsStore = (/0.,0.,0.,0./)
    prtCountsStoreAll = 0.
    prtCountsPrevStore = (/0.,0.,0.,0./)
    storedPrtCounts = (/0.,0.,0.,0./)
    storedPrtCountsAll = 0.    

    IF (PRESENT(doNotCopy)) THEN
       do_copy = .NOT. doNotCopy
    ELSE
       do_copy = .TRUE.
    END IF

    ! allow bad time ?
    IF( .not.PRESENT(bad_data) )THEN
       bad_data_in%allow_bad_top = .FALSE.
       bad_data_in%allow_bad_nav = .TRUE.
       bad_data_in%use_bad_nav = .FALSE.
       bad_data_in%allow_bad_calib = .FALSE.
       bad_data_in%allow_bad_time = .TRUE.
    ELSE
       bad_data_in = bad_data
    ENDIF
    IF( PRESENT(getTiePoints) &
         .OR. PRESENT(getLandMask) &
         .OR. PRESENT(getLocation) )THEN
       bad_data_in%allow_bad_top = .FALSE.
       bad_data_in%allow_bad_nav = .FALSE.
       bad_data_in%use_bad_nav = .FALSE.
       bad_data_in%allow_bad_calib = .TRUE.
       bad_data_in%allow_bad_time = .FALSE.
    ENDIF
    use_new_calibration_in = .FALSE.
    IF( PRESENT(use_new_calibration) )THEN
       use_new_calibration_in = use_new_calibration
    ENDIF
    calc_new_calibration_in = .FALSE.
    IF( PRESENT(calc_new_calibration) )THEN
       calc_new_calibration_in = calc_new_calibration
    ENDIF
    time_yearstart = 1975
    IF( PRESENT(use_2000_date) )THEN
       if( use_2000_date )then
          time_yearstart = 2000
       endif
    ENDIF
    copy_there = .FALSE.
    IF( PRESENT(copy_to_output) )THEN
       copy_there = copy_to_output
    ENDIF
    IF( PRESENT(outputData) )THEN
       IF( copy_there )THEN
          outData => outDataStore
          DataForFiltering => outDataStore
       ELSE
          outData => outputData
          DataForFiltering => outputData
       ENDIF
    ELSE
       outData => outDataStore
       DataForFiltering => outDataStore
    ENDIF
    out_radiances = .TRUE.
    IF( PRESENT(output_radiances) )THEN
       out_radiances = output_radiances
    ENDIF
    use_old_37_calibration_in = .FALSE.
    IF( PRESENT(use_old_37_calibration) )THEN
       use_old_37_calibration_in = use_old_37_calibration
    ENDIF
    load_landmask = .FALSE.
    IF( PRESENT(no_landmask) )THEN
       load_landmask = .not.no_landmask
    ENDIF
    IF( PRESENT(correct_solar_simple) )THEN
       corr_solar_simple = correct_solar_simple
    ELSE
       corr_solar_simple = .TRUE.
    ENDIF
    IF( PRESENT(new_vis_cal) )THEN
       new_vis_calibration = new_vis_cal
    ELSE
       new_vis_calibration = .TRUE.
    ENDIF
    IF( PRESENT(calibration_coef_overlap) )THEN
       cal_coef_overlap = calibration_coef_overlap
    ELSE
       cal_coef_overlap = 1.
    ENDIF
    IF( PRESENT(nedt_scale) )THEN
       nedt_scl = nedt_scale
    ELSE
       nedt_scl = 1.
    ENDIF
    IF( PRESENT(use_old_solar) )THEN
       use_old_solar_contam = use_old_solar
    ELSE
       use_old_solar_contam = .TRUE.
    ENDIF
    IF( PRESENT(use_iasi_calibration) )THEN
       use_iasi_cal = use_iasi_calibration
    ELSE
       use_iasi_cal = .FALSE.
    ENDIF
    IF( PRESENT(use_walton) )THEN
       use_walton_in = use_walton
    ELSE
       use_walton_in = .TRUE.
    ENDIF

    call DeAllocate_OutData_pntr( outData )
    call DeAllocate_OutData_pntr( DataForFiltering )
    AVHRR_Setup = .FALSE.

    ! Set Clavr mask to FALSE as default
    outData%Clavr_There = .FALSE.
    outData%walton_there = .FALSE.

    !Filename = TRIM(IMG%DataDir)//TRIM(IMG%DataFile)
    Filename=TRIM(IMG%DataFile)
    print*, "filename", Filename
#ifdef USE_GZFILE
    Unit = gzopen(Filename,'rb')
    IF(Unit.EQ.0)THEN
       CALL Check_IOS(1,'Opening File','Load_Imagery',&
            'NOAA_LoadAVHRRLevel1B.f90',.TRUE.,Filename)
    END IF
#else
    Unit = New_File_Unit()
    OPEN(Unit,FILE=Filename,STATUS='OLD',FORM='Unformatted',&
         ACCESS='STREAM',IOSTAT=Ios)
    CALL Check_IOS(Ios,'Opening File','Load_Imagery',&
         'NOAA_LoadAVHRRLevel1B.f90',.TRUE.,Filename)
#endif

    ! Read Headers
    CALL Check_Header( Unit, VersionNo, sizeOfRecords )
    SELECT CASE(VersionNo)
    CASE(0,1)
       write(*,'('' Reading in data in POD format'')')
       CALL Read_Data_HeaderV1(Unit,Coefs,DataType,sizeOfRecords,AVHRR_No,&
            have_12micron,headerNScans)
    CASE(2)
       write(*,'('' Reading in data in KLM version 2'')')
       CALL Read_Data_HeaderV2(Unit,Coefs,DataType,sizeOfRecords,AVHRR_No,&
            headerNScans)
    CASE(3,5)
       write(*,'('' Reading in data in KLM version 3+'')')
       CALL Read_Data_HeaderV3(Unit,Coefs,DataType,sizeOfRecords,AVHRR_No,&
           headerNScans)
    CASE DEFAULT
       call out_ERROR('Version Number not recognised',&
            'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END SELECT

    !
    ! Copy instrument coefficients from header if needed
    !
    IF( PRESENT(out_instr_coefs) )THEN
       out_instr_coefs%nuC = Coefs%nuC
       out_instr_coefs%Aval = Coefs%Aval
       out_instr_coefs%Bval = Coefs%Bval
       out_instr_coefs%F = Coefs%F
       out_instr_coefs%w = Coefs%w
       out_instr_coefs%nuc_numTemps = Coefs%nuC_numTemps
       out_instr_coefs%nuc_minT = Coefs%nuC_minT
       out_instr_coefs%nuc_maxT = Coefs%nuC_maxT
       out_instr_coefs%nuc_Array = Coefs%nuC_Array
       out_instr_coefs%Nspace = Coefs%Nspace
       out_instr_coefs%nonLinearCoefs3 = Coefs%nonLinearCoefs3
       out_instr_coefs%nonLinearCoefs4 = Coefs%nonLinearCoefs4
       out_instr_coefs%nonLinearCoefs5 = Coefs%nonLinearCoefs5
       out_instr_coefs%prttempCoefs1 = Coefs%prttempCoefs1
       out_instr_coefs%prttempCoefs2 = Coefs%prttempCoefs2
       out_instr_coefs%prttempCoefs3 = Coefs%prttempCoefs3
       out_instr_coefs%prttempCoefs4 = Coefs%prttempCoefs4
       out_instr_coefs%prttempCoefs5 = Coefs%prttempCoefs5
       out_instr_coefs%prttempCoefs6 = Coefs%prttempCoefs6
       out_instr_coefs%temperatureCoefs1 = Coefs%temperatureCoefs1
       out_instr_coefs%temperatureCoefs2 = Coefs%temperatureCoefs2
       out_instr_coefs%temperatureCoefs3 = Coefs%temperatureCoefs3
       out_instr_coefs%temperatureCoefs4 = Coefs%temperatureCoefs4
       out_instr_coefs%temperatureCoefs5 = Coefs%temperatureCoefs5
       out_instr_coefs%temperatureCoefs6 = Coefs%temperatureCoefs6
       out_instr_coefs%temperatureCoefs7 = Coefs%temperatureCoefs7
       out_instr_coefs%motorCurrentCoefs = Coefs%motorCurrentCoefs
       out_instr_coefs%number_Target_Temps = Coefs%number_Target_Temps
       out_instr_coefs%number_Scene_Temps = Coefs%number_Scene_Temps
       out_instr_coefs%Scene_Temp = Coefs%Scene_Temp
       out_instr_coefs%Target_Temp = Coefs%Target_Temp
       out_instr_coefs%Correction_Factor4 = Coefs%Correction_Factor4
       out_instr_coefs%Correction_Factor5 = Coefs%Correction_Factor5
     
    ENDIF

    ! use new calibration only for the ones that are valid (all except NOAA-06 and TIROS-N)
    IF ( use_new_calibration_in ) THEN
       SELECT CASE( AVHRR_No )
       CASE(1) ! TIROS-N
          use_new_calibration_in = .TRUE.
          CALL Gbcs_Warning(.TRUE.,&
               'Note new calibration coefficiens taken from NOAA-07 so AVHRR mismatch',&
               'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90')
       CASE(6) ! NOAA-6
          use_new_calibration_in = .TRUE.
          CALL Gbcs_Warning(.TRUE.,&
               'Note new calibration coefficiens taken from NOAA-07 so AVHRR mismatch',&
               'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90')
       CASE(7) ! NOAA-7
          use_new_calibration_in = .TRUE.
       CASE(8) ! NOAA-8
          use_new_calibration_in = .TRUE.
       CASE(9) ! NOAA-9
          use_new_calibration_in = .TRUE.
       CASE(10) ! NOAA-10
          use_new_calibration_in = .TRUE.
       CASE(11) ! NOAA-11
          use_new_calibration_in = .TRUE.
       CASE(12) ! NOAA-12
          use_new_calibration_in = .TRUE.
       CASE(13) ! NOAA-13
          use_new_calibration_in = .TRUE.
       CASE(14) ! NOAA-14
          use_new_calibration_in = .TRUE.
       CASE(15) ! NOAA-15
          use_new_calibration_in = .TRUE.
       CASE(16) ! NOAA-16
          use_new_calibration_in = .TRUE.
       CASE(17) ! NOAA-17
          use_new_calibration_in = .TRUE.
       CASE(18) ! NOAA-18
          use_new_calibration_in = .TRUE.
       CASE(19) ! NOAA-19
          use_new_calibration_in = .TRUE.
       CASE(-1) ! Metop-A == Metop-02
          use_new_calibration_in = .TRUE.
       CASE DEFAULT
          use_new_calibration_in = .FALSE.
       END SELECT
    END IF

    if( .not. use_new_calibration_in .and. .not. calc_new_calibration_in )then
       if( use_walton_in )then
          CALL out_ERROR('To activate Walton calibration nust also use_new_calibration',&
               'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
    endif


    ! the start/end times in the header are the frame times, without any
    ! clock corrections,  which may or may not be the same as the scanline
    ! lines.  But the difference is less than a second and this is acceptable
    ! I think for these gross filters.
    !
    ! check dates when things are valid, or mostly messed up
    CALL checkDate( start_time%year, start_time%month, start_time%day, &
         AVHRR_No, badDayS, validSolZaDecimalS, validClockDriftInfoS, &
         validClockCorrectionS, validSBBCS )
    CALL checkDate( end_time%year, end_time%month, end_time%day, AVHRR_No, &
         badDayE, validSolZaDecimalE, validClockDriftInfoE, &
         validClockCorrectionE, validSBBCE )

    badDay = badDayS.OR.badDayE
    validSolZaDecimal = validSolZaDecimalS.AND.validSolZaDecimalE
    ! for info only
    validClockDriftInfo = validClockDriftInfoS.AND.validClockDriftInfoE
    useL1BClock = validClockCorrectionS.AND.validClockCorrectionE
    validSBBC = validSBBCS.AND.validSBBCE

    IF(badDay)THEN
       ! nothing to be done, error message and stop
       CALL out_ERROR('bad day (see list of dates)',&
            'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

#ifdef USE_HDF5
    IF(.NOT.useL1BClock)THEN
       ! generally should use the time-corrected locations from the CLAVR-x
       ! runs but these may not be available.
       IF(PRESENT(clavrNavPathname))THEN
          IF(clavrNavPathname.NE.'')THEN
             CALL Read_CLAVR_Nav(clavrNavPathname, headerNScans)
          ELSE
             useL1BClock = .TRUE.
          END IF
       ELSE
          useL1BClock = .TRUE.
       END IF
    END IF
#else
    useL1BClock = .TRUE.
#endif

    outData%AVHRR_No = AVHRR_No
    CALL setupScanAngles( dataType, AVHRR_No )

    outData%satelliteAlt_There = .FALSE.

    IF( PRESENT(avhrr_rad_coefs) )THEN
       avhrr_rad_coefs%nuc = Coefs%nuc
       avhrr_rad_coefs%aVal = Coefs%aVal
       avhrr_rad_coefs%bVal = Coefs%bVal
    ENDIF

    ! Allocate for each scan line
    ALLOCATE(record(sizeOfRecords),STAT=Ios)
    CALL Check_IOS(Ios,'Allocating record','Load_Imagery',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    scanLinePos = 0
    ! Read in scan line data
    ReadLoop: DO scan=1,headerNScans
#ifdef USE_GZFILE
       Ios = gzread(Unit,record)
       IF( Ios .LT. 0 )THEN ! -1 error
          CALL Check_IOS(1,'Reading L1B records','Load_Imagery',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       END IF
       IF( 0 .EQ. Ios)THEN ! 0 EOF
          exit ReadLoop
       ENDIF
#else
       READ(Unit,IOSTAT=Ios)record
       IF( 0 .ne. Ios )THEN
          exit ReadLoop
       ENDIF
#endif
       scanLinePos = scanLinePos + 1
       SELECT CASE(VersionNo)
       CASE(0,1)
          CALL Parse_RecordV1(Coefs,DataType,scanLinePos,sizeOfRecords,&
               record,outData,bad_data_in,out_radiances,have_12micron,&
               getTiePoints=getTiePoints,getLandMask=getLandMask,&
               getLocation=getLocation,printAllErrors=printAllErrors,&
               ignore_solar_contamination=ignore_solar_contamination)
       CASE(2)
          CALL Parse_RecordV2(Coefs,DataType,scanLinePos,sizeOfRecords,&
               record,outData,bad_data_in,out_radiances,&
               getTiePoints=getTiePoints,getLandMask=getLandMask,&
               getLocation=getLocation,printAllErrors=printAllErrors,&
               ignore_solar_contamination=ignore_solar_contamination)
       CASE(3,5)
          CALL Parse_RecordV3(Coefs,DataType,scanLinePos,sizeOfRecords,&
               record,outData,bad_data_in,out_radiances,&
               getTiePoints=getTiePoints,getLandMask=getLandMask,&
               getLocation=getLocation,printAllErrors=printAllErrors,&
               ignore_solar_contamination=ignore_solar_contamination)
       CASE DEFAULT
          call out_ERROR('Version Number not recognised',&
               'Load_Imagery','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       END SELECT
!------Filtering outliers
       CALL filtering_sp_bb(scanLinePos,outData,DataForFiltering)
       CALL filtering_prt(scanLinePos,outData,DataForFiltering)

         
    END DO ReadLoop

#ifdef USE_GZFILE
    Ios = gzclose(Unit)
#else
    CLOSE(Unit)
#endif
    nScans = scanLinePos
    print*, "nscans", nScans
    ! for earlier AVHRRs there are a few (100 per year) cases with fewer lines in the file than described in the header.

    CALL Gbcs_Warning(nScans /= headerNScans,&
         "Scan lines in file (" //STR(nScans)// ") does not match header (" &
         //STR(headerNScans)// ")",ROUTINE_NAME,MODULE_NAME)

    ! Copy to correct sizing
    CALL Reallocate_Final_outData( outData, nScans )
    CALL Reallocate_Final_outData( DataForFiltering, nScans )

    ! Set Datatype
    IF( AVHRR_GAC .eq. DataType )THEN
       outData%isGAC = .TRUE.
    ELSE
       outData%isGAC = .FALSE.
    ENDIF
  
    IF( .NOT. ( PRESENT( getTiePoints ) &
         .OR. PRESENT( getLandMask ) &
         .OR. PRESENT( getLocation ) ) )THEN
       ! Recalibrate if needed
       if( use_new_calibration_in .or. calc_new_calibration_in )then

          CALL Recalibrate_AVHRR( outData, out_radiances,&
               use_night_only=use_night_only, &
               correct_solar_contamination=correct_solar_contamination, &
               correct_solar_contamination_all=correct_solar_contamination_all,&
               correct_earthshine=correct_earthshine, &
               moon_events=moon_events, &
               ignore_solar_contamination=ignore_solar_contamination,&
               avhrr_rad_coefs=avhrr_rad_coefs,&
               correct_solar_simple=corr_solar_simple,&
               new_vis_cal=new_vis_calibration,cal_overlap=cal_coef_overlap,&
               nedt_scl=nedt_scl,use_old_solar=use_old_solar_contam,&
               use_iasi_calibration=use_iasi_cal)
!!*******************************************
! Modif Marine 10 fevrier 2017
!!********************************************          
          print*, "coefs%nuc",coefs%nuc
          print*, "coefs%aval",coefs%Aval
          print*, "coefs%bval", coefs%bval
          print*, "coefs%nspace",coefs%nspace
          print*, "orbital temp",outdata%orbital_temperature
          call Get_Calib_Coefficients_CCI( outData%time(outdata%start_valid), outData%AVHRR_No,&
          coefs1, coefs2, coefs3, ncoefs, coefs_frac1, coefs_frac2, &
          cal_coef_overlap, use_iasi_cal )
         DO scanlinepos=1, Nscans
           outdata%rict_c3(scanlinepos)=NAN_R 
           outdata%rict_c4(scanlinepos)=NAN_R 
           outdata%rict_c5(scanlinepos)=NAN_R 
           if ((outData%smoothPrt(scanlinepos) .ne. NAN_R ) &
.and. (outdata%smoothPRT(scanlinepos) .gt. 0)) THEN     
             outdata%rict_c3(scanlinepos)=convertRadiance(outData%smoothPrt(scanlinepos),&
coefs%nuc(1),coefs%Aval(1),coefs%bval(1))
             outdata%rict_c4(scanlinepos)=convertRadiance(outData%smoothPrt(scanlinepos),&
coefs%nuc(2),coefs%Aval(2),coefs%bval(2))
             outdata%rict_c5(scanlinepos)=convertRadiance(outData%smoothPrt(scanlinepos),&
coefs%nuc(3),coefs%Aval(3),coefs%bval(3))
           end if 
         END DO
         print*, "NAN_R", NAN_R
         print*, "min,max,lon", minval(outdata%lon), maxval(outdata%lon)
         print*, "min,max rict", minval(outdata%rict_c3), maxval(outdata%rict_c3)
         print*, "min,max c3", minval(outdata%counts3), maxval(outdata%counts3)
         print*, "min,max sp3", minval(outdata%smoothsp3), maxval(outdata%smoothsp3)
         print*, "min,max bb3", minval(outdata%smoothbb3), maxval(outdata%smoothbb3)
         print*, "min,max rict", minval(outdata%rict_c4), maxval(outdata%rict_c4)
         print*, "min,max c4", minval(outdata%counts4), maxval(outdata%counts4)
         print*, "min,max sp4", minval(outdata%smoothsp4), maxval(outdata%smoothsp4)
         print*, "min,max bb4", minval(outdata%smoothbb4), maxval(outdata%smoothbb4)
         print*, "min,max ric5", minval(outdata%rict_c5), maxval(outdata%rict_c5)
         print*, "min,max c5", minval(outdata%counts5), maxval(outdata%counts5)
         print*, "min,max sp5", minval(outdata%smoothsp5), maxval(outdata%smoothsp5)
         print*, "min,max bb5", minval(outdata%smoothbb5), maxval(outdata%smoothbb5)
          call average_per_orbite(outData,NScans)
          print*,"ucs/icit", Outdata%ucs3, outdata%ucict3
          print*, Outdata%ucs4, outdata%ucict4
          print*, Outdata%ucs5, outdata%ucict5
!---------Deuxieme boucle siur les scanlines pour calculer les uncertitudes
          print*, "deuxieme boucle pr calculer sensibilites et uncertitudes"
          DO scanLinePos=1, Nscans
           CALL calcul_sensitivites(scanLinePos,outdata,coefs1,coefs2,coefs3)
    ! CALL calcul_urict(scanLinePos,outdata,coefs1,coefs2,coefs3)
           CALL radiance_uncertainties(scanLinePos,outData,coefs1, coefs2, coefs3)
           END DO
          print*, "prt", minval(outdata%smoothPrt), maxval(outdata%smoothPrt)
          print*, "sensitivities"     
          print*, minval(outdata%dre_over_dce3), maxval(outdata%dre_over_dce3)
          print*, minval(outdata%dre_over_dcs3), maxval(outdata%dre_over_dcs3)
          print*, minval(outdata%dre_over_dcict3),maxval(outdata%dre_over_dcict3)
          print*, minval(outdata%dre_over_db03),maxval(outdata%dre_over_db03)
          print*, minval(outdata%dre_over_da13),maxval(outdata%dre_over_da13)
          print*, minval(outdata%dre_over_da23),maxval(outdata%dre_over_da23)
          print*, minval(outdata%dre_over_db13),maxval(outdata%dre_over_db13)
          print*, minval(outdata%dre_over_dtinstr3),maxval(outdata%dre_over_dtinstr3)
          print*, minval(outdata%dre_over_drict3),maxval(outdata%dre_over_drict3)
          
          print*, "4"
          print*, minval(outdata%dre_over_dce4), maxval(outdata%dre_over_dce4)
          print*, minval(outdata%dre_over_dcs4), maxval(outdata%dre_over_dcs4)
          print*, minval(outdata%dre_over_dcict4),maxval(outdata%dre_over_dcict4)
          print*, minval(outdata%dre_over_db04),maxval(outdata%dre_over_db04)
          print*, minval(outdata%dre_over_da14),maxval(outdata%dre_over_da14)
          print*, minval(outdata%dre_over_da24),maxval(outdata%dre_over_da24)
          print*, minval(outdata%dre_over_db14),maxval(outdata%dre_over_db14)
          print*, minval(outdata%dre_over_dtinstr4),maxval(outdata%dre_over_dtinstr4)
          print*, minval(outdata%dre_over_drict4),maxval(outdata%dre_over_drict4)
          print*, "5" 
          print*, minval(outdata%dre_over_dce5), maxval(outdata%dre_over_dce5)
          print*, minval(outdata%dre_over_dcs5), maxval(outdata%dre_over_dcs5)
          print*, minval(outdata%dre_over_dcict5),maxval(outdata%dre_over_dcict5)
          print*, minval(outdata%dre_over_db05),maxval(outdata%dre_over_db05)
          print*, minval(outdata%dre_over_da15),maxval(outdata%dre_over_da15)
          print*, minval(outdata%dre_over_da25),maxval(outdata%dre_over_da25)
          print*, minval(outdata%dre_over_db15),maxval(outdata%dre_over_db15)
          print*, minval(outdata%dre_over_dtinstr5),maxval(outdata%dre_over_dtinstr5)
          print*, minval(outdata%dre_over_drict5),maxval(outdata%dre_over_drict5)
           
           
!*************************************
! end modif marine
!*************************************  
          !
          ! If want to apply Walton et al. uniformly to operational arrays
          ! Note - uses recalibration solar flags if available.
          !
          IF( use_walton_in )THEN
             CALL Calibration_Walton(outData,out_radiances,nedt_scl)
          ENDIF

          ! make sure that channel 3B (new calibration) hasn't overflowed
#ifdef USE_IEEE
          WHERE(.NOT.ieee_is_finite(outData%new_array3B)) outData%new_array3B=&
               NAN_R
#elif defined(USE_G95)
          WHERE(isNAN(outData%new_array3B)) outData%new_array3B=NAN_R
#endif
       endif
       ! make sure that channel 3B hasn't overflowed
#ifdef USE_IEEE
       WHERE(.NOT.ieee_is_finite(outData%array3B)) outData%array3B = NAN_R
#elif defined(USE_G95)
       WHERE(isNAN(outData%array3B)) outData%array3B = NAN_R
#endif

#ifdef USE_HDF5
       ! get the CLAVR-x cloud mask (only for GAC, LAC has mask = 7
       ! (== unprocessed))
       CALL Read_CLAVR_Cld(outData, clavrCldPathname, headerNScans)

       ! get the CLAVR-x cloud probability (only for GAC, LAC has prob = NAN_R)
       CALL Read_CLAVR_Prb(outData, clavrPrbPathname, headerNScans)
#endif

       IF (do_copy) THEN
       ! copy to IMG structure
          CALL Copy_To_IMG( Coefs, outData, IMG, use_new_calibration_in, &
               use_old_37_calibration_in, start_valid, stop_valid, &
               load_landmask)
       END IF
    endif
    print*, "Recalibration faite" 
 DO scanlinepos=1,Nscans
             if (outdata%time(scanlinepos) .eq.-1.e+30) then
                outdata%time(scanlinepos)=real_fillvalue
             end if
             if (outdata%scanlinenumber(scanlinepos) .eq.-1.e+30) then
                outdata%scanlinenumber(scanlinepos)=real_fillvalue
             end if
 if (outdata%rict_c3(scanlinepos) .eq.-1.e+30) then
                outdata%rict_c3(scanlinepos)=real_fillvalue
             end if
 if (outdata%rict_c4(scanlinepos) .eq.-1.e+30) then
                outdata%rict_c4(scanlinepos)=real_fillvalue
             end if
              if (outdata%rict_c3(scanlinepos) .eq.-1.e+30) then
                outdata%rict_c4(scanlinepos)=real_fillvalue
             end if
             DO j=1,409
                if (outdata%lat(j,scanlinepos) .eq. -1.e+30)then
                  outdata%lat(j,scanlinepos)=integer_fillvalue
                end if 
                if (outdata%lon(j,scanlinepos) .eq. -1.e+30) then
                 outdata%lon(j,scanlinepos)=integer_fillvalue
                end if
                if ((outdata%lat(j,scanlinepos) .lt. -90) .or. &
                   (outdata%lon(j,scanlinepos) .gt. 90)) then
                 outdata%lat(j,scanlinepos)=integer_fillvalue
                end if
                if ((outdata%lon(j,scanlinepos) .lt. 0) .or. &
                   (outdata%lon(j,scanlinepos) .gt. 360)) then
                 outdata%lon(j,scanlinepos)=integer_fillvalue
                end if
             
                if (outdata%solza(j,scanlinepos) .eq. -1.e+30) then
                  outdata%solza(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%satza(j,scanlinepos) .eq. -1.e+30) then
                  outdata%satza(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%relaz(j,scanlinepos) .eq. -1.e+30) then
                  outdata%relaz(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%array1(j,scanlinepos) .eq. -1.e+30) then
                  outdata%array1(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%array2(j,scanlinepos) .eq. -1.e+30) then
                  outdata%array2(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%array3A(j,scanlinepos) .eq. -1.e+30) then
                  outdata%array3B(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%counts1(j,scanlinepos) .eq. -1.e+30) then
                  outdata%counts1(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%counts2(j,scanlinepos) .eq. -1.e+30) then
                  outdata%counts2(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%counts3(j,scanlinepos) .eq. -1.e+30) then
                  outdata%counts3(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%counts4(j,scanlinepos) .eq. -1.e+30) then
                  outdata%counts4(j,scanlinepos)=real_fillvalue
                end if
                if (outdata%counts5(j,scanlinepos) .eq. -1.e+30) then
                  outdata%counts5(j,scanlinepos)=real_fillvalue
                end if
             
              END DO
END DO

    ! Copy to output if needed
    IF( copy_there )THEN
       IF( PRESENT(outputData) )THEN
          write(*,'('' Copying to outputData'')')
          CALL Copy_OutData(outData,outputData,start_valid,stop_valid)
       ENDIF
  Filename=TRIM(IMG%DataFile)
       ! Deallocate outdata
       call DeAllocate_OutData_pntr( outData )
    ENDIF

    IF( .NOT. PRESENT(outputData) )THEN
       ! Deallocate outdata if it is local to this module in outDataStore
       CALL DeAllocate_OutData_pntr( outData )
       ! this doesn't appear to free up the memory according to VmSize
       ! in /proc/self/status.  At least not instantly: VmSize isn't as
       ! large later in the program---perhaps the VM isn't returned
       ! immediately?
    ENDIF

    ! deallocate CLAVR-x lon/lat if used
    IF(ALLOCATED(clavrLonAnchor))DEALLOCATE(clavrLonAnchor)
    IF(ALLOCATED(clavrLatAnchor))DEALLOCATE(clavrLatAnchor)

    write(*,'('' Read in AVHRR data '')')

  END SUBROUTINE Load_Imagery

#ifdef USE_GZFILE
  SUBROUTINE Get_Land_Mask( IMG, outputData, landMaskDataArg )

    INTEGER,PARAMETER::nlons=43200+2
    INTEGER,PARAMETER::nlats=21600+2
    INTEGER,PARAMETER::ilon0=0-1
    INTEGER,PARAMETER::ilon1=43200-1+1
    INTEGER,PARAMETER::ilat0=-10800-1
    INTEGER,PARAMETER::ilat1=10800-1+1

    TYPE(Imagery), INTENT(IN) :: IMG
    TYPE(AVHRR_Data), target :: outputData
    TYPE(AVHRR_Land_Mask), target :: landMaskDataArg

    INTEGER(kind=GbcsInt1),DIMENSION(NDATA_GAC),TARGET::toleranceGAC &
         =(/24,23,22,22,21,20,20,19,19,18,18,18,17,17,16,16,16,15,15,15, &
         15,14,14,14,14,13,13,13,13,13,12,12,12,12,12,12,11,11,11,11, &
         11,11,11,10,10,10,10,10,10,10,10,10, 9, 9, 9, 9, 9, 9, 9, 9, &
         9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10, &
         10,10,10,10,10,10,10,11,11,11,11,11,11,11,12,12,12,12,12,12, &
         13,13,13,13,13,14,14,14,14,15,15,15,16,16,16,17,17,17,18,18, &
         19,19,20,20,21,21,22,22,23/)
    INTEGER(kind=GbcsInt1),DIMENSION(NDATA_LAC),TARGET::toleranceLAC &
         =(/16,16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,15,15, &
         15,15,14,14,14,14,14,14,14,14,14,14,14,14,14,13,13,13,13,13, &
         13,13,13,13,13,13,13,13,13,13,13,12,12,12,12,12,12,12,12,12, &
         12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,11,11,11, &
         11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,10,10,10,10,10, &
         10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
         10,10,10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
         9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
         9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
/         6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
         9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
         9, 9, 9, 9, 9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
         10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11, &
         11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12, &
         12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13, &
         13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14, &
         14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,16,16, &
         16,16,16,16,16,16,16,16/)

    INTEGER(kind=GbcsInt1),DIMENSION(:),POINTER::tolerance => NULL()

    TYPE(AVHRR_Data), POINTER :: outData => NULL()
    TYPE(AVHRR_Land_Mask), POINTER :: landMaskData => NULL()

    integer :: nelems, nlines, line, scanLine

    INTEGER(kind=GbcsInt1),DIMENSION(:,:),ALLOCATABLE::distance
    REAL(kind=GbcsDble),DIMENSION(ilon0:ilon1)::lon
    REAL(kind=GbcsDble),DIMENSION(ilat0:ilat1)::lat
    REAL(kind=GbcsDble)::dx,dy,dx2,dy2,mx,my
    INTEGER::i
    INTEGER::e,l

#ifdef USE_GZFILE
    INTEGER(kind=c_int) :: Unit
#else
    INTEGER :: Unit
#endif
    INTEGER :: Ios

    outData => outputData
    landMaskData => landMaskDataArg

    IF ( outData%nelem .EQ. NDATA_GAC ) THEN
       tolerance => toleranceGAC
    ELSE IF ( outData%nelem .EQ. NDATA_LAC ) THEN
       tolerance => toleranceLAC
    ELSE
       CALL out_ERROR('Data not GAC or LAC',&
            'Get_Land_Mask','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF

    ! set up the orbit land mask
    nelems = outData%nelem
    nlines = outData%arraySize

    landMaskData%arraySize = nlines
    landMaskData%nelem = nelems
    ALLOCATE(landMaskData%mask(nelems,nlines),stat=Ios)
    CALL Check_IOS(Ios,'Allocating land mask','Get_Land_Mask',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    ! set up the lons and lats of the global mask
    mx=120.0d0
    my=mx
    dx=1.0d0/mx
    dy=dx
    dx2=dx/2.0d0
    dy2=dx2
    lon=0.0d0+(/(i-1,i=1,nlons)/)*dx-0.5*dx
    lat=-90.0d0+(/(i-1,i=1,nlats)/)*dx-0.5*dx

    ! read the global land mask (distance to land) file
    ALLOCATE(distance(ilon0:ilon1,ilat0:ilat1))
    CALL Check_IOS(Ios,'Allocating distance','Get_Land_Mask',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    Unit = gzopen(IMG%LandSeaMaskFile,'rb')
    IF(Unit.EQ.0)THEN
       CALL Check_IOS(1,'Opening File','Get_Land_Mask',&
            'NOAA_LoadAVHRRLevel1B.f90',.TRUE.,IMG%LandSeaMaskFile)
    END IF
    Ios = gzread(Unit,distance)
    IF(Ios.LE.0)THEN
       CALL Check_IOS(1,'Reading land mask','Get_Land_Mask',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#ifdef USE_GZFILE
    Ios = gzclose(Unit)
#else
    CLOSE(Unit)
#endif

    ! find the indices of the orbit pixels in the global mask.
    ! Should work with rounding errors (lons 0.0-e, 360.0+e,
    ! lats -90.0-e,9 0.0+e) but still need to remove any bad values first.
    ! Compare distance and tolerance to derive an all-sea=Gbcs_Surface_Sea,
    ! some-land=Gbcs_Surface_Land mask.
    ! can't do in a oner:     ilon=FLOOR(mx*outData%Lon)
    !                         ilat=FLOOR(my*outData%Lat)
    !                         d(e,l)=distance(ilon(e,l),ilat(e,l)), etc

    landMaskData%mask=Gbcs_Surface_Space
    DO l=1,nlines
       DO e=1,nelems
          IF(outData%Lon(e,l).GT.-dx2 &
               .AND.outData%Lon(e,l).LT.360.0d0-dx2 &
               .AND.outData%lat(e,l).GT.-90.0d0-dy2 &
               .AND.outData%lat(e,l).LT.90.0d0+dy2)THEN
             ! better than          IF(outData%Lon(e,l).NE.NAN_R.AND.outData%lat(e,l).NE.NAN_R)THEN
             IF(tolerance(e).LE.distance(FLOOR(mx*outData%Lon(e,l)),FLOOR(my*outData%Lat(e,l))))THEN
                landMaskData%mask(e,l)=Gbcs_Surface_Sea
             ELSE
                landMaskData%mask(e,l)=Gbcs_Surface_Land
             END IF
          END IF
       END DO
    END DO

    IF(ALLOCATED(distance))DEALLOCATE(distance)

  end subroutine Get_Land_Mask
#endif

  subroutine Get_Tie_Points( outputData, tiePointDataArg )
  ! survey the image: is it ok to get tie points from?
  ! tie point start and step across track and along track
  ! put the tiepoints into a smaller AVHRR_Data structure

    TYPE(AVHRR_Data), target :: outputData, tiePointDataArg

    TYPE(AVHRR_Data), POINTER :: outData => NULL()
    TYPE(AVHRR_Data), POINTER :: tiePointData => NULL()

    INTEGER :: nelem, ntie, tie, line, line1, line2, scanLine, scanLineFinish
    REAL(GbcsDble) :: defaultTime
    logical :: goodLocation

    outData => outputData
    tiePointData => tiePointDataArg

    if ( outData%nelem .ne. NDATA_GAC ) then
       call out_ERROR('Data not GAC',&
            'Get_Tie_Points','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    endif

    nelem = ( NDATA_GAC - scanElemStart ) / scanElemStep + 1
    IF (outData%arraySize.LT.scanLineStart) THEN
       ntie = 0
       scanLineFinish = outData%arraySize
    ELSE
       ntie = ( outData%arraySize - scanLineStart ) / scanLineStep + 1
       scanLineFinish = scanLineStart + (ntie - 1) * scanLineStep
    END IF

    call Allocate_OutData( nelem, tiePointData )
    do while ( ntie .gt. tiePointData%arraySize )
       call Reallocate_outData( tiePointData, MEMORY_ALLOC_STEP )
    enddo
    CALL Reallocate_Final_outData( tiePointData, ntie )

    tiePointData%badTime = outData%badTime( scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%badNavigation = outData%badNavigation( scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%Lon = outData%Lon( scanElemStart::scanElemStep, scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%Lat = outData%Lat( scanElemStart::scanElemStep, scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%time = outData%time( scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%satZA = outData%satZa( scanElemStart::scanElemStep, scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%solZA = outData%solZA( scanElemStart::scanElemStep, scanLineStart:scanLineFinish:scanLineStep )
    tiePointData%relAz = outData%relAz( scanElemStart::scanElemStep, scanLineStart:scanLineFinish:scanLineStep )

    ! bad tie points get 0.0 for all the navigation values
    ! lon = 0, lat = 0 is a reasonable location, but
    ! something sensible for the time to make sure nwpgen.py
    ! is OK, choose the first valid.
    defaultTime = 0.0d0
    DO tie = 1, ntie
       IF (.NOT.tiePointData%badTime(tie)) THEN
          defaultTime = tiePointData%time(tie)
          EXIT
       END IF
    END DO
    ! if default time == 0.0 then nwpgen.py will probably fail

    ! dummy data in the bad tie points
    DO tie = 1, ntie
       IF (tiePointData%badTime(tie)) THEN
          tiePointData%time(tie) = defaultTime
       END IF
       IF (tiePointData%badNavigation(tie)) THEN
          tiePointData%Lon(:,tie) = 0.0
          tiePointData%Lat(:,tie) = 0.0
          tiePointData%satZA(:,tie) = 0.0
          tiePointData%solZA(:,tie) = 0.0
          tiePointData%relAz(:,tie) = 0.0
       END IF
    END DO

    ! make scan lines bad unless there are two good tie points
    ! surrounding them and there are no bad navigation points
    ! in between
    ! one tie point => all bad
    IF (ntie.LT.2)THEN
       ! all bad
       outData%badNavigation = .TRUE.
    ELSE
       ! lines before first tie point are bad
       DO line = 1, scanLineStart-1
          outData%badNavigation(line) =.TRUE.
       END DO
       ! lines after last tie point are bad (can have no lines)
       DO line = scanLineFinish+1, outData%arraySize
          outData%badNavigation(line) =.TRUE.
       END DO
       DO tie = 1, ntie-1
          line1 = scanLineStart + (tie - 1) * scanLineStep + 1
          line2 = scanLineStart + (tie+1 - 1) * scanLineStep - 1
          IF (tiePointData%badNavigation(tie) &
               .OR.tiePointData%badNavigation(tie+1) &
               .OR.ANY(outData%badNavigation(line1:line2))) THEN
             ! don't include the tie point lines
             outData%badNavigation(line1:line2) = .TRUE.
          END IF
       END DO
    END IF

    ! remove any lines (at tie points) that are good surrounded by bad
    ! this is done after testing each tie point section (see above)
    IF (ntie.LT.2)THEN
       ! all bad
       outData%badNavigation = .TRUE.
    ELSE
       ! first tie point
       line = scanLineStart
       IF (outData%badNavigation(line+1)) THEN
          outData%badNavigation(line) = .TRUE.
       ENDIF
       ! last tie point
       line = scanLineFinish
       IF (outData%badNavigation(line-1)) THEN
          outData%badNavigation(line) = .TRUE.
       ENDIF
       ! other tie points
       DO tie = 2, ntie-1
          line = scanLineStart + (tie - 1) * scanLineStep
          IF (outData%badNavigation(line-1) &
               .AND.outData%badNavigation(line+1)) THEN
             outData%badNavigation(line) = .TRUE.
          END IF
       END DO
    END IF

  end subroutine Get_Tie_Points

  ! Uses code from svdfit from Numerical Recipes (1992)
  INTEGER FUNCTION Poly_Fit( ndata, X, Y, ma, coefs )

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: X(ndata)
    REAL, INTENT(IN) :: Y(ndata)
    INTEGER, INTENT(IN) :: ma
    REAL, INTENT(OUT) :: coefs(ma)

    ! Local variables
    REAL, ALLOCATABLE :: covar(:,:)
    REAL, ALLOCATABLE :: afunc(:)
    REAL, ALLOCATABLE :: beta(:)
    INTEGER :: i,j,k,l,m
    REAL :: wt,ym
    INTEGER :: STAT

    Poly_Fit = 0
    ALLOCATE(covar(ma,ma),beta(ma),afunc(ma),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       call out_ERROR('Cannot allocate covar,beta,afunc',&
            'Poly_Fit','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    do j=1,ma
       do k=1,ma
          covar(j,k)=0.
       enddo
       beta(j)=0.
    enddo

    do i=1,ndata
       call funcs(x(i),afunc,ma)
       ym=y(i)
       do l=1,ma
          wt=afunc(l)
          do m=1,l
             covar(l,m)=covar(l,m)+wt*afunc(m)
          enddo
          beta(l)=beta(l)+ym*wt
       enddo
    enddo
    do j=2,ma
       do k=1,j-1
          covar(k,j)=covar(j,k)
       enddo
    enddo
    stat = gaussj(covar,ma,beta)
    IF( 0 .ne. stat )THEN
       Poly_Fit = 1
       RETURN
    ENDIF

    coefs = beta

    ! Deallocate memory
    DEALLOCATE(covar,beta,afunc)

  END FUNCTION Poly_Fit

  INTEGER FUNCTION gaussj(a,n,b)

    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: a(n,n)
    REAL, INTENT(INOUT) :: b(n)

    ! Local variables
    INTEGER, PARAMETER :: NMAX=50
    INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
    REAL :: big,dum,pivinv
    gaussj = 0
    do j=1,n
       ipiv(j)=0
    enddo
    do i=1,n
       big = 0.
       do j=1,n
          if(ipiv(j).ne.1)then
             do k=1,n
                if (ipiv(k).eq.0) then
                   if (abs(a(j,k)).ge.big)then
                      big=abs(a(j,k))
                      irow=j
                      icol=k
                   endif
                endif
             enddo
          endif
       enddo
       ipiv(icol)=ipiv(icol)+1
       if (irow.ne.icol) then
          do l=1,n
             dum=a(irow,l)
             a(irow,l)=a(icol,l)
             a(icol,l)=dum
          enddo
          dum=b(irow)
          b(irow)=b(icol)
          b(icol)=dum
       endif
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol).eq.0.) then
          call out_WARNING('Singular matrix',&
               'gaussj','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          gaussj = 1
          RETURN
       endif
       pivinv=1./a(icol,icol)
       a(icol,icol)=1.
       do l=1,n
          a(icol,l)=a(icol,l)*pivinv
       enddo
       b(icol)=b(icol)*pivinv
       do ll=1,n
          if(ll.ne.icol)then
             dum=a(ll,icol)
             a(ll,icol)=0.
             do l=1,n
                a(ll,l)=a(ll,l)-a(icol,l)*dum
             enddo
             b(ll)=b(ll)-b(icol)*dum
          endif
       enddo
    enddo
    do l=n,1,-1
       if(indxr(l).ne.indxc(l))then
          do k=1,n
             dum=a(k,indxr(l))
             a(k,indxr(l))=a(k,indxc(l))
             a(k,indxc(l))=dum
          enddo
       endif
    enddo
    return

  END FUNCTION gaussj

  SUBROUTINE funcs( x, afunc, ma )

    REAL, INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ma
    REAL, INTENT(OUT) :: afunc(ma)

    ! Local variables
    INTEGER :: I

    afunc(1) = 1
    DO I=2,ma
       afunc(I) = afunc(I-1)*x
    END DO

  END SUBROUTINE funcs

  ! Routine to fit a polynomial (quadratic) to gain^2 ws ICT radiance
  ! removes outliers to fit if needed
  LOGICAL FUNCTION Poly_Fit_Outliers(ndata,X,Y,coefs,stdev,sigma)RESULT(ok)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: X(ndata)
    REAL, INTENT(IN) :: Y(ndata)
    REAL, INTENT(OUT) :: coefs(3)
    REAL, INTENT(OUT) :: stdev
    REAL, INTENT(IN), OPTIONAL :: sigma

    ! Local variables
    INTEGER :: I, J
    INTEGER :: STAT
    REAL :: mean, sqval
    REAL :: standard_dev
    REAL, ALLOCATABLE :: result_val(:), newX(:), newY(:)

    ok = .TRUE.
    IF( 10 .ge. ndata )THEN
       ok = .FALSE.
       RETURN
    ENDIF

    stat = Poly_Fit(ndata,X,Y,3,coefs)

    ! Get standard deviation
    IF( ALLOCATED(result_val) )THEN
       DEALLOCATE(result_val)
    ENDIF
    ALLOCATE(result_val(ndata),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       call out_ERROR('Cannot allocate result_val',&
            'Poly_Fit_Outliers','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    mean = 0.
    sqval = 0.
    DO I=1,ndata
       result_val(I) = coefs(1) + coefs(2)*X(I) + coefs(3)*X(I)*X(I) - Y(I)
       mean = mean + result_val(I)
       sqval = sqval + result_val(I)*result_val(I)
    END DO
    mean = mean/ndata
    stdev = sqrt((sqval-ndata*mean*mean)/(ndata-1))
    standard_dev = sigma * stdev

    ! If want to sigma clip
    IF( PRESENT(sigma) )THEN

       IF( ALLOCATED(newX) )THEN
          DEALLOCATE(newX)
       ENDIF
       IF( ALLOCATED(newY) )THEN
          DEALLOCATE(newY)
       ENDIF
       ALLOCATE(newX(ndata),newY(ndata),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          call out_ERROR('Cannot allocate newX,newY',&
               'Poly_Fit_Outliers','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF

       ! Make new arrays
       J = 0
       DO I=1,ndata
          IF( standard_dev .gt. ABS(Y(I)-result_val(I)) )THEN
             J=J+1
             newX(J) = X(I)
             newY(J) = Y(I)
          ENDIF
       END DO

       ! Get sigma clipped coeffs
       IF( 10 .lt. J )THEN
          ! Get new coefficients
          stat = Poly_Fit(J,newX,newY,3,coefs)

          ! Get new standard deviation
          mean = 0.
          sqval = 0.
          DO I=1,J
             result_val(I) = coefs(1) + coefs(2)*newX(I) + &
                  coefs(3)*newX(I)*newX(I) - newY(I)
             mean = mean + result_val(I)
             sqval = sqval + result_val(I)*result_val(I)
          END DO
          mean = mean/J
          stdev = sqrt((sqval-ndata*mean*mean)/(J-1))

       ENDIF

       ! Deallocate memory
       DEALLOCATE(newX,newY)

       IF( 10 .ge. J )THEN
          ok = .FALSE.
          RETURN
       ENDIF

    ENDIF

    DEALLOCATE(result_val)

    RETURN

  END FUNCTION Poly_Fit_Outliers

  !
  ! Get CLAVRX coefficients for the older AVHRRs
  !
  ! CLAVR-X Version Nov. 2015 from CSPP CLAVRX_2_0
  !
  SUBROUTINE Get_CLAVRX_Coefs( outData, start_pos,&
       ch1_g_1, ch1_switch, ch1_g_2,&
       ch2_g_1, ch2_switch, ch2_g_2,&
       ch3_g_1, ch3_switch, ch3_g_2)

    USE GbcsDateTime

    TYPE(AVHRR_Data), INTENT(IN) :: outData
    INTEGER, INTENT(IN) :: start_pos
    REAL, INTENT(OUT) :: ch1_g_1, ch1_switch, ch1_g_2
    REAL, INTENT(OUT) :: ch2_g_1, ch2_switch, ch2_g_2
    REAL, INTENT(OUT) :: ch3_g_1, ch3_switch, ch3_g_2
    
    ! Local
    INTEGER :: I
    INTEGER :: new_start_pos
    REAL :: year_launch
    REAL :: Year
    REAL :: ch1_g0_low, ch1_g1_low, ch1_g2_low
    REAL :: ch1_g0_high, ch1_g1_high, ch1_g2_high
    REAL :: ch2_g0_low, ch2_g1_low, ch2_g2_low
    REAL :: ch2_g0_high, ch2_g1_high, ch2_g2_high
    REAL :: ch3_g0_low, ch3_g1_low, ch3_g2_low
    REAL :: ch3_g0_high, ch3_g1_high, ch3_g2_high
    REAL :: time_since_launch

    WRITE(*,'(''  Coefficient Version Date : Nov. 2015 taken from CSPP'')')

    ! Check year
    new_start_pos = -1
    Loop: DO I=start_pos,outData%arraySize
       IF( 0 .lt. outData%Year(I) )THEN
          new_start_pos = I
          EXIT Loop
       ENDIF
    END DO Loop
    IF( -1 .eq. new_start_pos )THEN
       CALL Gbcs_Critical(.TRUE.,'No year available for vis chan calibration',&
            'Get_CLAVRX_Coefs','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF

    ! Get year and fraction of year
    IF( Is_Leap_Year(outData%Year(new_start_pos)) )THEN
       Year = outData%Year(new_start_pos) + &
            (1.*outData%dayNo(new_start_pos))/366. +&
            outData%hours(new_start_pos)/(366.*24.)
    ELSE
       Year = outData%Year(new_start_pos) + &
            (1.*outData%dayNo(new_start_pos))/365. +&
            outData%hours(new_start_pos)/(365.*24.)
    ENDIF

    SELECT CASE(outData%AVHRR_No)
       CASE(1)
          year_launch=1978.783
          ch1_g0_low = 0.105
          ch1_g1_low = 27.015
          ch1_g2_low = -12.876
          ch1_g0_high = 0.105
          ch1_g1_high = 27.015
          ch1_g2_high = -12.876
          ch2_g0_low = 0.121
          ch2_g1_low = 10.7095
          ch2_g2_low = -0.643
          ch2_g0_high = 0.121
          ch2_g1_high = 10.7095
          ch2_g2_high = -0.643
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500
       CASE(6)
          year_launch=1979.490
          ch1_g0_low = 0.145
          ch1_g1_low = -28.117
          ch1_g2_low = 9.103
          ch1_g0_high = 0.145
          ch1_g1_high = -28.117
          ch1_g2_high = 9.103
          ch2_g0_low = 0.179
          ch2_g1_low = -30.385
          ch2_g2_low = 9.412
          ch2_g0_high = 0.179
          ch2_g1_high = -30.385
          ch2_g2_high = 9.412
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500          
       CASE(7)
          year_launch=1981.4764
          ch1_g0_low = 0.114
          ch1_g1_low = 5.492
          ch1_g2_low = -0.585
          ch1_g0_high = 0.114
          ch1_g1_high = 5.492
          ch1_g2_high = -0.585
          ch2_g0_low = 0.126
          ch2_g1_low = 6.801
          ch2_g2_low = -0.761
          ch2_g0_high = 0.126
          ch2_g1_high = 6.801
          ch2_g2_high = -0.761
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500                    
       CASE(8)
          year_launch=1983.241
          ch1_g0_low = 0.118
          ch1_g1_low = 14.334
          ch1_g2_low = -3.061
          ch1_g0_high = 0.118
          ch1_g1_high = 14.334
          ch1_g2_high = -3.061
          ch2_g0_low = 0.140
          ch2_g1_low = 17.711
          ch2_g2_low = -4.086
          ch2_g0_high = 0.140
          ch2_g1_high = 17.711
          ch2_g2_high = -4.086
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500                                    
       CASE(9)
          year_launch=1984.9480
          ch1_g0_low = 0.108
          ch1_g1_low = 6.657
          ch1_g2_low = -0.082
          ch1_g0_high = 0.108
          ch1_g1_high = 6.657
          ch1_g2_high = -0.082
          ch2_g0_low = 0.120
          ch2_g1_low = 5.340
          ch2_g2_low = -0.473
          ch2_g0_high = 0.120
          ch2_g1_high = 5.340
          ch2_g2_high = -0.473
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500                
       CASE(10)
          year_launch=1986.712
          ch1_g0_low = 0.103
          ch1_g1_low = 12.183
          ch1_g2_low = -2.010
          ch1_g0_high = 0.103
          ch1_g1_high = 12.183
          ch1_g2_high = -2.010
          ch2_g0_low = 0.134
          ch2_g1_low = 3.679
          ch2_g2_low = -0.270
          ch2_g0_high = 0.134
          ch2_g1_high = 3.679
          ch2_g2_high = -0.270
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500
       CASE(11)
          year_launch=1988.7064
          ch1_g0_low = 0.112
          ch1_g1_low = -0.876
          ch1_g2_low = 0.270
          ch1_g0_high = 0.112
          ch1_g1_high = -0.876
          ch1_g2_high = 0.270
          ch2_g0_low = 0.118
          ch2_g1_low = -0.238
          ch2_g2_low = 0.225
          ch2_g0_high = 0.118
          ch2_g1_high = -0.238
          ch2_g2_high = 0.225
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 999
          ch2_switch = 999
          ch3_switch = 999          
       CASE(12)
          year_launch=1991.3669
          ch1_g0_low = 0.129
          ch1_g1_low = 1.593
          ch1_g2_low = -0.130
          ch1_g0_high = 0.129
          ch1_g1_high = 1.593
          ch1_g2_high = -0.130
          ch2_g0_low = 0.151
          ch2_g1_low = 2.958
          ch2_g2_low = -0.255
          ch2_g0_high = 0.151
          ch2_g1_high = 2.958
          ch2_g2_high = -0.255
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500
       CASE(14)
          year_launch=1994.9699
          ch1_g0_low = 0.118
          ch1_g1_low = 5.493
          ch1_g2_low = -0.549
          ch1_g0_high = 0.118
          ch1_g1_high = 5.493
          ch1_g2_high = -0.549
          ch2_g0_low = 0.150
          ch2_g1_low = 2.349
          ch2_g2_low = -0.044
          ch2_g0_high = 0.150
          ch2_g1_high = 2.349
          ch2_g2_high = -0.044
          ch3_g0_low = 0.             
          ch3_g1_low = 0.             
          ch3_g2_low = 0.             
          ch3_g0_high = 0.            
          ch3_g1_high = 0.            
          ch3_g2_high = 0.            
          ch1_switch = 999
          ch2_switch = 999
          ch3_switch = 999
       CASE(15)
          year_launch=1998.3641
          ch1_g0_low = 0.061
          ch1_g1_low = -0.633
          ch1_g2_low = 0.031
          ch1_g0_high = 0.184
          ch1_g1_high = -0.633
          ch1_g2_high = 0.031
          ch2_g0_low = 0.070
          ch2_g1_low = 0.248
          ch2_g2_low = -0.010
          ch2_g0_high = 0.209
          ch2_g1_high = 0.248
          ch2_g2_high = -0.010
          ch3_g0_low = 0.            ! Don't trust this...  
          ch3_g1_low = 0.            !                      
          ch3_g2_low = 0.            !                      
          ch3_g0_high = 0.           !                      
          ch3_g1_high = 0.           !                      
          ch3_g2_high = 0.           !                      
          ch1_switch = 500
          ch2_switch = 500
          ch3_switch = 500
       CASE(16)
          year_launch=2000.7228
          ch1_g0_low = 0.0560
          ch1_g1_low = 0.427
          ch1_g2_low = -0.003
          ch1_g0_high = 0.169
          ch1_g1_high = 0.417
          ch1_g2_high = -0.003
          ch2_g0_low = 0.0580
          ch2_g1_low = 2.392
          ch2_g2_low = -0.269
          ch2_g0_high = 0.175
          ch2_g1_high = 2.392
          ch2_g2_high = -0.269
          ch3_g0_low = 0.030
          ch3_g1_low = -12.115
          ch3_g2_low = 4.462
          ch3_g0_high = 0.211
          ch3_g1_high = -12.115
          ch3_g2_high = 4.462
          ch1_switch = 498.96
          ch2_switch = 500.17
          ch3_switch = 499.43
       CASE(17)
          year_launch=2002.47912
          ch1_g0_low = 0.059
          ch1_g1_low = 0.199
          ch1_g2_low = 0.093
          ch1_g0_high = 0.176
          ch1_g1_high = 0.199
          ch1_g2_high = 0.093
          ch2_g0_low = 0.066
          ch2_g1_low = 3.635
          ch2_g2_low = -0.322
          ch2_g0_high = 0.199
          ch2_g1_high = 3.635
          ch2_g2_high = -0.322
          ch3_g0_low = 0.031
          ch3_g1_low = 2.040
          ch3_g2_low = -0.123
          ch3_g0_high = 0.217
          ch3_g1_high = 2.040
          ch3_g2_high = -0.123
          ch1_switch = 501.12
          ch2_switch = 500.73
          ch3_switch = 501.37
       CASE(18)
          year_launch = 2005.18891
          ch1_g0_low = 0.056
          ch1_g1_low = 2.058
          ch1_g2_low = -0.121
          ch1_g0_high = 0.167
          ch1_g1_high = 2.058
          ch1_g2_high = -0.121
          ch2_g0_low = 0.062
          ch2_g1_low = 2.655
          ch2_g2_low = -0.155
          ch2_g0_high = 0.186
          ch2_g1_high = 2.655
          ch2_g2_high = -0.155
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 500.54
          ch2_switch = 500.40
          ch3_switch = 500.56
       CASE(19)
          year_launch=2009.096
          ch1_g0_low = 0.054
          ch1_g1_low = 2.247
          ch1_g2_low = -0.317
          ch1_g0_high = 0.162
          ch1_g1_high = 2.247
          ch1_g2_high = -0.317
          ch2_g0_low = 0.061
          ch2_g1_low = 1.693
          ch2_g2_low = -0.132
          ch2_g0_high = 0.182
          ch2_g1_high = 1.693
          ch2_g2_high = -0.132
          ch3_g0_low = 0.
          ch3_g1_low = 0.
          ch3_g2_low = 0.
          ch3_g0_high = 0.
          ch3_g1_high = 0.
          ch3_g2_high = 0.
          ch1_switch = 496.43
          ch2_switch = 500.37
          ch3_switch = 496.11
       CASE(-1)
          year_launch=2006.77
          ch1_g0_low = 0.056
          ch1_g1_low = 0.906
          ch1_g2_low = -0.024
          ch1_g0_high = 0.167
          ch1_g1_high = 0.906
          ch1_g2_high = -0.024
          ch2_g0_low = 0.066
          ch2_g1_low = 0.814
          ch2_g2_low = 0.025
          ch2_g0_high = 0.199
          ch2_g1_high = 0.814
          ch2_g2_high = 0.025
          ch3_g0_low = 0.031
          ch3_g1_low = 2.020
          ch3_g2_low = -0.115
          ch3_g0_high = 0.220
          ch3_g1_high = 2.020
          ch3_g2_high = -0.115
          ch1_switch = 501.00
          ch2_switch = 500.00
          ch3_switch = 502.00
       CASE(-2)
          year_launch=2012.77
          ch1_g0_low = 0.060
          ch1_g1_low = -7.442
          ch1_g2_low = 2.828
          ch1_g0_high = 0.179
          ch1_g1_high = -7.442
          ch1_g2_high = 2.828
          ch2_g0_low = 0.058
          ch2_g1_low = 10.222
          ch2_g2_low = -3.271
          ch2_g0_high = 0.174
          ch2_g1_high = 10.222
          ch2_g2_high = -3.271
          ch3_g0_low = 0.035
          ch3_g1_low = -27.970
          ch3_g2_low = 10.925
          ch3_g0_high = 0.247
          ch3_g1_high = -27.970
          ch3_g2_high = 10.925
          ch1_switch = 501.12
          ch2_switch = 500.82
          ch3_switch = 501.32
       CASE DEFAULT
          call out_ERROR('Accessing CLAVR-X coeffs beyond MetOp-B',&
               'Get_CLAVX_Coeffs','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END SELECT

    time_since_launch = Year - year_launch
    ch1_g_1 = ch1_g0_low*(100.0 + ch1_g1_low*time_since_launch + &
         ch1_g2_low*time_since_launch**2)/100.
    ch1_g_2 = ch1_g0_high*(100.0 + ch1_g1_high*time_since_launch + &
         ch1_g2_high*time_since_launch**2)/100.

    ch2_g_1 = ch2_g0_low*(100.0 + ch2_g1_low*time_since_launch + &
         ch2_g2_low*time_since_launch**2)/100.
    ch2_g_2 = ch2_g0_high*(100.0 + ch2_g1_high*time_since_launch + &
         ch2_g2_high*time_since_launch**2)/100.

    ch3_g_1 = ch3_g0_low*(100.0 + ch3_g1_low*time_since_launch + &
         ch3_g2_low*time_since_launch**2)/100.
    ch3_g_2 = ch3_g0_high*(100.0 + ch3_g1_high*time_since_launch + &
         ch3_g2_high*time_since_launch**2)/100.

  END SUBROUTINE Get_CLAVRX_Coefs

  !
  ! Use CLAVR-X coefficients for pre 2004 vischan data
  ! After 2004 assume operational calibration OK
  ! as it should have been updated by then
  !
  SUBROUTINE Update_VisChan( outData, start_pos, stop_pos, noise, nedt_scale,&
       noise_all_data )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: start_pos, stop_pos
    REAL, INTENT(IN) :: noise(:,:)
    REAL, INTENT(IN) :: nedt_scale
    LOGICAL, INTENT(IN) :: noise_all_data

    ! Local variabiles
    INTEGER :: I,J,K,KK
    INTEGER :: minx, maxx
    INTEGER :: nmean1, nmean2, nmean3
    REAL :: mean1, mean2, mean3
    REAL :: ch1_g_1, ch1_switch, ch1_g_2
    REAL :: ch2_g_1, ch2_switch, ch2_g_2
    REAL :: ch3_g_1, ch3_switch, ch3_g_2

!    ! Check if after 01-01-2004
!    IF( outData%year(start_pos) .ge. VISCHAN_OPER_YEAR_MIN )THEN
!       !
!       ! Copy operational calibration over as should be OK
!       !
!       outData%new_array1 = outData%array1
!       outData%new_array2 = outData%array2
!       outData%new_array3a = outData%array3a
!       WRITE(*,'(''...Using operational VISCHAN calibration...'')')
!       RETURN
!    ENDIF

    outData%new_array1 = NAN_R
    outData%new_array2 = NAN_R
    outData%new_array3a = NAN_R
    outData%new_array1_error = NAN_R
    outData%new_array2_error = NAN_R
    outData%new_array3a_error = NAN_R

    !
    ! For older sensors use the CLAVR-X degradation model
    !
    CALL Get_CLAVRX_Coefs( outData, start_pos,&
         ch1_g_1, ch1_switch, ch1_g_2,&
         ch2_g_1, ch2_switch, ch2_g_2,&
         ch3_g_1, ch3_switch, ch3_g_2)
    !
    ! Loop round data
    !
    DO I=start_pos, stop_pos
       if( noise_all_data )then
          KK=1
       ELSE
          KK=I
       ENDIF
       !
       ! Get smoothed space counts
       !
       minx = MAX(1,I-NPIXEL_PRT_SMOOTH)
       maxx = MIN(outData%arraySize,I+NPIXEL_PRT_SMOOTH)
       nmean1 = 0
       mean1 = 0.
       nmean2 = 0
       mean2 = 0.
       nmean3 = 0
       mean3 = 0.
       DO J=minx,maxx
          DO K=1,NBBODY_SAMPLES
             IF( 0 .lt. outData%spaceFilter1(K,J) )THEN
                nmean1 = nmean1 + 1
                mean1 = mean1 + outData%spaceFilter1(K,J)
             ENDIF
             IF( 0 .lt. outData%spaceFilter2(K,J) )THEN
                nmean2 = nmean2 + 1
                mean2 = mean2 + outData%spaceFilter2(K,J)
             ENDIF
             IF( 0 .lt. outData%spaceFilter3a(K,J) )THEN
                nmean3 = nmean3 + 1
                mean3 = mean3 + outData%spaceFilter3a(K,J)
             ENDIF
          END DO
       END DO

       ! Now get albedo
       if( nmean1 .gt. 0 )then
          mean1 = mean1 / nmean1
          DO J=1,outData%nelem
             IF( outData%Counts1(J,I) .le. ch1_switch .and. &
                  outData%Counts1(J,I) .ge. 0. )THEN
                IF( outData%Counts1(J,I) > mean1 )THEN
                   outData%new_array1(J,I) = &
                        ch1_g_1*(outData%Counts1(J,I)-mean1)
                   outData%new_array1_error(J,I) = ABS(ch1_g_1*noise(1,KK))
                ELSE
                   outData%new_array1(J,I) = 0.0
                   outData%new_array1_error(J,I) = 0.0
                ENDIF
             ELSE IF(outData%Counts1(J,I) .gt. ch1_switch)THEN
                outData%new_array1(J,I) = ch1_g_1*(ch1_switch-mean1) + &
                     ch1_g_2*(outData%Counts1(J,I)-ch1_switch)
                outData%new_array1_error(J,I) = ABS(ch1_g_2*noise(1,KK))
             ELSE
                outData%new_array1(J,I) = NAN_R
                outData%new_array1_error(J,I) = NAN_R
             ENDIF
             IF( outData%new_array1(J,I) .ge. 0 )THEN
                outData%new_array1(J,I) = outData%new_array1(J,I)/100.
                outData%new_array1_error(J,I) = &
                     outData%new_array1_error(J,I)/100.
             ENDIF
             IF( outData%isGAC )THEN
                IF( 0 .lt. outData%new_array1_error(J,I) )THEN
                   IF( noise(1,KK) .ge. 3 )THEN
                      outData%new_array1_error(J,I) = &
                           outData%new_array1_error(J,I)/2.
                   ELSE
                      outData%new_array1_error(J,I) = &
                           outData%new_array1_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             ENDIF
          END DO          
       endif
       if( nmean2 .gt. 0 )then
          mean2 = mean2 / nmean2
          DO J=1,outData%nelem
             IF( outData%Counts2(J,I) .le. ch2_switch .and. &
                  outData%Counts2(J,I) .ge. 0. )THEN
                IF( outData%Counts2(J,I) > mean2 )THEN
                   outData%new_array2(J,I) = &
                     ch2_g_1*(outData%Counts2(J,I)-mean2)
                   outData%new_array2_error(J,I) = ABS(ch2_g_1*noise(2,KK))
                ELSE
                   outData%new_array2(J,I) = 0.0
                   outData%new_array2_error(J,I) = 0.0
                ENDIF
             ELSE IF( outData%Counts2(J,I) .gt. ch2_switch )THEN
                outData%new_array2(J,I) = ch2_g_1*(ch2_switch-mean2) + &
                     ch2_g_2*(outData%Counts2(J,I)-ch2_switch)
                outData%new_array2_error(J,I) = ABS(ch2_g_2*noise(2,KK))
             ELSE
                outData%new_array2(J,I) = NAN_R
                outData%new_array2_error(J,I) = NAN_R
             ENDIF
             IF( outData%new_array2(J,I) .ge. 0 )THEN
                outData%new_array2(J,I) = outData%new_array2(J,I)/100.
                outData%new_array2_error(J,I) = &
                     outData%new_array2_error(J,I)/100.
             ENDIF
             IF( outData%isGAC )THEN
                IF( 0 .lt. outData%new_array2_error(J,I) )THEN
                   IF( noise(2,KK) .ge. 3 )THEN
                      outData%new_array2_error(J,I) = &
                           outData%new_array2_error(J,I)/2.
                   ELSE
                      outData%new_array2_error(J,I) = &
                           outData%new_array2_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             ENDIF
          END DO
       endif
       if( nmean3 .gt. 0 )then
          mean3 = mean3 / nmean3
          DO J=1,outData%nelem
             IF( outData%Counts3(J,I) .le. ch3_switch .and. &
                  outData%Counts3(J,I) .ge. 0. )THEN
                IF( outData%Counts3(J,I) > mean3 )THEN
                   outData%new_array3a(J,I) = &
                        ch3_g_1*(outData%Counts3(J,I)-mean3)
                   outData%new_array3a_error(J,I) = ABS(ch3_g_1*noise(3,KK))
                ELSE 
                   outData%new_array3a(J,I) = 0.0
                   outData%new_array3a_error(J,I) = 0.0
                ENDIF
             ELSE IF( outData%Counts3(J,I) .gt. ch3_switch )THEN
                outData%new_array3a(J,I) = ch3_g_1*(ch3_switch-mean3) + &
                     ch3_g_2*(outData%Counts3(J,I)-ch3_switch)
                outData%new_array3a_error(J,I) = ABS(ch2_g_2*noise(3,KK))
             ELSE
                outData%new_array3a(J,I) = NAN_R
                outData%new_array3a_error(J,I) = NAN_R
             ENDIF
             IF( outData%new_array3a(J,I) .ge. 0 )THEN
                outData%new_array3a(J,I) = outData%new_array3a(J,I)/100.
                outData%new_array3a_error(J,I) = &
                     outData%new_array3a_error(J,I)/100.
             ENDIF
             IF( outData%isGAC )THEN
                IF( 0 .lt. outData%new_array3a_error(J,I) )THEN
                   IF( noise(3,KK) .ge. 3 )THEN
                      outData%new_array3a_error(J,I) = &
                           outData%new_array2_error(J,I)/2.
                   ELSE
                      outData%new_array3a_error(J,I) = &
                           outData%new_array3a_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             ENDIF
          END DO
       endif

    END DO
    
    WRITE(*,'(''...Using CLAVRX VISCHAN calibration...'')')

  END SUBROUTINE Update_VisChan

  ! Get noise estimates
  SUBROUTINE Get_Allan_deviation( outData, start_pos, stop_pos, noise, all )

    TYPE(AVHRR_Data), pointer :: outData
    INTEGER, INTENT(IN) :: start_pos
    INTEGER, INTENT(IN) :: stop_pos
    REAL, ALLOCATABLE, INTENT(OUT) :: noise(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all

    ! Local variables
    INTEGER :: I,J,K
    INTEGER :: startpos,stoppos
    INTEGER :: STAT
    INTEGER :: minx, maxx
    INTEGER :: nallan3, nallan4, nallan5
    REAL :: allan3, allan4, allan5
    INTEGER :: nallan1, nallan2, nallan3a
    REAL :: allan1, allan2, allan3a
    LOGICAL :: all_data

    IF( PRESENT(all) )THEN
       alL_data=all
    ELSE
       all_data=.FALSE.
    ENDIF

    IF( all_data )THEN
       ALLOCATE(noise(6,1),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          call out_ERROR('Cannot allocate noise',&
               'Get_Allan_Deviation','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       startpos=1
       stoppos=1
    ELSE
       ALLOCATE(noise(6,outData%arraySize),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          call out_ERROR('Cannot allocate noise',&
               'Get_Allan_Deviation','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       startpos=start_pos
       stoppos=stop_pos
    ENDIF
    noise = -1

    outLoop: DO I=startpos,stoppos
       IF( all_data )THEN
          minx = start_pos
          maxx = stop_pos
       ELSE
          minx = MAX(1,I-NPIXEL_PRT_SMOOTH)
          maxx = MIN(outData%arraySize,I+NPIXEL_PRT_SMOOTH)
       ENDIF
       allan1 = 0.
       nallan1 = 0
       allan2 = 0.
       nallan2 = 0
       allan3a = 0.
       nallan3a = 0
       allan3 = 0.
       nallan3 = 0
       allan4 = 0.
       nallan4 = 0
       allan5 = 0.
       nallan5 = 0
       DO J=minx,maxx
          DO K=1,NBBODY_SAMPLES-1
             IF( outData%spaceFilter1(K,J) .gt. 0 .and. &
                  outData%spaceFilter1(K+1,J) .gt. 0 )THEN
                allan1 = allan1 + (outData%spaceFilter1(K+1,J) - &
                     outData%spaceFilter1(K,J))**2
                nallan1 = nallan1 + 1
             ENDIF
             IF( outData%spaceFilter2(K,J) .gt. 0 .and. &
                  outData%spaceFilter2(K+1,J) .gt. 0 )THEN
                allan2 = allan2 + (outData%spaceFilter2(K+1,J) - &
                     outData%spaceFilter2(K,J))**2
                nallan2 = nallan2 + 1
             ENDIF
             IF( outData%transition3A3B(J) )THEN
                IF( outData%spaceFilter3(K,J) .gt. 0 .and. &
                     outData%spaceFilter3(K+1,J) .gt. 0 )THEN
                   allan3 = allan3 + (outData%spaceFilter3(K+1,J) - &
                        outData%spaceFilter3(K,J))**2
                   nallan3 = nallan3 + 1
                ENDIF
             ELSE
                IF( outData%spaceFilter3a(K,J) .gt. 900 .and. &
                     outData%spaceFilter3a(K+1,J) .gt. 900 )THEN
                   allan3a = allan3a + (outData%spaceFilter3a(K+1,J) - &
                        outData%spaceFilter3a(K,J))**2
                   nallan3a = nallan3a + 1
                ENDIF
             ENDIF
             IF( outData%spaceFilter4(K,J) .gt. 900 .and. &
                  outData%spaceFilter4(K+1,J) .gt. 900 )THEN
                allan4 = allan4 + (outData%spaceFilter4(K+1,J) - &
                     outData%spaceFilter4(K,J))**2
                nallan4 = nallan4 + 1
             ENDIF
             IF( outData%spaceFilter5(K,J) .gt. 900 .and. &
                  outData%spaceFilter5(K+1,J) .gt. 900 )THEN
                allan5 = allan5 + (outData%spaceFilter5(K+1,J) - &
                     outData%spaceFilter5(K,J))**2
                nallan5 = nallan5 + 1
             ENDIF
          END DO
       END DO
       IF( nallan1 .gt. 10 )THEN
          noise(1,I) = SQRT(allan1/(2*(nallan1-1)))
          IF( 0 .eq. noise(1,I) )THEN
             Loop1: DO J=I,start_pos,-1
                IF( noise(1,J) .gt. 0 )THEN
                   noise(1,I) = noise(1,J)
                   EXIT Loop1
                ENDIF
             END DO Loop1
          ENDIF
       ENDIF
       IF( nallan2 .gt. 10 )THEN
          noise(2,I) = SQRT(allan2/(2*(nallan2-1)))
          IF( 0 .eq. noise(2,I) )THEN
             Loop2: DO J=I,start_pos,-1
                IF( noise(2,J) .gt. 0 )THEN
                   noise(2,I) = noise(2,J)
                   EXIT Loop2
                ENDIF
             END DO Loop2
          ENDIF
       ENDIF
       IF( nallan3a .gt. 10 )THEN
          noise(3,I) = SQRT(allan3a/(2*(nallan3a-1)))
          IF( 0 .eq. noise(3,I) )THEN
             Loop3: DO J=I,start_pos,-1
                IF( noise(3,J) .gt. 0 )THEN
                   noise(3,I) = noise(3,J)
                   EXIT Loop3
                ENDIF
             END DO Loop3
          ENDIF
       ENDIF
       IF( nallan3 .gt. 10 )THEN
          noise(4,I) = SQRT(allan3/(2*(nallan3-1)))
          IF( 0 .eq. noise(4,I) )THEN
             Loop4: DO J=I,start_pos,-1
                IF( noise(4,J) .gt. 0 )THEN
                   noise(4,I) = noise(4,J)
                   EXIT Loop4
                ENDIF
             END DO Loop4
          ENDIF
       ENDIF
       IF( nallan4 .gt. 10 )THEN
          noise(5,I) = SQRT(allan4/(2*(nallan4-1)))
          IF( 0 .eq. noise(5,I) )THEN
             Loop5: DO J=I,start_pos,-1
                IF( noise(5,J) .gt. 0 )THEN
                   noise(5,I) = noise(5,J)
                   EXIT Loop5
                ENDIF
             END DO Loop5
          ENDIF
       ENDIF
       IF( nallan5 .gt. 10 )THEN
          noise(6,I) = SQRT(allan5/(2*(nallan5-1)))
          IF( 0 .eq. noise(6,I) )THEN
             Loop6: DO J=I,start_pos,-1
                IF( noise(6,J) .gt. 0 )THEN
                   noise(6,I) = noise(6,J)
                   EXIT Loop6
                ENDIF
             END DO Loop6
          ENDIF
       ENDIF
    END DO outLoop

  END SUBROUTINE GET_ALLAN_DEVIATION

  ! Do the AVHRR re-calibration.  This needs to be updated to unset
  ! outData%badCalibration for lines that had bad old calibration
  ! but good new calibration.
  SUBROUTINE Recalibrate_AVHRR( outData, out_radiances, use_night_only, &
       correct_solar_contamination, correct_earthshine, moon_events, &
       ignore_solar_contamination, correct_solar_contamination_all,&
       avhrr_rad_coefs,cal_overlap,correct_solar_simple,new_vis_cal,&
       nedt_scl,use_old_solar,noise_orbit,use_iasi_calibration)

    TYPE(AVHRR_Data), POINTER :: outData
    LOGICAL, INTENT(IN) :: out_radiances
    LOGICAL, INTENT(IN), OPTIONAL :: use_night_only
    LOGICAL, INTENT(IN), OPTIONAL :: correct_solar_contamination
    LOGICAL, INTENT(IN), OPTIONAL :: correct_earthshine
    LOGICAL, INTENT(IN), OPTIONAL :: moon_events
    LOGICAL, INTENT(IN), OPTIONAL :: ignore_solar_contamination
    LOGICAL, INTENT(IN), OPTIONAL :: correct_solar_contamination_all
    TYPE(AVHRR_Radiative_Coefs), OPTIONAL, INTENT(INOUT) :: avhrr_rad_coefs
    REAL, INTENT(IN), OPTIONAL :: cal_overlap
    LOGICAL, INTENT(IN), OPTIONAL :: correct_solar_simple
    LOGICAL, INTENT(IN), OPTIONAL :: new_vis_cal
    REAL, INTENT(IN), OPTIONAL :: nedt_scl
    LOGICAL, INTENT(IN), OPTIONAL :: use_old_solar
    LOGICAL, INTENT(IN), OPTIONAL :: noise_orbit
    LOGICAL, INTENT(IN), OPTIONAL :: use_iasi_calibration

    INTEGER :: nlines_removed
    INTEGER :: I, J, K, KK, JJ, K2, LL
    INTEGER :: J3, J4, J5

    INTEGER :: pos,pos_2,modpoint
    INTEGER :: start_pos
    INTEGER :: stop_pos
    LOGICAL :: interp_daynight
    LOGICAL :: correct_solar_contam
    LOGICAL :: correct_solar_contam_all
    LOGICAL :: correct_earthshine_data
    LOGICAL :: get_moon_events
    LOGICAL :: ignore_solar
    LOGICAL :: noise_all_data

    INTEGER :: minx, maxx
    REAL :: newa, newb
    INTEGER :: ntotal
    INTEGER :: ntotal_cnts
    INTEGER :: ntotal0
    REAL :: total1, total2, total3, total4
    REAL :: total1_cnts, total2_cnts, total3_cnts, total4_cnts
    REAL :: gain
    REAL :: R_ICT
    REAL :: counts
    REAL :: coefs1(8,2)
    REAL :: coefs2(8,2)
    REAL :: coefs3(8,2)
    REAL :: radiance

    INTEGER :: new_ndata
    REAL :: poly_coefs3(3), poly_coefs4(3), poly_coefs5(3)
    REAL :: stdev_val3, stdev_val4, stdev_val5
    REAL :: maxdev3, maxdev4, maxdev5
    REAL :: sigma_thresh
    LOGICAL :: use_poly_coefs3, use_poly_coefs4, use_poly_coefs5
    REAL :: stdev3, stdev4, stdev5

    INTEGER :: STAT

    REAL, ALLOCATABLE :: gain_array3(:,:), rad_array3(:)
    REAL, ALLOCATABLE :: gain_array4(:,:), rad_array4(:)
    REAL, ALLOCATABLE :: gain_array5(:,:), rad_array5(:)
    REAL, ALLOCATABLE :: rad_array3_correl(:), rad_array4_correl(:)
    REAL, ALLOCATABLE :: rad_array5_correl(:)
    LOGICAL, ALLOCATABLE :: good_data(:), solar_data(:)
    LOGICAL :: night_ok
    REAL :: radiance_error, error1, error2
    INTEGER :: mean3, mean4, mean5
    INTEGER :: sign3, sign4, sign5
    INTEGER :: solarStat
    INTEGER :: solarStat3,solarStat4,solarStat5
    REAL :: correlPRT
    LOGICAL :: printSolarWarn
    INTEGER :: ncoefs
    REAL, ALLOCATABLE :: new_calib3(:,:,:)
    REAL, ALLOCATABLE :: new_calib4(:,:,:)
    REAL, ALLOCATABLE :: new_calib5(:,:,:)
    REAL :: radiance2
    REAL :: coefs_frac1, coefs_frac2
    INTEGER :: norbital
    LOGICAL :: corr_solar_simple
    INTEGER :: n
    REAL :: mean
    REAL :: delta
    INTEGER :: n2
    REAL :: M2
    REAL :: mean_2
    REAL :: M2_2
    REAL :: gain_val,gain_val_2
    LOGICAL :: ok,ok2
    LOGICAL :: reset_gain
    REAL :: cal_coef_overlap
    LOGICAL :: operational_solar
    LOGICAL :: found
    REAL :: mean_gain
    LOGICAL, ALLOCATABLE :: solza_ascend(:)
    INTEGER :: flat_start, flat_stop
    INTEGER :: minI1, maxI1
    INTEGER :: minI2, maxI2
    INTEGER :: nmean1, nmean2
    REAL :: mean1, mean2
    LOGICAL :: new_vis_calibration
    REAL, ALLOCATABLE :: noise(:,:)
    REAL :: nedt_scale
    REAL :: max_stdev,max_stdev_2
    REAL, ALLOCATABLE :: smooth_solza(:)
    LOGICAL :: main_event
    REAL :: minza,maxza
    INTEGER :: midpoint
    INTEGER, ALLOCATABLE :: Used(:)
    INTEGER :: nUsed(100), UsedPos(100)
    LOGICAL :: use_old_solar_contam
    LOGICAL :: use_iasi_cal

    interp_daynight = .FALSE.
    IF( PRESENT(use_night_only) )THEN
       interp_daynight = use_night_only
       IF( interp_daynight )THEN
          write(*,'(''    Interpolating over daytime'')')
       ENDIF
    ENDIF
    correct_solar_contam = .FALSE.
    IF( PRESENT(correct_solar_contamination) )THEN
       correct_solar_contam = correct_solar_contamination
       IF( correct_solar_contam )THEN
          write(*,'(''    Correcting for solar contamination if present'')')
       ENDIF
    ENDIF
    correct_solar_contam_all = .FALSE.
    IF( PRESENT(correct_solar_contamination_all) )THEN
       correct_solar_contam_all = correct_solar_contamination_all
       IF( correct_solar_contam )THEN
          write(*,&
       '(''    Correcting for solar contamination (All Chans) if present'')')
       ENDIF
    ENDIF
    correct_earthshine_data = .FALSE.
    IF( PRESENT(correct_earthshine) )THEN
       correct_earthshine_data = correct_earthshine
       IF( correct_earthshine_data )THEN
          write(*,'(''    Correcting for Earthshine in 3.7mu channel'')')
       ENDIF
    ENDIF
    get_moon_events = .FALSE.
    IF( PRESENT(moon_events) )THEN
       get_moon_events = moon_events
       IF( get_moon_events )THEN
          write(*,'(''    Removing Moon events if present'')')
       ENDIF
    ENDIF
    ignore_solar = .FALSE.
    IF( PRESENT(ignore_solar_contamination) )THEN
       ignore_solar = ignore_solar_contamination
       IF( ignore_solar )THEN
          write(*,'(''    Removing Solar contaminated lines if present'')')
       ENDIF
    ENDIF
    IF( PRESENT(cal_overlap) )THEN
       cal_coef_overlap = cal_overlap
    ELSE
       cal_coef_overlap = 1.
    ENDIF
    IF( PRESENT(correct_solar_simple) )THEN
       corr_solar_simple = correct_solar_simple
    ELSE
       corr_solar_simple = .TRUE.
    ENDIF
    IF( PRESENT(new_vis_cal) )THEN
       new_vis_calibration = new_vis_cal
    ELSE
       new_vis_calibration = .TRUE.
    ENDIF
    IF( PRESENT(nedt_scl) )THEN
       nedt_scale = nedt_scl
    ELSE
       nedt_scale = 1.
    ENDIF
    IF( PRESENT(use_old_solar) )THEN
       use_old_solar_contam = use_old_solar
    ELSE
       use_old_solar_contam = .TRUE.
    ENDIF
    IF( PRESENT(noise_orbit) )THEN
       noise_all_data=noise_orbit
    ELSE
       noise_all_data=.TRUE.
    ENDIF
    IF( PRESENT(use_iasi_calibration) )THEN
       use_iasi_cal = use_iasi_calibration
       IF( use_iasi_cal .and. outData%AVHRR_No .eq. -1 )THEN
          write(*,'(''  Recalibrating AVHRR (Against IASI)'')')
       ELSE
          write(*,'(''  Recalibrating AVHRR'')')
       ENDIF
    ELSE
       use_iasi_cal = .FALSE.
       write(*,'(''  Recalibrating AVHRR'')')
    ENDIF

    ! Allocate new calibration outputs
    if( outData%dataFilled )then
       if( outData%newCalibration_There )then
          DEALLOCATE(outData%new_calib3,&
               outData%new_calib4,&
               outData%new_calib5,&
               outData%smoothPrt1,&
               outData%smoothPrt2,&
               outData%smoothPrt3,&
               outData%smoothPrt4,&
               outData%smoothPrt1Cnts,&
               outData%smoothPrt2Cnts,&
               outData%smoothPrt3Cnts,&
               outData%smoothPrt4Cnts,&
               outData%smoothPrt,&
               outData%smoothBB3,&
               outData%smoothBB4,&
               outData%smoothBB5,&
               outData%smoothSp3,&
               outData%smoothSp4,&
               outData%smoothSp5,&
               outData%Interpolated,&
               outData%solar_contamination_3B,&
               outData%solar_contamination_4,&
               outData%solar_contamination_5,&
               outData%moon_contamination,&
               outData%new_array1,&
               outData%new_array2,&
               outData%new_array3A,&
               outData%new_array3B,&
               outData%new_array4,&
               outData%new_array5,&
               outData%new_array1_error,&
               outData%new_array2_error,&
               outData%new_array3A_error,&
               outData%new_array3B_error,&
               outData%new_array4_error,&
               outData%new_array5_error,&
               STAT=STAT)
          CALL Check_IOS(STAT,'De-Allocating new calibration structure',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
    endif
    ALLOCATE(outData%new_calib3(3,outData%arraySize),&
         outData%new_calib4(3,outData%arraySize),&
         outData%new_calib5(3,outData%arraySize),&
         outData%smoothPrt1(outdata%arraySize),&
         outData%smoothPrt2(outdata%arraySize),&
         outData%smoothPrt3(outdata%arraySize),&
         outData%smoothPrt4(outdata%arraySize),&
         outData%smoothPrt1Cnts(outdata%arraySize),&
         outData%smoothPrt2Cnts(outdata%arraySize),&
         outData%smoothPrt3Cnts(outdata%arraySize),&
         outData%smoothPrt4Cnts(outdata%arraySize),&
         outData%smoothPrt(outdata%arraySize),&
         outData%smoothBB3(outdata%arraySize),&
         outData%smoothBB4(outdata%arraySize),&
         outData%smoothBB5(outdata%arraySize),&
         outData%smoothSp3(outdata%arraySize),&
         outData%smoothSp4(outdata%arraySize),&
         outData%smoothSp5(outdata%arraySize),&
         outData%Interpolated(outdata%arraySize),&
         outData%solar_contamination_3B(outdata%arraySize),&
         outData%solar_contamination_4(outdata%arraySize),&
         outData%solar_contamination_5(outdata%arraySize),&
         outData%moon_contamination(outdata%arraySize),&
         outData%new_array1(outData%nelem,outdata%arraySize),&
         outData%new_array2(outData%nelem,outdata%arraySize),&
         outData%new_array3A(outData%nelem,outdata%arraySize),&
         outData%new_array3B(outData%nelem,outdata%arraySize),&
         outData%new_array4(outData%nelem,outdata%arraySize),&
         outData%new_array5(outData%nelem,outdata%arraySize),&
         outData%new_array1_error(outData%nelem,outdata%arraySize),&
         outData%new_array2_error(outData%nelem,outdata%arraySize),&
         outData%new_array3A_error(outData%nelem,outdata%arraySize),&
         outData%new_array3B_error(outData%nelem,outdata%arraySize),&
         outData%new_array4_error(outData%nelem,outdata%arraySize),&
         outData%new_array5_error(outData%nelem,outdata%arraySize),&
         STAT=STAT)
    CALL Check_IOS(STAT,'Allocating new calibration structure',&
         'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outData%newCalibration_There = .TRUE.

    ! First see if Level1B has reported solar contamination
    printSolarWarn=.TRUE.
    DO I=1,outData%arraySize
       IF( outData%orig_solar_contamination_3B(I) .or. &
            outData%orig_solar_contamination_4(I) .or. &
            outData%orig_solar_contamination_5(I) )THEN
          IF( printSolarWarn )THEN
             WRITE(*,&
         '(''  WARNING: '')')
             WRITE(*,&
         '(''  WARNING: Solar contamination flagged in original Level 1B'')')
             WRITE(*,&
         '(''  WARNING: '')')
             printSolarWarn=.FALSE.
          ENDIF
       ENDIF
    END DO

    ! Set all calibration to good for the moment

    ! In case of level 1B remove first/last NPIXEL_PRT_SMOOTH lines because
    ! need +/- lines
    nlines_removed = 0
    I = 0
    remove_lines1: DO WHILE(.TRUE.)
       I=I+1
       if( I .gt. outData%arraySize )then
          call out_WARNING('Cannot find any valid AVHRR data',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          DO J=1,outData%arraySize
             outData%badCalibration(J) = .TRUE.
             call set_calib_bad( outData, J )
          END DO
          RETURN
       endif
       IF( .NOT. outData%badCalibration(I) )THEN
          nlines_removed = nlines_removed + 1
          if( nlines_removed .gt. NPIXEL_PRT_SMOOTH )then
             start_pos = I+1
             outData%badTime(1:start_pos) = .TRUE.
             exit remove_lines1
          endif
       endif
    END DO remove_lines1

    nlines_removed = 0
    I = outData%arraySize+1
    remove_lines2: DO WHILE(.TRUE.)
       I=I-1
       if( I .le. 0 )then
          call out_WARNING('Cannot find any valid AVHRR data',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          DO J=1,outData%arraySize
             outData%badCalibration(J) = .TRUE.
             call set_calib_bad( outData, J )
          END DO
          RETURN
       endif
       if( .NOT. outData%badCalibration(I) )then
          nlines_removed = nlines_removed + 1
          if( nlines_removed .gt. NPIXEL_PRT_SMOOTH )then
             stop_pos = I-1
             outData%badTime(stop_pos:outData%arraySize) = .TRUE.
             exit remove_lines2
          endif
       endif
    END DO remove_lines2

    ! Get noise estimates
    CALL Get_Allan_deviation( outData, start_pos, stop_pos, noise, all=noise_all_data )
    ALLOCATE(outData%noise_cnts(6,outData%arraySize),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate noise_cnts','Recalibrate_AVHRR',&
            'NOAA_LoadAVHRRLevel1B.f90')
    ENDIF
    IF( noise_all_data )THEN
       outData%noise_cnts = -1
       DO I=start_pos,stop_pos
          outData%noise_cnts(:,I) = noise(:,1)
       END DO
    ELSE
       outData%noise_cnts(:,start_pos:stop_pos) = noise
    ENDIF

    ! Now get smoothed PRT data and counts for re-calibration
    norbital = 0
    outData%orbital_temperature = 0.
    outData%start_valid = start_pos
    outData%stop_valid = stop_pos
    DO I=start_pos,stop_pos
       if( .NOT. outData%badCalibration(I) )then
          if( NAN_R .ne. outData%prt1(I) .and. &
               NAN_R .ne. outData%prt2(I) .and. &
               NAN_R .ne. outData%prt3(I) .and. &
               NAN_R .ne. outData%prt4(I) )then
             norbital = norbital + 1
             outData%orbital_temperature = outData%orbital_temperature + &
                  (outData%prt1(I) + outData%prt2(I) + outData%prt3(I) + outData%prt4(I))/4.
          endif
       endif
       minx = MAX(1,I-NPIXEL_PRT_SMOOTH)
       maxx = MIN(outData%arraySize,I+NPIXEL_PRT_SMOOTH)
       ntotal = 0
       total1 = 0.
       total2 = 0.
       total3 = 0.
       total4 = 0.
       ntotal_cnts = 0
       total1_cnts = 0.
       total2_cnts = 0.
       total3_cnts = 0.
       total4_cnts = 0.
       DO J=minx,maxx
          if( .NOT. outData%badCalibration(J) )then
             if( NAN_R .ne. outData%prt1(J) .and. &
                  NAN_R .ne. outData%prt2(J) .and. &
                  NAN_R .ne. outData%prt3(J) .and. &
                  NAN_R .ne. outData%prt4(J) )then
                total1 = total1 + outData%prt1(J)
                total2 = total2 + outData%prt2(J)
                total3 = total3 + outData%prt3(J)
                total4 = total4 + outData%prt4(J)
                ntotal = ntotal+1
             endif
             if( NAN_R .ne. outData%prt1(J) .and. &
                  NAN_R .ne. outData%prt2(J) .and. &
                  NAN_R .ne. outData%prt3(J) .and. &
                  NAN_R .ne. outData%prt4(J) )then
                total1_cnts = total1_cnts + SUM(outData%prt1CountsAll(:,J))
                total2_cnts = total2_cnts + SUM(outData%prt2CountsAll(:,J))
                total3_cnts = total3_cnts + SUM(outData%prt3CountsAll(:,J))
                total4_cnts = total4_cnts + SUM(outData%prt4CountsAll(:,J))
                ntotal_cnts = ntotal_cnts+3
             endif
          endif
       END DO
       if( 0 .eq. ntotal )then
          outData%smoothPrt1(I) = NAN_R
          outData%smoothPrt2(I) = NAN_R
          outData%smoothPrt3(I) = NAN_R
          outData%smoothPrt4(I) = NAN_R
       else
          outData%smoothPrt1(I) = total1 / ntotal
          outData%smoothPrt2(I) = total2 / ntotal
          outData%smoothPrt3(I) = total3 / ntotal
          outData%smoothPrt4(I) = total4 / ntotal
          outData%smoothPrt(I) = (outData%smoothPrt1(i) + &
               outData%smoothPrt2(i) + &
               outData%smoothPrt3(i) + &
               outData%smoothPrt4(i))/4.
       endif
       if( 0 .eq. ntotal_cnts )then
          outData%smoothPrt1Cnts(I) = NAN_R
          outData%smoothPrt2Cnts(I) = NAN_R
          outData%smoothPrt3Cnts(I) = NAN_R
          outData%smoothPrt4Cnts(I) = NAN_R
       else
          outData%smoothPrt1Cnts(I) = total1_cnts / ntotal_cnts
          outData%smoothPrt2Cnts(I) = total2_cnts / ntotal_cnts
          outData%smoothPrt3Cnts(I) = total3_cnts / ntotal_cnts
          outData%smoothPrt4Cnts(I) = total4_cnts / ntotal_cnts
       endif
    END DO
    IF( 5000 .gt. norbital )THEN
       write(*,*)'Norbital = ',norbital
       call out_WARNING('Not really enough data to determine the orbital temperature',&
            'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    outData%orbital_temperature = outData%orbital_temperature / norbital

    ! Now get smoothed Counts
    DO I=start_pos,stop_pos
       minx = MAX(1,I-NPIXEL_PRT_SMOOTH)
       maxx = MIN(outData%arraySize,I+NPIXEL_PRT_SMOOTH)
       ntotal = 0
       ntotal0 = 0
       total1 = 0.
       total2 = 0.
       total3 = 0.
       DO J=minx,maxx
          if( .NOT. outData%badCalibration(J) )then
             if( NAN_R .ne. outData%bb3(J) )then
                total1 = total1 + outData%bb3(J)
                ntotal0 = ntotal0 + 1
             endif
             if( NAN_R .ne. outData%bb4(J) .and. &
                  NAN_R .ne. outData%bb5(J) )then
                total2 = total2 + outData%bb4(J)
                total3 = total3 + outData%bb5(J)
                ntotal = ntotal+1
             endif
          endif
       END DO
       if( 0 .eq. ntotal0 )then
          outData%smoothBB3(I) = NAN_R
       else
          outData%smoothBB3(I) = total1/ntotal0
       endif
       if( 0 .eq. ntotal )then
          outData%smoothBB4(I) = NAN_R
          outData%smoothBB5(I) = NAN_R
       else
          outData%smoothBB4(I) = total2/ntotal
          outData%smoothBB5(I) = total3/ntotal
       endif
       ntotal = 0
       ntotal0 = 0
       total1 = 0.
       total2 = 0.
       total3 = 0.
       DO J=minx,maxx
          if( .NOT. outData%badCalibration(J) )then
             if( NAN_R .ne. outData%sp3(J) )then
                total1 = total1 + outData%sp3(J)
                ntotal0 = ntotal0 + 1
             endif
             if( NAN_R .ne. outData%sp4(J) .and. &
                  NAN_R .ne. outData%sp5(J) )then
                total2 = total2 + outData%sp4(J)
                total3 = total3 + outData%sp5(J)
                ntotal = ntotal+1
             endif
          endif
       END DO
       if( 0 .eq. ntotal0 )then
          outData%smoothSp3(I) = NAN_R
       else
          outData%smoothSp3(I) = total1/ntotal0
       endif
       if( 0 .eq. ntotal )then
          outData%smoothSp4(I) = NAN_R
          outData%smoothSp5(I) = NAN_R
       else
          outData%smoothSp4(I) = total2/ntotal
          outData%smoothSp5(I) = total3/ntotal
       endif
    END DO

    ! Find moon events and remove from calibration

    IF( get_moon_events )THEN
       call Find_Moon_Events( start_pos, stop_pos, outData )
    ENDIF

    ! Now with smoothed data we can do the calibration

    ! Get raw calibration coefficients
    call Get_Calib_Coefficients_CCI( outData%time(start_pos), outData%AVHRR_No,&
         coefs1, coefs2, coefs3, ncoefs, coefs_frac1, coefs_frac2, &
         cal_coef_overlap, use_iasi_cal )

    ! Note A and B band coefficients in the level1B files are stored
    ! for conversion of Radiance->temperature as
    !
    !       temperature = C2*nuc/log(1+(c1*nuc3/radiance)
    !   and
    !       output_temperature = A + B*temperature
    !
    ! Recalibration coefficients are stored in the temperature->radiance
    ! as the default sense so have to be converted as
    !
    !       new_a = -a/b
    !       new_b = 1/b
    !
    ! so to keep consistent with level 1B (KLM) parameters need nuc/a/b
    ! to work as
    !
    !       Tstar = (temp - a)/b
    !       Radiance = C1*nuc3/(exp(C2*nuc/Tstar)-1.)
    !

    DO I=1,ncoefs
       newa = -coefs1(6,I)/coefs1(7,I)
       newb = 1./coefs1(7,I)
       coefs1(6,I) = newa
       coefs1(7,I) = newb
       newa = -coefs2(6,I)/coefs2(7,I)
       newb = 1./coefs2(7,I)
       coefs2(6,I) = newa
       coefs2(7,I) = newb
       newa = -coefs3(6,I)/coefs3(7,I)
       newb = 1./coefs3(7,I)
       coefs3(6,I) = newa
       coefs3(7,I) = newb
    END DO

    ! Copy to output - take the first one for output only
    outData%new_cal_coefs3 = coefs1(:,1)
    outData%new_cal_coefs4 = coefs2(:,1)
    outData%new_cal_coefs5 = coefs3(:,1)

    ! Copy to radiative coef entry if present
    IF( PRESENT(avhrr_rad_coefs) )THEN
       avhrr_rad_coefs%nuc = (/coefs1(5,1),coefs2(5,1),coefs3(5,1)/)
       avhrr_rad_coefs%aVal = (/coefs1(6,1),coefs2(6,1),coefs3(6,1)/)
       avhrr_rad_coefs%bVal = (/coefs1(7,1),coefs2(7,1),coefs3(7,1)/)
    
   
    end if 
   
    ! Now get new calibration coefficients
    outData%Interpolated(start_pos:stop_pos) = .FALSE.

    ! Allocate gain calculations
    IF( ALLOCATED(gain_array3) )THEN
       DEALLOCATE(gain_array3,gain_array4,gain_array5)
    ENDIF
    IF( ALLOCATED(rad_array3) )THEN
       DEALLOCATE(rad_array3,rad_array4,rad_array5)
    ENDIF
    IF( ALLOCATED(rad_array3_correl) )THEN
       DEALLOCATE(rad_array3_correl,rad_array4_correl,rad_array5_correl)
    ENDIF
    IF( ALLOCATED(good_data) )THEN
       DEALLOCATE(good_data,solar_data)
    ENDIF
    new_ndata = stop_pos-start_pos+1
    IF( 0 .ge. new_ndata )THEN
       call out_WARNING('Cannot allocate gain_array/rad_array - ndata==0',&
            'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       DO J=1,outData%arraySize
          ! this assumes that we want to use the new calibration if we call
          ! this routine
          outData%badCalibration(J) = .TRUE.
          call set_calib_bad( outData, J )
       END DO
       RETURN
    ENDIF
    ALLOCATE( gain_array3(new_ndata,2), &
         gain_array4(new_ndata,2), &
         gain_array5(new_ndata,2), &
         rad_array3(new_ndata), &
         rad_array4(new_ndata), &
         rad_array5(new_ndata), &
         rad_array3_correl(new_ndata), &
         rad_array4_correl(new_ndata), &
         rad_array5_correl(new_ndata), &
         good_data(new_ndata),&
         solar_data(new_ndata),&
         new_calib3(3,outData%arraySize,2),&
         new_calib4(3,outData%arraySize,2),&
         new_calib5(3,outData%arraySize,2),&
         STAT=STAT )
    IF( 0 .ne. STAT )THEN
       call out_ERROR('Cannot allocate gain_array/rad_array',&
            'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    gain_array3 = NAN_R
    gain_array4 = NAN_R
    gain_array5 = NAN_R
    rad_array3 = 0.
    rad_array4 = 0.
    rad_array5 = 0.
    rad_array3_correl = 0.
    rad_array4_correl = 0.
    rad_array5_correl = 0.
    good_data = .TRUE.
    solar_data = .TRUE.

    ! Set calibration parameters to bad to begin
    ! Note that new_calib* has extra dimensions to deal with overlap
    ! in time between different calibration parameters
    outData%new_calib3 = NAN_R
    outData%new_calib4 = NAN_R
    outData%new_calib5 = NAN_R
    new_calib3 = NAN_R
    new_calib4 = NAN_R
    new_calib5 = NAN_R

    ! Get 'raw' gains.  If not correction this is the answer
    J3 = 0
    J4 = 0
    J5 = 0
    mean3 = 0.
    mean4 = 0.
    mean5 = 0.
    DO I=start_pos,stop_pos
       ! Check if this is data we want to derive model with
       if( NAN_R .ne. outData%smoothPrt(I) )then
!          minPRT = 1e30
!          minPRT = MIN(minPRT,outData%smoothPrt1(i))
!          minPRT = MIN(minPRT,outData%smoothPrt2(i))
!          minPRT = MIN(minPRT,outData%smoothPrt3(i))
!          minPRT = MIN(minPRT,outData%smoothPrt4(i))
          correlPRT = outData%smoothPrt(i)
          if( NAN_R .ne. outData%smoothBB3(I) .and. &
               NAN_R .ne. outData%smoothSp3(I) .and. &
               0 .lt. outData%smoothBB3(I) .and. &
               0 .lt. outData%smoothSp3(I) )then
             R_ICT = convertRadiance( outData%smoothPrt(i),coefs1(5,1),&
                  coefs1(6,1),coefs1(7,1))
             !print*, "rict", r_ict
             counts = (outData%smoothSp3(I) - outData%smoothBB3(I))
             J3 = J3+1
             POS = I-start_pos+1
             DO J=1,ncoefs
                gain_array3(POS,J) = ((eta_ict+coefs1(2,J))*R_ICT - &
                     coefs1(3,J) - coefs1(4,J)*counts*counts)/&
                     counts
                IF( J .eq. 1 )THEN
                   mean3 = mean3 + gain_array3(POS,J)
                   rad_array3(POS) = R_ICT
                   R_ICT = convertRadiance( correlPRT,coefs1(5,1),&
                        coefs1(6,1),coefs1(7,1))
                   rad_array3_correl(POS) = R_ICT
                ENDIF
                ! Convert to level 1B type coefficients
                new_calib3(1,I,J)=&
                     (coefs1(1,J)+gain_array3(POS,J)*outData%smoothSp3(i)+&
                     coefs1(4,J)*outData%smoothSp3(i)*outData%smoothSp3(i)) + &
                     coefs1(8,J)*outData%orbital_temperature
                new_calib3(2,I,J) = -(gain_array3(POS,J) + &
                     2*coefs1(4,J)*outData%smoothSp3(i))
                new_calib3(3,I,J) = coefs1(4,J)
             END DO
          else
             new_calib3(:,I,:) = NAN_R
          endif
          if( NAN_R .ne. outData%smoothBB4(I) .and. &
               NAN_R .ne. outData%smoothSp4(I) .and. &
               0 .lt. outData%smoothBB4(I) .and. &
               0 .lt. outData%smoothSp4(I) )then
             R_ICT = convertRadiance( outData%smoothPrt(i),coefs2(5,1),&
                  coefs2(6,1),coefs2(7,1))             
             counts = (outData%smoothSp4(I) - outData%smoothBB4(I))
             J4 = J4+1
             POS = I-start_pos+1
             DO J=1,ncoefs
                gain_array4(POS,J) = ((eta_ict+coefs2(2,J))*R_ICT - &
                     coefs2(3,J) - coefs2(4,J)*counts*counts)/&
                     counts
                IF( J .eq. 1)THEN
                   mean4 = mean4 + gain_array4(POS,J)
                   rad_array4(POS) = R_ICT
                   R_ICT = convertRadiance( correlPRT,coefs2(5,1),&
                        coefs2(6,1),coefs2(7,1))
                   rad_array4_correl(POS) = R_ICT
                ENDIF
                ! Convert to level 1B type coefficients
                new_calib4(1,I,J)=&
                     (coefs2(1,J)+gain_array4(POS,J)*outData%smoothSp4(i)+&
                     coefs2(4,J)*outData%smoothSp4(i)*outData%smoothSp4(i)) + &
                     coefs2(8,J)*outData%orbital_temperature
                new_calib4(2,I,J) = -(gain_array4(POS,J) + &
                     2*coefs2(4,J)*outData%smoothSp4(i))
                new_calib4(3,I,J) = coefs2(4,J)
             END DO
          else
             new_calib4(:,I,:) = NAN_R
          endif
          if( NAN_R .ne. outData%smoothBB5(I) .and. &
               NAN_R .ne. outData%smoothSp5(I) .and. &
               0 .lt. outData%smoothBB5(I) .and. &
               0 .lt. outData%smoothSp5(I) )then

             R_ICT = convertRadiance( outData%smoothPrt(i),coefs3(5,1),&
                  coefs3(6,1),coefs3(7,1) )
             counts = (outData%smoothSp5(I) - outData%smoothBB5(I))
             J5 = J5+1
             POS = I-start_pos+1
             DO J=1,ncoefs
                gain_array5(POS,J) = ((eta_ict+coefs3(2,J))*R_ICT - &
                     coefs3(3,J)-coefs3(4,J)*counts*counts)/&
                     counts
                IF( J .eq. 1 )THEN
                   mean5 = mean5 + gain_array5(POS,J)
                   rad_array5(POS) = R_ICT
                   R_ICT = convertRadiance( correlPRT,coefs3(5,1),&
                        coefs3(6,1),coefs3(7,1))
                   rad_array5_correl(POS) = R_ICT
                ENDIF
                ! Convert to level 1B type coefficients
                new_calib5(1,I,J)=&
                     (coefs3(1,J)+gain_array5(POS,J)*outData%smoothSp5(i)+&
                     coefs3(4,J)*outData%smoothSp5(i)*outData%smoothSp5(i)) + &
                     coefs3(8,J)*outData%orbital_temperature
                new_calib5(2,I,J) = -(gain_array5(POS,J) + &
                     2*coefs3(4,J)*outData%smoothSp5(i))
                new_calib5(3,I,J) = coefs3(4,J)
             END DO
          else
             new_calib5(:,I,:) = NAN_R
          endif
       else
          new_calib3(:,I,:) = NAN_R
          new_calib4(:,I,:) = NAN_R
          new_calib5(:,I,:) = NAN_R
       endif
    END DO
    IF( 0 .lt. J3 )THEN
       mean3 = mean3 / J3
    ENDIF
    IF( 0 .lt. J4 )THEN
       mean4 = mean4 / J4
    ENDIF
    IF( 0 .lt. J5 )THEN
       mean5 = mean5 / J5
    ENDIF
    if( 0 .gt. mean3 )then
       sign3 = -1
    else
       sign3 = 1
    endif
    if( 0 .gt. mean4 )then
       sign4 = -1
    else
       sign4 = 1
    endif
    if( 0 .gt. mean5 )then
       sign5 = -1
    else
       sign5 = 1
    endif

    !
    ! make sure that solar contamination (new) is set to FALSE
    !
    outData%solar_contamination_3B = .FALSE.
    outData%solar_contamination_4 = .FALSE.
    outData%solar_contamination_5 = .FALSE.

    ! If want to interpolate daytime from nighttime and/or
    ! correct for any solar contamination
    pos = outData%nelem/2
    IF( interp_daynight .or. correct_solar_contam .or. &
         correct_solar_contam_all)THEN

       night_ok = .TRUE.
       SolarStat = 0
       outData%solar_contamination_failure = .FALSE.
       IF( correct_solar_contam .or. correct_solar_contam_all )THEN
          ! Do this for 3.9 micron channel since this is most sensitive
          ! to it
          SolarStat = Remove_Solar_entries( stop_pos-start_pos+1,&
               outData%solZA(pos,start_pos:stop_pos),&
               gain_array3(:,1), rad_array3_correl,&
               outData%solar_contamination_3b, good_data, night_ok )
          IF( 0 .eq. SolarStat )THEN
             solar_data = good_data
          ELSE
             ! Set possible solar contamiation data to bad -problem
             ! with Remove_Solar_Entries
             DO I=start_pos,stop_pos
                IF( .not. good_data(I-start_pos+1) )THEN
                   outData%badCalibration(I) = .TRUE.
                   call set_calib_bad( outData, I )
                ENDIF
             END DO
          ENDIF
          IF( ignore_solar )THEN
             printSolarWarn=.TRUE.
             DO I=start_pos,stop_pos
                IF( .not. good_data(I-start_pos+1) )THEN
                   outData%badCalibration(I) = .TRUE.
                   call set_calib_bad( outData, I )
                   IF( printSolarWarn )THEN
                      WRITE(*,&
         '(''  WARNING: '')')
                      WRITE(*,&
         '(''  WARNING: Solar contamination removed from data'')')
                      WRITE(*,&
         '(''  WARNING: '')')
                      printSolarWarn=.FALSE.
                   ENDIF
                ENDIF
             END DO
          ENDIF
       endif
       ! If we have enough data to use nighttime only
       if( night_ok .and. interp_daynight )then
          do i=start_pos,stop_pos
             if(outData%solZA(pos,I) .le. 90)then
                good_data(i-start_pos+1) = .FALSE.
             endif
          end do
       endif

       ! Correct earthshine
       if( correct_earthshine_data .and. 0 .eq. SolarStat )then
          ! Get 3.7 mu radiance correction factor use in gain calculation
          CALL Earthshine_Correct( outData, coefs1, new_ndata, good_data, &
               start_pos, stop_pos, night_ok, gain_array3 )
       endif

       IF( 0 .eq. SolarStat .and. .not.ignore_solar )THEN
          SolarStat3 = Get_Poly_Coefs( new_ndata, rad_array3, gain_array3, &
               good_data, stdev3, maxdev3, poly_coefs3 )
          IF( correct_solar_contam_all )THEN
             SolarStat4 = Get_Poly_Coefs( new_ndata, rad_array4, gain_array4, &
                  good_data, stdev4, maxdev4, poly_coefs4 )
             SolarStat5 = Get_Poly_Coefs( new_ndata, rad_array5, gain_array5, &
                  good_data, stdev5, maxdev5, poly_coefs5 )
          ELSE
             SolarStat4 = 0
             SolarStat5 = 0
          ENDIF

          ! Only recalculate if OK
          IF( 0 .eq. SolarStat3 .and. 0 .eq. SolarStat4 .and. &
               0 .eq. SolarStat5 )THEN

             ! If using nighttime only then use interpolated data for all data
             if( night_ok .and. interp_daynight )then
                good_data = .FALSE.
             endif

             outData%gain_stdev(1) = stdev3
             outData%gain_maxdev(1) = maxdev3
             outData%poly_coefs3 = poly_coefs3
             IF( correct_solar_contam_all )THEN
                outData%gain_stdev(2) = stdev4
                outData%gain_stdev(3) = stdev5
                outData%gain_maxdev(2) = maxdev4
                outData%gain_maxdev(3) = maxdev5
                outData%poly_coefs4 = poly_coefs4
                outData%poly_coefs5 = poly_coefs5
             ELSE
                outData%gain_stdev(2) = NAN_R
                outData%gain_stdev(3) = NAN_R
                outData%gain_maxdev(2) = NAN_R
                outData%gain_maxdev(3) = NAN_R
                outData%poly_coefs4 = NAN_R
                outData%poly_coefs5 = NAN_R
             ENDIF

             outData%solar_contamination_3b = .FALSE.
             outData%solar_contamination_4 = .FALSE.
             outData%solar_contamination_5 = .FALSE.
             DO I=start_pos,stop_pos
                if( NAN_R .ne. outData%smoothPrt(I) .and. &
                     .not.good_data(i-start_pos+1) )then
                   IF( .not.solar_data(i-start_pos+1) )THEN
                      outData%solar_contamination_3b(I) = .TRUE.
                      IF( correct_solar_contam_all )THEN
                         outData%solar_contamination_4(I) = .TRUE.
                         outData%solar_contamination_5(I) = .TRUE.
                      ENDIF
                   ENDIF
                   outData%interpolated(I) = .TRUE.
                   ! 3.7 micron
                   R_ICT = convertRadiance( outData%smoothPrt(i),coefs1(5,1),&
                        coefs1(6,1),coefs1(7,1))
                   gain = sign3*sqrt(poly_coefs3(1) + poly_coefs3(2)*R_ICT + &
                        poly_coefs3(3)*R_ICT*R_ICT)
                   ! Convert to level 1B type coefficients
                   new_calib3(1,I,1)=&
                        (coefs1(1,1)+gain*outData%smoothSp3(i)+&
                        coefs1(4,1)*outData%smoothSp3(i)*outData%smoothSp3(i))
                   new_calib3(2,I,1) = -(gain + &
                        2*coefs1(4,1)*outData%smoothSp3(i))
                   new_calib3(3,I,1) = coefs1(4,1)
                   ! 11 micron
                   R_ICT = convertRadiance( outData%smoothPrt(i),coefs2(5,1),&
                        coefs2(6,1),coefs2(7,1))
                   IF( correct_solar_contam_all )THEN
                      gain = sign4*sqrt(poly_coefs4(1) + &
                           poly_coefs4(2)*R_ICT + &
                           poly_coefs4(3)*R_ICT*R_ICT)
                      ! Convert to level 1B type coefficients
                      new_calib4(1,I,1)=&
                           (coefs2(1,1)+gain*outData%smoothSp4(i)+&
                           coefs2(4,1)*outData%smoothSp4(i)*outData%smoothSp4(i))
                      new_calib4(2,I,1) = -(gain + &
                           2*coefs2(4,1)*outData%smoothSp4(i))
                      new_calib4(3,I,1) = coefs2(4,1)
                      ! 12 micron
                      R_ICT = convertRadiance( outData%smoothPrt(i),coefs3(5,1),&
                           coefs3(6,1),coefs3(7,1))
                      gain = sign5*sqrt(poly_coefs5(1) + &
                           poly_coefs5(2)*R_ICT + &
                           poly_coefs5(3)*R_ICT*R_ICT)
                      ! Convert to level 1B type coefficients
                      new_calib5(1,I,1)=&
                           (coefs3(1,1)+gain*outData%smoothSp5(i)+&
                           coefs3(4,1)*outData%smoothSp5(i)*outData%smoothSp5(i))
                      new_calib5(2,I,1) = -(gain + &
                           2*coefs3(4,1)*outData%smoothSp5(i))
                      new_calib5(3,I,1) = coefs3(4,1)
                   ENDIF
                ENDIF
             END DO
          ELSE
             CALL out_Warning('Cannot do solar contamination correction',&
                  'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ENDIF
       ENDIF

       ! Correct earthshine
    ELSE IF( correct_earthshine_data )THEN
       ! Get 3.7 mu radiance correction factor use in gain calculation
       CALL Earthshine_Correct( outData, coefs1(:,1), new_ndata, good_data, &
            start_pos, stop_pos, .TRUE., gain_array3(:,1) )
    ENDIF

    !
    ! If there are operational tests done use those if new methods doesn't find any...
    !
    operational_solar=.FALSE.
    DO I=start_pos,stop_pos
       IF( outData%orig_solar_contamination_3B(I) )THEN
          operational_solar=.TRUE.
       ENDIF
    END DO
    !
    ! If just a simple test for solar contamination do it here
    !
    reset_gain=.FALSE.
    IF( corr_solar_simple .and. .not. correct_solar_contam .and. &
         .not. correct_solar_contam_all )THEN
       WRITE(*,'(''... Trying to find solar contamination (new) ...'')')

       !
       ! First try and use new algo to get solar contamination
       !
       ! Smooth solza array
       !
       ALLOCATE(smooth_solza(outData%arraySize),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL out_Error('Cannot allocate smooth_solza',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       pos = outData%nelem/2
       DO I=start_pos,stop_pos
          minI1 = MAX(I-1000,1)
          maxI1 = MIN(i+1000,stop_pos)
          smooth_solza(I) = 0.
          nmean1 = 0
          DO J=minI1,maxI1
             IF( 0 .le. outData%solZA(pos,J) )THEN
                smooth_solza(I) = smooth_solza(I) + outData%solZA(pos,J)
                nmean1 = nmean1+1
             ENDIF
          END DO
          IF( 0 .lt. nmean1 )THEN
             smooth_solza(I) = smooth_solza(I) / nmean1
          ELSE
             smooth_solza(I) = -1e30
          ENDIF
       END DO
       ! 
       !
       ! Find ascending (increasing) sections of solZA
       !
       ALLOCATE(solza_ascend(outData%arraySize),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL out_Error('Cannot allocate solza_ascend',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       solza_ascend(:) = .FALSE.
       flat_start = -1
       flat_stop = -1
       DO I=start_pos,stop_pos
          IF( I-20 .lt. start_pos )THEN
             minI1 = start_pos
             maxI1 = start_pos+19
             minI2 = start_pos+21
             maxI2 = start_pos+40
          ELSE IF( I+20 .gt. stop_pos )THEN
             minI1 = stop_pos-40
             maxI1 = stop_pos-21
             minI2 = stop_pos-19
             maxI2 = stop_pos
          ELSE
             minI1 = I-20
             maxI1 = I-1
             minI2 = I+1
             maxI2 = I+20
          ENDIF
          mean1 = 0.
          nmean1 = 0
          DO J=minI1,maxI1
             IF( 0 .le. smooth_solZA(J) )THEN
                nmean1 = nmean1 + 1
                mean1 = mean1 + smooth_solZA(J)
             ENDIF
          END DO
          IF( nmean1 .gt. 0 )THEN
             mean1 = mean1/nmean1
          ENDIF
          mean2 = 0.
          nmean2 = 0
          DO J=minI2,maxI2
             IF( 0 .le. smooth_solZA(J) )THEN
                nmean2 = nmean2 + 1
                mean2 = mean2 + smooth_solZA(J)
             ENDIF
          END DO
          IF( nmean2 .gt. 0 )THEN
             mean2 = mean2/nmean2
          ENDIF
          IF( nmean1 .gt. 0 .and. nmean2 .gt. 0 )THEN
             IF( mean2 .gt. mean1 )THEN
                solza_ascend(I) = .TRUE.
             ELSE IF( mean2 .eq. mean1 .and. mean1 .gt. 125.0 )THEN
                IF( -1 .eq. flat_start )THEN
                   flat_start = I
                   flat_stop = I
                ELSE IF( flat_stop .eq. I-1 )THEN
                   flat_stop = I
                ENDIF
             ENDIF
          ENDIF
       END DO
       IF( flat_start .ne. -1 .and. flat_stop-flat_start .gt. 100 )THEN
          flat_stop = (flat_start+flat_stop)/2
          DO I=flat_start,flat_stop
             solza_ascend(I) = .TRUE.
          ENDDO
       ENDIF

       !
       ! From daytime data get mean 3.7 micron channel gain
       !
       ! Daytime defined as solZA < 80 degrees here to miss limb data
       !
       !
       ! Find minimum solza as start and increment by 5 degrees
       ! and check next 5 degrees to look for event
       !
       ALLOCATE(Used(outData%arraySize),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL out_Error('Cannot allocate Used',&
               'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       Used = 0
       nUsed = 0
       UsedPos = -1
       ok=.TRUE.
       J=0
       DO WHILE(ok)
          POS = -1
          minza = 1e30
          DO I=start_pos,stop_pos
             IF( .not. solza_ascend(I) )THEN
                IF( smooth_solZA(I) .lt. minza .and. Used(I) .eq. 0 )THEN
                   minza = smooth_solZA(I)
                   POS = I
                ENDIF
             ENDIF
          END DO
          IF( 0 .lt. POS )THEN
             J=J+1
             IF( J .gt. 100 )THEN
                ok = .FALSE.
                J=100
             ELSE
                !
                ! Track back to change/start
                !
                UsedLoop: DO I=POS,start_pos,-1
                   IF( .not. solza_ascend(I) )THEN
                      Used(I) = POS
                      nUsed(J) = nUsed(J) + 1
                      UsedPos(J) = POS
                   ELSE
                      EXIT UsedLoop
                   ENDIF
                END DO UsedLoop                
             ENDIF
          ELSE
             ok = .FALSE.
          ENDIF
       END DO
       !
       ! Find longest stretch of data and use that as POS
       !
       n=0
       POS = -1
       DO I=1,J
          IF( nUsed(I) .gt. n )THEN
             n = nUsed(I)
             POS = UsedPos(I)
          ENDIF
       END DO
       IF( 0 .ge. POS )THEN
          CALL Gbcs_Warning(.TRUE.,'SolarContam: No minza found',&
              'Reallocate_AVHRR','NOAA_LoadAVHRRLevel1B.f90') 
          DEALLOCATE(Used)
       ELSE
          minza = smooth_solza(POS)
          !
          ! Make sure there is distance to 80 degrees 
          !
          IF( minza .le. 50 )THEN
             LL = 4
          ELSE IF( minza .le. 60 )THEN
             LL = 3
          ELSE
             LL = 2
          ENDIF
          FindZALoop: DO KK=LL,14
             JJ=KK
             n = 0
             mean = 0.
             M2 = 0.
             n2 = 0
             mean_2 = 0.
             M2_2 = 0.
             maxza = minza + KK*5.
             IF( maxza .gt. 100 )THEN
                EXIT FindZALoop
             ENDIF
             FindZALoop1: DO I=POS,start_pos,-1
                IF( .not. solza_ascend(I) )THEN
                   IF( smooth_solZA(I) .lt. maxza )THEN                   
                      K = I-start_pos+1
                      IF( gain_array3(K,1) .gt. 0. )THEN
                         n=n+1
                         delta = gain_array3(K,1) - mean
                         mean = mean + delta/n
                         M2 = M2 + delta*(gain_array3(K,1)-mean)
                      ENDIF
                   ELSE
                      EXIT FindZALoop1
                   ENDIF
                ENDIF
             END DO FindZALoop1
             FindZALoop2: DO I=POS,start_pos,-1
                IF( .not. solza_ascend(I) )THEN
                   IF( smooth_solZA(I) .gt. maxza .and. &
                        smooth_solZA(I) .lt. maxza+5. )THEN
                      K = I-start_pos+1
                      IF( gain_array3(K,1) .gt. 0. )THEN
                         n2=n2+1
                         delta = gain_array3(K,1) - mean_2
                         mean_2 = mean_2 + delta/n2
                         M2_2 = M2_2 + delta*(gain_array3(K,1)-mean_2)
                      ENDIF
                   ENDIF
                ENDIF
             END DO FindZALoop2
             IF( n .gt. 0 .and. n2 .gt. 0 )THEN
                M2 = sqrt(M2/(n-1))
                IF( 3 .lt. (mean-mean_2)/M2 )THEN
                   EXIT FindZALoop
                ENDIF
             ENDIF
          END DO FindZALoop
          WRITE(*,'(''...minZA/maxZA = '',F6.2,''/'',F6.2)')minza,maxza
          IF( JJ .eq. 14 )THEN
             maxza=80.
          ENDIF
          !
          ! Now get maxima (or minima actually)
          !
          midpoint = outData%nelem/2
          DO J=1,ncoefs
             n = 0
             mean = 0.
             M2 = 0.
             n2 = 0
             mean_2 = 0.
             M2_2 = 0.
             DO I=start_pos,stop_pos
                K=I-start_pos+1
                IF( solza_ascend(I) )THEN
                   ! Use nighttime data
                   IF( smooth_solZA(I) .ge. 90. )THEN
                      !
                      ! remove possible 2nd spike
                      !
                      IF( smooth_solZA(I) .lt. 100. .or. &
                           smooth_solZA(I) .gt. 120 )THEN
                         IF( gain_array3(K,J) .gt. 0. )THEN
                            n=n+1
                            delta = gain_array3(K,J) - mean
                            mean = mean + delta/n
                            M2 = M2 + delta*(gain_array3(K,J)-mean)
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
                IF( .not. solza_ascend(I) )THEN
                   ! Use daytime data
                   IF( smooth_solZA(I) .gt. minza .and. &
                        smooth_solZA(I) .lt. maxza)THEN
                      IF( gain_array3(K,J) .gt. 0. )THEN
                         n2=n2+1
                         delta = gain_array3(K,J) - mean_2
                         mean_2 = mean_2 + delta/n2
                         M2_2 = M2_2 + delta*(gain_array3(K,J)-mean_2)
                      ENDIF
                   ENDIF
                ENDIF
             END DO
             reset_gain = .FALSE.
             main_event = .TRUE.
             IF( (n .gt. 1000 .and. M2 .ne. 0.) .or. (n2 .gt. 1000 .and. M2_2 .ne. 0.) )THEN
                !
                ! Now detect if there is solar contamination by looking for 10 stdev values +ve
                !
                max_stdev = 0.
                max_stdev_2 = 0.
                M2 = sqrt(M2/(n-1))             
                M2_2 = sqrt(M2_2/(n2-1))
                IF( M2 .gt. 0 )THEN
                   WRITE(*,'(''...Mean = '',E15.8,'' M2='',E15.8)')mean,M2
                ENDIF
                IF( M2_2 .gt. 0 )THEN
                   WRITE(*,'(''...2) Mean = '',E15.8,'' M2='',E15.8)')mean_2,M2_2
                ENDIF
                ok = .TRUE.
                DO WHILE(ok)
                   gain_val = 1e30
                   pos = -1
                   gain_val_2 = 1e30
                   pos_2 = -1
                   DO I=start_pos,stop_pos
                      K=I-start_pos+1
                      IF( gain_array3(K,J) .gt. 0. )THEN
                         IF( solza_ascend(I) .and. outData%solza(midpoint,I) .gt. 90.&
                              .and. n .gt. 1000 )THEN
                            IF( (mean-gain_array3(K,J))/M2 .gt. max_stdev )THEN
                               max_stdev = (mean-gain_array3(K,J))/M2
                            ENDIF
                            IF( gain_array3(K,J) .lt. mean - 10*M2 )THEN
                               reset_gain = .TRUE.
                               IF( gain_array3(K,J) .lt. gain_val )THEN
                                  gain_val = gain_array3(K,J)
                                  pos=I
                               ENDIF
                            ENDIF
                         ELSE IF( n2 .gt. 1000 )THEN
                            IF( (mean_2-gain_array3(K,J))/M2_2 .gt. max_stdev_2 )THEN
                               max_stdev_2 = (mean_2-gain_array3(K,J))/M2_2
                            ENDIF
                            IF( gain_array3(K,J) .lt. mean_2 - 6*M2_2 )THEN
                               reset_gain = .TRUE.
                               IF( gain_array3(K,J) .lt. gain_val_2 )THEN
                                  gain_val_2 = gain_array3(K,J)
                                  pos_2=I
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDIF
                   END DO
                   IF( 0 .lt. pos )THEN
                      WRITE(*,'(''...Found Max stdev from mean (solar)='',F10.6,''...'')')&
                           max_stdev
                   ENDIF
                   IF( 0 .lt. pos_2 )THEN
                      WRITE(*,'(''...Found Max stdev (2) from mean (solar)='',F10.6,''...'')')&
                           max_stdev_2
                   ENDIF
                   
                   IF( n .lt. 1000 )THEN
                      gain_val = 1e30
                      max_stdev = 0.
                      pos = -1
                      DO I=start_pos,stop_pos
                         K=I-start_pos+1
                         IF( gain_array3(K,J) .gt. 0. )THEN
                            IF( solza_ascend(I) .and. outData%solza(midpoint,I) .gt. 90. )THEN
                               IF( (mean_2-gain_array3(K,J))/M2_2 .gt. max_stdev )THEN
                                  max_stdev = (mean_2-gain_array3(K,J))/M2_2
                               ENDIF
                               IF( gain_array3(K,J) .lt. mean - 10*M2_2 )THEN
                                  reset_gain = .TRUE.
                                  IF( gain_array3(K,J) .lt. gain_val )THEN
                                     gain_val = gain_array3(K,J)
                                     pos=I
                                  ENDIF
                               ENDIF
                            ENDIF
                         ENDIF
                      END DO
                   ENDIF
                   IF( n2 .lt. 1000 )THEN
                      gain_val = 1e30
                      pos_2 = -1
                      max_stdev_2 = 0.
                      DO I=start_pos,stop_pos
                         K=I-start_pos+1
                         IF( gain_array3(K,J) .gt. 0. )THEN
                            IF( .not. (solza_ascend(I) .and. outData%solza(midpoint,I) .gt. 90.) )THEN
                               IF( (mean-gain_array3(K,J))/M2 .gt. max_stdev )THEN
                                  max_stdev = (mean-gain_array3(K,J))/M2
                               ENDIF
                               IF( gain_array3(K,J) .lt. mean - 5*M2 )THEN
                                  reset_gain = .TRUE.
                                  IF( gain_array3(K,J) .lt. gain_val )THEN
                                     gain_val = gain_array3(K,J)
                                     pos=I
                                  ENDIF
                               ENDIF
                            ENDIF
                         ENDIF
                      END DO
                   ENDIF
                   
                   !
                   ! If detected peak
                   !
                   IF( pos_2 .gt. 0 .and. n2 .gt. 1000)THEN
                      !
                      ! Main event
                      ! If we have a detection then work downward to 2nd mean
                      ! or an upturn < mean3 - 3*M2_2
                      !
                      !
                      ! Make sure original detection is reset
                      !
                      reset_gain = .TRUE.
                      gain_array3(pos_2-start_pos+1,J) = mean_2 - M2_2
                      outData%solar_contamination_3B(pos_2) = .TRUE.
                      !
                      ! Now move up/down away from pos to track to mean
                      ! or upturn within 3 sigma of mean2
                      !
                      Loop1: DO I=pos_2-1,start_pos,-1
                         IF( .not. solza_ascend(I) )THEN
                            K=I-start_pos+1
                            K2=K-20
                            IF( K2 .lt. 1 )K2 = 1
                            IF( gain_array3(K,J) .gt. 0. )THEN
                               IF( gain_array3(K,J) .lt. mean_2 )THEN
                                  !
                                  ! Only look for upturn in < mean - 3*M2
                                  !
                                  IF( gain_array3(K,J) .gt. mean_2 - 3*M2_2 )THEN
                                     ok2 = .FALSE.
                                     ! Look for upturn but only with valid data
                                     inLoop1: DO KK=K2,start_pos,-1
                                        ! Trap no further data
                                        IF( gain_array3(KK,J) .gt. 0. )THEN
                                           ok2 = .TRUE.
                                           IF( gain_array3(K,J) .gt. gain_array3(KK,J) )THEN
                                              EXIT Loop1
                                           ELSE
                                              gain_array3(K,J) = mean_2 - M2_2
                                              outData%solar_contamination_3B(I) = .TRUE.
                                              EXIT inLoop1
                                           ENDIF
                                        ENDIF
                                     END DO inLoop1
                                     IF( .not. ok2 )THEN
                                        EXIT Loop1
                                     ENDIF
                                  ELSE
                                     gain_array3(K,J) = mean_2 - M2_2
                                     outData%solar_contamination_3B(I) = .TRUE.
                                  ENDIF
                               ELSE
                                  EXIT Loop1
                               ENDIF
                            ENDIF
                         ELSE
                            EXIT Loop1
                         ENDIF
                      END DO Loop1
                      Loop2: DO I=pos_2+1,stop_pos
                         IF( .not. solza_ascend(I) )THEN
                            K=I-start_pos+1
                            K2=K+20
                            IF( K2 .gt. stop_pos )K2 = stop_pos
                            IF( gain_array3(K,J) .gt. 0. )THEN
                               IF( gain_array3(K,J) .lt. mean_2 )THEN
                                  !
                                  ! Only look for upturn in < mean - 4*M2
                                  !
                                  IF( gain_array3(K,J) .gt. mean_2 - 3*M2_2 )THEN
                                     ok2 = .FALSE.
                                     inLoop2: DO KK=K2,stop_pos
                                        IF( gain_array3(KK,J) .gt. 0. )THEN
                                           ok2 = .TRUE.
                                           IF( gain_array3(K,J) .gt. gain_array3(KK,J) )THEN
                                              EXIT Loop2
                                           ELSE
                                              gain_array3(K,J) = mean_2 - M2_2
                                              outData%solar_contamination_3B(I) = .TRUE.
                                              EXIT inLoop2
                                           ENDIF
                                        ENDIF
                                     END DO inLoop2
                                     IF( .not. ok2 )THEN
                                        EXIT Loop2
                                     ENDIF
                                  ELSE
                                     gain_array3(K,J) = mean_2 - M2_2
                                     outData%solar_contamination_3B(I) = .TRUE.
                                  ENDIF
                               ELSE
                                  EXIT Loop2
                               ENDIF
                            ENDIF
                         ELSE
                            EXIT Loop2
                         ENDIF
                      END DO Loop2
                   ELSE IF( pos .gt. 0 .and. n .gt. 1000 )THEN
                      reset_gain = .TRUE.
                      gain_array3(pos-start_pos+1,J) = mean
                      outData%solar_contamination_3B(pos) = .TRUE.
                      !
                      ! Now move up/down away from pos to track to mean
                      ! or upturn within 3 sigma of mean2
                      !
                      Loop1_2: DO I=pos-1,start_pos,-1
                         IF( solza_ascend(I) .and. smooth_solZA(I) .gt. 90. )THEN
                            K=I-start_pos+1
                            K2=K-20
                            IF( K2 .lt. 1 )K2 = 1
                            IF( gain_array3(K,J) .gt. 0. )THEN
                               IF( gain_array3(K,J) .lt. mean )THEN
                                  !
                                  ! Only look for upturn in < mean - 3*M2
                                  !
                                  IF( gain_array3(K,J) .gt. mean - 3*M2 )THEN
                                     ok2 = .FALSE.
                                     ! Look for upturn but only with valid data
                                     inLoop1_2: DO KK=K2,start_pos,-1
                                        ! Trap no further data
                                        IF( gain_array3(KK,J) .gt. 0. )THEN
                                           ok2 = .TRUE.
                                           IF( gain_array3(K,J) .gt. gain_array3(KK,J) )THEN
                                              EXIT Loop1_2
                                           ELSE
                                              gain_array3(K,J) = mean
                                              outData%solar_contamination_3B(I) = .TRUE.
                                              EXIT inLoop1_2
                                           ENDIF
                                        ENDIF
                                     END DO inLoop1_2
                                     IF( .not. ok2 )THEN
                                        EXIT Loop1_2
                                     ENDIF
                                  ELSE
                                     gain_array3(K,J) = mean
                                     outData%solar_contamination_3B(I) = .TRUE.
                                  ENDIF
                               ELSE
                                  EXIT Loop1_2
                               ENDIF
                            ENDIF
                         ELSE
                            EXIT Loop1_2
                         ENDIF
                      END DO Loop1_2
                      Loop2_2: DO I=pos+1,stop_pos
                         IF( solza_ascend(I) .and. smooth_solZA(I) .gt. 90. )THEN
                            K=I-start_pos+1
                            K2=K+20
                            IF( K2 .gt. stop_pos )K2 = stop_pos
                            IF( gain_array3(K,J) .gt. 0. )THEN
                               IF( gain_array3(K,J) .lt. mean )THEN
                                  !
                                  ! Only look for upturn in < mean - 4*M2
                                  !
                                  IF( gain_array3(K,J) .gt. mean - 3*M2 )THEN
                                     ok2 = .FALSE.
                                     inLoop2_2: DO KK=K2,stop_pos
                                        IF( gain_array3(KK,J) .gt. 0. )THEN
                                           ok2 = .TRUE.
                                           IF( gain_array3(K,J) .gt. gain_array3(KK,J) )THEN
                                              EXIT Loop2_2
                                           ELSE
                                              gain_array3(K,J) = mean
                                              outData%solar_contamination_3B(I) = .TRUE.
                                              EXIT inLoop2_2
                                           ENDIF
                                        ENDIF
                                     END DO inLoop2_2
                                     IF( .not. ok2 )THEN
                                        EXIT Loop2_2
                                     ENDIF
                                  ELSE
                                     gain_array3(K,J) = mean 
                                     outData%solar_contamination_3B(I) = .TRUE.
                                  ENDIF
                               ELSE
                                  EXIT Loop2_2
                               ENDIF
                            ENDIF
                         ELSE
                            EXIT Loop2_2
                         ENDIF
                      END DO Loop2_2
                   ELSE
                      ok = .FALSE.
                   ENDIF
                END DO
             ELSE
                outData%solar_contamination_failure = .TRUE.
                print *,'Number of solZA ascending data:',n,n2
                CALL out_Warning(&
                     'There is not enough good data to determine the solar contamination',&
                     'Recalibrate_AVHRR','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
             ENDIF
          END DO
       ENDIF
       DEALLOCATE(solza_ascend)
       DEALLOCATE(smooth_solZA)
    ENDIF
    IF( reset_gain )THEN
       WRITE(*,'(''...Detected solar contamination using Simple code (V1.0)...'')')
    ELSE IF( operational_solar .and. use_old_solar_contam )THEN
       !
       ! Reset gain for operational case of solar contamination detection
       !
       reset_gain = .FALSE.
       DO K=1,ncoefs
          I=start_pos
          ok = .TRUE.
          DO WHILE(ok)
             IF( outData%orig_solar_contamination_3B(I) )THEN
                found=.FALSE.
                solar_loop:DO J=I,stop_pos
                   IF( .not. outData%orig_solar_contamination_3B(J) )THEN
                      found=.TRUE.
                      IF(I .gt. start_pos )THEN
                         mean_gain = (gain_array3(I-start_pos,K)+gain_array3(J-start_pos+1,K))/2.
                      ELSE
                         mean_gain = gain_array3(J-start_pos+1,K)
                      ENDIF
                      reset_gain = .TRUE.
                      DO JJ=I,J-1
                         gain_array3(JJ-start_pos+1,K) = mean_gain
                         outData%solar_contamination_3B(JJ) = .TRUE.
                      END DO
                      I=J
                      EXIT solar_loop
                   ENDIF
                END DO solar_loop
                IF( .not. found )THEN
                   IF(I .gt. start_pos )THEN
                      mean_gain = gain_array3(I-start_pos+1,K)
                   ELSE
                      mean_gain = gain_array3(1,K)
                   ENDIF
                   reset_gain = .TRUE.
                   DO JJ=I,stop_pos
                      gain_array3(JJ-start_pos+1,K) = mean_gain
                      outData%solar_contamination_3B(JJ) = .TRUE.
                   END DO
                   ok = .FALSE.
                ENDIF
             ELSE
                I=I+1
                IF( I .gt. stop_pos )THEN
                   ok = .FALSE.
                ENDIF
             ENDIF
          END DO
       END DO
       IF( reset_gain )THEN
          WRITE(*,'(''...Using original Solar contamination flags with new calibration...'')')
       ENDIF
    ENDIF
    !
    ! if we have reset the gain then recalulate 3.7 micron coeffiecnts
    !
    IF( reset_gain )THEN
       DO J=1,ncoefs
          DO I=start_pos,stop_pos
             IF( new_calib3(1,I,J) .ne. NAN_R )THEN
                K=I-start_pos+1
                new_calib3(1,I,J)=&
                     (coefs1(1,J)+gain_array3(K,J)*outData%smoothSp3(i)+&
                     coefs1(4,J)*outData%smoothSp3(i)*outData%smoothSp3(i)) + &
                     coefs1(8,J)*outData%orbital_temperature
                new_calib3(2,I,J) = -(gain_array3(K,J) + &
                     2*coefs1(4,J)*outData%smoothSp3(i))
                new_calib3(3,I,J) = coefs1(4,J)
             ENDIF
          END DO
       END DO
    ELSE
       WRITE(*,'(''...No Solar contamination found...'')')
    ENDIF

    ! Copy first (prime) calibration coefficients over
    outData%new_calib3(:,:) = new_calib3(:,:,1)
    outData%new_calib4(:,:) = new_calib4(:,:,1)
    outData%new_calib5(:,:) = new_calib5(:,:,1)
   
    ! Now calibrate data
    outData%new_array3B = NAN_R
    outData%new_array4 = NAN_R
    outData%new_array5 = NAN_R
    outData%new_array3B_error = NAN_R
    outData%new_array4_error = NAN_R
    outData%new_array5_error = NAN_R
    DO I=start_pos,stop_pos
       if( noise_all_data )then
          K=1
       ELSE
          K=I
       ENDIF
       if( new_calib3(1,I,1) .ne. NAN_R )then
          DO J=1,outData%nelem
             counts = outData%counts3(J,I)
             if ( counts .ne. NAN_R ) then
                radiance = new_calib3(1,I,1) + &
                     new_calib3(2,I,1) * counts + &
                     new_calib3(3,I,1) * counts*counts
                if ((radiance .lt. 0) .and. (counts .ne. 992)) then                         
                !     print*, "r neg", radiance, counts, new_calib3(1,I,1),new_calib3(2,I,1), new_calib3(3,I,1) 
                end if    
                ! Interpolate two calibrations
                IF( 2 .eq. ncoefs )THEN
                   radiance2 = new_calib3(1,I,2) + &
                        new_calib3(2,I,2) * counts + &
                        new_calib3(3,I,2) * counts*counts
                   radiance = coefs_frac1*radiance + coefs_frac2*radiance2
                ENDIF
                IF( outData%interpolated(I) )THEN
                   IF( 0 .lt. noise(1,K) )THEN
                      radiance_error = sqrt((new_calib3(2,I,1) + &
                           2*new_calib3(3,I,1)*counts)**2*(noise(4,K)**2) + &
                           (stdev3 * counts)**2)
                   ELSE
                      radiance_error=NAN_R
                   ENDIF
                ELSE
                   IF( 0 .lt. noise(1,K) )THEN
                      radiance_error = (new_calib3(2,I,1) + &
                           2*new_calib3(3,I,1)*counts)*noise(4,K)
                   ELSE
                      radiance_error=NAN_R
                   ENDIF
                ENDIF
                IF( out_radiances )THEN
                   !print*,"rad"
                   outData%new_array3B(J,I) = radiance
                   outData%new_array3b_error(J,I) = radiance_error
                ELSE
                   outData%new_array3B(J,I) = convertBT(radiance,coefs1(5,1),&
                        coefs1(6,1),coefs1(7,1))
#ifdef USE_IEEE
                   IF( .NOT.ieee_is_finite(outData%new_array3B(J,I)) )THEN
                      outData%new_array3B(J,I) = NAN_R
                   END IF
#elif defined(USE_G95)
                   IF( isNAN(outData%new_array3B(J,I)) )THEN
                      outData%new_array3B(J,I) = NAN_R
                   END IF
#endif
                   IF( NAN_R .ne. radiance_error )THEN
                      if( 0. .lt. radiance-radiance_error )then
                         error1 = ABS(outData%new_array3b(J,I) - &
                              convertBT(radiance-radiance_error,coefs1(5,1),&
                              coefs1(6,1),coefs1(7,1)))
                      else
                         error1 = 10
                      endif
                      error2 = ABS(outData%new_array3b(J,I) - &
                           convertBT(radiance+radiance_error,coefs1(5,1),&
                           coefs1(6,1),coefs1(7,1)))
                      outData%new_array3B_error(J,I) = MAX(error1,error2)
#ifdef USE_IEEE
                      IF( .NOT.ieee_is_finite(outData%new_array3B_error(J,I)) )THEN
                         outData%new_array3B_error(J,I) = NAN_R
                      END IF
#elif defined(USE_G95)
                      IF( isNAN(outData%new_array3B_error(J,I)) )THEN
                         outData%new_array3B_error(J,I) = NAN_R
                      END IF
#endif
                   ELSE
                      outData%new_array3B_error(J,I) = NAN_R
                   ENDIF
                ENDIF
                IF( outData%isGAC )THEN
                   !
                   ! For GAC average nedt - but only if we're in the Gaussian
                   ! regime (counts > 3). If not use scale (default 1.)
                   !
                   IF( noise(4,K) .ge. 3 )THEN
                      !
                      ! Gaussian regime
                      !
                      outData%new_array3b_error(J,I) = &
                           outData%new_array3b_error(J,I)/2.
                   ELSE
                      outData%new_array3b_error(J,I) = &
                           outData%new_array3b_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             else
                outData%new_array3B(J,I) = NAN_R
             endif
               
          END DO
       endif
       if( new_calib4(1,I,1) .ne. NAN_R )then
          DO J=1,outData%nelem
             counts = outData%counts4(J,I)
             if ( counts .ne. NAN_R ) then
                radiance = new_calib4(1,I,1) + &
                     new_calib4(2,I,1) * counts + &
                     new_calib4(3,I,1) * counts*counts
                ! Interpolate two calibrations
                IF( 2 .eq. ncoefs )THEN
                   radiance2 = new_calib4(1,I,2) + &
                        new_calib4(2,I,2) * counts + &
                        new_calib4(3,I,2) * counts*counts
                   radiance = coefs_frac1*radiance + coefs_frac2*radiance2
                ENDIF
                IF( outData%interpolated(I) )THEN
                   IF( 0 .lt. noise(2,K) )THEN
                      radiance_error = sqrt((new_calib4(2,I,1) + &
                        2*new_calib4(3,I,1)*counts)**2*(noise(5,K)**2) + &
                        (stdev4 * counts)**2)
                   ELSE
                      radiance_error = NAN_R
                   ENDIF
                ELSE
                   IF( 0 .lt. noise(5,K) )THEN
                      radiance_error = (new_calib4(2,I,1) + &
                           2*new_calib4(3,I,1)*counts)*noise(5,K)
                   ELSE
                      radiance_error = NAN_R
                   ENDIF
                ENDIF
                IF( out_radiances )THEN
                   outData%new_array4(J,I) = radiance
                   outData%new_array4_error(J,I) = radiance_error
                ELSE
                   outData%new_array4(J,I) = convertBT(radiance,coefs2(5,1),&
                        coefs2(6,1),coefs2(7,1))
#ifdef USE_IEEE
                   IF( .NOT.ieee_is_finite(outData%new_array4(J,I)) )THEN
                      outData%new_array4(J,I) = NAN_R
                   END IF
#elif defined(USE_G95)
                   IF( isNAN(outData%new_array4(J,I)) )THEN
                      outData%new_array4(J,I) = NAN_R
                   END IF
#endif
                   IF( NAN_R .ne. radiance_error )THEN
                      if( 0. .lt. radiance-radiance_error )then
                         error1 = ABS(outData%new_array4(J,I) - &
                              convertBT(radiance-radiance_error,coefs2(5,1),&
                              coefs2(6,1),coefs2(7,1)))
                      else
                         error1 = 10
                      endif
                      error2 = ABS(outData%new_array4(J,I) - &
                           convertBT(radiance+radiance_error,coefs2(5,1),&
                           coefs2(6,1),coefs2(7,1)))
                      outData%new_array4_error(J,I) = MAX(error1,error2)
#ifdef USE_IEEE
                      IF( .NOT.ieee_is_finite(outData%new_array4_error(J,I)) )THEN
                         outData%new_array4_error(J,I) = NAN_R
                      END IF
#elif defined(USE_G95)
                      IF( isNAN(outData%new_array4_error(J,I)) )THEN
                         outData%new_array4_error(J,I) = NAN_R
                      END IF
#endif
                   ELSE
                      outData%new_array4_error(J,I) = NAN_R
                   ENDIF
                endif
                IF( outData%isGAC )THEN
                   !
                   ! For GAC average nedt - but only if we're in the Gaussian
                   ! regime (counts > 3). If not use scale (default 1.)
                   !
                   IF( noise(5,K) .ge. 3 )THEN
                      outData%new_array4_error(J,I) = &
                           outData%new_array4_error(J,I)/2.
                   ELSE
                      outData%new_array4_error(J,I) = &
                           outData%new_array4_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             else
                outData%new_array4(J,I) = NAN_R
             endif
          END DO
       endif
       if( new_calib5(1,I,1) .ne. NAN_R )then
          DO J=1,outData%nelem
             counts = outData%counts5(J,I)
             if ( counts .ne. NAN_R ) then
                radiance = new_calib5(1,I,1) + &
                     new_calib5(2,I,1) * counts + &
                     new_calib5(3,I,1) * counts*counts
                ! Interpolate two calibrations
                IF( 2 .eq. ncoefs )THEN
                   radiance2 = new_calib5(1,I,2) + &
                        new_calib5(2,I,2) * counts + &
                        new_calib5(3,I,2) * counts*counts
                   radiance = coefs_frac1*radiance + coefs_frac2*radiance2
                ENDIF
                IF( outData%interpolated(I) )THEN
                   IF( 0 .lt. noise(3,K) )THEN
                      radiance_error = sqrt((new_calib5(2,I,1) + &
                           2*new_calib5(3,I,1)*counts)**2*(noise(6,K))**2 + &
                           (stdev5 * counts)**2)
                   ELSE
                      radiance_error = NAN_R
                   ENDIF
                ELSE
                   IF( 0 .lt. noise(6,K) )THEN
                      radiance_error = (new_calib5(2,I,1) + &
                           2*new_calib5(3,I,1)*counts)*noise(6,K)
                   ELSE
                      radiance_error = NAN_R
                   ENDIF
                ENDIF
                IF( out_radiances )THEN
                   outData%new_array5(J,I) = radiance
                   outData%new_array5_error(J,I) = radiance_error
                ELSE
                   outData%new_array5(J,I) = convertBT(radiance,coefs3(5,1),&
                        coefs3(6,1),coefs3(7,1))
#ifdef USE_IEEE
                   IF( .NOT.ieee_is_finite(outData%new_array5(J,I)) )THEN
                      outData%new_array5(J,I) = NAN_R
                   END IF
#elif defined(USE_G95)
                   IF( isNAN(outData%new_array5(J,I)) )THEN
                      outData%new_array5(J,I) = NAN_R
                   END IF
#endif
                   IF( NAN_R .ne. radiance_error )THEN
                      if( 0. .lt. radiance-radiance_error )then
                         error1 = ABS(outData%new_array5(J,I) - &
                              convertBT(radiance-radiance_error,coefs3(5,1),&
                              coefs3(6,1),coefs3(7,1)))
                      else
                         error1 = 10
                      endif
                      error2 = ABS(outData%new_array5(J,I) - &
                           convertBT(radiance+radiance_error,coefs3(5,1),&
                           coefs3(6,1),coefs3(7,1)))
                      outData%new_array5_error(J,I) = MAX(error1,error2)
#ifdef USE_IEEE
                      IF( .NOT.ieee_is_finite(outData%new_array5_error(J,I)) )THEN
                         outData%new_array5_error(J,I) = NAN_R
                      END IF
#elif defined(USE_G95)
                      IF( isNAN(outData%new_array5_error(J,I)) )THEN
                         outData%new_array5_error(J,I) = NAN_R
                      END IF
#endif
                   ELSE
                      outData%new_array5_error(J,I) = NAN_R
                   ENDIF
                ENDIF
                IF( outData%isGAC )THEN
                   !
                   ! For GAC average nedt - but only if we're in the Gaussian
                   ! regime (counts > 3). If not use scale (default 1.)
                   !
                   IF( noise(6,K) .ge. 3 )THEN
                      outData%new_array5_error(J,I) = &
                           outData%new_array5_error(J,I)/2.
                   ELSE
                      outData%new_array5_error(J,I) = &
                           outData%new_array5_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             else
                outData%new_array5(J,I) = NAN_R
             endif
          END DO
       endif

    END DO

    ! Do new visible channel (CLAVR-X) if requested
    IF( new_vis_calibration )THEN
       CALL Update_Vischan( outData, start_pos, stop_pos, noise, nedt_scale, &
            noise_all_data )
    ENDIF

    ! Set non-calibrated data to bad
    DO I=1,start_pos-1
       outData%badCalibration(I) = .TRUE.
       call set_calib_bad( outData, I )
    END DO

    DO I=stop_pos+1,outData%arraySize
       outData%badCalibration(I) = .TRUE.
       call set_calib_bad( outData, I )
    END DO

    DEALLOCATE( gain_array3,gain_array4, &
         gain_array5,rad_array3,rad_array4, &
         rad_array5,good_data,solar_data,&
         rad_array3_correl,rad_array4_correl,&
         rad_array5_correl,new_calib3,new_calib4,&
         new_calib5,noise,&
         STAT=STAT )

  END SUBROUTINE Recalibrate_AVHRR

  SUBROUTINE Find_Moon_Events( start_pos, stop_pos, outData )

    INTEGER, INTENT(IN) :: start_pos
    INTEGER, INTENT(IN) :: stop_pos
    TYPE(AVHRR_Data), INTENT(INOUT) :: outData

    ! Local variables
    INTEGER :: STAT
    INTEGER(GbcsInt1), ALLOCATABLE :: OK37(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: OK11(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: OK12(:)
    INTEGER :: I,J,K
    INTEGER :: start, end
    INTEGER :: minpos, maxpos
    INTEGER :: width
    LOGICAL :: printOutput
    LOGICAL :: test_moon

    ! Allocate OK array
    ALLOCATE(OK37(outData%arraySize),&
         OK11(outData%arraySize),&
         OK12(outData%arraySize),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot allocate OK array','Find_Moon_Events',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    OK37 = 0
    OK37(start_pos:stop_pos) = 1
    OK11 = 0
    OK11(start_pos:stop_pos) = 1
    OK12 = 0
    OK12(start_pos:stop_pos) = 1
    outData%moon_contamination = .FALSE.

    ! Note don't use smoothed versions for this

    ! 3.7 micron channel
    test_moon = .TRUE.
    IF( 1 .eq. Filter_Moon_Events(start_pos,stop_pos,OK37,outData%sp3,3.,&
         test_moon) )THEN
       RETURN
    ENDIF
    STAT = Filter_Moon_Events(start_pos,stop_pos,OK37,outData%sp3,3.,&
         test_moon)

    ! 11 micron channel
    test_moon = .TRUE.
    IF( 1 .eq. Filter_Moon_Events(start_pos,stop_pos,OK11,outData%sp4,3.,&
         test_moon) )THEN
       RETURN
    ENDIF
    STAT = Filter_Moon_Events(start_pos,stop_pos,OK11,outData%sp4,3.,&
         test_moon)

    ! 12 micron channel
    test_moon = .TRUE.
    IF( 1 .eq. Filter_Moon_Events(start_pos,stop_pos,OK12,outData%sp3,3.,&
         test_moon) )THEN
       RETURN
    ENDIF
    STAT = Filter_Moon_Events(start_pos,stop_pos,OK12,outData%sp3,3.,&
         test_moon)

    ! Add width of 10%
    start = -1
    end = -1
    I = start_pos
    DO WHILE(I .le. stop_pos)
       IF( OK37(I) .eq. 0 .or. OK11(I) .eq. 0 .or. OK12(I) .eq. 0)THEN
          start = I
          Loop: DO J=I,stop_pos
             IF( OK37(J) .eq. 1 .and. OK11(J) .eq. 1 .and. OK12(J) .eq. 1)THEN
                end = J-1
                width = end-start
                IF( 10 .lt. width )THEN
                   minpos = MAX(1,I-width/10)
                   maxpos = MIN(stop_pos,J+width/10)
                   DO K=minpos,maxpos
                      OK37(K) = 0
                      OK11(K) = 0
                      OK12(K) = 0
                   END DO
                   I=maxpos-1
                ELSE
                   I=J
                ENDIF
                EXIT Loop
             ENDIF
          END DO LOOP
       END IF
       I=I+1
    END DO

    ! Filter out moon events
    printOutput = .TRUE.
    DO I=start_pos,stop_pos
       IF( 0 .eq. OK37(I) .or. 0 .eq. OK11(I) .or. 0 .eq. OK12(I) )THEN
          outData%moon_contamination(I) = .TRUE.
          outData%smoothPrt(I) = NAN_R
          outData%smoothsp3(I) = NAN_R
          outData%smoothsp4(I) = NAN_R
          outData%smoothsp5(I) = NAN_R
          outData%smoothbb3(I) = NAN_R
          outData%smoothbb4(I) = NAN_R
          outData%smoothbb5(I) = NAN_R
          IF( printOutput )THEN
             call out_WARNING('Moon events detected','Find_Moon_Events',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
             printOutput = .FALSE.
          ENDIF
       ENDIF
    END DO

    ! Deallocate OK array
    DEALLOCATE(OK37,OK11,OK12)

  END SUBROUTINE Find_Moon_Events

  INTEGER FUNCTION Filter_Moon_Events(start_pos, stop_pos, OK, array, &
       limit, test_moon ) RESULT(STATUS)

    INTEGER, INTENT(IN) :: start_pos
    INTEGER, INTENT(IN) :: stop_pos
    INTEGER(GbcsInt1), INTENT(INOUT) :: OK(:)
    REAL, INTENT(IN) :: array(:)
    REAL, INTENT(IN) :: limit
    LOGICAL, INTENT(INOUT) :: test_moon

    ! Local variables
    INTEGER :: I
    INTEGER :: ndata
    REAL :: maxVal
    REAL :: meanVal
    REAL :: stdevVal

    STATUS=0
    ndata = 0
    meanVal = 0.
    maxVal = 0.
    stdevVal = 0.
    DO I=start_pos,stop_pos
       IF( 1 .eq. OK(I) )THEN
          IF( 900 .lt. array(I) .and. 1000 .gt. array(I) )THEN
             ndata = ndata + 1
             meanVal = meanVal + array(I)
             IF( array(I) .gt. maxVal )THEN
                maxVal = array(I)
             ENDIF
          ENDIF
       ENDIF
    END DO
    IF( 100 .gt. ndata )THEN
       CALL out_WARNING('Not enough data to find moon events',&
            'Filter_Moon_Events','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       STATUS=1
       RETURN
    ENDIF
    meanVal = meanVal / ndata

    ! Coarse check to make sure we don't over do it...
    IF( 2 .gt. maxVal-meanVal .and. test_moon )THEN
       ! No moon event
       STATUS=1
       RETURN
    ENDIF
    test_moon = .FALSE.

    DO I=start_pos,stop_pos
       IF( 1 .eq. OK(I) )THEN
          IF( 900 .lt. array(I) .and. 1000 .gt. array(I) )THEN
             stdevVal = stdevVal + (array(I)-meanVal)**2
          ENDIF
       ENDIF
    END DO
    stdevVal = sqrt(stdevVal/(ndata-1))

    DO I=start_pos,stop_pos
       IF( meanVal+limit*stdevVal .lt. array(I) )THEN
          OK(I) = 0
       ENDIF
    END DO

  END FUNCTION Filter_Moon_Events

  SUBROUTINE Earthshine_Correct( outData, coefs, ndata, good_data, &
       start_pos, stop_pos, night_ok, gain_array )

    TYPE(AVHRR_Data), INTENT(INOUT) :: outData
    REAL, INTENT(IN) :: coefs(7)
    INTEGER, INTENT(IN) :: ndata
    LOGICAL, INTENT(IN) :: good_data(ndata)
    INTEGER, INTENT(IN) :: start_pos, stop_pos
    LOGICAL, INTENT(IN) :: night_ok
    REAL, INTENT(INOUT) :: gain_array(ndata)

    ! Local variables
    INTEGER :: I,J
    INTEGER :: STAT
    REAL, ALLOCATABLE :: radiance(:)
    LOGICAL, ALLOCATABLE :: good(:)
    REAL, ALLOCATABLE :: gain(:)
    INTEGER :: ndata_new
    REAL :: eta
    INTEGER :: pos

    ALLOCATE( radiance(ndata), good(ndata), gain(ndata), STAT=STAT )
    CALL Check_IOS(STAT,'Cannot allocate radiance/good array',&
         'Earthshine_Correct','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    ! Make arrays for fitting eta to
    IF( night_ok )THEN
       pos = outData%nelem/2
       good = .FALSE.
       J=0
       ndata_new = 0
       DO I=start_pos,stop_pos
          ! Check if this is data we want to derive model with
          if( NAN_R .ne. outData%smoothPrt(I) )then
             if( NAN_R .ne. outData%smoothBB3(I) .and. &
                  NAN_R .ne. outData%smoothSp3(I) .and. &
                  0 .lt. outData%smoothBB3(I) .and. &
                  0 .lt. outData%smoothSp3(I) )then
                ndata_new=ndata_new+1
                if( 90 .le. outData%solZA(pos,I) )then
                   good(ndata_new) = .TRUE.
                endif
             endif
          endif
       END DO
       if( ndata_new .ne. ndata )then
         call out_ERROR('ndata_new != ndata','Earthshine_Correct',&
              'NOAA_LoadAVHRRLevel1B.f90',.TRUE.,' ')
      endif
    ELSE
       ! Only ignore solar contaminated bits
       good = good_data
    ENDIF

    ! Get initial value of radiance
    CALL Get_gain_with_eta(  outData, coefs, ndata, ndata_new, start_pos, &
         stop_pos, gain_array, good_data, radiance, 0., get_radiance=.TRUE. )

    CALL Newton_Rapson(outData, coefs, ndata, start_pos, stop_pos, &
         gain_array, good_data, radiance, eta)

    ! Get final answer on input good_data
    write(*,'(''    Correcting for Eartshine'')')
    CALL Get_gain_with_eta(  outData, coefs, ndata, ndata_new, start_pos, &
         stop_pos, gain_array, good_data, radiance, eta )
    outData%earthshine_eta = eta
    write(*,'(''    eta = '',e12.5)')eta

    DEALLOCATE(radiance,good,gain)

  END SUBROUTINE Earthshine_Correct

  SUBROUTINE Newton_Rapson(outData, coefs, ndata, start_pos, stop_pos, &
       gain_array, good_data, radiance, eta)

    TYPE(AVHRR_Data), INTENT(INOUT) :: outData
    REAL, INTENT(IN) :: coefs(7)
    INTEGER, INTENT(IN) :: ndata
    INTEGER, INTENT(IN) :: start_pos, stop_pos
    REAL, INTENT(INOUT) :: gain_array(ndata)
    LOGICAL, INTENT(IN) :: good_data(ndata)
    REAL, INTENT(INOUT) :: radiance(ndata)
    REAL, INTENT(OUT) :: eta

    INTEGER :: I
    REAL :: previous_eta
    REAL, PARAMETER :: TOL = 1e-5
    LOGICAL :: start
    INTEGER, PARAMETER :: Max_Iter=20
    REAL :: correl
    REAL :: derivative
    REAL :: dx
    LOGICAL :: ok

    eta = 0.005
    start = .TRUE.

    ok = .FALSE.
    Loop: DO I=1,Max_Iter
       call Get_Function_Deriv( outData, coefs, ndata, start_pos, stop_pos, &
            gain_array, good_data, radiance, eta, correl, derivative )
       dx = correl/derivative
       eta = eta-dx
       if( start )then
          start = .FALSE.
       else
          if( ABS(dx/eta) .lt. TOL )then
             ok = .TRUE.
             EXIT Loop
          end if
       endif
    END DO Loop

    IF( .not.ok )THEN
       CALL out_WARNING('Exceeded number of iterations in Newton-Rapson',&
            'Newton_Rapson','NOAA_LoadAVHRRLevel1B.f90',.TRUE.,' ')
       eta = 0.
    ENDIF

  END SUBROUTINE Newton_Rapson

  SUBROUTINE Get_Function_Deriv( outData, coefs, ndata, start_pos, stop_pos, &
       gain_array, good_data, radiance, eta, correl, derivative )

    TYPE(AVHRR_Data), INTENT(INOUT) :: outData
    REAL, INTENT(IN) :: coefs(7)
    INTEGER, INTENT(IN) :: ndata
    INTEGER, INTENT(IN) :: start_pos, stop_pos
    LOGICAL, INTENT(IN) :: good_data(ndata)
    REAL, INTENT(INOUT) :: gain_array(ndata)
    REAL, INTENT(INOUT) :: radiance(ndata)
    REAL, INTENT(IN) :: eta
    REAL, INTENT(OUT) :: correl
    REAL, INTENT(OUT) :: derivative

    REAL :: delta_eta
    REAL :: correl2

    INTEGER :: I
    INTEGER :: ndata_new

    CALL Get_gain_with_eta(  outData, coefs, ndata, ndata_new, start_pos, &
         stop_pos, gain_array, good_data, radiance, eta, get_radiance=.FALSE. )
    correl = Get_Correlation( ndata_new, radiance, gain_array )

    DO I=1,ndata_new
       WRITE(101,*)radiance(I),gain_array(I)
    END DO
    STOP 1

    delta_eta = 0.0001

    CALL Get_gain_with_eta(  outData, coefs, ndata, ndata_new, start_pos, &
         stop_pos, gain_array, good_data, radiance, eta+delta_eta, &
         get_radiance=.TRUE. )
    correl2 = Get_Correlation( ndata_new, radiance, gain_array )

    derivative = (correl2-correl)/delta_eta

  END SUBROUTINE Get_Function_Deriv

  REAL FUNCTION Get_Correlation( ndata, X, Y )RESULT(correl)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: X(ndata)
    REAL, INTENT(IN) :: Y(ndata)

    INTEGER :: I
    REAL :: xmean, ymean
    REAL :: xsq, ysq, total

    xmean = 0.
    ymean = 0.
    DO I=1,ndata
       xmean = xmean + X(I)
       ymean = ymean + Y(I)
    END DO
    xmean = xmean / ndata
    ymean = ymean / ndata

    total = 0.
    xsq = 0.
    ysq = 0.
    DO I=1,ndata
       total = total + (X(I)-xmean)*(Y(I)-ymean)
       xsq = xsq + (X(I)-xmean)*(X(I)-xmean)
       ysq = ysq + (Y(I)-ymean)*(Y(I)-ymean)
    END DO
    correl = total/(sqrt(xsq)*sqrt(ysq))

  END FUNCTION Get_Correlation

  SUBROUTINE Get_gain_with_eta(  outData, coefs, ndata, ndata_new, &
       start_pos, stop_pos, gain_array, good_data, radiance, eta, &
       get_radiance )

    TYPE(AVHRR_Data), INTENT(INOUT) :: outData
    REAL, INTENT(IN) :: coefs(7)
    INTEGER, INTENT(IN) :: ndata
    INTEGER, INTENT(OUT) :: ndata_new
    INTEGER, INTENT(IN) :: start_pos, stop_pos
    REAL, INTENT(INOUT) :: gain_array(ndata)
    LOGICAL, INTENT(IN) :: good_data(ndata)
    REAL, INTENT(INOUT) :: radiance(ndata)
    REAL, INTENT(IN) :: eta
    LOGICAL, OPTIONAL :: get_radiance

    INTEGER :: I, J, J3
    REAL :: R_ICT
    REAL :: counts
    LOGICAL :: return_radiance
    REAL :: new_calib3_1
    REAL :: new_calib3_2
    REAL :: new_calib3_3
    INTEGER :: ndata_line

    return_radiance = .FALSE.
    IF( PRESENT(get_radiance) )THEN
       return_radiance = get_radiance
    ENDIF

    J3 = 0
    DO I=start_pos,stop_pos
       ! Check if this is data we want to derive model with
       if( NAN_R .ne. outData%smoothPrt(I) .and. &
            good_data(I-start_pos+1) )then
          if( NAN_R .ne. outData%smoothBB3(I) .and. &
               NAN_R .ne. outData%smoothSp3(I) .and. &
               0 .lt. outData%smoothBB3(I) .and. &
               0 .lt. outData%smoothSp3(I) )then
             R_ICT = convertRadiance( outData%smoothPrt(i),coefs(4),&
                  coefs(5),coefs(6))
             counts = (outData%smoothSp3(I) - outData%smoothBB3(I))
             J3 = J3+1
             gain_array(J3) = (eta*radiance(J3) + &
                  (eta_ict+coefs(2))*R_ICT - coefs(3)*counts*counts)/&
                  counts
             IF( return_radiance )THEN
                ! Convert to level 1B type coefficients
                new_calib3_1 = coefs(1)
                new_calib3_2 = gain_array(J3)
                new_calib3_3 = coefs(3)
                ndata_line = 0
                radiance(J3) = 0.
                ! Get average radiance across scan line
                DO J=1,outData%nelem
                   IF( NAN_R .ne. outData%counts3(J,i) )THEN
                      radiance(J3) = radiance(J3) + new_calib3_1 + &
                           new_calib3_2 * &
                           (outData%smoothSp3(i) - outData%counts3(J,i)) + &
                           new_calib3_3 * &
                           (outData%smoothSp3(i) - outData%counts3(J,i))**2
                      ndata_line = ndata_line + 1
                   ENDIF
                END DO
                radiance(J3) = radiance(J3) / ndata_line
             endif
          endif
       endif
    END DO
    IF( J3 .gt. ndata )THEN
       CALL out_ERROR('J3 > ndata','Get_gain_with_eta',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    ndata_new = J3

  END SUBROUTINE Get_gain_with_eta

  INTEGER FUNCTION Get_Poly_Coefs( ndata, rad_array, gain_array, good_data, &
       stdev, maxdev, Coefs )

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: rad_array(ndata)
    REAL, INTENT(IN) :: gain_array(ndata)
    LOGICAL, INTENT(IN) :: good_data(ndata)
    REAL, INTENT(OUT) :: stdev
    REAL, INTENT(OUT) :: maxdev
    REAL, INTENT(OUT) :: Coefs(3)

    REAL, ALLOCATABLE :: prt_fit(:)
    REAL, ALLOCATABLE :: gain_fit(:)
    REAL, ALLOCATABLE :: gainval(:)

    INTEGER :: STAT
    INTEGER :: I, J
    REAL :: meanval, signval, value

    ALLOCATE(prt_fit(ndata),gain_fit(ndata),gainval(ndata),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating prt_fit,gain_fit',&
         'Get_Poly_Coefs','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    Get_Poly_Coefs = 0

    J = 0
    meanval = 0.
    DO I=1,ndata
       IF( good_data(I) )THEN
          J=J+1
          prt_fit(J) = rad_array(I)
          gain_fit(J) = gain_array(I)*gain_array(I)
          meanval = meanval + gain_array(I)
       ENDIF
    END DO
    IF( 1 .ge. J )THEN
       Get_Poly_Coefs = 1
       CALL out_WARNING('No good data',&
         'Get_Poly_Coefs','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       RETURN
    ENDIF
    meanval = meanval / J
    IF( 0 .gt. meanval )then
       signval = -1
    ELSE
       signval = 1
    ENDIF

    STAT = Poly_Fit(J,prt_fit,gain_fit,3,Coefs)
    IF( 0 .eq. STAT )THEN
       stdev = 0.
       maxdev = 0.
       DO I=1,J
          gainval(I) = Coefs(1)+&
               Coefs(2)*prt_fit(I)+&
               Coefs(3)*prt_fit(I)*prt_fit(I)
          gainval(I) = signval * sqrt(gainval(I))
          value = (sqrt(gain_fit(I)) - gainval(I))
          maxdev = MAX(maxdev,ABS(value))
          stdev = stdev + value*value
       END DO
       stdev = sqrt(stdev/(J-1))
    ELSE
       Get_Poly_Coefs = 1
    ENDIF

    DEALLOCATE(prt_fit,gain_fit,gainval)

  END FUNCTION Get_Poly_Coefs

  ! Get times of solar contamination
  ! 2 step process - first use large range and get function
  ! use gain function to determine when actually to exclude data
  ! Also checks against the level1B stored contamination flag
  INTEGER FUNCTION Remove_Solar_Entries( ndata, solar_za, gain, prt, &
       solar_flag, good_values, night_ok )RESULT(ReturnStat)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: solar_za(ndata)
    REAL, INTENT(IN) :: gain(ndata)
    REAL, INTENT(IN) :: prt(ndata)
    LOGICAL, INTENT(IN) :: solar_flag(ndata)
    LOGICAL, INTENT(INOUT) :: good_values(ndata)
    LOGICAL, INTENT(INOUT) :: night_ok

    ! Parameters to start - solar ZA determined by NOAA-14 data in 2000
    REAL, PARAMETER :: start1=110., stop1=121.
    REAL, PARAMETER :: start2=50., stop2=121.
    REAL, PARAMETER :: sigma_thresh = 2.
    REAL, PARAMETER :: sigma_thresh_first = 4.
    REAL, PARAMETER :: sigma_thresh_large = 10.
    INTEGER :: start1_pos, stop1_pos
    INTEGER :: start2_pos, stop2_pos
    INTEGER :: start3_pos, stop3_pos
    INTEGER :: start1_pos_final, stop1_pos_final
    INTEGER :: start2_pos_final, stop2_pos_final
    INTEGER :: start3_pos_final, stop3_pos_final
    ! Local variables
    INTEGER :: I
    INTEGER :: minPos
    INTEGER :: maxPos
    REAL :: minValue
    REAL :: maxValue

    INTEGER :: new_ndata
    REAL, ALLOCATABLE :: good_prt(:), good_gain(:), good_gain_sq(:)
    LOGICAL, ALLOCATABLE :: bad(:), new_good_values(:)
    REAL, ALLOCATABLE :: gain_value(:)

    INTEGER, PARAMETER :: ncoefs = 3
    REAL :: coefs(ncoefs)
    REAL :: meangain
    REAL :: meanval, value, stdev, maxdev
    REAL :: min_night, max_night, min_day, max_day

    INTEGER :: STAT

    INTEGER :: bval, gval
    LOGICAL :: OK

    ALLOCATE(good_prt(ndata),good_gain(ndata),good_gain_sq(ndata),bad(ndata),&
         new_good_values(ndata),gain_value(ndata),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       call out_ERROR('Cannot allocate good_prt,good_gain',&
            'Remove_Solar_Entries','NOAA_LoadAVHRRLevel1B,f90',.FALSE.,' ')
    ENDIF

    ! First determine if we are on an up or down slope in solarZA
    ! determined by the position of the min and max values
    minPos = 1
    maxPos = 1
    minValue = 1e20
    maxValue = -1e20
    DO I=1,ndata
       IF( NAN_R .ne. solar_za(I) )THEN
          IF( solar_za(I) .lt. minValue )THEN
             minValue = solar_za(I)
             minPos = I
          ELSE IF( solar_za(I) .ge. maxValue )THEN
             maxValue = solar_za(I)
             maxPos = I
          ENDIF
       ENDIF
    END DO

    start1_pos = ndata
    stop1_pos = 0
    start2_pos = ndata
    stop2_pos = 0
    start3_pos = ndata
    stop3_pos = 0
    new_ndata = 0
    meanval = 0.
    bad = .FALSE.
    IF( maxPos .lt. minPos )THEN
       ! Maximum comes first
       DO I=1,maxPos
          IF( solar_za(I) .ge. start1 .and. solar_za(i) .le. stop1 )THEN
             bad(I) = .TRUE.
             start1_pos = MIN(start1_pos,I)
             stop1_pos = MAX(stop1_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
       DO I=maxPos+1,minPos
          IF( solar_za(I) .ge. start2 .and. solar_za(i) .le. stop2 )THEN
             bad(I) = .TRUE.
             start2_pos = MIN(start2_pos,I)
             stop2_pos = MAX(stop2_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
       DO I=minPos+1,ndata
          IF( solar_za(I) .ge. start1 .and. solar_za(i) .le. stop1 )THEN
             bad(I) = .TRUE.
             start3_pos = MIN(start3_pos,I)
             stop3_pos = MAX(stop3_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
    ELSE
       ! Minimum comes first
       DO I=1,minPos
          IF( solar_za(I) .ge. start2 .and. solar_za(i) .le. stop2 )THEN
             bad(I) = .TRUE.
             start1_pos = MIN(start1_pos,I)
             stop1_pos = MAX(stop1_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
       DO I=minPos+1,maxPos
          IF( solar_za(I) .ge. start1 .and. solar_za(i) .le. stop1 )THEN
             bad(I) = .TRUE.
             start2_pos = MIN(start2_pos,I)
             stop2_pos = MAX(stop2_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
       DO I=maxPos+1,ndata
          IF( solar_za(I) .ge. start2 .and. solar_za(i) .le. stop2 )THEN
             bad(I) = .TRUE.
             start3_pos = MIN(start3_pos,I)
             stop3_pos = MAX(stop3_pos,I)
          ELSE
             IF( NAN_R .ne. prt(I) .or. NAN_R .ne. gain(I) )THEN
                new_ndata = new_ndata + 1
                good_prt(new_ndata) = prt(I)
                good_gain(new_ndata) = gain(I)
                good_gain_sq(new_ndata) = gain(I)*gain(I)
                meanval = meanval + gain(I)
             ENDIF
          ENDIF
       END DO
    ENDIF

    IF( 1 .ge. new_ndata )THEN
       call out_ERROR('No good data for solar contamination',&
            'Remove_Solar_Entries','NOAA_LoadAVHRRLevel1B,f90',.FALSE.,' ')
    ENDIF
    meangain = meanval / new_ndata

    ! Now fit quadratic to prt vs gain^2 to get first model
    ! this is used to determine where the bad data really are
    stat = Poly_Fit( new_ndata, good_prt, good_gain_sq, ncoefs, coefs )

    ! Create new gain values and get standard deviation
    ReturnStat = 0
    IF( 0 .eq. stat )THEN
       stdev = 0.
       DO I=1,new_ndata
          gain_value(I) = coefs(1) + coefs(2)*good_prt(I) + &
               coefs(3)*good_prt(I)*good_prt(I)
          gain_value(I) = sqrt(gain_value(I))
          if( 0 .gt. meangain )then
             gain_value(I) = -gain_value(I)
          endif
          value = (good_gain(I)-gain_value(I))
          stdev = stdev + value*value
       END DO
       stdev = sqrt(stdev/(new_ndata-1))

       ! Only look for solar contamination if at least one point is at > 10
       ! sigma - if not then exit
       OK = .FALSE.
       Loop: DO I=1,ndata
          IF( bad(I) )THEN
             IF( NAN_R .ne. gain(I) )THEN
                gain_value(I) = coefs(1) + coefs(2)*prt(I) + &
                     coefs(3)*prt(I)*prt(I)
                gain_value(I) = sqrt(gain_value(I))
                if( 0 .gt. meangain )then
                   gain_value(I) = -gain_value(I)
                endif
                IF( ABS(gain_value(I)-gain(I)) .ge. sigma_thresh_large*stdev )THEN
                   OK = .TRUE.
                   EXIT Loop
                ENDIF
             ENDIF
          ENDIF
       END DO Loop

       IF( .not.OK )THEN
          RETURN
       ELSE
          write(*,'(''   **********************************'')')
          write(*,'(''      Removing Solar Contamination'')')
          write(*,'(''   **********************************'')')
       ENDIF

       ! Now have the standard deviation, find those entries within
       ! allowed ranges which are above 2*standard_deviation
       start1_pos_final = ndata
       stop1_pos_final = 0
       start2_pos_final = ndata
       stop2_pos_final = 0
       start3_pos_final = ndata
       stop3_pos_final = 0
       DO I=1,ndata
          IF( bad(I) )THEN
             IF( NAN_R .ne. gain(I) )THEN
                gain_value(I) = coefs(1) + coefs(2)*prt(I) + &
                     coefs(3)*prt(I)*prt(I)
                gain_value(I) = sqrt(gain_value(I))
                if( 0 .gt. meangain )then
                   gain_value(I) = -gain_value(I)
                endif
                IF( ABS(gain_value(I)-gain(I)) .ge. sigma_thresh_first*stdev )THEN
                   IF( I .ge. start1_pos .and. I .le. stop1_pos )THEN
                      start1_pos_final = MIN(start1_pos_final,I)
                      stop1_pos_final = MAX(stop1_pos_final,I)
                   ELSE IF( I .ge. start2_pos .and. I .le. stop2_pos )THEN
                      start2_pos_final = MIN(start2_pos_final,I)
                      stop2_pos_final = MAX(stop2_pos_final,I)
                   ELSE IF( I .ge. start3_pos .and. I .le. stop3_pos )THEN
                      start3_pos_final = MIN(start3_pos_final,I)
                      stop3_pos_final = MAX(stop3_pos_final,I)
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       END DO

       new_good_values = good_values
       DO I=1,ndata
          IF( I .ge. start1_pos_final .and. I .le. stop1_pos_final )THEN
             new_good_values(i) = .FALSE.
          ELSE IF( I .ge. start2_pos_final .and. I .le. stop2_pos_final )THEN
             new_good_values(i) = .FALSE.
          ELSE IF( I .ge. start3_pos_final .and. I .le. stop3_pos_final )THEN
             new_good_values(i) = .FALSE.
          ENDIF
       END DO

!    do i=1,ndata
!       if( NAN_R .ne. gain(I) )then
!          gval = 0
!          bval = 0
!          if( new_good_values(I) )THEN
!             gval = 1
!          endif
!          if( bad(I) )THEN
!             bval = 1
!          endif
!          write(101,*)i,gval,bval,prt(i),gain(i),sqrt(coefs(1) + coefs(2)*prt(I) + &
!               coefs(3)*prt(I)*prt(I)),stdev,solar_za(i)
!       endif
!    end do
!    stop 1

       ! Now get new coefs and re-evaluate goodvalues
       stat = Get_Poly_Coefs( ndata, prt, gain, new_good_values, stdev, &
            maxdev, coefs )

       ! Get min/max values in window to set to bad for final case
       if( 0 .eq. stat )then
          start1_pos_final = ndata
          stop1_pos_final = 0
          start2_pos_final = ndata
          stop2_pos_final = 0
          start3_pos_final = ndata
          stop3_pos_final = 0
          DO I=1,ndata
             IF( bad(I) )THEN
                IF( NAN_R .ne. gain(I) )THEN
                   gain_value(I) = coefs(1) + coefs(2)*prt(I) + &
                        coefs(3)*prt(I)*prt(I)
                   gain_value(I) = sqrt(gain_value(I))
                   if( 0 .gt. meangain )then
                      gain_value(I) = -gain_value(I)
                   endif
                   IF( ABS(gain_value(I)-gain(I)) .ge. sigma_thresh*stdev )THEN
                      IF( I .ge. start1_pos .and. I .le. stop1_pos )THEN
                         start1_pos_final = MIN(start1_pos_final,I)
                         stop1_pos_final = MAX(stop1_pos_final,I)
                      ELSE IF( I .ge. start2_pos .and. I .le. stop2_pos )THEN
                         start2_pos_final = MIN(start2_pos_final,I)
                         stop2_pos_final = MAX(stop2_pos_final,I)
                      ELSE IF( I .ge. start3_pos .and. I .le. stop3_pos )THEN
                         start3_pos_final = MIN(start3_pos_final,I)
                         stop3_pos_final = MAX(stop3_pos_final,I)
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
          END DO

          ! Adjust widths to make sure we are missing all bad data by 200 scan lines
          ! workaround to miss everything...
          if( ndata .gt. start1_pos_final )then
             start1_pos_final = MAX(1,start1_pos_final-150)
             stop1_pos_final = MIN(ndata,stop1_pos_final+150)
          endif
          if( ndata .gt. start2_pos_final )then
             start2_pos_final = MAX(1,start2_pos_final-150)
             stop2_pos_final = MIN(ndata,stop2_pos_final+150)
          endif
          if( ndata .gt. start3_pos_final )then
             start3_pos_final = MAX(1,start3_pos_final-150)
             stop3_pos_final = MIN(ndata,stop3_pos_final+150)
          endif

          DO I=1,ndata
             IF( I .ge. start1_pos_final .and. I .le. stop1_pos_final )THEN
                good_values(i) = .FALSE.
             ELSE IF( I .ge. start2_pos_final .and. I .le. stop2_pos_final )THEN
                good_values(i) = .FALSE.
             ELSE IF( I .ge. start3_pos_final .and. I .le. stop3_pos_final )THEN
                good_values(i) = .FALSE.
             ENDIF
          END DO

       ELSE
          returnStat = stat
       ENDIF
    ELSE
       returnStat = stat
    ENDIF

    IF( 0 .ne. returnStat )THEN
       ! Set full bad range to bad and print warning
       DO I=1,ndata
          IF( bad(I) )THEN
             good_values = .FALSE.
          ENDIF
       END DO
       call out_Warning(&
            'Problem in fitting gain vs prt for solar contamination',&
            'Remove_Solar_Entries','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       call out_Warning(&
            'Setting likely range to bad',&
            'Remove_Solar_Entries','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       night_ok = .FALSE.
       return
    ENDIF

!    do i=1,ndata
!       if( NAN_R .ne. gain(I) )then
!          gval = 0
!          bval = 0
!          if( good_values(I) )THEN
!             gval = 1
!          endif
!          if( bad(I) )THEN
!             bval = 1
!          endif
!          write(101,*)i,gval,bval,prt(i),gain(i),sqrt(coefs(1) + coefs(2)*prt(I) + &
!               coefs(3)*prt(I)*prt(I)),stdev,solar_za(i)
!       endif
!    end do
!    stop 1

    ! Now check which fraction of the night time is still OK
    min_night = 1e30
    max_night = -1e30
    min_day = 1e30
    max_day = -1e30
    DO I=1,ndata
       IF( solar_za(I) .gt. 90 .and. good_values(i) )THEN
          IF( min_night .gt. prt(I) )THEN
             min_night = prt(I)
          ENDIF
          IF( max_night .lt. prt(I) )THEN
             max_night = prt(I)
          ENDIF
       ELSE IF( solar_za(I) .le. 90 .and. good_values(i) )THEN
          IF( min_day .gt. prt(I) )THEN
             min_day = prt(I)
          ENDIF
          IF( max_day .lt. prt(I) )THEN
             max_day = prt(I)
          ENDIF
       ENDIF
    END DO



    ! Must have at least 60% of the range
    night_ok = .TRUE.
    IF( 0.5 .gt. (max_night-min_night)/(max_day-min_day) )THEN
       write(*,'(''    Cannot use nighttime only for prediction'')')
       night_ok = .FALSE.
    ENDIF

    ! Deallocate memory
    DEALLOCATE(good_gain,good_gain_sq,good_prt,bad,new_good_values,gain_value)

  END FUNCTION Remove_Solar_Entries

  ! Calibration coefficients based on
  ! CCI work with double difference clear sky matches from CCI
  ! NOAA-06 -> NOAA-11 just AVHRR-to-AVHRR
  ! NOAA-12 -> MetOp-A AVHRR-to-ATSR
  !         Created 10 Aug 2016
  SUBROUTINE Get_Calib_Coefficients_CCI( time, AVHRR_No, coefs1, coefs2, &
       coefs3, ncoefs, coefs_frac1, coefs_frac2, coef_window, use_iasi )

    REAL(GbcsDble), INTENT(IN) :: time
    INTEGER, INTENT(IN) :: AVHRR_No
    REAL, INTENT(OUT) :: coefs1(8,2)
    REAL, INTENT(OUT) :: coefs2(8,2)
    REAL, INTENT(OUT) :: coefs3(8,2)
    INTEGER, INTENT(OUT) :: ncoefs
    REAL, INTENT(OUT) :: coefs_frac1
    REAL, INTENT(OUT) :: coefs_frac2
    REAL, INTENT(IN) :: coef_window
    LOGICAL, INTENT(IN) :: use_iasi

    ! Local variables
    INTEGER :: I
    INTEGER :: ipos
    INTEGER :: newpos
    INTEGER :: ntime
    REAL(GbcsDble) :: s_time(2,4)
    REAL(GbcsReal) :: s_coefs1(8,4)
    REAL(GbcsReal) :: s_coefs2(8,4)
    REAL(GbcsReal) :: s_coefs3(8,4)
    REAL(GbcsDble) :: coefwindow
    REAL(GbcsReal) :: f

    ! WIndow to interpolated case in seconds
    coefwindow = coef_window*30.*86400.

    ! Numbers are
    !
    !             1  =  bias term
    !             2  =  eta/rho correction term
    !             3  =  non-linear term
    !             4  =  nuc
    !             5  =  aVal (band coefficient)
    !             6  =  bVal (band coefficient)
    !             7  = orbital temperature term
    !
    SELECT CASE( AVHRR_No )
    CASE(1)
       ! Made up values
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 0.04475183, -0.0962706196346, 0.0, 0.0, 2683.18, 1.73290, &
            0.997424, 0. /)
       s_coefs2(:,1) = (/ -6.42783284, -0.04252789, 0.0, 1.67568E-5, 927.508, 0.405709, &
            0.998741, 0.0385670 /)
       s_coefs3(:,1) = (/ -28.22137717,-0.06568356, 0.0, 1.063571E-5, 841.403, 0.392737, &
            0.998845, 0.12536896 /)
    CASE(6) 
       ! Made up values
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 0.04475183, -0.0962706196346, 0.0, 0.0, 2683.18, 1.73290, &
            0.997424, 0. /)
       s_coefs2(:,1) = (/ -6.42783284, -0.04252789, 0.0, 1.67568E-5, 927.508, 0.405709, &
            0.998741, 0.0385670 /)
       s_coefs3(:,1) = (/ -28.22137717,-0.06568356, 0.0, 1.063571E-5, 841.403, 0.392737, &
            0.998845, 0.12536896 /)
    CASE(7) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 0.04475183, -0.0962706196346, 0.0, 0.0, 2683.18, 1.73290, &
            0.997424, 0. /)
       s_coefs2(:,1) = (/ -6.42783284, -0.04252789, 0.0, 1.67568E-5, 927.508, 0.405709, &
            0.998741, 0.0385670 /)
       s_coefs3(:,1) = (/ -28.22137717,-0.06568356, 0.0, 1.063571E-5, 841.403, 0.392737, &
            0.998845, 0.12536896 /)
    CASE(8) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -2.42738951, 0.0466958517689, 0.0, 0.0, 2671.12, 2.31640, &
            1.00146, 0.00838107 /)
       s_coefs2(:,1) = (/ -1.33257068e+02, 2.06298458e-02, 0.0, 1.67568E-5, 915.281, 0.497544, &
            0.998654, 0.46258945 /)
       s_coefs3(:,1) = (/ 0.,0.,0.,0.,0.,0.,0.,0. /)
    CASE(9) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -1.08549431, -0.0507639896943, 0.0, 0.0, 2685.84, 1.74231, &
            0.99635, 0.00375667 /)
       s_coefs2(:,1) = (/ -3.44051742e+01, -5.82090651e-03, 0.0, 1.67568E-5, 929.574, 0.347473, &
            0.998839, 0.12352714 /)
       s_coefs3(:,1) = (/ -2.34917485, -0.04441365, 0.0, 1.063571E-5, 846.292, 0.493176, &
            0.998698, 0.0305383 /)
    CASE(10) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 0.08105582, 0.0444195463436, 0.0, 0.0, 2679.11, 1.93677, &
            0.998787, -0.00028523 /)
       s_coefs2(:,1) = (/ -10.43759699, 0.03732374, 0.0, 1.67568E-5, 910.515, 0.457244, &
            0.998768, 0.02932454 /)
       s_coefs3(:,1) = (/ 0.,0.,0.,0.,0.,0.,0.,0. /)
    CASE(11) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 0.00641791, 0.00367610714174, 0.0, 0.0, 2682.85, 1.68218, &
            0.997630, 0. /)
       s_coefs2(:,1) = (/ -9.3096892, 0.04557533, 0.0, 1.67568E-5, 928.047, 0.393074, &
            0.998759, 0.02040797 /)
       s_coefs3(:,1) = (/ -2.22947502e+01, 0.0, 1.65224150e-02, 1.063571E-5, 842.869, 0.406868, &
            0.998803, 0.07508146 /)
    CASE(12) 
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -0.6812111, 0.011175775445, 0.0, 0.0, 2671.81, 2.42853, &
            1.00102, 0.00237008 /)
       s_coefs2(:,1) = (/ -4.55986265e+01, -2.12927637e-03, 0.0, 1.67568E-5, 921.925, 0.560212, &
            0.998415, 0.16380201 /)
       s_coefs3(:,1) = (/ -5.58716406e+01, -1.06646194e-02, 0.0, 1.063571E-5, 838.003, 0.398581, &
            0.99883, 0.20447873 /)
    CASE(14)
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -0.22671186, 0.00411023343761, 0.0, 0.0, 2674.39, 2.39486, &
            1.00004, 0.00079423 /)
       s_coefs2(:,1) = (/ -2.49531168e+01, -1.06136789e-02, 0.0, 1.67568E-5, 930.521, 0.621590, &
            0.998294, 0.09372974 /)
       s_coefs3(:,1) = (/ -2.41356541e+01, -3.19547251e-03, 0.0, 1.063571E-5, 835.544, 0.419563, &
            0.998769, 0.09049913 /)
    CASE(15)
       ntime = 3
       s_time(:,1) = (/ 157766400.0d0, 883612800.0d0 /)
       s_time(:,2) = (/ 883612800.0d0, 1009843200.0d0 /)
       s_time(:,3) = (/ 1009843200.0d0, 3944678400.0d0 /)
       s_time(:,3) = (/ 1009843200.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -0.35925326, -0.00391622918693, 0.0, 0.0, 2693.67, 1.79289, &
            0.998704, 0.0012398 /)
       s_coefs2(:,1) = (/ -1.46969422e+01, -3.78504627e-03, 0.0, 1.6756799999999999e-05, 930.725, &
            0.533396, 0.998533, 0.05492564 /)
       s_coefs3(:,1) = (/ -1.72589079e+01, 4.23971234e-03, 0.0, 1.063571e-05, 838.477, 0.366515, &
            0.998930, 0.06061977 /)
       s_coefs1(:,2) = (/ -0.15675434, -0.0139010000165, 0.0, 0.0, 2693.67, 1.79289, &
            0.998704, 0.00056521 /)
       s_coefs2(:,2) = (/ -9.05642510e+00, -7.01951004e-03, 0.0, 1.6756799999999999e-05, 930.725, &
            0.533396, 0.998533, 0.03685354 /)
       s_coefs3(:,2) = (/ -8.80624181e+00, -1.71202854e-03, 0.0, 1.063571e-05, 838.477, 0.366515, &
            0.998930, 0.03411435 /)
       s_coefs1(:,3) = (/ -0.14610468, 0.0174068629097, 0.0, 0.0, 2693.67, 1.79289, &
            0.998704, 0.00049053 /)
       s_coefs2(:,3) = (/ -1.11019815e+01, 3.24968094e-03, 0.0, 1.6756799999999999e-05, 930.725, &
            0.533396, 0.998533, 0.04092865 /)
       s_coefs3(:,3) = (/ -1.08739078e+01, 8.09022450e-03, 0.0, 1.063571e-05, 838.477, 0.366515, &
            0.998930, 0.03795704 /)
    CASE(16)
       ntime = 4
       s_time(:,1) = (/ 157766400.0d0, 907200000.0d0 /)
       s_time(:,2) = (/ 907200000.0d0, 1022889600.0d0 /)
       s_time(:,3) = (/ 1022889600.0d0, 1073001600.0d0 /)
       s_time(:,4) = (/ 1073001600.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 1.49155438, -0.000902624969744, 0.0, 0.0, 2700.63, 1.61173, &
            0.997916, -0.00517501 /)
       s_coefs2(:,1) = (/ 8.03442772e+01, -2.91784882e-02, 0.0, 1.6756799999999999e-05, 918.956, &
            0.613784, 0.998183, -0.26635246 /)
       s_coefs3(:,1) = (/ 1.21091471e+02, -9.50955067e-03, 0.0, 1.063571e-05, 834.208, 0.359189, &
            0.998856, -0.41316913 /)
       s_coefs1(:,2) = (/ -0.10608244, 0.00151247187782, 0.0, 0.0, 2700.63, 1.61173, &
            0.997916, 0.00036288 /)
       s_coefs2(:,2) = (/ -22.55569832, -0.02482788, 0.0, 1.6756799999999999e-05, 918.956, &
            0.613784, 0.998183, 0.08896217 /)
       s_coefs3(:,2) = (/ -1.82549018e+01, -5.93550258e-03, 0.0, 1.063571e-05, 834.208, 0.359189, &
            0.998856, 0.0687013 /)
       s_coefs1(:,3) = (/ -0.34916101, 0.012558837755, 0.0, 0.0, 2700.63, 1.61173, &
            0.997916, 0.00117254 /)
       s_coefs2(:,3) = (/ -2.45476698e+01, -2.33725553e-02, 0.0, 1.6756799999999999e-05, 918.956, &
            0.613784, 0.998183, 0.09549935 /)
       s_coefs3(:,3) = (/ -2.06872378e+01, -5.61192256e-03, 0.0, 1.063571e-05, 834.208, 0.359189, &
            0.998856, 0.0770247 /)
       s_coefs1(:,4) = (/ -0.25389726, -0.00589812377425, 0.0, 0.0, 2700.63, 1.61173, &
            0.997916, 0.00087232 /)
       s_coefs2(:,4) = (/ -23.40616806, -0.03262973, 0.0, 1.6756799999999999e-05, 918.956, &
            0.613784, 0.998183, 0.0943693 /)
       s_coefs3(:,4) = (/ -1.75677971e+01, -1.40932079e-02, 0.0, 1.063571e-05, 834.208, 0.359189, &
            0.998856, 0.06938199 /)
    CASE(17)
       ntime = 2
       s_time(:,1) = (/ 157766400.0d0, 1112313600.0d0 /)
       s_time(:,2) = (/ 1112313600.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -0.12485158, 0.00204537373786, 0.0, 0.0, 2665.02, 1.92299, &
            0.998811, 0.00044127 /)
       s_coefs2(:,1) = (/ 3.72410471, 0.03647545, 0.0, 1.6756799999999999e-05, 942.392, 0.513380, &
            0.998565, -0.02063181 /)
       s_coefs3(:,1) = (/ -5.39759017, 0.00895862, 0.0, 1.063571e-05, 840.489, 0.290894, &
            0.999043, 0.01924988 /)
       s_coefs1(:,2) = (/ 0.30596653, 0.00289277795388, 0.0, 0.0, 2665.02, 1.92299, &
            0.998811, -0.00103927 /)
       s_coefs2(:,2) = (/ 11.96393153, 0.03662023, 0.0, 1.6756799999999999e-05, 942.392, 0.513380, &
            0.998565, -0.04903655 /)
       s_coefs3(:,2) = (/ -2.64244377, 0.01028158, 0.0, 1.063571e-05, 840.489, 0.290894, &
            0.999043, 0.00919356 /)
    CASE(18)
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -0.1313061, 0.00967957245062, 0.0, 0.0, 2658.37, 1.79343, &
            0.997407, 0.00045377 /)
       s_coefs2(:,1) = (/ -0.06209035, 0.01253979, 0.0, 1.67568E-5, 934.436, 0.536755, &
            0.998524, 0.0 /)
       s_coefs3(:,1) = (/ 9.98239547e-01, 4.26128628e-04, 0.0, 1.063571E-5, 834.548, 0.324453, &
            0.998947, 0.0 /)
    CASE(19)
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ -1.18171985, 0.0189763521632, 0.0, 0.0, 2670.00, 1.673956, &
            0.9973640, 0.00409561 /)
       s_coefs2(:,1) = (/ -3.22790849e+01, 1.20640747e-02, 0.0, 1.67568E-5, 928.90, 0.5395879, &
            0.9985342, 0.11224999 /)
       s_coefs3(:,1) = (/ -2.10300447e+01, 1.32912830e-02, 0.0, 1.063571E-5, 831.90, 0.3606376, &
            0.9989132, 0.07220871 /)
    CASE(-1)
       ntime = 1
       s_time(:,1) = (/ 157766400.0d0, 3944678400.0d0 /)
       s_coefs1(:,1) = (/ 1.11362134, 0.00179079409535, 0.0, 0.0, 2684.29, 2.27349, &
            0.997935, -0.00388403 /)
       IF( use_iasi )THEN
          !
          ! Note here we are using the IASI/MetOp-A calibration derived in 
          ! Mittaz et al. 2011
          ! We only have IASI/AVURR coefficients for the 11/12 micron channels
          !
          s_coefs2(:,1) = (/ 1.37261, 0.0403441, -0.00564703, 1.82018e-05, 927.200, 0.551260, &
               0.998509, 0. /)
          s_coefs3(:,1) = (/ 1.11838, -0.00479986, 0.00111601, 1.10115E-5, 837.7, 0.347160, &
               0.998947, 0. /)
       ELSE
          !
          ! Here are the equivilent values from the AATSR double difference which
          ! we are not using in CCI
          !
          s_coefs2(:,1) = (/ 6.88329791e+01, -6.01531058e-03, 0.0, 1.67568E-5, 926.566, 0.456849, &
               0.998596, -0.23486149 /)
          s_coefs3(:,1) = (/ 5.09659390e+01, 7.05761507e-04, 0.0, 1.063571E-5, 836.344, 0.136754, &
               0.999260, -0.17470683 /)
       ENDIF
    CASE DEFAULT
       call out_ERROR(&
            'Cannot match AVHRR number with allowed cases (7-19,MetOp-A)',&
            'Get_Calib_Coefficients_CCI','NOAA_LoadAVHRRLevel1B.f90',&
            .FALSE.,' ')
    END SELECT

    if( ntime .eq. 1 )then
       ncoefs = 1
       coefs1(:,1) = s_coefs1(:,1)
       coefs2(:,1) = s_coefs2(:,1)
       coefs3(:,1) = s_coefs3(:,1)
    else
       !
       ! Find which time zone we are in
       !
       ipos = -1
       do I=1,ntime
          if( time .ge. s_time(1,I) .and. time .le. s_time(2,I) )then
             ipos = i
          endif
       end do
       if( -1 .eq. ipos )then
          call out_ERROR('Cannot match input time to coefficients','Get_Calib_Coefficients_CCI',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
       
       coefs1(:,1) = s_coefs1(:,ipos)
       coefs2(:,1) = s_coefs2(:,ipos)
       coefs3(:,1) = s_coefs3(:,ipos)
       !
       ! See if we are within the window to try and interpolate Radiances
       !
       f = -1e30
       IF( ipos .gt. 1 )THEN
          IF( (time - s_time(1,ipos)) .lt. coefwindow )THEN
             f = (coefwindow+(time - s_time(1,ipos)))/(2.*coefwindow)
             newpos = ipos-1
          ENDIF
       ENDIF
       IF( ipos .lt. ntime )THEN
          IF( (s_time(2,ipos)-time) .lt. coefwindow )THEN
             f = (coefwindow+(s_time(2,ipos)-time))/(2.*coefwindow)
             newpos = ipos+1
          ENDIF
       ENDIF
       if( f .gt. -1e20 )then
          coefs_frac1 = f
          coefs_frac2 = (1-f)
          coefs1(:,2) = s_coefs1(:,newpos)
          coefs2(:,2) = s_coefs2(:,newpos)
          coefs3(:,2) = s_coefs3(:,newpos)
          ncoefs=2
       else
          ncoefs=1
       endif

    endif

  END SUBROUTINE Get_Calib_Coefficients_CCI

  ! Calibration coefficients based on
  ! Pre-launch : Mittaz, Harris & Sullivan (2008)
  ! In-Orbit : Mittaz & Harris, 2nd MERIS/ATSR Users Meeting, Frascati, Italy,
  !            22-26 Sept 2008

  SUBROUTINE Get_Calib_Coefficients_OLD( AVHRR_No, coefs1, coefs2, coefs3, &
       pre_launch )

    INTEGER, INTENT(IN) :: AVHRR_No
    REAL, INTENT(OUT) :: coefs1(7)
    REAL, INTENT(OUT) :: coefs2(7)
    REAL, INTENT(OUT) :: coefs3(7)
    LOGICAL, OPTIONAL :: pre_launch

    LOGICAL :: pre_launch_data

    if( PRESENT(pre_launch) )then
       pre_launch_data = pre_launch
    else
       pre_launch_data = .FALSE.
    endif

    ! Numbers are
    !
    !             1  =  2nd order (gamma) term
    !             2  =  rho (extra scattered radiation from instrument body)
    !             3  =  nuc
    !             4  =  aVal (band coefficient)
    !             5  =  bVal (band coefficient)
    !             6  =  constant offset (ECT) (alpha)
    !             7  =  ICT constant offset   (alpha')
    !
    IF( pre_launch_data )THEN

       SELECT CASE( AVHRR_No )
       CASE(15) ! NOAA-15
          coefs1 = (/ -5.0000000e-09, 0.014021210, 2693.67, 1.79289, &
               0.998704, 0.00061680884, 0.0049695020 /)
          coefs2 = (/ 1.7445034e-05, 0.029412672, 930.725, 0.533396, &
               0.998533, 0.55784984, 2.2423612 /)
          coefs3 = (/ 1.1596027e-05, 0.023182540, 838.477, 0.366515, &
               0.998930, 0.47486476, 1.7638984 /)
       CASE(16) ! NOAA-16
          coefs1 = (/ -5.0000000e-09, 0.014300421, 2700.63, 1.61173, &
               0.997916, 0.00020755702, 0.0049187551 /)
          coefs2 = (/ 8.7544046e-06, 0.030398189, 921.989, 0.661622, &
               0.998165, 0.57927427, 2.5610334 /)
          coefs3 = (/ 3.7901654e-06, 0.028376910, 831.384, 0.408960, &
               0.998784, 0.47407687, 2.4094642 /)
       CASE(17) ! NOAA-17
          coefs1 = (/ -5.0000000e-09, 0.019350447, 2665.02, 1.92299, &
               0.998811, 0.00048679, 0.0055925654 /)
          coefs2 = (/ 3.2261222e-05, 0.0297, 942.392, 0.513380, &
               0.998565, 0.705993, 2.4285780 /)
          coefs3 = (/ 1.3657313e-05, 0.0234, 840.489, 0.290894, &
               0.999043, 0.61184829, 1.8537181 /)
       CASE(18) ! NOAA-18
          coefs1 = (/ -5.0000000e-09, 0.012621663, 2658.37, 1.79343, &
               0.997407, 0.00095855, 0.0024209759 /)
          coefs2 = (/ 2.0571179e-05, 0.019886473, 934.436, 0.536755, &
               0.998524, 0.77172879, 1.5294463 /)
          coefs3 = (/ 7.6080007e-06, 0.017329559, 834.548, 0.324453, &
               0.998947, 0.5186498, 1.0483270 /)
       CASE(19) ! NOAA-19
          coefs1 = (/ 0.0, 0.0, 2670.00, 1.673956, &
               0.9973640, 0.0, 0.0 /)
          coefs2 = (/ 0.0, 0.0, 928.90, 0.5395879, &
               0.9985342, 0.0, 0.0 /)
          coefs3 = (/ 0.0, 0.0, 831.90, 0.3606376, &
               0.9989132, 0.0, 0.0 /)
       CASE(-1) ! Metop-A
          coefs1 = (/ -5.0000000e-09, 0.021368549, 2684.29, 2.27349, &
               0.997935, 0.00095782979, 0.0057592815 /)
          coefs2 = (/ 1.8444566e-05, 0.027956921, 931.826, 0.397824, &
               0.999694, 0.84437795, 2.2689228 /)
          coefs3 = (/ 1.1196305E-05, 0.026021140, 838.374, 0.268040, &
               0.999062, 0.62847434, 2.0601701 /)
       CASE DEFAULT
          call out_ERROR(&
               'Cannot match AVHRR number with allowed cases (15-19,MetOp-A)',&
               'Get_Calib_Coefficients','NOAA_LoadAVHRRLevel1B.f90',&
               .FALSE.,' ')
       END SELECT
    ELSE
       ! Set alpha' and rho to zero (NOAA-17 ATSR fits)
       ! Set alpha to 0.9 for NOAA-17 (filter shift implied by NOAA-17 ATSR fits)
       SELECT CASE( AVHRR_No )
       CASE(14) ! NOAA-14 - Made up numbers for now (apart from nuC,A,B)
          coefs1 = (/ -5.0000000e-09, 0.0, 2683.14, 2.65776, &
               1.00169, 0.0014, 0.0 /) ! guess at calibration
          coefs2 = (/ 1.8307023e-05, -0.00465605, 930.514, 0.627327, &
               0.998278, 1.23221666, 0.0371985 /)
          coefs3 = (/ 1.0808404e-05, -0.000478083, 835.840, 0.468771, &
               0.998695, 1.082536, 0.00109675 /)
       CASE(15) ! NOAA-15
          coefs1 = (/ -5.0000000e-09, 0.0, 2693.67, 1.79289, &
               0.998704, 0.00061680884, 0.0 /)
          coefs2 = (/ 1.7445034e-05, 0.0, 930.725, 0.533396, &
               0.998533, 0.55784984, 0.0 /)
          coefs3 = (/ 1.1596027e-05, 0.0, 838.477, 0.366515, &
               0.998930, 0.47486476, 0.0 /)
       CASE(16) ! NOAA-16
          coefs1 = (/ -5.0000000e-09, 0.0, 2700.63, 1.61173, &
               0.997916, 0.00020755702, 0.0 /)
          coefs2 = (/ 8.6071114e-06, -0.0056470302, 918.956, 0.613784, &
               0.998183, 1.29416, 0.0 /)
          coefs3 = (/ 4.0256298e-06, 0.0, 834.208, 0.359189, &
               0.998856, 1.16071, 0.0 /)
       CASE(17) ! NOAA-17
          coefs1 = (/ -5.0000000e-09, 0.0, 2665.02, 1.92299, &
               0.998811, 0.00048679, 0.0 /)
          coefs2 = (/ 3.2261222e-05, 0.0, 942.392, 0.513380, &
               0.998565, 0.9, 0.0 /)
          coefs3 = (/ 1.3657313e-05, 0.0, 840.489, 0.290894, &
               0.999043, 0.9, 0.0 /)
       CASE(18) ! NOAA-18
          coefs1 = (/ -5.0000000e-09, 0.0, 2658.37, 1.79343, &
               0.997407, 0.00095855, 0.0 /)
          coefs2 = (/ 2.0571179e-05, 0.0, 934.436, 0.536755, &
               0.998524, 0.77172879, 0.0 /)
          coefs3 = (/ 7.6080007e-06, 0.0, 834.548, 0.324453, &
               0.998947, 0.5186498, 0.0 /)
       CASE(19) ! NOAA-19
          coefs1 = (/ 0.0, 0.0, 2670.00, 1.673956, &
               0.9973640, 0.0, 0.0 /)
          coefs2 = (/ 0.0, 0.0, 928.90, 0.5395879, &
               0.9985342, 0.0, 0.0 /)
          coefs3 = (/ 0.0, 0.0, 831.90, 0.3606376, &
               0.9989132, 0.0, 0.0 /)
       CASE(-1) ! Metop-A
!          coefs1 = (/ -5.0000000e-09, 0.0, 2684.29, 2.27349, &
!               0.997935, 0.00095782979, 0.0 /)
!          coefs2 = (/ 1.8444566e-05, 0.0, 931.826, 0.397824, &
!               0.999694, 0.84437795, 0.0 /)
!          coefs3 = (/ 1.1196305E-05, 0.0, 838.374, 0.268040, &
!               0.999062, 0.62847434, 0.0 /)
          ! Values for 11/12 micron derived from IASI analysis
          ! Modified to reflect IASI fits = May 7 2010
          ! 3.7 mu =>
          !
          coefs1 = (/ -5.0000000e-09, 0.0, 2684.29, 2.27349, &
               0.997935, 0.0014, 0.0 /) ! guess at calibration
          ! 11 mu => alpha_zero=0.84879466 nu=931.604 a=0.403139 b=0.998687
          !       => non-linear=1.8307023e-05
          ! aval => 0.383422,0.0371985,-0.00469005,5.30000
          ! Note - changed variables to zero shift case so can match
          ! with RTMs
          coefs2 = (/ 1.82018e-05, -0.00564703, 927.200, 0.551260, &
               0.998509, 1.37261, 0.0403441 /)
          ! 12 mu => alpha_zero=0.65076100 nu=837.453 a=0.282462 b=0.999034
          !       => non-linear=1.0808404e-05
          ! aval => 0.431775,0.00109675,-0.000478083,0.287042
          ! Note - changed variables to zero shift case so can match
          ! with RTMs
          coefs3 = (/ 1.10115e-05, 0.00111601, 837.700,  0.347160, &
               0.998947, 1.11838, -0.00479986 /)
       CASE DEFAULT
          call out_ERROR(&
               'Cannot match AVHRR number with allowed cases (14-19,MetOp-A)',&
               'Get_Calib_Coefficients','NOAA_LoadAVHRRLevel1B.f90',&
               .FALSE.,' ')
       END SELECT
    ENDIF

  END SUBROUTINE Get_Calib_Coefficients_OLD

  SUBROUTINE Parse_DayNo( year, dayno, month, day )

    INTEGER, INTENT(IN) :: year
    INTEGER, INTENT(IN) :: dayno
    INTEGER, INTENT(OUT) :: month
    INTEGER, INTENT(OUT) :: day

    INTEGER :: i

    IF( (0 .eq. MOD(year,4) .and. 0 .ne. MOD(year,100)) .or. &
         (0 .eq. MOD(year,400)) )then
       leaploop: do i=2,13
          if( dayno .le.  dayNumberLeap(i) )then
             month = i-1
             day = dayno - dayNumberLeap(i-1)
             exit leaploop
          endif
       end do leaploop
    ELSE
       leaploop2: do i=2,13
          if( dayno .le. dayNumber(i) )then
             month = i-1
             day = dayno - dayNumber(i-1)
             exit leaploop2
          endif
       end do leaploop2
    ENDIF

  END SUBROUTINE Parse_DayNo

  SUBROUTINE Copy_To_IMG( Coefs, outData, IMG, new_calibration, &
       use_old_37_calibration, start_valid, stop_valid, &
       load_landmask,new_vis_calibration)

#ifdef NOAA_Compile
    USE MCIDAS_LoadImagery, ONLY: Read_In_Total_Landmask, &
         Add_To_Image_Start_Time
#endif

    TYPE(AVHRR_Instrument_Coefs), INTENT(IN) :: Coefs
    TYPE(AVHRR_Data), POINTER :: outData
    TYPE(Imagery), INTENT(INOUT) :: IMG
    LOGICAL, INTENT(IN) :: new_calibration
    LOGICAL, INTENT(IN) :: use_old_37_calibration
    INTEGER, INTENT(OUT) :: start_valid
    INTEGER, INTENT(OUT) :: stop_valid
    LOGICAL, INTENT(IN) :: load_landmask
    LOGICAL, INTENT(IN), OPTIONAL :: new_vis_calibration

    INTEGER :: nelem, nline, line
    INTEGER :: STAT

    INTEGER :: ii
    INTEGER :: FIndx
    integer :: inthour
    real :: val
    integer :: intmin
    integer :: intsec
    integer :: intsec1000
    real :: image_time_secs
    LOGICAL :: new_vis_cal

    IF( PRESENT(new_vis_calibration) )THEN
       IF( outData%newCalibration_There )THEN
          new_vis_cal = new_vis_calibration
       ELSE
          new_vis_cal = .FALSE.
       ENDIF
    ELSE
       IF( outData%newCalibration_There )THEN
          new_vis_cal = .TRUE.
       ELSE
          new_vis_cal = .FALSE.
       ENDIF
    ENDIF
    if( new_calibration .and. outData%newCalibration_There )THEN
       WRITE(*,'(''  Copying recalibrated radiances/BTs to IMG'')')
    else
       WRITE(*,'(''  Copying NOAA Level1B radiances/BTs to IMG'')')
    endif

    IF( new_vis_cal .and. outData%newCalibration_There )THEN
       WRITE(*,'(''  Copying recalibrated visible radiances to IMG'')')
    else
       WRITE(*,'(''  Copying NOAA Level1B visible radiances to IMG'')')
    endif
       

    if( .not. outData%dataFilled )then
       call out_WARNING('outData/IMG not filled with data',&
            'Copy_To_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       IMG%Dims%elem = 0
       IMG%Dims%line = 0
       IMG%Npts%elem = 0
       IMG%Npts%line = 0
       outData%start_valid = -1
       outData%stop_valid = -1
       return
    endif

    nelem = outData%nelem
    nline = outData%arraySize

    IMG%Dims%elem = outData%nelem
    IMG%Npts%elem = outData%nelem
    IMG%Dims%line = outData%arraySize
    IMG%Npts%line = outData%arraySize
    IMG%Origin%elem = 1
    IMG%Origin%line = 1
    IMG%Scan_Start_End_Available = .TRUE.

    start_valid = -1
    firstloop: do ii=1,outData%arraySize
       IF( .NOT. outData%badTime(ii) )THEN
          IMG%Start_Time%year = outData%year(ii)
          IMG%Start_Time%month = outData%month(ii)
          IMG%Start_Time%day = outData%day(ii)
          inthour = INT(outData%hours(ii))
          val = (outData%hours(ii) - inthour)*60.
          intmin = INT(val)
          val = (val - intmin)*60
          intsec = INT(val)
          val = val - intsec
          intsec1000 = NINT(val*1000)
          IMG%Start_Time%hour = inthour
          IMG%Start_Time%minute = intmin
          IMG%Start_Time%seconds = intsec
          IMG%Start_Time%sec1000 = intsec1000
          start_valid = ii
          exit firstloop
       endif
    end do firstloop

    stop_valid = -1
    secondloop: do ii=outData%arraySize,1,-1
       if( .NOT. outData%badTime(ii) )then
          IMG%End_Time%year = outData%year(ii)
          IMG%End_Time%month = outData%month(ii)
          IMG%End_Time%day = outData%day(ii)
          inthour = INT(outData%hours(ii))
          val = (outData%hours(ii) - inthour)*60.
          intmin = INT(val)
          val = (val - intmin)*60.
          intsec = INT(val)
          val = val - intsec
          intsec1000 = NINT(val*1000)
          IMG%End_Time%hour = inthour
          IMG%End_Time%minute = intmin
          IMG%End_Time%seconds = intsec
          IMG%End_Time%sec1000 = intsec1000
          stop_valid = ii
          exit secondloop
       endif
    end do secondloop

    outData%start_valid = start_valid
    outData%stop_valid = stop_valid

#ifdef NOAA_Compile
    ! Set mean time
    IF( start_valid .ne. -1 )THEN
       ! stop valid .ne. -1 also
       image_time_secs = outData%time(stop_valid) - outData%time(start_valid)
       call add_to_image_start_time( IMG%Start_Time, IMG%Mean_Time, &
            Image_Time_Secs/2. )
    ENDIF
#endif

    ! Reset start/stop valid to we match the right things
    start_valid = 1
    stop_valid = outData%ArraySize

    SELECT CASE( outData%AVHRR_No )
    CASE(1)
       IMG%AVHRR_name = 'TIROS-N'
    CASE(6)
       IMG%AVHRR_name = 'NOAA-6'
    CASE(7)
       IMG%AVHRR_name = 'NOAA-7'
    CASE(8)
       IMG%AVHRR_name = 'NOAA-8'
    CASE(9)
       IMG%AVHRR_name = 'NOAA-9'
    CASE(10)
       IMG%AVHRR_name = 'NOAA-10'
    CASE(11)
       IMG%AVHRR_name = 'NOAA-11'
    CASE(12)
       IMG%AVHRR_name = 'NOAA-12'
    CASE(13)
       IMG%AVHRR_name = 'NOAA-13'
    CASE(14)
       IMG%AVHRR_name = 'NOAA-14'
    CASE(15)
       IMG%AVHRR_name = 'NOAA-15'
    CASE(16)
       IMG%AVHRR_name = 'NOAA-16'
    CASE(17)
       IMG%AVHRR_name = 'NOAA-17'
    CASE(18)
       IMG%AVHRR_name = 'NOAA-18'
    CASE(19)
       IMG%AVHRR_name = 'NOAA-19'
    CASE(-1)
       IMG%AVHRR_name = 'Metop-02'
    CASE DEFAULT
       IMG%AVHRR_name = 'Unknown'
    END SELECT

    ! Allocate and fill
    if( ALLOCATED(IMG%ICT_temp) )then
       DEALLOCATE(IMG%ICT_temp)
    endif
    ALLOCATE(IMG%ICT_temp(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%ICT_temp',&
            'Copy_To_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    ! all AVHRRs (so far (still to see Metop-01 and -03) have
    ! the PRT weighting = 0.25 for all PRTs
    WHERE( NAN_R .NE. outData%prt1 .AND. &
         NAN_R .NE. outData%prt2 .AND. &
         NAN_R .NE. outData%prt3 .AND. &
         NAN_R .NE. outData%prt4 )
       IMG%ICT_temp = 0.25 * (outData%prt1 + outData%prt2 &
            + outData%prt3 + outData%prt4 )
    ELSEWHERE
       IMG%ICT_temp = NAN_R
    END WHERE

    if( ALLOCATED(IMG%scanLineNumber) )then
       DEALLOCATE(IMG%scanLineNumber)
    endif
    ALLOCATE(IMG%scanLineNumber(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%scanLineNumber',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IMG%scanLineNumber = outdata%scanLineNumber

    if( ALLOCATED(IMG%badTime) )then
       DEALLOCATE(IMG%badTime)
    endif
    ALLOCATE(IMG%badTime(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%badTime',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IMG%badTime = outdata%badTime

    if( ALLOCATED(IMG%badNavigation) )then
       DEALLOCATE(IMG%badNavigation)
    endif
    ALLOCATE(IMG%badNavigation(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%badNavigation',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IMG%badNavigation = outdata%badNavigation

    if( ALLOCATED(IMG%badCalibration) )then
       DEALLOCATE(IMG%badCalibration)
    endif
    ALLOCATE(IMG%badCalibration(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%badCalibration',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IMG%badCalibration = outdata%badCalibration

    if( ALLOCATED(IMG%scanline_time) )then
       DEALLOCATE(IMG%scanline_time)
    endif
    ALLOCATE(IMG%scanline_time(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%scanline_time',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    ! this can be inverted with invert_AVHRR_Date
    IMG%scanline_time = AVHRR_Date(outdata%year, outdata%dayNo, outdata%UTC_msecs)

    if( ALLOCATED(IMG%scanline_dtime) )then
       DEALLOCATE(IMG%scanline_dtime)
    endif
    ALLOCATE(IMG%scanline_dtime(nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%scanline_dtime',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IF( outData%start_valid .ne. -1 )THEN
       IMG%Reference_Time = IMG%Start_Time
       DO ii=1,outData%arraySize
          IF( outData%time(ii) .gt. 0)THEN
             IMG%scanline_dtime(ii) = SNGL(outdata%time(ii) - outdata%time(outData%start_valid))
          ENDIF
       END DO
       WHERE (outData%badTime) IMG%scanline_dtime = NAN_R
    END IF

    if( ALLOCATED(IMG%Lat) )then
       DEALLOCATE(IMG%Lat)
    endif
    ALLOCATE(IMG%Lat(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%Lat',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF

    IMG%Lat = outData%Lat

    if( ALLOCATED(IMG%Lon) )then
       DEALLOCATE(IMG%Lon)
    endif
    ALLOCATE(IMG%Lon(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%Lon',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%Lon = outData%Lon

    if( ALLOCATED(IMG%satZA) )then
       DEALLOCATE(IMG%satZA)
    endif
    ALLOCATE(IMG%satZA(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%satZA',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%satZA = outData%satZA
    if( ALLOCATED(IMG%solZA) )then
       DEALLOCATE(IMG%solZA)
    endif
    ALLOCATE(IMG%solZA(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%solZA',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%solZA = outData%solZA

    if( ALLOCATED(IMG%relAz) )then
       DEALLOCATE(IMG%relAz)
    endif
    ALLOCATE(IMG%relAz(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%relAz',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%relAz = outData%relAz

    ! Get land mask
    IF( load_landmask )THEN
#ifdef NOAA_Compile
       CALL Read_In_Total_Landmask( IMG, IMG%GbcsDataPath )
#endif
    ENDIF


    ! Now copy radiance etc. data
    DO ii = 1,IMG%N_Chans

       FIndx = IMG%Gbcs_Map(ii)

       SELECT CASE(FIndx)
       CASE(Gbcs_VIS_006)
          IMG%PReflToRad(Gbcs_VIS_006) = Coefs%F(1) / (100.0 * M_PI * Coefs%w(1))
          if( ALLOCATED(IMG%Chan_Data(Gbcs_VIS_006)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_VIS_006)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_VIS_006)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(VIS_006)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          IF( new_vis_cal )THEN
             IMG%Chan_Data(Gbcs_VIS_006)%d = outData%new_array1
          ELSE
             IMG%Chan_Data(Gbcs_VIS_006)%d = outData%array1
          ENDIF
          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_VIS_006)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_VIS_006)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_VIS_006)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(VIS_006)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_VIS_006)
          ENDIF
       CASE(Gbcs_VIS_008)
          IMG%PReflToRad(Gbcs_VIS_008) = Coefs%F(2) / (100.0 * M_PI * Coefs%w(2))
          if( ALLOCATED(IMG%Chan_Data(Gbcs_VIS_008)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_VIS_008)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_VIS_008)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(VIS_008)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          IF( new_vis_cal )THEN
             IMG%Chan_Data(Gbcs_VIS_008)%d = outData%new_array2
          ELSE
             IMG%Chan_Data(Gbcs_VIS_008)%d = outData%array2
          ENDIF
          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_VIS_008)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_VIS_008)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_VIS_008)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(NIR_01)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_VIS_008)
          ENDIF
       CASE(Gbcs_NIR_01x)
          IMG%PReflToRad(Gbcs_NIR_01x) = Coefs%F(3) / (100.0 * M_PI * Coefs%w(3))
          if( ALLOCATED(IMG%Chan_Data(Gbcs_NIR_01x)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_NIR_01x)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_NIR_01x)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(NIR_01)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          IF( new_vis_cal )THEN
             IMG%Chan_Data(Gbcs_NIR_01x)%d = outData%new_array3A
          ELSE
             IMG%Chan_Data(Gbcs_NIR_01x)%d = outData%array3A
          ENDIF

          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_NIR_01x)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_NIR_01x)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_NIR_01x)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(NIR_01)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_NIR_01x)
          ENDIF

       CASE(Gbcs_IR_03x)
          if( ALLOCATED(IMG%Chan_Data(Gbcs_IR_03x)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_IR_03x)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_IR_03x)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(IR_03x)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          ! If new calibration there and want to use it for 3.7 mu channel
          if( new_calibration .and. outData%newCalibration_There &
               .and. .not.use_old_37_calibration)then
             IMG%Chan_Data(Gbcs_IR_03x)%d = outData%new_array3B
          else
             IMG%Chan_Data(Gbcs_IR_03x)%d = outData%array3B
          endif
          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_IR_03x)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_IR_03x)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_IR_03x)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(IR_03x)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_IR_03x)
          ENDIF
       CASE(Gbcs_IR_10x)
          if( ALLOCATED(IMG%Chan_Data(Gbcs_IR_10x)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_IR_10x)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_IR_10x)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(IR_10x)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          if( new_calibration .and. outData%newCalibration_There )then
             IMG%Chan_Data(Gbcs_IR_10x)%d = outData%new_array4
          else
             IMG%Chan_Data(Gbcs_IR_10x)%d = outData%array4
          endif

          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_IR_10x)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_IR_10x)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_IR_10x)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(IR_10x)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_IR_10x)
          ENDIF

       CASE(Gbcs_IR_12x)
          if( ALLOCATED(IMG%Chan_Data(Gbcs_IR_12x)%d) )then
             DEALLOCATE(IMG%Chan_Data(Gbcs_IR_12x)%d)
          endif
          ALLOCATE(IMG%Chan_Data(Gbcs_IR_12x)%d(nelem,nline),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             CALL out_ERROR('Cannot ALLOCATE IMG%Chan_data(IR_12x)%d',&
                  'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
          ENDIF
          if( new_calibration .and. outData%newCalibration_There )then
             IMG%Chan_Data(Gbcs_IR_12x)%d = outData%new_array5
          else
             IMG%Chan_Data(Gbcs_IR_12x)%d = outData%array5
          endif
          IF( load_landmask )THEN
             if( ALLOCATED(IMG%Chan_LSD(Gbcs_IR_12x)%d) )then
                DEALLOCATE(IMG%Chan_LSD(Gbcs_IR_12x)%d)
             endif
             ALLOCATE(IMG%Chan_LSD(Gbcs_IR_12x)%d(nelem,nline),STAT = STAT)
             IF( 0 .ne. STAT )THEN
                CALL out_ERROR('Cannot ALLOCATE IMG%Chan_LSD(IR_12x)%d',&
                     'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
             ENDIF
             CALL AVHRR_Image_Calculate_LSD(IMG, Gbcs_IR_12x)
          ENDIF
       CASE DEFAULT
          CALL out_ERROR('Cannot get requested filter from AVHRR file',&
               'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
       END SELECT

    END DO

    ! if no clavr then the mask is the default of 3 (cloudy).  Otherwise
    ! subsequent processing fails (no 'clavr_there' in the IMG structure).
    IF( ALLOCATED(IMG%ClearingFlags) )THEN
       DEALLOCATE(IMG%ClearingFlags)
    ENDIF
    ALLOCATE(IMG%ClearingFlags(nelem,nline),STAT=STAT)
    IF( 0 .NE. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%CloudMask',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IF ( outData%Clavr_There ) THEN
       IMG%ClearingFlags = outData%clavr_mask
    ELSE
       IMG%ClearingFlags = 3
    END IF

#ifdef USE_HDF5
    ! if no clavr-x then the mask is the default of 7 (unprocessed).
    IF( ALLOCATED(IMG%CLAVRxMask) )THEN
       DEALLOCATE(IMG%CLAVRxMask)
    ENDIF
    ALLOCATE(IMG%CLAVRxMask(nelem,nline),STAT=STAT)
    IF( 0 .NE. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%CLAVRxMask',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%CLAVRxMask = IBITS(outData%clavrx_mask,0,8)

    ! if no clavr-x then the probability is the default of NAN_R (unprocessed).
    IF( ALLOCATED(IMG%CLAVRxPrb) )THEN
       DEALLOCATE(IMG%CLAVRxPrb)
    ENDIF
    ALLOCATE(IMG%CLAVRxPrb(nelem,nline),STAT=STAT)
    IF( 0 .NE. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%CLAVRxPrb',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    IMG%CLAVRxPrb = IBITS(outData%clavrx_prb,0,8) ! sign is preserved
#endif

    ! this maybe should be tightened up (NAN_Rs are used for the
    ! individual channels)
    if( ALLOCATED(IMG%Bad_Data) )then
       DEALLOCATE(IMG%Bad_Data)
    endif
    ALLOCATE(IMG%Bad_Data(nelem,nline),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL out_ERROR('Cannot ALLOCATE IMG%Bad_Data',&
            'Copy_to_IMG','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,'')
    ENDIF
    !
    DO line=1,nline
       IF (IMG%badTime(line) &
            .OR.IMG%badNavigation(line) &
            .OR.IMG%badCalibration(line)) THEN
          IMG%Bad_Data(:,line) = 1
       ELSE
          IMG%Bad_Data(:,line) = 0
       END IF
    END DO
  END SUBROUTINE Copy_To_IMG

  ! My version of LSD calculation from Marks GbcsImageUtil
  SUBROUTINE AVHRR_Image_Calculate_LSD(IMG, chan)

    USE GbcsPixelLoaders

    TYPE(Imagery), INTENT(INOUT) :: IMG
    INTEGER, INTENT(IN) :: chan

    ! Local variables
    INTEGER :: line, elem
    INTEGER :: numb

    DO line=1,IMG%Dims%line
       DO elem=1,IMG%Dims%elem
          IMG%Chan_LSD(chan)%d(elem,line) = &
               Calculate_LSD(IMG,chan,elem,line,numb)
       END DO
    END DO

  END SUBROUTINE AVHRR_Image_Calculate_LSD

  SUBROUTINE Parse_RecordV3(Coefs,DataType,scanLinePos,sizeOfRecords,&
       record,outData,bad_data,out_radiances,getTiePoints,getLandMask,&
       getLocation,printAllErrors,ignore_solar_contamination)

    TYPE(AVHRR_Instrument_Coefs), INTENT(IN) :: Coefs
    INTEGER, INTENT(IN) :: DataType
    INTEGER, INTENT(IN) :: scanLinePos
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER(GbcsInt1), INTENT(IN) :: record(sizeOfRecords)
    TYPE(AVHRR_Data), POINTER :: outData
    TYPE(AVHRR_Bad_Data), INTENT(IN) :: bad_data
    LOGICAL, INTENT(IN) :: out_radiances
    LOGICAL, OPTIONAL :: getTiePoints
    LOGICAL, OPTIONAL :: getLandMask
    LOGICAL, OPTIONAL :: getLocation
    LOGICAL, OPTIONAL :: printAllErrors
    LOGICAL, OPTIONAL :: ignore_solar_contamination

    INTEGER :: i,j,k
    LOGICAL :: badData
    LOGICAL :: bad_Top
    LOGICAL :: bad_Navigation
    LOGICAL :: bad_Time
    LOGICAL :: bad_Calibration
    INTEGER(GbcsInt2) :: shortVal
    INTEGER :: scanLineNumber
    INTEGER(GbcsInt4) :: qualityTop
    INTEGER(GbcsInt4) :: earthLocation
    INTEGER(GbcsInt2) :: qualityCalib(3)
    INTEGER(GbcsInt1) :: qualityCal
    INTEGER(GbcsInt1) :: qualityTime
    INTEGER(GbcsInt1) :: qualityEarth
    INTEGER(GbcsInt4) :: syncErrorCount
    INTEGER :: scanLineBitField
    LOGICAL :: timeCorrected
    INTEGER(GbcsInt2) :: clockDrift_msecs
    INTEGER :: frameDayno, frameUTC_msecs

    LOGICAL :: filter3a
    INTEGER :: utc
    INTEGER :: intVal
    INTEGER :: stat

    REAL :: calib1(2)
    REAL :: calib1_2(2)
    REAL :: calib1_intercept
    REAL :: calib2(2)
    REAL :: calib2_2(2)
    REAL :: calib2_intercept
    REAL :: calib3A(2)
    REAL :: calib3A_2(2)
    REAL :: calib3A_intercept
    REAL :: calib3(3)
    REAL :: calib4(3)
    REAL :: calib5(3)

    REAL :: satelliteAltitude
    REAL :: solarAngleIn(51)
    REAL :: satelliteAngleIn(51)
    REAL :: relativeAngleIn(51)
    REAL :: latitudeIn(51)
    REAL :: longitudeIn(51)

    REAL :: solarAngle(NDATA_MAX)
    REAL :: sateliteAngle(NDATA_MAX)
    REAL :: relativeAngle(NDATA_MAX)
    REAL :: latitude(NDATA_MAX)
    REAL :: longitude(NDATA_MAX)

    INTEGER :: prtPosition = 0

    INTEGER :: word
    INTEGER :: earthCounts(5,NDATA_MAX)
    INTEGER :: counts
    REAL :: radiance

    REAL :: minVal
    INTEGER :: minPos
    REAL :: minPosArray(4)

    REAL :: diff

    INTEGER :: nMean3
    INTEGER :: nMean4
    INTEGER :: nMean5

    REAL :: prtCounts
    REAL :: prtCounts3
    REAL :: prtCounts4
    REAL :: prtCounts5

    LOGICAL :: solar_3b
    LOGICAL :: solar_4
    LOGICAL :: solar_5
    LOGICAL :: print_All_Errors
    LOGICAL :: transition3A3B
    LOGICAL :: ignore_solar_contam

    ! Stores 10 bit words from minor frame
    INTEGER(GbcsInt2) :: minor_frame(103)

    print_All_Errors = .FALSE.
    IF( PRESENT(printAllErrors) )THEN
       print_All_Errors = printAllErrors
    ENDIF

    ignore_solar_contam = .FALSE.
    IF( PRESENT(ignore_solar_contamination) )THEN
       ignore_solar_contam = ignore_solar_contamination
    ENDIF

    if( .not. AVHRR_setup )then
       if( AVHRR_GAC .eq. DataType )then
          ndata = NDATA_GAC
          ndata_counts = 682
          do i=1,4
             prtRotate(i) = prtRotate_GAC(i)
          end do
       else if( AVHRR_LAC .eq. DataType .or. AVHRR_FRAC .eq. DataType )then
          ndata = NDATA_LAC
          ndata_counts = 3414
          do i=1,4
             prtRotate(i) = prtRotate_LAC(i)
          end do
       else
          call out_ERROR('Data not GAC/LAC',&
               'Parse_RecordV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
       AVHRR_Setup = .TRUE.
       ! Allocate output arrays
       call Allocate_OutData( ndata, outData )
       outData%badTime = .TRUE.
       outData%badNavigation = .TRUE.
       outData%badCalibration = .TRUE.
    endif

    ! check to see if we have enought storage
    if( scanLinePos .gt. outData%arraySize )then
       call Reallocate_outData( outData, MEMORY_ALLOC_STEP )
    endif

    ! assume bad to start with
    outData%badTime(scanLinePos) = .TRUE.
    outData%badNavigation(scanLinePos) = .TRUE.
    outData%badCalibration(scanLinePos) = .TRUE.

    ! Get scan line number
    scanLineNumber = header2ushort( 0, record )
    outData%scanLineNumber(scanLinePos) = scanLineNumber

    ! Check quality
    qualityTop = header2int( 24, record )
    DO i=1,103
       ! (can use header2short if just testing bits)
       minor_frame(i) = IBITS(header2ushort(1056-2+i*2,Record),0,10)
    END DO
    syncErrorCount = header2ushort(38,Record)
    qualityCalib(1) = header2short( 32, record )
    qualityCalib(2) = header2short( 34, record )
    qualityCalib(3) = header2short( 36, record )
    qualityTime = record(30)
    qualityCal = record(31)
    qualityEarth = record(32)
    earthLocation = header2int( 312, record )
    scanLineBitField = header2ushort( 12, record )
    transition3A3B=IBITS(scanLineBitField,0,2).eq.2
    outData%transition3a3b(scanLinePos) = .TRUE.
    ! this field is for info only: after late 2002 there are no corrections
    ! because the spacecraft time is correct
    timeCorrected = BTEST(scanLineBitField,14)
    ! clockDrift_msecs = spacecraft time - UTC
    clockDrift_msecs = header2short( 6, record)
    ! N15 and N16 have timeCorrected=T but the scanline time is the same as the
    ! frame time.  And N17,18,19 should have no corrections required (corrected in
    ! spacecraft) but have timeCorrected=T and scanline time is the same as the
    ! frame time.  Checking what happens on the orbit where the corrections are
    ! switched on, it is clear by looking at the lon/lat that the frame time and
    ! scanline time have both been corrected in the L1B file (even though one would
    ! not expect the frame time to be correct until after Aug 2002). This is different
    ! from v1.

    badData = .FALSE.
    bad_Top = .FALSE.
    bad_Navigation = .FALSE.
    bad_Time = .FALSE.
    bad_Calibration = .FALSE.

    ! Get time
    outData%year(scanLinePos) = header2ushort(2,record)
    outData%dayno(scanLinePos) = header2ushort(4,record)
    CALL Parse_DayNo( outData%year(scanLinePos), outData%dayno(scanLinePos),&
         outData%month(scanLinePos), outData%day(scanLinePos) )
    utc = header2int(8,record)
    outData%UTC_msecs(scanLinePos) = utc
    outData%hours(scanLinePos) = (utc*1e-3/3600.)
    outData%time(scanLinePos) = get_date_time( outData%year(scanLinePos),&
         outData%dayno(scanLinePos),outData%hours(scanLinePos), time_yearstart )

!!$    PRINT'(2i6,2b2,4i3,4b9.8,1b11.10,2b9.8)',scanLinePos,scanLineNumber,&
!!$         IBITS( earthLocation, 18, 1 ),& ! manoeuvre in progress
!!$         IBITS( earthLocation, 17, 1 ),& ! location within tolerance
!!$         IBITS( earthLocation, 12, 4 ),& ! up-to-date earth location
!!$         IBITS( earthLocation, 8, 4 ),& ! in YGC or nominal with good attitude
!!$         IBITS( earthLocation, 4, 4 ),& ! in nominal; or undocumented value
!!$         IBITS( earthLocation, 0, 4 ),& ! in nominal, no tests; or undocumented value
!!$         IBITS( qualityTop, 24, 8 ),&
!!$         IBITS( qualityTop, 16, 8 ),&
!!$         IBITS( qualityTop, 8, 8 ),&
!!$         IBITS( qualityTop, 0, 8 ),&
!!$         IBITS( minor_frame(7), 0, 10 ),&
!!$         qualityTime, qualityEarth

    bad_Top = CheckTopV2_V3(qualityTop, minor_frame, syncErrorCount)
    bad_Navigation = CheckNavigationV3(qualityTop, qualityTime, qualityEarth, &
         earthLocation, outData%AVHRR_No)
    bad_Time = CheckTimeV2_V3(qualityTop, qualityTime)
    bad_Calibration = CheckCalibrationV2_V3(qualityTop, qualityCal, solar_3b, &
         solar_4, solar_5)

    ! too bad to use at all
    IF( bad_Top )THEN
       IF( 0 .EQ. printOnce(1) .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (top/frame errors) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       ENDIF
       printOnce(1) = 1
       IF( bad_data%allow_bad_top )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_outdata_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    ! bad time means bad navigation
    if( bad_Time )then
       if( 0 .eq. printOnce(3) .or. print_All_Errors )then
          write(*,'(''Quality Problems (time) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       endif
       printOnce(3) = 1
       if( bad_data%allow_bad_time )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_outdata_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badTime(scanLinePos) = .FALSE.
    endif

    IF( bad_Time )THEN
       CALL set_time_bad( outData, scanLinePos )
    ENDIF

    ! CheckNavigationV3 also checks the time
    if( bad_Navigation )then
       if( 0 .eq. printOnce(2) .or. print_All_Errors )then
          write(*,'(''Quality Problems (navigation) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
          if( bad_data%use_bad_nav )then
             write(*,'(''  Using potentially bad navigation '')')
          endif
       endif
       printOnce(2) = 1
       if( bad_data%allow_bad_nav .or. bad_data%use_bad_nav )then
          badData = .TRUE.
          if( bad_data%use_bad_nav )then
             ! If use navigation then it is not considered 'bad'
             bad_navigation = .FALSE.
             outdata%badNavigation(scanLinePos) = .FALSE.
          endif
       else
          prtNumber = -1
          ! everything except time
          call set_nav_bad( outData, scanLinePos )
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badNavigation(scanLinePos) = .FALSE.
    endif

    if( .not. bad_Navigation )then
       ! Satellite angle in degrres
       do i=1,51
          ! Note that the start numbers are offset bu i*6 or i*8 from the
          ! record numbers
          solarAngleIn(i) = header2short((322+i*6),Record)/1e2
          satelliteAngleIn(i) = header2short((324+i*6),Record)/1e2
          relativeAngleIn(i) = header2short((326+i*6),Record)/1e2
          latitudeIn(i) = header2int((632+i*8),Record)/1e4
          longitudeIn(i) = header2int((636+i*8),Record)/1e4
       end do

       CALL InterpolateData(solarAngleIn,solarAngle,dataType)
       CALL InterpolateAngles(relativeAngleIn,satelliteAngleIn,relativeAngle,sateliteAngle,dataType)

       CALL InterpolateLonLat(longitudeIn,latitudeIn,longitude(1:ndata),latitude(1:ndata),dataType)
       outData%Lon(:,scanLinePos) = longitude(1:ndata)
       outData%Lat(:,scanLinePos) = latitude(1:ndata)
       outData%satZA(:,scanLinePos) = sateliteAngle(1:ndata)
       outData%solZA(:,scanLinePos) = solarAngle(1:ndata)
       outData%relAz(:,scanLinePos) = relativeAngle(1:ndata)
    else
       call set_nav_bad( outData, scanLinePos )
    endif

    IF (PRESENT(getTiePoints)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLandMask)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLocation)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF( bad_Calibration )THEN
       IF( 0 .EQ. printonce(4)  .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (calibration) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       printonce(4) = 1
       IF( bad_data%allow_bad_calib )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    ! Are we dealing with 3B
    shortVal = header2short(1068,record)
    if( BTEST(shortVal,0) )then
       filter3a = .TRUE.
    else
       filter3a = .FALSE.
    endif

    IF( .NOT.filter3a &
         .AND. (0 .NE. qualityCalib(1) .OR. transition3A3B) )THEN
       if( 0 .eq. printonce(5)  .or. print_All_Errors )then
          WRITE(*,'(''Quality Problems (Calib 1) or 3A/B transition, scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(5) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( 0 .ne. qualityCalib(2) )then
       if( 0 .eq. printonce(6)  .or. print_All_Errors )then
          write(*,'(''Quality Problems (Calib 2) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(6) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( 0 .ne. qualityCalib(3) )then
       if( 0 .eq. printonce(7)  .or. print_All_Errors )then
          write(*,'(''Quality Problems (Calib 3) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(7) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( BTEST(qualityCal,5) )then
       if( 0 .eq. printonce(8)  .or. print_All_Errors )then
          write(*,'(''PRTs not calibrated on scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(8) = 1
       prtNumber = -1
       call set_calib_bad( outData, scanLinePos )
       return
    endif

    if( BTEST(qualityCal,4) )then
       if( 0 .eq. printonce(9)  .or. print_All_Errors )then
          write(*,'(''PRTs marginally calibrated on scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(9) = 1
       prtNumber = -1
       call set_calib_bad( outData, scanLinePos )
       return
    endif

    if( .not. badData )then
       badDataLoop: do i=1,10
          if( 1 .eq. printOnce(i) )then
             printOnce = (/0,0,0,0,0,0,0,0,0,0/)
             write(*,'(''End of quality problems on scan line, position'',2i7)')&
                  scanLineNumber, scanLinePos
             exit badDataLoop
          endif
       end do badDataLoop
    endif

    ! Solar flags
    outData%orig_solar_contamination_3b(scanLinePos) = solar_3b
    outData%orig_solar_contamination_4(scanLinePos) = solar_4
    outData%orig_solar_contamination_5(scanLinePos) = solar_5

    ! Get calibration values
    calib1(2) = header2int(48,Record)/1e7
    calib1(1) = header2int(52,Record)/1e6
    calib1_2(2) = header2int(56,Record)/1e7
    calib1_2(1) = header2int(60,Record)/1e6
    calib1_intercept = header2int(64,Record)
    ! Channel 2
    calib2(2) = header2int(108,Record)/1e7
    calib2(1) = header2int(112,Record)/1e6
    calib2_2(2) = header2int(116,Record)/1e7
    calib2_2(1) = header2int(120,Record)/1e6
    calib2_intercept = header2int(124,Record)
    ! Infrared
    if( .not. filter3a )then
       calib3(1) = header2int(228,Record)/1e6
       calib3(2) = header2int(232,Record)/1e6
       calib3(3) = header2int(236,Record)/1e6
       ! Channel 3A
       calib3A(2) = 0.
       calib3A(1) = 0.
       calib3A_2(2) = 0.
       calib3A_2(1) = 0.
       calib3A_intercept = 0
    else
       calib3(1) = 0.
       calib3(2) = 0.
       calib3(3) = 0.
       ! Channel 3A
       calib3A(2) = header2int(168,Record)/1e7
       calib3A(1) = header2int(172,Record)/1e6
       calib3A_2(2) = header2int(176,Record)/1e7
       calib3A_2(1) = header2int(180,Record)/1e6
       calib3A_intercept = header2int(184,Record)
    endif
    calib4(1) = header2int(252,Record)/1e6
    calib4(2) = header2int(256,Record)/1e6
    calib4(3) = header2int(260,Record)/1e7
    calib5(1) = header2int(276,Record)/1e6
    calib5(2) = header2int(280,Record)/1e6
    calib5(3) = header2int(284,Record)/1e7

    outData%calib1(:,scanLinePos) = calib1
    outData%calib2(:,scanLinePos) = calib2
    outData%calib1_2(:,scanLinePos) = calib1_2
    outData%calib2_2(:,scanLinePos) = calib2_2
    outData%calib1_intercept(scanLinePos) = calib1_intercept
    outData%calib2_intercept(scanLinePos) = calib2_intercept
    outData%calib3A(:,scanLinePos) = calib3A
    outData%calib3A_2(:,scanLinePos) = calib3A_2
    outData%calib3A_intercept(scanLinePos) = calib3A_intercept
    outData%calib3(:,scanLinePos) = calib3
    outData%calib4(:,scanLinePos) = calib4
    outData%calib5(:,scanLinePos) = calib5

!  roll = header2short(320,Record)/1e3
!  pitch = header2short(322,Record)/1e3
!  yaw = header2short(324,Record)/1e3

    satelliteAltitude = header2ushort(326,Record)/10.
    outData%satelliteAlt_There = .TRUE.
    outData%satelliteAltitude(scanLinePos) = satelliteAltitude

    if( AVHRR_GAC .eq. dataType )then
       outData%radiator(scanLinePos) = convertTemperatureV3(record(4024),&
            Coefs%temperatureCoefs1,0)
       outData%electronics(scanLinePos) = convertTemperatureV3(record(4032),&
            Coefs%temperatureCoefs2,1)
       outData%cooler(scanLinePos) = convertTemperatureV3(record(4033),&
            Coefs%temperatureCoefs3,1)
       outData%baseplate(scanLinePos) = convertTemperatureV3(record(4034),&
            Coefs%temperatureCoefs4,1)
       outData%motor(scanLinePos) = convertTemperatureV3(record(4035),&
            Coefs%temperatureCoefs5,1)
       outData%a_d_conv(scanLinePos) = convertTemperatureV3(record(4036),&
            Coefs%temperatureCoefs6,1)
       outData%patch(scanLinePos) = convertTemperatureV3(record(4022),&
            Coefs%patchCoef,0)
       outData%patchExtended(scanLinePos) = convertTemperatureV3(record(4023),&
            Coefs%PatchCoefExt,0)

       ! Get Motor current
       outData%motorCurrent(scanLinePos) = convertCurrentV3(record(4030),&
            Coefs%motorCurrentCoefs)

       ! Get CLAVR flags
       intVal = header2int(4048,record)
       if( btest(intVal, 0) )then
          outData%Clavr_There = .TRUE.
          CALL getClavrFlags( AVHRR_GAC, outData%clavr_mask(:,scanLinePos), &
               4056, record )
       else
          outData%clavr_mask(:,scanLinePos) = 3 ! cloudy
       endif
    else if( AVHRR_LAC .eq. dataType )then
       outData%radiator(scanLinePos) = convertTemperatureV3(record(14952),&
            Coefs%temperatureCoefs1,0)
       outData%electronics(scanLinePos) = convertTemperatureV3(record(14960),&
            Coefs%temperatureCoefs2,1)
       outData%cooler(scanLinePos) = convertTemperatureV3(record(14961),&
            Coefs%temperatureCoefs3,1)
       outData%baseplate(scanLinePos) = convertTemperatureV3(record(14962),&
            Coefs%temperatureCoefs4,1)
       outData%motor(scanLinePos) = convertTemperatureV3(record(14963),&
            Coefs%temperatureCoefs5,1)
       outData%a_d_conv(scanLinePos) = convertTemperatureV3(record(14964),&
            Coefs%temperatureCoefs6,1)
       outData%patch(scanLinePos) = convertTemperatureV3(record(14949),&
            Coefs%patchCoef,0)
       outData%patchExtended(scanLinePos) = convertTemperatureV3(record(14950),&
            Coefs%PatchCoefExt,0)

       ! Get motor current
       outData%motorCurrent(scanLinePos) = convertCurrentV3(record(14958),&
            Coefs%motorCurrentCoefs)

       ! Get CLAVR flags
       intVal = header2int(14976,record)
       if( btest(intVal, 0) )then
          outData%Clavr_There = .TRUE.
          CALL getClavrFlags( AVHRR_LAC, outData%clavr_mask(:,scanLinePos), &
               14984, record )
       else
          outData%clavr_mask(:,scanLinePos) = 3 ! cloudy
       endif
    endif

    prtCounts3 = header2ushort(1090,record)
    prtCounts4 = header2ushort(1092,record)
    prtCounts5 = header2ushort(1094,record)
!---Modif du 12 nov---
!on veut imprimer les prt counts pour implementer la methode de Trishchenko
!    outData%prtcounts(1,scanLinePos)=prtCounts3
!    outData%prtcounts(2,scanLinePos)=prtCounts4
!    outData%prtcounts(3,scanLinePos)=prtCounts5
    !print*,outData%prtcounts(3,scanLinePos)
!------------------- 
    ! Need to get which PRT we are at
    if( -1 .eq. prtNumber )then
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
       endif
       storedPrtTemp = 0.
    else
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtCountsPrevStore(i) = prtCountsStore(i)
             else
                prtCountsPrevStore(i) = 0
             endif
          end do
       else
          prtCounts = (prtCounts3 + prtCounts4 + prtCounts5)/3.
          if( 1 .gt. prtNumber .or. 4 .lt. prtNumber )then
             ! Rejig array

             ! To be honest I can't remember why I wrote the following code
             ! should never really be used...

             if( 0 .lt. prtCountsPrevStore(1) .and. &
                  0 .lt. prtCountsPrevStore(2) .and. &
                  0 .lt. prtCountsPrevStore(3) .and. &
                  0 .lt.  prtCountsPrevStore(4) )then
                minPos = 0
                minVal = 1e32
                do j=1,4
                   diff = abs(prtCounts - prtCountsPrevStore(j))
                   if( minVal > diff )then
                      minVal = diff
                      minPos = j
                   endif
                end do
                minPosArray(4) = minPos
                do i=1,4
                   minPos = 0
                   minVal = 1e32
                   do j=1,4
                      diff = abs(prtCountsStore(i) - prtCountsPrevStore(j))
                      if( minVal > diff )then
                         minVal = diff
                         minPos = j
                      endif
                   end do
                   minPosArray(i) = minPos
                end do

                ! Now have array of indexes to matching entries - remap these
                ! ones
                ! Find repeat
                minPos = -1
                do i=2,4
                   if( minPosArray(i-1) .eq. minPosArray(i) )then
                      minPos = i
                   endif
                end do

                if( -1 .eq. minPos )then
                   ! Assume line is bad (for some reason) and set bad
                   prtNumber = -1
                   call set_calib_bad( outData, scanLinePos )
                   return
                else

                   ! Shift down
                   do i=minPos,3
                      minPosArray(i) = minPosArray(i+1)
                   end do

                   ! Set extra value to last value
                   minPosArray(4) = prtCounts
                endif

                write(*,*)'Resetting prtCounts array'

             else
                write(*,*)'No previous prt record - just taking first 4'
                prtNumber=1
             endif

          else
             prtCountsStore(prtNumber) = prtCounts
             prtCountsStoreAll(1,prtNumber) = prtCounts3
             prtCountsStoreAll(2,prtNumber) = prtCounts4
             prtCountsStoreAll(3,prtNumber) = prtCounts5
             prtNumber = prtNumber + 1
          endif

          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtPosition = prtRotate(i)
                if( 1 .gt. prtPosition .or. 4 .lt. prtPosition )then
                   call out_ERROR('prtPosition out of range',&
                        'Parse_RecordV3','NOAA_LoadAVHRRLevel1B.f90',&
                        .FALSE.,' ')
                endif
                storedPrtCounts(prtPosition) = prtCountsStore(i)
                storedPrtCountsAll(:,prtPosition) = prtCountsStoreAll(:,i)
                storedPrtTemp(prtPosition) = &
                     Coefs%prtTempCoefs1(prtPosition) + &
                     Coefs%prtTempCoefs2(prtPosition)*prtCountsStore(i) + &
                     Coefs%prtTempCoefs3(prtPosition)*(prtCountsStore(i)**2) + &
                     Coefs%prtTempCoefs4(prtPosition)*(prtCountsStore(i)**3) + &
                     Coefs%prtTempCoefs5(prtPosition)*(prtCountsStore(i)**4)
             endif
          end do
       end if
    end if

    if( storedPrtTemp(1) .gt. 200 .and. &
         storedPrtTemp(2) .gt. 200 .and. &
         storedPrtTemp(3) .gt. 200 .and. &
         storedPrtTemp(4) .gt. 200 )then
       outData%prt1CountsAll(:,scanLinePos) = storedPrtCountsAll(:,1)
       outData%prt2CountsAll(:,scanLinePos) = storedPrtCountsAll(:,2)
       outData%prt3CountsAll(:,scanLinePos) = storedPrtCountsAll(:,3)
       outData%prt4CountsAll(:,scanLinePos) = storedPrtCountsAll(:,4)
       outData%prt1Counts(scanLinePos) = storedPrtCounts(1)
       outData%prt2Counts(scanLinePos) = storedPrtCounts(2)
       outData%prt3Counts(scanLinePos) = storedPrtCounts(3)
       outData%prt4Counts(scanLinePos) = storedPrtCounts(4)
       outData%prt1(scanLinePos) = storedPrtTemp(1)
       outData%prt2(scanLinePos) = storedPrtTemp(2)
       outData%prt3(scanLinePos) = storedPrtTemp(3)
       outData%prt4(scanLinePos) = storedPrtTemp(4)
    else
       outData%prt1CountsAll(:,scanLinePos) = NAN_R
       outData%prt2CountsAll(:,scanLinePos) = NAN_R
       outData%prt3CountsAll(:,scanLinePos) = NAN_R
       outData%prt4CountsAll(:,scanLinePos) = NAN_R
       outData%prt1Counts(scanLinePos) = NAN_R
       outData%prt2Counts(scanLinePos) = NAN_R
       outData%prt3Counts(scanLinePos) = NAN_R
       outData%prt4Counts(scanLinePos) = NAN_R
       outData%prt1(scanLinePos) = NAN_R
       outData%prt2(scanLinePos) = NAN_R
       outData%prt3(scanLinePos) = NAN_R
       outData%prt4(scanLinePos) = NAN_R
    endif

    ! Get BB counts
    outData%bbodyFilter3(:,scanLinePos) = 0
    outData%bbodyFilter4(:,scanLinePos) = 0
    outData%bbodyFilter5(:,scanLinePos) = 0
    outData%bb3(scanLinePos) = 0.
    outData%bb4(scanLinePos) = 0.
    outData%bb5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       if( .not.filter3a )then
          outData%bbodyFilter3(i,scanLinePos) = &
               header2ushort((1100+(i-1)*6),record)
       endif
       outData%bbodyFilter4(i,scanLinePos) = &
            header2ushort((1102+(i-1)*6),record)
       outData%bbodyFilter5(i,scanLinePos) = &
            header2ushort((1104+(i-1)*6),record)
       ! in fact outData%bbodyFilter? are always +ve
       if( .not. filter3a )then
          if( 0 .lt. outData%bbodyFilter3(i,scanLinePos) )then
             outData%bb3(scanLinePos) = outData%bb3(scanLinePos) + &
                  outData%bbodyFilter3(i,scanLinePos)
             nMean3 = nMean3 + 1
          endif
       endif
       if( 0 .lt. outData%bbodyFilter4(i,scanLinePos) )then
          outData%bb4(scanLinePos) = outData%bb4(scanLinePos) + &
               outData%bbodyFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 0 .lt. outData%bbodyFilter5(i,scanLinePos) )then
          outData%bb5(scanLinePos) = outData%bb5(scanLinePos) + &
               outData%bbodyFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( .not. filter3a )then
       if( NBBODY_SAMPLES .ne. nMean3 )then
          write(*,*)'Dropped Black-body (3.7) counts'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       endif
       if( 0 .eq. nMean3 )then
          outData%bb3(scanLinePos) = NAN_R
          write(*,*)'Black-body counts (3.7) are zero'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       else
          outData%bb3(scanLinePos) = outData%bb3(scanLinePos) / nMean3
       end if
    else
       outData%bb3(scanLinePos) = NAN_R
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Black-body (11) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean4 )then
       outData%bb4(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb4(scanLinePos) = outData%bb4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Black-body (12) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean5 )then
       outData%bb5(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb5(scanLinePos) = outData%bb5(scanLinePos) / nMean5
    end if

    ! Get space counts
    outData%spaceFilter1(:,scanLinePos) = 0
    outData%spaceFilter2(:,scanLinePos) = 0
    outData%spaceFilter3a(:,scanLinePos) = 0
    outData%spaceFilter3(:,scanLinePos) = 0
    outData%spaceFilter4(:,scanLinePos) = 0
    outData%spaceFilter5(:,scanLinePos) = 0
    outData%sp3(scanLinePos) = 0.
    outData%sp4(scanLinePos) = 0.
    outData%sp5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       outData%spaceFilter1(i,scanLinePos) = &
            header2ushort((1160+(i-1)*10),record)
       outData%spaceFilter2(i,scanLinePos) = &
            header2ushort((1162+(i-1)*10),record)
       if( .not.filter3a )then
          outData%spaceFilter3(i,scanLinePos) = &
               header2ushort((1164+(i-1)*10),record)
          outData%spaceFilter3a(i,scanLinePos) = NAN_I
       else
          outData%spaceFilter3a(i,scanLinePos) = &
               header2ushort((1164+(i-1)*10),record)
          outData%spaceFilter3(i,scanLinePos) = NAN_I
       endif
       outData%spaceFilter4(i,scanLinePos) = &
            header2ushort((1166+(i-1)*10),record)
       outData%spaceFilter5(i,scanLinePos) = &
            header2ushort((1168+(i-1)*10),record)
       ! in fact outData%spaceFilter? are always +ve
       if( .not. filter3a )then
          if( 900 .lt. outData%spaceFilter3(i,scanLinePos) )then
             outData%sp3(scanLinePos) = outData%sp3(scanLinePos) + &
                  outData%spaceFilter3(i,scanLinePos)
             nMean3 = nMean3 + 1
          endif
       endif
       if( 900 .lt. outData%spaceFilter4(i,scanLinePos) )then
          outData%sp4(scanLinePos) = outData%sp4(scanLinePos) + &
               outData%spaceFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 900 .lt. outData%spaceFilter5(i,scanLinePos) )then
          outData%sp5(scanLinePos) = outData%sp5(scanLinePos) + &
               outData%spaceFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( .not. filter3a )then
       if( NBBODY_SAMPLES .ne. nMean3 )then
          write(*,*)'Dropped Space (3.7) counts'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       endif
       if( 0 .eq. nMean3 )then
          outData%sp3(scanLinePos) = NAN_R
          write(*,*)'Space counts (3.7) are zero'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       else
          outData%sp3(scanLinePos) = outData%sp3(scanLinePos) / nMean3
       end if
    else
       outData%sp3(scanLinePos) = NAN_R
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Space (11) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean4 )then
       outData%sp4(scanLinePos) = NAN_R
       write(*,*)'Space counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp4(scanLinePos) = outData%sp4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Space (12) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean5 )then
       outData%sp5(scanLinePos) = NAN_R
       write(*,*)'Space counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp5(scanLinePos) = outData%sp5(scanLinePos) / nMean5
    end if

    ! L1B file calibration is OK but if ignoring solar contam data, ignore
    IF( ignore_solar_contam .and. (solar_3b .or. solar_4 .or. solar_5) )THEN
       outData%badCalibration(scanLinePos) = .TRUE.
    ELSE
       outData%badCalibration(scanLinePos) = .FALSE.
    ENDIF

    if( .not. bad_Navigation .and. &
         .not.outData%badCalibration(scanLinePos) )then

       ! Get earth observations
       j = 0
       k = 0
       do i=0,ndata_counts-1
          ! Data is interlaced
          word = header2int((1264+i*4),record)
          ! Get bits 20-29
          earthCounts(k+1,j+1) = IBITS(word,20,10)
          k = k+1
          if( 5 .eq. k)then
             j = j + 1
             k = 0
          end if
          if( 682 .eq. ndata_counts )then
             ! Get bits 10-19
             earthCounts(k+1,j+1) = IBITS(word,10,10)
             k = k+1
             if( 5 .eq. k)then
                j = j + 1
                k = 0
             end if
             ! Get bits 0-9 if i != 681
             if( ndata_counts-1 .ne. i )then
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          else
             ! Get bits 10-19 && 0-9 if i != 3413
             if( ndata_counts-1 .ne. i )then
                ! Get bits 10-19
                earthCounts(k+1,j+1) = IBITS(word,10,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
                ! Get bits 0-9
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k+1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          end if
       end do

       if( 0 .ne. k .and. j .ne. ndata )then
          call out_ERROR('Problems parsing earth observations',&
               'Parse_RecordV3','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif

       do i=1,ndata
          ! Visible
          ! Channel 1
          counts = earthCounts(1,i)
          outData%counts1(i,scanLinePos) = counts
          if( counts .ge. calib1_intercept )then
             outData%array1(i,scanLinePos) = calib1_2(1) + &
                  calib1_2(2) * counts
          else
             outData%array1(i,scanLinePos) = calib1(1) + &
                  calib1(2) * counts
          endif
          if( 0 .gt. outData%array1(i,scanLinePos) )then
             outData%array1(i,scanLinePos) = 0.
          endif
          outData%array1(i,scanLinePos) &
               = outData%array1(i,scanLinePos) / 100.0 ! percent to fraction

          ! Channel 2
          counts = earthCounts(2,i)
          outData%counts2(i,scanLinePos) = counts
          if( counts .ge. calib2_intercept )then
             outData%array2(i,scanLinePos) = calib2_2(1) + &
                  calib2_2(2) * counts
          else
             outData%array2(i,scanLinePos) = calib2(1) + &
                  calib2(2) * counts
          endif
          if( 0 .gt. outData%array2(i,scanLinePos) )then
             outData%array2(i,scanLinePos) = 0.
          endif
          outData%array2(i,scanLinePos) &
               = outData%array2(i,scanLinePos) / 100.0 ! percent to fraction

          ! Channel 3 - either 3A/3B
          IF( filter3a )THEN

             IF( transition3A3B )THEN
                outData%array3A(i,scanLinePos) = NAN_R
                outData%transition3a3b(scanLinePos) = .TRUE.
             ELSE
                outData%transition3a3b(scanLinePos) = .FALSE.
                counts = earthCounts(3,i)
                IF( counts .GE. calib3A_intercept )THEN
                   outData%array3A(i,scanLinePos) = calib3A_2(1) + &
                        calib3A_2(2) * counts
                ELSE
                   outData%array3A(i,scanLinePos) = calib3A(1) + &
                        calib3A(2) * counts
                ENDIF
                IF( 0 .GT. outData%array3A(i,scanLinePos) )THEN
                   outData%array3A(i,scanLinePos) = 0.
                ENDIF
                outData%array3A(i,scanLinePos) &
                     = outData%array3A(i,scanLinePos) / 100.0 ! percent to fraction
             END IF

             ! it would be better to have a outData%filter3a
             ! flag but counts3 is *currently*) only used
             ! for the new calibration of 3B
             outData%counts3(i,scanLinePos) = NAN_R
             outData%array3B(i,scanLinePos) = NAN_R

          else

             IF( transition3A3B )THEN
                outData%array3B(i,scanLinePos) = NAN_R
                outData%counts3(i,scanLinePos) = NAN_R
                outData%transition3a3b(scanLinePos) = .FALSE.
             ELSE
                outData%transition3a3b(scanLinePos) = .TRUE.
                counts = earthCounts(3,i)
                outData%counts3(i,scanLinePos) = counts
                radiance = calib3(1) + calib3(2)*counts + &
                     calib3(3)*counts*counts
                IF( out_radiances )THEN
                   outData%array3B(i,scanLinePos) = radiance
                ELSE
                   outData%array3B(i,scanLinePos) = convertBT( radiance, &
                        Coefs%nuC(1), Coefs%Aval(1), Coefs%Bval(1) )
#ifdef USE_IEEE
                   IF( .NOT.ieee_is_finite(outData%array3B(i,scanLinePos)) )THEN
                      outData%array3B(i,scanLinePos) = NAN_R
                   END IF
#elif defined(USE_G95)
                   IF( isNAN(outData%array3B(i,scanLinePos)) )THEN
                      outData%array3B(i,scanLinePos) = NAN_R
                   END IF
#endif
                ENDIF
             END IF

             outData%array3A(i,scanLinePos) = NAN_R

          endif

          counts = earthCounts(4,i)
          outData%counts4(i,scanLinePos) = counts
          radiance = calib4(1) + calib4(2)*counts + &
               calib4(3)*counts*counts
          if( out_radiances )then
             outData%array4(i,scanLinePos) = radiance
          else
             outData%array4(i,scanLinePos) = convertBT( radiance, &
                  Coefs%nuC(2), Coefs%Aval(2), Coefs%Bval(2) )
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#endif
          endif

          counts = earthCounts(5,i)
          outData%counts5(i,scanLinePos) = counts
          radiance = calib5(1) + calib5(2)*counts + &
               calib5(3)*counts*counts
          if( out_radiances )then
             outData%array5(i,scanLinePos) = radiance
          else
             outData%array5(i,scanLinePos) = convertBT( radiance, &
                  Coefs%nuC(3), Coefs%Aval(3), Coefs%Bval(3) )
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array5(i,scanLinePos)) )THEN
                outData%array5(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array5(i,scanLinePos)) )THEN
                outData%array5(i,scanLinePos) = NAN_R
             END IF
#endif
          endif

       end do

    else

       ! set Earth observations to NAN because we don't know where we are
       outData%counts3(:,scanLinePos) = NAN_R
       outData%counts4(:,scanLinePos) = NAN_R
       outData%counts5(:,scanLinePos) = NAN_R
       outData%array1(:,scanLinePos) = NAN_R
       outData%array2(:,scanLinePos) = NAN_R
       outData%array3A(:,scanLinePos) = NAN_R
       outData%array3B(:,scanLinePos) = NAN_R
       outData%array4(:,scanLinePos) = NAN_R
       outData%array5(:,scanLinePos) = NAN_R

    endif
    outData%dataFilled = .TRUE.

  END SUBROUTINE Parse_RecordV3

  REAL FUNCTION convertBT( radiance, nuC, Aval, Bval )

    REAL, INTENT(IN) :: radiance
    REAL, INTENT(IN) :: nuC
    REAL, INTENT(IN) :: Aval
    REAL, INTENT(IN) :: Bval

    REAL :: constant

    IF( radiance .gt. 0 )THEN
       constant = 1. + (C1*nuC*nuC*nuC/radiance)
       constant = C2*nuC/log(constant)
       convertBT = (aVal+bVal*constant)
    ELSE
       convertBT = NAN_R
    ENDIF
    return

  END FUNCTION convertBT

  ! Use look up table radaince->temp conversion
  REAL FUNCTION convertBTold( radiance, numT, minT, maxT, nuC )

    REAL, INTENT(IN) :: radiance
    INTEGER(GbcsInt1), INTENT(IN) :: numT
    REAL(GbcsReal), INTENT(IN) :: minT(:)
    REAL(GbcsReal), INTENT(IN) :: maxT(:)
    REAL(GbcsReal), INTENT(IN) :: nuC(:)

    INTEGER :: I
    REAL :: constant
    REAL :: test_temp
    REAL :: nuC_Value

    IF( radiance .gt. 0 )THEN
       ! Work out which entry we want
       ! Assume SST temperature - 3rd entry (always exists)    
       constant = 1. + (C1*nuC(3)*nuC(3)*nuC(3)/radiance)
       constant = C2*nuC(3)/log(constant)
#ifdef USE_IEEE
       IF( .NOT.ieee_is_finite(constant) )THEN
          convertBTold = ieee_value(NAN_R,ieee_quiet_nan) ! NaN
          RETURN
       END IF
#elif defined(USE_G95)
       IF( isNAN(constant) )THEN
          convertBTold = NAN_R
          RETURN
       END IF
#endif
    ELSE
       convertBTold = NAN_R
       RETURN
    ENDIF

    nuC_Value = 0
    Loop: DO I=1,numT
       IF( constant .ge. minT(I) .and. constant .le. maxT(I) )THEN
          nuC_Value = nuC(I)
          EXIT Loop
       ENDIF
    END DO Loop
    IF( 0 .eq. nuC_Value )THEN
       ! This deals with the n/a for temps below 225 for 3.7 channel AVHRR-7
       IF( 0 .lt. nuC(1) )THEN
          IF( constant .lt. minT(1) )THEN
             nuC_Value = nuC(1)
          ENDIF
       ELSE
          IF( constant .lt. minT(2) )THEN
             nuC_Value = nuC(2)
          ENDIF
       ENDIF
       IF( constant .gt. maxT(numT) )THEN
          nuC_Value = nuC(numT)
       ENDIF
    ENDIF
    constant = 1. + (C1*nuC_Value*nuC_Value*nuC_Value/radiance)
    convertBTold = C2*nuC_Value/LOG(constant) ! might not be finite

    return

  END FUNCTION convertBTold

  REAL FUNCTION convertRadiance( temp, nuC, Aval, Bval )

    REAL, INTENT(IN) :: temp
    REAL, INTENT(IN) :: nuC
    REAL, INTENT(IN) :: Aval
    REAL, INTENT(IN) :: Bval

    REAL :: newtemp

    newtemp = (temp - aVal)/bval
    convertRadiance = C1*nuC*nuC*nuC/(exp(C2*nuC/newtemp) - 1.)

    return

  END FUNCTION convertRadiance

  SUBROUTINE setupScanAngles( data_type, AVHRR_No )

    INTEGER, INTENT(IN) :: data_type
    INTEGER, INTENT(IN) :: AVHRR_No

    INTEGER :: i = 0
    INTEGER :: j = 0
    INTEGER :: k = 0
    INTEGER :: l = 0

    INTEGER :: pos = 0

    REAL :: angle = 0.
    ! from http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm
    ! the start angle = 55.37 ~= 0.0541*1024 - 0.0541/2, i.e. angle from nadir
    ! to centre of edge LAC FOV.  For NOAA-16 start angle = 55.25 but the actual
    ! value just shifts the whole interpolation, so OK.  Define (in this code
    ! angle direction as righthand around direction of spacecraft motion, with
    ! zero at nadir.  Thus centre of first LAC of a scanline is at -55.37,
    ! here defined as -(1024-0.5)*0.0541
    REAL,PARAMETER :: stepAngle = 0.0541*Deg2Rad
    REAL :: startAngle = -(1024-0.5)*stepAngle

    if( .not. scanAngleSetup )then
       ! set up angles of the location points.  Locations for GAC and LAC
       ! are the same
       DO k=1,51
          i=(k-1)*40+25
          scanAngleIn(k)=startAngle+(i-1)*stepAngle
       END DO

       ! set up angles of the data points
       IF(AVHRR_GAC.EQ.data_type)THEN
          ! Setup scan angles
          DO j=1,409 ! = NDATA_GAC
             ! NOAA docs and email from Jeff Robel Fri, 15 Jan 2010
             ! 14:58:56 -0500 says to use the latest docs, which specifically
             ! say the GAC pixels start at the 1st LAC pixel in the scan line.
             ! I will need to check this for the pre-KLM AVHRRs.  I have checked
             ! this for Metop-2 and NOAA-18 and using 1st matches FRAC with GAC.
             ! So: GAC samples over 5 LAC's starting with the 1st and skipping
             ! the 5th, averaging over the 4 remaining LAC's.

             ! GAC gap lac index
             i=5*j
             IF (1.LE.AVHRR_No .AND. AVHRR_No.LE.14) THEN
                ! pre-KLM
                ! GAC centre is 0.5 LACs before the gap
                scanAngle(j)=startAngle+(i-1-0.5)*stepAngle
             ELSE ! Metops are -ve numbered here
                ! KLM
                ! GAC centre is 2.5 LACs before the gap
                scanAngle(j)=startAngle+(i-1-2.5)*stepAngle
             END IF
          END DO
       ELSE
          ! Setup scan angles for LAC
          DO i=1,2048
             scanAngle(i)=startAngle+(i-1)*stepAngle
          END DO
       END IF
    end if
    scanAngleSetup = .TRUE.

  end SUBROUTINE setupScanAngles

  SUBROUTINE InterpolateLonLat(longitudeIn, latitudeIn, longitude, latitude, dataType)

    REAL, DIMENSION(51), INTENT(IN) :: longitudeIn
    REAL, DIMENSION(51), INTENT(IN) :: latitudeIn
    REAL, DIMENSION(:), INTENT(OUT) :: longitude
    REAL, DIMENSION(:), INTENT(OUT) :: latitude
    INTEGER, INTENT(IN) :: dataType

    REAL, DIMENSION(51) :: thetaIn, phiIn, xIn, yIn, zIn
    REAL, DIMENSION(SIZE(longitude)) :: theta, phi, x, y, z

    thetaIn=(90.0-latitudeIn)*Deg2Rad
    phiIn=longitudeIn*Deg2Rad
    xIn=COS(phiIn)*SIN(thetaIn)
    yIn=SIN(phiIn)*SIN(thetaIn)
    zIn=COS(thetaIn)
    CALL InterpolateData(xIn,x,dataType)
    CALL InterpolateData(yIn,y,dataType)
    CALL InterpolateData(zIn,z,dataType)
    ! 1st arg is |sin(theta)|, and 0<=theta<=pi so sin(theta)>=0 OK
    ! should never have both arguments = 0
    ! cos(theta) = 0 => theta = pi/2 => sin(theta) = 1.0 and v.v.
    ! cos(phi)**2+sin(phi)**2 is never zero
    theta=ATAN2(SQRT(x**2+y**2),z)
    IF(ANY(x.EQ.0.0 .AND. y.EQ.0.0)) PRINT*,'x=y=0'
    WHERE(x.EQ.0.0 .AND. y.EQ.0.0)
       phi=0.0
    ELSEWHERE
       ! in ifort, atan2(0.0,0.0) = 0.0
       phi=ATAN2(y,x)
    END WHERE
    longitude=Rad2Deg*phi
    latitude=90.0-Rad2Deg*theta

    ! Make sure that latitude is 90S-90N
    IF(ANY(latitude.LT.-90.0)) PRINT*,'l<-90'
    IF(ANY(latitude.GT.90.0)) PRINT*,'l>90'
    WHERE(latitude.LT.-90.0) latitude=-90.0
    WHERE(latitude.GT.90.0) latitude=90.0

    ! Make sure that longitude is 0-360 E
    WHERE(longitude.LT.0.0) longitude = 360.0 + longitude

  END SUBROUTINE InterpolateLonLat

  SUBROUTINE InterpolateAngles(azimuthIn, zenithIn, azimuth, zenith, dataType)

    REAL, DIMENSION(51), INTENT(IN) :: azimuthIn
    REAL, DIMENSION(51), INTENT(IN) :: zenithIn
    REAL, DIMENSION(:), INTENT(OUT) :: azimuth
    REAL, DIMENSION(:), INTENT(OUT) :: zenith
    INTEGER, INTENT(IN) :: dataType

    REAL, DIMENSION(51) :: thetaIn, phiIn, xIn, yIn, zIn
    REAL, DIMENSION(SIZE(azimuth)) :: theta, phi, x, y, z

    thetaIn=zenithIn*Deg2Rad
    phiIn=azimuthIn*Deg2Rad
    xIn=COS(phiIn)*SIN(thetaIn)
    yIn=SIN(phiIn)*SIN(thetaIn)
    zIn=COS(thetaIn)
    CALL InterpolateData(xIn,x,dataType)
    CALL InterpolateData(yIn,y,dataType)
    CALL InterpolateData(zIn,z,dataType)

    theta=ACOS(z)
    phi=ATAN2(y,x)

    azimuth=Rad2Deg*phi
    zenith=Rad2Deg*theta

  END SUBROUTINE InterpolateAngles

  SUBROUTINE InterpolateData(AngleIn, Angle, dataType)

    REAL, INTENT(INOUT) :: AngleIn(51)
    REAL, INTENT(OUT) :: Angle(:)
    INTEGER, INTENT(IN) :: dataType

    INTEGER :: I

    if( .not.scanAngleSetup )then
       call out_ERROR('scan angle arrays are not setup',&
            'setupScanAngles',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    endif

    ! Check the all numbers are between -180 and +180
    do i=1,51
       if( -180.0 .gt. AngleIn(i) )then
          AngleIn(i) = AngleIn(i) + 360.
       else if( 180.0 .lt. AngleIn(i) )then
          AngleIn(i) = AngleIn(i) - 360.
       endif
    end do

    if( AVHRR_GAC .eq. dataType )then
       do i=1,NDATA_GAC
          Angle(I) = doInterp_GAC(i,AngleIn)
       end do
    else if(AVHRR_LAC .eq. dataType)then
       do i=1,NDATA_LAC
          Angle(I) = doInterp_LAC(i,AngleIn)
       end do
    endif

  END SUBROUTINE InterpolateData

  REAL FUNCTION doInterp_GAC(pos,arrayIn)

    INTEGER, INTENT(IN) :: pos
    REAL, INTENT(IN) :: arrayIn(:)
    INTEGER :: offset = 0

    ! Extrapolate one end and the other
    ! If in middle place data in middle bound
    if( pos .lt. 21 )then
       offset = 0
    else if(pos .gt. 389)then
       offset = 46
    else
       offset = (pos-5)/8 - 2
    endif

    doInterp_GAC = doInterp(pos,arrayIn,offset)
    return

  END FUNCTION doInterp_GAC

  REAL FUNCTION doInterp_LAC(pos,arrayIn)

    INTEGER, INTENT(IN) :: pos
    REAL, INTENT(IN) :: arrayIn(:)
    INTEGER :: offset = 0

    ! Extrapolate one end and the other
    ! If in middle place data in middle bound
    if( pos .lt. 105 )then
       offset = 0
    else if(pos .gt. 1945)then
       offset = 46
    else
       offset = (pos-25)/40 - 2
    endif

    doInterp_LAC = doInterp(pos,arrayIn,offset)
    return

  END FUNCTION doInterp_LAC

  REAL FUNCTION doInterp(pos,arrayIn,offset)

    INTEGER, INTENT(IN) :: pos
    REAL, INTENT(IN) :: arrayIn(:)
    INTEGER, INTENT(IN) :: offset

    INTEGER :: i = 0
    INTEGER :: j = 0

    REAL :: Lvalue = 0.
    REAL :: Sum = 0.

    REAL :: five_point(5)

    INTEGER :: positive = 0
    INTEGER :: date_line = 0

    ! Offset is value to set start of 5 point array to
    ! Try and match pos to middle point

    ! Check to see if we need to interpolate */
    do i=1,5
       if( scanAngle(pos) .eq. scanAngleIn(i+offset) )then
          doInterp = arrayIn(i+offset)
          return
       endif
    end do

    ! Copy over to array that can be modified
    date_line = 0
    positive = 0
    five_point = 0
    five_point(1) = ArrayIn(offset+1)
    do i=2,5
       five_point(i) = ArrayIn(i+offset)
       ! Check to see if we're crossing the date line
       if( five_point(i-1) - five_point(i) .gt. 180 )then
          positive = 1
       else if( five_point(i-1) - five_point(i) .lt. -180 )then
          positive = -1
       endif

       if( 1 .eq. positive )then
          five_point(i) = five_point(i) + 360.
          date_line = 1
       else if( -1 .eq. positive )then
          five_point(i) = five_point(i) - 360.
          date_line = 1
       endif
    end do

    Sum = 0.
    do i=1,5
       Lvalue = 1.
       do j=1,5
          if( j .ne. i )then
             Lvalue = Lvalue * (scanAngle(pos) - scanAngleIn(j+offset))/ &
                  (scanAngleIn(i+offset) - scanAngleIn(j+offset))
          endif
       end do
       Sum = Sum + Lvalue * five_point(i)
    end do

    ! If we've crossed the date line make sure that we check value to fit
    ! within +/- 180 - also takes care of case near the poles
    if( Sum .gt. 180. )then
       Sum = Sum - 360
    else if( Sum .lt. -180. )then
       Sum = Sum + 360
    endif
    doInterp = Sum

    return

  END FUNCTION doInterp

  SUBROUTINE getClavrFlags( dataType, clavr_mask, pos, record )

    INTEGER, INTENT(IN) :: dataType
    INTEGER(GbcsInt2), INTENT(OUT) :: clavr_mask(:)
    INTEGER, INTENT(IN) :: pos
    INTEGER(GbcsInt1), INTENT(IN) :: record(:)

    INTEGER :: nwords
    INTEGER(GbcsInt2) :: word
    INTEGER :: i, j, k
    INTEGER :: loc
    INTEGER :: temp

    IF( AVHRR_GAC .eq. dataType )THEN
       nwords = 52
    ELSE IF( AVHRR_LAC .eq. dataType )THEN
       nwords = 256
    ELSE
       call out_ERROR('dataType incorrect',&
            'getClavrFlags','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    ! Need to check if BIT order is correct
    ! ie. cloud is really cloud...
    ! yes it is correct: looking at images of the flag,
    ! areas go from 00 (clear) to 01 (probably clear) to 10 (probably cloudy)
    ! to 11 (cloudy)

    if( 52 .eq. nwords )then
       j = 0
       DO I=1,nwords-1
          word = header2short(pos+(i-1)*2,record)
          DO k=8,1,-1
             j = j+1
             loc = (k-1)*2
             clavr_mask(j) = IBITS(word,loc,2)
          END DO
       END DO
       word = header2short(pos+(nwords-1)*2,record)
       j = j+1
       clavr_mask(j) = IBITS(word,14,2)
    else
       j = 0
       DO I=1,nwords
          word = header2short(pos+(i-1)*2,record)
          DO k=8,1,-1
             j = j+1
             loc = (k-1)*2
             clavr_mask(j) = IBITS(word,loc,2)
          END DO
       END DO
    endif

  END SUBROUTINE getClavrFlags

  REAL FUNCTION convertTemperatureV3( intVal, Coefs, type )

    INTEGER(GbcsInt1), INTENT(IN) :: intVal
    REAL, INTENT(IN) :: Coefs(:)
    INTEGER, INTENT(IN) :: type

    INTEGER :: I
    INTEGER :: intValConv
    REAL :: volts
    REAL :: outTemp

    ! make sure that we're doing unsigned integers
    IF( 0 .gt. intVal )THEN
       intValConv = intVal + 256
    ELSE
       intValConv = intVal
    ENDIF
    volts = intValConv * 0.02
    outTemp = Coefs(1)
    do i=2,6
       outTemp = outTemp + Coefs(i) * (volts**(i-1))
    end do
    if( 1 .eq. type )then
       outTemp = outTemp + 273.15
    endif
    convertTemperatureV3 = outTemp

    return
  END FUNCTION convertTemperatureV3

  REAL FUNCTION convertCurrentV3( intVal, Coefs )

    INTEGER(GbcsInt1), INTENT(IN) :: intVal
    REAL, INTENT(IN) :: Coefs(:)

    INTEGER :: I
    INTEGER :: intValConv
    REAL :: volts
    REAL :: outCurr

    ! make sure that we're doing unsigned integers
    IF( 0 .gt. intVal )THEN
       intValConv = intVal + 256
    ELSE
       intValConv = intVal
    ENDIF
    volts = intValConv * 0.02
    outCurr = Coefs(1)
    do i=2,6
       outCurr = outCurr + Coefs(i) * (volts**(i-1))
    end do
    convertCurrentV3 = outCurr

    return
  END FUNCTION convertCurrentV3

  REAL FUNCTION convertTemperatureV2( intVal, Coefs, type )

    INTEGER(GbcsInt1), INTENT(IN) :: intVal
    REAL, INTENT(IN) :: Coefs(:)
    INTEGER, INTENT(IN) :: type

    INTEGER :: I
    INTEGER :: intValConv
    REAL :: volts
    REAL :: outTemp

    ! make sure that we're doing unsigned integers
    IF( 0 .gt. intVal )THEN
       intValConv = intVal + 256
    ELSE
       intValConv = intVal
    ENDIF
    volts = intValConv * 0.02
    outTemp = Coefs(1)
    do i=2,6
       outTemp = outTemp + Coefs(i) * (volts**(i-1))
    end do
    if( 1 .eq. type )then
       outTemp = outTemp + 273.15
    endif
    convertTemperatureV2 = outTemp

    return
  END FUNCTION convertTemperatureV2

  REAL FUNCTION convertCurrentV2( intVal, Coefs )

    INTEGER(GbcsInt1), INTENT(IN) :: intVal
    REAL, INTENT(IN) :: Coefs(:)

    INTEGER :: I
    INTEGER :: intValConv
    REAL :: volts
    REAL :: outCurr

    ! make sure that we're doing unsigned integers
    IF( 0 .gt. intVal )THEN
       intValConv = intVal + 256
    ELSE
       intValConv = intVal
    ENDIF
    volts = intValConv * 0.02

    outCurr = Coefs(1)
    do i=2,6
       outCurr = outCurr + Coefs(i) * (volts**(i-1))
    end do
    convertCurrentV2 = outCurr

    return
  END FUNCTION convertCurrentV2

  SUBROUTINE Parse_RecordV2(Coefs,DataType,scanLinePos,sizeOfRecords,&
       record,outData,bad_data,out_radiances,getTiePoints,getLandMask,&
       getLocation,printAllErrors,ignore_solar_contamination)

    TYPE(AVHRR_Instrument_Coefs), INTENT(IN) :: Coefs
    INTEGER, INTENT(IN) :: DataType
    INTEGER, INTENT(IN) :: scanLinePos
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER(GbcsInt1), INTENT(IN) :: record(sizeOfRecords)
    TYPE(AVHRR_Data), POINTER :: outData
    TYPE(AVHRR_Bad_Data), INTENT(IN) :: bad_data
    LOGICAL, INTENT(IN) :: out_radiances
    LOGICAL, OPTIONAL :: getTiePoints
    LOGICAL, OPTIONAL :: getLandMask
    LOGICAL, OPTIONAL :: getLocation
    LOGICAL, OPTIONAL :: printAllErrors
    LOGICAL, OPTIONAL :: ignore_solar_contamination

    INTEGER :: i,j,k
    LOGICAL :: badData
    LOGICAL :: bad_Top
    LOGICAL :: bad_Navigation
    LOGICAL :: bad_Time
    LOGICAL :: bad_Calibration
    INTEGER(GbcsInt2) :: shortVal
    INTEGER :: scanLineNumber
    INTEGER(GbcsInt4) :: qualityTop
    INTEGER(GbcsInt4) :: earthLocation
    INTEGER(GbcsInt2) :: qualityCalib(3)
    INTEGER(GbcsInt4) :: scanLineQuality
    INTEGER(GbcsInt1) :: qualityCal
    INTEGER(GbcsInt1) :: qualityTime
    INTEGER(GbcsInt1) :: qualityEarth
    INTEGER(GbcsInt4) :: syncErrorCount
    INTEGER :: scanLineBitField
    LOGICAL :: timeCorrected
    INTEGER(GbcsInt2) :: clockDrift_msecs
    INTEGER :: frameDayno, frameUTC_msecs

    LOGICAL :: filter3a
    INTEGER :: utc
    INTEGER :: intVal
    INTEGER :: stat

    REAL :: calib1(2)
    REAL :: calib1_2(2)
    REAL :: calib1_intercept
    REAL :: calib2(2)
    REAL :: calib2_2(2)
    REAL :: calib2_intercept
    REAL :: calib3A(2)
    REAL :: calib3A_2(2)
    REAL :: calib3A_intercept
    REAL :: calib3(3)
    REAL :: calib4(3)
    REAL :: calib5(3)

    REAL :: satelliteAltitude
    REAL :: solarAngleIn(51)
    REAL :: satelliteAngleIn(51)
    REAL :: relativeAngleIn(51)
    REAL :: latitudeIn(51)
    REAL :: longitudeIn(51)
    REAL :: d

    REAL :: solarAngle(NDATA_MAX)
    REAL :: sateliteAngle(NDATA_MAX)
    REAL :: relativeAngle(NDATA_MAX)
    REAL :: latitude(NDATA_MAX)
    REAL :: longitude(NDATA_MAX)

    INTEGER :: prtPosition = 0

    INTEGER :: word
    INTEGER :: earthCounts(5,NDATA_MAX)
    INTEGER :: counts
    REAL :: radiance

    REAL :: minVal
    INTEGER :: minPos
    REAL :: minPosArray(4)

    REAL :: diff

    INTEGER :: nMean3
    INTEGER :: nMean4
    INTEGER :: nMean5

    REAL :: prtCounts
    REAL :: prtCounts3
    REAL :: prtCounts4
    REAL :: prtCounts5

    LOGICAL :: solar_3b
    LOGICAL :: solar_4
    LOGICAL :: solar_5
    LOGICAL :: print_All_Errors
    LOGICAL :: transition3A3B
    LOGICAL :: ignore_solar_contam

    ! Stores 10 bit words from minor frame
    INTEGER(GbcsInt2) :: minor_frame(103)

    print_All_Errors = .FALSE.
    IF( PRESENT(printAllErrors) )THEN
       print_All_Errors = printAllErrors
    ENDIF

    ignore_solar_contam = .FALSE.
    IF( PRESENT(ignore_solar_contamination) )THEN
       ignore_solar_contam = ignore_solar_contamination
    ENDIF

    if( .not. AVHRR_setup )then
       if( AVHRR_GAC .eq. DataType )then
          ndata = NDATA_GAC
          ndata_counts = 682
          do i=1,4
             prtRotate(i) = prtRotate_GAC(i)
          end do
       else if( AVHRR_LAC .eq. DataType .or. AVHRR_FRAC .eq. DataType )then
          ndata = NDATA_LAC
          ndata_counts = 3414
          do i=1,4
             prtRotate(i) = prtRotate_LAC(i)
          end do
       else
          call out_ERROR('Data not GAC/LAC',&
               'Parse_RecordV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
       AVHRR_Setup = .TRUE.
       ! Allocate output arrays
       call Allocate_OutData( ndata, outData )
       outData%badTime = .TRUE.
       outData%badNavigation = .TRUE.
       outData%badCalibration = .TRUE.
    endif

    ! check to see if we have enought storage
    if( scanLinePos .gt. outData%arraySize )then
       call Reallocate_outData( outData, MEMORY_ALLOC_STEP )
    endif

    ! assume bad to start with
    outData%badTime(scanLinePos) = .TRUE.
    outData%badNavigation(scanLinePos) = .TRUE.
    outData%badCalibration(scanLinePos) = .TRUE.

    ! Get scan line number
    scanLineNumber = header2ushort( 0, record )
    outData%scanLineNumber(scanLinePos) = scanLineNumber

    ! Check quality
    qualityTop = header2int( 24, record )
    DO i=1,103
       ! (can use header2short if just testing bits)
       minor_frame(i) = IBITS(header2ushort(1056-2+i*2,Record),0,10)
    END DO
    syncErrorCount = header2ushort(38,Record)
    qualityCalib(1) = header2short( 32, record )
    qualityCalib(2) = header2short( 34, record )
    qualityCalib(3) = header2short( 36, record )
    scanLineQuality = header2int( 28, record )
    qualityTime = IBITS(scanLineQuality,16,8)
    qualityCal = IBITS(scanLineQuality,8,8)
    qualityEarth = IBITS(scanLineQuality,0,8)
    earthLocation = header2int( 312, record )
    scanLineBitField = header2ushort( 12, record )
    transition3A3B=IBITS(scanLineBitField,0,2).eq.2
    outData%transition3a3b(scanLinePos) = .TRUE.
    ! this field is for info only: after late 2002 there are no corrections
    ! because the spacecraft time is correct
    timeCorrected = BTEST(scanLineBitField,14)
    ! clockDrift_msecs = spacecraft time - UTC
    clockDrift_msecs = header2short( 6, record)
    ! N15 and N16 have timeCorrected=T but the scanline time is the same as the
    ! frame time.  And N17,18,19 should have no corrections required (corrected in
    ! spacecraft) but have timeCorrected=T and scanline time is the same as the
    ! frame time.  Checking what happens on the orbit where the corrections are
    ! switched on, it is clear by looking at the lon/lat that the frame time and
    ! scanline time have both been corrected in the L1B file (even though one would
    ! not expect the frame time to be correct until after Aug 2002). This is different
    ! from v1.

    badData = .FALSE.
    bad_Top = .FALSE.
    bad_Navigation = .FALSE.
    bad_Time = .FALSE.
    bad_Calibration = .FALSE.

    ! Get time
    outData%year(scanLinePos) = header2ushort(2,record)
    outData%dayno(scanLinePos) = header2ushort(4,record)
    CALL Parse_DayNo( outData%year(scanLinePos), outData%dayno(scanLinePos),&
         outData%month(scanLinePos), outData%day(scanLinePos) )
    utc = header2int(8,record)
    outData%UTC_msecs(scanLinePos) = utc
    outData%hours(scanLinePos) = (utc*1e-3/3600.)
    outData%time(scanLinePos) = get_date_time( outData%year(scanLinePos),&
         outData%dayno(scanLinePos),outData%hours(scanLinePos), time_yearstart )

!!$    PRINT'(2i6,2b2,4i3,4b9.8,1b11.10,2b9.8)',scanLinePos,scanLineNumber,&
!!$         IBITS( earthLocation, 18, 1 ),& ! manoeuvre in progress
!!$         IBITS( earthLocation, 17, 1 ),& ! location within tolerance
!!$         IBITS( earthLocation, 12, 4 ),& ! up-to-date earth location
!!$         IBITS( earthLocation, 8, 4 ),& ! in YGC or nominal with good attitude
!!$         IBITS( earthLocation, 4, 4 ),& ! in nominal; or undocumented value
!!$         IBITS( earthLocation, 0, 4 ),& ! in nominal, no tests; or undocumented value
!!$         IBITS( qualityTop, 24, 8 ),&
!!$         IBITS( qualityTop, 16, 8 ),&
!!$         IBITS( qualityTop, 8, 8 ),&
!!$         IBITS( qualityTop, 0, 8 ),&
!!$         IBITS( minor_frame(7), 0, 10 ),&
!!$         qualityTime, qualityEarth

    bad_Top = CheckTopV2_V3(qualityTop, minor_frame, syncErrorCount)
    bad_Navigation = CheckNavigationV2_V3(qualityTop, qualityTime, &
         qualityEarth, earthLocation, outData%AVHRR_No)
    ! if using clavr-x lon/lat then need to check those
    IF (.NOT.useL1BClock) THEN
       IF (ANY(clavrLonAnchor(:,scanLinePos).EQ.NAN_R) &
            .OR.ANY(clavrLatAnchor(:,scanLinePos).EQ.NAN_R))THEN
          bad_Navigation = .TRUE.
       ELSE
          ! check the nadir element
          i = 26
          latitudeIn(i) = header2int((632+i*8),Record)/1e4
          longitudeIn(i) = header2int((636+i*8),Record)/1e4
          CALL greatCircleDistance(longitudeIn(i),latitudeIn(i),clavrLonAnchor(i,scanLinePos),clavrLatAnchor(i,scanLinePos),d)
          IF(d.GT.maxTimingErrorDistance) THEN
             bad_Navigation = .TRUE.
          END IF
       END IF
    END IF
    bad_Time = CheckTimeV2_V3(qualityTop, qualityTime)
    bad_Calibration = CheckCalibrationV2_V3(qualityTop, qualityCal, solar_3b, solar_4, solar_5)

    ! too bad to use at all
    IF( bad_Top )THEN
       IF( 0 .EQ. printOnce(1) .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (top/frame errors) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       ENDIF
       printOnce(1) = 1
       IF( bad_data%allow_bad_top )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_outdata_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    ! bad time means bad navigation
    if( bad_Time )then
       if( 0 .eq. printOnce(3) .or. print_All_Errors )then
          write(*,'(''Quality Problems (time) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       endif
       printOnce(3) = 1
       if( bad_data%allow_bad_time )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_outdata_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badTime(scanLinePos) = .FALSE.
    endif

    IF( bad_Time )THEN
       CALL set_time_bad( outData, scanLinePos )
    ENDIF

    ! CheckNavigationV2_V3 also checks the time
    if( bad_Navigation )then
       if( 0 .eq. printOnce(2) .or. print_All_Errors )then
          write(*,&
               '(''Quality Problems (navigation) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
          if( bad_data%use_bad_nav )then
             write(*,'(''  Using potentially bad navigation '')')
          endif
       endif
       printOnce(2) = 1
       if( bad_data%allow_bad_nav .or. bad_data%use_bad_nav )then
          badData = .TRUE.
          if( bad_data%use_bad_nav )then
             ! If use navigation then it is not considered 'bad'
             bad_navigation = .FALSE.
             outdata%badNavigation(scanLinePos) = .FALSE.
          endif
       else
          prtNumber = -1
          ! everything except time
          call set_nav_bad( outData, scanLinePos )
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badNavigation(scanLinePos) = .FALSE.
    endif

    if( .not. bad_Navigation )then
       ! Satellite angle in degrres
       do i=1,51
          ! Note that the start numbers are offset bu i*6 or i*8 from the
          ! record numbers
          solarAngleIn(i) = header2short((322+i*6),Record)/1e2
          satelliteAngleIn(i) = header2short((324+i*6),Record)/1e2
          relativeAngleIn(i) = header2short((326+i*6),Record)/1e2
          latitudeIn(i) = header2int((632+i*8),Record)/1e4
          longitudeIn(i) = header2int((636+i*8),Record)/1e4
       end do

       IF (.NOT.useL1BClock) THEN
          ! use the clavr-x nav lon/lat that has been read in
          longitudeIn = clavrLonAnchor(:,scanLinePos)
          latitudeIn = clavrLatAnchor(:,scanLinePos)
       ENDIF

       CALL InterpolateData(solarAngleIn,solarAngle,dataType)
       CALL InterpolateAngles(relativeAngleIn,satelliteAngleIn,relativeAngle,sateliteAngle,dataType)

       CALL InterpolateLonLat(longitudeIn,latitudeIn,longitude(1:ndata),latitude(1:ndata),dataType)
       outData%Lon(:,scanLinePos) = longitude(1:ndata)
       outData%Lat(:,scanLinePos) = latitude(1:ndata)
       outData%satZA(:,scanLinePos) = sateliteAngle(1:ndata)
       outData%solZA(:,scanLinePos) = solarAngle(1:ndata)
       outData%relAz(:,scanLinePos) = relativeAngle(1:ndata)
    else
       call set_nav_bad( outData, scanLinePos )
    endif

    IF (PRESENT(getTiePoints)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLandMask)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLocation)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF( bad_Calibration )THEN
       IF( 0 .EQ. printonce(4)  .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (calibration) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       printonce(4) = 1
       IF( bad_data%allow_bad_calib )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    ! Are we dealing with 3B
    shortVal = header2short(1068,record)
    if( BTEST(shortVal,0) )then
       filter3a = .TRUE.
    else
       filter3a = .FALSE.
    endif

    IF( .NOT.filter3a &
         .AND. (0 .NE. qualityCalib(1) .OR. transition3A3B) )THEN
       if( 0 .eq. printonce(5)  .or. print_All_Errors )then
          WRITE(*,'(''Quality Problems (Calib 1) or 3A/B transition, scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(5) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( 0 .ne. qualityCalib(2) )then
       if( 0 .eq. printonce(6)  .or. print_All_Errors )then
          WRITE(*,'(''Quality Problems (Calib 2) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(6) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( 0 .ne. qualityCalib(3) )then
       if( 0 .eq. printonce(7)  .or. print_All_Errors )then
          WRITE(*,'(''Quality Problems (Calib 2) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(7) = 1
       if( bad_data%allow_bad_calib )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    endif

    if( BTEST(qualityCal,5) )then
       if( 0 .eq. printonce(8)  .or. print_All_Errors )then
          write(*,'(''PRTs not calibrated on scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(8) = 1
       prtNumber = -1
       call set_calib_bad( outData, scanLinePos )
       return
    endif

    if( BTEST(qualityCal,4) )then
       if( 0 .eq. printonce(9)  .or. print_All_Errors )then
          write(*,'(''PRTs marginally calibrated on scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       endif
       printonce(9) = 1
       prtNumber = -1
       call set_calib_bad( outData, scanLinePos )
       return
    endif

    if( .not. badData )then
       badDataLoop: do i=1,10
          if( 1 .eq. printOnce(i) )then
             printOnce = (/0,0,0,0,0,0,0,0,0,0/)
             write(*,'(''End of quality problems on scan line, position'',2i7)')&
                  scanLineNumber, scanLinePos
             exit badDataLoop
          endif
       end do badDataLoop
    endif

    ! Solar flags
    outData%orig_solar_contamination_3b(scanLinePos) = solar_3b
    outData%orig_solar_contamination_4(scanLinePos) = solar_4
    outData%orig_solar_contamination_5(scanLinePos) = solar_5

    ! Get calibration values
    calib1(2) = header2int(48,Record)/1e7
    calib1(1) = header2int(52,Record)/1e6
    calib1_2(2) = header2int(56,Record)/1e7
    calib1_2(1) = header2int(60,Record)/1e6
    calib1_intercept = header2int(64,Record)
    ! Channel 2
    calib2(2) = header2int(108,Record)/1e7
    calib2(1) = header2int(112,Record)/1e6
    calib2_2(2) = header2int(116,Record)/1e7
    calib2_2(1) = header2int(120,Record)/1e6
    calib2_intercept = header2int(124,Record)
    ! Infrared
    if( .not. filter3a )then
       calib3(1) = header2int(228,Record)/1e6
       calib3(2) = header2int(232,Record)/1e6
       calib3(3) = header2int(236,Record)/1e6
       ! Channel 3A
       calib3A(2) = 0.
       calib3A(1) = 0.
       calib3A_2(2) = 0.
       calib3A_2(1) = 0.
       calib3A_intercept = 0
    else
       calib3(1) = 0.
       calib3(2) = 0.
       calib3(3) = 0.
       ! Channel 3A
       calib3A(2) = header2int(168,Record)/1e7
       calib3A(1) = header2int(172,Record)/1e6
       calib3A_2(2) = header2int(176,Record)/1e7
       calib3A_2(1) = header2int(180,Record)/1e6
       calib3A_intercept = header2int(184,Record)
    endif
    calib4(1) = header2int(252,Record)/1e6
    calib4(2) = header2int(256,Record)/1e6
    calib4(3) = header2int(260,Record)/1e6
    calib5(1) = header2int(276,Record)/1e6
    calib5(2) = header2int(280,Record)/1e6
    calib5(3) = header2int(284,Record)/1e6

    outData%calib1(:,scanLinePos) = calib1
    outData%calib2(:,scanLinePos) = calib2
    outData%calib1_2(:,scanLinePos) = calib1_2
    outData%calib2_2(:,scanLinePos) = calib2_2
    outData%calib1_intercept(scanLinePos) = calib1_intercept
    outData%calib2_intercept(scanLinePos) = calib2_intercept
    outData%calib3A(:,scanLinePos) = calib3A
    outData%calib3A_2(:,scanLinePos) = calib3A_2
    outData%calib3A_intercept(scanLinePos) = calib3A_intercept
    outData%calib3(:,scanLinePos) = calib3
    outData%calib4(:,scanLinePos) = calib4
    outData%calib5(:,scanLinePos) = calib5

!  roll = header2short(320,Record)/1e3
!  pitch = header2short(322,Record)/1e3
!  yaw = header2short(324,Record)/1e3

    satelliteAltitude = header2ushort(326,Record)/10.
    outData%satelliteAlt_There = .TRUE.
    outData%satelliteAltitude(scanLinePos) = satelliteAltitude

    if( AVHRR_GAC .eq. dataType )then
       outData%radiator(scanLinePos) = convertTemperatureV2(record(4024),&
            Coefs%temperatureCoefs1,0)
       outData%electronics(scanLinePos) = convertTemperatureV2(record(4032),&
            Coefs%temperatureCoefs2,1)
       outData%cooler(scanLinePos) = convertTemperatureV2(record(4033),&
            Coefs%temperatureCoefs3,1)
       outData%baseplate(scanLinePos) = convertTemperatureV2(record(4034),&
            Coefs%temperatureCoefs4,1)
       outData%motor(scanLinePos) = convertTemperatureV2(record(4035),&
            Coefs%temperatureCoefs5,1)
       outData%a_d_conv(scanLinePos) = convertTemperatureV2(record(4036),&
            Coefs%temperatureCoefs6,1)
       outData%patch(scanLinePos) = convertTemperatureV2(record(4021),&
            Coefs%patchCoef,0)
       outData%patchExtended(scanLinePos) = convertTemperatureV2(record(4022),&
            Coefs%PatchCoefExt,0)

       ! Get motor current
       outData%motorCurrent(scanLinePos) = convertCurrentV2(record(4030),&
            Coefs%motorCurrentCoefs)

       ! Get CLAVR flags
       intVal = header2int(4048,record)
       if( btest(intVal, 0) )then
          outData%Clavr_There = .TRUE.
          CALL getClavrFlags( AVHRR_GAC, outData%clavr_mask(:,scanLinePos), &
               4056, record )
       else
          outData%clavr_mask(:,scanLinePos) = 3
       endif
    else if( AVHRR_LAC .eq. dataType )then
       outData%radiator(scanLinePos) = convertTemperatureV2(record(14952),&
            Coefs%temperatureCoefs1,0)
       outData%electronics(scanLinePos) = convertTemperatureV2(record(14960),&
            Coefs%temperatureCoefs2,1)
       outData%cooler(scanLinePos) = convertTemperatureV2(record(14961),&
            Coefs%temperatureCoefs3,1)
       outData%baseplate(scanLinePos) = convertTemperatureV2(record(14962),&
            Coefs%temperatureCoefs4,1)
       outData%motor(scanLinePos) = convertTemperatureV2(record(14963),&
            Coefs%temperatureCoefs5,1)
       outData%a_d_conv(scanLinePos) = convertTemperatureV2(record(14964),&
            Coefs%temperatureCoefs6,1)
       outData%patch(scanLinePos) = convertTemperatureV2(record(14949),&
            Coefs%patchCoef,0)
       outData%patchExtended(scanLinePos) = convertTemperatureV2(record(14950),&
            Coefs%PatchCoefExt,0)
       ! Get motor current
       outData%motorCurrent(scanLinePos) = convertCurrentV2(record(14958),&
            Coefs%motorCurrentCoefs)

       ! Get CLAVR flags
       intVal = header2int(14976,record)
       if( btest(intVal, 0) )then
          outData%Clavr_There = .TRUE.
          CALL getClavrFlags( AVHRR_LAC, outData%clavr_mask(:,scanLinePos), &
               14984, record )
       else
          outData%clavr_mask(:,scanLinePos) = 3
       endif
    endif

    prtCounts3 = header2ushort(1090,record)
    prtCounts4 = header2ushort(1092,record)
    prtCounts5 = header2ushort(1094,record)
!---Modif du 12 nov---
!on veut imprimer les prt counts pour implementer la methode de Trishchenko
    !outData%prtcounts(1,scanLinePos)=prtCounts3
    !outData%prtcounts(2,scanLinePos)=prtCounts4
    !outData%prtcounts(3,scanLinePos)=prtCounts5
    
!-------------------

    ! Need to get which PRT we are at
    if( -1 .eq. prtNumber )then
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
       endif
       storedPrtTemp = 0.
    else
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtCountsPrevStore(i) = prtCountsStore(i)
             else
                prtCountsPrevStore(i) = 0
             endif
          end do
       else
          prtCounts = (prtCounts3 + prtCounts4 + prtCounts5)/3.
          if( 1 .gt. prtNumber .or. 4 .lt. prtNumber )then
             ! Rejig array

             ! To be honest I can't remember why I wrote the following code
             ! should never really be used...

             if( 0 .lt. prtCountsPrevStore(1) .and. &
                  0 .lt. prtCountsPrevStore(2) .and. &
                  0 .lt. prtCountsPrevStore(3) .and. &
                  0 .lt.  prtCountsPrevStore(4) )then
                minPos = 0
                minVal = 1e32
                do j=1,4
                   diff = abs(prtCounts - prtCountsPrevStore(j))
                   if( minVal > diff )then
                      minVal = diff
                      minPos = j
                   endif
                end do
                minPosArray(4) = minPos
                do i=1,4
                   minPos = 0
                   minVal = 1e32
                   do j=1,4
                      diff = abs(prtCountsStore(i) - prtCountsPrevStore(j))
                      if( minVal > diff )then
                         minVal = diff
                         minPos = j
                      endif
                   end do
                   minPosArray(i) = minPos
                end do

                ! Now have array of indexes to matching entries - remap these
                ! ones
                ! Find repeat
                minPos = -1
                do i=2,4
                   if( minPosArray(i-1) .eq. minPosArray(i) )then
                      minPos = i
                   endif
                end do

                if( -1 .eq. minPos )then
                   ! Assume line is bad (for some reason) and set bad
                   prtNumber = -1
                   call set_calib_bad( outData, scanLinePos )
                   return
                else

                   ! Shift down
                   do i=minPos,3
                      minPosArray(i) = minPosArray(i+1)
                   end do

                   ! Set extra value to last value
                   minPosArray(4) = prtCounts
                endif

                write(*,*)'Resetting prtCounts array'

             else
                write(*,*)'No previous prt record - just taking first 4'
                prtNumber=1
             endif

          else
             prtCountsStore(prtNumber) = prtCounts
             prtCountsStoreAll(1,prtNumber) = prtCounts3
             prtCountsStoreAll(2,prtNumber) = prtCounts4
             prtCountsStoreAll(3,prtNumber) = prtCounts5
             prtNumber = prtNumber + 1
          endif

          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtPosition = prtRotate(i)
                if( 1 .gt. prtPosition .or. 4 .lt. prtPosition )then
                   call out_ERROR('prtPosition out of range',&
                        'Parse_RecordV2','NOAA_LoadAVHRRLevel1B.f90',&
                        .FALSE.,' ')
                endif
                storedPrtCounts(prtPosition) = prtCountsStore(i)
                storedPrtCountsAll(:,prtPosition) = prtCountsStoreAll(:,i)
                storedPrtTemp(prtPosition) = &
                     Coefs%prtTempCoefs1(prtPosition) + &
                     Coefs%prtTempCoefs2(prtPosition)*prtCountsStore(i) + &
                     Coefs%prtTempCoefs3(prtPosition)*(prtCountsStore(i)**2) + &
                     Coefs%prtTempCoefs4(prtPosition)*(prtCountsStore(i)**3) + &
                     Coefs%prtTempCoefs5(prtPosition)*(prtCountsStore(i)**4)
             endif
          end do
       end if
    end if

    if( storedPrtTemp(1) .gt. 200 .and. &
         storedPrtTemp(2) .gt. 200 .and. &
         storedPrtTemp(3) .gt. 200 .and. &
         storedPrtTemp(4) .gt. 200 )then
       outData%prt1CountsAll(:,scanLinePos) = storedPrtCountsAll(:,1)
       outData%prt2CountsAll(:,scanLinePos) = storedPrtCountsAll(:,2)
       outData%prt3CountsAll(:,scanLinePos) = storedPrtCountsAll(:,3)
       outData%prt4CountsAll(:,scanLinePos) = storedPrtCountsAll(:,4)
       outData%prt1Counts(scanLinePos) = storedPrtCounts(1)
       outData%prt2Counts(scanLinePos) = storedPrtCounts(2)
       outData%prt3Counts(scanLinePos) = storedPrtCounts(3)
       outData%prt4Counts(scanLinePos) = storedPrtCounts(4)
       outData%prt1(scanLinePos) = storedPrtTemp(1)
       outData%prt2(scanLinePos) = storedPrtTemp(2)
       outData%prt3(scanLinePos) = storedPrtTemp(3)
       outData%prt4(scanLinePos) = storedPrtTemp(4)
    else
       outData%prt1CountsAll(:,scanLinePos) = NAN_R
       outData%prt2CountsAll(:,scanLinePos) = NAN_R
       outData%prt3CountsAll(:,scanLinePos) = NAN_R
       outData%prt4CountsAll(:,scanLinePos) = NAN_R
       outData%prt1Counts(scanLinePos) = NAN_R
       outData%prt2Counts(scanLinePos) = NAN_R
       outData%prt3Counts(scanLinePos) = NAN_R
       outData%prt4Counts(scanLinePos) = NAN_R
       outData%prt1(scanLinePos) = NAN_R
       outData%prt2(scanLinePos) = NAN_R
       outData%prt3(scanLinePos) = NAN_R
       outData%prt4(scanLinePos) = NAN_R
    endif

    ! Get BB counts
    outData%bbodyFilter3(:,scanLinePos) = 0
    outData%bbodyFilter4(:,scanLinePos) = 0
    outData%bbodyFilter5(:,scanLinePos) = 0
    outData%bb3(scanLinePos) = 0.
    outData%bb4(scanLinePos) = 0.
    outData%bb5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       if( .not.filter3a )then
          outData%bbodyFilter3(i,scanLinePos) = &
               header2ushort((1100+(i-1)*6),record)
       endif
       outData%bbodyFilter4(i,scanLinePos) = &
            header2ushort((1102+(i-1)*6),record)
       outData%bbodyFilter5(i,scanLinePos) = &
            header2ushort((1104+(i-1)*6),record)
       ! in fact outData%bbodyFilter? are always +ve
       if( .not. filter3a )then
          if( 0 .lt. outData%bbodyFilter3(i,scanLinePos) )then
             outData%bb3(scanLinePos) = outData%bb3(scanLinePos) + &
                  outData%bbodyFilter3(i,scanLinePos)
             nMean3 = nMean3 + 1
          endif
       endif
       if( 0 .lt. outData%bbodyFilter4(i,scanLinePos) )then
          outData%bb4(scanLinePos) = outData%bb4(scanLinePos) + &
               outData%bbodyFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 0 .lt. outData%bbodyFilter5(i,scanLinePos) )then
          outData%bb5(scanLinePos) = outData%bb5(scanLinePos) + &
               outData%bbodyFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( .not. filter3a )then
       if( NBBODY_SAMPLES .ne. nMean3 )then
          write(*,*)'Dropped Black-body (3.7) counts'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       endif
       if( 0 .eq. nMean3 )then
          outData%bb3(scanLinePos) = NAN_R
          write(*,*)'Black-body counts (3.7) are zero'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       else
          outData%bb3(scanLinePos) = outData%bb3(scanLinePos) / nMean3
       end if
    else
       outData%bb3(scanLinePos) = NAN_R
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Black-body (11) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean4 )then
       outData%bb4(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb4(scanLinePos) = outData%bb4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Black-body (12) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean5 )then
       outData%bb5(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb5(scanLinePos) = outData%bb5(scanLinePos) / nMean5
    end if

    ! Get space counts
    outData%spaceFilter1(:,scanLinePos) = 0
    outData%spaceFilter2(:,scanLinePos) = 0
    outData%spaceFilter3a(:,scanLinePos) = 0
    outData%spaceFilter3(:,scanLinePos) = 0
    outData%spaceFilter4(:,scanLinePos) = 0
    outData%spaceFilter5(:,scanLinePos) = 0
    outData%sp3(scanLinePos) = 0.
    outData%sp4(scanLinePos) = 0.
    outData%sp5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       outData%spaceFilter1(i,scanLinePos) = &
            header2ushort((1160+(i-1)*10),record)
       outData%spaceFilter2(i,scanLinePos) = &
            header2ushort((1162+(i-1)*10),record)
       if( .not.filter3a )then
          outData%spaceFilter3(i,scanLinePos) = &
               header2ushort((1164+(i-1)*10),record)
          outData%spaceFilter3a(i,scanLinePos) = NAN_I
       else
          outData%spaceFilter3a(i,scanLinePos) = &
               header2ushort((1164+(i-1)*10),record)
          outData%spaceFilter3(i,scanLinePos) = NAN_I
       endif
       outData%spaceFilter4(i,scanLinePos) = &
            header2ushort((1166+(i-1)*10),record)
       outData%spaceFilter5(i,scanLinePos) = &
            header2ushort((1168+(i-1)*10),record)
       ! in fact outData%spaceFilter? are always +ve
       if( .not. filter3a )then
          if( 900 .lt. outData%spaceFilter3(i,scanLinePos) )then
             outData%sp3(scanLinePos) = outData%sp3(scanLinePos) + &
                  outData%spaceFilter3(i,scanLinePos)
             nMean3 = nMean3 + 1
          endif
       endif
       if( 900 .lt. outData%spaceFilter4(i,scanLinePos) )then
          outData%sp4(scanLinePos) = outData%sp4(scanLinePos) + &
               outData%spaceFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 900 .lt. outData%spaceFilter5(i,scanLinePos) )then
          outData%sp5(scanLinePos) = outData%sp5(scanLinePos) + &
               outData%spaceFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( .not. filter3a )then
       if( NBBODY_SAMPLES .ne. nMean3 )then
          write(*,*)'Dropped Space (3.7) counts'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       endif
       if( 0 .eq. nMean3 )then
          outData%sp3(scanLinePos) = NAN_R
          write(*,*)'Space counts (3.7) are zero'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       else
          outData%sp3(scanLinePos) = outData%sp3(scanLinePos) / nMean3
       end if
    else
       outData%sp3(scanLinePos) = NAN_R
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Space (11) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean4 )then
       outData%sp4(scanLinePos) = NAN_R
       write(*,*)'Space counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp4(scanLinePos) = outData%sp4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Space (12) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean5 )then
       outData%sp5(scanLinePos) = NAN_R
       write(*,*)'Space counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp5(scanLinePos) = outData%sp5(scanLinePos) / nMean5
    end if

    ! L1B file calibration is OK but if ignoring solar contam data, ignore
    IF( ignore_solar_contam .and. (solar_3b .or. solar_4 .or. solar_5) )THEN
       outData%badCalibration(scanLinePos) = .TRUE.
    ELSE
       outData%badCalibration(scanLinePos) = .FALSE.
    ENDIF

    if( .not. bad_Navigation .and. &
         .not.outData%badCalibration(scanLinePos) )then

       ! Get earth observations
       j = 0
       k = 0
       do i=0,ndata_counts-1
          ! Data is interlaced
          word = header2int((1264+i*4),record)
          ! Get bits 20-29
          earthCounts(k+1,j+1) = IBITS(word,20,10)
          k = k+1
          if( 5 .eq. k)then
             j = j + 1
             k = 0
          end if
          if( 682 .eq. ndata_counts )then
             ! Get bits 10-19
             earthCounts(k+1,j+1) = IBITS(word,10,10)
             k = k+1
             if( 5 .eq. k)then
                j = j + 1
                k = 0
             end if
             ! Get bits 0-9 if i != 681
             if( ndata_counts-1 .ne. i )then
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          else
             ! Get bits 10-19 && 0-9 if i != 3413
             if( ndata_counts-1 .ne. i )then
                ! Get bits 10-19
                earthCounts(k+1,j+1) = IBITS(word,10,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
                ! Get bits 0-9
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k+1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          end if
       end do

       if( 0 .ne. k .and. j .ne. ndata )then
          call out_ERROR('Problems parsing earth observations',&
               'Parse_RecordV2','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif

       do i=1,ndata
          ! Visible
          ! Channel 1
          counts = earthCounts(1,i)
          outData%counts1(i,scanLinePos) = counts
          if( counts .ge. calib1_intercept )then
             outData%array1(i,scanLinePos) = calib1_2(1) + &
                  calib1_2(2) * counts
          else
             outData%array1(i,scanLinePos) = calib1(1) + &
                  calib1(2) * counts
          endif
          if( 0 .gt. outData%array1(i,scanLinePos) )then
             outData%array1(i,scanLinePos) = 0.
          endif
          outData%array1(i,scanLinePos) &
               = outData%array1(i,scanLinePos) / 100.0 ! percent to fraction

          ! Channel 2
          counts = earthCounts(2,i)
          outData%counts2(i,scanLinePos) = counts
          if( counts .ge. calib2_intercept )then
             outData%array2(i,scanLinePos) = calib2_2(1) + &
                  calib2_2(2) * counts
          else
             outData%array2(i,scanLinePos) = calib2(1) + &
                  calib2(2) * counts
          endif
          if( 0 .gt. outData%array2(i,scanLinePos) )then
             outData%array2(i,scanLinePos) = 0.
          endif
          outData%array2(i,scanLinePos) &
               = outData%array2(i,scanLinePos) / 100.0 ! percent to fraction

          ! Channel 3 - either 3A/3B
          if( filter3a )then

             IF( transition3A3B )THEN
                outData%array3A(i,scanLinePos) = NAN_R
             ELSE
                counts = earthCounts(3,i)
                IF( counts .GE. calib3A_intercept )THEN
                   outData%array3A(i,scanLinePos) = calib3A_2(1) + &
                        calib3A_2(2) * counts
                ELSE
                   outData%array3A(i,scanLinePos) = calib3A(1) + &
                        calib3A(2) * counts
                ENDIF
                IF( 0 .GT. outData%array3A(i,scanLinePos) )THEN
                   outData%array3A(i,scanLinePos) = 0.
                ENDIF
                outData%array3A(i,scanLinePos) &
                     = outData%array3A(i,scanLinePos) / 100.0 ! percent to fraction
             END IF

             ! it would be better to have a outData%filter3a
             ! flag but counts3 is *currently*) only used
             ! for the new calibration of 3B
             outData%counts3(i,scanLinePos) = NAN_R
             outData%array3B(i,scanLinePos) = NAN_R
             outData%transition3A3B(scanLinePos) = .FALSE.

          else

             IF( transition3A3B )THEN
                outData%array3B(i,scanLinePos) = NAN_R
                outData%counts3(i,scanLinePos) = NAN_R
                outData%transition3A3B(scanLinePos) = .FALSE.
             ELSE
                outData%transition3A3B(scanLinePos) = .TRUE.
                counts = earthCounts(3,i)
                outData%counts3(i,scanLinePos) = counts
                radiance = calib3(1) + calib3(2)*counts + &
                     calib3(3)*counts*counts
                IF( out_radiances )THEN
                   outData%array3B(i,scanLinePos) = radiance
                ELSE
                   outData%array3B(i,scanLinePos) = convertBT( radiance, &
                        Coefs%nuC(1), Coefs%Aval(1), Coefs%Bval(1) )
#ifdef USE_IEEE
                   IF( .NOT.ieee_is_finite(outData%array3B(i,scanLinePos)) )THEN
                      outData%array3B(i,scanLinePos) = NAN_R
                   END IF
#elif defined(USE_G95)
                   IF( isNAN(outData%array3B(i,scanLinePos)) )THEN
                      outData%array3B(i,scanLinePos) = NAN_R
                   END IF
#endif
                ENDIF
             END IF

             outData%array3A(i,scanLinePos) = NAN_R

          endif

          counts = earthCounts(4,i)
          outData%counts4(i,scanLinePos) = counts
          radiance = calib4(1) + calib4(2)*counts + &
               calib4(3)*counts*counts
          if( out_radiances )then
             outData%array4(i,scanLinePos) = radiance
          else
             outData%array4(i,scanLinePos) = convertBT( radiance, &
                  Coefs%nuC(2), Coefs%Aval(2), Coefs%Bval(2) )
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#endif
          endif

          counts = earthCounts(5,i)
          outData%counts5(i,scanLinePos) = counts
          radiance = calib5(1) + calib5(2)*counts + &
               calib5(3)*counts*counts
          if( out_radiances )then
             outData%array5(i,scanLinePos) = radiance
          else
             outData%array5(i,scanLinePos) = convertBT( radiance, &
                  Coefs%nuC(3), Coefs%Aval(3), Coefs%Bval(3) )
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array5(i,scanLinePos)) )THEN
                outData%array5(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array5(i,scanLinePos)) )THEN
                outData%array5(i,scanLinePos) = NAN_R
             END IF
#endif
          endif

       end do

    else

       ! set Earth observations to NAN because we don't know where we are
       outData%counts3(:,scanLinePos) = NAN_R
       outData%counts4(:,scanLinePos) = NAN_R
       outData%counts5(:,scanLinePos) = NAN_R
       outData%array1(:,scanLinePos) = NAN_R
       outData%array2(:,scanLinePos) = NAN_R
       outData%array3A(:,scanLinePos) = NAN_R
       outData%array3B(:,scanLinePos) = NAN_R
       outData%array4(:,scanLinePos) = NAN_R
       outData%array5(:,scanLinePos) = NAN_R

    endif
    outData%dataFilled = .TRUE.

  END SUBROUTINE Parse_RecordV2

  SUBROUTINE Setup_Walton(outData,walton_str)

    TYPE(AVHRR_Data), INTENT(IN) :: outData
    TYPE(Walton_Struct), INTENT(OUT) :: walton_str

    INTEGER :: I,J,K
    INTEGER :: STAT
    INTEGER :: miny,maxy
    INTEGER :: ntotal3B
    INTEGER :: ntotal4
    INTEGER :: ntotal5
    INTEGER :: ntotalprt
    REAL :: prt1
    REAL :: prt2
    REAL :: prt3
    REAL :: prt4
    REAL :: total3B
    REAL :: total4
    REAL :: total5

    ALLOCATE(walton_str%sm_bb(3,outData%arraySize),&
         walton_str%sm_sp(3,outData%arraySize),&
         walton_str%sm_prt(4,outData%arraySize),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate arrays',&
            'Setup_Walton','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF
    walton_str%sm_bb = NAN_R
    walton_str%sm_sp = NAN_R
    walton_str%sm_prt = NAN_R

    DO I=NPIXEL_PRT_SMOOTH+1,outData%arraySize-NPIXEL_PRT_SMOOTH
       miny = I-NPIXEL_PRT_SMOOTH
       maxy = I+NPIXEL_PRT_SMOOTH
       ntotalprt = 0
       prt1=0.
       prt2=0.
       prt3=0.
       prt4=0.
       DO J=miny,maxy
          IF( 0 .lt. outData%prt1(J) .and. &
                0. .lt. outData%prt2(J) .and. &
                0. .lt. outData%prt3(J) .and. &
                0. .lt. outData%prt4(J) )THEN
             ntotalprt = ntotalprt + 1
             prt1 = prt1 + outData%prt1(J)
             prt2 = prt2 + outData%prt2(J)
             prt3 = prt3 + outData%prt3(J)
             prt4 = prt4 + outData%prt4(J)
          ENDIF
       END DO
       IF( 0 .lt. ntotalprt )THEN
          prt1 = prt1/ntotalprt
          prt2 = prt2/ntotalprt
          prt3 = prt3/ntotalprt
          prt4 = prt4/ntotalprt
          walton_str%sm_prt(1,I) = prt1
          walton_str%sm_prt(2,I) = prt2
          walton_str%sm_prt(3,I) = prt3
          walton_str%sm_prt(4,I) = prt4
       ENDIF

       ntotal3B = 0
       ntotal4 = 0
       ntotal5 = 0
       total3B = 0.
       total4 = 0.
       total5 = 0.
       DO J=miny,maxy
          DO K=1,10
             IF( outData%spaceFilter3(K,I) .gt. 0 )THEN
                total3B = total3B + outData%spaceFilter3(K,I)
                ntotal3B = ntotal3B + 1
             ENDIF
             IF( outData%spaceFilter4(K,I) .gt. 0 )THEN
                total4 = total4 + outData%spaceFilter4(K,I)
                ntotal4 = ntotal4 + 1
             ENDIF
             IF( outData%spaceFilter5(K,I) .gt. 0 )THEN
                total5 = total5 + outData%spaceFilter5(K,I)
                ntotal5 = ntotal5 + 1
             ENDIF
          END DO
       END DO
       IF( 0 .lt. ntotal3B )THEN
          walton_str%sm_sp(1,I) = total3B/ntotal3B
       ENDIF
       IF( 0 .lt. ntotal4 )THEN
          walton_str%sm_sp(2,I) = total4/ntotal4
       ENDIF
       IF( 0 .lt. ntotal5 )THEN
          walton_str%sm_sp(3,I) = total5/ntotal5
       ENDIF

       ntotal3B = 0
       ntotal4 = 0
       ntotal5 = 0
       total3B = 0.
       total4 = 0.
       total5 = 0.
       DO J=miny,maxy
          DO K=1,10
             IF( outData%bbodyFilter3(K,I) .gt. 0 )THEN
                total3B = total3B + outData%bbodyFilter3(K,I)
                ntotal3B = ntotal3B + 1
             ENDIF
             IF( outData%bbodyFilter4(K,I) .gt. 0 )THEN
                total4 = total4 + outData%bbodyFilter4(K,I)
                ntotal4 = ntotal4 + 1
             ENDIF
             IF( outData%bbodyFilter5(K,I) .gt. 0 )THEN
                total5 = total5 + outData%bbodyFilter5(K,I)
                ntotal5 = ntotal5 + 1
             ENDIF
          END DO
       END DO
       IF( 0 .lt. ntotal3B )THEN
          walton_str%sm_bb(1,I) = total3B/ntotal3B
       ENDIF
       IF( 0 .lt. ntotal4 )THEN
          walton_str%sm_bb(2,I) = total4/ntotal4
       ENDIF
       IF( 0 .lt. ntotal5 )THEN
          walton_str%sm_bb(3,I) = total5/ntotal5
       ENDIF

    END DO

    !
    ! These are the coefficients in the files provided from CSPP
    ! as it contains all AVHRRs
    !
    SELECT CASE(outData%AVHRR_No)
    CASE(1)
       walton_str%nuc(1) = 2655.7409
       walton_str%aval(1) = -1.6485446
       walton_str%bval(1) = 1.0020894
       walton_str%Nsp(1) = -0.0039
       walton_str%c(:,1) = (/0.00195,-0.015,0.011/)
       walton_str%nuc(2) = 913.05397
       walton_str%aval(2) = -0.53135445
       walton_str%bval(2) = 1.0014343
       walton_str%Nsp(2) = -8.130
       walton_str%c(:,2) = (/6.13,-0.131942,0.000673193/)
       walton_str%nuc(3) = 0.
       walton_str%aval(3) = 0.
       walton_str%bval(3) = 0.
       walton_str%Nsp(3) = 0.
       walton_str%c(:,3) = (/0.,0.,0./)
    CASE(6)
       walton_str%nuc(1) = 2671.5433
       walton_str%aval(1) = -1.7667110
       walton_str%bval(1) = 1.0024428
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 913.46088
       walton_str%aval(2) = -0.50395970
       walton_str%bval(2) = 1.0013592
       walton_str%Nsp(2) = -3.26
       walton_str%c(:,2) = (/2.24,-0.03964,0.00016925/)
       walton_str%nuc(3) = 0.
       walton_str%aval(3) = 0.
       walton_str%bval(3) = 0.
       walton_str%Nsp(3) = 0.
       walton_str%c(:,3) = (/0.,0.,0./)
    CASE(7)
       walton_str%nuc(1) = 2684.5233
       walton_str%aval(1) = -1.9488269
       walton_str%bval(1) = 1.0029260
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 928.23757
       walton_str%aval(2) = -0.52807997
       walton_str%bval(2) = 1.0014039
       walton_str%Nsp(2) = -5.16
       walton_str%c(:,2) = (/5.25,-0.10217,0.0004819/)
       walton_str%nuc(3) = 841.52137
       walton_str%aval(3) = -0.40557027
       walton_str%bval(3) = 1.0011789
       walton_str%Nsp(3) = -4.28
       walton_str%c(:,3) = (/3.93,-0.06317,0.0002425/)
    CASE(8)
       walton_str%nuc(1) = 2651.3776
       walton_str%aval(1) = -1.7764105
       walton_str%bval(1) = 1.0024260
       walton_str%Nsp(1) = 0.0
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 915.30330
       walton_str%aval(2) = -0.50017997
       walton_str%bval(2) = 1.0013460
       walton_str%Nsp(2) = -3.26
       walton_str%c(:,2) = (/2.24,-0.03964,0.00016925/)
       walton_str%nuc(3) = 0.
       walton_str%aval(3) = 0.
       walton_str%bval(3) = 0.
       walton_str%Nsp(3) = 0.
       walton_str%c(:,3) = (/0.,0.,0./)
    CASE(9)
       walton_str%nuc(1) = 2690.0451
       walton_str%aval(1) = -1.8832662
       walton_str%bval(1) = 1.0028978
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 930.50230
       walton_str%aval(2) = -0.5115335
       walton_str%bval(2) = 1.001357
       walton_str%Nsp(2) = -5.53
       walton_str%c(:,2) = (/5.24,-0.1136,0.0006033/)
       walton_str%nuc(3) = 845.750
       walton_str%aval(3) = -0.388215
       walton_str%bval(3) = 1.001121
       walton_str%Nsp(3) = -3.06
       walton_str%c(:,3) = (/2.42,-0.0469,0.0002198/)
    CASE(10)
       walton_str%nuc(1) = 2672.6164
       walton_str%aval(1) = -1.7986926
       walton_str%bval(1) = 1.0026326
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 910.49626
       walton_str%aval(2) = -0.45707063
       walton_str%bval(2) = 1.0012272
       walton_str%Nsp(2) = -7.29
       walton_str%c(:,2) = (/5.76,-0.1157,0.0005882/)
       walton_str%nuc(3) = 0.
       walton_str%aval(3) = 0.
       walton_str%bval(3) = 0.
       walton_str%Nsp(3) = 0.
       walton_str%c(:,3) = (/0.,0.,0./)
    CASE(11)
       walton_str%nuc(1) = 2680.05
       walton_str%aval(1) = -1.738973
       walton_str%bval(1) = 1.003354
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 927.462
       walton_str%aval(2) = -0.321199
       walton_str%bval(2) = 1.001213
       walton_str%Nsp(2) = -8.055
       walton_str%c(:,2) = (/7.21,-0.1588,0.0008739/)
       walton_str%nuc(3) = 840.746
       walton_str%aval(3) = -0.048652
       walton_str%bval(3) = 1.000664
       walton_str%Nsp(3) = -3.51
       walton_str%c(:,3) = (/2.92,-0.0540,0.0002504/)
    CASE(12)
       walton_str%nuc(1) = 2651.7708
       walton_str%aval(1) = -1.9052739
       walton_str%bval(1) = 1.003010
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 922.36261
       walton_str%aval(2) = -0.63404209
       walton_str%bval(2) = 1.0017076
       walton_str%Nsp(2) = -5.51
       walton_str%c(:,2) = (/5.11,-0.1107,0.0005968/)
       walton_str%nuc(3) = 838.02678
       walton_str%aval(3) = -0.41086587
       walton_str%bval(3) = 1.001201
       walton_str%Nsp(3) = -2.51
       walton_str%c(:,3) = (/1.91,-0.0370,0.0001775/)
    CASE(14)
       walton_str%nuc(1) = 2654.25
       walton_str%aval(1) = -1.885330
       walton_str%bval(1) = 1.003839
       walton_str%Nsp(1) = 0.0069
       walton_str%c(:,1) = (/-0.0031,0.00359,0.0000000/)
       walton_str%nuc(2) = 928.349
       walton_str%aval(2) = -0.308384
       walton_str%bval(2) = 1.001443
       walton_str%Nsp(2) = -4.05
       walton_str%c(:,2) = (/3.72,-0.07622,0.0003822/)
       walton_str%nuc(3) = 833.040
       walton_str%aval(3) = 0.022171
       walton_str%bval(3) = 1.000538
       walton_str%Nsp(3) = -2.29
       walton_str%c(:,3) = (/2.00,-0.03806,0.000174/)
    CASE(15)
       walton_str%nuc(1) = 2695.9743
       walton_str%aval(1) = -1.624481
       walton_str%bval(1) = 1.001989
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 925.4075
       walton_str%aval(2) = -0.338243
       walton_str%bval(2) = 1.001283
       walton_str%Nsp(2) = -4.50
       walton_str%c(:,2) = (/4.76,-0.0932,0.0004524/)
       walton_str%nuc(3) = 839.8979
       walton_str%aval(3) = -0.304856
       walton_str%bval(3) = 1.000977
       walton_str%Nsp(3) = -3.61
       walton_str%c(:,3) = (/3.83,-0.0659,0.0002811/)
    CASE(16)
       walton_str%nuc(1) = 2681.2540
       walton_str%aval(1) = -1.6774586
       walton_str%bval(1) = 1.0017316
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 922.34790
       walton_str%aval(2) = -0.55636216
       walton_str%bval(2) = 1.0014921
       walton_str%Nsp(2) = -2.467
       walton_str%c(:,2) = (/2.96,-0.05411,0.00024532/)
       walton_str%nuc(3) = 834.61814
       walton_str%aval(3) = 834.61814
       walton_str%bval(3) = 1.0012166
       walton_str%Nsp(3) = -2.009
       walton_str%c(:,3) = (/2.25,-0.03665,0.00014854/)
    CASE(17)
       walton_str%nuc(1) = 2669.1414
       walton_str%aval(1) = -1.7002941
       walton_str%bval(1) = 1.0026724
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 928.29959
       walton_str%aval(2) = -0.56634758
       walton_str%bval(2) = 1.0015205
       walton_str%Nsp(2) = -8.55
       walton_str%c(:,2) = (/8.22,-0.15795,0.00075579/)
       walton_str%nuc(3) = 840.20289
       walton_str%aval(3) = -0.37264803
       walton_str%bval(3) = 1.0010841
       walton_str%Nsp(3) = -3.97
       walton_str%c(:,3) = (/4.31,-0.07318,0.00030976/)
    CASE(18)
       walton_str%nuc(1) = 2660.6468
       walton_str%aval(1) = -1.7222650
       walton_str%bval(1) = 1.0028633
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 928.73452
       walton_str%aval(2) = -0.54696239
       walton_str%bval(2) = 1.0014581
       walton_str%Nsp(2) = -5.53
       walton_str%c(:,2) = (/5.82,-0.11069,0.00052337/)
       walton_str%nuc(3) = 834.08306
       walton_str%aval(3) = -0.39938376
       walton_str%bval(3) = 1.0011724
       walton_str%Nsp(3) = -2.22
       walton_str%c(:,3) = (/2.67,-0.04360,0.00017715/)
    CASE(19)
       walton_str%nuc(1) = 2670.2425
       walton_str%aval(1) = -1.6863857
       walton_str%bval(1) = 1.0025955
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 927.92374
       walton_str%aval(2) = -0.39419031
       walton_str%bval(2) = 1.0013299
       walton_str%Nsp(2) = -5.49
       walton_str%c(:,2) = (/5.70,-0.11187,0.00054668/)
       walton_str%nuc(3) = 831.28619
       walton_str%aval(3) = -0.26364620
       walton_str%bval(3) = 1.0009546
       walton_str%Nsp(3) = -3.39
       walton_str%c(:,3) = (/3.58,-0.05991,0.00024985/)
    CASE(-1)
       walton_str%nuc(1) = 2687.0392
       walton_str%aval(1) = -2.0653147
       walton_str%bval(1) = 1.0034418
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 927.27630
       walton_str%aval(2) = -0.56503332
       walton_str%bval(2) = 1.0015090
       walton_str%Nsp(2) = -4.98
       walton_str%c(:,2) = (/5.44,-0.10152,0.00046964/)
       walton_str%nuc(3) = 837.80762
       walton_str%aval(3) = -0.38472766
       walton_str%bval(3) = 1.0011264
       walton_str%Nsp(3) = -3.40
       walton_str%c(:,3) = (/3.84,-0.06249,0.00025239/)
    CASE(-2)
       walton_str%nuc(1) = 2664.3384
       walton_str%aval(1) = -1.7711318
       walton_str%bval(1) = 1.0029931
       walton_str%Nsp(1) = 0.
       walton_str%c(:,1) = (/0.,0.,0./)
       walton_str%nuc(2) = 933.71521
       walton_str%aval(2) = -0.51860807
       walton_str%bval(2) = 1.0013778
       walton_str%Nsp(2) = -4.98
       walton_str%c(:,2) = (/5.44,-0.10152,0.00046964/)
       walton_str%nuc(3) = 839.72764
       walton_str%aval(3) = -0.40059787
       walton_str%bval(3) = 1.0011702
       walton_str%Nsp(3) = -3.40
       walton_str%c(:,3) = (/3.84,-0.06249,0.00025239/)
    CASE DEFAULT
       CALL Gbcs_Critical(.TRUE.,'Cannot find AVHRR_No',&
            'Setup_Walton','NOAA_LoadAVHRRLevel1B.f90')
    END SELECT

  END SUBROUTINE Setup_Walton

  SUBROUTINE Deallocate_Walton(walton_str)

    TYPE(Walton_Struct), INTENT(INOUT) :: walton_str

    IF( ALLOCATED(walton_str%sm_bb) )THEN
       DEALLOCATE(walton_str%sm_bb)
    ENDIF
    IF( ALLOCATED(walton_str%sm_sp) )THEN
       DEALLOCATE(walton_str%sm_sp)
    ENDIF
    IF( ALLOCATED(walton_str%sm_prt) )THEN
       DEALLOCATE(walton_str%sm_prt)
    ENDIF

  END SUBROUTINE Deallocate_Walton

  REAL FUNCTION Get_Walton(counts,inchan,pos,walton_str,out_gain)

    REAL, INTENT(IN) :: counts
    INTEGER, INTENT(IN) :: inchan
    INTEGER, INTENT(IN) :: pos
    TYPE(Walton_Struct), INTENT(IN) :: walton_str
    REAL, INTENT(OUT) :: out_gain

    INTEGER :: chan

    REAL :: R_Lin
    REAL :: R_ICT
    REAL :: ICT_T
    REAL :: Rad
    REAL :: Tstar

    ! Make to correct index (from 3->1), 4->2, 5->3)
    chan = inchan-2

    IF( walton_str%sm_prt(1,pos) .lt. -1e20 .or. &
         walton_str%sm_prt(2,pos) .lt. -1e20 .or. &
         walton_str%sm_prt(3,pos) .lt. -1e20 .or. &
         walton_str%sm_prt(4,pos) .lt. -1e20 .or. &
         walton_str%sm_sp(chan,pos) .lt. 0 .or. &
         walton_str%sm_bb(chan,pos) .lt. 0 )THEN
       Get_Walton = NAN_R
       RETURN
    ENDIF

    ! Calculate the ICT temperature from the average
    ICT_T = (walton_str%sm_prt(1,pos) + walton_str%sm_prt(1,pos) + &
         walton_str%sm_prt(1,pos) + walton_str%sm_prt(1,pos))/4.

    ! Calculate the ICT radiance using band coefficients
    Tstar = (ICT_T-walton_str%aval(chan))/walton_str%bval(chan)
    R_ICT = C1*(walton_str%nuc(chan)**3)/(exp(C2*walton_str%nuc(chan)/Tstar)-1.)

    ! Calculate linear radiance term
    out_gain = (R_ICT-walton_str%Nsp(chan))/&
         (walton_str%sm_sp(chan,pos)-walton_str%sm_bb(chan,pos))
    R_Lin = walton_str%Nsp(chan)+out_gain*(walton_str%sm_sp(chan,pos)-counts)

    ! Calculate final radiance including any non-linearity
    Rad =  R_Lin + walton_str%c(1,chan) + walton_str%c(2,chan)*R_Lin + &
         walton_str%c(3,chan)*R_Lin*R_Lin

    Get_Walton=Rad
    RETURN
    
  END FUNCTION Get_Walton

  REAL FUNCTION Get_Walton_BT(Rad,inchan,walton_str)

    REAL, INTENT(IN) :: Rad
    INTEGER, INTENT(IN) :: inchan
    TYPE(Walton_Struct), INTENT(IN) :: walton_str

    INTEGER :: chan

    REAL :: Tstar

    IF( Rad .le. 0 )THEN
       Get_Walton_BT = -1e30
       RETURN
    ENDIF

    ! Make to correct index (from 3->1), 4->2, 5->3)
    chan = inchan-2

    ! Convert back to BT using band coefficients
    Tstar = C2*walton_str%nuc(chan)/LOG((C1*(walton_str%nuc(chan)**3)/Rad)+1.)

    Get_Walton_BT = walton_str%aval(chan) + walton_str%bval(chan)*Tstar
    RETURN

  END FUNCTION Get_Walton_BT

  SUBROUTINE Calibration_Walton(outData,out_radiances,nedt_scale)

    TYPE(AVHRR_Data), INTENT(INOUT) :: outData
    LOGICAL, INTENT(IN) :: out_radiances
    REAL, INTENT(IN) :: nedt_scale

    INTEGER :: I,J
    INTEGER :: STAT
    INTEGER :: start_pos
    INTEGER :: stop_pos
    REAL :: radiance
    REAL :: ICT_T
    REAL :: Tstar
    REAL :: R_ICT
    REAL :: R_Lin
    REAL :: newGain
    TYPE(Walton_Struct) :: walton_str
    REAL, ALLOCATABLE :: Gain3(:)
    REAL :: Gain4
    REAL :: Gain5
    REAL, ALLOCATABLE :: Coefs3(:,:)
    REAL, ALLOCATABLE :: Coefs4(:,:)
    REAL, ALLOCATABLE :: Coefs5(:,:)
    LOGICAL :: replaced_data
    REAL :: radiance_error
    REAL :: error1, error2
    REAL :: newrad

    !
    ! Allocate NeDT arrays for array version
    !
    ALLOCATE(outData%array3b_error(outData%nelem,outData%arraySize),&
         outData%array4_error(outData%nelem,outData%arraySize),&
         outData%array5_error(outData%nelem,outData%arraySize),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate outData%arrayN_error',&
            'Calibration_Walton','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF
    outData%array3b_error = NAN_R
    outData%array4_error = NAN_R
    outData%array5_error = NAN_R

    !
    ! These are local arrays
    !
    ALLOCATE(Gain3(outData%arraySize),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate gain arrays',&
            'Calibration_Walton','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF
    ALLOCATE(Coefs3(3,outData%arraySize),&
         Coefs4(3,outData%arraySize),&
         Coefs5(3,outData%arraySize),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate coefs arrays',&
            'Calibration_Walton','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF

    WRITE(*,'('' ======= Using Walton calibration consistently ======='')')

    CALL Setup_Walton(outData,walton_str)
    
    !
    ! For 3.7 micron channel - need to deal with solar contamination
    !
    ! Get gain first
    !
    Gain3(:) = -1e30
    DO I=1,outData%arraySize
       IF( walton_str%sm_prt(1,I) .gt. 0 .and. &
            walton_str%sm_prt(2,I) .gt. 0 .and. &
            walton_str%sm_prt(3,I) .gt. 0 .and. &
            walton_str%sm_prt(4,I) .gt. 0 .and. &
            walton_str%sm_sp(1,I) .gt. 0 .and. &
            walton_str%sm_bb(1,I) .gt. 0 )THEN
          ! Calculate the ICT temperature from the average
          ICT_T = (walton_str%sm_prt(1,I) + walton_str%sm_prt(1,I) + &
               walton_str%sm_prt(1,I) + walton_str%sm_prt(1,I))/4.

          ! Calculate the ICT radiance using band coefficients
          Tstar = (ICT_T-walton_str%aval(1))/walton_str%bval(1)
          R_ICT = C1*(walton_str%nuc(1)**3)/(exp(C2*walton_str%nuc(1)/Tstar)-1.)
          Gain3(I) = (R_ICT-walton_str%Nsp(1))/&
               (walton_str%sm_sp(1,I)-walton_str%sm_bb(1,I))
       ENDIF
    END DO

    ! 
    ! Now use relevant solar contamination flags to reset gain if needed
    !
    IF( outData%newCalibration_There )THEN
       !
       ! Use new solar flags and rewrite operational flags
       !
       replaced_data=.FALSE.
       start_pos=-1
       stop_pos=-1
       DO I=1,outData%arraySize
          IF( outData%solar_contamination_3B(I) )THEN
             IF( start_pos .eq. -1 )THEN
                start_pos = I
                stop_pos=I
             ELSE
                stop_pos = I
             ENDIF
          ELSE
             !
             ! If finished a window of contamination
             !
             IF( 0 .lt. start_pos )THEN
                replaced_data=.TRUE.
                !
                ! Get average gain and replace
                !
                newGain = (Gain3(start_pos) + Gain3(stop_pos))/2.
                DO J=start_pos,stop_pos
                   Gain3(J) = newGain
                END DO
             ENDIF
             !
             ! Reset counters
             !
             start_pos=-1
             stop_pos=-1
          ENDIF
          outData%orig_solar_contamination_3B(I) = &
               outData%solar_contamination_3B(I)
       END DO
       IF( replaced_data )THEN
          WRITE(*,'('' === Corrected Walton Gain for solar contamination ==='')')
          WRITE(*,'('' === Replaced operational solar flags (3B) with new =='')')
       ENDIF
    ELSE
       !
       ! Use operational flags
       !
       start_pos=-1
       stop_pos=-1
       DO I=1,outData%arraySize
          IF( outData%orig_solar_contamination_3B(I) )THEN
             IF( start_pos .eq. -1 )THEN
                start_pos = I
                stop_pos = I
             ELSE
                stop_pos = I
             ENDIF
          ELSE
             !
             ! If finished a window of contamination
             !
             IF( 0 .lt. start_pos )THEN
                !
                ! Get average gain and replace
                !
                newGain = (Gain3(start_pos) + Gain3(stop_pos))/2.
                DO J=start_pos,stop_pos
                   Gain3(J) = newGain
                END DO
             ENDIF
             !
             ! Reset counters
             !
             start_pos=-1
             stop_pos=-1
          ENDIF
       END DO
    ENDIF

    !
    ! Now use gain to get BTs
    !
    DO I=1,outData%arraySize
       IF( -1e20 .lt. Gain3(I) )THEN
          DO J=1,outData%nelem
             IF( outData%array3b(J,I) .ge. 0 )THEN
                R_Lin = walton_str%Nsp(1)+Gain3(I)*&
                     (walton_str%sm_sp(1,I)-outData%counts3(J,I))
                radiance =  R_Lin + walton_str%c(1,1) + &
                     walton_str%c(2,1)*R_Lin + &
                     walton_str%c(3,1)*R_Lin*R_Lin
                !
                ! Calculate operational coefficients
                !
                Coefs3(1,I) = (walton_str%c(1,1)+&
                     (1+walton_str%c(2,1))*walton_str%Nsp(1)+&
                     walton_str%c(2,1)*walton_str%Nsp(1)**2) + &
                     (1+walton_str%c(2,1)+2.*walton_str%c(3,1)*&
                     walton_str%Nsp(1))*Gain3(I)*walton_str%sm_sp(1,I) + &
                     walton_str%c(3,1)*Gain3(I)**2*walton_str%sm_sp(1,I)**2
                Coefs3(2,I) = -((1+walton_str%c(2,1)+2.*walton_str%c(3,1)*&
                     walton_str%Nsp(1))*Gain3(I) + &
                     2*walton_str%c(3,1)*Gain3(I)**2*&
                     walton_str%sm_sp(1,I))
                Coefs3(3,I) = walton_str%c(3,1)*Gain3(I)**2
                if( out_radiances )then
                   outData%array3B(J,I) = radiance                
                   outData%array3B_error(J,I) = ABS(Coefs3(2,I)+&
                        2*Coefs3(3,I)*outData%counts3(J,I))*&
                        outData%noise_cnts(4,I)
                else
                   outData%array3B(J,I) = &
                        get_walton_BT(radiance,3,walton_str)
                   radiance_error = ABS(Coefs3(2,I)+&
                        2*Coefs3(3,I)*outData%counts3(J,I))*&
                        outData%noise_cnts(4,I)
                   newrad = radiance-radiance_error
                   if( 0. .lt. newrad )then
                      error1 = ABS(outData%array3B(J,I) - &
                           get_walton_BT(newrad,3,walton_str))
                   else
                      error1 = 10
                   endif
                   newrad = radiance+radiance_error
                   error2 = ABS(outData%array3B(J,I) - &
                        get_walton_BT(newrad,3,walton_str))
                   outData%array4_error(J,I) = MAX(error1,error2)
                endif
                IF( outData%isGAC )THEN
                   IF( 0 .lt. outData%array3B_error(J,I) )THEN
                      IF( outData%counts3(J,I) .ge. 3 )THEN
                         outData%array3B_error(J,I) = &
                              outData%array3B_error(J,I)/2.
                      ELSE
                         outData%array3B_error(J,I) = &
                              outData%array3B_error(J,I)/nedt_scale
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
          END DO
       ENDIF
    END DO

    !
    ! Now 11/12 micron - no solar contamination correction
    !
    DO I=1,outData%arraySize
       DO J=1,outData%nelem

          ! Channel 4
          IF( outData%array4(J,I) .ge. 0 )THEN
             radiance = Get_Walton(outData%counts4(J,I),&
                  4,I,walton_str,Gain4)
             !
             ! Calculate operational coefficients
             !
             Coefs4(1,I) = (walton_str%c(1,2)+&
                  (1+walton_str%c(2,2))*walton_str%Nsp(2)+&
                  walton_str%c(2,2)*walton_str%Nsp(2)**2) + &
                  (1+walton_str%c(2,2)+2.*walton_str%c(3,2)*&
                  walton_str%Nsp(2))*Gain4*walton_str%sm_sp(2,I) + &
                  walton_str%c(3,2)*Gain4**2*walton_str%sm_sp(2,I)**2
             Coefs4(2,I) = -((1+walton_str%c(2,2)+2.*walton_str%c(3,2)*&
                  walton_str%Nsp(2))*Gain4 + &
                  2*walton_str%c(3,2)*Gain4**2*&
                  walton_str%sm_sp(2,I))
             Coefs4(3,I) = walton_str%c(3,2)*Gain4**2
             if( out_radiances )then
                outData%array4(J,I) = radiance                
                outData%array4_error(J,I) = ABS(Coefs4(2,I)+&
                     2*Coefs4(3,I)*outData%counts4(J,I))*&
                     outData%noise_cnts(5,I)
             else
                outData%array4(J,I) = get_walton_BT(radiance,4,walton_str)
                radiance_error = ABS(Coefs4(2,I)+&
                     2*Coefs4(3,I)*outData%counts4(J,I))*&                     
                     outData%noise_cnts(5,I)
                newrad = radiance-radiance_error
                if( 0. .lt. newrad )then
                   error1 = ABS(outData%array4(J,I) - &
                        get_walton_BT(newrad,4,walton_str))
                else
                   error1 = 10
                endif
                newrad = radiance+radiance_error
                error2 = ABS(outData%array4(J,I) - &
                     get_walton_BT(newrad,4,walton_str))
                newrad = get_walton_BT(radiance+radiance_error,4,walton_str)
                newrad = get_walton_BT(radiance-radiance_error,4,walton_str)
                outData%array4_error(J,I) = MAX(error1,error2)
             endif
             IF( outData%isGAC )THEN
                IF( 0 .lt. outData%array4_error(J,I) )THEN
                   IF( outData%counts4(J,I) .ge. 3 )THEN
                      outData%array4_error(J,I) = &
                           outData%array4_error(J,I)/2.
                   ELSE
                      outData%array4_error(J,I) = &
                           outData%array4_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             ENDIF
          ENDIF

          ! Check if we have 12 micron data
          IF( outData%array5(J,I) .ge. 0 )THEN
             radiance = Get_Walton(outData%counts5(J,I),5,i,&
                  walton_str,Gain5)
             Coefs5(1,I) = (walton_str%c(1,3)+&
                  (1+walton_str%c(2,3))*walton_str%Nsp(3)+&
                  walton_str%c(2,3)*walton_str%Nsp(3)**2) + &
                  (1+walton_str%c(2,3)+2.*walton_str%c(3,3)*&
                  walton_str%Nsp(3))*Gain5*walton_str%sm_sp(3,I) + &
                  walton_str%c(3,3)*Gain5**2*walton_str%sm_sp(3,I)**2
             Coefs5(2,I) = -((1+walton_str%c(2,3)+2.*walton_str%c(3,3)*&
                  walton_str%Nsp(3))*Gain5 + &
                  2*walton_str%c(3,3)*Gain5**2*&
                  walton_str%sm_sp(3,I))
             Coefs5(3,I) = walton_str%c(3,3)*Gain5**2
             if( out_radiances )then
                outData%array5(J,I) = radiance                
                outData%array5_error(J,I) = ABS(Coefs5(2,I)+&
                     2*Coefs5(3,I)*outData%counts5(J,I))*&
                     outData%noise_cnts(6,I)
             else
                outData%array5(J,I) = get_walton_BT(radiance,5,walton_str)
                radiance_error = ABS(Coefs5(2,I)+&
                     2*Coefs5(3,I)*outData%counts5(J,I))*&
                     outData%noise_cnts(6,I)
                newrad = radiance-radiance_error
                if( 0. .lt. newrad )then
                   error1 = ABS(outData%array5(J,I) - &
                        get_walton_BT(newrad,5,walton_str))
                else
                   error1 = 10
                endif
                newrad = radiance+radiance_error
                error2 = ABS(outData%array5(J,I) - &
                     get_walton_BT(newrad,5,walton_str))
                outData%array5_error(J,I) = MAX(error1,error2)
             ENDIF
             IF( outData%isGAC )THEN
                IF( 0 .lt. outData%array5_error(J,I) )THEN
                   IF( outData%counts5(J,I) .ge. 3 )THEN
                      outData%array5_error(J,I) = &
                           outData%array5_error(J,I)/2.
                   ELSE
                      outData%array5_error(J,I) = &
                           outData%array5_error(J,I)/nedt_scale
                   ENDIF
                ENDIF
             ENDIF
          ENDIF

       END DO
    END DO

    CALL Deallocate_Walton(walton_str)
    DEALLOCATE(Gain3)
    DEALLOCATE(Coefs3,Coefs4,Coefs5)

    outData%walton_there = .TRUE.

  END SUBROUTINE Calibration_Walton

  SUBROUTINE Parse_RecordV1(Coefs,DataType,scanLinePos,sizeOfRecords,&
       record,outData,bad_data,out_radiances,have_12micron,&
       getTiePoints,getLandMask,getLocation,printAllErrors,&
       ignore_solar_contamination)

    TYPE(AVHRR_Instrument_Coefs), INTENT(IN) :: Coefs
    INTEGER, INTENT(IN) :: DataType
    INTEGER, INTENT(IN) :: scanLinePos
    INTEGER, INTENT(IN) :: sizeOfRecords
    INTEGER(GbcsInt1), INTENT(IN) :: record(sizeOfRecords)
    TYPE(AVHRR_Data), POINTER :: outData
    TYPE(AVHRR_Bad_Data), INTENT(IN) :: bad_data
    LOGICAL, INTENT(IN) :: out_radiances
    LOGICAL, INTENT(IN) :: have_12micron
    LOGICAL, OPTIONAL :: getTiePoints
    LOGICAL, OPTIONAL :: getLandMask
    LOGICAL, OPTIONAL :: getLocation
    LOGICAL, OPTIONAL :: printAllErrors
    LOGICAL, OPTIONAL :: ignore_solar_contamination

    INTEGER :: i,j,k
    LOGICAL :: badData
    LOGICAL :: bad_Top
    LOGICAL :: bad_Navigation
    LOGICAL :: bad_Time
    LOGICAL :: bad_Calibration
    INTEGER(GbcsInt2) :: shortVal
    INTEGER :: scanLineNumber
    INTEGER(GbcsInt1) :: qualityTop(4)
    INTEGER(GbcsInt2) :: clockDriftField
    LOGICAL :: timeCorrected
    INTEGER(GbcsInt2) :: clockDrift_msecs
    INTEGER :: frameDayno, frameUTC_msecs

    INTEGER :: utc
    INTEGER :: intVal
    INTEGER :: stat

    REAL :: calib1(2)
    REAL :: calib1_2(2)
    REAL :: calib1_intercept
    REAL :: calib2(2)
    REAL :: calib2_2(2)
    REAL :: calib2_intercept
    REAL :: calib3A(2)
    REAL :: calib3A_2(2)
    REAL :: calib3A_intercept
    REAL :: calib3(3)
    REAL :: calib4(3)
    REAL :: calib5(3)

    REAL :: satelliteAltitude
    REAL :: solarAngleIn(51)
    REAL :: satelliteAngleIn(51)
    REAL :: relativeAngleIn(51)
    REAL :: latitudeIn(51)
    REAL :: longitudeIn(51)
    REAL :: d

    INTEGER :: pos
    REAL :: satellite_long
    REAL :: satellite_lat
    REAL :: solarAngle(NDATA_MAX)
    REAL :: relativeAngle(NDATA_MAX)
    REAL :: latitude(NDATA_MAX)
    REAL :: longitude(NDATA_MAX)

    INTEGER :: prtPosition = 0

    INTEGER :: word
    INTEGER :: earthCounts(5,NDATA_MAX)
    INTEGER :: counts
    REAL :: radiance

    REAL :: minVal
    INTEGER :: minPos
    REAL :: minPosArray(4)

    REAL :: diff

    INTEGER :: nMean3
    INTEGER :: nMean4
    INTEGER :: nMean5

    REAL :: prtCounts
    REAL :: prtCounts3
    REAL :: prtCounts4
    REAL :: prtCounts5
    REAL :: prtMean

    INTEGER :: solar_start
    INTEGER :: solar_dec

    LOGICAL :: solar_3b
    LOGICAL :: solar_4
    LOGICAL :: solar_5
    LOGICAL :: print_All_Errors

    LOGICAL :: ignore_solar_contam

    ! Stores 10 bit words from minor frame
    INTEGER(GbcsInt2) :: minor_frame(103)

    print_All_Errors = .FALSE.
    IF( PRESENT(printAllErrors) )THEN
       print_All_Errors = printAllErrors
    ENDIF

    ignore_solar_contam = .FALSE.
    IF( PRESENT(ignore_solar_contamination) )THEN
       ignore_solar_contam = ignore_solar_contamination
    ENDIF

    if( .not. AVHRR_setup )then
       if( AVHRR_GAC .eq. DataType )then
          ndata = NDATA_GAC
          ndata_counts = 682
          do i=1,4
             prtRotate(i) = prtRotate_GAC(i)
          end do
       else if( AVHRR_LAC .eq. DataType .or. AVHRR_FRAC .eq. DataType )then
          ndata = NDATA_LAC
          ndata_counts = 3414
          do i=1,4
             prtRotate(i) = prtRotate_LAC(i)
          end do
       else
          call out_ERROR('Data not GAC/LAC',&
               'Parse_RecordV1','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif
       AVHRR_Setup = .TRUE.
       ! Allocate output arrays
       call Allocate_OutData( ndata, outData )
       outData%badTime = .TRUE.
       outData%badNavigation = .TRUE.
       outData%badCalibration = .TRUE.
    endif

    ! check to see if we have enought storage
    if( scanLinePos .gt. outData%arraySize )then
       call Reallocate_outData( outData, MEMORY_ALLOC_STEP )
    endif

    ! assume bad to start with
    outData%badTime(scanLinePos) = .TRUE.
    outData%badNavigation(scanLinePos) = .TRUE.
    outData%badCalibration(scanLinePos) = .TRUE.

    ! Get scan line number
    scanLineNumber = header2ushort( 0, record )
    outData%scanLineNumber(scanLinePos) = scanLineNumber

    ! Check quality
    qualityTop = record(9:12)
    minor_frame = get_hrpt_frame(record(309:448))

    ! this field is for info only and is only valid on/after Nov 16 1994: test
    ! validClockDriftInfo (assigned by checkDate after this subroutine is called)
    ! no drift info was recorded nor corrections applied before Aug 29 2000
    ! but I've found that N14 1996d181, 1997d001, ...,
    ! has timeCorrected going ...TFTFTFTF... changing at each scan line and
    ! with clockDrift_msecs rather large = 4265, but the frame time is the same
    ! as the scanline time so no correction applied .
    clockDriftField = header2short( 3196, record ) ! clock drift * 2 + 0/1
    timeCorrected = BTEST(clockDriftField,0)
    ! clockDrift_msecs = UTC - spacecraft time (the other way round from the docs)
    ! i.e. this is what you add to the spacecraft clock to get correct time
    clockDrift_msecs = clockDriftField/2 ! integer division removes last bit
    ! note that v1 data has frame time as the 'raw' spacecraft time and the scanline
    ! time as the corrected time (cf v2 and v3!)

    badData = .FALSE.
    bad_Top = .FALSE.
    bad_Navigation = .FALSE.
    bad_Time = .FALSE.
    bad_Calibration = .FALSE.
    bad_Navigation = .FALSE.
    outData%transition3a3b(scanLinePos) = .TRUE.

    ! Get time
    CALL Parse_Date_Time( record, 3, outData%year(scanLinePos), &
         outData%dayno(scanLinePos), outData%hours(scanLinePos), outData%UTC_msecs(scanLinePos) )
    CALL Parse_DayNo( outData%year(scanLinePos), outData%dayno(scanLinePos),&
         outData%month(scanLinePos), outData%day(scanLinePos) )
    outData%time(scanLinePos) = get_date_time( outData%year(scanLinePos),&
        outData%dayno(scanLinePos),outData%hours(scanLinePos), time_yearstart )

!!$    PRINT'(2i6,1i20)',scanLinePos,scanLineNumber,outData%UTC_msecs(scanLinePos)
!!$    PRINT'(2i6,4b9.8,1b11.10)',scanLinePos,scanLineNumber,&
!!$         qualityTop,&
!!$         IBITS( minor_frame(7), 0, 10 )

    bad_Top = CheckTopV1(qualityTop, minor_frame)
    bad_Navigation = CheckNavigationV1(qualityTop)
    ! this should be caught by bad_Top but it may not always do it
    IF( 51 .NE. record(53) ) bad_Navigation = .TRUE.
    ! if using clavr-x lon/lat then need to check those
    IF (.NOT.useL1BClock) THEN
       IF (ANY(clavrLonAnchor(:,scanLinePos).EQ.NAN_R) &
            .OR.ANY(clavrLatAnchor(:,scanLinePos).EQ.NAN_R))THEN
          bad_Navigation = .TRUE.
       ELSE
          ! check the nadir element
          i = 26
          latitudeIn(i) = header2short((100+i*4),Record)/128.
          longitudeIn(i) = header2short((102+i*4),Record)/128.
          CALL greatCircleDistance(longitudeIn(i),latitudeIn(i),clavrLonAnchor(i,scanLinePos),clavrLatAnchor(i,scanLinePos),d)
          IF(d.GT.maxTimingErrorDistance) THEN
             bad_Navigation = .TRUE.
          END IF
       END IF
    END IF
    bad_Time = CheckTimeV1(qualityTop)
    bad_Calibration = CheckCalibrationV1(qualityTop, solar_3b, solar_4, solar_5)

    ! too bad to use at all
    IF( bad_Top )THEN
       IF( 0 .EQ. printOnce(1) .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (top/frame errors) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       ENDIF
       printOnce(1) = 1
       IF( bad_data%allow_bad_top )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_outdata_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    ! bad time means bad navigation
    if( bad_Time )then
       if( 0 .eq. printOnce(3) .or. print_All_Errors )then
          write(*,'(''Quality Problems (time) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
       endif
       printOnce(3) = 1
       if( bad_data%allow_bad_time )then
          badData = .TRUE.
       else
          prtNumber = -1
          call set_outdata_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badTime(scanLinePos) = .FALSE.
    endif

    IF( bad_Time )THEN
       CALL set_time_bad( outData, scanLinePos )
    ENDIF

    ! CheckNavigationV1 also checks the time
    if( bad_Navigation )then
       if( 0 .eq. printOnce(2) .or. print_All_Errors )then
          write(*,'(''Quality Problems (navigation) scan line, position '',2i7)')&
               scanLineNumber,scanLinePos
          if( bad_data%use_bad_nav )then
             write(*,'(''  Using potentially bad navigation '')')
          endif
       endif
       printOnce(2) = 1
       if( bad_data%allow_bad_nav .or. bad_data%use_bad_nav )then
          badData = .TRUE.
          if( bad_data%use_bad_nav )then
             ! If use navigation then it is not considered 'bad'
             bad_navigation = .FALSE.
             outdata%badNavigation(scanLinePos) = .FALSE.
          endif
       else
          prtNumber = -1
          ! everything except time
          call set_nav_bad( outData, scanLinePos )
          call set_calib_bad( outData, scanLinePos )
          return
       endif
    else
       outdata%badNavigation(scanLinePos) = .FALSE.
    endif

    if( .not. bad_Navigation )then
       ! Position and angles
       if( AVHRR_GAC .eq. DataType )then
          solar_start = 54
          solar_dec = 3177
       else
          solar_start = 54
          solar_dec = 6705
       endif
       do i=1,51
          ! Deal with decimal part which is stored elsewhere
          solarAngleIn(i) = Get_Solar_Angle(Record,solar_start,solar_dec,i)
          ! Note that the start numbers are offset bu i*6 or i*8 from the
          ! record numbers
          latitudeIn(i) = header2short((100+i*4),Record)/128.
          longitudeIn(i) = header2short((102+i*4),Record)/128.
       end do

       IF (.NOT.useL1BClock) THEN
          ! use the clavr-x nav lon/lat that has been read in
          longitudeIn = clavrLonAnchor(:,scanLinePos)
          latitudeIn = clavrLatAnchor(:,scanLinePos)
       ENDIF

       CALL InterpolateData(solarAngleIn,solarAngle,dataType)

       CALL InterpolateLonLat(longitudeIn,latitudeIn,longitude(1:ndata),latitude(1:ndata),dataType)

       pos = ndata/2 ! probably OK for pre-KLM with a GAC pixel centred on nadir
       satellite_long = longitude(pos)
       satellite_lat = latitude(pos)
       do i=1,ndata
          ! Get satellite ZA
          outData%satZA(i,scanLinePos) = get_satellite_za( longitude(i),&
               latitude(i),satellite_long,satellite_lat)
       end do

       ! Copy to output
       outData%Lon(:,scanLinePos) = longitude(1:ndata)
       outData%Lat(:,scanLinePos) = latitude(1:ndata)
       outData%solZA(:,scanLinePos) = solarAngle(1:ndata)
       outData%relAz(:,scanLinePos) = NAN_R
    else
       call set_nav_bad( outData, scanLinePos )
    endif

    IF (PRESENT(getTiePoints)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLandMask)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF (PRESENT(getLocation)) THEN
       IF( ANY(printOnce.EQ.1) )THEN
          printOnce = (/0,0,0,0,0,0,0,0,0,0/)
          WRITE(*,'(''End of quality problems on scan line, position'',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       outData%dataFilled = .TRUE.
       RETURN
    ENDIF

    IF( bad_Calibration )THEN
       IF( 0 .EQ. printonce(4)  .OR. print_All_Errors )THEN
          WRITE(*,'(''Quality Problems (calibration) scan line, position '',2i7)')&
               scanLineNumber, scanLinePos
       ENDIF
       printonce(4) = 1
       IF( bad_data%allow_bad_calib )THEN
          badData = .TRUE.
       ELSE
          prtNumber = -1
          CALL set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    ENDIF

    if( .not. badData )then
       badDataLoop: do i=1,10
          if( 1 .eq. printOnce(i) )then
             printOnce = (/0,0,0,0,0,0,0,0,0,0/)
             write(*,'(''End of quality problems on scan line, position'',2i7)')&
                  scanLineNumber, scanLinePos
             exit badDataLoop
          endif
       end do badDataLoop
    endif

    ! Solar flags
    outData%orig_solar_contamination_3b(scanLinePos) = solar_3b
    outData%orig_solar_contamination_4(scanLinePos) = solar_4
    outData%orig_solar_contamination_5(scanLinePos) = solar_5

    ! Get calibration values - slope/intercept
    calib1(2) = header2int(12,Record)/(2.**30)
    calib1(1) = header2int(16,Record)/(2.**22)
    ! Channel 2
    calib2(2) = header2int(20,Record)/(2.**30)
    calib2(1) = header2int(24,Record)/(2.**22)
    ! Infrared - native Level 1B calibration
    calib3(2) = header2int(28,Record)/(2.**30)
    calib3(1) = header2int(32,Record)/(2.**22)
    calib4(2) = header2int(36,Record)/(2.**30)
    calib4(1) = header2int(40,Record)/(2.**22)
    calib5(2) = header2int(44,Record)/(2.**30)
    calib5(1) = header2int(48,Record)/(2.**22)

    outData%calib1(:,scanLinePos) = calib1
    outData%calib2(:,scanLinePos) = calib2
    outData%calib3(:,scanLinePos) = calib3
    outData%calib4(:,scanLinePos) = calib4
    outData%calib5(:,scanLinePos) = calib5

    ! Nominal satelliteAltitude
    outData%satelliteAltitude(scanLinePos) = 833

    ! Now get calibration information for re-calibration
    ! From HRPT Minor frame
    prtCounts3 = minor_frame(18)
    prtCounts4 = minor_frame(19)
    prtCounts5 = minor_frame(20)
!---Modif du 12 nov---
!on veut imprimer les prt counts pour implementer la methode de Trishchenko
    !outData%prtcounts(1,scanLinePos)=prtCounts3
    !outData%prtcounts(2,scanLinePos)=prtCounts4
    !outData%prtcounts(3,scanLinePos)=prtCounts5
    
!-------------------

    ! Need to get which PRT we are at
    if( -1 .eq. prtNumber )then
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
       endif
       storedPrtTemp = 0.
    else
       if( 10 .gt. prtCounts3 .and. 10 .gt. prtCounts4 .and. &
            10 .gt. prtCounts5 )then
          prtNumber = 1
          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtCountsPrevStore(i) = prtCountsStore(i)
             else
                prtCountsPrevStore(i) = 0
             endif
          end do
       else
          prtCounts = (prtCounts3 + prtCounts4 + prtCounts5)/3.
          if( 1 .gt. prtNumber .or. 4 .lt. prtNumber )then
             ! Rejig array

             ! To be honest I can't remember why I wrote the following code
             ! should never really be used...

             if( 0 .lt. prtCountsPrevStore(1) .and. &
                  0 .lt. prtCountsPrevStore(2) .and. &
                  0 .lt. prtCountsPrevStore(3) .and. &
                  0 .lt.  prtCountsPrevStore(4) )then
                minPos = 0
                minVal = 1e32
                do j=1,4
                   diff = abs(prtCounts - prtCountsPrevStore(j))
                   if( minVal > diff )then
                      minVal = diff
                      minPos = j
                   endif
                end do
                minPosArray(4) = minPos
                do i=1,4
                   minPos = 0
                   minVal = 1e32
                   do j=1,4
                      diff = abs(prtCountsStore(i) - prtCountsPrevStore(j))
                      if( minVal > diff )then
                         minVal = diff
                         minPos = j
                      endif
                   end do
                   minPosArray(i) = minPos
                end do

                ! Now have array of indexes to matching entries - remap these
                ! ones
                ! Find repeat
                minPos = -1
                do i=2,4
                   if( minPosArray(i-1) .eq. minPosArray(i) )then
                      minPos = i
                   endif
                end do

                if( -1 .eq. minPos )then
                   ! Assume line is bad (for some reason) and set bad
                   prtNumber = -1
                   call set_calib_bad( outData, scanLinePos )
                   return
                else

                   ! Shift down
                   do i=minPos,3
                      minPosArray(i) = minPosArray(i+1)
                   end do

                   ! Set extra value to last value
                   minPosArray(4) = prtCounts
                endif

                write(*,*)'Resetting prtCounts array'

             else
                write(*,*)'No previous prt record - just taking first 4'
                prtNumber=1
             endif

          else
             prtCountsStore(prtNumber) = prtCounts
             prtCountsStoreAll(1,prtNumber) = prtCounts3
             prtCountsStoreAll(2,prtNumber) = prtCounts4
             prtCountsStoreAll(3,prtNumber) = prtCounts5
             prtNumber = prtNumber + 1
          endif

          do i=1,4
             if( 0 .lt. prtCountsStore(i) )then
                prtPosition = prtRotate(i)
                if( 1 .gt. prtPosition .or. 4 .lt. prtPosition )then
                   call out_ERROR('prtPosition out of range',&
                        'Parse_RecordV1','NOAA_LoadAVHRRLevel1B.f90',&
                        .FALSE.,' ')
                endif
                storedPrtTemp(prtPosition) = &
                     Coefs%prtTempCoefs1(prtPosition) + &
                     Coefs%prtTempCoefs2(prtPosition)*prtCountsStore(i) + &
                     Coefs%prtTempCoefs3(prtPosition)*(prtCountsStore(i)**2) + &
                     Coefs%prtTempCoefs4(prtPosition)*(prtCountsStore(i)**3) + &
                     Coefs%prtTempCoefs5(prtPosition)*(prtCountsStore(i)**4)
             endif
          end do
       end if
    end if

    if( storedPrtTemp(1) .gt. 200 .and. &
         storedPrtTemp(2) .gt. 200 .and. &
         storedPrtTemp(3) .gt. 200 .and. &
         storedPrtTemp(4) .gt. 200 )then
       outData%prt1CountsAll(:,scanLinePos) = storedPrtCountsAll(:,1)
       outData%prt2CountsAll(:,scanLinePos) = storedPrtCountsAll(:,2)
       outData%prt3CountsAll(:,scanLinePos) = storedPrtCountsAll(:,3)
       outData%prt4CountsAll(:,scanLinePos) = storedPrtCountsAll(:,4)
       outData%prt1Counts(scanLinePos) = storedPrtCounts(1)
       outData%prt2Counts(scanLinePos) = storedPrtCounts(2)
       outData%prt3Counts(scanLinePos) = storedPrtCounts(3)
       outData%prt4Counts(scanLinePos) = storedPrtCounts(4)
       outData%prt1(scanLinePos) = storedPrtTemp(1)
       outData%prt2(scanLinePos) = storedPrtTemp(2)
       outData%prt3(scanLinePos) = storedPrtTemp(3)
       outData%prt4(scanLinePos) = storedPrtTemp(4)
       prtMean = (storedPrtTemp(1) + storedPrtTemp(2) + &
            storedPrtTemp(3) + storedPrtTemp(4))/4.
    else
       outData%prt1CountsAll(:,scanLinePos) = NAN_R
       outData%prt2CountsAll(:,scanLinePos) = NAN_R
       outData%prt3CountsAll(:,scanLinePos) = NAN_R
       outData%prt4CountsAll(:,scanLinePos) = NAN_R
       outData%prt1Counts(scanLinePos) = NAN_R
       outData%prt2Counts(scanLinePos) = NAN_R
       outData%prt3Counts(scanLinePos) = NAN_R
       outData%prt4Counts(scanLinePos) = NAN_R
       outData%prt1(scanLinePos) = NAN_R
       outData%prt2(scanLinePos) = NAN_R
       outData%prt3(scanLinePos) = NAN_R
       outData%prt4(scanLinePos) = NAN_R
       prtMean = NAN_R
    endif

    ! Get BB counts
    outData%bbodyFilter3(:,scanLinePos) = 0
    outData%bbodyFilter4(:,scanLinePos) = 0
    outData%bbodyFilter5(:,scanLinePos) = 0
    outData%bb3(scanLinePos) = 0.
    outData%bb4(scanLinePos) = 0.
    outData%bb5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       ! Stored as 10 buit numbers
       outData%bbodyFilter3(i,scanLinePos) = minor_frame(23+(i-1)*3)
       outData%bbodyFilter4(i,scanLinePos) = minor_frame(24+(i-1)*3)
       outData%bbodyFilter5(i,scanLinePos) = minor_frame(25+(i-1)*3)
       ! in fact outData%bbodyFilter? are always +ve
       if( 0 .lt. outData%bbodyFilter3(i,scanLinePos) )then
          outData%bb3(scanLinePos) = outData%bb3(scanLinePos) + &
               outData%bbodyFilter3(i,scanLinePos)
          nMean3 = nMean3 + 1
       endif
       if( 0 .lt. outData%bbodyFilter4(i,scanLinePos) )then
          outData%bb4(scanLinePos) = outData%bb4(scanLinePos) + &
               outData%bbodyFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 0 .lt. outData%bbodyFilter5(i,scanLinePos) )then
          outData%bb5(scanLinePos) = outData%bb5(scanLinePos) + &
               outData%bbodyFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( NBBODY_SAMPLES .ne. nMean3 )then
       write(*,*)'Dropped Black-body (3.7) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean3 )then
       outData%bb3(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (3.7) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb3(scanLinePos) = outData%bb3(scanLinePos) / nMean3
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Black-body (11) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean4 )then
       outData%bb4(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb4(scanLinePos) = outData%bb4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Black-body (12) counts'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    endif
    if( 0 .eq. nMean5 )then
       outData%bb5(scanLinePos) = NAN_R
       write(*,*)'Black-body counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%bb5(scanLinePos) = outData%bb5(scanLinePos) / nMean5
    end if

    ! Get space counts
    outData%spaceFilter1(:,scanLinePos) = 0
    outData%spaceFilter2(:,scanLinePos) = 0
    outData%spaceFilter3a(:,scanLinePos) = 0
    outData%spaceFilter3(:,scanLinePos) = 0
    outData%spaceFilter4(:,scanLinePos) = 0
    outData%spaceFilter5(:,scanLinePos) = 0
    outData%sp3(scanLinePos) = 0.
    outData%sp4(scanLinePos) = 0.
    outData%sp5(scanLinePos) = 0.
    nMean3 = 0
    nMean4 = 0
    nMean5 = 0
    do i=1,NBBODY_SAMPLES
       outData%spaceFilter1(i,scanLinePos) = minor_frame(53+(i-1)*5)
       outData%spaceFilter2(i,scanLinePos) = minor_frame(54+(i-1)*5)
       outData%spaceFilter3(i,scanLinePos) = minor_frame(55+(i-1)*5)
       outData%spaceFilter3a(i,scanLinePos) = NAN_I
       outData%spaceFilter4(i,scanLinePos) = minor_frame(56+(i-1)*5)
       outData%spaceFilter5(i,scanLinePos) = minor_frame(57+(i-1)*5)
       ! in fact outData%spaceFilter? are always +ve
       if( 900 .lt. outData%spaceFilter3(i,scanLinePos) )then
          outData%sp3(scanLinePos) = outData%sp3(scanLinePos) + &
               outData%spaceFilter3(i,scanLinePos)
          nMean3 = nMean3 + 1
       endif
       if( 900 .lt. outData%spaceFilter4(i,scanLinePos) )then
          outData%sp4(scanLinePos) = outData%sp4(scanLinePos) + &
               outData%spaceFilter4(i,scanLinePos)
          nMean4 = nMean4 + 1
       endif
       if( 900 .lt. outData%spaceFilter5(i,scanLinePos) )then
          outData%sp5(scanLinePos) = outData%sp5(scanLinePos) + &
               outData%spaceFilter5(i,scanLinePos)
          nMean5 = nMean5 + 1
       endif
    end do
    if( NBBODY_SAMPLES .ne. nMean3 )then
       write(*,*)'Dropped Space (3.7) counts'
       IF( nMean3 .lt. 6 )THEN
          write(*,*)'Number Dropped Space (3.7) counts > 4 - setting to bad'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    endif
    if( 0 .eq. nMean3 )then
       outData%sp3(scanLinePos) = NAN_R
       write(*,*)'Space counts (3.7) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp3(scanLinePos) = outData%sp3(scanLinePos) / nMean3
    end if
    if( NBBODY_SAMPLES .ne. nMean4 )then
       write(*,*)'Dropped Space (11) counts'
       IF( nMean4 .lt. 6 )THEN
          write(*,*)'Number Dropped Space (11) counts > 4 - setting to bad'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    endif
    if( 0 .eq. nMean4 )then
       outData%sp4(scanLinePos) = NAN_R
       write(*,*)'Space counts (11) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp4(scanLinePos) = outData%sp4(scanLinePos) / nMean4
    end if
    if( NBBODY_SAMPLES .ne. nMean5 )then
       write(*,*)'Dropped Space (12) counts'
       IF( nMean5 .lt. 6 )THEN
          write(*,*)'Number Dropped Space (12) counts > 4 - setting to bad'
          call set_calib_bad( outData, scanLinePos )
          RETURN
       ENDIF
    endif
    if( 0 .eq. nMean5 )then
       outData%sp5(scanLinePos) = NAN_R
       write(*,*)'Space counts (12) are zero'
       call set_calib_bad( outData, scanLinePos )
       RETURN
    else
       outData%sp5(scanLinePos) = outData%sp5(scanLinePos) / nMean5
    end if

    ! L1B file calibration is OK but if ignoring solar contam data, ignore
    IF( ignore_solar_contam .and. (solar_3b .or. solar_4 .or. solar_5) )THEN
       outData%badCalibration(scanLinePos) = .TRUE.
    ELSE
       outData%badCalibration(scanLinePos) = .FALSE.
    ENDIF

    if( .not. bad_Navigation .and. &
         .not.outData%badCalibration(scanLinePos) )then

       ! Get earth observations
       j = 0
       k = 0
       do i=0,ndata_counts-1
          ! Data is interlaced
          word = header2int((448+i*4),record)
          ! Get bits 20-29
          earthCounts(k+1,j+1) = IBITS(word,20,10)
          k = k+1
          if( 5 .eq. k)then
             j = j + 1
             k = 0
          end if
          if( 682 .eq. ndata_counts )then
             ! Get bits 10-19
             earthCounts(k+1,j+1) = IBITS(word,10,10)
             k = k+1
             if( 5 .eq. k)then
                j = j + 1
                k = 0
             end if
             ! Get bits 0-9 if i != 681
             if( ndata_counts-1 .ne. i )then
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          else
             ! Get bits 10-19 && 0-9 if i != 3413
             if( ndata_counts-1 .ne. i )then
                ! Get bits 10-19
                earthCounts(k+1,j+1) = IBITS(word,10,10)
                k = k + 1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
                ! Get bits 0-9
                earthCounts(k+1,j+1) = IBITS(word,0,10)
                k = k+1
                if( 5 .eq. k)then
                   j = j + 1
                   k = 0
                end if
             end if
          end if
       end do

       if( 0 .ne. k .and. j .ne. ndata )then
          call out_ERROR('Problems parsing earth observations',&
               'Parse_RecordV1','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif

       do i=1,ndata
          ! Visible
          ! Channel 1
          counts = earthCounts(1,i)
          outData%counts1(i,scanLinePos) = counts
          outData%array1(i,scanLinePos) = calib1(1) + &
               calib1(2) * counts
          if( 0 .gt. outData%array1(i,scanLinePos) )then
             outData%array1(i,scanLinePos) = 0.
          endif
          outData%array1(i,scanLinePos) &
               = outData%array1(i,scanLinePos) / 100.0 ! percent to fraction

          ! Channel 2
          counts = earthCounts(2,i)
          outData%counts2(i,scanLinePos) = counts
          outData%array2(i,scanLinePos) = calib2(1) + &
               calib2(2) * counts
          if( 0 .gt. outData%array2(i,scanLinePos) )then
             outData%array2(i,scanLinePos) = 0.
          endif
          outData%array2(i,scanLinePos) &
               = outData%array2(i,scanLinePos) / 100.0 ! percent to fraction
          ! Channel 3
          counts = earthCounts(3,i)
          outData%counts3(i,scanLinePos) = counts
          IF( 0. .ne. Coefs%nonlinearCoefs3(1) )THEN
             ! 'new' radiance conversion (NOAA-14)
             radiance = calib3(1) + calib3(2)*counts
             radiance = Coefs%nonlinearCoefs3(1) + &
                  Coefs%nonlinearCoefs3(2)*radiance
          ELSE
             radiance = calib3(1) + calib3(2)*counts
          ENDIF
          if( out_radiances )then
             outData%array3B(i,scanLinePos) = radiance
          else
             ! Use older radiance->temperature conversion
             outData%array3B(i,scanLinePos) = convertBTold( radiance, &
                  Coefs%nuC_numTemps, Coefs%nuC_minT, &
                  Coefs%nuC_maxT, Coefs%nuC_array(:,1) )
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array3B(i,scanLinePos)) )THEN
                outData%array3B(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array3B(i,scanLinePos)) )THEN
                outData%array3B(i,scanLinePos) = NAN_R
             END IF
#endif
          endif
          outData%array3A(i,scanLinePos) = NAN_R
          
          counts = earthCounts(4,i)
          outData%counts4(i,scanLinePos) = counts
          IF( 0. .ne. Coefs%nonlinearCoefs4(1) )THEN
             ! 'new' radiance conversion (NOAA-13/NOAA-14)
             radiance = calib4(1) + calib4(2)*counts
             radiance = Coefs%nonlinearCoefs4(1) + &
                  Coefs%nonlinearCoefs4(2)*radiance + &
                  Coefs%nonlinearCoefs4(3)*radiance*radiance
          ELSE
             radiance = calib4(1) + calib4(2)*counts
          ENDIF
          if( out_radiances )then
             outData%array4(i,scanLinePos) = radiance
          else
             ! Use older radiance->temperature conversion
             outData%array4(i,scanLinePos) = convertBTold( radiance, &
                  Coefs%nuC_numTemps, Coefs%nuC_minT, &
                  Coefs%nuC_maxT, Coefs%nuC_array(:,2) )
             IF( 0. .eq. Coefs%nonlinearCoefs4(1) )THEN
                ! Use lookup table correction
                outData%array4(i,scanLinePos) = &
                     correctBT_Table(4,outData%array4(i,scanLinePos),&
                     prtMean,Coefs)
             ENDIF
#ifdef USE_IEEE
             IF( .NOT.ieee_is_finite(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#elif defined(USE_G95)
             IF( isNAN(outData%array4(i,scanLinePos)) )THEN
                outData%array4(i,scanLinePos) = NAN_R
             END IF
#endif
          endif
          
          ! Check if we have 12 micron data
          IF( Have_12micron )THEN
             counts = earthCounts(5,i)
             outData%counts5(i,scanLinePos) = counts
             IF( 0. .ne. Coefs%nonlinearCoefs5(1) )THEN
                ! 'new' radiance conversion (NOAA-13/NOAA-14)
                radiance = calib5(1) + calib5(2)*counts
                radiance = Coefs%nonlinearCoefs5(1) + &
                     Coefs%nonlinearCoefs5(2)*radiance + &
                     Coefs%nonlinearCoefs5(3)*radiance*radiance
             ELSE
                radiance = calib5(1) + calib5(2)*counts
             ENDIF
             if( out_radiances )then
                outData%array5(i,scanLinePos) = radiance
             else
                ! Use older radiance->temperature conversion
                outData%array5(i,scanLinePos) = convertBTold( radiance, &
                     Coefs%nuC_numTemps, Coefs%nuC_minT, &
                     Coefs%nuC_maxT, Coefs%nuC_array(:,3) )
                IF( 0. .eq. Coefs%nonlinearCoefs5(1) )THEN
                   ! Use lookup table correction
                   outData%array5(i,scanLinePos) = &
                        correctBT_Table(5,outData%array5(i,scanLinePos),&
                        prtMean,Coefs)
                ENDIF
#ifdef USE_IEEE
                IF( .NOT.ieee_is_finite(outData%array5(i,scanLinePos)) )THEN
                   outData%array5(i,scanLinePos) = NAN_R
                END IF
#elif defined(USE_G95)
                IF( isNAN(outData%array5(i,scanLinePos)) )THEN
                   outData%array5(i,scanLinePos) = NAN_R
                END IF
#endif
             endif
          ELSE
             outData%counts5(i,scanLinePos) = NAN_R
             outData%array5(i,scanLinePos) = NAN_R
          ENDIF
       end do

    else

       ! set Earth observations to NAN because we don't know where we are
       outData%counts3(:,scanLinePos) = NAN_R
       outData%counts4(:,scanLinePos) = NAN_R
       outData%counts5(:,scanLinePos) = NAN_R
       outData%array1(:,scanLinePos) = NAN_R
       outData%array2(:,scanLinePos) = NAN_R
       outData%array3A(:,scanLinePos) = NAN_R
       outData%array3B(:,scanLinePos) = NAN_R
       outData%array4(:,scanLinePos) = NAN_R
       outData%array5(:,scanLinePos) = NAN_R

    endif
    outData%dataFilled = .TRUE.

  END SUBROUTINE Parse_RecordV1

  REAL FUNCTION correctBT_Table( Channel, BT, prtMean_in, Coefs )RESULT(BTout)

    INTEGER, INTENT(IN) :: channel
    REAL, INTENT(IN) :: BT
    REAL, INTENT(IN) :: prtMean_in
    TYPE(AVHRR_Instrument_Coefs), INTENT(IN) :: Coefs

    ! Local variables
    INTEGER :: I
    INTEGER :: First, Second
    REAL :: prtMean
    REAL :: BT_diff(NLOOKUP_TABLE_SCENE)
    REAL :: BT_diff2(NLOOKUP_TABLE_SCENE)
    REAL :: Ratio1, Ratio2

#ifdef USE_IEEE
    IF( .NOT.ieee_is_finite(BT) .OR. NAN_R .EQ. prtMean_in )THEN
       BTout = ieee_value(NAN_R,ieee_quiet_nan) ! NaN
       RETURN
    ENDIF
#elif defined(USE_G95)
    IF( isNAN(BT) .OR. NAN_R .EQ. prtMean_in )THEN
       BTout = NAN_R
       RETURN
    ENDIF
#endif

    ! Convert to Centigrade
    prtMean = prtMean_in - 273.15

    ! Get scene temperature dependent values
    IF( 1 .eq. Coefs%number_Target_Temps )THEN
       IF( 4 .eq. Channel )THEN
          BT_Diff(1:Coefs%number_Scene_Temps) = &
               Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,1)
       ELSE IF( 5 .eq. Channel )THEN
          BT_Diff(1:Coefs%number_Scene_Temps) = &
               Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,1)
       ELSE
          CALL out_ERROR('Channel No out of range','correctBT_Table',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
    ELSE
       ! Map BT_Diff by instrument temperature
       IF( prtMean .le. Coefs%Target_Temp(1) )THEN
          IF( 4 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,1)
          ELSE IF( 5 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,1)
          ELSE
             CALL out_ERROR('Channel No out of range','correctBT_Table',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ENDIF
       ELSE IF( prtMean .ge. Coefs%Target_Temp(Coefs%number_Target_Temps) )THEN
          IF( 4 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
                  Coefs%number_Target_Temps)
          ELSE IF( 5 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
                  Coefs%number_Target_Temps)
          ELSE
             CALL out_ERROR('Channel No out of range','correctBT_Table',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ENDIF
       ELSE
          Loop: DO I=1,Coefs%number_Target_Temps-1
             IF( prtMean .ge. Coefs%Target_Temp(I) .and. &
                  prtMean .le. Coefs%Target_Temp(I+1) )THEN
                First = I
                Second = I+1
                Ratio1 = (Coefs%Target_Temp(I+1)-prtMean)/&
                     (Coefs%Target_Temp(I+1)-Coefs%Target_Temp(I))
                Ratio2 = (prtMean - Coefs%Target_Temp(I))/&
                     (Coefs%Target_Temp(I+1)-Coefs%Target_Temp(I))
                EXIT Loop
             ENDIF
          END DO Loop
          IF( 4 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Ratio1*Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
                  First) + &
                  Ratio2*Coefs%Correction_Factor4(1:Coefs%number_Scene_Temps,&
                  Second)
          ELSE IF( 5 .eq. Channel )THEN
             BT_Diff(1:Coefs%number_Scene_Temps) = &
                  Ratio1*Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
                  First) + &
                  Ratio2*Coefs%Correction_Factor5(1:Coefs%number_Scene_Temps,&
                  Second)
          ELSE
             CALL out_ERROR('Channel No out of range','correctBT_Table',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ENDIF

       ENDIF
    ENDIF

    ! Now have in BT_Diff array to interpolate relative to scene temperature
    IF( BT .ge. Coefs%Scene_Temp(1) )THEN
       BTout = BT_Diff(1)
       RETURN
    ELSE IF( BT .le. Coefs%Scene_Temp(Coefs%number_Scene_Temps) )THEN
       BTout = BT_Diff(Coefs%number_Scene_Temps)
       RETURN
    ELSE
       DO I=1,Coefs%number_Scene_Temps-1
          IF( BT .le. Coefs%Scene_Temp(I) .and. &
               BT .ge. Coefs%Scene_Temp(I+1) )THEN
             Ratio1 = (Coefs%Scene_Temp(I+1)-BT)/&
                  (Coefs%Scene_Temp(I+1)-Coefs%Scene_Temp(I))
             Ratio2 = (BT - Coefs%Scene_Temp(I))/&
                  (Coefs%Scene_Temp(I+1)-Coefs%Scene_Temp(I))
             BTout = BT + BT_Diff(I)*Ratio1 + BT_Diff(I+1)*Ratio2
             RETURN
          ENDIF
       END DO
    ENDIF

    if( .not. print_missed_matching )then
       CALL out_WARNING('Missed matching BTs in lookup Table','correctBT_Table',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       print_missed_matching = .TRUE.
    endif
    BTout = NAN_R

  END FUNCTION correctBT_Table

  ! Convert HRPT minor frame to 2byte numbers
  FUNCTION Get_HRPT_Frame(record)

    INTEGER(GbcsInt1), INTENT(IN) :: record(140)
    INTEGER(GbcsInt2) :: Get_HRPT_Frame(103)

    INTEGER :: I,J
    INTEGER(GbcsInt4) :: word
    INTEGER(GbcsInt2) :: twobyte

    J = 0
    DO I=1,140,4
       word = header2int((I-1),record)
       J=J+1
       IF( 103 .ge. J )THEN
          Get_HRPT_Frame(J) = IBITS(word,20,10)
          J=J+1
          IF( 103 .ge. J )THEN
             Get_HRPT_Frame(J) = IBITS(word,10,10)
             IF( 103 .ge. J )THEN
                J=J+1
                Get_HRPT_Frame(J) = IBITS(word,0,10)
             ENDIF
          ENDIF
       ENDIF
    END DO

  END FUNCTION Get_HRPT_Frame

  ! Use extra precision parts for solar angle
  REAL FUNCTION Get_Solar_Angle( Record, pos1, pos2, index )

    INTEGER(GbcsInt1), INTENT(IN) :: Record(:)
    INTEGER, INTENT(IN) :: pos1
    INTEGER, INTENT(IN) :: pos2
    INTEGER, INTENT(IN) :: index

    ! Local variables
    INTEGER :: pos
    INTEGER :: byte_pos
    INTEGER :: start_bit
    INTEGER(GbcsInt2) :: value
    INTEGER(GbcsInt1) :: part_no
    INTEGER(GbcsInt1) :: array(2)
    INTEGER(GbcsInt2) :: twobyte
    REAL :: solar_angle_large

    ! Main part of solar angle
    value = Record(pos1+index-1)
    if( 0 .gt. value )then
       value = 256+value
    endif
    solar_angle_large = value/2.

    IF(validSolZaDecimal)THEN
       ! valid on/after Oct 22 1992
       ! Get decimal part
       ! work out which byte and offset for 3 bit part
       pos = (index-1)*60
       byte_pos = INT(pos/160.)
       ! bit start
       start_bit = pos-byte_pos*8
       IF( 6 .GE. start_bit )THEN
          part_no = IBITS(Record(pos2+byte_pos),start_bit,3)
       ELSE
          ! Cross word boundary - need to take into account Endianess
          IF( Little_Endian )THEN
             array(1) = Record(pos2+byte_pos+1)
             array(2) = Record(pos2+byte_pos)
          ELSE
             array(1) = Record(pos2+byte_pos)
             array(2) = Record(pos2+byte_pos+1)
          ENDIF
          twobyte = TRANSFER(array,twobyte)
          part_no = IBITS(twobyte,start_bit,3)
       ENDIF
    ELSE
       part_no = 0
    ENDIF

    Get_Solar_Angle = solar_angle_large + part_no/10.
    RETURN

  END FUNCTION Get_Solar_Angle

  real function get_satellite_za( longitude, latitude, satellite_long, &
       satellite_lat )

    real, intent(in) :: longitude
    real, intent(in) :: latitude
    real, intent(in) :: satellite_long
    real, intent(in) :: satellite_lat

    double precision, parameter :: M_PI_2 = M_PI/2.0d0
    double precision, parameter :: DEG_TO_RAD = M_PI / 180.d0
    double precision, parameter :: semi_major_axis = 6378.137d0
    double precision, parameter :: f_minus = 289.257222101d0
    double precision, parameter :: e_squared = (2./f_minus) - (1./(f_minus*f_minus))
    DOUBLE PRECISION, PARAMETER :: r = 7220.d0 ! should have different for N11+14 vs N12
    double precision, parameter :: h = 0.

    ! geodetic height of point of interest
    ! height which we'll approximate by 0.

    double precision :: lambda
    double precision :: phi
    double precision :: lambda_s
    double precision :: phi_s

    double precision :: W
    double precision :: value
    double precision :: N

    double precision :: xp, yp, zp
    double precision :: xs, ys, zs

    double precision :: input(3)
    double precision :: output(3)

    double precision :: nu

    lambda = longitude*DEG_TO_RAD
    phi = latitude*DEG_TO_RAD
    lambda_s = satellite_long*DEG_TO_RAD
    phi_s = satellite_lat*DEG_TO_RAD

    value = sin(phi)
    W = sqrt(1 - e_squared * value * value)  ! sqrt( 1 - e^2 * sin^2(lat))

    N = semi_major_axis/W

    value = (N+h)
    xp = value*cos(lambda)*cos(phi)
    yp = value*sin(lambda)*cos(phi)
    zp = (N*(1-e_squared)+h)*sin(phi)

    xs = r*cos(lambda_s)*cos(phi_s)
    ys = r*sin(lambda_s)*cos(phi_s)
    zs = r*sin(phi_s)

    input(1) = xs - xp
    input(2) = ys - yp
    input(3) = zs - zp

    call rotate_matrix( lambda, phi, input, output )

    !alpha = atan2(output(1),output(2))
    nu = atan2(output(3),sqrt(output(1)*output(1)+output(2)*output(2)))

    get_satellite_za = 90. - nu/DEG_TO_RAD
    return

  END function get_satellite_za

  subroutine rotate_matrix( lambda, phi, input, output )

    double precision, intent(in) :: lambda
    double precision, intent(in) :: phi
    double precision, intent(in) :: input(:)
    double precision, intent(out) :: output(:)

    integer :: i, j
    double precision :: R(3,3)

    if( 3 .ne. size(input,DIM=1) )then
       CALL out_ERROR( 'input array not 1x3', &
            'Rotate_Matrix', &
            'NOAA/MCIDAS_LoadImagery.f90', .FALSE., ' ' )
    endif
    if( 3 .ne. size(output,DIM=1) )then
       CALL out_ERROR( 'output array not 1x3', &
            'Rotate_Matrix', &
            'NOAA/MCIDAS_LoadImagery.f90', .FALSE., ' ' )
    endif

    ! Rotation matrix
    R(1,1) = -sin(lambda)
    R(2,1) = cos(lambda)
    R(3,1) = 0.
    R(1,2) = -sin(phi)*cos(lambda)
    R(2,2) = -sin(phi)*sin(lambda)
    R(3,2) = cos(phi)
    R(1,3) = cos(phi)*cos(lambda)
    R(2,3) = cos(phi)*sin(lambda)
    R(3,3) = sin(phI)

    do i=1,3
       output(i) = 0.
       do j=1,3
          output(i) = output(i) + R(j,i)*input(j)
       end do
    end do

  end subroutine rotate_matrix

  ! This can be tightened up later to check for the actual
  ! orbits on messed-up days
  SUBROUTINE checkDate( year, month, day, AVHRR_No, badDay, &
       validSolZaDecimal, validClockDriftInfo, validClockCorrection, validSBBC )

    INTEGER, INTENT(IN) :: year
    INTEGER, INTENT(IN) :: month
    INTEGER, INTENT(IN) :: day
    INTEGER, INTENT(IN) :: AVHRR_No
    LOGICAL, INTENT(OUT) :: badDay
    LOGICAL, INTENT(OUT) :: validSolZaDecimal
    LOGICAL, INTENT(OUT) :: validClockDriftInfo
    LOGICAL, INTENT(OUT) :: validClockCorrection
    LOGICAL, INTENT(OUT) :: validSBBC

    badDay = .FALSE.

    ! Additional decimal portion of 51 Solar Zenith angles was
    ! switched on/off between Sep and Oct 1992, sure to be OK
    ! on/after Oct 22 1992
    IF(year.GT.1992 &
         .OR.(year.EQ.1992.AND.month.GT.10)&
         .OR.(year.EQ.1992.AND.month.EQ.10.AND.day.GT.21))THEN
       validSolZaDecimal = .TRUE.
    ELSE
       validSolZaDecimal = .FALSE.
    END IF

    ! Clock drift delta in milliseconds x 2 + indicator(0 = no time adjustment,
    ! 1 = time adjustment) sure to be valid on/after Nov 16 1994
    IF(year.GT.1994 &
         .OR.(year.EQ.1994.AND.month.GT.11)&
         .OR.(year.EQ.1994.AND.month.EQ.11.AND.day.GT.15))THEN
       validClockDriftInfo = .TRUE.
    ELSE
       validClockDriftInfo = .FALSE.
    END IF

    validClockCorrection = .FALSE.
    ! on Sep 7 1994 an aborted attempt was made to include the clock correction
    ! field, the bad orbits have been already deleted in the NOAA archive itself
    ! corrections off until/on Aug 28 2000
    ! N12, N14 switched on Aug 29 2000
    IF(AVHRR_No.EQ.12.OR.AVHRR_No.EQ.14)THEN
       ! do not use Aug 29 2000 (D00242)
       IF(year.EQ.2000.AND.month.EQ.8.AND.day.EQ.29)THEN
          badDay = .TRUE.
       END IF
       ! corrections on on/after Aug 30 2000
       IF(year.GT.2000 &
            .OR.(year.EQ.2000.AND.month.GT.8)&
            .OR.(year.EQ.2000.AND.month.EQ.8.AND.day.GT.29))THEN
          validClockCorrection = .TRUE.
          ! N12,N14 only were messed up in early 2001
          ! do not use Jan 1,2 2001
          IF(year.EQ.2001.AND.month.EQ.1.AND.day.EQ.1)THEN
             badDay = .TRUE.
          END IF
          IF(year.EQ.2001.AND.month.EQ.1.AND.day.EQ.2)THEN
             badDay = .TRUE.
          END IF
          ! corrections off Jan 3 to Jan 23 2001
          IF(year.EQ.2001.AND.month.EQ.1.AND.day.GT.2.AND.day.LT.24)THEN
             validClockCorrection = .FALSE.
          END IF
          ! do not use Jan 24 2001
          IF(year.EQ.2001.AND.month.EQ.1.AND.day.EQ.24)THEN
             badDay = .TRUE.
          END IF
          ! corrections on on/after Jan 25 2001
       END IF
    END IF
    ! N15,N16 switched on Dec 14 2000
    IF(AVHRR_No.EQ.15.OR.AVHRR_No.EQ.16)THEN
       ! do not use Dec 14 2000
       IF(year.EQ.2000.AND.month.EQ.12.AND.day.EQ.14)THEN
          badDay = .TRUE.
       END IF
       ! corrections on on/after Dec 15 2000
       IF(year.GT.2000 &
            .OR.(year.EQ.2000.AND.month.GT.12)&
            .OR.(year.EQ.2000.AND.month.EQ.12.AND.day.GT.14))THEN
          validClockCorrection = .TRUE.
       END IF
    END IF
    ! N17,N18,N19,M?? all have clock corrections on, either in L1B processing
    ! or the spacescraft clock is correct
    IF(AVHRR_No.GE.17.OR.AVHRR_No.LE.-1)THEN
       validClockCorrection = .TRUE.
    END IF

    ! navigation error of 1 second, do not use the following data
    ! N12: NSS.GHRR.ND.D96088* Mar 28 1996
    IF(AVHRR_No.EQ.12)THEN
       IF(year.EQ.1996.AND.month.EQ.3.AND.day.EQ.28)THEN
          badDay = .TRUE.
       END IF
    END IF

    ! the solar blackbody contamination (SBBC) indicator should not be used
    ! (I think they are not guaranteed to be zero) before Nov 16 1994,
    ! possibly OK on/after Dec 2 1993 (see pod guide historical record appendix E)
    IF(year.GT.1994 &
         .OR.(year.EQ.1994.AND.month.GT.11)&
         .OR.(year.EQ.1994.AND.month.EQ.11.AND.day.GT.15))THEN
       validSBBC = .TRUE.
    ELSE
       validSBBC = .FALSE.
    END IF

    RETURN

  END SUBROUTINE checkDate

  LOGICAL FUNCTION checkTopV1( qualityTop, minor_frame )&
       RESULT(output)

    INTEGER(GbcsInt1), INTENT(IN) :: qualityTop(4)
    INTEGER(GbcsInt2), INTENT(IN) :: minor_frame(103)

    output = .FALSE.

    ! FATAL FLAG - Data should not be used for product generation
    IF( BTEST(qualityTop(1),7) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! DATA JITTER - Resync occurred on this frame (should be OK but mark as bad)
    ! I think this and word 10 bits 7-3 (qualityTop(2)(3:7)) mean that I don't
    ! have to look at the HRPT frame info, specifically word 1 bit 2: 0 = frame
    ! stable; 1 = frame resync occurred
    IF( BTEST(qualityTop(1),4) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! P/N STATUS - Pseudo Noise (P/N) occurred (=1) on the frame, data not used
    ! for calibration computations
    IF( BTEST(qualityTop(1),0) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! BIT SYNC STATUS - Drop lock during frame
    IF( BTEST(qualityTop(2),7) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! SYNC ERROR - Frame Sync word error greater than zero
    IF( BTEST(qualityTop(2),6) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! FRAME SYNC LOCK - Frame Sync previously dropped lock
    ! this should not be a bad frame
    IF( BTEST(qualityTop(2),5) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! FLYWHEELING - Flywheeling detected during this frame
    IF( BTEST(qualityTop(2),4) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! BIT SLIPPAGE - Bit slippage detected during this frame
    IF( BTEST(qualityTop(2),3) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! TIP PARITY - bit 7 is 1st minor frame, bit 3 is 5th
    ! TIP parity should not be important for the AVHRR data but
    ! it does detect gross frame errors.
    IF( IBITS(qualityTop(3),3,5).ne.0 )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! SYNC ERRORS - Number of bit errors in frame sync
    ! this should be flagged in SYNC ERROR above but not always
    ! for really bad frames.
    IF( IBITS(qualityTop(4),2,6).NE.0 )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! frame sync, 60 bits
    IF( 644 .NE. minor_frame(1) .OR. &
         367 .NE. minor_frame(2) .OR. &
         860 .NE. minor_frame(3) .OR. &
         413 .NE. minor_frame(4) .OR. &
         527 .NE. minor_frame(5) .OR. &
         149 .NE. minor_frame(6) )THEN
       ! I hope this never occurs unless it it already marked bad
       ! in qualityTop(2)
       output = .TRUE.
       RETURN
    ENDIF
    ! Frame resync occurred
    ! I hope this never occurs unless it it already marked bad
    ! in qualityTop(1)
    IF( BTEST(minor_frame(7),2) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    RETURN

  END FUNCTION checkTopV1

  LOGICAL FUNCTION checkTopV2_V3( qualityTop, minor_frame, syncErrorCount )&
       RESULT(output)

    INTEGER(GbcsInt4), INTENT(IN) :: qualityTop
    INTEGER(GbcsInt2), INTENT(IN) :: minor_frame(103)
    INTEGER(GbcsInt4), INTENT(IN) :: syncErrorCount

    output = .FALSE.

    ! do not use scan for product generation
    IF( BTEST(qualityTop,31) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! instrument status changed with this scan
    ! probably ok
    IF( BTEST(qualityTop,25) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! sync lock dropped during this frame
    IF( BTEST(qualityTop,24) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! frame sync word error greater than zero
    IF( BTEST(qualityTop,23) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! frame sync previously dropped lock
    ! probably ok
    IF( BTEST(qualityTop,22) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! flywheeling detected during this frame
    ! or frame sync word not valid
    ! what is flywheeling?
    IF( BTEST(qualityTop,21) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! bit slippage detected during this frame
    ! what is bit slippage?
    IF( BTEST(qualityTop,20) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! TIP parity error detected
    ! TIP parity should not be important for the AVHRR data
    ! but this will detect gross errors
    IF( BTEST(qualityTop,8) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! resync occurred on this frame
    ! this is probably OK
    IF( BTEST(qualityTop,1) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! pseudo noise occurred on this frame
    ! could be whole frame, not just AVHRR data
    IF( BTEST(qualityTop,0) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! count of bit errors in frame sync
    ! should be caught by BTEST(qualityTop,23)
    IF( syncErrorCount.GT.0 )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! frame sync, 60 bits
    IF( 644 .NE. minor_frame(1) .OR. &
         367 .NE. minor_frame(2) .OR. &
         860 .NE. minor_frame(3) .OR. &
         413 .NE. minor_frame(4) .OR. &
         527 .NE. minor_frame(5) .OR. &
         149 .NE. minor_frame(6) )THEN
       ! I hope this never occurs unless it it already marked bad
       ! in qualityTop
       output = .TRUE.
       RETURN
    ENDIF
    ! Frame resync occurred
    ! I hope this never occurs unless it it already marked bad
    ! in qualityTop
    IF( BTEST(minor_frame(7),2) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! AVHRR input (0=pseudonoise; 1=normal)
    ! I hope this never occurs unless it it already marked bad
    ! in qualityTop
    IF( .not.BTEST(minor_frame(7),1) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    RETURN

  END FUNCTION checkTopV2_V3

  LOGICAL FUNCTION checkNavigationV1( qualityTop )&
       RESULT(output)

    INTEGER(GbcsInt1), INTENT(IN) :: qualityTop(4)

    output = .FALSE.

    ! TIME ERROR - A time sequence error was detected while processing this
    ! frame
    ! safest to remove these
    IF( BTEST(qualityTop(1),6) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! DATA GAP - A gap precedes this frame, should be ok but
    ! safest to remove these
    IF( BTEST(qualityTop(1),5) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! NO EARTH LOCATION - Earth location data not available
    IF( BTEST(qualityTop(1),2) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    RETURN

  END FUNCTION checkNavigationV1

  LOGICAL FUNCTION checkNavigationV2_V3( qualityTop, qualityTime, qualityEarth, earthLocation, AVHRR_No ) &
       RESULT(output)

    INTEGER(GbcsInt4), INTENT(IN) :: qualityTop
    INTEGER(GbcsInt1), INTENT(IN) :: qualityTime
    INTEGER(GbcsInt1), INTENT(IN) :: qualityEarth
    INTEGER(GbcsInt4), INTENT(IN) :: earthLocation
    INTEGER, INTENT(IN) :: AVHRR_No

    output = .FALSE.

    ! time sequence error detected within this scan
    ! safest to remove these
    IF( BTEST(qualityTop,30) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! data gap precedes this scan (should be ok but be safe)
    IF( BTEST(qualityTop,29) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! earth location data not available
    IF( BTEST(qualityTop,27) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! first good time following a clock update (nominally 0)
    ! probably OK but be safe
    IF( BTEST(qualityTop,26) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    IF( 0 .NE. qualityTime )THEN
       output = .TRUE.
       RETURN
    ENDIF
    IF( 0 .NE. qualityEarth )THEN
       output = .TRUE.
       RETURN
    ENDIF

    IF( AVHRR_No .GE. 1 )THEN ! NOAA
       ! up-to-date earth location
       IF( IBITS( earthLocation, 12, 4 ) .NE. 0 )THEN
          output = .TRUE.
          RETURN
       ENDIF
       ! in YGC or nominal or other mode with good attitude
       IF( IBITS( earthLocation, 8, 4 ) .NE. 0 &
            .AND. IBITS( earthLocation, 8, 4 ) .NE. 1 )THEN
          output = .TRUE.
          RETURN
       ENDIF
    ELSE ! MetOp
       ! manoeuvre not in progress
       IF( BTEST( earthLocation, 18 ) )THEN
          output = .TRUE.
          RETURN
       ENDIF
       ! not sure if the following are used by Metop
       ! up-to-date earth location
       IF( IBITS( earthLocation, 12, 4 ) .NE. 0 )THEN
          output = .TRUE.
          RETURN
       ENDIF
       ! in YGC or nominal or other mode with good attitude
       IF( IBITS( earthLocation, 8, 4 ) .NE. 0 &
            .AND. IBITS( earthLocation, 8, 4 ) .NE. 1 )THEN
          output = .TRUE.
          RETURN
       ENDIF
       ! not sure if other parts of location are useful for Metop
    ENDIF
    ! other parts of location are not useful to test: older NOAAs
    ! sometimes seem to use these as sequence markers (perhaps a test
    ! sequence but there are just too many).

    RETURN

  END FUNCTION checkNavigationV2_V3

  ! an extra test for V3
  LOGICAL FUNCTION checkNavigationV3( qualityTop, qualityTime, qualityEarth, earthLocation, AVHRR_No ) &
       RESULT(output)

    INTEGER(GbcsInt4), INTENT(IN) :: qualityTop
    INTEGER(GbcsInt1), INTENT(IN) :: qualityTime
    INTEGER(GbcsInt1), INTENT(IN) :: qualityEarth
    INTEGER(GbcsInt4), INTENT(IN) :: earthLocation
    INTEGER, INTENT(IN) :: AVHRR_No

    output = .FALSE.

    IF( checkNavigationV2_V3( qualityTop, qualityTime, qualityEarth, &
         earthLocation, AVHRR_No ))THEN
       output = .TRUE.
       RETURN
    ENDIF

    IF( AVHRR_No .GE. 1 )THEN ! NOAA
       ! this is not in V2, and should be covered by the qualityEarth
       ! test in checkNavigationV2_V3?
       ! location within tolerance
       IF( .NOT.BTEST( earthLocation, 17 ) )THEN
          output = .TRUE.
          RETURN
       ENDIF
    ELSE ! MetOp
       ! not sure if the following are used by Metop
       ! this is not in V2, and should be covered by the qualityEarth
       ! test in checkNavigationV2_V3?
       ! location within tolerance
       IF( .NOT.BTEST( earthLocation, 17 ) )THEN
          output = .TRUE.
          RETURN
       ENDIF
    ENDIF

    RETURN

  END FUNCTION checkNavigationV3

  LOGICAL FUNCTION checkTimeV1( qualityTop )&
       RESULT(output)

    INTEGER(GbcsInt1), INTENT(IN) :: qualityTop(4)

    output = .FALSE.

    ! TIME ERROR - A time sequence error was detected while processing this frame
    IF( BTEST(qualityTop(1),6) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! DATA GAP - A gap precedes this frame, should be ok but
    ! safest to remove these
    IF( BTEST(qualityTop(1),5) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    RETURN

  END FUNCTION checkTimeV1

  LOGICAL FUNCTION checkTimeV2_V3( qualityTop, qualityTime ) &
       RESULT(output)

    INTEGER(GbcsInt4), INTENT(IN) :: qualityTop
    INTEGER(GbcsInt1), INTENT(IN) :: qualityTime

    output = .FALSE.

    ! time sequence error detected within this scan
    IF( BTEST(qualityTop,30) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! data gap precedes this scan (should be ok but be safe)
    IF( BTEST(qualityTop,29) )THEN
       output = .TRUE.
       RETURN
    ENDIF
    ! first good time following a clock update (nominally 0)
    ! probably OK but be safe
    IF( BTEST(qualityTop,26) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    IF( 0 .NE. qualityTime )THEN
       output = .TRUE.
       RETURN
    ENDIF

    RETURN

  END FUNCTION checkTimeV2_V3

  LOGICAL FUNCTION checkCalibrationV1( qualityTop, solar_3b, solar_4, solar_5 )&
       RESULT(output)

    INTEGER(GbcsInt1), INTENT(IN) :: qualityTop(4)
    LOGICAL, INTENT(OUT) :: solar_3B
    LOGICAL, INTENT(OUT) :: solar_4
    LOGICAL, INTENT(OUT) :: solar_5

    INTEGER :: intval

    output = .FALSE.
    solar_3b = .FALSE.
    solar_4 = .FALSE.
    solar_5= .FALSE.

    ! CALIBRATION - Insufficient data for calibration
    IF( BTEST(qualityTop(1),3) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    IF(validSBBC)THEN
       ! these are only valid on/after Nov 16 1994
       IF( BTEST(qualityTop(2),2) )THEN
          solar_3b = .TRUE.
       ELSE
          solar_3b = .FALSE.
       ENDIF
       IF( BTEST(qualityTop(2),1) )THEN
          solar_4 = .TRUE.
       ELSE
          solar_4 = .FALSE.
       ENDIF
       IF( BTEST(qualityTop(2),0) )THEN
          solar_5 = .TRUE.
       ELSE
          solar_5 = .FALSE.
       ENDIF
    ENDIF

    RETURN

  END FUNCTION checkCalibrationV1

  LOGICAL FUNCTION checkCalibrationV2_V3(qualityTop, qualityCal, solar_3b, solar_4, solar_5 )&
       RESULT(output)

    INTEGER(GbcsInt4), INTENT(IN) :: qualityTop
    INTEGER(GbcsInt1), INTENT(IN) :: qualityCal
    LOGICAL, INTENT(OUT) :: solar_3B
    LOGICAL, INTENT(OUT) :: solar_4
    LOGICAL, INTENT(OUT) :: solar_5

    INTEGER :: intval

    output = .FALSE.
    solar_3b = .FALSE.
    solar_4 = .FALSE.
    solar_5= .FALSE.

    ! insufficient data for calibration
    IF( BTEST(qualityTop,28) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! scan line not calibrated: all IR channels failed calibration
    ! or because of bad time
    IF( BTEST(qualityCal,7) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! No visible calibration due to either the presence of MIRP psendonoise
    ! in place of AVHRR data (NOAA only) or calibration processing turned off
    ! not tested, mainly interested in IR channels

    ! scan line was not calibrated because of satellite maneuver (MetOp only)
    ! NOAA has zero fill
    IF( BTEST(qualityCal,0) )THEN
       output = .TRUE.
       RETURN
    ENDIF

    ! test for PRT and individual channels in parse routine

    ! Check solar contamination
    intval = IBITS(qualityTop,6,2)
    IF( 0 .ne. intval )THEN
       solar_3b = .TRUE.
    ELSE
       solar_3b = .FALSE.
    ENDIF
    intval = IBITS(qualityTop,4,2)
    IF( 0 .ne. intval )THEN
       solar_4 = .TRUE.
    ELSE
       solar_4 = .FALSE.
    ENDIF
    intval = IBITS(qualityTop,2,2)
    IF( 0 .ne. intval )THEN
       solar_5 = .TRUE.
    ELSE
       solar_5 = .FALSE.
    ENDIF

    RETURN

  END FUNCTION checkCalibrationV2_V3

  ! Set bad data
  SUBROUTINE set_outdata_bad( outData, scanLinePos )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: scanLinePos

    CALL set_time_bad( outData, scanLinePos )
    CALL set_nav_bad( outData, scanLinePos )
    CALL set_calib_bad( outData, scanLinePos )

  END SUBROUTINE set_outdata_bad

  SUBROUTINE set_time_bad( outData, scanLinePos )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: scanLinePos

    outData%year(scanLinePos) = NAN_I
    outData%month(scanLinePos) = NAN_I
    outData%day(scanLinePos) = NAN_I
    outData%dayNo(scanLinePos) = NAN_I
    outData%hours(scanLinePos) = NAN_R
    outData%UTC_msecs(scanLinePos) = NAN_I
    outData%time(scanLinePos) = NAN_D

  END SUBROUTINE set_time_bad

  SUBROUTINE set_nav_bad( outData, scanLinePos )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: scanLinePos

    ! don't set outData%scanLineNumber(scanLinePos) to NAN_I
    outData%Lon(:,scanLinePos) = NAN_R
    outData%Lat(:,scanLinePos) = NAN_R
    outData%satZA(:,scanLinePos) = NAN_R
    outData%solZA(:,scanLinePos) = NAN_R
    outData%relAz(:,scanLinePos) = NAN_R

  END SUBROUTINE set_nav_bad

  SUBROUTINE set_calib_bad( outData, scanLinePos )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: scanLinePos

    outData%Counts1(:,scanLinePos) = NAN_R
    outData%Counts2(:,scanLinePos) = NAN_R
    outData%Counts3(:,scanLinePos) = NAN_R
    outData%Counts4(:,scanLinePos) = NAN_R
    outData%Counts5(:,scanLinePos) = NAN_R
    outData%array1(:,scanLinePos) = NAN_R
    outData%array2(:,scanLinePos) = NAN_R
    outData%array3A(:,scanLinePos) = NAN_R
    outData%array3B(:,scanLinePos) = NAN_R
    outData%array4(:,scanLinePos) = NAN_R
    outData%array5(:,scanLinePos) = NAN_R
    outData%prt1(scanLinePos) = NAN_R
    outData%prt2(scanLinePos) = NAN_R
    outData%prt3(scanLinePos) = NAN_R
    outData%prt4(scanLinePos) = NAN_R
    outData%prt1Counts(scanLinePos) = NAN_R
    outData%prt2Counts(scanLinePos) = NAN_R
    outData%prt3Counts(scanLinePos) = NAN_R
    outData%prt4Counts(scanLinePos) = NAN_R
    outData%prt1CountsAll(:,scanLinePos) = NAN_R
    outData%prt2CountsAll(:,scanLinePos) = NAN_R
    outData%prt3CountsAll(:,scanLinePos) = NAN_R
    outData%prt4CountsAll(:,scanLinePos) = NAN_R
    outData%bb3(scanLinePos) = NAN_R
    outData%bb4(scanLinePos) = NAN_R
    outData%bb5(scanLinePos) = NAN_R
    outData%sp3(scanLinePos) = NAN_R
    outData%sp4(scanLinePos) = NAN_R
    outData%sp5(scanLinePos) = NAN_R
    outData%bbodyFilter3(:,scanLinePos) = NAN_I
    outData%bbodyFilter4(:,scanLinePos) = NAN_I
    outData%bbodyFilter5(:,scanLinePos) = NAN_I
    outData%spaceFilter1(:,scanLinePos) = NAN_I
    outData%spaceFilter2(:,scanLinePos) = NAN_I
    outData%spaceFilter3a(:,scanLinePos) = NAN_I
    outData%spaceFilter3(:,scanLinePos) = NAN_I
    outData%spaceFilter4(:,scanLinePos) = NAN_I
    outData%spaceFilter5(:,scanLinePos) = NAN_I
    outData%patch(scanLinePos) = NAN_R
    outData%patchExtended(scanLinePos) = NAN_R
    outData%Radiator(scanLinePos) = NAN_R
    outData%Cooler(scanLinePos) = NAN_R
    outData%a_d_conv(scanLinePos) = NAN_R
    outData%motor(scanLinePos) = NAN_R
    outData%electronics(scanLinePos) = NAN_R
    outData%baseplate(scanLinePos) = NAN_R
    outData%calib1(:,scanLinePos) = NAN_R
    outData%calib1_2(:,scanLinePos) = NAN_R
    outData%calib1_intercept(scanLinePos) = NAN_R
    outData%calib2(:,scanLinePos) = NAN_R
    outData%calib2_2(:,scanLinePos) = NAN_R
    outData%calib2_intercept(scanLinePos) = NAN_R
    outData%calib3A(:,scanLinePos) = NAN_R
    outData%calib3A_2(:,scanLinePos) = NAN_R
    outData%calib3A_intercept(scanLinePos) = NAN_R
    outData%calib3(:,scanLinePos) = NAN_R
    outData%calib4(:,scanLinePos) = NAN_R
    outData%calib5(:,scanLinePos) = NAN_R
    outData%orig_solar_contamination_3B(scanLinePos) = .FALSE.
    outData%orig_solar_contamination_4(scanLinePos) = .FALSE.
    outData%orig_solar_contamination_5(scanLinePos) = .FALSE.
    outData%clavr_mask(:,scanLinePos) = 3 ! default cloudy
    outData%clavrx_mask(:,scanLinePos) = 7 ! default unprocessed and cloudy
    outData%clavrx_prb(:,scanLinePos) = -128 ! default unprocessed
    IF( outData%walton_there )THEN
       outData%array3B_error(:,scanLinePos) = NAN_R
       outData%array4_error(:,scanLinePos) = NAN_R
       outData%array5_error(:,scanLinePos) = NAN_R
    ENDIF
    IF( outData%newCalibration_There )THEN
       outData%new_array1(:,scanLinePos) = NAN_R
       outData%new_array2(:,scanLinePos) = NAN_R
       outData%new_array3A(:,scanLinePos) = NAN_R
       outData%new_array3B(:,scanLinePos) = NAN_R
       outData%new_array4(:,scanLinePos) = NAN_R
       outData%new_array5(:,scanLinePos) = NAN_R
       outData%new_array1_error(:,scanLinePos) = NAN_R
       outData%new_array2_error(:,scanLinePos) = NAN_R
       outData%new_array3A_error(:,scanLinePos) = NAN_R
       outData%new_array3B_error(:,scanLinePos) = NAN_R
       outData%new_array4_error(:,scanLinePos) = NAN_R
       outData%new_array5_error(:,scanLinePos) = NAN_R
       outData%solar_contamination_3B(scanLinePos) = .FALSE.
       outData%solar_contamination_4(scanLinePos) = .FALSE.
       outData%solar_contamination_5(scanLinePos) = .FALSE.
       outData%new_calib3(:,scanLinePos) = NAN_R
       outData%new_calib4(:,scanLinePos) = NAN_R
       outData%new_calib5(:,scanLinePos) = NAN_R
    ENDIF
  END SUBROUTINE set_calib_bad

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
    WRITE(*,'(''Start Time:'',i4.4,''/'',i2.2,''/'',i2.2,'':'',i2.2,'':'',i2.2,'':'',i2.2)')&
         start_time%year,start_time%month,start_time%day,start_time%hour,&
         start_time%minute,start_time%seconds
    end_time   = AVHRR_Date(header2ushort(96,header), header2ushort(98,header), header2int(100,header))
    WRITE(*,'(''End Time:'',i4.4,''/'',i2.2,''/'',i2.2,'':'',i2.2,'':'',i2.2,'':'',i2.2)')&
         end_time%year,end_time%month,end_time%day,end_time%hour,&
         end_time%minute,end_time%seconds

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

  ! Read data header
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
    WRITE(*,'(''Start Time:'',i4.4,''/'',i2.2,''/'',i2.2,'':'',i2.2,'':'',i2.2,'':'',i2.2)')&
         start_time%year,start_time%month,start_time%day,start_time%hour,&
         start_time%minute,start_time%seconds
    end_time   = AVHRR_Date(header2ushort(96,header), header2ushort(98,header), header2int(100,header))
    WRITE(*,'(''End Time:'',i4.4,''/'',i2.2,''/'',i2.2,'':'',i2.2,'':'',i2.2,'':'',i2.2)')&
         end_time%year,end_time%month,end_time%day,end_time%hour,&
         end_time%minute,end_time%seconds

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

  ! Read data header
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
    WRITE(*,'(''Start Time:'',i4.4,''/'',i3.3,'':'',f6.3)')year,dayno,hours
    CALL Parse_Date_Time( header, 11, year, dayno, hours, UTC_msecs )
    end_time = AVHRR_Date(year, dayno, UTC_msecs)
    WRITE(*,'(''End Time:'',i4.4,''/'',i3.3,'':'',f6.3)')year,dayno,hours


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

  ! Routine to check if we have an ARVHIVE header and which version
  ! of file we have
  SUBROUTINE Check_Header( Unit, VersionNo, sizeOfRecords )

#ifdef USE_GZFILE
    INTEGER(c_int), INTENT(IN) :: Unit
#else
    INTEGER, INTENT(IN) :: Unit
#endif
    INTEGER, INTENT(OUT) :: VersionNo
    INTEGER, INTENT(OUT) :: sizeOfRecords

    ! Header size in bytes
    INTEGER(GbcsInt1) :: Header(ARCHIVE_HEADER_SIZE)
    INTEGER(GbcsInt1) :: Header_POD(ARCHIVE_HEADER_SIZE_POD)
    CHARACTER(LEN=2) :: tempString2
    CHARACTER(LEN=3) :: tempString3
    CHARACTER(LEN=4) :: tempString4
    CHARACTER(LEN=6) :: tempString6
    CHARACTER(LEN=15) :: tempString15
    CHARACTER(LEN=3) :: all1
    CHARACTER(LEN=3) :: all2
    CHARACTER(LEN=3) :: all3
    CHARACTER(LEN=3) :: all4
    INTEGER :: version
    INTEGER :: Ios
    INTEGER :: tenbit
    INTEGER :: type
    INTEGER :: year
    INTEGER :: dayno
    REAL :: hours
    INTEGER :: UTC_msecs
#ifndef USE_GZFILE
    INTEGER :: STAT
#endif

#ifdef USE_GZFILE
    Ios = gzread(Unit,Header)
    IF(Ios.LE.0)THEN
       CALL Check_IOS(1,'Reading ARCHIVE Header','Check_Header',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    END IF
#else
    READ(Unit,IOSTAT=STAT)Header
    CALL Check_IOS(STAT,'Reading ARCHIVE Header','Check_Header',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

    !
    ! Look for 'ALL' at 76, 79, 83, 87 (beg/end Latitude/Longitude)
    ! All means it's an archive heasder - this is true of both KLM (AVHRR/3) and POD (AVHRR/2s)
    !
    CALL Convert_Int1_String(75,3,Header,all1)
    CALL Convert_Int1_String(78,3,Header,all2)
    CALL Convert_Int1_String(81,3,Header,all3)
    CALL Convert_Int1_String(86,3,Header,all4)
    IF( 'ALL' .NE. all1 .and. 'ALL' .NE. all2 .and. &
         'ALL' .NE. all3 .and. 'ALL' .NE. all4 )THEN
!    IF( ICHAR('A') .NE. Header(76) .AND. &
!         ICHAR('L') .NE. Header(77) .AND. &
!         ICHAR('L') .NE. Header(78) .AND. &
!         ICHAR('A') .NE. Header(79) .AND. &
!         ICHAR('L') .NE. Header(80) .AND. &
!         ICHAR('L') .NE. Header(81) .AND. &
!         ICHAR('A') .NE. Header(82) .AND. &
!         ICHAR('L') .NE. Header(83) .AND. &
!         ICHAR('L') .NE. Header(84) .AND. &
!         ICHAR('A') .NE. Header(86) .AND. &
!         ICHAR('L') .NE. Header(87) .AND. &
!         ICHAR('L') .NE. Header(88) )THEN
       !
       ! Not from the CLASS Archive/Tape header
       !
       CALL Convert_Int1_String(0,3,Header,tempString3)
       IF(  'CMS' .eq. tempString3 .or. &
            'DSS' .eq. tempString3 .or. &
            'NSS' .eq. tempString3 .or. &
            'UKM' .eq. tempString3 )THEN
          ! Get version number
          version = Header2ushort( 4, Header )
          SELECT CASE(version)
          CASE(1)
             ! Some Noaa-15 data is listed as version 1 which shouldn't exist.
             VersionNo = 2
          CASE(2)
             VersionNo = 2
          CASE(3)
             VersionNo = 3
          CASE(4)
             VersionNo = 3
          CASE(5)
             VersionNo = 5
          CASE DEFAULT
             WRITE(*,*)' Version = ',Version
             call out_ERROR('Version Number out of range',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END SELECT
          type = Header2ushort( 76 , Header )
          SELECT CASE(type)
          case(2)
             sizeOfRecords = 4608
          case(1,3,13)
             sizeOfRecords = 15872
          CASE DEFAULT
             call out_ERROR('Data types out of range',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END SELECT
       ELSE
          ! Either Version 0 (TIROS-N) or 1 data - POD
          ! Have to check date to distinguish between TIROS-N and NOAA-11
          ! Check date
          IF( 1 .eq. Header(1) )THEN
             CALL Parse_Date_Time( Header, 3, year, dayno, hours, UTC_msecs )
             IF( 1990 .lt. year )THEN
                VersionNo = 1
             ELSE
                VersionNo = 0
             ENDIF
          ELSE
             VersionNo = 1
          ENDIF
          ! Get data type
          type = Parse_Data_Type( Header(2) )
          SELECT CASE(type)
          case(2)
             sizeOfRecords = 3220
          case(1,3)
             sizeOfRecords = 7400*2  !
          CASE DEFAULT
             call out_ERROR('Data types out of range',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END SELECT
       ENDIF

       ! Rewind file
#ifdef USE_GZFILE
       Ios = gzrewind(Unit)
#else
       rewind(Unit)
#endif

    ELSE

       !
       ! Check for NSS.GHRR in filename - from byte 31
       ! If there this is KLM archive header not POD (AVHRR/2) Tape header
       IF( ICHAR('N') .eq. Header(31) .and. &
            ICHAR('S') .eq. Header(32) .and. &
            ICHAR('S') .eq. Header(33) .and. &
            ICHAR('.') .eq. Header(34) .and. &
            ICHAR('G') .eq. Header(35) .and. &
            ICHAR('H') .eq. Header(36) .and. &
            ICHAR('R') .eq. Header(37) .and. &
            ICHAR('R') .eq. Header(38) )THEN

          ! Class archive header - can get version no and size from here
          CALL Convert_Int1_String(117,2,Header,tempString2)
          READ(tempString2,*,IOSTAT=Ios)tenbit
          call check_IOS(Ios,'Cannot parse format bit no',&
               'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          if( 10 .ne. tenbit )then
             call out_ERROR('Data file not in 10 bit format',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          endif

          CALL Convert_Int1_String(161,15,Header,tempString15)
          if( 'NOAA Level 1b v' .ne. tempString15 )then
             call out_ERROR('Not NOAA Level 1b data',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          endif

          CALL Convert_Int1_String(176,4,Header,tempString4)
          READ(tempString4,*,IOSTAT=Ios)version
          call check_IOS(Ios,'Cannot parse version number',&
               'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

          SELECT CASE(version)
          CASE(2)
             VersionNo = 2
          CASE(3)
             VersionNo = 3
          CASE(4)
             VersionNo = 3
          CASE(5)
             VersionNo = 5
          CASE DEFAULT
             call out_ERROR('Version Number out of range',&
                  'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END SELECT

          ! Get size of record
          CALL Convert_Int1_String(181,6,Header,tempString6)
          READ(tempString6,*,IOSTAT=Ios)sizeOfRecords
          call check_IOS(Ios,'Cannot parse size of records',&
               'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

       ELSE

          ! Rewind file
#ifdef USE_GZFILE
          Ios = gzrewind(Unit)
#else
          rewind(Unit)
#endif

          !
          ! Read in 122 byte POD header
          !
#ifdef USE_GZFILE
          Ios = gzread(Unit,Header_POD)
          IF(Ios.LE.0)THEN
             CALL Check_IOS(1,'Reading POD ARCHIVE Header','Check_Header',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END IF
#else
          READ(Unit,IOSTAT=STAT)Header_POD
          CALL Check_IOS(STAT,'Reading POD ARCHIVE Header','Check_Header',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

#ifdef USE_GZFILE
          Ios = gzread(Unit,Header)
          IF(Ios.LE.0)THEN
             CALL Check_IOS(1,'Reading POD ARCHIVE Header 2','Check_Header',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END IF
#else
          READ(Unit,IOSTAT=STAT)Header
          CALL Check_IOS(STAT,'Reading POD ARCHIVE Header 2','Check_Header',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

          CALL Convert_Int1_String(0,3,Header,tempString3)
          IF(  'CMS' .eq. tempString3 .or. &
               'DSS' .eq. tempString3 .or. &
               'NSS' .eq. tempString3 .or. &
               'UKM' .eq. tempString3 )THEN
             ! Get version number
             version = Header2ushort( 4, Header )
             SELECT CASE(version)
             CASE(1)
                ! Some Noaa-15 data is listed as version 1 which shouldn't exist.
                VersionNo = 2
             CASE(2)
                VersionNo = 2
             CASE(3)
                VersionNo = 3
             CASE(4)
                VersionNo = 3
             CASE(5)
                VersionNo = 5
             CASE DEFAULT
                WRITE(*,*)' Version = ',Version
                call out_ERROR('Version Number out of range',&
                     'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
             END SELECT
             type = Header2ushort( 76 , Header )
             SELECT CASE(type)
             case(2)
                sizeOfRecords = 4608
             case(1,3,13)
                sizeOfRecords = 15872
             CASE DEFAULT
                call out_ERROR('Data types out of range',&
                     'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
             END SELECT
          ELSE
             ! Either Version 0 (TIROS-N) or 1 data - POD
             ! Have to check date to distinguish between TIROS-N and NOAA-11
             ! Check date
             IF( 1 .eq. Header(1) )THEN
                CALL Parse_Date_Time( Header, 3, year, dayno, hours, UTC_msecs )
                IF( 1990 .lt. year )THEN
                   VersionNo = 1
                ELSE
                   VersionNo = 0
                ENDIF
             ELSE
                VersionNo = 1
             ENDIF
             ! Get data type
             type = Parse_Data_Type( Header(2) )
             SELECT CASE(type)
             case(2)
                sizeOfRecords = 3220
             case(1,3)
                sizeOfRecords = 7400*2  !
             CASE DEFAULT
                call out_ERROR('Data types out of range',&
                     'Check_Header','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
             END SELECT
          ENDIF

          ! Rewind file
#ifdef USE_GZFILE
          Ios = gzrewind(Unit)
#else
          rewind(Unit)
#endif
          !
          ! Strip 122 byte POD header
          !
#ifdef USE_GZFILE
          Ios = gzread(Unit,Header_POD)
          IF(Ios.LE.0)THEN
             CALL Check_IOS(1,'Reading POD ARCHIVE Header','Check_Header',&
                  'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          END IF
#else
          READ(Unit,IOSTAT=STAT)Header_POD
          CALL Check_IOS(STAT,'Reading POD ARCHIVE Header','Check_Header',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
#endif

       ENDIF
    ENDIF

  END SUBROUTINE Check_Header

#ifdef USE_HDF5
  SUBROUTINE Read_CLAVR_Nav(navPathname, nScans)

    CHARACTER(LEN=512), INTENT(IN) :: navPathname
    INTEGER, INTENT(IN) :: nscans

    INTEGER :: I,J
    INTEGER :: pos
    INTEGER :: stat
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: dspace_id
    INTEGER(HID_T) :: attr_id
    INTEGER(HSIZE_T), DIMENSION(2) :: dims
    INTEGER(HSIZE_T), DIMENSION(2) :: maxdims
    INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
    INTEGER :: rank
    INTEGER :: hdferr
    REAL :: missingValue
    LOGICAL :: exists
    INTEGER :: line
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_clavrLonAnchor, &
         temp_clavrLatAnchor

    attr_dims = 0

    ! use the CLAVR-x timing corrections if not done in L1B
    ! GAC only
    INQUIRE(file=TRIM(navPathname), exist=exists)
    IF (.NOT. exists) THEN
       ! no CLAVR-x file, L1B is as good as we get
       CALL out_WARNING('No CLAVR-x nav file',&
            'Read_CLAVR_Nav','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       useL1BClock = .TRUE.
    ELSE
       CALL h5open_f(hdferr)
       CALL h5fopen_f(navPathname,H5F_ACC_RDONLY_F,file_id,hdferr)

       ! assume lon and lat have same dimensions
       CALL h5dopen_f(file_id,'/pixel_longitude',dset_id,hdferr)
       CALL h5dget_space_f(dset_id,dspace_id,hdferr)
       CALL h5sget_simple_extent_ndims_f(dspace_id,rank,hdferr)
       IF (rank.NE.2) CALL out_ERROR('CLAVR-x input not 2D',&
            'Read_CLAVR_Nav','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       CALL h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,hdferr)
       CALL h5sclose_f(dspace_id,hdferr)
       ! nScans along track: dims(2)
       IF (nScans.NE.dims(2)) CALL out_ERROR('CLAVR-x does not have same number of scans as L1B',&
            'Read_CLAVR_Nav','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ! should be 51 across track : dims(1)
       IF (51.NE.dims(1)) THEN
          ALLOCATE(temp_clavrLonAnchor(dims(1),dims(2)),&
               clavrLonAnchor(51,dims(2)),stat=stat)
          CALL Check_IOS(stat,'Allocating clavrLonAnchor','Read_CLAVR_Nav',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ALLOCATE(temp_clavrLatAnchor(dims(1),dims(2)),&
               clavrLatAnchor(51,dims(2)),stat=stat)
          CALL Check_IOS(stat,'Allocating clavrLatAnchor','Read_CLAVR_Nav',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

          ! the scaling of the CLAVR-x output has been set to no scaling,
          ! the values are reals

          ! missing value (-999 usually in CLAVR-x)
          CALL h5aopen_f(dset_id,'RANGE_MISSING',attr_id,hdferr)
          CALL h5aread_f(attr_id,H5T_NATIVE_REAL,missingValue,attr_dims,hdferr)
          CALL h5aclose_f(attr_id,hdferr)

          CALL h5dread_f(dset_id,H5T_NATIVE_REAL,temp_clavrLonAnchor,dims,hdferr)
          CALL h5dclose_f(dset_id,hdferr)

          !
          ! Get tie points
          !
          ! Note have to scale longitude by 0.005493332 with missing values at
          ! -32768
          !
          if( dims(1) .eq. 2048 )then
             DO J=1,dims(2)
                DO I=1,51
                   POS = 25 + (I-1)*40
                   IF( temp_clavrLonAnchor(POS,J) .ne. -32768 )THEN
                      clavrLonAnchor(I,J) = temp_clavrLonAnchor(POS,J)*0.005493332
                   ELSE
                      clavrLonAnchor(I,J) = NAN_R
                   ENDIF
                END DO
             END DO
          else if( dims(1) .eq. 409 )then
             DO J=1,dims(2)
                DO I=1,51
                   POS = 5 + (I-1)*8
                   IF( temp_clavrLonAnchor(POS,J) .ne. -32768 )THEN
                      clavrLonAnchor(I,J) = temp_clavrLonAnchor(POS,J)*0.005493332
                   ELSE
                      clavrLonAnchor(I,J) = NAN_R
                   ENDIF
                END DO
             END DO
          else
             CALL out_ERROR('CLAVR-x does not have 409 or 2048 elements',&
                  'Read_CLAVR_Nav','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          endif
          !
          ! Deallocate temp arrays
          !
          DEALLOCATE(temp_clavrLonAnchor)


          ! bad scanlines have the lon and lat set to missing value
          WHERE(clavrLonAnchor.EQ.missingValue) clavrLonAnchor = NAN_R
          ! clavr-x also outputs NaNs, and more than one latitude == 0 in a scan line
          ! (this seems to be only at the end of the orbit and very occasionally),
          ! and also more than one latitude == 90, all meaning reposnx has failed,
          ! probably because the frame contains rubbish.
          ! Don't test for more than one longitude == 0 since that can be valid for
          ! near the pole and the scanline along the Greenwich meridian
#ifdef USE_IEEE
          WHERE(ieee_is_nan(clavrLonAnchor)) clavrLonAnchor = NAN_R
#elif defined(USE_G95)
          WHERE(isNAN(clavrLonAnchor)) clavrLonAnchor = NAN_R
#endif

          ! shift to correct range for reader is done in parse record V1 or V2

          CALL h5dopen_f(file_id,'/pixel_latitude',dset_id,hdferr)

          ! missing value (-999 usually in CLAVR-x)
          CALL h5aopen_f(dset_id,'RANGE_MISSING',attr_id,hdferr)
          CALL h5aread_f(attr_id,H5T_NATIVE_REAL,missingValue,attr_dims,hdferr)
          CALL h5aclose_f(attr_id,hdferr)

          CALL h5dread_f(dset_id,H5T_NATIVE_REAL,temp_clavrLatAnchor,dims,hdferr)
          CALL h5dclose_f(dset_id,hdferr)

          !
          ! Get tie points
          !
          if( dims(1) .eq. 2048 )then
             DO J=1,dims(2)
                DO I=1,51
                   POS = 25 + (I-1)*40
                   IF( temp_clavrLatAnchor(POS,J) .ne. -32768 )THEN
                      clavrLatAnchor(I,J) = temp_clavrLatAnchor(POS,J)*0.002746666
                   ELSE
                      clavrLatAnchor(I,J) = NAN_R
                   ENDIF
                END DO
             END DO
          else if( dims(1) .eq. 409 )then
             DO J=1,dims(2)
                DO I=1,51
                   POS = 5 + (I-1)*8
                   IF( temp_clavrLatAnchor(POS,J) .ne. -32768 )THEN
                      clavrLatAnchor(I,J) = temp_clavrLatAnchor(POS,J)*0.002746666
                   ELSE
                      clavrLatAnchor(I,J) = NAN_R
                   ENDIF
                END DO
             END DO
          else
             CALL out_ERROR('CLAVR-x does not have 409 or 2048 elements',&
                  'Read_CLAVR_Nav','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          endif
          !
          ! Deallocate temp arrays
          !
          DEALLOCATE(temp_clavrLatAnchor)

          ! bad scanlines have the lon and lat set to missing value
          WHERE(clavrLatAnchor.EQ.missingValue) clavrLatAnchor = NAN_R
          ! clavr-x also outputs NaNs, and more than one latitude == 0 in a scan line
          ! (this seems to be only at the end of the orbit and very occasionally),
          ! and also more than one latitude == 90, all meaning reposnx has failed,
          ! probably because the frame contains rubbish.
          ! Don't test for more than one longitude == 0 since that can be valid for
          ! near the pole and the scanline along the Greenwich meridian
#ifdef USE_IEEE
          WHERE(ieee_is_nan(clavrLatAnchor)) clavrLatAnchor = NAN_R
#elif defined(USE_G95)
          WHERE(isNAN(clavrLatAnchor)) clavrLatAnchor = NAN_R
#endif
       ELSE

          ALLOCATE(clavrLonAnchor(dims(1),dims(2)),stat=stat)
          CALL Check_IOS(stat,'Allocating clavrLonAnchor','Read_CLAVR_Nav',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
          ALLOCATE(clavrLatAnchor(dims(1),dims(2)),stat=stat)
          CALL Check_IOS(stat,'Allocating clavrLatAnchor','Read_CLAVR_Nav',&
               'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

          ! the scaling of the CLAVR-x output has been set to no scaling,
          ! the values are reals

          ! missing value (-999 usually in CLAVR-x)
          CALL h5aopen_f(dset_id,'RANGE_MISSING',attr_id,hdferr)
          CALL h5aread_f(attr_id,H5T_NATIVE_REAL,missingValue,attr_dims,hdferr)
          CALL h5aclose_f(attr_id,hdferr)

          CALL h5dread_f(dset_id,H5T_NATIVE_REAL,clavrLonAnchor,dims,hdferr)
          CALL h5dclose_f(dset_id,hdferr)

          ! bad scanlines have the lon and lat set to missing value
          WHERE(clavrLonAnchor.EQ.missingValue) clavrLonAnchor = NAN_R
          ! clavr-x also outputs NaNs, and more than one latitude == 0 in a scan line
          ! (this seems to be only at the end of the orbit and very occasionally),
          ! and also more than one latitude == 90, all meaning reposnx has failed,
          ! probably because the frame contains rubbish.
          ! Don't test for more than one longitude == 0 since that can be valid for
          ! near the pole and the scanline along the Greenwich meridian
#ifdef USE_IEEE
          WHERE(ieee_is_nan(clavrLonAnchor)) clavrLonAnchor = NAN_R
#elif defined(USE_G95)
          WHERE(isNAN(clavrLonAnchor)) clavrLonAnchor = NAN_R
#endif

          ! shift to correct range for reader is done in parse record V1 or V2

          CALL h5dopen_f(file_id,'/pixel_latitude',dset_id,hdferr)

          ! missing value (-999 usually in CLAVR-x)
          CALL h5aopen_f(dset_id,'RANGE_MISSING',attr_id,hdferr)
          CALL h5aread_f(attr_id,H5T_NATIVE_REAL,missingValue,attr_dims,hdferr)
          CALL h5aclose_f(attr_id,hdferr)

          CALL h5dread_f(dset_id,H5T_NATIVE_REAL,clavrLatAnchor,dims,hdferr)
          CALL h5dclose_f(dset_id,hdferr)

          ! bad scanlines have the lon and lat set to missing value
          WHERE(clavrLatAnchor.EQ.missingValue) clavrLatAnchor = NAN_R
          ! clavr-x also outputs NaNs, and more than one latitude == 0 in a scan line
          ! (this seems to be only at the end of the orbit and very occasionally),
          ! and also more than one latitude == 90, all meaning reposnx has failed,
          ! probably because the frame contains rubbish.
          ! Don't test for more than one longitude == 0 since that can be valid for
          ! near the pole and the scanline along the Greenwich meridian
#ifdef USE_IEEE
          WHERE(ieee_is_nan(clavrLatAnchor)) clavrLatAnchor = NAN_R
#elif defined(USE_G95)
          WHERE(isNAN(clavrLatAnchor)) clavrLatAnchor = NAN_R
#endif
       ENDIF
       DO line=1,dims(2)
          IF(COUNT(clavrLatAnchor(:,line).EQ.0.0).GT.1 &
               .OR.COUNT(clavrLatAnchor(:,line).EQ.90.0).GT.1 &
               .OR.COUNT(clavrLatAnchor(:,line).EQ.-90.0).GT.1)THEN
             clavrLatAnchor(:,line) = NAN_R
          END IF
       END DO

       CALL h5fclose_f(file_id,hdferr)
       CALL h5close_f(hdferr)

       ! the maximum timing correction magnitude is 1.75 s in the data that is in
       ! clavr-x (which is the same as the data from UMiami) RESTRICTED to the
       ! dates when there is AVHRR data available.  So expect the shift to be
       ! less than 1.75 * 6 * 1.1 = time/s * lines/s * km/line = 11.50 km, say 20km
       ! as a check but a bit tricky

    END IF

  END SUBROUTINE Read_CLAVR_Nav

  SUBROUTINE Read_CLAVR_Cld(outData, cldPathname, nScans)
!--- the output is the cloud mask 0-3, plus a bit set to zero
!--- to say if it's OK, which *here* means OK and is over sea
!--- If not OK then the mask is set to 3, so that you can
!--- simply test bits 1,0 if you're only interested in clear
!--- not in the statistics of the cloud.
!--- Default cld_mask has been changed to 3 to make sure we
!--- don't use unprocessed pixels.
!--- cld_mask is Bayesian
!    cld_mask = 2 byte cloud mask for level 1b output
!		  cloud mask codes:
!			0 - clear
!			1 - probably clear
!			2 - probably cloudy
!			3 - cloud
!
! format is
! bit 2 : 0 cloud mask valid, 1 cloud mask invalid
! bit 1:0 : cld_mask

!
! nScans is from the L1B header and should match CLAVR-x
! the actual size of the L1B data is given by outdata%arraySize
! and this is used to transfer the data

    TYPE(AVHRR_Data), POINTER :: outData
    CHARACTER(LEN=512), INTENT(IN) :: cldPathname
    INTEGER, INTENT(IN) :: nscans

    INTEGER :: pos
    INTEGER :: stat
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: dspace_id
    INTEGER(HID_T) :: attr_id
    INTEGER(HSIZE_T), DIMENSION(2) :: dims
    INTEGER(HSIZE_T), DIMENSION(2) :: maxdims
    INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
    INTEGER :: rank
    INTEGER :: hdferr
    REAL :: missingValue
    LOGICAL :: exists
    INTEGER :: line

    attr_dims = 0

    outData%Clavrx_There = .FALSE.
    outdata%clavrx_mask = 7 ! default is unprocessed and cloudy

    IF (.NOT. outData%isGac ) THEN
       CALL out_WARNING('CLAVR-x cloud only pre-calculated for GAC',&
            'Read_CLAVR_Cld','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       RETURN
    END IF

    INQUIRE(file=TRIM(cldPathname), exist=exists)
    IF (.NOT. exists) THEN
       CALL out_WARNING('No CLAVR-x cloud file',&
            'Read_CLAVR_Cld','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       RETURN
    END IF

    CALL h5open_f(hdferr)
    CALL h5fopen_f(cldPathname,H5F_ACC_RDONLY_F,file_id,hdferr)

    CALL h5dopen_f(file_id,'/clavr-x_cld_mask_results',dset_id,hdferr)
    CALL h5dget_space_f(dset_id,dspace_id,hdferr)
    CALL h5sget_simple_extent_ndims_f(dspace_id,rank,hdferr)
    IF (rank.NE.2) CALL out_ERROR('CLAVR-x input not 2D',&
         'Read_CLAVR_Cld','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    CALL h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,hdferr)
    CALL h5sclose_f(dspace_id,hdferr)
    ! should be 409 across track : dims(1)
    IF (409.NE.dims(1)) CALL out_ERROR('CLAVR-x data not 409 across track',&
         'Read_CLAVR_Cld','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ! and nScans along track: dims(2)
    IF (nScans.NE.dims(2)) CALL out_ERROR('CLAVR-x does not have same number of scans as L1B',&
         'Read_CLAVR_Cld','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    dims(2) = outdata%arraySize ! change to actual number scans in file

    CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,outdata%clavrx_mask,dims,hdferr)
    CALL h5dclose_f(dset_id,hdferr)

    CALL h5fclose_f(file_id,hdferr)
    CALL h5close_f(hdferr)

    outData%Clavrx_There = .TRUE.

  END SUBROUTINE Read_CLAVR_Cld

  SUBROUTINE Read_CLAVR_Prb(outData, prbPathname, nScans)
!--- probability is stored as 0:1 scaled to -127:127
!--- -128 is 'missing value'
!
! nScans is from the L1B header and should match CLAVR-x
! the actual size of the L1B data is given by outdata%arraySize
! and this is used to transfer the data

    TYPE(AVHRR_Data), POINTER :: outData
    CHARACTER(LEN=512), INTENT(IN) :: prbPathname
    INTEGER, INTENT(IN) :: nscans

    INTEGER :: pos
    INTEGER :: stat
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: dspace_id
    INTEGER(HID_T) :: attr_id
    INTEGER(HSIZE_T), DIMENSION(2) :: dims
    INTEGER(HSIZE_T), DIMENSION(2) :: maxdims
    INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
    INTEGER :: rank
    INTEGER :: hdferr
    LOGICAL :: exists
    INTEGER :: line

    attr_dims = 0

    outData%Clavrx_There = .FALSE.
    outdata%clavrx_prb = -128 ! default is unprocessed

    IF (.NOT. outData%isGac ) THEN
       CALL out_WARNING('CLAVR-x cloud only pre-calculated for GAC',&
            'Read_CLAVR_Prb','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       RETURN
    END IF

    INQUIRE(file=TRIM(prbPathname), exist=exists)
    IF (.NOT. exists) THEN
       CALL out_WARNING('No CLAVR-x cloud probability file',&
            'Read_CLAVR_Prb','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       RETURN
    END IF

    CALL h5open_f(hdferr)
    CALL h5fopen_f(prbPathname,H5F_ACC_RDONLY_F,file_id,hdferr)

    CALL h5dopen_f(file_id,'/cloud_probability',dset_id,hdferr)
    CALL h5dget_space_f(dset_id,dspace_id,hdferr)
    CALL h5sget_simple_extent_ndims_f(dspace_id,rank,hdferr)
    IF (rank.NE.2) CALL out_ERROR('CLAVR-x input not 2D',&
         'Read_CLAVR_Prb','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    CALL h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,hdferr)
    CALL h5sclose_f(dspace_id,hdferr)
    ! should be 409 across track : dims(1)
    IF (409.NE.dims(1)) CALL out_ERROR('CLAVR-x data not 409 across track',&
         'Read_CLAVR_Prb','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ! and nScans along track: dims(2)
    IF (nScans.NE.dims(2)) CALL out_ERROR('CLAVR-x does not have same number of scans as L1B',&
         'Read_CLAVR_Prb','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    dims(2) = outdata%arraySize ! change to actual number scans in file

    CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,outData%clavrx_prb,dims,hdferr)
    CALL h5dclose_f(dset_id,hdferr)

    CALL h5fclose_f(file_id,hdferr)
    CALL h5close_f(hdferr)

    outData%Clavrx_There = .TRUE.

  END SUBROUTINE Read_CLAVR_Prb
#endif

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

  SUBROUTINE Reallocate_outData( outData, extraLines )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: extraLines

    INTEGER :: numberOfLines

    CALL Reallocate_Arrays( outData%scanLineNumber, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%badTime, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%badNavigation, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%badCalibration, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%transition3A3B, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Lon, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Lat, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%satZA, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%solZA, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%relAz, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%counts1, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%counts2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%counts3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%counts4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%counts5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array1, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array3A, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array3B, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%array5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%year, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%month, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%day, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dayNo, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%hours, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%UTC_msecs, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%time, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt1, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt1Counts, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt2Counts, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt3Counts, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt4Counts, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt1CountsAll, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt2CountsAll, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt3CountsAll, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prt4CountsAll, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bb3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bb4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bb5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%sp3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%sp4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%sp5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bbodyFilter3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bbodyFilter4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%bbodyFilter5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter1, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter3a, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%spaceFilter5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%patch, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%patchExtended, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Radiator, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Cooler, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%a_d_conv, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%motor, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%motorCurrent, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%electronics, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%baseplate, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib1, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib1_2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib1_intercept, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%calib2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib2_2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib2_intercept, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%calib3A, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib3A_2, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib3A_intercept, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%calib3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%calib5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%clavr_mask, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%clavrx_mask, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%clavrx_prb, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%orig_solar_contamination_3B, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%orig_solar_contamination_4, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%orig_solar_contamination_5, extraLines, &
         numberOfLines )
    CALL Reallocate_Arrays( outData%satelliteAltitude, extraLines, &
         numberOfLines )
    outData%arraySize = numberOfLines

    IF( outData%walton_there )THEN
       CALL Reallocate_Arrays( outData%array3B_error, extraLines, numberOfLines )
       CALL Reallocate_Arrays( outData%array4_error, extraLines, numberOfLines )
       CALL Reallocate_Arrays( outData%array5_error, extraLines, numberOfLines )
    ENDIF
CALL Reallocate_Arrays( outData%prtNumber, extraLines, numberOfLines )
!***************
   
    CALL Reallocate_Arrays( outData%prtmean, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prtsigma, extraLines, numberOfLines )
!-----------
!-----------
    CALL Reallocate_Arrays( outData%btf3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%btf4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%btf5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%nef3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%nef4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%nef5, extraLines, numberOfLines )
 
    CALL Reallocate_Arrays( outData%ue3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ue4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ue5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%uce3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%uce4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%uce5, extraLines, numberOfLines )
 
    CALL Reallocate_Arrays( outData%ur3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ur4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ur5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%us3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%us4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%us5, extraLines, numberOfLines )


    CALL Reallocate_Arrays( outData%ustr3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ustr4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ustr5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%utict, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%utstar3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%utstar4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%utstar5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%urict3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%urict4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%urict5, extraLines, numberOfLines )
   
    CALL Reallocate_Arrays( outData%ve3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ve4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ve5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%vstr3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%vstr4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%vstr5, extraLines, numberOfLines )

!************
    CALL Reallocate_Arrays( outData%dre_over_dce3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dce4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dce5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_dcs3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dcs4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dcs5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_dcict3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dcict4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dcict5, extraLines, numberOfLines )  

    CALL Reallocate_Arrays( outData%dre_over_db03, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_db04, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_db05, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_da13, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_da14, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_da15, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_da23, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_da24, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_da25, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_db13, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_db14, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_db15, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_dtinstr3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dtinstr4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_dtinstr5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dre_over_drict3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_drict4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dre_over_drict5, extraLines, numberOfLines )

    CALL Reallocate_Arrays( outData%dtstar_over_daval3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_daval4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_daval5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dbval3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dbval4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dbval5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dnuc3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dnuc4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%dtstar_over_dnuc5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dtstar3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dtstar4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dtstar5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dnuc3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dnuc4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%drict_over_dnuc5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%prtcounts, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ramp_c3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ramp_c4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%ramp_c5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%gain_c3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%gain_c4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%gain_c5, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Rict_c3, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Rict_c4, extraLines, numberOfLines )
    CALL Reallocate_Arrays( outData%Rict_c5, extraLines, numberOfLines )
  END SUBROUTINE Reallocate_outData

  SUBROUTINE Reallocate_Final_outData( outData, finalNlines )

    TYPE(AVHRR_Data), POINTER :: outData
    INTEGER, INTENT(IN) :: finalNlines

    CALL Reallocate_Final_Arrays( outData%scanLineNumber, finalNlines )
    CALL Reallocate_Final_Arrays( outData%badTime, finalNlines )
    CALL Reallocate_Final_Arrays( outData%badNavigation, finalNlines )
    CALL Reallocate_Final_Arrays( outData%badCalibration, finalNlines )
    CALL Reallocate_Final_Arrays( outData%transition3A3B, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Lon, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Lat, finalNlines )
    CALL Reallocate_Final_Arrays( outData%satZA, finalNlines )
    CALL Reallocate_Final_Arrays( outData%solZA, finalNlines )
    CALL Reallocate_Final_Arrays( outData%relAz, finalNlines )
    CALL Reallocate_Final_Arrays( outData%counts1, finalNlines )
    CALL Reallocate_Final_Arrays( outData%counts2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%counts3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%counts4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%counts5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array1, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array3A, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array3B, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%array5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%year, finalNlines )
    CALL Reallocate_Final_Arrays( outData%month, finalNlines )
    CALL Reallocate_Final_Arrays( outData%day, finalNlines )
    CALL Reallocate_Final_Arrays( outData%dayNo, finalNlines )
    CALL Reallocate_Final_Arrays( outData%hours, finalNlines )
    CALL Reallocate_Final_Arrays( outData%UTC_msecs, finalNlines )
    CALL Reallocate_Final_Arrays( outData%time, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt1, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt1Counts, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt2Counts, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt3Counts, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt4Counts, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt1CountsAll, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt2CountsAll, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt3CountsAll, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prt4CountsAll, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bb3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bb4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bb5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%sp3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%sp4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%sp5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bbodyFilter3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bbodyFilter4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%bbodyFilter5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter1, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter3a, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%spaceFilter5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%patch, finalNlines )
    CALL Reallocate_Final_Arrays( outData%patchExtended, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Radiator, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Cooler, finalNlines )
    CALL Reallocate_Final_Arrays( outData%a_d_conv, finalNlines )
    CALL Reallocate_Final_Arrays( outData%motor, finalNlines )
    CALL Reallocate_Final_Arrays( outData%motorCurrent, finalNlines )
    CALL Reallocate_Final_Arrays( outData%electronics, finalNlines )
    CALL Reallocate_Final_Arrays( outData%baseplate, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib1, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib1_2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib1_intercept, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib2_2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib2_intercept, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib3A, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib3A_2, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib3A_intercept, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%calib5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%clavr_mask, finalNlines )
    CALL Reallocate_Final_Arrays( outData%clavrx_mask, finalNlines )
    CALL Reallocate_Final_Arrays( outData%clavrx_prb, finalNlines )
    CALL Reallocate_Final_Arrays( outData%orig_solar_contamination_3B, &
         finalNlines )
    CALL Reallocate_Final_Arrays( outData%orig_solar_contamination_4, &
         finalNlines )
    CALL Reallocate_Final_Arrays( outData%orig_solar_contamination_5, &
         finalNlines )
    CALL Reallocate_Final_Arrays( outData%satelliteAltitude, &
         finalNlines )

    IF( outData%walton_there )THEN
       CALL Reallocate_Final_Arrays( outData%array3B_error, finalNlines )
       CALL Reallocate_Final_Arrays( outData%array4_error, finalNlines )
       CALL Reallocate_Final_Arrays( outData%array5_error, finalNlines )
    ENDIF
  CALL Reallocate_Final_Arrays( outData%prtNumber, finalNlines )
!****************
  
    CALL Reallocate_Final_Arrays( outData%prtmean, finalNlines )
    CALL Reallocate_Final_Arrays( outData%prtsigma, finalNlines )
!----------------
    CALL Reallocate_Final_Arrays( outData%btf3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%btf4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%btf5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%nef3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%nef4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%nef5,finalNlines )
 
    CALL Reallocate_Final_Arrays( outData%ue3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ue4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ue5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%ur3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ur4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ur5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%us3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%us4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%us5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%uce3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%uce4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%uce5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%ustr3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ustr4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ustr5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%utict,finalNlines )

    CALL Reallocate_Final_Arrays( outData%utstar3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%utstar4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%utstar5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%urict3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%urict4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%urict5,finalNlines )
   
    CALL Reallocate_Final_Arrays( outData%ve3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ve4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%ve5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%vstr3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%vstr4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%vstr5,finalNlines )
!************
    CALL Reallocate_Final_Arrays( outData%dre_over_dce3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dce4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dce5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_dcs3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dcs4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dcs5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_dcict3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dcict4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dcict5,finalNlines )  

    CALL Reallocate_Final_Arrays( outData%dre_over_db03,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_db04,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_db05,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_da13,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_da14,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_da15,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_da23,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_da24,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_da25,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_db13,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_db14,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_db15,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_dtinstr3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dtinstr4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_dtinstr5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dre_over_drict3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_drict4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dre_over_drict5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dtstar_over_daval3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_daval4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_daval5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dtstar_over_dbval3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_dbval4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_dbval5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%dtstar_over_dnuc3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_dnuc4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%dtstar_over_dnuc5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%drict_over_dtstar3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%drict_over_dtstar4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%drict_over_dtstar5,finalNlines )

    CALL Reallocate_Final_Arrays( outData%drict_over_dnuc3,finalNlines )
    CALL Reallocate_Final_Arrays( outData%drict_over_dnuc4,finalNlines )
    CALL Reallocate_Final_Arrays( outData%drict_over_dnuc5,finalNlines )
!-----------
!--Modif du 12 nov
    CALL Reallocate_Final_Arrays( outData%prtcounts, finalNlines )
!-----------------
!----------------
!--Modif du 14 oct
    CALL Reallocate_Final_Arrays( outData%ramp_c3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%ramp_c4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%ramp_c5, finalNlines )
    CALL Reallocate_Final_Arrays( outData%gain_c3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%gain_c4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%gain_c5, finalNlines )
!-----------------
    CALL Reallocate_Final_Arrays( outData%Rict_c3, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Rict_c4, finalNlines )
    CALL Reallocate_Final_Arrays( outData%Rict_c5, finalNlines )
  
    outData%arraySize = finalNlines

  END SUBROUTINE Reallocate_Final_outData

  SUBROUTINE Allocate_OutData( ndata, outData )

    INTEGER, INTENT(IN) :: ndata
    TYPE(AVHRR_Data), POINTER :: outData

    INTEGER :: STAT

    outData%arraySize = MEMORY_ALLOC_STEP
    outData%nelem = ndata

    IF( ALLOCATED(outData%scanLineNumber) )THEN
       write(*,*)'scanLineNumber'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%badTime) )THEN
       write(*,*)'badTime'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%badNavigation) )THEN
       write(*,*)'badNavigation'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%badCalibration) )THEN
       write(*,*)'badCalibration'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%transition3a3b) )THEN
       write(*,*)'transition3A3B'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Lon) )THEN
       write(*,*)'Lon'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Lat) )THEN
       write(*,*)'Lat'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%satZA) )THEN
       write(*,*)'satZA'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%solZA) )THEN
       write(*,*)'solZA'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%relAz) )THEN
       write(*,*)'relAz'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Counts1) )THEN
       write(*,*)'Counts1'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Counts2) )THEN
       write(*,*)'Counts2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Counts3) )THEN
       write(*,*)'Counts3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Counts4) )THEN
       write(*,*)'Counts4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Counts5) )THEN
       write(*,*)'Counts5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array1) )THEN
       write(*,*)'array1'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array2) )THEN
       write(*,*)'array2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array3A) )THEN
       write(*,*)'array3A'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array3B) )THEN
       write(*,*)'array3B'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array4) )THEN
       write(*,*)'array4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%array5) )THEN
       write(*,*)'array5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%year) )THEN
       write(*,*)'year'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%month) )THEN
       write(*,*)'month'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%day) )THEN
       write(*,*)'day'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%dayNo) )THEN
       write(*,*)'dayNo'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%hours) )THEN
       write(*,*)'hours'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%UTC_msecs) )THEN
       write(*,*)'UTC_msecs'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%time) )THEN
       write(*,*)'time'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt1) )THEN
       write(*,*)'prt1'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt2) )THEN
       write(*,*)'prt2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt3) )THEN
       write(*,*)'prt3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt4) )THEN
       write(*,*)'prt4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt1Counts) )THEN
       write(*,*)'prt1Counts'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt2Counts) )THEN
       write(*,*)'prt2Counts'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt3Counts) )THEN
       write(*,*)'prt3Counts'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt4Counts) )THEN
       write(*,*)'prt4Counts'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt1CountsAll) )THEN
       write(*,*)'prt1CountsAll'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt2CountsAll) )THEN
       write(*,*)'prt2CountsAll'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt3CountsAll) )THEN
       write(*,*)'prt3CountsAll'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%prt4CountsAll) )THEN
       write(*,*)'prt4CountsAll'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bb3) )THEN
       write(*,*)'bb3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bb4) )THEN
       write(*,*)'bb4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bb5) )THEN
       write(*,*)'bb5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%sp3) )THEN
       write(*,*)'sp3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%sp4) )THEN
       write(*,*)'sp4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%sp5) )THEN
       write(*,*)'sp5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bbodyFilter3) )THEN
       write(*,*)'bbodyFilter3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bbodyFilter4) )THEN
       write(*,*)'bbodyFilter4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%bbodyFilter5) )THEN
       write(*,*)'bbodyFilter5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter1) )THEN
       write(*,*)'spaceFilter1'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter2) )THEN
       write(*,*)'spaceFilter2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter3a) )THEN
       write(*,*)'spaceFilter3a'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter3) )THEN
       write(*,*)'spaceFilter3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter4) )THEN
       write(*,*)'spaceFilter4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%spaceFilter5) )THEN
       write(*,*)'spaceFilter5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%patch) )THEN
       write(*,*)'patch'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%patchExtended) )THEN
       write(*,*)'patchExtended'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Radiator) )THEN
       write(*,*)'Radiator'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%Cooler) )THEN
       write(*,*)'Cooler'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%a_d_conv) )THEN
       write(*,*)'a_d_conv'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%motor) )THEN
       write(*,*)'motor'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%motorCurrent) )THEN
       write(*,*)'motorCurrent'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%electronics) )THEN
       write(*,*)'electronics'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%baseplate) )THEN
       write(*,*)'baseplate'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib1) )THEN
       write(*,*)'calib1'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib1_2) )THEN
       write(*,*)'calib1_2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib1_intercept) )THEN
       write(*,*)'calib1_intercept'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib2) )THEN
       write(*,*)'calib2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib2_2) )THEN
       write(*,*)'calib2_2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib2_intercept) )THEN
       write(*,*)'calib2_intercept'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib3A) )THEN
       write(*,*)'calib3A'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib3A_2) )THEN
       write(*,*)'calib3A_2'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib3A_intercept) )THEN
       write(*,*)'calib3A_intercept'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib3) )THEN
       write(*,*)'calib3'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib4) )THEN
       write(*,*)'calib4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%calib5) )THEN
       write(*,*)'calib5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%clavr_mask) )THEN
       write(*,*)'clavr_mask'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%clavrx_mask) )THEN
       write(*,*)'clavrx_mask'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%clavrx_prb) )THEN
       write(*,*)'clavrx_prb'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%orig_solar_contamination_3B) )THEN
       write(*,*)'orig_solar_contamination_3B'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%orig_solar_contamination_4) )THEN
       write(*,*)'orig_solar_contamination_4'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%orig_solar_contamination_5) )THEN
       write(*,*)'orig_solar_contamination_5'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF
    IF( ALLOCATED(outData%satelliteAltitude) )THEN
       write(*,*)'satelliteAltitude'
       CALL out_ERROR('Array should not be allocated',&
            'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    IF( outData%walton_there )THEN
       IF( ALLOCATED(outData%array3B_error) )THEN
          write(*,*)'array3B_error'
          CALL out_ERROR('Array should not be allocated',&
               'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       IF( ALLOCATED(outData%array4_error) )THEN
          write(*,*)'array4_error'
          CALL out_ERROR('Array should not be allocated',&
               'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       IF( ALLOCATED(outData%array5_error) )THEN
          write(*,*)'array5_error'
          CALL out_ERROR('Array should not be allocated',&
               'Allocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
    ENDIF

    ALLOCATE(outData%scanLineNumber(MEMORY_ALLOC_STEP),&
         outData%badTime(MEMORY_ALLOC_STEP),&
         outData%badNavigation(MEMORY_ALLOC_STEP),&
         outData%badCalibration(MEMORY_ALLOC_STEP),&
         outData%transition3A3B(MEMORY_ALLOC_STEP),&
         outData%Lon(ndata,MEMORY_ALLOC_STEP),&
         outData%Lat(ndata,MEMORY_ALLOC_STEP),&
         outData%satZA(ndata,MEMORY_ALLOC_STEP),&
         outData%solZA(ndata,MEMORY_ALLOC_STEP),&
         outData%relAz(ndata,MEMORY_ALLOC_STEP),&
         outData%Counts1(ndata,MEMORY_ALLOC_STEP),&
         outData%Counts2(ndata,MEMORY_ALLOC_STEP),&
         outData%Counts3(ndata,MEMORY_ALLOC_STEP),&
         outData%Counts4(ndata,MEMORY_ALLOC_STEP),&
         outData%Counts5(ndata,MEMORY_ALLOC_STEP),&
         outData%array1(ndata,MEMORY_ALLOC_STEP),&
         outData%array2(ndata,MEMORY_ALLOC_STEP),&
         outData%array3A(ndata,MEMORY_ALLOC_STEP),&
         outData%array3B(ndata,MEMORY_ALLOC_STEP),&
         outData%array4(ndata,MEMORY_ALLOC_STEP),&
         outData%array5(ndata,MEMORY_ALLOC_STEP),&
         outData%year(MEMORY_ALLOC_STEP),&
         outData%month(MEMORY_ALLOC_STEP),&
         outData%day(MEMORY_ALLOC_STEP),&
         outData%dayNo(MEMORY_ALLOC_STEP),&
         outData%hours(MEMORY_ALLOC_STEP),&
         outData%UTC_msecs(MEMORY_ALLOC_STEP),&
         outData%time(MEMORY_ALLOC_STEP),&
         outData%prt1(MEMORY_ALLOC_STEP),&
         outData%prt2(MEMORY_ALLOC_STEP),&
         outData%prt3(MEMORY_ALLOC_STEP),&
         outData%prt4(MEMORY_ALLOC_STEP),&
         outData%prt1Counts(MEMORY_ALLOC_STEP),&
         outData%prt2Counts(MEMORY_ALLOC_STEP),&
         outData%prt3Counts(MEMORY_ALLOC_STEP),&
         outData%prt4Counts(MEMORY_ALLOC_STEP),&
         outData%prt1CountsAll(NPRT_IND,MEMORY_ALLOC_STEP),&
         outData%prt2CountsAll(NPRT_IND,MEMORY_ALLOC_STEP),&
         outData%prt3CountsAll(NPRT_IND,MEMORY_ALLOC_STEP),&
         outData%prt4CountsAll(NPRT_IND,MEMORY_ALLOC_STEP),&
         outData%bb3(MEMORY_ALLOC_STEP),&
         outData%bb4(MEMORY_ALLOC_STEP),&
         outData%bb5(MEMORY_ALLOC_STEP),&
         outData%sp3(MEMORY_ALLOC_STEP),&
         outData%sp4(MEMORY_ALLOC_STEP),&
         outData%sp5(MEMORY_ALLOC_STEP),&
         outData%bbodyFilter3(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%bbodyFilter4(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%bbodyFilter5(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter1(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter2(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter3a(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter3(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter4(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%spaceFilter5(NBBODY_SAMPLES,MEMORY_ALLOC_STEP),&
         outData%patch(MEMORY_ALLOC_STEP),&
         outData%patchExtended(MEMORY_ALLOC_STEP),&
         outData%Radiator(MEMORY_ALLOC_STEP),&
         outData%Cooler(MEMORY_ALLOC_STEP),&
         outData%a_d_conv(MEMORY_ALLOC_STEP),&
         outData%motor(MEMORY_ALLOC_STEP),&
         outData%motorCurrent(MEMORY_ALLOC_STEP),&
         outData%electronics(MEMORY_ALLOC_STEP),&
         outData%baseplate(MEMORY_ALLOC_STEP),&
         outData%calib1(2,MEMORY_ALLOC_STEP),&
         outData%calib1_2(2,MEMORY_ALLOC_STEP),&
         outData%calib1_intercept(MEMORY_ALLOC_STEP),&
         outData%calib2(2,MEMORY_ALLOC_STEP),&
         outData%calib2_2(2,MEMORY_ALLOC_STEP),&
         outData%calib2_intercept(MEMORY_ALLOC_STEP),&
         outData%calib3A(2,MEMORY_ALLOC_STEP),&
         outData%calib3A_2(2,MEMORY_ALLOC_STEP),&
         outData%calib3A_intercept(MEMORY_ALLOC_STEP),&
         outData%calib3(3,MEMORY_ALLOC_STEP),&
         outData%calib4(3,MEMORY_ALLOC_STEP),&
         outData%calib5(3,MEMORY_ALLOC_STEP),&
         outData%clavr_mask(ndata,MEMORY_ALLOC_STEP),&
         outData%clavrx_mask(ndata,MEMORY_ALLOC_STEP),&
         outData%clavrx_prb(ndata,MEMORY_ALLOC_STEP),&
         outData%orig_solar_contamination_3B(MEMORY_ALLOC_STEP),&
         outData%orig_solar_contamination_4(MEMORY_ALLOC_STEP),&
         outData%orig_solar_contamination_5(MEMORY_ALLOC_STEP),&
         outData%satelliteAltitude(MEMORY_ALLOC_STEP),&
outData%prtNumber(MEMORY_ALLOC_STEP),&
!**************
       
!*************************
! modif du 20 aout 2015
         outData%prtmean(MEMORY_ALLOC_STEP),&
         outData%prtsigma(MEMORY_ALLOC_STEP),&
!*************************
!*************************
         outData%btf3(ndata,MEMORY_ALLOC_STEP),&
         outData%btf4(ndata,MEMORY_ALLOC_STEP),&
         outData%btf5(ndata,MEMORY_ALLOC_STEP),&
         outData%nef3(ndata,MEMORY_ALLOC_STEP),&
         outData%nef4(ndata,MEMORY_ALLOC_STEP),&
         outData%nef5(ndata,MEMORY_ALLOC_STEP),&

!*****uncertainties
         outData%ue3(ndata, MEMORY_ALLOC_STEP),&
         outData%ue4(ndata, MEMORY_ALLOC_STEP),&
         outData%ue5(ndata, MEMORY_ALLOC_STEP),&

         outData%ur3(ndata, MEMORY_ALLOC_STEP),&
         outData%ur4(ndata, MEMORY_ALLOC_STEP),&
         outData%ur5(ndata, MEMORY_ALLOC_STEP),&

         outData%us3(ndata, MEMORY_ALLOC_STEP),&
         outData%us4(ndata, MEMORY_ALLOC_STEP),&
         outData%us5(ndata, MEMORY_ALLOC_STEP),&
  
         outData%uce3(ndata,MEMORY_ALLOC_STEP),&
         outData%uce4(ndata,MEMORY_ALLOC_STEP),&
         outData%uce5(ndata,MEMORY_ALLOC_STEP),&
  
         outData%ustr3(ndata,MEMORY_ALLOC_STEP),&
         outData%ustr4(ndata,MEMORY_ALLOC_STEP),&
         outData%ustr5(ndata,MEMORY_ALLOC_STEP),&
         
         outData%utict(MEMORY_ALLOC_STEP),&
         
         outData%utstar3(MEMORY_ALLOC_STEP),&
         outData%utstar4(MEMORY_ALLOC_STEP),&
         outData%utstar5(MEMORY_ALLOC_STEP),&

         outData%urict3(MEMORY_ALLOC_STEP),&
         outData%urict4(MEMORY_ALLOC_STEP),&
         outData%urict5(MEMORY_ALLOC_STEP),&
         
         outData%ve3(ndata,MEMORY_ALLOC_STEP),&
         outData%ve4(ndata,MEMORY_ALLOC_STEP),&
         outData%ve5(ndata,MEMORY_ALLOC_STEP),&
         outData%vstr3(ndata,MEMORY_ALLOC_STEP),&
         outData%vstr4(ndata,MEMORY_ALLOC_STEP),&
         outData%vstr5(ndata,MEMORY_ALLOC_STEP),&
       
!******
         outData%dre_over_dce3(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dce4(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dce5(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_dcs3(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dcs4(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dcs5(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_dcict3(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dcict4(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dcict5(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_db03(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_db04(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_db05(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_da13(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_da14(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_da15(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_da23(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_da24(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_da25(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_db13(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_db14(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_db15(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_dtinstr3(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dtinstr4(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_dtinstr5(ndata,MEMORY_ALLOC_STEP),&

         outData%dre_over_drict3(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_drict4(ndata,MEMORY_ALLOC_STEP),&
         outData%dre_over_drict5(ndata,MEMORY_ALLOC_STEP),&
        
         outData%dtstar_over_daval3(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_daval4(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_daval5(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dbval3(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dbval4(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dbval5(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dnuc3(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dnuc4(MEMORY_ALLOC_STEP),&
         outData%dtstar_over_dnuc5(MEMORY_ALLOC_STEP),&
         outData%drict_over_dnuc3(MEMORY_ALLOC_STEP), &
         outData%drict_over_dnuc4(MEMORY_ALLOC_STEP), &
         outData%drict_over_dnuc5(MEMORY_ALLOC_STEP), &
         outData%drict_over_dtstar3(MEMORY_ALLOC_STEP),&
         outData%drict_over_dtstar4(MEMORY_ALLOC_STEP),&
         outData%drict_over_dtstar5(MEMORY_ALLOC_STEP), &
!**********************
!**********************
! modif du 12 nov 2015
         outData%prtcounts(3,MEMORY_ALLOC_STEP),&
!************************
!**********************
! modif du 14 oct 2014
         outData%ramp_c3(MEMORY_ALLOC_STEP),&
         outData%ramp_c4(MEMORY_ALLOC_STEP),&
         outData%ramp_c5(MEMORY_ALLOC_STEP),&
         outData%gain_c3(MEMORY_ALLOC_STEP),&
         outData%gain_c4(MEMORY_ALLOC_STEP),&
         outData%gain_c5(MEMORY_ALLOC_STEP),&
!*************************
! modif du 27 aout 2015
         outData%Rict_c3(MEMORY_ALLOC_STEP),&
         outData%Rict_c4(MEMORY_ALLOC_STEP),&
         outData%Rict_c5(MEMORY_ALLOC_STEP),&
         STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outData','Allocate_OutData',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    IF( outData%walton_there )THEN
       ALLOCATE(outData%array3B(ndata,MEMORY_ALLOC_STEP),&
            outData%array4(ndata,MEMORY_ALLOC_STEP),&
            outData%array5(ndata,MEMORY_ALLOC_STEP),&
            STAT=STAT)
       CALL Check_IOS(STAT,'Allocating outData (walton)','Allocate_OutData',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

  END SUBROUTINE Allocate_OutData

  SUBROUTINE Deallocate_OutData( outData )

    TYPE(AVHRR_Data), target, INTENT(INOUT) :: outData

    TYPE(AVHRR_Data), POINTER :: outData_ptr => NULL()

    outData_Ptr => outData
    CALL Deallocate_OutData_pntr( outData_Ptr )

  END SUBROUTINE Deallocate_OutData

  SUBROUTINE Deallocate_OutData_pntr( outData )

    TYPE(AVHRR_Data), POINTER :: outData

    INTEGER :: STAT

    outData%dataFilled = .FALSE.
    outData%arraySize = 0
    IF( ALLOCATED(outData%scanLineNumber) )THEN
       DEALLOCATE(outData%scanLineNumber,&
            outData%badTime,&
            outData%badNavigation,&
            outData%badCalibration,&
            outData%transition3A3B,&
            outData%Lon,&
            outData%Lat,&
            outData%satZA,&
            outData%solZA,&
            outData%relAz,&
            outData%Counts1,&
            outData%Counts2,&
            outData%Counts3,&
            outData%Counts4,&
            outData%Counts5,&
            outData%array1,&
            outData%array2,&
            outData%array3A,&
            outData%array3B,&
            outData%array4,&
            outData%array5,&
            outData%year,&
            outData%month,&
            outData%day,&
            outData%dayNo,&
            outData%hours,&
            outData%UTC_msecs,&
            outData%time,&
            outData%prt1,&
            outData%prt2,&
            outData%prt3,&
            outData%prt4,&
            outData%prt1Counts,&
            outData%prt2Counts,&
            outData%prt3Counts,&
            outData%prt4Counts,&
            outData%prt1CountsAll,&
            outData%prt2CountsAll,&
            outData%prt3CountsAll,&
            outData%prt4CountsAll,&
            outData%bb3,&
            outData%bb4,&
            outData%bb5,&
            outData%sp3,&
            outData%sp4,&
            outData%sp5,&
            outData%bbodyFilter3,&
            outData%bbodyFilter4,&
            outData%bbodyFilter5,&
            outData%spaceFilter1,&
            outData%spaceFilter2,&
            outData%spaceFilter3a,&
            outData%spaceFilter3,&
            outData%spaceFilter4,&
            outData%spaceFilter5,&
            outData%patch,&
            outData%patchExtended,&
            outData%Radiator,&
            outData%Cooler,&
            outData%a_d_conv,&
            outData%motor,&
            outData%motorCurrent,&
            outData%electronics,&
            outData%baseplate,&
            outData%calib1,&
            outData%calib1_2,&
            outData%calib1_intercept,&
            outData%calib2,&
            outData%calib2_2,&
            outData%calib2_intercept,&
            outData%calib3A,&
            outData%calib3A_2,&
            outData%calib3A_intercept,&
            outData%calib3,&
            outData%calib4,&
            outData%calib5,&
            outData%clavr_mask,&
            outData%clavrx_mask,&
            outData%clavrx_prb,&
            outData%orig_solar_contamination_3B,&
            outData%orig_solar_contamination_4,&
            outData%orig_solar_contamination_5,&
            outData%satelliteAltitude,&
  outData%prtNumber,&
!*****************
           
!************************
! modif du 14 oct 2015
            outData%ramp_c3,&
            outData%ramp_c4,&
            outData%ramp_c5,&
            outData%gain_c3,&
            outData%gain_c4,&
            outData%gain_c5,&
!*************************
!************************
! modif du 12 nov 2015
            outData%prtcounts,&
!*************************
! modif du 20 aout 2015
            outData%prtmean,&
            outData%prtsigma,&
!*************************
!*************************
            outData%btf3,&
            outData%btf4,&
            outData%btf5,&
            outData%nef3,&
            outData%nef4,&
            outData%nef5,&

!*****uncertainties
            outData%ue3,&
            outData%ue4,&
            outData%ue5,&

            outData%us3,&
            outData%us4,&
            outData%us5,&
            outData%ur3,&
            outData%ur4,&
            outData%ur5,&
  
            outData%uce3,&
            outData%uce4,&
            outData%uce5,&
  
            outData%ustr3,&
            outData%ustr4,&
            outData%ustr5,&
         
            outData%utict,&
         
            outData%utstar3,&
            outData%utstar4,&
            outData%utstar5,&

            outData%urict3,&
            outData%urict4,&
            outData%urict5,&
         
            outData%ve3,&
            outData%ve4,&
            outData%ve5,&
            outData%vstr3,&
            outData%vstr4,&
            outData%vstr5,&
!******
            outData%dre_over_dce3,&
            outData%dre_over_dce4,&
            outData%dre_over_dce5,&

            outData%dre_over_dcs3,&
            outData%dre_over_dcs4,&
            outData%dre_over_dcs5,&

            outData%dre_over_dcict3,&
            outData%dre_over_dcict4,&
            outData%dre_over_dcict5,&

            outData%dre_over_db03,&
            outData%dre_over_db04,&
            outData%dre_over_db05,&

            outData%dre_over_da13,&
            outData%dre_over_da14,&
            outData%dre_over_da15,&

            outData%dre_over_da23,&
            outData%dre_over_da24,&
            outData%dre_over_da25,&

            outData%dre_over_db13,&
            outData%dre_over_db14,&
            outData%dre_over_db15,&

            outData%dre_over_dtinstr3,&
            outData%dre_over_dtinstr4,&
            outData%dre_over_dtinstr5,&

            outData%dre_over_drict3,&
            outData%dre_over_drict4,&
            outData%dre_over_drict5,&
        
            outData%dtstar_over_daval3,&
            outData%dtstar_over_daval4,&
            outData%dtstar_over_daval5,&
            outData%dtstar_over_dbval3,&
            outData%dtstar_over_dbval4,&
            outData%dtstar_over_dbval5,&
            outData%dtstar_over_dnuc3,&
            outData%dtstar_over_dnuc4,&
            outData%dtstar_over_dnuc5,&

            outData%drict_over_dnuc3, &
            outData%drict_over_dnuc4, &
            outData%drict_over_dnuc5, &
            outData%drict_over_dtstar3,&
            outData%drict_over_dtstar4,&
            outData%drict_over_dtstar5, &
        
!*************************
!*************************
! modif du 27 aout 2015
            outData%Rict_c3,&
            outData%Rict_c4,&
            outData%Rict_c5,&
!*************************
          
            STAT=STAT)

       CALL Check_IOS(STAT,'DeAllocating outData','DeAllocate_OutData',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    !
    ! Walton calibration
    !
    IF( ALLOCATED(outData%array4_error) )THEN
       DEALLOCATE(outData%array3b_error,outData%array4_error,&
            outData%array5_error)
    ENDIF

    if( outData%newCalibration_There )then
       IF( ALLOCATED(outData%new_calib3) )THEN
          DEALLOCATE(outData%new_calib3,&
               outData%new_calib4,&
               outData%new_calib5,&
               outData%smoothPrt1,&
               outData%smoothPrt2,&
               outData%smoothPrt3,&
               outData%smoothPrt4,&
               outData%smoothPrt1Cnts,&
               outData%smoothPrt2Cnts,&
               outData%smoothPrt3Cnts,&
               outData%smoothPrt4Cnts,&
               outData%smoothPrt,&
               outData%smoothBB3,&
               outData%smoothBB4,&
               outData%smoothBB5,&
               outData%smoothSp3,&
               outData%smoothSp4,&
               outData%smoothSp5,&
               outData%Interpolated,&
               outData%solar_contamination_3B,&
               outData%solar_contamination_4,&
               outData%solar_contamination_5,&
               outData%moon_contamination,&
               outData%new_array1,&
               outData%new_array2,&
               outData%new_array3A,&
               outData%new_array3B,&
               outData%new_array4,&
               outData%new_array5,&
               outData%new_array1_error,&
               outData%new_array2_error,&
               outData%new_array3A_error,&
               outData%new_array3B_error,&
               outData%new_array4_error,&
               outData%new_array5_error,&
               outData%noise_cnts,&

               STAT=STAT)
          CALL Check_IOS(STAT,'DeAllocating outData (recal)',&
               'Deallocate_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       ENDIF
       outData%newCalibration_There = .FALSE.
    endif

  END SUBROUTINE Deallocate_OutData_pntr

  SUBROUTINE Copy_OutData( outData, outputData, start_valid, stop_valid )

    TYPE(AVHRR_Data), POINTER :: outData
    TYPE(AVHRR_Data), INTENT(OUT) :: outputData
    INTEGER, INTENT(IN) :: start_valid
    INTEGER, INTENT(IN) :: stop_valid

    INTEGER :: STAT
    LOGICAL :: walton_there

    if( .not.outData%dataFilled )then
       call out_WARNING('outData not filled with data',&
            'Copy_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       outputData%arraySize = 0
       return
    endif

    IF( outputData%dataFilled )THEN

       outputData%arraySize = 0
       DEALLOCATE(outputData%scanLineNumber,&
            outputData%badTime,&
            outputData%badNavigation,&
            outputData%badCalibration,&
            outputData%transition3A3B,&
            outputData%Lon,&
            outputData%Lat,&
            outputData%satZA,&
            outputData%solZA,&
            outputData%relAz,&
            outputData%Counts1,&
            outputData%Counts2,&
            outputData%Counts3,&
            outputData%Counts4,&
            outputData%Counts5,&
            outputData%array1,&
            outputData%array2,&
            outputData%array3A,&
            outputData%array3B,&
            outputData%array4,&
            outputData%array5,&
            outputData%year,&
            outputData%month,&
            outputData%day,&
            outputData%dayNo,&
            outputData%hours,&
            outputData%UTC_msecs,&
            outputData%time,&
            outputData%prt1,&
            outputData%prt2,&
            outputData%prt3,&
            outputData%prt4,&
            outputData%prt1Counts,&
            outputData%prt2Counts,&
            outputData%prt3Counts,&
            outputData%prt4Counts,&
            outputData%prt1CountsAll,&
            outputData%prt2CountsAll,&
            outputData%prt3CountsAll,&
            outputData%prt4CountsAll,&
            outputData%bb3,&
            outputData%bb4,&
            outputData%bb5,&
            outputData%sp3,&
            outputData%sp4,&
            outputData%sp5,&
            outputData%bbodyFilter3,&
            outputData%bbodyFilter4,&
            outputData%bbodyFilter5,&
            outputData%spaceFilter1,&
            outputData%spaceFilter2,&
            outputData%spaceFilter3a,&
            outputData%spaceFilter3,&
            outputData%spaceFilter4,&
            outputData%spaceFilter5,&
            outputData%patch,&
            outputData%patchExtended,&
            outputData%Radiator,&
            outputData%Cooler,&
            outputData%a_d_conv,&
            outputData%motor,&
            outputData%motorCurrent,&
            outputData%electronics,&
            outputData%baseplate,&
            outputData%calib1,&
            outputData%calib1_2,&
            outputData%calib1_intercept,&
            outputData%calib2,&
            outputData%calib2_2,&
            outputData%calib2_intercept,&
            outputData%calib3A,&
            outputData%calib3A_2,&
            outputData%calib3A_intercept,&
            outputData%calib3,&
            outputData%calib4,&
            outputData%calib5,&
            outputData%clavr_mask,&
            outputData%clavrx_mask,&
            outputData%clavrx_prb,&
            outputData%orig_solar_contamination_3B,&
            outputData%orig_solar_contamination_4,&
            outputData%orig_solar_contamination_5,&
            outData%satelliteAltitude,&
 outputData%ramp_c3,&
            outputData%ramp_c4,&
            outputData%ramp_c5,&
            outputData%gain_c3,&
            outputData%gain_c4,&
            outputData%gain_c5,&
!***********************
!*************************
! modif du 20 aout 2015
            outData%prtmean, &
            outData%prtsigma, &
!**************************
!**************************
            outData%btf3, &
            outData%btf4, &
            outData%btf5, &
            outData%nef3, &
            outData%nef4, &
            outData%nef5, &
!*****uncertainties
            outData%ue3,&
            outData%ue4,&
            outData%ue5,&

            outData%ur3,&
            outData%ur4,&
            outData%ur5,&
            outData%us3,&
            outData%us4,&
            outData%us5,&
  
            outData%uce3,&
            outData%uce4,&
            outData%uce5,&
  
            outData%ustr3,&
            outData%ustr4,&
            outData%ustr5,&
         
            outData%utict,&
         
            outData%utstar3,&
            outData%utstar4,&
            outData%utstar5,&

            outData%urict3,&
            outData%urict4,&
            outData%urict5,&
         
            outData%ve3,&
            outData%ve4,&
            outData%ve5,&
            outData%vstr3,&
            outData%vstr4,&
            outData%vstr5,&
  !******
            outData%dre_over_dce3,&
            outData%dre_over_dce4,&
            outData%dre_over_dce5,&

            outData%dre_over_dcs3,&
            outData%dre_over_dcs4,&
            outData%dre_over_dcs5,&

            outData%dre_over_dcict3,&
            outData%dre_over_dcict4,&
            outData%dre_over_dcict5,&

            outData%dre_over_db03,&
            outData%dre_over_db04,&
            outData%dre_over_db05,&

            outData%dre_over_da13,&
            outData%dre_over_da14,&
            outData%dre_over_da15,&

            outData%dre_over_da23,&
            outData%dre_over_da24,&
            outData%dre_over_da25,&

            outData%dre_over_db13,&
            outData%dre_over_db14,&
            outData%dre_over_db15,&

            outData%dre_over_dtinstr3,&
            outData%dre_over_dtinstr4,&
            outData%dre_over_dtinstr5,&

            outData%dre_over_drict3,&
            outData%dre_over_drict4,&
            outData%dre_over_drict5,&
        
            outData%dtstar_over_daval3,&
            outData%dtstar_over_daval4,&
            outData%dtstar_over_daval5,&
            outData%dtstar_over_dbval3,&
            outData%dtstar_over_dbval4,&
            outData%dtstar_over_dbval5,&
            outData%dtstar_over_dnuc3,&
            outData%dtstar_over_dnuc4,&
            outData%dtstar_over_dnuc5,&

            outData%drict_over_dnuc3, &
            outData%drict_over_dnuc4, &
            outData%drict_over_dnuc5, &
            outData%drict_over_dtstar3,&
            outData%drict_over_dtstar4,&
            outData%drict_over_dtstar5, &
!*************************          
! modif du 27 aout 2015
            outData%Rict_c3,&
            outData%Rict_c4,&
            outData%Rict_c5,&
            STAT=STAT)
       CALL Check_IOS(STAT,'DeAllocating outputData','Copy_OutData',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       !
       ! Walton calibration
       !
       walton_there = .FALSE.
       IF( ALLOCATED(outData%array4_error) )THEN
          walton_there = .TRUE.
          DEALLOCATE(outData%array3b_error,outData%array4_error,&
               outData%array5_error)
       ENDIF

       if( outputData%newCalibration_There )then
          DEALLOCATE(outputData%new_calib3,&
               outputData%new_calib4,&
               outputData%new_calib5,&
               outputData%smoothPrt1,&
               outputData%smoothPrt2,&
               outputData%smoothPrt3,&
               outputData%smoothPrt4,&
               outputData%smoothPrt1Cnts,&
               outputData%smoothPrt2Cnts,&
               outputData%smoothPrt3Cnts,&
               outputData%smoothPrt4Cnts,&
               outputData%smoothPrt,&
               outputData%smoothBB3,&
               outputData%smoothBB4,&
               outputData%smoothBB5,&
               outputData%smoothSp3,&
               outputData%smoothSp4,&
               outputData%smoothSp5,&
               outputData%Interpolated,&
               outputData%solar_contamination_3B,&
               outputData%solar_contamination_4,&
               outputData%solar_contamination_5,&
               outputData%moon_contamination,&
               outputData%new_array1,&
               outputData%new_array2,&
               outputData%new_array3A,&
               outputData%new_array3B,&
               outputData%new_array4,&
               outputData%new_array5,&
               outputData%new_array1_error,&
               outputData%new_array2_error,&
               outputData%new_array3A_error,&
               outputData%new_array3B_error,&
               outputData%new_array4_error,&
               outputData%new_array5_error,&
               STAT=STAT)
          CALL Check_IOS(STAT,'DeAllocating outputData (recal)',&
               'Copy_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
       endif

    ENDIF

    outputData%arraySize = stop_valid-start_valid+1
    outputData%nelem = outData%nelem
    outputData%Clavr_There = outData%Clavr_There
    outputData%Clavrx_There = outData%Clavrx_There
    outputData%start_valid = outData%start_valid
    outputData%stop_valid = outData%stop_valid
    ALLOCATE(outputData%scanLineNumber(outputData%arraySize),&
         outputData%badTime(outputData%arraySize),&
         outputData%badNavigation(outputData%arraySize),&
         outputData%badCalibration(outputData%arraySize),&
         outputData%transition3A3B(outputData%arraySize),&
         outputData%Lon(ndata,outputData%arraySize),&
         outputData%Lat(ndata,outputData%arraySize),&
         outputData%satZA(ndata,outputData%arraySize),&
         outputData%solZA(ndata,outputData%arraySize),&
         outputData%relAz(ndata,outputData%arraySize),&
         outputData%Counts1(ndata,outputData%arraySize),&
         outputData%Counts2(ndata,outputData%arraySize),&
         outputData%Counts3(ndata,outputData%arraySize),&
         outputData%Counts4(ndata,outputData%arraySize),&
         outputData%Counts5(ndata,outputData%arraySize),&
         outputData%array1(ndata,outputData%arraySize),&
         outputData%array2(ndata,outputData%arraySize),&
         outputData%array3A(ndata,outputData%arraySize),&
         outputData%array3B(ndata,outputData%arraySize),&
         outputData%array4(ndata,outputData%arraySize),&
         outputData%array5(ndata,outputData%arraySize),&
         outputData%year(outputData%arraySize),&
         outputData%month(outputData%arraySize),&
         outputData%day(outputData%arraySize),&
         outputData%dayNo(outputData%arraySize),&
         outputData%hours(outputData%arraySize),&
         outputData%UTC_msecs(outputData%arraySize),&
         outputData%time(outputData%arraySize),&
         outputData%prt1(outputData%arraySize),&
         outputData%prt2(outputData%arraySize),&
         outputData%prt3(outputData%arraySize),&
         outputData%prt4(outputData%arraySize),&
         outputData%prt1Counts(outputData%arraySize),&
         outputData%prt2Counts(outputData%arraySize),&
         outputData%prt3Counts(outputData%arraySize),&
         outputData%prt4Counts(outputData%arraySize),&
         outputData%prt1CountsAll(NPRT_IND,outputData%arraySize),&
         outputData%prt2CountsAll(NPRT_IND,outputData%arraySize),&
         outputData%prt3CountsAll(NPRT_IND,outputData%arraySize),&
         outputData%prt4CountsAll(NPRT_IND,outputData%arraySize),&
         outputData%bb3(outputData%arraySize),&
         outputData%bb4(outputData%arraySize),&
         outputData%bb5(outputData%arraySize),&
         outputData%sp3(outputData%arraySize),&
         outputData%sp4(outputData%arraySize),&
         outputData%sp5(outputData%arraySize),&
         outputData%bbodyFilter3(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%bbodyFilter4(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%bbodyFilter5(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter1(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter2(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter3a(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter3(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter4(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%spaceFilter5(NBBODY_SAMPLES,outputData%arraySize),&
         outputData%patch(outputData%arraySize),&
         outputData%patchExtended(outputData%arraySize),&
         outputData%Radiator(outputData%arraySize),&
         outputData%Cooler(outputData%arraySize),&
         outputData%a_d_conv(outputData%arraySize),&
         outputData%motor(outputData%arraySize),&
         outputData%motorCurrent(outputData%arraySize),&
         outputData%electronics(outputData%arraySize),&
         outputData%baseplate(outputData%arraySize),&
         outputData%calib1(2,outputData%arraySize),&
         outputData%calib1_2(2,outputData%arraySize),&
         outputData%calib1_intercept(outputData%arraySize),&
         outputData%calib2(2,outputData%arraySize),&
         outputData%calib2_2(2,outputData%arraySize),&
         outputData%calib2_intercept(outputData%arraySize),&
         outputData%calib3A(2,outputData%arraySize),&
         outputData%calib3A_2(2,outputData%arraySize),&
         outputData%calib3A_intercept(outputData%arraySize),&
         outputData%calib3(3,outputData%arraySize),&
         outputData%calib4(3,outputData%arraySize),&
         outputData%calib5(3,outputData%arraySize),&
         outputData%clavr_mask(ndata,outputData%arraySize),&
         outputData%clavrx_mask(ndata,outputData%arraySize),&
         outputData%clavrx_prb(ndata,outputData%arraySize),&
         outputData%orig_solar_contamination_3B(outputData%arraySize),&
         outputData%orig_solar_contamination_4(outputData%arraySize),&
         outputData%orig_solar_contamination_5(outputData%arraySize),&
         outData%satelliteAltitude(outputData%arraySize),&
outData%prtmean(outputData%arraySize),&
         outData%prtsigma(outputData%arraySize),&
!*************************
!*************************
         outputData%btf3(ndata,outputData%arraySize),&
         outputData%btf4(ndata,outputData%arraySize),&
         outputData%btf5(ndata,outputData%arraySize),&
         outputData%nef3(ndata,outputData%arraySize),&
         outputData%nef4(ndata,outputData%arraySize),&
         outputData%nef5(ndata,outputData%arraySize),&
!*****uncertainties
         outData%ue3(ndata,outputData%arraySize),&
         outData%ue4(ndata,outputData%arraySize),&
         outData%ue5(ndata,outputData%arraySize),&
 
         outData%ur3(ndata,outputData%arraySize),&
         outData%ur4(ndata,outputData%arraySize),&
         outData%ur5(ndata,outputData%arraySize),&
 
         outData%us3(ndata,outputData%arraySize),&
         outData%us4(ndata,outputData%arraySize),&
         outData%us5(ndata,outputData%arraySize),&
 
         outData%uce3(ndata,outputData%arraySize),&
         outData%uce4(ndata,outputData%arraySize),&
         outData%uce5(ndata,outputData%arraySize),&
  
         outData%ustr3(ndata,outputData%arraySize),&
         outData%ustr4(ndata,outputData%arraySize),&
         outData%ustr5(ndata,outputData%arraySize),&
         
         outData%utict(outputData%arraySize),&
         
         outData%utstar3(outputData%arraySize),&
         outData%utstar4(outputData%arraySize),&
         outData%utstar5(outputData%arraySize),&

         outData%urict3(outputData%arraySize),&
         outData%urict4(outputData%arraySize),&
         outData%urict5(outputData%arraySize),&
         
         outData%ve3(ndata,outputData%arraySize),&
         outData%ve4(ndata,outputData%arraySize),&
         outData%ve5(ndata,outputData%arraySize),&
         outData%vstr3(ndata,outputData%arraySize),&
         outData%vstr4(ndata,outputData%arraySize),&
         outData%vstr5(ndata,outputData%arraySize),&
  !******
         outData%dre_over_dce3(ndata,outputData%arraySize),&
         outData%dre_over_dce4(ndata,outputData%arraySize),&
         outData%dre_over_dce5(ndata,outputData%arraySize),&

         outData%dre_over_dcs3(ndata,outputData%arraySize),&
         outData%dre_over_dcs4(ndata,outputData%arraySize),&
         outData%dre_over_dcs5(ndata,outputData%arraySize),&

         outData%dre_over_dcict3(ndata,outputData%arraySize),&
         outData%dre_over_dcict4(ndata,outputData%arraySize),&
         outData%dre_over_dcict5(ndata,outputData%arraySize),&

         outData%dre_over_db03(ndata,outputData%arraySize),&
         outData%dre_over_db04(ndata,outputData%arraySize),&
         outData%dre_over_db05(ndata,outputData%arraySize),&

         outData%dre_over_da13(ndata,outputData%arraySize),&
         outData%dre_over_da14(ndata,outputData%arraySize),&
         outData%dre_over_da15(ndata,outputData%arraySize),&

         outData%dre_over_da23(ndata,outputData%arraySize),&
         outData%dre_over_da24(ndata,outputData%arraySize),&
         outData%dre_over_da25(ndata,outputData%arraySize),&

         outData%dre_over_db13(ndata,outputData%arraySize),&
         outData%dre_over_db14(ndata,outputData%arraySize),&
         outData%dre_over_db15(ndata,outputData%arraySize),&

         outData%dre_over_dtinstr3(ndata,outputData%arraySize),&
         outData%dre_over_dtinstr4(ndata,outputData%arraySize),&
         outData%dre_over_dtinstr5(ndata,outputData%arraySize),&

         outData%dre_over_drict3(ndata,outputData%arraySize),&
         outData%dre_over_drict4(ndata,outputData%arraySize),&
         outData%dre_over_drict5(ndata,outputData%arraySize),&
        
         outData%dtstar_over_daval3(outputData%arraySize),&
         outData%dtstar_over_daval4(outputData%arraySize),&
         outData%dtstar_over_daval5(outputData%arraySize),&
         outData%dtstar_over_dbval3(outputData%arraySize),&
         outData%dtstar_over_dbval4(outputData%arraySize),&
         outData%dtstar_over_dbval5(outputData%arraySize),&
         outData%dtstar_over_dnuc3(outputData%arraySize),&
         outData%dtstar_over_dnuc4(outputData%arraySize),&
         outData%dtstar_over_dnuc5(outputData%arraySize),&

         outData%drict_over_dnuc3(outputData%arraySize), &
         outData%drict_over_dnuc4(outputData%arraySize), &
         outData%drict_over_dnuc5(outputData%arraySize), &
         outData%drict_over_dtstar3(outputData%arraySize),&
         outData%drict_over_dtstar4(outputData%arraySize),&
         outData%drict_over_dtstar5(outputData%arraySize), &
!*************************
! modif du 27 aout 2015
         outData%Rict_c3(outputData%arraySize),&
         outData%Rict_c4(outputData%arraySize),&
         outData%Rict_c5(outputData%arraySize),&
!*************************
!*************************
! modif du 12 nov 2015
         outData%prtcounts(3,outputData%arraySize),&
!************************
!*************************
! modif du 14 oct 2015
         outData%ramp_c3(outputData%arraySize),&
         outData%ramp_c4(outputData%arraySize),&
         outData%ramp_c5(outputData%arraySize),&
         outData%gain_c3(outputData%arraySize),&
         outData%gain_c4(outputData%arraySize),&
         outData%gain_c5(outputData%arraySize),&

         STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outData','Copy_OutData',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    IF( walton_there )THEN
       ALLOCATE(outputData%array3B_error(ndata,outputData%arraySize),&
            outputData%array4_error(ndata,outputData%arraySize),&
            outputData%array5_error(ndata,outputData%arraySize),&
         STAT=STAT)
       CALL Check_IOS(STAT,'Allocating outData (walton)','Copy_OutData',&
            'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    ENDIF

    if( outData%newCalibration_There )then
       ALLOCATE(outputData%new_calib3(3,outputData%arraySize),&
            outputData%new_calib4(3,outputData%arraySize),&
            outputData%new_calib5(3,outputData%arraySize),&
            outputData%smoothPrt1(outputData%arraySize),&
            outputData%smoothPrt2(outputData%arraySize),&
            outputData%smoothPrt3(outputData%arraySize),&
            outputData%smoothPrt4(outputData%arraySize),&
            outputData%smoothPrt1Cnts(outputData%arraySize),&
            outputData%smoothPrt2Cnts(outputData%arraySize),&
            outputData%smoothPrt3Cnts(outputData%arraySize),&
            outputData%smoothPrt4Cnts(outputData%arraySize),&
            outputData%smoothPrt(outputData%arraySize),&
            outputData%smoothBB3(outputData%arraySize),&
            outputData%smoothBB4(outputData%arraySize),&
            outputData%smoothBB5(outputData%arraySize),&
            outputData%smoothSp3(outputData%arraySize),&
            outputData%smoothSp4(outputData%arraySize),&
            outputData%smoothSp5(outputData%arraySize),&
            outputData%Interpolated(outputData%arraySize),&
            outputData%new_array1(ndata,outputData%arraySize),&
            outputData%new_array2(ndata,outputData%arraySize),&
            outputData%new_array3A(ndata,outputData%arraySize),&
            outputData%new_array3B(ndata,outputData%arraySize),&
            outputData%new_array4(ndata,outputData%arraySize),&
            outputData%new_array5(ndata,outputData%arraySize),&
            outputData%new_array1_error(ndata,outputData%arraySize),&
            outputData%new_array2_error(ndata,outputData%arraySize),&
            outputData%new_array3A_error(ndata,outputData%arraySize),&
            outputData%new_array3B_error(ndata,outputData%arraySize),&
            outputData%new_array4_error(ndata,outputData%arraySize),&
            outputData%new_array5_error(ndata,outputData%arraySize),&
            STAT=STAT)
       CALL Check_IOS(STAT,'Allocating outputData (recal)',&
            'Copy_OutData','NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')
    endif

    outputData%scanLineNumber = outData%scanLineNumber(start_valid:stop_valid)
    outputData%badTime = outData%badTime(start_valid:stop_valid)
    outputData%badNavigation = outData%badNavigation(start_valid:stop_valid)
    outputData%badCalibration = outData%badCalibration(start_valid:stop_valid)
    outputData%transition3A3B = outData%transition3A3B(start_valid:stop_valid)
    outputData%Lon = outData%Lon(:,start_valid:stop_valid)
    outputData%Lat = outData%Lat(:,start_valid:stop_valid)
    outputData%satZA = outData%satZA(:,start_valid:stop_valid)
    outputData%solZA = outData%solZA(:,start_valid:stop_valid)
    outputData%relAz = outData%relAz(:,start_valid:stop_valid)
    outputData%Counts1 = outData%Counts1(:,start_valid:stop_valid)
    outputData%Counts2 = outData%Counts2(:,start_valid:stop_valid)
    outputData%Counts3 = outData%Counts3(:,start_valid:stop_valid)
    outputData%Counts4 = outData%Counts4(:,start_valid:stop_valid)
    outputData%Counts5 = outData%Counts5(:,start_valid:stop_valid)
    outputData%array1 = outData%array1(:,start_valid:stop_valid)
    outputData%array2 = outData%array2(:,start_valid:stop_valid)
    outputData%array3A = outData%array3A(:,start_valid:stop_valid)
    outputData%array3B = outData%array3B(:,start_valid:stop_valid)
    outputData%array4 = outData%array4(:,start_valid:stop_valid)
    outputData%array5 = outData%array5(:,start_valid:stop_valid)
    outputData%year = outData%year(start_valid:stop_valid)
    outputData%month = outData%month(start_valid:stop_valid)
    outputData%day = outData%day(start_valid:stop_valid)
    outputData%dayNo = outData%dayNo(start_valid:stop_valid)
    outputData%hours = outData%hours(start_valid:stop_valid)
    outputData%UTC_msecs = outData%UTC_msecs(start_valid:stop_valid)
    outputData%time = outData%time(start_valid:stop_valid)
    outputData%prt1 = outData%prt1(start_valid:stop_valid)
    outputData%prt2 = outData%prt2(start_valid:stop_valid)
    outputData%prt3 = outData%prt3(start_valid:stop_valid)
    outputData%prt4 = outData%prt4(start_valid:stop_valid)
    outputData%prt1Counts = outData%prt1Counts(start_valid:stop_valid)
    outputData%prt2Counts = outData%prt2Counts(start_valid:stop_valid)
    outputData%prt3Counts = outData%prt3Counts(start_valid:stop_valid)
    outputData%prt4Counts = outData%prt4Counts(start_valid:stop_valid)
    outputData%prt1CountsAll = outData%prt1CountsAll(:,start_valid:stop_valid)
    outputData%prt2CountsAll = outData%prt2CountsAll(:,start_valid:stop_valid)
    outputData%prt3CountsAll = outData%prt3CountsAll(:,start_valid:stop_valid)
    outputData%prt4CountsAll = outData%prt4CountsAll(:,start_valid:stop_valid)
    outputData%bb3 = outData%bb3(start_valid:stop_valid)
    outputData%bb4 = outData%bb4(start_valid:stop_valid)
    outputData%bb5 = outData%bb5(start_valid:stop_valid)
    outputData%sp3 = outData%sp3(start_valid:stop_valid)
    outputData%sp4 = outData%sp4(start_valid:stop_valid)
    outputData%sp5 = outData%sp5(start_valid:stop_valid)
    outputData%bbodyFilter3 = outData%bbodyFilter3(:,start_valid:stop_valid)
    outputData%bbodyFilter4 = outData%bbodyFilter4(:,start_valid:stop_valid)
    outputData%bbodyFilter5 = outData%bbodyFilter5(:,start_valid:stop_valid)
    outputData%spaceFilter1 = outData%spaceFilter1(:,start_valid:stop_valid)
    outputData%spaceFilter2 = outData%spaceFilter2(:,start_valid:stop_valid)
    outputData%spaceFilter3a = outData%spaceFilter3a(:,start_valid:stop_valid)
    outputData%spaceFilter3 = outData%spaceFilter3(:,start_valid:stop_valid)
    outputData%spaceFilter4 = outData%spaceFilter4(:,start_valid:stop_valid)
    outputData%spaceFilter5 = outData%spaceFilter5(:,start_valid:stop_valid)
    outputData%patch = outData%patch(start_valid:stop_valid)
    outputData%patchExtended = outData%patchExtended(start_valid:stop_valid)
    outputData%Radiator = outData%Radiator(start_valid:stop_valid)
    outputData%Cooler = outData%Cooler(start_valid:stop_valid)
    outputData%a_d_conv = outData%a_d_conv(start_valid:stop_valid)
    outputData%motor = outData%motor(start_valid:stop_valid)
    outputData%motorCurrent = outData%motorCurrent(start_valid:stop_valid)
    outputData%electronics = outData%electronics(start_valid:stop_valid)
    outputData%baseplate = outData%baseplate(start_valid:stop_valid)
    outputData%calib1 = outData%calib1(:,start_valid:stop_valid)
    outputData%calib1_2 = outData%calib1_2(:,start_valid:stop_valid)
    outputData%calib1_intercept = &
         outData%calib1_intercept(start_valid:stop_valid)
    outputData%calib2 = outData%calib2(:,start_valid:stop_valid)
    outputData%calib2_2 = outData%calib2_2(:,start_valid:stop_valid)
    outputData%calib2_intercept = &
         outData%calib2_intercept(start_valid:stop_valid)
    outputData%calib3A = outData%calib3A(:,start_valid:stop_valid)
    outputData%calib3A_2 = outData%calib3A_2(:,start_valid:stop_valid)
    outputData%calib3A_intercept = &
         outData%calib3A_intercept(start_valid:stop_valid)
    outputData%calib3 = outData%calib3(:,start_valid:stop_valid)
    outputData%calib4 = outData%calib4(:,start_valid:stop_valid)
    outputData%calib5 = outData%calib5(:,start_valid:stop_valid)
    outputData%clavr_mask = outData%clavr_mask(:,start_valid:stop_valid)
    outputData%clavrx_mask = outData%clavrx_mask(:,start_valid:stop_valid)
    outputData%clavrx_prb = outData%clavrx_prb(:,start_valid:stop_valid)
    outputData%orig_solar_contamination_3B = &
         outData%orig_solar_contamination_3B(start_valid:stop_valid)
    outputData%orig_solar_contamination_4 = &
         outData%orig_solar_contamination_4(start_valid:stop_valid)
    outputData%orig_solar_contamination_5 = &
         outData%orig_solar_contamination_5(start_valid:stop_valid)
    outputData%satelliteAltitude = &
         outData%satelliteAltitude(start_valid:stop_valid)
! modif du 16 nov
    outputData%prtNumber = outData%prtNumber(start_valid:stop_valid)
!****************
  
!*************************
! modif du 20 aout 2015
    outputData%prtmean = outData%prtmean(start_valid:stop_valid)
    outputData%prtsigma = outData%prtmean(start_valid:stop_valid)
!*************************
    outputData%btf3 = outData%btf3(:,start_valid:stop_valid)
    outputData%btf4 = outData%btf4(:,start_valid:stop_valid)
    outputData%btf5 = outData%btf5(:,start_valid:stop_valid)
    outputData%nef3 = outData%nef3(:,start_valid:stop_valid)
    outputData%nef4 = outData%nef4(:,start_valid:stop_valid)
    outputData%nef5 = outData%nef5(:,start_valid:stop_valid)
!*****uncertainties
    outputData%ue3 = outData%ue3(:,start_valid:stop_valid)
    outputData%ue4 = outData%ue4(:,start_valid:stop_valid)
    outputData%ue5 = outData%ue5(:,start_valid:stop_valid)

    outputData%ur3 = outData%ur3(:,start_valid:stop_valid)
    outputData%ur4 = outData%ur4(:,start_valid:stop_valid)
    outputData%ur5 = outData%ur5(:,start_valid:stop_valid)
    
    outputData%us3 = outData%us3(:,start_valid:stop_valid)
    outputData%us4 = outData%us4(:,start_valid:stop_valid)
    outputData%us5 = outData%us5(:,start_valid:stop_valid)
 
    outputData%uce3 = outData%uce3(:,start_valid:stop_valid)
    outputData%uce4 = outData%uce4(:,start_valid:stop_valid)
    outputData%uce5 = outData%uce5(:,start_valid:stop_valid)
  
    outputData%ustr3 = outData%ustr3(:,start_valid:stop_valid)
    outputData%ustr4 = outData%ustr4(:,start_valid:stop_valid)
    outputData%ustr5 = outData%ustr5(:,start_valid:stop_valid)
         
    outputData%utict = outData%utict(start_valid:stop_valid)
         
    outputData%utstar3 = outData%utstar3(start_valid:stop_valid)
    outputData%utstar4 = outData%utstar4(start_valid:stop_valid)
    outputData%utstar5 = outData%utstar5(start_valid:stop_valid)

    outputData%urict3 = outData%urict3(start_valid:stop_valid)
    outputData%urict4 = outData%urict4(start_valid:stop_valid)
    outputData%urict5 = outData%urict5(start_valid:stop_valid)
         
    outputData%ve3 = outData%ve3(:,start_valid:stop_valid)
    outputData%ve4 = outData%ve4(:,start_valid:stop_valid)
    outputData%ve5 = outData%ve5(:,start_valid:stop_valid)
    outputData%vstr3 = outData%vstr3(:,start_valid:stop_valid)
    outputData%vstr4 = outData%vstr4(:,start_valid:stop_valid)
    outputData%vstr5 = outData%vstr5(:,start_valid:stop_valid)

!******
    outData%dre_over_dce3 = outData%dre_over_dce3(:,start_valid:stop_valid)
    outData%dre_over_dce4 = outData%dre_over_dce4(:,start_valid:stop_valid)
    outData%dre_over_dce5 = outData%dre_over_dce5(:,start_valid:stop_valid)

    outData%dre_over_dcs3 = outData%dre_over_dcs3(:,start_valid:stop_valid)
    outData%dre_over_dcs4 = outData%dre_over_dcs4(:,start_valid:stop_valid)
    outData%dre_over_dcs5 = outData%dre_over_dcs5(:,start_valid:stop_valid)

    outData%dre_over_dcict3 = outData%dre_over_dcict3(:,start_valid:stop_valid)
    outData%dre_over_dcict4 = outData%dre_over_dcict4(:,start_valid:stop_valid)
    outData%dre_over_dcict5 = outData%dre_over_dcict5(:,start_valid:stop_valid)

    outData%dre_over_db03 = outData%dre_over_db03(:,start_valid:stop_valid)
    outData%dre_over_db04 = outData%dre_over_db04(:,start_valid:stop_valid)
    outData%dre_over_db05 = outData%dre_over_db05(:,start_valid:stop_valid)

    outData%dre_over_da13 = outData%dre_over_da13(:,start_valid:stop_valid)
    outData%dre_over_da14 = outData%dre_over_da14(:,start_valid:stop_valid)
    outData%dre_over_da15 = outData%dre_over_da15(:,start_valid:stop_valid)

    outData%dre_over_da23 = outData%dre_over_da23(:,start_valid:stop_valid)
    outData%dre_over_da24 = outData%dre_over_da24(:,start_valid:stop_valid)
    outData%dre_over_da25 = outData%dre_over_da25(:,start_valid:stop_valid)

    outData%dre_over_db13 = outData%dre_over_db13(:,start_valid:stop_valid)
    outData%dre_over_db14 = outData%dre_over_db14(:,start_valid:stop_valid)
    outData%dre_over_db15 = outData%dre_over_db15(:,start_valid:stop_valid)

    outData%dre_over_dtinstr3 = outData%dre_over_dtinstr3(:,start_valid:stop_valid)
    outData%dre_over_dtinstr4 = outData%dre_over_dtinstr4(:,start_valid:stop_valid)
    outData%dre_over_dtinstr4 = outData%dre_over_dtinstr5(:,start_valid:stop_valid)

    outData%dre_over_drict3 = outData%dre_over_drict3(:,start_valid:stop_valid)
    outData%dre_over_drict4 = outData%dre_over_drict4(:,start_valid:stop_valid)
    outData%dre_over_drict5 = outData%dre_over_drict5(:,start_valid:stop_valid)
        
    outData%dtstar_over_daval3 = outData%dtstar_over_daval3(start_valid:stop_valid)
    outData%dtstar_over_daval4 = outData%dtstar_over_daval4(start_valid:stop_valid)
    outData%dtstar_over_daval5 = outData%dtstar_over_daval5(start_valid:stop_valid)
    outData%dtstar_over_dbval3 = outData%dtstar_over_dbval3(start_valid:stop_valid)
    outData%dtstar_over_dbval4 = outData%dtstar_over_dbval4(start_valid:stop_valid)
    outData%dtstar_over_dbval5 = outData%dtstar_over_dbval5(start_valid:stop_valid)
    outData%dtstar_over_dnuc3 = outData%dtstar_over_dnuc3(start_valid:stop_valid)
    outData%dtstar_over_dnuc4 = outData%dtstar_over_dnuc4(start_valid:stop_valid)
    outData%dtstar_over_dnuc5 = outData%dtstar_over_dnuc5(start_valid:stop_valid)

    outData%drict_over_dnuc3 = outData%drict_over_dnuc3(start_valid:stop_valid)
    outData%drict_over_dnuc4 = outData%drict_over_dnuc4(start_valid:stop_valid)
    outData%drict_over_dnuc5 = outData%drict_over_dnuc5(start_valid:stop_valid)
    outData%drict_over_dtstar3 = outData%drict_over_dtstar3(start_valid:stop_valid)
    outData%drict_over_dtstar4 = outData%drict_over_dtstar4(start_valid:stop_valid)
    outData%drict_over_dtstar5 = outData%drict_over_dtstar5(start_valid:stop_valid)

!*************************
! modif du 27 aout 2015
    outputData%Rict_c3 = outData%Rict_c3(start_valid:stop_valid)
    outputData%Rict_c4 = outData%Rict_c4(start_valid:stop_valid)
    outputData%Rict_c5 = outData%Rict_c5(start_valid:stop_valid)
!*************************
!*************************
! modif du 14 oct 2015
    outputData%Ramp_c3 = outData%Ramp_c3(start_valid:stop_valid)
    outputData%Ramp_c4 = outData%Ramp_c4(start_valid:stop_valid)
    outputData%Ramp_c5 = outData%Ramp_c5(start_valid:stop_valid)
    outputData%gain_c3 = outData%Ramp_c3(start_valid:stop_valid)
    outputData%gain_c4 = outData%Ramp_c4(start_valid:stop_valid)
    outputData%gain_c5 = outData%Ramp_c5(start_valid:stop_valid)

    if( walton_there )then
       outputData%array3b_error = outData%array3b_error(:,start_valid:stop_valid)
       outputData%array4_error = outData%array4_error(:,start_valid:stop_valid)
       outputData%array5_error = outData%array5_error(:,start_valid:stop_valid)
    endif

    if( outData%newCalibration_There )then
       outputData%newCalibration_There = .TRUE.
       outputData%new_calib3 = outData%new_calib3(:,start_valid:stop_valid)
       outputData%new_calib4 = outData%new_calib4(:,start_valid:stop_valid)
       outputData%new_calib5 = outData%new_calib5(:,start_valid:stop_valid)
       outputData%smoothPrt1 = outData%smoothPrt1(start_valid:stop_valid)
       outputData%smoothPrt2 = outData%smoothPrt2(start_valid:stop_valid)
       outputData%smoothPrt3 = outData%smoothPrt3(start_valid:stop_valid)
       outputData%smoothPrt4 = outData%smoothPrt4(start_valid:stop_valid)
       outputData%smoothPrt1Cnts = outData%smoothPrt1Cnts(start_valid:stop_valid)
       outputData%smoothPrt2Cnts = outData%smoothPrt2Cnts(start_valid:stop_valid)
       outputData%smoothPrt3Cnts = outData%smoothPrt3Cnts(start_valid:stop_valid)
       outputData%smoothPrt4Cnts = outData%smoothPrt4Cnts(start_valid:stop_valid)
       outputData%smoothPrt = outData%smoothPrt(start_valid:stop_valid)
       outputData%smoothBB3 = outData%smoothBB3(start_valid:stop_valid)
       outputData%smoothBB4 = outData%smoothBB4(start_valid:stop_valid)
       outputData%smoothBB5 = outData%smoothBB5(start_valid:stop_valid)
       outputData%smoothSp3 = outData%smoothSp3(start_valid:stop_valid)
       outputData%smoothSp4 = outData%smoothSp4(start_valid:stop_valid)
       outputData%smoothSp5 = outData%smoothSp5(start_valid:stop_valid)
       outputData%Interpolated = outData%Interpolated(start_valid:stop_valid)
       outputData%solar_contamination_3B = &
            outData%solar_contamination_3B(start_valid:stop_valid)
       outputData%solar_contamination_4 = &
            outData%solar_contamination_4(start_valid:stop_valid)
       outputData%solar_contamination_5 = &
            outData%solar_contamination_5(start_valid:stop_valid)
       outputData%moon_contamination = &
            outData%moon_contamination(start_valid:stop_valid)
       outputData%new_array1 = outData%new_array1(:,start_valid:stop_valid)
       outputData%new_array2 = outData%new_array2(:,start_valid:stop_valid)
       outputData%new_array3A = outData%new_array3A(:,start_valid:stop_valid)
       outputData%new_array3B = outData%new_array3B(:,start_valid:stop_valid)
       outputData%new_array4 = outData%new_array4(:,start_valid:stop_valid)
       outputData%new_array5 = outData%new_array5(:,start_valid:stop_valid)
       outputData%new_array1_error = &
            outData%new_array1_error(:,start_valid:stop_valid)
       outputData%new_array2_error = &
            outData%new_array2_error(:,start_valid:stop_valid)
       outputData%new_array3A_error = &
            outData%new_array3A_error(:,start_valid:stop_valid)
       outputData%new_array3B_error = &
            outData%new_array3B_error(:,start_valid:stop_valid)
       outputData%new_array4_error = &
            outData%new_array4_error(:,start_valid:stop_valid)
       outputData%new_array5_error = &
            outData%new_array5_error(:,start_valid:stop_valid)
    endif
    outputData%dataFilled = .TRUE.

  END SUBROUTINE Copy_OutData

  SUBROUTINE Reallocate_Arrays_Real2D( inArray, numberToAdd, numberOfLines )

    REAL, ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    REAL, DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = Dims(2)+numberToAdd

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Real2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Real2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1),1:Dims(2)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(2)

  END SUBROUTINE Reallocate_Arrays_Real2D

  SUBROUTINE Reallocate_Final_Arrays_Real2D( inArray, numberOfLines )

    REAL, ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    REAL, DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = numberOfLines

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Real2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Real2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1),1:outDims(2))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Real2D

  SUBROUTINE Reallocate_Arrays_Int_2D( inArray, numberToAdd, numberOfLines )

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = Dims(2)+numberToAdd

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1),1:Dims(2)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(2)

  END SUBROUTINE Reallocate_Arrays_Int_2D

  SUBROUTINE Reallocate_Final_Arrays_Int_2D( inArray, numberOfLines )

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = numberOfLines

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1),1:outDims(2))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Int_2D

  SUBROUTINE Reallocate_Arrays_Int2_2D( inArray, numberToAdd, numberOfLines )

    INTEGER(GbcsInt2), ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    INTEGER(GbcsInt2), DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = Dims(2)+numberToAdd

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1),1:Dims(2)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(2)

  END SUBROUTINE Reallocate_Arrays_Int2_2D

  SUBROUTINE Reallocate_Final_Arrays_Int2_2D( inArray, numberOfLines )

    INTEGER(GbcsInt2), ALLOCATABLE, DIMENSION(:,:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    INTEGER(GbcsInt2), DIMENSION(:,:), ALLOCATABLE :: outArray
    INTEGER :: Dims(2)
    INTEGER :: outDims(2)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims = Dims
    outDims(2) = numberOfLines

    ALLOCATE(outArray(Dims(1),Dims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1),outDims(2)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int_2D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1),1:outDims(2))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Int2_2D

  SUBROUTINE Reallocate_Arrays_Real1D( inArray, numberToAdd, numberOfLines )

    REAL, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    REAL, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = Dims(1)+numberToAdd

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Real1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Real1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(1)

  END SUBROUTINE Reallocate_Arrays_Real1D

  SUBROUTINE Reallocate_Final_Arrays_Real1D( inArray, numberOfLines )

    REAL, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    REAL, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = numberOfLines

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Real1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Real1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Real1D

  SUBROUTINE Reallocate_Arrays_Int1D( inArray, numberToAdd, numberOfLines )

    INTEGER, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    INTEGER, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = Dims(1)+numberToAdd

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Int1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(1)

  END SUBROUTINE Reallocate_Arrays_Int1D

  SUBROUTINE Reallocate_Final_Arrays_Int1D( inArray, numberOfLines )

    INTEGER, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    INTEGER, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = numberOfLines

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Int1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Int1D

  SUBROUTINE Reallocate_Arrays_Logical1D( inArray, numberToAdd, numberOfLines )

    LOGICAL, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    LOGICAL, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = Dims(1)+numberToAdd

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Logical1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Logical1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(1)

  END SUBROUTINE Reallocate_Arrays_Logical1D

  SUBROUTINE Reallocate_Final_Arrays_Logical1D( inArray, numberOfLines )

    LOGICAL, ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    LOGICAL, DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = numberOfLines

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Logical1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Logical1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Logical1D

  SUBROUTINE Reallocate_Arrays_Double1D( inArray, numberToAdd, numberOfLines )

    REAL(GbcsDble), ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberToAdd
    INTEGER, INTENT(OUT) :: numberOfLines

    REAL(GbcsDble), DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = Dims(1)+numberToAdd

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Double1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray','Reallocate_Arrays_Double1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray(1:Dims(1)) = outArray
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)
    numberOfLines = outDims(1)

  END SUBROUTINE Reallocate_Arrays_Double1D

  SUBROUTINE Reallocate_Final_Arrays_Double1D( inArray, numberOfLines )

    REAL(GbcsDble), ALLOCATABLE, DIMENSION(:) :: inArray
    INTEGER, INTENT(IN) :: numberOfLines

    REAL(GbcsDble), DIMENSION(:), ALLOCATABLE :: outArray
    INTEGER :: Dims(1)
    INTEGER :: outDims(1)
    INTEGER :: STAT

    Dims = SHAPE(inArray)
    outDims(1) = numberOfLines

    ALLOCATE(outArray(Dims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Double1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    outArray = inArray
    IF(ALLOCATED(inArray))DEALLOCATE(inArray)
    ALLOCATE(inArray(outDims(1)),STAT=STAT)
    CALL Check_IOS(STAT,'Allocating outArray',&
         'Reallocate_Final_Arrays_Double1D',&
         'NOAA_LoadAVHRRLevel1B.f90',.FALSE.,' ')

    inArray = outArray(1:outDims(1))
    IF(ALLOCATED(outArray))DEALLOCATE(outArray)

  END SUBROUTINE Reallocate_Final_Arrays_Double1D

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


END MODULE NOAA_LoadAVHRRLevel1B
