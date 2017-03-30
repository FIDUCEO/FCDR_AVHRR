MODULE type

  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil

implicit none

PUBLIC
!########################################################
 ! AVHRR data type
  INTEGER, PUBLIC,PARAMETER :: AVHRR_GAC = 1
  INTEGER,  PUBLIC,PARAMETER :: AVHRR_LAC = 2
  INTEGER,  PUBLIC,PARAMETER :: AVHRR_HRPT = 3
  INTEGER,  PUBLIC,PARAMETER :: AVHRR_FRAC = 4

!########################################################
  ! Define Coefficients Structure
  INTEGER,  PUBLIC,PARAMETER :: MAX_NUC_TEMP = 4
  INTEGER, PARAMETER :: NFILTERS = 3
  INTEGER, PARAMETER :: NCALIB_COEF = 3
  INTEGER, PARAMETER :: NPRT = 4
  INTEGER, PARAMETER :: NPRT_IND = 3
  INTEGER, PARAMETER :: NTEMP_COEFS = 6
  INTEGER, PARAMETER :: NBBODY_SAMPLES = 10
  INTEGER, PARAMETER :: MEMORY_ALLOC_STEP = 20000
  INTEGER, PARAMETER :: NPIXEL_PRT_SMOOTH = 100
  INTEGER, PARAMETER :: NPIXEL_CNTS_SMOOTH = 27
  INTEGER, PARAMETER :: NLOOKUP_TABLE_SCENE = 14
  INTEGER, PARAMETER :: NLOOKUP_TABLE_TARGET = 4

  INTEGER, PARAMETER :: scanLineStep = 20 ! approx 65 km
  INTEGER, PARAMETER :: scanElemStep = 3 ! approx 65 km at scan edge
  ! always skip line 1 since it often has location errors (though
  ! many of these were perhaps from the 'two logical records per
  ! physical record': reduced bad tie-point orbits from 4% to 2% (for M2)
  ! often have location errors at end of orbit (N18)
  INTEGER, PARAMETER :: scanLineStart = 2
  INTEGER, PARAMETER :: scanElemStart = 1
  ! the maximum timing correction magnitude is 1.75 s in the data that is in
  ! clavr-x (which is the same as the data from UMiami) RESTRICTED to the
  ! dates when there is AVHRR data available.  So expect the shift to be
  ! less than 1.75 * 6 * 1.1 = time/s * lines/s * km/line = 11.50 km, say 20km
  ! as a check.
  REAL, PARAMETER    :: maxTimingErrorDistance = 20.0

!########################################################  
  ! Constants
  REAL, PARAMETER :: C1 = 1.1910427E-5
  REAL, PARAMETER :: C2 = 1.4387752
  REAL, PARAMETER :: eta_ict = 0.985140
  double precision, parameter :: M_PI = 3.14159265358979323846d0 

  INTEGER, PARAMETER, DIMENSION(13) :: dayNumber = &
       (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
  INTEGER, PARAMETER, DIMENSION(13) :: dayNumberLeap = &
       (/0,31,60,91,121,152,182,213,244,274,305,335,366/)


!######################################################## 
  TYPE AVHRR_Instrument_Coefs
     ! nuC for temperature -> radiance conversion 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: nuC
     ! A for temperature -> radiance conversion 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: Aval
     ! B for temperature -> radiance conversion 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: Bval
     ! F for reflectance -> radiance conversion 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: F
     ! w for reflectance -> radiance conversion 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: w
     ! Older style radiance->temperarure coefficients (pre KLM)
     ! nuC for temperature -> radiance conversion 
     ! Number of entries
     INTEGER(GbcsInt1) :: nuC_numTemps
     ! minimum temperature range
     REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP) :: nuC_minT
     ! maximum temperature range
     REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP) :: nuC_maxT
     ! nuC for temperature -> radiance conversion (old way)
     REAL(GbcsReal), DIMENSION(MAX_NUC_TEMP,NFILTERS) :: nuC_Array
     ! Space radiance for counts -> radiance 
     REAL(GbcsReal), DIMENSION(NFILTERS) :: Nspace
     ! detector non-linear coeficients filter 3B
     REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs3
     ! detector non-linear coeficients filter 3B
     REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs4
     ! detector non-linear coeficients filter 3B
     REAL(GbcsReal), DIMENSION(NCALIB_COEF) :: nonLinearCoefs5
     ! Counts -> temperature coefs (0) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs1
     ! Counts -> temperature coefs (1) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs2
     ! Counts -> temperature coefs (2) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs3
     ! Counts -> temperature coefs (3) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs4
     ! Counts -> temperature coefs (4) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs5
     ! Counts -> temperature coefs (4) for NPRTs 
     REAL(GbcsReal), DIMENSION(NPRT) :: prtTempCoefs6
     ! Counts -> temperature patch 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: patchCoef
     ! Counts -> temperature extended patch 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: patchCoefExt
     ! Counts -> temperature (Radiator) */
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs1
     ! Counts -> temperature (Electronics) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs2
     ! Counts -> temperature (Cooler) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs3
     ! Counts -> temperature (Baseplate) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs4
     ! Counts -> temperature (Motor) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs5
     ! Counts -> temperature (A/D convertor) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs6
     ! Counts -> temperature (Patch) 
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: temperatureCoefs7
     ! Counts->Volts->Current (Motor)
     REAL(GbcsReal), DIMENSION(NTEMP_COEFS) :: motorCurrentCoefs
     ! Lookup table for non-linear correction (old method)
     INTEGER :: number_Target_Temps
     INTEGER :: number_Scene_Temps
     REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE) :: Scene_Temp
     REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_TARGET) :: Target_Temp
     REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE,NLOOKUP_TABLE_TARGET) :: &
          Correction_Factor4
     REAL(GbcsReal), DIMENSION(NLOOKUP_TABLE_SCENE,NLOOKUP_TABLE_TARGET) :: &
          Correction_Factor5
  END TYPE AVHRR_Instrument_Coefs 
 
!########################################################  
  TYPE AVHRR_thresholds_coarse
     REAL, DIMENSION(3) :: Low_sp,High_sp
     REAL, DIMENSION(3) :: Low_bb,High_bb
     REAL, DIMENSION(4) :: Low_prt,High_prt
  END TYPE AVHRR_thresholds_coarse

!########################################################
  TYPE AVHRR_Data
     LOGICAL :: isGAC
     LOGICAL :: dataFilled
     INTEGER :: AVHRR_No
     INTEGER :: nelem
     INTEGER :: arraySize
     LOGICAL :: filter3a
     INTEGER :: start_valid
     INTEGER :: stop_valid
!*************************
!modif du 27 aout 2015
!et du 17 septembre 2015
!et du 23 novembre 2015 
!et du 16 fevrier 2016
     REAL ::  mean_patchExtended_orbite
     REAL ::  mean_patch_orbite
     REAL :: mean_space_counts_orbite_c3 
     REAL :: mean_bb_counts_orbite_c3
     REAL :: mean_rict_orbite_c3

     REAL :: mean_space_counts_orbite_c4 
     REAL :: mean_bb_counts_orbite_c4
     REAL :: mean_rict_orbite_c4
    
     REAL :: mean_space_counts_orbite_c5 
     REAL :: mean_bb_counts_orbite_c5
     REAL :: mean_rict_orbite_c5

     REAL :: mean_prt_orbite

     REAL, DIMENSION(3) :: mean_calib3_orbite
     REAL, DIMENSION(3) :: mean_calib4_orbite
     REAL, DIMENSION(3) :: mean_calib5_orbite

!*************************
     REAL, DIMENSION(3)::Nspace
     INTEGER, ALLOCATABLE, DIMENSION(:) :: scanLineNumber
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: badTime
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: badNavigation
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: badCalibration
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: transition3A3B
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Lon
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Lat
     REAL, ALLOCATABLE, DIMENSION(:,:) :: satZA
     REAL, ALLOCATABLE, DIMENSION(:,:) :: solZA
    REAL, ALLOCATABLE, DIMENSION(:,:) :: relAz
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts1
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts2
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts3
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts4
     REAL, ALLOCATABLE, DIMENSION(:,:) :: Counts5
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array1
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array2
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array3A
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array3B
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array4
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array5
     LOGICAL :: walton_there
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array3B_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array4_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: array5_error
     INTEGER, ALLOCATABLE, DIMENSION(:) :: year
     INTEGER, ALLOCATABLE, DIMENSION(:) :: month
     INTEGER, ALLOCATABLE, DIMENSION(:) :: day
     INTEGER, ALLOCATABLE, DIMENSION(:) :: dayNo
     REAL, ALLOCATABLE, DIMENSION(:) :: hours
     INTEGER, ALLOCATABLE, DIMENSION(:) :: UTC_msecs
     REAL(GbcsDble), ALLOCATABLE, DIMENSION(:) :: time
!****************************
! modif du 16 nov
     INTEGER, ALLOCATABLE, DIMENSION(:) :: prtNumber
!***********************
     REAL, ALLOCATABLE, DIMENSION(:) :: prt1
     REAL, ALLOCATABLE, DIMENSION(:) :: prt2
     REAL, ALLOCATABLE, DIMENSION(:) :: prt3
     REAL, ALLOCATABLE, DIMENSION(:) :: prt4
     REAL, ALLOCATABLE, DIMENSION(:) :: prt1Counts
     REAL, ALLOCATABLE, DIMENSION(:) :: prt2Counts
     REAL, ALLOCATABLE, DIMENSION(:) :: prt3Counts
     REAL, ALLOCATABLE, DIMENSION(:) :: prt4Counts
     REAL, ALLOCATABLE, DIMENSION(:,:) :: prt1CountsAll
     REAL, ALLOCATABLE, DIMENSION(:,:) :: prt2CountsAll
     REAL, ALLOCATABLE, DIMENSION(:,:) :: prt3CountsAll
     REAL, ALLOCATABLE, DIMENSION(:,:) :: prt4CountsAll
!*************************
!modif du 14 oct 2015

     REAL, ALLOCATABLE, DIMENSION(:) :: ramp_c3
     REAL, ALLOCATABLE, DIMENSION(:) :: ramp_c4
     REAL, ALLOCATABLE, DIMENSION(:) :: ramp_c5
     REAL, ALLOCATABLE, DIMENSION(:) :: gain_c3
     REAL, ALLOCATABLE, DIMENSION(:) :: gain_c4
     REAL, ALLOCATABLE, DIMENSION(:) :: gain_c5
!************************
!*************************
!modif du 12 nov 2015

     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: prtcounts
    
!************************
!*************************
!modif du 20 aout 2015

     REAL, ALLOCATABLE, DIMENSION(:) :: prtmean
     REAL, ALLOCATABLE, DIMENSION(:) :: prtsigma
!*************************
!*************************
!modif du 27 aout 2015

     REAL, ALLOCATABLE, DIMENSION(:) :: Rict_c3
     REAL, ALLOCATABLE, DIMENSION(:) :: Rict_c4
     REAL, ALLOCATABLE, DIMENSION(:) :: Rict_c5
!*************************
     REAL, ALLOCATABLE, DIMENSION(:) :: bb3
     REAL, ALLOCATABLE, DIMENSION(:) :: bb4
     REAL, ALLOCATABLE, DIMENSION(:) :: bb5
     REAL, ALLOCATABLE, DIMENSION(:) :: sp3
     REAL, ALLOCATABLE, DIMENSION(:) :: sp4
     REAL, ALLOCATABLE, DIMENSION(:) :: sp5
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sp3_filtered
     REAL, ALLOCATABLE, DIMENSION(:) :: sp4_filtered
     REAL, ALLOCATABLE, DIMENSION(:) :: sp5_filtered
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter3
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter4
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bbodyFilter5
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter1
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter2
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter3a
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter3
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter4
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spaceFilter5
     REAL, ALLOCATABLE, DIMENSION(:) :: patch
     REAL, ALLOCATABLE, DIMENSION(:) :: patchExtended
     REAL, ALLOCATABLE, DIMENSION(:) :: Radiator
     REAL, ALLOCATABLE, DIMENSION(:) :: Cooler
     REAL, ALLOCATABLE, DIMENSION(:) :: a_d_conv
     REAL, ALLOCATABLE, DIMENSION(:) :: motor
     REAL, ALLOCATABLE, DIMENSION(:) :: electronics
     REAL, ALLOCATABLE, DIMENSION(:) :: baseplate
     REAL, ALLOCATABLE, DIMENSION(:) :: motorCurrent
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib1
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib1_2
     REAL, ALLOCATABLE, DIMENSION(:) :: calib1_intercept
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib2
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib2_2
     REAL, ALLOCATABLE, DIMENSION(:) :: calib2_intercept
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3A
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3A_2
     REAL, ALLOCATABLE, DIMENSION(:) :: calib3A_intercept
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib3
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib4
     REAL, ALLOCATABLE, DIMENSION(:,:) :: calib5
     LOGICAL :: Clavr_There
     INTEGER(GbcsInt2), ALLOCATABLE, DIMENSION(:,:) :: clavr_mask
     LOGICAL :: Clavrx_There
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: clavrx_mask ! HDF5 wants integer to read into, even if the
                                                         ! file is byte
     INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: clavrx_prb
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_3B
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_4
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: orig_solar_contamination_5
     LOGICAL :: newCalibration_There
     REAL :: orbital_temperature
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib3
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib4
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_calib5
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt1
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt2
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt3
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt4
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt1Cnts
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt2Cnts
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt3Cnts
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt4Cnts
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothPrt
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB3
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB4
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothBB5
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp3
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp4
     REAL, ALLOCATABLE, DIMENSION(:) :: smoothSp5
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: Interpolated
     LOGICAL :: solar_contamination_failure
     LOGICAl, ALLOCATABLE, DIMENSION(:) :: solar_contamination_3B
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: solar_contamination_4
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: solar_contamination_5
     LOGICAL, ALLOCATABLE, DIMENSION(:) :: moon_contamination
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array1
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array2
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3A
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3B
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array4
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array5
     REAL, DIMENSION(8) :: new_cal_coefs3
     REAL, DIMENSION(8) :: new_cal_coefs4
     REAL, DIMENSION(8) :: new_cal_coefs5
     REAL, DIMENSION(3) :: gain_stdev
!***************************
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array1_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array2_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3A_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array3B_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array4_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: new_array5_error
     REAL, ALLOCATABLE, DIMENSION(:,:) :: noise_cnts
    
     REAL :: earthshine_eta 
     REAL :: poly_coefs3(3)
     REAL :: poly_coefs4(3)
     REAL :: poly_coefs5(3)
     REAL :: gain_maxdev(3)
     LOGICAL :: satelliteAlt_There
     REAL, ALLOCATABLE, DIMENSION(:) :: satelliteAltitude
!************************   
!MODIF 5 juillet 2016
!*************************
     REAL, ALLOCATABLE,DIMENSION(:)::tstar3,tstar4,tstar5
     REAL, ALLOCATABLE,DIMENSION(:,:)::btf3,btf4,btf5
     REAL, ALLOCATABLE,DIMENSION(:,:)::nef3,nef4,nef5

     REAL, ALLOCATABLE,DIMENSION(:,:)::ue3,ue4,ue5, &
                                       uce3,uce4,uce5, &
                                       ustr3,ustr4,ustr5, &
                                       ur3, ur4, ur5, &
                                       us3, us4, us5
                                       
     REAL                         :: ucs3,ucs4,ucs5, &
                                     ucict3,ucict4,ucict5
     REAL, ALLOCATABLE,DIMENSION(:):: utict, &
                                      utstar3, utstar4, utstar5,& 
                                      urict3,urict4,urict5
     
     REAL, ALLOCATABLE,DIMENSION(:,:):: ve3,ve4,ve5, &
                                        vstr3,vstr4,vstr5

     REAL, ALLOCATABLE,DIMENSION(:,:):: dre_over_dce3,dre_over_dce4,dre_over_dce5, &
                                        dre_over_dcs3, dre_over_dcs4,dre_over_dcs5, &
                                        dre_over_dcict3,dre_over_dcict4,dre_over_dcict5, &
                                        dre_over_db03, dre_over_db04, dre_over_db05, &
                                        dre_over_da13, dre_over_da14, dre_over_da15, &
                                        dre_over_da23, dre_over_da24, dre_over_da25, &
                                        dre_over_db13, dre_over_db14, dre_over_db15, &
                                        dre_over_dtinstr3, dre_over_dtinstr4, dre_over_dtinstr5, &
                                        dre_over_drict3, dre_over_drict4, dre_over_drict5, &
                                        dre_over_depsilon 

     REAL, ALLOCATABLE,DIMENSION(:):: dtstar_over_daval3,dtstar_over_daval4,dtstar_over_daval5, &
                                        dtstar_over_dbval3,dtstar_over_dbval4,dtstar_over_dbval5,&
                                        dtstar_over_dnuc3,dtstar_over_dnuc4,dtstar_over_dnuc5

     REAL, ALLOCATABLE,DIMENSION(:)::   drict_over_dnuc3,  drict_over_dnuc4, drict_over_dnuc5, &
                                        drict_over_dtstar3,drict_over_dtstar4, drict_over_dtstar5
                                       
 
    
  END TYPE AVHRR_Data

!########################################################
  TYPE AVHRR_Bad_Data
     LOGICAL :: allow_bad_top
     LOGICAL :: allow_bad_nav
     LOGICAL :: use_bad_nav
     LOGICAL :: allow_bad_time
     LOGICAL :: allow_bad_calib
  END TYPE AVHRR_Bad_Data

!########################################################
  TYPE AVHRR_Radiative_Coefs
     REAL :: nuc(NFILTERS)
     REAL :: aVal(NFILTERS)
     REAL :: bVal(NFILTERS)
  END TYPE AVHRR_Radiative_Coefs

!########################################################
  TYPE AVHRR_Land_Mask
     INTEGER :: nelem
     INTEGER :: arraySize
     INTEGER(GbcsInt1), ALLOCATABLE, DIMENSION(:,:) :: mask
  END TYPE AVHRR_Land_Mask

!########################################################
 ! PRT rotations
  INTEGER, PARAMETER :: prtRotate_GAC(4) = (/3,1,4,2/)
  INTEGER, PARAMETER :: prtRotate_LAC(4) = (/1,2,3,4/)

!########################################################
  ! Some global variables
  LOGICAL, SAVE :: AVHRR_Setup = .FALSE.
  INTEGER :: ndata 
  INTEGER :: ndata_counts
  INTEGER :: prtRotate(4)
  INTEGER :: printOnce(10)
  INTEGER, PARAMETER :: NDATA_GAC = 409
  INTEGER, PARAMETER :: NDATA_LAC = 2048
  INTEGER, PARAMETER :: NDATA_MAX = 2048
  INTEGER :: prtNumber
  REAL :: storedPrtTemp(4)
  REAL :: prtCountsStore(4)
  REAL :: prtCountsPrevStore(4)
  REAL :: storedPrtCounts(4)
  REAL :: prtCountsStoreAll(3,4)
  REAL :: storedPrtCountsAll(3,4)
  LOGICAL :: scanAngleSetup = .FALSE.
 
  REAL :: scanAngle(NDATA_MAX)
  REAL :: scanAngleIn(51)
  INTEGER :: time_yearstart
  TYPE(AVHRR_Data), TARGET, SAVE :: outDataStore
  TYPE(DateTime), SAVE :: start_time, end_time
  LOGICAL :: badDay
  LOGICAL :: validSolZaDecimal
  LOGICAL :: validClockDriftInfo
  LOGICAL :: useL1BClock
  LOGICAL :: validSBBC
  REAL, ALLOCATABLE, DIMENSION(:,:) :: clavrLonAnchor, clavrLatAnchor
  
  LOGICAL :: Little_Endian = .TRUE.
!########################################################
! Noise
  REAL                             :: prt_accuracy=0.1, prt_noise=0.015, &
                                      uaval=0, ubval=0, unuc=0
!#########################
!COUNTS TO RAD EARTH COEF 
! Coef calcules par l'harmonisation (version marine novembre 2016)
!########################
!1:bias term 
!2:ICT emissivity correction
!3:non-linear term 
!4:instrument temperature slope term
REAL , PARAMETER :: &
counts_to_rad_earth_coef_m02(3,4)=reshape((/1.11362134,0.00179079409535,0.,-0.00388403, &
6.88329791E+01, -6.01531058E-03,1.67568E-5,-0.23486149, &
5.09659390E+01, 7.05761507E-04, 1.063571E-5, -0.17470683/),(/3,4/)) 

REAL , PARAMETER :: &
!counts_to_rad_earth_coef_n19(3,4)=reshape((/1.15609851e+00, 6.10582553e-01, 1.58704597e-05, -4.01954808e-03, &
!5.02984886e+02, 5.45186283e-02, 4.78411936e-05, -1.74207901e+00, &
!5.84787339e+02, 4.93335608e-02, 4.90686745e-05, -2.02904891e+00/),(/3,4/))

counts_to_rad_earth_coef_n19(3,4)=reshape((/ &
-2.15239130E+00, 5.50695297E+00, 1.55046118E-04, 7.48275207E-03, &
2.86135070E+02, 5.48844112E-02, 4.46246537E-05, -9.88694358E-01, &
4.45629521E+02, 5.11796041E-02, 4.82689120E-05, -1.54555024E+00/),(/3,4/))

REAL, PARAMETER :: &
covariance_matrix_ch3(4,4)=reshape((/&
3.29033234e-11, 1.68572212e-11, 4.24054635e-16, -1.14412813e-13, &
1.68572212e-11, 5.25821800e-11, 1.35145044e-15, -5.86007060e-14, &
4.24054635e-16, 1.35145044e-15, 3.48812895e-20, -1.47407944e-18, &
-1.14412813e-13, -5.86007060e-14, -1.47407944e-18, 3.97841046e-16/),(/4,4/))

REAL, PARAMETER :: &
covariance_matrix_ch4(4,4)=reshape((/&
3.11354272e-02, 4.97622793e-08, -2.92995252e-10, -1.08284746e-04, &
4.97622793e-08, 1.64634254e-11, 1.62816551e-14, -1.71954874e-10, &
-2.92995252e-10, 1.62816551e-14, 4.07123831e-17, 1.02368988e-12, &
-1.08284746e-04, -1.71954874e-10, 1.02368988e-12, 3.76600213e-07/),(/4,4/))

REAL, PARAMETER :: &
covariance_matrix_ch5(4,4)=reshape((/&
4.80426207e-02, 1.62072681e-07, -3.69145208e-10, -1.67091016e-04, &
1.62072681e-07, 1.31727286e-11, 1.15136937e-14, -5.63190203e-10, &
-3.69145208e-10, 1.15136937e-14, 6.00894832e-17, 1.29281208e-12,&
-1.67091016e-04, -5.63190203e-10, 1.29281208e-12, 5.81139990e-07/),(/4,4/))

end module type
