! * Copyright (C) 2017 J.Mittaz University of Reading
! * This code was developed for the EC project Fidelity and Uncertainty in
! * Climate Data Records from Earth Observations (FIDUCEO).
! * Grant Agreement: 638822
! *
! * This program is free software; you can redistribute it and/or modify it
! * under the terms of the GNU General Public License as published by the Free
! * Software Foundation; either version 3 of the License, or (at your option)
! * any later version.
! * This program is distributed in the hope that it will be useful, but WITHOUT
! * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
! * more details.
! *
! * A copy of the GNU General Public License should have been supplied along
! * with this program; if not, see http://www.gnu.org/licenses/
!
!
! This code is a simple rewrite of Marines original code to split out
! the FIDUCEO uncertainties and sensitivities into a different structure
! and then write it out
!
! Only science changes are using the new readers filtering/Allan deviation
! calculations plus adding in a different systematic component to the 
! visible channels
!
! Written by J,Mittaz University of Reading 
! Under the FIDUCEO project (www.fiduceo.eu)
!
! Original 21-07-2017 (v0.1pre)
!
! MT: 11-11-2017: Define temp variables us1,us2,us3a to store structured uncertainties on the reflectance channels
! MT: 11-11-2017: fix problem of value not filling array with 0.03 for u_structured_Ch1
! MT: 11-11-2017: fix problem of value not filling array with 0.05 for u_structured_Ch2
! MT: 11-11-2017: fix problem of value not filling array with 0.05 for u_structured_Ch3a
! MT: 13-11-2017: allocated nsmoothBB3,4,5 and nsmoothSp3,4,5 to AVHRRout data structure in combine_orbits.f90 so that the calculations don't fail 
! MT: 19-12-2017: v0.3pre 
! MT: 09-03-2018: v0.5beta 
!
! Note: Coefs data from CCI are ordered as
!
!    1    = a0   (bias term)
!    2    = a1   (offset to emissivity)
!    3    =      (ICT bias term for IAsI calibration
!    4    = a2   (non-linear term)
!    5    = nuc  (SRF)
!    6    = aval   "
!    7    = bval   "
!    8    = a4   (Tinstrument term)

MODULE fiduceo_uncertainties

  ! From Marines code base
  USE GbcsConstants
  USE GbcsErrorHandler
  USE NOAA_LoadAVHRRLevel1B
  USE NETCDF

  IMPLICIT NONE

  TYPE FIDUCEO_Data
     REAL, ALLOCATABLE :: Rict_c3(:)
     REAL, ALLOCATABLE :: Rict_c4(:)
     REAL, ALLOCATABLE :: Rict_c5(:)
     REAL, ALLOCATABLE :: dtstar_over_dT3(:)
     REAL, ALLOCATABLE :: dtstar_over_dT4(:)
     REAL, ALLOCATABLE :: dtstar_over_dT5(:)
     REAL, ALLOCATABLE :: dtstar_over_daval3(:)
     REAL, ALLOCATABLE :: dtstar_over_daval4(:)
     REAL, ALLOCATABLE :: dtstar_over_daval5(:)
     REAL, ALLOCATABLE :: dtstar_over_dbval3(:)
     REAL, ALLOCATABLE :: dtstar_over_dbval4(:)
     REAL, ALLOCATABLE :: dtstar_over_dbval5(:)
     REAL, ALLOCATABLE :: dtstar_over_dnuc3(:)
     REAL, ALLOCATABLE :: dtstar_over_dnuc4(:)
     REAL, ALLOCATABLE :: dtstar_over_dnuc5(:)
     REAL, ALLOCATABLE :: drict_over_dnuc3(:)
     REAL, ALLOCATABLE :: drict_over_dnuc4(:)
     REAL, ALLOCATABLE :: drict_over_dnuc5(:)
     REAL, ALLOCATABLE :: drict_over_dtstar3(:)
     REAL, ALLOCATABLE :: drict_over_dtstar4(:)
     REAL, ALLOCATABLE :: drict_over_dtstar5(:)
     REAL, ALLOCATABLE :: dre_over_db03(:,:)
     REAL, ALLOCATABLE :: dre_over_db04(:,:)
     REAL, ALLOCATABLE :: dre_over_db05(:,:)
     REAL, ALLOCATABLE :: dre_over_db13(:,:)
     REAL, ALLOCATABLE :: dre_over_db14(:,:)
     REAL, ALLOCATABLE :: dre_over_db15(:,:)
     REAL, ALLOCATABLE :: dre_over_da13(:,:)
     REAL, ALLOCATABLE :: dre_over_da14(:,:)
     REAL, ALLOCATABLE :: dre_over_da15(:,:)
     REAL, ALLOCATABLE :: dre_over_da23(:,:)
     REAL, ALLOCATABLE :: dre_over_da24(:,:)
     REAL, ALLOCATABLE :: dre_over_da25(:,:)
     REAL, ALLOCATABLE :: dre_over_dtinstr3(:,:)
     REAL, ALLOCATABLE :: dre_over_dtinstr4(:,:)
     REAL, ALLOCATABLE :: dre_over_dtinstr5(:,:)
     REAL, ALLOCATABLE :: dre_over_dce3(:,:)
     REAL, ALLOCATABLE :: dre_over_dce4(:,:)
     REAL, ALLOCATABLE :: dre_over_dce5(:,:)
     REAL, ALLOCATABLE :: dre_over_dcict3(:,:)
     REAL, ALLOCATABLE :: dre_over_dcict4(:,:)
     REAL, ALLOCATABLE :: dre_over_dcict5(:,:)
     REAL, ALLOCATABLE :: dre_over_drict3(:,:)
     REAL, ALLOCATABLE :: dre_over_drict4(:,:)
     REAL, ALLOCATABLE :: dre_over_drict5(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs3(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs4(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs5(:,:)
     REAL, ALLOCATABLE :: ue3(:,:)
     REAL, ALLOCATABLE :: ue4(:,:)
     REAL, ALLOCATABLE :: ue5(:,:)
     REAL, ALLOCATABLE :: ur3(:,:)
     REAL, ALLOCATABLE :: ur4(:,:)
     REAL, ALLOCATABLE :: ur5(:,:)
     REAL, ALLOCATABLE :: us3(:,:)
     REAL, ALLOCATABLE :: us4(:,:)
     REAL, ALLOCATABLE :: us5(:,:)
     REAL, ALLOCATABLE :: btf3(:,:)
     REAL, ALLOCATABLE :: btf4(:,:)
     REAL, ALLOCATABLE :: btf5(:,:)
     REAL, ALLOCATABLE :: uce3(:,:)
     REAL, ALLOCATABLE :: uce4(:,:)
     REAL, ALLOCATABLE :: uce5(:,:)
     REAL, ALLOCATABLE :: urict3_r(:)
     REAL, ALLOCATABLE :: urict4_r(:)
     REAL, ALLOCATABLE :: urict5_r(:)
     REAL, ALLOCATABLE :: urict3_s(:)
     REAL, ALLOCATABLE :: urict4_s(:)
     REAL, ALLOCATABLE :: urict5_s(:)
     REAL, ALLOCATABLE :: ucict3(:)
     REAL, ALLOCATABLE :: ucict4(:)
     REAL, ALLOCATABLE :: ucict5(:)
     REAL, ALLOCATABLE :: ucs3(:)
     REAL, ALLOCATABLE :: ucs4(:)
     REAL, ALLOCATABLE :: ucs5(:)
     
     INTEGER, ALLOCATABLE :: flag_no_detection(:,:)
     INTEGER(GbcsInt1), ALLOCATABLE :: quality_channel_bitmask(:,:)
     INTEGER(GbcsInt1), ALLOCATABLE :: quality_scanline_bitmask(:)
  END TYPE FIDUCEO_Data

  ! Constants
  REAL, PARAMETER :: c1 = 1.1910427E-5
  REAL, PARAMETER :: c2 = 1.4387752
  REAL, PARAMETER :: eta_ict = 0.985140

  !
  ! This is where the FIDUCEO software version number is defined
  !
  ! MT: 19-12-2017: v0.3pre
  ! MT: 09-03-2018: v0.5beta
  CHARACTER(LEN=6) :: software_version = '0.5pre'

! MT: 11-11-2017: Define temp variables to store structured uncertainties on the reflectance channels
  REAL, ALLOCATABLE :: us1(:,:)
  REAL, ALLOCATABLE :: us2(:,:)
  REAL, ALLOCATABLE :: us3a(:,:)
 
  PRIVATE
  PUBLIC :: FIDUCEO_Data
  PUBLIC :: Add_FIDUCEO_Uncert

CONTAINS

  SUBROUTINE Add_FIDUCEO_Uncert(AVHRR,uuid_in,filename_nc,use_iasi_calibration)

    TYPE(AVHRR_Data), INTENT(INOUT) :: AVHRR
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    CHARACTER(LEN=*), INTENT(IN) :: filename_nc
    LOGICAL, OPTIONAL :: use_iasi_calibration

    ! Local variables
    INTEGER :: I
    INTEGER :: ncoefs
    REAL :: coefs1(8,2)
    REAL :: coefs2(8,2)
    REAL :: coefs3(8,2)
    REAL :: coefs_frac1
    REAL :: coefs_frac2
    REAL :: cal_coef_overlap
    CHARACTER(LEN=10) :: height
    CHARACTER(LEN=512) :: command_fcdr
    CHARACTER(LEN=256) :: temp_file
    LOGICAL :: use_iasi_cal
    LOGICAL :: twelve_micron_there
    TYPE(FIDUCEO_Data) :: FCDR

    IF( .not. AVHRR%valid_data_there )THEN
       RETURN
    ENDIF
    IF( PRESENT(use_iasi_calibration) )THEN
       use_iasi_cal = use_iasi_calibration
    ELSE
       use_iasi_cal = .FALSE.
    ENDIF

    call Get_Calib_Coefficients_CCI( AVHRR%time(AVHRR%start_valid), &
         AVHRR%AVHRR_No,&
         coefs1, coefs2, coefs3, ncoefs, coefs_frac1, coefs_frac2, &
         cal_coef_overlap, use_iasi_cal, twelve_micron_there )

    !
    ! Do FIDUCEO Uncertainties
    !
    DO I=1,AVHRR%arraySize
       CALL Get_Sensitivities(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)
       CALL calculate_urict(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)
       CALL radiance_uncertainties(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)
    END DO
    !
    ! Get Quality flags
    !
    CALL Get_Quality_Flags(AVHRR,FCDR)
    !
    ! write to NetCDF
    !
    temp_file=TRIM(uuid_in)//'.nc'
    CALL Write_Temp_NETCDF(temp_file,AVHRR,FCDR,twelve_micron_there)
!    CALL Rescale(AVHRR,FCDR)
    IF( 'None' .eq. filename_nc )THEN
       command_fcdr ='python2.7 write_easy_fcdr_from_netcdf.py '//TRIM(temp_file)
    ELSE
       command_fcdr ='python2.7 write_easy_fcdr_from_netcdf.py '//TRIM(temp_file)//' '//TRIM(filename_nc)
    ENDIF
    call SYSTEM(TRIM(command_fcdr))
       command_fcdr = 'rm -f '//TRIM(temp_file) !MT: 05-11-2017: comment to keep temp netcdf files
    call SYSTEM(TRIM(command_fcdr))
    print*, "remplissage"
!   Which is French for "filling"
!    call fill_netcdf(filename_nc,AVHRR,FCDR)

  END SUBROUTINE Add_FIDUCEO_Uncert

  !
  ! Work out quality flags from data in AVHRR_Data structure
  !
  SUBROUTINE Get_Quality_Flags(AVHRR,FCDR)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(FIDUCEO_Data), INTENT(INOUT) :: FCDR

    ! Local variables
    INTEGER :: I,J
    INTEGER :: STAT

    !
    ! Allocate and set quality bitmasks
    !
    ALLOCATE(FCDR%quality_channel_bitmask(6,AVHRR%arraySize),&
         FCDR%quality_scanline_bitmask(AVHRR%arraySize),&
         STAT=STAT)
    FCDR%quality_channel_bitmask = 0
    FCDR%quality_scanline_bitmask = 0

    !
    ! Scanline level flags
    !
    DO I=1,AVHRR%arraySize
       IF( AVHRR%badTop(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),0)
       ENDIF
       IF( AVHRR%badTime(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),1)
       ENDIF
       IF( AVHRR%badNavigation(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),2)
       ENDIF
       IF( AVHRR%badCalibration(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),3)
       ENDIF
       !
       ! Channel 3A present
       !
       IF( ANY(AVHRR%array3A(:,I) .ge. 0) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),4)
       ENDIF       
       !
       ! Solar contamination failure
       !
       IF( AVHRR%solar_contamination_failure(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),5)
       ENDIF
       !
       ! Solar contamination
       !
       IF( AVHRR%solar_contamination_3B(I) )THEN
          FCDR%quality_scanline_bitmask(I) = &
               IBSET(FCDR%quality_scanline_bitmask(I),6)
       ENDIF
    END DO

    !
    ! Channel level quality flags
    !
    DO I=1,AVHRR%arraySize
       IF( ALL(AVHRR%array1(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(1,I) = &
               IBSET(FCDR%quality_channel_bitmask(1,I),0)
       ENDIF
       IF( ALL(AVHRR%array2(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(2,I) = &
               IBSET(FCDR%quality_channel_bitmask(2,I),0)
       ENDIF
       IF( ALL(AVHRR%array3A(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(3,I) = &
               IBSET(FCDR%quality_channel_bitmask(3,I),0)
       ENDIF
       IF( ALL(AVHRR%array3B(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(4,I) = &
               IBSET(FCDR%quality_channel_bitmask(4,I),0)
       ENDIF
       IF( ALL(AVHRR%array4(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(5,I) = &
               IBSET(FCDR%quality_channel_bitmask(5,I),0)
       ENDIF
       IF( ALL(AVHRR%array5(:,I) .lt. 0) )THEN
          FCDR%quality_channel_bitmask(6,I) = &
               IBSET(FCDR%quality_channel_bitmask(6,I),0)
       ENDIF
       DO J=4,6
          IF( 1 .eq. FCDR%flag_no_detection(J-3,I) )THEN
             FCDR%quality_channel_bitmask(J,I) = &
                  IBSET(FCDR%quality_channel_bitmask(J,I),1)
          ENDIF
       END DO
    END DO

  END SUBROUTINE Get_Quality_Flags

  !
  ! Write a tempory netcdf file for python to then convert
  ! Different from Marine's original version
  !
  SUBROUTINE Write_Temp_NETCDF(temp_file,AVHRR,FCDR,twelve_micron_there)

    CHARACTER(LEN=*), INTENT(IN) :: temp_file
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(FIDUCEO_Data), INTENT(IN) :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there

    ! Local variables
    INTEGER :: ncid
    INTEGER :: latitude_varid
    INTEGER :: longitude_varid
    INTEGER :: time_varid
    INTEGER :: year_varid
    INTEGER :: month_varid
    INTEGER :: day_varid
    INTEGER :: hours_varid
    INTEGER :: satza_varid
    INTEGER :: solza_varid
    INTEGER :: relaz_varid
    INTEGER :: ch1_varid
    INTEGER :: ch2_varid
    INTEGER :: ch3a_varid
    INTEGER :: ch3b_varid
    INTEGER :: ch4_varid
    INTEGER :: ch5_varid
    INTEGER :: ch1_random_varid
    INTEGER :: ch2_random_varid
    INTEGER :: ch3a_random_varid
    INTEGER :: ch3b_random_varid
    INTEGER :: ch4_random_varid
    INTEGER :: ch5_random_varid
    INTEGER :: ch1_non_random_varid
    INTEGER :: ch2_non_random_varid
    INTEGER :: ch3a_non_random_varid
    INTEGER :: ch3b_non_random_varid
    INTEGER :: ch4_non_random_varid
    INTEGER :: ch5_non_random_varid
    INTEGER :: scan_qual_varid
    INTEGER :: chan_qual_varid
    INTEGER :: stat

    INTEGER :: dimid_nx
    INTEGER :: dimid_ny
    INTEGER :: dimid_ir
    INTEGER :: dims1(1)
    INTEGER :: dims2(2)

    CHARACTER(LEN=20) :: noaa_string
    
    stat = NF90_CREATE(temp_file,IOR(NF90_HDF5,NF90_CLOBBER),ncid)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nx',AVHRR%nelem,dimid_nx)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'ny',AVHRR%arraySize,dimid_ny)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nir',6,dimid_ir)
    call check(stat)

    dims1(1) = dimid_ny
    dims2(1) = dimid_nx
    dims2(2) = dimid_ny
    
    stat = NF90_DEF_VAR(ncid,'latitude',NF90_FLOAT,dims2,latitude_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, latitude_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, latitude_varid, 0, -1e30)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'longitude',NF90_FLOAT,dims2,longitude_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, longitude_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, longitude_varid, 0, -1e30)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'time',NF90_DOUBLE,dims1,time_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, time_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, time_varid, 0, -1d30)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'year',NF90_INT,dims1,year_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, year_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, year_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'month',NF90_INT,dims1,month_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, month_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, month_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'day',NF90_INT,dims1,day_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, day_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, day_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'hours',NF90_FLOAT,dims1,hours_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, hours_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, hours_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'satza',NF90_FLOAT,dims2,satza_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, satza_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, satza_varid, 0, -1e30)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'solza',NF90_FLOAT,dims2,solza_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, solza_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, solza_varid, 0, -1e30)
    call check(stat)

    IF( ALLOCATED(AVHRR%relaz) )THEN
       stat = NF90_DEF_VAR(ncid,'relaz',NF90_FLOAT,dims2,relaz_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, relaz_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, relaz_varid, 0, -1e30)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1',NF90_FLOAT,dims2,ch1_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2',NF90_FLOAT,dims2,ch2_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_varid, 0, -1e30)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a',NF90_FLOAT,dims2,ch3a_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_varid, 0, -1e30)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b',NF90_FLOAT,dims2,ch3b_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4',NF90_FLOAT,dims2,ch4_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_varid, 0, -1e30)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5',NF90_FLOAT,dims2,ch5_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_varid, 0, -1e30)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1_random',NF90_FLOAT,dims2,ch1_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_random_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2_random',NF90_FLOAT,dims2,ch2_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_random_varid, 0, -1e30)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a_random',NF90_FLOAT,dims2,ch3a_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_random_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_random_varid, 0, -1e30)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b_random',NF90_FLOAT,dims2,ch3b_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_random_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4_random',NF90_FLOAT,dims2,ch4_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_random_varid, 0, -1e30)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5_random',NF90_FLOAT,dims2,ch5_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_random_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_random_varid, 0, -1e30)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1_non_random',NF90_FLOAT,dims2,ch1_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_non_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_non_random_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2_non_random',NF90_FLOAT,dims2,ch2_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_non_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_non_random_varid, 0, -1e30)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a_non_random',NF90_FLOAT,dims2,ch3a_non_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_non_random_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_non_random_varid, 0, -1e30)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b_non_random',NF90_FLOAT,dims2,ch3b_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_non_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_non_random_varid, 0, -1e30)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4_non_random',NF90_FLOAT,dims2,ch4_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_non_random_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_non_random_varid, 0, -1e30)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5_non_random',NF90_FLOAT,dims2,ch5_non_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_non_random_varid, 1, 1, 9)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_non_random_varid, 0, -1e30)
       call check(stat)
    ENDIF

    dims1(1) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'quality_scanline_bitmask',NF90_UBYTE,dims1,&
         scan_qual_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, scan_qual_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,scan_qual_varid,'long_name',&
         'Bitmask for quality per scanline')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,scan_qual_varid,'flag_masks',&
         '1,2,4,8,16,32,64')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,scan_qual_varid,'flag_meanings',&
         'DO_NOT_USE, BAD_TIME, BAD_NAVIGATION, BAD_CALIBRATION, CHANNEL3A_PRESENT,SOLAR_CONTAMINATION_FAILURE,SOLAR_CONTAMINATION')
    call check(stat)

    dims2(1) = dimid_ir
    dims2(2) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'quality_channel_bitmask',NF90_UBYTE,dims2,&
         chan_qual_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, chan_qual_varid, 1, 1, 9)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,chan_qual_varid,'long_name',&
         'Bitmask for quality per channel/scanline')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,chan_qual_varid,'flag_masks',&
         '1,2')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,chan_qual_varid,'flag_meanings',&
         'BAD_CHANNEL, SOME_PIXELS_NOT_DETECTED_2SIGMA')
    call check(stat)

    !
    ! Define two global attributes
    !
    IF( AVHRR%AVHRR_No .eq. 1 )THEN
       noaa_string = 'TIROSN'
    ELSE IF( AVHRR%AVHRR_No .eq. 6 )THEN
       noaa_string = 'NOAA06'
    ELSE IF( AVHRR%AVHRR_No .eq. 7 )THEN
       noaa_string = 'NOAA07'
    ELSE IF( AVHRR%AVHRR_No .eq. 8 )THEN
       noaa_string = 'NOAA08'
    ELSE IF( AVHRR%AVHRR_No .eq. 9 )THEN
       noaa_string = 'NOAA09'
    ELSE IF( AVHRR%AVHRR_No .eq. 10 )THEN
       noaa_string = 'NOAA10'
    ELSE IF( AVHRR%AVHRR_No .eq. 11 )THEN
       noaa_string = 'NOAA11'
    ELSE IF( AVHRR%AVHRR_No .eq. 12 )THEN
       noaa_string = 'NOAA12'
    ELSE IF( AVHRR%AVHRR_No .eq. 14 )THEN
       noaa_string = 'NOAA14'
    ELSE IF( AVHRR%AVHRR_No .eq. 15 )THEN
       noaa_string = 'NOAA15'
    ELSE IF( AVHRR%AVHRR_No .eq. 16 )THEN
       noaa_string = 'NOAA16'
    ELSE IF( AVHRR%AVHRR_No .eq. 17 )THEN
       noaa_string = 'NOAA17'
    ELSE IF( AVHRR%AVHRR_No .eq. 18 )THEN
       noaa_string = 'NOAA18'
    ELSE IF( AVHRR%AVHRR_No .eq. 19 )THEN
       noaa_string = 'NOAA19'
    ELSE IF( AVHRR%AVHRR_No .eq. -1 )THEN
       noaa_string = 'METOPA'
    ELSE IF( AVHRR%AVHRR_No .eq. -2 )THEN
       noaa_string = 'METOPB'
    ELSE IF( AVHRR%AVHRR_No .eq. -3 )THEN
       noaa_string = 'METOPC'
    ENDIF
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'noaa_string',TRIM(noaa_string))
    call check(stat)    

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'version',TRIM(software_version))
    call check(stat)    

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'spatial_correlation_scale',NPIXEL_PRT_SMOOTH)
    call check(stat) 

    stat = NF90_ENDDEF(ncid)
    call check(stat)

    !
    ! Now write data
    !
!    WRITE(*,*)'Latitude writing'
    stat = NF90_PUT_VAR(ncid, latitude_varid, AVHRR%Lat)
    call check(stat)

!    WRITE(*,*)'Longitude writing'
    stat = NF90_PUT_VAR(ncid, longitude_varid, AVHRR%Lon)
    call check(stat)

!    WRITE(*,*)'Time writing'
    stat = NF90_PUT_VAR(ncid, time_varid, AVHRR%time)
    call check(stat)

!    WRITE(*,*)'year writing'
    stat = NF90_PUT_VAR(ncid, year_varid, AVHRR%year)
    call check(stat)

!    WRITE(*,*)'Month writing'
    stat = NF90_PUT_VAR(ncid, month_varid, AVHRR%month)
    call check(stat)

!    WRITE(*,*)'Day writing'
    stat = NF90_PUT_VAR(ncid, day_varid, AVHRR%day)
    call check(stat)

!    WRITE(*,*)'Hours writing'
    stat = NF90_PUT_VAR(ncid, hours_varid, AVHRR%hours)
    call check(stat)

!    WRITE(*,*)'Satza writing'
    stat = NF90_PUT_VAR(ncid, satza_varid, AVHRR%satza)
    call check(stat)

!    WRITE(*,*)'Solza writing'
    stat = NF90_PUT_VAR(ncid, solza_varid, AVHRR%solza)
    call check(stat)

    IF( ALLOCATED(AVHRR%relaz) )THEN
!       WRITE(*,*)'Relaz writing'
       stat = NF90_PUT_VAR(ncid, relaz_varid, AVHRR%relaz)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch1 writing'
    stat = NF90_PUT_VAR(ncid, ch1_varid, AVHRR%new_array1)
    call check(stat)

!    WRITE(*,*)'Ch2 writing'
    stat = NF90_PUT_VAR(ncid, ch2_varid, AVHRR%new_array2)
    call check(stat)

    IF( ALLOCATED(AVHRR%new_array3a) )THEN
!       WRITE(*,*)'Ch3a writing'
       stat = NF90_PUT_VAR(ncid, ch3a_varid, AVHRR%new_array3a)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch3b writing'
    stat = NF90_PUT_VAR(ncid, ch3b_varid, FCDR%btf3)
    call check(stat)

!    WRITE(*,*)'Ch4 writing'
    stat = NF90_PUT_VAR(ncid, ch4_varid, FCDR%btf4)
    call check(stat)

    IF( twelve_micron_there )THEN
!       WRITE(*,*)'Ch5 writing'
       stat = NF90_PUT_VAR(ncid, ch5_varid, FCDR%btf5)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch1 (Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch1_random_varid, AVHRR%new_array1_error)
    call check(stat)

!    WRITE(*,*)'Ch2 (Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch2_random_varid, AVHRR%new_array2_error)
    call check(stat)

    IF( ALLOCATED(AVHRR%new_array3a) )THEN
!       WRITE(*,*)'Ch3a (Rand) writing'
       stat = NF90_PUT_VAR(ncid, ch3a_random_varid, AVHRR%new_array3a_error)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch3b (Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch3b_random_varid, FCDR%ur3)
    call check(stat)

!    WRITE(*,*)'Ch4 (Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch4_random_varid, FCDR%ur4)
    call check(stat)

    IF( twelve_micron_there )THEN
!       WRITE(*,*)'Ch5 (Rand) writing'
       stat = NF90_PUT_VAR(ncid, ch5_random_varid, FCDR%ur5)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch1 (Non-Rand) writing'
!MT: 11-11-2017: fix problem of value not filling array     
    ALLOCATE(us1(AVHRR%nelem,AVHRR%arraySize), STAT=STAT)
!    DO I = 1,AVHRR%nelem
!       DO J = 1,AVHRR%arraySize
!          IF (ch1_varid(I,J).ne.-1e30) THEN
!             us1(I,J) = 0.03
!          ENDIF
!       ENDDO
!   ENDDO
    us1 = -1e30
    WHERE(AVHRR%new_array1.gt.-1e20)
       us1 = 0.03
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array1) )THEN
!       stat = NF90_PUT_VAR(ncid, ch1_non_random_varid, 0.03*AVHRR%new_array1_error)
       stat = NF90_PUT_VAR(ncid, ch1_non_random_varid, us1)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch2 (Non-Rand) writing'
!MT: 11-11-2017: fix problem of value not filling array     
    ALLOCATE(us2(AVHRR%nelem,AVHRR%arraySize), STAT=STAT)
    us2 = -1e30
    WHERE(AVHRR%new_array2.gt.-1e20)
       us2 = 0.05
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array2) )THEN
!       stat = NF90_PUT_VAR(ncid, ch2_non_random_varid, 0.05*AVHRR%new_array2_error)
       stat = NF90_PUT_VAR(ncid, ch2_non_random_varid, us2)
       call check(stat)
    ENDIF

!       WRITE(*,*)'Ch3a (Non-Rand) writing'
!MT: 11-11-2017: fix problem of value not filling array     
    ALLOCATE(us3a(AVHRR%nelem,AVHRR%arraySize), STAT=STAT)
    us3a = -1e30
    WHERE(AVHRR%new_array3a.gt.-1e20)
       us3a = 0.05
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array3a) )THEN
!       stat = NF90_PUT_VAR(ncid, ch3a_non_random_varid, 0.05*AVHRR%new_array3A_error)
       stat = NF90_PUT_VAR(ncid, ch3a_non_random_varid, us3a)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch3b (Non-Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch3b_non_random_varid, FCDR%us3)
    call check(stat)

!    WRITE(*,*)'Ch4 (Non-Rand) writing'
    stat = NF90_PUT_VAR(ncid, ch4_non_random_varid, FCDR%us4)
    call check(stat)

    IF( twelve_micron_there )THEN
!       WRITE(*,*)'Ch5 (Non-Rand) writing'
       stat = NF90_PUT_VAR(ncid, ch5_non_random_varid, FCDR%us5)
       call check(stat)
    ENDIF

!    stat = NF90_PUT_VAR(ncid, flag_no_detection_varid, FCDR%flag_no_detection)
!    call check(stat)

    stat = NF90_PUT_VAR(ncid, scan_qual_varid, &
         FCDR%quality_scanline_bitmask)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, chan_qual_varid, &
         FCDR%quality_channel_bitmask)
    call check(stat)

    stat = NF90_CLOSE(ncid)
    call check(stat)

  END SUBROUTINE Write_Temp_NETCDF

  !
  ! Based on Marines code calcul_urict
  !
  ! Corrected for coef update as well as science errors
  !
  SUBROUTINE calculate_urict(i,outData,coefs1,coefs2,coefs3,FCDR,&
       twelve_micron_there)

    INTEGER, INTENT(IN)                  :: i
    TYPE(AVHRR_Data), INTENT (IN)        :: outData
    REAL, DIMENSION(8,2), INTENT(IN)     :: coefs1,coefs2,coefs3
    TYPE(FIDUCEO_Data), INTENT(INOUT)    :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there

    ! Local variables
    ! prt_accuracy systematic, prt_noise random
    REAL                                 :: prt_accuracy=0.1, prt_noise=0., &
         urict, tstar, &
         aval,bval,nuc, &
         uc1=0,uc2=0,unuc=0,uaval=0,ubval=0
    INTEGER :: STAT
    REAL :: prtmean
    REAL :: utict_r
    REAL :: utict_s
    REAL :: prtsigma
    REAL :: utstar3_r
    REAL :: utstar4_r
    REAL :: utstar5_r
    REAL :: utstar3_s
    REAL :: utstar4_s
    REAL :: utstar5_s
    
    IF( .not. ALLOCATED(FCDR%urict4_r) )THEN
       ALLOCATE(FCDR%urict3_r(outdata%arraySize),&
            FCDR%urict4_r(outdata%arraySize),&
            FCDR%urict5_r(outdata%arraySize),&
            FCDR%urict3_s(outdata%arraySize),&
            FCDR%urict4_s(outdata%arraySize),&
            FCDR%urict5_s(outdata%arraySize),&
            STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate arrays',&
               'calculate_urict','fiduceo_uncertainties.f90')
       ENDIF
       FCDR%urict3_r = -1e30
       FCDR%urict4_r = -1e30
       FCDR%urict5_r = -1e30
       FCDR%urict3_s = -1e30
       FCDR%urict4_s = -1e30
       FCDR%urict5_s = -1e30
    ENDIF

    !Coefs : We use pre-launch values, for the moment the uncertainty =0
    !
    !             3  =  nuc
    !             4  =  aVal (band coefficient)
    !             5  =  bVal (band coefficient)


    !--On calcule u(Tict)
    ! Tict=1/4*Somme (Tprti)
    ! Tprti=a+b*Ti
    ! with a: systematic component 0.1 K
    !      b: random component 0.015K
    !Following the GUM: (u(tict))**2= u(a)**2+Tict**2 * u(b)**2 + u(b)**2 +sigmaTict**2 / 4

    prtmean = -1e30
    if ((outData%prt1(i) .ne. NAN_R ) &
         .and. (outData%prt2(i) .ne. NAN_R )&
         .and. (outData%prt3(i) .ne. NAN_R )&
         .and. (outData%prt4(i) .ne. NAN_R )) then
       prtmean = (outData%prt1(i)+outData%prt2(i)+&
            outData%prt3(i)+outData%prt4(i))/4.
       !
       ! Use dispersion of 4 PRTs as estimate of representivity error
       !
       ! On short timescales will be systematic
       !
       prtsigma=sqrt(((outData%prt1(i)-prtmean)**2&
            +(outData%prt2(i)-prtmean)**2&
            +(outData%prt3(i)-prtmean)**2&
            +(outData%prt4(i)-prtmean)**2)/3.)

       !
       ! Tict = (T1+T2+T3+T3)/4.
       !
       ! Random component from prt noise
       !
       ! u(Tict)**2 = (dTict/dT1)**2*u(T1)+(dTict/dT1)**2*u(T1)+(dTict/dT1)**2*u(T1)+(dTict/dT1)**2*u(T1))
       !
       ! u(Tict) **2= (1/16)*(4*u(T)**2) = (1/4)*u(T)**2
       utict_r = 0.5*prt_noise
       utict_s = sqrt(prt_accuracy**2 + prtsigma**2)

    end if


    !Ch3 
    nuc=coefs1(5,1)
    aval=coefs1(6,1)
    bval=coefs1(7,1)
    unuc=0

    !--On calcule u(Rict)
    ! T*=(Tict-A)/B
    ! Following the GUM: u(T*)**2= u(A)**2*(1/B)**2 +u(b)**2*(-(Tict-A)/B**2)**2 +(1/B)**2* u(Tict)**2
    ! u(T*)**2=(1/B)**2* u(Tict)**2
    if (prtmean .ne.  -1.00000002E+30 ) then 
       tstar=(prtmean-aval)/bval

       FCDR%Rict_c3(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)

       utstar3_s=sqrt(uaval**2*FCDR%dtstar_over_daval3(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval3(i)**2 &
            +utict_s**2*FCDR%dtstar_over_dT3(i)**2)
       utstar3_r=sqrt(utict_r**2*FCDR%dtstar_over_dT3(i)**2)
       !print*, outdata%utstar3(i)  
       FCDR%urict3_r(i)=sqrt(utstar3_r**2*FCDR%drict_over_dtstar3(i)**2)
       FCDR%urict3_s(i)=sqrt(utstar3_s**2*FCDR%drict_over_dtstar3(i)**2 + &
            unuc**2*FCDR%drict_over_dnuc3(i)**2)
       !print*, outdata%rict_c3(i), utict,outdata%urict3(i)
       !print*, "urict", outdata%urict3(i)
    end if
    !Ch4
    nuc=coefs2(5,1)
    aval=coefs2(6,1)
    bval=coefs2(7,1)
    unuc=0
    uaval=0
    ubval=0
    !--On calcule u(Rict)
    if (prtmean .ne. -1e30 )then
       tstar=(prtmean-aval)/bval

       FCDR%Rict_c4(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)

       utstar4_s=sqrt(uaval**2*FCDR%dtstar_over_daval4(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval4(i)**2 &
            +utict_s**2*FCDR%dtstar_over_dT4(i)**2)
       utstar4_r=sqrt(utict_r**2*FCDR%dtstar_over_dT4(i)**2)
       !print*, outdata%utstar3(i)  
       FCDR%urict4_r(i)=sqrt(utstar4_r**2*FCDR%drict_over_dtstar4(i)**2)
       FCDR%urict4_s(i)=sqrt(utstar4_s**2*FCDR%drict_over_dtstar4(i)**2 + &
            unuc**2*FCDR%drict_over_dnuc4(i)**2)
       !print*, outdata%rict_c4(i), utict,outdata%urict4(i)
    end if
    IF( twelve_micron_there )THEN
       !Ch5 
       nuc=coefs3(5,1)
       aval=coefs3(6,1)
       bval=coefs3(7,1)
       unuc=0
       if (prtmean .ne. -1e30 )then
          tstar=(prtmean-aval)/bval
          FCDR%Rict_c5(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)
          !print*, tstar, outdata%rict_c5(i)
          
          utstar5_s=sqrt(uaval**2*FCDR%dtstar_over_daval5(i)**2 &
               +ubval**2*FCDR%dtstar_over_dbval5(i)**2 &
               +utict_s**2*FCDR%dtstar_over_dT5(i)**2)
          utstar5_r=sqrt(utict_r**2*FCDR%dtstar_over_dT5(i)**2)
          !print*, outdata%utstar3(i)  
          FCDR%urict5_r(i)=sqrt(utstar5_r**2*FCDR%drict_over_dtstar5(i)**2)
          FCDR%urict5_s(i)=sqrt(utstar5_s**2*FCDR%drict_over_dtstar5(i)**2 + &
               unuc**2*FCDR%drict_over_dnuc5(i)**2)
       end if
    ENDIF

  END SUBROUTINE calculate_urict

  !
  ! Modified copy of Marines code
  !
  ! Corrected for coef update as well as science errors
  !
  SUBROUTINE Get_Sensitivities(i,outData,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)

    INTEGER, INTENT(IN) :: i
    TYPE(AVHRR_Data), INTENT(IN) :: outData
    REAL, INTENT(IN) :: coefs1(8,2)
    REAL, INTENT(IN) :: coefs2(8,2)
    REAL, INTENT(IN) :: coefs3(8,2)
    TYPE(FIDUCEO_Data), INTENT(INOUT) :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there

    ! Local variables
    INTEGER :: STAT
    INTEGER :: j
    REAL :: cs_cict_3
    REAL :: cs_cict_4
    REAL :: cs_cict_5
    REAL :: tstar
    REAL :: u

    IF( .not. ALLOCATED(FCDR%dtstar_over_daval3) )THEN
       ALLOCATE(FCDR%Rict_c3(outdata%arraySize),&
            FCDR%Rict_c4(outdata%arraySize),&
            FCDR%Rict_c5(outdata%arraySize),&
            FCDR%dtstar_over_dT3(outdata%arraySize),&
            FCDR%dtstar_over_dT4(outdata%arraySize),&
            FCDR%dtstar_over_dT5(outdata%arraySize),&
            FCDR%dtstar_over_daval3(outdata%arraySize),&
            FCDR%dtstar_over_daval4(outdata%arraySize),&
            FCDR%dtstar_over_daval5(outdata%arraySize),&
            FCDR%dtstar_over_dbval3(outdata%arraySize),&
            FCDR%dtstar_over_dbval4(outdata%arraySize),&
            FCDR%dtstar_over_dbval5(outdata%arraySize),&
            FCDR%dtstar_over_dnuc3(outdata%arraySize),&
            FCDR%dtstar_over_dnuc4(outdata%arraySize),&
            FCDR%dtstar_over_dnuc5(outdata%arraySize),&
            FCDR%drict_over_dtstar3(outdata%arraySize),&
            FCDR%drict_over_dtstar4(outdata%arraySize),&
            FCDR%drict_over_dtstar5(outdata%arraySize),&
            FCDR%drict_over_dnuc3(outdata%arraySize),&
            FCDR%drict_over_dnuc4(outdata%arraySize),&
            FCDR%drict_over_dnuc5(outdata%arraySize),&
            FCDR%dre_over_db03(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_db04(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_db05(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_db13(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_db14(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_db15(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da13(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da14(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da15(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da23(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da24(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da25(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dtinstr3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dtinstr4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dtinstr5(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dce3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dce4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dce5(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcict3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcict4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcict5(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_drict3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_drict4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_drict5(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs5(outdata%nelem,outdata%arraySize),&
            STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate arrays',&
               'Get_Sensitivities','fiduceo_uncertainties.f90')
       ENDIF
       !
       ! Set all values to bad initially
       !
       FCDR%Rict_c3 = NAN_R
       FCDR%Rict_c4 = NAN_R
       FCDR%Rict_c5 = NAN_R
       FCDR%dtstar_over_dT3 = NAN_R
       FCDR%dtstar_over_dT4 = NAN_R
       FCDR%dtstar_over_dT5 = NAN_R
       FCDR%dtstar_over_daval3 = NAN_R
       FCDR%dtstar_over_daval4 = NAN_R
       FCDR%dtstar_over_daval5 = NAN_R
       FCDR%dtstar_over_dbval3 = NAN_R
       FCDR%dtstar_over_dbval4 = NAN_R
       FCDR%dtstar_over_dbval5 = NAN_R
       FCDR%dtstar_over_dnuc3 = NAN_R
       FCDR%dtstar_over_dnuc4 = NAN_R
       FCDR%dtstar_over_dnuc5 = NAN_R
       FCDR%drict_over_dtstar3 = NAN_R
       FCDR%drict_over_dtstar4 = NAN_R
       FCDR%drict_over_dtstar5 = NAN_R
       FCDR%drict_over_dnuc3 = NAN_R
       FCDR%drict_over_dnuc4 = NAN_R
       FCDR%drict_over_dnuc5 = NAN_R
       FCDR%dre_over_db03 = NAN_R
       FCDR%dre_over_db04 = NAN_R
       FCDR%dre_over_db05 = NAN_R
       FCDR%dre_over_db13 = NAN_R
       FCDR%dre_over_db14 = NAN_R
       FCDR%dre_over_db15 = NAN_R
       FCDR%dre_over_da13 = NAN_R
       FCDR%dre_over_da14 = NAN_R
       FCDR%dre_over_da15 = NAN_R
       FCDR%dre_over_da23 = NAN_R
       FCDR%dre_over_da24 = NAN_R
       FCDR%dre_over_da25 = NAN_R
       FCDR%dre_over_dtinstr3 = NAN_R
       FCDR%dre_over_dtinstr4 = NAN_R
       FCDR%dre_over_dtinstr5 = NAN_R
       FCDR%dre_over_dce3 = NAN_R
       FCDR%dre_over_dce4 = NAN_R
       FCDR%dre_over_dce5 = NAN_R
       FCDR%dre_over_dcict3 = NAN_R
       FCDR%dre_over_dcict4 = NAN_R
       FCDR%dre_over_dcict5 = NAN_R
       FCDR%dre_over_drict3 = NAN_R
       FCDR%dre_over_drict4 = NAN_R
       FCDR%dre_over_drict5 = NAN_R
       FCDR%dre_over_dcs3 = NAN_R
       FCDR%dre_over_dcs4 = NAN_R
       FCDR%dre_over_dcs5 = NAN_R
    ENDIF

    !
    ! Calculate Rict radiances which was originally in Marines 
    ! average_per_orbite code
    !
    IF( outData%smoothPRT(i) .ne. NAN_R .and. outData%smoothPRT(i) .gt. 0. )THEN
       FCDR%Rict_c3(i) = convertRadiance(outData%smoothPRT(i), coefs1(5,1), &
            coefs1(6,1), coefs1(7,1))
       FCDR%Rict_c4(i) = convertRadiance(outData%smoothPRT(i), coefs2(5,1), &
            coefs2(6,1), coefs2(7,1))
       IF( twelve_micron_there )THEN
          FCDR%Rict_c5(i) = convertRadiance(outData%smoothPRT(i), coefs3(5,1), &
               coefs3(6,1), coefs3(7,1))
       ENDIF
    ENDIF

    !
    ! Taken direct from Marines code calcul_sensitivites
    !

    ! We define values that are the same for all the pixels of a scanline 
    cs_cict_3 = NAN_R
    if ( (outdata%smoothsp3(i) .ne. NAN_R) .and. (outData%smoothbb3(i) .ne. NAN_R)&
         .and.(outdata%smoothsp3(i) .gt. 0) .and. (outData%smoothbb3(i) .gt. 0))  THEN
       cs_cict_3=(outdata%smoothsp3(i)-outData%smoothbb3(i))
    end if
    cs_cict_4 = NAN_R
    if ( (outdata%smoothsp4(i) .ne. NAN_R) .and. (outData%smoothbb4(i) .ne.NAN_R)&
         .and.(outdata%smoothsp4(i) .gt. 0) .and. (outData%smoothbb4(i) .gt. 0))  THEN
       cs_cict_4=(outdata%smoothsp4(i)-outData%smoothbb4(i))
    end if
    cs_cict_5 = NAN_R
    IF( twelve_micron_there )THEN
       if ( (outdata%smoothsp5(i) .ne. NAN_R) .and. (outData%smoothbb5(i) .ne. NAN_R)&
            .and.(outdata%smoothsp5(i) .gt. 0) .and. (outData%smoothbb5(i) .gt. 0))  THEN
          cs_cict_5=(outdata%smoothsp5(i)-outData%smoothbb5(i))
       end if
    ENDIF

    !## ch3
    !###  sensitivities of Tstar ############
    FCDR%dtstar_over_dT3(i) = 1./coefs1(7,1)
    FCDR%dtstar_over_dT4(i) = 1./coefs2(7,1)
    IF( twelve_micron_there )THEN
       FCDR%dtstar_over_dT5(i) = 1./coefs3(7,1)
    ENDIF

    FCDR%dtstar_over_daval3(i)=-1./coefs1(7,1)
    FCDR%dtstar_over_daval4(i)=-1./coefs2(7,1)
    IF( twelve_micron_there )THEN
       FCDR%dtstar_over_daval5(i)=-1./coefs3(7,1)
    ENDIF

    if ((outdata%smoothprt(i) .ne. NAN_R) .and. (outdata%smoothprt(i) .gt. 0))then  
       FCDR%dtstar_over_dbval3(i)=-(outdata%smoothprt(i)-coefs1(6,1))/coefs1(7,1**2)
       FCDR%dtstar_over_dbval4(i)=-(outdata%smoothprt(i)-coefs2(6,1))/coefs2(7,1**2)
       IF( twelve_micron_there )THEN
          FCDR%dtstar_over_dbval5(i)=-(outdata%smoothprt(i)-coefs3(6,1))/coefs3(7,1**2)

       ENDIF
    end if

! This doesn't exist as Tstar = (T-a)/b with no nuc
!    FCDR%dtstar_over_dnuc3(i)=1./coefs1(5,1)
!    FCDR%dtstar_over_dnuc4(i)=1./coefs2(5,1)
!    FCDR%dtstar_over_dnuc5(i)=1./coefs3(5,1)
    FCDR%dtstar_over_dnuc3(i)=1.
    FCDR%dtstar_over_dnuc4(i)=1.
    IF( twelve_micron_there )THEN
       FCDR%dtstar_over_dnuc5(i)=1.
    ENDIF

    !###  sensitivities of Rict ############
    if ((outdata%smoothprt(i) .ne. NAN_R) .and. (outdata%smoothprt(i) .gt. 0)) then 
       tstar=(outdata%smoothprt(i)-coefs1(6,1))/coefs1(7,1)
       u=c2*coefs1(5,1)/tstar
       FCDR%drict_over_dnuc3(i)=(c1*coefs1(5,1)**2/(exp(u)-1)) &
            *(3-coefs1(5,1)*(c2/tstar)*exp(u)/(exp(u)-1))
       FCDR%drict_over_dtstar3(i)=c1*coefs1(5,1)**3*&
            (c2*coefs1(5,1)/tstar**2)&
            *exp(u)/(exp(u)-1)**2

       tstar=(outdata%smoothprt(i)-coefs2(6,1))/coefs2(7,1)
       u=c2*coefs2(5,1)/tstar
       FCDR%drict_over_dnuc4(i)=(c1*coefs2(5,1)**2/(exp(u)-1)) &
            *(3-coefs2(5,1)*(c2/tstar)*exp(u)/(exp(u)-1))
       FCDR%drict_over_dtstar4(i)=c1*coefs2(5,1)**3*&
            (c2*coefs2(5,1)/tstar**2)&
            *exp(u)/(exp(u)-1)**2

       IF( twelve_micron_there )THEN
          tstar=(outdata%smoothprt(i)-coefs3(6,1))/coefs3(7,1)
          u=c2*coefs3(5,1)/tstar
          FCDR%drict_over_dnuc5(i)=(c1*coefs3(5,1)**2/(exp(u)-1)) &
               *(3-coefs3(5,1)*(c2/tstar)*exp(u)/(exp(u)-1))
          FCDR%drict_over_dtstar5(i)=c1*coefs3(5,1)**3*&
               (c2*coefs3(5,1)/tstar**2)&
               *exp(u)/(exp(u)-1)**2
       ENDIF
    end if

    !--We calculate sensitivities that change from a pixel to another 
    !-Loop over the 409 pixels of a scan line
    do j=1,outdata%nelem
       !---Ch 3
       !----Sensitivities of earth radiance
       if ((outdata%smoothsp3(i) .ne. NAN_R) .and. &
            (outdata%smoothsp3(i) .gt. 0) .and. &
            (FCDR%rict_c3(i) .ne. NAN_R) .and. &
            (FCDR%rict_c3(i) .gt. 0) .and. &
            (outData%Counts3(j,i).ne. NAN_R) .and. &
            (outdata%Counts3(j,i) .gt. 0) .and. &
            (outdata%orbital_temperature .ne. NAN_R) .and. &
            (outdata%orbital_temperature .gt. 0) .and. &
            (cs_cict_3 .ne. NAN_R) .and. (cs_cict_3 .gt. 0)) then 

          FCDR%dre_over_db03(j,i)=1

          FCDR%dre_over_da13(j,i)=FCDR%Rict_c3(i)*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3

          FCDR%dre_over_da23(j,i)=-cs_cict_3*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               +(outdata%smoothsp3(i)-outData%Counts3(j,i))**2

          FCDR%dre_over_db13(j,i)=outdata%orbital_temperature
          FCDR%dre_over_dtinstr3(j,i)=coefs1(8,1)

          FCDR%dre_over_dce3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               /cs_cict_3 &
               -2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))

          FCDR%dre_over_dcict3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3**2&
               +2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))

          FCDR%dre_over_dcs3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3**2&
               +((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               /cs_cict_3

          FCDR%dre_over_drict3(j,i)=(eta_ict+coefs1(2,1))*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3
       end if

       !---Ch4
       !----Sensitivities  of earth radiance
       if ((outdata%smoothsp4(i) .ne. NAN_R) .and. &
            (outdata%smoothsp4(i) .gt. 0) .and. &
            (FCDR%rict_c4(i) .ne. NAN_R) .and. &
            (FCDR%rict_c4(i) .gt. 0) .and. &
            (outData%Counts4(j,i).ne. NAN_R) .and. &
            (outdata%Counts4(j,i) .gt. 0) .and. &
            (outdata%orbital_temperature .ne. NAN_R) .and. &
            (outdata%orbital_temperature .gt. 0) .and. &
            (cs_cict_4 .ne. NAN_R) .and. (cs_cict_4 .gt. 0)) then 

          FCDR%dre_over_db04(j,i)=1

          FCDR%dre_over_da14(j,i)=FCDR%Rict_c4(i)*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               /cs_cict_4

          FCDR%dre_over_da24(j,i)=-cs_cict_4*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               +(outdata%smoothsp4(i)-outData%Counts4(j,i))**2

          FCDR%dre_over_db14(j,i)=outdata%orbital_temperature
          FCDR%dre_over_dtinstr4(j,i)=coefs2(8,1)

          FCDR%dre_over_dce4(j,i)=-((eta_ict+coefs2(2,1))*FCDR%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
               /cs_cict_4 &
               -2*coefs2(4,1)*(outdata%smoothsp4(i)-outData%Counts4(j,i))

          FCDR%dre_over_dcict4(j,i)=-((eta_ict+coefs2(2,1))*FCDR%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
               *(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               /cs_cict_4**2&
               +2*coefs2(4,1)*(outdata%smoothsp4(i)-outData%Counts4(j,i))

          FCDR%dre_over_dcs4(j,i)=-((eta_ict+coefs2(2,1))*FCDR%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
               *(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               /cs_cict_4**2&
               +((eta_ict+coefs2(2,1))*FCDR%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
               /cs_cict_4

          FCDR%dre_over_drict4(j,i)=(eta_ict+coefs2(2,1))*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               /cs_cict_4
       end if

       IF( twelve_micron_there )THEN

          !---Ch5
          !----Sensitivities  of earth radiance
          if ((outdata%smoothsp5(i) .ne. NAN_R) .and. &
               (outdata%smoothsp5(i) .gt. 0) .and. &
               (FCDR%rict_c5(i) .ne. NAN_R) .and. &
               (FCDR%rict_c5(i) .gt. 0) .and. &
               (outData%Counts5(j,i).ne. NAN_R) .and. &
               (outdata%Counts5(j,i) .gt. 0) .and. &
               (outdata%orbital_temperature .ne. NAN_R) .and. &
               (outdata%orbital_temperature .gt. 0) .and. &
               (cs_cict_5 .ne. NAN_R) .and. (cs_cict_5 .gt. 0)) then 
             
             FCDR%dre_over_db05(j,i)=1

             FCDR%dre_over_da15(j,i)=FCDR%Rict_c5(i)*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5

             FCDR%dre_over_da25(j,i)=-cs_cict_5*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  +(outdata%smoothsp5(i)-outData%Counts5(j,i))**2

             FCDR%dre_over_db15(j,i)=outdata%orbital_temperature
             FCDR%dre_over_dtinstr5(j,i)=coefs3(8,1)

             FCDR%dre_over_dce5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
                  /cs_cict_5 &
                  -2*coefs2(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))

             FCDR%dre_over_dcict5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2) &
                  *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5**2&
                  +2*coefs2(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))
             
             FCDR%dre_over_dcs5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
                  *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5**2&
                  +((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
                  /cs_cict_5
             
             FCDR%dre_over_drict5(j,i)=(eta_ict+coefs3(2,1))*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5
          end if
       ENDIF
    end do 

  END SUBROUTINE Get_Sensitivities

  !
  ! Get radiance uncertainties
  !
  ! Note that only 2 sigma detections for flux are kept
  ! If fainter than 2 sigma then pixel set to NAN_R
  !
  SUBROUTINE radiance_uncertainties(i,outData,coefs1,coefs2,coefs3,FCDR,&
       twelve_micron_there)


    INTEGER, INTENT(IN) :: i
    TYPE(AVHRR_Data), INTENT(IN) :: outData
    REAL, INTENT(IN) :: coefs1(8,2)
    REAL, INTENT(IN) :: coefs2(8,2)
    REAL, INTENT(IN) :: coefs3(8,2)
    TYPE(FIDUCEO_Data), INTENT(INOUT) :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there

    ! Local variables
    INTEGER :: j
    INTEGER :: STAT
    REAL :: cs_cict_3
    REAL :: cs_cict_4
    REAL :: cs_cict_5

    REAL :: nrmax5,ur5, nrmin5, trmin5, trmax5, &
         nsmax5,us5, nsmin5, tsmin5, tsmax5, &
         nrmax3,ur3, nrmin3, trmin3, trmax3, &
         nsmax3,us3, nsmin3, tsmin3, tsmax3, &
         nrmax4,ur4, nrmin4, trmin4, trmax4, &
         nsmax4,us4, nsmin4, tsmin4, tsmax4
    REAL :: two_sigma

    !
    ! Allocate arrays
    !
    IF( .not. ALLOCATED(FCDR%ue3) )THEN
       ALLOCATE(FCDR%ue3(outData%nelem,outData%arraySize),&
            FCDR%ue4(outData%nelem,outData%arraySize),&
            FCDR%ue5(outData%nelem,outData%arraySize),&
            FCDR%ur3(outData%nelem,outData%arraySize),&
            FCDR%ur4(outData%nelem,outData%arraySize),&
            FCDR%ur5(outData%nelem,outData%arraySize),&
            FCDR%us3(outData%nelem,outData%arraySize),&
            FCDR%us4(outData%nelem,outData%arraySize),&
            FCDR%us5(outData%nelem,outData%arraySize),&
            FCDR%btf3(outData%nelem,outData%arraySize),&
            FCDR%btf4(outData%nelem,outData%arraySize),&
            FCDR%btf5(outData%nelem,outData%arraySize),&
            FCDR%uce3(outData%nelem,outData%arraySize),&
            FCDR%uce4(outData%nelem,outData%arraySize),&
            FCDR%uce5(outData%nelem,outData%arraySize),&
            FCDR%ucict3(outdata%arraySize),&
            FCDR%ucict4(outdata%arraySize),&
            FCDR%ucict5(outdata%arraySize),&
            FCDR%ucs3(outdata%arraySize),&
            FCDR%ucs4(outdata%arraySize),&
            FCDR%ucs5(outdata%arraySize),&
            FCDR%flag_no_detection(3,outData%arraySize),&
       STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate arrays',&
               'radiance_uncertainties','fiduceo_uncertainties.f90')
       ENDIF
       FCDR%ue3 = NAN_R
       FCDR%ue4 = NAN_R
       FCDR%ue5 = NAN_R
       FCDR%ur3 = NAN_R
       FCDR%ur4 = NAN_R
       FCDR%ur5 = NAN_R
       FCDR%us3 = NAN_R
       FCDR%us4 = NAN_R
       FCDR%us5 = NAN_R
       FCDR%btf3 = NAN_R
       FCDR%btf4 = NAN_R
       FCDR%btf5 = NAN_R
       FCDR%uce3 = NAN_R
       FCDR%uce4 = NAN_R
       FCDR%uce5 = NAN_R
       FCDR%ucict3 = NAN_R
       FCDR%ucict4 = NAN_R
       FCDR%ucict5 = NAN_R
       FCDR%ucs3 = NAN_R
       FCDR%ucs4 = NAN_R
       FCDR%ucs5 = NAN_R
       FCDR%flag_no_detection = 0
    ENDIF

    !
    ! Taken directly from Marines code radiance_uncertainties
    !

    !-On definit des grandeurs qui sont utiles et identiques pour tous les pixels d'une scan ligne
    cs_cict_3 = NAN_R
    if ( (outData%smoothsp3(i) .ne. NAN_R) .and. (outData%smoothbb3(i) .ne. NAN_R) &
         .and.(outData%smoothsp3(i) .gt. 0) .and. (outData%smoothbb3(i) .gt.  0))  THEN
       cs_cict_3=(outData%smoothsp3(i)-outData%smoothbb3(i))
    end if
    cs_cict_4 = NAN_R
    if ( (outDAta%smoothsp4(i) .ne. NAN_R) .and. (outData%smoothbb4(i) .ne. NAN_R)&
         .and. (outData%smoothsp4(i) .gt. 0) .and. (outData%smoothbb4(i) .gt. 0))  THEN
       cs_cict_4=(outData%smoothsp4(i)-outData%smoothbb4(i))
    end if
    cs_cict_5 = NAN_R
    IF( twelve_micron_there )THEN
       if ( (outData%smoothsp5(i) .ne. NAN_R) .and. (outData%smoothbb5(i) .ne. NAN_R)&
            .and. (outData%smoothsp5(i) .gt.0) .and. (outData%smoothbb5(i) .gt. 0))  THEN
          cs_cict_5=(outData%smoothsp5(i)-outData%smoothbb5(i))
       end if
    ENDIF

!    print *,'cs_cict_3=',cs_cict_3
!    print *,'cs_cict_4=',cs_cict_4
!    print *,'cs_cict_5=',cs_cict_5

    !--ON cacul les incertitudes qui changent d'un pixel à l'autre           
    do j=1,409 

       FCDR%ue3(j,i)=NAN_R
       FCDR%ue4(j,i)=NAN_R
       FCDR%ue5(j,i)=NAN_R
       FCDR%ur3(j,i)=NAN_R
       FCDR%ur4(j,i)=NAN_R
       FCDR%ur5(j,i)=NAN_R
       FCDR%us3(j,i)=NAN_R
       FCDR%us4(j,i)=NAN_R
       FCDR%us5(j,i)=NAN_R
       FCDR%btf3(j,i)=NAN_R
       FCDR%btf4(j,i)=NAN_R
       FCDR%btf5(j,i)=NAN_R
       FCDR%uce3(j,i)=outdata%noise_cnts(4,i)
       FCDR%uce4(j,i)=outdata%noise_cnts(5,i)
       IF( twelve_micron_there )THEN
          FCDR%uce5(j,i)=outdata%noise_cnts(6,i)
       ENDIF
!       FCDR%ucict3=outdata%noise_cnts_cal(4,i)
!       FCDR%ucict4=outdata%noise_cnts_cal(5,i)
!       IF( twelve_micron_there )THEN
!          FCDR%ucict5=outdata%noise_cnts_cal(6,i)
!       ENDIF
!       FCDR%ucs3=outdata%noise_cnts_cal(4,i)
!       FCDR%ucs4=outdata%noise_cnts_cal(5,i)
!       IF( twelve_micron_there )THEN
!          FCDR%ucs5=outdata%noise_cnts_cal(6,i)
!       ENDIF
       
!MT: 13-11-2017: allocated nsmoothBB3,4,5 and nsmoothSp3,4,5 to AVHRRout data structure in combine_orbits.f90 so that the calculations don't fail 
       IF( outdata%nsmoothBB3(i) .gt. 0 )THEN
          FCDR%ucict3(i)=outdata%noise_cnts(4,i)/SQRT(1.*outdata%nsmoothBB3(i))
       ELSE
          FCDR%ucict3(i)=NAN_R
       ENDIF
       IF( outdata%nsmoothBB4(i) .gt. 0 )THEN
          FCDR%ucict4(i)=outdata%noise_cnts(5,i)/SQRT(1.*outdata%nsmoothBB4(i))
       ELSE
          FCDR%ucict4(i)=NAN_R
       ENDIF
       IF( twelve_micron_there )THEN
          IF( outdata%nsmoothBB5(i) .gt. 0 )THEN
             FCDR%ucict5(i)=outdata%noise_cnts(6,i)/SQRT(1.*outdata%nsmoothBB5(i))
          ELSE
             FCDR%ucict5(i)=NAN_R
          ENDIF
       ENDIF
       IF( outdata%nsmoothSp3(i) .gt. 0 )THEN
          FCDR%ucs3(i)=outdata%noise_cnts(4,i)/SQRT(1.*outdata%nsmoothSp3(i))
       ELSE
          FCDR%ucs3(i)=NAN_R
       ENDIF
       IF( outdata%nsmoothSp4(i) .gt. 0 )THEN
          FCDR%ucs4(i)=outdata%noise_cnts(5,i)/SQRT(1.*outdata%nsmoothSp4(i))
       ELSE
          FCDR%ucs4(i)=NAN_R
       ENDIF
       IF( twelve_micron_there )THEN
          IF( outdata%nsmoothSp5(i) .gt. 0 )THEN
             FCDR%ucs5(i)=outdata%noise_cnts(6,i)/SQRT(1.*outdata%nsmoothSp5(i))
          ELSE
             FCDR%ucs5(i)=NAN_R
          ENDIF
       ENDIF
       ur3 = NAN_R
       us3 = NAN_R
      !---Ch 3
       if ((FCDR%uce3(j,i) .ne. NAN_R) &
            .and. (FCDR%uce3(j,i) .gt. 0) &
            .and. (FCDR%dre_over_dce3(j,i) .ne. NAN_R) &
            .and. (FCDR%urict3_r(i) .ge. 0) )then 
          ur3=sqrt(FCDR%dre_over_dce3(j,i)**2*FCDR%uce3(j,i)**2 + &
               FCDR%dre_over_drict3(j,i)**2*FCDR%urict3_r(i)**2)
       end if

       if ((FCDR%ucict3(i) .ne. NAN_R) &
            .and. (FCDR%ucict3(i) .gt. 0) &
            .and. (FCDR%ucs3(i) .ne. NAN_R) &
            .and. (FCDR%ucs3(i) .gt. 0) &
            .and. (FCDR%dre_over_drict3(j,i) .ne. NAN_R) & 
            .and. (FCDR%dre_over_dcs3(j,i) .ne. NAN_R) & 
            .and. (FCDR%dre_over_dcict3(j,i) .ne. NAN_R) &
            .and. (FCDR%urict3_s(i) .ge. 0) )then 
          IF( outData%walton_bias_correction )THEN
             ! Add in Walton bias correct uncertainty as well
             us3=sqrt((FCDR%dre_over_drict3(j,i)**2*FCDR%urict3_s(i)**2) &
                  +(FCDR%dre_over_dcs3(j,i)**2*FCDR%ucs3(i)**2) &
                  +(FCDR%dre_over_dcict3(j,i)**2*FCDR%ucict3(i)**2) &
                  +(outData%walton_bias_corr_uncert(1)**2)) 
          ELSE
             us3=sqrt((FCDR%dre_over_drict3(j,i)**2*FCDR%urict3_s(i)**2) &
                  +(FCDR%dre_over_dcs3(j,i)**2*FCDR%ucs3(i)**2) &
                  +(FCDR%dre_over_dcict3(j,i)**2*FCDR%ucict3(i)**2)) 
          ENDIF
       end if

       !---Ch 4
       ur4 = NAN_R
       us4 = NAN_R
       if ((FCDR%uce4(j,i) .ne. NAN_R) &
            .and. (FCDR%uce4(j,i) .gt. 0) &
            .and. (FCDR%dre_over_dce4(j,i) .ne. NAN_R) & 
            .and. (FCDR%urict4_r(i) .ge. 0) )then 
          ur4=sqrt(FCDR%dre_over_dce4(j,i)**2*FCDR%uce4(j,i)**2 + &
               FCDR%dre_over_drict4(j,i)**2*FCDR%urict4_r(i)**2)
       end if

       if ((FCDR%ucict4(i) .ne. NAN_R) &
            .and. (FCDR%ucict4(i) .gt. 0) &
            .and. (FCDR%ucs4(i) .ne. NAN_R) &
            .and. (FCDR%ucs4(i) .gt. 0) &
            .and. (FCDR%dre_over_drict4(j,i) .ne. NAN_R) & 
            .and. (FCDR%dre_over_dcs4(j,i) .ne. NAN_R) & 
            .and. (FCDR%dre_over_dcict4(j,i) .ne. NAN_R) &
            .and. (FCDR%urict4_s(i) .ge. 0) )then 
          IF( outData%walton_bias_correction )THEN
             us4=sqrt((FCDR%dre_over_drict4(j,i)**2*FCDR%urict4_s(i)**2) &
                  +(FCDR%dre_over_dcs4(j,i)**2*FCDR%ucs4(i)**2) &
                  +(FCDR%dre_over_dcict4(j,i)**2*FCDR%ucict4(i)**2) & 
                  +(outData%walton_bias_corr_uncert(2)**2)) 
          ELSE
             us4=sqrt((FCDR%dre_over_drict4(j,i)**2*FCDR%urict4_s(i)**2) &
                  +(FCDR%dre_over_dcs4(j,i)**2*FCDR%ucs4(i)**2) &
                  +(FCDR%dre_over_dcict4(j,i)**2*FCDR%ucict4(i)**2)) 
          ENDIF
       end if

!MT: 06-12-2017: allocate ur5 and us5 first
!       ur5 = NAN_R
!       us5 = NAN_R

       IF( twelve_micron_there )THEN
          ur5 = NAN_R
          us5 = NAN_R
          !---Ch 5
          if ((FCDR%uce5(j,i) .ne. NAN_R) &
               .and. (FCDR%uce5(j,i) .gt. 0) &
               .and. (FCDR%dre_over_dce5(j,i) .ne. NAN_R) &
               .and. (FCDR%urict5_r(i) .ge. 0) )then 
             ur5=sqrt(FCDR%dre_over_dce5(j,i)**2*FCDR%uce5(j,i)**2 + &
                  FCDR%dre_over_drict5(j,i)**2*FCDR%urict5_r(i)**2)
          end if

          if ((FCDR%ucict5(i) .ne. NAN_R) &
               .and. (FCDR%ucict5(i) .gt. 0) &
               .and. (FCDR%ucs5(i) .ne. NAN_R) &
               .and. (FCDR%ucs5(i) .gt. 0) &
               .and. (FCDR%dre_over_drict5(j,i) .ne. NAN_R) & 
               .and. (FCDR%dre_over_dcs5(j,i) .ne. NAN_R) & 
               .and. (FCDR%dre_over_dcict5(j,i) .ne. NAN_R) &
               .and. (FCDR%urict5_s(i) .ge. 0) )then 
             IF( outData%walton_bias_correction )THEN
                us5=sqrt((FCDR%dre_over_drict5(j,i)**2*FCDR%urict5_s(i)**2) &
                     +(FCDR%dre_over_dcs5(j,i)**2*FCDR%ucs5(i)**2) &
                     +(FCDR%dre_over_dcict5(j,i)**2*FCDR%ucict5(i)**2) & 
                     +(outData%walton_bias_corr_uncert(3)**2)) 
             ELSE
                us5=sqrt((FCDR%dre_over_drict5(j,i)**2*FCDR%urict5_s(i)**2) &
                     +(FCDR%dre_over_dcs5(j,i)**2*FCDR%ucs5(i)**2) &
                     +(FCDR%dre_over_dcict5(j,i)**2*FCDR%ucict5(i)**2)) 
             ENDIF
          end if
       ENDIF

       !---FIDUCEO : on convertit les radiances en BT
       if  ((outData%new_array3B(j,i) .gt. 0).and. (outData%new_array3B(j,i) .ne. NAN_R) &
            .and. ur3 .ne. NAN_R .and. us3 .ne. NAN_R ) then
          FCDR%btf3(j,i)=convertBT(outData%new_array3B(j,i),coefs1(5,1), coefs1(6,1), coefs1(7,1))
          nrmax3=outData%new_array3B(j,i)+ur3
          trmax3=convertBT(nrmax3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
          nrmin3=outData%new_array3B(j,i)-ur3
          two_sigma = outData%new_array3B(j,i)-2*ur3
          if( two_sigma .gt. 0. )then
             trmin3=convertBT(nrmin3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
             FCDR%ur3(j,i)=(trmax3-trmin3)/2.
          else
!             FCDR%ur3(j,i) = trmax3-FCDR%btf3(j,i)
             ! Set data to bad as no 1 sigma detection available
             FCDR%flag_no_detection(1,i) = 1
             FCDR%ur3(j,i) = NAN_R
             FCDR%us3(j,i) = NAN_R
             FCDR%btf3(j,i) = NAN_R
          endif

          nsmax3=outData%new_array3B(j,i)+us3
          nsmin3=outData%new_array3B(j,i)-us3
          two_sigma = outData%new_array3B(j,i)-2*us3
          tsmax3=convertBT(nsmax3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
          if( two_sigma .gt. 0. .and. NAN_R .ne. FCDR%ur3(j,i) )then
             tsmin3=convertBT(nsmin3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
             FCDR%us3(j,i)=(tsmax3-tsmin3)/2.
          else
!             FCDR%us3(j,i) = tsmax3-FCDR%btf3(j,i)
             ! Set data to bad as no 1 sigma detection available
             FCDR%flag_no_detection(1,i) = 1
             FCDR%ur3(j,i) = NAN_R
             FCDR%us3(j,i) = NAN_R
             FCDR%btf3(j,i) = NAN_R
          endif
       end if

       if  ((outData%new_array4(j,i) .gt. 0).and.(outData%new_array4(j,i) .ne. NAN_R) &
            .and. ur4 .ne. NAN_R .and. us4 .ne. NAN_R ) then
          FCDR%btf4(j,i)=convertBT(outdata%new_array4(j,i),coefs2(5,1), coefs2(6,1), coefs2(7,1))
          nrmax4=outdata%new_array4(j,i)+ur4 
          trmax4=convertBT(nrmax4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
          nrmin4=outdata%new_array4(j,i)-ur4
          two_sigma = outdata%new_array4(j,i)-2*ur4
          if( two_sigma .gt. 0 )then
             trmin4=convertBT(nrmin4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
             FCDR%ur4(j,i)=(trmax4-trmin4)/2.
          else
!             FCDR%ur4(j,i)=trmax4-FCDR%btf4(j,i)
             ! Set data to bad as no 1 sigma detection available
             FCDR%flag_no_detection(2,i) = 1
             FCDR%ur4(j,i) = NAN_R
             FCDR%us4(j,i) = NAN_R
             FCDR%btf4(j,i) = NAN_R
          endif
          nsmax4=outdata%new_array4(j,i)+us4
          tsmax4=convertBT(nsmax4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
          nsmin4=outdata%new_array4(j,i)-us4
          two_sigma = outdata%new_array4(j,i)-2*us4
          if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur4(j,i) )then
             tsmin4=convertBT(nsmin4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
             FCDR%us4(j,i)=(tsmax4-tsmin4)/2.
          else
!             FCDR%us4(j,i)=tsmax4-FCDR%btf4(j,i)
             ! Set data to bad as no 1 sigma detection available
             FCDR%flag_no_detection(2,i) = 1
             FCDR%ur4(j,i) = NAN_R
             FCDR%us4(j,i) = NAN_R
             FCDR%btf4(j,i) = NAN_R
          endif
       end if

       IF( twelve_micron_there )THEN
          if  ((outData%new_array5(j,i) .gt. 0).and. (outData%new_array5(j,i) .ne. NAN_R) &
               .and. ur5 .ne. NAN_R .and. us5 .ne. NAN_R ) then
             FCDR%btf5(j,i)=convertBT(outdata%new_array5(j,i),coefs3(5,1), coefs3(6,1), coefs3(7,1))
             nrmax5=outdata%new_array5(j,i)+ur5
             trmax5=convertBT(nrmax5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
             nrmin5=outdata%new_array5(j,i)-ur5
             trmin5=convertBT(nrmin5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
             two_sigma = outdata%new_array5(j,i)-2*ur5
             if( two_sigma .gt. 0 )then
                trmin5=convertBT(nrmin5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
                FCDR%ur5(j,i)=(trmax5-trmin5)/2.
             else
!                FCDR%ur5(j,i)=trmax5-FCDR%btf5(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(3,i) = 1
                FCDR%ur5(j,i) = NAN_R
                FCDR%us5(j,i) = NAN_R
                FCDR%btf5(j,i) = NAN_R
             endif
             nsmax5=outdata%new_array5(j,i)+us5
             tsmax5=convertBT(nsmax5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
             nsmin5=outdata%new_array5(j,i)-us5
             two_sigma = outdata%new_array5(j,i)-2*us5
             if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur5(j,i) )then
                tsmin5=convertBT(nsmin5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
                FCDR%us5(j,i)=(tsmax5-tsmin5)/2.
             else
!                FCDR%us5(j,i)=tsmax5-FCDR%btf5(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(3,i) = 1
                FCDR%ur5(j,i) = NAN_R
                FCDR%us5(j,i) = NAN_R
                FCDR%btf5(j,i) = NAN_R
             endif
          end if
       ENDIF

    end do 

  END SUBROUTINE radiance_uncertainties

  !
  ! Rescale data for writing
  !
  SUBROUTINE Rescale(AVHRR, FCDR)

    TYPE(AVHRR_Data), INTENT(INOUT) :: AVHRR
    TYPE(FIDUCEO_Data), INTENT(INOUT) :: FCDR
    
    ! Local variables
    INTEGER                        :: i,j

    REAL                           :: real_fillvalue=9.96921e+36
    INTEGER                        :: integer_fillvalue=-32767
    real, parameter                :: angle_add_offset = 0
    real, parameter                :: angle_scale_factor = 0.01

    real, parameter                :: bt_add_offset = 273.15
    real, parameter                :: bt_scale_factor = 0.01
    integer, parameter             :: bt_valid_min = -20000
    integer, parameter             :: bt_valid_max = 10000

    real, parameter                :: ref_add_offset = 0
    real, parameter                :: ref_scale_factor = 1e-4
    integer, parameter             :: ref_valid_min = 0
    integer, parameter             :: ref_valid_max = 15000

    real, parameter                :: u_add_offset = 0
    real, parameter                :: u_scale_factor = 1
    integer, parameter             :: u_valid_min = 0
    integer, parameter             :: u_valid_max = 1000

    !
    ! Taken directly from Marines code rescale
    !

    if (AVHRR%orbital_temperature .eq. NAN_R) then
       AVHRR%orbital_temperature=NAN_I
    else
       AVHRR%orbital_temperature=anint((AVHRR%orbital_temperature-bt_add_offset)/bt_scale_factor)
       if ((AVHRR%orbital_temperature .lt. -20000)  .or. (AVHRR%orbital_temperature .gt. 10000)) then
          AVHRR%orbital_temperature=NAN_I
       end if
    end if

    Do i=1, AVHRR%arraySize
       if (AVHRR%time(i) .eq. -1e30) then
          AVHRR%time(i)=real_fillvalue
       end if
       if (AVHRR%scanlinenumber(i) .eq.-1.e+30) then
          AVHRR%scanlinenumber(i)=real_fillvalue
       end if
       if (FCDR%rict_c3(i) .eq.-1.e+30) then
          FCDR%rict_c3(i)=real_fillvalue
       end if
       if (FCDR%rict_c4(i) .eq.-1.e+30) then
          FCDR%rict_c4(i)=real_fillvalue
       end if
       if (FCDR%rict_c3(i) .eq.-1.e+30) then
          FCDR%rict_c4(i)=real_fillvalue
       end if
       DO j=1,409
          if ( (AVHRR%lat(j,i) .eq. -1.e+30) .or. &
               (AVHRR%lat(j,i) .lt. -90) .or. &
               (AVHRR%lon(j,i) .gt. 90)) then
             AVHRR%lat(j,i)=integer_fillvalue
          end if
          if ((AVHRR%lon(j,i) .eq. -1.e+30) .or. &
               (AVHRR%lon(j,i) .lt. 0) .or. &
               (AVHRR%lon(j,i) .gt. 360)) then
             AVHRR%lon(j,i)=integer_fillvalue
          end if
          if (AVHRR%relaz(j,i) .eq. -1.e+30) then
             AVHRR%relaz(j,i)=real_fillvalue
          end if

          if (AVHRR%counts1(j,i) .eq. -1.e+30) then
             AVHRR%counts1(j,i)=real_fillvalue
          end if
          if (AVHRR%counts2(j,i) .eq. -1.e+30) then
             AVHRR%counts2(j,i)=real_fillvalue
          end if
          if (AVHRR%counts3(j,i) .eq. -1.e+30) then
             AVHRR%counts3(j,i)=real_fillvalue
          end if
          if (AVHRR%counts4(j,i) .eq. -1.e+30) then
             AVHRR%counts4(j,i)=real_fillvalue
          end if
          if (AVHRR%counts5(j,i) .eq. -1.e+30) then
             AVHRR%counts5(j,i)=real_fillvalue
          end if

          if (FCDR%btf3(j,i) .eq. real_fillvalue) then
             FCDR%btf3(j,i)=integer_fillvalue
          else
             !print*, AVHRR%btf3(j,i),(AVHRR%btf3(j,i)-bt_add_offset)/bt_scale_factor
             FCDR%btf3(j,i)=anint((FCDR%btf3(j,i)-bt_add_offset)/bt_scale_factor)
             !print*, AVHRR%btf3(j,i)
             if  ((FCDR%btf3(j,i) .lt. -20000)  .or. (FCDR%btf3(j,i) .gt. 10000)) then
                FCDR%btf3(j,i)=integer_fillvalue
             end if
          end if

          if (FCDR%btf4(j,i) .eq. real_fillvalue) then
             FCDR%btf4(j,i)=integer_fillvalue
          else
             FCDR%btf4(j,i)=anint((FCDR%btf4(j,i)-bt_add_offset)/bt_scale_factor)
             if  ((FCDR%btf4(j,i) .lt. -20000)  .or. (FCDR%btf4(j,i) .gt. 10000)) then
                FCDR%btf4(j,i)=integer_fillvalue
             end if
          end if
          if (FCDR%btf5(j,i) .eq. real_fillvalue) then
             FCDR%btf5(j,i)=integer_fillvalue
          else
             FCDR%btf5(j,i)=anint((FCDR%btf5(j,i)-bt_add_offset)/bt_scale_factor)
             if  ((FCDR%btf5(j,i) .lt. -20000)  .or. (FCDR%btf5(j,i) .gt. 10000)) then
                FCDR%btf5(j,i)=integer_fillvalue
             end if
          end if

          if (AVHRR%new_array3A(j,i) .eq. -1.e+30) then
             AVHRR%new_array3A(j,i)=integer_fillvalue
          else
             AVHRR%new_array3A(j,i)=anint((AVHRR%new_array3A(j,i)-ref_add_offset)/ref_scale_factor)
             if  ((AVHRR%new_array3A(j,i) .lt. -20000)  .or. (AVHRR%new_array3A(j,i) .gt. 10000)) then
                AVHRR%new_array3A(j,i)=integer_fillvalue
             end if
          end if

          if (AVHRR%new_array1(j,i) .eq. -1.e+30) then
             AVHRR%new_array1(j,i)=integer_fillvalue
          else
             AVHRR%new_array1(j,i)=anint((AVHRR%new_array1(j,i)-ref_add_offset)/ref_scale_factor)
             if  ((AVHRR%new_array1(j,i) .lt. -20000)  .or. (AVHRR%new_array1(j,i) .gt. 10000)) then
                AVHRR%new_array1(j,i)=integer_fillvalue
             end if
          end if
          if (AVHRR%new_array2(j,i) .eq. -1.e+30) then
             AVHRR%new_array2(j,i)=integer_fillvalue
          else
             AVHRR%new_array2(j,i)=anint((AVHRR%new_array2(j,i)-ref_add_offset)/ref_scale_factor)
             if  ((AVHRR%new_array2(j,i) .lt. -20000)  .or. (AVHRR%new_array2(j,i) .gt. 10000)) then
                AVHRR%new_array2(j,i)=integer_fillvalue
             end if
          end if
          if (AVHRR%satZA(j,i) .eq. -1.e+30) then
             AVHRR%satZA(j,i)=integer_fillvalue
          else
             AVHRR%satZA(j,i)=anint((AVHRR%satZA(j,i)-angle_add_offset)/angle_scale_factor)
             if  ((AVHRR%satZA(j,i) .lt. 0)  .or. (AVHRR%satZA(j,i) .gt. 9000)) then
                AVHRR%satZA(j,i)=integer_fillvalue
             end if
          end if

          if (AVHRR%solZA(j,i) .eq. -1.e+30) then
             AVHRR%solZA(j,i)=integer_fillvalue
          else
             AVHRR%solZA(j,i)=anint((AVHRR%solZA(j,i)-angle_add_offset)/angle_scale_factor)
             if  ((AVHRR%solZA(j,i) .lt. 0)  .or. (AVHRR%solZA(j,i) .gt. 18000)) then
                AVHRR%solZA(j,i)=integer_fillvalue
             end if
          end if

       END do
    end do

  END SUBROUTINE Rescale

  !
  ! Write Easy FCDR - taken straight from Marines code
  !
  SUBROUTINE fill_netcdf(filename,AVHRR,FCDR)
    character (len = *), INTENT(IN)         :: filename
    TYPE(AVHRR_Data),INTENT(IN)             :: AVHRR
    TYPE(FIDUCEO_Data),INTENT(IN)           :: FCDR

    integer                                 :: ncid

    integer, parameter                      :: NDIMS = 2 
    integer,parameter                       :: NI = 409, N10=10, N1=1,N3=3,N4=4

    character (len = *), parameter          :: I_NAME = "ni",&
         J_NAME = "nj", &
         TEN_NAME = "n10", &
         ONE_NAME = "n1", &
         THREE_NAME = "n3", &
         FOUR_NAME = "n4"

    integer                                 :: ni_dimid, nj_dimid,n10_dimid, &
         n1_dimid,n3_dimid,n4_dimid

    integer,dimension(2)                    ::   start_2 = (/ 1, 1/)
    !--ID par variable
    integer                                  :: scan_line_varid, lat_varid,lon_varid,&
         satza_varid,solza_varid,relaz_varid, &
         tict_varid, time_varid, &
         ch1_ref_varid, ch2_ref_varid,ch3a_ref_varid, &
         ch3b_bt_varid, ch4_bt_varid,ch5_bt_varid, &   
         ch1_ur_varid, ch2_ur_varid,ch3a_ur_varid, &
         ch3b_ur_varid, ch4_ur_varid,ch5_ur_varid, &
         ch1_us_varid, ch2_us_varid,ch3a_us_varid, &
         ch3b_us_varid, ch4_us_varid,ch5_us_varid
    integer  :: integer_fillvalue=-32767

    !
    ! Taken from Marines code fill_netcdf
    !

    print*, "ouverture"
    call check(nf90_open(filename, nf90_write,ncid)) 
    print*, "longitude" 
    call check(nf90_inq_varid(ncid,"longitude",lon_varid))
    call check( nf90_put_var(ncid, lon_varid,  AVHRR%lon) )
    print*, "latitude"  
    call check(nf90_inq_varid(ncid,"latitude",lat_varid))
    call check( nf90_put_var(ncid, lat_varid,  AVHRR%Lat) )
    print*, "scanline"
    call check(nf90_inq_varid(ncid,"scanline",scan_line_varid))
    call check( nf90_put_var(ncid, scan_line_varid,  AVHRR%scanlinenumber) )
    print*, "satza"
    call check(nf90_inq_varid(ncid,"satellite_zenith_angle",satza_varid))
    call check( nf90_put_var(ncid, satza_varid,  AVHRR%satZA) )
    print*, "relaz"
    call check(nf90_inq_varid(ncid,"satellite_azimuth_angle",relaz_varid))
    call check( nf90_put_var(ncid, relaz_varid,  AVHRR%relaz) )
    print*, "solaz"
    call check(nf90_inq_varid(ncid,"solar_zenith_angle",solza_varid))
    call check( nf90_put_var(ncid, solza_varid,  AVHRR%solZA) )
    print*, "time"
!MT: 19-12-2017: "Time" changed to "time" in writer
    call check(nf90_inq_varid(ncid,"time",time_varid))
    call check( nf90_put_var(ncid, time_varid,  AVHRR%time))! start = start_1, count = count_1) )
    print*, "radiances and BT"
    call check(nf90_inq_varid(ncid,"Ch1_Ref",ch1_ref_varid))
    call check(nf90_inq_varid(ncid,"Ch2_Ref",ch2_ref_varid))
    call check(nf90_inq_varid(ncid,"Ch3a_Ref",ch3a_ref_varid))
    call check(nf90_inq_varid(ncid,"Ch3b_Bt",ch3b_bt_varid))
    call check(nf90_inq_varid(ncid,"Ch4_Bt",ch4_bt_varid))
    call check(nf90_inq_varid(ncid,"Ch5_Bt",ch5_bt_varid))

    call check( nf90_put_var(ncid, ch1_ref_varid,  AVHRR%new_array1, start = start_2))! count = count_pixel) )
    call check( nf90_put_var(ncid, ch2_ref_varid,  AVHRR%new_array2, start = start_2))! count = count_pixel) ) 
    call check( nf90_put_var(ncid, ch3a_ref_varid,  AVHRR%new_array3A, start = start_2))! count = count_pixel) ) 
    call check( nf90_put_var(ncid, ch3b_bt_varid, FCDR%btf3, start =start_2))! count = count_pixel) )
    call check( nf90_put_var(ncid, ch4_bt_varid,  FCDR%btf4, start = start_2))! count = count_pixel) )
    call check( nf90_put_var(ncid, ch5_bt_varid,  FCDR%btf5, start = start_2))! count = count_pixel) )
    !print*, "tict" 
    call check(nf90_inq_varid(ncid,"T_ICT",tict_varid))
    call check(nf90_put_att(ncid, tict_varid, "UNITS","K"))
    !print * ,"24"  
    call check( nf90_put_var(ncid, tict_varid, AVHRR%orbital_temperature))! count = count_pixel) )
    !print * ,"25"  

    call check(nf90_inq_varid(ncid,"u_random_Ch1",ch1_ur_varid))
    !print * ,"26"  
    call check(nf90_inq_varid(ncid,"u_random_Ch2",ch2_ur_varid))
    !print * ,"27"  
    call check(nf90_inq_varid(ncid,"u_random_Ch3a",ch3a_ur_varid))
    !print * ,"28"  
    call check(nf90_inq_varid(ncid,"u_random_Ch3b",ch3b_ur_varid))
    !print * ,"29"  
    call check(nf90_inq_varid(ncid,"u_random_Ch4",ch4_ur_varid))
    !print * ,"30"  
    call check(nf90_inq_varid(ncid,"u_random_Ch5",ch5_ur_varid))
    !print * ,"31"  

    call check( nf90_put_var(ncid, ch1_ur_varid,  AVHRR%new_array1_error, start = start_2))! count = count_pixel) )
    !print * ,"32"  
    call check( nf90_put_var(ncid, ch2_ur_varid,  AVHRR%new_array2_error, start = start_2))! count = count_pixel) )
    !print * ,"33"  
    call check( nf90_put_var(ncid, ch3a_ur_varid,  AVHRR%new_array3A_error, start = start_2))! count = count_pixel) )
    !print * ,"34"  
    call check( nf90_put_var(ncid, ch3b_ur_varid,  FCDR%ur3, start = start_2))! count = count_pixel) )
    !print * ,"35"  
    call check( nf90_put_var(ncid, ch4_ur_varid,  FCDR%ur4, start = start_2))! count = count_pixel) )
    !print * ,"36"  
    call check( nf90_put_var(ncid, ch5_ur_varid,  FCDR%ur5, start = start_2))! count = count_pixel) )
    !print * ,"37"  

    call check(nf90_inq_varid(ncid,"u_non_random_Ch1",ch1_us_varid))
    !print * ,"38"  
    call check(nf90_inq_varid(ncid,"u_non_random_Ch2",ch2_us_varid))
    !print * ,"39"  
    call check(nf90_inq_varid(ncid,"u_non_random_Ch3a",ch3a_us_varid))
    !print * ,"40"  
    call check(nf90_inq_varid(ncid,"u_non_random_Ch3b",ch3b_us_varid))
    !print * ,"41"  
    call check(nf90_inq_varid(ncid,"u_non_random_Ch4",ch4_us_varid))
    !print * ,"42"  
    call check(nf90_inq_varid(ncid,"u_non_random_Ch5",ch5_us_varid))
    !print * ,"43"  

!MT: 07-11-2017: structured uncertainties on the reflectances are wrong
!    call check( nf90_put_var(ncid, ch1_us_varid, 0.03*AVHRR%new_array1, start = start_2))! count = count_pixel) )
    !print * ,"44"  
!    call check( nf90_put_var(ncid, ch2_us_varid,  0.05*AVHRR%new_array2, start = start_2))! count = count_pixel) )
    !print * ,"45"  
!    call check( nf90_put_var(ncid, ch3a_us_varid,  0.05*AVHRR%new_array3A, start = start_2))! count = count_pixel) )
    !print * ,"47"  
    call check( nf90_put_var(ncid, ch1_us_varid, us1, start = start_2))! count = count_pixel) )
    !print * ,"44"  
    call check( nf90_put_var(ncid, ch2_us_varid,  us2, start = start_2))! count = count_pixel) )
    !print * ,"45"  
    call check( nf90_put_var(ncid, ch3a_us_varid,  us3a, start = start_2))! count = count_pixel) )
    !print * ,"47"  
    call check( nf90_put_var(ncid, ch3b_us_varid, FCDR%us3, start = start_2))! count = count_pixel) )
    !print * ,"44"  
    call check( nf90_put_var(ncid, ch4_us_varid,  FCDR%us4, start = start_2))! count = count_pixel) )
    !print * ,"45"  
    call check( nf90_put_var(ncid, ch5_us_varid,  FCDR%us5, start = start_2))! count = count_pixel) )
    !print * ,"47"  

    !-Close the netcdf file. 
    call check( nf90_close(ncid) )
    !write(*,*) "file netcdf close"    

!    DEALLOCATE(us1)
!    DEALLOCATE(us2)
!    DEALLOCATE(us3a)

  END SUBROUTINE fill_netcdf

  !
  ! NetCDF check code - from Marines code
  !
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check

END MODULE fiduceo_uncertainties
