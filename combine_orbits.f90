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
! *  Module to combine three AVHRR orbits to create one 'good' orbit
! * equator to equator
!
! MODIFIED VERSION: M.Taylor University of Reading
! MT: 20-10-2017: 'between file3 and file4' --> 'between file4 and file5'
! MT: 24-10-2017: fix reversed logic in Resize_orbit_equator
! MT: 11-11-2017: added allocation of nmoothBB3,4,5 to fix error caused by 
!                 their absence in fiduceo_uncertainties.f90
! MT: 11-11-2017: added allocation of nmoothSp3,4,5 to fix error caused by 
!                 their absence in fiduceo_uncertainties.f90
! MT: 08-12-2017: added allocation of nsmoothPRT1,2,3,4 to fix error caused by 
!                 their absence in fiduceo_uncertainties.f90
! MT: 11-11-2017: write nmoothBB3,4,5 to AVHRRout data structure to fix error 
!                 caused by their absence in fiduceo_uncertainties.f90
! MT: 11-11-2017: write nmoothSp3,4,5 to AVHRRout data structure to fix error 
!                 caused by their absence in fiduceo_uncertainties.f90
! MT: 08-12-2017: write nmoothPrt1,2,3,4 to AVHRRout data structure to fix 
!                 error caused by their absence in fiduceo_uncertainties.f90
! MT: 11-12-2017: fix case of sub-orbit segments in Resize_orbit_equator
! MT: 09-03-2018: PyGAC geolocation
! JM: 27-04-2018: Rewrite pyGAC geolocation routines

MODULE Combine_Orbits
  
  USE GbcsKinds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsErrorHandler
  USE GbcsDateTime
  USE NOAA_LoadAVHRRLevel1B
  USE fiduceo_uncertainties
  USE fiduceo_calibration
  USE monte_carlo

! MT: routines for handling PyGAC .h5 output files containing updated geolocation
  USE NETCDF
  USE HDF5

  IMPLICIT NONE

  INTERFACE isNaN
     MODULE PROCEDURE isNaN_Real
     MODULE PROCEDURE isNaN_Dble
  END INTERFACE

  PRIVATE
  PUBLIC :: read_all_data 

CONTAINS

  LOGICAL FUNCTION Check_Scanline(AVHRR,POS)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    INTEGER, INTENT(IN) :: POS

    ! Local variables
    
    !
    ! For calibration information doesn't matter if navigation is off
    !
    IF( AVHRR%badTime(POS) .or. AVHRR%badCalibration(POS) )THEN
       Check_Scanline=.FALSE.
       RETURN
    ENDIF

    IF( (0. .gt. AVHRR%prt1(POS)) .or. (0. .gt. AVHRR%prt2(POS)) .or. &
         (0. .gt. AVHRR%prt3(POS)) .or. (0. .gt. AVHRR%prt4(POS)) )THEN
       Check_Scanline=.FALSE.
       RETURN
    ENDIF

    IF( 0. .gt. AVHRR%bb4(POS) )THEN
       Check_Scanline=.FALSE.
       RETURN
    ENDIF
    IF( ALLOCATED(AVHRR%bb5) )THEN
       IF( 0. .gt. AVHRR%bb5(POS) )THEN
          Check_Scanline=.FALSE.
          RETURN
       ENDIF
    ENDIF

    IF( 0. .gt. AVHRR%sp4(POS) )THEN
       Check_Scanline=.FALSE.
       RETURN
    ENDIF
    IF( ALLOCATED(AVHRR%sp5) )THEN
       IF( 0. .gt. AVHRR%sp5(POS) )THEN
          Check_Scanline=.FALSE.
          RETURN
       ENDIF
    ENDIF

    Check_Scanline=.TRUE.

  END FUNCTION Check_Scanline

  SUBROUTINE Copy_Scanline(AVHRR,POS1,AVHRR_New,POS2)

    TYPE(AVHRR_Data), INTENT(INOUT) :: AVHRR
    INTEGER, INTENT(IN) :: POS1
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR_New
    INTEGER, INTENT(IN) :: POS2

    AVHRR%ch3a_there(POS1) = AVHRR_New%ch3a_there(POS2)
    AVHRR%scnline_l1b(POS1) = AVHRR_new%scnline_l1b(POS2)
    AVHRR%scanLineNumber(POS1) = &
         AVHRR_New%scanLineNumber(POS2)
    AVHRR%badTop(POS1) = &
         AVHRR_New%badTop(POS2)
    AVHRR%badTime(POS1) = &
         AVHRR_New%badTime(POS2)
    AVHRR%badNavigation(POS1) = &
         AVHRR_New%badNavigation(POS2)
    AVHRR%badCalibration(POS1) = &
         AVHRR_New%badCalibration(POS2)
    AVHRR%missingLines(POS1) = &
         AVHRR_New%missingLines(POS2)
    AVHRR%Lon(:,POS1) = &
         AVHRR_New%Lon(:,POS2)
    AVHRR%Lat(:,POS1) = &
         AVHRR_New%Lat(:,POS2)
    AVHRR%satZA(:,POS1) = &
         AVHRR_New%satZA(:,POS2)
    AVHRR%solZA(:,POS1) = &
         AVHRR_New%solZA(:,POS2)
    IF( ALLOCATED(AVHRR%relAz) )THEN
       AVHRR%relAz(:,POS1) = &
            AVHRR_New%relAz(:,POS2)
    ENDIF
    AVHRR%counts1(:,POS1) = &
         AVHRR_New%counts1(:,POS2)
    AVHRR%counts2(:,POS1) = &
         AVHRR_New%counts2(:,POS2)
    AVHRR%counts3(:,POS1) = &
         AVHRR_New%counts3(:,POS2)
    AVHRR%counts4(:,POS1) = &
         AVHRR_New%counts4(:,POS2)
    IF( ALLOCATED(AVHRR%counts5) )THEN
       AVHRR%counts5(:,POS1) = &
            AVHRR_New%counts5(:,POS2)
    ENDIF
    AVHRR%array1(:,POS1) = &
         AVHRR_New%array1(:,POS2)
    AVHRR%array2(:,POS1) = &
         AVHRR_New%array2(:,POS2)
    AVHRR%array3A(:,POS1) = &
         AVHRR_New%array3A(:,POS2)
    AVHRR%array3B(:,POS1) = &
         AVHRR_New%array3B(:,POS2)
    AVHRR%array4(:,POS1) = &
         AVHRR_New%array4(:,POS2)
    AVHRR%array5(:,POS1) = &
         AVHRR_New%array5(:,POS2)
    AVHRR%year(POS1) = &
         AVHRR_New%year(POS2)
    AVHRR%month(POS1) = &
         AVHRR_New%month(POS2)
    AVHRR%day(POS1) = &
         AVHRR_New%day(POS2)
    AVHRR%dayNo(POS1) = &
         AVHRR_New%dayNo(POS2)
    AVHRR%hours(POS1) = &
         AVHRR_New%hours(POS2)
    AVHRR%UTC_msecs(POS1) = &
         AVHRR_New%UTC_msecs(POS2)
    AVHRR%time(POS1) = &
         AVHRR_New%time(POS2)
    AVHRR%prt1(POS1) = &
         AVHRR_New%prt1(POS2)
    AVHRR%prt2(POS1) = &
         AVHRR_New%prt2(POS2)
    AVHRR%prt3(POS1) = &
         AVHRR_New%prt3(POS2)
    AVHRR%prt4(POS1) = &
         AVHRR_New%prt4(POS2)
    AVHRR%prt1Counts(POS1) = &
         AVHRR_New%prt1Counts(POS2)
    AVHRR%prt2Counts(POS1) = &
         AVHRR_New%prt2Counts(POS2)
    AVHRR%prt3Counts(POS1) = &
         AVHRR_New%prt3Counts(POS2)
    AVHRR%prt4Counts(POS1) = &
         AVHRR_New%prt4Counts(POS2)
    AVHRR%prt1CountsAll(:,POS1) = &
         AVHRR_New%prt1CountsAll(:,POS2)
    AVHRR%prt2CountsAll(:,POS1) = &
         AVHRR_New%prt2CountsAll(:,POS2)
    AVHRR%prt3CountsAll(:,POS1) = &
         AVHRR_New%prt3CountsAll(:,POS2)
    AVHRR%prt4CountsAll(:,POS1) = &
         AVHRR_New%prt4CountsAll(:,POS2)
    AVHRR%bb3(POS1) = &
         AVHRR_New%bb3(POS2)
    AVHRR%bb4(POS1) = &
         AVHRR_New%bb4(POS2)
    IF( ALLOCATED(AVHRR%bb5) )THEN
       AVHRR%bb5(POS1) = &
            AVHRR_New%bb5(POS2)
    ENDIF
    AVHRR%sp3(POS1) = &
         AVHRR_New%sp3(POS2)
    AVHRR%sp4(POS1) = &
         AVHRR_New%sp4(POS2)
    IF( ALLOCATED(AVHRR%sp5) )THEN
       AVHRR%sp5(POS1) = &
            AVHRR_New%sp5(POS2)
    ENDIF
    AVHRR%bbodyFilter3(:,POS1) = &
         AVHRR_New%bbodyFilter3(:,POS2)
    AVHRR%bbodyFilter4(:,POS1) = &
         AVHRR_New%bbodyFilter4(:,POS2)
    IF( ALLOCATED(AVHRR%bbodyFilter5) )THEN
       AVHRR%bbodyFilter5(:,POS1) = &
            AVHRR_New%bbodyFilter5(:,POS2)
    ENDIF
    AVHRR%spaceFilter1(:,POS1) = &
         AVHRR_New%spaceFilter1(:,POS2)
    AVHRR%spaceFilter2(:,POS1) = &
         AVHRR_New%spaceFilter2(:,POS2)
    AVHRR%spaceFilter3a(:,POS1) = &
         AVHRR_New%spaceFilter3a(:,POS2)
    AVHRR%spaceFilter3(:,POS1) = &
         AVHRR_New%spaceFilter3(:,POS2)
    AVHRR%spaceFilter4(:,POS1) = &
         AVHRR_New%spaceFilter4(:,POS2)
    IF( ALLOCATED(AVHRR%spaceFilter5) )THEN
       AVHRR%spaceFilter5(:,POS1) = &
            AVHRR_New%spaceFilter5(:,POS2)
    ENDIF
    AVHRR%patch(POS1) = &
         AVHRR_New%patch(POS2)
    IF( ALLOCATED(AVHRR%patchExtended) )THEN
       AVHRR%patchExtended(POS1) = &
            AVHRR_New%patchExtended(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%Radiator) )THEN
       AVHRR%Radiator(POS1) = &
            AVHRR_New%Radiator(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%Cooler) )THEN
       AVHRR%Cooler(POS1) = &
            AVHRR_New%Cooler(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%a_d_conv) )THEN
       AVHRR%a_d_conv(POS1) = &
            AVHRR_New%a_d_conv(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%motor) )THEN
       AVHRR%motor(POS1) = &
            AVHRR_New%motor(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%motorCurrent) )THEN
       AVHRR%motorCurrent(POS1) = &
            AVHRR_New%motorCurrent(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%electronics) )THEN
       AVHRR%electronics(POS1) = &
            AVHRR_New%electronics(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%baseplate) )THEN
       AVHRR%baseplate(POS1) = &
            AVHRR_New%baseplate(POS2)
    ENDIF
    AVHRR%calib1(:,POS1) = &
         AVHRR_New%calib1(:,POS2)
    AVHRR%calib1_2(:,POS1) = &
         AVHRR_New%calib1_2(:,POS2)
    AVHRR%calib1_intercept(POS1) = &
         AVHRR_New%calib1_intercept(POS2)
    AVHRR%calib2(:,POS1) = &
         AVHRR_New%calib2(:,POS2)
    AVHRR%calib2_2(:,POS1) = &
         AVHRR_New%calib2_2(:,POS2)
    AVHRR%calib2_intercept(POS1) = &
         AVHRR_New%calib2_intercept(POS2)
    IF( ALLOCATED(AVHRR%calib3A) )THEN
       AVHRR%calib3A(:,POS1) = &
            AVHRR_New%calib3A(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%calib3A_2) )THEN
       AVHRR%calib3A_2(:,POS1) = &
            AVHRR_New%calib3A_2(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%calib3A_intercept) )THEN
       AVHRR%calib3A_intercept(POS1) = &
            AVHRR_New%calib3A_intercept(POS2)
    ENDIF
    AVHRR%calib3(:,POS1) = &
         AVHRR_New%calib3(:,POS2)
    AVHRR%calib4(:,POS1) = &
         AVHRR_New%calib4(:,POS2)
    IF( ALLOCATED(AVHRR%calib5) )THEN
       AVHRR%calib5(:,POS1) = &
            AVHRR_New%calib5(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%clavr_mask) )THEN
       AVHRR%clavr_mask(:,POS1) = &
            AVHRR_New%clavr_mask(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%clavrx_mask) )THEN
       AVHRR%clavrx_mask(:,POS1) = &
            AVHRR_New%clavrx_mask(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%clavrx_prb) )THEN
       AVHRR%clavrx_prb(:,POS1) = &
            AVHRR_New%clavrx_prb(:,POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%orig_solar_contamination_3B) )THEN
       AVHRR%orig_solar_contamination_3B(POS1) = &
            AVHRR_New%orig_solar_contamination_3B(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%orig_solar_contamination_4) )THEN
       AVHRR%orig_solar_contamination_4(POS1) = &
            AVHRR_New%orig_solar_contamination_4(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%orig_solar_contamination_5) )THEN
       AVHRR%orig_solar_contamination_5(POS1) = &
            AVHRR_New%orig_solar_contamination_5(POS2)
    ENDIF
    IF( ALLOCATED(AVHRR%satelliteAltitude) )THEN
       AVHRR%satelliteAltitude(POS1) = &
            AVHRR_New%satelliteAltitude(POS2)
    ENDIF

  END SUBROUTINE Copy_Scanline
  
  LOGICAL FUNCTION corr_diff_time(time,I,size,newtime,J,outJ)

    REAL(GbcsDble), INTENT(IN) :: time(:)
    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(IN) :: size
    REAL(GbcsDble), INTENT(IN) :: newtime(:)
    INTEGER, INTENT(IN) :: J
    INTEGER, INTENT(OUT) :: outJ

    ! Local variables
    INTEGER :: K
    REAL(GbcsDble) :: min_diff
    REAL(GbcsDble) :: diff

    corr_diff_time = .TRUE.
    min_diff = 1d30
    outJ = -1
    FindLoop: DO K=-3,3 
       IF( J+K .gt. 0 .and. J+K .le. size )THEN
          IF( time(I) .gt. 0 .and. newtime(J+K) .gt. 0 )THEN
             diff = ABS(time(I)-newtime(J+K))
             IF( diff .gt. 2 )THEN
                EXIT FindLoop
             ENDIF
             IF( diff .lt. min_diff )THEN
                outJ = J+K
                min_diff = diff
             ENDIF
          ENDIF
       ENDIF
    END DO FindLoop
    IF( -1 .eq. outJ )THEN
       corr_diff_time = .FALSE.
    ENDIF
    IF( min_diff .gt. 0.5 )THEN
       corr_diff_time = .FALSE.
    ENDIF

  END FUNCTION corr_diff_time

  SUBROUTINE Merge_AVHRR(AVHRR,AVHRR_New,test)

    TYPE(AVHRR_Data), INTENT(INOUT), TARGET :: AVHRR
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR_new
    INTEGER, OPTIONAL :: test

    ! Local variables
    INTEGER :: I
    INTEGER :: J
    INTEGER :: outJ
    INTEGER :: stat
    INTEGER :: start_pos
    INTEGER :: extra_lines
    INTEGER :: first_position
    INTEGER :: first_position_new
    INTEGER :: last_position
    INTEGER :: last_position_new
    INTEGER :: temp_position
    INTEGER :: back_pos
    INTEGER :: back_pos2
    INTEGER :: offset
    INTEGER :: test_print

    REAL(GbcsDble) :: start_time 

    LOGICAL, ALLOCATABLE :: good_data(:)

    TYPE(AVHRR_Data), POINTER :: AVHRR_Ptr

    IF( PRESENT(test) )THEN
       test_print = test
    ELSE
       test_print = 0
    ENDIF

    AVHRR%valid_data_there = .TRUE.

    !
    ! Setup AVHRR pointer
    !
    AVHRR_Ptr => AVHRR

    !
    ! Get first valid scanline
    !
    start_time = -1
    start_pos = -1
    I=0
    DO WHILE(-1 .eq. start_time)
       I = I+1
       IF( I .gt. AVHRR_New%arraySize )THEN
          CALL Gbcs_Warning(.TRUE.,'No good data found','Merge_AVHRR',&
               'extract_l1b_data.f90')
          RETURN
       ENDIF
       IF( .not. isNaN(AVHRR_New%time(I)) )THEN
          IF( AVHRR_New%time(I) .gt. 0 )THEN
             start_time = AVHRR_New%time(I)
             start_pos = I
          ENDIF
       ENDIF
    END DO

    !
    ! Start comparing data and fill out missing data at end of original set
    !
    ALLOCATE(good_data(AVHRR%arraySize),STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate good_data',&
            'Merge_AVHRR','extract_l1b_data.f90')
    ENDIF
    good_data=.TRUE.
    DO I=1,AVHRR%arraySize
       !
       ! Check that line has a time
       !
       IF( .not. isNaN(AVHRR%time(I)) .and. 0 .lt. AVHRR%time(I) )THEN
          ! 
          ! Check if in overlap
          !
          IF( ABS(AVHRR%time(I)) .ge. start_time )THEN
             !
             ! Check for bad navigation/flags etc.
             !
             IF( .not. Check_Scanline(AVHRR,I) )THEN
                !
                ! Bad data in first section
                !
                FindLoop2: DO J=start_pos,AVHRR_New%arraySize
                   !
                   ! As times don't exactly match have to look for
                   !
                   IF( corr_diff_time(AVHRR%time,I,AVHRR_New%arraySize,&
                        AVHRR_New%time,J,outJ) )THEN
                      ! 
                      ! Found same line
                      !
                      IF( Check_Scanline(AVHRR_New,outJ) )THEN
                         IF( test_print .gt. 0 )THEN
                            WRITE(test_print,*)I,outJ,&
                                 AVHRR%Time(I),AVHRR_New%Time(outJ)
                         ENDIF
                         CALL Copy_Scanline(AVHRR,I,AVHRR_New,outJ)
                      ENDIF
                      EXIT FindLoop2
                   ENDIF
                END DO FindLoop2
             ENDIF
          ENDIF
       ELSE
          !
          ! See if we are at a point which overlaps with new data
          !
          ! If not time present fill on the basis of the last good time
          !
          !
          ! Look back to find last good time to work out how to fill
          !
          back_pos=-1
          FindBack: DO J=I,1,-1
             IF( 0. .lt. AVHRR%time(J) )THEN
                back_pos=J
                EXIT findBack
             ENDIF
          END DO FindBack
          !
          ! Found previous time
          !
          IF( 0 .lt. back_pos )THEN
             !
             ! Check it is in overlap
             ! 
             IF( AVHRR%time(back_pos) .ge. start_time )THEN
                ! 
                ! Calculate offset from last time
                !
                offset = I-back_pos
                !
                ! Find overlap position for last time
                !
                back_pos2 = -1
                FindLoop: DO J=start_pos,AVHRR_New%arraySize
                   !
                   ! As times don't exactly match have to look for
                   !
                   IF( corr_diff_time(AVHRR%time,back_pos,&
                        AVHRR_New%arraySize,&
                        AVHRR_New%time,J,outJ) )THEN
                      back_pos2=outJ
                      EXIT FindLoop
                   ENDIF
                END DO FindLoop
                IF( 0 .lt. back_pos2 )THEN
                   !
                   ! Now check new scanline with offset
                   !
                   IF( Check_Scanline(AVHRR_New,back_pos2+offset) )THEN
                      CALL Copy_Scanline(AVHRR,back_pos+offset,&
                           AVHRR_New,back_pos2+offset)
                   ELSE
                      good_data(I) = .FALSE.
                   ENDIF
                ELSE
                   good_data(I) = .FALSE.
                ENDIF
             ELSE
                good_data(I) = .FALSE.
             ENDIF
          ELSE
             good_data(I) = .FALSE.
          ENDIF
       ENDIF
    END DO

    !
    ! Now AVHRR structure should be as full as it can be find last good line
    !
    last_position=-1
    FindLast: DO I=AVHRR%arraySize,1,-1
       IF( good_data(I) .and. AVHRR%Time(I) .gt. 0 )THEN
          last_position=I
          EXIT FindLast
       ENDIF
    END DO FindLast
    IF( -1 .eq. last_position )THEN
       CALL Gbcs_Critical(.TRUE.,'No good data found in first input structure',&
            'Merge_AVHRR','extract_l1b_data.f90')
    ENDIF

    ! 
    ! Find equivalent point in new input
    !
    first_position=-1
    FindFirst: DO I=1,AVHRR_New%arraySize
       IF( corr_diff_time(AVHRR%time,last_position,AVHRR_New%arraySize,AVHRR_New%time,I,outJ) )THEN
          first_position = outJ+1
          IF( first_position .gt. AVHRR_New%arraySize )THEN
             first_position = -2
          ENDIF
          EXIT FindFirst
       ENDIF
    END DO FindFirst
    !
    ! Find first/last valid point in new data stream so we can check if 
    ! there is a gap from previous or we are in overlap mode
    !
    first_position_new=-1
    FindFirst2: DO I=1,AVHRR_New%arraySize
       IF( .not. isNaN(AVHRR_New%time(I)) .and. AVHRR_New%time(I) .gt. 0 )THEN
             first_position_new = I
          EXIT FindFirst2
       ENDIF
    END DO FindFirst2
    last_position_new=-1
    FindFirst3: DO I=AVHRR_New%arraySize,1,-1
       IF( .not. isNaN(AVHRR_New%time(I)) .and. AVHRR_New%time(I) .gt. 0 )THEN
             last_position_new = I
          EXIT FindFirst3
       ENDIF
    END DO FindFirst3
    !
    ! Need to check if full overlap (since we will have checked/copied all
    ! new data) or we have a gap (need to add in new data)
    !
    IF( -1 .eq. first_position .and. 0 .lt. first_position_new )THEN
       !
       ! Valid last position
       !
       IF( last_position .le. AVHRR%arraySize )THEN
          !
          ! Is data contained?
          !
          IF( corr_diff_time(AVHRR%time,last_position,AVHRR_New%arraySize,&
               AVHRR_New%time,temp_position,outJ) )THEN
             ! No gap
             continue
          ELSE
             !
             ! Look to see if we have a gap
             !
             IF( AVHRR%time(last_position) .lt. &
                  AVHRR_New%time(first_position_new) )THEN
                first_position = first_position_new
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    !
    ! Check to see if we are full overlap (last new within AVHRR structure)
    ! If so don't append
    !
    IF( last_position .gt. 0 .and. last_position_new .gt. 0 )THEN
       IF( corr_diff_time(AVHRR%time,last_position,last_position_new,&
            AVHRR_New%time,temp_position,outJ) )THEN
          first_position = -1
       ELSE IF( AVHRR%time(last_position) .ge. &
            AVHRR_New%time(last_position_new) )THEN
          first_position = -1
       ENDIF
    ENDIF

    !
    ! Append new data to old structure
    !
    IF( first_position .gt. 0 )THEN
!       print *,first_position,last_position
!       print *,first_position_new,last_position_new
!       print *,AVHRR%time(first_position),AVHRR%time(last_position)
!       print *,AVHRR_New%time(first_position_new),AVHRR_New%time(last_position_new)
       extra_lines = AVHRR_New%arraySize - first_position - &
            (AVHRR%arraySize-last_position)
       CALL Reallocate_outData(AVHRR_Ptr,extra_lines)

       !
       ! Now add in extra data - original raw data
       !
       AVHRR%ch3a_there(last_position:AVHRR%arraySize) = &
            AVHRR_New%ch3a_there(first_position:AVHRR_New%arraySize)
       AVHRR%scnline_l1b(last_position:AVHRR%arraySize) = &
            AVHRR_New%scnline_l1b(first_position:AVHRR_New%arraySize)
       AVHRR%scanLineNumber(last_position:AVHRR%arraySize) = &
            AVHRR_New%scanLineNumber(first_position:AVHRR_New%arraySize)
       AVHRR%badTop(last_position:AVHRR%arraySize) = &
            AVHRR_New%badTop(first_position:AVHRR_New%arraySize)
       AVHRR%badTime(last_position:AVHRR%arraySize) = &
            AVHRR_New%badTime(first_position:AVHRR_New%arraySize)
       AVHRR%badNavigation(last_position:AVHRR%arraySize) = &
            AVHRR_New%badNavigation(first_position:AVHRR_New%arraySize)
       AVHRR%badCalibration(last_position:AVHRR%arraySize) = &
            AVHRR_New%badCalibration(first_position:AVHRR_New%arraySize)
       AVHRR%missingLines(last_position:AVHRR%arraySize) = &
            AVHRR_New%missingLines(first_position:AVHRR_New%arraySize)
       AVHRR%Lon(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%Lon(:,first_position:AVHRR_New%arraySize)
       AVHRR%Lat(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%Lat(:,first_position:AVHRR_New%arraySize)
       AVHRR%satZA(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%satZA(:,first_position:AVHRR_New%arraySize)
       AVHRR%solZA(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%solZA(:,first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%relAz) )THEN
          AVHRR%relAz(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%relAz(:,first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%counts1(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%counts1(:,first_position:AVHRR_New%arraySize)
       AVHRR%counts2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%counts2(:,first_position:AVHRR_New%arraySize)
       AVHRR%counts3(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%counts3(:,first_position:AVHRR_New%arraySize)
       AVHRR%counts4(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%counts4(:,first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%counts5) )THEN
          AVHRR%counts5(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%counts5(:,first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%array1(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array1(:,first_position:AVHRR_New%arraySize)
       AVHRR%array2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array2(:,first_position:AVHRR_New%arraySize)
       AVHRR%array3A(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array3A(:,first_position:AVHRR_New%arraySize)
       AVHRR%array3B(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array3B(:,first_position:AVHRR_New%arraySize)
       AVHRR%array4(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array4(:,first_position:AVHRR_New%arraySize)
       AVHRR%array5(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%array5(:,first_position:AVHRR_New%arraySize)
       AVHRR%year(last_position:AVHRR%arraySize) = &
            AVHRR_New%year(first_position:AVHRR_New%arraySize)
       AVHRR%month(last_position:AVHRR%arraySize) = &
            AVHRR_New%month(first_position:AVHRR_New%arraySize)
       AVHRR%day(last_position:AVHRR%arraySize) = &
            AVHRR_New%day(first_position:AVHRR_New%arraySize)
       AVHRR%dayNo(last_position:AVHRR%arraySize) = &
            AVHRR_New%dayNo(first_position:AVHRR_New%arraySize)
       AVHRR%hours(last_position:AVHRR%arraySize) = &
            AVHRR_New%hours(first_position:AVHRR_New%arraySize)
       AVHRR%UTC_msecs(last_position:AVHRR%arraySize) = &
            AVHRR_New%UTC_msecs(first_position:AVHRR_New%arraySize)
       AVHRR%time(last_position:AVHRR%arraySize) = &
            AVHRR_New%time(first_position:AVHRR_New%arraySize)
       AVHRR%prt1(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt1(first_position:AVHRR_New%arraySize)
       AVHRR%prt2(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt2(first_position:AVHRR_New%arraySize)
       AVHRR%prt3(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt3(first_position:AVHRR_New%arraySize)
       AVHRR%prt4(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt4(first_position:AVHRR_New%arraySize)
       AVHRR%prt1Counts(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt1Counts(first_position:AVHRR_New%arraySize)
       AVHRR%prt2Counts(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt2Counts(first_position:AVHRR_New%arraySize)
       AVHRR%prt3Counts(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt3Counts(first_position:AVHRR_New%arraySize)
       AVHRR%prt4Counts(last_position:AVHRR%arraySize) = &
            AVHRR_New%prt4Counts(first_position:AVHRR_New%arraySize)
       AVHRR%prt1CountsAll(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%prt1CountsAll(:,first_position:AVHRR_New%arraySize)
       AVHRR%prt2CountsAll(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%prt2CountsAll(:,first_position:AVHRR_New%arraySize)
       AVHRR%prt3CountsAll(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%prt3CountsAll(:,first_position:AVHRR_New%arraySize)
       AVHRR%prt4CountsAll(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%prt4CountsAll(:,first_position:AVHRR_New%arraySize)
       AVHRR%bb3(last_position:AVHRR%arraySize) = &
            AVHRR_New%bb3(first_position:AVHRR_New%arraySize)
       AVHRR%bb4(last_position:AVHRR%arraySize) = &
            AVHRR_New%bb4(first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%bb5) )THEN
          AVHRR%bb5(last_position:AVHRR%arraySize) = &
               AVHRR_New%bb5(first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%sp3(last_position:AVHRR%arraySize) = &
            AVHRR_New%sp3(first_position:AVHRR_New%arraySize)
       AVHRR%sp4(last_position:AVHRR%arraySize) = &
            AVHRR_New%sp4(first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%sp5) )THEN
          AVHRR%sp5(last_position:AVHRR%arraySize) = &
               AVHRR_New%sp5(first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%bbodyFilter3(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%bbodyFilter3(:,first_position:AVHRR_New%arraySize)
       AVHRR%bbodyFilter4(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%bbodyFilter4(:,first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%bbodyFilter5) )THEN
          AVHRR%bbodyFilter5(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%bbodyFilter5(:,first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%spaceFilter1(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%spaceFilter1(:,first_position:AVHRR_New%arraySize)
       AVHRR%spaceFilter2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%spaceFilter2(:,first_position:AVHRR_New%arraySize)
       AVHRR%spaceFilter3a(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%spaceFilter3a(:,first_position:AVHRR_New%arraySize)
       AVHRR%spaceFilter3(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%spaceFilter3(:,first_position:AVHRR_New%arraySize)
       AVHRR%spaceFilter4(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%spaceFilter4(:,first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%spaceFilter5) )THEN
          AVHRR%spaceFilter5(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%spaceFilter5(:,first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%patch(last_position:AVHRR%arraySize) = &
            AVHRR_New%patch(first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%patchExtended) )THEN
          AVHRR%patchExtended(last_position:AVHRR%arraySize) = &
               AVHRR_New%patchExtended(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%Radiator) )THEN
          AVHRR%Radiator(last_position:AVHRR%arraySize) = &
               AVHRR_New%Radiator(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%Cooler) )THEN
          AVHRR%Cooler(last_position:AVHRR%arraySize) = &
               AVHRR_New%Cooler(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%a_d_conv) )THEN
          AVHRR%a_d_conv(last_position:AVHRR%arraySize) = &
               AVHRR_New%a_d_conv(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%motor) )THEN
          AVHRR%motor(last_position:AVHRR%arraySize) = &
               AVHRR_New%motor(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%motorCurrent) )THEN
          AVHRR%motorCurrent(last_position:AVHRR%arraySize) = &
               AVHRR_New%motorCurrent(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%electronics) )THEN
          AVHRR%electronics(last_position:AVHRR%arraySize) = &
               AVHRR_New%electronics(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%baseplate) )THEN
          AVHRR%baseplate(last_position:AVHRR%arraySize) = &
               AVHRR_New%baseplate(first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%calib1(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib1(:,first_position:AVHRR_New%arraySize)
       AVHRR%calib1_2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib1_2(:,first_position:AVHRR_New%arraySize)
       AVHRR%calib1_intercept(last_position:AVHRR%arraySize) = &
            AVHRR_New%calib1_intercept(first_position:AVHRR_New%arraySize)
       AVHRR%calib2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib2(:,first_position:AVHRR_New%arraySize)
       AVHRR%calib2_2(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib2_2(:,first_position:AVHRR_New%arraySize)
       AVHRR%calib2_intercept(last_position:AVHRR%arraySize) = &
            AVHRR_New%calib2_intercept(first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%calib3A) )THEN
          AVHRR%calib3A(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%calib3A(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%calib3A_2) )THEN
          AVHRR%calib3A_2(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%calib3A_2(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%calib3A_intercept) )THEN
          AVHRR%calib3A_intercept(last_position:AVHRR%arraySize) = &
               AVHRR_New%calib3A_intercept(first_position:AVHRR_New%arraySize)
       ENDIF
       AVHRR%calib3(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib3(:,first_position:AVHRR_New%arraySize)
       AVHRR%calib4(:,last_position:AVHRR%arraySize) = &
            AVHRR_New%calib4(:,first_position:AVHRR_New%arraySize)
       IF( ALLOCATED(AVHRR%calib5) )THEN
          AVHRR%calib5(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%calib5(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%clavr_mask) )THEN
          AVHRR%clavr_mask(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%clavr_mask(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%clavrx_mask) )THEN
          AVHRR%clavrx_mask(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%clavrx_mask(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%clavrx_prb) )THEN
          AVHRR%clavrx_prb(:,last_position:AVHRR%arraySize) = &
               AVHRR_New%clavrx_prb(:,first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%orig_solar_contamination_3B) )THEN
          AVHRR%orig_solar_contamination_3B(last_position:AVHRR%arraySize) = &
               AVHRR_New%orig_solar_contamination_3B(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%orig_solar_contamination_4) )THEN
          AVHRR%orig_solar_contamination_4(last_position:AVHRR%arraySize) = &
               AVHRR_New%orig_solar_contamination_4(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%orig_solar_contamination_5) )THEN
          AVHRR%orig_solar_contamination_5(last_position:AVHRR%arraySize) = &
               AVHRR_New%orig_solar_contamination_5(first_position:AVHRR_New%arraySize)
       ENDIF
       IF( ALLOCATED(AVHRR%satelliteAltitude) )THEN
          AVHRR%satelliteAltitude(last_position:AVHRR%arraySize) = &
               AVHRR_New%satelliteAltitude(first_position:AVHRR_New%arraySize)
       ENDIF
    ENDIF
    AVHRR%newCalibration_there = .FALSE.
    AVHRR%walton_there = .FALSE.
        
  END SUBROUTINE Merge_AVHRR

  LOGICAL FUNCTION check_overlap(file1,file2)

    CHARACTER(LEN=*), INTENT(IN) :: file1
    CHARACTER(LEN=*), INTENT(IN) :: file2

    ! Local variables
    INTEGER POS_Dir1
    INTEGER POS_Dir2
    INTEGER POS_Year1
    INTEGER POS_Year2
    INTEGER POS_Start1
    INTEGER POS_Start2
    INTEGER POS_End1
    INTEGER POS_End2

    REAL(GbcsDble) :: jday1_start
    REAL(GbcsDble) :: jday1_end
    REAL(GbcsDble) :: jday2_start
    REAL(GbcsDble) :: jday2_end

    !
    ! If directory included find position of last /
    ! Note if '/' not there then INDEX returns zero which is OK
    !
    POS_Dir1 = INDEX(file1,'/',.TRUE.)
    POS_Dir2 = INDEX(file2,'/',.TRUE.)

    !
    ! Check filelength - must have at least 30 characters for NOAA Level 1B data
    !
    IF( LEN_TRIM(file1)-POS_Dir1 .lt. 30 )THEN
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (length)',&
            'check_overlap','extract_l1b_data.f90')
    ENDIF
    IF( LEN_TRIM(file2)-POS_Dir2 .lt. 30 )THEN
       CALL Gbcs_Critical(.TRUE.,'File2 not a NOAA Level 1B file (length)',&
            'check_overlap','extract_l1b_data.f90')
    ENDIF

    !
    ! Check this is an AVHRR file (NOAA Level 1B)
    ! If it is, get locations of year, day, hour, min
    !
    IF( 'NSS.GHRR.' .ne. file1(POS_Dir1+1:POS_Dir1+9) )THEN
       WRITE(*,*)TRIM(file1)
       WRITE(*,*)file1(POS_Dir1+1:POS_Dir1+9)
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (NSS.)',&
            'check_overlap','extract_l1b_data.f90')
    ENDIF
    IF( 'NSS.GHRR.' .ne. file2(POS_Dir2+1:POS_Dir2+9) )THEN
       WRITE(*,*)TRIM(file2)
       WRITE(*,*)file2(POS_Dir2+1:POS_Dir2+9)
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (NSS.)',&
            'check_overlap','extract_l1b_data.f90')
    ENDIF

    IF( file1(POS_Dir1+12:POS_Dir1+13) .ne. '.D' )THEN
       WRITE(*,*)TRIM(file1)
       WRITE(*,*)file1(POS_Dir1+12:POS_Dir1+13)
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (.D)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_Year1 = POS_Dir1+14
    ENDIF

    IF( file2(POS_Dir2+12:POS_Dir2+13) .ne. '.D' )THEN
       WRITE(*,*)TRIM(file2)
       WRITE(*,*)file2(POS_Dir2+12:POS_Dir2+13)
       CALL Gbcs_Critical(.TRUE.,'File2 not a NOAA Level 1B file (.D)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_Year2 = POS_Dir2+14
    ENDIF

    IF( file1(POS_Dir1+19:POS_Dir1+20) .ne. '.S' )THEN
       WRITE(*,*)TRIM(file1)
       WRITE(*,*)file1(POS_Dir1+19:POS_Dir1+20)
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (.S)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_Start1 = POS_Dir1+21
    ENDIF

    IF( file2(POS_Dir2+19:POS_Dir2+20) .ne. '.S' )THEN
       WRITE(*,*)TRIM(file2)
       WRITE(*,*)file2(POS_Dir2+19:POS_Dir2+20)
       CALL Gbcs_Critical(.TRUE.,'File2 not a NOAA Level 1B file (.S)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_Start2 = POS_Dir2+21
    ENDIF

    IF( file1(POS_Dir1+25:POS_Dir1+26) .ne. '.E' )THEN
       WRITE(*,*)TRIM(file1)
       WRITE(*,*)file1(POS_Dir1+25:POS_Dir1+26)
       CALL Gbcs_Critical(.TRUE.,'File1 not a NOAA Level 1B file (.E)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_End1 = POS_Dir1+27
    ENDIF

    IF( file2(POS_Dir2+25:POS_Dir2+26) .ne. '.E' )THEN
       WRITE(*,*)TRIM(file2)
       WRITE(*,*)file2(POS_Dir2+25:POS_Dir2+26)
       CALL Gbcs_Critical(.TRUE.,'File2 not a NOAA Level 1B file (.E)',&
            'check_overlap','extract_l1b_data.f90')
    ELSE
       POS_End2 = POS_Dir2+27
    ENDIF

    !
    ! Get start/end times from filename - convert to julian day to get round day/year 
    ! boundaries
    ! Need end time of file1 and start time of file2
    !
    CALL Get_Jul_Day(file1,POS_Year1,POS_Start1,POS_End1,jday1_start,jday1_end)
    CALL Get_Jul_Day(file2,POS_Year2,POS_Start2,POS_End2,jday2_start,jday2_end)

    check_overlap = .FALSE.
    IF( jday1_end .ge. jday2_start )THEN
       check_overlap = .TRUE.
    ELSE IF( ABS(jday1_end - jday2_start) .lt. 0.0007 )THEN
       !
       ! Check against a minute and a half (rounding errors)
       !
       check_overlap = .TRUE.
    ENDIF

    !
    ! Force check to true as check losing 10% files
    !
    check_overlap = .TRUE.

  END FUNCTION check_overlap

  SUBROUTINE Get_Jul_Day(file,POS_Year,POS_Start,POS_End,jday_start,jday_end)

    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: POS_Year
    INTEGER, INTENT(IN) :: POS_Start
    INTEGER, INTENT(IN) :: POS_End
    REAL(GbcsDble), INTENT(OUT) :: jday_start
    REAL(GbcsDble), INTENT(OUT) :: jday_end

    ! Local variables
    INTEGER :: IOS
    INTEGER :: year
    INTEGER :: dayno
    INTEGER :: month
    INTEGER :: day
    INTEGER :: hour
    INTEGER :: minute
    INTEGER :: hour_start
    INTEGER :: minute_start
    TYPE(DateTime) :: datestr

    READ(file(POS_Year:POS_Year+1),'(i2)',IOSTAT=IOS)year
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse year from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    READ(file(POS_Year+2:POS_Year+4),'(i3)',IOSTAT=IOS)dayno
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse dayno from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    CALL Day_Number_to_Date(year,dayno,month,day)
    READ(file(POS_Start:POS_Start+1),'(i2)',IOSTAT=IOS)hour
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse hour from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    READ(file(POS_Start+2:POS_Start+3),'(i2)',IOSTAT=IOS)minute
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse minute from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    datestr%year = year
    datestr%month = month
    datestr%day = day
    datestr%hour = hour
    datestr%minute = minute
    datestr%seconds = 0
    datestr%SEC1000 = 0
    jday_start = Date_to_JD(datestr)
    READ(file(POS_End:POS_End+1),'(i2)',IOSTAT=IOS)hour
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse hour from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    READ(file(POS_End+2:POS_End+3),'(i2)',IOSTAT=IOS)minute
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse minute from filename',&
            'Get_Jul_Day','extract_l1b_data.f90')
    ENDIF
    datestr%hour = hour
    datestr%minute = minute
    jday_end = Date_to_JD(datestr)

    !
    ! Check to see if we haven't gone over day boundary
    !
    IF( jday_end .lt. jday_start )THEN
       jday_end = jday_end+1.0D0
    ENDIF

  END SUBROUTINE Get_Jul_Day

  INTEGER FUNCTION get_instr(file)

    CHARACTER(LEN=*), INTENT(IN) :: file

    INTEGER :: POS_Dir
    CHARACTER(LEN=2) :: instr_str

    POS_Dir = INDEX(file,'/',.TRUE.)

    instr_str = file(POS_Dir+10:POS_Dir+11)

    SELECT CASE(instr_str)
    CASE('TN')
       get_instr = 1
    CASE('NA')
       get_instr = 6
    CASE('NC')
       get_instr = 7
    CASE('NE')
       get_instr = 8
    CASE('NF')
       get_instr = 9
    CASE('NG')
       get_instr = 10
    CASE('NH')
       get_instr = 11
    CASE('ND')
       get_instr = 12
    CASE('NJ')
       get_instr = 14
    CASE('NK')
       get_instr = 15
    CASE('NL')
       get_instr = 16
    CASE('NM')
       get_instr = 17
    CASE('NN')
       get_instr = 18
    CASE('NP')
       get_instr = 19
    CASE('M2')
       get_instr = -1
    CASE('M1')
       get_instr = -2
    END SELECT

  END FUNCTION get_instr

  SUBROUTINE read_all_data(mc_harm,&
       nFile,file1,file2,file3,file4,file5,uuid_in,&
       AVHRRout,year1,month1,day1,hour1,minute1,year2,month2,day2,&
       hour2,minute2,output_filename,walton_cal,split_single_file,&
       pygac1,pygac2,pygac3,pygac4,pygac5,gbcs_l1c_output,gbcs_l1c_cal,&
       walton_only,keep_temp,write_fcdr,output_solar_temp,montecarlo,&
       ocean_only)

    USE normal_generator

    TYPE(mc_harm_str), INTENT(IN) :: mc_harm
    INTEGER, INTENT(IN) :: nFile
    CHARACTER(LEN=*), INTENT(IN) :: file1
    CHARACTER(LEN=*), INTENT(IN) :: file2
    CHARACTER(LEN=*), INTENT(IN) :: file3
    CHARACTER(LEN=*), INTENT(IN) :: file4
    CHARACTER(LEN=*), INTENT(IN) :: file5
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    TYPE(AVHRR_Data), INTENT(OUT), TARGET :: AVHRRout
    INTEGER, INTENT(IN) :: year1
    INTEGER, INTENT(IN) :: month1
    INTEGER, INTENT(IN) :: day1
    INTEGER, INTENT(IN) :: hour1
    INTEGER, INTENT(IN) :: minute1
    INTEGER, INTENT(IN) :: year2
    INTEGER, INTENT(IN) :: month2
    INTEGER, INTENT(IN) :: day2
    INTEGER, INTENT(IN) :: hour2
    INTEGER, INTENT(IN) :: minute2
    CHARACTER(LEN=*), INTENT(IN) :: output_filename
    LOGICAL, INTENT(IN) :: walton_cal
    LOGICAL, INTENT(IN) :: split_single_file
    CHARACTER(LEN=*), INTENT(IN) :: pygac1
    CHARACTER(LEN=*), INTENT(IN) :: pygac2
    CHARACTER(LEN=*), INTENT(IN) :: pygac3
    CHARACTER(LEN=*), INTENT(IN) :: pygac4
    CHARACTER(LEN=*), INTENT(IN) :: pygac5
    LOGICAL, INTENT(IN), OPTIONAL :: gbcs_l1c_output
    LOGICAL, INTENT(IN), OPTIONAL :: gbcs_l1c_cal
    LOGICAL, INTENT(IN), OPTIONAL :: walton_only
    LOGICAL, INTENT(IN), OPTIONAL :: keep_temp
    LOGICAL, INTENT(IN), OPTIONAL :: write_fcdr
    INTEGER, INTENT(IN), OPTIONAL :: output_solar_temp
    LOGICAL, INTENT(IN), OPTIONAL :: montecarlo
    LOGICAL, INTENT(IN), OPTIONAL :: ocean_only

    ! Local variables
    TYPE(Imagery) :: IMG
    TYPE(AVHRR_Data), TARGET :: AVHRR
    TYPE(AVHRR_Data), TARGET :: AVHRRtmp
    TYPE(AVHRR_Data), TARGET :: AVHRR_MC
    TYPE(AVHRR_Data), TARGET :: AVHRR_Total
    TYPE(AVHRR_Instrument_Coefs) :: instr_coefs
    INTEGER :: start_pos, stop_pos
    TYPE(AVHRR_Data), POINTER :: pAVHRR=>NULL()
    TYPE(AVHRR_Data), POINTER :: pAVHRR_MC=>NULL()
    LOGICAL :: out_radiances = .FALSE.
    LOGICAL :: trim_data
    LOGICAL :: correctict
    LOGICAL :: applyict
    LOGICAL :: walton_biascorr
    LOGICAL :: apply_scenet_bias
    INTEGER :: trim_low
    INTEGER :: trim_high
    INTEGER :: I
    TYPE(Walton_Struct) :: walton_str
    LOGICAL :: walt_only
    INTEGER :: output_solartemp
    LOGICAL :: use_montecarlo
    TYPE(mc_delta_str) :: delta_rad
    INTEGER :: seedval
    LOGICAL :: oceanonly
    INTEGER :: start_valid
    INTEGER :: stop_valid

    IF( PRESENT(walton_only) )THEN
       walt_only = walton_only
    ELSE
       walt_only = .FALSE.
    ENDIF
    IF( PRESENT(output_solar_temp) )THEN
       output_solartemp = output_solar_temp
    ELSE
       output_solartemp = -1
    ENDIF
    IF( PRESENT(montecarlo) )THEN
       use_montecarlo = montecarlo
    ELSE
       use_montecarlo = .FALSE.
    ENDIF

    IF( use_montecarlo )THEN
       CALL Set_Random_Seed(seedval)
    ELSE
       seedval = -1
    ENDIF

    IF( PRESENT(ocean_only) )THEN
       oceanonly = ocean_only
    ELSE
       oceanonly = .FALSE.
    ENDIF

    !
    ! Setup GbcsDataPath
    !
    IMG%GbcsDataPath = '/group_workspaces/jasmin2/nceo_uor/users/jmittaz/AVHRR/&
         &code/git/gbcs_new_calibration/dat_cci/'

    AVHRR%newCalibration_there = .FALSE.
    AVHRR_Total%newCalibration_there = .FALSE.
!    instr = get_instr(file1)
    AVHRR_Total%arraySize = 0
    AVHRR_Total%norig_l1b = 0
    IF( nFile .eq. 1 )THEN
       CALL read_file(file1,AVHRR_Total,uuid_in,pygac1,1,&
            out_instr_coefs=instr_coefs)
       IF( .not. AVHRR_Total%valid_data_there )THEN
          RETURN
       ENDIF
       AVHRR_Total%orig_l1b(1) = file1(1:42)
       AVHRR_Total%norig_l1b = 1
    ELSE IF( nFile .eq. 2 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,pygac1,1,&
            out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          AVHRR_Total%orig_l1b(1) = file1(1:42)
          AVHRR_Total%norig_l1b = 1
          CALL read_file(file2,AVHRR,uuid_in,pygac2,2)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(2) = file2(1:42)
             AVHRR_Total%norig_l1b = 2
          ENDIF
       ELSE
          RETURN
       ENDIF
       CALL Deallocate_OutData(AVHRR)
    ELSE IF( nFile .eq. 3 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,pygac1,1,&
            out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          AVHRR_Total%orig_l1b(1) = file1(1:42)
          AVHRR_Total%norig_l1b = 1
          CALL read_file(file2,AVHRR,uuid_in,pygac2,2)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(2) = file2(1:42)
             AVHRR_Total%norig_l1b = 2
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
               'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in,pygac3,3)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(3) = file3(1:42)
             AVHRR_Total%norig_l1b = 3
          ENDIF
       ENDIF
       CALL Deallocate_OutData(AVHRR)
    ELSE IF( nfile .eq. 4 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,pygac1,1,&
            out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          AVHRR_Total%orig_l1b(1) = file1(1:42)
          AVHRR_Total%norig_l1b = 1
          CALL read_file(file2,AVHRR,uuid_in,pygac2,2)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(2) = file2(1:42)
             AVHRR_Total%norig_l1b = 2
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in,pygac3,3)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(3) = file3(1:42)
             AVHRR_Total%norig_l1b = 3
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file3,file4) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file3 and file4',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file4,AVHRR,uuid_in,pygac4,4)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(4) = file4(1:42)
             AVHRR_Total%norig_l1b = 4
          ENDIF
          CALL Deallocate_OutData(AVHRR)
       ENDIF
    ELSE IF( nfile .eq. 5 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,pygac1,1,&
            out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          AVHRR_Total%orig_l1b(1) = file1(1:42)
          AVHRR_Total%norig_l1b = 1
          CALL read_file(file2,AVHRR,uuid_in,pygac2,2)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(2) = file2(1:42)
             AVHRR_Total%norig_l1b = 2
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in,pygac3,3)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(3) = file3(1:42)
             AVHRR_Total%norig_l1b = 3
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file3,file4) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file3 and file4',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file4,AVHRR,uuid_in,pygac4,4)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(4) = file4(1:42)
             AVHRR_Total%norig_l1b = 4
          ENDIF
          CALL Deallocate_OutData(AVHRR)
!MT: 20-10-2017: 'between file3 and file4' --> 'between file4 and file5'
          IF( .not. check_overlap(file4,file5) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file4 and file5',&
                  'read_all_data','extract_l1b_data.f90') 
          ENDIF
          CALL read_file(file5,AVHRR,uuid_in,pygac5,5)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
             AVHRR_Total%orig_l1b(5) = file5(1:42)
             AVHRR_Total%norig_l1b = 5
          ENDIF
          CALL Deallocate_OutData(AVHRR)
       ENDIF
    ENDIF

    IF( AVHRR_Total%arraySize .le. 0 )THEN
       CALL Gbcs_Critical(.TRUE.,'No data available',&
            'read_all_data','combine_orbits.f90')
    ENDIF
    trim_data = .FALSE.
    IF( nfile .gt. 1 )THEN
       !
       ! Resize to a single orbit which is the assumption of the recalibrate
       ! routines
       !
       CALL Resize_Orbit_Equator(AVHRR_Total,AVHRR,start_pos,stop_pos,&
            year1,month1,day1,hour1,minute1,&
            year1,month2,day2,hour2,minute2,&
            trim_data,trim_low,trim_high)
       CALL Deallocate_OutData(AVHRR_Total)
       IF( AVHRR%arraySize .le. 0 )THEN
          CALL Gbcs_Critical(.TRUE.,'No data available',&
               'read_all_data','combine_orbits.f90')
       ENDIF

       !
       ! Copy top part of AVHRR arrays
       !
       IF( use_montecarlo )THEN
          pAVHRR_MC => AVHRR_MC
          CALL Allocate_OldCal(AVHRR,pAVHRR_MC,1,AVHRR%arraySize)
          CALL Copy_Top_Data(AVHRR,AVHRR_MC,1,AVHRR%arraySize)
          DO I=1,AVHRR%arraySize
             CALL Copy_All_Scan(AVHRR,AVHRR_MC,I,I,.FALSE.)
          END DO
       ENDIF

       !
       ! Make sure radiances are output as Marines code expects this
       !
       IMG%GbcsDataPath = '/group_workspaces/jasmin2/nceo_uor/users/jmittaz/&
            &AVHRR/code/git/gbcs_new_calibration/dat_cci/'
       IF( walton_cal )THEN
          CALL Setup_Walton( IMG, AVHRR, walton_str, srfonly=.TRUE.)
       ENDIF
       pAVHRR => AVHRR
       IF( walton_cal )THEN
          IF( walt_only )THEN
             WRITE(*,'('' ====== Walton et al. only calibration ====='')')
             correctict = .FALSE.
             applyict = .FALSE.
             walton_biascorr = .FALSE.
             apply_scenet_bias = .FALSE.
          ELSE
             WRITE(*,'('' ====== CORRECTED Walton et al. only calibration ====='')')
             correctict = .TRUE.
             applyict = .TRUE.
             walton_biascorr = .TRUE.
             apply_scenet_bias = .TRUE.
          ENDIF
       ELSE
          correctict = .TRUE.
          applyict = .TRUE.
          walton_biascorr = .FALSE.
          apply_scenet_bias = .FALSE.
       ENDIF
       IF( 102 .eq. output_solartemp )THEN
          CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,.TRUE.,walton_str,&
               moon_events=.TRUE.,use_old_solar=.FALSE.,&
               correct_solar_simple=.FALSE.,new_vis_cal=.TRUE.,&
               noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
               dig_noise=.TRUE.,all_noise=.TRUE.,&
               correctict=.TRUE.,applyict=.FALSE.,walton=walton_cal,&
               tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
               bad_tiny=.TRUE.,noise_ict=.TRUE.)
       ELSE IF( 101 .eq. output_solartemp )THEN
          CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,.TRUE.,walton_str,&
               moon_events=.TRUE.,use_old_solar=.FALSE.,&
               correct_solar_simple=.FALSE.,new_vis_cal=.TRUE.,&
               noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
               dig_noise=.TRUE.,all_noise=.TRUE.,&
               correctict=.FALSE.,applyict=.FALSE.,walton=walton_cal,&
               tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
               bad_tiny=.TRUE.,noise_ict=.TRUE.)
       ELSE IF( 100 .eq. output_solartemp )THEN 
          CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,.TRUE.,walton_str,&
               moon_events=.TRUE.,use_old_solar=.FALSE.,&
               correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
               noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
               dig_noise=.TRUE.,all_noise=.TRUE.,&
               correctict=correctict,applyict=applyict,walton=walton_cal,&
               tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
               bad_tiny=.TRUE.,noise_ict=.TRUE.)
       ELSE
          CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,.TRUE.,walton_str,&
               moon_events=.TRUE.,use_old_solar=.FALSE.,&
               correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
               noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
               dig_noise=.TRUE.,all_noise=.TRUE.,&
               correctict=correctict,applyict=applyict,walton=walton_cal,&
               tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
               bad_tiny=.TRUE.,noise_ict=.TRUE.)
       ENDIF

       IF( walton_cal )THEN
          CALL Setup_Walton( IMG, AVHRR, walton_str,&
               scenet_all_instr=.TRUE.,&
               apply_scenet_bias=apply_scenet_bias)
          CALL Calibration_Walton(IMG,AVHRR,.TRUE.,1.,&
               walton_str,walton_biascorr,applyict,ict_tinstr=.FALSE.,&
               apply_scenet_bias=apply_scenet_bias)
       ENDIF

       IF( use_montecarlo )THEN
          IF( walton_cal )THEN
             CALL Gbcs_Critical(.TRUE.,'Cannot do MonteCarlo with Walton',&
                  'read_all_data','combine_orbits.f90')
          ENDIF
          IF( oceanonly )THEN
             ! copy to IMG structure
             CALL Copy_To_IMG( instr_coefs, pAVHRR, IMG, .TRUE., &
                  .FALSE., start_valid, stop_valid, &
                  .FALSE.)
          ENDIF
          CALL Run_MonteCarlo(mc_harm,IMG,instr_coefs,&
               pAVHRR,AVHRR_MC,.TRUE.,walton_str,&
               delta_rad,moon_events=.TRUE.,use_old_solar=.FALSE.,&
               correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
               noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
               dig_noise=.TRUE.,all_noise=.TRUE.,&
               correctict=correctict,applyict=applyict,walton=walton_cal,&
               tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
               bad_tiny=.TRUE.,ocean_only=oceanonly)

       ENDIF

       !
       ! Resize to output
       !
       CALL Resize_Orbit(AVHRR,AVHRRtmp,start_pos,stop_pos,all=.TRUE.,&
            montecarlo=use_montecarlo,delta_rad=delta_rad)

       CALL fill_missing_lines(AVHRRtmp,AVHRRout,&
            montecarlo=use_montecarlo,delta_rad=delta_rad)

       CALL Deallocate_OutData(AVHRR)
       CALL Deallocate_OutData(AVHRRtmp)

!       CALL check_latlon(AVHRRout,47.609,-29.85500,usenew=.TRUE.)
       !
       ! Add in FIDUCEO uncertainties
       !
       CALL Add_FIDUCEO_Uncert(IMG,AVHRRout,uuid_in,output_filename,&
            gbcs_l1c_output=gbcs_l1c_output,&
            gbcs_l1c_cal=gbcs_l1c_cal,use_walton=walton_cal,&
            keep_temp=keep_temp,write_fcdr=write_fcdr,&
            monte_carlo=use_montecarlo,&
            delta_radiance=delta_rad,seedval=seedval,ocean_only=oceanonly)
       CALL Deallocate_OutData(AVHRRout)    
    ELSE
       !
       ! Calibrate whole orbit - it gets trimmed later
       !
       !
       ! Copy top part of AVHRR arrays for monte-carlo
       !
       IF( use_montecarlo )THEN
          pAVHRR_MC => AVHRR_MC
          CALL Allocate_OldCal(AVHRR,pAVHRR_MC,1,AVHRR%arraySize)
          CALL Copy_Top_Data(AVHRR,AVHRR_MC,1,AVHRR%arraySize)
          DO I=1,AVHRR%arraySize
             CALL Copy_All_Scan(AVHRR,AVHRR_MC,I,I,.FALSE.)
          END DO
       ENDIF
       !
       ! Make sure radiances are output as Marines code expects this
       !
       IF( walton_cal )THEN
          CALL Setup_Walton( IMG, AVHRR, walton_str, srfonly=.TRUE.)
       ENDIF
       pAVHRR => AVHRR_Total
       IF( walt_only )THEN
          correctict = .FALSE.
          applyict = .FALSE.
          walton_biascorr = .FALSE.
          apply_scenet_bias = .FALSE.
       ELSE
          correctict = .TRUE.
          applyict = .TRUE.
          walton_biascorr = .TRUE.
          apply_scenet_bias = .TRUE.
       ENDIF
       CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,.TRUE.,walton_str,&
            moon_events=.TRUE.,&
            correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
            noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
            dig_noise=.TRUE.,all_noise=.TRUE.,&
            correctict=correctict,applyict=applyict,walton=walton_cal,&
            tinstrapply=.FALSE.,noise_ict=.TRUE.)

       IF( walton_cal )THEN
          CALL Setup_Walton( IMG, AVHRR, walton_str,&
               scenet_all_instr=.TRUE.,&
               apply_scenet_bias=apply_scenet_bias)
          CALL Calibration_Walton(IMG,AVHRR,.TRUE.,1.,&
               walton_str,walton_biascorr,applyict,ict_tinstr=.FALSE.,&
               apply_scenet_bias=apply_scenet_bias)
       ENDIF
       IF( split_single_file )THEN
          IF( use_montecarlo )THEN
             IF( walton_cal )THEN
                CALL Gbcs_Critical(.TRUE.,'Cannot do MonteCarlo with Walton',&
                     'read_all_data','combine_orbits.f90')
             ENDIF
             CALL Run_MonteCarlo(mc_harm,IMG,instr_coefs,&
                  pAVHRR,AVHRR_MC,.TRUE.,walton_str,&
                  delta_rad,moon_events=.TRUE.,use_old_solar=.FALSE.,&
                  correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
                  noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
                  dig_noise=.TRUE.,all_noise=.TRUE.,&
                  correctict=correctict,applyict=applyict,walton=walton_cal,&
                  tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
                  bad_tiny=.TRUE.,ocean_only=oceanonly)
          ENDIF

          !
          ! Now split out sections of data
          !
          CALL Resize_Orbit_Equator(AVHRR_Total,AVHRRout,start_pos,stop_pos,&
               year1,month1,day1,hour1,minute1,&
               year1,month2,day2,hour2,minute2,&
               trim_data,trim_low,trim_high,&       
               montecarlo=use_montecarlo,delta_rad=delta_rad)
          CALL Resize_Orbit(AVHRRout,AVHRRtmp,start_pos,stop_pos,all=.TRUE.,&
               montecarlo=use_montecarlo,delta_rad=delta_rad)
          CALL fill_missing_lines(AVHRRtmp,AVHRR,&
               montecarlo=use_montecarlo,delta_rad=delta_rad)
          CALL Deallocate_OutData(AVHRRout)
          CALL Deallocate_OutData(AVHRRtmp)

          !
          ! Add in FIDUCEO uncertainties
          !
          CALL Add_FIDUCEO_Uncert(IMG,AVHRR,uuid_in,output_filename,&
               gbcs_l1c_output=gbcs_l1c_output,&
               gbcs_l1c_cal=gbcs_l1c_cal,use_walton=walton_cal,&
               keep_temp=keep_temp,write_fcdr=write_fcdr,&
               monte_carlo=use_montecarlo,&
               delta_radiance=delta_rad,seedval=seedval,ocean_only=oceanonly)
          CALL Deallocate_OutData(AVHRR)
       ELSE
          !
          ! Monte Carlo
          !
          IF( use_montecarlo )THEN
             IF( walton_cal )THEN
                CALL Gbcs_Critical(.TRUE.,'Cannot do MonteCarlo with Walton',&
                     'read_all_data','combine_orbits.f90')
             ENDIF
             pAVHRR => AVHRR_Total
             CALL Run_MonteCarlo(mc_harm,IMG,instr_coefs,&
                  pAVHRR,AVHRR_MC,.TRUE.,walton_str,&
                  delta_rad,moon_events=.TRUE.,use_old_solar=.FALSE.,&
                  correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
                  noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
                  dig_noise=.TRUE.,all_noise=.TRUE.,&
                  correctict=correctict,applyict=applyict,walton=walton_cal,&
                  tinstrapply=.FALSE.,output_solar_temp=output_solartemp,&
                  bad_tiny=.TRUE.,ocean_only=oceanonly)
          ENDIF
          !
          ! Output complete orbit with no equator splitting
          !
          CALL fill_missing_lines(AVHRR_Total,AVHRRout,&
               montecarlo=use_montecarlo,delta_rad=delta_rad)
          CALL Deallocate_OutData(AVHRR_Total)

          !
          ! Add in FIDUCEO uncertainties
          !
          CALL Add_FIDUCEO_Uncert(IMG,AVHRRout,uuid_in,output_filename,&
               gbcs_l1c_output=gbcs_l1c_output,&
               gbcs_l1c_cal=gbcs_l1c_cal,use_walton=walton_cal,&
               keep_temp=keep_temp,write_fcdr=write_fcdr,&
               monte_carlo=use_montecarlo,&
               delta_radiance=delta_rad,seedval=seedval,ocean_only=oceanonly)
          CALL Deallocate_OutData(AVHRRout)
       ENDIF
    ENDIF

    IF( use_montecarlo )THEN
       CALL Deallocate_OutData(AVHRR_MC)
    ENDIF

  END SUBROUTINE read_all_data

  LOGICAL FUNCTION get_equ_cross_type(AVHRR)RESULT(ascending)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR

    ascending = .TRUE.
    IF( AVHRR%AVHRR_No .eq. 1 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 6 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 7 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 8 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 9 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 10 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 11 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 12 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 14 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 15 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 16 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 17 )THEN
       ascending = .FALSE.
    ELSE IF( AVHRR%AVHRR_No .eq. 18 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. 19 )THEN
       ascending = .TRUE.
    ELSE IF( AVHRR%AVHRR_No .eq. -1 )THEN
       ascending = .FALSE.
    ELSE
       CALL Gbcs_Critical(.TRUE.,'Cannot match AVHRR No',&
            'get_equ_cross_type','NOAA_LoadAVHRRLevel1B.f90')
    ENDIF

    IF( ascending )THEN
       WRITE(*,*)'Equator crossing type (daytime) : ascending'
    ELSE
       WRITE(*,*)'Equator crossing type (daytime) : descending'
    ENDIF
    RETURN

  END FUNCTION get_equ_cross_type

  LOGICAL FUNCTION equator_crossing_test(lat1,lat2,ascending_type)RESULT(ok)

    REAL, INTENT(IN) :: lat1
    REAL, INTENT(IN) :: lat2
    LOGICAL, INTENT(IN) :: ascending_type

    ok = .FALSE.
    IF( ascending_type )THEN
       IF( lat1 .le. 0. .and. lat2 .ge. 0. )THEN
          ok = .TRUE.
       ENDIF
    ELSE
       IF( lat1 .ge. 0. .and. lat2 .le. 0. )THEN
          ok = .TRUE.
       ENDIF
    ENDIF

  END FUNCTION equator_crossing_test

  SUBROUTINE Resize_Orbit_Equator(AVHRR,AVHRRout,&
       out_start_pos,out_stop_pos,&
       year1,month1,day1,hour1,minute1,&
       year2,month2,day2,hour2,minute2,&
       trim_data,trim_low,trim_high,nosmooth,&
       montecarlo,delta_rad)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(AVHRR_Data), INTENT(OUT) :: AVHRRout
    INTEGER, INTENT(OUT) :: out_start_pos
    INTEGER, INTENT(OUT) :: out_stop_pos
    INTEGER, INTENT(IN) :: year1
    INTEGER, INTENT(IN) :: month1
    INTEGER, INTENT(IN) :: day1
    INTEGER, INTENT(IN) :: hour1
    INTEGER, INTENT(IN) :: minute1
    INTEGER, INTENT(IN) :: year2
    INTEGER, INTENT(IN) :: month2
    INTEGER, INTENT(IN) :: day2
    INTEGER, INTENT(IN) :: hour2
    INTEGER, INTENT(IN) :: minute2
    LOGICAL, INTENT(OUT) :: trim_data
    INTEGER, INTENT(OUT) :: trim_low
    INTEGER, INTENT(OUT) :: trim_high
    LOGICAL, INTENT(IN), OPTIONAL :: nosmooth
    LOGICAL, INTENT(IN), OPTIONAL :: montecarlo
    TYPE(mc_delta_str), INTENT(INOUT), OPTIONAL :: delta_rad

    ! Local variables
    INTEGER :: I,J
    INTEGER :: start_second
    INTEGER :: centre
    INTEGER :: start_pos
    INTEGER :: stop_pos
    INTEGER :: first_equ
    INTEGER :: last_equ
    INTEGER :: nlines
    LOGICAL :: no_smooth
    LOGICAL :: upwards
    LOGICAL :: first_seg
    INTEGER :: dayno
    REAL :: hours
    REAL(GbcsDble) :: time1
    REAL(GbcsDble) :: time2
    INTEGER :: diff
    LOGICAL :: make_orbit1
    LOGICAL :: make_orbit2
    LOGICAL :: ascending_type
    LOGICAL :: monte_carlo
    TYPE(mc_delta_str) :: delta_rad_tmp

    IF( PRESENT(nosmooth) )THEN
       no_smooth = nosmooth
    ELSE
       no_smooth = .FALSE.
    ENDIF

    IF( PRESENT(montecarlo) )THEN
       IF( .not. PRESENT(delta_rad) )THEN
          CALL Gbcs_Critical(.TRUE.,'delta_rad required',&
               'Resize_Orbit_Equator','combine_orbits.f90')
       ENDIF
       monte_carlo = montecarlo
    ELSE
       monte_carlo = .FALSE.
    ENDIF

    AVHRRout%arraySize=00
    !
    ! From Satellite No work out if need ascending or descending
    !
    ascending_type = get_equ_cross_type(AVHRR)

    ! Get times of predicted equation crossing times
    dayno = Day_Of_Year(year1,month1,day1)
    hours = hour1 + minute1/60.
    time1 = get_date_time(year1,dayno,hours,1975)
    dayno = Day_Of_Year(year2,month2,day2)
    hours = hour2 + minute2/60.
    time2 = get_date_time(year2,dayno,hours,1975)

    ! Find scanlines at time1/time2 and see if equator crossing is there
    make_orbit1=.FALSE.
    make_orbit2=.FALSE.
    upwards = .TRUE.
    first_equ = -1
    last_equ = -1
    centre = AVHRR%nelem/2
    Loop2: DO I=1,AVHRR%arraySize
       IF( 120. .gt. ABS(AVHRR%time(I)-time1) )THEN
          Loop1: DO J=I,AVHRR%arraySize
             IF( equator_crossing_test(AVHRR%lat(centre,I-1),&
                  AVHRR%lat(centre,I),ascending_type) )THEN
!             IF( AVHRR%lat(centre,I-1) .le. 0. .and. &
!                  AVHRR%lat(centre,I) .ge. 0. )THEN
                first_equ=I
                upwards = .TRUE.
                EXIT Loop2
             ENDIF
          END DO Loop1
       ENDIF
    END DO Loop2
    IF(first_equ .eq. -1)THEN
       checkLoop1: DO I=1,AVHRR%arraySize
          IF( AVHRR%time(I) .gt. 0 )THEN
             IF( AVHRR%time(I) .gt. time1 )THEN
                first_equ=I
                IF( AVHRR%time(1) .gt. time1 )THEN
                   make_orbit1=.TRUE.
                ENDIF
                EXIT checkLoop1
             ENDIF
          ENDIF
       END DO checkLoop1
       IF( first_equ .eq. -1 )THEN
          CALL Gbcs_Critical(.TRUE.,&
               'Cannot find equator where it should be (first)',&
               'Resize orbit equator','combine_orbits.f90')
       ENDIF
    ENDIF

    !
    ! Find equator in other direction
    !
    Loop3: DO I=first_equ+1,AVHRR%arraySize
       IF( 120. .gt. ABS(AVHRR%time(I)-time2) )THEN
          IF( equator_crossing_test(AVHRR%lat(centre,I-1),&
               AVHRR%lat(centre,I),ascending_type) )THEN
!          IF( AVHRR%lat(centre,I-1) .le. 0. .and. &
!               AVHRR%lat(centre,I) .ge. 0. )THEN
             last_equ=I-1
             EXIT Loop3
          ENDIF
       ENDIF
    END DO Loop3
    IF(last_equ .eq. -1)THEN
       checkLoop2: DO I=AVHRR%arraySize,first_equ,-1
          IF( AVHRR%time(I) .gt. 0 )THEN
             IF( AVHRR%time(I) .lt. time2 )THEN
                last_equ=I
                IF( AVHRR%time(AVHRR%arraySize) .lt. time2 )THEN
                   make_orbit2=.TRUE.
                ENDIF
                EXIT checkLoop2
             ENDIF
          ENDIF
       END DO checkLoop2
       IF( last_equ .eq. -1 )THEN
          CALL Gbcs_Critical(.TRUE.,&
               'Cannot find equator where it should be (first)',&
               'Resize orbit equator','combine_orbits.f90')
       ENDIF
    ENDIF

    IF( make_orbit1 )THEN
       CALL Gbcs_Warning(.TRUE.,'Actual first equator crossing not found',&
               'Resize orbit equator','combine_orbits.f90')            
    ENDIF
    IF( make_orbit2 )THEN
       CALL Gbcs_Warning(.TRUE.,'Actual second equator crossing not found',&
               'Resize orbit equator','combine_orbits.f90')            
    ENDIF

    ! Make choices on the basis of coverage. We have
    !
    !     not make_orbit1/not make_orbit2
    !        - Normal case where we are going equator to equator
    !     make_orbit1/not make_orbit2
    !        - found second equator crossing but data does not include
    !          first one. Need to make at least complete orbit and then
    !          trim later
    !     not make_orbit1/make_orbit2
    !        - found first equator crossing but data does not include
    !          second one. Need to make at least complete orbit (backward
    !          in time) and then trim later
    !     make_orbit1/make_orbit2
    !        - didn't find either - basically have from 1,arraySize 
    !          and no trimming
    IF( .not. make_orbit1 .and. .not. make_orbit2 )THEN
       trim_data = .FALSE.
    ELSE IF( make_orbit1 .and. .not. make_orbit2 )THEN
       trim_data = .TRUE.
       trim_low = 1
       trim_high = last_equ
       IF( first_equ + 15000 .gt. last_equ )THEN
          IF( first_equ + 15000 .gt. AVHRR%arraySize )THEN
             last_equ = AVHRR%arraySize
          ELSE
             last_equ = first_equ + 15000
          ENDIF
       ENDIF
!MT: 24-10-2017: fix reversed logic in Resize_orbit_equator
!    ELSE IF( make_orbit1 .and. .not. make_orbit2 )THEN
    ELSE IF( .not. make_orbit1 .and. make_orbit2 )THEN
       trim_data = .TRUE.
       trim_low = first_equ
       trim_high = last_equ
       IF( last_equ - 15000 .lt. first_equ )THEN
          IF( last_equ-15000 .lt. 1 )THEN
             first_equ = 1             
          ELSE
             trim_low = first_equ - (last_equ-15000)
             first_equ = last_equ - 15000
          ENDIF
       ENDIF
    ELSE 
       ! Only have section - all we have
!MT: 11-12-2017: reset to FALSE
       trim_data = .FALSE.
    ENDIF

!    IF( -1 .eq. first_equ )THEN
!       start_second = 1
!    ELSE
!       start_second = first_equ+1       
!    ENDIF
!    IF( start_second .lt. AVHRR%arraySize )THEN
!       Loop2: DO I=start_second+1,AVHRR%arraySize
!          IF( upwards )THEN
!             IF( AVHRR%lat(centre,I-1) .le. 0. .and. &
!                  AVHRR%lat(centre,I) .ge. 0. )THEN
!                last_equ=I
!                EXIT Loop2
!             ENDIF
!          ELSE
!             IF( AVHRR%lat(centre,I-1) .ge. 0. .and. &
!                  AVHRR%lat(centre,I) .le. 0. )THEN
!                last_equ=I-1
!                EXIT Loop2
!             ENDIF
!          ENDIF
!       END DO Loop2
!    ENDIF
!
!    IF( -1 .eq. first_equ .and. -1 .eq. last_equ )THEN
!       !
!       ! Case where we don't have the a full orbit and are within equator boun!ds
!       ! Just use everything
!       !
!       start_pos = 1
!       stop_pos = AVHRR%arraySize
!       out_start_pos = start_pos
!       out_stop_pos = stop_pos
!    ELSE IF( -1 .eq. last_equ )THEN
!       !
!       ! Again use all data for recalibration (get solar contamination)
!       !
!       start_pos = 1
!       stop_pos = AVHRR%arraySize
!       !
!       ! Output part depends on if you want first segment or last segment
!       !       
!       IF( first_seg )THEN
!          out_start_pos = 1
!          out_stop_pos = first_equ
!       ELSE
!          out_start_pos = first_equ
!          out_stop_pos = AVHRR%arraySize
!       ENDIF
!    ELSE
!       start_pos = first_equ
!       stop_pos = last_equ
!    ENDIF
    start_pos = first_equ
    stop_pos = last_equ

    IF( .not. no_smooth )THEN
       !
       ! Shift start/end to take into account smoothing
       !
       nlines = (stop_pos-start_pos) + 1
       out_start_pos = start_pos
       start_pos = MAX(1,start_pos-NPIXEL_PRT_SMOOTH)
       out_start_pos = out_start_pos-start_pos+1
       out_stop_pos = stop_pos
       stop_pos = MIN(AVHRR%arraySize,stop_pos+NPIXEL_PRT_SMOOTH)
       out_stop_pos = out_start_pos+nlines
!       if( trim_data )then
!          out_start_pos = trim_low
!          out_stop_pos = trim_high
!       endif
    ELSE
       out_start_pos = start_pos
       out_stop_pos = stop_pos
    ENDIF

    !
    ! Now resize to do recalibration
    !
    CALL Resize_Orbit(AVHRR,AVHRRout,start_pos,stop_pos)
    IF( monte_carlo )THEN
       CALL Copy_Delta_Rad(delta_rad,delta_rad_tmp,startpos=start_pos,&
            endpos=stop_pos)
       CALL Deallocate_DeltaRad(delta_rad)
       CALL Copy_Delta_Rad(delta_rad_tmp,delta_rad)
       CALL Deallocate_DeltaRad(delta_rad_tmp)
    ENDIF

  END SUBROUTINE Resize_Orbit_Equator

  SUBROUTINE Copy_DeltaRad_Ind(delta_rad_tmp,delta_rad,I,startpos,K)

    TYPE(mc_delta_str), INTENT(INOUT) :: delta_rad_tmp
    TYPE(mc_delta_str), INTENT(IN) :: delta_rad
    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(IN), OPTIONAL :: startpos
    INTEGER, INTENT(IN), OPTIONAL :: K

    ! Local variables
    INTEGER :: J

    IF( PRESENT(startpos) )THEN
       DO J=1,delta_rad%nelem
          delta_rad_tmp%ch1(J,I-startpos+1,:) = delta_rad%ch1(J,I,:)
          delta_rad_tmp%ch2(J,I-startpos+1,:) = delta_rad%ch2(J,I,:)
          delta_rad_tmp%ch3a(J,I-startpos+1,:) = delta_rad%ch3a(J,I,:)
          delta_rad_tmp%ch3(J,I-startpos+1,:) = delta_rad%ch3(J,I,:)
          delta_rad_tmp%ch4(J,I-startpos+1,:) = delta_rad%ch4(J,I,:)
          delta_rad_tmp%ch5(J,I-startpos+1,:) = delta_rad%ch5(J,I,:)
       END DO
    ELSE IF( PRESENT(K) )THEN
       DO J=1,delta_rad%nelem
          delta_rad_tmp%ch1(J,K,:) = delta_rad%ch1(J,I,:)
          delta_rad_tmp%ch2(J,K,:) = delta_rad%ch2(J,I,:)
          delta_rad_tmp%ch3a(J,K,:) = delta_rad%ch3a(J,I,:)
          delta_rad_tmp%ch3(J,K,:) = delta_rad%ch3(J,I,:)
          delta_rad_tmp%ch4(J,K,:) = delta_rad%ch4(J,I,:)
          delta_rad_tmp%ch5(J,K,:) = delta_rad%ch5(J,I,:)
       END DO
    ENDIF

  END SUBROUTINE Copy_DeltaRad_Ind

  SUBROUTINE Copy_Delta_Rad(delta_rad,delta_rad_tmp,&
       startpos,endpos)

    TYPE(mc_delta_str), INTENT(IN) :: delta_rad
    TYPE(mc_delta_str), INTENT(OUT) :: delta_rad_tmp
    INTEGER, INTENT(IN), OPTIONAL :: startpos
    INTEGER, INTENT(IN), OPTIONAL :: endpos

    ! Local variables
    INTEGER :: I
    INTEGER :: J
    INTEGER :: ndata
    LOGICAL :: resize

    IF( PRESENT(startpos) )THEN
       IF( .not. PRESENT(endpos) )THEN
          CALL Gbcs_Critical(.TRUE.,'Require endpos','Copy_Delta_Rad',&
               'combine_orbits.f90')
       ENDIF
       resize = .TRUE.
    ELSE IF( PRESENT(endpos) )THEN
       CALL Gbcs_Critical(.TRUE.,'Require endpos','Copy_Delta_Rad',&
            'combine_orbits.f90')
    ELSE
       resize = .FALSE.
    ENDIF

    IF( resize )THEN
       ndata = endpos-startpos+1
       CALL Allocate_DeltaRad(delta_rad%nelem,ndata,delta_rad%nMC,&
            delta_rad_tmp)
       DO I=startpos,endpos
          CALL Copy_DeltaRad_Ind(delta_rad_tmp,delta_rad,I,startpos=startpos)
       END DO
    ELSE
       CALL Allocate_DeltaRad(delta_rad%nelem,delta_rad%nscan,delta_rad%nMC,&
            delta_rad_tmp)
       delta_rad_tmp%ch1 = delta_rad%ch1
       delta_rad_tmp%ch2 = delta_rad%ch2
       delta_rad_tmp%ch3a = delta_rad%ch3a
       delta_rad_tmp%ch3 = delta_rad%ch3
       delta_rad_tmp%ch4 = delta_rad%ch4
       delta_rad_tmp%ch5 = delta_rad%ch5
    ENDIF

  END SUBROUTINE Copy_Delta_Rad

  SUBROUTINE Resize_Orbit(AVHRR,AVHRRout,startpos,endpos,all,&
       montecarlo,delta_rad)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(AVHRR_Data), INTENT(OUT), TARGET :: AVHRRout
    INTEGER, INTENT(IN) :: startpos
    INTEGER, INTENT(IN) :: endpos
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL, INTENT(IN), OPTIONAL :: montecarlo
    TYPE(mc_delta_str), INTENT(INOUT), OPTIONAL :: delta_rad

    ! Local variables
    INTEGER :: I,J,K
    INTEGER :: nsize
    INTEGER :: STAT
    LOGICAL :: alldata
    TYPE(AVHRR_Data), POINTER :: pAVHRRout
    INTEGER :: mem_scale
    TYPE(mc_delta_str) :: delta_rad_tmp
    LOGICAL :: monte_carlo

    IF( PRESENT(montecarlo) )THEN
       monte_carlo = montecarlo
       IF( .not. PRESENT(delta_rad) )THEN
          CALL Gbcs_Critical(.TRUE.,'Must include delta_rad for MC',&
               'Resize_Orbit','combine_orbits.f90')
       ENDIF
    ELSE
       monte_carlo = .FALSE.
    ENDIF

    IF( monte_carlo )THEN
       CALL Copy_Delta_Rad(delta_rad,delta_rad_tmp)
    ENDIF

    pAVHRRout => AVHRRout

    IF( PRESENT(all) )THEN
       alldata = all
    ELSE
       alldata = .TRUE.
    ENDIF

    !
    ! Allocate right size of output
    !
    AVHRRout%walton_there = .FALSE.

    CALL Allocate_OldCal(AVHRR,pAVHRRout,startpos,endpos)
    !
    ! Set all data to bad
    !
    DO I=1,AVHRRout%arraySize
       CALL INIT_OutData_Scanline(AVHRRout,I)
    END DO

    IF( alldata .and. AVHRR%newCalibration_There )THEN
       !
       ! Allocate newcal arrays and fill
       !
       CALL Allocate_NewCal( AVHRR, AVHRRout, alldata )

    ENDIF

    CALL Copy_Top_Data(AVHRR,AVHRRout,startpos,endpos)

    K = 0
    DO I=startpos,endpos
       K=K+1
       CALL Copy_All_Scan(AVHRR,AVHRRout,I,K,alldata)
    END DO

    IF( monte_carlo )THEN
       CALL Deallocate_DeltaRad(delta_rad)
       CALL Copy_Delta_Rad(delta_rad_tmp,delta_rad,startpos=startpos,&
            endpos=endpos)
       CALL Deallocate_DeltaRad(delta_rad_tmp)
    ENDIF

  END SUBROUTINE Resize_Orbit

  !
  ! TEST TEST TEST TEST TEST TEST TEST TEST
  !
  SUBROUTINE Check_latlon(AVHRR,lat,lon,usenew)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    REAL, INTENT(IN) :: lat
    REAL, INTENT(IN) :: lon
    LOGICAL, OPTIONAL :: usenew

    INTEGER :: I,J
    REAL :: dist
    REAL :: min_dist
    INTEGER :: IPOS, JPOS
    LOGICAL :: use_new

    IF( PRESENT(usenew) )THEN
       use_new = usenew
    ELSE
       use_new = .FALSE.
    ENDIF

    IPOS = -1
    JPOS = -1
    min_dist = 1e30

    DO I=1,AVHRR%arraySize
       DO J=1,AVHRR%nelem
          IF( AVHRR%lat(J,I) .gt. -100 .and. AVHRR%lon(J,I) .gt. -200. )THEN
             dist = acos(sin(AVHRR%lat(J,I)*Deg2Rad)*&
                  sin(lat*Deg2Rad) + &
                  cos(AVHRR%lat(J,I)*Deg2Rad)*&
                  cos(lat*Deg2Rad)*&
                  cos((AVHRR%lon(J,I)-lon)*Deg2Rad))*6371.
             IF( dist .lt. min_dist )THEN
                min_dist = dist
                IPOS = I
                JPOS = J
             ENDIF
          ENDIF
       END DO
    END DO

    IF( -1 .ne. IPOS .and. -1 .ne. JPOS )THEN
       WRITE(*,'(''======================================================'')')
       WRITE(*,*)'min dist (km) = ',min_dist,JPOS,IPOS
       WRITE(*,*)'Lat/Lon = ',AVHRR%lat(JPOS,IPOS),AVHRR%lon(JPOS,IPOS)
       IF( use_new )THEN
          WRITE(*,*)'Rad = ',AVHRR%array3b(JPOS,IPOS),&
               AVHRR%array4(JPOS,IPOS),&
               AVHRR%array5(JPOS,IPOS)
          WRITE(*,*)'BT = ',convertBT(AVHRR%array3b(JPOS,IPOS),&
               DBLE(AVHRR%nuc(1)),DBLE(AVHRR%aval(1)),DBLE(AVHRR%bval(1))),&
               convertBT(AVHRR%array4(JPOS,IPOS),&
               DBLE(AVHRR%nuc(2)),DBLE(AVHRR%aval(2)),DBLE(AVHRR%bval(2))),&
               convertBT(AVHRR%array5(JPOS,IPOS),&
               DBLE(AVHRR%nuc(3)),DBLE(AVHRR%aval(3)),DBLE(AVHRR%bval(3)))
       ELSE
          WRITE(*,*)'BT = ',AVHRR%array3b(JPOS,IPOS),&
               AVHRR%array4(JPOS,IPOS),&
               AVHRR%array5(JPOS,IPOS)
       ENDIF
       WRITE(*,'(''======================================================'')')
    ENDIF

  END SUBROUTINE Check_latlon
  !
  ! TEST TEST TEST TEST TEST TEST TEST TEST
  !
  SUBROUTINE read_file(infile,AVHRR,uuid_in,pygac_stem,fileno,out_instr_coefs)

    USE IFPORT

    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE(AVHRR_Data), TARGET, INTENT(OUT) :: AVHRR
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    CHARACTER(LEN=*), INTENT(IN) :: pygac_stem
    INTEGER, INTENT(IN) :: fileno
    TYPE(AVHRR_Instrument_Coefs), INTENT(out), OPTIONAL :: out_instr_coefs

    ! Local variables
    INTEGER :: pid
    INTEGER :: ostat
    INTEGER :: POS

    CHARACTER(LEN=256) :: inDirectory
    CHARACTER(LEN=256) :: infilename
    CHARACTER(LEN=256) :: temp_file
    CHARACTER(LEN=256) :: command_str

    TYPE(Imagery) :: IMG
    TYPE(AVHRR_Data), POINTER :: pAVHRR
    REAL :: coefs1(7)
    REAL :: coefs2(7)
    REAL :: coefs3(7)
    LOGICAL :: remove_file

    ! Check to see if we need to uncompress the data
    remove_file=.FALSE.
    IF( 0 .ne. INDEX(infile,'.gz') )THEN
       WRITE(command_str,'(''cp -f '',a,'' temp_file.'',a,''.gz'')')&
            TRIM(infile),TRIM(uuid_in)
       ostat = SYSTEM(command_str)
       CALL Gbcs_Critical(ostat.ne.0,'Cannot copy file to tempfile','Convert',&
            'extract_l1b_data.f90')
       WRITE(command_str,'(''gunzip -f temp_file.'',a,''.gz'')')&
            TRIM(uuid_in)
       ostat = SYSTEM(command_str)
       CALL Gbcs_Critical(ostat.ne.0,'Cannot gunzip tempfile','Convert',&
            'extract_l1b_data.f90')
       inDirectory = './'
       WRITE(inFilename,'(''temp_file.'',a)')TRIM(uuid_in)
       remove_file=.TRUE.
    ELSE
       POS=INDEX(infile,'/',.TRUE.)
       IF( 0 .ne. POS )THEN
          inDirectory = infile(1:POS)
       ELSE
          inDirectory = './'
       ENDIF
       inFilename = infile(POS+1:LEN(infile))
    ENDIF

    IMG%DataDir = TRIM(inDirectory)
    IMG%DataFile = TRIM(inFilename)

    IMG%GbcsDataPath = '/group_workspaces/jasmin2/nceo_uor/users/jmittaz/AVHRR/&
         &code/git/gbcs_new_calibration/dat_cci/'
    CALL Load_Imagery(IMG,outputData=AVHRR,use_new_calibration=.FALSE.,&
         out_instr_coefs=out_instr_coefs,use_walton=.FALSE.)
    AVHRR%scnline_l1b = fileno

!    CALL check_latlon(AVHRR,47.609,-29.85500)

    pAVHRR => AVHRR
    CALL Overlay_PyGAC_Data(pygac_stem,pAVHRR)

!    print *,'================= AFTER PyGAC ======================='
!    CALL check_latlon(AVHRR,47.609,-29.85500)

    IF( remove_file )THEN
       WRITE(command_str,'(''rm -f temp_file.'',a)')&
            TRIM(uuid_in)
       ostat = SYSTEM(command_str)
       CALL Gbcs_Critical(ostat.ne.0,'Cannot remove tempfile','Convert',&
            'extract_l1b_data.f90')
    ENDIF

  END SUBROUTINE read_file

  !
  ! Read pyGAC data and overwrite AVHRR information
  !
  SUBROUTINE Overlay_PyGAC_Data(pygac_stem,AVHRR)

    USE GbcsDateTime

    CHARACTER(LEN=*), INTENT(IN) :: pygac_stem
    TYPE(AVHRR_Data), POINTER, INTENT(INOUT) :: AVHRR

    ! Local data
    INTEGER :: I,J,K
    INTEGER :: pos
    INTEGER :: ncid
    INTEGER :: ostat
    INTEGER :: lat_id
    INTEGER :: grp_lat_id
    INTEGER :: nx
    INTEGER :: ny
    INTEGER :: qnx
    INTEGER :: qny
    INTEGER :: type
    INTEGER :: nmissing
    INTEGER :: scanline_step
    LOGICAL :: found
    REAL :: gain
    REAL :: offset
    REAL, ALLOCATABLE :: lat(:,:)
    REAL, ALLOCATABLE :: lon(:,:)
    REAL, ALLOCATABLE :: satza(:,:)
    REAL, ALLOCATABLE :: solza(:,:)
    REAL, ALLOCATABLE :: relaz(:,:)
    REAL, ALLOCATABLE :: ch1(:,:)
    REAL, ALLOCATABLE :: ch2(:,:)
    REAL, ALLOCATABLE :: ch3a(:,:)
    REAL, ALLOCATABLE :: ch3b(:,:)
    REAL, ALLOCATABLE :: ch4(:,:)
    REAL, ALLOCATABLE :: ch5(:,:)
    REAL, ALLOCATABLE :: temp_output(:,:)
    INTEGER(GbcsInt2), ALLOCATABLE :: qualflags(:,:)
    CHARACTER(LEN=256) :: filename
    CHARACTER(LEN=15) :: timestr(4)
    TYPE(DateTime) :: start_time
    TYPE(DateTime) :: stop_time
    REAL(GbcsDble) :: start_time_sec
    REAL(GbcsDble) :: time_delta
    REAL(GbcsDble) :: start_time_gbcs
    REAL(GbcsDble) :: new_time
    TYPE(DateTime) :: new_time_gbcs
    REAL :: hours
    LOGICAL :: ch3a_there
    LOGICAL :: ch5_there

    !
    ! Input name will be of form ECC_GAC_*_noaa*_start/end.h5
    ! where first * can be avhrr, qualflags,sunsatangles
    !
    ! lat/lon in avhrr version
    ! scanline numbering in qualflags
    ! angles in sunsatangles
    !
    ! Find _GAC_sunsatangles_ for lat/lon/angles
    !
    pos = INDEX(pygac_stem,'_GAC_avhrr_')
    IF( 0 .eq. pos )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot find _GAC_avhrr_ in pyGAC filename',&
            'Overlay_PyGAC_Data','combine_orbits.f90')
    END IF

    !
    ! Get Refl/BTs so we can check for fill values as there is extra 
    ! pyGAC filtering for scan motor errors where data is set to fill
    !
    filename = pygac_stem(1:pos-1)//'_GAC_avhrr_'//&
         &pygac_stem(pos+11:LEN_TRIM(pygac_stem))

    ostat = NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid)
    CALL check(ostat)

    CALL Get_Data_Grps(ncid,'None','Channel 1 reflectance',ch1,nx,ny,&
         tag_check=.FALSE.)
    CALL Get_Data_Grps(ncid,'None','Channel 2 reflectance',ch2,nx,ny,&
         tag_check=.FALSE.)
    CALL Get_Data_Grps(ncid,'None','Channel 3a reflectance',ch3a,nx,ny,&
         tag_check=.FALSE.)
    ch3a_there = .FALSE.
    IF( ALLOCATED(ch3a) )THEN
       ch3a_there = .TRUE.
    ENDIF
    CALL Get_Data_Grps(ncid,'None','Channel 3b brightness temperature',&
         ch3b,nx,ny,tag_check=.FALSE.)
    CALL Get_Data_Grps(ncid,'None','Channel 4 brightness temperature',&
         ch4,nx,ny,tag_check=.FALSE.)
    CALL Get_Data_Grps(ncid,'None','Channel 5 brightness temperature',&
         ch5,nx,ny,tag_check=.FALSE.)
    ch5_there = .FALSE.
    IF( ALLOCATED(ch5) )THEN
       ch5_there = .TRUE.
    ENDIF

    ostat = NF90_CLOSE(ncid)
    CALL check(ostat)

    !
    ! Filename to get lat/lon
    !
    filename = pygac_stem(1:pos-1)//'_GAC_sunsatangles_'//&
         &pygac_stem(pos+11:LEN_TRIM(pygac_stem))

    ostat = NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid)
    CALL check(ostat)

    !
    ! Get variables
    !
    CALL Get_Data_Grps(ncid,'where','Latitude',lat,nx,ny,gettime=timestr)    

    !
    ! Get rest of variables in this file
    !
    CALL Get_Data_Grps(ncid,'where','Longitude',lon,nx,ny)
    CALL Get_Data_Grps(ncid,'None','Satellite zenith angle',satza,nx,ny)
    CALL Get_Data_Grps(ncid,'None','Solar zenith angle',solza,nx,ny)
    CALL Get_Data_Grps(ncid,'None','Relative satellite-sun azimuth angle',&
         relaz,nx,ny)

    !
    ! Close file
    !
    ostat = NF90_CLOSE(ncid)
    CALL check(ostat)

    !
    ! Read in quality flags
    !
    !
    ! Filename to get lat/lon
    !
    filename = pygac_stem(1:pos-1)//'_GAC_qualflags_'//&
         &pygac_stem(pos+11:LEN_TRIM(pygac_stem))

    ostat = NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid)
    CALL check(ostat)

    !
    ! Get variables
    !
    CALL Get_Data_Grps(ncid,'None','Scanline quality flags',&
         temp_output,qnx,qny,short_output=qualflags,short_type=.TRUE.)

    !
    ! Close file
    !
    ostat = NF90_CLOSE(ncid)
    CALL check(ostat)

    !
    ! Setup time in pygac frame - uses first and last scanline number from
    ! qualflags
    !
    CALL Setup_PyGAC_Time(AVHRR%time_yearstart,timestr,&
         qualflags(1,1),qualflags(1,ny),start_time,&
         stop_time,start_time_sec,time_delta,start_time_gbcs)
    !
    ! Copy over lat/lon/angles
    ! Use first entry of qualflags as scanline number
    !
    AVHRR%scan_line_delta_time = time_delta
    nmissing = 0
    DO I=1,AVHRR%arraySize
       found = .FALSE.
       InnerLoop: DO J=1,qny
          IF( qualflags(1,J) .eq. AVHRR%scanLineNumber(I) )THEN
             !
             ! Found matching line
             !
             found = .TRUE.
             !
             ! Copy over data including updated time
             !
             AVHRR%Lon(:,I) = lon(:,J)
             AVHRR%Lat(:,I) = lat(:,J)
             AVHRR%satZA(:,I) = satza(:,J)
             AVHRR%solZA(:,I) = solza(:,J)
             AVHRR%relAz(:,I) = relaz(:,J)
             scanline_step = (qualflags(1,J)-qualflags(1,1))
             AVHRR%Time(I) = start_time_sec + scanline_step*time_delta
             new_time = start_time_gbcs + scanline_step*time_delta/86400.d0
             new_time_gbcs = JD_to_Date(new_time)
             AVHRR%year(I) = new_time_gbcs%year
             AVHRR%month(I) = new_time_gbcs%month
             AVHRR%day(I) = new_time_gbcs%day
             hours = new_time_gbcs%hour + new_time_gbcs%minute/60. + &
                  new_time_gbcs%seconds/3600. + &
                  new_time_gbcs%sec1000/3600000.
             AVHRR%hours(I) = hours
             AVHRR%UTC_msecs(I) = INT((hours/24.)*86400000.)
             !
             ! Check to see if there is bad data in pyGAC values
             ! They do extra flagging
             !
             DO K=1,AVHRR%nelem
                IF( ch1(K,J) .eq. NAN_R )THEN
                   AVHRR%counts1(K,I) = NAN_R
                ENDIF
                IF( ch2(K,J) .eq. NAN_R )THEN
                   AVHRR%counts2(K,I) = NAN_R
                ENDIF
                IF( ch3a_there )THEN
                   IF( AVHRR%array3A(K,I) .gt. 0. )THEN
                      IF( ch3a(K,J) .eq. NAN_R )THEN
                         AVHRR%counts3(K,I) = NAN_R
                      ENDIF
                   ENDIF
                ENDIF
                IF( AVHRR%array3B(K,I) .gt. 0. )THEN
                   IF( ch3b(K,J) .eq. NAN_R )THEN
                      AVHRR%counts3(K,I) = NAN_R
                   ENDIF
                ENDIF
                IF( ch4(K,J) .eq. NAN_R )THEN
                   AVHRR%counts4(K,I) = NAN_R
                ENDIF
                IF( ch5_there )THEN
                   IF( ch4(K,J) .eq. NAN_R )THEN
                      AVHRR%counts4(K,I) = NAN_R
                   ENDIF
                ENDIF
             END DO
          ENDIF
       END DO InnerLoop
       IF( .not. found )THEN
          nmissing = nmissing+1
          !
          ! Set time to bad to be used as a check when gap filling
          !
          AVHRR%Time(I) = -1e30
          AVHRR%year(I) = -1
          AVHRR%month(I) = -1
          AVHRR%day(I) = -1
          AVHRR%hours(I) = -1
          AVHRR%UTC_msecs(I) = -1
          AVHRR%Lon(:,I) = NAN_R
          AVHRR%Lat(:,I) = NAN_R
       ENDIF
    END DO
    WRITE(*,*)'PyGAC comparison Number of missing lines : ',nmissing

    DEALLOCATE(lon,lat,solza,satza,relaz,ch1,ch2,ch3b,ch4)
    IF( ch3a_there )THEN
       DEALLOCATE(ch3a)
    ENDIF
    IF( ch5_there )THEN
       DEALLOCATE(ch5)
    ENDIF
        
  END SUBROUTINE Overlay_PyGAC_Data

  SUBROUTINE Get_Data_Grps(ncid,top_name,tag_name,output,nx,ny,short_output,&
       short_type,gettime,tag_check)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: top_name
    CHARACTER(LEN=*), INTENT(IN) :: tag_name
    REAL, ALLOCATABLE, INTENT(OUT) :: output(:,:)
    INTEGER, INTENT(OUT) :: nx
    INTEGER, INTENT(OUT) :: ny
    INTEGER(GbcsInt2), ALLOCATABLE, INTENT(OUT), OPTIONAL :: short_output(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: short_type
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: gettime(4)
    LOGICAL, INTENT(IN), OPTIONAL :: tag_check

    ! Local variables
    INTEGER :: I,J
    INTEGER :: ostat
    INTEGER :: grp_id
    INTEGER :: var_id
    INTEGER :: type
    INTEGER(GbcsInt2), ALLOCATABLE :: var_short(:,:)
    INTEGER(GbcsInt4), ALLOCATABLE :: var_int(:,:)
    REAL :: gain
    REAL :: offset
    INTEGER(GbcsInt2) :: missingdata_short
    INTEGER(GbcsInt2) :: nodata_short
    INTEGER(GbcsInt4) :: missingdata_int
    INTEGER(GbcsInt4) :: nodata_int
    REAL :: missingdata_real
    REAL :: nodata_real
    LOGICAL :: shorttype
    LOGICAL :: ok
    LOGICAL :: check_ok

    IF( PRESENT(short_type) )THEN
       shorttype = short_type
    ELSE
       shorttype = .FALSE.
    ENDIF

    IF( PRESENT(tag_check) )THEN
       check_ok = tag_check
    ELSE
       check_ok = .TRUE.
    ENDIF
    CALL Get_Array_Grps(ncid,top_name,tag_name,grp_id,var_id,nx,ny,&
         type,gain,offset,missingdata_short,nodata_short,&
         missingdata_int,nodata_int,missingdata_real,nodata_real,&
         ok,gettime=gettime)
    IF( check_ok .and. .not. ok )THEN
       WRITE(*,'(''tag_name = '',a)')TRIM(tag_name)
       CALL Gbcs_Critical(.TRUE.,'Cannot find tag_name','Get_Data_Grps',&
            'combin_orbits.f90')
    ELSE IF( .not. ok )THEN
       RETURN
    ENDIF
    IF( .not. shorttype )THEN
       ALLOCATE(output(nx,ny),STAT=ostat)
    ELSE
       IF( .not. PRESENT(short_output) )THEN
          CALL Gbcs_Critical(.TRUE.,'Optional short_output not present',&
               'Get_Data_Grps','combine_orbits.f90')
       ENDIF
       ALLOCATE(short_output(nx,ny),STAT=ostat)
    ENDIF
    IF( 0 .ne. ostat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate main array',&
            'Get_Data_Grps','combine_orbits.f90')
    ENDIF
    !
    ! Depending on type of variable, read in and scale to float
    !
    IF( type .eq. NF90_SHORT )THEN
       ALLOCATE(var_short(nx,ny),STAT=ostat)
       IF( 0 .ne. ostat )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate var_short array',&
               'Get_Data_Grps','combine_orbits.f90')
       ENDIF
       ostat = NF90_GET_VAR(grp_id,var_id,var_short)
       CALL check(ostat)
       IF( shorttype )THEN
          DO I=1,ny
             DO J=1,nx
                IF( var_short(J,I) .eq. missingdata_short .or. &
                     var_short(J,I) .eq. nodata_short )THEN
                   short_output(J,I) = NAN_I
                ELSE
                   short_output(J,I) = var_short(J,I)
                ENDIF
             END DO
          END DO
       ELSE
          DO I=1,ny
             DO J=1,nx
                IF( var_short(J,I) .eq. missingdata_short .or. &
                     var_short(J,I) .eq. nodata_short )THEN
                   output(J,I) = NAN_R
                ELSE
                   output(J,I) = offset + gain*var_short(J,I)
                ENDIF
             END DO
          END DO
       ENDIF
       DEALLOCATE(var_short)
    ELSE IF( type .eq. NF90_INT )THEN
       ALLOCATE(var_int(nx,ny),STAT=ostat)
       IF( 0 .ne. ostat )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate var_short array',&
               'Get_Data_Grps','combine_orbits.f90')
       ENDIF
       ostat = NF90_GET_VAR(grp_id,var_id,var_int)
       CALL check(ostat)
       DO I=1,ny
          DO J=1,nx
             IF( var_int(J,I) .eq. missingdata_int .or. &
                  var_int(J,I) .eq. nodata_int )THEN
                output(J,I) = NAN_R
             ELSE
                output(J,I) = offset + gain*var_int(J,I)
             ENDIF
          END DO
       END DO       
       DEALLOCATE(var_int)
    ELSE IF( type .eq. NF90_FLOAT )THEN
       ostat = NF90_GET_VAR(grp_id,var_id,output)
       CALL check(ostat)
       DO I=1,ny
          DO J=1,nx
             IF( output(J,I) .eq. missingdata_real .or. &
                  output(J,I) .eq. nodata_real )THEN
                output(J,I) = NAN_R
             ENDIF
          END DO
       END DO       
    ENDIF

  END SUBROUTINE Get_Data_Grps

  !
  ! Note stored in groups so need group based reader
  !
  SUBROUTINE Get_Array_Grps(ncid,top_name,tag_name,grp_id,data_id,nx,ny,type,&
       gain,offset,missingdata_short,nodata_short,&
       missingdata_int,nodata_int,missingdata_real,nodata_real,ok,gettime)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: top_name
    CHARACTER(LEN=*), INTENT(IN) :: tag_name
    INTEGER, INTENT(OUT) :: grp_id
    INTEGER, INTENT(OUT) :: data_id
    INTEGER, INTENT(OUT) :: nx
    INTEGER, INTENT(OUT) :: ny
    INTEGER, INTENT(OUT) :: type
    REAL, INTENT(OUT) :: gain
    REAL, INTENT(OUT) :: offset
    INTEGER(GbcsInt2), INTENT(OUT) :: missingdata_short
    INTEGER(GbcsInt2), INTENT(OUT) :: nodata_short
    INTEGER(GbcsInt4), INTENT(OUT) :: missingdata_int
    INTEGER(GbcsInt4), INTENT(OUT) :: nodata_int
    REAL, INTENT(OUT) :: missingdata_real
    REAL, INTENT(OUT) :: nodata_real
    LOGICAL, INTENT(OUT) :: ok
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: gettime(4)

    ! Local variables
    INTEGER :: I,J
    INTEGER :: ostat
    INTEGER, PARAMETER :: max_ngrps=100
    INTEGER :: top_group_id
    INTEGER :: low_id
    INTEGER :: ntop
    INTEGER :: ntop_low
    INTEGER :: ngrps
    INTEGER :: n_attr
    INTEGER :: top_ncids(max_ngrps)
    INTEGER :: low_ncids(max_ngrps)
    INTEGER :: attr_ncid(max_ngrps)
    INTEGER :: include_parents
    CHARACTER(LEN=NF90_MAX_NAME) :: name
    CHARACTER(LEN=NF90_MAX_NAME) :: values

    ok = .TRUE.
    !
    ! Get list of groups
    !
    ostat = NF90_INQ_GRPS(ncid,ntop,top_ncids)
    CALL check(ostat)

    IF( top_name .ne. 'None' )THEN
       !
       ! Loop round groups to find top input tag
       !
       top_group_id=-1
       DO I=1,ntop
          ostat = NF90_INQ_GRPNAME(top_ncids(I),name) 
          CALL check(ostat)
          IF( TRIM(name) .eq. TRIM(top_name) )THEN
             top_group_id = top_ncids(I)
             EXIT
          ENDIF
       END DO
       IF( -1 .eq. top_group_id )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot top_group_id',&
               'Get_Array_Grps','combine_orbits.f90')
       ENDIF
    ELSE
       top_group_id = -1
    ENDIF

    !
    ! Get lower level ids - below where if present
    !
    IF( -1 .ne. top_group_id )THEN
       ostat = NF90_INQ_GRPS(top_group_id,ntop_low,low_ncids)
       CALL check(ostat)
    ELSE
       ntop_low = ntop
       low_ncids = top_ncids
    ENDIF

    !
    ! Loop round groups to find lower level input tag
    !
    data_id=-1
    nx = -1
    ny = -1
    low_loop: DO I=1,ntop_low
       ! 
       ! Get name from what group - dataset name
       !
       ostat = NF90_INQ_GRPS(low_ncids(I),n_attr,attr_ncid)
       CALL check(ostat)
       IF( 0 .eq. n_attr )THEN
          CALL Get_Grp_Attr(low_ncids(I),NF90_GLOBAL,tag_name,&
               gain,offset,type,grp_id,&
               data_id,missingdata_short,nodata_short,&
               missingdata_int,nodata_int,&
               missingdata_real,nodata_real,nx,ny,gettime=gettime)
          IF( -1 .ne. data_id )THEN
             EXIT low_loop
          ENDIF
       ELSE
          DO J=1,n_attr
             ostat = NF90_INQ_GRPNAME(attr_ncid(J),name) 
             CALL check(ostat)
             IF( name .eq. 'what' )THEN
                CALL Get_Grp_Attr(low_ncids(I),attr_ncid(J),tag_name,&
                     gain,offset,type,grp_id,&
                     data_id,missingdata_short,nodata_short,&
                     missingdata_int,nodata_int,&
                     missingdata_real,nodata_real,nx,ny,gettime=gettime)
                IF( -1 .ne. data_id )THEN
                   EXIT low_loop
                ENDIF
             ENDIF
          END DO
       ENDIF
    END DO low_loop
    IF( -1 .eq. data_id )THEN
       ok = .FALSE.
!       CALL Gbcs_Critical(.TRUE.,'Cannot get data_id',&
!            'Get_Array_Grps','combine_orbits.f90')
    ENDIF    

  END SUBROUTINE Get_Array_Grps
  
  SUBROUTINE Get_Grp_Attr(top_ncid,ncid,tag_name,gain,offset,type,grp_id,&
       data_id,missingdata_short,nodata_short,missingdata_int,nodata_int,&
       missingdata_real,nodata_real,nx,ny,gettime)

    INTEGER, INTENT(IN) :: top_ncid
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: tag_name
    REAL, INTENT(OUT) :: gain
    REAL, INTENT(OUT) :: offset
    INTEGER, INTENT(OUT) :: type
    INTEGER, INTENT(OUT) :: grp_id
    INTEGER, INTENT(OUT) :: data_id
    INTEGER(GbcsInt2), INTENT(OUT) :: missingdata_short
    INTEGER(GbcsInt2), INTENT(OUT) :: nodata_short
    INTEGER(GbcsInt4), INTENT(OUT) :: missingdata_int
    INTEGER(GbcsInt4), INTENT(OUT) :: nodata_int
    REAL, INTENT(OUT) :: missingdata_real
    REAL, INTENT(OUT) :: nodata_real
    INTEGER, INTENT(OUT) :: nx
    INTEGER, INTENT(OUT) :: ny
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: gettime(4)

    ! Local variables
    INTEGER, PARAMETER :: max_ngrps=100
    INTEGER :: ndims
    INTEGER :: nc_id
    INTEGER :: ostat
    INTEGER :: dim_ids(max_ngrps)
    INTEGER :: itemp
    CHARACTER(LEN=NF90_MAX_NAME) :: values
    CHARACTER(LEN=NF90_MAX_NAME) :: dim_name

    missingdata_short = -32001
    nodata_short = -32001
    missingdata_int = -999999
    nodata_int = -999999
    missingdata_real = NAN_R
    nodata_real = NAN_R

    grp_id = top_ncid
    data_id = -1

    IF( ncid .ne. NF90_GLOBAL )THEN
       ostat = NF90_GET_ATT(ncid, NF90_GLOBAL, 'dataset_name', values)
       nc_id = ncid
    ELSE
       ostat = NF90_GET_ATT(top_ncid, NF90_GLOBAL, 'dataset_name', values)
       nc_id = top_ncid
    ENDIF
    !
    ! If not found move one
    !
    IF( NF90_NOERR .ne. ostat )THEN
       RETURN
    ENDIF
    IF( values .eq. tag_name )THEN
       !
       ! Get gain and offset attributes
       !
       ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
            'gain', gain)
       CALL check(ostat)
       ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
            'offset', offset)
       CALL check(ostat)
       !
       ! Get pointer to data array and group_id
       !
       ostat = NF90_INQ_VARID(top_ncid,'data',data_id)
       CALL check(ostat)

       !
       ! Get data type and dimensions for data
       !
       ostat = NF90_INQUIRE_VARIABLE(top_ncid,data_id,xtype=type,&
            ndims=ndims,dimids=dim_ids)
       CALL check(ostat)

       IF( ndims .eq. 2 )THEN
          !
          ! Get dimension sizes
          !
          ! Note x/y switch due to file definition
          !
          ostat = NF90_INQUIRE_DIMENSION(top_ncid,dim_ids(1),dim_name,ny)
          CALL check(ostat)
          ostat = NF90_INQUIRE_DIMENSION(top_ncid,dim_ids(2),dim_name,nx)
          CALL check(ostat)
          !
          ! Sometimes seems we have a switch
          !
          IF( nx .ne. 409 )THEN
             itemp = nx
             nx = ny
             ny = itemp
          ENDIF
          !
          ! Get relevant missing data numbers
          !
          IF( type .eq. NF90_SHORT )THEN
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'missingdata', missingdata_short)
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'nodata', nodata_short)
             CALL check(ostat)
          ELSE IF( type .eq. NF90_INT )THEN
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'missingdata', missingdata_int)
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'nodata', nodata_int)
             CALL check(ostat)
          ELSE IF( type .eq. NF90_FLOAT )THEN
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'missingdata', missingdata_real)
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'nodata', nodata_real)
             CALL check(ostat)
          ENDIF
          
          !
          ! If requested get the time strings
          !
          IF( PRESENT(gettime) )THEN
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'startdate', gettime(1))
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'starttime', gettime(2))
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'enddate', gettime(3))
             CALL check(ostat)
             ostat = NF90_GET_ATT(nc_id, NF90_GLOBAL, &
                  'endtime', gettime(4))
             CALL check(ostat)
          ENDIF
       ELSE
          CALL Gbcs_Critical(.TRUE.,'Only 2d data arrays supported',&
               'Get_Grp_Attr','combine_orbits.f90')
       ENDIF
    END IF

  END SUBROUTINE Get_Grp_Attr

  SUBROUTINE fill_missing_lines(AVHRR,AVHRRout,montecarlo,delta_rad)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(AVHRR_Data), INTENT(OUT), TARGET :: AVHRRout
    LOGICAL, INTENT(IN), OPTIONAL :: montecarlo
    TYPE(mc_delta_str), INTENT(INOUT), OPTIONAL :: delta_rad

    ! Local variables
    INTEGER :: I
    INTEGER :: K
    INTEGER :: KK
    INTEGER :: II
    INTEGER :: start_line
    INTEGER :: stop_line
    INTEGER :: ntotal
    INTEGER :: nmissing
    INTEGER :: nmissing_total
    TYPE(AVHRR_Data), POINTER :: pAVHRRout
    LOGICAL :: too_many
    LOGICAL :: monte_carlo
    TYPE(mc_delta_str) :: delta_rad_tmp
    
    IF( PRESENT(montecarlo) )THEN
       IF( .not. PRESENT(delta_rad) )THEN
          CALL Gbcs_Critical(.TRUE.,'Must include delta_rad',&
               'fill_missing_lines','combine_orbits.f90')
       ENDIF
       monte_carlo = montecarlo
    ELSE
       monte_carlo = .FALSE.
    ENDIF

    ! Find how many lines needed for output
    start_line = -1
    startLoop: DO I=1,AVHRR%arraySize
       IF( AVHRR%year(I) .gt. 0 )THEN
          start_line = I
          EXIT startLoop
       END IF
    END DO startLoop
    IF( start_line .lt. 1 )THEN
       CALL Gbcs_Critical(.TRUE.,'No good data found','fill_missing_lines',&
            'combine_orbits.f90')
    ENDIF

    stop_line = -1
    startLoop2: DO I=AVHRR%arraySize,1,-1
       IF( AVHRR%year(I) .gt. 0 )THEN
          stop_line = I
          EXIT startLoop2
       END IF
    END DO startLoop2
    IF( stop_line .lt. 1 )THEN
       CALL Gbcs_Critical(.TRUE.,'No good data found','fill_missing_lines',&
            'combine_orbits.f90')
    ENDIF

    ! How many lines needed in output
    too_many = .FALSE.
    ntotal = NINT(((AVHRR%time(stop_line)-AVHRR%time(start_line))/&
             AVHRR%scan_line_delta_time))+1
    IF( ntotal .lt. stop_line-start_line+1 )THEN
       too_many = .TRUE.
       pAVHRRout => AVHRRout
       CALL Allocate_OldCal( AVHRR, pAVHRRout,start_line,stop_line)
       IF( monte_carlo )THEN
          CALL Allocate_DeltaRad( delta_rad%nelem,stop_line-start_line+1,&
               delta_rad%nMC,delta_rad_tmp)
       ENDIF
    ELSE
       ! Allocate output structure
       pAVHRRout => AVHRRout
       CALL Allocate_OldCal( AVHRR, pAVHRRout,start_line,start_line+ntotal)
       IF( monte_carlo )THEN
          CALL Allocate_DeltaRad( delta_rad%nelem,ntotal+1,&
               delta_rad%nMC,delta_rad_tmp)
       ENDIF
    ENDIF
    !
    ! Set all data to bad
    !
    DO I=1,AVHRRout%arraySize
       CALL INIT_OutData_Scanline(AVHRRout,I)
    END DO
    AVHRRout%missingLines = .TRUE.
    !
    ! New cal if needed
    !
    IF( AVHRR%newCalibration_There )THEN
       CALL Allocate_NewCal(AVHRR,AVHRRout,.TRUE.)
    ENDIF

    CALL Copy_Top_Data(AVHRR,AVHRRout,1,ntotal)

    !
    ! If ntotal < stop_line-start_line+1 then there is a problem with the 
    ! times and we have too much data. So we copy it all
    !
    nmissing_total = 0
    IF( too_many )THEN
       CALL Gbcs_Warning(.TRUE.,&
            'Problems with timing data ntotal < stop_line-start_line+1 - copying all data',&
            'fill_missing_lines',&
            'combine_orbits.f90')
       I=start_line
       K=1
       DO WHILE(I .le. stop_line)
          CALL Copy_All_Scan(AVHRR,AVHRRout,I,K,.TRUE.)
          IF( monte_carlo )THEN
             CALL Copy_DeltaRad_Ind(delta_rad_tmp,delta_rad,I,K=K)
          ENDIF
          AVHRRout%missingLines(K) = .FALSE.
          I=I+1
          K=K+1
       END DO
    ELSE
       !
       ! We have to make up missing lines
       !
       
       K=0
       I=start_line
       DO WHILE(I .le. stop_line)
          !
          ! If pyGAC copied data there, copy line and see if there is a 
          ! data gap
          !
          IF( AVHRR%year(I) .gt. 0 )THEN
             K=K+1
             IF( K .gt. AVHRRout%arraySize )THEN
                CALL Gbcs_Critical(.TRUE.,'K > arraySize',&
                     'fill_missing_lines','combine_orbits.f90')
             ENDIF
             CALL Copy_All_Scan(AVHRR,AVHRRout,I,K,.TRUE.)
             IF( monte_carlo )THEN
                CALL Copy_DeltaRad_Ind(delta_rad_tmp,delta_rad,I,K=K)
             ENDIF
             AVHRRout%missingLines(K) = .FALSE.
             !
             ! Check time to next scan line
             !
             II=I
             SearchLoop: DO WHILE(I .lt. stop_line )
                IF( AVHRR%year(I+1) .gt. 0 )THEN
                   nmissing = NINT((AVHRR%time(I+1)-AVHRR%time(II))/&
                        AVHRR%scan_line_delta_time)-1
                   IF( 0 .lt. nmissing )THEN
                      K=K+nmissing
                      nmissing_total = nmissing_total + nmissing
                   ENDIF
                   EXIT SearchLoop
                ELSE
                   I=I+1
                ENDIF
             END DO SearchLoop
          ElSE
             K=K+1
          ENDIF
          I=I+1
       END DO
    ENDIF

    WRITE(*,'(''  Number of missing scanlines added = '',i5)')nmissing_total

  END SUBROUTINE fill_missing_lines

  SUBROUTINE Setup_PyGAC_Time(time_yearstart,timestr,start_scnline,&
       stop_scnline,start_time,stop_time,start_time_sec,time_delta,&
       start_time_gbcs)

    INTEGER, INTENT(IN) :: time_yearstart
    CHARACTER(LEN=*), INTENT(IN) :: timestr(4)
    INTEGER(GbcsInt2), INTENT(IN) :: start_scnline
    INTEGER(GbcsInt2), INTENT(IN) :: stop_scnline
    TYPE(DateTime), INTENT(OUT) :: start_time
    TYPE(DateTime), INTENT(OUT) :: stop_time
    REAL(GbcsDble), INTENT(OUT) :: start_time_sec
    REAL(GbcsDble), INTENT(OUT) :: time_delta
    REAL(GbcsDble), INTENT(OUT) :: start_time_gbcs

    ! Local variables
    INTEGER :: IOS
    INTEGER :: dayno
    INTEGER :: microsec
    REAL :: hours
    REAL(GbcsDble) :: stop_time_sec
    CHARACTER(LEN=15) :: timestr_entry
    
    timestr_entry = timestr(1)
    READ(timestr_entry(1:4),*,IOSTAT=IOS)start_time%year
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start year',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(5:6),*,IOSTAT=IOS)start_time%month
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start month',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(7:8),*,IOSTAT=IOS)start_time%day
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start day',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    timestr_entry = timestr(2)
    READ(timestr_entry(1:2),*,IOSTAT=IOS)start_time%hour
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start hour',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(3:4),*,IOSTAT=IOS)start_time%minute
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start minute',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(5:6),*,IOSTAT=IOS)start_time%seconds
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start second',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(8:13),*,IOSTAT=IOS)microsec
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start microsecond',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    start_time%sec1000 = INT((microsec/1e6)*1e3)
    start_time%utc_offset = 0

    timestr_entry = timestr(3)
    READ(timestr_entry(1:4),*,IOSTAT=IOS)stop_time%year
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop year',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(5:6),*,IOSTAT=IOS)stop_time%month
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop month',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(7:8),*,IOSTAT=IOS)stop_time%day
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop day',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    timestr_entry = timestr(4)
    READ(timestr_entry(1:2),*,IOSTAT=IOS)stop_time%hour
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop hour',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(3:4),*,IOSTAT=IOS)stop_time%minute
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop minute',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(5:6),*,IOSTAT=IOS)stop_time%seconds
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse stop seconds',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    READ(timestr_entry(8:13),*,IOSTAT=IOS)microsec
    IF( 0 .ne. IOS )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot parse start microsecond',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF
    stop_time%sec1000 = INT((microsec/1e6)*1e3)
    stop_time%utc_offset = 0

    dayno = Day_Of_Year(start_time%year,start_time%month,start_time%day)
    hours = start_time%hour + start_time%minute/60. + &
         start_time%seconds/3600.+start_time%sec1000/3600000.
    start_time_sec = get_date_time(start_time%year,dayno,hours,time_yearstart)

    dayno = Day_Of_Year(stop_time%year,stop_time%month,stop_time%day)
    hours = stop_time%hour + stop_time%minute/60. + &
         stop_time%seconds/3600.+stop_time%sec1000/3600000.
    stop_time_sec = get_date_time(stop_time%year,dayno,hours,time_yearstart)

    time_delta = (stop_time_sec - start_time_sec)/&
         (1.*(stop_scnline-start_scnline-1))
    !
    ! Make sure time delta is 0.5 seconds
    !
    IF( ABS(INT(1000*time_delta+0.5)-500) .gt. 1 )THEN
       WRITE(*,*)'NY:',(stop_scnline-start_scnline)+1
       WRITE(*,*)'start_time:',start_time
       WRITE(*,*)' stop_time:',stop_time
       WRITE(*,*)'Time delta = ',time_delta
       CALL Gbcs_Critical(.TRUE.,'time delta != 0.5 secs',&
            'Setup_PyGAC_Time','combine_orbits.f90')
    ENDIF

    !
    ! Force time delta to be 0.5 seconds
    !
    time_delta = 0.5

    start_time_gbcs = Date_to_JD(start_time)

  END SUBROUTINE Setup_PyGAC_Time

  LOGICAL FUNCTION isNaN_Dble( inval )

    REAL(GbcsDBle), INTENT(IN) :: inval

    IF( inval .ne. inval )THEN
       isNaN_Dble = .TRUE.
    ELSE
       isNaN_Dble = .FALSE.
    ENDIF

  END FUNCTION isNaN_Dble

  LOGICAL FUNCTION isNaN_Real( inval )

    REAL(GbcsReal), INTENT(IN) :: inval

    IF( inval .ne. inval )THEN
       isNaN_Real = .TRUE.
    ELSE
       isNaN_Real = .FALSE.
    ENDIF

  END FUNCTION isNaN_Real

  SUBROUTINE check(ostat)
    
    INTEGER, INTENT (IN) :: ostat
    IF(ostat /= nf90_noerr) THEN
      WRITE(*,*)TRIM(nf90_strerror(ostat))
      stop "Stopped"
    END IF

  END SUBROUTINE check

END MODULE Combine_Orbits


