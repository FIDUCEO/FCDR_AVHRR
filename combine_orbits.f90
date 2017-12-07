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
! Module to combine three AVHRR orbits to create one 'good' orbit
! equator to equator
!
MODULE Combine_Orbits
  
  USE GbcsKinds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsErrorHandler
  USE GbcsDateTime
  USE NOAA_LoadAVHRRLevel1B
  USE fiduceo_uncertainties

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

    AVHRR%scanLineNumber(POS1) = &
         AVHRR_New%scanLineNumber(POS2)
    AVHRR%badTop(POS1) = &
         AVHRR_New%badTop(POS2)
    AVHRR%badTime(POS1) = &
         AVHRR_New%badTime(POS2)
    AVHRR%badNavigation(POS1) = &
         AVHRR_New%badNavigation(POS2)
    AVHRR%badCalibration(POS2) = &
         AVHRR_New%badCalibration(POS2)
    AVHRR%Lon(:,POS2) = &
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
          diff = ABS(time(I)-newtime(J+K))
          IF( diff .gt. 2 )THEN
             EXIT FindLoop
          ENDIF
          IF( diff .lt. min_diff )THEN
             outJ = J+K
             min_diff = diff
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

  SUBROUTINE Merge_AVHRR(AVHRR,AVHRR_New)

    TYPE(AVHRR_Data), INTENT(INOUT), TARGET :: AVHRR
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR_new
    
    ! Local variables
    INTEGER :: I
    INTEGER :: J
    INTEGER :: outJ
    INTEGER :: stat
    INTEGER :: start_pos
    INTEGER :: extra_lines
    INTEGER :: first_position
    INTEGER :: last_position
    INTEGER :: back_pos
    INTEGER :: back_pos2
    INTEGER :: offset
    
    REAL(GbcsDble) :: start_time 

    LOGICAL, ALLOCATABLE :: good_data(:)

    TYPE(AVHRR_Data), POINTER :: AVHRR_Ptr

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
       IF( good_data(I) )THEN
          last_position=I
          EXIT FindLast
       ENDIF
    END DO FindLast
    IF( -1 .eq. last_position )THEN
       CALL Gbcs_Critical(.TRUE.,'No good data found in first input structure',&
            'Merge_AVHRR','extract_l1b_data.f90')
    ENDIF
    IF( AVHRR%arraySize .lt. last_position )THEN
       last_position=last_position+1
    ELSE
       last_position=AVHRR%arraySize
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
    IF( -1 .eq. first_position )THEN
       !
       ! Gap so just add first of new one as new line
       !
       first_position = 1
    ENDIF

    !
    ! Append new data to old structure
    !
    IF( first_position .gt. 0 )THEN
       extra_lines = AVHRR_New%arraySize - first_position - &
            (AVHRR%arraySize-last_position)
       CALL Reallocate_outData(AVHRR_Ptr,extra_lines)

       !
       ! Now add in extra data - original raw data
       !
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
       ! Check against a minute and a half (rounding errors
       !
       check_overlap = .TRUE.
    ENDIF

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

  SUBROUTINE read_all_data(nFile,file1,file2,file3,file4,file5,uuid_in,&
       AVHRRout,year1,month1,day1,hour1,minute1,year2,month2,day2,&
       hour2,minute2,output_filename,walton_cal,split_single_file)

    INTEGER, INTENT(IN) :: nFile
    CHARACTER(LEN=*), INTENT(IN) :: file1
    CHARACTER(LEN=*), INTENT(IN) :: file2
    CHARACTER(LEN=*), INTENT(IN) :: file3
    CHARACTER(LEN=*), INTENT(IN) :: file4
    CHARACTER(LEN=*), INTENT(IN) :: file5
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    TYPE(AVHRR_Data), INTENT(OUT) :: AVHRRout
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

    ! Local variables
    TYPE(Imagery) :: IMG
    TYPE(AVHRR_Data), TARGET :: AVHRR
    TYPE(AVHRR_Data), TARGET :: AVHRR_Total
    TYPE(AVHRR_Instrument_Coefs) :: instr_coefs
    INTEGER :: start_pos, stop_pos
    TYPE(AVHRR_Data), POINTER :: pAVHRR=>NULL()
    LOGICAL :: out_radiances = .FALSE.
    LOGICAL :: trim_data
    INTEGER :: trim_low
    INTEGER :: trim_high
    INTEGER :: I
    TYPE(Walton_Struct) :: walton_str

    !
    ! Setup GbcsDataPath
    !
    IMG%GbcsDataPath = '/group_workspaces/cems/nceo_uor/users/jmittaz/AVHRR/code/git/gbcs_new_calibration/dat_cci/'

    AVHRR%newCalibration_there = .FALSE.
    AVHRR_Total%newCalibration_there = .FALSE.
!    instr = get_instr(file1)
    AVHRR_Total%arraySize = 0
    IF( nFile .eq. 1 )THEN
       CALL read_file(file1,AVHRR_Total,uuid_in,out_instr_coefs=instr_coefs)
       IF( .not. AVHRR_Total%valid_data_there )THEN
          RETURN
       ENDIF
    ELSE IF( nFile .eq. 2 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          CALL read_file(file2,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
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
       CALL read_file(file1,AVHRR_Total,uuid_in,out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          CALL read_file(file2,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
               'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
       ENDIF
       CALL Deallocate_OutData(AVHRR)
    ELSE IF( nfile .eq. 4 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          CALL read_file(file2,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file3,file4) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file3 and file4',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file4,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
       ENDIF
    ELSE IF( nfile .eq. 5 )THEN
       IF( .not. check_overlap(file1,file2) )THEN
          CALL Gbcs_Critical(.TRUE.,'No overlap between file1 and file2',&
               'read_all_data','extract_l1b_data.f90')
       ENDIF
       CALL read_file(file1,AVHRR_Total,uuid_in,out_instr_coefs=instr_coefs)
       IF( AVHRR_Total%valid_data_there )THEN
          CALL read_file(file2,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file2,file3) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file2 and file3',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file3,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file3,file4) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file3 and file4',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file4,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
          ENDIF
          CALL Deallocate_OutData(AVHRR)
          IF( .not. check_overlap(file4,file5) )THEN
             CALL Gbcs_Critical(.TRUE.,'No overlap between file3 and file4',&
                  'read_all_data','extract_l1b_data.f90')
          ENDIF
          CALL read_file(file5,AVHRR,uuid_in)
          IF( AVHRR%valid_data_there )THEN
             CALL merge_avhrr(AVHRR_Total,AVHRR)
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
       ! Make sure radiances are output as Marines code expects this
       !
       pAVHRR => AVHRR
       CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,out_radiances=.TRUE.,&
            moon_events=.TRUE.,&
            correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
            noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
            dig_noise=.TRUE.,all_noise=.TRUE.)
       
       IF( walton_cal )THEN
          CALL Calibration_Walton(IMG,AVHRR,.TRUE.,1.,walton_str,walton_cal)
       ENDIF
       !
       ! Resize to output
       !
       CALL Resize_Orbit(AVHRR,AVHRRout,start_pos,stop_pos,all=.TRUE.)
       CALL Deallocate_OutData(AVHRR)
       !
       ! Add in FIDUCEO uncertainties
       !
       CALL Add_FIDUCEO_Uncert(AVHRRout,uuid_in,output_filename)
       CALL Deallocate_OutData(AVHRRout)    
    ELSE
       !
       ! Calibrate whole orbit - it gets trimmed later
       !
       !
       ! Make sure radiances are output as Marines code expects this
       !
       pAVHRR => AVHRR_Total
       CALL Recalibrate_AVHRR(IMG,instr_coefs,pAVHRR,out_radiances=.TRUE.,&
            moon_events=.TRUE.,&
            correct_solar_simple=.TRUE.,new_vis_cal=.TRUE.,&
            noise_orbit=.TRUE.,filter_counts=.TRUE.,filter_prt=.TRUE.,&
            dig_noise=.TRUE.,all_noise=.TRUE.)       
       IF( walton_cal )THEN
          CALL Calibration_Walton(IMG,AVHRR_Total,.TRUE.,1.,walton_str,&
               walton_cal)
       ENDIF
       IF( split_single_file )THEN
          !
          ! Now split out sections of data
          !
          CALL Resize_Orbit_Equator(AVHRR_Total,AVHRRout,start_pos,stop_pos,&
               year1,month1,day1,hour1,minute1,&
               year1,month2,day2,hour2,minute2,&
               trim_data,trim_low,trim_high)       
          CALL Resize_Orbit(AVHRRout,AVHRR,start_pos,stop_pos,all=.TRUE.)
          !
          ! Add in FIDUCEO uncertainties
          !
          CALL Add_FIDUCEO_Uncert(AVHRR,uuid_in,output_filename)
          CALL Deallocate_OutData(AVHRR)
       ELSE
          !
          ! Output complete orbit with no equator splitting
          !
          ! Add in FIDUCEO uncertainties
          !
          CALL Add_FIDUCEO_Uncert(AVHRR_Total,uuid_in,output_filename)
          CALL Deallocate_OutData(AVHRR_Total)
       ENDIF
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
       trim_data,trim_low,trim_high,nosmooth)

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

    IF( PRESENT(nosmooth) )THEN
       no_smooth = nosmooth
    ELSE
       no_smooth = .FALSE.
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

    !
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
    !     not make_orbit1/not make_orbit2
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
    ELSE IF( make_orbit1 .and. .not. make_orbit2 )THEN
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
       out_stop_pos = stop_pos
       start_pos = MAX(1,start_pos-NPIXEL_PRT_SMOOTH)
       stop_pos = MIN(AVHRR%arraySize,stop_pos+NPIXEL_PRT_SMOOTH)
       if( trim_data )then
          out_start_pos = trim_low
          out_stop_pos = trim_high
       else
          out_start_pos = (out_start_pos-start_pos)+1
          out_stop_pos = out_start_pos + nlines - 1
       endif
    ELSE
       out_start_pos = start_pos
       out_stop_pos = stop_pos
    ENDIF

    !
    ! Now resize to do recalibration
    !
    CALL Resize_Orbit(AVHRR,AVHRRout,start_pos,stop_pos)

  END SUBROUTINE Resize_Orbit_Equator

  SUBROUTINE Resize_Orbit(AVHRR,AVHRRout,startpos,endpos,all)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(AVHRR_Data), INTENT(OUT), TARGET :: AVHRRout
    INTEGER, INTENT(IN) :: startpos
    INTEGER, INTENT(IN) :: endpos
    LOGICAL, INTENT(IN), OPTIONAL :: all

    ! Local variables
    INTEGER :: I,J,K
    INTEGER :: nsize
    INTEGER :: STAT
    LOGICAL :: alldata
    TYPE(AVHRR_Data), POINTER :: pAVHRRout
    INTEGER :: mem_scale

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
    CALL Allocate_OutData(AVHRR%nelem,pAVHRRout)
    nsize = endpos-startpos+1
    IF( nsize .gt. MEMORY_ALLOC_STEP )THEN
       CALL Gbcs_Warning(.TRUE.,&
            'Having to increase memory usage in output structure',&
            'Resize_Orbit','combine_orbits.f90')
       mem_scale = (nsize/MEMORY_ALLOC_STEP)*MEMORY_ALLOC_STEP
       CALL Reallocate_OutData(pAVHRRout,mem_scale)
    ENDIF
    CALL Reallocate_Final_outData(pAVHRRout,nsize)

    IF( alldata .and. AVHRR%newCalibration_There )THEN

       AVHRRout%newCalibration_There = AVHRR%newCalibration_There
       AVHRRout%orbital_temperature = AVHRR%orbital_temperature
       AVHRRout%new_cal_coefs3 = AVHRRout%new_cal_coefs3
       AVHRRout%new_cal_coefs4 = AVHRRout%new_cal_coefs4
       AVHRRout%new_cal_coefs5 = AVHRRout%new_cal_coefs5
       AVHRRout%gain_stdev = AVHRRout%gain_stdev
       AVHRRout%earthshine_eta = AVHRRout%earthshine_eta
       AVHRRout%poly_coefs3 = AVHRRout%poly_coefs3
       AVHRRout%poly_coefs4 = AVHRRout%poly_coefs4
       AVHRRout%poly_coefs5 = AVHRRout%poly_coefs5
       AVHRRout%gain_maxdev = AVHRRout%gain_maxdev

       ALLOCATE(AVHRRout%new_calib3(3,AVHRRout%arraySize),&
            AVHRRout%new_calib4(3,AVHRRout%arraySize),&
            AVHRRout%new_calib5(3,AVHRRout%arraySize),&
            AVHRRout%smoothPrt1(AVHRRout%arraySize),&
            AVHRRout%smoothPrt2(AVHRRout%arraySize),&
            AVHRRout%smoothPrt3(AVHRRout%arraySize),&
            AVHRRout%smoothPrt4(AVHRRout%arraySize),&
            AVHRRout%nsmoothPrt1(AVHRRout%arraySize),&
            AVHRRout%nsmoothPrt2(AVHRRout%arraySize),&
            AVHRRout%nsmoothPrt3(AVHRRout%arraySize),&
            AVHRRout%nsmoothPrt4(AVHRRout%arraySize),&
            AVHRRout%smoothPrt1Cnts(AVHRRout%arraySize),&
            AVHRRout%smoothPrt2Cnts(AVHRRout%arraySize),&
            AVHRRout%smoothPrt3Cnts(AVHRRout%arraySize),&
            AVHRRout%smoothPrt4Cnts(AVHRRout%arraySize),&
            AVHRRout%smoothPrt(AVHRRout%arraySize),&
            AVHRRout%smoothBB3(AVHRRout%arraySize),&
            AVHRRout%smoothBB4(AVHRRout%arraySize),&
            AVHRRout%smoothBB5(AVHRRout%arraySize),&
            AVHRRout%smoothSp3(AVHRRout%arraySize),&
            AVHRRout%smoothSp4(AVHRRout%arraySize),&
            AVHRRout%smoothSp5(AVHRRout%arraySize),&
            AVHRRout%nsmoothBB3(AVHRRout%arraySize),&
            AVHRRout%nsmoothBB4(AVHRRout%arraySize),&
            AVHRRout%nsmoothBB5(AVHRRout%arraySize),&
            AVHRRout%nsmoothSp3(AVHRRout%arraySize),&
            AVHRRout%nsmoothSp4(AVHRRout%arraySize),&
            AVHRRout%nsmoothSp5(AVHRRout%arraySize),&
            AVHRRout%Interpolated(AVHRRout%arraySize),&
            AVHRRout%solar_contamination_failure(AVHRRout%arraySize),&
            AVHRRout%solar_contamination_3B(AVHRRout%arraySize),&
            AVHRRout%solar_contamination_4(AVHRRout%arraySize),&
            AVHRRout%solar_contamination_5(AVHRRout%arraySize),&
            AVHRRout%moon_contamination(AVHRRout%arraySize),&
            AVHRRout%new_array1(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array2(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array3A(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array3B(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array4(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array5(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array1_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array2_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array3A_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array3B_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array4_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%new_array5_error(AVHRR%nelem,AVHRRout%arraySize),&
            AVHRRout%noise_cnts(6,AVHRRout%arraySize),&
            AVHRRout%noise_cnts_cal(6,AVHRRout%arraySize),&
            STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL Gbcs_Critical(.TRUE.,'Allocating outputData (recal)',&
               'Resize_Orbits','combine_orbits.f90')
       ENDIF
       IF( AVHRR%walton_there )THEN
          ALLOCATE(AVHRRout%array3B_error(AVHRR%nelem,AVHRRout%arraySize),&
               AVHRRout%array4_error(AVHRR%nelem,AVHRRout%arraySize),&
               AVHRRout%array5_error(AVHRR%nelem,AVHRRout%arraySize),STAT=STAT)

          IF( 0 .ne. STAT )THEN
             CALL Gbcs_Critical(.TRUE.,'Allocating outputData (walton)',&
                  'Resize_Orbits','combine_orbits.f90')
          ENDIF
       ENDIF
    ENDIF

    AVHRRout%isGAC = AVHRR%isGAC
    AVHRRout%dataFilled = AVHRR%dataFilled
    AVHRRout%AVHRR_No = AVHRR%AVHRR_No
    AVHRRout%nelem = AVHRR%nelem
    AVHRRout%filter3a = AVHRR%filter3a
    AVHRRout%start_valid = 1
    AVHRRout%stop_valid = (endpos-startpos)+1
    AVHRRout%valid_data_there = AVHRR%valid_data_there
    AVHRRout%walton_there = AVHRR%walton_there
    K = 0
    DO I=startpos,endpos
       K=K+1
       AVHRRout%scanLineNumber(K) = AVHRR%scanLineNumber(I)
       AVHRRout%badTop(K) = AVHRR%badTop(I)
       AVHRRout%badTime(K) = AVHRR%badTime(I)
       AVHRRout%badNavigation(K) = AVHRR%badNavigation(I)
       AVHRRout%badCalibration(K) = AVHRR%badCalibration(I)
       AVHRRout%transition3A3B(K) = AVHRR%transition3A3B(I)
       AVHRRout%Lon(:,K) = AVHRR%Lon(:,I)
       AVHRRout%Lat(:,K) = AVHRR%Lat(:,I)
       AVHRRout%satZA(:,K) = AVHRR%satZA(:,I)
       AVHRRout%solZA(:,K) = AVHRR%solZA(:,I)
       AVHRRout%relAz(:,K) = AVHRR%relAz(:,I)
       AVHRRout%Counts1(:,K) = AVHRR%Counts1(:,I)
       AVHRRout%Counts2(:,K) = AVHRR%Counts2(:,I)
       AVHRRout%Counts3(:,K) = AVHRR%Counts3(:,I)
       AVHRRout%Counts4(:,K) = AVHRR%Counts4(:,I)
       AVHRRout%Counts5(:,K) = AVHRR%Counts5(:,I)
       AVHRRout%array1(:,K) = AVHRR%array1(:,I)
       AVHRRout%array2(:,K) = AVHRR%array2(:,I)
       AVHRRout%array3A(:,K) = AVHRR%array3A(:,I)
       AVHRRout%array3B(:,K) = AVHRR%array3B(:,I)
       AVHRRout%array4(:,K) = AVHRR%array4(:,I)
       AVHRRout%array5(:,K) = AVHRR%array5(:,I)
       AVHRRout%year(K) = AVHRR%year(I)
       AVHRRout%month(K) = AVHRR%month(I)
       AVHRRout%day(K) = AVHRR%day(I)
       AVHRRout%dayno(K) = AVHRR%dayno(I)
       AVHRRout%hours(K) = AVHRR%hours(I)
       AVHRRout%UTC_msecs(K) = AVHRR%UTC_msecs(I)
       AVHRRout%time(K) = AVHRR%time(I)
       AVHRRout%prt1(K) = AVHRR%prt1(I)
       AVHRRout%prt2(K) = AVHRR%prt2(I)
       AVHRRout%prt3(K) = AVHRR%prt3(I)
       AVHRRout%prt4(K) = AVHRR%prt4(I)
       AVHRRout%prt1Counts(K) = AVHRR%prt1Counts(I)
       AVHRRout%prt2Counts(K) = AVHRR%prt2Counts(I)
       AVHRRout%prt3Counts(K) = AVHRR%prt3Counts(I)
       AVHRRout%prt4Counts(K) = AVHRR%prt4Counts(I)
       AVHRRout%prt1CountsAll(:,K) = AVHRR%prt1CountsAll(:,I)
       AVHRRout%prt2CountsAll(:,K) = AVHRR%prt2CountsAll(:,I)
       AVHRRout%prt3CountsAll(:,K) = AVHRR%prt3CountsAll(:,I)
       AVHRRout%prt4CountsAll(:,K) = AVHRR%prt4CountsAll(:,I)
       AVHRRout%bb3(K) = AVHRR%bb3(I)
       AVHRRout%bb4(K) = AVHRR%bb4(I)
       AVHRRout%bb5(K) = AVHRR%bb5(I)
       AVHRRout%sp3(K) = AVHRR%sp3(I)
       AVHRRout%sp4(K) = AVHRR%sp4(I)
       AVHRRout%sp5(K) = AVHRR%sp5(I)
       AVHRRout%bbodyFilter3(:,K) = AVHRR%bbodyFilter3(:,I)
       AVHRRout%bbodyFilter4(:,K) = AVHRR%bbodyFilter4(:,I)
       AVHRRout%bbodyFilter5(:,K) = AVHRR%bbodyFilter5(:,I)
       AVHRRout%spaceFilter1(:,K) = AVHRR%spaceFilter1(:,I)
       AVHRRout%spaceFilter2(:,K) = AVHRR%spaceFilter2(:,I)
       AVHRRout%spaceFilter3a(:,K) = AVHRR%spaceFilter3a(:,I)
       AVHRRout%spaceFilter3(:,K) = AVHRR%spaceFilter3(:,I)
       AVHRRout%spaceFilter4(:,K) = AVHRR%spaceFilter4(:,I)
       AVHRRout%spaceFilter5(:,K) = AVHRR%spaceFilter5(:,I)
       AVHRRout%patch(K) = AVHRR%patch(I)
       AVHRRout%patchExtended(K) = AVHRR%patchExtended(I)
       AVHRRout%Radiator(K) = AVHRR%Radiator(I)
       AVHRRout%Cooler(K) = AVHRR%Cooler(I)
       AVHRRout%a_d_conv(K) = AVHRR%a_d_conv(I)
       AVHRRout%motor(K) = AVHRR%motor(I)
       AVHRRout%motorCurrent(K) = AVHRR%motorCurrent(I)
       AVHRRout%electronics(K) = AVHRR%electronics(I)
       AVHRRout%baseplate(K) = AVHRR%baseplate(I)
       AVHRRout%calib1(:,K) = AVHRR%calib1(:,I)
       AVHRRout%calib1_2(:,K) = AVHRR%calib1_2(:,I)
       AVHRRout%calib1_intercept(K) = AVHRR%calib1_intercept(I)
       AVHRRout%calib2(:,K) = AVHRR%calib2(:,I)
       AVHRRout%calib2_2(:,K) = AVHRR%calib2_2(:,I)
       AVHRRout%calib2_intercept(K) = AVHRR%calib2_intercept(I)
       AVHRRout%calib3A(:,K) = AVHRR%calib3A(:,I)
       AVHRRout%calib3A_2(:,K) = AVHRR%calib3A_2(:,I)
       AVHRRout%calib3A_intercept(K) = AVHRR%calib3A_intercept(I)
       AVHRRout%calib3(:,K) = AVHRR%calib3(:,I)
       AVHRRout%calib4(:,K) = AVHRR%calib4(:,I)
       AVHRRout%calib5(:,K) = AVHRR%calib5(:,I)
       AVHRRout%clavr_mask(:,K) = AVHRR%clavr_mask(:,I)
       AVHRRout%clavrx_mask(:,K) = AVHRR%clavrx_mask(:,I)
       AVHRRout%clavrx_prb(:,K) = AVHRR%clavrx_prb(:,I)
       AVHRRout%orig_solar_contamination_3B(K) = &
            AVHRR%orig_solar_contamination_3B(I)
       AVHRRout%orig_solar_contamination_4(K) = &
            AVHRR%orig_solar_contamination_4(I)
       AVHRRout%orig_solar_contamination_5(K) = &
            AVHRR%orig_solar_contamination_5(I)
       AVHRRout%satelliteAltitude(K) = AVHRR%satelliteAltitude(I)

       IF( AVHRR%walton_there )THEN
          AVHRRout%array3B_error(:,K) = AVHRR%array3B_error(:,I)
          AVHRRout%array4_error(:,K) = AVHRR%array4_error(:,I)
          AVHRRout%array5_error(:,K) = AVHRR%array5_error(:,I)
       ENDIF
       
       !
       ! If new calbration there
       !
       IF( alldata .and. AVHRR%newCalibration_There )THEN

          AVHRRout%new_calib3(:,K) = AVHRR%new_calib3(:,I)
          AVHRRout%new_calib4(:,K) = AVHRR%new_calib4(:,I)
          AVHRRout%new_calib5(:,K) = AVHRR%new_calib5(:,I)
          AVHRRout%smoothPrt1(K) = AVHRR%smoothPrt1(I)
          AVHRRout%smoothPrt2(K) = AVHRR%smoothPrt2(I)
          AVHRRout%smoothPrt3(K) = AVHRR%smoothPrt3(I)
          AVHRRout%smoothPrt4(K) = AVHRR%smoothPrt4(I)
          AVHRRout%nsmoothPrt1(K) = AVHRR%nsmoothPrt1(I)
          AVHRRout%nsmoothPrt2(K) = AVHRR%nsmoothPrt2(I)
          AVHRRout%nsmoothPrt3(K) = AVHRR%nsmoothPrt3(I)
          AVHRRout%nsmoothPrt4(K) = AVHRR%nsmoothPrt4(I)
          AVHRRout%smoothPrt1Cnts(K) = AVHRR%smoothPrt1Cnts(I)
          AVHRRout%smoothPrt2Cnts(K) = AVHRR%smoothPrt2Cnts(I)
          AVHRRout%smoothPrt3Cnts(K) = AVHRR%smoothPrt3Cnts(I)
          AVHRRout%smoothPrt4Cnts(K) = AVHRR%smoothPrt4Cnts(I)
          AVHRRout%smoothPrt(K) = AVHRR%smoothPrt(I)
          AVHRRout%smoothBB3(K) = AVHRR%smoothBB3(I)
          AVHRRout%smoothBB4(K) = AVHRR%smoothBB4(I)
          AVHRRout%smoothBB5(K) = AVHRR%smoothBB5(I)
          AVHRRout%smoothSp3(K) = AVHRR%smoothSp3(I)
          AVHRRout%smoothSp4(K) = AVHRR%smoothSp4(I)
          AVHRRout%smoothSp5(K) = AVHRR%smoothSp5(I)
          AVHRRout%nsmoothBB3(K) = AVHRR%nsmoothBB3(I)
          AVHRRout%nsmoothBB4(K) = AVHRR%nsmoothBB4(I)
          AVHRRout%nsmoothBB5(K) = AVHRR%nsmoothBB5(I)
          AVHRRout%nsmoothSp3(K) = AVHRR%nsmoothSp3(I)
          AVHRRout%nsmoothSp4(K) = AVHRR%nsmoothSp4(I)
          AVHRRout%nsmoothSp5(K) = AVHRR%nsmoothSp5(I)
          AVHRRout%Interpolated(K) = AVHRR%Interpolated(I)
          AVHRRout%solar_contamination_failure(K) = &
               AVHRR%solar_contamination_failure(I)
          AVHRRout%solar_contamination_3B(K) = AVHRR%solar_contamination_3B(I)
          AVHRRout%solar_contamination_4(K) = AVHRR%solar_contamination_4(I)
          AVHRRout%solar_contamination_5(K) = AVHRR%solar_contamination_5(I)
          AVHRRout%moon_contamination(K) = AVHRR%moon_contamination(I)
          AVHRRout%new_array1(:,K) = AVHRR%new_array1(:,I)
          AVHRRout%new_array2(:,K) = AVHRR%new_array2(:,I)
          AVHRRout%new_array3A(:,K) = AVHRR%new_array3A(:,I)
          AVHRRout%new_array3B(:,K) = AVHRR%new_array3B(:,I)
          AVHRRout%new_array4(:,K) = AVHRR%new_array4(:,I)
          AVHRRout%new_array5(:,K) = AVHRR%new_array5(:,I)
          AVHRRout%new_array1_error(:,K) = AVHRR%new_array1_error(:,I)
          AVHRRout%new_array2_error(:,K) = AVHRR%new_array2_error(:,I)
          AVHRRout%new_array3A_error(:,K) = AVHRR%new_array3A_error(:,I)
          AVHRRout%new_array3B_error(:,K) = AVHRR%new_array3B_error(:,I)
          AVHRRout%new_array4_error(:,K) = AVHRR%new_array4_error(:,I)
          AVHRRout%new_array5_error(:,K) = AVHRR%new_array5_error(:,I)
          AVHRRout%noise_cnts(:,K) = AVHRR%noise_cnts(:,I)
          AVHRRout%noise_cnts_cal(:,K) = AVHRR%noise_cnts_cal(:,I)

       ENDIF

    END DO

  END SUBROUTINE Resize_Orbit

  SUBROUTINE read_file(infile,AVHRR,uuid_in,out_instr_coefs)

    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE(AVHRR_Data), INTENT(OUT) :: AVHRR
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    TYPE(AVHRR_Instrument_Coefs), INTENT(out), OPTIONAL :: out_instr_coefs

    ! Local variables
    INTEGER :: pid
    INTEGER :: stat
    INTEGER :: POS

    CHARACTER(LEN=256) :: inDirectory
    CHARACTER(LEN=256) :: infilename
    CHARACTER(LEN=256) :: temp_file
    CHARACTER(LEN=256) :: command_str

    TYPE(Imagery) :: IMG
    REAL :: coefs1(7)
    REAL :: coefs2(7)
    REAL :: coefs3(7)
    LOGICAL :: remove_file

    !
    ! Check to see if we need to uncompress the data
    !
    remove_file=.FALSE.
    IF( 0 .ne. INDEX(infile,'.gz') )THEN
       WRITE(command_str,'(''cp -f '',a,'' temp_file.'',a,''.gz'')')&
            TRIM(infile),TRIM(uuid_in)
       CALL SYSTEM(command_str,STATUS=stat)
       CALL Gbcs_Critical(stat.ne.0,'Cannot copy file to tempfile','Convert',&
            'extract_l1b_data.f90')
       WRITE(command_str,'(''gunzip -f temp_file.'',a,''.gz'')')&
            TRIM(uuid_in)
       CALL SYSTEM(command_str,STATUS=stat)
       CALL Gbcs_Critical(stat.ne.0,'Cannot gunzip tempfile','Convert',&
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

    CALL Load_Imagery(IMG,outputData=AVHRR,use_new_calibration=.FALSE.,&
         out_instr_coefs=out_instr_coefs,use_walton=.FALSE.)

    IF( remove_file )THEN
       WRITE(command_str,'(''rm -f temp_file.'',a)')&
            TRIM(uuid_in)
       CALL SYSTEM(command_str,STATUS=stat)
       CALL Gbcs_Critical(stat.ne.0,'Cannot remove tempfile','Convert',&
            'extract_l1b_data.f90')
    ENDIF

  END SUBROUTINE read_file
  
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

END MODULE Combine_Orbits

