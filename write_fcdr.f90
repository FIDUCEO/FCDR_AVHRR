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
! * ---------------------------------------------------------------------------

! * MT: 16-11-2017: fix for excess nArgs due to repeat orbits now applied in equator_to_equator.py

PROGRAM Extract_L1b_Data

  USE GbcsErrorHandler
  USE NOAA_LoadAVHRRLevel1B 
  USE Combine_Orbits

  IMPLICIT NONE

  CALL TopLevel()

CONTAINS

  SUBROUTINE TopLevel()

    INTEGER :: nArgs
    INTEGER :: nFiles
    INTEGER :: stat
    CHARACTER(LEN=256) :: file1
    CHARACTER(LEN=256) :: file2
    CHARACTER(LEN=256) :: file3
    CHARACTER(LEN=256) :: file4
    CHARACTER(LEN=256) :: file5
    CHARACTER(LEN=256) :: ofile
    CHARACTER(LEN=256) :: uuid_in
    CHARACTER(LEN=256) :: instr_str
    CHARACTER(LEN=256) :: y1
    CHARACTER(LEN=256) :: m1
    CHARACTER(LEN=256) :: d1
    CHARACTER(LEN=256) :: h1
    CHARACTER(LEN=256) :: min1
    CHARACTER(LEN=256) :: y2
    CHARACTER(LEN=256) :: m2
    CHARACTER(LEN=256) :: d2
    CHARACTER(LEN=256) :: h2
    CHARACTER(LEN=256) :: min2
    CHARACTER(LEN=256) :: sngle_split
    INTEGER :: year1, month1, day1, hour1, minute1
    INTEGER :: year2, month2, day2, hour2, minute2
    TYPE(AVHRR_Data) :: AVHRR
    INTEGER :: instr
    TYPE(AVHRR_Radiative_Coefs) :: avhrr_rad_coefs
    CHARACTER(LEN=1) :: first_segment_input
    LOGICAL :: first_segment
    LOGICAL :: new_filename
    LOGICAL :: split_single_file

!MT: 26-10-2017: extract list of input arguments: 
!    INTEGER :: i
!    CHARACTER(LEN=256) :: temp
!    nArgs = COMMAND_ARGUMENT_COUNT()
!    write(*,*)'nArgs=',nArgs
!    DO i = 1,nArgs
!        CALL GET_COMMAND_ARGUMENT(i,temp,STATUS=stat)
!        write(*,*),'inArgs=',temp
!    ENDDO

    nArgs = COMMAND_ARGUMENT_COUNT()
    IF( 14 .gt. nArgs .or. 19 .le. nArgs )THEN
       write(*,*)'nArgs=',nArgs
       CALL Gbcs_Critical(.TRUE.,&
            'USAGE: ./extract_l1b_data.exe uuid outfile eq_year1 eq_month1 &
            &eq_day1 eq_hour1 eq_min1 eq_year2 eq_month2 eq_day2 eq_hour2 &
            &eq_min2 split_single file1 (file2) (file3) (file4) (file5)',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(1,uuid_in,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get first command line argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(2,ofile,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get instr argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(3,y1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get year1 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(4,m1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get month1 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(5,d1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get day1 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(6,h1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get hour1 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(7,min1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get minute1 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(8,y2,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get year2 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(9,m2,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get month2 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(10,d2,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get day2 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(11,h2,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get hour2 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(12,min2,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get minute2 argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(13,sngle_split,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get sngle_split argument',&
            'Main','extract_l1b_data.f90')
    ENDIF

    CALL GET_COMMAND_ARGUMENT(14,file1,STATUS=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get file1 command line argument',&
            'Main','extract_l1b_data.f90')
    ENDIF
    nfiles = 1

    IF( nArgs .gt. 14 )THEN
       CALL GET_COMMAND_ARGUMENT(15,file2,STATUS=stat)
       IF( 0 .ne. stat )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot get file2 command line argument',&
               'Main','extract_l1b_data.f90')
       ENDIF
       nfiles = 2
       IF( nArgs .gt. 15 )THEN
          CALL GET_COMMAND_ARGUMENT(16,file3,STATUS=stat)
          IF( 0 .ne. stat )THEN
             CALL Gbcs_Critical(.TRUE.,'Cannot get file3 command line argument',&
                  'Main','extract_l1b_data.f90')
          ENDIF
          nfiles = 3
          IF( nArgs .gt. 16 )THEN
             CALL GET_COMMAND_ARGUMENT(17,file4,STATUS=stat)
             IF( 0 .ne. stat )THEN
                CALL Gbcs_Critical(.TRUE.,'Cannot get file4 command line argument',&
                     'Main','extract_l1b_data.f90')
             ENDIF
             nfiles = 4
             IF( nArgs .gt. 17 )THEN
                CALL GET_COMMAND_ARGUMENT(18,file5,STATUS=stat)
                IF( 0 .ne. stat )THEN
                   CALL Gbcs_Critical(.TRUE.,'Cannot get file5 command line argument',&
                        'Main','extract_l1b_data.f90')
                ENDIF
                nfiles = 5
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    READ(y1,'(i4)')year1
    READ(m1,'(i4)')month1
    READ(d1,'(i4)')day1
    READ(h1,'(i4)')hour1
    READ(min1,'(i4)')minute1

    READ(y2,'(i4)')year2
    READ(m2,'(i4)')month2
    READ(d2,'(i4)')day2
    READ(h2,'(i4)')hour2
    READ(min2,'(i4)')minute2

    IF( 'Y' .eq. sngle_split .or. 'y' .eq. sngle_split )THEN
       split_single_file = .TRUE.
    ELSE
       split_single_file = .FALSE.
    ENDIF

    new_filename = .FALSE.
    IF( ofile .eq. 'TIROSN' )then
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA06' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA07' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA08' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA09' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA10' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA11' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA12' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA14' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA15' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA16' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA17' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA18' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'NOAA19' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'METOPA' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ELSE IF( ofile .eq. 'METOPB' )THEN
       instr_str = ofile
       new_filename = .TRUE.
    ENDIF
!FIDUCEO_FCDR_L1C_AVHRR1_NOAA08_19840106013100_19840106E03200_EASY_v0.1_fv0.1.nc
    IF( new_filename )THEN
!       WRITE(ofile,'(''FIDUCEO_FCDR_L1C_'',a,''_'',i4.4,i2.2,i2.2,i2.2,i2.2,&
!            &i2.2,''_'',i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,''_EASY_v0.1_fv01.nc'')')&
!            TRIM(instr_str),year1,month1,day1,hour1,minute1,0,&
!            year2,month2,day2,hour2,minute2,0
       ofile = 'None'
    ENDIF
    !
    ! Note walton calibration for pre-beta
    !
    CALL read_all_data(nfiles,file1,file2,file3,file4,file5,uuid_in,&
         AVHRR,year1,month1,day1,hour1,minute1,year2,month2,day2,hour2,&
         minute2,ofile,.TRUE.,split_single_file)
    !
    ! Deallocate structure
    !
    CALL Deallocate_OutData(AVHRR)

  END SUBROUTINE TopLevel

END PROGRAM EXTRACT_L1B_DATA

