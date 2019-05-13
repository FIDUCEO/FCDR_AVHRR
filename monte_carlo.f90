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
! Module to get montecarlo inputs to calibration and run calibration
! output delta_BTs which still need to be stripped
!
MODULE Monte_Carlo

  USE GbcsKinds
  USE GbcsTypes
  USE NOAA_LoadAVHRRLevel1B
  USE GbcsErrorHandler
  USE fiduceo_uncertainties
  USE sampled_gaussian
  USE fiduceo_calibration

  IMPLICIT NONE

  INTERFACE Get_Random_EVEN
     MODULE PROCEDURE Get_Real_Random_EVEN
     MODULE PROCEDURE Get_Integer_Random_EVEN
  END INTERFACE

  INTERFACE Get_Random_ODD
     MODULE PROCEDURE Get_Real_Random_ODD
     MODULE PROCEDURE Get_Integer_Random_ODD
  END INTERFACE

  PRIVATE
  PUBLIC :: Run_MonteCarlo
  PUBLIC :: Allocate_DeltaRad
  PUBLIC :: Deallocate_DeltaRad

CONTAINS

  SUBROUTINE Allocate_DeltaRad(nelem,nscan,noMC,delta_rad)

    INTEGER, INTENT(IN) :: nelem
    INTEGER, INTENT(IN) :: nscan
    INTEGER, INTENT(IN) :: noMC
    TYPE(mc_delta_str), INTENT(INOUT) :: delta_rad

    ! Local variables
    INTEGER :: STAT
    
    delta_rad%nelem = nelem
    delta_rad%nscan = nscan
    delta_rad%nmc = noMC
    ALLOCATE(delta_rad%ch1(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         delta_rad%ch2(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         delta_rad%ch3a(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         delta_rad%ch3(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         delta_rad%ch4(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         delta_rad%ch5(delta_rad%nelem,delta_rad%nscan,delta_rad%nmc),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate delta_rad',&
            'Allocate_DeltaRad','monte_carlo.f90')
    ENDIF
    delta_rad%ch1 = NAN_R
    delta_rad%ch2 = NAN_R
    delta_rad%ch3a = NAN_R
    delta_rad%ch3 = NAN_R
    delta_rad%ch4 = NAN_R
    delta_rad%ch5 = NAN_R

  END SUBROUTINE Allocate_DeltaRad

  SUBROUTINE Deallocate_DeltaRad(delta_rad)

    TYPE(mc_delta_str), INTENT(INOUT) :: delta_rad

    DEALLOCATE(delta_rad%ch1,delta_rad%ch2,delta_rad%ch3a,&
         delta_rad%ch3,delta_rad%ch4,delta_rad%ch5)

  END SUBROUTINE Deallocate_DeltaRad

  SUBROUTINE Run_MonteCarlo(harm_mc,IMG,instr_coefs,AVHRR,AVHRR_MC,&
       out_radiances,walton_str,&
       delta_rad,moon_events,use_old_solar,&
       correct_solar_simple,new_vis_cal,&
       noise_orbit,filter_counts,filter_prt,&
       dig_noise,all_noise,&
       correctict,applyict,walton,&
       tinstrapply,output_solar_temp,&
       bad_tiny,digitize)

    TYPE(mc_harm_str), INTENT(IN) :: harm_mc
    TYPE(Imagery), INTENT(IN) :: IMG
    TYPE(AVHRR_Instrument_Coefs), INTENT(in) :: instr_coefs
    TYPE(AVHRR_Data), POINTER :: AVHRR
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR_MC
    LOGICAL, INTENT(IN) :: out_radiances
    TYPE(Walton_Struct), INTENT(IN) :: walton_str
    TYPE(mc_delta_str), INTENT(OUT) :: delta_rad
    LOGICAL, INTENT(IN), OPTIONAL :: moon_events
    LOGICAL, INTENT(IN), OPTIONAL :: correct_solar_simple
    LOGICAL, INTENT(IN), OPTIONAL :: new_vis_cal
    LOGICAL, INTENT(IN), OPTIONAL :: use_old_solar
    LOGICAL, INTENT(IN), OPTIONAL :: noise_orbit
    LOGICAL, INTENT(IN), OPTIONAL :: filter_counts
    LOGICAL, INTENT(IN), OPTIONAL :: filter_prt
    LOGICAL, INTENT(IN), OPTIONAL :: all_noise
    LOGICAL, INTENT(IN), OPTIONAL :: dig_noise
    LOGICAL, INTENT(IN), OPTIONAL :: correctict
    LOGICAL, INTENT(IN), OPTIONAL :: applyict
    LOGICAL, INTENT(IN), OPTIONAL :: tinstrapply
    LOGICAL, INTENT(IN), OPTIONAL :: walton
    INTEGER, INTENT(IN), OPTIONAL :: output_solar_temp
    LOGICAL, INTENT(IN), OPTIONAL :: bad_tiny
    LOGICAL, INTENT(IN), OPTIONAL :: digitize
    
    ! Local variables
    INTEGER :: noMC
    INTEGER :: I,J,K
    INTEGER :: STAT
    REAL :: prt_uncert
    REAL, ALLOCATABLE :: prt_adjust(:)
    REAL, ALLOCATABLE :: random_array(:)
    TYPE(AVHRR_Data), TARGET :: newAVHRR
    TYPE(AVHRR_Data), POINTER :: pnewAVHRR
    REAL, ALLOCATABLE :: Counts1(:,:,:)
    REAL, ALLOCATABLE :: Counts2(:,:,:)
    REAL, ALLOCATABLE :: Counts3(:,:,:)
    REAL, ALLOCATABLE :: Counts4(:,:,:)
    REAL, ALLOCATABLE :: Counts5(:,:,:)
    REAL, ALLOCATABLE :: bb3(:,:,:)
    REAL, ALLOCATABLE :: bb4(:,:,:)
    REAL, ALLOCATABLE :: bb5(:,:,:)
    REAL, ALLOCATABLE :: sp1(:,:,:)
    REAL, ALLOCATABLE :: sp2(:,:,:)
    REAL, ALLOCATABLE :: sp3a(:,:,:)
    REAL, ALLOCATABLE :: sp3(:,:,:)
    REAL, ALLOCATABLE :: sp4(:,:,:)
    REAL, ALLOCATABLE :: sp5(:,:,:)
    LOGICAL :: digitize_data
    REAL :: noise_cnts(6)
    LOGICAL :: ok
    
    IF( PRESENT(digitize) )THEN
       digitize_data = digitize
    ELSE
       digitize_data = .TRUE.
    ENDIF

    noMC = harm_mc%noMC
 
    write(*,'('' Running Monte-Carlo'')')
    !
    ! Get new values with MC perturabation
    !
    ok=.FALSE.
    noiseLoop: DO I=1,AVHRR%arraySize
       !
       ! Check channels which should always be there
       !
       IF( AVHRR%noise_cnts(1,I) .gt. 0 .and. &
            AVHRR%noise_cnts(2,I) .gt. 0 .and. &
            AVHRR%noise_cnts(4,I) .gt. 0 .and. &
            AVHRR%noise_cnts(5,I) .gt. 0 )THEN
          noise_cnts(:) = AVHRR%noise_cnts(:,I)
          ok=.TRUE.
          EXIT noiseLoop
       ENDIF
    END DO noiseLoop
    IF( .not. ok )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot find good noise',&
            'Run_MonteCarlo','monte_carlo.f90')
    ENDIF
    CALL Get_MC_Data(AVHRR_MC%arraySize,AVHRR_MC%nelem,noMC,AVHRR_MC%ncal,&
         AVHRR_MC%Counts1,AVHRR_MC%Counts2,AVHRR_MC%Counts3,&
         AVHRR_MC%Counts4,AVHRR_MC%Counts5,AVHRR_MC%bb3,&
         AVHRR_MC%bb4,AVHRR_MC%bb5,AVHRR_MC%bbodyFilter3,&
         AVHRR_MC%bbodyFilter4,AVHRR_MC%bbodyFilter5,AVHRR_MC%spaceFilter1,&
         AVHRR_MC%spaceFilter2,AVHRR_MC%spaceFilter3a,AVHRR_MC%sp3,&
         AVHRR_MC%sp4,AVHRR_MC%sp5,&
         AVHRR_MC%spaceFilter3,AVHRR_MC%spaceFilter4,AVHRR_MC%spaceFilter5,&
         noise_cnts,AVHRR%ch3a_there,&
         Counts1,Counts2,Counts3,Counts4,Counts5,bb3,bb4,bb5,&
         sp1,sp2,sp3a,sp3,sp4,sp5,digitize_data)
    !
    ! Get ICT temperature offset
    !
    ALLOCATE(prt_adjust(noMC),random_array(noMC),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate prt_adjust',&
            'Run_MonteCarlo','monte_carlo.f90')
    ENDIF
    prt_uncert = sqrt(prt_accuracy**2+AVHRR%ict_plane_uncert**2)
    IF( MOD(noMC,2) .eq. 0 )THEN
       CALL Get_Random_EVEN(noMC,0.,prt_uncert,prt_adjust,random_array,&
            .FALSE.)
    ELSE
       CALL Get_Random_ODD(noMC,0.,prt_uncert,prt_adjust,random_array,&
            .FALSE.)
    ENDIF

    !
    ! Allocate delta arrays
    !
    CALL Allocate_DeltaRad(AVHRR%nelem,AVHRR%arraySize,noMC,delta_rad)

    !
    ! Allocate new structure
    !
    pnewAVHRR => newAVHRR
    CALL Allocate_OldCal(AVHRR,pnewAVHRR,1,AVHRR%arraySize)
    CALL Allocate_NewCal(AVHRR,newAVHRR,.TRUE.)

    !
    ! Copy over data to new Structure - only original data
    !
    CALL Copy_Top_Data(AVHRR,newAVHRR,1,AVHRR%arraySize)

    DO I=1,noMC
       !
       ! Reset everything
       !
       DO J=1,AVHRR%arraySize
          CALL Copy_All_Scan(AVHRR_MC,newAVHRR,J,J,.FALSE.)
       END DO
       !
       ! Copy over new counts
       !
       newAVHRR%Counts1 = Counts1(:,:,I)
       newAVHRR%Counts2 = Counts2(:,:,I)
       newAVHRR%Counts3 = Counts3(:,:,I)
       newAVHRR%Counts4 = Counts4(:,:,I)
       newAVHRR%Counts5 = Counts5(:,:,I)
       newAVHRR%bbodyFilter3 = bb3(:,:,I)
       newAVHRR%bbodyFilter4 = bb4(:,:,I)
       newAVHRR%bbodyFilter5 = bb5(:,:,I)
       newAVHRR%spaceFilter1 = sp1(:,:,I)
       newAVHRR%spaceFilter2 = sp2(:,:,I)
       newAVHRR%spaceFilter3a = sp3a(:,:,I)
       newAVHRR%spaceFilter3 = sp3(:,:,I)
       newAVHRR%spaceFilter4 = sp4(:,:,I)
       newAVHRR%spaceFilter5 = sp5(:,:,I)       

       !
       ! Rerun Recalibration adding in PRT offset
       !
       CALL Recalibrate_AVHRR(IMG,instr_coefs,&
            pnewAVHRR,out_radiances,walton_str,&
            moon_events=moon_events,use_old_solar=use_old_solar,&
            correct_solar_simple=correct_solar_simple,&
            new_vis_cal=new_vis_cal,&
            noise_orbit=noise_orbit,filter_counts=filter_counts,&
            filter_prt=filter_prt,&
            dig_noise=dig_noise,all_noise=all_noise,&
            correctict=correctict,applyict=applyict,walton=walton,&
            tinstrapply=tinstrapply,output_solar_temp=output_solar_temp,&
            bad_tiny=bad_tiny,montecarlo=.TRUE.,harm_mc=harm_mc,harm_index=I,&
            mc_prt_offset=prt_adjust(I))

       !
       ! Now get delta Radiances
       !
       DO J=1,AVHRR%arraySize
          DO K=1,AVHRR%nelem
             IF( newAVHRR%new_array1(K,J) .gt. 0 .and. &
                  AVHRR%new_array1(K,J) .gt. 0 )THEN
                delta_rad%ch1(K,J,I) = newAVHRR%new_array1(K,J) - &
                     AVHRR%new_array1(K,J)
             ENDIF
             IF( newAVHRR%new_array2(K,J) .gt. 0 .and. &
                  AVHRR%new_array2(K,J) .gt. 0 )THEN
                delta_rad%ch2(K,J,I) = newAVHRR%new_array2(K,J) - &
                     AVHRR%new_array2(K,J)
             ENDIF
             IF( newAVHRR%new_array3a(K,J) .gt. 0 .and. &
                  AVHRR%new_array3a(K,J) .gt. 0 )THEN
                delta_rad%ch3a(K,J,I) = newAVHRR%new_array3a(K,J) - &
                     AVHRR%new_array3a(K,J)
             ENDIF
             IF( newAVHRR%new_array3b(K,J) .gt. 0 .and. &
                  AVHRR%new_array3b(K,J) .gt. 0 )THEN
                delta_rad%ch3(K,J,I) = newAVHRR%new_array3b(K,J) - &
                     AVHRR%new_array3b(K,J)
             ENDIF
             IF( newAVHRR%new_array4(K,J) .gt. 0 .and. &
                  AVHRR%new_array4(K,J) .gt. 0 )THEN
                delta_rad%ch4(K,J,I) = newAVHRR%new_array4(K,J) - &
                     AVHRR%new_array4(K,J)
             ENDIF
             IF( newAVHRR%new_array5(K,J) .gt. 0 .and. &
                  AVHRR%new_array5(K,J) .gt. 0 )THEN
                delta_rad%ch5(K,J,I) = newAVHRR%new_array5(K,J) - &
                     AVHRR%new_array5(K,J)
             ENDIF
          END DO
       END DO
       
    END DO
    DEALLOCATE(prt_adjust,Counts1,Counts2,Counts3,Counts4,Counts5,bb3,bb4,bb5,&
         sp1,sp2,sp3a,sp3,sp4,sp5,random_array)
    CALL Deallocate_OutData(newAVHRR)

  END SUBROUTINE Run_MonteCarlo

  !
  ! Get new values with MC perturabation
  !
  SUBROUTINE Get_MC_Data(ndata,nelem,noMC,ncal,&
       InCounts1,InCounts2,InCounts3,&
       InCounts4,InCounts5,Inbb3,&
       Inbb4,Inbb5,InbbFilter3,InbbFilter4,InbbFilter5,Inspace1,&
       Inspace2,Inspace3a,Inspace3,&
       Inspace4,Inspace5,InspaceFilter3,InspaceFilter4,InspaceFilter5,&
       noise,ch3a_there,&
       Counts1,Counts2,Counts3,Counts4,Counts5,bb3,bb4,bb5,&
       sp1,sp2,sp3a,sp3,sp4,sp5,digitize)

    INTEGER, INTENT(IN) :: ndata
    INTEGER, INTENT(IN) :: nelem
    INTEGER, INTENT(IN) :: noMC
    INTEGER, INTENT(IN) :: ncal
    REAL, INTENT(IN) :: InCounts1(nelem,ndata)
    REAL, INTENT(IN) :: InCounts2(nelem,ndata)
    REAL, INTENT(IN) :: InCounts3(nelem,ndata)
    REAL, INTENT(IN) :: InCounts4(nelem,ndata)
    REAL, INTENT(IN) :: InCounts5(nelem,ndata)
    REAL, INTENT(IN) :: Inbb3(ndata)
    REAL, INTENT(IN) :: Inbb4(ndata)
    REAL, INTENT(IN) :: Inbb5(ndata)
    INTEGER, INTENT(IN) :: InbbFilter3(ncal,ndata)
    INTEGER, INTENT(IN) :: InbbFilter4(ncal,ndata)
    INTEGER, INTENT(IN) :: InbbFIlter5(ncal,ndata)
    INTEGER, INTENT(IN) :: Inspace1(ncal,ndata)
    INTEGER, INTENT(IN) :: Inspace2(ncal,ndata)
    INTEGER, INTENT(IN) :: Inspace3a(ncal,ndata)
    REAL, INTENT(IN) :: Inspace3(ndata)
    REAL, INTENT(IN) :: Inspace4(ndata)
    REAL, INTENT(IN) :: Inspace5(ndata)
    INTEGER, INTENT(IN) :: InspaceFilter3(ncal,ndata)
    INTEGER, INTENT(IN) :: InspaceFilter4(ncal,ndata)
    INTEGER, INTENT(IN) :: InspaceFIlter5(ncal,ndata)
    REAL, INTENT(IN) :: noise(6)
    INTEGER, INTENT(IN) :: ch3a_there(ndata)
    REAL, INTENT(OUT), ALLOCATABLE :: Counts1(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: Counts2(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: Counts3(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: Counts4(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: Counts5(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: bb3(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: bb4(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: bb5(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp1(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp2(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp3a(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp3(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp4(:,:,:)
    REAL, INTENT(OUT), ALLOCATABLE :: sp5(:,:,:)
    LOGICAL, INTENT(IN) :: digitize

    ! Local variables
    INTEGER :: STAT
    INTEGER :: I,J
    REAL, ALLOCATABLE :: random_array(:)
    REAL, ALLOCATABLE :: out_real(:)

    !
    ! Allocate random number array
    !
    ALLOCATE(random_array(noMC),&
         out_real(noMC),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate MC random/cal',&
            'Get_MC_Data','monte_carlo.f90')
    ENDIF

    !
    ! Allocate output arrays
    !
    ALLOCATE(Counts1(nelem,ndata,noMC),&
         Counts2(nelem,ndata,noMC),&
         Counts3(nelem,ndata,noMC),&
         Counts4(nelem,ndata,noMC),&
         Counts5(nelem,ndata,noMC),&
         bb3(ncal,ndata,noMC),&
         bb4(ncal,ndata,noMC),&
         bb5(ncal,ndata,noMC),&
         sp1(ncal,ndata,noMC),&
         sp2(ncal,ndata,noMC),&
         sp3a(ncal,ndata,noMC),&
         sp3(ncal,ndata,noMC),&
         sp4(ncal,ndata,noMC),&
         sp5(ncal,ndata,noMC),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate MC output arrays',&
            'Get_MC_Data','monte_carlo.f90')
    ENDIF
   
    IF( 0 .eq. MOD(noMC,2) )THEN
       DO I=1,ndata
          IF( 0 .eq. MOD(I,1000) )THEN
             WRITE(*,'('' MC EVEN Scanline : '',i7)')I
          ENDIF
          !
          ! Earth counts arrays
          !
          IF( ch3a_there(I) .eq. 1 )THEN
             DO J=1,nelem
                CALL Get_Random_EVEN(noMC,InCounts1(J,I),noise(1),&
                     out_real,random_array,digitize)
                Counts1(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts2(J,I),noise(2),&
                     out_real,random_array,digitize)
                Counts2(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts3(J,I),noise(3),&
                     out_real,random_array,digitize)
                Counts3(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts4(J,I),noise(5),&
                     out_real,random_array,digitize)
                Counts4(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts5(J,I),noise(6),&
                     out_real,random_array,digitize)
                Counts5(J,I,:) = out_real
             END DO
          ELSE
             DO J=1,nelem
                CALL Get_Random_EVEN(noMC,InCounts1(J,I),noise(1),&
                     out_real,random_array,digitize)
                Counts1(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts2(J,I),noise(2),&
                     out_real,random_array,digitize)
                Counts2(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts3(J,I),noise(4),&
                     out_real,random_array,digitize)
                Counts3(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts4(J,I),noise(5),&
                     out_real,random_array,digitize)
                Counts4(J,I,:) = out_real
                CALL Get_Random_EVEN(noMC,InCounts5(J,I),noise(6),&
                     out_real,random_array,digitize)
                Counts5(J,I,:) = out_real
             END DO
          ENDIF
          !
          ! ICT counts
          ! Deal with smoothing window case
          !
          DO J=1,ncal
             IF( InBB3(I) .gt. 0 )THEN
                CALL Get_Random_EVEN(noMC,Inbb3(I),noise(4),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InbbFilter3(J,I),noise(4),&
                     out_real,random_array,digitize)
             ENDIF
             bb3(J,I,:) = out_real
             IF( InBB4(I) .gt. 0 )THEN
                CALL Get_Random_EVEN(noMC,Inbb4(I),noise(5),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InbbFilter4(J,I),noise(5),&
                     out_real,random_array,digitize)
             ENDIF
             bb4(J,I,:) = out_real
             IF( InBB5(I) .gt. 0 )THEN
                CALL Get_Random_EVEN(noMC,Inbb5(I),noise(6),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InbbFilter5(J,I),noise(6),&
                     out_real,random_array,digitize)
             ENDIF
             bb5(J,I,:) = out_real
          END DO

          !
          ! Space counts
          !
          DO J=1,ncal
             CALL Get_Random_EVEN(noMC,Inspace1(J,I),noise(1),&
                  out_real,random_array,digitize)
             sp1(J,I,:) = out_real
             CALL Get_Random_EVEN(noMC,Inspace2(J,I),noise(2),&
                  out_real,random_array,digitize)
             sp2(J,I,:) = out_real
             CALL Get_Random_EVEN(noMC,Inspace3a(J,I),noise(3),&
                  out_real,random_array,digitize)
             sp3a(J,I,:) = out_real
             IF( 0 .lt. Inspace3(I) )THEN
                CALL Get_Random_EVEN(noMC,Inspace3(I),noise(4),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InspaceFilter3(J,I),noise(4),&
                     out_real,random_array,digitize)
             ENDIF
             sp3(J,I,:) = out_real
             IF( 0 .lt. Inspace4(I) )THEN
                CALL Get_Random_EVEN(noMC,Inspace4(I),noise(5),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InspaceFilter4(J,I),noise(5),&
                     out_real,random_array,digitize)
             ENDIF
             sp4(J,I,:) = out_real
             IF( 0 .lt. Inspace5(I) )THEN
                CALL Get_Random_EVEN(noMC,Inspace5(I),noise(6),&
                     out_real,random_array,digitize)
             ELSE
                CALL Get_Random_EVEN(noMC,InspaceFilter5(J,I),noise(6),&
                     out_real,random_array,digitize)
             ENDIF
             sp5(J,I,:) = out_real
          END DO
       END DO
    ELSE
       DO I=1,ndata
          IF( 0 .eq. MOD(I,1000) )THEN
             WRITE(*,'('' MC EVEN Scanline : '',i7)')I
          ENDIF
          !
          ! Earth counts arrays
          !
          IF( ch3a_there(I) .eq. 1 )THEN
             DO J=1,nelem
                CALL Get_Random_ODD(noMC,InCounts1(J,I),noise(1),&
                     out_real,random_array,digitize)
                Counts1(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts2(J,I),noise(2),&
                     out_real,random_array,digitize)
                Counts2(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts3(J,I),noise(3),&
                     out_real,random_array,digitize)
                Counts3(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts4(J,I),noise(5),&
                     out_real,random_array,digitize)
                Counts4(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts5(J,I),noise(6),&
                     out_real,random_array,digitize)
                Counts5(J,I,:) = out_real
             END DO
          ELSE
             DO J=1,nelem
                CALL Get_Random_ODD(noMC,InCounts1(J,I),noise(1),&
                     out_real,random_array,digitize)
                Counts1(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts2(J,I),noise(2),&
                     out_real,random_array,digitize)
                Counts2(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts3(J,I),noise(4),&
                     out_real,random_array,digitize)
                Counts3(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts4(J,I),noise(5),&
                     out_real,random_array,digitize)
                Counts4(J,I,:) = out_real
                CALL Get_Random_ODD(noMC,InCounts5(J,I),noise(6),&
                     out_real,random_array,digitize)
                Counts5(J,I,:) = out_real
             END DO
          ENDIF
          !
          ! ICT counts
          !
          DO J=1,ncal
             CALL Get_Random_ODD(noMC,Inbb3(I),noise(4),&
                  out_real,random_array,digitize)
             bb3(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inbb4(I),noise(5),&
                  out_real,random_array,digitize)
             bb4(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inbb5(I),noise(6),&
                  out_real,random_array,digitize)
             bb5(J,I,:) = out_real
          END DO
          !
          ! Space counts
          !
          DO J=1,ncal
             CALL Get_Random_ODD(noMC,Inspace1(J,I),noise(1),&
                  out_real,random_array,digitize)
             sp1(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inspace2(J,I),noise(2),&
                  out_real,random_array,digitize)
             sp2(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inspace3a(J,I),noise(3),&
                  out_real,random_array,digitize)
             sp3a(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inspace3(I),noise(4),&
                  out_real,random_array,digitize)
             sp3(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inspace4(I),noise(5),&
                  out_real,random_array,digitize)
             sp4(J,I,:) = out_real
             CALL Get_Random_ODD(noMC,Inspace5(I),noise(6),&
                  out_real,random_array,digitize)
             sp5(J,I,:) = out_real
          END DO
       END DO
    ENDIF

    !
    ! Deallocate
    !
    DEALLOCATE(random_array,out_real)

  END SUBROUTINE Get_MC_Data

  SUBROUTINE Make_Cal_Data(noMC,ncal,cal_counts,output)

    INTEGER, INTENT(IN) :: noMC
    INTEGER, INTENT(IN) :: ncal
    REAL, INTENT(IN) :: cal_counts(ncal,noMC)
    REAL, INTENT(INOUT) :: output(noMC)

    ! Local variables
    INTEGER :: J
    INTEGER :: K
    INTEGER :: ndata

    DO J=1,noMC
       output(J) = 0.
       ndata = 0
       DO K=1,ncal
          IF( cal_counts(K,J) .gt. 0 )THEN
             output(J) = output(J) + cal_counts(K,J)
             ndata = ndata + 1
          ENDIF
       END DO
       IF( 0 .eq. ndata )THEN
          output(J) = NAN_R
       ELSE
          output(J) = output(J) / ndata
       ENDIF
    END DO

  END SUBROUTINE Make_Cal_Data
    
  SUBROUTINE Get_Real_Random_EVEN(noMC,InCounts,noise,Counts,random_array,&
       digitize)

    INTEGER, INTENT(IN) :: noMC
    REAL, INTENT(IN) :: InCounts
    REAL, INTENT(IN) :: noise
    REAL, INTENT(OUT) :: Counts(noMC)
    REAL, INTENT(INOUT) :: random_array(noMC)
    LOGICAL, INTENT(IN) :: digitize

    ! Local variables
    INTEGER :: K

    IF( InCounts .lt. 0 )THEN
       Counts = NAN_R
       RETURN
    ENDIF
    CALL Get_Gaussian_Eq_Area_Even_Switch(noMC,random_array,&
         shuffle=.TRUE.)
    IF( digitize )THEN
       DO K=1,noMC
          Counts(K) = INT(InCounts + noise*random_array(K))
       END DO
    ELSE
       DO K=1,noMC
          Counts(K) = InCounts + noise*random_array(K)
       END DO
    ENDIF

  END SUBROUTINE Get_Real_Random_EVEN

  SUBROUTINE Get_Integer_Random_EVEN(noMC,InCounts,noise,Counts,random_array,&
       digitize)

    INTEGER, INTENT(IN) :: noMC
    INTEGER, INTENT(IN) :: InCounts
    REAL, INTENT(IN) :: noise
    REAL, INTENT(OUT) :: Counts(noMC)
    REAL, INTENT(INOUT) :: random_array(noMC)
    LOGICAL, INTENT(IN) :: digitize

    ! Local variables
    INTEGER :: K

    IF( InCounts .lt. 0 )THEN
       Counts = NAN_R
       RETURN
    ENDIF
    CALL Get_Gaussian_Eq_Area_Even_Switch(noMC,random_array,&
         shuffle=.TRUE.)
    IF( digitize )THEN
       DO K=1,noMC
          Counts(K) = INT(InCounts + noise*random_array(K))
       END DO
    ELSE
       DO K=1,noMC
          Counts(K) = InCounts + noise*random_array(K)
       END DO
    ENDIF

  END SUBROUTINE Get_Integer_Random_EVEN

  SUBROUTINE Get_Real_Random_ODD(noMC,InCounts,noise,Counts,random_array,&
       digitize)

    INTEGER, INTENT(IN) :: noMC
    REAL, INTENT(IN) :: InCounts
    REAL, INTENT(IN) :: noise
    REAL, INTENT(OUT) :: Counts(noMC)
    REAL, INTENT(INOUT) :: random_array(noMC)
    LOGICAL, INTENT(IN) :: digitize

    ! Local variables
    INTEGER :: K

    IF( InCounts .lt. 0 )THEN
       Counts = NAN_R
       RETURN
    ENDIF
    CALL Get_Gaussian_Eq_Area_ODD_Switch(noMC,random_array,&
         shuffle=.TRUE.)
    IF( digitize )THEN
       DO K=1,noMC
          Counts(K) = INT(InCounts + noise*random_array(K))
       END DO
    ELSE
       DO K=1,noMC
          Counts(K) = InCounts + noise*random_array(K)
       END DO
    ENDIF

  END SUBROUTINE Get_Real_Random_ODD

  SUBROUTINE Get_Integer_Random_ODD(noMC,InCounts,noise,Counts,random_array,&
       digitize)

    INTEGER, INTENT(IN) :: noMC
    INTEGER, INTENT(IN) :: InCounts
    REAL, INTENT(IN) :: noise
    REAL, INTENT(OUT) :: Counts(noMC)
    REAL, INTENT(INOUT) :: random_array(noMC)
    LOGICAL, INTENT(IN) :: digitize

    ! Local variables
    INTEGER :: K

    IF( InCounts .lt. 0 )THEN
       Counts = NAN_R
       RETURN
    ENDIF
    CALL Get_Gaussian_Eq_Area_ODD_Switch(noMC,random_array,&
         shuffle=.TRUE.)
    IF( digitize )THEN
       DO K=1,noMC
          Counts(K) = INT(InCounts + noise*random_array(K))
       END DO
    ELSE
       DO K=1,noMC
          Counts(K) = InCounts + noise*random_array(K)
       END DO
    ENDIF

  END SUBROUTINE Get_Integer_Random_ODD

END MODULE Monte_Carlo
