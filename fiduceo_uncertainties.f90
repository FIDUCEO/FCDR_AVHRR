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
! MT: 11-11-2017: Define temp variables us1,us2,us3a to store structured 
!     uncertainties on the reflectance channels
! MT: 11-11-2017: fix problem of value not filling array with 0.03 for 
!     u_structured_Ch1
! MT: 11-11-2017: fix problem of value not filling array with 0.05 for 
!     u_structured_Ch2
! MT: 11-11-2017: fix problem of value not filling array with 0.05 for 
!     u_structured_Ch3a
! MT: 13-11-2017: allocated nsmoothBB3,4,5 and nsmoothSp3,4,5 to AVHRRout data 
!     structure in combine_orbits.f90 so that the calculations don't fail 
! MT: 19-12-2017: v0.3pre 
! MT: 09-03-2018: v0.5beta 
!
! JM: 14-05-2018: Added GBCS CCI L1C output routine
! JM: 12-07-2018: Added in FIDUCEO measurement equations
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
  USE GbcsTypes
  USE GbcsConstants
  USE GbcsErrorHandler
  USE NOAA_LoadAVHRRLevel1B
  USE NETCDF
  USE GbcsDateTime
  USE FIDUCEO_Calibration

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
     REAL, ALLOCATABLE :: dre_over_da03(:,:)
     REAL, ALLOCATABLE :: dre_over_da04(:,:)
     REAL, ALLOCATABLE :: dre_over_da05(:,:)
     REAL, ALLOCATABLE :: dre_over_da13(:,:)
     REAL, ALLOCATABLE :: dre_over_da14(:,:)
     REAL, ALLOCATABLE :: dre_over_da15(:,:)
     REAL, ALLOCATABLE :: dre_over_da23(:,:)
     REAL, ALLOCATABLE :: dre_over_da24(:,:)
     REAL, ALLOCATABLE :: dre_over_da25(:,:)
     REAL, ALLOCATABLE :: dre_over_da33(:,:)
     REAL, ALLOCATABLE :: dre_over_da34(:,:)
     REAL, ALLOCATABLE :: dre_over_da35(:,:)
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
     REAL, ALLOCATABLE :: dre_over_dcs1(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs2(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs3a(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs3(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs4(:,:)
     REAL, ALLOCATABLE :: dre_over_dcs5(:,:)
     REAL, ALLOCATABLE :: uc3(:,:)
     REAL, ALLOCATABLE :: uc4(:,:)
     REAL, ALLOCATABLE :: uc5(:,:)
     REAL, ALLOCATABLE :: ur3(:,:)
     REAL, ALLOCATABLE :: ur4(:,:)
     REAL, ALLOCATABLE :: ur5(:,:)
     REAL, ALLOCATABLE :: us1(:,:)
     REAL, ALLOCATABLE :: us2(:,:)
     REAL, ALLOCATABLE :: us3a(:,:)
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
     !
     ! extra CURUC variables
     REAL, ALLOCATABLE :: dBT_over_dT3(:,:)
     REAL, ALLOCATABLE :: dBT_over_dT4(:,:)
     REAL, ALLOCATABLE :: dBT_over_dT5(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcs3(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcs4(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcs5(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcict3(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcict4(:,:)
     REAL, ALLOCATABLE :: dBT_over_dcict5(:,:)
     REAL, ALLOCATABLE :: hu3(:,:)
     REAL, ALLOCATABLE :: hu4(:,:)
     REAL, ALLOCATABLE :: hu5(:,:)
     !
     ! Harmonisation covar
     !
     REAL, ALLOCATABLE :: Cmatrix3(:)
     REAL, ALLOCATABLE :: Omatrix3(:)
     REAL, ALLOCATABLE :: Cmatrix(:)
     REAL, ALLOCATABLE :: Omatrix(:)

     REAL :: nuc(3)
     REAL :: aval(3)
     REAL :: bval(3)
     
     INTEGER, ALLOCATABLE :: flag_no_detection(:,:)
     INTEGER(GbcsInt1), ALLOCATABLE :: quality_channel_bitmask(:,:)
     INTEGER(GbcsInt1), ALLOCATABLE :: quality_scanline_bitmask(:)
  END TYPE FIDUCEO_Data

  !
  ! Structure to hold radiance/bt deltas from MC runs
  !
  TYPE mc_delta_str
     INTEGER :: nscan
     INTEGER :: nelem
     INTEGER :: nmc
     REAL, ALLOCATABLE :: ch1(:,:,:)
     REAL, ALLOCATABLE :: ch2(:,:,:)
     REAL, ALLOCATABLE :: ch3a(:,:,:)
     REAL, ALLOCATABLE :: ch3(:,:,:)
     REAL, ALLOCATABLE :: ch4(:,:,:)
     REAL, ALLOCATABLE :: ch5(:,:,:)
  END type mc_delta_str

  ! Constants
  REAL, PARAMETER :: c1 = 1.1910427E-5
  REAL, PARAMETER :: c2 = 1.4387752
  REAL, PARAMETER :: eta_ict = 0.985140
  REAL, PARAMETER :: prt_accuracy=0.1, prt_noise=0.
  !
  ! This is where the FIDUCEO software version number is defined
  !
  ! MT: 19-12-2017: v0.3pre
  ! MT: 09-03-2018: v0.5beta
  ! JM: 12/04/2019: v0.2Bet (Beta)
  CHARACTER(LEN=6) :: software_version = '0.2Bet'

! MT: 11-11-2017: Define temp variables to store structured uncertainties on the reflectance channels
  REAL, ALLOCATABLE :: us1(:,:)
  REAL, ALLOCATABLE :: us2(:,:)
  REAL, ALLOCATABLE :: us3a(:,:)
 
  PRIVATE
  PUBLIC :: FIDUCEO_Data
  PUBLIC :: Add_FIDUCEO_Uncert
  PUBLIC :: mc_delta_str
  PUBLIC :: prt_accuracy

CONTAINS

  SUBROUTINE Add_FIDUCEO_Uncert(IMG,AVHRR,uuid_in,filename_nc,&
       use_iasi_calibration,gbcs_l1c_output,&
       gbcs_l1c_cal,use_walton,keep_temp,write_fcdr,monte_carlo,&
       delta_radiance,seedval,ocean_only)

    TYPE(Imagery), INTENT(IN) :: IMG
    TYPE(AVHRR_Data), INTENT(INOUT) :: AVHRR
    CHARACTER(LEN=*), INTENT(IN) :: uuid_in
    CHARACTER(LEN=*), INTENT(IN) :: filename_nc
    LOGICAL, OPTIONAL :: use_iasi_calibration
    LOGICAL, OPTIONAL :: gbcs_l1c_output
    LOGICAL, OPTIONAL :: gbcs_l1c_cal
    LOGICAL, OPTIONAL :: use_walton
    LOGICAL, OPTIONAL :: keep_temp
    LOGICAL, OPTIONAL :: write_fcdr
    LOGICAL, OPTIONAL :: monte_carlo
    TYPE(mc_delta_str), OPTIONAL, INTENT(INOUT) :: delta_radiance
    INTEGER, OPTIONAL :: seedval
    LOGICAL, OPTIONAL :: ocean_only

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
    LOGICAL :: l1c_output
    LOGICAL :: l1c_output_cal
    TYPE(FIDUCEO_Data) :: FCDR    
    LOGICAL :: usewalton
    LOGICAL :: keeptemp
    LOGICAL :: writefcdr
    REAL, ALLOCATABLE :: covar_tot(:,:,:)
    REAL, ALLOCATABLE :: covar3(:,:)
    REAL, ALLOCATABLE :: covar4(:,:)
    REAL, ALLOCATABLE :: covar5(:,:)
    INTEGER :: nparam
    INTEGER :: nparam3
    INTEGER :: stat
    LOGICAL :: montecarlo
    TYPE(mc_delta_str) :: delta_bts
    LOGICAL :: oceanonly

    IF( .not. AVHRR%valid_data_there )THEN
       RETURN
    ENDIF
    IF( PRESENT(use_iasi_calibration) )THEN
       use_iasi_cal = use_iasi_calibration
    ELSE
       use_iasi_cal = .FALSE.
    ENDIF
    IF( PRESENT(gbcs_l1c_output) )THEN
       l1c_output = gbcs_l1c_output
    ELSE
       l1c_output = .FALSE.
    ENDIF
    IF( PRESENT(gbcs_l1c_cal) )THEN
       l1c_output_cal = gbcs_l1c_cal
    ELSE
       l1c_output_cal = .FALSE.
    ENDIF
    IF( PRESENT(use_walton) )THEN
       usewalton = use_walton
    ELSE
       usewalton = .FALSE.
    ENDIF
    IF( PRESENT(keep_temp) )THEN
       keeptemp = keep_temp
    ELSE
       keeptemp = .FALSE.
    ENDIF
    IF( PRESENT(write_fcdr) )THEN
       writefcdr = write_fcdr
    ELSE
       writefcdr = .TRUE.
    ENDIF
    IF( PRESENT(monte_carlo) )THEN
       montecarlo = monte_carlo
    ELSE
       montecarlo = .FALSE.
    ENDIF
    IF( PRESENT(ocean_only) )THEN
       oceanonly = ocean_only
    ELSE
       oceanonly = .FALSE.
    ENDIF

    IF( montecarlo .and. .not. PRESENT(delta_radiance) .and. &
         .not. PRESENT(seedval) )THEN
       CALL Gbcs_Critical(.TRUE.,'delta_radiance needs to be set for montecarlo',&
            'Add_FIDUCEO_Uncert','fiduceo_uncertainties.f90')            
    ENDIF

    call Get_Calib_Coefficients_FIDUCEO( IMG, AVHRR%time(AVHRR%start_valid), &
         AVHRR%AVHRR_No, coefs1, coefs2, coefs3, ncoefs, coefs_frac1, &
         coefs_frac2, twelve_micron_there, covar=covar_tot, nparam=nparam, &
         nparam3=nparam3)
    ALLOCATE(covar3(nparam3,nparam3),covar4(nparam,nparam),&
         covar5(nparam,nparam),stat=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate covar3/4/5',&
            'Add_FIDUCEO_Uncert','fiduceo_uncertainties.f90')       
    ENDIF
    covar3 = covar_tot(1:nparam3,1:nparam3,1)
    covar4 = covar_tot(1:nparam,1:nparam,2)
    covar5 = covar_tot(1:nparam,1:nparam,3)

!    FCDR%nuc = (/coefs1(5,1),coefs2(5,1),coefs3(5,1)/)
!    FCDR%aval = (/coefs1(6,1),coefs2(6,1),coefs3(6,1)/)
!    FCDR%bval = (/coefs1(7,1),coefs2(7,1),coefs3(7,1)/)
    
    IF( usewalton )THEN
       coefs1(5,1) = AVHRR%nuc(1)
       coefs1(6,1) = AVHRR%aval(1)
       coefs1(7,1) = AVHRR%bval(1)
       coefs2(5,1) = AVHRR%nuc(2)
       coefs2(6,1) = AVHRR%aval(2)
       coefs2(7,1) = AVHRR%bval(2)
       coefs3(5,1) = AVHRR%nuc(3)
       coefs3(6,1) = AVHRR%aval(3)
       coefs3(7,1) = AVHRR%bval(3)
    ENDIF
    !
    ! make sure coefs are correct for the form of Rad->Temp and Temp-Rad
    !
    IF( coefs1(6,1) .gt. 0 )THEN
       coefs1(6,1) = -coefs1(6,1)/coefs1(7,1)
       coefs1(7,1) = 1./coefs1(7,1)
       coefs2(6,1) = -coefs2(6,1)/coefs2(7,1)
       coefs2(7,1) = 1./coefs2(7,1)
       coefs3(6,1) = -coefs3(6,1)/coefs3(7,1)
       coefs3(7,1) = 1./coefs3(7,1)
    ENDIF

    FCDR%nuc = AVHRR%nuc
    FCDR%aval = AVHRR%aval
    FCDR%bval = AVHRR%bval

    !
    ! Recalculate radiances using potentially more complex models than
    ! available in the Level 1B reader
    !
    

    !
    ! Do FIDUCEO Uncertainties
    !
    DO I=1,AVHRR%arraySize
       CALL Get_Sensitivities(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)
       CALL calculate_urict(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there)
       CALL radiance_uncertainties(I,AVHRR,coefs1,coefs2,coefs3,FCDR,&
            twelve_micron_there,nparam3,nparam,covar3,covar4,covar5)
    END DO
    !
    ! Get Quality flags
    !
    CALL Get_Quality_Flags(AVHRR,FCDR)

    !
    ! If monte-carlo data, then calculate delta BTs
    ! Forces use of same coefficients
    !
    IF( montecarlo )THEN
       CALL Get_Delta_BTs(AVHRR,delta_radiance,coefs1,coefs2,coefs3,&
            delta_bts,twelve_micron_there)
       CALL Deallocate_delta(delta_radiance)
    ENDIF

    !
    ! If GBCS output - write out
    !
    IF( l1c_output )THEN
       CALL Write_GBCS_L1C(AVHRR,FCDR,uuid_in,twelve_micron_there,&
            l1c_output_cal)
    ELSE
       !
       ! write to NetCDF
       !
       temp_file=TRIM(uuid_in)//'.nc'
       CALL Write_Temp_NETCDF(temp_file,AVHRR,FCDR,twelve_micron_there,&
            montecarlo,delta_bts,seedval)
!    CALL Rescale(AVHRR,FCDR)
       IF( writefcdr )THEN
          !
          ! Because of Gerrit's CURUC routines needing Python 3 but pyGAC
          ! needed Python 2 we have to run a script here that runs the 
          ! write easy python code after switching environments
          !
          IF( 'None' .eq. filename_nc )THEN
!             command_fcdr ='python2.7 write_easy_fcdr_from_netcdf.py '//TRIM(temp_file)
             IF( oceanonly )THEN
                command_fcdr = './write_easy_fcdr.sh '//TRIM(temp_file)//' --ocean'
             ELSE
                command_fcdr = './write_easy_fcdr.sh '//TRIM(temp_file)
             ENDIF
          ELSE
!             command_fcdr ='python2.7 write_easy_fcdr_from_netcdf.py '//TRIM(temp_file)//' '//TRIM(filename_nc)
             IF( oceanonly )THEN
                command_fcdr = './write_easy_fcdr.sh '//TRIM(temp_file)//' --output '//TRIM(filename_nc)//' --ocean'
             ELSE
                command_fcdr = './write_easy_fcdr.sh '//TRIM(temp_file)//' --output '//TRIM(filename_nc)
             ENDIF
          ENDIF
          call SYSTEM(TRIM(command_fcdr))
       ENDIF
       IF( .not. keeptemp )THEN
          command_fcdr = 'rm -f '//TRIM(temp_file) !MT: 05-11-2017: comment to keep temp netcdf files
          call SYSTEM(TRIM(command_fcdr))
       ENDIF
!       print*, "remplissage"
!   Which is French for "filling"
!    call fill_netcdf(filename_nc,AVHRR,FCDR)
    ENDIF

    DEALLOCATE(covar3,covar4,covar5)

  END SUBROUTINE Add_FIDUCEO_Uncert

  !
  ! Deallocate delta arrays
  !
  SUBROUTINE Deallocate_delta(delta)

    TYPE(mc_delta_str), INTENT(INOUT) :: delta

    ! Local variables
    INTEGER :: STAT

    DEALLOCATE(delta%ch1,delta%ch2,delta%ch3a,delta%ch3,delta%ch4,delta%ch5,&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot DEALLOCATE delta',&
            'Deallocate_delta','fiduceo_uncertainties.f90')
    ENDIF

  END SUBROUTINE Deallocate_delta

  !
  ! Get deltaBTs from Monte-Carlo runs
  ! Note we have delta radiances so use the derivative of the Planck
  ! function to get delta B
  !
  SUBROUTINE Get_Delta_BTs(AVHRR,delta_radiance,coefs1,coefs2,coefs3,&
            delta_bts,twelve_micron_there)
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(mc_delta_str), INTENT(IN) :: delta_radiance
    REAL, INTENT(IN) :: coefs1(8,2)
    REAL, INTENT(IN) :: coefs2(8,2)
    REAL, INTENT(IN) :: coefs3(8,2)
    TYPE(mc_delta_str), INTENT(OUT) :: delta_bts
    LOGICAL, INTENT(IN) :: twelve_micron_there

    ! Local variables
    INTEGER :: I,J,K
    INTEGER :: STAT
    REAL :: nrmin
    REAL :: nrmax
    REAL :: trmin
    REAL :: trmax
    REAL :: newBT37
    REAL :: newBT11
    REAL :: newBT12
    
    !
    ! Check AVHRR/delta_radiance inputs
    !
    IF( AVHRR%arraySize .ne. delta_radiance%nscan .or. &
         AVHRR%nelem .ne. delta_radiance%nelem )THEN
       print *,AVHRR%nelem,delta_radiance%nelem
       print *,AVHRR%arraySize,delta_radiance%nscan
       CALL Gbcs_Critical(.TRUE.,'Mismatch in size',&
            'Get_Delta_BTs','fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Allocate BT arrays
    !
    delta_bts%nelem = delta_radiance%nelem
    delta_bts%nscan = delta_radiance%nscan
    delta_bts%nmc = delta_radiance%nmc
    ALLOCATE(delta_bts%ch1(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         delta_bts%ch2(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         delta_bts%ch3a(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         delta_bts%ch3(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         delta_bts%ch4(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         delta_bts%ch5(delta_bts%nelem,delta_bts%nscan,delta_bts%nmc),&
         STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Allocation error delta_bts',&
            'Get_Delta_BTs','fiduceo_uncertainties.f90')
    ENDIF
    
    !
    ! Copy over vis chan deltas (needed in radiance space)
    !
    delta_bts%ch1 = delta_radiance%ch1
    delta_bts%ch2 = delta_radiance%ch2
    delta_bts%ch3a = delta_radiance%ch3a

    !
    ! Get delta BTs
    !
    DO J=1,delta_bts%nscan
       DO K=1,delta_bts%nelem
          !
          ! As this is constant over all MC runs calculate here
          !
          IF( AVHRR%new_array3B(k,j) .gt. 0 )THEN
             newBT37 = convertBT(AVHRR%new_array3B(k,j),DBLE(coefs1(5,1)), &
                  DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
          ENDIF
          IF( AVHRR%new_array4(k,j) .gt. 0 )THEN
             newBT11 = convertBT(AVHRR%new_array4(k,j),DBLE(coefs2(5,1)), &
                  DBLE(coefs2(6,1)), DBLE(coefs2(7,1)))
          ENDIF
          IF( AVHRR%new_array5(k,j) .gt. 0 )THEN
             newBT12 = convertBT(AVHRR%new_array5(k,j),DBLE(coefs3(5,1)), &
                  DBLE(coefs3(6,1)), DBLE(coefs3(7,1)))
          ENDIF
          !
          ! Loop round MC runs
          !
          DO I=1,delta_bts%nmc
             IF( -1e20 .lt. delta_radiance%ch3(k,j,i) )THEN
                nrmax=AVHRR%new_array3B(k,j)+delta_radiance%ch3(k,j,i)
                IF( nrmax .gt. 0 )THEN
                   delta_bts%ch3(k,j,i)=convertBT(nrmax,DBLE(coefs1(5,1)), &
                        DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))-&
                        newBT37
                ELSE
                   delta_bts%ch3(k,j,i) = NAN_R
                ENDIF
             ELSE
                delta_bts%ch3(k,j,i) = NAN_R
             ENDIF

             IF( -1e20 .lt. delta_radiance%ch4(k,j,i) )THEN
                nrmax=AVHRR%new_array4(k,j)+delta_radiance%ch4(k,j,i)
                IF( nrmax .gt. 0 )THEN
                   delta_bts%ch4(k,j,i)=convertBT(nrmax,DBLE(coefs2(5,1)), &
                        DBLE(coefs2(6,1)), DBLE(coefs2(7,1)))-&
                        newBT11
                ELSE
                   delta_bts%ch4(k,j,i) = NAN_R
                ENDIF
             ELSE
                delta_bts%ch4(k,j,i) = NAN_R
             ENDIF

             IF( twelve_micron_there )THEN
                IF( -1e20 .lt. delta_radiance%ch5(k,j,i) )THEN
                   nrmax=AVHRR%new_array5(k,j)+delta_radiance%ch5(k,j,i)
                   IF( nrmax .gt. 0 )THEN
                      delta_bts%ch5(k,j,i)=convertBT(nrmax,DBLE(coefs3(5,1)), &
                           DBLE(coefs3(6,1)), DBLE(coefs3(7,1)))-&
                           newBT12
                   ELSE
                      delta_bts%ch5(k,j,i) = NAN_R
                   ENDIF
                ELSE
                   delta_bts%ch5(k,j,i) = NAN_R
                ENDIF
             ELSE
                delta_bts%ch5(k,j,i) = NAN_R
             ENDIF
          END DO
       END DO
    END DO

  END SUBROUTINE Get_Delta_BTs
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
  ! Return a string for sensor
  !
  SUBROUTINE noaa_name(AVHRR,noaa_string,noaa)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    CHARACTER(LEN=*), INTENT(OUT) :: noaa_string
    LOGICAL, OPTIONAL :: noaa

    ! Local variables
    LOGICAL :: noaa_type

    IF( PRESENT(noaa) )THEN
       noaa_type = noaa
    ELSE
       noaa_type = .TRUE.
    ENDIF

    IF( noaa_type )THEN
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
       ELSE
          CALL Gbcs_Critical(.TRUE.,'Cannot match AVHRR_No','noaa_name',&
               'fiduceo_uncertainties.f90')
       ENDIF
    ELSE
       IF( AVHRR%AVHRR_No .eq. 1 )THEN
          noaa_string = 'AVHRRTN_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 6 )THEN
          noaa_string = 'AVHRR06_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 7 )THEN
          noaa_string = 'AVHRR07_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 8 )THEN
          noaa_string = 'AVHRR08_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 9 )THEN
          noaa_string = 'AVHRR09_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 10 )THEN
          noaa_string = 'AVHRR10_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 11 )THEN
          noaa_string = 'AVHRR11_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 12 )THEN
          noaa_string = 'AVHRR12_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 14 )THEN
          noaa_string = 'AVHRR14_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 15 )THEN
          noaa_string = 'AVHRR15_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 16 )THEN
          noaa_string = 'AVHRR16_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 17 )THEN
          noaa_string = 'AVHRR17_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 18 )THEN
          noaa_string = 'AVHRR18_G'
       ELSE IF( AVHRR%AVHRR_No .eq. 19 )THEN
          noaa_string = 'AVHRR19_G'
       ELSE IF( AVHRR%AVHRR_No .eq. -1 )THEN
          noaa_string = 'AVHRRMTA_G'
       ELSE IF( AVHRR%AVHRR_No .eq. -2 )THEN
          noaa_string = 'AVHRRMTB_G'
       ELSE IF( AVHRR%AVHRR_No .eq. -3 )THEN
          noaa_string = 'AVHRRMTC_G'
       ELSE
          CALL Gbcs_Critical(.TRUE.,'Cannot match AVHRR_No','noaa_name',&
               'fiduceo_uncertainties.f90')
       ENDIF       
    ENDIF

  END SUBROUTINE noaa_name

  !
  ! Write a tempory netcdf file for python to then convert
  ! Different from Marine's original version
  !
  SUBROUTINE Write_Temp_NETCDF(temp_file,AVHRR,FCDR,twelve_micron_there,&
       monte_carlo,delta_bts,seedval)

    CHARACTER(LEN=*), INTENT(IN) :: temp_file
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(FIDUCEO_Data), INTENT(IN) :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there
    LOGICAL, INTENT(IN) :: monte_carlo
    TYPE(mc_delta_str), INTENT(INOUT) :: delta_bts
    INTEGER, INTENT(IN) :: seedval

    ! Local variables
    INTEGER :: I
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
    INTEGER :: ch1_common_varid
    INTEGER :: ch2_common_varid
    INTEGER :: ch3a_common_varid
    INTEGER :: ch3b_common_varid
    INTEGER :: ch4_common_varid
    INTEGER :: ch5_common_varid
    INTEGER :: scan_qual_varid
    INTEGER :: chan_qual_varid
    INTEGER :: dBT3_over_dT_varid
    INTEGER :: dBT4_over_dT_varid
    INTEGER :: dBT5_over_dT_varid
    INTEGER :: dRe1_over_dCS_varid
    INTEGER :: dRe2_over_dCS_varid
    INTEGER :: dRe3a_over_dCS_varid
    INTEGER :: dBT3_over_dCS_varid
    INTEGER :: dBT4_over_dCS_varid
    INTEGER :: dBT5_over_dCS_varid
    INTEGER :: dBT3_over_dCICT_varid
    INTEGER :: dBT4_over_dCICT_varid
    INTEGER :: dBT5_over_dCICT_varid
    INTEGER :: cal_cnts_varid
    INTEGER :: earth_cnts_varid
    INTEGER :: nuc_varid
    INTEGER :: aval_varid
    INTEGER :: bval_varid
    INTEGER :: scanline_varid
    INTEGER :: oscanline_varid
    INTEGER :: smoothprt_varid
    INTEGER :: badNavigation_varid
    INTEGER :: badCalibration_varid
    INTEGER :: badTime_varid
    INTEGER :: missingLines_varid
    INTEGER :: solar3_varid
    INTEGER :: solar4_varid
    INTEGER :: solar5_varid
    INTEGER :: ch3a_there_varid
    INTEGER :: ch3b_harm_varid
    INTEGER :: ch4_harm_varid
    INTEGER :: ch5_harm_varid
    INTEGER :: ch1_MC_varid
    INTEGER :: ch2_MC_varid
    INTEGER :: ch3a_MC_varid
    INTEGER :: ch3_MC_varid
    INTEGER :: ch4_MC_varid
    INTEGER :: ch5_MC_varid
    INTEGER :: stat

    INTEGER :: dimid_nx
    INTEGER :: dimid_ny
    INTEGER :: dimid_ir
    INTEGER :: dimid_band
    INTEGER :: dimid_mc
    INTEGER :: dims1(1)
    INTEGER :: dims2(2)
    INTEGER :: dims3(3)
    INTEGER :: dims_band(1)
    INTEGER(GbcsInt1), ALLOCATABLE :: badNav(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: badCal(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: badTime(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: missingLines(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: solar3(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: solar4(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: solar5(:)

    REAL :: noise_cnts(6)
    REAL :: earth_noise_cnts(6)

    INTEGER :: compress_level=5
    CHARACTER(LEN=512) :: noaa_string
    
    stat = NF90_CREATE(temp_file,IOR(NF90_HDF5,NF90_CLOBBER),ncid)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nx',AVHRR%nelem,dimid_nx)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'ny',AVHRR%arraySize,dimid_ny)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nir',6,dimid_ir)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nband_coef',3,dimid_band)
    call check(stat)

    if( monte_carlo )THEN
       stat = NF90_DEF_DIM(ncid,'n_montecarlo',delta_bts%nmc,dimid_mc)
       call check(stat)
    ENDIF

    dims1(1) = dimid_ny
    dims2(1) = dimid_nx
    dims2(2) = dimid_ny
    
    stat = NF90_DEF_VAR(ncid,'latitude',NF90_FLOAT,dims2,latitude_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, latitude_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, latitude_varid, 0, NAN_R)
    call check(stat)

    !
    ! Write variables including define - ensures valid min/max held
    !
    stat = NF90_DEF_VAR(ncid,'longitude',NF90_FLOAT,dims2,longitude_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, longitude_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, longitude_varid, 0, NAN_R)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'time',NF90_DOUBLE,dims1,time_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, time_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, time_varid, 0, NAN_R)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'year',NF90_INT,dims1,year_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, year_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, year_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'month',NF90_INT,dims1,month_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, month_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, month_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'day',NF90_INT,dims1,day_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, day_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, day_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'hours',NF90_FLOAT,dims1,hours_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, hours_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, hours_varid, 0, -1)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'satza',NF90_FLOAT,dims2,satza_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, satza_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, satza_varid, 0, NAN_R)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'solza',NF90_FLOAT,dims2,solza_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, solza_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, solza_varid, 0, NAN_R)
    call check(stat)

    IF( ALLOCATED(AVHRR%relaz) )THEN
       stat = NF90_DEF_VAR(ncid,'relaz',NF90_FLOAT,dims2,relaz_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, relaz_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, relaz_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1',NF90_FLOAT,dims2,ch1_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2',NF90_FLOAT,dims2,ch2_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_varid, 0, NAN_R)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a',NF90_FLOAT,dims2,ch3a_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_varid, 0, NAN_R)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b',NF90_FLOAT,dims2,ch3b_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4',NF90_FLOAT,dims2,ch4_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_varid, 0, NAN_R)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5',NF90_FLOAT,dims2,ch5_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1_random',NF90_FLOAT,dims2,ch1_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_random_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2_random',NF90_FLOAT,dims2,ch2_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_random_varid, 0, NAN_R)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a_random',NF90_FLOAT,dims2,ch3a_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_random_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_random_varid, 0, NAN_R)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b_random',NF90_FLOAT,dims2,ch3b_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_random_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4_random',NF90_FLOAT,dims2,ch4_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_random_varid, 0, NAN_R)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5_random',NF90_FLOAT,dims2,ch5_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_random_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_random_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1_non_random',NF90_FLOAT,dims2,ch1_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_non_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_non_random_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2_non_random',NF90_FLOAT,dims2,ch2_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_non_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_non_random_varid, 0, NAN_R)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a_non_random',NF90_FLOAT,dims2,ch3a_non_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_non_random_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_non_random_varid, 0, NAN_R)
       call check(stat)
    ENDIF
    
    stat = NF90_DEF_VAR(ncid,'ch3b_non_random',NF90_FLOAT,dims2,ch3b_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_non_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_non_random_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4_non_random',NF90_FLOAT,dims2,ch4_non_random_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_non_random_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_non_random_varid, 0, NAN_R)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5_non_random',NF90_FLOAT,dims2,ch5_non_random_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_non_random_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_non_random_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch1_common',NF90_FLOAT,dims2,ch1_common_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_common_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch1_common_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch2_common',NF90_FLOAT,dims2,ch2_common_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_common_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch2_common_varid, 0, NAN_R)
    call check(stat)
    
    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'ch3a_common',NF90_FLOAT,dims2,ch3a_common_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_common_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_common_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    stat = NF90_DEF_VAR(ncid,'ch3b_common',NF90_FLOAT,dims2,ch3b_common_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_common_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_common_varid, 0, NAN_R)
    call check(stat)
    
    stat = NF90_DEF_VAR(ncid,'ch4_common',NF90_FLOAT,dims2,ch4_common_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_common_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_common_varid, 0, NAN_R)
    call check(stat)
    
    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'ch5_common',NF90_FLOAT,dims2,ch5_common_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_common_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_common_varid, 0, NAN_R)
       call check(stat)
    ENDIF

    dims1(1) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'quality_scanline_bitmask',NF90_UBYTE,dims1,&
         scan_qual_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, scan_qual_varid, 1, 1, compress_level)
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
    stat = NF90_DEF_VAR_DEFLATE(ncid, chan_qual_varid, 1, 1, compress_level)
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
    ! Now write fields needed for CURUC process
    !
    dims2(1) = dimid_nx
    dims2(2) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'dBT3_over_dT',NF90_FLOAT,dims2,dBT3_over_dT_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT3_over_dT_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT3_over_dT_varid, 0, NAN_R)
    call check(stat)    

    stat = NF90_DEF_VAR(ncid,'dBT4_over_dT',NF90_FLOAT,dims2,dBT4_over_dT_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT4_over_dT_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT4_over_dT_varid, 0, NAN_R)
    call check(stat)    

    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'dBT5_over_dT',NF90_FLOAT,dims2,dBT5_over_dT_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, dBT5_over_dT_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, dBT5_over_dT_varid, 0, NAN_R)
       call check(stat)    
    ENDIF

    stat = NF90_DEF_VAR(ncid,'dRe1_over_dCS',NF90_FLOAT,dims2,&
         dRe1_over_dCS_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dRe1_over_dCS_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dRe1_over_dCS_varid, 0, NAN_R)
    call check(stat)    

    stat = NF90_DEF_VAR(ncid,'dRe2_over_dCS',NF90_FLOAT,dims2,&
         dRe2_over_dCS_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dRe2_over_dCS_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dRe2_over_dCS_varid, 0, NAN_R)
    call check(stat)    

    IF( ALLOCATED(AVHRR%array3a) )THEN
       stat = NF90_DEF_VAR(ncid,'dRe3a_over_dCS',NF90_FLOAT,dims2,&
            dRe3a_over_dCS_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, dRe3a_over_dCS_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, dRe3a_over_dCS_varid, 0, NAN_R)
       call check(stat)    
    ENDIF

    stat = NF90_DEF_VAR(ncid,'dBT3_over_dCS',NF90_FLOAT,dims2,&
         dBT3_over_dCS_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT3_over_dCS_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT3_over_dCS_varid, 0, NAN_R)
    call check(stat)    

    stat = NF90_DEF_VAR(ncid,'dBT4_over_dCS',NF90_FLOAT,dims2,&
         dBT4_over_dCS_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT4_over_dCS_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT4_over_dCS_varid, 0, NAN_R)
    call check(stat)    

    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'dBT5_over_dCS',NF90_FLOAT,dims2,&
            dBT5_over_dCS_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, dBT5_over_dCS_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, dBT5_over_dCS_varid, 0, NAN_R)
       call check(stat)    
    ENDIF

    stat = NF90_DEF_VAR(ncid,'dBT3_over_dCICT',NF90_FLOAT,dims2,&
         dBT3_over_dCICT_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT3_over_dCICT_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT3_over_dCICT_varid, 0, NAN_R)
    call check(stat)    

    stat = NF90_DEF_VAR(ncid,'dBT4_over_dCICT',NF90_FLOAT,dims2,&
         dBT4_over_dCICT_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, dBT4_over_dCICT_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, dBT4_over_dCICT_varid, 0, NAN_R)
    call check(stat)    

    IF( twelve_micron_there )THEN
       stat = NF90_DEF_VAR(ncid,'dBT5_over_dCICT',NF90_FLOAT,dims2,&
            dBT5_over_dCICT_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, dBT5_over_dCICT_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, dBT5_over_dCICT_varid, 0, NAN_R)
       call check(stat)    
    ENDIF

    !
    ! ICT Temperature
    !
    dims1(1) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'smoothPRT',NF90_FLOAT,dims1,&
         smoothprt_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, smoothprt_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, smoothprt_varid, 0, NAN_R)
    call check(stat)    

    !
    ! AVHRR Flags for L1C
    !
    stat = NF90_DEF_VAR(ncid,'badNavigation',NF90_BYTE,dims1,badNavigation_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, badNavigation_varid, 1, 1, compress_level)
    call check(stat)

    stat = NF90_DEF_VAR(ncid,'badCalibration',NF90_BYTE,dims1,badCalibration_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, badCalibration_varid, 1, 1, compress_level)
    call check(stat)
 
    stat = NF90_DEF_VAR(ncid,'badTime',NF90_BYTE,dims1,badTime_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, badTime_varid, 1, 1, compress_level)
    call check(stat)
 
    stat = NF90_DEF_VAR(ncid,'missingLines',NF90_BYTE,dims1,missingLines_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, missingLines_varid, 1, 1, compress_level)
    call check(stat)
 
    stat = NF90_DEF_VAR(ncid,'solar_contam_3b',NF90_BYTE,dims1,solar3_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, solar3_varid, 1, 1, compress_level)
    call check(stat)
 
    stat = NF90_DEF_VAR(ncid,'solar_contam_4',NF90_BYTE,dims1,solar4_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, solar4_varid, 1, 1, compress_level)
    call check(stat)
 
    stat = NF90_DEF_VAR(ncid,'solar_contam_5',NF90_BYTE,dims1,solar5_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, solar5_varid, 1, 1, compress_level)
    call check(stat)

    !
    ! Write band coefficients
    !
    dims_band(1) = dimid_band
    stat = NF90_DEF_VAR(ncid,'nuc',NF90_FLOAT,dims_band,nuc_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, nuc_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, nuc_varid, 0, NAN_R)
    call check(stat)    
    
    stat = NF90_DEF_VAR(ncid,'aval',NF90_FLOAT,dims_band,aval_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, aval_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, aval_varid, 0, NAN_R)
    call check(stat)    
    
    stat = NF90_DEF_VAR(ncid,'bval',NF90_FLOAT,dims_band,bval_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, bval_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, bval_varid, 0, NAN_R)
    call check(stat)    
    
    dims1(1) = dimid_ir
    stat = NF90_DEF_VAR(ncid,'cal_cnts_noise',NF90_FLOAT,dims1,cal_cnts_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, cal_cnts_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, cal_cnts_varid, 0, NAN_R)
    call check(stat)    
    
    stat = NF90_DEF_VAR(ncid,'cnts_noise',NF90_FLOAT,dims1,earth_cnts_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, earth_cnts_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, earth_cnts_varid, 0, NAN_R)
    call check(stat)    
    
    dims1(1) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'ch3a_there',NF90_INT,dims1,ch3a_there_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_there_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3a_there_varid, 0, NAN_I)
    call check(stat)    
    
    stat = NF90_DEF_VAR(ncid,'scanline',NF90_INT,dims1,scanline_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, scanline_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, scanline_varid, 0, NAN_I)
    call check(stat)    
    
    dims1(1) = dimid_ny
    stat = NF90_DEF_VAR(ncid,'orig_scanline',NF90_INT,dims1,oscanline_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, oscanline_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, oscanline_varid, 0, NAN_I)
    call check(stat)    

    !
    ! For Harmonisation common uncertainty (needed for channel to channel
    ! correlation matrix) write out total Harmonisation uncertainty
    ! Will assume sensitivity of 1 for this term in CURUC calculation
    !
    stat = NF90_DEF_VAR(ncid,'ch3b_harm_uncertainty',NF90_FLOAT,&
         dims2,ch3b_harm_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch3b_harm_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch3b_harm_varid, 0, NAN_R)
    call check(stat)        
    
    stat = NF90_DEF_VAR(ncid,'ch4_harm_uncertainty',NF90_FLOAT,&
         dims2,ch4_harm_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_harm_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch4_harm_varid, 0, NAN_R)
    call check(stat)        
    
    stat = NF90_DEF_VAR(ncid,'ch5_harm_uncertainty',NF90_FLOAT,&
         dims2,ch5_harm_varid)
    call check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_harm_varid, 1, 1, compress_level)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, ch5_harm_varid, 0, NAN_R)
    call check(stat)        

    !
    ! If monte_carlo the add varids
    !
    IF( monte_carlo )THEN
       dims3(1) = dimid_nx
       dims3(2) = dimid_ny
       dims3(3) = dimid_mc
       stat = NF90_DEF_VAR(ncid,'ch1_MC',NF90_FLOAT,&
            dims3,ch1_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch1_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch1_MC_varid, 0, NAN_R)
       call check(stat)        
       stat = NF90_DEF_VAR(ncid,'ch2_MC',NF90_FLOAT,&
            dims3,ch2_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch2_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch2_MC_varid, 0, NAN_R)
       call check(stat)        
       stat = NF90_DEF_VAR(ncid,'ch3a_MC',NF90_FLOAT,&
            dims3,ch3a_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3a_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3a_MC_varid, 0, NAN_R)
       call check(stat)        
       stat = NF90_DEF_VAR(ncid,'ch3_MC',NF90_FLOAT,&
            dims3,ch3_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch3_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch3_MC_varid, 0, NAN_R)
       call check(stat)        
       stat = NF90_DEF_VAR(ncid,'ch4_MC',NF90_FLOAT,&
            dims3,ch4_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch4_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch4_MC_varid, 0, NAN_R)
       call check(stat)        
       stat = NF90_DEF_VAR(ncid,'ch5_MC',NF90_FLOAT,&
            dims3,ch5_MC_varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, ch5_MC_varid, 1, 1, compress_level)
       call check(stat)
       stat = NF90_DEF_VAR_FILL(ncid, ch5_MC_varid, 0, NAN_R)
       call check(stat)               
    ENDIF
    
    !
    ! Define global attributes
    !
    CALL NOAA_Name(AVHRR,noaa_string)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'noaa_string',TRIM(noaa_string))
    call check(stat)    

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'version',TRIM(software_version))
    call check(stat)    

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'spatial_correlation_scale',&
         NPIXEL_PRT_SMOOTH)
    call check(stat) 

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'ICT_Temperature_Uncertainty',&
         AVHRR%ict_plane_uncert)
    call check(stat) 

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'PRT_Uncertainty',&
         prt_accuracy)
    call check(stat) 

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'orbital_temperature',&
         AVHRR%orbital_temperature)
    call check(stat) 

    DO I=1,AVHRR%norig_l1b
       IF( 1 .eq. I )THEN
          noaa_string = TRIM(AVHRR%orig_l1b(I))
       ELSE
          noaa_string = TRIM(noaa_string)//','//TRIM(AVHRR%orig_l1b(I))
       ENDIF
    END DO
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'sources',TRIM(noaa_string))
    call check(stat) 

    IF( monte_carlo )THEN
       stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'montecarlo_seed',&
            seedval)
       call check(stat) 
    ENDIF
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
!       us1 = 0.03
       us1 = CSPP_Uncertainty(1)
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array1) )THEN
!       stat = NF90_PUT_VAR(ncid, ch1_non_random_varid, 0.03*AVHRR%new_array1_error)
       stat = NF90_PUT_VAR(ncid, ch1_non_random_varid, FCDR%us1)
       call check(stat)
       stat = NF90_PUT_VAR(ncid, ch1_common_varid, us1)
       call check(stat)
    ENDIF

!    WRITE(*,*)'Ch2 (Non-Rand) writing'
!MT: 11-11-2017: fix problem of value not filling array     
    ALLOCATE(us2(AVHRR%nelem,AVHRR%arraySize), STAT=STAT)
    us2 = -1e30
    WHERE(AVHRR%new_array2.gt.-1e20)
!       us2 = 0.05
       us2 = CSPP_Uncertainty(2)
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array2) )THEN
!       stat = NF90_PUT_VAR(ncid, ch2_non_random_varid, 0.05*AVHRR%new_array2_error)
       stat = NF90_PUT_VAR(ncid, ch2_non_random_varid, FCDR%us2)
       call check(stat)
       stat = NF90_PUT_VAR(ncid, ch2_common_varid, us2)
       call check(stat)
    ENDIF

!       WRITE(*,*)'Ch3a (Non-Rand) writing'
!MT: 11-11-2017: fix problem of value not filling array     
    ALLOCATE(us3a(AVHRR%nelem,AVHRR%arraySize), STAT=STAT)
    us3a = -1e30
    WHERE(AVHRR%new_array3a.gt.-1e20)
!       us3a = 0.05
       us3a = CSPP_Uncertainty(3)
    ENDWHERE
    IF( ALLOCATED(AVHRR%new_array3a) )THEN
!       stat = NF90_PUT_VAR(ncid, ch3a_non_random_varid, 0.05*AVHRR%new_array3A_error)
       stat = NF90_PUT_VAR(ncid, ch3a_non_random_varid, FCDR%us3a)
       call check(stat)
       stat = NF90_PUT_VAR(ncid, ch3a_common_varid, us3a)
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

    stat = NF90_PUT_VAR(ncid, ch3b_common_varid, FCDR%uc3)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, ch4_common_varid, FCDR%uc4)
    call check(stat)

    IF( twelve_micron_there )THEN
       stat = NF90_PUT_VAR(ncid, ch5_common_varid, FCDR%uc5)
       call check(stat)
    ENDIF

    ! CURUC variables
    stat = NF90_PUT_VAR(ncid, dBT3_over_dT_varid, FCDR%dBT_over_dT3)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, dBT4_over_dT_varid, FCDR%dBT_over_dT4)
    call check(stat)

    IF( twelve_micron_there )THEN
       stat = NF90_PUT_VAR(ncid, dBT5_over_dT_varid, FCDR%dBT_over_dT5)
       call check(stat)
    ENDIF

    stat = NF90_PUT_VAR(ncid, dRe1_over_dCS_varid, FCDR%dRe_over_dcs1)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, dRe2_over_dCS_varid, FCDR%dRe_over_dcs2)
    call check(stat)

    IF( ALLOCATED(AVHRR%new_array3a) )THEN
       stat = NF90_PUT_VAR(ncid, dRe3a_over_dCS_varid, FCDR%dRe_over_dcs3a)
       call check(stat)
    ENDIF

    stat = NF90_PUT_VAR(ncid, dBT3_over_dCS_varid, FCDR%dBT_over_dcs3)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, dBT4_over_dCS_varid, FCDR%dBT_over_dcs4)
    call check(stat)

    IF( twelve_micron_there )THEN
       stat = NF90_PUT_VAR(ncid, dBT5_over_dCS_varid, FCDR%dBT_over_dcs5)
       call check(stat)
    ENDIF

    stat = NF90_PUT_VAR(ncid, dBT3_over_dCICT_varid, FCDR%dBT_over_dcict3)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, dBT4_over_dCICT_varid, FCDR%dBT_over_dcict4)
    call check(stat)

    IF( twelve_micron_there )THEN
       stat = NF90_PUT_VAR(ncid, dBT5_over_dCICT_varid, FCDR%dBT_over_dcict5)
       call check(stat)
    ENDIF

    stat = NF90_PUT_VAR(ncid, smoothprt_varid, AVHRR%smoothPRT)
    call check(stat)

    noise_cnts = AVHRR%noise_cnts_cal(:,1)
    stat = NF90_PUT_VAR(ncid, cal_cnts_varid, noise_cnts)
    call check(stat)

    earth_noise_cnts = AVHRR%noise_cnts(:,1)
    stat = NF90_PUT_VAR(ncid, earth_cnts_varid, earth_noise_cnts)
    call check(stat)

    ALLOCATE(badNav(AVHRR%arraySize),&
         badCal(AVHRR%arraySize),&
         badTime(AVHRR%arraySize),&
         missingLines(AVHRR%arraySize),&
         solar3(AVHRR%arraySize),&
         solar4(AVHRR%arraySize),&
         solar5(AVHRR%arraySize),&
         STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate badNav etc',&
            'Write_TEMP_NetCDF','fiduceo_uncertainties.f90')
    ENDIF
    DO I=1,AVHRR%arraySize
       IF( AVHRR%badNavigation(I) )THEN
          badNav(I) = 1
       ELSE
          badNav(I) = 0
       ENDIF
       IF( AVHRR%badCalibration(I) )THEN
          badCal(I) = 1
       ELSE
          badCal(I) = 0
       ENDIF
       IF( AVHRR%badTime(I) )THEN
          badTime(I) = 1
       ELSE
          badTime(I) = 0
       ENDIF
       IF( AVHRR%missingLines(I) )THEN
          missingLines(I) = 1
       ELSE
          missingLines(I) = 0
       ENDIF
       IF( AVHRR%solar_contamination_3b(I) )THEN
          solar3(I) = 1
       ELSE
          solar3(I) = 0
       ENDIF
       IF( AVHRR%solar_contamination_4(I) )THEN
          solar4(I) = 1
       ELSE
          solar4(I) = 0
       ENDIF
       IF( AVHRR%solar_contamination_5(I) )THEN
          solar5(I) = 1
       ELSE
          solar5(I) = 0
       ENDIF
    END DO
    stat = NF90_PUT_VAR(ncid, badNavigation_varid, badNav)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, badNavigation_varid, badCal)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, badTime_varid, badTime)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, missingLines_varid, missingLines)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, solar3_varid, solar3)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, solar4_varid, solar4)
    call check(stat)
    stat = NF90_PUT_VAR(ncid, solar5_varid, solar5)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, nuc_varid, FCDR%nuc)
    call check(stat)
    
    stat = NF90_PUT_VAR(ncid, aval_varid, FCDR%aval)
    call check(stat)
    
    stat = NF90_PUT_VAR(ncid, bval_varid, FCDR%bval)
    call check(stat)
    
!    stat = NF90_PUT_VAR(ncid, flag_no_detection_varid, FCDR%flag_no_detection)
!    call check(stat)

    stat = NF90_PUT_VAR(ncid, scan_qual_varid, &
         FCDR%quality_scanline_bitmask)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, chan_qual_varid, &
         FCDR%quality_channel_bitmask)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, ch3a_there_varid, AVHRR%ch3a_there)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, scanline_varid, AVHRR%scanlinenumber)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, oscanline_varid, AVHRR%scnline_l1b)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, ch3b_harm_varid, FCDR%hu3)
    call check(stat)

    stat = NF90_PUT_VAR(ncid, ch4_harm_varid, FCDR%hu4)
    call check(stat)

    IF( twelve_micron_there )THEN
       stat = NF90_PUT_VAR(ncid, ch5_harm_varid, FCDR%hu5)
       call check(stat)
    ENDIF

    IF( monte_carlo )THEN
       stat = NF90_PUT_VAR(ncid, ch1_MC_varid, delta_bts%ch1)
       call check(stat)
       stat = NF90_PUT_VAR(ncid, ch2_MC_varid, delta_bts%ch2)
       call check(stat)
       IF( ALLOCATED(AVHRR%new_array3a) )THEN
          stat = NF90_PUT_VAR(ncid, ch3a_MC_varid, delta_bts%ch3a)
          call check(stat)
       ENDIF
       stat = NF90_PUT_VAR(ncid, ch3_MC_varid, delta_bts%ch3)
       call check(stat)
       stat = NF90_PUT_VAR(ncid, ch4_MC_varid, delta_bts%ch4)
       call check(stat)
       IF( twelve_micron_there )THEN
          stat = NF90_PUT_VAR(ncid, ch5_MC_varid, delta_bts%ch5)
          call check(stat)
       ENDIF
    ENDIF

    stat = NF90_CLOSE(ncid)
    call check(stat)

    DEALLOCATE(badNav,badCal,badTime,missingLines,solar3,solar4,solar5)
    IF( monte_carlo )THEN
       CALL Deallocate_delta(delta_bts)
    ENDIF

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
    REAL:: urict, tstar, &
         aval,bval,nuc, &
         uc1=0,uc2=0,unuc=0,uaval=0,ubval=0
    INTEGER :: STAT
    REAL :: prtmean
    REAL :: utict_r
    REAL :: utict_s
    REAL :: utict_2s
    REAL :: prtsigma
    REAL :: utstar3_r
    REAL :: utstar4_r
    REAL :: utstar5_r
    REAL :: utstar3_s
    REAL :: utstar4_s
    REAL :: utstar5_s
    REAL :: utstar3_2s
    REAL :: utstar4_2s
    REAL :: utstar5_2s
    
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
    if ((outData%smoothprt1(i) .ne. NAN_R ) &
         .and. (outData%smoothprt2(i) .ne. NAN_R )&
         .and. (outData%smoothprt3(i) .ne. NAN_R )&
         .and. (outData%smoothprt4(i) .ne. NAN_R )) then
       prtmean = outData%smoothprt(i)
       IF( outData%prt_correction(i) .gt. -1e20 )THEN
          prtmean = prtmean + outData%prt_correction(i)
       ENDIF
       !
       ! Use dispersion of 4 PRTs as estimate of representivity error
       !
       ! On short timescales will be systematic
       !
       ! Replace prtsigma with estimate from new ICT correction model
       !
!       prtsigma=sqrt(((outData%prt1(i)-prtmean)**2&
!            +(outData%prt2(i)-prtmean)**2&
!            +(outData%prt3(i)-prtmean)**2&
!            +(outData%prt4(i)-prtmean)**2)/3.)
       prtsigma = outData%ict_plane_uncert

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
       utict_2s = sqrt(prt_accuracy**2 + (2.*prtsigma)**2)

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
    if (prtmean .gt.  -1.E+20 ) then 
       tstar=(prtmean-aval)/bval

       FCDR%Rict_c3(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)

       utstar3_s=sqrt(uaval**2*FCDR%dtstar_over_daval3(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval3(i)**2 &
            +utict_s**2*FCDR%dtstar_over_dT3(i)**2)
       utstar3_2s=sqrt(uaval**2*FCDR%dtstar_over_daval3(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval3(i)**2 &
            +utict_2s**2*FCDR%dtstar_over_dT3(i)**2)
       utstar3_r=sqrt(utict_r**2*FCDR%dtstar_over_dT3(i)**2)
       !print*, outdata%utstar3(i)  
       FCDR%urict3_r(i)=sqrt(utstar3_r**2*FCDR%drict_over_dtstar3(i)**2)
       !
       ! If solar contamination double ICT temperature correction uncertainty
       !
       IF( outData%solar_contamination_3b(i) )THEN
          FCDR%urict3_s(i)=sqrt(utstar3_2s**2*FCDR%drict_over_dtstar3(i)**2 + &
               unuc**2*FCDR%drict_over_dnuc3(i)**2)
       ELSE
          FCDR%urict3_s(i)=sqrt(utstar3_s**2*FCDR%drict_over_dtstar3(i)**2 + &
               unuc**2*FCDR%drict_over_dnuc3(i)**2)
       ENDIF
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
    if (prtmean .gt. -1e20 )then
       tstar=(prtmean-aval)/bval

       FCDR%Rict_c4(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)

       utstar4_s=sqrt(uaval**2*FCDR%dtstar_over_daval4(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval4(i)**2 &
            +utict_s**2*FCDR%dtstar_over_dT4(i)**2)
       utstar4_2s=sqrt(uaval**2*FCDR%dtstar_over_daval4(i)**2 &
            +ubval**2*FCDR%dtstar_over_dbval4(i)**2 &
            +utict_2s**2*FCDR%dtstar_over_dT4(i)**2)
       utstar4_r=sqrt(utict_r**2*FCDR%dtstar_over_dT4(i)**2)
       !print*, outdata%utstar3(i)  
       FCDR%urict4_r(i)=sqrt(utstar4_r**2*FCDR%drict_over_dtstar4(i)**2)
       IF( outData%solar_contamination_3b(i) )THEN
          FCDR%urict4_s(i)=sqrt(utstar4_2s**2*FCDR%drict_over_dtstar4(i)**2 + &
               unuc**2*FCDR%drict_over_dnuc4(i)**2)
       ELSE
          FCDR%urict4_s(i)=sqrt(utstar4_s**2*FCDR%drict_over_dtstar4(i)**2 + &
               unuc**2*FCDR%drict_over_dnuc4(i)**2)
       ENDIF
       !print*, outdata%rict_c4(i), utict,outdata%urict4(i)
    end if
    IF( twelve_micron_there )THEN
       !Ch5 
       nuc=coefs3(5,1)
       aval=coefs3(6,1)
       bval=coefs3(7,1)
       unuc=0
       if (prtmean .gt. -1e30 )then
          tstar=(prtmean-aval)/bval
          FCDR%Rict_c5(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)
          !print*, tstar, outdata%rict_c5(i)
          
          utstar5_s=sqrt(uaval**2*FCDR%dtstar_over_daval5(i)**2 &
               +ubval**2*FCDR%dtstar_over_dbval5(i)**2 &
               +utict_s**2*FCDR%dtstar_over_dT5(i)**2)
          utstar5_2s=sqrt(uaval**2*FCDR%dtstar_over_daval5(i)**2 &
               +ubval**2*FCDR%dtstar_over_dbval5(i)**2 &
               +utict_2s**2*FCDR%dtstar_over_dT5(i)**2)
          utstar5_r=sqrt(utict_r**2*FCDR%dtstar_over_dT5(i)**2)
          !print*, outdata%utstar3(i)  
          FCDR%urict5_r(i)=sqrt(utstar5_r**2*FCDR%drict_over_dtstar5(i)**2)
          IF( outData%solar_contamination_3b(i) )THEN
             FCDR%urict5_s(i)=sqrt(utstar5_2s**2*FCDR%drict_over_dtstar5(i)**2 + &
                  unuc**2*FCDR%drict_over_dnuc5(i)**2)
          ELSE
             FCDR%urict5_s(i)=sqrt(utstar5_s**2*FCDR%drict_over_dtstar5(i)**2 + &
                  unuc**2*FCDR%drict_over_dnuc5(i)**2)
          ENDIF
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
    REAL :: factor
    REAL :: dBT_over_dre3
    REAL :: dBT_over_dre4
    REAL :: dBT_over_dre5
    REAL :: smoothprt

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
            FCDR%dre_over_da03(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da04(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da05(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da13(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da14(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da15(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da23(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da24(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da25(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da33(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da34(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_da35(outdata%nelem,outdata%arraySize),&
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
            FCDR%dre_over_dcs1(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs2(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs3a(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs3(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs4(outdata%nelem,outdata%arraySize),&
            FCDR%dre_over_dcs5(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dT3(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dT4(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dT5(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcs3(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcs4(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcs5(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcict3(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcict4(outdata%nelem,outdata%arraySize),&
            FCDR%dBT_over_dcict5(outdata%nelem,outdata%arraySize),&
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
       FCDR%dre_over_da03 = NAN_R
       FCDR%dre_over_da04 = NAN_R
       FCDR%dre_over_da05 = NAN_R
       FCDR%dre_over_da13 = NAN_R
       FCDR%dre_over_da14 = NAN_R
       FCDR%dre_over_da15 = NAN_R
       FCDR%dre_over_da23 = NAN_R
       FCDR%dre_over_da24 = NAN_R
       FCDR%dre_over_da25 = NAN_R
       FCDR%dre_over_da33 = NAN_R
       FCDR%dre_over_da34 = NAN_R
       FCDR%dre_over_da35 = NAN_R
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
       FCDR%dre_over_dcs1 = NAN_R
       FCDR%dre_over_dcs2 = NAN_R
       FCDR%dre_over_dcs3a = NAN_R
       FCDR%dre_over_dcs3 = NAN_R
       FCDR%dre_over_dcs4 = NAN_R
       FCDR%dre_over_dcs5 = NAN_R
       FCDR%dBT_over_dT3 = NAN_R
       FCDR%dBT_over_dT4 = NAN_R
       FCDR%dBT_over_dT5 = NAN_R
       FCDR%dBT_over_dcs3 = NAN_R
       FCDR%dBT_over_dcs4 = NAN_R
       FCDR%dBT_over_dcs5 = NAN_R
       FCDR%dBT_over_dcict3 = NAN_R
       FCDR%dBT_over_dcict4 = NAN_R
       FCDR%dBT_over_dcict5 = NAN_R
    ENDIF

    !
    ! Get smooth PRT corrected if necessary
    !
    smoothprt = NAN_R
    IF( outData%smoothPRT(i) .gt. 0. )THEN
       IF( outData%prt_correction(i) .gt. -1e20 )THEN
          smoothprt = outData%smoothPRT(I) + outData%prt_correction(i)
       ELSE
          smoothprt = outData%smoothPRT(I)
       ENDIF
    ENDIF
    !
    ! Calculate Rict radiances which was originally in Marines 
    ! average_per_orbite code
    !
    IF( smoothprt .gt. -1e20 )THEN
       FCDR%Rict_c3(i) = convertRadiance(smoothPRT, DBLE(coefs1(5,1)), &
            DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
       FCDR%Rict_c4(i) = convertRadiance(smoothPRT, DBLE(coefs2(5,1)), &
            DBLE(coefs2(6,1)), DBLE(coefs2(7,1)))
       IF( twelve_micron_there )THEN
          FCDR%Rict_c5(i) = convertRadiance(smoothPRT, DBLE(coefs3(5,1)), &
               DBLE(coefs3(6,1)), DBLE(coefs3(7,1)))
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

    if (smoothprt .gt. 0)then  
       FCDR%dtstar_over_dbval3(i)=-(smoothprt-coefs1(6,1))/coefs1(7,1**2)
       FCDR%dtstar_over_dbval4(i)=-(smoothprt-coefs2(6,1))/coefs2(7,1**2)
       IF( twelve_micron_there )THEN
          FCDR%dtstar_over_dbval5(i)=-(smoothprt-coefs3(6,1))/coefs3(7,1**2)

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
    if (smoothprt .gt. 0) then 
       tstar=(smoothprt-coefs1(6,1))/coefs1(7,1)
       u=c2*coefs1(5,1)/tstar
       FCDR%drict_over_dnuc3(i)=(c1*coefs1(5,1)**2/(exp(u)-1)) &
            *(3-coefs1(5,1)*(c2/tstar)*exp(u)/(exp(u)-1))
       FCDR%drict_over_dtstar3(i)=c1*coefs1(5,1)**3*&
            (c2*coefs1(5,1)/tstar**2)&
            *exp(u)/(exp(u)-1)**2

       tstar=(smoothprt-coefs2(6,1))/coefs2(7,1)
       u=c2*coefs2(5,1)/tstar
       FCDR%drict_over_dnuc4(i)=(c1*coefs2(5,1)**2/(exp(u)-1)) &
            *(3-coefs2(5,1)*(c2/tstar)*exp(u)/(exp(u)-1))
       FCDR%drict_over_dtstar4(i)=c1*coefs2(5,1)**3*&
            (c2*coefs2(5,1)/tstar**2)&
            *exp(u)/(exp(u)-1)**2

       IF( twelve_micron_there )THEN
          tstar=(smoothprt-coefs3(6,1))/coefs3(7,1)
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
            (outdata%Counts3(j,i) .ge. 0) .and. &
            (outdata%orbital_temperature .ne. NAN_R) .and. &
            (outdata%orbital_temperature .gt. 0) .and. &
            (cs_cict_3 .ne. NAN_R) .and. (cs_cict_3 .gt. 0)) then 

          FCDR%dre_over_da03(j,i)=1

          FCDR%dre_over_da13(j,i)=FCDR%Rict_c3(i)*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3

          FCDR%dre_over_da23(j,i)=-cs_cict_3*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               +(outdata%smoothsp3(i)-outData%Counts3(j,i))**2

          FCDR%dre_over_da33(j,i)=FIDUCEO_Tinstr_Model(1,outdata%orbital_temperature)
          FCDR%dre_over_dtinstr3(j,i)=coefs1(8,1)

          FCDR%dre_over_dce3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               /cs_cict_3 &
               -2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))

          FCDR%dre_over_dcict3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3**2&
               +2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))

          FCDR%dre_over_dcs1(:,i) = outdata%new_array1_dsp(i)
          FCDR%dre_over_dcs2(:,i) = outdata%new_array2_dsp(i)
          FCDR%dre_over_dcs3a(:,i) = outdata%new_array3a_dsp(i)

          FCDR%dre_over_dcs3(j,i)=-((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3**2&
               +((eta_ict+coefs1(2,1))*FCDR%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
               /cs_cict_3

          FCDR%dre_over_drict3(j,i)=(eta_ict+coefs1(2,1))*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
               /cs_cict_3
          
          IF( 0 .le. outData%new_array3b(j,i) )THEN
             factor = C1*(FCDR%nuc(1)**3)/outData%new_array3b(j,i)+1.
             dBT_over_dre3 = C1*C2*FCDR%nuc(1)**4/log(factor)**2/factor/&
                  FCDR%bval(1)/(outData%new_array3b(j,i)**2)

             FCDR%dBT_over_dT3(j,i) = dBT_over_dre3*FCDR%dre_over_drict3(j,i)*&
                  FCDR%drict_over_dtstar3(i)*FCDR%dtstar_over_dT3(i)

             FCDR%dBT_over_dcs3(j,i) = dBT_over_dre3*FCDR%dre_over_dcs3(j,i)
             FCDR%dBT_over_dcict3(j,i) = dBT_over_dre3*FCDR%dre_over_dcict3(j,i)
          ELSE
             FCDR%dBT_over_dT3(j,i) = NAN_R
             FCDR%dBT_over_dcs3(j,i) = NAN_R
             FCDR%dBT_over_dcict3(j,i) = NAN_R
          ENDIF
       end if

       !---Ch4
       !----Sensitivities  of earth radiance
       if ((outdata%smoothsp4(i) .ne. NAN_R) .and. &
            (outdata%smoothsp4(i) .gt. 0) .and. &
            (FCDR%rict_c4(i) .ne. NAN_R) .and. &
            (FCDR%rict_c4(i) .gt. 0) .and. &
            (outData%Counts4(j,i).ne. NAN_R) .and. &
            (outdata%Counts4(j,i) .ge. 0) .and. &
            (outdata%orbital_temperature .ne. NAN_R) .and. &
            (outdata%orbital_temperature .gt. 0) .and. &
            (cs_cict_4 .ne. NAN_R) .and. (cs_cict_4 .gt. 0)) then 

          FCDR%dre_over_da04(j,i)=1

          FCDR%dre_over_da14(j,i)=FCDR%Rict_c4(i)*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               /cs_cict_4

          FCDR%dre_over_da24(j,i)=-cs_cict_4*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
               +(outdata%smoothsp4(i)-outData%Counts4(j,i))**2

          FCDR%dre_over_da34(j,i)=FIDUCEO_Tinstr_Model(2,outdata%orbital_temperature)
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

          IF( 0 .le. outData%new_array4(j,i) )THEN
             factor = C1*(FCDR%nuc(2)**3)/outData%new_array4(j,i)+1.
             dBT_over_dre4 = C1*C2*FCDR%nuc(2)**4/log(factor)**2/factor/&
                  FCDR%bval(2)/(outData%new_array4(j,i)**2)

             FCDR%dBT_over_dT4(j,i) = dBT_over_dre4*FCDR%dre_over_drict4(j,i)*&
                  FCDR%drict_over_dtstar4(i)*FCDR%dtstar_over_dT4(i)

             FCDR%dBT_over_dcs4(j,i) = dBT_over_dre4*FCDR%dre_over_dcs4(j,i)
             FCDR%dBT_over_dcict4(j,i) = dBT_over_dre4*FCDR%dre_over_dcict4(j,i)
          ELSE
             FCDR%dBT_over_dT4(j,i) = NAN_R
             FCDR%dBT_over_dcs4(j,i) = NAN_R
             FCDR%dBT_over_dcict4(j,i) = NAN_R
          ENDIF
       end if

       IF( twelve_micron_there )THEN

          !---Ch5
          !----Sensitivities  of earth radiance
          if ((outdata%smoothsp5(i) .ne. NAN_R) .and. &
               (outdata%smoothsp5(i) .gt. 0) .and. &
               (FCDR%rict_c5(i) .ne. NAN_R) .and. &
               (FCDR%rict_c5(i) .gt. 0) .and. &
               (outData%Counts5(j,i).ne. NAN_R) .and. &
               (outdata%Counts5(j,i) .ge. 0) .and. &
               (outdata%orbital_temperature .ne. NAN_R) .and. &
               (outdata%orbital_temperature .gt. 0) .and. &
               (cs_cict_5 .ne. NAN_R) .and. (cs_cict_5 .gt. 0)) then 
             
             FCDR%dre_over_da05(j,i)=1

             FCDR%dre_over_da15(j,i)=FCDR%Rict_c5(i)*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5

             FCDR%dre_over_da25(j,i)=-cs_cict_5*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  +(outdata%smoothsp5(i)-outData%Counts5(j,i))**2

             FCDR%dre_over_da35(j,i)=FIDUCEO_Tinstr_Model(3,outdata%orbital_temperature)
             FCDR%dre_over_dtinstr5(j,i)=coefs3(8,1)

             FCDR%dre_over_dce5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs3(4,1)*cs_cict_5**2)&
                  /cs_cict_5 &
                  -2*coefs3(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))

             FCDR%dre_over_dcict5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs3(4,1)*cs_cict_5**2) &
                  *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5**2&
                  +2*coefs3(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))
             
             FCDR%dre_over_dcs5(j,i)=-((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs3(4,1)*cs_cict_5**2)&
                  *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5**2&
                  +((eta_ict+coefs3(2,1))*FCDR%Rict_c5(i)-coefs3(4,1)*cs_cict_5**2)&
                  /cs_cict_5
             
             FCDR%dre_over_drict5(j,i)=(eta_ict+coefs3(2,1))*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                  /cs_cict_5

             IF( 0 .le. outData%new_array5(j,i) )THEN
                factor = C1*(FCDR%nuc(3)**3)/outData%new_array5(j,i)+1.
                dBT_over_dre5 = C1*C2*FCDR%nuc(3)**4/log(factor)**2/factor/&
                     FCDR%bval(3)/(outData%new_array5(j,i)**2)
                
                FCDR%dBT_over_dT5(j,i) = dBT_over_dre5*FCDR%dre_over_drict5(j,i)*&
                     FCDR%drict_over_dtstar5(i)*FCDR%dtstar_over_dT5(i)
                
                FCDR%dBT_over_dcs5(j,i) = dBT_over_dre5*FCDR%dre_over_dcs5(j,i)
                FCDR%dBT_over_dcict5(j,i) = dBT_over_dre5*FCDR%dre_over_dcict5(j,i)
             ELSE
                FCDR%dBT_over_dT5(j,i) = NAN_R
                FCDR%dBT_over_dcs5(j,i) = NAN_R
                FCDR%dBT_over_dcict5(j,i) = NAN_R
             ENDIF
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
       twelve_micron_there,nparams3,nparams,covar3,covar4,covar5)

    INTEGER, INTENT(IN) :: i
    TYPE(AVHRR_Data), INTENT(IN) :: outData
    REAL, INTENT(IN) :: coefs1(8,2)
    REAL, INTENT(IN) :: coefs2(8,2)
    REAL, INTENT(IN) :: coefs3(8,2)
    TYPE(FIDUCEO_Data), INTENT(INOUT) :: FCDR
    LOGICAL, INTENT(IN) :: twelve_micron_there
    INTEGER :: nparams3
    INTEGER :: nparams
    REAL :: covar3(1:nparams3,1:nparams3)
    REAL :: covar4(1:nparams,1:nparams)
    REAL :: covar5(1:nparams,1:nparams)

    ! Local variables
    INTEGER :: j
    INTEGER :: STAT
    REAL :: cs_cict_3
    REAL :: cs_cict_4
    REAL :: cs_cict_5
    REAL :: noise

    REAL :: nrmax5,ur5, nrmin5, trmin5, trmax5, &
         nsmax5,us5, nsmin5, tsmin5, tsmax5, &
         nrmax3,ur3, nrmin3, trmin3, trmax3, &
         nsmax3,us3, nsmin3, tsmin3, tsmax3, &
         nrmax4,ur4, nrmin4, trmin4, trmax4, &
         nsmax4,us4, nsmin4, tsmin4, tsmax4, &
         uc3, uc4, uc5
    REAL :: two_sigma

    !
    ! Allocate arrays
    !
    IF( .not. ALLOCATED(FCDR%uc3) )THEN
       ALLOCATE(FCDR%uc3(outData%nelem,outData%arraySize),&
            FCDR%uc4(outData%nelem,outData%arraySize),&
            FCDR%uc5(outData%nelem,outData%arraySize),&
            FCDR%ur3(outData%nelem,outData%arraySize),&
            FCDR%ur4(outData%nelem,outData%arraySize),&
            FCDR%ur5(outData%nelem,outData%arraySize),&
            FCDR%us1(outData%nelem,outData%arraySize),&
            FCDR%us2(outData%nelem,outData%arraySize),&
            FCDR%us3a(outData%nelem,outData%arraySize),&
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
            FCDR%Cmatrix3(nparams3),&
            FCDR%Omatrix3(nparams3),&
            FCDR%Cmatrix(nparams),&
            FCDR%Omatrix(nparams),&
            FCDR%hu3(outData%nelem,outData%arraySize),&
            FCDR%hu4(outData%nelem,outData%arraySize),&
            FCDR%hu5(outData%nelem,outData%arraySize),&
       STAT=STAT)
       IF( 0 .ne. STAT )THEN
          CALL Gbcs_Critical(.TRUE.,'Cannot allocate arrays',&
               'radiance_uncertainties','fiduceo_uncertainties.f90')
       ENDIF
       FCDR%uc3 = NAN_R
       FCDR%uc4 = NAN_R
       FCDR%uc5 = NAN_R
       FCDR%ur3 = NAN_R
       FCDR%ur4 = NAN_R
       FCDR%ur5 = NAN_R
       FCDR%us1 = NAN_R
       FCDR%us2 = NAN_R
       FCDR%us3a = NAN_R
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
       FCDR%hu3 = NAN_R
       FCDR%hu4 = NAN_R
       FCDR%hu5 = NAN_R
    ENDIF

    !
    ! Vis channel averaging noise. Inlcudes digitisation noise
    ! Note digitisation noise does not average down
    !
    IF( outData%new_array1_dsp(i) .ne. NAN_R )THEN
       FCDR%us1(:,i)=ABS(outData%new_array1_dsp(i)*outData%new_array1_spnoise(i))
    ENDIF
    IF( outData%new_array2_dsp(i) .ne. NAN_R )THEN
       FCDR%us2(:,i)=ABS(outData%new_array2_dsp(i)*outData%new_array2_spnoise(i))
    ENDIF
    IF( outData%new_array3a_dsp(i) .gt. -1e20 .and. &
         outData%new_array3a_spnoise(i) .gt. -1e20 )THEN
       FCDR%us3a(:,i)=ABS(outData%new_array3a_dsp(i)*outData%new_array3a_spnoise(i))
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

       FCDR%uc3(j,i)=NAN_R
       FCDR%uc4(j,i)=NAN_R
       FCDR%uc5(j,i)=NAN_R
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
       
!MT: 13-11-2017: allocated nsmoothBB3,4,5 and nsmoothSp3,4,5 to AVHRRout data 
!    structure in combine_orbits.f90 so that the calculations don't fail 
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
       uc3 = NAN_R
       ur3 = NAN_R
       us3 = NAN_R
       !---Ch 3
       if ( (FCDR%dre_over_da03(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da13(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da23(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da33(j,i) .ne. NAN_R) )THEN
          !
          ! For calibration parameters
          !
          ! Use covariance matrix
          !
          FCDR%Cmatrix3(1) = FCDR%dre_over_da03(j,i)
          FCDR%Cmatrix3(2) = FCDR%dre_over_da13(j,i)
          if( nparams3 .eq. 3 )then
             FCDR%Cmatrix3(3) = FCDR%dre_over_da33(j,i)
          endif
          FCDR%Omatrix3 = MATMUL(FCDR%Cmatrix3,Covar3)
          if( nparams3 .eq. 3 )then
             uc3 = SQRT(FCDR%Omatrix3(1)*FCDR%Cmatrix3(1)+&
                  FCDR%Omatrix3(2)*FCDR%Cmatrix3(2)+&
                  FCDR%Omatrix3(3)*FCDR%Cmatrix3(3))
          else
             uc3 = SQRT(FCDR%Omatrix3(1)*FCDR%Cmatrix3(1)+&
                  FCDR%Omatrix3(2)*FCDR%Cmatrix3(2))
          endif
          !
          ! Store for CURUC
          !
          FCDR%hu3(j,i) = uc3
          !
          ! Add in ICT T error term
          !
          uc3 = SQRT(uc3*uc3+(FCDR%dre_over_drict3(j,i)**2*FCDR%urict3_s(i)**2))
       ENDIF

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
             ! Don't add in ICT T terms here - make it common
             us3=sqrt((FCDR%dre_over_dcs3(j,i)**2*FCDR%ucs3(i)**2) &
                  +(FCDR%dre_over_dcict3(j,i)**2*FCDR%ucict3(i)**2)) 
!             us3=sqrt((FCDR%dre_over_drict3(j,i)**2*FCDR%urict3_s(i)**2) &
!                  +(FCDR%dre_over_dcs3(j,i)**2*FCDR%ucs3(i)**2) &
!                  +(FCDR%dre_over_dcict3(j,i)**2*FCDR%ucict3(i)**2)) 
          ENDIF
       end if

       !---Ch 4
       uc4 = NAN_R
       ur4 = NAN_R
       us4 = NAN_R
       if ( (FCDR%dre_over_da04(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da14(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da24(j,i) .ne. NAN_R) &
            .and. (FCDR%dre_over_da34(j,i) .ne. NAN_R) )THEN
          !
          ! For calibration parameters
          !
          ! Use covariance matrix
          !
          FCDR%Cmatrix(1) = FCDR%dre_over_da04(j,i)
          FCDR%Cmatrix(2) = FCDR%dre_over_da14(j,i)
          FCDR%Cmatrix(3) = FCDR%dre_over_da14(j,i)
          if( nparams .eq. 4 )then
             FCDR%Cmatrix(4) = FCDR%dre_over_da34(j,i)
          endif
          FCDR%Omatrix = MATMUL(FCDR%Cmatrix,Covar4)
          if( nparams .eq. 4 )then
             uc4 = SQRT(FCDR%Omatrix(1)*FCDR%Cmatrix(1)+&
                  FCDR%Omatrix(2)*FCDR%Cmatrix(2)+&
                  FCDR%Omatrix(3)*FCDR%Cmatrix(3)+&
                  FCDR%Omatrix(4)*FCDR%Cmatrix(4))
          else
             uc4 = SQRT(FCDR%Omatrix(1)*FCDR%Cmatrix(1)+&
                  FCDR%Omatrix(2)*FCDR%Cmatrix(2)+&
                  FCDR%Omatrix(3)*FCDR%Cmatrix(3))
          endif
          !
          ! Store for CURUC
          !
          FCDR%hu4(j,i) = uc4
          !
          ! Add in ICT T error term
          !
          uc4 = SQRT(uc4*uc4+(FCDR%dre_over_drict4(j,i)**2*FCDR%urict4_s(i)**2))
       ENDIF

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
             us4=sqrt((FCDR%dre_over_dcs4(j,i)**2*FCDR%ucs4(i)**2) &
                  +(FCDR%dre_over_dcict4(j,i)**2*FCDR%ucict4(i)**2)) 
!             us4=sqrt((FCDR%dre_over_drict4(j,i)**2*FCDR%urict4_s(i)**2) &
!                  +(FCDR%dre_over_dcs4(j,i)**2*FCDR%ucs4(i)**2) &
!                  +(FCDR%dre_over_dcict4(j,i)**2*FCDR%ucict4(i)**2)) 
          ENDIF
       end if

!MT: 06-12-2017: allocate ur5 and us5 first
!       ur5 = NAN_R
!       us5 = NAN_R

       IF( twelve_micron_there )THEN
          uc5 = NAN_R
          ur5 = NAN_R
          us5 = NAN_R
          if ( (FCDR%dre_over_da05(j,i) .ne. NAN_R) &
               .and. (FCDR%dre_over_da15(j,i) .ne. NAN_R) &
               .and. (FCDR%dre_over_da25(j,i) .ne. NAN_R) &
               .and. (FCDR%dre_over_da35(j,i) .ne. NAN_R) )THEN
             !
             ! For calibration parameters
             !
             ! Use covariance matrix
             !
             FCDR%Cmatrix(1) = FCDR%dre_over_da05(j,i)
             FCDR%Cmatrix(2) = FCDR%dre_over_da15(j,i)
             FCDR%Cmatrix(3) = FCDR%dre_over_da15(j,i)
             if( nparams .eq. 4 )then
                FCDR%Cmatrix(4) = FCDR%dre_over_da35(j,i)
             endif
             FCDR%Omatrix = MATMUL(FCDR%Cmatrix,Covar5)
             if( nparams .eq. 4 )then
                uc5 = SQRT(FCDR%Omatrix(1)*FCDR%Cmatrix(1)+&
                     FCDR%Omatrix(2)*FCDR%Cmatrix(2)+&
                     FCDR%Omatrix(3)*FCDR%Cmatrix(3)+&
                     FCDR%Omatrix(4)*FCDR%Cmatrix(4))
             else
                uc5 = SQRT(FCDR%Omatrix(1)*FCDR%Cmatrix(1)+&
                     FCDR%Omatrix(2)*FCDR%Cmatrix(2)+&
                     FCDR%Omatrix(3)*FCDR%Cmatrix(3))
             endif
             !
             ! Store for CURUC
             !
             FCDR%hu5(j,i) = uc5
             !
             ! Add in ICT T error term
             !
             uc5 = SQRT(uc5*uc5+(FCDR%dre_over_drict5(j,i)**2*&
                  FCDR%urict5_s(i)**2))
          ENDIF
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
                us5=sqrt((FCDR%dre_over_dcs5(j,i)**2*FCDR%ucs5(i)**2) &
                     +(FCDR%dre_over_dcict5(j,i)**2*FCDR%ucict5(i)**2)) 
!                us5=sqrt((FCDR%dre_over_drict5(j,i)**2*FCDR%urict5_s(i)**2) &
!                     +(FCDR%dre_over_dcs5(j,i)**2*FCDR%ucs5(i)**2) &
!                     +(FCDR%dre_over_dcict5(j,i)**2*FCDR%ucict5(i)**2)) 
             ENDIF
          end if
       ENDIF

       !---FIDUCEO : on convertit les radiances en BT
       if  ((outData%new_array3B(j,i) .gt. 0).and. &
            (outData%new_array3B(j,i) .ne. NAN_R) &
            .and. ur3 .ne. NAN_R .and. us3 .ne. NAN_R .and. &
            uc3 .ne. NAN_R ) then
          IF( outData%new_array3B(j,i) .eq. TINY )THEN
             FCDR%flag_no_detection(1,i) = 1
             FCDR%uc3(j,i) = 400.
             FCDR%ur3(j,i) = 400.
             FCDR%us3(j,i) = 400.
             FCDR%btf3(j,i) = TINY
             FCDR%dBT_over_dT3(j,i) = TINY
             FCDR%dBT_over_dcs3(j,i) = TINY
             FCDR%dBT_over_dcict3(j,i) = TINY
          ELSE
             FCDR%btf3(j,i)=convertBT(outData%new_array3B(j,i),&
                  DBLE(coefs1(5,1)), &
                  DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
             nrmax3=outData%new_array3B(j,i)+ur3
             trmax3=convertBT(nrmax3,DBLE(coefs1(5,1)), &
                  DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
             nrmin3=outData%new_array3B(j,i)-ur3
             two_sigma = outData%new_array3B(j,i)-2*ur3
             if( two_sigma .gt. 0. )then
                trmin3=convertBT(nrmin3,DBLE(coefs1(5,1)), &
                     DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
                FCDR%ur3(j,i)=(trmax3-trmin3)/2.
             else
                !             FCDR%ur3(j,i) = trmax3-FCDR%btf3(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(1,i) = 1
                FCDR%uc3(j,i) = 400.
                FCDR%ur3(j,i) = 400.
                FCDR%us3(j,i) = 400.
                FCDR%btf3(j,i) = TINY
                FCDR%dBT_over_dT3(j,i) = TINY
                FCDR%dBT_over_dcs3(j,i) = TINY
                FCDR%dBT_over_dcict3(j,i) = TINY
             ENDIF

             nsmax3=outData%new_array3B(j,i)+us3
             nsmin3=outData%new_array3B(j,i)-us3
             two_sigma = outData%new_array3B(j,i)-2*us3
             tsmax3=convertBT(nsmax3,DBLE(coefs1(5,1)), DBLE(coefs1(6,1)), &
                  DBLE(coefs1(7,1)))
             if( two_sigma .gt. 0. .and. NAN_R .ne. FCDR%ur3(j,i) )then
                tsmin3=convertBT(nsmin3,DBLE(coefs1(5,1)), DBLE(coefs1(6,1)), &
                     DBLE(coefs1(7,1)))
                FCDR%us3(j,i)=(tsmax3-tsmin3)/2.
             else
                !             FCDR%us3(j,i) = tsmax3-FCDR%btf3(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(1,i) = 1
                FCDR%uc3(j,i) = 400.
                FCDR%ur3(j,i) = 400.
                FCDR%us3(j,i) = 400.
                FCDR%btf3(j,i) = TINY
                FCDR%dBT_over_dT3(j,i) = TINY
                FCDR%dBT_over_dcs3(j,i) = TINY
                FCDR%dBT_over_dcict3(j,i) = TINY
             endif

             nsmax3=outData%new_array3B(j,i)+uc3
             nsmin3=outData%new_array3B(j,i)-uc3
             two_sigma = outData%new_array3B(j,i)-2*uc3
             tsmax3=convertBT(nsmax3,DBLE(coefs1(5,1)), DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
             if( two_sigma .gt. 0. .and. NAN_R .ne. FCDR%ur3(j,i) )then
                tsmin3=convertBT(nsmin3,DBLE(coefs1(5,1)), DBLE(coefs1(6,1)), DBLE(coefs1(7,1)))
                FCDR%uc3(j,i)=(tsmax3-tsmin3)/2.
             else
                !             FCDR%us3(j,i) = tsmax3-FCDR%btf3(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(1,i) = 1
                FCDR%uc3(j,i) = 400.
                FCDR%ur3(j,i) = 400.
                FCDR%us3(j,i) = 400.
                FCDR%btf3(j,i) = TINY
                FCDR%dBT_over_dT3(j,i) = TINY
                FCDR%dBT_over_dcs3(j,i) = TINY
                FCDR%dBT_over_dcict3(j,i) = TINY
             endif
          ENDIF
       end if

       if  ((outData%new_array4(j,i) .gt. 0).and.&
            (outData%new_array4(j,i) .ne. NAN_R) &
            .and. ur4 .ne. NAN_R .and. us4 .ne. NAN_R .and. &
            uc4 .ne. NAN_R ) then
          IF( outData%new_array4(j,i) .eq. TINY )THEN
             FCDR%flag_no_detection(2,i) = 1
             FCDR%uc4(j,i) = 400.
             FCDR%ur4(j,i) = 400.
             FCDR%us4(j,i) = 400.
             FCDR%btf4(j,i) = TINY
             FCDR%dBT_over_dT4(j,i) = TINY
             FCDR%dBT_over_dcs4(j,i) = TINY
             FCDR%dBT_over_dcict4(j,i) = TINY
          ELSE
             FCDR%btf4(j,i)=convertBT(outdata%new_array4(j,i),&
                  DBLE(coefs2(5,1)), &
                  DBLE(coefs2(6,1)), DBLE(coefs2(7,1)))
             nrmax4=outdata%new_array4(j,i)+ur4 
             trmax4=convertBT(nrmax4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                  DBLE(coefs2(7,1)))
             nrmin4=outdata%new_array4(j,i)-ur4
             two_sigma = outdata%new_array4(j,i)-2*ur4
             if( two_sigma .gt. 0 )then
                trmin4=convertBT(nrmin4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                     DBLE(coefs2(7,1)))
                FCDR%ur4(j,i)=(trmax4-trmin4)/2.
             else
                !             FCDR%ur4(j,i)=trmax4-FCDR%btf4(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(2,i) = 1
                FCDR%uc4(j,i) = 400.
                FCDR%ur4(j,i) = 400.
                FCDR%us4(j,i) = 400.
                FCDR%btf4(j,i) = TINY
                FCDR%dBT_over_dT4(j,i) = TINY
                FCDR%dBT_over_dcs4(j,i) = TINY
                FCDR%dBT_over_dcict4(j,i) = TINY
             endif

             nsmax4=outdata%new_array4(j,i)+us4
             tsmax4=convertBT(nsmax4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                  DBLE(coefs2(7,1)))
             nsmin4=outdata%new_array4(j,i)-us4
             two_sigma = outdata%new_array4(j,i)-2*us4
             if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur4(j,i) )then
                tsmin4=convertBT(nsmin4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                     DBLE(coefs2(7,1)))
                FCDR%us4(j,i)=(tsmax4-tsmin4)/2.
             else
                !             FCDR%us4(j,i)=tsmax4-FCDR%btf4(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(2,i) = 1
                FCDR%uc4(j,i) = 400.
                FCDR%ur4(j,i) = 400.
                FCDR%us4(j,i) = 400.
                FCDR%btf4(j,i) = TINY
                FCDR%dBT_over_dT4(j,i) = TINY
                FCDR%dBT_over_dcs4(j,i) = TINY
                FCDR%dBT_over_dcict4(j,i) = TINY
             endif
             
             nsmax4=outdata%new_array4(j,i)+uc4
             tsmax4=convertBT(nsmax4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                  DBLE(coefs2(7,1)))
             nsmin4=outdata%new_array4(j,i)-uc4
             two_sigma = outdata%new_array4(j,i)-2*uc4
             if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur4(j,i) )then
                tsmin4=convertBT(nsmin4,DBLE(coefs2(5,1)), DBLE(coefs2(6,1)), &
                     DBLE(coefs2(7,1)))
                FCDR%uc4(j,i)=(tsmax4-tsmin4)/2.
             else
                !             FCDR%us4(j,i)=tsmax4-FCDR%btf4(j,i)
                ! Set data to bad as no 1 sigma detection available
                FCDR%flag_no_detection(2,i) = 1
                FCDR%uc4(j,i) = 400.
                FCDR%ur4(j,i) = 400.
                FCDR%us4(j,i) = 400.
                FCDR%btf4(j,i) = TINY
                FCDR%dBT_over_dT4(j,i) = TINY
                FCDR%dBT_over_dcs4(j,i) = TINY
                FCDR%dBT_over_dcict4(j,i) = TINY
             endif
          ENDIF
       end if

       IF( twelve_micron_there )THEN
          if  ((outData%new_array5(j,i) .gt. 0).and. &
               (outData%new_array5(j,i) .ne. NAN_R) &
               .and. ur5 .ne. NAN_R .and. us5 .ne. NAN_R ) then
             IF( outData%new_array5(j,i) .eq. TINY )THEN
                FCDR%flag_no_detection(3,i) = 1
                FCDR%uc5(j,i) = 400.
                FCDR%ur5(j,i) = 400.
                FCDR%us5(j,i) = 400.
                FCDR%btf5(j,i) = TINY
                FCDR%dBT_over_dT5(j,i) = TINY
                FCDR%dBT_over_dcs5(j,i) = TINY
                FCDR%dBT_over_dcict5(j,i) = TINY
             ELSE
                FCDR%btf5(j,i)=convertBT(outdata%new_array5(j,i),&
                     DBLE(coefs3(5,1)), DBLE(coefs3(6,1)), DBLE(coefs3(7,1)))
                nrmax5=outdata%new_array5(j,i)+ur5
                trmax5=convertBT(nrmax5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)), &
                     DBLE(coefs3(7,1)))
                nrmin5=outdata%new_array5(j,i)-ur5
                trmin5=convertBT(nrmin5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)), &
                     DBLE(coefs3(7,1)))
                two_sigma = outdata%new_array5(j,i)-2*ur5
                if( two_sigma .gt. 0 )then
                   trmin5=convertBT(nrmin5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)),&
                        DBLE(coefs3(7,1)))
                   FCDR%ur5(j,i)=(trmax5-trmin5)/2.
                else
                   !                FCDR%ur5(j,i)=trmax5-FCDR%btf5(j,i)
                   ! Set data to bad as no 1 sigma detection available
                   FCDR%flag_no_detection(3,i) = 1
                   FCDR%uc5(j,i) = 400.
                   FCDR%ur5(j,i) = 400.
                   FCDR%us5(j,i) = 400.
                   FCDR%btf5(j,i) = TINY
                   FCDR%dBT_over_dT5(j,i) = TINY
                   FCDR%dBT_over_dcs5(j,i) = TINY
                   FCDR%dBT_over_dcict5(j,i) = TINY
                endif
                nsmax5=outdata%new_array5(j,i)+us5
                tsmax5=convertBT(nsmax5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)),&
                     DBLE(coefs3(7,1)))
                nsmin5=outdata%new_array5(j,i)-us5
                two_sigma = outdata%new_array5(j,i)-2*us5
                if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur5(j,i) )then
                   tsmin5=convertBT(nsmin5,DBLE(coefs3(5,1)), &
                        DBLE(coefs3(6,1)), DBLE(coefs3(7,1)))
                   FCDR%us5(j,i)=(tsmax5-tsmin5)/2.
                else
                   !                FCDR%us5(j,i)=tsmax5-FCDR%btf5(j,i)
                   ! Set data to bad as no 1 sigma detection available
                   FCDR%flag_no_detection(3,i) = 1
                   FCDR%uc5(j,i) = 400.
                   FCDR%ur5(j,i) = 400.
                   FCDR%us5(j,i) = 400.
                   FCDR%btf5(j,i) = TINY
                   FCDR%dBT_over_dT5(j,i) = TINY
                   FCDR%dBT_over_dcs5(j,i) = TINY
                   FCDR%dBT_over_dcict5(j,i) = TINY
                endif
                nsmax5=outdata%new_array5(j,i)+uc5
                tsmax5=convertBT(nsmax5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)),&
                     DBLE(coefs3(7,1)))
                nsmin5=outdata%new_array5(j,i)-uc5
                two_sigma = outdata%new_array5(j,i)-2*us5
                if( two_sigma .gt. 0  .and. NAN_R .ne. FCDR%ur5(j,i) )then
                   tsmin5=convertBT(nsmin5,DBLE(coefs3(5,1)), DBLE(coefs3(6,1)),&
                        DBLE(coefs3(7,1)))
                   FCDR%uc5(j,i)=(tsmax5-tsmin5)/2.
                else
                   !                FCDR%us5(j,i)=tsmax5-FCDR%btf5(j,i)
                   ! Set data to bad as no 1 sigma detection available
                   FCDR%flag_no_detection(3,i) = 1
                   FCDR%uc5(j,i) = 400.
                   FCDR%ur5(j,i) = 400.
                   FCDR%us5(j,i) = 400.
                   FCDR%btf5(j,i) = TINY
                   FCDR%dBT_over_dT5(j,i) = TINY
                   FCDR%dBT_over_dcs5(j,i) = TINY
                   FCDR%dBT_over_dcict5(j,i) = TINY
                endif
             ENDIF
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
  ! Write GBCS dtime to netcdf
  !
  SUBROUTINE Write_GBCS_time(ncid,dims_time,stime,units)
    
    INTEGER, INTENT(IN) :: ncid
    INTEGER, INTENT(IN) :: dims_time
    INTEGER, INTENT(IN) :: stime
    CHARACTER(LEN=*), INTENT(IN) :: units

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims1(1)
    INTEGER :: varid
    INTEGER :: data(1)

    !
    ! Fill with values 
    !
    data(1) = stime    

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims1 = (/dims_time/)
    stat = NF90_DEF_VAR(ncid,'time',NF90_INT,dims1,varid)
    CALL check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name','reference time of sst file')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'standard_name','time')
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'calendar','gregorian')
    call check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

  END SUBROUTINE Write_GBCS_time

  !
  ! Write GBCS dtime to netcdf
  !
  SUBROUTINE Write_GBCS_dtime(ncid,dim_nx,dim_ny,dim_time,nx,ny,stime,indata,&
       long_name,units,coordinates,fill_value)
    
    INTEGER, INTENT(IN) :: ncid
    INTEGER, INTENT(IN) :: dim_nx
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: dim_time
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    INTEGER, INTENT(IN) :: stime
    REAL, INTENT(IN) :: indata(ny)
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    CHARACTER(LEN=*), INTENT(IN) :: coordinates
    REAL, INTENT(IN) :: fill_value

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims3(3)
    INTEGER :: varid
    REAL, ALLOCATABLE :: data(:,:,:)

    ALLOCATE(data(nx,ny,1),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_dtime',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with values 
    !
    DO I=1,ny
       DO J=1,nx
          IF( indata(I) .gt. -1e20 )THEN
             data(J,I,1) = indata(I)-stime
          ELSE
             data(J,I,1) = fill_value
          ENDIF
       END DO
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims3 = (/dim_nx,dim_ny,dim_time/)
    stat = NF90_DEF_VAR(ncid,'dtime',NF90_FLOAT,dims3,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'coordinates',TRIM(coordinates))
    call check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_dtime
 
  !
  ! Write GBCS float to netcdf
  !
  SUBROUTINE Write_GBCS_Float_1d(ncid,name,dim_ny,ny,indata,fill_value,&
       long_name,standard_name,units,valid_min,valid_max)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: ny
    REAL, INTENT(IN) :: indata(ny)
    REAL, INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: standard_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    REAL, INTENT(IN) :: valid_min
    REAL, INTENT(IN) :: valid_max

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims1(1)
    INTEGER :: varid
    REAL, ALLOCATABLE :: data(:)

    ALLOCATE(data(ny),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_Float',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    DO I=1,ny
       IF( indata(I) .lt. valid_min .or. indata(I) .gt. valid_max )THEN
          data(I) = fill_value
       ELSE
          data(I) = indata(I)
       ENDIF
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims1 = (/dim_ny/)
    stat = NF90_DEF_VAR(ncid,name,NF90_FLOAT,dims1,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'standard_name',TRIM(standard_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_min',valid_min)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_max',valid_max)
    call check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float_1d

  !
  ! Write GBCS float to netcdf
  !
  SUBROUTINE Write_GBCS_Float(ncid,name,dim_nx,dim_ny,nx,ny,indata,fill_value,&
       long_name,standard_name,units,valid_min,valid_max,reference_datum)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_nx
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    REAL, INTENT(IN) :: indata(nx,ny)
    REAL, INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: standard_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    REAL, INTENT(IN) :: valid_min
    REAL, INTENT(IN) :: valid_max
    CHARACTER(LEN=*), INTENT(IN) :: reference_datum

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims2(2)
    INTEGER :: varid
    REAL, ALLOCATABLE :: data(:,:)

    ALLOCATE(data(nx,ny),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_Float',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    DO I=1,ny
       DO J=1,nx
          IF( indata(J,I) .lt. valid_min .or. indata(J,I) .gt. valid_max )THEN
             data(J,I) = fill_value
          ELSE
             data(J,I) = indata(J,I)
          ENDIF
       END DO
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims2 = (/dim_nx,dim_ny/)
    stat = NF90_DEF_VAR(ncid,name,NF90_FLOAT,dims2,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'standard_name',TRIM(standard_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_min',valid_min)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_max',valid_max)
    call check(stat)
    IF( '' .ne. reference_datum )THEN
       stat = NF90_PUT_ATT(ncid,varid,'reference_datum',TRIM(reference_datum))
       call check(stat)
    ENDIF

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float

  !
  ! Write GBCS float to netcdf
  !
  SUBROUTINE Write_GBCS_Float_Time(ncid,name,dim_ny,dim_time,&
       ny,indata,fill_value,long_name,units,valid_min,valid_max,scale,offset)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: dim_time
    INTEGER, INTENT(IN) :: ny
    REAL, INTENT(IN) :: indata(ny)
    INTEGER(GbcsInt2), INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(GbcsInt2), INTENT(IN) :: valid_min
    INTEGER(GbcsInt2), INTENT(IN) :: valid_max
    REAL, INTENT(IN) :: scale
    REAL, INTENT(IN) :: offset

    ! Local variables
    INTEGER :: I
    INTEGER :: stat
    INTEGER :: dims2(2)
    INTEGER :: varid
    REAL :: valid_min_f
    REAL :: valid_max_f
    INTEGER(GbcsInt2), ALLOCATABLE :: data(:,:)

    ALLOCATE(data(ny,1),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_Float',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    valid_min_f = valid_min*scale+offset
    valid_max_f = valid_max*scale+offset
    DO I=1,ny
       IF( indata(I) .lt. valid_min_f .or. indata(I) .gt. valid_max_f )THEN
          data(I,1) = fill_value
       ELSE
          data(I,1) = INT((indata(I)-offset)/scale,KIND=GbcsInt2)
       ENDIF
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims2 = (/dim_ny,dim_time/)
    stat = NF90_DEF_VAR(ncid,name,NF90_SHORT,dims2,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_min',valid_min)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_max',valid_max)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'add_offset',offset)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'scale_factor',scale)
    call check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float_Time

  SUBROUTINE Write_GBCS_Float_Time_Int(ncid,name,dim_ny,dim_time,&
       ny,indata,fill_value,long_name,valid_min)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: dim_time
    INTEGER, INTENT(IN) :: ny
    INTEGER, INTENT(IN) :: indata(ny)
    INTEGER(GbcsInt2), INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    INTEGER(GbcsInt2), INTENT(IN) :: valid_min

    ! Local variables
    INTEGER :: I
    INTEGER :: stat
    INTEGER :: dims2(2)
    INTEGER :: varid
    INTEGER(GbcsInt2), ALLOCATABLE :: data(:,:)

    ALLOCATE(data(ny,1),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_Float',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    DO I=1,ny
       IF( indata(I) .lt. valid_min )THEN
          data(I,1) = fill_value
       ELSE
          data(I,1) = INT(indata(I),KIND=GbcsInt2)
       ENDIF
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims2 = (/dim_ny,dim_time/)
    stat = NF90_DEF_VAR(ncid,name,NF90_SHORT,dims2,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_min',valid_min)
    call check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float_Time_Int

  !
  ! Write GBCS scaled int (to float) to netcdf
  !
  SUBROUTINE Write_GBCS_Float_to_Int(ncid,name,dim_nx,dim_ny,nx,ny,&
       indata,fill_value,long_name,standard_name,units,valid_min,valid_max,&
       coordinates,scale,offset,comment)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_nx
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    REAL, INTENT(IN) :: indata(nx,ny)
    INTEGER(GbcsInt2), INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: standard_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(GbcsInt2), INTENT(IN) :: valid_min
    INTEGER(GbcsInt2), INTENT(IN) :: valid_max
    CHARACTER(LEN=*), INTENT(IN) :: coordinates
    REAL, INTENT(IN) :: scale
    REAL, INTENT(IN) :: offset
    CHARACTER(LEN=*), INTENT(IN) :: comment

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims2(2)
    INTEGER :: varid
    REAL :: valid_min_f, valid_max_f
    INTEGER(GbcsInt2), ALLOCATABLE :: data(:,:)

    ALLOCATE(data(nx,ny),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data','Write_GBCS_Float',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    valid_min_f = valid_min*scale+offset
    valid_max_f = valid_max*scale+offset
    DO I=1,ny
       DO J=1,nx
          IF( indata(J,I) .lt. valid_min_f .or. &
               indata(J,I) .gt. valid_max_f )THEN
             data(J,I) = fill_value
          ELSE
             data(J,I) = INT((indata(J,I)-offset)/scale,KIND=GbcsInt2)
          ENDIF
       END DO
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims2 = (/dim_nx,dim_ny/)
    stat = NF90_DEF_VAR(ncid,name,NF90_SHORT,dims2,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    IF( "" .ne. standard_name )THEN
       stat = NF90_PUT_ATT(ncid,varid,'standard_name',TRIM(standard_name))
       call check(stat)
    ENDIF
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'add_offset',offset)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'scale_factor',scale)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_min',valid_min)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_max',valid_max)
    call check(stat)
    IF( "" .ne. coordinates )THEN
       stat = NF90_PUT_ATT(ncid,varid,'coordinates',TRIM(coordinates))
       call check(stat)
    ENDIF
    IF( "" .ne. comment )THEN
       stat = NF90_PUT_ATT(ncid,varid,'comment',TRIM(comment))
       call check(stat)
    ENDIF
    
    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float_to_Int

  !
  ! Write GBCS scaled int (to float) to netcdf
  !
  SUBROUTINE Write_GBCS_Float_to_Int_3D(ncid,name,dim_nx,dim_ny,dim_time,nx,ny,&
       indata,fill_value,long_name,standard_name,units,valid_min,valid_max,&
       coordinates,scale,offset,comment)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_nx
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: dim_time
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    REAL, INTENT(IN) :: indata(nx,ny)
    INTEGER(GbcsInt2), INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: standard_name
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(GbcsInt2), INTENT(IN) :: valid_min
    INTEGER(GbcsInt2), INTENT(IN) :: valid_max
    CHARACTER(LEN=*), INTENT(IN) :: coordinates
    REAL, INTENT(IN) :: scale
    REAL, INTENT(IN) :: offset
    CHARACTER(LEN=*), INTENT(IN) :: comment

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims3(3)
    INTEGER :: varid
    REAL :: valid_min_f, valid_max_f
    INTEGER(GbcsInt2), ALLOCATABLE :: data(:,:,:)

    ALLOCATE(data(nx,ny,1),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data',&
            'Write_GBCS_Float_to_Int_3D',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Fill with fill values (GBCS)
    !
    valid_min_f = valid_min*scale+offset
    valid_max_f = valid_max*scale+offset
    DO I=1,ny
       DO J=1,nx
          IF( indata(J,I) .lt. valid_min_f .or. &
               indata(J,I) .gt. valid_max_f )THEN
             data(J,I,1) = fill_value
          ELSE
             data(J,I,1) = INT((indata(J,I)-offset)/scale,KIND=GbcsInt2)
          ENDIF
       END DO
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims3 = (/dim_nx,dim_ny,dim_time/)
    stat = NF90_DEF_VAR(ncid,name,NF90_SHORT,dims3,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    IF( "" .ne. standard_name )THEN
       stat = NF90_PUT_ATT(ncid,varid,'standard_name',TRIM(standard_name))
       call check(stat)
    ENDIF
    stat = NF90_PUT_ATT(ncid,varid,'units',TRIM(units))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'add_offset',offset)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'scale_factor',scale)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'valid_max',valid_max)
    call check(stat)
    IF( "" .ne. coordinates )THEN
       stat = NF90_PUT_ATT(ncid,varid,'coordinates',TRIM(coordinates))
       call check(stat)
    ENDIF
    IF( "" .ne. comment )THEN
       stat = NF90_PUT_ATT(ncid,varid,'comment',TRIM(comment))
       call check(stat)
    ENDIF
    
    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Float_to_Int_3D

  SUBROUTINE Write_GBCS_Byte(ncid,name,dim_nx,dim_ny,dim_time,nx,ny,&
       indata,fill_value,long_name,flag_meanings,nflags,flag_values,&
       valid_min,valid_max,comment,scale,offset)
    
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: dim_nx
    INTEGER, INTENT(IN) :: dim_ny
    INTEGER, INTENT(IN) :: dim_time
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    INTEGER(GbcsInt1), INTENT(IN) :: indata(ny)
    INTEGER(GbcsInt1), INTENT(IN) :: fill_value
    CHARACTER(LEN=*), INTENT(IN) :: long_name
    CHARACTER(LEN=*), INTENT(IN) :: flag_meanings
    INTEGER, INTENT(IN) :: nflags
    INTEGER(GbcsInt1), INTENT(IN) :: flag_values(nflags)
    INTEGER(GbcsInt1), INTENT(IN) :: valid_min
    INTEGER(GbcsInt1), INTENT(IN) :: valid_max
    CHARACTER(LEN=*), INTENT(IN) :: comment
    REAL, INTENT(IN), OPTIONAL :: scale
    REAL, INTENT(IN), OPTIONAL :: offset

    ! Local variables
    INTEGER :: I,J
    INTEGER :: stat
    INTEGER :: dims3(3)
    INTEGER :: varid
    INTEGER(GbcsInt1), ALLOCATABLE :: data(:,:,:)
    
    ALLOCATE(data(nx,ny,1),STAT=STAT)
    IF( 0 .ne. STAT )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate data',&
            'Write_GBCS_Byte',&
            'fiduceo_uncertainties.f90')
    ENDIF

    DO I=1,ny
       data(:,I,1) = indata(I)
    END DO

    !
    ! Go into define mode
    !
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    !
    ! Define variable
    !
    dims3 = (/dim_nx,dim_ny,dim_time/)
    stat = NF90_DEF_VAR(ncid,name,NF90_BYTE,dims3,varid)
    CALL check(stat)
    stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
    call check(stat)
    stat = NF90_DEF_VAR_FILL(ncid, varid, 0, fill_value)
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'long_name',TRIM(long_name))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'flag_meaings',TRIM(flag_meanings))
    call check(stat)
    stat = NF90_PUT_ATT(ncid,varid,'flag_values',flag_values)
    call check(stat)
    IF( "" .ne. comment )THEN
       stat = NF90_PUT_ATT(ncid,varid,'comment',TRIM(comment))
       call check(stat)
    ENDIF

    IF( PRESENT(scale) )THEN
       stat = NF90_PUT_ATT(ncid,varid,'scale_factor',scale)
       call check(stat)
    ENDIF
    IF( PRESENT(offset) )THEN
       stat = NF90_PUT_ATT(ncid,varid,'add_offset',offset)
       call check(stat)
    ENDIF
    
    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_VAR(ncid,varid,data)
    CALL check(stat)

    DEALLOCATE(data)

  END SUBROUTINE Write_GBCS_Byte

  SUBROUTINE Get_DateTime(AVHRR,stime,time_secs,start_year)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    INTEGER, INTENT(OUT) :: stime
    REAL, INTENT(OUT), ALLOCATABLE :: time_secs(:)
    INTEGER, INTENT(IN) :: start_year

    ! Local variables
    INTEGER :: I
    INTEGER :: stat
    TYPE(DateTime) :: date_time
    REAL(GbcsDble) :: jd
    REAL(GbcsDble) :: jd_start

    ALLOCATE(time_secs(AVHRR%arraySize),STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate time_secs','Get_DateTime',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Get since time
    !
    date_time%year = start_year
    date_time%month = 1
    date_time%day = 1
    date_time%hour = 0
    date_time%minute = 0
    date_time%seconds = 0
    date_time%sec1000 = 0
    date_time%utc_offset = 0
    jd_start = Date_to_Jd(date_time)

    !
    ! Get start time
    !
    CALL Convert_to_Datetime(AVHRR,1,date_time)
    jd = Date_to_Jd(date_time)
    stime = INT((jd-jd_start)*86400.d0)

    !
    ! Now loop round for times
    !
    DO I=1,AVHRR%arraySize
       IF( AVHRR%year(I) .gt. 0 .and. AVHRR%month(I) .gt. 0 .and. &
            AVHRR%day(i) .gt. 0 )THEN
          CALL Convert_to_Datetime(AVHRR,I,date_time)
          jd = Date_to_Jd(date_time)
          time_secs(I) = INT((jd-jd_start)*86400.d0)
       ELSE
          time_secs(I) = NAN_R
       ENDIF
    END DO

  END SUBROUTINE Get_DateTime

  SUBROUTINE Convert_to_DateTime(AVHRR,I,date_time)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    INTEGER, INTENT(IN) :: I
    TYPE(DateTime), INTENT(OUT) :: date_time

    ! Local variables
    REAL :: temp

    date_time%year = AVHRR%year(I)
    date_time%month = AVHRR%month(I)
    date_time%day = AVHRR%day(I)
    
    date_time%hour = INT(AVHRR%hours(I))
    temp = (AVHRR%hours(I)-date_time%hour)*60.
    date_time%minute = INT(temp)
    temp = (temp - date_time%minute)*60.
    date_time%seconds = INT(temp)
    date_time%sec1000 = INT((temp - date_time%seconds)*1000.)
    date_time%utc_offset = 0

  END SUBROUTINE Convert_to_DateTime

  !
  ! Make flags
  !
  SUBROUTINE Make_Flags(AVHRR,nflags,flags,flags_meanings,flag_values,&
       valid_min,valid_max)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    INTEGER, INTENT(OUT) :: nflags
    INTEGER(GbcsInt1), ALLOCATABLE, INTENT(OUT) :: flags(:)
    CHARACTER(LEN=*), INTENT(OUT) :: flags_meanings
    INTEGER(GbcsInt1), ALLOCATABLE, INTENT(OUT) :: flag_values(:)
    INTEGER(GbcsInt1), INTENT(OUT) :: valid_min
    INTEGER(GbcsInt1), INTENT(OUT) :: valid_max

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT

    nflags = 7
    flags_meanings = 'bad_navigation bad_calibration bad_timing missing_line &
         &solar_contamination_3B solar_contamination_4 solar_contamination_5'
    valid_min = 0_GbcsInt1
    ALLOCATE(flags(AVHRR%arraySize),flag_values(nflags),STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate flag_values',&
            'Make_Flags','fiduceo_uncertainties.f90')
    ENDIF
    DO I=1,nflags
       flag_values(I) = IBSET(0,I-1)       
    END DO

    valid_max = flag_values(nflags)

    flags = 0
    DO I=1,AVHRR%arraySize
       IF( AVHRR%badNavigation(I) )THEN
          flags(I) = IBSET(flags(I),0)
       ENDIF
       IF( AVHRR%badCalibration(I) )THEN
          flags(I) = IBSET(flags(I),1)
       ENDIF
       IF( AVHRR%badTime(I) )THEN
          flags(I) = IBSET(flags(I),2)
       ENDIF
       IF( AVHRR%missingLines(I) )THEN
          flags(I) = IBSET(flags(I),3)
       ENDIF
       IF( AVHRR%solar_contamination_3b(I) )THEN
          flags(I) = IBSET(flags(I),4)
       ENDIF
       IF( AVHRR%solar_contamination_4(I) )THEN
          flags(I) = IBSET(flags(I),5)
       ENDIF
       IF( AVHRR%solar_contamination_5(I) )THEN
          flags(I) = IBSET(flags(I),6)
       ENDIF
    END DO

  END SUBROUTINE Make_Flags

  !
  ! Strip out just GAC filename from full path
  !
  SUBROUTINE get_filename(infile,ofile)
    
    CHARACTER(LEN=*), INTENT(IN) :: infile
    CHARACTER(LEN=*), INTENT(OUT) :: ofile

    ! Local variables
    INTEGER :: I
    INTEGER :: POS
    
    POS = -1
    Loop: DO I=LEN_TRIM(infile),1,-1
       IF( '/' .eq. infile(I:I) )THEN
          POS=I+1
          EXIT Loop
       ENDIF
    END DO Loop

    IF( -1 .eq. POS )THEN
       ofile = infile(1:42)
    ELSE
       ofile = infile(POS:POS+41)
    ENDIF

  END SUBROUTINE get_filename

  !
  ! Write global attributes
  !
  SUBROUTINE Write_GBCS_L1C_Global(AVHRR,ncid,uuid,noaa_string)
    
    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: uuid
    CHARACTER(LEN=*), INTENT(IN) :: noaa_string
    
    ! Local variables
    INTEGER :: I
    INTEGER :: stat
    CHARACTER(LEN=256) :: string
    CHARACTER(LEN=512) :: ofile
    TYPE(DateTime) :: date_time
    TYPE(DateTime) :: date_time2
    
    WRITE(*,'('' ****** Writing CCI L1C format ******'')')
    stat = NF90_REDEF(ncid)
    CALL check(stat)

    WRITE(string,&
         '(''AVHRR Pre-Processing: '',a,'' L1C product (FIDUCEO based)'')')&
         TRIM(noaa_string)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',TRIM(string))
    CALL check(stat)

    WRITE(string,&
         '(a,''-ESACCI-L1C-vFIDUCEO'')')TRIM(noaa_string)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'id',TRIM(string))
    CALL check(stat)
    
    WRITE(string,&
         '(a,'' L1C product from the FIDUCEO project'')')TRIM(noaa_string)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'summary',TRIM(string))
    CALL check(stat)
    
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'references','http://www.fiduceo,eu')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution','FIDUCEO')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'history','Created using FIDUCEO code')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'license','This dataset is released &
         &for use under CC-BY licence &
         &(https://creativecommons.org/licenses/by/4.0/) and was developed &
         &in the EC FIDUCEO project \"Fidelity and Uncertainty in Climate &
         &Data Records from Earth Observations\". Grant Agreement: 638822.')
    CALL check(stat)

    stat = NF90_ENDDEF(ncid)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'product_version','v0.5 pre-beta')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'uuid',TRIM(uuid))
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'tracking_id',TRIM(uuid))
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'netcdf_version_id',&
         TRIM(nf90_inq_libvers()))
    CALL check(stat)

    CALL Get_Time(date_time)
    CALL write_isodate(string,date_time)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'date_created',&
         TRIM(string))
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'file_quality_level',3)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'spatial_resolution',&
         '4.0 km at nadir')
    CALL check(stat)

    CALL Convert_to_Datetime(AVHRR,1,date_time)
    CALL write_isodate(string,date_time)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'start_time',TRIM(string))
    CALL check(stat)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'time_coverage_start',TRIM(string))
    CALL check(stat)

    CALL Convert_to_Datetime(AVHRR,AVHRR%arraySize,date_time2)
    CALL write_isodate(string,date_time2)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'stop_time',TRIM(string))
    CALL check(stat)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'time_coverage_end',TRIM(string))
    CALL check(stat)

    CALL write_isoduration(string,date_time,date_time2)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'time_coverage_duration',TRIM(string))
    CALL check(stat)

    date_time%year = 2010
    date_time%month = 1
    date_time%day = 1
    date_time%hour = 0
    date_time%minute = 0
    date_time%seconds = 0
    date_time%sec1000 = 0
    date_time%utc_offset = 0
    date_time2%year = 2010
    date_time2%month = 1
    date_time2%day = 1
    date_time2%hour = 0
    date_time2%minute = 0
    date_time2%seconds = 1
    date_time2%sec1000 = 0
    date_time2%utc_offset = 0
    CALL write_isoduration(string,date_time,date_time2)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'time_coverage_resolution',&
         TRIM(string))
    CALL check(stat)

    WRITE(string,'(a,''-NOAA-L1-v1'')')TRIM(noaa_string)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'source',&
         TRIM(string))
    CALL check(stat)

    CALL noaa_name(AVHRR,string,noaa=.TRUE.)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'platform',&
         TRIM(string))
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'sensor','AVHRR_GAC')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'Metadata_Conventions',&
         'Unidata Dataset Discovery v1.0')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'metadata_link',&
         'http://www.esa-cci.org')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'keywords',&
         'AVHRR > L1C')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'standard_name_vocabularly',&
         'NetCDF Climate and Forecast (CF) Metadata Convention')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lat_units',&
         'degrees_north')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lat_resolution',&
         0.04)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lon_units',&
         'degrees_east')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lon_resolution',&
         0.04)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'acknowledgement',&
         'NOAA GAC Data. Processing funded by H2020 (EC)')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'creator_name',&
         'FIDUCEO')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'creator_email',&
         'Fiduceo-coordinator@lists.reading.ac.uk')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'creator_url',&
         'http://www.fiduceo.eu')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'creator_processing_institution',&
         'These data were produced on the JASMIN infrastructure at STFC as &
         &part of the H2020 FIDUCEO project')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'project',&
         'EC H2020 FIDUCEO project')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'northernmost_latitude',90.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'southernmost_latitude',90.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'easternmost_longitude',180.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'westernmost_longitude',-180.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lat_min',-90.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lat_max',90.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lon_min',-180.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'geospatial_lon_max',180.)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'processing_level','L1C')
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'cdm_data_type','swath')
    CALL check(stat)

    string = 'merged files: '
    DO I=1,AVHRR%norig_l1b
       string = TRIM(string)//' '//TRIM(AVHRR%orig_l1b(I))
    END DO
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'source_file',TRIM(string))
    CALL check(stat)
    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'gac_file',TRIM(string))
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'lines_truncated_start',1)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'lines_truncated_end',AVHRR%arraySize)
    CALL check(stat)

    stat = NF90_PUT_ATT(ncid,NF90_GLOBAL,'orbital_temperature',&
         AVHRR%orbital_temperature)
    CALL check(stat)    

  END SUBROUTINE Write_GBCS_L1C_Global

  !
  ! Routine to write out the GBCS CCI L1C format for use with the GBCS
  !
  SUBROUTINE Write_GBCS_L1C(AVHRR,FCDR,uuid,twelve_micron_there,&
       output_cal)

    TYPE(AVHRR_Data), INTENT(IN) :: AVHRR
    TYPE(FIDUCEO_Data) :: FCDR
    CHARACTER(LEN=*), INTENT(IN) :: uuid
    LOGICAL, INTENT(IN) :: twelve_micron_there
    LOGICAL, INTENT(IN) :: output_cal

    ! Local variables
    INTEGER :: stat
    INTEGER :: hour
    INTEGER :: minute
    INTEGER :: second
    INTEGER :: ncid
    INTEGER :: dimid_nx
    INTEGER :: dimid_ny
    INTEGER :: dimid_time
    INTEGER :: dimid_band
    INTEGER :: dims1(1)
    INTEGER :: dims2(2)
    INTEGER :: stime
    INTEGER :: nx
    INTEGER :: ny
    INTEGER :: nflags
    INTEGER :: varid
    INTEGER(GbcsInt1), ALLOCATABLE :: flags(:)
    INTEGER(GbcsInt1), ALLOCATABLE :: flag_values(:)
    INTEGER(GbcsInt1) :: valid_min
    INTEGER(GbcsInt1) :: valid_max
    REAL :: temp
    REAL, ALLOCATABLE :: noise(:,:)
    REAL, ALLOCATABLE :: time_secs(:)
    CHARACTER(LEN=10) :: noaa_string
    CHARACTER(LEN=256) :: output_filename
    CHARACTER(LEN=256) :: flags_meaning
    REAL, ALLOCATABLE :: lon(:,:)
    
    nx = AVHRR%nelem
    ny = AVHRR%arraySize
    ALLOCATE(noise(AVHRR%nelem,AVHRR%arraySize),STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate noise','Write_GBCS_L1C',&
            'fiduceo_uncertainties.f90')
    ENDIF

    !
    ! Parse start time/AVHRR_No to filename
    !
    CALL NOAA_Name(AVHRR,noaa_string,noaa=.FALSE.)
    IF( 0 .gt. AVHRR%year(1) .or. 0 .gt. AVHRR%month(1) .or. &
         0 .gt. AVHRR%day(1) .or. 0 .gt. AVHRR%hours(1) )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot get start time','Write_GBCS_L1C',&
            'fiduceo_uncertainties.f90')
    ENDIF
    hour = INT(AVHRR%hours(1))
    temp = (AVHRR%hours(1)-hour)*60.
    minute = INT(temp)
    second = INT((temp-minute)*60.)
    WRITE(output_filename,'(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,&
         &''-ESACCI-L1C-'',a,''-fv01.0.nc'')')AVHRR%year(1),&
         AVHRR%month(1),AVHRR%day(1),hour,minute,second,&
         TRIM(noaa_string)

    !
    ! Open file
    !
    stat = NF90_CREATE(output_filename,IOR(NF90_HDF5,NF90_CLOBBER),ncid)
    call check(stat)
    
    stat = NF90_DEF_DIM(ncid,'ni',AVHRR%nelem,dimid_nx)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'nj',AVHRR%arraySize,dimid_ny)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'time',1,dimid_time)
    call check(stat)

    stat = NF90_DEF_DIM(ncid,'band',3,dimid_band)
    call check(stat)

    !
    ! Make sure we are out of define mode
    !
    stat = NF90_ENDDEF(ncid)
    call check(stat)

    !
    ! Write variables
    !
    CALL Write_GBCS_Float(ncid,'lat',dimid_nx,dimid_ny,nx,ny,AVHRR%lat,&
         -32768.,'Latitude coordinates','latitude','degrees_north',&
         -90.,90.,'geographical coordinates, WGS84 projection')
    !
    ! make sure -180,180
    !
    ALLOCATE(lon(nx,ny),STAT=stat)
    IF( 0 .ne. stat )THEN
       CALL Gbcs_Critical(.TRUE.,'Cannot allocate lon','Write_GBCS_L1C',&
            'fiduceo_uncertainties.f90')
    ENDIF
    WHERE( AVHRR%lon .gt. 180. )
       lon = AVHRR%lon-360.
    ELSEWHERE
       lon = AVHRR%lon
    ENDWHERE
    CALL Write_GBCS_Float(ncid,'lon',dimid_nx,dimid_ny,nx,ny,lon,&
         -32768.,'Longitude coordinates','longitude','degrees_east',&
         -180.,180.,'geographical coordinates, WGS84 projection')
    DEALLOCATE(lon)

    !
    ! Get time in right units
    !
    CALL Get_DateTime(AVHRR,stime,time_secs,1981)

    !
    ! Write out time
    !
    CALL Write_GBCS_time(ncid,dimid_time,stime,'seconds since 1981-01-01 00:00:00')
    CALL Write_GBCS_dtime(ncid,dimid_nx,dimid_ny,dimid_time,nx,ny,stime,&
         time_secs,'scanline time difference from start time',&
         'seconds','lon lat',-32768.)

    !
    ! Now channels
    !
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch1',dimid_nx,dimid_ny,dimid_time,&
         nx,ny,AVHRR%new_array1,-32768_GbcsInt2,'Channel 1 Reflectance','',&
         'reflectance',0_GbcsInt2,15000_GbcsInt2,&
         'lon lat',0.0001,0.,'')
    
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch2',dimid_nx,dimid_ny,dimid_time,&
         nx,ny,AVHRR%new_array2,-32768_GbcsInt2,'Channel 2 Reflectance','',&
         'reflectance',0_GbcsInt2,15000_GbcsInt2,&
         'lon lat',0.0001,0.,'')
    
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch3b',dimid_nx,dimid_ny,dimid_time,&
         nx,ny,FCDR%btf3,-32768_GbcsInt2,&
         'Channel 3b Brightness Temperature','',&
         'kelvin',-20000_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.01,273.15,'')

    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch4',dimid_nx,dimid_ny,dimid_time,&
         nx,ny,FCDR%btf4,-32768_GbcsInt2,'Channel 4 Brightness Temperature','',&
         'kelvin',-20000_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.01,273.15,'')
    
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch5',dimid_nx,dimid_ny,dimid_time,&
         nx,ny,FCDR%btf5,-32768_GbcsInt2,'Channel 5 Brightness Temperature','',&
         'kelvin',-20000_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.01,273.15,'')
    
    !
    ! Now noise estimates - not here we combine independent and structured
    ! from FIDUCEO information
    !
    noise = NAN_R
    WHERE(AVHRR%new_array1 .gt. 0. .and. AVHRR%new_array1_error .gt. 0.)
       noise = SQRT((CSPP_Uncertainty(1)*AVHRR%new_array1)**2+&
            AVHRR%new_array1_error**2)
    ENDWHERE
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch1_noise',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,noise,-32768_GbcsInt2,'Channel 1 noise estimate','',&
         'reflectance',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',1.e-5,0.,'')
    
    noise = NAN_R
    WHERE(AVHRR%new_array2 .gt. 0. .and. AVHRR%new_array2_error .gt. 0.)
       noise = SQRT((CSPP_Uncertainty(2)*AVHRR%new_array2)**2+&
            AVHRR%new_array2_error**2)
    ENDWHERE
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch2_noise',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,noise,-32768_GbcsInt2,'Channel 2 noise estimate','',&
         'reflectance',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',1.e-5,0.,'')
    
    noise = NAN_R
    WHERE( FCDR%ur3 .gt. 0. .and. FCDR%us3 .gt. 0 )
       noise = SQRT(FCDR%ur3**2 + FCDR%us3**2 + FCDR%uc3**2)
    ENDWHERE
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch3b_nedt',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,noise,-32768_GbcsInt2,'Channel 3b noise estimate','',&
         'kelvin',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.001,0.,'')
    
    noise = NAN_R
    WHERE( FCDR%ur4 .gt. 0. .and. FCDR%us4 .gt. 0 )
       noise = SQRT(FCDR%ur4**2 + FCDR%us4**2 + FCDR%uc4**2)
    ENDWHERE
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch4_nedt',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,noise,-32768_GbcsInt2,'Channel 4 noise estimate','',&
         'kelvin',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.001,0.,'')
    
    noise = NAN_R
    WHERE( FCDR%ur5 .gt. 0. .and. FCDR%us5 .gt. 0 )
       noise = SQRT(FCDR%ur5**2 + FCDR%us5**2 + FCDR%uc5**2)
    ENDWHERE
    CALL Write_GBCS_Float_to_Int_3D(ncid,'ch5_nedt',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,noise,-32768_GbcsInt2,'Channel 5 noise estimate','',&
         'kelvin',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.001,0.,'')
    !
    ! Angles
    !
    CALL Write_GBCS_Float_to_Int_3D(ncid,'satellite_zenith_angle',&
         dimid_nx,dimid_ny,dimid_time,nx,ny,AVHRR%satza,-32768_GbcsInt2,&
         'satellite zenith angle','zenith angle',&
         'angular_degree',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.01,0.,&
         'The satellite zenith angle at time of the observations')
    
    CALL Write_GBCS_Float_to_Int_3D(ncid,'solar_zenith_angle',&
         dimid_nx,dimid_ny,dimid_time,nx,ny,AVHRR%solza,-32768_GbcsInt2,&
         'solar zenith angle','zenith angle',&
         'angular_degree',0_GbcsInt2,18000_GbcsInt2,&
         'lon lat',0.01,0.,&
         'The solar zenith angle at time of the observations')
    
    CALL Write_GBCS_Float_to_Int_3D(ncid,'relative_azimuth_angle',&
         dimid_nx,dimid_ny,dimid_time,nx,ny,AVHRR%relaz,-32768_GbcsInt2,&
         'relative azimuth angle','zenith angle',&
         'angular_degree',0_GbcsInt2,10000_GbcsInt2,&
         'lon lat',0.01,0.,&
         'The relative azimuth angle at time of the observations')
    
    !
    ! ICT temperature
    !
    CALL Write_GBCS_Float_Time(ncid,'ict_temp',dimid_ny,dimid_time,&
         ny,AVHRR%smoothPRT,-32768_GbcsInt2,&
         'Temperature of internal calibration target',&
         'kelvin',-20000_GbcsInt2,10000_GbcsInt2,0.01,273.15)

    !
    ! Quality flags
    !
    CALL Make_Flags(AVHRR,nflags,flags,flags_meaning,flag_values,&
         valid_min,valid_max)
    CALL Write_GBCS_Byte(ncid,'qual_flags',dimid_nx,dimid_ny,dimid_time,nx,ny,&
       flags,-128_GbcsInt1,'Quality Flags',flags_meaning,nflags,flag_values,&
       valid_min,valid_max,'Bad data is flagged when any one of the other &
         &categories (bad navigation, calibration, timing or missing line) &
         &are flagged')

    !
    ! Cloud masks
    !
    flags = 0
    flag_values(1) = 0_GbcsInt2
    CALL Write_GBCS_Byte(ncid,'cloud_mask',dimid_nx,dimid_ny,dimid_time,nx,ny,&
         flags,-128_GbcsInt1,'No Cloud Mask','None',1,flag_values,&
         0_GbcsInt1,0_GbcsInt1,'')
    CALL Write_GBCS_Byte(ncid,'cloud_probability',dimid_nx,dimid_ny,&
         dimid_time,nx,ny,flags,-128_GbcsInt1,'No Cloud Mask',&
         'None',1,flag_values,0_GbcsInt1,0_GbcsInt1,'',scale=0.,offset=0.)
    DEALLOCATE(flags,flag_values)

    !
    ! Line numbers
    !
    CALL Write_GBCS_Float_Time_Int(ncid,'l1b_line_number',dimid_ny,dimid_time,&
         ny,AVHRR%scanLineNumber,-32768_GbcsInt2,&
         'FIDUCEO Level 1 line number',0_GbcsInt2)

    !
    ! If calibration data wanted, then output it here
    !
    IF( output_cal )THEN
       !
       ! Write out - smoothPRT, prt_correction, ICT counts, Space counts
       ! Earth counts, nuc, aval, bval all at full resolution
       !
       CALL Write_GBCS_Float_1d(ncid,'smoothprt',dimid_ny,ny,&
            AVHRR%smoothPRT,NAN_R,'Smoothed average PRT temperature',&
            'ICT_T','Kelvin',0.,400.)
       
       CALL Write_GBCS_Float_1d(ncid,'smoothprt1',dimid_ny,ny,&
            AVHRR%smoothPRT1,NAN_R,'Smoothed PRT1 temperature',&
            'ICT_T','Kelvin',0.,400.)
       
       CALL Write_GBCS_Float_1d(ncid,'smoothprt2',dimid_ny,ny,&
            AVHRR%smoothPRT2,NAN_R,'Smoothed PRT2 temperature',&
            'ICT_T','Kelvin',0.,400.)
       
       CALL Write_GBCS_Float_1d(ncid,'smoothprt3',dimid_ny,ny,&
            AVHRR%smoothPRT3,NAN_R,'Smoothed PRT3 temperature',&
            'ICT_T','Kelvin',0.,400.)
       
       CALL Write_GBCS_Float_1d(ncid,'smoothprt4',dimid_ny,ny,&
            AVHRR%smoothPRT4,NAN_R,'Smoothed PRT4 temperature',&
            'ICT_T','Kelvin',0.,400.)
       
       CALL Write_GBCS_Float_1d(ncid,'prt_correction',dimid_ny,ny,&
            AVHRR%prt_correction,NAN_R,'Correction to ICT temperature',&
            'ICT_T_Corr','Kelvin',-100.,100.)
       
       CALL Write_GBCS_Float_1d(ncid,'bb3',dimid_ny,ny,&
            AVHRR%bb3,NAN_R,'Blackbody counts 3.7 micron',&
            'BB3','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float_1d(ncid,'bb4',dimid_ny,ny,&
            AVHRR%bb4,NAN_R,'Blackbody counts 11 micron',&
            'BB4','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float_1d(ncid,'bb5',dimid_ny,ny,&
            AVHRR%bb5,NAN_R,'Blackbody counts 12 micron',&
            'BB5','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float_1d(ncid,'sp3',dimid_ny,ny,&
            AVHRR%sp3,NAN_R,'Space counts 3.7 micron',&
            'Sp3','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float_1d(ncid,'sp4',dimid_ny,ny,&
            AVHRR%sp4,NAN_R,'Space counts 11 micron',&
            'Sp4','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float_1d(ncid,'sp5',dimid_ny,ny,&
            AVHRR%sp5,NAN_R,'Space counts 12 micron',&
            'Sp5','Counts',-10.,1000.)
              
       CALL Write_GBCS_Float(ncid,'counts3',dimid_nx,dimid_ny,nx,ny,&
            AVHRR%counts3,NAN_R,'Earth Counts (3.7 mu)',&
            'CE_3','Counts',&
            -10.,1000.,'')

       CALL Write_GBCS_Float(ncid,'counts4',dimid_nx,dimid_ny,nx,ny,&
            AVHRR%counts4,NAN_R,'Earth Counts (11 mu)',&
            'CE_4','Counts',&
            -10.,1000.,'')

       CALL Write_GBCS_Float(ncid,'counts5',dimid_nx,dimid_ny,nx,ny,&
            AVHRR%counts5,NAN_R,'Earth Counts (12 mu)',&
            'CE_5','Counts',&
            -10.,1000.,'')

       !
       ! Band coefficients
       !
       stat = NF90_REDEF(ncid)
       call check(stat)

       dims1 = dimid_band
       stat = NF90_DEF_VAR(ncid,'nuc',NF90_FLOAT,dims1,varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
       call check(stat)
       stat = NF90_PUT_VAR(ncid,varid,AVHRR%nuc)
       CALL check(stat)

       stat = NF90_DEF_VAR(ncid,'aval',NF90_FLOAT,dims1,varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
       call check(stat)
       stat = NF90_PUT_VAR(ncid,varid,AVHRR%aval)
       CALL check(stat)

       stat = NF90_DEF_VAR(ncid,'bval',NF90_FLOAT,dims1,varid)
       call check(stat)
       stat = NF90_DEF_VAR_DEFLATE(ncid, varid, 1, 1, 9)
       call check(stat)
       stat = NF90_PUT_VAR(ncid,varid,AVHRR%bval)
       CALL check(stat)

    ENDIF
    
    !
    ! Write Global attributes
    !
    CALL Write_GBCS_L1C_Global(AVHRR,ncid,uuid,noaa_string)

    !
    ! Close file
    !
    stat = NF90_CLOSE(ncid)
    call check(stat)

  END SUBROUTINE Write_GBCS_L1C

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
