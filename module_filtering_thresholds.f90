MODULE module_filtering_thresholds
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil
  USE type

implicit none
PUBLIC

CONTAINS

   SUBROUTINE filtering_sp_bb(scanLinePos,outData,DataForFiltering)

   
    INTEGER, INTENT(IN)                          :: scanLinePos
    TYPE(AVHRR_Data), INTENT (INOUT)             :: outData
    TYPE(AVHRR_DAta)                             :: DataForFiltering 

    TYPE(AVHRR_thresholds_coarse)                :: Threshold
   
    INTEGER                                      ::i,j,z=3,chanel,indice_chanel 

    INTEGER, DIMENSION(3)                        :: n_allan_sp,n_allan_bb, &
                                                    n_points_for_mean_sp, &
                                                    n_points_for_mean_bb, &
                                                    n_mean_bb, n_mean_sp
                                                   
    DOUBLE PRECISION, DIMENSION(3)               :: allan_variance_sp, &
                                                    allan_variance_bb, &
                                                    mean_sp, &
                                                    mean_bb, &
                                                    somme_sp, &
                                                    somme_bb, &
                                                    somme_for_mean_sp, &
                                                    somme_for_mean_bb, &
                                                    thresh_low_sp, &
                                                    thresh_high_sp, &
                                                    thresh_low_bb, &
                                                    thresh_high_bb
  !print*, "************* FILTERING***********"
  !print*, "scanLine",scanLinePos
  !print*, DataForFiltering%spaceFilter3(:,scanLinePos)  
  DO j=1,10
    DataForFiltering%spaceFilter3(j,scanLinePos)=outData%spaceFilter3(j,scanLinePos)
    DataForFiltering%spaceFilter4(j,scanLinePos)=outData%spaceFilter4(j,scanLinePos)
    DataForFiltering%spaceFilter5(j,scanLinePos)=outData%spaceFilter5(j,scanLinePos)
    DataForFiltering%BbodyFilter3(j,scanLinePos)=outData%BbodyFilter3(j,scanLinePos)
    DataForFiltering%BbodyFilter4(j,scanLinePos)=outData%BbodyFilter4(j,scanLinePos)
    DataForFiltering%BbodyFilter5(j,scanLinePos)=outData%BbodyFilter5(j,scanLinePos)
  END DO
         
  DataForFiltering%bb3(scanLinePos) = outData%bb3(scanLinePos)
  DataForFiltering%bb4(scanLinePos) = outData%bb4(scanLinePos)
  DataForFiltering%bb5(scanLinePos) = outData%bb5(scanLinePos)
  DataForFiltering%sp3(scanLinePos) = outData%sp3(scanLinePos)
  DataForFiltering%sp4(scanLinePos) = outData%sp4(scanLinePos)
  DataForFiltering%sp5(scanLinePos) = outData%sp5(scanLinePos)

  SELECT CASE( outData%AVHRR_No )
    CASE(1)
      !print*,'TIROS-N'
      Threshold%Low_sp=(/900,990,980/)
      Threshold%High_sp=(/1025,1005,1005/)
      Threshold%Low_bb=(/400,340,400/)
      Threshold%High_bb=(/850,550,510/)
    
    CASE(6)
      !print*,'NOAA-6'
      Threshold%Low_sp=(/950,980,985/)
      Threshold%High_sp=(/1025,1005,1000/)
      Threshold%Low_bb=(/500,300,300/)
      Threshold%High_bb=(/830,520,520/)
   
    CASE(7)
      !print*,'NOAA-7'
      Threshold%Low_sp=(/900,985,985/)
      Threshold%High_sp=(/1025,1000,1000/)
      Threshold%Low_bb=(/600,200,300/)
      Threshold%High_bb=(/930,420,520/)
    CASE(8)
      !print*,'NOAA-8'
      Threshold%Low_sp=(/940,980,980/)
      Threshold%High_sp=(/1025,1005,1005/)
      Threshold%Low_bb=(/600,300,300/)
      Threshold%High_bb=(/850,420,420/)

    CASE(9)
      !print*, 'NOAA-9'
      Threshold%Low_sp=(/960,985,985/)
      Threshold%High_sp=(/1025,1000,1000/)
      Threshold%Low_bb=(/500,200,300/)
      Threshold%High_bb=(/850,450,500/)
  
    CASE(10)
      !print*,'NOAA-10'
      Threshold%Low_sp=(/900,985,985/)
      Threshold%High_sp=(/1100,1015,1015/)
      Threshold%Low_bb=(/400,200,200/)
      Threshold%High_bb=(/950,450,450/)

    CASE(11)
      !print*,'NOAA-11'
      Threshold%Low_sp=(/960,980,980/)
      Threshold%High_sp=(/1025,1015,1015/)
      Threshold%Low_bb=(/500,300,200/)
      Threshold%High_bb=(/850,500,500/)
 
    CASE(12)
      !print*, 'NOAA-12'
      Threshold%Low_sp=(/950,980,980/)
      Threshold%High_sp=(/1025,1000,1015/)
      Threshold%Low_bb=(/400,200,200/)
      Threshold%High_bb=(/850,450,450/)

    CASE(14)
      !print*, 'NOAA-14'
      Threshold%Low_sp=(/970,980,980/)
      Threshold%High_sp=(/1025,1025,1020/)
      Threshold%Low_bb=(/600,300,300/)
      Threshold%High_bb=(/850,500,500/)

    CASE(15)
      !print*,'NOAA-15'
      Threshold%Low_sp=(/920,920,920/)
      Threshold%High_sp=(/1030,1030,1030/)
      Threshold%Low_bb=(/700,300,300/)
      Threshold%High_bb=(/950,530,530/)
    CASE(16)
      !print*,'NOAA-16'
      Threshold%Low_sp=(/900,980,980/)
      Threshold%High_sp=(/1100,1010,1000/)
      Threshold%Low_bb=(/700,300,350/)
      Threshold%High_bb=(/950,500,550/)
    CASE(17)
      !print*,'NOAA-17'
      Threshold%Low_sp=(/900,940,970/)
      Threshold%High_sp=(/1020,1025,1020/)
      Threshold%Low_bb=(/0,400,350/)
      Threshold%High_bb=(/1000,600,550/)
    CASE(18)
      !print*,'NOAA-18'
      Threshold%Low_sp=(/960,985,985/)
      Threshold%High_sp=(/1050,1000,1000/)
      Threshold%Low_bb=(/800,420,320/)
      Threshold%High_bb=(/950,550,550/)
    CASE(19)
      !print*,'NOAA-19'
      Threshold%Low_sp=(/900,985,985/)
      Threshold%High_sp=(/1030,1000,1000/)
      Threshold%Low_bb=(/780,430,430/)
      Threshold%High_bb=(/900,540,530/)
    CASE(-1)
      !print*,'Metop-02' 
       Threshold%Low_sp=(/900,985,985/)
       Threshold%High_sp=(/1030,1000,1000/)
       Threshold%Low_bb=(/780,430,450/)
       Threshold%High_bb=(/900,550,530/)

    CASE(-2) 
      !print*, "Metop-01"

  END SELECT
 
  
  !print*, "****** Premier Tri ******"
  !print*, "sp counts 3 ",OutData%spaceFilter3(:,scanLinePos)
  !print*, "sp counts 4 ",OutData%spaceFilter4(:,scanLinePos)
  !print*, "sp counts 5 ",OutData%spaceFilter5(:,scanLinePos)

  !print*, "Threholds sp 3", Threshold%Low_sp(1),Threshold%High_sp(1)
  !print*, "Threholds sp 4", Threshold%Low_sp(2),Threshold%High_sp(2)
  !print*, "Threholds sp 5", Threshold%Low_sp(3),Threshold%High_sp(3)

  !print*, "bb counts 3 ",OutData%bb3(scanLinePos),OutData%BbodyFilter3(:,scanLinePos)
  !print*, "bb counts 4 ",OutData%BbodyFilter4(:,scanLinePos)
  !print*, "bb counts 5 ",OutData%BbodyFilter5(:,scanLinePos)

  !print*, "Threholds bb 3", Threshold%Low_bb(1),Threshold%High_bb(1)
  !print*, "Threholds bb 4", Threshold%Low_bb(2),Threshold%High_bb(2)
  !print*, "Threholds bb 5", Threshold%Low_bb(3),Threshold%High_bb(3)

  Do j=1,10
    !print*, "***** Space counts *******"
    IF (( DataForFiltering%spaceFilter3(j,scanLinePos) .lt. Threshold%Low_sp(1) ) &
     .or. ( DataForFiltering%spaceFilter3(j,scanLinePos) .gt. Threshold%High_sp(1)) )  THEN
      DataForFiltering%spaceFilter3(j,scanLinePos)=-999                 
    ELSE
      DataForFiltering%spaceFilter3(j,scanLinePos)=outData%spaceFilter3(j,scanLinePos)
    END IF

    IF (( DataForFiltering%spaceFilter4(j,scanLinePos) .lt. Threshold%Low_sp(2) ) &
     .or. ( DataForFiltering%spaceFilter4(j,scanLinePos) .gt. Threshold%High_sp(2)) )  THEN
      DataForFiltering%spaceFilter4(j,scanLinePos)=-999               
    ELSE
      DataForFiltering%spaceFilter4(j,scanLinePos)=outData%spaceFilter4(j,scanLinePos)
    END IF

    IF (( DataForFiltering%spaceFilter5(j,scanLinePos) .lt. Threshold%Low_sp(3) ) &
     .or. ( DataForFiltering%spaceFilter5(j,scanLinePos) .gt. Threshold%High_sp(3)) )  THEN
      DataForFiltering%spaceFilter5(j,scanLinePos)=-999              
    ELSE
      DataForFiltering%spaceFilter5(j,scanLinePos)=outData%spaceFilter5(j,scanLinePos)
    END IF

    !print*, "***** ICT counts *******"
    IF (( DataForFiltering%BbodyFilter3(j,scanLinePos) .lt. Threshold%Low_bb(1) ) &
     .or. ( DataForFiltering%BbodyFilter3(j,scanLinePos) .gt. Threshold%High_bb(1)) )  THEN
      DataForFiltering%BbodyFilter3(j,scanLinePos)=-999              
    ELSE
      DataForFiltering%bbodyFilter3(j,scanLinePos)=outData%BbodyFilter3(j,scanLinePos)
    END IF

    IF (( DataForFiltering%BbodyFilter4(j,scanLinePos) .lt. Threshold%Low_bb(2) ) &
     .or. ( DataForFiltering%BbodyFilter4(j,scanLinePos) .gt. Threshold%High_bb(2)) )  THEN
      DataForFiltering%BbodyFilter4(j,scanLinePos)=-999                 
    ELSE
      DataForFiltering%BbodyFilter4(j,scanLinePos)=outData%BbodyFilter4(j,scanLinePos)
    END IF

    IF (( DataForFiltering%BbodyFilter5(j,scanLinePos) .lt. Threshold%Low_bb(3) ) &
     .or. ( DataForFiltering%BbodyFilter5(j,scanLinePos) .gt. Threshold%High_bb(3)) )  THEN
      DataForFiltering%BbodyFilter5(j,scanLinePos)=-999                
    ELSE
      DataForFiltering%BbodyFilter5(j,scanLinePos)=outData%BbodyFilter5(j,scanLinePos)
    END IF

  END DO
  !print*, "sp counts filtres 3",DataForFiltering%spaceFilter3(:,scanLinePos)
  !print*, "sp counts filtres 4",DataForFiltering%spaceFilter4(:,scanLinePos)
  !print*, "sp counts filtres 5",DataForFiltering%spaceFilter5(:,scanLinePos)

  !print*, "bb counts filtres 3",DataForFiltering%BbodyFilter3(:,scanLinePos)
  !print*, "bb counts filtres 4",DataForFiltering%BbodyFilter4(:,scanLinePos)
  !print*, "bb counts filtres 5",DataForFiltering%BbodyFilter5(:,scanLinePos)

!print*, " ****** Deuxieme tri  *******"
!-On calcule l'Allan deviation par paquets de 25 lignes (fenetre glissante)
!--Initialisation des tableaux
   n_points_for_mean_sp=(/0,0,0/)
   n_points_for_mean_bb=(/0,0,0/)
   n_allan_sp=(/0,0,0/)
   n_allan_bb=(/0,0,0/)
  
   somme_for_mean_sp=(/0,0,0/)
   somme_for_mean_bb=(/0,0,0/)
   somme_sp=(/0,0,0/)
   somme_bb=(/0,0,0/)

   allan_variance_sp=(/NAN_R,NAN_R,NAN_R/)
   allan_variance_bb=(/NAN_R,NAN_R,NAN_R/)
   
   mean_sp=(/NAN_R,NAN_R,NAN_R/)
   mean_bb=(/NAN_R,NAN_R,NAN_R/)
   
   thresh_low_sp=(/NAN_R,NAN_R,NAN_R/)
   thresh_high_sp=(/NAN_R,NAN_R,NAN_R/)
   thresh_low_bb=(/NAN_R,NAN_R,NAN_R/)
   thresh_high_bb=(/NAN_R,NAN_R,NAN_R/)
                                                 
   if ((scanLinePos>12) .and. (scanlinePos<outdata%arraysize-12)) then
!--On fait une boucle sur les n lignes qui entourent la ligne que l'on regarde
    DO i=scanLinePos-12,scanLinePos+12
!-----Boucle sur les 10 valeurs de chaque ligne     
      DO j=1,10
        IF ((DataForFiltering%spaceFilter3(j,i) .ne. NAN_I) .and. (DataForFiltering%spaceFilter3(j,i) .gt. 0))then
          somme_for_mean_sp(1)=somme_for_mean_sp(1)+DataForFiltering%spaceFilter3(j,i)
          n_points_for_mean_sp(1)=n_points_for_mean_sp(1)+1
        END IF 
  
        IF ((DataForFiltering%spaceFilter4(j,i) .ne. NAN_I) .and. (DataForFiltering%spaceFilter4(j,i) .gt. 0)) then
          somme_for_mean_sp(2)=somme_for_mean_sp(2)+DataForFiltering%spaceFilter4(j,i)
          n_points_for_mean_sp(2)=n_points_for_mean_sp(2)+1
        END IF 
  
        IF ((DataForFiltering%spaceFilter5(j,i) .ne. NAN_I) .and. (DataForFiltering%spaceFilter5(j,i) .gt. 0)) then
          somme_for_mean_sp(3)=somme_for_mean_sp(3)+DataForFiltering%spaceFilter5(j,i)
          n_points_for_mean_sp(3)=n_points_for_mean_sp(3)+1
        END IF 

        IF ((DataForFiltering%BbodyFilter3(j,i) .ne. NAN_I) .and. (DataForFiltering%BbodyFilter3(j,i) .gt. 0))then
          somme_for_mean_bb(1)=somme_for_mean_bb(1)+DataForFiltering%BbodyFilter3(j,i)
          n_points_for_mean_bb(1)=n_points_for_mean_bb(1)+1
        END IF

        IF ((DataForFiltering%BbodyFilter4(j,i) .ne. NAN_I) .and. (DataForFiltering%BbodyFilter4(j,i) .gt. 0))then
          somme_for_mean_bb(2)=somme_for_mean_bb(2)+DataForFiltering%BbodyFilter4(j,i)
          n_points_for_mean_bb(2)=n_points_for_mean_bb(2)+1
        END IF

        IF ((DataForFiltering%BbodyFilter5(j,i) .ne. NAN_I) .and. (DataForFiltering%BbodyFilter5(j,i) .gt. 0))then
          somme_for_mean_bb(3)=somme_for_mean_bb(3)+DataForFiltering%BbodyFilter5(j,i)
          n_points_for_mean_bb(3)=n_points_for_mean_bb(3)+1
        END IF 
      End DO !j=1,10

      DO j=1,9
         !print*, "j",j
        IF ((DataForFiltering%spaceFilter3(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter3(j+1,i) .ne. -999) &
      .and. (DataForFiltering%spaceFilter3(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter3(j,i) .ne. -999)) THEN
          somme_sp(1)=somme_sp(1)+(DataForFiltering%spaceFilter3(j+1,i)-DataForFiltering%spaceFilter3(j,i))**2
          n_allan_sp(1)=n_allan_sp(1)+1
        END IF

        IF ((DataForFiltering%spaceFilter4(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter4(j+1,i) .ne. -999) &
      .and. (DataForFiltering%spaceFilter4(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter4(j,i) .ne. -999)) THEN
          somme_sp(2)=somme_sp(2)+(DataForFiltering%spaceFilter4(j+1,i)-DataForFiltering%spaceFilter4(j,i))**2
          n_allan_sp(2)=n_allan_sp(2)+1
        END IF

         IF ((DataForFiltering%spaceFilter5(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter5(j+1,i) .ne. -999) &
      .and. (DataForFiltering%spaceFilter5(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%spaceFilter5(j,i) .ne. -999)) THEN
          somme_sp(3)=somme_sp(3)+(DataForFiltering%spaceFilter5(j+1,i)-DataForFiltering%spaceFilter5(j,i))**2
          n_allan_sp(3)=n_allan_sp(3)+1
        END IF

        IF ((DataForFiltering%BbodyFilter3(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter3(j+1,i) .ne. -999) &
      .and. (DataForFiltering%BbodyFilter3(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter3(j,i) .ne. -999)) THEN
          somme_bb(1)=somme_bb(1)+(DataForFiltering%BbodyFilter3(j+1,i)-DataForFiltering%BbodyFilter3(j,i))**2
          n_allan_bb(1)=n_allan_bb(1)+1
        END IF

         IF ((DataForFiltering%BbodyFilter4(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter4(j+1,i) .ne. -999) &
      .and. (DataForFiltering%BbodyFilter4(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter4(j,i) .ne. -999)) THEN
          somme_bb(2)=somme_bb(2)+(DataForFiltering%BbodyFilter4(j+1,i)-DataForFiltering%BbodyFilter4(j,i))**2
          n_allan_bb(2)=n_allan_bb(2)+1
          !print*, DataForFiltering%BbodyFilter4(j+1,i),DataForFiltering%BbodyFilter4(j,i),somme_bb(2), n_allan_bb(2)
        END IF

        IF ((DataForFiltering%BbodyFilter5(j+1,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter5(j+1,i) .ne. -999) &
      .and. (DataForFiltering%BbodyFilter5(j,i) .ne. NAN_I) &
      .and. (DataForFiltering%BbodyFilter5(j,i) .ne. -999)) THEN
          somme_bb(3)=somme_bb(3)+(DataForFiltering%BbodyFilter5(j+1,i)-DataForFiltering%BbodyFilter5(j,i))**2
          n_allan_bb(3)=n_allan_bb(3)+1
        END IF
      END DO !j=1,9
    END DO !Boucle sur les 25 scanlines

    Do chanel=3,5
      indice_chanel=chanel-2
!-----On en deduit les Allan deviation    
      IF (n_points_for_mean_sp(indice_chanel) .gt. 0) THEN
        mean_sp(indice_chanel)=somme_for_mean_sp(indice_chanel)/n_points_for_mean_sp(indice_chanel)
      END IF
      IF (somme_sp(indice_chanel) .eq. 0) THEN
        allan_variance_sp(indice_chanel)= 0 
      ELSE IF (somme_sp(indice_chanel) .gt. 0) THEN
        !allan_variance_sp(indice_chanel)= sqrt(somme_sp(indice_chanel)/(2.*n_allan_sp(indice_chanel))) 
        allan_variance_sp(indice_chanel)= sqrt(somme_sp(indice_chanel)/(2.*9)) 
      ELSE
        allan_variance_sp(indice_chanel)=NAN_R
      END IF
      
      IF (n_points_for_mean_bb(indice_chanel) .gt. 0) THEN
        mean_bb(indice_chanel)=somme_for_mean_bb(indice_chanel)/n_points_for_mean_bb(indice_chanel)
      END IF
      IF (somme_bb(indice_chanel) .eq. 0) THEN
        allan_variance_bb(indice_chanel)= 0 
      ELSE IF (somme_bb(indice_chanel) .gt. 0) THEN
        !allan_variance_bb(indice_chanel)= sqrt(somme_bb(indice_chanel)/(2.*n_allan_bb(indice_chanel)))
        allan_variance_bb(indice_chanel)= sqrt(somme_bb(indice_chanel)/(2.*18))
      ELSE
        allan_variance_bb(indice_chanel)=NAN_R
      END IF
!-----On calcule les seuils
      IF (allan_variance_sp(indice_chanel) .ne. NAN_R) THEN
         thresh_high_sp(indice_chanel)=mean_sp(indice_chanel)+z*allan_variance_sp(indice_chanel)
         thresh_low_sp(indice_chanel)=mean_sp(indice_chanel)-z*allan_variance_sp(indice_chanel)
      END IF
      IF (allan_variance_bb(indice_chanel) .ne. NAN_R) THEN
         thresh_high_bb(indice_chanel)=mean_bb(indice_chanel)+z*allan_variance_bb(indice_chanel)
         thresh_low_bb(indice_chanel)=mean_bb(indice_chanel)-z*allan_variance_bb(indice_chanel)
      END IF
    ENd DO !chanel
   

!---On trie les valeurs  et on calcule les moyennes par scanline  
    n_mean_bb=(/0,0,0/)
    n_mean_sp=(/0,0,0/)
   
    DataForFiltering%sp3(scanLinePos)=0
    DataForFiltering%sp4(scanLinePos)=0
    DataForFiltering%sp5(scanLinePos)=0
    DataForFiltering%bb3(scanLinePos)=0
    DataForFiltering%bb4(scanLinePos)=0
    DataForFiltering%bb5(scanLinePos)=0
    DO j=1,10
!*******Space counts********
!******* ch 3******
        IF (allan_variance_sp(1) .ne. NAN_R) THEN
           IF ( (DataForFiltering%spaceFilter3(j,scanLinePos) .lt. thresh_low_sp(1)) &
           .or. (DataForFiltering%spaceFilter3(j,scanLinePos) .gt. thresh_high_sp(1))) THEN
             DataForFiltering%spaceFilter3(j,scanLinePos)=-999
           ELSE IF (DataForFiltering%spaceFilter3(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%sp3(scanLinePos)=DataForFiltering%sp3(scanLinePos)&
                  +DataForFiltering%spaceFilter3(j,scanLinePos)
             n_mean_sp(1)=n_mean_sp(1)+1
           END IF
        ELSE
           DataForFiltering%spaceFilter3(j,scanLinePos)=-999
        END IF
!******* ch 4******
        IF (allan_variance_sp(2) .ne. NAN_R) THEN
           IF ( (DataForFiltering%spaceFilter4(j,scanLinePos) .lt. thresh_low_sp(2)) &
           .or. (DataForFiltering%spaceFilter4(j,scanLinePos) .gt. thresh_high_sp(2))) THEN
             DataForFiltering%spaceFilter4(j,scanLinePos)=-999
           ELSE IF (DataForFiltering%spaceFilter4(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%sp4(scanLinePos)=DataForFiltering%sp4(scanLinePos)&
                  +DataForFiltering%spaceFilter4(j,scanLinePos)
             n_mean_sp(2)=n_mean_sp(2)+1
           END IF
        ELSE
           DataForFiltering%spaceFilter4(j,scanLinePos)=-999
        END IF
!******* ch 5******
        IF (allan_variance_sp(3) .ne. NAN_R) THEN
           IF ( (DataForFiltering%spaceFilter5(j,scanLinePos) .lt. thresh_low_sp(3)) &
           .or. (DataForFiltering%spaceFilter5(j,scanLinePos) .gt. thresh_high_sp(3))) THEN
             DataForFiltering%spaceFilter5(j,scanLinePos)=-999
           ELSE IF (DataForFiltering%spaceFilter5(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%sp5(scanLinePos)=DataForFiltering%sp5(scanLinePos)&
                  +DataForFiltering%spaceFilter5(j,scanLinePos)
             n_mean_sp(3)=n_mean_sp(3)+1
           END IF
        ELSE
           DataForFiltering%spaceFilter5(j,scanLinePos)=-999
        END IF

!*******ICT ***************
        IF (allan_variance_bb(1) .ne. NAN_R) THEN
!**********ch 3******
          IF ( (DataForFiltering%bbodyFilter3(j,scanLinePos) .lt. thresh_low_bb(1)) &
          .or. (DataForFiltering%BbodyFilter3(j,scanLinePos) .gt. thresh_high_bb(1))) THEN
            DataForFiltering%BbodyFilter3(j,scanLinePos)=-999
          ELSE IF (DataForFiltering%BbodyFilter3(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%bb3(scanLinePos)=DataForFiltering%bb3(scanLinePos)&
                  +DataForFiltering%BbodyFilter3(j,scanLinePos)
             n_mean_bb(1)=n_mean_bb(1)+1
          END IF
        ELSE
           DataForFiltering%BbodyFilter3(j,scanLinePos)=-999
        END IF

!******* ch 4******
        IF (allan_variance_bb(2) .ne. NAN_R) THEN
          IF ( (DataForFiltering%bbodyFilter4(j,scanLinePos) .lt. thresh_low_bb(2)) &
          .or. (DataForFiltering%BbodyFilter4(j,scanLinePos) .gt. thresh_high_bb(2))) THEN
            DataForFiltering%BbodyFilter4(j,scanLinePos)=-999
          ELSE IF (DataForFiltering%BbodyFilter4(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%bb4(scanLinePos)=DataForFiltering%bb4(scanLinePos)&
                  +DataForFiltering%BbodyFilter4(j,scanLinePos)
             n_mean_bb(2)=n_mean_bb(2)+1
          END IF
        ELSE
           DataForFiltering%BbodyFilter4(j,scanLinePos)=-999
        END IF

!******* ch 5******
       
        IF (allan_variance_bb(3) .ne. NAN_R) THEN
          IF ( (DataForFiltering%bbodyFilter5(j,scanLinePos) .lt. thresh_low_bb(3)) &
          .or. (DataForFiltering%BbodyFilter5(j,scanLinePos) .gt. thresh_high_bb(3))) THEN
            DataForFiltering%BbodyFilter5(j,scanLinePos)=-999
          ELSE  IF (DataForFiltering%BbodyFilter5(j,scanLinePos) .ne. NAN_I) THEN
             DataForFiltering%bb5(scanLinePos)=DataForFiltering%bb5(scanLinePos)&
                  +DataForFiltering%BbodyFilter5(j,scanLinePos)
             n_mean_bb(3)=n_mean_bb(3)+1
          END IF
        ELSE
           DataForFiltering%BbodyFilter5(j,scanLinePos)=-999
        END IF
    END DO ! boucle sur les 10 valeurs 
    if( n_mean_sp(1)< 8)then
       DataForFiltering%sp3(scanLinePos) = -999
    else 
       
       DataforFiltering%sp3(scanLinePos) = DataForFiltering%sp3(scanLinePos) / n_mean_sp(1)
    end if
    if( n_mean_sp(2)< 8)then
       DataForFiltering%sp4(scanLinePos) = -999
    else 
       DataforFiltering%sp4(scanLinePos) = DataForFiltering%sp4(scanLinePos) / n_mean_sp(2)
    end if
    if( n_mean_sp(3)< 8)then
       DataForFiltering%sp5(scanLinePos) = -999
    else 
       DataforFiltering%sp5(scanLinePos) = DataForFiltering%sp5(scanLinePos) / n_mean_sp(3)
    end if

    if( n_mean_bb(1)< 8)then
       DataForFiltering%bb3(scanLinePos) = -999
    else 
       !print*, DataForFiltering%bb3(scanLinePos),n_mean_bb(1)
       DataforFiltering%bb3(scanLinePos) = DataForFiltering%bb3(scanLinePos) / n_mean_bb(1)
       !print*, DataForFiltering%bb3(scanLinePos)
    end if
    if( n_mean_bb(2)< 8)then
       DataForFiltering%bb4(scanLinePos) = -999
    else 
       DataforFiltering%bb4(scanLinePos) = DataForFiltering%bb4(scanLinePos) / n_mean_bb(2)
    end if
    if( n_mean_bb(3)< 8)then
       DataForFiltering%bb5(scanLinePos) = -999
    else 
       DataforFiltering%bb5(scanLinePos) = DataForFiltering%bb5(scanLinePos) / float(n_mean_bb(3))
    end if
    
    !print*,"sp3", DataForFiltering%spaceFilter3(:,scanLinePos)
    !print*,"sp 4", DataForFiltering%spaceFilter4(:,scanLinePos)
    !print*,"sp 5", DataForFiltering%spaceFilter5(:,scanLinePos)
   
    !print*,"bb 3", DataForFiltering%BbodyFilter3(:,scanLinePos)
    !print*,"bb 5", DataForFiltering%BbodyFilter5(:,scanLinePos)
    !print*,"sp3m", DataForFiltering%sp3(scanLinePos)
    !print*,"sp m4", DataForFiltering%sp4(scanLinePos)
    !print*,"sp m5", DataForFiltering%sp5(scanLinePos)
    !print*,"bbm3", DataForFiltering%Bb3(scanLinePos)
    !print*,"bb m4", DataForFiltering%Bb4(scanLinePos)
    !print*,"bb m5", DataForFiltering%Bb5(scanLinePos)
   end if ! scanline pos
   
DO j=1,10
  outData%spaceFilter3(j,scanLinePos)=DataForFiltering%spaceFilter3(j,scanLinePos)
  outData%spaceFilter4(j,scanLinePos)=DataForFiltering%spaceFilter4(j,scanLinePos)
  outData%spaceFilter5(j,scanLinePos)=DataForFiltering%spaceFilter5(j,scanLinePos)
  outData%BbodyFilter3(j,scanLinePos)=DataForFiltering%BbodyFilter3(j,scanLinePos)
  outData%BbodyFilter4(j,scanLinePos)=DataForFiltering%BbodyFilter4(j,scanLinePos)
  outData%BbodyFilter5(j,scanLinePos)=DataForFiltering%BbodyFilter5(j,scanLinePos)
END DO
outData%sp3(scanLinePos)=DataForFiltering%sp3(scanLinePos)
outData%sp4(scanLinePos)=DataForFiltering%sp4(scanLinePos)
outData%sp5(scanLinePos)=DataForFiltering%sp5(scanLinePos)
outData%Bb3(scanLinePos)=DataForFiltering%Bb3(scanLinePos)
outData%Bb4(scanLinePos)=DataForFiltering%Bb4(scanLinePos)
outData%Bb5(scanLinePos)=DataForFiltering%Bb5(scanLinePos)
  
END SUBROUTINE filtering_sp_bb 
 

END MODULE module_filtering_thresholds

