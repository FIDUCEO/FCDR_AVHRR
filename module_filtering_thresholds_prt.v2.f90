MODULE module_filtering_thresholds_prt
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

   SUBROUTINE filtering_prt(scanLinePos,outData,DataForFiltering,prtcounts3,prtcounts4,prtcounts5)

   
    INTEGER, INTENT(IN)                          :: scanLinePos,prtcounts3,prtcounts4,prtcounts5
    INTEGER,DIMENSION(4)                         :: number_prt=(/2,3,4,5/) 
    TYPE(AVHRR_Data), INTENT (INOUT)             :: outData
    INTEGEr, DIMENSION(3)                        :: prtcounts

    TYPE(AVHRR_thresholds_coarse)                :: Threshold
   
    INTEGER                                      ::i,j,z=3,i_prt,i_lecture

    INTEGER, DIMENSION(4)                        :: n_allan_prt, &
                                                    n_points_for_mean_prt, &
                                                    n_mean_prt
                                                   
    DOUBLE PRECISION, DIMENSION(4)               :: allan_variance_prt, &
                                                    mean_prt, &
                                                    somme_prt, &
                                                    somme_for_mean_prt, &
                                                    thresh_low_prt, &
                                                    thresh_high_prt
  !print*, "************* FILTERING PRT***********"
  !print*, "scanLine",scanLinePos
  !print*, DataForFiltering%spaceFilter3(:,scanLinePos)  
 
  prtcounts(1)=prtcounts3
  prtcounts(2)=prtcounts4
  prtcounts(3)=prtcounts5
 
  SELECT CASE( outData%AVHRR_No )
    CASE(1)
      !print*,'TIROS-N'
      Threshold%low_prt=(/250,250,250,250/)
      Threshold%High_prt=(/450,450,450,450/)

    CASE(6)
      !print*,'NOAA-6'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/400,400,400,450/)

    CASE(7)
      !print*,'NOAA-7'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/350,400,400,400/)
    
    CASE(8)
      !print*,'NOAA-8'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/300,350,350,350/)

    CASE(9)
      !print*, 'NOAA-9'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/500,500,500,500/)
  
    CASE(10)
      !print*,'NOAA-10'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/500,500,500,500/)

    CASE(11)
      !print*,'NOAA-11'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/500,500,500,500/)
 
    CASE(12)
      !print*, 'NOAA-12'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/550,550,550,550/)

    CASE(14)
      !print*, 'NOAA-14'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/600,600,600,600/)

    CASE(15)
      !print*,'NOAA-15'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/500,500,500,500/)

    CASE(16)
      !print*,'NOAA-16'
      Threshold%Low_prt=(/200,150,150,150/)
      Threshold%High_prt=(/550,550,550,550/)

    CASE(17)
      !print*,'NOAA-17'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/450,450,450,450/)
    
    CASE(18)
      !print*,'NOAA-18'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/300,300,300,300/)

    CASE(19)
      !print*,'NOAA-19'
      Threshold%Low_prt=(/150,150,150,150/)
      Threshold%High_prt=(/250,250,250,300/)

    CASE(-1)
      !print*,'Metop-02' 
       Threshold%Low_prt=(/150,150,150,150/)
       Threshold%High_prt=(/200,200,250,250/)

    CASE(-2) 
      !print*, "Metop-01"

  END SELECT
 
  
  !print*, "****** Premier Tri ******"
  
  DO i_prt=1,4
    if (outdata%prtnumber(scanLinePos) .eq. number_prt(i_prt)) THEN
        DO i_lecture=1,3
           IF (( prtcounts(i_lecture) .lt. Threshold%Low_prt(i_prt) ) &
     .or. ( prtcounts(i_lecture) .gt. Threshold%High_prt(i_prt)) )  THEN
              prtcounts(i_lecture)=-999                 
           ELSE
              prtcounts(i_lecture)=prtcounts(i_lecture)
           END IF
        END DO
    ELSE
       cycle
    END IF
  END DO
 
!print*, " ****** Deuxieme tri  *******"
!-On calcule l'Allan deviation par paquets de 25 lignes (fenetre glissante)
!--Initialisation des tableaux
   n_points_for_mean_prt=(/0,0,0,0/)
   n_allan_prt=(/0,0,0,0/)
   somme_for_mean_prt=(/0,0,0,0/)
   somme_prt=(/0,0,0,0/)
   allan_variance_prt=(/NAN_R,NAN_R,NAN_R,NAN_R/)
   mean_prt=(/NAN_R,NAN_R,NAN_R,NAN_R/)
   thresh_low_prt=(/NAN_R,NAN_R,NAN_R,NAN_R/)
   thresh_high_prt=(/NAN_R,NAN_R,NAN_R,NAN_R/)
                                                 
   if ((scanLinePos>12) .and. (scanlinePos<outdata%arraysize-12)) then
!--On fait une boucle sur les n lignes qui entourent la ligne que l'on regarde
    DO i=scanLinePos-12,scanLinePos+12
       Do i_prt=1,4
!-------Boucle sur les 3 valeurs de chaque ligne  
          if (outdata%prtnumber(i_prt) .eq.number_prt(i_prt) ) then
             DO j=1,3
                IF ((DataForFiltering%prtcounts(j,i) .ne. NAN_I) .and. (DataForFiltering%prtcounts(j,i) .gt. 0))then
                   somme_for_mean_prt(i_prt)=somme_for_mean_prt(i_prt)+DataForFiltering%prtcounts(j,i)
                   n_points_for_mean_prt(i_prt)=n_points_for_mean_prt(i_prt)+1
                END IF
             END DO   
       
             DO j=1,2
         !print*, "j",j
                IF ((DataForFiltering%prtcounts(j+1,i) .ne. NAN_I) &
                     .and. (DataForFiltering%prtcounts(j+1,i) .ne. -999) &
                     .and. (DataForFiltering%prtcounts(j,i) .ne. NAN_I) &
                     .and. (DataForFiltering%prtcounts(j,i) .ne. -999)) THEN
                   somme_prt(i_prt)=somme_prt(i_prt)+(DataForFiltering%prtcounts(j+1,i)-DataForFiltering%prtcounts(j,i))**2
                   n_allan_prt(i_prt)=n_allan_prt(i_prt)+1
                END IF
             END DO
          END IF!prtnumber
       END DO!i_prt
    END DO ! scanlines
    Do i_prt=1,4
!-----On en deduit les Allan deviation    
      IF (n_points_for_mean_prt(i_prt) .gt. 0) THEN
        mean_prt(i_prt)=somme_for_mean_prt(i_prt)/n_points_for_mean_prt(i_prt)
      END IF
      IF (somme_prt(i_prt) .eq. 0) THEN
        allan_variance_prt(i_prt)= 0 
      ELSE IF (somme_prt(i_prt) .gt. 0) THEN
        allan_variance_prt(i_prt)= sqrt(somme_prt(i_prt)/(2.*n_points_for_mean_prt(i_prt))) 
      ELSE
        allan_variance_prt(i_prt)=NAN_R
      END IF
     
!-----On calcule les seuils
      IF (allan_variance_prt(i_prt) .ne. NAN_R) THEN
         thresh_high_prt(i_prt)=mean_prt(i_prt)+z*allan_variance_prt(i_prt)
         thresh_low_prt(i_prt)=mean_prt(i_prt)-z*allan_variance_prt(i_prt)
      END IF
     
    ENd DO !prt
   !---On trie les valeurs  et on calcule les moyennes par scanline  
    n_mean_prt=(/0,0,0,0/)
    DO i_prt=1,4
       if (DataForFiltering%prtnumber(i_prt) .eq.number_prt(i_prt) ) then
          DO j=1,3
             IF (allan_variance_prt(i_prt) .ne. NAN_R) THEN
                IF ( (DataForFiltering%prtcounts(j,scanLinePos) .lt. thresh_low_prt(i_prt)) &
                 .or. (DataForFiltering%prtcounts(j,scanLinePos) .gt. thresh_high_prt(i_prt))) THEN
                   DataForFiltering%prtcounts(j,scanLinePos)=-999
                ELSE
                   n_mean_prt(i_prt)=n_mean_prt(i_prt)+1
                END IF
             ELSE
               DataForFiltering%prtcounts(j,scanLinePos)=-999
             END IF
          END DO
       end if
    END DO ! boucle sur les 10 valeurs 4 prt
end if ! scanline pos
   
outData%prtcounts(1,scanLinePos)=DataForFiltering%prtcounts(1,scanLinePos)
outData%prtcounts(2,scanLinePos)=DataForFiltering%prtcounts(2,scanLinePos)
outData%prtcounts(3,scanLinePos)=DataForFiltering%prtcounts(3,scanLinePos)
 
END SUBROUTINE filtering_prt
 

END MODULE module_filtering_thresholds_prt

