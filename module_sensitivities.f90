MODULE module_sensitivities
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

 SUBROUTINE calcul_sensitivites(i,outData,coefs1, coefs2, coefs3)
!This subroutine is called every scanline

     INTEGER, INTENT(IN)                          :: i
     real,dimension(8,2),INTENT(IN)                 :: coefs1,coefs2,coefs3
     TYPE(AVHRR_Data), INTENT (INOUT)             :: outData
 
     REAL                                         :: cs_cict_3, cs_cict_4,cs_cict_5, &
                                                     u,tstar ! variable intermediaire
     INTEGER                                      :: j !indice de boucle
     REAL, DIMENSION(3)                           :: b0,a1,a2,b1


!####################
!Coefs : We use pre-launch values, for the moment the uncertainty =0
!
!             3  =  nuc
!             4  =  aVal (band coefficient)
!             5  =  bVal (band coefficient)
!####################
!Reminder of the names I used 
! sp[3-4-5] : average of the 10 space counts value per scanline
! bb[3-4-5] : average of the 10 ICT(internal calibration target) counts value per scanline   
! prtmean : average PRT temperature of the scanline
! nuc: centroid wavenumber
! Rict : radiance of the ICT 
! Counts[3-4-5] : Earth counts, 1 value per pixel 
!#####################
!Initialisation des sensibilites a NAN_R
!outdata%dre_over_dce3(:,:)=NAN_R 
!outdata%dre_over_dce4(:,:)=NAN_R 
!outdata%dre_over_dce5(:,:)=NAN_R 
!outdata%dre_over_dcs3(:,:)=NAN_R 
!outdata%dre_over_dcs4(:,:)=NAN_R 
!outdata%dre_over_dcs5(:,:)=NAN_R 
!outdata%dre_over_dcict3(:,:)=NAN_R 
!outdata%dre_over_dcict4(:,:)=NAN_R 
!outdata%dre_over_dcict5(:,:)=NAN_R 
!outdata%dre_over_db03(:,:)=NAN_R 
!outdata%dre_over_db04(:,:)=NAN_R 
!outdata%dre_over_db05(:,:)=NAN_R 
!outdata%dre_over_da13(:,:)=NAN_R 
!outdata%dre_over_da14(:,:)=NAN_R 
!outdata%dre_over_da15(:,:)=NAN_R 
!outdata%dre_over_da23(:,:)=NAN_R 
!outdata%dre_over_da24(:,:)=NAN_R 
!outdata%dre_over_da25(:,:)=NAN_R 
!outdata%dre_over_db13(:,:)=NAN_R 
!outdata%dre_over_db14(:,:)=NAN_R 
!outdata%dre_over_db15(:,:)=NAN_R 
!outdata%dre_over_dtinstr3(:,:)=NAN_R 
!outdata%dre_over_dtinstr4(:,:)=NAN_R 
!outdata%dre_over_dtinstr5(:,:)=NAN_R 
!outdata%dre_over_drict3(:,:)=NAN_R 
!outdata%dre_over_drict4(:,:)=NAN_R 
!outdata%dre_over_drict5(:,:)=NAN_R 
!outdata%dre_over_depsilon(:,:)=NAN_R 
!outdata%dtstar_over_daval3(:)=NAN_R 
!outdata%dtstar_over_daval4(:)=NAN_R 
!outdata%dtstar_over_daval5(:)=NAN_R 
!outdata%dtstar_over_dbval3(:)=NAN_R 
!outdata%dtstar_over_dbval4(:)=NAN_R 
!outdata%dtstar_over_dbval5(:)=NAN_R 
!outdata%dtstar_over_dnuc3(:)=NAN_R 
!outdata%dtstar_over_dnuc4(:)=NAN_R 
!outdata%dtstar_over_dnuc5(:)=NAN_R 
!outdata%drict_over_dnuc3(:)=NAN_R 
!outdata%drict_over_dnuc4(:)=NAN_R 
!outdata%drict_over_dnuc5(:)=NAN_R 
!outdata%drict_over_dtstar3(:)=NAN_R 
!outdata%drict_over_dtstar4(:)=NAN_R 
!outdata%drict_over_dtstar5(:)=NAN_R
 
!-We define values that are the same for all the pixels of a scanline 
 if ( (outdata%smoothsp3(i) .ne. NAN_R) .and. (outData%smoothbb3(i) .ne. NAN_R)&
.and.(outdata%smoothsp3(i) .gt. 0) .and. (outData%smoothbb3(i) .gt. 0))  THEN
     cs_cict_3=(outdata%smoothsp3(i)-outData%smoothbb3(i))
  end if
  if ( (outdata%smoothsp4(i) .ne. NAN_R) .and. (outData%smoothbb4(i) .ne.NAN_R)&
.and.(outdata%smoothsp4(i) .gt. 0) .and. (outData%smoothbb4(i) .gt. 0))  THEN
     cs_cict_4=(outdata%smoothsp4(i)-outData%smoothbb4(i))
  end if
  if ( (outdata%smoothsp5(i) .ne. NAN_R) .and. (outData%smoothbb5(i) .ne. NAN_R)&
.and.(outdata%smoothsp5(i) .gt. 0) .and. (outData%smoothbb5(i) .gt. 0))  THEN
     cs_cict_5=(outdata%smoothsp5(i)-outData%smoothbb5(i))
  end if
!########

!## ch3
!###  sensitivities of Tstar ############
outdata%dtstar_over_daval3(i)=1./coefs1(5,1)
outdata%dtstar_over_daval4(i)=1./coefs2(5,1)
outdata%dtstar_over_daval5(i)=1./coefs3(5,1)

if ((outdata%smoothprt(i) .ne. NAN_R) .and. (outdata%smoothprt(i) .gt. 0))then  
  outdata%dtstar_over_dbval3(i)=(outdata%smoothprt(i)-coefs1(4,1))/coefs1(5,1)
  outdata%dtstar_over_dbval4(i)=(outdata%smoothprt(i)-coefs2(4,1))/coefs2(5,1)
  outdata%dtstar_over_dbval5(i)=(outdata%smoothprt(i)-coefs3(4,1))/coefs3(5,1)
end if 

outdata%dtstar_over_dnuc3(i)=1./coefs1(5,1)
outdata%dtstar_over_dnuc4(i)=1./coefs2(5,1)
outdata%dtstar_over_dnuc5(i)=1./coefs3(5,1)

!###  sensibilities of Rict ############
if ((outdata%smoothprt(i) .ne. NAN_R) .and. (outdata%smoothprt(i) .gt. 0)) then 
  tstar=(outdata%smoothprt(i)-coefs1(4,1))/coefs1(5,1)
  u=c2*coefs1(3,1)/tstar
  outdata%drict_over_dnuc3(i)=c1*coefs1(3,1)**3/(exp(u)-1) &
*(3-coefs1(3,1)*(c2/tstar)**2*exp(u))
  outdata%drict_over_dtstar3(i)=(c2*coefs1(3,1)/tstar**2)**2&
*exp(u)/(exp(u)-1)**2

  tstar=(outdata%smoothprt(i)-coefs2(4,1))/coefs2(5,1)
  u=c2*coefs2(3,1)/tstar
  outdata%drict_over_dnuc4(i)=c1*coefs2(3,1)**3/(exp(u)-1) &
*(3-coefs2(3,1)*(c2/tstar)**2*exp(u))
  outdata%drict_over_dtstar4(i)=(c2*coefs2(3,1)/tstar**2)**2&
*exp(u)/(exp(u)-1)**2

  tstar=(outdata%smoothprt(i)-coefs3(4,1))/coefs3(5,1)
  u=c2*coefs3(3,1)/tstar
  outdata%drict_over_dnuc5(i)=c1*coefs3(3,1)**3/(exp(u)-1) &
*(3-coefs3(3,1)*(c2/tstar)**2*exp(u))
  outdata%drict_over_dtstar5(i)=(c2*coefs3(3,1)/tstar**2)**2&
*exp(u)/(exp(u)-1)**2
end if
!###################################### 

!--We calculate sensitivities that change from a pixel to another 
!-Loop over the 409 pixels of a scan line
  do j=1,409
!---Ch 3
!----Sensitivities of earth radiance
     if ((outdata%smoothsp3(i) .ne. NAN_R) .and. &
(outdata%smoothsp3(i) .gt. 0) .and. &
(outdata%rict_c3(i) .ne. NAN_R) .and. &
(outdata%rict_c3(i) .gt. 0) .and. &
(outData%Counts3(j,i).ne. NAN_R) .and. &
(outdata%Counts3(j,i) .gt. 0) .and. &
(outdata%orbital_temperature .ne. NAN_R) .and. &
(outdata%orbital_temperature .gt. 0) .and. &
(cs_cict_3 .ne. NAN_R) .and. (cs_cict_3 .gt. 0)) then 
    
        outdata%dre_over_db03(j,i)=1
      
        outdata%dre_over_da13(j,i)=outdata%Rict_c3(i)*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
           /cs_cict_3
      
        outdata%dre_over_da23(j,i)=-cs_cict_3*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
           +(outdata%smoothsp3(i)-outData%Counts3(j,i))**2

        outdata%dre_over_db13(j,i)=outdata%orbital_temperature
        outdata%dre_over_dtinstr3(j,i)=coefs1(8,1)
                     
        outdata%dre_over_dce3(j,i)=-((eta_ict+coefs1(2,1))*outdata%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
           /cs_cict_3 &
           -2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))
      
       outdata%dre_over_dcict3(j,i)=-((eta_ict+coefs1(2,1))*outdata%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
            *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
           /cs_cict_3**2&
           +2*coefs1(4,1)*(outdata%smoothsp3(i)-outData%Counts3(j,i))
      
       outdata%dre_over_dcs3(j,i)=-((eta_ict+coefs1(2,1))*outData%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
           *(outdata%smoothsp3(i)-outData%Counts3(j,i))&
           /cs_cict_3**2&
           +((eta_ict+coefs1(2,1))*outData%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2)&
           /cs_cict_3
      
      outdata%dre_over_drict3(j,i)=(eta_ict+coefs1(2,1))*(outdata%smoothsp3(i)-outData%Counts3(j,i))&
                                   /cs_cict_3
    end if 

!---Ch4
!----Sensitivities  of earth radiance
     if ((outdata%smoothsp4(i) .ne. NAN_R) .and. &
(outdata%smoothsp4(i) .gt. 0) .and. &
(outdata%rict_c4(i) .ne. NAN_R) .and. &
(outdata%rict_c4(i) .gt. 0) .and. &
(outData%Counts4(j,i).ne. NAN_R) .and. &
(outdata%Counts4(j,i) .gt. 0) .and. &
(outdata%orbital_temperature .ne. NAN_R) .and. &
(outdata%orbital_temperature .gt. 0) .and. &
(cs_cict_4 .ne. NAN_R) .and. (cs_cict_4 .gt. 0)) then 
     
        outdata%dre_over_db04(j,i)=1
      
        outdata%dre_over_da14(j,i)=outdata%Rict_c4(i)*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
           /cs_cict_4
      
        outdata%dre_over_da24(j,i)=-cs_cict_4*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
           +(outdata%smoothsp4(i)-outData%Counts4(j,i))**2

        outdata%dre_over_db14(j,i)=outdata%orbital_temperature
        outdata%dre_over_dtinstr4(j,i)=coefs2(8,1)
                     
        outdata%dre_over_dce4(j,i)=-((eta_ict+coefs2(2,1))*outdata%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
           /cs_cict_4 &
           -2*coefs2(4,1)*(outdata%smoothsp4(i)-outData%Counts4(j,i))
      
       outdata%dre_over_dcict4(j,i)=-((eta_ict+coefs2(2,1))*outdata%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
            *(outdata%smoothsp4(i)-outData%Counts4(j,i))&
           /cs_cict_4**2&
           +2*coefs2(4,1)*(outdata%smoothsp4(i)-outData%Counts4(j,i))
      
       outdata%dre_over_dcs4(j,i)=-((eta_ict+coefs2(2,1))*outData%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
           *(outdata%smoothsp4(i)-outData%Counts4(j,i))&
           /cs_cict_4**2&
           +((eta_ict+coefs2(2,1))*outData%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2)&
           /cs_cict_4
      
      outdata%dre_over_drict4(j,i)=(eta_ict+coefs2(2,1))*(outdata%smoothsp4(i)-outData%Counts4(j,i))&
                                   /cs_cict_4
    end if 

!---Ch5
!----Sensitivities  of earth radiance
     if ((outdata%smoothsp5(i) .ne. NAN_R) .and. &
(outdata%smoothsp5(i) .gt. 0) .and. &
(outdata%rict_c5(i) .ne. NAN_R) .and. &
(outdata%rict_c5(i) .gt. 0) .and. &
(outData%Counts5(j,i).ne. NAN_R) .and. &
(outdata%Counts5(j,i) .gt. 0) .and. &
(outdata%orbital_temperature .ne. NAN_R) .and. &
(outdata%orbital_temperature .gt. 0) .and. &
(cs_cict_5 .ne. NAN_R) .and. (cs_cict_5 .gt. 0)) then 
     
        outdata%dre_over_db05(j,i)=1
      
        outdata%dre_over_da15(j,i)=outdata%Rict_c5(i)*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
           /cs_cict_5
      
        outdata%dre_over_da25(j,i)=-cs_cict_5*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
           +(outdata%smoothsp5(i)-outData%Counts5(j,i))**2

        outdata%dre_over_db15(j,i)=outdata%orbital_temperature
        outdata%dre_over_dtinstr5(j,i)=coefs3(8,1)
                     
        outdata%dre_over_dce5(j,i)=-((eta_ict+coefs3(2,1))*outdata%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
           /cs_cict_5 &
           -2*coefs2(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))
      
       outdata%dre_over_dcict5(j,i)=-((eta_ict+coefs3(2,1))*outdata%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2) &
            *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
           /cs_cict_5**2&
           +2*coefs2(4,1)*(outdata%smoothsp5(i)-outData%Counts5(j,i))
      
       outdata%dre_over_dcs5(j,i)=-((eta_ict+coefs3(2,1))*outData%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
           *(outdata%smoothsp5(i)-outData%Counts5(j,i))&
           /cs_cict_5**2&
           +((eta_ict+coefs3(2,1))*outData%Rict_c5(i)-coefs2(4,1)*cs_cict_5**2)&
           /cs_cict_5
      
      outdata%dre_over_drict5(j,i)=(eta_ict+coefs3(2,1))*(outdata%smoothsp5(i)-outData%Counts5(j,i))&
                                   /cs_cict_5
    end if 
  end do !409
END SUBROUTINE calcul_sensitivites
 
END MODULE module_sensitivities

