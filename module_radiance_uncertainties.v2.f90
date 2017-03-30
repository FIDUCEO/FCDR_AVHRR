MODULE module_radiance_uncertainties
  USE GbcsKInds
  USE GbcsConstants
  USE GbcsTypes
  USE GbcsSystemTools
  USE GbcsErrorHandler
  USE GbcsStringUtil
  USE type
  use module_sensitivities

implicit none
PUBLIC

CONTAINS

   SUBROUTINE radiance_uncertainties(i,outData,coefs1, coefs2, coefs3)

INTEGER, INTENT(IN)                  :: i
REAL,dimension(8,2),INTENT(IN)       :: coefs1,coefs2,coefs3
TYPE(AVHRR_Data), INTENT (INOUT)     :: outData
  
REAL                                 :: somme_allan_sp3,somme_allan_sp4,somme_allan_sp5, &
                                        somme_allan_bb3, somme_allan_bb4,somme_allan_bb5, &
                                        n_sp3, n_sp4, n_sp5, &
                                        n_bb3, n_bb4, n_bb5, &
                                        somme_ve3, somme_ve4, somme_ve5, &
                                        somme_vstr3,somme_vstr4,somme_vstr5, &
                                        cs_cict_3, cs_cict_4,cs_cict_5

INTEGER                              :: g,j, k,l
REAL                                 :: prt_accuracy=0.1, prt_noise=0.015, &
                                        uaval=0, ubval=0, unuc=0,&
                                        u

REAL                                 :: nrmax5,ur5, nrmin5, trmin5, trmax5, &
                                        nsmax5,us5, nsmin5, tsmin5, tsmax5, &
                                        nrmax3,ur3, nrmin3, trmin3, trmax3, &
                                        nsmax3,us3, nsmin3, tsmin3, tsmax3, &
                                        nrmax4,ur4, nrmin4, trmin4, trmax4, &
                                        nsmax4,us4, nsmin4, tsmin4, tsmax4
 
REAL, DIMENSION(3)                   :: u_b0, u_a1, u_a2,u_b1
REAL                                 :: real_fillvalue=9.96921e+36
INTEGER                              :: integer_fillvalue=-32767

u_b0(1)=covariance_matrix_ch3(1,1)
u_b0(2)=covariance_matrix_ch4(1,1)
u_b0(3)=covariance_matrix_ch5(1,1)

u_a1(1)=covariance_matrix_ch3(2,2)
u_a1(2)=covariance_matrix_ch4(2,2)
u_a1(3)=covariance_matrix_ch5(2,2)

u_a2(1)=covariance_matrix_ch3(3,3)
u_a2(2)=covariance_matrix_ch3(3,3)
u_a2(3)=covariance_matrix_ch3(3,3)

u_b1(1)=covariance_matrix_ch4(4,4)
u_b1(2)=covariance_matrix_ch3(4,4)
u_b1(3)=covariance_matrix_ch4(4,4)


!-On definit des grandeurs qui sont utiles et identiques pour tous les pixels d'une scan ligne
 if ( (outDAta%smoothsp3(i) .ne. NAN_R) .and. (outData%smoothbb3(i) .ne. NAN_R) &
.and.(outDAta%smoothsp3(i) .gt. 0) .and. (outData%smoothbb3(i) .gt.  0))  THEN
     cs_cict_3=(outData%smoothsp3(i)-outData%smoothbb3(i))
  end if
!print*, cs_cict_3
  if ( (outDAta%smoothsp4(i) .ne. NAN_R) .and. (outData%smoothbb4(i) .ne. NAN_R)&
.and. (outDAta%smoothsp4(i) .gt. 0) .and. (outData%smoothbb4(i) .gt. 0))  THEN
     cs_cict_4=(outData%smoothsp4(i)-outData%smoothbb4(i))
  end if
  if ( (outDAta%smoothsp5(i) .ne. NAN_R) .and. (outData%smoothbb5(i) .ne. NAN_R)&
.and. (outDAta%smoothsp5(i) .gt.0) .and. (outData%smoothbb5(i) .gt. 0))  THEN
     cs_cict_5=(outData%smoothsp5(i)-outData%smoothbb5(i))
  end if

!--ON cacul les incertitudes qui changent d'un pixel à l'autre           
  do j=1,409
     outData%ue3(j,i)=real_fillvalue
     outData%ue4(j,i)=real_fillvalue
     outData%ue5(j,i)=real_fillvalue
     outData%ur3(j,i)=real_fillvalue
     outData%ur4(j,i)=real_fillvalue
     outData%ur5(j,i)=real_fillvalue
     outData%us3(j,i)=real_fillvalue
     outData%us4(j,i)=real_fillvalue
     outData%us5(j,i)=real_fillvalue
     outData%btf3(j,i)=integer_fillvalue
     outData%btf4(j,i)=integer_fillvalue
     outData%btf5(j,i)=integer_fillvalue
     outData%uce3(j,i)=0.3
     outData%uce4(j,i)=0.3
     outData%uce5(j,i)=0.3
!---Ch 3
   !if ((outdata%ucs3 .ne. NAN_R) .and. (outData%ucict3 .ne. NAN_R)) then
   !    outdata%uce3(j,i)=(outData%ucs3+outData%ucict3)/2.
   ! end if
    !print*, "uce", outdata%uce3(j,i)
   if ((outData%smoothsp3(i) .ne. NAN_R) .and. (outdata%prtmean(i) .ne. NAN_R) .and. &
       (outData%Counts3(j,i).ne. NAN_R) .and. &
       (cs_cict_3 .ne. NAN_R)) then  
        
         ! outData%nef3(j,i)=coefs1(1,1)&
         !               +((eta_ict+coefs1(1,1))*outData%Rict_c3(i)-coefs1(4,1)*cs_cict_3**2) &
         !               *(outData%smoothsp3(i)-outData%Counts3(j,i)) &
         !               /cs_cict_3 &
         !               +coefs1(4,1)*((outData%smoothsp3(i)-outData%Counts3(j,i))**2) &
         !               +coefs1(8,1)* outData%orbital_temperature
      
       ur3=sqrt(outData%dre_over_dce3(j,i)**2*outData%uce3(j,i)**2)
       !print*,"sensitivites",outData%dre_over_drict3(j,i), outData%dre_over_dcs3(j,i), outData%dre_over_dcict3(j,i)
       !print*,"ucs,ucict",outData%ucs3, outData%ucict3           
       us3=sqrt(&!(outData%dre_over_drict3(j,i)**2*outData%urict3(i)**2) &
                       (outData%dre_over_dcs3(j,i)**2*outData%ucs3**2) &
                       +(outData%dre_over_dcict3(j,i)**2*outData%ucict3**2)) 
       !print*, "us3",us3
   end if

!---Ch 4
!   if ((outdata%ucs4 .ne. NAN_R) .and. (outData%ucict4 .ne. NAN_R)) then
!       outdata%uce4(j,i)=(outData%ucs4+outData%ucict4)/2.
!   end if

 !print*, eta_ict,outData%Rict_c4(i), cs_cict_4,outData%sp4(i),outData%Counts4(j,i),&
!                  outData%orbital_temperature
       !outData%nef4(j,i)=coefs2(1,1)&
       !                 +((eta_ict+coefs2(1,1))*outData%Rict_c4(i)-coefs2(4,1)*cs_cict_4**2) &
       !                 *(outData%smoothsp4(i)-outData%Counts4(j,i)) &
       !                 /cs_cict_4 &
       !                 +coefs2(4,1)*((outData%smoothsp4(i)-outData%Counts4(j,i))**2) &
       !                 +coefs2(8,1)*outData%orbital_temperature
     if ((outdata%uce4(j,i) .ne. NAN_R) &
.and. (outdata%uce4(j,i) .gt. 0) &
.and. (outdata%dre_over_dce4(j,i) .ne. NAN_R)) then 
       ur4=sqrt(outData%dre_over_dce4(j,i)**2*outData%uce4(j,i)**2)
 end if
                     
     if ((outdata%ucict4 .ne. NAN_R) &
.and. (outdata%ucict4 .gt. 0) &
.and. (outdata%ucs4 .ne. NAN_R) &
.and. (outdata%ucs4 .gt. 0) &
.and. (outdata%dre_over_drict4(j,i) .ne. NAN_R) & 
.and. (outdata%dre_over_dcs4(j,i) .ne. NAN_R) & 
.and. (outdata%dre_over_dcict4(j,i) .ne. NAN_R)) then 
       us4=sqrt((outData%dre_over_drict4(j,i)**2*outData%urict4(i)**2) &
                       +(outData%dre_over_dcs4(j,i)**2*outData%ucs4**2) &
                       +(outData%dre_over_dcict4(j,i)**2*outData%ucict4**2)) 
   end if

!---Ch 5
   !if ((outdata%ucs5 .ne. NAN_R) .and. (outData%ucict5 .ne. NAN_R)) then
   !    outdata%uce5(j,i)=(outData%ucs5+outData%ucict5)/2.
   !end if

   if ((outdata%smoothsp5(i) .ne. NAN_R) .and. (outdata%prtmean(i) .ne. NAN_R) .and. &
       (outData%Counts5(j,i).ne. NAN_R) .and. &
       (cs_cict_5 .ne. NAN_R)) then  
 
       !outData%nef5(j,i)=coefs3(1,1)&
       !                 +((eta_ict+coefs3(1,1))*outData%Rict_c5(i)-coefs3(4,1)*cs_cict_5**2) &
       !                 *(outdata%smoothsp5(i)-outData%Counts5(j,i)) &
       !                 /cs_cict_5 &
       !                 +coefs3(4,1)*((outdata%smoothsp5(i)-outData%Counts5(j,i))**2) &
       !                 +coefs3(8,1)*outData%orbital_temperature
      
       ur5=sqrt(outData%dre_over_dce5(j,i)**2*outData%uce5(j,i)**2)
                      
       us5=sqrt((outData%dre_over_drict5(j,i)**2*outData%urict5(i)**2) &
                       +(outData%dre_over_dcs5(j,i)**2*outData%ucs5**2) &
                       +(outData%dre_over_dcict5(j,i)**2*outData%ucict5**2))   
            
    end if
     
!---FIDUCEO : on convertit les radiances en BT
    if  ((outData%new_array3B(j,i) .gt. 0).and. (outData%new_array3B(j,i) .ne. NAN_R)) then
      outdata%btf3(j,i)=convertBT(outdata%new_array3B(j,i),coefs1(5,1), coefs1(6,1), coefs1(7,1))
      nrmax3=outdata%new_array3B(j,i)+ur3
      nrmin3=outdata%new_array3B(j,i)-ur3
      trmin3=convertBT(nrmin3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
      trmax3=convertBT(nrmax3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
      outdata%ur3(j,i)=(trmax3-trmin3)/2.
     
      nsmax3=outdata%new_array3B(j,i)+us3
      nsmin3=outdata%new_array3B(j,i)-us3
    
      tsmin3=convertBT(nsmin3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
      tsmax3=convertBT(nsmax3,coefs1(5,1), coefs1(6,1), coefs1(7,1))
      outdata%us3(j,i)=(tsmax3-tsmin3)/2.
    end if
      
    if  ((outData%new_array4(j,i) .gt. 0).and.(outData%new_array4(j,i) .ne. NAN_R)) then
      outdata%btf4(j,i)=convertBT(outdata%new_array4(j,i),coefs2(5,1), coefs2(6,1), coefs2(7,1))
      nrmax4=outdata%new_array4(j,i)+ur4
      nrmin4=outdata%new_array4(j,i)-ur4
      trmin4=convertBT(nrmin4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
      trmax4=convertBT(nrmax4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
      if (trmin4 .gt. trmax4) then
        print *, ur4, outdata%new_array4(j,i), outdata%btf4(j,i), trmin4, trmax4
      end if
      outdata%ur4(j,i)=(trmax4-trmin4)/2.
      nsmax4=outdata%new_array4(j,i)+us4
      nsmin4=outdata%new_array4(j,i)-us4
      tsmin4=convertBT(nsmin4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
      tsmax4=convertBT(nsmax4,coefs2(5,1), coefs2(6,1), coefs2(7,1))
      outdata%us4(j,i)=(tsmax4-tsmin4)/2.
      !print*, "ch4",outdata%new_array4(j,i),outData%btf4(j,i), outdata%ur4(j,i), outdata%us4(j,i)
    end if
      if  ((outData%new_array5(j,i) .gt. 0).and. (outData%new_array5(j,i) .ne. NAN_R)) then
      outdata%btf5(j,i)=convertBT(outdata%new_array5(j,i),coefs3(5,1), coefs3(6,1), coefs3(7,1))
      nrmax5=outdata%new_array5(j,i)+ur5
      nrmin5=outdata%new_array5(j,i)-ur5
      trmin5=convertBT(nrmin5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
      trmax5=convertBT(nrmax5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
      outdata%ur5(j,i)=(trmax5-trmin5)/2.
      nsmax5=outdata%new_array5(j,i)+us5
      nsmin5=outdata%new_array5(j,i)-us5
      tsmin5=convertBT(nsmin5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
      tsmax5=convertBT(nsmax5,coefs3(5,1), coefs3(6,1), coefs3(7,1))
      outdata%us5(j,i)=(tsmax5-tsmin5)/2.
      !print*, "ch5",outdata%new_array5(j,i),outData%btf5(j,i), outdata%ur5(j,i), outdata%us5(j,i)
    end if

  end do !409
!print*, minval(outdata%btf3), maxval(outdata%btf3)
contains 
 REAL FUNCTION convertBT( radiance, nuC, Aval, Bval )

    REAL, INTENT(IN) :: radiance
    REAL, INTENT(IN) :: nuC
    REAL, INTENT(IN) :: Aval
    REAL, INTENT(IN) :: Bval

    REAL :: constant

    constant = 1. + (C1*nuC*nuC*nuC/radiance)
    constant = C2*nuC/log(constant)
    convertBT = (aVal+bVal*constant)  
    return

  END FUNCTION convertBT
 END SUBROUTINE radiance_uncertainties
END MODULE module_radiance_uncertainties

