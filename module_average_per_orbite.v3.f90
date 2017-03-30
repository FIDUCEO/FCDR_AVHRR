MODULE module_average_per_orbite

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

   SUBROUTINE average_per_orbite(outData,NScans) 

 TYPE(AVHRR_Data), INTENT (INOUT)             :: outData
 TYPE(AVHRR_Instrument_Coefs)                 :: Coefs
 INTEGER, INTENT(IN)                          :: NScans
   
 real                                         :: somme_patch, somme_patchExtended, &
                                                 npatch, npatchExtended
 real                                         :: gain_c4, gain_c5, gain_c3
 integer                                      :: i,j
 REAL                                         :: somme_allan_sp3,somme_allan_sp4,somme_allan_sp5, &
                                                 somme_allan_bb3, somme_allan_bb4,somme_allan_bb5, &
                                                 n_sp3, n_sp4, n_sp5, &
                                                 n_bb3, n_bb4, n_bb5

 real                                         :: somme_c3, stdev_c3, somme_bb_c3,stdev_bb_c3, &
                                              somme_space_counts_c3, somme_bb_counts_c3, somme_rict_c3, &
                                              nmean_space_counts_c3, nmean_bb_counts_c3, nmean_rict_c3, &
                                              somme_2_c3, n_2_c3

 real                                         :: somme_c4, stdev_c4, somme_bb_c4,stdev_bb_c4, &
                                              somme_space_counts_c4, somme_bb_counts_c4, somme_rict_c4, &
                                              nmean_space_counts_c4, nmean_bb_counts_c4, nmean_rict_c4, &
                                              somme_2_c4, n_2_c4

 real                                         :: somme_c5, stdev_c5, somme_bb_c5,stdev_bb_c5, &
                                              somme_space_counts_c5, somme_bb_counts_c5, somme_rict_c5, &
                                              nmean_space_counts_c5, nmean_bb_counts_c5, nmean_rict_c5, &
                                              somme_2_c5, n_2_c5
                                                

 real                                         :: somme_prtmean,nmean_prtmean
 real                                         :: somme_calib3_2, somme_calib4_2, somme_calib5_2
 real                                         :: nmean_calib3_2, nmean_calib4_2, nmean_calib5_2
   

somme_patch=0
somme_patchExtended=0
npatch=0
npatchExtended=0

gain_c3=-32768
somme_rict_c3=0
nmean_rict_c3=0
somme_space_counts_c3=0
nmean_space_counts_c3=0
somme_bb_counts_c3=0
nmean_bb_counts_c3=0 
somme_2_c3=0
n_2_c3=0

gain_c4=-32768
somme_rict_c4=0
nmean_rict_c4=0
somme_space_counts_c4=0
nmean_space_counts_c4=0
somme_bb_counts_c4=0
nmean_bb_counts_c4=0 
somme_2_c4=0
n_2_c4=0

gain_c5=-32768
somme_rict_c5=0
nmean_rict_c5=0
somme_space_counts_c5=0
nmean_space_counts_c5=0
somme_bb_counts_c5=0
nmean_bb_counts_c5=0  
somme_2_c5=0
n_2_c5=0
           
somme_prtmean=0
nmean_prtmean=0 

somme_calib3_2=0
nmean_calib3_2=0

somme_calib4_2=0
nmean_calib4_2=0
           
somme_calib5_2=0
nmean_calib5_2=0   

somme_allan_sp3=0
somme_allan_sp4=0
somme_allan_sp5=0
somme_allan_bb3=0
somme_allan_bb4=0
somme_allan_bb5=0
n_sp3=0
n_sp4=0
n_sp5=0
n_bb3=0
n_bb4=0
n_bb5=0  

outdata%ucs3=NAN_R
outdata%ucs4=NAN_R
outdata%ucs5=NAN_R
outdata%ucict3=NAN_R
outdata%ucict4=NAN_R
outdata%ucict5=NAN_R

!calcule les moyennes par orbite
print*, "Nscans", Nscans
Do i=1, NScans
  do j=1,9
    if ((outData%spaceFilter3(j,i) .ne. -32768) .and. (outData%spaceFilter3(j+1,i) .ne. -32768)) then
      somme_allan_sp3=somme_allan_sp3+(outData%spaceFilter3(j+1,i)-outData%spaceFilter3(j,i))**2
      n_sp3=n_sp3+1
    end if   
    if ((outData%spaceFilter4(j,i) .ne. -32768) .and. (outData%spaceFilter4(j+1,i) .ne. -32768)) then
      somme_allan_sp4=somme_allan_sp4+(outData%spaceFilter4(j+1,i)-outData%spaceFilter4(j,i))**2
      n_sp4=n_sp4+1
    end if 
    if ((outData%spaceFilter5(j,i) .ne. -32768) .and. (outData%spaceFilter5(j+1,i) .ne. -32768)) then
      somme_allan_sp5=somme_allan_sp5+(outData%spaceFilter5(j+1,i)-outData%spaceFilter5(j,i))**2
      n_sp5=n_sp5+1
    end if
    if ((outData%bbodyFilter3(j,i) .ne. -32768) .and. (outData%bbodyFilter3(j+1,i) .ne. -32768)) then
      somme_allan_bb3=somme_allan_sp3+(outData%bbodyFilter3(j+1,i)-outData%bbodyFilter3(j,i))**2
      n_bb3=n_bb3+1
    end if 
    if ((outData%bbodyFilter4(j,i) .ne. -32768) .and. (outData%bbodyFilter4(j+1,i) .ne. -32768)) then
      somme_allan_bb4=somme_allan_sp4+(outData%bbodyFilter4(j+1,i)-outData%bbodyFilter4(j,i))**2
      n_bb4=n_bb4+1
    end if  
    if ((outData%bbodyFilter5(j,i) .ne. -32768) .and. (outData%bbodyFilter5(j+1,i) .ne. -32768)) then
      somme_allan_bb5=somme_allan_sp5+(outData%bbodyFilter5(j+1,i)-outData%bbodyFilter5(j,i))**2
      n_bb5=n_bb5+1
    end if
  end do !j
 
  if (outData%patch(i) .ne. NAN_R) then
    somme_patch=somme_patch+outData%patch(i)
    npatch=npatch+1
  end if
  if (outData%patchExtended(i) .ne. NAN_R) then
    somme_patchExtended=somme_patchExtended+outData%patchExtended(i)
    npatchExtended=npatchExtended+1
  end if
  if  (outData%prtmean(i) .ne. NAN_R) then 
      somme_prtmean=somme_prtmean+outData%smoothprt(i)
      nmean_prtmean=nmean_prtmean+1
  end if

  if  ((outData%smoothsp3(i) .ne. NAN_R) .and. (outData%smoothbb3(i) .ne. NAN_R)) then 
    somme_space_counts_c3=somme_space_counts_c3+outData%smoothsp3(i)
    nmean_space_counts_c3=nmean_space_counts_c3+1
    somme_bb_counts_c3=somme_bb_counts_c3+outData%smoothbb3(i)
    nmean_bb_counts_c3=nmean_bb_counts_c3+1
    do j=1,10
      somme_2_c3=somme_2_c3+outData%spaceFilter3(j,i)
      n_2_c3=n_2_c3+1
    end do
  end if 
  if  (outData%Rict_c3(i) .ne. NAN_R) then 
    somme_rict_c3=somme_rict_c3+outData%rict_c3(i)
    nmean_rict_c3=nmean_rict_c3+1
  end if
   
  if  ((outData%smoothsp4(i) .ne. NAN_R) .and. (outData%smoothbb4(i) .ne. NAN_R)) then 
    somme_space_counts_c4=somme_space_counts_c4+outData%smoothsp4(i)
    nmean_space_counts_c4=nmean_space_counts_c4+1
    somme_bb_counts_c4=somme_bb_counts_c4+outData%smoothbb4(i)
    nmean_bb_counts_c4=nmean_bb_counts_c4+1
    do j=1,10
      somme_2_c4=somme_2_c4+outData%spaceFilter4(j,i)
      n_2_c4=n_2_c4+1
    end do
  end if 
  if  (outData%Rict_c4(i) .ne. NAN_R) then 
    somme_rict_c4=somme_rict_c4+outData%rict_c4(i)
    nmean_rict_c4=nmean_rict_c4+1
  end if
  
!--Channel 5
  if  ((outData%smoothsp5(i) .ne. NAN_R) .and. (outData%smoothbb5(i) .ne. NAN_R)) then 
    somme_space_counts_c5=somme_space_counts_c5+outData%smoothsp5(i)
    nmean_space_counts_c5=nmean_space_counts_c5+1
    somme_bb_counts_c5=somme_bb_counts_c5+outData%smoothbb5(i)
    nmean_bb_counts_c5=nmean_bb_counts_c5+1
    do j=1,10
      somme_2_c5=somme_2_c5+outData%spaceFilter5(j,i)
      n_2_c5=n_2_c5+1
    end do
  end if 
  if  (outData%Rict_c5(i) .ne. NAN_R) then 
    somme_rict_c5=somme_rict_c5+outData%rict_c5(i)
    nmean_rict_c5=nmean_rict_c5+1
  end if
end do !scanlines
                    
if ( (somme_patch .gt. 0) .and. (npatch .gt. 0) ) then 
   outdata%mean_patch_orbite=somme_patch/npatch
end if 
         
if ( (somme_patchExtended .gt. 0) .and. (npatchExtended .gt. 0) ) then 
   outdata%mean_patchExtended_orbite=somme_patchExtended/npatchExtended
end if 

if ( (somme_prtmean .gt. 0) .and. (nmean_prtmean .gt. 0) ) then 
   outdata%mean_prt_orbite=somme_prtmean/nmean_prtmean
end if 

if ( (somme_rict_c3 .gt. 0) .and. (nmean_rict_c3 .gt. 0) ) then 
   outdata%mean_rict_orbite_c3=somme_rict_c3/nmean_rict_c3
end if 
          
if ( (somme_space_counts_c3 .gt. 0) .and. (nmean_space_counts_c3 .gt. 0) ) then 
   outData%mean_space_counts_orbite_c3=somme_space_counts_c3/nmean_space_counts_c3
end if 
       
if ( (somme_bb_counts_c3 .gt. 0) .and. (nmean_bb_counts_c3 .gt. 0) ) then 
   outdata%mean_bb_counts_orbite_c3=somme_bb_counts_c3/nmean_bb_counts_c3
end if
         
if ( (somme_rict_c4 .gt. 0) .and. (nmean_rict_c4 .gt. 0) ) then 
   outdata%mean_rict_orbite_c4=somme_rict_c4/nmean_rict_c4
end if 
           
if ( (somme_space_counts_c4 .gt. 0) .and. (nmean_space_counts_c4 .gt. 0) ) then 
   outData%mean_space_counts_orbite_c4=somme_space_counts_c4/nmean_space_counts_c4
end if 
   
if ( (somme_bb_counts_c4 .gt. 0) .and. (nmean_bb_counts_c4 .gt. 0) ) then 
   outdata%mean_bb_counts_orbite_c4=somme_bb_counts_c4/nmean_bb_counts_c4
end if 


if ( (somme_rict_c5 .gt. 0) .and. (nmean_rict_c5 .gt. 0) ) then 
   outdata%mean_rict_orbite_c5=somme_rict_c5/nmean_rict_c5
end if 

if ( (somme_space_counts_c5 .gt. 0) .and. (nmean_space_counts_c5 .gt. 0) ) then 
   outData%mean_space_counts_orbite_c5=somme_space_counts_c5/nmean_space_counts_c5           
end if

if ( (somme_bb_counts_c5 .gt. 0) .and. (nmean_bb_counts_c5 .gt. 0) ) then 
   outdata%mean_bb_counts_orbite_c5=somme_bb_counts_c5/nmean_bb_counts_c5
end if 

if (n_bb3>10) then
  print*, somme_allan_bb3, n_bb3
  outdata%ucict3= SQRT(somme_allan_bb3/(2*(n_bb3-1)))
  print*, outdata%ucict3
end if  
if (n_bb4>10) then
  outdata%ucict4= SQRT(somme_allan_bb4/(2*(n_bb5-1)))
end if  
if (n_bb5>10) then
  outdata%ucict5= SQRT(somme_allan_bb5/(2*(n_bb5-1)))
end if  
if (n_sp3>10) then
  outdata%ucs3= SQRT(somme_allan_sp3/(2*(n_sp3-1)))
end if  
if (n_sp4>10) then
  outdata%ucs4= SQRT(somme_allan_sp4/(2*(n_sp4-1)))
end if  
if (n_sp5>10) then
  outdata%ucs5= SQRT(somme_allan_sp5/(2*(n_sp5-1)))
end if

END SUBROUTINE average_per_orbite
 

END MODULE module_average_per_orbite

