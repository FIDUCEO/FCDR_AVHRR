MODULE module_average_per_scanline

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

   SUBROUTINE average_per_scanline(scanLinePos,outData, Coefs) 

TYPE(AVHRR_Data), INTENT (INOUT)              :: outData
TYPE(AVHRR_Instrument_Coefs), INTENT(INOUT)   :: Coefs
INTEGER, INTENT(IN)                           :: scanLinePos

INTEGER                                       :: j ! indice de boucle
REAL                                          :: somme_c3, stdev_c3, somme_bb_c3,stdev_bb_c3, &
                                                 somme_c4, stdev_c4, somme_bb_c4,stdev_bb_c4, &
                                                 somme_c5, stdev_c5,somme_bb_c5, stdev_bb_c5
                                  

somme_c3=-32768
stdev_c3=-32768
somme_c4=-32768
stdev_c4=-32768
somme_c5=-32768
stdev_c5=-32768
somme_bb_c3=-32768
stdev_bb_c3=-32768
somme_bb_c4=-32768
stdev_bb_c4=-32768
somme_bb_c5=-32768
stdev_bb_c5=-32768  

!Channel 3B
if ((outData%Rict_c3(scanLinePos) .ne. NAN_R) &
.and. (outData%sp3(scanLinePos) .ne. NAN_R) &
.and. (outData%bb3(scanLinePos) .ne. NAN_R)) then
   outData%gain_c3(scanLinePos)=(outData%Rict_c3(scanLinePos)-Coefs%Nspace(1))/ &
                 (outData%bb3(scanLinePos)-outData%sp3(scanLinePos))
end if

if ((outData%sp3(scanLinePos) .ne. NAN_R) .and. (outData%bb3(scanLinePos) .ne. NAN_R) &
.and. (outData%sp3(scanLinePos) .ne. -999) .and. (outData%bb3(scanLinePos) .ne. -999)) then 
!--on calcule la standard deviation de Space counts ligne par ligne 
   somme_c3=0  
   do j=1,10
      somme_c3=somme_c3+(outData%spaceFilter3(j,scanLinePos)-outData%sp3(scanLinePos))**2
   end do
!  print*, scanLinePos, outData%spaceFilter5(:,scanLinePos), outData%sp5(scanLinePos), somme_c3 
   if (somme_c3 .ge. 8) then
      stdev_c3= sqrt(somme_c3/10) 
   end if
!--On calcule la standard deviation de BB counts ligne par ligne      
   somme_bb_c3=0  
   do j=1,10
      somme_bb_c3=somme_c3+(outData%bbodyFilter3(j,scanLinePos)-outData%bb3(scanLinePos))**2
   end do
!--print*, i, outData%spaceFilter5(:,i), outData%sp5(i), somme     
   if (somme_bb_c3 .ge. 8) then
      stdev_bb_c3= sqrt(somme_bb_c3/10) 
   end if
end if 

!-----Channel 4
if ((outData%Rict_c4(scanLinePos) .ne. NAN_R) &
.and. (outData%sp4(scanLinePos) .ne. NAN_R) &
.and. (outData%bb4(scanLinePos) .ne. NAN_R)) then
   outData%gain_c4(scanLinePos)=(outData%Rict_c4(scanLinePos)-Coefs%Nspace(2))&
              /(outData%bb4(scanLinePos)-outData%sp4(scanLinePos))
end if 

if ((outData%sp4(scanLinePos) .ne. NAN_R) .and. (outData%bb4(scanLinePos) .ne. NAN_R) &
.and. (outData%sp4(scanLinePos) .ne. -999) .and. (outData%bb3(4) .ne. -999)) then 
!--on calcule la standard deviation de Space counts ligne par ligne
   somme_c4=0
   do j=1,10
      somme_c4=somme_c4+(outData%spaceFilter4(j,scanLinePos)-outData%sp4(scanLinePos))**2
   end do   
!  print*, i, outData%spaceFilter5(:,i), outData%sp5(i), somme     
   if (somme_c4 .ge. 8) then
      stdev_c4= sqrt(somme_c4/10)  
   end if
!--on calcule la standard deviation de BB counts ligne par ligne      
   somme_bb_c4=0  
   do j=1,10
      somme_bb_c4=somme_c4+(outData%bbodyFilter4(j,scanLinePos)-outData%bb4(scanLinePos))**2
   end do   
!  print*, i, outData%spaceFilter5(:,i), outData%sp5(i), somme     
   if (somme_bb_c4 .ge. 8) then
      stdev_bb_c4= sqrt(somme_bb_c4/10) 
   end if
end if 

!-----Channel 5
if ((outData%Rict_c5(scanLinePos) .ne. NAN_R) &
.and. (outData%sp5(scanLinePos) .ne. NAN_R) &
.and. (outData%bb5(scanLinePos) .ne. NAN_R)) then
  outData%gain_c5(scanLinePos) =(outData%Rict_c5(scanLinePos)-Coefs%Nspace(3))/&
                 (outData%bb5(scanLinePos)-outData%sp5(scanLinePos))
end if 
     
if  ((outData%sp5(scanLinePos) .ne. NAN_R) .and. (outData%bb5(scanLinePos) .ne. NAN_R) &
.and. (outData%sp5(scanLinePos) .ne. -999) .and. (outData%bb5(scanLinePos) .ne. -999)) then 
   somme_c5=0 
   do j=1,10
      somme_c5=somme_c5+(outData%spaceFilter5(j,scanLinePos)-outData%sp5(scanLinePos))**2
   end do   
!  print*, i, outData%spaceFilter5(:,i), outData%sp5(i), somme     
   if (somme_c5 .ge. 8) then
      stdev_c5= sqrt(somme_c5/10)
   end if
!--on calcule la standard deviation de BB counts ligne par ligne      
   somme_bb_c5=0  
   do j=1,10
      somme_bb_c5=somme_bb_c5+(outData%bbodyFilter5(j,scanLinePos)-outData%bb5(scanLinePos))**2
   end do   
!  print*, i, outData%spaceFilter5(:,i), outData%sp5(i), somme     
   if (somme_bb_c5 .ge. 8) then
      stdev_bb_c5= sqrt(somme_bb_c5/10)
   end if
end if 
           

END SUBROUTINE average_per_scanline
END MODULE module_average_per_scanline

