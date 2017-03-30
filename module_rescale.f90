module module_rescale

  use netcdf
  use type

implicit none
PUBLIC

CONTAINS

SUBROUTINE rescale(AVHRR,Nscans)
  INTEGER, INTENT(IN)            :: Nscans
  TYPE(AVHRR_Data),INTENT(INOUT)    :: AVHRR
  
 
  REAL                           :: real_fillvalue=9.96921e+36
  INTEGER                        :: integer_fillvalue=-32767
  INTEGER                        ::i,j



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
print*, "rescaling"
if (AVHRR%orbital_temperature .eq. real_fillvalue) then
							
							
         AVHRR%orbital_temperature=integer_fillvalue
       else
         AVHRR%orbital_temperature=anint((AVHRR%orbital_temperature-bt_add_offset)/bt_scale_factor)
         if  ((AVHRR%orbital_temperature .lt. -20000)  .or. (AVHRR%orbital_temperature .gt. 10000)) then
          AVHRR%orbital_temperature=integer_fillvalue
         end if
       end if
 print*,  AVHRR%orbital_temperature
  Do i=1, Nscans
     Do j=1,409
       if (AVHRR%btf3(j,i) .eq. real_fillvalue) then
         AVHRR%btf3(j,i)=integer_fillvalue
       else
         AVHRR%btf3(j,i)=anint((AVHRR%btf3(j,i)-bt_add_offset)/bt_scale_factor)
         if  ((AVHRR%btf3(j,i) .lt. -20000)  .or. (AVHRR%btf3(j,i) .gt. 10000)) then
          AVHRR%btf3(j,i)=integer_fillvalue
         end if
       end if
    if (AVHRR%btf4(j,i) .eq. real_fillvalue) then
         AVHRR%btf4(j,i)=integer_fillvalue
       else
         AVHRR%btf4(j,i)=anint((AVHRR%btf4(j,i)-bt_add_offset)/bt_scale_factor)
         if  ((AVHRR%btf4(j,i) .lt. -20000)  .or. (AVHRR%btf4(j,i) .gt. 10000)) then
          AVHRR%btf4(j,i)=integer_fillvalue
         end if
       end if
    if (AVHRR%btf5(j,i) .eq. real_fillvalue) then
         AVHRR%btf5(j,i)=integer_fillvalue
    else
         AVHRR%btf5(j,i)=anint((AVHRR%btf5(j,i)-bt_add_offset)/bt_scale_factor)
         if  ((AVHRR%btf5(j,i) .lt. -20000)  .or. (AVHRR%btf5(j,i) .gt. 10000)) then
          AVHRR%btf5(j,i)=integer_fillvalue
         end if
    end if
 !IF( PRESENT(output_radiances) )THEN
      
    if (AVHRR%array3A(j,i) .eq. real_fillvalue) then
         AVHRR%array3A(j,i)=integer_fillvalue
    else
         AVHRR%array3A(j,i)=anint((AVHRR%array3A(j,i)-ref_add_offset)/ref_scale_factor)
         if  ((AVHRR%array3A(j,i) .lt. -20000)  .or. (AVHRR%array3A(j,i) .gt. 10000)) then
          AVHRR%array3A(j,i)=integer_fillvalue
         end if
    end if
    if (AVHRR%array1(j,i) .eq. real_fillvalue) then
         AVHRR%array1(j,i)=integer_fillvalue
    else
         AVHRR%array1(j,i)=anint((AVHRR%array1(j,i)-ref_add_offset)/ref_scale_factor)
         if  ((AVHRR%array1(j,i) .lt. -20000)  .or. (AVHRR%array1(j,i) .gt. 10000)) then
          AVHRR%array1(j,i)=integer_fillvalue
         end if
    end if
    if (AVHRR%array2(j,i) .eq. real_fillvalue) then
         AVHRR%array2(j,i)=integer_fillvalue
    else
         AVHRR%array2(j,i)=anint((AVHRR%array2(j,i)-ref_add_offset)/ref_scale_factor)
         if  ((AVHRR%array2(j,i) .lt. -20000)  .or. (AVHRR%array2(j,i) .gt. 10000)) then
          AVHRR%array2(j,i)=integer_fillvalue
         end if
    end if
    if (AVHRR%satZA(j,i) .eq. real_fillvalue) then
         AVHRR%satZA(j,i)=integer_fillvalue
    else
         AVHRR%satZA(j,i)=anint((AVHRR%satZA(j,i)-angle_add_offset)/angle_scale_factor)
         if  ((AVHRR%satZA(j,i) .lt. 0)  .or. (AVHRR%satZA(j,i) .gt. 9000)) then
          AVHRR%satZA(j,i)=integer_fillvalue
         end if
    end if
   if (AVHRR%solZA(j,i) .eq. real_fillvalue) then
         AVHRR%solZA(j,i)=integer_fillvalue
    else
         AVHRR%solZA(j,i)=anint((AVHRR%solZA(j,i)-angle_add_offset)/angle_scale_factor)
         if  ((AVHRR%solZA(j,i) .lt. 0)  .or. (AVHRR%solZA(j,i) .gt. 18000)) then
          AVHRR%solZA(j,i)=integer_fillvalue
         end if
    end if
   

     END do
 end do

  
 
 
    
end subroutine rescale


end module module_rescale
