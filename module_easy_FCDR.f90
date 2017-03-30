module module_FCDR

  use netcdf
  use type

implicit none
PUBLIC

CONTAINS

SUBROUTINE write_FCDR(FILE_NAME,AVHRR)
  character (len = *), INTENT(IN):: FILE_NAME
  TYPE(AVHRR_Data),INTENT(IN)    :: AVHRR
  
  integer                        :: ncid

  REAL                           :: real_fillvalue=9.96921e+36
  real*4                         ::  real_fillvalue_2=-1.00000002E+30  
  INTEGER                           :: integer_fillvalue=-32767
! We are writing 3D data
  integer, parameter             :: NDIMS = 2 
  integer, parameter             :: NTIME = 1, NI = 409, NJ = 13639, N10=10
  character (len = *), parameter :: I_NAME = "ni"
  character (len = *), parameter :: J_NAME = "nj"
  character (len = *), parameter :: TEN_NAME = "n10"

  integer                        :: time_dimid, ni_dimid, nj_dimid,n10_dimid
  integer                        :: dimids_1(1),dimids_sp_ict(2),dimids_2(2)

! The start and count arrays will tell the netCDF library where to write our data.
  integer                        :: start(NDIMS),start_10(2),start_1(1), &
                                    count(NDIMS),count_10(2),count_1(1)

  integer                                 :: chunk_size_1_per_line(1),&
                                             chunk_size_10_per_line(2),&
                                             chunk_size_pixel(2),  &
                                             chunk_size_3_per_line(2),&
                                             chunk_size_4_per_line(2)

  character (len = *), parameter :: TIME_NAME = "time"
  character (len = *), parameter :: SCAN_LINE_NUMBER_NAME="scanLineNumber"
  character (len = *), parameter :: LAT_NAME="Lat"
  character (len = *), parameter :: LON_NAME="Lon"

  character (len = *), parameter :: satza_NAME="satZA"
  character (len = *), parameter :: solza_NAME="solZA"
  character (len = *), parameter :: relaz_NAME="relAZ"

  character (len = *), parameter :: CH3B_UR_NAME="ch3b_ur"
  character (len = *), parameter :: CH4_UR_NAME="ch4_ur"
  character (len = *), parameter :: CH5_UR_NAME="ch5_ur"

  character (len = *), parameter :: CH3B_US_NAME="ch3b_us"
  character (len = *), parameter :: CH4_US_NAME="ch4_us"
  character (len = *), parameter :: CH5_US_NAME="ch5_us"

  character (len = *), parameter :: CH1_COUNTS_NAME="ch1_counts"
  character (len = *), parameter :: CH2_COUNTS_NAME="ch2_counts"
  character (len = *), parameter :: CH3_COUNTS_NAME="ch3_counts"
  character (len = *), parameter :: CH4_COUNTS_NAME="ch4_counts"
  character (len = *), parameter :: CH5_COUNTS_NAME="ch5_counts"

  character (len = *), parameter :: CH1_NAME="ch1_ref"
  character (len = *), parameter :: CH2_NAME="ch2_ref"
  character (len = *), parameter :: CH3A_NAME="ch3a_ref"

  character (len = *), parameter :: CH3B_NAME="ch3b_bt"
  character (len = *), parameter :: CH4_NAME="ch4_bt"
  character (len = *), parameter :: CH5_NAME="ch5_bt"

  character (len = *), parameter :: T_ICT_NAME="t_ict"

  integer                        :: lat_varid,lon_varid, time_varid,&
                                    scanline_varid, &
                                    satza_varid,solza_varid,relaz_varid, &
                                    ch1_counts_varid, ch2_counts_varid,ch3_counts_varid,&
                                    ch4_counts_varid,ch5_counts_varid, &
                                    ch1_varid, ch2_varid,ch3a_varid,&
                                    ch3b_varid, ch4_varid,ch5_varid, &
                                    ch3b_ur_varid, ch4_ur_varid,ch5_ur_varid, & 
                                    ch3b_us_varid, ch4_us_varid,ch5_us_varid, & 
                                    tict_varid
          
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: FillValue = "FillValue"
  character (len = *), parameter :: add_offset = "add_offset"
  character (len = *), parameter :: scale_factor = "scale_factor"
  character (len = *), parameter :: valid_min = "valid_min"
  character (len = *), parameter :: valid_max = "valid_max"
  character (len = *), parameter :: coordinates = "coordinates"
  character (len = *), parameter :: long_name = "long_name"
  character (len = *), parameter :: standard_name = "standard_name"
  character (len = *), parameter :: ancilliaty_variables = "ancilliary_variables"
  character (len = *), parameter :: reference_datum = "reference_datum"
  character (len = *), parameter :: calendar = "calendar"

  character (len = *), parameter :: lat_standard_name= "Latitude"
  character (len = *), parameter :: lat_UNITS = "degrees_north"
  real, parameter                :: lat_valid_min = -90.
  real, parameter                :: lat_valid_max = 90.
  character (len = *), parameter :: lat_reference_datum = "geographical coordinates, WGS84 projection"

  character (len = *), parameter :: lon_standard_name= "Longitude"
  character (len = *), parameter :: lon_UNITS = "degrees_east"
  real, parameter                :: lon_valid_min = 0.
  real, parameter                :: lon_valid_max = 360.
  character (len = *), parameter :: lon_reference_datum = "geographical coordinates, WGS84 projection"
  
  character (len = *), parameter :: time_standard_name= "time"
  character (len = *), parameter :: time_UNITS = "seconds since 1970-01-01 00:00:00"
  character (len = *), parameter :: time_calendar = "gregorian"

  character (len = *), parameter :: scanline_UNITS = "s"

  character (len = *), parameter :: solza_standard_name= "solZA"
  character (len = *), parameter :: satza_standard_name= "satZA"
  character (len = *), parameter :: relaz_standard_name= "relAZ"
  character (len = *), parameter :: angle_UNITS = "degrees"
  real, parameter                :: angle_add_offset = 0.0
  real, parameter                :: angle_scale_factor = 0.01
  
  character (len = *), parameter :: bt_UNITS = "kelvin"
  real, parameter                :: bt_add_offset = 273.15
  real, parameter                :: bt_scale_factor = 0.01
  integer, parameter             :: bt_valid_min = 0
  integer, parameter             :: bt_valid_max = 1000

  character (len = *), parameter :: u_UNITS = "kelvin"
  real, parameter                :: u_add_offset = 0
  real, parameter                :: u_scale_factor = 1
  integer, parameter             :: u_valid_min = 0
  integer, parameter             :: u_valid_max = 1000

!Long names
  character (len = *), parameter  :: lat_long_name = "Latitude coordinates"
  character (len = *), parameter  :: lon_long_name = "Longitude coordinates"
  character (len = *), parameter  :: time_long_name = "reference time of sst file"
  character (len = *), parameter  :: satza_long_name = "Satellite zenith angle"
  character (len = *), parameter  :: solza_long_name = "Solar zenith angle"
  character (len = *), parameter  :: relaz_long_name = "Relative satellite angle"

!Earth counts 
  character (len = *), parameter  :: ch1_counts_long_name = "Channel 1 Earth counts"
  character (len = *), parameter  :: ch2_counts_long_name = "Channel 2  Earth counts"
  character (len = *), parameter  :: ch3_counts_long_name = "Channel 3A  Earth counts"
  character (len = *), parameter  :: ch4_counts_long_name = "Channel 4  Earth counts"
  character (len = *), parameter  :: ch5_counts_long_name = "Channel 5  Earth counts"

!Reflectances
  character (len = *), parameter  :: ch1_long_name = "Channel 1 Reflectance"
  character (len = *), parameter  :: ch2_long_name = "Channel 2 Reflectance"
  character (len = *), parameter  :: ch3a_long_name = "Channel 3A Reflectance"

!BT
  character (len = *), parameter  :: ch3b_long_name = "Channel 3B Brightness Temperature"
  character (len = *), parameter  :: ch4_long_name = "Channel 4 Brightness Temperature"
  character (len = *), parameter  :: ch5_long_name = "Channel 5 Brightness Temperature"

!uncertaintiees
  character (len = *), parameter  :: ch3b_ur_long_name = "Random component of the uncertainty on ch.3B BT"
  character (len = *), parameter  :: ch4_ur_long_name = "Random component of the uncertainty on ch.4 BT"
  character (len = *), parameter  :: ch5_ur_long_name =  "Random component of the uncertainty on ch.5 BT"   

  character (len = *), parameter  :: ch3b_us_long_name = "Non random component of the uncertainty on ch.3B BT"
  character (len = *), parameter  :: ch4_us_long_name = "Non random component of the uncertainty on ch.4 BT"
  character (len = *), parameter  :: ch5_us_long_name = "Non random component of the uncertainty on ch.5 BT"     
! Loop indices
 ! integer :: lat, lon, time, i
print*, "writing easy FCDR"
print*, "lat",minval(AVHRR%lat), maxval(AVHRR%lat)
print*, "lon",minval(AVHRR%lon), maxval(AVHRR%lon)
print*, "bt3",minval(AVHRR%btf3), maxval(AVHRR%btf3)
print*, "bt4",minval(AVHRR%btf4), maxval(AVHRR%btf4)
print*, "bt5", minval(AVHRR%btf5), maxval(AVHRR%btf5)
print*, "time", minval(AVHRR%time),  minval(AVHRR%time)  
print*, "scan", minval(AVHRR%scanlinenumber), maxval(AVHRR%scanlinenumber)
print*, "satza", minval(AVHRR%satza), maxval(AVHRR%satza)
print*, "bt5", minval(AVHRR%solza) , maxval(AVHRR%solza) 
print*, "relaz", minval(AVHRR%relaz), maxval(AVHRR%relaz)
print*, "tict",AVHRR%orbital_temperature

print*, "C1", minval(AVHRR%Counts1), maxval(AVHRR%Counts1)
print*, "C2", minval(AVHRR%Counts2), maxval(AVHRR%Counts2)
print*, "c3", minval(AVHRR%Counts3), maxval(AVHRR%Counts3)
print*, "C4", minval(AVHRR%Counts4), maxval(AVHRR%Counts4)
print*, "C5", minval(AVHRR%Counts5), maxval(AVHRR%Counts5)
print*, "ur3", minval(AVHRR%ur3), maxval(AVHRR%ur3)
print*, "ur4", minval(AVHRR%ur4), maxval(AVHRR%ur4)
print*, "ur5", minval(AVHRR%ur5), maxval(AVHRR%ur5)
print*, "us3", minval(AVHRR%us3), maxval(AVHRR%us3)
print*, "us4", minval(AVHRR%us4), maxval(AVHRR%us4)
print*, "us5", minval(AVHRR%us5), maxval(AVHRR%us5)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create the file. 
!!!!!!!!!!!!!!!!!
  !call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )
 call check( nf90_create(FILE_NAME, cmode=IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid=ncid) )
  print*, file_name,"created"

! Define the dimensions. The record dimension is defined to have unlimited length 
!- it can grow as needed. In this example it is the time dimension.)
  !call check( nf90_def_dim(ncid, TIME_NAME, NTIME,time_dimid) )
  call check( nf90_def_dim(ncid, J_NAME, NJ, nj_dimid) )
  call check( nf90_def_dim(ncid, I_NAME, NI, ni_dimid) )
  call check( nf90_def_dim(ncid, TEN_NAME, N10, n10_dimid) )
  print*, "dimensions defined"

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! Both of the netCDF variables we are creating share the same four dimensions. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
  dimids_2 = (/ ni_dimid, nj_dimid/)
  dimids_sp_ict = (/ n10_dimid, nj_dimid/)
  dimids_1 = (/ nj_dimid/)
  print*, "dimids defined"
  chunk_size_pixel = (/ 409,1000 /)
  chunk_size_10_per_line = (/ 10, 1000/)
  chunk_size_1_per_line = (/ 1000/)
  chunk_size_3_per_line = (/3,1000/)
  chunk_size_4_per_line = (/4,1000/)

!-Define the netCDF variables 
!-Assign units attributes to the netCDF variables.
!--Lat
   call check( nf90_def_var(ncid, LAT_NAME, NF90_FLOAT, dimids_2, lat_varid, &
  chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_put_att(ncid, lat_varid, "_FillValue", real_fillvalue) )
   call check( nf90_put_att(ncid, lat_varid, long_name, lat_long_name) )
   call check( nf90_put_att(ncid, lat_varid, standard_name, lat_standard_name) )
   call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
   call check( nf90_put_att(ncid, lat_varid, valid_min, lat_valid_min) )
   call check( nf90_put_att(ncid, lat_varid, valid_max, lat_valid_max) )
   call check( nf90_put_att(ncid, lat_varid, reference_datum,lat_reference_datum)) 
print*, "Lat defined"
!--Lon
   call check( nf90_def_var(ncid, LON_NAME, NF90_FLOAT, dimids_2, lon_varid, &
   chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, lon_varid, "_FillValue", real_fillvalue) )
   call check( nf90_put_att(ncid, lon_varid, long_name, lon_long_name) )
   call check( nf90_put_att(ncid, lon_varid, standard_name, lon_standard_name) )
   call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )
   call check( nf90_put_att(ncid, lon_varid, valid_min, lon_valid_min) )
   call check( nf90_put_att(ncid, lon_varid, valid_max, lon_valid_max) )
   call check( nf90_put_att(ncid, lon_varid, reference_datum,lon_reference_datum)) 
print*, "Lon defined"
!--Time
   call check( nf90_def_var(ncid, TIME_NAME, NF90_FLOAT, dimids_1, time_varid, &
    chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, time_varid, long_name, time_long_name) )
   call check( nf90_put_att(ncid, time_varid, standard_name, time_standard_name) )
   call check( nf90_put_att(ncid, time_varid, UNITS, time_units) )
   call check( nf90_put_att(ncid, time_varid, calendar, time_calendar) )
print*, "time defined"
   call check( nf90_def_var(ncid, SCAN_LINE_NUMBER_NAME, NF90_INT, dimids_1, scanline_varid, &
    chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, scanline_varid, "_FillValue", integer_fillvalue) )
print*, "scan line defined"
   call check( nf90_def_var(ncid, SATZA_NAME, NF90_INT, dimids_2, satza_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, satza_varid, "_FillValue", integer_fillvalue) )
   !call check( nf90_put_att(ncid, satza_varid, scale_factor, 0.01) )
   !call check( nf90_put_att(ncid, satza_varid, add_offset, 0) )
   call check( nf90_put_att(ncid, satza_varid, long_name, satza_long_name) )
   call check( nf90_put_att(ncid, satza_varid, standard_name, satza_standard_name) )
   call check( nf90_put_att(ncid, satza_varid, UNITS, angle_units) )
   call check( nf90_put_att(ncid, satza_varid, valid_min, 0) )
   call check( nf90_put_att(ncid, satza_varid, valid_max,9000) )
print*, "satza defined"
   call check( nf90_def_var(ncid, SOLZA_NAME, NF90_INT, dimids_2, solza_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, solza_varid, "_FillValue", integer_fillvalue) )
   !call check( nf90_put_att(ncid, solza_varid, scale_factor, 0.01) )
   !call check( nf90_put_att(ncid, solza_varid, add_offset, 0) )
   call check( nf90_put_att(ncid, solza_varid, long_name, satza_long_name) )
   call check( nf90_put_att(ncid, solza_varid, standard_name, satza_standard_name) )
   call check( nf90_put_att(ncid, solza_varid, UNITS, angle_units) )
   call check( nf90_put_att(ncid, solza_varid, valid_min, 0) )
   call check( nf90_put_att(ncid, solza_varid, valid_max,18000) )
print*, "solza defined" 
   call check( nf90_def_var(ncid, RELAZ_NAME, NF90_FLOAT, dimids_2, relaz_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_put_att(ncid, relaz_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, relaz_varid, long_name, relaz_long_name) )
   call check( nf90_put_att(ncid, relaz_varid, standard_name, relaz_standard_name) )
   call check( nf90_put_att(ncid, relaz_varid, UNITS, angle_units) )
print*, "relaz defined"
!--Earth counts 
   call check( nf90_def_var(ncid, CH1_COUNTS_NAME, NF90_FLOAT,dimids_2 , ch1_counts_varid) )
   call check( nf90_def_var(ncid, CH2_COUNTS_NAME, NF90_FLOAT, dimids_2, ch2_counts_varid) )
  call check( nf90_def_var(ncid, CH3_COUNTS_NAME, NF90_FLOAT, dimids_2, ch3_counts_varid) )
   call check( nf90_def_var(ncid, CH4_COUNTS_NAME, NF90_FLOAT, dimids_2, ch4_counts_varid) )
   call check( nf90_def_var(ncid, CH5_COUNTS_NAME, NF90_FLOAT, dimids_2, ch5_counts_varid) )

   call check( nf90_put_att(ncid, ch1_counts_varid, "_FillValue", real_FillValue) )
   call check( nf90_put_att(ncid, ch1_counts_varid, long_name,ch1_counts_long_name)) 

  call check( nf90_put_att(ncid, ch2_counts_varid, "_FillValue", real_FillValue) )
   call check( nf90_put_att(ncid, ch2_counts_varid, long_name,ch2_counts_long_name))

  call check( nf90_put_att(ncid, ch3_counts_varid, "_FillValue", real_FillValue) )
   call check( nf90_put_att(ncid, ch3_counts_varid, long_name,ch3_counts_long_name))

   call check( nf90_put_att(ncid, ch4_counts_varid, "_FillValue", real_FillValue) )
   call check( nf90_put_att(ncid, ch4_counts_varid, long_name,ch4_counts_long_name))

  call check( nf90_put_att(ncid, ch5_counts_varid, "_FillValue", real_FillValue) )
   call check( nf90_put_att(ncid, ch5_counts_varid, long_name,ch5_counts_long_name))

!--Reflectances
   call check( nf90_def_var(ncid, CH1_NAME, NF90_INT,dimids_2 , ch1_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH2_NAME, NF90_INT, dimids_2, ch2_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH3A_NAME, NF90_INT, dimids_2, ch3a_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))


   call check( nf90_put_att(ncid, ch1_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch1_varid, long_name,ch1_long_name)) 
   call check( nf90_put_att(ncid, ch2_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch2_varid, long_name,ch2_long_name))

   call check( nf90_put_att(ncid, ch3a_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch3a_varid, long_name,ch3a_long_name))
print*, "reflectances defined"  
!--BT 
   call check( nf90_def_var(ncid, CH3B_NAME, NF90_INT,dimids_2 , ch3b_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
  call check( nf90_def_var(ncid, CH4_NAME, NF90_INT, dimids_2, ch4_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_NAME, NF90_INT, dimids_2, ch5_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   
   call check( nf90_put_att(ncid, ch3b_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch3b_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch3b_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch3b_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch3b_varid, long_name,ch3b_long_name)) 

   call check( nf90_put_att(ncid, ch4_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch4_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch4_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch4_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch4_varid, long_name,ch4_long_name))

   call check( nf90_put_att(ncid, ch5_varid, "_FillValue", integer_FillValue) )
   call check( nf90_put_att(ncid, ch5_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch5_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch5_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch5_varid, long_name,ch5_long_name))
print*, "bt defined"
!---Uncertainty on earth BT
   call check( nf90_def_var(ncid, CH3B_UR_NAME, NF90_FLOAT,dimids_2 , ch3b_ur_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_UR_NAME, NF90_FLOAT, dimids_2, ch4_ur_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_UR_NAME, NF90_FLOAT, dimids_2, ch5_ur_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))

   call check( nf90_def_var(ncid, CH3B_US_NAME, NF90_FLOAT,dimids_2 , ch3b_us_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_US_NAME, NF90_FLOAT, dimids_2, ch4_us_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_US_NAME, NF90_FLOAT, dimids_2, ch5_us_varid, &
    chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))

  call check( nf90_put_att(ncid, ch3b_ur_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch3b_ur_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch3b_ur_varid, long_name,ch3b_ur_long_name))

  call check( nf90_put_att(ncid, ch4_ur_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch4_ur_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch4_ur_varid, long_name,ch4_ur_long_name))

  call check( nf90_put_att(ncid, ch5_ur_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch5_ur_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch5_ur_varid, long_name,ch5_ur_long_name))

  call check( nf90_put_att(ncid, ch3b_us_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch3b_us_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch3b_us_varid, long_name,ch3b_us_long_name))

  call check( nf90_put_att(ncid, ch4_us_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch4_us_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch4_us_varid, long_name,ch4_us_long_name))

  call check( nf90_put_att(ncid, ch5_us_varid, "_FillValue", real_FillValue) )
  call check( nf90_put_att(ncid, ch5_us_varid, UNITS, u_UNITS) )
  call check( nf90_put_att(ncid, ch5_us_varid, long_name,ch5_us_long_name))
 
  call check( nf90_def_var(ncid, T_ICT_NAME, NF90_INT,tict_varid)) 
   ! chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
  call check( nf90_put_att(ncid, tict_varid, "_FillValue", integer_FillValue) )
  call check( nf90_put_att(ncid, tict_varid, UNITS, bt_UNITS) )
  !call check( nf90_put_att(ncid, tict_varid, scale_factor, bt_scale_factor) )
  !call check( nf90_put_att(ncid, tict_varid, add_offset, bt_add_offset) )
  call check( nf90_put_att(ncid, tict_varid, valid_min, bt_valid_min) )
  call check( nf90_put_att(ncid, tict_varid, valid_max, bt_valid_max) )


   print*, "units, fillvalue... attribues"
 
!--End define mode.
   call check( nf90_enddef(ncid) )
   print*, "fin de la definition"  
  
!--Write the pretend data. This will write our surface pressure and surface temperature data. 
!--The arrays only hold one timestep worth of data. We will just rewrite the same data for each timestep. 
!--In a real :: application, the data would change between timesteps.
  
  count = (/ NI,NJ/)
  count_10 = (/ N10,NJ/)
  count_1 = (/ NJ/)  
  start = (/ 1, 1/)
  start_1 = (/ 1/)
  start_10=(/ N10,NJ/) 

  print*, "lon", maxval(AVHRR%lon), minval(AVHRR%lon)
  call check( nf90_put_var(ncid, lon_varid,  AVHRR%lon, start =start, count = count) )
  print*, "lat"
  call check( nf90_put_var(ncid, lat_varid,  AVHRR%lat, start = start, count = count) )
  call check( nf90_put_var(ncid, time_varid,  AVHRR%time, start = start_1, count = count_1) )
  print*, "line"
  call check( nf90_put_var(ncid, scanline_varid,  AVHRR%scanlinenumber, start =start_1, count = count_1) )
  print*, "satza"

  call check( nf90_put_var(ncid, satza_varid,  AVHRR%satza*100, start = start, count = count) )
  print*, "solza"
  call check( nf90_put_var(ncid, solza_varid,  AVHRR%solza, start = start, count = count) )
  print*, "relaz"
  call check( nf90_put_var(ncid, relaz_varid,  AVHRR%relaz, start = start, count = count) )
  print*, "tict", AVHRR%orbital_temperature
  !call check( nf90_put_var(ncid, tict_varid,  AVHRR%orbital_temperature))

  print*, "Earth counts"
  call check( nf90_put_var(ncid, ch1_counts_varid, AVHRR%Counts1, start =start, count = count) )
  call check( nf90_put_var(ncid, ch2_counts_varid,  AVHRR%Counts2, start = start, count = count) )
  call check( nf90_put_var(ncid, ch3_counts_varid, AVHRR%Counts3, start =start, count = count) )
  call check( nf90_put_var(ncid, ch4_counts_varid,  AVHRR%Counts4, start = start, count = count) )
  call check( nf90_put_var(ncid, ch5_counts_varid,  AVHRR%Counts5, start = start, count = count) )
  
  print*, "bt"

 call check( nf90_put_var(ncid, ch3b_varid, AVHRR%btf3, start =start, count = count) )
  call check( nf90_put_var(ncid, ch4_varid,  AVHRR%btf4, start = start, count = count) )
  call check( nf90_put_var(ncid, ch5_varid,  AVHRR%btf5, start = start, count = count) )

  !print*, "reflectances"
  !call check( nf90_put_var(ncid, ch1_varid, anint(AVHRR%array1), start =start, count = count) )
  !call check( nf90_put_var(ncid, ch2_varid,  AVHRR%array2, start = start, count = count) )
  !call check( nf90_put_var(ncid, ch3a_varid,  AVHRR%array3a, start = start, count = count) )
 
 
  print*, "Ne3"
  call check( nf90_put_var(ncid, ch3b_ur_varid, AVHRR%ur3, start =start, count = count) )
  call check( nf90_put_var(ncid, ch4_ur_varid,  AVHRR%ur4, start = start, count = count) )
  call check( nf90_put_var(ncid, ch5_ur_varid,  AVHRR%ur5, start = start, count = count) ) 
 
  call check( nf90_put_var(ncid, ch3b_us_varid, AVHRR%us3, start =start, count = count) )
  call check( nf90_put_var(ncid, ch4_us_varid,  AVHRR%us4, start = start, count = count) )
  call check( nf90_put_var(ncid, ch5_us_varid,  AVHRR%us5, start = start, count = count) ) 

         
!-Close the netcdf file. 
call check( nf90_close(ncid) )
write(*,*) "file netcdf close"    
end subroutine write_FCDR

!***************************************************
subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  
!********************************************

end module module_FCDR
