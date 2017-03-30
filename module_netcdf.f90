module module_netcdf

  use netcdf
  use type

implicit none
PUBLIC

CONTAINS

SUBROUTINE write_netcdf(FILE_NAME,AVHRR)
  
  character (len = *), INTENT(IN)         :: FILE_NAME
  TYPE(AVHRR_Data),INTENT(IN)             :: AVHRR

  integer                                 :: ncid

  integer, parameter                      :: NDIMS = 2 
  integer,parameter                      :: NI = 409, N10=10, N1=1,N3=3,N4=4
  
  character (len = *), parameter          :: I_NAME = "ni",&
                                             J_NAME = "nj", &
                                             TEN_NAME = "n10", &
                                             ONE_NAME = "n1", &
                                             THREE_NAME = "n3", &
                                             FOUR_NAME = "n4"
   
  integer                                 :: ni_dimid, nj_dimid,n10_dimid, &
                                             n1_dimid,n3_dimid,n4_dimid

  integer,dimension(2)                    ::   start_2 = (/ 1, 1/)
                                             !count_pixel = (/ NI,NF90_UNLIMITED/), &
                                             !count_10 = (/ N10,NF90_UNLIMITED/), &
                                             !count_3=(/3,NF90_UNLIMITED/)
                                                                      
  integer,dimension(1)                    :: count_nspace=(/3/), & 
                                             start_nspace=(/1/), &
                                             start_1=(/1/), count_1_value=(/1/)
                                             !count_1 = (/NF90_UNLIMITED/), &
                                   
  integer                                 :: dimids_1_per_line(1),&
                                             dimids_10_per_line(2),&
                                             dimids_pixel(2),  &
                                             dimids_3_per_line(2),&
                                             dimids_4_per_line(2), &
                                             dimids_nspace(1)
  
  integer                                 :: chunk_size_1_per_line(1),&
                                             chunk_size_10_per_line(2),&
                                             chunk_size_pixel(2),  &
                                             chunk_size_3_per_line(2),&
                                             chunk_size_4_per_line(2)
  
  character (len = *), parameter          :: SCAN_LINE_NAME="scanline", &
                                             LAT_NAME="Lat", &
                                             LON_NAME="Lon", &
                                             SATZA_NAME="SatZa", &
                                             SOLZA_NAME="SolZa", &
                                             TIME_NAME = "time", &
                                             YEAR_NAME = "year",&
                                             MONTH_NAME = "month",&
                                             DAY_NAME = "day",&
                                             DAY_NO_NAME = "day_no",&
                                             HOUR_NAME = "hour",&
                                             UTC_NAME = "utc",&
                                             CH3B_COUNTS_NAME="ch3b_counts",&
                                             CH4_COUNTS_NAME="ch4_counts",&
                                             CH5_COUNTS_NAME="ch5_counts",&
                                             CH3B_SPACE_COUNTS_NAME="ch3b_space_counts",&
                                             CH4_SPACE_COUNTS_NAME="ch4_space_counts",&
                                             CH5_SPACE_COUNTS_NAME="ch5_space_counts",&
                                             CH3B_ICT_COUNTS_NAME="ch3b_ict_counts",&
                                             CH4_ICT_COUNTS_NAME="ch4_ict_counts",&
                                             CH5_ICT_COUNTS_NAME="ch5_ict_counts",&
                                             CH3B_NE_NAME="ch3b_ne",&
                                             CH4_NE_NAME="ch4_ne",&
                                             CH5_NE_NAME="ch5_ne",&
                                             CH3B_BT_NAME="ch3b_bt",&
                                             CH4_BT_NAME="ch4_bt",&
                                             CH5_BT_NAME="ch5_bt",&
                                             CH3B_RICT_NAME="ch3b_rict",&
                                             CH4_RICT_NAME="ch4_rict",&
                                             CH5_RICT_NAME="ch5_rict",&
                                             CH3B_GAIN_NAME="ch3b_gain",&
                                             CH4_GAIN_NAME="ch4_gain",&
                                             CH5_GAIN_NAME="ch5_gain",&
                                             NSPACE_NAME="nspace",&
                                             CH3B_MEAN_SP_NAME="ch3b_mean_space",&
                                             CH4_MEAN_SP_NAME="ch4_mean_space",&
                                             CH5_MEAN_SP_NAME="ch5_mean_space",&
                                             CH3B_MEAN_BB_NAME="ch3b_mean_ict",&
                                             CH4_MEAN_BB_NAME="ch4_mean_ict",&
                                             CH5_MEAN_BB_NAME="ch5_mean_ict",&
                                             CH3B_RAMP_NAME="ch3b_ramp",&
                                             CH4_RAMP_NAME="ch4_ramp",&
                                             CH5_RAMP_NAME="ch5_ramp",&
                                             CH3B_COEF_CALIB_NAME="ch3b_coef_calib",&
                                             CH4_COEF_CALIB_NAME="ch4_coef_calib",&
                                             CH5_COEF_CALIB_NAME="ch5_coef_calib",&
                                             PRT_TEMP_1_NAME="prt_temp_1", &
                                             PRT_TEMP_2_NAME="prt_temp_2", &
                                             PRT_TEMP_3_NAME="prt_temp_3", &
                                             PRT_TEMP_4_NAME="prt_temp_4",&
                                             PRT_MEAN_TEMP_NAME="prt_mean_temp",&
                                             PRT_SIGMA_TEMP_NAME="prt_sigma_temp",& 
                                             PRT_COUNTS_NAME="prt_counts",&
                                              PRT_NUMBER_NAME="prt_number", &        
                                             RADIATOR_NAME="radiator",&
                                             ELECTRONICS_NAME="electronics",&
                                             COOLER_NAME="cooler",&
                                             BASEPLATE_NAME="baseplate",&
                                             MOTOR_NAME="motor",&
                                             A_D_CONV_NAME="a_d_conv",&
                                             PATCH_NAME="patch",&
                                             PATCH_EXTENDED_NAME="patch_extended",&  
                                             CH3B_MEAN_BB_ORBITE_NAME="ch3b_mean_bb_orbite",&
                                             CH3B_MEAN_SP_ORBITE_NAME="ch3b_mean_sp_orbite",&
                                             CH3B_MEAN_RICT_ORBITE_NAME="ch3b_mean_rict_orbite",&
                                             CH4_MEAN_BB_ORBITE_NAME="ch4_mean_bb_orbite",&
                                             CH4_MEAN_SP_ORBITE_NAME="ch4_mean_sp_orbite",&
                                             CH4_MEAN_RICT_ORBITE_NAME="ch4_mean_rict_orbite",&
                                             CH5_MEAN_BB_ORBITE_NAME="ch5_mean_bb_orbite",&
                                             CH5_MEAN_SP_ORBITE_NAME="ch5_mean_sp_orbite",&
                                             CH5_MEAN_RICT_ORBITE_NAME="ch5_mean_rict_orbite",&
                                             MEAN_PRT_ORBITE_NAME="mean_prt_orbite",&
                                             MEAN_PATCH_ORBITE_NAME="mean_patch_orbite",&
                                             MEAN_PATCH_EXTENDED_ORBITE_NAME="mean_patch_extended_orbite"  

!--ID par variable
  integer                                  :: scan_line_varid, lat_varid,lon_varid,&
                                              satza_varid,solza_varid, &
                                              time_varid ,year_varid,month_varid,&
                                              day_varid, day_no_varid,hour_varid,utc_varid, &
                                              ch3b_counts_varid, ch4_counts_varid,ch5_counts_varid, &
                                              ch3b_space_counts_varid, ch4_space_counts_varid,ch5_space_counts_varid, &
                                              ch3b_ict_counts_varid, ch4_ict_counts_varid,ch5_ict_counts_varid, &
                                              ch3b_ne_varid, ch4_ne_varid,ch5_ne_varid, &
                                              ch3b_bt_varid, ch4_bt_varid,ch5_bt_varid,&
                                              ch3b_rict_varid, ch4_rict_varid,ch5_rict_varid, &
                                              ch3b_gain_varid, ch4_gain_varid, ch5_gain_varid, &
                                              nspace_varid,&
                                              ch3b_mean_space_varid, ch4_mean_space_varid, ch5_mean_space_varid, &
                                              ch3b_mean_ict_varid, ch4_mean_ict_varid, ch5_mean_ict_varid, &  
                                              ch3b_coef_calib_varid, ch4_coef_calib_varid,ch5_coef_calib_varid, &
                                              ch3b_ramp_varid, ch4_ramp_varid,ch5_ramp_varid, &
                                              prt_temp_1_varid,prt_temp_2_varid, &
                                              prt_temp_3_varid, prt_temp_4_varid,&
                                              prt_mean_temp_varid, prt_sigma_temp_varid, &
                                              prt_counts_varid, prt_number_varid, &
                                              radiator_varid, electronics_varid, cooler_varid, baseplate_varid, &
                                              motor_varid, a_d_conv_varid, patch_varid,patch_extended_varid, &
                                              ch3b_mean_bb_orbite_varid, &
                                              ch3b_mean_sp_orbite_varid, &
                                              ch3b_mean_rict_orbite_varid,&
                                              ch4_mean_bb_orbite_varid, &
                                              ch4_mean_sp_orbite_varid, &
                                              ch4_mean_rict_orbite_varid,&
                                              ch5_mean_bb_orbite_varid, &
                                              ch5_mean_sp_orbite_varid, &
                                              ch5_mean_rict_orbite_varid,&
                                              mean_prt_orbite_varid, &
                                              mean_patch_orbite_varid, &
                                              mean_patch_extended_orbite_varid 
                                  
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

  real, parameter                :: real_FillValue = -1.00000002e+30
  integer, parameter             :: counts_FillValue = -32768

  character (len = *), parameter :: lat_standard_name= "Latitude"
  character (len = *), parameter :: lat_UNITS = "degrees_north"
  real, parameter                :: lat_valid_min = -90.
  real, parameter                :: lat_valid_max = 90.
  character (len = *), parameter :: lat_reference_datum = "geographical coordinates, WGS84 projection"

 
  character (len = *), parameter :: lon_standard_name= "Longitude"
  character (len = *), parameter :: lon_UNITS = "degrees_east"
  real, parameter                :: lon_valid_min = -180.
  real, parameter                :: lon_valid_max = 180.
  character (len = *), parameter :: lon_reference_datum = "geographical coordinates, WGS84 projection"
  
 
  character (len = *), parameter :: time_standard_name= "time"
  character (len = *), parameter :: time_UNITS = "seconds since 1981-01-01 00:00:00"
  character (len = *), parameter :: time_calendar = "gregorian"
  

  character (len = *), parameter :: bt_UNITS = "kelvin"
  integer, parameter             :: bt_valid_min = 0
  integer, parameter             :: bt_valid_max = 1000

  character (len = *), parameter :: ne_UNITS = "W.sr-1.m-2"
  real, parameter                :: ne_add_offset = 0
  real, parameter                :: ne_scale_factor = 1
  integer, parameter             :: ne_valid_min = 0
  integer, parameter             :: ne_valid_max = 1000

 
  character (len = *), parameter :: counts_UNITS = "counts"
  integer, parameter             :: counts_valid_min = 0
  integer, parameter             :: counts_valid_max = 1500


!Long names
  character (len = *), parameter  :: lat_long_name = "Latitude coordinates"
  character (len = *), parameter  :: lon_long_name = "Longitude coordinates"
  character (len = *), parameter  :: time_long_name = "reference time of sst file"
!BT
  character (len = *), parameter  :: ch3b_long_name = "Channel 3b Brightness Temperature"
  character (len = *), parameter  :: ch4_long_name = "Channel 4 Brightness Temperature"
  character (len = *), parameter  :: ch5_long_name = "Channel 5 Brightness Temperature"
!radiance
  character (len = *), parameter  :: ch3b_ne_long_name = "Channel 3B Earth radiance"
  character (len = *), parameter  :: ch4_ne_long_name = "Channel 4 Earth radiance"
  character (len = *), parameter  :: ch5_ne_long_name = "Channel 5 Earth radiance"

  character (len = *), parameter  :: ch3b_rict_long_name = "Channel 3B ICT Radiance"
  character (len = *), parameter  :: ch4_rict_long_name = "Channel 4 ICT Radiance"
  character (len = *), parameter  :: ch5_rict_long_name = "Channel 5 ICT Radiance"
!counts
  character (len = *), parameter  :: ch3b_ce_long_name = "Channel 3B Earth counts"
  character (len = *), parameter  :: ch4_ce_long_name = "Channel 4 Earth counts"
  character (len = *), parameter  :: ch5_ce_long_name = "Channel 5 Earth counts"

  character (len = *), parameter  :: ch3b_cict_long_name = "Channel 3B ICT counts"
  character (len = *), parameter  :: ch4_cict_long_name = "Channel 4 ICT counts"
  character (len = *), parameter  :: ch5_cict_long_name = "Channel 5 ICT counts"

  character (len = *), parameter  :: ch3b_csp_long_name = "Channel 3B Space counts"
  character (len = *), parameter  :: ch4_csp_long_name = "Channel 4 Space counts"
  character (len = *), parameter  :: ch5_csp_long_name = "Channel 5 Space counts"

 
  ! Loop indices
  integer :: lat, lon, time, i

!-2-!
!-Create the file. 
  call check( nf90_create(FILE_NAME, cmode=IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid=ncid) )
  !print*, "file created"

!-3-
!-Define the dimensions. The record dimension is defined to have unlimited length 
!- it can grow as needed. In this example it is the time dimension.)
  !call check( nf90_def_dim(ncid, TIME_NAME, NTIME,time_dimid) )
  call check( nf90_def_dim(ncid, J_NAME, NF90_UNLIMITED, nj_dimid) )
  call check( nf90_def_dim(ncid, I_NAME, NI, ni_dimid) )
  call check( nf90_def_dim(ncid, TEN_NAME, N10, n10_dimid) )
  call check( nf90_def_dim(ncid, THREE_NAME, N3, n3_dimid) )
  call check( nf90_def_dim(ncid, ONE_NAME, N1, n1_dimid) )
  call check( nf90_def_dim(ncid, FOUR_NAME, N4, n4_dimid) )
  !print*, "dimensions defined"
  
! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! Both of the netCDF variables we are creating share the same four dimensions. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
  dimids_pixel = (/ ni_dimid, nj_dimid/)
  dimids_10_per_line = (/ n10_dimid, nj_dimid/)
  dimids_1_per_line = (/ nj_dimid/)
  dimids_3_per_line = (/n3_dimid,nj_dimid/)
  dimids_4_per_line = (/n4_dimid,nj_dimid/)
  dimids_nspace=(/n3_dimid/)

  chunk_size_pixel = (/ 409,1000 /)
  chunk_size_10_per_line = (/ 10, 1000/)
  chunk_size_1_per_line = (/ 1000/)
  chunk_size_3_per_line = (/3,1000/)
  chunk_size_4_per_line = (/4,1000/)
!-4-
!-Define the netCDF variables 
!--Lat, Lon, time
   !print*, "lat lon time"
   call check( nf90_def_var(ncid, SCAN_LINE_NAME, NF90_INT, dimids_1_per_line, scan_line_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, LAT_NAME, NF90_FLOAT, dimids_pixel, lat_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, LON_NAME, NF90_FLOAT, dimids_pixel, lon_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   
   call check( nf90_def_var(ncid, YEAR_NAME, NF90_INT, dimids_1_per_line, year_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, MONTH_NAME, NF90_INT, dimids_1_per_line, month_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, DAY_NAME, NF90_INT, dimids_1_per_line, day_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, DAY_NO_NAME, NF90_INT, dimids_1_per_line, day_no_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, HOUR_NAME, NF90_FLOAT, dimids_1_per_line, hour_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, UTC_NAME, NF90_INT, dimids_1_per_line, utc_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, dimids_1_per_line, time_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   
   call check( nf90_def_var(ncid, SATZA_NAME, NF90_FLOAT, dimids_pixel, satza_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, SOLZA_NAME, NF90_FLOAT, dimids_pixel, solza_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9)) 

!--EARTH BT 
   !print*, "bt"  
   call check( nf90_def_var(ncid, CH3B_BT_NAME, NF90_FLOAT,dimids_pixel , ch3b_bt_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_BT_NAME, NF90_FLOAT, dimids_pixel, ch4_bt_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_BT_NAME, NF90_FLOAT, dimids_pixel, ch5_bt_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
!--Earth Radainces
  !print*, "radiances"
   call check( nf90_def_var(ncid, CH3B_NE_NAME, NF90_FLOAT,dimids_pixel , ch3b_ne_varid, &
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_NE_NAME, NF90_FLOAT, dimids_pixel, ch4_ne_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_NE_NAME, NF90_FLOAT, dimids_pixel, ch5_ne_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))

!--Earth Counts  
!print*, "counts"
   call check( nf90_def_var(ncid, CH3B_COUNTS_NAME, NF90_FLOAT,dimids_pixel , ch3b_counts_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_COUNTS_NAME, NF90_FLOAT, dimids_pixel, ch4_counts_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_COUNTS_NAME, NF90_FLOAT, dimids_pixel, ch5_counts_varid,&
        chunksizes=chunk_size_pixel, shuffle=.TRUE., deflate_level=9))
!--SPACE Counts
!print*, "sp"  
   call check( nf90_def_var(ncid, CH3B_SPACE_COUNTS_NAME,NF90_FLOAT, dimids_10_per_line, ch3b_space_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_SPACE_COUNTS_NAME, NF90_FLOAT, dimids_10_per_line, ch4_space_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_SPACE_COUNTS_NAME, NF90_FLOAT, dimids_10_per_line, ch5_space_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))
!--ICT Counts
!print*, "ict"  
   call check( nf90_def_var(ncid, CH3B_ICT_COUNTS_NAME, NF90_FLOAT,dimids_10_per_line , ch3b_ict_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_ICT_COUNTS_NAME, NF90_FLOAT, dimids_10_per_line, ch4_ict_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_ICT_COUNTS_NAME, NF90_FLOAT, dimids_10_per_line, ch5_ict_counts_varid,&
        chunksizes=chunk_size_10_per_line, shuffle=.TRUE., deflate_level=9))

!--Mean SP Counts
!print*, "mean sp"
   call check( nf90_def_var(ncid,CH3B_MEAN_SP_NAME ,NF90_FLOAT, dimids_1_per_line, ch3b_mean_space_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,CH4_MEAN_SP_NAME, NF90_FLOAT, dimids_1_per_line, ch4_mean_space_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,CH5_MEAN_SP_NAME, NF90_FLOAT, dimids_1_per_line, ch5_mean_space_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))

!--Mean ICT counts
!print*, "mean ict"
   call check( nf90_def_var(ncid,CH3B_MEAN_BB_NAME, NF90_FLOAT, dimids_1_per_line, ch3b_mean_ict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,CH4_MEAN_BB_NAME, NF90_FLOAT, dimids_1_per_line, ch4_mean_ict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,CH5_MEAN_BB_NAME, NF90_FLOAT, dimids_1_per_line, ch5_mean_ict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 


!--ICT Radainces
!print*, "ict radiance"
   call check( nf90_def_var(ncid, CH3B_RICT_NAME, NF90_FLOAT,dimids_1_per_line , ch3b_rict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH4_RICT_NAME, NF90_FLOAT, dimids_1_per_line, ch4_rict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid, CH5_RICT_NAME, NF90_FLOAT, dimids_1_per_line, ch5_rict_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))

!--Gain
!print*, "gain"
   call check( nf90_def_var(ncid, CH3B_GAIN_NAME, NF90_FLOAT, dimids_1_per_line, ch3b_gain_varid)) 
   call check( nf90_def_var(ncid,  CH4_GAIN_NAME, NF90_FLOAT, dimids_1_per_line, ch4_gain_varid)) 
   call check( nf90_def_var(ncid,   CH5_GAIN_NAME, NF90_FLOAT, dimids_1_per_line, ch5_gain_varid)) 
!--Nspace
!print*, "nspace"
    call check( nf90_def_var(ncid, NSPACE_NAME, NF90_FLOAT,dimids_nspace, nspace_varid)) 
!--Ramp
!print*, "ramp"
   call check( nf90_def_var(ncid,CH3B_RAMP_NAME, NF90_FLOAT, dimids_1_per_line, ch3b_ramp_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,CH4_RAMP_NAME, NF90_FLOAT, dimids_1_per_line, ch4_ramp_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,CH5_RAMP_NAME, NF90_FLOAT, dimids_1_per_line, ch5_ramp_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))

!--Coef Calibration
!print*, "coef calib"
!   call check( nf90_def_var(ncid, CH3B_COEF_CALIB_NAME,NF90_FLOAT, dimids_3_per_line, ch3b_coef_calib_varid,&
!        chunksizes=chunk_size_3_per_line, shuffle=.TRUE., deflate_level=9)) 
!   call check( nf90_def_var(ncid, CH4_COEF_CALIB_NAME,NF90_FLOAT, dimids_3_per_line, ch4_coef_calib_varid,&
!        chunksizes=chunk_size_3_per_line, shuffle=.TRUE., deflate_level=9)) 
!   call check( nf90_def_var(ncid, CH5_COEF_CALIB_NAME,NF90_FLOAT, dimids_3_per_line, ch5_coef_calib_varid,&
!        chunksizes=chunk_size_3_per_line, shuffle=.TRUE., deflate_level=9)) 

!--PRT TEMPERATURE
 !  print*, "prt temp"
   call check( nf90_def_var(ncid,PRT_TEMP_1_NAME,NF90_FLOAT, dimids_1_per_line, prt_temp_1_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,PRT_TEMP_2_NAME,NF90_FLOAT, dimids_1_per_line, prt_temp_2_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,PRT_TEMP_3_NAME,NF90_FLOAT, dimids_1_per_line, prt_temp_3_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,PRT_TEMP_4_NAME,NF90_FLOAT, dimids_1_per_line, prt_temp_4_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,PRT_MEAN_TEMP_NAME,NF90_FLOAT, dimids_1_per_line, prt_mean_temp_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,PRT_SIGMA_TEMP_NAME,NF90_FLOAT, dimids_1_per_line, prt_sigma_temp_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))  

!--PRT COUNTS
!print*, "prt counts"
   call check( nf90_def_var(ncid,  PRT_COUNTS_NAME,NF90_INT,dimids_3_per_line, prt_counts_varid,&
        chunksizes=chunk_size_3_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,PRT_NUMBER_NAME,NF90_INT, dimids_1_per_line, prt_number_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))      

!--TEMPERATURES DIVERSES
!print*, "temp divers"
   call check( nf90_def_var(ncid, RADIATOR_NAME,NF90_FLOAT, dimids_1_per_line, radiator_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,ELECTRONICS_NAME, NF90_FLOAT, dimids_1_per_line, electronics_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))
   call check( nf90_def_var(ncid,COOLER_NAME, NF90_FLOAT, dimids_1_per_line, cooler_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, BASEPLATE_NAME, NF90_FLOAT, dimids_1_per_line, baseplate_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,MOTOR_NAME, NF90_FLOAT, dimids_1_per_line, motor_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,A_D_CONV_NAME, NF90_FLOAT, dimids_1_per_line, a_d_conv_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid, PATCH_NAME, NF90_FLOAT, dimids_1_per_line, patch_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9)) 
   call check( nf90_def_var(ncid,PATCH_EXTENDED_NAME, NF90_FLOAT, dimids_1_per_line, patch_extended_varid,&
        chunksizes=chunk_size_1_per_line, shuffle=.TRUE., deflate_level=9))  

!--MOYENNES ORBITALES DIVERSES
!print*, "moyennes orbitales diverses"
   call check( nf90_def_var(ncid, CH3B_MEAN_SP_ORBITE_NAME, NF90_FLOAT,ch3b_mean_sp_orbite_varid))
   call check( nf90_def_var(ncid, CH4_MEAN_SP_ORBITE_NAME, NF90_FLOAT, ch4_mean_sp_orbite_varid)) 
   call check( nf90_def_var(ncid, CH5_MEAN_SP_ORBITE_NAME, NF90_FLOAT, ch5_mean_sp_orbite_varid)) 
   call check( nf90_def_var(ncid, CH3B_MEAN_BB_ORBITE_NAME, NF90_FLOAT, ch3b_mean_bb_orbite_varid))
   call check( nf90_def_var(ncid, CH4_MEAN_BB_ORBITE_NAME, NF90_FLOAT, ch4_mean_bb_orbite_varid))
   call check( nf90_def_var(ncid, CH5_MEAN_BB_ORBITE_NAME, NF90_FLOAT, ch5_mean_bb_orbite_varid))
   call check( nf90_def_var(ncid, CH3B_MEAN_RICT_ORBITE_NAME, NF90_FLOAT, ch3b_mean_rict_orbite_varid ))
   call check( nf90_def_var(ncid, CH4_MEAN_RICT_ORBITE_NAME, NF90_FLOAT, ch4_mean_rict_orbite_varid ))
   call check( nf90_def_var(ncid, CH5_MEAN_RICT_ORBITE_NAME, NF90_FLOAT, ch5_mean_rict_orbite_varid ))
   call check( nf90_def_var(ncid, MEAN_PRT_ORBITE_NAME, NF90_FLOAT, mean_prt_orbite_varid))
   call check( nf90_def_var(ncid, MEAN_PATCH_ORBITE_NAME,NF90_FLOAT, mean_patch_orbite_varid))
   call check( nf90_def_var(ncid, MEAN_PATCH_EXTENDED_ORBITE_NAME,NF90_FLOAT, mean_patch_extended_orbite_varid))


!*********************************************
!--5-
!-Assign units attributes to the netCDF variables.
!***********************************************
!--Lat
   call check( nf90_put_att(ncid, lat_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, lat_varid, long_name, lat_long_name) )
   call check( nf90_put_att(ncid, lat_varid, standard_name, lat_standard_name) )
   call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
   call check( nf90_put_att(ncid, lat_varid, valid_min, lat_valid_min) )
   call check( nf90_put_att(ncid, lat_varid, valid_max, lat_valid_max) )
   call check( nf90_put_att(ncid, lat_varid, reference_datum,lat_reference_datum)) 
!--Lon
   call check( nf90_put_att(ncid, lon_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, lon_varid, long_name, lon_long_name) )
   call check( nf90_put_att(ncid, lon_varid, standard_name, lon_standard_name) )
   call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )
   call check( nf90_put_att(ncid, lon_varid, valid_min, lon_valid_min) )
   call check( nf90_put_att(ncid, lon_varid, valid_max, lon_valid_max) )
   call check( nf90_put_att(ncid, lon_varid, reference_datum,lon_reference_datum)) 
!--Time
   call check( nf90_put_att(ncid, time_varid, long_name, time_long_name) )
   call check( nf90_put_att(ncid, time_varid, standard_name, time_standard_name) )
   call check( nf90_put_att(ncid, time_varid, UNITS, time_units) )
   call check( nf90_put_att(ncid, time_varid, calendar, time_calendar) )

!--Counts
   call check( nf90_put_att(ncid, ch3b_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch3b_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch3b_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch3b_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch3b_counts_varid, long_name,ch3b_ce_long_name))

   call check( nf90_put_att(ncid, ch4_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch4_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch4_counts_varid,valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch4_counts_varid,valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch4_counts_varid, long_name,ch4_ce_long_name))
  
   call check( nf90_put_att(ncid, ch5_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch5_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch5_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch5_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch5_counts_varid, long_name,ch5_ce_long_name))
!--space Counts
   call check( nf90_put_att(ncid, ch3b_space_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch3b_space_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch3b_space_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch3b_space_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch3b_space_counts_varid, long_name,ch3b_csp_long_name))

   call check( nf90_put_att(ncid, ch4_space_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch4_space_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch4_space_counts_varid,valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch4_space_counts_varid,valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch4_space_counts_varid, long_name,ch4_csp_long_name))
  
   call check( nf90_put_att(ncid, ch5_space_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch5_space_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch5_space_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch5_space_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch5_space_counts_varid, long_name,ch5_csp_long_name))

!--ICT Counts
   call check( nf90_put_att(ncid, ch3b_ict_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch3b_ict_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch3b_ict_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch3b_ict_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch3b_ict_counts_varid, long_name,ch3b_cict_long_name))

   call check( nf90_put_att(ncid, ch4_ict_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch4_ict_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch4_ict_counts_varid,valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch4_ict_counts_varid,valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch4_ict_counts_varid, long_name,ch4_cict_long_name))
  
   call check( nf90_put_att(ncid, ch5_ict_counts_varid, FillValue, counts_FillValue) )
   call check( nf90_put_att(ncid, ch5_ict_counts_varid, UNITS, counts_UNITS) )
   call check( nf90_put_att(ncid, ch5_ict_counts_varid, valid_min, counts_valid_min) )
   call check( nf90_put_att(ncid, ch5_ict_counts_varid, valid_max, counts_valid_max) )
   call check( nf90_put_att(ncid, ch5_ict_counts_varid, long_name,ch5_cict_long_name))
!--Radainces
   call check( nf90_put_att(ncid, ch3b_rict_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch3b_rict_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch3b_rict_varid, long_name,ch3b_rict_long_name))

   call check( nf90_put_att(ncid, ch4_rict_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch4_rict_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch4_rict_varid, long_name,ch4_rict_long_name))
  
   call check( nf90_put_att(ncid, ch5_rict_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch5_rict_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch5_rict_varid, long_name,ch5_rict_long_name))

   call check( nf90_put_att(ncid, ch3b_ne_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch3b_ne_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch3b_ne_varid, valid_min, ne_valid_min) )
   call check( nf90_put_att(ncid, ch3b_ne_varid, valid_max, ne_valid_max) )
   call check( nf90_put_att(ncid, ch3b_ne_varid, long_name,ch3b_ne_long_name))

   call check( nf90_put_att(ncid, ch4_ne_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch4_ne_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch4_ne_varid,valid_min, ne_valid_min) )
   call check( nf90_put_att(ncid, ch4_ne_varid,valid_max, ne_valid_max) )
   call check( nf90_put_att(ncid, ch4_ne_varid, long_name,ch4_ne_long_name))
  
   call check( nf90_put_att(ncid, ch5_ne_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch5_ne_varid, UNITS, ne_UNITS) )
   call check( nf90_put_att(ncid, ch5_ne_varid, valid_min, ne_valid_min) )
   call check( nf90_put_att(ncid, ch5_ne_varid, valid_max, ne_valid_max) )
   call check( nf90_put_att(ncid, ch5_ne_varid, long_name,ch5_ne_long_name))
!--BT
   call check( nf90_put_att(ncid, ch3b_bt_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch3b_bt_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch3b_bt_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch3b_bt_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch3b_bt_varid, long_name,ch3b_long_name)) 

   call check( nf90_put_att(ncid, ch4_bt_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch4_bt_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch4_bt_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch4_bt_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch4_bt_varid, long_name,ch4_long_name))

   call check( nf90_put_att(ncid, ch5_bt_varid, FillValue, real_FillValue) )
   call check( nf90_put_att(ncid, ch5_bt_varid, UNITS, bt_UNITS) )
   call check( nf90_put_att(ncid, ch5_bt_varid, valid_min, bt_valid_min) )
   call check( nf90_put_att(ncid, ch5_bt_varid, valid_max, bt_valid_max) )
   call check( nf90_put_att(ncid, ch5_bt_varid, long_name,ch5_long_name))

 !  print*, "units, fillvalue... attribues"
 
!--End define mode.
   call check( nf90_enddef(ncid) )
!   print*, "fin de la definition"  

!-6-  
 ! print*, "ecriture data"
 ! print*, "lon"
  call check( nf90_put_var(ncid, lon_varid,  AVHRR%lon) )
 ! print*, "lat"
  call check( nf90_put_var(ncid, lat_varid,  AVHRR%lat) )
  call check( nf90_put_var(ncid, scan_line_varid,  AVHRR%scanlinenumber) )
  call check( nf90_put_var(ncid, satza_varid,  AVHRR%satZA) )
 ! print*, "lat"
  call check( nf90_put_var(ncid, solza_varid,  AVHRR%solZA) )
!  print*, "time"

 
  call check( nf90_put_var(ncid, year_varid,  AVHRR%year))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, month_varid,  AVHRR%month))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, day_varid,  AVHRR%day))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, day_no_varid,  AVHRR%dayNo))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, hour_varid,  AVHRR%hours))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, utc_varid,  AVHRR%UTC_msecs))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, time_varid,  AVHRR%time))! start = start_1, count = count_1) )
 
!print*, "gain"
  call check( nf90_put_var(ncid, ch3b_gain_varid, AVHRR%gain_c3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch4_gain_varid, AVHRR%gain_c4))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch5_gain_varid, AVHRR%gain_c5))! start = start_1, count = count_1) )
!--Nspace
!  print*, "nspace"
  call check( nf90_put_var(ncid, nspace_varid, AVHRR%Nspace,start = start_nspace, count = count_nspace) )
!--Ramp
!  print*, "ramp"
  call check( nf90_put_var(ncid, ch3b_ramp_varid, AVHRR%ramp_c3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch4_ramp_varid, AVHRR%ramp_c3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch5_ramp_varid, AVHRR%ramp_c3))! start = start_1, count = count_1) )
!--Mean sp+ict
!  print*, "ramp"
  call check( nf90_put_var(ncid, ch3b_mean_space_varid,AVHRR%sp3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch4_mean_space_varid,AVHRR%sp4))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch5_mean_space_varid, AVHRR%sp5))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch3b_mean_ict_varid,AVHRR%bb3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch4_mean_ict_varid,AVHRR%bb4))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch5_mean_ict_varid,AVHRR%bb5))! start = start_1, count = count_1) )
!--Coef Calibration
 ! print*, "coef calib"
 ! call check( nf90_put_var(ncid, ch3b_coef_calib_varid,AVHRR%calib3, start = start_2, count = count_3) )
 ! call check( nf90_put_var(ncid, ch4_coef_calib_varid,AVHRR%calib4, start = start_2, count = count_3) )
 ! call check( nf90_put_var(ncid, ch5_coef_calib_varid,AVHRR%calib5, start = start_2, count = count_3) )
!--PRT TEMPERATURE
!  print*, "prt temp"
  call check( nf90_put_var(ncid, prt_temp_1_varid,AVHRR%prt1))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, prt_temp_2_varid,AVHRR%prt2))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, prt_temp_3_varid,AVHRR%prt3))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, prt_temp_4_varid,AVHRR%prt4))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, prt_mean_temp_varid,AVHRR%prtmean))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, prt_mean_temp_varid,AVHRR%prtsigma))! start = start_1, count = count_1) )
!--PRT COUNTS
!  print*, "prt counts"
  call check( nf90_put_var(ncid, prt_counts_varid, AVHRR%prtcounts))! start = start_2, count = count_3) )
  call check( nf90_put_var(ncid, prt_number_varid, AVHRR%prtNumber))! start = start_1, count = count_1) )                      
!  print*, "rict"
  call check( nf90_put_var(ncid, ch3b_rict_varid, AVHRR%rict_c3))! start =start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch4_rict_varid,  AVHRR%rict_c4))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, ch5_rict_varid,  AVHRR%rict_c5))! start = start_1, count = count_1) )

!  print*, "counts"
  call check( nf90_put_var(ncid, ch3b_counts_varid, AVHRR%Counts3))! start =start_2, count = count_pixel) )
  call check( nf90_put_var(ncid, ch4_counts_varid,  AVHRR%Counts4))! start = start_2, count = count_pixel) )
  call check( nf90_put_var(ncid, ch5_counts_varid,  AVHRR%Counts5))! start = start_2, count = count_pixel) ) 
!  print*, "sp counts"
  call check( nf90_put_var(ncid, ch3b_space_counts_varid, AVHRR%spaceFilter3))!start = start_2, count = count_10) )
  call check( nf90_put_var(ncid, ch4_space_counts_varid,  AVHRR%spaceFilter4))!start = start_2, count = count_10) )
  call check( nf90_put_var(ncid, ch5_space_counts_varid,  AVHRR%spaceFilter5))!start = start_2, count = count_10) )
!  print*, "ict counts"
  call check( nf90_put_var(ncid, ch3b_ict_counts_varid, AVHRR%bbodyFilter3))! start =start_2, count = count_10) )
  call check( nf90_put_var(ncid, ch4_ict_counts_varid,  AVHRR%bbodyFilter4))! start = start_2, count = count_10) )
  call check( nf90_put_var(ncid, ch5_ict_counts_varid,  AVHRR%bbodyFilter5))! start = start_2, count = count_10) )  
!  print*, "Ne"
  call check( nf90_put_var(ncid, ch3b_ne_varid, AVHRR%nef3, start =start_2))! count = count_pixel) )
  call check( nf90_put_var(ncid, ch4_ne_varid,  AVHRR%nef4, start = start_2))! count = count_pixel) )
  call check( nf90_put_var(ncid, ch5_ne_varid,  AVHRR%nef5, start = start_2))! count = count_pixel) ) 
!  print*, "bt"
  call check( nf90_put_var(ncid, ch3b_bt_varid, AVHRR%array3B, start =start_2))! count = count_pixel) )
  call check( nf90_put_var(ncid, ch4_bt_varid,  AVHRR%array4, start = start_2))! count = count_pixel) )
  call check( nf90_put_var(ncid, ch5_bt_varid,  AVHRR%array5, start = start_2))! count = count_pixel) )

!  print*, "divers temperatures"
  call check( nf90_put_var(ncid, radiator_varid,  AVHRR%Radiator))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, electronics_varid,  AVHRR%electronics))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, cooler_varid,  AVHRR%Cooler))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, baseplate_varid,  AVHRR%baseplate))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, motor_varid,  AVHRR%motor))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, a_d_conv_varid,  AVHRR%a_d_conv))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, patch_varid,  AVHRR%patch))! start = start_1, count = count_1) )
  call check( nf90_put_var(ncid, patch_extended_varid,  AVHRR%patchExtended))!start = start_1, count = count_1) )

! print*, "divers moyennes orbitales"
  call check( nf90_put_var(ncid, mean_prt_orbite_varid, AVHRR%mean_prt_orbite) )
  call check( nf90_put_var(ncid, mean_patch_orbite_varid, AVHRR%mean_patch_orbite) )
  call check( nf90_put_var(ncid, mean_patch_extended_orbite_varid, AVHRR%mean_patchExtended_orbite) )
  call check( nf90_put_var(ncid, ch3b_mean_sp_orbite_varid, AVHRR%mean_space_counts_orbite_c3) )
  call check( nf90_put_var(ncid, ch4_mean_sp_orbite_varid, AVHRR%mean_space_counts_orbite_c4) )
  call check( nf90_put_var(ncid, ch5_mean_sp_orbite_varid, AVHRR%mean_space_counts_orbite_c5) )
  call check( nf90_put_var(ncid, ch3b_mean_bb_orbite_varid, AVHRR%mean_bb_counts_orbite_c3) )
  call check( nf90_put_var(ncid, ch4_mean_bb_orbite_varid, AVHRR%mean_bb_counts_orbite_c4) )
  call check( nf90_put_var(ncid, ch5_mean_bb_orbite_varid, AVHRR%mean_bb_counts_orbite_c5) )
  call check( nf90_put_var(ncid, ch3b_mean_rict_orbite_varid, AVHRR%mean_rict_orbite_c3) )
  call check( nf90_put_var(ncid, ch4_mean_rict_orbite_varid, AVHRR%mean_rict_orbite_c4) )
  call check( nf90_put_var(ncid, ch5_mean_rict_orbite_varid, AVHRR%mean_rict_orbite_c5) )

         
!-Close the netcdf file. 
call check( nf90_close(ncid) )
!write(*,*) "file netcdf close"    
end subroutine write_netcdf

!***************************************************
subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  
!********************************************

end module module_netcdf
