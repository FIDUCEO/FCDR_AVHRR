 !SUBROUTINE MAIN0
PROGRAM Get_PRT_Data

  ! Reads in a AVHRR level1b and writes out the PRT data
!***********************
! Version du 9 novembre 2016
! programme qui sápplique à une seule orbite que l'on donne en argument d'entree
! Par conséquent, les stat mensuelles ne sont pas calculees ici.
!************************
! version du 1 nov  2015
! la difference par rapport a la version 6 est que l'on a change le nom des repertoires et des fichiers stat
! pour inclure le nom de línstrument

!************************
! version du 14 oct 2015
! Dans cette version on lit le ramp
!***********************
! 30 oct
! Dans cette version, on imprime egalement les prt counts dont on a besoin pour la methode de Trishchenko
! 6 nov
! On imprime aussi la temperature du detecteur
!***********************
!***********************
! Version du 19 novembre
! CORRECTION D'UN BUG
! LA RADIANCE DE L'ESPACE CHANGE EN FONCTION DU CANAL ET DE L'INSTRUMENET, PRECE!DEMMENT CETTE VALEUR ÉTAIT FIXEE A -2.22
!**********************
!***********************
! Version du 23 novembre
! on imprime aussi la temp moyenne du detecteur sur une orbite
!**********************


  USE GbcsKinds
  USE GbcsTypes
  USE GbcsErrorHandler
  USE GbcsCleanUp
  Use type
  USE NOAA_LoadAVHRRLevel1B
  USE PRT_UTILS
  USE module_netcdf
  USE module_FCDR
  USE module_rescale

  IMPLICIT NONE
   
  CHARACTER(LEN=512)          :: file,f,g
  CHARACTER(LEN=512)          :: output_type
  CHARACTER(LEN=13)           :: date
  CHARACTER(LEN=12)           :: date_2
  CHARACTER(LEN=1000)         :: Filename
  CHARACTER(LEN=512)          :: Command
  CHARACTER(LEN=512)          :: Temp_Filename,temp_file
  CHARACTER(LEN=512)          :: Temp_Dir="/work/scratch/mdesmons/temp/"
  real                        :: t1,t2
  TYPE(AVHRR_Data)            :: AVHRR
  TYPE(Imagery)               :: IMG
  TYPE(Satellite)             :: SAT

  INTEGER                     :: Ios
  INTEGER                     :: I,J,k,sizer
  character(len=256)          :: mkdirCmd1,mkdirCmd2
    logical                   :: dirExists1,dirExists2
  REAL(GbcsReal), ALLOCATABLE :: min_gradient(:)
  REAL(GbcsReal), ALLOCATABLE :: max_gradient(:)
  REAL(GbcsReal), ALLOCATABLE :: mean_gradient(:)
  REAL(GbcsReal)              :: gradients(6)
  REAL(GbcsReal)              :: mean_temp, mean_patch
  REAL(GbcsReal)              :: mean_grad
  REAL(GbcsReal)              :: min_temp, max_temp
  REAL(GbcsReal), ALLOCATABLE :: temperature(:)
  REAL(GbcsReal)              :: max_delta_cal(3)
  INTEGER                     :: temp1, temp2
  INTEGER                     :: POS
  LOGICAL                     :: remove_data
  CHARACTER(LEN=256)          :: outFile
  CHARACTER(LEN=256)          :: path="/group_workspaces/cems2/fiduceo/Users/mdesmons/avhrr_l1b/",&
                                 dir_filtre="essai_filtre/", &
!                                           1234567890123456789
                                 dir_stats="stats_xx/xxxx/xx/", &
!                                            1234567890123456789012
                                 dir_values="values_xx/xxxx/xx/xx/",&
                                 dir_lists_day_hours="listes_day_hours/"

  CHARACTER(LEN=356)          :: filelistedayhours="liste_day_hours_xx_xx_xxxx.data", & 
!                                                  1234567890123456789012345678901234        
                                 file_coef_calib3="coef_calib_c3_xx_xxxxxxxxxxxx.data", &
                                 file_coef_calib4="coef_calib_c4_xx_xxxxxxxxxxxx.data",&
                                 file_coef_calib5="coef_calib_c5_xx_xxxxxxxxxxxx.data", & 
!                                            1234567890123456789012345678     
                                 filestatc3="stat_c3_xx_xxxxxxxxxxxx.data", &
                                 filestatc4="stat_c4_xx_xxxxxxxxxxxx.data", &
                                 filestatc5="stat_c5_xx_xxxxxxxxxxxx.data", &
!                                             12345678901234567890123456789
                                 filestatprt="stat_prt_xx_xxxxxxxxxxxx.data",& 
!                                                  12345678901234567890123456789
                                 filestatdetector="stat_det_xx_xxxxxxxxxxxx.data" 

  CHARACTER(LEN=2)            :: ins
  CHARACTER(LEN=10)           :: instrument,m,y,d
  INTEGER                     :: year, month,day
  CHARACTER(len=2)            :: dday
  CHARACTER(len=6)            ::hms
                                            
  INTEGER                     :: ndata37_day
  INTEGER                     :: ndata37_night
  INTEGER                     :: ndata11
  INTEGER                     :: ndata12
  INTEGER                     :: nan37
  INTEGER                     :: nan11
  INTEGER                     :: nan12
  REAL                        :: min37, max37
  REAL                        :: min11, max11
  REAL                        :: min12, max12
  REAL                        :: nedt_c3, gain_c3, &
                                 nmean_sp_c3, somme_sp_c3, stdev_sp_c3, &
                                 nmean_bb_c3, somme_bb_c3, stdev_bb_c3
  REAL                        :: nedt_c4, gain_c4, &
                                 nmean_sp_c4, somme_sp_c4, stdev_sp_c4, & 
                                 nmean_bb_c4, somme_bb_c4, stdev_bb_c4
  REAL                        :: nedt_c5, gain_c5, &
                                 nmean_sp_c5, somme_sp_c5, stdev_sp_c5, &
                                 nmean_bb_c5, somme_bb_c5, stdev_bb_c5
  REAL                        :: somme_prt, nmean_prt, stdev_prt
  REAL                        :: somme_hours, nhours, mean_hours
  INTEGER                     :: STAT
!-------------------------------------------1234567890123456789012345678901234
  character (len = 34)           :: file_nc="yyyymmddhhmm00-ESACCI-AVHHRxx_G.nc"
  integer                        :: ncid

  call CPU_TIME(t1)
  
!-On recupere les arguments de la ligne de commande
  call getarg(1,instrument)
  call getarg(2, y)
  read(y,"(i4)")year
  call getarg(3, m)
  read(m,"(i2)")month
  CALL getarg(4, d)
  READ(d,"(i2.2)")day
  CALL getarg(5, f)
  READ(f,"(a131)",IOSTAT=Ios)output_type
  CALL getarg(6, g)
  READ(g,"(a131)",IOSTAT=Ios)Filename
  

!-On determine l'abreviation du nom de l'instrument
  if (instrument .eq."AVHRRTN_G") then 
    ins="TN"
  end if 
  if (instrument .eq. "AVHRR06_G") then
    ins="NA"
  end if 
  if (instrument .eq. "AVHRR07_G") then 
    ins="NC"
   end if 
  if (instrument .eq. "AVHRR08_G") then 
    ins="NE"
  end if 
  if (instrument .eq. "AVHRR09_G") then 
    ins="NF"
  end if 
  if (instrument .eq. "AVHRR10_G") then 
    ins="NG"
  end if 
  if (instrument .eq. "AVHRR11_G") then 
    ins="NH"
  end if 
  if (instrument .eq. "AVHRR12_G") then 
    ins="ND"
  end if 
  if (instrument .eq. "AVHRR14_G") then 
    ins="NJ"
  end if 
  if (instrument .eq. "AVHRR15_G") then 
    ins="NK"
  end if 
  if (instrument .eq. "AVHRR16_G") then 
    ins="NL"
  end if 
  if (instrument .eq. "AVHRR17_G") then 
    ins="NM"
  end if 
  if (instrument .eq. "AVHRR18_G") then 
    ins="NN"
  end if 
  if (instrument .eq. "AVHRR19_G") then 
    ins="NP"
  end if 
  if (instrument .eq. "AVHRRMTB_G") then 
    ins="M1"
  end if 
  if (instrument .eq. "AVHRRMTA_G") then 
    ins="M2"
  end if 
 
if (output_type .eq."easy_fcdr") then   
    dir_filtre="FCDR/"
end if
 
!-On completer le nom des repertoires et des fichiers
  write(dir_values(8:9),"(a2)")ins
  write(dir_values(11:14),"(i4)")year
  write(dir_values(16:17),"(i2.2)")month
  write(dir_values(19:20),"(i2.2)")day

  write(dir_stats(7:8),"(a2)")ins
  write(dir_stats(10:13),"(i4)")year
  write(dir_stats(15:16),"(i2.2)")month

  inquire(file=trim(path)//trim(dir_filtre)//trim(dir_values)//"/.", exist=dirExists1 ) 
  if (dirExists1) then
    write (*,*) "Directory values already exists: "
  else
    mkdirCmd1 = 'mkdir -p '//trim(path)//trim(dir_filtre)//trim(dir_values)
    write(*,'(a)') "Creating dir values "
    call system( mkdirCmd1 )
  endif

  inquire(file=trim(path)//trim(dir_filtre)//trim(dir_stats)//'/.', exist=dirExists2 ) 
  if (dirExists2) then
      write (*,*) "Directory already stats exists: "
  else
    mkdirCmd2 = 'mkdir -p '//trim(path)//trim(dir_filtre)//trim(dir_stats)
    write(*,'(a)') "Creating dir stats "
    call system( mkdirCmd2 )
  endif
 
  if ( (ins .eq. "M2") .or. (ins .eq. "M1")) then
    date=Filename(89:101)
    date_2=filename(89:100)
  else
    date=Filename(88:100)
    date_2=Filename(88:99)
  end if
!-On definit le nom du fichier netCDF  
        dday=filename(73:74)
        hms=filename(96:99)
        write(file_nc(1:4), "(i4.4)") year
        write(file_nc(5:6), "(i2.2)") month
        write(file_nc(7:8), "(a2)") dday 
        write(file_nc(9:12), "(a4)") hms 
        write(file_nc(23:31),"(a9)")instrument
        print*, "file_nc ", file_nc
     
  write(file_coef_calib3(15:16), "(a2)") ins
  write(file_coef_calib4(15:16), "(a2)") ins
  write(file_coef_calib5(15:16), "(a2)") ins
  write(file_coef_calib3(18:29), "(a12)") date
  write(file_coef_calib4(18:29), "(a12)") date
  write(file_coef_calib5(18:29), "(a12)") date
   
  write(filestatc4(9:10),"(a2)")ins
  write(filestatc5(9:10),"(a2)")ins
  write(filestatc3(9:10),"(a2)")ins
  write(filestatprt(10:11),"(a2)")ins
  write(filestatdetector(10:11),"(a2)")ins
  
  write(filestatc4(12:23),"(a12)") date
  write(filestatc5(12:23),"(a12)") date
  write(filestatc3(12:23),"(a12)") date
  write(filestatprt(13:24),"(a12)") date
  write(filestatdetector(13:24),"(a12)") date
  
  open(Unit=100, file=trim(path)//trim(dir_lists_day_hours)//trim(filelistedayhours), status="replace")
  open(Unit=400, file=trim(path)//trim(dir_filtre)//trim(dir_stats)//trim(filestatc4), status="replace")
  open(Unit=500, file=trim(path)//trim(dir_filtre)//trim(dir_stats)//trim(filestatc5), status="replace")
  open(Unit=300, file=trim(path)//trim(dir_filtre)//trim(dir_stats)//trim(filestatc3), status="replace")
  open(Unit=600, file=trim(path)//trim(dir_filtre)//trim(dir_stats)//trim(filestatprt), status="replace")
  open(Unit=700, file=trim(path)//trim(dir_filtre)//trim(dir_stats)//trim(filestatdetector), status="replace") 

!--Get unique filename based on process ID
  WRITE(Temp_File,'(''avhrr_data.'',i8.8)')GETPID()
  temp_filename=trim(temp_dir)//trim(temp_file)
  print*, "temp filename ", temp_filename

  remove_data = .TRUE.
  Ios = 0
  POS = INDEX(Filename,'.gz')
  print*, "POS ",POS
  IF( 0 .ne. POS )THEN
    write(command,'(''cp -f '',a,'' '',a,''.gz'')')&
    TRIM(Filename),TRIM(Temp_Filename)
    !print*, "commande ",command
    call SYSTEM(command)
    write(command,'(''gunzip -f '',a,''.gz'')')TRIM(Temp_Filename)
    !print*, "command",command
    call SYSTEM(command)
    IMG%DataFile = TRIM(Temp_Filename)
    remove_data = .TRUE.
  ELSE
!   Get filename
    POS = -1
    Find_loop: DO I=LEN(Filename),1,-1
      IF( '/' .eq. Filename(I:I) )THEN
        POS = I+1
        EXIT Find_Loop
      ENDIF
    END DO Find_loop
    IF( -1 .eq. POS )THEN
      Temp_Filename = TRIM(Filename)
      print*, "temp filename ", Temp_Filename
      IMG%DataDir = TRIM(Temp_Dir)
      IMG%DataFile = TRIM(Temp_Filename)
    ELSE
      Temp_Filename = Filename(POS:LEN(Filename))
      print*, "temp filename ", Temp_Filename
      Temp_Dir = Filename(1:POS-1)
!     IMG%DataDir = TRIM(Temp_Dir)
      IMG%DataFile = TRIM(Temp_Filename)
    ENDIF
  ENDIF !POS ne 0
   
  IMG%GbcsDataPath = '/net/orbit033l/home/jmittaz/GBCS_Oper/branch_phys/dat/'
  IMG%LandSeaMaskFile = ' '
  IMG%N_Chans = 3
  IMG%Gbcs_Map(1) = Gbcs_IR_03x
  IMG%Gbcs_Map(2) = Gbcs_IR_10x
  IMG%Gbcs_Map(3) = Gbcs_IR_12x
     
!-Read in Level 1B data
  CALL Load_Imagery( IMG, outputData=AVHRR, &
                      use_new_calibration=.TRUE., use_walton=.FALSE.,no_landmask=.TRUE.)
   CALL getarg(5, f)
  READ(f,"(a131)",IOSTAT=Ios)output_type
  CALL getarg(6, g)
  READ(g,"(a131)",IOSTAT=Ios)Filename
  
  if (output_type .eq."data") then 
    print*, "writing data"
    CALL write_netcdf(trim(path)//trim(dir_filtre)//trim(dir_values)//trim(file_nc), AVHRR) 
  else if (output_type .eq."easy_fcdr") then   
    print*, "Creation FCDR"
  !CALL rescale(AVHRR,AVHRR%arraysize)
  CALL write_FCDR(trim(path)//trim(dir_filtre)//trim(dir_values)//trim(file_nc), AVHRR) 
  end if 
  

!-Remove file
 write(command,'(''rm -rf '',a)')TRIM(Temp_Filename)

 IF ( AVHRR%start_valid .lt. AVHRR%stop_valid ) then
   !call average_per_orbite(AVHRR)        
    !gain_c3=(AVHRR%mean_rict_orbite_c3+AVHRR%nspace(1))/(AVHRR%mean_bb_counts_orbite_c3-AVHRR%mean_space_counts_orbite_c3)
    !nedt_c3=-1*stdev_sp_c3*gain_c3
    !gain_c4=(AVHRR%mean_rict_orbite_c4+AVHRR%nspace(2))/(AVHRR%mean_bb_counts_orbite_c4-AVHRR%mean_space_counts_orbite_c4)
    !nedt_c4=-1*stdev_sp_c4*gain_c4
    !gain_c5=(AVHRR%mean_rict_orbite_c5+AVHRR%nspace(3))/(AVHRR%mean_bb_counts_orbite_c5-AVHRR%mean_space_counts_orbite_c5)
    !nedt_c5=-1*stdev_sp_c5*gain_c5
   
    write(300,*),Filename(85:86),"     ",Filename(88:105),&
    AVHRR%year(1),AVHRR%month(1),AVHRR%day(1), AVHRR%dayNo(1), &
   AVHRR%hours(1),AVHRR%utc_msecs(1),AVHRR%time(1), &
    AVHRR%mean_rict_orbite_c3, AVHRR%mean_bb_counts_orbite_c3, stdev_bb_c3, &
    AVHRR%mean_space_counts_orbite_c3,stdev_sp_c3,  &
    gain_c3, nedt_c3, AVHRR%mean_calib3_orbite(2), &
    -1*stdev_sp_c3*AVHRR%mean_calib3_orbite(2)

    write(400,*)Filename(85:86),"     ",Filename(88:105),AVHRR%year(1),AVHRR%month(1),AVHRR%day(1), AVHRR%dayNo(1), &
    AVHRR%hours(1),AVHRR%utc_msecs(1),AVHRR%time(1), &
    AVHRR%mean_rict_orbite_c4, AVHRR%mean_bb_counts_orbite_c4, stdev_bb_c4, &
    AVHRR%mean_space_counts_orbite_c4,stdev_sp_c4, &
    gain_c4, nedt_c4, AVHRR%mean_calib4_orbite(2), &
    -1*stdev_sp_c4*AVHRR%mean_calib4_orbite(2)

    write(500,*)Filename(85:86),"     ",Filename(88:105),AVHRR%year(1),AVHRR%month(1),AVHRR%day(1), AVHRR%dayNo(1), &
    AVHRR%hours(1),AVHRR%utc_msecs(1),AVHRR%time(1), &
    AVHRR%mean_rict_orbite_c5, AVHRR%mean_bb_counts_orbite_c5, stdev_bb_c5, &
    AVHRR%mean_space_counts_orbite_c5,stdev_sp_c5,  &
    gain_c5, nedt_c5,AVHRR%mean_calib5_orbite(2), &
    -1*stdev_sp_c5*AVHRR%mean_calib5_orbite(2)

    write(600,*)Filename(85:86),"     ",Filename(88:105), &
    AVHRR%year(1),AVHRR%month(1),AVHRR%day(1), AVHRR%dayNo(1), &
    AVHRR%hours(1),AVHRR%utc_msecs(1),AVHRR%time(1), &
    AVHRR%mean_prt_orbite, stdev_prt 

   write(700,*)Filename(85:86),"     ",Filename(88:105), &
    AVHRR%year(1),AVHRR%month(1),AVHRR%day(1), AVHRR%dayNo(1), &
    AVHRR%hours(1),AVHRR%utc_msecs(1),AVHRR%time(1), &
    AVHRR%mean_patch_orbite, AVHRR%mean_patchExtended_orbite
  end if

  Close(100)
  close(300)
  close(400)
  close(500)
  close(600)
  close(700)
  
  if( remove_data )then
     print*, "enlever fichier"
     !print*, command
     call SYSTEM(command)
  endif
  
  106 Continue

  STAT = Release_Imagery(IMG)
  call DeAllocate_OutData( AVHRR )
  CALL CPU_TIME(t2)
  print*, "duree d'execution",  t2-t1, (t2-t1)/3600.

END PROGRAM Get_PRT_Data
!END SUBROUTINE MAIN0

      
