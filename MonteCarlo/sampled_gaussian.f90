!
! Generate gaussian distributed data sampled in nbin equal area bins to 
! conserve the gaussian sigma even with very small numbers of samples.
!
! * Copyright (C) 2019 J.Mittaz University of Reading
! * This code was developed for the EC project Fidelity and Uncertainty in
! * Climate Data Records from Earth Observations (FIDUCEO).
! * Grant Agreement: 638822
! *
! * This program is free software; you can redistribute it and/or modify it
! * under the terms of the GNU General Public License as published by the Free
! * Software Foundation; either version 3 of the License, or (at your option)
! * any later version.
! * This program is distributed in the hope that it will be useful, but WITHOUT
! * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
! * more details.
! *
! * A copy of the GNU General Public License should have been supplied along
! * with this program; if not, see http://www.gnu.org/licenses/
! * ---------------------------------------------------------------------------

MODULE SAMPLED_GAUSSIAN

  USE normal_generator

  IMPLICIT NONE
  
  !
  ! Constants used in NDTRI
  !
  ! sqrt(2pi)
  REAL(kind=8), PARAMETER :: s2pi = 2.50662827463100050242D0

  ! sqrt(2.)
  REAL(kind=8), PARAMETER :: sqrt_2 = SQRT(2.d0)

  ! approximation for 0 <= |y - 0.5| <= 3/8 
  REAL(kind=8), PARAMETER ::  P0(5) = (/&
       -5.99633501014107895267D1,&
       9.80010754185999661536D1,&
       -5.66762857469070293439D1,&
       1.39312609387279679503D1,&
       -1.23916583867381258016D0/)

  REAL(kind=8), PARAMETER :: Q0(8) = (/&
       1.95448858338141759834D0,&
       4.67627912898881538453D0,&
       8.63602421390890590575D1,&
       -2.25462687854119370527D2,&
       2.00260212380060660359D2,&
       -8.20372256168333339912D1,&
       1.59056225126211695515D1,&
       -1.18331621121330003142D0/)

  ! Approximation for interval z = sqrt(-2 log y ) between 2 and 8
  ! i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
  REAL(kind=8), PARAMETER :: P1(9) = (/&
       4.05544892305962419923D0,&
       3.15251094599893866154D1,&
       5.71628192246421288162D1,&
       4.40805073893200834700D1,&
       1.46849561928858024014D1,&
       2.18663306850790267539D0,&
       -1.40256079171354495875D-1,&
       -3.50424626827848203418D-2,&
       -8.57456785154685413611D-4/)

  REAL(kind=8), PARAMETER :: Q1(8) = (/&
       1.57799883256466749731D1,&
       4.53907635128879210584D1,&
       4.13172038254672030440D1,&
       1.50425385692907503408D1,&
       2.50464946208309415979D0,&
       -1.42182922854787788574D-1,&
       -3.80806407691578277194D-2,&
       -9.33259480895457427372D-4/)

  ! Approximation for interval z = sqrt(-2 log y ) between 8 and 64
  ! i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
  REAL(kind=8), PARAMETER :: P2(9) = (/&
       3.23774891776946035970D0,&
       6.91522889068984211695D0,&
       3.93881025292474443415D0,&
       1.33303460815807542389D0,&
       2.01485389549179081538D-1,&
       1.23716634817820021358D-2,&
       3.01581553508235416007D-4,&
       2.65806974686737550832D-6,&
       6.23974539184983293730D-9/)

  REAL(kind=8), PARAMETER :: Q2(8) = (/&
       6.02427039364742014255E0,&
       3.67983563856160859403E0,&
       1.37702099489081330271E0,&
       2.16236993594496635890E-1,&
       1.34204006088543189037E-2,&
       3.28014464682127739104E-4,&
       2.89247864745380683936E-6,&
       6.79019408009981274425E-9/)

  PRIVATE
  PUBLIC :: Get_Gaussian_Eq_Area
  PUBLIC :: Get_Gaussian_Eq_Area_EVEN
  PUBLIC :: Get_Gaussian_Eq_Area_ODD
  PUBLIC :: Get_Gaussian_Eq_Area_EVEN_Switch
  PUBLIC :: Get_Gaussian_Eq_Area_ODD_Switch
  PUBLIC :: Get_Gaussian_Eq_Area_EVEN_Three
  PUBLIC :: Get_Gaussian_Eq_Area_ODD_Three

CONTAINS

  !
  ! Top level wrapper which works out if ODD/EVEN
  ! For speed recommended to work out which routine to use 
  ! before outer loops so don't use this routine in this case
  !
  SUBROUTINE Get_Gaussian_Eq_Area(in_ntrials,array,shuffle,switch)

    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN), OPTIONAL :: shuffle
    LOGICAL, INTENT(IN), OPTIONAL :: switch

    ! Local variables
    LOGICAL :: in_shuffle
    LOGICAL :: in_switch

    IF( PRESENT(shuffle) )THEN
       in_shuffle = shuffle
    ELSE
       in_shuffle = .TRUE.
    ENDIF
    IF( PRESENT(switch) )THEN
       in_switch = switch
    ELSE
       in_switch = .TRUE.
    ENDIF

    IF( in_switch )THEN
       IF( 0 .eq. MOD(in_ntrials,2) )THEN
          CALL Get_Gaussian_Eq_Area_EVEN_Switch(in_ntrials,array,in_shuffle)
       ELSE
          CALL Get_Gaussian_Eq_Area_ODD_Switch(in_ntrials,array,in_shuffle)
       ENDIF
    ELSE
       IF( 0 .eq. MOD(in_ntrials,2) )THEN
          CALL Get_Gaussian_Eq_Area_EVEN(in_ntrials,array,in_shuffle)
       ELSE
          CALL Get_Gaussian_Eq_Area_ODD(in_ntrials,array,in_shuffle)
       ENDIF
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area

  !
  ! Get standard deviation/mean using one pass method
  !
  SUBROUTINE STANDARD_DEV(ndata,array,mean,stdev)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(IN) :: array(ndata)
    REAL, INTENT(OUT) :: mean
    REAL, INTENT(OUT) :: stdev

    ! Local variables
    INTEGER :: I
    INTEGER :: n
    REAL :: M2
    REAL :: delta
    REAL :: delta2
    
    n = 0
    M2 = 0.
    mean = 0.
    DO I=1,ndata
       n=n+1
       delta = array(I) - mean
       mean = mean + delta/n
       delta2 = array(I) - mean
       M2 = M2 + delta*delta2
    END DO

    stdev = SQRT(M2/(n-1))

  END SUBROUTINE STANDARD_DEV

  SUBROUTINE Get_Gaussian_Eq_Area_EVEN_Three(in_ntrials,array,shuffle)

    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT
    LOGICAL, SAVE :: alloc_array=.FALSE.
    INTEGER, SAVE :: narray=-1
    REAL, ALLOCATABLE, SAVE :: three_array(:,:)
    REAL :: dist(3)
    REAL :: mean
    REAL :: stdev
    INTEGER :: POS
    REAL :: MinDist
    
    IF( .not. alloc_array )THEN
       ALLOCATE(three_array(in_ntrials,3),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: three_EVEN: Cannot allocate array'')')
          stop
       ENDIF
       alloc_array = .TRUE.
       narray = in_ntrials
    ELSE IF( in_ntrials .ne. narray )THEN
       DEALLOCATE(three_array)
       ALLOCATE(three_array(in_ntrials,3),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: three_EVEN: Cannot allocate array'')')
          stop
       ENDIF
       narray = in_ntrials
    ENDIF

    DO I=1,3
       CALL Get_Gaussian_Eq_Area_EVEN(in_ntrials,array,.FALSE.)
       three_array(:,I) = array
       CALL Standard_Dev(in_ntrials,array,mean,stdev)
       dist(I) = mean**2+(1-stdev)**2
    END DO

    !
    ! Get minimum distance
    !
    POS = -1
    MinDist = 1e30
    DO I=1,3
       IF( dist(I) .lt. MinDist )THEN
          MinDist = dist(I)
          POS = I
       ENDIF
    END DO
    array = three_array(:,POS)

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_EVEN_Three

  SUBROUTINE Get_Gaussian_Eq_Area_ODD_Three(in_ntrials,array,shuffle)

    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT
    LOGICAL, SAVE :: alloc_array=.FALSE.
    INTEGER, SAVE :: narray=-1
    REAL, ALLOCATABLE, SAVE :: three_array(:,:)
    REAL :: dist(3)
    REAL :: mean
    REAL :: stdev
    INTEGER :: POS
    REAL :: MinDist
    
    IF( .not. alloc_array )THEN
       ALLOCATE(three_array(in_ntrials,3),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: three_EVEN: Cannot allocate array'')')
          stop
       ENDIF
       alloc_array = .TRUE.
       narray = in_ntrials
    ELSE IF( in_ntrials .ne. narray )THEN
       DEALLOCATE(three_array)
       ALLOCATE(three_array(in_ntrials,3),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: three_EVEN: Cannot allocate array'')')
          stop
       ENDIF
       narray = in_ntrials
    ENDIF

    DO I=1,3
       CALL Get_Gaussian_Eq_Area_ODD(in_ntrials,array,.FALSE.)
       three_array(:,I) = array
       CALL Standard_Dev(in_ntrials,array,mean,stdev)
       dist(I) = mean**2+(1-stdev)**2
    END DO

    !
    ! Get minimum distance
    !
    POS = -1
    MinDist = 1e30
    DO I=1,3
       IF( dist(I) .lt. MinDist )THEN
          MinDist = dist(I)
          POS = I
       ENDIF
    END DO
    array = three_array(:,POS)

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_ODD_Three

  SUBROUTINE ShuffleArray(in_ntrials,array)

    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(INOUT) :: array(in_ntrials)

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT
    REAL :: temp
    INTEGER :: itemp
    LOGICAL :: ok

    !
    ! Saved variables to save reallocations
    !
    LOGICAL, SAVE :: random_allocated=.FALSE.
    INTEGER, SAVE :: stored_ntrials=-1
    REAL, ALLOCATABLE, SAVE :: random_array(:)
    INTEGER, ALLOCATABLE, SAVE :: random_index(:)

    !
    ! If we need to allocate shuffle arrays
    !
    IF( .not. random_allocated )THEN
       ALLOCATE(random_array(in_ntrials),random_index(in_ntrials),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: sampled_gaussian: Cannot allocate eqArea arrays'')')
          stop
       ENDIF
       stored_ntrials = in_ntrials
       random_allocated = .TRUE.
    ELSE IF( stored_ntrials .ne. in_ntrials )THEN
       DEALLOCATE(random_array,random_index)
       ALLOCATE(random_array(in_ntrials),random_index(in_ntrials),STAT=STAT)
       IF( 0 .ne. STAT )THEN
          WRITE(*,&
               '(''ERROR: sampled_gaussian: Cannot allocate eqArea arrays'')')
          stop
       ENDIF
       stored_ntrials = in_ntrials
    ENDIF
    !
    ! Now get random array to sort
    !
    DO I=1,in_ntrials
       random_array(I) = uniform_random()
       random_index(I) = I
    END DO
    !
    ! Bubble sort
    !
    ok = .FALSE.
    DO WHILE(.not.ok)
       ok=.TRUE.
       DO I=1,in_ntrials-1
          IF( random_array(I+1) .lt. random_array(I) )THEN
             temp = random_array(I)
             random_array(I) = random_array(I+1)
             random_array(I+1) = temp
             itemp = random_index(I)
             random_index(I) = random_index(I+1)
             random_index(I+1) = itemp
             ok = .FALSE.
          ENDIF
       END DO
    END DO
    !
    ! Now resort output array based on index
    !
    random_array = array
    DO I=1,in_ntrials
       array(I) = random_array(random_index(I))
    END DO
    
  END SUBROUTINE ShuffleArray

  !
  ! Randomly generate nsamples from a Gaussian distribution using
  ! equal area bins to try and maintain the variance
  ! This case when we have a even number of trials
  !
  SUBROUTINE Get_Gaussian_Eq_Area_EVEN(in_ntrials,array,shuffle)
    
    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle

    ! Local variables
    INTEGER :: I
    INTEGER :: ntrials
    REAL(kind=8) :: rno
    REAL(kind=8) :: step_size

    !
    ! Loop round number of 1/2 bins
    ! As this is the even case, no offset start
    !
    ntrials = in_ntrials/2
    step_size = 1./ntrials
    !
    ! Positive case
    !
    DO I=1,ntrials
       rno = (I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = (I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I) = sqrt_2*inv_errfunc(rno)
    END DO
    !
    ! Negative case
    !
    DO I=1,ntrials
       rno = (I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = (I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I+ntrials) = -sqrt_2*inv_errfunc(rno)
    END DO

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_EVEN

  !
  ! Randomly generate nsamples from a Gaussian distribution using
  ! equal area bins to try and maintain the variance
  ! This case when we have a even number of trials
  !
  SUBROUTINE Get_Gaussian_Eq_Area_ODD(in_ntrials,array,shuffle)
    
    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle

    ! Local variables
    INTEGER :: I
    INTEGER :: ntrials
    REAL :: offset
    REAL(kind=8) :: rno
    REAL(kind=8) :: step_size

    !
    ! Central bin needs a different treatment for the ODD case
    ! Do 1/2 step size and randomly allocate +/-
    !
    ntrials = in_ntrials/2
    step_size = 1./(ntrials+0.5)
    offset = 0.5*step_size
    !
    ! central block
    !
    rno = offset*uniform_random_double()
    !
    ! Make sure rno isn't zero or one
    !
    DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
       rno = offset*uniform_random_double()
    END DO
    array(1) = sqrt_2*inv_errfunc(rno)
    IF( uniform_random() .gt. 0.5 )THEN
       array(1) = -array(1)
    ENDIF
    !
    ! Positive case
    !
    DO I=1,ntrials
       rno = offset + (I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = offset + (I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I+1) = sqrt_2*inv_errfunc(rno)
    END DO
    !
    ! Negative case
    !
    DO I=1,ntrials
       rno = offset+(I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = offset+(I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I+ntrials+1) = -sqrt_2*inv_errfunc(rno)
    END DO

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_ODD

  !
  ! Routine to make half of ntrials -ve
  !
  SUBROUTINE Get_Signval_ODD(ndata,signval)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(INOUT), ALLOCATABLE :: signval(:)

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT
    INTEGER :: n_data
    LOGICAL, SAVE :: array_alloc=.FALSE.
    INTEGER, SAVE :: stored_ndata=-1
    REAL, ALLOCATABLE, SAVE :: rand(:)    
    INTEGER, ALLOCATABLE, SAVE :: indexval(:)
    LOGICAL :: ok
    REAL :: temp
    INTEGER :: itemp

    IF( .not. array_alloc )THEN
    
       ALLOCATE(rand(ndata-1),indexval(ndata-1),STAT=STAT)
       IF( STAT .ne. 0 )THEN
          WRITE(*,&
               '(''ERROR: get_signval_odd: Cannot allocate signval arrays'')')
          STOP 1
       ENDIF
       array_alloc = .TRUE.
       stored_ndata = ndata

    ELSE IF( ndata .ne. stored_ndata )THEN

       DEALLOCATE(rand,indexval)
       ALLOCATE(rand(ndata-1),indexval(ndata-1),STAT=STAT)
       IF( STAT .ne. 0 )THEN
          WRITE(*,&
               '(''ERROR: get_signval_odd: Cannot allocate signval arrays'')')
          STOP 1
       ENDIF
       stored_ndata = ndata

    ENDIF
    signval = 0.
    n_data = ndata-1

    !
    ! Get random numbers
    !
    DO I=1,n_data
       rand(I) = uniform_random()
       indexval(I) = I
    END DO
    
    !
    ! Bubble sort
    !
    ok = .FALSE.
    DO WHILE(.not.ok)
       ok=.TRUE.
       DO I=1,n_data-1
          IF( rand(I+1) .lt. rand(I) )THEN
             temp = rand(I)
             rand(I) = rand(I+1)
             rand(I+1) = temp
             itemp = indexval(I)
             indexval(I) = indexval(I+1)
             indexval(I+1) = itemp
             ok = .FALSE.
          ENDIF
       END DO
    END DO

    !
    ! Take top half as negative
    !
    DO I=1,n_data/2
       signval(indexval(I)+1) = -1
    END DO

  END SUBROUTINE Get_Signval_ODD

  SUBROUTINE Get_Signval_EVEN(ndata,signval)

    INTEGER, INTENT(IN) :: ndata
    REAL, INTENT(INOUT), ALLOCATABLE :: signval(:)

    ! Local variables
    INTEGER :: I
    INTEGER :: STAT
    INTEGER :: n_data
    LOGICAL, SAVE :: array_alloc=.FALSE.
    INTEGER, SAVE :: stored_ndata=-1
    REAL, ALLOCATABLE, SAVE :: rand(:)    
    INTEGER, ALLOCATABLE, SAVE :: indexval(:)
    LOGICAL :: ok
    REAL :: temp
    INTEGER :: itemp

    IF( .not. array_alloc )THEN
    
       ALLOCATE(rand(ndata),indexval(ndata),STAT=STAT)
       IF( STAT .ne. 0 )THEN
          WRITE(*,&
               '(''ERROR: get_signval_even: Cannot allocate signval arrays'')')
          STOP 1
       ENDIF
       array_alloc = .TRUE.
       stored_ndata = ndata

    ELSE IF( ndata .ne. stored_ndata )THEN

       DEALLOCATE(rand,indexval)
       ALLOCATE(rand(ndata),indexval(ndata),STAT=STAT)
       IF( STAT .ne. 0 )THEN
          WRITE(*,&
               '(''ERROR: get_signval_even: Cannot allocate signval arrays'')')
          STOP 1
       ENDIF
       stored_ndata = ndata

    ENDIF
    signval = 0.
    n_data = ndata

    !
    ! Get random numbers
    !
    DO I=1,n_data
       rand(I) = uniform_random()
       indexval(I) = I
    END DO
    
    !
    ! Bubble sort
    !
    ok = .FALSE.
    DO WHILE(.not.ok)
       ok=.TRUE.
       DO I=1,n_data-1
          IF( rand(I+1) .lt. rand(I) )THEN
             temp = rand(I)
             rand(I) = rand(I+1)
             rand(I+1) = temp
             itemp = indexval(I)
             indexval(I) = indexval(I+1)
             indexval(I+1) = itemp
             ok = .FALSE.
          ENDIF
       END DO
    END DO

    !
    ! Take top half as negative
    !
    DO I=1,n_data/2
       signval(indexval(I)) = -1
    END DO

  END SUBROUTINE Get_Signval_EVEN

  !
  ! Randomly generate nsamples from a Gaussian distribution using
  ! equal area bins to try and maintain the variance
  ! This case when we have a even number of trials
  ! This is where the number of trials is over 1/2 gaussian and then the
  ! sigmas are randomly assigned +ve and -ve
  ! Has higher sigma resolution but doesn't technically fully cover 
  ! complete Gaussian in the same way as the arrays above
  !
  SUBROUTINE Get_Gaussian_Eq_Area_EVEN_Switch(in_ntrials,array,shuffle,half)
    
    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle
    LOGICAL, INTENT(IN), OPTIONAL :: half

    ! Local variables
    INTEGER :: I
    INTEGER :: ntrials
    REAL(kind=8) :: rno
    REAL(kind=8) :: step_size
    LOGICAL :: force_half
    INTEGER :: STAT
    LOGICAL, SAVE :: alloc_array=.FALSE.
    INTEGER, SAVE :: n_signval = -1
    REAL, ALLOCATABLE, SAVE :: signval(:)

    IF( PRESENT(half) )THEN
       force_half=half
    ELSE
       force_half=.TRUE.
    ENDIF

    !
    ! Loop round number of bins on one side of Gaussian only
    ! As this is the even case, no offset start
    !
    ntrials = in_ntrials
    step_size = 1./ntrials
    !
    ! If force half data to one side
    !
    IF( force_half )THEN
       IF( .not. alloc_array )THEN
          ALLOCATE(signval(ntrials),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             WRITE(*,&
                  '(''ERROR: get_eq_area_even: Cannot allocate signval arrays'')')
             STOP 1
          ENDIF
          alloc_array = .TRUE.
          n_signval = ntrials
       ELSE IF( n_signval .ne. ntrials )THEN
          DEALLOCATE(signval)
          ALLOCATE(signval(ntrials),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             WRITE(*,&
                  '(''ERROR: get_eq_area_even: Cannot allocate signval arrays'')')
             STOP 1
          ENDIF
          n_signval = ntrials
       ENDIF
       CALL Get_Signval_EVEN(ntrials,signval)
    ENDIF
    !
    ! Positive case
    !
    DO I=1,ntrials
       rno = (I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = (I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I) = sqrt_2*inv_errfunc(rno)
       !
       ! Move +/- direction randomly
       !
       IF( force_half )THEN
          IF( signval(I) .lt. 0 )THEN
             array(I) = -array(I)
          ENDIF
       ELSE
          IF( uniform_random() .gt. 0.5 )THEN
             array(I) = -array(I)
          ENDIF
       ENDIF
    END DO

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_EVEN_Switch

  !
  ! Randomly generate nsamples from a Gaussian distribution using
  ! equal area bins to try and maintain the variance
  ! This case when we have a even number of trials
  !
  SUBROUTINE Get_Gaussian_Eq_Area_ODD_Switch(in_ntrials,array,shuffle,half)
    
    INTEGER, INTENT(IN) :: in_ntrials
    REAL, INTENT(OUT) :: array(in_ntrials)
    LOGICAL, INTENT(IN) :: shuffle
    LOGICAL, INTENT(IN), OPTIONAL :: half

    ! Local variables
    INTEGER :: I
    INTEGER :: ntrials
    REAL :: offset
    REAL(kind=8) :: rno
    REAL(kind=8) :: step_size
    INTEGER :: STAT
    LOGICAL :: force_half
    LOGICAL, SAVE :: alloc_array=.FALSE.
    INTEGER, SAVE :: n_signval = -1
    REAL, ALLOCATABLE, SAVE :: signval(:)

    IF( PRESENT(half) )THEN
       force_half=half
    ELSE
       force_half=.TRUE.
    ENDIF

    ! Loop round number of bins on one side of Gaussian only
    !
    ! Central bin needs a different treatment for the ODD case
    ! Do start 1/2 step size and randomly allocate +/-
    !
    ntrials = in_ntrials
    step_size = 1./(ntrials+0.5)
    offset = 0.5*step_size
    !
    ! central block
    !
    rno = offset*uniform_random_double()
    !
    ! Make sure rno isn't zero or one
    !
    DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
       rno = offset*uniform_random_double()
    END DO
    array(1) = sqrt_2*inv_errfunc(rno)
    IF( uniform_random() .gt. 0.5 )THEN
       array(1) = -array(1)
    ENDIF

    IF( force_half )THEN
       IF( .not. alloc_array )THEN
          ALLOCATE(signval(ntrials),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             WRITE(*,&
                  '(''ERROR: get_eq_area_even: Cannot allocate signval arrays'')')
             STOP 1
          ENDIF
          alloc_array = .TRUE.
          n_signval = ntrials
       ELSE IF( n_signval .ne. ntrials )THEN
          DEALLOCATE(signval)
          ALLOCATE(signval(ntrials),STAT=STAT)
          IF( 0 .ne. STAT )THEN
             WRITE(*,&
                  '(''ERROR: get_eq_area_even: Cannot allocate signval arrays'')')
             STOP 1
          ENDIF
          n_signval = ntrials
       ENDIF
       CALL Get_Signval_ODD(ntrials,signval)
    ENDIF
    !
    ! One side only and randomly assign +/-
    !
    DO I=1,ntrials
       rno = offset + (I-1)*step_size+step_size*uniform_random_double()
       !
       ! Make sure rno isn't zero or one
       !
       DO WHILE( rno .lt. 0. .or. rno .ge. 1. )
          rno = offset+(I-1)*step_size+step_size*uniform_random_double()
       END DO
       array(I+1) = sqrt_2*inv_errfunc(rno)
       IF( force_half )THEN
          IF( signval(I+1) .lt. 0 )THEN
             array(I+1) = -array(I+1)
          ENDIF
       ELSE
          IF( uniform_random() .gt. 0.5 )THEN
             array(I+1) = -array(I+1)
          ENDIF
       ENDIF
    END DO

    !
    ! Shuffle within array if needed
    IF( shuffle )THEN
       CALL ShuffleArray(in_ntrials,array)
    ENDIF

  END SUBROUTINE Get_Gaussian_Eq_Area_ODD_Switch

  !
  ! Calculate the inverse error function (double precision)
  !
  ! Code based on C code in the scipy Python library which itself seems to use
  ! the CEPHES library:
  !
  !/*
  ! * Cephes Math Library Release 2.1:  January, 1989
  ! * Copyright 1984, 1987, 1989 by Stephen L. Moshier
  ! * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
  ! */
  !
  REAL(kind=8) FUNCTION inv_errfunc(y)

    REAL(kind=8), INTENT(IN) :: y

    inv_errfunc = NDTRI((y+1.D0)/2.0D0)/sqrt_2
    RETURN

  END FUNCTION inv_errfunc

  !
  ! Polynomial functions
  !
  ! * DESCRIPTION:
  ! *
  ! * Evaluates polynomial of degree N:
  ! *
  ! *                     2          N
  ! * y  =  C  + C x + C x  +...+ C x
  ! *        0    1     2          N
  ! *
  ! * Coefficients are stored in reverse order:
  ! *
  ! * coef[0] = C  , ..., coef[N] = C  .
  ! *            N                   0
  ! *
  ! *  The function p1evl() assumes that coef[N] = 1.0 and is
  ! * omitted from the array.  Its calling arguments are
  ! * otherwise the same as polevl().
  REAL(KIND=8) FUNCTION POLEVL(x,coef,n)

    REAL(kind=8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), INTENT(IN) :: coef(n)

    ! local variables
    REAL(kind=8) :: ans
    INTEGER :: i

    ans = coef(1)
    DO i=2,n
       ans = ans * x + coef(i)
    END DO

    POLEVL = ans
    RETURN

  END FUNCTION POLEVL

  REAL(KIND=8) FUNCTION P1EVL(x,coef,n)

    REAL(kind=8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), INTENT(IN) :: coef(n)

    ! local variables
    REAL(kind=8) :: ans
    INTEGER :: i
    
    ans = x+coef(1)
    DO i=2,n-1
       ans = ans * x + coef(i)
    END DO

    P1EVL = ans
    RETURN

  END FUNCTION P1EVL

  !
  ! Inverse of Normal distribution function
  ! Kept as close to original C as possible
  !
  ! * DESCRIPTION:
  ! *
  ! * Returns the argument, x, for which the area under the
  ! * Gaussian probability density function (integrated from
  ! * minus infinity to x) is equal to y.
  ! *
  ! *
  ! * For small arguments 0 < y < exp(-2), the program computes
  ! * z = sqrt( -2.0 * log(y) );  then the approximation is
  ! * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
  ! * There are two rational functions P/Q, one for 0 < y < exp(-32)
  ! * and the other for y up to exp(-2).  For larger arguments,
  ! * w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
  ! *
  ! *
  ! * ACCURACY:
  ! *
  ! *                      Relative error:
  ! * arithmetic   domain        # trials      peak         rms
  ! *    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
  ! *    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
  ! *
  ! 
  REAL(kind=8) FUNCTION NDTRI(y0)

    REAL(kind=8), INTENT(IN) :: y0

    ! Local variables   
    REAL(kind=8) :: x, y, z, y2, x0, x1
    INTEGER :: code

    !
    ! Check input
    !
    IF( y0 .le. 0. .or. y0 .ge. 1. )THEN
       WRITE(*,'(''NDTRI: ERROR: y0 <= 0 or y0 >= 1.'')')
       STOP 1
    ENDIF

    code = 1
    y = y0
    !
    ! 0.135... = exp(-2)
    !
    IF( y .gt. (1.0d0 - 0.13533528323661269189d0) )THEN
       y = 1.0D0 - y
       code = 0
    ENDIF

    IF( y .gt. 0.13533528323661269189d0 )THEN
       y = y - 0.5d0
       y2 = y*y
       x = y + y * (y2*POLEVL(y2,P0,5)/P1EVL(y2,Q0,9))
       NDTRI = x * s2pi
       RETURN
    ENDIF

    x = SQRT(-2.0*LOG(y))
    x0 = x - LOG(x)/x

    z = 1.0D0/x
    IF( x .lt. 8.0D0 )THEN
       x1 = z * POLEVL(z,P1,9) / P1EVL(z,Q1,9) 
    ELSE
       x1 = z * POLEVL(z,P2,9) / P1EVL(z,Q2,9)
    ENDIF

    x = x0 - x1

    IF( code .ne. 0 )THEN
       x = -x
    ENDIF

    NDTRI = x
    RETURN

  END FUNCTION NDTRI

END MODULE SAMPLED_GAUSSIAN
