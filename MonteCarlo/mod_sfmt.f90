module MOD_SFMT
  implicit none
  public
!/*-----------------
!  BASIC DEFINITIONS
!  -----------------*/
!/** Mersenne Exponent. The period of the sequence 
! *  is a multiple of 2^MEXP-1.
! * #define MEXP 19937 */
!/** SFMT generator has an internal state array of 128-bit integers,
! * and N is its size. */
    integer,parameter:: MEXP = 19937
    integer,parameter:: N = (MEXP / 128 + 1)
!/** N32 is the size of internal state array when regarded as an array
! * of 32-bit integers.*/
    integer,parameter:: N32 = (N * 4)
!/** N64 is the size of internal state array when regarded as an array
! * of 64-bit integers.*/
    integer,parameter:: N64 = (N * 2)
    


    
    integer(kind=4):: psfmt32(0:N32-1)
    integer(kind=8):: psfmt64(0:N64-1)
    integer:: idx
    integer:: initialized = 0


!====for MEXP == 19937
    integer,parameter:: POS1 = 122
    integer,parameter:: SL1 = 18
    integer,parameter:: SL2 = 1
    integer,parameter:: SR1 = 11
    integer,parameter:: SR2 = 1
!   integer,parameter:: MSK1 =      ! = 0x001fffefU = 0xdfffffefU
!   integer,parameter:: MSK2 =      ! = 0x001ecb7fU = 0xddfecb7fU
!   integer,parameter:: MSK3 =      ! = 0x001affffU = 0xbffaffffU
!   integer,parameter:: MSK4 =      ! = 0x001ffff6U = 0xbffffff6U
    integer(kind=8),parameter:: MSK12 = 8667995624701935_8    ! = 0x001ecb7f001fffefU
    integer(kind=8),parameter:: MSK34 = 9007156306837503_8    ! = 0x001ffff6001affffU
    integer(kind=8),parameter:: MSK_L = 70364449226751_8  ! = 0x00003fff00003fffU
    integer(kind=4),parameter:: PARITY1 = 1
    integer(kind=4),parameter:: PARITY2 = 0
    integer(kind=4),parameter:: PARITY3 = 0
    integer(kind=4),parameter:: PARITY4 = 331998852
!#define ALTI_SL2_PERM \
!(vector unsigned char)(1,2,3,23,5,6,7,0,9,10,11,4,13,14,15,8)
!#define ALTI_SL2_PERM64 \
!(vector unsigned char)(1,2,3,4,5,6,7,31,9,10,11,12,13,14,15,0)
!#define ALTI_SR2_PERM \
!(vector unsigned char)(7,0,1,2,11,4,5,6,15,8,9,10,17,12,13,14)
!#define ALTI_SR2_PERM64 \
!(vector unsigned char)(15,0,1,2,3,4,5,6,17,8,9,10,11,12,13,14)
!#define IDSTR  "SFMT-19937:122-18-1-11-1:dfffffef-ddfecb7f-bffaffff-bffffff6"
    
    integer(kind=4):: parity(0:3)= (/PARITY1, PARITY2, PARITY3, PARITY4/)

end module
