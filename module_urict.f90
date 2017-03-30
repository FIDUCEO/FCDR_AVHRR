MODULE module_urict
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

   SUBROUTINE calcul_urict(i,outData,coefs1,coefs2,coefs3)

INTEGER, INTENT(IN)                  :: i
TYPE(AVHRR_Data), INTENT (INOUT)     :: outData
REAL, DIMENSION(7), INTENT(IN)       :: coefs1,coefs2,coefs3
  
REAL                                 :: prt_accuracy=0.1, prt_noise=0., &
                                        urict, tstar, &
                                        aval,bval,nuc, &
                                        uc1=0,uc2=0,unuc=0,uaval=0,ubval=0

!Coefs : We use pre-launch values, for the moment the uncertainty =0
!
!             3  =  nuc
!             4  =  aVal (band coefficient)
!             5  =  bVal (band coefficient)


!--On calcule u(Tict)
! Tict=1/4*Somme (Tprti)
! Tprti=a+b*Ti
! with a: systematic component 0.1 K
!      b: random component 0.015K
!Following the GUM: (u(tict))**2= u(a)**2+Tict**2 * u(b)**2 + u(b)**2 +sigmaTict**2 / 4

if ((outData%prtmean(i) .ne. NAN_R )&
.and. (outData%prt1(i) .ne. NAN_R ) &
.and. (outData%prt2(i) .ne. NAN_R )&
.and. (outData%prt3(i) .ne. NAN_R )&
.and. (outData%prt4(i) .ne. NAN_R )) then
    outData%prtsigma(i)=((outData%prt1(i)-outData%prtmean(i))**2&
                +(outData%prt2(i)-outData%prtmean(i))**2&
                +(outData%prt3(i)-outData%prtmean(i))**2&
                +(outData%prt4(i)-outData%prtmean(i))**2)/3.
  
    outData%utict(i)=sqrt(prt_accuracy**2+prt_noise**2*outData%prtmean(i)**2+prt_noise**2*outData%prtsigma(i)**2)
 
end if


!Ch3 
nuc=coefs1(3)
aval=coefs1(4)
bval=coefs1(5)
unuc=0

!--On calcule u(Rict)
! T*=(Tict-A)/B
! Following the GUM: u(T*)**2= u(A)**2*(1/B)**2 +u(b)**2*(-(Tict-A)/B**2)**2 +(1/B)**2* u(Tict)**2
! u(T*)**2=(1/B)**2* u(Tict)**2
if (outData%prtmean(i) .ne.  -1.00000002E+30 ) then 
  tstar=(outData%prtmean(i)-aval)/bval

  outdata%Rict_c3(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)
 
  outdata%utstar3(i)=sqrt(uaval**2*outdata%dtstar_over_daval3(i)**2 &
                      +ubval**2*outdata%dtstar_over_dbval3(i)**2 &
                      +outdata%utict(i)**2*outdata%dtstar_over_dnuc3(i)**2)
  !print*, outdata%utstar3(i)  
  outData%urict3(i)=sqrt(unuc**2*outdata%drict_over_dnuc3(i)**2 &
                     +outdata%utstar3(i)**2*outdata%drict_over_dtstar3(i)**2)
!print*, outdata%rict_c3(i), outdata%utict(i),outdata%urict3(i)
  !print*, "urict", outdata%urict3(i)
end if
!Ch4
nuc=coefs2(3)
aval=coefs2(4)
bval=coefs2(5)
unuc=0
uaval=0
ubval=0
!--On calcule u(Rict)
if (outData%prtmean(i) .ne. NAN_R )then
tstar=(outData%prtmean(i)-aval)/bval

outdata%Rict_c4(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)

outdata%urict4(i)=sqrt(c1*nuc**3*(c2*nuc)**2/tstar**4*1/(exp(c2*nuc/tstar)-1))
 
outdata%utstar4(i)=sqrt(uaval**2*outdata%dtstar_over_daval4(i)**2 &
                      +ubval**2*outdata%dtstar_over_dbval4(i)**2 &
                      +outdata%utict(i)**2*outdata%dtstar_over_dnuc4(i)**2)
  
outData%urict4(i)=sqrt(unuc**2*outdata%drict_over_dnuc4(i)**2 &
                     +outdata%utstar4(i)**2*outdata%drict_over_dtstar4(i)**2)
!print*, outdata%rict_c4(i), outdata%utict(i),outdata%urict4(i)
end if
!Ch5 
nuc=coefs3(3)
aval=coefs3(4)
bval=coefs3(5)
unuc=0
if (outData%prtmean(i) .ne. NAN_R )then
tstar=(outData%prtmean(i)-aval)/bval
outdata%Rict_c5(i)=c1*nuc**3/(exp(c2*nuc/tstar)-1.)
!print*, tstar, outdata%rict_c5(i)

outdata%urict5(i)=sqrt(c1*nuc**3*(c2*nuc)**2/tstar**4*1/(exp(c2*nuc/tstar)-1))
 
outdata%utstar5(i)=sqrt(uaval**2*outdata%dtstar_over_daval5(i)**2 &
                      +ubval**2*outdata%dtstar_over_dbval5(i)**2 &
                      +outdata%utict(i)**2*outdata%dtstar_over_dnuc5(i)**2)
  
outData%urict5(i)=sqrt(unuc**2*outdata%drict_over_dnuc5(i)**2 &
                     +outdata%utstar5(i)**2*outdata%drict_over_dtstar5(i)**2)

end if


END SUBROUTINE calcul_urict
 

END MODULE module_urict

