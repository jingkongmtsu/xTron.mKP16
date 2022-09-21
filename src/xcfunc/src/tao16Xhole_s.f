*  Deck tao16Xhole_s
       SUBROUTINE tao16Xhole_s(hx,RA,RB,RLAPA,RLAPB,D1RA,D1RB,
     $ TA,TB,NGrid,NDEN)
c
c    ******************************************************************
c    *                                                                *
c    *  reference:                                                    *
c    *  J. Tao, Y. Mo, Accurate semilocal density functional for      *
c    *  condensed-matter physics and quantum  chemistry,              *
c    *  Phys. Rev. Lett. 117, 073001 (2016)                           *
c    * Evaluating the BR89 hx on a grid                               *
c    ******************************************************************
c
      IMPLICIT NONE 
#include "fderiv1.inc"
      REAL*8 hx(NGrid,*),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
     $  RLAPA(NGrid),RLAPB(NGrid),TA(NGrid),TB(NGrid),D1RB(NGrid,3)
      REAL*8 D1,D2,DD,GR,DLAPA,DLAPB,TAUA,TAUB,gra,grb,Tol,tt,ttt,tau0
      REAL*8 pi,ks1,SJ1,SJ3,Lh,Mh,hh1,hh2,hh3
      REAL*8 ff,ss,la,third,sixthrd,eighthrd,tenth
c
      INTEGER NGrid
      INTEGER i
      INTEGER NA,NB
      INTEGER NDEN
      parameter(pi=3.1415926535897930d0)
      parameter(third=1.0d0/3.0d0)
      parameter(sixthrd= 16.0d0/3.0d0)
      parameter(eighthrd = 8.0d0/3.0d0)
      parameter(tenth = 1.0d0/10.0d0)
      Tol = 1.0D-12
      
      ! initilizing the d1vars
c
      DO i = 1,NGrid
        hx(i,1) = 0.0d0
        hx(i,2) = 0.0d0
        SJ1=0.0d0
        SJ3=0.0d0
        ks1= 0.0d0
        ff=0.0d0
        Lh=0.0d0
        Mh=0.0d0
c
       IF((RA(i).gt.Tol).OR.(RB(i).gt.Tol)) THEN
c
         IF (RA(i).gt.Tol) then
            D1 = RA(i)
            gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
            TAUA = TA(i)
            DLAPA = RLAPA(i)
         else
            D1 = 0.0d0
            gra=0.0d0
            TAUA = 0.0d0
            DLAPA = 0.0d0
         endif
c
         IF (RB(i).gt.Tol) then
            D2 = RB(i)
            grb = D1RB(i,1)**2 + D1RB(i,2)**2 + D1RB(i,3)**2
            TAUB = TB(i)
            DLAPB = RLAPB(i)
         else
           D2 = 0.0d0
           grb=0.0d0
           TAUB = 0.0d0
           DLAPB = 0.0d0
         endif
            DD = D1+D2
            GR = 4.0d0*gra
c
         ss = 0.50d0
         ff=(1.0d0+ 0.0010571871787754341d0*(GR**2)/(DD**sixthrd)+ 
     $      0.09432139629937843d0*GR/(DD**eighthrd) )**tenth
          tt = (3.0d0*(pi**2)*DD)**third
          ks1=ff*tt
          ttt = tt**2
          tau0 = 3.0d0*DD*ttt/10.0d0
          la=0.68660d0
          Lh=3.0d0*(la**2-la-0.50d0)*(TAUA-tau0-GR/(72.0d0*DD))
     $    -TAUA+tau0+7.0d0*((2.0d0*la-1.0d0))**2*GR/(18.0d0*DD)
c
          Mh = (2.0d0*la-1.0d0)**2*GR/DD
        IF((ks1.gt.Tol).and.D1.gt.Tol) then
          SJ1 = -(Cos(ks1*ss)/(ks1*ss))+Sin(ks1*ss)/((ks1**2)*ss**2)
          SJ3 = (-15.0d0+(ks1**2)*ss**2)*Cos(ks1*ss)/((ks1**3)*ss**3)-
     $   (-15.0d0 + 6.0d0*(ks1**2)*ss**2)*Sin(ks1*ss)/((ks1**4)*ss**4)
           hh1 = 9.0d0*(SJ1**2)*DD/(ks1*ks1*ss**2)/2.0d0
           hh2= Lh*105.0d0*SJ1*SJ3/(ks1**4*(ss**2))
           hh3= 3675.0d0*Mh*SJ3*SJ3/(8.0d0*(ks1**6)*ss**4)
           hx(i,1) = -(hh1+hh2+hh3)
           hx(i,2) = hx(i,1)
        else
           hx(i,1) = 0.0d0
           hx(i,2) = 0.0d0
        endif  ! ks1
c
         IF(NDEN .ne. 1) then
           write(6,*) " Tao16 hole is coded only for closed shells"
           RETURN
         endif
c
c           write(6,*) "rho_alpha = ", D1 
c           write(6,*) "hx_alpha = ", hx(i,1)
c           write(6,*) "hx_beta = ", hx(i,2)
c Here D1F includes only the derivatives of b2.
c
        ENDIF   !  RA,RB vs Tol
c
      ENDDO  !  NGRID
        RETURN
        END
