*  Deck tao17Xhole_s
       SUBROUTINE tao17Xhole_s(hx,RA,RB,RLAPA,RLAPB,D1RA,D1RB,
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
      REAL*8 D1,D2,DT,GR,DLAPA,DLAPB,DLAPD,TAUA,TAUB,gra,grb,Tol
      REAL*8 pi,Zw,kf,Zueg,tau0,lamb,tauw
      REAL*8 LT,LP,LL,HH,HHs1,Hpbe,Hz1,eff,bet,Cc,Gg,Kk
      REAL*8 third,sixthrd,eighthrd,tenth,hh1,hh2,hh3,dd0,dd1,dd2,dd3
      REAL*8 p1,p2,p3,p4,p5,p6,Aa,Bb,Dd,Sr,ss,sf
      REAL*8 hx1,hx2,hx3,ttt,svntwo,nine2,hh0,Fxx,eunifx
      DOUBLE PRECISION Ftps,D11F(N_FUNC_DERIV_1)
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
      parameter(hh0=0.0060d0)
      parameter(hh1=2.89160d0)
      parameter(hh2=0.77680d0)
      parameter(hh3=2.08760d0)
      parameter(dd0=13.6950d0)
      parameter(dd1= -0.22190d0)
      parameter(dd2= 4.99170d0)
      parameter(dd3= 0.79720d0)
      parameter(p1=0.03020d0)
      parameter(p2= -0.10350d0)
      parameter(p3= 0.12720d0)
      parameter(p4= 0.12030d0)
      parameter(p5= 0.48590d0)
      parameter(p6= 0.10080d0)
      parameter(Aa=0.7572110d0)
      parameter(Bb= -0.1063640d0)
      parameter(Dd= 0.609650d0)
      parameter(svntwo = 7.0d0/2.0d0)
      parameter(nine2 = 9.0d0/2.0d0)
      Tol = 1.0D-12
      
      ! initilizing the d1vars
c
        DO i = 1,NGrid
          hx(i,1) = 0.0d0
          hx(i,2) = 0.0d0
          kf= 0.0d0
          Sr = 0.0d0
          Zw = 0.0d0
          Zueg = 0.0d0
        IF((RA(i).gt.Tol).OR.(RB(i).gt.Tol)) THEN
c
        IF (RA(i).gt.Tol) then
           D1 = RA(i)
           gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
           TAUA = TA(i)
           DLAPA = RLAPA(i)
         else
           D1 = 0.0d0
           gra =0.0d0
           TAUA = 0.0d0
           DLAPA = 0.0d0
         endif
c
           D2=D1
           grb=gra
           TAUB=TAUA
           DLAPB=DLAPA
           DT = D1+D2
c  This routine accepts only closed shell cases
           IF(NDEN .ne. 1) then
           write(6,*) " Tao17 hole is coded only for closed shells"
           RETURN
         endif
           GR = 4.0d0*gra 
           DLAPD=DLAPA+DLAPB
           kf = (3.0d0*DT*pi**2)**third
           Sr= 2.0d0*DSQRT(gra)/(2.0d0*kf*DT)
           ss = 0.50d0
c         sf here is s in the mathematica file and uf in the Tao paper
c          ss is rr in the math graphs
           sf= kf*ss
           tauw = GR/(8.0d0*DT)
           Zw= tauw/TAUA
           ttt = kf**2
           tau0 = 3.0d0*DT*ttt/10.0d0
           Zueg = TAUA/tau0
           eff = Erfc((Sr**2-36.0d0)/6.0d0)
           LT = -(3.0d0*Zueg/10.0d0-0.90d0+5.*(Sr**2)/6.0d0)/3.0d0
           LP = 0.20d0 - 2.0d0*(Sr**2)/27.0d0
           LL = 0.5*eff*LT + (1.0d0-eff)*LP
c
         Hz1 = (hh0 + hh1*Sr**2 + hh2*Sr**4 + hh3*Sr**6)/
     $  (dd0 + dd1*Sr**2 + dd2*Sr**4 + dd3*Sr**6)
c
         Hpbe = (p1*Sr**2 + p2*Sr**4 + p3*Sr**6)/
     $   (1.0d0 + p4*Sr**2 + p5*Sr**4 + p6*Sr**6)
c 
         HHs1=Hpbe*(1.0d0-eff/2.0d0) + Hz1*eff/2.0d0
         HH = HHs1*Zw**3
c
         bet = Aa + HH
         lamb = Dd + HH
c
c        Cc = (4.0d0*LL + 3.0d0*Aa**3 + 9.0d0*Aa*Aa*HH - 9.0d0*Aa*Dd**2 
c    $      - 18.0d0*Aa*Dd*HH + 8.0d0*Bb*lamb)/8.0d0
         Cc= 0.1250d0*(-1.2304319884513730d0
     $     -3.1490898640110010d0*HH + 4.0d0*LL - 0.8509120d0*lamb)
c
c   need to find and define Fxx before that:
c
       CALL functional_tpssx_pw_close
     &(TOL,D1,D1RA(i,1),D1RA(i,2),D1RA(i,3),TAUA,Ftps,D11F)
       eunifx = -3.0d0*(3.0d0*D1/(4.0d0*pi))**third/2.0d0
       IF((D1*eunifx).lt.-Tol) then
         Fxx = Ftps/(D1*eunifx)
       else
         Fxx = 0.0d0
       endif
c
       IF(HH.gt.Tol) then
       Gg=-63.0d0*lamb**3*(Fxx+Aa*dlog(bet/lamb)+HH*dlog(bet/HH))/8.0d0
     $    -24.0d0*(lamb**svntwo)*(3.0d0*Aa/(DSQRT(HH)
     $    + DSQRT(bet))-DSQRT(pi))/5.0d0 +603.0d0*Aa*(lamb**3)/40.0d0
     $    - 19.0d0*Bb*lamb**2/10.0d0-11.0d0*Cc*lamb/10.0d0
       else
        Gg = 0.0d0
       endif
c
       IF((HH.gt.Tol).and.bet.gt.Tol) then
        Kk=(2.0d0*lamb*(-15.0d0*Gg-6.0d0*Cc*lamb-4.0d0*Bb*(lamb**2) 
     $    - 18.0d0*Aa*(lamb**3)
     $    + 12.0d0*(lamb**svntwo)*((3.0d0*Aa)/(DSQRT(bet) 
     $    + DSQRT(HH))- DSQRT(pi))))/105.0d0
       else
        Kk = 0.0d0
       endif
c
       IF(sf.gt.Tol) then
        hx1= (-9.0d0*(1.0d0 - DEXP(-(Aa*(sf**2)))))/
     $        (4.0d0*DEXP(HH*(sf**2))*sf**4)
       hx2= DEXP(-(Dd*sf**2)-HH*sf**2)*(Bb+(9.0d0*Aa)/(4.0d0*sf**2))
       hx3= DEXP(-(Dd*sf**2)-HH*sf**2)*(Cc*sf**2 +Gg*sf**4 +Kk*sf**6)
       hx(i,1)= hx1+hx2+hx3
       else
        hx(i,1) = 0.0d0
       endif
        hx(i,2) = hx(i,1)
      
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
