*  Deck hirao00Xhole_s 
       SUBROUTINE hirao00Xhole_s(hx,ss,RA,RB,D1RA,D1RB,RLAPA,RLAPB,
     $ TA,TB,NGrid,NDEN)
c
c    ******************************************************************
c    *                                                                *
c    *  reference:                                                    *
c    * T. Tsubeda, K. Hirao, Parameter-free exchange functional       *
c    * Phys. Rev.B  62, 15527 (200).                                  *
c    * Evaluating the BR89 hx on a grid                               *
c    ******************************************************************
c
      IMPLICIT NONE 
#include "fderiv1.inc"
      REAL*8 hx(NGrid,*),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
     $   D1RB(NGrid,3),RLAPA(NGrid),RLAPB(NGrid),TA(NGrid),TB(NGrid)
      REAL*8 D1,D2,DLAPA,DLAPB,TAUA,TAUB,gra,grb,Tol
      REAL*8 pi,pi13,pi23,ba,bb,aa,ab,gg1,gg2,ss
      REAL*8 SJ1a,SJ3a,SJ1b,SJ3b,ks1,ks2,gam1,gam2
      REAL*8 third
c
      INTEGER NGrid
      INTEGER i
      INTEGER NA,NB
      INTEGER NDEN
      parameter(pi=3.1415926535897930d0)
      parameter(third=1.0d0/3.0d0)
      Tol = 1.0D-12
      
      ! initilizing the d1vars
c
      DO i = 1,NGrid
        hx(i,1) = 0.0d0
        hx(i,2) = 0.0d0
c
       IF((RA(i).gt.Tol).OR.(RB(i).gt.Tol)) THEN
c
         IF (RA(i).gt.Tol) then
            D1 = RA(i)
            gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
            TAUA = TA(i)/2.0d0
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
            TAUB = TB(i)/2.0d0
            DLAPB = RLAPB(i)
         else
           D2 = 0.0d0
           grb=0.0d0
           TAUB = 0.0d0
           DLAPB = 0.0d0
         endif
c
         ks1 = (6.0d0*(pi**2)*D1)**third 
         IF((ks1.gt.Tol).and.D1.gt.Tol) then
          SJ1a = -(Cos(ks1*ss)/(ks1*ss))+Sin(ks1*ss)/((ks1**2)*ss**2)
          SJ3a = ((-15.0d0+(ks1**2)*ss**2)*Cos(ks1*ss))/
     $    ((ks1**3)*ss**3)- 
     $   ((-15.0d0 + 6.0d0*(ks1**2)*ss**2)*Sin(ks1*ss))/((ks1**4)*ss**4)
          gg1= DLAPA/4.0d0 + (3.0d0*(ks1**2)*D1)/5.0d0 - TAUA
          gam1 = (3.0d0*D1*SJ1a)/(ks1*ss)+(35.0d0*SJ3a*gg1)/
     $    (2.0d0*(ks1**3)*ss)
          hx(i,1) = - gam1*gam1/D1
         else
           SJ1a = 0.0d0
           SJ3a = 0.0d0
           hx(i,1) = 0.0d0
         endif
c
           IF(NDEN .eq. 1) then
            hx(i,2) = hx(i,1)
         else
           ks2 = (6.0d0*(pi**2)*D2)**third 
          IF((ks2.gt.Tol).and.D2.gt.Tol) then
           SJ1b = -(Cos(ks2*ss)/(ks2*ss))+Sin(ks2*ss)/((ks2**2)*ss**2)
         SJ3b = ((-15.0d0+(ks2**2)*ss**2)*Cos(ks2*ss))/((ks2**3)*ss**3)-
     $   ((-15.0d0 + 6.0d0*(ks2**2)*ss**2)*Sin(ks2*ss))/((ks2**4)*ss**4)
          gg2= DLAPB/4.0d0 + (3.0d0*(ks2**2)*D2)/5.0d0 - TAUB
          gam2 = (3.0d0*D2*SJ1b)/(ks2*ss)+(35.0d0*SJ3b*gg2)/
     $    (2.0d0*(ks2**3)*ss)
          hx(i,2) = - gam2*gam2/D2
         else
          SJ1b=0.0d0
          SJ3b = 0.0d0
          hx(i,2) = 0.0d0
         endif
         endif
c
c           write(6,*) "b_alpha = ", ba, "  b_beta= ", bb
c           write(6,*) "a_alpha = ", aa, "a_beta = ", ab
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
