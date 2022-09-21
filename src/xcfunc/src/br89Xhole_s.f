*  Deck br89Xhole_s 
*  this functional is contributing from Emil Proynov
*
       SUBROUTINE br89Xhole_s(hx,sVal,RA,RB,D1RA,D1RB,RLAPA,RLAPB,
     $ TA,TB,NGrid,NDEN)
c
c    ******************************************************************
c    *                                                                *
c    *  reference: A. D. Becke and E. R. Johnson, J. Chem. Phys. 122, *
c    *             154104 (2005).                                     *
c    * Evaluating the BR89 hx on a grid                               *
c    ******************************************************************
c
      IMPLICIT NONE 
#include "fderiv1.inc"
c     REAL*8 hx(NGrid,*),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
      REAL*8 hx(NGrid),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
     $   RLAPA(NGrid),RLAPB(NGrid),TA(NGrid),TB(NGrid),D1RB(NGrid,3)
      REAL*8 D1,D2,DLAPA,DLAPB,TAUA,TAUB,gra,grb,Tol,h1a,h1b,ddem,Tol2
      REAL*8 pi,pi13,pi23,ba,bb,aa,ab,sVal,ee1,ee2,ee3,ee4,h2a,h2b
c
      INTEGER NGrid
      INTEGER i
      INTEGER NA,NB
      INTEGER NDEN
      parameter(pi=3.1415926535897930d0)
      Tol = 1.0D-10
      Tol2 = 1.0D-07
 
c     write(*,*) "In br89 hole function"     
      ! initilizing the d1vars
c
      DO i = 1,NGrid
        ba = 0.0d0
        bb = 0.0d0
        aa= 0.0d0
        ab= 0.0d0
        hx(i) = 0.0d0
        gra = 0.0d0
        TAUA = 0.0d0
        DLAPA = 0.0d0
        grb = 0.0d0
        TAUB = 0.0d0
        DLAPB = 0.0d0
c
c       IF((RA(i).gt.Tol).OR.(RB(i).gt.Tol)) THEN
c
           IF (RA(i).gt.Tol) then
              D1 = RA(i)
              gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
              TAUA = TA(i)
              DLAPA = RLAPA(i)
c          CALL brbrxb_hx(ba,aa,D1,gra,DLAPA,TAUA)
c          endif
c
c          IF (RB(i).gt.Tol) then
c             D2 = RB(i)
c             grb = D1RB(i,1)**2 + D1RB(i,2)**2 + D1RB(i,3)**2
c             TAUB = TB(i)
c             DLAPB = RLAPB(i)
c          CALL brbrxb_hx(bb,ab,D2,grb,DLAPB,TAUB)
c          endif
c
           CALL brbrxb_hx(ba,aa,D1,gra,DLAPA,TAUA)
c
            ee1= DEXP(-aa*DSQRT((ba+sVal)**2))
            ee2= DEXP(-aa*DSQRT((ba-sVal)**2))
c
            h1a = ee2*(1.0d0+DSQRT((ba-sVal)**2))
            h2a= ee1*(1.0d0+aa*DSQRT((ba+sVal)**2))
            ddem = 16.0d0*pi*ba*sVal
c           IF(DABS(ddem).gt.Tol) then
            IF(DSQRT(ddem**2).gt.Tol2) then
              hx(i) = (h2a-h1a)*aa/ddem
            else
              hx(i) = 0.0d0
            endif
c           supress the spurous positive values. Normally you
c           leave the positive values there unless
c           you really want to supress them.
c           IF(hx(i).gt.-Tol*100.0d0) hx(i) = 0.0d0
c
c           write(6,*) "b_alpha = ", ba
c           write(6,*) "a_alpha = ", aa, "a_beta = ", ab
c           write(6,*) "rho_alpha = ",  RA(i)
c           write(6,*) "hx_alpha = ", hx(i)
c           write(6,*) "SValue = ", sVal 
c           write(6,*) "gra = ", gra
c           write(6,*) "tau = ", TAUA
c           write(6,*) "LAP = ", DLAPA
c           write(6,*) "hx_beta = ", hx(i,2)
c Here D1F includes only the derivatives of b2.
c
        ENDIF   !  RA,RB vs Tol
c
      ENDDO  !  NGRID
        RETURN
        END
c
       SUBROUTINE brbrxb_hx(b,aaa,DN,GR,DLAP,TAU)
c
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 Tol,DS
       REAL*8 DN,DLAP,TAU,GR
       REAL*8 x_y, g_y, P12_y, Q12_y, yy, QS, dlam
       INTEGER i,j,k,n,m,l
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      DN is the electron density for the input spin direction c
c      (i.e. either D1, or D2, but not the sum of the two);    c
c      GR is the dot product of the  density gradient with     c
c      iself (per sin direction), NOT the grad.modulus;        c
c      TAU is the electron kinetic energy density for the      c
c      given spin direction                                    c
c      DLAP is the Laplacian of the given input spin density   c
c                                                              c
c     OUTPUT                                                   c
c      EX exchange energy for the current spin direction       c
c      dfdr is the funct. derivative with resp. to rho_a(b)    c
c      dfdg is the funct. derivative with resp. to grad_a(b)   c
c      dfdt is the funct. derivative with resp. to tau_a(b)    c
c      dfdl is the funct. derivative with resp. to Lap_a(b)    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
       parameter(five=5.0d0,six=6.0d0,three=3.0d0)
       parameter(third=1.0d0/3.0d0)
       parameter(twthr=2.d0/3.d0,fiftr=5.d0/3.d0)
       parameter(pi=3.141592653589793d0)
c
       parameter(a1=1.525525181200953d0,a2=0.4576575543602858d0,
     $ a3=.4292036732051034d0)
c
       parameter(aa1=-31.887610400665510D0,aa2=2.50D0,
     $ aa3=6.250D0,aa4=27.18969924908615D0,aa5=32.23619130191664D0)
c
       parameter(c0= 0.7566445420735584D0,c1=-2.636397787137096D0,
     $ c2= 5.474515996423288D0,c3=-12.65730812710829D0,
     $ c4= 4.125058472512136D0,c5=-30.42513395716384D0)
c
       parameter(b0= 0.4771976183772063d0,b1=-1.779981349455627D0,
     $ b2=3.843384186230215D0,b3=-9.591205088051849D0,
     $ b4= 2.173018028591672D0,b5=-30.42513385160366D0)
c
       parameter(d0 = 0.00004435009886795587D0)
       parameter(d1 = 0.5812865360445791D0)
       parameter(d2 = 66.74276451594061D0)
       parameter(d3 = 434.2678089722977D0)
       parameter(d4 = 824.7765766052239D0)
       parameter(d5 = 1657.965273158212D0)
c
       parameter(e0 = 0.00003347285060926091D0)
       parameter(e1 = 0.4791793102397135D0)
       parameter(e2 = 62.39226833857424D0)
       parameter(e3 = 463.1481642793812D0)
       parameter(e4 = 785.2360350104029D0)
       parameter(e5 = 1657.962968223273D0)
c
       parameter(q0 =  0.08873042794156832D0)
       parameter(q1 = -0.5912993446836108D0)
       parameter(q2 =  4.039440258145228D0)
       parameter(q3 = -4.215855548751412D0)
       parameter(q4 =  41.53059614758940D0)
       parameter(q5 = -96.31041316013488D0)
       parameter(q6 =  136.8224866906499D0)
       parameter(q7 = -186.5731887073031D0)
       parameter(q8 =  59.39109156063976D0)
c
       parameter(s0 = -0.000001794332402415239D0)
       parameter(s1 = -0.001066065367041282D0)
       parameter(s2 = -4.304047130400146D0)
       parameter(s3 = -122.2881565440333D0)
       parameter(s4 = -4000.619157339993D0)
       parameter(s5 = -2575.440153194542D0)
       parameter(s6 =  19352.57170922590D0)
       parameter(s7 =  95767.25426932056D0)
       parameter(s8 = -65554.94378973276D0)
       parameter(dlam=1.0D0)
       parameter(AXC= 1.24070098179880D0)
       parameter(detol=1.0D-12,dtol=1.0D-8,dtol2=1.D-6)
       Tol = 1.0D-12
c
c Initializations start
          b = zero
          aaa=zero
          yy=zero
          x_y=zero
          g_y=zero
          P12_y=zero
          Q12_y=zero
c Initializations end
c
      IF((DN.GT.detol).and.(TAU.gt.detol)) then
c
c
           DN13 = DN**third
           DN43= DN13*DN
           DN53= DN*(DN13**2)
           DS = TAU - GR/(DN*four)
           QS = DLAP/six - dlam*DS/three
c          write(6,*) "Qs =    ", QS
           pi13=pi**third
           pi23 = pi13**2
c  yy is the argument y in the BR model being certain
c  function of rho, and Qs:
           IF(dabs(QS).gt.dtol) then
               yy = twthr*pi23*DN53/QS
           else
               yy=dtol2
           endif
c analytic part starts from here
        IF(yy.gt.dtol2) then
c  The upper branch of x_y : 0 < y <= +infinity
            t0 = yy**2
            t33 = aa3 + aa4*t0
            t3 = dsqrt(t33)
            t5 = dlog(aa2+t3)
            t7 = dlog(yy)
              g_y = aa1 + aa5 + t5 - t7
c
            t2 = yy**2
            t4 = t2*yy
            t6 = t2**2
            t8 = t6*yy
            t9 = e0+e1*yy+e2*t2+e3*t4+e4*t6+e5*t8
            P12_y = (d0+d1*yy+d2*t2+d3*t4+d4*t6+d5*t8)/t9
             x_y = g_y*P12_y
c
        elseif (yy.le.(-dtol)) then
c  The lower branch of x_y : -infinity <= y <= 0
            g_y = -atan(a1*yy+a2) +a3
            t2 = yy**2
            t4 = t2*yy
            t6 = t2**2
            t8 = t6*yy
            t9 = b0+b1*yy+b2*t2+b3*t4+b4*t6+b5*t8
            P12_y = (c0+c1*yy+c2*t2+c3*t4+c4*t6+c5*t8)/t9
             x_y=g_y*P12_y
c
        else
           x_y=zero
           yy=zero
        endif  !  yy
c
           cc = two*pi13*DN43
           cc1=twthr*pi13*DN13
           cc2=two*pi13*DN13
           XEXP=DEXP(-x_y/3.0d0)
c
c analytic part ends here
           b = 0.50d0*x_y*XEXP/pi13/DN13
           IF(DABS(b).gt.Tol) then
            aaa = x_y/b
           else
            aaa = 0.0d0
           ENDIF
c
      ENDIF    !  DN and TAU vs detol
      RETURN
      END
