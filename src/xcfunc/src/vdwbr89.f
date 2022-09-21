*  Deck br89
*  this functional is contributing from Emil Proynov
*
       SUBROUTINE br89b2_vdwx(b2,D1F,INFOR,RA,RB,D1RA,D1RB,RLAPA,RLAPB,
     $ TA,TB,NGrid,NDEN)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the WdV exchange model (the b2 term) and its        *
c    *  derivatives on a grid using a home-made analytical            *
c    * interpolation of the x(y) function.                            *
c    *                                                                *
c    *  reference: A. D. Becke and E. R. Johnson, J. Chem. Phys. 122, *
c    *             154104 (2005).                                     *
c    *  the analytic derivatives of b2 term is from the below         *
c    *  reference:                                                    *
c    *  Analytical representation of the Becke-Roussel exchange       *
c    *  functional, Chem Phys Lett. 2008 March 31; vol 455,           *
c    *  issue 1-3, page 103–109.                                      *
c    *  Emil Proynov, Zhenting Gan, Jing Kong                         *
c    *                                                                *
c    *  OUTPUT:                                                       *
c    *     F      - Functional values                                 *
c    *     D1F    - First derivatives                                 *
c    *                                                                *
c    *  INPUT:                                                        *
c    *     RA,B   - Spin densities                                    *
c    *     D1RA,B   - Spin densitiy gradients                         *
c    *     TA,B   - Spin kinetic energy densities                     *
c    *     DLAPA,B - Laplacian of the spin densities                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
#include "fderiv1.inc"
      REAL*8 b2(NGrid,*),D1F(NGrid,*),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
     $D1RB(NGrid,3),RLAPA(NGrid),RLAPB(NGrid),TA(NGrid),TB(NGrid)
      REAL*8 dfdra,dfdrb,dfdga,dfdgb,dfdta,dfdtb,dfdla,dfdlb
      REAL*8 EX1,EX2,D1,D2,DLAPA,DLAPB,TAUA,TAUB,gra,grb,Tol
      REAL*8 pi,pi13,pi23,ba,bb
c
      INTEGER NGrid
      INTEGER i
      INTEGER NA,NB
      INTEGER NDEN
      INTEGER INFOR(*)
      INTEGER D1VARS(N_FUNC_DERIV_1)
      Tol = 1.0D-12
      
      ! initilizing the d1vars
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)

      DO i = 1,NGrid
        ba = 0.0d0
        bb = 0.0d0
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
c  Call first for the spin-up components
           CALL brbrxb_vdwx(ba,D1,gra,DLAPA,TAUA,dfdra,
     $               dfdga,dfdta,dfdla)
c
           IF(NDEN .eq. 1) then
             dfdrb=dfdra
             dfdgb=dfdga
             dfdtb=dfdta
             dfdlb=dfdla
             bb=ba
           else
c  Call next for the spin-down component
            CALL brbrxb_vdwx(bb,D2,grb,DLAPB,TAUB,dfdrb,
     $      dfdgb,dfdtb,dfdlb)
           endif
c
            b2(i,1) = ba*ba
            b2(i,2) = bb*bb
c Here D1F includes only the derivatives of b2.
c
            ID_RA_POS         = D1VARS(ID_RA)
            ID_RB_POS         = D1VARS(ID_RB)
            D1F(i,ID_RA_POS)  = dfdra
            D1F(i,ID_RB_POS)  = dfdrb
c
            ID_GAA_POS        = D1VARS(ID_GAA)
            ID_GAB_POS        = D1VARS(ID_GAB)
            ID_GBB_POS        = D1VARS(ID_GBB)
            D1F(i,ID_GAA_POS) = dfdga
            D1F(i,ID_GBB_POS) = dfdgb
            D1F(i,ID_GAB_POS) = 0.0d0
c
            ID_TA_POS         = D1VARS(ID_TA)
            ID_TB_POS         = D1VARS(ID_TB)
            ID_LA_POS         = D1VARS(ID_LA)
            ID_LB_POS         = D1VARS(ID_LB)
            D1F(i,ID_TA_POS)  = dfdta
            D1F(i,ID_TB_POS)  = dfdtb
            D1F(i,ID_LA_POS)  = dfdla
            D1F(i,ID_LB_POS)  = dfdlb
c
        ENDIF   !  RA,RB vs Tol
c
      ENDDO  !  NGRID
        RETURN
        END
c
       SUBROUTINE brbrxb_vdwx(b,DN,GR,DLAP,TAU,dfdr,dfdg,dfdt,dfdl)
c
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 dfdr,dfdg,dfdt,dfdl,Tol,DS,dxdy,dydr,Ux, dUdx
       REAL*8 DN,DLAP,TAU,GR,dUdr,dydg,dydt,dydl
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
c   Two possible choices of Lambda:
       parameter(dlam=1.0D0)
c      parameter(dlam=0.80D0)
       parameter(AXC= 1.24070098179880D0)
       parameter(detol=1.0D-12,dtol=1.0D-8,dtol2=1.D-6)
       Tol = 1.0D-12
c
c Initializations start
          dfdr=zero
          dfdg=zero
          dfdt=zero
          dfdl=zero
          b = zero
          yy=zero
          x_y=zero
          g_y=zero
          P12_y=zero
          Q12_y=zero
          dxdy=zero
          Ux=zero
          dUdx=zero
          dUdr = zero
          dydr=zero
          dydg=zero
          dydt=zero
          dydl=zero
          U12 = zero
          V23 = zero
          U3 = zero
c Initializations end
c
      IF((DN.GT.detol).and.(TAU.gt.detol)) then
c
c    First calculate the BR exchange energy density
c
           DN13 = DN**third
           DN43= DN13*DN
           DN53= DN*(DN13**2)
           DS = two*TAU - GR/(DN*four)
           QS = DLAP/six - dlam*DS/three
           pi13=pi**third
           pi23 = pi13**2
c  yy is the argument y in the BR model being certain
c  function of rho, and Qs:
           IF(abs(QS).gt.dtol) then
               yy = twthr*pi23*DN53/QS
           else
               yy=dtol2
           endif
c analytic part starts from here
        IF(yy.gt.dtol2) then
c  The upper branch of x_y : 0 < y <= +infinity
            t0 = yy**2
            t33 = aa3 + aa4*t0
            t3 = sqrt(t33)
            t5 = log(aa2+t3)
            t7 = log(yy)
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
            t10 = yy*t3
            dxx = aa2/t10
            tp2 = yy**2
            tp4 = tp2*yy
            tp6 = tp2**2
            tp8 = tp6*yy
            tp14 = tp6**2
            tp23 = t9**2
            Q12_y = (s0+s1*yy+s2*tp2+s3*tp4+s4*tp6+s5*tp8+s6*tp6*tp2
     $             + s7*tp6*tp4+s8*tp14)/tp23
            dxdy= -dxx*P12_y + g_y*Q12_y
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
            t92 = (a1*yy + a2)** 2
            ddx = a1/(one + t92)
            ts2 = yy**2
            ts4 = ts2*yy
            ts6 = ts2**2
            ts8 = ts6*yy
            ts14 = ts6**2
            ts23 = t9**2
            Q12_y = (q0+q1*yy+q2*ts2+q3*ts4+q4*ts6+q5*ts8+q6*ts6*ts2
     $            + q7*ts6*ts4+q8*ts14)/ts23
            dxdy= -ddx*P12_y + g_y*Q12_y
        else
           x_y=zero
           dxdy=zero
           yy=zero
        endif  !  yy
c
           cc = two*pi13*DN43
           cc1=twthr*pi13*DN13
           cc2=two*pi13*DN13
           XEXP=exp(-x_y/3.0d0)
c
c analytic part ends here
           WDV = XEXP*(x_y-three)/pi13/(six*DN13)
           b = 0.50d0*x_y*XEXP/pi13/DN13
c
c  Now start with the FBA of the functional derivative
              tr1 = pi13
              tr2 = tr1**2
              DN23 = DN13*DN13
        IF((x_y.gt.dtol).and.(abs(QS).gt.dtol)) then
c
              QS2 = QS**2
              dydr = 10.0D0*tr2*DN23/(QS*9.0d0)
     $             + tr2*dlam*GR/(18.0d0*DN13*QS2)
c
              dydg = - dlam*pi23*DN23/(QS2*18.0d0)
              dydt = four*dlam*pi23*DN53/(QS2*9.0d0)
              dydl = - pi23*DN53/(QS2*9.0d0)
              dbdr = -XEXP*x_y/(pi13*six*DN43)-WDV*dxdy*dydr
              dfdr = two*b*dbdr
c
              dbdg = -WDV*dxdy*dydg
              dfdg = two*b*dbdg
c
              dbdt = -WDV*dxdy*dydt
              dfdt = two*b*dbdt
c
              dbdl = -WDV*dxdy*dydl
              dfdl = two*b*dbdl
c
        endif   !  abs(QS)
c
c     This expression is simplified above
c       b = 0.125D0*(x_y**3)*exp(-x_y)/(PI*DN)
c       b = b**(1D0/3D0)

      ENDIF    !  DN and TAU vs detol
      RETURN
      END
