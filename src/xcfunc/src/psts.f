!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
!
! This file describes the second term for the PSTS functional, which is:
! \int [1-a(r)]*[e^tpss_x(r) - e^ex_x(r) dr]
! a(r) is some complicated functional, that's why we use several
! subroutines to describe it.
! a1 is used to  evaluate the abnormal regions that exchange 
! can dominate; the a2 is used to characterizes the abnormal 
! region that the hole density does not integrate to -1      
! finally the a(r) is:
! 1-a(r) = (1-a1)*(1-a2)      
! see the paper below for more information:      
! Density functional with full exact exchange, balanced nonlocality 
! of correlation, and constraint satisfaction
! Perdew, J.P. and Staroverov, V.N. and Tao, J. and Scuseria, G.E.
!
!
! note:
!      
! psts functional uses both small tau(tau with 1/2) and small exchange
! energy density (exa and exb are with 1/2). However, our program
! always use the big tau and big exchange energy density (both without
! 1/2). for tau variable, since PSTS kernel is TPSS it's already been
! handled over there. For the exchange energy density, you can see that
! we handled it in the a2(r) functional calculation and we also multiple
! 1/2 to the input ex density in the input of psts functional. here we
! just make an note on that.
!
!
!
!

      subroutine pstsa1(NDEN,TOL,RA,RB,DrhoAX,DrhoAY,DrhoAZ,
     &                  DrhoBX,DrhoBY,DrhoBZ,TA,TB,F,D1F)
      IMPLICIT NONE
!--------------------------------------------------------------------
! This is the a1(r) functional      
! INPUT :
! RA         : the alpha electron density
! RB         : the beta  electron density
! DRhoA      : the alpha rho'
! DRhoB      : the beta  rho'
! TA         : the alpha kinetic energy density
! TB         : the beta  kinetic energy density
! NDEN       : number of density
! TOL        : tolerance value
! OUTPUT:
! F          : functional values
! D1F        : the first  order functional derivatives
!--------------------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 
      ! for the a1(r), we do not need to input the infor arary
      ! we will define it inside
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)

      ! for the functional a1(r), it's a META-GGA functional
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS
      INTEGER  I

      ! input and output
      INTEGER NDEN
      REAL*8 RA,RB
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 TA,TB
      REAL*8 F, D1F(N_FUNC_DERIV_1)

      ! pointers for tpss correlation functional
      REAL*8 F_tpss
      REAL*8 D1F_tpss(N_FUNC_DERIV_1)

      ! functional value and 1st derivatives for lsd
      REAL*8 l,lra,lrb

      ! functional value and 1st derivatives for tpss
      REAL*8 p,pra,prb,pgaa,pgab,pgbb,pta,ptb

      ! functional value and 1st derivatives for u
      ! u = e_gltpssc/e_lsd
      REAL*8 u,ura,urb,ugaa,ugab,ugbb,uta,utb

      ! variables and constants
      REAL*8 A,B,F1,ZERO
      REAL*8 t,t1,t2
      REAL*8 TOL,TOL2

      ! constants
      A    = 3.74D0
      B    = 167D0
      F1   = 1.0D0
      ZERO = 0.0D0
      TOL2 = sqrt(TOL)

      ! firstly initilize variable position information
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR(i) = -1
      END DO
      VAR_INFOR(ID_RHO)   = 1
      VAR_INFOR(ID_GAMMA) = 1
      VAR_INFOR(ID_TAU)   = 1
      CALL INIT_FUNC_DERIV_1(VAR_INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_TB_POS  = D1VARS(ID_TB)

      ! initialize the functional derivatives 
      F = 0.0D0
      DO I = 1, N_FUNC_DERIV_1
         D1F(I) = 0.0D0
      END DO

      !----------------------------------------------------------------
      ! real work begins.....
      ! we note that RA > TOL and RB < TOL or RA < TOL and RB > TOL
      ! have already been considered in the LSD functional and 
      ! GL_TPSS functional, so we do not consider it here
      !----------------------------------------------------------------
      IF (RA .GT. TOL .OR. RB .GT. TOL) THEN

         ! get the lsd and tpss functional infor
         F_tpss  = 0.0D0
         CALL lsd_pointwise(TOL,RA,RB,l,lra,lrb)
         CALL tpssc_gllimit(NDEN,TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,TA,TB,F_tpss,D1F_tpss)
         p    = F_tpss
         pra  = D1F_tpss(ID_RA_POS)
         prb  = D1F_tpss(ID_RB_POS)
         pgaa = D1F_tpss(ID_GAA_POS)
         pgab = D1F_tpss(ID_GAB_POS)
         pgbb = D1F_tpss(ID_GBB_POS)
         pta  = D1F_tpss(ID_TA_POS)
         ptb  = D1F_tpss(ID_TB_POS)

         ! to calculate u
         ! here we should have more protection for the l
C         write(6,*)"p is", p
C         write(6,*)"l is", l
         IF (abs(l) .GT. TOL) THEN
            u    = p/l 
            ura  = (pra*l - lra*p)/(l*l) 
            ugaa = pgaa/l
            ugab = pgab/l
            uta  = pta/l

            ! consider the beta part of u
            IF (NDEN .EQ. 1) THEN
               urb   = ura
               ugbb  = ugaa
               utb   = uta
            ELSE
               urb  = (prb*l - lrb*p)/(l*l) 
               ugbb = pgbb/l
               utb  = ptb/l
            END IF

            t   =  F1 + A*dlog(abs(F1+B*u))
            IF (abs(t) .GT. TOL2) THEN
               t1  =  -F1*(F1/(t*t))
               t2  =  A*B/(F1+B*u)
               F   =  F + F1/t
               D1F(ID_RA_POS)  = D1F(ID_RA_POS)  + t1*t2*ura 
               D1F(ID_GAA_POS) = D1F(ID_GAA_POS) + t1*t2*ugaa 
               D1F(ID_GAB_POS) = D1F(ID_GAB_POS) + t1*t2*ugab
               D1F(ID_TA_POS)  = D1F(ID_TA_POS)  + t1*t2*uta

               ! consider the beta part 
               IF (NDEN .EQ. 1) THEN
                  D1F(ID_RB_POS)  = D1F(ID_RA_POS) 
                  D1F(ID_GBB_POS) = D1F(ID_GAA_POS)
                  D1F(ID_TB_POS)  = D1F(ID_TA_POS) 
               ELSE
                  D1F(ID_RB_POS)  = D1F(ID_RB_POS)  + t1*t2*urb 
                  D1F(ID_GBB_POS) = D1F(ID_GBB_POS) + t1*t2*ugbb 
                  D1F(ID_TB_POS)  = D1F(ID_TB_POS)  + t1*t2*utb
               END IF

            ELSE

               write(6,*)" the vale of t is:", t
               write(6,*)" the t value is too small in functional_psts"
               write(6,*)" therefore we just set all output values to 0"
               F               =  0.0D0
               D1F(ID_RA_POS)  =  0.0D0
               D1F(ID_GAA_POS) =  0.0D0
               D1F(ID_GAB_POS) =  0.0D0
               D1F(ID_TA_POS)  =  0.0D0
               D1F(ID_RB_POS)  =  0.0D0
               D1F(ID_GBB_POS) =  0.0D0
               D1F(ID_TB_POS)  =  0.0D0

            END IF  ! t^2 > tol

         ELSE

            F               =  0.0D0
            D1F(ID_RA_POS)  =  0.0D0
            D1F(ID_GAA_POS) =  0.0D0
            D1F(ID_GAB_POS) =  0.0D0
            D1F(ID_TA_POS)  =  0.0D0
            D1F(ID_RB_POS)  =  0.0D0
            D1F(ID_GBB_POS) =  0.0D0
            D1F(ID_TB_POS)  =  0.0D0

         END IF  ! l > TOL


C      write(6,*)"a1 is:",F
C      write(6,*)"a1RA is:",D1F(ID_RA_POS)
C      write(6,*)"a1RB is:",D1F(ID_RB_POS)
C      write(6,*)"a1GAA is:",D1F(ID_GAA_POS)
C      write(6,*)"a1GAB is:",D1F(ID_GAB_POS)
C      write(6,*)"a1GBB is:",D1F(ID_GBB_POS)
C      write(6,*)"a1TA  is:",D1F(ID_TA_POS)
C      write(6,*)"a1TB  is:",D1F(ID_TB_POS)

      ELSE

         F               =  0.0D0
         D1F(ID_RA_POS)  =  0.0D0
         D1F(ID_GAA_POS) =  0.0D0
         D1F(ID_GAB_POS) =  0.0D0
         D1F(ID_TA_POS)  =  0.0D0
         D1F(ID_RB_POS)  =  0.0D0
         D1F(ID_GBB_POS) =  0.0D0
         D1F(ID_TB_POS)  =  0.0D0

      END IF ! RA > TOL or RB > TOL

C      write(6,*)"a1 is:",F
C      write(6,*)"a1RA is:",D1F(ID_RA_POS)
C      write(6,*)"a1RB is:",D1F(ID_RB_POS)
C      write(6,*)"a1GAA is:",D1F(ID_GAA_POS)
C      write(6,*)"a1GAB is:",D1F(ID_GAB_POS)
C      write(6,*)"a1GBB is:",D1F(ID_GBB_POS)
C      write(6,*)"a1TA  is:",D1F(ID_TA_POS)
C      write(6,*)"a1TB  is:",D1F(ID_TB_POS)
      END

      SUBROUTINE lsd_pointwise(TOL,RA,RB,F,D1FRA,D1FRB)
c     ******************************************************************
c     *  This subroutine is used to calculate the LSD exchange energy  *
c     *  per electron. See the paper of psts.                          *
c     *                                                                *
c     *                                                                *
c     *  OUTPUT:                                                       *
c     *     F      - Functional values                                 *
c     *     D1F    - First derivatives                                 *
c     *                                                                *
c     *  INPUT:                                                        *
c     *     RA,B   - Spin densities                                    *
c     *                                                                *
c     *                                                                *
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 F,D1FRA,D1FRB,RA,RB,TOL,PI

C     constants
      PI  = 4.D0*DATAN(1.D0)
      F1  = 1.0D0
      F2  = 2.0D0
      F3  = 3.0D0
      F4  = 4.0D0
      F13 = F1/F3
      F43 = F4/F3
      F34 = F3/F4
      F12 = F1/F2
      Cs  = -(F34)*F12*(F3/PI)**F13

c     Do the real work
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN
         RA13  = (F2*RA)**F13
         RB13  = (F2*RB)**F13
         RA43  = RA13*(F2*RA)
         RB43  = RB13*(F2*RB)
         RAB   = RA+RB
         F     = Cs*(RA43+RB43)/RAB
         D1FRA = Cs*(F43*F2*RA13*RAB-(RA43+RB43))/(RAB**2)
         D1FRB = Cs*(F43*F2*RB13*RAB-(RA43+RB43))/(RAB**2)
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
         RA13  = (F2*RA)**F13
         RA43  = RA13*(F2*RA)
         RAB   = RA
         F     = Cs*RA43/RAB
         D1FRA = Cs*(F43*F2*RA13*RAB-RA43)/(RAB**2)
         D1FRB = 0.0D0
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
         RB13  = (F2*RB)**F13
         RB43  = RB13*(F2*RB)
         RAB   = RB
         F     = Cs*RB43/RAB
         D1FRA = 0.0D0
         D1FRB = Cs*(F43*F2*RB13*RAB-RB43)/(RAB**2)
      ELSE
         F     = 0.0D0
         D1FRA = 0.0D0
         D1FRB = 0.0D0
      END IF
      RETURN
      END

      subroutine pstsa2(NDEN,TOL,RA,RB,DrhoAX,DrhoAY,DrhoAZ,
     &                  DrhoBX,DrhoBY,DrhoBZ,TA,TB,
     &                  exa,exb,D1F_excha,D1F_exchb,
     &                  F,D1F)
      IMPLICIT NONE
!--------------------------------------------------------------------
! This is the a2(r) for the PSTS functional      
! INPUT :
! RA         : the alpha electron density
! RB         : the beta  electron density
! DRhoA      : the alpha rho'
! DRhoB      : the beta  rho'
! TA         : the alpha kinetic energy density
! TB         : the beta  kinetic energy density
! exa        : the alpha exchange energy density      
! exb        : the beta  exchange energy density      
! OUTPUT:
! F          : functional values
! D1F        : the first  order functional derivatives
!
! here it's worhty to make a note. The derivatives for exa and exb
! is passed out explicitly, therefore we assuem the a2(r) functional
! is only a "META-GGA" functional. This is the trick we need to make
! a reminder.
!
!--------------------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 
      ! for the a1(r), we do not need to input the infor arary
      ! we will define it inside
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)

      ! for the functional a1(r), it's a META-GGA functional
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS
      INTEGER  I

      ! input and output
      INTEGER NDEN
      REAL*8 RA,RB
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 TA,TB
      REAL*8 exa,exb
      REAL*8 D1F_excha, D1F_exchb
      REAL*8 F, D1F(*)

      REAL*8 ex  ! exa + exb

      ! pointers for tpss exchange functional
      REAL*8 F_tpss
      REAL*8 D1F_tpss(N_FUNC_DERIV_1)

      ! functional value and 1st derivatives for tpss
      REAL*8 p,pra,prb,pgaa,pgab,pgbb,pta,ptb

      ! functional value and 1st derivatives for v
      ! v = e^ex_x/e^tpss_x
      ! vexa and vexb is the derivatives for exa and exb
      REAL*8 v,vra,vrb,vgaa,vgab,vgbb,vta,vtb,vexa,vexb

      ! functional value and 1st derivatives for f(v)
      REAL*8 fv,fvv,fvra,fvrb,fvgaa,fvgab,fvgbb,fvta,fvtb,fvexa,fvexb
      REAL*8 p1,pc

      ! functinal and its derivatives for g(x)
      REAL*8 g,gra,grb

      ! variables and constants
      REAL*8 C,D,E,FF,F1,F2,F3
      REAL*8 TOL,TOL2

      !----------------------------------------------------------------
      ! real work begins, preparation
      !----------------------------------------------------------------
      ! constants
      F1  = 1.0D0
      F2  = 2.0D0
      F3  = 3.0D0
      C   = 0.909D0
      D   = 7.100D0
      E   = 9.610D0
      FF  = -(F3/F2)*(F1/dlog((F1-C)/F2))
      TOL2 = sqrt(TOL)

      ! firstly initilize variable position information
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR(i) = -1
      END DO
      VAR_INFOR(ID_RHO)   = 1
      VAR_INFOR(ID_GAMMA) = 1
      VAR_INFOR(ID_TAU)   = 1
      CALL INIT_FUNC_DERIV_1(VAR_INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_TB_POS  = D1VARS(ID_TB)

      ! initialize the functional derivatives 
      F = 0.0D0
      DO I = 1, N_FUNC_DERIV_1
         D1F(I) = 0.0D0
      END DO

      !----------------------------------------------------------------
      ! real work begins.....
      !----------------------------------------------------------------
      IF (RA .GT. TOL .OR. RB .GT. TOL) THEN


         !--------------------------------------------------------------
         ! get the lsd and tpss functional infor
         !--------------------------------------------------------------
         ex = exa + exb
         F_tpss  = 0.0D0
         IF (NDEN .EQ. 1) THEN
            CALL functional_tpssx_pw_close(TOL,RA,DRhoAX,DRhoAY,
     &           DRhoAZ,TA,F_tpss,D1F_tpss)
         ELSE
            CALL functional_tpssx_pw_open(TOL,RA,RB,DRhoAX,DRhoAY,
     &           DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,TA,TB,F_tpss,D1F_tpss)
         END IF
         p    = F_tpss
         pra  = D1F_tpss(ID_RA_POS)
         prb  = D1F_tpss(ID_RB_POS)
         pgaa = D1F_tpss(ID_GAA_POS)
         pgab = D1F_tpss(ID_GAB_POS)
         pgbb = D1F_tpss(ID_GBB_POS)
         pta  = D1F_tpss(ID_TA_POS)
         ptb  = D1F_tpss(ID_TB_POS)

         !--------------------------------------------------------------
         ! calculate g(x) and its derivatives
         !--------------------------------------------------------------
         CALL functional_pstsa2gx(TOL,D,E,RA,RB,g,gra,grb)

         !--------------------------------------------------------------
         ! calculate v
         ! pay attention to the trap here!!!
         ! usually, the ex and p should be very close, so we test the v
         !--------------------------------------------------------------
         IF (abs(p) .GT. TOL) THEN
            v    = ex/p 
            IF (v .GE. F1) THEN
               ! in this case, the f(v)
               ! is zero, hence a2 is zero
               F               = 0.0D0
               D1F(ID_RA_POS)  = 0.0D0
               D1F(ID_RB_POS)  = 0.0D0
               D1F(ID_GAA_POS) = 0.0D0
               D1F(ID_GAB_POS) = 0.0D0
               D1F(ID_GBB_POS) = 0.0D0
               D1F(ID_TA_POS)  = 0.0D0
               D1F(ID_TB_POS)  = 0.0D0
               D1F_excha       = 0.0D0
               D1F_exchb       = 0.0D0
               go to 10
            ELSE IF (v .LE. C) THEN
               ! in this case, f(v) = 1 hence it turns out to be 
               ! a constant,consider the g(x) function; all of 
               ! other derivatives will be zero except the density
               fv   =  F1
               fvra =  0.0D0
               fvrb =  0.0D0
               fvgaa=  0.0D0 
               fvgab=  0.0D0 
               fvgbb=  0.0D0 
               fvta =  0.0D0 
               fvtb =  0.0D0 
               fvexa=  0.0D0 
               fvexb=  0.0D0 
            ELSE
               ! for this case, we have to calculate the v which is
               ! used in the f(v). Here we have p approx ex
               vra  = -(ex/(p*p))*pra
               vgaa = -(ex/(p*p))*pgaa
               vgab = -(ex/(p*p))*pgab
               vta  = -(ex/(p*p))*pta
               vexa =  F1/p

               ! consider the beta part of v
               IF (NDEN .EQ. 1) THEN
                  vrb   = vra
                  vgbb  = vgaa
                  vtb   = vta
                  vexb  = vexa
               ELSE
                  vrb   = -(ex/(p*p))*prb
                  vgbb  = -(ex/(p*p))*pgbb
                  vtb   = -(ex/(p*p))*ptb
                  vexb  =  F1/p
               END IF

               ! next, let us calculate f(v)
               p1    = F1/(F1-v)**FF
               pc    = F1/(v-C)**FF
               IF (p1 .LE. pc) THEN
                  CALL pstsa2_fv1(C,FF,v,fv,fvv)
               ELSE
                  CALL pstsa2_fv2(C,FF,v,fv,fvv)
               END IF
               fvra  = fvv*vra
               fvrb  = fvv*vrb
               fvgaa = fvv*vgaa
               fvgab = fvv*vgab
               fvgbb = fvv*vgbb
               fvta  = fvv*vta
               fvtb  = fvv*vtb
               fvexa = fvv*vexa
               fvexb = fvv*vexb
            END IF  ! if v<=c ...
         ELSE

            ! here p is too small, we consider not to do the 
            ! following calculation anymore
            F               = 0.0D0
            D1F(ID_RA_POS)  = 0.0D0
            D1F(ID_RB_POS)  = 0.0D0
            D1F(ID_GAA_POS) = 0.0D0
            D1F(ID_GAB_POS) = 0.0D0
            D1F(ID_GBB_POS) = 0.0D0
            D1F(ID_TA_POS)  = 0.0D0
            D1F(ID_TB_POS)  = 0.0D0
            D1F_excha       = 0.0D0
            D1F_exchb       = 0.0D0
            go to 10
         END IF ! if p < TOL

         !--------------------------------------------------------------
         ! Finally, calculate the F and D1F
         ! F = g(x)*f(v)
         ! if we step here, we must have v <=C or C<v<1
         !--------------------------------------------------------------
         F               = g*fv
         D1F(ID_RA_POS)  = gra*fv + g*fvra
         D1F(ID_RB_POS)  = grb*fv + g*fvrb
         D1F(ID_GAA_POS) = g*fvgaa 
         D1F(ID_GAB_POS) = g*fvgab
         D1F(ID_GBB_POS) = g*fvgbb
         D1F(ID_TA_POS)  = g*fvta
         D1F(ID_TB_POS)  = g*fvtb
         D1F_excha       = g*fvexa
         D1F_exchb       = g*fvexb

      ELSE

         F               = 0.0D0
         D1F(ID_RA_POS)  = 0.0D0
         D1F(ID_RB_POS)  = 0.0D0
         D1F(ID_GAA_POS) = 0.0D0
         D1F(ID_GAB_POS) = 0.0D0
         D1F(ID_GBB_POS) = 0.0D0
         D1F(ID_TA_POS)  = 0.0D0
         D1F(ID_TB_POS)  = 0.0D0
         D1F_excha       = 0.0D0
         D1F_exchb       = 0.0D0

      END IF  ! if RA > TOL and RB > TOL


10    continue
C      write(6,*)"a2 is:",F
C      write(6,*)"a2RA is:",D1F(ID_RA_POS)
C      write(6,*)"a2RB is:",D1F(ID_RB_POS)
C      write(6,*)"a2GAA is:",D1F(ID_GAA_POS)
C      write(6,*)"a2GAB is:",D1F(ID_GAB_POS)
C      write(6,*)"a2GBB is:",D1F(ID_GBB_POS)
C      write(6,*)"a2TA  is:",D1F(ID_TA_POS)
C      write(6,*)"a2TB  is:",D1F(ID_TB_POS)
C      write(6,*)"exa2  is:",D1F_excha
C      write(6,*)"exb2  is:",D1F_exchb
      END

      subroutine functional_pstsa2gx(TOL,D,E,rhoA,rhoB,F,D1FRA,D1FRB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION rhoA
      DOUBLE PRECISION RA, RB 
      DOUBLE PRECISION rhoB
      DOUBLE PRECISION F,D1FRA,D1FRB
      DOUBLE PRECISION TOL ! tolerance
      DOUBLE PRECISION PI
      DOUBLE PRECISION D
      DOUBLE PRECISION E
      
      ! firstly deal with the constants
      PI  = 4.D0*DATAN(1.D0)

      ! clear the result array
      F     = 0.0D0
      D1FRA = 0.0D0
      D1FRB = 0.0D0
      
      ! the real calculation begins
      RA = rhoA
      RB = rhoB
      IF (RB.LE.TOL) THEN
         RB = 0.0D0 
      END IF 
      IF (RA.LE.TOL) THEN
         RA = 0.0D0 
      END IF 
      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN
         RETURN 
      END IF 
             
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN
         t3 = (RA-1.D0*RB)**2
         t5 = RA+RB
         t6 = t5**2
         t9 = (PI*t5)**(1.D0/3.D0)
         t10 = 1/t6*t9
         t18 = 0.1100642416298209D1*D*t3*t10/(1.D0+0.1100642416298209D1
     &*E*t3*t10)
         F =   t18 

         t2 = RA-1.D0*RB
         t4 = RA+RB
         t5 = t4**2
         t6 = 1/t5
         t8 = (PI*t4)**(1.D0/3.D0)
         t9 = t6*t8
         t10 = t2**2
         t11 = E*t10
         t14 = 1.D0+0.1100642416298209D1*t11*t9
         t15 = 1/t14
         t19 = D*t10
         t22 = 1/t5/t4*t8
         t26 = t19*t6
         t27 = t8**2
         t28 = 1/t27
         t33 = t14**2
         t49 = 0.2201284832596418D1*D*t2*t9*t15-0.2201284832596418D1
     &*t19*t22*t15+0.3668808054327363D0*t26*t28*t15*PI
     &-0.1100642416298209D1*t26
     &*t8/t33*(0.2201284832596418D1*E*t2*t9-0.2201284832596418D1*t11*t22
     &+0.3668808054327363D0*t11*t6*t28*PI)
         D1FRA =    t49 

         t2 = RA-1.D0*RB
         t4 = RA+RB
         t5 = t4**2
         t6 = 1/t5
         t8 = (PI*t4)**(1.D0/3.D0)
         t9 = t6*t8
         t10 = t2**2
         t11 = E*t10
         t14 = 1.D0+0.1100642416298209D1*t11*t9
         t15 = 1/t14
         t19 = D*t10
         t22 = 1/t5/t4*t8
         t26 = t19*t6
         t27 = t8**2
         t28 = 1/t27
         t33 = t14**2
         t49 = -0.2201284832596418D1*D*t2*t9*t15
     &-0.2201284832596418D1*t19*t22*t15+0.3668808054327363D0*t26*t28*
     &t15*PI-0.1100642416298209D1*t2
     &6*t8/t33*(-0.2201284832596418D1*E*t2*t9-0.2201284832596418D1*t11*t
     &22+0.3668808054327363D0*t11*t6*t28*PI)
         D1FRB =   t49 

      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
         t2 = (PI*RA)**(1.D0/3.D0)
         t9 = 0.1100642416298209D1*D*t2/(1.D0+0.1100642416298209D1*E*t2)
         F =     t9 

         t2 = (PI*RA)**(1.D0/3.D0)
         t3 = t2**2
         t8 = 1.D0+0.1100642416298209D1*E*t2
         t15 = t8**2
         t21 = 0.3668808054327363D0*D/t3/t8*PI-0.4038045761849199D0*D/t2
     &/t15*E*PI
         D1FRA =   t21 

      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
         t2 = (PI*RB)**(1.D0/3.D0)
         t9 = 0.1100642416298209D1*D*t2/(1.D0+0.1100642416298209D1*E*t2)
         F =   t9 

         t2 = (PI*RB)**(1.D0/3.D0)
         t3 = t2**2
         t8 = 1.D0+0.1100642416298209D1*E*t2
         t15 = t8**2
         t21 = 0.3668808054327363D0*D/t3/t8*PI-0.4038045761849199D0*D/t2
     &/t15*E*PI
         D1FRB =  t21 
      END IF
      RETURN
      END

      subroutine pstsa2_fv1(C,F,V,FV,FVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------------------
! This subroutine is used to calculate the f(v) in the a2(r)
! the expression is the first one is the B2 (appendix B)      
! It calculates two terms, on is f(v), the other is derivatives of
! f(v) with respect to the v      
!--------------------------------------------------------------------
      REAL*8 C,F,V,FV,FVV


      t3 = (1.D0-1.D0*v)**F
      t7 = (v-1.D0*C)**F
      t11 = dexp(1/t3-1.D0/t7)
      t13 = 1/(1.D0+t11)
      FV  = t13
      t2 = 1.D0-1.D0*v
      t3 = t2**F
      t4 = 1/t3
      t6 = v-1.D0*C
      t7 = t6**F
      t8 = 1/t7
      t11 = dexp(t4-1.D0*t8)
      t13 = (1.D0+t11)**2
      t25 = -1.D0/t13*(t4*F/t2+t8*F/t6)*t11
      FVV = t25
      RETURN
      END


      subroutine pstsa2_fv2(C,F,V,FV,FVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------------------
! This subroutine is used to calculate the f(v) in the a2(r)
! the expression is the second one is the B2 (appendix B)      
! It calculates two terms, on is f(v), the other is derivatives of
! f(v) with respect to the v      
!--------------------------------------------------------------------
      REAL*8 C,F,V,FV,FVV

      t3 = (v-1.D0*C)**F
      t7 = (1.D0-1.D0*v)**F
      t11 = dexp(1/t3-1.D0/t7)
      t14 = t11/(1.D0+t11)
      FV  = t14
      t2 = v-1.D0*C
      t3 = t2**F
      t4 = 1/t3
      t10 = 1.D0-1.D0*v
      t11 = t10**F
      t12 = 1/t11
      t17 = -1.D0*t4*F/t2-1.D0*t12*F/t10
      t20 = dexp(t4-1.D0*t12)
      t22 = 1.D0+t20
      t25 = t20**2
      t26 = t22**2
      t31 = t17*t20/t22-1.D0*t25/t26*t17
      fvv = t31
      RETURN
      END

      subroutine functional_psts
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DrhoA,DrhoB,TauA,TauB,
     &EXCHA,EXCHB,F,D1F)
      IMPLICIT NONE
!--------------------------------------------------------------------
! INPUT :
! NG         : the number of grid points
! rhoA       : the alpha electron density
! rhoB       : the beta  electron density
! DRhoA      : the alpha rho'
! DRhoB      : the beta  rho'
! TauA       : the alpha kinetic energy density
! TauB       : the beta  kinetic energy density
! OUTPUT:
! F          : functional values
! D1F        : the first order functional derivatives
!--------------------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 

      ! these are the psts functional position variable
      INTEGER  INFOR(*)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS,             ID_TB_POS
      INTEGER  ID_EXA_POS,ID_EXB_POS

      ! for TPSS functional, it's META-GGA and we do not have
      ! exchange energy density
      ! so we use different variables to represent the position
      ! information here
      INTEGER  VAR_INFOR_TPSS(MAX_VAR_TYPE)
      INTEGER  D1VARS_TPSS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS_TPSS, ID_GAA_POS_TPSS, ID_GAB_POS_TPSS
      INTEGER  ID_RB_POS_TPSS, ID_GBB_POS_TPSS 
      INTEGER  ID_TA_POS_TPSS, ID_TB_POS_TPSS 

      ! input and output
      INTEGER NG
      INTEGER NDEN ! number of density
      REAL*8 rhoA(NG),rhoB(NG)
      REAL*8 DRhoA(NG,3),DRhoB(NG,3)
      REAL*8 TauA(NG),TauB(NG)
      REAL*8 EXCHA(NG),EXCHB(NG)
      REAL*8 F(NG), D1F(NG,*)

      ! ordinary variables
      REAL*8 RA,RB
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 TA,TB

      ! exchange energy density variable for given grid point
      REAL*8 exa,exb,ex

      ! pointers for functional and its derivatives
      REAL*8 TF
      REAL*8 TD1F(N_FUNC_DERIV_1)

      ! functional value and 1st derivatives for tpss exchange
      REAL*8 t,tra,trb,tgaa,tgab,tgbb,tta,ttb

      ! functional value and 1st derivatives for a1(r) 
      REAL*8 a1,a1ra,a1rb,a1gaa,a1gab,a1gbb,a1ta,a1tb

      ! functional value and 1st derivatives for a2(r)
      ! a2exa and a2exb is the derivatives for a2 with respect to the 
      ! alpha or beta exchange variables
      REAL*8 a2,a2ra,a2rb,a2gaa,a2gab,a2gbb,a2ta,a2tb,a2exa,a2exb

      ! functional value and 1st derivatives for a(r)
      ! aexa and aexb is the derivatives with respect to the 
      ! alpha or beta exchange variables
      REAL*8 a,ara,arb,agaa,agab,agbb,ata,atb,aexa,aexb

      ! variables and constants
      REAL*8 F1
      REAL*8 F12
      PARAMETER(F1  =  1.0D0)
      PARAMETER(F12 =  0.5D0)
      REAL*8 TOL
      INTEGER I

      ! init position variables
      ! this is for PSTS functional
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_TB_POS  = D1VARS(ID_TB)
      ID_EXA_POS = D1VARS(ID_EXA)
      ID_EXB_POS = D1VARS(ID_EXB)

      ! initilize variable position information for META-GGA functional
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR_TPSS(I) = -1
      END DO
      VAR_INFOR_TPSS(ID_RHO)   = 1
      VAR_INFOR_TPSS(ID_GAMMA) = 1
      VAR_INFOR_TPSS(ID_TAU)   = 1
      CALL INIT_FUNC_DERIV_1(VAR_INFOR_TPSS,D1VARS_TPSS)
      ID_RA_POS_TPSS  = D1VARS_TPSS(ID_RA)
      ID_RB_POS_TPSS  = D1VARS_TPSS(ID_RB)
      ID_GAA_POS_TPSS = D1VARS_TPSS(ID_GAA)
      ID_GAB_POS_TPSS = D1VARS_TPSS(ID_GAB)
      ID_GBB_POS_TPSS = D1VARS_TPSS(ID_GBB)
      ID_TA_POS_TPSS  = D1VARS_TPSS(ID_TA)
      ID_TB_POS_TPSS  = D1VARS_TPSS(ID_TB)

      !----------------------------------------------------------------
      ! real work begins for both restricted and un-restricted cases
      !----------------------------------------------------------------
      DO I = 1, NG

         ! input variable for each grid point
         RA = rhoA(i)
         DRhoAX = DRhoA(i,1)
         DRhoAY = DRhoA(i,2)
         DRhoAZ = DRhoA(i,3)
         TA     = TauA(i)
         exa    = F12*EXCHA(i)

         ! for debugging purpose
C         IF (I .eq. 1) THEN
C            exa    = -327.269085669823D0;
C         ELSE IF (I .eq. 2) THEN
C            exa    = -327.269085703956D0;
C         ELSE IF (I .eq. 3) THEN
C            exa    = -327.269085686892D0;
C         ELSE IF (I .eq. 4) THEN
C            exa    = -327.269085686890D0;
C         END IF

         IF (NDEN .EQ. 1) THEN
            RB     = RA 
            DRhoBX = DRhoAX
            DRhoBY = DRhoAY
            DRhoBZ = DRhoAZ
            TB     = TA
            exb    = exa
         ELSE
            RB     = rhoB(i) 
            DRhoBX = DRhoB(i,1)
            DRhoBY = DRhoB(i,2)
            DRhoBZ = DRhoB(i,3)
            TB     = TauB(i)
            exb    = F12*EXCHB(i)
         END IF
         ex  = exa + exb

C         write(6,*)"----------------------------"
C         write(6,*)"for grid ", i
C         write(6,*)"----------------------------"
C         write(6,*)"alpha rho", RA
C         write(6,*)"alpha DRhoAX", DRhoAX
C         write(6,*)"alpha DRhoAY", DRhoAY
C         write(6,*)"alpha DRhoAZ", DRhoAZ
C         write(6,*)"alpha tau", TA
C         write(6,*)"alpha EXA", exa

         ! consider the small density separately for alpha and
         ! beta... 
         IF (RA .LT. TOL .AND. RB .LT. TOL) CYCLE

         ! get the a1(r)
         TF  = 0.0D0
         CALL pstsa1(NDEN,TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,TA,TB,TF,TD1F)
         a1    = TF
         a1ra  = TD1F(ID_RA_POS_TPSS)
         a1rb  = TD1F(ID_RB_POS_TPSS)
         a1gaa = TD1F(ID_GAA_POS_TPSS)
         a1gab = TD1F(ID_GAB_POS_TPSS)
         a1gbb = TD1F(ID_GBB_POS_TPSS)
         a1ta  = TD1F(ID_TA_POS_TPSS)
         a1tb  = TD1F(ID_TB_POS_TPSS)

         ! get the a2(r)
         ! see the functional definition above
         ! we already state that a2(r) is like META-GGA functional
         ! so there's no deriv in terms of exa and exb in TD1F
         TF    = 0.0D0
         a2exa = 0.0D0
         a2exb = 0.0D0
         CALL pstsa2(NDEN,TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,TA,TB,exa,exb,
     &        a2exa,a2exb,TF,TD1F)
         a2    = TF
         a2ra  = TD1F(ID_RA_POS_TPSS)
         a2rb  = TD1F(ID_RB_POS_TPSS)
         a2gaa = TD1F(ID_GAA_POS_TPSS)
         a2gab = TD1F(ID_GAB_POS_TPSS)
         a2gbb = TD1F(ID_GBB_POS_TPSS)
         a2ta  = TD1F(ID_TA_POS_TPSS)
         a2tb  = TD1F(ID_TB_POS_TPSS)

         ! finally, calculate the 1-a(r), and its derivatives
         a    = (F1-a1)*(F1-a2) 
         ara  = -a1ra *(F1-a2) -(F1-a1)*a2ra 
         agaa = -a1gaa*(F1-a2) -(F1-a1)*a2gaa 
         agab = -a1gab*(F1-a2) -(F1-a1)*a2gab 
         ata  = -a1ta *(F1-a2) -(F1-a1)*a2ta 
         aexa = -(F1-a1)*a2exa 

         ! consider the beta part 
         IF (NDEN .EQ. 1) THEN
            arb  = ara
            agbb = agaa
            atb  = ata
            aexb = aexa 
         ELSE
            arb  = -a1rb *(F1-a2) -(F1-a1)*a2rb 
            agbb = -a1gbb*(F1-a2) -(F1-a1)*a2gbb 
            atb  = -a1tb *(F1-a2) -(F1-a1)*a2tb 
            aexb = -(F1-a1)*a2exb 
         END IF

         ! get the functional and its derivatives for tpss exchange
         TF = 0.0D0
         IF (NDEN .EQ. 1) THEN
            CALL functional_tpssx_pw_close(TOL,RA,DRhoAX,DRhoAY,
     &           DRhoAZ,TA,TF,TD1F)
         ELSE
            CALL functional_tpssx_pw_open(TOL,RA,RB,DRhoAX,DRhoAY,
     &           DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,TA,TB,TF,TD1F)
         END IF
         t    = TF
         tra  = TD1F(ID_RA_POS_TPSS)
         trb  = TD1F(ID_RB_POS_TPSS)
         tgaa = TD1F(ID_GAA_POS_TPSS)
         tgab = TD1F(ID_GAB_POS_TPSS)
         tgbb = TD1F(ID_GBB_POS_TPSS)
         tta  = TD1F(ID_TA_POS_TPSS)
         ttb  = TD1F(ID_TB_POS_TPSS)

         ! finally let us calculate the results
         F(i) = F(i) + a*(t-ex) 
         D1F(i,ID_RA_POS)  = D1F(i,ID_RA_POS)  
     &                  + ara*(t-ex)  + a*tra 
         D1F(i,ID_GAA_POS) = D1F(i,ID_GAA_POS)  
     &                  + agaa*(t-ex) + a*tgaa 
         D1F(i,ID_GAB_POS) = D1F(i,ID_GAB_POS) 
     &                  + agab*(t-ex) + a*tgab
         D1F(i,ID_TA_POS)  = D1F(i,ID_TA_POS) 
     &                  + ata*(t-ex)  + a*tta 
         D1F(i,ID_EXA_POS) = D1F(i,ID_EXA_POS)  
     &                  + aexa*(t-ex) - a 

         ! beta part
         IF (NDEN .EQ. 1) THEN
            D1F(i,ID_RB_POS)  = D1F(i,ID_RA_POS)
            D1F(i,ID_GBB_POS) = D1F(i,ID_GAA_POS)
            D1F(i,ID_TB_POS)  = D1F(i,ID_TA_POS)
            D1F(i,ID_EXB_POS) = D1F(i,ID_EXA_POS)  
         ELSE
            D1F(i,ID_RB_POS)  = D1F(i,ID_RB_POS)  
     &                  + arb*(t-ex)  + a*trb 
            D1F(i,ID_GBB_POS) = D1F(i,ID_GBB_POS)  
     &                  + agbb*(t-ex) + a*tgbb 
            D1F(i,ID_TB_POS)  = D1F(i,ID_TB_POS) 
     &                  + atb*(t-ex)  + a*ttb 
            D1F(i,ID_EXB_POS) = D1F(i,ID_EXB_POS)  
     &                  + aexb*(t-ex) - a 
         END IF

C         write(6,*)"derivatives for alpha rho", ara*(t-ex)  + a*tra
C         write(6,*)"derivatives for GAA", agaa*(t-ex) + a*tgaa
C         write(6,*)"derivatives for GAB", agab*(t-ex) + a*tgab
C         write(6,*)"derivatives for GBB", agbb*(t-ex) + a*tgbb
C         write(6,*)"derivatives for alpha tau", ata*(t-ex)  + a*tta
C         write(6,*)"derivatives for alpha EXA", aexa*(t-ex) - a

      END DO
      END



