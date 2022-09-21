!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! Here it's the TPSS correlation functional, see      
! Climbing the density functional ladder: Nonempirical 
! meta-generalized gradient approximation designed for 
! molecules and solids
! J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
! Phys. Rev. Lett., 91 (2003) 146401
!      
! The lambda-dependent TPSS functional is shown in the 
! MCY functional, see:
! Mori-Sanchez, P. and Cohen, A.J. and Yang, W.
! Self-interaction-free exchange-correlation functional 
! for thermochemistry and kinetics
! The Journal of chemical physics, 2006,124, 091102  
!
! This is the Gorling-Levy limit of TPSS correlation energy per
! electrion
! see the paper of PSTS functional, the Appendix A:
! Density functional with full exact exchange, balanced nonlocality 
! of correlation, and constraint satisfaction
! Perdew, J.P. and Staroverov, V.N. and Tao, J. and Scuseria, G.E.
!
! note: 
!
! we note that TPSS correlation functional uses the "small tau"(
! tau variable is with 1/2 factor), however, in this program we always
! use the "big tau"(tau variable is without 1/2). Such transformation
! is already been handled inside. For TPSS correlation, the tau
! is used in the function of wt.f(please see the code over there),
! over there in using the tau variable we multiply 1/2.
!
!

      subroutine tpssc_gllimit(NDEN,TOL,RA,RB,DrhoAX,DrhoAY,DrhoAZ,
     &DrhoBX,DrhoBY,DrhoBZ,TA,TB,F,D1F)
      IMPLICIT NONE
!--------------------------------------------------------------------
! This is the Gorling-Levy limit of TPSS correlation energy per
! electrion
! It's calculated on a given grid point      
! INPUT :
! RA         : the alpha electron density
! RB         : the beta  electron density
! DRhoA      : the alpha rho'
! DRhoB      : the beta  rho'
! TA         : the alpha kinetic energy density
! TB         : the beta  kinetic energy density
! OUTPUT:
! F     : functional values
! D1F   : the first  order functional derivatives
!--------------------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 
      INTEGER  I
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS

      ! input and output
      INTEGER NDEN ! number of density
      REAL*8 RA,RB
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 TA,TB
      REAL*8 F, D1F(*)

      ! WT functional and REVPKZB functional
      REAL*8 F_pkzb,F_wt
      REAL*8 D1F_pkzb(N_FUNC_DERIV_1),D1F_wt(N_FUNC_DERIV_1)

      ! functional value and 1st derivatives for WT
      REAL*8 w,wra,wrb,wgaa,wgab,wgbb,wta,wtb

      ! functional value and 1st derivatives for REVPKZB
      REAL*8 p,pra,prb,pgaa,pgab,pgbb,pta,ptb

      ! variables and constants
      REAL*8 d,F2,F3
      REAL*8 TOL
      INTEGER IS_GLLIMIT,USE_LAMBDA

      !----------------------------------------------------------------
      ! preparation
      !----------------------------------------------------------------

      ! constants
      d  = 2.8D0
      F2 = 2.0D0
      F3 = 3.0D0

      ! in the PKZB functional calculation, we should calculate
      ! the GL LIMIT case, not lambda-dependent functional
      IS_GLLIMIT = 1
      USE_LAMBDA = 0

      ! initilize variable position information
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR(I) = -1
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
      ! real work begins, we do the restricted and un-restricted 
      ! together
      !----------------------------------------------------------------

      ! consider the small density
      ! here we do not specify that RA > TOL or RA < TOL ...
      ! cases, since they are specified in the PKZB and wt
      ! functionals. So we only care about the small density
      ! situation
      IF (RA .LT. TOL .AND. RB .LT. TOL) go to 10

      ! get the WT and PKZB functional infor
      F_wt   = 0.0D0
      F_pkzb = 0.0D0
      IF (NDEN .EQ. 1) THEN
         CALL functional_wt_close(TOL,RA,DRhoAX,DRhoAY,DRhoAZ,
     &        TA,F_wt,D1F_wt)
      ELSE
         CALL functional_wt_open(TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,TA,TB,F_wt,D1F_wt)
      END IF
      CALL functional_revpkzb(NDEN,TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,F_wt,D1F_wt,F_pkzb,D1F_pkzb,
     &        USE_LAMBDA,IS_GLLIMIT)
      w    = F_wt
      wra  = D1F_wt(ID_RA_POS)
      wrb  = D1F_wt(ID_RB_POS)
      wgaa = D1F_wt(ID_GAA_POS)
      wgab = D1F_wt(ID_GAB_POS)
      wgbb = D1F_wt(ID_GBB_POS)
      wta  = D1F_wt(ID_TA_POS)
      wtb  = D1F_wt(ID_TB_POS)
      p    = F_pkzb
      pra  = D1F_pkzb(ID_RA_POS)
      prb  = D1F_pkzb(ID_RB_POS)
      pgaa = D1F_pkzb(ID_GAA_POS)
      pgab = D1F_pkzb(ID_GAB_POS)
      pgbb = D1F_pkzb(ID_GBB_POS)
      pta  = D1F_pkzb(ID_TA_POS)
      ptb  = D1F_pkzb(ID_TB_POS)
C      write(6,*)"ourtau is:",w
C      write(6,*)"ourtaura is:",wra
C      write(6,*)"ourtaurb is:",wrb
C      write(6,*)"ourtaugaa is:",wgaa
C      write(6,*)"ourtaugab is:",wgab
C      write(6,*)"ourtaugbb is:",wgbb
C      write(6,*)"ourtauta is:",wta
C      write(6,*)"ourtautb is:",wtb
C       write(6,*)"ourpkzb is:",p
C      write(6,*)"ourpkzbra is:",pra
C      write(6,*)"ourpkzbrb is:",prb
C      write(6,*)"ourpkzbgaa is:",pgaa
C      write(6,*)"ourpkzbgbb is:",pgbb
C      write(6,*)"ourpkzbgab is:",pgab
C      write(6,*)"ourpkzbta is:",pta
C      write(6,*)"ourpkzbtb is:",ptb

      ! real calculation....
      F = p + d*p*p*w*w*w
      D1F(ID_RA_POS)  = D1F(ID_RA_POS)  
     &             + pra 
     &             + F2*d*p*pra*w*w*w 
     &             + F3*d*p*p*w*w*wra       
      D1F(ID_GAA_POS) = D1F(ID_GAA_POS)  
     &             + pgaa  
     &             + F2*d*p*pgaa*w*w*w 
     &             + F3*d*p*p*w*w*wgaa      
      D1F(ID_GAB_POS) = D1F(ID_GAB_POS) 
     &             + pgab  
     &             + F2*d*p*pgab*w*w*w 
     &             + F3*d*p*p*w*w*wgab      
      D1F(ID_TA_POS)  = D1F(ID_TA_POS) 
     &             + pta  
     &             + F2*d*p*pta*w*w*w 
     &             + F3*d*p*p*w*w*wta      


      ! beta is same with alpha for 1st derivatives
      IF (NDEN .EQ. 1) THEN
         D1F(ID_RB_POS)  = D1F(ID_RA_POS)
         D1F(ID_GBB_POS) = D1F(ID_GAA_POS)
         D1F(ID_TB_POS)  = D1F(ID_TA_POS)
      ELSE
         D1F(ID_RB_POS)  = D1F(ID_RB_POS)  
     &                + prb 
     &                + F2*d*p*prb*w*w*w 
     &                + F3*d*p*p*w*w*wrb       
         D1F(ID_GBB_POS) = D1F(ID_GBB_POS)  
     &                + pgbb  
     &                + F2*d*p*pgbb*w*w*w 
     &                + F3*d*p*p*w*w*wgbb      
         D1F(ID_TB_POS)  = D1F(ID_TB_POS) 
     &                + ptb  
     &                + F2*d*p*ptb*w*w*w 
     &                + F3*d*p*p*w*w*wtb      
      END IF

C        write(6,*)"ourtpss is:",F(i)
C        write(6,*)"ourtpssRA is:",D1F(i,ID_RA_POS)
C        write(6,*)"ourtpssRB is:",D1F(i,ID_RB_POS)
C        write(6,*)"ourtpssGAA is:",D1F(i,ID_GAA_POS)
C        write(6,*)"ourtpssGAB is:",D1F(i,ID_GAB_POS)
C        write(6,*)"ourtpssGBB is:",D1F(i,ID_GBB_POS)
C        write(6,*)"ourtpssTA  is:",D1F(i,ID_TA_POS)
C        write(6,*)"ourtpssTB  is:",D1F(i,ID_TB_POS)

 10   continue     
      END


      subroutine functional_tpssc
     &(INFOR,NG,NDEN,TOL,USE_LAMBDA,rhoA,rhoB,DrhoA,DrhoB,
     &TauA,TauB,F,D1F)
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
! use_lambda : whether it's the lambda-dependent functional  
!              if the value is 1 then it's lambda-TPSS correlation       
! OUTPUT:
! F     : functional values
! D1F   : the first  order functional derivatives
!--------------------------------------------------------------------
#include "fderiv1.inc"
      INTEGER INFOR(*)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS

      ! input and output
      INTEGER NG
      INTEGER USE_LAMBDA
      INTEGER NDEN ! number of density
      REAL*8 rhoA(NG),rhoB(NG)
      REAL*8 DRhoA(NG,3),DRhoB(NG,3)
      REAL*8 TauA(NG),TauB(NG)
      REAL*8 F(NG), D1F(NG,*)

      ! pointers for WT functional and REVPKZB functional
      REAL*8 F_pkzb,F_wt
      REAL*8 D1F_pkzb(N_FUNC_DERIV_1),D1F_wt(N_FUNC_DERIV_1)

      ! functional value and 1st derivatives for WT
      REAL*8 w,wra,wrb,wgaa,wgab,wgbb,wta,wtb

      ! functional value and 1st derivatives for REVPKZB
      REAL*8 p,pra,prb,pgaa,pgab,pgbb,pta,ptb

      ! variables and constants
      REAL*8 RA,RB,RT
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 TA,TB
      REAL*8 d,F2,F3,Fac
      REAL*8 TOL
      INTEGER I

      ! GL choice
      INTEGER IS_GLLIMIT


      !----------------------------------------------------------------
      ! preparation
      !----------------------------------------------------------------

      ! constants
      d  = 2.8D0
      F2 = 2.0D0
      F3 = 3.0D0
      IS_GLLIMIT = 0

      ! we have to remember that in lambda-TPSS, it is 2*Ec
      Fac = 1.0D0
      IF (USE_LAMBDA .EQ. 1) THEN
         Fac = 2.0D0
      END IF

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_TB_POS  = D1VARS(ID_TB)
      
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
         IF (NDEN .EQ. 1) THEN
            RB     = RA 
            DRhoBX = DRhoAX
            DRhoBY = DRhoAY
            DRhoBZ = DRhoAZ
            TB     = TA
         ELSE
            RB     = rhoB(i) 
            DRhoBX = DRhoB(i,1)
            DRhoBY = DRhoB(i,2)
            DRhoBZ = DRhoB(i,3)
            TB     = TauB(i)
         END IF

         ! consider the small density separately for alpha and
         ! beta... we are thinking that when the functional
         ! value should be zero 
         IF (RA .LT. TOL .AND. RB .LT. TOL) CYCLE

         ! get the WT and PKZB functional infor
         F_wt   = 0.0D0
         F_pkzb = 0.0D0
         IF (NDEN .EQ. 1) THEN
            CALL functional_wt_close(TOL,RA,DRhoAX,DRhoAY,DRhoAZ,
     &           TA,F_wt,D1F_wt)
         ELSE
            CALL functional_wt_open(TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &           DRhoBX,DRhoBY,DRhoBZ,TA,TB,F_wt,D1F_wt)
         END IF
         CALL functional_revpkzb(NDEN,TOL,RA,RB,DRhoAX,DRhoAY,DRhoAZ,
     &        DRhoBX,DRhoBY,DRhoBZ,F_wt,D1F_wt,F_pkzb,D1F_pkzb,
     &        USE_LAMBDA,IS_GLLIMIT)
         w    = F_wt
         wra  = D1F_wt(ID_RA_POS)
         wrb  = D1F_wt(ID_RB_POS)
         wgaa = D1F_wt(ID_GAA_POS)
         wgab = D1F_wt(ID_GAB_POS)
         wgbb = D1F_wt(ID_GBB_POS)
         wta  = D1F_wt(ID_TA_POS)
         wtb  = D1F_wt(ID_TB_POS)
         p    = F_pkzb
         pra  = D1F_pkzb(ID_RA_POS)
         prb  = D1F_pkzb(ID_RB_POS)
         pgaa = D1F_pkzb(ID_GAA_POS)
         pgab = D1F_pkzb(ID_GAB_POS)
         pgbb = D1F_pkzb(ID_GBB_POS)
         pta  = D1F_pkzb(ID_TA_POS)
         ptb  = D1F_pkzb(ID_TB_POS)
C         write(6,*)"ourtau is:",w
C         write(6,*)"ourtaura is:",wra
C         write(6,*)"ourtaurb is:",wrb
C         write(6,*)"ourtaugaa is:",wgaa
C         write(6,*)"ourtaugab is:",wgab
C         write(6,*)"ourtaugbb is:",wgbb
C         write(6,*)"ourtauta is:",wta
C         write(6,*)"ourtautb is:",wtb
C         write(6,*)"ourpkzb is:",p
C         write(6,*)"ourpkzbra is:",pra
C         write(6,*)"ourpkzbrb is:",prb
C         write(6,*)"ourpkzbgaa is:",pgaa
C         write(6,*)"ourpkzbgbb is:",pgbb
C         write(6,*)"ourpkzbgab is:",pgab
C         write(6,*)"ourpkzbta is:",pta
C         write(6,*)"ourpkzbtb is:",ptb

         ! real calculation....
         RT   = RA+RB
         F(i) = F(i)    + Fac*(RT*p     + d*RT*p*p*w*w*w)
         D1F(i,ID_RA_POS)  = D1F(i,ID_RA_POS)  
     &                  + Fac*(p        +  RT*pra   + d*p*p*w*w*w 
     &                  + F2*d*RT*p*pra*w*w*w +
     &                    F3*d*RT*p*p*w*w*wra)
         D1F(i,ID_GAA_POS) = D1F(i,ID_GAA_POS)  
     &                  + Fac*(RT*pgaa  
     &                  + F2*d*RT*p*pgaa*w*w*w 
     &                  + F3*d*RT*p*p*w*w*wgaa)      
         D1F(i,ID_GAB_POS) = D1F(i,ID_GAB_POS) 
     &                  + Fac*(RT*pgab  
     &                  + F2*d*RT*p*pgab*w*w*w 
     &                  + F3*d*RT*p*p*w*w*wgab)      
         D1F(i,ID_TA_POS)  = D1F(i,ID_TA_POS) 
     &                  + Fac*(RT*pta  
     &                  + F2*d*RT*p*pta*w*w*w 
     &                  + F3*d*RT*p*p*w*w*wta)      


         IF (NDEN .EQ. 1) THEN
            D1F(i,ID_RB_POS)  = D1F(i,ID_RA_POS)
            D1F(i,ID_GBB_POS) = D1F(i,ID_GAA_POS)
            D1F(i,ID_TB_POS)  = D1F(i,ID_TA_POS)
         ELSE
            D1F(i,ID_RB_POS)  = D1F(i,ID_RB_POS)  
     &                     + Fac*(p        +  RT*prb   + d*p*p*w*w*w 
     &                     + F2*d*RT*p*prb*w*w*w 
     &                     + F3*d*RT*p*p*w*w*wrb)       
            D1F(i,ID_GBB_POS) = D1F(i,ID_GBB_POS)  
     &                     + Fac*(RT*pgbb  
     &                     + F2*d*RT*p*pgbb*w*w*w 
     &                     + F3*d*RT*p*p*w*w*wgbb)      
            D1F(i,ID_TB_POS)  = D1F(i,ID_TB_POS) 
     &                     + Fac*(RT*ptb  
     &                     + F2*d*RT*p*ptb*w*w*w 
     &                     + F3*d*RT*p*p*w*w*wtb)      
         END IF



C        write(6,*)"ourtpss is:",F(i)
C        write(6,*)"ourtpssRA is:",D1F(i,ID_RA_POS)
C        write(6,*)"ourtpssRB is:",D1F(i,ID_RB_POS)
C        write(6,*)"ourtpssGAA is:",D1F(i,ID_GAA_POS)
C        write(6,*)"ourtpssGAB is:",D1F(i,ID_GAB_POS)
C        write(6,*)"ourtpssGBB is:",D1F(i,ID_GBB_POS)
C        write(6,*)"ourtpssTA  is:",D1F(i,ID_TA_POS)
C        write(6,*)"ourtpssTB  is:",D1F(i,ID_TB_POS)

      END DO

      END


      subroutine functional_revpkzb
     &(NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,
     & F_wt,D1F_wt,F,D1F,USE_LAMBDA,IS_GLLIMIT)
      IMPLICIT NONE
!--------------------------------------------------------------------
! This is the revised PKZB functional used in TPSS
! WT functional is tau^w/tau, which is the pieces inside the PKZB(
! see their original paper for more information)      
! INPUT     :
! rhoA      :  alpha density for a given grid point
! rhoB      :  beta  density for a given grid point
! DRhoAX    :  gradient density for alpha on X direction      
! DRhoAY    :  gradient density for alpha on Y direction      
! DRhoAZ    :  gradient density for alpha on Z direction      
! DRhoBX    :  gradient density for beta  on X direction      
! DRhoBY    :  gradient density for beta  on Y direction      
! DRhoBZ    :  gradient density for beta  on Z direction  
! F_wt      :  functional value for WT functional    
! D1F_wt    :  1st functional derivatives array for given grid point 
! USE_LAMBDA:  whether to use lambda-dependent functional  
! IS_GLLIMiT:  Whether it's Gorling-Levy limit form      
! OUTPUT:
! F     :  functional value for given grid point
! D1F   :  1st functional derivatives array for given grid point           
!--------------------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 

      ! this is for the PBE functional, store it's variable information
      INTEGER  VAR_PBE_INFOR(MAX_VAR_TYPE)
      INTEGER  D1VARS_PBE(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS_PBE, ID_GAA_POS_PBE, ID_GAB_POS_PBE
      INTEGER  ID_RB_POS_PBE, ID_GBB_POS_PBE 

      ! this is for the all of other functionals, META-GGA functionals
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS,             ID_TB_POS

      ! input and output 
      REAL*8 rhoA,rhoB
      REAL*8 DRhoAX,DRhoAY,DRhoAZ
      REAL*8 DRhoBX,DRhoBY,DRhoBZ
      REAL*8 F_wt,D1F_wt(*)
      REAL*8 F,D1F(*)
      INTEGER USE_LAMBDA
      INTEGER IS_GLLIMIT
      INTEGER NDEN

      ! constants 
      INTEGER  I
      REAL*8 F12,F1,F2
      REAL*8 TOL

      ! variables
      REAL*8 RA,RB,RhoT,xia,xib,xia_a,xia_b,xib_a,xib_b

      ! PBE functional and its 1st derivatives
      ! s indicates it is the spin-dependent
      REAL*8 F_pbes
      REAL*8 F_pbe
      REAL*8 D1F_pbe(N_FUNC_DERIV_1),D1F_pbes(N_FUNC_DERIV_1)
      REAL*8 p,pra,prb,pgaa,pgbb,pgab

      ! C(zeta,xi) functional and its 1st derivatives
      REAL*8 F_czeta,D1F_czeta(N_FUNC_DERIV_1)
      REAL*8 c,cra,crb,cgaa,cgab,cgbb

      ! WT functional and its 1st derivatives
      real*8 w,wra,wrb,wgaa,wgab,wgbb,wta,wtb

      !----------------------------------------------------------------
      ! preparation work
      !----------------------------------------------------------------

      ! initilize variable position information for PBE functional
      DO I = 1, MAX_VAR_TYPE
         VAR_PBE_INFOR(I) = -1
      END DO
      VAR_PBE_INFOR(ID_RHO)   = 1
      VAR_PBE_INFOR(ID_GAMMA) = 1
      CALL INIT_FUNC_DERIV_1(VAR_PBE_INFOR,D1VARS_PBE)
      ID_RA_POS_PBE  = D1VARS_PBE(ID_RA)
      ID_RB_POS_PBE  = D1VARS_PBE(ID_RB)
      ID_GAA_POS_PBE = D1VARS_PBE(ID_GAA)
      ID_GAB_POS_PBE = D1VARS_PBE(ID_GAB)
      ID_GBB_POS_PBE = D1VARS_PBE(ID_GBB)

      ! initilize variable position information for META-GGA functional
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

      ! arrange output for the WT functional
      w    = F_wt
      wra  = D1F_wt(ID_RA_POS)
      wrb  = D1F_wt(ID_RB_POS)
      wgaa = D1F_wt(ID_GAA_POS)
      wgab = D1F_wt(ID_GAB_POS)
      wgbb = D1F_wt(ID_GBB_POS)
      wta  = D1F_wt(ID_TA_POS)
      wtb  = D1F_wt(ID_TB_POS)

      ! arrange output for the C(zeta,xi) functional
      ! we note that this functional is same with PBE
      ! also the GGA functional, so use the same arragement of variable
      ! position
C      IF (USE_LAMBDA .EQ. 1 .OR. IS_GLLIMIT .EQ. 1) THEN
      IF (USE_LAMBDA .EQ. 1 ) THEN

         ! here we make to make additional explanation for the C(zeta,xi)
         ! function in the lambda-dependent case.
         ! The C can be 0.53 or 0, which is actually depends on the xi
         ! because C is the explicit function of xi. If xi is zero, even
         ! here it may be spin-unpolarized case, then the c is 0.53 rather
         ! than 0. Here 
         ! we note that for the lambda-dependent case, 
         ! because of the lambda all of its derivatives should be 
         ! very small, so they are zero
         IF (dabs((rhoA-rhoB)/(rhoA+rhoB)) .LT. 1.0D-12) THEN
            c    =  0.53D0
            cra  =  0.0D0
            crb  =  0.0D0
            cgaa =  0.0D0
            cgab =  0.0D0
            cgbb =  0.0D0
         ELSE
            c    =  0.0D0
            cra  =  0.0D0
            crb  =  0.0D0
            cgaa =  0.0D0
            cgab =  0.0D0
            cgbb =  0.0D0
         END IF

      ELSE

         F_czeta   = 0.0D0
         IF (NDEN .EQ. 1) THEN
            call functional_czeta_close(TOL,rhoA,DRhoAX,DRhoAY, 
     &            DRhoAZ, F_czeta, D1F_czeta)
            c    =  F_czeta
            cra  =  D1F_czeta(ID_RA_POS_PBE)
            crb  =  D1F_czeta(ID_RB_POS_PBE)
            cgaa =  D1F_czeta(ID_GAA_POS_PBE)
            cgab =  D1F_czeta(ID_GAB_POS_PBE)
            cgbb =  D1F_czeta(ID_GBB_POS_PBE)
         ELSE
            call functional_czeta_open(TOL,rhoA,rhoB,DRhoAX,DRhoAY, 
     &            DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,F_czeta,D1F_czeta)
            c    =  F_czeta
            cra  =  D1F_czeta(ID_RA_POS_PBE)
            crb  =  D1F_czeta(ID_RB_POS_PBE)
            cgaa =  D1F_czeta(ID_GAA_POS_PBE)
            cgab =  D1F_czeta(ID_GAB_POS_PBE)
            cgbb =  D1F_czeta(ID_GBB_POS_PBE)
         END IF

      END IF
C      write(6,*)"ourpkzbc is:",c
C      write(6,*)"cra is:",cra
C      write(6,*)"crb is:",crb
C      write(6,*)"cgaa is:",cgaa
C      write(6,*)"cgab is:",cgab
C      write(6,*)"cgbb is:",cgbb


      !constants 
      F12   = 0.5D0
      F1    = 1.0D0
      F2    = 2.0D0

      !----------------------------------------------------------------
      ! real work begins.....
      !----------------------------------------------------------------

      ! restricted case......
      IF (NDEN .EQ. 1) THEN

         ! here we consider the small density
         IF (rhoA .GT. TOL) THEN

            !-------------------------------------------------------
            ! deal with alpha + beta part
            !-------------------------------------------------------
            ! get the PBE value
            F_pbe = 0.0D0
            IF (USE_LAMBDA .EQ. 1) THEN
               call functional_lambda_pbe_close(TOL,rhoA,DRhoAX,DRhoAY,
     &            DRhoAZ, F_pbe, D1F_pbe)
            ELSE IF (IS_GLLIMIT .EQ. 1) THEN
               call functional_glpbe_pw_close(TOL,rhoA,DRhoAX,DRhoAY,
     &            DRhoAZ, F_pbe, D1F_pbe)
            ELSE
               call functional_pbecp_close(TOL,rhoA,DRhoAX,DRhoAY,
     &            DRhoAZ, F_pbe, D1F_pbe)
            END IF
            p    =  F_pbe
            pra  =  D1F_pbe(ID_RA_POS_PBE)
            prb  =  D1F_pbe(ID_RB_POS_PBE)
            pgaa =  D1F_pbe(ID_GAA_POS_PBE)
            pgab =  D1F_pbe(ID_GAB_POS_PBE)
            pgbb =  D1F_pbe(ID_GBB_POS_PBE)
C           write(6,*)"ourpbe is:",p
C           write(6,*)"ourpkzbpbe  is:",p
C           write(6,*)"ourpbera is:",pra
C           write(6,*)"ourpberb is:",prb
C           write(6,*)"ourpbegaa is:",pgaa
C           write(6,*)"ourpbegab is:",pgab
C           write(6,*)"ourpbegbb is:",pgbb

            ! this is  (1 + c*wt^2)*F_pbe
            F            = F            +    p*(F1 + c*w*w)
            D1F(ID_RA_POS)  = D1F(ID_RA_POS)  + 
     &                     pra*(F1+c*w*w)  + p*cra*w*w  + F2*p*c*w*wra
            D1F(ID_GAA_POS) = D1F(ID_GAA_POS) + 
     &                     pgaa*(F1+c*w*w) + p*cgaa*w*w + F2*p*c*w*wgaa
            D1F(ID_GAB_POS) = D1F(ID_GAB_POS) + 
     &                     pgab*(F1+c*w*w) + p*cgab*w*w + F2*p*c*w*wgab
            D1F(ID_TA_POS)  = D1F(ID_TA_POS)  +  F2*p*c*w*wta 


            !-------------------------------------------------------
            !deal with spin polarized part, with constraint of RA = RB
            !-------------------------------------------------------
            !this is -(1+c)*wt^2*sum(srho_a*pbea + srho_b*pbeb)
            !since RA = RB etc., so srho_a = srho_b = 1/2
            !pbea = pbeb
            !Hence finally the expression is:
            !-(1+c)*wt^2*pbea 

            F_pbes   = 0.0D0
            IF (USE_LAMBDA .EQ. 1) THEN
               call functional_lambda_pbe_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE IF (IS_GLLIMIT .EQ. 1) THEN
               call functional_glpbe_pw_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE
               call functional_pbecp_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            END IF

            !ensure that for the spin-polarized part, EC <= 0
            IF (F_pbes .LT. F_pbe) THEN
               p     = F_pbe
               pra   = D1F_pbe(ID_RA_POS_PBE)
               prb   = D1F_pbe(ID_RB_POS_PBE)
               pgaa  = D1F_pbe(ID_GAA_POS_PBE)
               pgab  = D1F_pbe(ID_GAB_POS_PBE)
               pgbb  = D1F_pbe(ID_GBB_POS_PBE)
            ELSE
               ! here we note that in this case, the contribution from
               ! beta part disppear, so we multiply 0.5
               p     = F_pbes
               pra   = F12*D1F_pbes(ID_RA_POS_PBE)
               prb   = 0.0D0
               pgaa  = F12*D1F_pbes(ID_GAA_POS_PBE)
               pgab  = 0.0D0
               pgbb  = 0.0D0
            END IF

C            write(6,*)"ouralphap is:",p
C            write(6,*)"ouralphapra is:",pra
C            write(6,*)"ouralphaprb is:",prb
C            write(6,*)"ouralphapgaa is:",pgaa
C            write(6,*)"ouralphapgab is:",pgab
C            write(6,*)"ouralphapgbb is:",pgbb

            F            = F            - (F1+c)*w*w*p
            D1F(ID_RA_POS)  = D1F(ID_RA_POS)  - 
     &                     cra*w*w*p    - (F1+c)*F2*w*wra*p  
     &                   - (F1+c)*w*w*pra
            D1F(ID_GAA_POS) = D1F(ID_GAA_POS) - 
     &                     cgaa*w*w*p   - (F1+c)*F2*w*wgaa*p 
     &                   - (F1+c)*w*w*pgaa
            D1F(ID_GAB_POS) = D1F(ID_GAB_POS) - 
     &                     cgab*w*w*p   - (F1+c)*F2*w*wgab*p 
     &                   - (F1+c)*w*w*pgab
            D1F(ID_TA_POS)  = D1F(ID_TA_POS)  - (F1+c)*F2*w*wta*p 
            D1F(ID_RB_POS)  = D1F(ID_RA_POS)
            D1F(ID_GBB_POS) = D1F(ID_GAA_POS)
            D1F(ID_TB_POS)  = D1F(ID_TA_POS)


         END IF ! close shell: if RA > TOL

      ELSE 

         ! now we are entering into the spin-polarized part

         ! differentiation term for rho_alpha/rhoT etc.
         RA    = rhoA
         RB    = rhoB
         RhoT  = RA + RB
         IF (RhoT .LT. TOL) go to 10
         xia   = RA/RhoT
         xib   = RB/RhoT
         xia_a = RB/(RhoT*RhoT)
         xia_b = -F1*RA/(RhoT*RhoT)
         xib_b = RA/(RhoT*RhoT)
         xib_a = -F1*RB/(RhoT*RhoT)


         !-------------------------------------------------------
         ! deal with alpha + beta part
         !-------------------------------------------------------
         ! get the PBE value
         F_pbe    = 0.0D0
         IF (USE_LAMBDA .EQ. 1) THEN
            call functional_lambda_pbe_open(TOL,rhoA,rhoB,DRhoAX,
     &           DRhoAY,DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,F_pbe,D1F_pbe)
         ELSE IF (IS_GLLIMIT .EQ. 1) THEN
            call functional_glpbe_pw_open(TOL,rhoA,rhoB,DRhoAX,
     &           DRhoAY,DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,F_pbe,D1F_pbe)
         ELSE
            call functional_pbecp_open(TOL,rhoA,rhoB,DRhoAX,
     &           DRhoAY,DRhoAZ,DRhoBX,DRhoBY,DRhoBZ,F_pbe,D1F_pbe)
         END IF
         p    =  F_pbe
         pra  =  D1F_pbe(ID_RA_POS_PBE)
         prb  =  D1F_pbe(ID_RB_POS_PBE)
         pgaa =  D1F_pbe(ID_GAA_POS_PBE)
         pgab =  D1F_pbe(ID_GAB_POS_PBE)
         pgbb =  D1F_pbe(ID_GBB_POS_PBE)
C         write(6,*)"ourpbe is:",p
C         write(6,*)"ourpkzbpbe  is:",p
C         write(6,*)"ourpbera is:",pra
C         write(6,*)"ourpberb is:",prb
C         write(6,*)"ourpbegaa is:",pgaa
C         write(6,*)"ourpbegab is:",pgab
C         write(6,*)"ourpbegbb is:",pgbb

         ! deal with the (alpha + beta) part
         ! this is  (1 + c*wt^2)*F_pbe
         F               = F            +    p*(F1 + c*w*w)
         D1F(ID_RA_POS)  = D1F(ID_RA_POS)  + 
     &                  pra*(F1+c*w*w)  + p*cra*w*w  + F2*p*c*w*wra
         D1F(ID_RB_POS)  = D1F(ID_RB_POS)  + 
     &                  prb*(F1+c*w*w)  + p*crb*w*w  + F2*p*c*w*wrb
         D1F(ID_GAA_POS) = D1F(ID_GAA_POS) + 
     &                  pgaa*(F1+c*w*w) + p*cgaa*w*w + F2*p*c*w*wgaa
         D1F(ID_GBB_POS) = D1F(ID_GBB_POS) + 
     &                  pgbb*(F1+c*w*w) + p*cgbb*w*w + F2*p*c*w*wgbb
         D1F(ID_GAB_POS) = D1F(ID_GAB_POS) + 
     &                  pgab*(F1+c*w*w) + p*cgab*w*w + F2*p*c*w*wgab
         D1F(ID_TA_POS)  = D1F(ID_TA_POS)  +    F2*p*c*w*wta 
         D1F(ID_TB_POS)  = D1F(ID_TB_POS)  +    F2*p*c*w*wtb 


         !-------------------------------------------------------
         ! deal with alpha part
         !-------------------------------------------------------
         !this is -(1+c)*wt^2*srho_a*pbea 

         IF (RA .GT. TOL) THEN

            ! get the PBE functional
            F_pbes   = 0.0D0
            IF (USE_LAMBDA .EQ. 1) THEN
               call functional_lambda_pbe_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE IF (IS_GLLIMIT .EQ. 1) THEN
               call functional_glpbe_pw_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE
               call functional_pbecp_open(TOL,rhoA,0.0D0,DRhoAX,
     &           DRhoAY,DRhoAZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            END IF
            IF (F_pbes .LT. F_pbe) THEN
               p     = F_pbe
               pra   = D1F_pbe(ID_RA_POS_PBE)
               prb   = D1F_pbe(ID_RB_POS_PBE)
               pgaa  = D1F_pbe(ID_GAA_POS_PBE)
               pgab  = D1F_pbe(ID_GAB_POS_PBE)
               pgbb  = D1F_pbe(ID_GBB_POS_PBE)
            ELSE
               p     = F_pbes
               pra   = D1F_pbes(ID_RA_POS_PBE)
               prb   = 0.0D0
               pgaa  = D1F_pbes(ID_GAA_POS_PBE)
               pgab  = 0.0D0
               pgbb  = 0.0D0
            END IF
C            write(6,*)"ouralphabe  is:",p
C            write(6,*)"ouralphapra is:",pra
C            write(6,*)"ouralphaprb is:",prb
C            write(6,*)"ouralphapgaa is:",pgaa
C            write(6,*)"ouralphapgab is:",pgab
C            write(6,*)"ouralphapgbb is:",pgbb

            F            = F                   - (F1+c)*w*w*xia*p
            D1F(ID_RA_POS)  = D1F(ID_RA_POS)          
     &                   - cra*w*w*xia*p       - (F1+c)*F2*w*wra*xia*p
     &                   - (F1+c)*w*w*xia_a*p  - (F1+c)*w*w*xia*pra
            D1F(ID_RB_POS)  = D1F(ID_RB_POS)          
     &                   - crb*w*w*xia*p       - (F1+c)*F2*w*wrb*xia*p
     &                   - (F1+c)*w*w*xia_b*p  - (F1+c)*w*w*xia*prb
            D1F(ID_GAA_POS) = D1F(ID_GAA_POS)          
     &                   - cgaa*w*w*xia*p      - (F1+c)*F2*w*wgaa*xia*p
     &                                         - (F1+c)*w*w*xia*pgaa
            D1F(ID_GBB_POS) = D1F(ID_GBB_POS)          
     &                   - cgbb*w*w*xia*p      - (F1+c)*F2*w*wgbb*xia*p
     &                                         - (F1+c)*w*w*xia*pgbb
            D1F(ID_GAB_POS) = D1F(ID_GAB_POS)  
     &                   - cgab*w*w*xia*p      - (F1+c)*F2*w*wgab*xia*p
     &                                         - (F1+c)*w*w*xia*pgab
            D1F(ID_TA_POS)  = D1F(ID_TA_POS)   - (F1+c)*F2*w*wta*xia*p 
            D1F(ID_TB_POS)  = D1F(ID_TB_POS)   - (F1+c)*F2*w*wtb*xia*p 

         END IF ! RA > TOL



         !-------------------------------------------------------
         ! deal with beta part
         !-------------------------------------------------------
         !this is -(1+c)*wt^2*srho_b*pbeb
         IF (RB .GT. TOL) THEN

            ! get the PBE functional
            F_pbes   = 0.0D0
            IF (USE_LAMBDA .EQ. 1) THEN
               call functional_lambda_pbe_open(TOL,rhoB,0.0D0,DRhoBX,
     &           DRhoBY,DRhoBZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE IF (IS_GLLIMIT .EQ. 1) THEN
               call functional_glpbe_pw_open(TOL,rhoB,0.0D0,DRhoBX,
     &           DRhoBY,DRhoBZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            ELSE
               call functional_pbecp_open(TOL,rhoB,0.0D0,DRhoBX,
     &           DRhoBY,DRhoBZ,0.0D0,0.0D0,0.0D0,F_pbes,D1F_pbes)
            END IF

            ! here we have to make some explanations. Since in the above
            ! codes we have taken the beta variable as alpha variable
            ! in the lambda_pbe, so here the assignment should be
            ! switched 
            IF (F_pbes .LT. F_pbe) THEN
               p     = F_pbe
               pra   = D1F_pbe(ID_RB_POS_PBE)
               prb   = D1F_pbe(ID_RA_POS_PBE)
               pgaa  = D1F_pbe(ID_GBB_POS_PBE)
               pgab  = D1F_pbe(ID_GAB_POS_PBE)
               pgbb  = D1F_pbe(ID_GAA_POS_PBE)
            ELSE
               p     = F_pbes
               pra   = 0.0D0
               prb   = D1F_pbes(ID_RA_POS_PBE)
               pgaa  = 0.0D0
               pgab  = 0.0D0
               pgbb  = D1F_pbes(ID_GAA_POS_PBE)
            END IF
C            write(6,*)"ourbetapbe is:",p
C            write(6,*)"ourbetapra is:",pra
C            write(6,*)"ourbetaprb is:",prb
C            write(6,*)"ourbetapgaa is:",pgaa
C            write(6,*)"ourbetapgab is:",pgab
C            write(6,*)"ourbetapgbb is:",pgbb

            F            = F                   - (F1+c)*w*w*xib*p
            D1F(ID_RA_POS)  = D1F(ID_RA_POS)          
     &                   - cra*w*w*xib*p       - (F1+c)*F2*w*wra*xib*p
     &                   - (F1+c)*w*w*xib_a*p  - (F1+c)*w*w*xib*pra
            D1F(ID_RB_POS)  = D1F(ID_RB_POS)          
     &                   - crb*w*w*xib*p       - (F1+c)*F2*w*wrb*xib*p
     &                   - (F1+c)*w*w*xib_b*p  - (F1+c)*w*w*xib*prb
            D1F(ID_GAA_POS) = D1F(ID_GAA_POS)          
     &                   - cgaa*w*w*xib*p      - (F1+c)*F2*w*wgaa*xib*p
     &                                         - (F1+c)*w*w*xib*pgaa
            D1F(ID_GBB_POS) = D1F(ID_GBB_POS)          
     &                   - cgbb*w*w*xib*p      - (F1+c)*F2*w*wgbb*xib*p
     &                                         - (F1+c)*w*w*xib*pgbb
            D1F(ID_GAB_POS) = D1F(ID_GAB_POS)  
     &                   - cgab*w*w*xib*p      - (F1+c)*F2*w*wgab*xib*p
     &                                         - (F1+c)*w*w*xib*pgab
            D1F(ID_TA_POS)  = D1F(ID_TA_POS)   - (F1+c)*F2*w*wta*xib*p 
            D1F(ID_TB_POS)  = D1F(ID_TB_POS)   - (F1+c)*F2*w*wtb*xib*p 


         END IF ! RB > TOL


      END IF ! nDen = 1

10    continue
      END
