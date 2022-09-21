
      SUBROUTINE b13coor_par(method,INFOR,NDEN,NA,NB,NG,THRESH,
     $VAL_P,VAL_Q,
     $HIRWTS,RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *   This is modified BR94 correlation functional                 *
c    *   which defined in the paper of:                               *
c    *   "Density functionals for static, dynamical, and strong       *
c    *    correlation "                                               *
c    *   A. D. Becke, J. Chem. Phys. 138, 074109 (2013)               *
c    *   this is for the parallel spin component                      *         
c    *   see equation 37,38,39 and 40                                 *         
c    *                                                                *
c    *   There are two suites of correlation functionals for B13      *
c    *   one is described by equation 37/38 (method = 1 above),       *
c    *   this one is lambda integrated form. The other one is         *
c    *   defined by equation 39/40, which is lambda=1 form            *     
c    *                                                                *
c    *   additionaly, it's worthy to note that for the parallel       *
c    *   part of b13coor, it depends on f_{sigma,sigma} which is      *
c    *   defined in equation 23 of above paper. This is related to    *
c    *   the static correlation part. We note, that this f is not     *
c    *   the f factor defined in becke05ex.f. The f_{sigma,sigma}     *
c    *   is the min(A1/A2,1) where A1/A2 are computed from the        *
c    *   becke05_A_sigma (see that function for more information).    *    
c    *                                                                *
c    *   the original Br94 correlation function is in paper below:    *         
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *                                                                *
c    *   input:                                                       *
c    *   method:     use equation 37/38(method == 1) or 39/40?        *         
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exchange energy density                          *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB,NE
      INTEGER method
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 VAL_P,VAL_Q
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_EXA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXB_POS

      ! now let's define the parameter
      REAL*8 CAA,CBB
      REAL*8 CAA2,CBB2

      ! variables for Neff and it's derivatives
      REAL*8 N1,N_RA,N_GAA,N_TA,N_LA,N_UA    
      REAL*8 N2,N_RB,N_GBB,N_TB,N_LB,N_UB    

      ! variables related to A1 and A2
      INTEGER iSpin
      REAL*8 A1
      REAL*8 DA1_RA,DA1_GAA,DA1_TA,DA1_LA,DA1_UA
      REAL*8 DA1_RB,DA1_GBB,DA1_TB,DA1_LB,DA1_UB
      REAL*8 A2,DA2DR,DA2DG,DA2DT,R2,R3

      ! variables for f and it's derivatives
      REAL*8 fac,f_RA,f_GAA,f_TA,f_LA,f_UA    
      REAL*8 f_RB,f_GBB,f_TB,f_LB,f_UB    
      REAL*8 A2_2

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAA,ZAA_RA,ZAA_GAA,ZAA_TA,ZAA_LA,ZAA_UA    
      REAL*8 ZBB,ZBB_RB,ZBB_GBB,ZBB_TB,ZBB_LB,ZBB_UB    
      REAL*8 ZAA2,ZAA3,ZAA4,ZAA5
      REAL*8 ZBB2,ZBB3,ZBB4,ZBB5
      REAL*8 FZ0,FZ1,FZ,FZ_Z,T,T2

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAA,ECAA_RA,ECAA_GAA,ECAA_TA,ECAA_LA,ECAA_UA  
      REAL*8 ECAA_RB,ECAA_GBB,ECAA_TB,ECAA_LB,ECAA_UB  
      REAL*8 ECBB,ECBB_RB,ECBB_GBB,ECBB_TB,ECBB_LB,ECBB_UB  
      REAL*8 ECBB_RA,ECBB_GAA,ECBB_TA,ECBB_LA,ECBB_UA

      ! constants
      REAL*8 ZERO,ONE,TWO,THREE,FOUR,FIVE,F12,F14,F13,F16,F112

C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! constants
      ZERO = 0.0D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      THREE= 3.0D0
      FOUR = 4.0D0
      FIVE = 5.0D0
      F14  = 0.25D0
      F12  = 0.5D0
      F13  = 1.0D0/3.0D0
      F16  = 1.0D0/6.0D0
      F112 = 1.0D0/12.0D0

      ! constants defined in the paper
      CAA   = -0.88D0
      CBB   = -0.88D0
      CAA2  = -0.005D0
      CBB2  = -0.005D0
      IF (method .eq. 1) THEN
         CAA2  = -0.01D0
         CBB2  = -0.01D0
      END IF

      ! set the total number of electrons
      NE = NA + NB

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! loop over NG
      DO I = 1,NG

        ! variables
        RA  = RhoA(i)
        RB  = RhoB(i)
        GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)
C         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
C         write(6,*)"D1RB(i,1)", DRB(i,1)
C         write(6,*)"D1RB(i,2)", DRB(i,2)
C         write(6,*)"D1RB(i,3)", DRB(i,3)
         
         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! EC_{alpha,alpha}
         ECAA     = ZERO 
         ECAA_RA  = ZERO 
         ECAA_GAA = ZERO 
         ECAA_TA  = ZERO 
         ECAA_LA  = ZERO 
         ECAA_UA  = ZERO 
         ECAA_RB  = ZERO 
         ECAA_GBB = ZERO 
         ECAA_TB  = ZERO 
         ECAA_LB  = ZERO 
         ECAA_UB  = ZERO 
         IF (NA > 1) THEN
            IF (RA > THRESH) THEN

               ! Neff for alpha
               N1       = ZERO 
               N_RA     = ZERO 
               N_GAA    = ZERO 
               N_TA     = ZERO 
               N_LA     = ZERO 
               N_UA     = ZERO 
               CALL becke05_neff(THRESH,RA,GAA,TA,LA,UA,HIRWT,N1,
     $ N_RA,N_GAA,N_TA,N_LA,N_UA)

               ! now it's ZAA
               ZAA       = ZERO 
               ZAA_RA    = ZERO 
               ZAA_GAA   = ZERO 
               ZAA_TA    = ZERO 
               ZAA_LA    = ZERO 
               ZAA_UA    = ZERO 
               IF (ABS(UA) > THRESH) THEN
                  ZAA     = TWO*CAA*N1*RA/UA
                  ZAA_RA  = TWO*CAA*(N_RA*RA+N1)/UA
                  ZAA_GAA = TWO*CAA*N_GAA*RA/UA
                  ZAA_TA  = TWO*CAA*N_TA*RA/UA
                  ZAA_LA  = TWO*CAA*N_LA*RA/UA
                  ZAA_UA  = TWO*CAA*RA*(N_UA/UA-N1/(UA*UA))
               END IF
c               write(6,*)"ZAA",ZAA
c               write(6,*)"ZAA_RA",ZAA_RA
   
               ! now it's function of ZAA
               FZ   = ZERO
               FZ_Z = ZERO
               IF (ABS(ZAA) > THRESH) THEN

C                 method is 1 means we use equation 37/38
C                 else we se equation 39/40                  
                  IF (method .eq. 1) THEN
                     ZAA2 = ZAA*ZAA
                     ZAA3 = ZAA2*ZAA
                     ZAA4 = ZAA3*ZAA
                     FZ0  = DLOG(ONE+ZAA/TWO)
                     FZ1  = ONE-(TWO/ZAA)*FZ0
                     FZ   = ZAA4*FZ1
                     FZ_Z = FOUR*ZAA3*FZ1 + ZAA4*(TWO*FZ0/ZAA2 -
     $                   ONE/(ZAA*(ONE+ZAA/TWO)))
                  ELSE
                     ZAA2 = ZAA*ZAA
                     ZAA3 = ZAA2*ZAA
                     ZAA4 = ZAA3*ZAA
                     ZAA5 = ZAA4*ZAA
                     T    = ONE+ZAA/TWO
                     T2   = T*T
                     FZ   = ZAA5/T
                     FZ_Z = FIVE*ZAA4/T-F12*ZAA5/T2
                  END IF
   
               END IF
c               write(6,*)"FZ for ZAA",FZ
   
               ! D and it's derivatives
               D    = ZERO
               DDDR = ZERO
               DDDG = ZERO
               DDDT = ONE
               IF (RA > THRESH .and. NA > 1) THEN
                  D    = TA - F14*GAA/RA
                  DDDR = F14*GAA/(RA*RA)
                  DDDG = -F14/RA
               END IF
C               write(6,*)"D",D
c               write(6,*)"D",D
c               write(6,*)"D_RA",DDDR

               ! do A2
               A2    = ZERO
               DA2DR = ZERO
               DA2DG = ZERO
               DA2DT = ZERO
               IF (RA > THRESH .and. NA > 1) THEN
                  A2    = F13*D/RA
                  DA2DT = F13/RA
                  R2    = RA*RA
                  IF (R2 > THRESH) THEN
                     DA2DR = -F13*TA/R2
                     DA2DG = -F112/R2
                     R3 = R2*RA
                     IF (R3 > THRESH) THEN
                        DA2DR = DA2DR + F16*GAA/R3
                     END IF
                  END IF
               END IF

               ! do A1
               iSpin = 1
               CALL becke05_A_sigma(iSpin,NA,NE,THRESH,
     $ VAL_P,VAL_Q,HIRWT,
     $ RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A1,
     $ DA1_RA,DA1_GAA,DA1_TA,DA1_LA,DA1_UA,
     $ DA1_RB,DA1_GBB,DA1_TB,DA1_LB,DA1_UB)

               ! f_sigma_sigma
               fac   = ZERO
               f_RA  = ZERO
               f_GAA = ZERO
               f_TA  = ZERO
               f_LA  = ZERO
               f_UA  = ZERO
               f_RB  = ZERO
               f_GBB = ZERO
               f_TB  = ZERO
               f_LB  = ZERO
               f_UB  = ZERO
               IF (DABS(A2) > THRESH) THEN
                  fac   = A1/A2
                  A2_2  = A2*A2
                  f_RA  = DA1_RA/A2  - A1*DA2DR/A2_2
                  f_GAA = DA1_GAA/A2 - A1*DA2DG/A2_2 
                  f_TA  = DA1_TA/A2  - A1*DA2DT/A2_2 
                  f_LA  = DA1_LA/A2
                  f_UA  = DA1_UA/A2 
                  f_RB  = DA1_RB/A2
                  f_GBB = DA1_GBB/A2
                  f_TB  = DA1_TB/A2
                  f_LB  = DA1_LB/A2
                  f_UB  = DA1_UB/A2
               END IF

               ! now let's assemble all of pices together
               ECAA      = CAA2*(ONE-fac)*RA*D*FZ
               ECAA_RA   = CAA2*(D*FZ+RA*DDDR*FZ+RA*D*FZ_Z*ZAA_RA)
     $*(ONE-fac) - CAA2*f_RA*RA*D*FZ
               ECAA_GAA  = CAA2*RA*(DDDG*FZ+D*FZ_Z*ZAA_GAA)*(ONE-fac)
     $ - CAA2*f_GAA*RA*D*FZ
               ECAA_TA   = CAA2*RA*(DDDT*FZ+D*FZ_Z*ZAA_TA)*(ONE-fac)
     $ - CAA2*f_TA*RA*D*FZ
               ECAA_LA   = CAA2*RA*D*FZ_Z*ZAA_LA*(ONE-fac)
     $ - CAA2*f_LA*RA*D*FZ
               ECAA_UA   = CAA2*RA*D*FZ_Z*ZAA_UA*(ONE-fac)
     $ - CAA2*f_UA*RA*D*FZ
               ECAA_RB   = -CAA2*f_RB*RA*D*FZ
               ECAA_GBB  = -CAA2*f_GBB*RA*D*FZ
               ECAA_TB   = -CAA2*f_TB*RA*D*FZ
               ECAA_LB   = -CAA2*f_LB*RA*D*FZ
               ECAA_UB   = -CAA2*f_UB*RA*D*FZ
c               write(6,*)"ECAA_GBB",ECAA_GBB
            END IF
         END IF

         ! EC_{beta,beta}
         ECBB     = ZERO
         ECBB_RB  = ZERO
         ECBB_GBB = ZERO
         ECBB_TB  = ZERO
         ECBB_LB  = ZERO
         ECBB_UB  = ZERO 
         ECBB_RA  = ZERO
         ECBB_GAA = ZERO
         ECBB_TA  = ZERO
         ECBB_LA  = ZERO
         ECBB_UA  = ZERO 
         IF (NB > 1) THEN
            IF (RB > THRESH) THEN

               ! compute NEff
               N2       = ZERO 
               N_RB     = ZERO
               N_GBB    = ZERO
               N_TB     = ZERO
               N_LB     = ZERO
               N_UB     = ZERO
               CALL becke05_neff(THRESH,RB,GBB,TB,LB,UB,HIRWT,N2,
     $ N_RB,N_GBB,N_TB,N_LB,N_UB)
c               write(6,*)"Neff GBB", N_GBB
   
               ! now it's ZBB
               ZBB       = ZERO 
               ZBB_RB    = ZERO 
               ZBB_GBB   = ZERO 
               ZBB_TB    = ZERO 
               ZBB_LB    = ZERO 
               ZBB_UB    = ZERO 
               IF (ABS(UB) > THRESH) THEN
                  ZBB     = TWO*CBB*N2*RB/UB
                  ZBB_RB  = TWO*CBB*(N_RB*RB+N2)/UB
                  ZBB_GBB = TWO*CBB*N_GBB*RB/UB
                  ZBB_TB  = TWO*CBB*N_TB*RB/UB
                  ZBB_LB  = TWO*CBB*N_LB*RB/UB
                  ZBB_UB  = TWO*CBB*RB*(N_UB/UB-N2/(UB*UB))
               END IF
c               write(6,*)"ZBB",ZBB
c               write(6,*)"ZBB_GBB",ZBB_GBB
   
               ! now it's function of ZBB
               FZ   = ZERO
               FZ_Z = ZERO
               IF (ABS(ZBB) > THRESH) THEN

                  IF (method .eq. 1) THEN
                     ZBB2 = ZBB*ZBB
                     ZBB3 = ZBB2*ZBB
                     ZBB4 = ZBB3*ZBB
                     FZ0  = DLOG(ONE+ZBB/TWO)
                     FZ1  = ONE-(TWO/ZBB)*FZ0
                     FZ   = ZBB4*FZ1
                     FZ_Z = FOUR*ZBB3*FZ1 + ZBB4*(TWO*FZ0/ZBB2 -
     $                   ONE/(ZBB*(ONE+ZBB/TWO)))
                  ELSE
                     ZBB2 = ZBB*ZBB
                     ZBB3 = ZBB2*ZBB
                     ZBB4 = ZBB3*ZBB
                     ZBB5 = ZBB4*ZBB
                     T    = ONE+ZBB/TWO
                     T2   = T*T
                     FZ   = ZBB5/T
                     FZ_Z = FIVE*ZBB4/T-F12*ZBB5/T2
                  END IF
               END IF
c               write(6,*)"FZ for ZBB",FZ
   
               ! D and it's derivatives
               D    = ZERO
               DDDR = ZERO
               DDDG = ZERO
               DDDT = ONE
               IF (RB > THRESH .and. NB > 1) THEN
                  D    = TB - F14*GBB/RB
                  DDDR = F14*GBB/(RB*RB)
                  DDDG = -F14/RB
               END IF
c              write(6,*)"D GBB",DDDG
   
               ! do A2
               A2    = ZERO
               DA2DR = ZERO
               DA2DG = ZERO
               DA2DT = ZERO
               IF (RB > THRESH .and. NB > 1) THEN
                  A2    = F13*D/RB
                  DA2DT = F13/RB
                  R2    = RB*RB
                  IF (R2 > THRESH) THEN
                     DA2DR = -F13*TB/R2
                     DA2DG = -F112/R2
                     R3 = R2*RB
                     IF (R3 > THRESH) THEN
                        DA2DR = DA2DR + F16*GBB/R3
                     END IF
                  END IF
               END IF
c               write(6,*)"A2 GBB",DA2DG

               ! do A1
               iSpin = 2
               CALL becke05_A_sigma(iSpin,NB,NE,THRESH,
     $ VAL_P,VAL_Q,HIRWT,
     $ RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A1,
     $ DA1_RA,DA1_GAA,DA1_TA,DA1_LA,DA1_UA,
     $ DA1_RB,DA1_GBB,DA1_TB,DA1_LB,DA1_UB)
c               write(6,*)"A1 GBB",DA1_GBB

               ! f_sigma_sigma
               fac   = ZERO
               f_RA  = ZERO
               f_GAA = ZERO
               f_TA  = ZERO
               f_LA  = ZERO
               f_UA  = ZERO
               f_RB  = ZERO
               f_GBB = ZERO
               f_TB  = ZERO
               f_LB  = ZERO
               f_UB  = ZERO
               IF (DABS(A2) > THRESH) THEN
                  fac   = A1/A2
                  A2_2  = A2*A2
                  f_RA  = DA1_RA/A2
                  f_GAA = DA1_GAA/A2
                  f_TA  = DA1_TA/A2
                  f_LA  = DA1_LA/A2
                  f_UA  = DA1_UA/A2 
                  f_RB  = DA1_RB/A2  - A1*DA2DR/A2_2
                  f_GBB = DA1_GBB/A2 - A1*DA2DG/A2_2 
                  f_TB  = DA1_TB/A2  - A1*DA2DT/A2_2 
                  f_LB  = DA1_LB/A2
                  f_UB  = DA1_UB/A2
               END IF
c               write(6,*)"fac GBB",f_GBB

               ! now let's assemble all of pices together
               ECBB      = CBB2*RB*D*FZ*(ONE-fac)
               ECBB_RB   = CBB2*(D*FZ+RB*DDDR*FZ+RB*D*FZ_Z*ZBB_RB)
     $*(ONE-fac) - CBB2*RB*D*FZ*f_RB
               ECBB_GBB  = CBB2*RB*(DDDG*FZ+D*FZ_Z*ZBB_GBB)*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_GBB
               ECBB_TB   = CBB2*RB*(DDDT*FZ+D*FZ_Z*ZBB_TB)*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_TB
               ECBB_LB   = CBB2*RB*D*FZ_Z*ZBB_LB*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_LB
               ECBB_UB   = CBB2*RB*D*FZ_Z*ZBB_UB*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_UB
               ECBB_RA   = -CBB2*RB*D*FZ*f_RA
               ECBB_GAA  = -CBB2*RB*D*FZ*f_GAA
               ECBB_TA   = -CBB2*RB*D*FZ*f_TA
               ECBB_LA   = -CBB2*RB*D*FZ*f_LA
               ECBB_UA   = -CBB2*RB*D*FZ*f_UA
c               write(6,*)"ECBB_GBB",ECBB_GBB
            END IF
         END IF

         ! collect energy
         F(i) = F(i) + ECAA + ECBB
c         write(6,*)"B13 EPAR", ECAA+ECBB

         ! collect all of terms for alpha derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAA_RA + ECBB_RA
c         write(6,*)"B13 PAR, RA", ECAA_RA + ECBB_RA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAA_GAA + ECBB_GAA
c         write(6,*)"B13 PAR, GAA", ECAA_GAA + ECBB_GAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAA_TA + ECBB_TA
c         write(6,*)"B13 PAR, TA", ECAA_TA + ECBB_TA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAA_LA + ECBB_LA
c         write(6,*)"B13 PAR, LA", ECAA_LA + ECBB_LA
         ID_EXA_POS=D1VARS(ID_EXA)
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + ECAA_UA + ECBB_UA
c         write(6,*)"B13 PAR, EXA", ECAA_UA + ECBB_UA
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECBB_RB + ECAA_RB
c         write(6,*)"B13 PAR, RB", ECAA_RB + ECBB_RB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECBB_GBB + ECAA_GBB
c         write(6,*)"B13 PAR, GBB", ECAA_GBB + ECBB_GBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECBB_TB + ECAA_TB
c         write(6,*)"B13 PAR, TB", ECAA_TB + ECBB_TB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECBB_LB + ECAA_LB
c         write(6,*)"B13 PAR, LB", ECAA_LB + ECBB_LB
         ID_EXB_POS=D1VARS(ID_EXB)
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + ECBB_UB + ECAA_UB
c         write(6,*)"B13 PAR, EXB", ECAA_UB + ECBB_UB

      END DO

      RETURN 
      END

      SUBROUTINE b13coor_opp(method,INFOR,NDEN,NA,NB,NG,THRESH,
     $VAL_P,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *   This is modified BR94 correlation functional                 *
c    *   which defined in the paper of:                               *
c    *   "Density functionals for static, dynamical, and strong       *
c    *    correlation "                                               *
c    *   A. D. Becke, J. Chem. Phys. 138, 074109 (2013)               *
c    *   this is for the opposite spin component                      *         
c    *   see equation 37,38,39 and 40                                 *         
c    *                                                                *
c    *   There are two suites of correlation functionals for B13      *
c    *   one is described by equation 37/38 (method = 1 above),       *
c    *   this one is lambda integrated form. The other one is         *
c    *   defined by equation 39/40, which is lambda=1 form            *     
c    *                                                                *
c    *   the original Br94 correlation function is in paper below:    *         
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *                                                                *
c    *   input:                                                       *
c    *   method:     use equation 37/38(method == 1) or 39/40?        *         
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exchange energy density                          *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER method
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 VAL_P
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! now let's define the parameter
      REAL*8 CAB,CAB2

      ! variables for Neff and it's derivatives
      REAL*8 N1,N_RA,N_GAA,N_TA,N_LA,N_UA    
      REAL*8 N2,N_RB,N_GBB,N_TB,N_LB,N_UB    

      ! variables for f and it's derivatives
      REAL*8 fac,f_RA,f_GAA,f_TA,f_LA,f_UA    
      REAL*8 f_RB,f_GBB,f_TB,f_LB,f_UB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAB,ZAB_RA,ZAB_GAA,ZAB_TA,ZAB_LA,ZAB_UA    
      REAL*8     ZAB_RB,ZAB_GBB,ZAB_TB,ZAB_LB,ZAB_UB    
      REAL*8 ZAB2,ZAB3,T,T2
      REAL*8 FZ0,FZ1,FZ,FZ_Z

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA,ECAB_UA  
      REAL*8      ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB,ECAB_UB  

      ! constant
      REAL*8 ZERO,ONE,TWO,THREE

C       
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C

       ! constants
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0

       ! constants defined in the paper
       ! we note that both equation 37/38 and 39/40 shares
       ! the same CAB2
       CAB   = -0.63D0
       CAB2  = -0.8D0

       ! firstly initilize variable position information
       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
       ! loop over NG
      DO I = 1,NG

         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)
c         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
         
         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! we can compute f in advance
         CALL becke05_f(THRESH,VAL_P,HIRWT,RA,GAA,TA,LA,UA,
     $ RB,GBB,TB,LB,UB,fac,f_RA,f_GAA,f_TA,f_LA,f_UA,
     $ f_RB,f_GBB,f_TB,f_LB,f_UB)
   
         ! N alpha
         N1       = ZERO 
         N_RA     = ZERO 
         N_GAA    = ZERO 
         N_TA     = ZERO 
         N_LA     = ZERO 
         N_UA     = ZERO 
         IF (RA > THRESH) THEN
            CALL becke05_neff(THRESH,RA,GAA,TA,LA,UA,HIRWT,N1,
     $ N_RA,N_GAA,N_TA,N_LA,N_UA)
         END IF

         ! N beta
         N2       = ZERO
         N_RB     = ZERO
         N_GBB    = ZERO
         N_TB     = ZERO
         N_LB     = ZERO
         N_UB     = ZERO
         IF (RB > THRESH) THEN
            CALL becke05_neff(THRESH,RB,GBB,TB,LB,UB,HIRWT,N2,
     $ N_RB,N_GBB,N_TB,N_LB,N_UB)
         END IF
   
         ! now it's ECAB
         ECAB     = ZERO
         ECAB_RA  = ZERO
         ECAB_GAA = ZERO
         ECAB_TA  = ZERO
         ECAB_LA  = ZERO
         ECAB_UA  = ZERO
         ECAB_RB  = ZERO
         ECAB_GBB = ZERO
         ECAB_TB  = ZERO
         ECAB_LB  = ZERO
         ECAB_UB  = ZERO
         IF (RA > THRESH .and. RB > THRESH .and. 
     $       NA > 0 .and. NB > 0) THEN

            ! now it's ZAB
            ZAB      = ZERO
            ZAB_RA   = ZERO
            ZAB_GAA  = ZERO
            ZAB_TA   = ZERO
            ZAB_LA   = ZERO
            ZAB_UA   = ZERO
            ZAB_RB   = ZERO
            ZAB_GBB  = ZERO
            ZAB_TB   = ZERO
            ZAB_LB   = ZERO
            ZAB_UB   = ZERO
            IF (ABS(UA) > THRESH) THEN
               ZAB     = CAB*N1*RA/UA
               ZAB_RA  = CAB*(N_RA*RA+N1)/UA
               ZAB_GAA = CAB*N_GAA*RA/UA
               ZAB_TA  = CAB*N_TA*RA/UA
               ZAB_LA  = CAB*N_LA*RA/UA
               ZAB_UA  = CAB*RA*(N_UA/UA-N1/(UA*UA))
            END IF
            IF (ABS(UB) > THRESH) THEN
               ZAB     = ZAB + CAB*N2*RB/UB
               ZAB_RB  = CAB*(N_RB*RB+N2)/UB
               ZAB_GBB = CAB*N_GBB*RB/UB
               ZAB_TB  = CAB*N_TB*RB/UB
               ZAB_LB  = CAB*N_LB*RB/UB
               ZAB_UB  = CAB*RB*(N_UB/UB-N2/(UB*UB))
            END IF
c            write(6,*)"ZAB",ZAB

            ! now it's function of ZAB
            FZ   = ZERO
            FZ_Z = ZERO
            IF (ABS(ZAB) > THRESH) THEN

               if (method .eq. 1) THEN
                  FZ0  = DLOG(ONE+ZAB)/ZAB
                  FZ   = ZAB*ZAB*(ONE-FZ0)
                  FZ_Z = TWO*ZAB*(ONE-FZ0) + ZAB*ZAB*(FZ0/ZAB -
     $                ONE/((ONE+ZAB)*ZAB))
               ELSE
                  T    = ONE+ZAB
                  T2   = T*T
                  ZAB2 = ZAB*ZAB
                  ZAB3 = ZAB*ZAB*ZAB
                  FZ   = ZAB3/T
                  FZ_Z = THREE*ZAB2/T-ZAB3/T2
               END IF

            END IF
c            write(6,*)"FZ for ZAB",FZ

            ! now let's assemble all of pices together
            ECAB      = CAB2*RA*RB*FZ*(ONE-fac)
            ECAB_RA   = CAB2*RB*(FZ+RA*FZ_Z*ZAB_RA)*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_RA
            ECAB_RB   = CAB2*RA*(FZ+RB*FZ_Z*ZAB_RB)*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_RB
            ECAB_GAA  = CAB2*RA*RB*FZ_Z*ZAB_GAA*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_GAA
            ECAB_GBB  = CAB2*RA*RB*FZ_Z*ZAB_GBB*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_GBB
            ECAB_TA   = CAB2*RA*RB*FZ_Z*ZAB_TA*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_TA
            ECAB_TB   = CAB2*RA*RB*FZ_Z*ZAB_TB*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_TB
            ECAB_LA   = CAB2*RA*RB*FZ_Z*ZAB_LA*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_LA
            ECAB_LB   = CAB2*RA*RB*FZ_Z*ZAB_LB*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_LB
            ECAB_UA   = CAB2*RA*RB*FZ_Z*ZAB_UA*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_UA
            ECAB_UB   = CAB2*RA*RB*FZ_Z*ZAB_UB*(ONE-fac)
     $ -CAB2*RA*RB*FZ*f_UB
         END IF
c         write(6,*)"ECAB",ECAB

         ! collect energy
         F(i) = F(i) + ECAB 

         ! collect all of terms for alpha derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAB_RA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAB_GAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAB_TA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAB_LA
         ID_EXA_POS=D1VARS(ID_EXA)
         D1F(i, ID_EXA_POS)  = D1F(i, ID_EXA_POS) + ECAB_UA
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECAB_RB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECAB_GBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECAB_TB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECAB_LB
         ID_EXB_POS=D1VARS(ID_EXB)
         D1F(i, ID_EXB_POS)  = D1F(i, ID_EXB_POS) + ECAB_UB

      END DO

      RETURN 
      END
