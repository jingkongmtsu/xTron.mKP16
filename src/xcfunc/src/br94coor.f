
      SUBROUTINE br94coor_par(INFOR,NDEN,NA,NB,NG,THRESH,
     $VAL_P,VAL_Q,HIRWTS,RhoA,RhoB,
     $DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *   This is BR94 correlation functional which defined in:        *
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *   this is for the parallel spin component                      *         
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
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
      INTEGER NA,NB,iSpin,NE
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 RA,RB,A2,A1,R2,DA2DR,DA2DT,DA2DG
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB,DA1_RA,DA1_GAA,DA1_TA,DA1_LA,DA1_UA
      REAL*8 DA1_RB,DA1_GBB,DA1_TB,DA1_LB,DA1_UB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS

      ! now let's define the parameter
      REAL*8 CAA,CBB,CAB,A2_2
      REAL*8 CAA2,CBB2,CAB2

      ! variables for U and it's derivatives
      REAL*8 UA,UA_RA,UA_GAA,UA_TA,UA_LA,UUA    
      REAL*8 UB,UB_RB,UB_GBB,UB_TB,UB_LB,UUB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAA,ZAA_RA,ZAA_GAA,ZAA_TA,ZAA_LA
      REAL*8 ZBB,ZBB_RB,ZBB_GBB,ZBB_TB,ZBB_LB
      REAL*8 ZAA2,ZAA3,ZAA4
      REAL*8 ZBB2,ZBB3,ZBB4
      REAL*8 ZAA_U,ZAB_U,ZBB_U
      REAL*8 FZ0,FZ1,FZ,FZ_Z,R3
      REAL*8 ZERO,ONE,TWO,THREE,FOUR,FIVE,F12,F14,F13,F16,F112
      REAL*8 ONEDU

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAA,ECAA_RA,ECAA_GAA,ECAA_TA,ECAA_LA,ECAA_UUA,ECAA_GBB
      REAL*8 ECBB,ECBB_RB,ECBB_GBB,ECBB_TB,ECBB_LB,ECBB_UUB,ECAA_LB  
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA,ECBB_UUA,ECAA_RB  
      REAL*8 ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB,ECAA_UUB,ECAA_TB  
      REAL*8 ECBB_GAA,ECBB_RA,ECBB_LA,ECBB_TA,VAL_P,VAL_Q
      REAL*8 fac,f_RA,f_GAA,f_TA,f_LA,f_UA,f_RB,f_GBB,f_TB,f_LB,f_UB

C       
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C

       ! constants defined in the paper
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
       CAA   = -0.88D0
       CBB   = -0.88D0
       CAB   = -0.63D0
       CAA2  = -0.01D0
       CBB2  = -0.01D0
       CAB2  = -0.8D0
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
         UUA  = EXA(i)
         UUB  = EXB(i)
         HIRWT = HIRWTS(i)
c         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
         
         ! EC_{alpha,alpha}
         UA       = 0.0D0
         UA_RA    = 0.0D0
         UA_GAA   = 0.0D0
         UA_TA    = 0.0D0
         UA_LA    = 0.0D0
         ECAA     = 0.0D0
         ECAA_RA  = 0.0D0
         ECAA_GAA = 0.0D0
         ECAA_TA  = 0.0D0
         ECAA_LA  = 0.0D0
         ECAA_UUA  = 0.0D0
         ECAA_UUB  = 0.0D0
         IF (NA > 1) THEN
            IF (RA > THRESH) THEN
   
               ! BR89 exchange hole
         CALL BR89HOLE(THRESH,RA,GAA,TA,LA,UA,UA_RA,UA_GAA,
     $                       UA_TA,UA_LA)
   
               ! now it's ZAA
               ZAA       = 0.0D0
               ZAA_RA    = 0.0D0
               ZAA_GAA   = 0.0D0
               ZAA_TA    = 0.0D0
               ZAA_LA    = 0.0D0
c               write(6,*)"UA",UA
               IF (ABS(UA) > THRESH) THEN
                  ONEDU   = 1.0D0/UA
                  ZAA     = CAA*2.0D0*ONEDU
                  ZAA_U   = -2.0D0*CAA/(UA*UA)
                  ZAA_RA  = ZAA_U*UA_RA
                  ZAA_GAA = ZAA_U*UA_GAA
                  ZAA_TA  = ZAA_U*UA_TA
                  ZAA_LA  = ZAA_U*UA_LA
               END IF
c               write(6,*)"ZAA",ZAA
c               write(6,*)"ZAA_RA",ZAA_RA
   
               ! now it's function of ZAA
               FZ   = 0.0D0
               FZ_Z = 0.0D0
               IF (ABS(ZAA) > THRESH) THEN
                  ZAA2 = ZAA*ZAA
                  ZAA3 = ZAA2*ZAA
                  ZAA4 = ZAA3*ZAA
                  FZ0  = DLOG(1.0D0+ZAA/2.0D0)
                  FZ1  = 1.0D0-(2.0D0/ZAA)*FZ0
                  FZ   = ZAA4*FZ1
                  FZ_Z = 4.0D0*ZAA3*FZ1 + ZAA4*(2.0D0*FZ0/ZAA2 -
     $                   1.0D0/(ZAA*(1.0D0+ZAA/2.0D0)))
   
               END IF
c               write(6,*)"FZ for ZAA",FZ
   
               ! D and it's derivatives
               D    = 0.0D0
               DDDR = 0.0D0
               DDDG = 0.0D0
               DDDT = 0.0D0
               IF (RA > THRESH) THEN
                  D    = TA - 0.25D0*GAA/RA
                  DDDR = 0.25D0*GAA/(RA*RA)
                  DDDG = -0.25D0/RA
                  DDDT = 1.0D0
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
     $ RA,GAA,TA,LA,UUA,RB,GBB,TB,LB,UUB,A1,
     $ DA1_RA,DA1_GAA,DA1_TA,DA1_LA,DA1_UA,
     $ DA1_RB,DA1_GBB,DA1_TB,DA1_LB,DA1_UB)
c
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
             IF(fac.le.ONE) then 
               ECAA      = CAA2*(ONE-fac)*RA*D*FZ
               ECAA_RA   = CAA2*(D*FZ+RA*DDDR*FZ+RA*D*FZ_Z*ZAA_RA)
     $*(ONE-fac) - CAA2*f_RA*RA*D*FZ
               ECAA_GAA  = CAA2*RA*(DDDG*FZ+D*FZ_Z*ZAA_GAA)*(ONE-fac)
     $ - CAA2*f_GAA*RA*D*FZ
               ECAA_TA   = CAA2*RA*(DDDT*FZ+D*FZ_Z*ZAA_TA)*(ONE-fac)
     $ - CAA2*f_TA*RA*D*FZ
               ECAA_LA   = CAA2*RA*D*FZ_Z*ZAA_LA*(ONE-fac)
     $ - CAA2*f_LA*RA*D*FZ
               ECAA_UUA  = -CAA2*f_UA*RA*D*FZ
               ECAA_RB   = -CAA2*f_RB*RA*D*FZ
               ECAA_GBB  = -CAA2*f_GBB*RA*D*FZ
               ECAA_TB   = -CAA2*f_TB*RA*D*FZ
               ECAA_LB   = -CAA2*f_LB*RA*D*FZ
               ECAA_UUB   = -CAA2*f_UB*RA*D*FZ
            endif
c               write(6,*)"ECAA_GBB",ECAA_GBB
c              ECAA_RA   = CAA2*(D*FZ+RA*DDDR*FZ+RA*D*FZ_Z*ZAA_RA)
c              ECAA_GAA  = CAA2*RA*(DDDG*FZ+D*FZ_Z*ZAA_GAA)
c              ECAA_TA   = CAA2*RA*(DDDT*FZ+D*FZ_Z*ZAA_TA)
c              ECAA_LA   = CAA2*RA*D*FZ_Z*ZAA_LA
c               write(6,*)"ECAA",ECAA
            END IF
         END IF

         ! EC_{beta,beta}
         UB       = 0.0D0
         UB_RB    = 0.0D0
         UB_GBB   = 0.0D0
         UB_TB    = 0.0D0
         UB_LB    = 0.0D0
         ECBB     = 0.0D0
         ECBB_RB  = 0.0D0
         ECBB_GBB = 0.0D0
         ECBB_TB  = 0.0D0
         ECBB_LB  = 0.0D0
         ECBB_UUB  = 0.0D0
         ECBB_UUA = 0.0D0
         IF (NB > 1) THEN
            IF (NDEN == 1) THEN
               UB       = UA
               UB_RB    = UA_RA
               UB_GBB   = UA_GAA
               UB_TB    = UA_TA
               UB_LB    = UA_LA
               ECBB     = ECAA
               ECBB_RB  = ECAA_RA
               ECBB_GBB = ECAA_GAA
               ECBB_TB  = ECAA_TA
               ECBB_LB  = ECAA_LA
               ECBB_UUB=ECAA_UUA
               ECBB_UUA=ECAA_UUB
            ELSE IF (RB > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RB,GBB,TB,LB,UB,UB_RB,UB_GBB,
     $                       UB_TB,UB_LB)
   
               ! now it's ZBB
               ZBB      = 0.0D0
               ZBB_RB   = 0.0D0
               ZBB_GBB  = 0.0D0
               ZBB_TB   = 0.0D0
               ZBB_LB   = 0.0D0
c               write(6,*)"UB",UB
               IF (ABS(UB) > THRESH) THEN
                  ONEDU   = 1.0D0/UB
                  ZBB     = CBB*2.0D0*ONEDU
                  ZBB_U   = -2.0D0*CBB/(UB*UB)
                  ZBB_RB  = ZBB_U*UB_RB
                  ZBB_GBB = ZBB_U*UB_GBB
                  ZBB_TB  = ZBB_U*UB_TB
                  ZBB_LB  = ZBB_U*UB_LB
               END IF
c               write(6,*)"ZBB",ZBB
   
               ! now it's function of ZBB
               FZ   = 0.0D0
               FZ_Z = 0.0D0
               IF (ABS(ZBB) > THRESH) THEN
                  ZBB2 = ZBB*ZBB
                  ZBB3 = ZBB2*ZBB
                  ZBB4 = ZBB3*ZBB
                  FZ0  = DLOG(1.0D0+ZBB/2.0D0)
                  FZ1  = 1.0D0-(2.0D0/ZBB)*FZ0
                  FZ   = ZBB4*FZ1
                  FZ_Z = 4.0D0*ZBB3*FZ1 + ZBB4*(2.0D0*FZ0/ZBB2 -
     $                   1.0D0/(ZBB*(1.0D0+ZBB/2.0D0)))
   
               END IF
c               write(6,*)"FZ for ZBB",FZ
   
               ! D and it's derivatives
               D    = TB - 0.25D0*GBB/RB
               DDDR = 0.25D0*GBB/(RB*RB)
               DDDG = -0.25D0/RB
               DDDT = 1.0D0
c              write(6,*)"D",D
   
               CALL becke05_A_sigma(iSpin,NB,NE,THRESH,
     $ VAL_P,VAL_Q,HIRWT,
     $ RA,GAA,TA,LA,UUA,RB,GBB,TB,LB,UUB,A1,
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
               ! now let's assemble all of pices together
                  IF(fac.le.ONE) then
               ECBB      = CBB2*RB*D*FZ*(ONE-fac)
               ECBB_RB   = CBB2*(D*FZ+RB*DDDR*FZ+RB*D*FZ_Z*ZBB_RB)
     $*(ONE-fac) - CBB2*RB*D*FZ*f_RB
               ECBB_GBB  = CBB2*RB*(DDDG*FZ+D*FZ_Z*ZBB_GBB)*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_GBB
               ECBB_TB   = CBB2*RB*(DDDT*FZ+D*FZ_Z*ZBB_TB)*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_TB
               ECBB_LB   = CBB2*RB*D*FZ_Z*ZBB_LB*(ONE-fac)
     $ -CBB2*RB*D*FZ*f_LB
               ECBB_UUB  = -CBB2*RB*D*FZ*f_UB
               ECBB_RA   = -CBB2*RB*D*FZ*f_RA
               ECBB_GAA  = -CBB2*RB*D*FZ*f_GAA
               ECBB_TA   = -CBB2*RB*D*FZ*f_TA
               ECBB_LA   = -CBB2*RB*D*FZ*f_LA
               ECBB_UUA  = -CBB2*RB*D*FZ*f_UA
             endif
            END IF
         END IF

         ! collect energy
         F(i) = F(i) + ECAA + ECBB

         ! collect all of terms for alpha derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAA_RA 
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAA_GAA 
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAA_TA 
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAA_LA 
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECBB_RB 
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECBB_GBB 
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECBB_TB 
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECBB_LB 

      END DO

      RETURN 
      END

      SUBROUTINE br94coor_op(INFOR,NDEN,NA,NB,NG,THRESH,
     $VAL_P,HIRWTS,RhoA,RhoB,
     $DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *   This is BR94 correlation functional which defined in:        *
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *   this is for the opposite spin component of the functional    *         
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
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
      INTEGER NA,NB,NE,iSpin
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB,VAL_P,VAL_Q
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS

      ! now let's define the parameter
      REAL*8 CAA,CBB,CAB
      REAL*8 CAA2,CBB2,CAB2

      ! variables for U and it's derivatives
      REAL*8 UA,UA_RA,UA_GAA,UA_TA,UA_LA    
      REAL*8 UB,UB_RB,UB_GBB,UB_TB,UB_LB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAB,ZAB_RA,ZAB_GAA,ZAB_TA,ZAB_LA    
      REAL*8 ZAB_RB,ZAB_GBB,ZAB_TB,ZAB_LB    
      REAL*8 ZAA_U,ZAB_U,ZBB_U,UUA,UUB
      REAL*8 FZ0,FZ1,FZ,FZ_Z
      REAL*8 ONEDU
      REAL*8 ZERO,ONE,TWO,THREE,FOUR,FIVE,F12,F14,F13,F16,F112

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA
      REAL*8 ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB,ECAB_UUA,ECAB_UUB  
      REAL*8 fac,f_RA,f_GAA,f_TA,f_LA,f_UA,f_RB,f_GBB,f_TB,f_LB,f_UB

C       
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       ! constants defined in the paper
       CAA   = -0.88D0
       CBB   = -0.88D0
       CAB   = -0.63D0
       CAA2  = -0.01D0
       CBB2  = -0.01D0
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
         UUA  = EXA(i)
         UUB  = EXB(i)
         HIRWT = HIRWTS(i)
c         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
         
         ! EC_{alpha,alpha}
         UA       = 0.0D0
         UA_RA    = 0.0D0
         UA_GAA   = 0.0D0
         UA_TA    = 0.0D0
         UA_LA    = 0.0D0
         IF (NA > 0) THEN
            IF (RA > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RA,GAA,TA,LA,UA,UA_RA,UA_GAA,
     $                       UA_TA,UA_LA)
c               write(6,*)"UA in ZAB", UA
            END IF
         END IF

         ! EC_{beta,beta}
         UB       = 0.0D0
         UB_RB    = 0.0D0
         UB_GBB   = 0.0D0
         UB_TB    = 0.0D0
         UB_LB    = 0.0D0
         IF (NB > 0) THEN
            IF (NDEN == 1) THEN
               UB       = UA
               UB_RB    = UA_RA
               UB_GBB   = UA_GAA
               UB_TB    = UA_TA
               UB_LB    = UA_LA
            ELSE IF (RB > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RB,GBB,TB,LB,UB,UB_RB,UB_GBB,
     $                       UB_TB,UB_LB)
c               write(6,*)"UB in ZAB", UB
            END IF
         END IF

         ! now it's ECAB
         ECAB     = 0.0D0
         ECAB_RA  = 0.0D0
         ECAB_GAA = 0.0D0
         ECAB_TA  = 0.0D0
         ECAB_LA  = 0.0D0
         ECAB_RB  = 0.0D0
         ECAB_GBB = 0.0D0
         ECAB_TB  = 0.0D0
         ECAB_LB  = 0.0D0
         ECAB_UUA = 0.0D0
         ECAB_UUB = 0.0D0
         IF (RA > THRESH .and. RB > THRESH .and. 
     $       NA > 0 .and. NB > 0) THEN

            ! now it's ZAB
            ZAB      = 0.0D0
            ZAB_RA   = 0.0D0
            ZAB_GAA  = 0.0D0
            ZAB_TA   = 0.0D0
            ZAB_LA   = 0.0D0
            ZAB_RB   = 0.0D0
            ZAB_GBB  = 0.0D0
            ZAB_TB   = 0.0D0
            ZAB_LB   = 0.0D0
            IF (ABS(UA) > THRESH) THEN
               ONEDU   = 1.0D0/UA
               ZAB     = CAB*ONEDU
               ZAB_U   = -CAB/(UA*UA)
               ZAB_RA  = ZAB_U*UA_RA
               ZAB_GAA = ZAB_U*UA_GAA
               ZAB_TA  = ZAB_U*UA_TA
               ZAB_LA  = ZAB_U*UA_LA
            END IF
            IF (NDEN == 1) THEN
               ZAB     = 2.0D0*ZAB 
               ZAB_RB  = ZAB_RA
               ZAB_GBB = ZAB_GAA
               ZAB_TB  = ZAB_TA
               ZAB_LB  = ZAB_LA
            ELSE IF (ABS(UB) > THRESH) THEN
               ONEDU   = 1.0D0/UB
               ZAB     = ZAB + CAB*ONEDU
               ZAB_U   = -CAB/(UB*UB)
               ZAB_RB  = ZAB_U*UB_RB
               ZAB_GBB = ZAB_U*UB_GBB
               ZAB_TB  = ZAB_U*UB_TB
               ZAB_LB  = ZAB_U*UB_LB
            END IF
c            write(6,*)"ZAB",ZAB

            ! now it's function of ZAB
            FZ   = 0.0D0
            FZ_Z = 0.0D0
            IF (ABS(ZAB) > THRESH) THEN
               FZ0  = DLOG(1.0D0+ZAB)/ZAB
               FZ   = ZAB*ZAB*(1.0D0-FZ0)
               FZ_Z = 2.0D0*ZAB*(1.0D0-FZ0) + ZAB*ZAB*(FZ0/ZAB -
     $                1.0D0/((1.0D0+ZAB)*ZAB))

            END IF
c            write(6,*)"FZ for ZAB",FZ
         CALL becke05_f(THRESH,VAL_P,HIRWT,RA,GAA,TA,LA,UUA,
     $ RB,GBB,TB,LB,UUB,fac,f_RA,f_GAA,f_TA,f_LA,f_UA,
     $ f_RB,f_GBB,f_TB,f_LB,f_UB)

            ! now let's assemble all of pices together
c           ECAB      = CAB2*RA*RB*FZ
c           ECAB_RA   = CAB2*RB*(FZ+RA*FZ_Z*ZAB_RA)
c           ECAB_RB   = CAB2*RA*(FZ+RB*FZ_Z*ZAB_RB)
c           ECAB_GAA  = CAB2*RA*RB*FZ_Z*ZAB_GAA
c           ECAB_GBB  = CAB2*RA*RB*FZ_Z*ZAB_GBB
c           ECAB_TA   = CAB2*RA*RB*FZ_Z*ZAB_TA
c           ECAB_TB   = CAB2*RA*RB*FZ_Z*ZAB_TB
c           ECAB_LA   = CAB2*RA*RB*FZ_Z*ZAB_LA
c           ECAB_LB   = CAB2*RA*RB*FZ_Z*ZAB_LB
        IF(fac.le.1.0D0) then
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
            ECAB_UUA   = - CAB2*RA*RB*FZ*f_UA 
            ECAB_UUB   = - CAB2*RA*RB*FZ*f_UB 
           endif
         END IF

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
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECAB_RB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECAB_GBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECAB_TB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECAB_LB

      END DO

      RETURN 
      END
