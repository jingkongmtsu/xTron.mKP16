C
C note by fenglai on Dec. 2015:
C
C we perform find difference comparison for KP14EC, it shows 
C that for some grid points (rho between 0.1 and 0.01) we have
C great difference between find difference and analytical ones.
C most of them are on beta parts.
C
C such comparison is with small p value in f factor calculating(p=0.01),
C we tested the B05 NDOP and B05 PAR they do not have such great
C difference. So the reason for these great difference should be 
C in the code of KP14C itself. Waiting for future solving.
C
C if you want to repeat the process, please use the following setting:
C 1 P in f factor calculating is set to 0.01 (original P value is not
C   stable for SCF convergence)
C 2 b = 1.355 alpha = 0.038 c_ndpar = 1.128
C please call the find difference fortran routine and inside call KP14EC
C
C Now. 28 2016: change the code for par_opt_kp14
C
C 1  set the parameter of 0.59 to 0.50, since 0.50 should be the correct
C    theory value
C
C 2  the CUTOFF1 value suggested by emil we still need to keep it. For 
C    deriving the Q_{ac} expression, it seems with smaller threshold value 
C    to determine the b*z the Q_{qc} will go wild and it will totally break 
C    the SCF result to make the energy goes to wild number. So we need to keep 
C    this cut off value to make the SCF stable.
C
C

      SUBROUTINE KP14EC(B05NDPARMethod,INFOR,NDEN,NG,NA,NB,b,alpha,
     $c_ndpar,c_ndpar_cap,THRESH,VAL_P,VAL_Q,
     $HIRWTS,RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,TF,D1F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the lambda dependent term correlation functional    *
c    *  this functional depends on three non-linear parameters,       *
c    *  namely b, alpha nad c_ndpar. See the paper for more details.  *         
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
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
      INTEGER NDEN,NG,I,J,method,B05NDPARMethod
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 b,alpha,PI,c_ndpar,c_ndpar_cap
      REAL*8 C
      REAL*8 THRESH
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 VAL_P,VAL_Q
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 TF(2*NG)
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

      ! static correlation for opposite spin
      ! it's derivatives
      REAL*8 usopp,usopp_RA,usopp_RB,usopp_GAA,usopp_GBB   
      REAL*8 usopp_TA,usopp_TB,usopp_LA,usopp_LB,usopp_UA,usopp_UB

      ! static correlation for parallel spin
      ! it's derivatives
      REAL*8 uspar,uspar_RA,uspar_RB,uspar_GAA,uspar_GBB   
      REAL*8 uspar_TA,uspar_TB,uspar_LA,uspar_LB,uspar_UA,uspar_UB

      ! static correlation for both of spin
      ! it's derivatives
      REAL*8 us,us_RA,us_RB,us_GAA,us_GBB   
      REAL*8 us_TA,us_TB,us_LA,us_LB,us_UA,us_UB

      ! static correlation for both of spin
      ! it's derivatives
      ! this will be formed by using c_ndpar
      ! so it's a new us, named as nus(new us)
      REAL*8 nus,nus_RA,nus_RB,nus_GAA,nus_GBB   
      REAL*8 nus_TA,nus_TB,nus_LA,nus_LB,nus_UA,nus_UB

      ! dynamic correlation for opposite spin
      ! it's derivatives
      REAL*8 udopp,udopp_RA,udopp_RB,udopp_GAA,udopp_GBB   
      REAL*8 udopp_TA,udopp_TB,udopp_LA,udopp_LB,udopp_UA,udopp_UB

      ! dynamic correlation for parallel spin
      ! it's derivatives
      REAL*8 udpar,udpar_RA,udpar_RB,udpar_GAA,udpar_GBB   
      REAL*8 udpar_TA,udpar_TB,udpar_LA,udpar_LB,udpar_UA,udpar_UB

      ! dynamic correlation for both of spin
      ! it's derivatives
      REAL*8 ud,ud_RA,ud_RB,ud_GAA,ud_GBB   
      REAL*8 ud_TA,ud_TB,ud_LA,ud_LB,ud_UA,ud_UB
      REAL*8 ud2

      ! z variable and it's derivatives
      REAL*8 z,z_RA,z_RB,z_GAA,z_GBB   
      REAL*8 z_TA,z_TB,z_LA,z_LB,z_UA,z_UB

      ! z0 variable and it's derivatives
      REAL*8 z0,z0_RA,z0_RB,z0_GAA,z0_GBB   
      REAL*8 z0_TA,z0_TB,z0_LA,z0_LB,z0_UA,z0_UB
      REAL*8 dzdz0

      !QAC
      REAL*8 LIMIT_QAC
      REAL*8 QAC, QAC_Z, bz, ebz, lebz, L2
      REAL*8 T0, T1, T2, T3

      ! DA and DB
      REAL*8 R89, R19, R179, R23, R53, D
      REAL*8 DA, DA_RA, DA_GAA, DA_TA
      REAL*8 DB, DB_RB, DB_GBB, DB_TB

      ! term describe spin error 
      ! it's modified based on us
      REAL*8 usm,usm_RA,usm_RB,usm_GAA,usm_GBB   
      REAL*8 usm_TA,usm_TB,usm_LA,usm_LB,usm_UA,usm_UB
      REAL*8 usm_opp, usm_par

      ! terms for the exp(-alpha*z2)
      REAL*8 eaz2, z3, eaz2_z

      ! tmp vactor for getting the derivatives
      ! totally it's 11 variables including GAB, see above
      ! so we just use the N_FUNC_DERIV_1
      REAL*8 TMP_DRA(3)
      REAL*8 TMP_DRB(3)
      REAL*8 TMP_D1F(N_FUNC_DERIV_1)

      ! constant
      REAL*8 ZERO
      REAL*8 ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE
      REAL*8 F12,F14,F13,F53,F23,F112

      ! the cut off 1 is related to Qac form evaluation
      ! if z is smaller than this value, through careful
      ! numerical test, we see that the expression
      ! of Qac goes to -1/2, and derivatives becomes 0
      !
      ! this can also be tested with maple program  
      REAL*8 CUTOFF1
      parameter(CUTOFF1 =  1.0D-8)

      ! set up some big number
      ! this is z's default value
      REAL*8 bignum
      parameter(bignum =  1.0D+8)

      ! set up the delta value for smoothing the z
      REAL*8 DELTA, INFINITY 
c      parameter(DELTA  =  0.0000005)
c      parameter(DELTA  =  0.005)
c      parameter(DELTA  =  0.05)
      parameter(DELTA=0.0000000001d0)  !
      parameter(INFINITY=1d6)  !

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
      EIGHT= 8.0D0
      NINE = 9.0D0
      F12  = 0.5D0
      F14  = 0.25D0
      F13  = 0.33333333333333333333333D0
      F53  = 1.66666666666666666666666D0
      F23  = 0.66666666666666666666666D0
      F112 = 0.08333333333333333333333D0
      PI   = 4.D0*DATAN(1.D0)

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! loop over NG
      DO I = 1,NG

c         write(6,*)"     "
c         write(6,*)"----------------------"
c         write(6,*)"for grid: ", i
          
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

         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! tmp dra and drb
         TMP_DRA(1) = DRA(i,1)
         TMP_DRA(2) = DRA(i,2)
         TMP_DRA(3) = DRA(i,3)
         TMP_DRB(1) = DRB(i,1)
         TMP_DRB(2) = DRB(i,2)
         TMP_DRB(3) = DRB(i,3)

         ! debug
c         write(6,*)"our RhoA: ",RA
c         write(6,*)"our GRhoA: ",GAA 
c         write(6,*)"our TA: ",TA
c         write(6,*)"our LA: ",LA 
c         write(6,*)"our UA: ",UA
c         write(6,*)"our RhoB: ",RB
c         write(6,*)"our GRhoB: ",GBB 
c         write(6,*)"our TB: ",TB
c         write(6,*)"our LB: ",LB 
c         write(6,*)"our UB: ",UB 
c         write(6,*)" B05NDPARMethod in kp14ec  ", B05NDPARMethod 

         ! get the b05 non-dynamic energy density for parallel
         uspar = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         if ( B05NDPARMethod .eq. 0 ) then
           CALL Becke05_NDPAR(INFOR,NDEN,1,NA,NB,THRESH,
     $                        VAL_P,VAL_Q,HIRWT,
     $                        RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $                        uspar,TMP_D1F)
         else if ( B05NDPARMethod .eq. 1 ) then
           CALL Becke05_NDPAR1(INFOR,NDEN,1,NA,NB,THRESH,c_ndpar_cap,
     $                        VAL_P,VAL_Q,HIRWT,
     $                        RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $                        uspar,TMP_D1F)
         else if ( B05NDPARMethod .eq. 2 ) then
           CALL Becke05_NDPAR2(INFOR,NDEN,1,NA,NB,THRESH,c_ndpar_cap,
     $                        VAL_P,VAL_Q,HIRWT,
     $                        RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $                        uspar,TMP_D1F)
         else
            write(*,*) "wrong choice of method for NDpar"
         endif
       
         ID_RA_POS  = D1VARS(ID_RA)
         uspar_RA   = TMP_D1F(ID_RA_POS)
         ID_GAA_POS = D1VARS(ID_GAA)
         uspar_GAA  = TMP_D1F(ID_GAA_POS)
         ID_TA_POS  = D1VARS(ID_TA)
         uspar_TA   = TMP_D1F(ID_TA_POS)
         ID_LA_POS  = D1VARS(ID_LA)
         uspar_LA   = TMP_D1F(ID_LA_POS)
         ID_EXA_POS = D1VARS(ID_EXA)
         uspar_UA   = TMP_D1F(ID_EXA_POS)
         ID_RB_POS  = D1VARS(ID_RB)
         uspar_RB   = TMP_D1F(ID_RB_POS)
         ID_GBB_POS = D1VARS(ID_GBB)
         uspar_GBB  = TMP_D1F(ID_GBB_POS)
         ID_TB_POS  = D1VARS(ID_TB)
         uspar_TB   = TMP_D1F(ID_TB_POS)
         ID_LB_POS  = D1VARS(ID_LB)
         uspar_LB   = TMP_D1F(ID_LB_POS)
         ID_EXB_POS = D1VARS(ID_EXB)
         uspar_UB   = TMP_D1F(ID_EXB_POS)

         ! get the b05 non-dynamic energy density for opposite
         usopp = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,VAL_P,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,usopp,TMP_D1F)
         usopp_RA   = TMP_D1F(ID_RA_POS)
         usopp_GAA  = TMP_D1F(ID_GAA_POS)
         usopp_TA   = TMP_D1F(ID_TA_POS)
         usopp_LA   = TMP_D1F(ID_LA_POS)
         usopp_UA   = TMP_D1F(ID_EXA_POS)
         usopp_RB   = TMP_D1F(ID_RB_POS)
         usopp_GBB  = TMP_D1F(ID_GBB_POS)
         usopp_TB   = TMP_D1F(ID_TB_POS)
         usopp_LB   = TMP_D1F(ID_LB_POS)
         usopp_UB   = TMP_D1F(ID_EXB_POS)

         ! now get the sum of them
         us = usopp+uspar
         us_RA   = uspar_RA  + usopp_RA 
         us_GAA  = uspar_GAA + usopp_GAA
         us_TA   = uspar_TA  + usopp_TA 
         us_LA   = uspar_LA  + usopp_LA 
         us_UA   = uspar_UA  + usopp_UA 
         us_RB   = uspar_RB  + usopp_RB 
         us_GBB  = uspar_GBB + usopp_GBB
         us_TB   = uspar_TB  + usopp_TB 
         us_LB   = uspar_LB  + usopp_LB 
         us_UB   = uspar_UB  + usopp_UB 

         ! derive the new us based on
         ! parameter c_ndpar
         nus = usopp+c_ndpar*uspar
         nus_RA   = c_ndpar*uspar_RA  + usopp_RA 
         nus_GAA  = c_ndpar*uspar_GAA + usopp_GAA
         nus_TA   = c_ndpar*uspar_TA  + usopp_TA 
         nus_LA   = c_ndpar*uspar_LA  + usopp_LA 
         nus_UA   = c_ndpar*uspar_UA  + usopp_UA 
         nus_RB   = c_ndpar*uspar_RB  + usopp_RB 
         nus_GBB  = c_ndpar*uspar_GBB + usopp_GBB
         nus_TB   = c_ndpar*uspar_TB  + usopp_TB 
         nus_LB   = c_ndpar*uspar_LB  + usopp_LB 
         nus_UB   = c_ndpar*uspar_UB  + usopp_UB 
c         write(6,*)"our us: ",us
c         write(6,*)"our nus: ",nus
c         write(6,*)"our us, opp: ",usopp
c         write(6,*)"our us, par: ",uspar
c         write(6,*)"nus deriv on RA: ", nus_RA
c         write(6,*)"nus deriv on RB: ", nus_RB
c         write(6,*)"nus deriv on GAA: ",nus_GAA
c         write(6,*)"nus deriv on GBB: ",nus_GBB
c         write(6,*)"nus deriv on TA: ", nus_TA
c         write(6,*)"nus deriv on TB: ", nus_TB
c         write(6,*)"nus deriv on LA: ", nus_LA
c         write(6,*)"nus deriv on LB: ", nus_LB
c         write(6,*)"nus deriv on EXA: ",nus_UA
c         write(6,*)"nus deriv on EXB: ",nus_UB
c         write(6,*)"us deriv on RA: ", us_RA
c         write(6,*)"us deriv on RB: ", us_RB
c         write(6,*)"us deriv on GAA: ",us_GAA
c         write(6,*)"us deriv on GBB: ",us_GBB
c         write(6,*)"us deriv on TA: ", us_TA
c         write(6,*)"us deriv on TB: ", us_TB
c         write(6,*)"us deriv on LA: ", us_LA
c         write(6,*)"us deriv on LB: ", us_LB
c         write(6,*)"us deriv on EXA: ",us_UA
c         write(6,*)"us deriv on EXB: ",us_UB

         ! get the dynamic energy density for parallel
         method = 2
         udpar = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL b13coor_par(method,INFOR,NDEN,NA,NB,1,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,udpar,TMP_D1F)
         udpar_RA   = TMP_D1F(ID_RA_POS)
         udpar_GAA  = TMP_D1F(ID_GAA_POS)
         udpar_TA   = TMP_D1F(ID_TA_POS)
         udpar_LA   = TMP_D1F(ID_LA_POS)
         udpar_UA   = TMP_D1F(ID_EXA_POS)
         udpar_RB   = TMP_D1F(ID_RB_POS)
         udpar_GBB  = TMP_D1F(ID_GBB_POS)
         udpar_TB   = TMP_D1F(ID_TB_POS)
         udpar_LB   = TMP_D1F(ID_LB_POS)
         udpar_UB   = TMP_D1F(ID_EXB_POS)

         ! get the dynamic energy density for opposite
         udopp = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL b13coor_opp(method,INFOR,NDEN,NA,NB,1,THRESH,
     $VAL_P,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,udopp,TMP_D1F)
         udopp_RA   = TMP_D1F(ID_RA_POS)
         udopp_GAA  = TMP_D1F(ID_GAA_POS)
         udopp_TA   = TMP_D1F(ID_TA_POS)
         udopp_LA   = TMP_D1F(ID_LA_POS)
         udopp_UA   = TMP_D1F(ID_EXA_POS)
         udopp_RB   = TMP_D1F(ID_RB_POS)
         udopp_GBB  = TMP_D1F(ID_GBB_POS)
         udopp_TB   = TMP_D1F(ID_TB_POS)
         udopp_LB   = TMP_D1F(ID_LB_POS)
         udopp_UB   = TMP_D1F(ID_EXB_POS)

         ! now get the sum of them
         ud      = udopp+udpar
         ud_RA   = udpar_RA  + udopp_RA 
         ud_GAA  = udpar_GAA + udopp_GAA
         ud_TA   = udpar_TA  + udopp_TA 
         ud_LA   = udpar_LA  + udopp_LA 
         ud_UA   = udpar_UA  + udopp_UA 
         ud_RB   = udpar_RB  + udopp_RB 
         ud_GBB  = udpar_GBB + udopp_GBB
         ud_TB   = udpar_TB  + udopp_TB 
         ud_LB   = udpar_LB  + udopp_LB 
         ud_UB   = udpar_UB  + udopp_UB 
c         write(6,*)"our ud: ",ud
c         write(6,*)"ud par deriv on RA: ", udpar_RA
c         write(6,*)"ud par deriv on RB: ", udpar_RB
c         write(6,*)"ud par deriv on GAA: ",udpar_GAA
c         write(6,*)"ud par deriv on GBB: ",udpar_GBB
c         write(6,*)"ud par deriv on TA: ", udpar_TA
c         write(6,*)"ud par deriv on TB: ", udpar_TB
c         write(6,*)"ud par deriv on LA: ", udpar_LA
c         write(6,*)"ud par deriv on LB: ", udpar_LB
c         write(6,*)"ud par deriv on EXA: ",udpar_UA
c         write(6,*)"ud par deriv on EXB: ",udpar_UB
c
c         write(6,*)"ud opp deriv on RA: ", udopp_RA
c         write(6,*)"ud opp deriv on RB: ", udopp_RB
c         write(6,*)"ud opp deriv on GAA: ",udopp_GAA
c         write(6,*)"ud opp deriv on GBB: ",udopp_GBB
c         write(6,*)"ud opp deriv on TA: ", udopp_TA
c         write(6,*)"ud opp deriv on TB: ", udopp_TB
c         write(6,*)"ud opp deriv on LA: ", udopp_LA
c         write(6,*)"ud opp deriv on LB: ", udopp_LB
c         write(6,*)"ud opp deriv on EXA: ",udopp_UA
c         write(6,*)"ud opp deriv on EXB: ",udopp_UB

         ! let's calculate original z 
         z0    = ZERO 
         z0_RA = ZERO 
         z0_GAA= ZERO 
         z0_TA = ZERO 
         z0_LA = ZERO 
         z0_UA = ZERO 
         z0_RB = ZERO 
         z0_GBB= ZERO 
         z0_TB = ZERO 
         z0_LB = ZERO 
         z0_UB = ZERO 
         IF (DABS(ud)>THRESH) THEN
            z0    = us/ud
            z0_RA = us_RA/ud  
            z0_GAA= us_GAA/ud 
            z0_TA = us_TA/ud  
            z0_LA = us_LA/ud  
            z0_UA = us_UA/ud  
            z0_RB = us_RB/ud  
            z0_GBB= us_GBB/ud 
            z0_TB = us_TB/ud  
            z0_LB = us_LB/ud  
            z0_UB = us_UB/ud  
            ud2   = ud*ud
            IF (ud2>THRESH) THEN
               z0_RA = z0_RA  - us*ud_RA/ud2
               z0_GAA= z0_GAA - us*ud_GAA/ud2
               z0_TA = z0_TA  - us*ud_TA/ud2
               z0_LA = z0_LA  - us*ud_LA/ud2
               z0_UA = z0_UA  - us*ud_UA/ud2
               z0_RB = z0_RB  - us*ud_RB/ud2
               z0_GBB= z0_GBB - us*ud_GBB/ud2
               z0_TB = z0_TB  - us*ud_TB/ud2
               z0_LB = z0_LB  - us*ud_LB/ud2
               z0_UB = z0_UB  - us*ud_UB/ud2
            END IF
         ELSE
            !This covers the case when there is no dynamic correlation,
            !typically only one orbital with less than an electron it
            !(not sure why ud is not zero for one electron).  Note, this
            !is regardless whether us is zero or not.  The latter could
            !happen, e.g. H(1/3alpha,1/3beta) for points that density
            !is zero. 
            z0 = INFINITY  
         END IF

         ! now let's smooth the z
         z    = ZERO
         z_RA = ZERO 
         z_GAA= ZERO 
         z_TA = ZERO 
         z_LA = ZERO 
         z_UA = ZERO 
         z_RB = ZERO 
         z_GBB= ZERO 
         z_TB = ZERO 
         z_LB = ZERO 
         z_UB = ZERO 
         IF (z0>=DELTA) THEN
            z    = z0    
            z_RA = z0_RA 
            z_GAA= z0_GAA
            z_TA = z0_TA 
            z_LA = z0_LA 
            z_UA = z0_UA 
            z_RB = z0_RB 
            z_GBB= z0_GBB
            z_TB = z0_TB 
            z_LB = z0_LB 
            z_UB = z0_UB 
         ELSE IF(DABS(z0)<DELTA)THEN
            z = (ONE/(FOUR*DELTA))*(z0+DELTA)*(z0+DELTA)
            dzdz0 = (ONE/(TWO*DELTA))*(z0+DELTA)
            z_RA  = dzdz0*z0_RA 
            z_GAA = dzdz0*z0_GAA
            z_TA  = dzdz0*z0_TA 
            z_LA  = dzdz0*z0_LA 
            z_UA  = dzdz0*z0_UA 
            z_RB  = dzdz0*z0_RB 
            z_GBB = dzdz0*z0_GBB
            z_TB  = dzdz0*z0_TB 
            z_LB  = dzdz0*z0_LB 
            z_UB  = dzdz0*z0_UB 
c            write(6,*)"ud or us cause abnormal z", us, ud
         END IF

         ! now it's the qac
         ! we need to see more carefully about the z
         ! the case we need to care about, is that when z is 
         ! very small, or z is some positive large number
         ! else we will do it in normal way
         bz = b*z
         QAC    = ZERO
         QAC_Z  = ZERO
         LIMIT_QAC = DLOG(ONE/THRESH)
         L2     = DLOG(TWO)
         IF (DABS(bz)<CUTOFF1) THEN
            QAC   = 0.50D0
            QAC_Z = ZERO
         ELSE IF (bz>LIMIT_QAC) THEN
            QAC   = -ONE + TWO*(bz-L2)/bz
            QAC_Z = TWO*L2*b/(bz*bz)
         ELSE
            ebz  = dexp(bz)
            lebz = DLOG(ONE+ebz)
            QAC  = (ONE+ebz)/(ONE-ebz)*(ONE-TWO*(lebz-L2)/bz)

            ! now it's derivatives
            T0    = ONE-TWO*(lebz-L2)/bz
            T1    = (ONE+ebz)/(ONE-ebz)
            T2    = -ebz/(ONE-ebz)
            QAC_Z = ebz*T0/(ONE-ebz)
            QAC_Z = QAC_Z - T1*T2*T0
            QAC_Z = QAC_Z + T1*(-TWO*ebz/((ONE+ebz)*bz)+
     $ TWO*(lebz-L2)/(bz*bz))
            QAC_Z = QAC_Z*b
         END IF
c        write(6,*)"our QAC_Z: ",QAC_Z

         ! now let's do DA 
         !
         ! let's remind ourself that this is related to
         ! the self spin error. It appears only NA > 1 
         !
         ! the DA calculation is in some strange way
         ! however, we found that if the T0 is some
         ! small negative number, for directly computing
         ! (D/R53)^(1/3) it some time returns NAN number 
         !
         ! be careful about the derivatives!
         ! now we add in protection like IF (DABS(R179*T0)>THRESH) THEN
         ! we found that if we do not have it for derivatives,
         ! it may go to NAN for very small input DFT variables
         !
         DA     = ZERO
         DA_RA  = ZERO
         DA_GAA = ZERO
         DA_TA  = ZERO
         IF (RA>THRESH .and. NA>1) THEN
            R23   = RA**F23
            R53   = RA*R23
            D     = TA-F14*GAA/RA
            T0    = D/R53
            DA    = DSIGN(DABS(T0)**F13,T0)

            ! now let's do derivative
            ! on rho
            T1     = TA*RA-F14*GAA
            T2     = T1*T1
            T0     = T2**F13
            R89    = RA**(EIGHT/NINE)
            R179   = RA*R89
            IF (DABS(T0)>THRESH) THEN
               DA_RA  = F13*(F14*GAA/R179-F53*D/R89)/T0
            END IF
c            write(6,*) "in DA_RA, R179*T0", R179*T0
c            IF (DABS(R179*T0)>THRESH) THEN
c               DA_RA = DA_RA + F13*(F14*GAA/R179)/T0
c            END IF
c            write(6,*) "in DA_RA, R89*T0", R89*T0
c            IF (DABS(R89*T0)>THRESH) THEN
c               DA_RA = DA_RA - F13*(F53*D/R89)/T0
c            END IF

            ! on tau and G
            R19    = RA/R89
            IF (DABS(T0)>THRESH) THEN
               DA_TA  = F13*R19/T0
            END IF
            IF (DABS(R89*T0)>THRESH) THEN
               DA_GAA = -F112/(T0*R89)
            END IF
         END IF
c         write(6,*)"our DA: ",DA
c         write(6,*)"our DA_RA: ",DA_RA
c         write(6,*)"our DA_GAA: ",DA_GAA
c         write(6,*)"our DA_TA: ",DA_TA

         ! now let's do DB 
         ! see comments for DA
         DB     = ZERO
         DB_RB  = ZERO
         DB_GBB = ZERO
         DB_TB  = ZERO
         IF (RB>THRESH .and. NB>1) THEN
            R23   = RB**F23
            R53   = RB*R23
            D     = TB-F14*GBB/RB
            T0    = D/R53
            DB    = DSIGN(DABS(T0)**F13,T0)

            ! now let's do derivative
            ! on rho
            T1     = TB*RB-F14*GBB
            T2     = T1*T1
            T0     = T2**F13
            R89    = RB**(EIGHT/NINE)
            R179   = RB*R89
            IF (DABS(T0)>THRESH) THEN
               DB_RB  = F13*(F14*GBB/R179-F53*D/R89)/T0
            END IF
c            IF (DABS(R179*T0)>THRESH) THEN
c               DB_RB = DB_RB + F13*(F14*GBB/R179)/T0
c            END IF
c            IF (DABS(R89*T0)>THRESH) THEN
c               DB_RB = DB_RB - F13*(F53*D/R89)/T0
c            END IF

            ! on tau and G
            R19    = RB/R89
            IF (DABS(T0)>THRESH) THEN
               DB_TB  = F13*R19/T0
            END IF
            IF (DABS(R89*T0)>THRESH) THEN
               DB_GBB = -F112/(T0*R89)
            END IF
         END IF
c         write(6,*)"our DB: ",DB
c         write(6,*)"our DB_RB: ",DB_RB
c         write(6,*)"our DB_GBB: ",DB_GBB
c         write(6,*)"our DB_TB: ",DB_TB

         ! now it's another term exp(-alpha*z2)
         eaz2 = ZERO
         IF (DABS(z)>THRESH) THEN
            eaz2 = DEXP(-alpha/(z*z))
         END IF
         z3 = z*z*z
         eaz2_z = ZERO
         IF (DABS(z3)>THRESH) THEN
            eaz2_z = eaz2*TWO*alpha/z3
         END IF
c         write(6,*)"our eaz2: ",eaz2
c         write(6,*)"our eaz2_z: ",eaz2_z

         ! now let's comput ethe spin error term
         C   = F12*DSQRT(alpha/PI)
         T0  = DA+DB
         usm = nus+C*nus*T0*eaz2 
         usm_opp = usopp*(ONE+C*T0*eaz2) 
         usm_par = uspar*(ONE+C*T0*eaz2) 
c         write(6,*)"DFF", T0
c         write(6,*)"S1",  DSQRT(alpha/PI)*eaz2
c         write(6,*)"alpha", alpha
c         write(6,*)"eshift", C*T0*eaz2
         usm_RA  = nus_RA   + C*nus_RA*T0*eaz2
     $ + C*nus*(DA_RA*eaz2  + T0*eaz2_z*z_RA)
         usm_RB  = nus_RB   + C*nus_RB*T0*eaz2
     $ + C*nus*(DB_RB*eaz2  + T0*eaz2_z*z_RB)
         usm_GAA = nus_GAA  + C*nus_GAA*T0*eaz2
     $ + C*nus*(DA_GAA*eaz2 + T0*eaz2_z*z_GAA)
         usm_GBB = nus_GBB  + C*nus_GBB*T0*eaz2
     $ + C*nus*(DB_GBB*eaz2 + T0*eaz2_z*z_GBB)
         usm_TA  = nus_TA   + C*nus_TA*T0*eaz2
     $ + C*nus*(DA_TA*eaz2  + T0*eaz2_z*z_TA)
         usm_TB  = nus_TB   + C*nus_TB*T0*eaz2
     $ + C*nus*(DB_TB*eaz2  + T0*eaz2_z*z_TB)
         usm_LA  = nus_LA   + C*nus_LA*T0*eaz2
     $ + C*nus*T0*eaz2_z*z_LA
         usm_LB  = nus_LB   + C*nus_LB*T0*eaz2
     $ + C*nus*T0*eaz2_z*z_LB
         usm_UA  = nus_UA   + C*nus_UA*T0*eaz2
     $ + C*nus*T0*eaz2_z*z_UA
         usm_UB  = nus_UB   + C*nus_UB*T0*eaz2
     $ + C*nus*T0*eaz2_z*z_UB
c         write(6,*)"usm deriv on RA: ", usm_RA
c         write(6,*)"usm deriv on RB: ", usm_RB
c         write(6,*)"usm deriv on GAA: ",usm_GAA
c         write(6,*)"usm deriv on GBB: ",usm_GBB
c         write(6,*)"usm deriv on TA: ", usm_TA
c         write(6,*)"usm deriv on TB: ", usm_TB
c         write(6,*)"usm deriv on LA: ", usm_LA
c         write(6,*)"usm deriv on LB: ", usm_LB
c         write(6,*)"usm deriv on EXA: ",usm_UA
c         write(6,*)"usm deriv on EXB: ",usm_UB
c
c         write(6,*)"d q dz: ", QAC_Z
c         write(6,*)"z deriv on RA: ", z_RA
c         write(6,*)"z deriv on RB: ", z_RB
c         write(6,*)"z deriv on GAA: ",z_GAA
c         write(6,*)"z deriv on GBB: ",z_GBB
c         write(6,*)"z deriv on TA: ", z_TA
c         write(6,*)"z deriv on TB: ", z_TB
c         write(6,*)"z deriv on LA: ", z_LA
c         write(6,*)"z deriv on LB: ", z_LB
c         write(6,*)"z deriv on EXA: ",z_UA
c         write(6,*)"z deriv on EXB: ",z_UB
c         write(6,*)"qac deriv on RA: ", QAC_Z*z_RA
c         write(6,*)"qac deriv on RB: ", QAC_Z*z_RB
c         write(6,*)"qac deriv on GAA: ",QAC_Z*z_GAA
c         write(6,*)"qac deriv on GBB: ",QAC_Z*z_GBB
c         write(6,*)"qac deriv on TA: ", QAC_Z*z_TA
c         write(6,*)"qac deriv on TB: ", QAC_Z*z_TB
c         write(6,*)"qac deriv on LA: ", QAC_Z*z_LA
c         write(6,*)"qac deriv on LB: ", QAC_Z*z_LB
c         write(6,*)"qac deriv on EXA: ",QAC_Z*z_UA
c         write(6,*)"qac deriv on EXB: ",QAC_Z*z_UB

         ! finally let's bring all of term together
         F(i) = F(i) + QAC*usm
         TF(i)    = TF(i) + QAC*usm_opp
         TF(i+ng) = TF(i+ng) + QAC*usm_par
c         write(6,*)"end0", QAC*nus
c        write(6,*)"KP14C energy", QAC*usm

         ! collect all of terms for derivatives
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS)  + QAC_Z*z_RA*usm 
     $ + QAC*usm_RA
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + QAC_Z*z_GAA*usm
     $ + QAC*usm_GAA
         ID_GAB_POS         = D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = ZERO
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS)  + QAC_Z*z_TA*usm
     $ + QAC*usm_TA
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS)  + QAC_Z*z_LA*usm
     $ + QAC*usm_LA
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + QAC_Z*z_UA*usm
     $ + QAC*usm_UA
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS)  + QAC_Z*z_RB*usm
     $ + QAC*usm_RB
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + QAC_Z*z_GBB*usm
     $ + QAC*usm_GBB
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS)  + QAC_Z*z_TB*usm
     $ + QAC*usm_TB
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS)  + QAC_Z*z_LB*usm
     $ + QAC*usm_LB
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + QAC_Z*z_UB*usm
     $ + QAC*usm_UB
c         write(6,*)"our deriv on RA: ",QAC_Z*z_RA*usm+QAC*usm_RA
c         write(6,*)"our deriv on RB: ",QAC_Z*z_RB*usm+QAC*usm_RB
c         write(6,*)"our deriv on GAA: ",QAC_Z*z_GAA*usm+QAC*usm_GAA
c         write(6,*)"our deriv on GBB: ",QAC_Z*z_GBB*usm+QAC*usm_GBB
c         write(6,*)"our deriv on TA: ",QAC_Z*z_TA*usm+QAC*usm_TA
c         write(6,*)"our deriv on TB: ",QAC_Z*z_TB*usm+QAC*usm_TB
c         write(6,*)"our deriv on LA: ",QAC_Z*z_LA*usm+QAC*usm_LA
c         write(6,*)"our deriv on LB: ",QAC_Z*z_LB*usm+QAC*usm_LB
c         write(6,*)"our deriv on EXA: ",QAC_Z*z_UA*usm+QAC*usm_UA
c         write(6,*)"our deriv on EXB: ",QAC_Z*z_UB*usm+QAC*usm_UB

      END DO
      RETURN
      END
