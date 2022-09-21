
      SUBROUTINE B13STRONG_AC2(INFOR,NDEN,NG,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the "strong" component for B13 Strong functional.   *
c    *  see the equations in section of VII in the following paper:   *
c    *  "Density functionals for static, dynamical, and strong        *
c    *   correlation"                                                 *
c    *   A. D. Becke  J. Chem. Phys.  138  074109                     *  
c    *                                                                *
c    *  this is the term with x^2*u.                                  *
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
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG
      INTEGER NA,NB
      INTEGER method
      INTEGER I,J
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 VAL_P,VAL_Q
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
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

      ! static correlation for opposite spin
      ! it's derivatives
      REAL*8 usopp,usopp_RA,usopp_RB,usopp_GAA,usopp_GBB   
      REAL*8 usopp_TA,usopp_TB,usopp_LA,usopp_LB,usopp_UA,usopp_UB

      ! static correlation for parallel spin
      ! it's derivatives
      REAL*8 uspar,uspar_RA,uspar_RB,uspar_GAA,uspar_GBB   
      REAL*8 uspar_TA,uspar_TB,uspar_LA,uspar_LB,uspar_UA,uspar_UB

      ! dynamic correlation for opposite spin
      ! it's derivatives
      REAL*8 udopp,udopp_RA,udopp_RB,udopp_GAA,udopp_GBB   
      REAL*8 udopp_TA,udopp_TB,udopp_LA,udopp_LB,udopp_UA,udopp_UB

      ! dynamic correlation for parallel spin
      ! it's derivatives
      REAL*8 udpar,udpar_RA,udpar_RB,udpar_GAA,udpar_GBB   
      REAL*8 udpar_TA,udpar_TB,udpar_LA,udpar_LB,udpar_UA,udpar_UB

      ! static correlation for both of spin
      ! it's derivatives
      REAL*8 us,us_RA,us_RB,us_GAA,us_GBB   
      REAL*8 us_TA,us_TB,us_LA,us_LB,us_UA,us_UB
      REAL*8 us2

      ! dynamic correlation for both of spin
      ! it's derivatives
      REAL*8 ud,ud_RA,ud_RB,ud_GAA,ud_GBB   
      REAL*8 ud_TA,ud_TB,ud_LA,ud_LB,ud_UA,ud_UB

      ! static correlation for uc in equation 52
      ! it's derivatives
      REAL*8 uc,uc_RA,uc_RB,uc_GAA,uc_GBB   
      REAL*8 uc_TA,uc_TB,uc_LA,uc_LB,uc_UA,uc_UB
      REAL*8 uc2

      ! the x defined in equation 52 
      REAL*8 x,x_RA,x_RB,x_GAA,x_GBB   
      REAL*8 x_TA,x_TB,x_LA,x_LB,x_UA,x_UB
      REAL*8 x2

      ! the x defined in equation 52, this is x^2*uc
      ! defined in equation 58
      REAL*8 x2uc,x2uc_RA,x2uc_RB,x2uc_GAA,x2uc_GBB   
      REAL*8 x2uc_TA,x2uc_TB,x2uc_LA,x2uc_LB,x2uc_UA,x2uc_UB

      ! tmp vactor for getting the derivatives
      ! totally it's 11 variables including GAB, see above
      ! so we just use the N_FUNC_DERIV_1
      REAL*8 TMP_DRA(3)
      REAL*8 TMP_DRB(3)
      REAL*8 TMP_D1F(N_FUNC_DERIV_1)

      ! constant
      REAL*8 ZERO
      REAL*8 TWO

C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! constants
      ZERO = 0.0D0
      TWO  = 2.0D0

      ! firstly initilize variable position information
      ! and we initialize all of position infor for deriv
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_LA_POS  = D1VARS(ID_LA)
      ID_EXA_POS = D1VARS(ID_EXA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TB_POS  = D1VARS(ID_TB)
      ID_LB_POS  = D1VARS(ID_LB)
      ID_EXB_POS = D1VARS(ID_EXB)
      ID_GAB_POS = D1VARS(ID_GAB)
      
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

         ! get the b05 non-dynamic energy density for parallel
         uspar = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL Becke05_NDPAR(INFOR,NDEN,1,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,uspar,TMP_D1F)
         uspar_RA   = TMP_D1F(ID_RA_POS)
         uspar_GAA  = TMP_D1F(ID_GAA_POS)
         uspar_TA   = TMP_D1F(ID_TA_POS)
         uspar_LA   = TMP_D1F(ID_LA_POS)
         uspar_UA   = TMP_D1F(ID_EXA_POS)
         uspar_RB   = TMP_D1F(ID_RB_POS)
         uspar_GBB  = TMP_D1F(ID_GBB_POS)
         uspar_TB   = TMP_D1F(ID_TB_POS)
         uspar_LB   = TMP_D1F(ID_LB_POS)
         uspar_UB   = TMP_D1F(ID_EXB_POS)
c         write(6,*)"uspar", uspar 

         ! get the b05 non-dynamic energy density for opposite
         usopp = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,
     $VAL_P,HIRWT,
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
c         write(6,*)"usopp", usopp 

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
c         write(6,*)"our us: ",us
c         write(6,*)"our nus: ",nus
c         write(6,*)"our us, opp: ",usopp
c         write(6,*)"our us, par: ",uspar
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
c         write(6,*)"udpar", udpar

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
c         write(6,*)"udopp", udopp

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

         ! let's calculate uc
         uc    = us + ud
         uc_RA = us_RA  + ud_RA
         uc_GAA= us_GAA + ud_GAA
         uc_TA = us_TA  + ud_TA
         uc_LA = us_LA  + ud_LA
         uc_UA = us_UA  + ud_UA
         uc_RB = us_RB  + ud_RB
         uc_GBB= us_GBB + ud_GBB
         uc_TB = us_TB  + ud_TB
         uc_LB = us_LB  + ud_LB
         uc_UB = us_UB  + ud_UB

         ! now it's x
         ! we only do these positive x values
         ! according to the  paper the negtive points are
         ! physically not meaningful
         x    = ZERO
         x_RA = ZERO 
         x_GAA= ZERO 
         x_TA = ZERO 
         x_LA = ZERO 
         x_UA = ZERO 
         x_RB = ZERO 
         x_GBB= ZERO 
         x_TB = ZERO 
         x_LB = ZERO 
         x_UB = ZERO 
         uc2  = uc*uc
         IF (DABS(uc)>THRESH .AND. us .LT. ZERO 
     $.AND. (us/uc) .GT. ZERO) THEN
            x    = us/uc
            x_RA = us_RA/uc
            x_GAA= us_GAA/uc
            x_TA = us_TA/uc
            x_LA = us_LA/uc
            x_UA = us_UA/uc
            x_RB = us_RB/uc
            x_GBB= us_GBB/uc
            x_TB = us_TB/uc
            x_LB = us_LB/uc
            x_UB = us_UB/uc

            ! now the other part of deriv
            IF (DABS(uc2)>THRESH) THEN
               x_RA = x_RA  - us*uc_RA/uc2
               x_GAA= x_GAA - us*uc_GAA/uc2
               x_TA = x_TA  - us*uc_TA/uc2
               x_LA = x_LA  - us*uc_LA/uc2
               x_UA = x_UA  - us*uc_UA/uc2
               x_RB = x_RB  - us*uc_RB/uc2
               x_GBB= x_GBB - us*uc_GBB/uc2
               x_TB = x_TB  - us*uc_TB/uc2
               x_LB = x_LB  - us*uc_LB/uc2
               x_UB = x_UB  - us*uc_UB/uc2
            END IF
         END IF
c         write(6,*)"x is x^2*uc", x
c         write(6,*)"uc is x^2*uc", uc
c         write(6,*)"x^2*uc", x*x*uc

         ! now it's x^2*uc
         x2      = x*x
         x2uc    = x2*uc
         x2uc_RA = TWO*x*x_RA*uc  + x2*uc_RA  
         x2uc_GAA= TWO*x*x_GAA*uc + x2*uc_GAA 
         x2uc_TA = TWO*x*x_TA*uc  + x2*uc_TA 
         x2uc_LA = TWO*x*x_LA*uc  + x2*uc_LA 
         x2uc_UA = TWO*x*x_UA*uc  + x2*uc_UA 
         x2uc_RB = TWO*x*x_RB*uc  + x2*uc_RB 
         x2uc_GBB= TWO*x*x_GBB*uc + x2*uc_GBB 
         x2uc_TB = TWO*x*x_TB*uc  + x2*uc_TB 
         x2uc_LB = TWO*x*x_LB*uc  + x2*uc_LB 
         x2uc_UB = TWO*x*x_UB*uc  + x2*uc_UB 

         ! finally let's bring all of term together
         F(i) = F(i) + x2uc

         ! collect all of terms for derivatives
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS)  + x2uc_RA
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + x2uc_GAA
         D1F(i, ID_GAB_POS) = ZERO
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS)  + x2uc_TA
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS)  + x2uc_LA
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + x2uc_UA
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS)  + x2uc_RB
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + x2uc_GBB
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS)  + x2uc_TB
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS)  + x2uc_LB
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + x2uc_UB

      END DO
      RETURN
      END

      SUBROUTINE B13STRONG_AC3(INFOR,NDEN,NG,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the "strong" component for B13 Strong functional.   *
c    *  see the equations in the following paper:                     *
c    *  "Communication: Calibration of a strong-correlation density   *
c    *   functional on transition-metal atoms"                        *
c    *   A. D. Becke  J. Chem. Phys.  138  161101                     *  
c    *                                                                *
c    *  this is the term with x^3*u.                                  *
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
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG
      INTEGER NA,NB
      INTEGER method
      INTEGER I,J
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 HIRWTS(NG)
      REAL*8 HIRWT
      REAL*8 VAL_P,VAL_Q
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
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

      ! static correlation for opposite spin
      ! it's derivatives
      REAL*8 usopp,usopp_RA,usopp_RB,usopp_GAA,usopp_GBB   
      REAL*8 usopp_TA,usopp_TB,usopp_LA,usopp_LB,usopp_UA,usopp_UB

      ! static correlation for parallel spin
      ! it's derivatives
      REAL*8 uspar,uspar_RA,uspar_RB,uspar_GAA,uspar_GBB   
      REAL*8 uspar_TA,uspar_TB,uspar_LA,uspar_LB,uspar_UA,uspar_UB

      ! dynamic correlation for opposite spin
      ! it's derivatives
      REAL*8 udopp,udopp_RA,udopp_RB,udopp_GAA,udopp_GBB   
      REAL*8 udopp_TA,udopp_TB,udopp_LA,udopp_LB,udopp_UA,udopp_UB

      ! dynamic correlation for parallel spin
      ! it's derivatives
      REAL*8 udpar,udpar_RA,udpar_RB,udpar_GAA,udpar_GBB   
      REAL*8 udpar_TA,udpar_TB,udpar_LA,udpar_LB,udpar_UA,udpar_UB

      ! static correlation for both of spin
      ! it's derivatives
      REAL*8 us,us_RA,us_RB,us_GAA,us_GBB   
      REAL*8 us_TA,us_TB,us_LA,us_LB,us_UA,us_UB
      REAL*8 us2,us3

      ! dynamic correlation for both of spin
      ! it's derivatives
      REAL*8 ud,ud_RA,ud_RB,ud_GAA,ud_GBB   
      REAL*8 ud_TA,ud_TB,ud_LA,ud_LB,ud_UA,ud_UB

      ! static correlation for uc 
      REAL*8 uc,uc_RA,uc_RB,uc_GAA,uc_GBB   
      REAL*8 uc_TA,uc_TB,uc_LA,uc_LB,uc_UA,uc_UB
      REAL*8 uc2

      ! the x defined in equation 52 
      REAL*8 x,x_RA,x_RB,x_GAA,x_GBB   
      REAL*8 x_TA,x_TB,x_LA,x_LB,x_UA,x_UB
      REAL*8 x2, x3

      ! the x defined in equation 52, see the above function
      ! this is x^3*uc
      REAL*8 x3uc,x3uc_RA,x3uc_RB,x3uc_GAA,x3uc_GBB   
      REAL*8 x3uc_TA,x3uc_TB,x3uc_LA,x3uc_LB,x3uc_UA,x3uc_UB

      ! tmp vactor for getting the derivatives
      ! totally it's 11 variables including GAB, see above
      ! so we just use the N_FUNC_DERIV_1
      REAL*8 TMP_DRA(3)
      REAL*8 TMP_DRB(3)
      REAL*8 TMP_D1F(N_FUNC_DERIV_1)

      ! constant
      REAL*8 ZERO
      REAL*8 THREE

C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! constants
      ZERO   = 0.0D0
      THREE  = 3.0D0

      ! firstly initilize variable position information
      ! and we initialize all of position infor for deriv
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_TA_POS  = D1VARS(ID_TA)
      ID_LA_POS  = D1VARS(ID_LA)
      ID_EXA_POS = D1VARS(ID_EXA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GBB_POS = D1VARS(ID_GBB)
      ID_TB_POS  = D1VARS(ID_TB)
      ID_LB_POS  = D1VARS(ID_LB)
      ID_EXB_POS = D1VARS(ID_EXB)
      ID_GAB_POS = D1VARS(ID_GAB)
      
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

         ! get the b05 non-dynamic energy density for parallel
         uspar = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL Becke05_NDPAR(INFOR,NDEN,1,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,uspar,TMP_D1F)
         uspar_RA   = TMP_D1F(ID_RA_POS)
         uspar_GAA  = TMP_D1F(ID_GAA_POS)
         uspar_TA   = TMP_D1F(ID_TA_POS)
         uspar_LA   = TMP_D1F(ID_LA_POS)
         uspar_UA   = TMP_D1F(ID_EXA_POS)
         uspar_RB   = TMP_D1F(ID_RB_POS)
         uspar_GBB  = TMP_D1F(ID_GBB_POS)
         uspar_TB   = TMP_D1F(ID_TB_POS)
         uspar_LB   = TMP_D1F(ID_LB_POS)
         uspar_UB   = TMP_D1F(ID_EXB_POS)

         ! get the b05 non-dynamic energy density for opposite
         usopp = ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,
     $VAL_P,HIRWT,
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
c         write(6,*)"our us: ",us
c         write(6,*)"our nus: ",nus
c         write(6,*)"our us, opp: ",usopp
c         write(6,*)"our us, par: ",uspar
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

         ! let's calculate uc
         uc    = us + ud
         uc_RA = us_RA  + ud_RA
         uc_GAA= us_GAA + ud_GAA
         uc_TA = us_TA  + ud_TA
         uc_LA = us_LA  + ud_LA
         uc_UA = us_UA  + ud_UA
         uc_RB = us_RB  + ud_RB
         uc_GBB= us_GBB + ud_GBB
         uc_TB = us_TB  + ud_TB
         uc_LB = us_LB  + ud_LB
         uc_UB = us_UB  + ud_UB

         ! now it's x
         ! we only do these positive x values
         ! according to the  paper the negtive points are
         ! physically not meaningful
         x    = ZERO
         x_RA = ZERO 
         x_GAA= ZERO 
         x_TA = ZERO 
         x_LA = ZERO 
         x_UA = ZERO 
         x_RB = ZERO 
         x_GBB= ZERO 
         x_TB = ZERO 
         x_LB = ZERO 
         x_UB = ZERO 
         uc2  = uc*uc
         IF (DABS(uc)>THRESH .AND. us .LT. ZERO 
     $.AND. (us/uc) .GT. ZERO) THEN
            x    = us/uc
            x_RA = us_RA/uc
            x_GAA= us_GAA/uc
            x_TA = us_TA/uc
            x_LA = us_LA/uc
            x_UA = us_UA/uc
            x_RB = us_RB/uc
            x_GBB= us_GBB/uc
            x_TB = us_TB/uc
            x_LB = us_LB/uc
            x_UB = us_UB/uc

            ! now the other part of deriv
            IF (DABS(uc2)>THRESH) THEN
               x_RA = x_RA  - us*uc_RA/uc2
               x_GAA= x_GAA - us*uc_GAA/uc2
               x_TA = x_TA  - us*uc_TA/uc2
               x_LA = x_LA  - us*uc_LA/uc2
               x_UA = x_UA  - us*uc_UA/uc2
               x_RB = x_RB  - us*uc_RB/uc2
               x_GBB= x_GBB - us*uc_GBB/uc2
               x_TB = x_TB  - us*uc_TB/uc2
               x_LB = x_LB  - us*uc_LB/uc2
               x_UB = x_UB  - us*uc_UB/uc2
            END IF
         END IF
c         write(6,*)"x is x^3*uc", x
c         write(6,*)"uc is x^3*uc", uc
c         write(6,*)"x^3*uc", x*x*x*uc

         ! now it's x^3*uc
         x2      = x*x
         x3      = x*x*x
         x3uc    = x3*uc
         x3uc_RA = THREE*x2*x_RA*uc  + x3*uc_RA  
         x3uc_GAA= THREE*x2*x_GAA*uc + x3*uc_GAA 
         x3uc_TA = THREE*x2*x_TA*uc  + x3*uc_TA 
         x3uc_LA = THREE*x2*x_LA*uc  + x3*uc_LA 
         x3uc_UA = THREE*x2*x_UA*uc  + x3*uc_UA 
         x3uc_RB = THREE*x2*x_RB*uc  + x3*uc_RB 
         x3uc_GBB= THREE*x2*x_GBB*uc + x3*uc_GBB 
         x3uc_TB = THREE*x2*x_TB*uc  + x3*uc_TB 
         x3uc_LB = THREE*x2*x_LB*uc  + x3*uc_LB 
         x3uc_UB = THREE*x2*x_UB*uc  + x3*uc_UB 

         ! finally let's bring all of term together
         F(i) = F(i) + x3uc

         ! collect all of terms for derivatives
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS)  + x3uc_RA
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + x3uc_GAA
         D1F(i, ID_GAB_POS) = ZERO
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS)  + x3uc_TA
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS)  + x3uc_LA
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + x3uc_UA
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS)  + x3uc_RB
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + x3uc_GBB
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS)  + x3uc_TB
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS)  + x3uc_LB
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + x3uc_UB

      END DO
      RETURN
      END
