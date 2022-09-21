      SUBROUTINE find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
c
c    ******************************************************************
c    *   this is the work function to test the function first         *
c    *   order of derivatives by comparing with find difference       *         
c    *   for function find_diff_check_deriv1_var                      *
c    *                                                                *
c    *   input:                                                       *
c    *   DERIV:      1 - work on rho, 2 - work on gamma etc.          *
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   delta:      the percentage of variation made on the variable *
c    *               for example, 0.01-> new var = var +(-) var*0.01  *         
c    *   criteria:   the criteria the see whether difference between  *
c    *               find difference and analytical one is significant*
c    *   stat:       return whether the comparison is successful      *
c    *                                                                *
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER I,NVAR
      INTEGER stat
      INTEGER NA,NB
      INTEGER DERIV
      REAL*8 THRESH
      REAL*8 delta
      REAL*8 CRITERIA
      REAL*8 RA,RB
      REAL*8 GAA,GBB,GAB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB

      ! this is used to store the analytical derivatives value 
      ! per each grid point
      REAL*8 DFDRA,DFDRB
      REAL*8 DFDGAA,DFDGAB,DFDGBB
      REAL*8 DFDTA,DFDTB
      REAL*8 DFDLA,DFDLB
      REAL*8 DFDUA,DFDUB
      REAL*8 origin

      ! this is the functional value for variable+delta
      REAL*8 FP
      REAL*8 DFPDRA,DFPDRB
      REAL*8 DFPDGAA,DFPDGAB,DFPDGBB
      REAL*8 DFPDTA,DFPDTB
      REAL*8 DFPDLA,DFPDLB
      REAL*8 DFPDUA,DFPDUB
      REAL*8 RA_P,RB_P
      REAL*8 GAA_P,GBB_P,GAB_P 
      REAL*8 TA_P,TB_P
      REAL*8 LA_P,LB_P
      REAL*8 UA_P,UB_P
      REAL*8 var_P

      ! this is the functional value for variable-delta
      REAL*8 FM
      REAL*8 DFMDRA,DFMDRB
      REAL*8 DFMDGAA,DFMDGAB,DFMDGBB
      REAL*8 DFMDTA,DFMDTB
      REAL*8 DFMDLA,DFMDLB
      REAL*8 DFMDUA,DFMDUB
      REAL*8 RA_M,RB_M
      REAL*8 GAA_M,GBB_M,GAB_M 
      REAL*8 TA_M,TB_M
      REAL*8 LA_M,LB_M
      REAL*8 UA_M,UB_M
      REAL*8 var_M

      ! this is used to store the find difference value
      REAL*8 diff

      ! constant
      REAL*8 ZERO
      PARAMETER(ZERO=0.0D0)

      ! set the number of var
      NVAR = 2
      IF (DERIV .eq. 2) THEN
         NVAR = 3
      END IF
      stat = 0

      ! set the NA, NB
      ! because right now the codes only care when
      ! NA/NB is 1 or not
      ! therefore we just set them to non-one case
      NA = 6
      NB = 6

      ! now let's test the alpha and beta deriv
      DO I = 1, NVAR

         ! set the original deriv
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               origin = DFDRA
            ELSE IF(DERIV .eq. 2) THEN
               origin = DFDGAA
            ELSE IF(DERIV .eq. 3) THEN
               origin = DFDTA
            ELSE IF(DERIV .eq. 4) THEN
               origin = DFDLA
            ELSE 
               origin = DFDUA
            END IF
         ELSE IF (I .eq. 2) THEN
            IF (DERIV .eq. 1) THEN
               origin = DFDRB
            ELSE IF(DERIV .eq. 2) THEN
               origin = DFDGBB
            ELSE IF(DERIV .eq. 3) THEN
               origin = DFDTB
            ELSE IF(DERIV .eq. 4) THEN
               origin = DFDLB
            ELSE 
               origin = DFDUB
            END IF
         ELSE
            origin = DFDGAB
         END IF

         ! increase var with delta
         RA_P  = RA 
         RB_P  = RB
         GAA_P = GAA 
         GAB_P = GAB 
         GBB_P = GBB 
         TA_P  = TA  
         TB_P  = TB  
         LA_P  = LA  
         LB_P  = LB  
         UA_P  = UA  
         UB_P  = UB  
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               RA_P  = RA + RA*DELTA 
               var_P = RA_P
            ELSE IF(DERIV .eq. 2) THEN
               GAA_P = GAA + GAA*DELTA 
               var_P = GAA_P
            ELSE IF(DERIV .eq. 3) THEN
               TA_P  = TA + TA*DELTA 
               var_P = TA_P
            ELSE IF(DERIV .eq. 4) THEN
               LA_P  = LA + LA*DELTA 
               var_P = LA_P
            ELSE 
               UA_P  = UA + UA*DELTA 
               var_P = UA_P
            END IF
         ELSE IF (I .eq. 2) THEN
            IF (DERIV .eq. 1) THEN
               RB_P  = RB + RB*DELTA 
               var_P = RB_P
            ELSE IF(DERIV .eq. 2) THEN
               GBB_P = GBB + GBB*DELTA 
               var_P = GBB_P
            ELSE IF(DERIV .eq. 3) THEN
               TB_P  = TB + TB*DELTA 
               var_P = TB_P
            ELSE IF(DERIV .eq. 4) THEN
               LB_P  = LB + LB*DELTA 
               var_P = LB_P
            ELSE 
               UB_P  = UB + UB*DELTA 
               var_P = UB_P
            END IF
         ELSE
            GAB_P  = GAB + GAB*DELTA 
            var_P  = GAB_P
         END IF

         ! now call the functional again
         CALL becke05_f(THRESH,RA_P,GAA_P,TA_P,LA_P,UA_P,
     $ RB_P,GBB_P,TB_P,LB_P,UB_P,FP,DFPDRA,DFPDGAA,DFPDTA,DFPDLA,
     $ DFPDUA,DFPDRB,DFPDGBB,DFPDTB,DFPDLB,DFPDUB)

         ! decrease var with delta
         RA_M  = RA 
         RB_M  = RB
         GAA_M = GAA 
         GAB_M = GAB 
         GBB_M = GBB 
         TA_M  = TA  
         TB_M  = TB  
         LA_M  = LA  
         LB_M  = LB  
         UA_M  = UA  
         UB_M  = UB  
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               RA_M  = RA - RA*DELTA 
               var_M = RA_M
            ELSE IF(DERIV .eq. 2) THEN
               GAA_M = GAA - GAA*DELTA 
               var_M = GAA_M
            ELSE IF(DERIV .eq. 3) THEN
               TA_M  = TA - TA*DELTA 
               var_M = TA_M
            ELSE IF(DERIV .eq. 4) THEN
               LA_M  = LA - LA*DELTA 
               var_M = LA_M
            ELSE 
               UA_M  = UA - UA*DELTA 
               var_M = UA_M
            END IF
         ELSE IF (I .eq. 2) THEN
            IF (DERIV .eq. 1) THEN
               RB_M  = RB - RB*DELTA 
               var_M = RB_M
            ELSE IF(DERIV .eq. 2) THEN
               GBB_M = GBB - GBB*DELTA 
               var_M = GBB_M
            ELSE IF(DERIV .eq. 3) THEN
               TB_M  = TB - TB*DELTA 
               var_M = TB_M
            ELSE IF(DERIV .eq. 4) THEN
               LB_M  = LB - LB*DELTA 
               var_M = LB_M
            ELSE 
               UB_M  = UB - UB*DELTA 
               var_M = UB_M
            END IF
         ELSE
            GAB_M  = GAB - GAB*DELTA 
            var_M  = GAB_M
         END IF

         ! now call the functional again
         CALL becke05_f(THRESH,RA_M,GAA_M,TA_M,LA_M,UA_M,
     $ RB_M,GBB_M,TB_M,LB_M,UB_M,FM,DFMDRA,DFMDGAA,DFMDTA,DFMDLA,
     $ DFMDUA,DFMDRB,DFMDGBB,DFMDTB,DFMDLB,DFMDUB)

         ! now let's see the result
         IF (DABS(var_P-var_M)>THRESH) THEN
            diff = (FP-FM)/(var_P-var_M)
            IF(DABS(diff-origin)>CRITERIA) THEN
               stat = 1
               IF (I .eq. 1) THEN
                  IF (DERIV .eq. 1) THEN
                     write(6,*)"find diff fail on RA"
                  ELSE IF(DERIV .eq. 2) THEN
                     write(6,*)"find diff fail on GAA"
                  ELSE IF(DERIV .eq. 3) THEN
                     write(6,*)"find diff fail on TA"
                  ELSE IF(DERIV .eq. 4) THEN
                     write(6,*)"find diff fail on LA"
                  ELSE
                     write(6,*)"find diff fail on UA"
                  END IF
               ELSE IF (I .eq. 2) THEN
                  IF (DERIV .eq. 1) THEN
                     write(6,*)"find diff fail on RB"
                  ELSE IF(DERIV .eq. 2) THEN
                     write(6,*)"find diff fail on GBB"
                  ELSE IF(DERIV .eq. 3) THEN
                     write(6,*)"find diff fail on TB"
                  ELSE IF(DERIV .eq. 4) THEN
                     write(6,*)"find diff fail on LB"
                  ELSE
                     write(6,*)"find diff fail on UB"
                  END IF
               ELSE
                  write(6,*)"find diff fail on GAB"
               END IF
               write(6,*)"on deriv, diff is, and percentage is", 
     $ DABS(diff-origin), DABS(diff-origin)/DABS(origin)
               write(6,*)"original  deriv value", origin
               write(6,*)"find diff deriv value", diff
            END IF
         END IF
      END DO
      RETURN
      END

      SUBROUTINE find_diff_check_deriv1_var(INFOR,NG,THRESH,Delta,
     $CRITERIA,RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB)
c
c    ******************************************************************
c    *   this is the driver function to test the function first       *
c    *   order of derivatives by comparing with find difference       *         
c    *   we note that this function is calling the testing functional *
c    *   in variable form input/output(that is to say, the rho, drho  *
c    *   and the output F,D1F etc. are not in array form)             * 
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   delta:      the percentage of variation made on the variable *
c    *               for example, 0.01-> new var = var +(-) var*0.01  *         
c    *   criteria:   the criteria the see whether difference between  *
c    *               find difference and analytical one is significant*
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NG,I,stat
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 DELTA
      REAL*8 CRITERIA
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 RA,RB
      REAL*8 GAA,GBB,GAB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB

      ! status check
      logical has_rho
      logical has_gamma
      logical has_tau
      logical has_lap
      logical has_exrho

      ! deriv status
      INTEGER deriv

      ! this is used to store the analytical derivatives value 
      ! per each grid point
      REAL*8 F
      REAL*8 DFDRA,DFDRB
      REAL*8 DFDGAA,DFDGAB,DFDGBB
      REAL*8 DFDTA,DFDTB
      REAL*8 DFDLA,DFDLB
      REAL*8 DFDUA,DFDUB

      ! constant
      REAL*8 ZERO
      PARAMETER(ZERO=0.0D0)

      ! set the NA, NB
      ! because right now the codes only care when
      ! NA/NB is 1 or not
      ! therefore we just set them to non-one case
      NA = 6
      NB = 6

      ! let's check what kind of variable we have
      has_rho    = .false.
      has_gamma  = .false.
      has_tau    = .false.
      has_lap    = .false.
      has_exrho  = .false.
      IF (INFOR(1) .gt. 0) THEN
         has_rho   = .true.
      END IF
      IF (INFOR(2) .gt. 0) THEN
         has_gamma = .true.
      END IF
      IF (INFOR(3) .gt. 0) THEN
         has_tau   = .true.
      END IF
      IF (INFOR(4) .gt. 0) THEN
         has_lap   = .true.
      END IF
      IF (INFOR(5) .gt. 0) THEN
         has_exrho = .true.
      END IF
      
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
         GAB = DRA(i,1)*DRB(i,1) + DRA(i,2)*DRB(i,2)
     &+ DRA(i,3)*DRB(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)

         ! debug
c         write(6,*)"our RhoA: ",RA
c         write(6,*)"our GAA: ",GAA 
c         write(6,*)"our GAB: ",GAB
c         write(6,*)"our TA: ",TA
c         write(6,*)"our LA: ",LA 
c         write(6,*)"our UA: ",UA
c         write(6,*)"our RhoB: ",RB
c         write(6,*)"our GBB: ",GBB 
c         write(6,*)"our TB: ",TB
c         write(6,*)"our LB: ",LB 
c         write(6,*)"our UB: ",UB 

         ! let's generate the analytical deriv
         F      = ZERO
         DFDRA  = ZERO
         DFDRB  = ZERO
         DFDGAA = ZERO
         DFDGAB = ZERO
         DFDGBB = ZERO
         DFDTA  = ZERO
         DFDTB  = ZERO
         DFDLA  = ZERO
         DFDLB  = ZERO
         DFDUA  = ZERO
         DFDUB  = ZERO

         ! here we need to add a function call
         ! we leave it to be blank
         CALL becke05_f(THRESH,RA,GAA,TA,LA,UA,
     $ RB,GBB,TB,LB,UB,F,DFDRA,DFDGAA,DFDTA,DFDLA,
     $ DFDUA,DFDRB,DFDGBB,DFDTB,DFDLB,DFDUB)

         ! now do rho
         IF (has_rho) THEN
            deriv = 1
            CALL find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do gamma
         IF (has_gamma) THEN
            deriv = 2
            CALL find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do tau
         IF (has_tau) THEN
            deriv = 3
            CALL find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do lap
         IF (has_lap) THEN
            deriv = 4
            CALL find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do exchange energy density
         IF (has_exrho) THEN
            deriv = 5
            CALL find_diff_check_deriv1_var_core(DERIV,THRESH,Delta,
     $CRITERIA,RA,RB,GAA,GAB,GBB,TA,TB,LA,LB,UA,UB,DFDRA,DFDRB,DFDGAA,
     $DFDGAB,DFDGBB,DFDTA,DFDTB,DFDLA,DFDLB,DFDUA,DFDUB,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

      END DO
      RETURN
      END


      SUBROUTINE find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,DRA,DRB,TA,TB,
     $LA,LB,UA,UB,D1F,HIRWT,stat)
c
c    ******************************************************************
c    *   this is the work function to test the function first         *
c    *   order of derivatives by comparing with find difference       *         
c    *   for function find_diff_check_deriv1_array                    *
c    *                                                                *
c    *   note:                                                        *
c    *   we can not test the derivatives for GAB in this function     *
c    *   form. Because increase/decrease the GAB by varying DRA/DRB,  *
c    *   will accordingly change the variable of GAA or GBB.          *
c    *   this will make the find difference not applicable.           *      
c    *                                                                *
c    *   input:                                                       *
c    *   DERIV:      1 - work on rho, 2 - work on gamma etc.          *
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   delta:      the percentage of variation made on the variable *
c    *               for example, 0.01-> new var = var +(-) var*0.01  *         
c    *   criteria:   the criteria the see whether difference between  *
c    *               find difference and analytical one is significant*
c    *   stat:       return whether the comparison is successful      *
c    *                                                                *
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE

      ! functional deriv position information
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS
      INTEGER INFOR(*)

      ! the variable input etc.
      INTEGER I,NVAR
      INTEGER stat
      INTEGER NA,NB,NDEN
      INTEGER DERIV
      REAL*8 THRESH
      REAL*8 delta
      REAL*8 CRITERIA
      REAL*8 RA,RB
      REAL*8 DRA(3)
      REAL*8 DRB(3)
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 D1F(N_FUNC_DERIV_1)

      ! this is used to store the analytical derivatives value 
      ! per each grid point
      REAL*8 DFDRA,DFDRB
      REAL*8 DFDGAA,DFDGAB,DFDGBB
      REAL*8 DFDTA,DFDTB
      REAL*8 DFDLA,DFDLB
      REAL*8 DFDUA,DFDUB
      REAL*8 origin

      ! this is the functional value for variable+delta
      REAL*8 FP
      REAL*8 D1FP(N_FUNC_DERIV_1)
      REAL*8 RA_P,RB_P
      REAL*8 DRA_P(3)
      REAL*8 DRB_P(3)
      REAL*8 GAA_P, GBB_P
      REAL*8 TA_P,TB_P
      REAL*8 LA_P,LB_P
      REAL*8 UA_P,UB_P
      REAL*8 var_P

      ! this is the functional value for variable-delta
      REAL*8 FM
      REAL*8 D1FM(N_FUNC_DERIV_1)
      REAL*8 RA_M,RB_M
      REAL*8 DRA_M(3)
      REAL*8 DRB_M(3)
      REAL*8 GAA_M, GBB_M
      REAL*8 TA_M,TB_M
      REAL*8 LA_M,LB_M
      REAL*8 UA_M,UB_M
      REAL*8 var_M

      ! this are parameters used by kp14c
      REAL*8 kp14c_b, kp14c_alpha, kp14c_ndpar

      ! this is used to store the find difference value
      REAL*8 diff

      ! the input hirshfeld weights
      REAL*8 HIRWT

      ! constant
      REAL*8 ZERO
      PARAMETER(ZERO=0.0D0)

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      ID_RA_POS=D1VARS(ID_RA)
      ID_GAA_POS=D1VARS(ID_GAA)
      ID_GAB_POS=D1VARS(ID_GAB)
      ID_TA_POS=D1VARS(ID_TA)
      ID_LA_POS=D1VARS(ID_LA)
      ID_EXA_POS=D1VARS(ID_EXA)
      ID_RB_POS=D1VARS(ID_RB)
      ID_GBB_POS=D1VARS(ID_GBB)
      ID_TB_POS=D1VARS(ID_TB)
      ID_LB_POS=D1VARS(ID_LB)
      ID_EXB_POS=D1VARS(ID_EXB)
      stat = 0

      ! now get the analytical deriv value
      ! initialize them first
      DFDRA  = ZERO
      DFDRB  = ZERO
      DFDGAA = ZERO
      DFDGAB = ZERO
      DFDGBB = ZERO
      DFDTA  = ZERO
      DFDTB  = ZERO
      DFDLA  = ZERO
      DFDLB  = ZERO
      DFDUA  = ZERO
      DFDUB  = ZERO

      ! rho
      IF (ID_RA_POS .gt. 0) THEN
         DFDRA  = D1F(ID_RA_POS)
      END IF
      IF (ID_RB_POS .gt. 0) THEN
         DFDRB  = D1F(ID_RB_POS)
      END IF

      ! gamma
      IF (ID_GAA_POS .gt. 0) THEN
         DFDGAA  = D1F(ID_GAA_POS)
      END IF
      IF (ID_GBB_POS .gt. 0) THEN
         DFDGBB  = D1F(ID_GBB_POS)
      END IF

      ! tau
      IF (ID_TA_POS .gt. 0) THEN
         DFDTA  = D1F(ID_TA_POS)
      END IF
      IF (ID_TB_POS .gt. 0) THEN
         DFDTB  = D1F(ID_TB_POS)
      END IF

      ! lap
      IF (ID_LA_POS .gt. 0) THEN
         DFDLA  = D1F(ID_LA_POS)
      END IF
      IF (ID_LB_POS .gt. 0) THEN
         DFDLB  = D1F(ID_LB_POS)
      END IF

      ! exchange energy density
      IF (ID_EXA_POS .gt. 0) THEN
         DFDUA  = D1F(ID_EXA_POS)
      END IF
      IF (ID_EXB_POS .gt. 0) THEN
         DFDUB  = D1F(ID_EXB_POS)
      END IF

      ! set the number of var
      NVAR = 2

      ! now let's test the alpha and beta deriv
      DO I = 1, NVAR

         ! set the original deriv
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               origin = DFDRA
            ELSE IF(DERIV .eq. 2) THEN
               origin = DFDGAA
            ELSE IF(DERIV .eq. 3) THEN
               origin = DFDTA
            ELSE IF(DERIV .eq. 4) THEN
               origin = DFDLA
            ELSE 
               origin = DFDUA
            END IF
         ELSE 
            IF (DERIV .eq. 1) THEN
               origin = DFDRB
            ELSE IF(DERIV .eq. 2) THEN
               origin = DFDGBB
            ELSE IF(DERIV .eq. 3) THEN
               origin = DFDTB
            ELSE IF(DERIV .eq. 4) THEN
               origin = DFDLB
            ELSE 
               origin = DFDUB
            END IF
         END IF

         ! increase var with delta
         RA_P  = RA 
         RB_P  = RB
         TA_P  = TA  
         TB_P  = TB  
         LA_P  = LA  
         LB_P  = LB  
         UA_P  = UA  
         UB_P  = UB  
         DRA_P(1) = DRA(1)
         DRA_P(2) = DRA(2)
         DRA_P(3) = DRA(3)
         DRB_P(1) = DRB(1)
         DRB_P(2) = DRB(2)
         DRB_P(3) = DRB(3)
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               RA_P  = RA + RA*DELTA 
               var_P = RA_P
            ELSE IF(DERIV .eq. 2) THEN
               DRA_P(1) = DRA(1) + DRA(1)*delta
               DRA_P(2) = DRA(2) + DRA(2)*delta
               DRA_P(3) = DRA(3) + DRA(3)*delta
               GAA_P    = DRA_P(1)*DRA_P(1)+DRA_P(2)*DRA_P(2)+
     $ DRA_P(3)*DRA_P(3)
               var_P = GAA_P
            ELSE IF(DERIV .eq. 3) THEN
               TA_P  = TA + TA*DELTA 
               var_P = TA_P
            ELSE IF(DERIV .eq. 4) THEN
               LA_P  = LA + LA*DELTA 
               var_P = LA_P
            ELSE 
               UA_P  = UA + UA*DELTA 
               var_P = UA_P
            END IF
         ELSE 
            IF (DERIV .eq. 1) THEN
               RB_P  = RB + RB*DELTA 
               var_P = RB_P
            ELSE IF(DERIV .eq. 2) THEN
               DRB_P(1) = DRB(1) + DRB(1)*delta
               DRB_P(2) = DRB(2) + DRB(2)*delta
               DRB_P(3) = DRB(3) + DRB(3)*delta
               GBB_P    = DRB_P(1)*DRB_P(1)+DRB_P(2)*DRB_P(2)+
     $ DRB_P(3)*DRB_P(3)
               var_P = GBB_P
            ELSE IF(DERIV .eq. 3) THEN
               TB_P  = TB + TB*DELTA 
               var_P = TB_P
            ELSE IF(DERIV .eq. 4) THEN
               LB_P  = LB + LB*DELTA 
               var_P = LB_P
            ELSE 
               UB_P  = UB + UB*DELTA 
               var_P = UB_P
            END IF
         END IF

         ! now call the functional again
         FP = ZERO
C         CALL Becke05_NDPAR(INFOR,NDEN,1,THRESH,
C     $RA_P,RB_P,DRA_P,DRB_P,TA_P,TB_P,LA_P,LB_P,UA_P,UB_P,FP,D1FP)
C         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,
C     $RA_P,RB_P,DRA_P,DRB_P,TA_P,TB_P,LA_P,LB_P,UA_P,UB_P,FP,D1FP)
C         CALL KP14EC(INFOR,NDEN,1,NA,NB,
C     $kp14c_b,kp14c_alpha,kp14c_ndpar,THRESH,
C     $RA_P,RB_P,DRA_P,DRB_P,TA_P,TB_P,LA_P,LB_P,UA_P,UB_P,FP,D1FP)
         CALL B13STRONG_AC2(INFOR,NDEN,1,NA,NB,THRESH,HIRWT,
     $RA_P,RB_P,DRA_P,DRB_P,TA_P,TB_P,LA_P,LB_P,UA_P,UB_P,FP,D1FP)

         ! decrease var with delta
         RA_M  = RA 
         RB_M  = RB
         TA_M  = TA  
         TB_M  = TB  
         LA_M  = LA  
         LB_M  = LB  
         UA_M  = UA  
         UB_M  = UB  
         DRA_M(1) = DRA(1)
         DRA_M(2) = DRA(2)
         DRA_M(3) = DRA(3)
         DRB_M(1) = DRB(1)
         DRB_M(2) = DRB(2)
         DRB_M(3) = DRB(3)
         IF (I .eq. 1) THEN
            IF (DERIV .eq. 1) THEN
               RA_M  = RA - RA*DELTA 
               var_M = RA_M
            ELSE IF(DERIV .eq. 2) THEN
               DRA_M(1) = DRA(1) - DRA(1)*delta
               DRA_M(2) = DRA(2) - DRA(2)*delta
               DRA_M(3) = DRA(3) - DRA(3)*delta
               GAA_M    = DRA_M(1)*DRA_M(1)+DRA_M(2)*DRA_M(2)+
     $ DRA_M(3)*DRA_M(3)
               var_M = GAA_M
            ELSE IF(DERIV .eq. 3) THEN
               TA_M  = TA - TA*DELTA 
               var_M = TA_M
            ELSE IF(DERIV .eq. 4) THEN
               LA_M  = LA - LA*DELTA 
               var_M = LA_M
            ELSE 
               UA_M  = UA - UA*DELTA 
               var_M = UA_M
            END IF
         ELSE 
            IF (DERIV .eq. 1) THEN
               RB_M  = RB - RB*DELTA 
               var_M = RB_M
            ELSE IF(DERIV .eq. 2) THEN
               DRB_M(1) = DRB(1) - DRB(1)*delta
               DRB_M(2) = DRB(2) - DRB(2)*delta
               DRB_M(3) = DRB(3) - DRB(3)*delta
               GBB_M    = DRB_M(1)*DRB_M(1)+DRB_M(2)*DRB_M(2)+
     $ DRB_M(3)*DRB_M(3)
               var_M = GBB_M
            ELSE IF(DERIV .eq. 3) THEN
               TB_M  = TB - TB*DELTA 
               var_M = TB_M
            ELSE IF(DERIV .eq. 4) THEN
               LB_M  = LB - LB*DELTA 
               var_M = LB_M
            ELSE 
               UB_M  = UB - UB*DELTA 
               var_M = UB_M
            END IF
         END IF

         ! now call the functional again
         FM = ZERO
C         CALL Becke05_NDPAR(INFOR,NDEN,1,THRESH,
C     $RA_M,RB_M,DRA_M,DRB_M,TA_M,TB_M,LA_M,LB_M,UA_M,UB_M,FM,D1FM)
C         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,
C     $RA_M,RB_M,DRA_M,DRB_M,TA_M,TB_M,LA_M,LB_M,UA_M,UB_M,FM,D1FM)
C         CALL KP14EC(INFOR,NDEN,1,NA,NB,
C     $kp14c_b,kp14c_alpha,kp14c_ndpar,THRESH,
C     $RA_M,RB_M,DRA_M,DRB_M,TA_M,TB_M,LA_M,LB_M,UA_M,UB_M,FM,D1FM)
         CALL B13STRONG_AC2(INFOR,NDEN,1,NA,NB,THRESH,HIRWT,
     $RA_M,RB_M,DRA_M,DRB_M,TA_M,TB_M,LA_M,LB_M,UA_M,UB_M,FM,D1FM)

         ! now let's see the result
         IF (DABS(var_P-var_M)>THRESH) THEN
            diff = (FP-FM)/(var_P-var_M)
            IF(DABS(diff-origin)>CRITERIA) THEN
               stat = 1
               IF (I .eq. 1) THEN
                  IF (DERIV .eq. 1) THEN
                     write(6,*)"find diff fail on RA"
                  ELSE IF(DERIV .eq. 2) THEN
                     write(6,*)"find diff fail on GAA"
                  ELSE IF(DERIV .eq. 3) THEN
                     write(6,*)"find diff fail on TA"
                  ELSE IF(DERIV .eq. 4) THEN
                     write(6,*)"find diff fail on LA"
                  ELSE
                     write(6,*)"find diff fail on UA"
                  END IF
               ELSE 
                  IF (DERIV .eq. 1) THEN
                     write(6,*)"find diff fail on RB"
                  ELSE IF(DERIV .eq. 2) THEN
                     write(6,*)"find diff fail on GBB"
                  ELSE IF(DERIV .eq. 3) THEN
                     write(6,*)"find diff fail on TB"
                  ELSE IF(DERIV .eq. 4) THEN
                     write(6,*)"find diff fail on LB"
                  ELSE
                     write(6,*)"find diff fail on UB"
                  END IF
               END IF
               write(6,*)"on deriv, diff is, and percentage is", 
     $ DABS(diff-origin), DABS(diff-origin)/DABS(origin)
               write(6,*)"original  deriv value", origin
               write(6,*)"find diff deriv value", diff
            END IF
         END IF
      END DO
      RETURN
      END

      SUBROUTINE find_diff_check_deriv1_array(INFOR,NG,THRESH,Delta,
     $CRITERIA,RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB)
c
c    ******************************************************************
c    *   this is the driver function to test the function first       *
c    *   order of derivatives by comparing with find difference       *         
c    *   we note that this function is calling the testing functional *
c    *   in array form input/output                                   *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   delta:      the percentage of variation made on the variable *
c    *               for example, 0.01-> new var = var +(-) var*0.01  *         
c    *   criteria:   the criteria the see whether difference between  *
c    *               find difference and analytical one is significant*
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   fenglai liu                                                  *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
#include "fderiv1.inc"
      INTEGER NG,I,J
      INTEGER NA,NB,NDEN
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 DELTA
      REAL*8 CRITERIA
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 RA,RB
      REAL*8 TMP_DRA(3)
      REAL*8 TMP_DRB(3)
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 GAA,GAB,GBB
      REAL*8 HIRWT
      REAL*8 F
      REAL*8 D1F(N_FUNC_DERIV_1)

      ! status check
      logical has_rho
      logical has_gamma
      logical has_tau
      logical has_lap
      logical has_exrho

      ! deriv status
      INTEGER deriv
      INTEGER stat

      ! this are parameters used by kp14c
      REAL*8 kp14c_b, kp14c_alpha, kp14c_ndpar

      ! constant
      REAL*8 ZERO
      PARAMETER(ZERO=0.0D0)

      ! set the NA, NB
      ! because right now the codes only care when
      ! NA/NB is 1 or not
      ! therefore we just set them to non-one case
      NA = 6
      NB = 6

      ! set the number of density
      NDEN = 2

      ! set the KP14C parameters
      ! if we call KP14C, you may want to change it here
      ! make sure you use the same value with the above
      ! core function
      kp14c_b     = 1.355D0
      kp14c_alpha = 0.038D0
      kp14c_ndpar = 1.128D0

      ! let's check what kind of variable we have
      has_rho    = .false.
      has_gamma  = .false.
      has_tau    = .false.
      has_lap    = .false.
      has_exrho  = .false.
      IF (INFOR(1) .gt. 0) THEN
         has_rho   = .true.
      END IF
      IF (INFOR(2) .gt. 0) THEN
         has_gamma = .true.
      END IF
      IF (INFOR(3) .gt. 0) THEN
         has_tau   = .true.
      END IF
      IF (INFOR(4) .gt. 0) THEN
         has_lap   = .true.
      END IF
      IF (INFOR(5) .gt. 0) THEN
         has_exrho = .true.
      END IF
      
      ! loop over NG
      DO I = 1,NG

c         write(6,*)"     "
c         write(6,*)"----------------------"
c         write(6,*)"for grid: ", i
          
         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)
         TMP_DRA(1) = DRA(i,1)
         TMP_DRA(2) = DRA(i,2)
         TMP_DRA(3) = DRA(i,3)
         TMP_DRB(1) = DRB(i,1)
         TMP_DRB(2) = DRB(i,2)
         TMP_DRB(3) = DRB(i,3)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)

         ! fake up the hirshfeld weight if needed
         HIRWT = 0.1E0

         ! debug
c         write(6,*)"our RhoA: ",RA
c         write(6,*)"our GAA: ", GAA
c         write(6,*)"our TA: ",TA
c         write(6,*)"our LA: ",LA 
c         write(6,*)"our UA: ",UA
c         write(6,*)"our RhoB: ",RB
c         write(6,*)"our GBB: ", GBB
c         write(6,*)"our TB: ",TB
c         write(6,*)"our LB: ",LB 
c         write(6,*)"our UB: ",UB 

         ! here we calculate the analytical derivatives
         F = ZERO
         DO J = 1, N_FUNC_DERIV_1
            D1F(J) = ZERO
         END DO
c         CALL Becke05_NDPAR(INFOR,NDEN,1,THRESH,
c     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,F,D1F)
c         CALL Becke05_NDOP(INFOR,NDEN,1,THRESH,
c     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,F,D1F)
c         CALL KP14EC(INFOR,NDEN,1,NA,NB,
c     $kp14c_b,kp14c_alpha,kp14c_ndpar,THRESH,
c     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,F,D1F)
         CALL B13STRONG_AC2(INFOR,NDEN,1,NA,NB,THRESH,HIRWT,
     $RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,F,D1F)

         ! now do rho
         IF (has_rho) THEN
            deriv = 1
            CALL find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $D1F,HIRWT,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do gamma
         IF (has_gamma) THEN
            deriv = 2
            CALL find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $D1F,HIRWT,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do tau
         IF (has_tau) THEN
            deriv = 3
            CALL find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $D1F,HIRWT,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do lap
         IF (has_lap) THEN
            deriv = 4
            CALL find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $D1F,HIRWT,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

         ! now do exchange energy density
         IF (has_exrho) THEN
            deriv = 5
            CALL find_diff_check_deriv1_array_core(INFOR,DERIV,
     $NDEN,NA,NB,kp14c_b,kp14c_alpha,kp14c_ndpar,
     $THRESH,Delta,CRITERIA,RA,RB,TMP_DRA,TMP_DRB,TA,TB,LA,LB,UA,UB,
     $D1F,HIRWT,stat)
            IF (stat .gt. 0) THEN
               write(6,*)"find difference failed on grid point", I
            END IF
         END IF

      END DO
      RETURN
      END

