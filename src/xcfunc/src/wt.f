!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! This is the WT term used in the TPSS functional, which is:
! WT = tau^w/tau;
! tau^w        = tau^w_alpha + tau^w_beta
! tau          = tau_alpha   + tau_beta
! tau^w_sigma  = (1/8)*GXX/RX
! tau_sigma    = tauX
!      
! see the original TPSS paper for more information
! Climbing the density functional ladder: Nonempirical 
! meta-generalized gradient approximation designed for 
! molecules and solids
! J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
! Phys. Rev. Lett., 91 (2003) 146401
!
!

      subroutine functional_wt_close
     &(TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives pointwisely.
! We note that this routine is used in the close shell case,
! that means, we have RA=RB, GAA=GAB=GBB,TA=TB and LA=LB.
! INTPUT :
! rhoA   : the alpha electron density
! rhoB   : the beta  electron density
! DRhoAX : the gradient of rho alpha on x direction
! DRhoAY : the gradient of rho alpha on y direction
! DRhoAZ : the gradient of rho alpha on z direction
! DRhoBX : the gradient of rho beta  on x direction
! DRhoBY : the gradient of rho beta  on y direction
! DRhoBZ : the gradient of rho beta  on z direction
! TauA   : the alpha kinetic energy density
! TauB   : the beta  kinetic energy density
! LapA   : the alpha laplacian density
! LapB   : the beta  laplacian density
! OUTPUT :
! F      : functional values
! D1F    : the first  order functional derivatives
! D2F    : the second order functional derivatives
! D3F    : the third  order functional derivatives
!--------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 
      INTEGER  I
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS
      DOUBLE PRECISION rhoA
      DOUBLE PRECISION RA, RB 
      DOUBLE PRECISION DRhoAX, DRhoAY, DRhoAZ
      DOUBLE PRECISION GAA, GAB, GBB 
      DOUBLE PRECISION TauA
      DOUBLE PRECISION TA, TB
      DOUBLE PRECISION F,D1F(N_FUNC_DERIV_1)
      DOUBLE PRECISION TOL ! tolerance
      
      ! initilize variable position information
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
      DO I = 1, N_FUNC_DERIV_1
         D1F(I) = 0.0D0
      END DO

      ! the real calculation begins
      RA = rhoA
      RB = RA 
      IF (RA.LE.TOL) RETURN
      GAA = DRhoAX*DRhoAX + DRhoAY*DRhoAY
     &+ DRhoAZ*DRhoAZ
      GBB = GAA
      GAB = GAA
      TA  = TauA
      TB  = TA
      t11 = 0.125D0*(GAA+GBB+2.D0*GAB)/(RA+RB)/(0.5D0*TA+0.5D0*TB)
      F   = t11 

      ! derivative
      t4 = (RA+RB)**2
      t13 = -0.125D0*(GAA+GBB+2.D0*GAB)/t4/(0.5D0*TA+0.5D0*TB)
      D1F( ID_RA_POS) = t13 
      t8 = 0.125D0/(RA+RB)/(0.5D0*TA+0.5D0*TB)
      D1F( ID_GAA_POS) = t8 
      t8 = 0.25D0/(RA+RB)/(0.5D0*TA+0.5D0*TB)
      D1F( ID_GAB_POS) = t8 
      t9 = (0.5D0*TA+0.5D0*TB)**2
      t13 = -0.625D-1*(GAA+GBB+2.D0*GAB)/(RA+RB)/t9
      D1F( ID_TA_POS) = t13 
      D1F( ID_TB_POS)= D1F( ID_TA_POS)
      D1F( ID_GBB_POS)= D1F( ID_GAA_POS)
      D1F( ID_RB_POS)= D1F( ID_RA_POS)
      RETURN
      END
                   
      subroutine functional_wt_open
     &(TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,
     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives pointwisely.
! We note that this routine is used in the open shell case.
! INTPUT :
! rhoA   : the alpha electron density
! rhoB   : the beta  electron density
! DRhoAX : the gradient of rho alpha on x direction
! DRhoAY : the gradient of rho alpha on y direction
! DRhoAZ : the gradient of rho alpha on z direction
! DRhoBX : the gradient of rho beta  on x direction
! DRhoBY : the gradient of rho beta  on y direction
! DRhoBZ : the gradient of rho beta  on z direction
! TauA   : the alpha kinetic energy density
! TauB   : the beta  kinetic energy density
! LapA   : the alpha laplacian density
! LapB   : the beta  laplacian density
! OUTPUT :
! F      : functional values
! D1F    : the first  order functional derivatives
! D2F    : the second order functional derivatives
! D3F    : the third  order functional derivatives
!--------------------------------------------------------
#include "fderiv1.inc"
#include "varlist.inc" 
      INTEGER  I
      INTEGER  VAR_INFOR(MAX_VAR_TYPE)
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS
      DOUBLE PRECISION rhoA
      DOUBLE PRECISION RA, RB 
      DOUBLE PRECISION DRhoAX, DRhoAY, DRhoAZ
      DOUBLE PRECISION GAA, GAB, GBB 
      DOUBLE PRECISION TauA
      DOUBLE PRECISION TA, TB
      DOUBLE PRECISION rhoB
      DOUBLE PRECISION DRhoBX, DRhoBY, DRhoBZ
      DOUBLE PRECISION TauB
      DOUBLE PRECISION F,D1F(N_FUNC_DERIV_1)
      DOUBLE PRECISION TOL ! tolerance
      
      ! initilize variable position information
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
      DO I = 1, N_FUNC_DERIV_1
         D1F(I) = 0.0D0
      END DO

      ! the real calculation begins
      RA = rhoA
      RB = rhoB
      GAA = DRhoAX*DRhoAX + DRhoAY*DRhoAY
     &+ DRhoAZ*DRhoAZ
      GBB = DRhoBX*DRhoBX + DRhoBY*DRhoBY
     &+ DRhoBZ*DRhoBZ
      GAB = DRhoAX*DRhoBX + DRhoAY*DRhoBY
     &+ DRhoAZ*DRhoBZ
      TA = TauA
      TB = TauB
      IF (RB.LE.TOL) THEN
         RB = 0.0D0 
         GBB = 0.0D0
         GAB = 0.0D0
         TB = 0.0D0
      END IF 
      IF (RA.LE.TOL) THEN
         RA = 0.0D0 
         GAA = 0.0D0
         GAB = 0.0D0
         TA = 0.0D0
      END IF 
      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN
         RETURN 
      END IF 
             
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN

         ! functional value
         t11 = 0.125D0*(GAA+GBB+2.D0*GAB)/(RA+RB)/(0.5D0*TA+0.5D0*TB)
         F = t11 

         ! derivative
         t4 = (RA+RB)**2
         t13 = -0.125D0*(GAA+GBB+2.D0*GAB)/t4/(0.5D0*TA+0.5D0*TB)
         D1F( ID_RA_POS) = t13 
         t4 = (RA+RB)**2
         t13 = -0.125D0*(GAA+GBB+2.D0*GAB)/t4/(0.5D0*TA+0.5D0*TB)
         D1F( ID_RB_POS) = t13 
         t8 = 0.125D0/(RA+RB)/(0.5D0*TA+0.5D0*TB)
         D1F( ID_GAA_POS)= t8 
         t8 = 0.125D0/(RA+RB)/(0.5D0*TA+0.5D0*TB)
         D1F( ID_GBB_POS)= t8 
         t8 = 0.25D0/(RA+RB)/(0.5D0*TA+0.5D0*TB)
         D1F( ID_GAB_POS)= t8 
         t9 = (0.5D0*TA+0.5D0*TB)**2
         t13 = -0.625D-1*(GAA+GBB+2.D0*GAB)/(RA+RB)/t9
         D1F( ID_TA_POS) = t13 
         t9 = (0.5D0*TA+0.5D0*TB)**2
         t13 = -0.625D-1*(GAA+GBB+2.D0*GAB)/(RA+RB)/t9
         D1F( ID_TB_POS) = t13 

      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN

         ! only RA exists
         t5 = 0.25D0*GAA/RA/TA
         F =  t5 
         t1 = RA**2
         t7 = -0.25D0*GAA/t1/TA
         D1F( ID_RA_POS) =  t7 
         t4 = 0.25D0/RA/TA
         D1F( ID_GAA_POS) = t4 
         t3 = TA**2
         t7 = -0.25D0*GAA/RA/t3
         D1F( ID_TA_POS) =  t7 
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN

         ! only RB exists
         t5 = 0.25D0*GBB/RB/TB
         F = t5 
         t1 = RB**2
         t7 = -0.25D0*GBB/t1/TB
         D1F( ID_RB_POS) =  t7 
         t4 = 0.25D0/RB/TB
         D1F( ID_GBB_POS) = t4 
         t3 = TB**2
         t7 = -0.25D0*GBB/RB/t3
         D1F( ID_TB_POS) =  t7 
      END IF
      RETURN
      END
