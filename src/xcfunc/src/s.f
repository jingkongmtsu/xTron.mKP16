!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! This is the first exchange functional
! see Slater's book:
! J. C. Slater, 
! "The Self-Consistent Field for Molecular and Solids, 
!  Quantum Theory of Molecular and Solids", Vol. 4 
! McGraw-Hill, New York, 1974
!
!



      subroutine functional_s
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)
      IMPLICIT NONE
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives values up to the third order.
! It's only the interface function, the real calculation
! work will be done by seperate routine.
! INPUT :
! INFOR : information related to the variables
! NG    : the number of grid points
! NDEN  : the number of densities
! TOL   : tolerance value for error
! rhoA  : the alpha electron density
! rhoB  : the beta  electron density
! DRhoA : the alpha rho'
! DRhoB : the beta  rho'
! TauA  : the alpha kinetic energy density
! TauB  : the beta  kinetic energy density
! LapA  : the alpha laplacian density
! LapB  : the beta  laplacian density
! OUTPUT:
! F     : functional values
! D1F   : the first  order functional derivatives
! D2F   : the second order functional derivatives
! D3F   : the third  order functional derivatives
!--------------------------------------------------------
      INTEGER INFOR(*)
      INTEGER NG,NDEN
      DOUBLE PRECISION TOL
      DOUBLE PRECISION rhoA(NG),rhoB(NG)
      DOUBLE PRECISION F(NG), D1F(NG,*)
      
      IF (NDEN .EQ. 1) THEN 
      CALL functional_s_close
     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)
      ELSE  
      CALL functional_s_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_s_close
     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives values up to the third order.
! We note that this routine is used in the close shell case,
! that means, we have RA=RB, GAA=GAB=GBB,TA=TB and LA=LB
! for each grid point.
!--------------------------------------------------------
#include "fderiv1.inc"
      INTEGER INFOR(*)
      INTEGER IDERIV, NG
      DOUBLE PRECISION rhoA(NG)
      DOUBLE PRECISION RA, RB ! rho vaule at each point
      DOUBLE PRECISION F(NG),D1F(NG,*)
      DOUBLE PRECISION TOL ! tolerance
      DOUBLE PRECISION PI
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      
      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! deal with the constants
      PI    = 4.D0*DATAN(1.D0)
      IDERIV    = 1

      ! the real calculation begins
      DO i = 1,NG
      RA = rhoA(i)
      RB = RA 
      IF (RA.LE.TOL) CYCLE
      t1 = RA**(1.D0/3.D0)
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = RB**(1.D0/3.D0)
      t11 = -0.1362840444624105D1*t1*RA*t4-0.1362840444624105D1*t7*RB*t4
      F(i) = F(i) +      t11 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = -0.181712059283214D1*t1*t3
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t6 
      END IF
      IF (IDERIV .GE. 1) THEN
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS)= D1F(i, ID_RA_POS)
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
      subroutine functional_s_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives values up to the third order.
! We note that this routine is used in the open shell case.
!--------------------------------------------------------
#include "fderiv1.inc"
      INTEGER INFOR(*)
      INTEGER IDERIV, NG
      DOUBLE PRECISION rhoA(NG)
      DOUBLE PRECISION RA, RB ! rho vaule at each point
      DOUBLE PRECISION rhoB(NG)
      DOUBLE PRECISION DRhoB(NG,3),TauB(NG),LapB(NG)
      DOUBLE PRECISION F(NG),D1F(NG,*)
      DOUBLE PRECISION TOL ! tolerance
      DOUBLE PRECISION PI
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      
      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! deal with the constants
      PI    = 4.D0*DATAN(1.D0)
      IDERIV    = 1

      ! the real calculation begins
      DO i = 1,NG
         RA = rhoA(i)
         RB = rhoB(i)
      IF (RB.LE.TOL) THEN
         RB = 0.0D0 
      END IF 
      IF (RA.LE.TOL) THEN
         RA = 0.0D0 
      END IF 
      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN
         CYCLE 
      END IF 
             
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN
      t1 = RA**(1.D0/3.D0)
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = RB**(1.D0/3.D0)
      t11 = -0.1362840444624105D1*t1*RA*t4-0.1362840444624105D1*t7*RB*t4
      F(i) = F(i) +      t11 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = -0.181712059283214D1*t1*t3
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t6 
      t1 = RB**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = -0.181712059283214D1*t1*t3
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t6 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = RA**(1.D0/3.D0)
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = -0.1362840444624105D1*t1*RA*t4
      F(i) = F(i) +      t7 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = -0.181712059283214D1*t1*t3
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t6 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = RB**(1.D0/3.D0)
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = -0.1362840444624105D1*t1*RB*t4
      F(i) = F(i) +      t7 
      IF (IDERIV .GE. 1) THEN
      t1 = RB**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = -0.181712059283214D1*t1*t3
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t6 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
