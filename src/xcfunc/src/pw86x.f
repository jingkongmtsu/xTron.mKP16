!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! J. P. Perdew, Y. Wang, Phys. Rev. B, 33, 8800 (1986) 
!
!



      subroutine functional_pw86x
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
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
      DOUBLE PRECISION DRhoA(NG,3),DRhoB(NG,3)
      DOUBLE PRECISION F(NG), D1F(NG,*)
      
      IF (NDEN .EQ. 1) THEN 
      CALL functional_pw86x_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
      ELSE  
      CALL functional_pw86x_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_pw86x_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
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
      DOUBLE PRECISION DRhoA(NG,3)
      DOUBLE PRECISION GAA, GAB, GBB ! the gamma value at each point
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
      GAA = DRhoA(i,1)*DRhoA(i,1) + DRhoA(i,2)*DRhoA(i,2)
     &+ DRhoA(i,3)*DRhoA(i,3)
      GBB = GAA
      GAB = GAA
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RA**2
      t12 = t3**2
      t18 = 1/t7/t6
      t19 = GAA**2
      t21 = t11**2
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t36 = (1.D0+0.9812451208562037D-1*t9*GAA/t12/t11+0.802551762691248
     &4D-1*t18*t19/t3/t21/RA+0.8680555571D-4*t28*t19*GAA/t31)**(1.D0/15.
     &D0)
      t39 = RB**(1.D0/3.D0)
      t43 = RB**2
      t44 = t39**2
      t49 = GBB**2
      t51 = t43**2
      t59 = t51**2
      t64 = (1.D0+0.9812451208562037D-1*t9*GBB/t44/t43+0.802551762691248
     &4D-1*t18*t49/t39/t51/RB+0.8680555571D-4*t28*t49*GBB/t59)**(1.D0/15
     &.D0)
      t67 = -0.1362840444624105D1*t2*t3*RA*t36-0.1362840444624105D1*t2*t
     &39*RB*t64
      F(i) = F(i) +      t67 
      IF (IDERIV .GE. 1) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t5 = PI**2
      t6 = t5**(1.D0/3.D0)
      t7 = t6**2
      t9 = 1/t7*GAA
      t10 = RA**2
      t11 = t3**2
      t18 = GAA**2
      t19 = 1/t6/t5*t18
      t20 = t10**2
      t26 = t5**2
      t29 = 1/t26*t18*GAA
      t30 = t20**2
      t35 = (1.D0+0.9812451208562037D-1*t9/t11/t10+0.8025517626912484D-1
     &*t19/t3/t20/RA+0.8680555571D-4*t29/t30)**(1.D0/15.D0)
      t40 = t35**2
      t41 = t40**2
      t43 = t41**2
      t64 = -0.181712059283214D1*t2*t3*t35-0.9085602964160698D-1*t2*t3*R
     &A/t43/t41/t40*(-0.2616653655616543D0*t9/t11/t10/RA-0.4280276067686
     &658D0*t19/t3/t20/t10-0.69444444568D-3*t29/t30/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t64 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RA**2
      t12 = t3**2
      t14 = 1/t12/t11
      t18 = 1/t7/t6
      t19 = GAA**2
      t21 = t11**2
      t24 = 1/t3/t21/RA
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t32 = 1/t31
      t36 = (1.D0+0.9812451208562037D-1*t9*GAA*t14+0.8025517626912484D-1
     &*t18*t19*t24+0.8680555571D-4*t28*t19*GAA*t32)**(1.D0/15.D0)
      t37 = t36**2
      t38 = t37**2
      t40 = t38**2
      t55 = -0.9085602964160698D-1*t2*t3*RA/t40/t38/t37*(0.9812451208562
     &037D-1*t9*t14+0.1605103525382497D0*t18*GAA*t24+0.26041666713D-3*t2
     &8*t19*t32)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t55 
      t1 = 0.D0
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t1 
      END IF
      IF (IDERIV .GE. 1) THEN
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS)= D1F(i, ID_GAA_POS)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS)= D1F(i, ID_RA_POS)
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
      subroutine functional_pw86x_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
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
      DOUBLE PRECISION DRhoA(NG,3)
      DOUBLE PRECISION GAA, GAB, GBB ! the gamma value at each point
      DOUBLE PRECISION rhoB(NG)
      DOUBLE PRECISION DRhoB(NG,3)
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
         GAA = DRhoA(i,1)*DRhoA(i,1) + DRhoA(i,2)*DRhoA(i,2)
     &+ DRhoA(i,3)*DRhoA(i,3)
         GBB = DRhoB(i,1)*DRhoB(i,1) + DRhoB(i,2)*DRhoB(i,2)
     &+ DRhoB(i,3)*DRhoB(i,3)
         GAB = DRhoA(i,1)*DRhoB(i,1) + DRhoA(i,2)*DRhoB(i,2)
     &+ DRhoA(i,3)*DRhoB(i,3)
      IF (RB.LE.TOL) THEN
         RB = 0.0D0 
         GBB = 0.0D0
         GAB = 0.0D0
      END IF 
      IF (RA.LE.TOL) THEN
         RA = 0.0D0 
         GAA = 0.0D0
         GAB = 0.0D0
      END IF 
      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN
         CYCLE 
      END IF 
             
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RA**2
      t12 = t3**2
      t18 = 1/t7/t6
      t19 = GAA**2
      t21 = t11**2
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t36 = (1.D0+0.9812451208562037D-1*t9*GAA/t12/t11+0.802551762691248
     &4D-1*t18*t19/t3/t21/RA+0.8680555571D-4*t28*t19*GAA/t31)**(1.D0/15.
     &D0)
      t39 = RB**(1.D0/3.D0)
      t43 = RB**2
      t44 = t39**2
      t49 = GBB**2
      t51 = t43**2
      t59 = t51**2
      t64 = (1.D0+0.9812451208562037D-1*t9*GBB/t44/t43+0.802551762691248
     &4D-1*t18*t49/t39/t51/RB+0.8680555571D-4*t28*t49*GBB/t59)**(1.D0/15
     &.D0)
      t67 = -0.1362840444624105D1*t2*t3*RA*t36-0.1362840444624105D1*t2*t
     &39*RB*t64
      F(i) = F(i) +      t67 
      IF (IDERIV .GE. 1) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t5 = PI**2
      t6 = t5**(1.D0/3.D0)
      t7 = t6**2
      t9 = 1/t7*GAA
      t10 = RA**2
      t11 = t3**2
      t18 = GAA**2
      t19 = 1/t6/t5*t18
      t20 = t10**2
      t26 = t5**2
      t29 = 1/t26*t18*GAA
      t30 = t20**2
      t35 = (1.D0+0.9812451208562037D-1*t9/t11/t10+0.8025517626912484D-1
     &*t19/t3/t20/RA+0.8680555571D-4*t29/t30)**(1.D0/15.D0)
      t40 = t35**2
      t41 = t40**2
      t43 = t41**2
      t64 = -0.181712059283214D1*t2*t3*t35-0.9085602964160698D-1*t2*t3*R
     &A/t43/t41/t40*(-0.2616653655616543D0*t9/t11/t10/RA-0.4280276067686
     &658D0*t19/t3/t20/t10-0.69444444568D-3*t29/t30/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t64 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RB**(1.D0/3.D0)
      t5 = PI**2
      t6 = t5**(1.D0/3.D0)
      t7 = t6**2
      t9 = 1/t7*GBB
      t10 = RB**2
      t11 = t3**2
      t18 = GBB**2
      t19 = 1/t6/t5*t18
      t20 = t10**2
      t26 = t5**2
      t29 = 1/t26*t18*GBB
      t30 = t20**2
      t35 = (1.D0+0.9812451208562037D-1*t9/t11/t10+0.8025517626912484D-1
     &*t19/t3/t20/RB+0.8680555571D-4*t29/t30)**(1.D0/15.D0)
      t40 = t35**2
      t41 = t40**2
      t43 = t41**2
      t64 = -0.181712059283214D1*t2*t3*t35-0.9085602964160698D-1*t2*t3*R
     &B/t43/t41/t40*(-0.2616653655616543D0*t9/t11/t10/RB-0.4280276067686
     &658D0*t19/t3/t20/t10-0.69444444568D-3*t29/t30/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t64 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RA**2
      t12 = t3**2
      t14 = 1/t12/t11
      t18 = 1/t7/t6
      t19 = GAA**2
      t21 = t11**2
      t24 = 1/t3/t21/RA
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t32 = 1/t31
      t36 = (1.D0+0.9812451208562037D-1*t9*GAA*t14+0.8025517626912484D-1
     &*t18*t19*t24+0.8680555571D-4*t28*t19*GAA*t32)**(1.D0/15.D0)
      t37 = t36**2
      t38 = t37**2
      t40 = t38**2
      t55 = -0.9085602964160698D-1*t2*t3*RA/t40/t38/t37*(0.9812451208562
     &037D-1*t9*t14+0.1605103525382497D0*t18*GAA*t24+0.26041666713D-3*t2
     &8*t19*t32)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t55 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RB**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RB**2
      t12 = t3**2
      t14 = 1/t12/t11
      t18 = 1/t7/t6
      t19 = GBB**2
      t21 = t11**2
      t24 = 1/t3/t21/RB
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t32 = 1/t31
      t36 = (1.D0+0.9812451208562037D-1*t9*GBB*t14+0.8025517626912484D-1
     &*t18*t19*t24+0.8680555571D-4*t28*t19*GBB*t32)**(1.D0/15.D0)
      t37 = t36**2
      t38 = t37**2
      t40 = t38**2
      t55 = -0.9085602964160698D-1*t2*t3*RB/t40/t38/t37*(0.9812451208562
     &037D-1*t9*t14+0.1605103525382497D0*t18*GBB*t24+0.26041666713D-3*t2
     &8*t19*t32)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t55 
      t1 = 0.D0
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t1 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t11 = RA**2
      t12 = t3**2
      t19 = GAA**2
      t21 = t11**2
      t27 = t6**2
      t31 = t21**2
      t36 = (1.D0+0.9812451208562037D-1/t8*GAA/t12/t11+0.802551762691248
     &4D-1/t7/t6*t19/t3/t21/RA+0.8680555571D-4/t27*t19*GAA/t31)**(1.D0/1
     &5.D0)
      t39 = -0.1362840444624105D1*t2*t3*RA*t36
      F(i) = F(i) +      t39 
      IF (IDERIV .GE. 1) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t5 = PI**2
      t6 = t5**(1.D0/3.D0)
      t7 = t6**2
      t9 = 1/t7*GAA
      t10 = RA**2
      t11 = t3**2
      t18 = GAA**2
      t19 = 1/t6/t5*t18
      t20 = t10**2
      t26 = t5**2
      t29 = 1/t26*t18*GAA
      t30 = t20**2
      t35 = (1.D0+0.9812451208562037D-1*t9/t11/t10+0.8025517626912484D-1
     &*t19/t3/t20/RA+0.8680555571D-4*t29/t30)**(1.D0/15.D0)
      t40 = t35**2
      t41 = t40**2
      t43 = t41**2
      t64 = -0.181712059283214D1*t2*t3*t35-0.9085602964160698D-1*t2*t3*R
     &A/t43/t41/t40*(-0.2616653655616543D0*t9/t11/t10/RA-0.4280276067686
     &658D0*t19/t3/t20/t10-0.69444444568D-3*t29/t30/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t64 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RA**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RA**2
      t12 = t3**2
      t14 = 1/t12/t11
      t18 = 1/t7/t6
      t19 = GAA**2
      t21 = t11**2
      t24 = 1/t3/t21/RA
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t32 = 1/t31
      t36 = (1.D0+0.9812451208562037D-1*t9*GAA*t14+0.8025517626912484D-1
     &*t18*t19*t24+0.8680555571D-4*t28*t19*GAA*t32)**(1.D0/15.D0)
      t37 = t36**2
      t38 = t37**2
      t40 = t38**2
      t55 = -0.9085602964160698D-1*t2*t3*RA/t40/t38/t37*(0.9812451208562
     &037D-1*t9*t14+0.1605103525382497D0*t18*GAA*t24+0.26041666713D-3*t2
     &8*t19*t32)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t55 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RB**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t11 = RB**2
      t12 = t3**2
      t19 = GBB**2
      t21 = t11**2
      t27 = t6**2
      t31 = t21**2
      t36 = (1.D0+0.9812451208562037D-1/t8*GBB/t12/t11+0.802551762691248
     &4D-1/t7/t6*t19/t3/t21/RB+0.8680555571D-4/t27*t19*GBB/t31)**(1.D0/1
     &5.D0)
      t39 = -0.1362840444624105D1*t2*t3*RB*t36
      F(i) = F(i) +      t39 
      IF (IDERIV .GE. 1) THEN
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RB**(1.D0/3.D0)
      t5 = PI**2
      t6 = t5**(1.D0/3.D0)
      t7 = t6**2
      t9 = 1/t7*GBB
      t10 = RB**2
      t11 = t3**2
      t18 = GBB**2
      t19 = 1/t6/t5*t18
      t20 = t10**2
      t26 = t5**2
      t29 = 1/t26*t18*GBB
      t30 = t20**2
      t35 = (1.D0+0.9812451208562037D-1*t9/t11/t10+0.8025517626912484D-1
     &*t19/t3/t20/RB+0.8680555571D-4*t29/t30)**(1.D0/15.D0)
      t40 = t35**2
      t41 = t40**2
      t43 = t41**2
      t64 = -0.181712059283214D1*t2*t3*t35-0.9085602964160698D-1*t2*t3*R
     &B/t43/t41/t40*(-0.2616653655616543D0*t9/t11/t10/RB-0.4280276067686
     &658D0*t19/t3/t20/t10-0.69444444568D-3*t29/t30/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t64 
      t2 = (1/PI)**(1.D0/3.D0)
      t3 = RB**(1.D0/3.D0)
      t6 = PI**2
      t7 = t6**(1.D0/3.D0)
      t8 = t7**2
      t9 = 1/t8
      t11 = RB**2
      t12 = t3**2
      t14 = 1/t12/t11
      t18 = 1/t7/t6
      t19 = GBB**2
      t21 = t11**2
      t24 = 1/t3/t21/RB
      t27 = t6**2
      t28 = 1/t27
      t31 = t21**2
      t32 = 1/t31
      t36 = (1.D0+0.9812451208562037D-1*t9*GBB*t14+0.8025517626912484D-1
     &*t18*t19*t24+0.8680555571D-4*t28*t19*GBB*t32)**(1.D0/15.D0)
      t37 = t36**2
      t38 = t37**2
      t40 = t38**2
      t55 = -0.9085602964160698D-1*t2*t3*RB/t40/t38/t37*(0.9812451208562
     &037D-1*t9*t14+0.1605103525382497D0*t18*GBB*t24+0.26041666713D-3*t2
     &8*t19*t32)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t55 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
