!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! A. D. Becke, 
! “Density-functional exchange-energy approximation 
!  with correct asymptotic-behavior,” 
! Phys. Rev. A, 38 (1988) 3098.
!
!



      subroutine functional_becke88
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
      CALL functional_becke88_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
      ELSE  
      CALL functional_becke88_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_becke88_close
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
      t1 = RA**(1.D0/3.D0)
      t2 = t1*RA
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = 1/t2
      t9 = dsqrt(GAA)
      t10 = t9*t7
      t11 = RA**2
      t12 = t1**2
      t17 = dsqrt(GAA/t12/t11+1)
      t19 = dlog(t10+t17)
      t26 = RB**(1.D0/3.D0)
      t27 = t26*RB
      t30 = 1/t27
      t32 = dsqrt(GBB)
      t33 = t32*t30
      t34 = RB**2
      t35 = t26**2
      t40 = dsqrt(GBB/t35/t34+1)
      t42 = dlog(t33+t40)
      t49 = -0.1362840444624105D1*t2*t4-0.42D-2*t7*GAA/(1.D0+0.252D-1*t1
     &0*t19)-0.1362840444624105D1*t27*t4-0.42D-2*t30*GBB/(1.D0+0.252D-1*
     &t33*t42)
      F(i) = F(i) +      t49 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = RA**2
      t8 = 1/t1/t6
      t10 = dsqrt(GAA)
      t12 = 1/t1/RA
      t13 = t10*t12
      t14 = t1**2
      t17 = GAA/t14/t6
      t19 = dsqrt(t17+1)
      t21 = dlog(t13+t19)
      t24 = 1.D0+0.252D-1*t13*t21
      t29 = t24**2
      t39 = dsqrt(t17+1.D0)
      t47 = -0.181712059283214D1*t1*t3+0.56D-2*t8*GAA/t24+0.42D-2*t12*GA
     &A/t29*(-0.336D-1*t10*t8*t21-0.336D-1*GAA/t14/t6/RA/t39)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t47 
      t1 = RA**(1.D0/3.D0)
      t3 = 1/t1/RA
      t4 = dsqrt(GAA)
      t5 = t4*t3
      t6 = RA**2
      t7 = t1**2
      t9 = 1/t7/t6
      t10 = GAA*t9
      t12 = dsqrt(t10+1)
      t14 = dlog(t5+t12)
      t17 = 1.D0+0.252D-1*t5*t14
      t22 = t17**2
      t29 = dsqrt(t10+1.D0)
      t37 = -0.42D-2*t3/t17+0.42D-2*t3*GAA/t22*(0.126D-1/t4*t3*t14+0.126
     &D-1*t9/t29)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t37 
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
                   
                   
                   
                   
      subroutine functional_becke88_open
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
      t1 = RA**(1.D0/3.D0)
      t2 = t1*RA
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = 1/t2
      t9 = dsqrt(GAA)
      t10 = t9*t7
      t11 = RA**2
      t12 = t1**2
      t17 = dsqrt(GAA/t12/t11+1)
      t19 = dlog(t10+t17)
      t26 = RB**(1.D0/3.D0)
      t27 = t26*RB
      t30 = 1/t27
      t32 = dsqrt(GBB)
      t33 = t32*t30
      t34 = RB**2
      t35 = t26**2
      t40 = dsqrt(GBB/t35/t34+1)
      t42 = dlog(t33+t40)
      t49 = -0.1362840444624105D1*t2*t4-0.42D-2*t7*GAA/(1.D0+0.252D-1*t1
     &0*t19)-0.1362840444624105D1*t27*t4-0.42D-2*t30*GBB/(1.D0+0.252D-1*
     &t33*t42)
      F(i) = F(i) +      t49 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = RA**2
      t8 = 1/t1/t6
      t10 = dsqrt(GAA)
      t12 = 1/t1/RA
      t13 = t10*t12
      t14 = t1**2
      t17 = GAA/t14/t6
      t19 = dsqrt(t17+1)
      t21 = dlog(t13+t19)
      t24 = 1.D0+0.252D-1*t13*t21
      t29 = t24**2
      t39 = dsqrt(t17+1.D0)
      t47 = -0.181712059283214D1*t1*t3+0.56D-2*t8*GAA/t24+0.42D-2*t12*GA
     &A/t29*(-0.336D-1*t10*t8*t21-0.336D-1*GAA/t14/t6/RA/t39)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t47 
      t1 = RB**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = RB**2
      t8 = 1/t1/t6
      t10 = dsqrt(GBB)
      t12 = 1/t1/RB
      t13 = t10*t12
      t14 = t1**2
      t17 = GBB/t14/t6
      t19 = dsqrt(t17+1)
      t21 = dlog(t13+t19)
      t24 = 1.D0+0.252D-1*t13*t21
      t29 = t24**2
      t39 = dsqrt(t17+1.D0)
      t47 = -0.181712059283214D1*t1*t3+0.56D-2*t8*GBB/t24+0.42D-2*t12*GB
     &B/t29*(-0.336D-1*t10*t8*t21-0.336D-1*GBB/t14/t6/RB/t39)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t47 
      t1 = RA**(1.D0/3.D0)
      t3 = 1/t1/RA
      t4 = dsqrt(GAA)
      t5 = t4*t3
      t6 = RA**2
      t7 = t1**2
      t9 = 1/t7/t6
      t10 = GAA*t9
      t12 = dsqrt(t10+1)
      t14 = dlog(t5+t12)
      t17 = 1.D0+0.252D-1*t5*t14
      t22 = t17**2
      t29 = dsqrt(t10+1.D0)
      t37 = -0.42D-2*t3/t17+0.42D-2*t3*GAA/t22*(0.126D-1/t4*t3*t14+0.126
     &D-1*t9/t29)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t37 
      t1 = RB**(1.D0/3.D0)
      t3 = 1/t1/RB
      t4 = dsqrt(GBB)
      t5 = t4*t3
      t6 = RB**2
      t7 = t1**2
      t9 = 1/t7/t6
      t10 = GBB*t9
      t12 = dsqrt(t10+1)
      t14 = dlog(t5+t12)
      t17 = 1.D0+0.252D-1*t5*t14
      t22 = t17**2
      t29 = dsqrt(t10+1.D0)
      t37 = -0.42D-2*t3/t17+0.42D-2*t3*GBB/t22*(0.126D-1/t4*t3*t14+0.126
     &D-1*t9/t29)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t37 
      t1 = 0.D0
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t1 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = RA**(1.D0/3.D0)
      t2 = t1*RA
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = 1/t2
      t9 = dsqrt(GAA)
      t10 = t9*t7
      t11 = RA**2
      t12 = t1**2
      t17 = dsqrt(GAA/t12/t11+1)
      t19 = dlog(t10+t17)
      t26 = -0.1362840444624105D1*t2*t4-0.42D-2*t7*GAA/(1.D0+0.252D-1*t1
     &0*t19)
      F(i) = F(i) +      t26 
      IF (IDERIV .GE. 1) THEN
      t1 = RA**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = RA**2
      t8 = 1/t1/t6
      t10 = dsqrt(GAA)
      t12 = 1/t1/RA
      t13 = t10*t12
      t14 = t1**2
      t17 = GAA/t14/t6
      t19 = dsqrt(t17+1)
      t21 = dlog(t13+t19)
      t24 = 1.D0+0.252D-1*t13*t21
      t29 = t24**2
      t39 = dsqrt(t17+1.D0)
      t47 = -0.181712059283214D1*t1*t3+0.56D-2*t8*GAA/t24+0.42D-2*t12*GA
     &A/t29*(-0.336D-1*t10*t8*t21-0.336D-1*GAA/t14/t6/RA/t39)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t47 
      t1 = RA**(1.D0/3.D0)
      t3 = 1/t1/RA
      t4 = dsqrt(GAA)
      t5 = t4*t3
      t6 = RA**2
      t7 = t1**2
      t9 = 1/t7/t6
      t10 = GAA*t9
      t12 = dsqrt(t10+1)
      t14 = dlog(t5+t12)
      t17 = 1.D0+0.252D-1*t5*t14
      t22 = t17**2
      t29 = dsqrt(t10+1.D0)
      t37 = -0.42D-2*t3/t17+0.42D-2*t3*GAA/t22*(0.126D-1/t4*t3*t14+0.126
     &D-1*t9/t29)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t37 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = RB**(1.D0/3.D0)
      t2 = t1*RB
      t4 = (1/PI)**(1.D0/3.D0)
      t7 = 1/t2
      t9 = dsqrt(GBB)
      t10 = t9*t7
      t11 = RB**2
      t12 = t1**2
      t17 = dsqrt(GBB/t12/t11+1)
      t19 = dlog(t10+t17)
      t26 = -0.1362840444624105D1*t2*t4-0.42D-2*t7*GBB/(1.D0+0.252D-1*t1
     &0*t19)
      F(i) = F(i) +      t26 
      IF (IDERIV .GE. 1) THEN
      t1 = RB**(1.D0/3.D0)
      t3 = (1/PI)**(1.D0/3.D0)
      t6 = RB**2
      t8 = 1/t1/t6
      t10 = dsqrt(GBB)
      t12 = 1/t1/RB
      t13 = t10*t12
      t14 = t1**2
      t17 = GBB/t14/t6
      t19 = dsqrt(t17+1)
      t21 = dlog(t13+t19)
      t24 = 1.D0+0.252D-1*t13*t21
      t29 = t24**2
      t39 = dsqrt(t17+1.D0)
      t47 = -0.181712059283214D1*t1*t3+0.56D-2*t8*GBB/t24+0.42D-2*t12*GB
     &B/t29*(-0.336D-1*t10*t8*t21-0.336D-1*GBB/t14/t6/RB/t39)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t47 
      t1 = RB**(1.D0/3.D0)
      t3 = 1/t1/RB
      t4 = dsqrt(GBB)
      t5 = t4*t3
      t6 = RB**2
      t7 = t1**2
      t9 = 1/t7/t6
      t10 = GBB*t9
      t12 = dsqrt(t10+1)
      t14 = dlog(t5+t12)
      t17 = 1.D0+0.252D-1*t5*t14
      t22 = t17**2
      t29 = dsqrt(t10+1.D0)
      t37 = -0.42D-2*t3/t17+0.42D-2*t3*GBB/t22*(0.126D-1/t4*t3*t14+0.126
     &D-1*t9/t29)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t37 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
