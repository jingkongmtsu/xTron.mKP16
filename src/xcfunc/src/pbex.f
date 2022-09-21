!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! J. P. Perdew, K. Burke, and M. Ernzerhof, 
! “Generalized gradient approximation made simple,” 
! Phys. Rev. Lett., 77 (1996) 3865-68.
!
!



      subroutine functional_pbex
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
      CALL functional_pbex_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
      ELSE  
      CALL functional_pbex_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_pbex_close
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
      t1 = PI**2
      t3 = (t1*RA)**(1.D0/3.D0)
      t5 = 1/PI
      t6 = t3**2
      t9 = RA**2
      t21 = (t1*RB)**(1.D0/3.D0)
      t23 = t21**2
      t26 = RB**2
      t37 = -0.1362840444624105D1*RA*t3*t5*(0.1804D1-0.804D0/(1.D0+0.206
     &7191009821618D-1*GAA/t6/t9))-0.1362840444624105D1*RB*t21*t5*(0.180
     &4D1-0.804D0/(1.D0+0.2067191009821618D-1*GBB/t23/t26))
      F(i) = F(i) +      t37 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = t3**2
      t7 = 1/t6
      t8 = GAA*t7
      t9 = RA**2
      t10 = 1/t9
      t13 = 1.D0+0.2067191009821618D-1*t8*t10
      t16 = 0.1804D1-0.804D0/t13
      t24 = t13**2
      t41 = -0.1362840444624105D1*t3*t4*t16-0.4542801482080349D0*RA*t7*P
     &I*t16-0.109572371747778D1*RA*t3*t4/t24*(-0.1378127339881079D-1*GAA
     &/t6/t2*t10*t1-0.4134382019643236D-1*t8/t9/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t41 
      t2 = PI**2
      t4 = (t2*RA)**(1.D0/3.D0)
      t8 = t4**2
      t11 = RA**2
      t16 = (1.D0+0.2067191009821618D-1*GAA/t8/t11)**2
      t21 = -0.226507021801839D-1/RA/t4/PI/t16
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t21 
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
                   
                   
                   
                   
      subroutine functional_pbex_open
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
      t1 = PI**2
      t3 = (t1*RA)**(1.D0/3.D0)
      t5 = 1/PI
      t6 = t3**2
      t9 = RA**2
      t21 = (t1*RB)**(1.D0/3.D0)
      t23 = t21**2
      t26 = RB**2
      t37 = -0.1362840444624105D1*RA*t3*t5*(0.1804D1-0.804D0/(1.D0+0.206
     &7191009821618D-1*GAA/t6/t9))-0.1362840444624105D1*RB*t21*t5*(0.180
     &4D1-0.804D0/(1.D0+0.2067191009821618D-1*GBB/t23/t26))
      F(i) = F(i) +      t37 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = t3**2
      t7 = 1/t6
      t8 = GAA*t7
      t9 = RA**2
      t10 = 1/t9
      t13 = 1.D0+0.2067191009821618D-1*t8*t10
      t16 = 0.1804D1-0.804D0/t13
      t24 = t13**2
      t41 = -0.1362840444624105D1*t3*t4*t16-0.4542801482080349D0*RA*t7*P
     &I*t16-0.109572371747778D1*RA*t3*t4/t24*(-0.1378127339881079D-1*GAA
     &/t6/t2*t10*t1-0.4134382019643236D-1*t8/t9/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t41 
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = t3**2
      t7 = 1/t6
      t8 = GBB*t7
      t9 = RB**2
      t10 = 1/t9
      t13 = 1.D0+0.2067191009821618D-1*t8*t10
      t16 = 0.1804D1-0.804D0/t13
      t24 = t13**2
      t41 = -0.1362840444624105D1*t3*t4*t16-0.4542801482080349D0*RB*t7*P
     &I*t16-0.109572371747778D1*RB*t3*t4/t24*(-0.1378127339881079D-1*GBB
     &/t6/t2*t10*t1-0.4134382019643236D-1*t8/t9/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t41 
      t2 = PI**2
      t4 = (t2*RA)**(1.D0/3.D0)
      t8 = t4**2
      t11 = RA**2
      t16 = (1.D0+0.2067191009821618D-1*GAA/t8/t11)**2
      t21 = -0.226507021801839D-1/RA/t4/PI/t16
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t21 
      t2 = PI**2
      t4 = (t2*RB)**(1.D0/3.D0)
      t8 = t4**2
      t11 = RB**2
      t16 = (1.D0+0.2067191009821618D-1*GBB/t8/t11)**2
      t21 = -0.226507021801839D-1/RB/t4/PI/t16
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t21 
      t1 = 0.D0
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t1 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = PI**2
      t3 = (t1*RA)**(1.D0/3.D0)
      t6 = t3**2
      t9 = RA**2
      t20 = -0.1362840444624105D1*RA*t3/PI*(0.1804D1-0.804D0/(1.D0+0.206
     &7191009821618D-1*GAA/t6/t9))
      F(i) = F(i) +      t20 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = t3**2
      t7 = 1/t6
      t8 = GAA*t7
      t9 = RA**2
      t10 = 1/t9
      t13 = 1.D0+0.2067191009821618D-1*t8*t10
      t16 = 0.1804D1-0.804D0/t13
      t24 = t13**2
      t41 = -0.1362840444624105D1*t3*t4*t16-0.4542801482080349D0*RA*t7*P
     &I*t16-0.109572371747778D1*RA*t3*t4/t24*(-0.1378127339881079D-1*GAA
     &/t6/t2*t10*t1-0.4134382019643236D-1*t8/t9/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t41 
      t2 = PI**2
      t4 = (t2*RA)**(1.D0/3.D0)
      t8 = t4**2
      t11 = RA**2
      t16 = (1.D0+0.2067191009821618D-1*GAA/t8/t11)**2
      t21 = -0.226507021801839D-1/RA/t4/PI/t16
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t21 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = PI**2
      t3 = (t1*RB)**(1.D0/3.D0)
      t6 = t3**2
      t9 = RB**2
      t20 = -0.1362840444624105D1*RB*t3/PI*(0.1804D1-0.804D0/(1.D0+0.206
     &7191009821618D-1*GBB/t6/t9))
      F(i) = F(i) +      t20 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = t3**2
      t7 = 1/t6
      t8 = GBB*t7
      t9 = RB**2
      t10 = 1/t9
      t13 = 1.D0+0.2067191009821618D-1*t8*t10
      t16 = 0.1804D1-0.804D0/t13
      t24 = t13**2
      t41 = -0.1362840444624105D1*t3*t4*t16-0.4542801482080349D0*RB*t7*P
     &I*t16-0.109572371747778D1*RB*t3*t4/t24*(-0.1378127339881079D-1*GBB
     &/t6/t2*t10*t1-0.4134382019643236D-1*t8/t9/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t41 
      t2 = PI**2
      t4 = (t2*RB)**(1.D0/3.D0)
      t8 = t4**2
      t11 = RB**2
      t16 = (1.D0+0.2067191009821618D-1*GBB/t8/t11)**2
      t21 = -0.226507021801839D-1/RB/t4/PI/t16
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t21 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
