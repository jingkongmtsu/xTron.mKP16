!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! LYP evaluates the Lee-Yang-Parr correlation functional     
! and its derivatives.
! C. Lee, W. Yang, and R. G. Parr, 
! “Development of the Colle-Salvetti correlation-energy 
!  formula into a functional of the electron density,” 
! Phys. Rev. B, 37 (1988) 785-89.
! However, this is the revised version of LYP:
! B. Miehlich, A. Savin, H. Stoll and H. Preuss, 
! “Results obtained with the correlation
!  energy density functionals of becke and lee, 
!  yang and parr”, 
! Chem. Phys. Lett.,157(3):200–206, 1989
!
!



      subroutine functional_lyp
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
      CALL functional_lyp_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
      ELSE  
      CALL functional_lyp_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_lyp_close
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
      t1 = RA*RB
      t2 = RA+RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t7 = 1/(1.D0+0.349D0*t4)
      t8 = 1/t2
      t13 = dexp(-0.2533D0*t4)
      t14 = t13*t7
      t15 = t2**2
      t17 = t3**2
      t19 = 1/t17/t15/t2
      t20 = t14*t19
      t21 = PI**2
      t22 = t21**(1.D0/3.D0)
      t23 = t22**2
      t24 = RA**2
      t25 = RA**(1.D0/3.D0)
      t26 = t25**2
      t28 = RB**2
      t29 = RB**(1.D0/3.D0)
      t30 = t29**2
      t38 = t4*t7
      t42 = GAA+GBB+2.D0*GAB
      t72 = 0.6666666666666667D0*t15
      t85 = -0.19672D0*t1*t7*t8-0.5144476616948204D-1*t20*t1*t23*(t26*t2
     &4+t30*t28)-0.649176D-2*t20*t1*(0.2611111111111111D1-0.9850555556D-
     &1*t4-0.1357222222D0*t38)*t42+0.649176D-2*t20*t1*(0.25D1-0.14072222
     &22D-1*t4-0.1938888889D-1*t38)*(GAA+GBB)+0.649176D-2*t20*t1*(0.2814
     &444444D-1*t4+0.3877777778D-1*t38-0.1222222222222222D1)*(RA*GAA*t8+
     &RB*GBB*t8)+0.432784D-2*t14/t17/t2*t42-0.649176D-2*t14*t19*(t72-1.D
     &0*t24)*GBB-0.649176D-2*t14*t19*(t72-1.D0*t28)*GAA
      F(i) = F(i) +      t85 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t7 = 1.D0+0.349D0*t3
      t8 = 1/t7
      t9 = t5*t8
      t10 = t1**2
      t11 = t10*t1
      t12 = t2**2
      t14 = 1/t12/t11
      t16 = 0.1333333333333333D1*RB
      t22 = t10**2
      t24 = 1/t22/t1
      t25 = t24*t5
      t26 = 0.6666666666666667D0*t10
      t27 = RA**2
      t29 = t26-1.D0*t27
      t34 = t7**2
      t35 = 1/t34
      t36 = t5*t35
      t47 = RB**2
      t49 = t26-1.D0*t47
      t59 = 1/t12/t22
      t68 = RA*RB
      t69 = 1/t10
      t76 = GAA+GBB+2.D0*GAB
      t85 = 1/t11
      t94 = 1/t1
      t97 = t9*t14
      t99 = t3*t8
      t101 = 0.2611111111111111D1-0.9850555556D-1*t3-0.1357222222D0*t99
      t106 = PI**2
      t107 = t106**(1.D0/3.D0)
      t108 = t107**2
      t110 = RA**(1.D0/3.D0)
      t111 = t110**2
      t112 = t111*t27
      t113 = RB**(1.D0/3.D0)
      t114 = t113**2
      t116 = t114*t47+t112
      t124 = -0.649176D-2*t9*t14*(-0.6666666666666667D0*RA+t16)*GBB-0.54
     &8120936D-3*t25*t8*t29*GBB-0.75520808D-3*t36*t24*t29*GBB-0.649176D-
     &2*t9*t14*(0.1333333333333333D1*RA+t16)*GAA-0.548120936D-3*t25*t8*t
     &49*GAA-0.75520808D-3*t36*t24*t49*GAA+0.2380312D-1*t9*t59*t29*GBB+0
     &.2380312D-1*t9*t59*t49*GAA+0.19672D0*t68*t8*t69-0.7213066666666667
     &D-2*t9/t12/t10*t76-0.2288509333333333D-1*t68*t35/t2/t10+0.36541395
     &73333333D-3*t85*t5*t8*t76+0.5034720533333333D-3*t36*t85*t76-0.1967
     &2D0*RB*t8*t94-0.649176D-2*t97*RB*t101*t76-0.5144476616948204D-1*t9
     &7*RB*t108*t116-0.1371860431186188D0*t97*t112*RB*t108
      t127 = 0.25D1-0.1407222222D-1*t3-0.1938888889D-1*t99
      t129 = GAA+GBB
      t135 = 0.2814444444D-1*t3+0.3877777778D-1*t99-0.1222222222222222D1
      t137 = RA*GAA
      t139 = RB*GBB
      t141 = t137*t94+t139*t94
      t146 = 1/t2/t1
      t148 = t146*t8
      t152 = 1/t12/t1*t35
      t159 = t25*t8
      t161 = t68*t101*t76
      t164 = t36*t24
      t168 = t68*t108*t116
      t173 = t9*t59
      t179 = t68*t127*t129
      t183 = t68*t135*t141
      s1 = 0.649176D-2*t97*RB*t127*t129+0.649176D-2*t97*RB*t135*t141-0.6
     &49176D-2*t97*t68*(0.3283518518666667D-1*t146+0.4524074073333333D-1
     &*t148-0.1578901851593333D-1*t152)*t76-0.548120936D-3*t159*t161-0.7
     &5520808D-3*t164*t161-0.4343653090243267D-2*t159*t168-0.59847411310
     &49744D-2*t164*t168+0.1886308092881008D0*t173*t168+0.2380312D-1*t17
     &3*t161
      t220 = s1-0.2380312D-1*t173*t179-0.2380312D-1*t173*t183+0.649176D-
     &2*t97*t68*(0.469074074D-2*t146+0.6462962963333333D-2*t148-0.225557
     &4074203333D-2*t152)*t129+0.548120936D-3*t159*t179+0.75520808D-3*t1
     &64*t179+0.649176D-2*t97*t68*(-0.938148148D-2*t146-0.12925925926666
     &67D-1*t148+0.4511148148406667D-2*t152)*t141+0.649176D-2*t97*t68*t1
     &35*(GAA*t94-1.D0*t137*t69-1.D0*t139*t69)+0.548120936D-3*t159*t183+
     &0.75520808D-3*t164*t183
      t221 = t124+t220
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t221 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t8 = 1/(1.D0+0.349D0*t3)
      t9 = t5*t8
      t10 = t1**2
      t12 = t2**2
      t14 = 1/t12/t10/t1
      t15 = t9*t14
      t16 = RA*RB
      t18 = t3*t8
      t30 = t10**2
      t34 = RA**2
      t47 = RB**2
      t53 = -0.649176D-2*t15*t16*(0.2611111111111111D1-0.9850555556D-1*t
     &3-0.1357222222D0*t18)+0.649176D-2*t15*t16*(0.25D1-0.1407222222D-1*
     &t3-0.1938888889D-1*t18)+0.649176D-2*t9/t12/t30*t34*RB*(0.281444444
     &4D-1*t3+0.3877777778D-1*t18-0.1222222222222222D1)+0.432784D-2*t9/t
     &12/t1-0.649176D-2*t9*t14*(0.6666666666666667D0*t10-1.D0*t47)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t53 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t8 = 1/(1.D0+0.349D0*t3)
      t9 = t5*t8
      t10 = t1**2
      t12 = t2**2
      t28 = -0.1298352D-1*t9/t12/t10/t1*RA*RB*(0.2611111111111111D1-0.98
     &50555556D-1*t3-0.1357222222D0*t3*t8)+0.865568D-2*t9/t12/t1
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t28 
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
                   
                   
                   
                   
      subroutine functional_lyp_open
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
      t1 = RA*RB
      t2 = RA+RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t7 = 1/(1.D0+0.349D0*t4)
      t8 = 1/t2
      t13 = dexp(-0.2533D0*t4)
      t14 = t13*t7
      t15 = t2**2
      t17 = t3**2
      t19 = 1/t17/t15/t2
      t20 = t14*t19
      t21 = PI**2
      t22 = t21**(1.D0/3.D0)
      t23 = t22**2
      t24 = RA**2
      t25 = RA**(1.D0/3.D0)
      t26 = t25**2
      t28 = RB**2
      t29 = RB**(1.D0/3.D0)
      t30 = t29**2
      t38 = t4*t7
      t42 = GAA+GBB+2.D0*GAB
      t72 = 0.6666666666666667D0*t15
      t85 = -0.19672D0*t1*t7*t8-0.5144476616948204D-1*t20*t1*t23*(t26*t2
     &4+t30*t28)-0.649176D-2*t20*t1*(0.2611111111111111D1-0.9850555556D-
     &1*t4-0.1357222222D0*t38)*t42+0.649176D-2*t20*t1*(0.25D1-0.14072222
     &22D-1*t4-0.1938888889D-1*t38)*(GAA+GBB)+0.649176D-2*t20*t1*(0.2814
     &444444D-1*t4+0.3877777778D-1*t38-0.1222222222222222D1)*(RA*GAA*t8+
     &RB*GBB*t8)+0.432784D-2*t14/t17/t2*t42-0.649176D-2*t14*t19*(t72-1.D
     &0*t24)*GBB-0.649176D-2*t14*t19*(t72-1.D0*t28)*GAA
      F(i) = F(i) +      t85 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t7 = 1.D0+0.349D0*t3
      t8 = 1/t7
      t9 = t5*t8
      t10 = t1**2
      t11 = t10*t1
      t12 = t2**2
      t14 = 1/t12/t11
      t16 = 0.1333333333333333D1*RB
      t22 = t10**2
      t24 = 1/t22/t1
      t25 = t24*t5
      t26 = 0.6666666666666667D0*t10
      t27 = RA**2
      t29 = t26-1.D0*t27
      t34 = t7**2
      t35 = 1/t34
      t36 = t5*t35
      t47 = RB**2
      t49 = t26-1.D0*t47
      t59 = 1/t12/t22
      t68 = RA*RB
      t69 = 1/t10
      t76 = GAA+GBB+2.D0*GAB
      t85 = 1/t11
      t94 = 1/t1
      t97 = t9*t14
      t99 = t3*t8
      t101 = 0.2611111111111111D1-0.9850555556D-1*t3-0.1357222222D0*t99
      t106 = PI**2
      t107 = t106**(1.D0/3.D0)
      t108 = t107**2
      t110 = RA**(1.D0/3.D0)
      t111 = t110**2
      t112 = t111*t27
      t113 = RB**(1.D0/3.D0)
      t114 = t113**2
      t116 = t114*t47+t112
      t124 = -0.649176D-2*t9*t14*(-0.6666666666666667D0*RA+t16)*GBB-0.54
     &8120936D-3*t25*t8*t29*GBB-0.75520808D-3*t36*t24*t29*GBB-0.649176D-
     &2*t9*t14*(0.1333333333333333D1*RA+t16)*GAA-0.548120936D-3*t25*t8*t
     &49*GAA-0.75520808D-3*t36*t24*t49*GAA+0.2380312D-1*t9*t59*t29*GBB+0
     &.2380312D-1*t9*t59*t49*GAA+0.19672D0*t68*t8*t69-0.7213066666666667
     &D-2*t9/t12/t10*t76-0.2288509333333333D-1*t68*t35/t2/t10+0.36541395
     &73333333D-3*t85*t5*t8*t76+0.5034720533333333D-3*t36*t85*t76-0.1967
     &2D0*RB*t8*t94-0.649176D-2*t97*RB*t101*t76-0.5144476616948204D-1*t9
     &7*RB*t108*t116-0.1371860431186188D0*t97*t112*RB*t108
      t127 = 0.25D1-0.1407222222D-1*t3-0.1938888889D-1*t99
      t129 = GAA+GBB
      t135 = 0.2814444444D-1*t3+0.3877777778D-1*t99-0.1222222222222222D1
      t137 = RA*GAA
      t139 = RB*GBB
      t141 = t137*t94+t139*t94
      t146 = 1/t2/t1
      t148 = t146*t8
      t152 = 1/t12/t1*t35
      t159 = t25*t8
      t161 = t68*t101*t76
      t164 = t36*t24
      t168 = t68*t108*t116
      t173 = t9*t59
      t179 = t68*t127*t129
      t183 = t68*t135*t141
      s1 = 0.649176D-2*t97*RB*t127*t129+0.649176D-2*t97*RB*t135*t141-0.6
     &49176D-2*t97*t68*(0.3283518518666667D-1*t146+0.4524074073333333D-1
     &*t148-0.1578901851593333D-1*t152)*t76-0.548120936D-3*t159*t161-0.7
     &5520808D-3*t164*t161-0.4343653090243267D-2*t159*t168-0.59847411310
     &49744D-2*t164*t168+0.1886308092881008D0*t173*t168+0.2380312D-1*t17
     &3*t161
      t220 = s1-0.2380312D-1*t173*t179-0.2380312D-1*t173*t183+0.649176D-
     &2*t97*t68*(0.469074074D-2*t146+0.6462962963333333D-2*t148-0.225557
     &4074203333D-2*t152)*t129+0.548120936D-3*t159*t179+0.75520808D-3*t1
     &64*t179+0.649176D-2*t97*t68*(-0.938148148D-2*t146-0.12925925926666
     &67D-1*t148+0.4511148148406667D-2*t152)*t141+0.649176D-2*t97*t68*t1
     &35*(GAA*t94-1.D0*t137*t69-1.D0*t139*t69)+0.548120936D-3*t159*t183+
     &0.75520808D-3*t164*t183
      t221 = t124+t220
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t221 
      t1 = RA+RB
      t2 = t1**2
      t3 = t2**2
      t5 = 1/t3/t1
      t6 = t1**(1.D0/3.D0)
      t7 = 1/t6
      t9 = dexp(-0.2533D0*t7)
      t10 = t5*t9
      t12 = 1.D0+0.349D0*t7
      t13 = 1/t12
      t14 = 0.6666666666666667D0*t2
      t15 = RA**2
      t17 = t14-1.D0*t15
      t22 = t12**2
      t23 = 1/t22
      t24 = t9*t23
      t29 = RB**2
      t31 = t14-1.D0*t29
      t40 = t9*t13
      t41 = t6**2
      t43 = 1/t41/t3
      t52 = t2*t1
      t54 = 1/t41/t52
      t55 = 0.1333333333333333D1*RA
      t68 = RA*RB
      t69 = 1/t2
      t76 = GAA+GBB+2.D0*GAB
      t85 = 1/t52
      t93 = t40*t54
      t94 = PI**2
      t95 = t94**(1.D0/3.D0)
      t96 = t95**2
      t98 = RA**(1.D0/3.D0)
      t99 = t98**2
      t101 = RB**(1.D0/3.D0)
      t102 = t101**2
      t103 = t102*t29
      t104 = t99*t15+t103
      t109 = t7*t13
      t111 = 0.2611111111111111D1-0.9850555556D-1*t7-0.1357222222D0*t109
      t118 = 0.25D1-0.1407222222D-1*t7-0.1938888889D-1*t109
      t120 = GAA+GBB
      t126 = 0.2814444444D-1*t7+0.3877777778D-1*t109-0.1222222222222222D
     &1
      t128 = RA*GAA
      t129 = 1/t1
      t131 = RB*GBB
      t133 = t128*t129+t131*t129
      t137 = -0.548120936D-3*t10*t13*t17*GBB-0.75520808D-3*t24*t5*t17*GB
     &B-0.548120936D-3*t10*t13*t31*GAA-0.75520808D-3*t24*t5*t31*GAA+0.23
     &80312D-1*t40*t43*t17*GBB+0.2380312D-1*t40*t43*t31*GAA-0.649176D-2*
     &t40*t54*(t55+0.1333333333333333D1*RB)*GBB-0.649176D-2*t40*t54*(t55
     &-0.6666666666666667D0*RB)*GAA+0.19672D0*t68*t13*t69-0.721306666666
     &6667D-2*t40/t41/t2*t76-0.2288509333333333D-1*t68*t23/t6/t2+0.36541
     &39573333333D-3*t85*t9*t13*t76+0.5034720533333333D-3*t24*t85*t76-0.
     &5144476616948204D-1*t93*RA*t96*t104-0.649176D-2*t93*RA*t111*t76+0.
     &649176D-2*t93*RA*t118*t120+0.649176D-2*t93*RA*t126*t133
      t143 = 1/t6/t1
      t145 = t143*t13
      t149 = 1/t41/t1*t23
      t156 = t10*t13
      t158 = t68*t111*t76
      t161 = t24*t5
      t165 = t68*t96*t104
      t170 = t40*t43
      t176 = t68*t118*t120
      t180 = t68*t126*t133
      s1 = -0.1371860431186188D0*t93*RA*t103*t96-0.649176D-2*t93*t68*(0.
     &3283518518666667D-1*t143+0.4524074073333333D-1*t145-0.157890185159
     &3333D-1*t149)*t76-0.548120936D-3*t156*t158-0.75520808D-3*t161*t158
     &-0.4343653090243267D-2*t156*t165-0.5984741131049744D-2*t161*t165+0
     &.1886308092881008D0*t170*t165+0.2380312D-1*t170*t158-0.2380312D-1*
     &t170*t176
      t220 = s1-0.2380312D-1*t170*t180+0.649176D-2*t93*t68*(0.469074074D
     &-2*t143+0.6462962963333333D-2*t145-0.2255574074203333D-2*t149)*t12
     &0+0.548120936D-3*t156*t176+0.75520808D-3*t161*t176+0.649176D-2*t93
     &*t68*(-0.938148148D-2*t143-0.1292592592666667D-1*t145+0.4511148148
     &406667D-2*t149)*t133+0.548120936D-3*t156*t180+0.75520808D-3*t161*t
     &180+0.649176D-2*t93*t68*t126*(-1.D0*t128*t69+GBB*t129-1.D0*t131*t6
     &9)-0.19672D0*RA*t13*t129
      t221 = t137+t220
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t221 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t8 = 1/(1.D0+0.349D0*t3)
      t9 = t5*t8
      t10 = t1**2
      t12 = t2**2
      t14 = 1/t12/t10/t1
      t15 = t9*t14
      t16 = RA*RB
      t18 = t3*t8
      t30 = t10**2
      t34 = RA**2
      t47 = RB**2
      t53 = -0.649176D-2*t15*t16*(0.2611111111111111D1-0.9850555556D-1*t
     &3-0.1357222222D0*t18)+0.649176D-2*t15*t16*(0.25D1-0.1407222222D-1*
     &t3-0.1938888889D-1*t18)+0.649176D-2*t9/t12/t30*t34*RB*(0.281444444
     &4D-1*t3+0.3877777778D-1*t18-0.1222222222222222D1)+0.432784D-2*t9/t
     &12/t1-0.649176D-2*t9*t14*(0.6666666666666667D0*t10-1.D0*t47)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t53 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t8 = 1/(1.D0+0.349D0*t3)
      t9 = t5*t8
      t10 = t1**2
      t12 = t2**2
      t14 = 1/t12/t10/t1
      t15 = t9*t14
      t16 = RA*RB
      t18 = t3*t8
      t30 = t10**2
      t34 = RB**2
      t47 = RA**2
      t53 = -0.649176D-2*t15*t16*(0.2611111111111111D1-0.9850555556D-1*t
     &3-0.1357222222D0*t18)+0.649176D-2*t15*t16*(0.25D1-0.1407222222D-1*
     &t3-0.1938888889D-1*t18)+0.649176D-2*t9/t12/t30*RA*t34*(0.281444444
     &4D-1*t3+0.3877777778D-1*t18-0.1222222222222222D1)+0.432784D-2*t9/t
     &12/t1-0.649176D-2*t9*t14*(0.6666666666666667D0*t10-1.D0*t47)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t53 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t8 = 1/(1.D0+0.349D0*t3)
      t9 = t5*t8
      t10 = t1**2
      t12 = t2**2
      t28 = -0.1298352D-1*t9/t12/t10/t1*RA*RB*(0.2611111111111111D1-0.98
     &50555556D-1*t3-0.1357222222D0*t3*t8)+0.865568D-2*t9/t12/t1
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t28 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t9 = t5/(1.D0+0.349D0*t3)
      t10 = t2**2
      t18 = t1**2
      t23 = RB**2
      t30 = 0.432784D-2*t9/t10/t1*(GAA+GBB+2.D0*GAB)-0.649176D-2*t9/t10/
     &t18/t1*(0.6666666666666667D0*t18-1.D0*t23)*GAA
      F(i) = F(i) +      t30 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = t1**2
      t3 = t2*t1
      t4 = 1/t3
      t5 = t1**(1.D0/3.D0)
      t6 = 1/t5
      t8 = dexp(-0.2533D0*t6)
      t11 = 1.D0+0.349D0*t6
      t12 = 1/t11
      t14 = GAA+GBB+2.D0*GAB
      t18 = t11**2
      t20 = t8/t18
      t24 = t8*t12
      t25 = t5**2
      t31 = t2**2
      t33 = 1/t31/t1
      t36 = RB**2
      t38 = 0.6666666666666667D0*t2-1.D0*t36
      t62 = 0.3654139573333333D-3*t4*t8*t12*t14+0.5034720533333333D-3*t2
     &0*t4*t14-0.7213066666666667D-2*t24/t25/t2*t14-0.548120936D-3*t33*t
     &8*t12*t38*GAA-0.75520808D-3*t20*t33*t38*GAA+0.2380312D-1*t24/t25/t
     &31*t38*GAA-0.649176D-2*t24/t25/t3*(0.1333333333333333D1*RA+0.13333
     &33333333333D1*RB)*GAA
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t62 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t9 = t5/(1.D0+0.349D0*t3)
      t10 = t2**2
      t15 = t1**2
      t20 = RB**2
      t26 = 0.432784D-2*t9/t10/t1-0.649176D-2*t9/t10/t15/t1*(0.666666666
     &6666667D0*t15-1.D0*t20)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t26 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t9 = t5/(1.D0+0.349D0*t3)
      t10 = t2**2
      t18 = t1**2
      t23 = RA**2
      t30 = 0.432784D-2*t9/t10/t1*(GAA+GBB+2.D0*GAB)-0.649176D-2*t9/t10/
     &t18/t1*(0.6666666666666667D0*t18-1.D0*t23)*GBB
      F(i) = F(i) +      t30 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = t1**2
      t3 = t2*t1
      t4 = 1/t3
      t5 = t1**(1.D0/3.D0)
      t6 = 1/t5
      t8 = dexp(-0.2533D0*t6)
      t11 = 1.D0+0.349D0*t6
      t12 = 1/t11
      t14 = GAA+GBB+2.D0*GAB
      t18 = t11**2
      t20 = t8/t18
      t24 = t8*t12
      t25 = t5**2
      t31 = t2**2
      t33 = 1/t31/t1
      t36 = RA**2
      t38 = 0.6666666666666667D0*t2-1.D0*t36
      t62 = 0.3654139573333333D-3*t4*t8*t12*t14+0.5034720533333333D-3*t2
     &0*t4*t14-0.7213066666666667D-2*t24/t25/t2*t14-0.548120936D-3*t33*t
     &8*t12*t38*GBB-0.75520808D-3*t20*t33*t38*GBB+0.2380312D-1*t24/t25/t
     &31*t38*GBB-0.649176D-2*t24/t25/t3*(0.1333333333333333D1*RA+0.13333
     &33333333333D1*RB)*GBB
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t62 
      t1 = RA+RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = dexp(-0.2533D0*t3)
      t9 = t5/(1.D0+0.349D0*t3)
      t10 = t2**2
      t15 = t1**2
      t20 = RA**2
      t26 = 0.432784D-2*t9/t10/t1-0.649176D-2*t9/t10/t15/t1*(0.666666666
     &6666667D0*t15-1.D0*t20)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t26 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
