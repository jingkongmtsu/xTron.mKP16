!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! S. H. Vosko, L. Wilk, and M. Nusair, 
! “Accurate spin-dependent electron liquid correlation 
!  energies for local spin density calculations: 
!  A critical analysis,”  the 1th functional
! Can. J. Phys., 58 (1980) 1200-11.
!                                                                
! This uses the Random Phase Approximation (RPA) parameters     
! rather than those for the Ceperley-Alder solution,        
! and is commonly used in the B3LYP hybrid method.   
!
!



      subroutine functional_vwn1rpa
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
      CALL functional_vwn1rpa_close
     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)
      ELSE  
      CALL functional_vwn1rpa_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_vwn1rpa_close
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
      t1 = RA+RB
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = t2**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t3
      t8 = 1/(0.1029581201158544D2*t4+t6+0.427198D2)
      t11 = dlog(0.6203504908994D0*t3*t8)
      t12 = 0.310907D-1*t11
      t13 = 0.1575246635799487D1*t4
      t17 = datan(0.4489988864D-1/(t13+0.13072D2))
      t18 = 0.2052197294D2*t17
      t19 = 0.7876233178997433D0*t4
      t21 = (t19+0.409286D0)**2
      t23 = dlog(t21*t8)
      t24 = 0.4431373769D-2*t23
      t27 = 1/(0.1584942278842832D2*t4+t6+0.101578D3)
      t30 = dlog(0.6203504908994D0*t3*t27)
      t35 = datan(0.1171685282D1/(t13+0.201231D2))
      t38 = (t19+0.743294D0)**2
      t40 = dlog(t38*t27)
      t45 = (RA-1.D0*RB)*t2
      t46 = 1.D0+t45
      t47 = t46**(1.D0/3.D0)
      t50 = 1.D0-1.D0*t45
      t51 = t50**(1.D0/3.D0)
      t57 = t1*(t12+t18+t24+0.1923661050931536D1*(0.1554535D-1*t30+0.618
     &8180274D0*t35+0.2667310006D-2*t40-t12-t18-t24)*(t47*t46+t51*t50-2.
     &D0))
      F(i) = F(i) +      t57 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = t2**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t3
      t7 = 0.1029581201158544D2*t4+t6+0.427198D2
      t8 = 1/t7
      t11 = dlog(0.6203504908994D0*t3*t8)
      t12 = 0.310907D-1*t11
      t13 = 0.1575246635799487D1*t4
      t14 = t13+0.13072D2
      t17 = datan(0.4489988864D-1/t14)
      t18 = 0.2052197294D2*t17
      t19 = 0.7876233178997433D0*t4
      t20 = t19+0.409286D0
      t21 = t20**2
      t23 = dlog(t21*t8)
      t24 = 0.4431373769D-2*t23
      t26 = 0.1584942278842832D2*t4+t6+0.101578D3
      t27 = 1/t26
      t30 = dlog(0.6203504908994D0*t3*t27)
      t32 = t13+0.201231D2
      t35 = datan(0.1171685282D1/t32)
      t37 = t19+0.743294D0
      t38 = t37**2
      t40 = dlog(t38*t27)
      t42 = 0.1554535D-1*t30+0.6188180274D0*t35+0.2667310006D-2*t40-t12-
     &t18-t24
      t44 = RA-1.D0*RB
      t45 = t44*t2
      t46 = 1.D0+t45
      t47 = t46**(1.D0/3.D0)
      t50 = 1.D0-1.D0*t45
      t51 = t50**(1.D0/3.D0)
      t53 = t47*t46+t51*t50-2.D0
      t56 = t3**2
      t57 = 1/t56
      t59 = t1**2
      t60 = 1/t59
      t63 = t7**2
      t64 = 1/t63
      t66 = t4**2
      t67 = t66**2
      t69 = 1/t67/t4
      t70 = t69*t60
      t73 = 0.2067834969664667D0*t57*t60
      t74 = -0.1715968668597574D1*t70-t73
      t78 = 1/t3
      t80 = (-0.2067834969664667D0*t57*t8*t60-0.6203504908994D0*t3*t64*t
     &74)*t78*t7
      t82 = t14**2
      t83 = 1/t82
      t89 = t83*t69*t60/(1.D0+0.2015999999884401D-2*t83)
      t101 = 0.4431373769D-2*(-0.2625411059665811D0*t20*t8*t70-1.D0*t21*
     &t64*t74)/t21*t7
      t105 = t26**2
      t106 = 1/t105
      t109 = -0.2641570464738054D1*t70-t73
      t116 = t32**2
      t117 = 1/t116
      t141 = t44*t60
      t155 = t12+t18+t24+0.1923661050931536D1*t42*t53+t1*(0.501179582447
     &3985D-1*t80+0.2419143801132913D0*t89+t101+0.1923661050931536D1*(0.
     &2505897912236993D-1*(-0.2067834969664667D0*t57*t27*t60-0.620350490
     &8994D0*t3*t106*t109)*t78*t26+0.190358047713073D0*t117*t69*t60/(1.D
     &0+0.137284640005542D1*t117)+0.2667310006D-2*(-0.2625411059665811D0
     &*t37*t27*t70-1.D0*t38*t106*t109)/t38*t26-0.5011795824473985D-1*t80
     &-0.2419143801132913D0*t89-t101)*t53+0.1923661050931536D1*t42*(0.13
     &33333333333333D1*t47*(t2-1.D0*t141)+0.1333333333333333D1*t51*(-1.D
     &0*t2+t141)))
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t155 
      END IF
      IF (IDERIV .GE. 1) THEN
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS)= D1F(i, ID_RA_POS)
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
      subroutine functional_vwn1rpa_open
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
      t1 = RA+RB
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = t2**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t3
      t8 = 1/(0.1029581201158544D2*t4+t6+0.427198D2)
      t11 = dlog(0.6203504908994D0*t3*t8)
      t12 = 0.310907D-1*t11
      t13 = 0.1575246635799487D1*t4
      t17 = datan(0.4489988864D-1/(t13+0.13072D2))
      t18 = 0.2052197294D2*t17
      t19 = 0.7876233178997433D0*t4
      t21 = (t19+0.409286D0)**2
      t23 = dlog(t21*t8)
      t24 = 0.4431373769D-2*t23
      t27 = 1/(0.1584942278842832D2*t4+t6+0.101578D3)
      t30 = dlog(0.6203504908994D0*t3*t27)
      t35 = datan(0.1171685282D1/(t13+0.201231D2))
      t38 = (t19+0.743294D0)**2
      t40 = dlog(t38*t27)
      t45 = (RA-1.D0*RB)*t2
      t46 = 1.D0+t45
      t47 = t46**(1.D0/3.D0)
      t50 = 1.D0-1.D0*t45
      t51 = t50**(1.D0/3.D0)
      t57 = t1*(t12+t18+t24+0.1923661050931536D1*(0.1554535D-1*t30+0.618
     &8180274D0*t35+0.2667310006D-2*t40-t12-t18-t24)*(t47*t46+t51*t50-2.
     &D0))
      F(i) = F(i) +      t57 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = t2**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t3
      t7 = 0.1029581201158544D2*t4+t6+0.427198D2
      t8 = 1/t7
      t11 = dlog(0.6203504908994D0*t3*t8)
      t12 = 0.310907D-1*t11
      t13 = 0.1575246635799487D1*t4
      t14 = t13+0.13072D2
      t17 = datan(0.4489988864D-1/t14)
      t18 = 0.2052197294D2*t17
      t19 = 0.7876233178997433D0*t4
      t20 = t19+0.409286D0
      t21 = t20**2
      t23 = dlog(t21*t8)
      t24 = 0.4431373769D-2*t23
      t26 = 0.1584942278842832D2*t4+t6+0.101578D3
      t27 = 1/t26
      t30 = dlog(0.6203504908994D0*t3*t27)
      t32 = t13+0.201231D2
      t35 = datan(0.1171685282D1/t32)
      t37 = t19+0.743294D0
      t38 = t37**2
      t40 = dlog(t38*t27)
      t42 = 0.1554535D-1*t30+0.6188180274D0*t35+0.2667310006D-2*t40-t12-
     &t18-t24
      t44 = RA-1.D0*RB
      t45 = t44*t2
      t46 = 1.D0+t45
      t47 = t46**(1.D0/3.D0)
      t50 = 1.D0-1.D0*t45
      t51 = t50**(1.D0/3.D0)
      t53 = t47*t46+t51*t50-2.D0
      t56 = t3**2
      t57 = 1/t56
      t59 = t1**2
      t60 = 1/t59
      t63 = t7**2
      t64 = 1/t63
      t66 = t4**2
      t67 = t66**2
      t69 = 1/t67/t4
      t70 = t69*t60
      t73 = 0.2067834969664667D0*t57*t60
      t74 = -0.1715968668597574D1*t70-t73
      t78 = 1/t3
      t80 = (-0.2067834969664667D0*t57*t8*t60-0.6203504908994D0*t3*t64*t
     &74)*t78*t7
      t82 = t14**2
      t83 = 1/t82
      t89 = t83*t69*t60/(1.D0+0.2015999999884401D-2*t83)
      t101 = 0.4431373769D-2*(-0.2625411059665811D0*t20*t8*t70-1.D0*t21*
     &t64*t74)/t21*t7
      t105 = t26**2
      t106 = 1/t105
      t109 = -0.2641570464738054D1*t70-t73
      t116 = t32**2
      t117 = 1/t116
      t141 = t44*t60
      t155 = t12+t18+t24+0.1923661050931536D1*t42*t53+t1*(0.501179582447
     &3985D-1*t80+0.2419143801132913D0*t89+t101+0.1923661050931536D1*(0.
     &2505897912236993D-1*(-0.2067834969664667D0*t57*t27*t60-0.620350490
     &8994D0*t3*t106*t109)*t78*t26+0.190358047713073D0*t117*t69*t60/(1.D
     &0+0.137284640005542D1*t117)+0.2667310006D-2*(-0.2625411059665811D0
     &*t37*t27*t70-1.D0*t38*t106*t109)/t38*t26-0.5011795824473985D-1*t80
     &-0.2419143801132913D0*t89-t101)*t53+0.1923661050931536D1*t42*(0.13
     &33333333333333D1*t47*(t2-1.D0*t141)+0.1333333333333333D1*t51*(-1.D
     &0*t2+t141)))
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t155 
      t1 = RA+RB
      t2 = 1/t1
      t3 = t2**(1.D0/3.D0)
      t4 = t2**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t3
      t7 = 0.1029581201158544D2*t4+t6+0.427198D2
      t8 = 1/t7
      t11 = dlog(0.6203504908994D0*t3*t8)
      t12 = 0.310907D-1*t11
      t13 = 0.1575246635799487D1*t4
      t14 = t13+0.13072D2
      t17 = datan(0.4489988864D-1/t14)
      t18 = 0.2052197294D2*t17
      t19 = 0.7876233178997433D0*t4
      t20 = t19+0.409286D0
      t21 = t20**2
      t23 = dlog(t21*t8)
      t24 = 0.4431373769D-2*t23
      t26 = 0.1584942278842832D2*t4+t6+0.101578D3
      t27 = 1/t26
      t30 = dlog(0.6203504908994D0*t3*t27)
      t32 = t13+0.201231D2
      t35 = datan(0.1171685282D1/t32)
      t37 = t19+0.743294D0
      t38 = t37**2
      t40 = dlog(t38*t27)
      t42 = 0.1554535D-1*t30+0.6188180274D0*t35+0.2667310006D-2*t40-t12-
     &t18-t24
      t44 = RA-1.D0*RB
      t45 = t44*t2
      t46 = 1.D0+t45
      t47 = t46**(1.D0/3.D0)
      t50 = 1.D0-1.D0*t45
      t51 = t50**(1.D0/3.D0)
      t53 = t47*t46+t51*t50-2.D0
      t56 = t3**2
      t57 = 1/t56
      t59 = t1**2
      t60 = 1/t59
      t63 = t7**2
      t64 = 1/t63
      t66 = t4**2
      t67 = t66**2
      t69 = 1/t67/t4
      t70 = t69*t60
      t73 = 0.2067834969664667D0*t57*t60
      t74 = -0.1715968668597574D1*t70-t73
      t78 = 1/t3
      t80 = (-0.2067834969664667D0*t57*t8*t60-0.6203504908994D0*t3*t64*t
     &74)*t78*t7
      t82 = t14**2
      t83 = 1/t82
      t89 = t83*t69*t60/(1.D0+0.2015999999884401D-2*t83)
      t101 = 0.4431373769D-2*(-0.2625411059665811D0*t20*t8*t70-1.D0*t21*
     &t64*t74)/t21*t7
      t105 = t26**2
      t106 = 1/t105
      t109 = -0.2641570464738054D1*t70-t73
      t116 = t32**2
      t117 = 1/t116
      t142 = t44*t60
      t155 = t12+t18+t24+0.1923661050931536D1*t42*t53+t1*(0.501179582447
     &3985D-1*t80+0.2419143801132913D0*t89+t101+0.1923661050931536D1*(0.
     &2505897912236993D-1*(-0.2067834969664667D0*t57*t27*t60-0.620350490
     &8994D0*t3*t106*t109)*t78*t26+0.190358047713073D0*t117*t69*t60/(1.D
     &0+0.137284640005542D1*t117)+0.2667310006D-2*(-0.2625411059665811D0
     &*t37*t27*t70-1.D0*t38*t106*t109)/t38*t26-0.5011795824473985D-1*t80
     &-0.2419143801132913D0*t89-t101)*t53+0.1923661050931536D1*t42*(0.13
     &33333333333333D1*t47*(-1.D0*t2-1.D0*t142)+0.1333333333333333D1*t51
     &*(t2+t142)))
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t155 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = 1/RA
      t2 = t1**(1.D0/3.D0)
      t3 = t1**(1.D0/6.D0)
      t7 = 1/(0.1584942278842832D2*t3+0.6203504908994D0*t2+0.101578D3)
      t10 = dlog(0.6203504908994D0*t2*t7)
      t16 = datan(0.1171685282D1/(0.1575246635799487D1*t3+0.201231D2))
      t20 = (0.7876233178997433D0*t3+0.743294D0)**2
      t22 = dlog(t20*t7)
      t25 = RA*(0.1554535D-1*t10+0.6188180274D0*t16+0.2667310006D-2*t22)
      F(i) = F(i) +      t25 
      IF (IDERIV .GE. 1) THEN
      t1 = 1/RA
      t2 = t1**(1.D0/3.D0)
      t3 = t1**(1.D0/6.D0)
      t6 = 0.1584942278842832D2*t3+0.6203504908994D0*t2+0.101578D3
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t13 = 0.1575246635799487D1*t3+0.201231D2
      t16 = datan(0.1171685282D1/t13)
      t19 = 0.7876233178997433D0*t3+0.743294D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t24 = t2**2
      t25 = 1/t24
      t27 = RA**2
      t28 = 1/t27
      t31 = t6**2
      t32 = 1/t31
      t34 = t3**2
      t35 = t34**2
      t37 = 1/t35/t3
      t38 = t37*t28
      t42 = -0.2641570464738054D1*t38-0.2067834969664667D0*t25*t28
      t50 = t13**2
      t51 = 1/t50
      t72 = 0.1554535D-1*t10+0.6188180274D0*t16+0.2667310006D-2*t22+RA*(
     &0.2505897912236993D-1*(-0.2067834969664667D0*t25*t7*t28-0.62035049
     &08994D0*t2*t32*t42)/t2*t6+0.190358047713073D0*t51*t37*t28/(1.D0+0.
     &137284640005542D1*t51)+0.2667310006D-2*(-0.2625411059665811D0*t19*
     &t7*t38-1.D0*t20*t32*t42)/t20*t6)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t72 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = 1/RB
      t2 = t1**(1.D0/3.D0)
      t3 = t1**(1.D0/6.D0)
      t7 = 1/(0.1584942278842832D2*t3+0.6203504908994D0*t2+0.101578D3)
      t10 = dlog(0.6203504908994D0*t2*t7)
      t16 = datan(0.1171685282D1/(0.1575246635799487D1*t3+0.201231D2))
      t20 = (0.7876233178997433D0*t3+0.743294D0)**2
      t22 = dlog(t20*t7)
      t25 = RB*(0.1554535D-1*t10+0.6188180274D0*t16+0.2667310006D-2*t22)
      F(i) = F(i) +      t25 
      IF (IDERIV .GE. 1) THEN
      t1 = 1/RB
      t2 = t1**(1.D0/3.D0)
      t3 = t1**(1.D0/6.D0)
      t6 = 0.1584942278842832D2*t3+0.6203504908994D0*t2+0.101578D3
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t13 = 0.1575246635799487D1*t3+0.201231D2
      t16 = datan(0.1171685282D1/t13)
      t19 = 0.7876233178997433D0*t3+0.743294D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t24 = t2**2
      t25 = 1/t24
      t27 = RB**2
      t28 = 1/t27
      t31 = t6**2
      t32 = 1/t31
      t34 = t3**2
      t35 = t34**2
      t37 = 1/t35/t3
      t38 = t37*t28
      t42 = -0.2641570464738054D1*t38-0.2067834969664667D0*t25*t28
      t50 = t13**2
      t51 = 1/t50
      t72 = 0.1554535D-1*t10+0.6188180274D0*t16+0.2667310006D-2*t22+RB*(
     &0.2505897912236993D-1*(-0.2067834969664667D0*t25*t7*t28-0.62035049
     &08994D0*t2*t32*t42)/t2*t6+0.190358047713073D0*t51*t37*t28/(1.D0+0.
     &137284640005542D1*t51)+0.2667310006D-2*(-0.2625411059665811D0*t19*
     &t7*t38-1.D0*t20*t32*t42)/t20*t6)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t72 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
