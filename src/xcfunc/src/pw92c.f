!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
!  J.P. Perdew and Y. Wang, 
!  “Accurate and simple analytic representation of the
!   electron-gas correlation energy”, 
!  Phys. Rev. B, 45(23):13244–13249, 1992
!
!



      subroutine functional_pw92c
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
      CALL functional_pw92c_close
     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)
      ELSE  
      CALL functional_pw92c_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_pw92c_close
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
      t2 = PI*t1
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t7 = t2**(1.D0/6.D0)
      t8 = 1/t7
      t11 = dsqrt(t2)
      t12 = 1/t11
      t14 = t3**2
      t15 = 1/t14
      t21 = dlog(1.D0+0.160819795D2/(0.7240101934316831D1*t8+0.325955091
     &9422292D1*t4+0.1418722816479667D1*t12+0.4069130045175293D0*t15))
      t23 = 0.621814D-1*(0.1941593353441141D0*t4+1.D0)*t21
      t34 = dlog(1.D0+0.2960874998D2/(0.9872129722569272D1*t8+0.32918048
     &09945063D1*t4+0.76232752193529D0*t12+0.4100250709496125D0*t15))
      t37 = RA-1.D0*RB
      t39 = t37/t1
      t40 = 1.D0+t39
      t41 = t40**(1.D0/3.D0)
      t44 = 1.D0-1.D0*t39
      t45 = t44**(1.D0/3.D0)
      t47 = t41*t40+t45*t44-2.D0
      t48 = t37**2
      t49 = t48**2
      t50 = t1**2
      t51 = t50**2
      t53 = t49/t51
      t69 = dlog(1.D0+0.32163959D2/(0.1345791371439445D2*t8+0.5630984149
     &097876D1*t4+0.2915214714219177D1*t12+0.5160664645478634D0*t15))
      t77 = t1*(-t23+0.37995525D-1*(0.1010773329762878D0*t4+1.D0)*t34*t4
     &7*(1.D0-1.D0*t53)+0.1923661050931536D1*(-0.310907D-1*(0.1866909697
     &07574D0*t4+1.D0)*t69+t23)*t47*t53)
      F(i) = F(i) +      t77 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = PI*t1
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t6 = 0.1941593353441141D0*t4+1.D0
      t7 = t2**(1.D0/6.D0)
      t8 = 1/t7
      t11 = dsqrt(t2)
      t12 = 1/t11
      t14 = t3**2
      t15 = 1/t14
      t17 = 0.7240101934316831D1*t8+0.3259550919422292D1*t4+0.1418722816
     &479667D1*t12+0.4069130045175293D0*t15
      t20 = 1.D0+0.160819795D2/t17
      t21 = dlog(t20)
      t23 = 0.621814D-1*t6*t21
      t25 = 0.1010773329762878D0*t4+1.D0
      t30 = 0.9872129722569272D1*t8+0.3291804809945063D1*t4+0.7623275219
     &3529D0*t12+0.4100250709496125D0*t15
      t33 = 1.D0+0.2960874998D2/t30
      t34 = dlog(t33)
      t35 = t25*t34
      t37 = RA-1.D0*RB
      t38 = 1/t1
      t39 = t37*t38
      t40 = 1.D0+t39
      t41 = t40**(1.D0/3.D0)
      t44 = 1.D0-1.D0*t39
      t45 = t44**(1.D0/3.D0)
      t47 = t41*t40+t45*t44-2.D0
      t48 = t37**2
      t49 = t48**2
      t50 = t1**2
      t51 = t50**2
      t52 = 1/t51
      t53 = t49*t52
      t55 = 1.D0-1.D0*t53
      t60 = 0.186690969707574D0*t4+1.D0
      t65 = 0.1345791371439445D2*t8+0.5630984149097876D1*t4+0.2915214714
     &219177D1*t12+0.5160664645478634D0*t15
      t68 = 1.D0+0.32163959D2/t65
      t69 = dlog(t68)
      t72 = -0.310907D-1*t60*t69+t23
      t73 = t72*t47
      t78 = 1/t3/t2*PI
      t79 = t78*t21
      t81 = t17**2
      t86 = 1/t7/t2*PI
      t89 = t11**2
      t92 = 1/t89/t11*PI
      t96 = 1/t14/t2*PI
      t102 = 0.10000000000813D1*t6/t81*(-0.1206683655719472D1*t86-0.1086
     &516973140764D1*t78-0.7093614082398337D0*t92-0.2712753363450195D0*t
     &96)/t20
      t107 = t30**2
      t122 = t37/t50
      t131 = 0.1333333333333333D1*t41*(t38-1.D0*t122)+0.1333333333333333
     &D1*t45*(-1.D0*t38+t122)
      t136 = t48*t37*t52
      t140 = t49/t51/t1
      t148 = t65**2
      t174 = -t23+0.37995525D-1*t35*t47*t55+0.1923661050931536D1*t73*t53
     &+t1*(0.4024366431588833D-2*t79+t102-0.1280162110677955D-2*t78*t34*
     &t47*t55-0.1125000000083839D1*t25/t107*(-0.1645354953761545D1*t86-0
     &.1097268269981688D1*t78-0.381163760967645D0*t92-0.2733500472997417
     &D0*t96)/t33*t47*t55+0.37995525D-1*t35*t131*t55+0.37995525D-1*t35*t
     &47*(-4.D0*t136+4.D0*t140)+0.1923661050931536D1*(0.1934784310629091
     &D-2*t78*t69+0.10000000000813D1*t60/t148*(-0.2242985619065741D1*t86
     &-0.1876994716365959D1*t78-0.1457607357109589D1*t92-0.3440443096985
     &756D0*t96)/t68-0.4024366431588833D-2*t79-t102)*t47*t53+0.192366105
     &0931536D1*t72*t131*t53+0.7694644203726145D1*t73*t136-0.76946442037
     &26145D1*t73*t140)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t174 
      END IF
      IF (IDERIV .GE. 1) THEN
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS)= D1F(i, ID_RA_POS)
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
      subroutine functional_pw92c_open
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
      t2 = PI*t1
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t7 = t2**(1.D0/6.D0)
      t8 = 1/t7
      t11 = dsqrt(t2)
      t12 = 1/t11
      t14 = t3**2
      t15 = 1/t14
      t21 = dlog(1.D0+0.160819795D2/(0.7240101934316831D1*t8+0.325955091
     &9422292D1*t4+0.1418722816479667D1*t12+0.4069130045175293D0*t15))
      t23 = 0.621814D-1*(0.1941593353441141D0*t4+1.D0)*t21
      t34 = dlog(1.D0+0.2960874998D2/(0.9872129722569272D1*t8+0.32918048
     &09945063D1*t4+0.76232752193529D0*t12+0.4100250709496125D0*t15))
      t37 = RA-1.D0*RB
      t39 = t37/t1
      t40 = 1.D0+t39
      t41 = t40**(1.D0/3.D0)
      t44 = 1.D0-1.D0*t39
      t45 = t44**(1.D0/3.D0)
      t47 = t41*t40+t45*t44-2.D0
      t48 = t37**2
      t49 = t48**2
      t50 = t1**2
      t51 = t50**2
      t53 = t49/t51
      t69 = dlog(1.D0+0.32163959D2/(0.1345791371439445D2*t8+0.5630984149
     &097876D1*t4+0.2915214714219177D1*t12+0.5160664645478634D0*t15))
      t77 = t1*(-t23+0.37995525D-1*(0.1010773329762878D0*t4+1.D0)*t34*t4
     &7*(1.D0-1.D0*t53)+0.1923661050931536D1*(-0.310907D-1*(0.1866909697
     &07574D0*t4+1.D0)*t69+t23)*t47*t53)
      F(i) = F(i) +      t77 
      IF (IDERIV .GE. 1) THEN
      t1 = RA+RB
      t2 = PI*t1
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t6 = 0.1941593353441141D0*t4+1.D0
      t7 = t2**(1.D0/6.D0)
      t8 = 1/t7
      t11 = dsqrt(t2)
      t12 = 1/t11
      t14 = t3**2
      t15 = 1/t14
      t17 = 0.7240101934316831D1*t8+0.3259550919422292D1*t4+0.1418722816
     &479667D1*t12+0.4069130045175293D0*t15
      t20 = 1.D0+0.160819795D2/t17
      t21 = dlog(t20)
      t23 = 0.621814D-1*t6*t21
      t25 = 0.1010773329762878D0*t4+1.D0
      t30 = 0.9872129722569272D1*t8+0.3291804809945063D1*t4+0.7623275219
     &3529D0*t12+0.4100250709496125D0*t15
      t33 = 1.D0+0.2960874998D2/t30
      t34 = dlog(t33)
      t35 = t25*t34
      t37 = RA-1.D0*RB
      t38 = 1/t1
      t39 = t37*t38
      t40 = 1.D0+t39
      t41 = t40**(1.D0/3.D0)
      t44 = 1.D0-1.D0*t39
      t45 = t44**(1.D0/3.D0)
      t47 = t41*t40+t45*t44-2.D0
      t48 = t37**2
      t49 = t48**2
      t50 = t1**2
      t51 = t50**2
      t52 = 1/t51
      t53 = t49*t52
      t55 = 1.D0-1.D0*t53
      t60 = 0.186690969707574D0*t4+1.D0
      t65 = 0.1345791371439445D2*t8+0.5630984149097876D1*t4+0.2915214714
     &219177D1*t12+0.5160664645478634D0*t15
      t68 = 1.D0+0.32163959D2/t65
      t69 = dlog(t68)
      t72 = -0.310907D-1*t60*t69+t23
      t73 = t72*t47
      t78 = 1/t3/t2*PI
      t79 = t78*t21
      t81 = t17**2
      t86 = 1/t7/t2*PI
      t89 = t11**2
      t92 = 1/t89/t11*PI
      t96 = 1/t14/t2*PI
      t102 = 0.10000000000813D1*t6/t81*(-0.1206683655719472D1*t86-0.1086
     &516973140764D1*t78-0.7093614082398337D0*t92-0.2712753363450195D0*t
     &96)/t20
      t107 = t30**2
      t122 = t37/t50
      t131 = 0.1333333333333333D1*t41*(t38-1.D0*t122)+0.1333333333333333
     &D1*t45*(-1.D0*t38+t122)
      t136 = t48*t37*t52
      t140 = t49/t51/t1
      t148 = t65**2
      t174 = -t23+0.37995525D-1*t35*t47*t55+0.1923661050931536D1*t73*t53
     &+t1*(0.4024366431588833D-2*t79+t102-0.1280162110677955D-2*t78*t34*
     &t47*t55-0.1125000000083839D1*t25/t107*(-0.1645354953761545D1*t86-0
     &.1097268269981688D1*t78-0.381163760967645D0*t92-0.2733500472997417
     &D0*t96)/t33*t47*t55+0.37995525D-1*t35*t131*t55+0.37995525D-1*t35*t
     &47*(-4.D0*t136+4.D0*t140)+0.1923661050931536D1*(0.1934784310629091
     &D-2*t78*t69+0.10000000000813D1*t60/t148*(-0.2242985619065741D1*t86
     &-0.1876994716365959D1*t78-0.1457607357109589D1*t92-0.3440443096985
     &756D0*t96)/t68-0.4024366431588833D-2*t79-t102)*t47*t53+0.192366105
     &0931536D1*t72*t131*t53+0.7694644203726145D1*t73*t136-0.76946442037
     &26145D1*t73*t140)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t174 
      t1 = RA+RB
      t2 = PI*t1
      t3 = t2**(1.D0/3.D0)
      t4 = 1/t3
      t6 = 0.1941593353441141D0*t4+1.D0
      t7 = t2**(1.D0/6.D0)
      t8 = 1/t7
      t11 = dsqrt(t2)
      t12 = 1/t11
      t14 = t3**2
      t15 = 1/t14
      t17 = 0.7240101934316831D1*t8+0.3259550919422292D1*t4+0.1418722816
     &479667D1*t12+0.4069130045175293D0*t15
      t20 = 1.D0+0.160819795D2/t17
      t21 = dlog(t20)
      t23 = 0.621814D-1*t6*t21
      t25 = 0.1010773329762878D0*t4+1.D0
      t30 = 0.9872129722569272D1*t8+0.3291804809945063D1*t4+0.7623275219
     &3529D0*t12+0.4100250709496125D0*t15
      t33 = 1.D0+0.2960874998D2/t30
      t34 = dlog(t33)
      t35 = t25*t34
      t37 = RA-1.D0*RB
      t38 = 1/t1
      t39 = t37*t38
      t40 = 1.D0+t39
      t41 = t40**(1.D0/3.D0)
      t44 = 1.D0-1.D0*t39
      t45 = t44**(1.D0/3.D0)
      t47 = t41*t40+t45*t44-2.D0
      t48 = t37**2
      t49 = t48**2
      t50 = t1**2
      t51 = t50**2
      t52 = 1/t51
      t53 = t49*t52
      t55 = 1.D0-1.D0*t53
      t60 = 0.186690969707574D0*t4+1.D0
      t65 = 0.1345791371439445D2*t8+0.5630984149097876D1*t4+0.2915214714
     &219177D1*t12+0.5160664645478634D0*t15
      t68 = 1.D0+0.32163959D2/t65
      t69 = dlog(t68)
      t72 = -0.310907D-1*t60*t69+t23
      t73 = t72*t47
      t78 = 1/t3/t2*PI
      t79 = t78*t21
      t81 = t17**2
      t86 = 1/t7/t2*PI
      t89 = t11**2
      t92 = 1/t89/t11*PI
      t96 = 1/t14/t2*PI
      t102 = 0.10000000000813D1*t6/t81*(-0.1206683655719472D1*t86-0.1086
     &516973140764D1*t78-0.7093614082398337D0*t92-0.2712753363450195D0*t
     &96)/t20
      t107 = t30**2
      t123 = t37/t50
      t131 = 0.1333333333333333D1*t41*(-1.D0*t38-1.D0*t123)+0.1333333333
     &333333D1*t45*(t38+t123)
      t136 = t48*t37*t52
      t140 = t49/t51/t1
      t148 = t65**2
      t174 = -t23+0.37995525D-1*t35*t47*t55+0.1923661050931536D1*t73*t53
     &+t1*(0.4024366431588833D-2*t79+t102-0.1280162110677955D-2*t78*t34*
     &t47*t55-0.1125000000083839D1*t25/t107*(-0.1645354953761545D1*t86-0
     &.1097268269981688D1*t78-0.381163760967645D0*t92-0.2733500472997417
     &D0*t96)/t33*t47*t55+0.37995525D-1*t35*t131*t55+0.37995525D-1*t35*t
     &47*(4.D0*t136+4.D0*t140)+0.1923661050931536D1*(0.1934784310629091D
     &-2*t78*t69+0.10000000000813D1*t60/t148*(-0.2242985619065741D1*t86-
     &0.1876994716365959D1*t78-0.1457607357109589D1*t92-0.34404430969857
     &56D0*t96)/t68-0.4024366431588833D-2*t79-t102)*t47*t53+0.1923661050
     &931536D1*t72*t131*t53-0.7694644203726145D1*t73*t136-0.769464420372
     &6145D1*t73*t140)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t174 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = PI*RA
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t7 = t1**(1.D0/6.D0)
      t11 = dsqrt(t1)
      t14 = t2**2
      t21 = dlog(1.D0+0.32163959D2/(0.1345791371439445D2/t7+0.5630984149
     &097876D1*t3+0.2915214714219177D1/t11+0.5160664645478634D0/t14))
      t24 = -0.310907D-1*RA*(0.186690969707574D0*t3+1.D0)*t21
      F(i) = F(i) +      t24 
      IF (IDERIV .GE. 1) THEN
      t1 = PI*RA
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = 0.186690969707574D0*t3+1.D0
      t6 = t1**(1.D0/6.D0)
      t10 = dsqrt(t1)
      t13 = t2**2
      t16 = 0.1345791371439445D2/t6+0.5630984149097876D1*t3+0.2915214714
     &219177D1/t10+0.5160664645478634D0/t13
      t19 = 1.D0+0.32163959D2/t16
      t20 = dlog(t19)
      t24 = 1/t2/t1
      t30 = t16**2
      t38 = t10**2
      t53 = -0.310907D-1*t5*t20+0.1934784310629091D-2*RA*t24*PI*t20+0.10
     &000000000813D1*RA*t5/t30*(-0.2242985619065741D1/t6/t1*PI-0.1876994
     &716365959D1*t24*PI-0.1457607357109589D1/t38/t10*PI-0.3440443096985
     &756D0/t13/t1*PI)/t19
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t53 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = PI*RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t7 = t1**(1.D0/6.D0)
      t11 = dsqrt(t1)
      t14 = t2**2
      t21 = dlog(1.D0+0.32163959D2/(0.1345791371439445D2/t7+0.5630984149
     &097876D1*t3+0.2915214714219177D1/t11+0.5160664645478634D0/t14))
      t24 = -0.310907D-1*RB*(0.186690969707574D0*t3+1.D0)*t21
      F(i) = F(i) +      t24 
      IF (IDERIV .GE. 1) THEN
      t1 = PI*RB
      t2 = t1**(1.D0/3.D0)
      t3 = 1/t2
      t5 = 0.186690969707574D0*t3+1.D0
      t6 = t1**(1.D0/6.D0)
      t10 = dsqrt(t1)
      t13 = t2**2
      t16 = 0.1345791371439445D2/t6+0.5630984149097876D1*t3+0.2915214714
     &219177D1/t10+0.5160664645478634D0/t13
      t19 = 1.D0+0.32163959D2/t16
      t20 = dlog(t19)
      t24 = 1/t2/t1
      t30 = t16**2
      t38 = t10**2
      t53 = -0.310907D-1*t5*t20+0.1934784310629091D-2*RB*t24*PI*t20+0.10
     &000000000813D1*RB*t5/t30*(-0.2242985619065741D1/t6/t1*PI-0.1876994
     &716365959D1*t24*PI-0.1457607357109589D1/t38/t10*PI-0.3440443096985
     &756D0/t13/t1*PI)/t19
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t53 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
