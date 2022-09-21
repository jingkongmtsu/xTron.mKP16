!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, M.R. Pederson, 
! D.J. Singh and C. Fiolhais, 
! “Atoms, molecules, solids and surfaces: Applications of the generalized
!  gradient approximation for exchange and correlation”, 
! Phys. Rev. B, 46(11):6671–6687, 1992
!
!



      subroutine functional_pw91x
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
      CALL functional_pw91x_close
     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)
      ELSE  
      CALL functional_pw91x_open
     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)
      END IF  
      RETURN  
      END     



      subroutine functional_pw91x_close
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
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t5 = 1/PI
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t15 = RA**2
      t16 = 1/t15
      t17 = GAA*t13*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t25 = 0.5405530067264707D-1*t8*t9*t22
      t27 = dexp(-0.7571335803467249D1*t17)
      t36 = GAA**2
      t40 = t15**2
      t49 = t1*RB
      t50 = t49**(1.D0/3.D0)
      t52 = dsqrt(GBB)
      t54 = t52/t50
      t55 = 1/RB
      t58 = t50**2
      t59 = 1/t58
      t61 = RB**2
      t62 = 1/t61
      t63 = GBB*t59*t62
      t66 = dsqrt(0.4601205204688956D1*t63+1)
      t68 = dlog(0.2145042005343708D1*t54*t55+t66)
      t71 = 0.5405530067264707D-1*t54*t55*t68
      t73 = dexp(-0.7571335803467249D1*t63)
      t82 = GBB**2
      t86 = t61**2
      t95 = -0.1362840444624105D1*RA*t3*t5*(1.D0+t25+0.7571335803467249D
     &-1*(0.2743D0-0.1508D0*t27)*GAA*t13*t16)/(1.D0+t25+0.22930050338078
     &5D-4*t36/t3/t2/t40)-0.1362840444624105D1*RB*t50*t5*(1.D0+t71+0.757
     &1335803467249D-1*(0.2743D0-0.1508D0*t73)*GBB*t59*t62)/(1.D0+t71+0.
     &229300503380785D-4*t82/t50/t49/t86)
      F(i) = F(i) +      t95 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t14 = GAA*t13
      t15 = RA**2
      t16 = 1/t15
      t17 = t14*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t23 = t9*t22
      t25 = 0.5405530067264707D-1*t8*t23
      t27 = dexp(-0.7571335803467249D1*t17)
      t30 = (0.2743D0-0.1508D0*t27)*GAA
      t34 = 1.D0+t25+0.7571335803467249D-1*t30*t13*t16
      t35 = GAA**2
      t37 = 1/t3/t2
      t38 = t35*t37
      t39 = t15**2
      t40 = 1/t39
      t43 = 1.D0+t25+0.229300503380785D-4*t38*t40
      t44 = 1/t43
      t53 = RA*t3
      t54 = t6*t37
      t57 = 0.1801843355754902D-1*t54*t23*t1
      t60 = 0.5405530067264707D-1*t8*t16*t22
      t70 = dsqrt(1.D0+0.4601205204688955D1*t17)
      t74 = 0.5405530067264707D-1*t8*t9*(-0.7150140017812359D0*t54*t9*t1
     &-0.2145042005343708D1*t8*t16)/t70
      t76 = 1/t12/t2
      t82 = 1/t15/RA
      t102 = t43**2
      t105 = t1**2
      t121 = -0.1362840444624105D1*t3*t4*t34*t44-0.4542801482080349D0*RA
     &*t13*PI*t34*t44-0.1362840444624105D1*t53*t4*(-t57-t60+t74-0.114175
     &7439162861D-1*(0.5047557202311499D1*GAA*t76*t16*t1+0.1514267160693
     &45D2*t14*t82)*t27*t17-0.5047557202311499D-1*t30*t76*t16*t1-0.15142
     &6716069345D0*t30*t13*t82)*t44+0.1362840444624105D1*t53*t4*t34/t102
     &*(-t57-t60+t74-0.3057340045077133D-4*t35/t3/t105/t15*t40*t1-0.9172
     &020135231398D-4*t38/t39/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t121 
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = RA*t3
      t5 = 1/PI
      t6 = dsqrt(GAA)
      t8 = 1/t3
      t10 = 1/RA
      t11 = t6*t8
      t14 = t3**2
      t15 = 1/t14
      t17 = RA**2
      t18 = 1/t17
      t19 = GAA*t15*t18
      t22 = dsqrt(0.4601205204688956D1*t19+1)
      t24 = dlog(0.2145042005343708D1*t11*t10+t22)
      t25 = t10*t24
      t27 = 0.2702765033632353D-1/t6*t8*t25
      t28 = t15*t18
      t31 = dsqrt(1.D0+0.4601205204688955D1*t19)
      t34 = 0.5797544527715597D-1*t28/t31
      t36 = 1/t3/t2
      t37 = t17**2
      t38 = 1/t37
      t41 = dexp(-0.7571335803467249D1*t19)
      t46 = 0.2743D0-0.1508D0*t41
      t53 = 0.5405530067264707D-1*t11*t25
      t54 = GAA**2
      t58 = 1.D0+t53+0.229300503380785D-4*t54*t36*t38
      t68 = t58**2
      t78 = -0.1362840444624105D1*t4*t5*(t27+t34+0.8644628978008849D-1*t
     &36*t38*t41*GAA+0.7571335803467249D-1*t46*t15*t18)/t58+0.1362840444
     &624105D1*t4*t5*(1.D0+t53+0.7571335803467249D-1*t46*GAA*t28)/t68*(t
     &27+t34+0.4586010067615699D-4*GAA*t36*t38)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t78 
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
                   
                   
                   
                   
      subroutine functional_pw91x_open
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
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t5 = 1/PI
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t15 = RA**2
      t16 = 1/t15
      t17 = GAA*t13*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t25 = 0.5405530067264707D-1*t8*t9*t22
      t27 = dexp(-0.7571335803467249D1*t17)
      t36 = GAA**2
      t40 = t15**2
      t49 = t1*RB
      t50 = t49**(1.D0/3.D0)
      t52 = dsqrt(GBB)
      t54 = t52/t50
      t55 = 1/RB
      t58 = t50**2
      t59 = 1/t58
      t61 = RB**2
      t62 = 1/t61
      t63 = GBB*t59*t62
      t66 = dsqrt(0.4601205204688956D1*t63+1)
      t68 = dlog(0.2145042005343708D1*t54*t55+t66)
      t71 = 0.5405530067264707D-1*t54*t55*t68
      t73 = dexp(-0.7571335803467249D1*t63)
      t82 = GBB**2
      t86 = t61**2
      t95 = -0.1362840444624105D1*RA*t3*t5*(1.D0+t25+0.7571335803467249D
     &-1*(0.2743D0-0.1508D0*t27)*GAA*t13*t16)/(1.D0+t25+0.22930050338078
     &5D-4*t36/t3/t2/t40)-0.1362840444624105D1*RB*t50*t5*(1.D0+t71+0.757
     &1335803467249D-1*(0.2743D0-0.1508D0*t73)*GBB*t59*t62)/(1.D0+t71+0.
     &229300503380785D-4*t82/t50/t49/t86)
      F(i) = F(i) +      t95 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t14 = GAA*t13
      t15 = RA**2
      t16 = 1/t15
      t17 = t14*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t23 = t9*t22
      t25 = 0.5405530067264707D-1*t8*t23
      t27 = dexp(-0.7571335803467249D1*t17)
      t30 = (0.2743D0-0.1508D0*t27)*GAA
      t34 = 1.D0+t25+0.7571335803467249D-1*t30*t13*t16
      t35 = GAA**2
      t37 = 1/t3/t2
      t38 = t35*t37
      t39 = t15**2
      t40 = 1/t39
      t43 = 1.D0+t25+0.229300503380785D-4*t38*t40
      t44 = 1/t43
      t53 = RA*t3
      t54 = t6*t37
      t57 = 0.1801843355754902D-1*t54*t23*t1
      t60 = 0.5405530067264707D-1*t8*t16*t22
      t70 = dsqrt(1.D0+0.4601205204688955D1*t17)
      t74 = 0.5405530067264707D-1*t8*t9*(-0.7150140017812359D0*t54*t9*t1
     &-0.2145042005343708D1*t8*t16)/t70
      t76 = 1/t12/t2
      t82 = 1/t15/RA
      t102 = t43**2
      t105 = t1**2
      t121 = -0.1362840444624105D1*t3*t4*t34*t44-0.4542801482080349D0*RA
     &*t13*PI*t34*t44-0.1362840444624105D1*t53*t4*(-t57-t60+t74-0.114175
     &7439162861D-1*(0.5047557202311499D1*GAA*t76*t16*t1+0.1514267160693
     &45D2*t14*t82)*t27*t17-0.5047557202311499D-1*t30*t76*t16*t1-0.15142
     &6716069345D0*t30*t13*t82)*t44+0.1362840444624105D1*t53*t4*t34/t102
     &*(-t57-t60+t74-0.3057340045077133D-4*t35/t3/t105/t15*t40*t1-0.9172
     &020135231398D-4*t38/t39/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t121 
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = dsqrt(GBB)
      t8 = t6/t3
      t9 = 1/RB
      t12 = t3**2
      t13 = 1/t12
      t14 = GBB*t13
      t15 = RB**2
      t16 = 1/t15
      t17 = t14*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t23 = t9*t22
      t25 = 0.5405530067264707D-1*t8*t23
      t27 = dexp(-0.7571335803467249D1*t17)
      t30 = (0.2743D0-0.1508D0*t27)*GBB
      t34 = 1.D0+t25+0.7571335803467249D-1*t30*t13*t16
      t35 = GBB**2
      t37 = 1/t3/t2
      t38 = t35*t37
      t39 = t15**2
      t40 = 1/t39
      t43 = 1.D0+t25+0.229300503380785D-4*t38*t40
      t44 = 1/t43
      t53 = RB*t3
      t54 = t6*t37
      t57 = 0.1801843355754902D-1*t54*t23*t1
      t60 = 0.5405530067264707D-1*t8*t16*t22
      t70 = dsqrt(1.D0+0.4601205204688955D1*t17)
      t74 = 0.5405530067264707D-1*t8*t9*(-0.7150140017812359D0*t54*t9*t1
     &-0.2145042005343708D1*t8*t16)/t70
      t76 = 1/t12/t2
      t82 = 1/t15/RB
      t102 = t43**2
      t105 = t1**2
      t121 = -0.1362840444624105D1*t3*t4*t34*t44-0.4542801482080349D0*RB
     &*t13*PI*t34*t44-0.1362840444624105D1*t53*t4*(-t57-t60+t74-0.114175
     &7439162861D-1*(0.5047557202311499D1*GBB*t76*t16*t1+0.1514267160693
     &45D2*t14*t82)*t27*t17-0.5047557202311499D-1*t30*t76*t16*t1-0.15142
     &6716069345D0*t30*t13*t82)*t44+0.1362840444624105D1*t53*t4*t34/t102
     &*(-t57-t60+t74-0.3057340045077133D-4*t35/t3/t105/t15*t40*t1-0.9172
     &020135231398D-4*t38/t39/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t121 
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = RA*t3
      t5 = 1/PI
      t6 = dsqrt(GAA)
      t8 = 1/t3
      t10 = 1/RA
      t11 = t6*t8
      t14 = t3**2
      t15 = 1/t14
      t17 = RA**2
      t18 = 1/t17
      t19 = GAA*t15*t18
      t22 = dsqrt(0.4601205204688956D1*t19+1)
      t24 = dlog(0.2145042005343708D1*t11*t10+t22)
      t25 = t10*t24
      t27 = 0.2702765033632353D-1/t6*t8*t25
      t28 = t15*t18
      t31 = dsqrt(1.D0+0.4601205204688955D1*t19)
      t34 = 0.5797544527715597D-1*t28/t31
      t36 = 1/t3/t2
      t37 = t17**2
      t38 = 1/t37
      t41 = dexp(-0.7571335803467249D1*t19)
      t46 = 0.2743D0-0.1508D0*t41
      t53 = 0.5405530067264707D-1*t11*t25
      t54 = GAA**2
      t58 = 1.D0+t53+0.229300503380785D-4*t54*t36*t38
      t68 = t58**2
      t78 = -0.1362840444624105D1*t4*t5*(t27+t34+0.8644628978008849D-1*t
     &36*t38*t41*GAA+0.7571335803467249D-1*t46*t15*t18)/t58+0.1362840444
     &624105D1*t4*t5*(1.D0+t53+0.7571335803467249D-1*t46*GAA*t28)/t68*(t
     &27+t34+0.4586010067615699D-4*GAA*t36*t38)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t78 
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = RB*t3
      t5 = 1/PI
      t6 = dsqrt(GBB)
      t8 = 1/t3
      t10 = 1/RB
      t11 = t6*t8
      t14 = t3**2
      t15 = 1/t14
      t17 = RB**2
      t18 = 1/t17
      t19 = GBB*t15*t18
      t22 = dsqrt(0.4601205204688956D1*t19+1)
      t24 = dlog(0.2145042005343708D1*t11*t10+t22)
      t25 = t10*t24
      t27 = 0.2702765033632353D-1/t6*t8*t25
      t28 = t15*t18
      t31 = dsqrt(1.D0+0.4601205204688955D1*t19)
      t34 = 0.5797544527715597D-1*t28/t31
      t36 = 1/t3/t2
      t37 = t17**2
      t38 = 1/t37
      t41 = dexp(-0.7571335803467249D1*t19)
      t46 = 0.2743D0-0.1508D0*t41
      t53 = 0.5405530067264707D-1*t11*t25
      t54 = GBB**2
      t58 = 1.D0+t53+0.229300503380785D-4*t54*t36*t38
      t68 = t58**2
      t78 = -0.1362840444624105D1*t4*t5*(t27+t34+0.8644628978008849D-1*t
     &36*t38*t41*GBB+0.7571335803467249D-1*t46*t15*t18)/t58+0.1362840444
     &624105D1*t4*t5*(1.D0+t53+0.7571335803467249D-1*t46*GBB*t28)/t68*(t
     &27+t34+0.4586010067615699D-4*GBB*t36*t38)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t78 
      t1 = 0.D0
      ID_GAB_POS=D1VARS(ID_GAB)
      D1F(i, ID_GAB_POS) = D1F(i, ID_GAB_POS) +      t1 
      END IF
      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t15 = RA**2
      t16 = 1/t15
      t17 = GAA*t13*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t25 = 0.5405530067264707D-1*t8*t9*t22
      t27 = dexp(-0.7571335803467249D1*t17)
      t36 = GAA**2
      t40 = t15**2
      t49 = -0.1362840444624105D1*RA*t3/PI*(1.D0+t25+0.7571335803467249D
     &-1*(0.2743D0-0.1508D0*t27)*GAA*t13*t16)/(1.D0+t25+0.22930050338078
     &5D-4*t36/t3/t2/t40)
      F(i) = F(i) +      t49 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = dsqrt(GAA)
      t8 = t6/t3
      t9 = 1/RA
      t12 = t3**2
      t13 = 1/t12
      t14 = GAA*t13
      t15 = RA**2
      t16 = 1/t15
      t17 = t14*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t23 = t9*t22
      t25 = 0.5405530067264707D-1*t8*t23
      t27 = dexp(-0.7571335803467249D1*t17)
      t30 = (0.2743D0-0.1508D0*t27)*GAA
      t34 = 1.D0+t25+0.7571335803467249D-1*t30*t13*t16
      t35 = GAA**2
      t37 = 1/t3/t2
      t38 = t35*t37
      t39 = t15**2
      t40 = 1/t39
      t43 = 1.D0+t25+0.229300503380785D-4*t38*t40
      t44 = 1/t43
      t53 = RA*t3
      t54 = t6*t37
      t57 = 0.1801843355754902D-1*t54*t23*t1
      t60 = 0.5405530067264707D-1*t8*t16*t22
      t70 = dsqrt(1.D0+0.4601205204688955D1*t17)
      t74 = 0.5405530067264707D-1*t8*t9*(-0.7150140017812359D0*t54*t9*t1
     &-0.2145042005343708D1*t8*t16)/t70
      t76 = 1/t12/t2
      t82 = 1/t15/RA
      t102 = t43**2
      t105 = t1**2
      t121 = -0.1362840444624105D1*t3*t4*t34*t44-0.4542801482080349D0*RA
     &*t13*PI*t34*t44-0.1362840444624105D1*t53*t4*(-t57-t60+t74-0.114175
     &7439162861D-1*(0.5047557202311499D1*GAA*t76*t16*t1+0.1514267160693
     &45D2*t14*t82)*t27*t17-0.5047557202311499D-1*t30*t76*t16*t1-0.15142
     &6716069345D0*t30*t13*t82)*t44+0.1362840444624105D1*t53*t4*t34/t102
     &*(-t57-t60+t74-0.3057340045077133D-4*t35/t3/t105/t15*t40*t1-0.9172
     &020135231398D-4*t38/t39/RA)
      ID_RA_POS=D1VARS(ID_RA)
      D1F(i, ID_RA_POS) = D1F(i, ID_RA_POS) +      t121 
      t1 = PI**2
      t2 = t1*RA
      t3 = t2**(1.D0/3.D0)
      t4 = RA*t3
      t5 = 1/PI
      t6 = dsqrt(GAA)
      t8 = 1/t3
      t10 = 1/RA
      t11 = t6*t8
      t14 = t3**2
      t15 = 1/t14
      t17 = RA**2
      t18 = 1/t17
      t19 = GAA*t15*t18
      t22 = dsqrt(0.4601205204688956D1*t19+1)
      t24 = dlog(0.2145042005343708D1*t11*t10+t22)
      t25 = t10*t24
      t27 = 0.2702765033632353D-1/t6*t8*t25
      t28 = t15*t18
      t31 = dsqrt(1.D0+0.4601205204688955D1*t19)
      t34 = 0.5797544527715597D-1*t28/t31
      t36 = 1/t3/t2
      t37 = t17**2
      t38 = 1/t37
      t41 = dexp(-0.7571335803467249D1*t19)
      t46 = 0.2743D0-0.1508D0*t41
      t53 = 0.5405530067264707D-1*t11*t25
      t54 = GAA**2
      t58 = 1.D0+t53+0.229300503380785D-4*t54*t36*t38
      t68 = t58**2
      t78 = -0.1362840444624105D1*t4*t5*(t27+t34+0.8644628978008849D-1*t
     &36*t38*t41*GAA+0.7571335803467249D-1*t46*t15*t18)/t58+0.1362840444
     &624105D1*t4*t5*(1.D0+t53+0.7571335803467249D-1*t46*GAA*t28)/t68*(t
     &27+t34+0.4586010067615699D-4*GAA*t36*t38)
      ID_GAA_POS=D1VARS(ID_GAA)
      D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) +      t78 
      END IF
      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t6 = dsqrt(GBB)
      t8 = t6/t3
      t9 = 1/RB
      t12 = t3**2
      t13 = 1/t12
      t15 = RB**2
      t16 = 1/t15
      t17 = GBB*t13*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t25 = 0.5405530067264707D-1*t8*t9*t22
      t27 = dexp(-0.7571335803467249D1*t17)
      t36 = GBB**2
      t40 = t15**2
      t49 = -0.1362840444624105D1*RB*t3/PI*(1.D0+t25+0.7571335803467249D
     &-1*(0.2743D0-0.1508D0*t27)*GBB*t13*t16)/(1.D0+t25+0.22930050338078
     &5D-4*t36/t3/t2/t40)
      F(i) = F(i) +      t49 
      IF (IDERIV .GE. 1) THEN
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = 1/PI
      t6 = dsqrt(GBB)
      t8 = t6/t3
      t9 = 1/RB
      t12 = t3**2
      t13 = 1/t12
      t14 = GBB*t13
      t15 = RB**2
      t16 = 1/t15
      t17 = t14*t16
      t20 = dsqrt(0.4601205204688956D1*t17+1)
      t22 = dlog(0.2145042005343708D1*t8*t9+t20)
      t23 = t9*t22
      t25 = 0.5405530067264707D-1*t8*t23
      t27 = dexp(-0.7571335803467249D1*t17)
      t30 = (0.2743D0-0.1508D0*t27)*GBB
      t34 = 1.D0+t25+0.7571335803467249D-1*t30*t13*t16
      t35 = GBB**2
      t37 = 1/t3/t2
      t38 = t35*t37
      t39 = t15**2
      t40 = 1/t39
      t43 = 1.D0+t25+0.229300503380785D-4*t38*t40
      t44 = 1/t43
      t53 = RB*t3
      t54 = t6*t37
      t57 = 0.1801843355754902D-1*t54*t23*t1
      t60 = 0.5405530067264707D-1*t8*t16*t22
      t70 = dsqrt(1.D0+0.4601205204688955D1*t17)
      t74 = 0.5405530067264707D-1*t8*t9*(-0.7150140017812359D0*t54*t9*t1
     &-0.2145042005343708D1*t8*t16)/t70
      t76 = 1/t12/t2
      t82 = 1/t15/RB
      t102 = t43**2
      t105 = t1**2
      t121 = -0.1362840444624105D1*t3*t4*t34*t44-0.4542801482080349D0*RB
     &*t13*PI*t34*t44-0.1362840444624105D1*t53*t4*(-t57-t60+t74-0.114175
     &7439162861D-1*(0.5047557202311499D1*GBB*t76*t16*t1+0.1514267160693
     &45D2*t14*t82)*t27*t17-0.5047557202311499D-1*t30*t76*t16*t1-0.15142
     &6716069345D0*t30*t13*t82)*t44+0.1362840444624105D1*t53*t4*t34/t102
     &*(-t57-t60+t74-0.3057340045077133D-4*t35/t3/t105/t15*t40*t1-0.9172
     &020135231398D-4*t38/t39/RB)
      ID_RB_POS=D1VARS(ID_RB)
      D1F(i, ID_RB_POS) = D1F(i, ID_RB_POS) +      t121 
      t1 = PI**2
      t2 = t1*RB
      t3 = t2**(1.D0/3.D0)
      t4 = RB*t3
      t5 = 1/PI
      t6 = dsqrt(GBB)
      t8 = 1/t3
      t10 = 1/RB
      t11 = t6*t8
      t14 = t3**2
      t15 = 1/t14
      t17 = RB**2
      t18 = 1/t17
      t19 = GBB*t15*t18
      t22 = dsqrt(0.4601205204688956D1*t19+1)
      t24 = dlog(0.2145042005343708D1*t11*t10+t22)
      t25 = t10*t24
      t27 = 0.2702765033632353D-1/t6*t8*t25
      t28 = t15*t18
      t31 = dsqrt(1.D0+0.4601205204688955D1*t19)
      t34 = 0.5797544527715597D-1*t28/t31
      t36 = 1/t3/t2
      t37 = t17**2
      t38 = 1/t37
      t41 = dexp(-0.7571335803467249D1*t19)
      t46 = 0.2743D0-0.1508D0*t41
      t53 = 0.5405530067264707D-1*t11*t25
      t54 = GBB**2
      t58 = 1.D0+t53+0.229300503380785D-4*t54*t36*t38
      t68 = t58**2
      t78 = -0.1362840444624105D1*t4*t5*(t27+t34+0.8644628978008849D-1*t
     &36*t38*t41*GBB+0.7571335803467249D-1*t46*t15*t18)/t58+0.1362840444
     &624105D1*t4*t5*(1.D0+t53+0.7571335803467249D-1*t46*GBB*t28)/t68*(t
     &27+t34+0.4586010067615699D-4*GBB*t36*t38)
      ID_GBB_POS=D1VARS(ID_GBB)
      D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) +      t78 
      END IF
      END IF
      END DO
      RETURN
      END
                   
                   
                   
                   
