!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! this is the Gorling-Levy limit of PBE correlation energy
! see the paper of PSTS functional, the Appendix A:
! Density functional with full exact exchange, balanced nonlocality 
! of correlation, and constraint satisfaction
! Perdew, J.P. and Staroverov, V.N. and Tao, J. and Scuseria, G.E.
! note:
! for the constants used in this program, for example; the gamma, beta and 
! omega etc., I just used the value in this Appendix A. We note that they
! may not accurate enough.
!
!

      subroutine functional_glpbe_pw_close
     &(TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives pointwisely.
! We note that this routine is used in the close shell case,
! that means, we have RA=RB, GAA=GAB=GBB,TA=TB and LA=LB.
! INTPUT :
! IDERIV : the order of functional derivatives, 1, 2, or 3
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
      INTEGER  ID_RA_POS,ID_GAA_POS,ID_GAB_POS,ID_GBB_POS,ID_RB_POS 
      DOUBLE PRECISION rhoA
      DOUBLE PRECISION RA, RB 
      DOUBLE PRECISION DRhoAX, DRhoAY, DRhoAZ
      DOUBLE PRECISION GAA, GAB, GBB 
      DOUBLE PRECISION F,D1F(N_FUNC_DERIV_1)
      DOUBLE PRECISION TOL ! tolerance
      DOUBLE PRECISION PI
      
      ! calcualte PI
      PI    = 4.D0*DATAN(1.D0)

      ! initilize variable position information
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR(i) = -1
      END DO
      VAR_INFOR(ID_RHO)   = 1
      VAR_INFOR(ID_GAMMA) = 1
      CALL INIT_FUNC_DERIV_1(VAR_INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      
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
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t41 = t1**2
      t44 = t25**2
      t46 = t28**2
      t47 = t29**2
      t52 = t18**2
      t61 = dlog(1.D0+1/(0.8561538578947555D-2*t1*t22*t25*t28/t29/t34/t1
     &8+0.7329994285711492D-4*t41*t21*t1*t44*t46/t47/t33/t32/t52))
      t64 = -0.3068528194400547D0/t1*t18*t17*t61
      F =    t64 

      ! deriv
      t1 = PI**2
      t2 = 1/t1
      t4 = RA-1.D0*RB
      t5 = RA+RB
      t6 = 1/t5
      t7 = t4*t6
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t20 = t1**(1.D0/3.D0)
      t21 = t20**2
      t22 = t21*t1
      t24 = dexp(-0.1520077282819692D0*t1)
      t25 = t22*t24
      t27 = GAA+GBB+2.D0*GAB
      t28 = t5**2
      t29 = 1/t28
      t30 = t27*t29
      t31 = t1*t5
      t32 = t31**(1.D0/3.D0)
      t33 = t32**2
      t34 = 1/t33
      t35 = 1/t18
      t36 = t34*t35
      t40 = t1**2
      t41 = t20*t1
      t42 = t40*t41
      t43 = t24**2
      t44 = t42*t43
      t45 = t27**2
      t46 = t28**2
      t47 = 1/t46
      t48 = t45*t47
      t50 = 1/t32/t31
      t51 = t18**2
      t52 = 1/t51
      t53 = t50*t52
      t57 = 0.8561538578947555D-2*t25*t30*t36+0.7329994285711492D-4*t44*
     &t48*t53
      t59 = 1.D0+1/t57
      t60 = dlog(t59)
      t62 = t4*t29
      t72 = 0.3333333333333333D0/t9*(t6-1.D0*t62)+0.3333333333333333D0/t
     &14*(-1.D0*t6+t62)
      t76 = t18*t17
      t78 = t57**2
      t133 = -0.9205584583201641D0*t2*t18*t60*t72+0.3068528194400547D0*t
     &2*t76/t78*(-0.1712307715789511D-1*t25*t27/t28/t5*t36-0.57076923859
     &65037D-2*t40*t21*t24*t30/t33/t31*t35-0.1712307715789511D-1*t22*t24
     &*t27*t29*t34/t76*t72-0.2931997714284597D-3*t44*t45/t46/t5*t53-0.97
     &73325714281989D-4*t1*t41*t43*t48/t32/t28*t52-0.2931997714284597D-3
     &*t42*t43*t45*t47*t50/t51/t17*t72)/t59
      D1F( ID_RA_POS) =    t133 
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t26 = t1*t22*t25
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t30 = 1/t29
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t35 = 1/t34
      t36 = 1/t18
      t41 = t1**2
      t44 = t25**2
      t45 = t41*t21*t1*t44
      t46 = t28**2
      t47 = t29**2
      t48 = 1/t47
      t52 = t18**2
      t54 = 1/t33/t32/t52
      t58 = 0.8561538578947555D-2*t26*t28*t30*t35*t36+0.7329994285711492
     &D-4*t45*t46*t48*t54
      t59 = t58**2
      t76 = 0.3068528194400547D0/t1*t18*t17/t59*(0.8561538578947555D-2*t
     &26*t30*t35*t36+0.1465998857142298D-3*t45*t28*t48*t54)/(1.D0+1/t58)
      D1F( ID_GAA_POS) =     t76 
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t26 = t1*t22*t25
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t30 = 1/t29
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t35 = 1/t34
      t36 = 1/t18
      t41 = t1**2
      t44 = t25**2
      t45 = t41*t21*t1*t44
      t46 = t28**2
      t47 = t29**2
      t48 = 1/t47
      t52 = t18**2
      t54 = 1/t33/t32/t52
      t58 = 0.8561538578947555D-2*t26*t28*t30*t35*t36+0.7329994285711492
     &D-4*t45*t46*t48*t54
      t59 = t58**2
      t76 = 0.3068528194400547D0/t1*t18*t17/t59*(0.1712307715789511D-1*t
     &26*t30*t35*t36+0.2931997714284597D-3*t45*t28*t48*t54)/(1.D0+1/t58)
      D1F( ID_GAB_POS) =  t76 
      D1F( ID_GBB_POS)= D1F( ID_GAA_POS)
      D1F( ID_RB_POS)= D1F( ID_RA_POS)
      RETURN
      END
                   
      subroutine functional_glpbe_pw_open
     &(TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,
     & DRhoBX,DRhoBY,DRhoBZ,F,D1F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!--------------------------------------------------------
! This routine is used to calculate the functional and
! functional derivatives pointwisely.
! We note that this routine is used in the open shell case.
! INTPUT :
! IDERIV : the order of functional derivatives, 1, 2, or 3
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
      INTEGER  ID_RA_POS,ID_GAA_POS,ID_GAB_POS,ID_GBB_POS,ID_RB_POS 
      DOUBLE PRECISION rhoA
      DOUBLE PRECISION RA, RB 
      DOUBLE PRECISION DRhoAX, DRhoAY, DRhoAZ
      DOUBLE PRECISION GAA, GAB, GBB 
      DOUBLE PRECISION rhoB
      DOUBLE PRECISION DRhoBX, DRhoBY, DRhoBZ
      DOUBLE PRECISION F,D1F(N_FUNC_DERIV_1)
      DOUBLE PRECISION TOL ! tolerance
      DOUBLE PRECISION PI
      
      ! calcualte PI
      PI    = 4.D0*DATAN(1.D0)

      ! initilize variable position information
      DO I = 1, MAX_VAR_TYPE
         VAR_INFOR(i) = -1
      END DO
      VAR_INFOR(ID_RHO)   = 1
      VAR_INFOR(ID_GAMMA) = 1
      CALL INIT_FUNC_DERIV_1(VAR_INFOR,D1VARS)
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      
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
         RETURN 
      END IF 
             
      IF (RA.GT.TOL .AND. RB.GT.TOL) THEN
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t41 = t1**2
      t44 = t25**2
      t46 = t28**2
      t47 = t29**2
      t52 = t18**2
      t61 = dlog(1.D0+1/(0.8561538578947555D-2*t1*t22*t25*t28/t29/t34/t1
     &8+0.7329994285711492D-4*t41*t21*t1*t44*t46/t47/t33/t32/t52))
      t64 = -0.3068528194400547D0/t1*t18*t17*t61
      F =     t64 

      ! deriv
      t1 = PI**2
      t2 = 1/t1
      t4 = RA-1.D0*RB
      t5 = RA+RB
      t6 = 1/t5
      t7 = t4*t6
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t20 = t1**(1.D0/3.D0)
      t21 = t20**2
      t22 = t21*t1
      t24 = dexp(-0.1520077282819692D0*t1)
      t25 = t22*t24
      t27 = GAA+GBB+2.D0*GAB
      t28 = t5**2
      t29 = 1/t28
      t30 = t27*t29
      t31 = t1*t5
      t32 = t31**(1.D0/3.D0)
      t33 = t32**2
      t34 = 1/t33
      t35 = 1/t18
      t36 = t34*t35
      t40 = t1**2
      t41 = t20*t1
      t42 = t40*t41
      t43 = t24**2
      t44 = t42*t43
      t45 = t27**2
      t46 = t28**2
      t47 = 1/t46
      t48 = t45*t47
      t50 = 1/t32/t31
      t51 = t18**2
      t52 = 1/t51
      t53 = t50*t52
      t57 = 0.8561538578947555D-2*t25*t30*t36+0.7329994285711492D-4*t44*
     &t48*t53
      t59 = 1.D0+1/t57
      t60 = dlog(t59)
      t62 = t4*t29
      t72 = 0.3333333333333333D0/t9*(t6-1.D0*t62)+0.3333333333333333D0/t
     &14*(-1.D0*t6+t62)
      t76 = t18*t17
      t78 = t57**2
      t133 = -0.9205584583201641D0*t2*t18*t60*t72+0.3068528194400547D0*t
     &2*t76/t78*(-0.1712307715789511D-1*t25*t27/t28/t5*t36-0.57076923859
     &65037D-2*t40*t21*t24*t30/t33/t31*t35-0.1712307715789511D-1*t22*t24
     &*t27*t29*t34/t76*t72-0.2931997714284597D-3*t44*t45/t46/t5*t53-0.97
     &73325714281989D-4*t1*t41*t43*t48/t32/t28*t52-0.2931997714284597D-3
     &*t42*t43*t45*t47*t50/t51/t17*t72)/t59
      D1F( ID_RA_POS) =    t133 
      t1 = PI**2
      t2 = 1/t1
      t4 = RA-1.D0*RB
      t5 = RA+RB
      t6 = 1/t5
      t7 = t4*t6
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t20 = t1**(1.D0/3.D0)
      t21 = t20**2
      t22 = t21*t1
      t24 = dexp(-0.1520077282819692D0*t1)
      t25 = t22*t24
      t27 = GAA+GBB+2.D0*GAB
      t28 = t5**2
      t29 = 1/t28
      t30 = t27*t29
      t31 = t1*t5
      t32 = t31**(1.D0/3.D0)
      t33 = t32**2
      t34 = 1/t33
      t35 = 1/t18
      t36 = t34*t35
      t40 = t1**2
      t41 = t20*t1
      t42 = t40*t41
      t43 = t24**2
      t44 = t42*t43
      t45 = t27**2
      t46 = t28**2
      t47 = 1/t46
      t48 = t45*t47
      t50 = 1/t32/t31
      t51 = t18**2
      t52 = 1/t51
      t53 = t50*t52
      t57 = 0.8561538578947555D-2*t25*t30*t36+0.7329994285711492D-4*t44*
     &t48*t53
      t59 = 1.D0+1/t57
      t60 = dlog(t59)
      t63 = t4*t29
      t72 = 0.3333333333333333D0/t9*(-1.D0*t6-1.D0*t63)+0.33333333333333
     &33D0/t14*(t6+t63)
      t76 = t18*t17
      t78 = t57**2
      t133 = -0.9205584583201641D0*t2*t18*t60*t72+0.3068528194400547D0*t
     &2*t76/t78*(-0.1712307715789511D-1*t25*t27/t28/t5*t36-0.57076923859
     &65037D-2*t40*t21*t24*t30/t33/t31*t35-0.1712307715789511D-1*t22*t24
     &*t27*t29*t34/t76*t72-0.2931997714284597D-3*t44*t45/t46/t5*t53-0.97
     &73325714281989D-4*t1*t41*t43*t48/t32/t28*t52-0.2931997714284597D-3
     &*t42*t43*t45*t47*t50/t51/t17*t72)/t59
      D1F( ID_RB_POS) =    t133 
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t26 = t1*t22*t25
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t30 = 1/t29
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t35 = 1/t34
      t36 = 1/t18
      t41 = t1**2
      t44 = t25**2
      t45 = t41*t21*t1*t44
      t46 = t28**2
      t47 = t29**2
      t48 = 1/t47
      t52 = t18**2
      t54 = 1/t33/t32/t52
      t58 = 0.8561538578947555D-2*t26*t28*t30*t35*t36+0.7329994285711492
     &D-4*t45*t46*t48*t54
      t59 = t58**2
      t76 = 0.3068528194400547D0/t1*t18*t17/t59*(0.8561538578947555D-2*t
     &26*t30*t35*t36+0.1465998857142298D-3*t45*t28*t48*t54)/(1.D0+1/t58)
      D1F( ID_GAA_POS) =     t76 
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t26 = t1*t22*t25
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t30 = 1/t29
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t35 = 1/t34
      t36 = 1/t18
      t41 = t1**2
      t44 = t25**2
      t45 = t41*t21*t1*t44
      t46 = t28**2
      t47 = t29**2
      t48 = 1/t47
      t52 = t18**2
      t54 = 1/t33/t32/t52
      t58 = 0.8561538578947555D-2*t26*t28*t30*t35*t36+0.7329994285711492
     &D-4*t45*t46*t48*t54
      t59 = t58**2
      t76 = 0.3068528194400547D0/t1*t18*t17/t59*(0.8561538578947555D-2*t
     &26*t30*t35*t36+0.1465998857142298D-3*t45*t28*t48*t54)/(1.D0+1/t58)
      D1F( ID_GBB_POS) =   t76 
      t1 = PI**2
      t5 = RA+RB
      t7 = (RA-1.D0*RB)/t5
      t9 = (1.D0+t7)**(1.D0/3.D0)
      t10 = t9**2
      t14 = (1.D0-1.D0*t7)**(1.D0/3.D0)
      t15 = t14**2
      t17 = 0.5D0*t10+0.5D0*t15
      t18 = t17**2
      t21 = t1**(1.D0/3.D0)
      t22 = t21**2
      t25 = dexp(-0.1520077282819692D0*t1)
      t26 = t1*t22*t25
      t28 = GAA+GBB+2.D0*GAB
      t29 = t5**2
      t30 = 1/t29
      t32 = t1*t5
      t33 = t32**(1.D0/3.D0)
      t34 = t33**2
      t35 = 1/t34
      t36 = 1/t18
      t41 = t1**2
      t44 = t25**2
      t45 = t41*t21*t1*t44
      t46 = t28**2
      t47 = t29**2
      t48 = 1/t47
      t52 = t18**2
      t54 = 1/t33/t32/t52
      t58 = 0.8561538578947555D-2*t26*t28*t30*t35*t36+0.7329994285711492
     &D-4*t45*t46*t48*t54
      t59 = t58**2
      t76 = 0.3068528194400547D0/t1*t18*t17/t59*(0.1712307715789511D-1*t
     &26*t30*t35*t36+0.2931997714284597D-3*t45*t28*t48*t54)/(1.D0+1/t58)
      D1F( ID_GAB_POS) =   t76 

      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t7 = dexp(-0.1520077282819692D0*t1)
      t9 = RA**2
      t12 = t1*RA
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t19 = t1**2
      t22 = t7**2
      t24 = GAA**2
      t25 = t9**2
      t36 = dlog(1.D0+1/(0.1359059534668767D-1*t1*t4*t7*GAA/t9/t14+0.184
     &7042819235409D-3*t19*t3*t1*t22*t24/t25/t13/t12))
      t39 = -0.1534264097200273D0/t1*t36
      F =     t39 

      ! deriv
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t7 = dexp(-0.1520077282819692D0*t1)
      t8 = t1*t4*t7
      t9 = RA**2
      t11 = GAA/t9
      t12 = t1*RA
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t15 = 1/t14
      t19 = t1**2
      t20 = t3*t1
      t22 = t7**2
      t23 = t19*t20*t22
      t24 = GAA**2
      t25 = t9**2
      t27 = t24/t25
      t29 = 1/t13/t12
      t33 = 0.1359059534668767D-1*t8*t11*t15+0.1847042819235409D-3*t23*t
     &27*t29
      t34 = t33**2
      t71 = 0.1534264097200273D0/t1/t34*(-0.2718119069337535D-1*t8*GAA/t
     &9/RA*t15-0.9060396897791782D-2*t19*t4*t7*t11/t14/t12-0.73881712769
     &41635D-3*t23*t24/t25/RA*t29-0.2462723758980545D-3*t1*t20*t22*t27/t
     &13/t9)/(1.D0+1/t33)
      D1F( ID_RA_POS) =    t71 
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t5 = t1*t4
      t7 = dexp(-0.1520077282819692D0*t1)
      t9 = RA**2
      t10 = 1/t9
      t12 = t1*RA
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t15 = 1/t14
      t19 = t1**2
      t22 = t7**2
      t23 = t19*t3*t1*t22
      t24 = GAA**2
      t25 = t9**2
      t26 = 1/t25
      t29 = 1/t13/t12
      t33 = 0.1359059534668767D-1*t5*t7*GAA*t10*t15+0.1847042819235409D-
     &3*t23*t24*t26*t29
      t34 = t33**2
      t51 = 0.1534264097200273D0/t1/t34*(0.1359059534668767D-1*t5*t7*t10
     &*t15+0.3694085638470818D-3*t23*GAA*t26*t29)/(1.D0+1/t33)
      D1F( ID_GAA_POS) =     t51 

      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t7 = dexp(-0.1520077282819692D0*t1)
      t9 = RB**2
      t12 = t1*RB
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t19 = t1**2
      t22 = t7**2
      t24 = GBB**2
      t25 = t9**2
      t36 = dlog(1.D0+1/(0.1359059534668767D-1*t1*t4*t7*GBB/t9/t14+0.184
     &7042819235409D-3*t19*t3*t1*t22*t24/t25/t13/t12))
      t39 = -0.1534264097200273D0/t1*t36
      F =       t39 

      ! deriv
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t7 = dexp(-0.1520077282819692D0*t1)
      t8 = t1*t4*t7
      t9 = RB**2
      t11 = GBB/t9
      t12 = t1*RB
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t15 = 1/t14
      t19 = t1**2
      t20 = t3*t1
      t22 = t7**2
      t23 = t19*t20*t22
      t24 = GBB**2
      t25 = t9**2
      t27 = t24/t25
      t29 = 1/t13/t12
      t33 = 0.1359059534668767D-1*t8*t11*t15+0.1847042819235409D-3*t23*t
     &27*t29
      t34 = t33**2
      t71 = 0.1534264097200273D0/t1/t34*(-0.2718119069337535D-1*t8*GBB/t
     &9/RB*t15-0.9060396897791782D-2*t19*t4*t7*t11/t14/t12-0.73881712769
     &41635D-3*t23*t24/t25/RB*t29-0.2462723758980545D-3*t1*t20*t22*t27/t
     &13/t9)/(1.D0+1/t33)
      D1F( ID_RB_POS) =    t71 
      t1 = PI**2
      t3 = t1**(1.D0/3.D0)
      t4 = t3**2
      t5 = t1*t4
      t7 = dexp(-0.1520077282819692D0*t1)
      t9 = RB**2
      t10 = 1/t9
      t12 = t1*RB
      t13 = t12**(1.D0/3.D0)
      t14 = t13**2
      t15 = 1/t14
      t19 = t1**2
      t22 = t7**2
      t23 = t19*t3*t1*t22
      t24 = GBB**2
      t25 = t9**2
      t26 = 1/t25
      t29 = 1/t13/t12
      t33 = 0.1359059534668767D-1*t5*t7*GBB*t10*t15+0.1847042819235409D-
     &3*t23*t24*t26*t29
      t34 = t33**2
      t51 = 0.1534264097200273D0/t1/t34*(0.1359059534668767D-1*t5*t7*t10
     &*t15+0.3694085638470818D-3*t23*GBB*t26*t29)/(1.D0+1/t33)
      D1F( ID_GBB_POS) =   t51 
      END IF
      RETURN
      END
                   
                   
                   
                   
