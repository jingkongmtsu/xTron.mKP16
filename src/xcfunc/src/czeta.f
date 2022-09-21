!--------------------------------------------------------------------------
!                 Reference section for the functional:
!--------------------------------------------------------------------------
! This is the C(zeta,xi) function used in TPSS functional
! Climbing the density functional ladder: Nonempirical 
! meta-generalized gradient approximation designed for 
! molecules and solids
! J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
! Phys. Rev. Lett., 91 (2003) 146401
!
!

      subroutine functional_czeta_close
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
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t52 = (1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)
     &*t6/t37*(0.5D0/t39/t26+0.5D0/t43/t23))**2
      t53 = t52**2
      t55 = (0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.226D1*t9*t3/t10/t5)/t53
       F =     t55 

      ! deriv
      t2 = RA-1.D0*RB
      t3 = RA+RB
      t4 = t3**2
      t5 = 1/t4
      t6 = t2*t5
      t8 = t2**2
      t9 = t4*t3
      t10 = 1/t9
      t14 = t4**2
      t15 = 1/t14
      t18 = t8**2
      t25 = 1/t14/t4
      t28 = t18*t8
      t34 = 1/t3
      t35 = t2*t34
      t37 = 1.D0-1.D0*t35
      t38 = t37**2
      t40 = 1.D0+t35
      t41 = t40**2
      t43 = GAB*t37
      t46 = GAA*t38+GBB*t41-2.D0*t43*t40
      t47 = t46*t5
      t48 = PI**2
      t49 = t48*t3
      t50 = t49**(1.D0/3.D0)
      t51 = t50**2
      t52 = 1/t51
      t53 = t40**(1.D0/3.D0)
      t57 = t37**(1.D0/3.D0)
      t61 = 0.5D0/t53/t40+0.5D0/t57/t37
      t62 = t52*t61
      t65 = 1.D0+0.120187464192284D0*t47*t62
      t66 = t65**2
      t67 = t66**2
      t82 = -1.D0*t34+t6
      t87 = t34-1.D0*t6
      t123 = (0.174D1*t6-0.174D1*t8*t10+0.2D1*t8*t2*t15-0.2D1*t18/t14/t3
     &+0.1356D2*t18*t2*t25-0.1356D2*t28/t14/t9)/t67-4.D0*(0.53D0+0.87D0*
     &t8*t5+0.5D0*t18*t15+0.226D1*t28*t25)/t67/t65*(0.120187464192284D0*
     &(2.D0*GAA*t37*t82+2.D0*GBB*t40*t87-2.D0*GAB*t82*t40-2.D0*t43*t87)*
     &t5*t62-0.2403749283845681D0*t46*t10*t62-0.8012497612818935D-1*t47/
     &t51/t49*t61*t48+0.120187464192284D0*t47*t52*(-0.6666666666666667D0
     &/t53/t41*t87-0.6666666666666667D0/t57/t38*t82))
       D1F( ID_RA_POS) =     t123 
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t38 = 1/t37
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t47 = 0.5D0/t39/t26+0.5D0/t43/t23
      t51 = 1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)*
     &t6*t38*t47
      t52 = t51**2
      t53 = t52**2
      t62 = -0.4807498567691361D0*(0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.22
     &6D1*t9*t3/t10/t5)/t53/t51*t24*t6*t38*t47
       D1F( ID_GAA_POS) =      t62 
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t48 = 1/t37*(0.5D0/t39/t26+0.5D0/t43/t23)
      t51 = 1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)*
     &t6*t48
      t52 = t51**2
      t53 = t52**2
      t61 = 0.9614997135382723D0*(0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.226
     &D1*t9*t3/t10/t5)/t53/t51*t23*t26*t6*t48
       D1F( ID_GAB_POS) =     t61 
      D1F( ID_GBB_POS)= D1F( ID_GAA_POS)
      D1F( ID_RB_POS)= D1F( ID_RA_POS)
      RETURN
      END
                   
                   
      subroutine functional_czeta_open
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
      DOUBLE PRECISION F,D1F(*)
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
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t52 = (1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)
     &*t6/t37*(0.5D0/t39/t26+0.5D0/t43/t23))**2
      t53 = t52**2
      t55 = (0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.226D1*t9*t3/t10/t5)/t53
      F =    t55 

      ! deriv
      t2 = RA-1.D0*RB
      t3 = RA+RB
      t4 = t3**2
      t5 = 1/t4
      t6 = t2*t5
      t8 = t2**2
      t9 = t4*t3
      t10 = 1/t9
      t14 = t4**2
      t15 = 1/t14
      t18 = t8**2
      t25 = 1/t14/t4
      t28 = t18*t8
      t34 = 1/t3
      t35 = t2*t34
      t37 = 1.D0-1.D0*t35
      t38 = t37**2
      t40 = 1.D0+t35
      t41 = t40**2
      t43 = GAB*t37
      t46 = GAA*t38+GBB*t41-2.D0*t43*t40
      t47 = t46*t5
      t48 = PI**2
      t49 = t48*t3
      t50 = t49**(1.D0/3.D0)
      t51 = t50**2
      t52 = 1/t51
      t53 = t40**(1.D0/3.D0)
      t57 = t37**(1.D0/3.D0)
      t61 = 0.5D0/t53/t40+0.5D0/t57/t37
      t62 = t52*t61
      t65 = 1.D0+0.120187464192284D0*t47*t62
      t66 = t65**2
      t67 = t66**2
      t82 = -1.D0*t34+t6
      t87 = t34-1.D0*t6
      t123 = (0.174D1*t6-0.174D1*t8*t10+0.2D1*t8*t2*t15-0.2D1*t18/t14/t3
     &+0.1356D2*t18*t2*t25-0.1356D2*t28/t14/t9)/t67-4.D0*(0.53D0+0.87D0*
     &t8*t5+0.5D0*t18*t15+0.226D1*t28*t25)/t67/t65*(0.120187464192284D0*
     &(2.D0*GAA*t37*t82+2.D0*GBB*t40*t87-2.D0*GAB*t82*t40-2.D0*t43*t87)*
     &t5*t62-0.2403749283845681D0*t46*t10*t62-0.8012497612818935D-1*t47/
     &t51/t49*t61*t48+0.120187464192284D0*t47*t52*(-0.6666666666666667D0
     &/t53/t41*t87-0.6666666666666667D0/t57/t38*t82))
       D1F( ID_RA_POS) =     t123 
      t2 = RA-1.D0*RB
      t3 = RA+RB
      t4 = t3**2
      t5 = 1/t4
      t6 = t2*t5
      t8 = t2**2
      t9 = t4*t3
      t10 = 1/t9
      t14 = t4**2
      t15 = 1/t14
      t18 = t8**2
      t25 = 1/t14/t4
      t28 = t18*t8
      t34 = 1/t3
      t35 = t2*t34
      t37 = 1.D0-1.D0*t35
      t38 = t37**2
      t40 = 1.D0+t35
      t41 = t40**2
      t43 = GAB*t37
      t46 = GAA*t38+GBB*t41-2.D0*t43*t40
      t47 = t46*t5
      t48 = PI**2
      t49 = t48*t3
      t50 = t49**(1.D0/3.D0)
      t51 = t50**2
      t52 = 1/t51
      t53 = t40**(1.D0/3.D0)
      t57 = t37**(1.D0/3.D0)
      t61 = 0.5D0/t53/t40+0.5D0/t57/t37
      t62 = t52*t61
      t65 = 1.D0+0.120187464192284D0*t47*t62
      t66 = t65**2
      t67 = t66**2
      t81 = t34+t6
      t87 = -1.D0*t34-1.D0*t6
      t123 = (-0.174D1*t6-0.174D1*t8*t10-0.2D1*t8*t2*t15-0.2D1*t18/t14/t
     &3-0.1356D2*t18*t2*t25-0.1356D2*t28/t14/t9)/t67-4.D0*(0.53D0+0.87D0
     &*t8*t5+0.5D0*t18*t15+0.226D1*t28*t25)/t67/t65*(0.120187464192284D0
     &*(2.D0*GAA*t37*t81+2.D0*GBB*t40*t87-2.D0*GAB*t81*t40-2.D0*t43*t87)
     &*t5*t62-0.2403749283845681D0*t46*t10*t62-0.8012497612818935D-1*t47
     &/t51/t49*t61*t48+0.120187464192284D0*t47*t52*(-0.6666666666666667D
     &0/t53/t41*t87-0.6666666666666667D0/t57/t38*t81))
       D1F( ID_RB_POS) =     t123 
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t38 = 1/t37
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t47 = 0.5D0/t39/t26+0.5D0/t43/t23
      t51 = 1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)*
     &t6*t38*t47
      t52 = t51**2
      t53 = t52**2
      t62 = -0.4807498567691361D0*(0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.22
     &6D1*t9*t3/t10/t5)/t53/t51*t24*t6*t38*t47
       D1F( ID_GAA_POS) =      t62 
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t38 = 1/t37
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t47 = 0.5D0/t39/t26+0.5D0/t43/t23
      t51 = 1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)*
     &t6*t38*t47
      t52 = t51**2
      t53 = t52**2
      t62 = -0.4807498567691361D0*(0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.22
     &6D1*t9*t3/t10/t5)/t53/t51*t27*t6*t38*t47
       D1F( ID_GBB_POS) =     t62 
      t2 = RA-1.D0*RB
      t3 = t2**2
      t4 = RA+RB
      t5 = t4**2
      t6 = 1/t5
      t9 = t3**2
      t10 = t5**2
      t21 = t2/t4
      t23 = 1.D0-1.D0*t21
      t24 = t23**2
      t26 = 1.D0+t21
      t27 = t26**2
      t34 = PI**2
      t36 = (t34*t4)**(1.D0/3.D0)
      t37 = t36**2
      t39 = t26**(1.D0/3.D0)
      t43 = t23**(1.D0/3.D0)
      t48 = 1/t37*(0.5D0/t39/t26+0.5D0/t43/t23)
      t51 = 1.D0+0.120187464192284D0*(GAA*t24+GBB*t27-2.D0*GAB*t23*t26)*
     &t6*t48
      t52 = t51**2
      t53 = t52**2
      t61 = 0.9614997135382723D0*(0.53D0+0.87D0*t3*t6+0.5D0*t9/t10+0.226
     &D1*t9*t3/t10/t5)/t53/t51*t23*t26*t6*t48
       D1F( ID_GAB_POS) =   t61 

      ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN
      t3 = (RA-1.D0*RB)**2
      t5 = (RA+RB)**2
      t9 = t3**2
      t10 = t5**2
      t19 = 0.53D0+0.87D0*t3/t5+0.5D0*t9/t10+0.226D1*t9*t3/t10/t5
       F =      t19 

      ! deriv
      t2 = RA-1.D0*RB
      t3 = RA+RB
      t4 = t3**2
      t8 = t2**2
      t9 = t4*t3
      t14 = t4**2
      t18 = t8**2
      t33 = 0.174D1*t2/t4-0.174D1*t8/t9+0.2D1*t8*t2/t14-0.2D1*t18/t14/t3
     &+0.1356D2*t18*t2/t14/t4-0.1356D2*t18*t8/t14/t9
       D1F( ID_RA_POS) =    t33 
      t1 = 0.D0
       D1F( ID_GAA_POS) =   t1 

      ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN
      t3 = (RA-1.D0*RB)**2
      t5 = (RA+RB)**2
      t9 = t3**2
      t10 = t5**2
      t19 = 0.53D0+0.87D0*t3/t5+0.5D0*t9/t10+0.226D1*t9*t3/t10/t5
       F =     t19 

      ! deriv
      t2 = RA-1.D0*RB
      t3 = RA+RB
      t4 = t3**2
      t8 = t2**2
      t9 = t4*t3
      t14 = t4**2
      t18 = t8**2
      t33 = -0.174D1*t2/t4-0.174D1*t8/t9-0.2D1*t8*t2/t14-0.2D1*t18/t14/t
     &3-0.1356D2*t18*t2/t14/t4-0.1356D2*t18*t8/t14/t9
       D1F( ID_RB_POS) =   t33 
      t1 = 0.D0
       D1F( ID_GBB_POS) =     t1 
      END IF
      RETURN
      END
                   
                   
                   
                   
