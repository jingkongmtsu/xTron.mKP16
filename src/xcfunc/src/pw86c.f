!
! This file is taken from Density functional repository
! ftp://ftp.dl.ac.uk/qcg/dft_library/contents.html
! this code directly use the GAA,GAB and GBB
! therefore we form it outside this function
!
!

      subroutine functional_pw86c
     &(INFOR,NG,NDEN,rhoA,rhoB,sigmaaa,sigmaab,sigmabb,F,D1F)
      IMPLICIT NONE
#include "fderiv1.inc"
      INTEGER INFOR(*)
      INTEGER NG,NDEN,ideriv,I
      DOUBLE PRECISION TOL
      DOUBLE PRECISION rhoA(NG),rhoB(NG)
      DOUBLE PRECISION sigmaaa(NG),sigmaab(NG),sigmabb(NG)
      DOUBLE PRECISION F(NG), D1F(NG,*)
      INTEGER  ID_RA_POS,ID_RB_POS,ID_GAA_POS,ID_GAB_POS,ID_GBB_POS
      INTEGER  D1VARS(N_FUNC_DERIV_1)

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! constants
      ideriv = 1

      ! now we need to figure out the pointer of fortran
      ID_RA_POS  = D1VARS(ID_RA)
      ID_RB_POS  = D1VARS(ID_RB)
      ID_GAA_POS = D1VARS(ID_GAA)
      ID_GAB_POS = D1VARS(ID_GAB)
      ID_GBB_POS = D1VARS(ID_GBB)
      
      IF (NDEN .EQ. 1) THEN 
      CALL rks_c_p86
     &(ideriv,NG,rhoA,sigmaaa,F(1),D1F(1,ID_RA_POS),D1F(1,ID_GAA_POS))
      ELSE  
      CALL uks_c_p86
     &(ideriv,NG,rhoA,rhoB,sigmaaa,sigmabb,sigmaab,
     & F(1),D1F(1,ID_RA_POS),D1F(1,ID_RB_POS),D1F(1,ID_GAA_POS),
     & D1F(1,ID_GBB_POS),D1F(1,ID_GAB_POS))
      END IF  

      ! finally, need to complete the derivatives array
      IF (NDEN .EQ. 1) THEN 
         DO I = 1, NG
            D1F(I,ID_RB_POS)  = D1F(I,ID_RA_POS)
            D1F(I,ID_GAB_POS) = D1F(I,ID_GAA_POS)
            D1F(I,ID_GBB_POS) = D1F(I,ID_GAA_POS)
         END DO
      END IF  

      RETURN  
      END     

c
c-----------------------------------------------------------------------
c
      real*8 function piecewise(o,a,b)
c
c     Fix for a Maple problem. The Maple Fortran generator is not
c     able to transform its piecewise command to Fortran. Therefore
c     this function was written to take care of this.
c
      implicit none
      logical o
      real*8 a, b
      if (o) then
         piecewise = a
      else
         piecewise = b
      endif
      return
      end
c:C_P86subrstart

c    Generated: Tue Mar  9 13:39:23 GMT 2004

      subroutine uks_c_p86
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      parameter(tol=1.0d-20)
      logical t7
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhob**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhob)
      t41 = rhob**(1.D0/3.D0)
      zk(i) = rhob*t18+0.7937005259840997D0*t39*t31*sigmabb/t41/rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhoa**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhoa)
      t41 = rhoa**(1.D0/3.D0)
      zk(i) = rhoa*t18+0.7937005259840997D0*t39*t31*sigmaaa/t41/rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t8 = t4**(1.D0/6.D0)
      t13 = 0.1423D0/(1.D0+0.8292885914166397D0*t8+0.20682485366586D0
     #*t5)
      t22 = (rhoa-1.D0*rhob)*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t46 = piecewise(1.D0 .le. t6,-t13+0.1923661050931536D1*(
     #-0.843D-1/(1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5)
     #+t13)*t30,0.311D-1*t34-0.48D-1+0.12407009817988D-2*t36
     #-0.719606569443304D-2*t5+0.1923661050931536D1*(-0.1555D-1*t34
     #+0.211D-1-0.80645563816922D-3*t36+0.421838333811592D-2*t5)*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t59 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t5
     #+0.2843543831490386D-5*t50)/(1.D0+0.5411317332115466D1*t5
     #+0.1816419932959077D0*t50+0.1763993811759022D-1*t4)
      t62 = rho**(1.D0/6.D0)
      t67 = dexp(-0.81290825D-3*t48/t59/t62/rho)
      t69 = rho**(1.D0/3.D0)
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t84 = dsqrt(t76*t74+t81*t79)
      zk(i) = rho*t46+0.7937005259840997D0*t67*t59*sigma/t69/rho/t84
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhob**(1.D0/6.D0)
      t36 = 1/t34/rhob
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhob**(1.D0/3.D0)
      t43 = 1/t41/rhob
      t44 = sigmabb*t43
      zk(i) = rhob*t18+0.7937005259840997D0*t40*t44
      vrhoa(i) = 0.D0
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhob**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhob(i) = t18+rhob*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmabb
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmabb/t41/t53
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t113 = dsqrt(rhob)
      vsigmabb(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhoa**(1.D0/6.D0)
      t36 = 1/t34/rhoa
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhoa**(1.D0/3.D0)
      t43 = 1/t41/rhoa
      t44 = sigmaaa*t43
      zk(i) = rhoa*t18+0.7937005259840997D0*t40*t44
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhoa**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhoa(i) = t18+rhoa*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmaaa
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmaaa/t41/t53
      vrhob(i) = 0.D0
      t113 = dsqrt(rhoa)
      vsigmaaa(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t7 = 1.D0 .le. t6
      t8 = t4**(1.D0/6.D0)
      t11 = 1.D0+0.8292885914166397D0*t8+0.20682485366586D0*t5
      t13 = 0.1423D0/t11
      t16 = 1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5
      t19 = -0.843D-1/t16+t13
      t21 = rhoa-1.D0*rhob
      t22 = t21*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t42 = -0.1555D-1*t34+0.211D-1-0.80645563816922D-3*t36
     #+0.421838333811592D-2*t5
      t46 = piecewise(t7,-t13+0.1923661050931536D1*t19*t30,0.311D-1
     #*t34-0.48D-1+0.12407009817988D-2*t36-0.719606569443304D-2*t5
     #+0.1923661050931536D1*t42*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t52 = 0.2568D-2+0.1443307452126544D-1*t5+0.2843543831490386D-5*t50
      t56 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t50
     #+0.1763993811759022D-1*t4
      t57 = 1/t56
      t59 = 0.1667D-2+t52*t57
      t61 = t48/t59
      t62 = rho**(1.D0/6.D0)
      t64 = 1/t62/rho
      t67 = dexp(-0.81290825D-3*t61*t64)
      t68 = t67*t59
      t69 = rho**(1.D0/3.D0)
      t71 = 1/t69/rho
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t83 = t76*t74+t81*t79
      t84 = dsqrt(t83)
      t85 = 1/t84
      t86 = sigma*t71*t85
      zk(i) = rho*t46+0.7937005259840997D0*t68*t86
      t93 = 0.1333333333333333D1*t24*t4-0.1333333333333333D1*t28*t4
      t98 = piecewise(t7,0.1923661050931536D1*t19
     #*t93,0.1923661050931536D1*t42*t93)
      t100 = t68*sigma
      t103 = t71/t84/t83
      t108 = 0.8333333333333333D0*t76*t4-0.8333333333333333D0*t81*t4
      t112 = t11**2
      t114 = t8**2
      t115 = t114**2
      t118 = rho**2
      t119 = 1/t118
      t120 = 1/t115/t8*t119
      t122 = 1/t50
      t123 = t122*t119
      t127 = 0.1423D0/t112*(-0.1382147652361066D0*t120
     #-0.6894161788861999D-1*t123)
      t128 = t16**2
      t144 = -0.1333333333333333D1*t24*t21*t119+0.1333333333333333D1
     #*t28*t21*t119
      t150 = t122*t34*t119
      t152 = t5*t4
      t165 = piecewise(t7,t127+0.1923661050931536D1*(0.843D-1/t128*(
     #-0.1835293601259385D0*t120-0.5399117105794445D-1*t123)-t127)*t30
     #+0.1923661050931536D1*t19*t144,-0.1036666666666667D-1*t4
     #-0.4135669939329333D-3*t150-0.4135669939329333D-3*t152
     #+0.2398688564811013D-2*t123+0.1923661050931536D1*
     #(0.5183333333333333D-2*t4+0.2688185460564067D-3*t150
     #+0.2688185460564067D-3*t152-0.1406127779371973D-2*t123)*t30
     #+0.1923661050931536D1*t42*t144)
      t166 = rho*t165
      t167 = t59**2
      t172 = 1/t5*t119
      t176 = t56**2
      t185 = (-0.4811024840421814D-2*t123-0.1895695887660258D-5*t172)
     #*t57-1.D0*t52/t176*(-0.1803772444038489D1*t123
     #-0.1210946621972718D0*t172-0.1763993811759022D-1*t119)
      t197 = 0.7937005259840997D0*(0.81290825D-3*t48/t167*t64*t185
     #+0.9483929583333333D-3*t61/t62/t118)*t67*t59*t86
      t200 = 0.7937005259840997D0*t67*t185*t86
      t206 = 0.10582673679788D1*t68*sigma/t69/t118*t85
      t216 = 0.3968502629920499D0*t100*t103*(-0.8333333333333333D0*t76
     #*t21*t119+0.8333333333333333D0*t81*t21*t119)
      vrhoa(i) = rho*t98-0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t217 = -t93
      t222 = piecewise(t7,0.1923661050931536D1*t19
     #*t217,0.1923661050931536D1*t42*t217)
      vrhob(i) = rho*t222+0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t228 = dsqrt(rho)
      t233 = t48/t228/t118*t67*t85
      t236 = t68*t71*t85
      vsigmaaa(i) = -0.322602852800907D-3*t233+0.7937005259840997D0*t236
      vsigmaab(i) = -0.645205705601814D-3*t233+0.1587401051968199D1*t236
      vsigmabb(i) = vsigmaaa(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_p86(ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      logical t5
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.1423D0/(1.D0
     #+0.8292885914166397D0*t6+0.20682485366586D0*t3),0.311D-1*t12
     #-0.48D-1+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rho**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rho)
      t41 = rho**(1.D0/3.D0)
      zk(i) = rho*t18+0.1D1*t39*t31*sigma/t41/rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t5,-0.1423D0/t9,0.311D-1*t12-0.48D-1
     #+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t49 = t9**2
      t51 = t6**2
      t52 = t51**2
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t73 = piecewise(t5,0.1423D0/t49*(-0.1382147652361066D0/t52/t6
     #*t56-0.6894161788861999D-1*t60),-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t59*t12*t56-0.4135669939329333D-3*t3*t2
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t80 = 1/t3*t56
      t84 = t28**2
      t93 = (-0.4811024840421814D-2*t60-0.1895695887660258D-5*t80)*t29
     #-1.D0*t24/t84*(-0.1803772444038489D1*t60-0.1210946621972718D0
     #*t80-0.1763993811759022D-1*t56)
      vrhoa(i) = rho*t47+t18+rho*t73+0.1D1*(0.81290825D-3*t20/t75*t36
     #*t93+0.9483929583333333D-3*t33/t34/t55)*t39*t31*sigma*t43+0.1D1
     #*t39*t93*t44-0.1333333333333333D1*t40*sigma/t41/t55
      t115 = dsqrt(rho)
      vsigmaaa(i) = -0.16258165D-2*t20/t115/t55*t39+0.4D1*t40*t43
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:C_P86subrend

c    Generated: Tue Mar  9 13:39:23 GMT 2004

      subroutine uks_c_p86_d2
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      logical t7
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhob**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhob)
      t41 = rhob**(1.D0/3.D0)
      zk(i) = rhob*t18+0.7937005259840997D0*t39*t31*sigmabb/t41/rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhoa**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhoa)
      t41 = rhoa**(1.D0/3.D0)
      zk(i) = rhoa*t18+0.7937005259840997D0*t39*t31*sigmaaa/t41/rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t8 = t4**(1.D0/6.D0)
      t13 = 0.1423D0/(1.D0+0.8292885914166397D0*t8+0.20682485366586D0
     #*t5)
      t22 = (rhoa-1.D0*rhob)*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t46 = piecewise(1.D0 .le. t6,-t13+0.1923661050931536D1*(
     #-0.843D-1/(1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5)
     #+t13)*t30,0.311D-1*t34-0.48D-1+0.12407009817988D-2*t36
     #-0.719606569443304D-2*t5+0.1923661050931536D1*(-0.1555D-1*t34
     #+0.211D-1-0.80645563816922D-3*t36+0.421838333811592D-2*t5)*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t59 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t5
     #+0.2843543831490386D-5*t50)/(1.D0+0.5411317332115466D1*t5
     #+0.1816419932959077D0*t50+0.1763993811759022D-1*t4)
      t62 = rho**(1.D0/6.D0)
      t67 = dexp(-0.81290825D-3*t48/t59/t62/rho)
      t69 = rho**(1.D0/3.D0)
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t84 = dsqrt(t76*t74+t81*t79)
      zk(i) = rho*t46+0.7937005259840997D0*t67*t59*sigma/t69/rho/t84
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhob**(1.D0/6.D0)
      t36 = 1/t34/rhob
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhob**(1.D0/3.D0)
      t43 = 1/t41/rhob
      t44 = sigmabb*t43
      zk(i) = rhob*t18+0.7937005259840997D0*t40*t44
      vrhoa(i) = 0.D0
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhob**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhob(i) = t18+rhob*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmabb
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmabb/t41/t53
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t113 = dsqrt(rhob)
      vsigmabb(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhoa**(1.D0/6.D0)
      t36 = 1/t34/rhoa
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhoa**(1.D0/3.D0)
      t43 = 1/t41/rhoa
      t44 = sigmaaa*t43
      zk(i) = rhoa*t18+0.7937005259840997D0*t40*t44
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhoa**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhoa(i) = t18+rhoa*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmaaa
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmaaa/t41/t53
      vrhob(i) = 0.D0
      t113 = dsqrt(rhoa)
      vsigmaaa(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t7 = 1.D0 .le. t6
      t8 = t4**(1.D0/6.D0)
      t11 = 1.D0+0.8292885914166397D0*t8+0.20682485366586D0*t5
      t13 = 0.1423D0/t11
      t16 = 1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5
      t19 = -0.843D-1/t16+t13
      t21 = rhoa-1.D0*rhob
      t22 = t21*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t42 = -0.1555D-1*t34+0.211D-1-0.80645563816922D-3*t36
     #+0.421838333811592D-2*t5
      t46 = piecewise(t7,-t13+0.1923661050931536D1*t19*t30,0.311D-1
     #*t34-0.48D-1+0.12407009817988D-2*t36-0.719606569443304D-2*t5
     #+0.1923661050931536D1*t42*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t52 = 0.2568D-2+0.1443307452126544D-1*t5+0.2843543831490386D-5*t50
      t56 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t50
     #+0.1763993811759022D-1*t4
      t57 = 1/t56
      t59 = 0.1667D-2+t52*t57
      t61 = t48/t59
      t62 = rho**(1.D0/6.D0)
      t64 = 1/t62/rho
      t67 = dexp(-0.81290825D-3*t61*t64)
      t68 = t67*t59
      t69 = rho**(1.D0/3.D0)
      t71 = 1/t69/rho
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t83 = t76*t74+t81*t79
      t84 = dsqrt(t83)
      t85 = 1/t84
      t86 = sigma*t71*t85
      zk(i) = rho*t46+0.7937005259840997D0*t68*t86
      t93 = 0.1333333333333333D1*t24*t4-0.1333333333333333D1*t28*t4
      t98 = piecewise(t7,0.1923661050931536D1*t19
     #*t93,0.1923661050931536D1*t42*t93)
      t100 = t68*sigma
      t103 = t71/t84/t83
      t108 = 0.8333333333333333D0*t76*t4-0.8333333333333333D0*t81*t4
      t112 = t11**2
      t114 = t8**2
      t115 = t114**2
      t118 = rho**2
      t119 = 1/t118
      t120 = 1/t115/t8*t119
      t122 = 1/t50
      t123 = t122*t119
      t127 = 0.1423D0/t112*(-0.1382147652361066D0*t120
     #-0.6894161788861999D-1*t123)
      t128 = t16**2
      t144 = -0.1333333333333333D1*t24*t21*t119+0.1333333333333333D1
     #*t28*t21*t119
      t150 = t122*t34*t119
      t152 = t5*t4
      t165 = piecewise(t7,t127+0.1923661050931536D1*(0.843D-1/t128*(
     #-0.1835293601259385D0*t120-0.5399117105794445D-1*t123)-t127)*t30
     #+0.1923661050931536D1*t19*t144,-0.1036666666666667D-1*t4
     #-0.4135669939329333D-3*t150-0.4135669939329333D-3*t152
     #+0.2398688564811013D-2*t123+0.1923661050931536D1*
     #(0.5183333333333333D-2*t4+0.2688185460564067D-3*t150
     #+0.2688185460564067D-3*t152-0.1406127779371973D-2*t123)*t30
     #+0.1923661050931536D1*t42*t144)
      t166 = rho*t165
      t167 = t59**2
      t172 = 1/t5*t119
      t176 = t56**2
      t185 = (-0.4811024840421814D-2*t123-0.1895695887660258D-5*t172)
     #*t57-1.D0*t52/t176*(-0.1803772444038489D1*t123
     #-0.1210946621972718D0*t172-0.1763993811759022D-1*t119)
      t197 = 0.7937005259840997D0*(0.81290825D-3*t48/t167*t64*t185
     #+0.9483929583333333D-3*t61/t62/t118)*t67*t59*t86
      t200 = 0.7937005259840997D0*t67*t185*t86
      t206 = 0.10582673679788D1*t68*sigma/t69/t118*t85
      t216 = 0.3968502629920499D0*t100*t103*(-0.8333333333333333D0*t76
     #*t21*t119+0.8333333333333333D0*t81*t21*t119)
      vrhoa(i) = rho*t98-0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t217 = -t93
      t222 = piecewise(t7,0.1923661050931536D1*t19
     #*t217,0.1923661050931536D1*t42*t217)
      vrhob(i) = rho*t222+0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t228 = dsqrt(rho)
      t233 = t48/t228/t118*t67*t85
      t236 = t68*t71*t85
      vsigmaaa(i) = -0.322602852800907D-3*t233+0.7937005259840997D0*t236
      vsigmaab(i) = -0.645205705601814D-3*t233+0.1587401051968199D1*t236
      vsigmabb(i) = vsigmaaa(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rhob**(1.D0/6.D0)
      t36 = 1/t34/rhob
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhob**(1.D0/3.D0)
      t43 = 1/t41/rhob
      t44 = sigmabb*t43
      zk(i) = rhob*t18+0.7937005259840997D0*t40*t44
      vrhoa(i) = 0.D0
      t47 = t9**2
      t48 = 1/t47
      t49 = t6**2
      t50 = t49**2
      t51 = t50*t6
      t52 = 1/t51
      t53 = rhob**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t60 = -0.1835293601259385D0*t52*t54-0.5399117105794445D-1*t58
      t64 = t57*t12
      t67 = t3*t2
      t71 = piecewise(t7,0.843D-1*t48*t60,-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t64*t54-0.1447484478765267D-3*t67
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t75 = t20/t73
      t77 = 1/t3
      t78 = t77*t54
      t80 = -0.4811024840421814D-2*t58-0.1895695887660258D-5*t78
      t82 = t28**2
      t83 = 1/t82
      t84 = t24*t83
      t88 = -0.1803772444038489D1*t58-0.1210946621972718D0*t78
     #-0.1763993811759022D-1*t54
      t91 = t80*t29-1.D0*t84*t88
      t96 = 1/t34/t53
      t99 = 0.81290825D-3*t75*t36*t91+0.9483929583333333D-3*t33*t96
      t100 = t99*t39
      t101 = t31*sigmabb
      t102 = t101*t43
      t105 = t39*t91
      t109 = 1/t41/t53
      t110 = sigmabb*t109
      vrhob(i) = t18+rhob*t71+0.7937005259840997D0*t100*t102
     #+0.7937005259840997D0*t105*t44-0.10582673679788D1*t40*t110
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t113 = dsqrt(rhob)
      t115 = 1/t113/t53
      vsigmabb(i) = -0.322602852800907D-3*t20*t115*t39
     #+0.7937005259840997D0*t40*t43
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t124 = t60**2
      t129 = t53**2
      t130 = 1/t129
      t133 = t53*rhob
      t134 = 1/t133
      t138 = 1/t22/t2
      t139 = t138*t130
      t141 = t57*t134
      t158 = piecewise(t7,-0.1686D0/t47/t9*t124+0.843D-1*t48*(
     #-0.1529411334382821D0/t51/t2*t130+0.367058720251877D0*t52*t134
     #-0.3599411403862963D-1*t139+0.1079823421158889D0*t141
     #),0.5183333333333333D-2*t54-0.9649896525101778D-4*t138*t12*t130
     #-0.1888622605627062D-2*t141+0.2894968957530533D-3*t64*t134
     #+0.1447484478765267D-3*t3*t54+0.6617071902926934D-3*t139)
      t163 = t91**2
      t173 = 1/t67*t130
      t175 = t77*t134
      t185 = t88**2
      t196 = (-0.3207349893614542D-2*t139+0.9622049680843627D-2*t141
     #-0.6318986292200858D-6*t173+0.3791391775320515D-5*t175)*t29-2.D0
     #*t80*t83*t88+2.D0*t24/t82/t28*t185-1.D0*t84*(
     #-0.1202514962692326D1*t139+0.3607544888076978D1*t141
     #-0.4036488739909061D-1*t173+0.2421893243945437D0*t175
     #+0.3527987623518044D-1*t134)
      t208 = t99**2
      v2rhob2(i) = 2.D0*t71+rhob*t158+0.7937005259840997D0*(
     #-0.16258165D-2*t20/t73/t31*t36*t163-0.1896785916666667D-2*t75
     #*t96*t91+0.81290825D-3*t75*t36*t196-0.2054851409722222D-2*t33
     #/t34/t133)*t39*t102+0.7937005259840997D0*t208*t39*t102
     #+0.1587401051968199D1*t100*t91*sigmabb*t43-0.2116534735957599D1
     #*t100*t101*t109+0.7937005259840997D0*t39*t196*t44
     #-0.2116534735957599D1*t105*t110+0.2469290525283866D1*t40*sigmabb
     #/t41/t133
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t233 = t41**2
      v2sigmabb2(i) = -0.4839042792013605D-3/t20*t115*t39
     #+0.1311232602576965D-6/t233/t133*t32*t39
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rhoa**(1.D0/6.D0)
      t36 = 1/t34/rhoa
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhoa**(1.D0/3.D0)
      t43 = 1/t41/rhoa
      t44 = sigmaaa*t43
      zk(i) = rhoa*t18+0.7937005259840997D0*t40*t44
      t47 = t9**2
      t48 = 1/t47
      t49 = t6**2
      t50 = t49**2
      t51 = t50*t6
      t52 = 1/t51
      t53 = rhoa**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t60 = -0.1835293601259385D0*t52*t54-0.5399117105794445D-1*t58
      t64 = t57*t12
      t67 = t3*t2
      t71 = piecewise(t7,0.843D-1*t48*t60,-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t64*t54-0.1447484478765267D-3*t67
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t75 = t20/t73
      t77 = 1/t3
      t78 = t77*t54
      t80 = -0.4811024840421814D-2*t58-0.1895695887660258D-5*t78
      t82 = t28**2
      t83 = 1/t82
      t84 = t24*t83
      t88 = -0.1803772444038489D1*t58-0.1210946621972718D0*t78
     #-0.1763993811759022D-1*t54
      t91 = t80*t29-1.D0*t84*t88
      t96 = 1/t34/t53
      t99 = 0.81290825D-3*t75*t36*t91+0.9483929583333333D-3*t33*t96
      t100 = t99*t39
      t101 = t31*sigmaaa
      t102 = t101*t43
      t105 = t39*t91
      t109 = 1/t41/t53
      t110 = sigmaaa*t109
      vrhoa(i) = t18+rhoa*t71+0.7937005259840997D0*t100*t102
     #+0.7937005259840997D0*t105*t44-0.10582673679788D1*t40*t110
      vrhob(i) = 0.D0
      t113 = dsqrt(rhoa)
      t115 = 1/t113/t53
      vsigmaaa(i) = -0.322602852800907D-3*t20*t115*t39
     #+0.7937005259840997D0*t40*t43
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t124 = t60**2
      t129 = t53**2
      t130 = 1/t129
      t133 = t53*rhoa
      t134 = 1/t133
      t138 = 1/t22/t2
      t139 = t138*t130
      t141 = t57*t134
      t158 = piecewise(t7,-0.1686D0/t47/t9*t124+0.843D-1*t48*(
     #-0.1529411334382821D0/t51/t2*t130+0.367058720251877D0*t52*t134
     #-0.3599411403862963D-1*t139+0.1079823421158889D0*t141
     #),0.5183333333333333D-2*t54-0.9649896525101778D-4*t138*t12*t130
     #-0.1888622605627062D-2*t141+0.2894968957530533D-3*t64*t134
     #+0.1447484478765267D-3*t3*t54+0.6617071902926934D-3*t139)
      t163 = t91**2
      t173 = 1/t67*t130
      t175 = t77*t134
      t185 = t88**2
      t196 = (-0.3207349893614542D-2*t139+0.9622049680843627D-2*t141
     #-0.6318986292200858D-6*t173+0.3791391775320515D-5*t175)*t29-2.D0
     #*t80*t83*t88+2.D0*t24/t82/t28*t185-1.D0*t84*(
     #-0.1202514962692326D1*t139+0.3607544888076978D1*t141
     #-0.4036488739909061D-1*t173+0.2421893243945437D0*t175
     #+0.3527987623518044D-1*t134)
      t208 = t99**2
      v2rhoa2(i) = 2.D0*t71+rhoa*t158+0.7937005259840997D0*(
     #-0.16258165D-2*t20/t73/t31*t36*t163-0.1896785916666667D-2*t75
     #*t96*t91+0.81290825D-3*t75*t36*t196-0.2054851409722222D-2*t33
     #/t34/t133)*t39*t102+0.7937005259840997D0*t208*t39*t102
     #+0.1587401051968199D1*t100*t91*sigmaaa*t43-0.2116534735957599D1
     #*t100*t101*t109+0.7937005259840997D0*t39*t196*t44
     #-0.2116534735957599D1*t105*t110+0.2469290525283866D1*t40*sigmaaa
     #/t41/t133
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t233 = t41**2
      v2sigmaaa2(i) = -0.4839042792013605D-3/t20*t115*t39
     #+0.1311232602576965D-6/t233/t133*t32*t39
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t7 = 1.D0 .le. t6
      t8 = t4**(1.D0/6.D0)
      t11 = 1.D0+0.8292885914166397D0*t8+0.20682485366586D0*t5
      t13 = 0.1423D0/t11
      t16 = 1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5
      t19 = -0.843D-1/t16+t13
      t21 = rhoa-1.D0*rhob
      t22 = t21*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t42 = -0.1555D-1*t34+0.211D-1-0.80645563816922D-3*t36
     #+0.421838333811592D-2*t5
      t46 = piecewise(t7,-t13+0.1923661050931536D1*t19*t30,0.311D-1
     #*t34-0.48D-1+0.12407009817988D-2*t36-0.719606569443304D-2*t5
     #+0.1923661050931536D1*t42*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t52 = 0.2568D-2+0.1443307452126544D-1*t5+0.2843543831490386D-5*t50
      t56 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t50
     #+0.1763993811759022D-1*t4
      t57 = 1/t56
      t59 = 0.1667D-2+t52*t57
      t60 = 1/t59
      t61 = t48*t60
      t62 = rho**(1.D0/6.D0)
      t64 = 1/t62/rho
      t67 = dexp(-0.81290825D-3*t61*t64)
      t68 = t67*t59
      t69 = rho**(1.D0/3.D0)
      t71 = 1/t69/rho
      t72 = sigma*t71
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t83 = t76*t74+t81*t79
      t84 = dsqrt(t83)
      t85 = 1/t84
      t86 = t72*t85
      zk(i) = rho*t46+0.7937005259840997D0*t68*t86
      t93 = 0.1333333333333333D1*t24*t4-0.1333333333333333D1*t28*t4
      t98 = piecewise(t7,0.1923661050931536D1*t19
     #*t93,0.1923661050931536D1*t42*t93)
      t100 = t68*sigma
      t102 = 1/t84/t83
      t103 = t71*t102
      t108 = 0.8333333333333333D0*t76*t4-0.8333333333333333D0*t81*t4
      t109 = t103*t108
      t112 = t11**2
      t113 = 1/t112
      t114 = t8**2
      t115 = t114**2
      t116 = t115*t8
      t117 = 1/t116
      t118 = rho**2
      t119 = 1/t118
      t120 = t117*t119
      t122 = 1/t50
      t123 = t122*t119
      t125 = -0.1382147652361066D0*t120-0.6894161788861999D-1*t123
      t127 = 0.1423D0*t113*t125
      t128 = t16**2
      t129 = 1/t128
      t132 = -0.1835293601259385D0*t120-0.5399117105794445D-1*t123
      t135 = 0.843D-1*t129*t132-t127
      t138 = t24*t21
      t141 = t28*t21
      t144 = -0.1333333333333333D1*t138*t119+0.1333333333333333D1*t141
     #*t119
      t149 = t122*t34
      t150 = t149*t119
      t152 = t5*t4
      t159 = 0.5183333333333333D-2*t4+0.2688185460564067D-3*t150
     #+0.2688185460564067D-3*t152-0.1406127779371973D-2*t123
      t165 = piecewise(t7,t127+0.1923661050931536D1*t135*t30
     #+0.1923661050931536D1*t19*t144,-0.1036666666666667D-1*t4
     #-0.4135669939329333D-3*t150-0.4135669939329333D-3*t152
     #+0.2398688564811013D-2*t123+0.1923661050931536D1*t159*t30
     #+0.1923661050931536D1*t42*t144)
      t166 = rho*t165
      t167 = t59**2
      t168 = 1/t167
      t169 = t48*t168
      t171 = 1/t5
      t172 = t171*t119
      t174 = -0.4811024840421814D-2*t123-0.1895695887660258D-5*t172
      t176 = t56**2
      t177 = 1/t176
      t178 = t52*t177
      t182 = -0.1803772444038489D1*t123-0.1210946621972718D0*t172
     #-0.1763993811759022D-1*t119
      t185 = t174*t57-1.D0*t178*t182
      t186 = t64*t185
      t190 = 1/t62/t118
      t193 = 0.81290825D-3*t169*t186+0.9483929583333333D-3*t61*t190
      t194 = t193*t67
      t195 = t194*t59
      t197 = 0.7937005259840997D0*t195*t86
      t198 = t67*t185
      t200 = 0.7937005259840997D0*t198*t86
      t202 = 1/t69/t118
      t204 = sigma*t202*t85
      t206 = 0.10582673679788D1*t68*t204
      t207 = t76*t21
      t210 = t81*t21
      t213 = -0.8333333333333333D0*t207*t119+0.8333333333333333D0*t210
     #*t119
      t214 = t103*t213
      t216 = 0.3968502629920499D0*t100*t214
      vrhoa(i) = rho*t98-0.3968502629920499D0*t100*t109+t46+t166+t197
     #+t200-t206-t216
      t217 = -t93
      t222 = piecewise(t7,0.1923661050931536D1*t19
     #*t217,0.1923661050931536D1*t42*t217)
      t224 = -t108
      t225 = t103*t224
      vrhob(i) = rho*t222-0.3968502629920499D0*t100*t225+t46+t166+t197
     #+t200-t206-t216
      t228 = dsqrt(rho)
      t230 = 1/t228/t118
      t231 = t48*t230
      t232 = t67*t85
      t233 = t231*t232
      t235 = t71*t85
      t236 = t68*t235
      vsigmaaa(i) = -0.322602852800907D-3*t233+0.7937005259840997D0*t236
      vsigmaab(i) = -0.645205705601814D-3*t233+0.1587401051968199D1*t236
      vsigmabb(i) = vsigmaaa(i)
      t240 = 2.D0*t165
      t244 = t185**2
      t252 = 1/t50/t4
      t253 = t118**2
      t254 = 1/t253
      t255 = t252*t254
      t257 = t118*rho
      t258 = 1/t257
      t259 = t122*t258
      t262 = 1/t152*t254
      t264 = t171*t258
      t274 = t182**2
      t285 = (-0.3207349893614542D-2*t255+0.9622049680843627D-2*t259
     #-0.6318986292200858D-6*t262+0.3791391775320515D-5*t264)*t57-2.D0
     #*t174*t177*t182+2.D0*t52/t176/t56*t274-1.D0*t178*(
     #-0.1202514962692326D1*t255+0.3607544888076978D1*t259
     #-0.4036488739909061D-1*t262+0.2421893243945437D0*t264
     #+0.3527987623518044D-1*t258)
      t297 = 0.7937005259840997D0*(-0.16258165D-2*t48/t167/t59*t64
     #*t244-0.1896785916666667D-2*t169*t190*t185+0.81290825D-3*t169
     #*t64*t285-0.2054851409722222D-2*t61/t62/t257)*t67*t59*t86
      t298 = t193**2
      t302 = 0.7937005259840997D0*t298*t67*t59*t86
      t306 = 0.7937005259840997D0*t195*t72*t102*t213
      t308 = 0.2116534735957599D1*t195*t204
      t311 = 0.1587401051968199D1*t194*t185*t86
      t314 = 0.7937005259840997D0*t67*t285*t86
      t315 = t198*sigma
      t317 = 0.7937005259840997D0*t315*t214
      t319 = 0.2116534735957599D1*t198*t204
      t320 = t24**2
      t321 = 1/t320
      t324 = t28**2
      t325 = 1/t324
      t328 = 0.4444444444444444D0*t321*t119+0.4444444444444444D0*t325
     #*t119
      t333 = piecewise(t7,0.1923661050931536D1*t19
     #*t328,0.1923661050931536D1*t42*t328)
      t334 = rho*t333
      t335 = t83**2
      t338 = t71/t84/t335
      t339 = t108**2
      t343 = 1/t75
      t346 = 1/t80
      t349 = 0.2777777777777778D0*t343*t119+0.2777777777777778D0*t346
     #*t119
      t352 = 0.3968502629920499D0*t100*t103*t349
      t353 = t240+t297+t302-t306-t308+t311+t314-t317-t319+t334
     #+0.5952753944880748D0*t100*t338*t339-t352
      t355 = t202*t102
      t358 = 0.10582673679788D1*t100*t355*t213
      t364 = 0.2469290525283866D1*t68*sigma/t69/t257*t85
      t375 = -0.2777777777777778D0*t343*t21*t258-0.8333333333333333D0
     #*t76*t119-0.2777777777777778D0*t346*t21*t258
     #+0.8333333333333333D0*t81*t119
      t377 = t100*t103*t375
      t379 = t315*t109
      t383 = t195*t72*t102*t108
      t397 = -0.4444444444444444D0*t321*t21*t258-0.1333333333333333D1
     #*t24*t119-0.4444444444444444D0*t325*t21*t258
     #+0.1333333333333333D1*t28*t119
      t406 = piecewise(t7,0.1923661050931536D1*t135*t93
     #+0.1923661050931536D1*t19*t397,0.1923661050931536D1*t159*t93
     #+0.1923661050931536D1*t42*t397)
      t407 = rho*t406
      t411 = t100*t338*t213*t108
      t414 = t100*t355*t108
      t416 = t213**2
      t419 = 0.5952753944880748D0*t100*t338*t416
      t420 = t21**2
      t434 = 0.3968502629920499D0*t100*t103*(0.2777777777777778D0*t343
     #*t420*t254+0.1666666666666667D1*t207*t258+0.2777777777777778D0
     #*t346*t420*t254-0.1666666666666667D1*t210*t258)
      t437 = t125**2
      t439 = 0.2846D0/t112/t11*t437
      t442 = 1/t116/t4*t254
      t444 = t117*t258
      t450 = 0.1423D0*t113*(-0.1151789710300888D0*t442
     #+0.2764295304722132D0*t444-0.4596107859241333D-1*t255
     #+0.13788323577724D0*t259)
      t453 = t132**2
      t478 = 0.4444444444444444D0*t321*t420*t254+0.2666666666666667D1
     #*t138*t258+0.4444444444444444D0*t325*t420*t254
     #-0.2666666666666667D1*t141*t258
      t484 = t252*t34*t254
      t487 = t149*t258
      t489 = t5*t119
      t506 = piecewise(t7,-t439+t450+0.1923661050931536D1*(-0.1686D0
     #/t128/t16*t453+0.843D-1*t129*(-0.1529411334382821D0*t442
     #+0.367058720251877D0*t444-0.3599411403862963D-1*t255
     #+0.1079823421158889D0*t259)+t439-t450)*t30+0.3847322101863073D1
     #*t135*t144+0.1923661050931536D1*t19*t478,0.1036666666666667D-1
     #*t119-0.2757113292886222D-3*t484-0.4521665800333405D-2*t259
     #+0.8271339878658667D-3*t487+0.4135669939329333D-3*t489
     #+0.1599125709874009D-2*t255+0.1923661050931536D1*(
     #-0.5183333333333333D-2*t119+0.1792123640376044D-3*t484
     #+0.2633043194706342D-2*t259-0.5376370921128133D-3*t487
     #-0.2688185460564067D-3*t489-0.9374185195813156D-3*t255)*t30
     #+0.3847322101863073D1*t159*t144+0.1923661050931536D1*t42*t478)
      t507 = rho*t506
      t508 = 2.D0*t98+t358+t364-0.7937005259840997D0*t377
     #-0.7937005259840997D0*t379-0.7937005259840997D0*t383+2.D0*t407
     #+0.119055078897615D1*t411+0.10582673679788D1*t414+t419-t434+t507
      v2rhoa2(i) = t353+t508
      t510 = t224**2
      t516 = -t397
      t525 = piecewise(t7,0.1923661050931536D1*t135*t217
     #+0.1923661050931536D1*t19*t516,0.1923661050931536D1*t159*t217
     #+0.1923661050931536D1*t42*t516)
      t526 = rho*t525
      t530 = -t100*t103*t375
      t532 = t315*t225
      t534 = t240+t297+2.D0*t222+0.5952753944880748D0*t100*t338*t510
     #+t302-t306-t308+t311+t314+2.D0*t526-0.7937005259840997D0*t530
     #-0.7937005259840997D0*t532
      t536 = t100*t355*t224
      t540 = t100*t338*t213*t224
      t544 = t195*t72*t102*t224
      t546 = 0.10582673679788D1*t536+0.119055078897615D1*t540
     #-0.7937005259840997D0*t544-t317-t319+t334-t352+t358+t364+t419
     #-t434+t507
      v2rhob2(i) = t534+t546
      t560 = -t328
      t565 = piecewise(t7,0.1923661050931536D1*t19
     #*t560,0.1923661050931536D1*t42*t560)
      t567 = t526-0.3968502629920499D0*t530-0.3968502629920499D0*t532
     #+0.5291336839893998D0*t536+0.5952753944880748D0*t540
     #-0.3968502629920499D0*t544+0.3968502629920499D0*t100*t103*t349
     #+0.5952753944880748D0*t100*t338*t108*t224+t98+rho*t565+t240+t222
     #+t297+t302-t306
      t573 = -t308+t311+t314-t319-t317+t358+t364-0.3968502629920499D0
     #*t377-0.3968502629920499D0*t379-0.3968502629920499D0*t383+t407
     #+0.5952753944880748D0*t411+t507+0.5291336839893998D0*t414-t434
     #+t419
      v2rhoab(i) = t567+t573
      t574 = t67*t102
      t576 = t231*t574*t108
      t578 = t68*t109
      t580 = 1/t48
      t590 = (0.406454125D-3*t580*t168*t186+0.4741964791666667D-3*t580
     #*t60*t190)*t67*t59*t86
      t591 = 0.7937005259840997D0*t590
      t595 = t193*t48*t230*t67*t85
      t596 = 0.322602852800907D-3*t595
      t599 = t194*t59*t71*t85
      t600 = 0.7937005259840997D0*t599
      t603 = t61*t230*t198*t85
      t604 = 0.322602852800907D-3*t603
      t605 = t198*t235
      t606 = 0.7937005259840997D0*t605
      t610 = t48/t228/t257*t232
      t611 = 0.430137137067876D-3*t610
      t613 = t68*t202*t85
      t614 = 0.10582673679788D1*t613
      t616 = t231*t574*t213
      t617 = 0.1613014264004535D-3*t616
      t618 = t68*t214
      t619 = 0.3968502629920499D0*t618
      v2rhoasigmaaa(i) = 0.1613014264004535D-3*t576
     #-0.3968502629920499D0*t578+t591-t596+t600-t604+t606+t611-t614
     #+t617-t619
      t622 = 0.1587401051968199D1*t590
      t623 = 0.645205705601814D-3*t595
      t624 = 0.1587401051968199D1*t599
      t625 = 0.645205705601814D-3*t603
      t626 = 0.1587401051968199D1*t605
      t627 = 0.8602742741357521D-3*t610
      t628 = 0.2116534735957599D1*t613
      t629 = 0.322602852800907D-3*t616
      t630 = 0.7937005259840997D0*t618
      v2rhoasigmaab(i) = 0.322602852800907D-3*t576
     #-0.7937005259840997D0*t578+t622-t623+t624-t625+t626+t627-t628
     #+t629-t630
      v2rhoasigmabb(i) = v2rhoasigmaaa(i)
      t632 = t231*t574*t224
      t634 = t68*t225
      v2rhobsigmaaa(i) = 0.1613014264004535D-3*t632
     #-0.3968502629920499D0*t634+t591-t596+t600-t604+t606+t611-t614
     #+t617-t619
      v2rhobsigmaab(i) = 0.322602852800907D-3*t632
     #-0.7937005259840997D0*t634+t622-t623+t624-t625+t626+t627-t628
     #+t629-t630
      v2rhobsigmabb(i) = v2rhobsigmaaa(i)
      t639 = t580*t230*t232
      t641 = t69**2
      t645 = 1/t641/t257*t60*t232
      v2sigmaaa2(i) = -0.4839042792013605D-3*t639
     #+0.1311232602576965D-6*t645
      v2sigmaaaab(i) = -0.9678085584027211D-3*t639
     #+0.2622465205153929D-6*t645
      v2sigmaaabb(i) = v2sigmaaa2(i)
      v2sigmaab2(i) = -0.1935617116805442D-2*t639
     #+0.5244930410307859D-6*t645
      v2sigmaabbb(i) = v2sigmaaaab(i)
      v2sigmabb2(i) = v2sigmaaabb(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_p86_d2
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      logical t5
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.1423D0/(1.D0
     #+0.8292885914166397D0*t6+0.20682485366586D0*t3),0.311D-1*t12
     #-0.48D-1+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rho**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rho)
      t41 = rho**(1.D0/3.D0)
      zk(i) = rho*t18+0.1D1*t39*t31*sigma/t41/rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t5,-0.1423D0/t9,0.311D-1*t12-0.48D-1
     #+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t49 = t9**2
      t51 = t6**2
      t52 = t51**2
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t73 = piecewise(t5,0.1423D0/t49*(-0.1382147652361066D0/t52/t6
     #*t56-0.6894161788861999D-1*t60),-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t59*t12*t56-0.4135669939329333D-3*t3*t2
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t80 = 1/t3*t56
      t84 = t28**2
      t93 = (-0.4811024840421814D-2*t60-0.1895695887660258D-5*t80)*t29
     #-1.D0*t24/t84*(-0.1803772444038489D1*t60-0.1210946621972718D0
     #*t80-0.1763993811759022D-1*t56)
      vrhoa(i) = rho*t47+t18+rho*t73+0.1D1*(0.81290825D-3*t20/t75*t36
     #*t93+0.9483929583333333D-3*t33/t34/t55)*t39*t31*sigma*t43+0.1D1
     #*t39*t93*t44-0.1333333333333333D1*t40*sigma/t41/t55
      t115 = dsqrt(rho)
      vsigmaaa(i) = -0.16258165D-2*t20/t115/t55*t39+0.4D1*t40*t43
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t11 = 0.1423D0/t9
      t12 = dlog(t4)
      t14 = t3*t12
      t18 = piecewise(t5,-t11,0.311D-1*t12-0.48D-1+0.12407009817988D-2
     #*t14-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t48 = rho*t47
      t49 = t9**2
      t50 = 1/t49
      t51 = t6**2
      t52 = t51**2
      t53 = t52*t6
      t54 = 1/t53
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t62 = -0.1382147652361066D0*t54*t56-0.6894161788861999D-1*t60
      t66 = t59*t12
      t69 = t3*t2
      t73 = piecewise(t5,0.1423D0*t50*t62,-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t66*t56-0.4135669939329333D-3*t69
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t76 = 1/t75
      t77 = t20*t76
      t79 = 1/t3
      t80 = t79*t56
      t82 = -0.4811024840421814D-2*t60-0.1895695887660258D-5*t80
      t84 = t28**2
      t85 = 1/t84
      t86 = t24*t85
      t90 = -0.1803772444038489D1*t60-0.1210946621972718D0*t80
     #-0.1763993811759022D-1*t56
      t93 = t82*t29-1.D0*t86*t90
      t94 = t36*t93
      t98 = 1/t34/t55
      t101 = 0.81290825D-3*t77*t94+0.9483929583333333D-3*t33*t98
      t102 = t101*t39
      t103 = t31*sigma
      t104 = t103*t43
      t107 = t39*t93
      t111 = 1/t41/t55
      t112 = sigma*t111
      vrhoa(i) = t48+t18+rho*t73+0.1D1*t102*t104+0.1D1*t107*t44
     #-0.1333333333333333D1*t40*t112
      t115 = dsqrt(rho)
      t117 = 1/t115/t55
      vsigmaaa(i) = -0.16258165D-2*t20*t117*t39+0.4D1*t40*t43
      t129 = (-0.843D-1/(1.D0+0.1101176160755631D1*t6
     #+0.1619735131738333D0*t3)+t11)*t56
      t135 = (-0.1555D-1*t12+0.211D-1-0.80645563816922D-3*t14
     #+0.421838333811592D-2*t3)*t56
      t137 = piecewise(t5,0.1709920934161366D1
     #*t129,0.1709920934161366D1*t135)
      t139 = t55*rho
      t147 = t62**2
      t152 = t55**2
      t153 = 1/t152
      t156 = 1/t139
      t160 = 1/t22/t2
      t161 = t160*t153
      t163 = t59*t156
      t180 = piecewise(t5,-0.2846D0/t49/t9*t147+0.1423D0*t50*(
     #-0.1151789710300888D0/t53/t2*t153+0.2764295304722132D0*t54*t156
     #-0.4596107859241333D-1*t161+0.13788323577724D0*t163
     #),0.1036666666666667D-1*t56-0.2757113292886222D-3*t160*t12*t153
     #-0.4521665800333405D-2*t163+0.8271339878658667D-3*t66*t156
     #+0.4135669939329333D-3*t3*t56+0.1599125709874009D-2*t161)
      t185 = piecewise(t5,-0.1709920934161366D1*t129,
     #-0.1709920934161366D1*t135)
      t193 = t93**2
      t203 = 1/t69*t153
      t205 = t79*t156
      t215 = t90**2
      t226 = (-0.3207349893614542D-2*t161+0.9622049680843627D-2*t163
     #-0.6318986292200858D-6*t203+0.3791391775320515D-5*t205)*t29-2.D0
     #*t82*t85*t90+2.D0*t24/t84/t28*t215-1.D0*t86*(
     #-0.1202514962692326D1*t161+0.3607544888076978D1*t163
     #-0.4036488739909061D-1*t203+0.2421893243945437D0*t205
     #+0.3527987623518044D-1*t156)
      t247 = t101**2
      v2rhoa2(i) = rho*t137+0.6222222222222222D1*t40*sigma/t41/t139
     #+2.D0*rho*t180+rho*t185-0.5333333333333333D1*t102*t103*t111
     #+0.2D1*(-0.16258165D-2*t20/t75/t31*t36*t193
     #-0.1896785916666667D-2*t77*t98*t93+0.81290825D-3*t77*t36*t226
     #-0.2054851409722222D-2*t33/t34/t139)*t39*t104+0.4D1*t102*t93
     #*sigma*t43+0.2D1*t39*t226*t44-0.5333333333333333D1*t107*t112
     #+0.2D1*t247*t39*t104+4.D0*t48+4.D0*t47+4.D0*t73
      t254 = 1/t20
      t266 = t117*t39
      v2rhoasigmaaa(i) = 0.4D1*(0.406454125D-3*t254*t76*t94
     #+0.4741964791666667D-3*t254*t32*t98)*t39*t104-0.16258165D-2*t101
     #*t20*t266+0.4D1*t102*t31*t43-0.16258165D-2*t33*t266*t93+0.4D1
     #*t107*t43+0.2167755333333333D-2*t20/t115/t139*t39
     #-0.5333333333333333D1*t40*t111
      t287 = t41**2
      v2sigmaaa2(i) = -0.9754899D-2*t254*t117*t39+0.264327929167225D-5
     #/t287/t139*t32*t39
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:C_P86subrend

