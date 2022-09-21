
      SUBROUTINE br89hole(THRESH,Rho,Gam,Tau,Lap,U,D1R,D1G,D1T,D1L)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the analytical BR89 exchange hole value point-wise  *
c    *  the function evaluate the Ux(r) A.1 in CPL paper (see below)  *
c    *  all of mathematical details can be found in the appendix of   *
c    *  this CPL paper                                                *
C    *                                                                *
c    *  reference: A. D. Becke and M. R. Roussel, Phys. Rev. A 39,    *
c    *             3761 (1989);                                       *
c    *   the paper below demonstrate an analytical hole function:     *
c    *             E. Proynov, Z. Gan and J. Kong, "Analytical        *
c    *   representation of the Becke-Roussel exchange functional"     *
c    *             Chem. Phys. Lett. 455, 103 (2008). paper CPL       *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *                                                                *
c    *   output:                                                      *
c    *   U: hole function value, see A.1                              *
c    *   D1R:  first derivative for U with Rho                        *
c    *   D1G:  first derivative for U with gamma                      *
c    *   D1T:  first derivative for U with tau                        *
c    *   D1L:  first derivative for U with Laplacian                  *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)

C
C    parameters defined in paper CPL, appendix
C       
       parameter(a1 = 1.525525181200953D0)
       parameter(a2 = 0.4576575543602858D0)
       parameter(a3 = 0.4292036732051034D0)

       parameter(B  = 2.085749716493756D0)

       parameter(c0 = 0.7566445420735584D0)
       parameter(c1 = -2.636397787137096D0)
       parameter(c2 = 5.474515996423288D0)
       parameter(c3 = -12.65730812710829D0)
       parameter(c4 = 4.125058472512136D0)
       parameter(c5 = -30.42513395716384D0)

       parameter(b0 = 0.4771976183772063D0)
       parameter(b1 = -1.779981349455627D0)
       parameter(b2 = 3.843384186230215D0)
       parameter(b3 = -9.591205088051849D0)
       parameter(b4 = 2.173018028591672D0)
       parameter(b5 = -30.42513385160366D0)
       
       parameter(d0 = 0.00004435009886795587D0)
       parameter(d1 = 0.5812865360445791D0)
       parameter(d2 = 66.74276451594061D0)
       parameter(d3 = 434.2678089722977D0)
       parameter(d4 = 824.7765766052239D0)
       parameter(d5 = 1657.965273158212D0)

       parameter(e0 = 0.00003347285060926091D0)
       parameter(e1 = 0.4791793102397135D0)
       parameter(e2 = 62.39226833857424D0)
       parameter(e3 = 463.1481642793812D0)
       parameter(e4 = 785.2360350104029D0)
       parameter(e5 = 1657.962968223273D0)

c
c      this is the cutoff threshold to determine
c      which branch of g(y) we really go
c      the implementation is sticky to the threshold value
c      provided by emil's old program
c
       parameter(CUT_OFF = 1.0D-10)

C
C      it's important to note that all of implementation here
C      follows the paper of CPL noted above
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F23  = 2.0D0/3.0D0
       F53  = 5.0D0/3.0D0
       F16  = 1.0D0/6.0D0
       F12  = 1.0D0/2.0D0
       F112 = 1.0D0/12.0D0
       F13  = 1.0D0/3.0D0
       PI13 = PI**F13
       PI23 = PI13*PI13

       ! let's build the y(x) in A.2
       Y    = 0.0D0
       DYDR = 0.0D0
       DYDL = 0.0D0
       DYDT = 0.0D0
       DYDG = 0.0D0
       IF (Rho > THRESH) THEN

          ! Q in equation 3 and 4
          ! and it's derivatives
          D    = Tau - 0.25D0*Gam/Rho
          Q    = F16*(Lap-2.0D0*D)
          DQDR = 0.0D0
          IF (Rho*Rho > THRESH) THEN
             DQDR = -F112*Gam/(Rho*Rho)
          END IF
          DQDL = F16
          DQDT = -F13
          DQDG = F112/Rho
C          write(6,*)"Rho",Rho
C          write(6,*)"Gam",Gam
C          write(6,*)"Tau",Tau
C          write(6,*)"Lap",Lap
C          write(6,*)"Q",Q
C          write(6,*)"DQDT",DQDT

          ! y(x) and it's derivatives
          Q2   = Q*Q;
          R23  = Rho**F23
          R53  = R23*Rho
          PREY = F23*PI23
          IF (ABS(Q) > THRESH) THEN
             Y    = PREY*R53/Q
             DYDQ = -Y/Q
             IF (Q2 > THRESH) THEN
                DYDR = PREY*(F53*R23*Q-DQDR*R53)/Q2
             END IF
             DYDL = DYDQ*DQDL
             DYDT = DYDQ*DQDT
             DYDG = DYDQ*DQDG
          END IF

       ENDIF

       ! now let's build x = x(y) function in A.2
       ! we note, that y > 0 is actually for y > threshold
       ! the rest of case will go to y<=0
       X     = 0.0D0
       DXDY  = 0.0D0
       IF (Y > CUT_OFF) THEN

          ! P1(y) and P2(y)
          Y2 = Y*Y
          Y3 = Y2*Y
          Y4 = Y3*Y
          Y5 = Y4*Y
          P1Y = d0 + d1*Y + d2*Y2 + d3*Y3 + d4*Y4 + d5*Y5
          DP1YDY = d1 + 2.0D0*d2*Y + 3.0D0*d3*Y2 + 4.0D0*d4*Y3 + 
     $5.0D0*d5*Y4 
          P2Y = e0 + e1*Y + e2*Y2 + e3*Y3 + e4*Y4 + e5*Y5
          DP2YDY = e1 + 2.0D0*e2*Y + 3.0D0*e3*Y2 + 4.0D0*e4*Y3 + 
     $5.0D0*e5*Y4 

          ! g(y) = arcsch(By) + 2
          BY      = B*Y
          ONEDBY  = 1.0D0/BY
          ONEDBY2 = ONEDBY*ONEDBY
          SQBY2   = DSQRT(1.0D0+ONEDBY2)
          GY      = DLOG(ONEDBY + SQBY2) + 2.0D0
          DGYDY   = -B/(BY*BY*SQBY2)

          ! now let's assemble all of them together
          X       = GY*P1Y/P2Y
          P2YSQ   = P2Y*P2Y
          IF (P2YSQ > THRESH) THEN
             DXDY    = (DGYDY*P1Y*P2Y+GY*DP1YDY*P2Y-GY*P1Y*DP2YDY)/P2YSQ
          END IF

       ELSE IF (Y < -CUT_OFF) THEN

          ! P1(y) and P2(y)
          Y2 = Y*Y
          Y3 = Y2*Y
          Y4 = Y3*Y
          Y5 = Y4*Y
          P1Y = c0 + c1*Y + c2*Y2 + c3*Y3 + c4*Y4 + c5*Y5
          DP1YDY = c1 + 2.0D0*c2*Y + 3.0D0*c3*Y2 + 4.0D0*c4*Y3 + 
     $5.0D0*c5*Y4 
          P2Y = b0 + b1*Y + b2*Y2 + b3*Y3 + b4*Y4 + b5*Y5
          DP2YDY = b1 + 2.0D0*b2*Y + 3.0D0*b3*Y2 + 4.0D0*b4*Y3 + 
     $5.0D0*b5*Y4 
c          write(6,*)"P12Y", P1Y/P2Y

          ! g(y) = -arctan(a1*y+a2) + a3
          A1Y     = a1*Y + a2
          GY      = -DATAN(A1Y) + a3
          DGYDY   = -a1/(1.0D0+A1Y*A1Y)
c          write(6,*)"GY", GY

          ! now let's assemble all of them together
          X       = GY*P1Y/P2Y
          P2YSQ   = P2Y*P2Y
          IF (P2YSQ > THRESH) THEN
             DXDY    = (DGYDY*P1Y*P2Y+GY*DP1YDY*P2Y-GY*P1Y*DP2YDY)/P2YSQ
          END IF
       ENDIF
c       write(6,*)"X", X, "Y", Y

       ! finally, let's do the BR89 exchange hole of U
       ! see the function of A.1
       U1    = 0.5D0  ! this is the limit of U1 as x is zero
       DU1DX = 1.0D0/6.0D0  ! this is the limit of the derivatives as x is zero
       IF (ABS(X) > THRESH) THEN

          ! firstly, do the U1 part, it only contains the x
          T1     = DEXP(-X)
          T2     = DEXP(X/3.0D0)
          T3     = 1.0D0-T1-0.5D0*X*T1
          U1     = T2*T3/X
          DU1DX  = F13*U1 - U1/X + T2*(F12*T1+F12*X*T1)/X
       END IF
C       write(6,*)"U", U1, "DU1DX", DU1DX
C       write(6,*)"DXDY", DXDY, "DYDR", DYDR

       ! now it's real U and it's derivatives
       CON    = -2.0D0*PI13
       R13    = Rho**F13
       U      = CON*R13*U1
       D1G    = CON*R13*DU1DX*DXDY*DYDG
       D1T    = CON*R13*DU1DX*DXDY*DYDT
       D1L    = CON*R13*DU1DX*DXDY*DYDL
C          write(6,*) "fenglai U12+U3", U1
C       write(6,*) "DYDT", "D1T", DYDT, D1T

       ! derivatives for rho, we need to see 
       D1R    = 0.0D0
       IF (Rho > THRESH) THEN
          MR23   = Rho**(-2.0D0/3.0D0)
          D1R    = CON*(F13*U1*MR23+R13*DU1DX*DXDY*DYDR)
       END IF
       
       RETURN 
       END
