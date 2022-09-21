C
C note by fenglai on Dec. 2015:
C 
C we have done a tough finite difference in comparing the first 
C order derivatives for the functionals. For B05, right now
C we have two issues:
C 1   in calculating Y, then X(Y) function; we basically omit
C     all of grids with density below the CUT_OFF3 (see the 
C     function of becke05_xy). For these grids Y is 0, X is 2
C     then accordingly N is 0, F factor is 0 etc. For these 
C     grids the non-dynamic contribution is zero. 
C
C 2   In calculating F factor, the large P value(115) brings 
C     trouble in comparing find difference with the analytical
C     ones. It's only when you decrease the P value (I tested 0.01)
C     then the F factor calculation can leading to reasonable 
C     match between find difference and analytical ones.
C     we believe, that large P could be a main reason leading 
C     SCF hard to converge in the open shell case
C
C changes made on Nov. 28 2016 for preparing new par_opt:
C
C 1   we do need to use cut_off3. we tried to use the thresh
C     input by the program instead, however it seems that the SCF
C     calculation becomes unstable. the reason for this is not known.
C     Originally the old way will bypass
C     a lot of grid points(because the original T=R^2*4*PI and
C     we compare T with 1.0E-8, so nearly all of R>1.0E-4 will
C     have X=2, Y=0, NEFF=0 so you will see that the functional
C     value and derivatives for these points are all zero accordingly).
C     it seems that for SCF stability such feature is needed to keep 
C     anyway. However, how unstable the SCF becomes we did not get
C     too much test. This is still a question hangling over there.
C     We tested Li atom(G3Large with svwn5 initial guess), the SCF 
C     calculation with tight threshold
C     value gets positive energy, however it's not going wild. 
C     therefore it maybe possible that we can improve it with 
C     advanced SCF algorithms like GDM.
C
C 2   now instead of using emil's three branch function, we use 
C     the newtion-raphon procedure to derive the dxy. This is the
C     more accurate way. Emil's original region III to derive the X value
C     through Y is not accurate. This can correct it and imrpove the precision.
C
C     further note: the region III using formula proposed by other groups.
C
C   
C   

      SUBROUTINE BR89_CURVATURE(THRESH,Rho,Gam,Tau,Lap,HIRWT,Q,
     $DQDR,DQDG,DQDT,DQDL)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the curvature for BR89 exchange hole model          *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   HIRWT: Hirshfeld weights using for B13 Strong functional     *
c    *          for the normal B05 functional etc. the input          *         
c    *          Hirshfeld weights should be zero every grid point     *         
c    *          For B13 strong functional it's value determined       *         
c    *          outside                                               *
C    *                                                                *
C    *   also it's worthy to note that hirshfeld weights only depends *
C    *   on the atomic density, therefore it's constant during SCF    *
C    *   process. We do not need to do derivatives for hirshfeld      *
C    *   weights.                                                     *
c    *                                                                *
c    *   output:                                                      *
c    *   Q   :  the curvature value                                   *
c    *   DQDR:  first derivative for Q with Rho                       *
c    *   DQDG:  first derivative for Q with gamma                     *
c    *   DQDT:  first derivative for Q with tau                       *
c    *   DQDL:  first derivative for Q with Laplacian                 *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F16  = 1.0D0/6.0D0
       F112 = 1.0D0/12.0D0
       F13  = 1.0D0/3.0D0
       F14  = 1.0D0/4.0D0
       F53  = 5.0D0/3.0D0
       F23  = 2.0D0/3.0D0
       ZERO = 0.0D0
       TWO  = 2.0D0

       ! build Q
       Q    = ZERO
       DQDR = ZERO
       DQDL = ZERO
       DQDT = ZERO
       DQDG = ZERO
       IF (Rho > THRESH) THEN
          D    = Tau - F14*Gam/Rho
          R23  = Rho**F23
          R53  = R23*Rho
          Q    = F16*(Lap-TWO*D) + HIRWT*R53
          DQDR = HIRWT*F53*R23
          IF (Rho*Rho > THRESH) THEN
             DQDR = DQDR -F112*Gam/(Rho*Rho) 
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
       END IF
       RETURN
       END

      SUBROUTINE NR_Becke05XY(THRESH,Y,X)
c    ******************************************************************
c    *                                                                *
c    *  evaluates the X(y) function for Becke05 exchange functional   *
c    *  through newton-raphson method                                 *
c    *                                                                *
c    *  this is for a given y value, we calculate it's corresponding  *
c    *  x value. the fomular is defined at equation 18 of paper below *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT REAL*8(A-Z)
      INTEGER  I
      INTEGER  LOOPS
      PARAMETER(LOOPS=1000)
      REAL*8 CRITERIA
      PARAMETER(CRITERIA=1.0D-10)

      ! constants
      ZERO = 0.0D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      THREE= 3.0D0
      FOUR = 4.0D0

      DO I = 1, LOOPS

         ! calculate ex
         ex = DEXP(X)
c         write(6,*)"in NR, cycle", I
c         write(6,*)"in NR, Y value", y
c         write(6,*)"in NR, initial X value", x

         ! compute fx
         xy = ZERO
         if (DABS(X*X)>THRESH) THEN
            xy = ((X-TWO)/(X*X))*(ex-ONE-X/TWO) - Y
         END IF

         ! compute the f'(x)
         dxy = ZERO
         if (DABS(X)>THRESH) THEN
            dxy = ex/X
         END IF
         if (DABS(X*X)>THRESH) THEN
            X2  = X*X
            dxy = dxy - THREE*ex/X2
         END IF
         if (DABS(X*X*X)>THRESH) THEN
            X3  = X*X*X
            dxy = dxy + FOUR*(ex-ONE)/X3
         END IF
c         write(6,*)"in NR, y(x)", xy
c         write(6,*)"in NR, dy(x)/dx", dxy

         ! let's see whether the difference between
         ! new X and old X is small enough
         diff = ONE
         IF (DABS(dxy)>THRESH) THEN
            diff = DABS(xy/dxy) 
         END IF
         IF (diff<CRITERIA) THEN
            return
         END IF
c         write(6,*)"in NR, diff", diff

         ! now derive the new X value
         NEW_X = X
         IF (DABS(dxy)>THRESH) THEN
            NEW_X = X - xy/dxy 
         END IF
c         write(6,*)"in NR, old x", x
c         write(6,*)"in NR, new x", new_x
         X = NEW_X
      END DO

      IF (I .ge. LOOPS) THEN
         stop "something wrong with the NR for solving X(y)"
      END IF

      RETURN 
      END

      SUBROUTINE becke05_xy(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,
     $X,D1R,D1G,D1T,D1L,D1U)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the X(y) function for Becke05 exchange functional   *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the x(y) functional is defined and analyzed in section III    *
c    *  of the above paper                                            *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   Ux:  the exact exchange energy density                       *
c    *                                                                *
c    *   output:                                                      *
c    *   X: the value of X(Y)                                         *
c    *   D1R:  first derivative for X with Rho                        *
c    *   D1G:  first derivative for X with gamma                      *
c    *   D1T:  first derivative for X with tau                        *
c    *   D1L:  first derivative for X with Laplacian                  *
c    *   D1U:  first derivative for X with exchange energy density    *         
c    *                                                                *
c    *                                                                *         
c    *  we note that Ux is not directly used. In constructing         *
c    *  the potential we actually use it's potential, which is        *
c    *  divided by rho                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)

C
C    parameters defined in paper J. Chem. Phys. 122, 64101 (2005), appendix
C       
C      this is for region I
       parameter(a1 =  0.9301841768238374D0 )
       parameter(a2 =  0.5485569431916153D0 )
c
       parameter(c0 = -5.9685280129072020D0 )
       parameter(c1 = -2.1837477426038480D0 )
       parameter(c2 = -4.9858864412437560D0 )
       parameter(c3 = -1.1341612120636830D0 )
       parameter(c4 = -1.6921426426199750D0 )
       parameter(c5 =  0.57089595383468940D0)
c
       parameter(b0 = -5.9685280130660880D0 )
       parameter(b1 = -2.030780232084790D0  )
       parameter(b2 = -4.6796750480012860D0 )
       parameter(b3 = -1.1188490577541170D0 )
       parameter(b4 = -1.8087055034029230D0 )
       parameter(b5 =  0.59226482161752910D0)
c
       parameter(q0 =   0.9129908719446399D0  )
       parameter(q1 =   3.655262558462426D0   )
       parameter(q2 =   0.1801828684494572D0  )
       parameter(q3 =  -3.062938667772561D0   )
       parameter(q4 =  -1.173405187745653D0   )
       parameter(q5 =  -1.662674088158794D0   )
       parameter(q6 =   0.6859613559654089D0  )
       parameter(q7 =   0.06595477584967326D0 )
       parameter(q8 =  -0.03038609318852905D0 )
       parameter(q9 =  -0.00000000000000077D0 )
c
c     THIS IS for region II
       parameter(d0 = 26.88413661379433D0  )
       parameter(d1 = 46.96693640414017D0  )
       parameter(d2 = 33.17800151829805D0  )
       parameter(d3 = 9.940088877152307D0  )
       parameter(d4 = 0.8786661786414733D0  )
       parameter(d5 = 0.01643722176146135D0  )
       parameter(d6 = 3.40424464772731D-5  )
c
       parameter(e0 = 3.360517076724290D0     )
       parameter(e1 = 4.623703278485152D0     )
       parameter(e2 = 2.688949840405010D0     )
       parameter(e3 = 0.6007166968496472D0    )
       parameter(e4 = 0.03922040006408070D0   )
       parameter(e5 = 0.0005438465669613952D0 )
       parameter(e6 = 0.00000078437439010087D0)
c
       parameter(k0 =   8.382230306773296D0         )
       parameter(k1 =   19.60257290836815D0         )
       parameter(k2 =   19.71894106502952D0         )
       parameter(k3 =   10.77146542063455D0         )
       parameter(k4 =   3.415698370189622D0         )
       parameter(k5 =   0.5813153924917321D0        )
       parameter(k6 =   0.05426084060061605D0       )
       parameter(k7 =   0.002299629631716270D0      )
       parameter(k8 =   0.00005119354330427682D0    )
       parameter(k9 =   0.000000322977561012273D0   )
       parameter(k10 =  0.000000001405232963383258D0)
C
C
C      this is for region III       
       parameter(l1 =  1.23767D0)
       parameter(l2 =  9.37587D0)
       parameter(l3 = 19.4777D0)
       parameter(l4 = 13.6816D0)
       parameter(l5 =  0.078655D0)
       parameter(l6 = 54.7264D0)
       parameter(l7 = 58.4331D0)
       parameter(l8 = 18.75174D0)
       parameter(l9 =  1.23767D0)

c
c      this is the cutoff threshold to determine
c      which branch of x(y) we really go
c
       parameter(CUT_OFF1 = 1.0D-13)
       parameter(CUT_OFF2 = 1000.0D0)
c
c
c      this is the cut off value emil suggest used 
c      in keeping the SCF stable       
       parameter(CUT_OFF3 = 1.0D-8)
c
c
c      this is the constant of e
c
c
       parameter(e = 2.7182818284590452353602875D0)
C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F23  = 2.0D0/3.0D0
       F34  = 3.0D0/4.0D0
       F53  = 5.0D0/3.0D0
       F16  = 1.0D0/6.0D0
       F12  = 1.0D0/2.0D0
       F112 = 1.0D0/12.0D0
       F13  = 1.0D0/3.0D0
       F14  = 1.0D0/4.0D0
       F83  = 8.0D0/3.0D0
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0
       FIVE = 5.0D0
       C24  = 24.0D0

       ! build Q
       Q    = ZERO
       DQDR = ZERO
       DQDG = ZERO
       DQDT = ZERO
       DQDL = ZERO
       CALL BR89_CURVATURE(THRESH,Rho,Gam,Tau,Lap,HIRWT,Q,
     $ DQDR,DQDG,DQDT,DQDL)
c       write(6,*)"our Q in building Y:", Q

       ! build Y
       ! here the use of cut_off3 is suggested
       ! by emil
       !
       ! we found that if we use smaller cut off value 
       ! than the cut_off3, the Y value and the corresponding 
       ! X value may result in derivatives calculation
       ! inaccurate. Therefore right now we use a loose
       ! cut off value so that to keep SCF stable.
       !
       ! basically, if less than the cut_off3 Y is zero.
       ! then X is 2. In this case the Neff is 0, too.
       !
       Y    = ZERO
       Y0   = ZERO
       DYDR = ZERO
       DYDG = ZERO
       DYDT = ZERO
       DYDL = ZERO
       DYDU = ZERO
       VX   = UX/Rho
       R2   = Rho*Rho
       C    = F34/PI
       T    = Rho*Rho*FOUR*PI
       IF (T > CUT_OFF3) THEN
          Y0   = -(C/R2)*VX
          Y    = Y0*Q
          DVDR = -UX/R2
          DYDR = -C*((DQDR*VX+Q*DVDR)/R2)
          DYDR = DYDR+C*(TWO*(Q/Rho)*(VX/R2))
          DYDU = -C*(Q/Rho)*(ONE/R2)
       END IF
       DYDG = Y0*DQDG
       DYDT = Y0*DQDT
       DYDL = Y0*DQDL
c       write(6,*)"Y: ", Y
c       write(6,*)"DYDR: ", DYDR
c       write(6,*)"DYDG: ", DYDG
c       write(6,*)"DYDT: ", DYDT
c       write(6,*)"DYDL: ", DYDL
c       write(6,*)"DYDU: ", DYDU

       IF (DABS(Y)<THRESH) THEN
          X    = TWO
          EX   = DEXP(X)
          DYDX = EX/X - THREE*EX/(X*X) + FOUR*(EX-ONE)/(X*X*X)
          DXDY = ONE/DYDX
          D1R = DXDY*DYDR
          D1G = DXDY*DYDG
          D1T = DXDY*DYDT
          D1L = DXDY*DYDL
          D1U = DXDY*DYDU
          RETURN 
       END IF

       ! now let's build x = x(y) 
       X     = ZERO
       DXDY  = ZERO
       IF (Y < -CUT_OFF1) THEN

          ! region I
          ! P1(y) and P2(y)
          Y2 = Y*Y
          Y3 = Y2*Y
          Y4 = Y3*Y
          Y5 = Y4*Y
          P1Y = c0 + c1*Y + c2*Y2 + c3*Y3 + c4*Y4 + c5*Y5
          P2Y = b0 + b1*Y + b2*Y2 + b3*Y3 + b4*Y4 + b5*Y5

          ! g(Y)
          T1  = TWO*DATAN(a2)+PI
          g_Y = TWO*(TWO*DATAN(a1*Y+a2)+PI)/T1

          ! now do diravative for g(y)
          DGY = FOUR*a1/((ONE+(a1*Y+a2)*(a1*Y+a2))*T1)

          ! dirivatives for P1(y)/P2(y)
          Y6 = Y5*Y
          Y7 = Y6*Y
          Y8 = Y7*Y
          Y9 = Y8*Y
          QY = q0 + q1*Y + q2*Y2 + q3*Y3 + q4*Y4 + q5*Y5 + q6*Y6 + q7*Y7
     $+ q8*Y8 + q9*Y9 
          DPYDY = ZERO
          IF (P2Y*P2Y>THRESH) THEN
             DPYDY =  QY/(P2Y*P2Y)
          END IF

          ! now it's x(y) and it's derivatives
          X    = g_Y*P1Y/P2Y
          DXDY = DGY*P1Y/P2Y + g_Y*DPYDY
c          write(6,*)"switch 1"
c          write(6,*)"X: ", X
c          write(6,*)"DXDY: ", DXDY

       ELSE IF (Y >= -CUT_OFF1 .and. Y <= CUT_OFF1 ) THEN

          ! this is asymptotic limit of region I defined 
          ! in equation 27 of the 2012 paper
          X    = F83 - e*e/THREE + 
     $ F13*DSQRT(FOUR-FOUR*e*e+e**4 + C24*Y)
          DXDY = FOUR/DSQRT(FOUR-FOUR*e*e+e**4 + C24*Y)

c          write(6,*)"switch 2"
c          write(6,*)"X: ", X
c          write(6,*)"DXDY: ", DXDY

       ELSE IF (Y > CUT_OFF1 .and. Y<CUT_OFF2) THEN

          ! region II
          ! P1(y) and P2(y)
          Y2 = Y*Y
          Y3 = Y2*Y
          Y4 = Y3*Y
          Y5 = Y4*Y
          Y6 = Y5*Y
          Y7 = Y6*Y
          Y8 = Y7*Y
          Y9 = Y8*Y
          Y10= Y9*Y
          R1Y = d0 + d1*Y + d2*Y2 + d3*Y3 + d4*Y4 + d5*Y5 + d6*Y6
          R2Y = e0 + e1*Y + e2*Y2 + e3*Y3 + e4*Y4 + e5*Y5 + e6*Y6
          QY  = k0 + k1*Y + k2*Y2 + k3*Y3 + k4*Y4 + k5*Y5 + k6*Y6
     $+ k7*Y7 + k8*Y8 + k9*Y9 + k10*Y10
          IF (DABS(R2Y) > THRESH) THEN
             X    = R1Y/(FOUR*R2Y)
             IF (DABS(R2Y*R2Y) > THRESH) THEN
                DXDY = QY/(R2Y*R2Y) 
             END IF
          END IF

c          write(6,*)"switch 3"
c          write(6,*)"X: ", X
c          write(6,*)"DXDY: ", DXDY

       ELSE IF (Y>CUT_OFF2) THEN

          ! region III
          LNY  = DLOG(Y)
          LNY1 = ONE/LNY
          LNY2 = LNY1/LNY
          LNY3 = LNY2/LNY
          LNY4 = LNY3/LNY
          T1   = LNY**(ONE+LNY1)
          X    = T1 + l1*LNY1 + l2*LNY2 - l3*LNY3 + l4*LNY4 - l5

          ! now dirivatives
          T2   = LNY**(LNY1)
          T3   = DLOG(LNY)
          T4   = ONE/(Y*LNY*LNY)
          DXDY = (T2/Y)*(-T3/LNY+ONE+LNY1) 
     $ + T4*(-l6*LNY3+l7*LNY2-l8*LNY1-l9)

c          write(6,*)"switch 4"
c          write(6,*)"X: ", X
c          write(6,*)"DXDY: ", DXDY

       ENDIF
C       write(6,*)"X", X, "Y", Y, "DXDY", DXDY

       ! finally, bring every thing for derivatives of x(y)
C       D1R = DXDY*DYDR
C       D1G = DXDY*DYDG
C       D1T = DXDY*DYDT
C       D1L = DXDY*DYDL
C       D1U = DXDY*DYDU

       ! this is the new way to calculate the DXDY
       ! comparing with emil's way, which has precision of 1.0D-4
       ! this way has precision of 1.0D-5
       xy = ZERO
       CALL NR_Becke05XY(THRESH,Y,X)
       IF (DABS(X*X)>THRESH) THEN
          xy = ((X-TWO)/(X*X))*(DEXP(X)-ONE-X/TWO)
       END IF
       dxy = ZERO
       If (DABS(X)>THRESH) THEN
          dxy = DEXP(X)/X
       END IF
       If (DABS(X*X)>THRESH) THEN
          X2  = X*X
          dxy = dxy - THREE*DEXP(X)/X2
       END IF
       If (DABS(X*X*X)>THRESH) THEN
          X3  = X*X*X
          dxy = dxy + FOUR*(DEXP(X)-ONE)/X3
       END IF
       If (DABS(dxy)>THRESH) THEN
          DXDY = ONE/dxy
       END IF
       D1R = DXDY*DYDR
       D1G = DXDY*DYDG
       D1T = DXDY*DYDT
       D1L = DXDY*DYDL
       D1U = DXDY*DYDU

       ! this is debug code for the x(y)
       ! we compare the xy with y
       ! also compare the calculated DXDY with emil's DXDY
c       if (DABS(xy-Y)> 0.000000001D0) THEN
c          write(6,*) "x(y) is different from Y", xy, Y
c          write(6,*) "diff is", DABS(xy-Y)
c       END IF
c       if (DABS(dxy)> 0.00000001D0 .and. 
c     $ DABS(ONE/dxy-DXDY)> 0.0001D0) THEN
c          write(6,*) "derivatives of x(y) is different from DXDY", 
c     $ ONE/dxy, DXDY 
c          write(6,*) "diff is", DABS(ONE/dxy-DXDY)
c       END IF

       ! do find difference for in terms of x(y)
c       Y1 = Y + Y*0.01
c       X1 = X
c       CALL NR_Becke05XY(THRESH,Y1,X1)
c       Y2 = Y - Y*0.01
c       X2 = X
c       CALL NR_Becke05XY(THRESH,Y2,X2)
c       IF (DABS(Y)>0.00000001D0) THEN
c          NEW_DXDY = (X1-X2)/(Y1-Y2)
c          IF (DABS(DXDY-NEW_DXDY)>0.0001D0) THEN
c             write(6,*)"find difference on the x(y) does not match"
c             write(6,*)"difference is", DABS(DXDY-NEW_DXDY)
c          END IF
c       END IF

       RETURN 
       END

      SUBROUTINE becke05_neff(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,N,
     $ D1R,D1G,D1T,D1L,D1U)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the N_{eff} for Becke05 exchange functional         *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   Ux:  the exact exchange energy density                       *
c    *                                                                *
c    *   output:                                                      *
c    *   N  :  N_{eff}                                                *
c    *   D1R:  first derivative for N with Rho                        *
c    *   D1G:  first derivative for N with gamma                      *
c    *   D1T:  first derivative for N with tau                        *
c    *   D1L:  first derivative for N with Laplacian                  *
c    *   D1U:  first derivative for N with exchange energy density    *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)
C
C      here the smooth constant defined for N_{eff}
C       
       parameter(DELTA     = 0.07D0)
       parameter(N_LIMIT   = 2.00D0)
C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F23  = 2.0D0/3.0D0
       F32  = 3.0D0/2.0D0
       F16  = 1.0D0/6.0D0
       F112 = 1.0D0/12.0D0
       F13  = 1.0D0/3.0D0
       F12  = 1.0D0/2.0D0
       F52  = 5.0D0/2.0D0
       F14  = 1.0D0/4.0D0
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0

       ! build Q
       Q    = ZERO
       DQDR = ZERO
       DQDG = ZERO
       DQDT = ZERO
       DQDL = ZERO
       CALL BR89_CURVATURE(THRESH,Rho,Gam,Tau,Lap,HIRWT,Q,
     $ DQDR,DQDG,DQDT,DQDL)
c       write(6,*)"Q: ", Q
c       write(6,*)"DQDR: ", DQDR
c       write(6,*)"DQDG: ", DQDG
c       write(6,*)"DQDT: ", DQDT
c       write(6,*)"DQDL: ", DQDL

       ! build X
       X    = ZERO
       DXDR = ZERO
       DXDG = ZERO
       DXDT = ZERO
       DXDL = ZERO
       DXDU = ZERO
       CALL becke05_xy(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,X,
     $ DXDR,DXDG,DXDT,DXDL,DXDU)
c       write(6,*)"X: ", X
c       write(6,*)"DXDR: ", DXDR
c       write(6,*)"DXDG: ", DXDG
c       write(6,*)"DXDT: ", DXDT
c       write(6,*)"DXDL: ", DXDL

       ! now build f(x), which reply on x and Q
       ! we note that, (X-TWO)/(X*Q) should be 
       ! put together
       ! this is because if X is smaller than 2,
       ! Q is negative therefore overall this value
       ! is always positive
       ! also, when Q is zero then X is two
       ! also when x is zero Q is some large negtive value
       ! therefore QX is never zero
       ! X itself should be always >=0
       FX   = ZERO
       DFDR = ZERO
       DFDG = ZERO
       DFDT = ZERO
       DFDL = ZERO
       DFDU = ZERO
       IF (DABS(Q) >THRESH .and. DABS((X-TWO)/(X*Q))>THRESH) THEN
          QX    = Q*X
          QX2   = QX*QX
          T     = (X-TWO)/(X*Q)
          T12   = T**F12  
          T32   = T12*T
          DTDR  = F32*T12*(DXDR/QX-(X-TWO)*(DQDR*X+DXDR*Q)/QX2)
          DTDG  = F32*T12*(DXDG/QX-(X-TWO)*(DQDG*X+DXDG*Q)/QX2)
          DTDT  = F32*T12*(DXDT/QX-(X-TWO)*(DQDT*X+DXDT*Q)/QX2)
          DTDL  = F32*T12*(DXDL/QX-(X-TWO)*(DQDL*X+DXDL*Q)/QX2)
          DTDU  = F32*T12*(DXDU/QX-(X-TWO)*(DXDU*Q)/QX2)

          CON   = F23**F32*PI
          S     = CON*DEXP(X)
          FX    = S*T32
          DFDR  = FX*DXDR + S*DTDR
          DFDG  = FX*DXDG + S*DTDG
          DFDT  = FX*DXDT + S*DTDT
          DFDL  = FX*DXDL + S*DTDL
          DFDU  = FX*DXDU + S*DTDU
       END IF

       ! now let's derive the Neff
       ! this is based on the original definition of N_{eff}
       ! in equation of 20 of 2012 JCP paper
       R12  = Rho**F12
       R32  = Rho*R12
       R52  = R32*Rho
       N1   = FX*R52
       TD1R = DFDR*R52+F52*FX*R32
       TD1G = DFDG*R52
       TD1T = DFDT*R52
       TD1L = DFDL*R52
       TD1U = DFDU*R52

       ! now we need to do a smooth for the N_{eff}
       ! for the SCF
       N    = ZERO
       D1R  = ZERO
       D1G  = ZERO
       D1T  = ZERO
       D1L  = ZERO
       D1U  = ZERO
       IF (N1 <= N_LIMIT-DELTA) THEN
          N    = N1
          D1R  = TD1R
          D1G  = TD1G
          D1T  = TD1T
          D1L  = TD1L
          D1U  = TD1U
       ELSE IF (N1 > N_LIMIT-DELTA .and. N1 < N_LIMIT+DELTA) THEN
          G0   = N1-N_LIMIT-DELTA
          N    = N_LIMIT-G0*G0/(FOUR*DELTA)
          DNDG = -G0/(TWO*DELTA)
          D1R  = TD1R*DNDG
          D1G  = TD1G*DNDG
          D1T  = TD1T*DNDG
          D1L  = TD1L*DNDG
          D1U  = TD1U*DNDG
       ELSE 
          N    = N_LIMIT+DELTA
          D1R  = ZERO 
          D1G  = ZERO 
          D1T  = ZERO 
          D1L  = ZERO 
          D1U  = ZERO 
       END IF
c       write(6,*)"N: ", N
c       write(6,*)"N1R: ", D1R
c       write(6,*)"N1G: ", D1G
c       write(6,*)"N1T: ", D1T
c       write(6,*)"N1L: ", D1L
c       write(6,*)"N1U: ", D1U

       RETURN 
       END


      SUBROUTINE becke05_M1(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,M1,
     $ D1R,D1G,D1T,D1L,D1U)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the M1 for Becke05 exchange functional              *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   Ux:  the exact exchange energy density                       *
c    *                                                                *
c    *   output:                                                      *
c    *   M1 :  equation 21 of 2012 paper                              *
c    *   D1R:  first derivative for M1 with Rho                       *
c    *   D1G:  first derivative for M1 with gamma                     *
c    *   D1T:  first derivative for M1 with tau                       *
c    *   D1L:  first derivative for M1 with Laplacian                 *
c    *   D1U:  first derivative for M1 with exchange energy density   *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)
C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F12  = 1.0D0/2.0D0
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       FOUR = 4.0D0
       SIX  = 6.0D0

       ! build Q
       Q    = ZERO
       DQDR = ZERO
       DQDG = ZERO
       DQDT = ZERO
       DQDL = ZERO
       CALL BR89_CURVATURE(THRESH,Rho,Gam,Tau,Lap,HIRWT,Q,
     $ DQDR,DQDG,DQDT,DQDL)

       ! build X
       X    = ZERO
       DXDR = ZERO
       DXDG = ZERO
       DXDT = ZERO
       DXDL = ZERO
       DXDU = ZERO
       CALL becke05_xy(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,X,
     $ DXDR,DXDG,DXDT,DXDL,DXDU)

       ! build Neff
       N    = ZERO
       DNDR = ZERO
       DNDG = ZERO
       DNDT = ZERO
       DNDL = ZERO
       DNDU = ZERO
       CALL becke05_neff(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,N,
     $ DNDR,DNDG,DNDT,DNDL,DNDU)

       ! now build f(x), which reply on x and Q
       ! we note that, (X-TWO)/(X*Q) should be 
       ! put together
       ! we need to be careful about this form
       ! some cases that X and Q are zero
       FX   = ZERO
       DFDR = ZERO
       DFDG = ZERO
       DFDT = ZERO
       DFDL = ZERO
       DFDU = ZERO
       IF (DABS(X*Q)>THRESH) THEN
          IF (DABS((X-TWO)/(X*Q))>THRESH) THEN

             ! get the F1 part, which only have X
             EXPX = DEXP(-X)
             F1   = ZERO
             DF1DX=ZERO
             IF (DABS(X)>THRESH) THEN
                X2   = X*X
                F1   = X + FOUR/X - (ONE+FOUR/X)*EXPX
                DF1DX= ONE-FOUR/X2+(FOUR/X2)*EXPX+(ONE+FOUR/X)*EXPX
             END IF

             ! now let's see the F2 part, depending on X and Q
             QX   = X*Q
             QX2  = QX*QX
             T    = ZERO
             IF (((X-TWO)/QX)>THRESH) THEN
                T = ((X-TWO)/QX)**F12
             END IF
             F2   = ZERO
             DF2DR= ZERO
             DF2DG= ZERO
             DF2DT= ZERO
             DF2DL= ZERO
             DF2DU= ZERO
             IF (DABS(T)>THRESH) THEN
                TM   = ONE/T
                F2   = T
                DF2DR= F12*TM*(DXDR/QX-(X-TWO)*(DQDR*X+DXDR*Q)/QX2)
                DF2DG= F12*TM*(DXDG/QX-(X-TWO)*(DQDG*X+DXDG*Q)/QX2)
                DF2DT= F12*TM*(DXDT/QX-(X-TWO)*(DQDT*X+DXDT*Q)/QX2)
                DF2DL= F12*TM*(DXDL/QX-(X-TWO)*(DQDL*X+DXDL*Q)/QX2)
                DF2DU= F12*TM*(DXDU/QX-(X-TWO)*(DXDU*Q)/QX2)
             END IF

             ! now bring all of stuff together
             C    = SIX**(-F12)
             FX   = C*F1*F2
             DFDR = C*DF1DX*DXDR*F2+C*F1*DF2DR
             DFDG = C*DF1DX*DXDG*F2+C*F1*DF2DG
             DFDT = C*DF1DX*DXDT*F2+C*F1*DF2DT
             DFDL = C*DF1DX*DXDL*F2+C*F1*DF2DL
             DFDU = C*DF1DX*DXDU*F2+C*F1*DF2DU
          END IF
       END IF

       ! now let's derive the M1
       M1   = ZERO 
       D1R  = ZERO 
       D1G  = ZERO 
       D1T  = ZERO 
       D1L  = ZERO 
       D1U  = ZERO 
       IF (Rho>THRESH) THEN
          R12  = Rho**F12
          M1   = N*FX*R12
          D1R  = DFDR*N*R12 + DNDR*FX*R12 
          D1G  = DFDG*N*R12 + DNDG*FX*R12 
          D1T  = DFDT*N*R12 + DNDT*FX*R12 
          D1L  = DFDL*N*R12 + DNDL*FX*R12 
          D1U  = DFDU*N*R12 + DNDU*FX*R12 
          IF (R12>THRESH) THEN
             D1R  = D1R + F12*N*FX/R12
          END IF
       END IF
c       write(6,*)"M1: ", M1
c       write(6,*)"DM1DR: ", D1R
c       write(6,*)"DM1DG: ", D1G
c       write(6,*)"DM1DT: ", D1T
c       write(6,*)"DM1DL: ", D1L
c       write(6,*)"DM1DU: ", D1U

       RETURN 
       END


      SUBROUTINE becke05_M2(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,M2,
     $ D1R,D1G,D1T,D1L,D1U)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the M2 for Becke05 exchange functional              *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   Ux:  the exact exchange energy density                       *
c    *                                                                *
c    *   output:                                                      *
c    *   M2 :  equation 22 of 2012 paper                              *
c    *   D1R:  first derivative for M2 with Rho                       *
c    *   D1G:  first derivative for M2 with gamma                     *
c    *   D1T:  first derivative for M2 with tau                       *
c    *   D1L:  first derivative for M2 with Laplacian                 *
c    *   D1U:  first derivative for M2 with exchange energy density   *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)
C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! constants
       PI   = 4.D0*DATAN(1.D0)
       F16  = 1.0D0/6.0D0
       F13  = 1.0D0/3.0D0
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       FOUR = 4.0D0
       C12  = 12.0D0

       ! build Q
       Q    = ZERO
       DQDR = ZERO
       DQDG = ZERO
       DQDT = ZERO
       DQDL = ZERO
       CALL BR89_CURVATURE(THRESH,Rho,Gam,Tau,Lap,HIRWT,Q,
     $ DQDR,DQDG,DQDT,DQDL)

       ! build X
       X    = ZERO
       DXDR = ZERO
       DXDG = ZERO
       DXDT = ZERO
       DXDL = ZERO
       DXDU = ZERO
       CALL becke05_xy(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,X,
     $ DXDR,DXDG,DXDT,DXDL,DXDU)

       ! build Neff
       N    = ZERO
       DNDR = ZERO
       DNDG = ZERO
       DNDT = ZERO
       DNDL = ZERO
       DNDU = ZERO
       CALL becke05_neff(THRESH,Rho,Gam,Tau,Lap,UX,HIRWT,N,
     $ DNDR,DNDG,DNDT,DNDL,DNDU)

       ! now build f(x), which only reply on x
       ! different from M1 and N_{eff}, here
       ! it's safe to seperate the (X-2)/X*Q
       ! into (X-2)/X and 1/Q two parts
       FX   = ZERO
       DFDX = ZERO
       IF (DABS(X)>THRESH) THEN
          X2   = X*X
          FX   = F16*(X-TWO)/X*(X2+C12)
          DFDX = F13*X-F13+FOUR/X2
       END IF

       ! now let's derive the M2
       M2   = ZERO
       D1R  = ZERO
       D1G  = ZERO
       D1T  = ZERO
       D1L  = ZERO
       D1U  = ZERO
       Q2   = Q*Q
       if (DABS(Q)>THRESH) THEN
          M2   = N*FX*Rho/Q
          D1R  = DFDX*DXDR*N*Rho/Q + DNDR*FX*Rho/Q + N*FX/Q
          D1G  = DFDX*DXDG*N*Rho/Q + DNDG*FX*Rho/Q 
          D1T  = DFDX*DXDT*N*Rho/Q + DNDT*FX*Rho/Q 
          D1L  = DFDX*DXDL*N*Rho/Q + DNDL*FX*Rho/Q 
          D1U  = DFDX*DXDU*N*Rho/Q + DNDU*FX*Rho/Q 
          IF (Q2>THRESH) THEN
             D1R  = D1R - N*FX*Rho*DQDR/Q2
             D1G  = D1G - N*FX*Rho*DQDG/Q2
             D1T  = D1T - N*FX*Rho*DQDT/Q2
             D1L  = D1L - N*FX*Rho*DQDL/Q2
          END IF
       END IF
C       write(6,*)"M2: ", M2
C       write(6,*)"DM2DR: ", D1R
C       write(6,*)"DM2DG: ", D1G
C       write(6,*)"DM2DT: ", D1T
C       write(6,*)"DM2DL: ", D1L
C       write(6,*)"DM2DU: ", D1U

       RETURN 
       END

      SUBROUTINE becke05_f(THRESH,VAL_P,HIRWT,Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,F,D1R1,D1G1,D1T1,D1L1,D1U1,
     $ D1R2,D1G2,D1T2,D1L2,D1U2)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the f factor for the exchange of Becke05 functional *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input(1 is for alpha spin components, 2 is for beta)         *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   UX:  the exact exchange energy density                       *
c    *                                                                *
c    *   output(1 is for alpha spin components, 2 is for beta)        *
c    *   F  :  this is the f factor defined in equation 13 and 46     *
c    *   D1R:  first derivative for F with Rho                        *
c    *   D1G:  first derivative for F with gamma                      *
c    *   D1T:  first derivative for F with tau                        *
c    *   D1L:  first derivative for F with Laplacian                  *
c    *   D1U:  first derivative for F with exchange energy density    *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    *   if you change the value of F limit below, be sure to check   *
c    *   all of formulas coded here! it may not work if you only      *
c    *   change the F_LIMIT value.                                    *         
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)

C
C      here it's related to the smooth technique of f factor
C       
       parameter(F_LIMIT = 1.0D0  )
       parameter(DELTA   = 0.05D0 )

C
C      the original parameter of P is replaced with the variable of P
C       
C       parameter(P       = 115.0D0)

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0
       FIVE = 5.0D0
       F12  = 0.50d0 
       P    = VAL_P

C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! build Neff for alpha
       N1    = ZERO
       DN1DR = ZERO
       DN1DG = ZERO
       DN1DT = ZERO
       DN1DL = ZERO
       DN1DU = ZERO
       CALL becke05_neff(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,N1,
     $ DN1DR,DN1DG,DN1DT,DN1DL,DN1DU)
c       write(6,*)"alpha neff: ", N1

       ! build Neff for beta
       N2    = ZERO
       DN2DR = ZERO
       DN2DG = ZERO
       DN2DT = ZERO
       DN2DL = ZERO
       DN2DU = ZERO
       CALL becke05_neff(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,N2,
     $ DN2DR,DN2DG,DN2DT,DN2DL,DN2DU)
c       write(6,*)"beta neff: ", N2

       ! now build f factor
       ! appear in equation 13, 14 of JCP 2012 paper
       ! we will include smooth inside according to equation 43
       ! alpha part first
       G1     = ZERO 
       DG1DN1 = ZERO
       DG1DN2 = ZERO
       IF (DABS(N2)>THRESH) THEN
          G1     = (F_LIMIT-N1)/N2
          DG1DN1 = -ONE/N2
          IF (N2*N2>THRESH) THEN
             DG1DN2 = (N1-F_LIMIT)/(N2*N2)
          END IF
       END IF
       F1     = ZERO
       DF1DN1 = ZERO
       DF1DN2 = ZERO
       IF (G1>F_LIMIT+DELTA) THEN
          F1 = F_LIMIT
       ELSE IF (G1>=F_LIMIT-DELTA .and. G1<=F_LIMIT+DELTA) THEN
          F1     = F_LIMIT-(G1-F_LIMIT-DELTA)*(G1-F_LIMIT-DELTA)/
     $(FOUR*DELTA)
          DF1DG  = -(G1-F_LIMIT-DELTA)/(TWO*DELTA)
          DF1DN1 = DF1DG*DG1DN1
          DF1DN2 = DF1DG*DG1DN2
       ELSE
          F1     = G1
          DF1DN1 = DG1DN1
          DF1DN2 = DG1DN2
       END IF

       ! now get the derivatives
       DF1DRA = DF1DN1*DN1DR
       DF1DGA = DF1DN1*DN1DG
       DF1DTA = DF1DN1*DN1DT
       DF1DLA = DF1DN1*DN1DL
       DF1DUA = DF1DN1*DN1DU
       DF1DRB = DF1DN2*DN2DR
       DF1DGB = DF1DN2*DN2DG
       DF1DTB = DF1DN2*DN2DT
       DF1DLB = DF1DN2*DN2DL
       DF1DUB = DF1DN2*DN2DU
c       write(6,*)"F1: ", F1
c       write(6,*)"F1 alpha rho: ", DF1DN1*DN1DR
c       write(6,*)"F1 alpha GAA: ", DF1DN1*DN1DG
c       write(6,*)"F1 alpha tau: ", DF1DN1*DN1DT
c       write(6,*)"F1 alpha lap: ", DF1DN1*DN1DL
c       write(6,*)"F1 alpha ex : ", DF1DN1*DN1DU
c       write(6,*)"F1 beta  rho: ", DF1DN2*DN2DR
c       write(6,*)"F1 beta  GAA: ", DF1DN2*DN2DG
c       write(6,*)"F1 beta  tau: ", DF1DN2*DN2DT
c       write(6,*)"F1 beta  lap: ", DF1DN2*DN2DL
c       write(6,*)"F1 beta  ex : ", DF1DN2*DN2DU

       ! now it's beta part
       G2     = ZERO
       DG2DN1 = ZERO
       DG2DN2 = ZERO
       IF (DABS(N1)>THRESH) THEN
          G2     = (F_LIMIT-N2)/N1
          DG2DN2 = -ONE/N1
          IF (N1*N1>THRESH) THEN
             DG2DN1 = (N2-F_LIMIT)/(N1*N1)
          END IF
       END IF
       F2     = ZERO
       DF2DN1 = ZERO
       DF2DN2 = ZERO
       IF (G2>F_LIMIT+DELTA) THEN
          F2 = F_LIMIT
       ELSE IF (G2>=F_LIMIT-DELTA .and. G2<=F_LIMIT+DELTA) THEN
          F2     = F_LIMIT-(G2-F_LIMIT-DELTA)*(G2-F_LIMIT-DELTA)/
     $(FOUR*DELTA)
          DF2DG  = -(G2-F_LIMIT-DELTA)/(TWO*DELTA)
          DF2DN1 = DF2DG*DG2DN1
          DF2DN2 = DF2DG*DG2DN2
       ELSE
          F2     = G2
          DF2DN1 = DG2DN1
          DF2DN2 = DG2DN2
       END IF

       ! now get the derivatives
       DF2DRA = DF2DN1*DN1DR
       DF2DGA = DF2DN1*DN1DG
       DF2DTA = DF2DN1*DN1DT
       DF2DLA = DF2DN1*DN1DL
       DF2DUA = DF2DN1*DN1DU
       DF2DRB = DF2DN2*DN2DR
       DF2DGB = DF2DN2*DN2DG
       DF2DTB = DF2DN2*DN2DT
       DF2DLB = DF2DN2*DN2DL
       DF2DUB = DF2DN2*DN2DU
c       write(6,*)"F2: ", F2
c       write(6,*)"F2 alpha rho: ", DF2DN1*DN1DR
c       write(6,*)"F2 alpha GAA: ", DF2DN1*DN1DG
c       write(6,*)"F2 alpha tau: ", DF2DN1*DN1DT
c       write(6,*)"F2 alpha lap: ", DF2DN1*DN1DL
c       write(6,*)"F2 alpha ex : ", DF2DN1*DN1DU
c       write(6,*)"F2 beta  rho: ", DF2DN2*DN2DR
c       write(6,*)"F2 beta  GAA: ", DF2DN2*DN2DG
c       write(6,*)"F2 beta  tau: ", DF2DN2*DN2DT
c       write(6,*)"F2 beta  lap: ", DF2DN2*DN2DL
c       write(6,*)"F2 beta  ex : ", DF2DN2*DN2DU

       ! now we step into the smooth function part
       ! let;s deal with Z first
       Z = ZERO
       T = F1*F1+F2*F2
       DZDRA = ZERO 
       DZDGA = ZERO 
       DZDTA = ZERO 
       DZDLA = ZERO 
       DZDUA = ZERO 
       DZDRB = ZERO 
       DZDGB = ZERO 
       DZDTB = ZERO 
       DZDLB = ZERO 
       DZDUB = ZERO 
       IF (T > THRESH) THEN
          Z = (F1-F2)/T
          DZDRA = (DF1DRA-DF2DRA)/T
          DZDGA = (DF1DGA-DF2DGA)/T
          DZDTA = (DF1DTA-DF2DTA)/T
          DZDLA = (DF1DLA-DF2DLA)/T
          DZDUA = (DF1DUA-DF2DUA)/T
          DZDRB = (DF1DRB-DF2DRB)/T
          DZDGB = (DF1DGB-DF2DGB)/T
          DZDTB = (DF1DTB-DF2DTB)/T
          DZDLB = (DF1DLB-DF2DLB)/T
          DZDUB = (DF1DUB-DF2DUB)/T
       END IF
       IF (T*T > THRESH) THEN
          T2    = T*T
          DZDRA = DZDRA - (F1-F2)*(TWO*F1*DF1DRA+TWO*F2*DF2DRA)/T2
          DZDGA = DZDGA - (F1-F2)*(TWO*F1*DF1DGA+TWO*F2*DF2DGA)/T2
          DZDTA = DZDTA - (F1-F2)*(TWO*F1*DF1DTA+TWO*F2*DF2DTA)/T2
          DZDLA = DZDLA - (F1-F2)*(TWO*F1*DF1DLA+TWO*F2*DF2DLA)/T2
          DZDUA = DZDUA - (F1-F2)*(TWO*F1*DF1DUA+TWO*F2*DF2DUA)/T2
          DZDRB = DZDRB - (F1-F2)*(TWO*F1*DF1DRB+TWO*F2*DF2DRB)/T2
          DZDGB = DZDGB - (F1-F2)*(TWO*F1*DF1DGB+TWO*F2*DF2DGB)/T2
          DZDTB = DZDTB - (F1-F2)*(TWO*F1*DF1DTB+TWO*F2*DF2DTB)/T2
          DZDLB = DZDLB - (F1-F2)*(TWO*F1*DF1DLB+TWO*F2*DF2DLB)/T2
          DZDUB = DZDUB - (F1-F2)*(TWO*F1*DF1DUB+TWO*F2*DF2DUB)/T2
       END IF

       ! finally, we apply the smooth function to get the final f
       ! the original f is the min(f1,f2)
       ! one thing to note, that PZ below could be very large then
       ! it may break the expression of HPZ
       ! therefore we use log to test it
       LIMIT_Z = DLOG(ONE/THRESH)
       HPZ   = ZERO
       DHPZDZ= ZERO
       PZ    = P*Z
       IF (DABS(PZ) < LIMIT_Z) THEN
          T   = DEXP(P*Z)
          HPZ = ONE/(ONE+T)
          DHPZDZ = -P*HPZ*HPZ*T
       ELSE IF (PZ < ZERO) THEN
          HPZ   = ONE
       ELSE IF (PZ > ZERO) THEN
          HPZ   = ZERO
       END IF
c       write(6,*)"Z value: ", Z
c       write(6,*)"PZ value: ", PZ
c       write(6,*)"HPZ value: ", HPZ
c       write(6,*)"HPZ deriv value: ", DHPZDZ

       ! find difference check on z
c       PZ1 = PZ + PZ*0.001D0
c       HPZ1   = ZERO
c       IF (DABS(PZ1) < LIMIT_Z) THEN
c          T    = DEXP(PZ1)
c          HPZ1 = ONE/(ONE+T)
c       ELSE IF (PZ1 < ZERO) THEN
c          HPZ1   = ONE
c       ELSE IF (PZ1 > ZERO) THEN
c          HPZ1   = ZERO
c       END IF
c       PZ2 = PZ - PZ*0.001D0
c       HPZ2   = ZERO
c       IF (DABS(PZ2) < LIMIT_Z) THEN
c          T    = DEXP(PZ2)
c          HPZ2 = ONE/(ONE+T)
c       ELSE IF (PZ2 < ZERO) THEN
c          HPZ2   = ONE
c       ELSE IF (PZ2 > ZERO) THEN
c          HPZ2   = ZERO
c       END IF
c       IF (DABS(PZ*0.001D0)>THRESH) THEN
c          diff_HPZ = (HPZ1-HPZ2)/(PZ1-PZ2)
c          IF (DABS(DABS(diff_HPZ*P)-DABS(DHPZDZ))>0.00001D0) THEN
c             write(6,*)"HPZ find difference", diff_HPZ*P, DHPZDZ
c          END IF
c       END IF

       ! now let's collect all of terms together
       F    = (F1-F2)*HPZ + F2
       D1R1 = (DF1DRA-DF2DRA)*HPZ+(F1-F2)*DHPZDZ*DZDRA+DF2DRA
       D1R2 = (DF1DRB-DF2DRB)*HPZ+(F1-F2)*DHPZDZ*DZDRB+DF2DRB
       D1G1 = (DF1DGA-DF2DGA)*HPZ+(F1-F2)*DHPZDZ*DZDGA+DF2DGA
       D1G2 = (DF1DGB-DF2DGB)*HPZ+(F1-F2)*DHPZDZ*DZDGB+DF2DGB
       D1T1 = (DF1DTA-DF2DTA)*HPZ+(F1-F2)*DHPZDZ*DZDTA+DF2DTA
       D1T2 = (DF1DTB-DF2DTB)*HPZ+(F1-F2)*DHPZDZ*DZDTB+DF2DTB
       D1L1 = (DF1DLA-DF2DLA)*HPZ+(F1-F2)*DHPZDZ*DZDLA+DF2DLA
       D1L2 = (DF1DLB-DF2DLB)*HPZ+(F1-F2)*DHPZDZ*DZDLB+DF2DLB
       D1U1 = (DF1DUA-DF2DUA)*HPZ+(F1-F2)*DHPZDZ*DZDUA+DF2DUA
       D1U2 = (DF1DUB-DF2DUB)*HPZ+(F1-F2)*DHPZDZ*DZDUB+DF2DUB
c       write(6,*)"F value: ", F
c       write(6,*)"F alpha rho: ", D1R1
c       write(6,*)"F alpha GAA: ", D1G1
c       write(6,*)"F alpha tau: ", D1T1
c       write(6,*)"F alpha lap: ", D1L1
c       write(6,*)"F alpha ex : ", D1U1
c       write(6,*)"F beta  rho: ", D1R2
c       write(6,*)"F beta  GAA: ", D1G2
c       write(6,*)"F beta  tau: ", D1T2
c       write(6,*)"F beta  lap: ", D1L2
c       write(6,*)"F beta  ex : ", D1U2

       RETURN 
       END

      SUBROUTINE ori_becke05_f(THRESH,HIRWT,Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the f factor for the exchange of Becke05 functional *
c    *  the inplementation is based on the Becke's original paper:    *
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *  because this paper use min function to evaluate the f factor, *
c    *  therefore no derivatives is available for this situation      *
c    *                                                                *
c    *   input(1 is for alpha spin components, 2 is for beta)         *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   UX:  the exact exchange energy density                       *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)

C
C      here it's related to the smooth technique of f factor
C       
       parameter(F_LIMIT = 1.0D0  )
       parameter(DELTA   = 0.05D0 )
       parameter(P       = 115.0D0)

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0
       FIVE = 5.0D0
       F12  = 0.50d0 

C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! build Neff for alpha
       N1    = ZERO
       DN1DR = ZERO
       DN1DG = ZERO
       DN1DT = ZERO
       DN1DL = ZERO
       DN1DU = ZERO
       CALL becke05_neff(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,N1,
     $ DN1DR,DN1DG,DN1DT,DN1DL,DN1DU)
       IF (DABS(N1)>F_LIMIT) THEN
          N1 = F_LIMIT
       END IF
c       write(6,*)"alpha neff: ", N1

       ! build Neff for beta
       N2    = ZERO
       DN2DR = ZERO
       DN2DG = ZERO
       DN2DT = ZERO
       DN2DL = ZERO
       DN2DU = ZERO
       CALL becke05_neff(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,N2,
     $ DN2DR,DN2DG,DN2DT,DN2DL,DN2DU)
       IF (DABS(N2)>F_LIMIT) THEN
          N2 = F_LIMIT
       END IF
c       write(6,*)"beta neff: ", N2

       ! now build f factor
       G1     = ONE
       G2     = ONE 
       IF (DABS(N2)>THRESH) THEN
          G1     = (F_LIMIT-N1)/N2
       END IF
       IF (DABS(N1)>THRESH) THEN
          G2     = (F_LIMIT-N2)/N1
       END IF

       ! now this is the final f factor
       F = dmin1(G1,G2,ONE)

       RETURN 
       END


      SUBROUTINE becke05_A_sigma(iSpin,NELE,NTOLELE,THRESH,
     $ VAL_P, VAL_Q,
     $ HIRWT,Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,A,D1R1,D1G1,D1T1,D1L1,D1U1,
     $ D1R2,D1G2,D1T2,D1L2,D1U2)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the A factor for the exchange of Becke05 functional *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input(1 is for spin alpha, 2 is for spin beta)               *
c    *   iSpin: the result A's spin state. 1 is for alpha, 2 for beta *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   UX:  the exact exchange energy density                       *
c    *                                                                *
c    *   output(see the explanation in the input above)               *
c    *   A  :  this is the A factor defined in equation 15,16 and 47  *
c    *   D1R:  first derivative for A with Rho                        *
c    *   D1G:  first derivative for A with gamma                      *
c    *   D1T:  first derivative for A with tau                        *
c    *   D1L:  first derivative for A with Laplacian                  *
c    *   D1U:  first derivative for A with exchange energy density    *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)
       INTEGER ISPIN
       INTEGER NELE,NTOLELE

C
C      here it's related to the smooth technique of A factor
C       

C
C      similarly with the parameter of P, now the Q value can also
C      be varied.
C
C       parameter(Q = 120.0D0)

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0
       FIVE = 5.0D0
       F12  = 0.5D0
       F14  = 0.25D0
       F16  = 1.0D0/6.0D0
       F13  = ONE/THREE
       F112 = 1.0D0/12.0D0
       F110 = 1.0D0/10.0D0
       F05 = 0.070D0
       Q    = VAL_Q

C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! build Neff for alpha
       N1    = ZERO
       DN1DR = ZERO
       DN1DG = ZERO
       DN1DT = ZERO
       DN1DL = ZERO
       DN1DU = ZERO
       CALL becke05_neff(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,N1,
     $ DN1DR,DN1DG,DN1DT,DN1DL,DN1DU)

       ! build Neff for beta
       N2    = ZERO
       DN2DR = ZERO
       DN2DG = ZERO
       DN2DT = ZERO
       DN2DL = ZERO
       DN2DU = ZERO
       CALL becke05_neff(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,N2,
     $ DN2DR,DN2DG,DN2DT,DN2DL,DN2DU)

       ! build f factor
       F     = ZERO
       DFR1  = ZERO
       DFG1  = ZERO
       DFT1  = ZERO
       DFL1  = ZERO
       DFU1  = ZERO
       DFR2  = ZERO
       DFG2  = ZERO
       DFT2  = ZERO
       DFL2  = ZERO
       DFU2  = ZERO
       CALL becke05_f(THRESH,VAL_P,HIRWT,Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,F,DFR1,DFG1,DFT1,DFL1,DFU1,
     $ DFR2,DFG2,DFT2,DFL2,DFU2)

       ! build M2_{sigma}
       M2     = ZERO
       DM2R   = ZERO
       DM2G   = ZERO
       DM2T   = ZERO
       DM2L   = ZERO
       DM2U   = ZERO
       IF (ISpin .eq. 1) THEN
          CALL becke05_M2(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,M2,
     $ DM2R,DM2G,DM2T,DM2L,DM2U)
       ELSE
          CALL becke05_M2(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,M2,
     $ DM2R,DM2G,DM2T,DM2L,DM2U)
       END IF

       ! now make A1
       A1     = ZERO
       DA1R1  = ZERO
       DA1G1  = ZERO
       DA1T1  = ZERO
       DA1L1  = ZERO
       DA1U1  = ZERO
       DA1R2  = ZERO
       DA1G2  = ZERO
       DA1T2  = ZERO
       DA1L2  = ZERO
       DA1U2  = ZERO
       IF (NTOLELE .gt. 1) THEN
          IF (ISpin .eq. 1) THEN
             T  = ONE - N1 - F*N2
             IF (DABS(M2) > THRESH) THEN
                A1     = T/M2
                DA1R1  = (-DN1DR-DFR1*N2)/M2
                DA1G1  = (-DN1DG-DFG1*N2)/M2
                DA1T1  = (-DN1DT-DFT1*N2)/M2
                DA1L1  = (-DN1DL-DFL1*N2)/M2
                DA1U1  = (-DN1DU-DFU1*N2)/M2
                DA1R2  = -(DFR2*N2+F*DN2DR)/M2
                DA1G2  = -(DFG2*N2+F*DN2DG)/M2
                DA1T2  = -(DFT2*N2+F*DN2DT)/M2
                DA1L2  = -(DFL2*N2+F*DN2DL)/M2
                DA1U2  = -(DFU2*N2+F*DN2DU)/M2
             END IF
             M2_2 = M2*M2
             IF ((M2_2) > THRESH) THEN
                DA1R1  = DA1R1 - DM2R*T/M2_2
                DA1G1  = DA1G1 - DM2G*T/M2_2
                DA1T1  = DA1T1 - DM2T*T/M2_2
                DA1L1  = DA1L1 - DM2L*T/M2_2
                DA1U1  = DA1U1 - DM2U*T/M2_2
             END IF
          ELSE

             ! now please remember, this is for beta spin of A1
             T  = ONE - N2 - F*N1
             IF (DABS(M2) > THRESH) THEN
                A1     = T/M2
                DA1R1  = -(DFR1*N1+F*DN1DR)/M2
                DA1G1  = -(DFG1*N1+F*DN1DG)/M2
                DA1T1  = -(DFT1*N1+F*DN1DT)/M2
                DA1L1  = -(DFL1*N1+F*DN1DL)/M2
                DA1U1  = -(DFU1*N1+F*DN1DU)/M2
                DA1R2  = (-DN2DR-DFR2*N1)/M2
                DA1G2  = (-DN2DG-DFG2*N1)/M2
                DA1T2  = (-DN2DT-DFT2*N1)/M2
                DA1L2  = (-DN2DL-DFL2*N1)/M2
                DA1U2  = (-DN2DU-DFU2*N1)/M2
             END IF
             M2_2 = M2*M2
             IF ((M2_2) > THRESH) THEN
                DA1R2  = DA1R2 - DM2R*T/M2_2
                DA1G2  = DA1G2 - DM2G*T/M2_2
                DA1T2  = DA1T2 - DM2T*T/M2_2
                DA1L2  = DA1L2 - DM2L*T/M2_2
                DA1U2  = DA1U2 - DM2U*T/M2_2
             END IF
          END IF
       END IF
c       write(6,*)"A1 value: ", A1
c       write(6,*)"A1 alpha rho: ", DA1R1
c       write(6,*)"A1 alpha GAA: ", DA1G1
c       write(6,*)"A1 alpha tau: ", DA1T1
c       write(6,*)"A1 alpha lap: ", DA1L1
c       write(6,*)"A1 alpha ex : ", DA1U1
c       write(6,*)"A1 beta  rho: ", DA1R2
c       write(6,*)"A1 beta  GAA: ", DA1G2
c       write(6,*)"A1 beta  tau: ", DA1T2
c       write(6,*)"A1 beta  lap: ", DA1L2
c       write(6,*)"A1 beta  ex : ", DA1U2

       ! now for the A2 part
       A2    = ZERO
       DA2DR = ZERO
       DA2DG = ZERO
       DA2DT = ZERO
       IF (NELE .gt. 1) THEN
          IF (ISpin .eq. 1) THEN
             IF (Rho1 > THRESH) THEN
                D     = Tau1 - F14*Gam1/Rho1
                A2    = F05*F13*D/Rho1
                DA2DT = F05*F13/Rho1
                R2    = Rho1*Rho1
                IF (R2 > THRESH) THEN
                   DA2DR = -F05*F13*Tau1/R2
                   DA2DG = -F05*F112/R2
                   R3 = R2*Rho1
                   IF (R3 > THRESH) THEN
                      DA2DR = DA2DR + F05*F16*Gam1/R3
                   END IF
                END IF
             END IF
          ELSE
             IF (Rho2 > THRESH) THEN
                D     = Tau2 - F14*Gam2/Rho2
                A2    = F05*F13*D/Rho2
                DA2DT = F05*F13/Rho2
                R2    = Rho2*Rho2 
                IF (R2 > THRESH) THEN
                   DA2DR = -F05*F13*Tau2/R2
                   DA2DG = -F05*F112/R2
                   R3 = R2*Rho2
                   IF (R3 > THRESH) THEN
                      DA2DR = DA2DR + F05*F16*Gam2/R3
                   END IF
                END IF
             END IF
          END IF
       END IF
c       write(6,*)"A2 value: ", A2
c       write(6,*)"A2 Rho: ", DA2DR
c       write(6,*)"A2 GRho: ", DA2DG
c       write(6,*)"A2 tau: ", DA2DT

       ! finally, let's try to get A by smoothing the A1 and A2
       ! this is done by similar expression for F factor
       !
       ! we note, that if the QT is large enough (both positive
       ! or negative enough), we have to consider them in seperate 
       ! way. however, their derivatives are all zero
       !
       HQT    = ZERO
       DHQTDT = ZERO
       T     = A1-A2
       QT    = Q*T
       LIMIT_QT = DLOG(ONE/THRESH)
       IF (DABS(QT)<LIMIT_QT) THEN

          ! expression for HQT
          EXQT = DEXP(QT)
          HQT  = ONE/(ONE+EXQT) 
          DHQTDT = -Q*HQT*HQT*EXQT
       ELSE IF (QT < ZERO) THEN
          HQT  = ONE
       ELSE IF (QT > ZERO) THEN
          HQT  = ZERO
       END IF
c       write(6,*)"T value: ", T
c       write(6,*)"QT value: ", QT
c       write(6,*)"HQT value: ", HQT
c       write(6,*)"HQT deriv value: ", DHQTDT
          
       ! derivatives
       DHR1  = ZERO
       DHG1  = ZERO
       DHT1  = ZERO
       DHL1  = ZERO
       DHU1  = ZERO
       DHR2  = ZERO
       DHG2  = ZERO
       DHT2  = ZERO
       DHL2  = ZERO
       DHU2  = ZERO
       IF (ISPIN .eq. 1) THEN
          DHR1  = DHQTDT*(DA1R1-DA2DR)
          DHG1  = DHQTDT*(DA1G1-DA2DG) 
          DHT1  = DHQTDT*(DA1T1-DA2DT) 
          DHL1  = DHQTDT*DA1L1 
          DHU1  = DHQTDT*DA1U1 
          DHR2  = DHQTDT*DA1R2 
          DHG2  = DHQTDT*DA1G2 
          DHT2  = DHQTDT*DA1T2 
          DHL2  = DHQTDT*DA1L2 
          DHU2  = DHQTDT*DA1U2 
       ELSE
          DHR1  = DHQTDT*DA1R1
          DHG1  = DHQTDT*DA1G1
          DHT1  = DHQTDT*DA1T1
          DHL1  = DHQTDT*DA1L1 
          DHU1  = DHQTDT*DA1U1 
          DHR2  = DHQTDT*(DA1R2-DA2DR) 
          DHG2  = DHQTDT*(DA1G2-DA2DG) 
          DHT2  = DHQTDT*(DA1T2-DA2DT) 
          DHL2  = DHQTDT*DA1L2 
          DHU2  = DHQTDT*DA1U2 
       END IF

       ! bring all of derivatives together
       D1R1  = ZERO
       D1G1  = ZERO
       D1T1  = ZERO
       D1L1  = ZERO
       D1U1  = ZERO
       D1R2  = ZERO
       D1G2  = ZERO
       D1T2  = ZERO
       D1L2  = ZERO
       D1U2  = ZERO
       A = (A1-A2)*HQT+A2
       IF (ISpin .eq. 1) THEN
          D1R1 = (DA1R1-DA2DR)*HQT+(A1-A2)*DHR1+DA2DR
          D1G1 = (DA1G1-DA2DG)*HQT+(A1-A2)*DHG1+DA2DG
          D1T1 = (DA1T1-DA2DT)*HQT+(A1-A2)*DHT1+DA2DT
          D1L1 = DA1L1*HQT+(A1-A2)*DHL1
          D1U1 = DA1U1*HQT+(A1-A2)*DHU1
          D1R2 = DA1R2*HQT+(A1-A2)*DHR2
          D1G2 = DA1G2*HQT+(A1-A2)*DHG2
          D1T2 = DA1T2*HQT+(A1-A2)*DHT2
          D1L2 = DA1L2*HQT+(A1-A2)*DHL2
          D1U2 = DA1U2*HQT+(A1-A2)*DHU2
       ELSE
          D1R1 = DA1R1*HQT+(A1-A2)*DHR1
          D1G1 = DA1G1*HQT+(A1-A2)*DHG1
          D1T1 = DA1T1*HQT+(A1-A2)*DHT1
          D1L1 = DA1L1*HQT+(A1-A2)*DHL1
          D1U1 = DA1U1*HQT+(A1-A2)*DHU1
          D1R2 = (DA1R2-DA2DR)*HQT+(A1-A2)*DHR2+DA2DR
          D1G2 = (DA1G2-DA2DG)*HQT+(A1-A2)*DHG2+DA2DG
          D1T2 = (DA1T2-DA2DT)*HQT+(A1-A2)*DHT2+DA2DT
          D1L2 = DA1L2*HQT+(A1-A2)*DHL2
          D1U2 = DA1U2*HQT+(A1-A2)*DHU2
       END IF
C       write(6,*)"A value: ", A
C       write(6,*)"A alpha rho: ", D1R1
C       write(6,*)"A alpha GAA: ", D1G1
C       write(6,*)"A alpha tau: ", D1T1
C       write(6,*)"A alpha lap: ", D1L1
C       write(6,*)"A alpha ex : ", D1U1
C       write(6,*)"A beta  rho: ", D1R2
C       write(6,*)"A beta  GAA: ", D1G2
C       write(6,*)"A beta  tau: ", D1T2
C       write(6,*)"A beta  lap: ", D1L2
C       write(6,*)"A beta  ex : ", D1U2

       RETURN 
       END

      SUBROUTINE ori_becke05_A_sigma(iSpin,NELE,NTOTELE,THRESH,HIRWT,
     $ Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,A)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the A factor for the exchange of Becke05 functional *
c    *  the inplementation is based on the Becke's original paper:    *
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input(1 is for spin alpha, 2 is for spin beta)               *
c    *   iSpin: the result A's spin state. 1 is for alpha, 2 for beta *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   UX:  the exact exchange energy density                       *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
       IMPLICIT REAL*8(A-Z)
       INTEGER ISPIN
       INTEGER NELE,NTOTELE

C
C      here it's related to the smooth technique of A factor
C       
       parameter(Q = 120.0D0)
       parameter(F_LIMIT = 1.0D0  )

       ! constants
       PI   = 4.D0*DATAN(1.D0)
       ZERO = 0.0D0
       ONE  = 1.0D0
       TWO  = 2.0D0
       THREE= 3.0D0
       FOUR = 4.0D0
       FIVE = 5.0D0
       F12  = 0.5D0
       F14  = 0.25D0
       F16  = 1.0D0/6.0D0
       F13  = ONE/THREE
       F112 = 1.0D0/12.0D0

C
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C
       ! build Neff for alpha
       N1    = ZERO
       DN1DR = ZERO
       DN1DG = ZERO
       DN1DT = ZERO
       DN1DL = ZERO
       DN1DU = ZERO
       CALL becke05_neff(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,N1,
     $ DN1DR,DN1DG,DN1DT,DN1DL,DN1DU)
       IF (DABS(N1)>F_LIMIT) THEN
          N1 = F_LIMIT
       END IF
c       write(6,*)"alpha neff: ", N1

       ! build Neff for beta
       N2    = ZERO
       DN2DR = ZERO
       DN2DG = ZERO
       DN2DT = ZERO
       DN2DL = ZERO
       DN2DU = ZERO
       CALL becke05_neff(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,N2,
     $ DN2DR,DN2DG,DN2DT,DN2DL,DN2DU)
       IF (DABS(N2)>F_LIMIT) THEN
          N2 = F_LIMIT
       END IF
c       write(6,*)"beta neff: ", N2

       ! build ori f factor
       F     = ZERO
       CALL ori_becke05_f(THRESH,HIRWT,Rho1,Gam1,Tau1,Lap1,UX1,
     $ Rho2,Gam2,Tau2,Lap2,UX2,F)

       ! build M2_{sigma}
       M2     = ZERO
       DM2R   = ZERO
       DM2G   = ZERO
       DM2T   = ZERO
       DM2L   = ZERO
       DM2U   = ZERO
       IF (ISpin .eq. 1) THEN
          CALL becke05_M2(THRESH,Rho1,Gam1,Tau1,Lap1,UX1,HIRWT,M2,
     $ DM2R,DM2G,DM2T,DM2L,DM2U)
       ELSE
          CALL becke05_M2(THRESH,Rho2,Gam2,Tau2,Lap2,UX2,HIRWT,M2,
     $ DM2R,DM2G,DM2T,DM2L,DM2U)
       END IF

       ! now make A1
       A1     = ZERO
       IF (NTOTELE .gt. 1) THEN
          IF (ISpin .eq. 1) THEN
             T  = ONE - N1 - F*N2
             IF (DABS(M2) > THRESH) THEN
                A1     = T/M2
             END IF
          ELSE
             ! now please remember, this is for beta spin of A1
             T  = ONE - N2 - F*N1
             IF (DABS(M2) > THRESH) THEN
                A1     = T/M2
             END IF
          END IF
       END IF

       ! now for the A2 part
       A2    = ZERO
       IF (NELE .gt. 1) THEN
          IF (ISpin .eq. 1) THEN
             IF (Rho1 > THRESH) THEN
                D     = Tau1 - F14*Gam1/Rho1
                A2    = F13*D/Rho1
             END IF
          ELSE
             IF (Rho2 > THRESH) THEN
                D     = Tau2 - F14*Gam2/Rho2
                A2    = F13*D/Rho2
             END IF
          END IF
       END IF
C       write(6,*)"A2 value: ", A2
C       write(6,*)"A2 Rho: ", DA2DR
C       write(6,*)"A2 GRho: ", DA2DG
C       write(6,*)"A2 tau: ", DA2DT

       A = DMIN1(A1,A2)
       RETURN 
       END

      SUBROUTINE Becke05_odd_electron(NG,NA,NB,THRESH,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,f1,f2,F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the odd electron population based on the paper:     *
c    *  Physical Review A 88 032510(2013)                             *
c    *                                                                *
c    *   input(1 is for spin alpha, 2 is for spin beta)               *
c    *   THRESH: threshold to determine the small rho etc.            *         
c    *   Rho: spin density (alpha rho or beta rho)                    *
c    *   Gam: |\nabla \rho|(spin polarized)                           *
c    *   Tau: kinetic density(without factor of 1/2)                  *
c    *   Lap: Laplacian of density                                    *
c    *   UX:  the exact exchange energy density                       *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NG,I
      INTEGER NA,NB,NE
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F(NG)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 N1,N2,DNDR,DNDG,DNDT,DNDL,DNDU
      REAL*8 M2A,M2B,DMR,DMG,DMT,DML,DMU
      REAL*8 A1,A2
      REAL*8 f1,f2
      REAL*8 OD1,OD2
      REAL*8 FAC
      REAL*8 HIRWT
      REAL*8 ZERO,ONE,TWO,FOUR

      REAL*8 F_LIMIT
      parameter(F_LIMIT = 1.0D0  )
      parameter(ZERO    = 0.0D0  )
      parameter(ONE     = 1.0D0  )
      parameter(TWO     = 2.0D0  )
      parameter(FOUR    = 4.0D0  )

      ! set up total number of electrons
      NE = NA + NB

      ! loop over NG
      DO I = 1,NG

C         write(6,*)"     "
C         write(6,*)"----------------------"
C         write(6,*)"for grid: ", i
          
         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)
         
         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! build Neff for alpha
         N1    = ZERO
         DNDR  = ZERO
         DNDG  = ZERO
         DNDT  = ZERO
         DNDL  = ZERO
         DNDU  = ZERO
         CALL becke05_neff(THRESH,RA,GAA,TA,LA,UA,HIRWT,N1,
     $   DNDR,DNDG,DNDT,DNDL,DNDU)
         IF (DABS(N1)>F_LIMIT) THEN
            N1 = F_LIMIT
         END IF
c         write(6,*)"alpha neff: ", N1

         ! build Neff for beta
         N2    = ZERO
         DNDR  = ZERO
         DNDG  = ZERO
         DNDT  = ZERO
         DNDL  = ZERO
         DNDU  = ZERO
         CALL becke05_neff(THRESH,RB,GBB,TB,LB,UB,HIRWT,N2,
     $   DNDR,DNDG,DNDT,DNDL,DNDU)
         IF (DABS(N2)>F_LIMIT) THEN
            N2 = F_LIMIT
         END IF
c         write(6,*)"beta neff: ", N2

         ! build ori f factor
         FAC     = ZERO
         CALL ori_becke05_f(THRESH,HIRWT,RA,GAA,TA,LA,UA,
     $   RB,GBB,TB,LB,UB,FAC)

         ! build M2_{sigma}
         M2A   = ZERO
         M2B   = ZERO
         DMR   = ZERO
         DMG   = ZERO
         DMT   = ZERO
         DML   = ZERO
         DMU   = ZERO
         CALL becke05_M2(THRESH,RA,GAA,TA,LA,UA,HIRWT,M2A,
     $   DMR,DMG,DMT,DML,DMU)
         CALL becke05_M2(THRESH,RB,GBB,TB,LB,UB,HIRWT,M2B,
     $   DMR,DMG,DMT,DML,DMU)

         ! build A1 and A2
         A1 = ZERO
         A2 = ZERO
         CALL ori_becke05_A_sigma(1,NA,NE,THRESH,HIRWT,RA,GAA,TA,LA,UA,
     $   RB,GBB,TB,LB,UB,A1)
         CALL ori_becke05_A_sigma(2,NB,NE,THRESH,HIRWT,RA,GAA,TA,LA,UA,
     $   RB,GBB,TB,LB,UB,A2)

         ! now it's the odd electron population now
         ! alpha
         OD1  = f1*FOUR*RA*FAC*N2 + f2*TWO*RA*A1*M2A
         OD2  = f1*FOUR*RB*FAC*N1 + f2*TWO*RB*A2*M2B
c    Change here for computing atomic charges: 
c    comment the two lines below and uncomment the two olines above for 
c    true odd el population
c         OD1 = RA
c         OD2 = RB
         F(i) = OD1 + OD2
      END DO
      RETURN 
      END

      SUBROUTINE Becke05_NDOP(INFOR,NDEN,NG,THRESH,VAL_P,HIRWTS,
     $      RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic opposite spin part of Becke05       *
c    *  exchange functional                                           *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 HIRWT
      REAL*8 VAL_P
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! tmp variables 
      ! the f factor and it's derivatives
      REAL*8 FAC   
      REAL*8 DFDRA 
      REAL*8 DFDGAA
      REAL*8 DFDTA 
      REAL*8 DFDLA 
      REAL*8 DFDUA 
      REAL*8 DFDRB 
      REAL*8 DFDGBB
      REAL*8 DFDTB 
      REAL*8 DFDLB 
      REAL*8 DFDUB 

      ! now for the real derivatives
      REAL*8 TEMP
      REAL*8 DRhoA 
      REAL*8 DGAA
      REAL*8 DTA 
      REAL*8 DLA 
      REAL*8 DUA 
      REAL*8 DRhoB 
      REAL*8 DGBB
      REAL*8 DTB 
      REAL*8 DLB 
      REAL*8 DUB 
      REAL*8 DGAB
      REAL*8 VA
      REAL*8 VB

      ! constant
      REAL*8 ZERO
      REAL*8 ONE
      REAL*8 F12
C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! constants
      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0

      ! loop over NG
      DO I = 1,NG

c         write(6,*)"     "
c         write(6,*)"----------------------"
c         write(6,*)"for grid: ", i
          
         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)

         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! debug
c         write(6,*)"our RhoA: ",RA
c         write(6,*)"our GRhoA: ",GAA 
c         write(6,*)"our TA: ",TA
c         write(6,*)"our LA: ",LA 
c         write(6,*)"our UA: ",UA
c         write(6,*)"our RhoB: ",RB
c         write(6,*)"our GRhoB: ",GBB 
c         write(6,*)"our TB: ",TB
c         write(6,*)"our LB: ",LB 
c         write(6,*)"our UB: ",UB 

         ! now derive the f factor
         FAC    = ZERO
         DFDRA  = ZERO
         DFDGAA = ZERO
         DFDTA  = ZERO
         DFDLA  = ZERO
         DFDUA  = ZERO
         DFDRB  = ZERO
         DFDGBB = ZERO
         DFDTB  = ZERO
         DFDLB  = ZERO
         DFDUB  = ZERO
         CALL becke05_f(THRESH,VAL_P,HIRWT,RA,GAA,TA,LA,UA,
     $RB,GBB,TB,LB,UB,
     $FAC,DFDRA,DFDGAA,DFDTA,DFDLA,DFDUA,DFDRB,DFDGBB,DFDTB,DFDLB,DFDUB)

         ! remember, here we use potential for exchange energy
         ! therefore it should be divided by rho
         VA = ZERO
         VB = ZERO
         IF (RA > THRESH) THEN
            VA = UA/RA
         END IF
         IF (RB > THRESH) THEN
            VB = UB/RB
         END IF

         ! now compose the finaly energy
         F(i) = F(i) + F12*FAC*(RA*VB+RB*VA)

         ! derivatives
         TEMP = RA*VB+RB*VA
         DGAA = F12*DFDGAA*TEMP
         DGBB = F12*DFDGBB*TEMP
         DTA  = F12*DFDTA*TEMP
         DTB  = F12*DFDTB*TEMP
         DLA  = F12*DFDLA*TEMP
         DLB  = F12*DFDLB*TEMP

         ! derivatives for UA and UB
         DUA  = F12*DFDUA*TEMP 
         DRhoA= F12*DFDRA*TEMP + F12*FAC*VB
         IF (RA > THRESH) THEN
            DUA  = DUA + F12*FAC*(RB/RA)
            IF (RA*RA > THRESH) THEN
               DRhoA= DRhoA - F12*FAC*(RB/RA)*(UA/RA)
            END IF
         END IF
         DUB  = F12*DFDUB*TEMP 
         DRhoB= F12*DFDRB*TEMP + F12*FAC*VA
         IF (RB > THRESH) THEN
            DUB  = DUB + F12*FAC*RA/RB
            IF (RB*RB > THRESH) THEN
               DRhoB= DRhoB - F12*FAC*(RA/RB)*(UB/RB)
            END IF
         END IF
c         write(6,*)"functional derivatives for NDOP rhoA:", DRhoA
c         write(6,*)"functional derivatives for NDOP rhoB:", DRhoB
c         write(6,*)"functional derivatives for NDOP GAA :", DGAA
c         write(6,*)"functional derivatives for NDOP GBB :", DGBB
c         write(6,*)"functional derivatives for NDOP DTA :", DTA
c         write(6,*)"functional derivatives for NDOP DTB :", DTB
c         write(6,*)"functional derivatives for NDOP DLA :", DLA
c         write(6,*)"functional derivatives for NDOP DLB :", DLB
c         write(6,*)"functional derivatives for NDOP DUA :", DUA
c         write(6,*)"functional derivatives for NDOP DUB :", DUB

         ! collect all of terms for derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + DRhoA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + DGAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = ZERO
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + DTA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + DLA
         ID_EXA_POS=D1VARS(ID_EXA)
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + DUA
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + DRhoB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + DGBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + DTB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + DLB
         ID_EXB_POS=D1VARS(ID_EXB)
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + DUB

      END DO
      RETURN 
      END
         
      SUBROUTINE Becke05_NDPAR1(INFOR,NDEN,NG,NA,NB,THRESH,c_ndpar_cap,
     $                          VAL_P,VAL_Q,HIRWTS,RhoA,
     $                          RhoB,DRA,DRB,TauA,TauB,
     $                          LapA,LapB,EXA,EXB,F,TMP_D1F)
c
c    ******************************************************************
c    *                                                                *
c     This is the same as _NDPAR except that the parallel contribution
c     is scaled on the total. See paper EP&JK jctc 2021.
c     It is called by kp14ec.
c
c
c    *  evaluates the non-dynamic parallel spin part of Becke05       *
c    *  exchange functional with capped from above NDpar              *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    * Jing Kong and E. Proynov                                       *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I,J
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH,c_ndpar_cap
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F,TMP_D1F(*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 HIRWT
      REAL*8 VAL_P
      REAL*8 VAL_Q
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! tmp variables 
      ! the A factor and it's derivatives
      REAL*8 TMP2_D1F(N_FUNC_DERIV_1)
      REAL*8 usopp,usopp_RA,usopp_RB,usopp_GAA,usopp_GBB
      REAL*8 usopp_TA,usopp_TB,usopp_LA,usopp_LB,usopp_UA,usopp_UB
      REAL*8 uspar,uspar_RA,uspar_RB,uspar_GAA,uspar_GBB
      REAL*8 uspar_TA,uspar_TB,uspar_LA,uspar_LB,uspar_UA,uspar_UB
      REAL*8 DF1DRA,DF2DRA,DF1DRB,DF2DRB,DF1DGA,DF2DGA,DF1DGB,DF2DGB
      REAL*8 DF1DTA,DF2DTA,DF1DTB,DF2DTB,DF1DLA,DF2DLA,DF1DLB,DF2DLB
      REAL*8 DF1DUA,DF2DUA,DF1DUB,DF2DUB,HPZ,DHPZDZ,DZDRA,DZDRB
      REAL*8 DZDGA,DZDGB,DZDTA,DZDTB,DZDLA,DZDLB,DZDUA,DZDUB
      REAL*8 P,Z,PZ,FF,F1,F2,T,T1,T2,D1R1,D1R2,D1G1,D1G2,D1T1,D1T2
      REAL*8 D1L1,D1L2,D1U1,D1U2,TP,DZ1,DZ2,LIMIT_Z

      ! constant
      REAL*8 ZERO
      REAL*8 ONE
      REAL*8 TWO 
      REAL*8 F12
      REAL*8 dtol
c
      dtol = 1.0d-14
      ZERO = 0.0D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      F12  = c_ndpar_cap 
      P    = 500.0D0
      FF = ZERO
      F1= ZERO
      F2= ZERO
      F = ZERO
c
        CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      DO i = 1, NG 
          usopp = ZERO
          DO J = 1, N_FUNC_DERIV_1
            TMP2_D1F(J) = ZERO
          END DO

        call Becke05_NDOP(INFOR,NDEN,1,THRESH,VAL_P,HIRWTS,
     $                  RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,
     $                  EXA,EXB,usopp,TMP2_D1F)
         ID_RA_POS  = D1VARS(ID_RA)
         ID_GAA_POS = D1VARS(ID_GAA)
         ID_GAB_POS = D1VARS(ID_GAB)
         ID_TA_POS  = D1VARS(ID_TA)
         ID_LA_POS  = D1VARS(ID_LA)
         ID_EXA_POS = D1VARS(ID_EXA)
         ID_RB_POS  = D1VARS(ID_RB)
         ID_GBB_POS = D1VARS(ID_GBB)
         ID_TB_POS  = D1VARS(ID_TB)
         ID_LB_POS  = D1VARS(ID_LB)
         ID_EXB_POS = D1VARS(ID_EXB)
c
         usopp_RA   = TMP2_D1F(ID_RA_POS)
         usopp_GAA  = TMP2_D1F(ID_GAA_POS)
         usopp_TA   = TMP2_D1F(ID_TA_POS)
         usopp_LA   = TMP2_D1F(ID_LA_POS)
         usopp_UA   = TMP2_D1F(ID_EXA_POS)
         usopp_RB   = TMP2_D1F(ID_RB_POS)
         usopp_GBB  = TMP2_D1F(ID_GBB_POS)
         usopp_TB   = TMP2_D1F(ID_TB_POS)
         usopp_LB   = TMP2_D1F(ID_LB_POS)
         usopp_UB   = TMP2_D1F(ID_EXB_POS)
c
         uspar=ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP2_D1F(J) = ZERO
         END DO
c
       call  Becke05_NDPAR(INFOR,NDEN,1,NA,NB,THRESH,
     $                     VAL_P,VAL_Q,HIRWTS,
     $                     RhoA,RhoB,DRA,DRB,TauA,TauB,
     $                     LapA,LapB,EXA,EXB,uspar,TMP2_D1F)
c
         uspar_RA   = TMP2_D1F(ID_RA_POS)
         uspar_GAA  = TMP2_D1F(ID_GAA_POS)
         uspar_TA   = TMP2_D1F(ID_TA_POS)
         uspar_LA   = TMP2_D1F(ID_LA_POS)
         uspar_UA   = TMP2_D1F(ID_EXA_POS)
         uspar_RB   = TMP2_D1F(ID_RB_POS)
         uspar_GBB  = TMP2_D1F(ID_GBB_POS)
         uspar_TB   = TMP2_D1F(ID_TB_POS)
         uspar_LB   = TMP2_D1F(ID_LB_POS)
         uspar_UB   = TMP2_D1F(ID_EXB_POS)
       ! let;s deal with Z first
         Z = ZERO
         F1= uspar
         F2 = F12*usopp
c        IF(usopp > ZERO) F2=ZERO 
c        Change here. It was found the line below helps
c        convergence a bit for stretched Ar2+.
c         IF(usopp > -dtol) F2= -dtol
c        IF (usopp < ZERO.and.uspar < ZERO ) then
         IF (usopp < ZERO ) then
         T = F1*F1+F2*F2
c        write(6,*)"T =   ", T
c
       DF1DRA = uspar_RA 
       DF1DGA = uspar_GAA
       DF1DTA = uspar_TA 
       DF1DLA = uspar_LA 
       DF1DUA = uspar_UA 
       DF1DRB = uspar_RB 
       DF1DGB = uspar_GBB
       DF1DTB = uspar_TB 
       DF1DLB = uspar_LB 
       DF1DUB = uspar_UB 
       DF2DRA = F12*usopp_RA 
       DF2DGA = F12*usopp_GAA
       DF2DTA = F12*usopp_TA 
       DF2DLA = F12*usopp_LA 
       DF2DUA = F12*usopp_UA 
       DF2DRB = F12*usopp_RB 
       DF2DGB = F12*usopp_GBB
       DF2DTB = F12*usopp_TB 
       DF2DLB = F12*usopp_LB 
       DF2DUB = F12*usopp_UB 
c
       DZDRA = ZERO 
       DZDGA = ZERO 
       DZDTA = ZERO 
       DZDLA = ZERO 
       DZDUA = ZERO 
       DZDRB = ZERO 
       DZDGB = ZERO 
       DZDTB = ZERO 
       DZDLB = ZERO 
       DZDUB = ZERO 
       D1R1 = ZERO
       D1R2 = ZERO
       D1G1 = ZERO
       D1G2 = ZERO
       D1T1 = ZERO
       D1T2 = ZERO
       D1L1 = ZERO
       D1L2 = ZERO
       D1U1 = ZERO
       D1U2 = ZERO
       FF = ZERO
       IF (T > THRESH) THEN
         T1= DSQRT(T)
         Z = (F2-F1)/T1
         T2 = T1**3
         DZ1= -F2*(F1+F2)/T2
         DZ2= F1*(F1+F2)/T2 
         DZDRA=DZ1*DF1DRA + DZ2*DF2DRA
         DZDGA=DZ1*DF1DGA + DZ2*DF2DGA
         DZDTA=DZ1*DF1DTA + DZ2*DF2DTA
         DZDLA=DZ1*DF1DLA + DZ2*DF2DLA
         DZDUA=DZ1*DF1DUA + DZ2*DF2DUA
         DZDRB=DZ1*DF1DRB + DZ2*DF2DRB
         DZDGB=DZ1*DF1DGB + DZ2*DF2DGB
         DZDTB=DZ1*DF1DTB + DZ2*DF2DTB
         DZDLB=DZ1*DF1DLB + DZ2*DF2DLB
         DZDUB=DZ1*DF1DUB + DZ2*DF2DUB
       END IF

       ! finally, we apply the smooth function to get the final f
       ! the original f is the min(f1,f2)
       ! one thing to note, that PZ below could be very large then
       ! it may break the expression of HPZ
       ! therefore we use log to test it
       LIMIT_Z = DLOG(ONE/THRESH)
c      LIMIT_Z = 99.0d0
       HPZ   = ZERO
       DHPZDZ= ZERO
       PZ    = P*Z
       IF (DABS(PZ) < LIMIT_Z) THEN
          TP   = DEXP(P*Z)
          HPZ = ONE/(ONE+TP)
          DHPZDZ = -P*HPZ*HPZ*TP
       ELSE IF (PZ < ZERO) THEN
          HPZ   = ONE
       ELSE IF (PZ > ZERO) THEN
          HPZ   = ZERO
       END IF
c
       FF    = (F1-F2)*HPZ + F2
c      if (.not.(abs(FF-uspar) < 1D-5 .or. abs(FF-F12*usopp) < 1D-5)) 
c    $ then
c          write(*,*) "catching this", FF, "  ", uspar, "  ", usopp
c      endif 

       ! now let's collect all of terms together
C         write(6,*)"uspar ",uspar, "
       D1R1 = (DF1DRA-DF2DRA)*HPZ+(F1-F2)*DHPZDZ*DZDRA+DF2DRA
       D1R2 = (DF1DRB-DF2DRB)*HPZ+(F1-F2)*DHPZDZ*DZDRB+DF2DRB
       D1G1 = (DF1DGA-DF2DGA)*HPZ+(F1-F2)*DHPZDZ*DZDGA+DF2DGA
       D1G2 = (DF1DGB-DF2DGB)*HPZ+(F1-F2)*DHPZDZ*DZDGB+DF2DGB
       D1T1 = (DF1DTA-DF2DTA)*HPZ+(F1-F2)*DHPZDZ*DZDTA+DF2DTA
       D1T2 = (DF1DTB-DF2DTB)*HPZ+(F1-F2)*DHPZDZ*DZDTB+DF2DTB
       D1L1 = (DF1DLA-DF2DLA)*HPZ+(F1-F2)*DHPZDZ*DZDLA+DF2DLA
       D1L2 = (DF1DLB-DF2DLB)*HPZ+(F1-F2)*DHPZDZ*DZDLB+DF2DLB
       D1U1 = (DF1DUA-DF2DUA)*HPZ+(F1-F2)*DHPZDZ*DZDUA+DF2DUA
       D1U2 = (DF1DUB-DF2DUB)*HPZ+(F1-F2)*DHPZDZ*DZDUB+DF2DUB
c
c      if (.not.(abs(D1R1-DF1DRA)< 1D-5.or.abs(D1R1-DF2DRA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1R1," ", DF1DRA," ",DF2DRA 
c      endif
c
c      if (.not.(abs(D1G1-DF1DGA)< 1D-5.or.abs(D1G1-DF2DGA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1G1," ", DF1DGA," ",DF2DGA
c      endif
c
c      if (.not.(abs(D1T1-DF1DGA)< 1D-5.or.abs(D1T1-DF2DTA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1T1," ", DF1DTA," ",DF2DTA
c      endif
c
c      if (.not.(abs(D1L1-DF1DLA)< 1D-5.or.abs(D1L1-DF2DLA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1L1," ", DF1DLA," ",DF2DLA
c      endif
c
c      if (.not.(abs(D1U1-DF1DUA)< 1D-5.or.abs(D1U1-DF2DUA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1U1," ", DF1DUA," ",DF2DUA
c      endif
c
c       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
c
c        ID_RA_POS  = D1VARS(ID_RA)
c        ID_GAA_POS = D1VARS(ID_GAA)
c        ID_GAB_POS = D1VARS(ID_GAB)
c        ID_TA_POS  = D1VARS(ID_TA)
c        ID_LA_POS  = D1VARS(ID_LA)
c        ID_EXA_POS = D1VARS(ID_EXA)
c        ID_RB_POS  = D1VARS(ID_RB)
c        ID_GBB_POS = D1VARS(ID_GBB)
c        ID_TB_POS  = D1VARS(ID_TB)
c        ID_LB_POS  = D1VARS(ID_LB)
c        ID_EXB_POS = D1VARS(ID_EXB)
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
          F= F + FF
         TMP_D1F(ID_RA_POS)  = TMP_D1F(ID_RA_POS) + D1R1 
         TMP_D1F(ID_GAA_POS) = TMP_D1F(ID_GAA_POS) + D1G1 
         TMP_D1F(ID_GAB_POS) = ZERO
         TMP_D1F(ID_TA_POS)  = TMP_D1F(ID_TA_POS) + D1T1
         TMP_D1F(ID_LA_POS)  = TMP_D1F(ID_LA_POS) + D1L1 
         TMP_D1F(ID_EXA_POS) = TMP_D1F(ID_EXA_POS) + D1U1 
         TMP_D1F(ID_RB_POS)  = TMP_D1F(ID_RB_POS) +  D1R2 
         TMP_D1F(ID_GBB_POS) = TMP_D1F(ID_GBB_POS) + D1G2 
         TMP_D1F(ID_TB_POS)  = TMP_D1F(ID_TB_POS) + D1T2 
         TMP_D1F(ID_LB_POS)  = TMP_D1F(ID_LB_POS) + D1L2 
         TMP_D1F(ID_EXB_POS) = TMP_D1F(ID_EXB_POS) + D1U2 
          ELSE IF (uspar < ZERO) then
c            F = F + 0.80d0*uspar
             F = F + uspar
          END IF  
       END DO
       RETURN
       END
          
      SUBROUTINE Becke05_NDPAR(INFOR,NDEN,NG,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,D1F)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic parallel spin part of Becke05       *
c    *  exchange functional                                           *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 HIRWT
      REAL*8 VAL_P
      REAL*8 VAL_Q
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! tmp variables 
      ! the A factor and it's derivatives
      REAL*8 A   
      REAL*8 DADRA 
      REAL*8 DADGAA
      REAL*8 DADTA 
      REAL*8 DADLA 
      REAL*8 DADUA 
      REAL*8 DADRB 
      REAL*8 DADGBB
      REAL*8 DADTB 
      REAL*8 DADLB 
      REAL*8 DADUB 

      ! tmp variables 
      ! the M1 and it's derivatives
      REAL*8 M1   
      REAL*8 DM1DRA 
      REAL*8 DM1DGAA
      REAL*8 DM1DTA 
      REAL*8 DM1DLA 
      REAL*8 DM1DUA 
      REAL*8 DM1DRB 
      REAL*8 DM1DGBB
      REAL*8 DM1DTB 
      REAL*8 DM1DLB 
      REAL*8 DM1DUB 

      ! now for the real derivatives
      REAL*8 DRhoA 
      REAL*8 DGAA
      REAL*8 DTA 
      REAL*8 DLA 
      REAL*8 DUA 
      REAL*8 DRhoB 
      REAL*8 DGBB
      REAL*8 DTB 
      REAL*8 DLB 
      REAL*8 DUB 
      REAL*8 DGAB

      ! constant
      REAL*8 ZERO
      REAL*8 ONE
      REAL*8 F12
      INTEGER ISPIN
      INTEGER NE
C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! constants
      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0
      NE   = NA+NB

      ! loop over NG
      DO I = 1,NG

C         write(6,*)"     "
C         write(6,*)"----------------------"
C         write(6,*)"for grid: ", i
          
         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)

         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! debug
C         write(6,*)"our RhoA: ",RA
C         write(6,*)"our GRhoA: ",GAA 
C         write(6,*)"our TA: ",TA
C         write(6,*)"our LA: ",LA 
C         write(6,*)"our UA: ",UA
C         write(6,*)"our RhoB: ",RB
C         write(6,*)"our GRhoB: ",GBB 
C         write(6,*)"our TB: ",TB
C         write(6,*)"our LB: ",LB 
C         write(6,*)"our UB: ",UB 

         ! now derive the alpha part first
         iSpin  = 1
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NA,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRA  = ZERO
         DM1DGAA = ZERO
         DM1DTA  = ZERO
         DM1DLA  = ZERO
         DM1DUA  = ZERO
         CALL becke05_M1(THRESH,RA,GAA,TA,LA,UA,HIRWT,M1,
     $DM1DRA,DM1DGAA,DM1DTA,DM1DLA,DM1DUA)

         ! now compose the finaly energy in terms of alpha
         F(i) = F(i) - F12*RA*A*M1

         ! derivatives
         DRhoA= -F12*A*M1-F12*RA*DADRA*M1-F12*RA*A*DM1DRA
         DGAA = -F12*RA*DADGAA*M1-F12*RA*A*DM1DGAA
         DTA  = -F12*RA*DADTA*M1-F12*RA*A*DM1DTA
         DLA  = -F12*RA*DADLA*M1-F12*RA*A*DM1DLA
         DUA  = -F12*RA*DADUA*M1-F12*RA*A*DM1DUA
         DRhoB= -F12*RA*DADRB*M1
         DGBB = -F12*RA*DADGBB*M1
         DTB  = -F12*RA*DADTB*M1
         DLB  = -F12*RA*DADLB*M1
         DUB  = -F12*RA*DADUB*M1

         ! collect all of terms for derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + DRhoA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + DGAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = ZERO
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + DTA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + DLA
         ID_EXA_POS=D1VARS(ID_EXA)
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + DUA
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + DRhoB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + DGBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + DTB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + DLB
         ID_EXB_POS=D1VARS(ID_EXB)
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + DUB

         ! now let's do beta part
         iSpin  = 2
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NB,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRB  = ZERO
         DM1DGBB = ZERO
         DM1DTB  = ZERO
         DM1DLB  = ZERO
         DM1DUB  = ZERO
         CALL becke05_M1(THRESH,RB,GBB,TB,LB,UB,HIRWT,M1,
     $DM1DRB,DM1DGBB,DM1DTB,DM1DLB,DM1DUB)

         ! now compose the finaly energy in terms of alpha
         F(i) = F(i) - F12*RB*A*M1

         ! derivatives
         DRhoA= -F12*RB*DADRA*M1
         DGAA = -F12*RB*DADGAA*M1
         DTA  = -F12*RB*DADTA*M1
         DLA  = -F12*RB*DADLA*M1
         DUA  = -F12*RB*DADUA*M1
         DRhoB= -F12*A*M1-F12*RB*DADRB*M1-F12*RB*A*DM1DRB
         DGBB = -F12*RB*DADGBB*M1-F12*RB*A*DM1DGBB
         DTB  = -F12*RB*DADTB*M1-F12*RB*A*DM1DTB
         DLB  = -F12*RB*DADLB*M1-F12*RB*A*DM1DLB
         DUB  = -F12*RB*DADUB*M1-F12*RB*A*DM1DUB

         ! collect all of terms for derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + DRhoA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + DGAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = ZERO
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + DTA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + DLA
         ID_EXA_POS=D1VARS(ID_EXA)
         D1F(i, ID_EXA_POS) = D1F(i, ID_EXA_POS) + DUA
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + DRhoB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + DGBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + DTB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + DLB
         ID_EXB_POS=D1VARS(ID_EXB)
         D1F(i, ID_EXB_POS) = D1F(i, ID_EXB_POS) + DUB
      END DO
      RETURN
      END

      SUBROUTINE Becke05_NDPAR2(INFOR,NDEN,NG,NA,NB,THRESH,c_ndpar_cap,
     $                          VAL_P,VAL_Q,HIRWTS,RhoA,
     $                          RhoB,DRA,DRB,TauA,TauB,
     $                          LapA,LapB,EXA,EXB,F,TMP_D1F)
c
c    ******************************************************************
c    *                                                                *
c    *                                                                *
c     This is the same as _NDPAR except that the parallel contribution
c     is scaled per spin. It did not work out as well as _NDPAR1 where
c     the scaling was on the total. So it is not used.
c
c
c    *  evaluates the non-dynamic parallel spin part of Becke05       *
c    *  exchange functional with capped from above NDpar              *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    * Jing Kong and E. Proynov                                       *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I,J
      INTEGER NA,NB,iSpin
      INTEGER INFOR(*)
      REAL*8 THRESH,c_ndpar_cap
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 HIRWT
      REAL*8 VAL_P
      REAL*8 VAL_Q
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! tmp variables 
      ! the A factor and it's derivatives
      REAL*8 TMP_D1F(*)
      REAL*8 TMP2_D1F(N_FUNC_DERIV_1)
      REAL*8 usopp,usopp_RA,usopp_RB,usopp_GAA,usopp_GBB
      REAL*8 usopp_TA,usopp_TB,usopp_LA,usopp_LB,usopp_UA,usopp_UB
      REAL*8 uspar,uspar_RA,uspar_RB,uspar_GAA,uspar_GBB,uspar2
      REAL*8 uspar_TA,uspar_TB,uspar_LA,uspar_LB,uspar_UA,uspar_UB
      REAL*8 DF1DRA,DF2DRA,DF1DRB,DF2DRB,DF1DGA,DF2DGA,DF1DGB,DF2DGB
      REAL*8 DF1DTA,DF2DTA,DF1DTB,DF2DTB,DF1DLA,DF2DLA,DF1DLB,DF2DLB
      REAL*8 DF1DUA,DF2DUA,DF1DUB,DF2DUB,HPZ,DHPZDZ,DZDRA,DZDRB
      REAL*8 DZDGA,DZDGB,DZDTA,DZDTB,DZDLA,DZDLB,DZDUA,DZDUB
      REAL*8 P,Z,PZ,FF,F1,F2,T,T1,T2,D1R1,D1R2,D1G1,D1G2,D1T1,D1T2
      REAL*8 D1L1,D1L2,D1U1,D1U2,TP,DZ1,DZ2,LIMIT_Z

      ! constant
      REAL*8 ZERO
      REAL*8 ONE
      REAL*8 TWO 
      REAL*8 F12
      REAL*8 dtol
c
      ZERO = 0.0D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      F12  = c_ndpar_cap
      P    = 250.0D0
      FF = ZERO
      F1= ZERO
      F2= ZERO
      F = ZERO
c
        CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
c     DO i = 1, NG 
          usopp = ZERO
          DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
          END DO

        call Becke05_NDOP(INFOR,NDEN,1,THRESH,VAL_P,HIRWTS,
     $               RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,
     $               EXA,EXB,usopp,TMP_D1F)
c
c        write(6,*)"usopp =   ", usopp
         ID_RA_POS  = D1VARS(ID_RA)
         ID_GAA_POS = D1VARS(ID_GAA)
         ID_GAB_POS = D1VARS(ID_GAB)
         ID_TA_POS  = D1VARS(ID_TA)
         ID_LA_POS  = D1VARS(ID_LA)
         ID_EXA_POS = D1VARS(ID_EXA)
         ID_RB_POS  = D1VARS(ID_RB)
         ID_GBB_POS = D1VARS(ID_GBB)
         ID_TB_POS  = D1VARS(ID_TB)
         ID_LB_POS  = D1VARS(ID_LB)
         ID_EXB_POS = D1VARS(ID_EXB)
c
         usopp_RA   = TMP_D1F(ID_RA_POS)
         usopp_GAA  = TMP_D1F(ID_GAA_POS)
         usopp_TA   = TMP_D1F(ID_TA_POS)
         usopp_LA   = TMP_D1F(ID_LA_POS)
         usopp_UA   = TMP_D1F(ID_EXA_POS)
         usopp_RB   = TMP_D1F(ID_RB_POS)
         usopp_GBB  = TMP_D1F(ID_GBB_POS)
         usopp_TB   = TMP_D1F(ID_TB_POS)
         usopp_LB   = TMP_D1F(ID_LB_POS)
         usopp_UB   = TMP_D1F(ID_EXB_POS)
c
         iSpin = 1
         uspar=ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP_D1F(J) = ZERO
         END DO
c
       call  Becke05_NDPARspin(INFOR,NDEN,1,NA,NB,THRESH,
     $                     VAL_P,VAL_Q,HIRWTS,
     $                     RhoA,RhoB,DRA,DRB,TauA,TauB,
     $                     LapA,LapB,EXA,EXB,uspar,TMP_D1F,1)
c
c        write(6,*)"uspar_a =   ", uspar
         uspar_RA   = TMP_D1F(ID_RA_POS)
         uspar_GAA  = TMP_D1F(ID_GAA_POS)
         uspar_TA   = TMP_D1F(ID_TA_POS)
         uspar_LA   = TMP_D1F(ID_LA_POS)
         uspar_UA   = TMP_D1F(ID_EXA_POS)
         uspar_RB   = TMP_D1F(ID_RB_POS)
         uspar_GBB  = TMP_D1F(ID_GBB_POS)
         uspar_TB   = TMP_D1F(ID_TB_POS)
         uspar_LB   = TMP_D1F(ID_LB_POS)
         uspar_UB   = TMP_D1F(ID_EXB_POS)
       ! let;s deal with Z first
         Z = ZERO
         F1= uspar
         F2 = c_ndpar_cap*usopp
         IF(usopp > ZERO) F2 = ZERO 
         T = F1*F1+F2*F2
c        write(6,*)"T =   ", T
c
       DF1DRA = uspar_RA 
       DF1DGA = uspar_GAA
       DF1DTA = uspar_TA 
       DF1DLA = uspar_LA 
       DF1DUA = uspar_UA 
       DF1DRB = uspar_RB 
       DF1DGB = uspar_GBB
       DF1DTB = uspar_TB 
       DF1DLB = uspar_LB 
       DF1DUB = uspar_UB 
       DF2DRA = F12*usopp_RA 
       DF2DGA = F12*usopp_GAA
       DF2DTA = F12*usopp_TA 
       DF2DLA = F12*usopp_LA 
       DF2DUA = F12*usopp_UA 
       DF2DRB = F12*usopp_RB 
       DF2DGB = F12*usopp_GBB
       DF2DTB = F12*usopp_TB 
       DF2DLB = F12*usopp_LB 
       DF2DUB = F12*usopp_UB 
c
       DZDRA = ZERO 
       DZDGA = ZERO 
       DZDTA = ZERO 
       DZDLA = ZERO 
       DZDUA = ZERO 
       DZDRB = ZERO 
       DZDGB = ZERO 
       DZDTB = ZERO 
       DZDLB = ZERO 
       DZDUB = ZERO 
       D1R1 = ZERO
       D1R2 = ZERO
       D1G1 = ZERO
       D1G2 = ZERO
       D1T1 = ZERO
       D1T2 = ZERO
       D1L1 = ZERO
       D1L2 = ZERO
       D1U1 = ZERO
       D1U2 = ZERO
       DZ1=ZERO
       DZ2=ZERO
       FF = ZERO
       IF (T > THRESH) THEN
         T1= DSQRT(T)
         Z = (F2-F1)/T1
         T2 = T1**3
       IF (T2 > THRESH) THEN
         DZ1= -F2*(F1+F2)/T2
         DZ2= F1*(F1+F2)/T2 
         DZDRA=DZ1*DF1DRA + DZ2*DF2DRA
         DZDGA=DZ1*DF1DGA + DZ2*DF2DGA
         DZDTA=DZ1*DF1DTA + DZ2*DF2DTA
         DZDLA=DZ1*DF1DLA + DZ2*DF2DLA
         DZDUA=DZ1*DF1DUA + DZ2*DF2DUA
         DZDRB=DZ1*DF1DRB + DZ2*DF2DRB
         DZDGB=DZ1*DF1DGB + DZ2*DF2DGB
         DZDTB=DZ1*DF1DTB + DZ2*DF2DTB
         DZDLB=DZ1*DF1DLB + DZ2*DF2DLB
         DZDUB=DZ1*DF1DUB + DZ2*DF2DUB
       END IF
       END IF

       ! finally, we apply the smooth function to get the final f
       ! the original f is the min(f1,f2)
       ! one thing to note, that PZ below could be very large then
       ! it may break the expression of HPZ
       ! therefore we use log to test it
       LIMIT_Z = DLOG(ONE/THRESH)
c      LIMIT_Z = 99.0d0
       HPZ   = ZERO
       DHPZDZ= ZERO
       PZ    = P*Z
       TP=ZERO
       IF (DABS(PZ) < LIMIT_Z) THEN
          TP   = DEXP(P*Z)
          HPZ = ONE/(ONE+TP)
          DHPZDZ = -P*HPZ*HPZ*TP
       ELSE IF (PZ < ZERO) THEN
          HPZ   = ONE
       ELSE IF (PZ > ZERO) THEN
          HPZ   = ZERO
       END IF
c
        FF  = (F1-F2)*HPZ + F2
c      if (.not.(abs(FF-uspar) < 1D-5 .or. abs(FF-F12*usopp) < 1D-5)) 
c    $ then
c          write(*,*) "catching this", FF, "  ", uspar, "  ", usopp
c      endif 

       ! now let's collect all of terms together
C         write(6,*)"uspar ",uspar, "
       D1R1 = (DF1DRA-DF2DRA)*HPZ+(F1-F2)*DHPZDZ*DZDRA+DF2DRA
       D1R2 = (DF1DRB-DF2DRB)*HPZ+(F1-F2)*DHPZDZ*DZDRB+DF2DRB
       D1G1 = (DF1DGA-DF2DGA)*HPZ+(F1-F2)*DHPZDZ*DZDGA+DF2DGA
       D1G2 = (DF1DGB-DF2DGB)*HPZ+(F1-F2)*DHPZDZ*DZDGB+DF2DGB
       D1T1 = (DF1DTA-DF2DTA)*HPZ+(F1-F2)*DHPZDZ*DZDTA+DF2DTA
       D1T2 = (DF1DTB-DF2DTB)*HPZ+(F1-F2)*DHPZDZ*DZDTB+DF2DTB
       D1L1 = (DF1DLA-DF2DLA)*HPZ+(F1-F2)*DHPZDZ*DZDLA+DF2DLA
       D1L2 = (DF1DLB-DF2DLB)*HPZ+(F1-F2)*DHPZDZ*DZDLB+DF2DLB
       D1U1 = (DF1DUA-DF2DUA)*HPZ+(F1-F2)*DHPZDZ*DZDUA+DF2DUA
       D1U2 = (DF1DUB-DF2DUB)*HPZ+(F1-F2)*DHPZDZ*DZDUB+DF2DUB
c
c      if (.not.(abs(D1R1-DF1DRA)< 1D-5.or.abs(D1R1-DF2DRA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1R1," ", DF1DRA," ",DF2DRA 
c      endif
c
c      if (.not.(abs(D1G1-DF1DGA)< 1D-5.or.abs(D1G1-DF2DGA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1G1," ", DF1DGA," ",DF2DGA
c      endif
c
c      if (.not.(abs(D1T1-DF1DGA)< 1D-5.or.abs(D1T1-DF2DTA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1T1," ", DF1DTA," ",DF2DTA
c      endif
c
c      if (.not.(abs(D1L1-DF1DLA)< 1D-5.or.abs(D1L1-DF2DLA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1L1," ", DF1DLA," ",DF2DLA
c      endif
c
c      if (.not.(abs(D1U1-DF1DUA)< 1D-5.or.abs(D1U1-DF2DUA) < 1D-5))
c    $ then
c          write(*,*) "catching this", D1U1," ", DF1DUA," ",DF2DUA
c      endif
c
c        CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
 
c        ID_RA_POS  = D1VARS(ID_RA)
c        ID_GAA_POS = D1VARS(ID_GAA)
c        ID_GAB_POS = D1VARS(ID_GAB)
c        ID_TA_POS  = D1VARS(ID_TA)
c        ID_LA_POS  = D1VARS(ID_LA)
c        ID_EXA_POS = D1VARS(ID_EXA)
c        ID_RB_POS  = D1VARS(ID_RB)
c        ID_GBB_POS = D1VARS(ID_GBB)
c        ID_TB_POS  = D1VARS(ID_TB)
c        ID_LB_POS  = D1VARS(ID_LB)
c        ID_EXB_POS = D1VARS(ID_EXB)
c
          F= F+FF
         TMP_D1F(ID_RA_POS)  = TMP_D1F(ID_RA_POS) + D1R1 
         TMP_D1F(ID_GAA_POS) = TMP_D1F(ID_GAA_POS) + D1G1 
         TMP_D1F(ID_GAB_POS) = ZERO
         TMP_D1F(ID_TA_POS)  = TMP_D1F(ID_TA_POS) + D1T1
         TMP_D1F(ID_LA_POS)  = TMP_D1F(ID_LA_POS) + D1L1 
         TMP_D1F(ID_EXA_POS) = TMP_D1F(ID_EXA_POS) + D1U1 
         TMP_D1F(ID_RB_POS)  = TMP_D1F(ID_RB_POS) +  D1R2 
         TMP_D1F(ID_GBB_POS) = TMP_D1F(ID_GBB_POS) + D1G2 
         TMP_D1F(ID_TB_POS)  = TMP_D1F(ID_TB_POS) + D1T2 
         TMP_D1F(ID_LB_POS)  = TMP_D1F(ID_LB_POS) + D1L2 
         TMP_D1F(ID_EXB_POS) = TMP_D1F(ID_EXB_POS) + D1U2 
c         
         iSpin = 2
c       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
c        ID_RA_POS  = D1VARS(ID_RA)
c        ID_GAA_POS = D1VARS(ID_GAA)
c        ID_GAB_POS = D1VARS(ID_GAB)
c        ID_TA_POS  = D1VARS(ID_TA)
c        ID_LA_POS  = D1VARS(ID_LA)
c        ID_EXA_POS = D1VARS(ID_EXA)
c        ID_RB_POS  = D1VARS(ID_RB)
c        ID_GBB_POS = D1VARS(ID_GBB)
c        ID_TB_POS  = D1VARS(ID_TB)
c        ID_LB_POS  = D1VARS(ID_LB)
c        ID_EXB_POS = D1VARS(ID_EXB)
c
         uspar2=ZERO
         FF = ZERO 
         F1= ZERO
         DO J = 1, N_FUNC_DERIV_1
            TMP2_D1F(J) = ZERO
         END DO
c
          uspar_RA   =  ZERO
          uspar_GAA  =  ZERO
          uspar_TA   =  ZERO
          uspar_LA   =  ZERO
          uspar_UA   =  ZERO
          uspar_RB   =  ZERO
          uspar_GBB  =  ZERO
          uspar_TB   =  ZERO
          uspar_LB   =  ZERO
          uspar_UB   =  ZERO
       call  Becke05_NDPARspin(INFOR,NDEN,1,NA,NB,THRESH,
     $                     VAL_P,VAL_Q,HIRWTS,
     $                     RhoA,RhoB,DRA,DRB,TauA,TauB,
     $                     LapA,LapB,EXA,EXB,uspar2,TMP2_D1F,2)
c
         uspar_RA   = TMP2_D1F(ID_RA_POS)
         uspar_GAA  = TMP2_D1F(ID_GAA_POS)
         uspar_TA   = TMP2_D1F(ID_TA_POS)
         uspar_LA   = TMP2_D1F(ID_LA_POS)
         uspar_UA   = TMP2_D1F(ID_EXA_POS)
         uspar_RB   = TMP2_D1F(ID_RB_POS)
         uspar_GBB  = TMP2_D1F(ID_GBB_POS)
         uspar_TB   = TMP2_D1F(ID_TB_POS)
         uspar_LB   = TMP2_D1F(ID_LB_POS)
         uspar_UB   = TMP2_D1F(ID_EXB_POS)
       ! let;s deal with Z first
         Z = ZERO
         F1= uspar2
         T = F1*F1+F2*F2
c        write(6,*)"T =   ", T
c
       DF1DRA = uspar_RA 
       DF1DGA = uspar_GAA
       DF1DTA = uspar_TA 
       DF1DLA = uspar_LA 
       DF1DUA = uspar_UA 
       DF1DRB = uspar_RB 
       DF1DGB = uspar_GBB
       DF1DTB = uspar_TB 
       DF1DLB = uspar_LB 
       DF1DUB = uspar_UB 
c
       DZDRA = ZERO 
       DZDGA = ZERO 
       DZDTA = ZERO 
       DZDLA = ZERO 
       DZDUA = ZERO 
       DZDRB = ZERO 
       DZDGB = ZERO 
       DZDTB = ZERO 
       DZDLB = ZERO 
       DZDUB = ZERO 
       D1R1 = ZERO
       D1R2 = ZERO
       D1G1 = ZERO
       D1G2 = ZERO
       D1T1 = ZERO
       D1T2 = ZERO
       D1L1 = ZERO
       D1L2 = ZERO
       D1U1 = ZERO
       D1U2 = ZERO
       FF = ZERO
       DZ1=ZERO
       DZ2=ZERO
       IF (T > THRESH) THEN
         T1= DSQRT(T)
         Z = (F2-F1)/T1
         T2 = T1**3
       IF (T2 > THRESH) THEN
         DZ1= -F2*(F1+F2)/T2
         DZ2= F1*(F1+F2)/T2 
         DZDRA=DZ1*DF1DRA + DZ2*DF2DRA
         DZDGA=DZ1*DF1DGA + DZ2*DF2DGA
         DZDTA=DZ1*DF1DTA + DZ2*DF2DTA
         DZDLA=DZ1*DF1DLA + DZ2*DF2DLA
         DZDUA=DZ1*DF1DUA + DZ2*DF2DUA
         DZDRB=DZ1*DF1DRB + DZ2*DF2DRB
         DZDGB=DZ1*DF1DGB + DZ2*DF2DGB
         DZDTB=DZ1*DF1DTB + DZ2*DF2DTB
         DZDLB=DZ1*DF1DLB + DZ2*DF2DLB
         DZDUB=DZ1*DF1DUB + DZ2*DF2DUB
       END IF
       END IF

       ! finally, we apply the smooth function to get the final f
       ! the original f is the min(f1,f2)
       ! one thing to note, that PZ below could be very large then
       ! it may break the expression of HPZ
       ! therefore we use log to test it
       LIMIT_Z = DLOG(ONE/THRESH)
c      LIMIT_Z = 99.0d0
       HPZ   = ZERO
       DHPZDZ= ZERO
       PZ    = P*Z
       TP=ZERO
       IF (DABS(PZ) < LIMIT_Z) THEN
          TP   = DEXP(P*Z)
          HPZ = ONE/(ONE+TP)
          DHPZDZ = -P*HPZ*HPZ*TP
       ELSE IF (PZ < ZERO) THEN
          HPZ   = ONE
       ELSE IF (PZ > ZERO) THEN
          HPZ   = ZERO
       END IF

       FF   = (F1-F2)*HPZ + F2
c      if (.not.(abs(FF-uspar) < 1D-5 .or. abs(FF-F12*usopp) < 1D-5)) 
c    $ then
c          write(*,*) "catching this", FF, "  ", uspar, "  ", usopp
c      endif 

       ! now let's collect all of terms together
C         write(6,*)"uspar ",uspar, "
       D1R1 = (DF1DRA-DF2DRA)*HPZ+(F1-F2)*DHPZDZ*DZDRA+DF2DRA
       D1R2 = (DF1DRB-DF2DRB)*HPZ+(F1-F2)*DHPZDZ*DZDRB+DF2DRB
       D1G1 = (DF1DGA-DF2DGA)*HPZ+(F1-F2)*DHPZDZ*DZDGA+DF2DGA
       D1G2 = (DF1DGB-DF2DGB)*HPZ+(F1-F2)*DHPZDZ*DZDGB+DF2DGB
       D1T1 = (DF1DTA-DF2DTA)*HPZ+(F1-F2)*DHPZDZ*DZDTA+DF2DTA
       D1T2 = (DF1DTB-DF2DTB)*HPZ+(F1-F2)*DHPZDZ*DZDTB+DF2DTB
       D1L1 = (DF1DLA-DF2DLA)*HPZ+(F1-F2)*DHPZDZ*DZDLA+DF2DLA
       D1L2 = (DF1DLB-DF2DLB)*HPZ+(F1-F2)*DHPZDZ*DZDLB+DF2DLB
       D1U1 = (DF1DUA-DF2DUA)*HPZ+(F1-F2)*DHPZDZ*DZDUA+DF2DUA
       D1U2 = (DF1DUB-DF2DUB)*HPZ+(F1-F2)*DHPZDZ*DZDUB+DF2DUB
c
c       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
 
c        ID_RA_POS  = D1VARS(ID_RA)
c        ID_GAA_POS = D1VARS(ID_GAA)
c        ID_GAB_POS = D1VARS(ID_GAB)
c        ID_TA_POS  = D1VARS(ID_TA)
c        ID_LA_POS  = D1VARS(ID_LA)
c        ID_EXA_POS = D1VARS(ID_EXA)
c        ID_RB_POS  = D1VARS(ID_RB)
c        ID_GBB_POS = D1VARS(ID_GBB)
c        ID_TB_POS  = D1VARS(ID_TB)
c        ID_LB_POS  = D1VARS(ID_LB)
c        ID_EXB_POS = D1VARS(ID_EXB)
c
         ID_RA_POS  = D1VARS(ID_RA)
         ID_GAA_POS = D1VARS(ID_GAA)
         ID_GAB_POS = D1VARS(ID_GAB)
         ID_TA_POS  = D1VARS(ID_TA)
         ID_LA_POS  = D1VARS(ID_LA)
         ID_EXA_POS = D1VARS(ID_EXA)
         ID_RB_POS  = D1VARS(ID_RB)
         ID_GBB_POS = D1VARS(ID_GBB)
         ID_TB_POS  = D1VARS(ID_TB)
         ID_LB_POS  = D1VARS(ID_LB)
         ID_EXB_POS = D1VARS(ID_EXB)
c
          F= F + FF
         TMP_D1F(ID_RA_POS)  = TMP_D1F(ID_RA_POS) + D1R1 
         TMP_D1F(ID_GAA_POS) = TMP_D1F(ID_GAA_POS) + D1G1 
         TMP_D1F(ID_GAB_POS) = ZERO
         TMP_D1F(ID_TA_POS)  = TMP_D1F(ID_TA_POS) + D1T1
         TMP_D1F(ID_LA_POS)  = TMP_D1F(ID_LA_POS) + D1L1 
         TMP_D1F(ID_EXA_POS) = TMP_D1F(ID_EXA_POS) + D1U1 
         TMP_D1F(ID_RB_POS)  = TMP_D1F(ID_RB_POS) +  D1R2 
         TMP_D1F(ID_GBB_POS) = TMP_D1F(ID_GBB_POS) + D1G2 
         TMP_D1F(ID_TB_POS)  = TMP_D1F(ID_TB_POS) + D1T2 
         TMP_D1F(ID_LB_POS)  = TMP_D1F(ID_LB_POS) + D1L2 
         TMP_D1F(ID_EXB_POS) = TMP_D1F(ID_EXB_POS) + D1U2 
c      END DO
       RETURN
       END

      SUBROUTINE Becke05_NDPARspin(INFOR,NDEN,NG,NA,NB,THRESH,
     $VAL_P,VAL_Q,HIRWTS,
     $RhoA,RhoB,DRA,DRB,TauA,TauB,LapA,LapB,EXA,EXB,F,TMP_D1F,iSpin)
c
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic parallel spin part of Becke05       *
c    *  exchange functional                                           *
c    *  the inplementation is based on the paper of:                  *
c    *  E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.              * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *  the original idea about Becke05 is coming from the paper:     *     
c    *  A. D. Becke, J. Chem. Phys. 122, 64101 (2005) etc.            *
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *   EXA,EXB:    exact exchange energy density                    *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB, iSpin
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 EXA(NG),EXB(NG)
      REAL*8 HIRWTS(NG)
      REAL*8 F
      REAL*8 TMP_D1F(*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
      REAL*8 UA,UB
      REAL*8 HIRWT
      REAL*8 VAL_P
      REAL*8 VAL_Q
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS
      INTEGER  ID_EXA_POS, ID_EXB_POS

      ! tmp variables 
      ! the A factor and it's derivatives
      REAL*8 A   
      REAL*8 DADRA 
      REAL*8 DADGAA
      REAL*8 DADTA 
      REAL*8 DADLA 
      REAL*8 DADUA 
      REAL*8 DADRB 
      REAL*8 DADGBB
      REAL*8 DADTB 
      REAL*8 DADLB 
      REAL*8 DADUB 

      ! tmp variables 
      ! the M1 and it's derivatives
      REAL*8 M1   
      REAL*8 DM1DRA 
      REAL*8 DM1DGAA
      REAL*8 DM1DTA 
      REAL*8 DM1DLA 
      REAL*8 DM1DUA 
      REAL*8 DM1DRB 
      REAL*8 DM1DGBB
      REAL*8 DM1DTB 
      REAL*8 DM1DLB 
      REAL*8 DM1DUB 

      ! now for the real derivatives
      REAL*8 DRhoA 
      REAL*8 DGAA
      REAL*8 DTA 
      REAL*8 DLA 
      REAL*8 DUA 
      REAL*8 DRhoB 
      REAL*8 DGBB
      REAL*8 DTB 
      REAL*8 DLB 
      REAL*8 DUB 
      REAL*8 DGAB

      ! constant
      REAL*8 ZERO
      REAL*8 ONE
      REAL*8 F12
c     INTEGER ISPIN
      INTEGER NE
C       
C     the tau used in this program is for "big tau", that 
C     is to say, without factor of 1/2
C

      ! firstly initilize variable position information
      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
      ! constants
      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0
      NE   = NA+NB

      ! loop over NG
      DO I = 1,NG

C         write(6,*)"     "
C         write(6,*)"----------------------"
C         write(6,*)"for grid: ", i
          
         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
         UA  = EXA(i)
         UB  = EXB(i)

         ! possible hirshfeld weights
         HIRWT = HIRWTS(i)

         ! debug
C         write(6,*)"our RhoA: ",RA
C         write(6,*)"our GRhoA: ",GAA 
C         write(6,*)"our TA: ",TA
C         write(6,*)"our LA: ",LA 
C         write(6,*)"our UA: ",UA
C         write(6,*)"our RhoB: ",RB
C         write(6,*)"our GRhoB: ",GBB 
C         write(6,*)"our TB: ",TB
C         write(6,*)"our LB: ",LB 
C         write(6,*)"our UB: ",UB 

         ! now derive the alpha part first
C        iSpin  = 1
         if ( iSpin .eq. 1 ) then
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NA,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRA  = ZERO
         DM1DGAA = ZERO
         DM1DTA  = ZERO
         DM1DLA  = ZERO
         DM1DUA  = ZERO
         CALL becke05_M1(THRESH,RA,GAA,TA,LA,UA,HIRWT,M1,
     $DM1DRA,DM1DGAA,DM1DTA,DM1DLA,DM1DUA)

         ! now compose the finaly energy in terms of alpha
         F = F - F12*RA*A*M1

         ! derivatives
         DRhoA= -F12*A*M1-F12*RA*DADRA*M1-F12*RA*A*DM1DRA
         DGAA = -F12*RA*DADGAA*M1-F12*RA*A*DM1DGAA
         DTA  = -F12*RA*DADTA*M1-F12*RA*A*DM1DTA
         DLA  = -F12*RA*DADLA*M1-F12*RA*A*DM1DLA
         DUA  = -F12*RA*DADUA*M1-F12*RA*A*DM1DUA
         DRhoB= -F12*RA*DADRB*M1
         DGBB = -F12*RA*DADGBB*M1
         DTB  = -F12*RA*DADTB*M1
         DLB  = -F12*RA*DADLB*M1
         DUB  = -F12*RA*DADUB*M1

         ! collect all of terms for derivatives
         ID_RA_POS=D1VARS(ID_RA)
         TMP_D1F(ID_RA_POS)  = TMP_D1F(ID_RA_POS) + DRhoA
         ID_GAA_POS=D1VARS(ID_GAA)
         TMP_D1F(ID_GAA_POS) = TMP_D1F(ID_GAA_POS) + DGAA
         ID_GAB_POS=D1VARS(ID_GAB)
         TMP_D1F(ID_GAB_POS) = ZERO
         ID_TA_POS=D1VARS(ID_TA)
         TMP_D1F(ID_TA_POS)  = TMP_D1F(ID_TA_POS) + DTA
         ID_LA_POS=D1VARS(ID_LA)
         TMP_D1F(ID_LA_POS)  = TMP_D1F(ID_LA_POS) + DLA
         ID_EXA_POS=D1VARS(ID_EXA)
         TMP_D1F(ID_EXA_POS) = TMP_D1F(ID_EXA_POS) + DUA
         ID_RB_POS=D1VARS(ID_RB)
         TMP_D1F(ID_RB_POS)  = TMP_D1F(ID_RB_POS) + DRhoB
         ID_GBB_POS=D1VARS(ID_GBB)
         TMP_D1F(ID_GBB_POS) = TMP_D1F(ID_GBB_POS) + DGBB
         ID_TB_POS=D1VARS(ID_TB)
         TMP_D1F(ID_TB_POS)  = TMP_D1F(ID_TB_POS) + DTB
         ID_LB_POS=D1VARS(ID_LB)
         TMP_D1F(ID_LB_POS)  = TMP_D1F(ID_LB_POS) + DLB
         ID_EXB_POS=D1VARS(ID_EXB)
         TMP_D1F(ID_EXB_POS) = TMP_D1F(ID_EXB_POS) + DUB

         ! now let's do beta part
c         iSpin  = 2
         else
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NB,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRB  = ZERO
         DM1DGBB = ZERO
         DM1DTB  = ZERO
         DM1DLB  = ZERO
         DM1DUB  = ZERO
         CALL becke05_M1(THRESH,RB,GBB,TB,LB,UB,HIRWT,M1,
     $DM1DRB,DM1DGBB,DM1DTB,DM1DLB,DM1DUB)

         ! now compose the finaly energy in terms of alpha
         F = F - F12*RB*A*M1

         ! derivatives
         DRhoA= -F12*RB*DADRA*M1
         DGAA = -F12*RB*DADGAA*M1
         DTA  = -F12*RB*DADTA*M1
         DLA  = -F12*RB*DADLA*M1
         DUA  = -F12*RB*DADUA*M1
         DRhoB= -F12*A*M1-F12*RB*DADRB*M1-F12*RB*A*DM1DRB
         DGBB = -F12*RB*DADGBB*M1-F12*RB*A*DM1DGBB
         DTB  = -F12*RB*DADTB*M1-F12*RB*A*DM1DTB
         DLB  = -F12*RB*DADLB*M1-F12*RB*A*DM1DLB
         DUB  = -F12*RB*DADUB*M1-F12*RB*A*DM1DUB

         ! collect all of terms for derivatives
         ID_RA_POS=D1VARS(ID_RA)
         TMP_D1F(ID_RA_POS)  = TMP_D1F(ID_RA_POS) + DRhoA
         ID_GAA_POS=D1VARS(ID_GAA)
         TMP_D1F(ID_GAA_POS) = TMP_D1F(ID_GAA_POS) + DGAA
         ID_GAB_POS=D1VARS(ID_GAB)
         TMP_D1F(ID_GAB_POS) = ZERO
         ID_TA_POS=D1VARS(ID_TA)
         TMP_D1F(ID_TA_POS)  = TMP_D1F(ID_TA_POS) + DTA
         ID_LA_POS=D1VARS(ID_LA)
         TMP_D1F(ID_LA_POS)  = TMP_D1F(ID_LA_POS) + DLA
         ID_EXA_POS=D1VARS(ID_EXA)
         TMP_D1F(ID_EXA_POS) = TMP_D1F(ID_EXA_POS) + DUA
         ID_RB_POS=D1VARS(ID_RB)
         TMP_D1F(ID_RB_POS)  = TMP_D1F(ID_RB_POS) + DRhoB
         ID_GBB_POS=D1VARS(ID_GBB)
         TMP_D1F(ID_GBB_POS) = TMP_D1F(ID_GBB_POS) + DGBB
         ID_TB_POS=D1VARS(ID_TB)
         TMP_D1F(ID_TB_POS)  = TMP_D1F(ID_TB_POS) + DTB
         ID_LB_POS=D1VARS(ID_LB)
         TMP_D1F(ID_LB_POS)  = TMP_D1F(ID_LB_POS) + DLB
         ID_EXB_POS=D1VARS(ID_EXB)
         TMP_D1F(ID_EXB_POS) = TMP_D1F(ID_EXB_POS) + DUB
      endif
      END DO
      RETURN
      END

