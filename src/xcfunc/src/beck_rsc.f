*  Deck beck_rsc
      SUBROUTINE beck_rsc(F,D1F,D1FEXCH,W,RA,RB,D1RA,D1RB,RLapA,RLapB,
     $            TA,TB,EX_HF_DENSA,EX_HF_DENSB,ACOP,ACPAR,
     $ODEL,ODELW,ODDSUM,NGrid,iterSCF,ICEN,NA,NB)
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic correlation energy by the           *
c    *  real-space post-HF model B05 of Becke with RI and fully SCF.  *
c    *                                                                *
c    *  reference: A. D. Becke, J. Chem. Phys. 122, 64101 (2005);     *
c    *  reference: E. Proynov, Y. Shao, and J. Kong, J. Chem. Phys.   * 
c    *             136, 034102 (2012)                                 *
c    *                                                                *
c    *                                                                *
c    *  OUTPUT:                                                       *
c    *     F    - Functional values                                   *
c    *     D1F  - First derivatives                                   *
c    *                                                                *
c    *  INPUT:                                                        *
c    *     RA,B   - Spin densities                                    *
c    *     D1RA,B   - Spin densitiy gradients                         *
c    *     TA,B   - Spin kinetic energy densities                     *
c    *     DLapA,B - Laplacian of the spin densities                  *
c    *     EX_HF_DENSA,B - exact HF-like exchange energy density      *
c    *                                                                *
c    ******************************************************************
c
        IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 F(*),D1F(NGrid,*),RA(NGrid),RB(NGrid),D1RA(NGrid,3),
     $  D1RB(NGrid,3),RLapA(NGrid),RLapB(NGrid),TA(NGrid),TB(NGrid),
     $  D1FEXCH(NGrid,*),W(NGrid),ODEL(NGrid),ODELW(NGrid)
        REAL*8 ODDSUM(*)

        REAL*8 fc1,fc2,UCA,UCB,DSA,DSB,DDA,DDB,ASA,ASB,ASSA,ASSB,sumsum
        REAL*8 EFNA,EFNB,UX1,UX2,tauA,tauB,VX1,VX2,DLapA,DLapB

        REAL*8 EX_HF_DENSA(NGrid)
        REAL*8 EX_HF_DENSB(NGrid)

c        REAL*8 FIT_DENSA(1)
c        REAL*8 FIT_DENSB(1)

       INTEGER NGrid,iterSCF,ICEN,NCENTER
       INTEGER NUMVAR
c
c        POINTER (pEX_HF_DENSA,EX_HF_DENSA)
c       POINTER (pEX_HF_DENSB,EX_HF_DENSB)
c
c        POINTER (pFIT_DENSA,FIT_DENSA)
c        POINTER (pFIT_DENSB,FIT_DENSB)

         PARAMETER(ID_RA = 1)
         PARAMETER(ID_RB = 2)
         PARAMETER(ID_GAA = 3)
         PARAMETER(ID_GAB = 4)
         PARAMETER(ID_GBB = 5)
         PARAMETER(ID_TA = 6)
         PARAMETER(ID_TB = 7)
         PARAMETER(ID_LA = 8)
         PARAMETER(ID_LB = 9)
         PARAMETER(ID_EXA = 10)
         PARAMETER(ID_EXB = 11)


        parameter(third=1.0d0/3.0d0,zero=0.0D0)
        parameter(dsmall=1.0D-15)
        parameter(dsmallb=1.0D-12)
        parameter(dsmall2=1.0D-14)
        parameter(detol=1.0D-08,dtol=1.0D-08)
        parameter(dtol2=1.0D-08)
c       parameter(ACOP=0.5140D0,ACPAR=0.6510D0)  ! Orig
c       parameter(ACOP=0.5150D0,ACPAR=0.6370D0)  ! DF07
        parameter(Tolbig=1.0d+15)
        parameter(delt=0.050d0)  ! smoothening of fa, fb
        parameter(smoth=115.0D0) ! unifies the definition of fcor
        parameter(smoth2=120.0D0)
        SAVE NCENTER
        DATA NCENTER /-1/

        Tol = 1.0D-08
           smth1 = 1.0d0 - delt
           smth2 = 1.0d0 + delt

C       fenglai: so far we only use file to tranfer the D1F for exchange
C       variable, so define the number of variables not including
C       exchange
        NUMVAR = 9           

c        Yihan and Emil, B05 adjustable parameters
c        CALL RemGet(maxiter,REM_MAXSCF)
c        IF(maxiter.lt.2) then
c          ACOP  = 0.5140D0
c          ACPAR = 0.6510D0
c        else

        ! fenglai, here we define the parameter in arbitray one
C        ACOP  = 0.5260D0 
C        ACPAR = 0.64670D0
C        ACOP2 = ACOP
C        ACPAR2 = ACPAR
c        endif
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 0, FM_BEG, ACOP2)
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 1, FM_BEG, ACPAR2)
c        IF (ABS(ACOP2).GT.0.00001) THEN
c          WRITE(*,1002) 'ACOP', ACOP, ACOP2
c1002       FORMAT(1x,'The default value of', A6, 
c    $ 'is changed from', F12.10, ' to', F12.10)
c           ACOP = ACOP2
c        ENDIF
c        IF (ABS(ACPAR2).GT.0.00001) THEN
c          WRITE(*,1002) 'ACPAR', ACPAR, ACPAR2
c           ACPAR = ACPAR2
c        ENDIF
c
c        CALL RemGet(NA,REM_NALPHA)
c        CALL RemGet(NB,REM_NBETA)
c        CALL RemGet(IUnres,REM_JUSTAL)
c        CALL RemGet(NAtoms,REM_NATOMS)
c       write(*,*)'IUnres=  ',IUnres
c
c

           NMO = 0
c        IF((NA.eq.NB).and.(IUnres.ne.1)) then
        IF(NA.eq.NB) then
          NMO = 1
        else
          NMO = 2
        endif
c
c        CALL ftnQAllocREAL8(pEX_HF_DENSA,NGrid)
c        CALL ftnQAllocREAL8(pEX_HF_DENSB,NGrid)
c        CALL VRLoad(EX_HF_DENSA,NGrid,0.0d0)
c        CALL VRLoad(EX_HF_DENSB,NGrid,0.0d0)

c        CALL ftnQAllocREAL8(pFIT_DENSA,NGrid)
c        CALL ftnQAllocREAL8(pFIT_DENSB,NGrid)
c        CALL VRLoad(FIT_DENSA,NGrid,0.0d0)
c        CALL VRLoad(FIT_DENSB,NGrid,0.0d0)

C        CALL RemGet(IMethod, REM_HF_EXCHANGE_WITH_GRID)
C       so far we only support the RI approximation        
        IMethod = 2
c         Using RI-approximated exact exchange

c          CALL FileMan(FM_READ,FILE_EX_HF_DENSITY,FM_DP,NGrid,
c    $                 0,FM_BEG,EX_HF_DENSA)
c          CALL FileMan(FM_READ,FILE_FITTED_DENSITY,FM_DP,NGrid,
c     $                 0,FM_BEG,FIT_DENSA)
C          call prtmat (EX_HF_DENSA, NGrid, 1, NGrid, 1,  6, 
C     $    'EX_HF_DENSA')
C          call prtmat (FIT_DENSA,   NGrid, 1, NGrid, 1,  6, 'FIT_DENSA')

c          IF (NMO.eq.2) THEN
c             CALL FileMan(FM_READ,FILE_EX_HF_DENSITY,FM_DP,NGrid,
c     $             NGrid,FM_BEG,EX_HF_DENSB)
c             CALL FileMan(FM_READ,FILE_FITTED_DENSITY,FM_DP,NGrid,
c     $             NGrid,FM_BEG,FIT_DENSB)
c          ENDIF
c
c         CALL VRLoad(F,NGrid,0.0d0)
c         CALL VRLoad(D1FEXCH,NGrid*2,0.0d0)
c         CALL VRLoad(ODEL,NGrid,0.0d0)
c         CALL VRLoad(ODELW,NGrid,0.0d0)
c         CALL VRLoad(D1F,NGrid*NUMVAR,0.0d0)  ! for the SCF version only
c ATTENTION! to change the allocation of D1F elsewhere too!!
c start initialization:
c
         fc1 = zero
         fc2 = zero
         fcor = zero
         UCA = zero
         UCB = zero
         DSA =dsmall
         DSB =dsmall
         DDA = zero
         DDB = zero
         ASA = zero
         ASSA = zero
         AASA = zero
         ASB = zero
         ASSB = zero
         BBSB = zero
         dASAdra = zero
         dASAdga = zero
         dASAdta = zero
         dASAdla = zero
         dASAdVxa = zero
         dASBdra = zero
         dASBdga = zero
         dASBdta = zero
         dASBdla = zero
         dASBdVxa =zero
         dASBdrb = zero
         dASBdgb = zero
         dASBdtb = zero
         dASBdlb = zero
         dASBdVxb =zero
         dASAdrb = zero
         dASAdgb = zero
         dASAdtb = zero
         dASAdlb = zero
         dASAdVxb =zero
         dASSAdra = zero
         dASSAdga = zero
         dASSAdta = zero
         dASSAdla = zero
         dASSAdVxa =zero
         dASSAdrb = zero
         dASSAdgb = zero
         dASSAdtb = zero
         dASSAdlb = zero
         dASSAdVxb =zero
         dASSBdrb = zero
         dASSBdgb = zero
         dASSBdtb = zero
         dASSBdlb = zero
         dASSBdVxb= zero
         dASSBdra = zero
         dASSBdga = zero
         dASSBdta = zero
         dASSBdla = zero
         dASSBdVxa =zero
         EFNA = 1.0d0
         EFN2A = 1.0d0
         EFNB = 1.0d0
         EFN2B = 1.0d0
         UX1= -zero
         UX2= -zero
         gra=zero
         grb=zero
         tauA =zero
         tauB =zero
         DLapA =zero
         DLapB =zero
         dfdra = zero
         dfdga = zero
         dfdta = zero
         dfdla = zero
         dfdrb = zero
         dfdgb = zero
         dfdtb = zero
         dfdlb = zero
         dfdVxa=zero
         dfdVxb=zero
         dfadra = zero
         dfadga = zero
         dfadta = zero
         dfadla = zero
         dfadrb = zero
         dfadgb = zero
         dfadtb = zero
         dfadlb = zero
         dfadVxa = zero
         dfadVxb = zero
         dfbdra = zero
         dfbdga = zero
         dfbdta = zero
         dfbdla = zero
         dfbdrb = zero
         dfbdgb = zero
         dfbdtb = zero
         dfbdlb = zero
         dfbdVxa = zero
         dfbdVxb = zero
         dfcdraa = zero
         dfcdgaa = zero
         dfcdtaa = zero
         dfcdlaa = zero
         dfcdVxaa =zero
         dfcdrbb = zero
         dfcdgbb = zero
         dfcdtbb = zero
         dfcdlbb = zero
         dfcdVxbb= zero
         dM1dra = zero
         dM1dta = zero
         dM1dla = zero
         dM1dga = zero
         dM1dVxa = zero
         dM2dra = zero
         dM2dta = zero
         dM2dla = zero
         dM2dga = zero
         dM2dVxa = zero
         dM1drb = zero
         dM1dtb = zero
         dM1dlb = zero
         dM1dgb = zero
         dM1dVxb = zero
         dM2drb = zero
         dM2dtb = zero
         dM2dlb = zero
         dM2dgb = zero
         dM2dVxb = zero
         dg3dra = zero
         dg3dta = zero
         dg3dga = zero
         dg3drb = zero
         dg3dtb = zero
         dg3dgb = zero
         dfc1g1 = zero
         dfc2g2 = zero
         VX1= -dsmall
         VX2= -dsmall
         SM1A = zero
         SM2A = zero
         SM1B = zero
         SM2B = zero
         ttgg = zero
         ttg1 = zero
         ttg2 = zero
         ttdB = zero
         ttdA = zero
         dfpardVxa = zero
         dfpardra = zero
         dfpardga = zero
         dfpardla = zero
         dfpardta = zero
         dfpardVxa = zero
         dfpardVxb = zero
         dfpardrb = zero
         dfpardgb = zero
         dfpardlb = zero
         dfpardtb = zero
         dfopdra =zero
         dfopdga =zero
         dfopdla =zero
         dfopdta =zero
         dfopdVxa =zero
         dfopdrb =zero
         dfopdgb =zero
         dfopdlb =zero
         dfopdtb =zero
         dfopdVxb =zero
         gc1 = zero
         gcc1 = zero
         ggg = zero
         ggp = zero
         gc2 = zero
         gcc2 = zero
         dNdra = zero
         dNdga = zero
         dNdta = zero
         dNdla = zero
         dNdVxa = zero
         dNdrb = zero
         dNdgb = zero
         dNdtb = zero
         dNdVxb = zero
         qqx = 0.0d0
         eey1=0.0d0
         eey2=0.0d0
         dNdVxa = zero
         dNdVxb = zero
         qqx = 0.0d0
         smtx = zero
         tx = 0.0d0
         qqfa = 0.0d0
         qqfb = 0.0d0
         xxx = zero
         yy1 = zero
         yy2 = zero
         VCPAR = zero
         Fop = zero
         Fpar = zero
         D1 = 0.0d0
         D2 = 0.0d0
         DF1 = 0.0d0
         DF2 = 0.0d0
         OD1 = 0.0d0
         OD2 = 0.0d0
         ds1 = 0.0d0
c      IF (iterSCF.gt.1) THEN
       IF (iterSCF.gt.(-1)) THEN
c
       DO i = 1,NGrid

c          write(6,*)"     "
c          write(6,*)"----------------------"
c          write(6,*)"for grid: ", i

        IF(NMO.eq.1) then
           EX_HF_DENSB(i)=EX_HF_DENSA(i)
c          RB(i) = RA(i)
c          TB(i) = TA(i)
         endif
          D1 = RA(i)
          D2 = RB(i)
          ds1 = abs(D1-D2)
c        IF (IMethod.eq.2) THEN
c        IF(NMO.eq.1) FIT_DENSB(i)=FIT_DENSA(i)
c          DF1 = FIT_DENSA(i)
c          DF2 = FIT_DENSB(i)
c        endif
c
cccccccccccccccccccc
c       IF(i.eq.5.or.i.eq.6) then
c         ndrv = 1
c  O
c        D1 = 4.5384906103506800d-4  + 0.0000005d0
c        D1 = 4.5384906103506800d-4  - 0.0000005d0
c        D1 = 1.1649710650678880d-3  + 0.000001d0
c        D1 = 1.1649710650678880d-3  - 0.000001d0
c
c  H
c
c        D1 = 0.170533691444758  + 0.0001d0
c        D1 = 0.170533691444758  - 0.0001d0
c        D1 = 1.0665313216417950d-5  + 0.00000005d0
c        D1 = 1.0665313216417950d-5  - 0.00000005d0
c      endif
ccccccccccccccccccc
         fc1 = zero
         fc2 = zero
         fcor = zero
         UCA = zero
         UCB = zero
         DSA =dsmall
         DSB =dsmall
         DDA = zero
         DDB = zero
         ASA = zero
         ASSA = zero
         ASB = zero
         ASSB = zero
         dASAdra = zero
         dASAdga = zero
         dASAdta = zero
         dASAdla = zero
         dASAdVxa = zero
         dASBdra = zero
         dASBdga = zero
         dASBdta = zero
         dASBdla = zero
         dASBdVxa =zero
         dASBdrb = zero
         dASBdgb = zero
         dASBdtb = zero
         dASBdlb = zero
         dASBdVxb =zero
         dASAdrb = zero
         dASAdgb = zero
         dASAdtb = zero
         dASAdlb = zero
         dASAdVxb =zero
         dASSAdra = zero
         dASSAdga = zero
         dASSAdta = zero
         dASSAdla = zero
         dASSAdVxa =zero
         dASSAdrb = zero
         dASSAdgb = zero
         dASSAdtb = zero
         dASSAdlb = zero
         dASSAdVxb =zero
         dASSBdrb = zero
         dASSBdgb = zero
         dASSBdtb = zero
         dASSBdlb = zero
         dASSBdVxb= zero
         dASSBdra = zero
         dASSBdga = zero
         dASSBdta = zero
         dASSBdla = zero
         dASSBdVxa =zero
         EFNA = 1.0d0
         EFN2A = 1.0d0
         EFNB = 1.0d0
         EFN2B = 1.0d0
         UX1= -zero
         UX2= -zero
         gra=zero
         grb=zero
         tauA =zero
         tauB =zero
         DLapA =zero
         DLapB =zero
         dfdra = zero
         dfdga = zero
         dfdta = zero
         dfdla = zero
         dfdrb = zero
         dfdgb = zero
         dfdtb = zero
         dfdlb = zero
         dfdVxa=zero
         dfdVxb=zero
         dfadra = zero
         dfadga = zero
         dfadta = zero
         dfadla = zero
         dfadrb = zero
         dfadgb = zero
         dfadtb = zero
         dfadlb = zero
         dfadVxa = zero
         dfadVxb = zero
         dfbdra = zero
         dfbdga = zero
         dfbdta = zero
         dfbdla = zero
         dfbdrb = zero
         dfbdgb = zero
         dfbdtb = zero
         dfbdlb = zero
         dfbdVxa = zero
         dfbdVxb = zero
         dfcdraa = zero
         dfcdgaa = zero
         dfcdtaa = zero
         dfcdlaa = zero
         dfcdVxaa =zero
         dfcdrbb = zero
         dfcdgbb = zero
         dfcdtbb = zero
         dfcdlbb = zero
         dfcdVxbb= zero
         dM1dra = zero
         dM1dta = zero
         dM1dla = zero
         dM1dga = zero
         dM1dVxa = zero
         dM2dra = zero
         dM2dta = zero
         dM2dla = zero
         dM2dga = zero
         dM2dVxa = zero
         dM1drb = zero
         dM1dtb = zero
         dM1dlb = zero
         dM1dgb = zero
         dM1dVxb = zero
         dM2drb = zero
         dM2dtb = zero
         dM2dlb = zero
         dM2dgb = zero
         dM2dVxb = zero
         dg3dra = zero
         dg3dta = zero
         dg3dga = zero
         dg3drb = zero
         dg3dtb = zero
         dg3dgb = zero
         dfc1g1 = zero
         dfc2g2 = zero
         SM1A = zero
         SM2A = zero
         SM1B = zero
         SM2B = zero
         ttgg = zero
         ttg1 = zero
         ttg2 = zero
         ttdB = zero
         ttdA = zero
         dfpardVxa = zero
         dfpardra = zero
         dfpardga = zero
         dfpardla = zero
         dfpardta = zero
         dfpardVxa = zero
         dfpardVxb = zero
         dfpardrb = zero
         dfpardgb = zero
         dfpardlb = zero
         dfpardtb = zero
         dfopdra =zero
         dfopdga =zero
         dfopdla =zero
         dfopdta =zero
         dfopdVxa =zero
         dfopdrb =zero
         dfopdgb =zero
         dfopdlb =zero
         dfopdtb =zero
         dfopdVxb =zero
         gc1 = zero
         gcc1 = zero
         ggg = zero
         ggp = zero
         gc2 = zero
         gcc2 = zero
         dNdra = zero
         dNdga = zero
         dNdta = zero
         dNdla = zero
         dNdVxa = zero
         dNdrb = zero
         dNdgb = zero
         dNdtb = zero
         dNdVxb = zero
         qqx = 0.0d0
         eey1=0.0d0
         eey2=0.0d0
         dNdVxa = zero
         dNdVxb = zero
         qqx = 0.0d0
         smtx = zero
         tx = 0.0d0
         qqfa = 0.0d0
         qqfb = 0.0d0
         xxx = zero
         yy1 = zero
         yy2 = zero
         fp1=zero
         fp2=zero
         VCPAR = zero
         Fop = zero
         Fpar = zero
         ndrvA = 2
         ndrvB = 2
c
             UX1 = -dsmall2
             UX2 = -dsmall2
             VX1=EX_HF_DENSA(i)
             VX2=EX_HF_DENSB(i)

        IF((VX1.ge.dsmall).or.(VX2.ge.dsmall)) GO TO 666
       IF((D1.gt.detol).OR.(D2.gt.detol)) then

         IF(D1.gt.detol) then
              UX1 = VX1/D1
cccccccccccccccccccccccccccc
c         if(EX_HF_DENSA(i).ge.-dsmall) then
c            write(*,*)'positive VxA=  ',EX_HF_DENSA(i)
c       CALL QCrash('positive VxA')
c            write(*,*) 'rhoA_fit =  ', DF1, 'i=  ',i
c          endif  ! EX
c
cccccccccccccccccccc
c         if(ndrv.eq.1) then
c             VX1 = -0.1643252735541850d0 + 0.0001d0
c             VX1 = -0.1643252735541850d0 - 0.0001d0
c             UX1= VX1/D1
c         endif
cccccccccccccccc
               gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
               tauA = TA(i)
               DLapA = RLapA(i)
c               write(6,*)"emil's RhoA: ",D1
c               write(6,*)"emil's GRhoA: ",gra
c               write(6,*)"emil's TauA: ",tauA 
c               write(6,*)"emil's LapA: ",DLAPA 
c               write(6,*)"emil's UA: ",VX1 
cccccccccccccc
c         if(ndrv.eq.1) tauA = 37.15155071715040d0 + 0.001d0
c         if(ndrv.eq.1) tauA = 37.15155071715040d0 - 0.001d0
c         if(ndrv.eq.1) tauA = 36.0147379088890d0 - 0.001d0
c
c         if(ndrv.eq.1) DLapA = 0.6556419974150890d0 + 0.0005d0
c         if(ndrv.eq.1) DLapA = 0.6556419974150890d0 - 0.00050d0
cccccccccccccc
c brsc is the routine doing the main RSC numerics
           CALL brsc(D1,DF1,gra,DLapA,tauA,VX1,UX1,EFNA,EFN2A,dNdra,
     $        dNdga,dNdta,dNdla,dNdVxa,dM1dra,dM2dra,dM1dta,dM2dta,
     $        dM1dla,dM2dla,dM1dga,dM2dga,dM1dVxa,dM2dVxa,dg3dra,
     $        dg3dta,dg3dga,SM1A,SM2A,DSA,ndrvA)
       else  ! D1
               SM1A = zero
               SM2A = zero
               DSA = dsmall
               gra = 0.0d0
               tauA = zero
               DLapA = zero
               UX1 = -dsmall
               EFNA = 1.0d0
               EFN2A = 1.0d0
               dNdra= zero
               dNdga= zero
               dNdta= zero
               dNdla= zero
               dNdVxa= zero
               dM1dra =zero
               dM2dra = zero
               dM1dta = zero
               dM2dta = zero
               dM1dla = zero
               dM2dla = zero
               dM1dga = zero
               dM2dga = zero
               dM1dVxa=zero
               dM2dVxa=zero
               dg3dra=zero
               dg3dta=zero
               dg3dga=zero
       endif ! D1, VX1
c
        IF(D2.gt.detol) then
              UX2 = VX2/D2
cccccccccccccccccccccccccccccccc
c          if((EX_HF_DENSB(i).ge.dsmall2).and.(NMO.eq.2)) then
c      write(*,*)'positive VxB=  ',EX_HF_DENSB(i)
c       CALL QCrash('positive VxB')
c           write(*,*) 'rhoB =  ', DF2, 'i=  ',i
c          endif  ! EX

c         if(ndrv.eq.1) then
c             VX2 = -72.75196487421990d0 + 0.005d0
c             VX2 = -72.75196487421990d0 - 0.005d0
c             UX2= VX2/D2
c         endif
               grb = D1RB(i,1)**2 + D1RB(i,2)**2 + D1RB(i,3)**2
               tauB = TB(i)
               DLapB = RLapB(i)

c               write(6,*)"emil's RhoB: ",D2
c               write(6,*)"emil's GRhoB: ",grb
c               write(6,*)"emil's TauB: ",tauB 
c               write(6,*)"emil's LapB: ",DLAPB 
c               write(6,*)"emil's UB: ",VX2 

           CALL brsc(D2,DF2,grb,DLapB,tauB,VX2,UX2,EFNB,EFN2B,dNdrb,
     $        dNdgb,dNdtb,dNdlb,dNdVxb,dM1drb,dM2drb,dM1dtb,dM2dtb,
     $        dM1dlb,dM2dlb,dM1dgb,dM2dgb,dM1dVxb,dM2dVxb,dg3drb,
     $        dg3dtb,dg3dgb,SM1B,SM2B,DSB,ndrvB)

           else  !  D2
               SM1B= zero
               SM2B= zero
               DSB= dsmall
               grb = zero
               tauB = zero
               DLapB = zero
               UX2 = -dsmall
               EFNB = 1.0d0
               EFN2B = 1.0d0
               dNdrb= zero
               dNdgb= zero
               dNdtb= zero
               dNdlb= zero
               dNdUxb= zero
               dM1drb =zero
               dM2drb = zero
               dM1dtb = zero
               dM2dtb = zero
               dM1dlb = zero
               dM2dlb = zero
               dM1dgb = zero
               dM2dgb = zero
               dM1dVxb=zero
               dM2dVxb=zero
               dg3drb=zero
               dg3dtb=zero
               dg3dgb=zero
       endif  ! VX2
c
                   fcor = 0.0d0
                   ffcor =0.0d0
                   Fop = 0.0d0
            IF((ndrvA.eq.0).or.(ndrvB.eq.0)) GO TO 777
c      If((NA.gt.0).and.(NB.gt.0)) then
         if((D1.gt.detol.and.ndrvA.eq.2).and.(D2.gt.detol
     $      .and.ndrvB.eq.2)) then
c
         IF(EFN2B.gt.dtol2) gcc1 = (1.0d0-EFN2A)/EFN2B
         IF(EFN2A.gt.dtol2) gcc2 = (1.0d0-EFN2B)/EFN2A
c           
         IF(EFNB.gt.dtol2) then
           gc1 = (1.0d0-EFNA)/EFNB
           ggg= gc1-1.0d0-delt
         IF((gc1.gt.smth1).and.(gc1.lt.smth2)) then
           fc1 = 1.0d0-ggg*ggg/(4.0D0*delt)
           dfc1g1 = -ggg/(2.0D0*delt)
        elseif (gc1.ge.smth2) then
           fc1 = 1.0d0
           dfc1g1 = 0.0d0
        else
            fc1 = gc1
            dfc1g1 = 1.0D0
         endif   ! gc1
        endif
c
c        write(6,*)"emil alpha neff: ", EFNA
c        write(6,*)"emil beta  neff: ", EFNB

       IF(EFNA.gt.dtol2) then
           gc2 = (1.0d0-EFNB)/EFNA
c          gc2 = 1.0d0/EFNA-EFNB/EFNA
           ggp= gc2-1.0d0-delt
        IF((gc2.gt.smth1).and.(gc2.lt.smth2)) then
           fc2 = 1.0d0-ggp*ggp/(4.0D0*delt)
           dfc2g2 = -ggp/(2.0D0*delt)
        elseif (gc2.ge.smth2) then
           fc2 = 1.0d0
           dfc2g2 = 0.0d0
        else
            fc2 = gc2
            dfc2g2 = 1.0D0
         endif   ! gc2
        endif
c
           tx = 1.0d0
           qqx = 0.50d0
           smtx = 0.0d0
           xxx=zero
c        fcor = fc1*0.50d0 + fc2*0.50d0
        if((abs(fc1).gt.detol).or.(abs(fc2).gt.detol)) then
             fp1 = fc1*fc1
             fp2 = fc2*fc2
           xxx = (fc1-fc2)/(fp1+fp2)
            smtx = smoth*xxx
c           tx = exp(smtx)
          if(abs(smtx).lt.99.0d0) then
            tx = exp(smtx)
            qqx = 1.0d0/(tx+1.0d0)
          elseif(smtx.lt.-99.0d0) then
             tx= zero
             qqx= 1.0d0
          else
             tx= zero
             qqx= zero
          endif  !smtx
         fcor = fc1*qqx + fc2*(1.0d0-qqx)
           if(gc1.le.-dsmall) gc1=0.0d0
           if(gcc1.le.-dsmall) gcc1=0.0d0
           if(gc2.le.-dsmall) gc2=0.0d0
           if(gcc2.le.-dsmall) gcc2=0.0d0
c        ffcor = min(gc1,gc2,1.0d0)
         ffcor = min(gcc1,gcc2,1.0d0)
        else  !  fc1,fc2
           xxx=zero
           qqx = 0.0d0
           tx = 0.0d0
           smtx=zero
          fcor = zero
         endif  !  fc1,fc2
c         write(6,*)"emil F1: ", fc1
c         write(6,*)"emil F2: ", fc2
c         write(6,*)"emil HPZ: ", qqx
c
          ttdA = zero
          ttdB = zero
        IF(D2.gt.dtol) then
          ttdB = UX2*D1/D2
        endif
        IF(D1.gt.dtol) then
          ttdA = UX1*D2/D1
        endif
          ttg2=D1*UX2
          ttg1=D2*UX1
           ttgg = ttg1+ttg2
         Fop=ACOP*0.50d0*fcor*ttgg
c         write(6,*)"emil F value: ", fcor
       else  ! NA
           gc1=zero
           gcc1 = zero
           fc1 = zero
           dfc1g1 = zero
           gc2=zero
           gcc2 = zero
           fc2 = zero
           dfc2g2 = zero
           xxx = zero
           tx = 0.0d0
           qqx = 0.0d0
           smtx = zero
           fcor = zero
           tx = 0.0d0
          Fop = 0.0d0
          ttdA=0.0d0
          ttdB=0.0d0
      endif   ! NA
cccccccccccc
 777  CONTINUE
         IF((ndrvA.eq.0).and.(ndrvB.eq.0)) GO TO 666
            yy1 = zero
            smtx2 = zero
            eey1 = 0.0d0
            ty1 = 0.0d0
            yy2 = zero
            smtx3 = zero
            eey2=0.0d0
            ty2 = 0.0d0
c
              ASSA =0.0d0
              A2SA = 0.0d0
              AASA = 0.0d0
      IF((D1.gt.detol).and.(ndrvA.eq.2)) then
          IF(abs(SM2A).gt.dtol) then
           ASA=(1.0d0-EFNA-fcor*EFNB)/SM2A
           A2SA=(1.0d0-EFN2A-ffcor*EFN2B)/SM2A
               endif
            DDA = DSA/(3.0d0*D1)
            yy1 = ASA-DDA
            smtx2 = smoth2*yy1
c            write(6,*)"emil's alpha A1:", ASA
c            write(6,*)"emil's alpha A2:", DDA
c           ty1 = exp(smtx2)
          if(abs(smtx2).lt.99.0d0) then
            ty1 = exp(smtx2)
            eey1 = 1.0d0/(ty1+1.0d0)
          elseif(smtx2.lt.-99.0d0) then
              ty1= zero
              eey1 = 1.0d0
          else
             ty1= zero
             eey1 = zero
          endif
            ASSA = ASA*eey1+DDA*(1.0d0-eey1)
            AASA = min(A2SA,DDA)
        endif  ! D1
ccccccccccccccccccccc
              ASSB=0.0d0
              A2SB = 0.0d0
              BBSB = 0.0d0
      IF((D2.gt.detol).and.(ndrvB.eq.2)) then
         IF(abs(SM2B).gt.dtol) then
            ASB=(1.0d0-EFNB-fcor*EFNA)/SM2B
            A2SB=(1.0d0-EFN2B-ffcor*EFN2A)/SM2B
               endif
            DDB = DSB/(3.0d0*D2)
            yy2 = ASB-DDB
            smtx3 = smoth2*yy2
c           ty2 = exp(smtx3)
          if(abs(smtx3).lt.99.0d0) then
            ty2 = exp(smtx3)
            eey2 = 1.0d0/(ty2+1.0d0)
          elseif(smtx3.lt.-99.0d0) then
              ty2= zero
              eey2 = 1.0d0
          else
              ty2= zero
              eey2 = zero
          endif
            ASSB = ASB*eey2+DDB*(1.0d0-eey2)
            BBSB = min(A2SB,DDB)
        endif   ! D2
c            write(6,*)"emil's beta A1:", ASB
c            write(6,*)"emil's beta A2:", DDB
ccccccccc
           UCA=0.0d0
           UCB=0.0d0
          if(ndrvA.eq.2)  UCA = - ASSA*SM1A
          if(ndrvB.eq.2)  UCB = - ASSB*SM1B
           VCPAR = (D1*UCA+D2*UCB)*0.50d0
c        else  !  NMO
c          VCPAR= zero
c        endif
ccccccccccccccc
         Fpar = ACPAR*VCPAR
         F(i)= Fop + Fpar
c          IF(NMO.eq.2) then
           apar1=ACPAR
c          else
c           apar1=0.0d0
c          endif
         OD1 = ACOP*4.0d0*D1*ffcor*EFN2B + apar1*2.0d0*D1*AASA*SM2A
         OD2 = ACOP*4.0d0*D2*ffcor*EFN2A + apar1*2.0d0*D2*BBSB*SM2B
          IF(OD1.le.(-dsmall)) OD1=0.0d0
          IF(OD2.le.(-dsmall)) OD2=0.0d0
            ODEL(i) = OD1+OD2
c           ODEL(i) = OD1+OD2 + ds1
            ODELW(i) = ODEL(i)*W(i)
        IF(ICEN.ne.NCENTER ) then
            NCENTER = ICEN
         endif
c          ODDSUM(NCENTER) = ODDSUM(NCENTER)+ODELW(i)
          
c       write(*,*) 'emil Fop = ', Fop
c       write(*,*) 'Fpar = ', Fpar
c
c Begin building the SCF potential:
        IF (IMethod.eq.2) THEN
c
            IF((ndrvA.eq.0).or.(ndrvB.eq.0)) GO TO 888
         If((D1.gt.detol.and.ndrvA.eq.2).and.(D2.gt.detol
     $      .and.ndrvB.eq.2)) then
       IF(EFNB.gt.dtol2) then
          dfadra = -dfc1g1*dNdra/EFNB
          dfadga = -dfc1g1*dNdga/EFNB
          dfadta = -dfc1g1*dNdta/EFNB
          dfadla = -dfc1g1*dNdla/EFNB
          dfadVxa = -dfc1g1*dNdVxa/EFNB
          dfadrb = -dfc1g1*dNdrb*(1.0d0-EFNA)/EFNB**2
          dfadgb = -dfc1g1*dNdgb*(1.0d0-EFNA)/EFNB**2
          dfadtb = -dfc1g1*dNdtb*(1.0d0-EFNA)/EFNB**2
          dfadlb = -dfc1g1*dNdlb*(1.0d0-EFNA)/EFNB**2
          dfadVxb = -dfc1g1*dNdVxb*(1.0d0-EFNA)/EFNB**2
       else
          dfadra = zero
          dfadga = zero
          dfadta = zero
          dfadla = zero
         dfadVxa = zero
          dfadrb = zero
          dfadgb = zero
          dfadtb = zero
          dfadlb = zero
         dfadVxb = zero
       endif
c        write(6,*)"emil's F1 alpha rho: ", dfadra
c        write(6,*)"emil's F1 alpha GAA: ", dfadga
c        write(6,*)"emil's F1 alpha tau: ", dfadta
c        write(6,*)"emil's F1 alpha lap: ", dfadla
c        write(6,*)"emil's F1 alpha ex : ", dfadVxa
c        write(6,*)"emil's F1 beta  rho: ", dfadrb
c        write(6,*)"emil's F1 beta  GAA: ", dfadgb
c        write(6,*)"emil's F1 beta  tau: ", dfadtb
c        write(6,*)"emil's F1 beta  lap: ", dfadlb
c        write(6,*)"emil's F1 beta  ex : ", dfadVxb
        IF(EFNA.gt.dtol2) then
          dfbdrb = -dfc2g2*dNdrb/EFNA
          dfbdgb = -dfc2g2*dNdgb/EFNA
          dfbdtb = -dfc2g2*dNdtb/EFNA
          dfbdlb = -dfc2g2*dNdlb/EFNA
          dfbdVxb = -dfc2g2*dNdVxb/EFNA
          dfbdra = -dfc2g2*dNdra*(1.0d0-EFNB)/EFNA**2
          dfbdga = -dfc2g2*dNdga*(1.0d0-EFNB)/EFNA**2
          dfbdta = -dfc2g2*dNdta*(1.0d0-EFNB)/EFNA**2
          dfbdla = -dfc2g2*dNdla*(1.0d0-EFNB)/EFNA**2
          dfbdVxa = -dfc2g2*dNdVxa*(1.0d0-EFNB)/EFNA**2
        else
          dfbdrb = zero
          dfbdgb = zero
          dfbdtb = zero
          dfbdlb = zero
         dfbdVxb = zero
          dfbdra = zero
          dfbdga = zero
          dfbdta = zero
          dfbdla = zero
         dfbdVxa = zero
        endif
c        write(6,*)"emil's F2 alpha rho: ", dfbdra
c        write(6,*)"emil's F2 alpha GAA: ", dfbdga
c        write(6,*)"emil's F2 alpha tau: ", dfbdta
c        write(6,*)"emil's F2 alpha lap: ", dfbdla
c        write(6,*)"emil's F2 alpha ex : ", dfbdVxa
c        write(6,*)"emil's F2 beta  rho: ", dfbdrb
c        write(6,*)"emil's F2 beta  GAA: ", dfbdgb
c        write(6,*)"emil's F2 beta  tau: ", dfbdtb
c        write(6,*)"emil's F2 beta  lap: ", dfbdlb
c        write(6,*)"emil's F2 beta  ex : ", dfbdVxb
c
        if ((abs(fc1).gt.dtol).or.(abs(fc2).gt.dtol)) then
                fp1= fc1*fc1
                fp2= fc2*fc2
          qqfa = smtx*qqx*qqx*tx/(fp1+fp2)
           fffa =(fc2*fc2-fc1*fc1 -2.0d0*fc1*fc2)
           fffb =(fc2*fc2-fc1*fc1 +2.0d0*fc1*fc2)
c
          dfcdraa = dfbdra+(dfadra-dfbdra)*qqx
     $    -fffa*dfbdra*qqfa-fffb*dfadra*qqfa

          dfcdgaa = dfbdga+(dfadga-dfbdga)*qqx
     $    -fffa*dfbdga*qqfa-fffb*dfadga*qqfa

          dfcdtaa = dfbdta+(dfadta-dfbdta)*qqx
     $    -fffa*dfbdta*qqfa-fffb*dfadta*qqfa

          dfcdlaa = dfbdla+(dfadla-dfbdla)*qqx
     $    -fffa*dfbdla*qqfa-fffb*dfadla*qqfa

          dfcdVxaa = dfbdVxa+(dfadVxa-dfbdVxa)*qqx
     $    -fffa*dfbdVxa*qqfa-fffb*dfadVxa*qqfa

          dfcdrbb = dfbdrb+(dfadrb-dfbdrb)*qqx
     $    -fffa*dfbdrb*qqfa-fffb*dfadrb*qqfa

          dfcdgbb = dfbdgb+(dfadgb-dfbdgb)*qqx
     $    -fffa*dfbdgb*qqfa-fffb*dfadgb*qqfa

          dfcdtbb = dfbdtb+(dfadtb-dfbdtb)*qqx
     $    -fffa*dfbdtb*qqfa-fffb*dfadtb*qqfa

          dfcdlbb = dfbdlb+(dfadlb-dfbdlb)*qqx
     $    -fffa*dfbdlb*qqfa-fffb*dfadlb*qqfa

          dfcdVxbb = dfbdVxb+(dfadVxb-dfbdVxb)*qqx
     $    -fffa*dfbdVxb*qqfa-fffb*dfadVxb*qqfa

        else
          dfcdraa = zero
          dfcdgaa = zero
          dfcdtaa = zero
          dfcdlaa = zero
         dfcdVxaa = zero
          dfcdrbb = zero
          dfcdgbb = zero
          dfcdtbb = zero
          dfcdlbb = zero
          dfcdVxbb= zero
          qqfa = 0.50d0
          qqx = 0.50d0
        endif
C        write(6,*)"emil's F alpha rho: ", dfcdraa
C        write(6,*)"emil's F alpha GAA: ", dfcdgaa
C        write(6,*)"emil's F alpha tau: ", dfcdtaa
C        write(6,*)"emil's F alpha lap: ", dfcdlaa
C        write(6,*)"emil's F alpha ex : ", dfcdVxaa
C        write(6,*)"emil's F beta  rho: ", dfcdrbb
C        write(6,*)"emil's F beta  GAA: ", dfcdgbb
C        write(6,*)"emil's F beta  tau: ", dfcdtbb
C        write(6,*)"emil's F beta  lap: ", dfcdlbb
C        write(6,*)"emil's F beta  ex : ", dfcdVxbb
          dfopdra = ACOP*0.50d0*(dfcdraa*ttgg+fcor*(UX2-ttdA))
          dfopdga = ACOP*0.50d0*dfcdgaa*ttgg
          dfopdta = ACOP*0.50d0*dfcdtaa*ttgg
          dfopdla = ACOP*0.50d0*dfcdlaa*ttgg
          dfopdVxa = ACOP*0.50d0*dfcdVxaa*ttgg
          dfopdVxa = dfopdVxa+ ACOP*0.50d0*fcor*D2/D1
          dfopdrb = ACOP*0.50d0*(dfcdrbb*ttgg+fcor*(UX1-ttdB))
          dfopdgb = ACOP*0.50d0*dfcdgbb*ttgg
          dfopdtb = ACOP*0.50d0*dfcdtbb*ttgg
          dfopdlb = ACOP*0.50d0*dfcdlbb*ttgg
          dfopdVxb = ACOP*0.50d0*dfcdVxbb*ttgg
          dfopdVxb = ACOP*0.50d0*(dfcdVxbb*ttgg+fcor*D1/D2)
      endif  ! NA vs NB
c=====up to here the opposite spin SCF part is done
c
 888      CONTINUE
            facta = zero
            factb=zero
c       IF((NMO.eq.2).and.(NA.ne.NB)) then
c
         if(abs(SM2A).gt.dtol2) facta= ASA/SM2A
       If((D1.gt.detol).and.(abs(SM2A).gt.dtol2).and.(ndrvA.eq.2)) then
            dASAdra = -(dNdra+EFNB*dfcdraa)/SM2A - dM2dra*facta
            dASAdga = -(dNdga+EFNB*dfcdgaa)/SM2A - dM2dga*facta
            dASAdta = -(dNdta+EFNB*dfcdtaa)/SM2A - dM2dta*facta
            dASAdla = -(dNdla+EFNB*dfcdlaa)/SM2A - dM2dla*facta
            dASAdVxa = -(dNdVxa+EFNB*dfcdVxaa)/SM2A - dM2dVxa*facta
         else
            dASAdra = zero
            dASAdga = zero
            dASAdta = zero
            dASAdla = zero
           dASAdVxa = zero
         endif
c
       if((D1.gt.detol).and.(abs(SM2B).gt.dtol).and.(D2.gt.detol).
     $ and.(ndrvB.eq.2)) then
c      if(D1.gt.dtol2.and.(abs(SM2B).gt.dtol2)) then
            dASBdra = - (fcor*dNdra+ EFNA*dfcdraa)/SM2B
            dASBdga = - (fcor*dNdga+ EFNA*dfcdgaa)/SM2B
            dASBdta = - (fcor*dNdta+ EFNA*dfcdtaa)/SM2B
            dASBdla = - (fcor*dNdla+ EFNA*dfcdlaa)/SM2B
          dASBdVxa = -(fcor*dNdVxa+ EFNA*dfcdVxaa)/SM2B
       else
            dASBdra = zero
            dASBdga = zero
            dASBdta = zero
            dASBdla = zero
           dASBdVxa = zero
       endif
c
          if(abs(SM2B).gt.dtol2) factb= ASB/SM2B
       IF((D2.gt.detol).and.(abs(SM2B).gt.dtol2).and.(ndrvB.eq.2)) then
            dASBdrb = -(dNdrb+EFNA*dfcdrbb)/SM2B - dM2drb*factb
            dASBdgb = -(dNdgb+EFNA*dfcdgbb)/SM2B - dM2dgb*factb
            dASBdtb = -(dNdtb+EFNA*dfcdtbb)/SM2B - dM2dtb*factb
            dASBdlb = -(dNdlb+EFNA*dfcdlbb)/SM2B - dM2dlb*factb
          dASBdVxb = -(dNdVxB+EFNA*dfcdVxbb)/SM2B-dM2dVxb*factb
        else
            dASBdrb = zero
            dASBdgb = zero
            dASBdtb = zero
            dASBdlb = zero
           dASBdVxb = zero
        endif
      IF((D2.gt.detol).and.(abs(SM2A).gt.dtol).and.(D1.gt.detol)
     $  .and.(ndrvA.eq.2)) then
            dASAdrb = - (fcor*dNdrb+ EFNB*dfcdrbb)/SM2A
            dASAdgb = - (fcor*dNdgb+ EFNB*dfcdgbb)/SM2A
            dASAdtb = - (fcor*dNdtb+ EFNB*dfcdtbb)/SM2A
            dASAdlb = - (fcor*dNdlb+ EFNB*dfcdlbb)/SM2A
         dASAdVxb = - (fcor*dNdVxb+ EFNB*dfcdVxbb)/SM2A
       else
            dASAdrb = zero
            dASAdgb = zero
            dASAdtb = zero
            dASAdlb = zero
           dASAdVxb = zero
       endif
c       write(6,*)"emil's alpha A1 alpha rho: ", dASAdra
c       write(6,*)"emil's alpha A1 alpha GAA: ", dASAdga
c       write(6,*)"emil's alpha A1 alpha tau: ", dASAdta
c       write(6,*)"emil's alpha A1 alpha lap: ", dASAdla
c       write(6,*)"emil's alpha A1 alpha ex : ", dASAdVxa
c       write(6,*)"emil's alpha A1 beta  rho: ", dASAdrb
c       write(6,*)"emil's alpha A1 beta  GAA: ", dASAdgb
c       write(6,*)"emil's alpha A1 beta  tau: ", dASAdtb
c       write(6,*)"emil's alpha A1 beta  lap: ", dASAdlb
c       write(6,*)"emil's alpha A1 beta  ex : ", dASAdVxb
c
c       ! beta part
c       write(6,*)"emil's beta A1 alpha rho: ", dASBdra
c       write(6,*)"emil's beta A1 alpha GAA: ", dASBdga
c       write(6,*)"emil's beta A1 alpha tau: ", dASBdta
c       write(6,*)"emil's beta A1 alpha lap: ", dASBdla
c       write(6,*)"emil's beta A1 alpha ex : ", dASBdVxa
c       write(6,*)"emil's beta A1 beta  rho: ", dASBdrb
c       write(6,*)"emil's beta A1 beta  GAA: ", dASBdgb
c       write(6,*)"emil's beta A1 beta  tau: ", dASBdtb
c       write(6,*)"emil's beta A1 beta  lap: ", dASBdlb
c       write(6,*)"emil's beta A1 beta  ex : ", dASBdVxb
c
          ppaa = eey1 -smoth2*eey1*eey1*ty1*yy1
          ppba = eey2 -smoth2*eey2*eey2*ty2*yy2
c
         IF(ndrvA.eq.2) then
           dASSAdra = dg3dra+(dASAdra-dg3dra)*ppaa
           dASSAdga = dg3dga+(dASAdga-dg3dga)*ppaa
           dASSAdta = dg3dta+(dASAdta-dg3dta)*ppaa
           dASSAdla = dASAdla*ppaa
           dASSAdVxa = dASAdVxa*ppaa
           dASSAdrb = dASAdrb*ppaa
           dASSAdgb = dASAdgb*ppaa
           dASSAdtb = dASAdtb*ppaa
           dASSAdlb = dASAdlb*ppaa
          dASSAdVxb = dASAdVxb*ppaa
         else
           dASSAdra =  0.0d0
           dASSAdga =  0.0d0
           dASSAdta =  0.0d0
           dASSAdla =  0.0d0
          dASSAdVxa =  0.0d0
           dASSAdrb =  0.0d0
           dASSAdgb =  0.0d0
           dASSAdtb =  0.0d0
           dASSAdlb =  0.0d0
          dASSAdVxb =  0.0d0
         endif
c        write(*,*) 'emil alpha A', AASA
c        write(6,*)"emil's alpha A alpha rho: ", dASSAdra
c        write(6,*)"emil's alpha A alpha GAA: ", dASSAdga
c        write(6,*)"emil's alpha A alpha tau: ", dASSAdta
c        write(6,*)"emil's alpha A alpha lap: ", dASSAdla
c        write(6,*)"emil's alpha A alpha ex : ", dASSAdVxa
c        write(6,*)"emil's alpha A beta  rho: ", dASSAdrb
c        write(6,*)"emil's alpha A beta  GAA: ", dASSAdgb
c        write(6,*)"emil's alpha A beta  tau: ", dASSAdtb
c        write(6,*)"emil's alpha A beta  lap: ", dASSAdlb
c        write(6,*)"emil's alpha A beta  ex : ", dASSAdVxb
c
         IF(ndrvB.eq.2) then
           dASSBdrb = dg3drb+ (dASBdrb- dg3drb)*ppba
           dASSBdgb = dg3dgb+ (dASBdgb- dg3dgb)*ppba
           dASSBdtb = dg3dtb+ (dASBdtb- dg3dtb)*ppba
           dASSBdlb = dASBdlb*ppba
           dASSBdVxb = dASBdVxb*ppba
           dASSBdra = dASBdra*ppba
           dASSBdga = dASBdga*ppba
           dASSBdta = dASBdta*ppba
           dASSBdla = dASBdla*ppba
           dASSBdVxa = dASBdVxa*ppba
         else
           dASSBdrb =  0.0d0
           dASSBdgb =  0.0d0
           dASSBdtb =  0.0d0
           dASSBdlb =  0.0d0
          dASSBdVxb =  0.0d0
           dASSBdra =  0.0d0
           dASSBdga =  0.0d0
           dASSBdta =  0.0d0
           dASSBdla =  0.0d0
          dASSBdVxa =  0.0d0
          endif
c        write(*,*) 'emil beta A', AASB
c        write(6,*)"emil's beta A alpha rho: ", dASSBdra
c        write(6,*)"emil's beta A alpha GAA: ", dASSBdga
c        write(6,*)"emil's beta A alpha tau: ", dASSBdta
c        write(6,*)"emil's beta A alpha lap: ", dASSBdla
c        write(6,*)"emil's beta A alpha ex : ", dASSBdVxa
c        write(6,*)"emil's beta A beta  rho: ", dASSBdrb
c        write(6,*)"emil's beta A beta  GAA: ", dASSBdgb
c        write(6,*)"emil's beta A beta  tau: ", dASSBdtb
c        write(6,*)"emil's beta A beta  lap: ", dASSBdlb
c        write(6,*)"emil's beta A beta  ex : ", dASSBdVxb
c
      IF ((D1.gt.detol).and.(ndrvA.eq.2)) then
            dfpardra = -ACPAR*0.50d0*(ASSA*SM1A+D1*dASSAdra*SM1A
     $               + D2*dASSBdra*SM1B + D1*ASSA*dM1dra)
            dfpardga = -ACPAR*0.50d0*(D1*dASSAdga*SM1A
     $               + D2*dASSBdga*SM1B + D1*ASSA*dM1dga)
            dfpardta = -ACPAR*0.50d0*(D1*dASSAdta*SM1A
     $               + D2*dASSBdta*SM1B + D1*ASSA*dM1dta)
            dfpardla = -ACPAR*0.50d0*(D1*dASSAdla*SM1A
     $               + D2*dASSBdla*SM1B + D1*ASSA*dM1dla)
            dfpardVxa = -ACPAR*0.50d0*(D1*dASSAdVxa*SM1A
     $                + D2*dASSBdVxa*SM1B + D1*ASSA*dM1dVxa)
         else
            dfpardra = zero
            dfpardga = zero
            dfpardta = zero
            dfpardla = zero
           dfpardVxa = zero
         endif
      IF ((D2.gt.dtol).and.(ndrvB.eq.2)) then
          dfpardrb = -ACPAR*0.50d0*(ASSB*SM1B+D2*dASSBdrb*SM1B
     $             + D1*dASSAdrb*SM1A + D2*ASSB*dM1drb)
          dfpardgb = -ACPAR*0.50d0*(D2*dASSBdgb*SM1B
     $             + D1*dASSAdgb*SM1A + D2*ASSB*dM1dgb)
          dfpardtb = -ACPAR*0.50d0*(D2*dASSBdtb*SM1B
     $             + D1*dASSAdtb*SM1A + D2*ASSB*dM1dtb)
          dfpardlb = -ACPAR*0.50d0*(D2*dASSBdlb*SM1B
     $             + D1*dASSAdlb*SM1A + D2*ASSB*dM1dlb)
          dfpardVxb = -ACPAR*0.50d0*(D2*dASSBdVxb*SM1B
     $               + D1*dASSAdVxb*SM1A + D2*ASSB*dM1dVxb)
       else
          dfpardrb = zero
          dfpardgb = zero
          dfpardtb = zero
          dfpardlb = zero
          dfpardVxb =zero
       endif
c
        dfdra = dfpardra+dfopdra
        dfdga = dfpardga+dfopdga
        dfdta = dfpardta+dfopdta
        dfdla = dfpardla+dfopdla
         dfdVxa = dfpardVxa+dfopdVxa
        dfdrb = dfpardrb+dfopdrb
        dfdgb = dfpardgb+dfopdgb
        dfdtb = dfpardtb+dfopdtb
        dfdlb = dfpardlb+dfopdlb
         dfdVxb = dfpardVxb+dfopdVxb
c        write(6,*)"emil for NDOP rhoA:", dfopdra
c        write(6,*)"emil for NDOP rhoB:", dfopdrb
c        write(6,*)"emil for NDOP GAA :", dfopdga
c        write(6,*)"emil for NDOP GBB :", dfopdgb
c        write(6,*)"emil for NDOP DTA :", dfopdta
c        write(6,*)"emil for NDOP DTB :", dfopdtb
c        write(6,*)"emil for NDOP DLA :", dfopdla
c        write(6,*)"emil for NDOP DLB :", dfopdlb
c        write(6,*)"emil for NDOP DUA :", dfopdVxa
c        write(6,*)"emil for NDOP DUB :", dfopdVxb
cccccccccccccccccccccccccccccccc
          D1F(i,ID_RA)  = dfdra
          D1F(i,ID_GAA) = dfdga
          D1F(i,ID_TA)  = dfdta
          D1F(i,ID_LA)  = dfdla
          D1FEXCH(i,1)  = dfdVxa
c
          D1F(i,ID_RB)  = dfdrb
          D1F(i,ID_GBB) = dfdgb
          D1F(i,ID_TB)  = dfdtb
          D1F(i,ID_LB)  = dfdlb
          D1FEXCH(i,2)  = dfdVxb
c
ccccccccccccccccccc
       else ! IMethod =1
cccccccccccccccccccccccccccccccccccc
c 188     Continue
          DNA13 = D1**third
          DNB13 = D2**third
          AXC = 1.24070098179880D0
          dfdra = -AXC*DNA13
          dfdga = zero
          dfdta = zero
          dfdla = zero

          dfdrb = -AXC*DNB13
          dfdgb = zero
          dfdtb = zero
          dfdlb = zero

             D1F(i,ID_RA)  = dfdra
             D1F(i,ID_GAA) = dfdga
             D1F(i,ID_TA)  = dfdta
             D1F(i,ID_LA)  = dfdla
             D1F(i,ID_RB)  = dfdrb
             D1F(i,ID_GAB) = 0.0d0
             D1F(i,ID_TB)  = dfdtb
             D1F(i,ID_LB)  = dfdlb
             D1F(i,ID_GBB) = dfdgb
          D1FEXCH(i,1)  = 0.0d0
          D1FEXCH(i,2)  = 0.0d0
cccccccccccccccccccccccccccccccc
        ENDIF  !IMethod
       ENDIF  ! RA,B
c       if (i.eq.5.or.i.eq.6) then
c             write(*,101) i, EX_HF_DENSA(i), fcor
c101           format(1x, 'iGrd', i3, 2f15.9)
c       endif
c          if(ndrv.eq.1) write(*,*) 'Fop = ', Fop,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'Fpar = ', Fpar,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'F(i) = ', F(i),'i = ',i
c          if(ndrv.eq.1) write(*,*) 'D1Fa = ', dfdra
c          if(ndrv.eq.1) write(*,*) 'D1Fb = ', dfdrb
c          if(ndrv.eq.1) write(*,*) 'DFdg = ', dfdga
c          if(ndrv.eq.1) write(*,*) 'DFdta = ', dfdta
c          if(ndrv.eq.1) write(*,*) 'DFdtb = ', dfdtb
c          if(ndrv.eq.1) write(*,*) 'DFdla = ', dfdla
c          if(ndrv.eq.1) write(*,*) 'DFdlb = ', dfdlb
c          if(ndrv.eq.1) write(*,*) 'DFdVxa = ', dfdVxa
c          if(ndrv.eq.1) write(*,*) ' dfopdra = ', dfopdra
c          if(ndrv.eq.1) write(*,*) ' dfpardra = ', dfpardra
c          if(ndrv.eq.1) write(*,*) ' dfadra = ', dfadra
c          if(ndrv.eq.1) write(*,*) ' dNdra = ', dNdra
c          if(ndrv.eq.1) write(*,*) ' dfopdrb = ', dfopdrb
c          if(ndrv.eq.1) write(*,*) 'fcor = ', fcor,'i = ',i
c          if(ndrv.eq.1) write(*,*) ' fc1 = ', fc1
c          if(ndrv.eq.1) write(*,*) ' fc2 = ', fc2
c          if(ndrv.eq.1) write(*,*) ' qqx = ', qqx
c          if((qqx.gt.0.10d0.and.qqx.lt.0.950d0).
c    $     and.qqx.ne.0.50d0) write(*,*) ' qqx = ', qqx,'i= ',i,
c    $     'iter  = ',iterSCF
c          if(eey1.gt.0.10d0.and.eey1.lt.0.950d0.
c    $     and.eey1.ne.0.50d0) write(*,*) ' eey1 = ', eey1,'i= ',i,
c    $     'iter  = ',iterSCF
c          if(eey2.gt.0.10d0.and.eey2.lt.0.950d0.
c    $     and.eey2.ne.0.50d0) write(*,*) ' eey2 = ', eey2,'i= ',i,
c    $     'iter  = ',iterSCF
c          if(ndrv.eq.1) write(*,*) ' gc1 = ', gc1
c          if(ndrv.eq.1) write(*,*) ' gc2 = ', gc2
c          if(ndrv.eq.1) write(*,*) ' dfpardrb = ', dfpardrb
c          if(ndrv.eq.1) write(*,*) ' M1A = ', SM1A
c          if(ndrv.eq.1) write(*,*) ' M2A = ', SM2A
c          if(ndrv.eq.1) write(*,*) ' VCPAR = ', VCPAR
c          if(ndrv.eq.1) write(*,*) ' ASA = ', ASA
c          if(ndrv.eq.1) write(*,*) ' ASB = ', ASB
c          if(ndrv.eq.1) write(*,*) ' ASSA = ', ASSA
c          if(ndrv.eq.1) write(*,*) ' ASSB = ', ASSB
c          if(ndrv.eq.1) write(*,*) ' dASSAdra = ', dASSAdra
c          if(ndrv.eq.1) write(*,*) ' dASSAdrb = ', dASSAdrb
c          if(ndrv.eq.1) write(*,*) ' dASAdra = ', dASAdra
c          if(ndrv.eq.1) write(*,*) ' dASAdrb = ', dASAdrb
c          if(ndrv.eq.1) write(*,*) ' DDA = ', DDA
c          if(ndrv.eq.1) write(*,*) ' DDB = ', DDB
c          if(ndrv.eq.1) write(*,*) ' EFNA = ', EFNA
c          if(ndrv.eq.1) write(*,*) ' EFNB = ', EFNB
c          if(ndrv.eq.1) write(*,*) ' dfcdraa = ', dfcdraa,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'LapA[i] = ', DLapA
c          if(ndrv.eq.1) write(*,*) 'LapB[i] = ', DLapB
c          if(ndrv.eq.1) write(*,*) 'tauA[i] = ', tauA
c          if(ndrv.eq.1) write(*,*) 'tauB[i] = ', tauB
c          if(ndrv.eq.1) write(*,*) 'rhoA_exact[i] = ', D1
c          if(ndrv.eq.1) write(*,*) 'rhoB[i] = ', D2
c          if(ndrv.eq.1) write(*,*) 'VXA[i] = ', VX1
c          if(ndrv.eq.1) write(*,*) 'gamB[i] = ', grb
 666      CONTINUE
        ENDDO  !  NGrid
       endif  ! iterSCF
c          write(*,*) 'ODDSUM =  ', ODDSUM(NCENTER), 'AIM  ', NCENTER 
c               IF(NCENTER.eq.NAtoms) then
c                 call VRtrace(sumsum, ODDSUM, NAtoms)
c          write(*,*) 'tot-sum =  ', sumsum, 'AIM  ', NCENTER 
c               endif
c
c      CALL ftnQFree(pEX_HF_DENSA)
c      CALL ftnQFree(pEX_HF_DENSB)
c      CALL ftnQFree(pFIT_DENSA)
c      CALL ftnQFree(pFIT_DENSB)
      RETURN
      END
c
      SUBROUTINE brsc(DN,DFN,GR,DLAP,tau,VX,UX,EFN,EFN2,dNdr,dNdg,dNdt,
     $ dNdl,
     $ dNdVx,dM1dr,dM2dr,dM1dt,dM2dt,dM1dl,dM2dl,dM1dg,dM2dg,dM1dVx,
     $ dM2dVx,dg3dn,dg3dt,dg3dg,SM1,SM2,DS,ndrv)
c
       IMPLICIT REAL*8(A-H,O-Z)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c      DN is the electron density for the input spin direction c
c      (i.e. either D1, or D2, but not the sum of the two);    c
c                                                              c
c      GR is the dot product of the  density gradient with     c
c      iself (per sin direction), NOT the grad.modulus;        c
c                                                              c
c      tau is the electron kinetic energy density for the      c
c      given spin direction                                    c
c                                                              c
c      DLAP is the density Laplacian                           c
c                                                              c
c      VX is the exact Slater exchange potential multiplied    c
c      by DN (for the given input spin direction)              c
c                                                              c
c     OUTPUT                                                   c
c     UX is the effective X Slater potential                   c
c                                                              c
c     EFN is the effective X hole normalizations               c
c                                                              c
c     EFN required to evaluate the functional derivativres     c
c     later on                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
       parameter(five=5.0d0,six=6.0d0,three=3.0d0)
       parameter(third=1.0d0/3.0d0)
       parameter(twthr=2.d0/3.d0, trihf=3.0d0/2.0d0)
       parameter(pi=3.141592653589793d0)
c
       parameter(aa1=0.9301841768238374D0,aa2=0.5485569431916153D0)
       parameter(ss1=0.9650036832623640d0)
       parameter(ss2=0.93018417682383740d0)
       parameter(ss3=0.54855694319161530d0)
       parameter(ss4=1.5158242410120670d0)
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
       parameter(d1 = 10.850190580831960d0  )
       parameter(d2 = 423.51564269756920d0  )
       parameter(d3 = 45.317325658682150d0  )
       parameter(d4 = 9.6414424525059270d0  )
       parameter(d5 = 1.9494258703373230d0  )
       parameter(d6 = 2.4210661522846020d0  )
       parameter(d7 = 2.1892374225965940d0  )
c
       parameter(e1 = 613.92073109477340d0     )
       parameter(e2 = 60.413933318777930d0     )
       parameter(e3 = 13.001906863306460d0     )
       parameter(e4 = 3.4757233899623390d0     )
       parameter(e5 = 2.5384493464526640d0     )
       parameter(e6 = 2.5561187910239140d0     )
       parameter(ee0 = 3.360517076724290D0     )
       parameter(ee1 = 4.623703278485152D0     )
       parameter(ee2 = 2.688949840405010D0     )
       parameter(ee3 = 0.6007166968496472D0    )
       parameter(ee4 = 0.03922040006408070D0   )
       parameter(ee5 = 0.0005438465669613952D0 )
       parameter(ee6 = 0.00000078437439010087D0)
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
       parameter(s0 =   8.382230306773296D0         )
       parameter(s1 =   19.60257290836815D0         )
       parameter(s2 =   19.71894106502952D0         )
       parameter(s3 =   10.77146542063455D0         )
       parameter(s4 =   3.415698370189622D0         )
       parameter(s5 =   0.5813153924917321D0        )
       parameter(s6 =   0.05426084060061605D0       )
       parameter(s7 =   0.002299629631716270D0      )
       parameter(s8 =   0.00005119354330427682D0    )
       parameter(s9 =   0.000000322977561012273D0   )
       parameter(s10 =  0.000000001405232963383258D0)
c
       parameter(a1 =  1.23767D0)
       parameter(a2 =  9.37587D0)
       parameter(a3 = 19.4777D0)
       parameter(a4 = 13.6816D0)
       parameter(a5 =  0.078655D0)
       parameter(a6 = 54.7264D0)
       parameter(a7 = 58.4331D0)
       parameter(a8 = 18.75174D0)
       parameter(a9 =  1.23767D0)
c
       parameter(dlam=1.0D0)
       parameter(AXC= 1.24070098179880D0)
       parameter(Tolbig=1.0d+10)
       parameter(dsmall=1.0D-13)
       parameter(dtol2=1.0D-08)
       parameter(dtol4=1.0D-06)
       parameter(dtol3=1.0D-08)
       parameter(detol0=1.0D-08)
       parameter(detol=1.0D-08,dtol=1.0D-08)
       parameter(delt=0.070d0)  ! smoothening of EFN

       REAL*8 xx1,xx2,yy,acc,X_y,EFN,DN,GR,DLAP,tau,SM1,SM2,DS
       INTEGER iter
ccccccccccccfor the numerical version onlyccccc
c      external FF
c      external rtsafe
cccccccccccccccccccccccccccccccccccccccccccccc
       Tol = 1.0D-08
c
           smth1 = 2.0d0 - delt
           smth2 = 2.0d0 + delt
c Initializations start
          dfdr=zero
          dfdt=zero
          dfdl=zero
          dfdVx = zero
          yy = dsmall
          X_y= dsmall
          g_y=zero
          P12_y=zero
          Q12_y=zero
          dxdy=zero
          dUdx=zero
          dUdr=zero
          dydn=zero
          dydg=zero
          dydt=zero
          dydl=zero
          dydVx=zero
          dQdn=zero
          dQdt=zero
          dQdl=zero
          dQdg=zero
          dNdr=zero
          dNdt=zero
          dNdl=zero
          dNdg=zero
          dNdVx=zero
          dM1dr = zero
          dM1dt = zero
          dM1dl = zero
          dM1dg = zero
          dM1dVx = zero
          dM2dr = zero
          dM2dt = zero
          dM2dl = zero
          dM2dg = zero
          dM2dVx = zero
          dg3dn = zero
          dg3dt = zero
          dg3dg = zero
          EFN = 1.0d0
          EFN1 = 1.0d0
          EFN2 = 1.0d0
          gc1=zero
          SM1 = 0.0d0
          SM2 = 0.0d0
          t1=zero
          t2=zero
          t3=zero
          t4=zero
          t5=zero
          t6=zero
          t7=zero
          t8=zero
          t9=zero
          t10=zero
          t11=zero
          t12=zero
          t13=zero
          t14=zero
          t15=zero
          t16=zero
          t17=zero
          tt9=zero
          ts14=zero
          ts15=zero
          ts16=zero
          ts23=zero
          tq1=zero
          tp6=zero
          tp7=zero
          t01=zero
          t02=zero
          tn1=zero
          tn2=zero
          tn3=zero
          tn4=zero
          tn5=zero
          tu1=zero
          tq1 = zero
          ndrv=2
c Initializations end
c          write(6,*)"in BRSC, emil's Rho: ",DN 
c          write(6,*)"in BRSC, emil's GRho: ",GR 
c          write(6,*)"in BRSC, emil's Tau: ",tau 
c          write(6,*)"in BRSC, emil's Lap: ",DLAP 
c          write(6,*)"in BRSC, emil's UX: ",UX 
c
c    First calculate X_y and the energy density
c
           DN13 = DN**third
           DN43= DN13*DN
           DN53= DN*(DN13*DN13)
           pi13=pi**third
           pi23 = pi13*pi13
c      IF(tau.gt.dsmall) then
           if(tau.lt.dsmall) tau=dsmall
c      IF(DN.gt.detol) then
           t1 = four*pi*DN
           tq1 = t1*DN
           DS = tau-GR/(DN*four)
           QS = DLAP/six -DS/three
C           write(6,*)"emil Q: ", QS
C           write(6,*)"emil UX in building Y: ", UX
c  To keep in mind that VX = Ux*DN:
          IF(tq1.gt.dtol2) then
             yy = -three*QS*UX/tq1
          endif
          if(abs(yy).lt.dsmall) ndrv =0
c            yy = -three*QS*UX/tq1
c       if(yy.gt.-Tolbig.and.yy.lt.Tolbig) then
cccccccccccfor the numerical version onlyccccccccccccccc
c           xx1=dtol
c           xx2 = 100.0d0
c        IF(yy.le.(-dsmall)) then
c           xx1 = detol
c           xx2 = 2.00000d0
c        elseif ((yy.ge.detol).and.(yy.lt.1000.0d0)) then
c           xx1 = 2.0000000d0
c           xx2 = 10.50d0
c        elseif ((yy.ge.1000.0d0).and.(yy.lt.Tolbig)) then
c           xx1 = 9.0d0
c           xx2 = 100.0d0
c        endif
c           acc = 1.0D-8
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           cc = two*pi13*DN43
           cc1=twthr*pi13*DN13
           cc2=two*pi13*DN13
cccccc set up for numerical evaluation of X_y cccccccccc
c          X_y = rtsafe(FF, yy, xx1, xx2, acc)
c            GO TO 115
cccccc end of set up for numerical eval. of X_y ccccccc
c        IF(yy.le.(-detol)) then
c           write(6,*)"emil Y: ", yy
         IF(yy.le.(-dsmall)) then
c        Region I : -infinity <= y <= 0
             t3 = atan(aa1*yy + aa2)
             t5 = atan(aa2)
             g_y = two*(two*t3 + pi)/(two*t5 + pi)
             t2 = yy*yy
             t4 = t2*yy
             t6 = t2*t2
             t8 = t6*yy
             tt9 = b0 +b1*yy +b2*t2 +b3*t4 +b4*t6 +b5*t8
             P12_y = dsmall
           IF(abs(tt9).gt.detol)
     $       P12_y = (c0+c1*yy+c2*t2+c3*t4+c4*t6+c5*t8)/tt9
             X_y=g_y*P12_y
c
             t5 = atan(aa2)
             tp1 = aa1*aa1
             tp7 = aa2*aa2
             dxx = 4.0d0*aa1/(one+tp1*t2+two*aa1*yy*aa2 + tp7)/
     $             (two*t5 + pi)
             ts14 = t6*t6
             ts15 = ts14*yy
             ts23 = tt9*tt9
             Q12_y = zero
           IF(ts23.gt.detol)
     $  Q12_y = (q0 +q1*yy +q2*t2 +q3*t4 +q4*t6 +q5*t8 +q6*t6*t2
     $        + q7*t6*t4 +q8*ts14 +q9*ts15)/ts23
             dxdy= dxx*P12_y + g_y*Q12_y
c
c          write(6,*)"switch 1"
c          write(6,*)"emil's X: ", X_y
c          write(6,*)"emil's DXDY: ", dxdy

         elseif((yy.ge.dsmall).AND.(yy.le.1000.0D0)) then
c          Region II : 0 < y < 1001
             t2 = yy*yy
             t4 = t2*yy
             t6 = t2*t2
             t8 = t6*yy
             t5 = t4*t4
             t9 = ee0 +ee1*yy +ee2*t2 +ee3*t4 +ee4*t6 +ee5*t8 +ee6*t5
             t10 = t2+e5*yy+e6
             X_y=dsmall
         IF(abs(t10).gt.detol)
     $      X_y = d1*(yy+d2)*(yy+d3)*(yy+d4)*(yy+d5)*(t2+d6*yy+d7)
     $         / (yy+e1)/(yy+e2)/(yy+e3)/(yy+e4)/t10
c
             ts14 = t6*t6
             ts15 = ts14*yy
             ts16 = t8*t8
             ts23 = t9*t9
              dxdy=zero
           IF(ts23.gt.detol)
     $       dxdy = (s0+s1*yy+s2*t2+s3*t4+s4*t6+s5*t8+s6*t6*t2
     $            + s7*t6*t4+s8*ts14+s9*ts15+s10*ts16)/ts23
c
c          write(6,*)"switch 3"
c          write(6,*)"emil's X: ", X_y
c          write(6,*)"emil's DXDY: ", dxdy

         elseif(yy.gt.1000.0D0) then
c            Region III : 1000 <= y
             t1 = log(yy)
             t2 = one/t1
c            t4 = t1**(one+t2)
             t22 = (t1+one)/t1
             t4 = t1**t22
             t6 = t1*t1
             t12 = t6*t6
        X_y = t4+a1/t1+a2/t6-a3/t6/t1+a4/t12-a5
             tp6 = log(t1)
             t7 = t1**t2
             dxdy1 = t7*(one - tp6/t1 + t2)/yy
             t8 = t2*t2/yy
             dxdy2 = t8*(-a9 - a8/t1 + a7/t6 - a6/t6/t1)
             dxdy = dxdy1+dxdy2

c             write(6,*)"switch 4"
c             write(6,*)"emil's X: ", X_y
c             write(6,*)"emil's DXDY: ", dxdy

        endif
            IF(abs(yy).lt.dsmall) X_y=2.000010d0
c           write(6,*)"emil X: ", X_y
c   end of calculating X(y)
c 115     CONTINUE   ! uncomment only for the numerical
c            twthr=2/3, trihf= 3/2
             t01 = sqrt(two)
             t02 = sqrt(three)
             tn0 = (twthr**trihf)*pi
             tn1 = DN*DN
             tn2 = sqrt(DN)
             tn4 = exp(X_y)
             tn5 = exp(-X_y)
             tu1 = X_y*X_y
             ts12 = 0.0d0
             tn10 = dsmall
             ts10 = dsmall
             xxxq = X_y*QS
c       IF((abs(X_y-2.0d0).gt.detol).and.X_y.gt.detol) then
        IF(X_y.gt.detol) then
             ts12 = X_y+4.0d0/X_y-(1.0d0+4.0d0/X_y)*tn5
          IF(abs(xxxq).gt.detol)
     $       tn10 = (X_y - two)*two/(xxxq*three)
           if(tn10.gt.dsmall) then
             ts10 = tn10*DN/4.0d0
c            tn11 = sqrt(abs(tn10))
             tn11 = sqrt(tn10)
c            ts11 = sqrt(abs(ts10))
             ts11 = sqrt(ts10)
             ALF3 = tn10*tn11
             EFN1 = ALF3*pi*tn4*DN*DN*sqrt(DN)
             EFN2 = EFN1
             IF(EFN2.gt.1.0d0) EFN2=1.0d0
           ggg= EFN1-2.0d0-delt
         IF((EFN1.gt.smth1).and.(EFN1.lt.smth2)) then
           EFN = 2.0d0-ggg*ggg/(4.0D0*delt)
           dNg1 = -ggg/(2.0D0*delt)
        elseif (EFN1.ge.smth2) then
           EFN = smth2
           dNg1 = 0.0d0
        else
            EFN = EFN1
           dNg1 = 1.0d0
         endif   ! EFN1
             SM1 = EFN*ts11*ts12
             SM2 = EFN*ts10*(tu1+12.0d0)
         else  !tn10
           EFN = 1.00d0
           EFN2=1.0d0
           SM1 = 0.0d0
           SM2 = 0.0d0
        endif   ! tn10
c        write(6,*)"emil N1: ", EFN1
c
c   Now prepare for the SCF potential:
c       IF(DN.gt.dtol2)  then
            t3 = DN*DN
            t4 = 9.0d0*QS+GR/(four*DN)
            dydn= t4*UX/(four*DN)/(pi*DN)/DN
            dQdn = -GR/(four*DN)/(three*DN)
            dydg = -UX/(DN*pi)/(four*DN)/(four*DN)
            dQdg = one/(DN*12.0D0)
            dydl = -UX/(pi*DN*two)/(four*DN)
            dydt = UX/(pi*DN)/(4.0D0*DN)
            dQdl = one/six
            dydVx = -three*QS/(four*DN)/(pi*DN)/DN
            dQdt = - third
            IF(DN.gt.dtol4) then
            dg3dn = -tau/DN/(three*DN)+GR/DN/(DN*six)/DN
            dg3dt = one/(three*DN)
            dg3dg = -one/(four*DN)/(three*DN)
            else
              dg3dn = 0.0d0
              dg3dt = 0.0d0
              dg3dg = 0.0d0
            endif
            dNdn = dNg1*5.0d0*EFN1/(two*DN)
            tt12 = (X_y-two)/(X_y*QS)
            dNdx = zero
c         write(6,*)"emil's DYDR: ", dydn
c         write(6,*)"emil's DYDG: ", dydg
c         write(6,*)"emil's DYDT: ", dydt
c         write(6,*)"emil's DYDL: ", dydl
c         write(6,*)"emil's DYDU: ", dydVx
C         write(6,*)"emil's A2 for rho: ", dg3dn
C         write(6,*)"emil's A2 for grho: ", dg3dg
C         write(6,*)"emil's A2 for tau: ", dg3dt
         if(tt12.gt.dsmall) then
            ts15 = sqrt(tt12)
c           ts15 = sqrt(abs(tt12))
           ts1 = sqrt(two)
            ts2 = sqrt(three)
            ts4 = DN*DN
            ts5 = sqrt(DN)
            ts9 = exp(X_y)
            ts17 = X_y*X_y
            dNdx = dNg1*(2.0d0/9.0d0)*ts1*ts2*pi*ts5*ts4*ts9*
     $      ts15*(ts17-2.0D0*X_y+3.0d0)/X_y/(X_y*QS)
          endif
            dNdQ = -dNg1*EFN1*three/(two*QS)
c
           t1 = DN*DN
           t6 = exp(X_y)
           t7 = X_y*X_y
           t12 = t7*t7
           t14 = QS*QS
           dM1dn = SM1*three/DN
           dM1dQ = -SM1*two/QS
           dM1dx = (two/9.0d0)*pi*t1*DN*(X_y-2.0d0)
     $     *(-t6*t7*X_y-12.0d0*t6*X_y+t6*t12+6.0d0*t6*t7+24.0d0*t6
     $     - 24.0d0)/t7/X_y/(QS*X_y)/QS
c
            t3 = sqrt(6.0d0)
            t4 = X_y - 2.0d0
            t5 = t4*t4
            t8 = X_y*X_y
            t10 = QS*QS
            tt13 = DN*t4/(X_y*QS)
            dM2dx = zero
           if(tt13.gt.(dtol)) then
            t18 = sqrt(tt13)
            t19 = exp(X_y)
            t23 = t8*t8
          dM2dx = pi*t1*DN*t3*t4/t8/X_y/QS
     $    *(t18*t19)/QS*(13.0d0*t8+60.0d0+t23-24.0d0*X_y)/27.0d0
        endif
           dM2dQ = -SM2*5.0d0/(two*QS)
           dM2dn=SM2*7.0d0/(two*DN)
c       else  !  dtol2
c         dydn = zero
c          dQdn = zero
c          dydg = zero
c          dQdg = zero
c          dydl = zero
c          dydt = zero
c          dQdl = zero
c          dydVx = zero
c          dQdt = zero
c          dg3dn = zero
c          dg3dt = zero
c          dg3dg = zero
c          dNdn = zero
c          dNdx = zero
c          dNdQ = zero
c          dM1dn=zero
c          dM1dx=zero
c          dM1dQ=zero
c          dM2dn=zero
c          dM2dx=zero
c          dM2dQ=zero
c      endif   ! dtol2
c    here dfdt == dNdt, etc.
         dNdr = dNdn + dNdx*dxdy*dydn +dNdQ*dQdn
         dNdt = dNdx*dxdy*dydt + dNdQ*dQdt
         dNdl = dNdx*dxdy*dydl+ dNdQ*dQdl
         dNdg = dNdx*dxdy*dydg + dNdQ*dQdg
         dNdVx = dNdx*dxdy*dydVx
c         write(6,*)"Q: ", QS
c         write(6,*)"DQDR: ", DQDn
c         write(6,*)"DQDG: ", DQDG
c         write(6,*)"DQDT: ", DQDT
c         write(6,*)"DQDL: ", DQDL
c         write(6,*)"X: ", X_y
c         write(6,*)"DXDR: ", dxdy*dydn
c         write(6,*)"DXDG: ", dxdy*dydg
c         write(6,*)"DXDT: ", dxdy*dydt
c         write(6,*)"DXDL: ", dxdy*dydl
c         write(6,*)"emil's N: ", EFN
c         write(6,*)"emil's N1R: ", DNDR
c         write(6,*)"emil's N1G: ", DNDG
c         write(6,*)"emil's N1T: ", DNDT
c         write(6,*)"emil's N1L: ", DNDL
c         write(6,*)"emil's N1U: ", DNDVX
         dM1dr = dM1dn + dM1dx*dxdy*dydn + dM1dQ*dQdn
         dM2dr = dM2dn + dM2dx*dxdy*dydn +dM2dQ*dQdn
         dM1dt = dM1dx*dxdy*dydt + dM1dQ*dQdt
         dM2dt = dM2dx*dxdy*dydt + dM2dQ*dQdt
         dM1dg = dM1dx*dxdy*dydg + dM1dQ*dQdg
         dM2dg = dM2dx*dxdy*dydg + dM2dQ*dQdg
         dM1dl = dM1dx*dxdy*dydl + dM1dQ*dQdl
         dM2dl = dM2dx*dxdy*dydl + dM2dQ*dQdl
         dM1dVx = dM1dx*dxdy*dydVx
         dM2dVx = dM2dx*dxdy*dydVx
C         write(6,*)"emil's M1: ", SM1
C         write(6,*)"emil's DM1DX: ", dM1dx
C         write(6,*)"emil's DM1DR: ", DM1DR
C         write(6,*)"emil's DM1DG: ", DM1DG
C         write(6,*)"emil's DM1DT: ", DM1DT
C         write(6,*)"emil's DM1DL: ", DM1DL
C         write(6,*)"emil's DM1DU: ", DM1DVX
C         write(6,*)"emil's M2: ", SM2
C         write(6,*)"emil's DM2DR: ", DM2DR
C         write(6,*)"emil's DM2DG: ", DM2DG
C         write(6,*)"emil's DM2DT: ", DM2DT
C         write(6,*)"emil's DM2DL: ", DM2DL
C         write(6,*)"emil's DM2DU: ", DM2DVX
       else   ! abs(X_y-2.0d0)
          dNdr = zero
          dNdt = zero
          dNdl = zero
          dNdg = zero
         dNdVx = zero
         dM1dr = zero
         dM2dr = zero
         dM1dt = zero
         dM2dt = zero
         dM1dg = zero
         dM2dg = zero
         dM1dl = zero
         dM2dl = zero
        dM1dVx = zero
        dM2dVx = zero
          EFN = 1.0d0
       endif  ! X_y
c      endif  ! Tolbig yy
c     ENDIF  ! DN and  tau
      RETURN
      END
ccccccccconly for the numerical version only ccccc
c     SUBROUTINE FF(X,FX,DFX,YY)
c     IMPLICIT REAL*8(A-H,O-Z)
c      REAL*8 X, FX, DFX, YY
c       t5 = exp(X)
c       t7 = X**2
c       t8 = t7*X
c       t2 = t5 - 1.0d0 -X*0.50d0
c       t3 = t2*(X-2.0d0)/(X*X)
c       FX = t3-YY
c       DFX = t5/X -0.3D1*t5/t7+(0.4D1*t5-0.4D1)/t8
c     return
c     end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
