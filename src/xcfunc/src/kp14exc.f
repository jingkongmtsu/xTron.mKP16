      SUBROUTINE kp14exc(iFunc,F,D1F,D1FEXCH,W,RA,RB,D1RA,D1RB,RLapA,
     $ RLapB,TA,TB,EX_HF_DENSA,EX_HF_DENSB,ODEL,ODELW,ODDSUM,
     $NGrid,iterSCF,ICEN,NA,NB)
c    ******************************************************************
c    *                                                                *
c    *  evaluates the non-dynamic correlation energy by the           *
c    *  real-space correlation model B13 of Becke with RI.            *
c    *                                                                *
c    *  reference: A. D. Becke, J. Chem. Phys. 122, 64101 (2005);     *
c    *  reference: A. D. Becke, J. Chem. Phys. 138, 074109 (2013);    *
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
        REAL*8 YSUM(1500)
        REAL*8 PKSUM(1500)
        REAL*8 AVNEF(1500)
        REAL*8 fc1,fc2,UCA,UCB,DSA,DSB,DDA,DDB,ASA,ASB,ASSA,ASSB,sumsum
        REAL*8 EFNA,EFNB,UX1,UX2,tauA,tauB,VX1,VX2,DLapA,DLapB

        REAL*8 EX_HF_DENSA(NGrid)
        REAL*8 EX_HF_DENSB(NGrid)

        ! fenglai: this is used to turn on/off KP14 correlation
        REAL*8 kp14fac
        INTEGER nDrvA,nDrvB

C        REAL*8 FIT_DENSA(1)
C        REAL*8 FIT_DENSB(1)

       INTEGER NGrid,iterSCF,ICEN,NCENTER
       INTEGER NUMVAR
       INTEGER icounttt, icounttp
       INTEGER POS_RA,POS_RB,POS_GAA,POS_GAB,POS_GBB
       INTEGER POS_TA,POS_TB,POS_LA,POS_LB
         PARAMETER(POS_RA = 1)
         PARAMETER(POS_RB = 2)
         PARAMETER(POS_GAA = 3)
         PARAMETER(POS_GAB = 4)
         PARAMETER(POS_GBB = 5)
         PARAMETER(POS_TA = 6)
         PARAMETER(POS_TB = 7)
         PARAMETER(POS_LA = 8)
         PARAMETER(POS_LB = 9)
         PARAMETER(POS_EXA = 10)
         PARAMETER(POS_EXB = 11)
c
C       POINTER (pEX_HF_DENSA,EX_HF_DENSA)
C        POINTER (pEX_HF_DENSB,EX_HF_DENSB)
c
C        POINTER (pFIT_DENSA,FIT_DENSA)
C        POINTER (pFIT_DENSB,FIT_DENSB)
        parameter(third=1.0d0/3.0d0,zero=0.0D0)
        parameter(fifthr= 5.0d0/3.0d0)
        parameter(twothr= 2.0d0/3.0d0)
        parameter(twthr= 2.0d0/3.0d0)
        parameter(dsmall=1.0D-12)
        parameter(dsmall3=1.0D-09)
        parameter(dsmall2=1.0D-14)
        parameter(detol=1.0D-08)
        parameter(dtol=5.0D-07)
        parameter(delt=0.050d0)  ! smoothening of fa, fb
        parameter(delt2=0.050d0)  !
        parameter(delt3=0.00000000050d0)  !
        parameter(smoth=115.0D0) ! unifies the definition of fcor
        parameter(smoth2=120.0D0)
c       parameter(aaa= 0.20D0)
        SAVE NCENTER
        DATA NCENTER /-1/
C       PI = ConvFac(VALUE_OF_PI)
        PI = 3.14159255358979323846D0
c       Tolbig = HUGE(0.0d0)*0.850d0
c        Tolbig = 1.0d+30
         Tolbig = 1.0d+14
        Tolbig2 = 1.0d+10
        smth1 = 1.0d0-delt
        smth2 = 1.0d0+delt
        smth3 = -delt3
        smth4 =  delt3

C       fenglai: so far we only use file to tranfer the D1F for exchange
C       variable, so define the number of variables not including
C       exchange
        NUMVAR = 9           

        ndrv = 5
        iqlamb = 10
        qlamb = 1.0d0
        ACOP = 1.0d0
        icorr = 5
        alf = 0.10d0
c       alf = 0.1590d0
        ialf= 0
        ialf= 0
c
c        CALL RemGet(icorr,5907)
        icorr = 5
        kp14fac = 1.0D0
c
c       IF(icorr.eq.2) then
c         CALL RemGet(iqlamb,5906)
c         qlamb = float(iqlamb)/10.0d0
c       endif
c       write(*,*) 'qlamb = ', qlamb, 'icorr = ', icorr
c          Bkp = 0.6990d0
c         ACOP  = 0.590D0
c         ACPAR = 1.0D0
c         DCOP = 0.50D0
c         DCPAR= 0.50D0
c         alf= 0.420d0
c
c        Yihan and Emil, B05 adjustable parameters
c         IF(ICEN.eq.1) ndrv=4
c        CALL RemGet(maxiter,REM_MAXSCF)
ccccc Parameters for KP14(stand,nofrac) with HF initial guess:
         if(icorr.eq.1) then
          Bkp = 0.6990d0
          ACOP  = 1.0d0
          ACPAR = 0.4740D0
          DCOP = 0.66D0
          DCPAR= 0.615D0
          alf= 0.420d0
         endif
ccccc Parameters for KP14(stand,nofrac) with LSD initial guess:
c         Bkp = 0.740d0
c         Ckp = 0.0d0
c         ACOP  = 1.0d0
c         ACPAR = 0.5120D0
c         DCOP = 0.6550D0
c         DCPAR= 0.615D0
c       endif
ccccc Parameters for KP14-YY(nofrac) with HF initial guess:
c         Bkp = 0.70d0
c         ACOP  = 1.0d0
c         ACPAR = 0.550D0
c         DCOP = 0.6590D0
c         DCPAR= 0.615D0
c cccc Parameters for KP14(stand,nofrac) combined with La averaged Ec(dyn):
        IF(icorr.eq.3) then
          Bkp = 1.810d0
          ACOP  = 1.0d0
          ACPAR = 0.50D0
          DCOP = 1.040D0
          DCPAR= 1.07D0
          alf= 0.1590d0
        endif
c
c cccc Parameters for KP14(frac2), only when icorr = 4
        IF(icorr.eq.4) then
          Bkp = 1.20d0
c         Bkp = 1.50d0
          ACOP  = 1.0d0
          ACPAR = 1.0D0
          DCOP = 0.50D0
          DCPAR= 0.50D0
          alf= 0.1590d0
c
c         Bkp = 0.7200d0
c         ACOP  = 1.0d0
c         ACPAR = 0.7730D0
c         DCOP = 0.590D0
c         DCPAR= 0.6150D0
c         alf = 0.420d0
        endif
        IF(icorr.eq.5) then
          kp14fac = 1.0D0
          Bkp = 1.355d0
          alf= 0.038d0
          ACOP  = 1.0d0
          ACPAR = 1.128D0
          DCOP = 1.0D0
          DCPAR= 1.0D0
        endif
c
         if(icorr.eq.6) then
          Bkp = 0.6990d0
          ACOP  = 1.0d0
          ACPAR = 1.0D0
          DCOP = 0.66D0
          DCPAR= 0.615D0
          alf= 0.420d0
         endif

c        fenglai: now let's choose functional
         if (ifunc .eq. 1) THEN
            DCOP = 0.0D0
            DCPAR= 0.0D0
            kp14fac = 1.0D0
         else if (ifunc .eq. 2) THEN
            kp14fac = 0.0D0
            DCPAR = 0.0D0
            DCOP = 1.0D0
         else if (ifunc .eq. 3) THEN
            kp14fac = 0.0D0
            DCOP = 0.0D0
            DCPAR = 1.0D0
         end if

c       write(*,*) 'icorr = ', icorr 
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 0, FM_BEG, Bkp2)
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 1, FM_BEG, ACPAR2)
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 2, FM_BEG, DCOP2)
c        CAll FileMan(FM_READ,FILE_B05_PARAMETERS,FM_DP, 1, 3, FM_BEG, DCPAR2)
c        IF (Bkp2.GT.-98.0) THEN
c          Bkp = Bkp2 
c        ENDIF
c        IF (ACPAR2.GT.-98.0) THEN
c           ACPAR = ACPAR2
c        ENDIF
c        IF (DCOP2.GT.-98.0) THEN
c           DCOP = DCOP2
c        ENDIF
c        IF (DCPAR2.GT.-98.0) THEN
c           DCPAR = DCPAR2
c        ENDIF
c       write(*,*) 'BKP = ', BKP 
c       write(*,*) 'ACPAR  = ', ACPAR 
c       write(*,*) 'DCPAR  = ', DCPAR 
c       write(*,*) 'DCOP  = ', DCOP 
c        IF((icorr.eq.4).or.(icorr.eq.5).or.(icorr.eq.6)) then
c          CALL RemGet(ialf,5908)
c           alf = float(ialf)/1000.0d0
c         endif
c       write(*,*) 'alf = ', alf

c        CALL RemGet(NA,REM_NALPHA)
c        CALL RemGet(NB,REM_NBETA)
c        CALL RemGet(IUnres,REM_JUSTAL)
c        CALL RemGet(NAtoms,REM_NATOMS)
c       write(*,*)'IUnres=  ',IUnres
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
c        CALL VRLoad(YSUM,NAtoms,0.0d0)
c        CALL VRLoad(PKSUM,NAtoms,0.0d0)
c        CALL VRLoad(AVNEF,NAtoms,0.0d0)
c
c        CALL ftnQAllocREAL8(pFIT_DENSA,NGrid)
c        CALL ftnQAllocREAL8(pFIT_DENSB,NGrid)
c        CALL VRLoad(FIT_DENSA,NGrid,0.0d0)
c        CALL VRLoad(FIT_DENSB,NGrid,0.0d0)

c        CALL RemGet(IMethod, REM_HF_EXCHANGE_WITH_GRID)
C       so far we only support the RI approximation        
         IMethod = 2
c         Using RI-approximated exact exchange

c          CALL FileMan(FM_READ,FILE_EX_HF_DENSITY,FM_DP,NGrid,
c     $                 0,FM_BEG,EX_HF_DENSA)
c          CALL FileMan(FM_READ,FILE_FITTED_DENSITY,FM_DP,NGrid,
c     $                 0,FM_BEG,FIT_DENSA)
C          call prtmat (EX_HF_DENSA, NGrid, 1, NGrid, 1,  6, 
C     $    'EX_HF_DENSA')

c          IF (NMO.eq.2) THEN
c             CALL FileMan(FM_READ,FILE_EX_HF_DENSITY,FM_DP,NGrid,
c     $             NGrid,FM_BEG,EX_HF_DENSB)
c            CALL FileMan(FM_READ,FILE_FITTED_DENSITY,FM_DP,NGrid,
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
         ffcor = 0.0d0
         fffcor = 0.0d0
         UCA = zero
         UCB = zero
         DSA =dsmall
         DSB =dsmall
         DDA = dsmall
         DFA = dsmall
         DFF = dsmall
         DFF0 = zero
         DDB = dsmall
         DFB = dsmall
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
         EFNA = 0.0d0
         EFN2A = 0.0d0
         EFNB = 0.0d0
         EFN2B = 0.0d0
         UX1=-zero
         UX2=-zero
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
         eshift = 0.0d0
         S1 = 0.0d0
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
         dg4dra = zero
         dg4dta = zero
         dg4dga = zero
         dg4drb = zero
         dg4dtb = zero
         dg4dgb = zero
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
         dfpar2dVxa = zero
         dfpardra = zero
         dfpar2dra = zero
         dfpardga = zero
         dfpar2dga = zero
         dfpardla = zero
         dfpar2dla = zero
         dfpardta = zero
         dfpar2dta = zero
         dfpardVxa = zero
         dfpar2dVxa = zero
         dfpardVxb = zero
         dfpar2dVxb = zero
         dfpardrb = zero
         dfpar2drb = zero
         dfpardgb = zero
         dfpar2dgb = zero
         dfpardlb = zero
         dfpar2dlb = zero
         dfpardtb = zero
         dfpar2dtb = zero
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
c         Zkp = Tolbig2
         Zkp = dsmall
         ZZkp = dsmall 
         PBkp = 0.5900d0
         dPdz = 0.0d0
         dPkpdra = 0.0d0
         dPkpdga = 0.0d0
         dPkpdta = 0.0d0
         dPkpdla = 0.0d0
         dPkpdVxa =0.0d0
         dEkpdra = 0.0d0
         dEkpdga = 0.0d0
         dEkpdla = 0.0d0
         dEkpdta = 0.0d0
         dEkpdVxa = 0.0d0 
         dPkpdrb = 0.0d0
         dPkpdgb = 0.0d0
         dPkpdtb = 0.0d0
         dPkpdlb = 0.0d0
         dPkpdVxb =0.0d0
         dEkpdrb = 0.0d0
         dEkpdgb = 0.0d0
         dEkpdlb = 0.0d0
         dEkpdtb = 0.0d0
         dEkpdVxb = 0.0d0
         dZdra = 0.0d0
         dZdga = 0.0d0
         dZdta = 0.0d0
         dZdla = 0.0d0
         dZdVxa = 0.0d0
         dZdrb = 0.0d0
         dZdgb = 0.0d0
         dZdtb = 0.0d0
         dZdlb = 0.0d0
         dZdVxb = 0.0d0
         ddffdra = 0.0d0
         ddffdga = 0.0d0
         ddffdta = 0.0d0
         ddffdrb = 0.0d0
         ddffdgb = 0.0d0
         ddffdtb = 0.0d0
c Ulamb is the lambda dependent prefactor without Und multiplied
         Ulamb = 1.0d0*qlamb 
         Ukp = 0.0d0 
         Ukp2 = 0.0d0 
         Ukp3 = 0.0d0 
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
c      write(*,*) 'NGrid for this batch = ', NGrid 
          Ebatch = 0.0d0
          icountUD = 0 
          icountUnd = 0 
          icountVX = 0
          icountNA = 0
          icountNB = 0
          icountD = 0
          icountNeff =0
          icounttt = 0
          icounttp = 0
c   LOOP OVER GRID STARTS
       DO i = 1,NGrid
c          write(*,*) 'Grid point No', i

         ! fenglai debug
c         write(6,*)"     "
c         write(6,*)"----------------------"
c         write(6,*)"for grid: ", i
          
        IF(NMO.eq.1) then
           EX_HF_DENSB(i)=EX_HF_DENSA(i)
c          RB(i) = RA(i)
c          TB(i) = TA(i)
         endif
          D1 = RA(i)
          D2 = RB(i)
          ds1 = abs(D1-D2)
        IF (IMethod.eq.2) THEN
c        IF(NMO.eq.1) FIT_DENSB(i)=FIT_DENSA(i)
c          DF1 = FIT_DENSA(i)
c          DF2 = FIT_DENSB(i)
        endif
c
cccccccccccccccccccc
c       IF((i.eq.8).or. (i.eq.15).or.(i.eq.20)) then

        ! fenglai: have a question about the ndrv??
c        IF(i.eq.10) then
c          ndrv = 1
c  N
c        D1 = 97.61423963505630d0 +0.050d0
c        D1 = 97.61423963505630d0 -0.050d0
c
c  H
c
c        D1 = 0.066093251507263210d0 + 0.000010d0
c        D1 = 0.066093251507263210d0 - 0.000010d0

c ch
c       D1 = 0.2002989991863230d0  + 0.000010d0
c       D1 = 0.2002989991863230d0  - 0.000010d0
c       endif
ccccccccccccccccccc
         fc1 = zero
         eshift = 0.0d0
         fc2 = zero
         fcor = zero
         ffcor = 0.0d0
         fffcor = 0.0d0
         PBkp = 0.59d0
         dPdz = 0.0d0
         Ulamb= qlamb 
         UCA = zero
         UCB = zero
         DSA =dsmall
         DSB =dsmall
         DDA = dsmall
         DFA = dsmall
         DFF = dsmall
         DFF0 = zero
         DDB = dsmall
         DFB = dsmall
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
         EFNA = 0.0d0
         EFN2A = 0.0d0
         EFNB = 0.0d0
         EFN2B = 0.0d0
         UDYN = zero
         ULADYN = zero
         EDYN = zero
         EAVDYN = 0.0d0
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
         dfpar2dVxa = zero
         dfpardra = zero
         dfpar2dra = zero
         dfpardga = zero
         dfpar2dga = zero
         dfpardla = zero
         dfpar2dla = zero
         dfpardta = zero
         dfpar2dta = zero
         dfpardVxa = zero
         dfpar2dVxa = zero
         dfpardVxb = zero
         dfpar2dVxb = zero
         dfpardrb = zero
         dfpar2drb = zero
         dfpardgb = zero
         dfpar2dgb = zero
         dfpardlb = zero
         dfpar2dlb = zero
         dfpardtb = zero
         dfpar2dtb = zero
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
         Ukp =  0.0d0 
         Ukp2 = 0.0d0 
         Ukp3 = 0.0d0 
         UUkp  = zero
c         Zkp = Tolbig2
         Zkp = dsmall
         ZZkp = dsmall
         Ekp14 = zero
         End0 = 0.0d0
         Elamb14 = zero
         Fpar = zero
         ndrvA = 2
         ndrvB = 2
c
             UX1 = -dsmall
             UX2 = -dsmall
             VX1=EX_HF_DENSA(i)
             VX2=EX_HF_DENSB(i)
c
         IF((VX1.gt.-dsmall).or.(VX2.gt.-dsmall)) GO TO 666
c        IF((D1.le.detol).OR.(D2.le.detol)) icountD = icountD + 1
         IF((D1.gt.detol).OR.(D2.gt.detol)) then

         IF(D1.gt.detol) then
              UX1 = VX1/D1
cccccccccccccccccccccccccccc
c         if(ndrv.eq.1) then
c  n2
c             VX1 = -599.5352554045510d0 + 0.001d0
c             VX1 = -599.5352554045510d0 - 0.001d0
c             UX1= VX1/D1

c  ch
c             VX1 = -0.1785385767938490d0 + 0.00001d0
c             VX1 = -0.1785385767938490d0 - 0.00001d0
c             UX1= VX1/D1
c         endif
cccccccccccccccc
               gra = D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
               tauA = TA(i)
               DLapA = RLapA(i)

               ! fenglai debug
c               write(6,*)"emil RhoA: ",D1
c               write(6,*)"emil GRhoA: ",gra
c               write(6,*)"emil TA: ",TauA
c               write(6,*)"emil LA: ",DLapA 
c               write(6,*)"emil UA: ",VX1
cccccccccccccc
c  ch
c         if(ndrv.eq.1) gra = 0.013861531360698990d0 + 0.000001d0
c         if(ndrv.eq.1) gra = 0.013861531360698990d0 - 0.000001d0

c n2
c         if(ndrv.eq.1) gra = 18.61259885126080d0 + 0.0010d0
c         if(ndrv.eq.1) gra = 18.61259885126080d0 - 0.0010d0

c  ch
c         if(ndrv.eq.1) tauA = 0.029284924703833880d0+0.000001d0
c         if(ndrv.eq.1) tauA = 0.029284924703833880d0-0.000001d0

c n2
c         if(ndrv.eq.1) tauA =  22.43602181678550d0 + 0.01d0
c         if(ndrv.eq.1) tauA =  22.43602181678550d0 - 0.01d0

c ch
c         if(ndrv.eq.1) DLapA =  -12.09368923688460d0 + 0.005d0
c         if(ndrv.eq.1) DLapA =  -12.09368923688460d0 - 0.005d0

c n2
c         if(ndrv.eq.1) DLapA =  -609401.9274966730d0 + 0.5d0
c         if(ndrv.eq.1) DLapA =  -609401.9274966730d0 - 0.5d0
cccccccccccccc
c brsc is the routine doing the main RSC numerics
           CALL brsc14(D1,DF1,gra,DLapA,tauA,VX1,UX1,EFNA,EFN2A,dNdra,
     $        dNdga,dNdta,dNdla,dNdVxa,dM1dra,dM2dra,dM1dta,dM2dta,
     $        dM1dla,dM2dla,dM1dga,dM2dga,dM1dVxa,dM2dVxa,dg3dra,
     $        dg3dta,dg3dga,dg4dra,dg4dta,dg4dga,SM1A,SM2A,DSA,detol,
     $        dsmall,ndrvA)
       else  ! D1
               SM1A = zero
               SM2A = zero
               DSA = dsmall
               gra = 0.0d0
               tauA = zero
               DLapA = zero
               UX1 = -dsmall
               EFNA = 0.00d0
               EFN1A = 0.00d0
               EFN2A = 0.00d0
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
               dg4dra=zero
               dg3dra=zero
               dg3dta=zero
               dg4dta=zero
               dg3dga=zero
               dg4dga=zero
       endif ! D1, VX1
c
        IF(D2.gt.detol) then
              UX2 = VX2/D2
cccccccccccccccccccccccccccccccc
c         if(ndrv.eq.1) then
c             VX2 = -72.75196487421990d0 + 0.005d0
c             VX2 = -72.75196487421990d0 - 0.005d0
c             UX2= VX2/D2
c         endif
               grb = D1RB(i,1)**2 + D1RB(i,2)**2 + D1RB(i,3)**2
               tauB = TB(i)
               DLapB = RLapB(i)

               ! fenglai debug
c               write(6,*)"emil RhoB: ",D2
c               write(6,*)"emil GRhoB: ",grb 
c               write(6,*)"emil TB: ",TauB
c               write(6,*)"emil LB: ",DLapB 
c               write(6,*)"emil UB: ",VX2
               
         CALL brsc14(D2,DF2,grb,DLapB,tauB,VX2,UX2,EFNB,EFN2B,dNdrb,
     $        dNdgb,dNdtb,dNdlb,dNdVxb,dM1drb,dM2drb,dM1dtb,dM2dtb,
     $        dM1dlb,dM2dlb,dM1dgb,dM2dgb,dM1dVxb,dM2dVxb,dg3drb,
     $        dg3dtb,dg3dgb,dg4drb,dg4dtb,dg4dgb,SM1B,SM2B,DSB,detol,
     $        dsmall,ndrvB)

           else  !  D2
               SM1B= zero
               SM2B= zero
               DSB= dsmall
               grb = zero
               tauB = zero
               DLapB = zero
               UX2 = -dsmall
               EFNB = 0.00d0
               EFN1B = 0.00d0
               EFN2B = 0.00d0
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
               dg4drb=zero
               dg4dtb=zero
               dg3dgb=zero
               dg4dgb=zero
       endif  ! VX2
c
                   fcor = 0.0d0
                   ffcor =0.0d0
                   ggcorr = zero 
                   fc1 = zero
                   fc2=zero
                   gcc1 = zero
                   gcc2 = zero
                   gc1=zero
                   gc2=zero
                   Fop = zero 
           IF((ndrvA.eq.0).or.(ndrvB.eq.0)) icountNeff = icountNeff+1
           IF((ndrvA.eq.0).or.(ndrvB.eq.0)) GO TO 777
         if((D1.gt.detol.and.ndrvA.eq.2).and.(D2.gt.detol
     $      .and.ndrvB.eq.2)) then
c
         IF(EFN2B.gt.detol) gcc1 = (1.0d0-EFN2A)/EFN2B
         IF(EFN2A.gt.detol) gcc2 = (1.0d0-EFN2B)/EFN2A
         IF(EFNB.gt.dsmall) then
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
        IF(EFNA.le.detol) icountNA = icountNA+1   
       IF(EFNA.gt.dsmall) then
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
c         write(6,*)"emil's alpha Neff", EFNA
c         write(6,*)"emil's beta  Neff", EFNB

           tx = 1.0d0
           qqx = 0.50d0
           smtx = 0.0d0
           xxx=zero
c        fcor = fc1*0.50d0 + fc2*0.50d0
        if((abs(fc1).gt.dsmall).or.(abs(fc2).gt.dsmall)) then
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
c          if(gcc1.lt.dsmall) gcc1=0.0d0
c          if(gcc2.lt.dsmall) gcc2=0.0d0
        else  !  fc1,fc2
           xxx=zero
           qqx = 0.0d0
           tx = 0.0d0
           smtx=zero
          fcor = zero
         endif  !  fc1,fc2
c         write(6,*)"emil's f1 factor", fc1
c         write(6,*)"emil's f2 factor", fc2
c         write(6,*)"emil's qqx factor", qqx
c         write(6,*)"emil's f factor", fcor
ccccccccccccccccccccc
c        fffcor = min(gcc1,gcc2,1.0d0)
         fffcor = min(gc1,gc2,1.0d0)
c        delt1 = 1.0d0 -aaa
c        delt2 = 1.0d0 +aaa
c        IF((ffcor.gt.delt1).and.(ffcor.lt.delt2)) then
c        IF((fffcor.gt.delt1).and.(fffcor.le.1.0d0)) then
c        ffcor = 3.0d0*fffcor*fffcor-2.0d0*(fffcor**3)
c        else
           ffcor = fffcor
c        endif
         ggcorr = min(2-gcc1-gcc2,1.0d0)
          ttdA = zero
          ttdB = zero
          ttgg = zero
          ttg1 = zero
          ttg2 = zero
c         ttdB = UX2*D1/D2
c         ttg2=D1*VX2/D2
c         ttdA = UX1*D2/D1
c         ttg1=D2*VX1/D1
cccccccccccccccccccccccccccccccc
        IF(D2.gt.detol) then
           ttdB = UX2*D1/D2
           ttg1=D2*UX1
         endif
         IF(D1.gt.detol) then
           ttdA = UX1*D2/D1
           ttg2=D1*UX2
         endif
          ttgg = ttg1+ttg2
c        Fop=ACOP*0.50d0*ffcor*ttgg
         Fop=0.50d0*fcor*ttgg
         Ukp = 0.50d0*fcor*ttgg
         UUkp = 0.50d0*ttgg
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
           ffcor = zero
           tx = 0.0d0
          Fop = 0.0d0
          Ukp = zero
          UUkp = zero
          ttdA=0.0d0
          ttdB=0.0d0
      endif 
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
            Ukp2 = 0.0d0
            Ukp3 = 0.0d0
c
              ASSA =0.0d0
              A2SA = 0.0d0
              AASA = 0.0d0
      IF((D1.gt.detol).and.(ndrvA.eq.2)) then
          IF(abs(SM2A).gt.detol) then
           ASA=(1.0d0-EFNA-fcor*EFNB)/SM2A
           A2SA=(1.0d0-EFN2A-fcor*EFN2B)/SM2A
               endif
            DDA = DSA/(3.0d0*D1)
c        IF(D1.gt.dtol.and.NA.gt.1) then
        IF(D1.gt.dtol.and.NA.gt.1.and.DSA.gt.dsmall3) then
            D113 = D1**third
            D153 = D1**fifthr
            DFA = DSA/D153
            D159 = D153**third
             DSA13 = (DSA)**third
            DFAA = DSA13/D159
             D123= D1**twthr
             D129= D123**third
            D1149=D159*D159*D129*D129
              DSA23=DSA13*DSA13
            ddffdra = dg4dra/(DSA23*3.0d0*D159)
     $              - 5.0d0*DSA13/(9.0d0*D1149)
            ddffdga = dg4dga/(DSA23*3.0d0*D159)
            ddffdta = dg4dta/(DSA23*3.0d0*D159)
        else
             DFAA = dsmall
             DFA= dsmall
            ddffdra = 0.0d0
            ddffdga = 0.0d0
            ddffdta = 0.0d0
        endif
c         write(6,*)"emil DA_RA: ", ddffdra
c         write(6,*)"emil DA_GAA: ",ddffdga
c         write(6,*)"emil DA_TA: ", ddffdta
            yy1 = ASA-DDA
            smtx2 = smoth2*yy1
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
c           ASSA = ASA*eey1+DDA*(1.0d0-eey1)
            ASSA = ASA*eey1+DDA-DDA*eey1
c           ASSA = min(ASA,DDA) 
            AASA = min(A2SA,DDA)
        endif  ! D1
ccccccccccccccccccccc
              ASSB=0.0d0
              A2SB = 0.0d0
              BBSB = 0.0d0
       IF((D2.gt.detol).and.(ndrvB.eq.2)) then
         IF(abs(SM2B).gt.detol) then
            ASB=(1.0d0-EFNB-fcor*EFNA)/SM2B
            A2SB=(1.0d0-EFN2B-fcor*EFN2A)/SM2B
               endif
            DDB = DSB/(3.0d0*D2)
c         IF(D2.gt.dtol.and.NB.gt.1) then
          IF(D2.gt.dtol.and.NB.gt.1.and.DSB.gt.dsmall3) then
            D213 = D2**third
            D253 = D2**fifthr
            D259 = D253**third
            DFB = DSB/D253
             DSB13 = (DSB)**third
            DFBB = DSB13/D259
            D223=D2**twothr
            D229 = D223**third
            D2149 = D259*D259*D229*D229
              DSB23=DSB13*DSB13
            ddffdrb = dg4drb/(DSB23*3.0d0*D259)
     $              - 5.0d0*DSB13/(9.0d0*D2149)
            ddffdgb = dg4dgb/(DSB23*3.0d0*D259)
            ddffdtb = dg4dtb/(DSB23*3.0d0*D259)
          else
              DFBB = dsmall
              DFB = dsmall
            ddffdrb = 0.0d0
            ddffdgb = 0.0d0
            ddffdtb = 0.0d0
          endif
c         write(6,*)"emil DA_RB: ", ddffdrb
c         write(6,*)"emil DA_GBB: ",ddffdgb
c         write(6,*)"emil DA_TB: ", ddffdtb
c
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
c           ASSB = ASB*eey2+DDB*(1.0d0-eey2)
            ASSB = ASB*eey2+DDB-DDB*eey2
c           ASSB = min(ASB,DDB) 
            BBSB = min(A2SB,DDB)
        endif   ! D2
ccccccccc
           UCA=0.0d0
           UCB=0.0d0
           Ukp2 = 0.0d0
           Ukp3 = 0.0d0
          if(ndrvA.eq.2)  UCA = - ASSA*SM1A
c         if(ndrvA.eq.2)  UCA = - AASA*SM1A
          if(ndrvB.eq.2)  UCB = - ASSB*SM1B
c         if(ndrvB.eq.2)  UCB = - BBSB*SM1B
           VCPAR = (D1*UCA+D2*UCB)*0.50d0
           Ukp2=Ukp+VCPAR
           Ukp3=Ukp+ACPAR*VCPAR

           ! fenglai debug
c           write(6,*)"emil's new us", Ukp3
c           write(6,*)"emil's us", Ukp2
c           write(6,*)"emil's par us", VCPAR
c           write(6,*)"emil's opp us", Ukp
c        else  !  NMO
c          VCPAR= zero
c        endif
ccccccccccccccc
         Fpar = ACPAR*VCPAR
          F(i)= Fpar
         OD1 = 2.0d0*D1*ffcor*EFN2B + D1*AASA*SM2A
         OD2 = 2.0d0*D2*ffcor*EFN2A + D2*BBSB*SM2B
          IF(OD1.le.(-dsmall)) OD1=0.0d0
          IF(OD2.le.(-dsmall)) OD2=0.0d0
c          ODEL(i) = OD1+OD2
c           ODEL(i) = OD1+OD2 + ds1
c            ODELW(i) = ODEL(i)*W(i)
c        IF(ICEN.ne.NCENTER ) then
c            NCENTER = ICEN
c         endif
c          ODDSUM(NCENTER) = ODDSUM(NCENTER)+ODELW(i)
          
c       write(*,*) 'Fop = ', Fop
c       write(*,*) 'Fpar = ', Fpar
c
c Begin building the SCF potential:
        IF (IMethod.eq.2) THEN
c
            IF((ndrvA.eq.0).or.(ndrvB.eq.0)) GO TO 888
         If((D1.gt.detol.and.ndrvA.eq.2).and.(D2.gt.detol
     $      .and.ndrvB.eq.2)) then
       IF(EFNB.gt.dsmall) then
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
        IF(EFNA.gt.dsmall) then
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
c
        if ((abs(fc1).gt.dsmall).or.(abs(fc2).gt.dsmall)) then
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

          dfcdVxaa = (dfbdVxa+(dfadVxa-dfbdVxa)*qqx
     $    -fffa*dfbdVxa*qqfa-fffb*dfadVxa*qqfa)
c    $    -fffa*dfbdVxa*qqfa-fffb*dfadVxa*qqfa)*0.50d0

          dfcdrbb = dfbdrb+(dfadrb-dfbdrb)*qqx
     $    -fffa*dfbdrb*qqfa-fffb*dfadrb*qqfa

          dfcdgbb = dfbdgb+(dfadgb-dfbdgb)*qqx
     $    -fffa*dfbdgb*qqfa-fffb*dfadgb*qqfa

          dfcdtbb = dfbdtb+(dfadtb-dfbdtb)*qqx
     $    -fffa*dfbdtb*qqfa-fffb*dfadtb*qqfa

          dfcdlbb = dfbdlb+(dfadlb-dfbdlb)*qqx
     $    -fffa*dfbdlb*qqfa-fffb*dfadlb*qqfa

          dfcdVxbb = (dfbdVxb+(dfadVxb-dfbdVxb)*qqx
     $    -fffa*dfbdVxb*qqfa-fffb*dfadVxb*qqfa)
c    $    -fffa*dfbdVxb*qqfa-fffb*dfadVxb*qqfa)*0.50d0

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
          dfopdra = 0.50d0*(dfcdraa*ttgg+fcor*(UX2-ttdA))
          dfopdga = 0.50d0*dfcdgaa*ttgg
          dfopdta = 0.50d0*dfcdtaa*ttgg
          dfopdla = 0.50d0*dfcdlaa*ttgg
          dfopdVxa = 0.50d0*dfcdVxaa*ttgg
          dfopdVxa = dfopdVxa+ 0.50d0*fcor*D2/D1
          dfopdrb = 0.50d0*(dfcdrbb*ttgg+fcor*(UX1-ttdB))
          dfopdgb = 0.50d0*dfcdgbb*ttgg
          dfopdtb = 0.50d0*dfcdtbb*ttgg
          dfopdlb = 0.50d0*dfcdlbb*ttgg
c         dfopdVxb = 0.50d0*dfcdVxbb*ttgg
          dfopdVxb = 0.50d0*(dfcdVxbb*ttgg+fcor*D1/D2)
      endif  ! NA vs NB
c=====up to here the opposite spin SCF part is done
c       write(6,*)"emil's f factor on RA",  dfcdraa
c       write(6,*)"emil's f factor on RB",  dfcdrbb
c       write(6,*)"emil's f factor on GAA", dfcdgaa
c       write(6,*)"emil's f factor on GBB", dfcdgbb
c       write(6,*)"emil's f factor on TA",  dfcdtaa
c       write(6,*)"emil's f factor on TB",  dfcdtbb
c       write(6,*)"emil's f factor on LA",  dfcdlaa
c       write(6,*)"emil's f factor on LB",  dfcdlbb
c       write(6,*)"emil's f factor on EXA", dfcdVxaa
c       write(6,*)"emil's f factor on EXB", dfcdVxbb

c
 888      CONTINUE
            facta = zero
            factb=zero
c       IF((NMO.eq.2).and.(NA.ne.NB)) then
c
         if(abs(SM2A).gt.dsmall) facta= ASA/SM2A
       If((D1.gt.detol).and.(abs(SM2A).gt.detol).and.(ndrvA.eq.2)) then
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
       if((D1.gt.detol).and.(abs(SM2B).gt.detol).and.(D2.gt.detol)
     $ .and.(ndrvB.eq.2)) then
c      if(D1.gt.detol.and.(abs(SM2B).gt.detol)) then
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
          if(abs(SM2B).gt.dsmall) factb= ASB/SM2B
       IF((D2.gt.detol).and.(abs(SM2B).gt.detol).and.(ndrvB.eq.2)) then
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
      IF((D2.gt.detol).and.(abs(SM2A).gt.detol).and.(D1.gt.detol)
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
c
      IF ((D1.gt.detol).and.(ndrvA.eq.2)) then
            dfpardra = -ACPAR*0.50d0*(ASSA*SM1A+D1*dASSAdra*SM1A
     $               + D2*dASSBdra*SM1B + D1*ASSA*dM1dra)
            dfpar2dra = -0.50d0*(ASSA*SM1A+D1*dASSAdra*SM1A
     $               + D2*dASSBdra*SM1B + D1*ASSA*dM1dra)
            dfpardga = -ACPAR*0.50d0*(D1*dASSAdga*SM1A
     $               + D2*dASSBdga*SM1B + D1*ASSA*dM1dga)
            dfpar2dga = -0.50d0*(D1*dASSAdga*SM1A
     $               + D2*dASSBdga*SM1B + D1*ASSA*dM1dga)
            dfpardta = -ACPAR*0.50d0*(D1*dASSAdta*SM1A
     $               + D2*dASSBdta*SM1B + D1*ASSA*dM1dta)
            dfpar2dta = -0.50d0*(D1*dASSAdta*SM1A
     $               + D2*dASSBdta*SM1B + D1*ASSA*dM1dta)
            dfpardla = -ACPAR*0.50d0*(D1*dASSAdla*SM1A
     $               + D2*dASSBdla*SM1B + D1*ASSA*dM1dla)
            dfpar2dla = -0.50d0*(D1*dASSAdla*SM1A
     $               + D2*dASSBdla*SM1B + D1*ASSA*dM1dla)
            dfpardVxa = -ACPAR*0.50d0*(D1*dASSAdVxa*SM1A
     $                + D2*dASSBdVxa*SM1B + D1*ASSA*dM1dVxa)
            dfpar2dVxa = -0.50d0*(D1*dASSAdVxa*SM1A
     $                + D2*dASSBdVxa*SM1B + D1*ASSA*dM1dVxa)
         else
            dfpardra = zero
            dfpar2dra = zero
            dfpardga = zero
            dfpar2dga = zero
            dfpardta = zero
            dfpar2dta = zero
            dfpardla = zero
            dfpar2dla = zero
           dfpardVxa = zero
           dfpar2dVxa = zero
         endif
      IF ((D2.gt.detol)) then
          dfpardrb = -ACPAR*0.50d0*(ASSB*SM1B+D2*dASSBdrb*SM1B
     $             + D1*dASSAdrb*SM1A + D2*ASSB*dM1drb)
          dfpar2drb = -0.50d0*(ASSB*SM1B+D2*dASSBdrb*SM1B
     $             + D1*dASSAdrb*SM1A + D2*ASSB*dM1drb)
          dfpardgb = -ACPAR*0.50d0*(D2*dASSBdgb*SM1B
     $             + D1*dASSAdgb*SM1A + D2*ASSB*dM1dgb)
          dfpar2dgb = -0.50d0*(D2*dASSBdgb*SM1B
     $             + D1*dASSAdgb*SM1A + D2*ASSB*dM1dgb)
          dfpardtb = -ACPAR*0.50d0*(D2*dASSBdtb*SM1B
     $             + D1*dASSAdtb*SM1A + D2*ASSB*dM1dtb)
          dfpar2dtb = -0.50d0*(D2*dASSBdtb*SM1B
     $             + D1*dASSAdtb*SM1A + D2*ASSB*dM1dtb)
          dfpardlb = -ACPAR*0.50d0*(D2*dASSBdlb*SM1B
     $             + D1*dASSAdlb*SM1A + D2*ASSB*dM1dlb)
          dfpar2dlb = -0.50d0*(D2*dASSBdlb*SM1B
     $             + D1*dASSAdlb*SM1A + D2*ASSB*dM1dlb)
          dfpardVxb = -ACPAR*0.50d0*(D2*dASSBdVxb*SM1B
     $               + D1*dASSAdVxb*SM1A + D2*ASSB*dM1dVxb)
          dfpar2dVxb = -0.50d0*(D2*dASSBdVxb*SM1B
     $               + D1*dASSAdVxb*SM1A + D2*ASSB*dM1dVxb)
       else
          dfpardrb = zero
          dfpar2drb = zero
          dfpardgb = zero
          dfpar2dgb = zero
          dfpardtb = zero
          dfpar2dtb = zero
          dfpardlb = zero
          dfpar2dlb = zero
          dfpardVxb =zero
          dfpar2dVxb =zero
       endif
c
         dfdra = dfpar2dra+dfopdra
         dfdga = dfpar2dga+dfopdga
         dfdta = dfpar2dta+dfopdta
         dfdla = dfpar2dla+dfopdla
         dfdVxa = dfpar2dVxa+dfopdVxa
         dfdrb = dfpar2drb+dfopdrb
         dfdgb = dfpar2dgb+dfopdgb
         dfdtb = dfpar2dtb+dfopdtb
         dfdlb = dfpar2dlb+dfopdlb
         dfdVxb = dfpar2dVxb+dfopdVxb
           dzxdgg=1.0d0
c

         dEopdn1=zero
         dEopdn2=zero
         dEopdg1=zero
         dEopdg2=zero
         dEopdt1=zero
         dEopdt2=zero
         dEopdl1=zero
         dEopdl2=zero
         dEopdVx1=zero
         dEopdVx2=zero
         dEpardn1=zero
         dEpardn2=zero
         dEpardg1=zero
         dEpardg2=zero
         dEpardt1=zero
         dEpardt2=zero
         dEpardl1=zero
         dEpardl2=zero
         dEpardVx1=zero
         dEpardVx2=zero
c
         duopdn1=zero
         duopdn2=zero
         duopdg1=zero
         duopdg2=zero
         duopdt1=zero
         duopdt2=zero
         duopdl1=zero
         duopdl2=zero
         duopdVx1=zero
         duopdVx2=zero
         dupardn1=zero
         dupardn2=zero
         dupardg1=zero
         dupardg2=zero
         dupardt1=zero
         dupardt2=zero
         dupardl1=zero
         dupardl2=zero
         dupardVx1=zero
         dupardVx2=zero
         EOPP = zero
         EP1=zero
         EP2=zero
         EPAR=zero
         UOPP = zero
         UPAR=zero
         UP1=zero
         UP2=zero
         ULAOPP = zero
         ULAP1=zero
         ULAP2 = zero
cccccccc kp14 dyn C begins
c       IF((D1.gt.detol).or.(D2.gt.detol)) then
         CALL B14dync(D1,VX1,EFNA,ASSA,dNdra,
     $   dNdga,dNdta,dNdla,dNdVxa,DDA,dg3dra,dg3dta,dg3dga,
     $   dfcdraa,dfcdgaa,dfcdtaa,dfcdlaa,dfcdVxaa,
     $   dASSAdra,dASSAdga,dASSAdta,dASSAdla,dASSAdVxa,
     $   dASSAdrb,dASSAdgb,dASSAdtb,dASSAdlb,dASSAdVxb,
     $   D2,VX2,EFNB,ASSB,
     $   dNdrb,dNdgb,dNdtb,dNdlb,dNdVxb,DDB,dg3drb,dg3dtb,dg3dgb,
     $   dfcdrbb,dfcdgbb, dfcdtbb,dfcdlbb,dfcdVxbb,
     $   dASSBdra,dASSBdga,dASSBdta,dASSBdla,dASSBdVxa,
     $   dASSBdrb,dASSBdgb,dASSBdtb,dASSBdlb,dASSBdVxb,fcor,
c
     $   dEopdn1,dEopdn2,dEopdg1,dEopdg2,dEopdt1,dEopdt2,
     $   dEopdl1,dEopdl2,dEopdVx1,dEopdVx2,dEpardn1,dEpardn2,
     $   dEpardg1,dEpardg2,dEpardt1,dEpardt2,dEpardl1,dEpardl2,
     $   dEpardVx1,dEpardVx2,
c
     $   duopdn1,duopdn2,duopdg1,duopdg2,duopdt1,duopdt2,
     $   duopdl1,duopdl2,duopdVx1,duopdVx2,dupardn1,dupardn2,
     $   dupardg1,dupardg2,dupardt1,dupardt2,dupardl1,dupardl2,
     $   dupardVx1,dupardVx2,UOPP,UPAR,UP1,UP2,ULAOPP,ULAP1,ULAP2,
     $   ULAPAR,EOPP,EP1,EP2,EPAR,detol,dsmall,ndrv,qlamb,NA,NB,icorr)
c       endif

cccccccccccccccccccccccccccccccc
c Initializing first
              PBkp = 0.59d0
               Ulamb = qlamb
                dPdz = 0.0d0
              YKPW = Tolbig2
              EFNW = 0.0d0
              Ekp14 = 0.0d0
              End0 = 0.0d0
              Elamb14=0.0d0
c             BXX = 0.0d0
          EDYN = DCOP*UOPP+DCPAR*UPAR
c   EAVDYN is the lambda integrated dynamic Ec
c   When EAVDYN is used then DCOP and DCPAR should be close to 1.0
          EAVDYN = DCOP*EOPP+DCPAR*EPAR
          UDYN = UOPP+UPAR

           ! fenglai debug
c           write(6,*)"emil's ud", UDYN
c   ULADYN is lambda dependent dynamic Uc
c ???      ULADYN = DCOP*ULAOPP+DCPAR*ULAPAR
           ULADYN = ULAOPP+ULAPAR
           dzxdgg=1.0d0
          IF(UDYN.ge.(dsmall)) then 
            icountUD = icountUD + 1
c          write(*,*) 'UDYN = ', UDYN,'i = ',i
          endif
        IF(Ukp2.ge.(dsmall)) then
              icountUnd = icountUnd + 1
c          write(*,*) 'Ukp2 = ', Ukp2,'i = ',i
          endif
        IF(UDYN.lt.(-dsmall)) then
         IF((Ukp3.lt.(-dsmall))) then
                  ggb = Ukp2/UDYN 
                  ggg = ggb+delt3
           IF((ggb.ge.smth3).and.(ggb.lt.smth4)) then
                Zkp = ggg*ggg/(4.0D0*delt2)
                dzxdgg = ggg/(2.0D0*delt2)
           elseif (ggb.lt.smth3) then
                Zkp = dsmall1
                dzxdgg = 0.0d0
           else
               Zkp = Ukp2/UDYN 
               dzxdgg = 1.0d0
           endif
c           write(6,*)"emil's z, originally defined", zkp
c                  YKPW = Zkp*W(i)
c                  EFNW = EFNA*W(i)
c                 BXX = Ukp/(Ukp+UDYN)
c          YSUM(NCENTER) = YSUM(NCENTER)+Zkp/NGrid
c          AVNEF(NCENTER) =  AVNEF(NCENTER) +EFNA/NGrid
ccccccccccccccccccccccccccccc
c   boost of the B05 f factor starts with Def BXX:
c        IF((BXX.gt.0.980d0).and.(ffcor.lt.0.990d0))   then
c                ffcor = 1.0d0
c                ffcor = fcor
c        endif
c        Ukp = 0.50d0*fcor*ttgg
c        UUkp = 0.50d0*ttgg
c         IF((UDYN.lt.(-detol)).and.(Ukp.lt.(-detol))) then
c                 Zkp = Ukp/UDYN
c           endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              Ulamb= 1.0d0*qlamb
              tt0 = Bkp*Zkp
              tlam0  = tt0*qlamb
              tte = exp(-tt0)
              tlame = exp(-tlam0)
              tt1 = tte + 1.0d0
              tt2 = tte -1.0d0
              tlam1 = tlame + 1.0d0
              tlam2 = tlame - 1.0d0
              tlog1 = 0.0d0
              dPdz = 0.0d0
              aakp = -1.0d0
              tt3 = 0.0d0
           IF(tt0.gt.detol) then
                  tt3 = log(tt1/2.0d0)
                  tlog1 = -1.0d0-2.0d0*tt3/tt0
               IF(abs(tt2).gt.dsmall) then
                   aakp = tt1/tt2
                else
                   aakp = -1.0d0
               endif
                  PBkp = aakp*tlog1
                  Ulamb = aakp*tlam2/tlam1
                  dPdz=0.0d0
                  t1 = Bkp*Zkp
                  t3 = exp(-t1)
                  t4 = Zkp*Zkp
                  t6 = Bkp*Bkp
                  t9 = t3*Zkp
                  t13 = log(t3 + 0.1D1)
                  t18 = exp(-0.2D1*t1)
                  t28 = (t3 - 0.1D1)**2
                  t29 = t28*t4*Bkp
             IF(abs(t29).gt.dsmall) then
       dPdz = (-0.2D1*t3*t4*t6+0.7725887222397820D0*t9*Bkp -
     $ 0.4D1*t9*Bkp*t13 + 0.2D1*t18*t1 - 0.1386294361119891D1*t18
     $ + 0.13862943611198910D1 + 0.2D1*t13*t18 - 0.2D1*t13)/t29
                     endif
                IF(EFNA.le.dsmall) then
                    PBkp = 0.590d0   
                    Ulamb = 1.0d0*qlamb 
                dPdz = 0.0d0
                endif
                IF(EFNB.le.dsmall) then
                    PBkp = 0.590d0   
                    Ulamb = 1.0d0*qlamb 
                dPdz = 0.0d0
                endif
ccccccccccccccc
          elseif(tt0.gt.dsmall) then   !  tt0
                tlog1 = -tt0/4.0d0+tt0**3/96.0d0
                aakp = -2.0d0/tt0-tt0/6.0d0+tt0**3/360.0d0
                  bp1 = 0.00002893518519d0
                  bp2 = 0.002430555556d0
                  bp3 = 0.02083333330d0
                  bp4 = 0.0001736111111d0
                  bp5 = 0.009722222224d0
                  bp6 = 0.041666666660d0
                 PBkp = bp1*tt0**6-bp2*tt0**4+bp3*tt0*tt0+0.50d0  
                 dPdz = bp4*(tt0**5)- bp5*tt0**4+bp6*tt0
             else
                    PBkp = 0.59d0   
                dPdz = 0.0d0
               endif  ! tt0
c
              IF(EFNA.le.dsmall) then
                   PBkp = 0.59d0
                   Ulamb = qlamb
               dPdz = 0.0d0
              endif
              IF(EFNB.le.dsmall) then
                   PBkp = 0.59d0
                   Ulamb = qlamb
               dPdz = 0.0d0
              endif
             else    !   Ukp3
              PBkp = 0.59d0
               Ulamb = qlamb 
                dPdz = 0.0d0
            endif   ! Ukp3
          else ! udyn
              PBkp = 1.0d0
                Ulamb = qlamb
                  dPdz = 0.0d0
          endif  ! udyn
             DNA13 = D1**third
             DA43 = DNA13*D1
             DNB13 = D2**third
             DB43 = DNB13*D2
             Ekp14 = Ukp2*PBkp
           If((icorr.eq.1).or.(icorr.eq.6)) then
             Ekp14 = Ukp3*PBkp
c           elseif(icorr.eq.6) then
c              Ekp14 = Ukp3*PBkp
           endif
              Elamb14 = Ukp*Ulamb
c
        ! fenglai debug
        IF((icorr.eq.4).or.(icorr.eq.2).or.(icorr.eq.5)) then
              alfpi = alf/PI
             IF(Zkp.gt.dtol) then
              ZZkp = 1.0d0/Zkp
             else
c              ZZkp = Tolbig2
c              ZZkp = 50.0d0 
               ZZkp = 1.0d0/detol
             endif
         DFF = (DFA)**third+(DFB)**third
c      DFF = dsign(dabs(DFA)**third,DFA)+dsign(dabs(DFB)**third,DFB)
ccccccccccc

              ! fenglai debug
c              write(6,*)"emil's DA", DFA**third
c              write(6,*)"emil's DB", DFB**third
c              write(6,*)"emil's QAC", PBkp
c              write(6,*)"emil's DFF", DFF
ccccccccccc
              es1 = exp(-alf*ZZkp*ZZkp)
              S1 = sqrt(alfpi)*es1
              eshift = 0.50d0*DFF*S1
c              write(6,*)"eaz2", es1
c              write(6,*)"S1", S1
             IF(Zkp.gt.dtol) then
              dS1dz = 2.0d0*sqrt(alfpi)*alf*es1/(Zkp*Zkp)/Zkp
              else
               dS1dz = 0.0d0
              endif
c
              End0 = Ukp3*PBkp
              Ekp14 = End0*(1.0d0+eshift)
              Elamb14 = Ukp3*Ulamb*(1.0d0 + eshift)

              ! fenglai
c              write(6,*)"emil, eshift", eshift
c              write(6,*)"emil, End0", end0
c              write(6,*)"emil's QAC", PBkp
c              write(6,*)"emil's EKP14", Ekp14
            endif   !   icorr 4
c
c          IF(icorr.eq.1) then
c             F(i) = F(i) + Ekp14 + EDYN
c
           if(icorr.eq.2) then
c            F(i) = VCPAR*qlamb+Elamb14+ULADYN
             F(i) = Elamb14+ULADYN
c
           elseif((icorr.eq.3).or.(icorr.eq.5).or.(icorr.eq.1)) then
            F(i) = kp14fac*Ekp14 + EAVDYN
c
           elseif((icorr.eq.4).or.(icorr.eq.6)) then
              F(i) = Ekp14 + EDYN
          endif
            Ebatch = Ebatch + F(i)*W(i)
c
c  addition to the SCF potential from the prefactor follows below:
           dZdra = 0.0d0
           dZdga = 0.0d0
           dZdta = 0.0d0
           dZdla = 0.0d0
           dZdVxa = 0.0d0
           dZdrb = 0.0d0
           dZdgb = 0.0d0
           dZdtb = 0.0d0
           dZdlb = 0.0d0
           dZdVxb = 0.0d0
c
           dS1dra = 0.0d0
           dS1dga = 0.0d0
           dS1dta = 0.0d0
           dS1dla = 0.0d0
           dS1dVxa = 0.0d0
           dS1drb = 0.0d0
           dS1dgb = 0.0d0
           dS1dtb = 0.0d0
           dS1dlb = 0.0d0
           dS1dVxb = 0.0d0
c
       IF(UDYN.lt.(-detol)) then
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdra,duopdn1,dupardn1,dZdra,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdga,duopdg1,dupardg1,dZdga,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdta,duopdt1,dupardt1,dZdta,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdla,duopdl1,dupardl1,dZdla,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdVxa,duopdVx1,
     $   dupardVx1,dZdVxa,detol)
c
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdrb,duopdn2,dupardn2,dZdrb,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdgb,duopdg2,dupardg2,dZdgb,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdtb,duopdt2,dupardt2,dZdtb,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdlb,duopdl2,dupardl2,dZdlb,detol)
      CALL ZDERIV(Zkp,dzxdgg,UDYN,dfdVxb,duopdVx2,dupardVx2,
     $            dZdVxb,detol)
       endif

       ! fenglai debug
c       write(6,*)"emil's us on RA",  dfdra
c       write(6,*)"emil's us on RB",  dfdrb
c       write(6,*)"emil's us on GAA", dfdga
c       write(6,*)"emil's us on GBB", dfdgb
c       write(6,*)"emil's us on TA",  dfdta
c       write(6,*)"emil's us on TB",  dfdtb
c       write(6,*)"emil's us on LA",  dfdla
c       write(6,*)"emil's us on LB",  dfdlb
c       write(6,*)"emil's us on EXA", dfdVxa
c       write(6,*)"emil's us on EXB", dfdVxb
c
c       write(6,*)"emil's ud par on RA",  dupardn1
c       write(6,*)"emil's ud par on RB",  dupardn2
c       write(6,*)"emil's ud par on GAA", dupardg1
c       write(6,*)"emil's ud par on GBB", dupardg2
c       write(6,*)"emil's ud par on TA",  dupardt1
c       write(6,*)"emil's ud par on TB",  dupardt2
c       write(6,*)"emil's ud par on LA",  dupardl1
c       write(6,*)"emil's ud par on LB",  dupardl2
c       write(6,*)"emil's ud par on EXA", dupardVx1
c       write(6,*)"emil's ud par on EXB", dupardVx2
c
c       write(6,*)"emil's ud opp on RA",  duopdn1
c       write(6,*)"emil's ud opp on RB",  duopdn2
c       write(6,*)"emil's ud opp on GAA", duopdg1
c       write(6,*)"emil's ud opp on GBB", duopdg2
c       write(6,*)"emil's ud opp on TA",  duopdt1
c       write(6,*)"emil's ud opp on TB",  duopdt2
c       write(6,*)"emil's ud opp on LA",  duopdl1
c       write(6,*)"emil's ud opp on LB",  duopdl2
c       write(6,*)"emil's ud opp on EXA", duopdVx1
c       write(6,*)"emil's ud opp on EXB", duopdVx2

c
         dPkpdra = dPdz*dZdra 
         dPkpdga = dPdz*dZdga 
         dPkpdta = dPdz*dZdta 
         dPkpdla = dPdz*dZdla 
         dPkpdVxa= dPdz*dZdVxa 
c
         dPkpdrb = dPdz*dZdrb
         dPkpdgb = dPdz*dZdgb
         dPkpdtb = dPdz*dZdtb
         dPkpdlb = dPdz*dZdlb
         dPkpdVxb= dPdz*dZdVxb
c
      IF(icorr.eq.1) then
         dEkpdra = dfopdra*PBkp + Ukp*dPkpdra
         dEkpdga = dfopdga*PBkp + Ukp*dPkpdga
         dEkpdla = dfopdla*PBkp + Ukp*dPkpdla
         dEkpdta = dfopdta*PBkp + Ukp*dPkpdta
         dEkpdVxa = dfopdVxa*PBkp + Ukp*dPkpdVxa
c
         dEkpdrb = dfopdrb*PBkp + Ukp*dPkpdrb
         dEkpdgb = dfopdgb*PBkp + Ukp*dPkpdgb
         dEkpdlb = dfopdlb*PBkp + Ukp*dPkpdlb
         dEkpdtb = dfopdtb*PBkp + Ukp*dPkpdtb
         dEkpdVxb = dfopdVxb*PBkp + Ukp*dPkpdVxb
c 
       D1F(i,POS_RA) = dEkpdra+dfpardra+DCOP*duopdn1+DCPAR*dupardn1
       D1F(i,POS_GAA)= dEkpdga+dfpardga+DCOP*duopdg1+DCPAR*dupardg1
       D1F(i,POS_TA) = dEkpdta+dfpardta+DCOP*duopdt1+DCPAR*dupardt1
       D1F(i,POS_LA) = dEkpdla+dfpardla+DCOP*duopdl1+DCPAR*dupardl1
      D1FEXCH(i,1)  = dEkpdVxa+dfpardVxa+DCOP*duopdVx1+DCPAR*dupardVx1
c
       D1F(i,POS_RB) = dEkpdrb+dfpardrb+DCOP*duopdn2+DCPAR*dupardn2
       D1F(i,POS_GBB)= dEkpdgb+dfpardgb+DCOP*duopdg2+DCPAR*dupardg2
       D1F(i,POS_TB) = (dEkpdtb+dfpardtb+DCOP*duopdt2+DCPAR*dupardt2)
       D1F(i,POS_LB) = dEkpdlb+dfpardlb+DCOP*duopdl2+DCPAR*dupardl2
       D1FEXCH(i,2)  = dEkpdVxb+dfpardVxb+DCOP*duopdVx2+DCPAR*dupardVx2
c
      elseif (icorr.eq.6) then
         dEkpdra = (dfopdra+dfpardra)*PBkp + Ukp3*dPkpdra
         dEkpdga = (dfopdga+dfpardga)*PBkp + Ukp3*dPkpdga
         dEkpdla = (dfopdla+dfpardla)*PBkp + Ukp3*dPkpdla
         dEkpdta = (dfopdta+dfpardta)*PBkp + Ukp3*dPkpdta
         dEkpdVxa = (dfopdVxa+dfpardVxa)*PBkp + Ukp3*dPkpdVxa
c
         dEkpdrb = (dfopdrb+dfpardrb)*PBkp + Ukp3*dPkpdrb
         dEkpdgb = (dfopdgb+dfpardgb)*PBkp + Ukp3*dPkpdgb
         dEkpdlb = (dfopdlb+dfpardlb)*PBkp + Ukp3*dPkpdlb
         dEkpdtb = (dfopdtb+dfpardtb)*PBkp + Ukp3*dPkpdtb
         dEkpdVxb = (dfopdVxb+dfpardVxb)*PBkp + Ukp3*dPkpdVxb
c
       D1F(i,POS_RA) = dEkpdra +DCOP*duopdn1+DCPAR*dupardn1
       D1F(i,POS_GAA)= dEkpdga +DCOP*duopdg1+DCPAR*dupardg1
       D1F(i,POS_TA) = dEkpdta +DCOP*duopdt1+DCPAR*dupardt1
       D1F(i,POS_LA) = dEkpdla +DCOP*duopdl1+DCPAR*dupardl1
      D1FEXCH(i,1)  = dEkpdVxa +DCOP*duopdVx1+DCPAR*dupardVx1
c
       D1F(i,POS_RB) = dEkpdrb +DCOP*duopdn2+DCPAR*dupardn2
       D1F(i,POS_GBB)= dEkpdgb +DCOP*duopdg2+DCPAR*dupardg2
       D1F(i,POS_TB) = dEkpdtb +DCOP*duopdt2+DCPAR*dupardt2
       D1F(i,POS_LB) = dEkpdlb +DCOP*duopdl2+DCPAR*dupardl2
      D1FEXCH(i,2)  = dEkpdVxb +DCOP*duopdVx2+DCPAR*dupardVx2
c
      elseif (icorr.eq.5) then
c   icorr=5 must set ACPAR=1, DCOP= 1, DCPAR = 1
         dS1dra = dS1dz*dZdra
         dS1dga = dS1dz*dZdga
         dS1dta = dS1dz*dZdta
         dS1dla = dS1dz*dZdla
         dS1dVxa= dS1dz*dZdVxa
c
         dS1drb = dS1dz*dZdrb
         dS1dgb = dS1dz*dZdgb
         dS1dtb = dS1dz*dZdtb
         dS1dlb = dS1dz*dZdlb
         dS1dVxb= dS1dz*dZdVxb
c
         dSnddra = 0.50d0*(ddffdra*S1+DFF*dS1dra)
         dSnddga = 0.50d0*(ddffdga*S1+DFF*dS1dga)
         dSnddta = 0.50d0*(ddffdta*S1+DFF*dS1dta)
         dSnddla = 0.50d0*DFF*dS1dla
         dSnddVxa =0.50d0*DFF*dS1dVxa
 
         dSnddrb = 0.50d0*(ddffdrb*S1+DFF*dS1drb)
         dSnddgb = 0.50d0*(ddffdgb*S1+DFF*dS1dgb)
         dSnddtb = 0.50d0*(ddffdtb*S1+DFF*dS1dtb)
         dSnddlb = 0.50d0*DFF*dS1dlb
         dSnddVxb =0.50d0*DFF*dS1dVxb
c

c         write(6,*)"nus deriv on RA: ", dfopdra+dfpardra
c         write(6,*)"nus deriv on RB: ", dfopdrb+dfpardrb
c         write(6,*)"nus deriv on GAA: ",dfopdga+dfpardga
c         write(6,*)"nus deriv on GBB: ",dfopdgb+dfpardgb
c         write(6,*)"nus deriv on TA: ", dfopdta+dfpardta
c         write(6,*)"nus deriv on TB: ", dfopdtb+dfpardtb
c         write(6,*)"nus deriv on LA: ", dfopdla+dfpardla
c         write(6,*)"nus deriv on LB: ", dfopdla+dfpardla
c         write(6,*)"nus deriv on EXA: ",dfopdVxa+dfpardVxa
c         write(6,*)"nus deriv on EXB: ",dfopdVxb+dfpardVxb
         dEkpdra = (dfopdra+dfpardra)*PBkp + Ukp3*dPkpdra
         dEkpdga = (dfopdga+dfpardga)*PBkp + Ukp3*dPkpdga
         dEkpdla = (dfopdla+dfpardla)*PBkp + Ukp3*dPkpdla
         dEkpdta = (dfopdta+dfpardta)*PBkp + Ukp3*dPkpdta
         dEkpdVxa =(dfopdVxa+dfpardVxa)*PBkp+Ukp3*dPkpdVxa
c
         dEkpdrb = (dfopdrb+dfpardrb)*PBkp + Ukp3*dPkpdrb
         dEkpdgb = (dfopdgb+dfpardgb)*PBkp + Ukp3*dPkpdgb
         dEkpdlb = (dfopdlb+dfpardlb)*PBkp + Ukp3*dPkpdlb
         dEkpdtb = (dfopdtb+dfpardtb)*PBkp + Ukp3*dPkpdtb
         dEkpdVxb = (dfopdVxb+dfpardVxb)*PBkp + Ukp3*dPkpdVxb
c
         dEnddra = dEkpdra*(1.0d0+eshift)+End0*dSnddra
         dEnddga = dEkpdga*(1.0d0+eshift)+End0*dSnddga
         dEnddta = dEkpdta*(1.0d0+eshift)+End0*dSnddta
         dEnddla = dEkpdla*(1.0d0+eshift)+End0*dSnddla
         dEnddVxa = dEkpdVxa*(1.0d0+eshift)+End0*dSnddVxa
c
         dEnddrb = dEkpdrb*(1.0d0+eshift)+End0*dSnddrb
         dEnddgb = dEkpdgb*(1.0d0+eshift)+End0*dSnddgb
         dEnddtb = dEkpdtb*(1.0d0+eshift)+End0*dSnddtb
         dEnddlb = dEkpdlb*(1.0d0+eshift)+End0*dSnddlb
         dEnddVxb = dEkpdVxb*(1.0d0+eshift)+End0*dSnddVxb

         ! fenglai debug
c       write(6,*)"emil's dQdz",  dPdz
c       write(6,*)"emil's z on RA",  dZdra
c       write(6,*)"emil's z on RB",  dZdrb
c       write(6,*)"emil's z on GAA", dZdga
c       write(6,*)"emil's z on GBB", dZdgb
c       write(6,*)"emil's z on TA",  dZdta
c       write(6,*)"emil's z on TB",  dZdtb
c       write(6,*)"emil's z on LA",  dZdla
c       write(6,*)"emil's z on LB",  dZdlb
c       write(6,*)"emil's z on EXA", dZdVxa
c       write(6,*)"emil's z on EXB", dZdVxb
c
c       write(6,*)"emil's qac on RA", dPkpdra
c       write(6,*)"emil's qac on RB", dPkpdrb
c       write(6,*)"emil's qac on GAA", dPkpdga
c       write(6,*)"emil's qac on GBB", dPkpdgb
c       write(6,*)"emil's qac on TA",  dPkpdta
c       write(6,*)"emil's qac on TB",  dPkpdtb
c       write(6,*)"emil's qac on LA",  dPkpdla
c       write(6,*)"emil's qac on LB",  dPkpdlb
c       write(6,*)"emil's qac on EXA", dPkpdVxa
c       write(6,*)"emil's qac on EXB", dPkpdVxb

       D1F(i,POS_RA) = kp14fac*dEnddra+DCOP*dEopdn1+DCPAR*dEpardn1
       D1F(i,POS_GAA)= kp14fac*dEnddga +DCOP*dEopdg1+DCPAR*dEpardg1
       D1F(i,POS_TA) = kp14fac*dEnddta +DCOP*dEopdt1+DCPAR*dEpardt1
       D1F(i,POS_LA) = kp14fac*dEnddla +DCOP*dEopdl1+DCPAR*dEpardl1
       D1FEXCH(i,1)  = kp14fac*dEnddVxa +DCOP*dEopdVx1+DCPAR*dEpardVx1
c
       D1F(i,POS_RB) = kp14fac*dEnddrb +DCOP*dEopdn2+DCPAR*dEpardn2
       D1F(i,POS_GBB)= kp14fac*dEnddgb +DCOP*dEopdg2+DCPAR*dEpardg2
       D1F(i,POS_TB) = kp14fac*dEnddtb +DCOP*dEopdt2+DCPAR*dEpardt2
       D1F(i,POS_LB) = kp14fac*dEnddlb +DCOP*dEopdl2+DCPAR*dEpardl2
       D1FEXCH(i,2)  = kp14fac*dEnddVxb +DCOP*dEopdVx2+DCPAR*dEpardVx2

       ! fenglai debug
c       write(6,*)"emil's functional derivatives on RA",
c     $ kp14fac*dEnddra+DCOP*dEopdn1+DCPAR*dEpardn1
c       write(6,*)"emil's functional derivatives on RB",
c     $ kp14fac*dEnddrb +DCOP*dEopdn2+DCPAR*dEpardn2
c       write(6,*)"emil's functional derivatives on GAA",
c     $ kp14fac*dEnddga +DCOP*dEopdg1+DCPAR*dEpardg1
c       write(6,*)"emil's functional derivatives on GBB",
c     $ kp14fac*dEnddgb +DCOP*dEopdg2+DCPAR*dEpardg2
c       write(6,*)"emil's functional derivatives on TA",
c     $ kp14fac*dEnddta +DCOP*dEopdt1+DCPAR*dEpardt1
c       write(6,*)"emil's functional derivatives on TB",
c     $ kp14fac*dEnddtb +DCOP*dEopdt2+DCPAR*dEpardt2
c       write(6,*)"emil's functional derivatives on LA",
c     $ kp14fac*dEnddla +DCOP*dEopdl1+DCPAR*dEpardl1
c       write(6,*)"emil's functional derivatives on LB",
c     $ kp14fac*dEnddlb +DCOP*dEopdl2+DCPAR*dEpardl2
c       write(6,*)"emil's functional derivatives on EXA",
c     $ kp14fac*dEnddVxa +DCOP*dEopdVx1+DCPAR*dEpardVx1
c       write(6,*)"emil's functional derivatives on EXB",
c     $ kp14fac*dEnddVxb +DCOP*dEopdVx2+DCPAR*dEpardVx2
      endif   ! icorr = 5
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

             D1F(i,POS_RA)  = dfdra
             D1F(i,POS_GAA) = dfdga
             D1F(i,POS_TA)  = dfdta
             D1F(i,POS_LA)  = dfdla
             D1F(i,POS_RB)  = dfdrb
             D1F(i,POS_GAB) = 0.0d0
             D1F(i,POS_TB)  = dfdtb
             D1F(i,POS_LB)  = dfdlb
             D1F(i,POS_GBB) = dfdgb
          D1FEXCH(i,1)  = 0.0d0
          D1FEXCH(i,2)  = 0.0d0
cccccccccccccccccccccccccccccccc
        ENDIF  !IMethod
       ENDIF  ! RA,B
c       if (i.eq.5.or.i.eq.6) then
c             write(*,101) i, EX_HF_DENSA(i), fcor
c101           format(1x, 'iGrd', i3, 2f15.9)
c       endif
c          IF(PBkp.lt.(dsmall)) then
c          write(*,*) 'PBkp = ', PBkp,'i = ',i
c          write(*,*) 'UDYN = ', UDYN,'i = ',i
c          endif
c          if(EFNA.lt.dsmall) write(*,*) 'EFNA = ', EFNA 
c          if(ndrv.eq.1) write(*,*) 'F(i) = ', F(i),'i = ',i
c          if(ndrv.eq.1) write(*,*) 'Ekp14 = ', Ekp14,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'DFAA = ', DFAA,'i = ',i
c          if(ndrv.eq.1) write(*,*) ' dEnddra = ',  dEnddra,'i = ',i
c          if(ndrv.eq.1) write(*,*) ' DSA = ',  DSA,'i = ',i
c          if(ndrv.eq.1) write(*,*) '  dg3dla = ', dg3dla,'i = ',i
c          if(ndrv.eq.1) write(*,*) '  ddffdra = ', ddffdra,'i = ',i
c          if(ndrv.eq.1) write(*,*) ' Fpar = ', Fpar,'i = ',i
c          if(ndrv.eq.1) write(*,*) ' dfpardVxa = ', dfpardVxa,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'PBkp = ', PBkp,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'Ukp = ', Ukp,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'D1Fa = ', D1F(i,POS_RA), 'i = ',i 
c          if(ndrv.eq.1) write(*,*) 'D1Fb = ', D1F(i,POS_RB), 'i = ',i 
c          if(ndrv.eq.1) write(*,*) 'gra = ', gra 
c          if(ndrv.eq.1) write(*,*) 'DFdga = ', D1F(i,POS_GAA) 
c          if(ndrv.eq.1) write(*,*) 'DFdgb = ', D1F(i,POS_GBB) 
c          if(ndrv.eq.1) write(*,*) 'DFdta = ', D1F(i,POS_TA) 
c          if(ndrv.eq.1) write(*,*) 'DFdtb = ', D1F(i,POS_TB) 
c          if(ndrv.eq.1) write(*,*) 'DFdla = ', D1F(i,POS_LA) 
c          if(ndrv.eq.1) write(*,*) 'DFdlb = ', D1F(i,POS_LB) 
c          if(ndrv.eq.1) write(*,*) 'DFdlb = ', dfdlb
c          if(ndrv.eq.1) write(*,*) 'DFdVxa = ', D1FEXCH(i,1) 
c          if(ndrv.eq.1) write(*,*) 'DFdVxb = ', D1FEXCH(i,2) 
c          if(ndrv.eq.1) write(*,*) ' Fop = ', Fop 
c          if(ndrv.eq.1) write(*,*) ' dfopdra = ', dfopdra
c          if(ndrv.eq.1) write(*,*) ' dfpardra = ', dfpardra
c          if(ndrv.eq.1) write(*,*) ' dfadra = ', dfadra
c          if(ndrv.eq.1) write(*,*) ' dNdra = ', dNdra
c          if(ndrv.eq.1) write(*,*) ' EFNA = ', EFNA 
c          if(ndrv.eq.1) write(*,*) ' dfopdrb = ', dfopdrb
c          write(*,*) ' fcor = ',  fcor,'i = ',i
c          write(*,*) ' fc1 = ', fc1
c          write(*,*) ' fc2 = ', fc2
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
c           write(*,*) ' EFNA = ', EFN2A
c          if(ndrv.eq.1) write(*,*) ' EFNB = ', EFNB
c          if(ndrv.eq.1) write(*,*) ' dfcdraa = ', dfcdraa,'i = ',i
c          if(ndrv.eq.1) write(*,*) 'LapA = ', DLapA
c          if(ndrv.eq.1) write(*,*) 'LapB[i] = ', DLapB
c          if(ndrv.eq.1) write(*,*) 'tauA = ', tauA
c          if(ndrv.eq.1) write(*,*) 'tauB[i] = ', tauB
c          if(ndrv.eq.1) write(*,*) 'rhoA[i] = ', D1, 'i = ',i
c          if(ndrv.eq.1) write(*,*) 'rhoB[i] = ', D2
 666      CONTINUE
c          if(ndrv.eq.1) write(*,*) 'VX1[i] = ', VX1
c          if(ndrv.eq.1) write(*,*) 'VX2[i] = ', VX2
c          if(ndrv.eq.1) write(*,*) 'gamB[i] = ', grb
c          write(*,*) 'Zkp = ', Zkp,'i = ',i
            ndrv=0
        ENDDO  !  NGrid
       endif  ! iterSCF
            ndrv=0
c       write(*,*) 'icountVX =  ', icountVX 
c       write(*,*) 'icountUD =  ', icountUD 
c       write(*,*) 'icountUnd =  ', icountUnd 
c       write(*,*) 'icountNA =  ', icountNA 
c       write(*,*) 'icountNB =  ', icountNB 
c       write(*,*) 'icountD =  ', icountD 
c       write(*,*) 'icountNeff =  ', icountNeff 
c       write(*,*) 'icounttt =  ', icounttt 
c       write(*,*) 'icounttp =  ', icounttp 
c       write(*,*) 'YSUM =  ', YSUM(NCENTER), 'CENTER  ', NCENTER 
c       write(*,*) 'PKSUM =  ', PKSUM(NCENTER), 'CENTER  ', NCENTER 
c       write(*,*) 'Ebatch =  ', Ebatch 
c       write(*,*) 'AVNEF =  ', AVNEF(NCENTER), 'CENTER  ', NCENTER 
c                TOTSUM = 0.0d0
c               IF(NCENTER.eq.NAtoms) then
c                 call VRtrace(TOTSUM, YSUM, NAtoms)
c           write(*,*) 'TOTSUM  =  ', TOTSUM 
c               endif
c
c      CALL ftnQFree(pEX_HF_DENSA)
c      CALL ftnQFree(pEX_HF_DENSB)
c      CALL ftnQFree(pFIT_DENSA)
c      CALL ftnQFree(pFIT_DENSB)
      RETURN
      END
c
      SUBROUTINE brsc14(DN,DFN,GR,DLAP,tau,VX,UX,EFN,EFN2,dNdr,dNdg,
     $ dNdt,dNdl,
     $ dNdVx,dM1dr,dM2dr,dM1dt,dM2dt,dM1dl,dM2dl,dM1dg,dM2dg,dM1dVx,
     $ dM2dVx,dg3dn,dg3dt,dg3dg,dg4dn,dg4dt,dg4dg,SM1,SM2,DS,
     $ detol,dsmall,ndrv)
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
c      parameter(dsmall=1.0D-13)
       parameter(dtol2=1.0D-08)
       parameter(dtol4=1.0D-06)
       parameter(dtol3=1.0D-08)
       parameter(dtol=1.0D-07)
       parameter(delt=0.070d0)  ! smoothening of EFN
        parameter(detol0=1.0D-08)
c      parameter(delt=0.0050d0)  ! smoothening of EFN

       REAL*8 xx1,xx2,yy,acc,X_y,EFN,DN,GR,DLAP,tau,SM1,SM2,DS
       INTEGER iter
ccccccccccccfor the numerical version onlyccccc
c      external FF
c      external rtsafe
cccccccccccccccccccccccccccccccccccccccccccccc
          Tol = 1.0d-08
           smth1 = 2.0d0  - delt
           smth2 = 2.0d0  + delt
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
          dg4dn = zero
          dg3dt = zero
          dg4dt = zero
          dg3dg = zero
          dg4dg = zero
          EFN = 0.0d0
          EFN1 = 0.0d0
          EFN2 = 0.0d0
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
          DS = dsmall
c Initializations end
c
c    First calculate X_y and the energy density
c
           DN13 = DN**third
           DN43= DN13*DN
           DN53= DN*(DN13*DN13)
           pi13=pi**third
           pi23 = pi13*pi13
           if(tau.lt.dsmall) tau=dsmall
           t1 = four*pi*DN
           tq1 = t1*DN
           DS = tau-GR/(DN*four)
           QS = DLAP/six -DS/three
c  To keep in mind that VX = Ux*DN:
          IF(tq1.gt.detol) then
             yy = -three*QS*UX/tq1
          endif
          if(abs(yy).lt.dsmall) ndrv =0
cccccccccccfor the numerical version onlyccccccccccccccc
c           xx1=detol
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
         IF(yy.lt.(-dsmall)) then
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
        endif
            IF(abs(yy).lt.dsmall) X_y=2.00001d0
c   end of calculating X(y)
c  115     CONTINUE   ! uncomment only for the numerical
c            twthr=2/3, trihf= 3/2
             t01 = sqrt(two)
             t02 = sqrt(three)
             tn0 = (twthr**trihf)*pi
             tn1 = DN*DN
             tn2 = sqrt(DN)
             tn4 = exp(X_y)
             tn5 = exp(-X_y)
             tu1 = X_y*X_y
             ts12 = dsmall
             tn10 = 0.0d0 
             ts10 = dsmall
             xxxq = X_y*QS
        IF(X_y.gt.detol) then
             ts12 = X_y+4.0d0/X_y-(1.0d0+4.0d0/X_y)*tn5
          IF(abs(xxxq).gt.detol)
     $      tn10 = (X_y - two)*two/(xxxq*three)
           if(tn10.gt.dsmall) then
             ts10 = tn10*DN/4.0d0
             tn11 = sqrt(tn10)
             ts11 = sqrt(ts10)
             ALF3 = tn10*tn11
c         write(*,*) 'ALF3 =  ', ALF3 
             EFN1 = ALF3*pi*tn4*DN*DN*sqrt(DN)
             EFN2 = EFN1
             IF(EFN2.gt.1.0d0) EFN2=1.0d0
c          ggg= EFN1-2.0d0-delt
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
           EFN = 0.0d0
           EFN2=0.00d0
           SM1 = 0.0d0
           SM2 = 0.0d0
        endif   ! tn10
c
c   Now prepare for the SCF potential:
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
            dg4dn = GR/(4.0d0*DN)/DN
            dg3dt = one/(three*DN)
            dg4dt = 1.0d0
            dg3dg = -one/(four*DN)/(three*DN)
            dg4dg = -1.0d0/(4.0d0*DN)
            else
              dg3dn = 0.0d0
              dg3dt = 0.0d0
              dg3dg = 0.0d0
              dg4dn = 0.0d0
              dg4dt = 0.0d0
              dg4dg = 0.0d0
            endif
            dNdn = dNg1*5.0d0*EFN1/(two*DN)
            tt12 = (X_y-two)/(X_y*QS)
            dNdx = zero
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
            if(tt13.gt.(detol)) then
            t18 = sqrt(tt13)
            t19 = exp(X_y)
            t23 = t8*t8
          dM2dx = pi*t1*DN*t3*t4/t8/X_y/QS
     $    *(t18*t19)/QS*(13.0d0*t8+60.0d0+t23-24.0d0*X_y)/27.0d0
        endif
           dM2dQ = -SM2*5.0d0/(two*QS)
           dM2dn=SM2*7.0d0/(two*DN)
c    here dfdt == dNdt, etc.
         dNdr = dNdn + dNdx*dxdy*dydn +dNdQ*dQdn
         dNdt = dNdx*dxdy*dydt + dNdQ*dQdt
         dNdl = dNdx*dxdy*dydl+ dNdQ*dQdl
         dNdg = dNdx*dxdy*dydg + dNdQ*dQdg
         dNdVx = dNdx*dxdy*dydVx
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
          EFN = 0.0d0
       endif  ! X_y
c      endif  ! Tolbig yy
c     ENDIF  ! DN and  tau
      RETURN
      END
c
       SUBROUTINE B14dync(D1,VX1,EFN1,ASSA,dNdra,
     $   dNdga,dNdta,dNdla,dNdVxa,DDA,dg3dra,dg3dta,dg3dga,
     $   dfopdra,dfopdga,dfopdta,dfopdla,dfopdVxa,
     $   dASSAdra,dASSAdga,dASSAdta,dASSAdla,dASSAdVxa,
     $   dASSAdrb,dASSAdgb,dASSAdtb,dASSAdlb,dASSAdVxb ,
     $   D2,VX2,EFN2,ASSB,
     $   dNdrb,dNdgb,dNdtb,dNdlb,dNdVxb,DDB,dg3drb,dg3dtb,dg3dgb,
     $   dfopdrb,dfopdgb,dfopdtb,dfopdlb,dfopdVxb,
     $   dASSBdra,dASSBdga,dASSBdta,dASSBdla,dASSBdVxa,
     $   dASSBdrb,dASSBdgb,dASSBdtb,dASSBdlb,dASSBdVxb,fopp,
c
     $   dEopdn1,dEopdn2,dEopdg1,dEopdg2,dEopdt1,dEopdt2,
     $   dEopdl1,dEopdl2,dEopdVx1,dEopdVx2,dEpardn1,dEpardn2,
     $   dEpardg1,dEpardg2,dEpardt1,dEpardt2,dEpardl1,dEpardl2,
     $   dEpardVx1,dEpardVx2,
c
     $   duopdn1,duopdn2,duopdg1,duopdg2,duopdt1,duopdt2,
     $   duopdl1,duopdl2,duopdVx1,duopdVx2,dupardn1,dupardn2,  
     $   dupardg1,dupardg2,dupardt1,dupardt2,dupardl1,dupardl2,
     $   dupardVx1,dupardVx2,UOPP,UPAR,UP1,UP2,ULAOPP,ULAP1,ULAP2,
     $   ULAPAR,EOPP,EP1,EP2,EPAR,detol,dsmall,ndrv,qlamb,NA,NB,icorr)
       IMPLICIT REAL*8(A-H,O-Z)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
       parameter(five=5.0d0,six=6.0d0,three=3.0d0)
       parameter(third=1.0d0/3.0d0)
       parameter(twthr=2.d0/3.d0, trihf=3.0d0/2.0d0)
       parameter(pi=3.141592653589793d0)
       parameter(Tolbig=1.0d+09)
       parameter(COPP= 0.630D0, CPAR=0.880D0)
       ZOPP= 0.0d0 
       Z1 = 0.0d0 
       Z2 = 0.0d0 
       UOPP = 0.0d0
       ULAOPP =0.0d0
       EOPP = 0.0d0
       UP1 = 0.0d0
       ULAP1= 0.0d0
       EP1 = 0.0d0
       UP2 = 0.0d0
       ULAP2= 0.0d0
       EP2 = 0.0d0
       UPAR= 0.0d0 
       ULAPAR= 0.0d0 
       EPAR = 0.0d0
       duopdn1=0.0D0
       duopdn2=0.0D0
       duopdg1=0.0D0
       duopdg2=0.0D0
       duopdt1=0.0D0
       duopdt2=0.0D0
       duopdl1=0.0D0
       duopdl2=0.0D0
       duopdVx1=0.0D0
       duopdVx2=0.0D0
       dZZdzop=0.0d0
c
       dupardn1=0.0D0
       dupardn2=0.0D0
       dupardg1=0.0D0
       dupardg2=0.0D0
       dupardt1=0.0D0
       dupardt2=0.0D0
       dupardl1=0.0D0
       dupardl2=0.0D0
       dupardVx1=0.0D0
       dupardVx2=0.0D0
c
       du1dn1=0.0D0
       du1dg1=0.0D0
       du1dt1=0.0D0
       du1dl1=0.0D0
       du1dVx1=0.0D0
c
       du2dn2=0.0D0
       du2dg2=0.0D0
       du2dt2=0.0D0
       du2dl2=0.0D0
       du2dVx2=0.0D0
c
c       write(6,*)"Neff RB", dNdrb
c       write(6,*)"Neff GBB", dNdgb
c       write(6,*)"A2 RB", dg3drb
c       write(6,*)"A2 GBB",dg3dgb
c       write(6,*)"A1 RB", dASSBdrb
c       write(6,*)"A1 GBB",dASSBdgb
c       write(6,*)"alpha A1 RB", dASSAdrb
c       write(6,*)"alpha A1 GBB",dASSAdgb
        IF(ndrv.eq.1) then
c       D1 = 0.199089068444670d0  + 0.0001d0
c       D1 = 0.199089068444670d0  - 0.0001d0
        endif
             tt1 = 0.0d0
             tt2 = 0.0d0
             ttd = 0.0d0
             tt3 = 0.0d0
             tla3 = 0.d0
             tlad = 0.0d0
             tlog2 = 0.0d0
             tlog3 = 0.0d0
             tlog4 = 0.0d0
             ZOPP = 0.0d0
       IF((VX1.gt.-detol).or.(VX2.gt.-detol)) GO TO 666
c        IF((D1.gt.detol).or.(D2.gt.detol)) then
           IF((D1.gt.detol).and.(D2.gt.detol)) then
            IF((VX1.lt.(-detol)).and.(VX2.lt.(-detol))) then
            tt1 = -EFN1*D1/VX1
            tt2 = -EFN2*D2/VX2
            ZOPP = COPP*(tt1 + tt2)
         IF(ndrv.eq.1) then
c          ZOPP = 1.415326759106320d0 + 0.00010d0 
c          ZOPP = 1.415326759106320d0 - 0.00010d0 
         endif
             ttd = 1.0d0+ZOPP
             tlad = 1.0d0+qlamb*ZOPP
             tt3 = ZOPP*ZOPP*ZOPP/ttd
c         if(ndrv.eq.1) write(*,*) 'tt3 = ', tt3 
             tla3 = ZOPP*ZOPP*ZOPP/tlad
             IF(ZOPP.gt.detol) then
              tlog2 = log(ttd)/ZOPP
              BOPP = (1.0d0-tlog2)
             dBopdZop = tlog2/ZOPP -1.0d0/(ttd*ZOPP)
             else
              BOPP = 1.0d0
              tlog2 = 1.0d0
             dBopdZop = 0.0d0
             endif
c            UOPP = -0.80d0*(1.0d0-fopp)*D1*D2*tt3
             UOPP = -0.80d0*D1*D2*tt3 + 0.80d0*fopp*D1*D2*tt3
          ULAOPP = -0.80d0*D1*D2*qlamb*tla3+0.80d0*fopp*D1*D2*qlamb*tla3
             EOPP = -0.80d0*(1.0d0-fopp)*D1*D2*ZOPP*ZOPP*BOPP
ccccc Now the derivatives of Uopp, in the order of: 
ccccc D1,gr1,tau1,lap1,Vx1 
             dZZdzop = 3.0d0*ZOPP*ZOPP/ttd - tt3/ttd
ccc
          IF(EFN1.gt.detol) then
             dzopdn1 = -COPP*(dNdra*D1 + EFN1)/VX1
           else 
              dzopdn1=0.0d0
           endif
          IF(EFN2.gt.detol) then
             dzopdn2 = -COPP*(dNdrb*D2+EFN2)/VX2
          else
             dzopdn2 = 0.0d0
          endif
c
       dduop = D2*tt3 + D1*D2*dzopdn1*dZZdzop
       duoppdn1 = dfopdra*D1*D2*tt3-(1.0d0-fopp)*dduop 
       duopdn1=duoppdn1*0.80d0
c
       duoppdn2 = dfopdrb*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*tt3+D1*D2*dzopdn2*dZZdzop)
       duopdn2=duoppdn2*0.80d0
cccc
            dzopdg1 = -COPP*(dNdga*D1)/VX1
            dzopdg2 = -COPP*(dNdgb*D2)/VX2
       duoppdg1 = dfopdga*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdg1*dZZdzop)
       duopdg1=duoppdg1*0.80d0
c
       duoppdg2 = dfopdgb*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdg2*dZZdzop)
       duopdg2=duoppdg2*0.80d0
cccc
            dzopdt1 = -COPP*(dNdta*D1)/VX1
            dzopdt2 = -COPP*(dNdtb*D2)/VX2
       duoppdt1 = dfopdta*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdt1*dZZdzop)
       duopdt1=duoppdt1*0.80d0
c
       duoppdt2 = dfopdtb*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdt2*dZZdzop)
       duopdt2=duoppdt2*0.80d0
cccc
            dzopdl1 = -COPP*(dNdla*D1)/VX1
            dzopdl2 = -COPP*(dNdlb*D2)/VX2
       duoppdl1 = dfopdla*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdl1*dZZdzop)
       duopdl1=duoppdl1*0.80d0
c
       duoppdl2 = dfopdlb*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdl2*dZZdzop)
       duopdl2=duoppdl2*0.80d0
cccc
         dzopdVx1 = -COPP*(-EFN1*D1/(VX1**2) +D1*dNdVxa/VX1)
         dzopdVx2 = -COPP*(-EFN2*D2/(VX2**2) +D2*dNdVxb/VX2)
c
       duoppdVx1 = dfopdVxa*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdVx1*dZZdzop)
       duopdVx1=duoppdVx1*0.80d0
c
       duoppdVx2 = dfopdVxb*D1*D2*tt3
     $         -(1.0d0-fopp)*(D1*D2*dzopdVx2*dZZdzop)
       duopdVx2 =duoppdVx2*0.80d0
c          endif   !   uopp
        endif     !  D1 and D2
          endif  ! VX1  and VX2
c           IF(UOPP.gt.-dsmall) UOPP = 0.0d0 
c          IF(icorr.eq.5) then
           IF((icorr.eq.5).or.(icorr.eq.1)) then
       ddeopr11 = -dfopdra*D1*D2*ZOPP*ZOPP*BOPP
       ddeopr12 =  D2*ZOPP*ZOPP*BOPP+2.0d0*D1*D2*ZOPP*dzopdn1*BOPP
       ddeopr13 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdn1
       dEopdn1 = -0.80d0*(ddeopr11+(1.0d0-fopp)*(ddeopr12+ddeopr13))
c
       ddeopr21 = -dfopdrb*D1*D2*ZOPP*ZOPP*BOPP
       ddeopr22 =  D1*ZOPP*ZOPP*BOPP+2.0d0*D1*D2*ZOPP*dzopdn2*BOPP
       ddeopr23 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdn2
       dEopdn2 = -0.80d0*(ddeopr21+(1.0d0-fopp)*(ddeopr22+ddeopr23))
c
       ddeopg11 = -dfopdga*D1*D2*ZOPP*ZOPP*BOPP
       ddeopg12 =  2.0d0*D1*D2*ZOPP*dzopdg1*BOPP
       ddeopg13 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdg1
       dEopdg1 = -0.80d0*(ddeopg11+(1.0d0-fopp)*(ddeopg12+ddeopg13))
c
       ddeopg21 = -dfopdgb*D1*D2*ZOPP*ZOPP*BOPP
       ddeopg22 =  2.0d0*D1*D2*ZOPP*dzopdg2*BOPP
       ddeopg23 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdg2
       dEopdg2 = -0.80d0*(ddeopg21+(1.0d0-fopp)*(ddeopg22+ddeopg23))
c
       ddeopt11 = -dfopdta*D1*D2*ZOPP*ZOPP*BOPP
       ddeopt12 =  2.0d0*D1*D2*ZOPP*dzopdt1*BOPP
       ddeopt13 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdt1
       dEopdt1 = -0.80d0*(ddeopt11+(1.0d0-fopp)*(ddeopt12+ddeopt13))
c
       ddeopt21 = -dfopdtb*D1*D2*ZOPP*ZOPP*BOPP
       ddeopt22 =  2.0d0*D1*D2*ZOPP*dzopdt2*BOPP
       ddeopt23 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdt2
       dEopdt2 = -0.80d0*(ddeopt21+(1.0d0-fopp)*(ddeopt22+ddeopt23))
c
       ddeopl11 = -dfopdla*D1*D2*ZOPP*ZOPP*BOPP
       ddeopl12 =  2.0d0*D1*D2*ZOPP*dzopdl1*BOPP
       ddeopl13 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdl1
       dEopdl1 = -0.80d0*(ddeopl11+(1.0d0-fopp)*(ddeopl12+ddeopl13))
c
       ddeopl21 = -dfopdlb*D1*D2*ZOPP*ZOPP*BOPP
       ddeopl22 =  2.0d0*D1*D2*ZOPP*dzopdl2*BOPP
       ddeopl23 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdl2
       dEopdl2 = -0.80d0*(ddeopl21+(1.0d0-fopp)*(ddeopl22+ddeopl23))
c
       ddeopX11 = -dfopdVxa*D1*D2*ZOPP*ZOPP*BOPP
       ddeopX12 =  2.0d0*D1*D2*ZOPP*dzopdVx1*BOPP
       ddeopX13 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdVx1
       dEopdVx1 = -0.80d0*(ddeopX11+(1.0d0-fopp)*(ddeopX12+ddeopX13))
c
       ddeopX21 = -dfopdVxb*D1*D2*ZOPP*ZOPP*BOPP
       ddeopX22 =  2.0d0*D1*D2*ZOPP*dzopdVx2*BOPP
       ddeopX23 = D1*D2*ZOPP*ZOPP*dBopdZop*dzopdVx2
       dEopdVx2 = -0.80d0*(ddeopX21+(1.0d0-fopp)*(ddeopX22+ddeopX23))
c
        endif ! icorr 5
cccccccc Now the parallel-spin component
              dzz1dz1=0.0d0
             pp1 = 0.0d0 
             pla1 = 0.0d0
             pp3 = 0.0d0 
             pp4 = 0.0d0 
             pp7 = 0.0d0 
             pp7b = 0.0d0 
             pp8 = 0.0d0 
             pp8b = 0.0d0 
             pp9= 0.0d0 
             pp9b= 0.0d0 
             pp10= 0.0d0
             pla10 = 0.0d0
             pp11 = 0.0d0
             pla11 = 0.0d0
ccccc D1,gr1,tau1,lap1,Vx1 
c
           IF(D1.gt.detol) then
             IF(VX1.le.(-detol)) then
             Z1 = -2.0d0*CPAR*EFN1*D1/VX1
             pp1 = 1.0d0+Z1/2.0d0
             pla1 = 1.0d0+Z1*qlamb/2.0d0
             pp3 = 5.0d0+2.0d0*Z1
             pp4 = pp1**2
             pp7 = Z1*Z1
             pp8 = pp7*pp7
             pp9= pp8*Z1
              app1 = abs(pp1)
           IF(app1.gt.dsmall) then
             pp10=pp9/pp1
           else
             pp10 =  0.0d0
           endif
              ala1 = abs(pla1)
           IF(ala1.gt.dsmall) then
             pla10=pp9/pla1
           else
             pla10 =  0.0d0
           endif
            if(pp4.gt.dsmall) then
             dzz1dz1 = (Z1**4)*pp3/pp4
            else
             dzz1dz1 = 0.0d0
            endif
            IF((abs(Z1).gt.detol).and.(pp1.gt.dsmall)) then
             tlog3 = 2.0d0*log(pp1)/Z1
             tlog5 = tlog3/Z1
            dB1dZ1 =  -1.0d0/pp1/Z1 + tlog5
            else
             tlog3 = 1.0d0
             tlog5 = 0.0d0
             dB1dZ1 = 0.0d0
            endif
           UP1 = -0.0050d0*(DDA-ASSA)*3.0d0*D1*D1*pp10
           ULAP1 = -0.0050d0*(DDA-ASSA)*3.0d0*qlamb*D1*D1*pla10
           B11 = (1.0d0-tlog3)
           EP1 = -0.030d0*(DDA-ASSA)*D1*D1*pp8*B11
c
            dz1dn1 = -2.0d0*CPAR*(D1*dNdra+EFN1)/VX1
            dz1dg1 = -2.0d0*CPAR*D1*dNdga/VX1
            dz1dt1 = -2.0d0*CPAR*D1*dNdta/VX1
            dz1dl1 = -2.0d0*CPAR*D1*dNdla/VX1
            dz1dVx1 = 2.0d0*CPAR*D1*(EFN1/(VX1*VX1) -dNdVxa/VX1)
c
         du1dnn1 = (dg3dra-dASSAdra)*3.0d0*D1*D1*pp10+(DDA-ASSA)
     $            *(pp10*6.0d0*D1+3.0d0*D1*D1*dz1dn1*dzz1dz1) 
         du1dn1 = -0.0050d0*du1dnn1
c
         du1dgg1 = (dg3dga-dASSAdga)*3.0d0*D1*D1*pp10
     $           +(DDA-ASSA)*3.0d0*D1*D1*dz1dg1*dzz1dz1
         du1dg1 = -0.0050d0*du1dgg1
c
         du1dtt1 = (dg3dta-dASSAdta)*3.0d0*D1*D1*pp10
     $           +(DDA-ASSA)*3.0d0*D1*D1*dz1dt1*dzz1dz1
         du1dt1 = -0.0050d0*du1dtt1
c
         du1dll1 = -dASSAdla*3.0d0*D1*D1*pp10
     $           +(DDA-ASSA)*3.0d0*D1*D1*dz1dl1*dzz1dz1
         du1dl1 = -0.0050d0*du1dll1
c
         du1dvv1 = -dASSAdVxa*3.0d0*D1*D1*pp10
     $           +(DDA-ASSA)*3.0d0*D1*D1*dz1dVx1*dzz1dz1
         du1dVx1 = -0.0050d0*du1dvv1
c
c           if(icorr.eq.5) then
            if((icorr.eq.5).or.(icorr.eq.1)) then
c
         du1dn11 = (dg3dra-dASSAdra)*D1*D1*pp8*B11 
         du1dn12 = 2.0d0*(DDA-ASSA)*D1*pp8*B11
         du1dn13 = (DDA-ASSA)*D1*D1*pp8*dB1dZ1*dz1dn1 
         du1dn14 = 4.0d0*(DDA-ASSA)*D1*D1*Z1*Z1*Z1*B11*dz1dn1
         dE1dn1 = -0.030d0*(du1dn11+du1dn12+ du1dn13+du1dn14)
c
         du1dg11 = (dg3dga-dASSAdga)*D1*D1*pp8*B11
         du1dg13 = (DDA-ASSA)*D1*D1*(pp8*dB1dZ1*dz1dg1
     $            + 4.0d0*Z1*Z1*Z1*B11*dz1dg1)
         dE1dg1 = -0.030d0*(du1dg11+du1dg13)
c
         du1dt11 = (dg3dta-dASSAdta)*D1*D1*pp8*B11
         du1dt13 = (DDA-ASSA)*D1*D1*(pp8*dB1dZ1*dz1dt1
     $            + 4.0d0*Z1*Z1*Z1*B11*dz1dt1)
         dE1dt1 = -0.030d0*(du1dt11+du1dt13)
c
         du1dl11 = (-dASSAdla)*D1*D1*pp8*B11
         du1dl13 = (DDA-ASSA)*D1*D1*(pp8*dB1dZ1*dz1dl1
     $            + 4.0d0*Z1*Z1*Z1*B11*dz1dl1)
         dE1dl1 = -0.030d0*(du1dl11+du1dl13)
c
         du1dv11 = (-dASSAdVxa)*D1*D1*pp8*B11
         du1dv13 = (DDA-ASSA)*D1*D1*(pp8*dB1dZ1*dz1dVx1
     $            + 4.0d0*Z1*Z1*Z1*B11*dz1dVx1)
         dE1dVx1 = -0.030d0*(du1dv11+du1dv13)
c
       endif  ! icorr 5
        endif    ! Vx1
       endif  !  D1
c
ccccc D2,gr2,tau2,lap2,Vx2 
               dzz2dz2 =0.0d0
         IF(D2.gt.detol) then
           IF(VX2.le.(-detol)) then
             Z2 = -2.0d0*CPAR*EFN2*D2/VX2
             pp2 = 1.0d0+Z2/2.0d0
             pla2 = 1.0d0+Z2*qlamb/2.0d0
             pp5 = 5.0d0+2.0d0*Z2
             pp6 = pp2**2
             pp7b = Z2*Z2
             pp8b = pp7b*pp7b
             pp9b = pp8b*Z2
             app2 = abs(pp2)
           IF(app2.gt.dsmall) then
             pp11=pp9b/pp2
           else
             pp11 = 0.0d0
           endif
               ala2 = abs(pla2)
             IF(ala2.gt.dsmall) then
             pla11=pp9b/pla2
           else
             pla11 = 0.0d0
           endif
            IF(pp6.gt.dsmall) then
             dzz2dz2 = (Z2**4)*pp5/pp6
            else
             dzz2dz2 = 0.0d0
            endif
c            write(6,*)"FZ_Z",dzz2dz2
            IF((abs(Z2).gt.detol).and.(pp2.gt.dsmall)) then
             tlog4 = 2.0d0*log(pp2)/Z2
             tlog6= tlog4/Z2
             dB2dZ2 =  -1.0d0/pp2/Z2 + tlog6 
            else
             tlog4 = 1.0d0
             tlog6 = 0.0d0
             dB2dZ2 = 0.0d0
            endif
           UP2 = -0.0050d0*(DDB-ASSB)*3.0d0*D2*D2*pp11
           ULAP2 = -0.0050d0*(DDB-ASSB)*3.0d0*qlamb*D2*D2*pla11
           B22 = (1.0d0-tlog4)
           EP2 = -0.030d0*(DDB-ASSB)*D2*D2*pp8b*B22
            dz2dn2 = -2.0d0*CPAR*(D2*dNdrb+EFN2)/VX2
            dz2dg2 = -2.0d0*CPAR*D2*dNdgb/VX2
            dz2dt2 = -2.0d0*CPAR*D2*dNdtb/VX2
            dz2dl2 = -2.0d0*CPAR*D2*dNdlb/VX2
            dz2dVx2 = 2.0d0*CPAR*D2*(EFN2/(VX2*VX2) -dNdVxb/VX2)
c            write(6,*)"ZBB", z2
c            write(6,*)"ZBB_RB", dz2dn2
c            write(6,*)"ZBB_GBB",dz2dg2
c
            IF(D1.gt.detol) then
             IF((icorr.eq.5).or.(icorr.eq.1)) then
              dE1dn1 = dE1dn1+0.030d0*dASSBdra*D2*D2*pp8b*B22
              dE1dg1 = dE1dg1+0.030d0*dASSBdga*D2*D2*pp8b*B22
              dE1dt1 = dE1dt1+0.030d0*dASSBdta*D2*D2*pp8b*B22
              dE1dl1 = dE1dl1+0.030d0*dASSBdla*D2*D2*pp8b*B22
              dE1dVx1 = dE1dVx1+0.030d0*dASSBdVxa*D2*D2*pp8b*B22
c             else
              du1dn1 = du1dn1+0.0050d0*dASSBdra*3.0d0*D2*D2*pp11
              du1dg1 = du1dg1+0.0050d0*dASSBdga*3.0d0*D2*D2*pp11
              du1dt1 = du1dt1+0.0050d0*dASSBdta*3.0d0*D2*D2*pp11
              du1dl1 = du1dl1+0.0050d0*dASSBdla*3.0d0*D2*D2*pp11
              du1dVx1 = du1dVx1+0.0050d0*dASSBdVxa*3.0d0*D2*D2*pp11
            endif   ! icorr 5
            endif   !  D1
c
         du2dnn2 = (dg3drb-dASSBdrb)*3.0d0*D2*D2*pp11+(DDB-ASSB)
     $            *(pp11*6.0d0*D2+3.0d0*D2*D2*dz2dn2*dzz2dz2)
         du2dn2 = -0.0050d0*du2dnn2
c
         du2dgg2 = (dg3dgb-dASSBdgb)*3.0d0*D2*D2*pp11
     $           +(DDB-ASSB)*3.0d0*D2*D2*dz2dg2*dzz2dz2
         du2dg2 = -0.0050d0*du2dgg2
c
         du2dtt2 = (dg3dtb-dASSBdtb)*3.0d0*D2*D2*pp11
     $           +(DDB-ASSB)*3.0d0*D2*D2*dz2dt2*dzz2dz2
         du2dt2 = -0.0050d0*du2dtt2
c
         du2dll2 = -dASSBdlb*3.0d0*D2*D2*pp11
     $           +(DDB-ASSB)*3.0d0*D2*D2*dz2dl2*dzz2dz2
         du2dl2 = -0.0050d0*du2dll2
c
         du2dvv2 = -dASSBdVxb*3.0d0*D2*D2*pp11
     $           +(DDB-ASSB)*3.0d0*D2*D2*dz2dVx2*dzz2dz2
         du2dVx2 = -0.0050d0*du2dvv2
c
            if((icorr.eq.5).or.(icorr.eq.1)) then
         du2dn21 = (dg3drb-dASSBdrb)*D2*D2*pp8b*B22
         du2dn22 = 2.0d0*(DDB-ASSB)*D2*pp8b*B22
         du2dn23 = (DDB-ASSB)*D2*D2*pp8b*dB2dZ2*dz2dn2
         du2dn24 = 4.0d0*(DDB-ASSB)*D2*D2*Z2*Z2*Z2*B22*dz2dn2
         dE2dn2 = -0.030d0*(du2dn21+du2dn22+ du2dn23+du2dn24)
c
         du2dg21 = (dg3dgb-dASSBdgb)*D2*D2*pp8b*B22
         du2dg23 = (DDB-ASSB)*D2*D2*(pp8b*dB2dZ2*dz2dg2
     $            + 4.0d0*Z2*Z2*Z2*B22*dz2dg2)
         dE2dg2 = -0.030d0*(du2dg21+du2dg23)
c
         du2dt21 = (dg3dtb-dASSBdtb)*D2*D2*pp8b*B22
         du2dt23 = (DDB-ASSB)*D2*D2*(pp8b*dB2dZ2*dz2dt2
     $            + 4.0d0*Z2*Z2*Z2*B22*dz2dt2)
         dE2dt2 = -0.030d0*(du2dt21+du2dt23)
c
         du2dl21 = (-dASSBdlb)*D2*D2*pp8b*B22
         du2dl23 = (DDB-ASSB)*D2*D2*(pp8b*dB2dZ2*dz2dl2
     $            + 4.0d0*Z2*Z2*Z2*B22*dz2dl2)
         dE2dl2 = -0.030d0*(du2dl21+du2dl23)
c
         du2dv21 = (-dASSBdVxb)*D2*D2*pp8b*B22
         du2dv23 = (DDB-ASSB)*D2*D2*(pp8b*dB2dZ2*dz2dVx2
     $            + 4.0d0*Z2*Z2*Z2*B22*dz2dVx2)
         dE2dVx2 = -0.030d0*(du2dv21+du2dv23)
c
           endif  ! icorr 5

c          write(6,*)"ECBB_RB",du2dn2
c          write(6,*)"ECBB_GBB",du2dg2

         IF(D1.gt.detol) then
           IF((icorr.eq.5).or.(icorr.eq.1)) then
            dE2dn2 = dE2dn2+0.030d0*dASSAdrb*D1*D1*pp8*B11
            dE2dg2 = dE2dg2+0.030d0*dASSAdgb*D1*D1*pp8*B11
            dE2dt2 = dE2dt2+0.030d0*dASSAdtb*D1*D1*pp8*B11
            dE2dl2 = dE2dl2+0.030d0*dASSAdlb*D1*D1*pp8*B11
            dE2dVx2 = dE2dVx2+0.030d0*dASSAdVxb*D1*D1*pp8*B11
           endif
            du2dn2 = du2dn2+0.0050d0*dASSAdrb*3.0d0*D1*D1*pp10
            du2dg2 = du2dg2+0.0050d0*dASSAdgb*3.0d0*D1*D1*pp10
            du2dt2 = du2dt2+0.0050d0*dASSAdtb*3.0d0*D1*D1*pp10
            du2dl2 = du2dl2+0.0050d0*dASSAdlb*3.0d0*D1*D1*pp10
            du2dVx2 = du2dVx2+0.0050d0*dASSAdVxb*3.0d0*D1*D1*pp10
          endif

c          write(6,*)"ECAA_RB",0.0050d0*dASSAdrb*3.0d0*D1*D1*pp10
c          write(6,*)"ECAA_GBB",0.0050d0*dASSAdgb*3.0d0*D1*D1*pp10
         endif  ! VX2
        endif   !   D2
c        IF(NA.eq.1) UP1=0.0d0
c        IF(NA.eq.1) EP1=0.0d0
c        IF(NB.eq.1) UP2=0.0d0
c        IF(NB.eq.1) EP2=0.0d0
            UPAR = UP1+UP2
            EPAR = EP1+EP2
            ULAPAR = ULAP1+ULAP2
          dupardn1 = du1dn1
          dupardg1 = du1dg1
          dupardt1 = du1dt1
          dupardl1 = du1dl1
          dupardVx1 = du1dVx1
c
          dupardn2 = du2dn2
          dupardg2 = du2dg2
          dupardt2 = du2dt2
          dupardl2 = du2dl2
          dupardVx2 = du2dVx2
c      endif !D1 or D2 
c         if(ndrv.eq.1) write(*,*) 'D1 = ', D1 
c         if(ndrv.eq.1) write(*,*) 'D2 = ', D2 
c         if(ndrv.eq.1) write(*,*) 'EOPP = ', EOPP
c         if(ndrv.eq.1) write(*,*) 'EP2 = ', EP2
c         if(ndrv.eq.1) write(*,*) 'EPAR = ', EPAR
c         if(ndrv.eq.1) write(*,*) 'ZOPP = ', ZOPP
c         if(ndrv.eq.1) write(*,*) 'dfcopdVxa = ', dfopdVxa 
c         if(ndrv.eq.1) write(*,*) 'dfcopdVxb = ', dfopdVxb 
c         if(ndrv.eq.1) write(*,*) 'dzopdt1 = ', dzopdt1
c         if(ndrv.eq.1) write(*,*) 'fcor = ', fopp 
c         if(ndrv.eq.1) write(*,*) 'BOPP = ', BOPP 
c         if(ndrv.eq.1) write(*,*) 'dBopdZop = ', dBopdZop 
c         if(ndrv.eq.1) write(*,*) 'duopdl1 = ', duopdl1
c        if(ndrv.eq.1) write(*,*) 'duopdl2 = ', duopdl2
c         if(ndrv.eq.1) write(*,*) 'duopdg1 = ', duopdg1
c         if(ndrv.eq.1) write(*,*) 'duopdg2 = ', duopdg2
c         if(ndrv.eq.1) write(*,*) 'duopdVx1 = ', duopdVx1
c         if(ndrv.eq.1) write(*,*) 'duopdVx2 = ', duopdVx2
c         if(ndrv.eq.1) write(*,*) 'duopdt1 = ', duopdn1
c         if(ndrv.eq.1) write(*,*) 'duopdt2 = ', duopdn2
c         if(ndrv.eq.1) write(*,*) 'dupardn1 = ', dupardn1
c         if(ndrv.eq.1) write(*,*) 'dupardn2 = ', dupardn2
c         if(ndrv.eq.1) write(*,*) 'dupardg1 = ', dupardg1
c         if(ndrv.eq.1) write(*,*) 'dupardg2 = ', dupardg2
c         if(ndrv.eq.1) write(*,*) 'dupardt1 = ', dupardt1
c         if(ndrv.eq.1) write(*,*) 'dupardt2 = ', dupardt2
c         if(ndrv.eq.1) write(*,*) 'dz1dn1 = ', dz1dn1 
c         if(ndrv.eq.1) write(*,*) 'dz2dn2 = ', dz2dn2 
c         if(ndrv.eq.1) write(*,*) 'pp10 = ', pp10 
c         if(ndrv.eq.1) write(*,*) 'pp11 = ', pp11 
c         if(ndrv.eq.1) write(*,*) 'dg3dra = ', dg3dra 
c         if(ndrv.eq.1) write(*,*) 'dg3drb = ', dg3drb 
c         if(ndrv.eq.1) write(*,*) 'DDA = ', DDA 
c         if(ndrv.eq.1) write(*,*) 'DDB = ', DDB 
c         if(ndrv.eq.1) write(*,*) 'dASSAdra = ', dASSAdra 
c         if(ndrv.eq.1) write(*,*) 'dASSBdrb = ', dASSBdrb 
c         if(ndrv.eq.1) write(*,*) 'B11 = ', B11 
c         if(ndrv.eq.1) write(*,*) 'B22 = ', B22 
c         if(ndrv.eq.1) write(*,*) 'dB1dZ1 = ', dB1dZ1 
c         if(ndrv.eq.1) write(*,*) 'dB2dZ2 = ', dB2dZ2 

 666   CONTINUE
      RETURN
      END
cccc
      SUBROUTINE ZDERIV(Z,dzdg,UDYN,dUnd,dudynop,dudynpar,dZds,detol)
        IMPLICIT REAL*8(A-H,O-Z)
        dZds = 0.0d0
        dUdyn = dudynop+dudynpar
        IF(UDYN.lt.(-detol)) then
           dZds = dzdg*(dUnd-dUdyn*Z)/UDYN
         else
            dZds = 0.0d0
        endif
        RETURN
        END
